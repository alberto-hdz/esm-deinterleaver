#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <random>
#include <vector>

#include "esm/pri_analyzer.hpp"

namespace {

struct Stream {
    std::vector<double> toa;
    std::vector<int> id;
};

// Interleave several pulse trains into one sorted stream, with TOA noise and dropped pulses.
// jitter is a fraction of the PRI (0.05 = +/-5%).
Stream make_stream(const std::vector<double>& pris, double drop, unsigned seed,
                   double T_us = 200000.0, const std::vector<double>& jitter = {}) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u01(0.0, 1.0);
    std::normal_distribution<double> toa_noise(0.0, 0.05);
    std::vector<std::pair<double, int>> all;
    for (std::size_t e = 0; e < pris.size(); ++e) {
        const double j = e < jitter.size() ? jitter[e] : 0.0;
        double t = u01(rng) * pris[e];
        while (t < T_us) {
            if (u01(rng) >= drop) all.push_back({t + toa_noise(rng), static_cast<int>(e)});
            t += pris[e] * (1.0 + j * (2.0 * u01(rng) - 1.0));
        }
    }
    std::sort(all.begin(), all.end());
    Stream s;
    for (const auto& [t, id] : all) {
        s.toa.push_back(t);
        s.id.push_back(id);
    }
    return s;
}

// Fraction of a train's pulses that belong to its most common truth id.
double purity(const esm::PulseTrain& tr, const Stream& s, int* majority = nullptr) {
    std::map<int, int> count;
    for (std::size_t j : tr.pulses) ++count[s.id[j]];
    auto best = std::max_element(count.begin(), count.end(),
                                 [](const auto& a, const auto& b) { return a.second < b.second; });
    if (majority) *majority = best->first;
    return static_cast<double>(best->second) / static_cast<double>(tr.pulses.size());
}

}  // namespace

TEST(Histogram, SingleTrainPeaksAtPri) {
    auto s = make_stream({500.0}, 0.0, 1);
    std::vector<bool> avail(s.toa.size(), true);
    esm::PriParams p;
    auto h = esm::difference_histogram(s.toa, avail, 1, p);
    auto peak = static_cast<std::size_t>(std::max_element(h.begin(), h.end()) - h.begin());
    EXPECT_NEAR(static_cast<double>(peak), 500.0, 1.0);
    EXPECT_GT(h[peak], esm::sdif_threshold(peak, 1, s.toa.size(), p));
}

TEST(SequenceSearch, ToleratesDroppedPulses) {
    // 50 pulses at 400 us with pulses 10, 11 and 30 missing
    std::vector<double> toa;
    for (int k = 0; k < 50; ++k) {
        if (k != 10 && k != 11 && k != 30) toa.push_back(400.0 * k);
    }
    std::vector<bool> avail(toa.size(), true);
    auto trains = esm::sequence_search(toa, avail, 400.0, {});
    ASSERT_EQ(trains.size(), 1u);
    EXPECT_EQ(trains[0].pulses.size(), toa.size());
}

TEST(SequenceSearch, SkipsUnavailablePulses) {
    std::vector<double> toa;
    for (int k = 0; k < 20; ++k) toa.push_back(400.0 * k);
    std::vector<bool> avail(toa.size(), true);
    for (std::size_t k = 0; k < toa.size(); k += 2) avail[k] = false;  // half already used
    // the pulses still available repeat every 800 us
    auto trains = esm::sequence_search(toa, avail, 800.0, {});
    ASSERT_EQ(trains.size(), 1u);
    EXPECT_EQ(trains[0].pulses.size(), 10u);
    for (std::size_t j : trains[0].pulses) EXPECT_TRUE(avail[j]);
}

TEST(Deinterleave, TwoInterleavedTrains) {
    auto s = make_stream({310.0, 457.0}, 0.10, 2);
    auto trains = esm::deinterleave(s.toa, {});
    ASSERT_EQ(trains.size(), 2u);
    for (const auto& tr : trains) EXPECT_GE(purity(tr, s), 0.95);
}

TEST(Deinterleave, FindsFundamentalNotHarmonic) {
    // with 30% of pulses dropped, 2x PRI shows up strong in the histogram
    auto s = make_stream({310.0}, 0.30, 3);
    auto trains = esm::deinterleave(s.toa, {});
    ASSERT_FALSE(trains.empty());
    EXPECT_NEAR(trains[0].pri_us, 310.0, 0.5);
}

TEST(Deinterleave, SixTrainsMostlyRecovered) {
    // Same case as my Python prototype: six constant-PRI emitters, 10% pulse loss.
    const std::vector<double> pris = {310.0, 457.0, 1023.0, 1250.0, 733.0, 2210.0};
    auto s = make_stream(pris, 0.10, 4);
    auto trains = esm::deinterleave(s.toa, {});
    std::map<int, std::size_t> recovered;  // truth id -> pulses in its cleanest train
    for (const auto& tr : trains) {
        int who = -1;
        if (purity(tr, s, &who) >= 0.95) recovered[who] = std::max(recovered[who], tr.pulses.size());
    }
    EXPECT_GE(recovered.size(), 5u);
}

TEST(Deinterleave, RandomPulsesGiveNoTrains) {
    // Poisson arrivals have no PRI, so nothing should come out
    std::mt19937 rng(5);
    std::exponential_distribution<double> gap(1.0 / 300.0);
    std::vector<double> toa;
    for (double t = 0.0; t < 200000.0; t += gap(rng)) toa.push_back(t);
    EXPECT_TRUE(esm::deinterleave(toa, {}).empty());
}

TEST(Deinterleave, TooFewPulses) {
    EXPECT_TRUE(esm::deinterleave({}, {}).empty());
    EXPECT_TRUE(esm::deinterleave({0.0, 100.0, 200.0}, {}).empty());
}

// This test started out disabled: a +/-5% jittered PRI spreads over ~80 histogram
// bins and plain SDIF never locked on. The jitter pass handles it now.
TEST(Deinterleave, FindsJitteredTrain) {
    auto s = make_stream({800.0, 310.0}, 0.10, 6, 200000.0, {0.05, 0.0});
    auto trains = esm::deinterleave(s.toa, {});
    bool found = false;
    for (const auto& tr : trains) {
        int who = -1;
        if (purity(tr, s, &who) >= 0.9 && who == 0 && tr.pulses.size() > 150) {
            found = true;
            EXPECT_EQ(tr.type, esm::PriType::Jitter);
            EXPECT_NEAR(tr.pri_us, 800.0, 8.0);
            EXPECT_NEAR(tr.jitter, 0.05, 0.015);
        }
    }
    EXPECT_TRUE(found);
}

TEST(Deinterleave, ConstantTrainIsNotCalledJitter) {
    auto s = make_stream({640.0}, 0.10, 7);
    auto trains = esm::deinterleave(s.toa, {});
    ASSERT_EQ(trains.size(), 1u);
    EXPECT_EQ(trains[0].type, esm::PriType::Constant);
}

TEST(JitterEstimate, MeasuresPriAndSpread) {
    auto s = make_stream({1000.0}, 0.10, 8, 400000.0, {0.08});
    std::vector<bool> avail(s.toa.size(), true);
    auto est = esm::estimate_jittered_pri(s.toa, avail, {});
    EXPECT_NEAR(est.pri_us, 1000.0, 10.0);
    EXPECT_NEAR(est.jitter, 0.08, 0.015);
}

TEST(JitterEstimate, RandomPulsesGiveNothing) {
    std::mt19937 rng(9);
    std::exponential_distribution<double> gap(1.0 / 300.0);
    std::vector<double> toa;
    for (double t = 0.0; t < 200000.0; t += gap(rng)) toa.push_back(t);
    std::vector<bool> avail(toa.size(), true);
    EXPECT_EQ(esm::estimate_jittered_pri(toa, avail, {}).pri_us, 0.0);
}

namespace {

// Staggered train: the intervals repeat in order, e.g. 300, 400, 500, 300, 400, ...
Stream make_stagger(const std::vector<double>& intervals, double drop, unsigned seed,
                    double T_us = 200000.0, Stream base = {}) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u01(0.0, 1.0);
    std::vector<std::pair<double, int>> all;
    for (std::size_t i = 0; i < base.toa.size(); ++i) all.push_back({base.toa[i], base.id[i]});
    const int id = base.toa.empty() ? 0 : 1 + *std::max_element(base.id.begin(), base.id.end());
    double t = u01(rng) * intervals[0];
    for (std::size_t k = 0; t < T_us; ++k) {
        if (u01(rng) >= drop) all.push_back({t, id});
        t += intervals[k % intervals.size()];
    }
    std::sort(all.begin(), all.end());
    Stream s;
    for (const auto& [tt, ii] : all) {
        s.toa.push_back(tt);
        s.id.push_back(ii);
    }
    return s;
}

// True if b is a rotation of a (stagger order has no fixed starting point) within tol.
bool same_cycle(const std::vector<double>& a, const std::vector<double>& b, double tol) {
    if (a.size() != b.size()) return false;
    for (std::size_t r = 0; r < a.size(); ++r) {
        bool ok = true;
        for (std::size_t i = 0; i < a.size() && ok; ++i) ok = std::abs(a[i] - b[(i + r) % b.size()]) <= tol;
        if (ok) return true;
    }
    return false;
}

}  // namespace

TEST(Deinterleave, StaggerBecomesOneEmitter) {
    auto s = make_stagger({300.0, 400.0, 500.0}, 0.10, 10);
    auto trains = esm::deinterleave(s.toa, {});
    ASSERT_EQ(trains.size(), 1u);
    EXPECT_EQ(trains[0].type, esm::PriType::Stagger);
    EXPECT_NEAR(trains[0].pri_us, 1200.0, 0.5);
    EXPECT_TRUE(same_cycle(trains[0].stagger_us, {300.0, 400.0, 500.0}, 1.0));
    EXPECT_GE(trains[0].pulses.size(), s.toa.size() * 95 / 100);
}

TEST(Deinterleave, StaggerNextToConstant) {
    auto s = make_stagger({300.0, 400.0, 500.0}, 0.10, 11, 200000.0, make_stream({457.0}, 0.10, 12));
    auto trains = esm::deinterleave(s.toa, {});
    ASSERT_EQ(trains.size(), 2u);
    int stagger = 0;
    int constant = 0;
    for (const auto& tr : trains) {
        EXPECT_GE(purity(tr, s), 0.95);
        stagger += tr.type == esm::PriType::Stagger;
        constant += tr.type == esm::PriType::Constant;
    }
    EXPECT_EQ(stagger, 1);
    EXPECT_EQ(constant, 1);
}

TEST(Deinterleave, LongStaggerFrame) {
    // 4-slot stagger with a ~9.8 ms frame. My first version couldn't see this at all
    // because max_pri_us was 3000.
    const std::vector<double> iv = {2631.0, 2288.6, 2919.5, 1989.1};
    auto s = make_stagger(iv, 0.05, 13, 400000.0);
    auto trains = esm::deinterleave(s.toa, {});
    ASSERT_EQ(trains.size(), 1u);
    EXPECT_EQ(trains[0].type, esm::PriType::Stagger);
    EXPECT_TRUE(same_cycle(trains[0].stagger_us, iv, 1.0));
}

TEST(Deinterleave, ConstantNextToJitterAtHighLoss) {
    // The case the identification step caught: at 30% pulse loss a constant train
    // breaks into pieces, and with the histogram PRI the pieces drifted apart and got
    // merged into a fake stagger. The least-squares PRI fit fixed it.
    auto s = make_stream({245.0, 383.0}, 0.30, 14, 200000.0, {0.0, 0.045});
    auto trains = esm::deinterleave(s.toa, {});
    bool found = false;
    for (const auto& tr : trains) {
        int who = -1;
        if (purity(tr, s, &who) >= 0.95 && who == 0) {
            found = true;
            EXPECT_EQ(tr.type, esm::PriType::Constant);
            EXPECT_NEAR(tr.pri_us, 245.0, 0.01);
        }
    }
    EXPECT_TRUE(found);
}
