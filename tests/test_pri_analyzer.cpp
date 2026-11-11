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

// Known limitation of plain SDIF with a fixed +/-1.5 us tolerance: a +/-5% jittered
// PRI spreads over ~80 bins and never clears the threshold. Phase 5 fixes this;
// I'll enable this test then.
TEST(Deinterleave, DISABLED_FindsJitteredTrain) {
    auto s = make_stream({800.0, 310.0}, 0.10, 6, 200000.0, {0.05, 0.0});
    auto trains = esm::deinterleave(s.toa, {});
    bool found = false;
    for (const auto& tr : trains) {
        int who = -1;
        if (purity(tr, s, &who) >= 0.9 && who == 0 && tr.pulses.size() > 150) found = true;
    }
    EXPECT_TRUE(found);
}
