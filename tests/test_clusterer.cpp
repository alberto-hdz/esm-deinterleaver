#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <random>
#include <set>
#include <vector>

#include "esm/clusterer.hpp"

namespace {

// n pulses around (rf, aoa, pw) with small measurement noise
void add_emitter(std::vector<esm::Pdw>& v, int id, int n, double rf, double aoa, double pw,
                 std::mt19937& rng) {
    std::normal_distribution<double> rf_n(0.0, 2.0), aoa_n(0.0, 3.0), pw_n(0.0, 0.05);
    for (int k = 0; k < n; ++k) {
        esm::Pdw p;
        p.toa_us = 100.0 * k;
        p.freq_mhz = rf + rf_n(rng);
        p.aoa_deg = std::fmod(aoa + aoa_n(rng) + 360.0, 360.0);
        p.pw_us = std::max(pw + pw_n(rng), 0.05);
        p.truth_id = id;
        v.push_back(p);
    }
}

// Every cluster should contain exactly one truth id and every truth id one cluster.
void expect_clean_partition(const std::vector<esm::Pdw>& v, const std::vector<int>& lab,
                            int expected_clusters) {
    ASSERT_EQ(esm::cluster_count(lab), expected_clusters);
    std::map<int, std::set<int>> truth_per_cluster;
    std::map<int, std::set<int>> cluster_per_truth;
    for (std::size_t i = 0; i < v.size(); ++i) {
        ASSERT_GE(lab[i], 0) << "pulse " << i << " ended up as noise";
        truth_per_cluster[lab[i]].insert(v[i].truth_id);
        cluster_per_truth[v[i].truth_id].insert(lab[i]);
    }
    for (const auto& [c, ids] : truth_per_cluster) EXPECT_EQ(ids.size(), 1u) << "cluster " << c;
    for (const auto& [t, cs] : cluster_per_truth) EXPECT_EQ(cs.size(), 1u) << "truth " << t;
}

}  // namespace

TEST(Aoa, DiffHandlesWraparound) {
    EXPECT_DOUBLE_EQ(esm::aoa_diff_deg(359.0, 1.0), 2.0);
    EXPECT_DOUBLE_EQ(esm::aoa_diff_deg(1.0, 359.0), 2.0);
    EXPECT_DOUBLE_EQ(esm::aoa_diff_deg(10.0, 350.0), 20.0);
    EXPECT_DOUBLE_EQ(esm::aoa_diff_deg(0.0, 180.0), 180.0);
    EXPECT_DOUBLE_EQ(esm::aoa_diff_deg(45.0, 45.0), 0.0);
}

TEST(Aoa, CircularMeanAcrossNorth) {
    // a plain average of 358 and 2 gives 180, which is backwards
    // compare with aoa_diff_deg since the answer can come back as 0 or 359.9999...
    EXPECT_LT(esm::aoa_diff_deg(esm::circular_mean_deg({358.0, 2.0}), 0.0), 1e-9);
    EXPECT_LT(esm::aoa_diff_deg(esm::circular_mean_deg({350.0, 10.0, 0.0}), 0.0), 1e-9);
    EXPECT_NEAR(esm::circular_mean_deg({80.0, 100.0}), 90.0, 1e-9);
}

TEST(Clusterer, SeparatesDistinctEmitters) {
    std::mt19937 rng(1);
    std::vector<esm::Pdw> v;
    add_emitter(v, 0, 100, 9000.0, 45.0, 1.0, rng);
    add_emitter(v, 1, 100, 3000.0, 200.0, 10.0, rng);
    add_emitter(v, 2, 100, 15000.0, 300.0, 0.5, rng);
    expect_clean_partition(v, esm::cluster_pdws(v, {}), 3);
}

TEST(Clusterer, EmitterAtNorthStaysOneCluster) {
    std::mt19937 rng(2);
    std::vector<esm::Pdw> v;
    add_emitter(v, 0, 200, 9000.0, 0.0, 1.0, rng);  // AOA noise straddles 0/360
    expect_clean_partition(v, esm::cluster_pdws(v, {}), 1);
}

TEST(Clusterer, SameRfDifferentAoaSplits) {
    std::mt19937 rng(3);
    std::vector<esm::Pdw> v;
    add_emitter(v, 0, 100, 9000.0, 30.0, 1.0, rng);
    add_emitter(v, 1, 100, 9000.0, 120.0, 1.0, rng);
    expect_clean_partition(v, esm::cluster_pdws(v, {}), 2);
}

TEST(Clusterer, SameRfAndAoaNeedsPwToSplit) {
    std::mt19937 rng(4);
    std::vector<esm::Pdw> v;
    add_emitter(v, 0, 100, 9000.0, 30.0, 1.0, rng);
    add_emitter(v, 1, 100, 9000.0, 30.0, 5.0, rng);
    expect_clean_partition(v, esm::cluster_pdws(v, {}), 2);
}

TEST(Clusterer, LongPulsesUseRelativePwTolerance) {
    // 40 vs 43 us (7.5%) should merge, 40 vs 50 us (25%) should not.
    std::vector<esm::Pdw> v;
    for (int k = 0; k < 10; ++k) v.push_back({100.0 * k, 9000.0, k % 2 ? 40.0 : 43.0, -40.0, 30.0, 0});
    EXPECT_EQ(esm::cluster_count(esm::cluster_pdws(v, {})), 1);
    for (auto& p : v) p.pw_us = (p.pw_us == 43.0) ? 50.0 : 40.0;
    EXPECT_EQ(esm::cluster_count(esm::cluster_pdws(v, {})), 2);
}

TEST(Clusterer, SparsePulsesAreNoise) {
    std::vector<esm::Pdw> v;
    for (int k = 0; k < 4; ++k) v.push_back({100.0 * k, 9000.0 + 1000.0 * k, 1.0, -40.0, 30.0, -1});
    const auto lab = esm::cluster_pdws(v, {});
    EXPECT_EQ(esm::cluster_count(lab), 0);
    EXPECT_TRUE(std::all_of(lab.begin(), lab.end(), [](int l) { return l == -1; }));
}

TEST(Clusterer, IgnoresAoaWhenDisabled) {
    std::mt19937 rng(5);
    std::vector<esm::Pdw> v;
    add_emitter(v, 0, 100, 9000.0, 30.0, 1.0, rng);
    add_emitter(v, 0, 100, 9000.0, 120.0, 1.0, rng);  // same id: should merge without AOA
    esm::ClusterParams p;
    p.use_aoa = false;
    expect_clean_partition(v, esm::cluster_pdws(v, p), 1);
}

TEST(Clusterer, InputOrderDoesNotChangeResult) {
    std::mt19937 rng(6);
    std::vector<esm::Pdw> v;
    add_emitter(v, 0, 100, 9000.0, 45.0, 1.0, rng);
    add_emitter(v, 1, 100, 9050.0, 45.0, 1.0, rng);
    add_emitter(v, 2, 100, 4000.0, 270.0, 2.0, rng);
    std::shuffle(v.begin(), v.end(), rng);
    expect_clean_partition(v, esm::cluster_pdws(v, {}), 3);
}

TEST(Clusterer, EmptyInput) {
    EXPECT_TRUE(esm::cluster_pdws({}, {}).empty());
}
