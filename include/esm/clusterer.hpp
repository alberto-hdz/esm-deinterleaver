#pragma once

#include <vector>

#include "esm/pdw.hpp"

namespace esm {

// Tolerances for deciding two pulses "look like" the same emitter.
// RF and AOA are absolute; PW uses whichever is bigger of the absolute and
// relative tolerance, since a 0.2 us window makes no sense for a 40 us pulse.
struct ClusterParams {
    double rf_tol_mhz = 10.0;
    double aoa_tol_deg = 10.0;
    double pw_tol_us = 0.2;
    double pw_tol_frac = 0.10;
    int min_pts = 5;        // DBSCAN: neighbors (including itself) needed to be a core pulse
    bool use_aoa = true;    // false for the single-antenna I/Q path (no AOA measurement)
};

// Smallest angle between two bearings, handles wraparound (359 vs 1 -> 2).
double aoa_diff_deg(double a, double b);

// Circular mean of bearings in degrees, result in [0, 360).
double circular_mean_deg(const std::vector<double>& angles);

// DBSCAN on (RF, AOA, PW). Returns one label per input pulse:
// 0..K-1 = cluster number, -1 = noise (didn't fit any cluster).
// I sort by RF first so each neighbor search only walks the pulses within
// rf_tol_mhz instead of comparing every pair.
std::vector<int> cluster_pdws(const std::vector<Pdw>& pdws, const ClusterParams& params);

// Number of clusters in a label vector (max label + 1).
int cluster_count(const std::vector<int>& labels);

}  // namespace esm
