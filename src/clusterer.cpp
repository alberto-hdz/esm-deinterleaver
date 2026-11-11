#include "esm/clusterer.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <vector>

namespace esm {
namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr int kUnvisited = -2;
constexpr int kNoise = -1;

bool similar(const Pdw& a, const Pdw& b, const ClusterParams& p) {
    if (std::abs(a.freq_mhz - b.freq_mhz) > p.rf_tol_mhz) return false;
    if (p.use_aoa && aoa_diff_deg(a.aoa_deg, b.aoa_deg) > p.aoa_tol_deg) return false;
    const double pw_tol = std::max(p.pw_tol_us, p.pw_tol_frac * std::max(a.pw_us, b.pw_us));
    return std::abs(a.pw_us - b.pw_us) <= pw_tol;
}

}  // namespace

double aoa_diff_deg(double a, double b) {
    double d = std::fmod(std::abs(a - b), 360.0);
    return d > 180.0 ? 360.0 - d : d;
}

double circular_mean_deg(const std::vector<double>& angles) {
    double s = 0.0;
    double c = 0.0;
    for (double a : angles) {
        s += std::sin(a * kPi / 180.0);
        c += std::cos(a * kPi / 180.0);
    }
    double m = std::atan2(s, c) * 180.0 / kPi;
    return m < 0.0 ? m + 360.0 : m;
}

std::vector<int> cluster_pdws(const std::vector<Pdw>& pdws, const ClusterParams& params) {
    const std::size_t n = pdws.size();

    // order[k] = index of the k-th pulse by RF; rank[i] = where pulse i sits in that order
    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        return pdws[a].freq_mhz < pdws[b].freq_mhz;
    });
    std::vector<std::size_t> rank(n);
    for (std::size_t k = 0; k < n; ++k) rank[order[k]] = k;

    // Walk left and right from pulse i in RF order, stopping once RF is out of tolerance.
    auto neighbors = [&](std::size_t i, std::vector<std::size_t>& out) {
        out.clear();
        const Pdw& p = pdws[i];
        for (std::size_t k = rank[i] + 1; k-- > 0;) {  // rank[i] down to 0
            const Pdw& q = pdws[order[k]];
            if (p.freq_mhz - q.freq_mhz > params.rf_tol_mhz) break;
            if (similar(p, q, params)) out.push_back(order[k]);
        }
        for (std::size_t k = rank[i] + 1; k < n; ++k) {
            const Pdw& q = pdws[order[k]];
            if (q.freq_mhz - p.freq_mhz > params.rf_tol_mhz) break;
            if (similar(p, q, params)) out.push_back(order[k]);
        }
    };

    std::vector<int> label(n, kUnvisited);
    std::vector<std::size_t> nbrs;
    std::vector<std::size_t> nbrs2;
    std::vector<std::size_t> queue;
    const auto min_pts = static_cast<std::size_t>(std::max(params.min_pts, 1));
    int next_cluster = 0;

    for (std::size_t i = 0; i < n; ++i) {
        if (label[i] != kUnvisited) continue;
        neighbors(i, nbrs);  // includes i itself
        if (nbrs.size() < min_pts) {
            label[i] = kNoise;  // may get upgraded to a border pulse later
            continue;
        }
        const int c = next_cluster++;
        label[i] = c;
        queue.clear();
        for (std::size_t j : nbrs) queue.push_back(j);
        for (std::size_t q = 0; q < queue.size(); ++q) {
            const std::size_t j = queue[q];
            if (label[j] == kNoise) label[j] = c;  // border pulse
            if (label[j] != kUnvisited) continue;
            label[j] = c;
            neighbors(j, nbrs2);
            if (nbrs2.size() >= min_pts) {  // j is a core pulse, keep expanding
                // only queue pulses not already claimed, otherwise a 1000-pulse
                // cluster puts ~1M duplicate entries on the queue
                for (std::size_t m : nbrs2) {
                    if (label[m] == kUnvisited || label[m] == kNoise) queue.push_back(m);
                }
            }
        }
    }
    return label;
}

int cluster_count(const std::vector<int>& labels) {
    int k = -1;
    for (int l : labels) k = std::max(k, l);
    return k + 1;
}

}  // namespace esm
