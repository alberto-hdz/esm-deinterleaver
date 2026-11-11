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

    const auto min_pts = static_cast<std::size_t>(std::max(params.min_pts, 1));

    // Core test: does pulse i have at least min_pts similar pulses (itself included)?
    // I walk out from i in RF order and stop as soon as I have enough. Inside a dense
    // cluster that's after a handful of pulses instead of the whole RF window.
    auto is_core = [&](std::size_t i) {
        const Pdw& p = pdws[i];
        std::size_t count = 0;
        std::size_t lo = rank[i];      // next candidate going left is lo - 1
        std::size_t hi = rank[i] + 1;  // next candidate going right is hi
        bool left = true;
        bool right = true;
        if (similar(p, p, params)) ++count;
        while ((left || right) && count < min_pts) {
            if (left) {
                if (lo == 0 || p.freq_mhz - pdws[order[lo - 1]].freq_mhz > params.rf_tol_mhz) {
                    left = false;
                } else if (similar(p, pdws[order[--lo]], params)) {
                    ++count;
                }
            }
            if (right && count < min_pts) {
                if (hi >= n || pdws[order[hi]].freq_mhz - p.freq_mhz > params.rf_tol_mhz) {
                    right = false;
                } else if (similar(p, pdws[order[hi++]], params)) {
                    ++count;
                }
            }
        }
        return count >= min_pts;
    };

    // "Open" pulses (not in a cluster yet) kept as two skip lists over the RF order,
    // so growing a cluster only visits pulses that could still join it. Without this,
    // every pulse in a big cluster re-scanned the whole RF window, which made the
    // clustering O(n^2): 15 s for 67k pulses. nxt[k] is the first open position >= k,
    // prv[k] is the last open position < k (shifted by one so 0 can mean "none").
    std::vector<std::size_t> nxt(n + 1);
    std::vector<std::size_t> prv(n + 1);
    std::iota(nxt.begin(), nxt.end(), std::size_t{0});
    std::iota(prv.begin(), prv.end(), std::size_t{0});
    auto find = [](std::vector<std::size_t>& link, std::size_t k) {
        std::size_t root = k;
        while (link[root] != root) root = link[root];
        while (link[k] != root) {  // path compression, so later lookups are O(1)
            const std::size_t up = link[k];
            link[k] = root;
            k = up;
        }
        return root;
    };
    auto close = [&](std::size_t pos) {
        nxt[pos] = pos + 1;
        prv[pos + 1] = pos;
    };

    // DBSCAN: a pulse with at least min_pts similar pulses around it is a "core" pulse.
    // Each new core pulse starts a cluster, and I grow it by pulling in its neighbors,
    // then their neighbors if they're core too, until nothing new gets added.
    // Pulses that never end up near a core pulse stay as noise; that's where the
    // spurious pulses from the simulator should land.
    std::vector<int> label(n, kUnvisited);
    std::vector<std::size_t> queue;
    int next_cluster = 0;

    for (std::size_t i = 0; i < n; ++i) {
        if (label[i] != kUnvisited) continue;
        if (!is_core(i)) {
            label[i] = kNoise;  // stays open: it can still become a border pulse later
            continue;
        }
        const int c = next_cluster++;
        label[i] = c;
        close(rank[i]);
        queue.assign(1, i);

        // A pulse goes on the queue once, when it's claimed. My first version only set
        // the label when a pulse came off the queue, so one pulse could be queued by
        // every one of its neighbors; at 336k pulses that ran out of memory (bad_alloc).
        for (std::size_t q = 0; q < queue.size(); ++q) {
            const std::size_t j = queue[q];
            if (q > 0 && !is_core(j)) continue;  // border pulse: in the cluster, but doesn't grow it
            const Pdw& p = pdws[j];
            auto claim = [&](std::size_t pos) {
                const std::size_t m = order[pos];
                if (!similar(p, pdws[m], params)) return;
                if (label[m] == kUnvisited) queue.push_back(m);  // noise pulses become border pulses
                label[m] = c;
                close(pos);
            };
            for (std::size_t pos = find(prv, rank[j]); pos > 0; pos = find(prv, pos - 1)) {
                if (p.freq_mhz - pdws[order[pos - 1]].freq_mhz > params.rf_tol_mhz) break;
                claim(pos - 1);
            }
            for (std::size_t pos = find(nxt, rank[j]); pos < n; pos = find(nxt, pos + 1)) {
                if (pdws[order[pos]].freq_mhz - p.freq_mhz > params.rf_tol_mhz) break;
                claim(pos);
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
