#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace esm {

// Settings for the PRI (pulse repetition interval) search. The SDIF threshold
// shape comes from Milojevic & Popovic (1992); I tuned the constants on my own
// simulated scenarios.
struct PriParams {
    double bin_us = 1.0;          // histogram bin width
    double max_pri_us = 12000.0;  // longest interval I look for (covers stagger frames)
    int max_level = 5;            // highest difference level c (gap of c pulses)
    double gamma = 0.15;          // threshold height
    double decay_us = 900.0;      // threshold falls off as exp(-tau / decay_us)
    double tol_us = 1.5;          // sequence search: how close a pulse must land to the prediction
    double tol_growth = 0.0;      // sequence search: widen tol by this much per missed pulse
    int max_misses = 3;           // sequence search: missed pulses allowed in a row
    int min_hits = 10;            // sequence search: shortest train I accept
    double min_fill = 0.6;        // sequence search: found pulses / expected pulses
    double min_merged_fill = 0.4; // same check after broken pieces are glued back together
    int min_pulses = 10;          // stop when fewer pulses than this are left
    double max_jitter = 0.12;     // jitter pass: widest jitter I try to follow (fraction of PRI)
    double jitter_share = 0.4;    // jitter pass: share of gaps that must agree on one PRI
};

enum class PriType : std::uint8_t { Constant, Stagger, Jitter };

// One pulse train pulled out of the stream (indices into the TOA vector).
// For a stagger, pri_us is the full frame period and stagger_us holds the
// intervals inside one frame.
struct PulseTrain {
    PriType type = PriType::Constant;
    double pri_us = 0.0;
    double jitter = 0.0;              // estimated +/- fraction, jitter trains only
    std::vector<double> stagger_us;   // stagger trains only
    std::vector<std::size_t> pulses;
};

// Histogram of toa[n + level] - toa[n], using only pulses where available[n] is true.
// Each bin is added to its two neighbors so a PRI sitting on a bin edge doesn't get
// split in half.
std::vector<double> difference_histogram(const std::vector<double>& toa,
                                         const std::vector<bool>& available, int level,
                                         const PriParams& p);

// SDIF threshold for bin k at difference level c with n pulses remaining.
double sdif_threshold(std::size_t k, int level, std::size_t n, const PriParams& p);

// Follow a train with the given PRI through the available pulses. Returns every
// train with at least min_hits pulses that fills at least min_fill of its expected
// slots. It does NOT mark them used; deinterleave() does that.
std::vector<PulseTrain> sequence_search(const std::vector<double>& toa,
                                        const std::vector<bool>& available, double pri_us,
                                        const PriParams& p);

// Several trains at the same PRI with different phases are one staggered emitter.
// Merges them and works out the intervals inside the frame.
PulseTrain merge_stagger(const std::vector<double>& toa, const std::vector<PulseTrain>& trains,
                         double frame_us);

// Jitter pass helper: if most gaps between the available pulses agree on one PRI
// within +/- max_jitter, returns {pri, jitter estimate}. Otherwise returns {0, 0}.
struct JitterEstimate {
    double pri_us = 0.0;
    double jitter = 0.0;
};
JitterEstimate estimate_jittered_pri(const std::vector<double>& toa,
                                     const std::vector<bool>& available, const PriParams& p);

// Full deinterleaver. toa must be sorted ascending.
// Pass 1 is the SDIF loop (constant and staggered PRIs).
// Pass 2 follows jittered PRIs in whatever pulses are left.
std::vector<PulseTrain> deinterleave(const std::vector<double>& toa, const PriParams& p);

const char* pri_type_name(PriType t);

}  // namespace esm
