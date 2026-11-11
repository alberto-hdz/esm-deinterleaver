#pragma once

#include <cstddef>
#include <vector>

namespace esm {

// Settings for the PRI (pulse repetition interval) search. The threshold
// constants gamma and theta come from the SDIF paper (Milojevic & Popovic, 1992);
// I tuned the starting values on my own simulated scenarios.
struct PriParams {
    double bin_us = 1.0;        // histogram bin width
    double max_pri_us = 3000.0; // longest PRI I look for
    int max_level = 5;          // highest difference level c (gap of c pulses)
    double gamma = 0.15;        // threshold height
    double theta = 0.3;         // how fast the threshold decays with tau
    double tol_us = 1.5;        // sequence search: how close a pulse must land to the prediction
    int max_misses = 3;         // sequence search: missed pulses allowed in a row
    int min_hits = 10;          // sequence search: shortest train I accept
    double min_fill = 0.6;      // sequence search: found pulses / expected pulses
    int min_pulses = 10;        // stop when fewer pulses than this are left
};

// One pulse train pulled out of the stream: its PRI and which pulses belong to it
// (indices into the TOA vector passed to deinterleave()).
struct PulseTrain {
    double pri_us = 0.0;
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

// Full SDIF loop. toa must be sorted ascending.
std::vector<PulseTrain> deinterleave(const std::vector<double>& toa, const PriParams& p);

}  // namespace esm
