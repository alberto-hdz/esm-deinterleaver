#pragma once

#include <complex>
#include <string>
#include <vector>

#include "esm/pdw.hpp"

namespace esm {

using Iq = std::complex<float>;

// Settings for turning raw I/Q samples into PDWs. Thresholds are relative to the
// noise floor, which the detector measures from the data itself.
struct DetectorParams {
    double fs_hz = 40e6;      // sample rate
    double fc_mhz = 9000.0;   // receiver center frequency (added to the measured offset)
    int smooth = 9;           // moving-average length on |x|^2, odd so it stays centered
    double on_db = 6.0;       // pulse starts when smoothed power is this far over the noise
    double off_db = 3.0;      // ...and ends when it drops below this (hysteresis)
    int min_samples = 10;     // shorter blips are noise, not pulses
};

struct DetectorStats {
    double noise_power = 0.0;  // mean noise power per sample
    std::size_t segments = 0;  // threshold crossings before the min-length filter
};

// Raw file of interleaved float32 I, Q, I, Q, ... (what numpy's complex64 tofile writes).
std::vector<Iq> read_iq_file(const std::string& path);

// Find pulses and measure TOA, PW, RF and amplitude for each one. There's only one
// antenna, so AOA can't be measured and comes back as 0.
std::vector<Pdw> detect_pulses(const std::vector<Iq>& x, const DetectorParams& p,
                               DetectorStats* stats = nullptr);

}  // namespace esm
