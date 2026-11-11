#include "esm/pulse_detector.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace esm {
namespace {

constexpr double kPi = 3.14159265358979323846;

// Linear interpolation of where the signal crosses `level` between samples i and i + 1.
// This is what gets TOA below one sample period (25 ns at 40 MSPS).
double crossing(const std::vector<double>& s, std::size_t i, double level) {
    const double a = s[i];
    const double b = s[i + 1];
    if (b == a) return static_cast<double>(i);
    return static_cast<double>(i) + (level - a) / (b - a);
}

}  // namespace

std::vector<Iq> read_iq_file(const std::string& path) {
    std::ifstream in(path, std::ios::binary | std::ios::ate);
    if (!in) throw std::runtime_error("cannot open " + path);
    const auto bytes = static_cast<std::size_t>(in.tellg());
    if (bytes % sizeof(Iq) != 0) {
        throw std::runtime_error(path + ": size is not a whole number of complex64 samples");
    }
    std::vector<Iq> x(bytes / sizeof(Iq));
    in.seekg(0);
    in.read(reinterpret_cast<char*>(x.data()), static_cast<std::streamsize>(bytes));
    return x;
}

std::vector<Pdw> detect_pulses(const std::vector<Iq>& x, const DetectorParams& p,
                               DetectorStats* stats) {
    std::vector<Pdw> out;
    const std::size_t n = x.size();
    const std::size_t w = static_cast<std::size_t>(std::max(p.smooth, 1));
    if (n < 2 * w + 2) return out;

    std::vector<double> pw(n);
    for (std::size_t i = 0; i < n; ++i) pw[i] = std::norm(x[i]);

    // Noise floor from the median instead of the mean, so the pulses themselves don't
    // pull it up. For complex Gaussian noise |x|^2 is exponential, and an exponential's
    // median is ln(2) times its mean.
    std::vector<double> tmp(pw);
    std::nth_element(tmp.begin(), tmp.begin() + static_cast<long>(n / 2), tmp.end());
    const double noise = tmp[n / 2] / std::log(2.0);

    // Centered moving average of the power. Per-sample power of noise jumps around a
    // lot (exponential distribution), so thresholding it directly would give constant
    // false alarms. Averaging 9 samples settles it down enough for a 6 dB threshold.
    std::vector<double> sm(n, 0.0);
    const std::size_t half = w / 2;
    double run = 0.0;
    for (std::size_t i = 0; i < w; ++i) run += pw[i];
    for (std::size_t c = half; c + half + 1 < n; ++c) {
        sm[c] = run / static_cast<double>(w);
        run += pw[c + half + 1] - pw[c - half];
    }

    const double on = noise * std::pow(10.0, p.on_db / 10.0);
    const double off = noise * std::pow(10.0, p.off_db / 10.0);
    std::size_t segments = 0;

    std::size_t i = half;
    while (i + half + 1 < n) {
        if (sm[i] < on) {
            ++i;
            continue;
        }
        // Found a start. Hysteresis: stay "in the pulse" until power drops below the
        // lower off threshold, so noise riding on top of a pulse can't chop it in two.
        const std::size_t start = i;
        while (i + half + 1 < n && sm[i] >= off) ++i;
        const std::size_t stop = i;  // first sample back below off
        ++segments;
        if (stop - start < static_cast<std::size_t>(p.min_samples)) continue;
        if (start <= half + 1 || stop + half + 2 >= n) continue;  // cut off by the capture edge

        // Plateau power: median of the smoothed power away from the edges.
        const std::size_t a = std::min(start + w, stop);
        const std::size_t b = std::max(stop > w ? stop - w : 0, a + 1);
        std::vector<double> mid(sm.begin() + static_cast<long>(a), sm.begin() + static_cast<long>(b));
        std::nth_element(mid.begin(), mid.begin() + static_cast<long>(mid.size() / 2), mid.end());
        const double plateau = mid[mid.size() / 2];

        // TOA and PW from the 50% points between the noise floor and the plateau.
        // The 9-sample average is good for deciding whether there's a pulse, but it
        // turns a one-sample edge into a 9-sample ramp, and noise on that ramp moved the
        // TOA by ~20 ns RMS at 12 dB SNR. For the edges I use a 3-sample average of the
        // amplitude instead, which keeps the edge sharp. e[k - ea] is that smoothed
        // amplitude at sample k, only around this pulse.
        const std::size_t ea = start > w ? start - w : 1;
        const std::size_t eb = std::min(stop + w, n - 1);
        std::vector<double> e(eb - ea, 0.0);
        for (std::size_t k = ea; k < eb; ++k) {
            e[k - ea] = (std::abs(x[k - 1]) + std::abs(x[k]) + std::abs(x[k + 1])) / 3.0;
        }
        const double noise_amp = std::sqrt(noise * kPi / 4.0);  // mean of a Rayleigh amplitude
        const double plateau_amp = std::sqrt(std::max(plateau - noise, 0.0));
        const double level = 0.5 * (plateau_amp + noise_amp);

        // Rising edge: scanning in from the left, the first sample that's above the
        // level and stays above it for the next sample too. Falling edge: the same from
        // the right. I first searched outward from the middle of the pulse, but at 8 dB
        // SNR the amplitude dips under the level inside the pulse now and then, and the
        // search stopped there: TOA came out over a microsecond late.
        const std::size_t len = e.size();
        std::size_t r = 1;
        while (r + 1 < len && !(e[r] >= level && e[r + 1] >= level)) ++r;
        std::size_t f = len - 1;
        while (f > r + 1 && !(e[f - 1] >= level && e[f - 2] >= level)) --f;
        if (r + 1 >= len || f <= r + 1) continue;
        const double t_rise = static_cast<double>(ea) + crossing(e, r - 1, level);
        const double t_fall = static_cast<double>(ea) + crossing(e, f - 1, level);
        r += ea;  // back to absolute sample indices for the frequency step below
        f += ea;

        // Frequency from the phase step between neighboring samples, summed over the
        // middle of the pulse: angle(sum x[k] * conj(x[k-1])) * fs / 2pi. Summing the
        // complex products first and taking one angle at the end is much less noisy
        // than averaging the angle of each step. I skip a few samples at each end,
        // where the pulse is still ramping.
        auto lo = static_cast<std::size_t>(std::ceil(t_rise)) + half;
        auto hi = static_cast<std::size_t>(std::floor(t_fall));
        hi = hi > half ? hi - half : 0;
        if (hi <= lo + 1) {  // very short pulse: use all of it
            lo = r;
            hi = f;
        }
        std::complex<double> acc(0.0, 0.0);
        for (std::size_t k = std::max(lo, std::size_t{1}); k < hi; ++k) {
            acc += std::complex<double>(x[k]) * std::conj(std::complex<double>(x[k - 1]));
        }
        const double f_off_hz = std::arg(acc) * p.fs_hz / (2.0 * kPi);

        Pdw d;
        d.toa_us = t_rise / p.fs_hz * 1e6;
        d.pw_us = (t_fall - t_rise) / p.fs_hz * 1e6;
        d.freq_mhz = p.fc_mhz + f_off_hz / 1e6;
        d.amp_db = 10.0 * std::log10(std::max(plateau - noise, 1e-30));
        d.aoa_deg = 0.0;
        out.push_back(d);
    }

    if (stats) {
        stats->noise_power = noise;
        stats->segments = segments;
    }
    return out;
}

}  // namespace esm
