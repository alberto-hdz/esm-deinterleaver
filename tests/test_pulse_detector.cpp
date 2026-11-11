#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <random>
#include <vector>

#include "esm/pulse_detector.hpp"

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kFs = 40e6;

struct Pulse {
    double toa_us;
    double pw_us;
    double offset_mhz;
    double snr_db;  // per-sample SNR
};

// Unit-power complex noise plus rectangular CW pulses with 25 ns linear edges.
// The envelope is evaluated at the real sample times, so a pulse can start
// between samples like it would in a real receiver.
std::vector<esm::Iq> make_capture(double ms, const std::vector<Pulse>& pulses, unsigned seed,
                                  bool noise = true) {
    const auto n = static_cast<std::size_t>(ms * 1e-3 * kFs);
    std::vector<std::complex<double>> x(n, {0.0, 0.0});
    std::mt19937 rng(seed);
    std::normal_distribution<double> g(0.0, std::sqrt(0.5));
    if (noise) {
        for (auto& v : x) v = {g(rng), g(rng)};
    }
    const double tr = 25e-9;
    for (const auto& p : pulses) {
        const double amp = std::sqrt(std::pow(10.0, p.snr_db / 10.0));
        const double t0 = p.toa_us * 1e-6;
        const double t1 = t0 + p.pw_us * 1e-6;
        const auto k0 = static_cast<std::size_t>(std::max(0.0, (t0 - tr) * kFs));
        const auto k1 = std::min(n, static_cast<std::size_t>((t1 + tr) * kFs) + 2);
        for (std::size_t k = k0; k < k1; ++k) {
            const double t = static_cast<double>(k) / kFs;
            // TOA is the 50% point of the rising ramp, PW is 50% to 50%
            const double env = std::clamp((t - t0) / tr + 0.5, 0.0, 1.0) *
                               std::clamp((t1 - t) / tr + 0.5, 0.0, 1.0);
            const double ph = 2.0 * kPi * p.offset_mhz * 1e6 * t;
            x[k] += amp * env * std::complex<double>(std::cos(ph), std::sin(ph));
        }
    }
    std::vector<esm::Iq> out(n);
    for (std::size_t k = 0; k < n; ++k) out[k] = esm::Iq(static_cast<float>(x[k].real()), static_cast<float>(x[k].imag()));
    return out;
}

esm::DetectorParams params() {
    esm::DetectorParams p;
    p.fs_hz = kFs;
    p.fc_mhz = 9000.0;
    return p;
}

}  // namespace

TEST(PulseDetector, MeasuresOnePulse) {
    auto x = make_capture(1.0, {{200.013, 2.0, 6.3, 15.0}}, 1);
    auto d = esm::detect_pulses(x, params());
    ASSERT_EQ(d.size(), 1u);
    EXPECT_NEAR(d[0].toa_us, 200.013, 0.025);   // within one sample period
    EXPECT_NEAR(d[0].pw_us, 2.0, 0.05);
    EXPECT_NEAR(d[0].freq_mhz, 9006.3, 0.3);
}

TEST(PulseDetector, NegativeFrequencyOffset) {
    auto x = make_capture(1.0, {{300.0, 1.0, -12.5, 15.0}}, 2);
    auto d = esm::detect_pulses(x, params());
    ASSERT_EQ(d.size(), 1u);
    EXPECT_NEAR(d[0].freq_mhz, 8987.5, 0.3);
}

TEST(PulseDetector, NoiseOnlyGivesNothing) {
    // 2 ms of noise is 80,000 chances for a false alarm
    auto x = make_capture(2.0, {}, 3);
    esm::DetectorStats st;
    auto d = esm::detect_pulses(x, params(), &st);
    EXPECT_TRUE(d.empty());
    EXPECT_NEAR(st.noise_power, 1.0, 0.05);  // the median-based floor finds the unit noise
}

TEST(PulseDetector, NoiseFloorIgnoresPulses) {
    // 10% duty cycle of strong pulses shouldn't move the noise estimate much
    std::vector<Pulse> ps;
    for (int k = 0; k < 20; ++k) ps.push_back({50.0 + 100.0 * k, 10.0, 3.0, 25.0});
    auto x = make_capture(2.1, ps, 4);
    esm::DetectorStats st;
    auto d = esm::detect_pulses(x, params(), &st);
    EXPECT_EQ(d.size(), 20u);
    EXPECT_LT(st.noise_power, 1.2);
}

TEST(PulseDetector, WeakPulseMissedStrongFound) {
    auto x = make_capture(1.0, {{100.0, 2.0, 0.0, 0.0}, {500.0, 2.0, 0.0, 12.0}}, 5);
    auto d = esm::detect_pulses(x, params());
    ASSERT_EQ(d.size(), 1u);
    EXPECT_NEAR(d[0].toa_us, 500.0, 0.025);
}

TEST(PulseDetector, TwoPulsesCloseTogether) {
    auto x = make_capture(1.0, {{100.0, 1.0, 2.0, 15.0}, {102.0, 1.0, 2.0, 15.0}}, 6);
    auto d = esm::detect_pulses(x, params());
    ASSERT_EQ(d.size(), 2u);
    EXPECT_NEAR(d[1].toa_us - d[0].toa_us, 2.0, 0.025);
}

TEST(PulseDetector, OverlappingPulsesMerge) {
    // Known limitation: one antenna and no channelizer, so two pulses that overlap in
    // time come out as one longer pulse.
    auto x = make_capture(1.0, {{100.0, 2.0, -8.0, 15.0}, {101.0, 2.0, 8.0, 15.0}}, 7);
    auto d = esm::detect_pulses(x, params());
    ASSERT_EQ(d.size(), 1u);
    EXPECT_NEAR(d[0].pw_us, 3.0, 0.1);
}

TEST(PulseDetector, PulseAtCaptureEdgeIsSkipped) {
    auto x = make_capture(1.0, {{0.0, 5.0, 0.0, 20.0}, {500.0, 2.0, 0.0, 20.0}}, 8);
    auto d = esm::detect_pulses(x, params());
    ASSERT_EQ(d.size(), 1u);
    EXPECT_NEAR(d[0].toa_us, 500.0, 0.025);
}

TEST(PulseDetector, AccuracyOverManyPulses) {
    // 200 pulses at random sub-sample start times, 12 dB SNR
    std::mt19937 rng(9);
    std::uniform_real_distribution<double> frac(0.0, 0.025);
    std::vector<Pulse> ps;
    for (int k = 0; k < 200; ++k) ps.push_back({20.0 + 20.0 * k + frac(rng), 1.5, 4.0, 12.0});
    auto x = make_capture(4.1, ps, 10);
    auto d = esm::detect_pulses(x, params());
    ASSERT_EQ(d.size(), ps.size());
    double toa_err2 = 0.0;
    double rf_err2 = 0.0;
    for (std::size_t k = 0; k < ps.size(); ++k) {
        toa_err2 += std::pow(d[k].toa_us - ps[k].toa_us, 2);
        rf_err2 += std::pow(d[k].freq_mhz - 9004.0, 2);
    }
    EXPECT_LT(std::sqrt(toa_err2 / 200.0) * 1000.0, 15.0);  // RMS TOA error under 15 ns (0.6 sample)
    EXPECT_LT(std::sqrt(rf_err2 / 200.0), 0.2);             // RMS RF error under 200 kHz
}

TEST(PulseDetector, ReadsComplex64File) {
    const std::string path = ::testing::TempDir() + "iq.c64";
    std::vector<esm::Iq> x = {{1.0f, -2.0f}, {0.5f, 0.25f}};
    std::ofstream(path, std::ios::binary).write(reinterpret_cast<const char*>(x.data()), 16);
    auto y = esm::read_iq_file(path);
    ASSERT_EQ(y.size(), 2u);
    EXPECT_EQ(y[0], x[0]);
    EXPECT_EQ(y[1], x[1]);
}
