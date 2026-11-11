#include <gtest/gtest.h>

#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "esm/emitter_library.hpp"

namespace {

esm::LibraryEntry entry(const std::string& name, double rf_lo, double rf_hi, double pw_lo,
                        double pw_hi, esm::PriType type, double pri_lo, double pri_hi,
                        int levels = 0) {
    esm::LibraryEntry e;
    e.name = name;
    e.rf_min_mhz = rf_lo;
    e.rf_max_mhz = rf_hi;
    e.pw_min_us = pw_lo;
    e.pw_max_us = pw_hi;
    e.pri_type = type;
    e.pri_min_us = pri_lo;
    e.pri_max_us = pri_hi;
    e.stagger_levels = levels;
    return e;
}

esm::Emitter measured(double rf, double pw, esm::PriType type, double pri,
                      std::vector<double> stagger = {}) {
    esm::Emitter e;
    e.freq_mhz = rf;
    e.pw_us = pw;
    e.pri_type = type;
    e.pri_us = pri;
    e.stagger_us = std::move(stagger);
    return e;
}

// Two entries in the same band with the same pulse width. Only PRI tells them apart,
// which is the case the library has to get right.
std::vector<esm::LibraryEntry> small_library() {
    return {entry("TRACK-CONST", 9000, 9400, 0.2, 0.5, esm::PriType::Constant, 150, 250),
            entry("TRACK-JIT", 9000, 9400, 0.2, 0.5, esm::PriType::Jitter, 300, 400),
            entry("SEARCH-STAG", 2700, 2900, 1.0, 2.0, esm::PriType::Stagger, 800, 900, 3)};
}

}  // namespace

TEST(Library, ExactMatch) {
    auto r = esm::identify(measured(9200, 0.3, esm::PriType::Constant, 200), small_library(), {});
    EXPECT_EQ(r.name, "TRACK-CONST");
    EXPECT_DOUBLE_EQ(r.score, 1.0);
    EXPECT_FALSE(r.ambiguous);
}

TEST(Library, PriSeparatesSameBand) {
    auto r = esm::identify(measured(9200, 0.3, esm::PriType::Jitter, 350), small_library(), {});
    EXPECT_EQ(r.name, "TRACK-JIT");
}

TEST(Library, SlightlyOutsideWindowStillMatches) {
    // 4 MHz past the top of the band: measurement noise, not a different emitter
    auto r = esm::identify(measured(9404, 0.3, esm::PriType::Constant, 200), small_library(), {});
    EXPECT_EQ(r.name, "TRACK-CONST");
    EXPECT_LT(r.score, 1.0);
}

TEST(Library, StaggerUsesMeanPri) {
    // frame 2550 us / 3 slots = 850 us mean PRI
    auto r = esm::identify(measured(2800, 1.5, esm::PriType::Stagger, 2550, {800, 850, 900}),
                           small_library(), {});
    EXPECT_EQ(r.name, "SEARCH-STAG");
}

TEST(Library, WrongStaggerLevelsLowersScore) {
    const auto lib = small_library();
    const auto three = measured(2800, 1.5, esm::PriType::Stagger, 2550, {800, 850, 900});
    const auto four = measured(2800, 1.5, esm::PriType::Stagger, 3400, {800, 850, 900, 850});
    EXPECT_GT(esm::match_score(three, lib[2], {}), esm::match_score(four, lib[2], {}));
}

TEST(Library, OutOfBandIsUnknown) {
    auto r = esm::identify(measured(12000, 0.3, esm::PriType::Constant, 200), small_library(), {});
    EXPECT_EQ(r.name, "UNKNOWN");
}

TEST(Library, NoPriIsNeverAConfidentMatch) {
    // RF and PW fit perfectly but there's no PRI: 0.6 * (1 + 1 + 0.5) / 3 = 0.5
    auto r = esm::identify(measured(9200, 0.3, esm::PriType::Constant, 0.0), small_library(), {});
    EXPECT_EQ(r.name, "UNKNOWN");
    EXPECT_NEAR(r.score, 0.5, 1e-9);
}

TEST(Library, FlagsAmbiguousMatch) {
    std::vector<esm::LibraryEntry> lib = {
        entry("A", 9000, 9400, 0.2, 0.5, esm::PriType::Constant, 150, 250),
        entry("B", 9000, 9400, 0.2, 0.5, esm::PriType::Constant, 200, 300)};
    auto r = esm::identify(measured(9200, 0.3, esm::PriType::Constant, 220), lib, {});
    EXPECT_TRUE(r.ambiguous);
}

TEST(Library, ReadsCsv) {
    const std::string path = ::testing::TempDir() + "lib.csv";
    std::ofstream(path) << "name,role,rf_min_mhz,rf_max_mhz,pw_min_us,pw_max_us,pri_type,"
                           "pri_min_us,pri_max_us,stagger_levels,jitter_max\n"
                           "X,search,2700,2900,1,2,stagger,800,900,3,0\n";
    const auto lib = esm::read_library_csv(path);
    ASSERT_EQ(lib.size(), 1u);
    EXPECT_EQ(lib[0].name, "X");
    EXPECT_EQ(lib[0].pri_type, esm::PriType::Stagger);
    EXPECT_EQ(lib[0].stagger_levels, 3);
}

TEST(Library, BadCsvThrows) {
    const std::string path = ::testing::TempDir() + "badlib.csv";
    std::ofstream(path) << "name,role,rf_min_mhz,rf_max_mhz,pw_min_us,pw_max_us,pri_type,"
                           "pri_min_us,pri_max_us,stagger_levels,jitter_max\n"
                           "X,search,2900,2700,1,2,constant,800,900,0,0\n";  // min > max
    EXPECT_THROW(esm::read_library_csv(path), std::runtime_error);
}

TEST(Library, ProjectLibraryLoads) {
    // the real config file should always parse; CI runs from the repo root
    const auto lib = esm::read_library_csv(ESM_SOURCE_DIR "/config/emitter_library.csv");
    EXPECT_EQ(lib.size(), 12u);
}
