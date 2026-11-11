#include <gtest/gtest.h>

#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "esm/pdw.hpp"

namespace {

std::string write_temp(const std::string& name, const std::string& text) {
    const std::string path = ::testing::TempDir() + name;
    std::ofstream(path) << text;
    return path;
}

}  // namespace

TEST(PdwIo, RoundTrip) {
    std::vector<esm::Pdw> in = {{10.5, 9000.0, 1.25, -40.0, 45.0, 0},
                                {20.25, 9001.5, 1.3, -41.5, 359.5, 1}};
    const std::string path = ::testing::TempDir() + "roundtrip.csv";
    esm::write_pdw_csv(path, in);
    const auto out = esm::read_pdw_csv(path);
    ASSERT_EQ(out.size(), in.size());
    for (std::size_t i = 0; i < in.size(); ++i) {
        EXPECT_NEAR(out[i].toa_us, in[i].toa_us, 1e-4);
        EXPECT_NEAR(out[i].freq_mhz, in[i].freq_mhz, 1e-3);
        EXPECT_NEAR(out[i].pw_us, in[i].pw_us, 1e-4);
        EXPECT_NEAR(out[i].amp_db, in[i].amp_db, 1e-2);
        EXPECT_NEAR(out[i].aoa_deg, in[i].aoa_deg, 1e-2);
        EXPECT_EQ(out[i].truth_id, in[i].truth_id);
    }
}

TEST(PdwIo, ColumnOrderDoesNotMatter) {
    const auto path = write_temp("reordered.csv",
                                 "aoa_deg,toa_us,pw_us,freq_mhz,amp_db\n"
                                 "90.0,5.0,2.0,3000.0,-50.0\n");
    const auto p = esm::read_pdw_csv(path);
    ASSERT_EQ(p.size(), 1u);
    EXPECT_DOUBLE_EQ(p[0].toa_us, 5.0);
    EXPECT_DOUBLE_EQ(p[0].freq_mhz, 3000.0);
    EXPECT_DOUBLE_EQ(p[0].aoa_deg, 90.0);
    EXPECT_EQ(p[0].truth_id, -1);  // no truth column -> stays -1
}

TEST(PdwIo, HandlesWindowsLineEndings) {
    const auto path = write_temp("crlf.csv",
                                 "toa_us,freq_mhz,pw_us,amp_db,aoa_deg,truth_id\r\n"
                                 "1.0,2000.0,0.5,-60.0,10.0,3\r\n");
    const auto p = esm::read_pdw_csv(path);
    ASSERT_EQ(p.size(), 1u);
    EXPECT_EQ(p[0].truth_id, 3);
}

TEST(PdwIo, MissingColumnThrows) {
    const auto path = write_temp("missing.csv", "toa_us,freq_mhz,pw_us,amp_db\n1,2,3,4\n");
    EXPECT_THROW(esm::read_pdw_csv(path), std::runtime_error);
}

TEST(PdwIo, BadValueReportsLineNumber) {
    const auto path = write_temp("bad.csv",
                                 "toa_us,freq_mhz,pw_us,amp_db,aoa_deg\n"
                                 "1,2,3,4,5\n"
                                 "1,abc,3,4,5\n");
    try {
        esm::read_pdw_csv(path);
        FAIL() << "expected an exception";
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find("line 3"), std::string::npos) << e.what();
    }
}

TEST(PdwIo, MissingFileThrows) {
    EXPECT_THROW(esm::read_pdw_csv("/no/such/file.csv"), std::runtime_error);
}
