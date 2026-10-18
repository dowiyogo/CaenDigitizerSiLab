// Tests de conversión de unidades y de decodificación de Status.
#include <gtest/gtest.h>

#include <nlohmann/json.hpp>

#include "daq/core/units.h"
#include "daq/hv/hv_types.h"
#include "daq/hv/status_bits.h"

using daq::hv::StatusBitMap;

TEST(Units, MicroampereToAmpere) {
  EXPECT_DOUBLE_EQ(daq::units::microampereToAmpere(15.0), 15e-6);
  EXPECT_DOUBLE_EQ(daq::units::ampereToMicroampere(2e-6), 2.0);
  EXPECT_DOUBLE_EQ(daq::units::volt(65.75), 65.75);
  EXPECT_NEAR(daq::units::applyDecimalExponent(-0.048, -6), -4.8e-8, 1e-20);
  EXPECT_DOUBLE_EQ(daq::units::applyDecimalExponent(65.75, 0), 65.75);
  EXPECT_NEAR(daq::units::removeDecimalExponent(1.5e-5, -6), 15.0, 1e-12);
}

TEST(Units, SymbolsAndPrefixes) {
  EXPECT_EQ(daq::units::unitSymbol(1), "A");
  EXPECT_EQ(daq::units::unitSymbol(2), "V");
  EXPECT_EQ(daq::units::unitSymbol(7), "V/s");
  EXPECT_EQ(daq::units::unitSymbol(8), "s");
  EXPECT_EQ(daq::units::decimalPrefix(-6), "u");
  EXPECT_EQ(daq::units::decimalPrefix(0), "");
}

TEST(Units, NumericValueConvertsAtBoundary) {
  daq::hv::ParamInfo imon;
  imon.name = "IMon";
  imon.type = daq::hv::ParamType::Numeric;
  imon.unitCode = 1;
  imon.exponent = -6;
  auto v = daq::hv::numericValue(imon, 15.0);
  EXPECT_DOUBLE_EQ(v.boardValue, 15.0);
  EXPECT_NEAR(v.si, 15e-6, 1e-18);
  EXPECT_EQ(imon.siUnit(), "A");
  EXPECT_EQ(imon.boardUnit(), "uA");
}

TEST(StatusBits, DefaultMapDecodesWrapperRev13Bits) {
  const auto m = StatusBitMap::wrapperRev13();
  EXPECT_TRUE(m.decode(0).empty());
  EXPECT_EQ(m.decode(0x0001), std::vector<std::string>({"on"}));
  EXPECT_EQ(m.decode(0x0003), std::vector<std::string>({"on", "ramp_up"}));
  EXPECT_EQ(m.decode(0x0200), std::vector<std::string>({"internal_trip"}));
  EXPECT_EQ(m.decode(0x0008 | 0x0001), std::vector<std::string>({"on", "overcurrent"}));
  EXPECT_EQ(m.decode(1u << 15), std::vector<std::string>({"bit15"}));
  EXPECT_TRUE(m.has(0x0209, "internal_trip"));
  EXPECT_FALSE(m.has(0x0209, "ramp_down"));
  EXPECT_FALSE(m.has(0xFFFF, "no_existe"));
  EXPECT_EQ(m.bitOf("unplugged"), 11);
}

TEST(StatusBits, FromJsonValidates) {
  auto ok = StatusBitMap::fromJson(nlohmann::json{{"on", 0}, {"trip", 9}});
  ASSERT_TRUE(ok.ok());
  EXPECT_EQ(ok.value().decode(0x0201), std::vector<std::string>({"on", "trip"}));
  EXPECT_FALSE(StatusBitMap::fromJson(nlohmann::json{{"a", 1}, {"b", 1}}).ok());
  EXPECT_FALSE(StatusBitMap::fromJson(nlohmann::json{{"a", 32}}).ok());
  EXPECT_FALSE(StatusBitMap::fromJson(nlohmann::json{{"a", "x"}}).ok());
  EXPECT_FALSE(StatusBitMap::fromJson(nlohmann::json::array()).ok());
}
