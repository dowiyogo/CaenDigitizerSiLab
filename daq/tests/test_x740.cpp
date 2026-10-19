// Tests del backend CAEN que no necesitan la placa: decodificación de buffers crudos del VX1740
// capturados en DAQ2-6 (comparados con CAEN_DGTZ_DecodeEvent), tabla de errores CAEN y
// desenrollado del TTT medido (8 ns, 31 bit, bit 31 = marca de vuelta).
#include <gtest/gtest.h>

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iterator>
#include <sstream>

#include <nlohmann/json.hpp>

#include "daq/daqrun/hv_client.h"
#include "daq/dig/x740_raw.h"

using namespace daq;
using namespace daq::dig;

namespace {
std::string slurp(const std::string& p) {
  std::ifstream in(p, std::ios::binary);
  return std::string((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
}
const std::string kData = DAQ_SOURCE_DIR "/tests/data/x740/";
}  // namespace

TEST(X740Raw, DecodesCapturedBlocksLikeCaenDecodeEvent) {
  for (const char* b : {"block_000", "block_001"}) {
    const std::string raw = slurp(kData + b + ".bin");
    const auto j = nlohmann::json::parse(slurp(kData + b + ".json"));
    ASSERT_EQ(raw.size(), j["bytes"].get<std::size_t>());
    auto spans = splitX740Buffer(raw.data(), raw.size());
    ASSERT_TRUE(spans.ok()) << spans.error().message;
    ASSERT_EQ(spans.value().size(), j["events"].size());
    const int n = j["record_length"];
    EXPECT_EQ(n, 144);
    for (std::size_t i = 0; i < spans.value().size(); ++i) {
      const auto& h = spans.value()[i].header;
      const auto& e = j["events"][i];
      EXPECT_EQ(h.sizeWords * 4, e["event_size"].get<std::uint32_t>());
      EXPECT_EQ(h.counter, e["counter"].get<std::uint32_t>());
      EXPECT_EQ(h.ttt, e["ttt"].get<std::uint32_t>());
      EXPECT_EQ(h.groupMask, 0xFFu);
      EXPECT_EQ(h.boardId, e["board_id"].get<std::uint32_t>());
      std::vector<std::uint16_t> s;
      std::vector<std::uint32_t> w(raw.size() / 4);
      std::memcpy(w.data(), raw.data(), raw.size());
      auto got = decodeX740Samples(w.data() + spans.value()[i].offsetWords, h.sizeWords, s);
      ASSERT_TRUE(got.ok()) << got.error().message;
      ASSERT_EQ(got.value(), n);
      int bad = 0;
      for (int c = 0; c < 64; ++c)
        for (int k = 0; k < n; ++k) bad += s[static_cast<std::size_t>(c * n + k)] != e["adc"][c * n + k].get<int>();
      EXPECT_EQ(bad, 0) << b;
    }
  }
  // Buffers dañados: errores, no basura.
  std::uint32_t junk[8] = {0x12345678, 0, 0, 0, 0, 0, 0, 0};
  EXPECT_FALSE(splitX740Buffer(junk, sizeof junk).ok());
  std::uint32_t truncated[4] = {0xA0000100, 0xFF, 1, 2};
  EXPECT_FALSE(splitX740Buffer(truncated, sizeof truncated).ok());
}

TEST(X740Raw, CaenErrorNames) {
  EXPECT_EQ(caenErrorName(0), "código 0 (CAEN_DGTZ_Success: operación completada)");
  EXPECT_NE(caenErrorName(-1).find("CAEN_DGTZ_CommError"), std::string::npos);
  EXPECT_NE(caenErrorName(-25).find("DigitizerAlreadyOpen"), std::string::npos);
  EXPECT_NE(caenErrorName(-18).find("Timeout"), std::string::npos);
  EXPECT_NE(caenErrorName(-12).find("desconocido"), std::string::npos);
  EXPECT_NE(caenErrorName(-99).find("NotYetImplemented"), std::string::npos);
}

TEST(TttUnwrapper, RealVx1740SequenceAcrossTwoWraps) {
  // Eventos 1500-3700 del run002 (SW trigger a 100 Hz): 15-37 s desde el inicio, dos vueltas.
  std::ifstream in(kData + "ttt_run002_1500_3700.csv");
  std::string line;
  std::vector<std::pair<std::uint32_t, std::int64_t>> ev;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    const auto c = line.find(',');
    ev.emplace_back(static_cast<std::uint32_t>(std::stoul(line.substr(0, c))), std::stoll(line.substr(c + 1)));
  }
  ASSERT_EQ(ev.size(), 2200u);
  int withFlag = 0;
  for (const auto& e : ev) withFlag += (e.first >> 31) & 1;
  EXPECT_GT(withFlag, 1500);  // después de la primera vuelta el bit 31 queda en 1
  // Origen aproximado: el primer evento está antes de la primera vuelta (15 s < 17.18 s).
  const std::int64_t origin = ev[0].second - static_cast<std::int64_t>((ev[0].first & 0x7FFFFFFF)) * 8 - 2000000;
  daqrun::TttUnwrapper u(origin, 8.0, 31);
  std::int64_t prev = -1;
  double maxLag = 0;
  for (const auto& e : ev) {
    const std::int64_t t = u.tttNs(e.first, e.second);
    EXPECT_GT(t, prev);
    prev = t;
    maxLag = std::max(maxLag, std::fabs(static_cast<double>(e.second - (origin + t))));
  }
  EXPECT_LT(maxLag, 5e6);  // < 5 ms entre el TTT desenrollado y la lectura del PC
  EXPECT_NEAR(static_cast<double>(prev) / 1e9, 37.0, 0.5);
  // Sin la máscara (32 bit) el desenrollado se rompe: lo que DAQ2-6 corrigió.
  daqrun::TttUnwrapper u32(origin, 8.0, 32);
  double maxLag32 = 0;
  for (const auto& e : ev) maxLag32 = std::max(maxLag32, std::fabs(static_cast<double>(e.second - (origin + u32.tttNs(e.first, e.second)))));
  EXPECT_GT(maxLag32, 1e9);
}
