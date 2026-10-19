// daqcampaign: abre, cierra y muestra campañas de daqrun (diseño 4.4).
//   daqcampaign [--sim] [--raw-root DIR] [--manifests DIR] [--operator NOMBRE]
//               [--authorize "<motivo>"] open|close|status NOMBRE
// Con --sim: /data/daq/sim/raw y /data/daq/sim/manifests (por omisión), para runs con digitizador
// simulado. Sin --sim: /data/daq/raw y /data/daq/manifests; abrir o cerrar una campaña real exige
// --authorize con el texto de la autorización (queda en el log y en el marcador).
#include <cstdio>
#include <cstdlib>
#include <string>

#include "daq/core/clock.h"
#include "daq/core/logger.h"
#include "daq/daqrun/campaign.h"
#include "daq/daqrun/config.h"

using namespace daq;
using namespace daq::daqrun;

namespace {
const char* kUsage =
    "Uso: daqcampaign [--sim] [--raw-root DIR] [--manifests DIR] [--operator NOMBRE] [--authorize \"motivo\"]\n"
    "                 open|close|status NOMBRE\n";
}

int main(int argc, char** argv) {
  Logger::setThreadName("daqcampaign");
  bool sim = false;
  std::string rawRoot, manifests, op, authorize, cmd, name;
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    if (a == "--sim") sim = true;
    else if (a == "--raw-root" && i + 1 < argc) rawRoot = argv[++i];
    else if (a == "--manifests" && i + 1 < argc) manifests = argv[++i];
    else if (a == "--operator" && i + 1 < argc) op = argv[++i];
    else if (a == "--authorize" && i + 1 < argc) authorize = argv[++i];
    else if (a == "--help" || a == "-h") {
      std::fputs(kUsage, stdout);
      return 0;
    } else if (cmd.empty()) cmd = a;
    else if (name.empty()) name = a;
    else {
      std::fprintf(stderr, "daqcampaign: argumento de más: %s\n%s", a.c_str(), kUsage);
      return 2;
    }
  }
  if (cmd.empty() || name.empty() || (cmd != "open" && cmd != "close" && cmd != "status")) {
    std::fputs(kUsage, stderr);
    return 2;
  }
  if (rawRoot.empty()) rawRoot = sim ? kSimRawRoot : kRealRawRoot;
  if (manifests.empty()) manifests = sim ? kSimManifests : kRealManifests;
  const CampaignLocation loc{rawRoot, manifests, name};
  SystemClock clock;
  if (cmd == "status") {
    std::printf("%s\n", campaignStatus(loc).dump(2).c_str());
    return 0;
  }
  const bool real = !sim;
  if (real && authorize.empty()) {
    std::fprintf(stderr, "daqcampaign: %s de una campaña real en %s exige --authorize \"<motivo>\"\n", cmd.c_str(),
                 rawRoot.c_str());
    return 2;
  }
  if (real) logWarn("campaña REAL ", name, ": ", cmd, " autorizado: ", authorize);
  Status st = cmd == "open" ? openCampaign(loc, op, sim, clock.now_utc_ns()) : closeCampaign(loc, clock.now_utc_ns());
  if (!st) {
    std::fprintf(stderr, "daqcampaign: %s\n", st.error().message.c_str());
    return 1;
  }
  std::printf("%s\n", campaignStatus(loc).dump(2).c_str());
  return 0;
}
