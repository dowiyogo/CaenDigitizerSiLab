// validate_hv_flags.C — validación 4 de DAQ2-5 (solo lectura).
// Run sintético con un grupo HV apagado por un hvmon simulado: comprueba en los segmentos que
//  - durante el apagado (lejos de los cambios) el bit 0 (HV_OFF) está exactamente en los canales del
//    digitizador del grupo según el mapeo JSON, y fuera del apagado en ninguno;
//  - el bit 9 (MAP_UNCONFIRMED) está en los canales cuyo hv_group o face no es "measurement";
//  - el bit 4 (HV_UNKNOWN) está en los canales sin grupo conocido (y en ningún otro con hvmon vivo);
//  - el bit 5 (HV_AMBIGUOUS) aparece solo cerca de los cambios.
// El intervalo de apagado sale del árbol hv del propio run (pw[hvCh]).
// Uso: root -l -b -q 'daq/tools/validation/validate_hv_flags.C("<dir>/<base>_seg*.root", 3, "<mapeo.json>")'
#include <TChain.h>

#include <cstdio>
#include <fstream>
#include <set>
#include <vector>

#include <nlohmann/json.hpp>

void validate_hv_flags(const char* pattern, int group = 3,
                       const char* mappingFile = "daq/config/mapping/hodo_mapping_2026_10.json", double marginS = 3.0) {
  std::ifstream in(mappingFile);
  const nlohmann::json m = nlohmann::json::parse(in);
  std::set<int> inGroup, unknown, unconfirmed;
  for (const auto& c : m["channels"]) {
    const int ch = c["dig_ch"];
    if (c["hv_group"].is_null()) unknown.insert(ch);
    else if (c["hv_group"].get<int>() == group) inGroup.insert(ch);
    const auto& cf = c["confirmed"];
    if (cf.value("hv_group", "") != "measurement" || cf.value("face", "") != "measurement") unconfirmed.insert(ch);
  }
  const int hvCh = m["hv_group_to_sy5527_channel"][std::to_string(group)];
  printf("grupo %d = CH%02d: canales del digitizador", group, hvCh);
  for (int c : inGroup) printf(" %d", c);
  printf("; sin grupo: %zu; sin confirmar (bit 9 esperado): %zu\n", unknown.size(), unconfirmed.size());

  // Intervalo de apagado desde el árbol hv.
  TChain hv("hv");
  hv.Add(pattern);
  Long64_t tResp = 0;
  UInt_t pw[12];
  hv.SetBranchAddress("t_resp_utc_ns", &tResp);
  hv.SetBranchAddress("pw", pw);
  Long64_t tOff = 0, tOn = 0;
  UInt_t prev = 1;
  for (Long64_t i = 0; i < hv.GetEntries(); ++i) {
    hv.GetEntry(i);
    if (prev == 1 && pw[hvCh] == 0 && tOff == 0) tOff = tResp;
    if (prev == 0 && pw[hvCh] == 1 && tOff != 0 && tOn == 0) tOn = tResp;
    prev = pw[hvCh];
  }
  printf("CH%02d apagado según el árbol hv: de t_resp %lld a %lld (%.1f s)\n", hvCh, tOff, tOn, (tOn - tOff) / 1e9);

  TChain ev("events");
  ev.Add(pattern);
  Long64_t t = 0;
  UShort_t f[64];
  ev.SetBranchAddress("t_utc_ns", &t);
  ev.SetBranchAddress("ch_flags", f);
  long nOff = 0, nOn = 0, nNear = 0, badOff = 0, badOn = 0, bad9 = 0, bad4 = 0, amb = 0, ambFar = 0;
  std::set<int> ambCh;
  for (Long64_t i = 0; i < ev.GetEntries(); ++i) {
    ev.GetEntry(i);
    const bool near = (tOff && std::abs(static_cast<double>(t - tOff)) / 1e9 < marginS + 5.0) ||
                      (tOn && std::abs(static_cast<double>(t - tOn)) / 1e9 < marginS + 5.0);
    const bool off = tOff && t > tOff && (tOn == 0 || t < tOn);
    std::set<int> b0;
    bool anyAmb = false;
    for (int c = 0; c < 64; ++c) {
      if (f[c] & (1u << 0)) b0.insert(c);
      if (((f[c] >> 9) & 1u) != (unconfirmed.count(c) ? 1u : 0u)) ++bad9;
      if (((f[c] >> 4) & 1u) != (unknown.count(c) ? 1u : 0u)) ++bad4;
      if (f[c] & (1u << 5)) {
        anyAmb = true;
        ambCh.insert(c);
      }
    }
    if (anyAmb) {
      ++amb;
      if (!near) ++ambFar;
    }
    if (near) {
      ++nNear;
      continue;
    }
    if (off) {
      ++nOff;
      if (b0 != inGroup) ++badOff;
    } else {
      ++nOn;
      if (!b0.empty()) ++badOn;
    }
  }
  printf("eventos: %lld; apagado (lejos de los cambios): %ld, con bit 0 distinto del grupo: %ld\n", ev.GetEntries(), nOff, badOff);
  printf("encendido (lejos de los cambios): %ld, con algún bit 0: %ld; cerca de los cambios (±%.0f s): %ld\n", nOn, badOn,
         marginS + 5.0, nNear);
  printf("bit 9 distinto de lo esperado: %ld (canal×evento); bit 4 distinto de lo esperado: %ld\n", bad9, bad4);
  printf("eventos con bit 5: %ld (lejos de los cambios: %ld); canales con bit 5:", amb, ambFar);
  for (int c : ambCh) printf(" %d", c);
  printf("\n");
}
