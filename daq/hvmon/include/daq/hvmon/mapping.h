// mapping.h
// Lectura del mapeo versionado (daq/config/mapping/hodo_mapping_2026_10.json, esquema
// daqv2.mapping/1): canal de digitizador -> estación, plano, tira, grupo HV, cara, y la procedencia
// de cada campo. Solo "measurement" cuenta como confirmado. El notificador lo usa para listar los
// canales del digitizador afectados por un canal HV.
// Tipo de valor, seguro para lectura concurrente.
#pragma once

#include <map>
#include <optional>
#include <string>
#include <vector>

#include "daq/core/result.h"

namespace daq::hvmon {

struct MappingEntry {
  int digCh = -1;
  std::string station;  // hod2018 | hod2019
  std::string plane;    // X | Y
  int strip = -1;
  std::optional<int> hvGroup;
  std::optional<int> hvCh;  // canal del SY5527 (= grupo - 1), si el grupo se conoce
  std::string face;
  std::string hvGroupConfirmed;  // measurement | cabling_record | derived | assumed | unknown
  std::string faceConfirmed;     // ídem para la cara física
  bool groupConfirmed() const { return hvGroupConfirmed == "measurement"; }
  bool faceIsConfirmed() const { return faceConfirmed == "measurement"; }
};

struct Mapping {
  std::string version;
  std::string sha256;  // del archivo leído (cabecera del run de daqrun)
  std::string path;
  std::vector<MappingEntry> channels;
  std::map<int, int> groupToHvCh;  // grupo HV -> canal del SY5527

  static Result<Mapping> fromFile(const std::string& path);
  std::optional<int> groupOfHvChannel(int hvCh) const;
  std::vector<MappingEntry> forGroup(int group) const;
  int unknownGroupCount(const std::string& station) const;
};

}  // namespace daq::hvmon
