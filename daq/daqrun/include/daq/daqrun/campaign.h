// campaign.h
// Campañas y runs (diseño 4.1, 4.4). Una campaña es el directorio <raw_root>/<nombre>/:
//  - daqcampaign open <nombre>: crea el directorio, el marcador CAMPAIGN_OPEN (JSON: nombre,
//    apertura UTC, operador, host, simulado) y el manifiesto <manifests_dir>/<nombre>.json.
//  - daqcampaign close <nombre>: rechaza si hay un run corriendo (lock); el marcador pasa a
//    CAMPAIGN_CLOSED (rename, nunca se borra) con el cierre UTC, la lista de runs y el SHA256 de
//    SHA256SUMS; manifiesto con status "closed"; el directorio queda sin permiso de escritura.
//  - Un run exige CAMPAIGN_OPEN y que no exista CAMPAIGN_CLOSED; toma el lock .run.lock del
//    directorio (un run a la vez por campaña) y el número siguiente al mayor existente.
// Nombres: <campaña>_<YYYYMMDDThhmmssZ>_run<NNN>_seg<MMM>.root; la campaña en el nombre los hace
// únicos entre campañas, y además se comprueba que ningún directorio de raw_root tenga ya ese
// nombre (hodo::resolveRawFile no debe encontrar dos).
// SHA256SUMS: "<sha256>  <archivo>\n", agregado con flock. Manifiesto: escritura atómica con flock.
// Seguridad entre hilos: funciones sin estado; la exclusión entre procesos es por flock.
#pragma once

#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "daq/core/instance_lock.h"
#include "daq/core/result.h"
#include "daq/core/time.h"

namespace daq::daqrun {

struct CampaignLocation {
  std::string rawRoot;
  std::string manifestsDir;
  std::string name;
  std::string dir() const { return rawRoot + "/" + name; }
  std::string manifestPath() const { return manifestsDir + "/" + name + ".json"; }
};

enum class CampaignState { Missing, Open, Closed, Invalid };
CampaignState campaignState(const CampaignLocation& loc);

Status openCampaign(const CampaignLocation& loc, const std::string& operatorName, bool simulated, UtcNs now);
Status closeCampaign(const CampaignLocation& loc, UtcNs now);
nlohmann::json campaignStatus(const CampaignLocation& loc);

// Reserva de un run: lock del directorio de la campaña + número + nombre base.
class RunReservation {
 public:
  static Result<RunReservation> reserve(const CampaignLocation& loc, UtcNs startUtc);
  int run() const { return run_; }
  // <campaña>_<YYYYMMDDThhmmssZ>_run<NNN>
  const std::string& baseName() const { return base_; }
  std::string segmentFile(int seg) const;  // nombre (sin directorio) del segmento seg (1, 2, ...)
  std::string indexFile() const { return base_ + ".index.json"; }
  const CampaignLocation& location() const { return loc_; }

 private:
  CampaignLocation loc_;
  InstanceLock lock_;
  int run_ = 0;
  std::string base_;
};

// true si algún directorio directo de rawRoot (o rawRoot mismo) ya tiene un archivo fileName.
bool fileNameUsedInAnyCampaign(const std::string& rawRoot, const std::string& fileName);

Status appendSha256Sum(const std::string& campaignDir, const std::string& fileName, const std::string& sha256);
// Agrega o reemplaza la entrada del run en el manifiesto (por número de run).
Status updateManifestRun(const CampaignLocation& loc, const nlohmann::json& runEntry);

}  // namespace daq::daqrun
