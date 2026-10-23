// run_info_dump.C — DAQ2-16 (solo lectura). Copia a un JSON las cabeceras de un segmento de run: run_info_json
// (incluye digitizer.config_written y digitizer.config_read: registros escritos y releídos por daqrun) y config_json.
// Salida: {"run_info": {...}, "config": {...}}.
// Uso: root -l -b -q 'daq/tools/biasscan/run_info_dump.C("<run>_seg001.root", "<salida.json>")'
#include <TFile.h>
#include <TNamed.h>

#include <cstdio>
#include <memory>
#include <string>

void run_info_dump(const char* segment, const char* out) {
  std::unique_ptr<TFile> f(TFile::Open(segment, "READ"));
  if (!f || f->IsZombie()) {
    printf("no se pudo abrir %s\n", segment);
    return;
  }
  auto* ri = f->Get<TNamed>("run_info_json");
  auto* cj = f->Get<TNamed>("config_json");
  FILE* o = fopen(out, "w");
  fprintf(o, "{\"run_info\": %s, \"config\": %s}\n", ri ? ri->GetTitle() : "null", cj ? cj->GetTitle() : "null");
  fclose(o);
}
