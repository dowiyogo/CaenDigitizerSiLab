// http_page.cpp
#include "http_page.h"

#include <TGraph.h>
#include <TH1D.h>
#include <TH2F.h>
#include <THttpServer.h>
#include <TNamed.h>

#include <sstream>

namespace daq::daqrun {

Result<std::unique_ptr<DaqrunPage>> DaqrunPage::create(const DaqrunConfig& cfg, const std::string& bindAddress) {
  using R = Result<std::unique_ptr<DaqrunPage>>;
  std::unique_ptr<DaqrunPage> p(new DaqrunPage());
  const std::string engine = "http:" + std::to_string(cfg.httpPort) + "?bind=" + bindAddress + "&thrds=2";
  p->server_ = std::make_unique<THttpServer>(engine.c_str());
  if (!p->server_->IsAnyEngine()) return R::failure(-1, "THttpServer no pudo iniciar " + engine);
  p->server_->SetReadOnly(kTRUE);
  p->server_->SetTimer(0, kTRUE);  // solo atiende en processRequests() (hilo principal)
  p->url_ = "http://" + bindAddress + ":" + std::to_string(cfg.httpPort) + "/";
  p->dt_ = std::make_unique<TH1D>("dt", "#Deltat entre triggers;log_{10}(#Deltat/s);eventos", MonitorView::kDtBins,
                                  MonitorView::kDtLogMin, MonitorView::kDtLogMax);
  p->oot_ = std::make_unique<TH1D>("oot_rate",
                                   "Tasa de pulsos fuera de tiempo sobre el umbral PROVISIONAL (hod2018 20 mV, hod2019 3 mV);"
                                   "canal;tasa [Hz]",
                                   64, -0.5, 63.5);
  const double lsb = dig::AdcConversion{cfg.digitizer.nBit, cfg.digitizer.maxVpp, 0x7FFF}.lsbMv();
  const int nb = std::max(1, static_cast<int>(cfg.ampHistMaxMv / lsb));
  p->amp_ = std::make_unique<TH2F>("amp", "Amplitud de pulsos fuera de tiempo (mV, sin calibración);canal;amplitud [mV]", 64,
                                   -0.5, 63.5, nb, 0.0, nb * lsb);
  p->chFlags_ = std::make_unique<TH1D>("ch_flags", "Eventos con el bit en algún canal;bit de ch_flags;eventos", 16, -0.5, 15.5);
  p->evtFlags_ = std::make_unique<TH1D>("evt_flags", "Eventos con el bit;bit de evt_flags;eventos", 9, -0.5, 8.5);
  p->rate_ = std::make_unique<TGraph>();
  p->rate_->SetName("trigger_rate");
  p->rate_->SetTitle("Tasa de trigger por período;horas UTC desde el inicio;tasa [Hz]");
  for (TObject* o : std::initializer_list<TObject*>{p->dt_.get(), p->rate_.get()}) p->server_->Register("/trigger", o);
  p->server_->Register("/pulses", p->oot_.get());
  p->server_->Register("/pulses", p->amp_.get());
  p->server_->Register("/flags", p->chFlags_.get());
  p->server_->Register("/flags", p->evtFlags_.get());
  p->table_ = std::make_unique<TNamed>("table", "");
  p->json_ = std::make_unique<TNamed>("state_json", "{}");
  p->server_->Register("/status", p->table_.get());
  p->server_->Register("/status", p->json_.get());
  return p;
}

DaqrunPage::~DaqrunPage() {
  if (server_) {
    for (TObject* o : std::initializer_list<TObject*>{dt_.get(), oot_.get(), amp_.get(), chFlags_.get(), evtFlags_.get(),
                                                      rate_.get(), table_.get(), json_.get()})
      if (o) server_->Unregister(o);
  }
}

void DaqrunPage::update(const MonitorView& v, const HvSample& hv, const nlohmann::json& status) {
  for (int b = 0; b < MonitorView::kDtBins; ++b) dt_->SetBinContent(b + 1, v.dtHist[static_cast<std::size_t>(b)]);
  for (int c = 0; c < 64; ++c) {
    const auto ci = static_cast<std::size_t>(c);
    oot_->SetBinContent(c + 1, v.ootTimeS[ci] > 0 ? static_cast<double>(v.ootCounts[ci]) / v.ootTimeS[ci] : 0.0);
    for (int b = 0; b < v.ampBins; ++b) amp_->SetBinContent(c + 1, b + 1, v.amp[static_cast<std::size_t>(c * v.ampBins + b)]);
  }
  const auto w = status.value("writer", nlohmann::json::object());
  const auto chf = w.value("ch_flag_events", nlohmann::json::array());
  const auto evf = w.value("evt_flag_events", nlohmann::json::array());
  for (std::size_t b = 0; b < chf.size() && b < 16; ++b) chFlags_->SetBinContent(static_cast<int>(b) + 1, chf[b].get<double>());
  for (std::size_t b = 0; b < evf.size() && b < 9; ++b) evtFlags_->SetBinContent(static_cast<int>(b) + 1, evf[b].get<double>());
  rate_->Set(0);
  const UtcNs t0 = v.firstEventUtcNs;
  for (const auto& r : v.rateHistory) rate_->AddPoint(static_cast<double>(r.first - t0) / 3.6e12, r.second);

  std::ostringstream os;
  os << "daqrun " << status.value("run_name", "") << (status.value("simulated", false) ? "  [SIMULADO]" : "") << "\n";
  if (!hv.alive) os << "AVISO: hvmon no disponible (" << hv.error << "): bit 4 (HV_UNKNOWN) en todos los canales\n";
  os << "calibración: none (umbrales PROVISIONALES en mV)\n";
  os << "eventos: " << v.events << "   tasa último período: " << v.last.triggerRateHz << " Hz\n";
  const auto rd = status.value("readout", nlohmann::json::object());
  os << "colas: escritura " << status.value("writer_queue", 0) << "/" << status.value("writer_capacity", 0)
     << " (máx " << rd.value("writer_max_occupancy", 0) << ", bloqueos " << rd.value("writer_blocked", 0) << "), monitor "
     << status.value("monitor_queue", 0) << "/" << status.value("monitor_capacity", 0) << " (descartados "
     << rd.value("monitor_dropped", 0) << ")\n";
  os << "segmento: " << w.value("current_segment", "") << "   cerrados: " << w.value("segments_closed", 0) << "\n";
  table_->SetTitle(os.str().c_str());
  json_->SetTitle(status.dump().c_str());
}

void DaqrunPage::processRequests() {
  if (server_) server_->ProcessRequests();
}

}  // namespace daq::daqrun
