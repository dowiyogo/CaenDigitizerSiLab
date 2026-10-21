// caen_digitizer.cpp
#include "daq/dig/caen_digitizer.h"

#include <CAENComm.h>
#include <CAENDigitizer.h>
#include <CAENVMElib.h>

#include <cstdio>
#include <cstring>
#include <fstream>

#include "daq/core/file_system.h"
#include "daq/core/logger.h"
#include "daq/dig/x740_raw.h"

namespace daq::dig {

namespace {
const char* trigModeName(CAEN_DGTZ_TriggerMode_t m) {
  switch (m) {
    case CAEN_DGTZ_TRGMODE_DISABLED: return "DISABLED";
    case CAEN_DGTZ_TRGMODE_EXTOUT_ONLY: return "EXTOUT_ONLY";
    case CAEN_DGTZ_TRGMODE_ACQ_ONLY: return "ACQ_ONLY";
    case CAEN_DGTZ_TRGMODE_ACQ_AND_EXTOUT: return "ACQ_AND_EXTOUT";
  }
  return "?";
}
std::string hex(std::uint32_t v) {
  char b[16];
  std::snprintf(b, sizeof b, "0x%08X", v);
  return b;
}
}  // namespace

CaenDigitizerBackend::CaenDigitizerBackend(Clock& clock) : clock_(clock) {}

CaenDigitizerBackend::~CaenDigitizerBackend() { close(); }

Status CaenDigitizerBackend::check(int code, const char* what) {
  if (code == CAEN_DGTZ_Success) return okStatus();
  lastError_ = std::string(what) + ": " + caenErrorName(code);
  logError("CAEN ", lastError_);
  return Status::failure(code, lastError_);
}

#define DAQ_CAEN(call) \
  if (auto s__ = check((call), #call); !s__) return s__

Status CaenDigitizerBackend::open(const DigitizerConfig& link) {
  cfg_ = link;
  char rel[64] = {0};
  CAEN_DGTZ_SWRelease(rel);
  info_["caen_digitizer_version"] = rel;
  char comm[64] = {0};
  CAENComm_SWRelease(comm);
  info_["caen_comm_version"] = comm;
  char vme[64] = {0};
  CAENVME_SWRelease(vme);
  info_["caen_vme_version"] = vme;
  info_["link"] = {{"type", link.linkType}, {"type_name", link.linkType == 5 ? "CAEN_DGTZ_USB_A4818" : "otro"},
                   {"link_num", link.linkNum}, {"conet_node", link.conetNode}, {"vme_base_address", hex(link.vmeBaseAddress)}};
  const int r = CAEN_DGTZ_OpenDigitizer(static_cast<CAEN_DGTZ_ConnectionType>(link.linkType), link.linkNum, link.conetNode,
                                        link.vmeBaseAddress, &handle_);
  DAQ_CAEN(r);
  open_ = true;
  CAEN_DGTZ_BoardInfo_t bi{};
  DAQ_CAEN(CAEN_DGTZ_GetInfo(handle_, &bi));
  info_["model"] = bi.ModelName;
  info_["model_code"] = bi.Model;
  info_["channels"] = bi.Channels;
  info_["form_factor"] = bi.FormFactor;
  info_["family_code"] = bi.FamilyCode;
  info_["roc_firmware"] = bi.ROC_FirmwareRel;
  info_["amc_firmware"] = bi.AMC_FirmwareRel;
  info_["serial"] = bi.SerialNumber;
  info_["pcb_revision"] = bi.PCB_Revision;
  info_["adc_nbits"] = bi.ADC_NBits;
  info_["comm_handle"] = bi.CommHandle;
  // A4818 (maestro CONET): número de serie y firmware por CAENComm_Info, si responde.
  char sn[64] = {0}, fw[64] = {0};
  const int e1 = CAENComm_Info(bi.CommHandle, CAENComm_PCI_Board_SN, sn);
  const int e2 = CAENComm_Info(bi.CommHandle, CAENComm_PCI_Board_FwRel, fw);
  info_["link"]["a4818_serial"] = e1 == 0 ? std::string(sn) : "error " + std::to_string(e1);
  info_["link"]["a4818_firmware"] = e2 == 0 ? std::string(fw) : "error " + std::to_string(e2);
  logInfo("digitizador abierto: ", bi.ModelName, " serie ", bi.SerialNumber, ", ROC ", bi.ROC_FirmwareRel, ", AMC ",
          bi.AMC_FirmwareRel, ", CAENDigitizer ", rel);
  return okStatus();
}

Result<std::uint32_t> CaenDigitizerBackend::readRegister(std::uint32_t addr) {
  std::uint32_t v = 0;
  const int r = CAEN_DGTZ_ReadRegister(handle_, addr, &v);
  if (r != CAEN_DGTZ_Success) return Result<std::uint32_t>::failure(r, "ReadRegister " + hex(addr) + ": " + caenErrorName(r));
  return v;
}

std::vector<std::pair<std::uint32_t, std::string>> CaenDigitizerBackend::standardRegisters() {
  std::vector<std::pair<std::uint32_t, std::string>> r{
      {0x8000, "group_configuration"}, {0x800C, "buffer_organization"},  {0x8020, "custom_size"},
      {0x8100, "acquisition_control"}, {0x8104, "acquisition_status"},   {0x810C, "global_trigger_mask"},
      {0x8110, "trigger_out_mask"},    {0x8114, "post_trigger"},         {0x811C, "front_panel_io_control"},
      {0x8120, "group_enable_mask"},   {0x8124, "roc_fpga_firmware"},    {0x812C, "event_stored"},
      {0x8140, "board_info"},          {0x814C, "event_size"},           {0x8178, "board_failure_status"},
      {0xEF00, "readout_control"},     {0xEF04, "readout_status"},       {0xEF08, "board_id"},
      {0xEF1C, "max_events_per_blt"}};
  for (std::uint32_t g = 0; g < 8; ++g) {
    const std::string p = "group" + std::to_string(g) + "_";
    r.push_back({0x1080 + 0x100 * g, p + "trigger_threshold"});
    r.push_back({0x1088 + 0x100 * g, p + "status"});
    r.push_back({0x108C + 0x100 * g, p + "amc_firmware"});
    r.push_back({0x1098 + 0x100 * g, p + "dc_offset"});
    r.push_back({0x10A8 + 0x100 * g, p + "channel_trigger_mask"});
  }
  return r;
}

nlohmann::json CaenDigitizerBackend::dumpRegisters(const std::vector<std::pair<std::uint32_t, std::string>>& regs) {
  nlohmann::json j = nlohmann::json::object();
  for (const auto& [addr, name] : regs) {
    auto v = readRegister(addr);
    j[hex(addr) + " " + name] = v ? hex(v.value()) : ("error " + v.error().message);
  }
  return j;
}

Status CaenDigitizerBackend::configure(const DigitizerConfig& cfg) {
  if (!open_) {
    if (auto s = open(cfg); !s) return s;
  }
  cfg_ = cfg;
  nlohmann::json written, readback;
  if (cfg.resetOnConfigure) DAQ_CAEN(CAEN_DGTZ_Reset(handle_));
  written["reset"] = cfg.resetOnConfigure;
  DAQ_CAEN(CAEN_DGTZ_SetRecordLength(handle_, static_cast<std::uint32_t>(cfg.recordLength)));
  written["record_length"] = cfg.recordLength;
  DAQ_CAEN(CAEN_DGTZ_SetGroupEnableMask(handle_, cfg.groupEnableMask));
  written["group_enable_mask"] = hex(cfg.groupEnableMask);
  for (std::uint32_t g = 0; g < 8; ++g)
    DAQ_CAEN(CAEN_DGTZ_SetGroupDCOffset(handle_, g, static_cast<std::uint32_t>(cfg.groupDcOffset[g])));
  written["group_dc_offset"] = cfg.groupDcOffset;
  DAQ_CAEN(CAEN_DGTZ_SetGroupSelfTrigger(handle_, CAEN_DGTZ_TRGMODE_DISABLED, 0xFF));
  written["group_self_trigger"] = "DISABLED";
  const auto& st = cfg.self;
  if (st.enabled) {
    // Auto-trigger por grupos (DAQ2-13): umbral 0x1n80 en cuentas ADC absolutas, canales que
    // participan 0x1nA8, polaridad negativa (0x8000 bit 6: dispara bajo el umbral).
    if (static_cast<int>(st.thresholdAdc.size()) != kGroups) return Status::failure(-1, "auto-trigger sin umbrales en ADC (falta la conversión)");
    nlohmann::json thr = nlohmann::json::object(), msk = nlohmann::json::object();
    for (std::uint32_t g = 0; g < 8; ++g) {
      if (!((st.groupMask >> g) & 1u)) continue;
      DAQ_CAEN(CAEN_DGTZ_SetGroupTriggerThreshold(handle_, g, static_cast<std::uint32_t>(st.thresholdAdc[g])));
      DAQ_CAEN(CAEN_DGTZ_SetChannelGroupMask(handle_, g, st.channelMask(static_cast<int>(g))));
      thr[std::to_string(g)] = st.thresholdAdc[g];
      msk[std::to_string(g)] = hex(st.channelMask(static_cast<int>(g)));
    }
    DAQ_CAEN(CAEN_DGTZ_SetGroupSelfTrigger(handle_, CAEN_DGTZ_TRGMODE_ACQ_ONLY, st.groupMask));
    DAQ_CAEN(CAEN_DGTZ_SetTriggerPolarity(handle_, 0, CAEN_DGTZ_TriggerOnFallingEdge));
    written["group_self_trigger"] = "ACQ_ONLY " + hex(st.groupMask);
    written["group_trigger_threshold_adc"] = thr;
    written["channel_group_mask"] = msk;
    written["trigger_polarity"] = "falling (bajo el umbral)";
  }
  const auto ext = cfg.externalTrigger ? CAEN_DGTZ_TRGMODE_ACQ_ONLY : CAEN_DGTZ_TRGMODE_DISABLED;
  DAQ_CAEN(CAEN_DGTZ_SetExtTriggerInputMode(handle_, ext));
  written["ext_trigger"] = trigModeName(ext);
  const auto io = cfg.ioLevel == "TTL" ? CAEN_DGTZ_IOLevel_TTL : CAEN_DGTZ_IOLevel_NIM;
  DAQ_CAEN(CAEN_DGTZ_SetIOLevel(handle_, io));
  written["io_level"] = cfg.ioLevel;
  const auto sw = cfg.softwareTrigger ? CAEN_DGTZ_TRGMODE_ACQ_ONLY : CAEN_DGTZ_TRGMODE_DISABLED;
  DAQ_CAEN(CAEN_DGTZ_SetSWTriggerMode(handle_, sw));
  written["sw_trigger"] = trigModeName(sw);
  if (st.enabled) {
    // Mayoría de grupos (0x810C[26:24] = nivel m, dispara con m + 1 grupos; [23:20] = ventana en
    // relojes de 8 ns). Se escribe después de los Set* que tocan 0x810C (grupos, TRG IN, SW).
    auto gm = readRegister(0x810C);
    if (!gm) return Status(gm.error());
    std::uint32_t v = gm.value() & ~((0x7u << 24) | (0xFu << 20));
    v |= (static_cast<std::uint32_t>(st.majorityLevel()) & 0x7u) << 24;
    v |= (static_cast<std::uint32_t>(st.windowClk) & 0xFu) << 20;
    DAQ_CAEN(CAEN_DGTZ_WriteRegister(handle_, 0x810C, v));
    written["global_trigger_mask_0x810C"] = hex(v);
    written["majority_level"] = st.majorityLevel();
    written["majority_window_clk"] = st.windowClk;
  }
  DAQ_CAEN(CAEN_DGTZ_SetAcquisitionMode(handle_, CAEN_DGTZ_SW_CONTROLLED));
  written["acquisition_mode"] = "SW_CONTROLLED";
  DAQ_CAEN(CAEN_DGTZ_SetMaxNumEventsBLT(handle_, static_cast<std::uint32_t>(cfg.maxEventsBlt)));
  written["max_events_blt"] = cfg.maxEventsBlt;
  // Post-trigger: como el DAQ antiguo, escribiendo 0x8114 (SetPostTriggerSize no sirve en el x740).
  DAQ_CAEN(CAEN_DGTZ_WriteRegister(handle_, 0x8114, static_cast<std::uint32_t>(cfg.postTriggerRegister)));
  written["post_trigger_register_0x8114"] = cfg.postTriggerRegister;
  // Relectura.
  std::uint32_t rl = 0;
  DAQ_CAEN(CAEN_DGTZ_GetRecordLength(handle_, &rl));
  recordLength_ = static_cast<int>(rl);
  readback["record_length"] = rl;
  std::uint32_t mask = 0;
  DAQ_CAEN(CAEN_DGTZ_GetGroupEnableMask(handle_, &mask));
  readback["group_enable_mask"] = hex(mask);
  std::vector<std::uint32_t> dc(8);
  for (std::uint32_t g = 0; g < 8; ++g) DAQ_CAEN(CAEN_DGTZ_GetGroupDCOffset(handle_, g, &dc[g]));
  readback["group_dc_offset"] = dc;
  CAEN_DGTZ_TriggerMode_t m{};
  DAQ_CAEN(CAEN_DGTZ_GetExtTriggerInputMode(handle_, &m));
  readback["ext_trigger"] = trigModeName(m);
  DAQ_CAEN(CAEN_DGTZ_GetSWTriggerMode(handle_, &m));
  readback["sw_trigger"] = trigModeName(m);
  CAEN_DGTZ_IOLevel_t lvl{};
  DAQ_CAEN(CAEN_DGTZ_GetIOLevel(handle_, &lvl));
  readback["io_level"] = lvl == CAEN_DGTZ_IOLevel_NIM ? "NIM" : "TTL";
  std::uint32_t blt = 0;
  DAQ_CAEN(CAEN_DGTZ_GetMaxNumEventsBLT(handle_, &blt));
  readback["max_events_blt"] = blt;
  CAEN_DGTZ_AcqMode_t am{};
  DAQ_CAEN(CAEN_DGTZ_GetAcquisitionMode(handle_, &am));
  readback["acquisition_mode"] = am == CAEN_DGTZ_SW_CONTROLLED ? "SW_CONTROLLED" : "otro";
  for (const auto& [addr, name] : std::vector<std::pair<std::uint32_t, const char*>>{
           {0x8114, "post_trigger_register_0x8114"}, {0x810C, "global_trigger_mask_0x810C"},
           {0x811C, "front_panel_io_0x811C"}, {0x8120, "group_enable_mask_0x8120"}, {0x8100, "acquisition_control_0x8100"},
           {0x800C, "buffer_organization_0x800C"}, {0x8020, "custom_size_0x8020"}, {0x8000, "group_configuration_0x8000"}}) {
    auto v = readRegister(addr);
    if (!v) return Status(v.error());
    readback[name] = addr == 0x8114 ? nlohmann::json(v.value()) : nlohmann::json(hex(v.value()));
  }
  if (st.enabled) {
    nlohmann::json thr = nlohmann::json::object(), msk = nlohmann::json::object(), mode = nlohmann::json::object(), reg = nlohmann::json::object();
    for (std::uint32_t g = 0; g < 8; ++g) {
      if (!((st.groupMask >> g) & 1u)) continue;
      std::uint32_t tv = 0, cm = 0;
      CAEN_DGTZ_TriggerMode_t gmode{};
      DAQ_CAEN(CAEN_DGTZ_GetGroupTriggerThreshold(handle_, g, &tv));
      DAQ_CAEN(CAEN_DGTZ_GetChannelGroupMask(handle_, g, &cm));
      DAQ_CAEN(CAEN_DGTZ_GetGroupSelfTrigger(handle_, g, &gmode));
      thr[std::to_string(g)] = tv;
      msk[std::to_string(g)] = hex(cm);
      mode[std::to_string(g)] = trigModeName(gmode);
      for (std::uint32_t a : {0x1080u + 0x100u * g, 0x10A8u + 0x100u * g}) {
        auto r = readRegister(a);
        if (!r) return Status(r.error());
        reg[hex(a)] = hex(r.value());
      }
    }
    CAEN_DGTZ_TriggerPolarity_t pol{};
    DAQ_CAEN(CAEN_DGTZ_GetTriggerPolarity(handle_, 0, &pol));
    readback["group_trigger_threshold_adc"] = thr;
    readback["channel_group_mask"] = msk;
    readback["group_self_trigger"] = mode;
    readback["trigger_polarity"] = pol == CAEN_DGTZ_TriggerOnFallingEdge ? "falling" : "rising";
    readback["group_registers"] = reg;
    auto gm = readRegister(0x810C);
    if (!gm) return Status(gm.error());
    readback["majority_level"] = (gm.value() >> 24) & 0x7u;
    readback["majority_window_clk"] = (gm.value() >> 20) & 0xFu;
    for (std::uint32_t g = 0; g < 8; ++g)
      if (((st.groupMask >> g) & 1u) && thr[std::to_string(g)].get<std::uint32_t>() != static_cast<std::uint32_t>(st.thresholdAdc[g]))
        return Status::failure(-1, "auto-trigger: el umbral releído del grupo " + std::to_string(g) + " no coincide con el escrito");
    if (((gm.value() >> 24) & 0x7u) != static_cast<std::uint32_t>(st.majorityLevel()) || ((gm.value() >> 20) & 0xFu) != static_cast<std::uint32_t>(st.windowClk))
      return Status::failure(-1, "auto-trigger: la mayoría releída de 0x810C no coincide con la escrita");
  }
  info_["config_written"] = written;
  info_["config_read"] = readback;
  if (recordLength_ != cfg.recordLength)
    logWarn("longitud de registro: pedida ", cfg.recordLength, ", la placa usa ", recordLength_);
  // Buffers de lectura y de evento decodificado.
  if (buffer_) CAEN_DGTZ_FreeReadoutBuffer(&buffer_);
  DAQ_CAEN(CAEN_DGTZ_MallocReadoutBuffer(handle_, &buffer_, &bufferSize_));
  info_["readout_buffer_bytes"] = bufferSize_;
  if (!event_) DAQ_CAEN(CAEN_DGTZ_AllocateEvent(handle_, &event_));
  if (!cfg.captureDir.empty()) PosixFileSystem().makeDirectories(cfg.captureDir, 0755);
  return okStatus();
}

Result<AcqOrigin> CaenDigitizerBackend::start() {
  AcqOrigin o;
  o.utcBeforeStart = clock_.now_utc_ns();
  const int r = CAEN_DGTZ_SWStartAcquisition(handle_);
  o.utcAfterStart = clock_.now_utc_ns();
  if (r != CAEN_DGTZ_Success) return Result<AcqOrigin>::failure(r, "SWStartAcquisition: " + caenErrorName(r));
  running_ = true;
  return o;
}

Status CaenDigitizerBackend::stop() {
  if (!running_) return okStatus();
  running_ = false;
  return check(CAEN_DGTZ_SWStopAcquisition(handle_), "CAEN_DGTZ_SWStopAcquisition");
}

Status CaenDigitizerBackend::sendSoftwareTrigger() { return check(CAEN_DGTZ_SendSWtrigger(handle_), "CAEN_DGTZ_SendSWtrigger"); }

Result<std::size_t> CaenDigitizerBackend::readBlock(EventBlock& out, std::size_t maxEvents) {
  using R = Result<std::size_t>;
  (void)maxEvents;  // la placa entrega hasta max_events_blt por ReadData
  out.recordLength = recordLength_;
  // Estado de adquisición antes de leer: memoria llena desde la lectura anterior.
  std::uint32_t st = 0;
  if (CAEN_DGTZ_ReadRegister(handle_, 0x8104, &st) == CAEN_DGTZ_Success) {
    ++statusReads_;
    lastStatus_ = st;
    if (st & memFullMask_) {
      ++memFullReads_;
      out.boardFull = true;
    }
  }
  std::uint32_t bsize = bufferSize_;
  int r = CAEN_DGTZ_ReadData(handle_, CAEN_DGTZ_SLAVE_TERMINATED_READOUT_MBLT, buffer_, &bsize);
  out.tReadUtcNs = clock_.now_utc_ns();
  out.tReadMonoNs = monotonicNowNs();
  ++readDataCalls_;
  if (r != CAEN_DGTZ_Success) return R::failure(r, "ReadData: " + caenErrorName(r));
  if (bsize == 0) return std::size_t{0};
  bytesRead_ += bsize;
  std::uint32_t nEv = 0;
  r = CAEN_DGTZ_GetNumEvents(handle_, buffer_, bsize, &nEv);
  if (r != CAEN_DGTZ_Success) return R::failure(r, "GetNumEvents: " + caenErrorName(r));
  nlohmann::json capture = nlohmann::json::array();
  const bool doCapture = !cfg_.captureDir.empty() && captured_ < cfg_.captureBlocks;
  for (std::uint32_t i = 0; i < nEv; ++i) {
    CAEN_DGTZ_EventInfo_t ei{};
    char* ptr = nullptr;
    r = CAEN_DGTZ_GetEventInfo(handle_, buffer_, bsize, static_cast<int32_t>(i), &ei, &ptr);
    if (r != CAEN_DGTZ_Success) return R::failure(r, "GetEventInfo: " + caenErrorName(r));
    r = CAEN_DGTZ_DecodeEvent(handle_, ptr, &event_);
    if (r != CAEN_DGTZ_Success) return R::failure(r, "DecodeEvent: " + caenErrorName(r));
    const auto* ev = static_cast<const CAEN_DGTZ_UINT16_EVENT_t*>(event_);
    Event e;
    e.counterHw = ei.EventCounter;
    e.tttRaw = ei.TriggerTimeTag;
    // 0 externo, 1 software, 2 auto-trigger, 3 auto-trigger o externo (no se distinguen en el evento)
    e.triggerKind = static_cast<std::uint8_t>(cfg_.self.enabled ? (cfg_.externalTrigger ? 3 : 2)
                                                                : (cfg_.softwareTrigger && !cfg_.externalTrigger ? 1 : 0));
    e.adc.assign(static_cast<std::size_t>(kChannels * recordLength_), 0);
    for (int c = 0; c < kChannels; ++c) {
      const int n = std::min<int>(static_cast<int>(ev->ChSize[c]), recordLength_);
      if (ev->DataChannel[c] && n > 0) std::memcpy(e.adc.data() + c * recordLength_, ev->DataChannel[c], static_cast<std::size_t>(n) * 2);
    }
    if (doCapture)
      capture.push_back({{"index", i}, {"event_size", ei.EventSize}, {"board_id", ei.BoardId}, {"pattern", ei.Pattern},
                         {"channel_mask", ei.ChannelMask}, {"counter", ei.EventCounter}, {"ttt", ei.TriggerTimeTag},
                         {"ch_size0", ev->ChSize[0]}, {"adc", e.adc}});
    out.events.push_back(std::move(e));
  }
  if (doCapture) {
    char base[64];
    std::snprintf(base, sizeof base, "/block_%03d", captured_);
    const std::string stem = cfg_.captureDir + base;
    std::ofstream(stem + ".bin", std::ios::binary).write(buffer_, bsize);
    std::ofstream(stem + ".json") << nlohmann::json{{"bytes", bsize}, {"record_length", recordLength_},
                                                    {"t_read_utc_ns", out.tReadUtcNs}, {"events", capture}}.dump() << "\n";
    ++captured_;
  }
  return static_cast<std::size_t>(nEv);
}

nlohmann::json CaenDigitizerBackend::info() const {
  nlohmann::json j = info_;
  j["record_length_actual"] = recordLength_;
  j["memory_full_mask_0x8104"] = hex(memFullMask_);
  j["acquisition_status_reads"] = statusReads_;
  j["memory_full_reads"] = memFullReads_;
  j["last_acquisition_status"] = hex(lastStatus_);
  j["read_data_calls"] = readDataCalls_;
  j["bytes_read"] = bytesRead_;
  j["last_error"] = lastError_;
  return j;
}

void CaenDigitizerBackend::close() {
  if (!open_) return;
  if (running_) stop();
  if (event_) CAEN_DGTZ_FreeEvent(handle_, &event_);
  event_ = nullptr;
  if (buffer_) CAEN_DGTZ_FreeReadoutBuffer(&buffer_);
  buffer_ = nullptr;
  check(CAEN_DGTZ_CloseDigitizer(handle_), "CAEN_DGTZ_CloseDigitizer");
  open_ = false;
  logInfo("digitizador cerrado");
}

}  // namespace daq::dig
