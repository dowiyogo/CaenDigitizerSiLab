// x740_raw.cpp
#include "daq/dig/x740_raw.h"

#include <cstring>
#include <map>

namespace daq::dig {

std::string caenErrorName(int code) {
  static const std::map<int, const char*> names{
      {0, "CAEN_DGTZ_Success: operación completada"},
      {-1, "CAEN_DGTZ_CommError: error de comunicación"},
      {-2, "CAEN_DGTZ_GenericError: error no especificado"},
      {-3, "CAEN_DGTZ_InvalidParam: parámetro inválido"},
      {-4, "CAEN_DGTZ_InvalidLinkType: tipo de enlace inválido"},
      {-5, "CAEN_DGTZ_InvalidHandle: handle inválido"},
      {-6, "CAEN_DGTZ_MaxDevicesError: demasiados dispositivos"},
      {-7, "CAEN_DGTZ_BadBoardType: operación no permitida en este tipo de placa"},
      {-8, "CAEN_DGTZ_BadInterruptLev: nivel de interrupción no permitido"},
      {-9, "CAEN_DGTZ_BadEventNumber: número de evento incorrecto"},
      {-10, "CAEN_DGTZ_ReadDeviceRegisterFail: no se pudo leer el registro"},
      {-11, "CAEN_DGTZ_WriteDeviceRegisterFail: no se pudo escribir el registro"},
      {-13, "CAEN_DGTZ_InvalidChannelNumber: canal inválido"},
      {-14, "CAEN_DGTZ_ChannelBusy: canal ocupado"},
      {-15, "CAEN_DGTZ_FPIOModeInvalid: modo FPIO inválido"},
      {-16, "CAEN_DGTZ_WrongAcqMode: modo de adquisición incorrecto"},
      {-17, "CAEN_DGTZ_FunctionNotAllowed: función no permitida en este módulo"},
      {-18, "CAEN_DGTZ_Timeout: timeout de comunicación"},
      {-19, "CAEN_DGTZ_InvalidBuffer: buffer inválido"},
      {-20, "CAEN_DGTZ_EventNotFound: evento no encontrado"},
      {-21, "CAEN_DGTZ_InvalidEvent: evento inválido"},
      {-22, "CAEN_DGTZ_OutOfMemory: sin memoria"},
      {-23, "CAEN_DGTZ_CalibrationError: no se pudo calibrar la placa"},
      {-24, "CAEN_DGTZ_DigitizerNotFound: no se pudo abrir el digitizador"},
      {-25, "CAEN_DGTZ_DigitizerAlreadyOpen: el digitizador ya está abierto"},
      {-26, "CAEN_DGTZ_DigitizerNotReady: el digitizador no está listo"},
      {-27, "CAEN_DGTZ_InterruptNotConfigured: IRQ no configurada"},
      {-28, "CAEN_DGTZ_DigitizerMemoryCorrupted: memoria flash corrupta"},
      {-29, "CAEN_DGTZ_DPPFirmwareNotSupported: firmware DPP no soportado"},
      {-30, "CAEN_DGTZ_InvalidLicense: licencia de firmware inválida"},
      {-31, "CAEN_DGTZ_InvalidDigitizerStatus: estado del digitizador corrupto"},
      {-32, "CAEN_DGTZ_UnsupportedTrace: traza no soportada"},
      {-33, "CAEN_DGTZ_InvalidProbe: sonda no soportada"},
      {-34, "CAEN_DGTZ_UnsupportedBaseAddress: dirección base no soportada"},
      {-99, "CAEN_DGTZ_NotYetImplemented: función no implementada"}};
  auto it = names.find(code);
  return "código " + std::to_string(code) + " (" + (it == names.end() ? "desconocido" : it->second) + ")";
}

Result<X740Header> decodeX740Header(const std::uint32_t* w, std::size_t nWords) {
  using R = Result<X740Header>;
  if (nWords < 4) return R::failure(-1, "evento x740: menos de 4 palabras");
  if ((w[0] >> 28) != 0xA) return R::failure(-1, "evento x740: palabra 0 sin la marca 0xA");
  X740Header h;
  h.sizeWords = w[0] & 0x0FFFFFFF;
  h.boardId = (w[1] >> 27) & 0x1F;
  h.boardFail = ((w[1] >> 26) & 1) != 0;
  h.pattern = (w[1] >> 8) & 0xFFFF;
  h.groupMask = w[1] & 0xFF;
  h.counter = w[2] & 0xFFFFFF;
  h.ttt = w[3];
  if (h.sizeWords < 4 || h.sizeWords > nWords) return R::failure(-1, "evento x740: tamaño fuera del buffer");
  return h;
}

Result<std::vector<X740EventSpan>> splitX740Buffer(const void* buffer, std::size_t bytes) {
  using R = Result<std::vector<X740EventSpan>>;
  std::vector<X740EventSpan> out;
  const std::size_t nWords = bytes / 4;
  std::vector<std::uint32_t> w(nWords);
  if (nWords) std::memcpy(w.data(), buffer, nWords * 4);
  std::size_t i = 0;
  while (i < nWords) {
    auto h = decodeX740Header(w.data() + i, nWords - i);
    if (!h) return R::failure(-1, h.error().message + " (palabra " + std::to_string(i) + ")");
    out.push_back(X740EventSpan{i, h.value()});
    i += h.value().sizeWords;
  }
  return out;
}

Result<int> decodeX740Samples(const std::uint32_t* event, std::size_t nWords, std::vector<std::uint16_t>& out) {
  auto h = decodeX740Header(event, nWords);
  if (!h) return Result<int>(h.error());
  int groups = 0;
  for (int g = 0; g < 8; ++g) groups += (h.value().groupMask >> g) & 1;
  if (groups == 0) return Result<int>::failure(-1, "evento x740 sin grupos");
  const std::size_t dataWords = h.value().sizeWords - 4;
  if (dataWords % static_cast<std::size_t>(groups) != 0) return Result<int>::failure(-1, "evento x740: datos no divisibles por grupo");
  const std::size_t perGroup = dataWords / static_cast<std::size_t>(groups);
  // Cada grupo: bloques de 9 palabras = 3 muestras consecutivas de los 8 canales (8 × 3 × 12 bit
  // = 288 bit). Dentro del bloque, los 24 valores de 12 bit van en orden canal 0 (3 muestras),
  // canal 1, ... como un flujo de bits continuo desde el bit 0 de la primera palabra.
  if (perGroup % 9 != 0) return Result<int>::failure(-1, "evento x740: datos de grupo no múltiplo de 9 palabras");
  const int n = static_cast<int>(perGroup / 9 * 3);
  out.assign(static_cast<std::size_t>(64 * n), 0);
  const std::uint32_t* p = event + 4;
  for (int g = 0; g < 8; ++g) {
    if (!((h.value().groupMask >> g) & 1)) continue;
    for (std::size_t blk = 0; blk < perGroup / 9; ++blk) {
      const std::uint32_t* b = p + blk * 9;
      for (int v = 0; v < 24; ++v) {
        const int bit = v * 12;
        const int wi = bit / 32, off = bit % 32;
        std::uint64_t two = b[wi];
        if (wi + 1 < 9) two |= static_cast<std::uint64_t>(b[wi + 1]) << 32;
        const auto val = static_cast<std::uint16_t>((two >> off) & 0xFFF);
        const int ch = v / 3, s = static_cast<int>(blk) * 3 + v % 3;
        out[static_cast<std::size_t>((g * 8 + ch) * n + s)] = val;
      }
    }
    p += perGroup;
  }
  return n;
}

}  // namespace daq::dig
