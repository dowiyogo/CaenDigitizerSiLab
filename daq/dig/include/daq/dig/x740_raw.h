// x740_raw.h
// Partes del backend CAEN que no necesitan la placa (tests sin hardware):
//  - caenErrorName(): código CAEN_DGTZ_ErrorCode -> nombre y texto (CAENDigitizer 2.17.3 no trae
//    CAEN_DGTZ_DecodeError; la tabla copia los comentarios de CAENDigitizerType.h).
//  - Decodificación del formato de evento del x740 en el buffer de ReadData: encabezado de 4
//    palabras (0xA en los bits 31-28 y tamaño en palabras, id de placa, patrón y máscara de grupos,
//    contador de eventos de 24 bit, TTT de 32 bit) y separación de los eventos del bloque.
//  - decodeX740Samples(): muestras de 12 bit de los grupos (empaquetado validado contra
//    CAEN_DGTZ_DecodeEvent con buffers capturados en DAQ2-6).
// Funciones puras: seguras entre hilos.
#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "daq/core/result.h"

namespace daq::dig {

std::string caenErrorName(int code);

struct X740Header {
  std::uint32_t sizeWords = 0;  // tamaño del evento en palabras de 32 bit (incluye el encabezado)
  std::uint32_t boardId = 0;
  bool boardFail = false;
  std::uint32_t pattern = 0;
  std::uint32_t groupMask = 0;
  std::uint32_t counter = 0;  // 24 bit
  std::uint32_t ttt = 0;      // tal como llega (32 bit)
};

Result<X740Header> decodeX740Header(const std::uint32_t* w, std::size_t nWords);

struct X740EventSpan {
  std::size_t offsetWords = 0;
  X740Header header;
};
// Recorre un buffer de ReadData (bytes) y devuelve los eventos por su tamaño de encabezado.
Result<std::vector<X740EventSpan>> splitX740Buffer(const void* buffer, std::size_t bytes);

// Muestras de un evento x740: out[c * nSamples + k] para los 64 canales (canales de grupos no
// habilitados quedan en 0). nSamples se deduce del tamaño del evento.
Result<int> decodeX740Samples(const std::uint32_t* event, std::size_t nWords, std::vector<std::uint16_t>& out);

}  // namespace daq::dig
