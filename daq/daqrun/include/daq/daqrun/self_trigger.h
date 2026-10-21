// self_trigger.h
// Umbrales del auto-trigger por grupos (DAQ2-13): de p.e. (trigger.self.threshold_pe) a cuentas ADC
// del registro 0x1n80, que es absoluto y común a los 8 canales del grupo.
//  - línea base de cada canal en ADC = pedestal de la calibración (mV) pasado a ADC con el DC offset
//    del grupo (AdcConversion, la fórmula del DAQ antiguo);
//  - por grupo: mediana de las líneas base y de G (mV/p.e.) de los canales que participan (sin G:
//    mediana de los canales de la estación con G); umbral = round(base − T·G/LSB);
//  - umbral efectivo por canal en p.e. = (base del canal − umbral)·LSB/G del canal: difiere entre
//    canales porque las líneas base del grupo no son iguales (hasta ~10 mV en c1). Va a la cabecera.
// Puro (sin hardware): lo usa daqrun antes de configurar el digitizador y lo prueban los tests.
#pragma once

#include "daq/core/result.h"
#include "daq/daqrun/calibration.h"
#include "daq/dig/digitizer.h"

namespace daq::daqrun {

Status applySelfTriggerThresholds(dig::DigitizerConfig& cfg, const PeCalibration& cal);

}  // namespace daq::daqrun
