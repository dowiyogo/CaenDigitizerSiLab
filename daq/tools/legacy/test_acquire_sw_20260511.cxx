// test_acquire_sw_20260511.cxx — PROGRAMA HEREDADO, solo como referencia. NO se compila ni se
// ejecuta (no está en ningún CMakeLists ni Makefile; abre el digitizador).
//
// Origen: blob colgante 31394dbdc72a803e98071620766a3b4ac7c9b1f0 del repositorio (ruta
// test_acquire_sw.cxx dentro del árbol colgante eccefd4ead41bd8b0fe493ade20467a28a645183),
// encontrado en la tarea DAQ2-3 y restaurado por René el 2026-10-19. Es la versión de
// test_acquire_sw.cxx que produjo logs/commissioning_2026-05/acquire_sw.log del 2026-05-11
// ("Starting software-trigger acquisition...", máscara de fuentes 0x810C = 0xC0000000, 10 eventos,
// "<rate>: 294.10 Hz", 1023 muestras).
//
// Trigger por software: readEventsTree(nevents, true, ...) con automatic = true, que hace que la
// biblioteca llame CAEN_DGTZ_SendSWtrigger en cada vuelta del bucle de lectura
// (CaenDigitizerSiLab.cxx, readEventsTree). Sirve de referencia para el modo de calibración de
// daqrun (docs/daq_v2_design.md, sección 10). El código siguiente es el del blob, sin cambios.
#include "CaenDigitizerSiLab.h"
#include "configuration.h"
#include <iostream>

int quit = 0;

int main(int argc, char *argv[])
{
    const char* cfg = "config_vx1740_swtrigger.conf";
    if (argc > 1) cfg = argv[1];

    std::cout << "Reading config file: " << cfg << std::endl;
    readconfig(cfg);

    std::cout << "CONFIG SUMMARY\n"
              << "model = " << model << "\n"
              << "connectionType = " << connectionType << "\n"
              << "linkNum = " << linkNum << "\n"
              << "conetNode = " << conetNode << "\n"
              << "vmeBaseAddress = 0x" << std::hex << vmeBaseAddress << std::dec << "\n"
              << "acqsamples = " << acqsamples << "\n"
              << "nevents = " << nevents << "\n"
              << "timeout = " << timeout << "\n";

    CaenDigitizerSiLab *dig = new CaenDigitizerSiLab();
    dig->setConnection(connectionType, linkNum, conetNode, vmeBaseAddress);
    dig->setModel(model);
    dig->setTriggerPolarity(triggerpolaritymode);

    std::cout << "Initializing digitizer..." << std::endl;
    int ret = dig->init();
    if (ret < 0) {
        std::cout << "Init failed" << std::endl;
        delete dig;
        return 1;
    }

    dig->getInfo();
    dig->setPolarizationType(polarization);
    dig->setNSamples(acqsamples, ptriggersize);
    dig->setTrigmV(vthreshold);
    dig->newFileTree("test_vx1740_sw.root");

    std::cout << "Starting software-trigger acquisition..." << std::endl;
    dig->readEventsTree(nevents, true, 0, timeout, triggerSource);

    dig->storeDataTree();
    dig->closeLastFile();
    dig->finish();
    delete dig;

    std::cout << "Done test_acquire_sw" << std::endl;
    return 0;
}
