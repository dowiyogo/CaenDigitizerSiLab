#include "CaenDigitizerSiLab.h"
#include "configuration.h"
#include <iostream>

int quit = 0;

int main(int argc, char *argv[])
{
    const char* cfg = "config_vx1740_autotrigger_or.conf";
    const char* outfile = "test_vx1740_autotrigger_or.root";
    if (argc > 1) cfg = argv[1];
    if (argc > 2) outfile = argv[2];

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
    dig->newFileTree(outfile);

    std::cout << "Starting VX1740 self-trigger OR acquisition..." << std::endl;
    dig->readEventsTree(nevents, false, 0, timeout, triggerSource);

    dig->storeDataTree();
    dig->closeLastFile();
    dig->finish();
    delete dig;

    std::cout << "Done test_acquire_sw" << std::endl;
    return 0;
}
