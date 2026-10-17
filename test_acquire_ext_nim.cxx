#include "CaenDigitizerSiLab.h"
#include "configuration.h"
#include <csignal>
#include <iostream>

int quit = 0;

void signalHandler(int)
{
    quit = 1;
}

int main(int argc, char *argv[])
{
    const char* cfg = "config_vx1740_ext_nim.conf";
    const char* outputFile = "test_ext_nim.root";

    if (argc > 1) cfg = argv[1];
    if (argc > 2) outputFile = argv[2];

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
              << "timeout = " << timeout << "\n"
              << "triggerSource = " << triggerSource << "\n";

    CaenDigitizerSiLab *dig = new CaenDigitizerSiLab();

    dig->setConnection(connectionType, linkNum, conetNode, vmeBaseAddress);
    dig->setModel(model);
    dig->setTriggerPolarity(triggerpolaritymode);

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

    dig->newFileTree(outputFile);

    std::cout << "Starting VX1740 external NIM trigger acquisition..." << std::endl;
    std::cout << "Waiting for pulses on TRG-IN NIM..." << std::endl;
    std::cout << "Press Ctrl+C once for a clean stop." << std::endl;

    std::signal(SIGINT, signalHandler);
    std::signal(SIGTERM, signalHandler);
    dig->readEventsTree(nevents, false, 0, timeout, triggerSource);

    if (quit)
        std::cout << "\nStop requested; writing and closing ROOT file..." << std::endl;

    dig->storeDataTree();
    dig->closeLastFile();
    dig->finish();
    delete dig;

    std::cout << "Done test_acquire_ext_nim" << std::endl;
    return 0;
}
