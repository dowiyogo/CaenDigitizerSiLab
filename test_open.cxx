#include "CaenDigitizerSiLab.h"
#include "configuration.h"
#include <iostream>

int quit = 0; // define here to satisfy references in library

int main(int argc, char *argv[])
{
    const char* cfg = "config.conf";
    if (argc > 1) cfg = argv[1];

    std::cout << "Reading config file: " << cfg << std::endl;
    if(readconfig(cfg)!=0){
        std::cout<<"Using default configuration\n";
    }

    std::cout << "CONFIG SUMMARY\n"
              << "model = " << model << "\n"
              << "connectionType = " << connectionType << "\n"
              << "linkNum = " << linkNum << "\n"
              << "conetNode = " << conetNode << "\n"
              << "vmeBaseAddress = 0x" << std::hex << vmeBaseAddress << std::dec << "\n";

    CaenDigitizerSiLab *dig = new CaenDigitizerSiLab();
    dig->setConnection(connectionType, linkNum, conetNode, vmeBaseAddress);
    dig->setModel(model);
    dig->setTriggerPolarity(triggerpolaritymode);

    std::cout<<"Initializing digitizer...\n";
    int ret = dig->init();
    if(ret<0){
        std::cout<<"Init failed\n";
        delete dig;
        return -1;
    }

    dig->getInfo();
    dig->finish();
    delete dig;
    std::cout<<"Done test_open\n";
    return 0;
}
