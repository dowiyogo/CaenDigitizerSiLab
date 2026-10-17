#include <iostream>
#include "CAENComm.h"

int main()
{
    const CAENComm_ConnectionType connectionType = CAENComm_USB_A4818;
    const int links[] = {0, 1, 62578};
    const int nodes[] = {0, 1, 2, 3};

    for (int linkNum : links) {
        for (int conetNode : nodes) {
            int handle = -1;
            char errorString[256] = {0};

            std::cout << "===== CAENComm_USB_A4818 linkNum=" << linkNum
                      << " conetNode=" << conetNode << " =====" << std::endl;

            CAENComm_ErrorCode ret = CAENComm_OpenDevice(connectionType,
                                                         linkNum,
                                                         conetNode,
                                                         0,
                                                         &handle);
            CAENComm_DecodeError(ret, errorString);

            std::cout << "ret = " << ret << std::endl;
            std::cout << "CAENComm error string: " << errorString << std::endl;

            if (ret == CAENComm_Success) {
                std::cout << "CAENComm open OK. handle=" << handle << std::endl;
                CAENComm_CloseDevice(handle);
            }
        }
    }

    return 0;
}
