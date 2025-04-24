#include "CaenDigitizerSiLab.h"

int main()
{
  CaenDigitizerSiLab *dig = new CaenDigitizerSiLab();
  dig->getInfo();
  dig->readTempAll();
  dig->init();
  return 0;
}
