# CaenDigitizerSiLab
C++ shared Library to use caen digitizers coupled with root.

To correctly use this library, edit your .bashrc file writting this line:

  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"dir"

Where 'dir' means the path to this folder.

Diagnostic VX1740 self-trigger example for hod2019 only:

  ./test_acquire_sw config_vx1740_autotrigger_hod2019_ch32_63.conf outputs/hod2019_diag.root
