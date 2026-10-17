#ifndef _CONFIGURATION_file_
#define _CONFIGURATION_file_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include "CaenDigitizerSiLab.h"


//extern int USB_HV;
//extern int USB_ard;
extern int model;
extern int polarization;
extern int vthreshold;
extern int acqsamples;
extern int nevents;
extern int bunches;
extern uint32_t ptriggersize;
extern int triggerpolaritymode;
extern int timeout;
extern int triggerSource;
extern uint32_t triggerGroupMask;

// New connection / addressing parameters
extern int connectionType; // 0=USB direct, 1=OpticalLink/CONET, 2=A4818, 3/4/5=A4818 bridge variants
extern int linkNum;        // link number or PID for A4818
extern int conetNode;      // CONET node (board position)
extern uint32_t vmeBaseAddress; // VME base address (hex allowed)

//extern char filename[200];


extern char buf[200];
int ReadUntil(int fd, char tok);
int readconfig(const char *file);

#endif
