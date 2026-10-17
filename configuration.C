#include "configuration.h"
#include "CaenDigitizerSiLab.h"

int model = 5730;
int polarization = 0;
int vthreshold = 50; //in mV
int acqsamples = 100;
int nevents = 100;
int bunches = 1;
uint32_t ptriggersize = 60;
int triggerpolaritymode = 1;
int timeout = 10;
int triggerSource = 0;
uint32_t triggerGroupMask = 0xFF;

// New connection defaults
int connectionType = 0; // 0=USB direct, 1=OpticalLink/CONET, 2=A4818, 3/4/5=A4818 bridge variants
int linkNum = 0;        // link number or PID for A4818
int conetNode = 0;      // CONET node (board position)
uint32_t vmeBaseAddress = 0;

char buf[200];
//char filename[200]="predeterminado";


int ReadUntil(int fd, char tok){
        if(fd<2){
                perror("ERROR fd to null file\n");exit(1);
        }
        int retread=1;
        int i=0;
        char data;
//      printf("leyendo archivo\n");
        while(true){
                if((retread=read(fd,&data,1))==-1){
                        perror("ERROR reading input\n"); exit(1);
                }
                if(retread==0) {
                        buf[i]=0;
                        return (i>0) ? i : -1;
                }
                if(data==tok)
                        break;
                if(i>=199){
                        perror("ERROR buf overflow\n"); exit(1);
                }
                buf[i++]=data;
        }
        buf[i]=0;
        return i;
}

int readconfig(const char *file){
        int fdin;
        char *param;
        if((fdin=open(file,O_RDONLY))==-1){
                perror("ERROR reading configuration file"); exit(1);
        }

        while(ReadUntil(fdin,'\n')>=0){
                param=strtok(buf," \t\r");
                if(param && param[0]!='#'){
                        char *value=strtok(NULL," \t\r");
                        if(!value)
                                continue;
//                      std::cout << param <<std::endl;
                        if(0==strcmp(param,"model"))
                                model=atoi(value);
                        if(0==strcmp(param,"polarization"))
                                polarization=atoi(value);
                        if(0==strcmp(param,"vthreshold"))
                                vthreshold=atoi(value);
                        if(0==strcmp(param,"acqsamples"))
                                acqsamples=atoi(value);
                        if(0==strcmp(param,"nevents"))
                                nevents=atoi(value);
                        if(0==strcmp(param,"bunches"))
                                bunches=atoi(value);
                        if(0==strcmp(param,"ptriggersize"))
                                ptriggersize=strtoul(value, NULL, 0);
                        if(0==strcmp(param,"triggerpolaritymode"))
                                triggerpolaritymode=atoi(value);
                        if(0==strcmp(param,"timeout"))
                                timeout=atoi(value);
                        if(0==strcmp(param,"triggerSource"))
                                triggerSource=atoi(value);
                        if(0==strcmp(param,"triggerGroupMask"))
                                triggerGroupMask=strtoul(value, NULL, 0);
                        if(0==strcmp(param,"connectionType"))
                                connectionType=atoi(value);
                        if(0==strcmp(param,"linkNum"))
                                linkNum=atoi(value);
                        if(0==strcmp(param,"conetNode"))
                                conetNode=atoi(value);
                        if(0==strcmp(param,"vmeBaseAddress"))
                                vmeBaseAddress=strtoul(value, NULL, 0);
                }
        }
        close(fdin);
        return 0;
}
