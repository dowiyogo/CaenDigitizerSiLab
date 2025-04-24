////////////////////////////////////////////////////////////////////////////////
////  Desarrollado por Orlando Soto.
////  Comentado, analizado y reparado por Jairo González.
////  Actualizado por Rene Rios.
////  
//// 
////
////  SiLab 2018 / X-PLUS 2025
////
////////////////////////////////////////////////////////////////////////////////


#include "CaenDigitizerSiLab.h"
#include "configuration.h"
#include "TBenchmark.h"
#include <signal.h>

CaenDigitizerSiLab *dig;
TBenchmark *bench;
char filename[50];

//Variables para control de flujo del programa (signal handler)
int quit=0;
int stage = 0;

//Funcion para cerrar correctamente el script con ctrl+c
void signalHandler(int sig)
{
  quit = 1;
  if (stage == 0){
    dig->storeData();
    dig->stopSWAcq();
    delete dig;
    bench->Show("example");
    delete bench;
    exit(sig);
  }
}

int main(int argc,char *argv[])
{
  struct sigaction act;
  act.sa_handler = signalHandler;
  sigemptyset(&act.sa_mask);
  act.sa_flags = 0;
  
  if (argc>2){		/*There are not arguments for input in this executable file*/
		printf("you just can insert only one parameter \n\
			example [<filename>]\n"); 
		return 0;
	}
	else if(argc>1){
		sprintf(filename,"%s.root",argv[1]);	
		}
	else {
		sprintf(filename,"data_from_digitizer.root");
	}

  printf("#############################################################\n"
         "#                                                           #\n"
         "#         C A E N    D I G I T I Z E R    X - P L U S       #\n"
         "#                          ULS                              #\n"
         "#                                                           #\n"
         "#############################################################\n\n"
         "Exit: Ctrl+c\n"
         "Force Exit: Ctrl+\\\n\n");

  //benchmark mide el desempeño del programa, no es vital
  bench = new TBenchmark();
  bench->Start("example");

  //Se lee el archivo config.conf para configurar el digitizer
  if(readconfig("config.conf"))
  	{
  		printf("Error al leer archivo de configuracion, tomando valores por defecto...");
  	}

  dig = new CaenDigitizerSiLab();

  dig->setModel(model);
  dig->setTriggerPolarity(triggerpolaritymode);
  printf("\nDigitizer initialization...\n"
          "----------------------------------------------------------\n");
  if (dig->init()>=0) {
    dig->getInfo();
    printf("----------------------------------------------------------\n\n");

    printf("Digitizer Configuration...\n"
          "----------------------------------------------------------\n");
    dig->setPolarizationType(polarization);//rango de polarizacion
    printf("polarization: %d\n",polarization);
    
    dig->setNSamples(acqsamples, ptriggersize);//samples por evento y post trigger size
    dig->setTrigmV(vthreshold);//threshold en milivolts

    int bunch_size=nevents;  //numero de eventos
    printf("Max Events: %d\n",nevents);
    int NBunch=bunches; //numero de grupos de eventos
    printf("----------------------------------------------------------\n\n");

    //Medición
    printf("Starting Data Adquisition...\n"
          "----------------------------------------------------------\n");
    dig->newFile(Form("%s",filename));
    sigaction(SIGINT, &act, 0); 
    for (int k=0;k<NBunch;k++)
    {
      std::cout<<"On bunch: "<<k<<std::endl;
      
      printf("start: triggerSource=%d\n",triggerSource);
      stage = 1;
      dig->readEvents(bunch_size,false,k*bunch_size,timeout,triggerSource); //lectura con selftrigger y timeout
      dig->storeData();
      dig->storeTempAll();
      printf("Done.\n----------------------------------------------------------\n\n");
      printf("Final Calibration...\n");
      dig->calibrate();

    }
    dig->closeLastFile();

    dig->finish();

    bench->Show("example");

    delete dig;
    delete bench;
  }
  else
  {
    bench->Show("example");
    delete dig;
    delete bench;
  }
  

  return 0;
}
