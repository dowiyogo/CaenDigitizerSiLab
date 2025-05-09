#include "CaenDigitizerSiLab.h"
using namespace std::chrono;

ClassImp(CaenDigitizerSiLab)

int32_t CaenDigitizerSiLab::init()
{
  if (kTempSupported)
  {
    meanTemp= new float[MaxNCh];
    varTemp= new float[MaxNCh];
    storedMeanTemp= new float[MaxNCh];
    storedVarTemp= new float[MaxNCh];
    tempFile.open("temp_measurements.txt");
    kNSamplesTemp=30; // number of samples to determine mean an var. of channel temperature.
    kDtTemp=1;// samplinbg time for temperature measurements.
  }


  NCh=MaxNCh;

  switch (kPolarizationType)
  {
  case 0:
    kOffset=0x7FFF;//-1.52590218966963675e-05
    break;
  case -1:
    kOffset=0x1000; // (1 - -1)*0x1000/0xffff - 1 = -8.74998092622262913e-01
    break;
  case 1:
    kOffset=0xefff;//8.74998092622262913e-01
    break;
  }
  triggermode = CAEN_DGTZ_TRGMODE_ACQ_ONLY;
  Evt = NULL;
  acqmode = CAEN_DGTZ_SW_CONTROLLED;
  buffer = NULL;
  evtptr = NULL;
  ret = CAEN_DGTZ_OpenDigitizer(CAEN_DGTZ_USB,0,0,0,&handle);
  if(ret != CAEN_DGTZ_Success) {
    printf("Can't open digitizer\n");
    return -1;
  }
  ret = CAEN_DGTZ_Reset(handle);// resetting Digitizer
  if (kDoCalibration)
  {
    calibrate();
  }
  //std::cout<<"number of samnples per acq: "<<kSamples<<std::endl;
  ret = CAEN_DGTZ_SetRecordLength(handle,kSamples); // samples a grabar por acquisition windows CAMBIAR es decir tamaño de buffer
  
  if (kModel==5740){
    ret = CAEN_DGTZ_SetGroupEnableMask(handle,0xF);//Todos habilitados
  }
  else {
    ret = CAEN_DGTZ_SetChannelEnableMask(handle,kEnableMask);
  }


  trigthresh=mV2adc(10);
  //std::cout<<"trigthreshold: "<< adc2mV(trigthresh)<<std::endl;

  uint32_t dat=0;
  uint32_t reg = 0x1088;
  if (kModel==5740){
    for (int32_t k=0;k<MaxNCh/8;k++)
    {
        ret = CAEN_DGTZ_ReadRegister(handle,(reg+k*0x100),&dat);
        ret = CAEN_DGTZ_SetGroupDCOffset(handle,k,kOffset);
        ret = CAEN_DGTZ_SetGroupTriggerThreshold(handle,k,trigthresh);                  /* Set selfTrigger threshold */
        //ret = CAEN_DGTZ_SetTriggerPolarity(handle, k, CAEN_DGTZ_TriggerOnFallingEdge);
        ret = CAEN_DGTZ_SetTriggerPolarity(handle, k, kTriggerpolaritymode);
    }

  }
  else{
    for (int32_t k=0;k<MaxNCh;k++)
    {
        ret = CAEN_DGTZ_ReadRegister(handle,(reg+k*0x100),&dat);
        ret = CAEN_DGTZ_SetChannelDCOffset(handle,k,kOffset);
        ret = CAEN_DGTZ_SetChannelTriggerThreshold(handle,k,trigthresh);                  /* Set selfTrigger threshold */
        //ret = CAEN_DGTZ_SetTriggerPolarity(handle, k, CAEN_DGTZ_TriggerOnFallingEdge);
        ret = CAEN_DGTZ_SetTriggerPolarity(handle, k, kTriggerpolaritymode);
    }
  }
  


  ret = CAEN_DGTZ_WriteRegister(handle, 0x811C, 0x000001); // Set Trig-In Trig-Out as TTL (not NIM)
  ret = CAEN_DGTZ_SetSWTriggerMode(handle, triggermode);//modo trigger
  ret = CAEN_DGTZ_SetPostTriggerSize(handle,60);
  //% a grabar por cada trigger del record length
  ret = CAEN_DGTZ_SetAcquisitionMode(handle, acqmode);			// modo de adquisicion
  ret = CAEN_DGTZ_SetMaxNumEventsBLT(handle,1);//numero maximo de eventos por transferencia

  //CAEN_DGTZ_SetChannelPairTriggerLogic(handle,6,7,CAEN_DGTZ_LOGIC_AND,2);
  if(ret != CAEN_DGTZ_Success) {
    printf("Errors during Digitizer Configuration.\n");
  }

  //TESTING LINE
  ret = CAEN_DGTZ_MallocReadoutBuffer(handle,&buffer,(uint32_t *)&size);


  //localiza buffer pointer en la memoria
  //the *buffer MUST be initialized to NULL
  //the size in byte of the buffer allocated
  return 0;
}

int32_t CaenDigitizerSiLab::newFile(const char* name)
{
  //New File
  TString branches;
  branches.Append(Form("Ch%d",0));
  for (int32_t k =1;k<NCh;k++)
  {
    branches.Append(":");
    branches.Append(Form("Ch%d",k));
  }
  branches.Append(":time:event:timestamp");
  ofile = new TFile(name,"recreate");
  data = new TNtuple("data","amp (adc ch) and time (nsample) time_stamp (us)",branches.Data() );

  return 0;
}

int32_t CaenDigitizerSiLab::closeLastFile()
{
  ofile->Close();
  return 0;
}

int32_t CaenDigitizerSiLab::setTrigADC(int32_t trigthresh)
{
  std::cout<<"trigthreshold: "<<adc2mV(trigthresh)<<std::endl;
  uint32_t dat = 0;
  if(kModel==5740){
    for (int32_t k=0;k<MaxNCh/8;k++)
    {
        ret = CAEN_DGTZ_ReadRegister(handle,(0x1088+k*0x100),&dat);
        ret = CAEN_DGTZ_SetGroupDCOffset(handle,k,kOffset);
        ret = CAEN_DGTZ_SetGroupTriggerThreshold(handle,k,trigthresh);                  /* Set selfTrigger threshold */
        //ret = CAEN_DGTZ_SetTriggerPolarity(handle, k, CAEN_DGTZ_TriggerOnFallingEdge);
        ret = CAEN_DGTZ_SetTriggerPolarity(handle, k, kTriggerpolaritymode);
    }

  }
  else{
    for (int32_t k=0;k<MaxNCh;k++)
    {
        ret = CAEN_DGTZ_ReadRegister(handle,(0x1088+k*0x100),&dat);
        ret = CAEN_DGTZ_SetChannelDCOffset(handle,k,kOffset);
        ret = CAEN_DGTZ_SetChannelTriggerThreshold(handle,k,trigthresh);                  /* Set selfTrigger threshold */
        //ret = CAEN_DGTZ_SetTriggerPolarity(handle, k, CAEN_DGTZ_TriggerOnFallingEdge);
        ret = CAEN_DGTZ_SetTriggerPolarity(handle, k, kTriggerpolaritymode);
    }
  }
  return 0;
}


int32_t CaenDigitizerSiLab::setTrigmV(float mV)
{
  trigthresh=mV2adc(mV);
  setTrigADC(trigthresh);
  return 0;
}


int32_t CaenDigitizerSiLab::getInfo()
{
  ret = CAEN_DGTZ_GetInfo(handle, &BoardInfo);
  printf("\nConnected to CAEN Digitizer Model %s\n", BoardInfo.ModelName);
  printf("\tROC FPGA Release is %s\n", BoardInfo.ROC_FirmwareRel);
  printf("\tAMC FPGA Release is %s\n\n", BoardInfo.AMC_FirmwareRel);
  return 0;
}

int32_t  CaenDigitizerSiLab::readEvents(int32_t events,bool automatic,int32_t start_event, uint32_t triggerSource)
{
  //printf("int32_t  CaenDigitizerSiLab::readEvents(int32_t events,bool automatic,int32_t start_event, uint32_t triggerSource)\n");
  int32_t count=0;
  uint32_t dat=0;

  //Se libera el buffer y se vuelve a ubicar, por si han habido cambios en las configuraciones de la ventana de adquisicion
  ret = CAEN_DGTZ_FreeReadoutBuffer(&buffer);
  ret = CAEN_DGTZ_MallocReadoutBuffer(handle,&buffer,(uint32_t *)&size);  //ret = CAEN_DGTZ_MallocReadoutBuffer(handle,&buffer,(uint32_t *)&size);
  
  if (!automatic){
    //ret = CAEN_DGTZ_SetChannelSelfTrigger(handle,CAEN_DGTZ_TRGMODE_ACQ_ONLY,(3<<6)); //Adjacent channels paired.
    //ret = CAEN_DGTZ_SetChannelSelfTrigger(handle,CAEN_DGTZ_TRGMODE_ACQ_ONLY,(0x3f<<2)); //Adjacent channels paired.
    if(triggerSource == 9){
      /* 
      ret = CAEN_DGTZ_SetChannelSelfTrigger(handle,CAEN_DGTZ_TRGMODE_DISABLED,0x00);
      ret = CAEN_DGTZ_SetExtTriggerInputMode(handle, CAEN_DGTZ_TRGMODE_ACQ_ONLY);
      printf("\n");
      */
      ret = CAEN_DGTZ_SetChannelSelfTrigger(handle,CAEN_DGTZ_TRGMODE_DISABLED,0xFF);
      ret = CAEN_DGTZ_SetExtTriggerInputMode(handle, CAEN_DGTZ_TRGMODE_ACQ_ONLY);

      CAEN_DGTZ_IOLevel_t myLevel={};
      CAEN_DGTZ_GetIOLevel (handle, &myLevel);
      
      CAEN_DGTZ_SetIOLevel (handle, CAEN_DGTZ_IOLevel_NIM);
      
      CAEN_DGTZ_GetIOLevel (handle, &myLevel);
    }
    else if ((0<=triggerSource)&&(triggerSource<=7)&&(kModel!=5740)){
      ret = CAEN_DGTZ_SetExtTriggerInputMode(handle, CAEN_DGTZ_TRGMODE_DISABLED);
      ret = CAEN_DGTZ_SetChannelSelfTrigger(handle,CAEN_DGTZ_TRGMODE_DISABLED,0xFF);
      ret = CAEN_DGTZ_SetChannelSelfTrigger(handle,CAEN_DGTZ_TRGMODE_ACQ_ONLY,0x1<<triggerSource);
    }
    else{
      if(kModel==5740){
        ret = CAEN_DGTZ_SetExtTriggerInputMode(handle, CAEN_DGTZ_TRGMODE_DISABLED);
        ret = CAEN_DGTZ_SetChannelSelfTrigger(handle,CAEN_DGTZ_TRGMODE_DISABLED,0xFF);
        ret = CAEN_DGTZ_SetGroupSelfTrigger(handle,CAEN_DGTZ_TRGMODE_ACQ_ONLY,(0x1<<0));
      }
      else{
        printf("Invalid Trigger Source, setting Ch0 as default source.\n");
        ret = CAEN_DGTZ_SetExtTriggerInputMode(handle, CAEN_DGTZ_TRGMODE_DISABLED);
        ret = CAEN_DGTZ_SetChannelSelfTrigger(handle,CAEN_DGTZ_TRGMODE_DISABLED,0xFF);
        ret = CAEN_DGTZ_SetChannelSelfTrigger(handle,CAEN_DGTZ_TRGMODE_ACQ_ONLY,(0x1<<0));
      }
    }
  }


  //setCoincidence(0);
  //setCoincidence(2);
  //setCoincidence(4);
  //setCoincidence(6);
  //setMajorCoincidence(0xe,1,0);
  startSWAcq();

  ret = CAEN_DGTZ_ReadRegister(handle,0x810C,&dat);
  //printf("\n hola line %d, reg: %x\n",__LINE__,dat);
  while (count < events)
  {
    if (automatic)
      ret = CAEN_DGTZ_SendSWtrigger(handle);

    ret = CAEN_DGTZ_ReadData(handle,CAEN_DGTZ_SLAVE_TERMINATED_READOUT_MBLT,buffer,(uint32_t *)&bsize);
    //entrega en numEvents el numero de eventos capturados
    uint32_t numEvents;
    ret = CAEN_DGTZ_GetNumEvents(handle,buffer,bsize,&numEvents);
    printf("\revents: %d",count);
    fflush(stdout);
    Float_t *data_arr = new Float_t[NCh+2];
    for (int32_t i=0;i< int32_t(numEvents);i++)
    {
      ret = CAEN_DGTZ_GetEventInfo(handle,buffer,bsize,i,&eventInfo,(char **)&evtptr);
      ret = CAEN_DGTZ_DecodeEvent(handle,evtptr,(void **)&Evt);

      for (int32_t j=0;j<int32_t(Evt->ChSize[0]);j++)
      {
	       for (int32_t k=0;k<NCh;k++)
	        {
	           data_arr[k] = (int32_t)Evt->DataChannel[k][j];
             data_arr[k] = adc2mV(data_arr[k]);
	        }
	        data_arr[NCh]=j*(kSamplingTime*1000000000); // Tiempo de sampling
	        data_arr[NCh+1]=count+i+start_event;
	        data->Fill(data_arr);
      }

      ret = CAEN_DGTZ_FreeEvent(handle,(void **)&Evt);
    }

    delete data_arr;
    count +=numEvents;
  }

  std::cout<<std::endl;
  stopSWAcq();
  return 0;
}

int32_t  CaenDigitizerSiLab::readEvents(int32_t maxEvents,bool automatic,int32_t start_event,double tm, uint32_t triggerSource) // tm= timeout
{
  timeval ti,tf;
  double time_elapsed = 0.0;
  int32_t count=0;
  uint32_t dat=0;

  //Se libera el buffer y se vuelve a ubicar, por si han habido cambios en las configuraciones de la ventana de adquisicion
  ret = CAEN_DGTZ_FreeReadoutBuffer(&buffer);
  ret = CAEN_DGTZ_MallocReadoutBuffer(handle,&buffer,(uint32_t *)&size);
  
  //ret = CAEN_DGTZ_MallocReadoutBuffer(handle,&buffer,(uint32_t *)&size);
  if (!automatic)
  {
    //ret = CAEN_DGTZ_SetChannelSelfTrigger(handle,CAEN_DGTZ_TRGMODE_ACQ_ONLY,(3<<6)); //Adjacent channels paired.
    //ret = CAEN_DGTZ_SetChannelSelfTrigger(handle,CAEN_DGTZ_TRGMODE_ACQ_ONLY,(0x3f<<2)); //Adjacent channels paired.
    if(triggerSource == 9)
    {
      /* 
      ret = CAEN_DGTZ_SetChannelSelfTrigger(handle,CAEN_DGTZ_TRGMODE_DISABLED,0x00);
      ret = CAEN_DGTZ_SetExtTriggerInputMode(handle, CAEN_DGTZ_TRGMODE_ACQ_ONLY);
      printf("\n");
      */
      ret = CAEN_DGTZ_SetExtTriggerInputMode(handle, CAEN_DGTZ_TRGMODE_ACQ_ONLY);

      CAEN_DGTZ_IOLevel_t myLevel={};
      CAEN_DGTZ_GetIOLevel (handle, &myLevel);
      
      CAEN_DGTZ_SetIOLevel (handle, CAEN_DGTZ_IOLevel_NIM);
      
      CAEN_DGTZ_GetIOLevel (handle, &myLevel);
      
    }
    else if ((0<=triggerSource)&&(triggerSource<=7)&&(kModel!=5740))
    {
      ret = CAEN_DGTZ_SetExtTriggerInputMode(handle, CAEN_DGTZ_TRGMODE_DISABLED);
      ret = CAEN_DGTZ_SetChannelSelfTrigger(handle,CAEN_DGTZ_TRGMODE_ACQ_ONLY,0x0<<triggerSource);
     //ret = CAEN_DGTZ_SetChannelSelfTrigger(handle,CAEN_DGTZ_TRGMODE_ACQ_ONLY,0x0F);
    }
    else
    {
      if(kModel==5740){
        ret = CAEN_DGTZ_SetExtTriggerInputMode(handle, CAEN_DGTZ_TRGMODE_DISABLED);
        ret = CAEN_DGTZ_SetGroupSelfTrigger(handle,CAEN_DGTZ_TRGMODE_ACQ_ONLY,(0x1<<0));
      }
      else
      {
        printf("Invalid Trigger Source, setting Ch0 as default source.\n");
        ret = CAEN_DGTZ_SetExtTriggerInputMode(handle, CAEN_DGTZ_TRGMODE_DISABLED);
        ret = CAEN_DGTZ_SetChannelSelfTrigger(handle,CAEN_DGTZ_TRGMODE_ACQ_ONLY,(0x1<<0));
      }
    }
  }

  //setCoincidence(0); // par 0 1 
  //setCoincidence(2); // par 2 3
  //setCoincidence(4); // par 4 5
  //setCoincidence(6); // par 6 7
  //setMajorCoincidence(0xe,1,0);
  //setCoincidence(0);
  //setCoincidence(2);
  //setAllCoincidencesToAND();
  //forceGlobalTriggerMask(0xFF,10);
  //enableAllChannelsORTrigger();
  //enableChannelPairsANDTrigger(8, 0);
  printTriggerConfiguration();
  //setMajorCoincidence(0x1, 0,0);
  startSWAcq();

  ret = CAEN_DGTZ_ReadRegister(handle,0x810C,&dat);
  gettimeofday(&ti,NULL);
  while ((count<maxEvents)&&(time_elapsed<tm)&&(quit!=1))
  {
    if (automatic)
      ret = CAEN_DGTZ_SendSWtrigger(handle);

    ret = CAEN_DGTZ_ReadData(handle,CAEN_DGTZ_SLAVE_TERMINATED_READOUT_MBLT,buffer,(uint32_t *)&bsize);
    //entrega en numEvents el numero de eventos capturados
    uint32_t numEvents;
    ret = CAEN_DGTZ_GetNumEvents(handle,buffer,bsize,&numEvents);
    printf("\revents: %d",count);
    fflush(stdout);
    auto now = steady_clock::now();
    double timestamp = duration_cast<microseconds>(now.time_since_epoch()).count();
    Float_t *data_arr = new Float_t[NCh+3];
    for (int32_t i=0;i<int32_t(numEvents);i++)
    {
      ret = CAEN_DGTZ_GetEventInfo(handle,buffer,bsize,i,&eventInfo,(char **)&evtptr);
      ret = CAEN_DGTZ_DecodeEvent(handle,evtptr,(void **)&Evt);

      for (int32_t j=0;j<int32_t(Evt->ChSize[0]);j++)
      {
	       for (int32_t k=0;k<NCh;k++)
	        {
	           data_arr[k] = (int32_t)Evt->DataChannel[k][j];
             data_arr[k] = adc2mV(data_arr[k]);
	        }
	        data_arr[NCh]=j*(kSamplingTime*1000000000); // Tiempo de sampling
	        data_arr[NCh+1]=count+i+start_event;
          data_arr[NCh+2]=timestamp;
	        data->Fill(data_arr);
      }

      ret = CAEN_DGTZ_FreeEvent(handle,(void **)&Evt);
    }

    delete data_arr;
    count +=numEvents;
    gettimeofday(&tf,NULL);
    time_elapsed = tf.tv_sec-ti.tv_sec + tf.tv_usec/1e6-ti.tv_usec/1e6;
  }
  std::cout<<std::endl;
  printf("time elapsed: %5.2f\n<rate>: %3.2f\n",time_elapsed,count/time_elapsed);
	printf("Retrieved %d Event\n",count);
  stopSWAcq();
  return 0;
}

int32_t CaenDigitizerSiLab::getTempMeanVar()
  {
    if (kTempSupported)
    {
      for (int32_t i =0; i<MaxNCh; i++)
      {
        meanTemp[i]=0;
        varTemp[i]=0;
      }

      for (int32_t k=0;k<kNSamplesTemp;k++)
      {
        for (int32_t i =0; i<MaxNCh; i++)
        {
  	       if ((kEnableMask>>i)&0x1)
           {
          	  readTemp(i);
          	  meanTemp[i]+=chTemp;
          	  varTemp[i]+=chTemp*chTemp;
  	       }
         }
         sleep(kDtTemp);
      }

      for (int32_t i =0; i<MaxNCh; i++)
      {
        meanTemp[i]/=kNSamplesTemp;
        varTemp[i]=(varTemp[i]/kNSamplesTemp-meanTemp[i]*meanTemp[i]);
      }
    }
    else
    {
      std::cout<<"Temperature measurement not supported"<<std::endl;
      return -1;
    }

    return 0;

  }


int32_t CaenDigitizerSiLab::storeData()
{
  ofile->cd();
  data->Write("",TObject::kOverwrite);
  //ofile->Close();
  // tempFile.close();
  return 0;
}

int32_t CaenDigitizerSiLab::setCoincidence(int32_t ch0)
{
  if (ch0 % 2 != 0)
  {
    printf("Coincidence not possible, use only with even channels\n");
    return -1;
  }

  int32_t ret = 0;
  uint32_t data = 0;
  uint32_t reg = 0x1000 + (ch0 << 8) + 0x84; // Dirección 0x1n84 para pareja ch0/ch0+1

  ret = CAEN_DGTZ_ReadRegister(handle, reg, &data);

  uint32_t opt = 0x0; // 00: AND
  data = (data & ~0x3) | opt;

  ret = CAEN_DGTZ_WriteRegister(handle, reg, data);
  ret = CAEN_DGTZ_ReadRegister(handle, reg, &data); // Verificación opcional

  return ret;
}


int32_t CaenDigitizerSiLab::setMajorCoincidence(int32_t blkmask, int32_t wd, int32_t level)
{
  uint32_t data = 0;
  uint32_t reg = 0x810C;

  // Leer valor actual
  ret = CAEN_DGTZ_ReadRegister(handle, reg, &data);

  // Limpiar bits [7:0], [23:20] y [26:24]
  data &= ~((0xFF) | (0xF << 20) | (0x7 << 24));

  // Establecer nuevos valores
  data |= (blkmask & 0xFF);        // parejas activadas
  data |= ((wd & 0xF) << 20);      // ventana de coincidencia
  data |= ((level & 0x7) << 24);   // nivel de mayoría

  // Escribir valor nuevo
  ret = CAEN_DGTZ_WriteRegister(handle, reg, data);

  // Leer de nuevo para verificar
  ret = CAEN_DGTZ_ReadRegister(handle, reg, &data);

  return ret;
}

int32_t CaenDigitizerSiLab::setAllCoincidencesToAND()
{
  int32_t ret = 0;
  uint32_t data = 0;
  uint32_t reg = 0;

  const int numChannels = 8; // Cambia esto si usas un modelo con más canales (por ejemplo 16)

  for (int ch = 0; ch < numChannels; ch += 2) {
    reg = 0x1000 + (ch << 8) + 0x84; // dirección 0x1n84
    ret = CAEN_DGTZ_ReadRegister(handle, reg, &data);

    data = (data & ~0x3) | 0x0; // Set bits [1:0] = 00 → AND

    ret = CAEN_DGTZ_WriteRegister(handle, reg, data);
    ret = CAEN_DGTZ_ReadRegister(handle, reg, &data); // lectura opcional para verificación
  }

  return ret;
}

int32_t CaenDigitizerSiLab::forceGlobalTriggerMask(uint8_t pair_mask, uint8_t coinc_window)
{
  int32_t ret = 0;
  uint32_t data = 0;
  uint32_t reg = 0x810C;

  // Leer valor actual
  ret |= CAEN_DGTZ_ReadRegister(handle, reg, &data);

  // Limpiar: bits [7:0] (pair mask), [23:20] (coinc window), [26:24] (majority)
  data &= ~((0xFF) | (0xF << 20) | (0x7 << 24));

  // Aplicar nueva configuración
  data |= (pair_mask & 0xFF);           // Qué parejas están activas
  data |= ((coinc_window & 0xF) << 20); // Ventana de coincidencia
  data |= (0 << 24);                    // Desactiva majority (nivel 0)

  // Escribir y verificar
  ret = CAEN_DGTZ_WriteRegister(handle, reg, data);
  ret = CAEN_DGTZ_ReadRegister(handle, reg, &data);

  // Tclk para DT5730
    constexpr double Tclk_ns = 8.0; // 8 ns por ciclo

  // Opcional: imprimir confirmación
  printf("Trigger mask set: 0x%02X, coincidence window: %.1f ns \n", pair_mask, coinc_window*Tclk_ns);

  return ret;
}

int32_t CaenDigitizerSiLab::finish()
{
  ret = CAEN_DGTZ_FreeReadoutBuffer(&buffer);	//Frees memory allocated by the CAEN_DGTZ_MallocReadoutBuffer
  std::cout<<"Free Readout Buffer. ret.: "<<ret<<std::endl;
	if (ret==0)
  {
    ret = CAEN_DGTZ_CloseDigitizer(handle);
    std::cout<<"Digitizer Comunication Closed. ret.: "<<ret<<std::endl;
    return ret;
  }
  return -1;
}

CaenDigitizerSiLab::~CaenDigitizerSiLab()
{
  if (kTempSupported)
  {
    delete [] meanTemp;
    delete [] varTemp;
    delete [] storedMeanTemp;
    delete [] storedVarTemp;
    tempFile.close();
  }
  //ofile->Close();

  //  CAEN_DGTZ_FreeEvent(handle,(void **)&Evt);
}

int32_t CaenDigitizerSiLab::configureGlobalTrigger(uint8_t pair_mask, uint8_t coinc_window, uint8_t majority_level)
{
    int32_t ret = 0;
    uint32_t data = 0;
    uint32_t reg = 0x810C;

    // Leer valor actual
    ret = CAEN_DGTZ_ReadRegister(handle, reg, &data);
    if (ret != 0) {
        printf("Error reading trigger register (0x810C)\n");
        return ret;
    }

    // Limpiar campos: [7:0]=parejas, [23:20]=ventana, [26:24]=mayoría
    data &= ~((0xFF) | (0xF << 20) | (0x7 << 24));

    // Construir nuevo valor
    data |= (pair_mask & 0xFF);                 // habilitación de parejas
    data |= ((coinc_window & 0xF) << 20);       // ventana de coincidencia
    data |= ((majority_level & 0x7) << 24);     // nivel de mayoría

    // Escribir configuración
    ret = CAEN_DGTZ_WriteRegister(handle, reg, data);
    if (ret != 0) {
        printf("Error writing trigger register (0x810C)\n");
        return ret;
    }

    // Leer nuevamente para verificar
    ret = CAEN_DGTZ_ReadRegister(handle, reg, &data);
    if (ret != 0) {
        printf("Error verifying trigger register (0x810C)\n");
        return ret;
    }

    // Tclk para DT5730
    constexpr double Tclk_ns = 8.0; // 8 ns por ciclo

    printf(">>> Global Trigger Configuration:\n");
    printf("  - Enabled Pairs Mask      = 0x%02X\n", pair_mask);
    printf("  - Coincidence Window      = %d clock cycles = %.1f ns\n", coinc_window, coinc_window * Tclk_ns);
    if (majority_level > 0)
        printf("  - Majority Level Required = %d channels\n", majority_level);
    else
        printf("  - Majority Trigger        = DISABLED\n");

    return 0;
}

void CaenDigitizerSiLab::printTriggerConfiguration() {
  if (handle < 0) {
      printf("Digitizer handle is invalid.\n");
      return;
  }

  uint32_t data = 0;
  int ret;

  printf("\n========== TRIGGER CONFIGURATION (DT5730) ==========\n");

  // 1. Trigger Thresholds por canal (CH0–CH7)
  printf("Trigger Thresholds per Channel:\n");
  for (int ch = 0; ch < 8; ++ch) {
      uint32_t addr = 0x1080 + (ch << 8);  // Threshold register
      ret = CAEN_DGTZ_ReadRegister(handle, addr, &data);
      if (ret == 0)
          printf("  - Channel %d: Threshold = %u\n", ch, data & 0x3FFF);
      else
          printf("  - Channel %d: Error reading threshold\n", ch);
  }

  // 2. Trigger Source Enable Mask (0x810C)
  ret = CAEN_DGTZ_ReadRegister(handle, 0x810C, &data);
  if (ret == 0) {
      printf("\nTrigger Source Enable Mask (0x810C): 0x%08X\n", data);

      uint8_t pair_mask      = data & 0xFF;
      uint8_t coinc_window   = (data >> 20) & 0xF;
      uint8_t majority_level = (data >> 24) & 0x7;
      constexpr double Tclk_ns = 8.0;

      printf("  - Enabled Pairs Mask      = 0x%02X\n", pair_mask);
      printf("  - Coincidence Window      = %u clock cycles = %.1f ns\n", coinc_window, coinc_window * Tclk_ns);
      if (majority_level > 0)
          printf("  - Majority Level Required = %u channels\n", majority_level);
      else
          printf("  - Majority Trigger        = DISABLED\n");

      // Bits 30: software trigger, 31: external trigger
      printf("  - Software Trigger        = %s\n", (data & (1U << 30)) ? "ENABLED" : "DISABLED");
      printf("  - External Trigger        = %s\n", (data & (1U << 31)) ? "ENABLED" : "DISABLED");

      // Canal individual habilitado
      printf("  - Channel Trigger Enable:\n");
      for (int ch = 0; ch < 8; ++ch) {
          bool enabled = data & (1 << ch);
          printf("    · CH%d: %s\n", ch, enabled ? "ENABLED" : "DISABLED");
      }

  } else {
      printf("\nTrigger Source Enable Mask: Error reading\n");
  }

  // 3. Coincidence Logic (even channels only)
  printf("\nCoincidence Logic Between Neighbor Channels:\n");
  for (int ch = 0; ch < 8; ch += 2) {
      uint32_t addr = 0x1084 + (ch << 8);
      ret = CAEN_DGTZ_ReadRegister(handle, addr, &data);
      if (ret == 0) {
          uint32_t mode = data & 0x3;
          const char* mode_str[] = {"AND", "Only N", "Only N+1", "OR"};
          printf("  - CH%d/CH%d: Mode = %s (0x%X)\n", ch, ch + 1, mode_str[mode], mode);
      } else {
          printf("  - CH%d/CH%d: Error reading coincidence mode\n", ch, ch + 1);
      }
  }

  printf("====================================================\n");
}


int32_t CaenDigitizerSiLab::enableAllChannelsORTrigger() {
  int32_t ret = 0;
  uint32_t data = 0;

  // 1. Configurar todos los pares en modo OR
  for (int ch = 0; ch < 8; ch += 2) {
      uint32_t addr = 0x1084 + (ch << 8); // Dirección del registro de coincidencia
      ret = CAEN_DGTZ_ReadRegister(handle, addr, &data);
      if (ret != 0) {
          printf("Error leyendo modo de coincidencia en CH%d/CH%d\n", ch, ch + 1);
          return ret;
      }

      // Set bits [1:0] = 0x3 (modo OR)
      data = (data & ~0x3) | 0x3;
      ret = CAEN_DGTZ_WriteRegister(handle, addr, data);
      if (ret != 0) {
          printf("Error configurando modo OR en CH%d/CH%d\n", ch, ch + 1);
          return ret;
      }
  }

  // 2. Configurar 0x810C: activar todos los pares, coincidencia de 4 ciclos, mayoría desactivada
  uint8_t pair_mask = 0xFF;          // CH0–CH7 (pares CH0/1, CH2/3, ...)
  uint8_t coinc_window = 4;          // 4 × 8 ns = 32 ns
  uint8_t majority_level = 0;        // desactivado

  // Leer valor actual de 0x810C
  ret = CAEN_DGTZ_ReadRegister(handle, 0x810C, &data);
  if (ret != 0) {
      printf("Error leyendo registro 0x810C\n");
      return ret;
  }

  // Limpiar bits [7:0], [23:20], [26:24]
  data &= ~((0xFF) | (0xF << 20) | (0x7 << 24));

  // Insertar nuevos valores
  data |= (pair_mask & 0xFF);
  data |= ((coinc_window & 0xF) << 20);
  data |= ((majority_level & 0x7) << 24);

  ret = CAEN_DGTZ_WriteRegister(handle, 0x810C, data);
  if (ret != 0) {
      printf("Error escribiendo registro 0x810C\n");
      return ret;
  }

  printf("Todos los pares configurados en modo OR, habilitados para trigger.\n");
  printf("Coincidence window = %u cycles (%.1f ns), Majority = DISABLED\n",
         coinc_window, coinc_window * 8.0);
  return 0;
}

int32_t CaenDigitizerSiLab::enableChannelPairsANDTrigger(uint8_t coinc_window, uint8_t majority_level) {
  int32_t ret = 0;
  uint32_t data = 0;

  // Validación básica
  if (coinc_window > 15 || majority_level > 7) {
      printf("Parámetros inválidos: coinc_window ≤ 15, majority_level ≤ 7\n");
      return -1;
  }

  // 1. Configurar todos los pares en modo AND
  for (int ch = 0; ch < 8; ch += 2) {
      uint32_t addr = 0x1084 + (ch << 8); // Dirección del registro de coincidencia CHch/CHch+1
      ret = CAEN_DGTZ_ReadRegister(handle, addr, &data);
      if (ret != 0) {
          printf("Error leyendo modo de coincidencia en CH%d/CH%d\n", ch, ch + 1);
          return ret;
      }

      // Modo AND: bits [1:0] = 0x0
      data = (data & ~0x3) | 0x0;
      ret = CAEN_DGTZ_WriteRegister(handle, addr, data);
      if (ret != 0) {
          printf("Error configurando modo AND en CH%d/CH%d\n", ch, ch + 1);
          return ret;
      }
  }

  // 2. Configurar el registro 0x810C
  uint8_t pair_mask = 0xFF;  // Todas las parejas habilitadas

  ret = CAEN_DGTZ_ReadRegister(handle, 0x810C, &data);
  if (ret != 0) {
      printf("Error leyendo registro 0x810C\n");
      return ret;
  }

  // Limpiar bits [7:0], [23:20], [26:24]
  data &= ~((0xFF) | (0xF << 20) | (0x7 << 24));

  // Insertar configuración
  data |= (pair_mask & 0xFF);
  data |= ((coinc_window & 0xF) << 20);
  data |= ((majority_level & 0x7) << 24);

  ret = CAEN_DGTZ_WriteRegister(handle, 0x810C, data);
  if (ret != 0) {
      printf("Error escribiendo registro 0x810C\n");
      return ret;
  }

  printf("Pares CH0/1 a CH6/7 configurados en modo AND y habilitados para trigger global.\n");
  printf("Ventana de coincidencia: %u ciclos (%.1f ns), Nivel de mayoría: %s\n",
         coinc_window, coinc_window * 8.0,
         majority_level ? std::to_string(majority_level).c_str() : "DESACTIVADO");

  return 0;
}
