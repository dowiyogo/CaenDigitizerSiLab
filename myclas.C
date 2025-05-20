#define myclas_cxx
#include "myclas.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void myclas::Loop()
{
   ROOT::EnableImplicitMT(); //Tell ROOT you want to go parallel
   ROOT::EnableThreadSafety();
//   In a ROOT session, you can do:
//      root> .L myclas.C
//      root> myclas t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   // Salida
   TFile *outFile = new TFile("trigger_analysis.root", "RECREATE");

   // Árbol de salida con branches originales + trigger
   TTree *newTree = new TTree("data_with_trigger", "Original data + trigger");

   // Variables de entrada
   Float_t in_Ch[8];
   Float_t in_time, in_event, in_timestamp;

   // Asociar branches originales
   fChain->SetBranchAddress("Ch0", &in_Ch[0]);
   fChain->SetBranchAddress("Ch1", &in_Ch[1]);
   fChain->SetBranchAddress("Ch2", &in_Ch[2]);
   fChain->SetBranchAddress("Ch3", &in_Ch[3]);
   fChain->SetBranchAddress("Ch4", &in_Ch[4]);
   fChain->SetBranchAddress("Ch5", &in_Ch[5]);
   fChain->SetBranchAddress("Ch6", &in_Ch[6]);
   fChain->SetBranchAddress("Ch7", &in_Ch[7]);
   fChain->SetBranchAddress("time", &in_time);
   fChain->SetBranchAddress("event", &in_event);
   fChain->SetBranchAddress("timestamp", &in_timestamp);

   // Variables de salida
   Float_t out_Ch[8];
   Float_t out_time, out_event, out_timestamp;
   UChar_t out_trigger = 0;

   // Asociar branches de salida
   newTree->Branch("Ch0", &out_Ch[0], "Ch0/F");
   newTree->Branch("Ch1", &out_Ch[1], "Ch1/F");
   newTree->Branch("Ch2", &out_Ch[2], "Ch2/F");
   newTree->Branch("Ch3", &out_Ch[3], "Ch3/F");
   newTree->Branch("Ch4", &out_Ch[4], "Ch4/F");
   newTree->Branch("Ch5", &out_Ch[5], "Ch5/F");
   newTree->Branch("Ch6", &out_Ch[6], "Ch6/F");
   newTree->Branch("Ch7", &out_Ch[7], "Ch7/F");
   newTree->Branch("time", &out_time, "time/F");
   newTree->Branch("event", &out_event, "event/F");
   newTree->Branch("timestamp", &out_timestamp, "timestamp/F");
   newTree->Branch("trigger", &out_trigger, "trigger/b"); // 8-bit mask

   // Mapa de máscaras por evento
   std::map<int, UChar_t> triggerMap;

   // Histograma 2D coincidencias canal-canal
   TH2I *h2 = new TH2I("hCoincidencias", "Coincidencias entre canales",
                       8, -0.5, 7.5, 8, -0.5, 7.5); // Canal vs Canal

   // Recorrido
   for (Long64_t jentry = 0; jentry < nentries; ++jentry) {
      fChain->GetEntry(jentry);

      int ev = static_cast<int>(in_event);
      if (triggerMap.find(ev) == triggerMap.end())
         triggerMap[ev] = 0;

      for (int ch = 0; ch < 8; ++ch) {
         if (in_Ch[ch] < -50)
            triggerMap[ev] |= (1 << ch);
      }
   }

   // Segunda pasada: llenar nuevo árbol y contar coincidencias
   for (Long64_t jentry = 0; jentry < nentries; ++jentry) {
      fChain->GetEntry(jentry);

      int ev = static_cast<int>(in_event);
      out_event = in_event;
      out_timestamp = in_timestamp;
      out_time = in_time;
      for (int i = 0; i < 8; ++i)
         out_Ch[i] = in_Ch[i];

      out_trigger = triggerMap[ev];
      newTree->Fill();

      // Contar coincidencias solo una vez por evento
      if (in_time == 0) {
         std::vector<int> active;
         for (int i = 0; i < 8; ++i) {
            if ((out_trigger >> i) & 1)
               active.push_back(i);
         }

         // Llenar histograma 2D para todas combinaciones
         for (size_t i = 0; i < active.size(); ++i)
            for (size_t j = i; j < active.size(); ++j)
               h2->Fill(active[i], active[j]);
      }
   }

   // Guardar todo
   outFile->cd();
   newTree->Write();
   h2->Write();
   outFile->Close();

   std::cout << "Análisis completado. Archivo 'trigger_analysis.root' generado.\n";
}

void myclas::ExportToStructuredTree(const char* outputFileName)
{
    // Abrir archivo de salida
    ROOT::EnableImplicitMT(); //Tell ROOT you want to go parallel
    ROOT::EnableThreadSafety();
    TFile* fout = new TFile(outputFileName, "RECREATE");
    TTree* tree = new TTree("Events", "Eventos estructurados");

    // Variables para almacenar los datos por evento
    std::vector<float> Ch[8];
    std::vector<float> vtime;
    float event_out;
    float timestamp_out;

    // Setear branches
    tree->Branch("Ch0", &Ch[0]);
    tree->Branch("Ch1", &Ch[1]);
    tree->Branch("Ch2", &Ch[2]);
    tree->Branch("Ch3", &Ch[3]);
    tree->Branch("Ch4", &Ch[4]);
    tree->Branch("Ch5", &Ch[5]);
    tree->Branch("Ch6", &Ch[6]);
    tree->Branch("Ch7", &Ch[7]);
    tree->Branch("time", &vtime);
    tree->Branch("event", &event_out);
    tree->Branch("timestamp", &timestamp_out);

    // Variables para leer
    float Chval[8], time, event, timestamp;

    // Set branch addresses (ajusta los nombres si es necesario)
    fChain->SetBranchAddress("Ch0", &Chval[0]);
    fChain->SetBranchAddress("Ch1", &Chval[1]);
    fChain->SetBranchAddress("Ch2", &Chval[2]);
    fChain->SetBranchAddress("Ch3", &Chval[3]);
    fChain->SetBranchAddress("Ch4", &Chval[4]);
    fChain->SetBranchAddress("Ch5", &Chval[5]);
    fChain->SetBranchAddress("Ch6", &Chval[6]);
    fChain->SetBranchAddress("Ch7", &Chval[7]);
    fChain->SetBranchAddress("time", &time);
    fChain->SetBranchAddress("event", &event);
    fChain->SetBranchAddress("timestamp", &timestamp);

    // Loop principal
    Long64_t nentries = fChain->GetEntries();
    float last_event = -1;
    for (Long64_t i=0; i<nentries; ++i)
    {
        fChain->GetEntry(i);
        if (i == 0 || event != last_event)
        {
            // Si no es el primer evento, guardar el anterior
            if (i != 0)
            {
                event_out = last_event;
                timestamp_out = timestamp; // toma el timestamp de la última muestra (o promedia si prefieres)
                tree->Fill();
            }
            // Limpiar los vectores
            for (int c=0; c<8; ++c) Ch[c].clear();
            vtime.clear();
            last_event = event;
        }
        // Llenar vectores con las muestras del evento actual
        for (int c=0; c<8; ++c) Ch[c].push_back(Chval[c]);
        vtime.push_back(time);
    }
    // Guardar el último evento
    event_out = last_event;
    timestamp_out = timestamp;
    tree->Fill();

    fout->Write();
    fout->Close();
}