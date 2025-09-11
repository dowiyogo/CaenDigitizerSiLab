#define myclass_cxx
#include "myclass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void myclass::Loop()
{
//   In a ROOT session, you can do:
//      root> .L myclass.C
//      root> myclass t
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

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
}
void myclass::PlotInterEventTimeDistribution()
{
  if (!fChain) {
    std::cerr << "fChain no está inicializado." << std::endl;
    return;
  }

  Long64_t nentries = fChain->GetEntries();
  if (nentries < 2) {
    std::cerr << "No hay suficientes eventos para calcular diferencias." << std::endl;
    return;
  }

  std::vector<float> timestamps;
  timestamps.reserve(nentries);

  for (Long64_t i = 0; i < nentries; ++i) {
    GetEntry(i);
    timestamps.push_back(timestamp);
  }

  // Ordenamos por si acaso (aunque debiera venir ordenado)
  std::sort(timestamps.begin(), timestamps.end());

  // Calculamos las diferencias entre timestamps consecutivos
  std::vector<float> deltas;
  deltas.reserve(nentries - 1);

  for (size_t i = 1; i < timestamps.size(); ++i) {
    float dt = timestamps[i] - timestamps[i - 1];
    if (dt > 0)
      deltas.push_back(dt);
  }

  // Crear histograma con delta t en microsegundos
  float max_dt = *std::max_element(deltas.begin(), deltas.end());
  int nbins = 100;
  TH1F* hDeltaT = new TH1F("hDeltaT", "Diferencias de tiempo entre eventos;#Delta t (#mus);Frecuencia", nbins, 0, max_dt);

  for (float dt : deltas) {
    hDeltaT->Fill(dt);
  }

  TCanvas* c1 = new TCanvas("c1", "Distribucion de Tiempos entre Eventos", 800, 600);
  hDeltaT->Draw();

  // Opción para ajustar a exponencial si se desea
  TF1* fitExp = new TF1("fitExp", "[0]*exp(-x/[1])", 0, max_dt);
  fitExp->SetParameters(hDeltaT->GetMaximum(), hDeltaT->GetMean());
  hDeltaT->Fit(fitExp, "R");

  std::cout << "Ajuste exponencial: y = A * exp(-x / tau)" << std::endl;
  std::cout << "A = " << fitExp->GetParameter(0)
            << ", tau = " << fitExp->GetParameter(1) << " microsegundos" << std::endl;
}

void myclass::PlotInterEventTimeDistribution(int ch, float threshold)
{
  if (!fChain) {
    std::cerr << "fChain no está inicializado." << std::endl;
    return;
  }

  if (ch < 0 || ch >= 8) {
    std::cerr << "Canal inválido. Debe estar entre 0 y 7." << std::endl;
    return;
  }

  Long64_t nentries = fChain->GetEntries();
  if (nentries < 2) {
    std::cerr << "No hay suficientes eventos para calcular diferencias." << std::endl;
    return;
  }

  std::vector<float> timestamps;
  timestamps.reserve(nentries);

  for (Long64_t i = 0; i < nentries; ++i) {
    GetEntry(i);

    bool pass = false;

    switch (ch) {
      case 0: for (float v : *Ch0) { if (v < threshold) { pass = true; break; } } break;
      case 1: for (float v : *Ch1) { if (v < threshold) { pass = true; break; } } break;
      case 2: for (float v : *Ch2) { if (v < threshold) { pass = true; break; } } break;
      case 3: for (float v : *Ch3) { if (v < threshold) { pass = true; break; } } break;
      case 4: for (float v : *Ch4) { if (v < threshold) { pass = true; break; } } break;
      case 5: for (float v : *Ch5) { if (v < threshold) { pass = true; break; } } break;
      case 6: for (float v : *Ch6) { if (v < threshold) { pass = true; break; } } break;
      case 7: for (float v : *Ch7) { if (v < threshold) { pass = true; break; } } break;
    }

    if (pass) {
      timestamps.push_back(timestamp);
    }
  }

  if (timestamps.size() < 2) {
    std::cerr << "Muy pocos eventos cumplen el criterio (valor < " << threshold << ")." << std::endl;
    return;
  }

  std::sort(timestamps.begin(), timestamps.end());

  std::vector<float> deltas;
  deltas.reserve(timestamps.size() - 1);

  for (size_t i = 1; i < timestamps.size(); ++i) {
    float dt = timestamps[i] - timestamps[i - 1];
    if (dt > 0) deltas.push_back(dt);
  }

  float max_dt = *std::max_element(deltas.begin(), deltas.end());
  int nbins = 100;
  TH1F* hDeltaT = new TH1F("hDeltaT",
    Form("Delta t para Ch%d < %.1f;#Delta t (#mus);Frecuencia", ch, threshold),
    nbins, 0, max_dt);

  for (float dt : deltas) {
    hDeltaT->Fill(dt);
  }

  TCanvas* c1 = new TCanvas("c1", "Distribucion de Tiempo entre Eventos (Filtrados)", 800, 600);
  hDeltaT->Draw();

  TF1* fitExp = new TF1("fitExp", "[0]*exp(-x/[1])", 0, max_dt);
  fitExp->SetParameters(hDeltaT->GetMaximum(), hDeltaT->GetMean());
  hDeltaT->Fit(fitExp, "R");

  std::cout << "Eventos seleccionados: " << timestamps.size() << std::endl;
  std::cout << "Ajuste exponencial: y = A * exp(-x / tau)" << std::endl;
  std::cout << "A = " << fitExp->GetParameter(0)
            << ", tau = " << fitExp->GetParameter(1) << " microsegundos" << std::endl;
}
