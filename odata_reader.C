#define odata_reader_cxx
#include "odata_reader.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <algorithm>
#include <vector>
#include <bitset>
#include <TH2I.h>

void odata_reader::Loop()
{
   if (!fChain) return;

   const int kFirst = 7500;         // inclusive
   const int kLast  = 9000;         // exclusivo
   const int kSave  = kLast - kFirst;   // 1500

   Long64_t nentries = fChain->GetEntriesFast();
   for (Long64_t ev = 0; ev < nentries; ++ev)
   {
      fChain->GetEntry(ev);           // carga TODA la traza del evento

      /* ── Copiar el sub-rango de muestras ───────────────────────────── */
      int idx = 0;
      for (int i = kFirst; i < kLast && i < size; ++i)   // ‘size’ viene del árbol
      {
         Ch0a[idx] = Ch0[i];
         Ch1a[idx] = Ch1[i];
         Ch2a[idx] = Ch2[i];
         Ch3a[idx] = Ch3[i];
         Ch4a[idx] = Ch4[i];
         Ch5a[idx] = Ch5[i];
         Ch6a[idx] = Ch6[i];
         Ch7a[idx] = Ch7[i];
         time_ns[idx] = i * 8;        // 8 ns por sample
         ++idx;
      }

      /* ── Guarda el evento recortado ────────────────────────────────── */
      oTree->Fill();

      if (!(ev % 1000))
         std::cout << Form("progreso: %5.1f %%\r",
                           100.*ev/nentries) << std::flush;
   }

   outFile->Write();
}

void odata_reader::PlotInterEventTimeDistribution(int ch,
                                                  float threshold,
                                                  int   nbins)
{
   /* ─── Comprobaciones básicas ─────────────────────────────────────────── */
   if (!fChain) {
      std::cerr << "fChain no está inicializado.\n";
      return;
   }
   if (ch < 0 || ch > 7) {
      std::cerr << "Canal inválido (debe ser 0-7).\n";
      return;
   }

   const Long64_t nentries = fChain->GetEntriesFast();
   if (nentries < 2) {
      std::cerr << "No hay suficientes eventos.\n";
      return;
   }

   /* ─── 1. Recorremos el árbol y seleccionamos eventos ─────────────────── */
   std::vector<float> ts;               // timestamps que cumplen el corte
   ts.reserve(nentries);

   for (Long64_t j = 0; j < nentries; ++j) {
      fChain->GetEntry(j);              // llena los miembros de la clase

      /* puntero al canal que interesa */
      const float* chan = nullptr;
      switch (ch) {
         case 0: chan = Ch0; break; case 1: chan = Ch1; break;
         case 2: chan = Ch2; break; case 3: chan = Ch3; break;
         case 4: chan = Ch4; break; case 5: chan = Ch5; break;
         case 6: chan = Ch6; break; case 7: chan = Ch7; break;
      }

      /* ¿Algún sample < threshold? */
      bool pass = false;
      for (int i = 0; i < size; ++i) {
         if (chan[i] < threshold) { pass = true; break; }
      }
      if (pass) ts.push_back(timestamp);   // timestamp es por evento
   }

   if (ts.size() < 2) {
      std::cerr << "Muy pocos eventos bajo el umbral " << threshold << ".\n";
      return;
   }

   /* ─── 2. Calculamos Δt entre timestamps sucesivos ────────────────────── */
   std::sort(ts.begin(), ts.end());

   std::vector<float> dt; dt.reserve(ts.size() - 1);
   for (std::size_t i = 1; i < ts.size(); ++i) {
      const float delta = ts[i] - ts[i - 1];
      if (delta > 0) dt.push_back(delta);
   }
   if (dt.empty()) {
      std::cerr << "No hay diferencias positivas de tiempo.\n";
      return;
   }

   /* ─── 3. Histograma y ajuste exponencial ─────────────────────────────── */
   const float max_dt = *std::max_element(dt.begin(), dt.end());
   TH1F* h = new TH1F("hDeltaT",
                      Form("#Delta t para Ch%d < %.1f ;#Delta t (#mu s);Conteo", ch, threshold),
                      nbins, 0, max_dt);

   for (float x : dt) h->Fill(x);

   TCanvas* c = new TCanvas("cInterEvt", "Distribucion #Delta t", 900, 600);
   c->SetLogy();               // escala log en Y: más cómoda para exp.
   h->Draw();

   TF1* fexp = new TF1("fexp", "[0]*exp(-x/[1])", 0, max_dt);
   fexp->SetParameters(h->GetMaximum(), h->GetMean());
   h->Fit(fexp, "RQ");         // R = rango definido, Q = sin printout

   /* ─── 4. Resumen en consola ──────────────────────────────────────────── */
   std::cout << "Eventos seleccionados: " << ts.size() << '\n'
             << "Ajuste exponencial  y = A·exp(-x/τ)\n"
             << "   A   = " << fexp->GetParameter(0) << '\n'
             << "   τ   = " << fexp->GetParameter(1) << " μs\n";
}

TH2I* odata_reader::MakeCoincidenceMap(float threshold)
{
   const Long64_t nentries = fChain->GetEntriesFast();
   if (nentries == 0) {
      std::cerr << "Árbol vacío.\n";
      return nullptr;
   }

   /* Histograma 8×8  */
   TH2I *hCoin = new TH2I("hCoin",
                          "Coincidencias entre canales;Canal i;Canal j",
                          8, -0.5, 7.5,   // eje X
                          8, -0.5, 7.5);  // eje Y

   /* Tabla de punteros a los 8 canales para recorrerlos con un bucle */
   const Float_t* chPtr[8] = { Ch0, Ch1, Ch2, Ch3, Ch4, Ch5, Ch6, Ch7 };

   int      lastEvent = -1;    // -1 → todavía no hemos procesado nada
   uint8_t  mask      = 0;     // bit-mask de canales activos en evento actual

   auto flushMask = [&](uint8_t m)
   {
      /* Inserta todos los pares (i,j) contenidos en la máscara */
      for (int i = 0; i < 8; ++i) if (m & (1u << i))
         for (int j = i; j < 8; ++j) if (m & (1u << j))
            hCoin->Fill(i, j);
   };

   /* ─── BARRIDO ÚNICO SOBRE EL ÁRBOL ─────────────────────────────────── */
   for (Long64_t ie = 0; ie < nentries; ++ie) {
      fChain->GetEntry(ie);

      int ev = static_cast<int>(event);   // ‘event’ viene del MakeClass

      /* Si cambiamos de evento, volcamos la máscara y la reiniciamos */
      if (ev != lastEvent && lastEvent != -1) {
         flushMask(mask);
         mask = 0;
      }
      lastEvent = ev;

      /* Actualizamos la máscara con los canales que superen el umbral  */
      for (int ch = 0; ch < 8; ++ch) {
         const Float_t* data = chPtr[ch];          // Arreglo ChX[size]
         for (int i = 0; i < size; ++i) {        // ‘size’ → nº samples
            if (data[i] < threshold) {           // ← criterio de disparo
               mask |= (1u << ch);
               break;                            // basta un sample
            }
         }
      }
   }

   /* Procesamos el último evento */
   if (lastEvent != -1) flushMask(mask);

   /* ─── Vista rápida ─────────────────────────────────────────────────── */
   auto c = new TCanvas("cCoin", "Coincidencias", 650, 600);
   c->SetRightMargin(0.14);   // que no tape la barra de colores
   hCoin->Draw("COLZ text");  // ‘text’ → cuenta visible en cada celda

   return hCoin;
}

TH1F* odata_reader::AllInterEventTimeDistribution(int nbins)
{
   /* ─── Validaciones básicas ─────────────────────────────────────────── */
   if (!fChain) {
      std::cerr << "fChain no está inicializado.\n";
      return nullptr;
   }
   const Long64_t nentries = fChain->GetEntriesFast();
   if (nentries < 2) {
      std::cerr << "Se necesitan al menos dos eventos.\n";
      return nullptr;
   }

   /* ─── 1. Copiamos todos los timestamps ─────────────────────────────── */
   std::vector<Float_t> ts; ts.reserve(nentries);
   for (Long64_t i = 0; i < nentries; ++i) {
      fChain->GetEntry(i);
      ts.push_back(timestamp);          // ‘timestamp’ es miembro autogenerado
   }

   std::sort(ts.begin(), ts.end());     // por si acaso

   /* ─── 2. Diferencias entre eventos sucesivos ───────────────────────── */
   std::vector<Float_t> dt; dt.reserve(ts.size() - 1);
   for (std::size_t i = 1; i < ts.size(); ++i) {
      const Float_t d = ts[i] - ts[i - 1];
      if (d > 0) dt.push_back(d);
   }
   if (dt.empty()) {
      std::cerr << "Todas las diferencias son 0 o negativas.\n";
      return nullptr;
   }

   /* ─── 3. Histograma Δt (μs) ─────────────────────────────────────────── */
   const Float_t max_dt = *std::max_element(dt.begin(), dt.end());
   const Float_t xmax   = (max_dt > 0) ? max_dt : 1.0e-6f;   // no rango 0

   TH1F *h = new TH1F("hDeltaT",
                      "Distribucion de #Delta t entre eventos;#Delta t (#mu s);Cuenta",
                      nbins, 0, xmax);

   for (Float_t x : dt) h->Fill(x);

   /* ─── 4. Dibujamos y (opcional) ajustamos exponencial ───────────────── */
   auto c = new TCanvas("cDeltaT", "Inter-Event #Delta t", 850, 600);
   c->SetLogy();          // escala log en Y facilita ver la cola
   h->Draw();

   if (dt.size() > 20) {  // ajustamos sólo si hay datos suficientes
      TF1 *fexp = new TF1("fexp", "[0]*exp(-x/[1])", 0, xmax);
      fexp->SetParameters(h->GetMaximum(), h->GetMean());
      h->Fit(fexp, "RQ");   // R = usa rango, Q = quiet (sin spam)
      std::cout << "Ajuste y = A · exp(-x/τ):  "
                << "A = "   << fexp->GetParameter(0)
                << ", τ = " << fexp->GetParameter(1) << " μs\n";
   }

   return h;   // el usuario decide si lo guarda
}

