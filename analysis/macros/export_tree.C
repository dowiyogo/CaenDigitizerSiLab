#define export_tree_cxx
#include "export_tree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <algorithm>
#include <vector>
#include <bitset>
#include <TH2I.h>

void export_tree::Loop()
{
   //ROOT::EnableImplicitMT(); //Tell ROOT you want to go parallel
   //ROOT::EnableThreadSafety();
//   In a ROOT session, you can do:
//      root> .L export_tree.C
//      root> export_tree t
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
   //Long64_t nentries = 100000;
   fChain->GetEntry(0);
   pre_event = event;

   Long64_t nbytes = 0, nb = 0, cnt=0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if (pre_event!=event){
         oTree->Fill();
         cnt=0;
      }
      pre_event=event;
      Ch0a[cnt] = Ch0;
      Ch1a[cnt] = Ch1;
      Ch2a[cnt] = Ch2;
      Ch3a[cnt] = Ch3;
      Ch4a[cnt] = Ch4;
      Ch5a[cnt] = Ch5;
      Ch6a[cnt] = Ch6;
      Ch7a[cnt] = Ch7;
      timens[cnt] = cnt*8; 
      cnt++;
      if(!(jentry%(10*SIZE)))std::cout<<Form("progress: %3.1f",float(jentry)/nentries*100.)<<"\r";
   }
   

   oTree->Fill();
   outFile->Write();
   //outFile->Close();
}
