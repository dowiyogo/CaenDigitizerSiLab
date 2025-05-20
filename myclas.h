//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon May 12 12:41:43 2025 by ROOT version 6.32.06
// from TTree data/amp (adc ch) and time (nsample) time_stamp (us)
// found on file: trigger_externo_50000_ch1-5NEW-ch6-8OLD.root
//////////////////////////////////////////////////////////

#ifndef myclas_h
#define myclas_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class myclas {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   //Float_t Chn[8][10000];

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         Ch0;
   Float_t         Ch1;
   Float_t         Ch2;
   Float_t         Ch3;
   Float_t         Ch4;
   Float_t         Ch5;
   Float_t         Ch6;
   Float_t         Ch7;
   Float_t         time;
   Float_t         event;
   Float_t         timestamp;

   // List of branches
   TBranch        *b_Ch0;   //!
   TBranch        *b_Ch1;   //!
   TBranch        *b_Ch2;   //!
   TBranch        *b_Ch3;   //!
   TBranch        *b_Ch4;   //!
   TBranch        *b_Ch5;   //!
   TBranch        *b_Ch6;   //!
   TBranch        *b_Ch7;   //!
   TBranch        *b_time;   //!
   TBranch        *b_event;   //!
   TBranch        *b_timestamp;   //!

   myclas(TTree *tree=0);
   virtual ~myclas();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual bool     Notify();
   virtual void     Show(Long64_t entry = -1);
   void ExportToStructuredTree(const char* outputFileName);
};

#endif

#ifdef myclas_cxx
myclas::myclas(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("trigger_externo_50000_ch1-5NEW-ch6-8OLD.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("trigger_externo_50000_ch1-5NEW-ch6-8OLD.root");
      }
      f->GetObject("data",tree);

   }
   Init(tree);
}

myclas::~myclas()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t myclas::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t myclas::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void myclas::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Ch0", &Ch0, &b_Ch0);
   fChain->SetBranchAddress("Ch1", &Ch1, &b_Ch1);
   fChain->SetBranchAddress("Ch2", &Ch2, &b_Ch2);
   fChain->SetBranchAddress("Ch3", &Ch3, &b_Ch3);
   fChain->SetBranchAddress("Ch4", &Ch4, &b_Ch4);
   fChain->SetBranchAddress("Ch5", &Ch5, &b_Ch5);
   fChain->SetBranchAddress("Ch6", &Ch6, &b_Ch6);
   fChain->SetBranchAddress("Ch7", &Ch7, &b_Ch7);
   fChain->SetBranchAddress("time", &time, &b_time);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("timestamp", &timestamp, &b_timestamp);
   Notify();
}

bool myclas::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

void myclas::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t myclas::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef myclas_cxx
