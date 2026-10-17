//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jun 11 18:03:08 2025 by ROOT version 6.32.10
// from TTree data/amp (adc ch) and time (nsample) time_stamp (us)
// found on file: trigger_externo_50000_ch1-5NEW-ch6-8OLD.root
//////////////////////////////////////////////////////////

#ifndef export_tree_h
#define export_tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#define SIZE 10000
// Header file for the classes stored in the TTree if any.

class export_tree {
public :
   TTree          *fChain, *oTree;   //!pointer to the analyzed TTree or TChain
   TFile          *outFile;
   Int_t           fCurrent; //!current Tree number in a TChain
   Int_t chsize;// = SIZE;
// Fixed size dimensions of array or collections stored in the TTree if any.
   Float_t         Ch0a[SIZE];
   Float_t         Ch1a[SIZE];
   Float_t         Ch2a[SIZE];
   Float_t         Ch3a[SIZE];
   Float_t         Ch4a[SIZE];
   Float_t         Ch5a[SIZE];
   Float_t         Ch6a[SIZE];
   Float_t         Ch7a[SIZE];
   Int_t         timens[SIZE];
   Float_t         pre_event;

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

   export_tree(TTree *tree=0);
   virtual ~export_tree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual bool     Notify();
   virtual void     Show(Long64_t entry = -1);

};

#endif

#ifdef export_tree_cxx
export_tree::export_tree(TTree *tree) : fChain(0) 
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

export_tree::~export_tree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t export_tree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t export_tree::LoadTree(Long64_t entry)
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

void export_tree::Init(TTree *tree)
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
   chsize = SIZE;
   outFile = new TFile("output_tree_file.root","RECREATE");
   oTree = new TTree("odata","Output Data");
   oTree->Branch("size",&chsize,"size/I");
   oTree->Branch("Ch0",Ch0a,"Ch0[size]/F");
   oTree->Branch("Ch1",Ch1a,"Ch1[size]/F");
   oTree->Branch("Ch2",Ch2a,"Ch2[size]/F");
   oTree->Branch("Ch3",Ch3a,"Ch3[size]/F");
   oTree->Branch("Ch4",Ch4a,"Ch4[size]/F");
   oTree->Branch("Ch5",Ch5a,"Ch5[size]/F");
   oTree->Branch("Ch6",Ch6a,"Ch6[size]/F");
   oTree->Branch("Ch7",Ch7a,"Ch7[size]/F");
   oTree->Branch("event",&event,"event/F");
   oTree->Branch("timens",timens,"timens[size]/I");
   oTree->Branch("timestamp",&timestamp,"timestamp/F");

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

bool export_tree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

void export_tree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t export_tree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef export_tree_cxx
