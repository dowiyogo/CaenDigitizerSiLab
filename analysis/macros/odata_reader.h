//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jun 12 18:41:43 2025 by ROOT version 6.32.10
// from TTree odata/Output Data
// found on file: output_tree_file.root
//////////////////////////////////////////////////////////

#ifndef odata_reader_h
#define odata_reader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#define SIZE 1500

// Header file for the classes stored in the TTree if any.

class odata_reader {
public :
   TTree          *fChain, *oTree;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   TFile          *outFile;

// Fixed size dimensions of array or collections stored in the TTree if any.

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
   Int_t           time_ns[SIZE];
   Float_t         pre_event;

   // Declaration of leaf types
   Int_t           size;
   Float_t         Ch0[10000];   //[size]
   Float_t         Ch1[10000];   //[size]
   Float_t         Ch2[10000];   //[size]
   Float_t         Ch3[10000];   //[size]
   Float_t         Ch4[10000];   //[size]
   Float_t         Ch5[10000];   //[size]
   Float_t         Ch6[10000];   //[size]
   Float_t         Ch7[10000];   //[size]
   Float_t         event;
   Int_t           timens[10000];   //[size]
   Float_t         timestamp;

   // List of branches
   TBranch        *b_size;   //!
   TBranch        *b_Ch0;   //!
   TBranch        *b_Ch1;   //!
   TBranch        *b_Ch2;   //!
   TBranch        *b_Ch3;   //!
   TBranch        *b_Ch4;   //!
   TBranch        *b_Ch5;   //!
   TBranch        *b_Ch6;   //!
   TBranch        *b_Ch7;   //!
   TBranch        *b_event;   //!
   TBranch        *b_timens;   //!
   TBranch        *b_timestamp;   //!

   odata_reader(TTree *tree=0);
   virtual ~odata_reader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual bool     Notify();
   virtual void     Show(Long64_t entry = -1);
   void PlotInterEventTimeDistribution(int ch,
      float threshold = -50.f,
      int   nbins      = 100);
   TH2I* MakeCoincidenceMap(float threshold = -50.f);
   TH1F* AllInterEventTimeDistribution(int nbins = 100);
};

#endif

#ifdef odata_reader_cxx
odata_reader::odata_reader(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("output_tree_file.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("output_tree_file.root");
      }
      f->GetObject("odata",tree);

   }
   Init(tree);
}

odata_reader::~odata_reader()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t odata_reader::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t odata_reader::LoadTree(Long64_t entry)
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

void odata_reader::Init(TTree *tree)
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
   outFile = new TFile("reduced_tree_file.root","RECREATE");
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
   oTree->Branch("timens",time_ns,"timens[size]/I");
   oTree->Branch("timestamp",&timestamp,"timestamp/F");

   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("size", &size, &b_size);
   fChain->SetBranchAddress("Ch0", Ch0, &b_Ch0);
   fChain->SetBranchAddress("Ch1", Ch1, &b_Ch1);
   fChain->SetBranchAddress("Ch2", Ch2, &b_Ch2);
   fChain->SetBranchAddress("Ch3", Ch3, &b_Ch3);
   fChain->SetBranchAddress("Ch4", Ch4, &b_Ch4);
   fChain->SetBranchAddress("Ch5", Ch5, &b_Ch5);
   fChain->SetBranchAddress("Ch6", Ch6, &b_Ch6);
   fChain->SetBranchAddress("Ch7", Ch7, &b_Ch7);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("timens", timens, &b_timens);
   fChain->SetBranchAddress("timestamp", &timestamp, &b_timestamp);
   Notify();
}

bool odata_reader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

void odata_reader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t odata_reader::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef odata_reader_cxx
