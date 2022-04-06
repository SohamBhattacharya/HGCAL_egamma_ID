#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TPad.h"
#include "TROOT.h"


int test()
{
    gStyle->SetOptStat(11111111);
    
    TChain *chain1 = new TChain("chain");
    
    TChain *inCh1 = new TChain("treeMaker/tree");
    inCh1->Add("ntupleTree_rerunTICL_modTICLeleWithRerunTICL_onRaw_2.root");
    inCh1->SetWeight(0.5);
    //inCh1->SetWeight(0.5, "global");
    
    chain1->Add(inCh1);
    
    //inCh1->Draw("gsfEleFromTICL_pT");
    chain1->Draw("gsfEleFromTICL_pT");
    
    return 0;
}
