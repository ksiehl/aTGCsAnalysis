#include <TFile.h>
#include <TTree.h>
#include <TTreeFormula.h>
#include <TBranch.h>
#include <TH1.h>
#include <TH2.h>
#include <THStack.h>
#include <TPad.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TDirectoryFile.h>
#include <iostream>
#include <iostream>


/*
 * Small macro to add weights corresponding to the cross-sections
 * 
 * Compile by: 
 * root -l addWeightSamples.cpp+
 */

double Nevents(std::string filename){
   TFile file(filename.c_str(), "READ");
  TTree * tree = (TTree*) file.Get("GenWeights/Tree");
  double gentWeight;
  tree -> SetBranchAddress("genWeight", &gentWeight);
  double sum = 0;
  for (unsigned int iEntry = 0; iEntry < tree -> GetEntries(); iEntry ++)
  {
    tree -> GetEntry(iEntry); 
    //sum += gentWeight/std::abs(gentWeight);
    sum += gentWeight;
  }
  return sum;
}
void addWeight(string FileName, float xsection, float lumi, std::string channel)
{
  double Nevents_ = Nevents(FileName);
  TFile file(FileName.c_str(), "UPDATE");
  TTree * tree = (TTree*) file.Get("treeDumper/BasicTree");
  double totWeight, topPtSF, btagWeight, totWeight_BTagUp, totWeight_BTagDown, totWeight_MistagUp, totWeight_MistagDown, totWeight_LeptonIDUp, totWeight_LeptonIDDown, triggerWeightHLTEle27NoER;
  tree -> SetBranchAddress("totWeight", &totWeight);
  tree -> SetBranchAddress("topPtSF", &topPtSF);
  tree -> SetBranchAddress("btagWeight", &btagWeight);
  tree -> SetBranchAddress("totWeight_BTagUp", &totWeight_BTagUp);
  tree -> SetBranchAddress("totWeight_BTagDown", &totWeight_BTagDown);
  tree -> SetBranchAddress("totWeight_MistagUp", &totWeight_MistagUp);
  tree -> SetBranchAddress("totWeight_MistagDown", &totWeight_MistagDown);
  tree -> SetBranchAddress("totWeight_LeptonIDUp", &totWeight_LeptonIDUp);
  tree -> SetBranchAddress("totWeight_LeptonIDDown", &totWeight_LeptonIDDown);
  if (channel == "ele") tree -> SetBranchAddress("triggerWeightHLTEle27NoER", &triggerWeightHLTEle27NoER);
  double totWeightWithLumi, totWeightWithLumiNoBtag, totWeightWithLumi_MistagUp, totWeightWithLumi_MistagDown, totWeightWithLumi_BTagUp, totWeightWithLumi_BTagDown, totWeightWithLumi_LeptonIDUp, totWeightWithLumi_LeptonIDDown;
  TBranch * br = tree -> Branch("totEventWeight", &totWeightWithLumi, "totEventWeight/D"); 
  TBranch * br_NoBtag = tree -> Branch("totEventWeightNoBtag", &totWeightWithLumiNoBtag, "totEventWeightNoBtag/D");
  TBranch * br_MistagUp = tree -> Branch("totEventWeight_MistagUp", &totWeightWithLumi_MistagUp, "totEventWeight_MistagUp/D"); 
  TBranch * br_MistagDown = tree -> Branch("totEventWeight_MistagDown", &totWeightWithLumi_MistagDown, "totEventWeight_MistagDown/D"); 
  TBranch * br_BTagUp = tree -> Branch("totEventWeight_BTagUp", &totWeightWithLumi_BTagUp, "totEventWeight_BTagUp/D"); 
  TBranch * br_BTagDown = tree -> Branch("totEventWeight_BTagDown", &totWeightWithLumi_BTagDown, "totEventWeight_BTagDown/D");
  TBranch * br_LeptonIDUp = tree -> Branch("totEventWeight_LeptonIDUp", &totWeightWithLumi_LeptonIDUp, "totEventWeight_LeptonIDUp/D");
  TBranch * br_LeptonIDDown = tree -> Branch("totEventWeight_LeptonIDDown", &totWeightWithLumi_LeptonIDDown, "totEventWeight_LeptonIDDown/D");
  std::cout << FileName << std::endl;
  std::cout << "Number of events (effective):" << Nevents_ << std::endl;

  // Calculate mean top Pt SF first
  double meanTopWeight = 0.;
  for (unsigned int iEntry = 0; iEntry < tree->GetEntries(); iEntry++) {
    tree->GetEntry(iEntry);
    meanTopWeight += topPtSF;
  }
  meanTopWeight /= tree->GetEntries();
  // need to divide by mean top Pt weight as it only affects shape, not total # events
  double weightLumi = (xsection*lumi)/(Nevents_*meanTopWeight);

  for (unsigned int iEntry = 0; iEntry < tree -> GetEntries(); iEntry ++)
  {
    tree -> GetEntry(iEntry); 
    // if (channel == "ele"){
    //   totWeightWithLumi = totWeight*weightLumi/triggerWeightHLTEle27NoER;
    //   totWeightWithLumi_BTagUp = totWeight_BTagUp*weightLumi/triggerWeightHLTEle27NoER;
    //   totWeightWithLumi_BTagDown = totWeight_BTagDown*weightLumi/triggerWeightHLTEle27NoER;
    //   totWeightWithLumi_MistagUp= totWeight_MistagUp*weightLumi/triggerWeightHLTEle27NoER;
    //   totWeightWithLumi_MistagDown = totWeight_MistagDown*weightLumi/triggerWeightHLTEle27NoER;
    // }
    // else { 
      totWeightWithLumi = totWeight*weightLumi;
      totWeightWithLumiNoBtag = totWeightWithLumi/btagWeight;
      totWeightWithLumi_BTagUp = totWeight_BTagUp*weightLumi;
      totWeightWithLumi_BTagDown = totWeight_BTagDown*weightLumi;
      totWeightWithLumi_MistagUp= totWeight_MistagUp*weightLumi;
      totWeightWithLumi_MistagDown = totWeight_MistagDown*weightLumi;
      totWeightWithLumi_LeptonIDUp= totWeight_LeptonIDUp*weightLumi;
      totWeightWithLumi_LeptonIDDown = totWeight_LeptonIDDown*weightLumi;
    // }
    br -> Fill();
    br_NoBtag -> Fill();
    br_MistagUp -> Fill();
    br_MistagDown -> Fill();
    br_BTagUp -> Fill();
    br_BTagDown -> Fill();
    br_LeptonIDUp -> Fill();
    br_LeptonIDDown -> Fill();
  }
  
  tree -> Write();

}

void addWeightSamples()
{
  double lumi = 35922.;
  std::string prefix = "/afs/cern.ch/work/m/maiqbal/private/aTGC/Samples_80X_Working/";

  //electron channel
  addWeight(prefix + "WW_ele.root", 49.997, lumi, "ele");
  addWeight(prefix + "WZ_ele.root", 11.46, lumi, "ele");
  addWeight(prefix + "s-ch_ele.root", 10.32*0.33, lumi, "ele");
  addWeight(prefix + "t-ch-top_ele.root", 136.02, lumi, "ele");
  addWeight(prefix + "t-ch-antitop_ele.root", 80.95, lumi, "ele");
  addWeight(prefix + "tW-ch-top_ele.root", 35.6, lumi, "ele");
  addWeight(prefix + "tW-ch-antitop_ele.root", 35.6, lumi, "ele");
  addWeight(prefix + "WJets_Ht100To200_ele.root", 1345.0*1.21, lumi, "ele");
  addWeight(prefix + "WJets_Ht200To400_ele.root", 359.7*1.21, lumi, "ele"); 
  addWeight(prefix + "WJets_Ht400To600_ele.root", 48.91*1.21, lumi, "ele"); 
  addWeight(prefix + "WJets_Ht600To800_ele.root", 12.05*1.21, lumi, "ele"); 
  addWeight(prefix + "WJets_Ht800To1200_ele.root", 5.501*1.21, lumi, "ele");
  addWeight(prefix + "WJets_Ht1200To2500_ele.root", 1.329*1.21, lumi, "ele"); 
  addWeight(prefix + "WJets_Ht2500ToInf_ele.root", 0.03216*1.21, lumi, "ele"); 
  addWeight(prefix + "ttbar_ele.root", 831.76, lumi, "ele");

  // The cross sections are from the GenXSecAnalyzer and multiplying factors are required to match the SM yield to the SM samples.
  addWeight(prefix + "WW-aTGC_MWW-600To800_ele.root", 0.1833 * 3.499, lumi, "ele");
  addWeight(prefix + "WW-aTGC_MWW-800ToInf_ele.root", 0.2366 * 2.816, lumi, "ele");
  addWeight(prefix + "WZ-aTGC_MWZ-600To800_ele.root", 0.06493 * 13.16, lumi, "ele");
  addWeight(prefix + "WZ-aTGC_MWZ-800ToInf_ele.root", 0.1012 * 8.81, lumi, "ele");


  //muon channel
  addWeight(prefix + "WW_mu.root", 49.997, lumi, "");
  addWeight(prefix + "WZ_mu.root", 11.46, lumi, "");
  addWeight(prefix + "s-ch_mu.root", 10.32*0.33, lumi,"");
  addWeight(prefix + "t-ch-top_mu.root", 136.02, lumi,""); 
  addWeight(prefix + "t-ch-antitop_mu.root", 80.95, lumi,""); 
  addWeight(prefix + "tW-ch-top_mu.root", 35.6, lumi,""); 
  addWeight(prefix + "tW-ch-antitop_mu.root", 35.6, lumi,""); 
  addWeight(prefix + "WJets_Ht100To200_mu.root", 1345.0*1.21, lumi,""); 
  addWeight(prefix + "WJets_Ht200To400_mu.root", 359.7*1.21, lumi,"");
  addWeight(prefix + "WJets_Ht400To600_mu.root", 48.91*1.21, lumi,"");
  addWeight(prefix + "WJets_Ht600To800_mu.root", 12.05*1.21, lumi,"");
  addWeight(prefix + "WJets_Ht800To1200_mu.root", 5.501*1.21, lumi,"");
  addWeight(prefix + "WJets_Ht1200To2500_mu.root", 1.329*1.21, lumi,""); 
  addWeight(prefix + "WJets_Ht2500ToInf_mu.root", 0.03216*1.21, lumi,"");
  addWeight(prefix + "ttbar_mu.root", 831.76, lumi,"");

  // The cross sections are from the GenXSecAnalyzer and multiplying factors are required to match the SM yield to the SM samples.
  addWeight(prefix + "WW-aTGC_MWW-600To800_mu.root", 0.1833 * 3.77, lumi, "");
  addWeight(prefix + "WW-aTGC_MWW-800ToInf_mu.root", 0.2366 * 3.179, lumi, "");
  addWeight(prefix + "WZ-aTGC_MWZ-600To800_mu.root", 0.06493 * 15.15, lumi, "");  
  addWeight(prefix + "WZ-aTGC_MWZ-800ToInf_mu.root", 0.1012 * 9.654, lumi, "");

}
