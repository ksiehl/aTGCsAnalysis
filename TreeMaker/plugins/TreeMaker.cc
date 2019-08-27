// system include files
#include <iostream>
#include <memory>
#include <map>
#include <vector>
#include <algorithm> 
#include "boost/algorithm/string.hpp"
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "FWCore/Framework/interface/EventSetup.h"


#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/BTauReco/interface/CATopJetTagInfo.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"

#include "TTree.h"
#include "TFile.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TRandom.h"

#include "DecayChannel.h"
#include "GenWUtils.h"
#include "DecayClass.h"
#include "Particle.h"
#include "SystematicsHelper.h"
#include "PU.h"
#include "PDFVariationMap.h"
#include "ScaleFactorHelper.h"
#include "BTagHelper.h"
#include "JetResolutionSmearer.h"
#include "Ele27WPLooseTrigTurnOn.h"

#include "aTGCsAnalysis/metZcalc/METzCalculator.cc"
#include "aTGCsAnalysis/metZcalc/METxyCorrector.cc"

//#define ANGLE_TESTING

namespace reco {
  typedef edm::Ptr<reco::Muon> MuonPtr;
  typedef edm::Ptr<reco::GsfElectron> ElectronPtr;
}

Double_t d_costheta1 = 0.0;
Double_t d_costheta2 = 0.0;
Double_t d_phi = 0.0;
Double_t d_costhetastar = 0.0;
Double_t d_phi1 = 0.0;
Double_t d_phi2 = 0.0;

void calculateAngles(TLorentzVector thep4M11, TLorentzVector thep4M12, TLorentzVector thep4M21, TLorentzVector thep4M22, Double_t& costheta1, Double_t& costheta2, Double_t& phi, Double_t& costhetastar, Double_t& phi1, Double_t& phi2);

const Int_t weights_number = 124;
const Bool_t METcorrect = false;
Int_t number_of_subjets = 0;
Int_t imaginary_neutrino = 0;

//Double_t lepton_xerror = -99.9;
//Double_t lepton_yerror = -99.9;
//Double_t lepton_zerror = -99.9;
//Double_t met_xerror = +99.9;
//Double_t met_yerror = +99.9;

#ifdef ANGLE_TESTING

Bool_t bad_angle = false;
//INTERMEDIATE STEPS VARIABLES
Double_t d_leptons_in_lep_px = 0.0;
Double_t d_leptons_in_lep_py = 0.0;
Double_t d_leptons_in_lep_pz = 0.0;

Double_t d_partons_in_lep_px = 0.0;
Double_t d_partons_in_lep_py = 0.0;
Double_t d_partons_in_lep_pz = 0.0;

Double_t d_parton1_in_lep_px = 0.0;
Double_t d_parton1_in_lep_py = 0.0;
Double_t d_parton1_in_lep_pz = 0.0;

Double_t d_parton2_in_lep_px = 0.0;
Double_t d_parton2_in_lep_py = 0.0;
Double_t d_parton2_in_lep_pz = 0.0;

Double_t d_lepton1_in_lep_px = 0.0;
Double_t d_lepton1_in_lep_py = 0.0;
Double_t d_lepton1_in_lep_pz = 0.0;

Double_t d_lepton1_dotted_x = 0.0;
Double_t d_lepton1_dotted_y = 0.0;
Double_t d_lepton1_dotted_z = 0.0;

Double_t d_leptons_in_had_px = 0.0;
Double_t d_leptons_in_had_py = 0.0;
Double_t d_leptons_in_had_pz = 0.0;

Double_t d_lepton1_in_had_px = 0.0;
Double_t d_lepton1_in_had_py = 0.0;
Double_t d_lepton1_in_had_pz = 0.0;

Double_t d_lepton2_in_had_px = 0.0;
Double_t d_lepton2_in_had_py = 0.0;
Double_t d_lepton2_in_had_pz = 0.0;

Double_t d_parton1_in_had_px = 0.0;
Double_t d_parton1_in_had_py = 0.0;
Double_t d_parton1_in_had_pz = 0.0;

Double_t d_parton1_dotted_x = 0.0;
Double_t d_parton1_dotted_y = 0.0;
Double_t d_parton1_dotted_z = 0.0;

Double_t d_complicated1_px = 0.0;
Double_t d_complicated1_py = 0.0;
Double_t d_complicated1_pz = 0.0;

Double_t d_complicated2_px = 0.0;
Double_t d_complicated2_py = 0.0;
Double_t d_complicated2_pz = 0.0;

Double_t d_lepton1WWframe_X = 0.0;
Double_t d_lepton1WWframe_Y = 0.0;
Double_t d_lepton1WWframe_Z = 0.0;

Double_t d_lepton_sumWWframe_X = 0.0;
Double_t d_lepton_sumWWframe_Y = 0.0;
Double_t d_lepton_sumWWframe_Z = 0.0;

Double_t d_parton1WWframe_X = 0.0;
Double_t d_parton1WWframe_Y = 0.0;
Double_t d_parton1WWframe_Z = 0.0;

Double_t d_parton_sumWWframe_X = 0.0;
Double_t d_parton_sumWWframe_Y = 0.0;
Double_t d_parton_sumWWframe_Z = 0.0;

Double_t d_boostWWframe_X = 0.0;
Double_t d_boostWWframe_Y = 0.0;
Double_t d_boostWWframe_Z = 0.0;

Double_t d_boostWlep_X = 0.0;
Double_t d_boostWlep_Y = 0.0;
Double_t d_boostWlep_Z = 0.0;

Double_t d_boostWhad_X = 0.0;
Double_t d_boostWhad_Y = 0.0;
Double_t d_boostWhad_Z = 0.0;

Double_t d_xdotx = 0.0;
Double_t d_xdoty = 0.0;
Double_t d_xdotz = 0.0;

Double_t d_ydotx = 0.0;
Double_t d_ydoty = 0.0;
Double_t d_ydotz = 0.0;

Double_t d_zdotx = 0.0;
Double_t d_zdoty = 0.0;
Double_t d_zdotz = 0.0;

Double_t d_lepton1WWframe_UX = 0.0; Double_t d_lepton1WWframe_UY = 0.0; Double_t d_lepton1WWframe_UZ = 0.0; Double_t d_lepton_sumWWframe_UX = 0.0; Double_t d_lepton_sumWWframe_UY = 0.0; Double_t d_lepton_sumWWframe_UZ = 0.0;

//////////////////////////////////////////////////////////////////

Int_t d_leptons_in_lep_px_good = 1;
Int_t d_leptons_in_lep_py_good = 1;
Int_t d_leptons_in_lep_pz_good = 1;

Int_t d_partons_in_lep_px_good = 1;
Int_t d_partons_in_lep_py_good = 1;
Int_t d_partons_in_lep_pz_good = 1;

Int_t d_parton1_in_lep_px_good = 1;
Int_t d_parton1_in_lep_py_good = 1;
Int_t d_parton1_in_lep_pz_good = 1;

Int_t d_parton2_in_lep_px_good = 1;
Int_t d_parton2_in_lep_py_good = 1;
Int_t d_parton2_in_lep_pz_good = 1;

Int_t d_lepton1_in_lep_px_good = 1;
Int_t d_lepton1_in_lep_py_good = 1;
Int_t d_lepton1_in_lep_pz_good = 1;

//Int_t d_lepton1_dotted_x_good = 1;
//Int_t d_lepton1_dotted_y_good = 1;
//Int_t d_lepton1_dotted_z_good = 1;

Int_t d_leptons_in_had_px_good = 1;
Int_t d_leptons_in_had_py_good = 1;
Int_t d_leptons_in_had_pz_good = 1;

Int_t d_lepton1_in_had_px_good = 1;
Int_t d_lepton1_in_had_py_good = 1;
Int_t d_lepton1_in_had_pz_good = 1;

Int_t d_lepton2_in_had_px_good = 1;
Int_t d_lepton2_in_had_py_good = 1;
Int_t d_lepton2_in_had_pz_good = 1;

Int_t d_parton1_in_had_px_good = 1;
Int_t d_parton1_in_had_py_good = 1;
Int_t d_parton1_in_had_pz_good = 1;

//Int_t d_parton1_dotted_x_good = 1;
//Int_t d_parton1_dotted_y_good = 1;
//Int_t d_parton1_dotted_z_good = 1;

//Int_t d_complicated1_px_good = 1;
//Int_t d_complicated1_py_good = 1;
//Int_t d_complicated1_pz_good = 1;

//Int_t d_complicated2_px_good = 1;
//Int_t d_complicated2_py_good = 1;
//Int_t d_complicated2_pz_good = 1;

Int_t d_lepton1WWframe_X_good = 1;
Int_t d_lepton1WWframe_Y_good = 1;
Int_t d_lepton1WWframe_Z_good = 1;

Int_t d_lepton_sumWWframe_X_good = 1;
Int_t d_lepton_sumWWframe_Y_good = 1;
Int_t d_lepton_sumWWframe_Z_good = 1;

Int_t d_parton1WWframe_X_good = 1;
Int_t d_parton1WWframe_Y_good = 1;
Int_t d_parton1WWframe_Z_good = 1;

Int_t d_parton_sumWWframe_X_good = 1;
Int_t d_parton_sumWWframe_Y_good = 1;
Int_t d_parton_sumWWframe_Z_good = 1;

//Int_t d_boostWWframe_X_good = 1;
//Int_t d_boostWWframe_Y_good = 1;
//Int_t d_boostWWframe_Z_good = 1;

//Int_t d_boostWlep_X_good = 1;
//Int_t d_boostWlep_Y_good = 1;
//Int_t d_boostWlep_Z_good = 1;

//Int_t d_boostWhad_X_good = 1;
//Int_t d_boostWhad_Y_good = 1;
//Int_t d_boostWhad_Z_good = 1;


//INTERMEDIATE STEPS
void intermediate_steps(TLorentzVector lepton1, TLorentzVector lepton2, TLorentzVector parton1, TLorentzVector parton2, Int_t neu_status, Double_t& leptons_in_lep_px, Double_t& leptons_in_lep_py, Double_t& leptons_in_lep_pz, Double_t& partons_in_lep_px, Double_t& partons_in_lep_py, Double_t& partons_in_lep_pz, Double_t& parton1_in_lep_px, Double_t& parton2_in_lep_px, Double_t& parton1_in_lep_py, Double_t& parton2_in_lep_py, Double_t& parton1_in_lep_pz, Double_t& parton2_in_lep_pz, Double_t& lepton1_in_lep_px, Double_t& lepton1_in_lep_py, Double_t& lepton1_in_lep_pz, Double_t& lepton1_dotted_x, Double_t& lepton1_dotted_y, Double_t& lepton1_dotted_z, Double_t& leptons_in_had_px, Double_t& leptons_in_had_py, Double_t& leptons_in_had_pz, Double_t& lepton1_in_had_px, Double_t& lepton1_in_had_py, Double_t& lepton1_in_had_pz, Double_t& lepton2_in_had_px, Double_t& lepton2_in_had_py, Double_t& lepton2_in_had_pz, Double_t& parton1_in_had_px, Double_t& parton1_in_had_py, Double_t& parton1_in_had_pz, Double_t& parton1_dotted_x, Double_t& parton1_dotted_y, Double_t& parton1_dotted_z, Double_t& complicated1_px, Double_t& complicated1_py, Double_t& complicated1_pz, Double_t& complicated2_px, Double_t& complicated2_py, Double_t& complicated2_pz, Double_t& lepton_sumWWframe_X, Double_t& lepton_sumWWframe_Y, Double_t& lepton_sumWWframe_Z, Double_t& lepton1WWframe_X, Double_t& lepton1WWframe_Y, Double_t& lepton1WWframe_Z, Double_t& parton_sumWWframe_X, Double_t& parton_sumWWframe_Y, Double_t& parton_sumWWframe_Z, Double_t& parton1WWframe_X, Double_t& parton1WWframe_Y, Double_t& parton1WWframe_Z, Double_t& costhetastar, Double_t& costheta1, Double_t& phi, Double_t& costheta2, Double_t& phi1, Double_t& phi2, Double_t& boostWWframe_X, Double_t& boostWWframe_Y, Double_t& boostWWframe_Z, Double_t& boostWlep_X, Double_t& boostWlep_Y, Double_t& boostWlep_Z, Double_t& boostWhad_X, Double_t& boostWhad_Y, Double_t& boostWhad_Z, Double_t& xdotx, Double_t& xdoty, Double_t& xdotz, Double_t& ydotx, Double_t& ydoty, Double_t& ydotz, Double_t& zdotx, Double_t& zdoty, Double_t& zdotz, Double_t& lepton1WWframe_UX, Double_t& lepton1WWframe_UY, Double_t& lepton1WWframe_UZ, Double_t& lepton_sumWWframe_UX, Double_t& lepton_sumWWframe_UY, Double_t& lepton_sumWWframe_UZ, Int_t& leptons_in_lep_px_good, Int_t& leptons_in_lep_py_good, Int_t& leptons_in_lep_pz_good, Int_t& partons_in_lep_px_good, Int_t& partons_in_lep_py_good, Int_t& partons_in_lep_pz_good, Int_t& parton1_in_lep_px_good, Int_t& parton2_in_lep_px_good, Int_t& parton1_in_lep_py_good, Int_t& parton2_in_lep_py_good, Int_t& parton1_in_lep_pz_good, Int_t& parton2_in_lep_pz_good, Int_t& lepton1_in_lep_px_good, Int_t& lepton1_in_lep_py_good, Int_t& lepton1_in_lep_pz_good, Int_t& leptons_in_had_px_good, Int_t& leptons_in_had_py_good, Int_t& leptons_in_had_pz_good, Int_t& lepton1_in_had_px_good, Int_t& lepton1_in_had_py_good, Int_t& lepton1_in_had_pz_good, Int_t& lepton2_in_had_px_good, Int_t& lepton2_in_had_py_good, Int_t& lepton2_in_had_pz_good, Int_t& parton1_in_had_px_good, Int_t& parton1_in_had_py_good, Int_t& parton1_in_had_pz_good, Int_t& lepton_sumWWframe_X_good, Int_t& lepton_sumWWframe_Y_good, Int_t& lepton_sumWWframe_Z_good, Int_t& lepton1WWframe_X_good, Int_t& lepton1WWframe_Y_good, Int_t& lepton1WWframe_Z_good, Int_t& parton_sumWWframe_X_good, Int_t& parton_sumWWframe_Y_good, Int_t& parton_sumWWframe_Z_good, Int_t& parton1WWframe_X_good, Int_t& parton1WWframe_Y_good, Int_t& parton1WWframe_Z_good);


 //INTERMEDIATE STEPS VARIABLES
  Double_t leptons_in_lep_px = -99.9;
  Double_t leptons_in_lep_py = -99.9;
  Double_t leptons_in_lep_pz = -99.9;
  
  Double_t partons_in_lep_px = -99.9;
  Double_t partons_in_lep_py = -99.9;
  Double_t partons_in_lep_pz = -99.9;
  
  Double_t parton1_in_lep_px = -99.9;
  Double_t parton1_in_lep_py = -99.9;  
  Double_t parton1_in_lep_pz = -99.9;

  Double_t parton2_in_lep_px = -99.9;
  Double_t parton2_in_lep_py = -99.9;
  Double_t parton2_in_lep_pz = -99.9;
  
  Double_t lepton1_in_lep_px = -99.9;
  Double_t lepton1_in_lep_py = -99.9;
  Double_t lepton1_in_lep_pz = -99.9;
  
  Double_t lepton1_dotted_x = -99.9;
  Double_t lepton1_dotted_y = -99.9;
  Double_t lepton1_dotted_z = -99.9;
  
  Double_t leptons_in_had_px = -99.9;
  Double_t leptons_in_had_py = -99.9;
  Double_t leptons_in_had_pz = -99.9;
  
  Double_t lepton1_in_had_px = -99.9;
  Double_t lepton1_in_had_py = -99.9;
  Double_t lepton1_in_had_pz = -99.9;
  
  Double_t lepton2_in_had_px = -99.9;
  Double_t lepton2_in_had_py = -99.9;
  Double_t lepton2_in_had_pz = -99.9;
  
  Double_t parton1_in_had_px = -99.9;
  Double_t parton1_in_had_py = -99.9;
  Double_t parton1_in_had_pz = -99.9;
  
  Double_t parton1_dotted_x = -99.9;
  Double_t parton1_dotted_y = -99.9;
  Double_t parton1_dotted_z = -99.9;
  
  Double_t complicated1_px = -99.9;
  Double_t complicated1_py = -99.9;
  Double_t complicated1_pz = -99.9;
  
  Double_t complicated2_px = -99.9;
  Double_t complicated2_py = -99.9;
  Double_t complicated2_pz = -99.9;
  
  Double_t lepton1WWframe_X = -99.9;
  Double_t lepton1WWframe_Y = -99.9;
  Double_t lepton1WWframe_Z = -99.9;
  
  Double_t lepton_sumWWframe_X = -99.9;
  Double_t lepton_sumWWframe_Y = -99.9;
  Double_t lepton_sumWWframe_Z = -99.9;
  
  Double_t parton1WWframe_X = -99.9;
  Double_t parton1WWframe_Y = -99.9;
  Double_t parton1WWframe_Z = -99.9;
  
  Double_t parton_sumWWframe_X = -99.9;
  Double_t parton_sumWWframe_Y = -99.9;
  Double_t parton_sumWWframe_Z = -99.9;

  Double_t boostWWframe_X = -99.9;
  Double_t boostWWframe_Y = -99.9;
  Double_t boostWWframe_Z = -99.9;
  
  Double_t boostWlep_X = -99.9;
  Double_t boostWlep_Y = -99.9;
  Double_t boostWlep_Z = -99.9;
  
  Double_t boostWhad_X = -99.9;
  Double_t boostWhad_Y = -99.9;
  Double_t boostWhad_Z = -99.9;

  Double_t xdotx = -99.9;
  Double_t xdoty = -99.9;
  Double_t xdotz = -99.9;

  Double_t ydotx = -99.9;
  Double_t ydoty = -99.9;
  Double_t ydotz = -99.9;

  Double_t zdotx = -99.9;
  Double_t zdoty = -99.9;
  Double_t zdotz = -99.9;

  Double_t lepton1WWframe_UX = -99.9; Double_t lepton1WWframe_UY = -99.9; Double_t lepton1WWframe_UZ = -99.9; Double_t lepton_sumWWframe_UX = -99.9; Double_t lepton_sumWWframe_UY = -99.9; Double_t lepton_sumWWframe_UZ = -99.9;

  ///////////////////////////////////////////////

  Int_t leptons_in_lep_px_good = -9;
  Int_t leptons_in_lep_py_good = -9;
  Int_t leptons_in_lep_pz_good = -9;
  
  Int_t partons_in_lep_px_good = -9;
  Int_t partons_in_lep_py_good = -9;
  Int_t partons_in_lep_pz_good = -9;
  
  Int_t parton1_in_lep_px_good = -9;
  Int_t parton1_in_lep_py_good = -9;  
  Int_t parton1_in_lep_pz_good = -9;

  Int_t parton2_in_lep_px_good = -9;
  Int_t parton2_in_lep_py_good = -9;
  Int_t parton2_in_lep_pz_good = -9;
  
  Int_t lepton1_in_lep_px_good = -9;
  Int_t lepton1_in_lep_py_good = -9;
  Int_t lepton1_in_lep_pz_good = -9;
  
  //Int_t lepton1_dotted_x_good = -9;
  //Int_t lepton1_dotted_y_good = -9;
  //Int_t lepton1_dotted_z_good = -9;
  
  Int_t leptons_in_had_px_good = -9;
  Int_t leptons_in_had_py_good = -9;
  Int_t leptons_in_had_pz_good = -9;
  
  Int_t lepton1_in_had_px_good = -9;
  Int_t lepton1_in_had_py_good = -9;
  Int_t lepton1_in_had_pz_good = -9;
  
  Int_t lepton2_in_had_px_good = -9;
  Int_t lepton2_in_had_py_good = -9;
  Int_t lepton2_in_had_pz_good = -9;
  
  Int_t parton1_in_had_px_good = -9;
  Int_t parton1_in_had_py_good = -9;
  Int_t parton1_in_had_pz_good = -9;
  
  //Int_t parton1_dotted_x_good = -9;
  //Int_t parton1_dotted_y_good = -9;
  //Int_t parton1_dotted_z_good = -9;
  
  //Int_t complicated1_px_good = -9;
  //Int_t complicated1_py_good = -9;
  //Int_t complicated1_pz_good = -9;
  
  //Int_t complicated2_px_good = -9;
  //Int_t complicated2_py_good = -9;
  //Int_t complicated2_pz_good = -9;
  
  Int_t lepton1WWframe_X_good = -9;
  Int_t lepton1WWframe_Y_good = -9;
  Int_t lepton1WWframe_Z_good = -9;
  
  Int_t lepton_sumWWframe_X_good = -9;
  Int_t lepton_sumWWframe_Y_good = -9;
  Int_t lepton_sumWWframe_Z_good = -9;
  
  Int_t parton1WWframe_X_good = -9;
  Int_t parton1WWframe_Y_good = -9;
  Int_t parton1WWframe_Z_good = -9;
  
  Int_t parton_sumWWframe_X_good = -9;
  Int_t parton_sumWWframe_Y_good = -9;
  Int_t parton_sumWWframe_Z_good = -9;

  //Int_t boostWWframe_X_good = -9;
  //Int_t boostWWframe_Y_good = -9;
  //Int_t boostWWframe_Z_good = -9;
  
  //Int_t boostWlep_X_good = -9;
  //Int_t boostWlep_Y_good = -9;
  //Int_t boostWlep_Z_good = -9;
  
  //Int_t boostWhad_X_good = -9;
  //Int_t boostWhad_Y_good = -9;
  //Int_t boostWhad_Z_good = -9;

  /////////////////////////////////////////////////////////////
  //*/
#endif
//
// class declaration
//

class TreeMaker : public edm::EDAnalyzer {
public:
  explicit TreeMaker(const edm::ParameterSet&);
  ~TreeMaker();
  //static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);
  
private:
  virtual void beginJob() override;
  virtual void beginRun(const edm::Run&, const edm::EventSetup&) override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endRun(edm::Run const& iEvent, edm::EventSetup const&) override;
  virtual void endJob() override;
  virtual float getPUPPIweight(float, float); 
  virtual float getSmearingFactor(float sf, float unc, float resolution, const pat::Jet & jet, const edm::View<reco::GenJet> & genJets, int variation, float drMax, float relResMax, bool usePuppiPt);
  // we need all these 3 overloaded methods as we use different 4-vector classes
  // math::XYZTLorentzVector is really a ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >
  virtual void saveDibosonMass(const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > & leptonicV_p4, const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > & hadronicV_p4, double & massVar);
  virtual void saveDibosonMass(math::XYZTLorentzVector & leptonicV_p4, const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > & hadronicV_p4, double & massVar);
  virtual void saveDibosonMass(math::XYZTLorentzVector & leptonicV_p4, math::XYZTLorentzVector & hadronicV_p4, double & massVar);
  virtual bool decaysHadronic(const reco::Candidate*);
 
  // ----------member data ---------------------------
  TTree* outTree_;

  Double_t costheta1 = -99.9;
  Double_t costheta2 = -99.9;
  Double_t costhetastar = -99.9;
  Double_t phi1 = -99.9;
  Double_t phi2 = -99.9;
  Double_t phi = -99.9;

  Double_t zerothsubjet_px = +99.9;
  Double_t zerothsubjet_py = +99.9;
  Double_t zerothsubjet_pz = +99.9;
  Double_t zerothsubjet_e  = +99.9;

  Double_t firstsubjet_px = +99.9;
  Double_t firstsubjet_py = +99.9;
  Double_t firstsubjet_pz = +99.9;
  Double_t firstsubjet_e  = +99.9;

  ////////////////////////////////

  Double_t zerothsubjet_pt = +99.9;
  Double_t zerothsubjet_phi = +99.9;
  Double_t zerothsubjet_eta = +99.9;
  Double_t zerothsubjet_m  = +99.9;

  Double_t firstsubjet_pt = +99.9;
  Double_t firstsubjet_phi = +99.9;
  Double_t firstsubjet_eta = +99.9;
  Double_t firstsubjet_m  = +99.9;

  Double_t gen_neutrino_pz = -99.9;
  Int_t found = 0;

  Double_t delta_neutrino0 = +99.9;
  Double_t delta_neutrino1 = +99.9;
  Double_t delta_neutrino2 = +99.9;
  Double_t delta_neutrino3 = +99.9;
  Double_t delta_neutrino4 = +99.9;
  Double_t delta_neutrino5 = +99.9;
  Double_t delta_neutrino6 = +99.9;
  Double_t delta_neutrino7 = +99.9;
  Double_t delta_neutrino8 = +99.9;
  
  //event info
  int nevent, run, lumi;
  
  //number of primary vertices
  int nPV;
  double gnPV;

  double PUweight;
  double btagWeight, btagWeight_BTagUp, btagWeight_BTagDown, btagWeight_MistagUp, btagWeight_MistagDown;
  double genWeight;
  double LeptonSF, LeptonSF_Up, LeptonSF_Down;
  double rho_;
  double totWeight, totWeight_BTagUp, totWeight_BTagDown, totWeight_MistagUp, totWeight_MistagDown, totWeight_LeptonIDUp, totWeight_LeptonIDDown;
  double VTagSF;
  double topPtSF;
  Particle Wboson_lep, METCand, Electron, Muon, Lepton;
  Particle subjet0, subjet1;
  double m_pruned;

  //Decay Info (gen level)
  DecayClass WDecayClass;
  //gen info

  double Wplus_gen_pt, Wplus_gen_eta, Wplus_gen_phi, Wplus_gen_mass;
  double Wminus_gen_pt, Wminus_gen_eta, Wminus_gen_phi, Wminus_gen_mass;
  ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> genHadronicVp4, genLeptonicVp4, genDibosonp4;
  double genPtV, genMWV;
  
  int N_had_W, N_lep_W;
  int N_had_Wgen, N_lep_Wgen;
  
  int nLooseEle, nLooseMu, nLep;

  //supercluster variables
  double sc_et, sc_eta;
  bool isEB;

  double tau1, tau2, tau3, tau21;
  
  double deltaR_LeptonWJet, deltaPhi_LeptonMet, deltaPhi_WJetMet, deltaPhi_WJetWlep;
  //systematics from MET
  double deltaPhi_LeptonMet_UnclEnUp, deltaPhi_LeptonMet_UnclEnDown, deltaPhi_LeptonMet_JECUp,deltaPhi_LeptonMet_JECDown, deltaPhi_LeptonMet_JERUp, deltaPhi_LeptonMet_JERDown, deltaPhi_LeptonMet_LeptonEnUp, deltaPhi_LeptonMet_LeptonEnDown;
  double deltaPhi_WJetMet_UnclEnUp, deltaPhi_WJetMet_UnclEnDown, deltaPhi_WJetMet_JECUp,deltaPhi_WJetMet_JECDown, deltaPhi_WJetMet_JERUp, deltaPhi_WJetMet_JERDown, deltaPhi_WJetMet_LeptonEnUp, deltaPhi_WJetMet_LeptonEnDown;
  double deltaPhi_WJetWlep_UnclEnUp, deltaPhi_WJetWlep_UnclEnDown, deltaPhi_WJetWlep_JECUp, deltaPhi_WJetWlep_JECDown, deltaPhi_WJetWlep_LeptonEnUp, deltaPhi_WJetWlep_LeptonEnDown;
  
  // Jet quantities
  int NAK8jet, NAK8jet_smearedUp, NAK8jet_smearedDown, NAK8jet_JECUp, NAK8jet_JECDown;
  int njets, njets_JERUp, njets_JERDown, njets_JECUp, njets_JECDown;
  int nbtag, nbtagMedium, nbtagLoose, nbtag_JERUp, nbtag_JERDown, nbtag_JECUp, nbtag_JECDown;

  double jet_pt, jet_eta, jet_phi, jet_mass, jet_mass_pruned, jet_mass_softdrop, jet_tau2tau1, jet_tau3tau2, jet_tau1, jet_tau2, jet_tau3;
  //PUPPI variables 
  double jet_pt_PUPPI, jet_eta_PUPPI, jet_phi_PUPPI, jet_mass_PUPPI, jet_tau1_PUPPI, jet_tau2_PUPPI, jet_tau3_PUPPI, jet_tau21_PUPPI, jet_tau32_PUPPI, jet_mass_softdrop_PUPPI, jet_tau21_DT;

  //gen info
  bool isMatched_;
  //JEC uncertainties
  double JECunc;
  double jet_pt_JECUp, jet_pt_JECDown, jet_mass_JECUp, jet_mass_JECDown, jet_mass_pruned_JECUp, jet_mass_pruned_JECDown, jet_mass_softdrop_JECUp, jet_mass_softdrop_JECDown, jet_mass_softdrop_PUPPI_JECUp, jet_mass_softdrop_PUPPI_JECDown;
  //JER uncerainties
  double jet_pt_JERUp, jet_pt_JERDown, jet_mass_JERUp, jet_mass_JERDown, jet_mass_softdrop_JERUp, jet_mass_softdrop_JERDown, jet_mass_pruned_JERUp, jet_mass_pruned_JERDown, jet_mass_softdrop_PUPPI_JERUp, jet_mass_softdrop_PUPPI_JERDown;
  //AK4 jets
  double jet2_pt, jet2_btag, jet3_pt, jet3_btag;
  double jet2_eta,jet2_phi, jet3_eta, jet3_phi;

  //additional info for AK4 jets
  std::vector<double> jetFlavours;

  std::vector<double> BgenjetStatus43_pt;
  std::vector<double> BgenjetStatus43_eta;
  std::vector<double> BgenjetStatus43_phi;
  std::vector<double> BgenjetStatus43_mass;
  std::vector<double> BgenjetStatus43_motherPDGID;

  std::vector<double> BgenjetStatus21_pt;
  std::vector<double> BgenjetStatus21_eta;
  std::vector<double> BgenjetStatus21_phi;
  std::vector<double> BgenjetStatus21_mass;
  std::vector<double> BgenjetStatus21_motherPDGID;

  //MET uncertainties
  double MET_UnclEnUp, MET_UnclEnDown, MET_JECUp, MET_JECDown, MET_JERUp, MET_JERDown,  MET_LeptonEnUp, MET_LeptonEnDown;
  //MET phi uncertainties
  double MET_Phi_UnclEnUp, MET_Phi_UnclEnDown, MET_Phi_JECUp, MET_Phi_JECDown, MET_Phi_JERUp, MET_Phi_JERDown, MET_Phi_LeptonEnUp, MET_Phi_LeptonEnDown;
  
  //m_lvj
  double m_lvj, m_lvj_SD;
  //m_lvj systematics
  double m_lvj_UnclEnUp, m_lvj_UnclEnDown, m_lvj_JECUp, m_lvj_JECDown, m_lvj_LeptonEnUp, m_lvj_LeptonEnDown, m_lvj_LeptonResUp, m_lvj_LeptonResDown, m_lvj_JERUp, m_lvj_JERDown;
  double m_lvj_SD_UnclEnUp, m_lvj_SD_UnclEnDown, m_lvj_SD_JECUp, m_lvj_SD_JECDown, m_lvj_SD_LeptonEnUp, m_lvj_SD_LeptonEnDown, m_lvj_SD_LeptonResUp, m_lvj_SD_LeptonResDown, m_lvj_SD_JERUp, m_lvj_SD_JERDown;

  double refXsec;
  //aTGC weights
  std::vector<double> aTGCWeights;
  double aTGCWeightUnitConv;

  int NominalPDF;

  std::vector<double> PDFWeights;
  std::vector<double> ScaleWeights;
  bool bit_HLT_Ele_105, bit_HLT_Ele_27, bit_HLT_Ele_45, bit_HLT_Ele_115, bit_HLT_Ele_30, bit_HLT_Ele_50_Jet_165, bit_BOTH_115_165;
  double triggerWeightHLTEle27NoER;
  
  
  //Defining Tokens
  edm::EDGetTokenT<std::vector< PileupSummaryInfo > > PUInfoToken_;
  edm::EDGetTokenT<edm::View<pat::MET> > metToken_;
  edm::EDGetTokenT<edm::View<reco::Candidate>> leptonicVToken_;
  edm::EDGetTokenT<edm::View<reco::Candidate>> genParticlesToken_;
  edm::EDGetTokenT<edm::View<pat::Jet>> fatJetsToken_, fatJetsSmearedUpToken_, fatJetsSmearedDownToken_;
  edm::EDGetTokenT<edm::View<reco::GenJet>> genJetsAK8Token_;
  edm::EDGetTokenT<edm::View<pat::Jet>> AK4JetsToken_, AK4JetsSmearedUpToken_, AK4JetsSmearedDownToken_, AK4JetsShiftedUpToken_, AK4JetsShiftedDownToken_;
  edm::EDGetTokenT<edm::View<reco::Vertex> > vertexToken_;
  edm::EDGetTokenT<edm::View<reco::Candidate>> looseMuToken_;
  edm::EDGetTokenT<edm::View<reco::Candidate>> looseEleToken_;
  edm::EDGetTokenT<edm::View<reco::Candidate>> leptonsToken_;
  edm::EDGetTokenT<edm::View<pat::Muon>> muonsToken_;
  edm::EDGetTokenT<edm::TriggerResults> TriggerResultsToken;
  bool isMC;
  edm::EDGetTokenT<double> rhoToken_;
  bool isSignal;
  std::string channel;
  edm::LumiReWeighting  LumiWeights_;
  edm::EDGetTokenT<GenEventInfoProduct> genInfoToken;
  edm::EDGetTokenT<LHEEventProduct> LHEEventProductTokenExternal;
  edm::EDGetTokenT<LHERunInfoProduct> lheProducerToken;
  SystematicsHelper SystematicsHelper_;
  MuonScaleFactor MuonScaleFactor_;
  ElectronScaleFactor ElectronScaleFactor_;
  JetResolutionSmearer<pat::Jet>JetResolutionSmearer_;
  BTagHelper BTagHelper_;
  // For PUPPI Softdrop Mass Correction
  TF1 *puppisd_corrGEN, *puppisd_corrRECO_cen, *puppisd_corrRECO_for;
  JetCorrectionUncertainty * jecUnc;
  double bTagDiscrCut, bTagDiscrCutMedium, bTagDiscrCutLoose;
};

//
// constructors and destructor
//
TreeMaker::TreeMaker(const edm::ParameterSet& iConfig):
  metToken_(consumes<edm::View<pat::MET>>(iConfig.getParameter<edm::InputTag>("metSrc"))),
  leptonicVToken_(consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>("leptonicVSrc"))),
  genParticlesToken_(mayConsume<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>("genSrc"))),
  fatJetsToken_(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("fatJetSrc"))),
  genJetsAK8Token_(mayConsume<edm::View<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("genJetsAK8Src"))),
  AK4JetsToken_(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("AK4JetSrc"))),
  vertexToken_(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertexSrc"))),
  looseMuToken_(consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>("looseMuSrc"))),
  looseEleToken_(consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>("looseEleSrc"))),
  leptonsToken_(consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>("leptonSrc"))),
  muonsToken_(mayConsume<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("leptonSrc"))),
  isMC(iConfig.getParameter<bool>("isMC")),
  rhoToken_(consumes<double> (iConfig.getParameter<edm::InputTag>("rho"))),
  isSignal(iConfig.getParameter<bool>("isSignal")),
  channel(iConfig.getParameter<std::string>("channel")),
  SystematicsHelper_(SystematicsHelper()),
  MuonScaleFactor_(
    "aTGCsAnalysis/TreeMaker/data/Muon_TrackingSF_RunBtoH.root",
    "aTGCsAnalysis/TreeMaker/data/Muon_IDSF_average_RunBtoH.root",
    "aTGCsAnalysis/TreeMaker/data/Muon_IsoSF_average_RunBtoH.root",
    "aTGCsAnalysis/TreeMaker/data/Muon_SingleLeptonTriggerSF_average_RunBtoH.root",
    "HighPtID",
    "tkLooseISO"
  ),
  ElectronScaleFactor_(
    "aTGCsAnalysis/TreeMaker/data/HEEP_SF.root",
    "aTGCsAnalysis/TreeMaker/data/egammaEffi.txt_EGM2D.root"
  ),
  JetResolutionSmearer_(iConfig.getParameter<bool>("isMC")),
  BTagHelper_(iConfig.getParameter<std::string>("BtagEffFile"), iConfig.getParameter<double>("BtagDiscrCut")),
  bTagDiscrCut(iConfig.getParameter<double>("BtagDiscrCut")),
  bTagDiscrCutMedium(iConfig.getParameter<double>("BtagDiscrCutMedium")),
  bTagDiscrCutLoose(iConfig.getParameter<double>("BtagDiscrCutLoose"))
{

  if ((channel != "mu") && (channel != "el")) {throw cms::Exception("InvalidValue") << "Invalid value for channel parameter, should be mu or el." << std::endl;}

  // For PUPPI Correction
  edm::FileInPath puppiCorr("aTGCsAnalysis/TreeMaker/data/puppiCorrSummer16.root");
  TFile* file = TFile::Open( puppiCorr.fullPath().c_str(),"READ");
  puppisd_corrGEN      = (TF1*)file->Get("puppiJECcorr_gen");
  puppisd_corrRECO_cen = (TF1*)file->Get("puppiJECcorr_reco_0eta1v3");
  puppisd_corrRECO_for = (TF1*)file->Get("puppiJECcorr_reco_1v3eta2v5");

  //loading PU and generator information for MC
   if (isMC)
     {
       SystematicsHelper_  = SystematicsHelper(channel, consumesCollector());

       // Use the AK8 jet label as basis
       fatJetsSmearedUpToken_ = consumes<edm::View<pat::Jet>>(edm::InputTag(iConfig.getParameter<edm::InputTag>("fatJetSrc").label() + "SmearedUp"));
       fatJetsSmearedDownToken_ = consumes<edm::View<pat::Jet>>(edm::InputTag(iConfig.getParameter<edm::InputTag>("fatJetSrc").label() + "SmearedDown"));

       // Use the AK4 jet label as basis
       AK4JetsSmearedUpToken_ = consumes<edm::View<pat::Jet>>(edm::InputTag(iConfig.getParameter<edm::InputTag>("AK4JetSrc").label() + "SmearedUp"));
       AK4JetsSmearedDownToken_ = consumes<edm::View<pat::Jet>>(edm::InputTag(iConfig.getParameter<edm::InputTag>("AK4JetSrc").label() + "SmearedDown"));

       AK4JetsShiftedUpToken_ = consumes<edm::View<pat::Jet>>(edm::InputTag(iConfig.getParameter<edm::InputTag>("AK4JetSrc").label() + "ShiftedUp"));
       AK4JetsShiftedDownToken_ = consumes<edm::View<pat::Jet>>(edm::InputTag(iConfig.getParameter<edm::InputTag>("AK4JetSrc").label() + "ShiftedDown"));

       PUInfoToken_ = consumes<std::vector< PileupSummaryInfo > >(iConfig.getParameter<edm::InputTag>("PUInfo"));

       //PU-reweighting
       LumiWeights_ = edm::LumiReWeighting(MC_dist(), data_dist());

       genInfoToken = mayConsume<GenEventInfoProduct> (iConfig.getParameter<edm::InputTag>( "genInfo" ) );
       LHEEventProductTokenExternal = mayConsume<LHEEventProduct> (iConfig.getParameter<edm::InputTag>( "LHEEventProductSrcExternal" ) );
       if(!isSignal)lheProducerToken = consumes< LHERunInfoProduct, edm::InRun >(edm::InputTag("externalLHEProducer"));
       else lheProducerToken = consumes< LHERunInfoProduct, edm::InRun >(edm::InputTag("externalLHEProducer"));
       VTagSF = iConfig.getParameter<double>("VTagSF");

     }

   jecUnc = nullptr;

   if (channel == "el") TriggerResultsToken = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggers"));

  //now do what ever initialization is needed
  edm::Service<TFileService> fs;
  outTree_ = fs->make<TTree>("BasicTree","BasicTree");

  //event info
  outTree_->Branch("event",	      &nevent,    	  "event/I"           );
  outTree_->Branch("lumi", 	      &lumi,   		  "lumi/I"  		);
  outTree_->Branch("run",	      &run,		  "run/I"  	       );

  outTree_->Branch("costheta1", &costheta1, "costheta1/D");
  outTree_->Branch("costheta2", &costheta2, "costheta2/D");
  outTree_->Branch("costhetastar", &costhetastar, "costhetastar/D");
  outTree_->Branch("phi1", &phi1, "phi1/D");
  outTree_->Branch("phi2", &phi2, "phi2/D");
  outTree_->Branch("phi", &phi, "phi/D");
  //outTree_->Branch("", &, "/D");
  outTree_->Branch("sj0_pt", &subjet0.pt, "sj0_pt/D");
  outTree_->Branch("sj0_phi", &subjet0.phi, "sj0_phi/D");
  outTree_->Branch("sj0_eta", &subjet0.eta, "sj0_eta/D");
  outTree_->Branch("sj0_mass", &subjet0.mass, "sj0_mass/D");

  outTree_->Branch("sj1_pt", &subjet1.pt, "sj1_pt/D");
  outTree_->Branch("sj1_phi", &subjet1.phi, "sj1_phi/D");
  outTree_->Branch("sj1_eta", &subjet1.eta, "sj1_eta/D");
  outTree_->Branch("sj1_mass", &subjet1.mass, "sj1_mass/D");

  outTree_->Branch("zeroth_subjet_px", &zerothsubjet_px, "zeroth_subjet_px/D");
  outTree_->Branch("zeroth_subjet_py", &zerothsubjet_py, "zeroth_subjet_py/D");
  outTree_->Branch("zeroth_subjet_pz", &zerothsubjet_pz, "zeroth_subjet_pz/D");
  outTree_->Branch("zeroth_subjet_e",  &zerothsubjet_e,  "zeroth_subjet_e/D");

  outTree_->Branch("first_subjet_px", &firstsubjet_px, "first_subjet_px/D");
  outTree_->Branch("first_subjet_py", &firstsubjet_py, "first_subjet_py/D");
  outTree_->Branch("first_subjet_pz", &firstsubjet_pz, "first_subjet_pz/D");
  outTree_->Branch("first_subjet_e",  &firstsubjet_e,  "first_subjet_e/D");


  ///////////////////////////////////////////////////////

  outTree_->Branch("zeroth_subjet_pt", &zerothsubjet_pt, "zeroth_subjet_pt/D");
  outTree_->Branch("zeroth_subjet_phi", &zerothsubjet_phi, "zeroth_subjet_phi/D");
  outTree_->Branch("zeroth_subjet_eta", &zerothsubjet_eta, "zeroth_subjet_eta/D");
  outTree_->Branch("zeroth_subjet_m",  &zerothsubjet_m,  "zeroth_subjet_m/D");

  outTree_->Branch("first_subjet_pt", &firstsubjet_pt, "first_subjet_pt/D");
  outTree_->Branch("first_subjet_phi", &firstsubjet_phi, "first_subjet_phi/D");
  outTree_->Branch("first_subjet_eta", &firstsubjet_eta, "first_subjet_eta/D");
  outTree_->Branch("first_subjet_m",  &firstsubjet_m,  "first_subjet_m/D");

  outTree_->Branch("gen_neutrino_pz", &gen_neutrino_pz, "gen_neutrino_pz/D");
  outTree_->Branch("found", &found, "found/I");

  outTree_->Branch("delta_neutrino0", &delta_neutrino0, "delta_neutrino0/D");
  outTree_->Branch("delta_neutrino1", &delta_neutrino1, "delta_neutrino1/D");
  outTree_->Branch("delta_neutrino2", &delta_neutrino2, "delta_neutrino2/D");
  outTree_->Branch("delta_neutrino3", &delta_neutrino3, "delta_neutrino3/D");
  outTree_->Branch("delta_neutrino4", &delta_neutrino4, "delta_neutrino4/D");
  outTree_->Branch("delta_neutrino5", &delta_neutrino5, "delta_neutrino5/D");
  outTree_->Branch("delta_neutrino6", &delta_neutrino6, "delta_neutrino6/D");
  outTree_->Branch("delta_neutrino7", &delta_neutrino7, "delta_neutrino7/D");
  outTree_->Branch("delta_neutrino8", &delta_neutrino8, "delta_neutrino8/D");

  outTree_->Branch("number_of_subjets", &number_of_subjets, "number_of_subjets/I");
  outTree_->Branch("imaginary_neutrino", &imaginary_neutrino, "imaginary_neutrino/I");
  //outTree_->Branch("lepton_xerror", &lepton_xerror, "lepton_xerror/D");
  //outTree_->Branch("lepton_yerror", &lepton_yerror, "lepton_yerror/D");
  //outTree_->Branch("lepton_zerror", &lepton_zerror, "lepton_zerror/D");
  //outTree_->Branch("met_xerror", &met_xerror, "met_xerror/D");
  //outTree_->Branch("met_yerror", &met_yerror, "met_yerror/D");
   
  //number of primary vertices
  outTree_->Branch("nPV",	      &nPV,		  "nPV/I"  	       );
  outTree_->Branch("rho",       &rho_,     "rho/D"          );
  //PUweight
  if (isMC)
    {
      outTree_->Branch("gnPV", &gnPV, "gnPV/D");

      outTree_->Branch("puweight",       &PUweight,     "puweight/D"          );
      outTree_->Branch("LeptonSF",       &LeptonSF,     "LeptonSF/D"          );
      outTree_->Branch("LeptonSF_Up",       &LeptonSF_Up,     "LeptonSF_Up/D"          );
      outTree_->Branch("LeptonSF_Down",       &LeptonSF_Down,     "LeptonSF_Down/D"          );
      outTree_->Branch("genweight",       &genWeight,     "genweight/D"          );
      outTree_->Branch("btagWeight",       &btagWeight,     "btagWeight/D"          );
      if(channel == "el") outTree_->Branch("triggerWeightHLTEle27NoER",       &triggerWeightHLTEle27NoER,     "triggerWeightHLTEle27NoER/D"          );
      outTree_->Branch("btagWeight_BTagUp",       &btagWeight_BTagUp,     "btagWeight_BTagUp/D"          );
      outTree_->Branch("btagWeight_BTagDown",       &btagWeight_BTagDown,     "btagWeight_BTagDown/D"          );
      outTree_->Branch("btagWeight_MistagUp",       &btagWeight_MistagUp,     "btagWeight_MistagUp/D"          );
      outTree_->Branch("btagWeight_MistagDown",       &btagWeight_MistagDown,     "btagWeight_MistagDown/D"          );
      outTree_->Branch("topPtSF",       &topPtSF,     "topPtSF/D"          );
      //total weights: central and systematics
      outTree_->Branch("totWeight",       &totWeight,     "totWeight/D"          );
     
      outTree_->Branch("totWeight_BTagUp",       &totWeight_BTagUp,     "totWeight_BTagUp/D"          );
      outTree_->Branch("totWeight_BTagDown",       &totWeight_BTagDown,     "totWeight_BTagDown/D"          );
     
      outTree_->Branch("totWeight_MistagUp",       &totWeight_MistagUp,     "totWeight_MistagUp/D"          );
      outTree_->Branch("totWeight_MistagDown",       &totWeight_MistagDown,     "totWeight_MistagDown/D"          );

      outTree_->Branch("totWeight_LeptonIDUp",       &totWeight_LeptonIDUp,     "totWeight_LeptonIDUp/D"          );
      outTree_->Branch("totWeight_LeptonIDDown",       &totWeight_LeptonIDDown,     "totWeight_LeptonIDDown/D"          );


      //PDF and scale weights: systematics
      outTree_->Branch("PDFWeights","std::vector<double>",&PDFWeights);
      outTree_->Branch("ScaleWeights","std::vector<double>",&ScaleWeights);
      //generator info about the decay of WW
      outTree_->Branch("WDecayClass",     &WDecayClass,    "WDecayClass/I"      );
      //generator W info
      outTree_->Branch("Wplus_gen_pt",     &Wplus_gen_pt,    "Wplus_gen_pt/D"      );
      outTree_->Branch("Wplus_gen_eta",     &Wplus_gen_eta,    "Wplus_gen_eta/D"      );
      outTree_->Branch("Wplus_gen_phi",     &Wplus_gen_phi,    "Wplus_gen_phi/D"      );
      outTree_->Branch("Wplus_gen_mass",     &Wplus_gen_mass,    "Wplus_gen_mass/D"      );

      outTree_->Branch("Wminus_gen_pt",     &Wminus_gen_pt,    "Wminus_gen_pt/D"      );
      outTree_->Branch("Wminus_gen_eta",     &Wminus_gen_eta,    "Wminus_gen_eta/D"      );
      outTree_->Branch("Wminus_gen_phi",     &Wminus_gen_phi,    "Wminus_gen_phi/D"      );
      outTree_->Branch("Wminus_gen_mass",     &Wminus_gen_mass,    "Wminus_gen_mass/D"      );

      outTree_->Branch("genPtV", &genPtV, "genPtV/D");
      outTree_->Branch("genMWV", &genMWV, "genMWV/D");
    };
  
  if (channel == "el")
    {
      outTree_->Branch("bit_HLT_Ele_105",       &bit_HLT_Ele_105,     "bit_HLT_Ele_105/B"          );
      outTree_->Branch("bit_HLT_Ele_27",       &bit_HLT_Ele_27,     "bit_HLT_Ele_27/B"          );
      outTree_->Branch("bit_HLT_Ele_45",       &bit_HLT_Ele_45,     "bit_HLT_Ele_45/B"          );
      outTree_->Branch("bit_HLT_Ele_115",       &bit_HLT_Ele_115,     "bit_HLT_Ele_115/B"          );
      outTree_->Branch("bit_HLT_Ele_30",       &bit_HLT_Ele_30,     "bit_HLT_Ele_30/B"          );
      outTree_->Branch("bit_HLT_Ele_50_Jet_165",       &bit_HLT_Ele_50_Jet_165,     "bit_HLT_Ele_50_Jet_165/B"          );
      outTree_->Branch("bit_BOTH_115_165",       &bit_BOTH_115_165,     "bit_BOTH_115_165/B"          );
    }
  
  //number of loose leptons
  outTree_->Branch("nLooseEle",      &nLooseEle, 	  "nLooseEle/I"       );
  outTree_->Branch("nLooseMu",      &nLooseMu, 	  "nLooseMu/I"       );
  
  // number of leptons for particular mu/ele channel (should be tight lepton)
  outTree_->Branch("nLep",            &nLep, 	      "nLep/I"       );
  
  //leptons (tight, ele or mu depends on the channel)
  outTree_->Branch("l_pt",	      &Lepton.pt,      "l_pt/D"         	);
  outTree_->Branch("l_eta",	      &Lepton.eta,     "l_eta/D"        	);
  outTree_->Branch("l_phi",	      &Lepton.phi,     "l_phi/D"        	);
  //lepton uncertainties
  if (isMC)
    {
      if(channel == "el")
	{
	  outTree_->Branch("sc_eta",       &sc_eta,     "sc_eta/D"          );  
	  outTree_->Branch("sc_et",       &sc_et,     "sc_et/D"          );  
	  outTree_->Branch("isEB",       &isEB,     "isEB/B"          );  
	}
      outTree_->Branch("l_pt_LeptonEnUp",       &Lepton.pt_LeptonEnUp,     "l_pt_LeptonEnUp/D"          );
      outTree_->Branch("l_pt_LeptonEnDown",     &Lepton.pt_LeptonEnDown,   "l_pt_LeptonEnDown/D"          );
      outTree_->Branch("l_pt_LeptonResUp",      &Lepton.pt_LeptonResUp,    "l_pt_LeptonResUp/D"          );
      outTree_->Branch("l_pt_LeptonResDown",    &Lepton.pt_LeptonResDown,  "l_pt_LeptonResDown/D"          );
    }
  //W leptonic observables
  outTree_->Branch("W_pt",	      &Wboson_lep.pt,     "W_pt/D"         );
  outTree_->Branch("W_eta",	      &Wboson_lep.eta,    "W_eta/D"        );
  outTree_->Branch("W_phi",	      &Wboson_lep.phi,    "W_phi/D"        );
  outTree_->Branch("W_mass",      &Wboson_lep.mass,   "W_mass/D"       );
  outTree_->Branch("W_mt",	      &Wboson_lep.mt,     "W_mt/D"         );

  if (isMC)
    {
      //pt
      outTree_->Branch("W_pt_LeptonEnUp",        &Wboson_lep.pt_LeptonEnUp,     "W_pt_LeptonEnUp/D"         );
      outTree_->Branch("W_pt_LeptonEnDown",      &Wboson_lep.pt_LeptonEnDown,   "W_pt_LeptonEnDown/D"       );
    
      outTree_->Branch("W_pt_LeptonResUp",        &Wboson_lep.pt_LeptonResUp,     "W_pt_LeptonResUp/D"         );
      outTree_->Branch("W_pt_LeptonResDown",      &Wboson_lep.pt_LeptonResDown,   "W_pt_LeptonResDown/D"       );

      outTree_->Branch("W_pt_JECUp",        &Wboson_lep.pt_JECUp,     "W_pt_JECUp/D"         );
      outTree_->Branch("W_pt_JECDown",      &Wboson_lep.pt_JECDown,   "W_pt_JECDown/D"       );

      outTree_->Branch("W_pt_UnclEnUp",        &Wboson_lep.pt_UnclEnUp,     "W_pt_UnclEnUp/D"         );
      outTree_->Branch("W_pt_UnclEnDown",      &Wboson_lep.pt_UnclEnDown,   "W_pt_UnclEnDown/D"       );

      //mass 
      outTree_->Branch("W_mass_LeptonEnUp",        &Wboson_lep.mass_LeptonEnUp,     "W_mass_LeptonEnUp/D"         );
      outTree_->Branch("W_mass_LeptonEnDown",      &Wboson_lep.mass_LeptonEnDown,   "W_mass_LeptonEnDown/D"       );
    
      outTree_->Branch("W_mass_LeptonResUp",        &Wboson_lep.mass_LeptonResUp,     "W_mass_LeptonResUp/D"         );
      outTree_->Branch("W_mass_LeptonResDown",      &Wboson_lep.mass_LeptonResDown,   "W_mass_LeptonResDown/D"       );

      outTree_->Branch("W_mass_JECUp",        &Wboson_lep.mass_JECUp,     "W_mass_JECUp/D"         );
      outTree_->Branch("W_mass_JECDown",      &Wboson_lep.mass_JECDown,   "W_mass_JECDown/D"       );

      outTree_->Branch("W_mass_UnclEnUp",        &Wboson_lep.mass_UnclEnUp,     "W_mass_UnclEnUp/D"         );
      outTree_->Branch("W_mass_UnclEnDown",      &Wboson_lep.mass_UnclEnDown,   "W_mass_UnclEnDown/D"       );

      //mt 
      outTree_->Branch("W_mt_LeptonEnUp",        &Wboson_lep.mt_LeptonEnUp,     "W_mt_LeptonEnUp/D"         );
      outTree_->Branch("W_mt_LeptonEnDown",      &Wboson_lep.mt_LeptonEnDown,   "W_mt_LeptonEnDown/D"       );
    
      outTree_->Branch("W_mt_LeptonResUp",        &Wboson_lep.mt_LeptonResUp,     "W_mt_LeptonResUp/D"         );
      outTree_->Branch("W_mt_LeptonResDown",      &Wboson_lep.mt_LeptonResDown,   "W_mt_LeptonResDown/D"       );

      outTree_->Branch("W_mt_JECUp",        &Wboson_lep.mt_JECUp,     "W_mt_JECUp/D"         );
      outTree_->Branch("W_mt_JECDown",      &Wboson_lep.mt_JECDown,   "W_mt_JECDown/D"       );

      outTree_->Branch("W_mt_UnclEnUp",        &Wboson_lep.mt_UnclEnUp,     "W_mt_UnclEnUp/D"         );
      outTree_->Branch("W_mt_UnclEnDown",      &Wboson_lep.mt_UnclEnDown,   "W_mt_UnclEnDown/D"       );   
    } 

  
  outTree_->Branch("charge_W_lep",    &Wboson_lep.charge, "charge_W_lep/D"     );
    
  outTree_->Branch("N_lep_W",	      &N_lep_W,		  "N_lep_W/I"          );
  
  outTree_->Branch("N_had_W_gen",     &N_had_Wgen,	  "N_had_W_gen/I"      );
  outTree_->Branch("N_lep_W_gen",     &N_lep_Wgen, 	  "N_lep_W_gen/I"      );

  //MET observables  
  outTree_->Branch("pfMET", 	      &METCand.pt, 	  "pfMET/D"              );
  outTree_->Branch("pfMETPhi",	      &METCand.phi, 	  "pfMETPhi/D"          );
  

  if (isMC)
    {
      //MET uncertainties
      //UnclEn
      outTree_->Branch("pfMET_UnclEnUp",         &MET_UnclEnUp,    "pfMET_UnclEnUp/D"              );
      outTree_->Branch("pfMET_UnclEnDown",         &MET_UnclEnDown,    "pfMET_UnclEnDown/D"              );
      //JER
      outTree_->Branch("pfMET_JERUp",         &MET_JERUp,    "pfMET_JERUp/D"              );
      outTree_->Branch("pfMET_JERDown",         &MET_JERDown,    "pfMET_JERDown/D"              );
      //JEC
      outTree_->Branch("pfMET_JECUp",         &MET_JECUp,    "pfMET_JECUp/D"              );
      outTree_->Branch("pfMET_JECDown",         &MET_JECDown,    "pfMET_JECDown/D"              );
      //Lepton energy energy scale
      outTree_->Branch("pfMET_LeptonEnUp",         &MET_LeptonEnUp,    "pfMET_LeptonEnUp/D"              );
      outTree_->Branch("pfMET_LeptonEnDown",         &MET_LeptonEnDown,    "pfMET_LeptonEnDown/D"              );
   
      //MET phi uncertainties
      //UnclEn
      outTree_->Branch("pfMETPhi_UnclEnUp",         &MET_Phi_UnclEnUp,    "pfMETPhi_UnclEnUp/D"              );
      outTree_->Branch("pfMETPhi_UnclEnDown",         &MET_Phi_UnclEnDown,    "pfMETPhi_UnclEnDown/D"              );
      //JER
      outTree_->Branch("pfMETPhi_JERUp",         &MET_Phi_JERUp,    "pfMETPhi_JERUp/D"              );
      outTree_->Branch("pfMETPhi_JERDown",         &MET_Phi_JERDown,    "pfMETPhi_JERDown/D"              );
      //JEC
      outTree_->Branch("pfMETPhi_JECUp",         &MET_Phi_JECUp,    "pfMETPhi_JECUp/D"              );
      outTree_->Branch("pfMETPhi_JECDown",         &MET_Phi_JECDown,    "pfMETPhi_JECDown/D"              );
      //Lepton energy scale
      outTree_->Branch("pfMETPhi_LeptonEnUp",         &MET_Phi_LeptonEnUp,    "pfMETPhi_LeptonEnUp/D"              );
      outTree_->Branch("pfMETPhi_LeptonEnDown",         &MET_Phi_LeptonEnDown,    "pfMETPhi_LeptonEnDown/D"              );
    }
  /// Other observables
  outTree_->Branch("deltaR_LeptonWJet",  &deltaR_LeptonWJet,	  "deltaR_LeptonWJet/D"   );
  outTree_->Branch("deltaPhi_LeptonMet", &deltaPhi_LeptonMet,	  "deltaPhi_LeptonMet/D"  );
  outTree_->Branch("deltaPhi_WJetMet",&deltaPhi_WJetMet,  "deltaPhi_WJetMet/D" );
  outTree_->Branch("deltaPhi_WJetWlep",&deltaPhi_WJetWlep,  "deltaPhi_WJetWlep/D" );
  if (isMC)
    {
      //UnclEn
      outTree_->Branch("deltaPhi_LeptonMet_UnclEnUp", &deltaPhi_LeptonMet_UnclEnUp,   "deltaPhi_LeptonMet_UnclEnUp/D"  );
      outTree_->Branch("deltaPhi_LeptonMet_UnclEnDown", &deltaPhi_LeptonMet_UnclEnDown,   "deltaPhi_LeptonMet_UnclEnDown/D"  );
      //JEC
      outTree_->Branch("deltaPhi_LeptonMet_JECUp", &deltaPhi_LeptonMet_JECUp,   "deltaPhi_LeptonMet_JECUp/D"  );
      outTree_->Branch("deltaPhi_LeptonMet_JECDown", &deltaPhi_LeptonMet_JECDown,   "deltaPhi_LeptonMet_JECDown/D"  );
      //Lepton energy scale
      outTree_->Branch("deltaPhi_LeptonMet_LeptonEnUp", &deltaPhi_LeptonMet_LeptonEnUp,   "deltaPhi_LeptonMet_LeptonEnUp/D"  );
      outTree_->Branch("deltaPhi_LeptonMet_LeptonEnDown", &deltaPhi_LeptonMet_LeptonEnDown,   "deltaPhi_LeptonMet_LeptonEnDown/D"  );
      //JER
      outTree_->Branch("deltaPhi_LeptonMet_JERUp", &deltaPhi_LeptonMet_JERUp,   "deltaPhi_LeptonMet_JERUp/D"  );
      outTree_->Branch("deltaPhi_LeptonMet_JERDown", &deltaPhi_LeptonMet_JERDown,   "deltaPhi_LeptonMet_LeptonEnDown/D"  );
      /////////////////////
      //UnclEn
      outTree_->Branch("deltaPhi_WJetMet_UnclEnUp", &deltaPhi_WJetMet_UnclEnUp,   "deltaPhi_WJetMet_UnclEnUp/D"  );
      outTree_->Branch("deltaPhi_WJetMet_UnclEnDown", &deltaPhi_WJetMet_UnclEnDown,   "deltaPhi_WJetMet_UnclEnDown/D"  );
      //JEC
      outTree_->Branch("deltaPhi_WJetMet_JECUp", &deltaPhi_WJetMet_JECUp,   "deltaPhi_WJetMet_JECUp/D"  );
      outTree_->Branch("deltaPhi_WJetMet_JECDown", &deltaPhi_WJetMet_JECDown,   "deltaPhi_WJetMet_JECDown/D"  );
      //Lepton energy scale
      outTree_->Branch("deltaPhi_WJetMet_LeptonEnUp", &deltaPhi_WJetMet_LeptonEnUp,   "deltaPhi_WJetMet_LeptonEnUp/D"  );
      outTree_->Branch("deltaPhi_WJetMet_LeptonEnDown", &deltaPhi_WJetMet_LeptonEnDown,   "deltaPhi_WJetMet_LeptonEnDown/D"  );
      //JER
      outTree_->Branch("deltaPhi_WJetMet_JERUp", &deltaPhi_WJetMet_JERUp,   "deltaPhi_WJetMet_JERUp/D"  );
      outTree_->Branch("deltaPhi_WJetMet_JERDown", &deltaPhi_WJetMet_JERDown,   "deltaPhi_WJetMet_LeptonEnDown/D"  );
      /////////////////////
      //UnclEn
      outTree_->Branch("deltaPhi_WJetWlep_UnclEnUp", &deltaPhi_WJetWlep_UnclEnUp,   "deltaPhi_WJetWlep_UnclEnUp/D"  );
      outTree_->Branch("deltaPhi_WJetWlep_UnclEnDown", &deltaPhi_WJetWlep_UnclEnDown,   "deltaPhi_WJetWlep_UnclEnDown/D"  );
      //JEC
      outTree_->Branch("deltaPhi_WJetWlep_JECUp", &deltaPhi_WJetWlep_JECUp,   "deltaPhi_WJetWlep_JECUp/D"  );
      outTree_->Branch("deltaPhi_WJetWlep_JECDown", &deltaPhi_WJetWlep_JECDown,   "deltaPhi_WJetWlep_JECDown/D"  );
      //Lepton energy scale
      outTree_->Branch("deltaPhi_WJetWlep_LeptonEnUp", &deltaPhi_WJetWlep_LeptonEnUp,   "deltaPhi_WJetWlep_LeptonEnUp/D"  );
      outTree_->Branch("deltaPhi_WJetWlep_LeptonEnDown", &deltaPhi_WJetWlep_LeptonEnDown,   "deltaPhi_WJetWlep_LeptonEnDown/D"  );
    }
  //Jet observables
  outTree_->Branch("NAK8jet",            &NAK8jet,              "NAK8jet/I"   );
  outTree_->Branch("jet_pt",  	      &jet_pt,	  	  "jet_pt/D"   );
  outTree_->Branch("jet_eta",  	      &jet_eta,	  	  "jet_eta/D"   );
  outTree_->Branch("jet_phi",  	      &jet_phi,	  	  "jet_phi/D"   );
  outTree_->Branch("jet_mass",         &jet_mass,       "jet_mass/D"   );
  outTree_->Branch("Mjpruned", &jet_mass_pruned,	  "Mjpruned/D"   );
  outTree_->Branch("jet_mass_softdrop",&jet_mass_softdrop,"jet_mass_softdrop/D"   );
  outTree_->Branch("jet_tau2tau1",    &jet_tau2tau1,	  "jet_tau2tau1/D"   );
  outTree_->Branch("jet_tau3tau2",    &jet_tau3tau2,    "jet_tau3tau2/D"   );
  outTree_->Branch("jet_tau1",    &jet_tau1,    "jet_tau1/D"   );
  outTree_->Branch("jet_tau2",    &jet_tau2,    "jet_tau2/D"   );
  outTree_->Branch("jet_tau3",    &jet_tau3,    "jet_tau3/D"   );
  //PUPPI variables
  outTree_->Branch("jet_pt_PUPPI",    &jet_pt_PUPPI,    "jet_pt_PUPPI/D"   );
  outTree_->Branch("jet_eta_PUPPI",    &jet_eta_PUPPI,    "jet_eta_PUPPI/D"   );
  outTree_->Branch("jet_phi_PUPPI",    &jet_phi_PUPPI,    "jet_phi_PUPPI/D"   );
  outTree_->Branch("jet_mass_PUPPI",    &jet_mass_PUPPI,    "jet_mass_PUPPI/D"   );
  outTree_->Branch("jet_tau1_PUPPI",    &jet_tau1_PUPPI,    "jet_tau1_PUPPI/D"   );
  outTree_->Branch("jet_tau2_PUPPI",    &jet_tau2_PUPPI,    "jet_tau2_PUPPI/D"   );
  outTree_->Branch("jet_tau3_PUPPI",    &jet_tau3_PUPPI,    "jet_tau3_PUPPI/D"   );
  outTree_->Branch("jet_tau21_PUPPI",    &jet_tau21_PUPPI,    "jet_tau21_PUPPI/D"   );
  outTree_->Branch("jet_tau32_PUPPI",    &jet_tau32_PUPPI,    "jet_tau32_PUPPI/D"   );
  outTree_->Branch("jet_mass_softdrop_PUPPI",    &jet_mass_softdrop_PUPPI,    "jet_mass_softdrop_PUPPI/D"   );
  outTree_->Branch("jet_tau21_DT",    &jet_tau21_DT,    "jet_tau21_DT/D"   );
  
  if (isMC)
    {
      //JEC uncertainties
      outTree_->Branch("JECunc",    &JECunc,    "JECunc/D"   ); 
      outTree_->Branch("jet_pt_JECUp",    &jet_pt_JECUp,    "jet_pt_JECUp/D"   ); 
      outTree_->Branch("jet_pt_JECDown",    &jet_pt_JECDown,    "jet_pt_JECDown/D"   );  
      outTree_->Branch("jet_mass_JECUp",    &jet_mass_JECUp,    "jet_mass_JECUp/D"   ); 
      outTree_->Branch("jet_mass_JECDown",    &jet_mass_JECDown,    "jet_mass_JECDown/D"   );  
      //JER uncertainties
      outTree_->Branch("jet_pt_JERUp",    &jet_pt_JERUp,    "jet_pt_JERUp/D"   ); 
      outTree_->Branch("jet_pt_JERDown",    &jet_pt_JERDown,    "jet_pt_JERDown/D"   );  
      outTree_->Branch("jet_mass_JERUp",    &jet_mass_JERUp,    "jet_mass_JERUp/D"   ); 
      outTree_->Branch("jet_mass_JERDown",    &jet_mass_JERDown,    "jet_mass_JERDown/D"   );  
      //JEC uncertainties
      outTree_->Branch("Mjpruned_JECUp",    &jet_mass_pruned_JECUp,    "Mjpruned_JECUp/D"   ); 
      outTree_->Branch("Mjpruned_JECDown",    &jet_mass_pruned_JECDown,    "Mjpruned_JECDown/D"   );  
      outTree_->Branch("jet_mass_softdrop_JECUp",    &jet_mass_softdrop_JECUp,    "jet_mass_softdrop_JECUp/D"   ); 
      outTree_->Branch("jet_mass_softdrop_JECDown",    &jet_mass_softdrop_JECDown,    "jet_mass_softdrop_JECDown/D"   );  
      outTree_->Branch("jet_mass_softdrop_PUPPI_JECUp",    &jet_mass_softdrop_PUPPI_JECUp,    "jet_mass_softdrop_PUPPI_JECUp/D"   );
      outTree_->Branch("jet_mass_softdrop_PUPPI_JECDown",    &jet_mass_softdrop_PUPPI_JECDown,    "jet_mass_softdrop_PUPPI_JECDown/D"   );
      //JER uncertainties
      outTree_->Branch("Mjpruned_JERUp",    &jet_mass_pruned_JERUp,    "Mjpruned_JERUp/D"   ); 
      outTree_->Branch("Mjpruned_JERDown",    &jet_mass_pruned_JERDown,    "Mjpruned_JERDown/D"   );  
      outTree_->Branch("jet_mass_softdrop_JERUp",    &jet_mass_softdrop_JERUp,    "jet_mass_softdrop_JERUp/D"   ); 
      outTree_->Branch("jet_mass_softdrop_JERDown",    &jet_mass_softdrop_JERDown,    "jet_mass_softdrop_JERDown/D"   );  
      outTree_->Branch("jet_mass_softdrop_PUPPI_JERUp",    &jet_mass_softdrop_PUPPI_JERUp,    "jet_mass_softdrop_PUPPI_JERUp/D"   );
      outTree_->Branch("jet_mass_softdrop_PUPPI_JERDown",    &jet_mass_softdrop_PUPPI_JERDown,    "jet_mass_softdrop_PUPPI_JERDown/D"   );
      outTree_->Branch("isMatched",    &isMatched_,    "isMatched/B"   ); 
      //add info for AK4 jets
      outTree_ -> Branch("jetFlavours",  &jetFlavours); 
    
      outTree_ -> Branch("BgenjetStatus21_pt",  &BgenjetStatus21_pt); 
      outTree_ -> Branch("BgenjetStatus21_eta",  &BgenjetStatus21_eta); 
      outTree_ -> Branch("BgenjetStatus21_phi",  &BgenjetStatus21_phi); 
      outTree_ -> Branch("BgenjetStatus21_mass",  &BgenjetStatus21_mass); 
      outTree_ -> Branch("BgenjetStatus21_motherPDGID",  &BgenjetStatus21_motherPDGID); 

      outTree_ -> Branch("BgenjetStatus43_pt",  &BgenjetStatus43_pt); 
      outTree_ -> Branch("BgenjetStatus43_eta",  &BgenjetStatus43_eta); 
      outTree_ -> Branch("BgenjetStatus43_phi",  &BgenjetStatus43_phi); 
      outTree_ -> Branch("BgenjetStatus43_mass",  &BgenjetStatus43_mass); 
      outTree_ -> Branch("BgenjetStatus43_motherPDGID",  &BgenjetStatus43_motherPDGID); 
    }
  outTree_->Branch("njets",         &njets,           "njets/I"   );
  outTree_->Branch("nbtag",         &nbtag,           "nbtag/I"   );
  outTree_->Branch("nbtagMedium",         &nbtagMedium,           "nbtagMedium/I"   );
  outTree_->Branch("nbtagLoose",         &nbtagLoose,           "nbtagLoose/I"   );

  if (isMC)
    {
      outTree_->Branch("njets_JERUp",         &njets_JERUp,           "njets_JERUp/I"   );
      outTree_->Branch("njets_JERDown",       &njets_JERDown,         "njets_JERDown/I"   );

      outTree_->Branch("nbtag_JERUp",         &nbtag_JERUp,           "nbtag_JERUp/I"   );
      outTree_->Branch("nbtag_JERDown",       &nbtag_JERDown,         "nbtag_JERDown/I"   );

      outTree_->Branch("njets_JECUp",         &njets_JECUp,           "njets_JECRUp/I"   );
      outTree_->Branch("njets_JECDown",       &njets_JECDown,         "njets_JECDown/I"   );

      outTree_->Branch("nbtag_JECUp",         &nbtag_JECUp,           "nbtag_JECUp/I"   );
      outTree_->Branch("nbtag_JECDown",       &nbtag_JECDown,         "nbtag_JECDown/I"   );

    }
  
  outTree_->Branch("jet2_pt",  	      &jet2_pt,	          "jet2_pt/D"   );
  outTree_->Branch("jet2_btag",       &jet2_btag,         "jet2_btag/D"   );
  outTree_->Branch("jet3_pt",  	      &jet3_pt,	          "jet3_pt/D"   );
  outTree_->Branch("jet3_btag",	      &jet3_btag,         "jet3_btag/D"   );
  
  outTree_->Branch("MWW",       &m_lvj,         "MWW/D"   );
  outTree_->Branch("MWW_SD",	      &m_lvj_SD,         "MWW_SD/D"   );
  if (isMC)
    {
      outTree_->Branch("MWW_UnclEnUp",       &m_lvj_UnclEnUp,         "MWW_UnclEnUp/D"   );
      outTree_->Branch("MWW_UnclEnDown",       &m_lvj_UnclEnDown,         "MWW_UnclEnDown/D"   );      
      outTree_->Branch("MWW_JECUp",       &m_lvj_JECUp,         "MWW_JECUp/D"   );
      outTree_->Branch("MWW_JECDown",       &m_lvj_JECDown,         "MWW_JECDown/D"   );  
      outTree_->Branch("MWW_LeptonEnUp",       &m_lvj_LeptonEnUp,         "MWW_LeptonEnUp/D"   );
      outTree_->Branch("MWW_LeptonEnDown",       &m_lvj_LeptonEnDown,         "MWW_LeptonEnDown/D"   );      
      outTree_->Branch("MWW_LeptonResUp",       &m_lvj_LeptonResUp,         "MWW_LeptonResUp/D"   );
      outTree_->Branch("MWW_LeptonResDown",       &m_lvj_LeptonResDown,         "MWW_LeptonResDown/D"   );   
      outTree_->Branch("MWW_JERUp",       &m_lvj_JERUp,         "MWW_JERUp/D"   );
      outTree_->Branch("MWW_JERDown",       &m_lvj_JERDown,         "MWW_JERDown/D"   );       

      outTree_->Branch("MWW_SD_UnclEnUp",       &m_lvj_SD_UnclEnUp,         "MWW_SD_UnclEnUp/D"   );
      outTree_->Branch("MWW_SD_UnclEnDown",       &m_lvj_SD_UnclEnDown,         "MWW_SD_UnclEnDown/D"   );
      outTree_->Branch("MWW_SD_JECUp",       &m_lvj_SD_JECUp,         "MWW_SD_JECUp/D"   );
      outTree_->Branch("MWW_SD_JECDown",       &m_lvj_SD_JECDown,         "MWW_SD_JECDown/D"   );
      outTree_->Branch("MWW_SD_LeptonEnUp",       &m_lvj_SD_LeptonEnUp,         "MWW_SD_LeptonEnUp/D"   );
      outTree_->Branch("MWW_SD_LeptonEnDown",       &m_lvj_SD_LeptonEnDown,         "MWW_SD_LeptonEnDown/D"   );
      outTree_->Branch("MWW_SD_LeptonResUp",       &m_lvj_SD_LeptonResUp,         "MWW_SD_LeptonResUp/D"   );
      outTree_->Branch("MWW_SD_LeptonResDown",       &m_lvj_SD_LeptonResDown,         "MWW_SD_LeptonResDown/D"   );
      outTree_->Branch("MWW_SD_JERUp",       &m_lvj_SD_JERUp,         "MWW_SD_JERUp/D"   );
      outTree_->Branch("MWW_SD_JERDown",       &m_lvj_SD_JERDown,         "MWW_SD_JERDown/D"   );

#ifdef ANGLE_TESTING
      //INTERMEDIATE STEPS VARIABLES
      outTree_->Branch("leptons_in_lep_px", &leptons_in_lep_px, "leptons_in_lep_px/D");
      outTree_->Branch("leptons_in_lep_py", &leptons_in_lep_py, "leptons_in_lep_py/D");
      outTree_->Branch("leptons_in_lep_pz", &leptons_in_lep_pz, "leptons_in_lep_pz/D");
  
      outTree_->Branch("partons_in_lep_px", &partons_in_lep_px, "partons_in_lep_px/D");
      outTree_->Branch("partons_in_lep_py", &partons_in_lep_py, "partons_in_lep_py/D");
      outTree_->Branch("partons_in_lep_pz", &partons_in_lep_pz, "partons_in_lep_pz/D");
  
      outTree_->Branch("parton1_in_lep_px", &parton1_in_lep_px, "parton1_in_lep_px/D");
      outTree_->Branch("parton1_in_lep_py", &parton1_in_lep_py, "parton1_in_lep_py/D"); 
      outTree_->Branch("parton1_in_lep_pz", &parton1_in_lep_pz, "parton1_in_lep_pz/D");

      outTree_->Branch("parton2_in_lep_px", &parton2_in_lep_px, "parton2_in_lep_px/D");
      outTree_->Branch("parton2_in_lep_py", &parton2_in_lep_py, "parton2_in_lep_py/D");
      outTree_->Branch("parton2_in_lep_pz", &parton2_in_lep_pz, "parton2_in_lep_pz/D");
  
      outTree_->Branch("lepton1_in_lep_px", &lepton1_in_lep_px, "lepton1_in_lep_px/D");
      outTree_->Branch("lepton1_in_lep_py", &lepton1_in_lep_py, "lepton1_in_lep_py/D");
      outTree_->Branch("lepton1_in_lep_pz", &lepton1_in_lep_pz, "lepton1_in_lep_pz/D");
  
      outTree_->Branch("lepton1_dotted_x", &lepton1_dotted_x, "lepton1_dotted_x/D");
      outTree_->Branch("lepton1_dotted_y", &lepton1_dotted_y, "lepton1_dotted_y/D");
      outTree_->Branch("lepton1_dotted_z", &lepton1_dotted_z, "lepton1_dotted_z/D");
  
      outTree_->Branch("leptons_in_had_px", &leptons_in_had_px, "leptons_in_had_px/D");
      outTree_->Branch("leptons_in_had_py", &leptons_in_had_py, "leptons_in_had_py/D");
      outTree_->Branch("leptons_in_had_pz", &leptons_in_had_pz, "leptons_in_had_pz/D");
  
      outTree_->Branch("lepton1_in_had_px", &lepton1_in_had_px, "lepton1_in_had_px/D");
      outTree_->Branch("lepton1_in_had_py", &lepton1_in_had_py, "lepton1_in_had_py/D");
      outTree_->Branch("lepton1_in_had_pz", &lepton1_in_had_pz, "lepton1_in_had_pz/D");
  
      outTree_->Branch("lepton2_in_had_px", &lepton2_in_had_px, "lepton2_in_had_px/D");
      outTree_->Branch("lepton2_in_had_py", &lepton2_in_had_py, "lepton2_in_had_py/D");
      outTree_->Branch("lepton2_in_had_pz", &lepton2_in_had_pz, "lepton2_in_had_pz/D");
  
      outTree_->Branch("parton1_in_had_px", &parton1_in_had_px, "parton1_in_had_px/D");
      outTree_->Branch("parton1_in_had_py", &parton1_in_had_py, "parton1_in_had_py/D");
      outTree_->Branch("parton1_in_had_pz", &parton1_in_had_pz, "parton1_in_had_pz/D");
  
      outTree_->Branch("parton1_dotted_x", &parton1_dotted_x, "parton1_dotted_x/D");
      outTree_->Branch("parton1_dotted_y", &parton1_dotted_y, "parton1_dotted_y/D");
      outTree_->Branch("parton1_dotted_z", &parton1_dotted_z, "parton1_dotted_z/D");
  
      outTree_->Branch("complicated1_px", &complicated1_px, "complicated1_px/D");
      outTree_->Branch("complicated1_py", &complicated1_py, "complicated1_py/D");
      outTree_->Branch("complicated1_pz", &complicated1_pz, "complicated1_pz/D");
  
      outTree_->Branch("complicated2_px", &complicated2_px, "complicated2_px/D");
      outTree_->Branch("complicated2_py", &complicated2_py, "complicated2_py/D");
      outTree_->Branch("complicated2_pz", &complicated2_pz, "complicated2_pz/D");
  
      outTree_->Branch("lepton_sumWWframe_X", &lepton_sumWWframe_X, "lepton_sumWWframe_X/D");
      outTree_->Branch("lepton_sumWWframe_Y", &lepton_sumWWframe_Y, "lepton_sumWWframe_Y/D");
      outTree_->Branch("lepton_sumWWframe_Z", &lepton_sumWWframe_Z, "lepton_sumWWframe_Z/D");
  
      outTree_->Branch("parton_sumWWframe_X", &parton_sumWWframe_X, "parton_sumWWframe_X/D");
      outTree_->Branch("parton_sumWWframe_Y", &parton_sumWWframe_Y, "parton_sumWWframe_Y/D");
      outTree_->Branch("parton_sumWWframe_Z", &parton_sumWWframe_Z, "parton_sumWWframe_Z/D");
  
      outTree_->Branch("lepton1WWframe_X", &lepton1WWframe_X, "lepton1WWframe_X/D");
      outTree_->Branch("lepton1WWframe_Y", &lepton1WWframe_Y, "lepton1WWframe_Y/D");
      outTree_->Branch("lepton1WWframe_Z", &lepton1WWframe_Z, "lepton1WWframe_Z/D");
  
      outTree_->Branch("parton1WWframe_X", &parton1WWframe_X, "parton1WWframe_X/D");
      outTree_->Branch("parton1WWframe_Y", &parton1WWframe_Y, "parton1WWframe_Y/D");
      outTree_->Branch("parton1WWframe_Z", &parton1WWframe_Z, "parton1WWframe_Z/D");

      outTree_->Branch("boostWWframe_X", &boostWWframe_X, "boostWWframe_X/D");
      outTree_->Branch("boostWWframe_Y", &boostWWframe_Y, "boostWWframe_Y/D");
      outTree_->Branch("boostWWframe_Z", &boostWWframe_Z, "boostWWframe_Z/D");
  
      outTree_->Branch("boostWlep_X", &boostWlep_X, "boostWlep_X/D");
      outTree_->Branch("boostWlep_Y", &boostWlep_Y, "boostWlep_Y/D");
      outTree_->Branch("boostWlep_Z", &boostWlep_Z, "boostWlep_Z/D");
  
      outTree_->Branch("boostWhad_X", &boostWhad_X, "boostWhad_X/D");
      outTree_->Branch("boostWhad_Y", &boostWhad_Y, "boostWhad_Y/D");
      outTree_->Branch("boostWhad_Z", &boostWhad_Z, "boostWhad_Z/D");

      outTree_->Branch("xdotx", &xdotx, "xdotx/D");
      outTree_->Branch("xdoty", &xdoty, "xdoty/D");
      outTree_->Branch("xdotz", &xdotz, "xdotz/D");

      outTree_->Branch("ydotx", &ydotx, "ydotx/D");
      outTree_->Branch("ydoty", &ydoty, "ydoty/D");
      outTree_->Branch("ydotz", &ydotz, "ydotz/D");

      outTree_->Branch("zdotx", &zdotx, "zdotx/D");
      outTree_->Branch("zdoty", &zdoty, "zdoty/D");
      outTree_->Branch("zdotz", &zdotz, "zdotz/D");

      outTree_->Branch("lepton1WWframe_UX", &lepton1WWframe_UX, "lepton1WWframe_UX/D");
      outTree_->Branch("lepton1WWframe_UY", &lepton1WWframe_UY, "lepton1WWframe_UY/D");
      outTree_->Branch("lepton1WWframe_UZ", &lepton1WWframe_UZ, "lepton1WWframe_UZ/D");
   
      outTree_->Branch("lepton_sumWWframe_UX", &lepton_sumWWframe_UX, "lepton_sumWWframe_UX/D");
      outTree_->Branch("lepton_sumWWframe_UY", &lepton_sumWWframe_UY, "lepton_sumWWframe_UY/D");
      outTree_->Branch("lepton_sumWWframe_UZ", &lepton_sumWWframe_UZ, "lepton_sumWWframe_UZ/D");

      /////////////////////////////////////////////////////////////////////////////////////////////
  
      outTree_->Branch("leptons_in_lep_px_good", &leptons_in_lep_px_good, "leptons_in_lep_px_good/I");
      outTree_->Branch("leptons_in_lep_py_good", &leptons_in_lep_py_good, "leptons_in_lep_py_good/I");
      outTree_->Branch("leptons_in_lep_pz_good", &leptons_in_lep_pz_good, "leptons_in_lep_pz_good/I");
  
      outTree_->Branch("partons_in_lep_px_good", &partons_in_lep_px_good, "partons_in_lep_px_good/I");
      outTree_->Branch("partons_in_lep_py_good", &partons_in_lep_py_good, "partons_in_lep_py_good/I");
      outTree_->Branch("partons_in_lep_pz_good", &partons_in_lep_pz_good, "partons_in_lep_pz_good/I");
  
      outTree_->Branch("parton1_in_lep_px_good", &parton1_in_lep_px_good, "parton1_in_lep_px_good/I");
      outTree_->Branch("parton1_in_lep_py_good", &parton1_in_lep_py_good, "parton1_in_lep_py_good/I"); 
      outTree_->Branch("parton1_in_lep_pz_good", &parton1_in_lep_pz_good, "parton1_in_lep_pz_good/I");

      outTree_->Branch("parton2_in_lep_px_good", &parton2_in_lep_px_good, "parton2_in_lep_px_good/I");
      outTree_->Branch("parton2_in_lep_py_good", &parton2_in_lep_py_good, "parton2_in_lep_py_good/I");
      outTree_->Branch("parton2_in_lep_pz_good", &parton2_in_lep_pz_good, "parton2_in_lep_pz_good/I");
  
      outTree_->Branch("lepton1_in_lep_px_good", &lepton1_in_lep_px_good, "lepton1_in_lep_px_good/I");
      outTree_->Branch("lepton1_in_lep_py_good", &lepton1_in_lep_py_good, "lepton1_in_lep_py_good/I");
      outTree_->Branch("lepton1_in_lep_pz_good", &lepton1_in_lep_pz_good, "lepton1_in_lep_pz_good/I");
  
      //outTree_->Branch("lepton1_dotted_x_good", &lepton1_dotted_x_good, "lepton1_dotted_x_good/I");
      //outTree_->Branch("lepton1_dotted_y_good", &lepton1_dotted_y_good, "lepton1_dotted_y_good/I");
      //outTree_->Branch("lepton1_dotted_z_good", &lepton1_dotted_z_good, "lepton1_dotted_z_good/I");
  
      outTree_->Branch("leptons_in_had_px_good", &leptons_in_had_px_good, "leptons_in_had_px_good/I");
      outTree_->Branch("leptons_in_had_py_good", &leptons_in_had_py_good, "leptons_in_had_py_good/I");
      outTree_->Branch("leptons_in_had_pz_good", &leptons_in_had_pz_good, "leptons_in_had_pz_good/I");
  
      outTree_->Branch("lepton1_in_had_px_good", &lepton1_in_had_px_good, "lepton1_in_had_px_good/I");
      outTree_->Branch("lepton1_in_had_py_good", &lepton1_in_had_py_good, "lepton1_in_had_py_good/I");
      outTree_->Branch("lepton1_in_had_pz_good", &lepton1_in_had_pz_good, "lepton1_in_had_pz_good/I");
  
      outTree_->Branch("lepton2_in_had_px_good", &lepton2_in_had_px_good, "lepton2_in_had_px_good/I");
      outTree_->Branch("lepton2_in_had_py_good", &lepton2_in_had_py_good, "lepton2_in_had_py_good/I");
      outTree_->Branch("lepton2_in_had_pz_good", &lepton2_in_had_pz_good, "lepton2_in_had_pz_good/I");
  
      outTree_->Branch("parton1_in_had_px_good", &parton1_in_had_px_good, "parton1_in_had_px_good/I");
      outTree_->Branch("parton1_in_had_py_good", &parton1_in_had_py_good, "parton1_in_had_py_good/I");
      outTree_->Branch("parton1_in_had_pz_good", &parton1_in_had_pz_good, "parton1_in_had_pz_good/I");
  
      //outTree_->Branch("parton1_dotted_x_good", &parton1_dotted_x_good, "parton1_dotted_x_good/I");
      //outTree_->Branch("parton1_dotted_y_good", &parton1_dotted_y_good, "parton1_dotted_y_good/I");
      //outTree_->Branch("parton1_dotted_z_good", &parton1_dotted_z_good, "parton1_dotted_z_good/I");
  
      //outTree_->Branch("complicated1_px_good", &complicated1_px_good, "complicated1_px_good/I");
      //outTree_->Branch("complicated1_py_good", &complicated1_py_good, "complicated1_py_good/I");
      //outTree_->Branch("complicated1_pz_good", &complicated1_pz_good, "complicated1_pz_good/I");
  
      //outTree_->Branch("complicated2_px_good", &complicated2_px_good, "complicated2_px_good/I");
      //outTree_->Branch("complicated2_py_good", &complicated2_py_good, "complicated2_py_good/I");
      //outTree_->Branch("complicated2_pz_good", &complicated2_pz_good, "complicated2_pz_good/I");
  
      outTree_->Branch("lepton_sumWWframe_X_good", &lepton_sumWWframe_X_good, "lepton_sumWWframe_X_good/I");
      outTree_->Branch("lepton_sumWWframe_Y_good", &lepton_sumWWframe_Y_good, "lepton_sumWWframe_Y_good/I");
      outTree_->Branch("lepton_sumWWframe_Z_good", &lepton_sumWWframe_Z_good, "lepton_sumWWframe_Z_good/I");
  
      outTree_->Branch("parton_sumWWframe_X_good", &parton_sumWWframe_X_good, "parton_sumWWframe_X_good/I");
      outTree_->Branch("parton_sumWWframe_Y_good", &parton_sumWWframe_Y_good, "parton_sumWWframe_Y_good/I");
      outTree_->Branch("parton_sumWWframe_Z_good", &parton_sumWWframe_Z_good, "parton_sumWWframe_Z_good/I");
  
      outTree_->Branch("lepton1WWframe_X_good", &lepton1WWframe_X_good, "lepton1WWframe_X_good/I");
      outTree_->Branch("lepton1WWframe_Y_good", &lepton1WWframe_Y_good, "lepton1WWframe_Y_good/I");
      outTree_->Branch("lepton1WWframe_Z_good", &lepton1WWframe_Z_good, "lepton1WWframe_Z_good/I");
  
      outTree_->Branch("parton1WWframe_X_good", &parton1WWframe_X_good, "parton1WWframe_X_good/I");
      outTree_->Branch("parton1WWframe_Y_good", &parton1WWframe_Y_good, "parton1WWframe_Y_good/I");
      outTree_->Branch("parton1WWframe_Z_good", &parton1WWframe_Z_good, "parton1WWframe_Z_good/I");

      //outTree_->Branch("boostWWframe_X_good", &boostWWframe_X_good, "boostWWframe_X_good/I");
      //outTree_->Branch("boostWWframe_Y_good", &boostWWframe_Y_good, "boostWWframe_Y_good/I");
      //outTree_->Branch("boostWWframe_Z_good", &boostWWframe_Z_good, "boostWWframe_Z_good/I");
  
      //outTree_->Branch("boostWlep_X_good", &boostWlep_X_good, "boostWlep_X_good/I");
      //outTree_->Branch("boostWlep_Y_good", &boostWlep_Y_good, "boostWlep_Y_good/I");
      //outTree_->Branch("boostWlep_Z_good", &boostWlep_Z_good, "boostWlep_Z_good/I");
  
      //outTree_->Branch("boostWhad_X_good", &boostWhad_X_good, "boostWhad_X_good/I");
      //outTree_->Branch("boostWhad_Y_good", &boostWhad_Y_good, "boostWhad_Y_good/I");
      //outTree_->Branch("boostWhad_Z_good", &boostWhad_Z_good, "boostWhad_Z_good/I");
#endif 
    }

  if (isSignal)
    {
      outTree_->Branch("aTGCWeights",  &aTGCWeights);
      outTree_->Branch("refXsec", &refXsec, "refXsec/D");
    }
  aTGCWeightUnitConv=1;

}


TreeMaker::~TreeMaker()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  nevent = iEvent.eventAuxiliary().event();
  run    = iEvent.eventAuxiliary().run();
  lumi   = iEvent.eventAuxiliary().luminosityBlock();

  edm::Handle< double > rhoH;
  iEvent.getByToken(rhoToken_,rhoH);
  rho_ = *rhoH;
      
  //Leptonic Ws
  edm::Handle<edm::View<reco::Candidate> > leptonicVs;
  iEvent.getByToken(leptonicVToken_, leptonicVs);
   
  //GenParticles
  edm::Handle<edm::View<reco::Candidate> > genParticles;
  if(isMC)iEvent.getByToken(genParticlesToken_, genParticles); 
   
  //Vertices 
  edm::Handle<edm::View<reco::Vertex> > vertices;
  iEvent.getByToken(vertexToken_, vertices);
   
  //Jets
  edm::Handle<edm::View<pat::Jet> > jets; 
  iEvent.getByToken(fatJetsToken_, jets);

  edm::Handle<edm::View<pat::Jet> > jetsSmearedUp;
  if(isMC) iEvent.getByToken(fatJetsSmearedUpToken_, jetsSmearedUp);

  edm::Handle<edm::View<pat::Jet> > jetsSmearedDown;
  if(isMC) iEvent.getByToken(fatJetsSmearedDownToken_, jetsSmearedDown);

  edm::Handle<edm::View<reco::GenJet>> genJetsAK8;
  if(isMC) iEvent.getByToken(genJetsAK8Token_, genJetsAK8);

  //AK4 Jets (for Btag veto )
  edm::Handle<edm::View<pat::Jet> > AK4Jets;
  iEvent.getByToken(AK4JetsToken_, AK4Jets);
  njets = AK4Jets->size();

  njets_JERUp = 0;
  njets_JERDown = 0;
  njets_JECUp = 0;
  njets_JECDown = 0;
  edm::Handle<edm::View<pat::Jet> > AK4JetsSmearedUp, AK4JetsSmearedDown, AK4JetsShiftedUp, AK4JetsShiftedDown;
  if(isMC)
    {
      iEvent.getByToken(AK4JetsSmearedUpToken_, AK4JetsSmearedUp);
      njets_JERUp = AK4JetsSmearedUp->size();
      iEvent.getByToken(AK4JetsSmearedDownToken_, AK4JetsSmearedDown);
      njets_JERDown = AK4JetsSmearedDown->size();

      iEvent.getByToken(AK4JetsShiftedUpToken_, AK4JetsShiftedUp);
      njets_JECUp = AK4JetsShiftedUp->size();
      iEvent.getByToken(AK4JetsShiftedDownToken_, AK4JetsShiftedDown);
      njets_JECDown = AK4JetsShiftedDown->size();
    }
   
  //MET
  edm::Handle<edm::View<pat::MET> > metHandle;
  iEvent.getByToken(metToken_, metHandle);
   
  //loose electrons
  edm::Handle<edm::View<reco::Candidate> > looseElectrons;
  iEvent.getByToken(looseEleToken_, looseElectrons);
   
  //loose muons
  edm::Handle<edm::View<reco::Candidate> > looseMuons;
  iEvent.getByToken(looseMuToken_, looseMuons); 
   
  //leptons (tight)
  edm::Handle<edm::View<reco::Candidate> > leptons;
  iEvent.getByToken(leptonsToken_, leptons); 

  std::map<std::string, math::XYZTLorentzVector>  SystMap; 
  std::map<std::string, math::XYZTLorentzVector>  LeptonSystMap;
  std::map<std::string, math::XYZTLorentzVector>  MetSystMap;
  if (isMC)
    {
      SystMap = SystematicsHelper_.getWLepSystematicsLoretzVectors(iEvent);
      LeptonSystMap = SystematicsHelper_.getLeptonSystematicsLoretzVectors(iEvent);
      MetSystMap = SystematicsHelper_.getMetSystematicsLoretzVectors(iEvent);
    }

  nPV = vertices->size();
   
  edm::Handle<std::vector< PileupSummaryInfo > >  PUInfo;
  edm::Handle <GenEventInfoProduct> genInfo; 
   
  if (isMC)
    {
      iEvent.getByToken(PUInfoToken_, PUInfo);
      std::vector<PileupSummaryInfo>::const_iterator PVI;
      float Tnpv = -1;
      for(PVI = PUInfo->begin(); PVI != PUInfo->end(); ++PVI)
	{
	  int BX = PVI->getBunchCrossing();
	  if(BX == 0)
	    { 
	      Tnpv = PVI->getTrueNumInteractions();
	      continue;
	    }
	}
      gnPV=Tnpv;
      PUweight = LumiWeights_.weight(Tnpv);

      iEvent.getByToken(genInfoToken , genInfo);
      genWeight = (genInfo -> weight());
      //btag weights
      if (njets > 0)
	{
	  btagWeight =  BTagHelper_.getEventWeight(AK4Jets);
	  if (btagWeight != btagWeight)
	    {
	      std::cout << "Event: " << nevent << " LS: " << lumi << std::endl;
	      throw std::runtime_error("Inf/NaN btagWeight for event ");
	    }
	  btagWeight_BTagUp =  BTagHelper_.getEventWeight(AK4Jets, UP, BTAG);
	  btagWeight_BTagDown =  BTagHelper_.getEventWeight(AK4Jets, DOWN, BTAG);
	  btagWeight_MistagUp =  BTagHelper_.getEventWeight(AK4Jets, UP, MISTAG);
	  btagWeight_MistagDown =  BTagHelper_.getEventWeight(AK4Jets, DOWN, MISTAG);
	}
      else
	{
	  btagWeight = 1.;
	  btagWeight_BTagUp = 1.;
	  btagWeight_BTagDown = 1.;
	  btagWeight_MistagUp = 1.;
	  btagWeight_MistagDown = 1.;
	}

      //PDF uncertainties
      edm::Handle<LHEEventProduct> LHEevtProductExternal;
      iEvent.getByToken(LHEEventProductTokenExternal, LHEevtProductExternal);
      range PDFRange = PDFVariationMap.at(NominalPDF);
      if(isSignal) PDFRange = PDFVariationMap.at(NominalPDF/100);
   
      //define number of PDF variations 
      unsigned int NPDFs = PDFRange.high - PDFRange.low + 1;
      PDFWeights.clear();
      ScaleWeights.clear();
      PDFWeights.resize(NPDFs);
      ScaleWeights.resize(9);

      unsigned int iPDF_ID = PDFRange.low;
      //if there are no weights for PDF uncertainties just fill with ones, that's the case for tW single top sample
      if (LHEevtProductExternal->weights().size() == 0) std::fill(PDFWeights.begin(), PDFWeights.end(), 1.);
      for (unsigned int i = 0; i < LHEevtProductExternal->weights().size(); i++)
	{
	  if (iPDF_ID > PDFRange.high) break;
	  if (LHEevtProductExternal->weights()[i].id == std::to_string(iPDF_ID))
	    {
	      unsigned int iPDF = iPDF_ID - PDFRange.low;
	      PDFWeights.at(iPDF) = (LHEevtProductExternal->weights()[i].wgt)/LHEevtProductExternal->originalXWGTUP();
	      iPDF_ID++;
	    }
	}	

      //scale variation uncertainties
      range RangeOfScaleVariation;
      if (NominalPDF == 263000) RangeOfScaleVariation = range(1,9);
      else RangeOfScaleVariation = range(1001, 1009);

      //if there are no weights for scale uncertainties just fill with ones, that's the case for tW single top  sample
      if (LHEevtProductExternal->weights().size() == 0 ) std::fill(ScaleWeights.begin(), ScaleWeights.end(), 1.);

      unsigned int iScale_ID = RangeOfScaleVariation.low;
      for (unsigned int i=0; i<LHEevtProductExternal->weights().size(); i++)
	{
	  if (iScale_ID > RangeOfScaleVariation.high) break;
	  if (LHEevtProductExternal->weights()[i].id == std::to_string(iScale_ID))
	    {
	      unsigned int iScale = iScale_ID - RangeOfScaleVariation.low;
	      ScaleWeights.at(iScale) = (LHEevtProductExternal->weights()[i].wgt)/LHEevtProductExternal->originalXWGTUP();
	      iScale_ID++;
	    }
	}


      if(isSignal)
	{
	  aTGCWeights.clear();
	  refXsec = LHEevtProductExternal -> originalXWGTUP();
	  int weightNumber = 1;
	  if(LHEevtProductExternal->weights().size())
	    {
	      aTGCWeightUnitConv = genWeight/(LHEevtProductExternal->weights().at(0).wgt);
	      for (unsigned int iwgt = 0; iwgt < LHEevtProductExternal->weights().size(); ++iwgt)
		{
		  const LHEEventProduct::WGT& wgt = LHEevtProductExternal->weights().at(iwgt);
		  if(boost::algorithm::contains(wgt.id, "rwgt_" + std::to_string(weightNumber)))
		    {
		      aTGCWeights.push_back(wgt.wgt);
		      weightNumber++;
		    }
		}
	    }
	}

      ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > Wplus_p4 = genWLorentzVector(genParticles, 1);
      ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > Wminus_p4 = genWLorentzVector(genParticles, -1);
      Wplus_gen_pt = Wplus_p4.Pt();
      Wplus_gen_eta = Wplus_p4.Eta();
      Wplus_gen_phi = Wplus_p4.Phi();
      Wplus_gen_mass = Wplus_p4.M();

      Wminus_gen_pt = Wminus_p4.Pt();
      Wminus_gen_eta = Wminus_p4.Eta();
      Wminus_gen_phi = Wminus_p4.Phi();
      Wminus_gen_mass = Wminus_p4.M();
    }

  // Store generator level vector boson pT and diboson mass
  genPtV = 0;
  genMWV = 0;
  if(isMC)
    {
      // Scan the list of gen particles, isolate the starting W or Z, identify the hadronically decaying one and store its pt.
      for(unsigned int iGen = 0; iGen<genParticles->size(); ++iGen)
	{
	  const reco::Candidate & p = genParticles->at(iGen);
	  /*std::cout<<p.pdgId()<<" "<<p.status();
	    genHadronicVp4.SetPt(p.pt());
	    genHadronicVp4.SetEta(p.eta());
	    genHadronicVp4.SetPhi(p.phi());
	    genHadronicVp4.SetM(p.mass());
	    std::cout<<" "<<genHadronicVp4.px()<<" "<<genHadronicVp4.py()<<" "<<genHadronicVp4.pz()<<" "<<genHadronicVp4.e()<<" "<<genHadronicVp4.M()<<std::endl;*/
	  if((fabs(p.pdgId()) == 23 || fabs(p.pdgId()) == 24) && p.status() == 22)
	    {
	      if(decaysHadronic(&p))
		{
		  genPtV=p.pt();

		  genHadronicVp4.SetPt(p.pt());
		  genHadronicVp4.SetEta(p.eta());
		  genHadronicVp4.SetPhi(p.phi());
		  genHadronicVp4.SetM(p.mass());
		  //std::cout<<"*Hadronically decaying W* "<<genHadronicVp4.px()<<" "<<genHadronicVp4.py()<<" "<<genHadronicVp4.pz()<<" "<<genHadronicVp4.e()<<" "<<genHadronicVp4.M()<<std::endl;
		}
	      else
		{
		  genLeptonicVp4.SetPt(p.pt());
		  genLeptonicVp4.SetEta(p.eta());
		  genLeptonicVp4.SetPhi(p.phi());
		  genLeptonicVp4.SetM(p.mass());
		  //std::cout<<"*Leptonically decaying W* "<<genLeptonicVp4.px()<<" "<<genLeptonicVp4.py()<<" "<<genLeptonicVp4.pz()<<" "<<genLeptonicVp4.e()<<" "<<genLeptonicVp4.M()<<std::endl;
		}
	    }
	}
      genDibosonp4 = genHadronicVp4+genLeptonicVp4;
      genMWV = genDibosonp4.M();
      //std::cout<<"Diboson mass: "<<genMWV<<std::endl<<std::endl;
    }

  if (isMC && jecUnc == nullptr)
    {
      // Only have to initialise once
      edm::ESHandle<JetCorrectorParametersCollection> JetCorParCollAK8;
      iSetup.get<JetCorrectionsRecord>().get("AK8PFchs",JetCorParCollAK8);
      JetCorrectorParameters const & JetCorPar = (*JetCorParCollAK8)["Uncertainty"];
      jecUnc = new JetCorrectionUncertainty(JetCorPar);
    }
  JECunc = 0.;

  //  Defining decay channel on the gen level
  N_had_Wgen  = 0, N_lep_Wgen = 0 ;

  if (isMC)
    {
      DefineDecayChannel(genParticles, N_lep_Wgen , N_had_Wgen);
      if (N_lep_Wgen == 0 && N_had_Wgen == 2) WDecayClass = Hadronic;
      else if (N_lep_Wgen == 1 && N_had_Wgen == 1) WDecayClass = Semileptonic;
      else if (N_lep_Wgen == 2 && N_had_Wgen == 0) WDecayClass = Leptonic;
      else WDecayClass = UnDefined;
    }

  N_lep_W = leptonicVs->size();
   
  //loose leptons
  nLooseEle = looseElectrons->size();
  nLooseMu = looseMuons->size();
   
  nLep = leptons->size();
   
  
  //electron channel
  if (( leptons->size()) > 0)
    {
      auto leptonPtr = leptons -> ptrAt(0);
      reco::MuonPtr asmuonPtr(leptonPtr);
      reco::ElectronPtr aselectronPtr(leptonPtr);
      if (channel == "mu")
	{
	  Lepton.pt =  asmuonPtr->tunePMuonBestTrack()->pt();
	  Lepton.eta = asmuonPtr->tunePMuonBestTrack()->eta();
	  Lepton.phi = asmuonPtr->tunePMuonBestTrack()->phi();
	}
      else if (channel == "el")
	{
	  Lepton.pt = (leptons->at(0)).pt();
	  Lepton.eta = (leptons->at(0)).eta();
	  Lepton.phi = (leptons->at(0)).phi();
	  reco::SuperClusterRef  superCluster = aselectronPtr->superCluster();
	  sc_eta = superCluster->eta();
	  sc_et = (superCluster->energy())* sin((aselectronPtr->superClusterPosition()).theta());
	  const reco::CaloClusterPtr& seed = aselectronPtr -> superCluster()->seed();
	  isEB = (seed->seed().subdetId() == EcalBarrel);
	  triggerWeightHLTEle27NoER  =  trigEle27NoER::turnOn(sc_et, sc_eta);	
	}
      if (isMC)
	{
	  Lepton.pt_LeptonEnUp = LeptonSystMap.at("LeptonEnUp").Pt();
	  Lepton.pt_LeptonEnDown = LeptonSystMap.at("LeptonEnDown").Pt();
	  Lepton.pt_LeptonResUp = LeptonSystMap.at("LeptonResUp").Pt();
	  Lepton.pt_LeptonResDown = LeptonSystMap.at("LeptonResDown").Pt();
	}
    }
  else
    {
      Lepton.pt = -99.;
      Lepton.eta = -99.;
      Lepton.phi = -99.;
      Lepton.pt_LeptonEnUp = -99.;
      Lepton.pt_LeptonEnDown = -99.;
      Lepton.pt_LeptonResUp = -99.;
      Lepton.pt_LeptonResDown = -99.;
    }

  if (channel == "mu")
    {
      // Muon scale factors apply cuts on their nVtx
      int goodNPV = 0;
      for (const auto & itr : *vertices) {if (fabs(itr.z()) <= 25 && itr.ndof() > 4 && fabs(itr.position().rho()) <= 2 && !itr.isFake()) {goodNPV++;}}

      LeptonSF = MuonScaleFactor_.getScaleFactor(Lepton.pt, std::abs(Lepton.eta), Lepton.phi, goodNPV);
      LeptonSF_Up = MuonScaleFactor_.getScaleFactor(Lepton.pt, std::abs(Lepton.eta), Lepton.phi, goodNPV, "up");
      LeptonSF_Down = MuonScaleFactor_.getScaleFactor(Lepton.pt, std::abs(Lepton.eta), Lepton.phi, goodNPV, "down");
    }
  else if (channel == "el")
    {
      LeptonSF = ElectronScaleFactor_.getScaleFactor(Lepton.pt, Lepton.eta, sc_eta, Lepton.phi, nPV);
      LeptonSF_Up = ElectronScaleFactor_.getScaleFactor(Lepton.pt, Lepton.eta, sc_eta, Lepton.phi, nPV, "up");
      LeptonSF_Down = ElectronScaleFactor_.getScaleFactor(Lepton.pt, Lepton.eta, sc_eta, Lepton.phi, nPV, "down");
    }
  //leptonically decaying W
  if (leptonicVs->size() > 0)
    {
      const reco::Candidate  & leptonicV = leptonicVs->at(0);   
      Wboson_lep.pt = leptonicV.pt();
      Wboson_lep.eta = leptonicV.eta();
      Wboson_lep.phi = leptonicV.phi();
      Wboson_lep.mass = leptonicV.mass();
      Wboson_lep.mt = leptonicV.mt();
      Wboson_lep.charge = leptonicV.charge();
      if (isMC)
	{
	  //pt
	  Wboson_lep.pt_LeptonEnUp = SystMap.at("LeptonEnUp").Pt();
	  Wboson_lep.pt_LeptonEnDown = SystMap.at("LeptonEnDown").Pt();
	  Wboson_lep.pt_LeptonResUp = SystMap.at("LeptonResUp").Pt();
	  Wboson_lep.pt_LeptonResDown = SystMap.at("LeptonResDown").Pt();
	  Wboson_lep.pt_JECUp = SystMap.at("JetEnUp").Pt();
	  Wboson_lep.pt_JECDown = SystMap.at("JetEnDown").Pt();
	  Wboson_lep.pt_UnclEnUp = SystMap.at("UnclusteredEnUp").Pt();
	  Wboson_lep.pt_UnclEnDown = SystMap.at("UnclusteredEnDown").Pt();

	  //mass
	  Wboson_lep.mass_LeptonEnUp = SystMap.at("LeptonEnUp").M();
	  Wboson_lep.mass_LeptonEnDown = SystMap.at("LeptonEnDown").M();
	  Wboson_lep.mass_LeptonResUp = SystMap.at("LeptonResUp").M();
	  Wboson_lep.mass_LeptonResDown = SystMap.at("LeptonResDown").M();
	  Wboson_lep.mass_JECUp = SystMap.at("JetEnUp").M();
	  Wboson_lep.mass_JECDown = SystMap.at("JetEnDown").M();
	  Wboson_lep.mass_UnclEnUp = SystMap.at("UnclusteredEnUp").M();
	  Wboson_lep.mass_UnclEnDown = SystMap.at("UnclusteredEnDown").M();


	  //mt
	  Wboson_lep.mt_LeptonEnUp = SystMap.at("LeptonEnUp").Mt();
	  Wboson_lep.mt_LeptonEnDown = SystMap.at("LeptonEnDown").Mt();
	  Wboson_lep.mt_LeptonResUp = SystMap.at("LeptonResUp").Mt();
	  Wboson_lep.mt_LeptonResDown = SystMap.at("LeptonResDown").Mt();
	  Wboson_lep.mt_JECUp = SystMap.at("JetEnUp").Mt();
	  Wboson_lep.mt_JECDown = SystMap.at("JetEnDown").Mt();
	  Wboson_lep.mt_UnclEnUp = SystMap.at("UnclusteredEnUp").Mt();
	  Wboson_lep.mt_UnclEnDown = SystMap.at("UnclusteredEnDown").Mt();

	}
    }
  
  else 
    {
      Wboson_lep.pt = -99.;
      Wboson_lep.eta = -99.;
      Wboson_lep.phi = -99.;
      Wboson_lep.mass = -99.;
      Wboson_lep.mt = -99.;
      Wboson_lep.charge = -99.;
      if (isMC)
	{
	  //pt
	  Wboson_lep.pt_LeptonEnUp = -99.;
	  Wboson_lep.pt_LeptonEnDown = -99.;
	  Wboson_lep.pt_LeptonResUp = -99.;
	  Wboson_lep.pt_LeptonResDown = -99.;
	  Wboson_lep.pt_JECUp = -99.;
	  Wboson_lep.pt_JECDown = -99.;
	  Wboson_lep.pt_UnclEnUp = -99.;
	  Wboson_lep.pt_UnclEnDown = -99.;

	  //mass
	  Wboson_lep.mass_LeptonEnUp = -99.;
	  Wboson_lep.mass_LeptonEnDown = -99.;
	  Wboson_lep.mass_LeptonResUp = -99.;
	  Wboson_lep.mass_LeptonResDown = -99.;
	  Wboson_lep.mass_JECUp = -99.;
	  Wboson_lep.mass_JECDown = -99.;
	  Wboson_lep.mass_UnclEnUp = -99.;
	  Wboson_lep.mass_UnclEnDown = -99.;

	  //mt 

	  Wboson_lep.mt_LeptonEnUp = -99.;
	  Wboson_lep.mt_LeptonEnDown = -99.;
	  Wboson_lep.mt_LeptonResUp = -99.;
	  Wboson_lep.mt_LeptonResDown = -99.;
	  Wboson_lep.mt_JECUp = -99.;
	  Wboson_lep.mt_JECDown = -99.;
	  Wboson_lep.mt_UnclEnUp = -99.;
	  Wboson_lep.mt_UnclEnDown = -99.;


	}
    }
  
  //MET quantities
  if (metHandle->size() > 0)
    {
      const pat::MET& metCand = metHandle->at(0);

      METCand.pt = metCand.pt();
      METCand.phi = metCand.phi();

      //MET uncertainties
      if (isMC)
	{
	  // MET uncertainties
	  MET_LeptonEnUp = MetSystMap.at("LeptonEnUp").Pt();
	  MET_LeptonEnDown = MetSystMap.at("LeptonEnDown").Pt();
	  MET_JECUp = MetSystMap.at("JetEnUp").Pt();
	  MET_JECDown = MetSystMap.at("JetEnDown").Pt();
	  MET_JERUp = MetSystMap.at("JetResUp").Pt();
	  MET_JERDown = MetSystMap.at("JetResDown").Pt();
	  MET_UnclEnUp = MetSystMap.at("UnclusteredEnUp").Pt();
	  MET_UnclEnDown = MetSystMap.at("UnclusteredEnDown").Pt();

	  //MET phi uncertainties
	  MET_Phi_LeptonEnUp = MetSystMap.at("LeptonEnUp").Phi();
	  MET_Phi_LeptonEnDown = MetSystMap.at("LeptonEnDown").Phi();
	  MET_Phi_JECUp = MetSystMap.at("JetEnUp").Phi();
	  MET_Phi_JECDown = MetSystMap.at("JetEnDown").Phi();
	  MET_Phi_JERUp = MetSystMap.at("JetResUp").Phi();
	  MET_Phi_JERDown = MetSystMap.at("JetResDown").Phi();
	  MET_Phi_UnclEnUp = MetSystMap.at("UnclusteredEnUp").Phi();
	  MET_Phi_UnclEnDown = MetSystMap.at("UnclusteredEnDown").Phi();
	}
    }
  else
    {
      METCand.pt = -99.;
      METCand.phi = -99.;

      //MET uncertainties
      //METUncl
      MET_UnclEnUp = -99. ;
      MET_UnclEnUp = -99.;
      //JER
      MET_JERUp = -99.;
      MET_JERDown = -99. ;
      //JEC
      MET_JECUp = -99. ;
      MET_JECDown = -99. ;
      //Lepton energy 
      MET_LeptonEnUp = -99.;
      MET_LeptonEnDown = -99. ;

      //MET phi uncertainties
      //METUncl
      MET_Phi_UnclEnUp = -99. ;
      MET_Phi_UnclEnUp = -99.;
      //JER
      MET_Phi_JERUp = -99.;
      MET_Phi_JERDown = -99. ;
      //JEC
      MET_Phi_JECUp = -99. ;
      MET_Phi_JECDown = -99. ;
      //Lepton energy 
      MET_Phi_LeptonEnUp = -99.;
      MET_Phi_LeptonEnDown = -99. ;
      
    }
   
  if (jets -> size() > 0 && leptonicVs -> size() > 0)
    {
      deltaR_LeptonWJet = deltaR(Lepton.eta,Lepton.phi,(jets -> at(0)).eta(), (jets -> at(0)).phi()); 
      deltaPhi_LeptonMet = deltaPhi(Lepton.phi, METCand.phi);
      deltaPhi_WJetMet = deltaPhi(jets->at(0).phi(), METCand.phi);
      deltaPhi_WJetWlep = deltaPhi(jets->at(0).phi(), Wboson_lep.phi);
      if (isMC)
	{
	  //Unclustered energy
	  deltaPhi_LeptonMet_UnclEnUp = deltaPhi(Lepton.phi, MET_Phi_UnclEnUp);
	  deltaPhi_LeptonMet_UnclEnDown = deltaPhi(Lepton.phi, MET_Phi_UnclEnDown);
	  //JEC
	  deltaPhi_LeptonMet_JECUp = deltaPhi(Lepton.phi, MET_Phi_JECUp);
	  deltaPhi_LeptonMet_JECDown = deltaPhi(Lepton.phi, MET_Phi_JECDown);
	  //lepton energy scale
	  deltaPhi_LeptonMet_LeptonEnUp = deltaPhi(Lepton.phi, MET_Phi_LeptonEnUp);
	  deltaPhi_LeptonMet_LeptonEnDown = deltaPhi(Lepton.phi, MET_Phi_LeptonEnDown);
	  //JER
	  deltaPhi_LeptonMet_JERUp = deltaPhi(Lepton.phi, MET_Phi_JERUp);
	  deltaPhi_LeptonMet_JERDown = deltaPhi(Lepton.phi, MET_Phi_JERDown);

	  //////////////////////
	  //Unclustered energy
	  deltaPhi_WJetMet_UnclEnUp = deltaPhi(jets->at(0).phi(), MET_Phi_UnclEnUp);
	  deltaPhi_WJetMet_UnclEnDown = deltaPhi(jets->at(0).phi(), MET_Phi_UnclEnDown);
	  //JEC
	  deltaPhi_WJetMet_JECUp = deltaPhi(jets->at(0).phi(), MET_Phi_JECUp);
	  deltaPhi_WJetMet_JECDown = deltaPhi(jets->at(0).phi(), MET_Phi_JECDown);
	  //lepton energy scale
	  deltaPhi_WJetMet_LeptonEnUp = deltaPhi(jets->at(0).phi(), MET_Phi_LeptonEnUp);
	  deltaPhi_WJetMet_LeptonEnDown = deltaPhi(jets->at(0).phi(), MET_Phi_LeptonEnDown);
	  //JER
	  deltaPhi_WJetMet_JERUp = deltaPhi(jets->at(0).phi(), MET_Phi_JERUp);
	  deltaPhi_WJetMet_JERDown = deltaPhi(jets->at(0).phi(), MET_Phi_JERDown);

	  //////////////////////
	  //Unclustered energy
	  deltaPhi_WJetWlep_UnclEnUp = deltaPhi(jets->at(0).phi(), SystMap.at("UnclusteredEnUp").Phi());
	  deltaPhi_WJetWlep_UnclEnDown = deltaPhi(jets->at(0).phi(), SystMap.at("UnclusteredEnDown").Phi());
	  //JEC
	  deltaPhi_WJetWlep_JECUp = deltaPhi(jets->at(0).phi(), SystMap.at("JetEnUp").Phi());
	  deltaPhi_WJetWlep_JECDown = deltaPhi(jets->at(0).phi(), SystMap.at("JetEnDown").Phi());
	  //lepton energy scale
	  deltaPhi_WJetWlep_LeptonEnUp = deltaPhi(jets->at(0).phi(), SystMap.at("LeptonEnUp").Phi());
	  deltaPhi_WJetWlep_LeptonEnDown = deltaPhi(jets->at(0).phi(), SystMap.at("LeptonEnDown").Phi());


	}
    }
  else 
    {
      deltaR_LeptonWJet = -99.; 
      deltaPhi_LeptonMet = -99.;
      deltaPhi_WJetMet = -99.;
      deltaPhi_WJetWlep = -99.;

      if (isMC)
	{
	  //Unclustered energy
	  deltaPhi_LeptonMet_UnclEnUp = -99.;
	  deltaPhi_LeptonMet_UnclEnDown = -99.;
	  //JEC
	  deltaPhi_LeptonMet_JECUp = -99.;
	  deltaPhi_LeptonMet_JECDown = -99.;
	  //lepton energy scale
	  deltaPhi_LeptonMet_LeptonEnUp = -99.;
	  deltaPhi_LeptonMet_LeptonEnDown = -99.;
	  //JER
	  deltaPhi_LeptonMet_JERUp = -99.;
	  deltaPhi_LeptonMet_JERDown = -99.;
	  /////////////
	  //Unclustered energy
	  deltaPhi_WJetMet_UnclEnUp = -99.;
	  deltaPhi_WJetMet_UnclEnDown = -99.;
	  //JEC
	  deltaPhi_WJetMet_JECUp = -99.;
	  deltaPhi_WJetMet_JECDown = -99.;
	  //lepton energy scale
	  deltaPhi_WJetMet_LeptonEnUp = -99.;
	  deltaPhi_WJetMet_LeptonEnDown = -99.;
	  //JER
	  deltaPhi_WJetMet_JERUp = -99.;
	  deltaPhi_WJetMet_JERDown = -99.;

	  //////////////////////
	  //Unclustered energy
	  deltaPhi_WJetWlep_UnclEnUp = -99.;
	  deltaPhi_WJetWlep_UnclEnDown = -99.;
	  //JEC
	  deltaPhi_WJetWlep_JECUp = -99.;
	  deltaPhi_WJetWlep_JECDown = -99.;
	  //lepton energy scale
	  deltaPhi_WJetWlep_LeptonEnUp = -99.;
	  deltaPhi_WJetWlep_LeptonEnDown = -99.;
	}

    }


  NAK8jet = jets -> size();
  if(isMC) NAK8jet_smearedUp = jetsSmearedUp -> size();
  if(isMC) NAK8jet_smearedDown = jetsSmearedDown -> size();
  JetResolutionSmearer_.setRhoAndSeed(rho_, iEvent);

  // Different types because the pat::Jet returns math::XYZTLorentzVector (uses E not M)
  // but for the SD ones we need to modify the mass, so need the M4D not E4D
  math::XYZTLorentzVector smearedJetUp, smearedJetDown;
  ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > smearedJetUp_SD, smearedJetDown_SD;

  if (jets -> size() > 0)
    {
      const pat::Jet & fatJet = jets->at(0);
      jet_tau2tau1 = fatJet.userFloat("NjettinessAK8:tau2")/fatJet.userFloat("NjettinessAK8:tau1");
      jet_tau3tau2 = fatJet.userFloat("NjettinessAK8:tau3")/fatJet.userFloat("NjettinessAK8:tau2");
      jet_tau1 = fatJet.userFloat("NjettinessAK8:tau1");
      jet_tau2 = fatJet.userFloat("NjettinessAK8:tau2");
      jet_tau3 = fatJet.userFloat("NjettinessAK8:tau3");
      // Need to manually apply correction factor to userFloat values - calculate by using
      // ratio of corrected PT to uncorrected pt
      double corr = fatJet.correctedP4(fatJet.currentJECLevel()).pt() / fatJet.correctedP4("Uncorrected").pt();
      jet_mass_pruned = corr * fatJet.userFloat("ak8PFJetsCHSPrunedMass");
      jet_mass_softdrop = corr * fatJet.userFloat("ak8PFJetsCHSSoftDropMass");

      //https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2016#Jets
      // Do these require correction factor?
      jet_pt_PUPPI = fatJet.userFloat("ak8PFJetsPuppiValueMap:pt");
      jet_eta_PUPPI = fatJet.userFloat("ak8PFJetsPuppiValueMap:eta");
      jet_phi_PUPPI = fatJet.userFloat("ak8PFJetsPuppiValueMap:phi");
      jet_mass_PUPPI = fatJet.userFloat("ak8PFJetsPuppiValueMap:mass");
      jet_tau1_PUPPI = fatJet.userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau1");
      jet_tau2_PUPPI = fatJet.userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau2");
      jet_tau3_PUPPI = fatJet.userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau3");
      jet_tau21_PUPPI = jet_tau2_PUPPI/jet_tau1_PUPPI;
      jet_tau32_PUPPI = jet_tau3_PUPPI/jet_tau2_PUPPI;


      TLorentzVector puppi_softdrop, puppi_softdrop_subjet;
      TLorentzVector zerothsubjet, firstsubjet;
      number_of_subjets = 0;
      auto const & sdSubjetsPuppi = fatJet.subjets("SoftDropPuppi");
      for ( auto const & it : sdSubjetsPuppi )
	{
	  puppi_softdrop_subjet.SetPtEtaPhiM(it->correctedP4(0).pt(),it->correctedP4(0).eta(),it->correctedP4(0).phi(),it->correctedP4(0).mass());
	  puppi_softdrop+=puppi_softdrop_subjet;
	  if (number_of_subjets == 0)
	    {
	      zerothsubjet.SetPtEtaPhiM(it->correctedP4(0).pt(),it->correctedP4(0).eta(),it->correctedP4(0).phi(),it->correctedP4(0).mass());
	      subjet0.pt = it->correctedP4(0).pt(); subjet0.eta = it->correctedP4(0).eta(); subjet0.phi = it->correctedP4(0).phi(); subjet0.mass = it->correctedP4(0).mass();
	    }
	  else if (number_of_subjets == 1)
	    {
	      firstsubjet.SetPtEtaPhiM(it->correctedP4(0).pt(),it->correctedP4(0).eta(),it->correctedP4(0).phi(),it->correctedP4(0).mass());
	      subjet1.pt = it->correctedP4(0).pt(); subjet1.eta = it->correctedP4(0).eta(); subjet1.phi = it->correctedP4(0).phi(); subjet1.mass = it->correctedP4(0).mass();
	    }
	  else std::cout << "More than two subjets counted.\n";
	  number_of_subjets++;
	}

      zerothsubjet_px = zerothsubjet.Px();
      zerothsubjet_py = zerothsubjet.Py();
      zerothsubjet_pz = zerothsubjet.Pz();
      zerothsubjet_e  = zerothsubjet.E();

      firstsubjet_px = firstsubjet.Px();
      firstsubjet_py = firstsubjet.Py();
      firstsubjet_pz = firstsubjet.Pz();
      firstsubjet_e  = firstsubjet.E();

      /////////////////////////////////////
      
      zerothsubjet_pt = zerothsubjet.Pt();
      zerothsubjet_phi = zerothsubjet.Phi();
      zerothsubjet_eta = zerothsubjet.Eta();
      zerothsubjet_m  = zerothsubjet.M();

      firstsubjet_pt = firstsubjet.Pt();
      firstsubjet_phi = firstsubjet.Phi();
      firstsubjet_eta = firstsubjet.Eta();
      firstsubjet_m  = firstsubjet.M();

      float puppiCorr= getPUPPIweight(jet_pt_PUPPI, jet_eta_PUPPI);
      jet_mass_softdrop_PUPPI = puppi_softdrop.M() * puppiCorr;
      jet_tau21_DT = jet_tau21_PUPPI + 0.063*std::log(jet_pt_PUPPI*jet_pt_PUPPI/jet_mass_PUPPI);

      jet_pt = fatJet.pt();
      jet_eta = fatJet.eta();
      jet_phi = fatJet.phi();
      jet_mass = fatJet.mass();

      if(isMC)
	{
	  isMatched_ = isMatchedToGenW(genParticles, fatJet);

	  //JEC uncertainty
	  jecUnc->setJetEta(jet_eta);
	  jecUnc->setJetPt(jet_pt);  // must be corrected PT
	  JECunc = jecUnc->getUncertainty(true);

	  jet_pt_JECDown = (1 - JECunc)*jet_pt;
	  jet_pt_JECUp   = (1 + JECunc)*jet_pt;
	  jet_mass_JECDown = (1 - JECunc)*jet_mass;
	  jet_mass_JECUp   = (1 + JECunc)*jet_mass;
	  jet_mass_pruned_JECDown = (1 - JECunc)*jet_mass_pruned;
	  jet_mass_pruned_JECUp = (1 + JECunc)*jet_mass_pruned;
	  jet_mass_softdrop_JECDown = (1 - JECunc)*jet_mass_softdrop;
	  jet_mass_softdrop_JECUp = (1 + JECunc)*jet_mass_softdrop;
	  jet_mass_softdrop_PUPPI_JECDown = (1 - JECunc)*jet_mass_softdrop_PUPPI;
	  jet_mass_softdrop_PUPPI_JECUp = (1 + JECunc)*jet_mass_softdrop_PUPPI;

	  // Numbers taken from JetWTagging twiki
	  float mSDSF = 1.0;
	  float mSDSFUnc = 0.2;
	  float mSDResolutionAbs = 10.1;
	  float mSDResolutionRel = mSDResolutionAbs / 80.; // FIXME! Need better number than 80

	  //JER uncertainty
	  if (jetsSmearedUp->size() > 0)
	    {
	      const pat::Jet & fatJetUp = jetsSmearedUp->at(0);
	      jet_pt_JERUp = fatJetUp.pt();
	      jet_mass_JERUp = fatJetUp.mass();
	      double JERUpCorrection = fatJetUp.pt()/jet_pt;
	      jet_mass_pruned_JERUp = JERUpCorrection*jet_mass_pruned;
	      jet_mass_softdrop_JERUp = JERUpCorrection*jet_mass_softdrop;
	      // For PUPPI SD mass, we treat it separately using official JMR SF & unc, and resolution.
	      // We don't have a gen level mass, so we'll use pT to calculate the factor for mass
	      // see https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetWtagging
	      float c = getSmearingFactor(mSDSF, mSDSFUnc, mSDResolutionRel, fatJet, *genJetsAK8, 1, 0.4, 99999., true);
	      jet_mass_softdrop_PUPPI_JERUp = c*jet_mass_softdrop_PUPPI;
	      smearedJetUp = fatJetUp.p4();
	      smearedJetUp_SD.SetPt(fatJetUp.pt());
	      smearedJetUp_SD.SetEta(fatJetUp.eta());
	      smearedJetUp_SD.SetPhi(fatJetUp.phi());
	      smearedJetUp_SD.SetM(jet_mass_softdrop_JERUp);
	    }
	  else
	    {
	      jet_pt_JERUp = -99;
	      jet_mass_JERUp = -99;
	      jet_mass_pruned_JERUp = -99;
	      jet_mass_softdrop_JERUp = -99;
	      jet_mass_softdrop_PUPPI_JERUp = -99;
	    }

	  if (jetsSmearedDown->size() > 0)
	    {
	      const pat::Jet & fatJetDown = jetsSmearedDown->at(0);
	      jet_pt_JERDown = fatJetDown.pt();
	      jet_mass_JERDown = fatJetDown.mass();
	      double JERDownCorrection = fatJetDown.pt()/jet_pt;
	      jet_mass_pruned_JERDown = JERDownCorrection*jet_mass_pruned;
	      jet_mass_softdrop_JERDown = JERDownCorrection*jet_mass_softdrop;
	      float c = getSmearingFactor(mSDSF, mSDSFUnc, mSDResolutionRel, fatJet, *genJetsAK8, -1, 0.4, 99999., true);
	      jet_mass_softdrop_PUPPI_JERDown = c*jet_mass_softdrop_PUPPI;
	      smearedJetDown = fatJetDown.p4();
	      smearedJetDown_SD.SetPt(fatJetDown.pt());
	      smearedJetDown_SD.SetEta(fatJetDown.eta());
	      smearedJetDown_SD.SetPhi(fatJetDown.phi());
	      smearedJetDown_SD.SetM(jet_mass_softdrop_JERDown);
	    }
	  else
	    {
	      jet_pt_JERDown = -99;
	      jet_mass_JERDown = -99;
	      jet_mass_pruned_JERDown = -99;
	      jet_mass_softdrop_JERDown = -99;
	      jet_mass_softdrop_PUPPI_JERDown = -99;
	    }

	}
    }
  
  else throw cms::Exception("InvalidValue") << "This shouldn't happen, we require at least 1 jet, but the size of the jet collection for this event is zero!" << std::endl; 
  
  //Loop over the collection of the AK4 jets which contain b-tagging information (to veto b-jets)
  if(isMC)
    {
      jetFlavours.clear();

      BgenjetStatus43_pt.clear();
      BgenjetStatus43_eta.clear();
      BgenjetStatus43_phi.clear();
      BgenjetStatus43_mass.clear();
      BgenjetStatus43_motherPDGID.clear();


      BgenjetStatus21_pt.clear();
      BgenjetStatus21_eta.clear();
      BgenjetStatus21_phi.clear();
      BgenjetStatus21_mass.clear();
      BgenjetStatus21_motherPDGID.clear();
  
      for (unsigned int iGen = 0; iGen < genParticles->size() && isMC; ++iGen)
	{
	  if (std::abs((genParticles->at(iGen)).pdgId()) == 5)
	    {
	      if ((genParticles -> at(iGen)).status() == 21)
		{
		  BgenjetStatus21_pt.push_back((genParticles -> at(iGen)).pt());
		  BgenjetStatus21_eta.push_back((genParticles -> at(iGen)).eta());
		  BgenjetStatus21_phi.push_back((genParticles -> at(iGen)).phi());
		  BgenjetStatus21_mass.push_back((genParticles -> at(iGen)).mass());
		  BgenjetStatus21_motherPDGID.push_back((genParticles -> at(iGen)).mother()->pdgId());
		}
	      if ((genParticles -> at(iGen)).status() == 43)
		{
		  BgenjetStatus43_pt.push_back((genParticles -> at(iGen)).pt());
		  BgenjetStatus43_eta.push_back((genParticles -> at(iGen)).eta());
		  BgenjetStatus43_phi.push_back((genParticles -> at(iGen)).phi());
		  BgenjetStatus43_mass.push_back((genParticles -> at(iGen)).mass());
		  BgenjetStatus43_motherPDGID.push_back((genParticles -> at(iGen)).mother()->pdgId());
		}
	    }
	}
    }
  // Count number of btagged AK4 jets
  nbtag = 0;
  nbtagMedium = 0;
  nbtagLoose = 0;
  nbtag_JECUp = 0;
  nbtag_JECDown = 0;
  nbtag_JERUp = 0;
  nbtag_JERDown = 0;
  for (const auto & itr : *AK4Jets)
    {
      //taken from: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80X#Supported_Algorithms_and_Operati
      if((itr.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")) > bTagDiscrCut){nbtag ++;}
      if((itr.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")) > bTagDiscrCutMedium){nbtagMedium ++;}
      if((itr.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")) > bTagDiscrCutLoose){nbtagLoose ++;}
      // std::cout << "Nominal " << itr.pt() << " : " << itr.eta() << " : " << itr.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") << std::endl;
      if(isMC) jetFlavours.push_back(itr.partonFlavour());
    }

  if (isMC)
    {
      for (const auto & itr : *AK4JetsSmearedUp)
	{
	  if (itr.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > bTagDiscrCut) {nbtag_JERUp++;}
	  // std::cout << "Smeared Up " << itr.pt() << " : " << itr.eta() << " : " << itr.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") << std::endl;
	}

      for (const auto & itr : *AK4JetsSmearedDown)
	{
	  if (itr.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > bTagDiscrCut) {nbtag_JERDown++;}
	  // std::cout << "Smeared Down " << itr.pt() << " : " << itr.eta() << " : " << itr.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") << std::endl;
	}

      for (const auto & itr : *AK4JetsShiftedUp)
	{
	  if (itr.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > bTagDiscrCut) {nbtag_JECUp++;}
	  // std::cout << "Shifted Up " << itr.pt() << " : " << itr.eta() << " : " << itr.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") << std::endl;
	}

      for (const auto & itr : *AK4JetsShiftedDown)
	{
	  if (itr.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > bTagDiscrCut) {nbtag_JECDown++;}
	  // std::cout << "Shifted Down " << itr.pt() << " : " << itr.eta() << " : " << itr.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") << std::endl;
	}
    }
 
  if (AK4Jets -> size() > 0)
    {
      jet2_pt = (AK4Jets -> at(0)).pt();
      jet2_eta = (AK4Jets -> at(0)).phi();
      jet2_phi = (AK4Jets -> at(0)).eta();
      jet2_btag = (AK4Jets -> at(0)).bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    }
  
  else 
    {
      jet2_pt = -99.;
      jet2_eta = -99.;
      jet2_phi = -99.;
      jet2_btag = -99.;
    }
  
  if (AK4Jets -> size() > 1)
    {
      jet3_pt = (AK4Jets -> at(1)).pt();
      jet3_eta = (AK4Jets -> at(1)).eta();
      jet3_phi = (AK4Jets -> at(1)).phi();
      jet3_btag = (AK4Jets -> at(1)).bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    } 
  else 
    {
      jet3_pt = -99.;
      jet3_eta = -99.;
      jet3_phi = -99.;
      jet3_btag = -99.;
    }
  
  //diboson mass
  ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > hadronicVp4, hadronicVp4_SD, leptonicVp4, lvj_p4, lvj_p4_SD;
  //hadronic V
  hadronicVp4.SetPt(jet_pt);
  hadronicVp4.SetEta(jet_eta);
  hadronicVp4.SetPhi(jet_phi);
  hadronicVp4.SetM(jet_mass);

  // hadronic V but with puppi SD mass
  hadronicVp4_SD.SetPt(jet_pt);
  hadronicVp4_SD.SetEta(jet_eta);
  hadronicVp4_SD.SetPhi(jet_phi);
  hadronicVp4_SD.SetM(jet_mass_softdrop_PUPPI);

  //std::cout<<"*Hadronically decaying W* "<<hadronicVp4.px()<<" "<<hadronicVp4.py()<<" "<<hadronicVp4.pz()<<" "<<hadronicVp4.e()<<" "<<hadronicVp4.M()<<std::endl;

  //leptonic W
  leptonicVp4.SetPt(Wboson_lep.pt);
  leptonicVp4.SetEta(Wboson_lep.eta);
  leptonicVp4.SetPhi(Wboson_lep.phi);
  leptonicVp4.SetM(Wboson_lep.mass);
  //std::cout<<"*Leptonically decaying W* "<<leptonicVp4.px()<<" "<<leptonicVp4.py()<<" "<<leptonicVp4.pz()<<" "<<leptonicVp4.e()<<" "<<leptonicVp4.M()<<std::endl;
   
  if (leptonicVs -> size() > 0)
    {
      saveDibosonMass(leptonicVp4, hadronicVp4, m_lvj);
      saveDibosonMass(leptonicVp4, hadronicVp4_SD, m_lvj_SD);
    }
  else
    {
      m_lvj = -99.;
      m_lvj_SD = -99.;
    }
  //std::cout<<"Diboson mass: "<<m_lvj<<std::endl<<std::endl;
  //std::cout<<"---------------------------------------"<<std::endl;
  //systematics
  if (isMC)
    {
      m_lvj_UnclEnUp = -99.;
      m_lvj_UnclEnDown = -99.;

      m_lvj_JECUp = -99.;
      m_lvj_JECDown = -99.;

      m_lvj_LeptonEnUp = -99.;
      m_lvj_LeptonEnDown = -99.;

      m_lvj_LeptonResUp = -99.;
      m_lvj_LeptonResDown = -99.;

      m_lvj_JERUp = -99.;
      m_lvj_JERDown = -99.;

      // SD versions
      m_lvj_SD_UnclEnUp = -99.;
      m_lvj_SD_UnclEnDown = -99.;

      m_lvj_SD_JECUp = -99.;
      m_lvj_SD_JECDown = -99.;

      m_lvj_SD_LeptonEnUp = -99.;
      m_lvj_SD_LeptonEnDown = -99.;

      m_lvj_SD_LeptonResUp = -99.;
      m_lvj_SD_LeptonResDown = -99.;

      m_lvj_SD_JERUp = -99.;
      m_lvj_SD_JERDown = -99.;

      if (leptonicVs -> size() > 0)
	{
	  //METUnclEn
	  saveDibosonMass(SystMap.at("UnclusteredEnUp"), hadronicVp4, m_lvj_UnclEnUp);
	  saveDibosonMass(SystMap.at("UnclusteredEnDown"), hadronicVp4, m_lvj_UnclEnDown);

	  saveDibosonMass(SystMap.at("UnclusteredEnUp"), hadronicVp4_SD, m_lvj_SD_UnclEnUp);
	  saveDibosonMass(SystMap.at("UnclusteredEnDown"), hadronicVp4_SD, m_lvj_SD_UnclEnDown);

	  //JEC
	  saveDibosonMass(SystMap.at("JetEnUp"), hadronicVp4*(1+JECunc), m_lvj_JECUp);
	  saveDibosonMass(SystMap.at("JetEnDown"), hadronicVp4*(1+JECunc), m_lvj_JECDown);

	  saveDibosonMass(SystMap.at("JetEnUp"), hadronicVp4_SD*(1+JECunc), m_lvj_SD_JECUp);
	  saveDibosonMass(SystMap.at("JetEnDown"), hadronicVp4_SD*(1+JECunc), m_lvj_SD_JECDown);

	  //lepton energy scale uncertainty
	  saveDibosonMass(SystMap.at("LeptonEnUp"), hadronicVp4, m_lvj_LeptonEnUp);
	  saveDibosonMass(SystMap.at("LeptonEnDown"), hadronicVp4, m_lvj_LeptonEnDown);

	  saveDibosonMass(SystMap.at("LeptonEnUp"), hadronicVp4_SD, m_lvj_SD_LeptonEnUp);
	  saveDibosonMass(SystMap.at("LeptonEnDown"), hadronicVp4_SD, m_lvj_SD_LeptonEnDown);

	  //lepton energy resolution uncertainty
	  saveDibosonMass(SystMap.at("LeptonResUp"), hadronicVp4, m_lvj_LeptonResUp);
	  saveDibosonMass(SystMap.at("LeptonResDown"), hadronicVp4, m_lvj_LeptonResDown);

	  saveDibosonMass(SystMap.at("LeptonResUp"), hadronicVp4_SD, m_lvj_SD_LeptonResUp);
	  saveDibosonMass(SystMap.at("LeptonResDown"), hadronicVp4_SD, m_lvj_SD_LeptonResDown);
	}

      //jet energy resolution uncertainty
      if (leptonicVs -> size() > 0 && jetsSmearedUp -> size() > 0)
	{
	  saveDibosonMass(SystMap.at("JetResUp"), smearedJetUp, m_lvj_JERUp);
	  saveDibosonMass(SystMap.at("JetResUp"), smearedJetUp_SD, m_lvj_SD_JERUp);
	}
      if (leptonicVs -> size() > 0 && jetsSmearedUp -> size() > 0)
	{
	  saveDibosonMass(SystMap.at("JetResDown"), smearedJetDown, m_lvj_JERDown);
	  saveDibosonMass(SystMap.at("JetResDown"), smearedJetDown_SD, m_lvj_SD_JERDown);
	}
    }

  edm::Handle<edm::TriggerResults> Triggers;
  if (channel == "el")
    {
      iEvent.getByToken(TriggerResultsToken, Triggers); 
      edm::TriggerNames names = iEvent.triggerNames(*Triggers);
      for (unsigned int iTrig = 0; iTrig < Triggers -> size(); iTrig ++)
	{
	  if( boost::algorithm::contains(names.triggerName(iTrig), "HLT_Ele105_CaloIdVT_GsfTrkIdT_v") ) bit_HLT_Ele_105 = Triggers -> accept(iTrig);
	  if( boost::algorithm::contains(names.triggerName(iTrig), "HLT_Ele27_WPLoose_Gsf_v") )  bit_HLT_Ele_27 =  Triggers -> accept(iTrig);
	  if( boost::algorithm::contains(names.triggerName(iTrig), "HLT_Ele45_WPLoose_Gsf_v") )  bit_HLT_Ele_45 =  Triggers -> accept(iTrig);
	  if( boost::algorithm::contains(names.triggerName(iTrig), "HLT_Ele115_CaloIdVT_GsfTrkIdT_v") )  bit_HLT_Ele_115 =  Triggers -> accept(iTrig);
	  if( boost::algorithm::contains(names.triggerName(iTrig), "HLT_Ele30_WPTight_Gsf_v") )  bit_HLT_Ele_30 =  Triggers -> accept(iTrig);
	  if( boost::algorithm::contains(names.triggerName(iTrig), "HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v") )  bit_HLT_Ele_50_Jet_165 =  Triggers -> accept(iTrig);
	}
    }
  if (bit_HLT_Ele_115 && bit_HLT_Ele_50_Jet_165) bit_BOTH_115_165 = 1;

  topPtSF=1.;
  bool topInGen=0, antitopInGen=0;
  double topSF=1., antitopSF=1.;
  //Double_t gen_neutrino_pz = 0.0;
  found = 0;
  if(isMC)
    {
      for(unsigned int iGen = 0; iGen < genParticles->size(); ++iGen)
	{
	  if((genParticles->at(iGen)).pdgId()==6 && (genParticles->at(iGen)).status()==22)
	    {
	      topInGen=1;
	      topSF=std::exp(0.0615-(0.0005*(genParticles->at(iGen)).pt()));
	    }
	  else if((genParticles->at(iGen)).pdgId()==-6 && (genParticles->at(iGen)).status()==22)
	    {
	      antitopInGen=1;
	      antitopSF=std::exp(0.0615-(0.0005*(genParticles->at(iGen)).pt()));
	    }
	  //if ((std::abs(genParticles->at(iGen).pdgId()) == 12 || std::abs(genParticles->at(iGen).pdgId()) == 14) && (genParticles->at(iGen).mother()->pdgId() == 23 || std::abs(genParticles->at(iGen).mother()->pdgId()) == 24 ))
	  if ((std::abs(genParticles->at(iGen).pdgId()) == 12 || std::abs(genParticles->at(iGen).pdgId()) == 14)) // && (genParticles->at(iGen).mother()->pdgId() == 23 || std::abs(genParticles->at(iGen).mother()->pdgId()) == 24 ))
	    {
	      //here is where we set the gen neutrino stuff
	      gen_neutrino_pz = genParticles->at(iGen).pz();
	      found = 1;
	    }
	}
      if(topInGen && antitopInGen)
	topPtSF = std::sqrt(topSF * antitopSF); 
      //	std::cout<<topSF<<" "<<antitopSF<<" "<<topPtSF<<std::endl;
    }

  double genWeightPosForaTGC=std::abs(genWeight);
  if(isMC && !isSignal) genWeightPosForaTGC=genWeight;

  if (isMC)
    {
      totWeight = PUweight*genWeightPosForaTGC*aTGCWeightUnitConv*LeptonSF*btagWeight*VTagSF*topPtSF;
      totWeight_BTagUp = PUweight*genWeightPosForaTGC*aTGCWeightUnitConv*LeptonSF*btagWeight_BTagUp*VTagSF*topPtSF;
      totWeight_BTagDown = PUweight*genWeightPosForaTGC*aTGCWeightUnitConv*LeptonSF*btagWeight_BTagDown*VTagSF*topPtSF;
      totWeight_MistagUp = PUweight*genWeightPosForaTGC*aTGCWeightUnitConv*LeptonSF*btagWeight_MistagUp*VTagSF*topPtSF;
      totWeight_MistagDown = PUweight*genWeightPosForaTGC*aTGCWeightUnitConv*LeptonSF*btagWeight_MistagDown*VTagSF*topPtSF;
      totWeight_LeptonIDUp = PUweight*genWeightPosForaTGC*aTGCWeightUnitConv*LeptonSF_Up*btagWeight*VTagSF*topPtSF;
      totWeight_LeptonIDDown = PUweight*genWeightPosForaTGC*aTGCWeightUnitConv*LeptonSF_Down*btagWeight*VTagSF*topPtSF;
    }
  //probably would leave it like that if we keep reweighting to data trigger efficiency in electron channel. In this case lepton ID & trigger scale factors are set to unity in the electron channel.
  /*if (isMC&&channel=="el"){
    totWeight *=  triggerWeightHLTEle27NoER;
    totWeight_BTagUp *= triggerWeightHLTEle27NoER;
    totWeight_BTagDown *= triggerWeightHLTEle27NoER;
    totWeight_MistagUp *= triggerWeightHLTEle27NoER; 
    totWeight_MistagDown *= triggerWeightHLTEle27NoER;
    }*/

  // ANGLES
  //   mm   mm   m   mmm  m      mmmmmm  mmmm 
  //   ##   #"m  # m"   " #      #      #"   "
  //  #  #  # #m # #   mm #      #mmmmm "#mmm 
  //  #mm#  #  # # #    # #      #          "#
  // #    # #   ##  "mmm" #mmmmm #mmmmm "mmm#"


  
  Double_t lepton_mass(-99.9), METpx(-99.0), METpy(-99.0);//, METx(0.0), METy(0.0);
  TLorentzVector lepton, met4vector, other_met4vector, jet0, jet1;
  jet0.SetPxPyPzE(zerothsubjet_px, zerothsubjet_py, zerothsubjet_pz, zerothsubjet_e);
  jet1.SetPxPyPzE(firstsubjet_px, firstsubjet_py, firstsubjet_pz, firstsubjet_e);

  /////////////////////////////////////////
  Double_t lepton_px(0.0), lepton_py(0.0), lepton_pz(0.0);
  
  if (channel == "el")
    {
      lepton_mass = 0.00051099;
      lepton_px = (leptons->at(0).px());
      lepton_py = (leptons->at(0).py());
      lepton_pz = (leptons->at(0).pz());
    }
   
  if (channel == "mu")
    {
      auto leptonPtr = leptons -> ptrAt(0);
      reco::MuonPtr asmuonPtr(leptonPtr);
      reco::ElectronPtr aselectronPtr(leptonPtr);
      lepton_mass = 0.105658367;
      lepton_px = asmuonPtr->tunePMuonBestTrack()->px();
      lepton_py = asmuonPtr->tunePMuonBestTrack()->py();
      lepton_pz = asmuonPtr->tunePMuonBestTrack()->pz();
    }

  
  if (metHandle->size() > 0)
    {
      const pat::MET& metCand = metHandle->at(0);

      METpx = metCand.px();
      METpy = metCand.py();
      //METx = metCand.px();
      //METy = metCand.py();
   }
  
  lepton.SetXYZM(lepton_px,lepton_py,lepton_pz,lepton_mass);
  /////////////////////////////////////////
    
  //METpx = METCand.pt * TMath::Cos(METCand.phi);
  //METpy = METCand.pt * TMath::Sin(METCand.phi);
  
  imaginary_neutrino = 0;
 
  //if (channel == "el") lepton_mass = 0.00051099;
  //if (channel == "mu") lepton_mass = 0.105658367;
  
 
  //lepton.SetPtEtaPhiM(Lepton.pt,Lepton.eta,Lepton.phi,lepton_mass);
  
  //lepton_xerror = lepton.Px() - lepton_px;//obsolete
  //lepton_yerror = lepton.Py() - lepton_py;
  //lepton_zerror = lepton.Pz() - lepton_pz;

  //met_xerror = METx - METpx;//obsolete
  //met_yerror = METy - METpy;

  METzCalculator calc;
  calc.SetMET(METpx, METpy);
  calc.SetLepton(lepton);
  calc.SetTruthInfo(gen_neutrino_pz);
  calc.SetLeptonType(channel);//technically already have lepton mass from 4-vector, but this is still needed
  calc.SetJets(jet0, jet1);

  Double_t met_pz = calc.Calculate(7);//currently supports values 0-7
  Double_t other_met_pz = calc.getOther();//must be run after Calculate method

  if (calc.IsComplex()) imaginary_neutrino = 1;

  if (imaginary_neutrino == 1 && METcorrect)
    {
      METxyCorrector corr;
      corr.setVocal(false);
      corr.initialMET(METpx, METpy);
      corr.setLepton(lepton);
      //corr.SetLeptonType(channel);

      METpx = corr.Correct(1);
      METpy = corr.Correct(2);
      Double_t newWmass = corr.Correct(0);

      calc.SetMET(METpx, METpy);
      calc.SetWmass(newWmass);
      met_pz = calc.Calculate(7);
    }

  met4vector.SetXYZM(METpx,METpy,met_pz,0.0);
  other_met4vector.SetXYZM(METpx,METpy,other_met_pz,0.0);
  
  if (Wboson_lep.charge < 0.0) calculateAngles(lepton, met4vector, jet0, jet1, d_costheta1, d_costheta2, d_phi, d_costhetastar, d_phi1, d_phi2);
  if (Wboson_lep.charge > 0.0) calculateAngles(met4vector, lepton, jet0, jet1, d_costheta1, d_costheta2, d_phi, d_costhetastar, d_phi1, d_phi2);
  if (Wboson_lep.charge == 0.0) {std::cout << "There's a zero charge W?"; exit(5);}

  if (metHandle->size() > 0 && number_of_subjets >= 1)
    {    
      costhetastar = d_costhetastar;
      costheta1 = d_costheta1;
      costheta2 = d_costheta2;
      phi = d_phi;
      phi1 = d_phi1;
      phi2 = d_phi2;
    }

  else
    {
      costhetastar = -99.9;
      costheta1 = -99.9;
      costheta2 = -99.9;
      phi = +99.9;
      phi1 = +99.9;
      phi2 = +99.9;
      return;
    }

  delta_neutrino0 =  calc.Calculate(0) - gen_neutrino_pz;
  delta_neutrino1 =  calc.Calculate(1) - gen_neutrino_pz;
  delta_neutrino2 =  calc.Calculate(2) - gen_neutrino_pz;
  delta_neutrino3 =  calc.Calculate(3) - gen_neutrino_pz;
  delta_neutrino4 =  calc.Calculate(4) - gen_neutrino_pz;
  delta_neutrino5 =  calc.Calculate(5) - gen_neutrino_pz;
  delta_neutrino6 =  calc.Calculate(6) - gen_neutrino_pz;
  delta_neutrino7 =  calc.Calculate(7) - gen_neutrino_pz;
  delta_neutrino8 =  calc.Calculate(8) - gen_neutrino_pz;

#ifdef ANGLE_TESTING  
   //INTERMEDIATE STEPS
  // 
  // MARKER: INTERMEDIATE
  //   _       _                               _ _       _       
  //  (_)_ __ | |_ ___ _ __ _ __ ___   ___  __| (_) __ _| |_ ___ 
  //  | | '_ \| __/ _ \ '__| '_ ` _ \ / _ \/ _` | |/ _` | __/ _ \ |
  //  | | | | | ||  __/ |  | | | | | |  __/ (_| | | (_| | ||  __/
  //  |_|_| |_|\__\___|_|  |_| |_| |_|\___|\__,_|_|\__,_|\__\___|
  //
    
 
   if (Wboson_lep.charge > 0.0) intermediate_steps(lepton, met4vector, jet0, jet1, imaginary_neutrino, d_leptons_in_lep_px, d_leptons_in_lep_py, d_leptons_in_lep_pz, d_partons_in_lep_px, d_partons_in_lep_py, d_partons_in_lep_pz, d_parton1_in_lep_px, d_parton2_in_lep_px, d_parton1_in_lep_py, d_parton2_in_lep_py, d_parton1_in_lep_pz, d_parton2_in_lep_pz, d_lepton1_in_lep_px, d_lepton1_in_lep_py, d_lepton1_in_lep_pz, d_lepton1_dotted_x, d_lepton1_dotted_y, d_lepton1_dotted_z, d_leptons_in_had_px, d_leptons_in_had_py, d_leptons_in_had_pz, d_lepton1_in_had_px, d_lepton1_in_had_py, d_lepton1_in_had_pz, d_lepton2_in_had_px, d_lepton2_in_had_py, d_lepton2_in_had_pz, d_parton1_in_had_px, d_parton1_in_had_py, d_parton1_in_had_pz, d_parton1_dotted_x, d_parton1_dotted_y, d_parton1_dotted_z, d_complicated1_px, d_complicated1_py, d_complicated1_pz, d_complicated2_px, d_complicated2_py, d_complicated2_pz, d_lepton_sumWWframe_X, d_lepton_sumWWframe_Y, d_lepton_sumWWframe_Z, d_lepton1WWframe_X, d_lepton1WWframe_Y, d_lepton1WWframe_Z, d_parton_sumWWframe_X, d_parton_sumWWframe_Y, d_parton_sumWWframe_Z, d_parton1WWframe_X, d_parton1WWframe_Y, d_parton1WWframe_Z, d_costhetastar, d_costheta1, d_phi, d_costheta2, d_phi1, d_phi2, d_boostWWframe_X, d_boostWWframe_Y, d_boostWWframe_Z, d_boostWlep_X, d_boostWlep_Y, d_boostWlep_Z, d_boostWhad_X, d_boostWhad_Y, d_boostWhad_Z, d_xdotx, d_xdoty, d_xdotz, d_ydotx, d_ydoty, d_ydotz, d_zdotx, d_zdoty, d_zdotz, d_lepton1WWframe_UX, d_lepton1WWframe_UY, d_lepton1WWframe_UZ, d_lepton_sumWWframe_UX, d_lepton_sumWWframe_UY, d_lepton_sumWWframe_UZ, d_leptons_in_lep_px_good, d_leptons_in_lep_py_good, d_leptons_in_lep_pz_good, d_partons_in_lep_px_good, d_partons_in_lep_py_good, d_partons_in_lep_pz_good, d_parton1_in_lep_px_good, d_parton2_in_lep_px_good, d_parton1_in_lep_py_good, d_parton2_in_lep_py_good, d_parton1_in_lep_pz_good, d_parton2_in_lep_pz_good, d_lepton1_in_lep_px_good, d_lepton1_in_lep_py_good, d_lepton1_in_lep_pz_good, d_leptons_in_had_px_good, d_leptons_in_had_py_good, d_leptons_in_had_pz_good, d_lepton1_in_had_px_good, d_lepton1_in_had_py_good, d_lepton1_in_had_pz_good, d_lepton2_in_had_px_good, d_lepton2_in_had_py_good, d_lepton2_in_had_pz_good, d_parton1_in_had_px_good, d_parton1_in_had_py_good, d_parton1_in_had_pz_good, d_lepton_sumWWframe_X_good, d_lepton_sumWWframe_Y_good, d_lepton_sumWWframe_Z_good, d_lepton1WWframe_X_good, d_lepton1WWframe_Y_good, d_lepton1WWframe_Z_good, d_parton_sumWWframe_X_good, d_parton_sumWWframe_Y_good, d_parton_sumWWframe_Z_good, d_parton1WWframe_X_good, d_parton1WWframe_Y_good, d_parton1WWframe_Z_good);
      
   if (Wboson_lep.charge < 0.0) intermediate_steps(met4vector, lepton, jet0, jet1, imaginary_neutrino, d_leptons_in_lep_px, d_leptons_in_lep_py, d_leptons_in_lep_pz, d_partons_in_lep_px, d_partons_in_lep_py, d_partons_in_lep_pz, d_parton1_in_lep_px, d_parton2_in_lep_px, d_parton1_in_lep_py, d_parton2_in_lep_py, d_parton1_in_lep_pz, d_parton2_in_lep_pz, d_lepton1_in_lep_px, d_lepton1_in_lep_py, d_lepton1_in_lep_pz, d_lepton1_dotted_x, d_lepton1_dotted_y, d_lepton1_dotted_z, d_leptons_in_had_px, d_leptons_in_had_py, d_leptons_in_had_pz, d_lepton1_in_had_px, d_lepton1_in_had_py, d_lepton1_in_had_pz, d_lepton2_in_had_px, d_lepton2_in_had_py, d_lepton2_in_had_pz, d_parton1_in_had_px, d_parton1_in_had_py, d_parton1_in_had_pz, d_parton1_dotted_x, d_parton1_dotted_y, d_parton1_dotted_z, d_complicated1_px, d_complicated1_py, d_complicated1_pz, d_complicated2_px, d_complicated2_py, d_complicated2_pz, d_lepton_sumWWframe_X, d_lepton_sumWWframe_Y, d_lepton_sumWWframe_Z, d_lepton1WWframe_X, d_lepton1WWframe_Y, d_lepton1WWframe_Z, d_parton_sumWWframe_X, d_parton_sumWWframe_Y, d_parton_sumWWframe_Z, d_parton1WWframe_X, d_parton1WWframe_Y, d_parton1WWframe_Z, d_costhetastar, d_costheta1, d_phi, d_costheta2, d_phi1, d_phi2, d_boostWWframe_X, d_boostWWframe_Y, d_boostWWframe_Z, d_boostWlep_X, d_boostWlep_Y, d_boostWlep_Z, d_boostWhad_X, d_boostWhad_Y, d_boostWhad_Z, d_xdotx, d_xdoty, d_xdotz, d_ydotx, d_ydoty, d_ydotz, d_zdotx, d_zdoty, d_zdotz, d_lepton1WWframe_UX, d_lepton1WWframe_UY, d_lepton1WWframe_UZ, d_lepton_sumWWframe_UX, d_lepton_sumWWframe_UY, d_lepton_sumWWframe_UZ, d_leptons_in_lep_px_good, d_leptons_in_lep_py_good, d_leptons_in_lep_pz_good, d_partons_in_lep_px_good, d_partons_in_lep_py_good, d_partons_in_lep_pz_good, d_parton1_in_lep_px_good, d_parton2_in_lep_px_good, d_parton1_in_lep_py_good, d_parton2_in_lep_py_good, d_parton1_in_lep_pz_good, d_parton2_in_lep_pz_good, d_lepton1_in_lep_px_good, d_lepton1_in_lep_py_good, d_lepton1_in_lep_pz_good, d_leptons_in_had_px_good, d_leptons_in_had_py_good, d_leptons_in_had_pz_good, d_lepton1_in_had_px_good, d_lepton1_in_had_py_good, d_lepton1_in_had_pz_good, d_lepton2_in_had_px_good, d_lepton2_in_had_py_good, d_lepton2_in_had_pz_good, d_parton1_in_had_px_good, d_parton1_in_had_py_good, d_parton1_in_had_pz_good, d_lepton_sumWWframe_X_good, d_lepton_sumWWframe_Y_good, d_lepton_sumWWframe_Z_good, d_lepton1WWframe_X_good, d_lepton1WWframe_Y_good, d_lepton1WWframe_Z_good, d_parton_sumWWframe_X_good, d_parton_sumWWframe_Y_good, d_parton_sumWWframe_Z_good, d_parton1WWframe_X_good, d_parton1WWframe_Y_good, d_parton1WWframe_Z_good);
    
   leptons_in_lep_px = d_leptons_in_lep_px;
   leptons_in_lep_py = d_leptons_in_lep_py;
   leptons_in_lep_pz = d_leptons_in_lep_pz;
  
   partons_in_lep_px = d_partons_in_lep_px;
   partons_in_lep_py = d_partons_in_lep_py;
   partons_in_lep_pz = d_partons_in_lep_pz;
  
   parton1_in_lep_px = d_parton1_in_lep_px;
   parton1_in_lep_py = d_parton1_in_lep_py;
   parton1_in_lep_pz = d_parton1_in_lep_pz;
 
   parton2_in_lep_px = d_parton2_in_lep_px;
   parton2_in_lep_py = d_parton2_in_lep_py;
   parton2_in_lep_pz = d_parton2_in_lep_pz;
  
   lepton1_in_lep_px = d_lepton1_in_lep_px;
   lepton1_in_lep_py = d_lepton1_in_lep_py;
   lepton1_in_lep_pz = d_lepton1_in_lep_pz;
  
   lepton1_dotted_x = d_lepton1_dotted_x;
   lepton1_dotted_y = d_lepton1_dotted_y;
   lepton1_dotted_z = d_lepton1_dotted_z;
  
   leptons_in_had_px = d_leptons_in_had_px;
   leptons_in_had_py = d_leptons_in_had_py;
   leptons_in_had_pz = d_leptons_in_had_pz;
  
   lepton1_in_had_px = d_lepton1_in_had_px;
   lepton1_in_had_py = d_lepton1_in_had_py;
   lepton1_in_had_pz = d_lepton1_in_had_pz;
  
   lepton2_in_had_px = d_lepton2_in_had_px;
   lepton2_in_had_py = d_lepton2_in_had_py;
   lepton2_in_had_pz = d_lepton2_in_had_pz;
  
   parton1_in_had_px = d_parton1_in_had_px;
   parton1_in_had_py = d_parton1_in_had_py;
   parton1_in_had_pz = d_parton1_in_had_pz;
  
   parton1_dotted_x = d_parton1_dotted_x;
   parton1_dotted_y = d_parton1_dotted_y;
   parton1_dotted_z = d_parton1_dotted_z;
  
   complicated1_px = d_complicated1_px;
   complicated1_py = d_complicated1_py;
   complicated1_pz = d_complicated1_pz;
  
   complicated2_px = d_complicated2_px;
   complicated2_py = d_complicated2_py;
   complicated2_pz = d_complicated2_pz;
  
   lepton_sumWWframe_X = d_lepton_sumWWframe_X;
   lepton_sumWWframe_Y = d_lepton_sumWWframe_Y;
   lepton_sumWWframe_Z = d_lepton_sumWWframe_Z;
  
   lepton1WWframe_X = d_lepton1WWframe_X;
   lepton1WWframe_Y = d_lepton1WWframe_Y;
   lepton1WWframe_Z = d_lepton1WWframe_Z;
  
   parton_sumWWframe_X = d_parton_sumWWframe_X;
   parton_sumWWframe_Y = d_parton_sumWWframe_Y;
   parton_sumWWframe_Z = d_parton_sumWWframe_Z;
  
   parton1WWframe_X = d_parton1WWframe_X;
   parton1WWframe_Y = d_parton1WWframe_Y;
   parton1WWframe_Z = d_parton1WWframe_Z;
  
   boostWWframe_X = d_boostWWframe_X;
   boostWWframe_Y = d_boostWWframe_Y;
   boostWWframe_Z = d_boostWWframe_Z;
  
   boostWlep_X = d_boostWlep_X;
   boostWlep_Y = d_boostWlep_Y;
   boostWlep_Z = d_boostWlep_Z;
  
   boostWhad_X = d_boostWhad_X;
   boostWhad_Y = d_boostWhad_Y;
   boostWhad_Z = d_boostWhad_Z;

   xdotx = d_xdotx;
   xdoty = d_xdoty;
   xdotz = d_xdotz;

   ydotx = d_ydotx;
   ydoty = d_ydoty;
   ydotz = d_ydotz;

   zdotx = d_zdotx;
   zdoty = d_zdoty;
   zdotz = d_zdotz;

   lepton1WWframe_UX = d_lepton1WWframe_UX;
   lepton1WWframe_UY = d_lepton1WWframe_UY;
   lepton1WWframe_UZ = d_lepton1WWframe_UZ;
   
   lepton_sumWWframe_UX = d_lepton_sumWWframe_UX;
   lepton_sumWWframe_UY = d_lepton_sumWWframe_UY;
   lepton_sumWWframe_UZ = d_lepton_sumWWframe_UZ;
    
   ///////////////////////////////////////////////////////////////////////

   leptons_in_lep_px_good = d_leptons_in_lep_px_good;
   leptons_in_lep_py_good = d_leptons_in_lep_py_good;
   leptons_in_lep_pz_good = d_leptons_in_lep_pz_good;
  
   partons_in_lep_px_good = d_partons_in_lep_px_good;
   partons_in_lep_py_good = d_partons_in_lep_py_good;
   partons_in_lep_pz_good = d_partons_in_lep_pz_good;
  
   parton1_in_lep_px_good = d_parton1_in_lep_px_good;
   parton1_in_lep_py_good = d_parton1_in_lep_py_good;
   parton1_in_lep_pz_good = d_parton1_in_lep_pz_good;

   parton2_in_lep_px_good = d_parton2_in_lep_px_good;
   parton2_in_lep_py_good = d_parton2_in_lep_py_good;
   parton2_in_lep_pz_good = d_parton2_in_lep_pz_good;

   lepton1_in_lep_px_good = d_lepton1_in_lep_px_good;
   lepton1_in_lep_py_good = d_lepton1_in_lep_py_good;
   lepton1_in_lep_pz_good = d_lepton1_in_lep_pz_good;
  
   //lepton1_dotted_x_good = d_lepton1_dotted_x_good;
   //lepton1_dotted_y_good = d_lepton1_dotted_y_good;
   //lepton1_dotted_z_good = d_lepton1_dotted_z_good;
  
   leptons_in_had_px_good = d_leptons_in_had_px_good;
   leptons_in_had_py_good = d_leptons_in_had_py_good;
   leptons_in_had_pz_good = d_leptons_in_had_pz_good;
  
   lepton1_in_had_px_good = d_lepton1_in_had_px_good;
   lepton1_in_had_py_good = d_lepton1_in_had_py_good;
   lepton1_in_had_pz_good = d_lepton1_in_had_pz_good;
  
   lepton2_in_had_px_good = d_lepton2_in_had_px_good;
   lepton2_in_had_py_good = d_lepton2_in_had_py_good;
   lepton2_in_had_pz_good = d_lepton2_in_had_pz_good;
  
   parton1_in_had_px_good = d_parton1_in_had_px_good;
   parton1_in_had_py_good = d_parton1_in_had_py_good;
   parton1_in_had_pz_good = d_parton1_in_had_pz_good;
  
   //parton1_dotted_x_good = d_parton1_dotted_x_good;
   //parton1_dotted_y_good = d_parton1_dotted_y_good;
   //parton1_dotted_z_good = d_parton1_dotted_z_good;
  
   //complicated1_px_good = d_complicated1_px_good;
   //complicated1_py_good = d_complicated1_py_good;
   //complicated1_pz_good = d_complicated1_pz_good;
  
   //complicated2_px_good = d_complicated2_px_good;
   //complicated2_py_good = d_complicated2_py_good;
   //complicated2_pz_good = d_complicated2_pz_good;
  
   lepton_sumWWframe_X_good = d_lepton_sumWWframe_X_good;
   lepton_sumWWframe_Y_good = d_lepton_sumWWframe_Y_good;
   lepton_sumWWframe_Z_good = d_lepton_sumWWframe_Z_good;
  
   lepton1WWframe_X_good = d_lepton1WWframe_X_good;
   lepton1WWframe_Y_good = d_lepton1WWframe_Y_good;
   lepton1WWframe_Z_good = d_lepton1WWframe_Z_good;
  
   parton_sumWWframe_X_good = d_parton_sumWWframe_X_good;
   parton_sumWWframe_Y_good = d_parton_sumWWframe_Y_good;
   parton_sumWWframe_Z_good = d_parton_sumWWframe_Z_good;
  
   parton1WWframe_X_good = d_parton1WWframe_X_good;
   parton1WWframe_Y_good = d_parton1WWframe_Y_good;
   parton1WWframe_Z_good = d_parton1WWframe_Z_good;
  
   //boostWWframe_X_good = d_boostWWframe_X_good;
   //boostWWframe_Y_good = d_boostWWframe_Y_good;
   //boostWWframe_Z_good = d_boostWWframe_Z_good;
  
   //boostWlep_X_good = d_boostWlep_X_good;
   //boostWlep_Y_good = d_boostWlep_Y_good;
   //boostWlep_Z_good = d_boostWlep_Z_good;
  
   //boostWhad_X_good = d_boostWhad_X_good;
   //boostWhad_Y_good = d_boostWhad_Y_good;
   //boostWhad_Z_good = d_boostWhad_Z_good;
   //*/ 
#endif
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  if(isSignal) {if(aTGCWeights.size() >= weights_number) outTree_->Fill();}
  else {outTree_->Fill();}

}

float TreeMaker::getPUPPIweight(float puppipt, float puppieta)
{

  float genCorr  = 1.;
  float recoCorr = 1.;
  float totalWeight = 1.;
        
  genCorr =  puppisd_corrGEN->Eval( puppipt );
  if( fabs(puppieta)  <= 1.3 ) {recoCorr = puppisd_corrRECO_cen->Eval( puppipt );}
  else {recoCorr = puppisd_corrRECO_for->Eval( puppipt );}

  if(isMC) totalWeight= genCorr * recoCorr;
  else totalWeight= recoCorr;

  return totalWeight;
}


float TreeMaker::getSmearingFactor(float sf, float unc, float resolution, const pat::Jet & jet, const edm::View<reco::GenJet> & genJets, int variation, float drMax, float relResMax, bool usePuppiPt) 
{
  // Calculate smearing factor using hybrid method
  // https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
  // i.e. rescale if match with gen-level object,
  // otherwise use stochastic method with scaling drawn from Gaussian
  // does matching with GenJets, using closest in dR, given some maximum dR
  // and relative resolution criteria
  //
  // sf is the scale factor, unc is its uncertainty
  // resolution is the relative pt resolution in simulation
  // jet is the reco jet
  // genJets are the collection of genJets to determine if there's a matching genJet
  // variation is whether to do nominal(0), up (1), or down (-1)
  // drMax is the maximum dR to count as a match with a genJet
  // relResMax is the maximum relative reoslution to count as a match
  // usePuppiPt is a bool to control whether to use Puppi Pt or normal pt for the reco jet
  if (!(variation==0 || abs(variation)==1)) {throw std::runtime_error("variation must be 0 (nominal) or +/-1");}

  float jet_pt = usePuppiPt ? jet.userFloat("ak8PFJetsPuppiValueMap:pt") : jet.pt();

  // First find if there's a match
  float ptGen = -1.;
  float dRBest = 9999;
  for (const auto & itr: genJets)
    {
      float dR = deltaR(jet, itr);
      float relRes = fabs(jet_pt - itr.pt())/jet_pt;
      if (dR < drMax && relRes < relResMax && dR < dRBest)
	{
	  dRBest = dR;
	  ptGen = itr.pt();
	}
    }
  // Now calc factor
  float this_sf = sf + (variation * unc);
  if (ptGen >= 0)
    {
      // scaling method
      // std::cout << "match" << std::endl;
      // std::cout << jet.pt() << " : " << jet.eta() << " : " << jet.phi() << std::endl;
      // std::cout << ptGen << " : " << etaGen << " : " << phiGen << std::endl;
      return 1 + ((this_sf-1)*(1-(ptGen/jet_pt)));
    }
  else
    {
      // std::cout << "no match" << std::endl;
      // stochastic method
      // initialise seed with reproducible number
      TRandom rand((int)(1000*jet.eta()));
      float random_gauss = rand.Gaus(0, resolution);
      return 1 + random_gauss * (sqrt(std::max((this_sf*this_sf) - 1, 0.0f)));
    }
}

void TreeMaker::saveDibosonMass(const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > & leptonicV_p4, const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > & hadronicV_p4, double & massVar)
{
  auto lvj_p4 = leptonicV_p4 + hadronicV_p4;
  massVar = lvj_p4.M();
}

void TreeMaker::saveDibosonMass(math::XYZTLorentzVector & leptonicV_p4, const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > & hadronicV_p4, double & massVar)
{
  auto lvj_p4 = leptonicV_p4 + hadronicV_p4;
  massVar = lvj_p4.M();
}

void TreeMaker::saveDibosonMass(math::XYZTLorentzVector & leptonicV_p4, math::XYZTLorentzVector & hadronicV_p4, double & massVar)
{
  auto lvj_p4 = leptonicV_p4 + hadronicV_p4;
  massVar = lvj_p4.M();
}


bool TreeMaker::decaysHadronic(const reco::Candidate* p)
{
  if (p!=NULL)
    {	
      //cout<<p->pdgId()<<" ";
      if (abs(p->pdgId())<7 || abs(p->pdgId())==21) return true;
      else
	{
	  for(size_t i = 0; i < p->numberOfDaughters(); ++i)
	    {
	      const reco::Candidate* d = (reco::Candidate*)p->daughter(i);
	      if(decaysHadronic(d)) return true;
	    }
	}
    }
  return false;	
}

void TreeMaker::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{


  edm::Handle<LHERunInfoProduct> run; 


  if(!isSignal && isMC) iRun.getByLabel( "externalLHEProducer", run);
  else if (isSignal)iRun.getByLabel( "externalLHEProducer", run);
  else return;
  std::cout << "Nominal : " << run->heprup().PDFSUP.first << std::endl;
  NominalPDF = run->heprup().PDFSUP.first;

}

void TreeMaker::endRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
}
// ------------ method called once each job just before starting event loop  ------------
void TreeMaker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void TreeMaker::endJob() {std::cout << "TreeMaker endJob()..." << std::endl;}

void calculateAngles(TLorentzVector thep4M11, TLorentzVector thep4M12, TLorentzVector thep4M21, TLorentzVector thep4M22, Double_t& costheta1, Double_t& costheta2, Double_t& phi, Double_t& costhetastar, Double_t& phi1, Double_t& phi2)
{


  TLorentzVector thep4H = thep4M11 + thep4M12 + thep4M21 + thep4M22;
  TLorentzVector thep4Z1 = thep4M11 + thep4M12;
  TLorentzVector thep4Z2 = thep4M21 + thep4M22;

  Double_t norm;

  TVector3 boostX = -(thep4H.BoostVector());
  TLorentzVector thep4Z1inXFrame(thep4Z1);
  TLorentzVector thep4Z2inXFrame(thep4Z2);      
  thep4Z1inXFrame.Boost(boostX);
  thep4Z2inXFrame.Boost(boostX);
  TVector3 theZ1X_p3 = TVector3(thep4Z1inXFrame.X(), thep4Z1inXFrame.Y(), thep4Z1inXFrame.Z());
  TVector3 theZ2X_p3 = TVector3(thep4Z2inXFrame.X(), thep4Z2inXFrame.Y(), thep4Z2inXFrame.Z());

  // calculate phi1, phi2, costhetastar
  ///phi1 = theZ1X_p3.Phi();
  ///phi2 = theZ2X_p3.Phi();

  ///////////////////////////////////////////////
  // check for z1/z2 convention, redefine all 4 vectors with convention
  /////////////////////////////////////////////// 
  TLorentzVector p4H, p4Z1, p4M11, p4M12, p4Z2, p4M21, p4M22;
  p4Z1 = thep4Z1; p4M11 = thep4M11; p4M12 = thep4M12;
  p4Z2 = thep4Z2; p4M21 = thep4M21; p4M22 = thep4M22;
  costhetastar = theZ1X_p3.CosTheta();

  // now helicity angles................................
  // ...................................................
  TVector3 boostZ1 = -(p4Z1.BoostVector());
  TLorentzVector p4Z2Z1(p4Z2);
  p4Z2Z1.Boost(boostZ1);
  //find the decay axis
  /////TVector3 unitx_1 = -Hep3Vector(p4Z2Z1);
  TVector3 unitx_1(-p4Z2Z1.X(), -p4Z2Z1.Y(), -p4Z2Z1.Z());
  norm = 1/(unitx_1.Mag());
  unitx_1*=norm;
  //boost daughters of z2
  TLorentzVector p4M21Z1(p4M21);
  TLorentzVector p4M22Z1(p4M22);
  p4M21Z1.Boost(boostZ1);
  p4M22Z1.Boost(boostZ1);
  //create z and y axes
  /////TVector3 unitz_1 = Hep3Vector(p4M21Z1).cross(Hep3Vector(p4M22Z1));
  TVector3 p4M21Z1_p3(p4M21Z1.X(), p4M21Z1.Y(), p4M21Z1.Z());
  TVector3 p4M22Z1_p3(p4M22Z1.X(), p4M22Z1.Y(), p4M22Z1.Z());
  TVector3 unitz_1 = p4M21Z1_p3.Cross(p4M22Z1_p3);
  norm = 1/(unitz_1.Mag());
  unitz_1 *= norm;
  TVector3 unity_1 = unitz_1.Cross(unitx_1);

  //caculate theta1
  TLorentzVector p4M11Z1(p4M11);
  p4M11Z1.Boost(boostZ1);
  TVector3 p3M11(p4M11Z1.X(), p4M11Z1.Y(), p4M11Z1.Z());
  TVector3 unitM11 = p3M11.Unit();
  Double_t x_m11 = unitM11.Dot(unitx_1); Double_t y_m11 = unitM11.Dot(unity_1); Double_t z_m11 = unitM11.Dot(unitz_1);
  TVector3 M11_Z1frame(y_m11, z_m11, x_m11);
  costheta1 = M11_Z1frame.CosTheta();
  //std::cout << "theta1: " << M11_Z1frame.Theta() << std::endl;
  //////-----------------------old way of calculating phi---------------/////////
  phi = M11_Z1frame.Phi();

  //set axes for other system
  TVector3 boostZ2 = -(p4Z2.BoostVector());
  TLorentzVector p4Z1Z2(p4Z1);
  p4Z1Z2.Boost(boostZ2);
  TVector3 unitx_2(-p4Z1Z2.X(), -p4Z1Z2.Y(), -p4Z1Z2.Z());
  norm = 1/(unitx_2.Mag());
  unitx_2 *= norm;
  //boost daughters of z2
  TLorentzVector p4M11Z2(p4M11);
  TLorentzVector p4M12Z2(p4M12);
  p4M11Z2.Boost(boostZ2);
  p4M12Z2.Boost(boostZ2);
  TVector3 p4M11Z2_p3(p4M11Z2.X(), p4M11Z2.Y(), p4M11Z2.Z());
  TVector3 p4M12Z2_p3(p4M12Z2.X(), p4M12Z2.Y(), p4M12Z2.Z());
  TVector3 unitz_2 = p4M11Z2_p3.Cross(p4M12Z2_p3);
  norm = 1/(unitz_2.Mag());
  unitz_2*=norm;
  TVector3 unity_2 = unitz_2.Cross(unitx_2);
  //calcuate theta2
  TLorentzVector p4M21Z2(p4M21);
  p4M21Z2.Boost(boostZ2);
  TVector3 p3M21(p4M21Z2.X(), p4M21Z2.Y(), p4M21Z2.Z());
  TVector3 unitM21 = p3M21.Unit();
  Double_t x_m21 = unitM21.Dot(unitx_2); Double_t y_m21 = unitM21.Dot(unity_2); Double_t z_m21 = unitM21.Dot(unitz_2);
  TVector3 M21_Z2frame(y_m21, z_m21, x_m21);
  costheta2 = M21_Z2frame.CosTheta();

  // calculate phi
  //calculating phi_n
  TLorentzVector n_p4Z1inXFrame(p4Z1);
  TLorentzVector n_p4M11inXFrame(p4M11);
  n_p4Z1inXFrame.Boost(boostX);
  n_p4M11inXFrame.Boost(boostX);        
  TVector3 n_p4Z1inXFrame_unit = n_p4Z1inXFrame.Vect().Unit();
  TVector3 n_p4M11inXFrame_unit = n_p4M11inXFrame.Vect().Unit();  
  TVector3 n_unitz_1(n_p4Z1inXFrame_unit);
  //// y-axis is defined by neg lepton cross z-axis
  //// the subtle part is here...
  //////////TVector3 n_unity_1 = n_p4M11inXFrame_unit.Cross(n_unitz_1);
  TVector3 n_unity_1 = n_unitz_1.Cross(n_p4M11inXFrame_unit);
  TVector3 n_unitx_1 = n_unity_1.Cross(n_unitz_1);

  TLorentzVector n_p4M21inXFrame(p4M21);
  n_p4M21inXFrame.Boost(boostX);
  TVector3 n_p4M21inXFrame_unit = n_p4M21inXFrame.Vect().Unit();
  //rotate into other plane
  TVector3 n_p4M21inXFrame_unitprime(n_p4M21inXFrame_unit.Dot(n_unitx_1), n_p4M21inXFrame_unit.Dot(n_unity_1), n_p4M21inXFrame_unit.Dot(n_unitz_1));

  ///////-----------------new way of calculating phi-----------------///////
  //Double_t phi_n =  n_p4M21inXFrame_unitprime.Phi();
  //
    //std::cout << "---------------------------" << std::endl;
    //std::cout << "phi: " << phi << std::endl;
    //std::cout << "phi_n: " << phi_n << std::endl;
    //std::cout << "phi + phi_n: " << (phi+phi_n) << std::endl;
  //
  /// and then calculate phi1
  TVector3 n_p4PartoninXFrame_unit(0.0, 0.0, 1.0);
  TVector3 n_p4PartoninXFrame_unitprime(n_p4PartoninXFrame_unit.Dot(n_unitx_1), n_p4PartoninXFrame_unit.Dot(n_unity_1), n_p4PartoninXFrame_unit.Dot(n_unitz_1));
  // negative sign is for arrow convention in paper
  phi1 = (n_p4PartoninXFrame_unitprime.Phi());

  // and the calculate phi2
  TLorentzVector n_p4Z2inXFrame(p4Z2);
  n_p4Z2inXFrame.Boost(boostX);
  TVector3 n_p4Z2inXFrame_unit = n_p4Z2inXFrame.Vect().Unit();
  ///////TLorentzVector n_p4M21inXFrame(p4M21);
  //////n_p4M21inXFrame.Boost(boostX);        
  ////TVector3 n_p4M21inXFrame_unit = n_p4M21inXFrame.Vect().Unit();  
  TVector3 n_unitz_2(n_p4Z2inXFrame_unit);
  //// y-axis is defined by neg lepton cross z-axis
  //// the subtle part is here...
  //////TVector3 n_unity_2 = n_p4M21inXFrame_unit.Cross(n_unitz_2);
  TVector3 n_unity_2 = n_unitz_2.Cross(n_p4M21inXFrame_unit);
  TVector3 n_unitx_2 = n_unity_2.Cross(n_unitz_2);
  TVector3 n_p4PartoninZ2PlaneFrame_unitprime(n_p4PartoninXFrame_unit.Dot(n_unitx_2), n_p4PartoninXFrame_unit.Dot(n_unity_2), n_p4PartoninXFrame_unit.Dot(n_unitz_2));
  phi2 = (n_p4PartoninZ2PlaneFrame_unitprime.Phi());

  //
    //Double_t phi12_0 = phi1 + phi2;
    //if (phi12_0 > TMath::Pi()) phi12 = phi12_0 - 2*TMath::Pi();
    //else if (phi12_0 < (-1.)*TMath::Pi()) phi12 = phi12_0 + 2*TMath::Pi();
    //else phi12 = phi12_0;
  //

}

#ifdef ANGLE_TESTING
//INTERMEDIATE STEPS
void intermediate_steps(TLorentzVector lepton1, TLorentzVector lepton2, TLorentzVector parton1, TLorentzVector parton2, Int_t neu_status, Double_t& leptons_in_lep_px, Double_t& leptons_in_lep_py, Double_t& leptons_in_lep_pz, Double_t& partons_in_lep_px, Double_t& partons_in_lep_py, Double_t& partons_in_lep_pz, Double_t& parton1_in_lep_px, Double_t& parton2_in_lep_px, Double_t& parton1_in_lep_py, Double_t& parton2_in_lep_py, Double_t& parton1_in_lep_pz, Double_t& parton2_in_lep_pz, Double_t& lepton1_in_lep_px, Double_t& lepton1_in_lep_py, Double_t& lepton1_in_lep_pz, Double_t& lepton1_dotted_x, Double_t& lepton1_dotted_y, Double_t& lepton1_dotted_z, Double_t& leptons_in_had_px, Double_t& leptons_in_had_py, Double_t& leptons_in_had_pz, Double_t& lepton1_in_had_px, Double_t& lepton1_in_had_py, Double_t& lepton1_in_had_pz, Double_t& lepton2_in_had_px, Double_t& lepton2_in_had_py, Double_t& lepton2_in_had_pz, Double_t& parton1_in_had_px, Double_t& parton1_in_had_py, Double_t& parton1_in_had_pz, Double_t& parton1_dotted_x, Double_t& parton1_dotted_y, Double_t& parton1_dotted_z, Double_t& complicated1_px, Double_t& complicated1_py, Double_t& complicated1_pz, Double_t& complicated2_px, Double_t& complicated2_py, Double_t& complicated2_pz, Double_t& lepton_sumWWframe_X, Double_t& lepton_sumWWframe_Y, Double_t& lepton_sumWWframe_Z, Double_t& lepton1WWframe_X, Double_t& lepton1WWframe_Y, Double_t& lepton1WWframe_Z, Double_t& parton_sumWWframe_X, Double_t& parton_sumWWframe_Y, Double_t& parton_sumWWframe_Z, Double_t& parton1WWframe_X, Double_t& parton1WWframe_Y, Double_t& parton1WWframe_Z, Double_t& costhetastar, Double_t& costheta1, Double_t& phi, Double_t& costheta2, Double_t& phi1, Double_t& phi2, Double_t& boostWWframe_X, Double_t& boostWWframe_Y, Double_t& boostWWframe_Z, Double_t& boostWlep_X, Double_t& boostWlep_Y, Double_t& boostWlep_Z, Double_t& boostWhad_X, Double_t& boostWhad_Y, Double_t& boostWhad_Z, Double_t& xdotx, Double_t& xdoty, Double_t& xdotz, Double_t& ydotx, Double_t& ydoty, Double_t& ydotz, Double_t& zdotx, Double_t& zdoty, Double_t& zdotz, Double_t& lepton1WWframe_UX, Double_t& lepton1WWframe_UY, Double_t& lepton1WWframe_UZ, Double_t& lepton_sumWWframe_UX, Double_t& lepton_sumWWframe_UY, Double_t& lepton_sumWWframe_UZ, Int_t& leptons_in_lep_px_good, Int_t& leptons_in_lep_py_good, Int_t& leptons_in_lep_pz_good, Int_t& partons_in_lep_px_good, Int_t& partons_in_lep_py_good, Int_t& partons_in_lep_pz_good, Int_t& parton1_in_lep_px_good, Int_t& parton2_in_lep_px_good, Int_t& parton1_in_lep_py_good, Int_t& parton2_in_lep_py_good, Int_t& parton1_in_lep_pz_good, Int_t& parton2_in_lep_pz_good, Int_t& lepton1_in_lep_px_good, Int_t& lepton1_in_lep_py_good, Int_t& lepton1_in_lep_pz_good, Int_t& leptons_in_had_px_good, Int_t& leptons_in_had_py_good, Int_t& leptons_in_had_pz_good, Int_t& lepton1_in_had_px_good, Int_t& lepton1_in_had_py_good, Int_t& lepton1_in_had_pz_good, Int_t& lepton2_in_had_px_good, Int_t& lepton2_in_had_py_good, Int_t& lepton2_in_had_pz_good, Int_t& parton1_in_had_px_good, Int_t& parton1_in_had_py_good, Int_t& parton1_in_had_pz_good, Int_t& lepton_sumWWframe_X_good, Int_t& lepton_sumWWframe_Y_good, Int_t& lepton_sumWWframe_Z_good, Int_t& lepton1WWframe_X_good, Int_t& lepton1WWframe_Y_good, Int_t& lepton1WWframe_Z_good, Int_t& parton_sumWWframe_X_good, Int_t& parton_sumWWframe_Y_good, Int_t& parton_sumWWframe_Z_good, Int_t& parton1WWframe_X_good, Int_t& parton1WWframe_Y_good, Int_t& parton1WWframe_Z_good)
{
  TVector3 x_hat(1.0,0.0,0.0);
  TVector3 y_hat(0.0,1.0,0.0);
  TVector3 z_hat(0.0,0.0,1.0);

  bad_angle = false;
  TLorentzVector fermion_sum = lepton1 + lepton2 + parton1 + parton2; //4-vector sum of all 4 particles XXX 
  TLorentzVector lepton_sum = lepton1 + lepton2; //4-vector sum of lepton and neutrino XXX
  TLorentzVector parton_sum = parton1 + parton2; //4-vector sum of partons XXX

  TString neu_output = "FIXTHIS";
  if (neu_status == 0) neu_output = "real";
  if (neu_status == 1) neu_output = "imaginary";

  //std::cout << "lepton sum mass is " << lepton_sum.M() << " and neutrino is " << neu_output << std::endl;
  //std::cout << "parton sum mass is " << parton_sum.M() << std::endl;
  //std::cout << "fermion sum mass is " << fermion_sum.M() << std::endl;

  Double_t norm;

  TVector3 boostWWframe = -(fermion_sum.BoostVector()); //boost to WW rest frame XXX
  boostWWframe_X = boostWWframe.X(); //newer
  boostWWframe_Y = boostWWframe.Y();
  boostWWframe_Z = boostWWframe.Z();

  if (boostWWframe_X != boostWWframe_X)
    {
      bad_angle = true;
      std::cout << "boostWWframe_X is bad.\n";
    }
  if (boostWWframe_Y != boostWWframe_Y)
    {
      bad_angle = true;
      std::cout << "boostWWframe_Y is bad.\n";
  
    }
  if (boostWWframe_Z != boostWWframe_Z)
    {
      bad_angle = true;
      std::cout << "boostWWframe_Z is bad.\n";
    }
  
  TLorentzVector lepton_sum_WWframe(lepton_sum); //4-vector sum of leptons boosted to WW rest frame (below) XXX  
  TLorentzVector parton_sum_WWframe(parton_sum); //4-vector sum of partons boosted to WW rest frame (below) XXX 
  lepton_sum_WWframe.Boost(boostWWframe);
  parton_sum_WWframe.Boost(boostWWframe);
  

  leptons_in_lep_px = lepton_sum_WWframe.X(); //new
  leptons_in_lep_py = lepton_sum_WWframe.Y();
  leptons_in_lep_pz = lepton_sum_WWframe.Z();
  
  if (leptons_in_lep_px != leptons_in_lep_px)
    {
      bad_angle = true;
      leptons_in_lep_px_good = 0;
       std::cout << "leptons_in_lep_px is bad" << std::endl;
    }
  
  if (leptons_in_lep_py != leptons_in_lep_py)
    {
      bad_angle = true;
      leptons_in_lep_py_good = 0;
      std::cout << "leptons_in_lep_py is bad" << std::endl;
    }
     
  if (leptons_in_lep_pz != leptons_in_lep_pz)
    {
      bad_angle = true;
      leptons_in_lep_pz_good = 0;
      std::cout << "leptons_in_lep_pz is bad" << std::endl;   
    }
  
 TVector3 lepton_sum_WWframe3vec = TVector3(leptons_in_lep_px, leptons_in_lep_py, leptons_in_lep_pz); //3-vector sum of leptons in WW frame XXX
  
  //////////////////////////////////////////////////////////////////
  
  // calculate phi1, phi2, costhetastar

  ///////////////////////////////////////////////
  // check for z1/z2 convention, redefine all 4 vectors with convention
  ///////////////////////////////////////////////
  
  //TLorentzVector fermion_sum, lepton_sum, lepton1, lepton2, parton_sum, parton1, parton2;
  
  costhetastar = lepton_sum_WWframe3vec.CosTheta(); //ANGLE THETA* 

  // now helicity angles................................
  // ...................................................
  TVector3 boostWlep = -(lepton_sum.BoostVector()); //boost to leptonic W rest frame XXX
  
  boostWlep_X = boostWlep.X(); //newer
  boostWlep_Y = boostWlep.Y();
  boostWlep_Z = boostWlep.Z();

  if (boostWlep_X != boostWlep_X)
    {
      bad_angle = true;
      std::cout << "boostWlep_X is bad" << std::endl;
    }
  if (boostWlep_Y != boostWlep_Y)
    {
      bad_angle = true;
      std::cout << "boostWlep_Y is bad" << std::endl;
    }
  if (boostWlep_Z != boostWlep_Z)
    {
      bad_angle = true;
      std::cout << "boostWlep_Z is bad" << std::endl;
    }
  
  TLorentzVector parton_sumWlep_frame(parton_sum); //4-vector sum of partons boosted to leptonic W rest frame (below) XXX
  parton_sumWlep_frame.Boost(boostWlep);
  //find the decay axis
 

  partons_in_lep_px = parton_sumWlep_frame.X(); //new
  partons_in_lep_py = parton_sumWlep_frame.Y();
  partons_in_lep_pz = parton_sumWlep_frame.Z();
  
  if (partons_in_lep_px != partons_in_lep_px)
    {
      bad_angle = true;
      partons_in_lep_px_good = 0;
      std::cout << "partons_in_lep_px is bad" << std::endl; 
    }
  if (partons_in_lep_py != partons_in_lep_py)
    {
      bad_angle = true;
      partons_in_lep_py_good = 0;
      std::cout << "partons_in_lep_py is bad" << std::endl;
    }
  if (partons_in_lep_pz != partons_in_lep_pz)
    {
      bad_angle = true;
      partons_in_lep_pz_good = 0;
      std::cout << "partons_in_lep_pz is bad" << std::endl;
    }
  
  TVector3 unitx_1(-partons_in_lep_px, -partons_in_lep_py, -partons_in_lep_pz); //unit 3-vector (below) of partons boosted to leptonic W rest frame XXX
  norm = 1.0 / (unitx_1.Mag());
  unitx_1 *= norm;

  if (unitx_1.Mag() != unitx_1.Mag() || unitx_1.Mag() == 0.0)
    {
      bad_angle = true;
      std::cout << "unit_x_1 is bad." << std::endl;
    }
  
  //boost daughters of z2
  TLorentzVector parton1Wlep_frame(parton1); //first parton boosted (below) to leptonic W rest frame XXX
  TLorentzVector parton2Wlep_frame(parton2); //second parton boosted (below) to leptonic W rest frame XXX
  parton1Wlep_frame.Boost(boostWlep);
  parton2Wlep_frame.Boost(boostWlep);
  
  //create z and y axes

  parton1_in_lep_px = parton1Wlep_frame.X(); //new
  parton2_in_lep_px = parton2Wlep_frame.X();
  
  parton1_in_lep_py = parton1Wlep_frame.Y();
  parton2_in_lep_py = parton2Wlep_frame.Y();
  
  parton1_in_lep_pz = parton1Wlep_frame.Z();
  parton2_in_lep_pz = parton2Wlep_frame.Z();
  
  if (parton1_in_lep_px != parton1_in_lep_px)
    {
       bad_angle = true;
      parton1_in_lep_px_good = 0;
      std::cout << "parton1_in_lep_px is bad." << std::endl;
    }
  if (parton1_in_lep_py != parton1_in_lep_py)
    {
       bad_angle = true;
      parton1_in_lep_py_good = 0;
      std::cout << "parton1_in_lep_py is bad." << std::endl;
    }
  if (parton1_in_lep_pz != parton1_in_lep_pz)
    {
      bad_angle = true;
      parton1_in_lep_pz_good = 0;
      std::cout << "parton1_in_lep_pz is bad." << std::endl;
    }

  if (parton2_in_lep_px != parton2_in_lep_px)
    {
      bad_angle = true;
      parton2_in_lep_px_good = 0;
      std::cout << "parton2_in_lep_px is bad." << std::endl;
    }
  if (parton2_in_lep_py != parton2_in_lep_py)
    {
      parton2_in_lep_py_good = 0;
      std::cout << "parton2_in_lep_px is bad." << std::endl;
    }
  if (parton2_in_lep_pz != parton2_in_lep_pz)
    {
      bad_angle = true;
      parton2_in_lep_pz_good = 0;
      std::cout << "parton2_in_lep_px is bad." << std::endl;
    }
  
  TVector3 parton1Wlep_frame3vec(parton1_in_lep_px, parton1_in_lep_py, parton1_in_lep_pz); //3-vector of first parton in leptonic frame XXX
  TVector3 parton2Wlep_frame3vec(parton2_in_lep_px, parton2_in_lep_py, parton2_in_lep_pz); //3-vector of second parton in leptonic frame XXX

  
  TVector3 unitz_1 = parton1Wlep_frame3vec.Cross(parton2Wlep_frame3vec); //unit 3-vector (below) of cross product of partons in leptonic frame XXX
  norm = 1.0 / (unitz_1.Mag());
  unitz_1 *= norm;

  if (unitz_1.Mag() != unitz_1.Mag() || unitz_1.Mag() == 0.0)
    {
      bad_angle = true;
      std::cout << "unitz_1 is bad.\n";
    }
   
  TVector3 unity_1 = unitz_1.Cross(unitx_1); //cross product of sum of partons with cross of partons, all in Wlep frame: (parton1 + parton2) X (parton1 X parton2) XXX

  if (unity_1.Mag() != unity_1.Mag() || unity_1.Mag() == 0.0)
    {
      bad_angle = true;
      std::cout << "unity_1 is bad.\n";
    }

  ///unit1: x = partons in lep ; z = cross product of partons and leptons ; y = z cross x

  //caculate theta1
  TLorentzVector lepton1Wlep_frame(lepton1); //first lepton boosted (below) to leptonic W rest frame XXX
  lepton1Wlep_frame.Boost(boostWlep);

  lepton1_in_lep_px = lepton1Wlep_frame.X(); // new
  lepton1_in_lep_py = lepton1Wlep_frame.Y();
  lepton1_in_lep_pz = lepton1Wlep_frame.Z();
  
  if (lepton1_in_lep_px != lepton1_in_lep_px) {lepton1_in_lep_px_good = 0; bad_angle = true; std::cout << "lepton1_in_lep_px is bad.\n";}
  if (lepton1_in_lep_py != lepton1_in_lep_py) {lepton1_in_lep_py_good = 0; bad_angle = true; std::cout << "lepton1_in_lep_py is bad.\n";}
  if (lepton1_in_lep_pz != lepton1_in_lep_pz) {lepton1_in_lep_pz_good = 0; bad_angle = true; std::cout << "lepton1_in_lep_pz is bad.\n";}
  
  TVector3 lepton1Wlep_frame3vec(lepton1_in_lep_px, lepton1_in_lep_py, lepton1_in_lep_pz); //3-vector of 1st lepton in leptonic frame XXX
  TVector3 lepton1Wlep_frame_unit3vec = lepton1Wlep_frame3vec.Unit(); //unit 3-vector of 1st lepton in leptonic frame XXX

  Double_t x_m11 = lepton1Wlep_frame_unit3vec.Dot(unitx_1);
  Double_t y_m11 = lepton1Wlep_frame_unit3vec.Dot(unity_1);
  Double_t z_m11 = lepton1Wlep_frame_unit3vec.Dot(unitz_1); //dot products of unit 3-vector of 1st lepton in lepton frame and parton units in lepton frame XXX

  lepton1_dotted_x = x_m11; // new
  lepton1_dotted_y = y_m11;
  lepton1_dotted_z = z_m11;

  if (lepton1_dotted_x != lepton1_dotted_x) {bad_angle = true; std::cout << "lepton1_dotted_x is bad.\n";}
  if (lepton1_dotted_y != lepton1_dotted_y) {bad_angle = true; std::cout << "lepton1_dotted_y is bad.\n";}
  if (lepton1_dotted_z != lepton1_dotted_z) {bad_angle = true; std::cout << "lepton1_dotted_z is bad.\n";}
  
  TVector3 M11_Z1frame(y_m11, z_m11, x_m11); //3-vector made of dot products of 1st lepton XXX


  ///////////////////////////////////////////////
  
  costheta1 = M11_Z1frame.CosTheta(); //ANGLE THETA1
 
  phi = M11_Z1frame.Phi(); ///ANGLE PHI

  //////////////////////////////////////////////////////////////////////// now on to last part

  //set axes for other system
  TVector3 boostWhad = -(parton_sum.BoostVector()); //boost to hadronic W rest frame XXX
  boostWhad_X = boostWhad.X(); //newer
  boostWhad_Y = boostWhad.Y();
  boostWhad_Z = boostWhad.Z();

  if (boostWhad_X != boostWhad_X) {bad_angle = true; std::cout << "boostWhad_X is bad.\n";}
  if (boostWhad_Y != boostWhad_Y) {bad_angle = true; std::cout << "boostWhad_Y is bad.\n";}
  if (boostWhad_Z != boostWhad_Z) {bad_angle = true; std::cout << "boostWhad_Z is bad.\n";}

  TLorentzVector lepton_sumWhad_frame(lepton_sum); //sum of leptonic vectors boosted to hadronic frame (below) XXX
  lepton_sumWhad_frame.Boost(boostWhad);

  leptons_in_had_px = lepton_sumWhad_frame.X(); //new
  leptons_in_had_py = lepton_sumWhad_frame.Y();
  leptons_in_had_pz = lepton_sumWhad_frame.Z();
  
  if (leptons_in_had_px != leptons_in_had_px) {leptons_in_had_px_good = 0; bad_angle = true; std::cout << "leptons_in_had_px is bad.\n";}
  if (leptons_in_had_py != leptons_in_had_py) {leptons_in_had_py_good = 0; bad_angle = true; std::cout << "leptons_in_had_py is bad.\n";}
  if (leptons_in_had_pz != leptons_in_had_pz) {leptons_in_had_pz_good = 0; bad_angle = true; std::cout << "leptons_in_had_pz is bad.\n";}
  
  TVector3 unitx_2(-leptons_in_had_px, -leptons_in_had_py, -leptons_in_had_pz); //unit 3-vector (below) of leptonic/hadronic boost XXX
  
  norm = 1.0 / (unitx_2.Mag());
  unitx_2 *= norm;

  if (unitx_2.Mag() != unitx_2.Mag() || unitx_2.Mag() == 0.0)
    {
      bad_angle = true;
      std::cout << "unitx_2 is bad.\n";
    }
  //boost daughters of z2 //////////////////
  TLorentzVector lepton1Whad_frame(lepton1); //first lepton boosted to hadronic frame XXX
  TLorentzVector lepton2Whad_frame(lepton2); //second lepton boosted to hadronic frame XXX
  lepton1Whad_frame.Boost(boostWhad);
  lepton2Whad_frame.Boost(boostWhad);

  lepton1_in_had_px = lepton1Whad_frame.X(); // new
  lepton2_in_had_px = lepton2Whad_frame.X();
  
  lepton1_in_had_py = lepton1Whad_frame.Y();
  lepton2_in_had_py = lepton2Whad_frame.Y();
  
  lepton1_in_had_pz = lepton1Whad_frame.Z();
  lepton2_in_had_pz = lepton2Whad_frame.Z();

  ////////////////
  if (lepton1_in_had_px != lepton1_in_had_px) {lepton1_in_had_px_good = 0; bad_angle = true; std::cout << "lepton1_in_had_px is bad.\n";}
  if (lepton1_in_had_py != lepton1_in_had_py) {lepton1_in_had_py_good = 0; bad_angle = true; std::cout << "lepton1_in_had_py is bad.\n";}
  if (lepton1_in_had_pz != lepton1_in_had_pz) {lepton1_in_had_pz_good = 0; bad_angle = true; std::cout << "lepton1_in_had_pz is bad.\n";}

  if (lepton2_in_had_px != lepton2_in_had_px) {lepton2_in_had_px_good = 0; bad_angle = true; std::cout << "lepton2_in_had_px is bad.\n";}
  if (lepton2_in_had_py != lepton2_in_had_py) {lepton2_in_had_py_good = 0; bad_angle = true; std::cout << "lepton2_in_had_py is bad.\n";}
  if (lepton2_in_had_pz != lepton2_in_had_pz) {lepton2_in_had_pz_good = 0; bad_angle = true; std::cout << "lepton2_in_had_pz is bad.\n";}
  
  TVector3 lepton1Whad_frame3vec(lepton1_in_had_px, lepton1_in_had_py, lepton1_in_had_pz); //3-vector of 1st lepton in hadronic frame XXX
  TVector3 lepton2Whad_frame3vec(lepton2_in_had_px, lepton2_in_had_py, lepton2_in_had_pz); //3-vector of 2nd lepton in hadronic frame XXX
  
  TVector3 unitz_2 = lepton1Whad_frame3vec.Cross(lepton2Whad_frame3vec); //unit of cross of leptons in hadronic frame XXX
  norm = 1.0 / (unitz_2.Mag());
  unitz_2 *= norm;
  TVector3 unity_2 = unitz_2.Cross(unitx_2); //cross of other units XXX

  if (unitz_2.Mag() != unitz_2.Mag() || unitz_2.Mag() == 0.0)
    {
      bad_angle = true;
      std::cout << "unitz_2 is bad.\n";
    }
  ///unit2: x = lepton sum in had; z = lepton1 cross lepton 2, all in had frame; y = z cross z
  
  //calcuate theta2
  TLorentzVector parton1Whad_frame(parton1); //first parton boosted to hadronic frame XXX
  parton1Whad_frame.Boost(boostWhad);

  parton1_in_had_px = parton1Whad_frame.X(); // new
  parton1_in_had_py = parton1Whad_frame.Y();
  parton1_in_had_pz = parton1Whad_frame.Z();
  
  if (parton1_in_had_px != parton1_in_had_px)
    {
    parton1_in_had_px_good = 0; bad_angle = true; std::cout << "parton1_in_had x is bad.\n";
    }
  if (parton1_in_had_py != parton1_in_had_py)
    {
    parton1_in_had_py_good = 0; bad_angle = true; std::cout << "parton1_in_had y is bad.\n";
    }
  if (parton1_in_had_pz != parton1_in_had_pz)
    {
    parton1_in_had_pz_good = 0; bad_angle = true; std::cout << "parton1_in_had z is bad.\n";
    }
  
  TVector3 parton1Whad_frame3vec(parton1_in_had_px, parton1_in_had_py, parton1_in_had_pz); //3-vector of 1st parton in hadronic frame XXX
  TVector3 parton1Whad_frame_unit3vec = parton1Whad_frame3vec.Unit(); //unit 3-vector of 1st parton in hadronic frame XXX

  
  Double_t x_m21 = parton1Whad_frame_unit3vec.Dot(unitx_2);
  Double_t y_m21 = parton1Whad_frame_unit3vec.Dot(unity_2);
  Double_t z_m21 = parton1Whad_frame_unit3vec.Dot(unitz_2); //dot products of unit 3-vector of 1st parton in hadronic frame and lepton units in parton frame XXX

  parton1_dotted_x = x_m21; // new
  parton1_dotted_y = y_m21;
  parton1_dotted_z = z_m21;

  if (parton1_dotted_x != parton1_dotted_x) {bad_angle = true; std::cout << "parton1 dotted x is bad" << std::endl;}
  if (parton1_dotted_y != parton1_dotted_y) {bad_angle = true; std::cout << "parton1 dotted y is bad" << std::endl;}
  if (parton1_dotted_z != parton1_dotted_z) {bad_angle = true; std::cout << "parton1 dotted z is bad" << std::endl;}
    
  TVector3 M21_Z2frame(y_m21, z_m21, x_m21); //3-vector made of dot products of 1st parton
  costheta2 = M21_Z2frame.CosTheta(); //ANGLE THETA2

  ////////////////////////// now a little bit more 
  
  // calculate phi
  //calculating phi_n
  TLorentzVector lepton_sumWWframe(lepton_sum); //sum of leptons, boosted to WW frame XXX
  TLorentzVector lepton1WWframe(lepton1); //1st lepton, boosted to WW frame XXX
  lepton_sumWWframe.Boost(boostWWframe);
  lepton1WWframe.Boost(boostWWframe);

  
  lepton_sumWWframe_X = lepton_sumWWframe.X(); //new
  lepton_sumWWframe_Y = lepton_sumWWframe.Y();
  lepton_sumWWframe_Z = lepton_sumWWframe.Z();
  
  lepton1WWframe_X = lepton1WWframe.X(); //new
  lepton1WWframe_Y = lepton1WWframe.Y();
  lepton1WWframe_Z = lepton1WWframe.Z();
  
  if (lepton_sumWWframe_X != lepton_sumWWframe_X) {lepton_sumWWframe_X_good = 0; bad_angle = true; std::cout << "lepton sum WWframe x is bad.\n";}
  if (lepton_sumWWframe_Y != lepton_sumWWframe_Y) {lepton_sumWWframe_Y_good = 0; bad_angle = true; std::cout << "lepton sum WWframe y is bad.\n";}
  if (lepton_sumWWframe_Z != lepton_sumWWframe_Z) {lepton_sumWWframe_Z_good = 0; bad_angle = true; std::cout << "lepton sum WWframe z is bad.\n";}

  if (lepton1WWframe_X != lepton1WWframe_X) {lepton1WWframe_X_good = 0; bad_angle = true; std::cout << "lepton1 WWframe x is bad.\n";}
  if (lepton1WWframe_Y != lepton1WWframe_Y) {lepton1WWframe_Y_good = 0; bad_angle = true; std::cout << "lepton1 WWframe y is bad.\n";}
  if (lepton1WWframe_Z != lepton1WWframe_Z) {lepton1WWframe_Z_good = 0; bad_angle = true; std::cout << "lepton1 WWframe z is bad.\n";}
  
  TVector3 lepton_sumWWframe_3vec(lepton_sumWWframe_X, lepton_sumWWframe_Y, lepton_sumWWframe_Z); //new
  TVector3 lepton1WWframe_3vec(lepton1WWframe_X, lepton1WWframe_Y, lepton1WWframe_Z);  //new
  
  //TVector3 lepton_sumWWframe_unit3vec = lepton_sumWWframe.Vect().Unit(); //unit 3-vector of sum of leptons, boosted to WW frame XXX
  TVector3 lepton_sumWWframe_unit3vec = lepton_sumWWframe_3vec.Unit(); //unit 3-vector of sum of leptons, boosted to WW frame XXX

  lepton_sumWWframe_UX = lepton_sumWWframe_unit3vec.X();
  lepton_sumWWframe_UY = lepton_sumWWframe_unit3vec.Y();
  lepton_sumWWframe_UZ = lepton_sumWWframe_unit3vec.Z();
  
  //TVector3 lepton1WWframe_unit3vec = lepton1WWframe.Vect().Unit(); //unit 3-vector of 1st lepton, boosted to WW frame XXX
  TVector3 lepton1WWframe_unit3vec = lepton1WWframe_3vec.Unit(); //unit 3-vector of 1st lepton, boosted to WW frame XXX

  lepton1WWframe_UX = lepton1WWframe_unit3vec.X();
  lepton1WWframe_UY = lepton1WWframe_unit3vec.Y();
  lepton1WWframe_UZ = lepton1WWframe_unit3vec.Z();
  
  TVector3 n_unitz_1(lepton_sumWWframe_unit3vec); //unit 3-vector of sum of leptons, boosted to WW frame XXX (yes, it's there twice)
  //// y-axis is defined by neg lepton cross z-axis
  //// the subtle part is here...
 
  TVector3 n_unity_1 = n_unitz_1.Cross(lepton1WWframe_unit3vec); //cross product of 1st lepton and sum of leptons, in WW frame XXX 
  TVector3 n_unitx_1 = n_unity_1.Cross(n_unitz_1); //cross product lepton sum and above cross product, all in WW frame XXX

  ////unit1n: z = lepton sum; y = z cross lepton1; x = y cross z

  //HERE?
  //Double_t& xdotx, Double_t& xdoty, Double_t& xdotz, Double_t& ydotx, Double_t& ydoty, Double_t& ydotz, Double_t& zdotx, Double_t& zdoty, Double_t& zdotz,
  //d_xdotx, d_xdoty, d_xdotz, d_ydotx, d_ydoty, d_ydotz, d_zdotx, d_zdoty, d_zdotz,

  //Double_t& lepton1WWframe_UX, Double_t& lepton1WWframe_UY, Double_t& lepton1WWframe_UZ, Double_t& lepton_sumWWframe_UX, Double_t& lepton_sumWWframe_UY, Double_t& lepton_sumWWframe_UZ,
  //d_lepton1WWframe_UX, d_lepton1WWframe_UY, d_lepton1WWframe_UZ, d_lepton_sumWWframe_UX, d_lepton_sumWWframe_UY, d_lepton_sumWWframe_UZ, 

  xdotx = n_unitx_1.Dot(x_hat);
  xdoty = n_unitx_1.Dot(y_hat);
  xdotz = n_unitx_1.Dot(z_hat);

  ydotx = n_unity_1.Dot(x_hat);
  ydoty = n_unity_1.Dot(y_hat);
  ydotz = n_unity_1.Dot(z_hat);

  zdotx = n_unitz_1.Dot(x_hat);
  zdoty = n_unitz_1.Dot(y_hat);
  zdotz = n_unitz_1.Dot(z_hat);

  TLorentzVector parton1WWframe(parton1); //1st parton, boosted to WW frame XXX
  parton1WWframe.Boost(boostWWframe);

  parton1WWframe_X = parton1WWframe.X(); //new
  parton1WWframe_Y = parton1WWframe.Y();
  parton1WWframe_Z = parton1WWframe.Z();
  
  if (parton1WWframe_X != parton1WWframe_X) {parton1WWframe_X_good = 0; bad_angle = true; std::cout << "parton1 WWframe x is bad.\n";}
  if (parton1WWframe_Y != parton1WWframe_Y) {parton1WWframe_Y_good = 0; bad_angle = true; std::cout << "parton1 WWframe y is bad.\n";}
  if (parton1WWframe_Z != parton1WWframe_Z) {parton1WWframe_Z_good = 0; bad_angle = true; std::cout << "parton1 WWframe y is bad.\n";}
  
  TVector3 parton1WWframe_3vec(parton1WWframe_X, parton1WWframe_Y, parton1WWframe_Z);  //new
  
  //TVector3 parton1WWframe_unit3vec = parton1WWframe.Vect().Unit(); //unit 3-vector of 1st parton, boosted to WW frame XXX
  TVector3 parton1WWframe_unit3vec = parton1WWframe_3vec.Unit(); //unit 3-vector of 1st parton, boosted to WW frame XXX
  //rotate into other plane
 

  ///////-----------------new way of calculating phi-----------------///////
  
  /// and then calculate phi1
  //TVector3 z_hat(0.0, 0.0, 1.0);
  TVector3 z_component_n1_units(z_hat.Dot(n_unitx_1), z_hat.Dot(n_unity_1), z_hat.Dot(n_unitz_1)); //n_units dotted into z-axis??? XXX PROBLEM
  // negative sign is for arrow convention in paper
  phi1 = (z_component_n1_units.Phi()); //ANGLE PHI1 PROBLEM

  complicated1_px = z_hat.Dot(n_unitx_1); //new
  complicated1_py = z_hat.Dot(n_unity_1);
  complicated1_pz = z_hat.Dot(n_unitz_1);

  // and the calculate phi2 
  TLorentzVector parton_sumWWframe(parton_sum); //sum of partons in WW frame XXX
  parton_sumWWframe.Boost(boostWWframe);

  parton_sumWWframe_X = parton_sumWWframe.X(); //new HERE NEXT
  parton_sumWWframe_Y = parton_sumWWframe.Y();
  parton_sumWWframe_Z = parton_sumWWframe.Z();
  
  if (parton_sumWWframe_X != parton_sumWWframe_X) {parton_sumWWframe_X_good = 0; bad_angle = true; std::cout << "parton sum WWframe x is bad.\n";}
  if (parton_sumWWframe_Y != parton_sumWWframe_Y) {parton_sumWWframe_Y_good = 0; bad_angle = true; std::cout << "parton sum WWframe y is bad.\n";}
  if (parton_sumWWframe_Z != parton_sumWWframe_Z) {parton_sumWWframe_Z_good = 0; bad_angle = true; std::cout << "parton sum WWframe z is bad.\n";}
  
  TVector3 parton_sumWWframe_3vec(parton_sumWWframe_X, parton_sumWWframe_Y, parton_sumWWframe_Z); //new
  
  //TVector3 parton_sumWWframe_unit3vec = parton_sumWWframe.Vect().Unit(); //unit 3-vector of sum of partons in WW frame XXX
  TVector3 parton_sumWWframe_unit3vec = parton_sumWWframe_3vec.Unit(); //unit 3-vector of sum of partons in WW frame XXX 
  
  TVector3 n_unitz_2(parton_sumWWframe_unit3vec); //unit 3-vector of sum of partons in WW frame XXX (yes, it's there twice)
  //// y-axis is defined by neg lepton cross z-axis
  //// the subtle part is here...
 
  TVector3 n_unity_2 = n_unitz_2.Cross(parton1WWframe_unit3vec); //sum of partons crossed into 2nd lepton, in WW frame XXX
  TVector3 n_unitx_2 = n_unity_2.Cross(n_unitz_2); //above cross product crossed into 2nd lepton, all in WW frame XXX

  //unit2n: z = parton sum; y = z cross parton1; x = y cross z

  TVector3 z_component_n2_units(z_hat.Dot(n_unitx_2), z_hat.Dot(n_unity_2), z_hat.Dot(n_unitz_2)); //above wierd thing in units of whatever XXX
  phi2 = (z_component_n2_units.Phi()); //ANGLE PHI2

  complicated2_px = z_hat.Dot(n_unitx_2); //new
  complicated2_py = z_hat.Dot(n_unity_2);
  complicated2_pz = z_hat.Dot(n_unitz_2);

  if (bad_angle == true)
    {
      std::cout << "=======================================================================================================\n";
      std::cout << "fermion_sum is X: " << fermion_sum.X() << " Y: " << fermion_sum.Y() << " Z: " << fermion_sum.Z() << " T: " << fermion_sum.T() << std::endl;
      std::cout << "lepton_sum is X: " << lepton_sum.X() << " Y: " << lepton_sum.Y() << " Z: " << lepton_sum.Z() << " T: " << lepton_sum.T() << std::endl;
      std::cout << "parton_sum is X: " << parton_sum.X() << " Y: " << parton_sum.Y() << " Z: " << parton_sum.Z() << " T: " << parton_sum.T() << std::endl;
      std::cout << "boostWWframe is X: " << boostWWframe.X() << " Y: " << boostWWframe.Y() << " Z: " << boostWWframe.Z() << std::endl;
      std::cout << "parton_sumWlep_frame is X: " << parton_sumWlep_frame.X() << " Y: " << parton_sumWlep_frame.Y() << " Z: " << parton_sumWlep_frame.Z() << " T: " << parton_sumWlep_frame.T() << std::endl;
      std::cout << "lepton_sum_WWframe is X: " << lepton_sum_WWframe.X() << " Y: " << lepton_sum_WWframe.Y() << " Z: " << lepton_sum_WWframe.Z() <<  " T: " << lepton_sum_WWframe.T() << std::endl;
      std::cout << "boostWlep is X: " << boostWlep.X() << " Y: " << boostWlep.Y() << " Z: " << boostWlep.Z() << std::endl;
      std::cout << "parton_sumWlep_frame is X: " << parton_sumWlep_frame.X() << " Y: " << parton_sumWlep_frame.Y() << " Z: " << parton_sumWlep_frame.Z() << "T: " << parton_sumWlep_frame.T() << std::endl;
      std::cout << "unitx_1 is X: " << unitx_1.X() << " Y: " << unitx_1.Y() << " Z: " << unitx_1.Z() << std::endl;
      std::cout << "parton1Wlep_frame is X: " << parton1Wlep_frame.X() << " Y: " << parton1Wlep_frame.Y() << " Z: " << parton1Wlep_frame.Z() << " T: " << parton1Wlep_frame.T() << std::endl;
      std::cout << "parton2Wlep_frame is X: " << parton2Wlep_frame.X() << " Y: " << parton2Wlep_frame.Y() << " Z: " << parton2Wlep_frame.Z() << " T: " << parton2Wlep_frame.T() << std::endl;
      std::cout << "unitz_1 is X: " << unitz_1.X() << " Y: " << unitz_1.Y() << " Z: " << unitx_1.Z() << std::endl;
      std::cout << "unity_1 is X: " << unity_1.X() << " Y: " << unity_1.Y() << " Z: " << unity_1.Z() << std::endl;
      std::cout << "lepton1Wlep_frame is X: " << lepton1Wlep_frame.X() << " Y: " << lepton1Wlep_frame.Y() << " Z: " << lepton1Wlep_frame.Z() << " T: " << lepton1Wlep_frame.T() << std::endl;
      std::cout << "lepton1 dotted is X: " << lepton1_dotted_x << " Y: " << lepton1_dotted_y << " Z: " << lepton1_dotted_z << std::endl;
      std::cout << "boostWhad is X: " << boostWhad.X() << " Y: " << boostWhad.Y() << " Z: " << boostWhad.Z() << std::endl;
      std::cout << "unitx_2 is X: " << unitx_2.X() << " Y: " << unitx_2.Y() << " Z: " << unitx_2.Z() << std::endl;
      std::cout << "lepton1Whad_frame is X: " << lepton1Whad_frame.X() << " Y: " << lepton1Whad_frame.Y() << " Z: " << lepton1Whad_frame.Z() << " T: " << lepton1Whad_frame.T() << std::endl;
      std::cout << "unitz_2 is X: " << unitz_2.X() << " Y: " << unitz_2.Y() << " Z: " << unitx_2.Z() << std::endl;
      std::cout << "unity_2 is X: " << unity_2.X() << " Y: " << unity_2.Y() << " Z: " << unity_2.Z() << std::endl;
      std::cout << "parton1Whad_frame is X: " << parton1Whad_frame.X() << " Y: " << parton1Whad_frame.Y() << " Z: " << parton1Whad_frame.Z() << " T: " << parton1Whad_frame.T() << std::endl;
      std::cout << "parton1 dotted is X: " << parton1_dotted_x << " Y: " << parton1_dotted_y << " Z: " << parton1_dotted_z << std::endl;
      std::cout << "lepton_sum_WWframe is X: " << lepton_sum_WWframe.X() << " Y: " << lepton_sum_WWframe.Y() << " Z: " << lepton_sum_WWframe.Z() <<  " T: " << lepton_sum_WWframe.T() << std::endl;
      std::cout << "lepton1WWframe is X: " << lepton1WWframe.X() << " Y: " << lepton1WWframe.Y() << " Z: " << lepton1WWframe.Z() <<  " T: " << lepton1WWframe.T() << std::endl;
      std::cout << "parton1WWframe is X: " << parton1WWframe.X() << " Y: " << parton1WWframe.Y() << " Z: " << parton1WWframe.Z() <<  " T: " << parton1WWframe.T() << std::endl;
      std::cout << "parton_sum_WWframe is X: " << parton_sum_WWframe.X() << " Y: " << parton_sum_WWframe.Y() << " Z: " << parton_sum_WWframe.Z() <<  " T: " << parton_sum_WWframe.T() << std::endl;
      
      std::cout << "n_unitx_1 is X: " << n_unitx_1.X() << " Y: " << n_unitx_1.Y() << " Z: " << n_unitx_1.Z() << std::endl;
      std::cout << "n_unitz_1 is X: " << n_unitz_1.X() << " Y: " << n_unitz_1.Y() << " Z: " << n_unitx_1.Z() << std::endl;
      std::cout << "n_unity_1 is X: " << n_unity_1.X() << " Y: " << n_unity_1.Y() << " Z: " << n_unity_1.Z() << std::endl;

      std::cout << "n_unitx_2 is X: " << n_unitx_2.X() << " Y: " << n_unitx_2.Y() << " Z: " << n_unitx_2.Z() << std::endl;
      std::cout << "n_unitz_2 is X: " << n_unitz_2.X() << " Y: " << n_unitz_2.Y() << " Z: " << n_unitx_2.Z() << std::endl;
      std::cout << "n_unity_2 is X: " << n_unity_2.X() << " Y: " << n_unity_2.Y() << " Z: " << n_unity_2.Z() << std::endl;
      std::cout << "=======================================================================================================\n";

      
      //std::cout << "something is X: " << something.X() << " Y: " << something.Y() << " Z: " << something.Z() << std::endl;
      //std::cout << "something is X: " << something.X() << " Y: " << something.Y() << " Z: " << something.Z() << " T: " << something.T() << std::endl;
    }
  bad_angle = false;

}
#endif

//define this as a plug-in
DEFINE_FWK_MODULE(TreeMaker);
