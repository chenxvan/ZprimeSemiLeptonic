#include <iostream>
#include <memory>

#include <UHH2/core/include/AnalysisModule.h>
#include <UHH2/core/include/Event.h>
#include <UHH2/core/include/Selection.h>
#include "UHH2/common/include/PrintingModules.h"

#include <UHH2/common/include/CleaningModules.h>
#include <UHH2/common/include/NSelections.h>
#include <UHH2/common/include/LumiSelection.h>
#include <UHH2/common/include/TriggerSelection.h>
#include <UHH2/common/include/JetCorrections.h>
#include <UHH2/common/include/ObjectIdUtils.h>
#include <UHH2/common/include/MuonIds.h>
#include <UHH2/common/include/ElectronIds.h>
#include <UHH2/common/include/JetIds.h>
#include <UHH2/common/include/TopJetIds.h>
#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/common/include/Utils.h>
#include <UHH2/common/include/AdditionalSelections.h>
#include "UHH2/common/include/LuminosityHists.h"
#include <UHH2/common/include/MCWeight.h>
#include <UHH2/common/include/MuonHists.h>
#include <UHH2/common/include/ElectronHists.h>
#include <UHH2/common/include/JetHists.h>
#include <UHH2/common/include/EventHists.h>
#include <UHH2/common/include/TopPtReweight.h>
#include <UHH2/common/include/CommonModules.h>

#include <UHH2/ZprimeSemiLeptonic/include/ModuleBASE.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicSelections.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicModules.h>
#include <UHH2/ZprimeSemiLeptonic/include/TTbarLJHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicGeneratorHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeCandidate.h>

//#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/common/include/TTbarReconstruction.h>
#include <UHH2/common/include/ReconstructionHypothesisDiscriminators.h>

using namespace std;
using namespace uhh2;

/*
██████  ███████ ███████ ██ ███    ██ ██ ████████ ██  ██████  ███    ██
██   ██ ██      ██      ██ ████   ██ ██    ██    ██ ██    ██ ████   ██
██   ██ █████   █████   ██ ██ ██  ██ ██    ██    ██ ██    ██ ██ ██  ██
██   ██ ██      ██      ██ ██  ██ ██ ██    ██    ██ ██    ██ ██  ██ ██
██████  ███████ ██      ██ ██   ████ ██    ██    ██  ██████  ██   ████
*/

class ZprimeBDTApply : public ModuleBASE {

public:
  explicit ZprimeBDTApply(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;
  void book_histograms(uhh2::Context&, vector<string>);
  void fill_histograms(uhh2::Event&, string);

protected:

  bool debug;
  
  // Cleaners (to make sure the pre-selection did everything right)
  std::unique_ptr<MuonCleaner> muon_cleaner;
  std::unique_ptr<ElectronCleaner> electron_cleaner;

  // Scale Factors -- Systematics
  unique_ptr<MCMuonScaleFactor> MuonID_module, MuonTrigger_module;
  //unique_ptr<MCElecScaleFactor> EleID_module, EleTrigger_module;

  // AnalysisModules
  unique_ptr<AnalysisModule> LumiWeight_module, PUWeight_module, printer_genparticles, BTagWeight_module, TopPtReweight_module, MCScale_module;

  // Taggers
  unique_ptr<AK8PuppiTopTagger> TopTaggerPuppi;

  // Mass reconstruction
  unique_ptr<ZprimeCandidateBuilder> CandidateBuilder;

  // Chi2 discriminator
  unique_ptr<ZprimeChi2Discriminator> Chi2DiscriminatorZprime;
  unique_ptr<ZprimeCorrectMatchDiscriminator> CorrectMatchDiscriminatorZprime;

  // Selections
  unique_ptr<Selection> Trigger1_selection, Trigger2_selection, NMuon1_selection, NMuon2_selection, NElectron_selection, TwoDCut_selection, Jet1_selection, Jet2_selection, Met_selection, Chi2_selection, TTbarMatchable_selection, Chi2CandidateMatched_selection, ZprimeTopTag_selection, BlindData_selection;
  std::unique_ptr<uhh2::Selection> met_sel;
  std::unique_ptr<uhh2::Selection> htlep_sel;
  std::unique_ptr<Selection> sel_1btag, sel_2btag;
  std::unique_ptr<Selection> TopJetBtagSubjet_selection;
  std::unique_ptr<Selection> HEM_selection;

  // NN variables handles
  //  unique_ptr<Variables_NN> Variables_module;

  //Handles
  Event::Handle<bool> h_is_zprime_reconstructed_chi2, h_is_zprime_reconstructed_correctmatch;
  Event::Handle<float> h_chi2;   Event::Handle<float> h_weight;
  Event::Handle<float> h_MET;   Event::Handle<int> h_NPV;
  Event::Handle<float> h_lep1_pt; Event::Handle<float> h_lep1_eta; 
  Event::Handle<float> h_ak4jet1_pt; Event::Handle<float> h_ak4jet1_eta; 
  Event::Handle<float> h_ak8jet1_pt; Event::Handle<float> h_ak8jet1_eta; 
  Event::Handle<float> h_Mttbar; 

  uhh2::Event::Handle<ZprimeCandidate*> h_BestZprimeCandidateChi2;

  // Configuration
  bool isMC, ispuppi, islooserselection;
  string Sys_MuonID, Sys_MuonTrigger, Sys_PU, Sys_btag; //Sys_EleID, Sys_EleTrigger
  TString sample;
  int runnr_oldtriggers = 299368;

  bool is2016v2, is2016v3, is2017v2, is2018;
  bool isMuon, isElectron;

  //BDT variables                                                                                                                   
  Event::Handle<int> h_ttagN;
  Event::Handle<int> h_btagN;
  Event::Handle<float> h_N_jets;
  Event::Handle<float> h_mass_jet;
  Event::Handle<float> h_mass_jet1;
  Event::Handle<float> h_mass_jet2;
  Event::Handle<float> h_ST;
  Event::Handle<float> h_csv_jet;
  Event::Handle<float> h_csv_jet1;
  Event::Handle<float> h_csv_jet2;
  Event::Handle<float> h_S33;
  Event::Handle<float> h_dR_mu_jet;
  Event::Handle<float> h_dRmin_mu_jet;
  Event::Handle<float> h_dRmin_mu_jet_scaled;
  Event::Handle<float> h_dR_ele_jet;
  Event::Handle<float> h_dRmin_ele_jet;
  Event::Handle<float> h_dRmin_ele_jet_scaled;
  Event::Handle<float> h_ptrel_mu_jet;
  Event::Handle<float> h_ptrel_ele_jet;
  Event::Handle<float> h_deepjet_topscore_jet;
  Event::Handle<float> h_deepjet_topscore_jet1;
  Event::Handle<float> h_deepjet_topscore_jet2;
  Event::Handle<float> h_deepjet_wscore_jet;
  Event::Handle<float> h_deepjet_wscore_jet1;
  Event::Handle<float> h_deepjet_wscore_jet2;
  Event::Handle<float> h_deepjet_zscore_jet;
  Event::Handle<float> h_deepjet_zscore_jet1;
  Event::Handle<float> h_deepjet_zscore_jet2;
  Event::Handle<float> h_deepjet_higgsscore_jet;
  Event::Handle<float> h_deepjet_higgsscore_jet1;
  Event::Handle<float> h_deepjet_higgsscore_jet2;
  Event::Handle<float> h_deepjet_qcdscore_jet;
  Event::Handle<float> h_deepjet_qcdscore_jet1;
  Event::Handle<float> h_deepjet_qcdscore_jet2;


};

void ZprimeBDTApply::book_histograms(uhh2::Context& ctx, vector<string> tags){
  for(const auto & tag : tags){
    string mytag = tag + "_Skimming";
    mytag = tag+"_General";
    book_HFolder(mytag, new ZprimeSemiLeptonicHists(ctx,mytag));
  }
}

void ZprimeBDTApply::fill_histograms(uhh2::Event& event, string tag){
  string mytag = tag + "_Skimming";
  mytag = tag+"_General";
  HFolder(mytag)->fill(event);
}

/*
█  ██████  ██████  ███    ██ ███████ ████████ ██████  ██    ██  ██████ ████████  ██████  ██████
█ ██      ██    ██ ████   ██ ██         ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
█ ██      ██    ██ ██ ██  ██ ███████    ██    ██████  ██    ██ ██         ██    ██    ██ ██████
█ ██      ██    ██ ██  ██ ██      ██    ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
█  ██████  ██████  ██   ████ ███████    ██    ██   ██  ██████   ██████    ██     ██████  ██   ██
*/

ZprimeBDTApply::ZprimeBDTApply(uhh2::Context& ctx){
  //  debug = true;
  debug = false;
  for(auto & kv : ctx.get_all()){
    cout << " " << kv.first << " = " << kv.second << endl;
  }
  // Configuration
  isMC = (ctx.get("dataset_type") == "MC");
  ispuppi = (ctx.get("is_puppi") == "true");
  TString mode = "chs";
  if(ispuppi) mode = "puppi";
  string tmp = ctx.get("dataset_version");
  sample = tmp;
  is2016v2 = (ctx.get("dataset_version").find("2016v2") != std::string::npos);
  is2016v3 = (ctx.get("dataset_version").find("2016v3") != std::string::npos);
  is2017v2 = (ctx.get("dataset_version").find("2017v2") != std::string::npos);
  is2018 = (ctx.get("dataset_version").find("2018") != std::string::npos);

  // Important selection values
  islooserselection = (ctx.get("is_looser_selection") == "true");
  double muon_pt(55.);
  double elec_pt(80.);
  double jet1_pt(150.);
  double jet2_pt(50.);
  double chi2_max(30.);
  double mtt_blind(3000.);
  int nmuon_min1, nmuon_max1;
  int nmuon_min2, nmuon_max2;
  int nele_min, nele_max;
  string trigger1,trigger2;
  double MET_cut, HT_lep_cut;
  isMuon = false; isElectron = false;
  if(ctx.get("channel") == "muon") isMuon = true;
  if(ctx.get("channel") == "electron") isElectron = true;

  if(isMuon){//semileptonic muon channel
    trigger1 = "HLT_Mu50_v*";
    if(is2016v2 || is2016v3)
      //trigger2 = "HLT_TkMu50_v*";
      trigger2 = "HLT_Mu50_v*"; //TkMu path does not exist in 2017/2018 and RunB 2016
    else
      trigger2 = "HLT_Mu50_v*"; //TkMu path does not exist in 2017/2018
    nmuon_min1 = 1, nmuon_max1 = -1;
    nmuon_min2 = 1, nmuon_max2 = 1;
    nele_min = 0; nele_max = 0;
    MET_cut = 50;
    jet1_pt = 150.;
    HT_lep_cut = 150;
  }
  if(isElectron){//semileptonic electron channel
    nmuon_min1 = 0; nmuon_max1 = 0;
    nmuon_min2 = 0; nmuon_max2 = 0;
    nele_min = 1; nele_max = 1;
    trigger1 = "HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v*";
    trigger2 = "HLT_Ele115_CaloIdVT_GsfTrkIdT_v*";
    MET_cut = 120;
    jet1_pt = 185.;
    HT_lep_cut = 0;
  }

 
  double TwoD_dr = 0.4, TwoD_ptrel = 25.;
  if(islooserselection){
    jet1_pt = 100.;
    TwoD_dr = 0.2;
    TwoD_ptrel = 10.;
  }
  const MuonId muonID(PtEtaCut(muon_pt, 2.4));
  const ElectronId electronID(PtEtaSCCut(elec_pt, 2.5));

  

  Sys_MuonID = ctx.get("Sys_MuonID");
  Sys_MuonTrigger = ctx.get("Sys_MuonTrigger");
//  Sys_EleID = ctx.get("Sys_EleID");
//  Sys_EleTrigger = ctx.get("Sys_EleTrigger");
  Sys_PU = ctx.get("Sys_PU");
  //  Sys_btag = ctx.get("Sys_BTagSF");

  BTag::algo btag_algo = BTag::DEEPJET;
  BTag::wp btag_wp_tight = BTag::WP_TIGHT;
  JetId id_btag = BTag(btag_algo, btag_wp_tight);

  double a_toppt = 0.0615; // par a TopPt Reweighting
  double b_toppt = -0.0005; // par b TopPt Reweighting 

  // Modules
  printer_genparticles.reset(new GenParticlesPrinter(ctx));
  muon_cleaner.reset(new MuonCleaner(muonID));
  electron_cleaner.reset(new ElectronCleaner(electronID));
  LumiWeight_module.reset(new MCLumiWeight(ctx));
  PUWeight_module.reset(new MCPileupReweight(ctx, Sys_PU));
  //BTagWeight_module.reset(new MCBTagDiscriminantReweighting(ctx, btag_algo, "jets", Sys_btag));
  TopPtReweight_module.reset(new TopPtReweight(ctx, a_toppt, b_toppt));
  MCScale_module.reset(new MCScaleVariation(ctx));


  if((is2016v3 || is2016v2) && isMuon){
    MuonID_module.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/deleokse/analysis/CMSSW_10_2_10/src/UHH2/common/data/2016/MuonID_EfficienciesAndSF_average_RunBtoH.root", "NUM_TightID_DEN_genTracks_eta_pt", 0., "MuonID", true, Sys_MuonID));
    MuonTrigger_module.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/deleokse/analysis/CMSSW_10_2_10/src/UHH2/common/data/2016/MuonTrigger_EfficienciesAndSF_average_RunBtoH.root", "IsoMu50_OR_IsoTkMu50_PtEtaBins", 0.5, "MuonTrigger", true, Sys_MuonTrigger));
  }
  if(is2017v2 && isMuon){
    MuonID_module.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/deleokse/RunII_102X_v2/CMSSW_10_2_17/src/UHH2/common/data/2017/MuonID_94X_RunBCDEF_SF_ID.root", "NUM_HighPtID_DEN_genTracks_pair_newTuneP_probe_pt_abseta", 0., "HighPtID", true, Sys_MuonID));
    MuonTrigger_module.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/deleokse/RunII_102X_v2/CMSSW_10_2_17/src/UHH2/common/data/2017/MuonTrigger_EfficienciesAndSF_RunBtoF_Nov17Nov2017.root", "Mu50_PtEtaBins/pt_abseta_ratio", 0.5, "Trigger", true, Sys_MuonTrigger));
  }
  if(is2018 && isMuon){
    MuonID_module.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/deleokse/RunII_102X_v2/CMSSW_10_2_17/src/UHH2/common/data/2018/Muon_ID_SF_RunABCD.root", "NUM_HighPtID_DEN_TrackerMuons_pair_newTuneP_probe_pt_abseta", 0., "HighPtID", true, Sys_MuonID));
    MuonTrigger_module.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/deleokse/RunII_102X_v2/CMSSW_10_2_17/src/UHH2/common/data/2018/Muon_Trigger_Eff_SF_AfterMuonHLTUpdate.root", "Mu50_OR_OldMu100_OR_TkMu100_PtEtaBins/pt_abseta_ratio", 0.5, "Trigger", true, Sys_MuonTrigger));
//    EleID_module.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/deleokse/RunII_102X_v2/CMSSW_10_2_16/src/UHH2/common/data/2018/2018_ElectronTight.root", 0., "TightID", Sys_EleID));
//    EleTrigger_module.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/deleokse/RunII_102X_v2/CMSSW_10_2_16/src/UHH2/common/data/2018/SF_2018.root", 0.5, "Trigger", Sys_EleTrigger));
  }

  // Selection modules
  Trigger1_selection.reset(new TriggerSelection(trigger1));
  Trigger2_selection.reset(new TriggerSelection(trigger2));
  NMuon1_selection.reset(new NMuonSelection(nmuon_min1, nmuon_max1));
  NMuon2_selection.reset(new NMuonSelection(nmuon_min2, nmuon_max2));
  NElectron_selection.reset(new NElectronSelection(nele_min, nele_max));
  TwoDCut_selection.reset(new TwoDCut(TwoD_dr, TwoD_ptrel));
  Jet1_selection.reset(new NJetSelection(1, -1, JetId(PtEtaCut(jet1_pt, 2.4))));
  Jet2_selection.reset(new NJetSelection(2, -1, JetId(PtEtaCut(jet2_pt, 2.4))));
  met_sel.reset(new METCut  (MET_cut   , uhh2::infinity));
  htlep_sel.reset(new HTlepCut(HT_lep_cut, uhh2::infinity));

  Chi2_selection.reset(new Chi2Cut(ctx, 0., chi2_max));
  TTbarMatchable_selection.reset(new TTbarSemiLepMatchableSelection());
  Chi2CandidateMatched_selection.reset(new Chi2CandidateMatchedSelection(ctx));
  ZprimeTopTag_selection.reset(new ZprimeTopTagSelection(ctx));
  BlindData_selection.reset(new BlindDataSelection(ctx, mtt_blind));

  HEM_selection.reset(new HEMSelection(ctx)); // HEM issue in 2018, veto on leptons and jets

  //  Variables_module.reset(new Variables_NN(ctx)); // variables for NN

  // Taggers
  TopTaggerPuppi.reset(new AK8PuppiTopTagger(ctx));

  // Zprime candidate builder
  CandidateBuilder.reset(new ZprimeCandidateBuilder(ctx, mode));

  // Zprime discriminators
  Chi2DiscriminatorZprime.reset(new ZprimeChi2Discriminator(ctx));
  h_is_zprime_reconstructed_chi2 = ctx.get_handle<bool>("is_zprime_reconstructed_chi2");
  CorrectMatchDiscriminatorZprime.reset(new ZprimeCorrectMatchDiscriminator(ctx));
  h_is_zprime_reconstructed_correctmatch = ctx.get_handle<bool>("is_zprime_reconstructed_correctmatch");
  h_BestZprimeCandidateChi2 = ctx.get_handle<ZprimeCandidate*>("ZprimeCandidateBestChi2");
  h_chi2 = ctx.declare_event_output<float> ("rec_chi2");
  h_MET = ctx.declare_event_output<float> ("met_pt");
  h_Mttbar = ctx.declare_event_output<float> ("Mttbar");
  h_lep1_pt = ctx.declare_event_output<float> ("lep1_pt");
  h_lep1_eta = ctx.declare_event_output<float> ("lep1_eta");
  h_ak4jet1_pt = ctx.declare_event_output<float> ("ak4jet1_pt");
  h_ak4jet1_eta = ctx.declare_event_output<float> ("ak4jet1_eta");
  h_ak8jet1_pt = ctx.declare_event_output<float> ("ak8jet1_pt");
  h_ak8jet1_eta = ctx.declare_event_output<float> ("ak8jet1_eta");

  h_NPV = ctx.declare_event_output<int> ("NPV");
  h_weight = ctx.declare_event_output<float> ("weight");

  sel_1btag.reset(new NJetSelection(1, 1, id_btag));
  sel_2btag.reset(new NJetSelection(2,-1, id_btag));


  TopJetBtagSubjet_selection.reset(new ZprimeBTagFatSubJetSelection(ctx));


  //variables for limit                                                                                                             
  h_ttagN = ctx.declare_event_output<int>("ttagN");
  h_btagN = ctx.declare_event_output<int>("btagN");
  //  h_wtagN = ctx.declare_event_output<int>("wtagN");


  //BDT variables                                                                                                                   
  h_N_jets = ctx.declare_event_output<float> ("N_jets");
  h_mass_jet = ctx.declare_event_output<float> ("mass_jet");
  h_mass_jet1 = ctx.declare_event_output<float> ("mass_jet1");
  h_mass_jet2 = ctx.declare_event_output<float> ("mass_jet2");
  h_ST = ctx.declare_event_output<float> ("st");
  h_csv_jet = ctx.declare_event_output<float> ("csv_jet");
  h_csv_jet1 = ctx.declare_event_output<float> ("csv_jet1");
  h_csv_jet2 = ctx.declare_event_output<float> ("csv_jet2");
  h_S33 = ctx.declare_event_output<float> ("s33");
  h_dR_mu_jet = ctx.declare_event_output<float> ("dR_mu_jet");
  h_dRmin_mu_jet = ctx.declare_event_output<float> ("dRmin_mu_jet");
  h_dRmin_mu_jet_scaled = ctx.declare_event_output<float> ("dRmin_mu_jet_scaled");
  h_dR_ele_jet = ctx.declare_event_output<float> ("dR_ele_jet");
  h_dRmin_ele_jet = ctx.declare_event_output<float> ("dRmin_ele_jet");
  h_dRmin_ele_jet_scaled = ctx.declare_event_output<float> ("dRmin_ele_jet_scaled");
  h_ptrel_mu_jet = ctx.declare_event_output<float> ("ptrel_mu_jet");
  h_ptrel_ele_jet = ctx.declare_event_output<float> ("ptrel_ele_jet");

  h_deepjet_topscore_jet = ctx.declare_event_output<float> ("deepjet_topscore_jet");
  h_deepjet_topscore_jet1 = ctx.declare_event_output<float> ("deepjet_topscore_jet1");
  h_deepjet_topscore_jet2 = ctx.declare_event_output<float> ("deepjet_topscore_jet2");
  h_deepjet_wscore_jet  = ctx.declare_event_output<float> ("deepjet_wscore_jet");
  h_deepjet_wscore_jet1 = ctx.declare_event_output<float> ("deepjet_wscore_jet1");
  h_deepjet_wscore_jet2 = ctx.declare_event_output<float> ("deepjet_wscore_jet2");
  h_deepjet_zscore_jet = ctx.declare_event_output<float> ("deepjet_zscore_jet");
  h_deepjet_zscore_jet1 = ctx.declare_event_output<float> ("deepjet_zscore_jet1");
  h_deepjet_zscore_jet2 = ctx.declare_event_output<float> ("deepjet_zscore_jet2");
  h_deepjet_higgsscore_jet = ctx.declare_event_output<float> ("deepjet_higgsscore_jet");
  h_deepjet_higgsscore_jet1  = ctx.declare_event_output<float> ("deepjet_higgsscore_jet1");
  h_deepjet_higgsscore_jet2 = ctx.declare_event_output<float> ("deepjet_higgsscore_jet2");
  h_deepjet_qcdscore_jet = ctx.declare_event_output<float> ("deepjet_qcdscore_jet");
  h_deepjet_qcdscore_jet1 = ctx.declare_event_output<float> ("deepjet_qcdscore_jet1");
  h_deepjet_qcdscore_jet2 = ctx.declare_event_output<float> ("deepjet_qcdscore_jet2");

  // Book histograms
  vector<string> histogram_tags = {"Weights", "Weights_MuonID", "Weights_PU", "Weights_Lumi", "Weights_TopPt", "Weights_MCScale", "Muon1", "TriggerMuon", "Muon2", "Electron1", "TriggerEle", "TwoDCut", "Jet1", "Jet2", "MET", "HTlep", "NNInputsBeforeReweight", "MatchableBeforeChi2Cut", "NotMatchableBeforeChi2Cut", "CorrectMatchBeforeChi2Cut", "NotCorrectMatchBeforeChi2Cut", "Chi2", "Matchable", "NotMatchable", "CorrectMatch", "NotCorrectMatch", "TopTagReconstruction", "NotTopTagReconstruction", "Btags2", "Btags1","TopJetBtagSubjet", "SR_1top", "SR_0top","CR_ttbar","CR_wjet","Test_1top"};
  book_histograms(ctx, histogram_tags);
}

/*
██████  ██████   ██████   ██████ ███████ ███████ ███████
██   ██ ██   ██ ██    ██ ██      ██      ██      ██
██████  ██████  ██    ██ ██      █████   ███████ ███████
██      ██   ██ ██    ██ ██      ██           ██      ██
██      ██   ██  ██████   ██████ ███████ ███████ ███████
*/

bool ZprimeBDTApply::process(uhh2::Event& event){

  if(debug)   cout << "++++++++++++ NEW EVENT ++++++++++++++" << endl;
  if(debug)   cout<<" run.event: "<<event.run<<". "<<event.event<<endl;
  // Initialize reco flags with false
  event.set(h_is_zprime_reconstructed_chi2, false);
  event.set(h_is_zprime_reconstructed_correctmatch, false);
  event.set(h_chi2,-100);
  event.set(h_MET,-100);
  event.set(h_Mttbar,-100);
  event.set(h_lep1_pt,-100);
  event.set(h_lep1_eta,-100);
  event.set(h_ak4jet1_pt,-100);
  event.set(h_ak4jet1_eta,-100);
  event.set(h_ak8jet1_pt,-100);
  event.set(h_ak8jet1_eta,-100);
  event.set(h_NPV,-100);
  event.set(h_weight,-100);


  event.set(h_ttagN,-100);
  //  event.set(h_wtagN,-100);                                                                                                      
  event.set(h_btagN,-100);


  event.set(h_N_jets,0);
  event.set(h_mass_jet,0);
  event.set(h_mass_jet1,0);
  event.set(h_mass_jet2,0);
  event.set(h_ST,0);
  event.set(h_csv_jet,0);
  event.set(h_csv_jet1,0);
  event.set(h_csv_jet2,0);
  event.set(h_S33,0);
  event.set(h_dR_mu_jet,0);
  event.set(h_dRmin_mu_jet,0);
  event.set(h_dRmin_mu_jet_scaled,0);
  event.set(h_dR_ele_jet,0);
  event.set(h_dRmin_ele_jet,0);
  event.set(h_dRmin_ele_jet_scaled,0);
  event.set(h_ptrel_mu_jet,0);
  event.set(h_ptrel_ele_jet,0);

  event.set(h_deepjet_topscore_jet,0);
  event.set(h_deepjet_topscore_jet1,0);
  event.set(h_deepjet_topscore_jet2,0);
  event.set(h_deepjet_wscore_jet,0);
  event.set(h_deepjet_wscore_jet1,0);
  event.set(h_deepjet_wscore_jet2,0);
  event.set(h_deepjet_zscore_jet,0);
  event.set(h_deepjet_zscore_jet1,0);
  event.set(h_deepjet_zscore_jet2,0);
  event.set(h_deepjet_higgsscore_jet,0);
  event.set(h_deepjet_higgsscore_jet1,0);
  event.set(h_deepjet_higgsscore_jet2,0);
  event.set(h_deepjet_qcdscore_jet,0);
  event.set(h_deepjet_qcdscore_jet1,0);
  event.set(h_deepjet_qcdscore_jet2,0);


  // Printing
  // if(!event.isRealData) printer_genparticles->process(event);

  // TODO Apply things that should've been done in the pre-selection already... Fix pre-selection and then remove these steps
  if(isMuon) muon_cleaner->process(event);
  if(isElectron) electron_cleaner->process(event);
  if(debug)  cout<<"Muon and Electron cleaner ok"<<endl;


  if(!HEM_selection->passes(event)){
    if(!isMC) return false;
    else event.weight = event.weight*(1-0.64774715284); // calculated following instructions ar https://twiki.cern.ch/twiki/bin/view/CMS/PdmV2018Analysis
  }

  // Run top-tagging
  TopTaggerPuppi->process(event);
  if(debug) cout<<"Top Tagger ok"<<endl;
  
  fill_histograms(event, "Weights");

  // Weight modules
  if(isMuon){
    MuonID_module->process(event);
    if(debug)  cout<<"MuonID ok"<<endl;
  }
  fill_histograms(event, "Weights_MuonID");
 // if(isElectron){
 //   EleID_module->process(event);
 //   if(debug)  cout<<"EleID ok"<<endl;
 // }
  PUWeight_module->process(event);
  if(debug)  cout<<"PUWeight ok"<<endl;
  fill_histograms(event, "Weights_PU");

  LumiWeight_module->process(event);
  if(debug)  cout<<"LumiWeight ok"<<endl;
  fill_histograms(event, "Weights_Lumi");

  TopPtReweight_module->process(event);
  fill_histograms(event, "Weights_TopPt");

  MCScale_module->process(event);
  fill_histograms(event, "Weights_MCScale");
 
  //BTagWeight_module->process(event);
 
  if(!(Trigger1_selection->passes(event)|| Trigger2_selection->passes(event))) return false;
  if(isMuon){
    if(!NMuon1_selection->passes(event)) return false;
    fill_histograms(event, "Muon1");
    MuonTrigger_module->process_onemuon(event, 0);
    fill_histograms(event, "TriggerMuon");
    if(!NMuon2_selection->passes(event)) return false;
    fill_histograms(event, "Muon2");
  }
  if(isElectron){
    if(!NElectron_selection->passes(event)) return false;
    fill_histograms(event, "Electron1");
    //EleTrigger_module->process(event);
    //fill_histograms(event, "TriggerEle");
  }
  if((event.muons->size()+event.electrons->size()) != 1) return false; //veto events without leptons or with too many 
  if(debug) cout<<"N leptons ok: Nelectrons="<<event.electrons->size()<<" Nmuons="<<event.muons->size()<<endl;
  if(!TwoDCut_selection->passes(event)) return false;
  fill_histograms(event, "TwoDCut");

  // Here, the Zprime must be reconstructed (we ensured to have >= 2 AK4 jets, >= 1 muon)
  // Only consider well-separated AK4 jets
  // //FixME: 
  // for(const TopJet & chsjet : *event.topjets){
  //   int n_separatedjets = 0;
  //   for(unsigned int k = 0; k < event.jets->size(); k++){
  //     if(deltaR(event.jets->at(k), chsjet) > 1.2) n_separatedjets++;
  //   }
  //   //    cout<<"n_separatedjets = "<<n_separatedjets<<endl;
  //   if(n_separatedjets<1) return false;
  // }

  CandidateBuilder->process(event);
  if(debug) cout<<"CandidateBuilder is ok"<<endl;
  Chi2DiscriminatorZprime->process(event);
  if(debug)  cout<<"Chi2DiscriminatorZprime is ok"<<endl;
  CorrectMatchDiscriminatorZprime->process(event);
  if(debug) cout<<"CorrectMatchDiscriminatorZprime is ok"<<endl;
  if(sample.Contains("_blinded")){
    if(!BlindData_selection->passes(event)) return false;
  }
  if(!Jet1_selection->passes(event)) return false;
  if(debug) cout<<"Jet1_selection is ok"<<endl;
  fill_histograms(event, "Jet1");

  if(!Jet2_selection->passes(event)) return false;
  if(debug) cout<<"Jet2_selection is ok"<<endl;
  fill_histograms(event, "Jet2");

  // MET selection
  if(!met_sel->passes(event)) return false;
  if(debug) cout<<"MET is ok"<<endl;
  fill_histograms(event, "MET");
  if(isMuon){
    if(!htlep_sel->passes(event)) return false;
    fill_histograms(event, "HTlep");
    if(debug) cout<<"HTlep is ok"<<endl;
  }

  //  Variables for NN 
  //  Variables_module->process(event);
  //  fill_histograms(event, "NNInputsBeforeReweight");

  if(TTbarMatchable_selection->passes(event)) fill_histograms(event, "MatchableBeforeChi2Cut");
  else fill_histograms(event, "NotMatchableBeforeChi2Cut");
  if(debug) cout<<"TTbarMatchable_selection is ok"<<endl;

  if(Chi2CandidateMatched_selection->passes(event)) fill_histograms(event, "CorrectMatchBeforeChi2Cut");
  else fill_histograms(event, "NotCorrectMatchBeforeChi2Cut");
  if(debug) cout<<"Chi2CandidateMatched_selection is ok"<<endl;

  if(!Chi2_selection->passes(event)) return false;
  fill_histograms(event, "Chi2");

  if(debug) cout<<"Chi2_selection is ok"<<endl;

  if(TTbarMatchable_selection->passes(event)) fill_histograms(event, "Matchable");
  else fill_histograms(event, "NotMatchable");
  if(debug) cout<<"TTbarMatchable_selection is ok"<<endl;

  if(Chi2CandidateMatched_selection->passes(event)) fill_histograms(event, "CorrectMatch");
  else fill_histograms(event, "NotCorrectMatch");
  if(debug) cout<<"Chi2CandidateMatched_selection is ok"<<endl;

  if(ZprimeTopTag_selection->passes(event)) fill_histograms(event, "TopTagReconstruction");
  else fill_histograms(event, "NotTopTagReconstruction");
  if(debug) cout<<"ZprimeTopTag_selection is ok"<<endl;

  /*
  if(event.bdt > 0.4){
    if() fill_histograms(event, "SR_02");
    else fill_histograms(event, "SR_01");
  }
  else{
    if 
      }    
  */
  //Test with b-tagging
  if(sel_2btag->passes(event)) fill_histograms(event, "Btags2");
  if(sel_1btag->passes(event)) fill_histograms(event, "Btags1");
  if(debug) cout<<"Btags1 is ok"<<endl;
  if(TopJetBtagSubjet_selection->passes(event)) fill_histograms(event, "TopJetBtagSubjet");
  if(debug) cout<<"TopJetBtagSubjet_selection is ok"<<endl;

  //Fill some vars for monitoring and comparison
  bool is_zprime_reconstructed_chi2 = event.get(h_is_zprime_reconstructed_chi2); 
  if(is_zprime_reconstructed_chi2){
    ZprimeCandidate* BestZprimeCandidate = event.get(h_BestZprimeCandidateChi2);
    //    float chi2 = BestZprimeCandidate->discriminator("chi2_total");
    event.set(h_chi2,BestZprimeCandidate->discriminator("chi2_total"));
    event.set(h_Mttbar,BestZprimeCandidate->Zprime_v4().M());
  }
  if(debug) cout<<"Set ttbar reconstruction vars for monitoring"<<endl;

  event.set(h_weight,event.weight);
  event.set(h_MET,event.met->pt());
  if(isMuon){
    event.set(h_lep1_pt,event.muons->at(0).pt());
    event.set(h_lep1_eta,event.muons->at(0).eta());
  }
  if(isElectron){
    event.set(h_lep1_pt,event.electrons->at(0).pt());
    event.set(h_lep1_eta,event.electrons->at(0).eta());
  }
  if(event.jets->size()>0){
    event.set(h_ak4jet1_pt,event.jets->at(0).pt());
    event.set(h_ak4jet1_eta,event.jets->at(0).eta());
  }
  if(event.topjets->size()>0){
    event.set(h_ak8jet1_pt,event.topjets->at(0).pt());
    event.set(h_ak8jet1_eta,event.topjets->at(0).eta());
  }
  event.set(h_NPV,event.pvs->size());
  if(debug) cout<<"Set some vars for monitoring"<<endl;

  //limit variables                                                                                                                 
  //  vector<TopJet> WTags = event.get(h_AK8WTags);                                                                                 
  //  event.set(h_wtagN,WTags.size());                                                                                              

  //  vector<TopJet> TopTags = event.toppuppijets;                                                                                  
  event.set(h_ttagN,event.toppuppijets->size());

  int btagN = 0;
  int jetbtagN = 0;
  int btag_subjet=0;

  if(event.toppuppijets->size() == 0){
    if(sel_2btag->passes(event)) jetbtagN = jetbtagN + 2;
    if(sel_1btag->passes(event)) jetbtagN = jetbtagN + 1;
  }

  if(event.toppuppijets->size() == 1){
    for(auto & topjet : *event.topjets){
      auto subjets = topjet.subjets();
      for (auto & subjet : subjets) {
        if(sel_1btag->passes(event)) btag_subjet++;
      }


      for(unsigned int j=0; j<event.jets->size(); j++){
        double dR_jet = deltaR(topjet, event.jets->at(j));
        if(dR_jet > 1.2) continue;
        if(sel_1btag->passes(event)) jetbtagN++;
      }
    }
  }

  btagN = jetbtagN + btag_subjet;
  event.set(h_btagN,btagN);


  //BDT
  event.set(h_N_jets,event.jets->size());
  for(unsigned int i=0; i<event.jets->size(); i++){
    event.set(h_mass_jet,event.jets->at(i).v4().M());
    event.set(h_csv_jet,event.jets->at(i).btag_combinedSecondaryVertex());
    if(isMuon){
      double dRmin_muon_jet = 99999;
      for(unsigned int j=0; j<event.muons->size(); j++){
        double dR_mu_jet = deltaR(event.jets->at(i), event.muons->at(j));
        if(dR_mu_jet < dRmin_muon_jet) dRmin_muon_jet = dR_mu_jet;
        event.set(h_dR_mu_jet,dR_mu_jet);
      }
    }
    if(isElectron){
      double dRmin_ele_jet = 99999;
      for(unsigned int j=0; j<event.electrons->size(); j++){
        double dR_ele_jet = deltaR(event.jets->at(i), event.electrons->at(j));
        if(dR_ele_jet < dRmin_ele_jet) dRmin_ele_jet = dR_ele_jet;
        event.set(h_dR_ele_jet,dR_ele_jet);
      }
    }
    if(i==0){
      event.set(h_mass_jet1,event.jets->at(i).v4().M());
      event.set(h_csv_jet1,event.jets->at(i).btag_combinedSecondaryVertex());
    }

    else if(i==1){
      event.set(h_mass_jet2,event.jets->at(i).v4().M());
      event.set(h_csv_jet2,event.jets->at(i).btag_combinedSecondaryVertex());
    }

  }

  double st = 0., st_jets = 0., st_lep = 0.;


  for(unsigned int i=0; i<event.jets->size(); i++){
    st_jets += event.jets->at(i).pt();
  }
  for(unsigned int i=0; i<event.electrons->size(); i++){
    st_lep += event.electrons->at(i).pt();
  }
  for(unsigned int i=0; i<event.muons->size(); i++){
    st_lep += event.muons->at(i).pt();
  }
  st = st_jets + st_lep + event.met->pt();

  event.set(h_ST,st);

  double s33 = -1., mag = -1.;
  for(const Jet jet : *event.jets){
    mag += (jet.v4().Px()*jet.v4().Px()+jet.v4().Py()*jet.v4().Py()+jet.v4().Pz()*jet.v4().Pz());
    s33 += jet.v4().Pz()*jet.v4().Pz();
  }
  s33 = s33 / mag;
  event.set(h_S33,s33);

  if(isMuon){

    if(event.muons->at(0).has_tag(Muon::twodcut_dRmin) && event.muons->at(0).has_tag(Muon::twodcut_pTrel)){
      event.set(h_dRmin_mu_jet,event.muons->at(0).get_tag(Muon::twodcut_dRmin));
      event.set(h_ptrel_mu_jet,event.muons->at(0).get_tag(Muon::twodcut_pTrel));
      event.set(h_dRmin_mu_jet_scaled, event.muons->at(0).get_tag(Muon::twodcut_dRmin)*event.jets->at(0).pt());
    }
  }

  if(isElectron){
    if(event.electrons->at(0).has_tag(Electron::twodcut_dRmin) && event.electrons->at(0).has_tag(Electron::twodcut_pTrel)){
      event.set(h_dRmin_ele_jet,event.electrons->at(0).get_tag(Electron::twodcut_dRmin));
      event.set(h_ptrel_ele_jet,event.electrons->at(0).get_tag(Electron::twodcut_pTrel));
      event.set(h_dRmin_ele_jet_scaled, event.electrons->at(0).get_tag(Electron::twodcut_dRmin)*event.jets->at(0).pt());
    }
  }


  vector<TopJet>* AK8Puppijets = event.toppuppijets;
  unsigned int NAK8Puppijets = 0;
  for(unsigned int i=0; i<AK8Puppijets->size(); i++){
    if(AK8Puppijets->at(i).numberOfDaughters()<2) continue;
    NAK8Puppijets++;

    event.set(h_deepjet_topscore_jet, AK8Puppijets->at(i).btag_DeepBoosted_raw_score_top());
    event.set(h_deepjet_wscore_jet, AK8Puppijets->at(i).btag_DeepBoosted_raw_score_w());
    event.set(h_deepjet_zscore_jet, AK8Puppijets->at(i).btag_DeepBoosted_raw_score_z());
    event.set(h_deepjet_higgsscore_jet, AK8Puppijets->at(i).btag_DeepBoosted_raw_score_h());
    event.set(h_deepjet_qcdscore_jet, AK8Puppijets->at(i).btag_DeepBoosted_raw_score_qcd());


    if(i==0){
      event.set(h_deepjet_topscore_jet1, AK8Puppijets->at(i).btag_DeepBoosted_raw_score_top());
      event.set(h_deepjet_wscore_jet1, AK8Puppijets->at(i).btag_DeepBoosted_raw_score_w());
      event.set(h_deepjet_zscore_jet1, AK8Puppijets->at(i).btag_DeepBoosted_raw_score_z());
      event.set(h_deepjet_higgsscore_jet1, AK8Puppijets->at(i).btag_DeepBoosted_raw_score_h());
      event.set(h_deepjet_qcdscore_jet1, AK8Puppijets->at(i).btag_DeepBoosted_raw_score_qcd());
    }
    if(i==1){
      event.set(h_deepjet_topscore_jet2, AK8Puppijets->at(i).btag_DeepBoosted_raw_score_top());
      event.set(h_deepjet_wscore_jet2, AK8Puppijets->at(i).btag_DeepBoosted_raw_score_w());
      event.set(h_deepjet_zscore_jet2, AK8Puppijets->at(i).btag_DeepBoosted_raw_score_z());
      event.set(h_deepjet_higgsscore_jet2, AK8Puppijets->at(i).btag_DeepBoosted_raw_score_h());
      event.set(h_deepjet_qcdscore_jet2, AK8Puppijets->at(i).btag_DeepBoosted_raw_score_qcd());
    }

  }


  if(event.bdt_score > 0.75){
    if(ZprimeTopTag_selection->passes(event)) fill_histograms(event, "SR_1top");
    else fill_histograms(event, "SR_0top");
    //    if (event.topjets->size() == 1) fill_histograms(event, "SR_1top");
    //    else if (event.topjets->size() == 0) fill_histograms(event, "SR_0top");
  }
  else {
    if(event.bdt_score < 0.15) fill_histograms(event, "CR_wjet");
    else if (event.bdt_score > 0.6) fill_histograms(event, "CR_ttbar");
  }

  if (event.topjets->size() == 1) fill_histograms(event, "Test_1top");
  return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(ZprimeBDTApply)
