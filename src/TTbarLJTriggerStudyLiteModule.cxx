#include <iostream>
#include <memory>

#include <UHH2/core/include/AnalysisModule.h>
#include <UHH2/core/include/Event.h>
#include <UHH2/core/include/Selection.h>
#include <UHH2/common/include/CleaningModules.h>

#include <UHH2/common/include/MCWeight.h>
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
#include <UHH2/common/include/TTbarReconstruction.h>
#include <UHH2/common/include/ReconstructionHypothesisDiscriminators.h>
#include <UHH2/common/include/Utils.h>
#include <UHH2/common/include/JetHists.h>
#include <UHH2/common/include/TopPtReweight.h>
#include "UHH2/common/include/LuminosityHists.h"

#include <UHH2/ZprimeSemiLeptonic/include/ModuleBASE.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicSelections.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicUtils.h>
#include <UHH2/ZprimeSemiLeptonic/include/TTbarLJTriggerEffHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/DileptonHists.h>
#include <UHH2/ZprimeSemiLeptonic/include/EffyTTbarRECOHists.h>

#include <UHH2/ZprimeSemiLeptonic/include/SF_elec.h>
#include <UHH2/ZprimeSemiLeptonic/include/SF_ttagging.h>
#include <UHH2/ZprimeSemiLeptonic/include/SF_WjetsREWGT.h>

#include <TMVA/Tools.h>
#include <TMVA/Reader.h>
#include <TMVA/MethodCuts.h>


using namespace std;

class TTbarLJTriggerStudyLiteModule : public ModuleBASE {

 public:
  explicit TTbarLJTriggerStudyLiteModule(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;

 protected:
  enum lepton { muon, elec };
  lepton channel_;

  // selections
  std::unique_ptr<uhh2::Selection> lumi_sel;

  std::unique_ptr<uhh2::Selection> jet2_sel;
  std::unique_ptr<uhh2::Selection> jet1_sel;
  std::unique_ptr<uhh2::Selection> trigger_sel;
  std::unique_ptr<uhh2::Selection> trigger2_sel;
  std::unique_ptr<uhh2::Selection> trigger3_sel;
  std::unique_ptr<uhh2::Selection> trigger4_sel;
  std::unique_ptr<uhh2::Selection> trigger5_sel;
  std::unique_ptr<uhh2::Selection> trigger6_sel;
  std::unique_ptr<uhh2::Selection> tag_trigger_sel;
  std::unique_ptr<uhh2::Selection> tag_trigger2_sel;
  std::unique_ptr<uhh2::Selection> tag_trigger3_sel;
  std::unique_ptr<uhh2::Selection> met_sel;
  std::unique_ptr<uhh2::Selection> htlep_sel;
  //  std::unique_ptr<uhh2::Selection> triangc_sel;
  std::unique_ptr<uhh2::Selection> twodcut_sel;
  // cleaners
  std::unique_ptr<MuonCleaner>     muoSR_cleaner;
  std::unique_ptr<ElectronCleaner> eleSR_cleaner;
  //  ElectronId eleID = ElectronID_Spring15_25ns_tight_noIso;
  ElectronId eleID =   ElectronID_Fall17_tight_noIso; //check OLD version
  //  ElectronId eleID = ElectronID_Spring16_medium_noIso;         // MVA version of Zprime
  //  ElectronId eleID = ElectronID_MVAGeneralPurpose_Spring16_loose; //Cut-based version of Zprime
  std::unique_ptr<JetCleaner>                      jet_IDcleaner;
  //  std::unique_ptr<JetCorrector>                    jet_corrector, jet_corrector_BCD, jet_corrector_EFearly, jet_corrector_FlateG, jet_corrector_H;
  std::unique_ptr<JetCorrector>                    jet_corrector_A,jet_corrector_B, jet_corrector_C, jet_corrector_D, jet_corrector_E, jet_corrector_F, jet_corrector_G,jet_corrector_H, jet_corrector_MC;
  std::unique_ptr<GenericJetResolutionSmearer>     jetER_smearer;
  //  std::unique_ptr<JetLeptonCleaner> jetlepton_cleaner, JLC_BCD, JLC_EFearly, JLC_FlateG, JLC_H;
  std::unique_ptr<JetLeptonCleaner_by_KEYmatching> JLC_A, JLC_B, JLC_C, JLC_D,JLC_E, JLC_F, JLC_G, JLC_H, JLC_MC;
  std::unique_ptr<JetCleaner>                      jet_cleaner1;
  std::unique_ptr<JetCleaner>                      jet_cleaner2;

  JetId btag_ID_;
  CSVBTag::wp b_working_point;

  bool use_ttagging_;
  TopJetId ttag_ID_;
  float    ttag_minDR_jet_;

  bool blind_DATA_;
  bool isTrigger_sel;
  bool store_PDF_weights_;

  //// Data/MC scale factors
  std::unique_ptr<uhh2::AnalysisModule> pileupSF;
  std::unique_ptr<uhh2::AnalysisModule> muonID_SF;
  std::unique_ptr<uhh2::AnalysisModule> muonHLT_SF;
  //  std::unique_ptr<uhh2::AnalysisModule> muonTRK_SF;
  std::unique_ptr<uhh2::AnalysisModule> elecID_SF;
  std::unique_ptr<uhh2::AnalysisModule> elecGsf_SF;

  //// VARS

  // event
  Event::Handle<int>   h_run;
  Event::Handle<int>   h_lumi_block;
  Event::Handle<int>   h_event;
  Event::Handle<int>   h_pvN;
  Event::Handle<int>   h_muoN;
  Event::Handle<int>   h_eleN;
 
  Event::Handle<float> tt_TMVA_response;// response of TMVA method
  
  ////

  //Runnumbers for applying different corrections
  constexpr static int s_runnr_B_2016  = 275376; //up to this one, including this one
  constexpr static int s_runnr_C_2016  = 276283; //up to this one, including this one
  constexpr static int s_runnr_D_2016 =  276811; //up to this one, including this one
  constexpr static int s_runnr_E_2016 =  277420; //up to this one, including this one
  constexpr static int s_runnr_F_2016 =  278801; //up to this one, including this one = Fearly
  constexpr static int s_runnr_G_2016 =  280385; //up to this one, including this one
  constexpr static int s_runnr_H_2016 =  284044; //up to this one, including this one

  constexpr static int s_runnr_B_2017  = 299329; //up to this one, including this one
  constexpr static int s_runnr_C_2017  = 302029; //up to this one, including this one
  constexpr static int s_runnr_D_2017 =  303434; //up to this one, including this one
  constexpr static int s_runnr_E_2017 =  304826; //up to this one, including this one
  constexpr static int s_runnr_F_2017  = 306462; //up to this one, including this one

  constexpr static int s_runnr_A_2018  = 316995; //up to this one, including this one
  constexpr static int s_runnr_B_2018  = 319310; //up to this one, including this one
  constexpr static int s_runnr_C_2018  = 320065; //up to this one, including this one
  constexpr static int s_runnr_D_2018 =  325175; //up to this one, including this one



  std::unique_ptr<Hists> lumihists_nosel; std::unique_ptr<Hists> lumihists_dilep;
  std::unique_ptr<Hists> lumihists_twodcut; std::unique_ptr<Hists> lumihists_jet2;
  std::unique_ptr<Hists> lumihists_jet1; std::unique_ptr<Hists> lumihists_met;
  std::unique_ptr<Hists> lumihists_tag; std::unique_ptr<Hists> lumihists_probe;


  //new for the 2017, 2018
  bool is2016v2, is2016v3, is2017, is2018;
  bool isMC, ispuppi;
};

TTbarLJTriggerStudyLiteModule::TTbarLJTriggerStudyLiteModule(uhh2::Context& ctx){

  //// CONFIGURATION
  const bool isMC = (ctx.get("dataset_type") == "MC");

  const std::string& channel = ctx.get("channel", "");
  is2016v2 = (ctx.get("dataset_version").find("2016v2") != std::string::npos);
  is2016v3 = (ctx.get("dataset_version").find("2016v3") != std::string::npos);
  is2017 = (ctx.get("dataset_version").find("2017") != std::string::npos);
  is2018 = (ctx.get("dataset_version").find("2018") != std::string::npos);



  use_ttagging_ = false;

  blind_DATA_ = ((ctx.get("dataset_version").find("BLINDED") != std::string::npos) && (ctx.get("dataset_type") == "DATA") && !isMC);
  isTrigger_sel = ((ctx.get("dataset_version").find("HLT") != std::string::npos) && (ctx.get("dataset_type") == "DATA") );
  const std::string& store_PDF_weights = ctx.get("store_PDF_weights", "");
  if     (store_PDF_weights == "true")  store_PDF_weights_ = true;
  else if(store_PDF_weights == "false") store_PDF_weights_ = false;
  else {

    std::string log("TTbarLJTriggerStudyLiteModule::TTbarLJTriggerStudyLiteModule -- ");
    log += "invalid argument for 'store_PDF_weights' key in xml file (must be 'true' or 'false'): \""+store_PDF_weights+"\"";

    throw std::runtime_error(log);
  }
  //

  const std::string& keyword = ctx.get("keyword");

  float jet1_pt(-1.), jet2_pt(-1.), MET(-1.), HT_lep(-1.);
  float muon1_pt (-1.),ele1_pt(-1.);

  use_ttagging_ = false;
  muon1_pt = 55.;
  //  ele1_pt = 50.; 
  //jet1_pt  = 185.;
    ele1_pt = 50.; 
    jet1_pt  = 185.;
  //  jet1_pt  = 160.;
  //  jet1_pt  = 200.;

  jet2_pt  = 50.; 
  //  jet2_pt  = 30.; 
  //  MET      =  50.;
  MET      =  120.;
  HT_lep   = 150.;  //used the value from the 2016 analysis note, might need to change it later


  //// COMMON MODULES

  /* CMS-certified luminosity sections */
  if(!isMC) lumi_sel.reset(new LumiSelection(ctx));

  /**********************************/

  ////

  //// OBJ CLEANING
  //  const     MuonId muoSR(AndId<Muon>    (PtEtaCut  (muon1_pt, 2.5), MuonIDTight()));
  //  const     MuonId muoSR((PtEtaCut(muon1_pt, 2.5), MuonID(Muon::CutBasedIdTight)));
  const MuonId muoSR(AndId<Muon> (PtEtaCut(muon1_pt, 2.5), MuonID(Muon::CutBasedIdTight)));
  const ElectronId eleSR(AndId<Electron>(PtEtaSCCut(ele1_pt, 2.4), eleID));
  muoSR_cleaner.reset(new     MuonCleaner(muoSR));
  eleSR_cleaner.reset(new ElectronCleaner(eleSR));
  //  const JetId jetID(JetPFID(JetPFID::WP_TIGHT_PUPPI));
  const JetPFID jetID_CHS(JetPFID::WP_TIGHT_CHS); 
  const JetPFID jetID_PUPPI(JetPFID::WP_TIGHT_PUPPI); 




  //Jet energy correction
  std::vector<std::string> JEC_AK4CHS_MC, JEC_AK8CHS_MC, JEC_AK4Puppi_MC, JEC_AK8Puppi_MC;
  std::vector<std::string> JEC_AK4CHS_A, JEC_AK4CHS_B, JEC_AK4CHS_C, JEC_AK4CHS_D, JEC_AK4CHS_E, JEC_AK4CHS_F, JEC_AK4CHS_G, JEC_AK4CHS_H;
  std::vector<std::string> JEC_AK4Puppi_A, JEC_AK4Puppi_B, JEC_AK4Puppi_C, JEC_AK4Puppi_D, JEC_AK4Puppi_E, JEC_AK4Puppi_F, JEC_AK4Puppi_G, JEC_AK4Puppi_H;
  std::vector<std::string> JEC_AK8CHS_A, JEC_AK8CHS_B, JEC_AK8CHS_C, JEC_AK8CHS_D, JEC_AK8CHS_E, JEC_AK8CHS_F, JEC_AK8CHS_G, JEC_AK8CHS_H;
  std::vector<std::string> JEC_AK8Puppi_A, JEC_AK8Puppi_B, JEC_AK8Puppi_C, JEC_AK8Puppi_D, JEC_AK8Puppi_E, JEC_AK8Puppi_F, JEC_AK8Puppi_G, JEC_AK8Puppi_H;

  if(is2017){
    cout<<"ZprimePreselectionModule uses JEC for 2017 data/MC"<<endl;
    JEC_AK4CHS_B       = JERFiles::Fall17_17Nov2017_V32_B_L123_AK4PFchs_DATA;
    JEC_AK4CHS_C       = JERFiles::Fall17_17Nov2017_V32_C_L123_AK4PFchs_DATA;
    JEC_AK4CHS_D       = JERFiles::Fall17_17Nov2017_V32_D_L123_AK4PFchs_DATA;
    JEC_AK4CHS_E       = JERFiles::Fall17_17Nov2017_V32_E_L123_AK4PFchs_DATA;
    JEC_AK4CHS_F       = JERFiles::Fall17_17Nov2017_V32_F_L123_AK4PFchs_DATA;
    JEC_AK4CHS_MC       = JERFiles::Fall17_17Nov2017_V32_L123_AK4PFchs_MC;

    JEC_AK8CHS_B       = JERFiles::Fall17_17Nov2017_V32_B_L123_AK8PFchs_DATA;
    JEC_AK8CHS_C       = JERFiles::Fall17_17Nov2017_V32_C_L123_AK8PFchs_DATA;
    JEC_AK8CHS_D       = JERFiles::Fall17_17Nov2017_V32_D_L123_AK8PFchs_DATA;
    JEC_AK8CHS_E       = JERFiles::Fall17_17Nov2017_V32_E_L123_AK8PFchs_DATA;
    JEC_AK8CHS_F       = JERFiles::Fall17_17Nov2017_V32_F_L123_AK8PFchs_DATA;
    JEC_AK8CHS_MC       = JERFiles::Fall17_17Nov2017_V32_L123_AK8PFchs_MC;

    JEC_AK4Puppi_B = JERFiles::Fall17_17Nov2017_V32_B_L123_AK4PFPuppi_DATA;
    JEC_AK4Puppi_C = JERFiles::Fall17_17Nov2017_V32_C_L123_AK4PFPuppi_DATA;
    JEC_AK4Puppi_D = JERFiles::Fall17_17Nov2017_V32_D_L123_AK4PFPuppi_DATA;
    JEC_AK4Puppi_E = JERFiles::Fall17_17Nov2017_V32_E_L123_AK4PFPuppi_DATA;
    JEC_AK4Puppi_F = JERFiles::Fall17_17Nov2017_V32_F_L123_AK4PFPuppi_DATA;
    JEC_AK4Puppi_MC = JERFiles::Fall17_17Nov2017_V32_L123_AK4PFPuppi_MC;

    JEC_AK8Puppi_B = JERFiles::Fall17_17Nov2017_V32_B_L123_AK8PFPuppi_DATA;
    JEC_AK8Puppi_C = JERFiles::Fall17_17Nov2017_V32_C_L123_AK8PFPuppi_DATA;
    JEC_AK8Puppi_D = JERFiles::Fall17_17Nov2017_V32_D_L123_AK8PFPuppi_DATA;
    JEC_AK8Puppi_E = JERFiles::Fall17_17Nov2017_V32_E_L123_AK8PFPuppi_DATA;
    JEC_AK8Puppi_F = JERFiles::Fall17_17Nov2017_V32_F_L123_AK8PFPuppi_DATA;
    JEC_AK8Puppi_MC = JERFiles::Fall17_17Nov2017_V32_L123_AK8PFPuppi_MC;
  }
  
  if(is2018){
    cout<<"ZprimePreselectionModule uses JEC for 2018 data/MC"<<endl;
    JEC_AK4CHS_A       = JERFiles::Autumn18_V8_A_L123_AK4PFchs_DATA;
    JEC_AK4CHS_B       = JERFiles::Autumn18_V8_B_L123_AK4PFchs_DATA;
    JEC_AK4CHS_C       = JERFiles::Autumn18_V8_C_L123_AK4PFchs_DATA;
    JEC_AK4CHS_D       = JERFiles::Autumn18_V8_D_L123_AK4PFchs_DATA;
    JEC_AK4CHS_MC       = JERFiles::Autumn18_V8_L123_AK4PFchs_MC;

    JEC_AK8CHS_A       = JERFiles::Autumn18_V8_A_L123_AK8PFchs_DATA;
    JEC_AK8CHS_B       = JERFiles::Autumn18_V8_B_L123_AK8PFchs_DATA;
    JEC_AK8CHS_C       = JERFiles::Autumn18_V8_C_L123_AK8PFchs_DATA;
    JEC_AK8CHS_D       = JERFiles::Autumn18_V8_D_L123_AK8PFchs_DATA;
    JEC_AK8CHS_MC       = JERFiles::Autumn18_V8_L123_AK8PFchs_MC;

    JEC_AK4Puppi_A = JERFiles::Autumn18_V8_A_L123_AK4PFPuppi_DATA;
    JEC_AK4Puppi_B = JERFiles::Autumn18_V8_B_L123_AK4PFPuppi_DATA;
    JEC_AK4Puppi_C = JERFiles::Autumn18_V8_C_L123_AK4PFPuppi_DATA;
    JEC_AK4Puppi_D = JERFiles::Autumn18_V8_D_L123_AK4PFPuppi_DATA;
    JEC_AK4Puppi_MC = JERFiles::Autumn18_V8_L123_AK4PFPuppi_MC;

    JEC_AK8Puppi_A = JERFiles::Autumn18_V8_A_L123_AK8PFPuppi_DATA;
    JEC_AK8Puppi_B = JERFiles::Autumn18_V8_B_L123_AK8PFPuppi_DATA;
    JEC_AK8Puppi_C = JERFiles::Autumn18_V8_C_L123_AK8PFPuppi_DATA;
    JEC_AK8Puppi_D = JERFiles::Autumn18_V8_D_L123_AK8PFPuppi_DATA;
    JEC_AK8Puppi_MC = JERFiles::Autumn18_V8_L123_AK8PFPuppi_MC;
  }


  if(!ispuppi){
    if(is2018) jet_corrector_A.reset(new JetCorrector(ctx, JEC_AK4CHS_A));
    jet_corrector_B.reset(new JetCorrector(ctx, JEC_AK4CHS_B));
    jet_corrector_C.reset(new JetCorrector(ctx, JEC_AK4CHS_C));
    jet_corrector_D.reset(new JetCorrector(ctx, JEC_AK4CHS_D));
    if(is2017 or is2016v2 or is2016v3){
      //      jet_corrector_D.reset(new JetCorrector(ctx, JEC_AK4CHS_D));//tmp
      jet_corrector_E.reset(new JetCorrector(ctx, JEC_AK4CHS_E));
      jet_corrector_F.reset(new JetCorrector(ctx, JEC_AK4CHS_F));
      if(is2016v2 or is2016v3){
	jet_corrector_G.reset(new JetCorrector(ctx, JEC_AK4CHS_G));
	jet_corrector_H.reset(new JetCorrector(ctx, JEC_AK4CHS_H));
      }
    }
    jet_corrector_MC.reset(new JetCorrector(ctx,            JEC_AK4CHS_MC));

    if(is2018) JLC_A.reset(new JetLeptonCleaner_by_KEYmatching(ctx, JEC_AK4CHS_A));
    JLC_B.reset(new JetLeptonCleaner_by_KEYmatching(ctx, JEC_AK4CHS_B));
    JLC_C.reset(new JetLeptonCleaner_by_KEYmatching(ctx, JEC_AK4CHS_C));
    JLC_D.reset(new JetLeptonCleaner_by_KEYmatching(ctx, JEC_AK4CHS_D));
    if(is2017 or is2016v2 or is2016v3){
      //  JLC_D.reset(new JetLeptonCleaner_by_KEYmatching(ctx, JEC_AK4CHS_D));//tmp
      JLC_E.reset(new JetLeptonCleaner_by_KEYmatching(ctx, JEC_AK4CHS_E));
      JLC_F.reset(new JetLeptonCleaner_by_KEYmatching(ctx, JEC_AK4CHS_F));
      if(is2016v2 or is2016v3){
	JLC_G.reset(new JetLeptonCleaner_by_KEYmatching(ctx, JEC_AK4CHS_G));
	JLC_H.reset(new JetLeptonCleaner_by_KEYmatching(ctx, JEC_AK4CHS_H));
      }
    }
    JLC_MC.reset(new JetLeptonCleaner_by_KEYmatching(ctx,            JEC_AK4CHS_MC));
  }
  
  else{
    if(is2018) jet_corrector_A.reset(new JetCorrector(ctx, JEC_AK4Puppi_A));
    jet_corrector_B.reset(new JetCorrector(ctx, JEC_AK4Puppi_B));
    jet_corrector_C.reset(new JetCorrector(ctx, JEC_AK4Puppi_C));
    jet_corrector_D.reset(new JetCorrector(ctx, JEC_AK4Puppi_D));
    if(is2017 or is2016v2 or is2016v3){
      //  jet_corrector_D.reset(new JetCorrector(ctx, JEC_AK4Puppi_D));//tmp
      jet_corrector_E.reset(new JetCorrector(ctx, JEC_AK4Puppi_E));
      jet_corrector_F.reset(new JetCorrector(ctx, JEC_AK4Puppi_F));
      if(is2016v2 or is2016v3){
	jet_corrector_G.reset(new JetCorrector(ctx, JEC_AK4Puppi_G));
	jet_corrector_H.reset(new JetCorrector(ctx, JEC_AK4Puppi_H));
      }
    }
    jet_corrector_MC.reset(new JetCorrector(ctx, JEC_AK4Puppi_MC));

    //Jet-lepton cleaning
    if(is2018) JLC_A.reset(new JetLeptonCleaner_by_KEYmatching(ctx, JEC_AK4Puppi_A));
    JLC_B.reset(new JetLeptonCleaner_by_KEYmatching(ctx, JEC_AK4Puppi_B));
    JLC_C.reset(new JetLeptonCleaner_by_KEYmatching(ctx, JEC_AK4Puppi_C));
    JLC_D.reset(new JetLeptonCleaner_by_KEYmatching(ctx, JEC_AK4Puppi_D));
    if(is2017 or is2016v2 or is2016v3){
      //  JLC_D.reset(new JetLeptonCleaner_by_KEYmatching(ctx, JEC_AK4Puppi_D));//tmp
      JLC_E.reset(new JetLeptonCleaner_by_KEYmatching(ctx, JEC_AK4Puppi_E));
      JLC_F.reset(new JetLeptonCleaner_by_KEYmatching(ctx, JEC_AK4Puppi_F));
      if(is2016v2 or is2016v3){
	JLC_G.reset(new JetLeptonCleaner_by_KEYmatching(ctx, JEC_AK4Puppi_G));
	JLC_H.reset(new JetLeptonCleaner_by_KEYmatching(ctx, JEC_AK4Puppi_H));
      }
    }
    JLC_MC.reset(new JetLeptonCleaner_by_KEYmatching(ctx,            JEC_AK4Puppi_MC));
  }




  //  jetlepton_cleaner.reset(new JetLeptonCleaner(ctx, JEC_AK4));
  jet_IDcleaner.reset(new JetCleaner(ctx, jetID_PUPPI));

 
  jet_cleaner1.reset(new JetCleaner(ctx, 15., 3.0));
  jet_cleaner2.reset(new JetCleaner(ctx, 30., 2.4));

  //// EVENT SELECTION
  jet2_sel.reset(new NJetSelection(2, -1, JetId(PtEtaCut(jet2_pt, 2.4))));
  jet1_sel.reset(new NJetSelection(1, -1, JetId(PtEtaCut(jet1_pt, 2.4))));

  const std::string& trigger = ctx.get("probe_trigger", "NULL");
  if(trigger != "NULL") trigger_sel.reset(new TriggerSelection(trigger));
  else                  trigger_sel.reset(new uhh2::AndSelection(ctx));
  const std::string& trigger2 = ctx.get("probe_trigger2", "NULL");
  if(trigger2 != "HLT_NULL_v*") trigger2_sel.reset(new TriggerSelection(trigger2));
  else                  trigger2_sel.reset(new uhh2::AndSelection(ctx));
  const std::string& trigger3 = ctx.get("probe_trigger3", "NULL");
  if(trigger3 != "HLT_NULL_v*") trigger3_sel.reset(new TriggerSelection(trigger3));
  else                  trigger3_sel.reset(new uhh2::AndSelection(ctx));
  const std::string& trigger4 = ctx.get("probe_trigger4", "NULL");
  if(trigger4 != "HLT_NULL_v*") trigger4_sel.reset(new TriggerSelection(trigger4));
  else                  trigger4_sel.reset(new uhh2::AndSelection(ctx));
  const std::string& trigger5 = ctx.get("probe_trigger5", "NULL");
  if(trigger5 != "HLT_NULL_v*") trigger5_sel.reset(new TriggerSelection(trigger5));
  else                  trigger5_sel.reset(new uhh2::AndSelection(ctx));
  const std::string& trigger6 = ctx.get("probe_trigger6", "NULL");
  if(trigger6 != "HLT_NULL_v*") trigger6_sel.reset(new TriggerSelection(trigger6));
  else                  trigger6_sel.reset(new uhh2::AndSelection(ctx));

  const std::string& tag_trigger = ctx.get("tag_trigger", "NULL");
  if(tag_trigger != "NULL") tag_trigger_sel.reset(new TriggerSelection(tag_trigger));
  else  tag_trigger_sel.reset(new uhh2::AndSelection(ctx));
  const std::string& tag_trigger2 = ctx.get("tag_trigger2", "NULL");
  if(tag_trigger2 != "NULL") tag_trigger2_sel.reset(new TriggerSelection(tag_trigger2));
  else  tag_trigger2_sel.reset(new uhh2::AndSelection(ctx));
  const std::string& tag_trigger3 = ctx.get("tag_trigger3", "NULL");
  if(tag_trigger3 != "NULL") tag_trigger3_sel.reset(new TriggerSelection(tag_trigger3));
  else  tag_trigger3_sel.reset(new uhh2::AndSelection(ctx));
  std::cout<<"tag_trigger = "<<tag_trigger<<" tag_trigger2 = "<<tag_trigger2<<std::endl;
  std::cout<<"probe_trigger = "<<trigger<<" probe_trigger2 = "<<trigger2<<" probe_trigger3 = "<<trigger3<<" probe_trigger4 = "<<trigger4<<" probe_trigger5 = "<<trigger5<<" probe_trigger6 = "<<trigger6<<std::endl;
  met_sel  .reset(new METCut  (MET   , uhh2::infinity));
  htlep_sel.reset(new HTlepCut(HT_lep, uhh2::infinity));
  twodcut_sel.reset(new TwoDCutALL(0.4, 20.)); //may need to change
  // twodcut_sel.reset(new TwoDCut1(0.4, 20.));

  //// HISTS
  std::vector<std::string> htags_1({
      "nosel",
      "dilep",
      "twodcut",
      "jet2",
      "jet1",
       "met",
      "tag",
      "tagNprobe",
  });

  for(const auto& tag : htags_1){

    book_HFolder(tag, new TTbarLJTriggerEffHists(ctx, tag, ttag_ID_, ttag_minDR_jet_));
  }
  //

  
  //Luminosity
  
  lumihists_nosel.reset(new LuminosityHists(ctx, "lumi_nosel"));
  lumihists_dilep.reset(new LuminosityHists(ctx, "lumi_dilep"));
  lumihists_twodcut.reset(new LuminosityHists(ctx, "lumi_twodcut"));
  lumihists_jet2.reset(new LuminosityHists(ctx, "lumi_jet2"));
  lumihists_jet1.reset(new LuminosityHists(ctx, "lumi_jet1"));
  lumihists_met.reset(new LuminosityHists(ctx, "lumi_met"));
  lumihists_tag.reset(new LuminosityHists(ctx, "lumi_tag"));
  lumihists_probe.reset(new LuminosityHists(ctx, "lumi_tagNprobe"));

  ////

  //// Data/MC scale factors

  //

  // muon-ID
  const std::string& muonID_SFac    = ctx.get("muonID_SF_file");
  const std::string& muonID_directory    = ctx.get("muonID_SF_directory");

  // muon-HLT
  const std::string& muonHLT_SFac   = ctx.get("muonHLT_SF_file");
  const std::string& muonHLT_directory   = ctx.get("muonHLT_SF_directory");

  // // muon-TRK
  //  const std::string& muonTRK_SFac   = ctx.get("muonTRK_SF_file");


  // elec-ID
  const std::string& elecID_SFac    = ctx.get("elecID_SF_file");
  const std::string& elecGsf_SFac    = ctx.get("elecGsf_SF_file");

  ////

  //// VARS
  ctx.undeclare_all_event_output();

  //pileup (define it after undeclaring all other variables to keep the weights in the output)
  pileupSF.reset(new MCPileupReweight(ctx));

  //muon scale factors
  muonID_SF.reset(new MCMuonScaleFactor(ctx, muonID_SFac, muonID_directory, 1.0, "ID"));
    muonHLT_SF.reset(new MCMuonScaleFactor(ctx, muonHLT_SFac, muonHLT_directory, 0.5, "HLT"));
  //  muonTRK_SF.reset(new MCMuonTrkScaleFactor(ctx, muonTRK_SFac, 0.0, "TRK"));

  //electron scale factors
  elecID_SF.reset(new MCElecScaleFactor(ctx, elecID_SFac, 0.0, "ID"));
  elecGsf_SF.reset(new MCElecScaleFactor(ctx, elecGsf_SFac, 0.0, "Gsf"));

  // //b-tagging scale factors
  // btagSF.reset(new MCBTagScaleFactor(ctx, b_working_point));

  // event
  h_run             = ctx.declare_event_output<int>("run");
  h_lumi_block      = ctx.declare_event_output<int>("lumi_block");
  h_event           = ctx.declare_event_output<int>("event");

  h_pvN             = ctx.declare_event_output<int>("pvN");

  h_muoN            = ctx.declare_event_output<int>("muoN");
  h_eleN            = ctx.declare_event_output<int>("eleN");

  //Dummy TMVA response
  tt_TMVA_response            = ctx.declare_event_output<float>("TMVA_response");

}

bool TTbarLJTriggerStudyLiteModule::process(uhh2::Event& event){
  // std::cout<<" --- New event! --- "<<std::endl;
  //  std::cout<<"Jet N = "<<event.jets->size()<<std::endl;
  int jet_n = event.jets->size();
  if(jet_n<2) return false;

  HFolder("nosel")->fill(event);
  if(event.isRealData) lumihists_nosel->fill(event);
  //double w_before = event.weight;
  //std::cout<<" event.weight [BEFORE] = "<<event.weight<<std::endl;
  event.set(tt_TMVA_response, 0);
  //// COMMON MODULES


  // if(!event.isRealData){

  //   ttgenprod->process(event);
  // }

  /* CMS-certified luminosity sections */
  if(event.isRealData){

    if(!lumi_sel->passes(event)) return false;
  }
  //  std::cout<<" Went through Lumi_sel "<<std::endl;
  ////

  //// Data/MC scale factors
  


  //pileup
  pileupSF->process(event);

  // //  std::cout<<" event.weight [after pileup] = "<<event.weight<<std::endl;
  // // // b-tagging
  // // btagSF->process(event);
  //  double w_after_pileup = event.weight;
  
  //// LEPTON SELECTION
  muoSR_cleaner->process(event);
  sort_by_pt<Muon>(*event.muons);
  eleSR_cleaner->process(event);
  sort_by_pt<Electron>(*event.electrons);


  //// LEPTON SFs
  //  double w_before = event.weight;
  //  std::cout<<" event.weight [BEFORE] = "<<event.weight<<std::endl;
  // muon SFs
  muonID_SF->process(event);
  muonHLT_SF->process(event);
  //  muonTRK_SF->process(event);

  // elec SFs
  elecID_SF->process(event);
  elecGsf_SF->process(event);
  // //  double w_after = event.weight;
  // //std::cout<<" event.weight [AFTER lepton SFs] = "<<event.weight<<std::endl;
  //  double w_after_lepSFs = event.weight;
  

  //di-lepton selection
  bool pass_dilep(0);
  pass_dilep = (event.muons->size()==1) && (event.electrons->size()==1) && (event.muons->at(0).charge()!=event.electrons->at(0).charge());  //electron-muon
  if(!pass_dilep) return false;
  //  if(!(event.muons->at(0).pt() > 55. && (abs(event.muons->at(0).eta()) <2.5))) return false;
  //  if(!(event.electrons->at(0).pt() > 50. && (abs(event.electrons->at(0).eta()) <2.4))) return false;

  HFolder("dilep")->fill(event);
  if(event.isRealData) lumihists_dilep->fill(event);


  if(event.electrons->at(0).Class()==4) return false;//TEST veto gap electrons  

  // // //// JET selection
  jet_IDcleaner->process(event);
  //  std::cout<<"Jet N = "<<event.jets->size()<<std::endl;
  //  std::cout<<" Went through jet_IDcleaner "<<std::endl;

  if(event.isRealData){
    bool apply_A = false;
    bool apply_B = false;
    bool apply_C = false;
    bool apply_D = false;
    bool apply_E = false;
    bool apply_F = false;
    bool apply_G = false;
    bool apply_H = false;
    if(is2016v2 || is2016v3){
      if(event.run <= s_runnr_B_2016)  apply_B = true;
      else if(event.run <= s_runnr_C_2016) apply_C = true;
      else if(event.run <= s_runnr_D_2016) apply_D = true;
      else if(event.run <= s_runnr_E_2016) apply_E = true;
      else if(event.run <= s_runnr_F_2016) apply_F = true;
      else if(event.run <= s_runnr_G_2016) apply_G = true;
      else if(event.run <= s_runnr_H_2016) apply_H = true;
      else throw std::runtime_error("run number not covered by if-statements in process-routine.");
    }

    if(is2017){
      if(event.run <= s_runnr_B_2017)  apply_B = true;
      else if(event.run <= s_runnr_C_2017) apply_C = true;
      else if(event.run <= s_runnr_D_2017) apply_D = true;
      else if(event.run <= s_runnr_E_2017) apply_E = true;
      else if(event.run <= s_runnr_F_2017) apply_F = true;
      else throw std::runtime_error("run number not covered by if-statements in process-routine.");
    }
    if(is2018){
      if(event.run <= s_runnr_A_2018)  apply_A = true;
      else if(event.run <= s_runnr_B_2018)  apply_B = true;
      else if(event.run <= s_runnr_C_2018) apply_C = true;
      else if(event.run <= s_runnr_D_2018) apply_D = true;
      else throw std::runtime_error("run number not covered by if-statements in process-routine.");
    }
    


    if(apply_A+apply_B+apply_C+apply_D+apply_E+apply_F+apply_G+apply_H != 1) throw std::runtime_error("In ZprimePreselectionModule.cxx: Sum of apply_* when applying JECs is not == 1. Fix this.");




  //apply proper JECs


  if(apply_A){
    JLC_A->process(event);
    jet_corrector_A->process(event);
    jet_corrector_A->correct_met(event);
  }
  if(apply_B){
    JLC_B->process(event);
    jet_corrector_B->process(event);
    jet_corrector_B->correct_met(event);
  }
  if(apply_C){
    JLC_C->process(event);
    jet_corrector_C->process(event);
    jet_corrector_C->correct_met(event);
  }
  if(apply_D){
    JLC_D->process(event);
    jet_corrector_D->process(event);
    jet_corrector_D->correct_met(event);
  }

  if(apply_E){
    JLC_E->process(event);
    jet_corrector_E->process(event);
    jet_corrector_E->correct_met(event);

  }
  if(apply_F){
    JLC_F->process(event);
    jet_corrector_F->process(event);
    jet_corrector_F->correct_met(event);
  }    
  if(apply_G){
    JLC_G->process(event);
    jet_corrector_G->process(event);
    jet_corrector_G->correct_met(event);
  }
  if(apply_H){
    JLC_H->process(event);
    jet_corrector_H->process(event);
    jet_corrector_H->correct_met(event);
  }
 }
 else{ //MC
   //    cout<<"JLC, JEC, JER ..."<<endl;
   JLC_MC->process(event);
   jet_corrector_MC->process(event);
   jet_corrector_MC->correct_met(event);

 }


// std::cout<<" Went through JEC "<<std::endl;
  jet_cleaner1->process(event);
  sort_by_pt<Jet>(*event.jets);



  //  const int 
  jet_n = event.jets->size();
  //  std::cout<<"Jet N = "<<event.jets->size()<<std::endl;
  if(jet_n<2) return false;

  /* lepton-2Dcut variables */



  const bool pass_twodcut = twodcut_sel->passes(event); {

    for(auto& muo : *event.muons){

      float    dRmin, pTrel;
      std::tie(dRmin, pTrel) = drmin_pTrel(muo, *event.jets);

      muo.set_tag(Muon::twodcut_dRmin, dRmin);
      muo.set_tag(Muon::twodcut_pTrel, pTrel);
    }

    for(auto& ele : *event.electrons){

      float    dRmin, pTrel;
      std::tie(dRmin, pTrel) = drmin_pTrel(ele, *event.jets);

      ele.set_tag(Electron::twodcut_dRmin, dRmin);
      ele.set_tag(Electron::twodcut_pTrel, pTrel);
    }
  }

  // //// LEPTON-2Dcut selection
   if(!pass_twodcut) return false;
   HFolder("twodcut")->fill(event);
   if(event.isRealData) lumihists_twodcut->fill(event);
   jet_cleaner2->process(event);                                                                                                   
  sort_by_pt<Jet>(*event.jets);

  // /* 2nd AK4 jet selection */
  const bool pass_jet2 = jet2_sel->passes(event);

  if(!pass_jet2) return false;
  HFolder("jet2")->fill(event);
  if(event.isRealData) lumihists_jet2->fill(event);

  // /* 1st AK4 jet selection */
  const bool pass_jet1 = jet1_sel->passes(event);
  if(!pass_jet1) return false;
  HFolder("jet1")->fill(event);
  if(event.isRealData) lumihists_jet1->fill(event);
  //// MET selection
  const bool pass_met = met_sel->passes(event);
  if(!pass_met) return false;
  HFolder("met")->fill(event);
  if(event.isRealData) lumihists_met->fill(event);
  //tag trigger
  bool pass_tag_trigger=false; bool pass_tag_trigger2=false; bool pass_tag_trigger3=false;
  pass_tag_trigger = tag_trigger_sel->passes(event);
  pass_tag_trigger2 = tag_trigger2_sel->passes(event);
  pass_tag_trigger3 = tag_trigger3_sel->passes(event);
  //  if(event.run>=274953) pass_tag_trigger2 = tag_trigger2_sel->passes(event);
  //  else pass_tag_trigger2=false;//TkMu50 becomes available after 274953

  if(!pass_tag_trigger && !pass_tag_trigger2 && !pass_tag_trigger3) return false;
  //if(!pass_tag_trigger) return false; //TEST
  HFolder("tag")->fill(event);
  //  double w_after = event.weight;
  //  std::cout<<" event.weight = "<<event.weight<<std::endl;
  // std::cout<<" event.weight [AFTER Pile-up]/[BEFORE] = "<<100.*(w_after_pileup/w_before)<<std::endl;
  // std::cout<<" event.weight [AFTER lep SFs]/[BEFORE] = "<<100.*(w_after_lepSFs/w_before)<<std::endl;
  // std::cout<<" event.weight [AFTER]/[BEFORE] = "<<100.*(w_after/w_before)<<std::endl;
  if(event.isRealData) lumihists_tag->fill(event);


  //probe trugger
  bool pass_trigger = false;  
  bool pass_trigger2 = false;
  //  bool pass_trigger3 = false;  
  //bool pass_trigger4 = false;
  // bool pass_trigger5 = false;  bool pass_trigger6 = false;

  pass_trigger = trigger_sel->passes(event);
  pass_trigger2 = trigger2_sel->passes(event);
  //  pass_trigger3 = trigger3_sel->passes(event);
  // pass_trigger4 = trigger4_sel->passes(event);
  // pass_trigger5 = trigger5_sel->passes(event);
  // pass_trigger6 = trigger6_sel->passes(event);
  // if(pass_trigger2 || pass_trigger3 || pass_trigger4 || pass_trigger5 || pass_trigger6)
  //   std::cout<<"pass_trigger2 || pass_trigger3 || pass_trigger4 || pass_trigger5 || pass_trigger6 "<<std::endl;
  //  if(pass_trigger6)  std::cout<<"pass_trigger6 = "<<pass_trigger6<<std::endl;
  //  if(!pass_trigger && !pass_trigger2 && !pass_trigger3 && !pass_trigger4 && !pass_trigger5 && !pass_trigger6) return false;
  //  if(!pass_trigger && !pass_trigger2 && !pass_trigger3 && !pass_trigger6) return false;//TEST
  //  if(!pass_trigger && !pass_trigger2 && !pass_trigger3) return false;//TEST
  //  std::cout<<"pass_tag_trigger = "<<pass_tag_trigger<<" pass_trigger = "<<pass_trigger<<std::endl;
  //  if(!pass_trigger) return false;//TEST with one trigger
  if(!pass_trigger && !pass_trigger2) return false;
  //  if(!pass_trigger && !pass_trigger2 && !pass_trigger3) return false;
  HFolder("tagNprobe")->fill(event);
  if(event.isRealData) lumihists_probe->fill(event);

  // // /*****************/

  //// VARS

  // event
  event.set(h_run       , event.run);
  event.set(h_lumi_block, event.luminosityBlock);
  event.set(h_event     , event.event);

  event.set(h_pvN  , event.pvs->size());

  event.set(h_muoN , event.muons    ->size());
  event.set(h_eleN , event.electrons->size());
  //  std::cout<<" --- END --- "<<std::endl;
  return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(TTbarLJTriggerStudyLiteModule)
