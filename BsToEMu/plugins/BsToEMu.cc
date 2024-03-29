#include <memory>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h" // muy importante para MiniAOD (in English: very important for MiniAOD)

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "RecoEgamma/EgammaElectronAlgos/interface/GsfElectronAlgo.h"
#include "DataFormats/ParticleFlowReco/interface/GsfPFRecTrack.h"
#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateTransform.h"
#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateMode.h"
//#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
//#include "DataFormats/EgammaCandidates/interface/GsfElectronCore.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
//#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexSorter.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoVertex/PrimaryVertexProducer/interface/VertexHigherPtSquared.h"

#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TLorentzVector.h"
#include "KVFitter.h"
#include "TFile.h"
#include "TTree.h"
//using reco::TrackCollection;
//
#include "Util.h"

class BsToEMu : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit BsToEMu(const edm::ParameterSet&);
  ~BsToEMu();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  KVFitter vtxfit;
  Utilities  utilityFuntions;
  // ----------member data ---------------------------
  //edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
  //edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
  //edm::EDGetTokenT<pat::MuonCollection> muonToken_;
  edm::EDGetTokenT<edm::View<pat::Electron>> electronToken_;
  edm::EDGetTokenT<edm::View<pat::Muon>> muonToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
  edm::EDGetTokenT<reco::VertexCollection> vertices_;
  edm::EDGetTokenT<reco::VertexCollection> ltrk_;
  edm::EDGetTokenT<edm::TriggerResults> trg_;
  //edm::EDGetTokenT<edm::View<pat::TriggerObjectStandAlone>> tobjs_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> tobjs_;

  bool isMC_, OnlyGen_;
  TTree *tree;
  unsigned int nB;

  int run,lumi,event, npv,ncands, hlt_path; //std::vector<int> *trgIndex;
  std::vector<float> *t_pt, *t_eta, *t_phi;
  std::vector<float> *mass, *mcorr, *pt, *eta, *phi;
  std::vector<int> *charge;
  std::vector<float> *dr_mu, *dr_ele;
  std::vector<float> *pv_x, *pv_y, *pv_z; float bs_x0, bs_y0, bs_z0; std::vector<float> *bs_x, *bs_y, *vx, *vy, *vz;
  std::vector<float> *vtx_chi2, *vtx_prob, *cos2d, *lxy, *lxy_err, *lxy_sig;
  std::vector<float> *mu_pt, *mu_eta, *mu_phi, *mu_e, *mu_mass;
  std::vector<int> *mu_charge;
  std::vector<float> *mu_id_loose, *mu_id_soft, *mu_id_medium, *mu_id_tight, *mu_id_soft_mva, *mu_dxy, *mu_dxy_e, *mu_dxy_sig, *mu_dz, *mu_dz_e, *mu_dz_sig, *mu_bs_dxy, *mu_bs_dxy_e, *mu_bs_dxy_sig, *mu_cov_pos_def, *mu_trkIsolation;
  std::vector<float> *ele_pt, *ele_eta, *ele_phi, *ele_e, *ele_mass;
  std::vector<int> *ele_charge;
  std::vector<float> *ele_id_iso_wp90, *ele_id_iso_wp80, *ele_id_iso_wpLoose, *ele_id_iso, *ele_id_noIso, *ele_id_noIso_wp90, *ele_id_noIso_wp80, *ele_id_noIso_wpLoose, *ele_dxy, *ele_dxy_e, *ele_dxy_sig, *ele_dz, *ele_dz_e, *ele_dz_sig, *ele_bs_dxy, *ele_bs_dxy_e, *ele_bs_dxy_sig, *ele_cov_pos_def, *ele_trkIsolation;
  std::vector<float> *mu_trk_chi2, *mu_trk_ndof, *mu_trk_prob;
  std::vector<float> *ele_sigmaietaieta, *ele_trk_chi2, *ele_trk_ndof, *ele_trk_prob;
  std::vector<float> *charge_pt;
  TH1F *h1_eff00 = new TH1F("h1_eff00", "hlt is passed", 10, 0, 10);
  TH1F *h1_eff01 = new TH1F("h1_eff01", "muon pt and eta", 10, 0, 10);
  TH1F *h1_eff02 = new TH1F("h1_eff02", "muon best track", 10, 0, 10);
  TH1F *h1_eff03 = new TH1F("h1_eff03", "pf muon", 10, 0, 10);
  TH1F *h1_eff04 = new TH1F("h1_eff04", "global muon", 10, 0, 10);
  TH1F *h1_eff05 = new TH1F("h1_eff05", "loose muon", 10, 0, 10);
  TH1F *h1_eff06 = new TH1F("h1_eff06", "trg_dr", 10, 0, 10);
  TH1F *h1_eff07 = new TH1F("h1_eff07", "trg_pt", 10, 0, 10);
  TH1F *h1_eff08 = new TH1F("h1_eff08", "trg matched", 10, 0, 10);
  TH1F *h1_eff09 = new TH1F("h1_eff09", "gsfTrack", 10, 0, 10);
  TH1F *h1_eff10 = new TH1F("h1_eff10", "ele-muon charge", 10, 0, 10);
  TH1F *h1_eff11 = new TH1F("h1_eff11", "ele pt and eta", 10, 0, 10);
  TH1F *h1_eff12 = new TH1F("h1_eff12", "ele iso wpLoose", 10, 0, 10);
  TH1F *h1_eff13 = new TH1F("h1_eff13", "vertex valid", 10, 0, 10);
  TH1F *h1_eff14 = new TH1F("h1_eff14", "cos2d", 10, 0, 10);
  TH1F *h1_eff15 = new TH1F("h1_eff15", "mass of B-cand", 10, 0, 10);

  TH1F *h1_mu_pt_1 = new TH1F("h1_mu_pt_1", "p_{T} (#mu) before", 100, 0, 50);
  TH1F *h1_mu_pt_2 = new TH1F("h1_mu_pt_2", "p_{T} (#mu) after", 100, 0, 50);
  TH1F *h1_mu_pt_3 = new TH1F("h1_mu_pt_3", "p_{T} (#mu) before Loose/PF", 100, 0, 50);
  TH1F *h1_mu_pt_4 = new TH1F("h1_mu_pt_4", "p_{T} (#mu) after Loose/PF", 100, 0, 50);
  TH1F *h1_mu_eta_1 = new TH1F("h1_mu_eta_1", "#eta (#mu) before", 200, -6, 6);
  TH1F *h1_mu_eta_2 = new TH1F("h1_mu_eta_2", "#eta (#mu) after", 100, -3, 3);
  TH1F *h1_mu_eta_3 = new TH1F("h1_mu_eta_3", "#eta (#mu) before Loose/PF", 100, -3, 3);
  TH1F *h1_mu_eta_4 = new TH1F("h1_mu_eta_4", "#eta (#mu) after Loose/PF", 100, -3, 3);
  TH1F *h1_delPt_b1 = new TH1F("h1_delPt_b1", "#Deltap_T(#mu, trg) before", 100, 0, 6);
  TH1F *h1_delPt_b2 = new TH1F("h1_delPt_b2", "#Deltap_T(#mu, trg) before", 100, 0, 0.5);
  TH1F *h1_delPt_a = new TH1F("h1_delPt_a", "#Deltap_T(#mu, trg) after", 100, 0, 0.5);
  TH1F *h1_delRtrg_b1 = new TH1F("h1_delRtrg_b1", "#DeltaR(#mu, trg) before", 100, 0, 20);
  TH1F *h1_delRtrg_b2 = new TH1F("h1_delRtrg_b2", "#DeltaR(#mu, trg) before", 100, 0, 0.5);
  TH1F *h1_delRtrg_a = new TH1F("h1_delRtrg_a", "#DeltaR(#mu, trg) after", 100, 0, 0.5);
  TH1F *h1_el_pt_1 = new TH1F("h1_el_pt_1", "p_{T} (e) before charge", 100, 0, 50);
  TH1F *h1_el_pt_2 = new TH1F("h1_el_pt_2", "p_{T} (e) after charge", 100, 0, 50);
  TH1F *h1_el_pt_3 = new TH1F("h1_el_pt_3", "p_{T} (e) before", 100, 0, 50);
  TH1F *h1_el_pt_4 = new TH1F("h1_el_pt_4", "p_{T} (e) after", 100, 0, 50);
  TH1F *h1_el_eta_1 = new TH1F("h1_el_eta_1", "#eta (e) before", 200, -3, 3);
  TH1F *h1_el_eta_2 = new TH1F("h1_el_eta_2", "#eta (e) after", 100, -3, 3);
  TH1F *h1_el_eta_3 = new TH1F("h1_el_eta_3", "#eta (e) before", 200, -3, 3);
  TH1F *h1_el_eta_4 = new TH1F("h1_el_eta_4", "#eta (e) after", 100, -3, 3);
  TH1F *h1_el_idIso_1 = new TH1F("h1_el_idIso_1", "electron id iso before", 100, -1, 1.01);
  TH1F *h1_el_idNoIso_1 = new TH1F("h1_el_idNoIso_1", "electron id noIso before", 100, -1, 1.01);
  TH1F *h1_el_idIso_2 = new TH1F("h1_el_idIso_2", "electron id iso after", 100, -1, 1.01);
  TH1F *h1_el_idNoIso_2 = new TH1F("h1_el_idNoIso_2", "electron id noIso after", 100, -1, 1.01);
  TH1F *h1_cos2d_1 = new TH1F("h1_cos2d_1", "cos#alpha_{2D} before", 100, 0.9, 1.001 );
  TH1F *h1_cos2d_2 = new TH1F("h1_cos2d_2", "cos#alpha_{2D} after", 100, 0.9, 1.001 );
  TH1F *h1_mass_1 = new TH1F("h1_mass_1", "mass of e#mu before", 100, 0, 10);
  TH1F *h1_mass_2 = new TH1F("h1_mass_2", "mass of e#mu after", 100, 0, 10);
};

BsToEMu::BsToEMu(const edm::ParameterSet& iConfig):
  //tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))){
  //electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
  //muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  electronToken_(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electrons"))),
  muonToken_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
  beamSpotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
  pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
  vertices_ (consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  ltrk_ (consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("ltrk"))),
  trg_ (consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("trg_res"))),
//tobjs_ (consumes<edm::View<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("tobjs"))),
  tobjs_ (consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("tobjs"))),
  
  isMC_(iConfig.getParameter<bool>("isMC")),
  OnlyGen_(iConfig.getParameter<bool>("OnlyGen")),
  //now do what ever initialization is needed
  tree(0),
  nB(0),
  run(0), lumi(0), event(0), npv(0), ncands(0), hlt_path(0), //trgIndex(0),
   t_pt(0), t_eta(0), t_phi(0),
  mass(0), mcorr(0), pt(0), eta(0), phi(0), charge(0), dr_mu(0), dr_ele(0), pv_x(0), pv_y(0), pv_z(0), bs_x0(0), bs_y0(0), bs_z0(0), bs_x(0), bs_y(0), vx(0), vy(0), vz(0), vtx_chi2(0), vtx_prob(0), cos2d(0), lxy(0), lxy_err(0), lxy_sig(0),
  mu_pt(0), mu_eta(0), mu_phi(0), mu_e(0), mu_mass(0), mu_charge(0), mu_id_loose(0), mu_id_soft(0), mu_id_medium(0), mu_id_tight(0), mu_id_soft_mva(0), mu_dxy(0), mu_dxy_e(0), mu_dxy_sig(0), mu_dz(0), mu_dz_e(0), mu_dz_sig(0), mu_bs_dxy(0), mu_bs_dxy_e(0), mu_bs_dxy_sig(0), mu_cov_pos_def(0), mu_trkIsolation(0), ele_pt(0), ele_eta(0), ele_phi(0), ele_e(0), ele_mass(0), ele_charge(0), ele_id_iso_wp90(0), ele_id_iso_wp80(0), ele_id_iso_wpLoose(0), ele_id_iso(0), ele_id_noIso(0), ele_id_noIso_wp90(0), ele_id_noIso_wp80(0), ele_id_noIso_wpLoose(0), ele_dxy(0), ele_dxy_e(0), ele_dxy_sig(0), ele_dz(0), ele_dz_e(0), ele_dz_sig(0), ele_bs_dxy(0), ele_bs_dxy_e(0), ele_bs_dxy_sig(0), ele_cov_pos_def(0), ele_trkIsolation(0),
  mu_trk_chi2(0), mu_trk_ndof(0), mu_trk_prob(0),
  ele_sigmaietaieta(0), ele_trk_chi2(0), ele_trk_ndof(0), ele_trk_prob(0),
  charge_pt(0)
{}

BsToEMu::~BsToEMu(){}

void BsToEMu::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  using std::vector;
  using namespace edm;
  using namespace reco;
  using namespace std;

  //for(const auto& track : iEvent.get(tracksToken_) ) {   }
  //edm::Handle<pat::MuonCollection> muons;
  //iEvent.getByToken(muonToken_, muons);
  //edm::Handle<pat::ElectronCollection> electrons;
  //iEvent.getByToken(electronToken_, electrons);
  edm::Handle<pat::PackedCandidateCollection> pfs;
  iEvent.getByToken(pfToken_, pfs);
  
  //std::vector<const reco::Candidate *> leptons;
  //for (const pat::Muon &mu : *muons) leptons.push_back(&mu);
  //for (const pat::Electron &el : *electrons) leptons.push_back(&el);
  
  edm::Handle<VertexCollection> primaryVertices;
  iEvent.getByToken(vertices_, primaryVertices);

  edm::Handle<VertexCollection> ltrk;
  iEvent.getByToken(ltrk_, ltrk);

  edm::Handle<edm::TriggerResults> trg_res;
  iEvent.getByToken(trg_, trg_res);
  
  //edm::Handle<pat::TriggerObjectStandAlone> tobjs;
  edm::Handle<pat::TriggerObjectStandAloneCollection> tobjs;
  iEvent.getByToken(tobjs_, tobjs);

  lumi = iEvent.id().luminosityBlock();
  run = iEvent.id().run();
  event = iEvent.id().event();
  npv = primaryVertices->size();//iEvent.id().size();

  edm::Handle< View<pat::Muon> > muons;
  iEvent.getByToken(muonToken_, muons );
  
  edm::Handle< View<pat::Electron> > electrons;
  iEvent.getByToken(electronToken_, electrons );
  
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  iEvent.getByToken(beamSpotToken_, recoBeamSpotHandle );
  //beamSpot = beamSpotHandle.product();

  bs_x0 = recoBeamSpotHandle -> x0();
  bs_y0 = recoBeamSpotHandle -> y0();
  bs_z0 = recoBeamSpotHandle -> z0();

  //cout << "something ...   " << recoBeamSpotHandle -> x0() << '\t' << recoBeamSpotHandle -> y0() << '\t' << recoBeamSpotHandle -> z0() <<endl;
  //bs_x  = beamSpot->x0();
  //bs_y  = beamSpot->x0();

  //edm::ESHandle<TransientTrackBuilder> theB; 
  //iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB); 

  //cout << "================================ Something started ======================================\n";
  unsigned int nn[20] = {0};
  int test_event=event;

  static int nhlt = 10;
  const char *paths[nhlt] = { "HLT_Mu7_IP4"     ,
			    "HLT_Mu8_IP3"     ,
			    "HLT_Mu8_IP5"     ,
			    "HLT_Mu8_IP6"     ,
			    "HLT_Mu8p5_IP3p5" ,
			    "HLT_Mu9_IP4"     ,
			    "HLT_Mu9_IP5"     ,
			    "HLT_Mu9_IP6"     ,
			    "HLT_Mu10p5_IP3p5",
			    "HLT_Mu12_IP6"
			    };

  //const edm::TriggerNames &names = iEvent.triggerNames(*trg_res);  int hlt_passed = 0;
  //for( unsigned int i=0, n = trg_res->size(); i < n; ++i) { if( names.triggerName(i) == "HLT_Mu7_IP4_part0_v2" ){ hlt_passed=1; break; }} //cout << names.triggerName(i) << endl; }
  const char *var1;
  auto trg_names = iEvent.triggerNames(*trg_res);
  int hlt_passed=0;
  for(int ipath=0; ipath < nhlt; ++ipath){
    for(unsigned iname = 0; iname < trg_res->size(); ++iname){ //.triggerNames()
      std::string name = trg_names.triggerName(iname);
      var1 = name.c_str();
      if(strstr(var1, paths[ipath]) ){ // && (strlen(var1) - strlen(paths[ipath]) < 5) ){
	if( trg_res -> accept(iname) && strstr(var1, "HLT_Mu7_IP4")) hlt_passed=1;
	if( trg_res -> accept(iname) && strstr(var1, "HLT_Mu8_IP3")) hlt_passed=2;
	if( trg_res -> accept(iname) && strstr(var1, "HLT_Mu8_IP5")) hlt_passed=3;
	if( trg_res -> accept(iname) && strstr(var1, "HLT_Mu8_IP6")) hlt_passed=4;
	if( trg_res -> accept(iname) && strstr(var1, "HLT_Mu8p5_IP3p5")) hlt_passed=5;
	if( trg_res -> accept(iname) && strstr(var1, "HLT_Mu9_IP4")) hlt_passed=6;
	if( trg_res -> accept(iname) && strstr(var1, "HLT_Mu9_IP5")) hlt_passed=7;
	if( trg_res -> accept(iname) && strstr(var1, "HLT_Mu9_IP6")) hlt_passed=8;
	if( trg_res -> accept(iname) && strstr(var1, "HLT_Mu10p5_IP3p5")) hlt_passed=9;
	if( trg_res -> accept(iname) && strstr(var1, "HLT_Mu12_IP6")) hlt_passed=10;
      }
    }
  }
  
  // Checking trigger condition, step one for efficiency calculation
  if( !hlt_passed ) return; 
  nn[0]++;
  hlt_path=hlt_passed;
  //cout << tobjs -> pt() << "  <--- Testing trigger mattching object size.\n"; // working but not correctly
  pat::TriggerObjectStandAlone good_tobjs;
  for(pat::TriggerObjectStandAlone obj : *tobjs) {
    //cout << "************************************** yaha to chal sakata a hai na *******************************************************\n";
    obj.unpackFilterLabels(iEvent, *trg_res);
    if(obj.hasFilterLabel("hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered7IP4Q") ||
       obj.hasFilterLabel("hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered8Q") ||
       obj.hasFilterLabel("hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered8IP5Q") ||
       obj.hasFilterLabel("hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered8IP6Q") ||
       obj.hasFilterLabel("hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered8p5Q") ||
       obj.hasFilterLabel("hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered9IP4Q") ||
       obj.hasFilterLabel("hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered9IP5Q") ||
       obj.hasFilterLabel("hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered9Q") ||
       obj.hasFilterLabel("hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered10p5Q") ||
       obj.hasFilterLabel("hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered12Q") ){
      //cout << "(((((((((((((((((((((((((((((((((((((((((((( if ke andar bhi aa to gya ))))))))))))))))))))))))))))))))))))))))))))\n";
      t_pt ->push_back(obj.pt());
      t_eta->push_back(obj.eta());
      t_phi->push_back(obj.phi());
      //cout <<")))))))))))))))))))))))))))))))))))))))))))) if ke bahar bhi nikal gya (((((((((((((((((((((((((((((((((((((((((((( \n";
    }
  }

  //cout <<  " ================================= Trigger objects are fine ================================================\n";
  std::vector <reco::Track> tracksTifit;
  ParticleMass muon_mass = 0.10565837; //pdg mass
  ParticleMass el_mass = 0.0005109989461; //pdg mass

  for(View<pat::Muon>::const_iterator iMu = muons->begin(); iMu != muons->end(); ++iMu){

    // muon selection 
    h1_mu_pt_1 -> Fill(iMu  -> pt());
    h1_mu_eta_1 -> Fill(iMu -> eta());
    if( iMu  -> pt() < 6.0 || fabs(iMu -> eta()) > 2.5 ) continue;
    nn[1]++;

    if( ! iMu -> bestTrack()) continue;
    nn[2]++;

    h1_mu_pt_2 -> Fill(iMu  -> pt());
    h1_mu_eta_2 -> Fill(iMu -> eta());

    // Checking muon whether it is triggerd, step two for efficiency calculation
    TLorentzVector  mu4, trg4, dif;
    double dpT = 999, t_dr = 999;
    mu4.SetXYZM(iMu->px(), iMu->py(), iMu->pz(), muon_mass);
    for(int trg_i = 0; trg_i < (int)t_eta -> size(); trg_i++){
      double deta = iMu -> eta() - t_eta -> at(0);
      double dphi = fabs(iMu -> phi() - t_phi->at(0));
      if(dphi > M_PI) dphi =  2*M_PI - dphi;
      double t_dr1 = sqrt(deta*deta+dphi*dphi);
      if( t_dr > t_dr1 ) t_dr = t_dr1;
      if( t_dr <= 0.02 ) nn[6]++;
      
      trg4.SetPtEtaPhiM(t_pt->at(trg_i), t_eta->at(trg_i), t_phi->at(trg_i), muon_mass);
      dif = mu4 - trg4;
      double dpT1 = dif.Pt()/iMu->pt();
      if( dpT > dpT1 ) dpT = dpT1;
      if( dpT <= 0.05 ) nn[7]++;
    }
    h1_delPt_b1 -> Fill(dpT);
    h1_delPt_b2 -> Fill(dpT);
    if( dpT > 0.05 ) continue;
    h1_delPt_a -> Fill(dpT);
    h1_delRtrg_b1 -> Fill(t_dr);
    h1_delRtrg_b2 -> Fill(t_dr);
    if( t_dr > 0.02 ) continue;
    h1_delRtrg_a -> Fill(t_dr);
    nn[8]++;

    //cout << "============================= Muons are matching with trigger =========================================\n";

    h1_mu_pt_3 -> Fill(iMu  -> pt());
    h1_mu_eta_3 -> Fill(iMu -> eta());

    if( ! iMu -> isPFMuon() ) continue;
    nn[3]++;
    //if( ! iMu -> isGlobalMuon() ) continue;
    //nn[4]++;
    //if( ! iMu -> isLooseMuon() ) continue;
    if( ! iMu -> isMediumMuon() ) continue;
    nn[5]++;
    h1_mu_pt_4 -> Fill(iMu  -> pt());
    h1_mu_eta_4 -> Fill(iMu -> eta());
    
    for(View<pat::Electron>::const_iterator iE = electrons->begin(); iE != electrons->end(); ++iE){
      
      // electron selection
      if( ! iE -> gsfTrack().isNonnull() ) continue;
      nn[9]++;
      h1_el_pt_1 -> Fill(iE -> pt());
      h1_el_eta_1 -> Fill(iE -> eta());
      //if( (iMu -> charge()) * (iE->charge()) == 1 ) continue;
      //nn[10]++;
      
      h1_el_pt_2 -> Fill(iE -> pt());
      h1_el_eta_2 -> Fill(iE -> eta());
      if( iE -> pt() < 4.0 || fabs(iE -> eta()) > 2.5) continue;
      nn[11]++;

      h1_el_pt_3 -> Fill(iE -> pt());
      h1_el_eta_3 -> Fill(iE -> eta());

      h1_el_idIso_1 -> Fill(iE->userFloat("ElectronMVAEstimatorRun2Fall17IsoV2Values"));
      h1_el_idNoIso_1 -> Fill(iE->userFloat("ElectronMVAEstimatorRun2Fall17NoIsoV2Values"));
      if( ! iE -> electronID("mvaEleID-Fall17-iso-V2-wpLoose") ) continue;
      nn[12]++;
      h1_el_pt_4 -> Fill(iE -> pt());
      h1_el_eta_4 -> Fill(iE -> eta());
      h1_el_idIso_2 -> Fill(iE->userFloat("ElectronMVAEstimatorRun2Fall17IsoV2Values"));
      h1_el_idNoIso_2 -> Fill(iE->userFloat("ElectronMVAEstimatorRun2Fall17NoIsoV2Values"));
      //cout << "=========================== Electrons are good ===========================================\n";


      // Checking vertex conditions, step three for efficiency calculation
      tracksTifit.clear();
      tracksTifit.push_back( *(iMu->bestTrack()) );
      tracksTifit.push_back( *(iE->gsfTrack() ) );
      auto emuVtxFit = vtxfit.Fit( tracksTifit );
      if( ! emuVtxFit.isValid()) continue;
      nn[13]++;

      auto pv = (*primaryVertices)[0];
      int pvIndex=0;
      double dz_mu0 = fabs(iMu -> bestTrack() -> dz(pv.position()));
      for(int ivtx = 1; ivtx < npv; ++ivtx ){
	double dz_mu1 = fabs(iMu -> bestTrack() -> dz((*primaryVertices)[ivtx].position()));
	if( dz_mu0 > dz_mu1 ){
	  dz_mu0 = dz_mu1;
	  pv = (*primaryVertices)[ivtx];
	  pvIndex = ivtx;
	}
      }
      
      //cout << "==============================\nPV searches are done...\n===================================\n";

      TLorentzVector ele4, bcand;
      ele4.SetXYZM(iE->px(), iE->py(), iE->pz(), el_mass);
      bcand=mu4+ele4;

      auto bs_point = Vertex::Point( recoBeamSpotHandle -> x(pv.z()), recoBeamSpotHandle -> y(pv.z()), recoBeamSpotHandle -> z0() );
      auto bs_error = recoBeamSpotHandle -> covariance3D();
      float chi2_temp=0, ndof_temp=0;
      auto bs_temp = Vertex(bs_point, bs_error, chi2_temp, ndof_temp, 3);
      
      TVector3 vect_lxy, vect_pt;
      vect_lxy.SetXYZ(emuVtxFit.position().x() - bs_temp.position().x(),
		      emuVtxFit.position().y() - bs_temp.position().y(),
		      0. );
      vect_pt.SetXYZ(bcand.Px(), bcand.Py(), 0. );
      
      double vtxCos = vect_pt.Dot(vect_lxy) / ( vect_pt.Mag() * vect_lxy.Mag() );

      h1_cos2d_1 -> Fill(vtxCos);
      if( vtxCos < 0.99 ) continue;
      nn[14]++;
      h1_cos2d_2 -> Fill(vtxCos);

      h1_mass_1 -> Fill(bcand.M());
      if ( bcand.M()<4. || bcand.M()>7. ) continue;
      nn[15]++;
      h1_mass_2 -> Fill(bcand.M());

      for(std::vector<pat::PackedCandidate>::const_iterator iPFS = pfs->begin(); iPFS != pfs->end(); ++iPFS){
	if( test_event != event ) break;
	if (iPFS -> charge() == 0 ) continue;
	if (iPFS -> pt()< 1.) continue;
	if (!iPFS -> hasTrackDetails()) continue;
	if( (int)iPFS -> vertexRef().key() != pvIndex ) continue;
	charge_pt -> push_back(iPFS -> pt());
      }
      test_event++;
      //cout << "*****************************************************\nOne round is done fully\n*******************************************\n";
      cos2d -> push_back( vtxCos );

      bs_x -> push_back(bs_temp.position().x());
      bs_y -> push_back(bs_temp.position().y());
      vx -> push_back(emuVtxFit.position().x());
      vy -> push_back(emuVtxFit.position().y());
      vz -> push_back(emuVtxFit.position().z());
      vtx_chi2 -> push_back(emuVtxFit.normalisedChiSquared());
      vtx_prob -> push_back(TMath::Prob(emuVtxFit.totalChiSquared(), emuVtxFit.degreesOfFreedom()));
      // TMath::Prob(fChisquare,fNpfits-fNpar); // https://root.cern.ch/root/roottalk/roottalk01/4648.html 

      auto llxy= VertexDistanceXY().distance(bs_temp, emuVtxFit.vertexState());
      lxy -> push_back(llxy.value());
      lxy_err -> push_back(llxy.error());
      lxy_sig -> push_back(llxy.significance());

      pv_x -> push_back(pv.position().x());
      pv_y -> push_back(pv.position().y());
      pv_z -> push_back(pv.position().z());

      //trgIndex -> push_back(trg_flag);
      mass  -> push_back( bcand.M());
      pt    -> push_back( bcand.Pt());
      eta   -> push_back( bcand.Eta());
      phi   -> push_back( bcand.Phi());
      dr_mu -> push_back( bcand.DeltaR(mu4));
      dr_ele-> push_back( bcand.DeltaR(ele4));

      mu_pt       -> push_back( iMu -> pt());
      mu_eta      -> push_back( iMu -> eta());
      mu_phi      -> push_back( iMu -> phi());
      mu_e        -> push_back( iMu -> energy());
      mu_mass     -> push_back( iMu -> mass());
      mu_charge   -> push_back( iMu -> charge());
      mu_id_loose -> push_back( iMu -> isLooseMuon());
      mu_id_soft  -> push_back( iMu -> isSoftMuon(pv));
      mu_id_medium-> push_back( iMu -> isMediumMuon());
      mu_id_tight -> push_back( iMu -> isTightMuon(pv));
      mu_id_soft_mva-> push_back( iMu -> softMvaValue());
      mu_dxy      -> push_back( iMu -> bestTrack() -> dxy(pv.position()));
      mu_dxy_e    -> push_back( iMu -> bestTrack() -> dxyError(pv.position(), pv.error()));
      mu_dxy_sig  -> push_back( iMu -> bestTrack() -> dxy(pv.position() ) / iMu -> bestTrack() -> dxyError(pv.position(), pv.error()) );
      mu_dz       -> push_back( iMu -> bestTrack() -> dz(pv.position()));
      mu_dz_e     -> push_back( iMu -> bestTrack() -> dzError());
      mu_dz_sig   -> push_back( iMu -> bestTrack() -> dz(pv.position()) / iMu -> bestTrack() -> dzError());
      mu_bs_dxy   -> push_back( iMu -> bestTrack() -> dxy(bs_temp.position()));
      mu_bs_dxy_e -> push_back( iMu -> bestTrack() -> dxyError(bs_temp.position(), bs_temp.error()));
      mu_bs_dxy_sig-> push_back( iMu -> bestTrack() -> dxy(bs_temp.position()) / iMu -> bestTrack() -> dxyError(bs_temp.position(), bs_temp.error()) );
      mu_trkIsolation -> push_back( utilityFuntions.computeTrkMuonIsolation( *iMu , *iE , *pfs , pvIndex ) );

      ele_pt -> push_back( iE -> pt());
      ele_eta -> push_back( iE -> eta());
      ele_phi -> push_back( iE -> phi());
      ele_e -> push_back( iE -> energy());
      ele_mass -> push_back( iE -> mass());
      ele_charge -> push_back( iE -> charge());
      ele_id_iso_wp90 -> push_back( iE -> electronID("mvaEleID-Fall17-iso-V2-wp90"));
      ele_id_iso_wp80 -> push_back( iE -> electronID("mvaEleID-Fall17-iso-V2-wp80"));
      ele_id_iso_wpLoose -> push_back( iE -> electronID("mvaEleID-Fall17-iso-V2-wpLoose"));
      ele_id_iso -> push_back(iE->userFloat("ElectronMVAEstimatorRun2Fall17IsoV2Values"));
      ele_id_noIso -> push_back(iE->userFloat("ElectronMVAEstimatorRun2Fall17NoIsoV2Values"));
      ele_id_noIso_wp90 -> push_back( iE -> electronID("mvaEleID-Fall17-noIso-V2-wp90"));
      ele_id_noIso_wp80 -> push_back( iE -> electronID("mvaEleID-Fall17-noIso-V2-wp80"));
      ele_id_noIso_wpLoose -> push_back( iE -> electronID("mvaEleID-Fall17-noIso-V2-wpLoose"));
      ele_dxy -> push_back( iE -> gsfTrack() -> dxy(pv.position()));
      ele_dxy_e -> push_back( iE -> gsfTrack() -> dxyError(pv.position(),  pv.error() ) );
      ele_dxy_sig -> push_back( iE -> gsfTrack() -> dxy(pv.position()) / iE -> gsfTrack() -> dxyError(pv.position(),  pv.error()));
      ele_dz -> push_back( iE -> gsfTrack() -> dz(pv.position()));
      ele_dz_e -> push_back( iE -> gsfTrack() -> dzError() ); //pv.position(),  pv.error() ) );                                                                                                                    
      ele_dz_sig -> push_back( iE -> gsfTrack() -> dz(pv.position()) / iE -> gsfTrack() -> dzError());
      ele_bs_dxy -> push_back( iE -> gsfTrack() -> dxy(bs_temp.position()));
      ele_bs_dxy_e -> push_back( iE -> gsfTrack() -> dxyError(bs_temp.position(), bs_temp.error()));
      ele_bs_dxy_sig -> push_back( iE -> gsfTrack() -> dxy(bs_temp.position()) / iE -> gsfTrack() -> dxyError(bs_temp.position(), bs_temp.error()));
      ele_trkIsolation -> push_back( utilityFuntions.computeTrkEleIsolation( *iE, *iMu, *pfs, pvIndex ));

      mu_trk_chi2 -> push_back( iMu->bestTrack()->chi2() );
      mu_trk_ndof -> push_back( iMu->bestTrack()->ndof() );
      //mu_trk_prob -> push_back( TMath::Prob(iMu->bestTrack()->totalChiSquared(), iMu->bestTrack()->degreesOfFreedom() ) );                                                                                       
      //mu_trk_prob -> push_back( TMath::Prob(iMu->bestTrack()->chi2(), iMu->bestTrack()->ndof() ) );
      
      ele_sigmaietaieta -> push_back( iE -> full5x5_sigmaIetaIeta() );
      ele_trk_chi2 -> push_back( iE->gsfTrack()->chi2() );
      ele_trk_ndof -> push_back( iE->gsfTrack()->ndof() );
      //ele_trk_prob -> push_back(TMath::Prob(iE->gsfTrack()->totalChiSquared(), iE->gsfTrack()->degreesOfFreedom() ) );
      //ele_trk_prob -> push_back(TMath::Prob(iE->gsfTrack()->chi2(), iE->gsfTrack()->ndof() ) );

      nB++;
    }
  }
  ncands = pt->size();
  if(nn[0] != 0 ) h1_eff00 -> Fill(nn[0]);
  if(nn[1] != 0 ) h1_eff01 -> Fill(nn[1]);
  if(nn[2] != 0 ) h1_eff02 -> Fill(nn[2]);
  if(nn[3] != 0 ) h1_eff03 -> Fill(nn[3]);
  if(nn[4] != 0 ) h1_eff04 -> Fill(nn[4]);
  if(nn[5] != 0 ) h1_eff05 -> Fill(nn[5]);
  if(nn[6] != 0 ) h1_eff06 -> Fill(nn[6]);
  if(nn[7] != 0 ) h1_eff07 -> Fill(nn[7]);
  if(nn[8] != 0 ) h1_eff08 -> Fill(nn[8]);
  if(nn[9] != 0 ) h1_eff09 -> Fill(nn[9]);
  if(nn[10] != 0 ) h1_eff10 -> Fill(nn[10]);
  if(nn[11] != 0 ) h1_eff11 -> Fill(nn[11]);
  if(nn[12] != 0 ) h1_eff12 -> Fill(nn[12]);
  if(nn[13] != 0 ) h1_eff13 -> Fill(nn[13]);
  if(nn[14] != 0 ) h1_eff14 -> Fill(nn[14]);
  if(nn[15] != 0 ) h1_eff15 -> Fill(nn[15]);

  //trgIndex -> clear();
  t_pt     -> clear();
  t_eta    -> clear();
  t_phi    -> clear();

  if (nB > 0 ){
    tree -> Fill();
    nB = 0;
  }

  mass     -> clear();
  pt       -> clear();
  eta      -> clear();
  phi      -> clear(); 
  dr_mu    -> clear();
  dr_ele   -> clear(); 
  pv_x     -> clear();
  pv_y     -> clear(); 
  pv_z     -> clear(); 
  bs_x     -> clear();
  bs_y     -> clear();
  vx       -> clear();
  vy       -> clear();
  vz       -> clear(); 
  vtx_chi2 -> clear();
  vtx_prob -> clear();
  cos2d    -> clear();
  lxy      -> clear();
  lxy_err  -> clear();
  lxy_sig  -> clear(); 
  mu_pt    -> clear();
  mu_eta   -> clear();
  mu_phi   -> clear();
  mu_e     -> clear();
  mu_mass  -> clear();
  mu_charge-> clear();
  mu_id_loose-> clear();
  mu_id_soft -> clear();
  mu_id_medium -> clear();
  mu_id_tight  -> clear();
  mu_id_soft_mva  -> clear();
  mu_dxy     -> clear();
  mu_dxy_e   -> clear();
  mu_dxy_sig -> clear();
  mu_dz      -> clear();
  mu_dz_e    -> clear();
  mu_dz_sig  -> clear();
  mu_bs_dxy  -> clear();
  mu_bs_dxy_e  -> clear();
  mu_bs_dxy_sig  -> clear();
  mu_cov_pos_def  -> clear();
  mu_trkIsolation  -> clear(); 
  ele_pt    -> clear();
  ele_eta   -> clear();
  ele_phi   -> clear();
  ele_e     -> clear();
  ele_mass  -> clear();
  ele_charge  -> clear();
  ele_id_iso_wp90  -> clear();
  ele_id_iso_wp80  -> clear();
  ele_id_iso_wpLoose  -> clear();
  ele_id_iso -> clear();
  ele_id_noIso -> clear();
  ele_id_noIso_wp90  -> clear();
  ele_id_noIso_wp80  -> clear();
  ele_id_noIso_wpLoose  -> clear();
  ele_dxy  -> clear();
  ele_dxy_e  -> clear();
  ele_dxy_sig  -> clear();
  ele_dz  -> clear();
  ele_dz_e  -> clear();
  ele_dz_sig  -> clear();
  ele_bs_dxy  -> clear();
  ele_bs_dxy_e  -> clear();
  ele_bs_dxy_sig  -> clear();
  ele_cov_pos_def  -> clear();
  ele_trkIsolation  -> clear(); 

  mu_trk_chi2-> clear();
  mu_trk_ndof-> clear();
  mu_trk_prob -> clear();
  ele_sigmaietaieta-> clear();
  ele_trk_chi2-> clear();
  ele_trk_ndof-> clear();
  ele_trk_prob -> clear();
  charge_pt ->clear();

  //#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  //ESHandle<SetupData> pSetup;
  //iSetup.get<SetupRecord>().get(pSetup);
  //#endif
}
  
void BsToEMu::beginJob() {
  std::cout << "Beginning analyzer with isMC = " << isMC_ << std::endl;
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("ntuple","Bs -> eMu ntuple");
  tree -> Branch("nB", &nB, "nB/i" );
  tree -> Branch("run", &run );
  tree -> Branch("lumi", &lumi );
  tree -> Branch("event", &event );
  tree -> Branch("npv", &npv );
  tree -> Branch("ncands", &ncands );

  tree -> Branch("hlt_path", &hlt_path);
  //tree -> Branch("trgIndex", &trgIndex );
  tree -> Branch("t_pt", &t_pt);
  tree -> Branch("t_eta", &t_eta);
  tree -> Branch("t_phi", &t_phi);  

  tree -> Branch("mass", &mass );
  tree -> Branch("pt", &pt );
  tree -> Branch("eta", &eta );
  tree -> Branch("phi", &phi );
  tree -> Branch("dr_mu", &dr_mu );
  tree -> Branch("dr_ele", &dr_ele );
  tree -> Branch("pv_x", &pv_x );
  tree -> Branch("pv_y", &pv_y );
  tree -> Branch("pv_z", &pv_z );
  tree -> Branch("bs_x0", &bs_x0 );
  tree -> Branch("bs_y0", &bs_y0 );
  tree -> Branch("bs_z0", &bs_z0 );
  tree -> Branch("bs_x", &bs_x );
  tree -> Branch("bs_y", &bs_y );
  tree -> Branch("vx", &vx );
  tree -> Branch("vy", &vy );
  tree -> Branch("vz", &vz );
  tree -> Branch("vtx_chi2", &vtx_chi2 );
  tree -> Branch("vtx_prob", &vtx_prob );
  tree -> Branch("cos2d", &cos2d );
  tree -> Branch("lxy", &lxy );
  tree -> Branch("lxy_err", &lxy_err );
  tree -> Branch("lxy_sig", &lxy_sig );
  tree -> Branch("mu_pt", &mu_pt );
  tree -> Branch("mu_eta", &mu_eta );
  tree -> Branch("mu_phi", &mu_phi );
  tree -> Branch("mu_e", &mu_e );
  tree -> Branch("mu_mass", &mu_mass );
  tree -> Branch("mu_charge", &mu_charge );
  tree -> Branch("mu_id_loose", &mu_id_loose );
  tree -> Branch("mu_id_soft", &mu_id_soft );
  tree -> Branch("mu_id_medium", &mu_id_medium );
  tree -> Branch("mu_id_tight", &mu_id_tight );
  tree -> Branch("mu_id_soft_mva", &mu_id_soft_mva );
  tree -> Branch("mu_dxy", &mu_dxy );
  tree -> Branch("mu_dxy_e", &mu_dxy_e );
  tree -> Branch("mu_dxy_sig", &mu_dxy_sig );
  tree -> Branch("mu_dz", &mu_dz );
  tree -> Branch("mu_dz_e", &mu_dz_e );
  tree -> Branch("mu_dz_sig", &mu_dz_sig );
  tree -> Branch("mu_bs_dxy", &mu_bs_dxy );
  tree -> Branch("mu_bs_dxy_e", &mu_bs_dxy_e );
  tree -> Branch("mu_bs_dxy_sig", &mu_bs_dxy_sig );
  tree -> Branch("mu_cov_pos_def", &mu_cov_pos_def );
  tree -> Branch("mu_trkIsolation", &mu_trkIsolation );
  tree -> Branch("ele_pt", &ele_pt );
  tree -> Branch("ele_eta", &ele_eta );
  tree -> Branch("ele_phi", &ele_phi );
  tree -> Branch("ele_e", &ele_e );
  tree -> Branch("ele_mass", &ele_mass );
  tree -> Branch("ele_charge", &ele_charge );
  tree -> Branch("ele_id_iso_wp90", &ele_id_iso_wp90 );
  tree -> Branch("ele_id_iso_wp80", &ele_id_iso_wp80 );
  tree -> Branch("ele_id_iso_wpLoose", &ele_id_iso_wpLoose );
  tree -> Branch("ele_id_iso", &ele_id_iso );
  tree -> Branch("ele_id_noIso", &ele_id_noIso );
  tree -> Branch("ele_id_noIso_wp90", &ele_id_noIso_wp90 );
  tree -> Branch("ele_id_noIso_wp80", &ele_id_noIso_wp80 );
  tree -> Branch("ele_id_noIso_wpLoose", &ele_id_noIso_wpLoose );
  tree -> Branch("ele_dxy", &ele_dxy );
  tree -> Branch("ele_dxy_e", &ele_dxy_e );
  tree -> Branch("ele_dxy_sig", &ele_dxy_sig );
  tree -> Branch("ele_dz", &ele_dz );
  tree -> Branch("ele_dz_e", &ele_dz_e );
  tree -> Branch("ele_dz_sig", &ele_dz_sig );
  tree -> Branch("ele_bs_dxy", &ele_bs_dxy );
  tree -> Branch("ele_bs_dxy_e", &ele_bs_dxy_e );
  tree -> Branch("ele_bs_dxy_sig", &ele_bs_dxy_sig );
  tree -> Branch("ele_cov_pos_def", &ele_cov_pos_def );
  tree -> Branch("ele_trkIsolation", &ele_trkIsolation );

  tree -> Branch("mu_trk_chi2", &mu_trk_chi2 );
  tree -> Branch("mu_trk_ndof", &mu_trk_ndof );
  tree -> Branch("mu_trk_prob", &mu_trk_prob );
  tree -> Branch("ele_sigmaietaieta", &ele_sigmaietaieta );
  tree -> Branch("ele_trk_chi2", &ele_trk_chi2 );
  tree -> Branch("ele_trk_ndof", &ele_trk_ndof );
  tree -> Branch("ele_trk_prob", &ele_trk_prob );

  tree -> Branch("charge_pt", &charge_pt );
}
void BsToEMu::endJob() {
  tree -> GetDirectory()->cd();
  tree -> Write();
  h1_eff00 -> Write();
  h1_eff01 -> Write();
  h1_eff02 -> Write();
  h1_eff03 -> Write();
  h1_eff04 -> Write();
  h1_eff05 -> Write();
  h1_eff06 -> Write();
  h1_eff07 -> Write();
  h1_eff08 -> Write();
  h1_eff09 -> Write();
  h1_eff10 -> Write();
  h1_eff11 -> Write();
  h1_eff12 -> Write();
  h1_eff13 -> Write();
  h1_eff14 -> Write();
  h1_eff15 -> Write();

  h1_mu_pt_1     -> Write();
  h1_mu_pt_2     -> Write();
  h1_mu_pt_3     -> Write();
  h1_mu_pt_4     -> Write();
  h1_mu_eta_1    -> Write();
  h1_mu_eta_2    -> Write();
  h1_mu_eta_3    -> Write();
  h1_mu_eta_4    -> Write();
  h1_delPt_b1    -> Write();
  h1_delPt_b2    -> Write();
  h1_delPt_a     -> Write();
  h1_delRtrg_b1  -> Write();
  h1_delRtrg_b2  -> Write();
  h1_delRtrg_a   -> Write();
  h1_el_pt_1     -> Write();
  h1_el_pt_2     -> Write();
  h1_el_pt_3     -> Write();
  h1_el_pt_4     -> Write();
  h1_el_eta_1    -> Write();
  h1_el_eta_2    -> Write();
  h1_el_eta_3    -> Write();
  h1_el_eta_4    -> Write();
  h1_el_idIso_1  -> Write();
  h1_el_idNoIso_1-> Write();
  h1_el_idIso_2  -> Write();
  h1_el_idNoIso_2-> Write();
  h1_cos2d_1     -> Write();
  h1_cos2d_2     -> Write();
  h1_mass_1      -> Write();
  h1_mass_2      -> Write();
}
void BsToEMu::fillDescriptions(edm::ConfigurationDescriptions& descriptions) { 
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(BsToEMu);
