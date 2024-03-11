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

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "TLorentzVector.h"
#include "KVFitter.h"
#include "TFile.h"
#include "TTree.h"
//using reco::TrackCollection;
#include "Util.h"

class Bu2JPsiK : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit Bu2JPsiK(const edm::ParameterSet&);
  ~Bu2JPsiK();

  bool IsTheSame(const pat::PackedCandidate &tk, const pat::Muon &mu);
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  void fillGenBranches( const edm::Event& );

  KVFitter vtxfit;

  // ----------member data ---------------------------
  //edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
  //edm::EDGetTokenT<edm::View<pat::Electron>> electronToken_;
  edm::EDGetTokenT<edm::View<pat::Muon>> muonToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
  edm::EDGetTokenT<reco::VertexCollection> vertices_;
  edm::EDGetTokenT<reco::VertexCollection> ltrk_;
  edm::EDGetTokenT<edm::TriggerResults> trg_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> tobjs_;

  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesCollection_;

  bool isMC_, OnlyGen_;
  TTree *tree, *tree_gen;
  unsigned int nB;

  int run,lumi,event, npv,ncands, hlt_path;
  std::vector<int>   *trgIndex;
  std::vector<float> *t_pt, *t_eta, *t_phi;
  std::vector<float> *mass, *mcorr, *pt, *eta, *phi;
  std::vector<int>   *charge;
  std::vector<float> *pv_x, *pv_y, *pv_z; float bs_x0, bs_y0, bs_z0; std::vector<float> *bs_x, *bs_y, *vx, *vy, *vz;
  std::vector<float> *J_mass, *J_pt, *J_eta, *J_phi;
  std::vector<float> *vtx_chi2, *vtx_prob, *cos2d, *lxy, *lxy_err, *lxy_sig, *l3d, *l3d_err, *l3d_sig;
  std::vector<float> *mu1_pt, *mu1_eta, *mu1_phi, *mu1_mass, *mu1_charge, *mu1_id_loose, *mu1_id_soft, *mu1_id_medium, *mu1_id_tight, *mu1_id_soft_mva, *mu1_dxy, *mu1_dxy_e, *mu1_dxy_sig, *mu1_dz, *mu1_dz_e, *mu1_dz_sig, *mu1_bs_dxy, *mu1_bs_dxy_e, *mu1_bs_dxy_sig, *mu1_cov_pos_def, *mu1_trk_chi2, *mu1_trk_ndof, *mu1_trk_prob, *dr_mu1;
  std::vector<float> *mu2_pt, *mu2_eta, *mu2_phi, *mu2_mass, *mu2_charge, *mu2_id_loose, *mu2_id_soft, *mu2_id_medium, *mu2_id_tight, *mu2_id_soft_mva, *mu2_dxy, *mu2_dxy_e, *mu2_dxy_sig, *mu2_dz, *mu2_dz_e, *mu2_dz_sig, *mu2_bs_dxy, *mu2_bs_dxy_e, *mu2_bs_dxy_sig, *mu2_cov_pos_def, *mu2_trk_chi2, *mu2_trk_ndof, *mu2_trk_prob, *dr_mu2;
  std::vector<float> *trk_pt, *trk_eta, *trk_phi;
  std::vector<int>   *trk_charge;
  std::vector<float> *trk_dxy, *trk_dz, *trk_dxy_sig, *trk_dz_sig, *trk_innerHits;
  std::vector<float> *cosAl;

  std::vector<float> *c_mu1_pt, *c_mu1_eta, *c_mu1_phi, *c_mu2_pt, *c_mu2_eta, *c_mu2_phi, *c_trk_pt, *c_trk_eta, *c_trk_phi, *c_J_mass, *c_J_pt, *c_J_eta, *c_J_phi, *c_mass, *c_pt, *c_eta, *c_phi;

  int gen_nB, gen_trks;
  std::vector<int>  *gen_B_pid, *gen_BMuM_pid, *gen_BMuP_pid, *gen_BK_pid, *gen_BJpsi_pid, *gen_otrk_pid;
  std::vector<float> *gen_B_pt, *gen_B_energy, *gen_B_eta, *gen_B_phi, *gen_B_pz;
  std::vector<float> *gen_BMuM_pt, *gen_BMuM_eta, *gen_BMuM_phi;
  std::vector<float> *gen_BMuP_pt, *gen_BMuP_eta, *gen_BMuP_phi;
  std::vector<float> *gen_BK_pt, *gen_BK_energy, *gen_BK_eta, *gen_BK_phi;
  std::vector<float> *gen_BJpsi_pt, *gen_BJpsi_energy, *gen_BJpsi_eta, *gen_BJpsi_phi;
  std::vector<float> *gen_otrk_pt, *gen_otrk_eta, *gen_otrk_phi;
};
Bu2JPsiK::Bu2JPsiK(const edm::ParameterSet& iConfig):
  muonToken_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
  beamSpotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
  pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
  vertices_ (consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  ltrk_ (consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("ltrk"))),
  trg_ (consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("trg_res"))),
  tobjs_ (consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("tobjs"))),
  
  genParticlesCollection_ (consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"))),

  isMC_(iConfig.getParameter<bool>("isMC")),
  OnlyGen_(iConfig.getParameter<bool>("OnlyGen")),
  //now do what ever initialization is needed
  tree(0), tree_gen(0),
  nB(0),
  run(0), lumi(0), event(0), npv(0), ncands(0), hlt_path(0), trgIndex(0),
  t_pt(0), t_eta(0), t_phi(0),
  mass(0), mcorr(0), pt(0), eta(0), phi(0), charge(0), 
  pv_x(0), pv_y(0), pv_z(0), bs_x0(0), bs_y0(0), bs_z0(0), bs_x(0), bs_y(0), vx(0), vy(0), vz(0),
  J_mass(0), J_pt(0), J_eta(0), J_phi(0),
  vtx_chi2(0), vtx_prob(0), cos2d(0), lxy(0), lxy_err(0), lxy_sig(0), l3d(0), l3d_err(0), l3d_sig(0),
  mu1_pt(0), mu1_eta(0), mu1_phi(0), mu1_mass(0), mu1_charge(0), mu1_id_loose(0), mu1_id_soft(0), mu1_id_medium(0), mu1_id_tight(0), mu1_id_soft_mva(0), mu1_dxy(0), mu1_dxy_e(0), mu1_dxy_sig(0), mu1_dz(0), mu1_dz_e(0), mu1_dz_sig(0), mu1_bs_dxy(0), mu1_bs_dxy_e(0), mu1_bs_dxy_sig(0), mu1_cov_pos_def(0), mu1_trk_chi2(0), mu1_trk_ndof(0), mu1_trk_prob(0), dr_mu1(0),
  mu2_pt(0), mu2_eta(0), mu2_phi(0), mu2_mass(0), mu2_charge(0), mu2_id_loose(0), mu2_id_soft(0), mu2_id_medium(0), mu2_id_tight(0), mu2_id_soft_mva(0), mu2_dxy(0), mu2_dxy_e(0), mu2_dxy_sig(0), mu2_dz(0), mu2_dz_e(0), mu2_dz_sig(0), mu2_bs_dxy(0), mu2_bs_dxy_e(0), mu2_bs_dxy_sig(0), mu2_cov_pos_def(0), mu2_trk_chi2(0), mu2_trk_ndof(0), mu2_trk_prob(0), dr_mu2(0),
  trk_pt(0), trk_eta(0), trk_phi(0), 
  trk_charge(0),
  trk_dxy(0), trk_dz(0), trk_dxy_sig(0), trk_dz_sig(0), trk_innerHits(0),
  cosAl(0),

  c_mu1_pt(0), c_mu1_eta(0), c_mu1_phi(0), c_mu2_pt(0), c_mu2_eta(0), c_mu2_phi(0), c_trk_pt(0), c_trk_eta(0), c_trk_phi(0),
  c_J_mass(0), c_J_pt(0), c_J_eta(0), c_J_phi(0), c_mass(0), c_pt(0), c_eta(0), c_phi(0),

  gen_nB(0), gen_trks(0),
  gen_B_pid(0), gen_BMuM_pid(0), gen_BMuP_pid(0), gen_BK_pid(0), gen_BJpsi_pid(0), gen_otrk_pid(0),
  gen_B_pt(0), gen_B_energy(0), gen_B_eta(0), gen_B_phi(0), gen_B_pz(0),
  gen_BMuM_pt(0), gen_BMuM_eta(0), gen_BMuM_phi(0),
  gen_BMuP_pt(0), gen_BMuP_eta(0), gen_BMuP_phi(0),
  gen_BK_pt(0), gen_BK_energy(0), gen_BK_eta(0), gen_BK_phi(0),
  gen_BJpsi_pt(0), gen_BJpsi_energy(0), gen_BJpsi_eta(0), gen_BJpsi_phi(0),
  gen_otrk_pt(0), gen_otrk_eta(0), gen_otrk_phi(0)

{}

Bu2JPsiK::~Bu2JPsiK(){}

void Bu2JPsiK::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  using std::vector;
  using namespace edm;
  using namespace reco;
  using namespace std;
  //std::cout<<__LINE__<<"\n"; 

  if( isMC_ ) {
    fillGenBranches(iEvent);
  }
  if( OnlyGen_ ) return;

  edm::Handle<pat::PackedCandidateCollection> pfs;
  iEvent.getByToken(pfToken_, pfs);
    
  edm::Handle<VertexCollection> primaryVertices;
  iEvent.getByToken(vertices_, primaryVertices);
    
  edm::Handle<VertexCollection> ltrk;
  iEvent.getByToken(ltrk_, ltrk);
    
  edm::Handle<edm::TriggerResults> trg_res;
  iEvent.getByToken(trg_, trg_res);
    
  edm::Handle<pat::TriggerObjectStandAloneCollection> tobjs;
  iEvent.getByToken(tobjs_, tobjs);
    
  edm::Handle< View<pat::Muon> > muons;
  iEvent.getByToken(muonToken_, muons );
    
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  iEvent.getByToken(beamSpotToken_, recoBeamSpotHandle );

  edm::ESHandle<TransientTrackBuilder> theB; 
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB); 
    
  lumi = iEvent.id().luminosityBlock();
  run = iEvent.id().run();
  event = iEvent.id().event();
  npv = primaryVertices->size();//iEvent.id().size();
    
  bs_x0 = recoBeamSpotHandle -> x0();
  bs_y0 = recoBeamSpotHandle -> y0();
  bs_z0 = recoBeamSpotHandle -> z0();
    
  static int nhlt = 10;
  const char *paths[nhlt]
    = { "HLT_Mu7_IP4"     ,
	"HLT_Mu8_IP3"     ,
	"HLT_Mu8_IP5"     ,
	"HLT_Mu8_IP6"     ,
	"HLT_Mu8p5_IP3p5" ,
	"HLT_Mu9_IP4"     ,
	"HLT_Mu9_IP5"     ,
	"HLT_Mu9_IP6"     ,
	"HLT_Mu10p5_IP3p5",
	"HLT_Mu12_IP6"  };
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
  if( !hlt_passed ) return; 
  hlt_path=hlt_passed;
  pat::TriggerObjectStandAlone good_tobjs;
  for(pat::TriggerObjectStandAlone obj : *tobjs) {
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
      t_pt ->push_back(obj.pt());
      t_eta->push_back(obj.eta());
      t_phi->push_back(obj.phi());
    }
  }
  std::vector <reco::Track> tracksTifit, jpsifit;
  for(View<pat::Muon>::const_iterator iMu = muons->begin(); iMu != muons->end(); ++iMu){
    if( iMu  -> pt() < 6.0 || fabs(iMu -> eta()) > 2.5 ) continue;
    if( ! iMu -> bestTrack()) continue;
    if( ! iMu -> isPFMuon() ) continue;
    if( ! iMu -> isGlobalMuon() ) continue;
    if( ! iMu -> isLooseMuon() ) continue;
    int trg_flag=-999;
    for(int trg_i = 0; trg_i < (int)t_eta -> size(); trg_i++){
      double deta = iMu -> eta() - t_eta -> at(0);
      double dphi = fabs(iMu -> phi() - t_phi->at(0));
      if(dphi > M_PI) dphi =  2*M_PI - dphi;
      double t_dr = sqrt(deta*deta+dphi*dphi);
      if( t_dr > 0.02 ) continue; 
	
      double dpT = fabs((iMu  -> pt() - t_pt->at(0))/iMu  -> pt());
      if( dpT > 0.05 ) continue;
      trg_flag=trg_i; break;
    }
    if(trg_flag < 0) continue;
      
    auto pv = (*primaryVertices)[0];
    //int pvIndex=0;
    double dz_mu0 = fabs(iMu -> bestTrack() -> dz(pv.position()));
    for(int ivtx = 1; ivtx < npv; ++ivtx ){
      double dz_mu1 = fabs(iMu -> bestTrack() -> dz((*primaryVertices)[ivtx].position()));
      if( dz_mu0 > dz_mu1 ){
	dz_mu0 = dz_mu1;
	pv = (*primaryVertices)[ivtx];
	//pvIndex = ivtx;
      }
    }
      
    for(View<pat::Muon>::const_iterator jMu = muons->begin(); jMu != muons->end(); ++jMu){
      if(jMu->pt() < 4. || fabs(jMu->eta()) > 2.5 ) continue;
      if(iMu->charge() + jMu->charge() != 0 ) continue;
	
      ParticleMass muon_mass = 0.10565837;
      ParticleMass psi_mass = 3.096916;
      float muon_sigma = muon_mass*1.e-6;
      TLorentzVector mu1L, mu2L, JPsiCand;
      mu1L.SetXYZM(iMu->px(), iMu->py(), iMu->pz(), muon_mass);
      mu2L.SetXYZM(jMu->px(), jMu->py(), jMu->pz(), muon_mass);
      jpsifit.clear();
      jpsifit.push_back( *(iMu->bestTrack()) );
      jpsifit.push_back( *(jMu->bestTrack()) );
      auto mu1mu2vtxFit = vtxfit.Fit(jpsifit);
      if( !mu1mu2vtxFit.isValid()) continue;
      
      JPsiCand = mu1L + mu2L;
      if( JPsiCand.M() < 2.9 || JPsiCand.M() > 3.3 ) continue;
      //reco::TransientTrack muon1TT((*theB).build(iMu->innerTrack().get()));
      //reco::TransientTrack muon2TT((*theB).build(jMu->innerTrack().get()));
      reco::TransientTrack muon1TT((*theB).build(iMu->bestTrack()));
      reco::TransientTrack muon2TT((*theB).build(jMu->bestTrack()));
      for(std::vector<pat::PackedCandidate>::const_iterator iTrack = pfs->begin(); iTrack != pfs -> end(); ++iTrack){
	if( iTrack->charge() == 0 ) continue;
	if( fabs(iTrack -> pdgId()) != 211 ) continue;
	if( iTrack->pt() < 0.95) continue;
	if(! iTrack -> trackHighPurity() ) continue;
	if(! iTrack -> hasTrackDetails() ) continue;
	if( IsTheSame(*iTrack, *iMu ) || IsTheSame(*iTrack, *jMu) ) continue;
	  
	//ParticleMass kaon_mass = 0.13957061; // pion mass
	ParticleMass kaon_mass = 0.493677; // kaon mass
	float kaon_sigma = kaon_mass*1.e-6;
	float chi = 0.;
	float ndf = 0.;

	TLorentzVector kaonL, bcand;
	kaonL.SetXYZM(iTrack->px(), iTrack->py(), iTrack->pz(), kaon_mass);
	if( (JPsiCand+kaonL).M()  < 4. || (JPsiCand+kaonL).M() > 7 ) continue;
	  
	KinematicParticleFactoryFromTransientTrack pFactory;
	reco::TransientTrack kaon1TT((*theB).build(iTrack->pseudoTrack()));
	vector<RefCountedKinematicParticle> vFitMCParticles;
	vFitMCParticles.clear();
	vFitMCParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
	vFitMCParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
	vFitMCParticles.push_back(pFactory.particle(kaon1TT,kaon_mass,chi,ndf,kaon_sigma));
	
	// JPsi mass constraint is applied in the final B+ fit,  
	MultiTrackKinematicConstraint *  j_psi_c = new  TwoTrackMassKinematicConstraint(psi_mass);
	KinematicConstrainedVertexFitter kcvFitter;
	RefCountedKinematicTree vertexFitTree = kcvFitter.fit(vFitMCParticles, j_psi_c);
	if (!vertexFitTree->isValid()) {
	  // std::cout << "caught an exception in the B vertex fit with MC" << std::endl;                                        
	  continue;
	}
	RefCountedKinematicParticle bCandMC = vertexFitTree->currentParticle();
	RefCountedKinematicVertex bDecayVertexMC = vertexFitTree->currentDecayVertex();
	if (!bDecayVertexMC->vertexIsValid()) continue; // B MC fit vertex is validity
	if(bCandMC->currentState().mass()<4.0 || bCandMC->currentState().mass()>7.) continue;
	if(bDecayVertexMC->chiSquared()<0 || bDecayVertexMC->chiSquared()>50 ) continue; //continue from negative chi2
	if(TMath::Prob(bDecayVertexMC->chiSquared(),(int)bDecayVertexMC->degreesOfFreedom() < 0.01)) continue; //std::cout << "B vtx probablity ...\n";
	vertexFitTree->movePointerToTheFirstChild();
	RefCountedKinematicParticle mu1CandMC = vertexFitTree->currentParticle();
	vertexFitTree->movePointerToTheNextChild();
	RefCountedKinematicParticle mu2CandMC = vertexFitTree->currentParticle();
	vertexFitTree->movePointerToTheNextChild();
	RefCountedKinematicParticle T1CandMC = vertexFitTree->currentParticle();
	KinematicParameters muon1KP = mu1CandMC->currentState().kinematicParameters();
	KinematicParameters muon2KP = mu2CandMC->currentState().kinematicParameters();   
	KinematicParameters kaon1KP = T1CandMC->currentState().kinematicParameters();
	TLorentzVector kaon1TL, muon1TL, muon2TL, bcanTL;
	muon1TL.SetXYZM(muon1KP.momentum().x(), muon1KP.momentum().y(), muon1KP.momentum().z(), muon_mass);
	muon2TL.SetXYZM(muon2KP.momentum().x(), muon2KP.momentum().y(), muon2KP.momentum().z(), muon_mass);
	kaon1TL.SetXYZM(kaon1KP.momentum().x(), kaon1KP.momentum().y(), kaon1KP.momentum().z(), kaon_mass);
	bcanTL.SetXYZM(bCandMC->currentState().globalMomentum().x(), bCandMC->currentState().globalMomentum().y(), bCandMC->currentState().globalMomentum().z(), bCandMC->currentState().mass());

	tracksTifit.clear();
	tracksTifit.push_back( *(iMu->bestTrack()) );
	tracksTifit.push_back( *(jMu->bestTrack()) );
	tracksTifit.push_back( *(iTrack->bestTrack() ) );
	auto mu1mu2KvtxFit = vtxfit.Fit( tracksTifit );
	if( ! mu1mu2KvtxFit.isValid()) continue;
	//mu1mu2KvtxFit.movePointerToTheFirstChild();
	//cout << "Testing... " << mu1mu2KvtxFit.vertexState().size() << endl;
	  
	bcand = mu1L + mu2L + kaonL;
	auto bs_point = Vertex::Point( recoBeamSpotHandle -> x(pv.z()), recoBeamSpotHandle -> y(pv.z()), recoBeamSpotHandle -> z0() );
	auto bs_error = recoBeamSpotHandle -> covariance3D();
	float chi2_temp=0, ndof_temp=0;
	auto bs_temp = Vertex(bs_point, bs_error, chi2_temp, ndof_temp, 3);
	TVector3 vect_lxy, vect_pt;
	vect_lxy.SetXYZ(mu1mu2KvtxFit.position().x() - bs_temp.position().x(),
			mu1mu2KvtxFit.position().y() - bs_temp.position().y(),
			0. );
	vect_pt.SetXYZ(bcand.Px(), bcand.Py(), 0. );
	double vtxCos = vect_pt.Dot(vect_lxy) / ( vect_pt.Mag() * vect_lxy.Mag() );
	  
	//if( vtxCos < 0.99 ) continue;
	if ( bcand.M()<4. || bcand.M()>7. ) continue;
	// Number of charged tracks is missing here.
	  
	// muon (triggered) momentum in rest frame of b-cand
	TVector3 b_boost;
	TLorentzVector boost_mu;
	boost_mu.SetXYZM(iMu->px(), iMu->py(), iMu->pz(), muon_mass);
	b_boost = -bcand.BoostVector();
	boost_mu.Boost(b_boost);
	// vector perpendicular to plane (b-cand momentum and beam axis)
	TVector3 perp_bpt, boost_mu_p;
	perp_bpt.SetXYZ(bcand.Py(), -bcand.Px(), 0. );
	boost_mu_p.SetXYZ(boost_mu.Px(), boost_mu.Py(), boost_mu.Pz());
	double cos_alpha_mu = perp_bpt.Dot(boost_mu_p) / (perp_bpt.Mag() * boost_mu_p.Mag());
	//double angle_mu = acos(cos_alpha_mu);
	//alpha_mu -> push_back(angle_mu);
	cosAl		-> push_back(cos_alpha_mu);
	cos2d		-> push_back( vtxCos );
	  
	pv_x		-> push_back(pv.position().x());
	pv_y		-> push_back(pv.position().y());
	pv_z		-> push_back(pv.position().z());
	  
	bs_x		-> push_back(bs_temp.position().x());
	bs_y		-> push_back(bs_temp.position().y());
	vx		-> push_back(mu1mu2KvtxFit.position().x());
	vy		-> push_back(mu1mu2KvtxFit.position().y());
	vz		-> push_back(mu1mu2KvtxFit.position().z());

	trgIndex	-> push_back(trg_flag);
	mass		-> push_back( bcand.M());
	//mcorr		-> push_back(); 
	pt		-> push_back( bcand.Pt());
	eta		-> push_back( bcand.Eta());
	phi		-> push_back( bcand.Phi() ); 
	//charge	-> push_back( bcand.Charge() );

	J_mass		-> push_back( JPsiCand.M() );
	J_pt		-> push_back( JPsiCand.Pt() );
	J_eta		-> push_back( JPsiCand.Eta() );
	J_phi		-> push_back( JPsiCand.Phi() );
	vtx_chi2	-> push_back( mu1mu2KvtxFit.normalisedChiSquared() );
	vtx_prob	-> push_back( TMath::Prob(mu1mu2KvtxFit.totalChiSquared(), mu1mu2KvtxFit.degreesOfFreedom()) );

	auto llxy= VertexDistanceXY().distance(bs_temp, mu1mu2KvtxFit.vertexState());
	lxy		-> push_back( llxy.value() );
	lxy_err		-> push_back( llxy.error() );
	lxy_sig		-> push_back( llxy.significance() );
	//l3d		-> push_back();
	//l3d_err		-> push_back();
	//l3d_sig		-> push_back();
	//cout << "(" << muon1TL.Pt() - iMu -> pt() << ", " << muon2TL.Pt() - jMu -> pt() << ")  pt difference of first muon.";
	mu1_pt		-> push_back( iMu -> pt() );
	mu1_eta		-> push_back( iMu -> eta() );
	mu1_phi		-> push_back( iMu -> phi() );
	mu1_mass	-> push_back( iMu -> mass() );
	mu1_charge	-> push_back( iMu -> charge() );
	mu1_id_loose	-> push_back( iMu -> isLooseMuon() );
	mu1_id_soft	-> push_back( iMu -> isSoftMuon(pv) );
	mu1_id_medium	-> push_back( iMu -> isMediumMuon() );
	mu1_id_tight	-> push_back( iMu -> isTightMuon(pv) );
	mu1_id_soft_mva	-> push_back( iMu -> softMvaValue() );
	mu1_dxy		-> push_back( iMu -> bestTrack() -> dxy(pv.position()) );
	mu1_dxy_e	-> push_back( iMu -> bestTrack() -> dxyError(pv.position(), pv.error()) );
	mu1_dxy_sig	-> push_back( iMu -> bestTrack() -> dxy(pv.position() ) / iMu -> bestTrack() -> dxyError(pv.position(), pv.error()) );
	mu1_dz		-> push_back( iMu -> bestTrack() -> dz(pv.position()) );
	mu1_dz_e	-> push_back( iMu -> bestTrack() -> dzError() );
	mu1_dz_sig	-> push_back( iMu -> bestTrack() -> dz(pv.position()) / iMu -> bestTrack() -> dzError() );
	mu1_bs_dxy	-> push_back( iMu -> bestTrack() -> dxy(bs_temp.position()) );
	mu1_bs_dxy_e	-> push_back( iMu -> bestTrack() -> dxyError(bs_temp.position(), bs_temp.error()) );
	mu1_bs_dxy_sig	-> push_back( iMu -> bestTrack() -> dxy(bs_temp.position()) / iMu -> bestTrack() -> dxyError(bs_temp.position(), bs_temp.error()) );
	//mu1_cov_pos_def	-> push_back();
	mu1_trk_chi2	-> push_back( iMu->bestTrack()->chi2() );
	mu1_trk_ndof	-> push_back( iMu->bestTrack()->ndof() );
	//mu1_trk_prob	-> push_back();
	dr_mu1		-> push_back( bcand.DeltaR(mu1L ) );
	mu2_pt		-> push_back( jMu -> pt() );
	mu2_eta		-> push_back( jMu -> eta() );
	mu2_phi		-> push_back( jMu -> phi() );
	mu2_mass	-> push_back( jMu -> mass() );
	mu2_charge	-> push_back( jMu -> charge() );
	mu2_id_loose	-> push_back( jMu -> isLooseMuon() );
	mu2_id_soft	-> push_back( jMu -> isSoftMuon(pv) );
	mu2_id_medium	-> push_back( jMu -> isMediumMuon() );
	mu2_id_tight	-> push_back( jMu -> isTightMuon(pv) );
	mu2_id_soft_mva	-> push_back( jMu -> softMvaValue() );
	mu2_dxy		-> push_back( jMu -> bestTrack() -> dxy(pv.position()) );
	mu2_dxy_e	-> push_back( jMu -> bestTrack() -> dxyError(pv.position(), pv.error()) );
	mu2_dxy_sig	-> push_back( jMu -> bestTrack() -> dxy(pv.position() ) / jMu -> bestTrack() -> dxyError(pv.position(), pv.error()) );
	mu2_dz		-> push_back( jMu -> bestTrack() -> dz(pv.position()) );
	mu2_dz_e	-> push_back( jMu -> bestTrack() -> dzError() );
	mu2_dz_sig	-> push_back( jMu -> bestTrack() -> dz(pv.position()) / jMu -> bestTrack() -> dzError() );
	mu2_bs_dxy	-> push_back( jMu -> bestTrack() -> dxy(bs_temp.position()) );
	mu2_bs_dxy_e	-> push_back( jMu -> bestTrack() -> dxyError(bs_temp.position(), bs_temp.error()) );
	mu2_bs_dxy_sig	-> push_back( jMu -> bestTrack() -> dxy(bs_temp.position()) / jMu -> bestTrack() -> dxyError(bs_temp.position(), bs_temp.error()) );
	//mu2_cov_pos_def	-> push_back();
	mu2_trk_chi2	-> push_back( jMu->bestTrack()->chi2() );
	mu2_trk_ndof	-> push_back( jMu->bestTrack()->ndof() );
	//mu2_trk_prob	-> push_back();
	dr_mu2		-> push_back( bcand.DeltaR(mu2L ) );
	trk_pt		-> push_back( iTrack -> pt() );
	trk_eta		-> push_back( iTrack -> eta() );
	trk_phi		-> push_back( iTrack -> phi() );
	trk_charge	-> push_back( iTrack -> charge() );
	trk_dxy		-> push_back( iTrack->dxy() );
	trk_dz		-> push_back( iTrack->dz() );
	trk_dxy_sig	-> push_back( iTrack -> dxy() / iTrack -> dxyError() );
	trk_dz_sig	-> push_back( iTrack -> dz() / iTrack -> dzError() );
	trk_innerHits	-> push_back( iTrack->lostInnerHits() );

	c_mu1_pt  -> push_back( muon1TL.Pt() );
	c_mu1_eta -> push_back( muon1TL.Eta() );
	c_mu1_phi -> push_back( muon1TL.Phi() );
	c_mu2_pt  -> push_back( muon2TL.Pt() );
	c_mu2_eta -> push_back( muon2TL.Eta() );
	c_mu2_phi -> push_back( muon2TL.Phi() );
	c_trk_pt  -> push_back( kaon1TL.Pt() );
	c_trk_eta -> push_back( kaon1TL.Eta() );
	c_trk_phi -> push_back( kaon1TL.Phi() );
	c_J_mass -> push_back( (muon1TL + muon2TL).M() );
	c_J_pt   -> push_back( (muon1TL + muon2TL).Pt() );
	c_J_eta  -> push_back( (muon1TL + muon2TL).Eta() );
	c_J_phi  -> push_back( (muon1TL + muon2TL).Phi() );
	c_mass   -> push_back( bcanTL.M() );
	c_pt     -> push_back( bcanTL.Pt() );
	c_eta    -> push_back( bcanTL.Eta() );
	c_phi    -> push_back( bcanTL.Phi() );

	nB++;
      }
    }
  }
  if( nB > 0 ){
    tree -> Fill();
    nB = 0;
  }
  trgIndex	-> clear();
  t_pt		-> clear();
  t_eta		-> clear();
  t_phi		-> clear();

  mass		-> clear(); 
  //mcorr		-> clear(); 
  pt		-> clear(); 
  eta		-> clear(); 
  phi		-> clear(); 
  //charge	-> clear();

  pv_x		-> clear(); 
  pv_y		-> clear(); 
  pv_z		-> clear();
  bs_x		-> clear();
  bs_y		-> clear();
  vx		-> clear();
  vy		-> clear();
  vz		-> clear();

  J_mass	-> clear();
  J_pt		-> clear();
  J_eta		-> clear();
  J_phi		-> clear();
  vtx_chi2	-> clear();
  vtx_prob	-> clear();
  cos2d		-> clear();
  lxy			-> clear();
  lxy_err		-> clear();
  lxy_sig		-> clear();
  //l3d			-> clear();
  //l3d_err		-> clear();
  //l3d_sig		-> clear();
  mu1_pt		-> clear();
  mu1_eta		-> clear();
  mu1_phi		-> clear();
  mu1_mass		-> clear();
  mu1_charge		-> clear();
  mu1_id_loose		-> clear();
  mu1_id_soft		-> clear();
  mu1_id_medium		-> clear();
  mu1_id_tight		-> clear();
  mu1_id_soft_mva	-> clear();
  mu1_dxy		-> clear();
  mu1_dxy_e		-> clear();
  mu1_dxy_sig		-> clear();
  mu1_dz		-> clear();
  mu1_dz_e		-> clear();
  mu1_dz_sig		-> clear();
  mu1_bs_dxy		-> clear();
  mu1_bs_dxy_e		-> clear();
  mu1_bs_dxy_sig	-> clear();
  //mu1_cov_pos_def	-> clear();
  mu1_trk_chi2		-> clear();
  mu1_trk_ndof		-> clear();
  //mu1_trk_prob		-> clear();
  dr_mu1		-> clear();
  mu2_pt		-> clear();
  mu2_eta		-> clear();
  mu2_phi		-> clear();
  mu2_mass		-> clear();
  mu2_charge		-> clear();
  mu2_id_loose		-> clear();
  mu2_id_soft		-> clear();
  mu2_id_medium		-> clear();
  mu2_id_tight		-> clear();
  mu2_id_soft_mva	-> clear();
  mu2_dxy		-> clear();
  mu2_dxy_e		-> clear();
  mu2_dxy_sig		-> clear();
  mu2_dz		-> clear();
  mu2_dz_e		-> clear();
  mu2_dz_sig		-> clear();
  mu2_bs_dxy		-> clear();
  mu2_bs_dxy_e		-> clear();
  mu2_bs_dxy_sig	-> clear();
  //mu2_cov_pos_def	-> clear();
  mu2_trk_chi2		-> clear();
  mu2_trk_ndof		-> clear();
  //mu2_trk_prob		-> clear();
  dr_mu2		-> clear();
  trk_pt		-> clear();
  trk_eta		-> clear();
  trk_phi		-> clear();
  trk_charge		-> clear();
  trk_dxy		-> clear();
  trk_dz		-> clear();
  trk_dxy_sig		-> clear();
  trk_dz_sig		-> clear();
  trk_innerHits		-> clear();
  cosAl			-> clear();

  c_mu1_pt   -> clear();
  c_mu1_eta  -> clear();
  c_mu1_phi  -> clear();
  c_mu2_pt   -> clear();
  c_mu2_eta  -> clear();
  c_mu2_phi  -> clear();
  c_trk_pt   -> clear();
  c_trk_eta  -> clear();
  c_trk_phi  -> clear();
  c_J_mass -> clear();
  c_J_pt   -> clear();
  c_J_eta  -> clear();
  c_J_phi  -> clear();
  c_mass   -> clear();
  c_pt     -> clear();
  c_eta    -> clear();
  c_phi    -> clear();

// #ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
//    ESHandle<SetupData> pSetup;
//    iSetup.get<SetupRecord>().get(pSetup);
// #endif
}
bool Bu2JPsiK::IsTheSame(const pat::PackedCandidate &tk, const pat::Muon &mu){
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}

void Bu2JPsiK::beginJob(){
  std::cout << "Beginning analyzer with isMC = " << isMC_ << std::endl;
  edm::Service<TFileService> fs;
  if( ! OnlyGen_ ){
    tree = fs->make<TTree>("ntuple","B+ -> JPsiK+ ntuple");
    tree -> Branch("nB",		&nB, "nB/i" );
    tree -> Branch("run",		&run );
    tree -> Branch("lumi",	&lumi );
    tree -> Branch("event",	&event );
    tree -> Branch("npv",		&npv );
    tree -> Branch("ncands",	&ncands );

    tree -> Branch("hlt_path",	&hlt_path);
    tree -> Branch("trgIndex",	&trgIndex );
    tree -> Branch("t_pt",	&t_pt);
    tree -> Branch("t_eta",	&t_eta);
    tree -> Branch("t_phi",	&t_phi);  

    tree -> Branch("mass",	&mass );
    tree -> Branch("pt",		&pt );
    tree -> Branch("eta",		&eta );
    tree -> Branch("phi",		&phi );
    //tree -> Branch( "charge",	&charge );

    tree -> Branch( "pv_x",	&pv_x );
    tree -> Branch( "pv_y",	&pv_y );
    tree -> Branch( "pv_z",	&pv_z );
    tree -> Branch( "bs_x0",	&bs_x0 );
    tree -> Branch( "bs_y0",	&bs_y0 );
    tree -> Branch( "bs_z0",	&bs_z0 );
    tree -> Branch( "bs_x",	&bs_x );
    tree -> Branch( "bs_y",	&bs_y );
    tree -> Branch( "vx",		&vx );
    tree -> Branch( "vy",		&vy );
    tree -> Branch( "vz",		&vz );

    tree -> Branch( "J_mass",		&J_mass );
    tree -> Branch( "J_pt",		&J_pt );
    tree -> Branch( "J_eta",		&J_eta );
    tree -> Branch( "J_phi",		&J_phi );
    tree -> Branch( "vtx_chi2",		&vtx_chi2 );
    tree -> Branch( "vtx_prob",		&vtx_prob );
    tree -> Branch( "cos2d",		&cos2d );
    tree -> Branch( "lxy",		&lxy );
    tree -> Branch( "lxy_err",		&lxy_err );
    tree -> Branch( "lxy_sig",		&lxy_sig );
    //tree -> Branch( "l3d",		&l3d );
    //tree -> Branch( "l3d_err",		&l3d_err );
    //tree -> Branch( "l3d_sig",		&l3d_sig );
    tree -> Branch( "mu1_pt",		&mu1_pt );
    tree -> Branch( "mu1_eta",		&mu1_eta );
    tree -> Branch( "mu1_phi",		&mu1_phi );
    tree -> Branch( "mu1_mass",		&mu1_mass );
    tree -> Branch( "mu1_charge",		&mu1_charge );
    tree -> Branch( "mu1_id_loose",	&mu1_id_loose );
    tree -> Branch( "mu1_id_soft",	&mu1_id_soft );
    tree -> Branch( "mu1_id_medium",	&mu1_id_medium );
    tree -> Branch( "mu1_id_tight",	&mu1_id_tight );
    tree -> Branch( "mu1_id_soft_mva",	&mu1_id_soft_mva );
    tree -> Branch( "mu1_dxy",		&mu1_dxy );
    tree -> Branch( "mu1_dxy_e",		&mu1_dxy_e );
    tree -> Branch( "mu1_dxy_sig",	&mu1_dxy_sig );
    tree -> Branch( "mu1_dz",		&mu1_dz );
    tree -> Branch( "mu1_dz_e",		&mu1_dz_e );
    tree -> Branch( "mu1_dz_sig",		&mu1_dz_sig );
    tree -> Branch( "mu1_bs_dxy",		&mu1_bs_dxy );
    tree -> Branch( "mu1_bs_dxy_e",	&mu1_bs_dxy_e );
    tree -> Branch( "mu1_bs_dxy_sig",	&mu1_bs_dxy_sig );
    //tree -> Branch( "mu1_cov_pos_def",	&mu1_cov_pos_def );
    tree -> Branch( "mu1_trk_chi2",	&mu1_trk_chi2 );
    tree -> Branch( "mu1_trk_ndof",	&mu1_trk_ndof );
    //tree -> Branch( "mu1_trk_prob",	&mu1_trk_prob );
    tree -> Branch( "dr_mu1",		&dr_mu1 );
    tree -> Branch( "mu2_pt",		&mu2_pt );
    tree -> Branch( "mu2_eta",		&mu2_eta );
    tree -> Branch( "mu2_phi",		&mu2_phi );
    tree -> Branch( "mu2_mass",		&mu2_mass );
    tree -> Branch( "mu2_charge",		&mu2_charge );
    tree -> Branch( "mu2_id_loose",	&mu2_id_loose );
    tree -> Branch( "mu2_id_soft",	&mu2_id_soft );
    tree -> Branch( "mu2_id_medium",	&mu2_id_medium );
    tree -> Branch( "mu2_id_tight",	&mu2_id_tight );
    tree -> Branch( "mu2_id_soft_mva",	&mu2_id_soft_mva );
    tree -> Branch( "mu2_dxy",		&mu2_dxy );
    tree -> Branch( "mu2_dxy_e",		&mu2_dxy_e );
    tree -> Branch( "mu2_dxy_sig",	&mu2_dxy_sig );
    tree -> Branch( "mu2_dz",		&mu2_dz );
    tree -> Branch( "mu2_dz_e",		&mu2_dz_e );
    tree -> Branch( "mu2_dz_sig",		&mu2_dz_sig );
    tree -> Branch( "mu2_bs_dxy",		&mu2_bs_dxy );
    tree -> Branch( "mu2_bs_dxy_e",	&mu2_bs_dxy_e );
    tree -> Branch( "mu2_bs_dxy_sig",	&mu2_bs_dxy_sig );
    //tree -> Branch( "mu2_cov_pos_def",	&mu2_cov_pos_def );
    tree -> Branch( "mu2_trk_chi2",	&mu2_trk_chi2 );
    tree -> Branch( "mu2_trk_ndof",	&mu2_trk_ndof );
    //tree -> Branch( "mu2_trk_prob",	&mu2_trk_prob );
    tree -> Branch( "dr_mu2",		&dr_mu2 );
    tree -> Branch( "trk_pt",		&trk_pt );
    tree -> Branch( "trk_eta",		&trk_eta );
    tree -> Branch( "trk_phi",		&trk_phi );
    tree -> Branch( "trk_charge",		&trk_charge );
    tree -> Branch( "trk_dxy",		&trk_dxy );
    tree -> Branch( "trk_dz",		&trk_dz );
    tree -> Branch( "trk_dxy_sig",	&trk_dxy_sig );
    tree -> Branch( "trk_dz_sig",		&trk_dz_sig );
    tree -> Branch( "trk_innerHits",	&trk_innerHits );
    tree -> Branch( "cosAl",		&cosAl );

    tree -> Branch( "c_mu1_pt",   &c_mu1_pt );
    tree -> Branch( "c_mu1_eta",  &c_mu1_eta );
    tree -> Branch( "c_mu1_phi",  &c_mu1_phi );
    tree -> Branch( "c_mu2_pt",   &c_mu2_pt );
    tree -> Branch( "c_mu2_eta",  &c_mu2_eta );
    tree -> Branch( "c_mu2_phi",  &c_mu2_phi );
    tree -> Branch( "c_trk_pt",   &c_trk_pt );
    tree -> Branch( "c_trk_eta",  &c_trk_eta );
    tree -> Branch( "c_trk_phi",  &c_trk_phi );
    tree -> Branch( "c_J_mass", &c_J_mass );
    tree -> Branch( "c_J_pt",   &c_J_pt );
    tree -> Branch( "c_J_eta",  &c_J_eta );
    tree -> Branch( "c_J_phi",  &c_J_phi );
    tree -> Branch( "c_mass",   &c_mass );
    tree -> Branch( "c_pt",     &c_pt );
    tree -> Branch( "c_eta",    &c_eta );
    tree -> Branch( "c_phi",    &c_phi );

  }
  if( isMC_ ) {
    tree_gen = fs->make<TTree>("gen_ntuple","B+ -> JPsiK+ Gen");
    tree_gen -> Branch( "gen_nB", &gen_nB );
    tree_gen -> Branch( "gen_trks", &gen_trks );
    tree_gen -> Branch( "gen_B_pt",		&gen_B_pt );                          
    tree_gen -> Branch( "gen_B_energy",		&gen_B_energy );         
    tree_gen -> Branch( "gen_B_eta",		&gen_B_eta );            
    tree_gen -> Branch( "gen_B_phi",		&gen_B_phi );            
    tree_gen -> Branch( "gen_B_pz",		&gen_B_pz );             
    tree_gen -> Branch( "gen_B_pid",		&gen_B_pid );          
    tree_gen -> Branch( "gen_BMuM_pt",		&gen_BMuM_pt );          
    tree_gen -> Branch( "gen_BMuM_eta",		&gen_BMuM_eta );         
    tree_gen -> Branch( "gen_BMuM_phi",		&gen_BMuM_phi );         
    tree_gen -> Branch( "gen_BMuM_pid",		&gen_BMuM_pid );         
    tree_gen -> Branch( "gen_BMuP_pt",		&gen_BMuP_pt );          
    tree_gen -> Branch( "gen_BMuP_eta",		&gen_BMuP_eta );         
    tree_gen -> Branch( "gen_BMuP_phi",		&gen_BMuP_phi );         
    tree_gen -> Branch( "gen_BMuP_pid",		&gen_BMuP_pid );         
    tree_gen -> Branch( "gen_BK_pt",		&gen_BK_pt );            
    tree_gen -> Branch( "gen_BK_energy",	&gen_BK_energy );        
    tree_gen -> Branch( "gen_BK_eta",		&gen_BK_eta );           
    tree_gen -> Branch( "gen_BK_phi",		&gen_BK_phi );           
    tree_gen -> Branch( "gen_BK_pid",		&gen_BK_pid );           
    tree_gen -> Branch( "gen_BJpsi_pt",		&gen_BJpsi_pt );         
    tree_gen -> Branch( "gen_BJpsi_energy",	&gen_BJpsi_energy );     
    tree_gen -> Branch( "gen_BJpsi_eta",	&gen_BJpsi_eta );        
    tree_gen -> Branch( "gen_BJpsi_phi",	&gen_BJpsi_phi );            
    tree_gen -> Branch( "gen_BJpsi_pid",	&gen_BJpsi_pid );            
    tree_gen -> Branch( "gen_otrk_pt",  &gen_otrk_pt );
    tree_gen -> Branch( "gen_otrk_eta", &gen_otrk_eta );
    tree_gen -> Branch( "gen_otrk_phi", &gen_otrk_phi );
    tree_gen -> Branch( "gen_otrk_pid", &gen_otrk_pid );
  }  
}

void Bu2JPsiK::endJob(){
  if( !OnlyGen_ ){
    tree -> GetDirectory() -> cd();
    tree -> Write();
  }
  if( isMC_ ){
    tree_gen -> GetDirectory() -> cd();
    tree_gen -> Write();
  }
}
void Bu2JPsiK::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
void Bu2JPsiK::fillGenBranches( const edm::Event& iEvent){
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticlesCollection_, genParticles);

  //static int xi1 =0; 
  //int xi1=0, xi2 = 0;
  //std::cout  << xi1++ << "\t";
  //std::cout<<__LINE__<<"\t" <<  genParticles -> size() << std::endl; 
  for(auto p = genParticles -> begin(); p != genParticles -> end(); ++p) {
    if( fabs(p -> pdgId()) == 521 ){ // B+
      //xi1++;
      //if( p -> mother() -> pdgId() == 523 || p -> mother() -> pdgId() == 521)
      //std::cout << xi2++ << "\t" << p-> mother() -> pdgId() << "\t" << p -> mother() -> mother() -> pdgId() << "\n";
      if( p -> numberOfDaughters() != 2 ) continue;
      if( fabs(p -> daughter(0) -> pdgId()) != 443 ) continue; // not Jpsi
      if( fabs(p -> daughter(1) -> pdgId()) != 321 ) continue; // not K+
      if( fabs(p -> daughter(0) -> daughter(0) -> pdgId()) != 13 ) continue; // not Mu-
      if( fabs(p -> daughter(0) -> daughter(1) -> pdgId()) != 13 ) continue; // not Mu+
      if( p -> daughter(0) -> daughter(0) -> pt() < 6 && p -> daughter(0) -> daughter(1) -> pt() < 6 ) continue;
      gen_nB++;
      gen_BMuM_pt      -> push_back( p -> daughter(0) -> daughter(0) -> pt() );
      gen_BMuM_eta     -> push_back( p -> daughter(0) -> daughter(0) -> eta() );
      gen_BMuM_phi     -> push_back( p -> daughter(0) -> daughter(0) -> phi() );
      gen_BMuM_pid     -> push_back( p -> daughter(0) -> daughter(0) -> pdgId() );
      gen_BMuP_pt      -> push_back( p -> daughter(0) -> daughter(1) -> pt() );
      gen_BMuP_eta     -> push_back( p -> daughter(0) -> daughter(1) -> eta() );
      gen_BMuP_phi     -> push_back( p -> daughter(0) -> daughter(1) -> phi() );
      gen_BMuP_pid     -> push_back( p -> daughter(0) -> daughter(1) -> pdgId() );
      gen_BK_pt        -> push_back( p -> daughter(1) -> pt() );
      gen_BK_energy    -> push_back( p -> daughter(1) -> energy() );
      gen_BK_eta       -> push_back( p -> daughter(1) -> eta() );
      gen_BK_phi       -> push_back( p -> daughter(1) -> phi() );
      gen_BK_pid       -> push_back( p -> daughter(1) -> pdgId() );
      gen_BJpsi_pt     -> push_back( p -> daughter(0) -> pt() );
      gen_BJpsi_energy -> push_back( p -> daughter(0) -> energy() );
      gen_BJpsi_eta    -> push_back( p -> daughter(0) -> eta() );
      gen_BJpsi_phi    -> push_back( p -> daughter(0) -> phi() );
      gen_BJpsi_pid    -> push_back( p -> daughter(0) -> pdgId() );
      gen_B_pt         -> push_back( p -> pt());
      gen_B_energy     -> push_back( p -> energy());
      gen_B_eta        -> push_back( p -> eta());
      gen_B_phi        -> push_back( p -> phi());
      gen_B_pz         -> push_back( p -> pz());
      gen_B_pid        -> push_back( p -> pdgId());
      //genKeta = p -> daughter(1) -> eta(); genKphi = p -> daughter(1) -> phi();
      //std::cout << "(" << genKeta << ", " << genKphi << ")\t";
    }
  }
  if(gen_nB > 0 ){
    for(auto p = genParticles -> begin(); p != genParticles -> end(); ++p) {
      if( fabs(p -> pdgId()) == 321 || fabs(p -> pdgId()) == 211 ){
	float deleta = fabs(p -> eta() - gen_BK_eta -> at(gen_nB -1 ));
	float delphi = fabs(p -> phi() - gen_BK_phi -> at(gen_nB -1 ));
	if(delphi > M_PI)delphi = 2*M_PI - delphi;
	float delr = sqrt(deleta*deleta+delphi*delphi);
	if(delr < 0.3  && delr > 1e-3 ){
	  gen_otrk_pt  -> push_back( p -> pt() );
	  gen_otrk_eta -> push_back( p -> eta() );
	  gen_otrk_phi -> push_back( p -> phi() );
	  gen_otrk_pid -> push_back( p -> pdgId() );
	  gen_trks++;
	}
      }
    }
  }
  //std::cout<<__LINE__<<"\n"; 
  if( gen_nB > 0 ) tree_gen -> Fill();
  //if( xi1 > 1 ) std::cout << xi1 << " <-- number of B+\n";
  //if( gen_nB > 1 ) std::cout << "# of B: " << gen_nB << ", CT: " << gen_trks << std::endl;
  gen_nB = 0; gen_trks = 0;
  gen_B_pt -> clear();
  gen_B_energy -> clear();
  gen_B_eta -> clear();
  gen_B_phi -> clear();
  gen_B_pz -> clear();
  gen_B_pid -> clear();
  gen_BMuM_pt -> clear();
  gen_BMuM_eta -> clear();
  gen_BMuM_phi -> clear();
  gen_BMuM_pid -> clear();
  gen_BMuP_pt -> clear();
  gen_BMuP_eta -> clear();
  gen_BMuP_phi -> clear();
  gen_BMuP_pid -> clear();
  gen_BK_pt -> clear();
  gen_BK_energy -> clear();
  gen_BK_eta -> clear();
  gen_BK_phi -> clear();
  gen_BK_pid -> clear();
  gen_BJpsi_pt -> clear();
  gen_BJpsi_energy -> clear();
  gen_BJpsi_eta -> clear();
  gen_BJpsi_phi -> clear();
  gen_BJpsi_pid -> clear();
  gen_otrk_pt  -> clear();
  gen_otrk_eta -> clear();
  gen_otrk_phi -> clear();
  gen_otrk_pid -> clear();
}

DEFINE_FWK_MODULE(Bu2JPsiK);
