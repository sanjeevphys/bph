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

#include "DataFormats/PatCandidates/interface/PackedCandidate.h" // muy importante para MiniAOD

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//#include "RecoEgamma/EgammaElectronAlgos/interface/GsfElectronAlgo.h"
//#include "DataFormats/ParticleFlowReco/interface/GsfPFRecTrack.h"
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
  unsigned int nB, nn;

  int run,lumi,event, npv,ncands, hlt_path;
  std::vector<float> *t_pt, *t_eta, *t_phi;
  std::vector<float> *mass, *pt, *eta, *rap, *phi;
  std::vector<float> *dr_mu, *dr_ele;
  std::vector<float> *pv_x, *pv_y, *pv_z; float bs_x0, bs_y0, bs_z0; 
  std::vector<float> *bs_x, *bs_y, *vx, *vy, *vz, *vtx_chi2, *vtx_prob, *cos2d, *lxy, *lxy_err, *lxy_sig;
  std::vector<float> *vx_ctf, *vy_ctf, *vz_ctf, *vtx_chi2_ctf, *vtx_prob_ctf, *cos2d_ctf, *lxy_ctf, *lxy_err_ctf, *lxy_sig_ctf;
  std::vector<int> *pfmu, *globalMu, *mu_charge, *ele_charge;
  std::vector<float> *mu_pt, *mu_eta, *mu_rap, *mu_phi, *mu_e, *mu_mass, *mu_id_loose, *mu_id_soft, *mu_id_medium, *mu_id_tight, *mu_id_soft_mva, *mu_dxy, *mu_dxy_e, *mu_dxy_sig, *mu_dz, *mu_dz_e, *mu_dz_sig, *mu_bs_dxy, *mu_bs_dxy_e, *mu_bs_dxy_sig, *mu_cov_pos_def, *mu_trkIsolation;
  std::vector<float> *ele_pt, *ele_eta, *ele_rap, *ele_phi, *ele_e, *ele_mass, *ele_id_iso_wp90, *ele_id_iso_wp80, *ele_id_iso_wpLoose, *ele_id_noIso_wp90, *ele_id_noIso_wp80, *ele_id_noIso_wpLoose, *ele_dxy, *ele_dxy_e, *ele_dxy_sig, *ele_dz, *ele_dz_e, *ele_dz_sig, *ele_bs_dxy, *ele_bs_dxy_e, *ele_bs_dxy_sig, *ele_cov_pos_def, *ele_trkIsolation;
};

BsToEMu::BsToEMu(const edm::ParameterSet& iConfig):
  electronToken_(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electrons"))),
  muonToken_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
  beamSpotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
  pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
  vertices_ (consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  ltrk_ (consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("ltrk"))),
  trg_ (consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("trg_res"))),
  tobjs_ (consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("tobjs"))),
  
  isMC_(iConfig.getParameter<bool>("isMC")),
  OnlyGen_(iConfig.getParameter<bool>("OnlyGen")),
  //now do what ever initialization is needed
  tree(0),
  nB(0), nn(0),
  run(0), lumi(0), event(0), npv(0), ncands(0), hlt_path(0), 
  t_pt(0), t_eta(0), t_phi(0),
  mass(0), pt(0), eta(0), rap(0), phi(0), dr_mu(0), dr_ele(0), pv_x(0), pv_y(0), pv_z(0), bs_x0(0), bs_y0(0), bs_z0(0), bs_x(0), bs_y(0), vx(0), vy(0), vz(0),
  vtx_chi2(0), vtx_prob(0), cos2d(0), lxy(0), lxy_err(0), lxy_sig(0),
  vx_ctf(0), vy_ctf(0), vz_ctf(0), vtx_chi2_ctf(0), vtx_prob_ctf(0), cos2d_ctf(0), lxy_ctf(0), lxy_err_ctf(0), lxy_sig_ctf(0),
  pfmu(0), globalMu(0), mu_charge(0), ele_charge(0),
  mu_pt(0), mu_eta(0), mu_rap(0), mu_phi(0), mu_e(0), mu_mass(0), mu_id_loose(0), mu_id_soft(0), mu_id_medium(0), mu_id_tight(0), mu_id_soft_mva(0), mu_dxy(0), mu_dxy_e(0), mu_dxy_sig(0), mu_dz(0), mu_dz_e(0), mu_dz_sig(0), mu_bs_dxy(0), mu_bs_dxy_e(0), mu_bs_dxy_sig(0), mu_cov_pos_def(0), mu_trkIsolation(0), ele_pt(0), ele_eta(0), ele_rap(0), ele_phi(0), ele_e(0), ele_mass(0), ele_id_iso_wp90(0), ele_id_iso_wp80(0), ele_id_iso_wpLoose(0), ele_id_noIso_wp90(0), ele_id_noIso_wp80(0), ele_id_noIso_wpLoose(0), ele_dxy(0), ele_dxy_e(0), ele_dxy_sig(0), ele_dz(0), ele_dz_e(0), ele_dz_sig(0), ele_bs_dxy(0), ele_bs_dxy_e(0), ele_bs_dxy_sig(0), ele_cov_pos_def(0), ele_trkIsolation(0)
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

  bs_x0 = recoBeamSpotHandle -> x0();
  bs_y0 = recoBeamSpotHandle -> y0();
  bs_z0 = recoBeamSpotHandle -> z0();

  //edm::ESHandle<TransientTrackBuilder> theB; 
  //iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB); 

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

  const char *var1;
  auto trg_names = iEvent.triggerNames(*trg_res);
  int hlt_passed=0;
  for(int ipath=0; ipath < nhlt; ++ipath){
    for(unsigned iname = 0; iname < trg_res->size(); ++iname){ //.triggerNames()
      std::string name = trg_names.triggerName(iname);
      var1 = name.c_str();
      if(strstr(var1, paths[ipath]) ){ // && (strlen(var1) - strlen(paths[ipath]) < 5) ){
	//cout << "This is not testing: ";
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
  nn++;
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
  std::vector <reco::Track> tracksTifit, ctfFit;
  for(View<pat::Muon>::const_iterator iMu = muons->begin(); iMu != muons->end(); ++iMu){
    for(View<pat::Electron>::const_iterator iE = electrons->begin(); iE != electrons->end(); ++iE){
      if( (iMu -> charge()) * (iE->charge()) == 1 ) continue;
      if( iMu  -> pt() < 7.0 || fabs(iMu -> eta()) > 2.5) continue;
      if( iE -> pt() < 4.0 || fabs(iE -> eta()) > 2.5) continue;

      if( ! iMu -> bestTrack()) continue;
      if( ! iE -> gsfTrack().isNonnull() ) continue;
      if( ! iE-> closestCtfTrackRef().isNonnull() ) continue;
       
      //cout << "============================ Loop has begin =========================================\n";
      tracksTifit.clear();
      ctfFit.clear();
      
      //cout << ".-.-.-.-.-.--.-.-.-.-.-.-.-.-. clear vertex fitter .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.- \n";
      tracksTifit.push_back( *(iMu->bestTrack()) );
      tracksTifit.push_back( *(iE->gsfTrack() ) );

      ctfFit.push_back( *(iMu->bestTrack()) );
      ctfFit.push_back( *(iE-> closestCtfTrackRef()) );
      
      //cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~  tracks are passed ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
      auto ctfVtxFit = vtxfit.Fit( ctfFit );
      auto emuVtxFit = vtxfit.Fit( tracksTifit );
      //cout << "Validity of vtx fit: " << (int) ctfVtxFit.isValid() << (int) emuVtxFit.isValid() << endl;

      if( ! ctfVtxFit.isValid()) continue;
      if( ! emuVtxFit.isValid()) continue;
      
      //cout << "###########################  validity checekd ###############################\n";
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
      
      //pat::CompositeCandidate EMuCand;
      //EMuCand.addDaughter(*iMu, "iMu");
      //EMuCand.addDaughter(*iE, "iE");
      
      //reco::TransientTrack MuTT((*theB).build(iMu->track()));
      //reco::TransientTrack ETT((*theB).build(iE->track()));

      ParticleMass muon_mass = 0.10565837; //pdg mass
      ParticleMass el_mass = 0.0005109989461; //pdg mass

      TLorentzVector mu4, ele4, bcand;
      mu4.SetXYZM(iMu->px(), iMu->py(), iMu->pz(), muon_mass);
      ele4.SetXYZM(iE->px(), iE->py(), iE->pz(), el_mass);
      bcand=mu4+ele4;

      auto bs_point = Vertex::Point( recoBeamSpotHandle -> x(pv.z()), recoBeamSpotHandle -> y(pv.z()), recoBeamSpotHandle -> z0() );
      auto bs_error = recoBeamSpotHandle -> covariance3D();
      float chi2_temp=0, ndof_temp=0;
      auto bs_temp = Vertex(bs_point, bs_error, chi2_temp, ndof_temp, 3);
      
      TVector3 vect_lxy, vect_pt, vect_lxy_ctf;

      vect_lxy.SetXYZ(emuVtxFit.position().x() - bs_temp.position().x(),
		      emuVtxFit.position().y() - bs_temp.position().y(),
		      0. );
      vect_lxy_ctf.SetXYZ(ctfVtxFit.position().x() - bs_temp.position().x(),
			  ctfVtxFit.position().y() - bs_temp.position().y(),
			  0. );
      vect_pt.SetXYZ(bcand.Px(), bcand.Py(), 0. );

      double vtxCos = vect_pt.Dot(vect_lxy) / ( vect_pt.Mag() * vect_lxy.Mag() );
      double vtxCos_ctf = vect_pt.Dot(vect_lxy_ctf) / ( vect_pt.Mag() * vect_lxy_ctf.Mag() );
      if( vtxCos < 0.9 ) continue;

      if ( bcand.M()<4. || bcand.M()>7. ) continue;
      cos2d -> push_back( vtxCos );
      cos2d_ctf -> push_back( vtxCos_ctf );

      //if( iMu -> isPFMuon() ) { cout << "This is PF Muon\n"; } else { cout << "Not a PF Muon\n"; }
      //if( iMu -> isGlobalMuon() ) { cout << "This is Global Muon\n"; } else { cout << "Not a Global Muon\n"; }
      bs_x -> push_back(bs_temp.position().x());
      bs_y -> push_back(bs_temp.position().y());
      vx -> push_back(emuVtxFit.position().x());
      vy -> push_back(emuVtxFit.position().y());
      vz -> push_back(emuVtxFit.position().z());
      vtx_chi2 -> push_back(emuVtxFit.normalisedChiSquared());
      vtx_prob -> push_back(TMath::Prob(emuVtxFit.totalChiSquared(), emuVtxFit.degreesOfFreedom()));
      // TMath::Prob(fChisquare,fNpfits-fNpar); // https://root.cern.ch/root/roottalk/roottalk01/4648.html 

      //cout << "__________________________________________________________________________________________\n";
      //vx_ctf -> push_back(ctfVtxFit.position().x());
      //vy_ctf -> push_back(ctfVtxFit.position().y());
      //vz_ctf -> push_back(ctfVtxFit.position().z());
      //cout << "..............................................................................................\n";
      vtx_chi2_ctf -> push_back(ctfVtxFit.normalisedChiSquared());
      vtx_prob_ctf -> push_back(TMath::Prob(ctfVtxFit.totalChiSquared(), ctfVtxFit.degreesOfFreedom()));


      //cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
      auto llxy= VertexDistanceXY().distance(bs_temp, emuVtxFit.vertexState());
      lxy -> push_back(llxy.value());
      lxy_err -> push_back(llxy.error());
      lxy_sig -> push_back(llxy.significance());

      auto llxy_ctf= VertexDistanceXY().distance(bs_temp, ctfVtxFit.vertexState());
      lxy_ctf -> push_back(llxy_ctf.value());
      lxy_err_ctf -> push_back(llxy_ctf.error());
      lxy_sig_ctf -> push_back(llxy_ctf.significance());

      //cout << "-----------------------------------------------------------------------------------------\n";
      //delz
      //delz_sig

      pv_x -> push_back(pv.position().x());
      pv_y -> push_back(pv.position().y());
      pv_z -> push_back(pv.position().z());
      
      mass  -> push_back( bcand.M());
      pt    -> push_back( bcand.Pt());
      eta   -> push_back( bcand.Eta());
      rap   -> push_back( bcand.Rapidity());
      phi   -> push_back( bcand.Phi());
      dr_mu -> push_back( bcand.DeltaR(mu4));
      dr_ele-> push_back( bcand.DeltaR(ele4));

      pfmu -> push_back(iMu->isPFMuon());
      globalMu -> push_back(iMu->isGlobalMuon());

      mu_pt       -> push_back( iMu -> pt());
      mu_eta      -> push_back( iMu -> eta());
      mu_rap      -> push_back( iMu -> rapidity());
      mu_phi      -> push_back( iMu -> phi());
      /*mu_e        -> push_back( iMu -> energy());
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
      //mu_cov_pos_def-> push_back( is_pos_def( convert_cov( iMu -> bestTrack() -> covariance() ) ) );
      //cout << "Muon Track Isolation : " << utilityFuntions.computeTrkMuonIsolation( *iMu , *iE , *pfs , pvIndex ) << endl;
      mu_trkIsolation -> push_back( utilityFuntions.computeTrkMuonIsolation( *iMu , *iE , *pfs , pvIndex ) );*/
      //cout << iE->pfIsolationVariables() << endl;
	      
      ele_pt -> push_back( iE -> pt());
      ele_eta -> push_back( iE -> eta());
      ele_rap -> push_back( iE -> rapidity());
      ele_phi -> push_back( iE -> phi());
      /*ele_e -> push_back( iE -> energy());
      ele_mass -> push_back( iE -> mass());
      ele_charge -> push_back( iE -> charge());
      ele_id_iso_wp90 -> push_back( iE -> electronID("mvaEleID-Fall17-iso-V2-wp90"));
      ele_id_iso_wp80 -> push_back( iE -> electronID("mvaEleID-Fall17-iso-V2-wp80"));
      ele_id_iso_wpLoose -> push_back( iE -> electronID("mvaEleID-Fall17-iso-V2-wpLoose"));
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
      //ele_cov_pos_def -> push_back( iE -> ());
      ele_trkIsolation -> push_back( utilityFuntions.computeTrkEleIsolation( *iE, *iMu, *pfs, pvIndex ));*/

      //ncands,mu_cov_pos_def,ele_cov_pos_def
      nB++;
      //cout << ",,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, Loop is over for this count ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\n";
    }
  }
  //cout << ":::::::::::::::::::::::::::::: Just outside loop :::::::::::::::::::::::::::::::::::::::::::::\n";
  ncands = pt->size();
  if (nB > 0 ){
    //std::cout << "filling tree for " << nB << " inputs..." << endl;
    tree -> Fill();
    nB = 0;
  }
  
  t_pt     -> clear();
  t_eta    -> clear();
  t_phi    -> clear();

  mass  -> clear();
  pt    -> clear();
  eta   -> clear();
  rap   -> clear();
  phi   -> clear();
  dr_mu -> clear();
  dr_ele-> clear();

  pv_x -> clear(); pv_y -> clear(); pv_z -> clear();
  //bs_x0 -> clear();
  //bs_y0 -> clear();
  //bs_z0 -> clear();
  bs_x -> clear();
  bs_y -> clear();
  vx -> clear();
  vy -> clear();
  vz -> clear();
  vtx_chi2 -> clear();
  vtx_prob -> clear();
  cos2d -> clear();
  lxy -> clear();
  lxy_err -> clear();
  lxy_sig -> clear();

  //vx_ctf -> clear();
  //vy_ctf -> clear();
  //vz_ctf -> clear();
  vtx_chi2_ctf -> clear();
  vtx_prob_ctf -> clear();
  cos2d_ctf -> clear();
  lxy_ctf -> clear();
  lxy_err_ctf -> clear();
  lxy_sig_ctf -> clear();

  //delz -> clear();
  //delz_sig -> clear();

  pfmu -> clear();
  globalMu -> clear();
  mu_pt->clear();
  mu_eta      -> clear();
  mu_rap      -> clear();
  mu_phi      -> clear();
  /*mu_e        -> clear();
  mu_mass     -> clear();
  mu_charge   -> clear();
  mu_id_loose -> clear();
  mu_id_soft  -> clear();
  mu_id_medium-> clear();
  mu_id_tight -> clear();
  mu_id_soft_mva-> clear();
  mu_dxy      -> clear();
  mu_dxy_e    -> clear();
  mu_dxy_sig  -> clear();
  mu_dz       -> clear();
  mu_dz_e     -> clear();
  mu_dz_sig   -> clear();
  mu_bs_dxy   -> clear();
  mu_bs_dxy_e -> clear();
  mu_bs_dxy_sig -> clear();
  mu_cov_pos_def -> clear();
  mu_trkIsolation -> clear();*/

  ele_pt -> clear();
  ele_eta -> clear();
  ele_rap -> clear();
  ele_phi -> clear();
  /*ele_e -> clear();
  ele_mass -> clear();
  ele_charge -> clear();
  ele_id_iso_wp90 -> clear();
  ele_id_iso_wp80 -> clear();
  ele_id_iso_wpLoose -> clear();
  ele_id_noIso_wp90 -> clear();
  ele_id_noIso_wp80 -> clear();
  ele_id_noIso_wpLoose -> clear();
  ele_dxy -> clear();
  ele_dxy_e -> clear();
  ele_dxy_sig -> clear();
  ele_dz -> clear();
  ele_dz_e -> clear();
  ele_dz_sig -> clear();
  ele_bs_dxy -> clear();
  ele_bs_dxy_e -> clear();
  ele_bs_dxy_sig -> clear();
  ele_cov_pos_def -> clear();
  ele_trkIsolation -> clear();*/
  //#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  //ESHandle<SetupData> pSetup;
  //iSetup.get<SetupRecord>().get(pSetup);
  //#endif
  //cout << "******************************************* At the END **********************************\n";
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
  tree -> Branch("t_pt", &t_pt);
  tree -> Branch("t_eta", &t_eta);
  tree -> Branch("t_phi", &t_phi);

  tree -> Branch("mass", &mass );
  tree -> Branch("pt", &pt );
  tree -> Branch("eta", &eta );
  tree -> Branch("rap", &rap );
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

  //tree -> Branch("vy_ctf", &vy_ctf );
  //tree -> Branch("vz_ctf", &vz_ctf );
  tree -> Branch("vtx_chi2_ctf", &vtx_chi2_ctf );
  tree -> Branch("vtx_prob_ctf", &vtx_prob_ctf );
  tree -> Branch("cos2d_ctf", &cos2d_ctf );
  tree -> Branch("lxy_ctf", &lxy_ctf );
  tree -> Branch("lxy_err_ctf", &lxy_err_ctf );
  tree -> Branch("lxy_sig_ctf", &lxy_sig_ctf );

  tree -> Branch("pfmu", &pfmu);
  tree -> Branch("globalMu", &globalMu);

  tree -> Branch("mu_pt", &mu_pt );
  tree -> Branch("mu_eta", &mu_eta );
  tree -> Branch("mu_rap", &mu_rap );
  tree -> Branch("mu_phi", &mu_phi );
  /*tree -> Branch("mu_e", &mu_e );
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
  tree -> Branch("mu_trkIsolation", &mu_trkIsolation );*/

  tree -> Branch("ele_pt", &ele_pt );
  tree -> Branch("ele_eta", &ele_eta ); 
  tree -> Branch("ele_rap", &ele_rap );
  tree -> Branch("ele_phi", &ele_phi );
  /*tree -> Branch("ele_e", &ele_e );
  tree -> Branch("ele_mass", &ele_mass );
  tree -> Branch("ele_charge", &ele_charge );
  tree -> Branch("ele_id_iso_wp90", &ele_id_iso_wp90 );
  tree -> Branch("ele_id_iso_wp80", &ele_id_iso_wp80 );
  tree -> Branch("ele_id_iso_wpLoose", &ele_id_iso_wpLoose );
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
  tree -> Branch("ele_trkIsolation", &ele_trkIsolation );*/
}
void BsToEMu::endJob() {
  tree -> GetDirectory()->cd();
  tree -> Write();
}
void BsToEMu::fillDescriptions(edm::ConfigurationDescriptions& descriptions) { 
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(BsToEMu);
