#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

class Utilities {
 private:
  
 public:
  float computeTrkMuonIsolation(const pat::Muon& muon, 
				 const pat::Electron& the_electron,
				 //edm::Handle<edm::View<pat::PackedCandidate>> pfCandHandle_,
				 std::vector<pat::PackedCandidate> pfCandHandle_,
				 unsigned int primaryVertexIndex,
				 float minPt=0.5, float dR=0.5,
				 std::vector<const pat::PackedCandidate*> ignoreTracks = 
				 std::vector<const pat::PackedCandidate*>());
  float computeTrkEleIsolation(const pat::Electron& the_electron, const pat::Muon& muon, std::vector<pat::PackedCandidate> pfCandHandle_,
			       unsigned int primaryVertexIndex, float minPt=0.5, float dR=0.5,
			       std::vector<const pat::PackedCandidate*> ignoreTracks = std::vector<const pat::PackedCandidate*>());
};

// https://github.com/drkovalskyi/Bmm5/blob/origin/RunII-NanoAODv8/NanoAOD/plugins/BxToMuMuProducer.cc#L861-L887 
float Utilities::computeTrkMuonIsolation(const pat::Muon& the_muon, 
					 const pat::Electron& the_electron, 
					 std::vector<pat::PackedCandidate> pfCandidates,
					 unsigned int primaryVertexIndex,
					 float minPt, float dR,
					 std::vector<const pat::PackedCandidate*> ignoreTracks){
  float sumPt(0);
  for (const auto& pfCand: pfCandidates){
    bool ignore_track = false;
    for (auto trk: ignoreTracks){
      if (trk==&pfCand){
	ignore_track = true;
	break;
      }
    }
    if (ignore_track) continue;
    if (deltaR(the_muon, pfCand) < 0.01 || deltaR(the_electron, pfCand) < 0.01) continue;
    if (pfCand.charge() == 0 ) continue;
    if (!pfCand.hasTrackDetails()) continue;
    if (pfCand.pt()<minPt) continue;
    if (pfCand.vertexRef().key()!=primaryVertexIndex) continue;
    if (deltaR(the_muon, pfCand) > dR) continue;
    sumPt += pfCand.pt();
  }
  return the_muon.pt()/(the_muon.pt()+sumPt);
}

float Utilities::computeTrkEleIsolation(const pat::Electron& the_electron, const pat::Muon& the_muon, 
					std::vector<pat::PackedCandidate> pfCandidates,
					unsigned int primaryVertexIndex,
					float minPt, float dR,
					std::vector<const pat::PackedCandidate*> ignoreTracks){
  float sumPt(0);
  for (const auto& pfCand: pfCandidates){
    bool ignore_track = false;
    for (auto trk: ignoreTracks){
      if (trk==&pfCand){
	ignore_track = true;
	break;
      }
    }
    if (ignore_track) continue;
    if (deltaR(the_muon, pfCand) < 0.01 || deltaR(the_electron, pfCand) < 0.01) continue;
    if (pfCand.charge() == 0 ) continue;
    if (!pfCand.hasTrackDetails()) continue;
    if (pfCand.pt()<minPt) continue;
    if (pfCand.vertexRef().key()!=primaryVertexIndex) continue;
    if (deltaR(the_electron, pfCand) > dR) continue;
    sumPt += pfCand.pt();
  }
  return the_electron.pt()/(the_electron.pt()+sumPt);
}
