#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "MagneticField/ParametrizedEngine/src/OAEParametrizedMagneticField.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

// RM: dirrrrrrty
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"


class KVFitter {

  public:
    KVFitter() {};
    virtual ~KVFitter() {};

    // constructed from reco::TrackRef
    reco::TransientTrack getTransientTrack(const reco::TrackRef& trackRef) {    
      reco::TransientTrack transientTrack(trackRef, paramField);
      return transientTrack;
    }

    // constructed from reco::Track
    reco::TransientTrack getTransientTrack(const reco::Track& track) {    
      reco::TransientTrack transientTrack(track, paramField);
      return transientTrack;
    }

    // constructed from reco::TrackRef
    TransientVertex Fit(const std::vector<reco::TrackRef> & tracks)
    {
      // do tau vertex fit
      std::vector<reco::TransientTrack> tks;
      for (std::vector<reco::TrackRef>::const_iterator itk = tracks.begin(); itk != tracks.end(); ++itk){
          tks.push_back(getTransientTrack(*itk));
      }
    
      KalmanVertexFitter kvf;
      TransientVertex tv = kvf.vertex(tks);
            
      return tv;    
    };

    // constructed from reco::Track
    TransientVertex Fit(const std::vector<reco::Track> & tracks)
    {
      // do tau vertex fit
      std::vector<reco::TransientTrack> tks;
      for (std::vector<reco::Track>::const_iterator itk = tracks.begin(); itk != tracks.end(); ++itk){
          tks.push_back(getTransientTrack(*itk));
      }
    
      KalmanVertexFitter kvf;
      TransientVertex tv = kvf.vertex(tks);
      
      return tv;          
    };

  private:
    OAEParametrizedMagneticField *paramField = new OAEParametrizedMagneticField("3_8T");

};
