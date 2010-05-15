// my includes
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/TrackReco/interface/TrackExtraBase.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "../interface/TrackerIsolation.h"

//CLHEP
//#include <CLHEP/Vector/LorentzVector.h>

#include <TLorentzVector.h>
#include <TVector3.h>

using namespace std;
TrackerIsolation::TrackerIsolation (){}

TrackerIsolation::TrackerIsolation (const Track *eleTrack, const TrackCollection *trackColl) : 
  _eleTrack(eleTrack),
  _tracks(trackColl)   
{
  _extRadius = 0.20;
  _intRadius = 0.015;
}

TrackerIsolation::~TrackerIsolation (){}

void TrackerIsolation::setExtRadius (float extRadius){_extRadius = extRadius; }

void TrackerIsolation::setIntRadius (float intRadius){_intRadius = intRadius; }

float TrackerIsolation::getPtTracks (bool relative) const {
  
  float dummyPt = 0 ;
  
  //  Hep3Vector elePAtVtx(_eleTrack->px(), _eleTrack->py(), _eleTrack->pz()); 
  TVector3 elePAtVtx;
  elePAtVtx.SetXYZ(_eleTrack->px(), _eleTrack->py(), _eleTrack->pz()); 
  float ele_pt = _eleTrack->pt();

  TrackCollection::const_iterator this_track;
  for(this_track = _tracks->begin(); this_track != _tracks->end(); this_track++ ){ 
    
    //    Hep3Vector trackPAtVtx(this_track->px(),this_track->py(),this_track->pz());
    TVector3  trackPAtVtx; 
    trackPAtVtx.SetXYZ(this_track->px(),this_track->py(),this_track->pz());
    float this_pt  = trackPAtVtx.Perp();

    // usually low pt tracks are fakes
    if ( this_pt < 1.0 ) continue;

    double dr = elePAtVtx.DeltaR(trackPAtVtx);
    if ( fabs(dr) < _extRadius && fabs(dr) > _intRadius ){ 
      if ( relative ) dummyPt += this_pt/ele_pt;
      else dummyPt += this_pt;
    } 
    
  } //end loop over tracks		       
  
  return dummyPt;
}

bool TrackerIsolation::isIsolated (float ptCut) const {

  bool dummyIsolation = true ;  
  if (TrackerIsolation::getPtTracks() > ptCut ) dummyIsolation = false ;  
  
  return dummyIsolation ;
}
