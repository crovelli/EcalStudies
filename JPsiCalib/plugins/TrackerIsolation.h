#ifndef TrackerIsolation_h
#define TrackerIsolation_h

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/TrackReco/interface/TrackExtraBase.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

using namespace edm;
using namespace std;
using namespace reco;

class TrackerIsolation{
 public:
  
  //constructors
  TrackerIsolation();
  TrackerIsolation(const Track *eleTrack, const TrackCollection *trackColl);

  //methods
  void setExtRadius (float extRadius);
  void setIntRadius (float intRadius);
  float getPtTracks (bool relative=true) const;
  bool isIsolated (float ptCut = 0.05) const;

  //destructor 
  ~TrackerIsolation();
  
 private:

  const Track *_eleTrack;  	  
  const TrackCollection *_tracks;
  
  float _extRadius;
  float _intRadius;
};

#endif
