#ifndef __SVTXPURITYSTUDY_H__
#define __SVTXPURITYSTUDY_H__

#include <fun4all/SubsysReco.h>
#include <string>
#include <TH1D.h>
#include <TH2D.h>

class PHCompositeNode;

class SvtxPurityStudy: public SubsysReco {

public: 

  SvtxPurityStudy(const std::string &name="SvtxPurityStudy");

  int Init(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  void set_nlayers(const unsigned int &nlayers) {_nlayers = nlayers;}
	void Detector(const std::string &d) {_detector = d;}

private:

  // event counter
  unsigned long long _event;

  // number of layers
  unsigned int _nlayers;

	// which calorimeter are we looking at
	std::string _detector;

  // output histograms ---------------------------------------------------------
  
  TH1D* _truept_particles_leaving7Hits;       // pattern reco eff baseline
  TH1D* _truept_track_found;                // cut eff denominator
	TH2D* _cal_e_over_recop_pt;
  TH2D* _truept_ecut;                  // e/p cut eff baseline
  TH2D* _recopt_tracks_all_ecut;                   // purity 
  TH2D* _recopt_tracks_recoWithExactHits_ecut;     // purity by nhit match
  TH2D* _recopt_tracks_recoWithin1Hit_ecut;
  TH2D* _recopt_tracks_recoWithin2Hits_ecut;

};

#endif // __SVTXPURITYSTUDY_H__
