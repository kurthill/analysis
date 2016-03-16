
#include "SvtxPurityStudy.h"

#include <phool/getClass.h>
#include <fun4all/Fun4AllServer.h>

#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4VtxPoint.h>

#include <g4hough/SvtxVertexMap.h>
#include <g4hough/SvtxVertex.h>
#include <g4hough/SvtxTrackMap.h>
#include <g4hough/SvtxTrack.h>

#include <g4eval/SvtxEvalStack.h>
#include <g4eval/SvtxTrackEval.h>
#include <g4eval/SvtxVertexEval.h>
#include <g4eval/SvtxTruthEval.h>

#include <TH1D.h>
#include <TH2D.h>

#include <math.h>
#include <iostream>

using namespace std;

SvtxPurityStudy::SvtxPurityStudy(const string &name)
  : SubsysReco(name),
    _event(0),
    _nlayers(7),
		_detector("HCALIN"),
    _truept_particles_leaving7Hits(NULL),
		_truept_track_found(NULL),
		_cal_e_over_recop_pt(NULL),
    _truept_ecut(NULL),
    _recopt_tracks_all_ecut(NULL),
    _recopt_tracks_recoWithExactHits_ecut(NULL),
    _recopt_tracks_recoWithin1Hit_ecut(NULL),
    _recopt_tracks_recoWithin2Hits_ecut(NULL) {
}

int SvtxPurityStudy::Init(PHCompositeNode *topNode) {

  // register histograms
  Fun4AllServer *se = Fun4AllServer::instance();

  _truept_particles_leaving7Hits = new TH1D("truept_particles_leaving7Hits",
					    "truept_particles_leaving7Hits", 20,0.0,10.0);

  _truept_track_found = new TH1D("truept_track_found",
					    "truept_track_found", 20,0.0,10.0);

  _cal_e_over_recop_pt = new TH2D("cal_e_over_recop_pt",
			"cal_e_over_recop_pt",100,0.0,10.0,12,0,12);

  _truept_ecut = new TH2D("truept_ecut",
			"truept_ecut",20,0.0,10.0,20,0.5,20.5);

  _recopt_tracks_all_ecut = new TH2D("recopt_tracks_all_ecut",
			"recopt_tracks_all_ecut",20,0.0,10.0,20,0.5,20.5);

  _recopt_tracks_recoWithExactHits_ecut = new TH2D("recopt_tracks_recoWithExactHits_ecut",
			"recopt_tracks_recoWithExactHits_ecut",20,0.0,10.0,20,0.5,20.5);

  _recopt_tracks_recoWithin1Hit_ecut = new TH2D("recopt_tracks_recoWithin1Hit_ecut",
			"recopt_tracks_recoWithin1Hit_10",20,0.0,10.0,20,0.5,20.5);
  
  _recopt_tracks_recoWithin2Hits_ecut = new TH2D("recopt_tracks_recoWithin2Hits_ecut",
			"recopt_tracks_recoWithin2Hits_10",20,0.0,10.0,20,0.5,20.5);

 
  se->registerHisto(_truept_particles_leaving7Hits);       
  se->registerHisto(_truept_track_found);       
  se->registerHisto(_cal_e_over_recop_pt);
  se->registerHisto(_truept_ecut);
  se->registerHisto(_recopt_tracks_all_ecut);
  se->registerHisto(_recopt_tracks_recoWithExactHits_ecut);
  se->registerHisto(_recopt_tracks_recoWithin1Hit_ecut);
  se->registerHisto(_recopt_tracks_recoWithin2Hits_ecut);

  return 0;
}

int SvtxPurityStudy::process_event(PHCompositeNode *topNode) {

  ++_event;
  
  // need things off of the DST...
  PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
  if (!truthinfo) {
    cerr << PHWHERE << " ERROR: Can't find G4TruthInfo" << endl;
    exit(-1);
  }

  SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode,"SvtxTrackMap");
  if (!trackmap) {
    cerr << PHWHERE << " ERROR: Can't find SvtxTrackMap" << endl;
    exit(-1);
  }

  SvtxVertexMap* vertexmap = findNode::getClass<SvtxVertexMap>(topNode,"SvtxVertexMap");
  if (!vertexmap) {
    cerr << PHWHERE << " ERROR: Can't find SvtxVertexMap" << endl;
    exit(-1);
  }

  // create SVTX eval stack
  SvtxEvalStack svtxevalstack(topNode);

  //SvtxVertexEval*   vertexeval = svtxevalstack.get_vertex_eval();
  SvtxTrackEval*     trackeval = svtxevalstack.get_track_eval();
  SvtxTruthEval*     trutheval = svtxevalstack.get_truth_eval();
  
  // loop over all truth particles
  PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange();
  for (PHG4TruthInfoContainer::ConstIterator iter = range.first; 
       iter != range.second; 
       ++iter) {
    PHG4Particle* g4particle = iter->second;
    
		if (trutheval->get_embed(g4particle) != 0) 
			continue;
    
    std::set<PHG4Hit*> g4hits = trutheval->all_truth_hits(g4particle);     
    float ng4hits = g4hits.size();  

    float truept = sqrt(pow(g4particle->get_px(),2)+pow(g4particle->get_py(),2));
		//cout << "NEW TRUTH TRACK" << endl;
    
    // examine truth particles that leave 7 detector hits
		//   (acceptance condition)
    if (ng4hits == _nlayers) {
      _truept_particles_leaving7Hits->Fill(truept);
    
      SvtxTrack* track = trackeval->best_track_from(g4particle);
    
      if (!track) {continue;}
      
      unsigned int nfromtruth = trackeval->get_nclusters_contribution(track,g4particle);
      unsigned int ndiff = abs((int)nfromtruth-(int)_nlayers);

			// Must be exact match, for now . . .
      if (ndiff != 0)
				continue;

      float recop = track->get_p();
			float cal_e;
			if (_detector == "CEMC")
				cal_e = track->get_cal_energy_3x3(SvtxTrack::CEMC);
			else if (_detector == "HCALIN")
				cal_e = track->get_cal_energy_3x3(SvtxTrack::HCALIN);
			else
			{
    		cerr << PHWHERE << " ERROR: Detector " << _detector 
						 << " not supported" << endl;
    		exit(-1);
			}

			// if there was no track projection cal_e will be NaN
			if ( isnan(cal_e) )
				continue;

			_cal_e_over_recop_pt->Fill( cal_e/recop, truept );

		//cout << " RECO P = " << recop << ", EMC E = " << cal_e << ", E/PT = " << cal_e/recop << endl;

			// Fill histograms for cut efficiency
			//
      _truept_track_found->Fill(truept);
			float low = 0.;
			int cutbin = 0;
			while (cutbin < 21)
			{
				cutbin++;
				low += 0.1;
				if (cal_e/recop > low)
					_truept_ecut->Fill(truept,cutbin);
			}

    }    
  }

  // loop over all reco particles
  for (SvtxTrackMap::Iter iter = trackmap->begin();
       iter != trackmap->end();
       ++iter) {
    
    SvtxTrack*    track      = iter->second;
      
    float recopt = track->get_pt();
    float recop = track->get_p();
		float cal_e;
		if (_detector == "CEMC")
			cal_e = track->get_cal_energy_3x3(SvtxTrack::CEMC);
		else if (_detector == "HCALIN")
			cal_e = track->get_cal_energy_3x3(SvtxTrack::HCALIN);
		else
		{
    	cerr << PHWHERE << " ERROR: Detector " << _detector 
					 << " not supported" << endl;
    	exit(-1);
		}

		// if there was no track projection cal_e will be NaN
		if ( isnan(cal_e) )
			continue;

		//cout << "NEW TRACK" << endl;
		//cout << " RECO P = " << recop << ", EMC E = " << cal_e << ", E/PT = " << cal_e/recop << endl;
        
    PHG4Particle* g4particle = trackeval->max_truth_particle_by_nclusters(track);
    if (!g4particle) {continue;}

    //float truept = sqrt(pow(g4particle->get_px(),2)+pow(g4particle->get_py(),2));

		// non-embedded results (purity measures)
    if (trutheval->get_embed(g4particle) != 0) 
			continue;

    unsigned int nfromtruth = trackeval->get_nclusters_contribution(track,g4particle);

    unsigned int ndiff = abs((int)nfromtruth-(int)_nlayers);

		// Fill histograms for purity
		//
		float low = 0.;
		int cutbin = 0;
		while (cutbin < 21)
		{
			cutbin++;
			low += 0.1;
			if (cal_e/recop > low)
			{
				_recopt_tracks_all_ecut->Fill(recopt,cutbin);
				if (ndiff <= 2)
					_recopt_tracks_recoWithin2Hits_ecut->Fill(recopt,cutbin);
				if (ndiff <= 1)
					_recopt_tracks_recoWithin1Hit_ecut->Fill(recopt,cutbin);
				if (ndiff <= 0)
					_recopt_tracks_recoWithExactHits_ecut->Fill(recopt,cutbin);
			}
		}
    
  }

  return 0;
}

int SvtxPurityStudy::End(PHCompositeNode *topNode) {
 return 0;
}
