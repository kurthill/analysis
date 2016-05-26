
#include "CaloResponse.h"

#include <phool/getClass.h>
#include <fun4all/Fun4AllServer.h>

#include <g4main/PHG4TruthInfoContainer.h>
#include <g4hough/SvtxTrack.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4VtxPoint.h>
#include <g4eval/SvtxEvalStack.h>
#include <g4eval/SvtxTrackEval.h>
#include <g4cemc/RawCluster.h>
#include <g4cemc/RawTowerGeomContainer.h>
#include <g4cemc/RawClusterContainer.h>

#include <TH1D.h>

#include <math.h>
#include <iostream>

using namespace std;

CaloResponse::CaloResponse(const string &name)
  : SubsysReco(name),
    _event(0),
    _h_calresponse_3x3(NULL),
    _h_calresponse_5x5(NULL),
    _h_calresponse_cluster(NULL),
    _h_calresponse_cluster_nocharge(NULL)
{
  _cal_names.push_back("CEMC");
  _cal_names.push_back("HCALIN");
  _cal_names.push_back("HCALOUT");
  _cal_types.push_back(SvtxTrack::CEMC);
  _cal_types.push_back(SvtxTrack::HCALIN);
  _cal_types.push_back(SvtxTrack::HCALOUT);
}

int CaloResponse::Init(PHCompositeNode *topNode) {

  // register histograms
  Fun4AllServer *se = Fun4AllServer::instance();

  _h_calresponse_3x3 = new TH1D(
      "h_calresponse_3x3","calresponse_3x3", 75,0.0,1.5);
  _h_calresponse_5x5 = new TH1D(
      "h_calresponse_5x5","calresponse_5x5", 75,0.0,1.5);
  _h_calresponse_cluster = new TH1D(
      "h_calresponse_cluster","calresponse_cluster", 75,0.0,1.5);

  _h_calresponse_cluster_nocharge = new TH1D(
      "h_calresponse_cluster_nocharge","calresponse_cluster_nocharge", 75,0.0,1.5);

  se->registerHisto(_h_calresponse_3x3);
  se->registerHisto(_h_calresponse_5x5);
  se->registerHisto(_h_calresponse_cluster);
  se->registerHisto(_h_calresponse_cluster_nocharge);

  return 0;
}

int CaloResponse::process_event(PHCompositeNode *topNode) {

  ++_event;
  
  // need things off of the DST...
  PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
  if (!truthinfo) {
    cerr << PHWHERE << " ERROR: Can't find G4TruthInfo" << endl;
    exit(-1);
  }

  // create SVTX eval stack
  SvtxEvalStack svtxevalstack(topNode);

  //SvtxVertexEval*   vertexeval = svtxevalstack.get_vertex_eval();
  SvtxTrackEval*     trackeval = svtxevalstack.get_track_eval();
  

  // loop over all truth particles
  PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange();
  for (PHG4TruthInfoContainer::ConstIterator iter = range.first; 
       iter != range.second; 
       ++iter) 
	{

    // Arrays to hold calorimeter energies
    double cal_e_3x3[_num_cal_layers];
    memset( cal_e_3x3, 0, _num_cal_layers*sizeof(double) );
    double cal_e_5x5[_num_cal_layers];
    memset( cal_e_5x5, 0, _num_cal_layers*sizeof(double) );
    double cal_e_clu[_num_cal_layers];
    memset( cal_e_clu, 0, _num_cal_layers*sizeof(double) );
    double cal_e_clu_nocharge[_num_cal_layers];
    memset( cal_e_clu_nocharge, 0, _num_cal_layers*sizeof(double) );


    PHG4Particle* g4particle = iter->second;
    double true_energy = g4particle->get_e();
    
    // get vertex
    PHG4VtxPoint* vtx = truthinfo->GetPrimaryVtx(g4particle->get_vtx_id());

    // get track
    SvtxTrack* track = trackeval->best_track_from(g4particle);
      

    // Loop over cal layers
    for (int i=0; i<_num_cal_layers; ++i) 
    {

      // pull the clusters
      string clusternodename = "CLUSTER_" + _cal_names[i];
      RawClusterContainer *clusterList = findNode::getClass<RawClusterContainer>(topNode,clusternodename.c_str());
      if (!clusterList) {
        cerr << PHWHERE << " ERROR: Can't find node " << clusternodename << endl;
        exit(-1);
      }    


      if (track) 
      {
			  cal_e_3x3[i] = track->get_cal_energy_3x3(_cal_types[i]);
			  cal_e_5x5[i] = track->get_cal_energy_5x5(_cal_types[i]);
			  cal_e_clu[i] = track->get_cal_cluster_e(_cal_types[i]);
      }

      // for no track (no charge), use straight line projection from initial 
      // momentum vector to calorimeter
      double px = g4particle->get_px();
      double py = g4particle->get_py();
      double pz = g4particle->get_pz();
      // shift origin to vertex
      px = px + vtx->get_x();
      py = py + vtx->get_y();
      pz = pz + vtx->get_z();
      
      double phi = atan2(py,px);
      double eta = asinh(pz/sqrt(px*px+py*py));

    
      // loop over all clusters and find nearest
      // I stole this code from PHG4SvtxTrackProjection
      double min_r = DBL_MAX;
      double min_index = -9999;
      double min_dphi = NAN;
      double min_deta = NAN;
      double min_e = NAN;
      for (unsigned int k = 0; k < clusterList->size(); ++k) 
      {
        RawCluster *cluster = clusterList->getCluster(k);

        double dphi = 
              atan2(sin(phi-cluster->get_phi()),cos(phi-cluster->get_phi()));
		    double deta = eta-cluster->get_eta();
		    double r = sqrt(pow(dphi,2)+pow(deta,2));

		    if (r < min_r) 
        {
		      min_index = k;
		      min_r = r;
		      min_dphi = dphi;
		      min_deta = deta;
		      min_e = cluster->get_energy();
		    }
      }

      if (min_index != -9999) 
      {
        cal_e_clu_nocharge[i] = min_e;

    	  if (verbosity > 1) 
        {
		      cout << " nearest cluster dphi = " << min_dphi << " deta = " << min_deta << " e = " << min_e << endl;
		    }
      }

    } // end calorimeter layer loop


    // Fill histograms
    if (track) 
    {
      double cal_e_3x3_tot = cal_e_3x3[0] + cal_e_3x3[1] + cal_e_3x3[2];
      _h_calresponse_3x3->Fill(cal_e_3x3_tot/true_energy);
    
      double cal_e_5x5_tot = cal_e_5x5[0] + cal_e_5x5[1] + cal_e_5x5[2];
      _h_calresponse_5x5->Fill(cal_e_5x5_tot/true_energy);
    
      double cal_e_clu_tot = cal_e_clu[0] + cal_e_clu[1] + cal_e_clu[2];
      _h_calresponse_cluster->Fill(cal_e_clu_tot/true_energy);
    }
    
    double cal_e_clu_nocharge_tot = cal_e_clu_nocharge[0] + 
                         cal_e_clu_nocharge[1] + 
                         cal_e_clu_nocharge[2];
    _h_calresponse_cluster_nocharge->Fill(cal_e_clu_nocharge_tot/true_energy);

  }// truth loop

  return 0;
}


int CaloResponse::End(PHCompositeNode *topNode) {
 return 0;
}
