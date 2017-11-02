#include "pFlowTreeMaker.h"

#include <phool/getClass.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/PHCompositeNode.h>

#include "TLorentzVector.h"
#include <iostream>

#include <calotrigger/CaloTriggerInfo.h>

#include <g4cemc/RawClusterContainer.h>
#include <g4cemc/RawCluster.h>

#include <g4cemc/RawTowerGeom.h>
#include <g4cemc/RawTower.h>
#include <g4cemc/RawTowerContainer.h>
#include <g4cemc/RawTowerGeomContainer_Cylinderv1.h>
#include <g4cemc/RawTowerGeomContainer.h>

#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4VtxPoint.h>

//#include "/opt/sphenix/core/Eigen/include/eigen3/Eigen/Core"
#include <g4hough/PHG4HoughTransform.h>
//#include <g4hough/PHG4SvtxTrackProjection.h>
//#include <g4hough/PHG4GenFitTrackProjection.h>
#include <g4hough/SvtxVertexMap.h>
#include <g4hough/SvtxVertex.h>
#include <g4hough/SvtxTrackMap.h>
#include <g4hough/SvtxTrack_v1.h>
#include <g4hough/SvtxTrackState_v1.h>

#include <g4eval/SvtxEvalStack.h>
#include <g4eval/SvtxTrackEval.h>
#include <g4eval/SvtxVertexEval.h>
#include <g4eval/SvtxTruthEval.h>



pFlowTreeMaker::pFlowTreeMaker(const std::string &name) : SubsysReco("pFlowTreeMaker"),
  _foutname(name),
  _magfield(1.5)
{
}

int pFlowTreeMaker::Init(PHCompositeNode *topNode)
{

  _f = new TFile( _foutname.c_str(), "RECREATE");

  _tree = new TTree("ttree","a succulent orange tree");

  _tree->Branch("vz",&_b_vz, "vz/F");

  _tree->Branch("tower_n",&_b_tower_n, "tower_n/I");
  _tree->Branch("tower_layer",&_b_tower_layer, "tower_layer[tower_n]/I");
  _tree->Branch("tower_E",&_b_tower_E, "tower_E[tower_n]/F");
  _tree->Branch("tower_eta",&_b_tower_eta, "tower_eta[tower_n]/F");
  _tree->Branch("tower_phi",&_b_tower_phi, "tower_phi[tower_n]/F");

  _tree->Branch("cluster_n", &_b_cluster_n,"cluster_n/I");
  _tree->Branch("cluster_id", &_b_cluster_id,"cluster_id[cluster_n]/i");
  _tree->Branch("cluster_chi2", _b_cluster_chi2,"cluster_chi2[cluster_n]/F");
  _tree->Branch("cluster_prob", _b_cluster_prob,"cluster_prob[cluster_n]/F");
  _tree->Branch("cluster_e", _b_cluster_e,"cluster_e[cluster_n]/F");
  _tree->Branch("cluster_pt", _b_cluster_pt,"cluster_pt[cluster_n]/F");
  _tree->Branch("cluster_eta",_b_cluster_eta,"cluster_eta[cluster_n]/F");
  _tree->Branch("cluster_phi",_b_cluster_phi,"cluster_phi[cluster_n]/F");
  _tree->Branch("cluster_layer",_b_cluster_layer,"cluster_layer[cluster_n]/I");


  _tree->Branch("particle_n", &_b_particle_n,"particle_n/I");
  _tree->Branch("particle_pid", &_b_particle_pid,"particle_pid[particle_n]/I");
  _tree->Branch("particle_track_id", &_b_particle_track_id,"particle_track_id[particle_n]/I");
  _tree->Branch("particle_p", _b_particle_p,"particle_p[particle_n]/F");
  _tree->Branch("particle_pt", _b_particle_pt,"particle_pt[particle_n]/F");
  _tree->Branch("particle_eta", _b_particle_eta,"particle_eta[particle_n]/F");
  _tree->Branch("particle_phi", _b_particle_phi,"particle_phi[particle_n]/F");
  _tree->Branch("particle_proj_eta_layer0", _b_particle_proj_eta_layer0,"particle_proj_eta_layer0[particle_n]/F");
  _tree->Branch("particle_proj_eta_layer1", _b_particle_proj_eta_layer1,"particle_proj_eta_layer1[particle_n]/F");
  _tree->Branch("particle_proj_eta_layer2", _b_particle_proj_eta_layer2,"particle_proj_eta_layer2[particle_n]/F");
  _tree->Branch("particle_proj_phi_layer0", _b_particle_proj_phi_layer0,"particle_proj_phi_layer0[particle_n]/F");
  _tree->Branch("particle_proj_phi_layer1", _b_particle_proj_phi_layer1,"particle_proj_phi_layer1[particle_n]/F");
  _tree->Branch("particle_proj_phi_layer2", _b_particle_proj_phi_layer2,"particle_proj_phi_layer2[particle_n]/F");

  _tree->Branch("track_n", &_b_track_n,"track_n/I");
  _tree->Branch("track_id", &_b_track_id,"track_id[track_n]/i");
  _tree->Branch("track_p", _b_track_p,"track_p[track_n]/F");
  _tree->Branch("track_pt", _b_track_pt,"track_pt[track_n]/F");
  _tree->Branch("track_eta", _b_track_eta,"track_eta[track_n]/F");
  _tree->Branch("track_phi", _b_track_phi,"track_phi[track_n]/F");
  _tree->Branch("track_cal_cluster_id_layer0", _b_track_cal_cluster_id_layer0,"track_cal_cluster_id_layer0[track_n]/i");
  _tree->Branch("track_cal_cluster_id_layer1", _b_track_cal_cluster_id_layer1,"track_cal_cluster_id_layer1[track_n]/i");
  _tree->Branch("track_cal_cluster_id_layer2", _b_track_cal_cluster_id_layer2,"track_cal_cluster_id_layer2[track_n]/i");
  _tree->Branch("track_deta_layer0", _b_track_deta_layer0,"track_deta_layer0[track_n]/F");
  _tree->Branch("track_deta_layer1", _b_track_deta_layer1,"track_deta_layer1[track_n]/F");
  _tree->Branch("track_deta_layer2", _b_track_deta_layer2,"track_deta_layer2[track_n]/F");
  _tree->Branch("track_dphi_layer0", _b_track_dphi_layer0,"track_dphi_layer0[track_n]/F");
  _tree->Branch("track_dphi_layer1", _b_track_dphi_layer1,"track_dphi_layer1[track_n]/F");
  _tree->Branch("track_dphi_layer2", _b_track_dphi_layer2,"track_dphi_layer2[track_n]/F");
  _tree->Branch("track_proj_eta_layer0", _b_track_proj_eta_layer0,"track_proj_eta_layer0[track_n]/F");
  _tree->Branch("track_proj_eta_layer1", _b_track_proj_eta_layer1,"track_proj_eta_layer1[track_n]/F");
  _tree->Branch("track_proj_eta_layer2", _b_track_proj_eta_layer2,"track_proj_eta_layer2[track_n]/F");
  _tree->Branch("track_proj_phi_layer0", _b_track_proj_phi_layer0,"track_proj_phi_layer0[track_n]/F");
  _tree->Branch("track_proj_phi_layer1", _b_track_proj_phi_layer1,"track_proj_phi_layer1[track_n]/F");
  _tree->Branch("track_proj_phi_layer2", _b_track_proj_phi_layer2,"track_proj_phi_layer2[track_n]/F");
  _tree->Branch("track_cal_energy_3x3_layer0", _b_track_cal_energy_3x3_layer0,"track_cal_energy_3x3_layer0[track_n]/F");
  _tree->Branch("track_cal_energy_3x3_layer1", _b_track_cal_energy_3x3_layer1,"track_cal_energy_3x3_layer1[track_n]/F");
  _tree->Branch("track_cal_energy_3x3_layer2", _b_track_cal_energy_3x3_layer2,"track_cal_energy_3x3_layer2[track_n]/F");
  _tree->Branch("track_cal_energy_5x5_layer0", _b_track_cal_energy_5x5_layer0,"track_cal_energy_5x5_layer0[track_n]/F");
  _tree->Branch("track_cal_energy_5x5_layer1", _b_track_cal_energy_5x5_layer1,"track_cal_energy_5x5_layer1[track_n]/F");
  _tree->Branch("track_cal_energy_5x5_layer2", _b_track_cal_energy_5x5_layer2,"track_cal_energy_5x5_layer2[track_n]/F");

  return 0;

}

int pFlowTreeMaker::process_event(PHCompositeNode *topNode)
{

  std::cout << "kurthill : at process_event, tree size is: " << _tree->GetEntries() << std::endl;

  RawTowerContainer *towersEM3 = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_CEMC");
  RawTowerContainer *towersIH3 = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_HCALIN");
  RawTowerContainer *towersOH3 = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_HCALOUT");

  RawTowerGeomContainer *geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  RawTowerGeomContainer *geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  RawTowerGeomContainer *geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");

  // Get the calorimeter radii for projections
  float cal_radii[3];
  cal_radii[0] = geomEM->get_radius();
  cal_radii[1] = geomIH->get_radius();
  cal_radii[2] = geomOH->get_radius();


  // --- SvtxTrackMap
  SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode,"SvtxTrackMap");
  if ( !trackmap )
    {
      cerr << PHWHERE << " ERROR: Can't find SvtxTrackMap" << endl;
      exit(-1);
    }


  // --- Create SVTX eval stack
  SvtxEvalStack svtxevalstack(topNode);
  // --- Get evaluator objects from the eval stack
  //SvtxVertexEval*   vertexeval = svtxevalstack.get_vertex_eval();
  SvtxTrackEval*     trackeval = svtxevalstack.get_track_eval();
  SvtxTruthEval*     trutheval = svtxevalstack.get_truth_eval();



  _b_tower_n = 0;

  cout << " Looping over towers" << endl;
  
  {
    RawTowerContainer::ConstRange begin_end = towersEM3->getTowers();
    for (RawTowerContainer::ConstIterator rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter) 
    {
      RawTower *tower = rtiter->second;
      RawTowerGeom *tower_geom = geomEM->get_tower_geometry(tower->get_key());
      _b_tower_layer[ _b_tower_n ] = 0;
      _b_tower_E[ _b_tower_n ] = tower->get_energy();
      _b_tower_eta[ _b_tower_n ] = tower_geom->get_eta();
      _b_tower_phi[ _b_tower_n ] = tower_geom->get_phi();
      _b_tower_n++;
    }
  }
  {
    RawTowerContainer::ConstRange begin_end = towersIH3->getTowers();
    for (RawTowerContainer::ConstIterator rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter) 
    {
      RawTower *tower = rtiter->second;
      RawTowerGeom *tower_geom = geomIH->get_tower_geometry(tower->get_key());
      _b_tower_layer[ _b_tower_n ] = 1;
      _b_tower_E[ _b_tower_n ] = tower->get_energy();
      _b_tower_eta[ _b_tower_n ] = tower_geom->get_eta();
      _b_tower_phi[ _b_tower_n ] = tower_geom->get_phi();
      _b_tower_n++;
    }
  }
  {
    RawTowerContainer::ConstRange begin_end = towersOH3->getTowers();
    for (RawTowerContainer::ConstIterator rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter) 
    {
      RawTower *tower = rtiter->second;
      RawTowerGeom *tower_geom = geomOH->get_tower_geometry(tower->get_key());
      _b_tower_layer[ _b_tower_n ] = 2;
      _b_tower_E[ _b_tower_n ] = tower->get_energy();
      _b_tower_eta[ _b_tower_n ] = tower_geom->get_eta();
      _b_tower_phi[ _b_tower_n ] = tower_geom->get_phi();
      _b_tower_n++;
    }
  }
 

  cout << " Looping over clusters" << endl;

  const char* l[] = {"CLUSTER_CEMC","CLUSTER_HCALIN","CLUSTER_HCALOUT"};

  _b_cluster_n = 0;
  for(int layer = 0; layer < 3; layer++)
  {
    RawClusterContainer *clusters = findNode::getClass<RawClusterContainer>(topNode,Form("%s",l[layer]));

    RawClusterContainer::ConstRange begin_end = clusters->getClusters();
    RawClusterContainer::ConstIterator rtiter;

    for (rtiter = begin_end.first; rtiter !=  begin_end.second; ++rtiter) 
    {
      RawCluster *cluster = rtiter->second;

      float pt = cluster->get_energy() / cosh(  cluster->get_eta() );

      if (pt < 0.1) continue;

      _b_cluster_id[ _b_cluster_n ] =  cluster->get_id();
      _b_cluster_chi2[ _b_cluster_n ] =  cluster->get_chi2();
      _b_cluster_prob[ _b_cluster_n ] =  cluster->get_prob();
      _b_cluster_e[ _b_cluster_n ] = cluster->get_energy();
      _b_cluster_pt[ _b_cluster_n ] = pt;
      _b_cluster_eta[ _b_cluster_n ] =  cluster->get_eta();
      _b_cluster_phi[ _b_cluster_n ] =  cluster->get_phi();
      _b_cluster_layer[ _b_cluster_n ] = layer;

      _b_cluster_n++;
    }
  }


  cout << " Looping over truth particles" << endl;

  {
    PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
    PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange();

    _b_particle_n = 0;
    for ( PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter ) 
    {
      PHG4Particle* g4particle = iter->second;

      TLorentzVector t; t.SetPxPyPzE( g4particle->get_px(), g4particle->get_py(), g4particle->get_pz(), g4particle->get_e() );
      
      if (t.Pt() < 0.05) continue;
      if (fabs(t.Eta()) > 1.1) continue;
      
      int pid = g4particle->get_pid();
      int vid = g4particle->get_vtx_id();
      PHG4VtxPoint* vtx = truthinfo->GetVtx(vid);
      
      _b_vz = vtx->get_z();
      
      _b_particle_p[ _b_particle_n ] = t.P();
      _b_particle_pt[ _b_particle_n ] = t.Pt();
      _b_particle_eta[ _b_particle_n ] = t.Eta();
      _b_particle_phi[ _b_particle_n ] = t.Phi();
      _b_particle_pid[ _b_particle_n ] = pid;

      
      // Try to project the truth particle to the calorimeter layers
      // currently not working
      //
      int c;
      if(pid == 22 || pid == 130 || 2112)
        c = 0;
      else if(pid < 0)
        c = -1;
      else if(pid > 0)
        c = 1;

      //SvtxTrack_v1 truthTrack;
      SvtxTrackState_v1 truthTrack;
      //truthTrack.set_charge(c);
      //truthTrack.set_dca(0);
      truthTrack.set_x(vtx->get_x());
      truthTrack.set_y(vtx->get_y());
      truthTrack.set_z(vtx->get_z());
      truthTrack.set_px(g4particle->get_px());
      truthTrack.set_py(g4particle->get_py());
      truthTrack.set_pz(g4particle->get_pz());
      for(int ii = 0; ii < 6; ii++)
        for(int jj = 0; jj < 6; jj++)
          truthTrack.set_error(ii,jj,0);

      std::vector<double> point_EM;
      std::vector<double> point_HI;
      std::vector<double> point_HO;
      
      PHG4HoughTransform hough;
      hough.projectToRadius(&truthTrack,c,_magfield,cal_radii[0],point_EM);
      hough.projectToRadius(&truthTrack,c,_magfield,cal_radii[1],point_HI);
      hough.projectToRadius(&truthTrack,c,_magfield,cal_radii[2],point_HO);


      double x = point_EM[0];
      double y = point_EM[1];
      double z = point_EM[2];
      if(c == 0)
      {
        _b_particle_proj_phi_layer0[ _b_particle_n ] = t.Phi();
        _b_particle_proj_eta_layer0[ _b_particle_n ] = t.Eta();
      }
      else
      {
        _b_particle_proj_phi_layer0[ _b_particle_n ] = atan2(y,x);
        _b_particle_proj_eta_layer0[ _b_particle_n ] = asinh(z/sqrt(x*x+y*y));
      }

      x = point_HI[0];
      y = point_HI[1];
      z = point_HI[2];
      if(c == 0)
      {
        _b_particle_proj_phi_layer1[ _b_particle_n ] = t.Phi();
        _b_particle_proj_eta_layer1[ _b_particle_n ] = t.Eta();
      }
      else
      {
        _b_particle_proj_phi_layer1[ _b_particle_n ] = atan2(y,x);
        _b_particle_proj_eta_layer1[ _b_particle_n ] = asinh(z/sqrt(x*x+y*y));
      }

      x = point_HO[0];
      y = point_HO[1];
      z = point_HO[2];
      if(c == 0)
      {
        _b_particle_proj_phi_layer2[ _b_particle_n ] = t.Phi();
        _b_particle_proj_eta_layer2[ _b_particle_n ] = t.Eta();
      }
      else
      {
        _b_particle_proj_phi_layer2[ _b_particle_n ] = atan2(y,x);
        _b_particle_proj_eta_layer2[ _b_particle_n ] = asinh(z/sqrt(x*x+y*y));
      }

      
      // --- Get the reconsructed SvtxTrack based on the best candidate from the truth info
      SvtxTrack* track = trackeval->best_track_from(g4particle);
      if(track != NULL)
        _b_particle_track_id[ _b_particle_n ] = track->get_id();
      else
        _b_particle_track_id[ _b_particle_n ] = -1;


      //cout << "get_track_id = " << g4particle->get_track_id() << ", phi = " << t.Phi() << endl;
      //cout << "get_id = " << track->get_id() << ", phi = " << track->get_phi() << ", get_cal_dphi = " << track->get_cal_dphi(SvtxTrack::CEMC) << endl << endl;

      _b_particle_n++;
    }
  }

  //cout << " Loop over reco tracks" << endl;

  // loop over all reco particles
  _b_track_n = 0;
  for ( SvtxTrackMap::Iter iter = trackmap->begin(); iter != trackmap->end(); ++iter )
  {
    // --- Get the StxTrack object (from the iterator)
    SvtxTrack* track = iter->second;
    cout << "get_id = " << track->get_id() << ", phi = " << track->get_phi() << ", get_cal_dphi = " << track->get_cal_dphi(SvtxTrack::CEMC) << endl;

    _b_track_id[ _b_track_n ] = track->get_id();
    _b_track_p[ _b_track_n ] = track->get_p();
    _b_track_pt[ _b_track_n ] = track->get_pt();
    _b_track_eta[ _b_track_n ] = track->get_eta();
    _b_track_phi[ _b_track_n ] = track->get_phi();
    _b_track_cal_cluster_id_layer0[ _b_track_n ] = track->get_cal_cluster_id(SvtxTrack::CEMC);
    _b_track_cal_cluster_id_layer1[ _b_track_n ] = track->get_cal_cluster_id(SvtxTrack::HCALIN);
    _b_track_cal_cluster_id_layer2[ _b_track_n ] = track->get_cal_cluster_id(SvtxTrack::HCALOUT);
    _b_track_deta_layer0[ _b_track_n ] = track->get_cal_deta(SvtxTrack::CEMC);
    _b_track_deta_layer1[ _b_track_n ] = track->get_cal_deta(SvtxTrack::HCALIN);
    _b_track_deta_layer2[ _b_track_n ] = track->get_cal_deta(SvtxTrack::HCALOUT);
    _b_track_dphi_layer0[ _b_track_n ] = track->get_cal_dphi(SvtxTrack::CEMC);
    _b_track_dphi_layer1[ _b_track_n ] = track->get_cal_dphi(SvtxTrack::HCALIN);
    _b_track_dphi_layer2[ _b_track_n ] = track->get_cal_dphi(SvtxTrack::HCALOUT);

    _b_track_cal_energy_3x3_layer0[ _b_track_n ] = track->get_cal_energy_3x3(SvtxTrack::CEMC);
    _b_track_cal_energy_3x3_layer1[ _b_track_n ] = track->get_cal_energy_3x3(SvtxTrack::HCALIN);
    _b_track_cal_energy_3x3_layer2[ _b_track_n ] = track->get_cal_energy_3x3(SvtxTrack::HCALOUT);
    _b_track_cal_energy_5x5_layer0[ _b_track_n ] = track->get_cal_energy_5x5(SvtxTrack::CEMC);
    _b_track_cal_energy_5x5_layer1[ _b_track_n ] = track->get_cal_energy_5x5(SvtxTrack::HCALIN);
    _b_track_cal_energy_5x5_layer2[ _b_track_n ] = track->get_cal_energy_5x5(SvtxTrack::HCALOUT);

    // find the track projections to the calorimeters
    //
    std::vector<double> point_EM;
    std::vector<double> point_HI;
    std::vector<double> point_HO;
    
    PHG4HoughTransform hough;
    hough.projectToRadius(track,_magfield,cal_radii[0],point_EM);
    hough.projectToRadius(track,_magfield,cal_radii[1],point_HI);
    hough.projectToRadius(track,_magfield,cal_radii[2],point_HO);

    double x = point_EM[0];
    double y = point_EM[1];
    double z = point_EM[2];
    _b_track_proj_phi_layer0[ _b_track_n ] = atan2(y,x);
    _b_track_proj_eta_layer0[ _b_track_n ] = asinh(z/sqrt(x*x+y*y));;

    x = point_HI[0];
    y = point_HI[1];
    z = point_HI[2];
    _b_track_proj_phi_layer1[ _b_track_n ] = atan2(y,x);
    _b_track_proj_eta_layer1[ _b_track_n ] = asinh(z/sqrt(x*x+y*y));;

    x = point_HO[0];
    y = point_HO[1];
    z = point_HO[2];
    _b_track_proj_phi_layer2[ _b_track_n ] = atan2(y,x);
    _b_track_proj_eta_layer2[ _b_track_n ] = asinh(z/sqrt(x*x+y*y));;

    _b_track_n++;
  }

/*
  cout " Looping over towers" << endl;
  
  {
    RawTowerContainer::ConstRange begin_end = towersEM3->getTowers();
    for (RawTowerContainer::ConstIterator rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter) 
    {

      RawTower *tower = rtiter->second;
      RawTowerGeom *tower_geom = geomEM->get_tower_geometry(tower->get_key());
      _b_tower_layer[ _b_tower_n ] = 0;
      _b_tower_E[ _b_tower_n ] = tower->get_energy();
      _b_tower_eta[ _b_tower_n ] = tower_geom->get_eta();
      _b_tower_phi[ _b_tower_n ] = tower_geom->get_phi();
      
      float deta = _b_tower_eta[ _b_tower_n ] - _b_particle_eta[ 0 ] ;
      float dphi = _b_tower_phi[ _b_tower_n ] - _b_particle_phi[ 0 ] ;
      if (dphi > +3.14159) dphi -= 2 * 3.14159;
      if (dphi < -3.14159) dphi += 2 * 3.14159;
      float dR = sqrt( pow( deta, 2 ) + pow( dphi, 2 ) );
      if (dR > 0.3) continue;

      _b_tower_n++;
    }
  }
  {
    RawTowerContainer::ConstRange begin_end = towersIH3->getTowers();
    for (RawTowerContainer::ConstIterator rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter) 
    {

      RawTower *tower = rtiter->second;
      RawTowerGeom *tower_geom = geomIH->get_tower_geometry(tower->get_key());
      _b_tower_layer[ _b_tower_n ] = 1;
      _b_tower_E[ _b_tower_n ] = tower->get_energy();
      _b_tower_eta[ _b_tower_n ] = tower_geom->get_eta();
      _b_tower_phi[ _b_tower_n ] = tower_geom->get_phi();

      float deta = _b_tower_eta[ _b_tower_n ] - _b_particle_eta[ 0 ] ;
      float dphi = _b_tower_phi[ _b_tower_n ] - _b_particle_phi[ 0 ] ;
      if (dphi > +3.14159) dphi -= 2 * 3.14159;
      if (dphi < -3.14159) dphi += 2 * 3.14159;
      float dR = sqrt( pow( deta, 2 ) + pow( dphi, 2 ) );
      if (dR > 0.3) continue;

      _b_tower_n++;
    }
  }
  {
    RawTowerContainer::ConstRange begin_end = towersOH3->getTowers();
    for (RawTowerContainer::ConstIterator rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter) 
    {

      RawTower *tower = rtiter->second;
      RawTowerGeom *tower_geom = geomOH->get_tower_geometry(tower->get_key());
      _b_tower_layer[ _b_tower_n ] = 2;
      _b_tower_E[ _b_tower_n ] = tower->get_energy();
      _b_tower_eta[ _b_tower_n ] = tower_geom->get_eta();
      _b_tower_phi[ _b_tower_n ] = tower_geom->get_phi();

      float deta = _b_tower_eta[ _b_tower_n ] - _b_particle_eta[ 0 ] ;
      float dphi = _b_tower_phi[ _b_tower_n ] - _b_particle_phi[ 0 ] ;
      if (dphi > +3.14159) dphi -= 2 * 3.14159;
      if (dphi < -3.14159) dphi += 2 * 3.14159;
      float dR = sqrt( pow( deta, 2 ) + pow( dphi, 2 ) );
      if (dR > 0.3) continue;

      _b_tower_n++;
    }
  }
    
  */


  _tree->Fill();

  return 0;
}



int pFlowTreeMaker::End(PHCompositeNode *topNode)
{

  _f->Write();
  _f->Close();

  delete _tree;
  delete _f;

  return 0;
}

