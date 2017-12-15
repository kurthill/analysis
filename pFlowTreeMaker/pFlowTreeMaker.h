#ifndef __PFLOWTREEMAKER_H__
#define __PFLOWTREEMAKER_H__

// --- need to check all these includes...
#include <fun4all/SubsysReco.h>
#include <vector>

#include "TTree.h"
#include "TFile.h"

class PHCompositeNode;

class pFlowTreeMaker: public SubsysReco
{

 public:

  pFlowTreeMaker(const std::string &name="pFlowTreeMaker.root");

  int Init(PHCompositeNode*);
  int process_event(PHCompositeNode*);
  int End(PHCompositeNode*);

 private:

  TFile *_f;

  TTree *_tree;

  std::string _foutname;
  float _magfield;

  float _b_vz;

  int _b_tower_n;
  int _b_tower_layer[10000];
  float _b_tower_E[10000];
  float _b_tower_eta[10000];
  float _b_tower_phi[10000];

  int _b_cluster_n;
  unsigned int _b_cluster_id[500];
  float _b_cluster_chi2[500];
  float _b_cluster_prob[500];
  float _b_cluster_e[500];
  float _b_cluster_pt[500];
  float _b_cluster_eta[500];
  float _b_cluster_phi[500];
  int _b_cluster_layer[500];

  int _b_particle_n;
  int _b_particle_pid[1000];
  int _b_particle_track_id[1000];
  float _b_particle_p[1000];
  float _b_particle_pt[1000];
  float _b_particle_eta[1000];
  float _b_particle_phi[1000];
  float _b_particle_proj_eta_layer0[1000];
  float _b_particle_proj_eta_layer1[1000];
  float _b_particle_proj_eta_layer2[1000];
  float _b_particle_proj_phi_layer0[1000];
  float _b_particle_proj_phi_layer1[1000];
  float _b_particle_proj_phi_layer2[1000];

  int _b_track_n;
  unsigned int _b_track_id[1000];
  float _b_track_p[1000];
  float _b_track_pt[1000];
  float _b_track_eta[1000];
  float _b_track_phi[1000];
  // track calorimeter projections
  unsigned int _b_track_cal_cluster_id_layer0[1000];
  unsigned int _b_track_cal_cluster_id_layer1[1000];
  unsigned int _b_track_cal_cluster_id_layer2[1000];
  float _b_track_deta_layer0[1000];
  float _b_track_deta_layer1[1000];
  float _b_track_deta_layer2[1000];
  float _b_track_dphi_layer0[1000];
  float _b_track_dphi_layer1[1000];
  float _b_track_dphi_layer2[1000];
  float _b_track_proj_eta_layer0[1000];
  float _b_track_proj_eta_layer1[1000];
  float _b_track_proj_eta_layer2[1000];
  float _b_track_proj_phi_layer0[1000];
  float _b_track_proj_phi_layer1[1000];
  float _b_track_proj_phi_layer2[1000];
  float _b_track_cal_energy_3x3_layer0[1000];
  float _b_track_cal_energy_3x3_layer1[1000];
  float _b_track_cal_energy_3x3_layer2[1000];
  float _b_track_cal_energy_5x5_layer0[1000];
  float _b_track_cal_energy_5x5_layer1[1000];
  float _b_track_cal_energy_5x5_layer2[1000];

  int _b_jet_n;
  float _b_jet_e[500];
  float _b_jet_pt[500];
  float _b_jet_eta[500];
  float _b_jet_phi[500];

};

#endif // __PFLOWTREEMAKER_H__
