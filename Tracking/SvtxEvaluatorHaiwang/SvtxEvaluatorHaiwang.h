#ifndef SVTXEVALUATOR_HAIWANG_H__
#define SVTXEVALUATOR_HAIWANG_H__

//===============================================
/// \file SvtxEvaluatorHaiwang.h
/// \brief Compares reconstructed tracks to truth particles
/// \author Michael P. McCumber (revised SVTX version)
//===============================================


#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;

class SvtxEvalStack;
class TFile;
class TNtuple;

/// \class SvtxEvaluatorHaiwang
///
/// \brief Compares reconstructed tracks to truth particles
///
/// Plan: This module will trace the reconstructed tracks back to
/// the greatest contributor Monte Carlo particle and then
/// test one against the other.
///
class SvtxEvaluatorHaiwang : public SubsysReco {
  
public:
 
  SvtxEvaluatorHaiwang(const std::string &name = "SVTXEVALUATOR",
                const std::string &filename = "g4eval.root");
  virtual ~SvtxEvaluatorHaiwang() {}
		
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  void set_strict(bool b) {_strict = b;}
  
  void do_vertex_eval(bool b) {_do_vertex_eval = b;}
  void do_gpoint_eval(bool b) {_do_gpoint_eval = b;}
  void do_g4hit_eval(bool b) {_do_g4hit_eval = b;}
  void do_hit_eval(bool b) {_do_hit_eval = b;}
  void do_cluster_eval(bool b) {_do_cluster_eval = b;}
  void do_gtrack_eval(bool b) {_do_gtrack_eval = b;}
  void do_track_eval(bool b) {_do_track_eval = b;}

  void scan_for_embedded(bool b) {_scan_for_embedded = b;}
  
 private:

  unsigned int _ievent;

  // eval stack
  SvtxEvalStack* _svtxevalstack;
  
  //----------------------------------
  // evaluator output ntuples

  bool _strict;
  unsigned int _errors;
  
  bool _do_vertex_eval;
  bool _do_gpoint_eval;
  bool _do_g4hit_eval;
  bool _do_hit_eval;
  bool _do_cluster_eval;
  bool _do_gtrack_eval;
  bool _do_track_eval;

  bool _scan_for_embedded;
  
  TNtuple *_ntp_vertex;
  TNtuple *_ntp_gpoint;
  TNtuple *_ntp_g4hit;
  TNtuple *_ntp_hit;
  TNtuple *_ntp_cluster;
  TNtuple *_ntp_gtrack;
  TNtuple *_ntp_track;

  // evaluator output file
  std::string _filename;
  TFile *_tfile;

  // output subroutines
  void fillOutputNtuples(PHCompositeNode* topNode); ///< dump the evaluator information into ntuple for external analysis
  void printInputInfo(PHCompositeNode* topNode);    ///< print out the input object information (debugging upstream components)
  void printOutputInfo(PHCompositeNode* topNode);   ///< print out the ancestry information for detailed diagnosis
};

#endif // SVTXEVALUATOR_HAIWANG_H__
