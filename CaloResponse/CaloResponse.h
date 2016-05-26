#ifndef __CaloResponse_H__
#define __CaloResponse_H__

#include <fun4all/SubsysReco.h>
#include <string>
#include <TH1D.h>
#include <g4hough/SvtxTrack.h>

#include <vector>
#include <string>

class PHCompositeNode;
class RawClusterContainer;

class CaloResponse: public SubsysReco {

public: 

  CaloResponse(const std::string &name="CaloResponse");

  int Init(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

private:

  // event counter
  unsigned long long _event;

  const int _num_cal_layers = 3;
  TH1D* _h_calresponse_3x3;
  TH1D* _h_calresponse_5x5;
  TH1D* _h_calresponse_cluster;
  TH1D* _h_calresponse_cluster_nocharge;

  std::vector<SvtxTrack::CAL_LAYER> _cal_types;
  std::vector<std::string> _cal_names;
};

#endif // __CaloResponse_H__
