#ifndef __SoftLeptonTaggingTruth_H__
#define __SoftLeptonTaggingTruth_H__

#include <fun4all/SubsysReco.h>
#include <phool/PHCompositeNode.h>

#include <string>
#include <set>
#include <map>
#include <utility>      // std::pair, std::make_pair
#include <stdint.h>
#ifndef __CINT__
#include <memory>
#endif

#include <TString.h>

class PHCompositeNode;
class Fun4AllHistoManager;
class TH1F;
class JetEvalStack;
class JetTruthEval;
class Jet;
class TAxis;
class PHG4TruthInfoContainer;

/// \class SoftLeptonTaggingTruth
class SoftLeptonTaggingTruth : public SubsysReco
{

public:

  enum enu_flags
  {
    //! spectrum of truth jets
    kProcessTruthSpectrum = 1 << 1,

    kDefaultFlag = kProcessTruthSpectrum
  };

  SoftLeptonTaggingTruth(const std::string & truth_jet, enu_flags flags =
      kDefaultFlag);
  virtual
  ~SoftLeptonTaggingTruth();

  //! add reco jet to the process list
  //! @return number of reco jet on list
  int
  add_reco_jet(const std::string & reco_jet)
  {
    _reco_jets.insert(reco_jet);
    return _reco_jets.size();
  }

  uint32_t
  get_flags() const
  {
    return _flags;
  }

  void
  set_flags(enu_flags flags)
  {
    _flags = (uint32_t) flags;
  }

  void
  set_flag(enu_flags flag)
  {
    _flags |= (uint32_t) flag;
  }

  bool
  flag(enu_flags flag)
  {
    return _flags & flag;
  }

  void
  reset_flag(enu_flags flag)
  {
    _flags &= ~(uint32_t) flag;
  }

  //! Energy ratio difference cut from 1 for matched jets
  double
  get_jet_match_E_Ratio() const
  {
    return _jet_match_E_Ratio;
  }

  //! Energy ratio difference cut from 1 for matched jets
  void
  set_jet_match_E_Ratio(double jetMatchDERatio)
  {
    _jet_match_E_Ratio = jetMatchDERatio;
  }

  //! Eta difference cut for matched jets
  double
  get_jet_match_dR() const
  {
    return _jet_match_dR;
  }

  //! Eta difference cut for matched jets
  void
  set_jet_match_dR(double jetMatchDEta)
  {
    _jet_match_dR = jetMatchDEta;
  }

  //! Phi difference cut for matched jets
  double
  get_jet_match_dca() const
  {
    return _jet_match_dca;
  }

  //! Phi difference cut for matched jets
  void
  set_jet_match_dca(double jetMatchDPhi)
  {
    _jet_match_dca = jetMatchDPhi;
  }

  //! set eta range
  void
  set_eta_range(double low, double high);

  int
  Init(PHCompositeNode *topNode);
  int
  InitRun(PHCompositeNode *topNode);
  int
  process_event(PHCompositeNode *topNode);
  int
  End(PHCompositeNode *topNode);

  //! Get a pointer to the default hist manager for QA modules
  static Fun4AllHistoManager *
  getHistoManager();

private:

  void
  useLogBins(TAxis * axis);

  int
  Init_Spectrum(PHCompositeNode *topNode, const std::string & jet_name);
  int
  process_Spectrum(PHCompositeNode *topNode, const std::string & jet_name,
      const bool is_reco_jet);

  //! common prefix for QA histograms
  std::string
  get_histo_prefix(const std::string & src_jet_name = "",
      const std::string & reco_jet_name = "");

#ifndef __CINT__
  //! cache the jet evaluation modules
  typedef std::map<std::string, std::shared_ptr<JetEvalStack>> jetevalstacks_map;
  jetevalstacks_map _jetevalstacks;
  std::shared_ptr<JetTruthEval> _jettrutheval;
#endif

  //! truth jet name
  std::string _truth_jet;

  //! list of reco jet
  std::set<std::string> _reco_jets;

  PHG4TruthInfoContainer* _truth_container;

  uint32_t _flags;

  //! eta range
  std::pair<double, double> eta_range;

  //! string description of eta range
  //! @return TString as ROOT likes
  TString
  get_eta_range_str(const char * eta_name = "#eta_{Jet}") const;

  //! acceptance cut on jet object
  bool
  jet_acceptance_cut(const Jet * jet) const;

  //! Eta difference cut for matched jets
  double _jet_match_dR;

  //! Phi difference cut for matched jets
  double _jet_match_dca;

  //! Energy ratio difference cut from 1 for matched jets
  double _jet_match_E_Ratio;
};

#endif // __CALOEVALUATOR_H__
