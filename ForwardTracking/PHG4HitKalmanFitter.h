/*!
 *  \file		PHG4HitKalmanFitter.h
 *  \brief		Kalman Filter based on smeared truth PHG4Hit
 *  \details	Kalman Filter based on smeared truth PHG4Hit
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef __PHG4HitKalmanFitter_H__
#define __PHG4HitKalmanFitter_H__

#include <fun4all/SubsysReco.h>
#include <g4main/PHG4HitContainer.h>
#include <GenFit/MeasuredStateOnPlane.h>
#include <phgenfit/Measurement.h>
#include <string>
#include <vector>

#include "/afs/rhic.bnl.gov/x8664_sl6/opt/sphenix/core/root-5.34.34/include/TVector3.h"

class PHG4Particle;
namespace PHGenFit {
class PlanarMeasurement;
} /* namespace PHGenFit */

namespace PHGenFit {
class Track;
} /* namespace PHGenFit */

namespace genfit {
class GFRaveVertexFactory;
} /* namespace genfit */

class SvtxTrack;
namespace PHGenFit {
class Fitter;
} /* namespace PHGenFit */

class SvtxTrackMap;
class SvtxVertexMap;
class SvtxVertex;
class PHCompositeNode;
class PHG4TruthInfoContainer;
class SvtxClusterMap;
class SvtxEvalStack;
class TFile;
class TTree;

class PHG4HitKalmanFitter: public SubsysReco {
public:

	//! Default constructor
	PHG4HitKalmanFitter(const std::string &name = "PHG4HitKalmanFitter");

	//! dtor
	~PHG4HitKalmanFitter();

	//!Initialization, called for initialization
	int Init(PHCompositeNode *);

	//!Initialization Run, called for initialization of a run
	int InitRun(PHCompositeNode *);

	//!Process Event, called for each event
	int process_event(PHCompositeNode *);

	//!End, write and close files
	int End(PHCompositeNode *);

	/// set verbosity
	void Verbosity(int verb) {
		verbosity = verb; // SubsysReco verbosity
	}

	bool is_do_evt_display() const {
		return _do_evt_display;
	}

	void set_do_evt_display(bool doEvtDisplay) {
		_do_evt_display = doEvtDisplay;
	}

	double get_FGEM_phi_resolution() const {
		return _FGEM_phi_resolution;
	}

	void set_FGEM_phi_resolution(double fgemPhiResolution) {
		_FGEM_phi_resolution = fgemPhiResolution;
	}

	double get_FGEM_r_resolution() const {
		return _FGEM_r_resolution;
	}

	void set_FGEM_r_resolution(double fgemRResolution) {
		_FGEM_r_resolution = fgemRResolution;
	}

	const std::string& get_fit_alg_name() const {
		return _fit_alg_name;
	}

	void set_fit_alg_name(const std::string& fitAlgName) {
		_fit_alg_name = fitAlgName;
	}

	const std::string& get_mag_field_file_name() const {
		return _mag_field_file_name;
	}

	void set_mag_field_file_name(const std::string& magFieldFileName) {
		_mag_field_file_name = magFieldFileName;
	}

	float get_mag_field_re_scaling_factor() const {
		return _mag_field_re_scaling_factor;
	}

	void set_mag_field_re_scaling_factor(float magFieldReScalingFactor) {
		_mag_field_re_scaling_factor = magFieldReScalingFactor;
	}

	bool is_reverse_mag_field() const {
		return _reverse_mag_field;
	}

	void set_reverse_mag_field(bool reverseMagField) {
		_reverse_mag_field = reverseMagField;
	}

	double get_pat_rec_hit_finding_eff() const {
		return _pat_rec_hit_finding_eff;
	}

	void set_pat_rec_hit_finding_eff(double patRecHitFindingEff) {
		if(!(patRecHitFindingEff>=0&&patRecHitFindingEff<=1)) {
			std::cout<<"ERROR: _pat_rec_hit_finding_eff out of range! \n";
		}
		_pat_rec_hit_finding_eff = patRecHitFindingEff;
	}

	double get_pat_rec_nosise_prob() const {
		return _pat_rec_nosise_prob;
	}

	void set_pat_rec_nosise_prob(double patRecNosiseProb) {
		if(!(patRecNosiseProb <= 1. && patRecNosiseProb >= 0)) {
			std::cout<<"ERROR: _pat_rec_nosise_prob out of range! \n";
			return;
		}
		_pat_rec_nosise_prob = patRecNosiseProb;
	}

private:

	/*!
	 * Create needed nodes.
	 */
	int CreateNodes(PHCompositeNode *);

	/*!
	 * Get all the all the required nodes off the node tree.
	 */
	int GetNodes(PHCompositeNode *);

	/*!
	 *
	 */
	int PseudoPatternRecognition(const PHG4Particle* particle,
			std::vector<PHGenFit::Measurement*> & meas_out, TVector3& seed_pos,
			TVector3& seed_mom, TMatrixDSym& seed_cov, const bool do_smearing = true);

	PHGenFit::PlanarMeasurement* PHG4HitToMeasurementVerticalPlane(const PHG4Hit* g4hit);

	PHGenFit::PlanarMeasurement* VertexMeasurement(const TVector3 &vtx, const double dr,
			const double dphi);

	/*!
	 * Make SvtxTrack from PHGenFit::Track
	 */
	SvtxTrack* MakeSvtxTrack(const PHGenFit::Track* phgf_track_in);

	//! Event counter
	int _event;


	//! Input Node pointers
	PHG4TruthInfoContainer* _truth_container;

	std::vector<PHG4HitContainer*> _phg4hits;
	std::vector<std::string> _phg4hits_names;

	//! Output Node pointers
	SvtxTrackMap* _trackmap_out;


	/*!
	 *	GenFit fitter interface
	 */
	PHGenFit::Fitter* _fitter;

	//!
	std::string _mag_field_file_name;

	//! rescale mag field, modify the original mag field read in
	float _mag_field_re_scaling_factor;

	//! Switch to reverse Magnetic field
	bool _reverse_mag_field;

	//!
	std::string _fit_alg_name;

	//!
	bool _do_evt_display;

	/*!
	 * For PseudoPatternRecognition function.
	 */

	double _FGEM_phi_resolution;

	double _FGEM_r_resolution;

	//!
	double _pat_rec_hit_finding_eff;

	//!
	double _pat_rec_nosise_prob;



};

#endif /*__PHG4HitKalmanFitter_H__*/





















