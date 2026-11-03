// -*- C++ -*-
#ifndef EFT2OBS_RIVETPLUGINS_WGQQ_H
#define EFT2OBS_RIVETPLUGINS_WGQQ_H

#include "HepMC3/GenParticle.h"
#include "Rivet/Analysis.hh"
#include <map>
#include <string>
#include <vector>

namespace Rivet {

int find_bins(double value, const std::vector<double> binning);

struct HadronicWDecay {
  HepMC3::ConstGenParticlePtr boson;
  HepMC3::ConstGenParticlePtr quark1;
  HepMC3::ConstGenParticlePtr quark2;
};

// Find terminal hadronic W decays, deduplicating same-signed W copy chains.
std::vector<HadronicWDecay> FindHadronicWDecays(
    const std::vector<HepMC3::ConstGenParticlePtr>& particles);

class WGQQ : public Analysis {
public:
  struct WGammaRivetVariables {
    bool w_flag;

    double fj_pt;
    double fj_eta;
    double fj_phi;
    double fj_M;
    double fj_MSD;

    double j0_pt;
    double j0_eta;
    double j0_phi;
    double j0_M;

    double j1_pt;
    double j1_eta;
    double j1_phi;
    double j1_M;

    double p0_pt;
    double p0_eta;
    double p0_phi;
    double p0_M;

    double phi_quark;
    double phi_sjet;
    double phi_sjet_matched; //may be different from phi_sjet due to the jets order
    double wg_M;

    double w_pt;

    WGammaRivetVariables();
    void resetVars();
  };

  struct WGSystem {

    FourMomentum wg_system;

    FourMomentum c_w_boson;
    FourMomentum c_subjet1;
    FourMomentum c_subjet2;
    FourMomentum c_photon;

    FourMomentum r_w_boson;
    FourMomentum r_subjet1;
    FourMomentum r_subjet2;
    FourMomentum r_photon;

    WGSystem(FourMomentum const &j1, FourMomentum const &j2, FourMomentum const &pho,
             bool verbose);

    double Phi();
    double SymPhi();
  };

  static constexpr double kFatJetPtCut = 200.;
  static constexpr double kFatJetAbsEtaCut = 4.7;
  static constexpr double kFatJetDR = 0.8;
  static constexpr double kJetDR = 0.4;
  static constexpr double kJetPtCut = 30.;

  static constexpr double kPhotonPtCut = 200.;
  static constexpr double kPhotonAbsEtaCut = 2.5;
  static constexpr double kPhotonPtMax = 2000.;

  std::vector<double> pta_bins_;
  std::vector<double> pta_full_bins_ = {200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1500, 2000};
  std::vector<double> mwa_bins_ = {400, 500, 700, 1000, 2000, 4000};
  std::vector<double> phi_bins_ = {
      -PI,     -PI * 4 / 5., -PI * 3 / 5., -PI * 2 / 5., -PI / 5., 0,
      PI / 5., PI * 2 / 5.,  PI * 3 / 5.,  PI * 4 / 5.,  PI};
  WGammaRivetVariables vars_;
  std::map<std::string, Histo1DPtr> _h;

  RIVET_DEFAULT_ANALYSIS_CTOR(WGQQ);

  void init() override;

  /// Perform the per-event analysis
  void analyze(const Event &event) override;

  void finalize() override;
};

} // namespace Rivet

#endif
