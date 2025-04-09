// -*- C++ -*-
#include "HepMC/IteratorRange.h"
#include "Rivet/Analysis.hh"
#include "Rivet/Event.hh"
#include "Rivet/Math/LorentzTrans.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "fastjet/contrib/SoftDrop.hh"
namespace Rivet {

class WGQQ : public Analysis {
public:
  struct WGammaRivetVariables {
    bool is_wg_gen;

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

    double phi;
    double phi_f;

    double wg_M;

    WGammaRivetVariables() { resetVars(); }
    void resetVars() {
      is_wg_gen = false;
      fj_pt = 0.;
      fj_eta = 0.;
      fj_phi = 0.;
      fj_M = 0.;
      fj_MSD = 0.;

      p0_pt = 0.;
      p0_eta = 0.;
      p0_phi = 0.;
      p0_M = 0.;

      j0_pt = 0.;
      j0_eta = 0.;
      j0_phi = 0.;
      j0_M = 0.;
      j1_pt = 0.;
      j1_eta = 0.;
      j1_phi = 0.;
      j1_M = 0.;

      phi = 0.;
      phi_f = 0.;
      
      wg_M = 0;
    }
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

    WGSystem(PseudoJet const &j1, PseudoJet const &j2, Particle const &pho,
             bool verbose);

    double Phi();
    double SymPhi();
  };

  double fatjet_pt_cut_ = 200.;
  double fatjet_abs_eta_cut_ = 4.7;
  double fatjet_dr_ = 0.8;

  double photon_pt_cut_ = 200.;
  double photon_abs_eta_cut_ = 2.5;

  WGammaRivetVariables vars_;
  map<string, Histo1DPtr> _h;

  RIVET_DEFAULT_ANALYSIS_CTOR(WGQQ);

  void init() {
    vars_.resetVars();

    FinalState fs;
    declare(fs, "FinalState");

    // Jets - all final state particles excluding neutrinos
    VetoedFinalState vfs;
    vfs.vetoNeutrinos();
    FastJets fastjets(vfs, FastJets::ANTIKT, fatjet_dr_);
    declare(fastjets, "FJets");

    // Photons
    IdentifiedFinalState photons(fs, 0);
    photons.acceptIdPair(PID::PHOTON);
    PromptFinalState prompt_photons(photons);
    // prompt_photons.acceptMuonDecays(true);
    // prompt_photons.acceptTauDecays(true);
    declare(prompt_photons, "Photons");

    // Booking of histograms
    book(_h["fjet_sdmass"], "fjet_sdmass", 18, 20, 200);
    book(_h["sjet_phi"], "sjet_phi", 10, -PI, PI);
    book(_h["sjet_phi_f"], "sjet_phi_f", 10, PI / 2., PI / 2.);
    book(_h["wg_mass"], "wg_mass", 50, 0, 2000);
    book(_h["p_pt"], "p_pt", 35, 200, 1500);
    book(_h["fjet_pt"], "fjet_pt", 35, 200, 1500);

  }

  /// Perform the per-event analysis
  void analyze(const Event &event) {
    vars_.resetVars();
    const Particles photons =
        applyProjection<FinalState>(event, "Photons").particlesByPt(Cuts:: pT > photon_pt_cut_*GeV && Cuts::abseta < photon_abs_eta_cut_);

    if (photons.size() == 0) {
      vetoEvent;
    }
    auto p0 = photons.at(0);

    // Filter jets on pT, eta and DR with lepton and photon
    const Jets fjets =
        applyProjection<FastJets>(event, "FJets").jetsByPt([&](Jet const &j) {
          return j.pt() > fatjet_pt_cut_ &&
                 std::abs(j.eta()) < fatjet_abs_eta_cut_ &&
                 deltaR(j, p0) > fatjet_dr_;
        });
    
    fastjet::contrib::SoftDrop sd(0.0, 0.1, fatjet_dr_);
    PseudoJet SD_fjet;
    PseudoJets SD_subjets;
    bool flag_subjet = false;
    for (auto const &fj : fjets) {
      SD_fjet = sd(fj);
      ClusterSequence subjet_sdcs(SD_fjet.constituents(), JetDefinition(fastjet::antikt_algorithm, 0.4));
      PseudoJets sdsubjets = fastjet::sorted_by_pt(subjet_sdcs.inclusive_jets(10.0));
      if (sdsubjets.size() >= 2) {
        SD_subjets = sdsubjets;
        flag_subjet = true;
        break;
      }
    }
    if (photons.size() > 0) {
      vars_.p0_pt = p0.pt();
      vars_.p0_eta = p0.eta();
      vars_.p0_phi = p0.phi(PhiMapping::MINUSPI_PLUSPI);
      vars_.p0_M = p0.mass();
      _h["p_pt"]->fill(vars_.p0_pt);
    }

    if (photons.size() >= 1 && flag_subjet) {
      // Populate variables
      vars_.is_wg_gen = true;
      vars_.j0_pt = SD_subjets[0].pt();
      vars_.j0_eta = SD_subjets[0].eta();
      vars_.j0_phi = SD_subjets[0].phi_std();
      vars_.j0_M = SD_subjets[0].m();
      vars_.j1_pt = SD_subjets[0].pt();
      vars_.j1_eta = SD_subjets[0].eta();
      vars_.j1_phi = SD_subjets[0].phi_std();
      vars_.j1_M = SD_subjets[0].m();

      vars_.fj_pt = fjets[0].pt();
      vars_.fj_eta = fjets[0].eta();
      vars_.fj_phi = fjets[0].phi(PhiMapping::MINUSPI_PLUSPI);
      vars_.fj_M = fjets[0].mass();

      vars_.wg_M = (p0.momentum() + fjets[0].momentum()).mass(); 


      // Now calculate EFT phi observables
      auto wg_system = WGSystem(SD_subjets[0], SD_subjets[1], p0, false);

      vars_.phi = wg_system.Phi();
      vars_.phi_f = wg_system.SymPhi();

      _h["fjet_sdmass"]->fill(SD_fjet.m());
      _h["sjet_phi"]->fill(vars_.phi);
      _h["sjet_phi_f"]->fill(vars_.phi_f);
      _h["wg_mass"]->fill(vars_.wg_M);
      _h["fjet_pt"]->fill(vars_.fj_pt);

    }
  }

  void finalize() {
    // Scale according to cross section
    for (std::string const &x :
         {"fjet_sdmass", "sjet_phi", "sjet_phi_f", "wg_mass", "p_pt", "fjet_pt"}) {
      if (crossSection() < 0.) {
        // Assume av. evt weight gives xsec
        scale(_h[x], 1.0 / femtobarn / numEvents());
      } else {
        scale(_h[x], crossSection() / femtobarn / sumOfWeights());
      }
    }
  }
};

WGQQ::WGSystem::WGSystem(PseudoJet const &j0, PseudoJet const &j1,
                         Particle const &pho, bool verbose) {

  wg_system += momentum(j0);
  wg_system += momentum(j1);
  wg_system += pho.momentum();

  if (verbose) {
    std::cout << "> wg_system: " << wg_system << "\n";
    std::cout << "> subjet1  : " << momentum(j0) << "\n";
    std::cout << "> subjet2  : " << momentum(j1) << "\n";
    std::cout << "> photon   : " << pho.momentum() << "\n";
  }

  auto boost = LorentzTransform::mkFrameTransformFromBeta(wg_system.betaVec());

  c_subjet1 = boost.transform(momentum(j0));
  c_subjet2 = boost.transform(momentum(j1));
  c_photon = boost.transform(pho.momentum());

  if (verbose) {
    std::cout << "> c_subjet1  : " << c_subjet1 << "\n";
    std::cout << "> c_subjet2  : " << c_subjet2 << "\n";
    std::cout << "> c_photon   : " << c_photon << "\n";
  }

  FourMomentum c_w_boson;
  c_w_boson += c_subjet1;
  c_w_boson += c_subjet2;

  auto r_uvec = wg_system.vector3().unit();

  if (verbose) {
    std::cout << "> c_w_boson : " << c_w_boson << "\n";
    std::cout << "> r_uvec   : " << r_uvec << "\n";
  }

  auto z_uvec = c_w_boson.vector3().unit();
  auto y_uvec = z_uvec.cross(r_uvec).unit();
  auto x_uvec = y_uvec.cross(z_uvec).unit();

  if (verbose) {
    std::cout << "> x_uvec   : " << x_uvec << "\n";
    std::cout << "> y_uvec   : " << y_uvec << "\n";
    std::cout << "> z_uvec   : " << z_uvec << "\n";
  }

  Matrix3 rot_matrix;
  rot_matrix.setRow(0, x_uvec).setRow(1, y_uvec).setRow(2, z_uvec);
  auto rotator = LorentzTransform();
  rotator = rotator.postMult(rot_matrix);
  if (verbose) {
    std::cout << "> rotator   : " << rotator << "\n";
  }

  r_w_boson = rotator.transform(c_w_boson);
  r_subjet1 = rotator.transform(c_subjet1);
  r_subjet2 = rotator.transform(c_subjet2);
  r_photon = rotator.transform(c_photon);

  if (verbose) {
    std::cout << "> r_subjet1  : " << r_subjet1 << "\n";
    std::cout << "> r_subjet2  : " << r_subjet2 << "\n";
    std::cout << "> r_photon   : " << r_photon << "\n";
    std::cout << "> r_w_boson  : " << r_w_boson << "\n";
  }
}

double WGQQ::WGSystem::Phi() {
  return r_subjet1.phi(PhiMapping::MINUSPI_PLUSPI);
}

double WGQQ::WGSystem::SymPhi() {
  double j0_phi = Phi();
  if (j0_phi > PI / 2.) {
    return PI - j0_phi;
  } else if (j0_phi < -1. * (PI / 2.)) {
    return -1. * (PI + j0_phi);
  } else {
    return j0_phi;
  }
}

RIVET_DECLARE_PLUGIN(WGQQ);

} // namespace Rivet
