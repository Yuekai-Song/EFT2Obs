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
#include "Rivet/Projections/WFinder.hh"
#include <cstdlib>
namespace Rivet {

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

    WGammaRivetVariables() { resetVars(); }
    void resetVars() {
      w_flag = true;
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

      phi_quark = 0.;
      phi_sjet = 0.;
      phi_sjet_matched = 0.;
      wg_M = 0;

      w_pt = 0.;
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

    WGSystem(FourMomentum const &j1, FourMomentum const &j2, FourMomentum const &pho,
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
    book(_h["wg_mass"], "wg_mass", 50, 0, 2000);
    book(_h["p_pt"], "p_pt", 35, 200, 1500);
    book(_h["p_pt_matched"], "p_pt_matched", 35, 200, 1500);
    book(_h["fjet_pt"], "fjet_pt", 35, 200, 1500);
    book(_h["w_pt"], "w_pt", 40, 0, 1500);
    book(_h["w_pt_matched"], "w_pt_matched", 40, 0, 1500);
    book(_h["w_flag"], "w_flag", 2, 0, 2);
    book(_h["phi_quark"], "phi_quark", 10, -PI, PI);
    book(_h["phi_quark_matched"], "phi_quark_matched", 10, -PI, PI);
    book(_h["phi_sjet"], "phi_sjet", 10, -PI, PI);
    book(_h["phi_sjet_matched"], "phi_sjet_matched", 10, -PI, PI);
    
  }

  /// Perform the per-event analysis
  void analyze(const Event &event) {
    vars_.resetVars();

    const Particles photons =
        applyProjection<FinalState>(event, "Photons").particlesByPt(Cuts:: pT > photon_pt_cut_*GeV && Cuts::abseta < photon_abs_eta_cut_);

    if (photons.size() > 0) {
      auto p0 = photons.at(0);
      vars_.p0_pt = p0.pt();
      vars_.p0_eta = p0.eta();
      vars_.p0_phi = p0.phi(PhiMapping::MINUSPI_PLUSPI);
      vars_.p0_M = p0.mass();
      _h["p_pt"]->fill(vars_.p0_pt / GeV);

      std::vector<ConstGenParticlePtr> w, w_quarks;
      for (ConstGenParticlePtr p : HepMCUtils::particles(event.genEvent())) {
        if (abs(p->pdg_id()) == PID::WPLUSBOSON && p->status() == 62) {
          w.push_back(p);
        }
      }
      if (w.size() != 1)
        vars_.w_flag = false;
      else {
        w_quarks = HepMCUtils::particles(w[0], HepMC::children);
        if (w_quarks.size() != 2)
          vars_.w_flag = false;
        else if (abs(w_quarks[0]->pdg_id() + w_quarks[1]->pdg_id()) != 1 ||
          abs(w_quarks[0]->pdg_id()) >= PID::BQUARK ||
          abs(w_quarks[1]->pdg_id()) >= PID::BQUARK ||
          w_quarks[0]->pdg_id() == 0 || w_quarks[1]->pdg_id() == 0)
          vars_.w_flag = false;
      }
      _h["w_flag"]->fill(vars_.w_flag);
      if (vars_.w_flag) {
        vars_.w_pt = w[0]->momentum().perp();
        if (w_quarks[0]->momentum().perp() < w_quarks[1]->momentum().perp())
          swap(w_quarks[0], w_quarks[1]);
        // unsigned ind = w_quarks[0]->momentum().perp() > w_quarks[1]->momentum().perp() ? 0 : 1;
        auto wg_system = WGSystem(w_quarks[0]->momentum(), w_quarks[1]->momentum(), p0.momentum(), false);
        vars_.phi_quark = wg_system.Phi();
        _h["w_pt"]->fill(vars_.w_pt / GeV);
        _h["phi_quark"]->fill(vars_.phi_quark);
      }

      // Filter jets on pT, eta and DR with lepton and photon
      const Jets fjets =
          applyProjection<FastJets>(event, "FJets").jetsByPt([&](Jet const &j) {
            return j.pt() > fatjet_pt_cut_ &&
                  std::abs(j.eta()) < fatjet_abs_eta_cut_ &&
                  deltaR(j, p0) > fatjet_dr_;
          });
      
      fastjet::contrib::SoftDrop sd(0.0, 0.1, fatjet_dr_);
      PseudoJet SD_fjet;
      PseudoJets SD_subjets, SD_subjets_matched;
      bool flag_subjet = false, flag_subjet_matched = false;
      for (auto const &fj : fjets) {
        SD_fjet = sd(fj);
        PseudoJets sdsubjets;
        if (SD_fjet.has_pieces())
          sdsubjets = SD_fjet.pieces();
        else if (SD_fjet.has_constituents())
          sdsubjets = SD_fjet.constituents();
        sdsubjets = fastjet::sorted_by_pt(sdsubjets);
        if (sdsubjets.size() >= 2) {
          if (!flag_subjet) {
            SD_subjets = sdsubjets;
            flag_subjet = true;
          }
          // Check if the subjet is matched to the quark
          if (!flag_subjet_matched && vars_.w_flag) {
            if (deltaR(momentum(sdsubjets[0]), FourMomentum(w_quarks[0]->momentum())) < 0.4 &&
                deltaR(momentum(sdsubjets[1]), FourMomentum(w_quarks[1]->momentum())) < 0.4) {
              SD_subjets_matched = {sdsubjets[0], sdsubjets[1]};
              flag_subjet_matched = true;
            } else if (deltaR(momentum(sdsubjets[1]), FourMomentum(w_quarks[0]->momentum())) < 0.4 &&
                       deltaR(momentum(sdsubjets[0]), FourMomentum(w_quarks[1]->momentum())) < 0.4) {
              SD_subjets_matched = {sdsubjets[1], sdsubjets[0]};
              flag_subjet_matched = true;
            }
          }
          if (flag_subjet && flag_subjet_matched)
            break;
        }
      }

      if (flag_subjet) {
        // Populate variables
        vars_.j0_pt = SD_subjets[0].pt();
        vars_.j0_eta = SD_subjets[0].eta();
        vars_.j0_phi = SD_subjets[0].phi_std();
        vars_.j0_M = SD_subjets[0].m();
        vars_.j1_pt = SD_subjets[1].pt();
        vars_.j1_eta = SD_subjets[1].eta();
        vars_.j1_phi = SD_subjets[1].phi_std();
        vars_.j1_M = SD_subjets[1].m();
        vars_.fj_pt = fjets[0].pt();
        vars_.fj_eta = fjets[0].eta();
        vars_.fj_phi = fjets[0].phi(PhiMapping::MINUSPI_PLUSPI);
        vars_.fj_M = fjets[0].mass();
        vars_.wg_M = (p0.momentum() + fjets[0].momentum()).mass();
        // Now calculate EFT phi observables
        auto wg_system = WGSystem(momentum(SD_subjets[0]), momentum(SD_subjets[1]), p0.momentum(), false);
        vars_.phi_sjet = wg_system.Phi();
        _h["fjet_sdmass"]->fill(SD_fjet.m() / GeV);
        _h["phi_sjet"]->fill(vars_.phi_sjet);
        _h["wg_mass"]->fill(vars_.wg_M / GeV);
        _h["fjet_pt"]->fill(vars_.fj_pt / GeV);
      }

      if (flag_subjet_matched) {
        _h["p_pt_matched"]->fill(vars_.p0_pt/ GeV);
        _h["w_pt_matched"]->fill(vars_.w_pt / GeV);
        _h["phi_quark_matched"]->fill(vars_.phi_quark);
        auto wg_system = WGSystem(momentum(SD_subjets_matched[0]), momentum(SD_subjets_matched[1]), p0.momentum(), false);
        vars_.phi_sjet_matched = wg_system.Phi();
        _h["phi_sjet_matched"]->fill(vars_.phi_sjet_matched);
      }
    }
  }

  void finalize() {
    // Scale according to cross section
    for (std::string const &x :
         {"fjet_sdmass", "wg_mass", "p_pt", "p_pt_matched", "fjet_pt", "w_pt", "w_pt_matched", "phi_sjet", "phi_sjet_matched", "phi_quark", "phi_quark_matched"}) {
      if (crossSection() < 0.) {
        // Assume av. evt weight gives xsec
        scale(_h[x], 1.0 / femtobarn / numEvents());
      } else {
        scale(_h[x], crossSection() / femtobarn / sumOfWeights());
      }
    }
  }
};

WGQQ::WGSystem::WGSystem(FourMomentum const &j0, FourMomentum const &j1,
                         FourMomentum const &pho, bool verbose) {

  wg_system += j0;
  wg_system += j1;
  wg_system += pho;

  if (verbose) {
    std::cout << "> wg_system: " << wg_system << "\n";
    std::cout << "> subjet1  : " << j0 << "\n";
    std::cout << "> subjet2  : " << j1 << "\n";
    std::cout << "> photon   : " << pho << "\n";
  }

  auto boost = LorentzTransform::mkFrameTransformFromBeta(wg_system.betaVec());

  c_subjet1 = boost.transform(j0);
  c_subjet2 = boost.transform(j1);
  c_photon = boost.transform(pho);

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
