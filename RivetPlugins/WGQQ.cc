// -*- C++ -*-
#include "WGQQ.h"
#include "Rivet/Event.hh"
#include "Rivet/Math/LorentzTrans.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Tools/ParticleName.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <iterator>
#include <set>

namespace Rivet {

int find_bins(double value, const std::vector<double> binning) {
  auto it = std::upper_bound(binning.begin(), binning.end(), value);
  if (it == binning.begin() || it == binning.end())
    return -1;
  else
    return std::distance(binning.begin(), it) - 1;
}

// use direct quarks before FSR, reject continuing W copies, and deduplicate
// copies by their earliest same-signed-PDG ancestor. No status-code cut.
std::vector<HadronicWDecay> FindHadronicWDecays(
    const std::vector<HepMC3::ConstGenParticlePtr>& particles) {
  std::vector<HadronicWDecay> result;
  std::set<HepMC3::ConstGenParticlePtr> used_roots;
  for (const auto& boson : particles) {
    if (!boson || std::abs(boson->pdg_id()) != 24) continue;
    auto root = boson;
    std::set<HepMC3::ConstGenParticlePtr> visited;
    bool cycle = false;
    while (root) {
      if (!visited.insert(root).second) {
        cycle = true;
        break;
      }
      HepMC3::ConstGenParticlePtr mother;
      for (const auto& parent : root->parents()) {
        if (parent->pdg_id() == root->pdg_id()) {
          mother = parent;
          break;
        }
      }
      if (!mother) break;
      root = mother;
    }
    if (cycle || used_roots.count(root)) continue;
    bool has_boson_copy = false;
    std::vector<HepMC3::ConstGenParticlePtr> quarks;
    for (const auto& child : boson->children()) {
      has_boson_copy |= child->pdg_id() == boson->pdg_id();
      const int id = std::abs(child->pdg_id());
      if (id >= PID::DQUARK && id <= PID::BQUARK) quarks.push_back(child);
    }
    if (has_boson_copy || quarks.size() != 2) continue;
    used_roots.insert(root);
    result.push_back({boson, quarks[0], quarks[1]});
  }
  return result;
}

WGQQ::WGammaRivetVariables::WGammaRivetVariables() { resetVars(); }

void WGQQ::WGammaRivetVariables::resetVars() {
  w_flag = false;
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

void WGQQ::init() {
  pta_bins_.clear();
  for (double x = 200; x <= 1200; x += 50) {
    pta_bins_.push_back(x);
  }
  pta_bins_.push_back(2000);
  vars_.resetVars();

  FinalState fs;
  declare(fs, "FinalState");

  // Jets - all final state particles excluding neutrinos
  VetoedFinalState vfs;
  vfs.vetoNeutrinos();
  FastJets fastjets(vfs, JetAlg::ANTIKT, kFatJetDR);
  declare(fastjets, "FJets");

  FastJets fastjets_ak4(vfs, JetAlg::ANTIKT, kJetDR);
  declare(fastjets_ak4, "Jets");

  // Photons
  IdentifiedFinalState photons(fs, 0);
  photons.acceptIdPair(PID::PHOTON);
  PromptFinalState prompt_photons(photons);
  // prompt_photons.acceptMuonDecays(true);
  // prompt_photons.acceptTauDecays(true);
  declare(prompt_photons, "Photons");

  // Booking of histograms
  book(_h["p_pt"], "p_pt", pta_full_bins_);
  book(_h["w_flag"], "w_flag", 2, 0, 2);
  book(_h["w_pt"], "w_pt", 40, 0, 1500);
  book(_h["phi_quark"], "phi_quark", phi_bins_);

  book(_h["fjet_sdmass"], "fjet_sdmass", 18, 20, 200);
  book(_h["wp_mass"], "wp_mass", mwa_bins_);
  book(_h["fjet_pt"], "fjet_pt", 35, 200, 1500);
  book(_h["phi_sjet"], "phi_sjet", phi_bins_);
  book(_h["phi_sjet_ajet"], "phi_sjet_ajet", phi_bins_);
  book(_h["phi_sjet_noajet"], "phi_sjet_noajet", phi_bins_);

  book(_h["phi_quark_matched"], "phi_quark_matched", phi_bins_);
  book(_h["phi_sjet_matched"], "phi_sjet_matched", phi_bins_);

  book(_h["ptg"], "ptg", pta_bins_);
  book(_h["ptg_ajet"], "ptg_ajet", pta_bins_);

  for (size_t i = 0; i < phi_bins_.size() - 1; i++) {
    book(_h["ptg_phi_" + to_string(i)], "ptg_phi_" + to_string(i), pta_bins_);
    book(_h["mwp_phi_" + to_string(i)], "mwp_phi_" + to_string(i), mwa_bins_);
    book(_h["ptg_phi_" + to_string(i) + "_ajet"], "ptg_phi_" + to_string(i) + "_ajet", pta_bins_);
    book(_h["mwp_phi_" + to_string(i) + "_ajet"], "mwp_phi_" + to_string(i) + "_ajet", mwa_bins_);
  }
}

void WGQQ::analyze(const Event &event) {
  vars_.resetVars();

  const Particles photons =
      apply<FinalState>(event, "Photons").particlesByPt(Cuts:: pT > kPhotonPtCut*GeV && Cuts::abseta < kPhotonAbsEtaCut);

  if (photons.size() > 0) {
    auto p0 = photons.at(0);
    // Same signal/OOA pT boundary as Acorn prepare.py, after candidate choice.
    if (p0.pt() >= kPhotonPtMax * GeV) return;
    vars_.p0_pt = p0.pt();
    vars_.p0_eta = p0.eta();
    vars_.p0_phi = p0.phi(PhiMapping::MINUSPI_PLUSPI);
    vars_.p0_M = p0.mass();
    _h["p_pt"]->fill(vars_.p0_pt / GeV);

    const auto hadronic_ws = FindHadronicWDecays(HepMCUtils::particles(event.genEvent()));
    // Match Acorn WGQQ's signal-kinematics definition: exactly one hadronic W.
    vars_.w_flag = hadronic_ws.size() == 1;
    std::vector<ConstGenParticlePtr> w_quarks;
    _h["w_flag"]->fill(vars_.w_flag);
    if (vars_.w_flag) {
      const auto& decay = hadronic_ws[0];
      w_quarks = {decay.quark1, decay.quark2};
      vars_.w_pt = decay.boson->momentum().perp();
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
        apply<FastJets>(event, "FJets").jetsByPt([&](Jet const &j) {
          return j.pt() > kFatJetPtCut &&
                std::abs(j.eta()) < kFatJetAbsEtaCut &&
                deltaR(j, p0) > kFatJetDR;
        });

    fastjet::contrib::SoftDrop sd(0.0, 0.1, kFatJetDR);
    PseudoJet SD_fjet;
    PseudoJets SD_subjets, SD_subjets_matched;
    Jet fjet;
    bool flag_subjet = false, flag_subjet_matched = false;
    for (auto const &fj : fjets) {
      SD_fjet = sd(fj);
      PseudoJets sdsubjets;
      if (SD_fjet.has_pieces())
        sdsubjets = SD_fjet.pieces();
      else if (SD_fjet.has_constituents())
        sdsubjets = SD_fjet.constituents();
      sdsubjets = fastjet::sorted_by_pt(sdsubjets);
      if (sdsubjets.size() == 2) {
        if (!flag_subjet) {
          SD_subjets = sdsubjets;
          fjet = fj;
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
      vars_.fj_pt = fjet.pt();
      vars_.fj_eta = fjet.eta();
      vars_.fj_phi = fjet.phi(PhiMapping::MINUSPI_PLUSPI);
      vars_.fj_M = fjet.mass();
      vars_.wg_M = (p0.momentum() + fjet.momentum()).mass();
      // Now calculate EFT phi observables
      auto wg_system = WGSystem(momentum(SD_subjets[0]), momentum(SD_subjets[1]), p0.momentum(), false);
      vars_.phi_sjet = wg_system.Phi();
      _h["fjet_sdmass"]->fill(SD_fjet.m() / GeV);
      _h["phi_sjet"]->fill(vars_.phi_sjet);
      _h["wp_mass"]->fill(vars_.wg_M / GeV);
      _h["fjet_pt"]->fill(vars_.fj_pt / GeV);
      const Jets jets = apply<FastJets>(event, "Jets").jetsByPt([&](Jet const &j) {
          return j.pt() > kJetPtCut && deltaR(j, p0) > kJetDR && deltaR(j, fjet) > kFatJetDR;
      });

      int phi_bin = find_bins(vars_.phi_sjet, phi_bins_);
      if (phi_bin >= 0 && jets.empty()) {
        _h["ptg_phi_" + to_string(phi_bin)]->fill(vars_.p0_pt / GeV);
        _h["mwp_phi_" + to_string(phi_bin)]->fill(vars_.wg_M / GeV);
      }
      if (phi_bin >= 0 && !jets.empty()) {
        _h["ptg_phi_" + to_string(phi_bin) + "_ajet"]->fill(vars_.p0_pt / GeV);
        _h["mwp_phi_" + to_string(phi_bin) + "_ajet"]->fill(vars_.wg_M / GeV);
      }
      if (jets.empty()) {
        _h["ptg"]->fill(vars_.p0_pt / GeV);
        _h["phi_sjet_noajet"]->fill(vars_.phi_sjet);
      }
      else {
        _h["ptg_ajet"]->fill(vars_.p0_pt / GeV);;
        _h["phi_sjet_ajet"]->fill(vars_.phi_sjet);
      }
    }

    if (flag_subjet_matched) {
      _h["phi_quark_matched"]->fill(vars_.phi_quark);
      auto wg_system = WGSystem(momentum(SD_subjets_matched[0]), momentum(SD_subjets_matched[1]), p0.momentum(), false);
      vars_.phi_sjet_matched = wg_system.Phi();
      _h["phi_sjet_matched"]->fill(vars_.phi_sjet_matched);
    }
  }
}

void WGQQ::finalize() {
  // Scale according to cross section
  for (auto &x : _h) {
    if (crossSection() < 0.) {
      // Assume av. evt weight gives xsec
      scale(x.second, 1.0 / femtobarn / numEvents());
    } else {
      scale(x.second, crossSection() / femtobarn / sumOfWeights());
    }
  }
}

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
