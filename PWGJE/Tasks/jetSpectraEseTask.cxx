// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file jetSpectraEseTask.cxx
/// \brief jet spectra analysis framework with ESE (19/08/2024)
///
/// \author Joachim C. K. B. Hansen, Lund University

#include <string>

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"

#include "Common/DataModel/EseTable.h"
#include "Common/DataModel/Qvectors.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"
#include "Framework/StaticFor.h"

struct JetSpectraEseTask {
  ConfigurableAxis binJetPt{"binJetPt", {200, 0., 200.}, ""};
  ConfigurableAxis bindPhi{"bindPhi", {100, -TMath::Pi() - 1, TMath::Pi() + 1}, ""};
  ConfigurableAxis binESE{"binESE", {100, 0, 100}, ""};
  ConfigurableAxis binCos{"binCos", {100, -1.05, 1.05}, ""};
  ConfigurableAxis binOccupancy{"binOccupancy", {5000, 0, 25000}, ""};

  Configurable<float> jetPtMin{"jetPtMin", 5.0, "minimum jet pT cut"};
  Configurable<float> jetR{"jetR", 0.2, "jet resolution parameter"};
  Configurable<float> vertexZCut{"vertexZCut", 10.0, "vertex z cut"};
  Configurable<std::vector<float>> CentRange{"CentRange", {30, 50}, "centrality region of interest"};
  Configurable<double> leadingJetPtCut{"fLeadingJetPtCut", 5.0, "leading jet pT cut"};

  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};

  Configurable<int> fColSwitch{"fColSwitch", 0, "collision switch"};

  Configurable<bool> cfgEvSelOccupancy{"cfgEvSelOccupancy", false, "Flag for occupancy cut"};
  Configurable<std::vector<int>> cfgCutOccupancy{"cfgCutOccupancy", {0, 500}, "Occupancy cut"};
  Configurable<std::vector<float>> cfgOccupancyPtCut{"cfgPtCut", {0, 100}, "pT cut"};

  AxisSpec jetPtAxis = {binJetPt, "#it{p}_{T,jet}"};
  AxisSpec dPhiAxis = {bindPhi, "#Delta#phi"};
  AxisSpec eseAxis = {binESE, "#it{q}_{2}"};
  AxisSpec cosAxis = {binCos, ""};
  AxisSpec occAxis = {binOccupancy, "Occupancy"};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  int eventSelection = -1;
  int trackSelection = -1;

  void init(o2::framework::InitContext&)
  {
    eventSelection = jetderiveddatautilities::initialiseEventSelection(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

    LOGF(info, "jetSpectraEse::init()");

    switch (fColSwitch) {
      case 0:
        LOGF(info, "JetSpectraEseTask::init() - using data");
        registry.add("hEventCounter", "event status;event status;entries", {HistType::kTH1F, {{10, 0.0, 10.0}}});
        registry.add("hJetPt", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{jetPtAxis}}});
        registry.add("hJetPt_bkgsub", "jet pT background sub;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{jetPtAxis}}});
        registry.add("hJetEta", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}});
        registry.add("hJetPhi", "jet #phi;#phi_{jet};entries", {HistType::kTH1F, {{80, -1.0, 7.}}});
        registry.add("hRho", ";#rho;entries", {HistType::kTH1F, {{100, 0, 200.}}});
        registry.add("hJetArea", ";area_{jet};entries", {HistType::kTH1F, {{100, 0, 10.}}});
        registry.add("hdPhi", "#Delta#phi;entries;", {HistType::kTH1F, {{dPhiAxis}}});
        registry.add("hCentJetPtdPhiq2", "", {HistType::kTHnSparseF, {{100,0,100},{jetPtAxis}, {dPhiAxis}, {eseAxis}}});
        registry.add("hPsi2FT0C", ";Centrality; #Psi_{2}", {HistType::kTH2F, {{100, 0, 100}, {150, -2.5, 2.5}}});
        registry.add("hPsi2FT0A", ";Centrality; #Psi_{2}", {HistType::kTH2F, {{100, 0, 100}, {150, -2.5, 2.5}}});
        registry.add("hPsi2FV0A", ";Centrality; #Psi_{2}", {HistType::kTH2F, {{100, 0, 100}, {150, -2.5, 2.5}}});
        registry.add("hPsi2TPCpos", ";Centrality; #Psi_{2}", {HistType::kTH2F, {{100, 0, 100}, {150, -2.5, 2.5}}});
        registry.add("hPsi2TPCneg", ";Centrality; #Psi_{2}", {HistType::kTH2F, {{100, 0, 100}, {150, -2.5, 2.5}}});
        registry.add("hCosPsi2FT0CmFT0A", ";Centrality;cos(2(#Psi_{2}^{FT0C}-#Psi_{2}^{FT0A}));#it{q}_{2}", {HistType::kTH3F, {{100, 0, 100}, {cosAxis}, {eseAxis}}});
        registry.add("hCosPsi2FT0CmFV0A", ";Centrality;cos(2(#Psi_{2}^{FT0C}-#Psi_{2}^{FV0A}));#it{q}_{2}", {HistType::kTH3F, {{100, 0, 100}, {cosAxis}, {eseAxis}}});
        registry.add("hCosPsi2FV0AmFT0A", ";Centrality;cos(2(#Psi_{2}^{FT0C}-#Psi_{2}^{FV0A}));#it{q}_{2}", {HistType::kTH3F, {{100, 0, 100}, {cosAxis}, {eseAxis}}});
        registry.add("hCosPsi2FT0AmFT0C", ";Centrality;cos(2(#Psi_{2}^{FT0A}-#Psi_{2}^{FT0C}));#it{q}_{2}", {HistType::kTH3F, {{100, 0, 100}, {cosAxis}, {eseAxis}}});
        registry.add("hCosPsi2FT0AmFV0A", ";Centrality;cos(2(#Psi_{2}^{FT0C}-#Psi_{2}^{FV0A}));#it{q}_{2}", {HistType::kTH3F, {{100, 0, 100}, {cosAxis}, {eseAxis}}});
        registry.add("hCosPsi2FV0AmFT0C", ";Centrality;cos(2(#Psi_{2}^{FV0A}-#Psi_{2}^{FT0C}));#it{q}_{2}", {HistType::kTH3F, {{100, 0, 100}, {cosAxis}, {eseAxis}}});
        registry.add("hCosPsi2TPCposmTPCneg", ";Centrality;cos(2(#Psi_{2}^{TPCpos}-#Psi_{2}^{TPCneg}));#it{q}_{2}", {HistType::kTH3F, {{100, 0, 100}, {cosAxis}, {eseAxis}}});
        registry.add("hCosPsi2TPCposmFV0A", ";Centrality;cos(2(#Psi_{2}^{TPCpos}-#Psi_{2}^{FV0A}));#it{q}_{2}", {HistType::kTH3F, {{100, 0, 100}, {cosAxis}, {eseAxis}}});
        registry.add("hCosPsi2TPCnegmFV0A", ";Centrality;cos(2(#Psi_{2}^{TPCneg}-#Psi_{2}^{FV0A}));#it{q}_{2}", {HistType::kTH3F, {{100, 0, 100}, {cosAxis}, {eseAxis}}});
        registry.add("hCosPsi2FT0AmTPCpos", ";Centrality;cos(2(#Psi_{2}^{FT0A}-#Psi_{2}^{TPCpos}));#it{q}_{2}", {HistType::kTH3F, {{100, 0, 100}, {cosAxis}, {eseAxis}}});
        registry.add("hCosPsi2FT0AmTPCneg", ";Centrality;cos(2(#Psi_{2}^{FT0A}-#Psi_{2}^{TPCneg}));#it{q}_{2}", {HistType::kTH3F, {{100, 0, 100}, {cosAxis}, {eseAxis}}});
        registry.add("hCosPsi2FT0CmTPCpos", ";Centrality;cos(2(#Psi_{2}^{FT0C}-#Psi_{2}^{TPCpos}));#it{q}_{2}", {HistType::kTH3F, {{100, 0, 100}, {cosAxis}, {eseAxis}}});
        registry.add("hCosPsi2FT0CmTPCneg", ";Centrality;cos(2(#Psi_{2}^{FT0C}-#Psi_{2}^{TPCneg}));#it{q}_{2}", {HistType::kTH3F, {{100, 0, 100}, {cosAxis}, {eseAxis}}});
        break;
      case 1:
        LOGF(info, "JetSpectraEseTask::init() - using MC");
        registry.add("h_mc_collisions", "event status;event status;entries", {HistType::kTH1F, {{10, 0.0, 10.0}}});
        registry.add("h_part_jet_pt", "particle level jet pT;#it{p}_{T,jet part} (GeV/#it{c});entries", {HistType::kTH1F, {{jetPtAxis}}});
        registry.add("h_part_jet_eta", "particle level jet #eta;#eta_{jet part};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}});
        registry.add("h_part_jet_phi", "particle level jet #phi;#phi_{jet part};entries", {HistType::kTH1F, {{80, -1.0, 7.}}});
        registry.add("h_part_jet_pt_match", "particle level jet pT;#it{p}_{T,jet part} (GeV/#it{c});entries", {HistType::kTH1F, {{jetPtAxis}}});
        registry.add("h_part_jet_eta_match", "particle level jet #eta;#eta_{jet part};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}});
        registry.add("h_part_jet_phi_match", "particle level jet #phi;#phi_{jet part};entries", {HistType::kTH1F, {{80, -1.0, 7.}}});
        registry.add("h_detector_jet_pt", "detector level jet pT;#it{p}_{T,jet det} (GeV/#it{c});entries", {HistType::kTH1F, {{jetPtAxis}}});
        registry.add("h_detector_jet_eta", "detector level jet #eta;#eta_{jet det};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}});
        registry.add("h_detector_jet_phi", "detector level jet #phi;#phi_{jet det};entries", {HistType::kTH1F, {{80, -1.0, 7.}}});
        registry.add("h_matched_jets_pt_delta", "#it{p}_{T,jet part}; det - part", {HistType::kTH2F, {{jetPtAxis}, {200, -20., 20.0}}});
        registry.add("h_matched_jets_eta_delta", "#eta_{jet part}; det - part", {HistType::kTH2F, {{100, -1.0, 1.0}, {200, -20.0, 20.0}}});
        registry.add("h_matched_jets_phi_delta", "#phi_{jet part}; det - part", {HistType::kTH2F, {{80, -1.0, 7.}, {200, -20.0, 20.}}});
        registry.add("h_response_mat_match", "#it{p}_{T, jet det}; #it{p}_{T, jet part}", HistType::kTH2F, {jetPtAxis, jetPtAxis});
        break;
      case 2:
        LOGF(info, "JetSpectraEseTask::init() - using Occupancy processing");
        registry.add("hEventCounterOcc", "event status;event status;entries", {HistType::kTH1F, {{10, 0.0, 10.0}}});
        registry.add("hTrackPt", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTHnSparseF, {{100,0,100}, {100,0,100},{eseAxis}, {occAxis}}});
        registry.add("hTrackEta", "track #eta;#eta_{track};entries", {HistType::kTH3F, {{100,0,100},{100, -1.0, 1.0},{occAxis}}});
        registry.add("hTrackPhi", "track #phi;#phi_{track};entries", {HistType::kTH3F, {{100,0,100},{80, -1.0, 7.},{occAxis}}});
        registry.add("hOccupancy", "Occupancy;Occupancy;entries", {HistType::kTH1F, {{occAxis}}});
        registry.add("hPsiOccupancy", "Occupancy;#Psi_{2};entries", {HistType::kTH3F, {{100,0,100}, {150,-2.5,2.5}, {occAxis}}});
        break;
    }
  }

  Filter jetCuts = aod::jet::pt > jetPtMin&& aod::jet::r == nround(jetR.node() * 100.0f) && nabs(aod::jet::eta) < 0.9f - jetR;
  Filter colFilter = nabs(aod::jcollision::posZ) < vertexZCut;
  Filter mcCollisionFilter = nabs(aod::jmccollision::posZ) < vertexZCut;

  void processESEDataCharged(soa::Join<aod::JetCollisions, aod::JCollisionPIs, aod::BkgChargedRhos>::iterator const& collision,
                             soa::Join<aod::Collisions, aod::CentFT0Cs, aod::QvectorFT0CVecs, aod::QvectorFT0AVecs, aod::QvectorFV0AVecs, aod::QvectorTPCposVecs, aod::QvectorTPCnegVecs, aod::QPercentileFT0Cs> const&,
                             soa::Filtered<aod::ChargedJets> const& jets,
                             aod::JetTracks const& tracks)
  {
    float counter{0.5f};
    registry.fill(HIST("hEventCounter"), counter++);
    const auto originalCollision = collision.collision_as<soa::Join<aod::Collisions, aod::CentFT0Cs, aod::QvectorFT0CVecs, aod::QvectorFT0AVecs, aod::QvectorFV0AVecs, aod::QvectorTPCposVecs, aod::QvectorTPCnegVecs, aod::QPercentileFT0Cs>>();
    registry.fill(HIST("hEventCounter"), counter++);

    const auto vPsi2 = procEP<true>(originalCollision);
    const auto qPerc = originalCollision.qPERCFT0C();
    if (qPerc[0] < 0)
      return;
    registry.fill(HIST("hEventCounter"), counter++);

    if (!jetderiveddatautilities::selectCollision(collision, eventSelection))
      return;

    registry.fill(HIST("hEventCounter"), counter++);

    if (!isAcceptedLeadTrack(tracks))
      return;

    if (cfgEvSelOccupancy) {
      auto occupancy = collision.trackOccupancyInTimeRange();
      if (occupancy < cfgCutOccupancy->at(0) || occupancy > cfgCutOccupancy->at(1))
        registry.fill(HIST("hEventCounter"), counter++);
        return;
    }

    registry.fill(HIST("hEventCounter"), counter++);
    registry.fill(HIST("hRho"), collision.rho());
    for (auto const& jet : jets) {
      float jetpT_bkgsub = jet.pt() - (collision.rho() * jet.area());
      registry.fill(HIST("hJetPt"), jet.pt());
      registry.fill(HIST("hJetPt_bkgsub"), jetpT_bkgsub);
      registry.fill(HIST("hJetEta"), jet.eta());
      registry.fill(HIST("hJetPhi"), jet.phi());
      registry.fill(HIST("hJetArea"), jet.area());

      float dPhi = RecoDecay::constrainAngle(jet.phi() - vPsi2, -o2::constants::math::PI);
      registry.fill(HIST("hdPhi"), dPhi);
      registry.fill(HIST("hCentJetPtdPhiq2"), originalCollision.centFT0C(), jetpT_bkgsub, dPhi, qPerc[0]);
    }
    registry.fill(HIST("hEventCounter"), counter++);

    if (originalCollision.centFT0C() < CentRange->at(0) || originalCollision.centFT0C() > CentRange->at(1)) /* for counting */
      return;
    registry.fill(HIST("hEventCounter"), counter++);
  }
  PROCESS_SWITCH(JetSpectraEseTask, processESEDataCharged, "process ese collisions", true);

  void processESEOccupancy( soa::Join<aod::JetCollisions, aod::JCollisionPIs>::iterator const& collision, 
                            soa::Join<aod::JetTracks, aod::JTrackPIs> const& tracks, 
                            soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection> const&,
                            soa::Join<aod::Collisions, aod::CentFT0Cs, aod::QPercentileFT0Cs, aod::QvectorFT0AVecs> const&)
  {
    float count{0.5};
    registry.fill(HIST("hEventCounterOcc"), count++);
    const auto orgCol = collision.collision_as<soa::Join<aod::Collisions, aod::CentFT0Cs, aod::QPercentileFT0Cs, aod::QvectorFT0AVecs>>();
    const auto vPsi2 = procEP<false>(orgCol);
    const auto qPerc = orgCol.qPERCFT0C();

    auto occupancy = collision.trackOccupancyInTimeRange();
    registry.fill(HIST("hPsiOccupancy"), orgCol.centFT0C(), vPsi2, occupancy);
    registry.fill(HIST("hOccupancy"), occupancy);
    
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) 
      return;
    registry.fill(HIST("hEventCounterOcc"), count++);

    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection))
        continue;

      
      registry.fill(HIST("hTrackPt"),orgCol.centFT0C(), track.pt(), qPerc[0], occupancy);
      if (track.pt() < cfgOccupancyPtCut->at(0) || track.pt() > cfgOccupancyPtCut->at(1))
        continue;
      registry.fill(HIST("hTrackEta"),orgCol.centFT0C(), track.eta(), occupancy);
      registry.fill(HIST("hTrackPhi"),orgCol.centFT0C(), track.phi(), occupancy);
    }

  }
  PROCESS_SWITCH(JetSpectraEseTask, processESEOccupancy, "process occupancy", false);

 
  void processMCParticleLevel(soa::Filtered<aod::ChargedMCParticleLevelJets>::iterator const& jet)
  {
    registry.fill(HIST("h_part_jet_pt"), jet.pt());
    registry.fill(HIST("h_part_jet_eta"), jet.eta());
    registry.fill(HIST("h_part_jet_phi"), jet.phi());
  }
  PROCESS_SWITCH(JetSpectraEseTask, processMCParticleLevel, "jets on particle level MC", false);

  using JetMCPTable = soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>>;
  void processMCChargedMatched(soa::Filtered<aod::JetCollisionsMCD>::iterator const& collision,
                               soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>> const& mcdjets,
                               JetMCPTable const&,
                               aod::JetTracks const&,
                               aod::JetParticles const&)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection))
      return;

    float counter{0.5f};
    registry.fill(HIST("h_mc_collisions"), counter++);
    for (const auto& mcdjet : mcdjets) {

      registry.fill(HIST("h_detector_jet_pt"), mcdjet.pt());
      registry.fill(HIST("h_detector_jet_eta"), mcdjet.eta());
      registry.fill(HIST("h_detector_jet_phi"), mcdjet.phi());
      for (auto& mcpjet : mcdjet.template matchedJetGeo_as<JetMCPTable>()) {

        registry.fill(HIST("h_part_jet_pt_match"), mcpjet.pt());
        registry.fill(HIST("h_part_jet_eta_match"), mcpjet.eta());
        registry.fill(HIST("h_part_jet_phi_match"), mcpjet.phi());

        registry.fill(HIST("h_matched_jets_pt_delta"), mcpjet.pt(), mcdjet.pt() - mcpjet.pt());
        registry.fill(HIST("h_matched_jets_phi_delta"), mcpjet.phi(), mcdjet.phi() - mcpjet.phi());
        registry.fill(HIST("h_matched_jets_eta_delta"), mcpjet.eta(), mcdjet.eta() - mcpjet.eta());

        registry.fill(HIST("h_response_mat_match"), mcdjet.pt(), mcpjet.pt());
      }
    }
    registry.fill(HIST("h_mc_collisions"), counter++);
  }
  PROCESS_SWITCH(JetSpectraEseTask, processMCChargedMatched, "jet matched mcp and mcd", false);

  template <typename T>
  bool isAcceptedLeadTrack(T const& tracks)
  {
    double leadingTrackpT{0.0};
    for (const auto& track : tracks) {
      if (track.pt() > leadingJetPtCut) {
        if (track.pt() > leadingTrackpT) {
          leadingTrackpT = track.pt();
        }
      }
    }
    if (leadingTrackpT == 0.0)
      return false;
    else
      return true;
  }
  static constexpr const char* cosList[] = {"hCosPsi2FT0CmFT0A","hCosPsi2FT0CmFV0A","hCosPsi2FV0AmFT0A","hCosPsi2FT0AmFT0C","hCosPsi2FT0AmFV0A","hCosPsi2FV0AmFT0C","hCosPsi2TPCposmTPCneg","hCosPsi2TPCposmFV0A","hCosPsi2TPCnegmFV0A","hCosPsi2FT0AmTPCpos","hCosPsi2FT0AmTPCneg","hCosPsi2FT0CmTPCpos","hCosPsi2FT0CmTPCneg"};
  template <bool fill, typename EPCol>
  float procEP(EPCol const& vec)
  { 
    constexpr std::array<float, 2> ampCut{1e-8, 0.0};
    auto computeEP = [&ampCut](auto sumAmpl, auto qImVec, auto qReVec, auto det) { return sumAmpl > ampCut[det] ? 0.5 * std::atan2(qImVec[0], qReVec[0]) : 0.0; };
    std::map<std::string, float> epMap;
    epMap["FT0A"] = computeEP(vec.sumAmplFT0A(), vec.qvecFT0AImVec(), vec.qvecFT0AReVec(), 0);
    if constexpr(fill){
      epMap["FV0A"] = computeEP(vec.sumAmplFV0A(), vec.qvecFV0AImVec(), vec.qvecFV0AReVec(), 0);
      epMap["FT0C"] = computeEP(vec.sumAmplFT0C(), vec.qvecFT0CImVec(), vec.qvecFT0CReVec(), 0);
      epMap["TPCpos"] = computeEP(vec.nTrkTPCpos(), vec.qvecTPCposImVec(), vec.qvecTPCposReVec(), 1);
      epMap["TPCneg"] = computeEP(vec.nTrkTPCneg(), vec.qvecTPCnegImVec(), vec.qvecTPCnegReVec(), 1);
      fillEP(std::make_index_sequence<5>{}, vec, epMap);

      auto cosPsi = [](float psiX, float psiY) { return (static_cast<double>(psiX) == 0.0 || static_cast<double>(psiY) == 0.0) ? 0.0f : std::cos(2.0 * (psiX - psiY)); };
      std::vector<float> epCorrContainer{};
      for (const auto& name : cosList) {
        const auto [x, y] = getName(name);
        epCorrContainer.push_back(cosPsi(epMap[std::string(x)], epMap[std::string(y)]));
      }
      fillEPCos(std::make_index_sequence<13>{}, vec, epCorrContainer);
    }
    return epMap["FT0A"];
  }
  template <std::size_t... Idx, typename collision>
  void fillEPCos(const std::index_sequence<Idx...>&, const collision& col, const std::vector<float>& Corr) {
    (registry.fill(HIST(cosList[Idx]), col.centFT0C(), Corr[Idx], col.qPERCFT0C()[0]), ...);
  }
  static constexpr const char* epList[] = {"hPsi2FT0A", "hPsi2FV0A", "hPsi2FT0C", "hPsi2TPCpos", "hPsi2TPCneg"};
  template <std::size_t... Idx, typename collision>
  void fillEP(const std::index_sequence<Idx...>&, const collision& col, const std::map<std::string, float>& epMap) { 
    (registry.fill(HIST(epList[Idx]), col.centFT0C(), epMap.at(std::string(removePrefix(epList[Idx])))), ...);
  }
  std::pair<std::string_view, std::string_view> getName(std::string_view str) {
    size_t pos = str.find("m");
    std::string_view x = str.substr(8, pos - 8);
    std::string_view y = str.substr(pos + 1);
    return {x, y};
  }
  constexpr std::string_view removePrefix(std::string_view str) {
    constexpr std::string_view prefix = "hPsi2";
    return str.substr(prefix.size());
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<JetSpectraEseTask>(cfgc)}; }
