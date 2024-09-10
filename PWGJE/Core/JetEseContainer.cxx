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

/*
Author: Joachim Hansen
*/

#include "JetEseContainer.h"


ClassImp(JetEseContainer);



JetEseContainer::JetEseContainer() : TNamed("", ""),
                               jet_data{nullptr},
                               ptA{nullptr},
                               eseA{nullptr},
                               dPhiA{nullptr}
{
}

JetEseContainer::JetEseContainer(const char* name) : TNamed(name, name),
                                             jet_data{nullptr},
                                             ptA{nullptr},
                                             eseA{nullptr},
                                             dPhiA{nullptr}
{
}

JetEseContainer::~JetEseContainer()
{
};

void JetEseContainer::Init(bool autoRange)
{

  jet_data.reset(new TObjArray());
  jet_data->SetName("jetESE_Data");
  jet_data->SetOwner(kTRUE);

  if (autoRange) {
    this->SetPtBinRange(200,0,200);
    this->SetEseRange(100,0,100);
    this->SetdPhiBinRange(100,-TMath::Pi()-1, TMath::Pi()+1);

    jet_data->Add(new TH1F(this->GetPtName("_dPhi"), ";#Delta#phi;entries", dPhiA->GetNbins(), dPhiA->GetXmin(), dPhiA->GetXmax()));
    // jet_data->Add(new TH1F(this->GetPtName("_pt"), ";#it{p}_{T,jet};entries", ptA->GetNbins(), ptA->GetXmin(), ptA->GetXmax()));
    // jet_data->Add(new TH2F(this->GetPtName("_pt_q2"), ";#it{p}_{T,jet};#it{q}_{2}", ptA->GetNbins(), ptA->GetXmin(), ptA->GetXmax(), eseA->GetNbins(), eseA->GetXmin(), eseA->GetXmax()));
    jet_data->Add(new TH3F(this->GetPtName("_pt_dPhi_q2"), ";#it{p}_{T,jet};#Delta#phi;#it{q}_{2}", ptA->GetNbins(), ptA->GetXmin(), ptA->GetXmax(), 
                                                                                    dPhiA->GetNbins(), dPhiA->GetXmin(), dPhiA->GetXmax(), 
                                                                                    eseA->GetNbins(), eseA->GetXmin(), eseA->GetXmax()));
  }
  
};

void JetEseContainer::SetPtBinRange(int nBins, int nMin, int nMax) {
  ptA = std::shared_ptr<TAxis>(new TAxis(nBins,nMin,nMax));
};

void JetEseContainer::SetEseRange(int Resolution, int min, int max) {
  eseA = std::shared_ptr<TAxis>(new TAxis(Resolution, min, max));
};

void JetEseContainer::SetdPhiBinRange(int nBins, int nMin, int nMax) {
  dPhiA = std::shared_ptr<TAxis>(new TAxis(nBins,nMin,nMax));
};

void JetEseContainer::Fill(const float &pt) {
  auto tar = jet_data;
  if (!tar)
    return;

  TH1F* th = reinterpret_cast<TH1F*>(tar->FindObject(Form("%s",this->GetPtName("_pt"))));
  if (!th) {
    tar->Add(new TH1F(this->GetPtName("_pt"), "#it{p}_{T,jet}", ptA->GetNbins(), ptA->GetXmin(), ptA->GetXmax()));
    th = reinterpret_cast<TH1F*>(tar->At(tar->GetEntries() - 1));
  }
  th->Fill(pt);
};

void JetEseContainer::Fill(const float &pt, const float &perc) {
  auto tar = jet_data;
  if (!tar)
    return;

  TH2F* th2 = reinterpret_cast<TH2F*>(tar->FindObject(Form("%s",this->GetPtName("_pt_q2"))));
  if (!th2) {
    tar->Add(new TH2F(this->GetPtName("_pt_q2"), "#it{p}_{T,jet};q_{2};", ptA->GetNbins(), ptA->GetXmin(), ptA->GetXmax(),eseA->GetNbins(), eseA->GetXmin(), eseA->GetXmax()));
    th2 = reinterpret_cast<TH2F*>(tar->At(tar->GetEntries() - 1));
  }
  th2->Fill(pt, perc);
};

void JetEseContainer::Fill(const float &pt, const float &perc, const float &phi, const float &vPsi) {
  auto tar = jet_data;
  if (!tar)
    return;

  float dphi = DeltaRPhi(phi,vPsi);
  this->FilldPhi(dphi);
  // this->Fill(pt);
  // this->Fill(pt,perc);

  TH3F* th3 = reinterpret_cast<TH3F*>(tar->FindObject(Form("%s",this->GetPtName("_pt_dPhi_q2"))));
  if (!th3) {
    tar->Add(new TH3F(this->GetPtName("_pt_dPhi_q2"), "#itp_{T,jet};#Delta#phi;q_{2};", ptA->GetNbins(), ptA->GetXmin(), ptA->GetXmax(), 
                                                                                    dPhiA->GetNbins(), dPhiA->GetXmin(), dPhiA->GetXmax(), 
                                                                                    eseA->GetNbins(), eseA->GetXmin(), eseA->GetXmax()));
    th3 = reinterpret_cast<TH3F*>(tar->At(tar->GetEntries() - 1));
  }
  th3->Fill(pt, dphi, perc);
};

void JetEseContainer::FilldPhi(const float &phi) {
  auto tar = jet_data;
  if (!tar)
    return;

  TH1F* thp = reinterpret_cast<TH1F*>(tar->FindObject(Form("%s",this->GetPtName("_dPhi"))));
  if (!thp) {
    return;
  }
  thp->Fill(phi);
};

float JetEseContainer::DeltaRPhi(Double_t mphi, Double_t vphi) {
  if (vphi < -1 * TMath::Pi()) vphi += (2 * TMath::Pi());
  else if (vphi > TMath::Pi()) vphi -= (2 * TMath::Pi());
  
  if (mphi < -1 * TMath::Pi()) mphi += (2 * TMath::Pi());
  else if (mphi > TMath::Pi()) mphi -= (2 * TMath::Pi());
  float dphi = mphi - vphi;
  
  if (dphi < -1 * TMath::Pi()) dphi += (2 * TMath::Pi());
  else if (dphi > TMath::Pi()) dphi -= (2 * TMath::Pi());
  
  return dphi; // dphi in [-Pi, Pi]                                                                                                                                                     
}
