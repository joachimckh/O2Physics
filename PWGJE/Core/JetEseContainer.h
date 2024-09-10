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
#ifndef PWGJE_CORE_JETESECONTAINER_H_
#define PWGJE_CORE_JETESECONTAINER_H_

#include <vector>
#include <string>
#include <utility>
#include <algorithm>
#include <memory>

#include "TNamed.h"
#include "TObjArray.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TFile.h"
#include "TCollection.h"
#include "TString.h"
#include "TMath.h"



class JetEseContainer : public TNamed {
  public:
    JetEseContainer();
    explicit JetEseContainer(const char* name);
    ~JetEseContainer();

    void Init(bool autoRange);
    void SetPtBinRange(int nBins, int nMin, int nMax);
    void SetEseRange(int Resolution, int min, int max);
    void SetdPhiBinRange(int nBins, int nMin, int nMax);
    
    void Fill(const float &pt);
    void Fill(const float &pt, const float &perc);
    void Fill(const float &pt, const float &perc, const float &phi, const float &vPsi);

    void FilldPhi(const float &phi);

    std::shared_ptr<TObjArray> GetDataArray() { return jet_data; }
    
    void MergeEsePt();
    float DeltaRPhi(Double_t mphi, Double_t vphi);

  private:
    std::shared_ptr<TObjArray> jet_data;
    std::shared_ptr<TAxis> ptA; //!
    std::shared_ptr<TAxis> eseA; //!
    std::shared_ptr<TAxis> dPhiA; //!


    
    const char* GetPtName(const char* pf = "_pt") {
      return Form("jet%s", pf);
    };
    
  ClassDef(JetEseContainer, 1); // calibration class

};

#endif // PWGJE_CORE_JETESECONTAINER_H_
