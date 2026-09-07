/*
 * LowECoherentFiller.h
 *
 * Histogram container/filler for the dedicated low-energy coherent-MC branch.
 * It books m(gamma gamma), DeltaE without a DeltaPhi cut, DeltaE-vs-DeltaPhi,
 * and DeltaE after the 8/10/12-degree cuts for the four low-energy bins defined
 * in LowECoherentCommon.h.
 *
 * This header is preserved functionally unchanged from the supplied analysis.
 */
#ifndef LOWE_COHERENT_FILLER_H
#define LOWE_COHERENT_FILLER_H

#include "LowECoherentCommon.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TDirectory.h"
#include "TString.h"

class LowECoherentFiller {
public:
    LowECoherentFiller()
    {
        for (int ie = 0; ie < LowECoh::kNEnergyBins; ++ie) {
            hMgg[ie] =
                new TH1D(Form("coh_lowE_mgg_E%d", ie),
                         Form("Coherent MC, %.0f-%.0f MeV;"
                              "m_{#gamma#gamma} (MeV);Events",
                              LowECoh::kEnergyLow[ie],
                              LowECoh::kEnergyHigh[ie]),
                         180, 0.0, 360.0);

            hDeltaENoPhi[ie] =
                new TH1D(Form("coh_lowE_deltaE_E%d_nodphi", ie),
                         Form("Coherent MC, %.0f-%.0f MeV, no #Delta#Phi cut;"
                              "#Delta E_{#pi^{0}}^{*} (MeV);Events",
                              LowECoh::kEnergyLow[ie],
                              LowECoh::kEnergyHigh[ie]),
                         130, -80.0, 50.0);

            hDeltaEVsDeltaPhi[ie] =
                new TH2D(Form("coh_lowE_deltaE_vs_deltaPhi_E%d", ie),
                         Form("Coherent MC, %.0f-%.0f MeV;"
                              "#Delta E_{#pi^{0}}^{*} (MeV);#Delta#Phi (deg)",
                              LowECoh::kEnergyLow[ie],
                              LowECoh::kEnergyHigh[ie]),
                         130, -80.0, 50.0,
                         120, -20.0, 40.0);

            hMgg[ie]->Sumw2();
            hDeltaENoPhi[ie]->Sumw2();
            hDeltaEVsDeltaPhi[ie]->Sumw2();

            for (int ic = 0; ic < LowECoh::kNCuts; ++ic) {
                const int cut =
                    (int)LowECoh::kDeltaPhiCutDeg[ic];

                hDeltaE[ie][ic] =
                    new TH1D(Form("coh_lowE_deltaE_E%d_dphi%d",
                                  ie, cut),
                             Form("Coherent MC, %.0f-%.0f MeV, "
                                  "#Delta#Phi < %d^{#circ};"
                                  "#Delta E_{#pi^{0}}^{*} (MeV);Events",
                                  LowECoh::kEnergyLow[ie],
                                  LowECoh::kEnergyHigh[ie],
                                  cut),
                             130, -80.0, 50.0);

                hDeltaE[ie][ic]->Sumw2();
            }
        }
    }

    void Fill(double eGammaMeV,
              const TLorentzVector& gamma1,
              const TLorentzVector& gamma2,
              double weight = 1.0,
              double mggMinMeV = 110.0,
              double mggMaxMeV = 155.0)
    {
        const int ie = LowECoh::FindEnergyBin(eGammaMeV);
        if (ie < 0) return;

        const TLorentzVector pi0 = gamma1 + gamma2;
        const double mgg = pi0.M();

        hMgg[ie]->Fill(mgg, weight);
        if (mgg < mggMinMeV || mgg > mggMaxMeV)
            return;

        const double deltaE =
            LowECoh::DeltaECoherentHypothesis(eGammaMeV, pi0);

        const double deltaPhi =
            LowECoh::DeltaPhiCoherentDeg(eGammaMeV,
                                         gamma1, gamma2);

        if (!TMath::Finite(deltaE) ||
            !TMath::Finite(deltaPhi))
            return;

        hDeltaENoPhi[ie]->Fill(deltaE, weight);
        hDeltaEVsDeltaPhi[ie]->Fill(deltaE, deltaPhi, weight);

        for (int ic = 0; ic < LowECoh::kNCuts; ++ic) {
            if (deltaPhi < LowECoh::kDeltaPhiCutDeg[ic])
                hDeltaE[ie][ic]->Fill(deltaE, weight);
        }
    }

    void Write(TDirectory* directory = 0)
    {
        TDirectory* old = gDirectory;
        if (directory) directory->cd();

        for (int ie = 0; ie < LowECoh::kNEnergyBins; ++ie) {
            hMgg[ie]->Write();
            hDeltaENoPhi[ie]->Write();
            hDeltaEVsDeltaPhi[ie]->Write();

            for (int ic = 0; ic < LowECoh::kNCuts; ++ic)
                hDeltaE[ie][ic]->Write();
        }

        if (old) old->cd();
    }

    TH1D* hMgg[LowECoh::kNEnergyBins];
    TH1D* hDeltaENoPhi[LowECoh::kNEnergyBins];
    TH2D* hDeltaEVsDeltaPhi[LowECoh::kNEnergyBins];
    TH1D* hDeltaE[LowECoh::kNEnergyBins][LowECoh::kNCuts];
};

#endif
