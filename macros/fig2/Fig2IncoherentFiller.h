/*
 * Fig2IncoherentFiller.h
 *
 * Purpose
 * -------
 * Shared histogram filler for incoherent/breakup Monte Carlo used in the
 * Fig. 2 studies. It uses the common Fig2 kinematics and creates, for each
 * of the four reference photon-energy bins, the m(gamma gamma), DeltaE_pi
 * before the DeltaPhi selection, DeltaE_pi-vs-DeltaPhi, and DeltaE_pi after
 * the standard 8, 10 and 12 degree DeltaPhi cuts.
 *
 * Important version detail
 * ------------------------
 * The v5 header uses the wide DeltaE_pi range -400..50 MeV (450 bins).
 * Earlier v4 copies of the calling macro can look identical while including
 * a narrower older filler. This repository keeps this v5 header together
 * with the Acqu template maker so their behaviour is unambiguous.
 *
 * This file is original analysis code. Only this documentation header has
 * been added; class behaviour and histogram definitions are unchanged.
 */

#ifndef FIG2_INCOHERENT_FILLER_H
#define FIG2_INCOHERENT_FILLER_H

#include "Fig2Common.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TDirectory.h"
#include "TString.h"

class Fig2IncoherentFiller {
public:
    Fig2IncoherentFiller()
    {
        for (int ie = 0; ie < Fig2::kNEnergyBins; ++ie) {
            hMgg[ie] = new TH1D(Form("inc_mgg_E%d", ie),
                                Form("Incoherent MC, E_{#gamma}=%.0f MeV;m_{#gamma#gamma} (MeV);Weighted events",
                                     Fig2::kEnergyCenter[ie]),
                                180, 0.0, 360.0);
            hDeltaEVsDeltaPhi[ie] =
                new TH2D(Form("inc_deltaE_vs_deltaPhi_E%d", ie),
                         Form("Incoherent MC, E_{#gamma}=%.0f MeV;#Delta E_{#pi^{0}}^{*} (MeV);#Delta#Phi (deg)",
                              Fig2::kEnergyCenter[ie]),
                         450, -400.0, 50.0, 100, -10.0, 40.0);
            hDeltaENoPhi[ie] =
                new TH1D(Form("inc_deltaE_E%d_nodphi", ie),
                         Form("Incoherent MC, %.0f MeV, no #Delta#Phi cut;#Delta E_{#pi^{0}}^{*} (MeV);Weighted events",
                              Fig2::kEnergyCenter[ie]),
                         450, -400.0, 50.0);

            hMgg[ie]->Sumw2();
            hDeltaEVsDeltaPhi[ie]->Sumw2();
            hDeltaENoPhi[ie]->Sumw2();

            for (int ic = 0; ic < Fig2::kNCuts; ++ic) {
                const int cutLabel = (int)Fig2::kDeltaPhiCutDeg[ic];
                hDeltaE[ie][ic] =
                    new TH1D(Form("inc_deltaE_E%d_dphi%d", ie, cutLabel),
                             Form("Incoherent MC, %.0f MeV, #Delta#Phi < %d^{#circ};#Delta E_{#pi^{0}}^{*} (MeV);Weighted events",
                                  Fig2::kEnergyCenter[ie], cutLabel),
                             450, -400.0, 50.0);
                hDeltaE[ie][ic]->Sumw2();
            }
        }
    }

    /*
     * Fill this directly from an existing event loop once the best gamma pair
     * has been selected and the invariant-mass cut has been applied.
     */
    void Fill(double eGammaMeV,
              const TLorentzVector& gamma1,
              const TLorentzVector& gamma2,
              double weight = 1.0,
              bool fillMassBeforeCut = true,
              double mggMinMeV = 110.0,
              double mggMaxMeV = 155.0)
    {
        const int ie = Fig2::FindEnergyBin(eGammaMeV);
        if (ie < 0) return;

        const TLorentzVector pi0 = gamma1 + gamma2;
        const double mgg = pi0.M();

        if (fillMassBeforeCut) hMgg[ie]->Fill(mgg, weight);
        if (mgg < mggMinMeV || mgg > mggMaxMeV) return;

        const double deltaE =
            Fig2::DeltaECoherentHypothesis(eGammaMeV, pi0);
        const double deltaPhi =
            Fig2::DeltaPhiCoherentDeg(eGammaMeV, gamma1, gamma2);

        if (!TMath::Finite(deltaE) || !TMath::Finite(deltaPhi)) return;

        hDeltaENoPhi[ie]->Fill(deltaE, weight);
        hDeltaEVsDeltaPhi[ie]->Fill(deltaE, deltaPhi, weight);

        for (int ic = 0; ic < Fig2::kNCuts; ++ic) {
            if (deltaPhi < Fig2::kDeltaPhiCutDeg[ic])
                hDeltaE[ie][ic]->Fill(deltaE, weight);
        }
    }

    void Write(TDirectory* directory = 0)
    {
        TDirectory* old = gDirectory;
        if (directory) directory->cd();

        for (int ie = 0; ie < Fig2::kNEnergyBins; ++ie) {
            hMgg[ie]->Write();
            hDeltaENoPhi[ie]->Write();
            hDeltaEVsDeltaPhi[ie]->Write();

            for (int ic = 0; ic < Fig2::kNCuts; ++ic)
                hDeltaE[ie][ic]->Write();
        }

        if (old) old->cd();
    }

    TH1D* hMgg[Fig2::kNEnergyBins];
    TH1D* hDeltaENoPhi[Fig2::kNEnergyBins];
    TH2D* hDeltaEVsDeltaPhi[Fig2::kNEnergyBins];
    TH1D* hDeltaE[Fig2::kNEnergyBins][Fig2::kNCuts];
};

#endif
