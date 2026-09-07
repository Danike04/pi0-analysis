/*
 * Fig2SaveInputs.h
 *
 * Helper used by quantify_fig2_cut_stability_root6.C to write the 4 x 3
 * grid of unscaled data and coherent-MC missing-energy histograms to a
 * standardized ROOT file. The output histogram names are used by later
 * Fig. 2 template-comparison macros.
 *
 * Depends on Fig2Common.h. Behaviour below is preserved from the original.
 */
#ifndef FIG2_SAVE_INPUTS_H
#define FIG2_SAVE_INPUTS_H

#include "Fig2Common.h"
#include "TFile.h"
#include "TH1.h"
#include "TString.h"
#include <iostream>

/*
 * Add this header to the macro that already creates the user's preferred
 * data/coherent plots. After all 12 histograms have been filled, call:
 *
 *   SaveFig2Inputs(hData, hCoherent, "fig2_data_coherent.root");
 *
 * The array order must be [energy index][DeltaPhi-cut index]:
 * energies = 224, 294, 320, 366 MeV
 * cuts     = 8, 10, 12 degrees
 */
inline bool SaveFig2Inputs(TH1* hData[Fig2::kNEnergyBins][Fig2::kNCuts],
                           TH1* hCoherent[Fig2::kNEnergyBins][Fig2::kNCuts],
                           const char* outputFile = "fig2_data_coherent.root")
{
    TFile out(outputFile, "RECREATE");
    if (out.IsZombie()) {
        std::cerr << "Cannot create " << outputFile << std::endl;
        return false;
    }

    for (int ie = 0; ie < Fig2::kNEnergyBins; ++ie) {
        for (int ic = 0; ic < Fig2::kNCuts; ++ic) {
            const int cut = (int)Fig2::kDeltaPhiCutDeg[ic];
            if (!hData[ie][ic] || !hCoherent[ie][ic]) {
                std::cerr << "Missing histogram at energy index " << ie
                          << ", cut index " << ic << std::endl;
                out.Close();
                return false;
            }

            TH1* d = (TH1*)hData[ie][ic]->Clone(
                Form("data_E%d_dphi%d", ie, cut));
            TH1* c = (TH1*)hCoherent[ie][ic]->Clone(
                Form("coh_E%d_dphi%d", ie, cut));

            d->SetDirectory(&out);
            c->SetDirectory(&out);
            d->Write();
            c->Write();
        }
    }

    out.Write();
    out.Close();
    std::cout << "Saved standardized Fig. 2 inputs to " << outputFile << std::endl;
    return true;
}

#endif
