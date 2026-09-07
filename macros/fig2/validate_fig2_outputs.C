/*
 * validate_fig2_outputs.C
 *
 * Purpose
 * -------
 * Small consistency checker for the standardized Fig. 2 files. It verifies
 * that, for all four photon-energy bins and DeltaPhi cuts 8/10/12 degrees,
 * the data, coherent and incoherent histograms expected by the overlay
 * workflow are present.
 *
 * Default inputs
 * --------------
 * fig2_data_coherent.root
 * fig2_incoherent_templates.root
 *
 * Expected histogram families
 * ---------------------------
 * data_E<bin>_dphi<cut>
 * coh_E<bin>_dphi<cut>
 * inc_deltaE_E<bin>_dphi<cut>
 *
 * Usage
 * -----
 * root -l -b -q validate_fig2_outputs.C
 *
 * No files are written; integrals/status are printed to stdout.
 * Original analysis code follows unchanged; only this header was added.
 */

#include "Fig2Common.h"

#include "TFile.h"
#include "TH1.h"
#include "TString.h"

#include <iostream>

void validate_fig2_outputs(
    const char* dataCoherentFile = "fig2_data_coherent.root",
    const char* incoherentFile = "fig2_incoherent_templates.root")
{
    TFile* dataCoherent =
        TFile::Open(dataCoherentFile, "READ");

    TFile* incoherent =
        TFile::Open(incoherentFile, "READ");

    if (!dataCoherent || dataCoherent->IsZombie()) {
        std::cerr << "Cannot open "
                  << dataCoherentFile
                  << std::endl;
        return;
    }

    if (!incoherent || incoherent->IsZombie()) {
        std::cerr << "Cannot open "
                  << incoherentFile
                  << std::endl;
        dataCoherent->Close();
        return;
    }

    bool ok = true;

    for (int energyIndex = 0;
         energyIndex < Fig2::kNEnergyBins;
         ++energyIndex) {
        for (int cutIndex = 0;
             cutIndex < Fig2::kNCuts;
             ++cutIndex) {
            const int cut =
                static_cast<int>(
                    Fig2::kDeltaPhiCutDeg[cutIndex]);

            const char* dataName =
                Form("data_E%d_dphi%d",
                     energyIndex,
                     cut);

            const char* coherentName =
                Form("coh_E%d_dphi%d",
                     energyIndex,
                     cut);

            const char* incoherentName =
                Form("inc_deltaE_E%d_dphi%d",
                     energyIndex,
                     cut);

            TH1* data =
                dynamic_cast<TH1*>(
                    dataCoherent->Get(dataName));

            TH1* coherent =
                dynamic_cast<TH1*>(
                    dataCoherent->Get(coherentName));

            TH1* inc =
                dynamic_cast<TH1*>(
                    incoherent->Get(incoherentName));

            if (!data || !coherent || !inc) {
                ok = false;
                std::cerr
                    << "Missing:"
                    << (!data ? Form(" %s", dataName) : "")
                    << (!coherent ? Form(" %s", coherentName) : "")
                    << (!inc ? Form(" %s", incoherentName) : "")
                    << std::endl;
                continue;
            }

            std::cout
                << "E" << energyIndex
                << ", dphi" << cut
                << ": data=" << data->Integral()
                << ", coherent=" << coherent->Integral()
                << ", incoherent=" << inc->Integral()
                << std::endl;
        }
    }

    dataCoherent->Close();
    incoherent->Close();

    if (ok)
        std::cout
            << "All 36 required histograms are present."
            << std::endl;
}
