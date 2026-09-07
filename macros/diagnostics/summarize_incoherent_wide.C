/*
 * summarize_incoherent_wide.C
 *
 * PURPOSE
 *   Produce a compact numerical summary of wide-range incoherent templates.
 *   For each main Fig.2 energy bin and DeltaPhi = 8/10/12 deg it reports event
 *   counts and survival both over the complete histogram range and over the
 *   Fig.2 window -60 to 40 MeV.
 *
 * INPUT HISTOGRAMS
 *   inc_deltaE_Ei_nodphi and inc_deltaE_Ei_dphi{8,10,12}.
 *
 * DEFAULT OUTPUT
 *   incoherent_4bins_summary_clean.txt.
 *
 * USAGE EXAMPLE
 *   root -l -b -q 'summarize_incoherent_wide.C("fig2_incoherent_4bins_1M_wide.root")'
 *
 * REPOSITORY STATUS
 *   Fig.2 numerical diagnostic; executable code unchanged.
 */
#include <TFile.h>
#include <TH1D.h>
#include <TAxis.h>
#include <TString.h>

#include <fstream>
#include <iomanip>
#include <iostream>

void summarize_incoherent_wide(
    const char* inputFile = "fig2_incoherent_4bins_1M_wide.root",
    const char* outputFile = "incoherent_4bins_summary_clean.txt")
{
    TFile file(inputFile, "READ");

    if (file.IsZombie()) {
        std::cerr << "Cannot open " << inputFile << std::endl;
        return;
    }

    const int cuts[3] = {8, 10, 12};
    const double eLow[4]  = {223, 283, 319, 356};
    const double eHigh[4] = {234, 294, 330, 366};

    std::ofstream out(outputFile);

    out << "# Global: complete range -400 to 50 MeV\n";
    out << "# Fig2: window -60 to 40 MeV\n";
    out << "# E low high cut "
        << "N_before_global N_after_global survival_global_percent "
        << "N_before_Fig2 N_after_Fig2 survival_Fig2_percent\n";

    out << std::fixed << std::setprecision(6);

    for (int ie = 0; ie < 4; ++ie) {
        TH1D* hBefore =
            (TH1D*)file.Get(Form("inc_deltaE_E%d_nodphi", ie));

        if (!hBefore) {
            std::cerr << "Missing before histogram for E" << ie << std::endl;
            continue;
        }

        const int nBins = hBefore->GetNbinsX();
        const int binLow =
            hBefore->GetXaxis()->FindBin(-60.0 + 1.0e-6);
        const int binHigh =
            hBefore->GetXaxis()->FindBin(40.0 - 1.0e-6);

        const double nBeforeGlobal =
            hBefore->Integral(0, nBins + 1);
        const double nBeforeFig2 =
            hBefore->Integral(binLow, binHigh);

        for (int ic = 0; ic < 3; ++ic) {
            TH1D* hAfter =
                (TH1D*)file.Get(
                    Form("inc_deltaE_E%d_dphi%d", ie, cuts[ic]));

            if (!hAfter) {
                std::cerr << "Missing after histogram E"
                          << ie << ", cut " << cuts[ic] << std::endl;
                continue;
            }

            const double nAfterGlobal =
                hAfter->Integral(0, hAfter->GetNbinsX() + 1);
            const double nAfterFig2 =
                hAfter->Integral(binLow, binHigh);

            const double survivalGlobal =
                nBeforeGlobal > 0.0
                ? 100.0 * nAfterGlobal / nBeforeGlobal
                : 0.0;

            const double survivalFig2 =
                nBeforeFig2 > 0.0
                ? 100.0 * nAfterFig2 / nBeforeFig2
                : 0.0;

            out << ie << " "
                << eLow[ie] << " "
                << eHigh[ie] << " "
                << cuts[ic] << " "
                << nBeforeGlobal << " "
                << nAfterGlobal << " "
                << survivalGlobal << " "
                << nBeforeFig2 << " "
                << nAfterFig2 << " "
                << survivalFig2 << "\n";
        }
    }

    out.close();

    std::cout << "Written: " << outputFile << std::endl;
}
