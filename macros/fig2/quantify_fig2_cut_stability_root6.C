/*
 * quantify_fig2_cut_stability_root6.C
 *
 * PURPOSE
 *   Quantify how the Fig. 2 missing-energy comparison changes when the
 *   opening-angle cut is varied between DeltaPhi < 8, 10, and 12 degrees.
 *
 * INPUTS
 *   By default this macro reads the three ROOT files created by running
 *   paper_fig2_missing_energy_root6.C three times and manually renaming the
 *   outputs:
 *     paper_fig2_missing_energy_cut8.root
 *     paper_fig2_missing_energy_cut10.root
 *     paper_fig2_missing_energy_cut12.root
 *
 * WHAT THE MACRO ACTUALLY DOES
 *   - Reads data and coherent-MC histograms after each DeltaPhi cut.
 *   - Normalizes the coherent MC to the data in the configurable core region
 *     (default: -10 to +15 MeV).
 *   - Measures the negative-tail excess in the configurable tail region
 *     (default: -60 to -20 MeV).
 *   - Calculates a chi2 comparison over the configurable comparison range.
 *   - Produces overlays and a text table of the resulting metrics.
 *   - Saves the unscaled data/coherent histograms in a standardized ROOT file
 *     (default: fig2_data_coherent.root) for later template comparisons.
 *
 * OUTPUTS
 *   fig2_cut_stability_overlays.pdf
 *   fig2_cut_stability_metrics.txt
 *   fig2_data_coherent.root
 *
 * DEPENDENCIES
 *   Fig2SaveInputs.h
 *   Fig2Common.h (included by Fig2SaveInputs.h)
 *
 * USAGE
 *   root -l -b -q 'quantify_fig2_cut_stability_root6.C()'
 *
 * NOTE FOR THE REPOSITORY
 *   The analysis logic, function signature, defaults and output names below
 *   are preserved from the original macro. Documentation has been added
 *   without changing its behaviour.
 */
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

#include "Fig2SaveInputs.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"

namespace {

struct Metrics {
    double scale;
    double dataCore;
    double mcCore;
    double dataTail;
    double mcTail;
    double tailExcess;
    double tailExcessErr;
    double chi2;
    int ndf;
    double pvalue;
};

double IntegralAndErrorRange(TH1D* h,
                             double xmin,
                             double xmax,
                             double& error)
{
    if (!h) {
        error = 0.0;
        return 0.0;
    }

    int firstBin = h->GetXaxis()->FindBin(xmin + 1e-9);
    int lastBin  = h->GetXaxis()->FindBin(xmax - 1e-9);

    firstBin = std::max(1, firstBin);
    lastBin  = std::min(h->GetNbinsX(), lastBin);

    if (lastBin < firstBin) {
        error = 0.0;
        return 0.0;
    }

    return h->IntegralAndError(firstBin, lastBin, error);
}

Metrics Evaluate(TH1D* data,
                 TH1D* mc,
                 double coreMin,
                 double coreMax,
                 double tailMin,
                 double tailMax,
                 double compareMin,
                 double compareMax)
{
    Metrics result;
    result.scale = 0.0;
    result.dataCore = 0.0;
    result.mcCore = 0.0;
    result.dataTail = 0.0;
    result.mcTail = 0.0;
    result.tailExcess = 0.0;
    result.tailExcessErr = 0.0;
    result.chi2 = 0.0;
    result.ndf = 0;
    result.pvalue = 0.0;

    double dataCoreErr = 0.0;
    double mcCoreErr = 0.0;

    result.dataCore =
        IntegralAndErrorRange(data, coreMin, coreMax, dataCoreErr);
    result.mcCore =
        IntegralAndErrorRange(mc, coreMin, coreMax, mcCoreErr);

    if (result.mcCore > 0.0)
        result.scale = result.dataCore / result.mcCore;

    double dataTailErr = 0.0;
    double mcTailErrRaw = 0.0;

    result.dataTail =
        IntegralAndErrorRange(data, tailMin, tailMax, dataTailErr);

    const double mcTailRaw =
        IntegralAndErrorRange(mc, tailMin, tailMax, mcTailErrRaw);

    result.mcTail = result.scale * mcTailRaw;
    const double mcTailErr = result.scale * mcTailErrRaw;

    result.tailExcess = result.dataTail - result.mcTail;
    result.tailExcessErr =
        std::sqrt(dataTailErr * dataTailErr + mcTailErr * mcTailErr);

    int firstBin = data->GetXaxis()->FindBin(compareMin + 1e-9);
    int lastBin  = data->GetXaxis()->FindBin(compareMax - 1e-9);

    firstBin = std::max(1, firstBin);
    lastBin  = std::min(data->GetNbinsX(), lastBin);

    int usedBins = 0;

    for (int bin = firstBin; bin <= lastBin; ++bin) {
        const double observed = data->GetBinContent(bin);
        const double expected = result.scale * mc->GetBinContent(bin);

        const double observedError = data->GetBinError(bin);
        const double expectedError = result.scale * mc->GetBinError(bin);
        const double variance =
            observedError * observedError + expectedError * expectedError;

        if (variance <= 0.0)
            continue;

        const double difference = observed - expected;
        result.chi2 += difference * difference / variance;
        ++usedBins;
    }

    // One fitted parameter: the coherent normalization in the core.
    result.ndf = usedBins > 1 ? usedBins - 1 : 0;

    if (result.ndf > 0)
        result.pvalue = TMath::Prob(result.chi2, result.ndf);

    return result;
}

void DeleteSavedInputs(
    TH1* data[Fig2::kNEnergyBins][Fig2::kNCuts],
    TH1* coherent[Fig2::kNEnergyBins][Fig2::kNCuts])
{
    for (int energyIndex = 0;
         energyIndex < Fig2::kNEnergyBins;
         ++energyIndex) {
        for (int cutIndex = 0;
             cutIndex < Fig2::kNCuts;
             ++cutIndex) {
            delete data[energyIndex][cutIndex];
            delete coherent[energyIndex][cutIndex];
            data[energyIndex][cutIndex] = 0;
            coherent[energyIndex][cutIndex] = 0;
        }
    }
}

} // namespace

void quantify_fig2_cut_stability_root6(
    const char* file8 = "paper_fig2_missing_energy_cut8.root",
    const char* file10 = "paper_fig2_missing_energy_cut10.root",
    const char* file12 = "paper_fig2_missing_energy_cut12.root",
    double coreMin = -10.0,
    double coreMax = 15.0,
    double tailMin = -60.0,
    double tailMax = -20.0,
    double compareMin = -60.0,
    double compareMax = 40.0,
    const char* outputPdf = "fig2_cut_stability_overlays.pdf",
    const char* outputMetrics = "fig2_cut_stability_metrics.txt",
    const char* outputInputs = "fig2_data_coherent.root")
{
    gStyle->SetOptStat(0);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetLegendBorderSize(0);

    const int cuts[Fig2::kNCuts] = {8, 10, 12};
    const char* inputFiles[Fig2::kNCuts] = {file8, file10, file12};

    const int energyLow[Fig2::kNEnergyBins] = {223, 283, 319, 356};
    const int energyHigh[Fig2::kNEnergyBins] = {234, 294, 330, 366};

    const char* labels[Fig2::kNEnergyBins] = {
        "224 MeV",
        "294 MeV",
        "320 MeV",
        "366 MeV"
    };

    TFile* input[Fig2::kNCuts] = {0, 0, 0};

    for (int cutIndex = 0;
         cutIndex < Fig2::kNCuts;
         ++cutIndex) {
        input[cutIndex] = TFile::Open(inputFiles[cutIndex], "READ");

        if (!input[cutIndex] || input[cutIndex]->IsZombie()) {
            std::cerr << "Cannot open "
                      << inputFiles[cutIndex]
                      << std::endl;

            for (int previous = 0; previous <= cutIndex; ++previous) {
                if (input[previous])
                    input[previous]->Close();
            }
            return;
        }
    }

    TH1* savedData[Fig2::kNEnergyBins][Fig2::kNCuts] = {{0}};
    TH1* savedCoherent[Fig2::kNEnergyBins][Fig2::kNCuts] = {{0}};

    std::ofstream metrics(outputMetrics);
    if (!metrics) {
        std::cerr << "Cannot create " << outputMetrics << std::endl;

        for (int cutIndex = 0;
             cutIndex < Fig2::kNCuts;
             ++cutIndex)
            input[cutIndex]->Close();

        return;
    }

    metrics << "# Quantitative stability test of the opening-angle cut\n";
    metrics << "# DeltaE_pi0^* is the CM missing energy defined by paper Eq. (5)\n";
    metrics << "# Coherent MC normalized to data in DeltaE core = ["
            << coreMin << "," << coreMax << "] MeV\n";
    metrics << "# Negative tail = ["
            << tailMin << "," << tailMax << "] MeV\n";
    metrics << "# Chi2 comparison range = ["
            << compareMin << "," << compareMax << "] MeV\n";
    metrics << "# cut E_low E_high scale data_core mc_core_raw "
               "data_tail mc_tail_scaled tail_excess tail_excess_err "
               "tail_significance chi2 ndf chi2_ndf pvalue\n";

    TCanvas* canvas =
        new TCanvas("c_fig2_stability",
                    "Fig. 2 cut stability",
                    1200,
                    900);

    canvas->Print(Form("%s[", outputPdf));

    bool allHistogramsFound = true;

    for (int cutIndex = 0;
         cutIndex < Fig2::kNCuts;
         ++cutIndex) {
        canvas->Clear();
        canvas->Divide(2, 2);

        for (int energyIndex = 0;
             energyIndex < Fig2::kNEnergyBins;
             ++energyIndex) {
            TH1D* sourceData = dynamic_cast<TH1D*>(
                input[cutIndex]->Get(
                    Form("data_deltaE_after_phi_E%d", energyIndex)));

            TH1D* sourceCoherent = dynamic_cast<TH1D*>(
                input[cutIndex]->Get(
                    Form("coh_mc_deltaE_after_phi_E%d", energyIndex)));

            if (!sourceData || !sourceCoherent) {
                std::cerr
                    << "Missing histogram(s) in "
                    << inputFiles[cutIndex]
                    << " for energy index "
                    << energyIndex
                    << std::endl;

                allHistogramsFound = false;

                canvas->cd(energyIndex + 1);
                TLatex missing;
                missing.SetNDC();
                missing.SetTextAlign(22);
                missing.SetTextSize(0.05);
                missing.DrawLatex(0.5, 0.5, "Missing input histogram");
                continue;
            }

            // These copies remain unscaled and are written to
            // fig2_data_coherent.root for the final three-template overlay.
            savedData[energyIndex][cutIndex] =
                dynamic_cast<TH1D*>(
                    sourceData->Clone(
                        Form("saved_data_cut%d_E%d",
                             cuts[cutIndex],
                             energyIndex)));

            savedCoherent[energyIndex][cutIndex] =
                dynamic_cast<TH1D*>(
                    sourceCoherent->Clone(
                        Form("saved_coherent_cut%d_E%d",
                             cuts[cutIndex],
                             energyIndex)));

            if (!savedData[energyIndex][cutIndex] ||
                !savedCoherent[energyIndex][cutIndex]) {
                std::cerr << "Cannot clone input histogram for energy index "
                          << energyIndex
                          << ", cut "
                          << cuts[cutIndex]
                          << std::endl;
                allHistogramsFound = false;
                continue;
            }

            savedData[energyIndex][cutIndex]->SetDirectory(0);
            savedCoherent[energyIndex][cutIndex]->SetDirectory(0);

            // Separate drawing copies. Only this coherent copy is scaled.
            TH1D* data = dynamic_cast<TH1D*>(
                savedData[energyIndex][cutIndex]->Clone(
                    Form("draw_data_cut%d_E%d",
                         cuts[cutIndex],
                         energyIndex)));

            TH1D* coherent = dynamic_cast<TH1D*>(
                savedCoherent[energyIndex][cutIndex]->Clone(
                    Form("draw_coherent_cut%d_E%d",
                         cuts[cutIndex],
                         energyIndex)));

            data->SetDirectory(0);
            coherent->SetDirectory(0);

            const Metrics result =
                Evaluate(data,
                         coherent,
                         coreMin,
                         coreMax,
                         tailMin,
                         tailMax,
                         compareMin,
                         compareMax);

            if (result.scale > 0.0)
                coherent->Scale(result.scale);

            const double significance =
                result.tailExcessErr > 0.0
                    ? result.tailExcess / result.tailExcessErr
                    : 0.0;

            const double chi2Ndf =
                result.ndf > 0
                    ? result.chi2 / result.ndf
                    : 0.0;

            metrics
                << cuts[cutIndex] << " "
                << energyLow[energyIndex] << " "
                << energyHigh[energyIndex] << " "
                << result.scale << " "
                << result.dataCore << " "
                << result.mcCore << " "
                << result.dataTail << " "
                << result.mcTail << " "
                << result.tailExcess << " "
                << result.tailExcessErr << " "
                << significance << " "
                << result.chi2 << " "
                << result.ndf << " "
                << chi2Ndf << " "
                << result.pvalue
                << "\n";

            canvas->cd(energyIndex + 1);
            gPad->SetTicks(1, 1);
            gPad->SetLeftMargin(0.13);
            gPad->SetBottomMargin(0.13);

            data->SetTitle(
                Form("%s, #Delta#Phi < %d^{#circ}",
                     labels[energyIndex],
                     cuts[cutIndex]));

            data->GetXaxis()->SetTitle(
                "#Delta E_{#pi^{0}}^{*} (MeV)");
            data->GetYaxis()->SetTitle("Counts");
            data->GetXaxis()->SetRangeUser(compareMin, compareMax);

            data->SetLineColor(kBlack);
            data->SetMarkerColor(kBlack);
            data->SetMarkerStyle(20);
            data->SetMarkerSize(0.45);
            data->SetLineWidth(2);

            coherent->SetLineColor(kRed + 1);
            coherent->SetLineWidth(2);
            coherent->SetFillStyle(0);

            const double maximum =
                std::max(data->GetMaximum(),
                         coherent->GetMaximum());

            data->SetMinimum(0.0);
            data->SetMaximum(maximum > 0.0 ? 1.22 * maximum : 1.0);

            data->Draw("E1");
            coherent->Draw("HIST SAME");
            data->Draw("E1 SAME");

            TLegend* legend =
                new TLegend(0.50, 0.72, 0.89, 0.88);
            legend->SetBorderSize(0);
            legend->SetFillStyle(0);
            legend->AddEntry(data, "Data after cut", "lep");
            legend->AddEntry(
                coherent,
                "Coherent MC (core scaled)",
                "l");
            legend->Draw();

            TLatex text;
            text.SetNDC();
            text.SetTextSize(0.035);
            text.DrawLatex(
                0.16,
                0.88,
                Form("#chi^{2}/ndf = %.2f", chi2Ndf));
            text.DrawLatex(
                0.16,
                0.83,
                Form("tail excess = %.1f #sigma", significance));
        }

        canvas->Print(outputPdf);
    }

    canvas->Print(Form("%s]", outputPdf));
    metrics.close();

    bool saved = false;
    if (allHistogramsFound) {
        saved = SaveFig2Inputs(
            savedData,
            savedCoherent,
            outputInputs);
    } else {
        std::cerr
            << "The stability PDF was produced, but "
            << outputInputs
            << " was not written because at least one input histogram "
               "was missing."
            << std::endl;
    }

    for (int cutIndex = 0;
         cutIndex < Fig2::kNCuts;
         ++cutIndex)
        input[cutIndex]->Close();

    DeleteSavedInputs(savedData, savedCoherent);

    delete canvas;

    std::cout
        << "Saved "
        << outputMetrics
        << " and "
        << outputPdf
        << std::endl;

    if (saved)
        std::cout << "Saved " << outputInputs << std::endl;
}
