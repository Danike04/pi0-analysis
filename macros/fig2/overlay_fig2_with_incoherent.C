/*
 * overlay_fig2_with_incoherent.C
 *
 * Purpose
 * -------
 * Overlay standardized experimental/coherent Fig. 2 histograms with a
 * standardized incoherent Monte Carlo template for each of the four energy
 * bins and DeltaPhi cuts 8, 10 and 12 degrees.
 *
 * Inputs
 * ------
 * Default data/coherent file: fig2_data_coherent.root
 *   produced by `quantify_fig2_cut_stability_root6.C`.
 * Default incoherent file: fig2_incoherent_templates.root
 *   produced, for Acqu MC, by `make_acqu_incoherent_fig2_templates.C`.
 *
 * Normalization modes
 * -------------------
 * mode 0 (default): scale coherent MC to the data core (-10..15 MeV), then
 *                   scale the incoherent template to the remaining negative
 *                   tail (-60..-20 MeV).
 * mode 1:            non-negative simultaneous least-squares fit of the two
 *                   templates in the comparison interval (-60..40 MeV).
 *
 * Output defaults
 * ---------------
 * See the function signature below for the PDF, ROOT and metrics filenames.
 *
 * This belongs to the pre-full/empty incoherent-template study chain. Later
 * full-minus-empty refits use dedicated breakup-channel files and a separate
 * refit macro.
 *
 * Original analysis code follows unchanged; only this header was added.
 */

#include "Fig2Common.h"

#include "TArrayD.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

namespace {

TH1* GetHistogram(TFile* file,
                  const char* primary,
                  const char* alternative1 = 0,
                  const char* alternative2 = 0,
                  const char* alternative3 = 0)
{
    const char* names[4] = {
        primary,
        alternative1,
        alternative2,
        alternative3
    };

    for (int index = 0; index < 4; ++index) {
        if (!names[index] || !names[index][0])
            continue;

        TH1* histogram =
            dynamic_cast<TH1*>(file->Get(names[index]));

        if (histogram)
            return histogram;
    }

    return 0;
}

TH1D* RemapToDataBinning(const TH1* source,
                         const TH1* data,
                         const char* newName)
{
    const int numberBins = data->GetNbinsX();
    const TArrayD* variableBins =
        data->GetXaxis()->GetXbins();

    TH1D* output = 0;

    if (variableBins &&
        variableBins->GetSize() == numberBins + 1) {
        output = new TH1D(
            newName,
            source->GetTitle(),
            numberBins,
            variableBins->GetArray());
    } else {
        output = new TH1D(
            newName,
            source->GetTitle(),
            numberBins,
            data->GetXaxis()->GetXmin(),
            data->GetXaxis()->GetXmax());
    }

    output->Sumw2();
    output->SetDirectory(0);

    for (int sourceBin = 1;
         sourceBin <= source->GetNbinsX();
         ++sourceBin) {
        const double x =
            source->GetXaxis()->GetBinCenter(sourceBin);

        const int outputBin = output->FindBin(x);

        if (outputBin < 1 ||
            outputBin > output->GetNbinsX())
            continue;

        const double oldContent =
            output->GetBinContent(outputBin);
        const double oldError =
            output->GetBinError(outputBin);

        const double addedContent =
            source->GetBinContent(sourceBin);
        const double addedError =
            source->GetBinError(sourceBin);

        output->SetBinContent(
            outputBin,
            oldContent + addedContent);

        output->SetBinError(
            outputBin,
            std::sqrt(oldError * oldError +
                      addedError * addedError));
    }

    return output;
}

double IntegralRange(const TH1* histogram,
                     double minimum,
                     double maximum)
{
    int firstBin =
        histogram->GetXaxis()->FindBin(minimum + 1e-9);
    int lastBin =
        histogram->GetXaxis()->FindBin(maximum - 1e-9);

    firstBin = std::max(1, firstBin);
    lastBin =
        std::min(histogram->GetNbinsX(), lastBin);

    if (lastBin < firstBin)
        return 0.0;

    return histogram->Integral(firstBin, lastBin);
}

double IntegralErrorRange(const TH1* histogram,
                          double minimum,
                          double maximum)
{
    int firstBin =
        histogram->GetXaxis()->FindBin(minimum + 1e-9);
    int lastBin =
        histogram->GetXaxis()->FindBin(maximum - 1e-9);

    firstBin = std::max(1, firstBin);
    lastBin =
        std::min(histogram->GetNbinsX(), lastBin);

    double variance = 0.0;

    for (int bin = firstBin; bin <= lastBin; ++bin) {
        const double error = histogram->GetBinError(bin);
        variance += error * error;
    }

    return std::sqrt(variance);
}

double ScaleCoherentInCore(const TH1* data,
                           const TH1* coherent,
                           double coreMin,
                           double coreMax)
{
    const double dataIntegral =
        IntegralRange(data, coreMin, coreMax);

    const double coherentIntegral =
        IntegralRange(coherent, coreMin, coreMax);

    return coherentIntegral > 0.0
        ? dataIntegral / coherentIntegral
        : 0.0;
}

double ScaleIncoherentToTailResidual(
    const TH1* data,
    const TH1* coherentScaled,
    const TH1* incoherent,
    double tailMin,
    double tailMax)
{
    const double residual =
        IntegralRange(data, tailMin, tailMax)
        - IntegralRange(
            coherentScaled,
            tailMin,
            tailMax);

    const double incoherentIntegral =
        IntegralRange(
            incoherent,
            tailMin,
            tailMax);

    if (incoherentIntegral <= 0.0 ||
        residual <= 0.0)
        return 0.0;

    return residual / incoherentIntegral;
}

void FitTwoNonNegativeTemplates(
    const TH1* data,
    const TH1* coherent,
    const TH1* incoherent,
    double fitMin,
    double fitMax,
    double& coherentScale,
    double& incoherentScale)
{
    double coherentCoherent = 0.0;
    double incoherentIncoherent = 0.0;
    double coherentIncoherent = 0.0;
    double coherentData = 0.0;
    double incoherentData = 0.0;

    for (int bin = 1;
         bin <= data->GetNbinsX();
         ++bin) {
        const double x =
            data->GetXaxis()->GetBinCenter(bin);

        if (x < fitMin || x > fitMax)
            continue;

        const double observed =
            data->GetBinContent(bin);

        const double coherentValue =
            coherent->GetBinContent(bin);

        const double incoherentValue =
            incoherent->GetBinContent(bin);

        double error = data->GetBinError(bin);
        if (error <= 0.0)
            error = std::sqrt(std::max(1.0, observed));

        const double weight = 1.0 / (error * error);

        coherentCoherent +=
            weight * coherentValue * coherentValue;

        incoherentIncoherent +=
            weight * incoherentValue * incoherentValue;

        coherentIncoherent +=
            weight * coherentValue * incoherentValue;

        coherentData +=
            weight * coherentValue * observed;

        incoherentData +=
            weight * incoherentValue * observed;
    }

    const double determinant =
        coherentCoherent * incoherentIncoherent
        - coherentIncoherent * coherentIncoherent;

    if (determinant > 0.0) {
        coherentScale =
            (coherentData * incoherentIncoherent
             - incoherentData * coherentIncoherent)
            / determinant;

        incoherentScale =
            (incoherentData * coherentCoherent
             - coherentData * coherentIncoherent)
            / determinant;
    } else {
        coherentScale = 0.0;
        incoherentScale = 0.0;
    }

    if (coherentScale < 0.0) {
        coherentScale = 0.0;
        incoherentScale =
            incoherentIncoherent > 0.0
                ? incoherentData / incoherentIncoherent
                : 0.0;
    }

    if (incoherentScale < 0.0) {
        incoherentScale = 0.0;
        coherentScale =
            coherentCoherent > 0.0
                ? coherentData / coherentCoherent
                : 0.0;
    }

    coherentScale = std::max(0.0, coherentScale);
    incoherentScale = std::max(0.0, incoherentScale);
}

double Chi2Ndf(const TH1* data,
               const TH1* model,
               double minimum,
               double maximum,
               int fittedParameters,
               int& ndf)
{
    double chi2 = 0.0;
    int usedBins = 0;

    for (int bin = 1;
         bin <= data->GetNbinsX();
         ++bin) {
        const double x =
            data->GetXaxis()->GetBinCenter(bin);

        if (x < minimum || x > maximum)
            continue;

        const double observed =
            data->GetBinContent(bin);

        const double expected =
            model->GetBinContent(bin);

        double variance =
            data->GetBinError(bin) *
            data->GetBinError(bin);

        variance +=
            model->GetBinError(bin) *
            model->GetBinError(bin);

        if (variance <= 0.0)
            variance = std::max(1.0, observed);

        const double difference = observed - expected;
        chi2 += difference * difference / variance;
        ++usedBins;
    }

    ndf = usedBins - fittedParameters;

    return ndf > 0 ? chi2 / ndf : 0.0;
}

double CoherentOnlyTailExcessSignificance(
    const TH1* data,
    const TH1* coherentScaled,
    double tailMin,
    double tailMax)
{
    const double excess =
        IntegralRange(data, tailMin, tailMax)
        - IntegralRange(
            coherentScaled,
            tailMin,
            tailMax);

    const double dataError =
        IntegralErrorRange(data, tailMin, tailMax);

    const double coherentError =
        IntegralErrorRange(
            coherentScaled,
            tailMin,
            tailMax);

    const double totalError =
        std::sqrt(dataError * dataError +
                  coherentError * coherentError);

    return totalError > 0.0
        ? excess / totalError
        : 0.0;
}

void StyleData(TH1* histogram)
{
    histogram->SetMarkerStyle(20);
    histogram->SetMarkerSize(0.45);
    histogram->SetMarkerColor(kBlack);
    histogram->SetLineColor(kBlack);
    histogram->SetLineWidth(2);
}

void StyleCoherent(TH1* histogram)
{
    histogram->SetLineColor(kRed + 1);
    histogram->SetLineWidth(2);
    histogram->SetFillStyle(0);
}

void StyleIncoherent(TH1* histogram)
{
    histogram->SetLineColor(kBlue + 1);
    histogram->SetLineWidth(2);
    histogram->SetFillStyle(0);
}

void StyleSum(TH1* histogram)
{
    histogram->SetLineColor(kGreen + 2);
    histogram->SetLineWidth(3);
    histogram->SetLineStyle(2);
    histogram->SetFillStyle(0);
}

} // namespace

/*
 * normalizationMode:
 *
 * 0:
 *   Coherent MC is normalized in the core.
 *   Incoherent MC is normalized to the remaining negative-tail yield.
 *
 * 1:
 *   Non-negative simultaneous two-template least-squares fit.
 */
void overlay_fig2_with_incoherent(
    const char* dataCoherentFile =
        "fig2_data_coherent.root",
    const char* incoherentFile =
        "fig2_incoherent_templates.root",
    const char* outputPdf =
        "fig2_with_incoherent.pdf",
    const char* outputRoot =
        "fig2_with_incoherent.root",
    int normalizationMode = 0,
    double coreMin = -10.0,
    double coreMax = 15.0,
    double tailMin = -60.0,
    double tailMax = -20.0,
    double fitMin = -60.0,
    double fitMax = 40.0,
    const char* outputMetrics =
        "fig2_with_incoherent_metrics.txt")
{
    gStyle->SetOptStat(0);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetLegendBorderSize(0);

    TFile* dataCoherent =
        TFile::Open(dataCoherentFile, "READ");

    TFile* incoherentInput =
        TFile::Open(incoherentFile, "READ");

    if (!dataCoherent ||
        dataCoherent->IsZombie()) {
        std::cerr
            << "Cannot open "
            << dataCoherentFile
            << std::endl;
        return;
    }

    if (!incoherentInput ||
        incoherentInput->IsZombie()) {
        std::cerr
            << "Cannot open "
            << incoherentFile
            << std::endl;
        dataCoherent->Close();
        return;
    }

    TFile output(outputRoot, "RECREATE");
    if (output.IsZombie()) {
        std::cerr
            << "Cannot create "
            << outputRoot
            << std::endl;

        dataCoherent->Close();
        incoherentInput->Close();
        return;
    }

    std::ofstream metrics(outputMetrics);
    if (!metrics) {
        std::cerr
            << "Cannot create "
            << outputMetrics
            << std::endl;

        output.Close();
        dataCoherent->Close();
        incoherentInput->Close();
        return;
    }

    metrics
        << "# data + coherent MC + incoherent MC comparison\n"
        << "# mode " << normalizationMode << "\n"
        << "# core [" << coreMin << "," << coreMax << "] MeV\n"
        << "# tail [" << tailMin << "," << tailMax << "] MeV\n"
        << "# comparison [" << fitMin << "," << fitMax << "] MeV\n"
        << "# cut energy_index E_low E_high coherent_scale "
           "incoherent_scale coherent_only_tail_sigma "
           "sum_chi2_ndf ndf\n";

    TString pdfName(outputPdf);

    TCanvas* canvas =
        new TCanvas(
            "c_fig2_with_incoherent",
            "Fig. 2 with incoherent MC",
            1200,
            900);

    canvas->Print(pdfName + "[");

    for (int cutIndex = 0;
         cutIndex < Fig2::kNCuts;
         ++cutIndex) {
        const int cut =
            static_cast<int>(
                Fig2::kDeltaPhiCutDeg[cutIndex]);

        canvas->Clear();
        canvas->Divide(2, 2);

        for (int energyIndex = 0;
             energyIndex < Fig2::kNEnergyBins;
             ++energyIndex) {
            canvas->cd(energyIndex + 1);

            gPad->SetTicks(1, 1);
            gPad->SetLeftMargin(0.13);
            gPad->SetBottomMargin(0.13);

            TH1* dataSource =
                GetHistogram(
                    dataCoherent,
                    Form("data_E%d_dphi%d",
                         energyIndex,
                         cut));

            TH1* coherentSource =
                GetHistogram(
                    dataCoherent,
                    Form("coh_E%d_dphi%d",
                         energyIndex,
                         cut));

            TH1* incoherentSource =
                GetHistogram(
                    incoherentInput,
                    Form("inc_deltaE_E%d_dphi%d",
                         energyIndex,
                         cut));

            if (!dataSource ||
                !coherentSource ||
                !incoherentSource) {
                std::cerr
                    << "Missing histogram for energy index "
                    << energyIndex
                    << ", DeltaPhi cut "
                    << cut
                    << std::endl;

                TLatex missing;
                missing.SetNDC();
                missing.SetTextAlign(22);
                missing.SetTextSize(0.05);
                missing.DrawLatex(
                    0.5,
                    0.5,
                    "Missing input histogram");
                continue;
            }

            TH1* data =
                dynamic_cast<TH1*>(
                    dataSource->Clone(
                        Form("data_scaled_E%d_dphi%d",
                             energyIndex,
                             cut)));

            TH1D* coherent =
                RemapToDataBinning(
                    coherentSource,
                    data,
                    Form("coh_scaled_E%d_dphi%d",
                         energyIndex,
                         cut));

            TH1D* incoherent =
                RemapToDataBinning(
                    incoherentSource,
                    data,
                    Form("inc_scaled_E%d_dphi%d",
                         energyIndex,
                         cut));

            data->SetDirectory(0);
            data->Sumw2();

            double coherentScale = 0.0;
            double incoherentScale = 0.0;

            if (normalizationMode == 0) {
                coherentScale =
                    ScaleCoherentInCore(
                        data,
                        coherent,
                        coreMin,
                        coreMax);

                coherent->Scale(coherentScale);

                incoherentScale =
                    ScaleIncoherentToTailResidual(
                        data,
                        coherent,
                        incoherent,
                        tailMin,
                        tailMax);

                incoherent->Scale(incoherentScale);
            } else {
                FitTwoNonNegativeTemplates(
                    data,
                    coherent,
                    incoherent,
                    fitMin,
                    fitMax,
                    coherentScale,
                    incoherentScale);

                coherent->Scale(coherentScale);
                incoherent->Scale(incoherentScale);
            }

            TH1* total =
                dynamic_cast<TH1*>(
                    coherent->Clone(
                        Form("sum_scaled_E%d_dphi%d",
                             energyIndex,
                             cut)));

            total->SetDirectory(0);
            total->Add(incoherent);

            StyleData(data);
            StyleCoherent(coherent);
            StyleIncoherent(incoherent);
            StyleSum(total);

            data->SetTitle(
                Form("%.0f MeV, #Delta#Phi < %d^{#circ}",
                     Fig2::kEnergyCenter[energyIndex],
                     cut));

            data->GetXaxis()->SetTitle(
                "#Delta E_{#pi^{0}}^{*} (MeV)");
            data->GetYaxis()->SetTitle("Counts");
            data->GetXaxis()->SetRangeUser(fitMin, fitMax);

            const double maximum =
                std::max(
                    data->GetMaximum(),
                    total->GetMaximum());

            data->SetMinimum(0.0);
            data->SetMaximum(
                maximum > 0.0
                    ? 1.24 * maximum
                    : 1.0);

            data->Draw("E1");
            coherent->Draw("HIST SAME");
            incoherent->Draw("HIST SAME");
            total->Draw("HIST SAME");
            data->Draw("E1 SAME");

            int ndf = 0;
            const double chi2Ndf =
                Chi2Ndf(
                    data,
                    total,
                    fitMin,
                    fitMax,
                    2,
                    ndf);

            const double tailSignificance =
                CoherentOnlyTailExcessSignificance(
                    data,
                    coherent,
                    tailMin,
                    tailMax);

            TLatex text;
            text.SetNDC();
            text.SetTextSize(0.035);
            text.DrawLatex(
                0.16,
                0.88,
                Form("#chi^{2}/ndf = %.2f",
                     chi2Ndf));
            text.DrawLatex(
                0.16,
                0.83,
                Form("coherent-only tail = %.1f #sigma",
                     tailSignificance));

            TLegend* legend =
                new TLegend(0.47, 0.64, 0.90, 0.88);
            legend->SetBorderSize(0);
            legend->SetFillStyle(0);
            legend->AddEntry(
                data,
                "Data after cut",
                "lep");
            legend->AddEntry(
                coherent,
                "Coherent MC",
                "l");
            legend->AddEntry(
                incoherent,
                "Incoherent MC",
                "l");
            legend->AddEntry(
                total,
                "Coherent + incoherent",
                "l");
            legend->Draw();

            metrics
                << cut << " "
                << energyIndex << " "
                << Fig2::kEnergyLow[energyIndex] << " "
                << Fig2::kEnergyHigh[energyIndex] << " "
                << coherentScale << " "
                << incoherentScale << " "
                << tailSignificance << " "
                << chi2Ndf << " "
                << ndf
                << "\n";

            output.cd();
            data->Write();
            coherent->Write();
            incoherent->Write();
            total->Write();
        }

        output.cd();
        canvas->Write(
            Form("canvas_dphi%d", cut));

        canvas->Print(pdfName);
    }

    canvas->Print(pdfName + "]");

    metrics.close();
    output.Write();
    output.Close();

    dataCoherent->Close();
    incoherentInput->Close();

    delete canvas;

    std::cout
        << "Wrote "
        << outputPdf
        << ", "
        << outputRoot
        << " and "
        << outputMetrics
        << std::endl;
}
