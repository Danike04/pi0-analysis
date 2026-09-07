/*
 * diagnose_incoherent_selection.C
 *
 * PURPOSE
 *   Quantify how the DeltaPhi selections affect the standardized incoherent
 *   Fig.2 templates. For each of the four Fig.2 energy bins it compares the
 *   no-cut distribution with the 8/10/12 degree selections and writes survival
 *   fractions, tail/core fractions, mean and RMS.
 *
 * INPUT
 *   Standardized template file produced by make_acqu_incoherent_fig2_templates.C
 *   (default `fig2_incoherent_templates.root`). Requires Fig2Common.h.
 *
 * OUTPUTS
 *   incoherent_before_after_deltaphi.pdf and
 *   incoherent_before_after_deltaphi.txt by default.
 *
 * USAGE EXAMPLE
 *   root -l -b -q 'diagnose_incoherent_selection.C()'
 *
 * REPOSITORY STATUS
 *   Fig.2 selection diagnostic; executable code unchanged.
 */
#include "Fig2Common.h"

#include "TCanvas.h"
#include "TFile.h"
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

double IntegralRange(const TH1D* h, double xmin, double xmax)
{
    if (!h) return 0.0;
    int b1 = std::max(1, h->GetXaxis()->FindBin(xmin + 1e-9));
    int b2 = std::min(h->GetNbinsX(), h->GetXaxis()->FindBin(xmax - 1e-9));
    return b2 >= b1 ? h->Integral(b1, b2) : 0.0;
}

TH1D* CloneDetached(TH1D* source, const char* name)
{
    if (!source) return 0;
    TH1D* h = dynamic_cast<TH1D*>(source->Clone(name));
    if (h) h->SetDirectory(0);
    return h;
}

void NormalizeShape(TH1D* h)
{
    const double integral = IntegralRange(h, -60.0, 40.0);
    if (integral > 0.0) h->Scale(1.0 / integral);
}

void StyleNoCut(TH1D* h)
{
    h->SetLineColor(kBlack);
    h->SetLineWidth(3);
    h->SetFillStyle(0);
}

void StyleCut(TH1D* h, int index)
{
    const int colors[Fig2::kNCuts] = {kRed + 1, kBlue + 1, kGreen + 2};
    const int styles[Fig2::kNCuts] = {1, 2, 7};
    h->SetLineColor(colors[index]);
    h->SetLineWidth(2);
    h->SetLineStyle(styles[index]);
    h->SetFillStyle(0);
}

void WriteMetrics(std::ofstream& out,
                  int energyIndex,
                  const char* sample,
                  const TH1D* h,
                  double noCutIntegral)
{
    const double total = IntegralRange(h, -60.0, 40.0);
    const double farLeft = IntegralRange(h, -60.0, -20.0);
    const double negative = IntegralRange(h, -60.0, 0.0);
    const double core = IntegralRange(h, -10.0, 15.0);
    const double positive = IntegralRange(h, 0.0, 40.0);

    out << energyIndex << " "
        << Fig2::kEnergyLow[energyIndex] << " "
        << Fig2::kEnergyHigh[energyIndex] << " "
        << sample << " "
        << total << " "
        << (noCutIntegral > 0.0 ? total / noCutIntegral : 0.0) << " "
        << (total > 0.0 ? farLeft / total : 0.0) << " "
        << (total > 0.0 ? negative / total : 0.0) << " "
        << (total > 0.0 ? core / total : 0.0) << " "
        << (total > 0.0 ? positive / total : 0.0) << " "
        << h->GetMean() << " "
        << h->GetRMS()
        << "\n";
}

} // namespace

void diagnose_incoherent_selection(
    const char* inputFile = "fig2_incoherent_templates.root",
    const char* outputPdf = "incoherent_before_after_deltaphi.pdf",
    const char* outputMetrics = "incoherent_before_after_deltaphi.txt")
{
    gStyle->SetOptStat(0);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetLegendBorderSize(0);

    TFile* input = TFile::Open(inputFile, "READ");
    if (!input || input->IsZombie()) {
        std::cerr << "Cannot open " << inputFile << std::endl;
        return;
    }

    std::ofstream out(outputMetrics);
    if (!out) {
        std::cerr << "Cannot create " << outputMetrics << std::endl;
        input->Close();
        return;
    }

    out << "# Incoherent MC before and after the DeltaPhi selection\n";
    out << "# E_index E_low E_high sample integral survival "
           "fraction_DeltaE_lt_minus20 fraction_DeltaE_lt_0 "
           "fraction_core_minus10_to_15 fraction_DeltaE_gt_0 mean rms\n";

    TH1D* noCut[Fig2::kNEnergyBins] = {0};
    TH1D* cut[Fig2::kNEnergyBins][Fig2::kNCuts] = {{0}};

    for (int ie = 0; ie < Fig2::kNEnergyBins; ++ie) {
        noCut[ie] = CloneDetached(
            dynamic_cast<TH1D*>(input->Get(Form("inc_deltaE_E%d_nodphi", ie))),
            Form("diagnostic_nocut_E%d", ie));

        if (!noCut[ie]) {
            std::cerr << "Missing inc_deltaE_E" << ie
                      << "_nodphi. Rebuild templates with v4." << std::endl;
            input->Close();
            return;
        }

        for (int ic = 0; ic < Fig2::kNCuts; ++ic) {
            const int dphi = static_cast<int>(Fig2::kDeltaPhiCutDeg[ic]);
            cut[ie][ic] = CloneDetached(
                dynamic_cast<TH1D*>(
                    input->Get(Form("inc_deltaE_E%d_dphi%d", ie, dphi))),
                Form("diagnostic_cut_E%d_dphi%d", ie, dphi));

            if (!cut[ie][ic]) {
                std::cerr << "Missing cut histogram for E" << ie
                          << ", DeltaPhi " << dphi << std::endl;
                input->Close();
                return;
            }
        }
    }

    TCanvas* canvas =
        new TCanvas("c_incoherent_diagnostic",
                    "Incoherent MC selection diagnostic",
                    1200, 900);

    TString pdf(outputPdf);
    canvas->Print(pdf + "[");

    // Page 1: absolute yields.
    canvas->Clear();
    canvas->Divide(2, 2);

    for (int ie = 0; ie < Fig2::kNEnergyBins; ++ie) {
        canvas->cd(ie + 1);
        gPad->SetTicks(1, 1);
        gPad->SetLeftMargin(0.13);
        gPad->SetBottomMargin(0.13);

        StyleNoCut(noCut[ie]);
        for (int ic = 0; ic < Fig2::kNCuts; ++ic)
            StyleCut(cut[ie][ic], ic);

        noCut[ie]->SetTitle(
            Form("%.0f MeV: ppnn before/after #Delta#Phi cut",
                 Fig2::kEnergyCenter[ie]));
        noCut[ie]->GetXaxis()->SetTitle("#Delta E_{#pi^{0}}^{*} (MeV)");
        noCut[ie]->GetYaxis()->SetTitle("MC events");
        noCut[ie]->GetXaxis()->SetRangeUser(-60.0, 40.0);
        noCut[ie]->SetMinimum(0.0);

        double ymax = noCut[ie]->GetMaximum();
        for (int ic = 0; ic < Fig2::kNCuts; ++ic)
            ymax = std::max(ymax, cut[ie][ic]->GetMaximum());
        noCut[ie]->SetMaximum(ymax > 0.0 ? 1.25 * ymax : 1.0);

        noCut[ie]->Draw("HIST");
        for (int ic = 0; ic < Fig2::kNCuts; ++ic)
            cut[ie][ic]->Draw("HIST SAME");

        TLegend* legend = new TLegend(0.48, 0.64, 0.89, 0.88);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->AddEntry(noCut[ie], "No #Delta#Phi cut", "l");
        for (int ic = 0; ic < Fig2::kNCuts; ++ic)
            legend->AddEntry(cut[ie][ic],
                             Form("#Delta#Phi < %.0f^{#circ}",
                                  Fig2::kDeltaPhiCutDeg[ic]),
                             "l");
        legend->Draw();

        const double base = IntegralRange(noCut[ie], -60.0, 40.0);

        TLatex text;
        text.SetNDC();
        text.SetTextSize(0.033);
        for (int ic = 0; ic < Fig2::kNCuts; ++ic) {
            const double selected =
                IntegralRange(cut[ie][ic], -60.0, 40.0);
            const double survival =
                base > 0.0 ? 100.0 * selected / base : 0.0;
            text.DrawLatex(0.16, 0.88 - 0.045 * ic,
                           Form("survival %.0f^{#circ}: %.1f%%",
                                Fig2::kDeltaPhiCutDeg[ic], survival));
        }

        WriteMetrics(out, ie, "no_cut", noCut[ie], base);
        for (int ic = 0; ic < Fig2::kNCuts; ++ic)
            WriteMetrics(out, ie,
                         Form("dphi%.0f", Fig2::kDeltaPhiCutDeg[ic]),
                         cut[ie][ic], base);
    }
    canvas->Print(pdf);

    // Page 2: normalized shapes.
    canvas->Clear();
    canvas->Divide(2, 2);

    for (int ie = 0; ie < Fig2::kNEnergyBins; ++ie) {
        canvas->cd(ie + 1);
        gPad->SetTicks(1, 1);
        gPad->SetLeftMargin(0.13);
        gPad->SetBottomMargin(0.13);

        TH1D* noCutShape =
            CloneDetached(noCut[ie], Form("shape_nocut_E%d", ie));
        TH1D* cutShape[Fig2::kNCuts] = {0};

        for (int ic = 0; ic < Fig2::kNCuts; ++ic)
            cutShape[ic] =
                CloneDetached(cut[ie][ic],
                              Form("shape_cut_E%d_%d", ie, ic));

        NormalizeShape(noCutShape);
        for (int ic = 0; ic < Fig2::kNCuts; ++ic)
            NormalizeShape(cutShape[ic]);

        StyleNoCut(noCutShape);
        for (int ic = 0; ic < Fig2::kNCuts; ++ic)
            StyleCut(cutShape[ic], ic);

        noCutShape->SetTitle(
            Form("%.0f MeV: normalized ppnn shapes",
                 Fig2::kEnergyCenter[ie]));
        noCutShape->GetXaxis()->SetTitle("#Delta E_{#pi^{0}}^{*} (MeV)");
        noCutShape->GetYaxis()->SetTitle("Unit-area shape");
        noCutShape->GetXaxis()->SetRangeUser(-60.0, 40.0);
        noCutShape->SetMinimum(0.0);

        double ymax = noCutShape->GetMaximum();
        for (int ic = 0; ic < Fig2::kNCuts; ++ic)
            ymax = std::max(ymax, cutShape[ic]->GetMaximum());
        noCutShape->SetMaximum(ymax > 0.0 ? 1.25 * ymax : 1.0);

        noCutShape->Draw("HIST");
        for (int ic = 0; ic < Fig2::kNCuts; ++ic)
            cutShape[ic]->Draw("HIST SAME");

        TLegend* legend = new TLegend(0.48, 0.64, 0.89, 0.88);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->AddEntry(noCutShape, "No #Delta#Phi cut", "l");
        for (int ic = 0; ic < Fig2::kNCuts; ++ic)
            legend->AddEntry(cutShape[ic],
                             Form("#Delta#Phi < %.0f^{#circ}",
                                  Fig2::kDeltaPhiCutDeg[ic]),
                             "l");
        legend->Draw();

        delete noCutShape;
        for (int ic = 0; ic < Fig2::kNCuts; ++ic)
            delete cutShape[ic];
    }

    canvas->Print(pdf);
    canvas->Print(pdf + "]");

    delete canvas;
    for (int ie = 0; ie < Fig2::kNEnergyBins; ++ie) {
        delete noCut[ie];
        for (int ic = 0; ic < Fig2::kNCuts; ++ic)
            delete cut[ie][ic];
    }

    out.close();
    input->Close();

    std::cout << "Wrote " << outputPdf
              << " and " << outputMetrics << std::endl;
}
