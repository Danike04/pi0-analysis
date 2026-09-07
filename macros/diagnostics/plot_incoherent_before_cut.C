/*
 * plot_incoherent_before_cut.C
 *
 * PURPOSE
 *   Plot the incoherent ppnn missing-energy templates BEFORE the DeltaPhi cut
 *   in the four main Fig.2 energy bins. Optionally normalizes each distribution
 *   to unit area and writes the displayed histograms to a ROOT file.
 *
 * INPUT HISTOGRAMS
 *   inc_deltaE_E0_nodphi ... inc_deltaE_E3_nodphi from the standardized
 *   incoherent template file.
 *
 * DEFAULT OUTPUTS
 *   incoherent_before_deltaphi.pdf and incoherent_before_deltaphi.root.
 *
 * USAGE EXAMPLE
 *   root -l -b -q 'plot_incoherent_before_cut.C("fig2_incoherent_templates.root")'
 *
 * REPOSITORY STATUS
 *   Fig.2 shape diagnostic; executable code unchanged.
 */
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"

#include <algorithm>
#include <iostream>

namespace {

double IntegralRange(const TH1D* h, double xmin, double xmax)
{
    if (!h) return 0.0;

    int b1 = h->GetXaxis()->FindBin(xmin + 1e-9);
    int b2 = h->GetXaxis()->FindBin(xmax - 1e-9);

    b1 = std::max(1, b1);
    b2 = std::min(h->GetNbinsX(), b2);

    return b2 >= b1 ? h->Integral(b1, b2) : 0.0;
}

} // namespace

/*
 * Mostra esclusivamente il MC incoerente ppnn prima del taglio DeltaPhi.
 *
 * Richiede un file prodotto dalla v4, contenente:
 *
 *   inc_deltaE_E0_nodphi
 *   inc_deltaE_E1_nodphi
 *   inc_deltaE_E2_nodphi
 *   inc_deltaE_E3_nodphi
 */
void plot_incoherent_before_cut(
    const char* inputFile = "fig2_incoherent_templates.root",
    const char* outputPdf = "incoherent_before_deltaphi.pdf",
    const char* outputRoot = "incoherent_before_deltaphi.root",
    bool normalizeToUnitArea = false,
    double xMin = -80.0,
    double xMax = 20.0)
{
    gStyle->SetOptStat(0);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetLegendBorderSize(0);

    const int nEnergyBins = 4;

    const double energyLow[nEnergyBins] = {
        223.0, 283.0, 319.0, 356.0
    };

    const double energyHigh[nEnergyBins] = {
        234.0, 294.0, 330.0, 366.0
    };

    const char* energyLabel[nEnergyBins] = {
        "223-234 MeV",
        "283-294 MeV",
        "319-330 MeV",
        "356-366 MeV"
    };

    TFile* input = TFile::Open(inputFile, "READ");

    if (!input || input->IsZombie()) {
        std::cerr << "Cannot open " << inputFile << std::endl;
        return;
    }

    TFile output(outputRoot, "RECREATE");

    if (output.IsZombie()) {
        std::cerr << "Cannot create " << outputRoot << std::endl;
        input->Close();
        return;
    }

    TCanvas* canvas =
        new TCanvas("c_incoherent_before",
                    "Incoherent ppnn before DeltaPhi cut",
                    1200,
                    900);

    canvas->Divide(2, 2);

    for (int ie = 0; ie < nEnergyBins; ++ie) {
        TH1D* source = dynamic_cast<TH1D*>(
            input->Get(Form("inc_deltaE_E%d_nodphi", ie)));

        if (!source) {
            std::cerr
                << "Missing histogram inc_deltaE_E"
                << ie
                << "_nodphi"
                << std::endl;

            canvas->cd(ie + 1);

            TLatex missing;
            missing.SetNDC();
            missing.SetTextAlign(22);
            missing.SetTextSize(0.05);
            missing.DrawLatex(0.5, 0.5, "Missing no-cut histogram");
            continue;
        }

        TH1D* h = dynamic_cast<TH1D*>(
            source->Clone(Form("inc_before_E%d", ie)));

        h->SetDirectory(0);
        h->Sumw2();

        const double rawIntegral =
            IntegralRange(h, xMin, xMax);

        if (normalizeToUnitArea && rawIntegral > 0.0)
            h->Scale(1.0 / rawIntegral);

        const double negativeFraction =
            rawIntegral > 0.0
                ? IntegralRange(h, xMin, 0.0)
                  / IntegralRange(h, xMin, xMax)
                : 0.0;

        const double farLeftFraction =
            rawIntegral > 0.0
                ? IntegralRange(h, -60.0, -40.0)
                  / IntegralRange(h, xMin, xMax)
                : 0.0;

        canvas->cd(ie + 1);

        gPad->SetTicks(1, 1);
        gPad->SetLeftMargin(0.13);
        gPad->SetBottomMargin(0.13);
        gPad->SetRightMargin(0.04);
        gPad->SetTopMargin(0.10);

        h->SetTitle(
            Form("%s, before #Delta#Phi cut",
                 energyLabel[ie]));

        h->GetXaxis()->SetTitle(
            "#Delta E_{#pi^{0}}^{*} (MeV)");

        h->GetYaxis()->SetTitle(
            normalizeToUnitArea
                ? "Unit-area shape"
                : "MC events");

        h->GetXaxis()->SetRangeUser(xMin, xMax);

        h->SetLineColor(kBlue + 1);
        h->SetLineWidth(3);
        h->SetFillColorAlpha(kBlue + 1, 0.18);
        h->SetFillStyle(1001);

        h->SetMinimum(0.0);
        h->SetMaximum(
            h->GetMaximum() > 0.0
                ? 1.28 * h->GetMaximum()
                : 1.0);

        h->Draw("HIST");

        TLine* zero =
            new TLine(0.0, 0.0, 0.0, h->GetMaximum());

        zero->SetLineColor(kGray + 2);
        zero->SetLineStyle(2);
        zero->SetLineWidth(2);
        zero->Draw("SAME");

        TLine* minusTwenty =
            new TLine(-20.0, 0.0, -20.0, h->GetMaximum());

        minusTwenty->SetLineColor(kRed + 1);
        minusTwenty->SetLineStyle(3);
        minusTwenty->SetLineWidth(2);
        minusTwenty->Draw("SAME");

        TLegend* legend =
            new TLegend(0.48, 0.72, 0.89, 0.88);

        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->AddEntry(
            h,
            "Incoherent ppnn MC",
            "lf");
        legend->AddEntry(
            minusTwenty,
            "#Delta E=-20 MeV",
            "l");
        legend->Draw();

        TLatex text;
        text.SetNDC();
        text.SetTextSize(0.034);

        text.DrawLatex(
            0.16,
            0.88,
            Form("Entries in range: %.0f",
                 rawIntegral));

        text.DrawLatex(
            0.16,
            0.83,
            Form("Mean = %.1f MeV",
                 h->GetMean()));

        text.DrawLatex(
            0.16,
            0.78,
            Form("Fraction #Delta E<0: %.1f%%",
                 100.0 * negativeFraction));

        text.DrawLatex(
            0.16,
            0.73,
            Form("Fraction -60<#Delta E<-40: %.1f%%",
                 100.0 * farLeftFraction));

        output.cd();
        h->Write();

        delete zero;
        delete minusTwenty;
    }

    output.cd();
    canvas->Write();
    output.Write();
    output.Close();

    canvas->Print(outputPdf);

    delete canvas;
    input->Close();

    std::cout
        << "Wrote "
        << outputPdf
        << " and "
        << outputRoot
        << std::endl;
}
