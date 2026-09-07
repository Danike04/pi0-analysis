/*
 * ps_tagging_eff_clean.C
 *
 * Purpose
 * -------
 * Calculate the channel-by-channel tagging efficiency measured with the
 * Pair Spectrometer (PS) from an ARHist ROOT file.
 *
 * Physics quantity
 * ----------------
 * For each FPD channel the macro evaluates
 *
 *     epsilon_tag^PS = correctionFactor * (B - C) / A,
 *
 * where
 *   A = FPD_ScalerAcc,
 *   B = PairSpec_SumGated,
 *   C = PairSpec_SumGatedDly.
 *
 * The default correction factor is 2500, as used in the 2026 analysis.
 * Statistical uncertainties are propagated treating A, B and C as
 * independent Poisson counts.
 *
 * Input
 * -----
 * A ROOT histogram file (typically ARHist_CBTagg_<run>.root) containing
 * the three histograms above. Histogram lookup is recursive, so the
 * histograms may live inside ROOT subdirectories.
 *
 * Output
 * ------
 * <outputBase>.txt   channel-by-channel numerical results
 * <outputBase>.root  input histograms plus derived efficiency histogram/graph
 * <outputBase>.pdf   plot
 * <outputBase>.png   plot
 *
 * Usage
 * -----
 * root -l -b -q 'ps_tagging_eff_clean.C("ARHist_CBTagg_RUN.root","PS_TaggEff_RUN",2500.0)'
 *
 * Notes
 * -----
 * ROOT bin 1 is written as FPD channel 0 in the text output.
 * The 20--25% band drawn on the plot is diagnostic only; it is not used
 * to accept or reject channels in the numerical output.
 *
 * The code below is the original analysis code. This header only documents
 * its behaviour; the calculation, defaults and outputs are unchanged.
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>

#include "TFile.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TClass.h"
#include "TH1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TBox.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "TString.h"

// Cerca ricorsivamente un oggetto con nome esatto in tutte le directory ROOT.
TObject* FindObjectRecursive(TDirectory* dir, const TString& wanted)
{
    if (!dir) return nullptr;

    TObject* direct = dir->Get(wanted);
    if (direct) return direct;

    TIter next(dir->GetListOfKeys());
    TKey* key = nullptr;

    while ((key = (TKey*)next())) {
        TString keyName = key->GetName();
        TString className = key->GetClassName();

        if (keyName == wanted)
            return key->ReadObj();

        TClass* cl = TClass::GetClass(className);
        if (cl && cl->InheritsFrom(TDirectory::Class())) {
            TDirectory* subdir = dynamic_cast<TDirectory*>(key->ReadObj());
            TObject* found = FindObjectRecursive(subdir, wanted);
            if (found) return found;
        }
    }

    return nullptr;
}

TH1* FindHistogram(TDirectory* dir,
                   const std::vector<TString>& candidateNames,
                   TString& matchedName)
{
    for (const TString& name : candidateNames) {
        TObject* obj = FindObjectRecursive(dir, name);
        TH1* hist = dynamic_cast<TH1*>(obj);

        if (hist) {
            matchedName = name;
            return hist;
        }
    }

    return nullptr;
}

void ps_tagging_eff_clean(const char* inputFile = "ARHist_CBTagg_RUN.root",
                    const char* outputBase = "PS_TaggEff_RUN",
                    double correctionFactor = 2500.0)
{
    // Stile pulito da presentazione/pubblicazione.
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetCanvasColor(kWhite);
    gStyle->SetPadColor(kWhite);
    gStyle->SetFrameFillColor(kWhite);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetLineWidth(2);
    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelSize(0.040, "XYZ");
    gStyle->SetTitleSize(0.050, "XYZ");
    gStyle->SetTitleOffset(1.05, "X");
    gStyle->SetTitleOffset(1.20, "Y");
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetEndErrorSize(2);
    TGaxis::SetMaxDigits(4);

    TFile* fin = TFile::Open(inputFile, "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "ERROR: cannot open " << inputFile << std::endl;
        return;
    }

    TString nameA, nameB, nameC;

    TH1* hA = FindHistogram(
        fin,
        {"FPD_ScalerAcc",
         "FPD_scalerAcc"},
        nameA
    );

    TH1* hB = FindHistogram(
        fin,
        {"PairSpec_SumGated",
         "Pairspec_SumGated",
         "PairSpec_Sumgated",
         "Pairspec_Sumgated"},
        nameB
    );

    TH1* hC = FindHistogram(
        fin,
        {"PairSpec_SumGatedDly",
         "Pairspec_SumGatedDly",
         "PairSpec_SumgatedDly",
         "Pairspec_SumgatedDly"},
        nameC
    );

    if (!hA || !hB || !hC) {
        std::cerr << "\nERROR: one or more histograms are missing.\n"
                  << "A FPD_ScalerAcc        : " << (hA ? "FOUND" : "MISSING") << "\n"
                  << "B PairSpec_SumGated    : " << (hB ? "FOUND" : "MISSING") << "\n"
                  << "C PairSpec_SumGatedDly : " << (hC ? "FOUND" : "MISSING") << "\n\n"
                  << "Open the file interactively and inspect it with .ls,\n"
                  << "or run: rootls -r " << inputFile << std::endl;
        fin->Close();
        return;
    }

    std::cout << "Found:\n"
              << " A = " << nameA << "  bins=" << hA->GetNbinsX() << "\n"
              << " B = " << nameB << "  bins=" << hB->GetNbinsX() << "\n"
              << " C = " << nameC << "  bins=" << hC->GetNbinsX() << "\n";

    const int nBins = std::min({hA->GetNbinsX(),
                                hB->GetNbinsX(),
                                hC->GetNbinsX()});

    if (hA->GetNbinsX() != hB->GetNbinsX() ||
        hA->GetNbinsX() != hC->GetNbinsX()) {
        std::cerr << "WARNING: histogram bin numbers differ. "
                  << "Using the first " << nBins << " bins only.\n";
    }

    const TArrayD* bins = hA->GetXaxis()->GetXbins();
    TH1D* hR = nullptr;

    if (bins && bins->GetSize() == nBins + 1) {
        hR = new TH1D("PS_TaggEff",
                      "Tagging efficiency from Pair Spectrometer;"
                      "FPD channel;#epsilon_{tag}^{PS}",
                      nBins, bins->GetArray());
    } else {
        hR = new TH1D("PS_TaggEff",
                      "Tagging efficiency from Pair Spectrometer;"
                      "FPD channel;#epsilon_{tag}^{PS}",
                      nBins,
                      hA->GetXaxis()->GetXmin(),
                      hA->GetXaxis()->GetXmax());
    }

    hR->SetDirectory(nullptr);
    hR->SetMarkerStyle(20);
    hR->SetMarkerSize(0.55);
    hR->SetLineWidth(1);

    std::ofstream txt(TString::Format("%s.txt", outputBase));
    txt << "# channel A_FPD_ScalerAcc B_PairSpec_SumGated "
           "C_PairSpec_SumGatedDly BminusC R dR\n";
    txt << "# R = " << correctionFactor << " * (B-C) / A\n";

    std::vector<double> sensibleValues;
    int nValid = 0;
    int nExpectedBand = 0;
    int nZeroA = 0;

    for (int bin = 1; bin <= nBins; ++bin) {
        const double A = hA->GetBinContent(bin);
        const double B = hB->GetBinContent(bin);
        const double C = hC->GetBinContent(bin);
        const double signal = B - C;

        double R = 0.0;
        double dR = 0.0;

        if (A > 0.0) {
            R = correctionFactor * signal / A;

            // B e C trattati come conteggi poissoniani indipendenti.
            // L'errore statistico di A è incluso come termine poissoniano.
            const double varianceSignal = std::max(0.0, B) + std::max(0.0, C);
            const double varianceA = A;

            const double termSignal =
                correctionFactor * correctionFactor * varianceSignal / (A * A);

            const double termA =
                correctionFactor * correctionFactor *
                signal * signal * varianceA / (A * A * A * A);

            dR = std::sqrt(termSignal + termA);
            ++nValid;

            if (R >= 0.20 && R <= 0.25)
                ++nExpectedBand;

            // Intervallo largo usato solo per una mediana diagnostica robusta.
            if (R > 0.0 && R < 0.60)
                sensibleValues.push_back(R);
        } else {
            ++nZeroA;
        }

        hR->SetBinContent(bin, R);
        hR->SetBinError(bin, dR);

        const char* label = hA->GetXaxis()->GetBinLabel(bin);
        if (label && label[0] != '\0')
            hR->GetXaxis()->SetBinLabel(bin, label);

        // Canale ROOT: bin 1 -> canale 0.
        txt << (bin - 1) << " "
            << A << " " << B << " " << C << " "
            << signal << " " << R << " " << dR << "\n";
    }

    double median = 0.0;
    if (!sensibleValues.empty()) {
        std::sort(sensibleValues.begin(), sensibleValues.end());
        const size_t n = sensibleValues.size();
        median = (n % 2)
            ? sensibleValues[n / 2]
            : 0.5 * (sensibleValues[n / 2 - 1] +
                     sensibleValues[n / 2]);
    }

    std::cout << "\nValid channels (A>0): " << nValid << "\n"
              << "Channels with A=0:       " << nZeroA << "\n"
              << "Channels in 0.20-0.25:   " << nExpectedBand << "\n"
              << "Median for 0<R<0.60:     " << median << "\n";

    // Per il disegno usiamo un TGraphErrors: i canali non validi
    // non vengono mostrati come punti artificiali a zero.
    TGraphErrors* gR = new TGraphErrors();
    gR->SetName("g_PS_TaggEff");
    gR->SetTitle("");

    int graphPoint = 0;
    double largestVisibleValue = 0.0;

    for (int bin = 1; bin <= nBins; ++bin) {
        const double A = hA->GetBinContent(bin);
        if (A <= 0.0)
            continue;

        const double channel = bin - 1;
        const double value = hR->GetBinContent(bin);
        const double error = hR->GetBinError(bin);

        gR->SetPoint(graphPoint, channel, value);
        gR->SetPointError(graphPoint, 0.0, error);
        ++graphPoint;

        if (std::isfinite(value + error))
            largestVisibleValue =
                std::max(largestVisibleValue, value + error);
    }

    gR->SetMarkerStyle(20);
    gR->SetMarkerSize(0.55);
    gR->SetMarkerColor(kBlue + 1);
    gR->SetLineColor(kBlue + 1);
    gR->SetLineWidth(1);

    TCanvas* c = new TCanvas("c_PS_TaggEff",
                             "PS tagging efficiency",
                             1400, 850);

    c->SetLeftMargin(0.115);
    c->SetRightMargin(0.035);
    c->SetBottomMargin(0.125);
    c->SetTopMargin(0.105);
    c->SetGridx();
    c->SetGridy();

    const double xmin = -1.0;
    const double xmax = nBins;
    const double ymin = -0.01;

    // Mantiene un asse leggibile senza tagliare eventuali punti alti.
    double ymax = 0.42;
    if (largestVisibleValue > ymax)
        ymax = std::min(0.60, 1.08 * largestVisibleValue);

    TH1D* frame = new TH1D("frame_PS_TaggEff", "",
                           nBins, xmin, xmax);

    frame->SetDirectory(nullptr);
    frame->SetMinimum(ymin);
    frame->SetMaximum(ymax);
    frame->GetXaxis()->SetTitle("FPD channel");
    frame->GetYaxis()->SetTitle(
        "PS tagging efficiency, #epsilon_{tag}^{PS}");
    frame->GetXaxis()->CenterTitle(false);
    frame->GetYaxis()->CenterTitle(false);
    frame->GetXaxis()->SetNdivisions(510);
    frame->GetYaxis()->SetNdivisions(508);
    frame->GetXaxis()->SetTitleOffset(1.08);
    frame->GetYaxis()->SetTitleOffset(1.12);
    frame->Draw("AXIS");

    // Fascia attesa 20%-25%, più leggibile delle sole due linee.
    TBox* expectedBand = new TBox(xmin, 0.20, xmax, 0.25);
    expectedBand->SetFillColorAlpha(kAzure - 9, 0.35);
    expectedBand->SetLineColor(kAzure - 4);
    expectedBand->SetLineWidth(1);
    expectedBand->Draw("same");

    TLine* l20 = new TLine(xmin, 0.20, xmax, 0.20);
    TLine* l25 = new TLine(xmin, 0.25, xmax, 0.25);

    l20->SetLineColor(kGray + 2);
    l25->SetLineColor(kGray + 2);
    l20->SetLineStyle(7);
    l25->SetLineStyle(7);
    l20->SetLineWidth(2);
    l25->SetLineWidth(2);
    l20->Draw("same");
    l25->Draw("same");

    // "PZ" mantiene punti e barre d'errore senza linee di collegamento.
    gR->Draw("PZ same");

    // Ridisegna gli assi sopra la fascia.
    frame->Draw("AXIS same");

    TLatex title;
    title.SetNDC();
    title.SetTextFont(42);
    title.SetTextSize(0.047);
    title.SetTextAlign(22);
    title.DrawLatex(
        0.54, 0.955,
        "Tagging efficiency from Pair Spectrometer");

    TLatex info;
    info.SetNDC();
    info.SetTextFont(42);
    info.SetTextSize(0.030);
    info.SetTextColor(kGray + 2);
    info.SetTextAlign(31);
    info.DrawLatex(
        0.955, 0.905,
        TString::Format("R = %.0f (B-C)/A", correctionFactor));

    TLegend* legend = new TLegend(0.145, 0.755, 0.405, 0.875);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.030);
    legend->AddEntry(gR, "PS tagging efficiency", "pe");
    legend->AddEntry(expectedBand, "Expected range: 20-25%", "f");
    legend->Draw();

    c->RedrawAxis();
    c->SaveAs(TString::Format("%s.pdf", outputBase));
    c->SaveAs(TString::Format("%s.png", outputBase));

    TFile* fout = TFile::Open(
        TString::Format("%s.root", outputBase), "RECREATE");

    hA->Write("A_FPD_ScalerAcc");
    hB->Write("B_PairSpec_SumGated");
    hC->Write("C_PairSpec_SumGatedDly");
    hR->Write();
    gR->Write();
    frame->Write();
    c->Write();
    fout->Close();

    txt.close();
    fin->Close();

    std::cout << "\nSaved:\n"
              << " " << outputBase << ".root\n"
              << " " << outputBase << ".pdf\n"
              << " " << outputBase << ".png\n"
              << " " << outputBase << ".txt\n";
}
