/*
 * fig2_breakup_shape_comparison_366.C
 *
 * Purpose
 * -------
 * Shape-only comparison of d+d, t+p and 3He+n breakup Monte Carlo at the
 * 366-MeV Fig. 2 energy point.  Each channel is mapped to a common binning and
 * normalized to equal area in the requested display range; means/RMS values
 * are also written.  No absolute breakup normalization is inferred here.
 *
 * Default inputs
 * --------------
 *   fig2_dd_E3_366_wide.root
 *   fig2_tp_E3_366_wide.root
 *   fig2_he3n_all_available_bins_wide.root
 *
 * Outputs
 * -------
 *   fig2_breakup_shape_comparison_366.root
 *   fig2_breakup_shape_comparison_366.txt
 *   fig2_breakup_shape_comparison_366.pdf
 *
 * Usage
 * -----
 *   root -l -b -q 'fig2_breakup_shape_comparison_366.C()'
 *
 * The analysis body below is unchanged from the supplied source.
 */
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"

namespace Fig2ShapeComparison366 {

TH1D* MapToCommonBinning(TH1* source, const char* name,
                         int nBins, double histMin, double histMax)
{
    if (!source) return 0;

    TH1D* out = new TH1D(name, "", nBins, histMin, histMax);
    out->SetDirectory(0);
    out->Sumw2();

    for (int b = 1; b <= source->GetNbinsX(); ++b) {
        const double x = source->GetXaxis()->GetBinCenter(b);
        if (x < histMin || x >= histMax) continue;

        const int bo = out->GetXaxis()->FindBin(x);
        const double c0 = out->GetBinContent(bo);
        const double e0 = out->GetBinError(bo);
        const double c1 = source->GetBinContent(b);
        const double e1 = source->GetBinError(b);

        out->SetBinContent(bo, c0 + c1);
        out->SetBinError(bo, std::sqrt(e0*e0 + e1*e1));
    }
    return out;
}

double IntegralInRange(TH1D* h, double xmin, double xmax)
{
    if (!h) return 0.0;
    const int b1 = h->GetXaxis()->FindBin(xmin + 1.0e-9);
    const int b2 = h->GetXaxis()->FindBin(xmax - 1.0e-9);
    return h->Integral(b1, b2);
}

bool NormalizeToUnitArea(TH1D* h, double xmin, double xmax)
{
    const double area = IntegralInRange(h, xmin, xmax);
    if (area <= 0.0) return false;
    h->Scale(1.0 / area);
    return true;
}

void MomentsInRange(TH1D* h, double xmin, double xmax,
                    double& mean, double& rms)
{
    mean = 0.0;
    rms = 0.0;
    if (!h) return;

    const int b1 = h->GetXaxis()->FindBin(xmin + 1.0e-9);
    const int b2 = h->GetXaxis()->FindBin(xmax - 1.0e-9);

    double sumW = 0.0;
    double sumWX = 0.0;
    double sumWX2 = 0.0;

    for (int b = b1; b <= b2; ++b) {
        const double w = h->GetBinContent(b);
        const double x = h->GetXaxis()->GetBinCenter(b);
        sumW += w;
        sumWX += w*x;
        sumWX2 += w*x*x;
    }

    if (sumW <= 0.0) return;

    mean = sumWX / sumW;
    const double variance = sumWX2/sumW - mean*mean;
    rms = variance > 0.0 ? std::sqrt(variance) : 0.0;
}

void StyleDD(TH1D* h)
{
    h->SetLineColor(kBlue + 1);
    h->SetLineWidth(3);
    h->SetLineStyle(1);
    h->SetFillStyle(0);
}

void StyleTP(TH1D* h)
{
    h->SetLineColor(kMagenta + 2);
    h->SetLineWidth(3);
    h->SetLineStyle(7);
    h->SetFillStyle(0);
}

void StyleHe3n(TH1D* h)
{
    h->SetLineColor(kGreen + 2);
    h->SetLineWidth(3);
    h->SetLineStyle(2);
    h->SetFillStyle(0);
}

} // namespace Fig2ShapeComparison366


void fig2_breakup_shape_comparison_366(
    const char* ddFile = "fig2_dd_E3_366_wide.root",
    const char* tpFile = "fig2_tp_E3_366_wide.root",
    const char* he3nFile = "fig2_he3n_all_available_bins_wide.root",
    double displayMin = -60.0,
    double displayMax = 40.0,
    double histMin = -120.0,
    double histMax = 120.0,
    int nBins = 240)
{
    using namespace Fig2ShapeComparison366;

    gStyle->SetOptStat(0);

    TFile* fDD = TFile::Open(ddFile, "READ");
    TFile* fTP = TFile::Open(tpFile, "READ");
    TFile* fHe3n = TFile::Open(he3nFile, "READ");

    if (!fDD || fDD->IsZombie()) {
        std::cerr << "Cannot open d+d file: " << ddFile << std::endl;
        return;
    }
    if (!fTP || fTP->IsZombie()) {
        std::cerr << "Cannot open t+p file: " << tpFile << std::endl;
        fDD->Close();
        return;
    }
    if (!fHe3n || fHe3n->IsZombie()) {
        std::cerr << "Cannot open 3He+n file: " << he3nFile << std::endl;
        fDD->Close();
        fTP->Close();
        return;
    }

    const int nPanels = 4;
    const int cuts[nPanels] = {-1, 8, 10, 12};
    const char* panelLabels[nPanels] = {
        "before #Delta#Phi cut",
        "#Delta#Phi < 8^{#circ}",
        "#Delta#Phi < 10^{#circ}",
        "#Delta#Phi < 12^{#circ}"
    };

    TFile* output = TFile::Open(
        "fig2_breakup_shape_comparison_366.root", "RECREATE");

    if (!output || output->IsZombie()) {
        std::cerr << "Cannot create output ROOT file." << std::endl;
        fDD->Close();
        fTP->Close();
        fHe3n->Close();
        return;
    }

    std::ofstream txt("fig2_breakup_shape_comparison_366.txt");
    txt << "# Each template is normalized independently to unit area in ["
        << displayMin << "," << displayMax << "] MeV.\n"
        << "# panel channel raw_visible_area mean rms\n";

    TCanvas* c = new TCanvas(
        "c_breakup_shape_comparison_366",
        "366 MeV breakup-shape comparison",
        1300, 900);
    c->Divide(2, 2, 0.002, 0.002);

    for (int ip = 0; ip < nPanels; ++ip) {
        const TString hname = cuts[ip] < 0
            ? "inc_deltaE_E3_nodphi"
            : Form("inc_deltaE_E3_dphi%d", cuts[ip]);

        TH1* ddRaw = dynamic_cast<TH1*>(fDD->Get(hname));
        TH1* tpRaw = dynamic_cast<TH1*>(fTP->Get(hname));
        TH1* he3nRaw = dynamic_cast<TH1*>(fHe3n->Get(hname));

        if (!ddRaw || !tpRaw || !he3nRaw) {
            std::cerr << "Missing histogram " << hname
                      << " in one or more files." << std::endl;
            continue;
        }

        TH1D* dd = MapToCommonBinning(
            ddRaw, Form("dd_shape_panel%d", ip),
            nBins, histMin, histMax);
        TH1D* tp = MapToCommonBinning(
            tpRaw, Form("tp_shape_panel%d", ip),
            nBins, histMin, histMax);
        TH1D* he3n = MapToCommonBinning(
            he3nRaw, Form("he3n_shape_panel%d", ip),
            nBins, histMin, histMax);

        if (!dd || !tp || !he3n) continue;

        const double ddRawArea = IntegralInRange(dd, displayMin, displayMax);
        const double tpRawArea = IntegralInRange(tp, displayMin, displayMax);
        const double he3nRawArea = IntegralInRange(he3n, displayMin, displayMax);

        if (!NormalizeToUnitArea(dd, displayMin, displayMax) ||
            !NormalizeToUnitArea(tp, displayMin, displayMax) ||
            !NormalizeToUnitArea(he3n, displayMin, displayMax)) {
            std::cerr << "Zero visible integral for " << hname << std::endl;
            continue;
        }

        double ddMean = 0.0, ddRms = 0.0;
        double tpMean = 0.0, tpRms = 0.0;
        double he3nMean = 0.0, he3nRms = 0.0;

        MomentsInRange(dd, displayMin, displayMax, ddMean, ddRms);
        MomentsInRange(tp, displayMin, displayMax, tpMean, tpRms);
        MomentsInRange(he3n, displayMin, displayMax, he3nMean, he3nRms);

        txt << ip << " dd " << ddRawArea << " "
            << ddMean << " " << ddRms << "\n";
        txt << ip << " tp " << tpRawArea << " "
            << tpMean << " " << tpRms << "\n";
        txt << ip << " he3n " << he3nRawArea << " "
            << he3nMean << " " << he3nRms << "\n";

        StyleDD(dd);
        StyleTP(tp);
        StyleHe3n(he3n);

        c->cd(ip + 1);
        gPad->SetTicks(1, 1);
        gPad->SetLeftMargin(0.13);
        gPad->SetBottomMargin(0.13);

        const double ymax = std::max(
            dd->GetMaximum(),
            std::max(tp->GetMaximum(), he3n->GetMaximum()));

        dd->SetTitle(Form("366 MeV, %s", panelLabels[ip]));
        dd->GetXaxis()->SetTitle("#Delta E_{#pi^{0}}^{*} (MeV)");
        dd->GetYaxis()->SetTitle("Arbitrary units (equal area)");
        dd->GetXaxis()->SetRangeUser(displayMin, displayMax);
        dd->SetMinimum(0.0);
        dd->SetMaximum(ymax > 0.0 ? 1.28*ymax : 1.0);

        dd->Draw("HIST");
        tp->Draw("HIST SAME");
        he3n->Draw("HIST SAME");

        TLegend* leg = new TLegend(0.56, 0.68, 0.89, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.035);
        leg->AddEntry(dd, "d+d", "l");
        leg->AddEntry(tp, "t+p", "l");
        leg->AddEntry(he3n, "^{3}He+n", "l");
        leg->Draw();

        TLatex label;
        label.SetNDC();
        label.SetTextSize(0.029);
        label.DrawLatex(
            0.15, 0.88,
            Form("unit area in [%.0f, %.0f] MeV",
                 displayMin, displayMax));

        output->cd();
        dd->Write();
        tp->Write();
        he3n->Write();
    }

    c->SaveAs("fig2_breakup_shape_comparison_366.pdf");
    output->cd();
    c->Write();

    txt.close();
    output->Close();
    fDD->Close();
    fTP->Close();
    fHe3n->Close();

    std::cout << "Saved:\n"
              << "  fig2_breakup_shape_comparison_366.pdf\n"
              << "  fig2_breakup_shape_comparison_366.root\n"
              << "  fig2_breakup_shape_comparison_366.txt\n";
}
