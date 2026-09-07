/*
 * diagnose_coherent_lowE.C
 *
 * Purpose
 * -------
 * Diagnostic plotting/summary macro for a ROOT file produced by the low-energy
 * coherent template maker.  It compares no-DeltaPhi and 8/10/12-degree
 * histograms and reports integrals together with m(gamma gamma) / DeltaE plots.
 *
 * Historical default input: coh_lowE_E0_135_155.root.
 * Outputs use the configurable outputPrefix and include *_summary.txt,
 * *_mgg.pdf and *_deltaE.pdf.
 *
 * Usage
 * -----
 *   root -l -b -q 'diagnose_coherent_lowE.C()'
 *
 * Keep LowECoherentCommon.h in the same directory.  The code body below is
 * unchanged from the supplied source.
 */
#include "LowECoherentCommon.h"

#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TString.h"

#include <fstream>
#include <iostream>
#include <algorithm>

namespace {

void StyleHistogram(TH1D* h, int color, int width)
{
    if (!h) return;
    h->SetLineColor(color);
    h->SetLineWidth(width);
    h->SetStats(0);
}

double IntegralAll(TH1D* h)
{
    return h ? h->Integral(0, h->GetNbinsX()+1) : 0.0;
}

} // namespace

void diagnose_coherent_lowE(
    const char* inputFile =
        "coh_lowE_E0_135_155.root",
    const char* outputPrefix =
        "coh_lowE_E0_135_155")
{
    TFile* f = TFile::Open(inputFile, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Cannot open "
                  << inputFile << std::endl;
        return;
    }

    gStyle->SetOptStat(0);

    TString txtName =
        TString::Format("%s_summary.txt", outputPrefix);
    std::ofstream txt(txtName.Data());

    txt << "# bin E_low E_high N_mgg N_before "
           "mean_deltaE rms_deltaE "
           "N_after8 eff8 "
           "N_after10 eff10 "
           "N_after12 eff12\n";

    TCanvas cMgg("cMgg", "mgg", 1200, 900);
    cMgg.Divide(2,2);

    TCanvas cDE("cDE", "DeltaE", 1200, 900);
    cDE.Divide(2,2);

    for (int ie = 0; ie < LowECoh::kNEnergyBins; ++ie) {
        TH1D* hm =
            (TH1D*)f->Get(Form("coh_lowE_mgg_E%d", ie));

        TH1D* h0 =
            (TH1D*)f->Get(
                Form("coh_lowE_deltaE_E%d_nodphi", ie));

        TH1D* h8 =
            (TH1D*)f->Get(
                Form("coh_lowE_deltaE_E%d_dphi8", ie));

        TH1D* h10 =
            (TH1D*)f->Get(
                Form("coh_lowE_deltaE_E%d_dphi10", ie));

        TH1D* h12 =
            (TH1D*)f->Get(
                Form("coh_lowE_deltaE_E%d_dphi12", ie));

        const double nm = IntegralAll(hm);
        const double n0 = IntegralAll(h0);
        const double n8 = IntegralAll(h8);
        const double n10 = IntegralAll(h10);
        const double n12 = IntegralAll(h12);

        txt << ie << " "
            << LowECoh::kEnergyLow[ie] << " "
            << LowECoh::kEnergyHigh[ie] << " "
            << nm << " "
            << n0 << " "
            << (h0 ? h0->GetMean() : 0.0) << " "
            << (h0 ? h0->GetRMS() : 0.0) << " "
            << n8 << " "
            << (n0 > 0.0 ? n8/n0 : 0.0) << " "
            << n10 << " "
            << (n0 > 0.0 ? n10/n0 : 0.0) << " "
            << n12 << " "
            << (n0 > 0.0 ? n12/n0 : 0.0)
            << "\n";

        cMgg.cd(ie+1);
        if (hm && nm > 0.0) {
            StyleHistogram(hm, kBlack, 2);
            hm->GetXaxis()->SetRangeUser(80.0, 190.0);
            hm->Draw("hist");

            TLatex latex;
            latex.SetNDC();
            latex.SetTextSize(0.045);
            latex.DrawLatex(
                0.16, 0.84,
                Form("%.0f-%.0f MeV",
                     LowECoh::kEnergyLow[ie],
                     LowECoh::kEnergyHigh[ie]));
            latex.DrawLatex(
                0.16, 0.77,
                Form("entries = %.0f", nm));
        } else {
            TLatex latex;
            latex.SetNDC();
            latex.SetTextSize(0.06);
            latex.DrawLatex(0.25, 0.50,
                            "No events in this bin");
        }

        cDE.cd(ie+1);
        if (h0 && n0 > 0.0) {
            StyleHistogram(h0, kBlack, 3);
            StyleHistogram(h8, kRed+1, 2);
            StyleHistogram(h10, kBlue+1, 2);
            StyleHistogram(h12, kGreen+2, 2);

            double ymax = h0->GetMaximum();
            if (h8) ymax = std::max(ymax, h8->GetMaximum());
            if (h10) ymax = std::max(ymax, h10->GetMaximum());
            if (h12) ymax = std::max(ymax, h12->GetMaximum());

            h0->SetMaximum(1.20*ymax);
            h0->GetXaxis()->SetRangeUser(-50.0, 40.0);
            h0->Draw("hist");
            if (h8) h8->Draw("hist same");
            if (h10) h10->Draw("hist same");
            if (h12) h12->Draw("hist same");

            TLegend leg(0.57, 0.62, 0.88, 0.87);
            leg.SetBorderSize(0);
            leg.SetFillStyle(0);
            leg.AddEntry(h0, "No #Delta#Phi cut", "l");
            leg.AddEntry(h8, "#Delta#Phi < 8^{#circ}", "l");
            leg.AddEntry(h10, "#Delta#Phi < 10^{#circ}", "l");
            leg.AddEntry(h12, "#Delta#Phi < 12^{#circ}", "l");
            leg.Draw();

            TLatex latex;
            latex.SetNDC();
            latex.SetTextSize(0.039);
            latex.DrawLatex(
                0.15, 0.86,
                Form("%.0f-%.0f MeV",
                     LowECoh::kEnergyLow[ie],
                     LowECoh::kEnergyHigh[ie]));
            latex.DrawLatex(
                0.15, 0.80,
                Form("#mu = %.2f MeV, RMS = %.2f MeV",
                     h0->GetMean(), h0->GetRMS()));
            latex.DrawLatex(
                0.15, 0.74,
                Form("#epsilon_{8/10/12} = "
                     "%.1f%% / %.1f%% / %.1f%%",
                     n0 > 0.0 ? 100.0*n8/n0 : 0.0,
                     n0 > 0.0 ? 100.0*n10/n0 : 0.0,
                     n0 > 0.0 ? 100.0*n12/n0 : 0.0));
        } else {
            TLatex latex;
            latex.SetNDC();
            latex.SetTextSize(0.06);
            latex.DrawLatex(0.25, 0.50,
                            "No events in this bin");
        }
    }

    TString mggPdf =
        TString::Format("%s_mgg.pdf", outputPrefix);
    TString dePdf =
        TString::Format("%s_deltaE.pdf", outputPrefix);

    cMgg.SaveAs(mggPdf.Data());
    cDE.SaveAs(dePdf.Data());

    txt.close();
    f->Close();

    std::cout << "Wrote " << mggPdf << std::endl;
    std::cout << "Wrote " << dePdf << std::endl;
    std::cout << "Wrote " << txtName << std::endl;
}
