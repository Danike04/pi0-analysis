/*
 * compare_data_incoherent_mm.C
 *
 * PURPOSE
 *   Compare data and an incoherent-MC missing-mass template globally and in
 *   paper energy bins. The incoherent template is normalized in a configurable
 *   high-sideband (default 45-120 MeV), with a peak-height fallback if the
 *   sideband cannot be used.
 *
 * EXPECTED INPUT HISTOGRAMS
 *   Data ROOT: h_mm_after_mgg_cut, h_mm_vs_egamma_after_mgg_cut.
 *   Incoherent ROOT: inc_mm_after_mgg, inc_mm_vs_egamma_after_mgg.
 *
 * HISTORICAL DEPENDENCY
 *   The default file is `incoherent_mc_diagnostics.root`. No producer of that
 *   exact diagnostic file/histogram pair was found in the supplied archive,
 *   so this macro is retained as a diagnostic with that dependency explicitly
 *   documented rather than presented as a self-contained workflow step.
 *
 * USAGE EXAMPLE
 *   root -l -b -q 'compare_data_incoherent_mm.C()'
 *
 * REPOSITORY STATUS
 *   Data/incoherent-MC diagnostic; executable code unchanged.
 */
#include <iostream>
#include <fstream>
#include <cmath>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TString.h"
#include "TStyle.h"
#include "TMath.h"

static double integral_range(TH1* h, double xmin, double xmax)
{
    if(!h) return 0.0;
    int b1 = h->GetXaxis()->FindBin(xmin + 1.0e-6);
    int b2 = h->GetXaxis()->FindBin(xmax - 1.0e-6);
    return h->Integral(b1, b2);
}

static double maximum_range(TH1* h, double xmin, double xmax)
{
    if(!h) return 0.0;
    int b1 = h->GetXaxis()->FindBin(xmin + 1.0e-6);
    int b2 = h->GetXaxis()->FindBin(xmax - 1.0e-6);
    double m = 0.0;
    for(int b=b1; b<=b2; b++) {
        if(h->GetBinContent(b) > m) m = h->GetBinContent(b);
    }
    return m;
}

static double choose_scale(TH1* hData, TH1* hInc,
                           double normMin, double normMax,
                           double drawMin, double drawMax,
                           bool& usedSideband)
{
    usedSideband = false;
    double dSide = integral_range(hData, normMin, normMax);
    double iSide = integral_range(hInc, normMin, normMax);

    if(dSide > 0.0 && iSide > 0.0) {
        usedSideband = true;
        return dSide / iSide;
    }

    double dMax = maximum_range(hData, drawMin, drawMax);
    double iMax = maximum_range(hInc, drawMin, drawMax);
    if(dMax > 0.0 && iMax > 0.0) return dMax / iMax;
    return 1.0;
}

static void style_data(TH1* h)
{
    h->SetLineColor(kBlack);
    h->SetMarkerColor(kBlack);
    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.65);
    h->SetLineWidth(1);
}

static void style_inc(TH1* h)
{
    h->SetLineColor(kRed+1);
    h->SetMarkerColor(kRed+1);
    h->SetLineWidth(2);
    h->SetFillStyle(0);
}

void compare_data_incoherent_mm(
    const char* dataFile="paper_style_dsigma_dOmega_cm_from_acqu.root",
    const char* incFile="incoherent_mc_diagnostics.root",
    double normMin=45.0,
    double normMax=120.0,
    double signalMin=-10.0,
    double signalMax=40.0)
{
    gStyle->SetOptStat(0);

    TFile* fd = TFile::Open(dataFile);
    TFile* fi = TFile::Open(incFile);
    if(!fd || fd->IsZombie()) {
        std::cout << "Cannot open data file: " << dataFile << std::endl;
        return;
    }
    if(!fi || fi->IsZombie()) {
        std::cout << "Cannot open incoherent MC file: " << incFile << std::endl;
        return;
    }

    TH1D* hDataGlobal = (TH1D*)fd->Get("h_mm_after_mgg_cut");
    TH1D* hIncGlobal  = (TH1D*)fi->Get("inc_mm_after_mgg");
    TH2D* hData2D = (TH2D*)fd->Get("h_mm_vs_egamma_after_mgg_cut");
    TH2D* hInc2D  = (TH2D*)fi->Get("inc_mm_vs_egamma_after_mgg");

    if(!hDataGlobal || !hIncGlobal || !hData2D || !hInc2D) {
        std::cout << "Missing histogram. Expected:" << std::endl;
        std::cout << "  data: h_mm_after_mgg_cut, h_mm_vs_egamma_after_mgg_cut" << std::endl;
        std::cout << "  MC:   inc_mm_after_mgg, inc_mm_vs_egamma_after_mgg" << std::endl;
        fd->ls();
        fi->ls();
        return;
    }

    const double drawMin = -100.0;
    const double drawMax = 150.0;

    TFile* fout = new TFile("data_vs_incoherent_mm.root", "RECREATE");
    std::ofstream txt("data_vs_incoherent_mm.txt");
    txt << "# Preliminary data vs incoherent ppnn MC comparison\n";
    txt << "# MC normalization sideband = [" << normMin << "," << normMax << "] MeV\n";
    txt << "# signal window = [" << signalMin << "," << signalMax << "] MeV\n";
    txt << "# IMPORTANT: sideband normalization is preliminary and assumes the selected sideband is described by this MC component.\n";
    txt << "# E_low E_high data_side inc_side scale used_sideband data_signal predicted_inc_signal inc_over_data_signal\n";

    // Global comparison.
    TH1D* dGlobal = (TH1D*)hDataGlobal->Clone("cmp_data_mm_global");
    TH1D* iGlobal = (TH1D*)hIncGlobal->Clone("cmp_inc_mm_global_scaled");
    dGlobal->SetDirectory(0);
    iGlobal->SetDirectory(0);

    bool globalSide = false;
    double globalScale = choose_scale(dGlobal, iGlobal, normMin, normMax,
                                      drawMin, drawMax, globalSide);
    iGlobal->Scale(globalScale);
    style_data(dGlobal);
    style_inc(iGlobal);

    TCanvas* cGlobal = new TCanvas("c_data_vs_inc_global", "Data vs incoherent MC", 1100, 750);
    dGlobal->SetTitle("Missing mass after m_{#gamma#gamma} cut;MM-M_{^{4}He} (MeV);Prompt-random counts / scaled MC");
    dGlobal->GetXaxis()->SetRangeUser(drawMin, drawMax);
    double ymax = TMath::Max(dGlobal->GetMaximum(), iGlobal->GetMaximum());
    dGlobal->SetMaximum(1.20*ymax);
    dGlobal->Draw("E");
    iGlobal->Draw("HIST SAME");

    TLine* l1 = new TLine(signalMin, 0.0, signalMin, 1.15*ymax);
    TLine* l2 = new TLine(signalMax, 0.0, signalMax, 1.15*ymax);
    l1->SetLineStyle(2); l2->SetLineStyle(2);
    l1->SetLineWidth(2); l2->SetLineWidth(2);
    l1->Draw(); l2->Draw();

    TLegend* leg = new TLegend(0.55,0.70,0.88,0.88);
    leg->SetBorderSize(0);
    leg->AddEntry(dGlobal, "Data: prompt - random", "lep");
    leg->AddEntry(iGlobal, "Incoherent ppnn MC (scaled)", "l");
    leg->AddEntry(l1, "Current MM selection", "l");
    leg->Draw();

    TLatex label;
    label.SetNDC();
    label.SetTextSize(0.032);
    if(globalSide)
        label.DrawLatex(0.14,0.86,Form("MC normalized in %.0f < MM < %.0f MeV",normMin,normMax));
    else
        label.DrawLatex(0.14,0.86,"Sideband invalid: MC normalized to the data maximum");

    cGlobal->SaveAs("data_vs_incoherent_mm_global.pdf");

    double dSideG = integral_range(dGlobal, normMin, normMax);
    double iSideRawG = integral_range(hIncGlobal, normMin, normMax);
    double dSigG = integral_range(dGlobal, signalMin, signalMax);
    double iSigG = integral_range(iGlobal, signalMin, signalMax);
    double fracG = (dSigG != 0.0) ? iSigG/dSigG : 0.0;

    std::cout << "Global sideband scale = " << globalScale
              << " (used sideband: " << (globalSide ? "yes" : "no") << ")" << std::endl;
    std::cout << "Global predicted incoherent counts in MM window = " << iSigG << std::endl;
    std::cout << "Global data counts in MM window = " << dSigG << std::endl;
    std::cout << "Preliminary incoherent/data ratio in MM window = " << fracG << std::endl;

    txt << "GLOBAL GLOBAL " << dSideG << " " << iSideRawG << " " << globalScale
        << " " << (globalSide ? 1 : 0) << " " << dSigG << " " << iSigG << " " << fracG << "\n";

    // Paper energy bins.
    const int nE = 17;
    double eLow[nE]  = {201,211,223,235,247,259,271,283,295,309,319,331,343,356,367,379,391};
    double eHigh[nE] = {210,222,234,246,258,270,282,294,308,318,330,342,355,366,378,390,401};

    TString pdfName = "data_vs_incoherent_mm_by_energy.pdf";
    TCanvas* cE = new TCanvas("c_data_vs_inc_energy", "Data vs incoherent MC by energy", 1300, 850);
    cE->Divide(3,2);
    cE->Print(pdfName + "[");

    for(int ie=0; ie<nE; ie++) {
        int bx1d = hData2D->GetXaxis()->FindBin(eLow[ie] + 1.0e-6);
        int bx2d = hData2D->GetXaxis()->FindBin(eHigh[ie] - 1.0e-6);
        int bx1i = hInc2D->GetXaxis()->FindBin(eLow[ie] + 1.0e-6);
        int bx2i = hInc2D->GetXaxis()->FindBin(eHigh[ie] - 1.0e-6);

        TH1D* hd = hData2D->ProjectionY(Form("cmp_data_mm_Ebin%d",ie), bx1d, bx2d, "e");
        TH1D* hi = hInc2D->ProjectionY(Form("cmp_inc_mm_Ebin%d_scaled",ie), bx1i, bx2i, "e");
        hd->SetDirectory(0);
        hi->SetDirectory(0);

        bool usedSide = false;
        double scale = choose_scale(hd, hi, normMin, normMax, drawMin, drawMax, usedSide);
        double incSideRaw = integral_range(hi, normMin, normMax);
        hi->Scale(scale);

        double dataSide = integral_range(hd, normMin, normMax);
        double dataSig = integral_range(hd, signalMin, signalMax);
        double incSig = integral_range(hi, signalMin, signalMax);
        double ratio = (dataSig != 0.0) ? incSig/dataSig : 0.0;

        txt << eLow[ie] << " " << eHigh[ie] << " "
            << dataSide << " " << incSideRaw << " " << scale << " "
            << (usedSide ? 1 : 0) << " " << dataSig << " " << incSig << " " << ratio << "\n";

        style_data(hd);
        style_inc(hi);
        hd->SetTitle(Form("%.0f < E_{#gamma} < %.0f MeV;MM-M_{^{4}He} (MeV);Counts",eLow[ie],eHigh[ie]));
        hd->GetXaxis()->SetRangeUser(drawMin, drawMax);
        double ym = TMath::Max(hd->GetMaximum(),hi->GetMaximum());
        if(ym <= 0.0) ym = 1.0;
        hd->SetMaximum(1.28*ym);

        int pad = (ie % 6) + 1;
        cE->cd(pad);
        hd->Draw("E");
        hi->Draw("HIST SAME");

        TLine* a = new TLine(signalMin,0.0,signalMin,1.18*ym);
        TLine* b = new TLine(signalMax,0.0,signalMax,1.18*ym);
        a->SetLineStyle(2); b->SetLineStyle(2);
        a->Draw(); b->Draw();

        TLatex t;
        t.SetNDC();
        t.SetTextSize(0.045);
        t.DrawLatex(0.13,0.86,Form("inc/data in cut = %.3g",ratio));
        t.SetTextSize(0.035);
        t.DrawLatex(0.13,0.80, usedSide ? "sideband scaled" : "maximum scaled");

        fout->cd();
        hd->Write();
        hi->Write();

        if(pad == 6 || ie == nE-1) {
            cE->Print(pdfName);
            if(ie != nE-1) {
                cE->Clear();
                cE->Divide(3,2);
            }
        }
    }
    cE->Print(pdfName + "]");

    fout->cd();
    dGlobal->Write();
    iGlobal->Write();
    cGlobal->Write();
    fout->Close();
    txt.close();

    std::cout << "Saved:" << std::endl;
    std::cout << "  data_vs_incoherent_mm_global.pdf" << std::endl;
    std::cout << "  data_vs_incoherent_mm_by_energy.pdf" << std::endl;
    std::cout << "  data_vs_incoherent_mm.root" << std::endl;
    std::cout << "  data_vs_incoherent_mm.txt" << std::endl;
    std::cout << "Interpret the sideband-scaled numbers as preliminary only." << std::endl;

    fd->Close();
    fi->Close();
}
