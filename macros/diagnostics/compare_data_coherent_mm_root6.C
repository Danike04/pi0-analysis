/*
 * compare_data_coherent_mm_root6.C
 *
 * PURPOSE
 *   Compare the reconstructed missing-mass distribution in data with coherent
 *   MC globally and in the 17 paper photon-energy bins.
 *
 * EXPECTED INPUT HISTOGRAMS
 *   Data ROOT: h_mm_after_mgg_cut, h_mm_vs_egamma_after_mgg_cut.
 *   Coherent MC ROOT: coh_mm_after_mgg, coh_mm_vs_egamma_after_mgg.
 *   The latter are produced by analyze_coherent_mc_root6.C.
 *
 * NORMALIZATION / REGIONS
 *   Uses configurable core, signal and high-side tail ranges to scale/compare
 *   shapes and reports integrals in text output.
 *
 * DEFAULT INPUTS
 *   paper_style_dsigma_dOmega_cm_from_acqu.root
 *   coherent_mc_diagnostics.root
 *
 * USAGE EXAMPLE
 *   root -l -b -q 'compare_data_coherent_mm_root6.C()'
 *
 * REPOSITORY STATUS
 *   Useful data/MC validation diagnostic; executable code unchanged.
 */
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TString.h"
#include "TStyle.h"

namespace {

TH1D* copy_to_common_axis(const TH1* src,
                          const char* name,
                          const char* title,
                          int nBins=250,
                          double xMin=-100.0,
                          double xMax=150.0)
{
    if(!src) return nullptr;

    TH1D* out = new TH1D(name, title, nBins, xMin, xMax);
    out->Sumw2();

    for(int b=1; b<=out->GetNbinsX(); ++b) {
        const double x = out->GetBinCenter(b);
        const int bs = src->GetXaxis()->FindBin(x);
        if(bs < 1 || bs > src->GetNbinsX()) continue;
        out->SetBinContent(b, src->GetBinContent(bs));
        out->SetBinError(b, src->GetBinError(bs));
    }

    return out;
}

double integral_range(const TH1* h, double x1, double x2)
{
    if(!h) return 0.0;
    const int b1 = h->GetXaxis()->FindBin(x1 + 1.0e-9);
    const int b2 = h->GetXaxis()->FindBin(x2 - 1.0e-9);
    return h->Integral(b1, b2);
}

void draw_vertical(double x, double ymax, int style, int width=2)
{
    TLine* l = new TLine(x, 0.0, x, ymax);
    l->SetLineStyle(style);
    l->SetLineWidth(width);
    l->Draw("same");
}

} // namespace

void compare_data_coherent_mm_root6(
    const char* dataFile="paper_style_dsigma_dOmega_cm_from_acqu.root",
    const char* coherentFile="coherent_mc_diagnostics.root",
    double coreMin=-10.0,
    double coreMax=15.0,
    double signalMin=-10.0,
    double signalMax=40.0,
    double tailMin=40.0,
    double tailMax=120.0)
{
    gStyle->SetOptStat(0);

    TFile* fd = TFile::Open(dataFile, "READ");
    if(!fd || fd->IsZombie()) {
        std::cerr << "Cannot open data file: " << dataFile << std::endl;
        return;
    }

    TFile* fc = TFile::Open(coherentFile, "READ");
    if(!fc || fc->IsZombie()) {
        std::cerr << "Cannot open coherent MC file: " << coherentFile << std::endl;
        fd->Close();
        return;
    }

    TH1* hDataRaw = dynamic_cast<TH1*>(fd->Get("h_mm_after_mgg_cut"));
    TH1* hCohRaw  = dynamic_cast<TH1*>(fc->Get("coh_mm_after_mgg"));
    TH2* hDataERaw = dynamic_cast<TH2*>(fd->Get("h_mm_vs_egamma_after_mgg_cut"));
    TH2* hCohERaw  = dynamic_cast<TH2*>(fc->Get("coh_mm_vs_egamma_after_mgg"));

    if(!hDataRaw || !hCohRaw || !hDataERaw || !hCohERaw) {
        std::cerr << "Missing required histograms. Expected:" << std::endl;
        std::cerr << "  data: h_mm_after_mgg_cut, h_mm_vs_egamma_after_mgg_cut" << std::endl;
        std::cerr << "  MC:   coh_mm_after_mgg, coh_mm_vs_egamma_after_mgg" << std::endl;
        fc->Close();
        fd->Close();
        return;
    }

    const double eLow[17]  = {201,211,223,235,247,259,271,283,295,309,319,331,343,356,367,379,391};
    const double eHigh[17] = {210,222,234,246,258,270,282,294,308,318,330,342,355,366,378,390,401};

    TH1D* hData = copy_to_common_axis(hDataRaw, "data_mm_common",
        "Data vs coherent MC after m_{#gamma#gamma} cut;MM-M_{^{4}He} (MeV);Counts / MeV");
    TH1D* hCoh = copy_to_common_axis(hCohRaw, "coherent_mm_common",
        "Coherent MC;MM-M_{^{4}He} (MeV);Counts / MeV");

    const double dataCore = integral_range(hData, coreMin, coreMax);
    const double cohCore = integral_range(hCoh, coreMin, coreMax);
    const double scale = (cohCore > 0.0) ? dataCore / cohCore : 0.0;
    hCoh->Scale(scale);

    const double dataSignal = integral_range(hData, signalMin, signalMax);
    const double cohSignal = integral_range(hCoh, signalMin, signalMax);
    const double dataTail = integral_range(hData, tailMin, tailMax);
    const double cohTail = integral_range(hCoh, tailMin, tailMax);
    const double tailExcess = dataTail - cohTail;

    std::ofstream txt("data_vs_coherent_mm.txt");
    txt << "# Data vs coherent MC shape comparison\n";
    txt << "# coherent MC normalized to data in core = [" << coreMin << "," << coreMax << "] MeV\n";
    txt << "# signal window = [" << signalMin << "," << signalMax << "] MeV\n";
    txt << "# positive tail = [" << tailMin << "," << tailMax << "] MeV\n";
    txt << "# This is a shape comparison, not an absolute MC normalization.\n";
    txt << "# E_low E_high data_core mc_core_raw scale data_signal mc_signal_scaled data_tail mc_tail_scaled tail_excess\n";
    txt << "GLOBAL GLOBAL " << dataCore << " " << cohCore << " " << scale << " "
        << dataSignal << " " << cohSignal << " " << dataTail << " " << cohTail << " " << tailExcess << "\n";

    TFile* fout = new TFile("data_vs_coherent_mm.root", "RECREATE");

    TCanvas* cGlobal = new TCanvas("c_data_vs_coherent_global", "Data vs coherent MC", 1000, 750);
    cGlobal->SetGrid();

    hData->SetMarkerStyle(20);
    hData->SetMarkerSize(0.7);
    hData->SetLineColor(kBlack);
    hData->SetMarkerColor(kBlack);

    hCoh->SetLineColor(kBlue+1);
    hCoh->SetLineWidth(3);

    const double ymaxGlobal = 1.18 * std::max(hData->GetMaximum(), hCoh->GetMaximum());
    hData->SetMaximum(ymaxGlobal);
    hData->GetXaxis()->SetRangeUser(-60.0, 120.0);
    hData->Draw("E");
    hCoh->Draw("HIST SAME");

    draw_vertical(signalMin, ymaxGlobal, 2);
    draw_vertical(signalMax, ymaxGlobal, 2);
    draw_vertical(coreMax, ymaxGlobal, 3, 1);

    TLegend* leg = new TLegend(0.54,0.68,0.89,0.89);
    leg->SetBorderSize(0);
    leg->AddEntry(hData, "Data: prompt - random", "lep");
    leg->AddEntry(hCoh, "Coherent MC, core-scaled", "l");
    leg->AddEntry((TObject*)0, Form("core: %.0f < MM < %.0f MeV",coreMin,coreMax), "");
    leg->Draw();

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.033);
    lat.DrawLatex(0.13,0.87,Form("Signal: data %.0f, coherent MC %.0f",dataSignal,cohSignal));
    lat.DrawLatex(0.13,0.82,Form("Tail %.0f--%.0f MeV: data %.0f, MC %.0f, excess %.0f",tailMin,tailMax,dataTail,cohTail,tailExcess));

    cGlobal->Print("data_vs_coherent_mm_global.pdf");
    hData->Write();
    hCoh->Write();
    cGlobal->Write();

    TCanvas* cBook = new TCanvas("c_data_vs_coherent_by_energy", "Data vs coherent MC by energy", 1400, 850);
    cBook->Print("data_vs_coherent_mm_by_energy.pdf[");

    for(int page=0; page<3; ++page) {
        cBook->Clear();
        cBook->Divide(3,2,0.001,0.001);

        for(int ip=0; ip<6; ++ip) {
            const int ie = 6*page + ip;
            if(ie >= 17) continue;

            cBook->cd(ip+1);
            gPad->SetGrid();

            const int xd1 = hDataERaw->GetXaxis()->FindBin(eLow[ie] + 1.0e-6);
            const int xd2 = hDataERaw->GetXaxis()->FindBin(eHigh[ie] - 1.0e-6);
            const int xc1 = hCohERaw->GetXaxis()->FindBin(eLow[ie] + 1.0e-6);
            const int xc2 = hCohERaw->GetXaxis()->FindBin(eHigh[ie] - 1.0e-6);

            TH1D* dRaw = hDataERaw->ProjectionY(Form("data_mm_E%d_raw",ie), xd1, xd2, "e");
            TH1D* cRaw = hCohERaw->ProjectionY(Form("coh_mm_E%d_raw",ie), xc1, xc2, "e");
            dRaw->SetDirectory(nullptr);
            cRaw->SetDirectory(nullptr);

            TH1D* d = copy_to_common_axis(dRaw, Form("data_mm_E%d",ie),
                Form("%.0f < E_{#gamma} < %.0f MeV;MM-M_{^{4}He} (MeV);Counts / MeV",eLow[ie],eHigh[ie]));
            TH1D* c = copy_to_common_axis(cRaw, Form("coh_mm_E%d",ie), "coherent MC");
            delete dRaw;
            delete cRaw;

            const double dCore = integral_range(d, coreMin, coreMax);
            const double cCoreRaw = integral_range(c, coreMin, coreMax);
            const double s = (cCoreRaw > 0.0) ? dCore / cCoreRaw : 0.0;
            c->Scale(s);

            const double dSig = integral_range(d, signalMin, signalMax);
            const double cSig = integral_range(c, signalMin, signalMax);
            const double dTail = integral_range(d, tailMin, tailMax);
            const double cTail = integral_range(c, tailMin, tailMax);
            const double excess = dTail - cTail;

            txt << eLow[ie] << " " << eHigh[ie] << " " << dCore << " " << cCoreRaw << " " << s << " "
                << dSig << " " << cSig << " " << dTail << " " << cTail << " " << excess << "\n";

            d->SetMarkerStyle(20);
            d->SetMarkerSize(0.45);
            d->SetLineColor(kBlack);
            d->SetMarkerColor(kBlack);
            c->SetLineColor(kBlue+1);
            c->SetLineWidth(2);

            const double ymax = 1.18 * std::max(d->GetMaximum(), c->GetMaximum());
            d->SetMaximum(ymax > 0.0 ? ymax : 1.0);
            d->GetXaxis()->SetRangeUser(-60.0,120.0);
            d->Draw("E");
            c->Draw("HIST SAME");
            draw_vertical(signalMin, d->GetMaximum(), 2, 1);
            draw_vertical(signalMax, d->GetMaximum(), 2, 1);

            TLatex l;
            l.SetNDC();
            l.SetTextSize(0.045);
            l.DrawLatex(0.12,0.86,Form("core scale = %.3g",s));
            l.DrawLatex(0.12,0.79,Form("tail excess = %.0f",excess));

            if(ie == 0) {
                TLegend* le = new TLegend(0.54,0.68,0.89,0.87);
                le->SetBorderSize(0);
                le->AddEntry(d,"Data","lep");
                le->AddEntry(c,"Coherent MC","l");
                le->Draw();
            }

            fout->cd();
            d->Write();
            c->Write();
        }

        cBook->Print("data_vs_coherent_mm_by_energy.pdf");
    }

    cBook->Print("data_vs_coherent_mm_by_energy.pdf]");
    cBook->Write();

    fout->Close();
    txt.close();
    fc->Close();
    fd->Close();

    std::cout << "Coherent MC normalized to data in core [" << coreMin << "," << coreMax << "] MeV" << std::endl;
    std::cout << "GLOBAL data signal = " << dataSignal << ", scaled coherent MC signal = " << cohSignal << std::endl;
    std::cout << "GLOBAL data tail = " << dataTail << ", scaled coherent MC tail = " << cohTail
              << ", excess = " << tailExcess << std::endl;
    std::cout << "Saved:" << std::endl;
    std::cout << "  data_vs_coherent_mm_global.pdf" << std::endl;
    std::cout << "  data_vs_coherent_mm_by_energy.pdf" << std::endl;
    std::cout << "  data_vs_coherent_mm.root" << std::endl;
    std::cout << "  data_vs_coherent_mm.txt" << std::endl;
}
