/*
 * add_table1_first6_overlay_to_root.C
 *
 * PURPOSE
 *   Post-process the ROOT output of paper_style_diff_xs_from_acqu.C by adding
 *   graphs from the published/reference Table-1 differential cross sections
 *   and producing comparison plots in the same CM energy/angle binning.
 *
 * INPUTS
 *   - paper_style_dsigma_dOmega_cm_from_acqu.root by default.
 *   - paper_table1_dsigma_dOmega.txt by default, with columns
 *       E_low E_high E_center theta_cm value uncertainty.
 *
 * BEHAVIOUR
 *   Despite the historical filename, the implementation contains
 *   add_table1_all_overlay_to_root(), which loops over all available paper
 *   energy bins. add_table1_first6_overlay_to_root() simply calls that routine.
 *   The input ROOT file is opened in UPDATE mode and new graphs/canvases are
 *   written into it.
 *
 * USAGE EXAMPLE
 *   root -l -b -q 'add_table1_first6_overlay_to_root.C("paper_style_dsigma_dOmega_cm_from_acqu.root","paper_table1_dsigma_dOmega.txt")'
 *
 * REPOSITORY STATUS
 *   Useful paper-comparison study. Executable code below is unchanged.
 */
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TObject.h"
#include "TString.h"

TGraphErrors* load_table1_graph_overlay(const char* tableFile, double eLow, double eHigh)
{
    std::ifstream in(tableFile);
    if(!in.is_open()) return 0;

    TGraphErrors* gr = new TGraphErrors();
    gr->SetName(Form("table1_dsigma_dOmega_cm_E_%.0f_%.0f", eLow, eHigh));
    gr->SetTitle(Form("Table 1 %.0f < E_{#gamma}^{lab} < %.0f MeV;#theta^{cm}_{#pi^{0}} (deg);d#sigma/d#Omega (#mub/sr)", eLow, eHigh));
    gr->SetMarkerStyle(24);
    gr->SetMarkerSize(0.9);
    gr->SetMarkerColor(kRed+1);
    gr->SetLineColor(kRed+1);

    std::string line;
    while(std::getline(in, line)) {
        if(line.empty() || line[0] == '#') continue;

        std::istringstream ss(line);
        double el = 0.0;
        double eh = 0.0;
        double ec = 0.0;
        double th = 0.0;
        double y = 0.0;
        double ey = 0.0;
        if(!(ss >> el >> eh >> ec >> th >> y >> ey)) continue;
        if(std::fabs(el - eLow) > 0.1 || std::fabs(eh - eHigh) > 0.1) continue;

        int ip = gr->GetN();
        gr->SetPoint(ip, th, y);
        gr->SetPointError(ip, 0.0, ey);
    }

    if(gr->GetN() == 0) {
        delete gr;
        return 0;
    }

    return gr;
}

TGraphErrors* graph_from_projection(TH1D* h, const char* name)
{
    if(!h) return 0;

    TGraphErrors* gr = new TGraphErrors();
    gr->SetName(name);
    gr->SetTitle(Form("%s;#theta^{cm}_{#pi^{0}} (deg);d#sigma/d#Omega (#mub/sr)", name));
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0.9);
    gr->SetMarkerColor(kBlue+1);
    gr->SetLineColor(kBlue+1);

    for(int ib=1; ib<=h->GetNbinsX(); ib++) {
        double y = h->GetBinContent(ib);
        double ey = h->GetBinError(ib);
        if(y == 0.0 && ey == 0.0) continue;

        int ip = gr->GetN();
        gr->SetPoint(ip, h->GetBinCenter(ib), y);
        gr->SetPointError(ip, 0.5 * h->GetBinWidth(ib), ey);
    }

    return gr;
}

void add_table1_all_overlay_to_root(const char* rootFile="paper_style_dsigma_dOmega_cm_from_acqu.root",
                                    const char* tableFile="paper_table1_dsigma_dOmega.txt")
{
    const double eLow[17]  = {201, 211, 223, 235, 247, 259, 271, 283, 295, 309, 319, 331, 343, 356, 367, 379, 391};
    const double eHigh[17] = {210, 222, 234, 246, 258, 270, 282, 294, 308, 318, 330, 342, 355, 366, 378, 390, 401};

    TFile* f = TFile::Open(rootFile, "UPDATE");
    if(!f || f->IsZombie()) {
        std::cout << "Cannot open " << rootFile << std::endl;
        return;
    }

    TH2D* hMap = (TH2D*)f->Get("h_dsigma_dOmega_cm_paper_bins");
    if(!hMap) {
        std::cout << "Missing h_dsigma_dOmega_cm_paper_bins in " << rootFile << std::endl;
        f->Close();
        return;
    }

    const int nCompareBins = std::min(17, hMap->GetNbinsX());
    TGraphErrors* grAcqu[17];
    TGraphErrors* grTable1[17];
    for(int i=0; i<17; i++) {
        grAcqu[i] = 0;
        grTable1[i] = 0;
    }

    TH2D* hTable1Map = (TH2D*)hMap->Clone("h_dsigma_dOmega_cm_table1_all");
    hTable1Map->Reset("ICES");
    hTable1Map->SetTitle("Table 1 d#sigma/d#Omega in CM, all energy bins;E_{#gamma}^{lab} bin;#theta^{cm}_{#pi^{0}} (deg)");

    TCanvas* c = new TCanvas("c_paper_style_dsigma_dOmega_cm_all_overlay",
                             "paper-style bins with Table 1 overlay", 1100, 700);
    c->Print("paper_style_dsigma_dOmega_cm_all_overlay.pdf[");

    for(int ie=0; ie<nCompareBins; ie++) {
        TH1D* hProj = hMap->ProjectionY(Form("projY_dsigma_dOmega_cm_Ebin%d", ie + 1), ie + 1, ie + 1);
        grAcqu[ie] = graph_from_projection(hProj, Form("acqu_dsigma_dOmega_cm_Ebin%d", ie + 1));
        grTable1[ie] = load_table1_graph_overlay(tableFile, eLow[ie], eHigh[ie]);

        if(grTable1[ie]) {
            for(int ip=0; ip<grTable1[ie]->GetN(); ip++) {
                double x = 0.0;
                double y = 0.0;
                grTable1[ie]->GetPoint(ip, x, y);
                int by = hTable1Map->GetYaxis()->FindBin(x);
                hTable1Map->SetBinContent(ie + 1, by, y);
                hTable1Map->SetBinError(ie + 1, by, grTable1[ie]->GetErrorY(ip));
            }
        }

        c->Clear();
        c->SetGrid();
        if(grAcqu[ie]) grAcqu[ie]->Draw("AP");
        if(grTable1[ie]) {
            if(grAcqu[ie]) grTable1[ie]->Draw("P SAME");
            else grTable1[ie]->Draw("AP");
        }

        TLegend* leg = new TLegend(0.52, 0.72, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        if(grAcqu[ie]) leg->AddEntry(grAcqu[ie], "ACQU", "p");
        if(grTable1[ie]) leg->AddEntry(grTable1[ie], "Paper Table 1", "p");
        leg->Draw();
        c->Print("paper_style_dsigma_dOmega_cm_all_overlay.pdf");
    }
    c->Print("paper_style_dsigma_dOmega_cm_all_overlay.pdf]");

    TCanvas* cMap = new TCanvas("c_paper_table1_dsigma_dOmega_cm_map", "Table 1 dsigma/dOmega CM map", 1200, 750);
    hTable1Map->GetXaxis()->LabelsOption("v");
    hTable1Map->Draw("COLZ");
    cMap->SaveAs("paper_table1_dsigma_dOmega_cm_map.pdf");

    hTable1Map->Write("", TObject::kOverwrite);
    for(int ie=0; ie<nCompareBins; ie++) {
        if(grAcqu[ie]) grAcqu[ie]->Write("", TObject::kOverwrite);
        if(grTable1[ie]) grTable1[ie]->Write("", TObject::kOverwrite);
    }
    c->Write("", TObject::kOverwrite);
    cMap->Write("", TObject::kOverwrite);
    f->Close();

    std::cout << "Updated " << rootFile << std::endl;
    std::cout << "Saved paper_style_dsigma_dOmega_cm_all_overlay.pdf" << std::endl;
    std::cout << "Saved paper_table1_dsigma_dOmega_cm_map.pdf" << std::endl;
}

void add_table1_first6_overlay_to_root(const char* rootFile="paper_style_dsigma_dOmega_cm_from_acqu.root",
                                       const char* tableFile="paper_table1_dsigma_dOmega.txt")
{
    add_table1_all_overlay_to_root(rootFile, tableFile);
}
