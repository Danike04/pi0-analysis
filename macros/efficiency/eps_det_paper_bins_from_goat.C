/*
 * eps_det_paper_bins_from_goat.C
 *
 * PURPOSE
 *   Calculate the paper-bin detection efficiency using histograms already
 *   produced by the GoAT MC analysis rather than reconstructing AcquRoot trees
 *   event by event.
 *
 * DEFAULT HISTOGRAMS
 *   denominator: TaggerBinning/theta_MC
 *   numerator:   TaggerBinning/theta_MM
 *   If theta_MM is absent, the original macro falls back to
 *   TaggerBinning/theta_all.
 *
 * ENERGY MAPPING
 *   Tagger channel -> photon energy is read from the supplied FPD file. The
 *   default energyColumn is 15 and the default number of channels is 328.
 *   These values describe the historical GoAT/FPD format and should be checked
 *   when using a different configuration.
 *
 * EFFICIENCY
 *   For each paper photon-energy bin and theta_cm interval the macro integrates
 *   the generated and reconstructed GoAT histograms and calculates
 *
 *       eps_det(E,theta) = Nreco / Ngen.
 *
 * OUTPUTS
 *   eps_det_paper_bins_from_goat.txt
 *   eps_det_paper_bins_from_goat.root
 *   eps_det_paper_bins_from_goat.pdf
 *
 * The text output uses the format accepted by
 * paper_style_diff_xs_from_acqu.C in paper-bin efficiency mode.
 *
 * NOTE FOR THE REPOSITORY
 *   Analysis logic and defaults below are preserved from the original macro.
 */
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH2.h"
#include "TH2D.h"
#include "TMath.h"
#include "TString.h"

bool load_goat_egamma_map(const char* fpdFile,
                          double egammaByCh[],
                          int nCh,
                          int energyColumn)
{
    for(int ch=0; ch<nCh; ch++) egammaByCh[ch] = -1.0;
    if(!fpdFile || TString(fpdFile).Length() == 0) return false;

    std::ifstream in(fpdFile);
    if(!in.is_open()) return false;

    std::string line;
    while(std::getline(in, line)) {
        if(line.empty() || line[0] == '#') continue;

        std::istringstream ss(line);
        std::vector<std::string> tok;
        std::string item;
        while(ss >> item) tok.push_back(item);
        if(tok.size() < 2) continue;

        int ch = -1;
        if(tok[0] == "Element:" && tok.size() > 1) {
            ch = std::atoi(tok[1].c_str());
        } else {
            ch = std::atoi(tok[0].c_str());
        }
        if(ch < 0 || ch >= nCh) continue;

        int col = energyColumn - 1;
        if(col < 0 || col >= int(tok.size())) continue;

        char* endptr = 0;
        double egamma = std::strtod(tok[col].c_str(), &endptr);
        if(endptr == tok[col].c_str()) continue;
        if(egamma > 0.0) egammaByCh[ch] = egamma;
    }

    for(int ch=0; ch<nCh; ch++) {
        if(egammaByCh[ch] > 0.0) return true;
    }
    return false;
}

int goat_channel_axis(TH2* h, int nCh)
{
    if(!h) return 0;
    int nx = h->GetXaxis()->GetNbins();
    int ny = h->GetYaxis()->GetNbins();
    double xmin = h->GetXaxis()->GetXmin();
    double xmax = h->GetXaxis()->GetXmax();
    double ymin = h->GetYaxis()->GetXmin();
    double ymax = h->GetYaxis()->GetXmax();

    bool xLooksChannel = (nx >= nCh-5 && nx <= nCh+30) || (xmin <= 0.5 && xmax >= nCh-0.5);
    bool yLooksChannel = (ny >= nCh-5 && ny <= nCh+30) || (ymin <= 0.5 && ymax >= nCh-0.5);

    if(xLooksChannel && !yLooksChannel) return 1;
    if(yLooksChannel && !xLooksChannel) return 2;
    return 1;
}

double goat_integral_channel_theta(TH2* h,
                                   int ch,
                                   int channelAxis,
                                   double thetaLow,
                                   double thetaHigh,
                                   double& err)
{
    err = 0.0;
    if(!h) return 0.0;

    if(channelAxis == 1) {
        int bx = h->GetXaxis()->FindBin(ch + 0.5);
        int y1 = h->GetYaxis()->FindBin(thetaLow);
        int y2 = h->GetYaxis()->FindBin(thetaHigh);
        return h->IntegralAndError(bx, bx, y1, y2, err);
    }

    int by = h->GetYaxis()->FindBin(ch + 0.5);
    int x1 = h->GetXaxis()->FindBin(thetaLow);
    int x2 = h->GetXaxis()->FindBin(thetaHigh);
    return h->IntegralAndError(x1, x2, by, by, err);
}

double goat_solid_angle(double thetaLow, double thetaHigh)
{
    double th1 = thetaLow * TMath::DegToRad();
    double th2 = thetaHigh * TMath::DegToRad();
    return 2.0 * TMath::Pi() * (std::cos(th1) - std::cos(th2));
}

void eps_det_paper_bins_from_goat(const char* fname="He4Pi0MC.root",
                                  const char* fpdFile="FPD_855_new.dat",
                                  int energyColumn=15,
                                  int nCh=328,
                                  double thetaCmMinDeg=5.0,
                                  double thetaCmMaxDeg=150.0,
                                  double thetaCmBinWidthDeg=5.0,
                                  const char* genHistName="TaggerBinning/theta_MC",
                                  const char* recHistName="TaggerBinning/theta_MM")
{
    std::vector<double> eBinLow;
    std::vector<double> eBinHigh;
    const double eLowPaper[]  = {201, 211, 223, 235, 247, 259, 271, 283, 295, 309, 319, 331, 343, 356, 367, 379, 391};
    const double eHighPaper[] = {210, 222, 234, 246, 258, 270, 282, 294, 308, 318, 330, 342, 355, 366, 378, 390, 401};
    const int nPaperEBins = 17;
    for(int i=0; i<nPaperEBins; i++) {
        eBinLow.push_back(eLowPaper[i]);
        eBinHigh.push_back(eHighPaper[i]);
    }
    int nEBins = int(eBinLow.size());

    if(thetaCmBinWidthDeg <= 0.0) thetaCmBinWidthDeg = 5.0;
    int nThetaBins = int(std::floor((thetaCmMaxDeg - thetaCmMinDeg) / thetaCmBinWidthDeg + 0.5));
    if(nThetaBins < 1) nThetaBins = 1;

    TFile* f = TFile::Open(fname);
    if(!f || f->IsZombie()) {
        std::cout << "Cannot open " << fname << std::endl;
        return;
    }

    TH2* hGen = (TH2*)f->Get(genHistName);
    TH2* hRec = (TH2*)f->Get(recHistName);
    if(!hRec) hRec = (TH2*)f->Get("TaggerBinning/theta_all");

    if(!hGen || !hRec) {
        std::cout << "Missing required histograms:" << std::endl;
        std::cout << "  gen: " << genHistName << " -> " << hGen << std::endl;
        std::cout << "  rec: " << recHistName << " -> " << hRec << std::endl;
        f->Close();
        return;
    }

    double egammaByCh[400];
    bool haveEgammaMap = load_goat_egamma_map(fpdFile, egammaByCh, 400, energyColumn);
    if(!haveEgammaMap) {
        std::cout << "Could not read E_gamma map from " << fpdFile << std::endl;
        f->Close();
        return;
    }

    int genAxis = goat_channel_axis(hGen, nCh);
    int recAxis = goat_channel_axis(hRec, nCh);

    TH2D* hEps = new TH2D("h_eps_det_paper_bins_from_goat",
                          "#epsilon_{det}(E_{#gamma}^{lab},#theta^{cm}) from GoAT;E_{#gamma}^{lab} bin;#theta_{#pi^{0}}^{cm} (deg)",
                          nEBins, 0.0, double(nEBins),
                          nThetaBins, thetaCmMinDeg, thetaCmMinDeg + nThetaBins * thetaCmBinWidthDeg);
    for(int ie=0; ie<nEBins; ie++) {
        hEps->GetXaxis()->SetBinLabel(ie + 1, Form("%.0f-%.0f", eBinLow[ie], eBinHigh[ie]));
    }

    std::ofstream out("eps_det_paper_bins_from_goat.txt");
    out << "# input: " << fname << "\n";
    out << "# denominator: " << genHistName << "\n";
    out << "# numerator: " << hRec->GetName() << "\n";
    out << "# E_low E_high E_center thetaCm_low thetaCm_high thetaCm_center DeltaOmega_sr Ngen_E Ngen_bin Nreco eps_det deps_det\n";

    for(int ie=0; ie<nEBins; ie++) {
        for(int ith=0; ith<nThetaBins; ith++) {
            double th1 = thetaCmMinDeg + ith * thetaCmBinWidthDeg;
            double th2 = th1 + thetaCmBinWidthDeg;
            double thCenter = 0.5 * (th1 + th2);
            double deltaOmega = goat_solid_angle(th1, th2);

            double nGen = 0.0;
            double nRec = 0.0;
            double dGen2 = 0.0;
            double dRec2 = 0.0;
            double nGenE = 0.0;

            for(int ch=0; ch<nCh; ch++) {
                double eg = egammaByCh[ch];
                if(eg < eBinLow[ie] || eg > eBinHigh[ie]) continue;

                double dg = 0.0;
                double dr = 0.0;
                double ngTheta = goat_integral_channel_theta(hGen, ch, genAxis, th1, th2, dg);
                double nrTheta = goat_integral_channel_theta(hRec, ch, recAxis, th1, th2, dr);
                double ngAll = goat_integral_channel_theta(hGen, ch, genAxis, thetaCmMinDeg, thetaCmMaxDeg, dg);

                nGen += ngTheta;
                nRec += nrTheta;
                dGen2 += dg * dg;
                dRec2 += dr * dr;
                nGenE += ngAll;
            }

            double eps = 0.0;
            double deps = 0.0;
            if(nGen > 0.0) {
                eps = nRec / nGen;
                if(eps >= 0.0 && eps <= 1.0) {
                    deps = std::sqrt(eps * (1.0 - eps) / nGen);
                } else {
                    double rel2 = 0.0;
                    if(nRec > 0.0) rel2 += dRec2 / (nRec * nRec);
                    if(nGen > 0.0) rel2 += dGen2 / (nGen * nGen);
                    deps = std::fabs(eps) * std::sqrt(rel2);
                }
                hEps->SetBinContent(ie + 1, ith + 1, eps);
                hEps->SetBinError(ie + 1, ith + 1, deps);
            }

            out << eBinLow[ie] << " "
                << eBinHigh[ie] << " "
                << 0.5 * (eBinLow[ie] + eBinHigh[ie]) << " "
                << th1 << " "
                << th2 << " "
                << thCenter << " "
                << deltaOmega << " "
                << nGenE << " "
                << nGen << " "
                << nRec << " "
                << eps << " "
                << deps << "\n";
        }
    }
    out.close();

    TCanvas* c = new TCanvas("c_eps_det_paper_bins_from_goat", "eps det paper bins from GoAT", 1200, 750);
    hEps->GetXaxis()->LabelsOption("v");
    hEps->Draw("COLZ");
    c->SaveAs("eps_det_paper_bins_from_goat.pdf");

    TFile* fout = new TFile("eps_det_paper_bins_from_goat.root", "RECREATE");
    hEps->Write();
    c->Write();
    fout->Close();

    std::cout << "Saved:" << std::endl;
    std::cout << "  eps_det_paper_bins_from_goat.txt" << std::endl;
    std::cout << "  eps_det_paper_bins_from_goat.pdf" << std::endl;
    std::cout << "  eps_det_paper_bins_from_goat.root" << std::endl;

    f->Close();
}
