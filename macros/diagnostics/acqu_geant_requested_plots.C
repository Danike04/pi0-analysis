/*
 * acqu_geant_requested_plots.C
 *
 * Purpose
 * -------
 * Produce reconstruction-level diagnostic plots from an Acqu/Geant coherent
 * He4(gamma,pi0)He4 Monte Carlo sample. It was used to inspect the same
 * reconstructed observables and cuts later used in the data analysis.
 *
 * Input
 * -----
 * An Acqu/Geant ROOT file containing `tracks`, `tagger` and
 * `setupParameters`. The default filename is Acqu_geant_He4pi0.root.
 *
 * Main default selections
 * -----------------------
 * cluster energy > 20 MeV
 * veto energy < 1 MeV
 * 110 <= m(gamma gamma) <= 155 MeV
 * -10 <= missing mass <= 40 MeV relative to the He4 mass
 * theta_pi0^CM binning: 5--150 degrees in 5-degree bins
 *
 * What is plotted
 * ---------------
 * Among other diagnostics, the macro produces m(gamma gamma), missing mass,
 * E_gamma-vs-angle maps, lab/CM pion-angle distributions, gamma-gamma
 * opening angle versus the two-body minimum opening angle, and reconstructed
 * counts in the paper energy/CM-angle bins before/after selections.
 *
 * Output
 * ------
 * A collection of PDF/PNG diagnostic plots and
 * `acqu_geant_requested_plots.root` containing the histograms/graphs.
 *
 * Usage example
 * -------------
 * root -l -b -q 'acqu_geant_requested_plots.C("Acqu_geant_He4pi0.root")'
 *
 * This is a Monte-Carlo validation/diagnostic macro, not the detector-
 * efficiency extraction itself. The original code below is unchanged; only
 * this documentation header has been added.
 */

#include <cmath>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TStyle.h"
#include "TTree.h"
#include "TVirtualPad.h"

int find_paper_energy_bin(double x, const std::vector<double>& low, const std::vector<double>& high)
{
    for(size_t i=0; i<low.size(); i++) {
        if(x >= low[i] && x <= high[i]) return int(i);
    }
    return -1;
}

double coherent_pion_lab_energy_geant(double egamma, double thetaCmDeg, double mpi0, double targetMass)
{
    double s = targetMass * targetMass + 2.0 * egamma * targetMass;
    double sqrtS = TMath::Sqrt(s);
    double ePiCm = (s + mpi0 * mpi0 - targetMass * targetMass) / (2.0 * sqrtS);
    double pPiCm2 = ePiCm * ePiCm - mpi0 * mpi0;
    if(pPiCm2 < 0.0) pPiCm2 = 0.0;
    double pPiCm = TMath::Sqrt(pPiCm2);

    double betaCm = egamma / (egamma + targetMass);
    double gammaCm = 1.0 / TMath::Sqrt(1.0 - betaCm * betaCm);
    double cosThetaCm = TMath::Cos(thetaCmDeg * TMath::DegToRad());

    return gammaCm * (ePiCm + betaCm * pPiCm * cosThetaCm);
}

double phi_min_opening_angle_geant(double ePiLab, double mpi0)
{
    if(ePiLab <= mpi0) return 180.0;

    double arg = TMath::Sqrt(ePiLab * ePiLab - mpi0 * mpi0) / ePiLab;
    if(arg < -1.0) arg = -1.0;
    if(arg > 1.0) arg = 1.0;

    return 2.0 * TMath::ACos(arg) * TMath::RadToDeg();
}

void acqu_geant_requested_plots(const char* fname="Acqu_geant_He4pi0.root",
                                Long64_t maxEvents=-1,
                                double eMin=20.0,
                                double vetoMax=1.0,
                                int maxNeutralForPairs=16,
                                Long64_t printEvery=50000,
                                double mggMin=110.0,
                                double mggMax=155.0,
                                double mmMin=-10.0,
                                double mmMax=40.0,
                                double thetaCmMinDeg=5.0,
                                double thetaCmMaxDeg=150.0,
                                double thetaCmBinWidthDeg=5.0)
{
    const int nChMax = 352;
    const int maxTracks = 512;
    const int maxTagged = 2048;

    const double beamEMax = 855.0;
    const double targetMass = 3727.38;
    const double mpi0 = 134.9768;

    TFile* f = TFile::Open(fname);
    if(!f || f->IsZombie()) {
        std::cout << "Cannot open " << fname << std::endl;
        return;
    }

    TTree* tracks = (TTree*)f->Get("tracks");
    TTree* tagger = (TTree*)f->Get("tagger");
    TTree* setup = (TTree*)f->Get("setupParameters");

    if(!tracks || !tagger || !setup) {
        std::cout << "Missing tracks/tagger/setupParameters" << std::endl;
        f->Close();
        return;
    }

    int nTracks = 0;
    double clusterEnergy[maxTracks];
    double theta[maxTracks];
    double phi[maxTracks];
    double vetoEnergy[maxTracks];

    tracks->SetBranchStatus("*", 0);
    tracks->SetBranchStatus("nTracks", 1);
    tracks->SetBranchStatus("clusterEnergy", 1);
    tracks->SetBranchStatus("theta", 1);
    tracks->SetBranchStatus("phi", 1);
    tracks->SetBranchStatus("vetoEnergy", 1);
    tracks->SetBranchAddress("nTracks", &nTracks);
    tracks->SetBranchAddress("clusterEnergy", clusterEnergy);
    tracks->SetBranchAddress("theta", theta);
    tracks->SetBranchAddress("phi", phi);
    tracks->SetBranchAddress("vetoEnergy", vetoEnergy);

    int nTagged = 0;
    int taggedChannel[maxTagged];

    tagger->SetBranchStatus("*", 0);
    tagger->SetBranchStatus("nTagged", 1);
    tagger->SetBranchStatus("taggedChannel", 1);
    tagger->SetBranchAddress("nTagged", &nTagged);
    tagger->SetBranchAddress("taggedChannel", taggedChannel);

    int nTagger = 0;
    double taggerPhotonEnergy[nChMax];

    setup->SetBranchStatus("*", 0);
    setup->SetBranchStatus("nTagger", 1);
    setup->SetBranchStatus("TaggerPhotonEnergy", 1);
    setup->SetBranchAddress("nTagger", &nTagger);
    setup->SetBranchAddress("TaggerPhotonEnergy", taggerPhotonEnergy);
    setup->GetEntry(0);

    int nCh = (nTagger < nChMax) ? nTagger : nChMax;

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
    if(thetaCmMaxDeg <= thetaCmMinDeg) {
        thetaCmMinDeg = 5.0;
        thetaCmMaxDeg = 150.0;
    }
    int nThetaBins = int(TMath::Floor((thetaCmMaxDeg - thetaCmMinDeg) / thetaCmBinWidthDeg + 0.5));
    if(nThetaBins < 1) nThetaBins = 1;

    TH1D* hMgg = new TH1D("h_mgg_best_pair",
                          "Best #gamma#gamma pair;m_{#gamma#gamma} (MeV);Counts",
                          400, 0.0, 400.0);

    TH1D* hMM = new TH1D("h_mm_best_pair",
                         "Missing mass using best #gamma#gamma pair;MM-M_{^{4}He} (MeV);Counts",
                         300, -150.0, 150.0);

    TH2D* hMMvsMgg = new TH2D("h_mm_vs_mgg",
                              "MM vs m_{#gamma#gamma};m_{#gamma#gamma} (MeV);MM-M_{^{4}He} (MeV)",
                              250, 0.0, 250.0,
                              300, -150.0, 150.0);

    TH2D* hEgammaVsMM = new TH2D("h_egamma_vs_mm",
                                 "E_{#gamma} vs MM;MM-M_{^{4}He} (MeV);E_{#gamma} (MeV)",
                                 300, -150.0, 150.0,
                                 352, 0.0, beamEMax);

    TH2D* hEgammaVsTheta = new TH2D("h_egamma_vs_theta",
                                    "E_{#gamma} vs #theta_{#pi^{0}};#theta_{#pi^{0}} (deg);E_{#gamma} (MeV)",
                                    180, 0.0, 180.0,
                                    352, 0.0, beamEMax);

    TH2D* hEgammaVsThetaCmBefore = new TH2D("h_egamma_vs_theta_cm_before_mgg_cut",
                                            "E_{#gamma} vs #theta^{cm}_{#pi^{0}} before m_{#gamma#gamma} cut;#theta^{cm}_{#pi^{0}} (deg);E_{#gamma} (MeV)",
                                            180, 0.0, 180.0,
                                            352, 0.0, beamEMax);

    TH2D* hEgammaVsThetaCmAfter = new TH2D("h_egamma_vs_theta_cm_after_mgg_cut",
                                           "E_{#gamma} vs #theta^{cm}_{#pi^{0}} after m_{#gamma#gamma} cut;#theta^{cm}_{#pi^{0}} (deg);E_{#gamma} (MeV)",
                                           180, 0.0, 180.0,
                                           352, 0.0, beamEMax);

    TH2D* hOpeningPhiMinBefore = new TH2D("h_opening_angle_vs_phi_min_before_mgg_cut",
                                          "GEANT #Phi_{#gamma#gamma} vs #Phi_{min} before m_{#gamma#gamma} cut;#Phi_{min} Eq. 2 (deg);measured #Phi_{#gamma#gamma} (deg)",
                                          180, 0.0, 180.0,
                                          180, 0.0, 180.0);

    TH2D* hOpeningPhiMinAfter = new TH2D("h_opening_angle_vs_phi_min_after_mgg_cut",
                                         "GEANT #Phi_{#gamma#gamma} vs #Phi_{min} after m_{#gamma#gamma} cut;#Phi_{min} Eq. 2 (deg);measured #Phi_{#gamma#gamma} (deg)",
                                         180, 0.0, 180.0,
                                         180, 0.0, 180.0);

    TH2D* hOpeningMinusPhiMinAfter = new TH2D("h_opening_minus_phi_min_vs_egamma_after_mgg_cut",
                                              "GEANT #Phi_{#gamma#gamma}-#Phi_{min} after m_{#gamma#gamma} cut;E_{#gamma} (MeV);#Phi_{#gamma#gamma}-#Phi_{min} (deg)",
                                              352, 0.0, beamEMax,
                                              240, -60.0, 60.0);

    TH2D* hOpeningEgammaBefore = new TH2D("h_geant_phigg_vs_egamma_before_mgg_cut",
                                          "GEANT #Phi_{#gamma#gamma} vs E_{#gamma} before m_{#gamma#gamma} cut;E_{#gamma} (MeV);#Phi_{#gamma#gamma} (deg)",
                                          352, 0.0, beamEMax,
                                          180, 0.0, 180.0);

    TH2D* hOpeningEgammaAfter = new TH2D("h_geant_phigg_vs_egamma_after_mgg_cut",
                                         "GEANT #Phi_{#gamma#gamma} vs E_{#gamma} after m_{#gamma#gamma} cut;E_{#gamma} (MeV);#Phi_{#gamma#gamma} (deg)",
                                         352, 0.0, beamEMax,
                                         180, 0.0, 180.0);

    TH2D* hOpeningEgammaAfterMM = new TH2D("h_geant_phigg_vs_egamma_after_mgg_mm_cut",
                                           "GEANT #Phi_{#gamma#gamma} vs E_{#gamma} after m_{#gamma#gamma} and MM cuts;E_{#gamma} (MeV);#Phi_{#gamma#gamma} (deg)",
                                           352, 0.0, beamEMax,
                                           180, 0.0, 180.0);

    TH2D* hRecoPaperBefore = new TH2D("h_geant_reco_counts_paper_bins_before_mgg_cut",
                                      "GEANT reco counts in paper bins before m_{#gamma#gamma} cut;E_{#gamma}^{lab} bin;#theta^{cm}_{#pi^{0}} (deg)",
                                      nEBins, 0.0, double(nEBins),
                                      nThetaBins, thetaCmMinDeg, thetaCmMinDeg + nThetaBins * thetaCmBinWidthDeg);

    TH2D* hRecoPaperAfter = new TH2D("h_geant_reco_counts_paper_bins_after_mgg_cut",
                                     "GEANT reco counts in paper bins after m_{#gamma#gamma} cut;E_{#gamma}^{lab} bin;#theta^{cm}_{#pi^{0}} (deg)",
                                     nEBins, 0.0, double(nEBins),
                                     nThetaBins, thetaCmMinDeg, thetaCmMinDeg + nThetaBins * thetaCmBinWidthDeg);

    TH2D* hRecoPaperAfterMM = new TH2D("h_geant_reco_counts_paper_bins_after_mgg_mm_cut",
                                       "GEANT reco counts in paper bins after m_{#gamma#gamma} and MM cuts;E_{#gamma}^{lab} bin;#theta^{cm}_{#pi^{0}} (deg)",
                                       nEBins, 0.0, double(nEBins),
                                       nThetaBins, thetaCmMinDeg, thetaCmMinDeg + nThetaBins * thetaCmBinWidthDeg);

    TH2D* hPhiMinMap = new TH2D("h_phi_min_eq2_cm_bins",
                                "GEANT #Phi_{min} from Eq. 2 using coherent two-body kinematics;E_{#gamma}^{lab} bin;#theta^{cm}_{#pi^{0}} (deg)",
                                nEBins, 0.0, double(nEBins),
                                nThetaBins, thetaCmMinDeg, thetaCmMinDeg + nThetaBins * thetaCmBinWidthDeg);

    std::vector<TGraphErrors*> grPhiMin(nEBins, 0);
    std::vector<TGraphErrors*> grGeantCounts(nEBins, 0);
    for(int ie=0; ie<nEBins; ie++) {
        grPhiMin[ie] = new TGraphErrors();
        grPhiMin[ie]->SetName(Form("phi_min_eq2_E_%.0f_%.0f", eBinLow[ie], eBinHigh[ie]));
        grPhiMin[ie]->SetTitle(Form("GEANT #Phi_{min}, %.0f < E_{#gamma}^{lab} < %.0f MeV;#theta^{cm}_{#pi^{0}} (deg);#Phi_{min} Eq. 2 (deg)", eBinLow[ie], eBinHigh[ie]));
        grPhiMin[ie]->SetMarkerStyle(21);
        grPhiMin[ie]->SetMarkerColor(kMagenta+2);
        grPhiMin[ie]->SetLineColor(kMagenta+2);

        grGeantCounts[ie] = new TGraphErrors();
        grGeantCounts[ie]->SetName(Form("geant_counts_after_mgg_mm_E_%.0f_%.0f", eBinLow[ie], eBinHigh[ie]));
        grGeantCounts[ie]->SetTitle(Form("GEANT %.0f < E_{#gamma}^{lab} < %.0f MeV;#theta^{cm}_{#pi^{0}} (deg);Counts", eBinLow[ie], eBinHigh[ie]));
        grGeantCounts[ie]->SetMarkerStyle(20);
        grGeantCounts[ie]->SetMarkerColor(kBlue+1);
        grGeantCounts[ie]->SetLineColor(kBlue+1);

        for(int ith=0; ith<nThetaBins; ith++) {
            double th1 = thetaCmMinDeg + ith * thetaCmBinWidthDeg;
            double th2 = th1 + thetaCmBinWidthDeg;
            double thCenter = 0.5 * (th1 + th2);
            double eCenter = 0.5 * (eBinLow[ie] + eBinHigh[ie]);
            double ePiLabCoh = coherent_pion_lab_energy_geant(eCenter, thCenter, mpi0, targetMass);
            double phiMinDeg = phi_min_opening_angle_geant(ePiLabCoh, mpi0);

            hPhiMinMap->SetBinContent(ie + 1, ith + 1, phiMinDeg);
            grPhiMin[ie]->SetPoint(ith, thCenter, phiMinDeg);
            grPhiMin[ie]->SetPointError(ith, 0.5 * thetaCmBinWidthDeg, 0.0);
        }
    }

    Long64_t nTrackEntries = tracks->GetEntries();
    Long64_t nTaggerEntries = tagger->GetEntries();
    Long64_t nentries = (nTrackEntries < nTaggerEntries) ? nTrackEntries : nTaggerEntries;
    if(maxEvents > 0 && maxEvents < nentries) nentries = maxEvents;

    const double deg = TMath::DegToRad();
    TLorentzVector target(0.0, 0.0, 0.0, targetMass);

    Long64_t nWithTag = 0;
    Long64_t nWithTwoNeutral = 0;
    Long64_t nFilled = 0;

    for(Long64_t ev=0; ev<nentries; ev++) {
        if(printEvery > 0 && (ev == 0 || ev % printEvery == 0)) {
            std::cout << "event " << ev << "/" << nentries
                      << " (" << 100.0 * double(ev) / double(nentries) << "%)"
                      << " filled=" << nFilled << std::endl;
        }

        tracks->GetEntry(ev);
        tagger->GetEntry(ev);

        bool seenCh[nChMax] = {false};
        int eventChannels[nChMax];
        int nEventChannels = 0;

        int ntag = (nTagged < maxTagged) ? nTagged : maxTagged;
        for(int it=0; it<ntag; it++) {
            int ch = taggedChannel[it];
            if(ch < 0 || ch >= nCh) continue;
            if(seenCh[ch]) continue;
            seenCh[ch] = true;
            eventChannels[nEventChannels] = ch;
            nEventChannels++;
        }

        if(nEventChannels == 0) continue;
        nWithTag++;

        int neutralIdx[maxTracks];
        int nNeutral = 0;
        int nt = (nTracks < maxTracks) ? nTracks : maxTracks;
        for(int i=0; i<nt; i++) {
            if(clusterEnergy[i] <= eMin) continue;
            if(vetoEnergy[i] >= vetoMax) continue;
            neutralIdx[nNeutral] = i;
            nNeutral++;
        }

        if(nNeutral < 2) continue;
        nWithTwoNeutral++;

        for(int a=0; a<nNeutral-1; a++) {
            for(int b=a+1; b<nNeutral; b++) {
                if(clusterEnergy[neutralIdx[b]] > clusterEnergy[neutralIdx[a]]) {
                    int tmp = neutralIdx[a];
                    neutralIdx[a] = neutralIdx[b];
                    neutralIdx[b] = tmp;
                }
            }
        }

        int nUse = (maxNeutralForPairs > 0)
                 ? ((nNeutral < maxNeutralForPairs) ? nNeutral : maxNeutralForPairs)
                 : nNeutral;

        double bestM = -1.0;
        double bestDist = 1.0e30;
        double bestOpeningDeg = -1.0;
        TLorentzVector bestPi0;

        for(int a=0; a<nUse; a++) {
            int i = neutralIdx[a];
            for(int b=a+1; b<nUse; b++) {
                int j = neutralIdx[b];

                double th1 = theta[i] * deg;
                double th2 = theta[j] * deg;
                double ph1 = phi[i] * deg;
                double ph2 = phi[j] * deg;

                TLorentzVector g1;
                TLorentzVector g2;
                g1.SetPxPyPzE(clusterEnergy[i]*TMath::Sin(th1)*TMath::Cos(ph1),
                              clusterEnergy[i]*TMath::Sin(th1)*TMath::Sin(ph1),
                              clusterEnergy[i]*TMath::Cos(th1),
                              clusterEnergy[i]);
                g2.SetPxPyPzE(clusterEnergy[j]*TMath::Sin(th2)*TMath::Cos(ph2),
                              clusterEnergy[j]*TMath::Sin(th2)*TMath::Sin(ph2),
                              clusterEnergy[j]*TMath::Cos(th2),
                              clusterEnergy[j]);

                TLorentzVector pi0 = g1 + g2;
                double mgg = pi0.M();
                if(mgg <= 0.0) continue;

                double dist = TMath::Abs(mgg - mpi0);
                if(dist < bestDist) {
                    bestDist = dist;
                    bestM = mgg;
                    bestPi0 = pi0;
                    bestOpeningDeg = g1.Angle(g2.Vect()) * TMath::RadToDeg();
                }
            }
        }

        if(bestM <= 0.0) continue;

        double thetaPi0Deg = bestPi0.Theta() * TMath::RadToDeg();
        bool passMggCut = (bestM >= mggMin && bestM <= mggMax);
        hMgg->Fill(bestM);

        for(int k=0; k<nEventChannels; k++) {
            int ch = eventChannels[k];
            double egamma = taggerPhotonEnergy[ch];
            if(egamma <= 0.0 || egamma > beamEMax + 5.0) continue;

            TLorentzVector beam(0.0, 0.0, egamma, egamma);
            TLorentzVector recoil = beam + target - bestPi0;
            double mm = recoil.M() - targetMass;
            bool passMMCut = (mm >= mmMin && mm <= mmMax);

            TLorentzVector pi0CM = bestPi0;
            TLorentzVector totalCM = beam + target;
            pi0CM.Boost(-totalCM.BoostVector());
            double thetaCmDeg = pi0CM.Theta() * TMath::RadToDeg();
            double ePiLabCoh = coherent_pion_lab_energy_geant(egamma, thetaCmDeg, mpi0, targetMass);
            double phiMinDeg = phi_min_opening_angle_geant(ePiLabCoh, mpi0);
            double openingMinusPhiMin = bestOpeningDeg - phiMinDeg;

            hMM->Fill(mm);
            hMMvsMgg->Fill(bestM, mm);
            hEgammaVsMM->Fill(mm, egamma);
            hEgammaVsTheta->Fill(thetaPi0Deg, egamma);

            hEgammaVsThetaCmBefore->Fill(thetaCmDeg, egamma);
            if(bestOpeningDeg >= 0.0) hOpeningPhiMinBefore->Fill(phiMinDeg, bestOpeningDeg);
            if(bestOpeningDeg >= 0.0) hOpeningEgammaBefore->Fill(egamma, bestOpeningDeg);

            int ie = find_paper_energy_bin(egamma, eBinLow, eBinHigh);
            if(ie >= 0 && thetaCmDeg >= thetaCmMinDeg && thetaCmDeg < thetaCmMaxDeg) {
                hRecoPaperBefore->Fill(double(ie) + 0.5, thetaCmDeg);
            }

            if(passMggCut) {
                hEgammaVsThetaCmAfter->Fill(thetaCmDeg, egamma);
                if(bestOpeningDeg >= 0.0) {
                    hOpeningPhiMinAfter->Fill(phiMinDeg, bestOpeningDeg);
                    hOpeningMinusPhiMinAfter->Fill(egamma, openingMinusPhiMin);
                    hOpeningEgammaAfter->Fill(egamma, bestOpeningDeg);
                    if(passMMCut) hOpeningEgammaAfterMM->Fill(egamma, bestOpeningDeg);
                }
                if(ie >= 0 && thetaCmDeg >= thetaCmMinDeg && thetaCmDeg < thetaCmMaxDeg) {
                    hRecoPaperAfter->Fill(double(ie) + 0.5, thetaCmDeg);
                    if(passMMCut) hRecoPaperAfterMM->Fill(double(ie) + 0.5, thetaCmDeg);
                }
            }

            nFilled++;
        }
    }

    for(int ie=0; ie<nEBins; ie++) {
        for(int ith=0; ith<nThetaBins; ith++) {
            double y = hRecoPaperAfterMM->GetBinContent(ie + 1, ith + 1);
            double ey = (y > 0.0) ? TMath::Sqrt(y) : 0.0;
            if(y <= 0.0 && ey <= 0.0) continue;

            double thCenter = hRecoPaperAfterMM->GetYaxis()->GetBinCenter(ith + 1);
            int ip = grGeantCounts[ie]->GetN();
            grGeantCounts[ie]->SetPoint(ip, thCenter, y);
            grGeantCounts[ie]->SetPointError(ip, 0.5 * thetaCmBinWidthDeg, ey);
        }
    }

    gStyle->SetOptStat(0);

    TCanvas* cMMgg = new TCanvas("c_mm_vs_mgg", "MM vs mgg", 900, 700);
    cMMgg->SetRightMargin(0.14);
    hMMvsMgg->Draw("COLZ");
    cMMgg->SaveAs("acqu_geant_mm_vs_mgg.png");

    TCanvas* cEMM = new TCanvas("c_egamma_vs_mm", "Egamma vs MM", 900, 700);
    cEMM->SetRightMargin(0.14);
    hEgammaVsMM->Draw("COLZ");
    cEMM->SaveAs("acqu_geant_egamma_vs_mm.png");

    TCanvas* cETheta = new TCanvas("c_egamma_vs_theta", "Egamma vs theta", 900, 700);
    cETheta->SetRightMargin(0.14);
    hEgammaVsTheta->Draw("COLZ");
    cETheta->SaveAs("acqu_geant_egamma_vs_theta.png");

    TCanvas* cEThetaCM = new TCanvas("c_egamma_vs_theta_cm_before_after", "Egamma vs theta CM before/after", 1200, 600);
    cEThetaCM->Divide(2, 1);
    cEThetaCM->cd(1);
    gPad->SetRightMargin(0.14);
    hEgammaVsThetaCmBefore->Draw("COLZ");
    cEThetaCM->cd(2);
    gPad->SetRightMargin(0.14);
    hEgammaVsThetaCmAfter->Draw("COLZ");
    cEThetaCM->SaveAs("acqu_geant_egamma_vs_theta_cm_before_after.pdf");

    TCanvas* cOpeningPhi = new TCanvas("c_opening_angle_vs_phi_min", "opening angle vs Phi_min", 1200, 600);
    cOpeningPhi->Divide(2, 1);
    cOpeningPhi->cd(1);
    gPad->SetRightMargin(0.14);
    hOpeningPhiMinBefore->Draw("COLZ");
    cOpeningPhi->cd(2);
    gPad->SetRightMargin(0.14);
    hOpeningPhiMinAfter->Draw("COLZ");
    cOpeningPhi->SaveAs("acqu_geant_opening_angle_vs_phi_min_before_after.pdf");

    TCanvas* cOpeningDelta = new TCanvas("c_opening_minus_phi_min_vs_egamma", "opening minus Phi_min vs Egamma", 900, 700);
    cOpeningDelta->SetRightMargin(0.14);
    hOpeningMinusPhiMinAfter->Draw("COLZ");
    cOpeningDelta->SaveAs("acqu_geant_opening_minus_phi_min_vs_egamma_after_mgg_cut.pdf");

    TCanvas* cOpeningEgamma = new TCanvas("c_geant_phigg_vs_egamma", "GEANT Phi_gg vs Egamma", 1500, 500);
    cOpeningEgamma->Divide(3, 1);
    cOpeningEgamma->cd(1);
    gPad->SetRightMargin(0.14);
    hOpeningEgammaBefore->Draw("COLZ");
    cOpeningEgamma->cd(2);
    gPad->SetRightMargin(0.14);
    hOpeningEgammaAfter->Draw("COLZ");
    cOpeningEgamma->cd(3);
    gPad->SetRightMargin(0.14);
    hOpeningEgammaAfterMM->Draw("COLZ");
    cOpeningEgamma->SaveAs("acqu_geant_phigg_vs_egamma_before_after.pdf");

    TCanvas* cPhiMinMap = new TCanvas("c_phi_min_eq2_cm_map", "Phi_min Eq. 2 CM map", 1200, 750);
    for(int ie=0; ie<nEBins; ie++) {
        hPhiMinMap->GetXaxis()->SetBinLabel(ie + 1, Form("%.0f-%.0f", eBinLow[ie], eBinHigh[ie]));
        hRecoPaperBefore->GetXaxis()->SetBinLabel(ie + 1, Form("%.0f-%.0f", eBinLow[ie], eBinHigh[ie]));
        hRecoPaperAfter->GetXaxis()->SetBinLabel(ie + 1, Form("%.0f-%.0f", eBinLow[ie], eBinHigh[ie]));
        hRecoPaperAfterMM->GetXaxis()->SetBinLabel(ie + 1, Form("%.0f-%.0f", eBinLow[ie], eBinHigh[ie]));
    }
    hPhiMinMap->GetXaxis()->LabelsOption("v");
    hPhiMinMap->Draw("COLZ");
    cPhiMinMap->SaveAs("acqu_geant_phi_min_eq2_cm_map.pdf");

    TCanvas* cRecoPaper = new TCanvas("c_geant_reco_counts_paper_bins", "GEANT reco counts paper bins", 1200, 750);
    hRecoPaperAfterMM->GetXaxis()->LabelsOption("v");
    hRecoPaperAfterMM->Draw("COLZ");
    cRecoPaper->SaveAs("acqu_geant_reco_counts_paper_bins_after_mgg_mm_cut.pdf");

    TCanvas* cCountsBook = new TCanvas("c_geant_counts_grouped6_book", "GEANT counts grouped by six", 1500, 900);
    std::vector<TCanvas*> cCountsPages;
    cCountsBook->Print("acqu_geant_counts_paper_bins_grouped6.pdf[");
    for(int first=0; first<nEBins; first+=6) {
        int group = first / 6;
        TCanvas* cPage = new TCanvas(Form("c_geant_counts_paper_bins_grouped6_page%d", group + 1),
                                     Form("GEANT counts paper bins page %d", group + 1),
                                     1500, 900);
        cPage->Divide(3, 2);
        for(int ipad=0; ipad<6; ipad++) {
            int ie = first + ipad;
            if(ie >= nEBins) break;
            cPage->cd(ipad + 1);
            gPad->SetGrid();
            if(grGeantCounts[ie]->GetN() > 0) {
                grGeantCounts[ie]->Draw("AP");
                grGeantCounts[ie]->GetXaxis()->SetLimits(0.0, 160.0);
            } else {
                TH1F* frame = gPad->DrawFrame(0.0, 0.0, 160.0, 1.0);
                frame->SetTitle(Form("GEANT %.0f < E_{#gamma}^{lab} < %.0f MeV;#theta^{cm}_{#pi^{0}} (deg);Counts",
                                      eBinLow[ie], eBinHigh[ie]));
            }
        }
        cPage->Print("acqu_geant_counts_paper_bins_grouped6.pdf");
        cCountsPages.push_back(cPage);
    }
    cCountsBook->Print("acqu_geant_counts_paper_bins_grouped6.pdf]");

    TCanvas* cPhiBook = new TCanvas("c_phi_min_eq2_grouped6_book", "Phi min grouped by six", 1500, 900);
    std::vector<TCanvas*> cPhiPages;
    cPhiBook->Print("acqu_geant_phi_min_eq2_grouped6.pdf[");
    for(int first=0; first<nEBins; first+=6) {
        int group = first / 6;
        TCanvas* cPage = new TCanvas(Form("c_phi_min_eq2_grouped6_page%d", group + 1),
                                     Form("GEANT Phi min Eq. 2 page %d", group + 1),
                                     1500, 900);
        cPage->Divide(3, 2);
        for(int ipad=0; ipad<6; ipad++) {
            int ie = first + ipad;
            if(ie >= nEBins) break;
            cPage->cd(ipad + 1);
            gPad->SetGrid();
            grPhiMin[ie]->Draw("ALP");
        }
        cPage->Print("acqu_geant_phi_min_eq2_grouped6.pdf");
        cPhiPages.push_back(cPage);
    }
    cPhiBook->Print("acqu_geant_phi_min_eq2_grouped6.pdf]");

    TCanvas* c1D = new TCanvas("c_mgg_mm_1d", "mgg and MM", 1200, 600);
    c1D->Divide(2, 1);
    c1D->cd(1);
    hMgg->Draw("HIST");
    c1D->cd(2);
    hMM->Draw("HIST");
    c1D->SaveAs("acqu_geant_mgg_mm_1d.png");

    TFile* fout = new TFile("acqu_geant_requested_plots.root", "RECREATE");
    hMgg->Write();
    hMM->Write();
    hMMvsMgg->Write();
    hEgammaVsMM->Write();
    hEgammaVsTheta->Write();
    hEgammaVsThetaCmBefore->Write();
    hEgammaVsThetaCmAfter->Write();
    hOpeningPhiMinBefore->Write();
    hOpeningPhiMinAfter->Write();
    hOpeningMinusPhiMinAfter->Write();
    hOpeningEgammaBefore->Write();
    hOpeningEgammaAfter->Write();
    hOpeningEgammaAfterMM->Write();
    hRecoPaperBefore->Write();
    hRecoPaperAfter->Write();
    hRecoPaperAfterMM->Write();
    hPhiMinMap->Write();
    for(int ie=0; ie<nEBins; ie++) {
        grPhiMin[ie]->Write();
        grGeantCounts[ie]->Write();
    }
    cMMgg->Write();
    cEMM->Write();
    cETheta->Write();
    cEThetaCM->Write();
    cOpeningPhi->Write();
    cOpeningDelta->Write();
    cOpeningEgamma->Write();
    cPhiMinMap->Write();
    cRecoPaper->Write();
    cCountsBook->Write();
    for(int i=0; i<int(cCountsPages.size()); i++) cCountsPages[i]->Write();
    cPhiBook->Write();
    for(int i=0; i<int(cPhiPages.size()); i++) cPhiPages[i]->Write();
    c1D->Write();
    fout->Close();

    std::cout << "Input file: " << fname << std::endl;
    std::cout << "Processed events = " << nentries << std::endl;
    std::cout << "Events with tag = " << nWithTag << std::endl;
    std::cout << "Events with >=2 neutral clusters = " << nWithTwoNeutral << std::endl;
    std::cout << "Filled tag/event combinations = " << nFilled << std::endl;
    std::cout << "Cuts: clusterEnergy > " << eMin
              << " MeV, vetoEnergy < " << vetoMax << std::endl;
    std::cout << "Pi0 cuts: mgg [" << mggMin << "," << mggMax << "] MeV, MM ["
              << mmMin << "," << mmMax << "] MeV" << std::endl;
    std::cout << "Saved:" << std::endl;
    std::cout << "  acqu_geant_mm_vs_mgg.png" << std::endl;
    std::cout << "  acqu_geant_egamma_vs_mm.png" << std::endl;
    std::cout << "  acqu_geant_egamma_vs_theta.png" << std::endl;
    std::cout << "  acqu_geant_egamma_vs_theta_cm_before_after.pdf" << std::endl;
    std::cout << "  acqu_geant_opening_angle_vs_phi_min_before_after.pdf" << std::endl;
    std::cout << "  acqu_geant_opening_minus_phi_min_vs_egamma_after_mgg_cut.pdf" << std::endl;
    std::cout << "  acqu_geant_phigg_vs_egamma_before_after.pdf" << std::endl;
    std::cout << "  acqu_geant_phi_min_eq2_cm_map.pdf" << std::endl;
    std::cout << "  acqu_geant_phi_min_eq2_grouped6.pdf" << std::endl;
    std::cout << "  acqu_geant_reco_counts_paper_bins_after_mgg_mm_cut.pdf" << std::endl;
    std::cout << "  acqu_geant_counts_paper_bins_grouped6.pdf" << std::endl;
    std::cout << "  acqu_geant_mgg_mm_1d.png" << std::endl;
    std::cout << "  acqu_geant_requested_plots.root" << std::endl;

    f->Close();
}
