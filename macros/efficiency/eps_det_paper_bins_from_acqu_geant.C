/*
 * eps_det_paper_bins_from_acqu_geant.C
 *
 * PURPOSE
 *   Estimate the pi0 reconstruction/detection efficiency in the same
 *   (Egamma_lab, theta_pi0_cm) bins used by paper_style_diff_xs_from_acqu.C,
 *   starting from an AcquRoot reconstruction of a coherent Geant MC sample.
 *
 * RECONSTRUCTED NUMERATOR
 *   Events are reconstructed using neutral clusters, the best gamma-gamma pair
 *   (closest to the pi0 mass), the configurable m(gamma gamma) window, and
 *   optionally a missing-mass cut. applyMMCut is FALSE by default.
 *
 * GENERATED DENOMINATOR — IMPORTANT ASSUMPTION
 *   The macro counts generated/tagged events in each photon-energy bin and
 *   estimates the generated number in each theta_cm bin as
 *
 *       Ngen(E,theta) = Ngen(E) * DeltaOmega(theta) / (4*pi).
 *
 *   This is valid when the MC was generated isotropically in the CM frame.
 *   The macro does not read a separate generated theta truth distribution.
 *
 * EFFICIENCY
 *       eps_det(E,theta) = Nreco(E,theta) / Ngen(E,theta)
 *
 * OUTPUTS
 *   eps_det_paper_bins_from_acqu_geant.txt
 *   eps_det_paper_bins_from_acqu_geant.root
 *   eps_det_paper_bins_from_acqu_geant.pdf
 *   eps_det_paper_bins_reco_counts.pdf
 *   eps_det_paper_bins_checks.pdf
 *
 * The text output format is directly readable by
 * paper_style_diff_xs_from_acqu.C when paper-bin efficiency mode is enabled.
 *
 * NOTE FOR THE REPOSITORY
 *   The numerical method and defaults below are preserved from the original.
 */
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TTree.h"

int find_paper_bin(double x, const std::vector<double>& low, const std::vector<double>& high)
{
    for(size_t i=0; i<low.size(); i++) {
        if(x >= low[i] && x <= high[i]) return int(i);
    }
    return -1;
}

double paper_solid_angle_bin(double thetaMinDeg, double thetaMaxDeg)
{
    double th1 = thetaMinDeg * TMath::DegToRad();
    double th2 = thetaMaxDeg * TMath::DegToRad();
    return 2.0 * TMath::Pi() * (std::cos(th1) - std::cos(th2));
}

void eps_det_paper_bins_from_acqu_geant(const char* fname="Acqu_geant_He4pi0.root",
                                        Long64_t maxEvents=-1,
                                        double eMin=20.0,
                                        double vetoMax=1.0,
                                        double mggMin=110.0,
                                        double mggMax=155.0,
                                        bool applyMMCut=false,
                                        double mmMin=-10.0,
                                        double mmMax=40.0,
                                        double thetaCmMinDeg=5.0,
                                        double thetaCmMaxDeg=150.0,
                                        double thetaCmBinWidthDeg=5.0,
                                        Long64_t printEvery=50000)
{
    const int nChMax = 352;
    const int maxTracks = 512;
    const int maxTagged = 2048;

    const double beamEMax = 855.0;
    const double targetMass = 3727.38;
    const double mpi0 = 134.9768;

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
    int nThetaBins = int(std::floor((thetaCmMaxDeg - thetaCmMinDeg) / thetaCmBinWidthDeg + 0.5));
    if(nThetaBins < 1) nThetaBins = 1;

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

    int nCh = std::min(nTagger, nChMax);

    std::vector<double> nGenE(nEBins, 0.0);
    std::vector< std::vector<double> > nGen(nEBins, std::vector<double>(nThetaBins, 0.0));
    std::vector< std::vector<double> > nDet(nEBins, std::vector<double>(nThetaBins, 0.0));
    std::vector< std::vector<double> > nReco(nEBins, std::vector<double>(nThetaBins, 0.0));

    TH1D* hMgg = new TH1D("h_mgg_best_pair_mc",
                          "MC best #gamma#gamma pair;m_{#gamma#gamma} (MeV);Counts",
                          400, 0.0, 400.0);
    TH1D* hMM = new TH1D("h_mm_best_pair_mc",
                         "MC missing mass after m_{#gamma#gamma} cut;MM-M_{^{4}He} (MeV);Counts",
                         250, -100.0, 150.0);
    TH2D* hRecoMap = new TH2D("h_reco_counts_paper_bins",
                              "Reco counts in paper bins;E_{#gamma}^{lab} bin;#theta_{#pi^{0}}^{cm} (deg)",
                              nEBins, 0.0, double(nEBins),
                              nThetaBins, thetaCmMinDeg, thetaCmMinDeg + nThetaBins * thetaCmBinWidthDeg);
    TH2D* hEps = new TH2D("h_eps_det_paper_bins",
                          "#epsilon_{det}(E_{#gamma}^{lab},#theta^{cm});E_{#gamma}^{lab} bin;#theta_{#pi^{0}}^{cm} (deg)",
                          nEBins, 0.0, double(nEBins),
                          nThetaBins, thetaCmMinDeg, thetaCmMinDeg + nThetaBins * thetaCmBinWidthDeg);

    for(int ie=0; ie<nEBins; ie++) {
        hRecoMap->GetXaxis()->SetBinLabel(ie + 1, Form("%.0f-%.0f", eBinLow[ie], eBinHigh[ie]));
        hEps->GetXaxis()->SetBinLabel(ie + 1, Form("%.0f-%.0f", eBinLow[ie], eBinHigh[ie]));
    }

    Long64_t nentries = std::min(tracks->GetEntries(), tagger->GetEntries());
    if(maxEvents > 0 && maxEvents < nentries) nentries = maxEvents;

    const double deg = TMath::DegToRad();
    TLorentzVector target(0.0, 0.0, 0.0, targetMass);

    Long64_t nEventsWithTag = 0;
    Long64_t nEventsWithTwoNeutral = 0;
    Long64_t nEventsReco = 0;

    for(Long64_t ev=0; ev<nentries; ev++) {
        if(printEvery > 0 && (ev == 0 || ev % printEvery == 0)) {
            std::cout << "event " << ev << "/" << nentries
                      << " (" << 100.0 * double(ev) / double(nentries) << "%)"
                      << " tag=" << nEventsWithTag
                      << " twoNeutral=" << nEventsWithTwoNeutral
                      << " reco=" << nEventsReco
                      << std::endl;
        }

        tracks->GetEntry(ev);
        tagger->GetEntry(ev);

        bool seenCh[nChMax] = {false};
        int eventChannels[nChMax];
        int nEventChannels = 0;

        int ntag = std::min(nTagged, maxTagged);
        for(int it=0; it<ntag; it++) {
            int ch = taggedChannel[it];
            if(ch < 0 || ch >= nCh) continue;
            if(seenCh[ch]) continue;
            seenCh[ch] = true;
            eventChannels[nEventChannels] = ch;
            nEventChannels++;

            int ie = find_paper_bin(taggerPhotonEnergy[ch], eBinLow, eBinHigh);
            if(ie >= 0) nGenE[ie] += 1.0;
        }

        if(nEventChannels == 0) continue;
        nEventsWithTag++;

        int neutralIdx[maxTracks];
        int nNeutral = 0;
        int nt = std::min(nTracks, maxTracks);
        for(int i=0; i<nt; i++) {
            if(clusterEnergy[i] <= eMin) continue;
            if(vetoEnergy[i] >= vetoMax) continue;
            neutralIdx[nNeutral] = i;
            nNeutral++;
        }

        if(nNeutral < 2) continue;
        nEventsWithTwoNeutral++;

        double bestM = -1.0;
        double bestDist = 1.0e30;
        TLorentzVector bestPi0;

        for(int a=0; a<nNeutral; a++) {
            int i = neutralIdx[a];
            for(int b=a+1; b<nNeutral; b++) {
                int j = neutralIdx[b];

                double th1 = theta[i] * deg;
                double th2 = theta[j] * deg;
                double ph1 = phi[i] * deg;
                double ph2 = phi[j] * deg;

                TLorentzVector g1;
                TLorentzVector g2;
                g1.SetPxPyPzE(clusterEnergy[i]*std::sin(th1)*std::cos(ph1),
                              clusterEnergy[i]*std::sin(th1)*std::sin(ph1),
                              clusterEnergy[i]*std::cos(th1),
                              clusterEnergy[i]);
                g2.SetPxPyPzE(clusterEnergy[j]*std::sin(th2)*std::cos(ph2),
                              clusterEnergy[j]*std::sin(th2)*std::sin(ph2),
                              clusterEnergy[j]*std::cos(th2),
                              clusterEnergy[j]);

                TLorentzVector pi0 = g1 + g2;
                double mgg = pi0.M();
                if(mgg <= 0.0) continue;

                double dist = std::fabs(mgg - mpi0);
                if(dist < bestDist) {
                    bestDist = dist;
                    bestM = mgg;
                    bestPi0 = pi0;
                }
            }
        }

        if(bestM <= 0.0) continue;
        hMgg->Fill(bestM);

        bool passMgg = (bestM >= mggMin && bestM <= mggMax);
        if(!passMgg) continue;
        nEventsReco++;

        for(int k=0; k<nEventChannels; k++) {
            int ch = eventChannels[k];
            double egamma = taggerPhotonEnergy[ch];
            int ie = find_paper_bin(egamma, eBinLow, eBinHigh);
            if(ie < 0) continue;

            TLorentzVector beam(0.0, 0.0, egamma, egamma);
            TLorentzVector recoil = beam + target - bestPi0;
            double mm = recoil.M() - targetMass;
            hMM->Fill(mm);

            if(applyMMCut && (mm < mmMin || mm > mmMax)) continue;

            TLorentzVector pi0CM = bestPi0;
            TLorentzVector totalCM = beam + target;
            pi0CM.Boost(-totalCM.BoostVector());
            double thetaCmDeg = pi0CM.Theta() * TMath::RadToDeg();

            if(thetaCmDeg < thetaCmMinDeg || thetaCmDeg >= thetaCmMaxDeg) continue;
            int ith = int((thetaCmDeg - thetaCmMinDeg) / thetaCmBinWidthDeg);
            if(ith < 0 || ith >= nThetaBins) continue;

            nReco[ie][ith] += 1.0;
            nDet[ie][ith] += 1.0;
            hRecoMap->Fill(double(ie) + 0.5, thetaCmDeg);
        }
    }

    for(int ie=0; ie<nEBins; ie++) {
        for(int ith=0; ith<nThetaBins; ith++) {
            double th1 = thetaCmMinDeg + ith * thetaCmBinWidthDeg;
            double th2 = th1 + thetaCmBinWidthDeg;
            double fracOmega = paper_solid_angle_bin(th1, th2) / (4.0 * TMath::Pi());
            nGen[ie][ith] = nGenE[ie] * fracOmega;

            if(nGen[ie][ith] <= 0.0) continue;

            double eps = nReco[ie][ith] / nGen[ie][ith];
            double deps = 0.0;
            if(eps >= 0.0 && eps <= 1.0) {
                deps = std::sqrt(eps * (1.0 - eps) / nGen[ie][ith]);
            } else {
                deps = std::sqrt(nReco[ie][ith]) / nGen[ie][ith];
            }
            hEps->SetBinContent(ie + 1, ith + 1, eps);
            hEps->SetBinError(ie + 1, ith + 1, deps);
        }
    }

    std::ofstream out("eps_det_paper_bins_from_acqu_geant.txt");
    out << "# input: " << fname << "\n";
    out << "# definition: eps_det(E_lab,theta_cm) = Nreco / Ngen_bin\n";
    out << "# denominator: Ngen_bin = Ntagged_in_Ebin * DeltaOmega_cm/(4*pi), valid if MC was generated isotropically in CM\n";
    out << "# reco selection: two neutral clusters E>" << eMin << ", veto<" << vetoMax
        << ", best pair " << mggMin << " < mgg < " << mggMax << " MeV\n";
    out << "# applyMMCut = " << applyMMCut << ", MM cut = [" << mmMin << "," << mmMax << "] MeV\n";
    out << "# E_low E_high E_center thetaCm_low thetaCm_high thetaCm_center DeltaOmega_sr Ngen_E Ngen_bin Nreco eps_det deps_det\n";

    for(int ie=0; ie<nEBins; ie++) {
        for(int ith=0; ith<nThetaBins; ith++) {
            double th1 = thetaCmMinDeg + ith * thetaCmBinWidthDeg;
            double th2 = th1 + thetaCmBinWidthDeg;
            double thCenter = 0.5 * (th1 + th2);
            double deltaOmega = paper_solid_angle_bin(th1, th2);
            double eps = 0.0;
            double deps = 0.0;
            if(nGen[ie][ith] > 0.0) {
                eps = nReco[ie][ith] / nGen[ie][ith];
                if(eps >= 0.0 && eps <= 1.0) {
                    deps = std::sqrt(eps * (1.0 - eps) / nGen[ie][ith]);
                } else {
                    deps = std::sqrt(nReco[ie][ith]) / nGen[ie][ith];
                }
            }

            out << eBinLow[ie] << " "
                << eBinHigh[ie] << " "
                << 0.5 * (eBinLow[ie] + eBinHigh[ie]) << " "
                << th1 << " "
                << th2 << " "
                << thCenter << " "
                << deltaOmega << " "
                << nGenE[ie] << " "
                << nGen[ie][ith] << " "
                << nReco[ie][ith] << " "
                << eps << " "
                << deps << "\n";
        }
    }
    out.close();

    TCanvas* c1 = new TCanvas("c_eps_det_paper_bins", "eps det paper bins", 1200, 750);
    hEps->GetXaxis()->LabelsOption("v");
    hEps->Draw("COLZ");
    c1->SaveAs("eps_det_paper_bins_from_acqu_geant.pdf");

    TCanvas* c2 = new TCanvas("c_eps_det_paper_reco_counts", "reco counts paper bins", 1200, 750);
    hRecoMap->GetXaxis()->LabelsOption("v");
    hRecoMap->Draw("COLZ");
    c2->SaveAs("eps_det_paper_bins_reco_counts.pdf");

    TCanvas* c3 = new TCanvas("c_eps_det_paper_checks", "checks", 1100, 500);
    c3->Divide(2, 1);
    c3->cd(1);
    hMgg->Draw();
    c3->cd(2);
    hMM->Draw();
    c3->SaveAs("eps_det_paper_bins_checks.pdf");

    TFile* fout = new TFile("eps_det_paper_bins_from_acqu_geant.root", "RECREATE");
    hEps->Write();
    hRecoMap->Write();
    hMgg->Write();
    hMM->Write();
    c1->Write();
    c2->Write();
    c3->Write();
    fout->Close();

    std::cout << "Processed events = " << nentries << std::endl;
    std::cout << "Events with tag = " << nEventsWithTag << std::endl;
    std::cout << "Events with >=2 neutral = " << nEventsWithTwoNeutral << std::endl;
    std::cout << "Events passing reco selection = " << nEventsReco << std::endl;
    std::cout << "Saved:" << std::endl;
    std::cout << "  eps_det_paper_bins_from_acqu_geant.txt" << std::endl;
    std::cout << "  eps_det_paper_bins_from_acqu_geant.pdf" << std::endl;
    std::cout << "  eps_det_paper_bins_reco_counts.pdf" << std::endl;
    std::cout << "  eps_det_paper_bins_checks.pdf" << std::endl;
    std::cout << "  eps_det_paper_bins_from_acqu_geant.root" << std::endl;

    f->Close();
}
