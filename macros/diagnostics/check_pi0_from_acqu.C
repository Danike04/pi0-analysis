/*
 * check_pi0_from_acqu.C
 *
 * Purpose
 * -------
 * Minimal Acqu reconstruction sanity check for pi0 -> gamma gamma. It loops
 * over neutral track pairs, applies the cluster-energy and veto selections,
 * fills the invariant mass of all accepted pairs, and separately fills the
 * pair closest to the nominal pi0 mass in each event.
 *
 * Input
 * -----
 * An Acqu ROOT file with a `tracks` tree. Default:
 * Acqu_CBTagg_31838.root.
 *
 * Default cuts
 * ------------
 * cluster energy >= 20 MeV
 * veto energy <= 1 MeV
 *
 * Output
 * ------
 * check_pi0_from_acqu.png
 * check_pi0_from_acqu.root
 *
 * Usage
 * -----
 * root -l -b -q 'check_pi0_from_acqu.C("Acqu_CBTagg_31838.root")'
 *
 * This is a reconstruction diagnostic only; no timing, tagger-energy or
 * cross-section normalization is applied.
 * Original analysis code follows unchanged; only this header was added.
 */

#include <algorithm>
#include <cmath>
#include <iostream>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TString.h"
#include "TTree.h"

void check_pi0_from_acqu(const char* fname="Acqu_CBTagg_31838.root",
                         Long64_t maxEvents=200000,
                         double eMin=20.0,
                         double vetoMax=1.0)
{
    TFile* f = TFile::Open(fname);
    if(!f || f->IsZombie()) {
        std::cout << "Cannot open " << fname << std::endl;
        return;
    }

    TTree* tracks = (TTree*)f->Get("tracks");
    if(!tracks) {
        std::cout << "Missing tree: tracks" << std::endl;
        f->Close();
        return;
    }

    const int maxTracks = 512;
    int nTracks = 0;
    double clusterEnergy[maxTracks];
    double theta[maxTracks];
    double phi[maxTracks];
    double vetoEnergy[maxTracks];
    int centralVeto[maxTracks];

    tracks->SetBranchAddress("nTracks", &nTracks);
    tracks->SetBranchAddress("clusterEnergy", clusterEnergy);
    tracks->SetBranchAddress("theta", theta);
    tracks->SetBranchAddress("phi", phi);
    tracks->SetBranchAddress("vetoEnergy", vetoEnergy);
    tracks->SetBranchAddress("centralVeto", centralVeto);

    TH1D* hAll = new TH1D("h_mgg_all_pairs",
                          "All neutral cluster pairs;m_{#gamma#gamma} (MeV);Counts",
                          400, 0.0, 400.0);
    TH1D* hBest = new TH1D("h_mgg_best_pair",
                           "Best pair per event, closest to m_{#pi^{0}};m_{#gamma#gamma} (MeV);Counts",
                           400, 0.0, 400.0);

    Long64_t nentries = tracks->GetEntries();
    if(maxEvents > 0 && maxEvents < nentries) nentries = maxEvents;

    const double deg = TMath::DegToRad();
    const double mpi0 = 134.9768;

    Long64_t usedEvents = 0;
    Long64_t usedPairs = 0;

    for(Long64_t ev=0; ev<nentries; ev++) {
        tracks->GetEntry(ev);

        double bestM = -1.0;
        double bestDist = 1.0e30;

        int nt = std::min(nTracks, maxTracks);
        for(int i=0; i<nt; i++) {
            if(clusterEnergy[i] < eMin) continue;
            if(vetoEnergy[i] > vetoMax) continue;

            for(int j=i+1; j<nt; j++) {
                if(clusterEnergy[j] < eMin) continue;
                if(vetoEnergy[j] > vetoMax) continue;

                double th1 = theta[i] * deg;
                double th2 = theta[j] * deg;
                double ph1 = phi[i] * deg;
                double ph2 = phi[j] * deg;

                double cosPsi = std::sin(th1)*std::sin(th2)*std::cos(ph1-ph2)
                              + std::cos(th1)*std::cos(th2);
                cosPsi = std::max(-1.0, std::min(1.0, cosPsi));

                double m2 = 2.0 * clusterEnergy[i] * clusterEnergy[j] * (1.0 - cosPsi);
                if(m2 <= 0.0) continue;

                double mgg = std::sqrt(m2);
                hAll->Fill(mgg);
                usedPairs++;

                double dist = std::fabs(mgg - mpi0);
                if(dist < bestDist) {
                    bestDist = dist;
                    bestM = mgg;
                }
            }
        }

        if(bestM > 0.0) {
            hBest->Fill(bestM);
            usedEvents++;
        }
    }

    TCanvas* c = new TCanvas("c_check_pi0_from_acqu", "pi0 check from Acqu", 1100, 500);
    c->Divide(2, 1);
    c->cd(1);
    hAll->Draw();
    c->cd(2);
    hBest->Draw();
    c->SaveAs("check_pi0_from_acqu.png");

    TFile* fout = new TFile("check_pi0_from_acqu.root", "RECREATE");
    hAll->Write();
    hBest->Write();
    c->Write();
    fout->Close();

    std::cout << "Input file: " << fname << std::endl;
    std::cout << "Processed events = " << nentries << std::endl;
    std::cout << "Events with at least one pair = " << usedEvents << std::endl;
    std::cout << "Accepted pairs = " << usedPairs << std::endl;
    std::cout << "Cuts: Ecluster > " << eMin << " MeV, vetoEnergy < " << vetoMax << " MeV" << std::endl;
    std::cout << "Saved: check_pi0_from_acqu.png, check_pi0_from_acqu.root" << std::endl;

    f->Close();
}
