/*
 * eff_det_from_acqu_tracks.C
 *
 * PURPOSE
 *   Early reconstructed-MC efficiency cross-check using only the Acqu `tracks`
 *   tree. It treats every tree entry as a generated event, selects neutral-like
 *   Crystal Ball clusters, forms all gamma-gamma pairs and counts events with
 *   a pi0 candidate in the invariant-mass window.
 *
 * DEFAULTS
 *   m_gg = 110-145 MeV; cluster threshold = 0 MeV.
 *
 * IMPORTANT LIMITATION
 *   N_generated is taken as tracks->GetEntries(). Therefore this macro is only
 *   meaningful for MC files in which the reconstructed tree entries have the
 *   intended one-to-one relation with generated events. The later paper-bin
 *   efficiency macros are the preferred production tools.
 *
 * OUTPUTS
 *   acqu_pi0_invariant_mass.pdf, acqu_nTracks.pdf,
 *   eff_det_from_acqu_tracks_output.root.
 *
 * USAGE EXAMPLE
 *   root -l -b -q 'eff_det_from_acqu_tracks.C("Acqu_Geant_4He_pi0.root")'
 *
 * REPOSITORY STATUS
 *   Efficiency cross-check/study; executable code unchanged.
 */
#include <iostream>
#include <cmath>
#include <vector>

void eff_det_from_acqu_tracks(const char* filename = "Acqu_Geant_4He_pi0.root")
{
    // ============================================================
    // Detection efficiency from Acqu reconstructed tracks
    // No MC truth information is used.
    // ============================================================

    TFile* f = TFile::Open(filename);

    if (!f || f->IsZombie()) {
        std::cout << "Error: cannot open file " << filename << std::endl;
        return;
    }

    TTree* t = (TTree*) f->Get("tracks");

    if (!t) {
        std::cout << "Error: tree 'tracks' not found." << std::endl;
        return;
    }

    const int MAXTRACKS = 100;

    Int_t nTracks;
    Double_t clusterEnergy[MAXTRACKS];
    Double_t theta[MAXTRACKS];
    Double_t phi[MAXTRACKS];
    Double_t vetoEnergy[MAXTRACKS];
    Int_t detectors[MAXTRACKS];
    Int_t clusterSize[MAXTRACKS];
    Int_t centralVeto[MAXTRACKS];

    t->SetBranchAddress("nTracks", &nTracks);
    t->SetBranchAddress("clusterEnergy", clusterEnergy);
    t->SetBranchAddress("theta", theta);
    t->SetBranchAddress("phi", phi);
    t->SetBranchAddress("vetoEnergy", vetoEnergy);
    t->SetBranchAddress("detectors", detectors);
    t->SetBranchAddress("clusterSize", clusterSize);
    t->SetBranchAddress("centralVeto", centralVeto);

    Long64_t Ngen = t->GetEntries();

    // Pi0 invariant mass cut
    const double mPi0L = 110.0;
    const double mPi0R = 145.0;

    // Optional cluster energy threshold.
    // Start with 0.0 if you want to reproduce the pure reconstructed histogram.
    // Try 20 or 50 MeV as systematic checks.
    const double E_thr = 0.0;

    TH1D* hMgg = new TH1D(
        "hMgg",
        "#pi^{0} invariant mass from Acqu tracks;m_{#gamma#gamma} [MeV];Counts",
        250,
        0,
        250
    );

    TH1D* hNTracks = new TH1D(
        "hNTracks",
        "Number of reconstructed tracks;nTracks;Events",
        20,
        0,
        20
    );

    Long64_t N_with_two_clusters = 0;
    Long64_t N_pi0 = 0;

    for (Long64_t i = 0; i < Ngen; i++) {

        t->GetEntry(i);

        hNTracks->Fill(nTracks);

        std::vector<int> good;

        for (int j = 0; j < nTracks; j++) {

            // detectors == 1 seems to be Crystal Ball from your Show(0).
            // centralVeto == -1 and vetoEnergy == 0 are neutral-like conditions.
            bool isCB = (detectors[j] == 1);
            bool isNeutral = (centralVeto[j] < 0 || vetoEnergy[j] == 0);
            bool aboveThr = (clusterEnergy[j] > E_thr);

            if (isCB && isNeutral && aboveThr) {
                good.push_back(j);
            }
        }

        if ((int)good.size() >= 2) {
            N_with_two_clusters++;
        }

        bool eventHasPi0 = false;

        // Loop over all possible photon pairs
        for (size_t a = 0; a < good.size(); a++) {
            for (size_t b = a + 1; b < good.size(); b++) {

                int i1 = good[a];
                int i2 = good[b];

                double E1 = clusterEnergy[i1]; // MeV
                double E2 = clusterEnergy[i2]; // MeV

                double th1 = theta[i1] * TMath::DegToRad();
                double th2 = theta[i2] * TMath::DegToRad();

                double ph1 = phi[i1] * TMath::DegToRad();
                double ph2 = phi[i2] * TMath::DegToRad();

                double cos12 =
                    std::sin(th1) * std::sin(th2) * std::cos(ph1 - ph2)
                    + std::cos(th1) * std::cos(th2);

                if (cos12 >  1.0) cos12 =  1.0;
                if (cos12 < -1.0) cos12 = -1.0;

                double mgg2 = 2.0 * E1 * E2 * (1.0 - cos12);

                if (mgg2 <= 0) continue;

                double mgg = std::sqrt(mgg2);

                hMgg->Fill(mgg);

                if (mgg > mPi0L && mgg < mPi0R) {
                    eventHasPi0 = true;
                }
            }
        }

        if (eventHasPi0) {
            N_pi0++;
        }
    }

    double eff_two_clusters = (double)N_with_two_clusters / (double)Ngen;
    double err_two_clusters =
        std::sqrt(eff_two_clusters * (1.0 - eff_two_clusters) / (double)Ngen);

    double eff_pi0 = (double)N_pi0 / (double)Ngen;
    double err_pi0 =
        std::sqrt(eff_pi0 * (1.0 - eff_pi0) / (double)Ngen);

    std::cout << "========================================" << std::endl;
    std::cout << "DETECTION EFFICIENCY FROM ACQU TRACKS" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Input file = " << filename << std::endl;
    std::cout << "Generated events = " << Ngen << std::endl;
    std::cout << "No MC truth information used." << std::endl;
    std::cout << "Cluster energy threshold = " << E_thr << " MeV" << std::endl;
    std::cout << "Pi0 mass cut = [" << mPi0L << ", " << mPi0R << "] MeV" << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Events with at least two neutral CB clusters = "
              << N_with_two_clusters << std::endl;
    std::cout << "Efficiency two clusters = "
              << eff_two_clusters * 100.0
              << " +/- "
              << err_two_clusters * 100.0
              << " %" << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Events with reconstructed pi0 = "
              << N_pi0 << std::endl;
    std::cout << "eps_det = "
              << eff_pi0
              << " = "
              << eff_pi0 * 100.0
              << " +/- "
              << err_pi0 * 100.0
              << " %" << std::endl;
    std::cout << "========================================" << std::endl;

    TCanvas* c1 = new TCanvas("c1", "Invariant mass", 900, 650);

    hMgg->SetLineColor(kBlack);
    hMgg->SetLineWidth(2);
    hMgg->Draw("E");

    double ymax = hMgg->GetMaximum();

    TLine* l1 = new TLine(mPi0L, 0, mPi0L, ymax);
    TLine* l2 = new TLine(mPi0R, 0, mPi0R, ymax);

    l1->SetLineColor(kRed);
    l2->SetLineColor(kRed);
    l1->SetLineWidth(2);
    l2->SetLineWidth(2);

    l1->Draw("same");
    l2->Draw("same");

    c1->SaveAs("acqu_pi0_invariant_mass.pdf");

    TCanvas* c2 = new TCanvas("c2", "Number of tracks", 900, 650);
    hNTracks->Draw("E");
    c2->SaveAs("acqu_nTracks.pdf");

    TFile* fout = new TFile("eff_det_from_acqu_tracks_output.root", "RECREATE");
    hMgg->Write();
    hNTracks->Write();
    fout->Close();

    std::cout << "Saved plots:" << std::endl;
    std::cout << "  acqu_pi0_invariant_mass.pdf" << std::endl;
    std::cout << "  acqu_nTracks.pdf" << std::endl;
    std::cout << "Saved histograms in eff_det_from_acqu_tracks_output.root" << std::endl;
}
