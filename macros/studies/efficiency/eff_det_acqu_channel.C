/*
 * eff_det_acqu_channel.C
 *
 * PURPOSE
 *   Early Acqu-reconstructed MC efficiency study restricted to selected MC
 *   tagger channels. Reconstructs neutral gamma-gamma pairs and measures the
 *   fraction passing the pi0 invariant-mass selection.
 *
 * HISTORICAL ASSUMPTIONS
 *   - N_generated is hard-coded to 100000.
 *   - MC channels 259 and 260 are hard-coded and both are enabled by default.
 *   - Default m_gg window is 110-145 MeV; cluster threshold is 0 MeV.
 *
 * Because of these dataset-specific assumptions this is NOT the production
 * detector-efficiency extractor. It is retained as a useful development
 * cross-check of the reconstruction.
 *
 * OUTPUTS
 *   Acqu invariant-mass/tagger/multiplicity PDFs and
 *   eff_det_acqu_channel_output.root.
 *
 * USAGE EXAMPLE
 *   root -l -b -q 'eff_det_acqu_channel.C("Acqu_Geant_4He_pi0.root")'
 *
 * REPOSITORY STATUS
 *   Historical efficiency study; executable code unchanged.
 */
#include <iostream>
#include <cmath>
#include <vector>

void eff_det_acqu_channel(const char* filename = "Acqu_Geant_4He_pi0.root")
{
    // ============================================================
    // Detection efficiency from Acqu reconstructed tracks
    // Treats MC as real data:
    // - uses reconstructed clusters
    // - uses tagger channel
    // - does NOT use MC truth
    // ============================================================

    TFile* f = TFile::Open(filename);

    if (!f || f->IsZombie()) {
        std::cout << "Error: cannot open file " << filename << std::endl;
        return;
    }

    TTree* tracks = (TTree*) f->Get("tracks");
    TTree* tagger = (TTree*) f->Get("tagger");

    if (!tracks) {
        std::cout << "Error: tree 'tracks' not found." << std::endl;
        return;
    }

    if (!tagger) {
        std::cout << "Error: tree 'tagger' not found." << std::endl;
        return;
    }

    Long64_t Ntracks = tracks->GetEntries();
    Long64_t Ntagger = tagger->GetEntries();

    if (Ntracks != Ntagger) {
        std::cout << "Warning: tracks and tagger have different entries!" << std::endl;
        std::cout << "tracks entries = " << Ntracks << std::endl;
        std::cout << "tagger entries = " << Ntagger << std::endl;
    }

    const int MAXTRACKS = 100;
    const int MAXTAGGED = 100;

    // Tracks branches
    Int_t nTracks;
    Double_t clusterEnergy[MAXTRACKS];
    Double_t theta[MAXTRACKS];
    Double_t phi[MAXTRACKS];
    Double_t vetoEnergy[MAXTRACKS];
    Int_t detectors[MAXTRACKS];
    Int_t clusterSize[MAXTRACKS];
    Int_t centralVeto[MAXTRACKS];

    tracks->SetBranchAddress("nTracks", &nTracks);
    tracks->SetBranchAddress("clusterEnergy", clusterEnergy);
    tracks->SetBranchAddress("theta", theta);
    tracks->SetBranchAddress("phi", phi);
    tracks->SetBranchAddress("vetoEnergy", vetoEnergy);
    tracks->SetBranchAddress("detectors", detectors);
    tracks->SetBranchAddress("clusterSize", clusterSize);
    tracks->SetBranchAddress("centralVeto", centralVeto);

    // Tagger branches
    Int_t nTagged;
    Int_t taggedChannel[MAXTAGGED];
    Double_t taggedTime[MAXTAGGED];

    tagger->SetBranchAddress("nTagged", &nTagged);
    tagger->SetBranchAddress("taggedChannel", taggedChannel);
    tagger->SetBranchAddress("taggedTime", taggedTime);

    // ============================================================
    // Analysis settings
    // ============================================================

    const double Ngen = 100000.0;

    // MC tagger channel.
    // From Show(0) you saw taggedChannel = 259.
    // Use 259, 260, or both depending on what you want.
    const int CH_MC_1 = 259;
    const int CH_MC_2 = 260;

    // Set this to true if you want to accept both 259 and 260
    const bool useTwoChannels = true;

    // Pi0 invariant mass cut
    const double mPi0L = 110.0;
    const double mPi0R = 145.0;

    // Optional cluster energy threshold in MeV.
    // Start from 0 MeV, then try 20 or 50 MeV as checks.
    const double E_thr = 0.0;

    // ============================================================
    // Histograms
    // ============================================================

    TH1D* hMgg_all = new TH1D(
        "hMgg_all",
        "m_{#gamma#gamma}, all events;m_{#gamma#gamma} [MeV];Counts",
        250, 0, 250
    );

    TH1D* hMgg_ch = new TH1D(
        "hMgg_ch",
        "m_{#gamma#gamma}, selected tagger channel;m_{#gamma#gamma} [MeV];Counts",
        250, 0, 250
    );

    TH1D* hTagg = new TH1D(
        "hTagg",
        "Tagged channel;tagged channel;Events",
        400, 0, 400
    );

    TH1D* hNTracks = new TH1D(
        "hNTracks",
        "Number of reconstructed tracks;nTracks;Events",
        20, 0, 20
    );

    // ============================================================
    // Counters
    // ============================================================

    Long64_t N_total_entries = tracks->GetEntries();
    Long64_t N_tagged_selected = 0;
    Long64_t N_two_clusters_selected = 0;
    Long64_t N_pi0_selected = 0;

    Long64_t N_pi0_all_channels = 0;

    // ============================================================
    // Event loop
    // ============================================================

    for (Long64_t ev = 0; ev < N_total_entries; ev++) {

        tracks->GetEntry(ev);
        tagger->GetEntry(ev);

        hNTracks->Fill(nTracks);

        bool channelSelected = false;

        for (int it = 0; it < nTagged; it++) {

            hTagg->Fill(taggedChannel[it]);

            if (useTwoChannels) {
                if (taggedChannel[it] == CH_MC_1 || taggedChannel[it] == CH_MC_2) {
                    channelSelected = true;
                }
            } else {
                if (taggedChannel[it] == CH_MC_2) {
                    channelSelected = true;
                }
            }
        }

        if (channelSelected) {
            N_tagged_selected++;
        }

        std::vector<int> good;

        for (int j = 0; j < nTracks; j++) {

            bool isCB = (detectors[j] == 1);
            bool isNeutral = (centralVeto[j] < 0 || vetoEnergy[j] == 0);
            bool aboveThr = (clusterEnergy[j] > E_thr);

            if (isCB && isNeutral && aboveThr) {
                good.push_back(j);
            }
        }

        if (channelSelected && (int)good.size() >= 2) {
            N_two_clusters_selected++;
        }

        bool hasPi0_all = false;
        bool hasPi0_ch = false;

        for (size_t a = 0; a < good.size(); a++) {
            for (size_t b = a + 1; b < good.size(); b++) {

                int i1 = good[a];
                int i2 = good[b];

                double E1 = clusterEnergy[i1];
                double E2 = clusterEnergy[i2];

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

                hMgg_all->Fill(mgg);

                if (channelSelected) {
                    hMgg_ch->Fill(mgg);
                }

                if (mgg > mPi0L && mgg < mPi0R) {
                    hasPi0_all = true;

                    if (channelSelected) {
                        hasPi0_ch = true;
                    }
                }
            }
        }

        if (hasPi0_all) {
            N_pi0_all_channels++;
        }

        if (hasPi0_ch) {
            N_pi0_selected++;
        }
    }

    // ============================================================
    // Efficiencies
    // ============================================================

    double eff_pi0_global = (double) N_pi0_all_channels / Ngen;
    double err_pi0_global = std::sqrt(eff_pi0_global * (1.0 - eff_pi0_global) / Ngen);

    double eff_pi0_ch = (double) N_pi0_selected / Ngen;
    double err_pi0_ch = std::sqrt(eff_pi0_ch * (1.0 - eff_pi0_ch) / Ngen);

    double eff_two_clusters_ch = (double) N_two_clusters_selected / Ngen;
    double err_two_clusters_ch =
        std::sqrt(eff_two_clusters_ch * (1.0 - eff_two_clusters_ch) / Ngen);

    // Efficiency conditional on selected tagger events.
    // This is NOT the usual absolute detection efficiency, but it is useful as a check.
    double eff_conditional = 0.0;
    double err_conditional = 0.0;

    if (N_tagged_selected > 0) {
        eff_conditional = (double) N_pi0_selected / (double) N_tagged_selected;
        err_conditional =
            std::sqrt(eff_conditional * (1.0 - eff_conditional) / (double) N_tagged_selected);
    }

    // ============================================================
    // Output
    // ============================================================

    std::cout << "========================================" << std::endl;
    std::cout << "DETECTION EFFICIENCY FROM ACQU TREE" << std::endl;
    std::cout << "========================================" << std::endl;

    std::cout << "Input file = " << filename << std::endl;
    std::cout << "Total entries = " << N_total_entries << std::endl;
    std::cout << "Generated events used for denominator = " << Ngen << std::endl;
    std::cout << "No MC truth information used." << std::endl;

    std::cout << "----------------------------------------" << std::endl;

    if (useTwoChannels) {
        std::cout << "Selected tagger channels = "
                  << CH_MC_1 << " and " << CH_MC_2 << std::endl;
    } else {
        std::cout << "Selected tagger channel = "
                  << CH_MC_2 << std::endl;
    }

    std::cout << "Pi0 mass cut = ["
              << mPi0L << ", " << mPi0R << "] MeV" << std::endl;

    std::cout << "Cluster energy threshold = "
              << E_thr << " MeV" << std::endl;

    std::cout << "----------------------------------------" << std::endl;

    std::cout << "Tagged events in selected channel(s) = "
              << N_tagged_selected << std::endl;

    std::cout << "Events with at least two neutral CB clusters in selected channel(s) = "
              << N_two_clusters_selected << std::endl;

    std::cout << "Events with reconstructed pi0 in selected channel(s) = "
              << N_pi0_selected << std::endl;

    std::cout << "----------------------------------------" << std::endl;

    std::cout << "Global pi0 efficiency, all channels:" << std::endl;
    std::cout << "Npi0_all = " << N_pi0_all_channels << std::endl;
    std::cout << "eps_pi0_global = "
              << eff_pi0_global * 100.0
              << " +/- "
              << err_pi0_global * 100.0
              << " %" << std::endl;

    std::cout << "----------------------------------------" << std::endl;

    std::cout << "Absolute pi0 efficiency for selected channel(s):" << std::endl;
    std::cout << "Npi0_selected = " << N_pi0_selected << std::endl;
    std::cout << "eps_det = "
              << eff_pi0_ch
              << " = "
              << eff_pi0_ch * 100.0
              << " +/- "
              << err_pi0_ch * 100.0
              << " %" << std::endl;

    std::cout << "----------------------------------------" << std::endl;

    std::cout << "Conditional efficiency within selected tagged events:" << std::endl;
    std::cout << "Npi0_selected / Ntagged_selected = "
              << eff_conditional * 100.0
              << " +/- "
              << err_conditional * 100.0
              << " %" << std::endl;

    std::cout << "========================================" << std::endl;

    // ============================================================
    // Plots
    // ============================================================

    gStyle->SetOptStat(0);

    TCanvas* c1 = new TCanvas("c1", "Invariant mass all", 900, 650);
    hMgg_all->SetLineColor(kBlack);
    hMgg_all->SetLineWidth(2);
    hMgg_all->Draw("E");

    double ymax1 = hMgg_all->GetMaximum();
    TLine* l1a = new TLine(mPi0L, 0, mPi0L, ymax1);
    TLine* l1b = new TLine(mPi0R, 0, mPi0R, ymax1);
    l1a->SetLineColor(kRed);
    l1b->SetLineColor(kRed);
    l1a->SetLineWidth(2);
    l1b->SetLineWidth(2);
    l1a->Draw("same");
    l1b->Draw("same");

    c1->SaveAs("acqu_mgg_all_channels.pdf");

    TCanvas* c2 = new TCanvas("c2", "Invariant mass selected channel", 900, 650);
    hMgg_ch->SetLineColor(kBlack);
    hMgg_ch->SetLineWidth(2);
    hMgg_ch->Draw("E");

    double ymax2 = hMgg_ch->GetMaximum();
    TLine* l2a = new TLine(mPi0L, 0, mPi0L, ymax2);
    TLine* l2b = new TLine(mPi0R, 0, mPi0R, ymax2);
    l2a->SetLineColor(kRed);
    l2b->SetLineColor(kRed);
    l2a->SetLineWidth(2);
    l2b->SetLineWidth(2);
    l2a->Draw("same");
    l2b->Draw("same");

    c2->SaveAs("acqu_mgg_selected_channel.pdf");

    TCanvas* c3 = new TCanvas("c3", "Tagger channels", 900, 650);
    hTagg->SetLineColor(kBlack);
    hTagg->SetLineWidth(2);
    hTagg->Draw("E");
    c3->SaveAs("acqu_tagger_channels.pdf");

    TCanvas* c4 = new TCanvas("c4", "nTracks", 900, 650);
    hNTracks->SetLineColor(kBlack);
    hNTracks->SetLineWidth(2);
    hNTracks->Draw("E");
    c4->SaveAs("acqu_nTracks.pdf");

    TFile* fout = new TFile("eff_det_acqu_channel_output.root", "RECREATE");
    hMgg_all->Write();
    hMgg_ch->Write();
    hTagg->Write();
    hNTracks->Write();
    fout->Close();

    std::cout << "Saved plots:" << std::endl;
    std::cout << "  acqu_mgg_all_channels.pdf" << std::endl;
    std::cout << "  acqu_mgg_selected_channel.pdf" << std::endl;
    std::cout << "  acqu_tagger_channels.pdf" << std::endl;
    std::cout << "  acqu_nTracks.pdf" << std::endl;
}
