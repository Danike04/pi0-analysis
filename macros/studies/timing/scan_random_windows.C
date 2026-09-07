/*
 * scan_random_windows.C
 *
 * PURPOSE
 *   Campaign-specific timing-sideband study for the GoAT reconstruction. The
 *   macro keeps one prompt window fixed and scans six alternative random
 *   timing windows to check how strongly the prompt-random-subtracted missing-
 *   energy yield depends on the random sideband choice.
 *
 * WHAT THE MACRO ACTUALLY DOES
 *   - Opens GoAT_CBTagg_31838.root by default.
 *   - Uses tagger channel 51 and Egamma = 212 MeV.
 *   - Requires 110 < m(pi0) < 145 MeV.
 *   - Uses a fixed prompt window 525-610 ns.
 *   - Scans random windows from 700-850 ns through 825-975 ns.
 *   - Reconstructs missing energy and compares prompt, scaled-random and
 *     prompt-random yields in the 17.5-67.5 MeV integration region.
 *   - Also prints a diagnostic cross-section normalization estimate using
 *     hard-coded eps_tag, eps_det and target thickness values.
 *
 * OUTPUT
 *   Diagnostic canvases and numerical information printed to stdout. The
 *   original macro does not write a ROOT/PDF summary file.
 *
 * USAGE
 *   root -l 'scan_random_windows.C'
 *   scan_random_windows();
 *
 * REPOSITORY STATUS
 *   Preserved as a historical timing/systematic study, not as a production
 *   selection macro. The executable code below is unchanged from the raw file.
 */
void scan_random_windows() {

    // ============================================================
    // RANDOM-WINDOW SCAN
    // Checks whether the chosen random timing sideband over-subtracts
    // ============================================================

    const char* filename = "GoAT_CBTagg_31838.root";

    const int CH = 51;
    const double Egamma = 212.0;       // MeV
    const double Mtarget = 3727.38;    // MeV, 4He target mass

    // pi0 invariant mass cut
    const double mPi0L = 110.0;
    const double mPi0R = 145.0;

    // Missing-energy signal window
    const double meL = 17.5;
    const double meR = 67.5;

    // Normalization constants
    const double eps_tag = 0.3;
    const double eps_det = 0.70;
    const double thickness = 0.940e-7; // microbarn^-1

    // Fixed prompt window
    const double promptL = 525.0;
    const double promptR = 610.0;

    // Random windows to scan
    const int NTEST = 6;

    double randomL_arr[NTEST] = {
       700.0, 725.0, 750.0, 775.0, 800.0, 825.0
    };

    double randomR_arr[NTEST] = {
        850.0, 875.0, 900.0, 925.0, 950.0, 975.0
    };

    // ============================================================
    // Open file and trees
    // ============================================================

    TFile *f = new TFile(filename);

    if (!f || f->IsZombie()) {
        std::cout << "Error: cannot open file " << filename << std::endl;
        return;
    }

    TTree *tagger = (TTree*)f->Get("tagger");
    TTree *pi0    = (TTree*)f->Get("neutralPions");
    TTree *pair   = (TTree*)f->Get("pairSpec");

    if (!tagger || !pi0) {
        std::cout << "Error: tagger or neutralPions tree not found." << std::endl;
        return;
    }

    // ============================================================
    // Branches
    // ============================================================

    const int MAXTAG = 5000;
    const int MAXPI0 = 200;

    Int_t nTagged;
    Int_t taggedChannel[MAXTAG];
    Double_t taggedTime[MAXTAG];

    Int_t nParticles;
    Double_t clusterEnergy[MAXPI0];
    Double_t theta[MAXPI0];
    Double_t phi[MAXPI0];
    Double_t mass[MAXPI0];
    Double_t time[MAXPI0];

    tagger->SetBranchAddress("nTagged", &nTagged);
    tagger->SetBranchAddress("taggedChannel", taggedChannel);
    tagger->SetBranchAddress("taggedTime", taggedTime);

    pi0->SetBranchAddress("nParticles", &nParticles);
    pi0->SetBranchAddress("clusterEnergy", clusterEnergy);
    pi0->SetBranchAddress("theta", theta);
    pi0->SetBranchAddress("phi", phi);
    pi0->SetBranchAddress("mass", mass);
    pi0->SetBranchAddress("time", time);

    // ============================================================
    // Electron flux Ne
    // ============================================================

    double Ne = 0.0;

    if (pair) {

        Double_t scalerOpen[328];
        pair->SetBranchAddress("scalerOpen", scalerOpen);

        Long64_t nPair = pair->GetEntries();

        for (Long64_t i = 0; i < nPair; i++) {
            pair->GetEntry(i);
            Ne += scalerOpen[CH];
        }
    }

    double denom = eps_tag * Ne * thickness * eps_det;
    double Y_for_150 = 150.0 * denom;

    std::cout << "========================================" << std::endl;
    std::cout << "RANDOM-WINDOW DIAGNOSTIC SETTINGS" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Input file = " << filename << std::endl;
    std::cout << "Tagger channel = " << CH << std::endl;
    std::cout << "Photon energy = " << Egamma << " MeV" << std::endl;
    std::cout << "Prompt window = " << promptL << " to " << promptR << std::endl;
    std::cout << "Missing-energy window = "
              << meL << " to " << meR << " MeV" << std::endl;
    std::cout << "Electron flux Ne = " << Ne << std::endl;
    std::cout << "Normalization denominator = " << denom << std::endl;
    std::cout << "Yield needed for sigma = 150 microbarn = "
              << Y_for_150 << std::endl;
    std::cout << "========================================" << std::endl;

    Long64_t nEntries = tagger->GetEntries();

    // ============================================================
    // Timing histogram after pi0 mass cut
    // ============================================================

    TH1D *hDT_all = new TH1D(
        "hDT_all",
        "Timing spectrum after pi0 mass cut; t_{tagger}-t_{#pi^{0}}; counts",
        1000, 0, 1000
    );

    for (Long64_t iev = 0; iev < nEntries; iev++) {

        tagger->GetEntry(iev);
        pi0->GetEntry(iev);

        for (int itag = 0; itag < nTagged; itag++) {

            if (taggedChannel[itag] != CH) continue;

            for (int ip = 0; ip < nParticles; ip++) {

                double m = mass[ip];

                if (m < mPi0L || m > mPi0R) continue;

                double dt = taggedTime[itag] - time[ip];

                hDT_all->Fill(dt);
            }
        }
    }

    TCanvas *cDT = new TCanvas("cDT", "Timing spectrum", 900, 600);
    hDT_all->Draw();

    // Prompt lines
    TLine *pL = new TLine(promptL, 0, promptL, hDT_all->GetMaximum());
    TLine *pR = new TLine(promptR, 0, promptR, hDT_all->GetMaximum());

    pL->SetLineColor(kBlue);
    pR->SetLineColor(kBlue);
    pL->SetLineWidth(3);
    pR->SetLineWidth(3);

    pL->Draw("same");
    pR->Draw("same");

    // ============================================================
    // Scan random windows
    // ============================================================

    std::cout << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "RANDOM-WINDOW SCAN" << std::endl;
    std::cout << "========================================" << std::endl;

    for (int itest = 0; itest < NTEST; itest++) {

        double randomL = randomL_arr[itest];
        double randomR = randomR_arr[itest];

        double promptWidth = promptR - promptL;
        double randomWidth = randomR - randomL;

        double randomWeight = -promptWidth / randomWidth;

        TH1D *hME_prompt = new TH1D(
            Form("hME_prompt_test%d", itest),
            Form("Prompt ME; E_{miss} (MeV); counts"),
            300, -150, 150
        );

        TH1D *hME_random_scaled = new TH1D(
            Form("hME_random_scaled_test%d", itest),
            Form("Scaled random ME; E_{miss} (MeV); counts"),
            300, -150, 150
        );

        TH1D *hME_sub = new TH1D(
            Form("hME_sub_test%d", itest),
            Form("Prompt-random ME; E_{miss} (MeV); counts"),
            300, -150, 150
        );

        double nPromptComb = 0.0;
        double nRandomComb = 0.0;

        double nPi0Prompt = 0.0;
        double nPi0Random = 0.0;

        for (Long64_t iev = 0; iev < nEntries; iev++) {

            tagger->GetEntry(iev);
            pi0->GetEntry(iev);

            for (int itag = 0; itag < nTagged; itag++) {

                if (taggedChannel[itag] != CH) continue;

                for (int ip = 0; ip < nParticles; ip++) {

                    double m = mass[ip];

                    if (m < mPi0L || m > mPi0R) continue;

                    double dt = taggedTime[itag] - time[ip];

                    double weight = 0.0;
                    bool isPrompt = false;
                    bool isRandom = false;

                    if (dt >= promptL && dt <= promptR) {
                        weight = +1.0;
                        isPrompt = true;
                        nPromptComb++;
                        nPi0Prompt++;
                    }
                    else if (dt >= randomL && dt <= randomR) {
                        weight = randomWeight;
                        isRandom = true;
                        nRandomComb++;
                        nPi0Random++;
                    }
                    else {
                        continue;
                    }

                    // pi0 kinematics:
                    // clusterEnergy treated as pion momentum magnitude
                    double p = clusterEnergy[ip];
                    double Epi = sqrt(p*p + m*m);

                    double th = theta[ip] * TMath::DegToRad();
                    double ph = phi[ip] * TMath::DegToRad();

                    double px = p * sin(th) * cos(ph);
                    double py = p * sin(th) * sin(ph);
                    double pz = p * cos(th);

                    double pxmiss = -px;
                    double pymiss = -py;
                    double pzmiss = Egamma - pz;

                    double pmiss2 = pxmiss*pxmiss + pymiss*pymiss + pzmiss*pzmiss;

                    double Erecoil = sqrt(Mtarget*Mtarget + pmiss2);
                    double missingEnergy = Egamma + Mtarget - Epi - Erecoil;

                    if (isPrompt) {
                        hME_prompt->Fill(missingEnergy);
                    }

                    if (isRandom) {
                        hME_random_scaled->Fill(missingEnergy, -randomWeight);
                    }

                    hME_sub->Fill(missingEnergy, weight);
                }
            }
        }

        double Y_prompt = hME_prompt->Integral(
            hME_prompt->FindBin(meL),
            hME_prompt->FindBin(meR)
        );

        double Y_random_scaled = hME_random_scaled->Integral(
            hME_random_scaled->FindBin(meL),
            hME_random_scaled->FindBin(meR)
        );

        double Y_sub = hME_sub->Integral(
            hME_sub->FindBin(meL),
            hME_sub->FindBin(meR)
        );

        double sigma = 0.0;

        if (denom > 0) {
            sigma = Y_sub / denom;
        }

        double peakPrompt = hME_prompt->GetBinCenter(
            hME_prompt->GetMaximumBin()
        );

        double peakSub = hME_sub->GetBinCenter(
            hME_sub->GetMaximumBin()
        );

        std::cout << std::endl;
        std::cout << "Random window = "
                  << randomL << " to " << randomR << std::endl;

        std::cout << "Random width = " << randomWidth << std::endl;
        std::cout << "Random weight = " << randomWeight << std::endl;

        std::cout << "Prompt pi0 candidates = " << nPi0Prompt << std::endl;
        std::cout << "Random pi0 candidates = " << nPi0Random << std::endl;

        std::cout << "Net pi0 candidates = "
                  << nPi0Prompt + randomWeight*nPi0Random << std::endl;

        std::cout << "Prompt ME yield = " << Y_prompt << std::endl;
        std::cout << "Scaled random ME yield = " << Y_random_scaled << std::endl;
        std::cout << "Prompt-random ME yield = " << Y_sub << std::endl;

        if (Y_prompt > 0) {
            std::cout << "Fraction subtracted in signal window = "
                      << Y_random_scaled / Y_prompt << std::endl;
        }

        std::cout << "Cross section estimate = "
                  << sigma << " microbarn" << std::endl;

        std::cout << "Prompt ME peak = " << peakPrompt << " MeV" << std::endl;
        std::cout << "Prompt-random ME peak = " << peakSub << " MeV" << std::endl;

        // ------------------------------------------------------------
        // Draw prompt vs scaled random for each random window
        // ------------------------------------------------------------

        TCanvas *c = new TCanvas(
            Form("c_random_test%d", itest),
            Form("Random window %.0f-%.0f", randomL, randomR),
            1000, 700
        );

        hME_prompt->SetLineColor(kBlack);
        hME_prompt->SetLineWidth(2);

        hME_random_scaled->SetLineColor(kRed);
        hME_random_scaled->SetLineWidth(2);

        hME_sub->SetLineColor(kBlue);
        hME_sub->SetLineWidth(2);

        hME_prompt->Draw("hist");
        hME_random_scaled->Draw("hist same");
        hME_sub->Draw("hist same");

        TLegend *leg = new TLegend(0.58, 0.70, 0.88, 0.88);
        leg->AddEntry(hME_prompt, "Prompt", "l");
        leg->AddEntry(hME_random_scaled, "Scaled random", "l");
        leg->AddEntry(hME_sub, "Prompt - random", "l");
        leg->Draw();

        TLine *meLineL = new TLine(meL, 0, meL, hME_prompt->GetMaximum());
        TLine *meLineR = new TLine(meR, 0, meR, hME_prompt->GetMaximum());

        meLineL->SetLineColor(kGreen+2);
        meLineR->SetLineColor(kGreen+2);
        meLineL->SetLineWidth(2);
        meLineR->SetLineWidth(2);

        meLineL->Draw("same");
        meLineR->Draw("same");
    }

    std::cout << std::endl;
    std::cout << "Done. Inspect the canvases and compare the random-window dependence." << std::endl;
}
