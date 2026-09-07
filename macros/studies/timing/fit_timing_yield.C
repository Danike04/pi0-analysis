/*
 * fit_timing_yield.C
 *
 * PURPOSE
 *   Historical cross-check that extracts the prompt pi0 yield by fitting the
 *   tagger-pi0 timing spectrum after invariant-mass and missing-energy cuts,
 *   instead of using a fixed prompt/random sideband subtraction.
 *
 * WHAT THE MACRO ACTUALLY DOES
 *   - Opens GoAT_CBTagg_31838.root by default.
 *   - Uses tagger channel 51 and Egamma = 212 MeV.
 *   - Requires 110 < m(pi0) < 145 MeV and missing energy 30-80 MeV.
 *   - Builds the timing spectrum, fits the prompt structure in 400-750 ns and
 *     obtains a fitted prompt yield.
 *   - Prints a cross-section estimate using hard-coded eps_tag=0.3,
 *     eps_det=0.70 and target thickness 0.940e-7 microbarn^-1.
 *
 * OUTPUT
 *   timing_fit_yield.pdf plus fit/yield information printed to stdout.
 *
 * USAGE
 *   root -l -b -q 'fit_timing_yield.C'
 *
 * REPOSITORY STATUS
 *   Kept as a timing/yield cross-check only. The fixed constants and single-
 *   channel setup make it unsuitable as a production cross-section macro.
 *   The executable code below is unchanged from the raw file.
 */
void fit_timing_yield() {

    const char* filename = "GoAT_CBTagg_31838.root";

    const int CH = 51;
    const double Egamma = 212.0;
    const double Mtarget = 3727.38;

    // Cuts
    const double mPi0L = 110.0;
    const double mPi0R = 145.0;

    const double meL = 30.0;
    const double meR = 80.0;

    // Normalization
    const double eps_tag = 0.3;
    const double eps_det = 0.70;
    const double thickness = 0.940e-7; // microbarn^-1

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

    TH1D *hDT_signal = new TH1D(
        "hDT_signal",
        "Timing after #pi^{0} mass and missing-energy cuts; t_{tagger}-t_{#pi^{0}}; counts",
        450, 300, 900
    );

    Long64_t nEntries = tagger->GetEntries();

    for (Long64_t iev = 0; iev < nEntries; iev++) {

        tagger->GetEntry(iev);
        pi0->GetEntry(iev);

        for (int it = 0; it < nTagged; it++) {

            if (taggedChannel[it] != CH) continue;

            for (int ip = 0; ip < nParticles; ip++) {

                double m = mass[ip];

                if (m < mPi0L || m > mPi0R) continue;

                // pi0 kinematics
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

                if (missingEnergy < meL || missingEnergy > meR) continue;

                double dt = taggedTime[it] - time[ip];

                hDT_signal->Fill(dt);
            }
        }
    }

    // Rebin only to stabilize the fit visually/statistically
    hDT_signal->Rebin(3);

    // Fit range: around the prompt region and nearby background
    const double fitL = 400.0;
    const double fitR = 750.0;

    TF1 *fit = new TF1(
        "fit",
        "gaus(0) + pol1(3)",
        fitL,
        fitR
    );

    // Initial guesses
    double peakX = hDT_signal->GetBinCenter(hDT_signal->GetMaximumBin());
    double peakY = hDT_signal->GetMaximum();

    fit->SetParameter(0, peakY);    // Gaussian amplitude
    fit->SetParameter(1, peakX);    // Gaussian mean
    fit->SetParameter(2, 25.0);     // Gaussian sigma
    fit->SetParameter(3, 20.0);     // background offset
    fit->SetParameter(4, 0.0);      // background slope

    fit->SetParLimits(1, 520.0, 610.0);
    fit->SetParLimits(2, 5.0, 80.0);

    hDT_signal->Fit(fit, "R");

    double A     = fit->GetParameter(0);
    double mean  = fit->GetParameter(1);
    double sig_t = fit->GetParameter(2);

    double binWidth = hDT_signal->GetBinWidth(1);

    // Gaussian area in counts:
    // Integral of A exp(...) dt divided by bin width
    double Y_fit = A * sig_t * sqrt(2.0*TMath::Pi()) / binWidth;

    // Electron flux Ne
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
    double sigma_xs = 0.0;

    if (denom > 0) {
        sigma_xs = Y_fit / denom;
    }

    std::cout << "========================================" << std::endl;
    std::cout << "TIMING FIT YIELD EXTRACTION" << std::endl;
    std::cout << "========================================" << std::endl;

    std::cout << "Input file = " << filename << std::endl;
    std::cout << "Tagger channel = " << CH << std::endl;
    std::cout << "Photon energy = " << Egamma << " MeV" << std::endl;

    std::cout << "Pi0 mass cut = "
              << mPi0L << " to " << mPi0R << " MeV" << std::endl;

    std::cout << "Missing-energy cut = "
              << meL << " to " << meR << " MeV" << std::endl;

    std::cout << "Fit range = "
              << fitL << " to " << fitR << std::endl;

    std::cout << "Gaussian mean = " << mean << std::endl;
    std::cout << "Gaussian sigma_t = " << sig_t << std::endl;
    std::cout << "Fitted prompt yield Y_fit = " << Y_fit << std::endl;

    std::cout << "Electron flux Ne = " << Ne << std::endl;
    std::cout << "Normalization denominator = " << denom << std::endl;

    std::cout << "Cross section from timing fit = "
              << sigma_xs << " microbarn" << std::endl;

    std::cout << "========================================" << std::endl;

    TCanvas *c = new TCanvas("c", "Timing fit", 1000, 700);

    hDT_signal->SetLineColor(kBlack);
    hDT_signal->SetLineWidth(2);
    hDT_signal->Draw("E");

    fit->SetLineColor(kRed);
    fit->SetLineWidth(3);
    fit->Draw("same");

    c->SaveAs("timing_fit_yield.pdf");
}
