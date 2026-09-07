/*
 * check_timing_cleaning.C
 *
 * PURPOSE
 *   Check how simple event-cleaning requirements modify the tagger-pi0 timing
 *   spectrum in a representative GoAT run.
 *
 * WHAT THE MACRO ACTUALLY DOES
 *   - Opens GoAT_CBTagg_31832.root by default.
 *   - Uses tagger channel 51 and Egamma = 212 MeV.
 *   - Requires 110 < m(pi0) < 145 MeV and 30 < theta_pi0 < 150 degrees.
 *   - Compares the timing spectrum after: mass cut only; mass + missing-energy
 *     cut; mass + low tagger multiplicity; and the combined cleaning.
 *   - The missing-energy window is 15-65 MeV.
 *   - The historical code declares maxTaggedCut=60 but the comparison plotted
 *     in the macro uses nTagged <= 10; this behaviour is preserved.
 *
 * OUTPUT
 *   timing_cleaning_check.pdf and event-count summaries printed to stdout.
 *
 * USAGE
 *   root -l -b -q 'check_timing_cleaning.C'
 *
 * REPOSITORY STATUS
 *   Historical selection diagnostic. It documents an exploratory cleaning
 *   study and is not part of the production pi0 extraction. Code unchanged.
 */
void check_timing_cleaning() {

    const char* filename = "GoAT_CBTagg_31832.root";

    const int CH = 51;
    const double Egamma = 212.0;
    const double Mtarget = 3727.38;

    // Current cuts
    const double mPi0L = 110.0;
    const double mPi0R = 145.0;

    const double meL = 15.0;
    const double meR = 65.0;

    const int maxTaggedCut = 60;

    TFile *f = new TFile(filename);

    if (!f || f->IsZombie()) {
        std::cout << "Error: cannot open file " << filename << std::endl;
        return;
    }

    TTree *tagger = (TTree*)f->Get("tagger");
    TTree *pi0    = (TTree*)f->Get("neutralPions");

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

    TH1D *hDT_mass = new TH1D(
        "hDT_mass",
        "Timing: mass cut only; t_{tagger}-t_{#pi^{0}}; counts",
        1000, 100, 1000
    );

    TH1D *hDT_mass_ME = new TH1D(
        "hDT_mass_ME",
        "Timing: mass cut + missing-energy cut; t_{tagger}-t_{#pi^{0}}; counts",
        1000, 100, 1000
    );

    TH1D *hDT_mass_lowMult = new TH1D(
        "hDT_mass_lowMult",
        "Timing: mass cut + nTagged <= 10; t_{tagger}-t_{#pi^{0}}; counts",
        1000, 100, 1000
    );

    TH1D *hDT_clean = new TH1D(
        "hDT_clean",
        "Timing: mass cut + missing-energy cut + nTagged <= 10; t_{tagger}-t_{#pi^{0}}; counts",
        1000, 100, 1000
    );

    Long64_t nEntries = tagger->GetEntries();

    for (Long64_t iev = 0; iev < nEntries; iev++) {

        tagger->GetEntry(iev);
        pi0->GetEntry(iev);

        for (int it = 0; it < nTagged; it++) {

            if (taggedChannel[it] != CH) continue;

            for (int ip = 0; ip < nParticles; ip++) {

                double m = mass[ip];

                // First physical cut: pi0 mass
                if (m < mPi0L || m > mPi0R) continue;
if (theta[ip] < 30.0 || theta[ip] > 150.0) continue;

                double dt = taggedTime[it] - time[ip];

                // Compute missing energy
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

                bool passME = (missingEnergy >= meL && missingEnergy <= meR);
                bool passLowMult = (nTagged <= maxTaggedCut);

                hDT_mass->Fill(dt);

                if (passME) {
                    hDT_mass_ME->Fill(dt);
                }

                if (passLowMult) {
                    hDT_mass_lowMult->Fill(dt);
                }

                if (passME && passLowMult) {
                    hDT_clean->Fill(dt);
                }
            }
        }
    }

    // Rebin for readability only
    hDT_mass->Rebin(4);
    hDT_mass_ME->Rebin(4);
    hDT_mass_lowMult->Rebin(4);
    hDT_clean->Rebin(4);

    TCanvas *c = new TCanvas("c", "Timing cleaning check", 1200, 900);
    c->Divide(2,2);

    c->cd(1);
    hDT_mass->Draw("hist");

    c->cd(2);
    hDT_mass_ME->Draw("hist");

    c->cd(3);
    hDT_mass_lowMult->Draw("hist");

    c->cd(4);
    hDT_clean->Draw("hist");

    c->SaveAs("timing_cleaning_check.pdf");

    std::cout << "Saved timing_cleaning_check.pdf" << std::endl;

    std::cout << "Entries:" << std::endl;
    std::cout << "Mass cut only = " << hDT_mass->Integral() << std::endl;
    std::cout << "Mass + ME cut = " << hDT_mass_ME->Integral() << std::endl;
    std::cout << "Mass + nTagged <= 10 = " << hDT_mass_lowMult->Integral() << std::endl;
    std::cout << "Mass + ME + nTagged <= 10 = " << hDT_clean->Integral() << std::endl;
}
