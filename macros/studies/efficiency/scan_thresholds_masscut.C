/*
 * scan_thresholds_masscut.C
 *
 * PURPOSE
 *   Scan the Geant deposited-energy threshold used to identify both pi0 decay
 *   photons and print the corresponding two-photon and m_gg-cut efficiencies.
 *
 * THRESHOLDS
 *   5, 10, 15, 20, 25, 30, 40, 50 and 60 MeV.
 *   The invariant-mass window is fixed to 110-145 MeV.
 *
 * INPUT
 *   Geant ROOT tree `h1` with nhits/ecryst/pcryst/dircos branches.
 *
 * OUTPUT
 *   Text summary on stdout only.
 *
 * USAGE EXAMPLE
 *   root -l -b -q 'scan_thresholds_masscut.C("Geant_4He_pi0.root")'
 *
 * REPOSITORY STATUS
 *   Useful detector-threshold systematic study; executable code unchanged.
 */
void scan_thresholds_masscut(const char* filename = "Geant_4He_pi0.root")
{
    TFile* f = TFile::Open(filename);
    TTree* t = (TTree*) f->Get("h12");

    const int MAXHITS = 1000;
    const int MAXPART = 10;

    Int_t nhits;
    Float_t ecryst[MAXHITS];
    Int_t pcryst[MAXHITS];

    Float_t dircos[MAXPART][3];

    t->SetBranchAddress("nhits", &nhits);
    t->SetBranchAddress("ecryst", ecryst);
    t->SetBranchAddress("pcryst", pcryst);
    t->SetBranchAddress("dircos", dircos);

    Long64_t Ngen = t->GetEntries();

    double thresholds[] = {0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.040, 0.050, 0.060};
    int nThr = sizeof(thresholds) / sizeof(double);

    std::cout << "Threshold scan with pi0 mass cut" << std::endl;
    std::cout << "Ngen = " << Ngen << std::endl;
    std::cout << "Mass cut: 110 < mgg < 145 MeV" << std::endl;
    std::cout << "----------------------------------------" << std::endl;

    for (int ithr = 0; ithr < nThr; ithr++) {

        double Ethr = thresholds[ithr];

        Long64_t Ntwo = 0;
        Long64_t Ntwo_masscut = 0;

        for (Long64_t i = 0; i < Ngen; i++) {

            t->GetEntry(i);

            double Edep1 = 0.0;
            double Edep2 = 0.0;

            for (int j = 0; j < nhits; j++) {
                if (pcryst[j] == 1) Edep1 += ecryst[j];
                if (pcryst[j] == 2) Edep2 += ecryst[j];
            }

            bool gamma1_seen = Edep1 > Ethr;
            bool gamma2_seen = Edep2 > Ethr;

            if (gamma1_seen && gamma2_seen) {

                Ntwo++;

                double cos12 =
                    dircos[0][0] * dircos[1][0] +
                    dircos[0][1] * dircos[1][1] +
                    dircos[0][2] * dircos[1][2];

                if (cos12 >  1.0) cos12 =  1.0;
                if (cos12 < -1.0) cos12 = -1.0;

                double mgg_GeV = sqrt(2.0 * Edep1 * Edep2 * (1.0 - cos12));
                double mgg_MeV = 1000.0 * mgg_GeV;

                if (mgg_MeV > 110.0 && mgg_MeV < 145.0) {
                    Ntwo_masscut++;
                }
            }
        }

        double eff_two = (double) Ntwo / (double) Ngen;
        double err_two = sqrt(eff_two * (1.0 - eff_two) / Ngen);

        double eff_mass = (double) Ntwo_masscut / (double) Ngen;
        double err_mass = sqrt(eff_mass * (1.0 - eff_mass) / Ngen);

        std::cout << "Ethr = " << Ethr * 1000.0 << " MeV"
                  << "   Ntwo = " << Ntwo
                  << "   eff_two = " << eff_two * 100.0
                  << " +/- " << err_two * 100.0 << " %"
                  << "   Nmass = " << Ntwo_masscut
                  << "   eff_mass = " << eff_mass * 100.0
                  << " +/- " << err_mass * 100.0 << " %"
                  << std::endl;
    }
}
