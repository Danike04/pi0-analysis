/*
 * eff_MC_pi0.C
 *
 * PURPOSE
 *   Early Geant-level pi0 efficiency study. Sums deposited energy associated
 *   with the two decay photons, applies a per-photon threshold and reconstructs
 *   m_gg from the deposited energies and generated photon directions.
 *
 * DEFAULTS
 *   Photon deposited-energy threshold = 20 MeV.
 *   Reconstructed m_gg window = 110-145 MeV.
 *
 * OUTPUTS
 *   Efficiency summary on stdout, diagnostic histograms/canvases and
 *   eff_MC_pi0_output.root.
 *
 * RELATION TO PRODUCTION
 *   This was used to understand the detector-threshold/mass-cut efficiency.
 *   It does not reproduce the full Acqu reconstruction and is therefore kept
 *   under studies rather than as the production efficiency macro.
 *
 * USAGE EXAMPLE
 *   root -l -b -q 'eff_MC_pi0.C("Geant_4He_pi0-140-740MeV.root")'
 *
 * REPOSITORY STATUS
 *   Geant efficiency study; executable code unchanged.
 */
#include <iostream>
#include <cmath>

void eff_MC_pi0(const char* filename = "Geant_4He_pi0-140-740MeV.root")
{
    TFile* f = TFile::Open(filename);

    if (!f || f->IsZombie()) {
        std::cerr << "Errore: non riesco ad aprire il file " << filename << std::endl;
        return;
    }

    TTree* t = (TTree*) f->Get("h12");

    if (!t) {
        std::cerr << "Errore: TTree h12 non trovato" << std::endl;
        return;
    }

    const int MAXHITS = 1000;
    const int MAXPART = 10;

    Int_t nhits;
    Int_t npart;

    Float_t ecryst[MAXHITS];
    Int_t   pcryst[MAXHITS];

    Float_t felab[MAXPART];
    Float_t fklab[MAXPART];
    Float_t dircos[MAXPART][3];

    t->SetBranchAddress("nhits",  &nhits);
    t->SetBranchAddress("npart",  &npart);
    t->SetBranchAddress("ecryst", ecryst);
    t->SetBranchAddress("pcryst", pcryst);
    t->SetBranchAddress("elab",   felab);
    t->SetBranchAddress("klab",   fklab);
    t->SetBranchAddress("dircos", dircos);

    Long64_t Ngen = t->GetEntries();

    // Soglia minima di energia depositata per considerare un gamma visto.
    // ecryst sembra essere in GeV.
    // 0.020 GeV = 20 MeV.
    double Ethr = 0.020;

    Long64_t N_one_gamma_seen = 0;
    Long64_t N_two_gamma_seen = 0;
    Long64_t N_pi0_mass_cut_truth = 0;
    Long64_t N_two_gamma_and_mass_cut = 0;

    TH1D* hEdep1 = new TH1D("hEdep1", "Deposited energy gamma 1;E_{dep} [GeV];Counts", 200, 0, 1.0);
    TH1D* hEdep2 = new TH1D("hEdep2", "Deposited energy gamma 2;E_{dep} [GeV];Counts", 200, 0, 1.0);
    TH1D* hEtot  = new TH1D("hEtot",  "Total deposited energy;E_{dep,tot} [GeV];Counts", 200, 0, 1.0);

    TH1D* hMgg_truth = new TH1D("hMgg_truth", "Generated #gamma#gamma invariant mass;m_{#gamma#gamma} [MeV];Counts", 200, 0, 250);
    TH1D* hMgg_recoE = new TH1D("hMgg_recoE", "m_{#gamma#gamma} using deposited energies;m_{#gamma#gamma} [MeV];Counts", 200, 0, 250);

    for (Long64_t i = 0; i < Ngen; i++) {

        t->GetEntry(i);

        double Edep1 = 0.0;
        double Edep2 = 0.0;

        for (int j = 0; j < nhits; j++) {

            if (pcryst[j] == 1) {
                Edep1 += ecryst[j];
            }

            if (pcryst[j] == 2) {
                Edep2 += ecryst[j];
            }
        }

        hEdep1->Fill(Edep1);
        hEdep2->Fill(Edep2);
        hEtot->Fill(Edep1 + Edep2);

        bool gamma1_seen = Edep1 > Ethr;
        bool gamma2_seen = Edep2 > Ethr;

        if (gamma1_seen || gamma2_seen) {
            N_one_gamma_seen++;
        }

        if (gamma1_seen && gamma2_seen) {
            N_two_gamma_seen++;
        }

        // Calcolo dell'angolo tra i due gamma usando le direzioni generate
        double cos12 =
            dircos[0][0] * dircos[1][0] +
            dircos[0][1] * dircos[1][1] +
            dircos[0][2] * dircos[1][2];

        // Protezione numerica
        if (cos12 >  1.0) cos12 =  1.0;
        if (cos12 < -1.0) cos12 = -1.0;

        // Massa invariante usando le energie generate.
        // felab è in GeV, quindi converto il risultato in MeV.
        double E1_gen = felab[0];
        double E2_gen = felab[1];

        double mgg_truth_GeV = std::sqrt(2.0 * E1_gen * E2_gen * (1.0 - cos12));
        double mgg_truth_MeV = 1000.0 * mgg_truth_GeV;

        hMgg_truth->Fill(mgg_truth_MeV);

        if (mgg_truth_MeV > 110.0 && mgg_truth_MeV < 145.0) {
            N_pi0_mass_cut_truth++;
        }

        // Massa "ricostruita" molto approssimata:
        // uso energia depositata nei cristalli ma direzione vera dei gamma.
        // Questa NON è una vera ricostruzione a cluster, ma è utile come prima stima.
        if (gamma1_seen && gamma2_seen) {

            double mgg_reco_GeV = std::sqrt(2.0 * Edep1 * Edep2 * (1.0 - cos12));
            double mgg_reco_MeV = 1000.0 * mgg_reco_GeV;

            hMgg_recoE->Fill(mgg_reco_MeV);

            if (mgg_reco_MeV > 110.0 && mgg_reco_MeV < 145.0) {
                N_two_gamma_and_mass_cut++;
            }
        }
    }

    double eff_one = (double) N_one_gamma_seen / (double) Ngen;
    double err_one = std::sqrt(eff_one * (1.0 - eff_one) / (double) Ngen);

    double eff_two = (double) N_two_gamma_seen / (double) Ngen;
    double err_two = std::sqrt(eff_two * (1.0 - eff_two) / (double) Ngen);

    double eff_mass = (double) N_two_gamma_and_mass_cut / (double) Ngen;
    double err_mass = std::sqrt(eff_mass * (1.0 - eff_mass) / (double) Ngen);

    std::cout << "========================================" << std::endl;
    std::cout << "File MC = " << filename << std::endl;
    std::cout << "Generated events Ngen = " << Ngen << std::endl;
    std::cout << "Energy threshold = " << Ethr * 1000.0 << " MeV" << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Events with at least one gamma seen = " << N_one_gamma_seen << std::endl;
    std::cout << "Efficiency one gamma seen = "
              << eff_one << " +/- " << err_one
              << "  =  " << eff_one * 100.0 << " +/- " << err_one * 100.0 << " %" << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Events with both gammas seen = " << N_two_gamma_seen << std::endl;
    std::cout << "Detection efficiency two gammas = "
              << eff_two << " +/- " << err_two
              << "  =  " << eff_two * 100.0 << " +/- " << err_two * 100.0 << " %" << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Generated events with 110 < mgg_truth < 145 MeV = "
              << N_pi0_mass_cut_truth << std::endl;
    std::cout << "Events with both gammas seen and 110 < mgg_recoE < 145 MeV = "
              << N_two_gamma_and_mass_cut << std::endl;
    std::cout << "Efficiency with pi0 mass cut = "
              << eff_mass << " +/- " << err_mass
              << "  =  " << eff_mass * 100.0 << " +/- " << err_mass * 100.0 << " %" << std::endl;
    std::cout << "========================================" << std::endl;

    TCanvas* c1 = new TCanvas("c1", "Deposited energies", 1200, 800);
    c1->Divide(2,2);
    c1->cd(1);
    hEdep1->Draw();
    c1->cd(2);
    hEdep2->Draw();
    c1->cd(3);
    hEtot->Draw();
    c1->cd(4);
    hMgg_recoE->Draw();

    TCanvas* c2 = new TCanvas("c2", "Invariant masses", 1000, 700);
    hMgg_truth->SetLineColor(kBlue);
    hMgg_recoE->SetLineColor(kRed);
    hMgg_truth->Draw();
    hMgg_recoE->Draw("same");

    TLegend* leg = new TLegend(0.60, 0.70, 0.88, 0.88);
    leg->AddEntry(hMgg_truth, "Generated energies", "l");
    leg->AddEntry(hMgg_recoE, "Deposited energies", "l");
    leg->Draw();

    TFile* fout = new TFile("eff_MC_pi0_output.root", "RECREATE");
    hEdep1->Write();
    hEdep2->Write();
    hEtot->Write();
    hMgg_truth->Write();
    hMgg_recoE->Write();
    fout->Close();

    std::cout << "Istogrammi salvati in eff_MC_pi0_output.root" << std::endl;
}
