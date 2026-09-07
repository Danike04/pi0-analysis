/*
 * analyze_coherent_mc_root6.C
 *
 * PURPOSE
 *   Detailed reconstruction diagnostic for coherent 4He(gamma,pi0)4He
 *   Acqu/Geant MC. It follows reconstructed neutral clusters through m_gg,
 *   missing mass, missing energy, opening-angle/phi-min variables and paper
 *   energy/CM-angle bins.
 *
 * DEFAULT SELECTION
 *   cluster E >= 20 MeV, vetoEnergy <= 1 MeV,
 *   110 < m_gg < 155 MeV, -10 < MM-M_4He < 40 MeV.
 *
 * OUTPUTS
 *   coherent_mc_diagnostics.txt/.root and several summary PDFs. The ROOT file
 *   contains, among others, `coh_mm_after_mgg` and
 *   `coh_mm_vs_egamma_after_mgg`, consumed by compare_data_coherent_mm_root6.C.
 *
 * USAGE EXAMPLE
 *   root -l -b -q 'analyze_coherent_mc_root6.C("Acqu_geant_He4pi0.root")'
 *
 * REPOSITORY STATUS
 *   General coherent-MC diagnostic; executable code unchanged.
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
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"

int inc_find_energy_bin(double egamma,
                        const double eLow[],
                        const double eHigh[],
                        int nBins)
{
    for(int i=0; i<nBins; i++) {
        if(egamma >= eLow[i] && egamma <= eHigh[i]) return i;
    }
    return -1;
}

double inc_coherent_pion_lab_energy(double egamma,
                                    double thetaCmDeg,
                                    double mpi0,
                                    double targetMass)
{
    const double s = targetMass * targetMass + 2.0 * egamma * targetMass;
    const double sqrtS = std::sqrt(s);
    const double ePiCm = (s + mpi0 * mpi0 - targetMass * targetMass) / (2.0 * sqrtS);

    double pPiCm2 = ePiCm * ePiCm - mpi0 * mpi0;
    if(pPiCm2 < 0.0) pPiCm2 = 0.0;
    const double pPiCm = std::sqrt(pPiCm2);

    const double betaCm = egamma / (egamma + targetMass);
    const double gammaCm = 1.0 / std::sqrt(1.0 - betaCm * betaCm);
    const double cosThetaCm = std::cos(thetaCmDeg * TMath::DegToRad());

    return gammaCm * (ePiCm + betaCm * pPiCm * cosThetaCm);
}

double inc_phi_min_deg(double ePiLab, double mpi0)
{
    if(ePiLab <= mpi0) return 180.0;

    double arg = std::sqrt(ePiLab * ePiLab - mpi0 * mpi0) / ePiLab;
    if(arg < -1.0) arg = -1.0;
    if(arg > 1.0) arg = 1.0;

    return 2.0 * std::acos(arg) * TMath::RadToDeg();
}

double inc_safe_fraction(Long64_t num, Long64_t den)
{
    if(den <= 0) return 0.0;
    return double(num) / double(den);
}

void analyze_coherent_mc_root6(const char* fname="Acqu_geant_He4pi0.root",
                                 Long64_t maxEvents=-1,
                                 double clusterEnergyMin=20.0,
                                 double vetoEnergyMax=1.0,
                                 double mggMin=110.0,
                                 double mggMax=155.0,
                                 double mmMin=-10.0,
                                 double mmMax=40.0,
                                 double thetaCmMinDeg=5.0,
                                 double thetaCmMaxDeg=150.0,
                                 double thetaCmBinWidthDeg=5.0)
{
    const int nCh = 352;
    const int maxTracks = 512;
    const int maxTagged = 2048;

    const double beamEMax = 855.0;
    const double targetMass = 3727.38;
    const double mpi0 = 134.9768;
    const double deg = TMath::DegToRad();

    const double eLowPaper[]  = {201,211,223,235,247,259,271,283,295,309,319,331,343,356,367,379,391};
    const double eHighPaper[] = {210,222,234,246,258,270,282,294,308,318,330,342,355,366,378,390,401};
    const int nEnergyBins = 17;

    if(thetaCmBinWidthDeg <= 0.0) thetaCmBinWidthDeg = 5.0;
    if(thetaCmMaxDeg <= thetaCmMinDeg) {
        thetaCmMinDeg = 5.0;
        thetaCmMaxDeg = 150.0;
    }
    int nThetaBins = int(std::floor((thetaCmMaxDeg-thetaCmMinDeg)/thetaCmBinWidthDeg + 0.5));
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
        std::cout << "Missing tracks, tagger or setupParameters tree" << std::endl;
        f->Close();
        return;
    }

    if(!tracks->GetBranch("nTracks") ||
       !tracks->GetBranch("clusterEnergy") ||
       !tracks->GetBranch("theta") ||
       !tracks->GetBranch("phi") ||
       !tracks->GetBranch("vetoEnergy") ||
       !tagger->GetBranch("nTagged") ||
       !tagger->GetBranch("taggedChannel") ||
       !setup->GetBranch("nTagger") ||
       !setup->GetBranch("TaggerPhotonEnergy")) {
        std::cout << "One or more required branches are missing" << std::endl;
        f->Close();
        return;
    }

    int nTracks = 0;
    double clusterEnergy[maxTracks];
    double theta[maxTracks];
    double phi[maxTracks];
    double vetoEnergy[maxTracks];

    tracks->SetBranchAddress("nTracks", &nTracks);
    tracks->SetBranchAddress("clusterEnergy", clusterEnergy);
    tracks->SetBranchAddress("theta", theta);
    tracks->SetBranchAddress("phi", phi);
    tracks->SetBranchAddress("vetoEnergy", vetoEnergy);

    int nTagged = 0;
    int taggedChannel[maxTagged];
    tagger->SetBranchAddress("nTagged", &nTagged);
    tagger->SetBranchAddress("taggedChannel", taggedChannel);

    int nTagger = 0;
    double taggerPhotonEnergy[nCh];
    setup->SetBranchAddress("nTagger", &nTagger);
    setup->SetBranchAddress("TaggerPhotonEnergy", taggerPhotonEnergy);
    setup->GetEntry(0);
    const int nValidCh = std::min(nTagger, nCh);

    Long64_t nentries = tracks->GetEntries();
    if(tagger->GetEntries() < nentries) nentries = tagger->GetEntries();
    if(maxEvents > 0 && maxEvents < nentries) nentries = maxEvents;

    TH1D* hNTracks = new TH1D("coh_ntracks", "Coherent MC reconstructed multiplicity;nTracks;Events", 20, -0.5, 19.5);
    TH1D* hMgg = new TH1D("coh_mgg_best", "Coherent MC best #gamma#gamma pair;m_{#gamma#gamma} (MeV);Events", 300, 0.0, 300.0);

    TH1D* hMMBefore = new TH1D("coh_mm_before_mgg", "Coherent MC missing mass before m_{#gamma#gamma} cut;MM-M_{^{4}He} (MeV);Tag combinations", 300, -150.0, 150.0);
    TH1D* hMMAfterMgg = new TH1D("coh_mm_after_mgg", "Coherent MC missing mass after m_{#gamma#gamma} cut;MM-M_{^{4}He} (MeV);Tag combinations", 300, -150.0, 150.0);
    TH1D* hMMAfterAll = new TH1D("coh_mm_after_mgg_mm", "Coherent MC after m_{#gamma#gamma} and MM cuts;MM-M_{^{4}He} (MeV);Tag combinations", 300, -150.0, 150.0);

    TH1D* hDeltaEBefore = new TH1D("coh_deltaE_before_mgg", "Coherent MC missing energy before m_{#gamma#gamma} cut;#DeltaE_{#pi^{0}} (MeV);Tag combinations", 300, -200.0, 100.0);
    TH1D* hDeltaEAfterMgg = new TH1D("coh_deltaE_after_mgg", "Coherent MC missing energy after m_{#gamma#gamma} cut;#DeltaE_{#pi^{0}} (MeV);Tag combinations", 300, -200.0, 100.0);
    TH1D* hDeltaEAfterAll = new TH1D("coh_deltaE_after_mgg_mm", "Coherent MC missing energy after m_{#gamma#gamma} and MM cuts;#DeltaE_{#pi^{0}} (MeV);Tag combinations", 300, -200.0, 100.0);

    TH1D* hDeltaPhiAfterMgg = new TH1D("coh_delta_phi_after_mgg", "Coherent MC opening-angle excess after m_{#gamma#gamma} cut;#Phi_{#gamma#gamma}-#Phi_{min}^{coh} (deg);Tag combinations", 220, -20.0, 90.0);
    TH1D* hDeltaPhiAfterAll = new TH1D("coh_delta_phi_after_mgg_mm", "Coherent MC opening-angle excess after m_{#gamma#gamma} and MM cuts;#Phi_{#gamma#gamma}-#Phi_{min}^{coh} (deg);Tag combinations", 220, -20.0, 90.0);

    TH2D* hMMDeltaEBefore = new TH2D("coh_mm_vs_deltaE_before_mgg", "Coherent MC MM vs missing energy before m_{#gamma#gamma} cut;#DeltaE_{#pi^{0}} (MeV);MM-M_{^{4}He} (MeV)", 300, -200.0, 100.0, 300, -150.0, 150.0);
    TH2D* hMMDeltaEAfterMgg = new TH2D("coh_mm_vs_deltaE_after_mgg", "Coherent MC MM vs missing energy after m_{#gamma#gamma} cut;#DeltaE_{#pi^{0}} (MeV);MM-M_{^{4}He} (MeV)", 300, -200.0, 100.0, 300, -150.0, 150.0);
    TH2D* hMMDeltaEAfterAll = new TH2D("coh_mm_vs_deltaE_after_mgg_mm", "Coherent MC MM vs missing energy after all current cuts;#DeltaE_{#pi^{0}} (MeV);MM-M_{^{4}He} (MeV)", 300, -200.0, 100.0, 300, -150.0, 150.0);

    TH2D* hOpeningPhiMinBefore = new TH2D("coh_opening_vs_phiMin_before_mgg", "Coherent MC opening angle vs coherent #Phi_{min};#Phi_{min}^{coh} (deg);#Phi_{#gamma#gamma} (deg)", 180, 0.0, 180.0, 180, 0.0, 180.0);
    TH2D* hOpeningPhiMinAfterMgg = new TH2D("coh_opening_vs_phiMin_after_mgg", "Coherent MC opening angle vs coherent #Phi_{min} after m_{#gamma#gamma};#Phi_{min}^{coh} (deg);#Phi_{#gamma#gamma} (deg)", 180, 0.0, 180.0, 180, 0.0, 180.0);
    TH2D* hOpeningPhiMinAfterAll = new TH2D("coh_opening_vs_phiMin_after_mgg_mm", "Coherent MC opening angle vs coherent #Phi_{min} after all current cuts;#Phi_{min}^{coh} (deg);#Phi_{#gamma#gamma} (deg)", 180, 0.0, 180.0, 180, 0.0, 180.0);

    TH2D* hDeltaEEgammaAfterMgg = new TH2D("coh_deltaE_vs_egamma_after_mgg", "Coherent MC missing energy vs E_{#gamma} after m_{#gamma#gamma};E_{#gamma} (MeV);#DeltaE_{#pi^{0}} (MeV)", 220, 0.0, beamEMax, 300, -200.0, 100.0);
    TH2D* hMMEgammaAfterMgg = new TH2D("coh_mm_vs_egamma_after_mgg", "Coherent MC missing mass vs E_{#gamma} after m_{#gamma#gamma};E_{#gamma} (MeV);MM-M_{^{4}He} (MeV)", 220, 0.0, beamEMax, 300, -150.0, 150.0);
    TH2D* hThetaEgammaAfterMgg = new TH2D("coh_theta_cm_vs_egamma_after_mgg", "Coherent MC reconstructed #theta_{cm} vs E_{#gamma} after m_{#gamma#gamma};E_{#gamma} (MeV);#theta_{#pi^{0}}^{cm} (deg)", 220, 0.0, beamEMax, 180, 0.0, 180.0);

    TH2D* hPaperMgg = new TH2D("coh_counts_after_mgg_paper_bins", "Coherent MC counts after m_{#gamma#gamma};E_{#gamma} bin;#theta_{#pi^{0}}^{cm} (deg)", nEnergyBins, 0.0, double(nEnergyBins), nThetaBins, thetaCmMinDeg, thetaCmMinDeg+nThetaBins*thetaCmBinWidthDeg);
    TH2D* hPaperAll = new TH2D("coh_counts_after_mgg_mm_paper_bins", "Coherent MC counts after m_{#gamma#gamma} and MM;E_{#gamma} bin;#theta_{#pi^{0}}^{cm} (deg)", nEnergyBins, 0.0, double(nEnergyBins), nThetaBins, thetaCmMinDeg, thetaCmMinDeg+nThetaBins*thetaCmBinWidthDeg);
    TH2D* hPaperLeak = new TH2D("coh_fraction_mm_given_mgg_paper_bins", "Fraction of reconstructed coherent MC surviving MM cut;E_{#gamma} bin;#theta_{#pi^{0}}^{cm} (deg)", nEnergyBins, 0.0, double(nEnergyBins), nThetaBins, thetaCmMinDeg, thetaCmMinDeg+nThetaBins*thetaCmBinWidthDeg);

    for(int ie=0; ie<nEnergyBins; ie++) {
        const TString label = Form("%.0f-%.0f", eLowPaper[ie], eHighPaper[ie]);
        hPaperMgg->GetXaxis()->SetBinLabel(ie+1, label.Data());
        hPaperAll->GetXaxis()->SetBinLabel(ie+1, label.Data());
        hPaperLeak->GetXaxis()->SetBinLabel(ie+1, label.Data());
    }

    std::vector<Long64_t> nPairByE(nEnergyBins, 0);
    std::vector<Long64_t> nMggByE(nEnergyBins, 0);
    std::vector<Long64_t> nAllByE(nEnergyBins, 0);

    Long64_t nEventsWithPair = 0;
    Long64_t nEventsPassMgg = 0;
    Long64_t nTagCombBefore = 0;
    Long64_t nTagCombPassMgg = 0;
    Long64_t nTagCombPassAll = 0;
    Long64_t nInvalidTagChannel = 0;

    TLorentzVector target(0.0, 0.0, 0.0, targetMass);

    for(Long64_t ev=0; ev<nentries; ev++) {
        tracks->GetEntry(ev);
        tagger->GetEntry(ev);
        hNTracks->Fill(nTracks);

        int nt = nTracks;
        if(nt > maxTracks) nt = maxTracks;

        int bestI = -1;
        int bestJ = -1;
        double bestM = -1.0;
        double bestDist = 1.0e30;
        double bestOpeningDeg = -1.0;
        TLorentzVector bestPi0;

        for(int i=0; i<nt; i++) {
            if(clusterEnergy[i] < clusterEnergyMin) continue;
            if(vetoEnergy[i] > vetoEnergyMax) continue;

            for(int j=i+1; j<nt; j++) {
                if(clusterEnergy[j] < clusterEnergyMin) continue;
                if(vetoEnergy[j] > vetoEnergyMax) continue;

                const double th1 = theta[i] * deg;
                const double th2 = theta[j] * deg;
                const double ph1 = phi[i] * deg;
                const double ph2 = phi[j] * deg;

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

                const TLorentzVector pi0 = g1 + g2;
                const double mgg = pi0.M();
                if(mgg <= 0.0) continue;

                const double dist = std::fabs(mgg-mpi0);
                if(dist < bestDist) {
                    bestDist = dist;
                    bestM = mgg;
                    bestI = i;
                    bestJ = j;
                    bestPi0 = pi0;
                    bestOpeningDeg = g1.Angle(g2.Vect()) * TMath::RadToDeg();
                }
            }
        }

        if(bestI < 0 || bestJ < 0) continue;
        nEventsWithPair++;
        hMgg->Fill(bestM);

        const bool passMgg = (bestM >= mggMin && bestM <= mggMax);
        if(passMgg) nEventsPassMgg++;

        int ntag = nTagged;
        if(ntag > maxTagged) ntag = maxTagged;

        for(int it=0; it<ntag; it++) {
            const int ch = taggedChannel[it];
            if(ch < 0 || ch >= nValidCh) {
                nInvalidTagChannel++;
                continue;
            }

            const double egamma = taggerPhotonEnergy[ch];
            if(egamma <= 0.0 || egamma > beamEMax+5.0) continue;

            const int ie = inc_find_energy_bin(egamma, eLowPaper, eHighPaper, nEnergyBins);

            TLorentzVector beam(0.0, 0.0, egamma, egamma);
            TLorentzVector recoil = beam + target - bestPi0;
            const double mm = recoil.M() - targetMass;

            TLorentzVector totalCM = beam + target;
            TLorentzVector pi0CM = bestPi0;
            pi0CM.Boost(-totalCM.BoostVector());
            const double thetaCmDeg = pi0CM.Theta() * TMath::RadToDeg();

            const double s = targetMass*targetMass + 2.0*egamma*targetMass;
            const double sqrtS = std::sqrt(s);
            const double ePiCmExpected = (s + mpi0*mpi0 - targetMass*targetMass)/(2.0*sqrtS);
            const double betaCm = egamma/(egamma+targetMass);
            const double gammaCm = 1.0/std::sqrt(1.0-betaCm*betaCm);
            const double ePiCmMeasured = gammaCm*(bestPi0.E() - betaCm*bestPi0.Pz());
            const double deltaE = ePiCmMeasured - ePiCmExpected;

            const double ePiLabCoh = inc_coherent_pion_lab_energy(egamma, thetaCmDeg, mpi0, targetMass);
            const double phiMinDeg = inc_phi_min_deg(ePiLabCoh, mpi0);
            const double deltaPhi = bestOpeningDeg - phiMinDeg;

            nTagCombBefore++;
            if(ie >= 0) nPairByE[ie]++;

            hMMBefore->Fill(mm);
            hDeltaEBefore->Fill(deltaE);
            hMMDeltaEBefore->Fill(deltaE, mm);
            hOpeningPhiMinBefore->Fill(phiMinDeg, bestOpeningDeg);

            if(!passMgg) continue;
            nTagCombPassMgg++;
            if(ie >= 0) nMggByE[ie]++;

            hMMAfterMgg->Fill(mm);
            hDeltaEAfterMgg->Fill(deltaE);
            hMMDeltaEAfterMgg->Fill(deltaE, mm);
            hOpeningPhiMinAfterMgg->Fill(phiMinDeg, bestOpeningDeg);
            hDeltaPhiAfterMgg->Fill(deltaPhi);
            hDeltaEEgammaAfterMgg->Fill(egamma, deltaE);
            hMMEgammaAfterMgg->Fill(egamma, mm);
            hThetaEgammaAfterMgg->Fill(egamma, thetaCmDeg);

            int ith = -1;
            if(thetaCmDeg >= thetaCmMinDeg && thetaCmDeg < thetaCmMaxDeg) {
                ith = int((thetaCmDeg-thetaCmMinDeg)/thetaCmBinWidthDeg);
            }
            if(ie >= 0 && ith >= 0 && ith < nThetaBins) {
                hPaperMgg->Fill(double(ie)+0.5, thetaCmDeg);
            }

            const bool passMM = (mm >= mmMin && mm <= mmMax);
            if(!passMM) continue;

            nTagCombPassAll++;
            if(ie >= 0) nAllByE[ie]++;

            hMMAfterAll->Fill(mm);
            hDeltaEAfterAll->Fill(deltaE);
            hMMDeltaEAfterAll->Fill(deltaE, mm);
            hOpeningPhiMinAfterAll->Fill(phiMinDeg, bestOpeningDeg);
            hDeltaPhiAfterAll->Fill(deltaPhi);

            if(ie >= 0 && ith >= 0 && ith < nThetaBins) {
                hPaperAll->Fill(double(ie)+0.5, thetaCmDeg);
            }
        }

        if(ev % 10000 == 0) {
            std::cout << "event " << ev << "/" << nentries << std::endl;
        }
    }

    for(int bx=1; bx<=nEnergyBins; bx++) {
        for(int by=1; by<=nThetaBins; by++) {
            const double den = hPaperMgg->GetBinContent(bx, by);
            const double num = hPaperAll->GetBinContent(bx, by);
            if(den <= 0.0) continue;
            const double frac = num/den;
            double err2 = frac*(1.0-frac)/den;
            if(err2 < 0.0) err2 = 0.0;
            hPaperLeak->SetBinContent(bx, by, frac);
            hPaperLeak->SetBinError(bx, by, std::sqrt(err2));
        }
    }

    std::ofstream out("coherent_mc_diagnostics.txt");
    out << "# Coherent MC diagnostic analysis\n";
    out << "# input = " << fname << "\n";
    out << "# processed events = " << nentries << "\n";
    out << "# clusterEnergyMin = " << clusterEnergyMin << " MeV\n";
    out << "# vetoEnergyMax = " << vetoEnergyMax << "\n";
    out << "# mgg cut = [" << mggMin << "," << mggMax << "] MeV\n";
    out << "# MM cut = [" << mmMin << "," << mmMax << "] MeV\n";
    out << "# events_with_pair = " << nEventsWithPair << "\n";
    out << "# events_pass_mgg = " << nEventsPassMgg << "\n";
    out << "# tag_combinations_before = " << nTagCombBefore << "\n";
    out << "# tag_combinations_pass_mgg = " << nTagCombPassMgg << "\n";
    out << "# tag_combinations_pass_mgg_mm = " << nTagCombPassAll << "\n";
    out << "# fraction_events_with_pair = " << inc_safe_fraction(nEventsWithPair,nentries) << "\n";
    out << "# fraction_mgg_given_pair = " << inc_safe_fraction(nEventsPassMgg,nEventsWithPair) << "\n";
    out << "# fraction_mm_given_mgg_tagcomb = " << inc_safe_fraction(nTagCombPassAll,nTagCombPassMgg) << "\n";
    out << "# IMPORTANT: this is a reconstructed-MC survival fraction, not an absolute incoherent contamination in data.\n";
    out << "# E_low E_high pair_tagcomb mgg_tagcomb mgg_mm_tagcomb frac_mm_given_mgg binomial_error\n";

    for(int ie=0; ie<nEnergyBins; ie++) {
        const double frac = inc_safe_fraction(nAllByE[ie], nMggByE[ie]);
        double err = 0.0;
        if(nMggByE[ie] > 0) {
            double err2 = frac*(1.0-frac)/double(nMggByE[ie]);
            if(err2 < 0.0) err2 = 0.0;
            err = std::sqrt(err2);
        }
        out << eLowPaper[ie] << " "
            << eHighPaper[ie] << " "
            << nPairByE[ie] << " "
            << nMggByE[ie] << " "
            << nAllByE[ie] << " "
            << frac << " "
            << err << "\n";
    }
    out.close();

    gStyle->SetOptStat(0);

    TCanvas* cSummary = new TCanvas("c_inc_summary", "Coherent MC summary", 1500, 900);
    cSummary->Divide(3,2);
    cSummary->cd(1); hMgg->Draw("HIST");
    cSummary->cd(2); hMMAfterMgg->Draw("HIST");
    cSummary->cd(3); hDeltaEAfterMgg->Draw("HIST");
    cSummary->cd(4); hMMDeltaEAfterMgg->Draw("COLZ");
    cSummary->cd(5); hOpeningPhiMinAfterMgg->Draw("COLZ");
    cSummary->cd(6); hDeltaPhiAfterMgg->Draw("HIST");
    cSummary->SaveAs("coherent_mc_summary.pdf");

    TCanvas* cAfterAll = new TCanvas("c_inc_after_all", "Coherent MC after current cuts", 1500, 500);
    cAfterAll->Divide(3,1);
    cAfterAll->cd(1); hMMAfterAll->Draw("HIST");
    cAfterAll->cd(2); hDeltaEAfterAll->Draw("HIST");
    cAfterAll->cd(3); hOpeningPhiMinAfterAll->Draw("COLZ");
    cAfterAll->SaveAs("coherent_mc_after_mgg_mm.pdf");

    TCanvas* cEnergy = new TCanvas("c_inc_energy_dependence", "Coherent MC energy dependence", 1500, 500);
    cEnergy->Divide(3,1);
    cEnergy->cd(1); hDeltaEEgammaAfterMgg->Draw("COLZ");
    cEnergy->cd(2); hMMEgammaAfterMgg->Draw("COLZ");
    cEnergy->cd(3); hThetaEgammaAfterMgg->Draw("COLZ");
    cEnergy->SaveAs("coherent_mc_energy_dependence.pdf");

    TCanvas* cPaper = new TCanvas("c_inc_paper_bins", "Coherent MC paper bins", 1600, 500);
    cPaper->Divide(3,1);
    cPaper->cd(1); hPaperMgg->GetXaxis()->LabelsOption("v"); hPaperMgg->Draw("COLZ");
    cPaper->cd(2); hPaperAll->GetXaxis()->LabelsOption("v"); hPaperAll->Draw("COLZ");
    cPaper->cd(3); hPaperLeak->GetXaxis()->LabelsOption("v"); hPaperLeak->SetMinimum(0.0); hPaperLeak->SetMaximum(1.0); hPaperLeak->Draw("COLZ TEXT");
    cPaper->SaveAs("coherent_mc_paper_bins.pdf");

    TFile* fout = new TFile("coherent_mc_diagnostics.root", "RECREATE");
    hNTracks->Write();
    hMgg->Write();
    hMMBefore->Write();
    hMMAfterMgg->Write();
    hMMAfterAll->Write();
    hDeltaEBefore->Write();
    hDeltaEAfterMgg->Write();
    hDeltaEAfterAll->Write();
    hDeltaPhiAfterMgg->Write();
    hDeltaPhiAfterAll->Write();
    hMMDeltaEBefore->Write();
    hMMDeltaEAfterMgg->Write();
    hMMDeltaEAfterAll->Write();
    hOpeningPhiMinBefore->Write();
    hOpeningPhiMinAfterMgg->Write();
    hOpeningPhiMinAfterAll->Write();
    hDeltaEEgammaAfterMgg->Write();
    hMMEgammaAfterMgg->Write();
    hThetaEgammaAfterMgg->Write();
    hPaperMgg->Write();
    hPaperAll->Write();
    hPaperLeak->Write();
    fout->Close();

    std::cout << "\nInput: " << fname << std::endl;
    std::cout << "Processed events = " << nentries << std::endl;
    std::cout << "Events with reconstructed gamma-gamma pair = " << nEventsWithPair << std::endl;
    std::cout << "Events passing mgg cut = " << nEventsPassMgg << std::endl;
    std::cout << "Tag combinations before cuts = " << nTagCombBefore << std::endl;
    std::cout << "Tag combinations passing mgg = " << nTagCombPassMgg << std::endl;
    std::cout << "Tag combinations passing mgg + MM = " << nTagCombPassAll << std::endl;
    std::cout << "Fraction with pair = " << inc_safe_fraction(nEventsWithPair,nentries) << std::endl;
    std::cout << "Fraction passing mgg among pair events = " << inc_safe_fraction(nEventsPassMgg,nEventsWithPair) << std::endl;
    std::cout << "Fraction passing MM among mgg tag combinations = " << inc_safe_fraction(nTagCombPassAll,nTagCombPassMgg) << std::endl;
    std::cout << "Invalid tag channels skipped = " << nInvalidTagChannel << std::endl;
    std::cout << "Saved: coherent_mc_diagnostics.root, coherent_mc_diagnostics.txt, and four PDFs" << std::endl;
    std::cout << "NOTE: no scaler/tagging-efficiency normalization is applied; this macro studies shapes and survival fractions only." << std::endl;

    f->Close();
}
