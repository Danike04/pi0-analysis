/*
 * paper_fig2_missing_energy_root6.C
 *
 * PURPOSE
 *   Reconstruct the pi0 -> gamma gamma candidate in AcquRoot data and in a
 *   coherent 4He(gamma,pi0)4He Monte Carlo sample, then compare the
 *   missing-energy distributions before and after an opening-angle cut.
 *
 * WHAT THE MACRO ACTUALLY DOES
 *   - Reads the trees: tracks, tagger, setupParameters.
 *   - Selects neutral clusters using clusterEnergy and vetoEnergy.
 *   - Among all accepted neutral-cluster pairs, keeps the pair whose
 *     invariant mass is closest to the pi0 mass.
 *   - Requires mggMin < m(gamma gamma) < mggMax.
 *   - For experimental data, performs prompt-random subtraction using the
 *     taggedTime windows passed to the main function.
 *   - Calculates DeltaE_pi0 in the gamma+4He centre-of-mass frame.
 *   - Calculates DeltaPhi = measured gamma-gamma opening angle - the minimum
 *     opening angle expected for coherent pi0 production.
 *   - Fills four photon-energy intervals: 223-234, 283-294, 319-330,
 *     and 356-366 MeV.
 *   - Produces histograms before and after DeltaPhi < deltaPhiCut.
 *   - Scales the coherent MC shape to the data only for plotting.
 *
 * IMPORTANT
 *   No missing-mass cut is applied in this macro.
 *
 * INPUT EXPECTATION
 *   dataFile and coherentMCFile must each contain the energy intervals used
 *   above. The default filenames are historical working-directory defaults.
 *
 * OUTPUTS (fixed names in the original workflow)
 *   paper_fig2_missing_energy.root
 *   paper_fig2_missing_energy.pdf
 *   paper_fig2_missing_energy.txt
 *
 * HISTORICAL WORKFLOW
 *   The macro was run separately for DeltaPhi = 8, 10, and 12 degrees. After
 *   each run the ROOT output was manually renamed to:
 *     paper_fig2_missing_energy_cut8.root
 *     paper_fig2_missing_energy_cut10.root
 *     paper_fig2_missing_energy_cut12.root
 *   These three files are inputs to quantify_fig2_cut_stability_root6.C.
 *
 * USAGE EXAMPLE
 *   root -l -b -q 'paper_fig2_missing_energy_root6.C(
 *       "DATA.root","COHERENT_MC.root",-1,-1,10.0)'
 *
 * NOTE FOR THE REPOSITORY
 *   The analysis logic, function signature, defaults, histogram names and
 *   output names below are preserved from the original macro. Documentation
 *   has been added without changing its behaviour.
 */
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"
#include "TVirtualPad.h"

namespace PaperFig2 {

const int kNCh = 352;
const int kMaxTracks = 512;
const int kMaxTagged = 2048;
const int kNBins = 4;

const double kELow[kNBins]  = {223.0, 283.0, 319.0, 356.0};
const double kEHigh[kNBins] = {234.0, 294.0, 330.0, 366.0};
const char* kLabels[kNBins] = {"224 MeV", "294 MeV", "320 MeV", "366 MeV"};

int FindEnergyBin(double egamma)
{
    for(int i=0; i<kNBins; ++i) {
        if(egamma >= kELow[i] && egamma <= kEHigh[i]) return i;
    }
    return -1;
}

double CoherentPionLabEnergy(double egamma,
                             double thetaCmDeg,
                             double mpi0,
                             double targetMass)
{
    const double s = targetMass*targetMass + 2.0*egamma*targetMass;
    const double sqrtS = std::sqrt(s);
    const double ePiCm = (s + mpi0*mpi0 - targetMass*targetMass)/(2.0*sqrtS);

    double pPiCm2 = ePiCm*ePiCm - mpi0*mpi0;
    if(pPiCm2 < 0.0) pPiCm2 = 0.0;
    const double pPiCm = std::sqrt(pPiCm2);

    const double betaCm = egamma/(egamma + targetMass);
    const double gammaCm = 1.0/std::sqrt(1.0 - betaCm*betaCm);
    const double cosThetaCm = std::cos(thetaCmDeg*TMath::DegToRad());

    return gammaCm*(ePiCm + betaCm*pPiCm*cosThetaCm);
}

double PhiMinDeg(double ePiLab, double mpi0)
{
    if(ePiLab <= mpi0) return 180.0;

    double arg = std::sqrt(ePiLab*ePiLab - mpi0*mpi0)/ePiLab;
    if(arg < -1.0) arg = -1.0;
    if(arg > 1.0) arg = 1.0;

    return 2.0*std::acos(arg)*TMath::RadToDeg();
}

bool BuildBestPi0(int nTracks,
                  const double clusterEnergy[],
                  const double theta[],
                  const double phi[],
                  const double vetoEnergy[],
                  double clusterEnergyMin,
                  double vetoEnergyMax,
                  double mpi0,
                  TLorentzVector& bestPi0,
                  double& bestMgg,
                  double& bestOpeningDeg)
{
    const int nt = std::min(nTracks, kMaxTracks);
    const double deg = TMath::DegToRad();

    double bestDist = 1.0e30;
    bool found = false;
    bestMgg = -1.0;
    bestOpeningDeg = -1.0;

    for(int i=0; i<nt; ++i) {
        if(clusterEnergy[i] < clusterEnergyMin) continue;
        if(vetoEnergy[i] > vetoEnergyMax) continue;

        for(int j=i+1; j<nt; ++j) {
            if(clusterEnergy[j] < clusterEnergyMin) continue;
            if(vetoEnergy[j] > vetoEnergyMax) continue;

            const double th1 = theta[i]*deg;
            const double th2 = theta[j]*deg;
            const double ph1 = phi[i]*deg;
            const double ph2 = phi[j]*deg;

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

            const TLorentzVector candidate = g1 + g2;
            const double mgg = candidate.M();
            if(mgg <= 0.0) continue;

            const double dist = std::fabs(mgg - mpi0);
            if(dist < bestDist) {
                bestDist = dist;
                bestPi0 = candidate;
                bestMgg = mgg;
                bestOpeningDeg = g1.Angle(g2.Vect())*TMath::RadToDeg();
                found = true;
            }
        }
    }

    return found;
}

bool CalculateKinematics(const TLorentzVector& pi0,
                         double openingDeg,
                         double egamma,
                         double targetMass,
                         double mpi0,
                         double& deltaE,
                         double& deltaPhi)
{
    if(egamma <= 0.0 || openingDeg < 0.0) return false;

    TLorentzVector target(0.0, 0.0, 0.0, targetMass);
    TLorentzVector beam(0.0, 0.0, egamma, egamma);
    TLorentzVector totalCM = beam + target;

    TLorentzVector pi0CM = pi0;
    pi0CM.Boost(-totalCM.BoostVector());
    const double thetaCmDeg = pi0CM.Theta()*TMath::RadToDeg();

    const double s = targetMass*targetMass + 2.0*egamma*targetMass;
    const double sqrtS = std::sqrt(s);
    const double ePiCmExpected =
        (s + mpi0*mpi0 - targetMass*targetMass)/(2.0*sqrtS);

    // Paper Eqs. (3)-(5): both energies are evaluated in the gamma+4He CM.
    // pi0CM is the measured gamma-gamma four-vector boosted from the lab to that CM.
    const double ePiCmMeasured = pi0CM.E();
    deltaE = ePiCmMeasured - ePiCmExpected;

    const double ePiLabCoh =
        CoherentPionLabEnergy(egamma, thetaCmDeg, mpi0, targetMass);
    const double phiMin = PhiMinDeg(ePiLabCoh, mpi0);
    deltaPhi = openingDeg - phiMin;

    return true;
}

void FillData(const char* dataFile,
              TH1D* before[],
              TH1D* after[],
              Long64_t maxEvents,
              double clusterEnergyMin,
              double vetoEnergyMax,
              double mggMin,
              double mggMax,
              double deltaPhiCut,
              double promptMin,
              double promptMax,
              double randomMin,
              double randomMax)
{
    const double targetMass = 3727.38;
    const double mpi0 = 134.9768;
    const double randomWeight = (promptMax-promptMin)/(randomMax-randomMin);

    TFile* f = TFile::Open(dataFile);
    if(!f || f->IsZombie()) {
        std::cout << "Cannot open experimental file " << dataFile << std::endl;
        return;
    }

    TTree* tracks = (TTree*)f->Get("tracks");
    TTree* tagger = (TTree*)f->Get("tagger");
    TTree* setup = (TTree*)f->Get("setupParameters");
    if(!tracks || !tagger || !setup) {
        std::cout << "Experimental file is missing tracks, tagger, or setupParameters" << std::endl;
        f->Close();
        return;
    }

    int nTracks = 0;
    double clusterEnergy[kMaxTracks];
    double theta[kMaxTracks];
    double phi[kMaxTracks];
    double vetoEnergy[kMaxTracks];
    tracks->SetBranchAddress("nTracks", &nTracks);
    tracks->SetBranchAddress("clusterEnergy", clusterEnergy);
    tracks->SetBranchAddress("theta", theta);
    tracks->SetBranchAddress("phi", phi);
    tracks->SetBranchAddress("vetoEnergy", vetoEnergy);

    int nTagged = 0;
    int taggedChannel[kMaxTagged];
    double taggedTime[kMaxTagged];
    tagger->SetBranchAddress("nTagged", &nTagged);
    tagger->SetBranchAddress("taggedChannel", taggedChannel);
    tagger->SetBranchAddress("taggedTime", taggedTime);

    int nTagger = 0;
    double taggerPhotonEnergy[kNCh];
    setup->SetBranchAddress("nTagger", &nTagger);
    setup->SetBranchAddress("TaggerPhotonEnergy", taggerPhotonEnergy);
    setup->GetEntry(0);
    const int nValidCh = std::min(nTagger, kNCh);

    Long64_t nentries = std::min(tracks->GetEntries(), tagger->GetEntries());
    if(maxEvents > 0 && maxEvents < nentries) nentries = maxEvents;

    Long64_t selectedPairs = 0;
    Long64_t promptTags = 0;
    Long64_t randomTags = 0;

    for(Long64_t ev=0; ev<nentries; ++ev) {
        tracks->GetEntry(ev);
        tagger->GetEntry(ev);

        TLorentzVector pi0;
        double mgg = -1.0;
        double openingDeg = -1.0;
        if(!BuildBestPi0(nTracks, clusterEnergy, theta, phi, vetoEnergy,
                         clusterEnergyMin, vetoEnergyMax, mpi0,
                         pi0, mgg, openingDeg)) continue;
        if(mgg < mggMin || mgg > mggMax) continue;
        selectedPairs++;

        const int ntags = std::min(nTagged, kMaxTagged);
        for(int it=0; it<ntags; ++it) {
            const int ch = taggedChannel[it];
            if(ch < 0 || ch >= nValidCh) continue;
            const double egamma = taggerPhotonEnergy[ch];
            const int ib = FindEnergyBin(egamma);
            if(ib < 0) continue;

            double weight = 0.0;
            const double tt = taggedTime[it];
            if(tt > promptMin && tt < promptMax) {
                weight = 1.0;
                promptTags++;
            } else if(tt > randomMin && tt < randomMax) {
                weight = -randomWeight;
                randomTags++;
            } else {
                continue;
            }

            double deltaE = 0.0;
            double deltaPhi = 0.0;
            if(!CalculateKinematics(pi0, openingDeg, egamma,
                                    targetMass, mpi0, deltaE, deltaPhi)) continue;

            before[ib]->Fill(deltaE, weight);
            if(deltaPhi < deltaPhiCut) after[ib]->Fill(deltaE, weight);
        }

        if(ev % 500000 == 0) {
            std::cout << "data event " << ev << "/" << nentries << std::endl;
        }
    }

    std::cout << "Experimental data: processed " << nentries
              << ", selected pi0 pairs " << selectedPairs
              << ", prompt tags " << promptTags
              << ", random tags " << randomTags
              << ", random weight " << randomWeight << std::endl;

    f->Close();
}

void FillCoherentMC(const char* mcFile,
                    TH1D* before[],
                    TH1D* after[],
                    Long64_t maxEvents,
                    double clusterEnergyMin,
                    double vetoEnergyMax,
                    double mggMin,
                    double mggMax,
                    double deltaPhiCut)
{
    const double targetMass = 3727.38;
    const double mpi0 = 134.9768;

    TFile* f = TFile::Open(mcFile);
    if(!f || f->IsZombie()) {
        std::cout << "Cannot open coherent MC file " << mcFile << std::endl;
        return;
    }

    TTree* tracks = (TTree*)f->Get("tracks");
    TTree* tagger = (TTree*)f->Get("tagger");
    TTree* setup = (TTree*)f->Get("setupParameters");
    if(!tracks || !tagger || !setup) {
        std::cout << "Coherent MC file is missing tracks, tagger, or setupParameters" << std::endl;
        f->Close();
        return;
    }

    int nTracks = 0;
    double clusterEnergy[kMaxTracks];
    double theta[kMaxTracks];
    double phi[kMaxTracks];
    double vetoEnergy[kMaxTracks];
    tracks->SetBranchAddress("nTracks", &nTracks);
    tracks->SetBranchAddress("clusterEnergy", clusterEnergy);
    tracks->SetBranchAddress("theta", theta);
    tracks->SetBranchAddress("phi", phi);
    tracks->SetBranchAddress("vetoEnergy", vetoEnergy);

    int nTagged = 0;
    int taggedChannel[kMaxTagged];
    tagger->SetBranchAddress("nTagged", &nTagged);
    tagger->SetBranchAddress("taggedChannel", taggedChannel);

    int nTagger = 0;
    double taggerPhotonEnergy[kNCh];
    setup->SetBranchAddress("nTagger", &nTagger);
    setup->SetBranchAddress("TaggerPhotonEnergy", taggerPhotonEnergy);
    setup->GetEntry(0);
    const int nValidCh = std::min(nTagger, kNCh);

    Long64_t nentries = std::min(tracks->GetEntries(), tagger->GetEntries());
    if(maxEvents > 0 && maxEvents < nentries) nentries = maxEvents;

    Long64_t selectedPairs = 0;
    Long64_t tagCombinations = 0;
    Long64_t tagCombinationsAfter = 0;

    for(Long64_t ev=0; ev<nentries; ++ev) {
        tracks->GetEntry(ev);
        tagger->GetEntry(ev);

        TLorentzVector pi0;
        double mgg = -1.0;
        double openingDeg = -1.0;
        if(!BuildBestPi0(nTracks, clusterEnergy, theta, phi, vetoEnergy,
                         clusterEnergyMin, vetoEnergyMax, mpi0,
                         pi0, mgg, openingDeg)) continue;
        if(mgg < mggMin || mgg > mggMax) continue;
        selectedPairs++;

        const int ntags = std::min(nTagged, kMaxTagged);
        for(int it=0; it<ntags; ++it) {
            const int ch = taggedChannel[it];
            if(ch < 0 || ch >= nValidCh) continue;
            const double egamma = taggerPhotonEnergy[ch];
            const int ib = FindEnergyBin(egamma);
            if(ib < 0) continue;

            double deltaE = 0.0;
            double deltaPhi = 0.0;
            if(!CalculateKinematics(pi0, openingDeg, egamma,
                                    targetMass, mpi0, deltaE, deltaPhi)) continue;

            tagCombinations++;
            before[ib]->Fill(deltaE);
            if(deltaPhi < deltaPhiCut) {
                after[ib]->Fill(deltaE);
                tagCombinationsAfter++;
            }
        }

        if(ev % 500000 == 0) {
            std::cout << "MC event " << ev << "/" << nentries << std::endl;
        }
    }

    std::cout << "Coherent MC: processed " << nentries
              << ", selected pi0 pairs " << selectedPairs
              << ", selected tag combinations " << tagCombinations
              << ", after DeltaPhi cut " << tagCombinationsAfter << std::endl;

    f->Close();
}

} // namespace PaperFig2

void paper_fig2_missing_energy_root6(
    const char* dataFile="Acqu_CBTagg_31837.root",
    const char* coherentMCFile="Acqu_geant_He4pi0.root",
    Long64_t maxDataEvents=-1,
    Long64_t maxMCEvents=-1,
    double deltaPhiCut=10.0,
    double clusterEnergyMin=20.0,
    double vetoEnergyMax=1.0,
    double mggMin=110.0,
    double mggMax=155.0,
    double promptMin=700.0,
    double promptMax=800.0,
    double randomMin=450.0,
    double randomMax=680.0,
    double plotMin=-60.0,
    double plotMax=40.0)
{
    using namespace PaperFig2;

    gStyle->SetOptStat(0);

    TH1D* dataBefore[kNBins];
    TH1D* dataAfter[kNBins];
    TH1D* mcBefore[kNBins];
    TH1D* mcAfter[kNBins];

    for(int i=0; i<kNBins; ++i) {
        dataBefore[i] = new TH1D(Form("data_deltaE_before_phi_E%d", i),
            Form("Data %.0f-%.0f MeV;#Delta E_{#pi^{0}}^{*} (MeV);Prompt-random counts",
                 kELow[i], kEHigh[i]),
            240, -120.0, 120.0);
        dataAfter[i] = new TH1D(Form("data_deltaE_after_phi_E%d", i),
            Form("Data after #Delta#Phi < %.1f deg %.0f-%.0f MeV;#Delta E_{#pi^{0}}^{*} (MeV);Prompt-random counts",
                 deltaPhiCut, kELow[i], kEHigh[i]),
            240, -120.0, 120.0);
        mcBefore[i] = new TH1D(Form("coh_mc_deltaE_before_phi_E%d", i),
            Form("Coherent MC %.0f-%.0f MeV;#Delta E_{#pi^{0}}^{*} (MeV);Counts",
                 kELow[i], kEHigh[i]),
            240, -120.0, 120.0);
        mcAfter[i] = new TH1D(Form("coh_mc_deltaE_after_phi_E%d", i),
            Form("Coherent MC after #Delta#Phi < %.1f deg %.0f-%.0f MeV;#Delta E_{#pi^{0}}^{*} (MeV);Counts",
                 deltaPhiCut, kELow[i], kEHigh[i]),
            240, -120.0, 120.0);

        dataBefore[i]->Sumw2();
        dataAfter[i]->Sumw2();
        mcBefore[i]->Sumw2();
        mcAfter[i]->Sumw2();
    }

    FillData(dataFile, dataBefore, dataAfter, maxDataEvents,
             clusterEnergyMin, vetoEnergyMax, mggMin, mggMax,
             deltaPhiCut, promptMin, promptMax, randomMin, randomMax);

    FillCoherentMC(coherentMCFile, mcBefore, mcAfter, maxMCEvents,
                   clusterEnergyMin, vetoEnergyMax, mggMin, mggMax,
                   deltaPhiCut);

    std::ofstream out("paper_fig2_missing_energy.txt");
    out << "# Reproduction of the missing-energy validation used in paper Fig. 2\n";
    out << "# DeltaE_pi0^* = E_pi0^*(gamma1,gamma2) - E_pi0,coh^*(Egamma); star denotes the gamma+4He CM\n";
    out << "# Phi_min uses the coherent pion energy in the laboratory, as in paper Eq. (2)\n";
    out << "# data = " << dataFile << "\n";
    out << "# coherent MC = " << coherentMCFile << "\n";
    out << "# cuts: Ecluster > " << clusterEnergyMin
        << " MeV, veto < " << vetoEnergyMax
        << ", " << mggMin << " < mgg < " << mggMax
        << " MeV, DeltaPhi < " << deltaPhiCut << " deg\n";
    out << "# No missing-mass cut is applied.\n";
    out << "# E_low E_high data_before data_after data_after_over_before mc_before mc_after mc_eff scale_mc_to_data_after\n";

    TCanvas* c = new TCanvas("c_paper_fig2_missing_energy",
                             "Paper Fig. 2 missing-energy validation",
                             1200, 900);
    c->Divide(2,2,0.002,0.002);

    for(int i=0; i<kNBins; ++i) {
        const int bminData = dataAfter[i]->FindBin(plotMin + 1e-6);
        const int bmaxData = dataAfter[i]->FindBin(plotMax - 1e-6);
        const int bminMC = mcAfter[i]->FindBin(plotMin + 1e-6);
        const int bmaxMC = mcAfter[i]->FindBin(plotMax - 1e-6);

        const double dataBeforeInt = dataBefore[i]->Integral(0, dataBefore[i]->GetNbinsX()+1);
        const double dataAfterInt = dataAfter[i]->Integral(0, dataAfter[i]->GetNbinsX()+1);
        const double mcBeforeInt = mcBefore[i]->Integral(0, mcBefore[i]->GetNbinsX()+1);
        const double mcAfterInt = mcAfter[i]->Integral(0, mcAfter[i]->GetNbinsX()+1);

        const double dataAfterPlot = dataAfter[i]->Integral(bminData, bmaxData);
        const double mcAfterPlot = mcAfter[i]->Integral(bminMC, bmaxMC);
        const double scale = (mcAfterPlot > 0.0) ? dataAfterPlot/mcAfterPlot : 0.0;

        TH1D* mcScaled = (TH1D*)mcAfter[i]->Clone(Form("coh_mc_scaled_E%d", i));
        if(scale > 0.0) mcScaled->Scale(scale);

        out << kELow[i] << " " << kEHigh[i] << " "
            << dataBeforeInt << " " << dataAfterInt << " "
            << ((dataBeforeInt != 0.0) ? dataAfterInt/dataBeforeInt : 0.0) << " "
            << mcBeforeInt << " " << mcAfterInt << " "
            << ((mcBeforeInt > 0.0) ? mcAfterInt/mcBeforeInt : 0.0) << " "
            << scale << "\n";

        c->cd(i+1);
        gPad->SetTicks(1,1);
        gPad->SetLeftMargin(0.13);
        gPad->SetBottomMargin(0.13);

        dataBefore[i]->SetLineColor(kGray+2);
        dataBefore[i]->SetLineWidth(3);
        dataBefore[i]->SetMarkerColor(kGray+2);
        dataBefore[i]->SetMarkerStyle(20);
        dataBefore[i]->SetMarkerSize(0.45);

        dataAfter[i]->SetLineColor(kBlack);
        dataAfter[i]->SetLineWidth(3);
        dataAfter[i]->SetMarkerColor(kBlack);
        dataAfter[i]->SetMarkerStyle(20);
        dataAfter[i]->SetMarkerSize(0.45);

        mcScaled->SetLineColor(kRed+1);
        mcScaled->SetLineWidth(2);
        mcScaled->SetFillStyle(0);

        dataBefore[i]->GetXaxis()->SetRangeUser(plotMin, plotMax);
        dataBefore[i]->GetXaxis()->SetTitle("#Delta E_{#pi^{0}}^{*} (MeV)");
        dataBefore[i]->GetYaxis()->SetTitle("Counts");
        dataBefore[i]->SetTitle(kLabels[i]);

        double ymax = std::max(dataBefore[i]->GetMaximum(), dataAfter[i]->GetMaximum());
        ymax = std::max(ymax, mcScaled->GetMaximum());
        if(ymax <= 0.0) ymax = 1.0;
        dataBefore[i]->SetMinimum(0.0);
        dataBefore[i]->SetMaximum(1.18*ymax);

        dataBefore[i]->Draw("HIST");
        dataAfter[i]->Draw("HIST SAME");
        mcScaled->Draw("HIST SAME");

        TLegend* leg = new TLegend(0.48,0.68,0.89,0.89);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.037);
        leg->AddEntry(dataBefore[i], "Data after m_{#gamma#gamma}", "l");
        leg->AddEntry(dataAfter[i], Form("Data after #Delta#Phi < %.0f^{#circ}", deltaPhiCut), "l");
        leg->AddEntry(mcScaled, "Coherent MC, shape scaled", "l");
        leg->Draw();

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.038);
        latex.DrawLatex(0.16,0.84,Form("%.0f < E_{#gamma} < %.0f MeV", kELow[i], kEHigh[i]));
    }

    c->Print("paper_fig2_missing_energy.pdf");

    TFile* fout = TFile::Open("paper_fig2_missing_energy.root", "RECREATE");
    for(int i=0; i<kNBins; ++i) {
        dataBefore[i]->Write();
        dataAfter[i]->Write();
        mcBefore[i]->Write();
        mcAfter[i]->Write();
    }
    c->Write();
    fout->Close();
    out.close();

    std::cout << "Saved:\n"
              << "  paper_fig2_missing_energy.pdf\n"
              << "  paper_fig2_missing_energy.root\n"
              << "  paper_fig2_missing_energy.txt\n"
              << "No missing-mass cut was applied.\n";
}
