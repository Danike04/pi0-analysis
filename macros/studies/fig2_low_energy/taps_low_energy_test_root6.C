/*
 * taps_low_energy_test_root6.C
 *
 * PURPOSE
 *   Specialized low-energy Fig.2/TAPS cross-check for the two photon-energy
 *   intervals 201-210 MeV and 211-222 MeV. It compares data and coherent MC in
 *   DeltaE_pi0 versus reconstructed pion angle and tests DeltaPhi cuts of
 *   8, 10 and 12 degrees.
 *
 * WHAT THE MACRO ACTUALLY DOES
 *   - Reconstructs pi0 -> gamma gamma from Acqu tracks.
 *   - Uses the gamma+4He centre-of-mass missing-energy definition.
 *   - Performs prompt-random subtraction for data.
 *   - Fills 2D DeltaE_pi0 vs theta distributions for data and coherent MC.
 *   - For each DeltaPhi cut, scales coherent MC in a defined core region and
 *     evaluates residual/tail behaviour with numerical uncertainties.
 *
 * DEFAULT INPUTS
 *   dataFile       = Acqu_CBTagg_31837.root
 *   coherentMCFile = Acqu_geant_He4pi0.root
 *
 * OUTPUTS (fixed historical names)
 *   taps_low_energy_test.root
 *   taps_low_energy_test.pdf
 *   taps_low_energy_test.txt
 *
 * USAGE EXAMPLE
 *   root -l -b -q 'taps_low_energy_test_root6.C("DATA.root","MC.root")'
 *
 * REPOSITORY STATUS
 *   Preserved as a specialized low-energy study under studies/fig2_low_energy.
 *   It is not required by the main Fig.2 workflow. Code below is unchanged.
 */
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"
#include "TVirtualPad.h"

namespace TapsLowE {

const int kNCh = 352;
const int kMaxTracks = 512;
const int kMaxTagged = 2048;
const int kNE = 2;
const int kNCuts = 3;

const double kELow[kNE]  = {201.0, 211.0};
const double kEHigh[kNE] = {210.0, 222.0};
const char* kLabels[kNE] = {"201-210 MeV", "211-222 MeV"};
const double kCuts[kNCuts] = {8.0, 10.0, 12.0};

int FindEnergyBin(double egamma)
{
    for(int i=0; i<kNE; ++i) {
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

    TLorentzVector target(0.0,0.0,0.0,targetMass);
    TLorentzVector beam(0.0,0.0,egamma,egamma);
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
        CoherentPionLabEnergy(egamma,thetaCmDeg,mpi0,targetMass);
    deltaPhi = openingDeg - PhiMinDeg(ePiLabCoh,mpi0);
    return true;
}

void FillData2D(const char* dataFile,
                TH2D* h2[],
                Long64_t maxEvents,
                double clusterEnergyMin,
                double vetoEnergyMax,
                double mggMin,
                double mggMax,
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
    tracks->SetBranchAddress("nTracks",&nTracks);
    tracks->SetBranchAddress("clusterEnergy",clusterEnergy);
    tracks->SetBranchAddress("theta",theta);
    tracks->SetBranchAddress("phi",phi);
    tracks->SetBranchAddress("vetoEnergy",vetoEnergy);

    int nTagged = 0;
    int taggedChannel[kMaxTagged];
    double taggedTime[kMaxTagged];
    tagger->SetBranchAddress("nTagged",&nTagged);
    tagger->SetBranchAddress("taggedChannel",taggedChannel);
    tagger->SetBranchAddress("taggedTime",taggedTime);

    int nTagger = 0;
    double taggerPhotonEnergy[kNCh];
    setup->SetBranchAddress("nTagger",&nTagger);
    setup->SetBranchAddress("TaggerPhotonEnergy",taggerPhotonEnergy);
    setup->GetEntry(0);
    const int nValidCh = std::min(nTagger,kNCh);

    Long64_t nentries = std::min(tracks->GetEntries(),tagger->GetEntries());
    if(maxEvents > 0 && maxEvents < nentries) nentries = maxEvents;

    for(Long64_t ev=0; ev<nentries; ++ev) {
        tracks->GetEntry(ev);
        tagger->GetEntry(ev);

        TLorentzVector pi0;
        double mgg = -1.0;
        double openingDeg = -1.0;
        if(!BuildBestPi0(nTracks,clusterEnergy,theta,phi,vetoEnergy,
                         clusterEnergyMin,vetoEnergyMax,mpi0,
                         pi0,mgg,openingDeg)) continue;
        if(mgg < mggMin || mgg > mggMax) continue;

        const int ntags = std::min(nTagged,kMaxTagged);
        for(int it=0; it<ntags; ++it) {
            const int ch = taggedChannel[it];
            if(ch < 0 || ch >= nValidCh) continue;
            const double egamma = taggerPhotonEnergy[ch];
            const int ie = FindEnergyBin(egamma);
            if(ie < 0) continue;

            double weight = 0.0;
            const double tt = taggedTime[it];
            if(tt > promptMin && tt < promptMax) weight = 1.0;
            else if(tt > randomMin && tt < randomMax) weight = -randomWeight;
            else continue;

            double deltaE = 0.0;
            double deltaPhi = 0.0;
            if(!CalculateKinematics(pi0,openingDeg,egamma,targetMass,mpi0,
                                    deltaE,deltaPhi)) continue;
            h2[ie]->Fill(deltaPhi,deltaE,weight);
        }
        if(ev % 500000 == 0) std::cout << "data event " << ev << "/" << nentries << std::endl;
    }
    f->Close();
}

void FillMC2D(const char* mcFile,
              TH2D* h2[],
              Long64_t maxEvents,
              double clusterEnergyMin,
              double vetoEnergyMax,
              double mggMin,
              double mggMax)
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
    tracks->SetBranchAddress("nTracks",&nTracks);
    tracks->SetBranchAddress("clusterEnergy",clusterEnergy);
    tracks->SetBranchAddress("theta",theta);
    tracks->SetBranchAddress("phi",phi);
    tracks->SetBranchAddress("vetoEnergy",vetoEnergy);

    int nTagged = 0;
    int taggedChannel[kMaxTagged];
    tagger->SetBranchAddress("nTagged",&nTagged);
    tagger->SetBranchAddress("taggedChannel",taggedChannel);

    int nTagger = 0;
    double taggerPhotonEnergy[kNCh];
    setup->SetBranchAddress("nTagger",&nTagger);
    setup->SetBranchAddress("TaggerPhotonEnergy",taggerPhotonEnergy);
    setup->GetEntry(0);
    const int nValidCh = std::min(nTagger,kNCh);

    Long64_t nentries = std::min(tracks->GetEntries(),tagger->GetEntries());
    if(maxEvents > 0 && maxEvents < nentries) nentries = maxEvents;

    for(Long64_t ev=0; ev<nentries; ++ev) {
        tracks->GetEntry(ev);
        tagger->GetEntry(ev);

        TLorentzVector pi0;
        double mgg = -1.0;
        double openingDeg = -1.0;
        if(!BuildBestPi0(nTracks,clusterEnergy,theta,phi,vetoEnergy,
                         clusterEnergyMin,vetoEnergyMax,mpi0,
                         pi0,mgg,openingDeg)) continue;
        if(mgg < mggMin || mgg > mggMax) continue;

        const int ntags = std::min(nTagged,kMaxTagged);
        for(int it=0; it<ntags; ++it) {
            const int ch = taggedChannel[it];
            if(ch < 0 || ch >= nValidCh) continue;
            const double egamma = taggerPhotonEnergy[ch];
            const int ie = FindEnergyBin(egamma);
            if(ie < 0) continue;

            double deltaE = 0.0;
            double deltaPhi = 0.0;
            if(!CalculateKinematics(pi0,openingDeg,egamma,targetMass,mpi0,
                                    deltaE,deltaPhi)) continue;
            h2[ie]->Fill(deltaPhi,deltaE);
        }
        if(ev % 500000 == 0) std::cout << "MC event " << ev << "/" << nentries << std::endl;
    }
    f->Close();
}

struct Metrics {
    double scale;
    double dataCore;
    double mcCore;
    double dataTail;
    double mcTail;
    double excess;
    double excessErr;
    double chi2;
    int ndf;
    double pvalue;
};

double IntegralAndErrorRange(TH1D* h,double xmin,double xmax,double& err)
{
    const int b1 = h->FindBin(xmin+1e-9);
    const int b2 = h->FindBin(xmax-1e-9);
    return h->IntegralAndError(b1,b2,err);
}

Metrics Evaluate(TH1D* data,TH1D* mc,
                 double coreMin,double coreMax,
                 double tailMin,double tailMax,
                 double compareMin,double compareMax)
{
    Metrics m;
    m.scale = m.dataCore = m.mcCore = 0.0;
    m.dataTail = m.mcTail = m.excess = m.excessErr = 0.0;
    m.chi2 = 0.0;
    m.ndf = 0;
    m.pvalue = 0.0;

    double edc=0.0, emc=0.0;
    m.dataCore = IntegralAndErrorRange(data,coreMin,coreMax,edc);
    m.mcCore = IntegralAndErrorRange(mc,coreMin,coreMax,emc);
    if(m.mcCore > 0.0) m.scale = m.dataCore/m.mcCore;

    double edt=0.0, emtRaw=0.0;
    m.dataTail = IntegralAndErrorRange(data,tailMin,tailMax,edt);
    const double mcTailRaw = IntegralAndErrorRange(mc,tailMin,tailMax,emtRaw);
    m.mcTail = m.scale*mcTailRaw;
    const double emt = m.scale*emtRaw;
    m.excess = m.dataTail-m.mcTail;
    m.excessErr = std::sqrt(edt*edt+emt*emt);

    const int b1 = data->FindBin(compareMin+1e-9);
    const int b2 = data->FindBin(compareMax-1e-9);
    int used=0;
    for(int b=b1; b<=b2; ++b) {
        const double d = data->GetBinContent(b);
        const double q = m.scale*mc->GetBinContent(b);
        const double ed = data->GetBinError(b);
        const double eq = m.scale*mc->GetBinError(b);
        const double var = ed*ed+eq*eq;
        if(var <= 0.0) continue;
        const double diff = d-q;
        m.chi2 += diff*diff/var;
        ++used;
    }
    m.ndf = used>1 ? used-1 : 0;
    if(m.ndf>0) m.pvalue = TMath::Prob(m.chi2,m.ndf);
    return m;
}

} // namespace TapsLowE

void taps_low_energy_test_root6(
    const char* dataFile="Acqu_CBTagg_31837.root",
    const char* coherentMCFile="Acqu_geant_He4pi0.root",
    Long64_t maxDataEvents=-1,
    Long64_t maxMCEvents=-1,
    double clusterEnergyMin=20.0,
    double vetoEnergyMax=1.0,
    double mggMin=110.0,
    double mggMax=155.0,
    double promptMin=700.0,
    double promptMax=800.0,
    double randomMin=450.0,
    double randomMax=680.0,
    double coreMin=-10.0,
    double coreMax=15.0,
    double tailMin=-60.0,
    double tailMax=-20.0,
    double compareMin=-60.0,
    double compareMax=40.0)
{
    using namespace TapsLowE;
    gStyle->SetOptStat(0);

    TH2D* data2D[kNE];
    TH2D* mc2D[kNE];
    for(int ie=0; ie<kNE; ++ie) {
        data2D[ie] = new TH2D(Form("data_deltaE_vs_deltaPhi_E%d",ie),
            Form("Data %.0f-%.0f MeV;#Delta#Phi (deg);#Delta E_{#pi^{0}}^{*} (MeV)",kELow[ie],kEHigh[ie]),
            360,-60.0,120.0,240,-120.0,120.0);
        mc2D[ie] = new TH2D(Form("coh_mc_deltaE_vs_deltaPhi_E%d",ie),
            Form("Coherent MC %.0f-%.0f MeV;#Delta#Phi (deg);#Delta E_{#pi^{0}}^{*} (MeV)",kELow[ie],kEHigh[ie]),
            360,-60.0,120.0,240,-120.0,120.0);
        data2D[ie]->Sumw2();
        mc2D[ie]->Sumw2();
    }

    FillData2D(dataFile,data2D,maxDataEvents,clusterEnergyMin,vetoEnergyMax,
               mggMin,mggMax,promptMin,promptMax,randomMin,randomMax);
    FillMC2D(coherentMCFile,mc2D,maxMCEvents,clusterEnergyMin,vetoEnergyMax,
             mggMin,mggMax);

    TFile* fout = TFile::Open("taps_low_energy_test.root","RECREATE");
    std::ofstream out("taps_low_energy_test.txt");
    out << "# Direct test of the two low-energy bins relevant to the suspected TAPS point\n";
    out << "# DeltaE_pi0^* = E_pi0^*(gamma1,gamma2) - E_pi0,coh^*(Egamma); star denotes the gamma+4He CM\n";
    out << "# Phi_min uses the coherent pion energy in the laboratory\n";
    out << "# data = " << dataFile << "\n";
    out << "# coherent MC = " << coherentMCFile << "\n";
    out << "# cuts: Ecluster > " << clusterEnergyMin << " MeV, veto < " << vetoEnergyMax
        << ", " << mggMin << " < mgg < " << mggMax << " MeV\n";
    out << "# MC normalized in DeltaE core = [" << coreMin << "," << coreMax << "] MeV\n";
    out << "# negative tail = [" << tailMin << "," << tailMax << "] MeV\n";
    out << "# chi2 range = [" << compareMin << "," << compareMax << "] MeV\n";
    out << "# cut E_low E_high data_after mc_after scale data_core mc_core_raw data_tail mc_tail_scaled tail_excess tail_excess_err tail_significance chi2 ndf chi2_ndf pvalue\n";

    TCanvas* c = new TCanvas("c_taps_low_energy","TAPS low-energy test",1200,600);
    c->Print("taps_low_energy_test.pdf[");

    for(int ic=0; ic<kNCuts; ++ic) {
        c->Clear();
        c->Divide(2,1,0.002,0.002);
        for(int ie=0; ie<kNE; ++ie) {
            const int xFirst = 1;
            const int xLast = data2D[ie]->GetXaxis()->FindBin(kCuts[ic]-1e-9);

            TH1D* data = data2D[ie]->ProjectionY(Form("data_after_cut%.0f_E%d",kCuts[ic],ie),xFirst,xLast,"e");
            TH1D* mc = mc2D[ie]->ProjectionY(Form("coh_mc_after_cut%.0f_E%d",kCuts[ic],ie),xFirst,xLast,"e");
            data->SetDirectory(0);
            mc->SetDirectory(0);

            Metrics m = Evaluate(data,mc,coreMin,coreMax,tailMin,tailMax,compareMin,compareMax);
            TH1D* mcScaled = (TH1D*)mc->Clone(Form("coh_mc_scaled_cut%.0f_E%d",kCuts[ic],ie));
            mcScaled->SetDirectory(0);
            if(m.scale > 0.0) mcScaled->Scale(m.scale);

            const double sig = m.excessErr>0.0 ? m.excess/m.excessErr : 0.0;
            const double chi2ndf = m.ndf>0 ? m.chi2/m.ndf : 0.0;
            out << kCuts[ic] << " " << kELow[ie] << " " << kEHigh[ie] << " "
                << data->Integral(0,data->GetNbinsX()+1) << " "
                << mc->Integral(0,mc->GetNbinsX()+1) << " "
                << m.scale << " " << m.dataCore << " " << m.mcCore << " "
                << m.dataTail << " " << m.mcTail << " " << m.excess << " "
                << m.excessErr << " " << sig << " " << m.chi2 << " " << m.ndf << " "
                << chi2ndf << " " << m.pvalue << "\n";

            c->cd(ie+1);
            gPad->SetTicks(1,1);
            gPad->SetLeftMargin(0.13);
            gPad->SetBottomMargin(0.13);

            data->SetTitle(Form("%s, #Delta#Phi < %.0f^{#circ}",kLabels[ie],kCuts[ic]));
            data->GetXaxis()->SetTitle("#Delta E_{#pi^{0}}^{*} (MeV)");
            data->GetYaxis()->SetTitle("Counts");
            data->GetXaxis()->SetRangeUser(compareMin,compareMax);
            data->SetLineColor(kBlack);
            data->SetMarkerColor(kBlack);
            data->SetMarkerStyle(20);
            data->SetMarkerSize(0.55);
            data->SetLineWidth(2);

            mcScaled->SetLineColor(kRed+1);
            mcScaled->SetLineWidth(2);
            mcScaled->SetFillStyle(0);

            double ymax = std::max(data->GetMaximum(),mcScaled->GetMaximum());
            data->SetMinimum(0.0);
            data->SetMaximum(ymax>0.0 ? 1.22*ymax : 1.0);
            data->Draw("E");
            mcScaled->Draw("HIST SAME");

            TLegend* leg = new TLegend(0.50,0.72,0.89,0.88);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->AddEntry(data,"Data after cut","lep");
            leg->AddEntry(mcScaled,"Coherent MC (core scaled)","l");
            leg->Draw();

            TLatex tx;
            tx.SetNDC();
            tx.SetTextSize(0.038);
            tx.DrawLatex(0.16,0.87,Form("#chi^{2}/ndf = %.2f",chi2ndf));
            tx.DrawLatex(0.16,0.81,Form("negative-tail excess = %.1f #sigma",sig));

            fout->cd();
            data->Write();
            mc->Write();
            mcScaled->Write();
        }
        c->Print("taps_low_energy_test.pdf");
    }
    c->Print("taps_low_energy_test.pdf]");

    fout->cd();
    for(int ie=0; ie<kNE; ++ie) {
        data2D[ie]->Write();
        mc2D[ie]->Write();
    }
    c->Write();
    fout->Close();
    out.close();

    std::cout << "Saved taps_low_energy_test.root, taps_low_energy_test.pdf, and taps_low_energy_test.txt" << std::endl;
}
