/*
 * fig2_empty_target_same_analysis.C
 *
 * Purpose
 * -------
 * Apply the same event reconstruction, prompt-random subtraction and
 * DeltaPhi selections used in the Fig. 2 Acqu-data studies to an EMPTY-target
 * data file.  The macro is a diagnostic of the empty-target missing-energy
 * contribution; it does not perform FULL-EMPTY photon-flux normalization.
 *
 * Historical default input: Acqu_CBTagg_32410.root.
 * Default selections: cluster energy 20 MeV, veto <= 1 MeV,
 * 110 < m(gamma gamma) < 155 MeV, prompt 700-800 ns,
 * random 450-680 ns; histograms are produced before and after 8/10/12 deg.
 *
 * Outputs
 * -------
 *   fig2_empty_<outputTag>.root/.txt
 *   fig2_empty_<outputTag>_before.pdf
 *   fig2_empty_<outputTag>_after.pdf
 *
 * Usage
 * -----
 *   root -l -b -q 'fig2_empty_target_same_analysis.C()'
 *
 * The code body below is unchanged from the supplied source.
 */
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"

namespace Fig2FromAcqu {

const int kNCh = 352;
const int kMaxTracks = 512;
const int kMaxTagged = 2048;
const int kNE = 4;
const int kNCuts = 3;

const double kELow[kNE]  = {223.0, 283.0, 319.0, 356.0};
const double kEHigh[kNE] = {234.0, 294.0, 330.0, 366.0};
const char* kLabels[kNE] = {"224 MeV", "294 MeV", "320 MeV", "366 MeV"};
const int kCuts[kNCuts] = {8, 10, 12};

struct Chi2Result {
    double chi2;
    int ndf;
    double pvalue;
};

struct TwoComponentFit {
    double ppnnScale;
    double he3nScale;
};

int FindEnergyBin(double egamma)
{
    for (int i = 0; i < kNE; ++i) {
        const bool inside =
            egamma >= kELow[i] &&
            (egamma < kEHigh[i] ||
             (i == kNE - 1 && egamma <= kEHigh[i]));
        if (inside) return i;
    }
    return -1;
}

double IntegralRange(TH1D* h, double xmin, double xmax, double* err = 0)
{
    if (err) *err = 0.0;
    if (!h) return 0.0;
    const int b1 = h->GetXaxis()->FindBin(xmin + 1.0e-9);
    const int b2 = h->GetXaxis()->FindBin(xmax - 1.0e-9);
    if (err) return h->IntegralAndError(b1, b2, *err);
    return h->Integral(b1, b2);
}

double IntegralGlobal(TH1* h)
{
    return h ? h->Integral(0, h->GetNbinsX() + 1) : 0.0;
}

TH1D* MapToReference(TH1* source, TH1D* reference, const char* name)
{
    if (!source || !reference) return 0;

    TH1D* out = dynamic_cast<TH1D*>(reference->Clone(name));
    if (!out) return 0;

    out->SetDirectory(0);
    out->Reset("ICES");
    out->Sumw2();

    const double xmin = out->GetXaxis()->GetXmin();
    const double xmax = out->GetXaxis()->GetXmax();

    for (int b = 1; b <= source->GetNbinsX(); ++b) {
        const double x = source->GetXaxis()->GetBinCenter(b);
        if (x < xmin || x >= xmax) continue;

        const int bo = out->GetXaxis()->FindBin(x);
        const double c0 = out->GetBinContent(bo);
        const double e0 = out->GetBinError(bo);
        const double c1 = source->GetBinContent(b);
        const double e1 = source->GetBinError(b);

        out->SetBinContent(bo, c0 + c1);
        out->SetBinError(bo, std::sqrt(e0*e0 + e1*e1));
    }
    return out;
}

double CoreScale(TH1D* data, TH1D* coherent,
                 double coreMin, double coreMax)
{
    const double d = IntegralRange(data, coreMin, coreMax);
    const double c = IntegralRange(coherent, coreMin, coreMax);
    return c > 0.0 ? d/c : 0.0;
}

double FitPPNNBefore(TH1D* data, TH1D* coherent,
                     TH1D* ppnn,
                     double tailMin, double tailMax)
{
    const int b1 = data->GetXaxis()->FindBin(tailMin + 1.0e-9);
    const int b2 = data->GetXaxis()->FindBin(tailMax - 1.0e-9);

    double num = 0.0;
    double den = 0.0;

    for (int b = b1; b <= b2; ++b) {
        const double p = ppnn->GetBinContent(b);
        if (p <= 0.0) continue;

        const double ed = data->GetBinError(b);
        const double ec = coherent->GetBinError(b);
        const double var = ed*ed + ec*ec;
        if (var <= 0.0) continue;

        const double w = 1.0/var;
        const double y =
            data->GetBinContent(b) -
            coherent->GetBinContent(b);

        num += w*p*y;
        den += w*p*p;
    }

    if (den <= 0.0) return 0.0;
    return std::max(0.0, num/den);
}

double TwoComponentSSE(TH1D* data, TH1D* coherent,
                       TH1D* ppnn, TH1D* he3n,
                       double pScale, double hScale,
                       double tailMin, double tailMax)
{
    const int b1 = data->GetXaxis()->FindBin(tailMin + 1.0e-9);
    const int b2 = data->GetXaxis()->FindBin(tailMax - 1.0e-9);

    double sse = 0.0;

    for (int b = b1; b <= b2; ++b) {
        const double ed = data->GetBinError(b);
        const double ec = coherent->GetBinError(b);
        const double var = ed*ed + ec*ec;
        if (var <= 0.0) continue;

        const double diff =
            data->GetBinContent(b) -
            coherent->GetBinContent(b) -
            pScale*ppnn->GetBinContent(b) -
            hScale*he3n->GetBinContent(b);

        sse += diff*diff/var;
    }
    return sse;
}

TwoComponentFit FitPPNNHe3nBefore(TH1D* data, TH1D* coherent,
                                  TH1D* ppnn, TH1D* he3n,
                                  double tailMin, double tailMax)
{
    TwoComponentFit result = {0.0, 0.0};

    const int b1 = data->GetXaxis()->FindBin(tailMin + 1.0e-9);
    const int b2 = data->GetXaxis()->FindBin(tailMax - 1.0e-9);

    double a11 = 0.0;
    double a12 = 0.0;
    double a22 = 0.0;
    double rhs1 = 0.0;
    double rhs2 = 0.0;

    for (int b = b1; b <= b2; ++b) {
        const double ed = data->GetBinError(b);
        const double ec = coherent->GetBinError(b);
        const double var = ed*ed + ec*ec;
        if (var <= 0.0) continue;

        const double w = 1.0/var;
        const double y =
            data->GetBinContent(b) -
            coherent->GetBinContent(b);

        const double p = ppnn->GetBinContent(b);
        const double h = he3n->GetBinContent(b);

        a11 += w*p*p;
        a12 += w*p*h;
        a22 += w*h*h;
        rhs1 += w*p*y;
        rhs2 += w*h*y;
    }

    double bestP = 0.0;
    double bestH = 0.0;
    double bestSSE = TwoComponentSSE(
        data, coherent, ppnn, he3n,
        bestP, bestH, tailMin, tailMax);

    if (a11 > 0.0) {
        const double p = std::max(0.0, rhs1/a11);
        const double sse = TwoComponentSSE(
            data, coherent, ppnn, he3n,
            p, 0.0, tailMin, tailMax);
        if (sse < bestSSE) {
            bestSSE = sse;
            bestP = p;
            bestH = 0.0;
        }
    }

    if (a22 > 0.0) {
        const double h = std::max(0.0, rhs2/a22);
        const double sse = TwoComponentSSE(
            data, coherent, ppnn, he3n,
            0.0, h, tailMin, tailMax);
        if (sse < bestSSE) {
            bestSSE = sse;
            bestP = 0.0;
            bestH = h;
        }
    }

    const double det = a11*a22 - a12*a12;
    if (det > 0.0) {
        const double p = (rhs1*a22 - rhs2*a12)/det;
        const double h = (rhs2*a11 - rhs1*a12)/det;

        if (p >= 0.0 && h >= 0.0) {
            const double sse = TwoComponentSSE(
                data, coherent, ppnn, he3n,
                p, h, tailMin, tailMax);

            if (sse < bestSSE) {
                bestSSE = sse;
                bestP = p;
                bestH = h;
            }
        }
    }

    result.ppnnScale = bestP;
    result.he3nScale = bestH;
    return result;
}

Chi2Result EvaluateChi2(TH1D* data, TH1D* model,
                        double xmin, double xmax,
                        int nFittedParameters)
{
    Chi2Result result = {0.0, 0, 0.0};

    const int b1 = data->GetXaxis()->FindBin(xmin + 1.0e-9);
    const int b2 = data->GetXaxis()->FindBin(xmax - 1.0e-9);

    int used = 0;

    for (int b = b1; b <= b2; ++b) {
        const double ed = data->GetBinError(b);
        const double em = model->GetBinError(b);
        const double var = ed*ed + em*em;
        if (var <= 0.0) continue;

        const double diff =
            data->GetBinContent(b) -
            model->GetBinContent(b);

        result.chi2 += diff*diff/var;
        ++used;
    }

    result.ndf = std::max(0, used - nFittedParameters);
    if (result.ndf > 0)
        result.pvalue = TMath::Prob(result.chi2, result.ndf);

    return result;
}

double TailResidualSigma(TH1D* data, TH1D* model,
                         double tailMin, double tailMax,
                         double& residual, double& residualError)
{
    double ed = 0.0;
    double em = 0.0;

    const double d = IntegralRange(data, tailMin, tailMax, &ed);
    const double m = IntegralRange(model, tailMin, tailMax, &em);

    residual = d-m;
    residualError = std::sqrt(ed*ed + em*em);

    return residualError > 0.0
        ? residual/residualError
        : 0.0;
}

double CoherentPionLabEnergy(double egamma,
                             double thetaCmDeg,
                             double mpi0,
                             double targetMass)
{
    const double s =
        targetMass*targetMass +
        2.0*egamma*targetMass;

    const double sqrtS = std::sqrt(s);

    const double ePiCm =
        (s + mpi0*mpi0 - targetMass*targetMass) /
        (2.0*sqrtS);

    const double pPiCm =
        std::sqrt(std::max(0.0, ePiCm*ePiCm - mpi0*mpi0));

    const double betaCm =
        egamma/(egamma + targetMass);

    const double gammaCm =
        1.0/std::sqrt(1.0 - betaCm*betaCm);

    return gammaCm *
        (ePiCm +
         betaCm*pPiCm*
         std::cos(thetaCmDeg*TMath::DegToRad()));
}

double PhiMinDeg(double ePiLab, double mpi0)
{
    if (ePiLab <= mpi0) return 180.0;

    double arg =
        std::sqrt(ePiLab*ePiLab - mpi0*mpi0) /
        ePiLab;

    arg = std::max(-1.0, std::min(1.0, arg));

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

    bool found = false;
    double bestDistance = 1.0e30;

    for (int i = 0; i < nt; ++i) {
        if (clusterEnergy[i] < clusterEnergyMin) continue;
        if (vetoEnergy[i] > vetoEnergyMax) continue;

        for (int j = i + 1; j < nt; ++j) {
            if (clusterEnergy[j] < clusterEnergyMin) continue;
            if (vetoEnergy[j] > vetoEnergyMax) continue;

            const double th1 = theta[i]*deg;
            const double th2 = theta[j]*deg;
            const double ph1 = phi[i]*deg;
            const double ph2 = phi[j]*deg;

            TLorentzVector g1, g2;

            g1.SetPxPyPzE(
                clusterEnergy[i]*std::sin(th1)*std::cos(ph1),
                clusterEnergy[i]*std::sin(th1)*std::sin(ph1),
                clusterEnergy[i]*std::cos(th1),
                clusterEnergy[i]);

            g2.SetPxPyPzE(
                clusterEnergy[j]*std::sin(th2)*std::cos(ph2),
                clusterEnergy[j]*std::sin(th2)*std::sin(ph2),
                clusterEnergy[j]*std::cos(th2),
                clusterEnergy[j]);

            const TLorentzVector candidate = g1 + g2;
            const double mgg = candidate.M();
            if (mgg <= 0.0) continue;

            const double distance = std::fabs(mgg - mpi0);

            if (distance < bestDistance) {
                bestDistance = distance;
                bestPi0 = candidate;
                bestMgg = mgg;
                bestOpeningDeg =
                    g1.Angle(g2.Vect())*TMath::RadToDeg();
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
    TLorentzVector target(0.0, 0.0, 0.0, targetMass);
    TLorentzVector beam(0.0, 0.0, egamma, egamma);
    TLorentzVector totalCM = beam + target;

    TLorentzVector pi0CM = pi0;
    pi0CM.Boost(-totalCM.BoostVector());

    const double thetaCmDeg =
        pi0CM.Theta()*TMath::RadToDeg();

    const double s =
        targetMass*targetMass +
        2.0*egamma*targetMass;

    const double sqrtS = std::sqrt(s);

    const double ePiCmExpected =
        (s + mpi0*mpi0 - targetMass*targetMass) /
        (2.0*sqrtS);

    const double betaCm =
        egamma/(egamma + targetMass);

    const double gammaCm =
        1.0/std::sqrt(1.0 - betaCm*betaCm);

    const double ePiCmMeasured =
        gammaCm*(pi0.E() - betaCm*pi0.Pz());

    deltaE = ePiCmMeasured - ePiCmExpected;

    const double coherentLabEnergy =
        CoherentPionLabEnergy(
            egamma, thetaCmDeg, mpi0, targetMass);

    deltaPhi =
        openingDeg -
        PhiMinDeg(coherentLabEnergy, mpi0);

    return true;
}

bool FillExperimentalData(
    const char* dataFile,
    TH1D* before[kNE],
    TH1D* after[kNCuts][kNE],
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
    const double randomWeight =
        (promptMax-promptMin)/(randomMax-randomMin);

    TFile* f = TFile::Open(dataFile);
    if (!f || f->IsZombie()) {
        std::cerr << "Cannot open " << dataFile << std::endl;
        return false;
    }

    TTree* tracks =
        dynamic_cast<TTree*>(f->Get("tracks"));
    TTree* tagger =
        dynamic_cast<TTree*>(f->Get("tagger"));
    TTree* setup =
        dynamic_cast<TTree*>(f->Get("setupParameters"));

    if (!tracks || !tagger || !setup) {
        std::cerr << "Missing trees in " << dataFile << std::endl;
        f->Close();
        return false;
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
    setup->SetBranchAddress(
        "TaggerPhotonEnergy",
        taggerPhotonEnergy);
    setup->GetEntry(0);

    const int nValidChannels =
        std::min(nTagger, kNCh);

    Long64_t entries =
        std::min(tracks->GetEntries(),
                 tagger->GetEntries());

    if (maxEvents > 0 && maxEvents < entries)
        entries = maxEvents;

    for (Long64_t event = 0; event < entries; ++event) {
        tracks->GetEntry(event);
        tagger->GetEntry(event);

        TLorentzVector pi0;
        double mgg = -1.0;
        double openingDeg = -1.0;

        if (!BuildBestPi0(
                nTracks,
                clusterEnergy,
                theta,
                phi,
                vetoEnergy,
                clusterEnergyMin,
                vetoEnergyMax,
                mpi0,
                pi0,
                mgg,
                openingDeg))
            continue;

        if (mgg < mggMin || mgg > mggMax)
            continue;

        const int nTags =
            std::min(nTagged, kMaxTagged);

        for (int iTag = 0; iTag < nTags; ++iTag) {
            const int channel = taggedChannel[iTag];
            if (channel < 0 || channel >= nValidChannels)
                continue;

            const double egamma =
                taggerPhotonEnergy[channel];

            const int energyBin =
                FindEnergyBin(egamma);

            if (energyBin < 0)
                continue;

            double weight = 0.0;
            const double time = taggedTime[iTag];

            if (time > promptMin && time < promptMax)
                weight = 1.0;
            else if (time > randomMin && time < randomMax)
                weight = -randomWeight;
            else
                continue;

            double deltaE = 0.0;
            double deltaPhi = 0.0;

            if (!CalculateKinematics(
                    pi0, openingDeg, egamma,
                    targetMass, mpi0,
                    deltaE, deltaPhi))
                continue;

            before[energyBin]->Fill(deltaE, weight);

            for (int ic = 0; ic < kNCuts; ++ic) {
                if (deltaPhi < kCuts[ic])
                    after[ic][energyBin]->Fill(
                        deltaE, weight);
            }
        }

        if ((event + 1) % 500000 == 0)
            std::cout << "Event "
                      << event + 1 << " / "
                      << entries << std::endl;
    }

    f->Close();
    return true;
}

void StyleData(TH1D* h)
{
    h->SetLineColor(kBlack);
    h->SetMarkerColor(kBlack);
    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.45);
    h->SetLineWidth(2);
}

void StyleCoherent(TH1D* h)
{
    h->SetLineColor(kRed + 1);
    h->SetLineWidth(2);
}

void StylePPNN(TH1D* h)
{
    h->SetLineColor(kBlue + 1);
    h->SetLineWidth(2);
    h->SetLineStyle(2);
}

void StyleHe3n(TH1D* h)
{
    h->SetLineColor(kMagenta + 2);
    h->SetLineWidth(2);
    h->SetLineStyle(7);
}

void StylePPNNTotal(TH1D* h)
{
    h->SetLineColor(kOrange + 7);
    h->SetLineWidth(2);
    h->SetLineStyle(9);
}

void StyleTotal(TH1D* h)
{
    h->SetLineColor(kGreen + 2);
    h->SetLineWidth(3);
}

} // namespace Fig2FromAcqu


void fig2_empty_target_same_analysis(
    const char* dataFile =
        "../Acqu_CBTagg_32410.root",
    const char* outputTag =
        "empty_32410_32414",
    Long64_t maxEvents = -1,
    double clusterEnergyMin = 20.0,
    double vetoEnergyMax = 1.0,
    double mggMin = 110.0,
    double mggMax = 155.0,
    double promptMin = 700.0,
    double promptMax = 800.0,
    double randomMin = 450.0,
    double randomMax = 680.0)
{
    using namespace Fig2FromAcqu;

    gStyle->SetOptStat(0);

    TString tag(outputTag);
    tag.ReplaceAll("/", "_");
    tag.ReplaceAll(" ", "_");

    TH1D* dataBefore[kNE];
    TH1D* dataAfter[kNCuts][kNE];

    for (int ie = 0; ie < kNE; ++ie) {
        dataBefore[ie] = new TH1D(
            Form("empty_data_before_E%d", ie),
            "", 240, -120.0, 120.0);
        dataBefore[ie]->SetDirectory(0);
        dataBefore[ie]->Sumw2();

        for (int ic = 0; ic < kNCuts; ++ic) {
            dataAfter[ic][ie] = new TH1D(
                Form("empty_data_after_cut%d_E%d",
                     kCuts[ic], ie),
                "", 240, -120.0, 120.0);
            dataAfter[ic][ie]->SetDirectory(0);
            dataAfter[ic][ie]->Sumw2();
        }
    }

    if (!FillExperimentalData(
            dataFile,
            dataBefore,
            dataAfter,
            maxEvents,
            clusterEnergyMin,
            vetoEnergyMax,
            mggMin,
            mggMax,
            promptMin,
            promptMax,
            randomMin,
            randomMax))
        return;

    const TString beforePdf =
        Form("fig2_empty_%s_before.pdf", tag.Data());
    const TString afterPdf =
        Form("fig2_empty_%s_after.pdf", tag.Data());
    const TString rootName =
        Form("fig2_empty_%s.root", tag.Data());
    const TString textName =
        Form("fig2_empty_%s.txt", tag.Data());

    TFile* output =
        TFile::Open(rootName, "RECREATE");

    if (!output || output->IsZombie()) {
        std::cerr << "Cannot create " << rootName << std::endl;
        return;
    }

    std::ofstream txt(textName.Data());
    txt << "# Empty-target analysis with the same event selection as full target\n";
    txt << "# Data file: " << dataFile << "\n";
    txt << "# energy_index energy_label panel integral_-60_40 error_-60_40 "
           "integral_all\n";

    auto prepareHistogram = [](TH1D* h,
                               const char* title) {
        if (!h) return;
        StyleData(h);
        h->SetTitle(title);
        h->GetXaxis()->SetTitle("#Delta E_{#pi^{0}}^{*} (MeV)");
        h->GetYaxis()->SetTitle("Counts");
        h->GetXaxis()->SetRangeUser(-60.0, 40.0);
        h->GetYaxis()->SetTitleOffset(1.25);

        double maximum = h->GetMaximum();
        double minimum = h->GetMinimum();
        if (maximum <= 0.0) maximum = 1.0;

        h->SetMaximum(1.18*maximum);
        if (minimum < 0.0)
            h->SetMinimum(1.18*minimum);
        else
            h->SetMinimum(0.0);
    };

    TCanvas* beforeCanvas =
        new TCanvas("c_empty_before",
                    "empty target before DeltaPhi cut",
                    1300, 900);
    beforeCanvas->Divide(2, 2, 0.002, 0.002);

    for (int ie = 0; ie < kNE; ++ie) {
        beforeCanvas->cd(ie + 1);
        gPad->SetTicks(1, 1);
        gPad->SetLeftMargin(0.12);
        gPad->SetBottomMargin(0.12);

        prepareHistogram(
            dataBefore[ie],
            Form("Empty target: %s, before #Delta#Phi cut",
                 kLabels[ie]));

        dataBefore[ie]->Draw("E1");

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.040);
        latex.DrawLatex(0.16, 0.84, "Empty target");
        latex.SetTextSize(0.032);
        latex.DrawLatex(
            0.16, 0.78,
            Form("Entries = %.0f",
                 dataBefore[ie]->GetEntries()));

        double error = 0.0;
        const double integral =
            IntegralRange(dataBefore[ie], -60.0, 40.0,
                          &error);

        txt << ie << " \"" << kLabels[ie]
            << "\" before "
            << integral << " "
            << error << " "
            << IntegralGlobal(dataBefore[ie])
            << "\n";

        output->cd();
        dataBefore[ie]->Write();
    }

    beforeCanvas->SaveAs(beforePdf);

    TCanvas* afterCanvas =
        new TCanvas("c_empty_after",
                    "empty target after DeltaPhi cuts",
                    1300, 900);

    afterCanvas->Print(afterPdf + "[");

    for (int ic = 0; ic < kNCuts; ++ic) {
        afterCanvas->Clear();
        afterCanvas->Divide(2, 2, 0.002, 0.002);

        for (int ie = 0; ie < kNE; ++ie) {
            afterCanvas->cd(ie + 1);
            gPad->SetTicks(1, 1);
            gPad->SetLeftMargin(0.12);
            gPad->SetBottomMargin(0.12);

            prepareHistogram(
                dataAfter[ic][ie],
                Form("Empty target: %s, #Delta#Phi < %d^{#circ}",
                     kLabels[ie], kCuts[ic]));

            dataAfter[ic][ie]->Draw("E1");

            TLatex latex;
            latex.SetNDC();
            latex.SetTextSize(0.040);
            latex.DrawLatex(0.16, 0.84, "Empty target");
            latex.SetTextSize(0.032);
            latex.DrawLatex(
                0.16, 0.78,
                Form("Entries = %.0f",
                     dataAfter[ic][ie]->GetEntries()));

            double error = 0.0;
            const double integral =
                IntegralRange(dataAfter[ic][ie],
                              -60.0, 40.0,
                              &error);

            txt << ie << " \"" << kLabels[ie]
                << "\" cut" << kCuts[ic] << " "
                << integral << " "
                << error << " "
                << IntegralGlobal(dataAfter[ic][ie])
                << "\n";

            output->cd();
            dataAfter[ic][ie]->Write();
        }

        afterCanvas->Print(afterPdf);
    }

    afterCanvas->Print(afterPdf + "]");

    output->cd();
    beforeCanvas->Write();
    afterCanvas->Write();
    output->Close();
    txt.close();

    std::cout
        << "\nCreated:\n"
        << "  " << beforePdf << "\n"
        << "  " << afterPdf << "\n"
        << "  " << rootName << "\n"
        << "  " << textName << "\n"
        << std::endl;
}
