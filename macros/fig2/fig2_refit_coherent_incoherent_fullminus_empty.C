/*
 * fig2_refit_coherent_incoherent_fullminus_empty.C
 *
 * PURPOSE
 *   Fit the FULL-EMPTY missing-energy distributions with coherent and breakup
 *   Monte Carlo templates, independently for each energy interval and for the
 *   selections before the DeltaPhi cut and after DeltaPhi = 8, 10, 12 deg.
 *
 * DATA INPUT
 *   dataHistFile must be a ROOT file produced by the FULL-EMPTY subtraction
 *   macro and must contain:
 *     he4_before_E0 ... he4_before_E3
 *     he4_cut8_E0  ... he4_cut8_E3
 *     he4_cut10_E0 ... he4_cut10_E3
 *     he4_cut12_E0 ... he4_cut12_E3
 *
 * COHERENT INPUT
 *   coherentFile8 / coherentFile10 / coherentFile12 are the ROOT files made
 *   by paper_fig2_missing_energy_root6.C for the three opening-angle cuts.
 *   The 'before cut' coherent histogram is taken from coherentFile8; the
 *   before-cut histogram is independent of the selected DeltaPhi threshold.
 *
 * BREAKUP INPUTS
 *   he3nFile supplies the 3He+n template.
 *   tpFile supplies the t+p template when a non-empty histogram is available.
 *
 * WHAT THE FIT DOES
 *   The macro performs two separate non-negative two-template fits:
 *     data = a_coh * coherent + a_inc * (3He+n)
 *   and, when available,
 *     data = a_coh * coherent + a_inc * (t+p)
 *
 *   These are independent comparisons. The macro does NOT perform one
 *   simultaneous coherent + 3He+n + t+p three-component fit.
 *
 *   Default fit interval:  -60 to +40 MeV
 *   Default tail interval: -60 to -20 MeV
 *
 * OUTPUTS
 *   Generated from outputTag:
 *     fig2_refit_coh_incoh_<tag>_before.pdf
 *     fig2_refit_coh_incoh_<tag>_after.pdf
 *     fig2_refit_coh_incoh_<tag>.txt
 *     fig2_refit_coh_incoh_<tag>.root
 *
 * NOTE FOR THE REPOSITORY
 *   Historical defaults such as 'full33_minus_empty5' are kept for backward
 *   compatibility. Pass explicit input paths and outputTag for new analyses.
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




namespace RefitAfterCuts {

using namespace Fig2FromAcqu;

const int kNPanels = 4;
const char* kPanelKey[kNPanels] = {
    "before", "cut8", "cut10", "cut12"
};

const char* kPanelTitle[kNPanels] = {
    "before #Delta#Phi cut",
    "#Delta#Phi < 8^{#circ}",
    "#Delta#Phi < 10^{#circ}",
    "#Delta#Phi < 12^{#circ}"
};

struct FitResult {
    double coherentScale;
    double incoherentScale;
    double chi2;
    int ndf;
    double chi2ndf;
    bool boundary;
};

double FitChi2(
    TH1D* data,
    TH1D* coherent,
    TH1D* incoherent,
    double coherentScale,
    double incoherentScale,
    double xmin,
    double xmax)
{
    if (!data || !coherent || !incoherent)
        return 1.0e300;

    const int b1 =
        data->GetXaxis()->FindBin(xmin + 1.0e-9);

    const int b2 =
        data->GetXaxis()->FindBin(xmax - 1.0e-9);

    double chi2 = 0.0;

    for (int b = b1; b <= b2; ++b) {
        const double error = data->GetBinError(b);

        if (error <= 0.0)
            continue;

        const double model =
            coherentScale*coherent->GetBinContent(b) +
            incoherentScale*incoherent->GetBinContent(b);

        const double difference =
            data->GetBinContent(b) - model;

        chi2 += difference*difference/(error*error);
    }

    return chi2;
}

FitResult FitCoherentIncoherent(
    TH1D* data,
    TH1D* coherent,
    TH1D* incoherent,
    double xmin,
    double xmax)
{
    FitResult result = {
        0.0, 0.0,
        0.0, 0, 0.0,
        true
    };

    if (!data || !coherent || !incoherent)
        return result;

    const int b1 =
        data->GetXaxis()->FindBin(xmin + 1.0e-9);

    const int b2 =
        data->GetXaxis()->FindBin(xmax - 1.0e-9);

    double a11 = 0.0;
    double a12 = 0.0;
    double a22 = 0.0;
    double rhs1 = 0.0;
    double rhs2 = 0.0;

    int usedBins = 0;

    for (int b = b1; b <= b2; ++b) {
        const double error = data->GetBinError(b);

        if (error <= 0.0)
            continue;

        const double weight = 1.0/(error*error);
        const double y = data->GetBinContent(b);
        const double c = coherent->GetBinContent(b);
        const double i = incoherent->GetBinContent(b);

        a11 += weight*c*c;
        a12 += weight*c*i;
        a22 += weight*i*i;
        rhs1 += weight*c*y;
        rhs2 += weight*i*y;

        ++usedBins;
    }

    double bestC = 0.0;
    double bestI = 0.0;

    double bestChi2 =
        FitChi2(
            data, coherent, incoherent,
            bestC, bestI,
            xmin, xmax);

    if (a11 > 0.0) {
        const double c =
            std::max(0.0, rhs1/a11);

        const double chi2 =
            FitChi2(
                data, coherent, incoherent,
                c, 0.0,
                xmin, xmax);

        if (chi2 < bestChi2) {
            bestChi2 = chi2;
            bestC = c;
            bestI = 0.0;
        }
    }

    if (a22 > 0.0) {
        const double i =
            std::max(0.0, rhs2/a22);

        const double chi2 =
            FitChi2(
                data, coherent, incoherent,
                0.0, i,
                xmin, xmax);

        if (chi2 < bestChi2) {
            bestChi2 = chi2;
            bestC = 0.0;
            bestI = i;
        }
    }

    const double determinant =
        a11*a22 - a12*a12;

    if (determinant > 0.0) {
        const double c =
            (rhs1*a22 - rhs2*a12)/determinant;

        const double i =
            (rhs2*a11 - rhs1*a12)/determinant;

        if (c >= 0.0 && i >= 0.0) {
            const double chi2 =
                FitChi2(
                    data, coherent, incoherent,
                    c, i,
                    xmin, xmax);

            if (chi2 < bestChi2) {
                bestChi2 = chi2;
                bestC = c;
                bestI = i;
            }
        }
    }

    const int activeParameters =
        (bestC > 0.0 ? 1 : 0) +
        (bestI > 0.0 ? 1 : 0);

    result.coherentScale = bestC;
    result.incoherentScale = bestI;
    result.chi2 = bestChi2;
    result.ndf =
        std::max(0, usedBins - activeParameters);

    result.chi2ndf =
        result.ndf > 0
        ? result.chi2/result.ndf
        : 0.0;

    result.boundary =
        !(bestC > 0.0 && bestI > 0.0);

    return result;
}

TH1D* ScaledClone(
    TH1D* source,
    double scale,
    const char* name)
{
    if (!source)
        return 0;

    TH1D* result =
        dynamic_cast<TH1D*>(
            source->Clone(name));

    if (!result)
        return 0;

    result->SetDirectory(0);
    result->Scale(scale);

    return result;
}

TH1D* SumModel(
    TH1D* coherent,
    TH1D* incoherent,
    const char* name)
{
    if (!coherent || !incoherent)
        return 0;

    TH1D* result =
        dynamic_cast<TH1D*>(
            coherent->Clone(name));

    if (!result)
        return 0;

    result->SetDirectory(0);
    result->Add(incoherent);

    return result;
}

void StyleFittedCoherent(TH1D* histogram)
{
    if (!histogram) return;

    histogram->SetLineColor(kRed + 1);
    histogram->SetLineWidth(2);
    histogram->SetLineStyle(2);
}

void StyleFittedHe3n(TH1D* histogram)
{
    if (!histogram) return;

    histogram->SetLineColor(kBlue + 1);
    histogram->SetLineWidth(2);
    histogram->SetLineStyle(7);
}

void StyleTotalHe3n(TH1D* histogram)
{
    if (!histogram) return;

    histogram->SetLineColor(kGreen + 2);
    histogram->SetLineWidth(3);
    histogram->SetLineStyle(1);
}

void StyleTotalTP(TH1D* histogram)
{
    if (!histogram) return;

    histogram->SetLineColor(kMagenta + 2);
    histogram->SetLineWidth(3);
    histogram->SetLineStyle(7);
}

double TailSigma(
    TH1D* data,
    TH1D* model,
    double tailMin,
    double tailMax)
{
    double residual = 0.0;
    double residualError = 0.0;

    return TailResidualSigma(
        data, model,
        tailMin, tailMax,
        residual, residualError);
}

} // namespace RefitAfterCuts


void fig2_refit_coherent_incoherent_fullminus_empty(
    const char* dataHistFile =
        "fig2_full_empty_full33_minus_empty5.root",
    const char* coherentFile8 =
        "../paper_fig2_missing_energy_cut8.root",
    const char* coherentFile10 =
        "../paper_fig2_missing_energy_cut10.root",
    const char* coherentFile12 =
        "../paper_fig2_missing_energy_cut12.root",
    const char* he3nFile =
        "fig2_he3n_all_available_bins_wide.root",
    const char* tpFile =
        "fig2_tp_E3_366_wide.root",
    const char* outputTag =
        "full33_minus_empty5",
    double fitMin = -60.0,
    double fitMax = 40.0,
    double tailMin = -60.0,
    double tailMax = -20.0)
{
    using namespace Fig2FromAcqu;
    using namespace RefitAfterCuts;

    gStyle->SetOptStat(0);

    TString tag(outputTag);
    tag.ReplaceAll("/", "_");
    tag.ReplaceAll(" ", "_");

    TFile* coherentFiles[kNCuts] = {
        TFile::Open(coherentFile8),
        TFile::Open(coherentFile10),
        TFile::Open(coherentFile12)
    };

    for (int ic = 0; ic < kNCuts; ++ic) {
        if (!coherentFiles[ic] ||
            coherentFiles[ic]->IsZombie()) {
            std::cerr
                << "Cannot open coherent template file "
                << (ic == 0 ? coherentFile8 :
                    ic == 1 ? coherentFile10 :
                              coherentFile12)
                << std::endl;
            return;
        }
    }

    TFile* fHe3n = TFile::Open(he3nFile);
    TFile* fTP = TFile::Open(tpFile);
    TFile* fData = TFile::Open(dataHistFile);

    if (!fHe3n || fHe3n->IsZombie()) {
        std::cerr
            << "Cannot open " << he3nFile
            << std::endl;
        return;
    }

    if (!fTP || fTP->IsZombie()) {
        std::cerr
            << "Cannot open " << tpFile
            << std::endl;
        return;
    }

    if (!fData || fData->IsZombie()) {
        std::cerr
            << "Cannot open full-empty histogram file "
            << dataHistFile << std::endl;
        return;
    }

    TH1D* dataBefore[kNE] = {0};
    TH1D* dataAfter[kNCuts][kNE] = {{0}};

    for (int ie = 0; ie < kNE; ++ie) {
        TH1D* sourceBefore =
            dynamic_cast<TH1D*>(
                fData->Get(
                    Form("he4_before_E%d", ie)));

        if (!sourceBefore) {
            std::cerr
                << "Missing he4_before_E" << ie
                << " in " << dataHistFile
                << std::endl;
            return;
        }

        dataBefore[ie] =
            dynamic_cast<TH1D*>(
                sourceBefore->Clone(
                    Form("refit_he4_before_E%d", ie)));

        dataBefore[ie]->SetDirectory(0);
        dataBefore[ie]->Sumw2();

        for (int ic = 0; ic < kNCuts; ++ic) {
            TH1D* sourceAfter =
                dynamic_cast<TH1D*>(
                    fData->Get(
                        Form("he4_cut%d_E%d",
                             kCuts[ic], ie)));

            if (!sourceAfter) {
                std::cerr
                    << "Missing he4_cut"
                    << kCuts[ic] << "_E" << ie
                    << " in " << dataHistFile
                    << std::endl;
                return;
            }

            dataAfter[ic][ie] =
                dynamic_cast<TH1D*>(
                    sourceAfter->Clone(
                        Form("refit_he4_cut%d_E%d",
                             kCuts[ic], ie)));

            dataAfter[ic][ie]->SetDirectory(0);
            dataAfter[ic][ie]->Sumw2();
        }
    }

    fData->Close();
    delete fData;

    const TString beforePdf =
        Form(
            "fig2_refit_coh_incoh_%s_before.pdf",
            tag.Data());

    const TString afterPdf =
        Form(
            "fig2_refit_coh_incoh_%s_after.pdf",
            tag.Data());

    const TString textName =
        Form(
            "fig2_refit_coh_incoh_%s.txt",
            tag.Data());

    const TString rootName =
        Form(
            "fig2_refit_coh_incoh_%s.root",
            tag.Data());

    std::ofstream textOutput(textName.Data());

    textOutput
        << "# Independent coherent+incoherent refit "
           "for every energy and every DeltaPhi selection\n"
        << "# Full-empty histogram file: " << dataHistFile << "\n"
        << "# Fit interval: ["
        << fitMin << "," << fitMax << "] MeV\n"
        << "# Tail interval: ["
        << tailMin << "," << tailMax << "] MeV\n"
        << "# panel energy model coherent_scale "
           "incoherent_scale chi2 ndf chi2ndf "
           "tail_sigma boundary\n";

    TFile* rootOutput =
        TFile::Open(rootName, "RECREATE");

    if (!rootOutput || rootOutput->IsZombie()) {
        std::cerr
            << "Cannot create " << rootName
            << std::endl;
        return;
    }

    TCanvas* beforeCanvas =
        new TCanvas(
            "c_refit_before",
            "independent refit before cut",
            1300, 900);

    beforeCanvas->Divide(
        2, 2, 0.002, 0.002);

    TCanvas* afterCanvas =
        new TCanvas(
            "c_refit_after",
            "independent refit after cuts",
            1300, 900);

    afterCanvas->Print(afterPdf + "[");

    for (int panel = 0;
         panel < kNPanels;
         ++panel) {

        if (panel > 0) {
            afterCanvas->Clear();
            afterCanvas->Divide(
                2, 2, 0.002, 0.002);
        }

        for (int ie = 0; ie < kNE; ++ie) {
            TH1D* data =
                panel == 0
                ? dataBefore[ie]
                : dataAfter[panel - 1][ie];

            TH1* coherentRaw = 0;
            TH1* he3nRaw = 0;
            TH1* tpRaw = 0;

            if (panel == 0) {
                coherentRaw =
                    dynamic_cast<TH1*>(
                        coherentFiles[0]->Get(
                            Form(
                                "coh_mc_deltaE_before_phi_E%d",
                                ie)));

                he3nRaw =
                    dynamic_cast<TH1*>(
                        fHe3n->Get(
                            Form(
                                "inc_deltaE_E%d_nodphi",
                                ie)));

                tpRaw =
                    dynamic_cast<TH1*>(
                        fTP->Get(
                            Form(
                                "inc_deltaE_E%d_nodphi",
                                ie)));
            } else {
                const int ic = panel - 1;

                coherentRaw =
                    dynamic_cast<TH1*>(
                        coherentFiles[ic]->Get(
                            Form(
                                "coh_mc_deltaE_after_phi_E%d",
                                ie)));

                he3nRaw =
                    dynamic_cast<TH1*>(
                        fHe3n->Get(
                            Form(
                                "inc_deltaE_E%d_dphi%d",
                                ie, kCuts[ic])));

                tpRaw =
                    dynamic_cast<TH1*>(
                        fTP->Get(
                            Form(
                                "inc_deltaE_E%d_dphi%d",
                                ie, kCuts[ic])));
            }

            if (!data ||
                !coherentRaw ||
                !he3nRaw)
                continue;

            TH1D* coherent =
                MapToReference(
                    coherentRaw,
                    data,
                    Form(
                        "coh_raw_%s_E%d",
                        kPanelKey[panel], ie));

            TH1D* he3n =
                MapToReference(
                    he3nRaw,
                    data,
                    Form(
                        "he3n_raw_%s_E%d",
                        kPanelKey[panel], ie));

            TH1D* tp =
                tpRaw &&
                IntegralGlobal(tpRaw) > 0.0
                ? MapToReference(
                    tpRaw,
                    data,
                    Form(
                        "tp_raw_%s_E%d",
                        kPanelKey[panel], ie))
                : 0;

            if (!coherent || !he3n)
                continue;

            const FitResult fitHe3n =
                FitCoherentIncoherent(
                    data,
                    coherent,
                    he3n,
                    fitMin,
                    fitMax);

            TH1D* coherentHe3n =
                ScaledClone(
                    coherent,
                    fitHe3n.coherentScale,
                    Form(
                        "coherent_he3n_%s_E%d",
                        kPanelKey[panel], ie));

            TH1D* fittedHe3n =
                ScaledClone(
                    he3n,
                    fitHe3n.incoherentScale,
                    Form(
                        "he3n_fitted_%s_E%d",
                        kPanelKey[panel], ie));

            TH1D* totalHe3n =
                SumModel(
                    coherentHe3n,
                    fittedHe3n,
                    Form(
                        "total_he3n_%s_E%d",
                        kPanelKey[panel], ie));

            const double tailSigmaHe3n =
                TailSigma(
                    data,
                    totalHe3n,
                    tailMin,
                    tailMax);

            textOutput
                << kPanelKey[panel] << " "
                << kLabels[ie] << " "
                << "he3n "
                << fitHe3n.coherentScale << " "
                << fitHe3n.incoherentScale << " "
                << fitHe3n.chi2 << " "
                << fitHe3n.ndf << " "
                << fitHe3n.chi2ndf << " "
                << tailSigmaHe3n << " "
                << (fitHe3n.boundary ? 1 : 0)
                << "\n";

            FitResult fitTP = {
                0.0, 0.0,
                0.0, 0, 0.0,
                true
            };

            TH1D* coherentTP = 0;
            TH1D* fittedTP = 0;
            TH1D* totalTP = 0;
            double tailSigmaTP = 0.0;

            if (tp) {
                fitTP =
                    FitCoherentIncoherent(
                        data,
                        coherent,
                        tp,
                        fitMin,
                        fitMax);

                coherentTP =
                    ScaledClone(
                        coherent,
                        fitTP.coherentScale,
                        Form(
                            "coherent_tp_%s_E%d",
                            kPanelKey[panel], ie));

                fittedTP =
                    ScaledClone(
                        tp,
                        fitTP.incoherentScale,
                        Form(
                            "tp_fitted_%s_E%d",
                            kPanelKey[panel], ie));

                totalTP =
                    SumModel(
                        coherentTP,
                        fittedTP,
                        Form(
                            "total_tp_%s_E%d",
                            kPanelKey[panel], ie));

                tailSigmaTP =
                    TailSigma(
                        data,
                        totalTP,
                        tailMin,
                        tailMax);

                textOutput
                    << kPanelKey[panel] << " "
                    << kLabels[ie] << " "
                    << "tp "
                    << fitTP.coherentScale << " "
                    << fitTP.incoherentScale << " "
                    << fitTP.chi2 << " "
                    << fitTP.ndf << " "
                    << fitTP.chi2ndf << " "
                    << tailSigmaTP << " "
                    << (fitTP.boundary ? 1 : 0)
                    << "\n";
            }

            StyleData(data);
            StyleFittedCoherent(coherentHe3n);
            StyleFittedHe3n(fittedHe3n);
            StyleTotalHe3n(totalHe3n);
            StyleTotalTP(totalTP);

            TCanvas* canvas =
                panel == 0
                ? beforeCanvas
                : afterCanvas;

            canvas->cd(ie + 1);

            gPad->SetTicks(1, 1);
            gPad->SetLeftMargin(0.13);
            gPad->SetBottomMargin(0.13);

            data->SetTitle(
                Form(
                    "%s, %s",
                    kLabels[ie],
                    kPanelTitle[panel]));

            data->GetXaxis()->SetTitle(
                "#Delta E_{#pi^{0}}^{*} (MeV)");

            data->GetYaxis()->SetTitle(
                "Counts");

            data->GetXaxis()->SetRangeUser(
                fitMin,
                fitMax);

            double maximum =
                std::max(
                    data->GetMaximum(),
                    totalHe3n->GetMaximum());

            if (totalTP) {
                maximum =
                    std::max(
                        maximum,
                        totalTP->GetMaximum());
            }

            data->SetMinimum(0.0);
            data->SetMaximum(
                maximum > 0.0
                ? 1.38*maximum
                : 1.0);

            data->Draw("E");
            coherentHe3n->Draw("HIST SAME");
            fittedHe3n->Draw("HIST SAME");
            totalHe3n->Draw("HIST SAME");

            if (totalTP)
                totalTP->Draw("HIST SAME");

            data->Draw("E SAME");

            TLegend* legend =
                new TLegend(
                    0.42, 0.54,
                    0.89, 0.88);

            legend->SetBorderSize(0);
            legend->SetFillStyle(0);
            legend->SetTextSize(0.026);

            legend->AddEntry(
                data,
                "Full - empty",
                "lep");

            legend->AddEntry(
                coherentHe3n,
                "Coherent MC (refitted)",
                "l");

            legend->AddEntry(
                fittedHe3n,
                "^{3}He+n MC (refitted)",
                "l");

            legend->AddEntry(
                totalHe3n,
                "Coherent + ^{3}He+n",
                "l");

            if (totalTP) {
                legend->AddEntry(
                    totalTP,
                    "Coherent + t+p (independent fit)",
                    "l");
            }

            legend->Draw();

            TLatex label;
            label.SetNDC();
            label.SetTextSize(0.026);

            label.DrawLatex(
                0.14, 0.91,
                Form(
                    "^{3}He+n: a_{coh}=%.3g, a_{inc}=%.3g",
                    fitHe3n.coherentScale,
                    fitHe3n.incoherentScale));

            label.DrawLatex(
                0.14, 0.86,
                Form(
                    "#chi^{2}/ndf=%.2f, tail=%.1f#sigma",
                    fitHe3n.chi2ndf,
                    tailSigmaHe3n));

            if (totalTP) {
                label.DrawLatex(
                    0.14, 0.81,
                    Form(
                        "t+p: #chi^{2}/ndf=%.2f, tail=%.1f#sigma",
                        fitTP.chi2ndf,
                        tailSigmaTP));
            }

            rootOutput->cd();

            data->Write(
                Form(
                    "data_%s_E%d",
                    kPanelKey[panel], ie));

            coherentHe3n->Write();
            fittedHe3n->Write();
            totalHe3n->Write();

            if (coherentTP)
                coherentTP->Write();

            if (fittedTP)
                fittedTP->Write();

            if (totalTP)
                totalTP->Write();
        }

        if (panel > 0)
            afterCanvas->Print(afterPdf);
    }

    beforeCanvas->SaveAs(beforePdf);
    afterCanvas->Print(afterPdf + "]");

    rootOutput->cd();
    beforeCanvas->Write();
    afterCanvas->Write();

    textOutput.close();
    rootOutput->Close();

    for (int ic = 0; ic < kNCuts; ++ic)
        coherentFiles[ic]->Close();

    fHe3n->Close();
    fTP->Close();

    std::cout
        << "\nIndependent post-cut refits completed.\n"
        << "Saved:\n"
        << "  " << beforePdf << "\n"
        << "  " << afterPdf << "\n"
        << "  " << textName << "\n"
        << "  " << rootName << "\n";
}
