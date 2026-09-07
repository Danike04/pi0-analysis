/*
 * fig2_breakup_normalization_study_all.C
 *
 * Purpose
 * -------
 * Campaign-specific Fig. 2 study comparing several normalization strategies
 * for breakup templates.  Experimental Acqu data are reconstructed with the
 * same pi0 / prompt-random / DeltaPhi selections used in the surrounding
 * Fig. 2 studies.  Coherent templates are taken from the three historical
 * DeltaPhi-cut ROOT files.  ppnn is combined in turn with 3He+n and t+p and
 * several ways of fixing the breakup normalizations are compared.
 *
 * Main inputs (historical defaults)
 * ---------------------------------
 *   dataFile       : Acqu_CBTagg_31832_31836-31850.root
 *   coherentFile8  : paper_fig2_missing_energy_cut8.root
 *   coherentFile10 : paper_fig2_missing_energy_cut10.root
 *   coherentFile12 : paper_fig2_missing_energy_cut12.root
 *   ppnnFile       : fig2_incoherent_4bins_1M_wide.root
 *   he3nFile       : fig2_he3n_all_available_bins_wide.root
 *   tpFile         : fig2_tp_E3_366_wide.root
 *
 * Default analysis regions
 * ------------------------
 *   coherent core : -10 to +15 MeV
 *   breakup tail  : -60 to -20 MeV
 *   comparison    : -60 to +40 MeV
 *
 * Outputs
 * -------
 * Files begin with "fig2_breakup_normstudy_<outputTag>" and include a ROOT
 * file, text summaries/scans, and PDF comparisons of the tested methods.
 *
 * Usage
 * -----
 *   root -l -b -q 'fig2_breakup_normalization_study_all.C()'
 *
 * Status
 * ------
 * Physics study / normalization cross-check, not a general-purpose production
 * step.  The code below is preserved exactly from the supplied analysis file.
 */
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>

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



namespace Fig2FromAcqu {

void StyleModelHe3n(TH1D* h)
{
    h->SetLineColor(kGreen + 2);
    h->SetLineWidth(3);
    h->SetLineStyle(1);
    h->SetFillStyle(0);
}

void StyleModelTP(TH1D* h)
{
    h->SetLineColor(kMagenta + 2);
    h->SetLineWidth(3);
    h->SetLineStyle(7);
    h->SetFillStyle(0);
}

} // namespace Fig2FromAcqu



namespace BreakupNormalizationStudy {

using namespace Fig2FromAcqu;

const int kNStudyPanels = 4;
const char* kPanelNames[kNStudyPanels] = {
    "before",
    "dphi8",
    "dphi10",
    "dphi12"
};

const char* kPanelTitles[kNStudyPanels] = {
    "before #Delta#Phi cut",
    "#Delta#Phi < 8^{#circ}",
    "#Delta#Phi < 10^{#circ}",
    "#Delta#Phi < 12^{#circ}"
};

struct OneFit {
    double scale;
    double error;
    double chi2;
    int ndf;
    int usedBins;
};

struct TwoFit {
    double ppnnScale;
    double secondaryScale;
    double ppnnError;
    double secondaryError;
    double correlation;
    double chi2;
    int ndf;
    int usedBins;
    int activeParameters;
    bool boundary;
};

struct ModelMetrics {
    double chi2;
    int ndf;
    double chi2ndf;
};

double FitVariance(TH1D* data, TH1D* coherent, int bin)
{
    if (!data || !coherent) return 0.0;

    const double ed = data->GetBinError(bin);
    const double ec = coherent->GetBinError(bin);

    return ed*ed + ec*ec;
}

OneFit FitOneTemplate(
    TH1D* data[kNStudyPanels],
    TH1D* coherent[kNStudyPanels],
    TH1D* templ[kNStudyPanels],
    int firstPanel,
    int lastPanel,
    double xmin,
    double xmax,
    double fixedPPNNScale = 0.0,
    TH1D* ppnn[kNStudyPanels] = 0)
{
    OneFit result = {0.0, 0.0, 0.0, 0, 0};

    double a11 = 0.0;
    double rhs = 0.0;

    for (int ip = firstPanel; ip <= lastPanel; ++ip) {
        if (!data[ip] || !coherent[ip] || !templ[ip])
            continue;

        const int b1 =
            data[ip]->GetXaxis()->FindBin(xmin + 1.0e-9);
        const int b2 =
            data[ip]->GetXaxis()->FindBin(xmax - 1.0e-9);

        for (int b = b1; b <= b2; ++b) {
            const double var =
                FitVariance(data[ip], coherent[ip], b);

            if (var <= 0.0)
                continue;

            double y =
                data[ip]->GetBinContent(b) -
                coherent[ip]->GetBinContent(b);

            if (ppnn && ppnn[ip]) {
                y -= fixedPPNNScale *
                     ppnn[ip]->GetBinContent(b);
            }

            const double x =
                templ[ip]->GetBinContent(b);

            if (x == 0.0)
                continue;

            const double w = 1.0/var;

            a11 += w*x*x;
            rhs += w*x*y;
            ++result.usedBins;
        }
    }

    if (a11 > 0.0) {
        result.scale =
            std::max(0.0, rhs/a11);

        result.error =
            1.0/std::sqrt(a11);
    }

    for (int ip = firstPanel; ip <= lastPanel; ++ip) {
        if (!data[ip] || !coherent[ip] || !templ[ip])
            continue;

        const int b1 =
            data[ip]->GetXaxis()->FindBin(xmin + 1.0e-9);
        const int b2 =
            data[ip]->GetXaxis()->FindBin(xmax - 1.0e-9);

        for (int b = b1; b <= b2; ++b) {
            const double var =
                FitVariance(data[ip], coherent[ip], b);

            if (var <= 0.0)
                continue;

            double model =
                coherent[ip]->GetBinContent(b) +
                result.scale*templ[ip]->GetBinContent(b);

            if (ppnn && ppnn[ip]) {
                model += fixedPPNNScale *
                         ppnn[ip]->GetBinContent(b);
            }

            const double diff =
                data[ip]->GetBinContent(b) - model;

            result.chi2 += diff*diff/var;
        }
    }

    result.ndf =
        std::max(0, result.usedBins -
                    (result.scale > 0.0 ? 1 : 0));

    return result;
}

double TwoFitSSE(
    TH1D* data[kNStudyPanels],
    TH1D* coherent[kNStudyPanels],
    TH1D* ppnn[kNStudyPanels],
    TH1D* secondary[kNStudyPanels],
    int firstPanel,
    int lastPanel,
    double xmin,
    double xmax,
    double ppnnScale,
    double secondaryScale)
{
    double chi2 = 0.0;

    for (int ip = firstPanel; ip <= lastPanel; ++ip) {
        if (!data[ip] || !coherent[ip] ||
            !ppnn[ip] || !secondary[ip])
            continue;

        const int b1 =
            data[ip]->GetXaxis()->FindBin(xmin + 1.0e-9);
        const int b2 =
            data[ip]->GetXaxis()->FindBin(xmax - 1.0e-9);

        for (int b = b1; b <= b2; ++b) {
            const double var =
                FitVariance(data[ip], coherent[ip], b);

            if (var <= 0.0)
                continue;

            const double model =
                coherent[ip]->GetBinContent(b) +
                ppnnScale*ppnn[ip]->GetBinContent(b) +
                secondaryScale*
                    secondary[ip]->GetBinContent(b);

            const double diff =
                data[ip]->GetBinContent(b) - model;

            chi2 += diff*diff/var;
        }
    }

    return chi2;
}

TwoFit FitTwoTemplates(
    TH1D* data[kNStudyPanels],
    TH1D* coherent[kNStudyPanels],
    TH1D* ppnn[kNStudyPanels],
    TH1D* secondary[kNStudyPanels],
    int firstPanel,
    int lastPanel,
    double xmin,
    double xmax)
{
    TwoFit result = {
        0.0, 0.0,
        0.0, 0.0,
        0.0,
        0.0, 0, 0, 0,
        true
    };

    double a11 = 0.0;
    double a12 = 0.0;
    double a22 = 0.0;
    double rhs1 = 0.0;
    double rhs2 = 0.0;

    for (int ip = firstPanel; ip <= lastPanel; ++ip) {
        if (!data[ip] || !coherent[ip] ||
            !ppnn[ip] || !secondary[ip])
            continue;

        const int b1 =
            data[ip]->GetXaxis()->FindBin(xmin + 1.0e-9);
        const int b2 =
            data[ip]->GetXaxis()->FindBin(xmax - 1.0e-9);

        for (int b = b1; b <= b2; ++b) {
            const double var =
                FitVariance(data[ip], coherent[ip], b);

            if (var <= 0.0)
                continue;

            const double w = 1.0/var;

            const double y =
                data[ip]->GetBinContent(b) -
                coherent[ip]->GetBinContent(b);

            const double p =
                ppnn[ip]->GetBinContent(b);

            const double s =
                secondary[ip]->GetBinContent(b);

            a11 += w*p*p;
            a12 += w*p*s;
            a22 += w*s*s;
            rhs1 += w*p*y;
            rhs2 += w*s*y;

            ++result.usedBins;
        }
    }

    double bestP = 0.0;
    double bestS = 0.0;

    double bestChi2 =
        TwoFitSSE(
            data, coherent, ppnn, secondary,
            firstPanel, lastPanel,
            xmin, xmax,
            bestP, bestS);

    if (a11 > 0.0) {
        const double p =
            std::max(0.0, rhs1/a11);

        const double chi2 =
            TwoFitSSE(
                data, coherent, ppnn, secondary,
                firstPanel, lastPanel,
                xmin, xmax,
                p, 0.0);

        if (chi2 < bestChi2) {
            bestChi2 = chi2;
            bestP = p;
            bestS = 0.0;
        }
    }

    if (a22 > 0.0) {
        const double s =
            std::max(0.0, rhs2/a22);

        const double chi2 =
            TwoFitSSE(
                data, coherent, ppnn, secondary,
                firstPanel, lastPanel,
                xmin, xmax,
                0.0, s);

        if (chi2 < bestChi2) {
            bestChi2 = chi2;
            bestP = 0.0;
            bestS = s;
        }
    }

    const double det =
        a11*a22 - a12*a12;

    if (det > 0.0) {
        const double p =
            (rhs1*a22 - rhs2*a12)/det;

        const double s =
            (rhs2*a11 - rhs1*a12)/det;

        if (p >= 0.0 && s >= 0.0) {
            const double chi2 =
                TwoFitSSE(
                    data, coherent,
                    ppnn, secondary,
                    firstPanel, lastPanel,
                    xmin, xmax,
                    p, s);

            if (chi2 < bestChi2) {
                bestChi2 = chi2;
                bestP = p;
                bestS = s;
            }
        }
    }

    result.ppnnScale = bestP;
    result.secondaryScale = bestS;
    result.chi2 = bestChi2;

    result.activeParameters =
        (bestP > 0.0 ? 1 : 0) +
        (bestS > 0.0 ? 1 : 0);

    result.ndf =
        std::max(
            0,
            result.usedBins -
            result.activeParameters);

    result.boundary =
        !(bestP > 0.0 && bestS > 0.0);

    if (!result.boundary && det > 0.0) {
        result.ppnnError =
            std::sqrt(a22/det);

        result.secondaryError =
            std::sqrt(a11/det);

        result.correlation =
            -a12/std::sqrt(a11*a22);
    } else if (bestP > 0.0 && a11 > 0.0) {
        result.ppnnError =
            1.0/std::sqrt(a11);
    } else if (bestS > 0.0 && a22 > 0.0) {
        result.secondaryError =
            1.0/std::sqrt(a22);
    }

    return result;
}

TH1D* BuildModel(
    TH1D* coherent,
    TH1D* ppnn,
    TH1D* secondary,
    double ppnnScale,
    double secondaryScale,
    const char* name)
{
    if (!coherent)
        return 0;

    TH1D* model =
        dynamic_cast<TH1D*>(
            coherent->Clone(name));

    if (!model)
        return 0;

    model->SetDirectory(0);

    if (ppnn && ppnnScale != 0.0)
        model->Add(ppnn, ppnnScale);

    if (secondary && secondaryScale != 0.0)
        model->Add(secondary, secondaryScale);

    return model;
}

ModelMetrics EvaluateAllPanels(
    TH1D* data[kNStudyPanels],
    TH1D* coherent[kNStudyPanels],
    TH1D* ppnn[kNStudyPanels],
    TH1D* secondary[kNStudyPanels],
    double ppnnScale,
    double secondaryScale,
    double xmin,
    double xmax,
    int fittedParameters)
{
    ModelMetrics result = {0.0, 0, 0.0};
    int usedBins = 0;

    for (int ip = 0; ip < kNStudyPanels; ++ip) {
        if (!data[ip] || !coherent[ip])
            continue;

        TH1D* model =
            BuildModel(
                coherent[ip],
                ppnn ? ppnn[ip] : 0,
                secondary ? secondary[ip] : 0,
                ppnnScale,
                secondaryScale,
                Form("tmp_eval_model_p%d", ip));

        if (!model)
            continue;

        const Chi2Result panel =
            EvaluateChi2(
                data[ip], model,
                xmin, xmax, 0);

        result.chi2 += panel.chi2;
        usedBins += panel.ndf;

        delete model;
    }

    result.ndf =
        std::max(0, usedBins - fittedParameters);

    result.chi2ndf =
        result.ndf > 0
        ? result.chi2/result.ndf
        : 0.0;

    return result;
}

void ScaleCoherentSet(
    TH1D* data[kNStudyPanels],
    TH1D* coherentUnscaled[kNStudyPanels],
    TH1D* coherentScaled[kNStudyPanels],
    int energyBin,
    double coreMin,
    double coreMax,
    const char* suffix)
{
    for (int ip = 0; ip < kNStudyPanels; ++ip) {
        coherentScaled[ip] = 0;

        if (!data[ip] || !coherentUnscaled[ip])
            continue;

        coherentScaled[ip] =
            dynamic_cast<TH1D*>(
                coherentUnscaled[ip]->Clone(
                    Form(
                        "coh_scaled_E%d_p%d_%s",
                        energyBin, ip, suffix)));

        coherentScaled[ip]->SetDirectory(0);

        const double scale =
            CoreScale(
                data[ip],
                coherentScaled[ip],
                coreMin, coreMax);

        coherentScaled[ip]->Scale(scale);
    }
}

void StyleSinglePPNN(TH1D* h)
{
    if (!h) return;
    h->SetLineColor(kBlue + 1);
    h->SetLineStyle(2);
    h->SetLineWidth(3);
}

void StyleSingleSecondary(TH1D* h)
{
    if (!h) return;
    h->SetLineColor(kGreen + 2);
    h->SetLineStyle(1);
    h->SetLineWidth(3);
}

void StyleSequential(TH1D* h)
{
    if (!h) return;
    h->SetLineColor(kBlue + 1);
    h->SetLineStyle(2);
    h->SetLineWidth(3);
}

void StyleBeforeFit(TH1D* h)
{
    if (!h) return;
    h->SetLineColor(kGreen + 2);
    h->SetLineStyle(1);
    h->SetLineWidth(3);
}

void StyleGlobalFit(TH1D* h)
{
    if (!h) return;
    h->SetLineColor(kMagenta + 2);
    h->SetLineStyle(7);
    h->SetLineWidth(3);
}

void WriteMethodLine(
    std::ofstream& txt,
    int energyBin,
    const char* secondaryName,
    const char* method,
    double ppnnScale,
    double secondaryScale,
    const ModelMetrics& metrics,
    double correlation,
    bool boundary)
{
    txt
        << kLabels[energyBin] << " "
        << secondaryName << " "
        << method << " "
        << ppnnScale << " "
        << secondaryScale << " "
        << metrics.chi2ndf << " "
        << metrics.chi2 << " "
        << metrics.ndf << " "
        << correlation << " "
        << (boundary ? 1 : 0)
        << "\n";
}

} // namespace BreakupNormalizationStudy


void fig2_breakup_normalization_study_all(
    const char* dataFile =
        "../Acqu_CBTagg_31832_31836-31850.root",
    const char* coherentFile8 =
        "../paper_fig2_missing_energy_cut8.root",
    const char* coherentFile10 =
        "../paper_fig2_missing_energy_cut10.root",
    const char* coherentFile12 =
        "../paper_fig2_missing_energy_cut12.root",
    const char* ppnnFile =
        "fig2_incoherent_4bins_1M_wide.root",
    const char* he3nFile =
        "fig2_he3n_all_available_bins_wide.root",
    const char* tpFile =
        "fig2_tp_E3_366_wide.root",
    const char* outputTag =
        "all16",
    Long64_t maxEvents = -1,
    double clusterEnergyMin = 20.0,
    double vetoEnergyMax = 1.0,
    double mggMin = 110.0,
    double mggMax = 155.0,
    double promptMin = 700.0,
    double promptMax = 800.0,
    double randomMin = 450.0,
    double randomMax = 680.0,
    double coreMin = -10.0,
    double coreMax = 15.0,
    double tailMin = -60.0,
    double tailMax = -20.0,
    double compareMin = -60.0,
    double compareMax = 40.0)
{
    using namespace Fig2FromAcqu;
    using namespace BreakupNormalizationStudy;

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
                << "Cannot open coherent file "
                << ic << std::endl;
            return;
        }
    }

    TFile* fPPNN = TFile::Open(ppnnFile);
    TFile* fHe3n = TFile::Open(he3nFile);
    TFile* fTP = TFile::Open(tpFile);

    if (!fPPNN || fPPNN->IsZombie() ||
        !fHe3n || fHe3n->IsZombie() ||
        !fTP || fTP->IsZombie()) {
        std::cerr
            << "Cannot open one or more breakup files."
            << std::endl;
        return;
    }

    TH1D* dataBefore[kNE];
    TH1D* dataAfter[kNCuts][kNE];

    for (int ie = 0; ie < kNE; ++ie) {
        dataBefore[ie] = new TH1D(
            Form("study_data_before_E%d", ie),
            "", 240, -120.0, 120.0);

        dataBefore[ie]->SetDirectory(0);
        dataBefore[ie]->Sumw2();

        for (int ic = 0; ic < kNCuts; ++ic) {
            dataAfter[ic][ie] = new TH1D(
                Form(
                    "study_data_after_cut%d_E%d",
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

    TH1D* data[kNStudyPanels][kNE];
    TH1D* coherentUnscaled[kNStudyPanels][kNE];
    TH1D* coherentDefault[kNStudyPanels][kNE];
    TH1D* ppnn[kNStudyPanels][kNE];
    TH1D* he3n[kNStudyPanels][kNE];
    TH1D* tp[kNStudyPanels][kNE];

    for (int ip = 0; ip < kNStudyPanels; ++ip) {
        for (int ie = 0; ie < kNE; ++ie) {
            data[ip][ie] = 0;
            coherentUnscaled[ip][ie] = 0;
            coherentDefault[ip][ie] = 0;
            ppnn[ip][ie] = 0;
            he3n[ip][ie] = 0;
            tp[ip][ie] = 0;
        }
    }

    for (int ie = 0; ie < kNE; ++ie) {
        data[0][ie] =
            dynamic_cast<TH1D*>(
                dataBefore[ie]->Clone(
                    Form("study_data_p0_E%d", ie)));

        data[0][ie]->SetDirectory(0);

        for (int ic = 0; ic < kNCuts; ++ic) {
            data[ic + 1][ie] =
                dynamic_cast<TH1D*>(
                    dataAfter[ic][ie]->Clone(
                        Form(
                            "study_data_p%d_E%d",
                            ic + 1, ie)));

            data[ic + 1][ie]->SetDirectory(0);
        }
    }

    for (int ie = 0; ie < kNE; ++ie) {
        for (int ip = 0; ip < kNStudyPanels; ++ip) {
            TH1* cRaw = 0;
            TH1* pRaw = 0;
            TH1* hRaw = 0;
            TH1* tRaw = 0;

            if (ip == 0) {
                cRaw = dynamic_cast<TH1*>(
                    coherentFiles[0]->Get(
                        Form(
                            "coh_mc_deltaE_before_phi_E%d",
                            ie)));

                pRaw = dynamic_cast<TH1*>(
                    fPPNN->Get(
                        Form(
                            "inc_deltaE_E%d_nodphi",
                            ie)));

                hRaw = dynamic_cast<TH1*>(
                    fHe3n->Get(
                        Form(
                            "inc_deltaE_E%d_nodphi",
                            ie)));

                tRaw = dynamic_cast<TH1*>(
                    fTP->Get(
                        Form(
                            "inc_deltaE_E%d_nodphi",
                            ie)));
            } else {
                const int ic = ip - 1;

                cRaw = dynamic_cast<TH1*>(
                    coherentFiles[ic]->Get(
                        Form(
                            "coh_mc_deltaE_after_phi_E%d",
                            ie)));

                pRaw = dynamic_cast<TH1*>(
                    fPPNN->Get(
                        Form(
                            "inc_deltaE_E%d_dphi%d",
                            ie, kCuts[ic])));

                hRaw = dynamic_cast<TH1*>(
                    fHe3n->Get(
                        Form(
                            "inc_deltaE_E%d_dphi%d",
                            ie, kCuts[ic])));

                tRaw = dynamic_cast<TH1*>(
                    fTP->Get(
                        Form(
                            "inc_deltaE_E%d_dphi%d",
                            ie, kCuts[ic])));
            }

            if (cRaw) {
                coherentUnscaled[ip][ie] =
                    MapToReference(
                        cRaw,
                        data[ip][ie],
                        Form(
                            "study_coh_unscaled_p%d_E%d",
                            ip, ie));
            }

            if (pRaw) {
                ppnn[ip][ie] =
                    MapToReference(
                        pRaw,
                        data[ip][ie],
                        Form(
                            "study_ppnn_p%d_E%d",
                            ip, ie));
            }

            if (hRaw && IntegralGlobal(hRaw) > 0.0) {
                he3n[ip][ie] =
                    MapToReference(
                        hRaw,
                        data[ip][ie],
                        Form(
                            "study_he3n_p%d_E%d",
                            ip, ie));
            }

            if (tRaw && IntegralGlobal(tRaw) > 0.0) {
                tp[ip][ie] =
                    MapToReference(
                        tRaw,
                        data[ip][ie],
                        Form(
                            "study_tp_p%d_E%d",
                            ip, ie));
            }
        }

        TH1D* dataEnergy[kNStudyPanels];
        TH1D* cohUnscaledEnergy[kNStudyPanels];
        TH1D* cohScaledEnergy[kNStudyPanels];

        for (int ip = 0; ip < kNStudyPanels; ++ip) {
            dataEnergy[ip] = data[ip][ie];
            cohUnscaledEnergy[ip] =
                coherentUnscaled[ip][ie];
            cohScaledEnergy[ip] = 0;
        }

        ScaleCoherentSet(
            dataEnergy,
            cohUnscaledEnergy,
            cohScaledEnergy,
            ie,
            coreMin,
            coreMax,
            "default");

        for (int ip = 0; ip < kNStudyPanels; ++ip)
            coherentDefault[ip][ie] =
                cohScaledEnergy[ip];
    }

    const TString rootName =
        Form(
            "fig2_breakup_normstudy_%s.root",
            tag.Data());

    const TString summaryName =
        Form(
            "fig2_breakup_normstudy_%s_summary.txt",
            tag.Data());

    const TString scanName =
        Form(
            "fig2_breakup_normstudy_%s_window_scan.csv",
            tag.Data());

    const TString singlePdf =
        Form(
            "fig2_breakup_normstudy_%s_single_channels.pdf",
            tag.Data());

    const TString he3nPdf =
        Form(
            "fig2_breakup_normstudy_%s_he3n_methods.pdf",
            tag.Data());

    const TString tpPdf =
        Form(
            "fig2_breakup_normstudy_%s_tp_methods.pdf",
            tag.Data());

    TFile* output =
        TFile::Open(rootName, "RECREATE");

    std::ofstream summary(summaryName.Data());
    std::ofstream scan(scanName.Data());

    summary
        << "# Complete breakup-normalization study\n"
        << "# Data file: " << dataFile << "\n"
        << "# Default coherent core: ["
        << coreMin << "," << coreMax << "] MeV\n"
        << "# Default breakup tail: ["
        << tailMin << "," << tailMax << "] MeV\n"
        << "# Evaluation interval: ["
        << compareMin << "," << compareMax << "] MeV\n"
        << "# Methods:\n"
        << "# single: one breakup template fitted before the cut\n"
        << "# sequential: ppnn fitted first before the cut; "
           "secondary fitted to the remaining tail\n"
        << "# before2: ppnn and secondary fitted simultaneously "
           "before the cut\n"
        << "# global2: shared ppnn and secondary scales fitted "
           "simultaneously before and after all DeltaPhi cuts\n"
        << "# energy secondary method ppnn_scale secondary_scale "
           "chi2ndf_all chi2_all ndf_all correlation boundary\n";

    scan
        << "energy,secondary,core_min,core_max,"
           "tail_min,tail_max,ppnn_scale,ppnn_error,"
           "secondary_scale,secondary_error,correlation,"
           "boundary,fit_chi2,fit_ndf,all_chi2ndf\n";

    OneFit singlePPNN[kNE];
    OneFit singleHe3n[kNE];
    OneFit singleTP[kNE];

    TwoFit sequentialHe3n[kNE];
    TwoFit beforeHe3n[kNE];
    TwoFit globalHe3n[kNE];

    TwoFit sequentialTP[kNE];
    TwoFit beforeTP[kNE];
    TwoFit globalTP[kNE];

    for (int ie = 0; ie < kNE; ++ie) {
        TH1D* d[kNStudyPanels];
        TH1D* c[kNStudyPanels];
        TH1D* p[kNStudyPanels];
        TH1D* h[kNStudyPanels];
        TH1D* t[kNStudyPanels];

        for (int ip = 0; ip < kNStudyPanels; ++ip) {
            d[ip] = data[ip][ie];
            c[ip] = coherentDefault[ip][ie];
            p[ip] = ppnn[ip][ie];
            h[ip] = he3n[ip][ie];
            t[ip] = tp[ip][ie];
        }

        singlePPNN[ie] =
            FitOneTemplate(
                d, c, p,
                0, 0,
                tailMin, tailMax);

        if (h[0]) {
            singleHe3n[ie] =
                FitOneTemplate(
                    d, c, h,
                    0, 0,
                    tailMin, tailMax);

            const OneFit seqSecondary =
                FitOneTemplate(
                    d, c, h,
                    0, 0,
                    tailMin, tailMax,
                    singlePPNN[ie].scale,
                    p);

            sequentialHe3n[ie].ppnnScale =
                singlePPNN[ie].scale;
            sequentialHe3n[ie].secondaryScale =
                seqSecondary.scale;
            sequentialHe3n[ie].ppnnError =
                singlePPNN[ie].error;
            sequentialHe3n[ie].secondaryError =
                seqSecondary.error;
            sequentialHe3n[ie].correlation = 0.0;
            sequentialHe3n[ie].boundary =
                sequentialHe3n[ie].ppnnScale <= 0.0 ||
                sequentialHe3n[ie].secondaryScale <= 0.0;
            sequentialHe3n[ie].activeParameters =
                (sequentialHe3n[ie].ppnnScale > 0.0 ? 1 : 0) +
                (sequentialHe3n[ie].secondaryScale > 0.0 ? 1 : 0);

            beforeHe3n[ie] =
                FitTwoTemplates(
                    d, c, p, h,
                    0, 0,
                    tailMin, tailMax);

            globalHe3n[ie] =
                FitTwoTemplates(
                    d, c, p, h,
                    0, kNStudyPanels - 1,
                    tailMin, tailMax);

            const ModelMetrics mSingle =
                EvaluateAllPanels(
                    d, c,
                    0, h,
                    0.0,
                    singleHe3n[ie].scale,
                    compareMin, compareMax,
                    1);

            const ModelMetrics mSequential =
                EvaluateAllPanels(
                    d, c,
                    p, h,
                    sequentialHe3n[ie].ppnnScale,
                    sequentialHe3n[ie].secondaryScale,
                    compareMin, compareMax,
                    sequentialHe3n[ie].activeParameters);

            const ModelMetrics mBefore =
                EvaluateAllPanels(
                    d, c,
                    p, h,
                    beforeHe3n[ie].ppnnScale,
                    beforeHe3n[ie].secondaryScale,
                    compareMin, compareMax,
                    beforeHe3n[ie].activeParameters);

            const ModelMetrics mGlobal =
                EvaluateAllPanels(
                    d, c,
                    p, h,
                    globalHe3n[ie].ppnnScale,
                    globalHe3n[ie].secondaryScale,
                    compareMin, compareMax,
                    globalHe3n[ie].activeParameters);

            WriteMethodLine(
                summary, ie,
                "he3n", "single",
                0.0,
                singleHe3n[ie].scale,
                mSingle,
                0.0,
                singleHe3n[ie].scale <= 0.0);

            WriteMethodLine(
                summary, ie,
                "he3n", "sequential",
                sequentialHe3n[ie].ppnnScale,
                sequentialHe3n[ie].secondaryScale,
                mSequential,
                0.0,
                sequentialHe3n[ie].boundary);

            WriteMethodLine(
                summary, ie,
                "he3n", "before2",
                beforeHe3n[ie].ppnnScale,
                beforeHe3n[ie].secondaryScale,
                mBefore,
                beforeHe3n[ie].correlation,
                beforeHe3n[ie].boundary);

            WriteMethodLine(
                summary, ie,
                "he3n", "global2",
                globalHe3n[ie].ppnnScale,
                globalHe3n[ie].secondaryScale,
                mGlobal,
                globalHe3n[ie].correlation,
                globalHe3n[ie].boundary);
        }

        if (t[0]) {
            singleTP[ie] =
                FitOneTemplate(
                    d, c, t,
                    0, 0,
                    tailMin, tailMax);

            const OneFit seqSecondary =
                FitOneTemplate(
                    d, c, t,
                    0, 0,
                    tailMin, tailMax,
                    singlePPNN[ie].scale,
                    p);

            sequentialTP[ie].ppnnScale =
                singlePPNN[ie].scale;
            sequentialTP[ie].secondaryScale =
                seqSecondary.scale;
            sequentialTP[ie].ppnnError =
                singlePPNN[ie].error;
            sequentialTP[ie].secondaryError =
                seqSecondary.error;
            sequentialTP[ie].correlation = 0.0;
            sequentialTP[ie].boundary =
                sequentialTP[ie].ppnnScale <= 0.0 ||
                sequentialTP[ie].secondaryScale <= 0.0;
            sequentialTP[ie].activeParameters =
                (sequentialTP[ie].ppnnScale > 0.0 ? 1 : 0) +
                (sequentialTP[ie].secondaryScale > 0.0 ? 1 : 0);

            beforeTP[ie] =
                FitTwoTemplates(
                    d, c, p, t,
                    0, 0,
                    tailMin, tailMax);

            globalTP[ie] =
                FitTwoTemplates(
                    d, c, p, t,
                    0, kNStudyPanels - 1,
                    tailMin, tailMax);

            const ModelMetrics mSingle =
                EvaluateAllPanels(
                    d, c,
                    0, t,
                    0.0,
                    singleTP[ie].scale,
                    compareMin, compareMax,
                    1);

            const ModelMetrics mSequential =
                EvaluateAllPanels(
                    d, c,
                    p, t,
                    sequentialTP[ie].ppnnScale,
                    sequentialTP[ie].secondaryScale,
                    compareMin, compareMax,
                    sequentialTP[ie].activeParameters);

            const ModelMetrics mBefore =
                EvaluateAllPanels(
                    d, c,
                    p, t,
                    beforeTP[ie].ppnnScale,
                    beforeTP[ie].secondaryScale,
                    compareMin, compareMax,
                    beforeTP[ie].activeParameters);

            const ModelMetrics mGlobal =
                EvaluateAllPanels(
                    d, c,
                    p, t,
                    globalTP[ie].ppnnScale,
                    globalTP[ie].secondaryScale,
                    compareMin, compareMax,
                    globalTP[ie].activeParameters);

            WriteMethodLine(
                summary, ie,
                "tp", "single",
                0.0,
                singleTP[ie].scale,
                mSingle,
                0.0,
                singleTP[ie].scale <= 0.0);

            WriteMethodLine(
                summary, ie,
                "tp", "sequential",
                sequentialTP[ie].ppnnScale,
                sequentialTP[ie].secondaryScale,
                mSequential,
                0.0,
                sequentialTP[ie].boundary);

            WriteMethodLine(
                summary, ie,
                "tp", "before2",
                beforeTP[ie].ppnnScale,
                beforeTP[ie].secondaryScale,
                mBefore,
                beforeTP[ie].correlation,
                beforeTP[ie].boundary);

            WriteMethodLine(
                summary, ie,
                "tp", "global2",
                globalTP[ie].ppnnScale,
                globalTP[ie].secondaryScale,
                mGlobal,
                globalTP[ie].correlation,
                globalTP[ie].boundary);
        }

        const ModelMetrics mPPNN =
            EvaluateAllPanels(
                d, c,
                p, 0,
                singlePPNN[ie].scale,
                0.0,
                compareMin, compareMax,
                1);

        WriteMethodLine(
            summary, ie,
            "none", "ppnn_single",
            singlePPNN[ie].scale,
            0.0,
            mPPNN,
            0.0,
            singlePPNN[ie].scale <= 0.0);
    }

    summary
        << "\n# Interpretation flags\n";

    for (int ie = 0; ie < kNE; ++ie) {
        if (he3n[0][ie]) {
            summary
                << "# " << kLabels[ie]
                << " he3n global correlation = "
                << globalHe3n[ie].correlation
                << ", boundary = "
                << (globalHe3n[ie].boundary ? 1 : 0)
                << "\n";
        }

        if (tp[0][ie]) {
            summary
                << "# " << kLabels[ie]
                << " tp global correlation = "
                << globalTP[ie].correlation
                << ", boundary = "
                << (globalTP[ie].boundary ? 1 : 0)
                << "\n";
        }
    }

    TCanvas* cSingle =
        new TCanvas(
            "c_normstudy_single",
            "single-channel study",
            1300, 900);

    cSingle->Print(singlePdf + "[");

    for (int ip = 0; ip < kNStudyPanels; ++ip) {
        cSingle->Clear();
        cSingle->Divide(2, 2, 0.002, 0.002);

        for (int ie = 0; ie < kNE; ++ie) {
            TH1D* modelPPNN =
                BuildModel(
                    coherentDefault[ip][ie],
                    ppnn[ip][ie],
                    0,
                    singlePPNN[ie].scale,
                    0.0,
                    Form(
                        "single_ppnn_p%d_E%d",
                        ip, ie));

            TH1D* modelHe3n =
                he3n[ip][ie]
                ? BuildModel(
                    coherentDefault[ip][ie],
                    0,
                    he3n[ip][ie],
                    0.0,
                    singleHe3n[ie].scale,
                    Form(
                        "single_he3n_p%d_E%d",
                        ip, ie))
                : 0;

            TH1D* modelTP =
                tp[ip][ie]
                ? BuildModel(
                    coherentDefault[ip][ie],
                    0,
                    tp[ip][ie],
                    0.0,
                    singleTP[ie].scale,
                    Form(
                        "single_tp_p%d_E%d",
                        ip, ie))
                : 0;

            StyleData(data[ip][ie]);
            StyleCoherent(coherentDefault[ip][ie]);
            StyleSinglePPNN(modelPPNN);
            StyleSingleSecondary(modelHe3n);

            if (modelTP) {
                modelTP->SetLineColor(kMagenta + 2);
                modelTP->SetLineStyle(7);
                modelTP->SetLineWidth(3);
            }

            cSingle->cd(ie + 1);

            gPad->SetTicks(1, 1);
            gPad->SetLeftMargin(0.13);
            gPad->SetBottomMargin(0.13);

            data[ip][ie]->SetTitle(
                Form(
                    "%s, %s",
                    kLabels[ie],
                    kPanelTitles[ip]));

            data[ip][ie]->GetXaxis()->SetTitle(
                "#Delta E_{#pi^{0}}^{*} (MeV)");

            data[ip][ie]->GetYaxis()->SetTitle(
                "Counts");

            data[ip][ie]->GetXaxis()->SetRangeUser(
                compareMin, compareMax);

            double ymax =
                std::max(
                    data[ip][ie]->GetMaximum(),
                    modelPPNN
                    ? modelPPNN->GetMaximum() : 0.0);

            if (modelHe3n)
                ymax = std::max(
                    ymax,
                    modelHe3n->GetMaximum());

            if (modelTP)
                ymax = std::max(
                    ymax,
                    modelTP->GetMaximum());

            data[ip][ie]->SetMinimum(0.0);
            data[ip][ie]->SetMaximum(
                ymax > 0.0 ? 1.35*ymax : 1.0);

            data[ip][ie]->Draw("E");
            coherentDefault[ip][ie]->Draw("HIST SAME");

            if (modelPPNN)
                modelPPNN->Draw("HIST SAME");

            if (modelHe3n)
                modelHe3n->Draw("HIST SAME");

            if (modelTP)
                modelTP->Draw("HIST SAME");

            data[ip][ie]->Draw("E SAME");

            TLegend* leg =
                new TLegend(
                    0.42, 0.60,
                    0.89, 0.88);

            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.027);

            leg->AddEntry(
                data[ip][ie],
                "Data", "lep");

            leg->AddEntry(
                coherentDefault[ip][ie],
                "Coherent MC", "l");

            if (modelPPNN)
                leg->AddEntry(
                    modelPPNN,
                    "Coherent + ppnn only",
                    "l");

            if (modelHe3n)
                leg->AddEntry(
                    modelHe3n,
                    "Coherent + ^{3}He+n only",
                    "l");

            if (modelTP)
                leg->AddEntry(
                    modelTP,
                    "Coherent + t+p only",
                    "l");

            leg->Draw();

            output->cd();

            if (modelPPNN)
                modelPPNN->Write();

            if (modelHe3n)
                modelHe3n->Write();

            if (modelTP)
                modelTP->Write();
        }

        cSingle->Print(singlePdf);
    }

    cSingle->Print(singlePdf + "]");
    output->cd();
    cSingle->Write();

    TCanvas* cHe3n =
        new TCanvas(
            "c_normstudy_he3n",
            "he3n normalization methods",
            1300, 900);

    cHe3n->Print(he3nPdf + "[");

    for (int ip = 0; ip < kNStudyPanels; ++ip) {
        cHe3n->Clear();
        cHe3n->Divide(2, 2, 0.002, 0.002);

        for (int ie = 0; ie < kNE; ++ie) {
            if (!he3n[ip][ie])
                continue;

            TH1D* modelSeq =
                BuildModel(
                    coherentDefault[ip][ie],
                    ppnn[ip][ie],
                    he3n[ip][ie],
                    sequentialHe3n[ie].ppnnScale,
                    sequentialHe3n[ie].secondaryScale,
                    Form(
                        "he3n_sequential_p%d_E%d",
                        ip, ie));

            TH1D* modelBefore =
                BuildModel(
                    coherentDefault[ip][ie],
                    ppnn[ip][ie],
                    he3n[ip][ie],
                    beforeHe3n[ie].ppnnScale,
                    beforeHe3n[ie].secondaryScale,
                    Form(
                        "he3n_before2_p%d_E%d",
                        ip, ie));

            TH1D* modelGlobal =
                BuildModel(
                    coherentDefault[ip][ie],
                    ppnn[ip][ie],
                    he3n[ip][ie],
                    globalHe3n[ie].ppnnScale,
                    globalHe3n[ie].secondaryScale,
                    Form(
                        "he3n_global2_p%d_E%d",
                        ip, ie));

            StyleData(data[ip][ie]);
            StyleCoherent(coherentDefault[ip][ie]);
            StyleSequential(modelSeq);
            StyleBeforeFit(modelBefore);
            StyleGlobalFit(modelGlobal);

            cHe3n->cd(ie + 1);

            gPad->SetTicks(1, 1);
            gPad->SetLeftMargin(0.13);
            gPad->SetBottomMargin(0.13);

            data[ip][ie]->SetTitle(
                Form(
                    "%s, %s",
                    kLabels[ie],
                    kPanelTitles[ip]));

            data[ip][ie]->GetXaxis()->SetTitle(
                "#Delta E_{#pi^{0}}^{*} (MeV)");

            data[ip][ie]->GetYaxis()->SetTitle(
                "Counts");

            data[ip][ie]->GetXaxis()->SetRangeUser(
                compareMin, compareMax);

            double ymax =
                data[ip][ie]->GetMaximum();

            ymax = std::max(ymax, modelSeq->GetMaximum());
            ymax = std::max(ymax, modelBefore->GetMaximum());
            ymax = std::max(ymax, modelGlobal->GetMaximum());

            data[ip][ie]->SetMinimum(0.0);
            data[ip][ie]->SetMaximum(
                ymax > 0.0 ? 1.35*ymax : 1.0);

            data[ip][ie]->Draw("E");
            coherentDefault[ip][ie]->Draw("HIST SAME");
            modelSeq->Draw("HIST SAME");
            modelBefore->Draw("HIST SAME");
            modelGlobal->Draw("HIST SAME");
            data[ip][ie]->Draw("E SAME");

            TLegend* leg =
                new TLegend(
                    0.41, 0.57,
                    0.89, 0.88);

            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.026);

            leg->AddEntry(
                data[ip][ie],
                "Data", "lep");

            leg->AddEntry(
                coherentDefault[ip][ie],
                "Coherent MC", "l");

            leg->AddEntry(
                modelSeq,
                "Sequential: ppnn then ^{3}He+n",
                "l");

            leg->AddEntry(
                modelBefore,
                "Two-component fit before cut",
                "l");

            leg->AddEntry(
                modelGlobal,
                "Global fit: before + all cuts",
                "l");

            leg->Draw();

            TLatex label;
            label.SetNDC();
            label.SetTextSize(0.026);

            label.DrawLatex(
                0.14, 0.90,
                Form(
                    "global: ppnn %.3g, ^{3}He+n %.3g, #rho %.2f",
                    globalHe3n[ie].ppnnScale,
                    globalHe3n[ie].secondaryScale,
                    globalHe3n[ie].correlation));

            output->cd();

            modelSeq->Write();
            modelBefore->Write();
            modelGlobal->Write();
        }

        cHe3n->Print(he3nPdf);
    }

    cHe3n->Print(he3nPdf + "]");
    output->cd();
    cHe3n->Write();

    bool anyTP = false;

    for (int ie = 0; ie < kNE; ++ie) {
        if (tp[0][ie])
            anyTP = true;
    }

    if (anyTP) {
        TCanvas* cTP =
            new TCanvas(
                "c_normstudy_tp",
                "tp normalization methods",
                1300, 900);

        cTP->Print(tpPdf + "[");

        for (int ip = 0; ip < kNStudyPanels; ++ip) {
            cTP->Clear();
            cTP->Divide(2, 2, 0.002, 0.002);

            for (int ie = 0; ie < kNE; ++ie) {
                cTP->cd(ie + 1);
                gPad->SetTicks(1, 1);
                gPad->SetLeftMargin(0.13);
                gPad->SetBottomMargin(0.13);

                if (!tp[ip][ie]) {
                    TH1D* empty =
                        dynamic_cast<TH1D*>(
                            data[ip][ie]->Clone(
                                Form(
                                    "tp_missing_p%d_E%d",
                                    ip, ie)));

                    empty->Reset("ICES");
                    empty->SetTitle(
                        Form(
                            "%s, %s",
                            kLabels[ie],
                            kPanelTitles[ip]));

                    empty->GetXaxis()->SetTitle(
                        "#Delta E_{#pi^{0}}^{*} (MeV)");

                    empty->GetYaxis()->SetTitle(
                        "Counts");

                    empty->GetXaxis()->SetRangeUser(
                        compareMin, compareMax);

                    empty->SetMinimum(0.0);
                    empty->SetMaximum(1.0);
                    empty->Draw();

                    TLatex missing;
                    missing.SetNDC();
                    missing.SetTextSize(0.045);
                    missing.DrawLatex(
                        0.22, 0.52,
                        "t+p template not available");

                    continue;
                }

                TH1D* modelSeq =
                    BuildModel(
                        coherentDefault[ip][ie],
                        ppnn[ip][ie],
                        tp[ip][ie],
                        sequentialTP[ie].ppnnScale,
                        sequentialTP[ie].secondaryScale,
                        Form(
                            "tp_sequential_p%d_E%d",
                            ip, ie));

                TH1D* modelBefore =
                    BuildModel(
                        coherentDefault[ip][ie],
                        ppnn[ip][ie],
                        tp[ip][ie],
                        beforeTP[ie].ppnnScale,
                        beforeTP[ie].secondaryScale,
                        Form(
                            "tp_before2_p%d_E%d",
                            ip, ie));

                TH1D* modelGlobal =
                    BuildModel(
                        coherentDefault[ip][ie],
                        ppnn[ip][ie],
                        tp[ip][ie],
                        globalTP[ie].ppnnScale,
                        globalTP[ie].secondaryScale,
                        Form(
                            "tp_global2_p%d_E%d",
                            ip, ie));

                StyleData(data[ip][ie]);
                StyleCoherent(coherentDefault[ip][ie]);
                StyleSequential(modelSeq);
                StyleBeforeFit(modelBefore);
                StyleGlobalFit(modelGlobal);

                data[ip][ie]->SetTitle(
                    Form(
                        "%s, %s",
                        kLabels[ie],
                        kPanelTitles[ip]));

                data[ip][ie]->GetXaxis()->SetTitle(
                    "#Delta E_{#pi^{0}}^{*} (MeV)");

                data[ip][ie]->GetYaxis()->SetTitle(
                    "Counts");

                data[ip][ie]->GetXaxis()->SetRangeUser(
                    compareMin, compareMax);

                double ymax =
                    data[ip][ie]->GetMaximum();

                ymax = std::max(
                    ymax, modelSeq->GetMaximum());

                ymax = std::max(
                    ymax, modelBefore->GetMaximum());

                ymax = std::max(
                    ymax, modelGlobal->GetMaximum());

                data[ip][ie]->SetMinimum(0.0);
                data[ip][ie]->SetMaximum(
                    ymax > 0.0 ? 1.35*ymax : 1.0);

                data[ip][ie]->Draw("E");
                coherentDefault[ip][ie]->Draw("HIST SAME");
                modelSeq->Draw("HIST SAME");
                modelBefore->Draw("HIST SAME");
                modelGlobal->Draw("HIST SAME");
                data[ip][ie]->Draw("E SAME");

                TLegend* leg =
                    new TLegend(
                        0.41, 0.57,
                        0.89, 0.88);

                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                leg->SetTextSize(0.026);

                leg->AddEntry(
                    data[ip][ie],
                    "Data", "lep");

                leg->AddEntry(
                    coherentDefault[ip][ie],
                    "Coherent MC", "l");

                leg->AddEntry(
                    modelSeq,
                    "Sequential: ppnn then t+p",
                    "l");

                leg->AddEntry(
                    modelBefore,
                    "Two-component fit before cut",
                    "l");

                leg->AddEntry(
                    modelGlobal,
                    "Global fit: before + all cuts",
                    "l");

                leg->Draw();

                TLatex label;
                label.SetNDC();
                label.SetTextSize(0.026);

                label.DrawLatex(
                    0.14, 0.90,
                    Form(
                        "global: ppnn %.3g, t+p %.3g, #rho %.2f",
                        globalTP[ie].ppnnScale,
                        globalTP[ie].secondaryScale,
                        globalTP[ie].correlation));

                output->cd();

                modelSeq->Write();
                modelBefore->Write();
                modelGlobal->Write();
            }

            cTP->Print(tpPdf);
        }

        cTP->Print(tpPdf + "]");
        output->cd();
        cTP->Write();
    }

    const double scanCoreMin[] = {
        -10.0, -10.0, -5.0
    };

    const double scanCoreMax[] = {
        15.0, 10.0, 15.0
    };

    const double scanTailMin[] = {
        -60.0, -60.0, -55.0, -50.0
    };

    const double scanTailMax[] = {
        -20.0, -15.0, -20.0, -20.0
    };

    const int nCoreScans = 3;
    const int nTailScans = 4;

    for (int ie = 0; ie < kNE; ++ie) {
        TH1D* d[kNStudyPanels];
        TH1D* cUnscaled[kNStudyPanels];
        TH1D* p[kNStudyPanels];
        TH1D* h[kNStudyPanels];
        TH1D* t[kNStudyPanels];

        for (int ip = 0; ip < kNStudyPanels; ++ip) {
            d[ip] = data[ip][ie];
            cUnscaled[ip] =
                coherentUnscaled[ip][ie];
            p[ip] = ppnn[ip][ie];
            h[ip] = he3n[ip][ie];
            t[ip] = tp[ip][ie];
        }

        for (int icore = 0;
             icore < nCoreScans;
             ++icore) {

            TH1D* cScan[kNStudyPanels];

            ScaleCoherentSet(
                d,
                cUnscaled,
                cScan,
                ie,
                scanCoreMin[icore],
                scanCoreMax[icore],
                Form("scan_core%d", icore));

            for (int itail = 0;
                 itail < nTailScans;
                 ++itail) {

                if (h[0]) {
                    const TwoFit fit =
                        FitTwoTemplates(
                            d, cScan, p, h,
                            0, kNStudyPanels - 1,
                            scanTailMin[itail],
                            scanTailMax[itail]);

                    const ModelMetrics metrics =
                        EvaluateAllPanels(
                            d, cScan, p, h,
                            fit.ppnnScale,
                            fit.secondaryScale,
                            compareMin,
                            compareMax,
                            fit.activeParameters);

                    scan
                        << kLabels[ie] << ",he3n,"
                        << scanCoreMin[icore] << ","
                        << scanCoreMax[icore] << ","
                        << scanTailMin[itail] << ","
                        << scanTailMax[itail] << ","
                        << fit.ppnnScale << ","
                        << fit.ppnnError << ","
                        << fit.secondaryScale << ","
                        << fit.secondaryError << ","
                        << fit.correlation << ","
                        << (fit.boundary ? 1 : 0) << ","
                        << fit.chi2 << ","
                        << fit.ndf << ","
                        << metrics.chi2ndf
                        << "\n";
                }

                if (t[0]) {
                    const TwoFit fit =
                        FitTwoTemplates(
                            d, cScan, p, t,
                            0, kNStudyPanels - 1,
                            scanTailMin[itail],
                            scanTailMax[itail]);

                    const ModelMetrics metrics =
                        EvaluateAllPanels(
                            d, cScan, p, t,
                            fit.ppnnScale,
                            fit.secondaryScale,
                            compareMin,
                            compareMax,
                            fit.activeParameters);

                    scan
                        << kLabels[ie] << ",tp,"
                        << scanCoreMin[icore] << ","
                        << scanCoreMax[icore] << ","
                        << scanTailMin[itail] << ","
                        << scanTailMax[itail] << ","
                        << fit.ppnnScale << ","
                        << fit.ppnnError << ","
                        << fit.secondaryScale << ","
                        << fit.secondaryError << ","
                        << fit.correlation << ","
                        << (fit.boundary ? 1 : 0) << ","
                        << fit.chi2 << ","
                        << fit.ndf << ","
                        << metrics.chi2ndf
                        << "\n";
                }
            }

            for (int ip = 0; ip < kNStudyPanels; ++ip)
                delete cScan[ip];
        }
    }

    output->cd();

    for (int ip = 0; ip < kNStudyPanels; ++ip) {
        for (int ie = 0; ie < kNE; ++ie) {
            if (data[ip][ie])
                data[ip][ie]->Write();

            if (coherentDefault[ip][ie])
                coherentDefault[ip][ie]->Write();

            if (ppnn[ip][ie])
                ppnn[ip][ie]->Write();

            if (he3n[ip][ie])
                he3n[ip][ie]->Write();

            if (tp[ip][ie])
                tp[ip][ie]->Write();
        }
    }

    summary.close();
    scan.close();
    output->Close();

    for (int ic = 0; ic < kNCuts; ++ic)
        coherentFiles[ic]->Close();

    fPPNN->Close();
    fHe3n->Close();
    fTP->Close();

    std::cout
        << "\nComplete normalization study finished.\n"
        << "Saved:\n"
        << "  " << singlePdf << "\n"
        << "  " << he3nPdf << "\n";

    if (anyTP)
        std::cout
            << "  " << tpPdf << "\n";

    std::cout
        << "  " << summaryName << "\n"
        << "  " << scanName << "\n"
        << "  " << rootName << "\n";
}
