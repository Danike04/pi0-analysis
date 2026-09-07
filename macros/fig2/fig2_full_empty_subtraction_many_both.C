/*
 * fig2_full_empty_subtraction_many_both.C
 *
 * PURPOSE
 *   Apply the same pi0 / missing-energy selection to many FULL-target and
 *   many EMPTY-target AcquRoot files, normalize the EMPTY contribution to the
 *   FULL photon flux, and save FULL - scaled(EMPTY) histograms for the four
 *   Fig. 2 photon-energy intervals and the DeltaPhi = 8, 10, 12 degree cuts.
 *
 * INPUT FILE LISTS
 *   fullListFile  : text file with one FULL-target ROOT path per line.
 *   emptyListFile : text file with one EMPTY-target ROOT path per line.
 *   Empty lines and lines beginning with '#' are ignored.
 *
 * EVENT SELECTION
 *   The event selection is implemented inside this macro and uses the
 *   tracks, tagger and setupParameters trees. The default values are:
 *     clusterEnergy > 20 MeV
 *     vetoEnergy <= 1 MeV
 *     110 < m(gamma gamma) < 155 MeV
 *     prompt: 700-800 ns
 *     random: 450-680 ns
 *     DeltaPhi cuts: 8, 10 and 12 degrees
 *
 * PHOTON-FLUX NORMALIZATION
 *   For each of the four energy intervals the macro integrates tagger
 *   electron scalers over all listed files and converts them to photon flux
 *   using the supplied tagging-efficiency file. The FPD configuration file
 *   is used to map tagger channels to scaler indices.
 *
 *   The real FPD map is required for production; the macro aborts if an
 *   explicitly supplied map cannot be read. A constant tagging efficiency
 *   is used only when no tagging-efficiency filename is requested.
 *
 *   The subtraction factor is calculated independently for each energy bin:
 *       alpha = Phi_full / Phi_empty
 *       H_sub = H_full - alpha * H_empty
 *   The uncertainty on alpha is propagated to the EMPTY-scaled and
 *   subtracted histogram errors.
 *
 * OUTPUT HISTOGRAM INTERFACE
 *   The FULL-EMPTY histograms used by the later refit are written as:
 *     he4_before_E0 ... he4_before_E3
 *     he4_cut8_E0  ... he4_cut8_E3
 *     he4_cut10_E0 ... he4_cut10_E3
 *     he4_cut12_E0 ... he4_cut12_E3
 *
 * OUTPUT FILES
 *   Names are generated from outputTag:
 *     fig2_full_empty_<tag>.root
 *     fig2_full_empty_<tag>.txt
 *     fig2_full_empty_<tag>_before.pdf
 *     fig2_full_empty_<tag>_after.pdf
 *
 * INPUT PATHS
 *   The default auxiliary filenames are relative to the current working
 *   directory. For production use, pass explicit paths to the validated FPD
 *   map and zero-based tagging-efficiency table.
 *
 * NOTE FOR THE REPOSITORY
 *   This is the version that accepts lists for BOTH FULL and EMPTY datasets.
 *   Final repository review also aligned the tagging-efficiency parser with
 *   the validated zero-based production table and removed the unsafe FPD
 *   fallback from the production normalization path.
 */

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TChain.h"
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
    const std::vector<std::string>& dataFiles,
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

    TChain* tracks = new TChain("tracks");
    TChain* tagger = new TChain("tagger");
    TChain* setup  = new TChain("setupParameters");

    if (dataFiles.empty()) {
        std::cerr << "No data files supplied." << std::endl;
        delete tracks;
        delete tagger;
        delete setup;
        return false;
    }

    int tracksAdded = 0;
    int taggerAdded = 0;
    int setupAdded = 0;

    for (std::size_t iFile = 0; iFile < dataFiles.size(); ++iFile) {
        tracksAdded += tracks->Add(dataFiles[iFile].c_str());
        taggerAdded += tagger->Add(dataFiles[iFile].c_str());
        setupAdded  += setup->Add(dataFiles[iFile].c_str());
    }

    const int expectedFiles = static_cast<int>(dataFiles.size());

    if (tracksAdded != expectedFiles ||
        taggerAdded != expectedFiles ||
        setupAdded != expectedFiles) {
        std::cerr
            << "Cannot add all required FULL-target trees. "
            << "Expected " << expectedFiles
            << " files, added tracks/tagger/setup = "
            << tracksAdded << "/"
            << taggerAdded << "/"
            << setupAdded
            << std::endl;
        delete tracks;
        delete tagger;
        delete setup;
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

    delete tracks;
    delete tagger;
    delete setup;
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


namespace FullEmptyAnalysis {

using namespace Fig2FromAcqu;

const int kScalerArraySize = 8706;

bool LoadFPDScalerMap(const char* fpdFile, int scalerIndex[], int nCh)
{
    for (int ch = 0; ch < nCh; ++ch)
        scalerIndex[ch] = -1;

    std::ifstream input(fpdFile);
    if (!input.is_open())
        return false;

    std::string line;
    int loaded = 0;

    while (std::getline(input, line)) {
        if (line.find("Element:") != 0)
            continue;

        std::istringstream stream(line);
        std::string label;
        int channel = -1;

        stream >> label >> channel;
        if (channel < 0 || channel >= nCh)
            continue;

        std::vector<std::string> fields;
        std::string token;

        while (stream >> token)
            fields.push_back(token);

        if (fields.empty())
            continue;

        const int index =
            std::atoi(fields.back().c_str());

        if (index < 0 || index >= kScalerArraySize)
            continue;

        scalerIndex[channel] = index;
        ++loaded;
    }

    return loaded > 0;
}

bool LoadTaggingEfficiency(
    const char* fileName,
    double efficiency[],
    double uncertainty[],
    int nCh)
{
    for (int ch = 0; ch < nCh; ++ch) {
        efficiency[ch] = -1.0;
        uncertainty[ch] = 0.0;
    }

    if (!fileName || std::string(fileName).empty())
        return false;

    std::ifstream input(fileName);
    if (!input.is_open())
        return false;

    std::string line;
    int loaded = 0;

    while (std::getline(input, line)) {
        if (line.empty() || line[0] == '#')
            continue;

        std::istringstream stream(line);

        int fileChannel = -1;
        double value = 0.0;
        double error = 0.0;

        if (!(stream >> fileChannel >> value))
            continue;

        if (!(stream >> error))
            error = 0.0;

        // Production tagging-efficiency tables use zero-based tagger channels.
        const int channel = fileChannel;

        if (channel < 0 || channel >= nCh)
            continue;

        efficiency[channel] = value;
        uncertainty[channel] = error;
        ++loaded;
    }

    return loaded > 0;
}

bool ReadTaggerSetup(
    const char* fileName,
    int& nValidChannels,
    double photonEnergy[kNCh])
{
    nValidChannels = 0;
    for (int ch = 0; ch < kNCh; ++ch)
        photonEnergy[ch] = 0.0;

    TFile file(fileName, "READ");
    if (file.IsZombie()) {
        std::cerr << "Cannot open " << fileName << std::endl;
        return false;
    }

    TTree* setup =
        dynamic_cast<TTree*>(file.Get("setupParameters"));

    if (!setup) {
        std::cerr << "Missing setupParameters in "
                  << fileName << std::endl;
        return false;
    }

    int nTagger = 0;
    setup->SetBranchAddress("nTagger", &nTagger);
    setup->SetBranchAddress(
        "TaggerPhotonEnergy", photonEnergy);

    if (setup->GetEntries() <= 0) {
        std::cerr << "Empty setupParameters in "
                  << fileName << std::endl;
        return false;
    }

    setup->GetEntry(0);
    nValidChannels = std::min(nTagger, kNCh);

    return nValidChannels > 0;
}

bool FillExperimentalDataOne(
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
        (promptMax - promptMin) /
        (randomMax - randomMin);

    TFile* file = TFile::Open(dataFile);

    if (!file || file->IsZombie()) {
        std::cerr << "Cannot open " << dataFile << std::endl;
        return false;
    }

    TTree* tracks =
        dynamic_cast<TTree*>(file->Get("tracks"));
    TTree* tagger =
        dynamic_cast<TTree*>(file->Get("tagger"));
    TTree* setup =
        dynamic_cast<TTree*>(
            file->Get("setupParameters"));

    if (!tracks || !tagger || !setup) {
        std::cerr << "Missing trees in "
                  << dataFile << std::endl;
        file->Close();
        delete file;
        return false;
    }

    int nTracks = 0;
    double clusterEnergy[kMaxTracks];
    double theta[kMaxTracks];
    double phi[kMaxTracks];
    double vetoEnergy[kMaxTracks];

    tracks->SetBranchAddress("nTracks", &nTracks);
    tracks->SetBranchAddress(
        "clusterEnergy", clusterEnergy);
    tracks->SetBranchAddress("theta", theta);
    tracks->SetBranchAddress("phi", phi);
    tracks->SetBranchAddress(
        "vetoEnergy", vetoEnergy);

    int nTagged = 0;
    int taggedChannel[kMaxTagged];
    double taggedTime[kMaxTagged];

    tagger->SetBranchAddress("nTagged", &nTagged);
    tagger->SetBranchAddress(
        "taggedChannel", taggedChannel);
    tagger->SetBranchAddress(
        "taggedTime", taggedTime);

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

    for (Long64_t event = 0;
         event < entries; ++event) {

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

        for (int iTag = 0;
             iTag < nTags; ++iTag) {

            const int channel =
                taggedChannel[iTag];

            if (channel < 0 ||
                channel >= nValidChannels)
                continue;

            const double egamma =
                taggerPhotonEnergy[channel];

            const int energyBin =
                FindEnergyBin(egamma);

            if (energyBin < 0)
                continue;

            const double time =
                taggedTime[iTag];

            double weight = 0.0;

            if (time > promptMin &&
                time < promptMax)
                weight = 1.0;
            else if (time > randomMin &&
                     time < randomMax)
                weight = -randomWeight;
            else
                continue;

            double deltaE = 0.0;
            double deltaPhi = 0.0;

            if (!CalculateKinematics(
                    pi0,
                    openingDeg,
                    egamma,
                    targetMass,
                    mpi0,
                    deltaE,
                    deltaPhi))
                continue;

            before[energyBin]->Fill(
                deltaE, weight);

            for (int ic = 0;
                 ic < kNCuts; ++ic) {
                if (deltaPhi < kCuts[ic])
                    after[ic][energyBin]->Fill(
                        deltaE, weight);
            }
        }

        if ((event + 1) % 500000 == 0)
            std::cout << "EMPTY event "
                      << event + 1 << " / "
                      << entries << std::endl;
    }

    file->Close();
    delete file;
    return true;
}

struct FluxResult {
    double photonFlux[kNE];
    double electronCounts[kNE];
    double photonFluxVariance[kNE];
    Long64_t scalerEntries;
    int validChannels;
    bool usedTaggingEfficiencyFile;
    bool usedFPDMap;
};

FluxResult CalculateFlux(
    const std::vector<std::string>& files,
    const char* fpdFile,
    const char* taggingEfficiencyFile,
    double fallbackTaggingEfficiency,
    double minimumTaggingEfficiency,
    bool skipFirstScalerEntry)
{
    FluxResult result;

    for (int ie = 0; ie < kNE; ++ie) {
        result.photonFlux[ie] = 0.0;
        result.electronCounts[ie] = 0.0;
        result.photonFluxVariance[ie] = 0.0;
    }

    result.scalerEntries = 0;
    result.validChannels = 0;
    result.usedTaggingEfficiencyFile = false;
    result.usedFPDMap = false;

    if (files.empty())
        return result;

    double photonEnergy[kNCh];
    int nValidChannels = 0;

    if (!ReadTaggerSetup(
            files.front().c_str(),
            nValidChannels,
            photonEnergy))
        return result;

    result.validChannels = nValidChannels;

    int scalerIndex[kNCh];

    const bool haveFPDMap =
        LoadFPDScalerMap(
            fpdFile,
            scalerIndex,
            nValidChannels);

    result.usedFPDMap = haveFPDMap;

    if (!haveFPDMap) {
        std::cerr
            << "ERROR: cannot read required FPD map from "
            << fpdFile << "\n"
            << "FULL/EMPTY flux normalization requires the real "
            << "non-contiguous FPD mapping; aborting."
            << std::endl;
        return result;
    }

    double taggingEfficiency[kNCh];
    double taggingEfficiencyError[kNCh];

    const bool haveTaggingEfficiency =
        LoadTaggingEfficiency(
            taggingEfficiencyFile,
            taggingEfficiency,
            taggingEfficiencyError,
            nValidChannels);

    result.usedTaggingEfficiencyFile =
        haveTaggingEfficiency;

    if (!haveTaggingEfficiency &&
        taggingEfficiencyFile &&
        std::string(taggingEfficiencyFile).size() > 0) {
        std::cerr
            << "ERROR: tagging-efficiency file was requested but cannot be read: "
            << taggingEfficiencyFile << std::endl;
        return result;
    }

    if (!haveTaggingEfficiency) {
        std::cerr
            << "No tagging-efficiency file requested; using constant epsilon_tag = "
            << fallbackTaggingEfficiency
            << ". The same value is used for "
            << "full and empty."
            << std::endl;
    }

    double electronByChannel[kNCh];
    double electronVarianceByChannel[kNCh];

    for (int ch = 0;
         ch < nValidChannels; ++ch) {
        electronByChannel[ch] = 0.0;
        electronVarianceByChannel[ch] = 0.0;
    }

    for (std::size_t iFile = 0;
         iFile < files.size(); ++iFile) {

        TFile file(files[iFile].c_str(), "READ");

        if (file.IsZombie()) {
            std::cerr << "Cannot open "
                      << files[iFile]
                      << " for scaler integration."
                      << std::endl;
            continue;
        }

        TTree* scalers =
            dynamic_cast<TTree*>(
                file.Get("scalers"));

        if (!scalers) {
            std::cerr << "Missing scalers in "
                      << files[iFile]
                      << std::endl;
            continue;
        }

        UInt_t scalerValues[kScalerArraySize];
        scalers->SetBranchAddress(
            "scalers", scalerValues);

        const Long64_t entries =
            scalers->GetEntries();

        result.scalerEntries += entries;

        for (Long64_t entry = 0;
             entry < entries; ++entry) {

            scalers->GetEntry(entry);

            if (skipFirstScalerEntry &&
                entry == 0)
                continue;

            for (int ch = 0;
                 ch < nValidChannels; ++ch) {

                const int index =
                    scalerIndex[ch];

                if (index < 0 ||
                    index >= kScalerArraySize)
                    continue;

                const double value =
                    scalerValues[index];

                electronByChannel[ch] += value;
                electronVarianceByChannel[ch] += value;
            }
        }
    }

    for (int ch = 0;
         ch < nValidChannels; ++ch) {

        const int energyBin =
            FindEnergyBin(photonEnergy[ch]);

        if (energyBin < 0)
            continue;

        double efficiency =
            fallbackTaggingEfficiency;
        double efficiencyError = 0.0;

        if (haveTaggingEfficiency &&
            taggingEfficiency[ch] >=
                minimumTaggingEfficiency) {
            efficiency =
                taggingEfficiency[ch];
            efficiencyError =
                taggingEfficiencyError[ch];
        }

        const double electrons =
            electronByChannel[ch];

        result.electronCounts[energyBin] +=
            electrons;

        result.photonFlux[energyBin] +=
            electrons * efficiency;

        // Poisson scaler uncertainty plus the
        // quoted tagging-efficiency uncertainty.
        result.photonFluxVariance[energyBin] +=
            electronVarianceByChannel[ch] *
                efficiency * efficiency +
            electrons * electrons *
                efficiencyError *
                efficiencyError;
    }

    return result;
}

TH1D* MakeHistogram(
    const char* name)
{
    TH1D* histogram =
        new TH1D(name, "", 240, -120.0, 120.0);

    histogram->SetDirectory(0);
    histogram->Sumw2();

    return histogram;
}

TH1D* CloneHistogram(
    TH1D* source,
    const char* name)
{
    if (!source)
        return 0;

    TH1D* clone =
        dynamic_cast<TH1D*>(
            source->Clone(name));

    if (!clone)
        return 0;

    clone->SetDirectory(0);
    return clone;
}

void StyleFull(TH1D* histogram)
{
    histogram->SetLineColor(kBlack);
    histogram->SetMarkerColor(kBlack);
    histogram->SetMarkerStyle(20);
    histogram->SetMarkerSize(0.42);
    histogram->SetLineWidth(2);
}

void StyleEmpty(TH1D* histogram)
{
    histogram->SetLineColor(kRed + 1);
    histogram->SetMarkerColor(kRed + 1);
    histogram->SetMarkerStyle(24);
    histogram->SetMarkerSize(0.40);
    histogram->SetLineWidth(2);
}

void StyleSubtracted(TH1D* histogram)
{
    histogram->SetLineColor(kBlue + 1);
    histogram->SetMarkerColor(kBlue + 1);
    histogram->SetMarkerStyle(21);
    histogram->SetMarkerSize(0.40);
    histogram->SetLineWidth(2);
}

void PrepareAxes(
    TH1D* full,
    TH1D* empty,
    TH1D* subtracted,
    const char* title)
{
    if (!full || !empty || !subtracted)
        return;

    full->SetTitle(title);
    full->GetXaxis()->SetTitle(
        "#Delta E_{#pi^{0}}^{*} (MeV)");
    full->GetYaxis()->SetTitle(
        "Counts at full-target exposure");
    full->GetXaxis()->SetRangeUser(-60.0, 40.0);
    full->GetYaxis()->SetTitleOffset(1.30);

    double maximum =
        std::max(full->GetMaximum(),
        std::max(empty->GetMaximum(),
                 subtracted->GetMaximum()));

    double minimum =
        std::min(full->GetMinimum(),
        std::min(empty->GetMinimum(),
                 subtracted->GetMinimum()));

    if (maximum <= 0.0)
        maximum = 1.0;

    full->SetMaximum(1.20 * maximum);

    if (minimum < 0.0)
        full->SetMinimum(1.20 * minimum);
    else
        full->SetMinimum(0.0);
}

void DrawPanel(
    TH1D* full,
    TH1D* emptyScaled,
    TH1D* subtracted,
    const char* title,
    double fullFlux,
    double emptyFlux,
    double alpha)
{
    PrepareAxes(
        full, emptyScaled, subtracted, title);

    StyleFull(full);
    StyleEmpty(emptyScaled);
    StyleSubtracted(subtracted);

    full->Draw("E1");
    emptyScaled->Draw("E1 SAME");
    subtracted->Draw("E1 SAME");

    TLegend legend(0.53, 0.65, 0.88, 0.87);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetTextSize(0.032);
    legend.AddEntry(full, "Full target", "lep");
    legend.AddEntry(
        emptyScaled,
        "Empty target, flux-scaled",
        "lep");
    legend.AddEntry(
        subtracted,
        "Full - empty", "lep");
    legend.DrawClone();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.030);
    latex.DrawLatex(
        0.15, 0.86,
        Form("#alpha = #Phi_{full}/#Phi_{empty} = %.3g",
             alpha));

    latex.SetTextSize(0.026);
    latex.DrawLatex(
        0.15, 0.81,
        Form("#Phi_{full}=%.4e, #Phi_{empty}=%.4e",
             fullFlux, emptyFlux));
}

void BuildScaledAndSubtracted(
    TH1D* full,
    TH1D* emptyRaw,
    TH1D* emptyScaled,
    TH1D* subtracted,
    double alpha,
    double alphaError)
{
    if (!full || !emptyRaw ||
        !emptyScaled || !subtracted)
        return;

    const int nBins = full->GetNbinsX();

    for (int bin = 0; bin <= nBins + 1; ++bin) {
        const double f = full->GetBinContent(bin);
        const double ef = full->GetBinError(bin);
        const double e = emptyRaw->GetBinContent(bin);
        const double ee = emptyRaw->GetBinError(bin);

        const double scaled = alpha * e;
        const double scaledVariance =
            alpha * alpha * ee * ee +
            e * e * alphaError * alphaError;

        emptyScaled->SetBinContent(bin, scaled);
        emptyScaled->SetBinError(
            bin, std::sqrt(std::max(0.0, scaledVariance)));

        subtracted->SetBinContent(bin, f - scaled);
        subtracted->SetBinError(
            bin,
            std::sqrt(std::max(
                0.0, ef * ef + scaledVariance)));
    }
}


} // namespace FullEmptyAnalysis


void fig2_full_empty_subtraction_many_both(
    const char* fullListFile =
        "full99_good_files.txt",
    const char* emptyListFile =
        "empty5_files.txt",
    const char* outputTag =
        "full99_minus_empty5",
    const char* fpdFile =
        "FPD_855_new.dat",
    const char* taggingEfficiencyFile =
        "ExpBkgSub_COPP_TaggEff_31834.dat",
    double fallbackTaggingEfficiency = 0.20,
    double minimumTaggingEfficiency = 0.02,
    bool skipFirstScalerEntry = true,
    Long64_t maxFullEvents = -1,
    Long64_t maxEmptyEvents = -1,
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
    using namespace FullEmptyAnalysis;

    gStyle->SetOptStat(0);

    TString tag(outputTag);
    tag.ReplaceAll("/", "_");
    tag.ReplaceAll(" ", "_");

    std::vector<std::string> fullFiles;
    {
        std::ifstream input(fullListFile);
        if (!input) {
            std::cerr
                << "Cannot open FULL-target file list "
                << fullListFile << std::endl;
            return;
        }

        std::string line;
        while (std::getline(input, line)) {
            if (line.empty() || line[0] == '#')
                continue;
            fullFiles.push_back(line);
        }
    }

    if (fullFiles.empty()) {
        std::cerr
            << "FULL-target file list is empty: "
            << fullListFile << std::endl;
        return;
    }

    std::cout
        << "Loaded " << fullFiles.size()
        << " FULL-target files from "
        << fullListFile << std::endl;

    std::vector<std::string> emptyFiles;
    {
        std::ifstream input(emptyListFile);
        if (!input) {
            std::cerr
                << "Cannot open EMPTY-target file list "
                << emptyListFile << std::endl;
            return;
        }

        std::string line;
        while (std::getline(input, line)) {
            if (line.empty() || line[0] == '#')
                continue;
            emptyFiles.push_back(line);
        }
    }

    if (emptyFiles.empty()) {
        std::cerr
            << "EMPTY-target file list is empty: "
            << emptyListFile << std::endl;
        return;
    }

    std::cout
        << "Loaded " << emptyFiles.size()
        << " EMPTY-target files from "
        << emptyListFile << std::endl;

    TH1D* fullBefore[kNE];
    TH1D* emptyBefore[kNE];
    TH1D* fullAfter[kNCuts][kNE];
    TH1D* emptyAfter[kNCuts][kNE];

    for (int ie = 0; ie < kNE; ++ie) {
        fullBefore[ie] =
            MakeHistogram(
                Form("full_before_E%d", ie));

        emptyBefore[ie] =
            MakeHistogram(
                Form("empty_before_E%d", ie));

        for (int ic = 0;
             ic < kNCuts; ++ic) {
            fullAfter[ic][ie] =
                MakeHistogram(
                    Form("full_cut%d_E%d",
                         kCuts[ic], ie));

            emptyAfter[ic][ie] =
                MakeHistogram(
                    Form("empty_cut%d_E%d",
                         kCuts[ic], ie));
        }
    }

    std::cout
        << "\n=== Filling FULL target from "
        << fullFiles.size() << " files ===\n"
        << std::endl;

    if (!FillExperimentalData(
            fullFiles,
            fullBefore,
            fullAfter,
            maxFullEvents,
            clusterEnergyMin,
            vetoEnergyMax,
            mggMin,
            mggMax,
            promptMin,
            promptMax,
            randomMin,
            randomMax))
        return;

    std::cout
        << "\n=== Filling EMPTY target from "
        << emptyFiles.size() << " files ===\n"
        << std::endl;

    if (!FillExperimentalData(
            emptyFiles,
            emptyBefore,
            emptyAfter,
            maxEmptyEvents,
            clusterEnergyMin,
            vetoEnergyMax,
            mggMin,
            mggMax,
            promptMin,
            promptMax,
            randomMin,
            randomMax))
        return;

    std::cout
        << "\n=== Integrating FULL photon flux ===\n";

    const FluxResult fullFlux =
        CalculateFlux(
            fullFiles,
            fpdFile,
            taggingEfficiencyFile,
            fallbackTaggingEfficiency,
            minimumTaggingEfficiency,
            skipFirstScalerEntry);

    std::cout
        << "\n=== Integrating EMPTY photon flux ===\n";

    const FluxResult emptyFlux =
        CalculateFlux(
            emptyFiles,
            fpdFile,
            taggingEfficiencyFile,
            fallbackTaggingEfficiency,
            minimumTaggingEfficiency,
            skipFirstScalerEntry);

    TString rootName =
        Form("fig2_full_empty_%s.root",
             tag.Data());

    TString textName =
        Form("fig2_full_empty_%s.txt",
             tag.Data());

    TString beforePdf =
        Form("fig2_full_empty_%s_before.pdf",
             tag.Data());

    TString afterPdf =
        Form("fig2_full_empty_%s_after.pdf",
             tag.Data());

    TFile* output =
        TFile::Open(rootName, "RECREATE");

    if (!output || output->IsZombie()) {
        std::cerr << "Cannot create "
                  << rootName << std::endl;
        return;
    }

    std::ofstream textOutput(textName.Data());

    textOutput
        << "# Full-empty subtraction\n"
        << "# Full file list: " << fullListFile << "\n"
        << "# Number of full files: " << fullFiles.size() << "\n"
        << "# Empty file list: " << emptyListFile << "\n"
        << "# Number of empty files: " << emptyFiles.size() << "\n"
        << "# FPD map: " << fpdFile << "\n"
        << "# Tagging efficiency: "
        << taggingEfficiencyFile << "\n"
        << "# skipFirstScalerEntry = "
        << skipFirstScalerEntry << "\n"
        << "# full scaler entries = "
        << fullFlux.scalerEntries << "\n"
        << "# empty scaler entries = "
        << emptyFlux.scalerEntries << "\n"
        << "# ie label panel Phi_full Phi_empty alpha alpha_error "
        << "full_int empty_raw_int empty_scaled_int "
        << "subtracted_int subtracted_error\n";

    TH1D* scaledEmptyBefore[kNE];
    TH1D* subtractedBefore[kNE];
    TH1D* scaledEmptyAfter[kNCuts][kNE];
    TH1D* subtractedAfter[kNCuts][kNE];

    double alpha[kNE];
    double alphaError[kNE];

    for (int ie = 0; ie < kNE; ++ie) {
        alpha[ie] = 0.0;
        alphaError[ie] = 0.0;

        if (fullFlux.photonFlux[ie] <= 0.0 ||
            emptyFlux.photonFlux[ie] <= 0.0) {
            std::cerr
                << "Invalid photon flux in "
                << kLabels[ie]
                << ": full="
                << fullFlux.photonFlux[ie]
                << ", empty="
                << emptyFlux.photonFlux[ie]
                << std::endl;
            output->Close();
            return;
        }

        alpha[ie] =
            fullFlux.photonFlux[ie] /
            emptyFlux.photonFlux[ie];

        const double relativeVariance =
            fullFlux.photonFluxVariance[ie] /
                (fullFlux.photonFlux[ie] *
                 fullFlux.photonFlux[ie]) +
            emptyFlux.photonFluxVariance[ie] /
                (emptyFlux.photonFlux[ie] *
                 emptyFlux.photonFlux[ie]);

        alphaError[ie] =
            alpha[ie] *
            std::sqrt(std::max(0.0, relativeVariance));

        scaledEmptyBefore[ie] =
            CloneHistogram(
                emptyBefore[ie],
                Form("empty_scaled_before_E%d", ie));

        subtractedBefore[ie] =
            CloneHistogram(
                fullBefore[ie],
                Form("he4_before_E%d", ie));

        BuildScaledAndSubtracted(
            fullBefore[ie],
            emptyBefore[ie],
            scaledEmptyBefore[ie],
            subtractedBefore[ie],
            alpha[ie],
            alphaError[ie]);

        for (int ic = 0;
             ic < kNCuts; ++ic) {

            scaledEmptyAfter[ic][ie] =
                CloneHistogram(
                    emptyAfter[ic][ie],
                    Form("empty_scaled_cut%d_E%d",
                         kCuts[ic], ie));

            subtractedAfter[ic][ie] =
                CloneHistogram(
                    fullAfter[ic][ie],
                    Form("he4_cut%d_E%d",
                         kCuts[ic], ie));

            BuildScaledAndSubtracted(
                fullAfter[ic][ie],
                emptyAfter[ic][ie],
                scaledEmptyAfter[ic][ie],
                subtractedAfter[ic][ie],
                alpha[ie],
                alphaError[ie]);
        }
    }

    TCanvas* canvas =
        new TCanvas(
            "c_full_empty_before",
            "full minus empty before cut",
            1400, 950);

    canvas->Divide(2, 2, 0.002, 0.002);

    for (int ie = 0; ie < kNE; ++ie) {
        canvas->cd(ie + 1);
        gPad->SetTicks(1, 1);
        gPad->SetLeftMargin(0.12);
        gPad->SetBottomMargin(0.12);

        DrawPanel(
            fullBefore[ie],
            scaledEmptyBefore[ie],
            subtractedBefore[ie],
            Form("%s, before #Delta#Phi cut",
                 kLabels[ie]),
            fullFlux.photonFlux[ie],
            emptyFlux.photonFlux[ie],
            alpha[ie]);

        double subtractedError = 0.0;

        const double fullIntegral =
            IntegralRange(
                fullBefore[ie],
                -60.0, 40.0);

        const double emptyRawIntegral =
            IntegralRange(
                emptyBefore[ie],
                -60.0, 40.0);

        const double emptyScaledIntegral =
            IntegralRange(
                scaledEmptyBefore[ie],
                -60.0, 40.0);

        const double subtractedIntegral =
            IntegralRange(
                subtractedBefore[ie],
                -60.0, 40.0,
                &subtractedError);

        textOutput
            << ie << " \"" << kLabels[ie]
            << "\" before "
            << fullFlux.photonFlux[ie] << " "
            << emptyFlux.photonFlux[ie] << " "
            << alpha[ie] << " "
            << alphaError[ie] << " "
            << fullIntegral << " "
            << emptyRawIntegral << " "
            << emptyScaledIntegral << " "
            << subtractedIntegral << " "
            << subtractedError << "\n";

        output->cd();
        fullBefore[ie]->Write();
        emptyBefore[ie]->Write();
        scaledEmptyBefore[ie]->Write();
        subtractedBefore[ie]->Write();
    }

    canvas->SaveAs(beforePdf);

    TCanvas* afterCanvas =
        new TCanvas(
            "c_full_empty_after",
            "full minus empty after cuts",
            1400, 950);

    afterCanvas->Print(afterPdf + "[");

    for (int ic = 0;
         ic < kNCuts; ++ic) {

        afterCanvas->Clear();
        afterCanvas->Divide(
            2, 2, 0.002, 0.002);

        for (int ie = 0;
             ie < kNE; ++ie) {

            afterCanvas->cd(ie + 1);
            gPad->SetTicks(1, 1);
            gPad->SetLeftMargin(0.12);
            gPad->SetBottomMargin(0.12);

            DrawPanel(
                fullAfter[ic][ie],
                scaledEmptyAfter[ic][ie],
                subtractedAfter[ic][ie],
                Form("%s, #Delta#Phi < %d^{#circ}",
                     kLabels[ie], kCuts[ic]),
                fullFlux.photonFlux[ie],
                emptyFlux.photonFlux[ie],
                alpha[ie]);

            double subtractedError = 0.0;

            const double fullIntegral =
                IntegralRange(
                    fullAfter[ic][ie],
                    -60.0, 40.0);

            const double emptyRawIntegral =
                IntegralRange(
                    emptyAfter[ic][ie],
                    -60.0, 40.0);

            const double emptyScaledIntegral =
                IntegralRange(
                    scaledEmptyAfter[ic][ie],
                    -60.0, 40.0);

            const double subtractedIntegral =
                IntegralRange(
                    subtractedAfter[ic][ie],
                    -60.0, 40.0,
                    &subtractedError);

            textOutput
                << ie << " \"" << kLabels[ie]
                << "\" cut" << kCuts[ic] << " "
                << fullFlux.photonFlux[ie] << " "
                << emptyFlux.photonFlux[ie] << " "
                << alpha[ie] << " "
                << alphaError[ie] << " "
                << fullIntegral << " "
                << emptyRawIntegral << " "
                << emptyScaledIntegral << " "
                << subtractedIntegral << " "
                << subtractedError << "\n";

            output->cd();
            fullAfter[ic][ie]->Write();
            emptyAfter[ic][ie]->Write();
            scaledEmptyAfter[ic][ie]->Write();
            subtractedAfter[ic][ie]->Write();
        }

        afterCanvas->Print(afterPdf);
    }

    afterCanvas->Print(afterPdf + "]");

    output->cd();
    canvas->Write();
    afterCanvas->Write();

    output->Close();
    textOutput.close();

    std::cout
        << "\n=== FLUX SUMMARY ===\n";

    for (int ie = 0; ie < kNE; ++ie) {
        std::cout
            << kLabels[ie]
            << ": Phi_full="
            << fullFlux.photonFlux[ie]
            << ", Phi_empty="
            << emptyFlux.photonFlux[ie]
            << ", alpha="
            << alpha[ie]
            << " +/- "
            << alphaError[ie]
            << std::endl;
    }

    std::cout
        << "\nCreated:\n"
        << "  " << beforePdf << "\n"
        << "  " << afterPdf << "\n"
        << "  " << rootName << "\n"
        << "  " << textName << "\n"
        << std::endl;
}
