/*
 * fig2_coherent_lowE_multibin_quantitative_final.C
 *
 * Purpose
 * -------
 * Low-energy coherent-MC closure / applicability study below the four main
 * Fig. 2 photon-energy bins.  Experimental Acqu data are reconstructed in
 * 135-155, 155-175, 175-195 and 195-215 MeV and compared with the dedicated
 * low-energy coherent templates before and after the 8/10/12 degree cuts.
 *
 * Default input data: Acqu_CBTagg_31837.root.
 * Default coherent templates: fig2_coherent_lowE_available_bins_wide.root.
 * Default core/tail/compare regions are [-10,15], [-60,-20], [-60,40] MeV.
 *
 * Usage
 * -----
 *   root -l -b -q 'fig2_coherent_lowE_multibin_quantitative_final.C()'
 *
 * This is a dedicated low-energy study, not the principal Fig. 2 production
 * macro.  The code body below is unchanged from the supplied source.
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
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"

namespace Fig2CoherentLowEFinal {

const int kNCh = 352;
const int kMaxTracks = 512;
const int kMaxTagged = 2048;

const int kNE = 4;
const int kNCuts = 3;

const double kELow[kNE]  = {135.0, 155.0, 175.0, 195.0};
const double kEHigh[kNE] = {155.0, 175.0, 195.0, 215.0};
const char* kLabels[kNE] = {
    "145 MeV", "165 MeV", "185 MeV", "205 MeV"
};

const int kCuts[kNCuts] = {8, 10, 12};

struct Chi2Result {
    double chi2;
    int ndf;
    double pvalue;
};

int FindEnergyBin(double egamma)
{
    for (int i = 0; i < kNE; ++i) {
        const bool inBin =
            egamma >= kELow[i] &&
            (egamma < kEHigh[i] ||
             (i == kNE - 1 && egamma <= kEHigh[i]));
        if (inBin) return i;
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

double IntegralGlobal(TH1D* h)
{
    return h ? h->Integral(0, h->GetNbinsX() + 1) : 0.0;
}

double CoreScale(TH1D* data, TH1D* coherent,
                 double coreMin, double coreMax)
{
    const double d = IntegralRange(data, coreMin, coreMax);
    const double c = IntegralRange(coherent, coreMin, coreMax);
    return c > 0.0 ? d / c : 0.0;
}

Chi2Result EvaluateChi2(TH1D* data, TH1D* model,
                        double xmin, double xmax,
                        int nFittedParameters)
{
    Chi2Result result = {0.0, 0, 0.0};
    if (!data || !model) return result;

    const int b1 = data->GetXaxis()->FindBin(xmin + 1.0e-9);
    const int b2 = data->GetXaxis()->FindBin(xmax - 1.0e-9);

    int used = 0;

    for (int b = b1; b <= b2; ++b) {
        const double d = data->GetBinContent(b);
        const double m = model->GetBinContent(b);
        const double ed = data->GetBinError(b);
        const double em = model->GetBinError(b);
        const double variance = ed * ed + em * em;

        if (variance <= 0.0) continue;

        const double diff = d - m;
        result.chi2 += diff * diff / variance;
        ++used;
    }

    result.ndf = used - nFittedParameters;
    if (result.ndf < 0) result.ndf = 0;

    if (result.ndf > 0)
        result.pvalue = TMath::Prob(result.chi2, result.ndf);

    return result;
}

double TailResidualSigma(TH1D* data, TH1D* model,
                         double tailMin, double tailMax,
                         double& residual,
                         double& residualError)
{
    double ed = 0.0;
    double em = 0.0;

    const double d = IntegralRange(data, tailMin, tailMax, &ed);
    const double m = IntegralRange(model, tailMin, tailMax, &em);

    residual = d - m;
    residualError = std::sqrt(ed * ed + em * em);

    return residualError > 0.0
        ? residual / residualError
        : 0.0;
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
    h->SetLineWidth(3);
    h->SetFillStyle(0);
}

double CoherentPionLabEnergy(double egamma,
                             double thetaCmDeg,
                             double mpi0,
                             double targetMass)
{
    const double s =
        targetMass * targetMass +
        2.0 * egamma * targetMass;

    const double sqrtS = std::sqrt(s);

    const double ePiCm =
        (s + mpi0 * mpi0 -
         targetMass * targetMass) /
        (2.0 * sqrtS);

    double pPiCm2 = ePiCm * ePiCm - mpi0 * mpi0;
    if (pPiCm2 < 0.0) pPiCm2 = 0.0;

    const double pPiCm = std::sqrt(pPiCm2);
    const double betaCm = egamma / (egamma + targetMass);
    const double gammaCm =
        1.0 / std::sqrt(1.0 - betaCm * betaCm);

    const double cosThetaCm =
        std::cos(thetaCmDeg * TMath::DegToRad());

    return gammaCm *
        (ePiCm + betaCm * pPiCm * cosThetaCm);
}

double PhiMinDeg(double ePiLab, double mpi0)
{
    if (ePiLab <= mpi0) return 180.0;

    double arg =
        std::sqrt(ePiLab * ePiLab - mpi0 * mpi0) /
        ePiLab;

    if (arg < -1.0) arg = -1.0;
    if (arg > 1.0) arg = 1.0;

    return 2.0 * std::acos(arg) * TMath::RadToDeg();
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

    double bestDistance = 1.0e30;
    bool found = false;

    bestMgg = -1.0;
    bestOpeningDeg = -1.0;

    for (int i = 0; i < nt; ++i) {
        if (clusterEnergy[i] < clusterEnergyMin) continue;
        if (vetoEnergy[i] > vetoEnergyMax) continue;

        for (int j = i + 1; j < nt; ++j) {
            if (clusterEnergy[j] < clusterEnergyMin) continue;
            if (vetoEnergy[j] > vetoEnergyMax) continue;

            const double th1 = theta[i] * deg;
            const double th2 = theta[j] * deg;
            const double ph1 = phi[i] * deg;
            const double ph2 = phi[j] * deg;

            TLorentzVector gamma1;
            TLorentzVector gamma2;

            gamma1.SetPxPyPzE(
                clusterEnergy[i] * std::sin(th1) * std::cos(ph1),
                clusterEnergy[i] * std::sin(th1) * std::sin(ph1),
                clusterEnergy[i] * std::cos(th1),
                clusterEnergy[i]);

            gamma2.SetPxPyPzE(
                clusterEnergy[j] * std::sin(th2) * std::cos(ph2),
                clusterEnergy[j] * std::sin(th2) * std::sin(ph2),
                clusterEnergy[j] * std::cos(th2),
                clusterEnergy[j]);

            const TLorentzVector candidate = gamma1 + gamma2;
            const double mgg = candidate.M();

            if (mgg <= 0.0) continue;

            const double distance = std::fabs(mgg - mpi0);

            if (distance < bestDistance) {
                bestDistance = distance;
                bestPi0 = candidate;
                bestMgg = mgg;
                bestOpeningDeg =
                    gamma1.Angle(gamma2.Vect()) *
                    TMath::RadToDeg();
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
    if (egamma <= 0.0 || openingDeg < 0.0)
        return false;

    TLorentzVector target(0.0, 0.0, 0.0, targetMass);
    TLorentzVector beam(0.0, 0.0, egamma, egamma);
    TLorentzVector totalCM = beam + target;

    TLorentzVector pi0CM = pi0;
    pi0CM.Boost(-totalCM.BoostVector());

    const double thetaCmDeg =
        pi0CM.Theta() * TMath::RadToDeg();

    const double s =
        targetMass * targetMass +
        2.0 * egamma * targetMass;

    const double sqrtS = std::sqrt(s);

    const double ePiCmExpected =
        (s + mpi0 * mpi0 -
         targetMass * targetMass) /
        (2.0 * sqrtS);

    const double betaCm =
        egamma / (egamma + targetMass);

    const double gammaCm =
        1.0 / std::sqrt(1.0 - betaCm * betaCm);

    const double ePiCmMeasured =
        gammaCm * (pi0.E() - betaCm * pi0.Pz());

    deltaE = ePiCmMeasured - ePiCmExpected;

    const double ePiLabCoherent =
        CoherentPionLabEnergy(
            egamma,
            thetaCmDeg,
            mpi0,
            targetMass);

    deltaPhi =
        openingDeg -
        PhiMinDeg(ePiLabCoherent, mpi0);

    return true;
}

void FillExperimentalData(
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
        std::cerr
            << "Cannot open experimental file "
            << dataFile << std::endl;
        return;
    }

    TTree* tracks =
        dynamic_cast<TTree*>(file->Get("tracks"));

    TTree* tagger =
        dynamic_cast<TTree*>(file->Get("tagger"));

    TTree* setup =
        dynamic_cast<TTree*>(
            file->Get("setupParameters"));

    if (!tracks || !tagger || !setup) {
        std::cerr
            << "Experimental file is missing "
               "tracks, tagger or setupParameters"
            << std::endl;
        file->Close();
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

    Long64_t selectedPairs = 0;
    Long64_t usedTags = 0;

    for (Long64_t event = 0;
         event < entries;
         ++event) {

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

        ++selectedPairs;

        const int nTags =
            std::min(nTagged, kMaxTagged);

        for (int iTag = 0;
             iTag < nTags;
             ++iTag) {

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

            double weight = 0.0;
            const double time = taggedTime[iTag];

            if (time > promptMin &&
                time < promptMax) {
                weight = 1.0;
            } else if (time > randomMin &&
                       time < randomMax) {
                weight = -randomWeight;
            } else {
                continue;
            }

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

            for (int iCut = 0;
                 iCut < kNCuts;
                 ++iCut) {

                if (deltaPhi < kCuts[iCut]) {
                    after[iCut][energyBin]->Fill(
                        deltaE, weight);
                }
            }

            ++usedTags;
        }

        if ((event + 1) % 500000 == 0) {
            std::cout
                << "Data event "
                << event + 1 << " / "
                << entries << std::endl;
        }
    }

    std::cout
        << "Experimental data: processed "
        << entries
        << ", selected pi0 pairs "
        << selectedPairs
        << ", used tag combinations "
        << usedTags
        << std::endl;

    file->Close();
}


TH1D* MapToReference(TH1D* source,
                     TH1D* reference,
                     const char* name)
{
    if (!source || !reference) return 0;

    TH1D* output =
        dynamic_cast<TH1D*>(
            reference->Clone(name));

    if (!output) return 0;

    output->SetDirectory(0);
    output->Reset("ICES");
    output->Sumw2();

    const double xmin =
        output->GetXaxis()->GetXmin();

    const double xmax =
        output->GetXaxis()->GetXmax();

    for (int bin = 1;
         bin <= source->GetNbinsX();
         ++bin) {

        const double x =
            source->GetXaxis()->GetBinCenter(bin);

        if (x < xmin || x >= xmax)
            continue;

        const int outputBin =
            output->GetXaxis()->FindBin(x);

        const double oldContent =
            output->GetBinContent(outputBin);

        const double oldError =
            output->GetBinError(outputBin);

        const double addContent =
            source->GetBinContent(bin);

        const double addError =
            source->GetBinError(bin);

        output->SetBinContent(
            outputBin,
            oldContent + addContent);

        output->SetBinError(
            outputBin,
            std::sqrt(
                oldError * oldError +
                addError * addError));
    }

    return output;
}

void DrawMissingPanel(int eLow,
                      int eHigh,
                      const char* message)
{
    TH1D* frame = new TH1D(
        Form("frame_%d_%d_%p",
             eLow,
             eHigh,
             (void*)gPad),
        "",
        100,
        -60.0,
        40.0);

    frame->SetDirectory(0);
    frame->SetStats(0);
    frame->SetMinimum(0.0);
    frame->SetMaximum(1.0);

    frame->GetXaxis()->SetTitle(
        "#Delta E_{#pi^{0}}^{*} (MeV)");

    frame->GetYaxis()->SetTitle("Counts");
    frame->Draw();

    TLatex text;
    text.SetNDC();
    text.SetTextAlign(22);
    text.SetTextSize(0.055);

    text.DrawLatex(
        0.50,
        0.57,
        Form("%d-%d MeV",
             eLow,
             eHigh));

    text.SetTextSize(0.040);
    text.DrawLatex(
        0.50,
        0.47,
        message);
}

} // namespace Fig2CoherentLowEFinal

void fig2_coherent_lowE_multibin_quantitative_final(
    const char* dataFile =
        "../Acqu_CBTagg_31837.root",
    const char* coherentTemplateFile =
        "fig2_coherent_lowE_available_bins_wide.root",
    Long64_t maxDataEvents = -1,
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
    using namespace Fig2CoherentLowEFinal;

    std::cout
        << "RUNNING coherent low-E multibin "
           "quantitative macro"
        << std::endl;

    gStyle->SetOptStat(0);

    TH1D* dataBefore[kNE];
    TH1D* dataAfter[kNCuts][kNE];

    for (int iEnergy = 0;
         iEnergy < kNE;
         ++iEnergy) {

        dataBefore[iEnergy] =
            new TH1D(
                Form(
                    "data_deltaE_before_phi_E%d",
                    iEnergy),
                Form(
                    "Data %.0f-%.0f MeV;"
                    "#Delta E_{#pi^{0}}^{*} (MeV);"
                    "Prompt-random counts",
                    kELow[iEnergy],
                    kEHigh[iEnergy]),
                240,
                -120.0,
                120.0);

        dataBefore[iEnergy]->Sumw2();

        for (int iCut = 0;
             iCut < kNCuts;
             ++iCut) {

            dataAfter[iCut][iEnergy] =
                new TH1D(
                    Form(
                        "data_deltaE_after_phi_cut%d_E%d",
                        kCuts[iCut],
                        iEnergy),
                    Form(
                        "Data, #Delta#Phi < %d^{#circ}, "
                        "%.0f-%.0f MeV;"
                        "#Delta E_{#pi^{0}}^{*} (MeV);"
                        "Prompt-random counts",
                        kCuts[iCut],
                        kELow[iEnergy],
                        kEHigh[iEnergy]),
                    240,
                    -120.0,
                    120.0);

            dataAfter[iCut][iEnergy]->Sumw2();
        }
    }

    FillExperimentalData(
        dataFile,
        dataBefore,
        dataAfter,
        maxDataEvents,
        clusterEnergyMin,
        vetoEnergyMax,
        mggMin,
        mggMax,
        promptMin,
        promptMax,
        randomMin,
        randomMax);

    TFile* coherentInput =
        TFile::Open(
            coherentTemplateFile,
            "READ");

    if (!coherentInput ||
        coherentInput->IsZombie()) {

        std::cerr
            << "Cannot open coherent template file "
            << coherentTemplateFile
            << std::endl;

        return;
    }

    TFile* output =
        TFile::Open(
            "fig2_coherent_lowE_multibin_quantitative.root",
            "RECREATE");

    if (!output ||
        output->IsZombie()) {

        std::cerr
            << "Cannot create output ROOT file"
            << std::endl;

        coherentInput->Close();
        return;
    }

    std::ofstream textOutput(
        "fig2_coherent_lowE_multibin_quantitative.txt");

    textOutput
        << "# Coherent MC normalized independently "
           "in the core ["
        << coreMin << ","
        << coreMax
        << "] MeV before and after each DeltaPhi cut\n";

    textOutput
        << "# comparison range = ["
        << compareMin << ","
        << compareMax
        << "] MeV\n";

    textOutput
        << "# negative-tail residual in ["
        << tailMin << ","
        << tailMax
        << "] MeV\n";

    textOutput
        << "# cut E_low E_high "
           "coh_scale_before coh_scale_after "
           "Ndata_before Ncoh_before_raw "
           "Ndata_after Ncoh_after_raw "
           "coh_survival_raw "
           "chi2ndf_coh p_coh "
           "tail_residual tail_error tail_sigma\n";

    double coherentScaleBefore[kNE] =
        {0.0, 0.0, 0.0, 0.0};

    double coherentBeforeRawGlobal[kNE] =
        {0.0, 0.0, 0.0, 0.0};

    bool coherentAvailable[kNE] =
        {false, false, false, false};

    TCanvas* beforeCanvas =
        new TCanvas(
            "c_coherent_lowE_before",
            "Low-energy coherent comparison before cut",
            1300,
            900);

    beforeCanvas->Divide(
        2,
        2,
        0.002,
        0.002);

    for (int iEnergy = 0;
         iEnergy < kNE;
         ++iEnergy) {

        TH1D* coherentRaw =
            dynamic_cast<TH1D*>(
                coherentInput->Get(
                    Form(
                        "coh_lowE_deltaE_E%d_nodphi",
                        iEnergy)));

        coherentBeforeRawGlobal[iEnergy] =
            IntegralGlobal(coherentRaw);

        coherentAvailable[iEnergy] =
            coherentRaw &&
            coherentBeforeRawGlobal[iEnergy] > 0.0;

        beforeCanvas->cd(iEnergy + 1);

        gPad->SetTicks(1, 1);
        gPad->SetLeftMargin(0.13);
        gPad->SetBottomMargin(0.13);

        if (!coherentAvailable[iEnergy]) {
            DrawMissingPanel(
                (int)kELow[iEnergy],
                (int)kEHigh[iEnergy],
                "coherent MC not available");
            continue;
        }

        TH1D* data =
            dynamic_cast<TH1D*>(
                dataBefore[iEnergy]->Clone(
                    Form(
                        "data_before_E%d",
                        iEnergy)));

        data->SetDirectory(0);

        TH1D* coherent =
            MapToReference(
                coherentRaw,
                data,
                Form(
                    "coh_before_E%d",
                    iEnergy));

        if (!coherent) {
            DrawMissingPanel(
                (int)kELow[iEnergy],
                (int)kEHigh[iEnergy],
                "cannot map coherent MC");
            continue;
        }

        coherentScaleBefore[iEnergy] =
            CoreScale(
                data,
                coherent,
                coreMin,
                coreMax);

        coherent->Scale(
            coherentScaleBefore[iEnergy]);

        const Chi2Result chi =
            EvaluateChi2(
                data,
                coherent,
                compareMin,
                compareMax,
                1);

        const double chi2ndf =
            chi.ndf > 0
            ? chi.chi2 / chi.ndf
            : 0.0;

        double tailResidual = 0.0;
        double tailError = 0.0;

        const double tailSigma =
            TailResidualSigma(
                data,
                coherent,
                tailMin,
                tailMax,
                tailResidual,
                tailError);

        StyleData(data);
        StyleCoherent(coherent);

        data->SetTitle(
            Form(
                "%s, before #Delta#Phi cut",
                kLabels[iEnergy]));

        data->GetXaxis()->SetTitle(
            "#Delta E_{#pi^{0}}^{*} (MeV)");

        data->GetYaxis()->SetTitle("Counts");

        data->GetXaxis()->SetRangeUser(
            compareMin,
            compareMax);

        const double maximum =
            std::max(
                data->GetMaximum(),
                coherent->GetMaximum());

        data->SetMinimum(0.0);

        data->SetMaximum(
            maximum > 0.0
            ? 1.30 * maximum
            : 1.0);

        data->Draw("E");
        coherent->Draw("HIST SAME");
        data->Draw("E SAME");

        TLegend* legend =
            new TLegend(
                0.47,
                0.69,
                0.89,
                0.88);

        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetTextSize(0.031);

        legend->AddEntry(
            data,
            "Data before cut",
            "lep");

        legend->AddEntry(
            coherent,
            "Coherent MC (core scaled)",
            "l");

        legend->Draw();

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.031);

        latex.DrawLatex(
            0.15,
            0.88,
            Form(
                "coh scale = %.4g",
                coherentScaleBefore[iEnergy]));

        latex.DrawLatex(
            0.15,
            0.83,
            Form(
                "#chi^{2}/ndf = %.2f",
                chi2ndf));

        latex.DrawLatex(
            0.15,
            0.78,
            Form(
                "tail residual = %.1f#sigma",
                tailSigma));

        output->cd();
        data->Write();
        coherent->Write();
    }

    beforeCanvas->SaveAs(
        "fig2_coherent_lowE_multibin_before_normalization.pdf");

    output->cd();
    beforeCanvas->Write();

    TCanvas* afterCanvas =
        new TCanvas(
            "c_coherent_lowE_after",
            "Low-energy coherent comparison after cut",
            1300,
            900);

    afterCanvas->Print(
        "fig2_coherent_lowE_multibin_after_prediction.pdf[");

    for (int iCut = 0;
         iCut < kNCuts;
         ++iCut) {

        afterCanvas->Clear();

        afterCanvas->Divide(
            2,
            2,
            0.002,
            0.002);

        for (int iEnergy = 0;
             iEnergy < kNE;
             ++iEnergy) {

            afterCanvas->cd(iEnergy + 1);

            gPad->SetTicks(1, 1);
            gPad->SetLeftMargin(0.13);
            gPad->SetBottomMargin(0.13);

            if (!coherentAvailable[iEnergy]) {
                DrawMissingPanel(
                    (int)kELow[iEnergy],
                    (int)kEHigh[iEnergy],
                    "coherent MC not available");
                continue;
            }

            TH1D* coherentRaw =
                dynamic_cast<TH1D*>(
                    coherentInput->Get(
                        Form(
                            "coh_lowE_deltaE_E%d_dphi%d",
                            iEnergy,
                            kCuts[iCut])));

            const double coherentAfterRawGlobal =
                IntegralGlobal(coherentRaw);

            if (!coherentRaw ||
                coherentAfterRawGlobal <= 0.0) {

                DrawMissingPanel(
                    (int)kELow[iEnergy],
                    (int)kEHigh[iEnergy],
                    "coherent cut template missing");
                continue;
            }

            TH1D* data =
                dynamic_cast<TH1D*>(
                    dataAfter[iCut][iEnergy]->Clone(
                        Form(
                            "data_after_cut%d_E%d",
                            kCuts[iCut],
                            iEnergy)));

            data->SetDirectory(0);

            TH1D* coherent =
                MapToReference(
                    coherentRaw,
                    data,
                    Form(
                        "coh_after_cut%d_E%d",
                        kCuts[iCut],
                        iEnergy));

            if (!coherent) {
                DrawMissingPanel(
                    (int)kELow[iEnergy],
                    (int)kEHigh[iEnergy],
                    "cannot map coherent MC");
                continue;
            }

            const double coherentScaleAfter =
                CoreScale(
                    data,
                    coherent,
                    coreMin,
                    coreMax);

            coherent->Scale(
                coherentScaleAfter);

            const Chi2Result chi =
                EvaluateChi2(
                    data,
                    coherent,
                    compareMin,
                    compareMax,
                    1);

            const double chi2ndf =
                chi.ndf > 0
                ? chi.chi2 / chi.ndf
                : 0.0;

            double tailResidual = 0.0;
            double tailError = 0.0;

            const double tailSigma =
                TailResidualSigma(
                    data,
                    coherent,
                    tailMin,
                    tailMax,
                    tailResidual,
                    tailError);

            const double coherentSurvival =
                coherentBeforeRawGlobal[iEnergy] > 0.0
                ? coherentAfterRawGlobal /
                  coherentBeforeRawGlobal[iEnergy]
                : 0.0;

            textOutput
                << kCuts[iCut] << " "
                << kELow[iEnergy] << " "
                << kEHigh[iEnergy] << " "
                << coherentScaleBefore[iEnergy] << " "
                << coherentScaleAfter << " "
                << IntegralGlobal(
                       dataBefore[iEnergy]) << " "
                << coherentBeforeRawGlobal[iEnergy] << " "
                << IntegralGlobal(
                       dataAfter[iCut][iEnergy]) << " "
                << coherentAfterRawGlobal << " "
                << coherentSurvival << " "
                << chi2ndf << " "
                << chi.pvalue << " "
                << tailResidual << " "
                << tailError << " "
                << tailSigma << "\n";

            StyleData(data);
            StyleCoherent(coherent);

            data->SetTitle(
                Form(
                    "%s, #Delta#Phi < %d^{#circ}",
                    kLabels[iEnergy],
                    kCuts[iCut]));

            data->GetXaxis()->SetTitle(
                "#Delta E_{#pi^{0}}^{*} (MeV)");

            data->GetYaxis()->SetTitle("Counts");

            data->GetXaxis()->SetRangeUser(
                compareMin,
                compareMax);

            const double maximum =
                std::max(
                    data->GetMaximum(),
                    coherent->GetMaximum());

            data->SetMinimum(0.0);

            data->SetMaximum(
                maximum > 0.0
                ? 1.32 * maximum
                : 1.0);

            data->Draw("E");
            coherent->Draw("HIST SAME");
            data->Draw("E SAME");

            TLegend* legend =
                new TLegend(
                    0.47,
                    0.71,
                    0.89,
                    0.88);

            legend->SetBorderSize(0);
            legend->SetFillStyle(0);
            legend->SetTextSize(0.031);

            legend->AddEntry(
                data,
                "Data after cut",
                "lep");

            legend->AddEntry(
                coherent,
                "Coherent MC (core scaled)",
                "l");

            legend->Draw();

            TLatex latex;
            latex.SetNDC();
            latex.SetTextSize(0.030);

            latex.DrawLatex(
                0.15,
                0.89,
                Form(
                    "coh scale: before %.4g, after %.4g",
                    coherentScaleBefore[iEnergy],
                    coherentScaleAfter));

            latex.DrawLatex(
                0.15,
                0.84,
                Form(
                    "#chi^{2}/ndf = %.2f",
                    chi2ndf));

            latex.DrawLatex(
                0.15,
                0.79,
                Form(
                    "tail residual = %.1f#sigma",
                    tailSigma));

            output->cd();
            data->Write();
            coherent->Write();
        }

        afterCanvas->Print(
            "fig2_coherent_lowE_multibin_after_prediction.pdf");
    }

    afterCanvas->Print(
        "fig2_coherent_lowE_multibin_after_prediction.pdf]");

    output->cd();

    for (int iEnergy = 0;
         iEnergy < kNE;
         ++iEnergy) {

        dataBefore[iEnergy]->Write(
            Form(
                "raw_data_before_E%d",
                iEnergy));

        for (int iCut = 0;
             iCut < kNCuts;
             ++iCut) {

            dataAfter[iCut][iEnergy]->Write(
                Form(
                    "raw_data_after_cut%d_E%d",
                    kCuts[iCut],
                    iEnergy));
        }
    }

    afterCanvas->Write();

    textOutput.close();
    output->Close();
    coherentInput->Close();

    std::cout
        << "Saved:\n"
        << "  fig2_coherent_lowE_multibin_before_normalization.pdf\n"
        << "  fig2_coherent_lowE_multibin_after_prediction.pdf\n"
        << "  fig2_coherent_lowE_multibin_quantitative.root\n"
        << "  fig2_coherent_lowE_multibin_quantitative.txt\n";
}
