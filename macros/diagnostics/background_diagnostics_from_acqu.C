/*
 * background_diagnostics_from_acqu.C
 *
 * Purpose
 * -------
 * Extended diagnostic of the experimental Acqu pi0 selection. The macro
 * reconstructs the best two-neutral-cluster pi0 candidate and studies the
 * correlations among invariant mass, missing mass, photon energy, pion
 * angles and coherent missing energy DeltaE_pi.
 *
 * Main default selections
 * -----------------------
 * neutral cluster energy >= 20 MeV
 * veto energy <= 1 MeV
 * 110 <= m(gamma gamma) <= 155 MeV
 * prompt tagger time: 700--800 ns
 * random tagger time: 450--680 ns
 * baseline missing mass: -10--40 MeV relative to the He4 mass
 *
 * The prompt-random subtraction uses the ratio of prompt/random time-window
 * widths. The macro also scans several missing-mass / DeltaE_pi selections
 * and a fitted diagonal band in the (MM, DeltaE_pi) plane.
 *
 * Inputs / normalization
 * ----------------------
 * The interface also contains the same scaler, tagging-efficiency,
 * detection-efficiency and target-thickness inputs used by the historical
 * cross-section checks. See the function arguments below for defaults.
 *
 * Important warning
 * -----------------
 * This is a DIAGNOSTIC macro, although it also calculates a cross-section
 * check. It writes files named `sigma_from_acqu.txt`, `sigma_from_acqu.pdf`
 * and `sigma_from_acqu_vs_egamma.pdf`, which can overwrite outputs produced
 * by `cross_section_from_acqu.C` if both macros are run in the same working
 * directory. Run this macro in a separate output directory.
 *
 * Representative outputs
 * ----------------------
 * background_cut_stability.txt
 * background_diagonal_band_scan.txt
 * background_diagnostics_from_acqu.root
 * mgg_best_from_acqu.pdf
 * missing_mass_before_after_from_acqu.pdf
 * missing_mass_vs_delta_epi_before_after_from_acqu.pdf
 * plus additional 2D/angle diagnostic PDFs.
 *
 * Usage example
 * -------------
 * root -l -b -q 'background_diagnostics_from_acqu.C("Acqu_CBTagg_31837.root")'
 *
 * The code below is the original analysis code. Only this documentation
 * header has been added; calculations, defaults and output names are unchanged.
 */

#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TLegend.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TProfile.h"
#include "TString.h"
#include "TTree.h"

bool load_fpd_scaler_map(const char* fpdFile, int scalerIndex[], int nCh)
{
    for(int ch=0; ch<nCh; ch++) scalerIndex[ch] = -1;

    std::ifstream in(fpdFile);
    if(!in.is_open()) return false;

    std::string line;
    while(std::getline(in, line)) {
        if(line.find("Element:") != 0) continue;

        std::istringstream ss(line);
        std::string label;
        int ch = -1;
        ss >> label >> ch;
        if(ch < 0 || ch >= nCh) continue;

        std::vector<std::string> fields;
        std::string tok;
        while(ss >> tok) fields.push_back(tok);
        if(fields.empty()) continue;

        scalerIndex[ch] = std::atoi(fields.back().c_str());
    }

    return true;
}

TGraphErrors* load_sigma_reference(const char* refFile,
                                   const char* graphName,
                                   int color,
                                   int marker,
                                   double scale=1.0)
{
    if(!refFile || std::string(refFile).empty()) return 0;

    std::ifstream in(refFile);
    if(!in.is_open()) {
        std::cout << "Warning: could not open reference file " << refFile << std::endl;
        return 0;
    }

    TGraphErrors* gr = new TGraphErrors();
    gr->SetName(graphName);
    gr->SetTitle(graphName);
    gr->SetLineColor(color);
    gr->SetMarkerColor(color);
    gr->SetMarkerStyle(marker);
    gr->SetMarkerSize(1.0);
    gr->SetLineWidth(1);

    std::string line;
    int ip = 0;
    while(std::getline(in, line)) {
        if(line.empty() || line[0] == '#') continue;

        std::istringstream ss(line);
        double e = 0.0;
        double sigma = 0.0;
        double dsigma = 0.0;

        if(!(ss >> e >> sigma)) continue;
        if(!(ss >> dsigma)) dsigma = 0.0;

        gr->SetPoint(ip, e, sigma * scale);
        gr->SetPointError(ip, 0.0, dsigma * scale);
        ip++;
    }

    if(ip == 0) {
        delete gr;
        return 0;
    }

    return gr;
}

bool load_eps_det_file(const char* epsFile,
                       double epsDetByCh[],
                       double dEpsDetByCh[],
                       int nCh)
{
    for(int ch=0; ch<nCh; ch++) {
        epsDetByCh[ch] = -1.0;
        dEpsDetByCh[ch] = 0.0;
    }

    if(!epsFile || std::string(epsFile).empty()) return false;

    std::ifstream in(epsFile);
    if(!in.is_open()) {
        std::cout << "Warning: could not open eps_det file " << epsFile << std::endl;
        return false;
    }

    std::string line;
    while(std::getline(in, line)) {
        if(line.empty() || line[0] == '#') continue;

        std::istringstream ss(line);

        int ch = -1;
        double egamma = 0.0;
        double ngen = 0.0;
        double ndet = 0.0;
        double nreco = 0.0;
        double epsDet = 0.0;
        double dEpsDet = 0.0;
        double epsReco = 0.0;
        double dEpsReco = 0.0;

        if(!(ss >> ch >> egamma >> ngen >> ndet >> nreco >> epsDet >> dEpsDet >> epsReco >> dEpsReco)) {
            continue;
        }

        if(ch < 0 || ch >= nCh) continue;

        // Use the reconstructed pi0 efficiency from eps_det_reco_fast_from_acqu_geant.txt.
        epsDetByCh[ch] = epsReco;
        dEpsDetByCh[ch] = dEpsReco;
    }

    return true;
}

bool load_tagging_eff_file(const char* taggEffFile,
                           double epsTagByCh[],
                           double dEpsTagByCh[],
                           int nCh)
{
    for(int ch=0; ch<nCh; ch++) {
        epsTagByCh[ch] = -1.0;
        dEpsTagByCh[ch] = 0.0;
    }

    if(!taggEffFile || std::string(taggEffFile).empty()) return false;

    std::ifstream in(taggEffFile);
    if(!in.is_open()) {
        std::cout << "Warning: could not open tagging efficiency file " << taggEffFile << std::endl;
        return false;
    }

    std::string line;
    while(std::getline(in, line)) {
        if(line.empty() || line[0] == '#') continue;

        std::istringstream ss(line);
        int chFile = -1;
        double eff = 0.0;
        double deff = 0.0;

        if(!(ss >> chFile >> eff >> deff)) continue;

        // The tagging-efficiency table used by this analysis is already
        // zero-based and uses the same channel numbering as the Acqu trees.
        int ch = chFile;
        if(ch < 0 || ch >= nCh) continue;

        epsTagByCh[ch] = eff;
        dEpsTagByCh[ch] = deff;
    }

    return true;
}

TString path_next_to_input(const char* inputFile, const char* neighborFile)
{
    TString in(inputFile);
    Ssiz_t slash = in.Last('/');
    if(slash == kNPOS) return TString(neighborFile);

    TString dir = in(0, slash + 1);
    return dir + neighborFile;
}

void background_diagnostics_from_acqu(const char* fname="Acqu_CBTagg_31837.root",
                             Long64_t maxEvents=-1,
                             double eps_tag=0.30,
                             double eps_det=0.7064,
                             double thickness=0.940e-7,
                             bool skipFirstScalerEntry=true,
                             double minScalerSum=1.0e7,
                             const char* fpdFile="FPD_855_new.dat",
                             double promptMin=700.0,
                             double promptMax=800.0,
                             double randomMin=450.0,
                             double randomMax=680.0,
                             double mmMin=-10.0,
                             double mmMax=40.0,
                             double mggMin=110.0,   //110
                             double mggMax=155.0,   //155
                             const char* referenceFile="",
                             double referenceScale=1.0,
                             const char* epsDetFile="",
                             double minEpsDet=0.02,
                             const char* taggEffFile="ExpBkgSub_COPP_TaggEff_31834.dat",
                             double minEpsTag=0.02)
{
    const int nCh = 352;
    const int maxTracks = 512;
    const int maxTagged = 2048;

    const double beamE = 855.0;       // MeV
    const double targetMass = 3727.38; // MeV, 4He
    const double mpi0 = 134.9768;

    const double eMin = 20.0;
    const double vetoMax = 1.0;
    const double randomWeight = (promptMax - promptMin) / (randomMax - randomMin);

    TFile* f = TFile::Open(fname);
    if(!f || f->IsZombie()) {
        std::cout << "Cannot open " << fname << std::endl;
        return;
    }

    TTree* tracks = (TTree*)f->Get("tracks");
    TTree* tagger = (TTree*)f->Get("tagger");
    TTree* scalers = (TTree*)f->Get("scalers");
    TTree* setup = (TTree*)f->Get("setupParameters");

    if(!tracks || !tagger || !scalers || !setup) {
        std::cout << "Missing one of the required trees: tracks, tagger, scalers, setupParameters" << std::endl;
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
    double taggedTime[maxTagged];

    tagger->SetBranchAddress("nTagged", &nTagged);
    tagger->SetBranchAddress("taggedChannel", taggedChannel);
    tagger->SetBranchAddress("taggedTime", taggedTime);

    int nTagger = 0;
    double taggerPhotonEnergy[nCh];
    setup->SetBranchAddress("nTagger", &nTagger);
    setup->SetBranchAddress("TaggerPhotonEnergy", taggerPhotonEnergy);
    setup->GetEntry(0);
    int nValidCh = std::min(nTagger, nCh);

    UInt_t sc[8706];
    scalers->SetBranchAddress("scalers", sc);

    int scalerIndex[nCh];
    bool haveFpdMap = load_fpd_scaler_map(fpdFile, scalerIndex, nCh);
    if(!haveFpdMap) {
        std::cout << "Warning: could not read FPD scaler map from " << fpdFile << std::endl;
        std::cout << "Falling back to scalerIndex[ch] = 2000 + ch" << std::endl;
        for(int ch=0; ch<nCh; ch++) scalerIndex[ch] = 2000 + ch;
    }

    TString epsDetPath(epsDetFile);
    if(epsDetPath.Length() == 0) {
        epsDetPath = path_next_to_input(fname, "eps_det_reco_fast_from_acqu_geant.txt");
    }

    TString referencePath(referenceFile);
    if(referencePath.Length() == 0) {
        referencePath = path_next_to_input(fname, "paper_sigma_4he_pi0.txt");
    }

    double epsDetByCh[nCh];
    double dEpsDetByCh[nCh];
    bool haveEpsDetFile = load_eps_det_file(epsDetPath.Data(), epsDetByCh, dEpsDetByCh, nCh);
    if(haveEpsDetFile) {
        std::cout << "Using channel-dependent eps_det from " << epsDetPath << std::endl;
    } else {
        std::cout << "Using constant eps_det = " << eps_det << std::endl;
    }

    double epsTagByCh[nCh];
    double dEpsTagByCh[nCh];
    bool haveTaggEffFile = load_tagging_eff_file(taggEffFile, epsTagByCh, dEpsTagByCh, nCh);
    if(haveTaggEffFile) {
        std::cout << "Using channel-dependent eps_tag from " << taggEffFile << std::endl;
    } else {
        std::cout << "Using constant eps_tag = " << eps_tag << std::endl;
    }

    double Ne[nCh] = {0.0};
    double dNe2[nCh] = {0.0};
    for(Long64_t is=0; is<scalers->GetEntries(); is++) {
        scalers->GetEntry(is);
        if(skipFirstScalerEntry && is == 0) continue;

        for(int ch=0; ch<nCh; ch++) {
            if(scalerIndex[ch] < 0 || scalerIndex[ch] >= 8706) continue;
            double val = sc[scalerIndex[ch]];
            Ne[ch] += val;
            dNe2[ch] += val;
        }
    }

    TH1D* hMMp[nCh];
    TH1D* hMMr[nCh];
    for(int ch=0; ch<nCh; ch++) {
        hMMp[ch] = new TH1D(Form("hMM_prompt_ch%d", ch),
                            Form("prompt MM channel %d;MM-M_{^{4}He} (MeV);Counts", ch),
                            250, -100.0, 150.0);
        hMMr[ch] = new TH1D(Form("hMM_random_ch%d", ch),
                            Form("random MM channel %d;MM-M_{^{4}He} (MeV);Counts", ch),
                            250, -100.0, 150.0);
    }

    TH1D* hMggBest = new TH1D("h_mgg_best_selected",
                              "Best #gamma#gamma pair before mass cut;m_{#gamma#gamma} (MeV);Counts",
                              400, 0.0, 400.0);

    TH1D* hMMBefore = new TH1D("h_mm_before_mgg_cut",
                               "Missing mass before m_{#gamma#gamma} cut;MM-M_{^{4}He} (MeV);Prompt - w Random",
                               250, -100.0, 150.0);
    TH1D* hMMAfter = new TH1D("h_mm_after_mgg_cut",
                              "Missing mass after m_{#gamma#gamma} cut;MM-M_{^{4}He} (MeV);Prompt - w Random",
                              250, -100.0, 150.0);

    TH2D* hMggMMBefore = new TH2D("h_mgg_vs_mm_before_mgg_cut",
                                  "m_{#gamma#gamma} vs missing mass before m_{#gamma#gamma} cut;m_{#gamma#gamma} (MeV);MM-M_{^{4}He} (MeV)",
                                  250, 0.0, 250.0,
                                  250, -100.0, 150.0);
    TH2D* hMggMMAfter = new TH2D("h_mgg_vs_mm_after_mgg_cut",
                                 "m_{#gamma#gamma} vs missing mass after m_{#gamma#gamma} cut;m_{#gamma#gamma} (MeV);MM-M_{^{4}He} (MeV)",
                                 250, 0.0, 250.0,
                                 250, -100.0, 150.0);
    TH2D* hMMEgammaBefore = new TH2D("h_mm_vs_egamma_before_mgg_cut",
                                     "Missing mass vs E_{#gamma} before m_{#gamma#gamma} cut;E_{#gamma} (MeV);MM-M_{^{4}He} (MeV)",
                                     352, 0.0, beamE,
                                     250, -100.0, 150.0);
    TH2D* hMMEgammaAfter = new TH2D("h_mm_vs_egamma_after_mgg_cut",
                                    "Missing mass vs E_{#gamma} after m_{#gamma#gamma} cut;E_{#gamma} (MeV);MM-M_{^{4}He} (MeV)",
                                    352, 0.0, beamE,
                                    250, -100.0, 150.0);
    TH2D* hEgammaThetaBefore = new TH2D("h_egamma_vs_theta_before_mgg_cut",
                                        "E_{#gamma} vs #theta_{#pi^{0}} before m_{#gamma#gamma} cut;#theta_{#pi^{0}} (deg);E_{#gamma} (MeV)",
                                        180, 0.0, 180.0,
                                        352, 0.0, beamE);
    TH2D* hEgammaThetaAfter = new TH2D("h_egamma_vs_theta_after_mgg_cut",
                                       "E_{#gamma} vs #theta_{#pi^{0}} after m_{#gamma#gamma} cut;#theta_{#pi^{0}} (deg);E_{#gamma} (MeV)",
                                       180, 0.0, 180.0,
                                       352, 0.0, beamE);

    TH1D* hDeltaEpiBefore = new TH1D("h_delta_epi_before_mgg_cut",
                                     "#DeltaE_{#pi}^{CM} before m_{#gamma#gamma} cut;E_{#pi}^{CM,meas} - E_{#pi}^{CM,coh} (MeV);Prompt - w Random",
                                     300, -150.0, 150.0);
    TH1D* hDeltaEpiAfter = new TH1D("h_delta_epi_after_mgg_cut",
                                    "#DeltaE_{#pi}^{CM} after m_{#gamma#gamma} cut;E_{#pi}^{CM,meas} - E_{#pi}^{CM,coh} (MeV);Prompt - w Random",
                                    300, -150.0, 150.0);
    TH2D* hMMDeltaBefore = new TH2D("h_mm_vs_delta_epi_before_mgg_cut",
                                    "Missing mass vs #DeltaE_{#pi}^{CM} before m_{#gamma#gamma} cut;#DeltaE_{#pi}^{CM} (MeV);MM-M_{^{4}He} (MeV)",
                                    300, -150.0, 150.0,
                                    250, -100.0, 150.0);
    TH2D* hMMDeltaAfter = new TH2D("h_mm_vs_delta_epi_after_mgg_cut",
                                   "Missing mass vs #DeltaE_{#pi}^{CM} after m_{#gamma#gamma} cut;#DeltaE_{#pi}^{CM} (MeV);MM-M_{^{4}He} (MeV)",
                                   300, -150.0, 150.0,
                                   250, -100.0, 150.0);
    TH2D* hMMDeltaPromptAfter = new TH2D("h_mm_vs_delta_epi_prompt_after_mgg_cut",
                                         "Prompt missing mass vs #DeltaE_{#pi}^{CM} after m_{#gamma#gamma} cut;#DeltaE_{#pi}^{CM} (MeV);MM-M_{^{4}He} (MeV)",
                                         300, -150.0, 150.0,
                                         250, -100.0, 150.0);
    TH2D* hMggDeltaBefore = new TH2D("h_mgg_vs_delta_epi_before_mgg_cut",
                                     "m_{#gamma#gamma} vs #DeltaE_{#pi}^{CM};m_{#gamma#gamma} (MeV);#DeltaE_{#pi}^{CM} (MeV)",
                                     250, 0.0, 250.0,
                                     300, -150.0, 150.0);
    TH2D* hEgammaDeltaAfter = new TH2D("h_egamma_vs_delta_epi_after_mgg_cut",
                                       "E_{#gamma} vs #DeltaE_{#pi}^{CM} after m_{#gamma#gamma} cut;E_{#gamma} (MeV);#DeltaE_{#pi}^{CM} (MeV)",
                                       352, 0.0, beamE,
                                       300, -150.0, 150.0);
    TH2D* hThetaCmDeltaAfter = new TH2D("h_theta_cm_vs_delta_epi_after_mgg_cut",
                                        "#theta_{CM} vs #DeltaE_{#pi}^{CM} after m_{#gamma#gamma} cut;#theta_{CM} (deg);#DeltaE_{#pi}^{CM} (MeV)",
                                        180, 0.0, 180.0,
                                        300, -150.0, 150.0);
    TH2D* hThetaLabCmAfter = new TH2D("h_theta_lab_vs_theta_cm_after_mgg_cut",
                                      "#theta_{lab} vs #theta_{CM} after m_{#gamma#gamma} cut;#theta_{CM} (deg);#theta_{lab} (deg)",
                                      180, 0.0, 180.0,
                                      180, 0.0, 180.0);
    TH1D* hBandResidualAfter = new TH1D("h_diagonal_band_residual_after_mgg_cut",
                                        "Distance from fitted MM-#DeltaE_{#pi}^{CM} band after m_{#gamma#gamma} cut;perpendicular distance from band (MeV);Prompt - w Random",
                                        300, -150.0, 150.0);

    Long64_t nentries = std::min(tracks->GetEntries(), tagger->GetEntries());
    if(maxEvents > 0 && maxEvents < nentries) nentries = maxEvents;

    const double deg = TMath::DegToRad();
    TLorentzVector target(0.0, 0.0, 0.0, targetMass);

    Long64_t nPairEvents = 0;
    Long64_t nPi0Events = 0;
    Long64_t nPromptTagsBefore = 0;
    Long64_t nRandomTagsBefore = 0;
    Long64_t nPromptTags = 0;
    Long64_t nRandomTags = 0;

    const int nMmScan = 5;
    const double mmLowScan[nMmScan] = {-20.0, -15.0, -10.0, -5.0, 0.0};
    const double mmHighScan[nMmScan] = {50.0, 45.0, 40.0, 35.0, 30.0};
    const int nDeScan = 5;
    const double deHalfScan[nDeScan] = {40.0, 30.0, 20.0, 15.0, 10.0};
    double cutYield[nMmScan][nDeScan] = {};
    double cutErr2[nMmScan][nDeScan] = {};
    double cutPrompt[nMmScan][nDeScan] = {};
    double cutRandom[nMmScan][nDeScan] = {};
    const int nBandScan = 6;
    const double bandHalfScan[nBandScan] = {50.0, 40.0, 30.0, 25.0, 20.0, 15.0};
    double bandYield[nBandScan] = {};
    double bandErr2[nBandScan] = {};
    double bandPrompt[nBandScan] = {};
    double bandRandom[nBandScan] = {};

    for(Long64_t ev=0; ev<nentries; ev++) {
        tracks->GetEntry(ev);
        tagger->GetEntry(ev);

        int nt = std::min(nTracks, maxTracks);
        int bestI = -1;
        int bestJ = -1;
        double bestM = -1.0;
        double bestDist = 1.0e30;
        TLorentzVector bestPi0;

        for(int i=0; i<nt; i++) {
            if(clusterEnergy[i] < eMin) continue;
            if(vetoEnergy[i] > vetoMax) continue;

            for(int j=i+1; j<nt; j++) {
                if(clusterEnergy[j] < eMin) continue;
                if(vetoEnergy[j] > vetoMax) continue;

                double th1 = theta[i] * deg;
                double th2 = theta[j] * deg;
                double ph1 = phi[i] * deg;
                double ph2 = phi[j] * deg;

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

                TLorentzVector pi0 = g1 + g2;
                double mgg = pi0.M();
                if(mgg <= 0.0) continue;

                double dist = std::fabs(mgg - mpi0);
                if(dist < bestDist) {
                    bestDist = dist;
                    bestM = mgg;
                    bestI = i;
                    bestJ = j;
                    bestPi0 = pi0;
                }
            }
        }

        if(bestI < 0 || bestJ < 0) continue;
        nPairEvents++;
        hMggBest->Fill(bestM);

        bool passMggCut = (bestM >= mggMin && bestM <= mggMax);
        if(passMggCut) nPi0Events++;
        double thetaPi0Deg = bestPi0.Theta() * TMath::RadToDeg();

        int ntag = std::min(nTagged, maxTagged);
        for(int it=0; it<ntag; it++) {
            int ch = taggedChannel[it];
            if(ch < 0 || ch >= nCh) continue;

            double egamma = taggerPhotonEnergy[ch];
            if(egamma <= 0.0 || egamma > beamE + 5.0) continue;

            TLorentzVector beam(0.0, 0.0, egamma, egamma);
            TLorentzVector recoil = beam + target - bestPi0;
            double mm = recoil.M() - targetMass;
            TLorentzVector total = beam + target;
            TLorentzVector pi0CM = bestPi0;
            pi0CM.Boost(-total.BoostVector());
            double wcm = total.M();
            // Coherent two-body reference: gamma + 4He -> pi0 + 4He.
            double ePiCohCM = (wcm*wcm + mpi0*mpi0 - targetMass*targetMass) / (2.0*wcm);
            double deltaEpi = pi0CM.E() - ePiCohCM;
            double thetaCmDeg = pi0CM.Theta() * TMath::RadToDeg();

            double tt = taggedTime[it];
            double weight = 0.0;
            if(tt > promptMin && tt < promptMax) {
                nPromptTagsBefore++;
                weight = 1.0;
            } else if(tt > randomMin && tt < randomMax) {
                nRandomTagsBefore++;
                weight = -randomWeight;
            } else {
                continue;
            }

            hMMBefore->Fill(mm, weight);
            hMggMMBefore->Fill(bestM, mm, weight);
            hMMEgammaBefore->Fill(egamma, mm, weight);
            hEgammaThetaBefore->Fill(thetaPi0Deg, egamma, weight);
            hDeltaEpiBefore->Fill(deltaEpi, weight);
            hMMDeltaBefore->Fill(deltaEpi, mm, weight);
            hMggDeltaBefore->Fill(bestM, deltaEpi, weight);

            if(!passMggCut) continue;

            hMMAfter->Fill(mm, weight);
            hMggMMAfter->Fill(bestM, mm, weight);
            hMMEgammaAfter->Fill(egamma, mm, weight);
            hEgammaThetaAfter->Fill(thetaPi0Deg, egamma, weight);
            hDeltaEpiAfter->Fill(deltaEpi, weight);
            hMMDeltaAfter->Fill(deltaEpi, mm, weight);
            if(weight > 0.0) hMMDeltaPromptAfter->Fill(deltaEpi, mm);
            hEgammaDeltaAfter->Fill(egamma, deltaEpi, weight);
            hThetaCmDeltaAfter->Fill(thetaCmDeg, deltaEpi, weight);
            hThetaLabCmAfter->Fill(thetaCmDeg, thetaPi0Deg, weight);

            for(int im=0; im<nMmScan; im++) {
                if(mm < mmLowScan[im] || mm > mmHighScan[im]) continue;
                for(int id=0; id<nDeScan; id++) {
                    if(std::fabs(deltaEpi) > deHalfScan[id]) continue;
                    cutYield[im][id] += weight;
                    cutErr2[im][id] += weight * weight;
                    if(weight > 0.0) cutPrompt[im][id] += 1.0;
                    else cutRandom[im][id] += 1.0;
                }
            }

            if(weight > 0.0) {
                nPromptTags++;
                hMMp[ch]->Fill(mm);
            } else {
                nRandomTags++;
                hMMr[ch]->Fill(mm);
            }
        }
    }

    TProfile* pMMDeltaBand = hMMDeltaPromptAfter->ProfileX("p_mm_vs_delta_epi_prompt_after_mgg_cut",
                                                           1, -1, "s");
    TF1* fMMDeltaBand = new TF1("f_mm_delta_epi_band",
                                "pol1",
                                -40.0,
                                40.0);
    pMMDeltaBand->Fit(fMMDeltaBand, "QNR");
    double bandP0 = fMMDeltaBand->GetParameter(0);
    double bandP1 = fMMDeltaBand->GetParameter(1);
    if(!std::isfinite(bandP0) || !std::isfinite(bandP1)) {
        bandP0 = 0.0;
        bandP1 = 0.0;
    }

    for(Long64_t ev=0; ev<nentries; ev++) {
        tracks->GetEntry(ev);
        tagger->GetEntry(ev);

        int nt = std::min(nTracks, maxTracks);
        int bestI = -1;
        int bestJ = -1;
        double bestM = -1.0;
        double bestDist = 1.0e30;
        TLorentzVector bestPi0;

        for(int i=0; i<nt; i++) {
            if(clusterEnergy[i] < eMin) continue;
            if(vetoEnergy[i] > vetoMax) continue;

            for(int j=i+1; j<nt; j++) {
                if(clusterEnergy[j] < eMin) continue;
                if(vetoEnergy[j] > vetoMax) continue;

                double th1 = theta[i] * deg;
                double th2 = theta[j] * deg;
                double ph1 = phi[i] * deg;
                double ph2 = phi[j] * deg;

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

                TLorentzVector pi0 = g1 + g2;
                double mgg = pi0.M();
                if(mgg <= 0.0) continue;

                double dist = std::fabs(mgg - mpi0);
                if(dist < bestDist) {
                    bestDist = dist;
                    bestM = mgg;
                    bestI = i;
                    bestJ = j;
                    bestPi0 = pi0;
                }
            }
        }

        if(bestI < 0 || bestJ < 0) continue;
        if(bestM < mggMin || bestM > mggMax) continue;

        int ntag = std::min(nTagged, maxTagged);
        for(int it=0; it<ntag; it++) {
            int ch = taggedChannel[it];
            if(ch < 0 || ch >= nCh) continue;

            double egamma = taggerPhotonEnergy[ch];
            if(egamma <= 0.0 || egamma > beamE + 5.0) continue;

            TLorentzVector beam(0.0, 0.0, egamma, egamma);
            TLorentzVector recoil = beam + target - bestPi0;
            double mm = recoil.M() - targetMass;
            TLorentzVector total = beam + target;
            TLorentzVector pi0CM = bestPi0;
            pi0CM.Boost(-total.BoostVector());
            double wcm = total.M();
            double ePiCohCM = (wcm*wcm + mpi0*mpi0 - targetMass*targetMass) / (2.0*wcm);
            double deltaEpi = pi0CM.E() - ePiCohCM;

            double tt = taggedTime[it];
            double weight = 0.0;
            if(tt > promptMin && tt < promptMax) {
                weight = 1.0;
            } else if(tt > randomMin && tt < randomMax) {
                weight = -randomWeight;
            } else {
                continue;
            }

            if(mm < mmMin || mm > mmMax) continue;

            double residual = (mm - (bandP0 + bandP1 * deltaEpi)) / std::sqrt(1.0 + bandP1 * bandP1);
            hBandResidualAfter->Fill(residual, weight);
            for(int ib=0; ib<nBandScan; ib++) {
                if(std::fabs(residual) > bandHalfScan[ib]) continue;
                bandYield[ib] += weight;
                bandErr2[ib] += weight * weight;
                if(weight > 0.0) bandPrompt[ib] += 1.0;
                else bandRandom[ib] += 1.0;
            }
        }
    }

    std::ofstream bandOut("background_diagonal_band_scan.txt");
    bandOut << "# input: " << fname << "\n";
    bandOut << "# all entries below include mgg cut [" << mggMin << "," << mggMax << "] MeV\n";
    bandOut << "# and baseline mm cut [" << mmMin << "," << mmMax << "] MeV\n";
    bandOut << "# fitted band: MM = p0 + p1*DeltaEpi, p0 = " << bandP0
            << ", p1 = " << bandP1 << "\n";
    bandOut << "# residual = (MM - (p0 + p1*DeltaEpi))/sqrt(1+p1^2)\n";
    bandOut << "# residual_max prompt random yield dY yield_over_dY\n";
    for(int ib=0; ib<nBandScan; ib++) {
        double y = bandYield[ib];
        double dy = std::sqrt(bandErr2[ib]);
        double signif = (dy > 0.0) ? y / dy : 0.0;
        bandOut << std::setw(10) << bandHalfScan[ib] << " "
                << std::setw(12) << bandPrompt[ib] << " "
                << std::setw(12) << bandRandom[ib] << " "
                << std::setw(14) << y << " "
                << std::setw(14) << dy << " "
                << std::setw(14) << signif << "\n";
    }
    bandOut.close();

    std::ofstream cutOut("background_cut_stability.txt");
    cutOut << "# input: " << fname << "\n";
    cutOut << "# all entries below include the mgg cut [" << mggMin << "," << mggMax << "] MeV\n";
    cutOut << "# yield = prompt - randomWeight*random, dY = sqrt(prompt + randomWeight^2*random)\n";
    cutOut << "# randomWeight = " << randomWeight << "\n";
    cutOut << "# mm_low mm_high abs_deltaEpi_max prompt random yield dY yield_over_dY\n";
    for(int im=0; im<nMmScan; im++) {
        for(int id=0; id<nDeScan; id++) {
            double y = cutYield[im][id];
            double dy = std::sqrt(cutErr2[im][id]);
            double signif = (dy > 0.0) ? y / dy : 0.0;
            cutOut << std::setw(10) << mmLowScan[im] << " "
                   << std::setw(10) << mmHighScan[im] << " "
                   << std::setw(10) << deHalfScan[id] << " "
                   << std::setw(12) << cutPrompt[im][id] << " "
                   << std::setw(12) << cutRandom[im][id] << " "
                   << std::setw(14) << y << " "
                   << std::setw(14) << dy << " "
                   << std::setw(14) << signif << "\n";
        }
        cutOut << "\n";
    }
    cutOut.close();

    std::ofstream out("sigma_from_acqu.txt");
    out << "# input: " << fname << "\n";
    out << "# processed events = " << nentries << "\n";
    out << "# fpd scaler map = " << fpdFile << "\n";
    out << "# haveFpdMap = " << haveFpdMap << "\n";
    out << "# skipFirstScalerEntry = " << skipFirstScalerEntry << "\n";
    out << "# prompt = [" << promptMin << "," << promptMax << "]\n";
    out << "# random = [" << randomMin << "," << randomMax << "]\n";
    out << "# randomWeight = " << randomWeight << "\n";
    out << "# mgg cut = [" << mggMin << "," << mggMax << "] MeV\n";
    out << "# mm cut = [" << mmMin << "," << mmMax << "] MeV\n";
    out << "# eps_tag = " << eps_tag << "\n";
    out << "# tagging efficiency file = " << taggEffFile << "\n";
    out << "# have tagging efficiency file = " << haveTaggEffFile << "\n";
    out << "# min eps_tag = " << minEpsTag << "\n";
    out << "# eps_det = " << eps_det << "\n";
    out << "# eps_det file = " << epsDetPath << "\n";
    out << "# have eps_det file = " << haveEpsDetFile << "\n";
    out << "# min eps_det = " << minEpsDet << "\n";
    out << "# thickness = " << thickness << " ub^-1\n";
    out << "# nTagger = " << nTagger << "\n";
    out << "# minScalerSum = " << minScalerSum << "\n";
    out << "# channel Egamma scalerIndex eps_tag deps_tag eps_det deps_det Y dY Ne sigma_ub dsigma_ub prompt random status\n";

    TGraphErrors* gr = new TGraphErrors();
    gr->SetName("sigma_from_acqu_vs_channel");
    gr->SetTitle("^{4}He(#gamma,#pi^{0})^{4}He from Acqu;Tagger channel;#sigma (#mub)");
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0.8);
    gr->SetMarkerColor(kBlue+1);
    gr->SetLineColor(kBlue+1);

    TGraphErrors* grE = new TGraphErrors();
    grE->SetName("sigma_from_acqu_vs_egamma");
    grE->SetTitle("^{4}He(#gamma,#pi^{0})^{4}He from Acqu;E_{#gamma} (MeV);#sigma (#mub)");
    grE->SetMarkerStyle(20);
    grE->SetMarkerSize(0.8);
    grE->SetMarkerColor(kBlue+1);
    grE->SetLineColor(kBlue+1);

    TGraphErrors* grRef = load_sigma_reference(referencePath.Data(),
                                               "paper_sigma_vs_egamma",
                                               kBlack,
                                               20,
                                               referenceScale);

    int ip = 0;
    for(int ch=0; ch<nValidCh; ch++) {
        int b1 = hMMp[ch]->GetXaxis()->FindBin(mmMin);
        int b2 = hMMp[ch]->GetXaxis()->FindBin(mmMax);

        double ep = 0.0;
        double er = 0.0;
        double p = hMMp[ch]->IntegralAndError(b1, b2, ep);
        double r = hMMr[ch]->IntegralAndError(b1, b2, er);

        double Y = p - randomWeight * r;
        double dY = std::sqrt(ep*ep + randomWeight*randomWeight*er*er);

        const char* status = "ok";
        if(Ne[ch] <= 0.0) status = "bad_scaler_zero";
        else if(Ne[ch] < minScalerSum) status = "bad_scaler_low";

        double epsDetUse = eps_det;
        double dEpsDetUse = 0.0;
        if(haveEpsDetFile) {
            epsDetUse = epsDetByCh[ch];
            dEpsDetUse = dEpsDetByCh[ch];
            if(epsDetUse <= 0.0) status = "bad_eps_det_zero";
            else if(epsDetUse < minEpsDet) status = "bad_eps_det_low";
        }

        double epsTagUse = eps_tag;
        double dEpsTagUse = 0.0;
        if(haveTaggEffFile) {
            epsTagUse = epsTagByCh[ch];
            dEpsTagUse = dEpsTagByCh[ch];
            if(epsTagUse <= 0.0) status = "bad_eps_tag_zero";
            else if(epsTagUse < minEpsTag) status = "bad_eps_tag_low";
        }

        double sigma = 0.0;
        double dsigma = 0.0;
        if(Ne[ch] > 0.0 && epsDetUse > 0.0 && epsTagUse > 0.0) {
            sigma = Y / (Ne[ch] * epsTagUse * epsDetUse * thickness);
        }
        double dNe = std::sqrt(dNe2[ch]);
        double rel2 = 0.0;
        if(Y > 0.0 && dY > 0.0) rel2 += (dY/Y) * (dY/Y);
        if(Ne[ch] > 0.0 && dNe > 0.0) rel2 += (dNe/Ne[ch]) * (dNe/Ne[ch]);
        if(epsTagUse > 0.0 && dEpsTagUse > 0.0) rel2 += (dEpsTagUse/epsTagUse) * (dEpsTagUse/epsTagUse);
        if(epsDetUse > 0.0 && dEpsDetUse > 0.0) rel2 += (dEpsDetUse/epsDetUse) * (dEpsDetUse/epsDetUse);
        dsigma = std::fabs(sigma) * std::sqrt(rel2);

        out << ch << " "
            << taggerPhotonEnergy[ch] << " "
            << scalerIndex[ch] << " "
            << epsTagUse << " "
            << dEpsTagUse << " "
            << epsDetUse << " "
            << dEpsDetUse << " "
            << Y << " "
            << dY << " "
            << Ne[ch] << " "
            << sigma << " "
            << dsigma << " "
            << p << " "
            << r << " "
            << status << "\n";

        if(Y > 0.0 && Ne[ch] >= minScalerSum && epsDetUse >= minEpsDet && epsTagUse >= minEpsTag) {
            gr->SetPoint(ip, ch, sigma);
            gr->SetPointError(ip, 0.0, dsigma);
            grE->SetPoint(ip, taggerPhotonEnergy[ch], sigma);
            grE->SetPointError(ip, 0.0, dsigma);
            ip++;
        }
    }
    out.close();

    TCanvas* c1 = new TCanvas("c_mgg_from_acqu", "mgg from Acqu", 800, 600);
    hMggBest->Draw();
    c1->SaveAs("mgg_best_from_acqu.pdf");

    TCanvas* cMM = new TCanvas("c_mm_before_after_from_acqu", "missing mass before/after", 1200, 600);
    cMM->Divide(2, 1);
    cMM->cd(1);
    hMMBefore->Draw("HIST");
    cMM->cd(2);
    hMMAfter->Draw("HIST");
    cMM->SaveAs("missing_mass_before_after_from_acqu.pdf");

    TCanvas* c2D = new TCanvas("c_mgg_vs_mm_before_after_from_acqu", "mgg vs missing mass before/after", 1200, 600);
    c2D->Divide(2, 1);
    c2D->cd(1);
    hMggMMBefore->Draw("COLZ");
    c2D->cd(2);
    hMggMMAfter->Draw("COLZ");
    c2D->SaveAs("mgg_vs_missing_mass_before_after_from_acqu.pdf");

    TCanvas* cMME = new TCanvas("c_mm_vs_egamma_before_after_from_acqu", "missing mass vs Egamma before/after", 1200, 600);
    cMME->Divide(2, 1);
    cMME->cd(1);
    hMMEgammaBefore->Draw("COLZ");
    cMME->cd(2);
    hMMEgammaAfter->Draw("COLZ");
    cMME->SaveAs("missing_mass_vs_egamma_before_after_from_acqu.pdf");

    TCanvas* cETheta = new TCanvas("c_egamma_vs_theta_before_after_from_acqu", "Egamma vs theta before/after", 1200, 600);
    cETheta->Divide(2, 1);
    cETheta->cd(1);
    hEgammaThetaBefore->Draw("COLZ");
    cETheta->cd(2);
    hEgammaThetaAfter->Draw("COLZ");
    cETheta->SaveAs("egamma_vs_theta_before_after_from_acqu.pdf");

    TCanvas* cDelta = new TCanvas("c_delta_epi_before_after_from_acqu", "Delta Epi before/after", 1200, 600);
    cDelta->Divide(2, 1);
    cDelta->cd(1);
    hDeltaEpiBefore->Draw("HIST");
    cDelta->cd(2);
    hDeltaEpiAfter->Draw("HIST");
    cDelta->SaveAs("delta_epi_before_after_from_acqu.pdf");

    TCanvas* cMMDelta = new TCanvas("c_mm_vs_delta_epi_before_after_from_acqu", "MM vs Delta Epi before/after", 1200, 600);
    cMMDelta->Divide(2, 1);
    cMMDelta->cd(1);
    hMMDeltaBefore->Draw("COLZ");
    cMMDelta->cd(2);
    hMMDeltaAfter->Draw("COLZ");
    cMMDelta->SaveAs("missing_mass_vs_delta_epi_before_after_from_acqu.pdf");

    TCanvas* cBandFit = new TCanvas("c_mm_delta_epi_band_fit_from_acqu", "MM vs Delta Epi band fit", 900, 700);
    hMMDeltaAfter->Draw("COLZ");
    fMMDeltaBand->SetLineColor(kRed);
    fMMDeltaBand->SetLineWidth(3);
    fMMDeltaBand->Draw("SAME");
    cBandFit->SaveAs("missing_mass_vs_delta_epi_band_fit_from_acqu.pdf");

    TCanvas* cBandResidual = new TCanvas("c_diagonal_band_residual_from_acqu", "diagonal band residual", 900, 700);
    hBandResidualAfter->Draw("HIST");
    cBandResidual->SaveAs("diagonal_band_residual_from_acqu.pdf");

    TCanvas* cMggDelta = new TCanvas("c_mgg_vs_delta_epi_from_acqu", "mgg vs Delta Epi", 900, 700);
    hMggDeltaBefore->Draw("COLZ");
    cMggDelta->SaveAs("mgg_vs_delta_epi_from_acqu.pdf");

    TCanvas* cEgammaDelta = new TCanvas("c_egamma_vs_delta_epi_from_acqu", "Egamma vs Delta Epi", 900, 700);
    hEgammaDeltaAfter->Draw("COLZ");
    cEgammaDelta->SaveAs("egamma_vs_delta_epi_after_mgg_from_acqu.pdf");

    TCanvas* cThetaDelta = new TCanvas("c_theta_cm_vs_delta_epi_from_acqu", "theta CM vs Delta Epi", 900, 700);
    hThetaCmDeltaAfter->Draw("COLZ");
    cThetaDelta->SaveAs("theta_cm_vs_delta_epi_after_mgg_from_acqu.pdf");

    TCanvas* cThetaLabCm = new TCanvas("c_theta_lab_vs_theta_cm_from_acqu", "theta lab vs theta CM", 900, 700);
    hThetaLabCmAfter->Draw("COLZ");
    cThetaLabCm->SaveAs("theta_lab_vs_theta_cm_after_mgg_from_acqu.pdf");

    TCanvas* c2 = new TCanvas("c_sigma_from_acqu", "sigma from Acqu", 1100, 700);
    c2->SetGrid();
    gr->Draw("AP");
    c2->SaveAs("sigma_from_acqu.pdf");

    TCanvas* c3 = new TCanvas("c_sigma_from_acqu_egamma", "sigma from Acqu vs Egamma", 1100, 700);
    c3->SetGrid();
    grE->Draw("AP");
    if(grRef) {
        grRef->Draw("P SAME");
        TLegend* leg = new TLegend(0.58, 0.72, 0.88, 0.88);
        leg->AddEntry(grE, "Acqu extraction", "p");
        leg->AddEntry(grRef, "TAPS paper", "p");
        leg->Draw();
    }
    c3->SaveAs("sigma_from_acqu_vs_egamma.pdf");

    if(grRef) {
        c3->SaveAs("sigma_from_acqu_vs_egamma_with_reference.pdf");
    }

    TFile* fout = new TFile("background_diagnostics_from_acqu.root", "RECREATE");
    hMggBest->Write();
    hMMBefore->Write();
    hMMAfter->Write();
    hMggMMBefore->Write();
    hMggMMAfter->Write();
    hMMEgammaBefore->Write();
    hMMEgammaAfter->Write();
    hEgammaThetaBefore->Write();
    hEgammaThetaAfter->Write();
    hDeltaEpiBefore->Write();
    hDeltaEpiAfter->Write();
    hMMDeltaBefore->Write();
    hMMDeltaAfter->Write();
    hMMDeltaPromptAfter->Write();
    hMggDeltaBefore->Write();
    hEgammaDeltaAfter->Write();
    hThetaCmDeltaAfter->Write();
    hThetaLabCmAfter->Write();
    hBandResidualAfter->Write();
    pMMDeltaBand->Write();
    fMMDeltaBand->Write();
    for(int ch=0; ch<nCh; ch++) {
        hMMp[ch]->Write();
        hMMr[ch]->Write();
    }
    gr->Write();
    grE->Write();
    if(grRef) grRef->Write();
    c1->Write();
    cMM->Write();
    c2D->Write();
    cMME->Write();
    cETheta->Write();
    cDelta->Write();
    cMMDelta->Write();
    cBandFit->Write();
    cBandResidual->Write();
    cMggDelta->Write();
    cEgammaDelta->Write();
    cThetaDelta->Write();
    cThetaLabCm->Write();
    c2->Write();
    c3->Write();
    fout->Close();

    std::cout << "Input file: " << fname << std::endl;
    std::cout << "Processed events = " << nentries << std::endl;
    std::cout << "Events with gamma-gamma pair = " << nPairEvents << std::endl;
    std::cout << "Events passing pi0 mass cut = " << nPi0Events << std::endl;
    std::cout << "Prompt tags before mgg cut = " << nPromptTagsBefore << std::endl;
    std::cout << "Random tags before mgg cut = " << nRandomTagsBefore << std::endl;
    std::cout << "Prompt tags used = " << nPromptTags << std::endl;
    std::cout << "Random tags used = " << nRandomTags << std::endl;
    std::cout << "Random weight = " << randomWeight << std::endl;
    std::cout << "Saved:" << std::endl;
    std::cout << "  background_cut_stability.txt" << std::endl;
    std::cout << "  background_diagonal_band_scan.txt" << std::endl;
    std::cout << "  delta_epi_before_after_from_acqu.pdf" << std::endl;
    std::cout << "  missing_mass_vs_delta_epi_before_after_from_acqu.pdf" << std::endl;
    std::cout << "  missing_mass_vs_delta_epi_band_fit_from_acqu.pdf" << std::endl;
    std::cout << "  diagonal_band_residual_from_acqu.pdf" << std::endl;
    std::cout << "  mgg_vs_delta_epi_from_acqu.pdf" << std::endl;
    std::cout << "  egamma_vs_delta_epi_after_mgg_from_acqu.pdf" << std::endl;
    std::cout << "  theta_cm_vs_delta_epi_after_mgg_from_acqu.pdf" << std::endl;
    std::cout << "  theta_lab_vs_theta_cm_after_mgg_from_acqu.pdf" << std::endl;
    std::cout << "  sigma_from_acqu.txt" << std::endl;
    std::cout << "  sigma_from_acqu.pdf" << std::endl;
    std::cout << "  sigma_from_acqu_vs_egamma.pdf" << std::endl;
    if(grRef) std::cout << "  sigma_from_acqu_vs_egamma_with_reference.pdf" << std::endl;
    std::cout << "  mgg_best_from_acqu.pdf" << std::endl;
    std::cout << "  missing_mass_before_after_from_acqu.pdf" << std::endl;
    std::cout << "  mgg_vs_missing_mass_before_after_from_acqu.pdf" << std::endl;
    std::cout << "  missing_mass_vs_egamma_before_after_from_acqu.pdf" << std::endl;
    std::cout << "  egamma_vs_theta_before_after_from_acqu.pdf" << std::endl;
    std::cout << "  background_diagnostics_from_acqu.root" << std::endl;

    f->Close();
}
