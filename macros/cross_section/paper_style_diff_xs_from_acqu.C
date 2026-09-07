/*
 * paper_style_diff_xs_from_acqu.C
 *
 * PURPOSE
 *   Extract d(sigma)/d(Omega) for 4He(gamma,pi0)4He using the 17 laboratory
 *   photon-energy intervals used for the paper-style comparison and angular
 *   bins in the pi0 centre-of-mass angle.
 *
 * DATA SELECTION
 *   - Reads tracks, tagger, scalers and setupParameters.
 *   - Cluster threshold: 20 MeV; vetoEnergy <= 1 MeV.
 *   - Chooses the gamma-gamma pair closest to the pi0 mass.
 *   - Default m(gamma gamma) window: 110-155 MeV.
 *   - Performs prompt-random subtraction using the configurable time windows.
 *   - Yield in every (Egamma, theta_cm) bin is obtained by integrating the
 *     missing-mass histogram over the configurable MM window (default
 *     -10 to 40 MeV).
 *
 * IMPORTANT: OPENING-ANGLE QUANTITY
 *   The macro calculates Phi_min and produces several opening-angle/Phi_min
 *   diagnostic plots, but the event loop DOES NOT reject events using a
 *   DeltaPhi/Phi_min cut. Therefore the differential cross section produced
 *   here uses the m(gamma gamma) and missing-mass selections, not the Fig. 2
 *   DeltaPhi coherent-selection cut.
 *
 * NORMALIZATION
 *   The photon-flux denominator is accumulated from tagger electron scalers,
 *   tagging efficiency and target thickness over all tagger channels inside
 *   each paper energy bin.
 *
 * DETECTION EFFICIENCY MODES
 *   1. Constant eps_det fallback.
 *   2. Channel-dependent efficiency from epsDetFile (historical
 *      eps_det_reco_fast_from_acqu_geant.txt format).
 *   3. If usePaperBinEps=true (or epsDetPaperFile is supplied), a detection
 *      efficiency specific to each (Egamma, theta_cm) bin is used. The files
 *      produced by eps_det_paper_bins_from_acqu_geant.C and
 *      eps_det_paper_bins_from_goat.C follow the expected text format. A
 *      requested paper-bin efficiency file must be readable; otherwise the
 *      macro aborts rather than silently changing normalization mode.
 *
 * BINNING
 *   Photon-energy bins (MeV):
 *     201-210, 211-222, 223-234, 235-246, 247-258, 259-270,
 *     271-282, 283-294, 295-308, 309-318, 319-330, 331-342,
 *     343-355, 356-366, 367-378, 379-390, 391-401.
 *   Default theta_cm range: 5-150 deg in 5-deg bins.
 *
 * CROSS SECTION
 *   d(sigma)/d(Omega) = Y /
 *       (photon-flux normalization * eps_det(E,theta) * DeltaOmega)
 *
 * MAIN OUTPUTS
 *   paper_style_dsigma_dOmega_cm_from_acqu.txt
 *   paper_style_dsigma_dOmega_cm_from_acqu.root
 *   paper_style_dsigma_dOmega_cm_by_energy.pdf
 *   paper_style_dsigma_dOmega_cm_map.pdf
 *   plus diagnostic and reference-comparison PDFs written by the macro.
 *
 * A convenience function, paper_style_diff_xs_from_acqu_goat_eps(), is also
 * defined at the bottom of this same original file and enables the GoAT
 * paper-bin efficiency file.
 *
 * NOTE FOR THE REPOSITORY
 *   Physics selections and the public function interface are preserved.
 *   Repository preparation removed machine-specific default paths and makes requested
 *   production normalization/efficiency inputs fail closed when unreadable.
 */
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TTree.h"
#include "TVirtualPad.h"

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

bool load_eps_det_paper_bins(const char* epsFile,
                             const std::vector<double>& eBinLow,
                             const std::vector<double>& eBinHigh,
                             double thetaCmMinDeg,
                             double thetaCmBinWidthDeg,
                             int nThetaBins,
                             std::vector< std::vector<double> >& epsDet,
                             std::vector< std::vector<double> >& dEpsDet)
{
    int nEBins = int(eBinLow.size());
    epsDet.assign(nEBins, std::vector<double>(nThetaBins, -1.0));
    dEpsDet.assign(nEBins, std::vector<double>(nThetaBins, 0.0));

    if(!epsFile || std::string(epsFile).empty()) return false;

    std::ifstream in(epsFile);
    if(!in.is_open()) {
        std::cout << "Warning: could not open paper-bin eps_det file " << epsFile << std::endl;
        return false;
    }

    std::string line;
    while(std::getline(in, line)) {
        if(line.empty() || line[0] == '#') continue;

        std::istringstream ss(line);
        double eLow = 0.0;
        double eHigh = 0.0;
        double eCenter = 0.0;
        double thLow = 0.0;
        double thHigh = 0.0;
        double thCenter = 0.0;
        double deltaOmega = 0.0;
        double nGenE = 0.0;
        double nGen = 0.0;
        double nReco = 0.0;
        double eps = 0.0;
        double deps = 0.0;

        if(!(ss >> eLow >> eHigh >> eCenter
                >> thLow >> thHigh >> thCenter
                >> deltaOmega >> nGenE >> nGen >> nReco >> eps >> deps)) {
            continue;
        }

        int ie = -1;
        for(int i=0; i<nEBins; i++) {
            if(std::fabs(eLow - eBinLow[i]) < 0.1 && std::fabs(eHigh - eBinHigh[i]) < 0.1) {
                ie = i;
                break;
            }
        }
        if(ie < 0) continue;

        int ith = int(std::floor((thCenter - thetaCmMinDeg) / thetaCmBinWidthDeg));
        if(ith < 0 || ith >= nThetaBins) continue;

        epsDet[ie][ith] = eps;
        dEpsDet[ie][ith] = deps;
    }

    int nLoaded = 0;
    for(int ie=0; ie<nEBins; ie++) {
        for(int ith=0; ith<nThetaBins; ith++) {
            if(epsDet[ie][ith] > 0.0) nLoaded++;
        }
    }

    return nLoaded > 0;
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

double estimate_energy_bin_width(const double energy[], int ch, int nValidCh)
{
    if(ch < 0 || ch >= nValidCh || energy[ch] <= 0.0) return 0.0;

    int prev = -1;
    int next = -1;

    for(int i=ch-1; i>=0; i--) {
        if(energy[i] > 0.0) {
            prev = i;
            break;
        }
    }

    for(int i=ch+1; i<nValidCh; i++) {
        if(energy[i] > 0.0) {
            next = i;
            break;
        }
    }

    if(prev >= 0 && next >= 0) return 0.5 * std::fabs(energy[next] - energy[prev]);
    if(prev >= 0) return std::fabs(energy[ch] - energy[prev]);
    if(next >= 0) return std::fabs(energy[next] - energy[ch]);

    return 0.0;
}

double solid_angle_bin(double thetaMinDeg, double thetaMaxDeg)
{
    double th1 = thetaMinDeg * TMath::DegToRad();
    double th2 = thetaMaxDeg * TMath::DegToRad();
    return 2.0 * TMath::Pi() * (std::cos(th1) - std::cos(th2));
}

int find_bin(double x, const std::vector<double>& low, const std::vector<double>& high)
{
    for(size_t i=0; i<low.size(); i++) {
        if(x >= low[i] && x <= high[i]) return int(i);
    }
    return -1;
}

TGraphErrors* make_table1_graph(int ie, double eLow, double eHigh)
{
    static const char* table1Data =
        "201 210 205.5 7.0 0.1 0.1\n"
        "211 222 216.5 7.0 0.3 0.1\n"
        "223 234 228.5 7.0 0.4 0.1\n"
        "235 246 240.5 7.0 0.4 0.2\n"
        "247 258 252.5 7.0 1.8 0.3\n"
        "259 270 264.5 7.0 1.3 0.3\n"
        "201 210 205.5 12.0 2.1 0.1\n"
        "211 222 216.5 12.0 2.9 0.1\n"
        "223 234 228.5 12.0 3.6 0.2\n"
        "235 246 240.5 12.0 5.4 0.3\n"
        "247 258 252.5 12.0 6.2 0.3\n"
        "259 270 264.5 12.0 7.5 0.3\n"
        "201 210 205.5 17.0 4.2 0.2\n"
        "211 222 216.5 17.0 5.5 0.2\n"
        "223 234 228.5 17.0 6.9 0.2\n"
        "235 246 240.5 17.0 9.2 0.3\n"
        "247 258 252.5 17.0 11.2 0.3\n"
        "259 270 264.5 17.0 13.2 0.3\n"
        "201 210 205.5 22.0 6.9 0.2\n"
        "211 222 216.5 22.0 8.3 0.2\n"
        "223 234 228.5 22.0 10.6 0.2\n"
        "235 246 240.5 22.0 14.3 0.3\n"
        "247 258 252.5 22.0 16.8 0.3\n"
        "259 270 264.5 22.0 20.3 0.3\n"
        "201 210 205.5 27.0 9.2 0.3\n"
        "211 222 216.5 27.0 11.5 0.3\n"
        "223 234 228.5 27.0 14.0 0.3\n"
        "235 246 240.5 27.0 18.6 0.4\n"
        "247 258 252.5 27.0 22.4 0.4\n"
        "259 270 264.5 27.0 26.6 0.4\n"
        "201 210 205.5 32.0 11.4 0.3\n"
        "211 222 216.5 32.0 13.3 0.3\n"
        "223 234 228.5 32.0 15.8 0.4\n"
        "235 246 240.5 32.0 22.5 0.4\n"
        "247 258 252.5 32.0 28.7 0.4\n"
        "259 270 264.5 32.0 33.2 0.4\n"
        "201 210 205.5 37.0 12.7 0.4\n"
        "211 222 216.5 37.0 15.9 0.4\n"
        "223 234 228.5 37.0 18.9 0.4\n"
        "235 246 240.5 37.0 26.4 0.5\n"
        "247 258 252.5 37.0 33.2 0.5\n"
        "259 270 264.5 37.0 39.7 0.6\n"
        "201 210 205.5 42.0 15.0 0.4\n"
        "211 222 216.5 42.0 17.7 0.4\n"
        "223 234 228.5 42.0 23.6 0.5\n"
        "235 246 240.5 42.0 29.9 0.6\n"
        "247 258 252.5 42.0 37.6 0.7\n"
        "259 270 264.5 42.0 43.2 0.7\n"
        "201 210 205.5 47.0 15.9 0.4\n"
        "223 234 228.5 62.0 31.5 0.7\n"
        "235 246 240.5 62.0 41.4 0.8\n"
        "247 258 252.5 62.0 48.0 0.8\n"
        "259 270 264.5 62.0 54.4 0.8\n"
        "211 222 216.5 67.0 24.7 0.7\n"
        "223 234 228.5 67.0 34.8 0.6\n"
        "235 246 240.5 67.0 40.7 0.5\n"
        "247 258 252.5 67.0 46.2 0.5\n"
        "259 270 264.5 67.0 49.2 0.5\n"
        "201 210 205.5 72.0 19.5 0.6\n"
        "211 222 216.5 72.0 25.7 0.6\n"
        "223 234 228.5 72.0 32.5 0.5\n"
        "235 246 240.5 72.0 39.0 0.5\n"
        "247 258 252.5 72.0 42.7 0.4\n"
        "259 270 264.5 72.0 45.3 0.4\n"
        "201 210 205.5 77.0 20.4 0.5\n"
        "211 222 216.5 77.0 24.4 0.4\n"
        "223 234 228.5 77.0 30.2 0.4\n"
        "235 246 240.5 77.0 36.5 0.4\n"
        "247 258 252.5 77.0 39.5 0.4\n"
        "259 270 264.5 77.0 41.0 0.4\n"
        "201 210 205.5 82.0 19.2 0.4\n"
        "211 222 216.5 82.0 23.8 0.4\n"
        "223 234 228.5 82.0 27.5 0.4\n"
        "235 246 240.5 82.0 32.7 0.4\n"
        "247 258 252.5 82.0 35.2 0.4\n"
        "259 270 264.5 82.0 35.8 0.3\n"
        "201 210 205.5 87.0 17.6 0.3\n"
        "211 222 216.5 87.0 20.1 0.3\n"
        "223 234 228.5 87.0 24.9 0.3\n"
        "235 246 240.5 87.0 29.1 0.4\n"
        "247 258 252.5 87.0 30.8 0.4\n"
        "259 270 264.5 87.0 30.3 0.3\n"
        "201 210 205.5 92.0 17.0 0.3\n"
        "211 222 216.5 92.0 19.1 0.3\n"
        "223 234 228.5 92.0 22.0 0.3\n"
        "235 246 240.5 92.0 26.0 0.4\n"
        "247 258 252.5 92.0 27.4 0.4\n"
        "259 270 264.5 92.0 26.6 0.4\n"
        "201 210 205.5 97.0 14.9 0.3\n"
        "211 222 216.5 97.0 16.7 0.3\n"
        "223 234 228.5 97.0 19.3 0.3\n"
        "235 246 240.5 97.0 21.9 0.4\n"
        "247 258 252.5 97.0 22.5 0.4\n"
        "259 270 264.5 97.0 21.7 0.5\n"
        "201 210 205.5 102.0 13.3 0.2\n"
        "211 222 216.5 102.0 14.1 0.2\n"
        "223 234 228.5 102.0 15.6 0.3\n"
        "235 246 240.5 102.0 17.1 0.3\n"
        "247 258 252.5 102.0 16.4 0.4\n"
        "259 270 264.5 102.0 15.4 0.5\n"
        "201 210 205.5 107.0 11.9 0.2\n"
        "211 222 216.5 107.0 12.7 0.2\n"
        "223 234 228.5 107.0 13.0 0.3\n"
        "235 246 240.5 107.0 14.3 0.3\n"
        "247 258 252.5 107.0 12.6 0.4\n"
        "259 270 264.5 107.0 12.0 0.4\n"
        "201 210 205.5 112.0 10.3 0.2\n"
        "211 222 216.5 112.0 11.1 0.2\n"
        "223 234 228.5 112.0 11.8 0.3\n"
        "235 246 240.5 112.0 11.9 0.3\n"
        "247 258 252.5 112.0 10.2 0.3\n"
        "259 270 264.5 112.0 8.9 0.3\n"
        "201 210 205.5 117.0 8.7 0.2\n"
        "211 222 216.5 117.0 8.8 0.2\n"
        "223 234 228.5 117.0 9.2 0.2\n"
        "235 246 240.5 117.0 9.6 0.2\n"
        "247 258 252.5 117.0 8.2 0.2\n"
        "259 270 264.5 117.0 6.3 0.2\n"
        "201 210 205.5 122.0 7.0 0.2\n"
        "211 222 216.5 122.0 7.8 0.2\n"
        "223 234 228.5 122.0 7.4 0.2\n"
        "235 246 240.5 122.0 7.0 0.2\n"
        "247 258 252.5 122.0 6.6 0.2\n"
        "259 270 264.5 122.0 4.8 0.2\n"
        "201 210 205.5 127.0 6.1 0.3\n"
        "211 222 216.5 127.0 6.0 0.2\n"
        "223 234 228.5 127.0 6.1 0.3\n"
        "235 246 240.5 127.0 6.4 0.2\n"
        "247 258 252.5 127.0 5.0 0.2\n"
        "259 270 264.5 127.0 3.6 0.2\n"
        "201 210 205.5 132.0 5.1 0.3\n"
        "211 222 216.5 132.0 5.1 0.2\n"
        "223 234 228.5 132.0 5.5 0.3\n"
        "235 246 240.5 132.0 4.9 0.3\n"
        "247 258 252.5 132.0 3.7 0.2\n"
        "259 270 264.5 132.0 2.5 0.1\n"
        "201 210 205.5 137.0 3.3 0.2\n"
        "211 222 216.5 137.0 3.9 0.2\n"
        "223 234 228.5 137.0 3.1 0.3\n"
        "235 246 240.5 137.0 3.2 0.3\n"
        "247 258 252.5 137.0 2.9 0.2\n"
        "259 270 264.5 137.0 1.9 0.2\n"
        "201 210 205.5 142.0 3.1 0.2\n"
        "211 222 216.5 142.0 3.5 0.2\n"
        "223 234 228.5 142.0 2.8 0.4\n"
        "235 246 240.5 142.0 0.8 0.3\n"
        "247 258 252.5 142.0 1.0 0.2\n"
        "259 270 264.5 142.0 0.5 0.2\n"
        "201 210 205.5 147.0 2.8 0.2\n"
        "211 222 216.5 147.0 2.3 0.2\n"
        "223 234 228.5 147.0 1.5 0.3\n"
        "271 282 276.5 7.0 2.2 0.3\n"
        "283 294 288.5 7.0 3.4 0.3\n"
        "295 308 301.5 7.0 2.9 0.2\n"
        "309 318 313.5 7.0 4.5 0.3\n"
        "319 330 324.5 7.0 4.2 0.3\n"
        "331 342 336.5 7.0 2.5 0.3\n"
        "271 282 276.5 12.0 8.8 0.3\n"
        "283 294 288.5 12.0 10.3 0.3\n"
        "295 308 301.5 12.0 10.3 0.3\n"
        "309 318 313.5 12.0 11.9 0.4\n"
        "319 330 324.5 12.0 11.6 0.3\n"
        "331 342 336.5 12.0 10.9 0.3\n"
        "271 282 276.5 17.0 15.5 0.3\n"
        "283 294 288.5 17.0 16.7 0.3\n"
        "295 308 301.5 17.0 19.0 0.3\n"
        "309 318 313.5 17.0 20.2 0.4\n"
        "319 330 324.5 17.0 20.5 0.3\n"
        "331 342 336.5 17.0 20.5 0.3\n"
        "271 282 276.5 22.0 23.0 0.4\n"
        "283 294 288.5 22.0 24.7 0.4\n"
        "295 308 301.5 22.0 27.3 0.4\n"
        "309 318 313.5 22.0 29.2 0.5\n"
        "319 330 324.5 22.0 29.3 0.4\n"
        "331 342 336.5 22.0 30.7 0.4\n"
        "271 282 276.5 27.0 30.0 0.4\n"
        "283 294 288.5 27.0 32.6 0.4\n"
        "295 308 301.5 27.0 35.6 0.4\n"
        "309 318 313.5 27.0 38.4 0.5\n"
        "319 330 324.5 27.0 37.9 0.4\n"
        "331 342 336.5 27.0 37.4 0.4\n"
        "271 282 276.5 32.0 36.9 0.4\n"
        "283 294 288.5 32.0 40.3 0.4\n"
        "295 308 301.5 32.0 42.7 0.4\n"
        "309 318 313.5 32.0 46.3 0.6\n"
        "319 330 324.5 32.0 44.6 0.4\n"
        "331 342 336.5 32.0 43.8 0.4\n"
        "271 282 276.5 37.0 44.1 0.5\n"
        "283 294 288.5 37.0 45.9 0.5\n"
        "295 308 301.5 37.0 49.3 0.5\n"
        "309 318 313.5 37.0 50.5 0.6\n"
        "319 330 324.5 37.0 50.0 0.5\n"
        "331 342 336.5 37.0 47.8 0.5\n"
        "271 282 276.5 42.0 47.7 0.7\n"
        "283 294 288.5 42.0 50.6 0.6\n"
        "295 308 301.5 42.0 53.0 0.6\n"
        "309 318 313.5 42.0 54.4 0.8\n"
        "319 330 324.5 42.0 51.9 0.6\n"
        "331 342 336.5 42.0 50.0 0.5\n"
        "271 282 276.5 47.0 52.2 1.2\n"
        "283 294 288.5 47.0 52.6 1.0\n"
        "295 308 301.5 47.0 54.2 0.9\n"
        "309 318 313.5 47.0 53.5 1.0\n"
        "319 330 324.5 47.0 51.7 0.7\n"
        "331 342 336.5 47.0 48.5 0.7\n"
        "271 282 276.5 51.0 51.4 3.3\n"
        "283 294 288.5 51.0 49.2 2.4\n"
        "295 308 301.5 51.0 55.4 2.0\n"
        "309 318 313.5 51.0 55.8 1.9\n"
        "319 330 324.5 51.0 49.6 1.2\n"
        "331 342 336.5 51.0 45.9 1.1\n"
        "271 282 276.5 58.0 50.0 1.6\n"
        "283 294 288.5 58.0 51.3 1.4\n"
        "295 308 301.5 58.0 51.2 1.4\n"
        "309 318 313.5 58.0 52.3 1.4\n"
        "319 330 324.5 58.0 48.8 1.0\n"
        "331 342 336.5 58.0 40.7 0.8\n"
        "271 282 276.5 62.0 52.9 0.7\n"
        "283 294 288.5 62.0 49.7 0.6\n"
        "295 308 301.5 62.0 50.7 0.6\n"
        "309 318 313.5 62.0 44.9 0.7\n"
        "319 330 324.5 62.0 41.0 0.5\n"
        "331 342 336.5 62.0 33.4 0.4\n"
        "271 282 276.5 67.0 49.3 0.5\n"
        "283 294 288.5 67.0 46.2 0.4\n"
        "295 308 301.5 67.0 44.2 0.4\n"
        "309 318 313.5 67.0 39.0 0.5\n"
        "319 330 324.5 67.0 33.9 0.3\n"
        "331 342 336.5 67.0 28.3 0.3\n"
        "271 282 276.5 72.0 44.0 0.4\n"
        "283 294 288.5 72.0 41.0 0.4\n"
        "295 308 301.5 72.0 37.8 0.3\n"
        "309 318 313.5 72.0 32.9 0.4\n"
        "319 330 324.5 72.0 27.9 0.3\n"
        "331 342 336.5 72.0 22.0 0.2\n"
        "271 282 276.5 77.0 39.4 0.3\n"
        "283 294 288.5 77.0 34.4 0.3\n"
        "295 308 301.5 77.0 31.0 0.3\n"
        "309 318 313.5 77.0 25.6 0.3\n"
        "319 330 324.5 77.0 20.8 0.2\n"
        "331 342 336.5 77.0 16.4 0.2\n"
        "271 282 276.5 82.0 33.9 0.3\n"
        "283 294 288.5 82.0 29.3 0.3\n"
        "295 308 301.5 82.0 24.7 0.2\n"
        "309 318 313.5 82.0 20.6 0.3\n"
        "319 330 324.5 82.0 15.5 0.2\n"
        "331 342 336.5 82.0 12.1 0.2\n"
        "271 282 276.5 87.0 28.2 0.3\n"
        "283 294 288.5 87.0 23.9 0.3\n"
        "295 308 301.5 87.0 19.0 0.2\n"
        "309 318 313.5 87.0 15.0 0.3\n"
        "319 330 324.5 87.0 10.9 0.2\n"
        "331 342 336.5 87.0 8.8 0.1\n"
        "271 282 276.5 92.0 24.3 0.4\n"
        "283 294 288.5 92.0 19.0 0.3\n"
        "295 308 301.5 92.0 14.4 0.2\n"
        "309 318 313.5 92.0 10.8 0.3\n"
        "319 330 324.5 92.0 7.8 0.2\n"
        "331 342 336.5 92.0 6.0 0.1\n"
        "271 282 276.5 97.0 20.1 0.5\n"
        "283 294 288.5 97.0 17.0 0.4\n"
        "295 308 301.5 97.0 10.8 0.3\n"
        "309 318 313.5 97.0 8.0 0.3\n"
        "319 330 324.5 97.0 5.2 0.2\n"
        "331 342 336.5 97.0 4.3 0.1\n"
        "271 282 276.5 101.0 12.1 0.6\n"
        "283 294 288.5 101.0 10.8 0.8\n"
        "295 308 301.5 101.0 7.1 0.5\n"
        "309 318 313.5 101.0 4.9 0.4\n"
        "319 330 324.5 101.0 3.7 0.2\n"
        "331 342 336.5 101.0 2.1 0.2\n"
        "271 282 276.5 107.0 7.5 0.5\n"
        "283 294 288.5 107.0 4.7 0.5\n"
        "295 308 301.5 107.0 1.6 0.3\n"
        "309 318 313.5 107.0 3.0 0.4\n"
        "319 330 324.5 107.0 1.9 0.2\n"
        "331 342 336.5 107.0 1.4 0.2\n"
        "271 282 276.5 112.0 6.4 0.3\n"
        "283 294 288.5 112.0 4.3 0.2\n"
        "295 308 301.5 112.0 2.6 0.2\n"
        "309 318 313.5 112.0 2.1 0.2\n"
        "319 330 324.5 112.0 1.7 0.1\n"
        "331 342 336.5 112.0 1.6 0.1\n"
        "271 282 276.5 117.0 4.1 0.2\n"
        "283 294 288.5 117.0 2.7 0.1\n"
        "295 308 301.5 117.0 2.0 0.1\n"
        "309 318 313.5 117.0 1.8 0.1\n"
        "319 330 324.5 117.0 1.5 0.1\n"
        "331 342 336.5 117.0 1.3 0.1\n"
        "271 282 276.5 122.0 3.2 0.1\n"
        "283 294 288.5 122.0 2.0 0.1\n"
        "295 308 301.5 122.0 1.6 0.1\n"
        "309 318 313.5 122.0 1.3 0.1\n"
        "319 330 324.5 122.0 1.4 0.1\n"
        "331 342 336.5 122.0 1.3 0.1\n"
        "271 282 276.5 127.0 2.2 0.1\n"
        "283 294 288.5 127.0 1.3 0.1\n"
        "295 308 301.5 127.0 1.3 0.1\n"
        "309 318 313.5 127.0 1.3 0.1\n"
        "319 330 324.5 127.0 1.2 0.1\n"
        "331 342 336.5 127.0 1.2 0.1\n"
        "271 282 276.5 132.0 1.4 0.1\n"
        "283 294 288.5 132.0 0.9 0.1\n"
        "295 308 301.5 132.0 1.0 0.1\n"
        "309 318 313.5 132.0 1.1 0.1\n"
        "319 330 324.5 132.0 1.1 0.1\n"
        "331 342 336.5 132.0 1.0 0.1\n"
        "271 282 276.5 137.0 1.0 0.1\n"
        "283 294 288.5 137.0 0.8 0.1\n"
        "295 308 301.5 137.0 0.9 0.1\n"
        "309 318 313.5 137.0 1.1 0.1\n"
        "319 330 324.5 137.0 1.1 0.1\n"
        "331 342 336.5 137.0 0.9 0.1\n"
        "343 355 349.0 8.0 2.4 0.3\n"
        "356 366 361.0 8.0 2.0 0.4\n"
        "379 390 384.5 8.0 9.9 3.5\n"
        "343 355 349.0 12.0 9.1 0.3\n"
        "356 366 361.0 12.0 10.7 0.4\n"
        "367 378 372.5 12.0 5.2 0.6\n"
        "379 390 384.5 12.0 7.4 0.8\n"
        "391 401 396.0 12.0 8.4 1.5\n"
        "343 355 349.0 17.0 19.9 0.3\n"
        "356 366 361.0 17.0 21.4 0.4\n"
        "367 378 372.5 17.0 18.2 0.6\n"
        "379 390 384.5 17.0 18.0 0.7\n"
        "391 401 396.0 17.0 17.0 0.7\n"
        "343 355 349.0 22.0 28.9 0.4\n"
        "356 366 361.0 22.0 30.9 0.4\n"
        "367 378 372.5 22.0 27.4 0.6\n"
        "379 390 384.5 22.0 26.7 0.6\n"
        "391 401 396.0 22.0 24.2 0.6\n"
        "343 355 349.0 27.0 36.1 0.4\n"
        "356 366 361.0 27.0 37.7 0.4\n"
        "367 378 372.5 27.0 34.1 0.6\n"
        "379 390 384.5 27.0 31.2 0.5\n"
        "391 401 396.0 27.0 29.2 0.5\n"
        "343 355 349.0 32.0 40.3 0.4\n"
        "356 366 361.0 32.0 43.3 0.5\n"
        "367 378 372.5 32.0 37.2 0.6\n"
        "379 390 384.5 32.0 34.3 0.5\n"
        "391 401 396.0 32.0 32.3 0.5\n"
        "343 355 349.0 37.0 45.1 0.4\n"
        "356 366 361.0 37.0 47.0 0.5\n"
        "367 378 372.5 37.0 40.2 0.6\n"
        "379 390 384.5 37.0 35.7 0.6\n"
        "391 401 396.0 37.0 33.1 0.6\n"
        "343 355 349.0 42.0 45.9 0.5\n"
        "356 366 361.0 42.0 46.5 0.6\n"
        "367 378 372.5 42.0 38.3 0.8\n"
        "379 390 384.5 42.0 35.7 0.8\n"
        "391 401 396.0 42.0 31.1 0.8\n"
        "343 355 349.0 47.0 43.1 0.7\n"
        "356 366 361.0 47.0 43.8 0.9\n"
        "367 378 372.5 47.0 35.6 1.1\n"
        "379 390 384.5 47.0 27.0 1.1\n"
        "391 401 396.0 47.0 28.7 1.2\n"
        "343 355 349.0 52.0 37.9 1.0\n"
        "356 366 361.0 52.0 36.8 1.3\n"
        "367 378 372.5 52.0 27.7 1.9\n"
        "379 390 384.5 52.0 27.4 2.2\n"
        "391 401 396.0 52.0 20.3 2.5\n"
        "343 355 349.0 58.0 36.4 0.8\n"
        "356 366 361.0 58.0 32.4 1.0\n"
        "367 378 372.5 58.0 20.9 1.4\n"
        "379 390 384.5 58.0 17.1 1.8\n"
        "391 401 396.0 58.0 15.0 2.8\n"
        "343 355 349.0 62.0 29.0 0.4\n"
        "356 366 361.0 62.0 27.1 0.5\n"
        "367 378 372.5 62.0 19.4 0.6\n"
        "379 390 384.5 62.0 15.5 0.7\n"
        "391 401 396.0 62.0 12.6 0.8\n"
        "343 355 349.0 67.0 23.5 0.3\n"
        "356 366 361.0 67.0 21.9 0.3\n"
        "367 378 372.5 67.0 16.6 0.4\n"
        "379 390 384.5 67.0 12.0 0.4\n"
        "391 401 396.0 67.0 10.3 0.4\n"
        "343 355 349.0 72.0 18.4 0.2\n"
        "356 366 361.0 72.0 16.7 0.2\n"
        "367 378 372.5 72.0 12.0 0.3\n"
        "379 390 384.5 72.0 9.9 0.3\n"
        "391 401 396.0 72.0 7.4 0.2\n"
        "343 355 349.0 77.0 13.3 0.2\n"
        "356 366 361.0 77.0 11.8 0.2\n"
        "367 378 372.5 77.0 8.5 0.2\n"
        "379 390 384.5 77.0 6.7 0.2\n"
        "391 401 396.0 77.0 5.4 0.2\n"
        "343 355 349.0 82.0 9.1 0.1\n"
        "356 366 361.0 82.0 7.9 0.1\n"
        "367 378 372.5 82.0 5.9 0.2\n"
        "379 390 384.5 82.0 4.7 0.2\n"
        "391 401 396.0 82.0 3.7 0.1\n"
        "343 355 349.0 87.0 6.6 0.1\n"
        "356 366 361.0 87.0 5.9 0.1\n"
        "367 378 372.5 87.0 4.1 0.1\n"
        "379 390 384.5 87.0 3.0 0.1\n"
        "391 401 396.0 87.0 2.6 0.1\n"
        "343 355 349.0 92.0 4.6 0.1\n"
        "356 366 361.0 92.0 4.0 0.1\n"
        "367 378 372.5 92.0 2.7 0.1\n"
        "379 390 384.5 92.0 2.2 0.1\n"
        "391 401 396.0 92.0 1.5 0.1\n"
        "343 355 349.0 97.0 3.0 0.1\n"
        "356 366 361.0 97.0 2.8 0.1\n"
        "367 378 372.5 97.0 1.4 0.1\n"
        "379 390 384.5 97.0 1.0 0.1\n"
        "391 401 396.0 97.0 0.8 0.1\n"
        "343 355 349.0 102.0 2.1 0.1\n"
        "356 366 361.0 102.0 1.9 0.1\n"
        "367 378 372.5 102.0 0.4 0.1\n"
        "343 355 349.0 107.0 0.2 0.1\n"
        "356 366 361.0 107.0 0.8 0.1\n"
        "367 378 372.5 107.0 0.3 0.1\n"
        "343 355 349.0 112.0 1.3 0.1\n"
        "356 366 361.0 112.0 1.0 0.1\n"
        "367 378 372.5 112.0 0.8 0.1\n"
        "343 355 349.0 117.0 1.1 0.1\n"
        "356 366 361.0 117.0 1.2 0.1\n"
        "367 378 372.5 117.0 0.4 0.1\n"
        "343 355 349.0 122.0 1.1 0.1\n"
        "356 366 361.0 122.0 1.0 0.1\n"
        "367 378 372.5 122.0 0.7 0.1\n"
        ;

    TGraphErrors* gr = new TGraphErrors();
    gr->SetName(Form("table1_dsigma_dOmega_cm_E_%.0f_%.0f", eLow, eHigh));
    gr->SetTitle(Form("Table 1 %.0f < E_{#gamma}^{lab} < %.0f MeV;#theta^{cm}_{#pi^{0}} (deg);d#sigma/d#Omega (#mub/sr)", eLow, eHigh));
    gr->SetMarkerStyle(24);
    gr->SetMarkerSize(0.9);
    gr->SetMarkerColor(kRed+1);
    gr->SetLineColor(kRed+1);
    gr->SetLineWidth(1);

    std::istringstream in(table1Data);
    std::string line;
    while(std::getline(in, line)) {
        if(line.empty()) continue;

        std::istringstream ss(line);
        double el = 0.0;
        double eh = 0.0;
        double ec = 0.0;
        double th = 0.0;
        double y = 0.0;
        double ey = 0.0;
        if(!(ss >> el >> eh >> ec >> th >> y >> ey)) continue;
        if(std::fabs(el - eLow) > 0.1 || std::fabs(eh - eHigh) > 0.1) continue;

        int ip = gr->GetN();
        gr->SetPoint(ip, th, y);
        gr->SetPointError(ip, 0.0, ey);
    }

    if(gr->GetN() == 0) {
        delete gr;
        return 0;
    }

    return gr;
}
double coherent_pion_lab_energy(double egamma, double thetaCmDeg, double mpi0, double targetMass)
{
    double s = targetMass * targetMass + 2.0 * egamma * targetMass;
    double sqrtS = std::sqrt(s);
    double ePiCm = (s + mpi0 * mpi0 - targetMass * targetMass) / (2.0 * sqrtS);
    double pPiCm2 = ePiCm * ePiCm - mpi0 * mpi0;
    if(pPiCm2 < 0.0) pPiCm2 = 0.0;
    double pPiCm = std::sqrt(pPiCm2);

    double betaCm = egamma / (egamma + targetMass);
    double gammaCm = 1.0 / std::sqrt(1.0 - betaCm * betaCm);
    double cosThetaCm = std::cos(thetaCmDeg * TMath::DegToRad());

    return gammaCm * (ePiCm + betaCm * pPiCm * cosThetaCm);
}

double phi_min_opening_angle_deg(double ePiLab, double mpi0)
{
    if(ePiLab <= mpi0) return 180.0;

    double arg = std::sqrt(ePiLab * ePiLab - mpi0 * mpi0) / ePiLab;
    if(arg < -1.0) arg = -1.0;
    if(arg > 1.0) arg = 1.0;

    return 2.0 * std::acos(arg) * TMath::RadToDeg();
}

void paper_style_diff_xs_from_acqu(const char* fname="Acqu_CBTagg_31837.root",
                             Long64_t maxEvents=-1,
                             double eps_tag=0.20,
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
                             double minEpsTag=0.02,
                             double thetaCmMinDeg=5.0,
                             double thetaCmMaxDeg=150.0,
                             double thetaCmBinWidthDeg=5.0,
                             bool usePaperBinEps=false,
                             const char* epsDetPaperFile="")
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
        std::cerr << "ERROR: could not read required FPD scaler map from "
                  << fpdFile << std::endl;
        std::cerr << "Production normalization requires the real non-contiguous "
                  << "FPD mapping; aborting instead of using 2000 + channel."
                  << std::endl;
        f->Close();
        return;
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
    } else if(taggEffFile && std::string(taggEffFile).size() > 0) {
        std::cerr << "ERROR: tagging-efficiency file was requested but could not be read: "
                  << taggEffFile << std::endl;
        f->Close();
        return;
    } else {
        std::cout << "Using explicitly requested constant eps_tag = " << eps_tag << std::endl;
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

    std::vector<double> eBinLow;
    std::vector<double> eBinHigh;
    const double eLowPaper[]  = {201, 211, 223, 235, 247, 259, 271, 283, 295, 309, 319, 331, 343, 356, 367, 379, 391};
    const double eHighPaper[] = {210, 222, 234, 246, 258, 270, 282, 294, 308, 318, 330, 342, 355, 366, 378, 390, 401};
    const int nPaperEBins = 17;
    for(int i=0; i<nPaperEBins; i++) {
        eBinLow.push_back(eLowPaper[i]);
        eBinHigh.push_back(eHighPaper[i]);
    }
    int nEBins = int(eBinLow.size());

    if(thetaCmBinWidthDeg <= 0.0) thetaCmBinWidthDeg = 5.0;
    if(thetaCmMaxDeg <= thetaCmMinDeg) {
        thetaCmMinDeg = 5.0;
        thetaCmMaxDeg = 150.0;
    }
    int nThetaBins = int(std::floor((thetaCmMaxDeg - thetaCmMinDeg) / thetaCmBinWidthDeg + 0.5));
    if(nThetaBins < 1) nThetaBins = 1;

    TString epsDetPaperPath(epsDetPaperFile);
    bool wantPaperBinEps = usePaperBinEps || epsDetPaperPath.Length() > 0;
    if(wantPaperBinEps && epsDetPaperPath.Length() == 0) {
        epsDetPaperPath = path_next_to_input(fname, "eps_det_paper_bins_from_acqu_geant.txt");
    }

    std::vector< std::vector<double> > epsDetPaper;
    std::vector< std::vector<double> > dEpsDetPaper;
    bool haveEpsDetPaperFile = false;
    if(wantPaperBinEps) {
        haveEpsDetPaperFile = load_eps_det_paper_bins(epsDetPaperPath.Data(),
                                                      eBinLow,
                                                      eBinHigh,
                                                      thetaCmMinDeg,
                                                      thetaCmBinWidthDeg,
                                                      nThetaBins,
                                                      epsDetPaper,
                                                      dEpsDetPaper);
    }
    if(haveEpsDetPaperFile) {
        std::cout << "Using paper-bin eps_det(E,theta_cm) from " << epsDetPaperPath << std::endl;
    } else if(wantPaperBinEps) {
        std::cerr << "ERROR: paper-bin eps_det was requested but no valid file was found: "
                  << epsDetPaperPath << std::endl;
        f->Close();
        return;
    } else {
        std::cout << "Using channel-dependent/fallback eps_det. Paper-bin eps_det is disabled by default." << std::endl;
    }

    TH1D* hMMp[nCh];
    TH1D* hMMr[nCh];
    std::vector< std::vector<TH1D*> > hMMpPaper;
    std::vector< std::vector<TH1D*> > hMMrPaper;
    hMMpPaper.resize(nEBins);
    hMMrPaper.resize(nEBins);

    for(int ch=0; ch<nCh; ch++) {
        hMMp[ch] = new TH1D(Form("hMM_prompt_ch%d", ch),
                            Form("prompt MM channel %d;MM-M_{^{4}He} (MeV);Counts", ch),
                            250, -100.0, 150.0);
        hMMr[ch] = new TH1D(Form("hMM_random_ch%d", ch),
                            Form("random MM channel %d;MM-M_{^{4}He} (MeV);Counts", ch),
                            250, -100.0, 150.0);

    }

    for(int ie=0; ie<nEBins; ie++) {
        hMMpPaper[ie].resize(nThetaBins);
        hMMrPaper[ie].resize(nThetaBins);
        for(int ith=0; ith<nThetaBins; ith++) {
            double th1 = thetaCmMinDeg + ith * thetaCmBinWidthDeg;
            double th2 = th1 + thetaCmBinWidthDeg;
            hMMpPaper[ie][ith] = new TH1D(Form("hMM_prompt_Ebin%d_th%d", ie, ith),
                                          Form("prompt MM %.0f-%.0f MeV, %.1f < theta_cm < %.1f deg;MM-M_{^{4}He} (MeV);Counts",
                                               eBinLow[ie], eBinHigh[ie], th1, th2),
                                          250, -100.0, 150.0);
            hMMrPaper[ie][ith] = new TH1D(Form("hMM_random_Ebin%d_th%d", ie, ith),
                                          Form("random MM %.0f-%.0f MeV, %.1f < theta_cm < %.1f deg;MM-M_{^{4}He} (MeV);Counts",
                                               eBinLow[ie], eBinHigh[ie], th1, th2),
                                          250, -100.0, 150.0);
        }
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
    TH1D* hMMAfterMM = new TH1D("h_mm_after_mgg_mm_cut",
                                "Missing mass after m_{#gamma#gamma} and MM cuts;MM-M_{^{4}He} (MeV);Prompt - w Random",
                                250, -100.0, 150.0);

    TH2D* hMggMMBefore = new TH2D("h_mgg_vs_mm_before_mgg_cut",
                                  "m_{#gamma#gamma} vs missing mass before m_{#gamma#gamma} cut;m_{#gamma#gamma} (MeV);MM-M_{^{4}He} (MeV)",
                                  250, 0.0, 250.0,
                                  250, -100.0, 150.0);
    TH2D* hMggMMAfter = new TH2D("h_mgg_vs_mm_after_mgg_cut",
                                 "m_{#gamma#gamma} vs missing mass after m_{#gamma#gamma} cut;m_{#gamma#gamma} (MeV);MM-M_{^{4}He} (MeV)",
                                 250, 0.0, 250.0,
                                 250, -100.0, 150.0);
    TH2D* hMggMMAfterMM = new TH2D("h_mgg_vs_mm_after_mgg_mm_cut",
                                   "m_{#gamma#gamma} vs missing mass after m_{#gamma#gamma} and MM cuts;m_{#gamma#gamma} (MeV);MM-M_{^{4}He} (MeV)",
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
    TH2D* hMMEgammaAfterMM = new TH2D("h_mm_vs_egamma_after_mgg_mm_cut",
                                      "Missing mass vs E_{#gamma} after m_{#gamma#gamma} and MM cuts;E_{#gamma} (MeV);MM-M_{^{4}He} (MeV)",
                                      352, 0.0, beamE,
                                      250, -100.0, 150.0);
    TH2D* hEgammaThetaBefore = new TH2D("h_egamma_vs_theta_before_mgg_cut",
                                        "E_{#gamma} vs #theta^{cm}_{#pi^{0}} before m_{#gamma#gamma} cut;#theta^{cm}_{#pi^{0}} (deg);E_{#gamma} (MeV)",
                                        180, 0.0, 180.0,
                                        352, 0.0, beamE);
    TH2D* hEgammaThetaAfter = new TH2D("h_egamma_vs_theta_after_mgg_cut",
                                       "E_{#gamma} vs #theta^{cm}_{#pi^{0}} after m_{#gamma#gamma} cut;#theta^{cm}_{#pi^{0}} (deg);E_{#gamma} (MeV)",
                                       180, 0.0, 180.0,
                                       352, 0.0, beamE);
    TH2D* hEgammaThetaAfterMM = new TH2D("h_egamma_vs_theta_after_mgg_mm_cut",
                                         "E_{#gamma} vs #theta^{cm}_{#pi^{0}} after m_{#gamma#gamma} and MM cuts;#theta^{cm}_{#pi^{0}} (deg);E_{#gamma} (MeV)",
                                         180, 0.0, 180.0,
                                         352, 0.0, beamE);
    TH2D* hOpeningPhiMinBefore = new TH2D("h_opening_angle_vs_phi_min_before_mgg_cut",
                                          "#Phi_{#gamma#gamma} vs #Phi_{min} before m_{#gamma#gamma} cut;#Phi_{min} Eq. 2 (deg);measured #Phi_{#gamma#gamma} (deg)",
                                          180, 0.0, 180.0,
                                          180, 0.0, 180.0);
    TH2D* hOpeningPhiMinAfter = new TH2D("h_opening_angle_vs_phi_min_after_mgg_cut",
                                         "#Phi_{#gamma#gamma} vs #Phi_{min} after m_{#gamma#gamma} cut;#Phi_{min} Eq. 2 (deg);measured #Phi_{#gamma#gamma} (deg)",
                                         180, 0.0, 180.0,
                                         180, 0.0, 180.0);
    TH2D* hOpeningPhiMinAfterMM = new TH2D("h_opening_angle_vs_phi_min_after_mgg_mm_cut",
                                           "#Phi_{#gamma#gamma} vs #Phi_{min} after m_{#gamma#gamma} and MM cuts;#Phi_{min} Eq. 2 (deg);measured #Phi_{#gamma#gamma} (deg)",
                                           180, 0.0, 180.0,
                                           180, 0.0, 180.0);
    TH2D* hOpeningEgammaBefore = new TH2D("h_phigg_vs_egamma_before_mgg_cut",
                                          "#Phi_{#gamma#gamma} vs E_{#gamma} before m_{#gamma#gamma} cut;E_{#gamma} (MeV);#Phi_{#gamma#gamma} (deg)",
                                          352, 0.0, beamE,
                                          180, 0.0, 180.0);
    TH2D* hOpeningEgammaAfter = new TH2D("h_phigg_vs_egamma_after_mgg_cut",
                                         "#Phi_{#gamma#gamma} vs E_{#gamma} after m_{#gamma#gamma} cut;E_{#gamma} (MeV);#Phi_{#gamma#gamma} (deg)",
                                         352, 0.0, beamE,
                                         180, 0.0, 180.0);
    TH2D* hOpeningEgammaAfterMM = new TH2D("h_phigg_vs_egamma_after_mgg_mm_cut",
                                           "#Phi_{#gamma#gamma} vs E_{#gamma} after m_{#gamma#gamma} and MM cuts;E_{#gamma} (MeV);#Phi_{#gamma#gamma} (deg)",
                                           352, 0.0, beamE,
                                           180, 0.0, 180.0);

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

    for(Long64_t ev=0; ev<nentries; ev++) {
        tracks->GetEntry(ev);
        tagger->GetEntry(ev);

        int nt = std::min(nTracks, maxTracks);
        int bestI = -1;
        int bestJ = -1;
        double bestM = -1.0;
        double bestDist = 1.0e30;
        double bestOpeningDeg = -1.0;
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
                    bestOpeningDeg = g1.Angle(g2.Vect()) * TMath::RadToDeg();
                }
            }
        }

        if(bestI < 0 || bestJ < 0) continue;
        nPairEvents++;
        hMggBest->Fill(bestM);

        bool passMggCut = (bestM >= mggMin && bestM <= mggMax);
        if(passMggCut) nPi0Events++;
        double thetaPi0LabDeg = bestPi0.Theta() * TMath::RadToDeg();

        int ntag = std::min(nTagged, maxTagged);
        for(int it=0; it<ntag; it++) {
            int ch = taggedChannel[it];
            if(ch < 0 || ch >= nCh) continue;

            double egamma = taggerPhotonEnergy[ch];
            if(egamma <= 0.0 || egamma > beamE + 5.0) continue;
            int ie = find_bin(egamma, eBinLow, eBinHigh);

            TLorentzVector beam(0.0, 0.0, egamma, egamma);
            TLorentzVector recoil = beam + target - bestPi0;
            double mm = recoil.M() - targetMass;

            TLorentzVector pi0CM = bestPi0;
            TLorentzVector totalCM = beam + target;
            pi0CM.Boost(-totalCM.BoostVector());
            double thetaCmDeg = pi0CM.Theta() * TMath::RadToDeg();
            double ePiLabCoh = coherent_pion_lab_energy(egamma, thetaCmDeg, mpi0, targetMass);
            double phiMinDeg = phi_min_opening_angle_deg(ePiLabCoh, mpi0);

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
            hEgammaThetaBefore->Fill(thetaCmDeg, egamma, weight);
            if(bestOpeningDeg >= 0.0) hOpeningPhiMinBefore->Fill(phiMinDeg, bestOpeningDeg, weight);
            if(bestOpeningDeg >= 0.0) hOpeningEgammaBefore->Fill(egamma, bestOpeningDeg, weight);

            if(!passMggCut) continue;

            hMMAfter->Fill(mm, weight);
            hMggMMAfter->Fill(bestM, mm, weight);
            hMMEgammaAfter->Fill(egamma, mm, weight);
            hEgammaThetaAfter->Fill(thetaCmDeg, egamma, weight);
            if(bestOpeningDeg >= 0.0) hOpeningPhiMinAfter->Fill(phiMinDeg, bestOpeningDeg, weight);
            bool passMMCut = (mm >= mmMin && mm <= mmMax);
            if(passMMCut) {
                hMMAfterMM->Fill(mm, weight);
                hMggMMAfterMM->Fill(bestM, mm, weight);
                hMMEgammaAfterMM->Fill(egamma, mm, weight);
                hEgammaThetaAfterMM->Fill(thetaCmDeg, egamma, weight);
            }
            if(bestOpeningDeg >= 0.0) {
                hOpeningEgammaAfter->Fill(egamma, bestOpeningDeg, weight);
                if(passMMCut) {
                    hOpeningPhiMinAfterMM->Fill(phiMinDeg, bestOpeningDeg, weight);
                    hOpeningEgammaAfterMM->Fill(egamma, bestOpeningDeg, weight);
                }
            }

            if(weight > 0.0) {
                nPromptTags++;
                hMMp[ch]->Fill(mm);
            } else {
                nRandomTags++;
                hMMr[ch]->Fill(mm);
            }

            if(ie >= 0 && thetaCmDeg >= thetaCmMinDeg && thetaCmDeg < thetaCmMaxDeg) {
                int ith = int((thetaCmDeg - thetaCmMinDeg) / thetaCmBinWidthDeg);
                if(ith >= 0 && ith < nThetaBins) {
                    if(weight > 0.0) {
                        hMMpPaper[ie][ith]->Fill(mm);
                    } else {
                        hMMrPaper[ie][ith]->Fill(mm);
                    }
                }
            }
        }
    }

    std::ofstream out("paper_style_dsigma_dOmega_cm_from_acqu.txt");
    out << "# input: " << fname << "\n";
    out << "# mode: paper-style differential cross section, theta in CM, photon energy bins in lab\n";
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
    out << "# eps_det paper-bin file = " << epsDetPaperPath << "\n";
    out << "# use eps_det paper-bin requested = " << wantPaperBinEps << "\n";
    out << "# have eps_det paper-bin file = " << haveEpsDetPaperFile << "\n";
    out << "# min eps_det = " << minEpsDet << "\n";
    out << "# thickness = " << thickness << " ub^-1\n";
    out << "# nTagger = " << nTagger << "\n";
    out << "# minScalerSum = " << minScalerSum << "\n";
    out << "# energy bins: same style as Table 1 of arXiv:nucl-ex/9907020, photon laboratory energy in MeV\n";
    out << "# theta_cm binning = [" << thetaCmMinDeg << "," << thetaCmMaxDeg << "] deg, width = " << thetaCmBinWidthDeg << " deg, nThetaBins = " << nThetaBins << "\n";
    out << "# E_low E_high E_center thetaCm_low thetaCm_high thetaCm_center DeltaOmega_cm_sr fluxDen eps_det deps_det Y dY dsigma_dOmega_ub_per_sr ddsigma_dOmega_ub_per_sr prompt random status\n";

    std::vector<double> normDenE(nEBins, 0.0);
    std::vector<double> normErr2E(nEBins, 0.0);

    for(int ch=0; ch<nValidCh; ch++) {
        int ie = find_bin(taggerPhotonEnergy[ch], eBinLow, eBinHigh);
        if(ie < 0) continue;

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

        double dNe = std::sqrt(dNe2[ch]);
        double relNorm2 = 0.0;
        if(Ne[ch] > 0.0 && dNe > 0.0) relNorm2 += (dNe/Ne[ch]) * (dNe/Ne[ch]);
        if(epsTagUse > 0.0 && dEpsTagUse > 0.0) relNorm2 += (dEpsTagUse/epsTagUse) * (dEpsTagUse/epsTagUse);
        if(!haveEpsDetPaperFile && epsDetUse > 0.0 && dEpsDetUse > 0.0) relNorm2 += (dEpsDetUse/epsDetUse) * (dEpsDetUse/epsDetUse);

        bool validNorm = (Ne[ch] >= minScalerSum &&
                          epsTagUse >= minEpsTag &&
                          thickness > 0.0);
        if(!haveEpsDetPaperFile) validNorm = validNorm && (epsDetUse >= minEpsDet);

        double normDen = 0.0;
        if(validNorm) {
            normDen = Ne[ch] * epsTagUse * thickness;
            if(!haveEpsDetPaperFile) normDen *= epsDetUse;
        }
        if(validNorm) {
            normDenE[ie] += normDen;
            normErr2E[ie] += normDen * normDen * relNorm2;
        }
    }

    std::vector<TGraphErrors*> grPaper(nEBins, 0);
    std::vector<TGraphErrors*> grPhiMin(nEBins, 0);
    TH2D* hPaperMap = new TH2D("h_dsigma_dOmega_cm_paper_bins",
                               "d#sigma/d#Omega in CM, paper-style bins;E_{#gamma}^{lab} bin;#theta^{cm}_{#pi^{0}} (deg)",
                               nEBins, 0.0, double(nEBins),
                               nThetaBins, thetaCmMinDeg, thetaCmMinDeg + nThetaBins * thetaCmBinWidthDeg);
    TH2D* hPhiMinMap = new TH2D("h_phi_min_eq2_cm_bins",
                                "#Phi_{min} from Eq. 2 using coherent two-body kinematics;E_{#gamma}^{lab} bin;#theta^{cm}_{#pi^{0}} (deg)",
                                nEBins, 0.0, double(nEBins),
                                nThetaBins, thetaCmMinDeg, thetaCmMinDeg + nThetaBins * thetaCmBinWidthDeg);

    for(int ie=0; ie<nEBins; ie++) {
        grPaper[ie] = new TGraphErrors();
        grPaper[ie]->SetName(Form("dsigma_dOmega_cm_E_%.0f_%.0f", eBinLow[ie], eBinHigh[ie]));
        grPaper[ie]->SetTitle(Form("%.0f < E_{#gamma}^{lab} < %.0f MeV;#theta^{cm}_{#pi^{0}} (deg);d#sigma/d#Omega (#mub/sr)", eBinLow[ie], eBinHigh[ie]));
        grPaper[ie]->SetMarkerStyle(20);
        grPaper[ie]->SetMarkerSize(0.9);
        grPaper[ie]->SetMarkerColor(kBlue+1);
        grPaper[ie]->SetLineColor(kBlue+1);
        grPhiMin[ie] = new TGraphErrors();
        grPhiMin[ie]->SetName(Form("phi_min_eq2_E_%.0f_%.0f", eBinLow[ie], eBinHigh[ie]));
        grPhiMin[ie]->SetTitle(Form("#Phi_{min}, %.0f < E_{#gamma}^{lab} < %.0f MeV;#theta^{cm}_{#pi^{0}} (deg);#Phi_{min} Eq. 2 (deg)", eBinLow[ie], eBinHigh[ie]));
        grPhiMin[ie]->SetMarkerStyle(21);
        grPhiMin[ie]->SetMarkerSize(0.8);
        grPhiMin[ie]->SetMarkerColor(kMagenta+2);
        grPhiMin[ie]->SetLineColor(kMagenta+2);

        for(int ith=0; ith<nThetaBins; ith++) {
            double th1 = thetaCmMinDeg + ith * thetaCmBinWidthDeg;
            double th2 = th1 + thetaCmBinWidthDeg;
            double thCenter = 0.5 * (th1 + th2);
            double deltaOmega = solid_angle_bin(th1, th2);
            double eCenter = 0.5 * (eBinLow[ie] + eBinHigh[ie]);
            double ePiLabCoh = coherent_pion_lab_energy(eCenter, thCenter, mpi0, targetMass);
            double phiMinDeg = phi_min_opening_angle_deg(ePiLabCoh, mpi0);
            hPhiMinMap->SetBinContent(ie + 1, ith + 1, phiMinDeg);
            grPhiMin[ie]->SetPoint(ith, thCenter, phiMinDeg);
            grPhiMin[ie]->SetPointError(ith, 0.5 * thetaCmBinWidthDeg, 0.0);

            int b1 = hMMpPaper[ie][ith]->GetXaxis()->FindBin(mmMin);
            int b2 = hMMpPaper[ie][ith]->GetXaxis()->FindBin(mmMax);

            double ep = 0.0;
            double er = 0.0;
            double p = hMMpPaper[ie][ith]->IntegralAndError(b1, b2, ep);
            double r = hMMrPaper[ie][ith]->IntegralAndError(b1, b2, er);

            double Y = p - randomWeight * r;
            double dY = std::sqrt(ep*ep + randomWeight*randomWeight*er*er);
            double dsigma_dOmega = 0.0;
            double ddsigma_dOmega = 0.0;
            double epsDetBin = 1.0;
            double dEpsDetBin = 0.0;
            const char* status = "ok";

            if(normDenE[ie] <= 0.0) status = "bad_energy_norm";
            else if(deltaOmega <= 0.0) status = "bad_deltaOmega";
            else if(haveEpsDetPaperFile) {
                epsDetBin = epsDetPaper[ie][ith];
                dEpsDetBin = dEpsDetPaper[ie][ith];
                if(epsDetBin <= 0.0) status = "bad_eps_det_paper_zero";
                else if(epsDetBin < minEpsDet) status = "bad_eps_det_paper_low";
            }

            if(status[0] == 'o' && status[1] == 'k' && status[2] == '\0') {
                double denom = normDenE[ie] * epsDetBin * deltaOmega;
                dsigma_dOmega = Y / denom;
                double statErr = dY / denom;
                double relNorm2 = normErr2E[ie] / (normDenE[ie] * normDenE[ie]);
                if(haveEpsDetPaperFile && epsDetBin > 0.0 && dEpsDetBin > 0.0) {
                    relNorm2 += (dEpsDetBin / epsDetBin) * (dEpsDetBin / epsDetBin);
                }
                double normErr = std::fabs(dsigma_dOmega) * std::sqrt(relNorm2);
                ddsigma_dOmega = std::sqrt(statErr*statErr + normErr*normErr);

                int ip = grPaper[ie]->GetN();
                grPaper[ie]->SetPoint(ip, thCenter, dsigma_dOmega);
                grPaper[ie]->SetPointError(ip, 0.5 * thetaCmBinWidthDeg, ddsigma_dOmega);
                hPaperMap->SetBinContent(ie + 1, ith + 1, dsigma_dOmega);
                hPaperMap->SetBinError(ie + 1, ith + 1, ddsigma_dOmega);
            }

            out << eBinLow[ie] << " "
                << eBinHigh[ie] << " "
                << 0.5 * (eBinLow[ie] + eBinHigh[ie]) << " "
                << th1 << " "
                << th2 << " "
                << thCenter << " "
                << deltaOmega << " "
                << normDenE[ie] << " "
                << epsDetBin << " "
                << dEpsDetBin << " "
                << Y << " "
                << dY << " "
                << dsigma_dOmega << " "
                << ddsigma_dOmega << " "
                << p << " "
                << r << " "
                << status << "\n";
        }
    }
    out.close();

    const int nCompareBins = nEBins;
    std::vector<TGraphErrors*> grTable1(nCompareBins, 0);
    for(int ie=0; ie<nCompareBins; ie++) {
        grTable1[ie] = make_table1_graph(ie, eBinLow[ie], eBinHigh[ie]);
    }

    TH2D* hTable1Map = new TH2D("h_dsigma_dOmega_cm_table1_all",
                                "Table 1 d#sigma/d#Omega in CM, all energy bins;E_{#gamma}^{lab} bin;#theta^{cm}_{#pi^{0}} (deg)",
                                nCompareBins, 0.0, double(nCompareBins),
                                nThetaBins, thetaCmMinDeg, thetaCmMinDeg + nThetaBins * thetaCmBinWidthDeg);
    for(int ie=0; ie<nCompareBins; ie++) {
        hTable1Map->GetXaxis()->SetBinLabel(ie + 1, Form("%.0f-%.0f", eBinLow[ie], eBinHigh[ie]));
        if(!grTable1[ie]) continue;
        for(int ip=0; ip<grTable1[ie]->GetN(); ip++) {
            double x = 0.0;
            double y = 0.0;
            grTable1[ie]->GetPoint(ip, x, y);
            int by = hTable1Map->GetYaxis()->FindBin(x);
            hTable1Map->SetBinContent(ie + 1, by, y);
            hTable1Map->SetBinError(ie + 1, by, grTable1[ie]->GetErrorY(ip));
        }
    }

    TCanvas* c1 = new TCanvas("c_mgg_from_acqu", "mgg from Acqu", 800, 600);
    hMggBest->Draw();
    c1->SaveAs("mgg_best_from_acqu.pdf");

    TCanvas* cMM = new TCanvas("c_mm_before_after_from_acqu", "missing mass before/after", 1500, 500);
    cMM->Divide(3, 1);
    cMM->cd(1);
    hMMBefore->Draw("HIST");
    cMM->cd(2);
    hMMAfter->Draw("HIST");
    cMM->cd(3);
    hMMAfterMM->Draw("HIST");
    cMM->SaveAs("missing_mass_before_after_from_acqu.pdf");

    TCanvas* c2D = new TCanvas("c_mgg_vs_mm_before_after_from_acqu", "mgg vs missing mass before/after", 1500, 500);
    c2D->Divide(3, 1);
    c2D->cd(1);
    hMggMMBefore->Draw("COLZ");
    c2D->cd(2);
    hMggMMAfter->Draw("COLZ");
    c2D->cd(3);
    hMggMMAfterMM->Draw("COLZ");
    c2D->SaveAs("mgg_vs_missing_mass_before_after_from_acqu.pdf");

    TCanvas* cMME = new TCanvas("c_mm_vs_egamma_before_after_from_acqu", "missing mass vs Egamma before/after", 1500, 500);
    cMME->Divide(3, 1);
    cMME->cd(1);
    hMMEgammaBefore->Draw("COLZ");
    cMME->cd(2);
    hMMEgammaAfter->Draw("COLZ");
    cMME->cd(3);
    hMMEgammaAfterMM->Draw("COLZ");
    cMME->SaveAs("missing_mass_vs_egamma_before_after_from_acqu.pdf");

    TCanvas* cETheta = new TCanvas("c_egamma_vs_theta_before_after_from_acqu", "Egamma vs theta before/after", 1500, 500);
    cETheta->Divide(3, 1);
    cETheta->cd(1);
    hEgammaThetaBefore->Draw("COLZ");
    cETheta->cd(2);
    hEgammaThetaAfter->Draw("COLZ");
    cETheta->cd(3);
    hEgammaThetaAfterMM->Draw("COLZ");
    cETheta->SaveAs("egamma_vs_theta_before_after_from_acqu.pdf");

    TCanvas* cPaper = new TCanvas("c_paper_style_dsigma_dOmega_cm", "paper-style dsigma/dOmega CM", 1100, 700);
    cPaper->Print("paper_style_dsigma_dOmega_cm_by_energy.pdf[");
    for(int ie=0; ie<nEBins; ie++) {
        cPaper->Clear();
        cPaper->SetGrid();
        grPaper[ie]->Draw("AP");
        cPaper->Print("paper_style_dsigma_dOmega_cm_by_energy.pdf");
    }
    cPaper->Print("paper_style_dsigma_dOmega_cm_by_energy.pdf]");

    TCanvas* cPaperMap = new TCanvas("c_paper_style_dsigma_dOmega_cm_map", "paper-style dsigma/dOmega CM map", 1200, 750);
    for(int ie=0; ie<nEBins; ie++) {
        hPaperMap->GetXaxis()->SetBinLabel(ie + 1, Form("%.0f-%.0f", eBinLow[ie], eBinHigh[ie]));
    }
    hPaperMap->GetXaxis()->LabelsOption("v");
    hPaperMap->Draw("COLZ");
    cPaperMap->SaveAs("paper_style_dsigma_dOmega_cm_map.pdf");

    std::vector<TCanvas*> cOverlayPages;
    TCanvas* cOverlayBook = new TCanvas("c_paper_style_dsigma_dOmega_cm_grouped6_overlay_book",
                                        "paper-style bins with Table 1 overlay, grouped by six", 1500, 900);
    cOverlayBook->Print("paper_style_dsigma_dOmega_cm_grouped6_overlay.pdf[");
    for(int first=0; first<nCompareBins; first+=6) {
        int group = first / 6;
        TCanvas* cOverlayPage = new TCanvas(Form("c_paper_style_dsigma_dOmega_cm_grouped6_overlay_page%d", group + 1),
                                            Form("paper-style bins with Table 1 overlay, page %d", group + 1),
                                            1500, 900);
        cOverlayPage->Divide(3, 2);
        for(int ipad=0; ipad<6; ipad++) {
            int ie = first + ipad;
            if(ie >= nCompareBins) break;

            cOverlayPage->cd(ipad + 1);
            gPad->SetGrid();
            grPaper[ie]->SetMarkerStyle(20);
            grPaper[ie]->SetMarkerColor(kBlue+1);
            grPaper[ie]->SetLineColor(kBlue+1);
            grPaper[ie]->Draw("AP");
            grPaper[ie]->GetXaxis()->SetLimits(0.0, 160.0);
            grPaper[ie]->GetYaxis()->SetTitleOffset(1.25);
            if(grTable1[ie]) grTable1[ie]->Draw("P SAME");
            TLegend* leg = new TLegend(0.52, 0.72, 0.88, 0.88);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->AddEntry(grPaper[ie], "ACQU", "p");
            if(grTable1[ie]) leg->AddEntry(grTable1[ie], "Paper Table 1", "p");
            leg->Draw();
        }
        cOverlayPage->Print("paper_style_dsigma_dOmega_cm_grouped6_overlay.pdf");
        cOverlayPages.push_back(cOverlayPage);
    }
    cOverlayBook->Print("paper_style_dsigma_dOmega_cm_grouped6_overlay.pdf]");

    TCanvas* cTable1Map = new TCanvas("c_paper_table1_dsigma_dOmega_cm_map", "Table 1 dsigma/dOmega CM map", 1200, 750);
    hTable1Map->GetXaxis()->LabelsOption("v");
    hTable1Map->Draw("COLZ");
    cTable1Map->SaveAs("paper_table1_dsigma_dOmega_cm_map.pdf");

    TCanvas* cPhiMinMap = new TCanvas("c_phi_min_eq2_cm_map", "Phi_min Eq. 2 CM map", 1200, 750);
    for(int ie=0; ie<nEBins; ie++) {
        hPhiMinMap->GetXaxis()->SetBinLabel(ie + 1, Form("%.0f-%.0f", eBinLow[ie], eBinHigh[ie]));
    }
    hPhiMinMap->GetXaxis()->LabelsOption("v");
    hPhiMinMap->Draw("COLZ");
    cPhiMinMap->SaveAs("phi_min_eq2_cm_map.pdf");

    TCanvas* cOpeningPhiMin = new TCanvas("c_opening_angle_vs_phi_min", "opening angle vs Phi_min", 1500, 500);
    cOpeningPhiMin->Divide(3, 1);
    cOpeningPhiMin->cd(1);
    hOpeningPhiMinBefore->Draw("COLZ");
    cOpeningPhiMin->cd(2);
    hOpeningPhiMinAfter->Draw("COLZ");
    cOpeningPhiMin->cd(3);
    hOpeningPhiMinAfterMM->Draw("COLZ");
    cOpeningPhiMin->SaveAs("opening_angle_vs_phi_min_before_after.pdf");

    TCanvas* cOpeningPhiMinAfterMM = new TCanvas("c_opening_angle_vs_phi_min_after_mgg_mm",
                                                 "opening angle vs Phi_min after mgg and MM cuts",
                                                 900, 700);
    cOpeningPhiMinAfterMM->SetRightMargin(0.14);
    hOpeningPhiMinAfterMM->Draw("COLZ");
    cOpeningPhiMinAfterMM->SaveAs("opening_angle_vs_phi_min_after_mgg_mm.pdf");

    TCanvas* cOpeningEgamma = new TCanvas("c_phigg_vs_egamma_from_acqu", "Phi_gg vs Egamma", 1500, 500);
    cOpeningEgamma->Divide(3, 1);
    cOpeningEgamma->cd(1);
    hOpeningEgammaBefore->Draw("COLZ");
    cOpeningEgamma->cd(2);
    hOpeningEgammaAfter->Draw("COLZ");
    cOpeningEgamma->cd(3);
    hOpeningEgammaAfterMM->Draw("COLZ");
    cOpeningEgamma->SaveAs("phigg_vs_egamma_before_after_from_acqu.pdf");

    TCanvas* cOpeningEgammaAfterMM = new TCanvas("c_phigg_vs_egamma_after_mgg_mm_from_acqu",
                                                 "Phi_gg vs Egamma after mgg and MM cuts",
                                                 900, 700);
    cOpeningEgammaAfterMM->SetRightMargin(0.14);
    hOpeningEgammaAfterMM->Draw("COLZ");
    cOpeningEgammaAfterMM->SaveAs("phigg_vs_egamma_after_mgg_mm_from_acqu.pdf");

    TCanvas* cPhiMin = new TCanvas("c_phi_min_eq2_by_energy", "Phi_min Eq. 2 by energy", 1100, 700);
    cPhiMin->Print("phi_min_eq2_by_energy.pdf[");
    for(int ie=0; ie<nEBins; ie++) {
        cPhiMin->Clear();
        cPhiMin->SetGrid();
        grPhiMin[ie]->Draw("ALP");
        cPhiMin->Print("phi_min_eq2_by_energy.pdf");
    }
    cPhiMin->Print("phi_min_eq2_by_energy.pdf]");

    TFile* fout = new TFile("paper_style_dsigma_dOmega_cm_from_acqu.root", "RECREATE");
    hMggBest->Write();
    hMMBefore->Write();
    hMMAfter->Write();
    hMMAfterMM->Write();
    hMggMMBefore->Write();
    hMggMMAfter->Write();
    hMggMMAfterMM->Write();
    hMMEgammaBefore->Write();
    hMMEgammaAfter->Write();
    hMMEgammaAfterMM->Write();
    hEgammaThetaBefore->Write();
    hEgammaThetaAfter->Write();
    hEgammaThetaAfterMM->Write();
    hOpeningPhiMinBefore->Write();
    hOpeningPhiMinAfter->Write();
    hOpeningPhiMinAfterMM->Write();
    hOpeningEgammaBefore->Write();
    hOpeningEgammaAfter->Write();
    hOpeningEgammaAfterMM->Write();
    for(int ch=0; ch<nCh; ch++) {
        hMMp[ch]->Write();
        hMMr[ch]->Write();
    }
    for(int ie=0; ie<nEBins; ie++) {
        for(int ith=0; ith<nThetaBins; ith++) {
            hMMpPaper[ie][ith]->Write();
            hMMrPaper[ie][ith]->Write();
        }
        grPaper[ie]->Write();
        grPhiMin[ie]->Write();
    }
    for(int ie=0; ie<nCompareBins; ie++) {
        if(grTable1[ie]) grTable1[ie]->Write();
    }
    hPaperMap->Write();
    hTable1Map->Write();
    hPhiMinMap->Write();
    c1->Write();
    cMM->Write();
    c2D->Write();
    cMME->Write();
    cETheta->Write();
    cPaper->Write();
    cPaperMap->Write();
    cOverlayBook->Write();
    for(int i=0; i<int(cOverlayPages.size()); i++) {
        cOverlayPages[i]->Write();
    }
    cTable1Map->Write();
    cPhiMinMap->Write();
    cOpeningPhiMin->Write();
    cOpeningPhiMinAfterMM->Write();
    cOpeningEgamma->Write();
    cOpeningEgammaAfterMM->Write();
    cPhiMin->Write();
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
    std::cout << "  paper_style_dsigma_dOmega_cm_from_acqu.txt" << std::endl;
    std::cout << "  paper_style_dsigma_dOmega_cm_by_energy.pdf" << std::endl;
    std::cout << "  paper_style_dsigma_dOmega_cm_map.pdf" << std::endl;
    std::cout << "  paper_style_dsigma_dOmega_cm_grouped6_overlay.pdf" << std::endl;
    std::cout << "  paper_table1_dsigma_dOmega_cm_map.pdf" << std::endl;
    std::cout << "  phi_min_eq2_cm_map.pdf" << std::endl;
    std::cout << "  phi_min_eq2_by_energy.pdf" << std::endl;
    std::cout << "  opening_angle_vs_phi_min_before_after.pdf" << std::endl;
    std::cout << "  opening_angle_vs_phi_min_after_mgg_mm.pdf" << std::endl;
    std::cout << "  phigg_vs_egamma_before_after_from_acqu.pdf" << std::endl;
    std::cout << "  phigg_vs_egamma_after_mgg_mm_from_acqu.pdf" << std::endl;
    std::cout << "  mgg_best_from_acqu.pdf" << std::endl;
    std::cout << "  missing_mass_before_after_from_acqu.pdf" << std::endl;
    std::cout << "  mgg_vs_missing_mass_before_after_from_acqu.pdf" << std::endl;
    std::cout << "  missing_mass_vs_egamma_before_after_from_acqu.pdf" << std::endl;
    std::cout << "  egamma_vs_theta_before_after_from_acqu.pdf" << std::endl;
    std::cout << "  paper_style_dsigma_dOmega_cm_from_acqu.root" << std::endl;

    f->Close();
}

void paper_style_diff_xs_from_acqu_goat_eps(const char* fname="Acqu_CBTagg_31837.root",
                                            const char* epsDetPaperFile="eps_det_paper_bins_from_goat.txt",
                                            const char* fpdFile="FPD_855_new.dat",
                                            const char* taggEffFile="ExpBkgSub_COPP_TaggEff_31834.dat")
{
    paper_style_diff_xs_from_acqu(fname,
                                  -1,
                                  0.20,
                                  0.7064,
                                  0.940e-7,
                                  true,
                                  1.0e7,
                                  fpdFile,
                                  700.0,
                                  800.0,
                                  450.0,
                                  680.0,
                                  -10.0,
                                  40.0,
                                  110.0,
                                  155.0,
                                  "",
                                  1.0,
                                  "",
                                  0.02,
                                  taggEffFile,
                                  0.02,
                                  5.0,
                                  150.0,
                                  5.0,
                                  true,
                                  epsDetPaperFile);
}
