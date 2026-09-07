/*
 * cross_section_from_acqu.C
 *
 * PURPOSE
 *   Extract the total 4He(gamma,pi0)4He cross section independently for each
 *   tagger channel from an AcquRoot file.
 *
 * EVENT SELECTION
 *   - Reads tracks, tagger, scalers and setupParameters.
 *   - Accepts clusters with E >= 20 MeV and vetoEnergy <= 1 MeV.
 *   - Chooses the neutral-cluster pair whose invariant mass is closest to the
 *     pi0 mass.
 *   - Applies the configurable gamma-gamma invariant-mass window (default
 *     110-155 MeV).
 *   - Fills prompt and random missing-mass histograms for each tagger channel.
 *   - Integrates the configurable missing-mass window (default -10 to 40 MeV).
 *
 * YIELD AND CROSS SECTION
 *   Y = N_prompt - w_random * N_random
 *   sigma = Y / (Ne * eps_tag * eps_det * target_thickness)
 *
 *   Ne is obtained from the tagger electron scalers. Channel-to-scaler mapping
 *   is read from the required FPD configuration file. The validated mapping is
 *   non-contiguous, so production extraction aborts if the map cannot be read.
 *
 * EFFICIENCIES
 *   Tagging efficiency:
 *     - channel dependent when taggEffFile is supplied;
 *     - if a supplied file cannot be read, the macro aborts;
 *     - pass an empty taggEffFile explicitly to use constant eps_tag.
 *
 *   Detection efficiency:
 *     - channel dependent when epsDetFile can be read;
 *     - otherwise the constant eps_det argument is used.
 *     - when epsDetFile is empty, the macro looks next to the input ROOT file
 *       for 'eps_det_reco_fast_from_acqu_geant.txt'. No generator for that
 *       exact historical file was found in the supplied macro archives, so
 *       this dependency must be supplied externally or the constant fallback
 *       is used.
 *
 * IMPORTANT
 *   This macro does NOT apply the DeltaPhi/opening-angle cut used in the
 *   separate Fig. 2 coherent-selection study.
 *
 * OUTPUTS
 *   sigma_from_acqu.txt
 *   sigma_from_acqu.pdf
 *   sigma_from_acqu_vs_egamma.pdf
 *   sigma_from_acqu_vs_egamma_with_reference.pdf (if reference data loaded)
 *   mgg_best_from_acqu.pdf
 *   cross_section_from_acqu.root
 *
 * USAGE EXAMPLE
 *   root -l -b -q 'cross_section_from_acqu.C("Acqu_CBTagg_31837.root")'
 *
 * NOTE FOR THE REPOSITORY
 *   Physics selections and the public function interface are preserved.
 *   Repository preparation removed machine-specific default paths and makes the real
 *   FPD map / explicitly requested tagging-efficiency input fail closed.
 */
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TMath.h"
#include "TLorentzVector.h"
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

void cross_section_from_acqu(const char* fname="Acqu_CBTagg_31837.root",
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

    Long64_t nentries = std::min(tracks->GetEntries(), tagger->GetEntries());
    if(maxEvents > 0 && maxEvents < nentries) nentries = maxEvents;

    const double deg = TMath::DegToRad();
    TLorentzVector target(0.0, 0.0, 0.0, targetMass);

    Long64_t nPairEvents = 0;
    Long64_t nPi0Events = 0;
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

        if(bestM < mggMin || bestM > mggMax) continue;
        nPi0Events++;

        int ntag = std::min(nTagged, maxTagged);
        for(int it=0; it<ntag; it++) {
            int ch = taggedChannel[it];
            if(ch < 0 || ch >= nCh) continue;

            double egamma = taggerPhotonEnergy[ch];
            if(egamma <= 0.0 || egamma > beamE + 5.0) continue;

            TLorentzVector beam(0.0, 0.0, egamma, egamma);
            TLorentzVector recoil = beam + target - bestPi0;
            double mm = recoil.M() - targetMass;

            double tt = taggedTime[it];
            if(tt > promptMin && tt < promptMax) {
                hMMp[ch]->Fill(mm);
                nPromptTags++;
            } else if(tt > randomMin && tt < randomMax) {
                hMMr[ch]->Fill(mm);
                nRandomTags++;
            }
        }
    }

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

    TFile* fout = new TFile("cross_section_from_acqu.root", "RECREATE");
    hMggBest->Write();
    for(int ch=0; ch<nCh; ch++) {
        hMMp[ch]->Write();
        hMMr[ch]->Write();
    }
    gr->Write();
    grE->Write();
    if(grRef) grRef->Write();
    c1->Write();
    c2->Write();
    c3->Write();
    fout->Close();

    std::cout << "Input file: " << fname << std::endl;
    std::cout << "Processed events = " << nentries << std::endl;
    std::cout << "Events with gamma-gamma pair = " << nPairEvents << std::endl;
    std::cout << "Events passing pi0 mass cut = " << nPi0Events << std::endl;
    std::cout << "Prompt tags used = " << nPromptTags << std::endl;
    std::cout << "Random tags used = " << nRandomTags << std::endl;
    std::cout << "Random weight = " << randomWeight << std::endl;
    std::cout << "Saved:" << std::endl;
    std::cout << "  sigma_from_acqu.txt" << std::endl;
    std::cout << "  sigma_from_acqu.pdf" << std::endl;
    std::cout << "  sigma_from_acqu_vs_egamma.pdf" << std::endl;
    if(grRef) std::cout << "  sigma_from_acqu_vs_egamma_with_reference.pdf" << std::endl;
    std::cout << "  mgg_best_from_acqu.pdf" << std::endl;
    std::cout << "  cross_section_from_acqu.root" << std::endl;

    f->Close();
}
