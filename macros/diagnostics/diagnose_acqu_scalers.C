/*
 * diagnose_acqu_scalers.C
 *
 * Purpose
 * -------
 * Diagnose how Acqu scaler entries should be summed for the tagger electron
 * counts. In particular, it compares sums including scaler-tree entry 0 with
 * sums that skip entry 0.
 *
 * Input
 * -----
 * An Acqu ROOT file containing the trees "scalers" and "setupParameters".
 * The macro assumes a contiguous block of tagger scalers starting at
 * firstScaler (default 2000) and inspects nCh channels (default 352).
 *
 * Output
 * ------
 * diagnose_acqu_scalers.txt   per-channel sums and neighbour ratios
 * diagnose_acqu_scalers.png   comparison plot (logarithmic y axis)
 * diagnose_acqu_scalers.root  graphs and canvas
 *
 * Usage
 * -----
 * root -l -b -q 'diagnose_acqu_scalers.C("Acqu_CBTagg_31838.root",2000,352)'
 *
 * This is a diagnostic macro, not a production normalization step.
 * The code below is the original analysis code; only this documentation
 * header has been added.
 */

#include <cmath>
#include <fstream>
#include <iostream>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TString.h"
#include "TTree.h"

void diagnose_acqu_scalers(const char* fname="Acqu_CBTagg_31838.root",
                           int firstScaler=2000,
                           int nCh=352)
{
    TFile* f = TFile::Open(fname);
    if(!f || f->IsZombie()) {
        std::cout << "Cannot open " << fname << std::endl;
        return;
    }

    TTree* scalers = (TTree*)f->Get("scalers");
    TTree* setup = (TTree*)f->Get("setupParameters");
    if(!scalers || !setup) {
        std::cout << "Missing scalers or setupParameters tree" << std::endl;
        f->Close();
        return;
    }

    UInt_t sc[8706];
    scalers->SetBranchAddress("scalers", sc);

    int nTagger = 0;
    double egamma[512];
    setup->SetBranchAddress("nTagger", &nTagger);
    setup->SetBranchAddress("TaggerPhotonEnergy", egamma);
    setup->GetEntry(0);

    double neIncl[512] = {0.0};
    double neExcl[512] = {0.0};

    for(Long64_t is=0; is<scalers->GetEntries(); is++) {
        scalers->GetEntry(is);
        for(int ch=0; ch<nCh; ch++) {
            double val = sc[firstScaler + ch];
            neIncl[ch] += val;
            if(is > 0) neExcl[ch] += val;
        }
    }

    std::ofstream out("diagnose_acqu_scalers.txt");
    out << "# input: " << fname << "\n";
    out << "# firstScaler = " << firstScaler << "\n";
    out << "# channel Egamma Ne_including_entry0 Ne_excluding_entry0 ratio_to_previous_excl\n";

    TGraph* gIncl = new TGraph();
    TGraph* gExcl = new TGraph();
    gIncl->SetName("Ne_including_entry0");
    gExcl->SetName("Ne_excluding_entry0");
    gIncl->SetTitle("Acqu scaler sums;Tagger channel;Scaler sum");
    gExcl->SetTitle("Acqu scaler sums;Tagger channel;Scaler sum");
    gIncl->SetLineColor(kRed+1);
    gIncl->SetMarkerColor(kRed+1);
    gIncl->SetMarkerStyle(24);
    gExcl->SetLineColor(kBlue+1);
    gExcl->SetMarkerColor(kBlue+1);
    gExcl->SetMarkerStyle(20);

    std::cout << "Suspicious jumps in Ne_excluding_entry0:" << std::endl;
    for(int ch=0; ch<nCh; ch++) {
        double ratio = 0.0;
        if(ch > 0 && neExcl[ch-1] > 0.0) ratio = neExcl[ch] / neExcl[ch-1];

        out << ch << " "
            << egamma[ch] << " "
            << neIncl[ch] << " "
            << neExcl[ch] << " "
            << ratio << "\n";

        gIncl->SetPoint(ch, ch, neIncl[ch]);
        gExcl->SetPoint(ch, ch, neExcl[ch]);

        if(ch > 0 && (ratio < 0.25 || ratio > 4.0)) {
            std::cout << "  ch " << ch-1 << " -> " << ch
                      << " : " << neExcl[ch-1]
                      << " -> " << neExcl[ch]
                      << "  ratio = " << ratio
                      << std::endl;
        }
    }
    out.close();

    TCanvas* c = new TCanvas("c_diagnose_acqu_scalers", "Acqu scaler diagnosis", 1100, 700);
    c->SetGrid();
    c->SetLogy();
    gIncl->Draw("APL");
    gExcl->Draw("PL SAME");

    TLegend* leg = new TLegend(0.58, 0.75, 0.88, 0.88);
    leg->AddEntry(gIncl, "including entry 0", "lp");
    leg->AddEntry(gExcl, "excluding entry 0", "lp");
    leg->Draw();

    c->SaveAs("diagnose_acqu_scalers.png");

    TFile* fout = new TFile("diagnose_acqu_scalers.root", "RECREATE");
    gIncl->Write();
    gExcl->Write();
    c->Write();
    fout->Close();

    std::cout << "Saved:" << std::endl;
    std::cout << "  diagnose_acqu_scalers.txt" << std::endl;
    std::cout << "  diagnose_acqu_scalers.png" << std::endl;
    std::cout << "  diagnose_acqu_scalers.root" << std::endl;

    f->Close();
}
