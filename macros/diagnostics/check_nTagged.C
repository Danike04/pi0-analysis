/*
 * check_nTagged.C
 *
 * Purpose
 * -------
 * Inspect the event-by-event tagger multiplicity nTagged in a GoAT ROOT
 * file. The macro fills the nTagged distribution and prints its mean, RMS
 * and cumulative event counts below several multiplicity thresholds.
 *
 * Input
 * -----
 * The filename is hard-coded in the original macro as
 * GoAT_CBTagg_31838.root. The file must contain the "tagger" tree with the
 * nTagged branch.
 *
 * Output
 * ------
 * nTagged_distribution.pdf plus summary values printed to stdout.
 *
 * Usage
 * -----
 * Place/run the macro where GoAT_CBTagg_31838.root is available, then use
 *
 *     root -l -b -q check_nTagged.C
 *
 * This is a historical diagnostic, not a production analysis step.
 * The code below is the original analysis code; only this documentation
 * header has been added.
 */

void check_nTagged() {

    const char* filename = "GoAT_CBTagg_31838.root";

    TFile *f = new TFile(filename);

    if (!f || f->IsZombie()) {
        std::cout << "Error: cannot open file " << filename << std::endl;
        return;
    }

    TTree *tagger = (TTree*)f->Get("tagger");

    if (!tagger) {
        std::cout << "Error: tagger tree not found." << std::endl;
        return;
    }

    Int_t nTagged;

    tagger->SetBranchAddress("nTagged", &nTagged);

    TH1D *hN = new TH1D(
        "hNTagged",
        "Tagger multiplicity per event; nTagged; events",
        500, 0, 500
    );

    Long64_t nEntries = tagger->GetEntries();

    for (Long64_t i = 0; i < nEntries; i++) {
        tagger->GetEntry(i);
        hN->Fill(nTagged);
    }

    TCanvas *c = new TCanvas("c", "nTagged distribution", 900, 600);
    hN->Draw();

    std::cout << "Entries = " << nEntries << std::endl;
    std::cout << "Mean nTagged = " << hN->GetMean() << std::endl;
    std::cout << "RMS nTagged = " << hN->GetRMS() << std::endl;

    std::cout << "Events with nTagged <= 10  = " << hN->Integral(1, hN->FindBin(10)) << std::endl;
    std::cout << "Events with nTagged <= 25  = " << hN->Integral(1, hN->FindBin(25)) << std::endl;
    std::cout << "Events with nTagged <= 50  = " << hN->Integral(1, hN->FindBin(50)) << std::endl;
    std::cout << "Events with nTagged <= 100 = " << hN->Integral(1, hN->FindBin(100)) << std::endl;
    std::cout << "Events with nTagged <= 150 = " << hN->Integral(1, hN->FindBin(150)) << std::endl;
    std::cout << "Events with nTagged <= 200 = " << hN->Integral(1, hN->FindBin(200)) << std::endl;

    c->SaveAs("nTagged_distribution.pdf");
}
