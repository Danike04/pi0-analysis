/*
 * check_ME_cut_visual.C
 *
 * PURPOSE
 *   Visual historical check of the chosen missing-energy integration window
 *   in a representative GoAT data run.
 *
 * WHAT THE MACRO ACTUALLY DOES
 *   - Opens GoAT_CBTagg_31832.root by default and selects tagger channel 51.
 *   - Uses Egamma = 212 MeV and requires 110 < m(pi0) < 145 MeV.
 *   - Performs prompt/random subtraction with prompt 760-790 ns and random
 *     450-700 ns.
 *   - Reconstructs missing energy and compares the fixed [-10,40] MeV window
 *     with an alternative peak +/-25 MeV integration interval.
 *
 * OUTPUT
 *   ME_cut_visual_ch59.pdf plus the corresponding integrals printed to stdout.
 *   Note that the historical output filename says "ch59" although the code
 *   selects CH=51; the inconsistent filename is intentionally preserved.
 *
 * USAGE
 *   root -l -b -q 'check_ME_cut_visual.C'
 *
 * REPOSITORY STATUS
 *   Historical cut-visualization study. It is retained because it documents
 *   how the missing-energy window was inspected, not as a production macro.
 *   The executable code below is unchanged from the raw file.
 */
void check_ME_cut_visual() {

    const char* filename = "GoAT_CBTagg_31832.root";

    const int CH = 51;

    // Cambia qui se vuoi provare altre energie
    const double Egamma = 212.0;

    const double Mtarget = 3727.38;

    // Timing attuale
    const double promptL = 760.0;
    const double promptR = 790.0;

    const double randomL = 450.0;
    const double randomR = 700.0;

    const double randomWeight = -(promptR - promptL)/(randomR - randomL);

    // Pi0 mass cut
    const double mPi0L = 110.0;
    const double mPi0R = 145.0;

    // Cut del PDF
    const double cutL = -10.0;
    const double cutR = 40.0;

    TFile *f = new TFile(filename);

    TTree *tagger = (TTree*)f->Get("tagger");
    TTree *pi0    = (TTree*)f->Get("neutralPions");

    const int MAXTAG = 5000;
    const int MAXPI0 = 200;

    Int_t nTagged;
    Int_t taggedChannel[MAXTAG];
    Double_t taggedTime[MAXTAG];

    Int_t nParticles;
    Double_t clusterEnergy[MAXPI0];
    Double_t theta[MAXPI0];
    Double_t phi[MAXPI0];
    Double_t mass[MAXPI0];
    Double_t time[MAXPI0];

    tagger->SetBranchAddress("nTagged", &nTagged);
    tagger->SetBranchAddress("taggedChannel", taggedChannel);
    tagger->SetBranchAddress("taggedTime", taggedTime);

    pi0->SetBranchAddress("nParticles", &nParticles);
    pi0->SetBranchAddress("clusterEnergy", clusterEnergy);
    pi0->SetBranchAddress("theta", theta);
    pi0->SetBranchAddress("phi", phi);
    pi0->SetBranchAddress("mass", mass);
    pi0->SetBranchAddress("time", time);

    TH1D *hME_prompt = new TH1D(
        "hME_prompt",
        "Missing Energy PROMPT only;E_{miss} (MeV);counts",
        300, -150, 150
    );

    TH1D *hME_sub = new TH1D(
        "hME_sub",
        "Missing Energy PROMPT - RANDOM;E_{miss} (MeV);counts",
        300, -150, 150
    );

    TH1D *hDT = new TH1D(
        "hDT",
        "Timing; t_{tagger}-t_{#pi^{0}};counts",
        1000, 300, 1100
    );

    Long64_t nEntries = tagger->GetEntries();

    for (Long64_t iev = 0; iev < nEntries; iev++) {

        tagger->GetEntry(iev);
        pi0->GetEntry(iev);

        for (int it = 0; it < nTagged; it++) {

            if (taggedChannel[it] != CH) continue;

            for (int ip = 0; ip < nParticles; ip++) {

                double m = mass[ip];

                if (m < mPi0L || m > mPi0R) continue;

                double dt = taggedTime[it] - time[ip];

                hDT->Fill(dt);

                double weight = 0.0;
                bool isPrompt = false;

                if (dt >= promptL && dt <= promptR) {
                    weight = +1.0;
                    isPrompt = true;
                }
                else if (dt >= randomL && dt <= randomR) {
                    weight = randomWeight;
                }
                else {
                    continue;
                }

                // Corretto: clusterEnergy trattata come p_pi0
                double p = clusterEnergy[ip];
                double Epi = sqrt(p*p + m*m);

                double th = theta[ip] * TMath::DegToRad();
                double ph = phi[ip] * TMath::DegToRad();

                double px = p * sin(th) * cos(ph);
                double py = p * sin(th) * sin(ph);
                double pz = p * cos(th);

                double pxmiss = -px;
                double pymiss = -py;
                double pzmiss = Egamma - pz;

                double pmiss2 = pxmiss*pxmiss + pymiss*pymiss + pzmiss*pzmiss;

                double Erecoil = sqrt(Mtarget*Mtarget + pmiss2);

                double missingEnergy = Egamma + Mtarget - Epi - Erecoil;

                if (isPrompt) {
                    hME_prompt->Fill(missingEnergy);
                }

                hME_sub->Fill(missingEnergy, weight);
            }
        }
    }

    double peakPrompt = hME_prompt->GetBinCenter(hME_prompt->GetMaximumBin());
    double peakSub = hME_sub->GetBinCenter(hME_sub->GetMaximumBin());

    double yPromptCut = hME_prompt->Integral(
        hME_prompt->FindBin(cutL),
        hME_prompt->FindBin(cutR)
    );

    double ySubCut = hME_sub->Integral(
        hME_sub->FindBin(cutL),
        hME_sub->FindBin(cutR)
    );

    double yPromptAroundPeak = hME_prompt->Integral(
        hME_prompt->FindBin(peakPrompt - 25),
        hME_prompt->FindBin(peakPrompt + 25)
    );

    double ySubAroundPeak = hME_sub->Integral(
        hME_sub->FindBin(peakSub - 25),
        hME_sub->FindBin(peakSub + 25)
    );

    std::cout << "==================================" << std::endl;
    std::cout << "Channel = " << CH << std::endl;
    std::cout << "Egamma = " << Egamma << " MeV" << std::endl;
    std::cout << "Prompt window = " << promptL << " to " << promptR << std::endl;
    std::cout << "Random window = " << randomL << " to " << randomR << std::endl;
    std::cout << "Random weight = " << randomWeight << std::endl;
    std::cout << "Peak prompt = " << peakPrompt << " MeV" << std::endl;
    std::cout << "Peak prompt-random = " << peakSub << " MeV" << std::endl;
    std::cout << "Prompt integral [-10,40] = " << yPromptCut << std::endl;
    std::cout << "Prompt integral peak +/-25 = " << yPromptAroundPeak << std::endl;
    std::cout << "Prompt-random integral [-10,40] = " << ySubCut << std::endl;
    std::cout << "Prompt-random integral peak +/-25 = " << ySubAroundPeak << std::endl;
    std::cout << "==================================" << std::endl;

    TCanvas *c = new TCanvas("c", "Visual ME cut check", 1200, 900);
    c->Divide(1,3);

    c->cd(1);
    hDT->Draw();

    c->cd(2);
    hME_prompt->Draw();

    double ymaxP = hME_prompt->GetMaximum();

    TLine *pCutL = new TLine(cutL, 0, cutL, ymaxP);
    TLine *pCutR = new TLine(cutR, 0, cutR, ymaxP);

    pCutL->SetLineColor(kRed);
    pCutR->SetLineColor(kRed);
    pCutL->SetLineWidth(3);
    pCutR->SetLineWidth(3);

    pCutL->Draw("same");
    pCutR->Draw("same");

    TLine *pPeakL = new TLine(peakPrompt - 25, 0, peakPrompt - 25, ymaxP);
    TLine *pPeakR = new TLine(peakPrompt + 25, 0, peakPrompt + 25, ymaxP);

    pPeakL->SetLineColor(kBlue);
    pPeakR->SetLineColor(kBlue);
    pPeakL->SetLineWidth(2);
    pPeakR->SetLineWidth(2);
    pPeakL->SetLineStyle(2);
    pPeakR->SetLineStyle(2);

    pPeakL->Draw("same");
    pPeakR->Draw("same");

    TLatex *txt1 = new TLatex();
    txt1->SetNDC();
    txt1->SetTextSize(0.04);
    txt1->DrawLatex(0.15, 0.85, "Red: PDF cut [-10,40]");
    txt1->DrawLatex(0.15, 0.78, "Blue dashed: peak #pm 25 MeV");

    c->cd(3);
    hME_sub->Draw();

    double ymaxS = hME_sub->GetMaximum();

    TLine *sCutL = new TLine(cutL, 0, cutL, ymaxS);
    TLine *sCutR = new TLine(cutR, 0, cutR, ymaxS);

    sCutL->SetLineColor(kRed);
    sCutR->SetLineColor(kRed);
    sCutL->SetLineWidth(3);
    sCutR->SetLineWidth(3);

    sCutL->Draw("same");
    sCutR->Draw("same");

    TLine *sPeakL = new TLine(peakSub - 25, 0, peakSub - 25, ymaxS);
    TLine *sPeakR = new TLine(peakSub + 25, 0, peakSub + 25, ymaxS);

    sPeakL->SetLineColor(kBlue);
    sPeakR->SetLineColor(kBlue);
    sPeakL->SetLineWidth(2);
    sPeakR->SetLineWidth(2);
    sPeakL->SetLineStyle(2);
    sPeakR->SetLineStyle(2);

    sPeakL->Draw("same");
    sPeakR->Draw("same");

    TLatex *txt2 = new TLatex();
    txt2->SetNDC();
    txt2->SetTextSize(0.04);
    txt2->DrawLatex(0.15, 0.85, "Red: PDF cut [-10,40]");
    txt2->DrawLatex(0.15, 0.78, "Blue dashed: peak #pm 25 MeV");

    c->SaveAs("ME_cut_visual_ch59.pdf");
}
