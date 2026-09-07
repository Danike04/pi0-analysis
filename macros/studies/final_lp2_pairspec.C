/*
 * final_lp2_pairspec.C
 *
 * Purpose
 * -------
 * Reproduce the 2026 run-by-run consistency study comparing the Ladder/P2
 * ratio with the Pair-Spectrometer response.
 *
 * This is NOT the general tagging-efficiency extraction macro. It is a
 * dedicated study for a fixed set of runs and was used to test whether the
 * Ladder/P2 observable tracks changes in tagging efficiency.
 *
 * Inputs
 * ------
 * 1. lp2_values.csv in the current working directory. Each readable row is
 *        run  LadderP2Ratio
 *    (commas are also accepted as separators).
 * 2. ARHist_CBTagg_<run>.root for each selected run, also expected in the
 *    current working directory.
 *
 * For each valid FPD channel the PairSpec quantity is
 *
 *     R_i = 2500 * (B_i - C_i) / A_i,
 *
 * with A=FPD scaler, B=gated PairSpec and C=delayed PairSpec. Channels with
 * A<=0 or non-finite/non-physical R_i (R_i<=0 or R_i>1) are rejected.
 * The macro sums R_i over channels, then studies
 *
 *     LadderP2Ratio * PairSpecSum
 *
 * versus run number.
 *
 * Historical run selection
 * ------------------------
 * The run list is hard-coded in the macro. Runs 32791--32796 are explicitly
 * excluded because they were quadrupole-current test runs; runs without an
 * entry in lp2_values.csv are skipped.
 *
 * Output
 * ------
 * final_lp2_pairspec_values.csv
 * LadderP2Ratio_PairSpec_consistency.pdf
 * LadderP2Ratio_PairSpec_consistency.png
 * LadderP2Ratio_PairSpec_consistency.root
 *
 * Usage
 * -----
 * root -l -b -q final_lp2_pairspec.C
 *
 * The code below is the original analysis code. This header only documents
 * its behaviour; the calculation, run list and outputs are unchanged.
 */

#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TClass.h>
#include <TROOT.h>
#include <TH1.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TString.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <cmath>


struct Result
{
    int run;

    double ladderP2Ratio;

    double pairSpecSum;
    double pairSpecError;

    double product;
    double productError;

    int usedChannels;
    int rejectedChannels;
};


TObject* FindObjectRecursive(
    TDirectory* directory,
    const char* objectName
)
{
    if (!directory)
        return nullptr;

    TObject* directObject =
        directory->Get(objectName);

    if (directObject)
        return directObject;

    TIter nextKey(directory->GetListOfKeys());

    while (TKey* key =
           dynamic_cast<TKey*>(nextKey())) {

        TClass* objectClass =
            gROOT->GetClass(key->GetClassName());

        if (!objectClass)
            continue;

        if (!objectClass->InheritsFrom(
                TDirectory::Class()
            )) {
            continue;
        }

        TDirectory* subDirectory =
            dynamic_cast<TDirectory*>(
                directory->Get(key->GetName())
            );

        if (!subDirectory)
            continue;

        TObject* result =
            FindObjectRecursive(
                subDirectory,
                objectName
            );

        if (result)
            return result;
    }

    return nullptr;
}


TH1* FindHistogram(
    TFile* file,
    const std::vector<std::string>& possibleNames
)
{
    if (!file)
        return nullptr;

    for (const std::string& name : possibleNames) {

        TObject* object =
            FindObjectRecursive(
                file,
                name.c_str()
            );

        TH1* histogram =
            dynamic_cast<TH1*>(object);

        if (histogram)
            return histogram;
    }

    return nullptr;
}


std::map<int, double> ReadLadderP2Ratios(
    const char* filename
)
{
    std::map<int, double> values;

    std::ifstream input(filename);

    if (!input.is_open()) {

        std::cerr
            << "Cannot open "
            << filename
            << std::endl;

        return values;
    }

    std::string line;

    while (std::getline(input, line)) {

        if (line.empty())
            continue;

        std::replace(
            line.begin(),
            line.end(),
            ',',
            ' '
        );

        std::stringstream stream(line);

        int run = 0;
        double value = 0.0;

        if (!(stream >> run >> value))
            continue;

        if (value > 0.0)
            values[run] = value;
    }

    return values;
}


bool IsQuadrupoleTestRun(int run)
{
    return run >= 32791 && run <= 32796;
}


void final_lp2_pairspec()
{
    /*
     * Selected interval.
     *
     * Run 32785 is kept in the list but will be skipped
     * automatically because its LadderP2Ratio is missing.
     *
     * Runs 32791-32796 are excluded because they were
     * quadrupole-current tests.
     */

    const int runArray[] = {
        32776,
        32780,
        32785,
        32786,
        32787,
        32788,
        32789,
        32790,
        32791,
        32792,
        32793,
        32794,
        32795,
        32796,
        32797,
        32798
    };

    const int numberOfRuns =
        sizeof(runArray) / sizeof(runArray[0]);

    /*
     * Read LadderP2Ratio values from the file that
     * you already filled using the elog.
     */

    const std::map<int, double> ladderP2Ratios =
        ReadLadderP2Ratios(
            "lp2_values.csv"
        );

    if (ladderP2Ratios.empty()) {

        std::cerr
            << "No LadderP2Ratio values found in "
            << "lp2_values.csv"
            << std::endl;

        return;
    }

    std::vector<Result> results;

    std::ofstream output(
        "final_lp2_pairspec_values.csv"
    );

    output
        << "run,"
        << "LadderP2Ratio,"
        << "PairSpecSum,"
        << "PairSpecError,"
        << "LadderP2Ratio_times_PairSpecSum,"
        << "productError,"
        << "usedChannels,"
        << "rejectedChannels\n";

    /*
     * Loop over selected runs.
     */

    for (
        int index = 0;
        index < numberOfRuns;
        ++index
    ) {
        const int run =
            runArray[index];

        if (IsQuadrupoleTestRun(run)) {

            std::cout
                << "Run " << run
                << " skipped: quadrupole-current test."
                << std::endl;

            continue;
        }

        std::map<int, double>::const_iterator lp2Iterator =
            ladderP2Ratios.find(run);

        if (lp2Iterator == ladderP2Ratios.end()) {

            std::cout
                << "Run " << run
                << " skipped: LadderP2Ratio unavailable."
                << std::endl;

            continue;
        }

        TString filename =
            Form(
                "ARHist_CBTagg_%d.root",
                run
            );

        TFile* file =
            TFile::Open(
                filename,
                "READ"
            );

        if (!file || file->IsZombie()) {

            std::cerr
                << "Cannot open "
                << filename
                << std::endl;

            if (file)
                delete file;

            continue;
        }

        TH1* hFPD =
            FindHistogram(
                file,
                {
                    "FPD_ScalerAcc",
                    "FPD_scalerAcc"
                }
            );

        TH1* hGated =
            FindHistogram(
                file,
                {
                    "PairSpec_SumGated"
                }
            );

        TH1* hDelayed =
            FindHistogram(
                file,
                {
                    "PairSpec_SumGatedDly"
                }
            );

        if (!hFPD || !hGated || !hDelayed) {

            std::cerr
                << "Missing histograms in run "
                << run
                << std::endl;

            file->Close();
            delete file;

            continue;
        }

        int numberOfBins =
            hFPD->GetNbinsX();

        numberOfBins =
            std::min(
                numberOfBins,
                hGated->GetNbinsX()
            );

        numberOfBins =
            std::min(
                numberOfBins,
                hDelayed->GetNbinsX()
            );

        double pairSpecSum = 0.0;
        double pairSpecVariance = 0.0;

        int usedChannels = 0;
        int rejectedChannels = 0;

        /*
         * Channel-by-channel PairSpec calculation:
         *
         * R_i = 2500 * (B_i - C_i) / A_i
         *
         * PairSpecSum = sum_i R_i
         */

        for (
            int bin = 1;
            bin <= numberOfBins;
            ++bin
        ) {
            const double A =
                hFPD->GetBinContent(bin);

            const double B =
                hGated->GetBinContent(bin);

            const double C =
                hDelayed->GetBinContent(bin);

            if (A <= 0.0) {

                ++rejectedChannels;
                continue;
            }

            const double R =
                2500.0
                * (B - C)
                / A;

            if (
                !std::isfinite(R)
                ||
                R <= 0.0
                ||
                R > 1.0
            ) {
                ++rejectedChannels;
                continue;
            }

            double RError = 0.0;

            if (B + C > 0.0) {

                RError =
                    2500.0
                    * std::sqrt(B + C)
                    / A;
            }

            pairSpecSum += R;

            pairSpecVariance +=
                RError * RError;

            ++usedChannels;
        }

        const double pairSpecError =
            std::sqrt(pairSpecVariance);

        const double ladderP2Ratio =
            lp2Iterator->second;

        const double product =
            ladderP2Ratio
            * pairSpecSum;

        /*
         * No uncertainty is available for LadderP2Ratio.
         * Only the PairSpec statistical uncertainty is used.
         */

        const double productError =
            ladderP2Ratio
            * pairSpecError;

        Result result;

        result.run =
            run;

        result.ladderP2Ratio =
            ladderP2Ratio;

        result.pairSpecSum =
            pairSpecSum;

        result.pairSpecError =
            pairSpecError;

        result.product =
            product;

        result.productError =
            productError;

        result.usedChannels =
            usedChannels;

        result.rejectedChannels =
            rejectedChannels;

        results.push_back(result);

        output
            << run << ","
            << std::setprecision(12)
            << ladderP2Ratio << ","
            << pairSpecSum << ","
            << pairSpecError << ","
            << product << ","
            << productError << ","
            << usedChannels << ","
            << rejectedChannels
            << "\n";

        std::cout
            << "Run " << run
            << " | LadderP2Ratio = "
            << ladderP2Ratio
            << " | PairSpec sum = "
            << pairSpecSum
            << " | product = "
            << product
            << std::endl;

        file->Close();
        delete file;
    }

    output.close();

    if (results.size() < 2) {

        std::cerr
            << "Not enough valid points."
            << std::endl;

        return;
    }

    /*
     * Calculate mean and run-to-run RMS.
     */

    double mean = 0.0;

    for (const Result& result : results)
        mean += result.product;

    mean /= results.size();

    double rms = 0.0;

    for (const Result& result : results) {

        rms +=
            std::pow(
                result.product - mean,
                2
            );
    }

    rms =
        std::sqrt(
            rms / results.size()
        );

    const double relativeRMS =
        100.0 * rms / mean;

    /*
     * Create clean graph.
     */

    TGraphErrors* graph =
        new TGraphErrors();

    graph->SetName(
        "gLadderP2RatioTimesPairSpec"
    );

    graph->SetTitle(
        "LadderP2Ratio #times PairSpec sum versus run number;"
        "Run number;"
        "LadderP2Ratio #times PairSpec sum"
    );

    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1.2);
    graph->SetLineWidth(1);

    double minimumY =
        results.front().product;

    double maximumY =
        results.front().product;

    for (
        unsigned int index = 0;
        index < results.size();
        ++index
    ) {
        const Result& result =
            results[index];

        graph->SetPoint(
            index,
            result.run,
            result.product
        );

        graph->SetPointError(
            index,
            0.0,
            result.productError
        );

        minimumY =
            std::min(
                minimumY,
                result.product
            );

        maximumY =
            std::max(
                maximumY,
                result.product
            );
    }

    const double yMargin =
        0.18 * (maximumY - minimumY);

    graph->SetMinimum(
        minimumY - yMargin
    );

    graph->SetMaximum(
        maximumY + yMargin
    );

    gStyle->SetOptStat(0);
    gStyle->SetTitleBorderSize(0);

    TCanvas* canvas =
        new TCanvas(
            "canvas",
            "LadderP2Ratio PairSpec consistency",
            1200,
            750
        );

    canvas->SetLeftMargin(0.13);
    canvas->SetRightMargin(0.05);
    canvas->SetBottomMargin(0.13);
    canvas->SetTopMargin(0.10);

    canvas->SetGridx();
    canvas->SetGridy();

    graph->Draw("AP");

    graph->GetXaxis()->SetNoExponent(true);
    graph->GetXaxis()->SetMaxDigits(6);
    graph->GetXaxis()->SetNdivisions(510);

    graph->GetXaxis()->SetTitleSize(0.045);
    graph->GetYaxis()->SetTitleSize(0.045);

    graph->GetXaxis()->SetLabelSize(0.036);
    graph->GetYaxis()->SetLabelSize(0.036);

    graph->GetXaxis()->SetTitleOffset(1.15);
    graph->GetYaxis()->SetTitleOffset(1.25);

    /*
     * Horizontal line at the arithmetic mean.
     */

    TLine* meanLine =
        new TLine(
            results.front().run - 0.5,
            mean,
            results.back().run + 0.5,
            mean
        );

    meanLine->SetLineStyle(2);
    meanLine->SetLineWidth(2);
    meanLine->Draw("SAME");

    /*
     * Only the useful information is shown.
     * No constant fit and no misleading chi-square.
     */

    TLatex text;

    text.SetNDC(true);
    text.SetTextSize(0.034);

    text.DrawLatex(
        0.16,
        0.84,
        Form(
            "Mean = %.2f",
            mean
        )
    );

    text.DrawLatex(
        0.16,
        0.79,
        Form(
            "Run-to-run RMS / mean = %.2f%%",
            relativeRMS
        )
    );

    text.SetTextSize(0.027);

    text.DrawLatex(
        0.16,
        0.735,
        "Runs 32791-32796 excluded: quadrupole-current tests"
    );

    canvas->SaveAs(
        "LadderP2Ratio_PairSpec_consistency.pdf"
    );

    canvas->SaveAs(
        "LadderP2Ratio_PairSpec_consistency.png"
    );

    TFile outputROOT(
        "LadderP2Ratio_PairSpec_consistency.root",
        "RECREATE"
    );

    graph->Write();
    meanLine->Write("meanLine");

    outputROOT.Close();

    std::cout
        << "\nFinal result:"
        << "\n  Number of normal runs = "
        << results.size()
        << "\n  Mean product = "
        << mean
        << "\n  RMS = "
        << rms
        << "\n  Relative RMS = "
        << relativeRMS
        << " %"
        << "\n\nCreated:"
        << "\n  LadderP2Ratio_PairSpec_consistency.pdf"
        << "\n  LadderP2Ratio_PairSpec_consistency.png"
        << "\n  LadderP2Ratio_PairSpec_consistency.root"
        << "\n  final_lp2_pairspec_values.csv"
        << std::endl;
}
