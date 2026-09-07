/*
 * fig2_ppnn_dd_E3_quantitative_final.C
 *
 * Purpose
 * -------
 * Quantitative Fig. 2 test of ppnn plus d+d at the E3 / 366-MeV point, with
 * the same before-cut normalization / after-cut prediction strategy used by
 * the related ppnn breakup studies.
 *
 * Default inputs include the three historical data/coherent cut files,
 * fig2_incoherent_4bins_1M_wide.root and fig2_dd_E3_366_wide.root.
 * Default regions: core [-10,15], tail [-60,-20], comparison [-60,40] MeV.
 *
 * Usage
 * -----
 *   root -l -b -q 'fig2_ppnn_dd_E3_quantitative_final.C()'
 *
 * The code body below is unchanged from the supplied source.
 */
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"

namespace Fig2BreakupDDFinal {

struct Chi2Result {
    double chi2;
    int ndf;
    double pvalue;
};

struct TwoComponentFit {
    double ppnnScale;
    double ddScale;
};

static double IntegralRange(TH1D* h, double xmin, double xmax, double* err = 0)
{
    if (err) *err = 0.0;
    if (!h) return 0.0;
    const int b1 = h->GetXaxis()->FindBin(xmin + 1.0e-9);
    const int b2 = h->GetXaxis()->FindBin(xmax - 1.0e-9);
    if (err) return h->IntegralAndError(b1, b2, *err);
    return h->Integral(b1, b2);
}

static TH1D* MapToReference(TH1D* source, TH1D* reference, const char* name)
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

static double CoreScale(TH1D* data, TH1D* coherent,
                        double coreMin, double coreMax)
{
    const double d = IntegralRange(data, coreMin, coreMax);
    const double c = IntegralRange(coherent, coreMin, coreMax);
    return c > 0.0 ? d/c : 0.0;
}

// Used for bins where only the ppnn template is available.
static double FitPPNNBefore(TH1D* dataBefore,
                            TH1D* coherentBeforeScaled,
                            TH1D* ppnnBeforeRaw,
                            double tailMin,
                            double tailMax)
{
    if (!dataBefore || !coherentBeforeScaled || !ppnnBeforeRaw) return 0.0;

    const int b1 = dataBefore->GetXaxis()->FindBin(tailMin + 1.0e-9);
    const int b2 = dataBefore->GetXaxis()->FindBin(tailMax - 1.0e-9);

    double numerator = 0.0;
    double denominator = 0.0;

    for (int b = b1; b <= b2; ++b) {
        const double d = dataBefore->GetBinContent(b);
        const double c = coherentBeforeScaled->GetBinContent(b);
        const double p = ppnnBeforeRaw->GetBinContent(b);
        if (p <= 0.0) continue;

        const double ed = dataBefore->GetBinError(b);
        const double ec = coherentBeforeScaled->GetBinError(b);
        const double variance = ed*ed + ec*ec;
        if (variance <= 0.0) continue;

        const double weight = 1.0/variance;
        numerator += weight*p*(d-c);
        denominator += weight*p*p;
    }

    if (denominator <= 0.0) return 0.0;
    const double scale = numerator/denominator;
    return scale > 0.0 ? scale : 0.0;
}

static double TwoComponentSSE(TH1D* dataBefore,
                              TH1D* coherentBeforeScaled,
                              TH1D* ppnnBeforeRaw,
                              TH1D* ddBeforeRaw,
                              double ppnnScale,
                              double ddScale,
                              double tailMin,
                              double tailMax)
{
    const int b1 = dataBefore->GetXaxis()->FindBin(tailMin + 1.0e-9);
    const int b2 = dataBefore->GetXaxis()->FindBin(tailMax - 1.0e-9);
    double sse = 0.0;

    for (int b = b1; b <= b2; ++b) {
        const double d = dataBefore->GetBinContent(b);
        const double c = coherentBeforeScaled->GetBinContent(b);
        const double p = ppnnBeforeRaw->GetBinContent(b);
        const double h = ddBeforeRaw->GetBinContent(b);
        const double ed = dataBefore->GetBinError(b);
        const double ec = coherentBeforeScaled->GetBinError(b);
        const double variance = ed*ed + ec*ec;
        if (variance <= 0.0) continue;

        const double diff = (d-c) - ppnnScale*p - ddScale*h;
        sse += diff*diff/variance;
    }
    return sse;
}

// Simultaneous non-negative weighted least-squares fit of ppnn and d+d.
// It is performed ONCE before the DeltaPhi cut for E3. The two scales are
// then frozen and propagated after every opening-angle cut.
static TwoComponentFit FitPPNNDDBefore(TH1D* dataBefore,
                                         TH1D* coherentBeforeScaled,
                                         TH1D* ppnnBeforeRaw,
                                         TH1D* ddBeforeRaw,
                                         double tailMin,
                                         double tailMax)
{
    TwoComponentFit result = {0.0, 0.0};
    if (!dataBefore || !coherentBeforeScaled ||
        !ppnnBeforeRaw || !ddBeforeRaw) return result;

    const int b1 = dataBefore->GetXaxis()->FindBin(tailMin + 1.0e-9);
    const int b2 = dataBefore->GetXaxis()->FindBin(tailMax - 1.0e-9);

    double a11 = 0.0;
    double a12 = 0.0;
    double a22 = 0.0;
    double rhs1 = 0.0;
    double rhs2 = 0.0;

    for (int b = b1; b <= b2; ++b) {
        const double d = dataBefore->GetBinContent(b);
        const double c = coherentBeforeScaled->GetBinContent(b);
        const double p = ppnnBeforeRaw->GetBinContent(b);
        const double h = ddBeforeRaw->GetBinContent(b);
        const double ed = dataBefore->GetBinError(b);
        const double ec = coherentBeforeScaled->GetBinError(b);
        const double variance = ed*ed + ec*ec;
        if (variance <= 0.0) continue;

        const double w = 1.0/variance;
        const double y = d-c;
        a11 += w*p*p;
        a12 += w*p*h;
        a22 += w*h*h;
        rhs1 += w*p*y;
        rhs2 += w*h*y;
    }

    // Candidate 0: both components fixed to zero.
    double bestP = 0.0;
    double bestH = 0.0;
    double bestSSE = TwoComponentSSE(dataBefore, coherentBeforeScaled,
                                     ppnnBeforeRaw, ddBeforeRaw,
                                     bestP, bestH, tailMin, tailMax);

    // Candidate 1: ppnn only.
    if (a11 > 0.0) {
        const double p = std::max(0.0, rhs1/a11);
        const double sse = TwoComponentSSE(dataBefore, coherentBeforeScaled,
                                           ppnnBeforeRaw, ddBeforeRaw,
                                           p, 0.0, tailMin, tailMax);
        if (sse < bestSSE) {
            bestSSE = sse;
            bestP = p;
            bestH = 0.0;
        }
    }

    // Candidate 2: d+d only.
    if (a22 > 0.0) {
        const double h = std::max(0.0, rhs2/a22);
        const double sse = TwoComponentSSE(dataBefore, coherentBeforeScaled,
                                           ppnnBeforeRaw, ddBeforeRaw,
                                           0.0, h, tailMin, tailMax);
        if (sse < bestSSE) {
            bestSSE = sse;
            bestP = 0.0;
            bestH = h;
        }
    }

    // Candidate 3: unconstrained two-component solution, accepted only if
    // both fitted scales are non-negative.
    const double det = a11*a22 - a12*a12;
    if (det > 0.0) {
        const double p = (rhs1*a22 - rhs2*a12)/det;
        const double h = (rhs2*a11 - rhs1*a12)/det;
        if (p >= 0.0 && h >= 0.0) {
            const double sse = TwoComponentSSE(dataBefore, coherentBeforeScaled,
                                               ppnnBeforeRaw, ddBeforeRaw,
                                               p, h, tailMin, tailMax);
            if (sse < bestSSE) {
                bestSSE = sse;
                bestP = p;
                bestH = h;
            }
        }
    }

    result.ppnnScale = bestP;
    result.ddScale = bestH;
    return result;
}

static Chi2Result EvaluateChi2(TH1D* data, TH1D* model,
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
        const double variance = ed*ed + em*em;
        if (variance <= 0.0) continue;

        const double diff = d-m;
        result.chi2 += diff*diff/variance;
        ++used;
    }

    result.ndf = used - nFittedParameters;
    if (result.ndf < 0) result.ndf = 0;
    if (result.ndf > 0) result.pvalue = TMath::Prob(result.chi2, result.ndf);
    return result;
}

static double TailResidualSigma(TH1D* data, TH1D* model,
                                double tailMin, double tailMax,
                                double& residual,
                                double& residualError)
{
    double ed = 0.0;
    double em = 0.0;
    const double d = IntegralRange(data, tailMin, tailMax, &ed);
    const double m = IntegralRange(model, tailMin, tailMax, &em);
    residual = d-m;
    residualError = std::sqrt(ed*ed + em*em);
    return residualError > 0.0 ? residual/residualError : 0.0;
}

static void StyleData(TH1D* h)
{
    h->SetLineColor(kBlack);
    h->SetMarkerColor(kBlack);
    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.45);
    h->SetLineWidth(2);
}

static void StyleCoherent(TH1D* h)
{
    h->SetLineColor(kRed+1);
    h->SetLineWidth(2);
    h->SetFillStyle(0);
}

static void StylePPNN(TH1D* h)
{
    h->SetLineColor(kBlue+1);
    h->SetLineWidth(2);
    h->SetLineStyle(2);
    h->SetFillStyle(0);
}

static void StyleDD(TH1D* h)
{
    h->SetLineColor(kMagenta+2);
    h->SetLineWidth(2);
    h->SetLineStyle(7);
    h->SetFillStyle(0);
}

static void StylePPNNTotal(TH1D* h)
{
    h->SetLineColor(kOrange+7);
    h->SetLineWidth(2);
    h->SetLineStyle(9);
    h->SetFillStyle(0);
}

static void StyleTotal(TH1D* h)
{
    h->SetLineColor(kGreen+2);
    h->SetLineWidth(3);
    h->SetLineStyle(1);
    h->SetFillStyle(0);
}

} // namespace Fig2BreakupDDFinal

void fig2_ppnn_dd_E3_quantitative_final(
    const char* file8="../paper_fig2_missing_energy_cut8.root",
    const char* file10="../paper_fig2_missing_energy_cut10.root",
    const char* file12="../paper_fig2_missing_energy_cut12.root",
    const char* ppnnFile="fig2_incoherent_4bins_1M_wide.root",
    const char* ddFile="fig2_dd_E3_366_wide.root",
    double coreMin=-10.0,
    double coreMax=15.0,
    double tailMin=-60.0,
    double tailMax=-20.0,
    double compareMin=-60.0,
    double compareMax=40.0)
{
    using namespace Fig2BreakupDDFinal;
    std::cout << "RUNNING ppnn + d+d E3 quantitative macro" << std::endl;
    gStyle->SetOptStat(0);

    const int nCuts = 3;
    const int cuts[nCuts] = {8,10,12};
    const char* files[nCuts] = {file8,file10,file12};

    const int nE = 4;
    const int eLow[nE] = {223,283,319,356};
    const int eHigh[nE] = {234,294,330,366};
    const char* labels[nE] = {"224 MeV","294 MeV","320 MeV","366 MeV"};

    TFile* fin[nCuts] = {0,0,0};
    for (int ic = 0; ic < nCuts; ++ic) {
        fin[ic] = TFile::Open(files[ic]);
        if (!fin[ic] || fin[ic]->IsZombie()) {
            std::cerr << "Cannot open " << files[ic] << std::endl;
            return;
        }
    }

    TFile* fPPNN = TFile::Open(ppnnFile);
    if (!fPPNN || fPPNN->IsZombie()) {
        std::cerr << "Cannot open " << ppnnFile << std::endl;
        return;
    }

    TFile* fDD = TFile::Open(ddFile);
    if (!fDD || fDD->IsZombie()) {
        std::cerr << "Cannot open " << ddFile << std::endl;
        return;
    }

    TFile* fOut = TFile::Open("fig2_ppnn_dd_E3_quantitative.root", "RECREATE");
    if (!fOut || fOut->IsZombie()) {
        std::cerr << "Cannot create output ROOT file" << std::endl;
        return;
    }

    std::ofstream txt("fig2_ppnn_dd_E3_quantitative.txt");
    txt << "# ppnn and d+d normalizations fitted ONCE before the DeltaPhi cut in ["
        << tailMin << "," << tailMax << "] MeV and reused after every cut\n";
    txt << "# Only E3 is expected to contain a non-empty d+d template; E0-E2 use ppnn only\n";
    txt << "# coherent MC is independently normalized in the core ["
        << coreMin << "," << coreMax << "] MeV before and after each cut\n";
    txt << "# The d+d input file is fig2_dd_E3_366_wide.root; availability is still detected from its non-zero E3 integral\n";
    txt << "# breakup normalizations are shape-only, not absolute cross-section predictions\n";
    txt << "# cut E_low E_high ppnn_scale dd_scale "
        << "Nppnn_before_global Nppnn_after_global ppnn_survival_global "
        << "Ndd_before_global Ndd_after_global dd_survival_global "
        << "coh_scale_after chi2ndf_coh chi2ndf_coh_ppnn chi2ndf_all "
        << "p_coh p_coh_ppnn p_all "
        << "tail_residual_coh_ppnn tail_error_coh_ppnn tail_sigma_coh_ppnn "
        << "tail_residual_all tail_error_all tail_sigma_all\n";

    double ppnnScale[nE] = {0.0,0.0,0.0,0.0};
    double ddScale[nE] = {0.0,0.0,0.0,0.0};
    double nPPNNBeforeGlobal[nE] = {0.0,0.0,0.0,0.0};
    double nDDBeforeGlobal[nE] = {0.0,0.0,0.0,0.0};

    // Fit breakup templates only BEFORE the DeltaPhi cut.
    TCanvas* cBefore = new TCanvas("c_breakup_before", "breakup normalization before cut", 1300, 900);
    cBefore->Divide(2,2,0.002,0.002);

    for (int ie = 0; ie < nE; ++ie) {
        TH1D* d0 = dynamic_cast<TH1D*>(fin[0]->Get(Form("data_deltaE_before_phi_E%d",ie)));
        TH1D* c0 = dynamic_cast<TH1D*>(fin[0]->Get(Form("coh_mc_deltaE_before_phi_E%d",ie)));
        TH1D* p0 = dynamic_cast<TH1D*>(fPPNN->Get(Form("inc_deltaE_E%d_nodphi",ie)));
        TH1D* h0 = dynamic_cast<TH1D*>(fDD->Get(Form("inc_deltaE_E%d_nodphi",ie)));
        const bool useDD = h0 && h0->Integral(0,h0->GetNbinsX()+1) > 0.0;

        if (!d0 || !c0 || !p0) {
            std::cerr << "Missing BEFORE histogram for E index " << ie << std::endl;
            fOut->Close();
            return;
        }

        TH1D* data = dynamic_cast<TH1D*>(d0->Clone(Form("data_before_E%d",ie)));
        TH1D* coh = dynamic_cast<TH1D*>(c0->Clone(Form("coh_before_E%d",ie)));
        data->SetDirectory(0);
        coh->SetDirectory(0);

        TH1D* ppnn = MapToReference(p0, data, Form("ppnn_before_E%d",ie));
        TH1D* dd = 0;
        if (useDD)
            dd = MapToReference(h0, data, Form("dd_before_E%d",ie));

        if (!ppnn || (useDD && !dd)) {
            std::cerr << "Cannot map breakup BEFORE histogram E" << ie << std::endl;
            fOut->Close();
            return;
        }

        nPPNNBeforeGlobal[ie] = p0->Integral(0,p0->GetNbinsX()+1);
        if (useDD)
            nDDBeforeGlobal[ie] = h0->Integral(0,h0->GetNbinsX()+1);

        const double cohScale = CoreScale(data, coh, coreMin, coreMax);
        coh->Scale(cohScale);

        if (useDD) {
            const TwoComponentFit fit = FitPPNNDDBefore(
                data, coh, ppnn, dd, tailMin, tailMax);
            ppnnScale[ie] = fit.ppnnScale;
            ddScale[ie] = fit.ddScale;
        } else {
            ppnnScale[ie] = FitPPNNBefore(data, coh, ppnn, tailMin, tailMax);
        }

        ppnn->Scale(ppnnScale[ie]);
        if (dd) dd->Scale(ddScale[ie]);

        TH1D* totalPPNN = dynamic_cast<TH1D*>(coh->Clone(Form("total_ppnn_before_E%d",ie)));
        totalPPNN->SetDirectory(0);
        totalPPNN->Add(ppnn);

        TH1D* totalAll = dynamic_cast<TH1D*>(totalPPNN->Clone(Form("total_all_before_E%d",ie)));
        totalAll->SetDirectory(0);
        if (dd) totalAll->Add(dd);

        StyleData(data);
        StyleCoherent(coh);
        StylePPNN(ppnn);
        if (dd) StyleDD(dd);
        StylePPNNTotal(totalPPNN);
        StyleTotal(totalAll);

        cBefore->cd(ie+1);
        gPad->SetTicks(1,1);
        gPad->SetLeftMargin(0.13);
        gPad->SetBottomMargin(0.13);

        data->SetTitle(Form("%s, before #Delta#Phi cut", labels[ie]));
        data->GetXaxis()->SetTitle("#Delta E_{#pi^{0}}^{*} (MeV)");
        data->GetYaxis()->SetTitle("Counts");
        data->GetXaxis()->SetRangeUser(compareMin, compareMax);

        double ymax = std::max(data->GetMaximum(), totalAll->GetMaximum());
        data->SetMinimum(0.0);
        data->SetMaximum(ymax > 0.0 ? 1.30*ymax : 1.0);

        data->Draw("E");
        coh->Draw("HIST SAME");
        ppnn->Draw("HIST SAME");
        if (dd) dd->Draw("HIST SAME");
        totalAll->Draw("HIST SAME");
        data->Draw("E SAME");

        TLegend* leg = new TLegend(0.47,0.61,0.89,0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.029);
        leg->AddEntry(data,"Data before cut","lep");
        leg->AddEntry(coh,"Coherent MC (core scaled)","l");
        leg->AddEntry(ppnn,"ppnn MC (fitted before cut)","l");
        if (dd) leg->AddEntry(dd,"d+d MC (fitted before cut)","l");
        leg->AddEntry(totalAll,dd ? "Coherent + ppnn + d+d" : "Coherent + ppnn","l");
        leg->Draw();

        TLatex tx;
        tx.SetNDC();
        tx.SetTextSize(0.032);
        tx.DrawLatex(0.15,0.88,Form("ppnn scale = %.4g",ppnnScale[ie]));
        if (dd)
            tx.DrawLatex(0.15,0.83,Form("d+d scale = %.4g",ddScale[ie]));

        fOut->cd();
        data->Write();
        coh->Write();
        ppnn->Write();
        if (dd) dd->Write();
        totalPPNN->Write();
        totalAll->Write();
    }

    cBefore->SaveAs("fig2_ppnn_dd_E3_before_normalization.pdf");
    fOut->cd();
    cBefore->Write();

    // Prediction AFTER each DeltaPhi cut with all breakup scales frozen.
    TCanvas* cAfter = new TCanvas("c_breakup_after", "breakup prediction after cut", 1300, 900);
    cAfter->Print("fig2_ppnn_dd_E3_after_prediction.pdf[");

    for (int ic = 0; ic < nCuts; ++ic) {
        cAfter->Clear();
        cAfter->Divide(2,2,0.002,0.002);

        for (int ie = 0; ie < nE; ++ie) {
            TH1D* d0 = dynamic_cast<TH1D*>(fin[ic]->Get(Form("data_deltaE_after_phi_E%d",ie)));
            TH1D* c0 = dynamic_cast<TH1D*>(fin[ic]->Get(Form("coh_mc_deltaE_after_phi_E%d",ie)));
            TH1D* p0 = dynamic_cast<TH1D*>(fPPNN->Get(Form("inc_deltaE_E%d_dphi%d",ie,cuts[ic])));
            const bool useDD = nDDBeforeGlobal[ie] > 0.0;
            TH1D* h0 = useDD
                ? dynamic_cast<TH1D*>(fDD->Get(Form("inc_deltaE_E%d_dphi%d",ie,cuts[ic])))
                : 0;

            if (!d0 || !c0 || !p0 || (useDD && !h0)) {
                std::cerr << "Missing AFTER histogram for cut " << cuts[ic]
                          << ", E index " << ie << std::endl;
                continue;
            }

            TH1D* data = dynamic_cast<TH1D*>(d0->Clone(Form("data_after_cut%d_E%d",cuts[ic],ie)));
            TH1D* coh = dynamic_cast<TH1D*>(c0->Clone(Form("coh_after_cut%d_E%d",cuts[ic],ie)));
            data->SetDirectory(0);
            coh->SetDirectory(0);

            TH1D* ppnn = MapToReference(p0, data, Form("ppnn_after_cut%d_E%d",cuts[ic],ie));
            TH1D* dd = 0;
            if (useDD)
                dd = MapToReference(h0, data, Form("dd_after_cut%d_E%d",cuts[ic],ie));
            if (!ppnn || (useDD && !dd)) continue;

            const double nPPNNAfterGlobal = p0->Integral(0,p0->GetNbinsX()+1);
            const double ppnnSurvivalGlobal = nPPNNBeforeGlobal[ie] > 0.0
                ? nPPNNAfterGlobal/nPPNNBeforeGlobal[ie] : 0.0;

            double nDDAfterGlobal = 0.0;
            double ddSurvivalGlobal = 0.0;
            if (dd) {
                nDDAfterGlobal = h0->Integral(0,h0->GetNbinsX()+1);
                ddSurvivalGlobal = nDDBeforeGlobal[ie] > 0.0
                    ? nDDAfterGlobal/nDDBeforeGlobal[ie] : 0.0;
            }

            const double cohScale = CoreScale(data, coh, coreMin, coreMax);
            coh->Scale(cohScale);
            ppnn->Scale(ppnnScale[ie]);
            if (dd) dd->Scale(ddScale[ie]);

            TH1D* totalPPNN = dynamic_cast<TH1D*>(coh->Clone(Form("total_ppnn_after_cut%d_E%d",cuts[ic],ie)));
            totalPPNN->SetDirectory(0);
            totalPPNN->Add(ppnn);

            TH1D* totalAll = dynamic_cast<TH1D*>(totalPPNN->Clone(Form("total_all_after_cut%d_E%d",cuts[ic],ie)));
            totalAll->SetDirectory(0);
            if (dd) totalAll->Add(dd);

            const Chi2Result chiCoh = EvaluateChi2(data, coh, compareMin, compareMax, 1);
            const Chi2Result chiPPNN = EvaluateChi2(data, totalPPNN, compareMin, compareMax, 1);
            const Chi2Result chiAll = EvaluateChi2(data, totalAll, compareMin, compareMax, 1);
            const double chi2ndfCoh = chiCoh.ndf > 0 ? chiCoh.chi2/chiCoh.ndf : 0.0;
            const double chi2ndfPPNN = chiPPNN.ndf > 0 ? chiPPNN.chi2/chiPPNN.ndf : 0.0;
            const double chi2ndfAll = chiAll.ndf > 0 ? chiAll.chi2/chiAll.ndf : 0.0;

            double tailResidualPPNN = 0.0;
            double tailErrorPPNN = 0.0;
            const double tailSigmaPPNN = TailResidualSigma(
                data, totalPPNN, tailMin, tailMax,
                tailResidualPPNN, tailErrorPPNN);

            double tailResidualAll = 0.0;
            double tailErrorAll = 0.0;
            const double tailSigmaAll = TailResidualSigma(
                data, totalAll, tailMin, tailMax,
                tailResidualAll, tailErrorAll);

            txt << cuts[ic] << " " << eLow[ie] << " " << eHigh[ie] << " "
                << ppnnScale[ie] << " " << ddScale[ie] << " "
                << nPPNNBeforeGlobal[ie] << " " << nPPNNAfterGlobal << " "
                << ppnnSurvivalGlobal << " "
                << nDDBeforeGlobal[ie] << " " << nDDAfterGlobal << " "
                << ddSurvivalGlobal << " " << cohScale << " "
                << chi2ndfCoh << " " << chi2ndfPPNN << " " << chi2ndfAll << " "
                << chiCoh.pvalue << " " << chiPPNN.pvalue << " " << chiAll.pvalue << " "
                << tailResidualPPNN << " " << tailErrorPPNN << " " << tailSigmaPPNN << " "
                << tailResidualAll << " " << tailErrorAll << " " << tailSigmaAll << "\n";

            StyleData(data);
            StyleCoherent(coh);
            StylePPNN(ppnn);
            if (dd) StyleDD(dd);
            StylePPNNTotal(totalPPNN);
            StyleTotal(totalAll);

            cAfter->cd(ie+1);
            gPad->SetTicks(1,1);
            gPad->SetLeftMargin(0.13);
            gPad->SetBottomMargin(0.13);

            data->SetTitle(Form("%s, #Delta#Phi < %d^{#circ}",labels[ie],cuts[ic]));
            data->GetXaxis()->SetTitle("#Delta E_{#pi^{0}}^{*} (MeV)");
            data->GetYaxis()->SetTitle("Counts");
            data->GetXaxis()->SetRangeUser(compareMin, compareMax);

            double ymax = std::max(data->GetMaximum(), totalAll->GetMaximum());
            data->SetMinimum(0.0);
            data->SetMaximum(ymax > 0.0 ? 1.32*ymax : 1.0);

            data->Draw("E");
            coh->Draw("HIST SAME");
            if (nPPNNAfterGlobal > 0.0) ppnn->Draw("HIST SAME");
            if (dd && nDDAfterGlobal > 0.0) dd->Draw("HIST SAME");
            if (dd) totalPPNN->Draw("HIST SAME");
            totalAll->Draw("HIST SAME");
            data->Draw("E SAME");

            TLegend* leg = new TLegend(0.45,0.57,0.89,0.88);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.027);
            leg->AddEntry(data,"Data after cut","lep");
            leg->AddEntry(coh,"Coherent MC (core scaled)","l");
            leg->AddEntry(ppnn,"ppnn prediction (fixed scale)","l");
            if (dd) {
                leg->AddEntry(dd,"d+d prediction (fixed scale)","l");
                leg->AddEntry(totalPPNN,"Coherent + ppnn","l");
                leg->AddEntry(totalAll,"Coherent + ppnn + d+d","l");
            } else {
                leg->AddEntry(totalAll,"Coherent + ppnn","l");
            }
            leg->Draw();

            TLatex tx;
            tx.SetNDC();
            tx.SetTextSize(0.029);
            if (dd) {
                tx.DrawLatex(0.15,0.89,Form("raw ppnn: %.0f #rightarrow %.0f; d+d: %.0f #rightarrow %.0f",
                                           nPPNNBeforeGlobal[ie],nPPNNAfterGlobal,
                                           nDDBeforeGlobal[ie],nDDAfterGlobal));
                tx.DrawLatex(0.15,0.84,Form("#chi^{2}/ndf: coh %.2f, +ppnn %.2f, +d+d %.2f",
                                           chi2ndfCoh,chi2ndfPPNN,chi2ndfAll));
                tx.DrawLatex(0.15,0.79,Form("tail residual: +ppnn %.1f#sigma, all %.1f#sigma",
                                           tailSigmaPPNN,tailSigmaAll));
            } else {
                tx.DrawLatex(0.15,0.89,Form("raw ppnn: %.0f #rightarrow %.0f",
                                           nPPNNBeforeGlobal[ie],nPPNNAfterGlobal));
                tx.DrawLatex(0.15,0.84,Form("#chi^{2}/ndf: %.2f #rightarrow %.2f",
                                           chi2ndfCoh,chi2ndfAll));
                tx.DrawLatex(0.15,0.79,Form("tail residual after total: %.1f#sigma",
                                           tailSigmaAll));
            }

            fOut->cd();
            data->Write();
            coh->Write();
            ppnn->Write();
            if (dd) dd->Write();
            totalPPNN->Write();
            totalAll->Write();
        }

        cAfter->Print("fig2_ppnn_dd_E3_after_prediction.pdf");
    }

    cAfter->Print("fig2_ppnn_dd_E3_after_prediction.pdf]");
    fOut->cd();
    cAfter->Write();

    txt.close();
    fOut->Close();
    for (int ic = 0; ic < nCuts; ++ic) if (fin[ic]) fin[ic]->Close();
    fPPNN->Close();
    fDD->Close();

    std::cout << "Saved:\n"
              << "  fig2_ppnn_dd_E3_before_normalization.pdf\n"
              << "  fig2_ppnn_dd_E3_after_prediction.pdf\n"
              << "  fig2_ppnn_dd_E3_quantitative.root\n"
              << "  fig2_ppnn_dd_E3_quantitative.txt\n";
}
