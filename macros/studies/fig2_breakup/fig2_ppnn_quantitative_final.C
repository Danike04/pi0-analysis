/*
 * fig2_ppnn_quantitative_final.C
 *
 * Purpose
 * -------
 * Quantitative ppnn-only test of the negative missing-energy tail in Fig. 2.
 * The ppnn normalization is determined once before the DeltaPhi cut in the
 * tail region and then reused after the 8/10/12 degree cuts.  Coherent MC is
 * independently normalized to the data core for each selection.
 *
 * Important interpretation
 * ------------------------
 * The ppnn scale is a shape normalization to the observed tail residual; it
 * is explicitly not an absolute cross-section prediction.
 *
 * Default regions: core [-10,15] MeV, tail [-60,-20] MeV,
 * comparison [-60,40] MeV.
 *
 * Usage
 * -----
 *   root -l -b -q 'fig2_ppnn_quantitative_final.C()'
 *
 * Fixed outputs include fig2_ppnn_quantitative.root/.txt and the before/after
 * PDF comparisons.  The code body below is unchanged from the supplied file.
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

namespace Fig2PPNNFinal {

struct Chi2Result {
    double chi2;
    int ndf;
    double pvalue;
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

// The ppnn normalization is determined ONLY from the distributions before
// the DeltaPhi cut. It is never refitted after the cut.
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

static void StyleTotal(TH1D* h)
{
    h->SetLineColor(kGreen+2);
    h->SetLineWidth(3);
    h->SetLineStyle(1);
    h->SetFillStyle(0);
}

} // namespace Fig2PPNNFinal

void fig2_ppnn_quantitative_final(
    const char* file8="../paper_fig2_missing_energy_cut8.root",
    const char* file10="../paper_fig2_missing_energy_cut10.root",
    const char* file12="../paper_fig2_missing_energy_cut12.root",
    const char* ppnnFile="fig2_incoherent_4bins_1M_wide.root",
    double coreMin=-10.0,
    double coreMax=15.0,
    double tailMin=-60.0,
    double tailMax=-20.0,
    double compareMin=-60.0,
    double compareMax=40.0)
{
    using namespace Fig2PPNNFinal;
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

    TFile* fOut = TFile::Open("fig2_ppnn_quantitative.root", "RECREATE");
    if (!fOut || fOut->IsZombie()) {
        std::cerr << "Cannot create output ROOT file" << std::endl;
        return;
    }

    std::ofstream txt("fig2_ppnn_quantitative.txt");
    txt << "# ppnn normalization fitted ONCE before the DeltaPhi cut in ["
        << tailMin << "," << tailMax << "] MeV and reused after every cut\n";
    txt << "# coherent MC is independently normalized in the core ["
        << coreMin << "," << coreMax << "] MeV before and after each cut\n";
    txt << "# ppnn normalization is shape-only, not an absolute cross-section prediction\n";
    txt << "# cut E_low E_high ppnn_scale Nppnn_before_global Nppnn_after_global survival_global "
        << "coh_scale_after chi2ndf_coh chi2ndf_total p_coh p_total "
        << "tail_residual_total tail_error_total tail_sigma_total\n";

    double ppnnScale[nE] = {0.0,0.0,0.0,0.0};
    double nPPNNBeforeGlobal[nE] = {0.0,0.0,0.0,0.0};

    // First and only ppnn fit: use the distributions BEFORE DeltaPhi cut.
    TCanvas* cBefore = new TCanvas("c_ppnn_before", "ppnn normalization before cut", 1300, 900);
    cBefore->Divide(2,2,0.002,0.002);

    for (int ie = 0; ie < nE; ++ie) {
        TH1D* d0 = dynamic_cast<TH1D*>(fin[0]->Get(Form("data_deltaE_before_phi_E%d",ie)));
        TH1D* c0 = dynamic_cast<TH1D*>(fin[0]->Get(Form("coh_mc_deltaE_before_phi_E%d",ie)));
        TH1D* p0 = dynamic_cast<TH1D*>(fPPNN->Get(Form("inc_deltaE_E%d_nodphi",ie)));

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
        if (!ppnn) {
            std::cerr << "Cannot map ppnn BEFORE histogram E" << ie << std::endl;
            fOut->Close();
            return;
        }

        nPPNNBeforeGlobal[ie] = p0->Integral(0,p0->GetNbinsX()+1);

        const double cohScale = CoreScale(data, coh, coreMin, coreMax);
        coh->Scale(cohScale);
        ppnnScale[ie] = FitPPNNBefore(data, coh, ppnn, tailMin, tailMax);
        ppnn->Scale(ppnnScale[ie]);

        TH1D* total = dynamic_cast<TH1D*>(coh->Clone(Form("total_before_E%d",ie)));
        total->SetDirectory(0);
        total->Add(ppnn);

        StyleData(data);
        StyleCoherent(coh);
        StylePPNN(ppnn);
        StyleTotal(total);

        cBefore->cd(ie+1);
        gPad->SetTicks(1,1);
        gPad->SetLeftMargin(0.13);
        gPad->SetBottomMargin(0.13);

        data->SetTitle(Form("%s, before #Delta#Phi cut", labels[ie]));
        data->GetXaxis()->SetTitle("#Delta E_{#pi^{0}}^{*} (MeV)");
        data->GetYaxis()->SetTitle("Counts");
        data->GetXaxis()->SetRangeUser(compareMin, compareMax);

        double ymax = std::max(data->GetMaximum(), total->GetMaximum());
        data->SetMinimum(0.0);
        data->SetMaximum(ymax > 0.0 ? 1.30*ymax : 1.0);

        data->Draw("E");
        coh->Draw("HIST SAME");
        ppnn->Draw("HIST SAME");
        total->Draw("HIST SAME");
        data->Draw("E SAME");

        TLegend* leg = new TLegend(0.49,0.65,0.89,0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.031);
        leg->AddEntry(data,"Data before cut","lep");
        leg->AddEntry(coh,"Coherent MC (core scaled)","l");
        leg->AddEntry(ppnn,"ppnn MC (tail fitted once)","l");
        leg->AddEntry(total,"Coherent + ppnn","l");
        leg->Draw();

        TLatex tx;
        tx.SetNDC();
        tx.SetTextSize(0.033);
        tx.DrawLatex(0.15,0.88,Form("ppnn scale = %.4g",ppnnScale[ie]));

        fOut->cd();
        data->Write();
        coh->Write();
        ppnn->Write();
        total->Write();
    }

    cBefore->SaveAs("fig2_ppnn_before_normalization.pdf");
    fOut->cd();
    cBefore->Write();

    // Prediction AFTER each DeltaPhi cut, with ppnnScale fixed above.
    TCanvas* cAfter = new TCanvas("c_ppnn_after", "ppnn prediction after cut", 1300, 900);
    cAfter->Print("fig2_ppnn_after_prediction.pdf[");

    for (int ic = 0; ic < nCuts; ++ic) {
        cAfter->Clear();
        cAfter->Divide(2,2,0.002,0.002);

        for (int ie = 0; ie < nE; ++ie) {
            TH1D* d0 = dynamic_cast<TH1D*>(fin[ic]->Get(Form("data_deltaE_after_phi_E%d",ie)));
            TH1D* c0 = dynamic_cast<TH1D*>(fin[ic]->Get(Form("coh_mc_deltaE_after_phi_E%d",ie)));
            TH1D* p0 = dynamic_cast<TH1D*>(fPPNN->Get(Form("inc_deltaE_E%d_dphi%d",ie,cuts[ic])));

            if (!d0 || !c0 || !p0) {
                std::cerr << "Missing AFTER histogram for cut " << cuts[ic]
                          << ", E index " << ie << std::endl;
                continue;
            }

            TH1D* data = dynamic_cast<TH1D*>(d0->Clone(Form("data_after_cut%d_E%d",cuts[ic],ie)));
            TH1D* coh = dynamic_cast<TH1D*>(c0->Clone(Form("coh_after_cut%d_E%d",cuts[ic],ie)));
            data->SetDirectory(0);
            coh->SetDirectory(0);

            TH1D* ppnn = MapToReference(p0, data, Form("ppnn_after_cut%d_E%d",cuts[ic],ie));
            if (!ppnn) continue;

            const double nAfterGlobal = p0->Integral(0,p0->GetNbinsX()+1);
            const double survivalGlobal = nPPNNBeforeGlobal[ie] > 0.0
                ? nAfterGlobal/nPPNNBeforeGlobal[ie] : 0.0;

            const double cohScale = CoreScale(data, coh, coreMin, coreMax);
            coh->Scale(cohScale);
            ppnn->Scale(ppnnScale[ie]); // fixed from BEFORE cut; no refit here

            TH1D* total = dynamic_cast<TH1D*>(coh->Clone(Form("total_after_cut%d_E%d",cuts[ic],ie)));
            total->SetDirectory(0);
            total->Add(ppnn);

            const Chi2Result chiCoh = EvaluateChi2(data, coh, compareMin, compareMax, 1);
            const Chi2Result chiTot = EvaluateChi2(data, total, compareMin, compareMax, 1);
            const double chi2ndfCoh = chiCoh.ndf > 0 ? chiCoh.chi2/chiCoh.ndf : 0.0;
            const double chi2ndfTot = chiTot.ndf > 0 ? chiTot.chi2/chiTot.ndf : 0.0;

            double tailResidual = 0.0;
            double tailError = 0.0;
            const double tailSigma = TailResidualSigma(
                data, total, tailMin, tailMax, tailResidual, tailError);

            txt << cuts[ic] << " " << eLow[ie] << " " << eHigh[ie] << " "
                << ppnnScale[ie] << " " << nPPNNBeforeGlobal[ie] << " "
                << nAfterGlobal << " " << survivalGlobal << " " << cohScale << " "
                << chi2ndfCoh << " " << chi2ndfTot << " "
                << chiCoh.pvalue << " " << chiTot.pvalue << " "
                << tailResidual << " " << tailError << " " << tailSigma << "\n";

            StyleData(data);
            StyleCoherent(coh);
            StylePPNN(ppnn);
            StyleTotal(total);

            cAfter->cd(ie+1);
            gPad->SetTicks(1,1);
            gPad->SetLeftMargin(0.13);
            gPad->SetBottomMargin(0.13);

            data->SetTitle(Form("%s, #Delta#Phi < %d^{#circ}",labels[ie],cuts[ic]));
            data->GetXaxis()->SetTitle("#Delta E_{#pi^{0}}^{*} (MeV)");
            data->GetYaxis()->SetTitle("Counts");
            data->GetXaxis()->SetRangeUser(compareMin, compareMax);

            double ymax = std::max(data->GetMaximum(), total->GetMaximum());
            data->SetMinimum(0.0);
            data->SetMaximum(ymax > 0.0 ? 1.30*ymax : 1.0);

            data->Draw("E");
            coh->Draw("HIST SAME");
            if (nAfterGlobal > 0.0) ppnn->Draw("HIST SAME");
            total->Draw("HIST SAME");
            data->Draw("E SAME");

            TLegend* leg = new TLegend(0.49,0.65,0.89,0.88);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.031);
            leg->AddEntry(data,"Data after cut","lep");
            leg->AddEntry(coh,"Coherent MC (core scaled)","l");
            leg->AddEntry(ppnn,"ppnn prediction (fixed scale)","l");
            leg->AddEntry(total,"Coherent + ppnn","l");
            leg->Draw();

            TLatex tx;
            tx.SetNDC();
            tx.SetTextSize(0.030);
            tx.DrawLatex(0.15,0.89,Form("raw ppnn: %.0f #rightarrow %.0f",nPPNNBeforeGlobal[ie],nAfterGlobal));
            tx.DrawLatex(0.15,0.84,Form("#chi^{2}/ndf: %.2f #rightarrow %.2f",chi2ndfCoh,chi2ndfTot));
            tx.DrawLatex(0.15,0.79,Form("tail residual after total: %.1f #sigma",tailSigma));

            fOut->cd();
            data->Write();
            coh->Write();
            ppnn->Write();
            total->Write();
        }

        cAfter->Print("fig2_ppnn_after_prediction.pdf");
    }

    cAfter->Print("fig2_ppnn_after_prediction.pdf]");
    fOut->cd();
    cAfter->Write();

    txt.close();
    fOut->Close();
    for (int ic = 0; ic < nCuts; ++ic) if (fin[ic]) fin[ic]->Close();
    fPPNN->Close();

    std::cout << "Saved:\n"
              << "  fig2_ppnn_before_normalization.pdf\n"
              << "  fig2_ppnn_after_prediction.pdf\n"
              << "  fig2_ppnn_quantitative.root\n"
              << "  fig2_ppnn_quantitative.txt\n";
}
