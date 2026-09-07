/*
 * fig2_breakup_three_channels_on_data_366.C
 *
 * Purpose
 * -------
 * Compare the shapes of the d+d, t+p and 3He+n breakup templates with the
 * 366-MeV (E3) experimental missing-energy distribution for the no-DeltaPhi,
 * 8, 10 and 12 degree selections.  This is a channel-identification / shape
 * study, not the main FULL-EMPTY refit macro.
 *
 * Inputs
 * ------
 *   paper_fig2_missing_energy_cut8.root
 *   paper_fig2_missing_energy_cut10.root
 *   paper_fig2_missing_energy_cut12.root
 *   fig2_dd_E3_366_wide.root
 *   fig2_tp_E3_366_wide.root
 *   fig2_he3n_all_available_bins_wide.root
 *
 * Outputs (fixed historical names)
 * --------------------------------
 *   fig2_breakup_three_channels_on_data_366.root
 *   fig2_breakup_three_channels_on_data_366.txt
 *   fig2_breakup_three_channels_on_data_366.pdf
 *
 * Usage
 * -----
 *   root -l -b -q 'fig2_breakup_three_channels_on_data_366.C()'
 *
 * The analysis body below is unchanged from the supplied source.
 */
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"

namespace Fig2ThreeChannelsOnData366 {

struct Chi2Result {
    double chi2;
    int ndf;
    double pvalue;
};

double IntegralRange(TH1* h, double xmin, double xmax, double* error = 0)
{
    if (error) *error = 0.0;
    if (!h) return 0.0;

    const int b1 = h->GetXaxis()->FindBin(xmin + 1.0e-9);
    const int b2 = h->GetXaxis()->FindBin(xmax - 1.0e-9);

    if (error)
        return h->IntegralAndError(b1, b2, *error);

    return h->Integral(b1, b2);
}

double IntegralGlobal(TH1* h)
{
    return h ? h->Integral(0, h->GetNbinsX() + 1) : 0.0;
}

TH1D* MapToReference(TH1* source, TH1D* reference, const char* name)
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
        out->SetBinError(
            bo,
            std::sqrt(e0 * e0 + e1 * e1));
    }

    return out;
}

double CoreScale(TH1D* data,
                 TH1D* coherent,
                 double coreMin,
                 double coreMax)
{
    const double d = IntegralRange(data, coreMin, coreMax);
    const double c = IntegralRange(coherent, coreMin, coreMax);

    return c > 0.0 ? d / c : 0.0;
}

double EqualTailScale(TH1D* data,
                      TH1D* coherentScaled,
                      TH1D* templateRaw,
                      double tailMin,
                      double tailMax)
{
    if (!data || !coherentScaled || !templateRaw)
        return 0.0;

    const double residual =
        IntegralRange(data, tailMin, tailMax) -
        IntegralRange(coherentScaled, tailMin, tailMax);

    const double templateArea =
        IntegralRange(templateRaw, tailMin, tailMax);

    if (residual <= 0.0 || templateArea <= 0.0)
        return 0.0;

    return residual / templateArea;
}

Chi2Result EvaluateChi2(TH1D* data,
                        TH1D* model,
                        double xmin,
                        double xmax,
                        int nFittedParameters)
{
    Chi2Result result = {0.0, 0, 0.0};

    if (!data || !model)
        return result;

    const int b1 =
        data->GetXaxis()->FindBin(xmin + 1.0e-9);
    const int b2 =
        data->GetXaxis()->FindBin(xmax - 1.0e-9);

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
    if (result.ndf < 0)
        result.ndf = 0;

    if (result.ndf > 0)
        result.pvalue = TMath::Prob(
            result.chi2,
            result.ndf);

    return result;
}

void StyleData(TH1D* h)
{
    h->SetLineColor(kBlack);
    h->SetMarkerColor(kBlack);
    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.50);
    h->SetLineWidth(2);
}

void StyleCoherent(TH1D* h)
{
    h->SetLineColor(kRed + 1);
    h->SetLineWidth(2);
    h->SetLineStyle(2);
    h->SetFillStyle(0);
}

void StyleTotalDD(TH1D* h)
{
    h->SetLineColor(kBlue + 1);
    h->SetLineWidth(3);
    h->SetLineStyle(1);
    h->SetFillStyle(0);
}

void StyleTotalTP(TH1D* h)
{
    h->SetLineColor(kMagenta + 2);
    h->SetLineWidth(3);
    h->SetLineStyle(7);
    h->SetFillStyle(0);
}

void StyleTotalHe3n(TH1D* h)
{
    h->SetLineColor(kGreen + 2);
    h->SetLineWidth(3);
    h->SetLineStyle(9);
    h->SetFillStyle(0);
}

} // namespace Fig2ThreeChannelsOnData366


void fig2_breakup_three_channels_on_data_366(
    const char* file8 =
        "../paper_fig2_missing_energy_cut8.root",
    const char* file10 =
        "../paper_fig2_missing_energy_cut10.root",
    const char* file12 =
        "../paper_fig2_missing_energy_cut12.root",
    const char* ddFile =
        "fig2_dd_E3_366_wide.root",
    const char* tpFile =
        "fig2_tp_E3_366_wide.root",
    const char* he3nFile =
        "fig2_he3n_all_available_bins_wide.root",
    double coreMin = -10.0,
    double coreMax = 15.0,
    double tailMin = -60.0,
    double tailMax = -20.0,
    double compareMin = -60.0,
    double compareMax = 40.0)
{
    using namespace Fig2ThreeChannelsOnData366;

    std::cout
        << "RUNNING three-channel comparison on data at 366 MeV"
        << std::endl;

    gStyle->SetOptStat(0);

    const int nPanels = 4;
    const int cuts[nPanels] = {-1, 8, 10, 12};

    const char* panelLabels[nPanels] = {
        "before #Delta#Phi cut",
        "#Delta#Phi < 8^{#circ}",
        "#Delta#Phi < 10^{#circ}",
        "#Delta#Phi < 12^{#circ}"
    };

    const char* dataFiles[nPanels] = {
        file8,
        file8,
        file10,
        file12
    };

    TFile* fData[nPanels] = {0, 0, 0, 0};

    for (int ip = 0; ip < nPanels; ++ip) {
        fData[ip] = TFile::Open(dataFiles[ip], "READ");

        if (!fData[ip] || fData[ip]->IsZombie()) {
            std::cerr
                << "Cannot open data file: "
                << dataFiles[ip]
                << std::endl;
            return;
        }
    }

    TFile* fDD = TFile::Open(ddFile, "READ");
    TFile* fTP = TFile::Open(tpFile, "READ");
    TFile* fHe3n = TFile::Open(he3nFile, "READ");

    if (!fDD || fDD->IsZombie()) {
        std::cerr
            << "Cannot open d+d file: "
            << ddFile
            << std::endl;
        return;
    }

    if (!fTP || fTP->IsZombie()) {
        std::cerr
            << "Cannot open t+p file: "
            << tpFile
            << std::endl;
        return;
    }

    if (!fHe3n || fHe3n->IsZombie()) {
        std::cerr
            << "Cannot open 3He+n file: "
            << he3nFile
            << std::endl;
        return;
    }

    // The arbitrary normalizations are determined once before the cut:
    // each breakup template is forced to have the same integral as
    // data - coherent MC in the negative tail [tailMin, tailMax].
    TH1D* dataBeforeRaw =
        dynamic_cast<TH1D*>(
            fData[0]->Get(
                "data_deltaE_before_phi_E3"));

    TH1D* cohBeforeRaw =
        dynamic_cast<TH1D*>(
            fData[0]->Get(
                "coh_mc_deltaE_before_phi_E3"));

    TH1* ddBeforeInput =
        dynamic_cast<TH1*>(
            fDD->Get(
                "inc_deltaE_E3_nodphi"));

    TH1* tpBeforeInput =
        dynamic_cast<TH1*>(
            fTP->Get(
                "inc_deltaE_E3_nodphi"));

    TH1* he3nBeforeInput =
        dynamic_cast<TH1*>(
            fHe3n->Get(
                "inc_deltaE_E3_nodphi"));

    if (!dataBeforeRaw ||
        !cohBeforeRaw ||
        !ddBeforeInput ||
        !tpBeforeInput ||
        !he3nBeforeInput) {

        std::cerr
            << "Missing one or more BEFORE histograms."
            << std::endl;
        return;
    }

    TH1D* dataBefore =
        dynamic_cast<TH1D*>(
            dataBeforeRaw->Clone(
                "data_before_for_scale"));

    TH1D* cohBefore =
        dynamic_cast<TH1D*>(
            cohBeforeRaw->Clone(
                "coh_before_for_scale"));

    dataBefore->SetDirectory(0);
    cohBefore->SetDirectory(0);

    TH1D* ddBefore =
        MapToReference(
            ddBeforeInput,
            dataBefore,
            "dd_before_for_scale");

    TH1D* tpBefore =
        MapToReference(
            tpBeforeInput,
            dataBefore,
            "tp_before_for_scale");

    TH1D* he3nBefore =
        MapToReference(
            he3nBeforeInput,
            dataBefore,
            "he3n_before_for_scale");

    if (!ddBefore || !tpBefore || !he3nBefore) {
        std::cerr
            << "Cannot map the BEFORE breakup templates."
            << std::endl;
        return;
    }

    const double coherentScaleBefore =
        CoreScale(
            dataBefore,
            cohBefore,
            coreMin,
            coreMax);

    cohBefore->Scale(
        coherentScaleBefore);

    const double targetTailArea =
        std::max(
            0.0,
            IntegralRange(
                dataBefore,
                tailMin,
                tailMax) -
            IntegralRange(
                cohBefore,
                tailMin,
                tailMax));

    const double ddScale =
        EqualTailScale(
            dataBefore,
            cohBefore,
            ddBefore,
            tailMin,
            tailMax);

    const double tpScale =
        EqualTailScale(
            dataBefore,
            cohBefore,
            tpBefore,
            tailMin,
            tailMax);

    const double he3nScale =
        EqualTailScale(
            dataBefore,
            cohBefore,
            he3nBefore,
            tailMin,
            tailMax);

    std::cout
        << "Arbitrary scales fixed BEFORE the cut:\n"
        << "  target tail area = "
        << targetTailArea << "\n"
        << "  d+d scale       = "
        << ddScale << "\n"
        << "  t+p scale       = "
        << tpScale << "\n"
        << "  3He+n scale     = "
        << he3nScale
        << std::endl;

    TFile* output =
        TFile::Open(
            "fig2_breakup_three_channels_on_data_366.root",
            "RECREATE");

    if (!output || output->IsZombie()) {
        std::cerr
            << "Cannot create output ROOT file."
            << std::endl;
        return;
    }

    std::ofstream txt(
        "fig2_breakup_three_channels_on_data_366.txt");

    txt
        << "# The three breakup templates are arbitrarily normalized ONCE "
           "before the DeltaPhi cut.\n"
        << "# Each template is scaled to the same data-minus-coherent "
           "integral in ["
        << tailMin << "," << tailMax << "] MeV.\n"
        << "# These three scales are then frozen after all DeltaPhi cuts.\n"
        << "# The coherent MC is normalized independently in the core ["
        << coreMin << "," << coreMax
        << "] MeV for every panel.\n"
        << "# target_tail_area "
        << targetTailArea << "\n"
        << "# dd_scale " << ddScale << "\n"
        << "# tp_scale " << tpScale << "\n"
        << "# he3n_scale " << he3nScale << "\n"
        << "# cut coh_scale raw_dd raw_tp raw_he3n "
           "chi2ndf_coh_dd chi2ndf_coh_tp chi2ndf_coh_he3n\n";

    TCanvas* canvas =
        new TCanvas(
            "c_breakup_three_channels_on_data_366",
            "Three breakup channels on data at 366 MeV",
            1350,
            920);

    canvas->Divide(
        2,
        2,
        0.002,
        0.002);

    for (int ip = 0; ip < nPanels; ++ip) {
        const bool before = cuts[ip] < 0;

        const TString dataName =
            before
            ? "data_deltaE_before_phi_E3"
            : "data_deltaE_after_phi_E3";

        const TString coherentName =
            before
            ? "coh_mc_deltaE_before_phi_E3"
            : "coh_mc_deltaE_after_phi_E3";

        const TString templateName =
            before
            ? "inc_deltaE_E3_nodphi"
            : Form(
                "inc_deltaE_E3_dphi%d",
                cuts[ip]);

        TH1D* dataInput =
            dynamic_cast<TH1D*>(
                fData[ip]->Get(dataName));

        TH1D* coherentInput =
            dynamic_cast<TH1D*>(
                fData[ip]->Get(coherentName));

        TH1* ddInput =
            dynamic_cast<TH1*>(
                fDD->Get(templateName));

        TH1* tpInput =
            dynamic_cast<TH1*>(
                fTP->Get(templateName));

        TH1* he3nInput =
            dynamic_cast<TH1*>(
                fHe3n->Get(templateName));

        if (!dataInput ||
            !coherentInput ||
            !ddInput ||
            !tpInput ||
            !he3nInput) {

            std::cerr
                << "Missing histogram for panel "
                << ip
                << " ("
                << panelLabels[ip]
                << ")."
                << std::endl;
            continue;
        }

        TH1D* data =
            dynamic_cast<TH1D*>(
                dataInput->Clone(
                    Form(
                        "data_panel%d",
                        ip)));

        TH1D* coherent =
            dynamic_cast<TH1D*>(
                coherentInput->Clone(
                    Form(
                        "coherent_panel%d",
                        ip)));

        data->SetDirectory(0);
        coherent->SetDirectory(0);

        TH1D* dd =
            MapToReference(
                ddInput,
                data,
                Form(
                    "dd_panel%d",
                    ip));

        TH1D* tp =
            MapToReference(
                tpInput,
                data,
                Form(
                    "tp_panel%d",
                    ip));

        TH1D* he3n =
            MapToReference(
                he3nInput,
                data,
                Form(
                    "he3n_panel%d",
                    ip));

        if (!dd || !tp || !he3n)
            continue;

        const double coherentScale =
            CoreScale(
                data,
                coherent,
                coreMin,
                coreMax);

        coherent->Scale(
            coherentScale);

        dd->Scale(
            ddScale);

        tp->Scale(
            tpScale);

        he3n->Scale(
            he3nScale);

        TH1D* totalDD =
            dynamic_cast<TH1D*>(
                coherent->Clone(
                    Form(
                        "coherent_plus_dd_panel%d",
                        ip)));

        TH1D* totalTP =
            dynamic_cast<TH1D*>(
                coherent->Clone(
                    Form(
                        "coherent_plus_tp_panel%d",
                        ip)));

        TH1D* totalHe3n =
            dynamic_cast<TH1D*>(
                coherent->Clone(
                    Form(
                        "coherent_plus_he3n_panel%d",
                        ip)));

        totalDD->SetDirectory(0);
        totalTP->SetDirectory(0);
        totalHe3n->SetDirectory(0);

        totalDD->Add(dd);
        totalTP->Add(tp);
        totalHe3n->Add(he3n);

        const Chi2Result chiDD =
            EvaluateChi2(
                data,
                totalDD,
                compareMin,
                compareMax,
                1);

        const Chi2Result chiTP =
            EvaluateChi2(
                data,
                totalTP,
                compareMin,
                compareMax,
                1);

        const Chi2Result chiHe3n =
            EvaluateChi2(
                data,
                totalHe3n,
                compareMin,
                compareMax,
                1);

        const double chi2ndfDD =
            chiDD.ndf > 0
            ? chiDD.chi2 / chiDD.ndf
            : 0.0;

        const double chi2ndfTP =
            chiTP.ndf > 0
            ? chiTP.chi2 / chiTP.ndf
            : 0.0;

        const double chi2ndfHe3n =
            chiHe3n.ndf > 0
            ? chiHe3n.chi2 / chiHe3n.ndf
            : 0.0;

        txt
            << cuts[ip] << " "
            << coherentScale << " "
            << IntegralGlobal(ddInput) << " "
            << IntegralGlobal(tpInput) << " "
            << IntegralGlobal(he3nInput) << " "
            << chi2ndfDD << " "
            << chi2ndfTP << " "
            << chi2ndfHe3n
            << "\n";

        StyleData(data);
        StyleCoherent(coherent);
        StyleTotalDD(totalDD);
        StyleTotalTP(totalTP);
        StyleTotalHe3n(totalHe3n);

        canvas->cd(
            ip + 1);

        gPad->SetTicks(
            1,
            1);

        gPad->SetLeftMargin(
            0.13);

        gPad->SetBottomMargin(
            0.13);

        data->SetTitle(
            Form(
                "366 MeV, %s",
                panelLabels[ip]));

        data->GetXaxis()->SetTitle(
            "#Delta E_{#pi^{0}}^{*} (MeV)");

        data->GetYaxis()->SetTitle(
            "Counts");

        data->GetXaxis()->SetRangeUser(
            compareMin,
            compareMax);

        double ymax =
            data->GetMaximum();

        ymax = std::max(
            ymax,
            totalDD->GetMaximum());

        ymax = std::max(
            ymax,
            totalTP->GetMaximum());

        ymax = std::max(
            ymax,
            totalHe3n->GetMaximum());

        data->SetMinimum(
            0.0);

        data->SetMaximum(
            ymax > 0.0
            ? 1.35 * ymax
            : 1.0);

        data->Draw(
            "E");

        coherent->Draw(
            "HIST SAME");

        totalDD->Draw(
            "HIST SAME");

        totalTP->Draw(
            "HIST SAME");

        totalHe3n->Draw(
            "HIST SAME");

        data->Draw(
            "E SAME");

        TLegend* legend =
            new TLegend(
                0.45,
                0.57,
                0.89,
                0.88);

        legend->SetBorderSize(
            0);

        legend->SetFillStyle(
            0);

        legend->SetTextSize(
            0.029);

        legend->AddEntry(
            data,
            "Data",
            "lep");

        legend->AddEntry(
            coherent,
            "Coherent MC (core scaled)",
            "l");

        legend->AddEntry(
            totalDD,
            "Coherent + d+d",
            "l");

        legend->AddEntry(
            totalTP,
            "Coherent + t+p",
            "l");

        legend->AddEntry(
            totalHe3n,
            "Coherent + ^{3}He+n",
            "l");

        legend->Draw();

        TLatex label;
        label.SetNDC();
        label.SetTextSize(
            0.028);

        label.DrawLatex(
            0.15,
            0.89,
            "equal tail area fixed before cut");

        label.DrawLatex(
            0.15,
            0.84,
            Form(
                "#chi^{2}/ndf: d+d %.2f, t+p %.2f, ^{3}He+n %.2f",
                chi2ndfDD,
                chi2ndfTP,
                chi2ndfHe3n));

        output->cd();

        data->Write();
        coherent->Write();
        dd->Write();
        tp->Write();
        he3n->Write();
        totalDD->Write();
        totalTP->Write();
        totalHe3n->Write();
    }

    canvas->SaveAs(
        "fig2_breakup_three_channels_on_data_366.pdf");

    output->cd();
    canvas->Write();

    txt.close();
    output->Close();

    for (int ip = 0; ip < nPanels; ++ip) {
        if (fData[ip])
            fData[ip]->Close();
    }

    fDD->Close();
    fTP->Close();
    fHe3n->Close();

    std::cout
        << "Saved:\n"
        << "  fig2_breakup_three_channels_on_data_366.pdf\n"
        << "  fig2_breakup_three_channels_on_data_366.root\n"
        << "  fig2_breakup_three_channels_on_data_366.txt\n";
}
