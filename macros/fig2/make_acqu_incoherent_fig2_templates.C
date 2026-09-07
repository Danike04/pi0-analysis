/*
 * make_acqu_incoherent_fig2_templates.C
 *
 * Purpose
 * -------
 * Convert an Acqu/Geant incoherent or breakup Monte Carlo ROOT file into the
 * standardized Fig. 2 template histograms produced by
 * `Fig2IncoherentFiller.h`.
 *
 * Input
 * -----
 * Default: Acqu_geant_ppnn.root. The file must contain `tracks`, `tagger`
 * and `setupParameters`. The macro reconstructs the best neutral-cluster
 * gamma-gamma pair, applies the pi0 invariant-mass selection, maps tagger
 * channels to photon energy, and fills the common Fig. 2 DeltaE/DeltaPhi
 * histograms.
 *
 * Main defaults
 * -------------
 * output: fig2_incoherent_templates.root
 * cluster threshold: 20 MeV
 * veto maximum: 0.10 MeV
 * m(gamma gamma): 110--155 MeV
 * tagger channel offset: configurable argument
 *
 * The output contains `inc_deltaE_E<bin>_dphi<cut>` for the 4 energy bins
 * and DeltaPhi cuts 8, 10 and 12 degrees, plus before-cut diagnostics and
 * processing counters.
 *
 * Usage example
 * -------------
 * root -l -b -q 'make_acqu_incoherent_fig2_templates.C("Acqu_geant_ppnn.root","fig2_incoherent_templates.root")'
 *
 * Dependencies
 * ------------
 * `Fig2Common.h` and `Fig2IncoherentFiller.h` must be in the same macro
 * directory (or otherwise visible to ROOT's include path).
 *
 * Original analysis code follows unchanged; only this header was added.
 */

#include "Fig2Common.h"
#include "Fig2IncoherentFiller.h"

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TParameter.h"
#include "TString.h"

#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <algorithm>

namespace {

class LeafArray {
public:
    LeafArray() : leaf(0) {}

    bool Bind(TTree* tree, const char* branchName, bool required = true)
    {
        TBranch* branch = tree ? tree->GetBranch(branchName) : 0;
        if (!branch) {
            if (required)
                std::cerr << "Missing branch '" << branchName
                          << "' in tree '" << (tree ? tree->GetName() : "null")
                          << "'" << std::endl;
            return !required;
        }

        leaf = tree->GetLeaf(branchName);
        if (!leaf && branch->GetListOfLeaves() &&
            branch->GetListOfLeaves()->GetEntries() > 0)
            leaf = (TLeaf*)branch->GetListOfLeaves()->At(0);

        if (!leaf) {
            if (required)
                std::cerr << "Cannot obtain leaf for branch '" << branchName
                          << "'" << std::endl;
            return !required;
        }

        std::cout << "  " << tree->GetName() << "/" << branchName
                  << " [" << leaf->GetTypeName() << "]" << std::endl;
        return true;
    }

    bool Valid() const { return leaf != 0; }

    int Size() const
    {
        if (!leaf) return 0;
        int n = leaf->GetNdata();
        if (n <= 0) n = leaf->GetLen();
        return n > 0 ? n : 1;
    }

    double At(int index) const
    {
        if (!leaf || index < 0 || index >= Size()) return 0.0;
        return leaf->GetValue(index);
    }

    double First() const { return Size() > 0 ? At(0) : 0.0; }

private:
    TLeaf* leaf;
};

double ToRadians(double value, int angleUnit)
{
    return angleUnit == 2 ? Fig2::DegToRad(value) : value;
}

} // namespace

/*
 * Macro specifico per l'output Acqu/AcquMC mostrato dal file:
 *
 *   tracks:
 *     nTracks, clusterEnergy, theta, phi, vetoEnergy
 *
 *   tagger:
 *     nTagged, taggedChannel
 *
 *   setupParameters:
 *     nTagger, TaggerPhotonEnergy
 *
 * I tree "tracks" e "tagger" sono letti con lo stesso entry index.
 * L'energia del fotone si ottiene da:
 *
 *   TaggerPhotonEnergy[taggedChannel + taggerChannelOffset]
 *
 * angleUnit:
 *   0 = automatico
 *   1 = radianti
 *   2 = gradi
 *
 * energyScale:
 *   1.0 se clusterEnergy e TaggerPhotonEnergy sono in MeV
 *   1000.0 se sono in GeV
 */
void make_acqu_incoherent_fig2_templates(
    const char* inputFile = "Acqu_geant_ppnn.root",
    const char* outputFile = "fig2_incoherent_templates.root",
    double energyScale = 1.0,
    int angleUnit = 0,
    int taggerChannelOffset = 0,
    double clusterThresholdMeV = 20.0,
    double vetoMaxMeV = 0.10,
    double mggMinMeV = 110.0,
    double mggMaxMeV = 155.0,
    Long64_t maxEvents = -1)
{
    TFile* input = TFile::Open(inputFile, "READ");
    if (!input || input->IsZombie()) {
        std::cerr << "Cannot open " << inputFile << std::endl;
        return;
    }

    TTree* tracks = dynamic_cast<TTree*>(input->Get("tracks"));
    TTree* tagger = dynamic_cast<TTree*>(input->Get("tagger"));
    TTree* setup  = dynamic_cast<TTree*>(input->Get("setupParameters"));

    if (!tracks || !tagger || !setup) {
        std::cerr << "Required trees not found. Need tracks, tagger and "
                     "setupParameters." << std::endl;
        input->Close();
        return;
    }

    std::cout << "Input: " << inputFile << std::endl;
    std::cout << "tracks entries = " << tracks->GetEntries() << std::endl;
    std::cout << "tagger entries = " << tagger->GetEntries() << std::endl;

    LeafArray nTracks, clusterEnergy, theta, phi, vetoEnergy;
    LeafArray nTagged, taggedChannel;
    LeafArray nTagger, taggerPhotonEnergy;

    bool ok = true;
    ok &= nTracks.Bind(tracks, "nTracks");
    ok &= clusterEnergy.Bind(tracks, "clusterEnergy");
    ok &= theta.Bind(tracks, "theta");
    ok &= phi.Bind(tracks, "phi");
    vetoEnergy.Bind(tracks, "vetoEnergy", false);

    ok &= nTagged.Bind(tagger, "nTagged");
    ok &= taggedChannel.Bind(tagger, "taggedChannel");

    ok &= nTagger.Bind(setup, "nTagger");
    ok &= taggerPhotonEnergy.Bind(setup, "TaggerPhotonEnergy");

    if (!ok) {
        input->Close();
        return;
    }

    setup->GetEntry(0);

    const int numberTaggerChannels = (int)nTagger.First();
    const int energyTableSize = taggerPhotonEnergy.Size();

    std::cout << "nTagger = " << numberTaggerChannels
              << ", TaggerPhotonEnergy entries = " << energyTableSize
              << std::endl;

    if (energyTableSize <= 0) {
        std::cerr << "Empty TaggerPhotonEnergy table." << std::endl;
        input->Close();
        return;
    }

    Fig2IncoherentFiller filler;

    Long64_t entries = std::min(tracks->GetEntries(), tagger->GetEntries());
    if (maxEvents >= 0 && maxEvents < entries) entries = maxEvents;

    Long64_t countTaggedEvents = 0;
    Long64_t countValidTaggerEnergy = 0;
    Long64_t countTwoNeutral = 0;
    Long64_t countBestPair = 0;
    Long64_t countMassPass = 0;
    Long64_t countFilledTags = 0;

    int resolvedAngleUnit = angleUnit;
    double minMappedEnergy = 1e99;
    double maxMappedEnergy = -1e99;
    int minChannelSeen = 999999;
    int maxChannelSeen = -999999;

    for (Long64_t entry = 0; entry < entries; ++entry) {
        tracks->GetEntry(entry);
        tagger->GetEntry(entry);

        const int nt = std::min((int)nTracks.First(),
                     std::min(clusterEnergy.Size(),
                     std::min(theta.Size(), phi.Size())));
        const int ng = std::min((int)nTagged.First(), taggedChannel.Size());

        if (ng <= 0) continue;
        ++countTaggedEvents;

        if (nt < 2) continue;

        if (resolvedAngleUnit == 0) {
            double maxAbsTheta = 0.0;
            double maxAbsPhi = 0.0;
            for (int i = 0; i < nt; ++i) {
                maxAbsTheta = std::max(maxAbsTheta, std::fabs(theta.At(i)));
                maxAbsPhi = std::max(maxAbsPhi, std::fabs(phi.At(i)));
            }
            resolvedAngleUnit =
                (maxAbsTheta > 3.5 || maxAbsPhi > 6.5) ? 2 : 1;
            std::cout << "Angle unit inferred as "
                      << (resolvedAngleUnit == 2 ? "degrees" : "radians")
                      << std::endl;
        }

        std::vector<int> neutral;
        for (int i = 0; i < nt; ++i) {
            const double energy = clusterEnergy.At(i)*energyScale;
            if (energy < clusterThresholdMeV) continue;

            if (vetoEnergy.Valid() && i < vetoEnergy.Size()) {
                const double veto = vetoEnergy.At(i)*energyScale;
                if (veto > vetoMaxMeV) continue;
            }
            neutral.push_back(i);
        }

        if (neutral.size() < 2) continue;
        ++countTwoNeutral;

        TLorentzVector bestGamma1;
        TLorentzVector bestGamma2;
        double bestDistance = 1e99;
        bool foundPair = false;

        for (size_t ia = 0; ia < neutral.size(); ++ia) {
            for (size_t ib = ia + 1; ib < neutral.size(); ++ib) {
                const int i = neutral[ia];
                const int j = neutral[ib];

                const double e1 = clusterEnergy.At(i)*energyScale;
                const double e2 = clusterEnergy.At(j)*energyScale;
                const double th1 = ToRadians(theta.At(i), resolvedAngleUnit);
                const double th2 = ToRadians(theta.At(j), resolvedAngleUnit);
                const double ph1 = ToRadians(phi.At(i), resolvedAngleUnit);
                const double ph2 = ToRadians(phi.At(j), resolvedAngleUnit);

                const TLorentzVector gamma1 =
                    Fig2::PhotonP4(e1, th1, ph1);
                const TLorentzVector gamma2 =
                    Fig2::PhotonP4(e2, th2, ph2);

                const double mass = (gamma1 + gamma2).M();
                const double distance =
                    std::fabs(mass - Fig2::kPi0MassMeV);

                if (distance < bestDistance) {
                    bestDistance = distance;
                    bestGamma1 = gamma1;
                    bestGamma2 = gamma2;
                    foundPair = true;
                }
            }
        }

        if (!foundPair) continue;
        ++countBestPair;

        const double bestMass = (bestGamma1 + bestGamma2).M();
        if (bestMass >= mggMinMeV && bestMass <= mggMaxMeV)
            ++countMassPass;

        for (int itag = 0; itag < ng; ++itag) {
            const int rawChannel = (int)taggedChannel.At(itag);
            const int lookupChannel = rawChannel + taggerChannelOffset;

            minChannelSeen = std::min(minChannelSeen, rawChannel);
            maxChannelSeen = std::max(maxChannelSeen, rawChannel);

            if (lookupChannel < 0 || lookupChannel >= energyTableSize)
                continue;

            const double eGamma =
                taggerPhotonEnergy.At(lookupChannel)*energyScale;

            minMappedEnergy = std::min(minMappedEnergy, eGamma);
            maxMappedEnergy = std::max(maxMappedEnergy, eGamma);
            ++countValidTaggerEnergy;

            if (Fig2::FindEnergyBin(eGamma) < 0) continue;

            filler.Fill(eGamma,
                        bestGamma1,
                        bestGamma2,
                        1.0,
                        true,
                        mggMinMeV,
                        mggMaxMeV);
            ++countFilledTags;
        }

        if ((entry + 1) % 100000 == 0)
            std::cout << "Processed " << entry + 1
                      << " / " << entries << std::endl;
    }

    TFile output(outputFile, "RECREATE");
    if (output.IsZombie()) {
        std::cerr << "Cannot create " << outputFile << std::endl;
        input->Close();
        return;
    }

    filler.Write(&output);

    TParameter<Long64_t>("n_entries", entries).Write();
    TParameter<Long64_t>("n_tagged_events", countTaggedEvents).Write();
    TParameter<Long64_t>("n_valid_tagger_energy", countValidTaggerEnergy).Write();
    TParameter<Long64_t>("n_two_neutral", countTwoNeutral).Write();
    TParameter<Long64_t>("n_best_pair", countBestPair).Write();
    TParameter<Long64_t>("n_mass_pass", countMassPass).Write();
    TParameter<Long64_t>("n_filled_tags", countFilledTags).Write();
    TParameter<int>("tagger_channel_offset", taggerChannelOffset).Write();
    TParameter<double>("cluster_threshold_MeV", clusterThresholdMeV).Write();
    TParameter<double>("mgg_min_MeV", mggMinMeV).Write();
    TParameter<double>("mgg_max_MeV", mggMaxMeV).Write();

    output.Write();
    output.Close();
    input->Close();

    std::cout << "\nWrote " << outputFile << std::endl;
    std::cout << "raw tagger channels seen: "
              << minChannelSeen << " ... " << maxChannelSeen << std::endl;

    if (countValidTaggerEnergy > 0)
        std::cout << "mapped photon energies: "
                  << minMappedEnergy << " ... "
                  << maxMappedEnergy << " MeV" << std::endl;

    std::cout << "tagged events     = " << countTaggedEvents << std::endl;
    std::cout << "valid tag energies= " << countValidTaggerEnergy << std::endl;
    std::cout << ">=2 neutral       = " << countTwoNeutral << std::endl;
    std::cout << "best pair         = " << countBestPair << std::endl;
    std::cout << "mgg pass          = " << countMassPass << std::endl;
    std::cout << "filled tags       = " << countFilledTags << std::endl;
}
