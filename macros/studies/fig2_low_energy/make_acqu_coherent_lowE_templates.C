/*
 * make_acqu_coherent_lowE_templates.C
 *
 * Purpose
 * -------
 * Convert a coherent Acqu/Geant pi0 MC sample into the low-energy Fig. 2
 * template format defined by LowECoherentCommon.h / LowECoherentFiller.h.
 * Energy bins are 135-155, 155-175, 175-195 and 195-215 MeV; the filler writes
 * m(gamma gamma), no-DeltaPhi and 8/10/12-degree missing-energy histograms.
 *
 * Historical default input: Acqu_geant_He4pi0_coh_135_155.root.
 * Default reconstruction cuts: cluster threshold 20 MeV, veto <= 0.10 MeV,
 * 110 < m(gamma gamma) < 155 MeV.
 *
 * Usage
 * -----
 *   root -l -b -q 'make_acqu_coherent_lowE_templates.C()'
 *
 * Keep LowECoherentCommon.h and LowECoherentFiller.h in the same directory.
 * The code body below is unchanged from the supplied source.
 */
#include "LowECoherentCommon.h"
#include "LowECoherentFiller.h"

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TParameter.h"
#include "TString.h"

#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>

namespace {

class LeafArray {
public:
    LeafArray() : leaf(0) {}

    bool Bind(TTree* tree,
              const char* branchName,
              bool required = true)
    {
        TBranch* branch =
            tree ? tree->GetBranch(branchName) : 0;

        if (!branch) {
            if (required)
                std::cerr << "Missing branch '" << branchName
                          << "' in tree '"
                          << (tree ? tree->GetName() : "null")
                          << "'" << std::endl;
            return !required;
        }

        leaf = tree->GetLeaf(branchName);

        if (!leaf && branch->GetListOfLeaves() &&
            branch->GetListOfLeaves()->GetEntries() > 0)
            leaf =
                (TLeaf*)branch->GetListOfLeaves()->At(0);

        if (!leaf) {
            if (required)
                std::cerr << "Cannot obtain leaf for branch '"
                          << branchName << "'" << std::endl;
            return !required;
        }

        std::cout << "  " << tree->GetName()
                  << "/" << branchName
                  << " [" << leaf->GetTypeName()
                  << "]" << std::endl;

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
        if (!leaf || index < 0 || index >= Size())
            return 0.0;
        return leaf->GetValue(index);
    }

    double First() const
    {
        return Size() > 0 ? At(0) : 0.0;
    }

private:
    TLeaf* leaf;
};

double ToRadians(double value, int angleUnit)
{
    return angleUnit == 2 ?
           LowECoh::DegToRad(value) : value;
}

} // namespace

void make_acqu_coherent_lowE_templates(
    const char* inputFile =
        "../Acqu_geant_He4pi0_coh_135_155.root",
    const char* outputFile =
        "coh_lowE_E0_135_155.root",
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

    TTree* tracks =
        dynamic_cast<TTree*>(input->Get("tracks"));
    TTree* tagger =
        dynamic_cast<TTree*>(input->Get("tagger"));
    TTree* setup =
        dynamic_cast<TTree*>(input->Get("setupParameters"));

    if (!tracks || !tagger || !setup) {
        std::cerr << "Required trees not found. Need tracks, "
                     "tagger and setupParameters." << std::endl;
        input->Close();
        return;
    }

    std::cout << "Input: " << inputFile << std::endl;
    std::cout << "tracks entries = "
              << tracks->GetEntries() << std::endl;
    std::cout << "tagger entries = "
              << tagger->GetEntries() << std::endl;

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
    ok &= taggerPhotonEnergy.Bind(setup,
                                  "TaggerPhotonEnergy");

    if (!ok) {
        input->Close();
        return;
    }

    setup->GetEntry(0);

    const int energyTableSize =
        taggerPhotonEnergy.Size();

    std::cout << "nTagger = "
              << (int)nTagger.First()
              << ", TaggerPhotonEnergy entries = "
              << energyTableSize << std::endl;

    LowECoherentFiller filler;

    Long64_t entries =
        std::min(tracks->GetEntries(),
                 tagger->GetEntries());

    if (maxEvents >= 0 && maxEvents < entries)
        entries = maxEvents;

    Long64_t countTaggedEvents = 0;
    Long64_t countValidTaggerEnergy = 0;
    Long64_t countTwoNeutral = 0;
    Long64_t countBestPair = 0;
    Long64_t countMassPass = 0;
    Long64_t countFilledTags = 0;

    int resolvedAngleUnit = angleUnit;
    double minMappedEnergy = 1e99;
    double maxMappedEnergy = -1e99;

    for (Long64_t entry = 0; entry < entries; ++entry) {
        tracks->GetEntry(entry);
        tagger->GetEntry(entry);

        const int nt =
            std::min((int)nTracks.First(),
            std::min(clusterEnergy.Size(),
            std::min(theta.Size(), phi.Size())));

        const int ng =
            std::min((int)nTagged.First(),
                     taggedChannel.Size());

        if (ng <= 0) continue;
        ++countTaggedEvents;

        if (nt < 2) continue;

        if (resolvedAngleUnit == 0) {
            double maxAbsTheta = 0.0;
            double maxAbsPhi = 0.0;

            for (int i = 0; i < nt; ++i) {
                maxAbsTheta =
                    std::max(maxAbsTheta,
                             std::fabs(theta.At(i)));
                maxAbsPhi =
                    std::max(maxAbsPhi,
                             std::fabs(phi.At(i)));
            }

            resolvedAngleUnit =
                (maxAbsTheta > 3.5 ||
                 maxAbsPhi > 6.5) ? 2 : 1;

            std::cout << "Angle unit inferred as "
                      << (resolvedAngleUnit == 2 ?
                          "degrees" : "radians")
                      << std::endl;
        }

        std::vector<int> neutral;

        for (int i = 0; i < nt; ++i) {
            const double energy =
                clusterEnergy.At(i)*energyScale;

            if (energy < clusterThresholdMeV)
                continue;

            if (vetoEnergy.Valid() &&
                i < vetoEnergy.Size()) {
                const double veto =
                    vetoEnergy.At(i)*energyScale;
                if (veto > vetoMaxMeV)
                    continue;
            }

            neutral.push_back(i);
        }

        if (neutral.size() < 2) continue;
        ++countTwoNeutral;

        TLorentzVector bestGamma1;
        TLorentzVector bestGamma2;
        double bestDistance = 1e99;
        bool foundPair = false;

        for (size_t ia = 0;
             ia < neutral.size(); ++ia) {
            for (size_t ib = ia + 1;
                 ib < neutral.size(); ++ib) {

                const int i = neutral[ia];
                const int j = neutral[ib];

                const double e1 =
                    clusterEnergy.At(i)*energyScale;
                const double e2 =
                    clusterEnergy.At(j)*energyScale;

                const double th1 =
                    ToRadians(theta.At(i),
                              resolvedAngleUnit);
                const double th2 =
                    ToRadians(theta.At(j),
                              resolvedAngleUnit);
                const double ph1 =
                    ToRadians(phi.At(i),
                              resolvedAngleUnit);
                const double ph2 =
                    ToRadians(phi.At(j),
                              resolvedAngleUnit);

                const TLorentzVector gamma1 =
                    LowECoh::PhotonP4(e1, th1, ph1);
                const TLorentzVector gamma2 =
                    LowECoh::PhotonP4(e2, th2, ph2);

                const double mass =
                    (gamma1 + gamma2).M();

                const double distance =
                    std::fabs(mass -
                              LowECoh::kPi0MassMeV);

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

        const double bestMass =
            (bestGamma1 + bestGamma2).M();

        if (bestMass >= mggMinMeV &&
            bestMass <= mggMaxMeV)
            ++countMassPass;

        for (int itag = 0; itag < ng; ++itag) {
            const int rawChannel =
                (int)taggedChannel.At(itag);

            const int lookupChannel =
                rawChannel + taggerChannelOffset;

            if (lookupChannel < 0 ||
                lookupChannel >= energyTableSize)
                continue;

            const double eGamma =
                taggerPhotonEnergy.At(lookupChannel)*
                energyScale;

            minMappedEnergy =
                std::min(minMappedEnergy, eGamma);
            maxMappedEnergy =
                std::max(maxMappedEnergy, eGamma);

            ++countValidTaggerEnergy;

            if (LowECoh::FindEnergyBin(eGamma) < 0)
                continue;

            filler.Fill(eGamma,
                        bestGamma1,
                        bestGamma2,
                        1.0,
                        mggMinMeV,
                        mggMaxMeV);

            ++countFilledTags;
        }

        if ((entry + 1) % 100000 == 0)
            std::cout << "Processed "
                      << entry + 1 << " / "
                      << entries << std::endl;
    }

    TFile output(outputFile, "RECREATE");

    if (output.IsZombie()) {
        std::cerr << "Cannot create "
                  << outputFile << std::endl;
        input->Close();
        return;
    }

    filler.Write(&output);

    TParameter<Long64_t>("n_entries", entries).Write();
    TParameter<Long64_t>("n_tagged_events",
                         countTaggedEvents).Write();
    TParameter<Long64_t>("n_valid_tagger_energy",
                         countValidTaggerEnergy).Write();
    TParameter<Long64_t>("n_two_neutral",
                         countTwoNeutral).Write();
    TParameter<Long64_t>("n_best_pair",
                         countBestPair).Write();
    TParameter<Long64_t>("n_mass_pass",
                         countMassPass).Write();
    TParameter<Long64_t>("n_filled_tags",
                         countFilledTags).Write();

    TParameter<double>("cluster_threshold_MeV",
                       clusterThresholdMeV).Write();
    TParameter<double>("mgg_min_MeV",
                       mggMinMeV).Write();
    TParameter<double>("mgg_max_MeV",
                       mggMaxMeV).Write();

    output.Write();
    output.Close();
    input->Close();

    std::cout << "\nWrote " << outputFile << std::endl;

    if (countValidTaggerEnergy > 0)
        std::cout << "mapped photon energies: "
                  << minMappedEnergy << " ... "
                  << maxMappedEnergy << " MeV"
                  << std::endl;

    std::cout << "tagged events      = "
              << countTaggedEvents << std::endl;
    std::cout << "valid tag energies = "
              << countValidTaggerEnergy << std::endl;
    std::cout << ">=2 neutral        = "
              << countTwoNeutral << std::endl;
    std::cout << "best pair          = "
              << countBestPair << std::endl;
    std::cout << "mgg pass           = "
              << countMassPass << std::endl;
    std::cout << "filled tags        = "
              << countFilledTags << std::endl;
}
