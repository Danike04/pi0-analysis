/*
 * inspect_scalers.C
 *
 * Purpose
 * -------
 * Low-level diagnostic for the Acqu "scalers" tree. It compares the first
 * and last scaler-tree entries for all 8706 scaler indices and prints the
 * non-zero/changing indices sorted by |last-first|.
 *
 * Input
 * -----
 * An Acqu ROOT file containing the "scalers" tree and its scalers[8706]
 * branch. The macro also reads eventNumber, eventID, nmr and curr.
 *
 * Output
 * ------
 * Text printed to stdout only (up to the first 200 changing scaler indices).
 *
 * Usage
 * -----
 * root -l -b -q 'inspect_scalers.C("Acqu_CBTagg_31836.root")'
 *
 * This macro is useful when identifying or checking scaler indices. It is
 * not part of the production cross-section calculation.
 *
 * The code below is the original analysis code; only this documentation
 * header has been added.
 */

#include <iostream>
#include <vector>
#include <algorithm>

void inspect_scalers(const char* filename = "Acqu_CBTagg_31836.root")
{
    TFile* f = TFile::Open(filename);

    if (!f || f->IsZombie()) {
        std::cout << "Error: cannot open file " << filename << std::endl;
        return;
    }

    TTree* t = (TTree*) f->Get("scalers");

    if (!t) {
        std::cout << "Error: tree scalers not found" << std::endl;
        return;
    }

    const int NSCAL = 8706;

    UInt_t scalers[NSCAL];
    Int_t eventNumber;
    Int_t eventID;
    Float_t nmr;
    Float_t curr;

    t->SetBranchAddress("eventNumber", &eventNumber);
    t->SetBranchAddress("eventID", &eventID);
    t->SetBranchAddress("scalers", scalers);
    t->SetBranchAddress("nmr", &nmr);
    t->SetBranchAddress("curr", &curr);

    Long64_t N = t->GetEntries();

    std::vector<UInt_t> first(NSCAL);
    std::vector<UInt_t> last(NSCAL);

    t->GetEntry(0);
    for (int i = 0; i < NSCAL; i++) {
        first[i] = scalers[i];
    }

    t->GetEntry(N - 1);
    for (int i = 0; i < NSCAL; i++) {
        last[i] = scalers[i];
    }

    struct Item {
        int idx;
        UInt_t first;
        UInt_t last;
        Long64_t diff;
    };

    std::vector<Item> changed;

    for (int i = 0; i < NSCAL; i++) {
        Long64_t diff = (Long64_t)last[i] - (Long64_t)first[i];

        if (diff != 0 || last[i] != 0 || first[i] != 0) {
            changed.push_back({i, first[i], last[i], diff});
        }
    }

    std::sort(changed.begin(), changed.end(),
              [](const Item& a, const Item& b) {
                  return std::abs(a.diff) > std::abs(b.diff);
              });

    std::cout << "========================================" << std::endl;
    std::cout << "Scaler inspection" << std::endl;
    std::cout << "File = " << filename << std::endl;
    std::cout << "Scaler entries = " << N << std::endl;
    std::cout << "Non-zero/changing scaler channels = " << changed.size() << std::endl;
    std::cout << "Showing first 200 sorted by |last-first|" << std::endl;
    std::cout << "----------------------------------------" << std::endl;

    int nshow = std::min((int)changed.size(), 200);

    for (int k = 0; k < nshow; k++) {
        std::cout << "idx = " << changed[k].idx
                  << "   first = " << changed[k].first
                  << "   last = " << changed[k].last
                  << "   diff = " << changed[k].diff
                  << std::endl;
    }

    std::cout << "========================================" << std::endl;
}
