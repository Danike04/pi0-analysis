/*
 * inspect_mc_tree.C
 *
 * PURPOSE
 *   Generic ROOT-file inspection utility used during the incoherent-MC work.
 *   Lists every object in the file and, for each TTree, prints entry count,
 *   branches and leaf types.
 *
 * INPUT/OUTPUT
 *   Input ROOT filename is configurable; output is printed to stdout only.
 *
 * USAGE EXAMPLE
 *   root -l -b -q 'inspect_mc_tree.C("ppnn_mc.root")'
 *
 * REPOSITORY STATUS
 *   Generic MC diagnostic utility; executable code unchanged.
 */
#include "TFile.h"
#include "TKey.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TClass.h"
#include "TObject.h"
#include <iostream>

void inspect_mc_tree(const char* inputFile = "ppnn_mc.root")
{
    TFile* f = TFile::Open(inputFile, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Cannot open " << inputFile << std::endl;
        return;
    }

    std::cout << "\nFile: " << inputFile << "\n" << std::endl;
    TIter nextKey(f->GetListOfKeys());
    TKey* key = 0;

    while ((key = (TKey*)nextKey())) {
        std::cout << key->GetClassName() << "  " << key->GetName() << std::endl;

        TClass* cl = TClass::GetClass(key->GetClassName());
        if (!cl || !cl->InheritsFrom(TTree::Class())) continue;

        TTree* t = (TTree*)f->Get(key->GetName());
        std::cout << "  entries: " << t->GetEntries() << std::endl;

        TObjArray* branches = t->GetListOfBranches();
        for (int ib = 0; ib < branches->GetEntries(); ++ib) {
            TBranch* b = (TBranch*)branches->At(ib);
            std::cout << "  branch: " << b->GetName();
            const char* cls = b->GetClassName();
            if (cls && cls[0]) std::cout << "   class=" << cls;

            TObjArray* leaves = b->GetListOfLeaves();
            if (leaves && leaves->GetEntries() > 0) {
                std::cout << "   leaves=";
                for (int il = 0; il < leaves->GetEntries(); ++il) {
                    TLeaf* leaf = (TLeaf*)leaves->At(il);
                    if (il) std::cout << ", ";
                    std::cout << leaf->GetName() << ":" << leaf->GetTypeName();
                }
            }
            std::cout << std::endl;
        }
    }

    f->Close();
}
