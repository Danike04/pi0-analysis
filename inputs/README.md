# Versioned auxiliary inputs

This directory contains small auxiliary inputs required to reproduce the
retained analysis. Large experimental and Monte-Carlo ROOT files are not
stored in Git.

## FPD map

    fpd/FPD_855_new.dat

This is the real non-contiguous FPD channel-to-scaler mapping used by the
production normalization.

## Tagging efficiency

    tagging_efficiency/ExpBkgSub_COPP_TaggEff_31834.dat
    tagging_efficiency/ExpBkgSub_COPP_TaggEff_32284.dat
    tagging_efficiency/ExpBkgSub_COPP_TaggEff_32313.dat

Format:

    channel  eps_tag  deps_tag

The channel numbering is zero-based.

See:

    ../config/RUN_TAGGEFF_MAP.md

for explicitly validated run-to-table assignments.

## Detection efficiency

    detection_efficiency/eps_det_paper_bins_final_highstat.txt

This is the paper-bin efficiency table used in the final independent-sample
cross-section check.

The table uses the broad coherent MC for the full 17-bin coverage, with
higher-statistics dedicated coherent MC replacing the corresponding rows for:

    223-234 MeV  -> 224coh
    283-294 MeV  -> 294coh
    319-330 MeV  -> 320coh
    356-366 MeV  -> 366coh

The large MC ROOT source files are external to the repository.
