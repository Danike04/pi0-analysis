# Tagging efficiency and beam normalization

## Versioned production tables

The cross-section workflow accepts a text table with rows:

```text
channel  eps_tag  deps_tag
```

The retained campaign tables are stored under:

```text
inputs/tagging_efficiency/
```

Currently versioned files include:

```text
ExpBkgSub_COPP_TaggEff_31834.dat
ExpBkgSub_COPP_TaggEff_32284.dat
ExpBkgSub_COPP_TaggEff_32313.dat
```

Validated run-to-table assignments used in retained cross-checks are recorded in:

```text
config/RUN_TAGGEFF_MAP.md
```

Do not infer an undocumented run range from the filename alone.

The production tables used in this analysis are **zero-based**. The parser therefore uses the file channel directly and does not subtract one.

The same convention is used by `fig2_full_empty_subtraction_many_both.C` for FULL/EMPTY photon-flux normalization.

## Pair-Spectrometer extractor

Canonical macro:

```text
macros/efficiency/ps_tagging_eff_clean.C
```

For each channel it reads:

```text
A = FPD_ScalerAcc  (FPD_scalerAcc also accepted)
B = PairSpec_SumGated
C = PairSpec_SumGatedDly
```

and evaluates

```text
eps_tag^PS = 2500 * (B-C) / A
```

by default.

Example:

```bash
root -l -b -q '/path/to/pi0-analysis/macros/efficiency/ps_tagging_eff_clean.C("ARHist_CBTagg_RUN.root","PS_TaggEff_RUN",2500.0)'
```

The PS method was sanity-checked on the dedicated normal-current runs 32608/32609 and produced the expected approximately 0.20--0.25 channel-by-channel scale. Treat this as a campaign sanity check rather than an independent absolute calibration of every production run.

## Ladder/P2 -- PairSpec study

`macros/studies/final_lp2_pairspec.C` is a campaign-specific correlation study, not the general tagging-efficiency extractor.

## FPD scaler diagnostics

`diagnose_acqu_scalers.C` retains the historical contiguous-scaler treatment as a diagnostic example. It is **not** the production mapping. Production code must use the FPD map:

```text
inputs/fpd/FPD_855_new.dat
```

## Normalization chain

```text
versioned FPD map -> scaler index -> electron counts Ne
selected eps_tag(channel) table
Ne * eps_tag -> photon flux
```

The FPD map and tagging-table channel convention must both be correct before a cross section or FULL-EMPTY normalization is interpreted physically.
