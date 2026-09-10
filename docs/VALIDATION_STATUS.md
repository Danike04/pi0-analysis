# Validation status

This file records what was exercised during the Mainz analysis and the final repository review. Execution, numerical sanity checks and physics validation are deliberately kept separate.

## Status labels

- **RUNS** — executes in the tested ROOT environment and writes the expected outputs.
- **SANITY CHECKED** — output also passed a limited numerical/physics consistency check.
- **VALIDATED** — compared with a trusted reference for the stated purpose.
- **QUESTIONABLE** — executable, but its intended physics use still has a known caveat.
- **NOT YET TESTED** — source inspected/documented only.

A macro producing files without errors is not automatically physics-validated.

## Reference reconstruction run

The primary reconstruction reference was `Acqu_CBTagg_31837.root`.

### `check_pi0_from_acqu.C`

**Status: RUNS / workflow PASS**

The macro produced the expected outputs. The best-pair invariant-mass peak was observed near 128.5 MeV and accepted as the working reconstruction check for this analysis pass.

## FPD mapping and electron counts

**Status: SANITY CHECKED / mapping PASS**

The validated FPD map (`FPD_855_new.dat`) is non-contiguous. The resulting electron counts were smooth over the tested range. Production normalization must use the map rather than a `2000 + channel` assumption.

The final repository version makes the cross-section and FULL-EMPTY production paths fail closed if the required FPD map cannot be read.

`diagnose_acqu_scalers.C` still contains the historical contiguous treatment and is retained only as a low-level diagnostic.

## Tagging-efficiency channel indexing

**Status: VALIDATED for the production table convention**

The production table `ExpBkgSub_COPP_TaggEff_31834.dat` is zero-based. The cross-section-family parser was corrected from `chFile - 1` to:

```cpp
int ch = chFile;
```

The correction was directly checked with source-table channels 54 and 77.

The corrected files include:

```text
macros/cross_section/cross_section_from_acqu.C
macros/cross_section/paper_style_diff_xs_from_acqu.C
macros/diagnostics/background_diagnostics_from_acqu.C
macros/studies/cross_section/cross_section_from_acqu_with_2d.C
macros/studies/cross_section/angular_cross_section_from_acqu.C
```

### Final repository review: Fig. 2 FULL-EMPTY parser

The same zero-based production table is read by `fig2_full_empty_subtraction_many_both.C`. Final review found that this parser still contained a historical `fileChannel - 1` shift. It has now been corrected to `fileChannel`.

**Status: software correction applied; FULL-EMPTY normalized outputs should be regenerated before final physics use.**

This correction is logically fixed by the already validated table convention, but the full long FULL/EMPTY production pass has not been rerun inside this review environment.

## Event-selection/background diagnostics

### `background_diagnostics_from_acqu.C`

**Status: SANITY CHECKED on full run 31837**

Recorded validation integrals included:

```text
h_mgg_best_selected             1.61394e6
h_mm_before_mgg_cut             1.13055e6
h_mm_after_mgg_cut              625657
h_delta_epi_before_mgg_cut      1.13055e6
h_delta_epi_after_mgg_cut       625657
```

The macro writes `sigma_from_acqu.*` as well as diagnostics, so it should be run in a separate output directory.

## Detection efficiency

### `eps_det_paper_bins_from_acqu_geant.C`

**Status: RUNS / SANITY CHECKED**

Reference coherent samples at 224, 294, 320 and 366 MeV were exercised. The broad `Acqu_geant_He4pi0.root` sample produced a complete 17-energy-bin efficiency table.

The denominator is

```text
Ngen(E,theta) = Ngen(E) * DeltaOmega/(4*pi)
```

and therefore assumes isotropic CM generation.

If the data yield uses a missing-mass selection, the MC efficiency must be regenerated with the same selection.

### Retained high-statistics hybrid table

**Status: SANITY CHECKED / reproducibility PASS**

The retained table is:

```text
inputs/detection_efficiency/eps_det_paper_bins_final_highstat.txt
```

It uses the broad coherent MC for the full 17-bin coverage, with the following energy bins replaced by dedicated higher-statistics coherent MC results:

```text
223-234 MeV -> 224coh
283-294 MeV -> 294coh
319-330 MeV -> 320coh
356-366 MeV -> 366coh
```

The broad and dedicated tables used for this replacement were generated with the same retained reconstruction selections, including the missing-mass cut `[-10,40] MeV`.

`scripts/build_hybrid_efficiency.py` reproduces the replacement deterministically and checks the angular-bin keys before replacing rows. Its numerical output was compared with the retained table and gave an exact match.

## Differential cross section

### `paper_style_diff_xs_from_acqu.C`

**Status: RUNS / SANITY CHECKED on full statistics; normalization consistency checked**

The chain has been exercised with the real FPD mapping, corrected zero-based tagging efficiency and a complete paper-bin efficiency table. All 17 photon-energy bins can be populated.

A consistency issue was found while comparing two non-overlapping data samples: channels rejected from the photon-flux normalization could still contribute to the event yield. The production macro now constructs one `validNormChannel` mask and applies it to both numerator and denominator.

The corrected workflow was reprocessed independently for:

```text
runs 31840-31849 -> ExpBkgSub_COPP_TaggEff_31834.dat
runs 32330-32339 -> ExpBkgSub_COPP_TaggEff_32313.dat
```

using the retained high-statistics hybrid detection-efficiency table.

For the 487 common `(Egamma,theta_cm)` points, the statistical comparison between the two extracted data samples gave:

```text
mean pull               = -0.0062
RMS pull                =  1.1303
chi2/N                  =  1.2750
median xs323/xs318      =  1.0050
```

For the 283--294 MeV bin, where the broad-MC efficiency had produced a strong common point-to-point structure, the local fluctuation correlation decreased from approximately `0.90` to approximately `0.23` after using the dedicated higher-statistics coherent MC.

These two run groups are statistically independent experimental event samples. The extracted cross sections still share common systematic corrections, including the detection efficiency, so this comparison is a consistency check rather than a proof of complete systematic independence.

The pipeline is technically operational. The final coherent-event MM selection remains a physics choice and should not be chosen by tuning agreement with the embedded paper reference.

The macro computes DeltaPhi/Phi-min diagnostics but does not apply a DeltaPhi cut to the differential-cross-section yield.

## Pair-Spectrometer tagging efficiency

### `ps_tagging_eff_clean.C`

**Status: SANITY CHECKED for the 2026 PS measurement method**

The macro uses

```text
eps_tag^PS = 2500 * (B-C) / A
```

with `FPD_ScalerAcc`/`FPD_scalerAcc`, `PairSpec_SumGated`, and `PairSpec_SumGatedDly`.

Dedicated normal-current runs 32608/32609 produced the expected approximately 0.20--0.25 channel-by-channel scale. This supports the method, but it should not be confused with an independent absolute calibration of every run in the production cross-section table.

## Fig. 2 / DeltaEpi

### Kinematics and prepared-template interface

**Status: SANITY CHECKED / operational**

The retained workflow uses four reference energy intervals and DeltaPhi cuts of 8, 10 and 12 degrees. Prepared coherent, 3He+n and E3 t+p templates were successfully used by the refit workflow.

The refit performs independent two-template comparisons:

```text
data = coherent + 3He+n
```

and, when available,

```text
data = coherent + t+p
```

It does not perform one simultaneous coherent + 3He+n + t+p fit.

Default fit range: `[-60,+40] MeV`.

### Precomputed FULL-EMPTY files

Prepared `he4_*` histogram files remain useful for fast code/interface checks. Because the FULL-EMPTY tagging-efficiency channel parser was corrected during final repository review, old flux-normalized products made with the shifted parser should be treated as historical and regenerated before final normalized physics use.

## Overall repository state

The retained package contains an operational analysis chain for reconstruction, FPD/scaler normalization, tagging efficiency, Acqu/Geant detection efficiency, differential cross section, background diagnostics and Fig. 2 template studies.

Small campaign inputs required for the retained production cross-section workflow are now versioned under `inputs/`, while large experimental and Monte-Carlo ROOT files remain external.

Known remaining physics caveats are explicit rather than hidden in path or parser assumptions.
