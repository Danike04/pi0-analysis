# Fig. 2 / DeltaEpi coherent-incoherent workflow

The four retained photon-energy intervals are:

```text
E0  223-234 MeV   (224 MeV label)
E1  283-294 MeV   (294 MeV label)
E2  319-330 MeV   (320 MeV label)
E3  356-366 MeV   (366 MeV label)
```

The missing-energy variable is evaluated in the gamma+He4 CM frame:

```text
DeltaE_pi0 = E_pi0,measured^CM - E_pi0,coherent^CM
```

## Data/coherent diagnostic branch

`paper_fig2_missing_energy_root6.C` reads one data file and one coherent-MC file, reconstructs the best neutral gamma-gamma pair, applies the invariant-mass selection and fills DeltaEpi before/after a configurable DeltaPhi cut.

Default core selections include:

```text
cluster energy >= 20 MeV
veto energy <= 1 MeV
110 < mgg < 155 MeV
prompt 700-800 ns
random 450-680 ns
DeltaPhi configurable, default 10 deg
```

No missing-mass cut is applied in this macro.

## Prepared coherent templates

The high-stat combined coherent template filenames used by the retained refit are:

```text
paper_fig2_missing_energy_cut8_coh224_294_320_366new.root
paper_fig2_missing_energy_cut10_coh224_294_320_366new.root
paper_fig2_missing_energy_cut12_coh224_294_320_366new.root
```

The corresponding single-energy coherent MC source files include:

```text
Acqu_geant_He4pi0_224coh.root
Acqu_geant_He4pi0_294coh.root
Acqu_geant_He4pi0_320coh.root
Acqu_geant_He4pi0_366coh.root
```

These large ROOT files and prepared template ROOT files are not versioned in Git.

## 3He+n breakup template

The retained high-stat prepared input is:

```text
fig2_he3n_all_available_bins_wide_highstat294_320_366.root
```

with histogram interface:

```text
inc_deltaE_E<0..3>_nodphi
inc_deltaE_E<0..3>_dphi8
inc_deltaE_E_E<0..3>_dphi10
inc_deltaE_E<0..3>_dphi12
```

Historical incoherent source files include:

```text
Acqu_geant_He4pi0_294_25M.root
Acqu_geant_He4pi0_320.root
Acqu_geant_He4pi0_366.root
```

The names without `coh` are the incoherent/breakup samples in this analysis.

## t+p template

The prepared file is:

```text
fig2_tp_E3_366_wide.root
```

It follows the E0--E3 naming interface but is intentionally populated only for E3.

## FULL-EMPTY preparation

Canonical macro:

```text
fig2_full_empty_subtraction_many_both.C
```

It accepts one list of FULL data files and one list of EMPTY data files, calculates the photon flux for each set and subtracts

```text
alpha = Phi_full / Phi_empty
H_sub = H_full - alpha * H_empty
```

The output interface consumed by the refit is:

```text
he4_before_E0..E3
he4_cut8_E0..E3
he4_cut10_E0..E3
he4_cut12_E0..E3
```

### Normalization requirements

The FULL-EMPTY macro requires the real non-contiguous FPD map and the appropriate tagging-efficiency table for production. Retained small campaign inputs are versioned under:

```text
inputs/fpd/
inputs/tagging_efficiency/
```

Run-to-tagging-efficiency assignments that were explicitly validated are recorded in:

```text
config/RUN_TAGGEFF_MAP.md
```

The tagging-efficiency tables are interpreted as **zero-based**, matching the validated cross-section convention.

During final repository review the historical `fileChannel - 1` shift in this FULL-EMPTY parser was corrected to `fileChannel`. Therefore precomputed FULL-EMPTY products made with the shifted parser should be regenerated before they are treated as final normalized physics outputs.

## Final retained refit

`fig2_refit_coherent_incoherent_fullminus_empty.C` reads the FULL-EMPTY data histograms and the coherent/breakup templates.

It performs independent non-negative two-template comparisons:

```text
data = a_coh * coherent + a_inc * (3He+n)
```

and, when the t+p histogram is populated,

```text
data = a_coh * coherent + a_inc * (t+p)
```

It does **not** perform one simultaneous coherent + 3He+n + t+p fit.

Default ranges:

```text
fit   [-60,+40] MeV
tail  [-60,-20] MeV
```

A portable example call is provided in `QUICKSTART.md`.

## Standardized-template branch

The repository also preserves:

```text
paper_fig2_missing_energy_root6.C
    -> quantify_fig2_cut_stability_root6.C
    -> fig2_data_coherent.root

make_acqu_incoherent_fig2_templates.C
    -> fig2_incoherent_templates.root

validate_fig2_outputs.C
overlay_fig2_with_incoherent.C
```

This branch is useful for diagnostic/cut-stability studies and uses `Fig2Common.h`, `Fig2IncoherentFiller.h` and `Fig2SaveInputs.h` for shared definitions.

## Connection to the cross-section selection

The DeltaEpi coherent/incoherent distributions are the appropriate place to study breakup contamination. They should inform a physics selection; the differential-cross-section selection should not be chosen simply to optimize agreement with the paper reference.
