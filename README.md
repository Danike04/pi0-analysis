# A2 He4 pi0 analysis macros

ROOT/C++ analysis code for the 2026 A2/MAMI study of neutral-pion photoproduction on He4.

The repository contains the retained production analysis, efficiency tools, Fig. 2 / DeltaEpi workflow, diagnostics, and clearly separated development studies. Large data and Monte-Carlo ROOT files are not part of the repository.

## Repository layout

```text
macros/
  cross_section/   production cross-section extraction
  efficiency/      detection and tagging-efficiency tools
  fig2/            DeltaEpi, FULL-EMPTY and template analysis
  diagnostics/     reconstruction and normalization diagnostics
  studies/         campaign-specific or development studies

inputs/
  fpd/             validated FPD scaler map
  tagging_efficiency/
  detection_efficiency/

config/
  RUN_TAGGEFF_MAP.md

scripts/
  build_hybrid_efficiency.py

docs/
  QUICKSTART.md
  INPUTS_AND_FORMATS.md
  CROSS_SECTION_WORKFLOW.md
  FIG2_WORKFLOW.md
  TAGGING_EFFICIENCY_WORKFLOW.md
  DIAGNOSTICS_AND_STUDIES.md
  VALIDATION_STATUS.md
```

## Requirements

- ROOT 6
- a C++ compiler usable by ROOT/ACLiC
- Acqu ROOT files with the trees required by the chosen macro
- external experimental and Monte-Carlo ROOT files required by the workflow
- the small campaign auxiliary inputs versioned under `inputs/`

Typical ROOT usage is:

```bash
root -l
```

then, for example:

```cpp
.L /path/to/pi0-analysis/macros/cross_section/paper_style_diff_xs_from_acqu.C+
```

Run production macros from a separate output directory because many historical macros intentionally keep fixed output filenames.

## Versioned auxiliary inputs

Small campaign inputs required by the retained analysis are stored directly in the repository:

```text
inputs/fpd/FPD_855_new.dat
inputs/tagging_efficiency/
inputs/detection_efficiency/eps_det_paper_bins_final_highstat.txt
```

Validated run-to-tagging-efficiency assignments are recorded in:

```text
config/RUN_TAGGEFF_MAP.md
```

Large experimental and Monte-Carlo ROOT files remain external.

The retained high-statistics paper-bin detection-efficiency table can be rebuilt from broad and dedicated coherent-MC efficiency tables using:

```text
scripts/build_hybrid_efficiency.py
```

## Main production chain

```text
Acqu DATA
   -> pi0 reconstruction and prompt/random subtraction
   -> real FPD scaler mapping
   -> channel-dependent tagging efficiency
   -> eps_det(Egamma,theta_cm) from coherent MC
   -> paper_style_diff_xs_from_acqu.C
   -> d(sigma)/dOmega
```

Core files:

```text
macros/cross_section/cross_section_from_acqu.C
macros/cross_section/paper_style_diff_xs_from_acqu.C
macros/efficiency/eps_det_paper_bins_from_acqu_geant.C
```

The differential-cross-section analysis uses 17 laboratory photon-energy bins between 201 and 401 MeV and, by default, 5-degree pion CM angular bins from 5 to 150 degrees.

The normalization is

```text
dsigma/dOmega = Y / (fluxDen * eps_det * DeltaOmega)
```

with prompt-random-subtracted yield `Y`.

## Important conventions and validated corrections

### FPD mapping

The production FPD scaler mapping is non-contiguous. The production cross-section and FULL-EMPTY normalization code require the real FPD map and do not use `2000 + channel` as a production fallback.

### Tagging-efficiency channel numbering

The production tagging-efficiency table used in this analysis is **zero-based**:

```text
channel  eps_tag  deps_tag
```

The retained parsers therefore use the table channel directly:

```cpp
int ch = chFile;
```

This applies both to the cross-section family and to the FULL-EMPTY Fig. 2 flux-normalization path.

### Detection efficiency

`eps_det_paper_bins_from_acqu_geant.C` uses

```text
Ngen(E,theta) = Ngen(E) * DeltaOmega/(4*pi)
```

and therefore assumes isotropic generation in the CM frame.

If a data yield uses a missing-mass cut, regenerate the efficiency with the same missing-mass selection.

### DeltaEpi

For the Fig. 2 workflow,

```text
DeltaE_pi0 = E_pi0,measured^CM - E_pi0,coherent^CM
```

in the gamma+He4 centre-of-mass frame.

## Where to start

Read, in order:

```text
docs/QUICKSTART.md
docs/INPUTS_AND_FORMATS.md
docs/VALIDATION_STATUS.md
```

Then use either:

```text
docs/CROSS_SECTION_WORKFLOW.md
docs/FIG2_WORKFLOW.md
```

## Production versus studies

The recommended production macros are under `macros/cross_section`, `macros/efficiency` and `macros/fig2`. Files under `macros/studies` are retained as useful cross-checks or campaign-specific development work and should not be substituted silently for the production workflow.

## Validation

A macro that executes successfully is not automatically a physics-validated result. `docs/VALIDATION_STATUS.md` distinguishes:

```text
RUNS
SANITY CHECKED
VALIDATED
QUESTIONABLE
NOT YET TESTED
```

Use that file before treating an output as final physics.

Just a test
