# Inputs and formats

Large experimental and Monte-Carlo ROOT files are not stored in this repository. Small campaign auxiliary inputs required by the retained production workflow are versioned under `inputs/`.

## Required inputs

Depending on the workflow, the analysis requires:

- Acqu ROOT data files
- coherent/incoherent Monte-Carlo ROOT files
- the real FPD scaler map
- the tagging-efficiency table
- for the retained differential-cross-section result, the paper-bin detection-efficiency table

Production Acqu files typically use the trees:

- `tracks`
- `tagger`
- `setupParameters`
- `scalers`

## Versioned auxiliary inputs

The retained campaign files are stored as:

```text
inputs/fpd/FPD_855_new.dat
inputs/tagging_efficiency/ExpBkgSub_COPP_TaggEff_31834.dat
inputs/tagging_efficiency/ExpBkgSub_COPP_TaggEff_32284.dat
inputs/tagging_efficiency/ExpBkgSub_COPP_TaggEff_32313.dat
inputs/detection_efficiency/eps_det_paper_bins_final_highstat.txt
```

Validated run-to-tagging-efficiency assignments are recorded in:

```text
config/RUN_TAGGEFF_MAP.md
```

Do not infer an undocumented run-range assignment from a tagging-efficiency filename alone.

## FPD mapping

The FPD scaler mapping is non-contiguous.

Production normalization must use the real FPD map and must not assume:

```text
scaler = 2000 + channel
```

The retained production map is:

```text
inputs/fpd/FPD_855_new.dat
```

## Tagging efficiency

The tagging-efficiency table has the format:

```text
channel  eps_tag  deps_tag
```

Channel numbering is zero-based. The channel read from the file is therefore used directly, without subtracting one.

The same convention is used by the differential-cross-section normalization and by the FULL-EMPTY Fig. 2 normalization path.

## Detection efficiency

The differential cross-section workflow uses:

```text
eps_det(Egamma, theta_cm)
```

generated with:

```text
macros/efficiency/eps_det_paper_bins_from_acqu_geant.C
```

For each coherent MC input, the generated denominator is estimated as:

```text
Ngen(E,theta) = Ngen(E) * DeltaOmega/(4*pi)
```

which assumes isotropic generation in the CM frame.

Data and reconstructed MC must use compatible selections. In particular, if the data use a missing-mass cut, the efficiency must be generated using the same missing-mass cut.

The retained high-statistics table is:

```text
inputs/detection_efficiency/eps_det_paper_bins_final_highstat.txt
```

It uses the broad coherent MC for complete 17-bin coverage and dedicated coherent MC in four energy intervals. The deterministic replacement step is implemented by:

```text
scripts/build_hybrid_efficiency.py
```

See `CROSS_SECTION_WORKFLOW.md` and `QUICKSTART.md` for the exact retained configuration.

## Combining multiple runs

Do not average the cross sections of different runs.

Combine the prompt-random-subtracted yields and photon-flux normalizations:

```text
Y_total       = sum_r Y_r
fluxDen_total = sum_r fluxDen_r
```

and calculate:

```text
dsigma/dOmega =
    Y_total / (fluxDen_total * eps_det * DeltaOmega)
```

This preserves the normalization of every run and avoids ambiguities caused by concatenating scaler trees with `hadd`.

## External data

Only the large event-level inputs remain external to the repository: experimental Acqu ROOT files, broad coherent MC, dedicated coherent MC and incoherent/breakup MC.

Their local paths are installation-dependent and should be supplied explicitly when running the macros. Small FPD, tagging-efficiency and retained detection-efficiency inputs should normally be taken from the versioned files under `inputs/` when reproducing the retained analysis.
