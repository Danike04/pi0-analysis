# Cross-section workflow

## Production macros

```text
macros/cross_section/cross_section_from_acqu.C
macros/cross_section/paper_style_diff_xs_from_acqu.C
macros/efficiency/eps_det_paper_bins_from_acqu_geant.C
macros/efficiency/eps_det_paper_bins_from_goat.C
```

`cross_section_from_acqu.C` extracts a channel-by-channel cross section. `paper_style_diff_xs_from_acqu.C` produces the CM angular differential cross section in the 17 paper-style photon-energy bins.

## Photon normalization

Production normalization requires:

1. the real FPD channel-to-scaler map;
2. the channel-dependent tagging-efficiency table;
3. the target thickness configured in the macro call.

The versioned campaign inputs are under:

```text
inputs/fpd/
inputs/tagging_efficiency/
```

Validated run-to-tagging-efficiency assignments used in retained cross-checks are recorded in:

```text
config/RUN_TAGGEFF_MAP.md
```

The FPD mapping is non-contiguous. The production cross-section macros abort if the map cannot be read rather than using the historical `2000 + channel` approximation.

The production tagging-efficiency table is zero-based:

```cpp
int ch = chFile;
```

### Valid normalization-channel mask

The differential yield and photon-flux normalization must use the same set of valid tagger channels.

`paper_style_diff_xs_from_acqu.C` therefore builds a normalization-channel mask from the scaler, tagging-efficiency and normalization requirements and applies the same mask to the event yield. A channel rejected from the photon-flux denominator cannot contribute to the numerator.

This consistency correction prevents the numerator and denominator from being formed from different tagger-channel sets.

## Differential-cross-section binning

Laboratory photon-energy bins (MeV):

```text
201-210, 211-222, 223-234, 235-246, 247-258, 259-270,
271-282, 283-294, 295-308, 309-318, 319-330, 331-342,
343-355, 356-366, 367-378, 379-390, 391-401
```

The reconstructed pion angle is evaluated in the gamma+He4 CM frame. The default angular interval is 5--150 degrees in 5-degree bins.

The yield is obtained from prompt-random-subtracted missing-mass histograms integrated over the configured MM window.

```text
dsigma/dOmega = Y / (fluxDen * eps_det(E,theta) * DeltaOmega)
```

The Phi-min / DeltaPhi quantity is diagnostic only in this macro; no Fig. 2 DeltaPhi cut is applied to the differential-cross-section yield.

## Detection efficiency from Acqu/Geant

`eps_det_paper_bins_from_acqu_geant.C` reconstructs coherent MC in the same `(Egamma,theta_cm)` bins.

The generated denominator is

```text
Ngen(E,theta) = Ngen(E) * DeltaOmega/(4*pi)
```

so the method assumes isotropic generation in the CM frame.

The broad coherent MC used for complete 17-bin coverage is:

```text
Acqu_geant_He4pi0.root
```

Dedicated coherent files used for the retained high-statistics efficiency table are:

```text
Acqu_geant_He4pi0_224coh.root
Acqu_geant_He4pi0_294coh.root
Acqu_geant_He4pi0_320coh.root
Acqu_geant_He4pi0_366coh.root
```

The retained hybrid table uses dedicated coherent MC for:

```text
223-234 MeV -> 224coh
283-294 MeV -> 294coh
319-330 MeV -> 320coh
356-366 MeV -> 366coh
```

and the broad coherent MC for the remaining 13 energy bins.

The 224coh sample is monoenergetic near 224.186 MeV and is used as the available proxy for the 223--234 MeV paper bin.

The retained table is versioned as:

```text
inputs/detection_efficiency/eps_det_paper_bins_final_highstat.txt
```

It can be reconstructed from the broad and dedicated efficiency tables using:

```text
scripts/build_hybrid_efficiency.py
```

The numerical rows produced by that script were checked to match the retained high-statistics table exactly.

## Selection consistency

The historical efficiency macro defaults to `applyMMCut=false`. If the data cross section uses an MM window, production use should explicitly enable the same MM cut in the MC efficiency.

The retained high-statistics table was produced with:

```text
110 < mgg < 155 MeV
MM in [-10,40] MeV
```

Do not choose the MM window by optimizing agreement with the embedded paper points. The coherent/breakup studies are the appropriate place to motivate the physics selection.

## Combining multiple runs

Do not average cross sections from separate runs and do not use an `hadd`-merged scaler tree for production normalization.

For each `(Egamma,theta_cm)` bin combine the run-level quantities as:

```text
Y_total       = sum_r Y_r
fluxDen_total = sum_r fluxDen_r
```

and then evaluate:

```text
dsigma/dOmega = Y_total / (fluxDen_total * eps_det * DeltaOmega)
```

This preserves the normalization associated with each run and each run-period tagging-efficiency table.

## Example

See `QUICKSTART.md` for a complete shell/ROOT example using the versioned FPD, tagging-efficiency and high-statistics paper-bin efficiency inputs.

## GoAT alternative

`eps_det_paper_bins_from_goat.C` maps GoAT generated/reconstructed histograms into the same text-table interface used by `paper_style_diff_xs_from_acqu.C`.

## Historical channel-efficiency fallback

`cross_section_from_acqu.C` can look next to the input ROOT file for `eps_det_reco_fast_from_acqu_geant.txt` when no channel-efficiency file is supplied. No producer of that exact historical file is part of this repository. For the paper-style differential cross section, the recommended production route is the explicit paper-bin efficiency table.
