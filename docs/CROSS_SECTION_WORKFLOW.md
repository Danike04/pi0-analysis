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

The FPD mapping is non-contiguous. The production cross-section macros abort if the map cannot be read rather than using the historical `2000 + channel` approximation.

The production tagging-efficiency table is zero-based:

```cpp
int ch = chFile;
```

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

The broad coherent MC historically used for the full 17-bin table is named:

```text
Acqu_geant_He4pi0.root
```

Single-energy coherent files used for reference checks include:

```text
Acqu_geant_He4pi0_224coh.root
Acqu_geant_He4pi0_294coh.root
Acqu_geant_He4pi0_320coh.root
Acqu_geant_He4pi0_366coh.root
```

## Selection consistency

The historical efficiency macro defaults to `applyMMCut=false`. If the data cross section uses an MM window, production use should explicitly enable the same MM cut in the MC efficiency.

Do not choose the MM window by optimizing agreement with the embedded paper points. The coherent/breakup studies are the appropriate place to motivate the physics selection.

## Example

See `QUICKSTART.md` for a complete shell/ROOT example using explicit data, MC, FPD, tagging-efficiency and paper-bin efficiency paths.

## GoAT alternative

`eps_det_paper_bins_from_goat.C` maps GoAT generated/reconstructed histograms into the same text-table interface used by `paper_style_diff_xs_from_acqu.C`.

## Historical channel-efficiency fallback

`cross_section_from_acqu.C` can look next to the input ROOT file for `eps_det_reco_fast_from_acqu_geant.txt` when no channel-efficiency file is supplied. No producer of that exact historical file is part of this repository. For the paper-style differential cross section, the recommended production route is the explicit paper-bin efficiency table.
