# Inputs and formats

Large experimental and Monte-Carlo ROOT files are not stored in this repository.

## Required inputs

Depending on the workflow, the analysis requires:

- Acqu ROOT data files
- coherent/incoherent Monte-Carlo ROOT files
- the real FPD scaler map
- the tagging-efficiency table

Production Acqu files typically use the trees:

- `tracks`
- `tagger`
- `setupParameters`
- `scalers`

## FPD mapping

The FPD scaler mapping is non-contiguous.

Production normalization must use the real FPD map and must not assume:

    scaler = 2000 + channel

## Tagging efficiency

The tagging-efficiency table has the format:

    channel  eps_tag  deps_tag

Channel numbering is zero-based. The channel read from the file is therefore
used directly, without subtracting one.

## Detection efficiency

The differential cross-section workflow uses:

    eps_det(Egamma, theta_cm)

generated with:

    macros/efficiency/eps_det_paper_bins_from_acqu_geant.C

For the retained broad-MC workflow:

    Ngen(E,theta) = Ngen(E) * DeltaOmega/(4*pi)

which assumes isotropic generation in the CM frame.

Data and reconstructed MC must use compatible selections. In particular, if
the data use a missing-mass cut, the efficiency must be regenerated using the
same missing-mass cut.

## Combining multiple runs

Do not average the cross sections of different runs.

Combine the prompt-random-subtracted yields and photon-flux normalizations:

    Y_total       = sum_r Y_r
    fluxDen_total = sum_r fluxDen_r

and calculate:

    dsigma/dOmega =
        Y_total / (fluxDen_total * eps_det * DeltaOmega)

This preserves the normalization of every run and avoids ambiguities caused by
concatenating scaler trees with `hadd`.

## External data

Paths to data, MC, FPD maps and tagging-efficiency files are installation- and
campaign-dependent and should be supplied explicitly when running the macros.
