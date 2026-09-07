# Diagnostics and studies

The repository separates production analysis from diagnostic and development
code.

## Production

The recommended production workflows are contained in:

    macros/cross_section/
    macros/efficiency/
    macros/fig2/

## Diagnostics

Files under:

    macros/diagnostics/

are intended for reconstruction checks, scaler inspection, background studies,
MC/data comparisons and other validation tasks.

They are not automatically part of the production cross-section chain.

## Studies

Files under:

    macros/studies/

contain historical studies, alternative approaches and campaign-specific
cross-checks.

They are retained for reference but must not be substituted silently for
production macros.

Some may contain special cuts, assumptions, input files or normalization
choices specific to the study for which they were written.

## Validation

See:

    docs/VALIDATION_STATUS.md

The retained status terminology is:

    RUNS
    SANITY CHECKED
    VALIDATED
    QUESTIONABLE
    NOT YET TESTED

A macro compiling or running successfully does not by itself imply physics
validation.
