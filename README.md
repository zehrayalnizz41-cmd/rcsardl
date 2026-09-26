# RCS-ARDL for Stata

**Version 0.1.5 — 26 September 2026**

`rcsardl` implements the validated point-estimation algorithm for the **Robust Cross-Sectionally Augmented Autoregressive Distributed Lag (RCS-ARDL)** approach developed in:

> Zehra Yalnız (2026), *Robust Estimation of Equilibrium Adjustment in Cross-Sectionally Dependent Dynamic Panels: The RCS-ARDL Approach*.

## What the command does

The command preserves the conventional CS-ARDL long-run coefficient vector and robustifies only the equilibrium-adjustment component.

The implemented sequence is:

1. Estimate conventional CS-ARDL with `xtdcce2`.
2. Retain the conventional CS-ARDL long-run coefficients unchanged.
3. Construct the equilibrium error from those long-run coefficients.
4. Estimate a robust error-correction regression separately for each panel unit using `rreg`.
5. Use cross-sectional medians of Δy, Δx, and lagged ECT as robust common controls.
6. Retain dynamically admissible unit-specific adjustment estimates satisfying `-2 < phi_i < 0`.
7. Compute the robust center as the median of admissible `phi_i`.
8. Scale dispersion as `1.4826 × MAD`.
9. Apply Huber-type weights with `c = 1.345`.
10. Aggregate admissible adjustment coefficients using their Huber-weighted mean.

## Requirements

- StataNow 19 (the current release was validated there)
- `xtdcce2`
- Stata built-in `rreg`

If `xtdcce2` is not installed, install it before using `rcsardl`.

## Installation from a local folder

Place `rcsardl.ado` and `rcsardl.sthlp` in a folder on Stata's adopath, or add the folder temporarily:

```stata
adopath ++ "C:\path\to\folder"
which rcsardl
help rcsardl
```

## Installation from GitHub

After this repository is published, installation will use the repository's raw-content path:

```stata
net install rcsardl, from("https://raw.githubusercontent.com/zehrayalnizz41-cmd/rcsardl/main/")
```

Replace `YOUR-GITHUB-USERNAME` with the actual GitHub username before publishing the installation instructions.

## Basic use

```stata
xtset id year
rcsardl lnco2 lnenergy lngdp, panel(id) time(year)
```

Useful stored results include:

```stata
ereturn list
matrix list e(longrun)
matrix list e(unit_adjustment)
```

## Empirical validation

The command was checked against the manuscript's balanced OECD application (38 countries, 1990–2020). Version 0.1.5 reproduces the reported point-estimation benchmark:

| Quantity | Validated result |
|---|---:|
| Conventional CS-ARDL adjustment | -0.6057704 |
| RCS-ARDL adjustment | -0.241517 |
| Admissible units | 34/38 |
| Median admissible phi | -0.194498 |
| MAD-based robust scale | 0.156503 |
| Huber cutoff | 0.210497 |
| Huber down-weighted units | 8 |
| Long-run energy | 0.8890433 |
| Long-run GDP | -0.0638206 |

## Bootstrap inference

Version 0.1.5 deliberately **does not implement bootstrap inference**.

The manuscript reports a 50-replication block bootstrap, but the exact resampling/block scheme is not present in the available replication `.do` files. To avoid silently substituting a different bootstrap design, this release contains only the point estimator that has been reproduced and benchmarked exactly.

## Versioning

- **0.1.5**: validated point estimator; named extraction of CS-ARDL coefficients; robust stage uses the full available panel conditional on the classical long-run vector.
- Bootstrap inference is reserved for a later release after the exact resampling procedure is documented and independently reproduced.

## Citation

Please cite the accompanying manuscript when using the estimator:

Yalnız, Z. (2026). *Robust Estimation of Equilibrium Adjustment in Cross-Sectionally Dependent Dynamic Panels: The RCS-ARDL Approach*. Manuscript submitted to *Empirical Economics*.

## License

MIT License. See `LICENSE`.
