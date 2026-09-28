RCS-ARDL for Stata
Version 0.4.0 — 28 September 2026
rcsardl implements the validated point-estimation algorithm for the Robust Cross-Sectionally Augmented Autoregressive Distributed Lag (RCS-ARDL) approach developed in:
Zehra Yalnız (2026), Robust Estimation of Equilibrium Adjustment in Cross-Sectionally Dependent Dynamic Panels: The RCS-ARDL Approach.

Version 0.4.0 preserves the Frozen-V3 point estimator and adds crlags() support and optional panel-synchronous block-bootstrap inference.
What the command does
The command preserves the conventional CS-ARDL long-run coefficient vector and robustifies only the equilibrium-adjustment component.
The implemented sequence is:
1. Estimate conventional CS-ARDL with xtdcce2.
2. Retain the conventional CS-ARDL long-run coefficients unchanged.
3. Construct the equilibrium error from those long-run coefficients.
4. Estimate a robust error-correction regression separately for each panel unit using rreg.
5. Use cross-sectional medians of Δy, Δx, and lagged ECT as robust common controls.
6. Retain dynamically admissible unit-specific adjustment estimates satisfying -2 < phi_i < 0.
7. Compute the robust center as the median of admissible phi_i.
8. Scale dispersion as 1.4826 × MAD.
9. Apply Huber-type weights with c = 1.345.
10. Aggregate admissible adjustment coefficients using their Huber-weighted mean.
Requirements
- Stata 17 or later
- xtdcce2
- Stata built-in rreg
The current release has been tested under StataNow 19.
If xtdcce2 is not installed, install it before using rcsardl.
Installation from a local folder
Place rcsardl.ado and rcsardl.sthlp in a folder on Stata's adopath, or add the folder temporarily:
adopath ++ "C:\path\to\folder"
which rcsardl
help rcsardl
Installation from GitHub
net install rcsardl, from("https://raw.githubusercontent.com/zehrayalnizz41-cmd/rcsardl/main/")
After installation:
which rcsardl
help rcsardl
Basic use
xtset id year
rcsardl lnco2 lnenergy lngdp, panel(id) time(year)
The number of lags of the cross-sectional averages can be specified with crlags():
rcsardl lnco2 lnenergy lngdp, panel(id) time(year) crlags(0)
or:
rcsardl lnco2 lnenergy lngdp, panel(id) time(year) crlags(1)
Useful stored results include:
ereturn list
matrix list e(longrun)
matrix list e(longrun_se)
matrix list e(unit_adjustment)
matrix list e(results)
Empirical validation
The Frozen-V3 point estimator was checked against the manuscript's balanced OECD application (38 countries, 1990–2020). The validated point-estimation benchmark is:
Quantity	Validated result
Conventional CS-ARDL adjustment	-0.6057704
RCS-ARDL adjustment	-0.241517
Admissible units	34/38
Median admissible phi	-0.194498
MAD-based robust scale	0.156503
Huber cutoff	0.210497
Huber down-weighted units	8
Long-run energy	0.8890433
Long-run GDP	-0.0638206


Version 0.4.0 preserves this Frozen-V3 point-estimation structure.
Bootstrap inference
Version 0.4.0 adds optional panel-synchronous block-bootstrap inference for the RCS-ARDL adjustment coefficient.
The supported bootstrap schemes are:
- sbb: stationary block bootstrap
- cbb: circular block bootstrap
- mbb: moving block bootstrap
- nbb: nonoverlapping block bootstrap
The same sampled time indices are applied simultaneously to all cross-sectional units in each replication.
This panel-synchronous design preserves the contemporaneous cross-sectional structure while resampling serially dependent time blocks.
Each successful bootstrap replication re-estimates the complete RCS-ARDL procedure.
Example:
rcsardl lnef lnurban lnenergy d_lngdppc, panel(i) time(year) crlags(1) bootstrap reps(499) type(sbb) lblock(auto) seed(12345)
The lblock(auto) option uses an automatic block-length procedure adapted from the Baum–Otero blockboot approach.
In the panel implementation, variable-specific automatic block lengths are calculated from cross-sectional median time series and combined into a common block length for synchronous resampling.
Bootstrap validation
The Version 0.4.0 SBB implementation was tested with:
rcsardl lnef lnurban lnenergy d_lngdppc, panel(i) time(year) crlags(1) bootstrap reps(499) type(sbb) lblock(auto) seed(12345)
The software test produced:
Quantity	Result
Requested bootstrap replications	499
Successful bootstrap replications	498
Common block length	4
RCS-ARDL adjustment	-0.3028
Bootstrap standard error	0.1130
95% percentile bootstrap CI	[-0.5313, -0.0559]


These values are reported only as a software-validation example.
Stored bootstrap results
When bootstrap inference is requested, additional stored results include:
ereturn list
matrix list e(boot_dist)
matrix list e(autoblock)
Important bootstrap scalars include:
e(boot_se)        Bootstrap standard error
e(boot_ll)        Lower percentile confidence limit
e(boot_ul)        Upper percentile confidence limit
e(boot_reps)      Requested bootstrap replications
e(boot_success)   Successful bootstrap replications
e(boot_lblock)    Common bootstrap block length
Versioning
- 0.4.0: adds crlags(), panel-synchronous block-bootstrap inference, SBB/CBB/MBB/NBB schemes, automatic block-length selection, bootstrap standard errors, percentile confidence intervals, and formatted estimation output.
- 0.1.5: validated Frozen-V3 point estimator; named extraction of CS-ARDL coefficients; robust stage uses the full available panel conditional on the conventional long-run vector.
Methodological attribution
The RCS-ARDL estimator was developed by Zehra Yalnız.
The block-resampling and automatic block-length logic in Version 0.4.0 is adapted from the Baum–Otero blockboot approach.
The RCS-ARDL panel implementation differs by applying common sampled time indices synchronously across all panel units and by re-estimating the complete RCS-ARDL procedure within each successful bootstrap replication.
Citation
Please cite the accompanying manuscript when using the estimator:
Yalnız, Z. (2026). Robust Estimation of Equilibrium Adjustment in Cross-Sectionally Dependent Dynamic Panels: The RCS-ARDL Approach. Manuscript submitted to Empirical Economics.
License
MIT License. See LICENSE.
