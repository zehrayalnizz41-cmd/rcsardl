{smcl}
{* *! version 0.4.0 28sep2026}{...}
{vieweralsosee "xtdcce2" "help xtdcce2"}{...}
{title:Title}

{phang}
{bf:rcsardl} {hline 2} Robust Cross-Sectionally Augmented ARDL estimator

{title:Syntax}

{p 8 17 2}
{cmd:rcsardl} {it:depvar indepvars} {cmd:,}
{opt panel(varname)}
{opt time(varname)}
[{opt crlags(#)}
{opt bootstrap}
{opt reps(#)}
{opt type(string)}
{opt lblock(string)}
{opt seed(#)}
{opt level(#)}]

{title:Description}

{pstd}
{cmd:rcsardl} implements the Robust Cross-Sectionally Augmented
Autoregressive Distributed Lag (RCS-ARDL) approach described by
Zehra Yalnız (2026).

{pstd}
The command first estimates a conventional CS-ARDL model using
{cmd:xtdcce2}. The conventional CS-ARDL long-run coefficient vector is
retained unchanged. Robustification is applied only to the
equilibrium-adjustment component.

{pstd}
Conditional on the conventional CS-ARDL long-run coefficients, the command
constructs the lagged equilibrium error and estimates a robust error-correction
regression separately for each panel unit using Stata's {cmd:rreg}. The robust
stage includes cross-sectional medians of the dependent-variable change,
explanatory-variable changes, and lagged equilibrium error.

{pstd}
Unit-specific robust adjustment coefficients are retained only when

{p 12 12 2}
-2 < phi_i < 0.

{pstd}
The retained coefficients are aggregated using a cross-sectional median,
MAD scaling by 1.4826, and Huber-type weights with tuning constant 1.345.
The resulting Huber-weighted mean is reported as the RCS-ARDL adjustment
coefficient.

{pstd}
Version 0.4.0 preserves the validated Frozen-V3 point-estimation structure
and adds {cmd:crlags()} support and optional panel-synchronous block-bootstrap
inference.

{title:Options}

{phang}
{opt panel(varname)} specifies the numeric panel identifier.

{phang}
{opt time(varname)} specifies the numeric time variable.

{phang}
{opt crlags(#)} specifies the number of lags of the cross-sectional averages
used in the conventional CS-ARDL stage. The specified value is passed to
{cmd:xtdcce2} through {cmd:cr_lags()} and is preserved in bootstrap
re-estimation.

{phang}
{opt bootstrap} requests panel-synchronous block-bootstrap inference for the
RCS-ARDL adjustment coefficient.

{phang}
{opt reps(#)} specifies the requested number of bootstrap replications.
The default is 499.

{phang}
{opt type(string)} specifies the block-bootstrap scheme. Supported types are
{cmd:sbb}, {cmd:cbb}, {cmd:mbb}, and {cmd:nbb}.

{phang}
{opt lblock(string)} specifies the block length. Use {cmd:lblock(auto)} for
automatic block-length selection or supply a positive integer.

{phang}
{opt seed(#)} sets the random-number seed for bootstrap resampling.
The default is -1, which leaves the current random-number state unchanged.

{phang}
{opt level(#)} specifies the confidence level used for reported confidence
intervals. The default is 95.

{title:Bootstrap inference}

{pstd}
The bootstrap is panel-synchronous: within each replication, the same sampled
time indices are applied simultaneously to all cross-sectional units. This
preserves the contemporaneous cross-sectional structure while resampling
serially dependent time blocks.

{pstd}
Each successful bootstrap replication re-estimates the complete RCS-ARDL
procedure, including conventional CS-ARDL estimation, long-run coefficient
extraction, equilibrium-error construction, unit-specific robust
error-correction regressions, dynamic-admissibility screening, MAD scaling,
Huber weighting, and RCS-ARDL aggregation.

{pstd}
Supported block-bootstrap schemes are:

{p 8 12 2}
{cmd:sbb}: stationary block bootstrap{break}
{cmd:cbb}: circular block bootstrap{break}
{cmd:mbb}: moving block bootstrap{break}
{cmd:nbb}: nonoverlapping block bootstrap

{pstd}
With {cmd:lblock(auto)}, variable-specific automatic block lengths are
calculated from cross-sectional median time series and combined into a common
block length for synchronous panel resampling. The block-resampling and
automatic block-length logic are adapted from the Baum-Otero {cmd:blockboot}
approach.

{pstd}
Bootstrap inference requires a balanced common-time panel over the estimation
sample.

{title:Dependencies}

{pstd}
{cmd:rcsardl} requires {cmd:xtdcce2}. The command also uses Stata's built-in
{cmd:rreg}. The ado file is written for Stata 17 or later. The current release
has been tested under StataNow 19.

{title:Stored results}

{pstd}
{cmd:rcsardl} stores the following in {cmd:e()}:

{synoptset 30 tabbed}{...}
{synopt:{cmd:e(ect_cs)}}conventional CS-ARDL adjustment coefficient{p_end}
{synopt:{cmd:e(ect_cs_se)}}standard error of conventional adjustment coefficient{p_end}
{synopt:{cmd:e(ect_rcs)}}RCS-ARDL robust adjustment coefficient{p_end}
{synopt:{cmd:e(n_total)}}number of panel units considered{p_end}
{synopt:{cmd:e(n_admissible)}}number satisfying -2 < phi_i < 0{p_end}
{synopt:{cmd:e(n_excluded)}}number excluded by dynamic admissibility{p_end}
{synopt:{cmd:e(n_downweighted)}}number receiving Huber weight below one{p_end}
{synopt:{cmd:e(median_phi)}}median admissible robust adjustment coefficient{p_end}
{synopt:{cmd:e(mad_raw)}}unscaled median absolute deviation{p_end}
{synopt:{cmd:e(mad_scale)}}1.4826 times MAD{p_end}
{synopt:{cmd:e(huber_cutoff)}}1.345 times the MAD-based scale{p_end}
{synopt:{cmd:e(huber_c)}}Huber tuning constant (1.345){p_end}
{synopt:{cmd:e(admiss_lower)}}lower admissibility bound (-2){p_end}
{synopt:{cmd:e(admiss_upper)}}upper admissibility bound (0){p_end}
{synopt:{cmd:e(longrun)}}row vector of conventional CS-ARDL long-run coefficients{p_end}
{synopt:{cmd:e(longrun_se)}}row vector of conventional long-run standard errors{p_end}
{synopt:{cmd:e(unit_adjustment)}}matrix of admissible panel IDs, robust phi_i, and Huber weights{p_end}
{synopt:{cmd:e(results)}}formatted estimation-results matrix{p_end}

{pstd}
When bootstrap inference is requested, additional stored results include:

{synoptset 30 tabbed}{...}
{synopt:{cmd:e(boot_se)}}bootstrap standard error of the RCS-ARDL adjustment coefficient{p_end}
{synopt:{cmd:e(boot_ll)}}lower percentile bootstrap confidence limit{p_end}
{synopt:{cmd:e(boot_ul)}}upper percentile bootstrap confidence limit{p_end}
{synopt:{cmd:e(boot_reps)}}requested number of bootstrap replications{p_end}
{synopt:{cmd:e(boot_success)}}number of successful bootstrap replications{p_end}
{synopt:{cmd:e(boot_lblock)}}common bootstrap block length{p_end}
{synopt:{cmd:e(boot_dist)}}bootstrap distribution of the RCS-ARDL adjustment coefficient{p_end}
{synopt:{cmd:e(autoblock)}}variable-specific automatic block lengths when {cmd:lblock(auto)} is used{p_end}

{title:Examples}

{phang2}{cmd:. xtset id year}

{phang2}{cmd:. rcsardl lnco2 lnenergy lngdp, panel(id) time(year)}

{phang2}{cmd:. rcsardl lnco2 lnenergy lngdp, panel(id) time(year) crlags(1)}

{phang2}{cmd:. rcsardl lnef lnurban lnenergy d_lngdppc, panel(i) time(year) crlags(1) bootstrap reps(499) type(sbb) lblock(auto) seed(12345)}

{phang2}{cmd:. ereturn list}

{phang2}{cmd:. matrix list e(unit_adjustment)}

{title:Validated OECD point-estimation benchmark}

{pstd}
Using the balanced 38-country OECD panel for 1990-2020 used in the manuscript,
the validated Frozen-V3 point estimator reproduces the following results:

{p 8 12 2}
Conventional CS-ARDL adjustment: -0.6057704{break}
RCS-ARDL adjustment: -0.241517 (rounded){break}
Admissible units: 34/38{break}
Median admissible phi: -0.194498{break}
MAD-based scale: 0.156503{break}
Huber cutoff: 0.210497{break}
Huber down-weighted units: 8{break}
Long-run energy coefficient: 0.8890433{break}
Long-run GDP coefficient: -0.0638206

{title:Bootstrap software test}

{pstd}
The Version 0.4.0 stationary block-bootstrap implementation was tested with
499 requested replications, {cmd:crlags(1)}, {cmd:lblock(auto)}, and seed
12345. The software test produced 498 successful replications, a common block
length of 4, an RCS-ARDL adjustment coefficient of -0.3028, a bootstrap
standard error of 0.1130, and a 95 percent percentile confidence interval of
[-0.5313, -0.0559].

{pstd}
These values are reported as a software-validation example and are not
intended as general empirical findings.

{title:Methodological attribution}

{pstd}
The RCS-ARDL estimator was developed by Zehra Yalnız.

{pstd}
The block-resampling and automatic block-length logic in Version 0.4.0 is
adapted from the Baum-Otero {cmd:blockboot} approach. The RCS-ARDL panel
implementation applies common sampled time indices synchronously across all
panel units and re-estimates the complete RCS-ARDL procedure within each
successful bootstrap replication.

{title:Reference}

{pstd}
Yalnız, Z. (2026). Robust Estimation of Equilibrium Adjustment in
Cross-Sectionally Dependent Dynamic Panels: The RCS-ARDL Approach.
Manuscript submitted to {it:Empirical Economics}.

{title:Author}

{pstd}
Zehra Yalnız

{title:License}

{pstd}
MIT License. See the LICENSE file distributed with the package.
