{smcl}
{* *! version 0.1.5 26sep2026}{...}
{vieweralsosee "xtdcce2" "help xtdcce2"}{...}
{title:Title}

{phang}
{bf:rcsardl} {hline 2} Robust Cross-Sectionally Augmented ARDL adjustment estimator

{title:Syntax}

{p 8 17 2}
{cmd:rcsardl} {it:depvar indepvars} {cmd:,}
{opt panel(varname)}
{opt time(varname)}

{title:Description}

{pstd}
{cmd:rcsardl} implements the point-estimation algorithm of the RCS-ARDL
approach described by Zehra Yalnız (2026). The command first estimates a
conventional CS-ARDL model using {cmd:xtdcce2}. The conventional CS-ARDL
long-run coefficient vector is retained unchanged. Robustification is applied
only to the equilibrium-adjustment component.

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

{title:Options}

{phang}
{opt panel(varname)} specifies the numeric panel identifier.

{phang}
{opt time(varname)} specifies the numeric time variable.

{title:Dependencies}

{pstd}
{cmd:rcsardl} requires {cmd:xtdcce2}. The command also uses Stata's built-in
{cmd:rreg}. The current release was validated in StataNow 19.

{title:Stored results}

{pstd}
{cmd:rcsardl} stores the following in {cmd:e()}:

{synoptset 28 tabbed}{...}
{synopt:{cmd:e(ect_cs)}}conventional CS-ARDL adjustment coefficient{p_end}
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
{synopt:{cmd:e(unit_adjustment)}}matrix of admissible panel IDs, robust phi_i, and Huber weights{p_end}

{title:Example}

{phang2}{cmd:. xtset id year}

{phang2}{cmd:. rcsardl lnco2 lnenergy lngdp, panel(id) time(year)}

{phang2}{cmd:. ereturn list}

{phang2}{cmd:. matrix list e(unit_adjustment)}

{title:Validated OECD benchmark}

{pstd}
Using the balanced 38-country OECD panel for 1990-2020 used in the manuscript,
the validated implementation reproduces the following point-estimation results:

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

{title:Important limitation of version 0.1.5}

{pstd}
Version 0.1.5 implements the validated point estimator only. It does not
implement the manuscript's 50-replication block-bootstrap inference because
the exact resampling/block scheme is not documented in the available
replication code. No alternative bootstrap scheme is substituted silently.

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
