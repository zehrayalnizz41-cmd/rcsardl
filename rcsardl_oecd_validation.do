*******************************************************
* rcsardl_oecd_validation.do
* RCS-ARDL Version 0.4.0
* Regression test for the manuscript's validated OECD
* Frozen-V3 point-estimation benchmark.
*
* Requires the validated OECD dataset already loaded
* in memory with variables:
* id year lnco2 lnenergy lngdp
*
* Note:
* The original validated benchmark used one lag of the
* cross-sectional averages. Version 0.4.0 therefore uses
* crlags(1) explicitly in this regression test.
*
* This file checks the point estimator only.
* Bootstrap inference is not part of this OECD benchmark.
*******************************************************

xtset id year

rcsardl lnco2 lnenergy lngdp, panel(id) time(year) crlags(1)

display "Checking validated OECD point-estimation benchmark..."

assert abs(e(ect_cs) - (-0.6057704)) < 0.000001
assert abs(e(ect_rcs) - (-0.241517)) < 0.00001
assert e(n_total) == 38
assert e(n_admissible) == 34
assert e(n_downweighted) == 8
assert abs(e(median_phi) - (-0.194498)) < 0.00001
assert abs(e(mad_scale) - 0.156503) < 0.00001
assert abs(e(huber_cutoff) - 0.210497) < 0.00001

matrix LR = e(longrun)
assert abs(LR[1,1] - 0.8890433) < 0.000001
assert abs(LR[1,2] - (-0.0638206)) < 0.000001

display "OECD Frozen-V3 point-estimation benchmark reproduced successfully."
