*******************************************************
* rcsardl_example.do
* RCS-ARDL Version 0.4.0
* Minimal usage and bootstrap examples
*******************************************************

* Load your panel dataset first.
* Required: numeric panel ID, numeric time variable,
* dependent variable, and at least two explanatory variables.

xtset id year

*******************************************************
* 1. Basic RCS-ARDL estimation
*******************************************************

rcsardl lnco2 lnenergy lngdp, panel(id) time(year)

*******************************************************
* 2. Specify cross-sectional-average lags
*******************************************************

rcsardl lnco2 lnenergy lngdp, panel(id) time(year) crlags(1)

*******************************************************
* 3. Panel-synchronous block-bootstrap inference
*******************************************************

rcsardl lnco2 lnenergy lngdp, panel(id) time(year) crlags(1) bootstrap reps(499) type(sbb) lblock(auto) seed(12345)

*******************************************************
* 4. Stored results
*******************************************************

ereturn list
matrix list e(longrun)
matrix list e(longrun_se)
matrix list e(unit_adjustment)
matrix list e(results)

* Bootstrap-specific stored results are available
* after estimation with the bootstrap option.

matrix list e(boot_dist)
matrix list e(autoblock)
