*******************************************************
* rcsardl_example.do
* Minimal usage example
*******************************************************

* Load your panel dataset first.
* Required: numeric panel ID, numeric time variable,
* dependent variable, and at least two explanatory variables.

xtset id year

rcsardl lnco2 lnenergy lngdp, panel(id) time(year)

ereturn list
matrix list e(longrun)
matrix list e(unit_adjustment)
