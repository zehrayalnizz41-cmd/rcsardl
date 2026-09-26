*! rcsardl.ado
*! Version 0.1.5  26sep2026
*! RCS-ARDL point estimator implementing the frozen V3 adjustment algorithm
*! Method: Zehra Yalnız (2026), "Robust Estimation of Equilibrium Adjustment in
*! Cross-Sectionally Dependent Dynamic Panels: The RCS-ARDL Approach"
*!
*! IMPORTANT:
*! - Long-run coefficients are taken unchanged from conventional CS-ARDL.
*! - Robustification applies only to the equilibrium-adjustment coefficient.
*! - Dynamic admissibility is fixed at -2 < phi_i < 0.
*! - MAD normal-consistency factor is fixed at 1.4826.
*! - Huber tuning constant is fixed at 1.345.
*! - Cross-sectional controls in the robust stage are medians.
*! - This version does NOT yet implement the empirical block bootstrap.

capture program drop rcsardl
program define rcsardl, eclass sortpreserve
    version 19.0

    syntax varlist(ts numeric min=3) [if] [in], PANEL(varname numeric) TIME(varname numeric)

    marksample touse

    gettoken depvar indepvars : varlist
    local K : word count `indepvars'

    capture which xtdcce2
    if _rc {
        di as error "rcsardl requires xtdcce2."
        di as error "Install xtdcce2 before running rcsardl."
        exit 499
    }

    quietly xtset `panel' `time'

    local lrlist "L.`depvar'"
    foreach x of local indepvars {
        local lrlist "`lrlist' `x' L.`x'"
    }

    capture noisily xtdcce2 `depvar' if `touse', ///
        lr(`lrlist') ///
        lr_options(ardl) ///
        cr(`depvar' `indepvars') cr_lags(1)

    if _rc {
        di as error "Conventional CS-ARDL estimation failed (xtdcce2 return code = " _rc ")."
        exit _rc
    }

    tempname B LR U
    matrix `B' = e(b)

    * ------------------------------------------------------------
    * Extract CS-ARDL coefficients BY NAME, not by matrix position.
    * This is safer across xtdcce2 versions.
    * ------------------------------------------------------------
    capture scalar __rcs_ect_cs = _b[lr_`depvar']
    if _rc {
        di as error "Could not identify adjustment coefficient lr_`depvar' in xtdcce2 results."
        exit 498
    }

    matrix `LR' = J(1,`K',.)
    local j = 0
    foreach x of local indepvars {
        local ++j
        capture scalar __rcs_theta = _b[lr_`x']
        if _rc {
            di as error "Could not identify long-run coefficient lr_`x' in xtdcce2 results."
            exit 498
        }
        matrix `LR'[1,`j'] = __rcs_theta
    }
    matrix colnames `LR' = `indepvars'

    * ------------------------------------------------------------
    * Frozen-V3 robust stage uses the FULL available panel after
    * conditioning on the classical CS-ARDL long-run coefficients.
    * Do NOT restrict ECT construction to xtdcce2 e(sample).
    * ------------------------------------------------------------
    quietly sort `panel' `time'

    tempvar ect_cs L_ect_cs d_y med_dy med_ect
    quietly gen double `ect_cs' = `depvar' if `touse'

    local j = 0
    foreach x of local indepvars {
        local ++j
        scalar __rcs_theta = `LR'[1,`j']
        quietly replace `ect_cs' = `ect_cs' - __rcs_theta*`x' if `touse'
    }

    quietly gen double `L_ect_cs' = L.`ect_cs'
    quietly gen double `d_y' = D.`depvar'

    * Generate ALL time-series transformed regressors while data are
    * still sorted by panel and time.  Do not use bysort until all L./D.
    * operators have been evaluated.
    local dxlist
    local meddxlist
    local j = 0
    foreach x of local indepvars {
        local ++j
        tempvar dx`j' meddx`j'
        quietly gen double `dx`j'' = D.`x'
        local dxlist "`dxlist' `dx`j''"
        local meddxlist "`meddxlist' `meddx`j''"
    }

    * Cross-sectional medians by time.
    quietly bysort `time': egen double `med_dy' = median(`d_y')
    quietly bysort `time': egen double `med_ect' = median(`L_ect_cs')

    local j = 0
    foreach x of local indepvars {
        local ++j
        quietly bysort `time': egen double `meddx`j'' = median(`dx`j'')
    }

    quietly levelsof `panel' if `touse', local(panelids)
    local NTOTAL : word count `panelids'

    tempfile units
    tempname unitpost
    postfile `unitpost' double uid phi using `units', replace

    foreach ii of local panelids {
        capture quietly rreg ///
            `d_y' `L_ect_cs' `dxlist' ///
            `med_dy' `meddxlist' `med_ect' ///
            if `panel'==`ii' & `touse'

        if _rc==0 {
            capture scalar __rcs_phi = _b[`L_ect_cs']
            if _rc==0 {
                if (__rcs_phi < 0 & __rcs_phi > -2) {
                    post `unitpost' (`ii') (__rcs_phi)
                }
            }
        }
    }
    postclose `unitpost'

    preserve
        quietly use `units', clear
        quietly drop if missing(phi)
        quietly count
        local NADM = r(N)

        if `NADM' < 1 {
            restore
            di as error "No dynamically admissible unit-specific adjustment coefficients were obtained."
            exit 498
        }

        quietly summarize phi, detail
        scalar __rcs_med_phi = r(p50)

        quietly gen double dev_phi = abs(phi-__rcs_med_phi)
        quietly summarize dev_phi, detail
        scalar __rcs_mad_raw = r(p50)

        if (__rcs_mad_raw<=0 | missing(__rcs_mad_raw)) {
            scalar __rcs_mad_raw = 0.000001
        }

        scalar __rcs_scale = 1.4826*__rcs_mad_raw
        scalar __rcs_cutoff = 1.345*__rcs_scale

        quietly gen double w_phi = 1
        quietly replace w_phi = __rcs_cutoff/dev_phi ///
            if dev_phi>__rcs_cutoff & dev_phi>0

        quietly count if w_phi < 1
        local NDOWN = r(N)

        quietly gen double wx_phi = w_phi*phi
        quietly summarize wx_phi, meanonly
        scalar __rcs_sum_wx = r(sum)
        quietly summarize w_phi, meanonly
        scalar __rcs_sum_w = r(sum)

        scalar __rcs_ect_rcs = __rcs_sum_wx/__rcs_sum_w

        mkmat uid phi w_phi, matrix(`U')
        matrix colnames `U' = panel_id phi_robust huber_weight
    restore

    ereturn clear
    ereturn scalar ect_cs = __rcs_ect_cs
    ereturn scalar ect_rcs = __rcs_ect_rcs
    ereturn scalar n_total = `NTOTAL'
    ereturn scalar n_admissible = `NADM'
    ereturn scalar n_excluded = `NTOTAL' - `NADM'
    ereturn scalar n_downweighted = `NDOWN'
    ereturn scalar median_phi = __rcs_med_phi
    ereturn scalar mad_raw = __rcs_mad_raw
    ereturn scalar mad_scale = __rcs_scale
    ereturn scalar huber_cutoff = __rcs_cutoff
    ereturn scalar huber_c = 1.345
    ereturn scalar admiss_lower = -2
    ereturn scalar admiss_upper = 0

    ereturn matrix longrun = `LR'
    ereturn matrix unit_adjustment = `U'

    ereturn local cmd "rcsardl"
    ereturn local depvar "`depvar'"
    ereturn local indepvars "`indepvars'"
    ereturn local panelvar "`panel'"
    ereturn local timevar "`time'"
    ereturn local method "Frozen V3 RCS-ARDL adjustment estimator; full robust-stage panel"

    di
    di as text "RCS-ARDL (Frozen V3 adjustment estimator)"
    di as text "Long-run coefficients retained from conventional CS-ARDL"
    di as text "{hline 68}"
    di as text "Conventional CS-ARDL adjustment coefficient" ///
        _col(52) as result %10.4f e(ect_cs)
    di as text "RCS-ARDL robust adjustment coefficient" ///
        _col(52) as result %10.4f e(ect_rcs)
    di as text "Admissible units (-2 < phi_i < 0)" ///
        _col(52) as result %6.0f e(n_admissible) "/" %6.0f e(n_total)
    di as text "Huber down-weighted units" ///
        _col(52) as result %10.0f e(n_downweighted)
    di as text "Median admissible phi" ///
        _col(52) as result %10.4f e(median_phi)
    di as text "MAD-based robust scale (1.4826 x MAD)" ///
        _col(52) as result %10.4f e(mad_scale)
    di as text "Huber cutoff (1.345 x scale)" ///
        _col(52) as result %10.4f e(huber_cutoff)
    di as text "{hline 68}"
    di as text "Long-run coefficient vector:"
    matrix list e(longrun), noheader format(%10.4f)
end
