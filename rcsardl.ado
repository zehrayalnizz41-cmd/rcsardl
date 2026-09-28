*! rcsardl.ado
*! Version 0.4.0  28sep2026
*! RCS-ARDL with panel-synchronous block-bootstrap inference
*! Core estimator: Zehra Yalnız (2026)
*! Bootstrap design adapted from the block-bootstrap logic of Baum and Otero.
*!
*! Core estimator is the frozen V3 algorithm:
*! - CS-ARDL long-run coefficients are retained unchanged.
*! - Robustification applies only to the equilibrium-adjustment coefficient.
*! - Dynamic admissibility: -2 < phi_i < 0.
*! - MAD consistency factor: 1.4826.
*! - Huber tuning constant: 1.345.
*! - Cross-sectional controls in the robust stage are medians.

capture program drop rcsardl
capture program drop _rcsardl_point

program define _rcsardl_point, rclass
    version 17.0

    syntax varlist(ts numeric min=3) [if] [in], ///
        PANEL(varname numeric) TIME(varname numeric) ///
        [CRLAGS(integer 0)]

    marksample touse

    if (`crlags' < 0) {
        exit 198
    }

    gettoken depvar indepvars : varlist
    local K : word count `indepvars'

    capture which xtdcce2
    if _rc {
        exit 499
    }

    quietly xtset `panel' `time'

    local lrlist "L.`depvar'"
    foreach x of local indepvars {
        local lrlist "`lrlist' `x' L.`x'"
    }

    capture quietly xtdcce2 `depvar' if `touse', ///
        lr(`lrlist') ///
        lr_options(ardl) ///
        cr(`depvar' `indepvars') ///
        cr_lags(`crlags')

    if _rc {
        exit _rc
    }

    tempname LR LRSE U

    capture scalar __rcs_ect_cs = _b[lr_`depvar']
    if _rc {
        exit 498
    }

    capture scalar __rcs_ect_cs_se = _se[lr_`depvar']
    if _rc {
        scalar __rcs_ect_cs_se = .
    }

    matrix `LR' = J(1,`K',.)
    matrix `LRSE' = J(1,`K',.)

    local j = 0
    foreach x of local indepvars {
        local ++j

        capture scalar __rcs_theta = _b[lr_`x']
        if _rc {
            exit 498
        }

        matrix `LR'[1,`j'] = __rcs_theta

        capture scalar __rcs_theta_se = _se[lr_`x']
        if _rc {
            scalar __rcs_theta_se = .
        }
        matrix `LRSE'[1,`j'] = __rcs_theta_se
    }

    matrix colnames `LR' = `indepvars'
    matrix colnames `LRSE' = `indepvars'

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
    postfile `unitpost' double uid phi using `units'

    foreach ii of local panelids {
        capture quietly rreg ///
            `d_y' `L_ect_cs' `dxlist' ///
            `med_dy' `meddxlist' `med_ect' ///
            if `panel'==`ii' & `touse'

        if (_rc==0) {
            capture scalar __rcs_phi = _b[`L_ect_cs']
            if (_rc==0) {
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

        if (`NADM' < 1) {
            restore
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

    return scalar ect_cs = __rcs_ect_cs
    return scalar ect_cs_se = __rcs_ect_cs_se
    return scalar ect_rcs = __rcs_ect_rcs
    return scalar n_total = `NTOTAL'
    return scalar n_admissible = `NADM'
    return scalar n_excluded = `NTOTAL' - `NADM'
    return scalar n_downweighted = `NDOWN'
    return scalar median_phi = __rcs_med_phi
    return scalar mad_raw = __rcs_mad_raw
    return scalar mad_scale = __rcs_scale
    return scalar huber_cutoff = __rcs_cutoff

    return matrix longrun = `LR'
    return matrix longrun_se = `LRSE'
    return matrix unit_adjustment = `U'
end


program define rcsardl, eclass
    version 17.0

    syntax varlist(ts numeric min=3) [if] [in], ///
        PANEL(varname numeric) TIME(varname numeric) ///
        [ CRLAGS(integer 0) BOOTstrap REPS(integer 499) ///
          TYPE(string) LBLOCK(string) SEED(integer -1) LEVEL(real 95) ]

    marksample touse

    if (`crlags' < 0) {
        di as error "crlags() must be a nonnegative integer."
        exit 198
    }

    if (`reps' < 1) {
        di as error "reps() must be positive."
        exit 198
    }

    if (`level' <= 0 | `level' >= 100) {
        di as error "level() must be between 0 and 100."
        exit 198
    }

    if ("`type'"=="") local type "sbb"
    local type = lower("`type'")

    if !inlist("`type'","sbb","cbb","mbb","nbb") {
        di as error "type() must be sbb, cbb, mbb, or nbb."
        exit 198
    }

    if ("`lblock'"=="") local lblock "auto"

    gettoken depvar indepvars : varlist
    local K : word count `indepvars'

    capture quietly _rcsardl_point `varlist' if `touse', ///
        panel(`panel') time(`time') crlags(`crlags')

    if _rc {
        di as error "RCS-ARDL point estimation failed (return code = " _rc ")."
        exit _rc
    }

    tempname LR LRSE U ABLK BOOTDIST RESULTS
    matrix `LR' = r(longrun)
    matrix `LRSE' = r(longrun_se)
    matrix `U' = r(unit_adjustment)

    scalar __main_ect_cs = r(ect_cs)
    scalar __main_ect_cs_se = r(ect_cs_se)
    scalar __main_ect_rcs = r(ect_rcs)
    scalar __main_ntotal = r(n_total)
    scalar __main_nadm = r(n_admissible)
    scalar __main_nexcl = r(n_excluded)
    scalar __main_ndown = r(n_downweighted)
    scalar __main_medphi = r(median_phi)
    scalar __main_madraw = r(mad_raw)
    scalar __main_scale = r(mad_scale)
    scalar __main_cutoff = r(huber_cutoff)

    scalar __boot_se = .
    scalar __boot_ll = .
    scalar __boot_ul = .
    scalar __boot_lblock = .
    scalar __boot_success = 0

    if ("`bootstrap'"!="") {

        if (`seed'!=-1) quietly set seed `seed'

        tempfile original base bootresults

        quietly save `original'

        quietly keep if `touse'

        capture isid `panel' `time'
        if _rc {
            quietly use `original', clear
            di as error "Bootstrap requires unique panel-time observations."
            exit 459
        }

        tempvar Ti Nt srcidx

        quietly bysort `panel': gen long `Ti' = _N
        quietly summarize `Ti', meanonly
        local Tmin = r(min)
        local Tmax = r(max)

        if (`Tmin' != `Tmax') {
            quietly use `original', clear
            di as error "Panel-synchronous bootstrap currently requires a balanced panel."
            exit 459
        }

        quietly levelsof `panel', local(panelids)
        local Npanel : word count `panelids'

        quietly bysort `time': gen long `Nt' = _N
        quietly summarize `Nt', meanonly

        if (r(min)!=`Npanel' | r(max)!=`Npanel') {
            quietly use `original', clear
            di as error "All panel units must share the same time periods."
            exit 459
        }

        local T = `Tmax'

        if (`T' < 4) {
            quietly use `original', clear
            di as error "Too few time observations for block bootstrap."
            exit 2001
        }

        quietly egen long `srcidx' = group(`time')
        quietly keep `panel' `depvar' `indepvars' `srcidx'
        quietly save `base'

        if (lower("`lblock'")=="auto") {

            quietly collapse (median) `depvar' `indepvars', by(`srcidx')
            quietly sort `srcidx'

            local nx = `K' + 1
            matrix `ABLK' = J(`nx',1,.)
            local rows
            local sumL = 0
            local j = 0

            foreach x of varlist `depvar' `indepvars' {
                local ++j
                mata: st_numscalar("__rcs_Lone", rcs_bstar("`x'","`type'"))
                local Lone = scalar(__rcs_Lone)
                matrix `ABLK'[`j',1] = `Lone'
                local sumL = `sumL' + `Lone'
                local rows "`rows' `x'"
            }

            matrix rownames `ABLK' = `rows'
            matrix colnames `ABLK' = Auto_Block

            local L = floor(`sumL'/`nx')
            if (`L' < 1) local L = 1
            if (`L' > `T') local L = `T'

            quietly use `base', clear
        }
        else {
            local L = real("`lblock'")

            if missing(`L') {
                quietly use `original', clear
                di as error "lblock() must be auto or a positive integer."
                exit 198
            }

            if (`L' < 1 | floor(`L') != `L') {
                quietly use `original', clear
                di as error "lblock() must be auto or a positive integer."
                exit 198
            }

            if (`L' > `T') {
                quietly use `original', clear
                di as error "Block length cannot exceed T."
                exit 198
            }
        }

        scalar __boot_lblock = `L'

        tempname bootpost
        postfile `bootpost' double phi_rcs using `bootresults'

        forvalues b = 1/`reps' {

            quietly clear
            quietly set obs `T'

            tempvar boottime bootsrc
            quietly gen long `boottime' = _n

            mata: rcs_make_indices(`T',`L',"`type'")
            quietly getmata (`bootsrc') = rcs_boot_indices
            quietly rename `bootsrc' `srcidx'

            capture quietly joinby `srcidx' using `base'

            if (_rc==0) {
                quietly sort `panel' `boottime'

                capture quietly _rcsardl_point `depvar' `indepvars', ///
                    panel(`panel') time(`boottime') crlags(`crlags')

                if (_rc==0) {
                    capture scalar __bp = r(ect_rcs)
                    if (_rc==0 & !missing(__bp)) {
                        post `bootpost' (__bp)
                    }
                }
            }
        }

        postclose `bootpost'

        quietly use `bootresults', clear
        quietly drop if missing(phi_rcs)
        quietly count
        local Bsuccess = r(N)
        scalar __boot_success = `Bsuccess'

        if (`Bsuccess' < 2) {
            quietly use `original', clear
            di as error "Fewer than two successful bootstrap replications were obtained."
            exit 498
        }

        quietly summarize phi_rcs
        scalar __boot_se = r(sd)

        local alpha = (100-`level')/2
        local upper = 100-`alpha'

        quietly centile phi_rcs, centile(`alpha' `upper')
        scalar __boot_ll = r(c_1)
        scalar __boot_ul = r(c_2)

        mkmat phi_rcs, matrix(`BOOTDIST')
        matrix colnames `BOOTDIST' = ect_rcs_boot

        quietly use `original', clear
    }

    scalar __zcrit = invnormal(1-(100-`level')/200)

    local NRES = `K' + 2
    matrix `RESULTS' = J(`NRES',6,.)
    matrix colnames `RESULTS' = Coef StdErr z Pvalue CI_L CI_U

    local rnames
    local j = 0

    foreach x of local indepvars {
        local ++j
        scalar __b = `LR'[1,`j']
        scalar __se = `LRSE'[1,`j']
        matrix `RESULTS'[`j',1] = __b
        matrix `RESULTS'[`j',2] = __se

        if (!missing(__se) & __se>0) {
            scalar __z = __b/__se
            scalar __p = 2*normal(-abs(__z))
            scalar __ll = __b-__zcrit*__se
            scalar __ul = __b+__zcrit*__se
            matrix `RESULTS'[`j',3] = __z
            matrix `RESULTS'[`j',4] = __p
            matrix `RESULTS'[`j',5] = __ll
            matrix `RESULTS'[`j',6] = __ul
        }

        local rn = strtoname("LR_`x'")
        local rnames "`rnames' `rn'"
    }

    local rowcs = `K' + 1
    matrix `RESULTS'[`rowcs',1] = __main_ect_cs
    matrix `RESULTS'[`rowcs',2] = __main_ect_cs_se

    if (!missing(__main_ect_cs_se) & __main_ect_cs_se>0) {
        scalar __z = __main_ect_cs/__main_ect_cs_se
        scalar __p = 2*normal(-abs(__z))
        scalar __ll = __main_ect_cs-__zcrit*__main_ect_cs_se
        scalar __ul = __main_ect_cs+__zcrit*__main_ect_cs_se
        matrix `RESULTS'[`rowcs',3] = __z
        matrix `RESULTS'[`rowcs',4] = __p
        matrix `RESULTS'[`rowcs',5] = __ll
        matrix `RESULTS'[`rowcs',6] = __ul
    }

    local rowrcs = `K' + 2
    matrix `RESULTS'[`rowrcs',1] = __main_ect_rcs

    if ("`bootstrap'"!="") {
        matrix `RESULTS'[`rowrcs',2] = __boot_se
        if (!missing(__boot_se) & __boot_se>0) {
            scalar __z = __main_ect_rcs/__boot_se
            scalar __p = 2*normal(-abs(__z))
            matrix `RESULTS'[`rowrcs',3] = __z
            matrix `RESULTS'[`rowrcs',4] = __p
        }
        matrix `RESULTS'[`rowrcs',5] = __boot_ll
        matrix `RESULTS'[`rowrcs',6] = __boot_ul
    }

    local rnames "`rnames' CS_ARDL_ECT RCS_ARDL_ECT"
    matrix rownames `RESULTS' = `rnames'

    di
    di as text "RCS-ARDL estimation results"
    di as text "Dependent variable: " as result "`depvar'"
    di as text "Cross-sectional average lags: " as result %4.0f `crlags'
    di as text "Admissible units: " as result %5.0f __main_nadm ///
        as text " / " as result %5.0f __main_ntotal
    di as text "Huber down-weighted units: " as result %5.0f __main_ndown

    if ("`bootstrap'"!="") {
        di as text "Bootstrap: " as result "`type'" ///
            as text ", common block length = " as result %4.0f __boot_lblock ///
            as text ", successful reps = " as result %5.0f __boot_success
    }

    di
    di as text "{hline 89}"
    di as text "Parameter" ///
        _col(24) "Coef." ///
        _col(36) "Std. Err." ///
        _col(49) "z" ///
        _col(59) "P>|z|" ///
        _col(70) "[" `level' "% Conf. Interval]"
    di as text "{hline 89}"

    local j = 0
    foreach x of local indepvars {
        local ++j

        scalar __db = `LR'[1,`j']
        scalar __dse = `LRSE'[1,`j']
        scalar __dz = .
        scalar __dp = .
        scalar __dll = .
        scalar __dul = .

        if (!missing(__dse) & __dse>0) {
            scalar __dz = __db/__dse
            scalar __dp = 2*normal(-abs(__dz))
            scalar __dll = __db-__zcrit*__dse
            scalar __dul = __db+__zcrit*__dse
        }

        di as text "LR: `x'" ///
            _col(23) as result %10.4f __db ///
            _col(35) %10.4f __dse ///
            _col(47) %9.2f __dz ///
            _col(58) %8.4f __dp ///
            _col(69) %9.4f __dll ///
            _col(79) %9.4f __dul
    }

    di as text "{hline 89}"

    scalar __db = __main_ect_cs
    scalar __dse = __main_ect_cs_se
    scalar __dz = .
    scalar __dp = .
    scalar __dll = .
    scalar __dul = .

    if (!missing(__dse) & __dse>0) {
        scalar __dz = __db/__dse
        scalar __dp = 2*normal(-abs(__dz))
        scalar __dll = __db-__zcrit*__dse
        scalar __dul = __db+__zcrit*__dse
    }

    di as text "CS-ARDL ECT" ///
        _col(23) as result %10.4f __db ///
        _col(35) %10.4f __dse ///
        _col(47) %9.2f __dz ///
        _col(58) %8.4f __dp ///
        _col(69) %9.4f __dll ///
        _col(79) %9.4f __dul

    scalar __db = __main_ect_rcs
    scalar __dse = .
    scalar __dz = .
    scalar __dp = .
    scalar __dll = .
    scalar __dul = .

    if ("`bootstrap'"!="") {
        scalar __dse = __boot_se
        scalar __dll = __boot_ll
        scalar __dul = __boot_ul

        if (!missing(__dse) & __dse>0) {
            scalar __dz = __db/__dse
            scalar __dp = 2*normal(-abs(__dz))
        }
    }

    di as text "RCS-ARDL ECT" ///
        _col(23) as result %10.4f __db ///
        _col(35) %10.4f __dse ///
        _col(47) %9.2f __dz ///
        _col(58) %8.4f __dp ///
        _col(69) %9.4f __dll ///
        _col(79) %9.4f __dul

    di as text "{hline 89}"
    di as text "Median admissible phi = " as result %9.4f __main_medphi ///
        as text "   MAD scale = " as result %9.4f __main_scale ///
        as text "   Huber cutoff = " as result %9.4f __main_cutoff

    if ("`bootstrap'"!="") {
        di as text "RCS-ARDL ECT: bootstrap SE; percentile bootstrap confidence interval."
    }

    ereturn clear

    ereturn scalar ect_cs = __main_ect_cs
    ereturn scalar ect_cs_se = __main_ect_cs_se
    ereturn scalar ect_rcs = __main_ect_rcs
    ereturn scalar cr_lags = `crlags'
    ereturn scalar n_total = __main_ntotal
    ereturn scalar n_admissible = __main_nadm
    ereturn scalar n_excluded = __main_nexcl
    ereturn scalar n_downweighted = __main_ndown
    ereturn scalar median_phi = __main_medphi
    ereturn scalar mad_raw = __main_madraw
    ereturn scalar mad_scale = __main_scale
    ereturn scalar huber_cutoff = __main_cutoff
    ereturn scalar huber_c = 1.345
    ereturn scalar admiss_lower = -2
    ereturn scalar admiss_upper = 0
    ereturn scalar level = `level'

    ereturn matrix longrun = `LR'
    ereturn matrix longrun_se = `LRSE'
    ereturn matrix unit_adjustment = `U'
    ereturn matrix results = `RESULTS'

    if ("`bootstrap'"!="") {
        ereturn scalar boot_se = __boot_se
        ereturn scalar boot_ll = __boot_ll
        ereturn scalar boot_ul = __boot_ul
        ereturn scalar boot_reps = `reps'
        ereturn scalar boot_success = __boot_success
        ereturn scalar boot_lblock = __boot_lblock
        ereturn matrix boot_dist = `BOOTDIST'
        if (lower("`lblock'")=="auto") {
            ereturn matrix autoblock = `ABLK'
        }
        ereturn local bootstrap "yes"
        ereturn local boot_type "`type'"
        ereturn local boot_lblock_method "`lblock'"
    }
    else {
        ereturn local bootstrap "no"
    }

    ereturn local cmd "rcsardl"
    ereturn local depvar "`depvar'"
    ereturn local indepvars "`indepvars'"
    ereturn local panelvar "`panel'"
    ereturn local timevar "`time'"
    ereturn local method "Frozen V3 RCS-ARDL; panel-synchronous block-bootstrap inference"
end


mata:

real scalar rcs_bstar(string scalar vname, string scalar btype)
{
    real colvector x, rhok, insignif, acov, kk, lam, x1, x2
    real scalar T, Kn, mmax, Bmax, c, rhokcrit
    real scalar mean_x, var_x, k, j, first_run, mhat, M
    real scalar Ghat, DCBhat, DSBhat, BstarSB, BstarCB, T1

    x = st_data(.,vname)
    x = select(x,x:<.)
    T = rows(x)

    if (T<4) return(1)

    Kn = max((5,ceil(log10(T))))
    mmax = ceil(sqrt(T))+Kn
    if (mmax>T-1) mmax = T-1

    Bmax = ceil(min((3*sqrt(T),T/3)))
    if (Bmax<1) Bmax = 1

    c = invnormal(0.975)
    rhokcrit = c*sqrt(log10(T)/T)

    mean_x = mean(x)
    var_x = mean((x:-mean_x):^2)

    if (var_x<=0 | missing(var_x)) return(1)

    rhok = J(mmax,1,.)
    insignif = J(mmax,1,.)

    for (k=1; k<=mmax; k++) {
        x1 = x[1..(T-k)]
        T1 = rows(x1)
        x2 = x[(1+k)..T]
        rhok[k] = (mean((x1:-mean_x):*(x2:-mean_x))/var_x)*(T1/T)
        insignif[k] = abs(rhok[k]) < rhokcrit
    }

    first_run = 0
    if (mmax>=Kn) {
        for (k=1; k<=mmax-Kn+1; k++) {
            if (sum(insignif[k..(k+Kn-1)])==Kn) {
                first_run = k
                break
            }
        }
    }

    if (first_run>0) mhat = first_run
    else mhat = mmax

    M = min((2*mhat,mmax))
    if (M<1) M = 1

    acov = J(2*M+1,1,.)
    kk = (-M..M)'

    for (j=1; j<=rows(kk); j++) {
        k = kk[j]

        if (k<0) {
            x1 = x[(1-k)..T]
            x2 = x[1..(T+k)]
            T1 = rows(x2)
        }
        else if (k>0) {
            x1 = x[1..(T-k)]
            x2 = x[(1+k)..T]
            T1 = rows(x1)
        }
        else {
            x1 = x
            x2 = x
            T1 = T
        }

        acov[j] = mean((x1:-mean_x):*(x2:-mean_x))*(T1/T)
    }

    lam = (abs(kk:/M):<0.5) :+ ///
          2:*(1:-abs(kk:/M)):*(abs(kk:/M):>=0.5):*(abs(kk:/M):<=1)

    Ghat = sum(lam:*abs(kk):*acov)
    DCBhat = (4/3)*sum(lam:*acov)^2
    DSBhat = 2*sum(lam:*acov)^2

    if (DSBhat<=0 | missing(DSBhat) | missing(Ghat)) {
        BstarSB = round(T^(1/3))
    }
    else {
        BstarSB = ((2*Ghat^2)/DSBhat)^(1/3)*T^(1/3)
    }

    if (DCBhat<=0 | missing(DCBhat) | missing(Ghat)) {
        BstarCB = ceil(T^(1/3))
    }
    else {
        BstarCB = ((2*Ghat^2)/DCBhat)^(1/3)*T^(1/3)
    }

    BstarSB = (BstarSB>Bmax ? Bmax : (BstarSB<1 ? 1 : round(BstarSB)))
    BstarCB = (BstarCB>Bmax ? Bmax : (BstarCB<1 ? 1 : ceil(BstarCB)))

    if (btype=="sbb") return(BstarSB)
    return(BstarCB)
}


void rcs_make_indices(real scalar T, real scalar L, string scalar btype)
{
    external real colvector rcs_boot_indices
    real scalar i, j, s, pos, nb, pick, p, rem
    real colvector starts

    rcs_boot_indices = J(T,1,.)

    if (btype=="sbb") {
        p = 1/L
        rcs_boot_indices[1] = runiformint(1,1,1,T)

        for (i=2; i<=T; i++) {
            if (runiform(1,1)<p) {
                rcs_boot_indices[i] = runiformint(1,1,1,T)
            }
            else {
                s = rcs_boot_indices[i-1]+1
                if (s>T) s = 1
                rcs_boot_indices[i] = s
            }
        }
        return
    }

    if (btype=="cbb") {
        pos = 1
        while (pos<=T) {
            s = runiformint(1,1,1,T)
            for (j=0; j<L-1; j++) {
                if (pos>T) break
                rcs_boot_indices[pos] = mod(s-1+j,T)+1
                pos++
            }
        }
        return
    }

    if (btype=="mbb") {
        pos = 1
        while (pos<=T) {
            if (L>=T) s = 1
            else s = runiformint(1,1,1,T-L+1)

            for (j=0; j<L-1; j++) {
                if (pos>T) break
                rcs_boot_indices[pos] = s+j
                pos++
            }
        }
        return
    }

    if (btype=="nbb") {
        nb = floor(T/L)

        if (nb<1) {
            rcs_boot_indices = (1..T)'
            return
        }

        starts = (1 :+ (0..nb-1):*L)'
        pos = 1

        while (pos<=T) {
            pick = runiformint(1,1,1,nb)
            s = starts[pick]
            rem = min((L,T-pos+1))

            for (j=0; j<rem; j++) {
                rcs_boot_indices[pos] = s+j
                pos++
            }
        }
        return
    }
}

end
