*! rcsardl.ado
*! Version 0.2.0  28sep2026
*! RCS-ARDL estimator with optional panel-synchronous block bootstrap
*!
*! Method:
*! Zehra Yalnız (2026)
*! "Robust Estimation of Equilibrium Adjustment in
*! Cross-Sectionally Dependent Dynamic Panels:
*! The RCS-ARDL Approach"
*!
*! CORE RCS-ARDL:
*! - Long-run coefficients are retained from conventional CS-ARDL.
*! - Robustification applies to the equilibrium-adjustment coefficient.
*! - Dynamic admissibility: -2 < phi_i < 0.
*! - MAD normal-consistency factor: 1.4826.
*! - Huber tuning constant: 1.345.
*! - Robust-stage cross-sectional controls are medians.
*!
*! NEW IN 0.2.0:
*! - crlags(#) option for cross-sectional-average lags.
*! - Optional panel-synchronous block bootstrap.
*! - sbb, cbb, mbb and nbb resampling schemes.
*! - lblock(auto) automatic common block length.
*! - Bootstrap SE and percentile confidence interval.
*!
*! Bootstrap design:
*! The same sampled time indices are applied simultaneously
*! to all cross-sectional units in each replication.
*! This preserves contemporaneous cross-sectional structure
*! while retaining serial dependence within sampled blocks.


capture program drop rcsardl
capture program drop _rcsardl_point


/******************************************************************
 INTERNAL POINT-ESTIMATION ENGINE
******************************************************************/

program define _rcsardl_point, rclass sortpreserve
    version 17.0

    syntax varlist(ts numeric min=3) [if] [in], ///
        PANEL(varname numeric) ///
        TIME(varname numeric) ///
        [CRLAGS(integer 0)]

    marksample touse

    if `crlags' < 0 {
        di as error "crlags() must be a nonnegative integer."
        exit 198
    }

    gettoken depvar indepvars : varlist
    local K : word count `indepvars'

    capture which xtdcce2
    if _rc {
        di as error "rcsardl requires xtdcce2."
        di as error "Install xtdcce2 before running rcsardl."
        exit 499
    }

    quietly xtset `panel' `time'

    /**************************************************************
     Conventional CS-ARDL long-run stage
    **************************************************************/

    local lrlist "L.`depvar'"

    foreach x of local indepvars {
        local lrlist "`lrlist' `x' L.`x'"
    }

    capture quietly xtdcce2 `depvar' if `touse', ///
        lr(`lrlist') ///
        lr_options(ardl) ///
        cr(`depvar' `indepvars') ///
        cr_lags(`crlags') ///
        reportconstant

    if _rc {
        exit _rc
    }

    tempname LR U

    /**************************************************************
     Extract conventional adjustment coefficient
    **************************************************************/

    capture scalar __rcs_ect_cs = _b[lr_`depvar']

    if _rc {
        exit 498
    }

    /**************************************************************
     Extract long-run coefficients by NAME
    **************************************************************/

    matrix `LR' = J(1,`K',.)

    local j = 0

    foreach x of local indepvars {

        local ++j

        capture scalar __rcs_theta = _b[lr_`x']

        if _rc {
            exit 498
        }

        matrix `LR'[1,`j'] = __rcs_theta
    }

    matrix colnames `LR' = `indepvars'


    /**************************************************************
     Frozen-V3 robust adjustment stage
    **************************************************************/

    quietly sort `panel' `time'

    tempvar ect_cs L_ect_cs d_y med_dy med_ect

    quietly gen double `ect_cs' = `depvar' if `touse'

    local j = 0

    foreach x of local indepvars {

        local ++j

        scalar __rcs_theta = `LR'[1,`j']

        quietly replace `ect_cs' = ///
            `ect_cs' - __rcs_theta*`x' ///
            if `touse'
    }


    /**************************************************************
     Error-correction component
    **************************************************************/

    quietly gen double `L_ect_cs' = L.`ect_cs'

    quietly gen double `d_y' = D.`depvar'


    /**************************************************************
     Differenced regressors
    **************************************************************/

    local dxlist
    local meddxlist

    local j = 0

    foreach x of local indepvars {

        local ++j

        tempvar dx`j' meddx`j'

        quietly gen double `dx`j'' = D.`x'

        local dxlist ///
            "`dxlist' `dx`j''"

        local meddxlist ///
            "`meddxlist' `meddx`j''"
    }


    /**************************************************************
     Cross-sectional medians
    **************************************************************/

    quietly bysort `time': ///
        egen double `med_dy' = median(`d_y')

    quietly bysort `time': ///
        egen double `med_ect' = median(`L_ect_cs')


    local j = 0

    foreach x of local indepvars {

        local ++j

        quietly bysort `time': ///
            egen double `meddx`j'' = median(`dx`j'')
    }


    /**************************************************************
     Panel identifiers
    **************************************************************/

    quietly levelsof `panel' if `touse', local(panelids)

    local NTOTAL : word count `panelids'


    /**************************************************************
     Unit-specific robust regressions
    **************************************************************/

    tempfile units

    tempname unitpost

    postfile `unitpost' ///
        double uid phi ///
        using `units', replace


    foreach ii of local panelids {

        capture quietly rreg ///
            `d_y' ///
            `L_ect_cs' ///
            `dxlist' ///
            `med_dy' ///
            `meddxlist' ///
            `med_ect' ///
            if `panel'==`ii' & `touse'


        if _rc==0 {

            capture scalar __rcs_phi = _b[`L_ect_cs']

            if _rc==0 {

                /*
                 Dynamic admissibility:
                 -2 < phi_i < 0
                */

                if (__rcs_phi < 0 & __rcs_phi > -2) {

                    post `unitpost' ///
                        (`ii') ///
                        (__rcs_phi)
                }
            }
        }
    }

    postclose `unitpost'


    /**************************************************************
     Robust cross-sectional aggregation
    **************************************************************/

    preserve

        quietly use `units', clear

        quietly drop if missing(phi)

        quietly count

        local NADM = r(N)


        if `NADM' < 1 {

            restore

            exit 498
        }


        /*
         Median
        */

        quietly summarize phi, detail

        scalar __rcs_med_phi = r(p50)


        /*
         MAD
        */

        quietly gen double dev_phi = ///
            abs(phi-__rcs_med_phi)

        quietly summarize dev_phi, detail

        scalar __rcs_mad_raw = r(p50)


        /*
         Numerical safeguard
        */

        if (__rcs_mad_raw<=0 | missing(__rcs_mad_raw)) {

            scalar __rcs_mad_raw = 0.000001
        }


        /*
         Normal-consistent MAD scale
        */

        scalar __rcs_scale = ///
            1.4826*__rcs_mad_raw


        /*
         Huber cutoff
        */

        scalar __rcs_cutoff = ///
            1.345*__rcs_scale


        /*
         Huber weights
        */

        quietly gen double w_phi = 1

        quietly replace w_phi = ///
            __rcs_cutoff/dev_phi ///
            if dev_phi>__rcs_cutoff & dev_phi>0


        quietly count if w_phi < 1

        local NDOWN = r(N)


        /*
         Huber weighted adjustment coefficient
        */

        quietly gen double wx_phi = ///
            w_phi*phi

        quietly summarize wx_phi, meanonly

        scalar __rcs_sum_wx = r(sum)


        quietly summarize w_phi, meanonly

        scalar __rcs_sum_w = r(sum)


        scalar __rcs_ect_rcs = ///
            __rcs_sum_wx/__rcs_sum_w


        /*
         Unit-specific results matrix
        */

        mkmat uid phi w_phi, matrix(`U')

        matrix colnames `U' = ///
            panel_id ///
            phi_robust ///
            huber_weight

    restore


    /**************************************************************
     Return point-estimation results
    **************************************************************/

    return scalar ect_cs = __rcs_ect_cs
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
    return matrix unit_adjustment = `U'

end



/******************************************************************
 MAIN USER COMMAND
******************************************************************/

program define rcsardl, eclass sortpreserve
    version 17.0

    syntax varlist(ts numeric min=3) [if] [in], ///
        PANEL(varname numeric) ///
        TIME(varname numeric) ///
        [ CRLAGS(integer 0) ///
          BOOTstrap ///
          REPS(integer 499) ///
          TYPE(string) ///
          LBLOCK(string) ///
          SEED(integer -1) ///
          LEVEL(real 95) ]


    marksample touse


    /**************************************************************
     Validate options
    **************************************************************/

    if `crlags' < 0 {

        di as error ///
            "crlags() must be a nonnegative integer."

        exit 198
    }


    if `reps' < 1 {

        di as error ///
            "reps() must be a positive integer."

        exit 198
    }


    if `level' <= 0 | `level' >= 100 {

        di as error ///
            "level() must be between 0 and 100."

        exit 198
    }


    if "`type'" == "" {

        local type "sbb"
    }


    local type = lower("`type'")


    if !inlist("`type'","sbb","cbb","mbb","nbb") {

        di as error ///
            "type() must be sbb, cbb, mbb, or nbb."

        exit 198
    }


    if "`lblock'" == "" {

        local lblock "auto"
    }


    gettoken depvar indepvars : varlist


    /**************************************************************
     Check xtdcce2
    **************************************************************/

    capture which xtdcce2

    if _rc {

        di as error ///
            "rcsardl requires xtdcce2."

        di as error ///
            "Install xtdcce2 before running rcsardl."

        exit 499
    }


    /**************************************************************
     ORIGINAL POINT ESTIMATE
    **************************************************************/

    capture quietly _rcsardl_point ///
        `varlist' if `touse', ///
        panel(`panel') ///
        time(`time') ///
        crlags(`crlags')


    if _rc {

        di as error ///
            "RCS-ARDL point estimation failed (return code = " ///
            _rc ")."

        exit _rc
    }


    /**************************************************************
     Save original results BEFORE bootstrap overwrites r()
    **************************************************************/

    tempname LR U

    matrix `LR' = r(longrun)

    matrix `U' = r(unit_adjustment)


    scalar __main_ect_cs = r(ect_cs)

    scalar __main_ect_rcs = r(ect_rcs)

    scalar __main_ntotal = r(n_total)

    scalar __main_nadm = r(n_admissible)

    scalar __main_nexcluded = r(n_excluded)

    scalar __main_ndown = r(n_downweighted)

    scalar __main_medphi = r(median_phi)

    scalar __main_madraw = r(mad_raw)

    scalar __main_madscale = r(mad_scale)

    scalar __main_cutoff = r(huber_cutoff)


    /**************************************************************
     Bootstrap defaults
    **************************************************************/

    scalar __boot_se = .
    scalar __boot_ll = .
    scalar __boot_ul = .

    scalar __boot_reps_requested = 0
    scalar __boot_reps_success = 0

    scalar __boot_lblock = .

    tempname BOOTDIST


    /**************************************************************
     PANEL-SYNCHRONOUS BLOCK BOOTSTRAP
    **************************************************************/

    if "`bootstrap'" != "" {

        if `seed' != -1 {

            quietly set seed `seed'
        }


        preserve

            /*
             Restrict bootstrap base sample
            */

            quietly keep if `touse'


            /*
             Unique panel-time observations required
            */

            capture isid `panel' `time'

            if _rc {

                restore

                di as error ///
                    "Bootstrap requires unique panel-time observations."

                exit 459
            }


            /*
             Bootstrap currently requires a balanced panel.
             This ensures that synchronous time resampling is
             well defined across all cross-sectional units.
            */

            tempvar __Ti __Nt __srcidx

            quietly bysort `panel': ///
                gen long `__Ti' = _N

            quietly summarize `__Ti', meanonly

            local Tmin = r(min)
            local Tmax = r(max)


            if `Tmin' != `Tmax' {

                restore

                di as error ///
                    "Panel-synchronous block bootstrap currently requires a balanced panel."

                exit 459
            }


            quietly levelsof `panel', local(__panels)

            local __N : word count `__panels'


            quietly bysort `time': ///
                gen long `__Nt' = _N

            quietly summarize `__Nt', meanonly


            if r(min) != `__N' | r(max) != `__N' {

                restore

                di as error ///
                    "Panel-synchronous block bootstrap requires common time periods across units."

                exit 459
            }


            local __T = `Tmax'


            if `__T' < 4 {

                restore

                di as error ///
                    "Too few time observations for block bootstrap."

                exit 2001
            }


            /*
             Create ordered source-time index
            */

            quietly egen long `__srcidx' = group(`time')


            tempfile __bootbase

            quietly keep ///
                `panel' ///
                `depvar' ///
                `indepvars' ///
                `__srcidx'

            quietly save `__bootbase', replace


            /********************************************************
             BLOCK LENGTH
            ********************************************************/

            if lower("`lblock'") == "auto" {

                /*
                 Common block length selected from the
                 cross-sectional median of the dependent variable.

                 This produces one common time-block length for
                 synchronous panel resampling.
                */

                tempvar __medy

                quietly collapse ///
                    (median) `__medy'=`depvar', ///
                    by(`__srcidx')

                quietly sort `__srcidx'


                mata: ///
                    st_numscalar( ///
                    "__rcs_auto_block", ///
                    rcs_auto_block( ///
                    st_data(.,"`__medy'"), ///
                    "`type'" ///
                    ))


                local __L = scalar(__rcs_auto_block)


                quietly use `__bootbase', clear
            }

            else {

                local __L = real("`lblock'")


                if missing(`__L') {

                    restore

                    di as error ///
                        "lblock() must contain a positive integer or auto."

                    exit 198
                }


                if `__L' < 1 | floor(`__L') != `__L' {

                    restore

                    di as error ///
                        "lblock() must contain a positive integer or auto."

                    exit 198
                }


                if `__L' > `__T' {

                    restore

                    di as error ///
                        "Block length cannot exceed the time dimension."

                    exit 198
                }
            }


            scalar __boot_lblock = `__L'


            /********************************************************
             Bootstrap storage
            ********************************************************/

            tempfile __bootresults

            tempname __bootpost


            postfile `__bootpost' ///
                double phi_rcs ///
                using `__bootresults', replace


            /********************************************************
             Bootstrap replications
            ********************************************************/

            forvalues __b = 1/`reps' {

                quietly clear

                quietly set obs `__T'


                tempvar __boottime __bootsrc


                quietly gen long `__boottime' = _n


                /*
                 Mata generates the source time index:
                 same index sequence for EVERY cross-sectional unit.
                */

                mata: ///
                    rcs_make_indices( ///
                    `__T', ///
                    `__L', ///
                    "`type'" ///
                    )


                quietly getmata ///
                    (`__bootsrc') = rcs_boot_indices


                /*
                 Rename source key to match base dataset
                */

                quietly rename ///
                    `__bootsrc' ///
                    `__srcidx'


                /*
                 Match selected times with all panel units.

                 Duplicate selected dates are retained and assigned
                 different bootstrap-time positions.
                */

                capture quietly joinby ///
                    `__srcidx' ///
                    using `__bootbase'


                if _rc==0 {

                    quietly sort ///
                        `panel' ///
                        `__boottime'


                    /*
                     Complete RCS-ARDL is re-estimated in every draw:
                     CS-ARDL -> LR coefficients -> ECT ->
                     unit phi -> admissibility -> MAD/Huber.
                    */

                    capture quietly _rcsardl_point ///
                        `depvar' `indepvars', ///
                        panel(`panel') ///
                        time(`__boottime') ///
                        crlags(`crlags')


                    if _rc==0 {

                        capture scalar ///
                            __boot_phi = r(ect_rcs)


                        if _rc==0 & !missing(__boot_phi) {

                            post `__bootpost' ///
                                (__boot_phi)
                        }
                    }
                }
            }


            postclose `__bootpost'


            /********************************************************
             Bootstrap distribution
            ********************************************************/

            quietly use `__bootresults', clear

            quietly drop if missing(phi_rcs)

            quietly count

            local __BSUCCESS = r(N)


            scalar __boot_reps_requested = `reps'

            scalar __boot_reps_success = `__BSUCCESS'


            if `__BSUCCESS' < 2 {

                restore

                di as error ///
                    "Fewer than two successful bootstrap replications were obtained."

                exit 498
            }


            quietly summarize phi_rcs

            scalar __boot_se = r(sd)


            /*
             Percentile confidence interval
            */

            local __alpha = ///
                (100-`level')/2

            local __upper = ///
                100-`__alpha'


            quietly centile phi_rcs, ///
                centile(`__alpha' `__upper')


            scalar __boot_ll = r(c_1)

            scalar __boot_ul = r(c_2)


            /*
             Save bootstrap draws
            */

            mkmat phi_rcs, matrix(`BOOTDIST')

            matrix colnames `BOOTDIST' = ///
                ect_rcs_boot


        restore
    }


    /**************************************************************
     POST e() RESULTS
    **************************************************************/

    ereturn clear


    /**************************************************************
     Point estimates
    **************************************************************/

    ereturn scalar ect_cs = ///
        __main_ect_cs

    ereturn scalar ect_rcs = ///
        __main_ect_rcs


    ereturn scalar cr_lags = ///
        `crlags'


    /**************************************************************
     Unit-level robust-stage diagnostics
    **************************************************************/

    ereturn scalar n_total = ///
        __main_ntotal

    ereturn scalar n_admissible = ///
        __main_nadm

    ereturn scalar n_excluded = ///
        __main_nexcluded

    ereturn scalar n_downweighted = ///
        __main_ndown


    ereturn scalar median_phi = ///
        __main_medphi

    ereturn scalar mad_raw = ///
        __main_madraw

    ereturn scalar mad_scale = ///
        __main_madscale


    ereturn scalar huber_cutoff = ///
        __main_cutoff

    ereturn scalar huber_c = ///
        1.345


    ereturn scalar admiss_lower = ///
        -2

    ereturn scalar admiss_upper = ///
        0


    /**************************************************************
     Main matrices
    **************************************************************/

    ereturn matrix longrun = `LR'

    ereturn matrix unit_adjustment = `U'


    /**************************************************************
     Bootstrap results
    **************************************************************/

    if "`bootstrap'" != "" {

        ereturn scalar boot_se = ///
            __boot_se

        ereturn scalar boot_ll = ///
            __boot_ll

        ereturn scalar boot_ul = ///
            __boot_ul

        ereturn scalar boot_reps = ///
            __boot_reps_requested

        ereturn scalar boot_success = ///
            __boot_reps_success

        ereturn scalar boot_lblock = ///
            __boot_lblock


        ereturn matrix boot_dist = ///
            `BOOTDIST'


        ereturn local bootstrap ///
            "yes"

        ereturn local boot_type ///
            "`type'"

        ereturn local boot_lblock_method ///
            "`lblock'"

        ereturn scalar level = ///
            `level'
    }

    else {

        ereturn local bootstrap ///
            "no"
    }


    /**************************************************************
     Identification macros
    **************************************************************/

    ereturn local cmd ///
        "rcsardl"

    ereturn local depvar ///
        "`depvar'"

    ereturn local indepvars ///
        "`indepvars'"

    ereturn local panelvar ///
        "`panel'"

    ereturn local timevar ///
        "`time'"

    ereturn local method ///
        "Frozen V3 RCS-ARDL adjustment estimator"


    /**************************************************************
     DISPLAY
    **************************************************************/

    di

    di as text ///
        "RCS-ARDL (Frozen V3 adjustment estimator)"

    di as text ///
        "Long-run coefficients retained from conventional CS-ARDL"

    di as text "{hline 72}"


    di as text ///
        "Cross-sectional average lags" ///
        _col(55) ///
        as result %10.0f e(cr_lags)


    di as text ///
        "Conventional CS-ARDL adjustment coefficient" ///
        _col(55) ///
        as result %10.4f e(ect_cs)


    di as text ///
        "RCS-ARDL robust adjustment coefficient" ///
        _col(55) ///
        as result %10.4f e(ect_rcs)


    di as text ///
        "Admissible units (-2 < phi_i < 0)" ///
        _col(55) ///
        as result ///
        %6.0f e(n_admissible) ///
        "/" ///
        %6.0f e(n_total)


    di as text ///
        "Huber down-weighted units" ///
        _col(55) ///
        as result %10.0f e(n_downweighted)


    di as text ///
        "Median admissible phi" ///
        _col(55) ///
        as result %10.4f e(median_phi)


    di as text ///
        "MAD-based robust scale (1.4826 x MAD)" ///
        _col(55) ///
        as result %10.4f e(mad_scale)


    di as text ///
        "Huber cutoff (1.345 x scale)" ///
        _col(55) ///
        as result %10.4f e(huber_cutoff)


    /**************************************************************
     Bootstrap display
    **************************************************************/

    if "`bootstrap'" != "" {

        di as text "{hline 72}"

        di as text ///
            "Panel-synchronous block bootstrap"


        di as text ///
            "Bootstrap type" ///
            _col(55) ///
            as result "`type'"


        di as text ///
            "Block length" ///
            _col(55) ///
            as result %10.0f e(boot_lblock)


        di as text ///
            "Requested bootstrap replications" ///
            _col(55) ///
            as result %10.0f e(boot_reps)


        di as text ///
            "Successful bootstrap replications" ///
            _col(55) ///
            as result %10.0f e(boot_success)


        di as text ///
            "Bootstrap SE of RCS adjustment coefficient" ///
            _col(55) ///
            as result %10.4f e(boot_se)


        di as text ///
            "`level'% percentile CI - lower" ///
            _col(55) ///
            as result %10.4f e(boot_ll)


        di as text ///
            "`level'% percentile CI - upper" ///
            _col(55) ///
            as result %10.4f e(boot_ul)
    }


    di as text "{hline 72}"

    di as text ///
        "Long-run coefficient vector:"


    matrix list e(longrun), ///
        noheader ///
        format(%10.4f)

end



/******************************************************************
 MATA ROUTINES
******************************************************************/

mata:


/******************************************************************
 Automatic block-length selector

 Adapted to provide one common block length for the
 panel-synchronous RCS-ARDL bootstrap.

 For SBB: stationary-bootstrap block-length expression.
 For CBB/MBB/NBB: fixed-block expression.
******************************************************************/

real scalar rcs_auto_block(
    real colvector x,
    string scalar btype
)
{
    real scalar T
    real scalar Kn
    real scalar mmax
    real scalar Bmax
    real scalar crit
    real scalar meanx
    real scalar varx

    real scalar k
    real scalar j
    real scalar first_run
    real scalar mhat
    real scalar M

    real scalar Ghat
    real scalar Dhat
    real scalar Bstar

    real colvector rho
    real colvector insignif
    real colvector acov
    real colvector kk
    real colvector lam

    real colvector x1
    real colvector x2


    /*
     Remove missing values
    */

    x = select(x, x :< .)

    T = rows(x)


    if (T < 4) {

        return(1)
    }


    Kn = max((5,ceil(log10(T))))

    mmax = ceil(sqrt(T)) + Kn

    mmax = min((mmax,T-1))


    Bmax = min((ceil(3*sqrt(T)),floor(T/3)))

    Bmax = max((1,Bmax))


    crit = ///
        invnormal(0.975)*sqrt(log10(T)/T)


    meanx = mean(x)

    varx = mean((x :- meanx):^2)


    /*
     Degenerate series safeguard
    */

    if (varx<=0 | missing(varx)) {

        return(1)
    }


    rho = J(mmax,1,.)

    insignif = J(mmax,1,0)


    /*
     Sample autocorrelations
    */

    for (k=1; k<=mmax; k++) {

        x1 = x[1..(T-k)]

        x2 = x[(1+k)..T]


        rho[k] = ///
            mean((x1:-meanx):*(x2:-meanx)) ///
            / varx ///
            * ((T-k)/T)


        insignif[k] = ///
            abs(rho[k]) < crit
    }


    /*
     Locate first run of Kn insignificant autocorrelations
    */

    first_run = 0


    if (mmax>=Kn) {

        for (k=1; k<=mmax-Kn+1; k++) {

            if (sum(insignif[k..(k+Kn-1)])==Kn) {

                first_run = k

                break
            }
        }
    }


    if (first_run>0) {

        mhat = first_run
    }

    else {

        mhat = mmax
    }


    M = min((2*mhat,mmax))

    M = max((1,M))


    /*
     Autocovariances from -M to M
    */

    acov = J(2*M+1,1,.)

    kk = (-M..M)'


    for (j=1; j<=rows(kk); j++) {

        k = kk[j]


        if (k<0) {

            x1 = x[(1-k)..T]

            x2 = x[1..(T+k)]
        }

        else if (k>0) {

            x1 = x[1..(T-k)]

            x2 = x[(1+k)..T]
        }

        else {

            x1 = x

            x2 = x
        }


        acov[j] = ///
            mean((x1:-meanx):*(x2:-meanx)) ///
            * (rows(x1)/T)
    }


    /*
     Flat-top weighting used in block-length calculation
    */

    lam = J(rows(kk),1,0)


    for (j=1; j<=rows(kk); j++) {

        if (abs(kk[j]/M)<0.5) {

            lam[j] = 1
        }

        else if (abs(kk[j]/M)<=1) {

            lam[j] = ///
                2*(1-abs(kk[j]/M))
        }
    }


    Ghat = ///
        sum(lam:*abs(kk):*acov)


    /*
     Stationary bootstrap
    */

    if (btype=="sbb") {

        Dhat = ///
            2*(sum(lam:*acov)^2)
    }

    /*
     Circular / moving / non-overlapping
    */

    else {

        Dhat = ///
            (4/3)*(sum(lam:*acov)^2)
    }


    /*
     Numerical fallback
    */

    if (Dhat<=0 | missing(Dhat) | missing(Ghat)) {

        Bstar = round(T^(1/3))
    }

    else {

        Bstar = ///
            ((2*Ghat^2)/Dhat)^(1/3) ///
            * T^(1/3)
    }


    if (missing(Bstar)) {

        Bstar = round(T^(1/3))
    }


    if (btype=="sbb") {

        Bstar = round(Bstar)
    }

    else {

        Bstar = ceil(Bstar)
    }


    Bstar = ///
        max((1,min((Bstar,Bmax))))


    return(Bstar)
}



/******************************************************************
 Generate one synchronous bootstrap time-index sequence
******************************************************************/

void rcs_make_indices(
    real scalar T,
    real scalar L,
    string scalar btype
)
{
    external real colvector rcs_boot_indices

    real scalar i
    real scalar j
    real scalar s
    real scalar pos
    real scalar nb
    real scalar pick
    real scalar p

    real colvector starts


    rcs_boot_indices = J(T,1,.)


    /**************************************************************
     STATIONARY BLOCK BOOTSTRAP
    **************************************************************/

    if (btype=="sbb") {

        p = 1/L


        rcs_boot_indices[1] = ///
            runiformint(1,1,1,T)


        for (i=2; i<=T; i++) {

            if (runiform(1,1)<p) {

                rcs_boot_indices[i] = ///
                    runiformint(1,1,1,T)
            }

            else {

                s = ///
                    rcs_boot_indices[i-1] + 1


                if (s>T) {

                    s = 1
                }


                rcs_boot_indices[i] = s
            }
        }

        return
    }


    /**************************************************************
     CIRCULAR BLOCK BOOTSTRAP
    **************************************************************/

    if (btype=="cbb") {

        pos = 1


        while (pos<=T) {

            s = ///
                runiformint(1,1,1,T)


            for (j=0; j<L-1; j++) {

                if (pos>T) break


                rcs_boot_indices[pos] = ///
                    mod(s-1+j,T)+1


                pos++
            }
        }

        return
    }


    /**************************************************************
     MOVING BLOCK BOOTSTRAP
    **************************************************************/

    if (btype=="mbb") {

        pos = 1


        while (pos<=T) {

            if (L>=T) {

                s = 1
            }

            else {

                s = ///
                    runiformint( ///
                    1,1,1,T-L+1 ///
                    )
            }


            for (j=0; j<L-1; j++) {

                if (pos>T) break


                rcs_boot_indices[pos] = ///
                    s+j


                pos++
            }
        }

        return
    }


    /**************************************************************
     NONOVERLAPPING BLOCK BOOTSTRAP
    **************************************************************/

    if (btype=="nbb") {

        nb = floor(T/L)


        if (nb<1) {

            rcs_boot_indices = ///
                (1..T)'

            return
        }


        starts = ///
            (1 :+ (0..nb-1):*L)'


        pos = 1


        while (pos<=T) {

            pick = ///
                runiformint(1,1,1,nb)


            s = starts[pick]


            for (j=0; j<L-1; j++) {

                if (pos>T) break

                if (s+j>T) break


                rcs_boot_indices[pos] = ///
                    s+j


                pos++
            }
        }


        /*
         Safety fill if T is not an exact multiple of L
        */

        while (pos<=T) {

            rcs_boot_indices[pos] = ///
                runiformint(1,1,1,T)

            pos++
        }

        return
    }

}

end
