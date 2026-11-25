*! btransport 0.1  25nov2026  Fernando Borraz
*! Bounded transport: combining a survey with census data when a discrete
*! variable is observed only in the survey. Sharp bounds, smallest discrepancy
*! compatible with the data, breakdown value and subsampling confidence bound.
*
* Syntax (see help btransport):
*   btransport y a [pw=w] [if] [in], census(file) cell(var) contrast(hi lo) class(#)
*       [cweight(name) coutcome(name) cif(exp) svar(var) samples(s1 [s2])
*        delta(numlist) ranks(5) bins(25) nomonotone metric(max|mean)
*        reps(0) frac(0.25) seed(#) lpmaxiter(1000)]
*
* Contrast: mean of y in class `class' of a, cell `hi' minus cell `lo'.
* With samples(s1 s2): change in the contrast between samples s1 and s2 (s1 - s2).

program btransport, rclass
    version 16
    syntax varlist(min=2 max=2 numeric) [if] [in] [pweight aweight], ///
        CENSus(string) CELL(varname numeric) CONtrast(numlist min=2 max=2 integer) ///
        CLass(real) [CWeight(name) COUTcome(name) CIF(string asis)              ///
        SVar(varname numeric) SAMPles(numlist min=1 max=2 integer)             ///
        DELta(numlist >=0 <=1) RANKs(integer 5) BINs(integer 25) NOMONotone    ///
        METric(string) REPs(integer 0) FRAC(real 0.25) SEED(string)            ///
        LPMAXiter(integer 1000)]

    tokenize `varlist'
    local y `1'
    local a `2'
    if "`coutcome'" == "" local coutcome `y'
    if "`metric'" == "" local metric max
    if !inlist("`metric'", "max", "mean") {
        di as err "metric() must be max or mean"
        exit 198
    }
    local imet = cond("`metric'" == "max", 1, 2)
    local mono = cond("`nomonotone'" == "", 1, 0)
    if ("`svar'" == "") != ("`samples'" == "") {
        di as err "svar() and samples() must be specified together"
        exit 198
    }
    if `ranks' < 2 | `bins' < `ranks' {
        di as err "ranks() must be at least 2 and bins() at least ranks()"
        exit 198
    }
    if mod(`bins', `ranks') != 0 {
        di as err "bins() must be a multiple of ranks(), so that census bins are nested in survey rank groups"
        exit 198
    }
    _btr_clean
    if `reps' > 0 & (`frac' <= 0 | `frac' >= 1) {
        di as err "frac() must lie in (0, 1)"
        exit 198
    }
    tokenize `contrast'
    local chi `1'
    local clo `2'
    if `chi' == `clo' {
        di as err "contrast() requires two different cells"
        exit 198
    }
    local nsamp 1
    if "`samples'" != "" {
        local nsamp : word count `samples'
        tokenize `samples'
        local s1 `1'
        local s2 `2'
    }

    * ---------------- survey ----------------
    marksample touse
    markout `touse' `cell' `svar'
    tempvar w ak sel1 sel2
    if "`weight'" == "" quietly gen double `w' = 1 if `touse'
    else                quietly gen double `w' `exp' if `touse'
    quietly replace `touse' = 0 if missing(`w') | `w' <= 0
    quietly replace `touse' = 0 if !inlist(`cell', `clo', `chi')
    if "`svar'" != "" {
        if `nsamp' == 1 quietly replace `touse' = 0 if `svar' != `s1'
        else            quietly replace `touse' = 0 if !inlist(`svar', `s1', `s2')
    }
    quietly count if `touse'
    if r(N) == 0 {
        di as err "no survey observations"
        exit 2000
    }
    quietly levelsof `a' if `touse', local(alev)
    local K : word count `alev'
    if `K' < 2 {
        di as err "`a' must take at least two values in the survey"
        exit 198
    }
    local kidx 0
    local i 0
    foreach v of local alev {
        local ++i
        if `v' == `class' local kidx `i'
    }
    if `kidx' == 0 {
        di as err "class(`class') is not a value of `a' in the survey"
        exit 198
    }
    quietly egen `ak' = group(`a') if `touse'
    forvalues i = 1/`nsamp' {
        foreach c in `clo' `chi' {
            local scond 1
            if "`svar'" != "" local scond `svar' == `s`i''
            quietly count if `touse' & `cell' == `c' & `scond'
            if r(N) == 0 {
                di as err "no survey observations in cell `c'" cond("`svar'" != "", " of sample `s`i''", "")
                exit 2000
            }
            local kk 0
            foreach v of local alev {
                local ++kk
                quietly count if `touse' & `cell' == `c' & `scond' & `ak' == `kk'
                if r(N) == 0 {
                    di as err "class `v' of `a' has no survey observations in cell `c'" ///
                        cond("`svar'" != "", " of sample `s`i''", "") "; its share is zero and the programs are not defined"
                    exit 2000
                }
            }
        }
        quietly gen byte `sel`i'' = `touse' & `scond'
        mata: __btr_X`i' = st_data(., ("`y'", "`ak'", "`w'", "`cell'"), "`sel`i''")
        quietly count if `sel`i''
        local Ns`i' = r(N)
    }

    * ---------------- census ----------------
    preserve
    forvalues i = 1/`nsamp' {
        local j 0
        foreach c in `clo' `chi' {
            local ++j
            local sopt
            if "`svar'" != "" local sopt svar(`svar') s(`s`i'')
            _btr_census, file(`"`census'"') y(`coutcome') w(`cweight') cell(`cell') c(`c') ///
                `sopt' cif(`cif') bins(`bins') ranks(`ranks') mname(__btr_S`i'`j')
            local Nc`i'`j' = r(N)
        }
    }
    restore

    * ---------------- computation ----------------
    if "`seed'" != "" set seed `seed'
    mata: btr_run(`nsamp', `kidx', `mono', `imet', `lpmaxiter', `ranks', ///
                  strtoreal(tokens("`delta'")), `reps', `frac', `clo', `chi')

    tempname RES
    matrix `RES' = __btr_res
    local D0   = `RES'[1, 1]
    local sg   = `RES'[1, 2]
    local dmin = `RES'[1, 3]
    local ds   = `RES'[1, 4]
    local bp   = `RES'[1, 5]
    local lb   = `RES'[1, 6]
    local dlo  = `RES'[1, 7]
    local dhi  = `RES'[1, 8]
    local nd   = `RES'[1, 10]
    local nc   = `RES'[1, 12]
    local nbf  = `RES'[1, 13]
    local nno  = `RES'[1, 14]

    * ---------------- output ----------------
    di _n as txt "Bounded transport" _col(50) "Classes of `a': K = " as res `K'
    di as txt "Contrast: mean of `y' for class `class' of `a', `cell' = `chi' minus `cell' = `clo'"
    if `nsamp' == 2 di as txt "          change between `svar' = `s1' and `svar' = `s2' (`s1' minus `s2')"
    else if "`svar'" != "" di as txt "          in `svar' = `s1'"
    di as txt "Metric: `metric' over ranks; monotonicity: " cond(`mono', "imposed", "not imposed") ///
       "; `ranks' survey rank groups, `bins' census bins"
    di as txt "{hline 64}"
    di as txt "Contrast under exact transport (delta = 0)" _col(50) as res %10.4f `D0'
    di as txt "delta_min (smallest compatible discrepancy)" _col(50) as res %10.4f `dmin'
    di as txt "delta* (restriction stops binding)" _col(50) as res %10.4f `ds'
    if `sg' == 0 {
        di as txt "No conclusion on the sign of the contrast holds at delta_min"
    }
    else {
        di as txt "Conclusion: contrast " cond(`sg' > 0, "> 0", "< 0")
        if missing(`bp') di as txt "Breakdown value: the conclusion is not overturned for any delta"
        else             di as txt "Breakdown value, delta_bp" _col(50) as res %10.4f `bp'
        if `reps' > 0 & !missing(`bp') {
            di as txt "95% lower confidence bound (subsampling)" _col(50) as res %10.4f `lb'
        }
        else if `reps' > 0 {
            di as txt "95% lower confidence bound (subsampling)" _col(50) as txt "not applicable"
        }
    }
    if `reps' > 0 {
        di as txt "90% interval for delta_min (subsampling)" _col(50) as res "[" %6.4f `dlo' ", " %6.4f `dhi' "]"
        local frs = cond(`frac' < 1, "0" + string(`frac'), string(`frac'))
        di as txt "Subsamples: `reps', fraction `frs'; dropped (shares not defined): " as res `nd' ///
           as txt "; dropped (no convergence): " as res `nc'
        if `nc' > 0 di as txt "(increase lpmaxiter() if many subsamples do not converge)"
        if `nno' > 0 & `nno' < . di as txt "Subsamples in which the conclusion cannot be overturned (entered as 1): " as res `nno'
    }
    if "`delta'" != "" {
        tempname BND
        matrix `BND' = __btr_bnd
        matrix colnames `BND' = delta lower upper
        di as txt "{hline 64}"
        di as txt "Sharp bounds on the contrast" _n "      delta       lower       upper"
        local anymiss 0
        forvalues r = 1/`=rowsof(`BND')' {
            di as res %11.4f `BND'[`r', 1] %12.4f `BND'[`r', 2] %12.4f `BND'[`r', 3]
            if missing(`BND'[`r', 2]) | missing(`BND'[`r', 3]) local anymiss 1
        }
        if `anymiss' & `nbf' == 0 di as txt "(missing: delta below delta_min, no feasible census distribution)"
        if `nbf' > 0 di as err "warning: `nbf' bound(s) missing because the linear program did not converge; increase lpmaxiter()"
        return matrix bounds = `BND'
    }
    di as txt "{hline 64}"
    _btr_clean

    return scalar K     = `K'
    return scalar D0    = `D0'
    return scalar sign  = `sg'
    return scalar dmin  = `dmin'
    return scalar dstar = `ds'
    return scalar bp    = `bp'
    if `reps' > 0 {
        return scalar lb      = `lb'
        return scalar dmin_lo = `dlo'
        return scalar dmin_hi = `dhi'
        return scalar reps    = `reps'
        return scalar ndrop   = `nd'
        return scalar nconv   = `nc'
        return scalar nnotover = `nno'
    }
    forvalues i = 1/`nsamp' {
        return scalar N_survey`i' = `Ns`i''
    }
    return local metric "`metric'"
end

* removes the command's global objects (Mata and Stata)
program _btr_clean
    capture mata: mata drop __btr_X1 __btr_X2 __btr_S11 __btr_S12 __btr_S21 __btr_S22 __btr_lpst __btr_lpcode __btr_lptext
    capture matrix drop __btr_res __btr_bnd
end

* census distribution of one cell: bins (weight s, mean m, survey rank group q)
program _btr_census, rclass
    syntax , file(string) y(name) [w(string)] cell(name) c(integer) ///
             [svar(name) s(integer -999) cif(string asis)] bins(integer) ranks(integer) mname(name)
    local cond `cell' == `c'
    if "`svar'" != "" local cond `cond' & `svar' == `s'
    if `"`cif'"' != "" local cond `cond' & (`cif')
    quietly use if `cond' using `"`file'"', clear
    if "`w'" == "" {
        tempvar one
        quietly gen double `one' = 1
        local w `one'
    }
    quietly keep `y' `w'
    quietly drop if missing(`y', `w') | `w' <= 0
    if _N == 0 {
        di as err "no census observations in cell `c'" cond("`svar'" != "", " of sample `s'", "")
        exit 2000
    }
    return scalar N = _N
    collapse (sum) __w = `w', by(`y')                      // tied values form one group
    sort `y'
    quietly gen double __cw = sum(__w)
    quietly egen double __tw = total(__w)
    quietly gen double __mr = (__cw - __w/2)/__tw           // mid-rank of the group
    quietly gen int __bin = min(`bins', floor(__mr*`bins') + 1)
    collapse (rawsum) __sw = __w (mean) __m = `y' __rbar = __mr [aw = __w], by(__bin)
    quietly egen double __tot = total(__sw)
    quietly gen double __s = __sw/__tot
    quietly gen int __q = min(`ranks', floor(__rbar*`ranks') + 1)
    mata: `mname' = st_data(., ("__s", "__m", "__q"))
end

version 16
mata:
mata set matastrict off

// share of class k by outcome rank group (K x R); X = (y, class 1..K, weight)
real matrix btr_Q(real matrix X, real scalar K, real scalar R)
{
    real scalar n, i, j, tot, cum, gw, q, k, den
    real matrix Y, Q
    real colvector mr, qt
    n = rows(X)
    Q = J(K, R, 0)
    if (n == 0) return(Q)
    Y = (X[., (1, 2, 3)], (1::n))
    Y = Y[order(Y, (1, 4)), .]
    tot = sum(Y[., 3])
    mr = J(n, 1, .)
    cum = 0
    i = 1
    while (i <= n) {
        j = i
        while (j < n) {
            if (Y[j + 1, 1] == Y[i, 1]) j++
            else break
        }
        gw = sum(Y[|i, 3 \ j, 3|])
        mr[|i \ j|] = J(j - i + 1, 1, (cum + gw/2)/tot)
        cum = cum + gw
        i = j + 1
    }
    qt = rowmin((J(n, 1, R), floor(mr*R) :+ 1))
    for (q = 1; q <= R; q++) {
        den = sum(Y[., 3] :* (qt :== q))
        if (den <= 0) {
            Q[., q] = J(K, 1, .)        // empty rank group: shares not defined
            continue
        }
        for (k = 1; k <= K; k++) Q[k, q] = sum(Y[., 3] :* (qt :== q) :* (Y[., 2] :== k))/den
    }
    return(Q)
}

// Q of each block: (sample 1 cell lo, sample 1 cell hi, [sample 2 lo, sample 2 hi])
pointer(real matrix) rowvector btr_Qs(pointer(real matrix) rowvector pX, real scalar nsamp,
                                      real scalar K, real scalar R, real scalar clo, real scalar chi)
{
    pointer(real matrix) rowvector pQ
    real scalar i
    real matrix X
    pQ = J(1, 2*nsamp, NULL)
    for (i = 1; i <= nsamp; i++) {
        X = *pX[i]
        pQ[2*i - 1] = &(btr_Q(select(X[., (1, 2, 3)], X[., 4] :== clo), K, R))
        pQ[2*i]     = &(btr_Q(select(X[., (1, 2, 3)], X[., 4] :== chi), K, R))
    }
    return(pQ)
}

// Linear program over rho (K x J per block), u = |rho - q| and t.
//   mode 0: min t (k = 0: delta_min ; k > 0: breakdown, adds sense * contrast <= 0)
//   mode 1: min (sense = 1) or max (sense = -1) of the contrast subject to t <= delta
//   metric 1: maximum total variation over ranks ; 2: weighted average
// Returns missing if the program is not defined or not solved; the reason is stored in
// __btr_lpst: 0 solved, 1 infeasible (LinearProgram error code 1), 2 not defined
// (empty rank group or zero class share), 3 no convergence or other error.
// For status 3, __btr_lpcode and __btr_lptext keep LinearProgram's errorcode() and errortext().
real scalar btr_lp(pointer(real matrix) rowvector pS, pointer(real matrix) rowvector pQ,
                   real rowvector cf, real scalar k, real scalar mode, real scalar sense,
                   real scalar delta, real scalar mono, real scalar metric, real scalar maxit)
{
    real scalar B, b, K, n, tix, Jb, off, uo, j, kk, a, val, neq, nin, ie, ii
    real rowvector s, m, sm, c, lb, ub, offs, Js
    real colvector pi, beq, bin
    real matrix SM, Q, R0, Aeq, Ain
    class LinearProgram scalar lp
    external real scalar __btr_lpst, __btr_lpcode
    external string scalar __btr_lptext

    __btr_lpst = 0
    __btr_lpcode = 0
    __btr_lptext = ""
    B = cols(pS)
    K = rows(*pQ[1])
    for (b = 1; b <= B; b++) {
        if (missing(*pQ[b]) > 0) {
            __btr_lpst = 2
            return(.)
        }
    }
    offs = J(1, B, .)
    Js = J(1, B, .)
    n = 0
    neq = 0
    nin = 0
    for (b = 1; b <= B; b++) {
        offs[b] = n
        Js[b] = rows(*pS[b])
        n = n + 2*K*Js[b]
        neq = neq + Js[b] + K - 1
        nin = nin + 2*K*Js[b] + (metric == 1 ? Js[b] : 1) + (mono ? K - 1 : 0)
    }
    tix = n + 1
    n = n + 1
    if (mode == 0 & k > 0) nin = nin + 1
    Aeq = J(neq, n, 0)
    beq = J(neq, 1, 0)
    Ain = J(nin, n, 0)
    bin = J(nin, 1, 0)
    c   = J(1, n, 0)
    ie = 0
    ii = 0
    for (b = 1; b <= B; b++) {
        SM = *pS[b]
        Q  = *pQ[b]
        Jb = Js[b]
        off = offs[b]
        uo  = off + K*Jb
        s  = SM[., 1]'
        m  = SM[., 2]'
        sm = s :* m
        R0 = J(K, Jb, .)
        for (j = 1; j <= Jb; j++) R0[., j] = Q[., SM[j, 3]]
        pi = R0*s'
        if (min(pi) < 1e-8) {
            __btr_lpst = 2
            return(.)
        }
        // (C1) sum_k rho_kj = 1
        for (j = 1; j <= Jb; j++) {
            ie++
            for (kk = 1; kk <= K; kk++) Aeq[ie, off + (kk - 1)*Jb + j] = 1
            beq[ie] = 1
        }
        // (C2) sum_j s_j rho_kj = pi_k, k = 1..K-1 (the K-th is implied)
        for (kk = 1; kk <= K - 1; kk++) {
            ie++
            Aeq[|ie, off + (kk - 1)*Jb + 1 \ ie, off + kk*Jb|] = s
            beq[ie] = pi[kk]
        }
        // u_kj >= |rho_kj - q_kj|
        for (kk = 1; kk <= K; kk++) {
            for (j = 1; j <= Jb; j++) {
                ii++
                Ain[ii, off + (kk - 1)*Jb + j] = 1
                Ain[ii, uo + (kk - 1)*Jb + j] = -1
                bin[ii] = R0[kk, j]
                ii++
                Ain[ii, off + (kk - 1)*Jb + j] = -1
                Ain[ii, uo + (kk - 1)*Jb + j] = -1
                bin[ii] = -R0[kk, j]
            }
        }
        // total variation: (1/2) sum_k u_kj <= t in every bin (max) or on average (mean)
        if (metric == 1) {
            for (j = 1; j <= Jb; j++) {
                ii++
                for (kk = 1; kk <= K; kk++) Ain[ii, uo + (kk - 1)*Jb + j] = 0.5
                Ain[ii, tix] = -1
            }
        }
        else {
            ii++
            for (kk = 1; kk <= K; kk++) Ain[|ii, uo + (kk - 1)*Jb + 1 \ ii, uo + kk*Jb|] = 0.5*s
            Ain[ii, tix] = -1
        }
        // monotonicity: mu_1 <= ... <= mu_K
        if (mono) {
            for (a = 1; a <= K - 1; a++) {
                ii++
                Ain[|ii, off + (a - 1)*Jb + 1 \ ii, off + a*Jb|] = sm/pi[a]
                Ain[|ii, off + a*Jb + 1 \ ii, off + (a + 1)*Jb|] = -sm/pi[a + 1]
            }
        }
        // contrast: sum_b cf_b mu_k(b)
        if (k > 0) c[|off + (k - 1)*Jb + 1 \ off + k*Jb|] = c[|off + (k - 1)*Jb + 1 \ off + k*Jb|] + cf[b]*sm/pi[k]
    }
    lb = J(1, n, 0)
    ub = J(1, n, 1)
    lp = LinearProgram()
    if (mode == 0) {
        if (k > 0) {
            ii++
            Ain[ii, .] = sense*c
        }
        c = J(1, n, 0)
        c[tix] = 1
        lp.setMaxOrMin("min")
    }
    else {
        ub[tix] = delta
        lp.setMaxOrMin(sense == 1 ? "min" : "max")
    }
    lp.setCoefficients(c)
    lp.setEquality(Aeq, beq)
    lp.setInequality(Ain, bin)
    lp.setBounds(lb, ub)
    lp.setMaxiter(maxit)
    val = lp.optimize()
    if (lp.errorcode() != 0) {
        __btr_lpst = (lp.errorcode() == 1 ? 1 : 3)
        __btr_lpcode = lp.errorcode()
        __btr_lptext = lp.errortext()
        return(.)
    }
    return(val)
}

// error in the full sample: message and exit
void btr_fail(real scalar st, string scalar what)
{
    external real scalar __btr_lpcode
    external string scalar __btr_lptext
    if (st == 3 & __btr_lpcode == 3) errprintf("%s: the linear program did not converge; increase lpmaxiter()\n", what)
    else if (st == 3) errprintf("%s: the linear program failed: %s\n", what, __btr_lptext)
    else if (st == 2) errprintf("%s: survey shares not defined (empty rank group or class)\n", what)
    else errprintf("%s: no census distribution satisfies the restrictions\n", what)
    exit(498)
}

// contrast under exact transport (without monotonicity)
real scalar btr_m0(real matrix SM, real matrix Q, real scalar k)
{
    real rowvector s, m, r
    real scalar j
    s = SM[., 1]'
    m = SM[., 2]'
    r = J(1, rows(SM), .)
    for (j = 1; j <= rows(SM); j++) r[j] = Q[k, SM[j, 3]]
    return(sum(s :* m :* r)/sum(s :* r))
}

// delta*: value at which the restriction stops binding
real scalar btr_dstar(pointer(real matrix) rowvector pS, pointer(real matrix) rowvector pQ, real scalar metric)
{
    real scalar b, j, ds, v
    real matrix SM, Q
    ds = 0
    for (b = 1; b <= cols(pS); b++) {
        SM = *pS[b]
        Q  = *pQ[b]
        if (metric == 1) ds = max((ds, 1 - min(Q)))
        else {
            v = 0
            for (j = 1; j <= rows(SM); j++) v = v + SM[j, 1]*(1 - min(Q[., SM[j, 3]]))
            ds = max((ds, v))
        }
    }
    return(ds)
}

real scalar btr_pctl(real colvector x, real scalar p)
{
    real colvector y
    y = select(x, x :< .)
    if (rows(y) == 0) return(.)
    y = sort(y, 1)
    return(y[max((1, ceil(p*rows(y))))])
}

void btr_run(real scalar nsamp, real scalar k, real scalar mono, real scalar metric,
             real scalar maxit, real scalar R, real rowvector dl, real scalar reps,
             real scalar frac, real scalar clo, real scalar chi)
{
    external real matrix __btr_X1, __btr_X2, __btr_S11, __btr_S12, __btr_S21, __btr_S22
    pointer(real matrix) rowvector pX, pS, pQ, pXb, pQb
    real rowvector cf
    real scalar K, i, b, D0, dmin, bpp, bpm, sg, bp, ds, lb, dlo, dhi, r, n, nm, nd, nc, nbf, st, stp, stm
    real scalar std, stb, nno
    external real scalar __btr_lpst
    real matrix BND, X
    real colvector bb, db, pp

    if (nsamp == 1) {
        pX = (&__btr_X1)
        pS = (&__btr_S11, &__btr_S12)
        cf = (-1, 1)
    }
    else {
        pX = (&__btr_X1, &__btr_X2)
        pS = (&__btr_S11, &__btr_S12, &__btr_S21, &__btr_S22)
        cf = (-1, 1, 1, -1)
    }
    K = 0
    for (i = 1; i <= nsamp; i++) K = max((K, max((*pX[i])[., 2])))
    pQ = btr_Qs(pX, nsamp, K, R, clo, chi)
    // empty rank groups in the full sample: explicit error
    for (b = 1; b <= cols(pQ); b++) {
        if (missing(*pQ[b]) > 0) {
            errprintf("survey rank group(s) with no observations in cell %g%s; reduce ranks()\n",
                      (mod(b, 2) == 1 ? clo : chi),
                      (nsamp == 2 ? " of sample " + strofreal(ceil(b/2)) + " in samples()" : ""))
            exit(498)
        }
    }

    D0 = 0
    for (b = 1; b <= cols(pS); b++) D0 = D0 + cf[b]*btr_m0(*pS[b], *pQ[b], k)
    dmin = btr_lp(pS, pQ, cf, 0, 0, 1, ., mono, metric, maxit)
    if (dmin >= .) btr_fail(__btr_lpst, "delta_min")
    bpp = btr_lp(pS, pQ, cf, k, 0,  1, ., mono, metric, maxit)      // overturns "contrast > 0"
    stp = __btr_lpst
    if (bpp >= . & stp != 1) btr_fail(stp, "breakdown value")
    bpm = btr_lp(pS, pQ, cf, k, 0, -1, ., mono, metric, maxit)      // overturns "contrast < 0"
    stm = __btr_lpst
    if (bpm >= . & stm != 1) btr_fail(stm, "breakdown value")
    // sign: the conclusion that requires the larger discrepancy to be overturned. An
    // infeasible program (code 1) means that conclusion is not overturned by any
    // discrepancy: its sign holds and bp is missing. No conclusion: sg = 0, bp missing.
    sg = 0
    bp = .
    if (bpp < . & bpm < .) {
        if (bpp > bpm + 1e-6) {
            sg = 1
            bp = bpp
        }
        else if (bpm > bpp + 1e-6) {
            sg = -1
            bp = bpm
        }
    }
    else if (bpp >= . & bpm < .) sg = 1
    else if (bpm >= . & bpp < .) sg = -1
    else btr_fail(1, "breakdown value")
    ds = btr_dstar(pS, pQ, metric)

    nbf = 0
    BND = J(cols(dl), 3, .)
    for (i = 1; i <= cols(dl); i++) {
        BND[i, 1] = dl[i]
        BND[i, 2] = btr_lp(pS, pQ, cf, k, 1,  1, dl[i], mono, metric, maxit)
        if (BND[i, 2] >= . & __btr_lpst == 3) nbf++
        BND[i, 3] = btr_lp(pS, pQ, cf, k, 1, -1, dl[i], mono, metric, maxit)
        if (BND[i, 3] >= . & __btr_lpst == 3) nbf++
    }

    lb = .
    dlo = .
    dhi = .
    nd = .
    nc = .
    nno = .
    if (reps > 0) {
        nd = 0
        nc = 0
        nno = 0
        bb = J(reps, 1, .)
        db = J(reps, 1, .)
        for (r = 1; r <= reps; r++) {
            pXb = J(1, nsamp, NULL)
            for (i = 1; i <= nsamp; i++) {
                X  = *pX[i]
                n  = rows(X)
                nm = floor(frac*n)
                pp = order(runiform(n, 1), 1)
                pXb[i] = &(X[pp[|1 \ nm|], .])
            }
            pQb = btr_Qs(pXb, nsamp, K, R, clo, chi)
            db[r] = btr_lp(pS, pQb, cf, 0, 0, 1, ., mono, metric, maxit)
            std = __btr_lpst
            stb = 0
            if (sg != 0 & bp < .) {
                bb[r] = btr_lp(pS, pQb, cf, k, 0, sg, ., mono, metric, maxit)
                stb = __btr_lpst
                // infeasible: the conclusion is not overturned in this subsample -> the maximum, 1
                if (bb[r] >= . & stb == 1) {
                    bb[r] = 1
                    nno++
                }
            }
            // dropped: delta_min or the breakdown program fails (counted once per subsample)
            if (db[r] >= . | (sg != 0 & bp < . & bb[r] >= .)) {
                st = (sg != 0 & bp < . & bb[r] >= . ? stb : std)
                if (st == 3) nc++
                else         nd++
            }
        }
        if (sg != 0 & bp < .) {
            if (bp > dmin + 1e-6) lb = max((0, bp - sqrt(frac)*btr_pctl(bb :- bp, 0.95)))
            else                  lb = 0
        }
        dlo = max((0, dmin - sqrt(frac)*btr_pctl(db :- dmin, 0.95)))
        dhi = dmin - sqrt(frac)*btr_pctl(db :- dmin, 0.05)
    }
    st_matrix("__btr_res", (D0, sg, dmin, ds, bp, lb, dlo, dhi, reps, nd, K, nc, nbf, nno))
    st_matrix("__btr_bnd", BND)
}
end
