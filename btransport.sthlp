{smcl}
{* *! version 0.1  25nov2026}{...}
{title:Title}

{phang}
{bf:btransport} {hline 2} Bounded transport: combining a survey with census data when a discrete variable is observed only in the survey


{title:Syntax}

{p 8 17 2}
{cmd:btransport} {it:outcome} {it:discvar} [{it:weight}] [{it:if}] [{it:in}]{cmd:,}
{opt cens:us(filename)} {opt cell(varname)} {opt con:trast(# #)} {opt cl:ass(#)}
[{it:options}]

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt cens:us(filename)}}census or administrative file with {it:outcome} and {it:cellvar}, without {it:discvar}{p_end}
{synopt:{opt cell(varname)}}cells observed in both sources (e.g. schooling groups); same name in both files{p_end}
{synopt:{opt con:trast(hi lo)}}contrast: mean outcome of the class in cell {it:hi} minus cell {it:lo}{p_end}
{synopt:{opt cl:ass(#)}}value of {it:discvar} whose mean outcome is contrasted{p_end}

{syntab:Census}
{synopt:{opt cw:eight(name)}}weight variable in the census file; default equal weights{p_end}
{synopt:{opt cout:come(name)}}outcome variable in the census file; default same name as {it:outcome}{p_end}
{synopt:{opt cif(exp)}}restriction applied to the census file{p_end}

{syntab:Samples}
{synopt:{opt sv:ar(varname)}}variable identifying samples, in both files{p_end}
{synopt:{opt samp:les(s1 [s2])}}one sample, or the change in the contrast between {it:s1} and {it:s2} ({it:s1} minus {it:s2}){p_end}

{syntab:Bounds and discretization}
{synopt:{opt del:ta(numlist)}}report sharp bounds on the contrast at these values of delta{p_end}
{synopt:{opt rank:s(#)}}outcome-rank groups in the survey; default 5{p_end}
{synopt:{opt bin:s(#)}}census rank bins; default 25; must be a multiple of {opt ranks()}{p_end}
{synopt:{opt nomon:otone}}do not impose that the mean outcome increases with {it:discvar}{p_end}
{synopt:{opt met:ric(max|mean)}}maximum (default) or weighted average of the total variation over ranks{p_end}

{syntab:Inference}
{synopt:{opt reps(#)}}subsamples for the confidence bound; default 0 (none){p_end}
{synopt:{opt frac(#)}}subsample fraction of each survey; default 0.25{p_end}
{synopt:{opt seed(#)}}random-number seed{p_end}
{synopt:{opt lpmax:iter(#)}}maximum iterations of each linear program; default 1000{p_end}
{synoptline}
{p 4 6 2}
{cmd:pweight}s and {cmd:aweight}s are allowed for the survey. Frequency and importance weights
are not allowed, because subsampling draws respondents, not weighted units.{p_end}


{title:Description}

{pstd}
{cmd:btransport} combines a small survey that records a discrete variable {it:discvar}
(for example, an ability class) with a large census or administrative source that records
the outcome but not {it:discvar}. Within each cell observed in both sources, bounded transport
allows the distribution of {it:discvar} at each rank of the outcome distribution to differ
between the census and the survey by at most delta in total variation. At delta = 0 the
survey's information carries over exactly; as delta grows, the identified set widens until it
reaches the Cross-Manski bounds. The share of each class within each cell is carried over exactly.

{pstd}
The command reports, for the contrast in the mean outcome of one class between two cells
(or the change in that contrast between two samples): the value under exact transport;
delta_min, the smallest discrepancy compatible with the data (with monotonicity); delta*, the
value at which the restriction stops binding; the sign of the conclusion and its breakdown
value, the smallest discrepancy that overturns it; optionally, a one-sided 95% lower confidence
bound for the breakdown value and a 90% interval for delta_min by subsampling; and sharp
bounds on the contrast at chosen values of delta. All quantities are values of linear programs.

{pstd}
{it:discvar} may take any number K >= 2 of values; classes are ordered by their values.
The census file is read once per cell and sample; it is not modified.


{title:Remarks}

{pstd}
The contrast under exact transport, {cmd:r(D0)}, does not impose monotonicity. All other
quantities impose that the mean outcome increases with {it:discvar} unless {opt nomonotone}
is specified.

{pstd}
Census bins are nested in survey rank groups; hence {opt bins()} must be a multiple of
{opt ranks()}. The command stops with an error if, in the full survey, a class has no
members in a cell or a rank group has no observations.

{pstd}
In the full sample, if a linear program does not converge the command stops with an error;
it never reports a conclusion from a program that was not solved. If the program that
overturns a conclusion is infeasible, the conclusion is not overturned for any discrepancy:
{cmd:r(sign)} gives its sign and {cmd:r(bp)} is missing.

{pstd}
Subsamples in which a program for delta_min or for the breakdown value cannot be solved are
dropped, each counted once: in {cmd:r(ndrop)} when the survey shares are not defined, mainly
because a class has no members in a cell of the subsample, and in {cmd:r(nconv)} when a
program does not converge. A subsample in which the
conclusion cannot be overturned enters with breakdown value 1.

{pstd}
The confidence bound requires that every class has positive share in every cell of the
population; with very small cells, use a larger {opt frac()}.


{title:Examples}

{pstd}Earnings from administrative records, health status (4 classes) from a survey, two age groups:{p_end}
{phang2}{cmd:. btransport lnearn health [pw=wt], census(admin.dta) cweight(wt) cell(agegroup) contrast(2 1) class(4)}{p_end}
{phang2}{cmd:. btransport lnearn health [pw=wt], census(admin.dta) cweight(wt) cell(agegroup) contrast(2 1) class(4) delta(0 .05 .1) reps(199) seed(1)}{p_end}

{pstd}Wages from the census, an ability class from a longitudinal survey, change in the college premium between two cohorts:{p_end}
{phang2}{cmd:. btransport lnw ability [pw=wt], census(census.dta) cweight(perwt) cell(educ) contrast(3 2) class(3) svar(cohort) samples(3 1) reps(199)}{p_end}

{pstd}
See {cmd:btransport_example.do} for a complete example with simulated data.


{title:Stored results}

{synoptset 16 tabbed}{...}
{p2col 5 16 20 2: Scalars}{p_end}
{synopt:{cmd:r(K)}}number of classes of {it:discvar}{p_end}
{synopt:{cmd:r(D0)}}contrast under exact transport, without monotonicity{p_end}
{synopt:{cmd:r(sign)}}sign of the conclusion (1, -1, or 0 if none holds at delta_min){p_end}
{synopt:{cmd:r(dmin)}}smallest discrepancy compatible with the data{p_end}
{synopt:{cmd:r(dstar)}}value at which the restriction stops binding{p_end}
{synopt:{cmd:r(bp)}}breakdown value; missing if no conclusion holds at delta_min or if the conclusion is never overturned{p_end}
{synopt:{cmd:r(lb)}}95% lower confidence bound for the breakdown value (with {opt reps()}){p_end}
{synopt:{cmd:r(dmin_lo)}, {cmd:r(dmin_hi)}}90% interval for delta_min (with {opt reps()}){p_end}
{synopt:{cmd:r(ndrop)}}subsamples dropped because shares are not defined (with {opt reps()}){p_end}
{synopt:{cmd:r(nconv)}}subsamples dropped because a program did not converge (with {opt reps()}){p_end}
{synopt:{cmd:r(nnotover)}}subsamples in which the conclusion cannot be overturned, entered with breakdown value 1 (with {opt reps()}){p_end}
{synopt:{cmd:r(N_survey1)}, {cmd:r(N_survey2)}}survey observations used{p_end}

{p2col 5 16 20 2: Matrices}{p_end}
{synopt:{cmd:r(bounds)}}delta, lower and upper bound on the contrast (with {opt delta()}){p_end}


{title:Reference}

{pstd}
Borraz, F. Bounded transport: combining survey and census data when a key variable is missing.
Working paper.


{title:Author}

{pstd}
Fernando Borraz, Universidad de la Republica.
