*==============================================================================
* btransport_example.do — example with simulated data (K = 4 classes)
*
* Generic application: administrative earnings records (large, without health
* information) combined with a health survey (small, with self-reported health
* status in 4 classes). Cells observed in both sources: two age groups.
* Question: what is the earnings gap between age groups for people in
* excellent health?
* Requires btransport.ado on the adopath (e.g. in the personal ado folder).
*==============================================================================
clear all
set more off
version 16

set seed 12345
set obs 60000
gen byte agegroup = 1 + (runiform() < 0.5)                        // 1 = 25-44, 2 = 45-64
gen double u = runiform() - 0.20*(agegroup == 2)                    // health worsens with age
gen byte health = 1 + (u > 0.15) + (u > 0.45) + (u > 0.75)          // 1 poor ... 4 excellent
gen double lnearn = 9.5 + 0.30*(agegroup == 2) + 0.12*health + rnormal(0, 0.60)
gen double wt = 1

* administrative records: 50,000 people, without health information
tempfile admin
preserve
keep in 1/50000
drop health u
save `admin'
restore

* health survey: 1,500 of the remaining people, with health status
keep in 50001/51500
drop u

* earnings gap between age groups for the excellent-health class:
* exact transport, delta_min, breakdown value with confidence bound, and bounds at several deltas
btransport lnearn health [pw = wt], census(`admin') cweight(wt) cell(agegroup) contrast(2 1) class(4) ///
    delta(0 0.05 0.10 0.20) reps(199) frac(0.25) seed(2026)
return list
