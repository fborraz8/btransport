# btransport

Stata command for **bounded transport**: combining a small survey that records a
discrete variable (for example, an ability class or a health status) with a large
census or administrative source that records the outcome but not that variable.

Within each cell observed in both sources, bounded transport allows the
distribution of the discrete variable at each rank of the outcome distribution to
differ between the census and the survey by at most δ in total variation. At δ = 0
the survey's information carries over exactly; as δ grows, the identified set
widens until it reaches the Cross–Manski bounds.

For the contrast in the mean outcome of one class between two cells, or the change
in that contrast between two samples, `btransport` reports:

- the value under exact transport;
- δ_min, the smallest discrepancy compatible with the data;
- the breakdown value, the smallest discrepancy that overturns the conclusion,
  with a one-sided confidence bound by subsampling;
- sharp bounds on the contrast at chosen values of δ.

All quantities are values of linear programs. The discrete variable may take any
number of values.

## Installation

Requires Stata 16 or later.

```stata
net install btransport, from("https://raw.githubusercontent.com/fborraz8/btransport/main/")
net get btransport, from("https://raw.githubusercontent.com/fborraz8/btransport/main/")
```

The second line downloads the example do-file to the current directory.

## Example

```stata
btransport lnearn health [pw = wt], census(admin.dta) cweight(wt) cell(agegroup) ///
    contrast(2 1) class(4) delta(0 0.05 0.10 0.20) reps(199) frac(0.25) seed(2026)
```

`btransport_example.do` runs this example on simulated data. See `help btransport`
for all options and stored results.

## Reference

Borraz, F. *Bounded Transport: Combining Survey and Census Data When a Key Variable
Is Missing*. Working paper, Universidad de la República.

## Author

Fernando Borraz, Universidad de la República.
