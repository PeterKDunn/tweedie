# dtweedie(): methods and regimes (patches 0028–0029)

Nov 4, 2026 · @Pete

dtweedie() no longer interpolates: it chooses, for each value, between a corrected saddlepoint approximation, the series and the Fourier inversion, and its relative error is now below about 4e-13 everywhere tested (before: up to about 3e-5).

## The regimes

Everything depends on p and on z = φ y^(p−2), the dispersion after rescaling y to 1 (by the scaling property, each method's relative accuracy depends on p and z only). For each y > 0 (p not 0, 1, 2 or 3):

| Order | Method | Used where | Accuracy there |
| --- | --- | --- | --- |
| 1 | Saddlepoint approximation with its asymptotic correction (new) | Its estimated relative error < 1e-13, i.e. z small: z < about 0.07 for p from 1.1 to 2.05; 0.03 at p = 2.5; 3e-3 at p = 4; 4e-4 at p = 8; 1e-4 at p = 12; 7e-3 at p = 1.01 | about 1e-13 or better |
| 2 | Series (rewritten for 1 < p < 2) | Otherwise, for 1 < p < 2 | about 1e-13 |
| 3 | Series | Otherwise, for p > 2, where its own error estimate is < 1e-12 (few terms: peak index 1/{z(p−2)} small) | about 1e-13 |
| 4 | Fourier inversion, method 3 (y rescaled to 1) | Otherwise (p > 2, middle values of z) | about 4e-13 |
| — | Fallback: series or saddlepoint, whichever has the smaller error estimate, with a warning | Where the inversion fails | never needed in any test |

Special values stay exact: closed forms at p = 1, 2, 3; Pr(Y = 0) = exp(−λ) for 1 < p < 2; y = Inf gives 0; y = NA gives NA. No interpolation is used anywhere.

The correction: f/f\_sp ≈ 1 + a₁(p) z + a₂(p) z² + … (eight terms), with a₁(p) = p(p−3)/24. The coefficients were derived in exact rational arithmetic from the cumulants. Checks: a\_m(3) = 0 for every m (the saddlepoint is exact for the inverse Gaussian), and a\_m(2) = −1/12, 1/288, 139/51840, … (Stirling's series, the gamma case). For p near 1 the density has ripples the expansion cannot see; their size, exp{−2π²(p−1)/(z(2−p)²)}, is added to the error estimate.

## Accuracy before and after

The new dtweedie() is accurate to about 1e-13 everywhere tested; the old one was off by up to 3e-5 at small dispersion, with no warning.

Worst relative error over a grid of 543 cases (z from 1e-8 to 100, the body and both sides of the mean), against references summed in Rmpfr or the inversion integral in 240-bit arithmetic:

| p | Old dtweedie() | New dtweedie() |
| --- | --- | --- |
| 1.01 | 3e-5 | 1e-13 |
| 1.05 | 8e-6 | 7e-15 |
| 1.2 | 1e-7 | 2e-14 |
| 1.5 | 4e-9 | 2e-14 |
| 1.8 | 5e-8 | 2e-15 |
| 1.95 | 4e-7 | 1e-14 |
| 1.99 | 2e-6 | 1e-14 |
| 2.001 | 1e-5 | 4e-13 |
| 2.01 | 7e-7 | 4e-14 |
| 2.05 | 1e-7 | 4e-14 |
| 2.2 | 1e-7 | 3e-14 |
| 2.5 | 3e-8 | 3e-14 |
| 3.5 | 1e-8 | 1e-13 |
| 4 | 8e-9 | 7e-14 |
| 6 | 9e-8 | 5e-15 |
| 8 | 3e-7 | 5e-15 |
| 12 | 2e-10 | 6e-14 |

Four further checks, none used while building it:

- The 1,614 cases of the regime tool (both tails to 1e-30): worst 2.7e-13 (old: 2.5e-7).
- 69 fresh random cases with new Rmpfr references: worst 1.2e-14.
- The density integrates to the distribution function (∫\_a^b f = F(b) − F(a)) to 1.7e-13 in 120 random cases, and to a total mass of 1 to rounding.
- A random sweep of 3,000 parameter sets (p from 1.001 to 20, φ from 1e-10 to 1e5, μ from 1e-4 to 1e4): no errors, NaNs, negative values or warnings. Where it disagreed with the inversion, Rmpfr references showed dtweedie() right and the inversion wrong.

## Bugs found and fixed

Checking each method on its own turned up these, all fixed in patch 0028:

| Function | Problem | Effect before | Fix |
| --- | --- | --- | --- |
| tweedie\_dev() | Computed as a difference of large terms | Lost most digits near y = μ (relative error up to 1e12 at y/μ − 1 = 1e-12); made the saddlepoint and inversion densities wrong at small φ | Rewritten with log(y/μ), expm1() and a series near y = μ: full accuracy (5e-15) for every p and y |
| dtweedie\_inversion() | The method argument was ignored; the method with the smallest multiplier was always used | Methods 1 and 2 lose all relative accuracy in the tails, and failed (returning 0 or 0.5) for tiny y when φ is large | method is now used; the default, 3, as documented; method = NULL gives the old automatic choice |
| dtweedie\_series(), 1 < p < 2 | One range of N for all values of y | dtweedie\_series(c(1e-10, 1e6), 1.5, 1, 0.01) needed 15 Gb and failed | Each value summed from its own peak |
| dtweedie\_series(), 1 < p < 2 | Each term a difference of numbers of size N log N, then λ subtracted | Relative errors about 1e-9 at λ = 1e4, 1e-7 at λ = 1e6 | Terms from dpois() and dgamma(): 5e-11 or better at any λ |
| dtweedie\_series(), p > 2 | Values with peak index up to 1e6 tried, though double precision cannot sum beyond about 23(p−1); error estimate used the largest term over all values | Huge matrices; values rejected for no reason when y was spread | Only summable values tried; the estimate counts only the terms that matter, and is returned (details = TRUE) |
| dtweedie() and special\_cases() | y = NA gave an error | — | NA gives NA |

One more was found but not fixed (outside the density): ptweedie() with y = NA also gives an error.

## Side effects on ptweedie(), and speed

ptweedie() gets more accurate, because for p > 2 it integrates the density from dtweedie\_inversion() and uses tweedie\_dev(). On the 1,055 reference cases for p > 2 (patch 0027 → with 0028):

| Measure | 0027 | 0027 + 0028 |
| --- | --- | --- |
| Worst unflagged error, 2 < p < 2.05, ξ from 1e4 to 1e7 | 8.6e-10 | 6.3e-10 |
| Worst unflagged error, 2.05 ≤ p < 2.5, ξ from 1e4 to 1e7 | 1.3e-10 | 5.6e-11 |
| Worst unflagged error, 2.05 ≤ p < 2.5, ξ from 1 to 1e3 | 4.6e-7 | 1.1e-12 |
| Errors above 1e-6 (all warned) | 7 | 2 |

One test value moved the other way: p = 2.001, ξ = 3e6, y in the lower tail, from 2e-11 to 1.1e-10 (still within the stated 1e-9); its test tolerance is now 3e-10. Table 4 and the Checks should be regenerated, since some values change in the last digits.

Speed, 10,000 values (seconds; y drawn with the right mean and variance; μ = 2):

| p | φ | Old | New |
| --- | --- | --- | --- |
| 1.1, 1.5, 1.9 | 0.01 to 10 | 0.00 to 0.82 | 0.01 to 0.51 |
| 2.5 | 0.01 / 1 / 10 | 0.01 / 0.01 / 0.04 | 0.02 / 8.5 / 3.9 |
| 4 | 0.01 / 1 / 10 | 0.01 / 0.03 / 0.03 | 16 / 2.5 / 0.56 |
| 8 | 0.01 / 1 / 10 | 0.10 / >120 / >120 | 2.9 / 0.13 / 0.02 |

For 1 < p < 2 the speed is about the same. For p > 2 the new version is much slower where the inversion is needed (middle values of z): about 1–2 ms a value, against almost nothing for the old interpolation. That is where stored accurate values (your planned interpolation) would pay off. The old version hung at p = 8 with φ ≥ 1.

## Patch 0029: the regimes in one file

Patch 0029 (optional, applied after 0028) moves every threshold that picks a method into one function, tweedie\_regime\_settings(), in a new file R/tweedie\_regimes.R. Results are identical to 0028, bit for bit, on all 4,227 test values.

- To change a threshold: edit its one line in tweedie\_regime\_settings(); dtweedie() and ptweedie() read it from there (for example d\_tol\_saddle = 1e-13, p\_xi\_saddle = 1e7).
- To see where values go: tweedie:::dtweedie\_regime(y, mu, phi, power) returns "exact", "saddlepoint", "series" or "series/inversion" for each value; tweedie:::ptweedie\_regime(q, mu, phi, power) returns "series", "inversion", "integration in tails", "integration" or "saddlepoint" (the rows of Table 4).
- Limit: for p > 2 the density's choice between series and inversion, and ptweedie()'s choice of tail or body, are known only after the series or the inversion has been computed, so those two labels name both possibilities.

## Open questions and next steps

- [ ] Apply 0027, 0028, 0029 in that order and run the tests (827 expectations, all pass here; R CMD check: no warnings or notes from the package).
- [ ] Regenerate the Tables and Checks (Table 4 values change in the last digits), including Check B with the skellam\_ref fix.
- [ ] Decide whether to export the corrected saddlepoint (now internal: tweedie:::dtweedie\_saddle\_corrected()).
- [ ] Decide whether dtweedie() should get a log argument: the saddlepoint and series already work on the log scale, and logLiktweedie() takes log(dtweedie()), which underflows in the far tails.
- [ ] Interpolation, later: build it from stored values of the corrected saddlepoint, series and inversion, for speed where the inversion is used (p > 2, middle values of z).
- [ ] ptweedie() with y = NA gives an error (not fixed).
- [ ] dtweedie\_interpolation.R and stored\_grids() are no longer called by dtweedie(); keep them for the interpolation work, or remove them.
