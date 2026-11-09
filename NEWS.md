tweedie 3.2.0 (Release date: )
==============

## Summary
Numerous changes, mostly internally. These changes substantially increase the
relative accuracy of the algorithms for dtweedie() and especially for ptweedie().

**Crucially, the interpolation is no longer used for dtweedie()** but will be
reinstated in a future update; the purpose of this update is improving accuracy.
The result is much greater relative accuracy, but some evaluations are slightly 
slower (though not painfully so). This will also mean that the faster
interpolation method will produce more accurate results when it is reinstated.

Some of the important changes for improving accuracy include: 
* new methods are introduced for computation (e.g., direct integration for the CDF); 
* improved (corrected) saddlepoint approximations;  
* more computations performed on the log-scale to prevent under- and over-flow; 
* a re-scaling identity now is used for the CDF.

dtweedie() and ptweedie() have new inputs for log-outputs. More warnings and
errors are reported when accuracy targets are not met.

These changes have all been prompted by preparing the Tweedie CDF paper. For
example, the direct integration method inspired the more accurate computation
of the density.


## Changed and fixed
* ptweedie_inversion(): with details = TRUE, now also returns small_tail, flagging
  tail probabilities below 1e-8. The inversion computes F = 1/2 -/+ I/pi, so its
  error is absolute (a few times 1e-16), and such small tails can have poor relative
  accuracy. The new Accuracy section of the help page explains this. No warning is
  given (as for the other methods), and the values returned are unchanged; ptweedie()
  already recomputes small tails by integrating the density.
* dtweedie_saddle(), dtweedie_series(), ptweedie_series(): new Accuracy sections in
  the help pages, describing where each method is less accurate.
* ptweedie() and dtweedie(): Numerous improvements that improve accuracy and precision
* ptweedie() and dtweedie(): Numerous new warnings in the code when accuracy is not reached or fails 
* ptweedie() and dtweedie(): Some refining of the inversion approach
* dtweedie() and ptweedie(): The thresholds that decide which method dtweedie() and ptweedie() use are now all in one place, tweedie_regime_settings() (R/tweedie_regimes.R); the internal functions dtweedie_regime() and ptweedie_regime() say which regime each value falls in. The results are unchanged
* ptweedie(): Use the re-scaling identity to simplify CDF computations
* ptweedie(): Increased a computation limit to allows greater precision
* ptweedie(): new log.p input; if TRUE, the log of the CDF is returned
* ptweedie(): Add lower.tail argument for ptweedie() functions and ptweedie_series()
* ptweedie(): New testthat tests for the log.p and lower.tail changes in ptweedie()
* ptweedie_series(): Changes in the ptweedie_series() approach as per paper
* ptweedie_series(): New input: max.terms, so that the truncation of the series can be controlled
* ptweedie(): using saddlepoint approx and direct integration, greatly improving accuracy
* ptweedie(): for p > 2, uses the Lugannani-Rice saddlepoint approximation when the dispersion is very small (xi = mu^(2-p)/{phi (p-2)} > 1e7, and in the body when xi > 1e6 and phi mu^(p-2) <= 1e-5), where the inversion and tail integration can fail
* ptweedie(): for p > 2 and xi > 1e4: tail probabilities below 1e-2 (not 1e-5) are computed by integrating the density, and the integrated value is used whenever the integration succeeds (before, between 1e-10 and 1e-5 the inversion value was kept if the two disagreed, which near p = 2 gave silent errors of up to 0.5%)
* ptweedie(): for p > 2 and 1e5 < xi <= 1e7: values in the body are also computed by integrating the density (the inversion's error there reached about 1e-8 near p = 2; now about 1e-9 or better), with the inversion where the integration returns no value; this is slower (about 0.1 s per value)
* ptweedie(): for p > 2: the saddlepoint approximation is also used at xi = 1e7 when phi mu^(p-2) <= 1e-5, where it is more accurate than the integration
* ptweedie(): for p > 2: where the inversion fails (near p = 2, for very small y and large dispersion it could return 0.5, or 0, whatever the true value), the density is integrated or, failing that, the saddlepoint approximation is used; these values are rough, and still flagged by a warning
* ptweedie(): y = NA gives NA (it gave an error)
* ptweedie_inversion(): a result is now flagged (exitstatus = 1, with a warning) when the pre-acceleration zone reaches its cap, since acceleration then starts where it is not valid
* dtweedie(): new argument log; if TRUE, the log of the density is returned, computed without underflow (each method works on the log scale), so it stays finite far into the tails where the density underflows to zero. tweedie_AIC(), logLiktweedie() and tweedie_profile() now use it, so an observation far in the tail no longer makes the log-likelihood -Inf. dtweedie_series() and dtweedie_inversion() return the log density too (details = TRUE: logf, logdensity)
* dtweedie(): no longer interpolates. For each value it uses the saddlepoint approximation with its asymptotic correction where z = phi y^(p-2) is small enough (its estimated relative error below 1e-13), otherwise the series (1 < p < 2; and p > 2 where its estimated error is below 1e-12), otherwise the Fourier inversion. The relative error is now below about 4e-13 everywhere tested (before, up to about 1e-5 for p near 1 at small dispersion, and about 2.5e-7 for 7 <= p <= 10)
* dtweedie_series() for 1 < p < 2: each term is computed with dpois() and dgamma(), and the terms are summed from each value's own peak, so it keeps full accuracy at large lambda (before, relative errors of about 1e-7 at lambda = 1e6) and no longer needs a huge matrix when the values of y are widely spread (it could exhaust the memory)
* dtweedie_series() for p > 2: only values that can be summed in double precision are tried, and the error estimate counts only the terms that matter; with details = TRUE it also returns the estimated relative error (relerr)
* dtweedie_inversion(): the method argument is now used (it was ignored, and the method with the smallest multiplier was always chosen; method = NULL does that); the default, method 3, keeps the relative accuracy in the tails. This also makes ptweedie() more accurate for p > 2 where it integrates the density
* dtweedie(): y = NA gives NA (it gave an error)
* tweedie_dev(): computed without cancellation, so it has full relative accuracy when y is close to mu (it lost most of its digits there); this also makes the saddlepoint and inversion densities more accurate at small dispersion
* Numerous new testthat tests 
* Minor coding fixes and improvements 
* Minor improvements and typo fixes in man pages throghout


## Fixed
* Fixed typo in README.Rmd.


tweedie 3.1.0 (Release date: 2026-05-20)
==============

## Summary
Minor fixes and improvements

## Changed
* Added  details  input to dtweedie_series(), that returns the  lo  and  hi  index for the summing of terms.
  Other necessary changes made downstream.
* Added (tweedie_integrand(), tweedie_profile()) and improved (tweedie_plot()) the passing of plotting parameters.

## Fixed
* Fixed tweedie_integrand: Plot 2 now plots cos( Im k(t) ) for the PDF (not sin).



tweedie 3.0.19 (Release date: 2026-04-26)
==============

## Summary
Minor fixes and improvements

## Changed
* Improvements to tweedie_plot() and passing plotting options

## Fixed
* Fix error with accelerate when returning early: return properly(Thanks Jeonghwan Lee)


tweedie 3.0.17 (Release date: 2026-02-26)
==============

Changes:

* Fixed an error with xi = 1 and phi \ne 1 in qtweedie
  (thanks Milan Bouchet-Valat).
* Added a test to prevent this again.
* Relocated hex.R and fixed the hex-producing script.
* Fix Makefile.

tweedie 3.0.14 (Release date: 2026-02-16)
==============

Changes:

* Improved the vignette.
* Some internal renaming.
* Fix some xi = 1 cases (thanks Milan Bouchet-Valat).
* Relocate some messages in tweedie_profile().
* Add poison example to vignette.


tweedie 3.0.12 (Release date: 2026-02-07)
==============

Changes:

* Trying to fix bugs that pop up (seemingly at random) with rhub etc. checks.



tweedie 3.0.5 (Release date: 2026-01-30)
==============

Changes:

* FORTRAN code restructured to make the similar flow in the three zones 
  (initial; pre-acceleration; acceleration) clearer
* Some fixes to documentation to pass tests.
* Some minor fixes to R code.



tweedie 3.0.4 (Release date: 2026-01-20)
==============

Changes:

* Some fixes to implementation of IGexact


tweedie 3.0.3 (Release date: 2025-11-29)
==============

Changes:

* Add IGexact for [dp]tweedie_inversion: whether to use exact values or inversion when p = 3.
* Fixed some comments



tweedie 3.0.2 (Release date: 2025-11-29)
==============

Changes:

* All code moved from FORTRAN77 to FORTRAN90. Almost no FORTRAN code remains from version < 3.
* PDF and CDF computations consolidated and code shared where possible, substantially reducing the amount of FORTRAN code.
* Separated FORTRAN code into different files for easier debugging.
* Improvements to the acceleration algorithm and root-finding algorithms, so should work better for more cases.
* Added verbose (shows what's happening behind the scenes) and details (reports on the fitting) as options for many user-facing R functions.
* Added  ptweedie_inversion()  to the man page for  dtweedie.
* Tidied the man pages; added examples.
* Removed the almost-never used dtweedie.stable() function.
* Changed function names (e.g,  tweedie.convert()  to  tweedie_convert()).
* Separated R functions into separate files depending on purpose (e.g., dtweedie.R  and  ptweedie.R).
* Moved the  tweedie_Extra  files into the main package.
* dtweedie.igrand()  (now tweedie_igrand()) to plot the integrand for the DF also.


tweedie 2.3.5 (Release date: 2022-08-17)
==============

Changes:

* Added outputs gamma.mean and gamma.phi to tweedie.convert()
* Added more error checks to  tweedie.convert()  to prevent a model being provided
* Minor edits

tweedie 2.3.1 (Release date: 2017-11-15)
==============

Changes:

* Updated ptweedie.series() to fix a bug (reported by Lu Yang), where incorrect answers could sometimes be returned.
* Other minor fixes

tweedie 2.3.0 (Release date: 2017-11-06)
==============

Changes:

* Fixed an issue with AICtweedie(), where the incorrect AIC was given when prior weights used (reported by David Scollnik)
* Fixed a compilation error (in the subroutine smallp(), where variables were declared as initialised (thanks to Iñaki Úcar <i.ucar86@gmail.com>)
* Other minor edits.

tweedie 2.2.10 (Release date: 2017-08-23)
==============

Changes:

* Kept it even more quiet

tweedie 2.2.6 (Release date: 2017-08-22)
==============

Changes:

* Fixed tweedie.f() to keep it quiet more often (sometimes, diagnostic reports meant for internal monitoring, were printed)
* Fixed a problem reported by Gustavo Lacerda, where  ptweedie()  returned NaN
* As a result, the series is now used in far more cases when  1<xi<2  for  ptweedie()
* Changed  rtweedie()  algorithm for the case 1 < xi < 2 (thanks to Carlos J. Gil Bellosta)
* New function  tweedie.convert()  added

tweedie 2.2.5 (Release date: 2016-12-19)
==============

Changes:

* Added CITATION file
* Fixed an issue where the Tweedie cdf could return a value greater than one (reported by Jeremie Juste)
* Minor tidy of FORTRAN code
* Minor tidies in R code

tweedie 2.1.9 (Release date: 2014-06-06)
==============

Changes:

* Some administrative fixes for CRAN

tweedie 2.1.8 (Release date: 2013-09-10)
==============

Changes:

* Fixed an issue where ptweedie() would fail for very small y; set this to 0 when y<1.0e-300 (based on a report by Johann Cuenin)

tweedie 2.1.7 (Release date: 2013-01-15)
==============

Changes:

* Admin release (e.g. .First.lib()  removed)
* Some minor edits to manual
* Added the  control  input to  tweedie.profile() (thanks to Giri Khageswor, DPI Victoria)
* Minor fixes in the manual

tweedie 2.1.5 (Release date: 2012-11-01)
==============

Changes:

* Changed the example in  tweedie-package  to execute faster (CRAN requirement)

tweedie 2.1.4 (Release date: 2012-10-31)
==============

Changes:

* Fixed an error in dtweedie() that reported NA in the case when power=1 and phi != 1 (thanks to Dina Farkas)
* [dqpr]stable now in package  stabledist  rather then fBasics; fixed
* Fixed some typos in the help for tweedie-package (thanks to Peng Yu)
* rtweedie()  reported an error if power = 1; fixed (thanks for Peng Yu)
* Edits to conform with new first argument of .Fortran (i.e. .NAME rather than name)
* Added NAMESPACE

tweedie 2.1.0 (Release date: 2011-06-09)
==============

Changes:

* Changed tweedie.profile() to ignore values of p/xi outside (1, 2) rather than report an error.
* In some unusual cases, when p/xi = 0 was used (with add0 = TRUE), the mle of p was between 0 and 1 (which is impossible).  We report a warning message to check the data and the call to tweedie.profile(), but then set the mle to the value of p/xi giving the larger value of the likelihood
* Made the functions usually called by users able to accept xi or power
* A few minor edits to FORTRAN code; some variables not declared

tweedie 2.0.8 (Release date: 2011-06-08)
==============

Changes:

* Minor fixes to documentation
* Fixed the  add0  input
* Some minor changes to the code to tidy up adding the zero
* If values of xi/p are given between 0 and 1, or less than 0, they are now omitted (with a message) rather than creating an error. If there are no values left after omitting the problem value, an error message is given.

tweedie 2.0.7 (Release date: 2010-09-30)
==============

Changes:

* Ensured tweedie.profile() does not use power = 1.  This case (power=1 and phi not equal to 1) is too hard for me to deal with at present.
* Fixed an error introduced in version 2.0.5, where the value of xi.vec/p.vec was set to 1.2 (y >= 0) or 1.5 (y > 0) when not explicitly specified
* Fixed an error that reported the wrong mle of phi when the mle occurred at an endpoint of the given  xi  values.

tweedie 2.0.5 (Release date: 2010-08-27)
==============

Changes:

* Change to dtweedie.inversion() to ensure density = 0 is returned when y < 0 (1 < p < 2) or p <- 0 (when p > 2)
* Change to tweedie.profile() to fix a problem that p = 1 returned an error
* Changed so that tweedie.profile() works with p/xi = 0 when add0=TRUE (default is FALSE)
* Fixed some minor outputting messages (when verbose == 2)  
* Location of CITATION file moved to correct location

tweedie 2.0.4 (Release date: 2010-07-12)
==============

Changes:

* Minor edits

tweedie 2.0.3 (Release date: 2009-12-18)
==============

Changes:

* Changed default p.vec: There were too many values
* Added the facility to refer to  p  as  xi  in line with GLMs text

tweedie 2.0.1 (Release date: 2009-11-17)
==============

Changes:

* Changed the default p.vec when 1 < p < 2 to seq(1.2, 1.8, by = 0.05) (it was  seq(1.2, 1.8, by = 0.1)  )
* Slightly changed default output (added sep = "" to some paste commands)
* Added  AICtweedie()  to compute AIC for Tweedie glms

tweedie 2.0.0 (Release date: 2009-08-10)
==============

Changes:

* Made  method = "inversion" the default (was "series")
* Slightly changed default output (added sep = "" to some paste commands)
* An error introduced earlier (unsure when exactly; prob v 1.6.1)
* In trying to identify and fix the error, tidied some of the FORTRAN code
* Made  do.smooth = TRUE  the default (was FALSE)
* If  p.vec  is not supplied,  tweedie.profile()  makes a sensible guess
* Made  verbose = FALSE the default  (was TRUE)
* Made minor changes to output when  verbose = FALSE

tweedie 1.6.8 (Release date: 2009-07-18)
==============

Changes:

* Correct error in dtweedie.saddle() when y = 0 for 1 < p < 2 (Thanks to Glenn Meyers <GMeyers@iso.com>)

tweedie 1.6.7 (Release date: 2009--06-30)
==============

Changes:

* Correct error in rtweedie() when power = 1 and phi \ne 1 (Thanks to Frederic Gosselin <frederic.gosselin@cemagref.fr>)

tweedie 1.6.6 (Release date: 2009-09-19)
==============

Changes:

* Corrected specification of GPL; no changes to functionality, code or documentation
* In dtweedie(), added power == 3 (inverse Gaussian) as a special case
* Minor change in documentation
* Added  data, weights, and  offset  input arguments
* Fixed a bug where "=>" was used rather than ">="
* Fixed a bug in tweedie.profile(): exact zeros still allowed p = 2 (but not p > 2)

tweedie 1.6.2 (Release date: 2009-02-16)
==============

Changes:

* Fixed a bug in  ptweedie.inversion()  where variable  cdf  was undefined
* Fixed a bug in  qtweedie()  where any(power) < 1  should have been  any(power < 1)

tweedie 1.6.1 (Release date: 2009-02-06)
==============

Changes:

* Fixed error in tweedie.f()  (twice) where a float was used in a do-loop
* Fixed a bug in ptweedie(): It now returns 0 when y < 0 (rather than reporting an error)
* Fixed a bug in ptweedie.inversion(): y.len was undefined
* Fixed a bug in dtweedie.inversion(): if( any(phi < 0) ) written as if (any(phi) < 0 )
* Built for R version 2.6.0
* Changed the tweedie.plot() to produce no default x and y labels

tweedie 1.5.2 (Release date: 2007-09-01)
==============

Changes:

* Built for R version 2.5.1
* Removed the inversion computation from the example for  dtweedie(); it was very slow and hence caused problems

tweedie 1.5.1 (Release date: 2007-05-03)
==============

* Built for R version 2.5.0
* Fixed citation details
* Added note in dtweedie() help page indicating where the methods are defined (Dunn and Smyth, 2007).
* Added tweedie-package.Rd
* Fixed a small bug in  dtweedie.inversion() (replaced  p  with   power)
* Built for R version 2.5
* Made change to use of pmatch() in function  tweedie()
* Swapped the order of methods 1 and 2 in the call of tweedie.inversion  to make consistent with a change in the associated paper Dunn and Smyth.
* Fixed an uninitialized variable error in the FORTRAN code which sometimes caused the inversion procedure to fail.
* Other minor changes (fixing spelling errors, updated Dunn and Smyth reference; etc.)
* Some very minor changes (some cosmetic) in the code
* Made the code more robust, exiting gracefully when Inf and NA appear in the computed log-likelihood
* Slightly changed the output plot in tweedie.profile(), to reflect that some value of  p.vec  may not have a corresponding likelihood computed accurately (or, indeed, at all)
* Added new function   tweedie.plot()  to make it easier to plot Tweedie densities

tweedie 1.2 (Release date: 2006-09-05)
==============

Changes:

* Built for R version 2.0.0
* Major fix: The mle estimate of phi was incorrect in tweedie.profile() (thanks Sarah Lennox, QDPI&F) * Some better warnings in  tweedie.profile()
* tweedie.profile() now allows the actual computed points to be added to the plot using the  do.points  option (default is TRUE)
* The tweedie.profile() function now smooth and computes confidence intervals even when the computed likelihood values contain Infs or NAs.

tweedie 1.01 (Release date: 2004-09-30)
==============

Changes:

* Built for R version 1.9
* Fixed a bug in rtweedie(): it failed for p = 2 (thanks Sarah Lennox, QDPI&F)
* Fixed a bug in qtweedie(): it failed in some cases (thanks Gordon Smyth, WEHI)
* Added a proper  .First.lib  function (!) (thanks to James Wettenhall and Gordon Smyth, WEHI)
* Other minor bugs fixed
* The  tweedie()  family function removed, and package made dependent on the  statmod  package
* Updated for R version 1.9
