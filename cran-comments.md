## Resubmission

I was told "package tweedie_2.3.1.tar.gz does not pass the incoming checks automatically"

The URL to which I was pointed has no errors or warnings showing, just one NOTE.  Hopefully this fixes.

## Test environments
* local OS X install, R 3.4.2
* ubuntu 12.04 (on travis-ci), R 3.4.2
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

* I have run R CMD check on the downstream dependencies.
  
Checked ChainLadder: 0 errors | 0 warnings | 0 notes
Checked cplm       : 0 errors | 0 warnings | 0 notes
Checked dsm        : 0 errors | 0 warnings | 1 note 
Checked HeritSeq   : 1 error  | 0 warnings | 0 notes
Checked mcglm      : 0 errors | 0 warnings | 0 notes
Checked mvabund    : 1 error  | 0 warnings | 0 notes
Checked raw        : 0 errors | 0 warnings | 0 notes
Checked statmod    : 0 errors | 0 warnings | 0 notes

The errors are because (other) packages are not available, not to do with tweedie

* revdep maintainers were not notified of the releases, as the change do not impact them.
