# msir 1.3.2 (2019-01 NOT ON CRAN)

* Update vignette.
* Added website for inclusion in mclust-org.
* Fix a bug in `msir.slices()`.
* Fix a bug on recycling array of length 1.
* **rgl** package is no longer a dependency but only suggested.

# msir 1.3.1 (2016-04)

* Fix issues related to `NAMESPACE`.

# msir 1.3 (2013-02)

* Removed **MASS** as suggested package since now it is not needed for examples in help pages.
* Fixed comments and notes from package check.

# msir 1.2 (2012-08)

* Compatibility bug fixes with new version 4.0 of **mclust** package.
* Add `verbose` argument to `summary.msir()`.

# msir 1.1 (2012-05)

* New `mvdnorm()` function which should be more stable in the case of near singularity of covariance matrix.
* Added byte-code compiling in the `DESCRIPTION` file.
* Corrected `msir.components()`, now it preserves the order of slices.
* Introduced the possibility to control EM algorithm.
* Bug correction and change in models' name in `msir.regularizedSigma()`.
* The `'msir'` object returns the normalized basis vectors in component `$basis` and raw eigenvectors in component `$evectors`.
* Update reference to CSDA paper.
* Bugs fix.

# msir 1.0 (2012-05)

* Initial release submitted to CRAN.
