# `spur`: A Stata Package around Spatial Unit Roots

This package implements methods for diagnosing and correcting spatial unit roots developed by Müller and Watson (2024). 

**When using this code, please cite [Becker, Boll and Voth (2025) ??????]().**

## Installation

Installation from GitHub:

    . net install spur, replace from(https://raw.githubusercontent.com/pdavidboll/spur/main/)

If you *also* want to download the example data and do-file (will be saved to the current directory), run the following *in addition* to the above:

    . net get spur, replace from(https://raw.githubusercontent.com/pdavidboll/spur/main/)

Then, you can run the example do-file, which reproduces Table 1 from [Müller and Watson (2024)](https://www.princeton.edu/~umueller/SPUR.pdf) based on the data from [Chetty et. al. (2014)](https://doi.org/10.1093/qje/qju022) (see [Becker, Boll and Voth (2025) ??????]() for details):

    . do example

## References

Becker, Sascha O., P. David Boll and Hans-Joachim Voth "Spatial Unit Roots in Regressions: A Practitioner's Guide and a Stata Package", 2025.

Müller, Ulrich K. and Mark W. Watson "Spatial Unit Roots and Spurious Regression", Econometrica 92 (2024), 1661–1695. https://www.princeton.edu/~umueller/SPUR.pdf.

