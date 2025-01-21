{smcl}
{* *! version 1.1.0 Jan 2025}{...}
{title:Title}

{pstd}
{hi:spurtest} {hline 2} Diagnostic tests for Spatial Unit Roots


{title:Syntax}

{phang}
Spatial I(1) test

{p 8 16 2}
{cmd:spurtest i1} {varname} {ifin} [{cmd:,} {it:options}]


{p 4 8 2}
Spatial I(0) test

{p 8 16 2}
{cmd:spurtest i0} {varname} {ifin} [{cmd:,} {it:options}]


{p 4 8 2}
Spatial I(1) test for residuals

{p 8 16 2}
{cmd:spurtest i1resid} {depvar} [{indepvars}] {ifin} [{cmd:,} {it:options}]


{p 4 8 2}
Spatial I(0) test for residuals

{p 8 16 2}
{cmd:spurtest i0resid} {depvar} [{indepvars}] {ifin} [{cmd:,} {it:options}]


{synoptset 11}{...}
{synopthdr}
{synoptline}
{synopt :{opt q(#)}}specifies the number of weighted averages to be used in the test. The default is {opt q(15)}.{p_end}
{synopt :{opt nrep(#)}}specifies the number of Monte Carlo draws to be used to simulate the distribution of the test statistic. The default is {opt nrep(100000)}.{p_end}
{synopt :{opt latlong}}specifies that the spatial coordinates are given in latitude (stored in {it:s_1}) and longitude (stored in {it:s_2}) (see below).{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
This command implements four tests developed by {help spurtest##mw2024:Müller and Watson (2024)} to diagnose unit roots in spatial processes. {cmd:spurtest i1resid} and {cmd:spurtest i0resid} test the null hypotheses of I(1) (unit root) and I(0) (stationarity), respectively, for an observed outcome variable {it:varname}. {cmd:spurtest i1} and {cmd:spurtest i0} test these hypotheses for the errors of a linear regression of {it:depvar} on {it:indepvars} (a constant is automatically included). 

{pstd}
{cmd:spurtest} is part of a package of commands that also includes {cmd:spurtransform} and {cmd:spurhalflife}. A practical guide to these methods is provided in a paper accompanying this implementation ({help spurtest##bbv2025:Becker et. al (2025)}), please cite this paper when using this code.

{pstd}
Note: For this command (and all other commands in this package) to work, the spatial coordinates must be stored in the variables {it:s_*}, where * is a positive integer. This is for consistency with the {cmd:scpc} command developed by {help spurtest##mw2022:Müller and Watson (2022, 2023)} (available from Ulrich Müller's website), which this package is designed to work alongside. If the option latlong is specified, {it:s_1} is interpreted as latitude and {it:s_2} as longitude, and no other {it:s_*} variables may be present. If the option is not
specified, the {it:p} {it:s_*} variables present are interpreted as coordinates in {it:p}−dimensional Euclidean space.


{marker options}{...}
{title:Options}

{phang}
{opt q(#)} specifies the number of weighted averages to be used in the test. The default is {opt q(15)}. See {help spurtest##mw2024:Müller and Watson (2024)} and {help spurtest##bbv2025:Becker et. al (2025)} for details.

{phang}
{opt nrep(#)} specifies the number of Monte Carlo draws to be used to simulate the distribution of the test statistic. The default is {opt nrep(100000)}. See {help spurtest##mw2024:Müller and Watson (2024)} and {help spurtest##bbv2025:Becker et. al (2025)} for details.

{phang}
{opt latlong} specifies that the spatial coordinates are given in latitude (stored in {it:s_1}) and longitude (stored in {it:s_2}) (see above).

{marker examples}{...}
{title:Examples}

{phang}{cmd:. sysuse auto}{p_end}
{phang}{cmd:. gen s_1=rnormal(0,1)} (generate random coordinates){p_end}
{phang}{cmd:. gen s_2=rnormal(0,1)}{p_end}
{phang}{cmd:. spurtest i1 mpg}{p_end}
{phang}{cmd:. spurtest i1 mpg, q(15) nrep(100000)}{space 6}(equivalent to above command){p_end}
{phang}{cmd:. spurtest i0resid mpg weight length}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:spurtest} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 23 26 2: Scalars}{p_end}
{synopt:{cmd:r(teststat)}}Test statistic (LFUR or LFST){p_end}
{synopt:{cmd:r(p)}}P-value of the test{p_end}
{synopt:{cmd:r(ha_param)}}Parameter for alternative hypothesis (c_a or g_a). See {help spurtest##mw2024:Müller and Watson (2024)} and {help spurtest##bbv2025:Becker et. al (2025)} for details. {p_end}

{p2col 5 15 19 2: Matrices}{p_end}
{synopt:{cmd:r(cv)}}Critical values at 1%, 5%, and 10% levels.{p_end}

{p2colreset}{...}


{marker authors}{...}
{title:Authors}

{pstd}
Sascha O. Becker{break}
University of Warwick{break}
Coventry, UK{break}
s.o.becker@warwick.ac.uk{break}

{pstd}
P. David Boll{break}
University of Warwick{break}
Coventry, UK{break}
david.boll@warwick.ac.uk{break}

{pstd}
Hans-Joachim Voth{break}
University of Zurich{break}
Zurich, Switzerland{break}
voth@econ.uzh.ch{break}

{pstd}
This command and all others in this package were written based on the Matlab replication files for {help spurtest##mw2024:Müller and Watson (2024)}, available from the authors' websites. We thank Ulrich Müller and Mark Watson for their support and useful conversations.


{marker disclaimer}{...}
{title:Disclaimer}

{pstd}
This software is provided "as is", without warranty of any kind.
If you have suggestions or want to report problems, please create a new issue in the {browse "https://github.com/pdavidboll/SPUR/issues":project repository} or contact the project maintainer.


{marker references}{...}
{title:References}

{marker bbv2025}{...}
{phang}Becker, Sascha O., P. David Boll and Hans-Joachim Voth "Spatial Unit Roots in Regressions: A Practitioner's Guide and a Stata Package", 2025.

{marker mw2022}{...}
{phang}Müller, Ulrich K. and Mark W. Watson "Spatial Correlation Robust Inference", Econometrica 90 (2022), 2901–2935. {browse "https://www.princeton.edu/~umueller/SHAR.pdf"}.

{marker mw2023}{...}
{phang}Müller, Ulrich K. and Mark W. Watson "Spatial Correlation Robust Inference in Linear Regression and Panel Models", Journal of Business and Economic Statistics 41 (2023), 1050–1064. {browse "https://www.princeton.edu/~umueller/SpatialRegression.pdf"}.

{marker mw2024}{...}
{phang}Müller, Ulrich K. and Mark W. Watson "Spatial Unit Roots and Spurious Regression", Econometrica 92 (2024), 1661–1695. {browse "https://www.princeton.edu/~umueller/SPUR.pdf"}.
{p_end}