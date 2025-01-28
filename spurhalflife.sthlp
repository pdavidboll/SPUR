{smcl}
{* *! version 1.1.0 Jan 2025}{...}
{title:Title}

{pstd}
{hi:spurhalflife} {hline 2} Confidence Sets for the Half-Life of Spatial Processes


{title:Syntax}

{p 8 16 2}
{cmd:spurhalflife} {varname} {ifin} [{cmd:,} {it:options}]


{synoptset 11}{...}
{synopthdr}
{synoptline}
{synopt :{opt q(#)}}specifies the number of weighted averages to be used in the test. The default is {opt q(15)}.{p_end}
{synopt :{opt nrep(#)}}specifies the number of Monte Carlo draws to be used to simulate the distribution of the test statistic. The default is {opt nrep(100000)}.{p_end}
{synopt :{opt l:evel(#)}} specifies the desired confidence level in percent. The default is {opt level(95)}.{p_end}
{synopt :{opt latlong}}specifies that the spatial coordinates are given in latitude (stored in {it:s_1}) and longitude (stored in {it:s_2}) (see below).{p_end}
{synopt :{opt norm:dist}}specifies that the results are to be returned as fractions of the maximum pairwise distance in the sample.{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
This command implements a method to construct confidence sets about the half-life of spatial processes developed by {help spurtest##mw2024:Müller and Watson (2024)}. This inference is based on the persistence parameter in the Spatial Local-to-Unity model, which maps directly into the half-life (i.e., the distance at which the correlation is 1/2).

{pstd}
{cmd:spurhalflife} is part of a package of commands that also includes {cmd:spurtest} and {cmd:spurtransform}. A practical guide to these methods is provided in a paper accompanying this implementation ({help spurtest##bbv2025:Becker et. al (2025)}), please cite this paper when using this code.

{pstd}
Note: For this command (and all other commands in this package) to work, the spatial coordinates must be stored in the variables {it:s_*}, where * is a positive integer. This is for consistency with the {cmd:scpc} command developed by {help spurtest##mw2022:Müller and Watson (2022, 2023)} (available from Ulrich Müller's website), which this package is designed to work alongside. 

{phang2}
If the option {opt latlong} is specified, {it:s_1} is interpreted as latitude and {it:s_2} as longitude, and no other {it:s_*} variables may be present. Latitude and longitude are expected in decimal degrees (i.e., lat lies in [-90, 90] and lon in [-180, 180]).

{phang2}
If the option is not specified, the {it:p} {it:s_*} variables present are interpreted as coordinates in {it:p}−dimensional Euclidean space. This is appropriate for projected coordinates.

{pstd}
Note: This command and all others in this package rely on the {cmd: moremata} package by Ben Jann, which can be installed from SSC: {stata "ssc install moremata, replace":auto-install moremata}.


{marker options}{...}
{title:Options}

{phang}
{opt q(#)} specifies the number of weighted averages to be used in the test. The default is {opt q(15)}. See {help spurtest##mw2024:Müller and Watson (2024)} and {help spurtest##bbv2025:Becker et. al (2025)} for details.

{phang}
{opt nrep(#)} specifies the number of Monte Carlo draws to be used to simulate the distribution of the test statistic. The default is {opt nrep(100000)}. See {help spurtest##mw2024:Müller and Watson (2024)} and {help spurtest##bbv2025:Becker et. al (2025)} for details.

{phang}
{opt l:evel(#)} specifies the desired confidence level in percent. The default is {opt level(95)}.

{phang}
{opt latlong} specifies that the spatial coordinates are given in latitude (stored in {it:s_1}) and longitude (stored in {it:s_2}) (see above).

{phang}
{opt norm:dist} specifies that the results are to be returned as fractions of the maximum pairwise distance in the sample. Otherwise, they are returned in meters (if {opt latlong}) or the units of the original Euclidean coordinates (if not {opt latlong}).

{marker examples}{...}
{title:Examples}

{phang}{cmd:. sysuse auto}{p_end}
{phang}{cmd:. gen s_1=rnormal(0,1)} (generate random coordinates){p_end}
{phang}{cmd:. gen s_2=rnormal(0,1)}{p_end}
{phang}{cmd:. spurhalflife mpg}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:spurhalflife} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 23 26 2: Scalars}{p_end}
{synopt:{cmd:r(ci_l)}}Lower bound of the confidence interval{p_end}
{synopt:{cmd:r(ci_U)}}Upper bound of the confidence interval{p_end}
{synopt:{cmd:r(max_dist)}}Maximum pairwise distance in the sample{p_end}
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