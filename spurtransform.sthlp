{smcl}
{* *! version 1.1.0 Jan 2025}{...}
{title:Title}

{pstd}
{hi:spurtransform} {hline 2} Differencing Transformations to Eliminate Spatial Unit Roots


{title:Syntax}

{p 8 16 2}
{cmd:spurtransform} {varlist} {ifin} {cmd:,} {opt pre:fix(string)} [{it: options}]


{synoptset 22}{...}
{synopthdr}
{synoptline}
{synopt :{opt pre:fix(string)}}specifies the prefix for the variable names under which the transformed data will be stored.{p_end}
{synopt :{opt transformation(string)}}specifies the type of transformation. Must be one of {opt nn}, {opt iso}, {opt cluster}, {opt lbmgls}. Defaults to {opt lbmgls}.{p_end}
{synopt :{opt radius(#)}}specifies the radius to be used for isotropic differencing. Only allowed with {opt transformation(iso)}.{p_end}
{synopt :{opt clustvar(varname)}}specifies the variable that is to be used for clustering. Only allowed with {opt transformation(cluster)}.{p_end}
{synopt :{opt latlong}}specifies that the spatial coordinates are given in latitude (stored in {it:s_1}) and longitude (stored in {it:s_2}) (see below).{p_end}
{synopt :{opt r:eplace}}allows the command to overwrite variables when storing the transformed data.{p_end}
{synopt :{opt separately}}executes the transformation separately for all variables in {it: varlist}.{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
This command implements four spatial differencing transformations proposed by {help spurtest##mw2024:Müller and Watson (2024)} to eliminate unit roots in spatial processes. By default, it applies the LBM-GLS transformation, favoured by {help spurtest##mw2024:Müller and Watson (2024)}. Other transformations can be selected using the {opt transformation()} option.

{pstd}
{cmd:spurtransform} is part of a package of commands that also includes {cmd:spurtest} and {cmd:spurhalflife}. A practical guide to these methods is provided in a paper accompanying this implementation ({help spurtest##bbv2025:Becker et. al (2025)}), please cite this paper when using this code.

{pstd}
Note: For this command (and all other commands in this package) to work, the spatial coordinates must be stored in the variables {it:s_*}, where * is a positive integer. This is for consistency with the {cmd:scpc} command developed by {help spurtest##mw2022:Müller and Watson (2022, 2023)} (available from Ulrich Müller's website), which this package is designed to work alongside. If the option latlong is specified, {it:s_1} is interpreted as latitude and {it:s_2} as longitude, and no other {it:s_*} variables may be present. If the option is not
specified, the {it:p} {it:s_*} variables present are interpreted as coordinates in {it:p}−dimensional Euclidean space.

{pstd}
Note: This command and all others in this package rely on the {cmd: moremata} package by Ben Jann, which can be installed from SSC: {stata "ssc install moremata, replace":auto-install moremata}.


{marker options}{...}
{title:Options}

{phang}
{opt pre:fix(string)} specifies the prefix for the variable names under which the transformed data will be stored. For example, when transforming the variable {it:x}, specifying {opt prefix("h_")} will result in the transformed variable to be stored as {it:h_x}. 

{phang}
{opt transformation(string)} specifies the type of transformation. Must be one of {opt nn}, {opt iso}, {opt cluster}, {opt lbmgls}. Defaults to {opt lbmgls}.

{phang2}
{opt transformation(nn)} requests the Nearest-Neighbor transformation, i.e. differencing with respect to each observation's nearest neighbor.

{phang2}
{opt transformation(iso)} requests the Isotropic transformation, i.e. differencing with respect to an unweighted average of all (other) observations within a {opt radius} meter radius.

{phang2}
{opt transformation(cluster)} requests the Cluster transformation, i.e. unweighted demeaning within a cluster of observations, where the cluster variable is specified by {opt clustvar}.

{phang2}
{opt transformation(lbmgls)} requests the LBM-GLS transformation, i.e. GLS transformation based on the covariance matrix of a Lévy-Brownian motion. This is the default option.

{phang2}
See {help spurtest##mw2024:Müller and Watson (2024)} and {help spurtest##bbv2025:Becker et. al (2025)} for details.

{phang}
{opt radius(#)} specifies the radius in metres (if {opt latlong}), or in the units of the original coordinates (if not {opt latlong}), which is to be used for isotropic differencing. Only allowed with {opt transformation(iso)}.

{phang}
{opt cluster(varname)} specifies the variable that is to be used for clustering. Only allowed with {opt transformation(cluster)}.

{phang}
{opt latlong} specifies that the spatial coordinates are given in latitude (stored in {it:s_1}) and longitude (stored in {it:s_2}) (see above).

{phang}
{opt r:eplace} allows the command to overwrite variables when storing the transformed data.

{phang}
{opt separately} executes the transformation separately for all variables in varlist. This leads to different results if there are missing observations in some variables, because the default behavior is to construct the transformation matrix based only on those observations for which all variables are non-missing.


{marker examples}{...}
{title:Examples}

{phang}{cmd:. sysuse auto}{p_end}
{phang}{cmd:. gen s_1=rnormal(0,1)} (generate random coordinates){p_end}
{phang}{cmd:. gen s_2=rnormal(0,1)}{p_end}
{phang}{cmd:. spurtransform mpg, prefix("h_")}{p_end}
{phang}{cmd:. spurtransform mpg, prefix("h_") transformation(lbmgls)}{space 6}(equivalent to above command){p_end}
{phang}{cmd:. spurtransform mpg weight length, prefix("h_")}{space 6}(applies transformation to all three variables){p_end}
{phang}{cmd:. spurtransform mpg, prefix("h_") transformation(iso) radius(0.2)}{space 6}(applies isotropic transformation within radius of 0.2){p_end}


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