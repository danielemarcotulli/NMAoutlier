NMAoutlier: Detecting Outliers in Network Meta-Analysis
================

[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/NMAoutlier)](https://cran.r-project.org/package=NMAoutlier)
[![CRAN\_time\_from\_release](https://www.r-pkg.org/badges/ago/NMAoutlier)](https://cran.r-project.org/package=NMAoutlier)
[![Monthly
Downloads](http://cranlogs.r-pkg.org/badges/NMAoutlier)](http://cranlogs.r-pkg.org/badges/NMAoutlier)
[![Total
Downloads](http://cranlogs.r-pkg.org/badges/grand-total/NMAoutlier)](http://cranlogs.r-pkg.org/badges/grand-total/NMAoutlier)
[![Build
Status](https://travis-ci.org/petropouloumaria/NMAoutlier.svg?branch=master)](https://travis-ci.org/petropouloumaria/NMAoutlier)

<img src="man/figures/NMAoutlier_logo.png" width=300 align="right" style="margin-left:20px; margin-right: 20px;"/>

## Description

A package that provides measures and methodologies for detecting
outlying and influential studies in network meta-analysis.

-   **1) Simply outlier and influential detection measures:** Raw,
    Standardized, Studentized residuals; Mahalanobis distance and
    leverage.
-   **2) Outlier and influential detection measures by considering a
    study deletion (Shift the mean):** Raw, Standardized, Studentized
    deleted residuals; Cook distance; COVRATIO; weight “leave one out”;
    leverage “leave one out”; heterogeneity “leave one out”; R
    heterogeneity; R Qtotal; R Qheterogeneity; R Qinconsistency and
    DFBETAS.
-   Plots for all the above outlier and influential detection measures
    (simple and deletion  
    measures) and Q-Q plot for network meta-analysis.
-   **3) Forward search algorithm in network meta-analysis (FS).**
-   Forward plots (fwdplot) for the monitoring measures in each step of
    forward search algorithm. Monitoring measures: P-scores; z-values
    for difference of direct and indirect evidence with back-calculation
    method; Standardized residuals; heterogeneity variance estimator;
    cook distance; ratio of variances; Q statistics.
-   Forward plot for summary estimates and their confidence intervals
    for each treatment in each step of forward search algorithm.  
-   **4) Random shift variance NMA model (RSVM NMA) (Shift the study
    variance).**
-   Plots for the monitoring measures for random shift variance model.
    Monitoring measures: P-scores; z-values for difference of for
    difference of direct and indirect evidence with back-calculation
    method; Standardized residuals; heterogeneity variance estimator;
    over-dispersion parameter; leverage; Q statistics; Likelihood Ratio
    Test (LRT).
-   Plots for the for summary estimates and their confidence intervals
    for each treatment estimated with the random shift variance model.

## Installation

You can install the **NMAoutlier** package from GitHub repository as
follows:

Installation using R package
**[devtools](https://cran.r-project.org/package=devtools)**:

``` r
install.packages("devtools")
devtools::install_github("petropouloumaria/NMAoutlier")
```

## Usage

Example of network meta-analysis comparing the relative effects of four
smoking cessation counseling programs, no contact (A), self-help (B),
individual counseling (C) and group counseling (D). The outcome is the
number of individuals with successful smoking cessation at 6 to 12
months. The data are in contrast format with odds ratio (OR) and its
standard error. Arm-level data can be found in Dias et al. (2013).

Reference:

Higgins D, Jackson JK, Barrett G, Lu G, Ades AE, and White IR. (2012):
Consistency and inconsistency in network meta-analysis: concepts and
models for multi-arm studies. Research Synthesis Methods 3(2):
98–110.

Dias S, Welton NJ, Sutton AJ, Caldwell DM, Lu G, and Ades AE. (2013): Evidence
Synthesis for Decision Making 4: Inconsistency in networks of evidence
based on randomized controlled trials. Medical Decision Making, 33:
641–656.

You can load the **NMAoutlier** library

``` r
library(NMAoutlier)
```

Load the dataset smoking cessation from **netmeta** package.

``` r
data(smokingcessation, package = "netmeta")
```

Transform data from arm-based to contrast-based format using the
function **pairwise** from **netmeta** package.

``` r
library(netmeta)
p1 <- pairwise(list(treat1, treat2, treat3),
              list(event1, event2, event3),
              list(n1, n2, n3),              
              data=smokingcessation,
              sm="OR")
```

**Part 1: Simply outlier detection measures**

You can calculate simply outlier and influential detection measures with
**NMAoutlier.measures** function as follows:

``` r
measures <- NMAoutlier.measures(p1)
```

You can see the Mahalanobis distance for each study

``` r
measures$Mahalanobis.distance
```

You can plot the Mahalanobis distance for each study with **measplot**
function as follows:

``` r
measplot(measures, "mah")
```

You can figure out the Q-Q plot for network meta-analysis with
**Qnetplot** function as follows:

``` r
Qnetplot(measures)
```

**Part 2: Outlier detection measures considered deletion (Shift the
mean)**

You can calculate outlier and influential detection measures considered
study deletion with **NMAoutlier,measures** function as follows:

``` r
deletion <- NMAoutlier,measures(p1, measure = "deletion")
```

You can see the standardized deleted residuals for each study

``` r
deletion$estand.deleted
```

You can see the COVRATIO for each study

``` r
deletion$Covratio
```

You can plot the R statistic for Qinconsistency with function
**measplot** as follows:

``` r
measplot(deletion, "rqinc", measure = "deletion")
```

**Part 3: Forward Search Algorithm - (Outlier detection Methodology)**

You can conduct the Forward Search algorithm with **NMAoutlier**
function as follows:

``` r
FSresult <- NMAoutlier(p1, small.values = "bad")
```

You can see the forward plots with **fwdplot** function for Cook
distance as follows:

``` r
fwdplot(FSresult,"cook")
```

<img src="man/figures/fwdplot-cook-distance.png" width=450 style="margin-left: auto; margin-right: auto; display: block;"/>

Or you can plot the Ratio of variances as follows:

``` r
fwdplot(FSresult,"ratio")
```

<img src="man/figures/fwdplot-ratio.png" width=450 style="margin-left: auto; margin-right: auto; display: block;"/>

You can plot the differences of direct and indirect estimates (z-values)
as follows:

``` r
fwdplot(FSresult,"nsplit")
```

<img src="man/figures/fwdplot-diff-direct-indirect.png" width=700 style="margin-left: auto; margin-right: auto; display: block;"/>

You can see the forward plots for summary relative treatment estimates
of B, C and D versus the reference A with **fwdplotest** function as
follows:

``` r
fwdplotest(FSresult)
```

<img src="man/figures/fwdplot-summary-estimates.png" style="margin-left: auto; margin-right: auto; display: block;"/>

**Part 4: Random Shift Variance Network Meta-analysis (RSV NMA) (Shift
the study variance)/Outlier detection methodology and sensitivity
analysis by downweighting outlier**

You can conduct the Random Shift Variance RSV NMA model by shifting the
variance for each study with **NMAoutlier.rsv** function as follows:

``` r
RSVresult <- NMAoutlier.rsv(p1, small.values = "bad")
```

You can see the Likelihood Ratio Test (LRT) with RSV NMA model for each
study

``` r
RSVresult$LRT 
```

Or you can see the over-dispersion parameter of RSV NMA model for each
study

``` r
RSVresult$over_disp
```

You can plot the Likelihood Ratio Test (LRT) of RSV NMA model for each
study with **rsvplot** function as follows:

``` r
rsvplot(RSVresult, "LRT")
```

You can see the plots for summary relative treatment estimates for B, C
and D versus the reference A with **rsvplotest** function as follows:

``` r
rsvplotest(RSVresult)
```
