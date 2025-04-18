% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NMAoutlier-package.R
\docType{package}
\name{NMAoutlier-package}
\alias{NMAoutlier-package}
\title{NMAoutlier: Brief overview of measures and methodologies for
detection of outlying and influential studies in network
meta-analysis.}
\description{
R package \bold{NMAoutlier} provides methods and tools to detect
outlier and influential studies in network meta-analysis.
}
\details{
R package \bold{NMAoutlier} is a tool to detect outliers (studies
with extreme results) and influential studies in network
meta-analysis (Petropoulou, 2020). The package can calculate:
simple outlier and influential measures; outlier and influential
measures considered study deletion (shift the mean); the outlier
detection methodology with Forward Search (FS) algorithm
(Petropoulou et al., 2021). All proposed outlier and influential
detection methods were fitted the frequentist NMA model by graph
theory introduced by Rücker (2012) and implemented inR package
\bold{netmeta}.

The \bold{NMAoutlier} package implements the following methods
described in Petropoulou (2020).
\itemize{
 \item \bold{Simple outlier and influential detection measures}
   (function \code{\link{NMAoutlier.measures}}):
 \enumerate{
   \item raw residuals,
   \item standardized residuals,
   \item studentized residuals,
   \item Mahalanobis distance,
   \item leverage;
  }
\item \bold{Outlier and influential detection measures considered
 study deletion (shift the mean)} (function
 \code{\link{NMAoutlier.measures}}):
 \enumerate{
   \item raw deleted residuals,
   \item standardized deleted residuals,
   \item studentized deleted residuals,
   \item Cook's distance,
   \item COVRATIO,
   \item weight leave one out,
   \item leverage leave one out,
   \item heterogeneity leave one out,
   \item R heterogeneity,
   \item R Qtotal,
   \item R Qheterogeneity,
   \item R Qinconsistency,
   \item DFBETAS;
 }
\item Plots of the several outlier and influential detection
  (simple and deletion) measures (function
  (\code{\link{measplot}}));
\item Q-Q plot for network meta-analysis (function
  \code{\link{Qnetplot}});
\item \bold{Forward Search algorithm in network meta-analysis}
  (function (\code{\link{NMAoutlier}})) based on Petropoulou et
  al. (2021);
\item forward plots (\code{\link{fwdplot}}) with monitoring
  statistics in each step of the FS algorithm:
\enumerate{
  \item P-scores (Rücker & Schwarzer, 2015),
  \item z-values for difference of direct and indirect evidence
    with back-calculation method (König et al., 2013; Dias et al.,
    2010),
  \item standardized residuals,
  \item heterogeneity variance estimator,
  \item Cook's distance,
  \item ratio of variances,
  \item Q statistics (Krahn et al., 2013);
  }
\item forward plots (\code{\link{fwdplotest}}) for summary
  treatment estimates in each iteration of the FS algorithm
  (Petropoulou et al., 2021).
}

Type \code{help(package = "NMAoutlier")} for a listing of R functions
available in \bold{NMAoutlier}.

Type \code{citation("NMAoutlier")} on how to cite \bold{NMAoutlier}
in publications.

To report problems and bugs, please send an email to Dr. Maria
Petropoulou \email{maria.petropoulou@uniklinik-freiburg.de}.

The development version of \bold{NMAoutlier} is available on GitHub
\url{https://github.com/petropouloumaria/NMAoutlier}.
}
\references{
Dias S, Welton NJ, Caldwell DM, Ades AE (2010):
Checking consistency in mixed treatment comparison meta-analysis.
\emph{Statistics in Medicine},
\bold{29}, 932--44

König J, Krahn U, Binder H (2013):
Visualizing the flow of evidence in network meta-analysis and
characterizing mixed treatment comparisons.
\emph{Statistics in Medicine},
\bold{32}, 5414--29

Krahn U, Binder H, König J (2013):
A graphical tool for locating inconsistency in network meta-analyses.
\emph{BMC Medical Research Methodology},
\bold{13}, 35

Petropoulou M (2020):
Exploring methodological challenges in network meta-analysis models
and developing methodology for outlier detection.
\emph{PhD dissertation}

Petropoulou M, Salanti G, Rücker G, Schwarzer G, Moustaki I,
Mavridis D (2021):
A forward search algorithm for detecting extreme study effects in
network meta-analysis.
\emph{Statistics in Medicine}

Rücker G (2012):
Network meta-analysis, electrical networks and graph theory.
\emph{Research Synthesis Methods},
\bold{3}, 312--24

Rücker G, Schwarzer G (2015):
Ranking treatments in frequentist network meta-analysis works
without resampling methods.
\emph{BMC Medical Research Methodology},
\bold{15}, 58
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/petropouloumaria/NMAoutlier}
}

}
\author{
Petropoulou Maria \email{maria.petropoulou@uniklinik-freiburg.de}.
}
