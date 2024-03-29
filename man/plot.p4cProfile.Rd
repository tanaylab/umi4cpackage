% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/p4c_package.R
\name{plot.p4cProfile}
\alias{plot.p4cProfile}
\alias{plotSingleProf}
\alias{plotCompProf}
\title{Plot 4C nearcis profile.}
\usage{
\method{plot}{p4cProfile}(p4c_obj1, p4c_obj2, trend_scale = "adaptive", png_fn = NULL, ...)

plotSingleProf(
  p4c_obj,
  trend_scale,
  png_fn,
  plot.colorbar,
  add.func,
  xlim,
  ylim,
  trend_only,
  main,
  sd,
  ...
)

plotCompProf(
  p4c_obj1,
  ref_p4c_obj,
  trend_scale,
  png_fn,
  col,
  min_win_cov,
  xlim,
  zlim,
  legend.text,
  ylim,
  dgram.method,
  main,
  sd,
  filename1,
  filename2,
  trend_only
)
}
\arguments{
\item{p4c_obj1}{p4cProfile object created by \code{p4cNewProfile}.}

\item{p4c_obj2}{Optional. A second p4cProfile object that will be compared to p4c_obj1.}

\item{trend_scale}{This parameter controls the smoothing method for the trendline. Either
'adaptive' (default), or, for fixed window trend, an integer with the desired size
of the window (in restriction fragments units).}

\item{png_fn}{File to save the plot to. The default (NULL) will result in plotting on the
current graphic device.}

\item{plot.colorbar}{Optional. Logical that defines whether a colorbar for the domainogram
should be plotted.}

\item{xlim}{Similar to R basic graphics. A two elements vector that controls the genomic
focus of the plot (X axis).}

\item{ylim}{Similar to R basic graphics. A two elements vector that controls the trendline
ylim (In units of molecules).}

\item{trend_only}{Optional. Logical that defines whether only the trend plot should be displayed.}

\item{main}{Plot title.}

\item{sd}{Controls the margins of the trendline. The Number of bionomial theoretical standard deviations.}

\item{col}{Vector with length 2 that controls the colors of each profile in a comparative plot.}

\item{min_win_cov}{If the default smoothing method is used ('adaptive'), this parameter
controls the window size by requiring that no less than \code{min_win_cov} molecules will be included
in that genomic window.}

\item{legend.text}{Character vector with length 2. Controls the text in the color legend of the trendline.
The default are the names of the profiles.}

\item{filename1}{Optional. File for the bedgraph with the trendline of `p4c_obj1`}

\item{filename2}{Optional. File for the bedgraph with the trendline of `p4c_obj2`}
}
\description{
Plot a nearcis profile in the defined genomic scope. Can be applied on a single
p4cProfile object or a two object for a comparison plot.
}
\details{
This function calls two methods depending on whether it is called on a single object
or two objects. When it is called on a single p4cProfile object it will produce a nearcis
plot with smoothed trendline, and a contact intensity domainogram.
Calling this function for two object will produce a comparison plot of where the
4C profiles are normalized for UMI coverage. The two p4cProfile objects must be
constructed for the \strong{same bait} with the \strong{same genomic scope}.
}
\examples{
\donttest{
# Create a p4cProfile object:
fc <- p4cNewProfile("umi4C_example_CMK_ANK1_TSS")

# Plot a profile on:
plot(fc, xlim = c(41554693, 41754693), ylim = c(0, 5), plot.colorbar = TRUE)


# Comparing two profiles:
# Create a second profile
fc_cond2 <- p4cNewProfile("umi4C_example_DND41_ANK1_TSS")

# Plot profiles
plot(fc, fc_cond2, col = c("red", "blue"), min_win_cov = 100, main = "4C figure", png_fn = "fig.png")
}

}
