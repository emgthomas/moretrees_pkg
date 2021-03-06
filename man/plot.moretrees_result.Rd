% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.moretrees_result.R
\name{plot.moretrees_result}
\alias{plot.moretrees_result}
\title{\code{plot.moretrees_result} plots the groups discovered by MOReTreeS
on the outcome tree.}
\usage{
\method{plot}{moretrees_result}(x, group.text.size = 4,
  group.text.offset = 0.1, legend.text.size = 10,
  layout = "rectangular", horizontal = TRUE, ...)
}
\arguments{
\item{x}{Output from \code{moretrees()}
An object of class "moretrees_result".}

\item{group.text.size}{Text size for the group labels}

\item{group.text.offset}{Offset of the group label from the 
leaves of the tree}

\item{legend.text.size}{Text size for legend}

\item{layout}{Layout for the tree, most likely "rectangular" (the default)
or "slanted", but see the \code{layout} option of \code{ggtree()} for more
possibilities}

\item{horizontal}{If TRUE (the default), the tree will be plotted with
the root node at the top and all other nodes below. If FALSE, the tree
will be plotted with the root node to the left and all other nodes to
the right.}

\item{...}{Not used.}
}
\value{
A plot showing the groups discovered by MOReTreeS on the original
outcome tree.
}
\description{
\code{plot.moretrees_result} plots the groups discovered by MOReTreeS
on the outcome tree.
}
\examples{
# See vignette
vignette("moretrees")
}
\seealso{
Other MOReTrees results: \code{\link{print.moretrees_result}},
  \code{\link{print.summary.moretrees_compact}},
  \code{\link{print.summary.moretrees_long}},
  \code{\link{summary.moretrees_result}}
}
\concept{MOReTrees results}
