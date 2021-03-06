% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spike_and_slab_logistic_moretrees.R
\name{spike_and_slab_logistic_moretrees}
\alias{spike_and_slab_logistic_moretrees}
\title{Description}
\usage{
spike_and_slab_logistic_moretrees(dsgn, vi_params_init, hyperparams_init,
  random_init, random_init_vals, tol, tol_hyper, max_iter, print_freq,
  update_hyper_freq, hyper_fixed)
}
\arguments{
\item{dsgn}{list of outputs from \code{moretrees_design_tree()}}

\item{vi_params_init, hyperparams_init, random_init, random_init_vals, tol, tol_hyper, max_iter, print_freq, update_hyper_freq, hyper_fixed}{see documentation for \code{moretrees()}}
}
\value{
A named list containing the following entried:
\describe{
\item{\code{vi_params}}{named list of final variational parameter estimates}
\item{\code{hyperparams}}{named list of final hyperparameter estimates}
\item{\code{hyper_fixed}}{named list of fixed hyperparameters}
\item{\code{ELBO_track}}{numeric vector containing the values of the objective function
(ELBO) at the end of every iteration}
}
}
\description{
\code{spike_and_slab_logistic_moretrees()} fits MOReTreeS models for
  binary data.
  The posterior is approximated via variational inference.
  This function returns the parameters of the variational approximation.
}
\seealso{
Other Internal VI functions: \code{\link{moretrees_init_logistic}},
  \code{\link{update_hyperparams_logistic_moretrees}},
  \code{\link{update_vi_params_logistic_moretrees}}
}
\concept{Internal VI functions}
