% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Deconvolution.R
\name{SCDC_prop}
\alias{SCDC_prop}
\title{Proportion estimation}
\usage{
SCDC_prop(
  bulk.eset,
  sc.eset,
  ct.varname,
  sample,
  ct.sub,
  iter.max = 1000,
  nu = 1e-04,
  epsilon = 0.01,
  truep = NULL,
  weight.basis = T,
  ct.cell.size = NULL,
  Transform_bisque = F,
  ...
)
}
\arguments{
\item{bulk.eset}{ExpressionSet object for bulk samples}

\item{sc.eset}{ExpressionSet object for single cell samples}

\item{ct.varname}{variable name for 'cell types'}

\item{sample}{variable name for subject/samples}

\item{ct.sub}{a subset of cell types that are selected to construct basis matrix}

\item{iter.max}{the maximum number of iteration in WNNLS}

\item{nu}{a small constant to facilitate the calculation of variance}

\item{epsilon}{a small constant number used for convergence criteria}

\item{truep}{true cell-type proportions for bulk samples if known}

\item{ct.cell.size}{default is NULL, which means the "library size" is calculated based on the data. Users can specify a vector of cell size factors corresponding to the ct.sub according to prior knowledge. The vector should be named: names(ct.cell.size input) should not be NULL.}

\item{Transform_bisque}{The bulk sample transformation from bisqueRNA. Aiming to reduce the systematic difference between single cells and bulk samples.}
}
\value{
Estimated proportion, basis matrix, predicted gene expression levels for bulk samples
}
\description{
Proportion estimation function for multi-subject case
}
