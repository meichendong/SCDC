% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Deconvolution.R
\name{SCDC_basis_ONE}
\alias{SCDC_basis_ONE}
\title{Basis matrix for single cells from one subject}
\usage{
SCDC_basis_ONE(x, ct.sub = NULL, ct.varname, sample, ct.cell.size = NULL)
}
\arguments{
\item{x}{ExpressionSet object for single cells}

\item{ct.sub}{a subset of cell types that are selected to construct basis matrix}

\item{ct.varname}{variable name for 'cell types'}

\item{sample}{variable name for subject/samples}

\item{ct.cell.size}{default is NULL, which means the "library size" is calculated based on the data. Users can specify a vector of cell size factors corresponding to the ct.sub according to prior knowledge. The vector should be named: names(ct.cell.size input) should not be NULL.}
}
\value{
a list of basis matrix, sum of cell-type-specific library size, sample variance matrix, basis matrix by mvw, mvw matrix.
}
\description{
Basis matrix construction for single cells from one subject
}
