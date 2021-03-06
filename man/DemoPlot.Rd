% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Visualization.R
\name{DemoPlot}
\alias{DemoPlot}
\title{Demographic figure for single cells}
\usage{
DemoPlot(
  eset,
  cluster,
  sample,
  select.ct,
  Palette = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
    "#D55E00", "#CC79A7")
)
}
\arguments{
\item{eset}{ExpressionSet object of the single cell data}

\item{cluster}{the variable name for "cell types" or "cluster name"}

\item{sample}{the variable name for subjects}

\item{select.ct}{the vector of cell types of interest, for example, c("alpha","beta")}

\item{Palette}{the color palette to be used in the demographic figure}
}
\value{
A figure showing the number of cells and percentage of cell types for single cells after clustering
}
\description{
Create a demographic figure for the input single cell dataset for each of the subject
}
