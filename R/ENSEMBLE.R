##################################
#####   ENSEMBLE PROCEDURE    #####
##################################
library(L1pack)

#Default methods for weight calculation includes minimize: RMSD/mAD of Y, Pearson

##################################################
#' ENSEMBLE function
#' @description ENSEMBLE function for deconvolution results output from SCDC_prop.
#' @name SCDC_ENSEMBLE
#' @param bulk.eset ExpressionSet object for bulk samples
#' @param sc.eset.list list of ExpressionSet objects for single cell reference datasets. Note that the variable names of multiple ExpressionSet objects should be consistent. Not required if prop.input is specified.
#' @param ct.sub vector. a subset of cell types that are selected to construct basis matrix. Not required if prop.input is specified.
#' @param ct.varname character string specifying the variable name for 'cell types'. Not required if prop.input is specified.
#' @param sample character string specifying the variable name for subject/samples. Not required if prop.input is specified.
#' @param grid.search logical. whether to allow grid search method to derive the ENSEMBLE weights.
#' @param search.length a number between 0 to 0.5. if using "Grid search", the step length used. Smaller search.length derives more accurate optimization results.
#' @param iter.max the maximum number of iteration in WNNLS
#' @param nu a small constant to facilitate the calculation of variance
#' @param epsilon a small constant number used for convergence criteria
#' @param truep true cell-type proportions for bulk samples if known
#' @param prop.input A list of SCDC_prop objects. Default is NULL. Users can use the SCDC function CreateSCDCpropObj to create the SCDC_prop objects first, and then put all objects into a list as this input. This allows users to ENSEMBLE results calculated from other deconvolution methods.
#' @param ct.cell.size default is NULL, which means the "library size" is calculated based on the data. Users can specify a vector of cell size factors corresponding to the ct.sub according to prior knowledge. The vector should be named: names(ct.cell.size input) should not be NULL.
#' @import L1pack
#' @import nnls
#' @import Biobase
#' @export
SCDC_ENSEMBLE <- function(bulk.eset, sc.eset.list = NULL, ct.varname, sample,
                          ct.sub, grid.search = F, search.length = 0.05,
                          iter.max = 2000, nu = 1e-04, epsilon = 0.001,
                          truep = NULL, weight.basis =T,
                          prop.input = NULL, Transform_bisque = F, ct.cell.size = NULL,
                          ...){
  # STEP 1: CALCULATE PROPORTION USING REF SEPARATELY.
  if (!is.null(prop.input)){
    message("Using user-input estimated proportion list ...")
    prop.list <- prop.input
  } else {
    prop.list <- lapply(sc.eset.list, function(zz){
      if (length(unique(zz@phenoData@data[,sample])) > 1){
        SCDC_prop(bulk.eset = bulk.eset, sc.eset = zz, ct.varname = ct.varname, sample = sample, truep = truep,
                  ct.sub = ct.sub, iter.max = iter.max, nu = nu, epsilon = epsilon, weight.basis = weight.basis,
                  Transform_bisque = Transform_bisque, ct.cell.size = ct.cell.size)
      } else {
        SCDC_prop_ONE(bulk.eset = bulk.eset, sc.eset = zz, ct.varname = ct.varname, sample = sample, truep = truep,
                      ct.sub = ct.sub, iter.max = iter.max, nu = nu, epsilon = epsilon, weight.basis = weight.basis,
                      ct.cell.size = ct.cell.size)
      }
    })
  }

  row.list <- sapply(1:length(prop.list), function(x){
    rownames(prop.list[[x]]$yhat)
  })
  gene.prop <- Reduce("intersect", row.list)
  gene.prop2 <- intersect(gene.prop, rownames(bulk.eset))
  subj.order <- colnames(bulk.eset)
  ycpm <- getCPM0(exprs(bulk.eset)[gene.prop2,subj.order])
  g.filter <- rowSums(ycpm) < quantile(rowSums(ycpm), 0.95) & rowSums(ycpm) > quantile(rowSums(ycpm), 0.15)
  gene.use <- gene.prop2[g.filter]
  length(gene.use)
  yv <- c(getCPM0(exprs(bulk.eset)[gene.use,subj.order]))*1e5 #vectorize y to make computing faster. scale to 100,000

  y.list <- do.call(cbind, lapply(prop.list, function(x){
    c(getCPM0(x$yhat[gene.use,subj.order]))*1e5
  }))


  sse <- function(x,y){
    sum((x-y)^2, na.rm = T)
  }
  sae <- function(x,y){
    sum(abs(x-y), na.rm = T)
  }
  rmsd <- function(x,y){
    sqrt(mean((x-y)^2, na.rm = T))
  }

  # STEP2: WEIGHT CALCULATION
  # -------------------------------
  # sum of squared errors
  message("Searching ENSEMBLE weight by Sum of Squared Errors or Sum of Abs Errors ......")
  sses <- apply(y.list, 2, function(x){
    sse(yv, x)
  })
  sse.wt <- 1/sses / sum(1/sses)
  # sum of absolute errors
  saes <- apply(y.list, 2, function(x){
    sae(yv, x)
  })
  sae.wt <- 1/saes / sum(1/saes)

  # RMSD
  rmsds <- apply(y.list, 2, function(x){
    rmsd(yv, x)
  })
  rmsd.wt <- 1/rmsds / sum(1/rmsds)


  message("Searching ENSEMBLE weight by LAD -- Minimizing mAD of Y measurement")
  w_lad <-NA
  dt <- data.frame(y = yv, y.list)
  fitlad <- L1pack::lad(y~.-1, data = dt, method = "BR")
  w_lad <- fitlad$coefficients
  w_lad[w_lad <0] <- 0
  w_lad <- w_lad/sum(w_lad)

  w_nnls <- NA
  message("Searching ENSEMBLE weight by NNLS -- Minimizing MSE of Y measurement")
  fitnnls <- nnls::nnls(A = as.matrix(y.list), b = yv)
  w_nnls <- fitnnls$x
  w_nnls <- w_nnls / sum(w_nnls)

  #-------------------------

  # combo of proportions according to selected weights
  combo <- lapply(prop.list, function(x){
    x$prop.est.mvw
  })

  # -----------------------------
  # grid search
  w_p_pearson <- NA; w_mad <- NA; w_rmsd <- NA; w_spearman <- NA; w_y_pearson <- NA
  gridres <- NULL
  if (grid.search){
    message("Grid search for ENSEMBLE weight ...")
    # GRID SEARCH MATRIX:
    gridmat <- getSearchGrid(lengthby = search.length, nparam = length(prop.list))
    ptm <- proc.time()
    search.prop <- NULL
    if (!is.null(truep)){
      message("Searching according to proportion--Pearson Correlation...")
      search.prop <- t(
        apply(gridmat, 1, function(x){
          temp <- wt_prop(x, proplist = combo)
          w_eval <- SCDC_peval(ptrue = as.matrix(truep), pest = temp,
                               pest.names = c("SCDC"),
                               select.ct = ct.sub)
          w_eval$evals.table
        })
      )
      colnames(search.prop) <- c("RMSD","mAD","R")
      w_p_pearson <- gridmat[which.max(search.prop[,3]),]
    }

    message("Searching according to bulk expression measurement...")
    search.y <- t(
      apply(gridmat, 1, function(x){
        temp <- wt_y(x, y.list = y.list)
        # mse, mad, pearson
        c(sqrt(mean((yv-temp)^2)), mean(abs(yv - temp)), cor(yv, temp, method = "pearson"), cor(yv, temp, method = "spearman"))
      })
    )
    colnames(search.y) <- c("RMSD_Y","mAD_Y","Pearson_Y","Spearman_Y")
    ptm2 <- proc.time() - ptm
    message("Grid search used", ptm2[3]," seconds.")
    # summarize weight results:
    w_rmsd <- gridmat[which.min(search.y[,1]),]
    w_mad <- gridmat[which.min(search.y[,2]),]
    w_y_pearson <- gridmat[which.max(search.y[,3]),]
    w_spearman <- gridmat[which.max(search.y[,4]),]
  }

  # ----------------------------------------------
  # summarize all weights and performances
  if (!is.null(truep)){
    weight.mat <- rbind(
      sse.wt,
      sae.wt,
      rmsd.wt,
      w_lad,
      w_nnls,
      w_p_pearson,
      w_mad,
      w_rmsd,
      w_spearman)
    rownames(weight.mat) <- c("inverse SSE", "inverse SAE","inverse RMSD","LAD","NNLS",
                              "Pearson_prop", "mAD_Y","RMSD_Y","Spearman_Y")
    eval.prop <- t(
      apply(weight.mat, 1, function(x){
        temp <- wt_prop(x, proplist = combo)
        w_eval <- SCDC_peval(ptrue = as.matrix(truep), pest = temp,
                             pest.names = c("SCDC"),
                             select.ct = ct.sub)
        w_eval$evals.table
      })
    )
    colnames(eval.prop) <- c("RMSD_prop","mAD_prop","R_prop")
    eval.y <- t(
      apply(weight.mat, 1, function(x){
        temp <- wt_y(x, y.list = y.list)
        # mse, mad, pearson
        c(sqrt(mean((yv-temp)^2)), mean(abs(yv - temp)), cor(yv, temp, method = "pearson"), cor(yv, temp, method = "spearman"))
      })
    )
    colnames(eval.y) <- c("RMSD_Y", "mAD_Y","Pearson_Y","Spearman_Y")
    if (grid.search){
      gridres <- cbind(search.prop, search.y, gridmat)
    }
    out <- cbind(round(weight.mat,2), eval.y, eval.prop)
  } else {
    weight.mat <- rbind(
      sse.wt,
      sae.wt,
      rmsd.wt,
      w_lad,
      w_nnls,
      w_mad,
      w_rmsd,
      w_spearman)
    rownames(weight.mat) <- c("inverse SSE", "inverse SAE","inverse RMSD","LAD","NNLS",
                              "mAD_Y","RMSD_Y","Spearman_Y")
    eval.y <- t(
      apply(weight.mat, 1, function(x){
        temp <- wt_y(x, y.list = y.list)
        # mse, mad, pearson
        c(sqrt(mean((yv-temp)^2)), mean(abs(yv - temp)), cor(yv, temp, method = "pearson"), cor(yv, temp, method = "spearman"))
      })
    )
    colnames(eval.y) <- c("RMSD_Y", "mAD_Y","Pearson_Y","Spearman_Y")
    if (grid.search){
      gridres <- cbind(search.y, gridmat)
    }
    out <- cbind(round(weight.mat,2), eval.y)
  }
  return(list(w_table = out, prop.list = prop.list, prop.only = combo, gridres = gridres))
}
