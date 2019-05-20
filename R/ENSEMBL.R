##################################
#####   ENSEMBL PROCEDURE    #####
##################################

##################################################
#' ENSEMBL function
#' @description ENSEMBL function for deconvolution results output from SCDC_prop
#' @name SCDC_ENSEMBLE
#' @param bulk.eset ExpressionSet object for bulk samples
#' @param sc.eset.list list of ExpressionSet objects for single cell reference datasets
#' @param ct.sub a subset of cell types that are selected to construct basis matrix
#' @param ct.varname variable name for 'cell types'
#' @param sample variable name for subject/samples
#' @param search.method method to derive the ENSEMBL weights, including "Grid search", "LAD", "NNLS".
#' @param search.length if using "Grid search", the step length used.
#' @param iter.max the maximum number of iteration in WNNLS
#' @param nu a small constant to facilitate the calculation of variance
#' @param epsilon a small constant number used for convergence criteria
#' @param truep true cell-type proportions for bulk samples if known
#' @export
SCDC_ENSEMBLE <- function(bulk.eset, sc.eset.list, ct.varname, sample,
                          ct.sub, search.method = c("Grid search", "LAD", "NNLS"), search.length = 0.05,
                          iter.max = 2000, nu = 1e-04, epsilon = 0.001, truep = NULL,
                          ...){
  prop.list <- lapply(sc.eset.list, function(zz){
    if (length(unique(zz@phenoData@data[,sample])) > 1){
      SCDC_prop(bulk.eset = bulk.eset, sc.eset = zz, ct.varname = ct.varname, sample = sample, truep = truep,
                ct.sub = ct.sub, iter.max = iter.max, nu = nu, epsilon = epsilon)
    } else {
      SCDC_propONE(bulk.eset = bulk.eset, sc.eset = zz, ct.varname = ct.varname, sample = sample, truep = truep,
                    ct.sub = ct.sub, iter.max = iter.max, nu = nu, epsilon = epsilon)
    }
  })

  row.list <- sapply(1:length(prop.list), function(x){
    rownames(prop.list[[x]]$yhat)
  })
  gene.use <- Reduce("intersect", row.list)
  yobs <- exprs(bulk.eset)[gene.use,]
  w_pearson <- NA; w_mad <- NA; w_rmsd <- NA; w_spearman <- NA
  if ("Grid search" %in% search.method){
    message("Searching ENSEMBLE weight by Grid search:")
    # GRID SEARCH MATRIX:
    gridmat <- getSearchGrid(lengthby = search.length, nparam = length(prop.list))
    testp <- NA
    ptm <- proc.time()
    if (!is.null(truep)){
      message("Searching according to proportion--Pearson Correlation...")
      testp <- sapply(1:dim(gridmat)[1], function(gg){
        temp <- matrix(0, ncol = length(ct.sub), nrow = dim(prop.list[[1]]$prop.est.mvw)[1])
        colnames(temp) <- ct.sub
        wtemp <- NULL
        for (z in 1:length(prop.list)){
          wtemp <- gridmat[gg,z]*prop.list[[z]]$prop.est.mvw
          mcol <- match(colnames(wtemp), colnames(temp))
          temp[,mcol] <- wtemp + temp[,mcol]
        }
        rownames(temp) <- rownames(truep)
        w_eval <- SCDC_peval(ptrue = as.matrix(truep), pest = temp,
                           pest.names = c("SCDC"),
                           select.ct = ct.sub)
        w_eval$evals.table[3]
      })
    }

    message("Searching according to bulk expression measurement...")
    testy <- sapply(1:dim(gridmat)[1], function(gg){
      temp <- 0
      for (z in 1:length(prop.list)){
        ytemp <- gridmat[gg,z]*prop.list[[z]]$yhat[gene.use,]
        temp <- ytemp + temp
      }
      w_eval <- SCDC_yeval(y = yobs, yest = temp, yest.names = "SCDC")
      w_eval$yevals.table
    })
    ptm2 <- proc.time() - ptm
    message("Grid search used", ptm2[3]," seconds.")
    # summarize results:
    w_spearman <- gridmat[which.max(t(testy)[,1]),]
    w_rmsd <- gridmat[which.min(t(testy)[,2]),]
    w_mad<- gridmat[which.min(t(testy)[,3]),]
    w_pearson <- gridmat[which.max(testp),]
  }

  xlist <- sapply(1:length(prop.list), function(z){
    c(getCPM0(prop.list[[z]]$yhat[gene.use,]))
  })
  colnames(xlist) <- names(prop.list)

  w_lad <-NA
  if ("LAD" %in% search.method){
    message("Searching ENSEMBLE weight by LAD -- Minimizing mAD of Y measurement")
    dt <- data.frame(y = c(getCPM0(yobs)),xlist)
    fitlad <- lad(y~.-1, data = dt, method = c("BR", "EM"))
    w_lad <- fitlad$coefficients
    w_lad[w_lad <0] <- 0
    w_lad <- w_lad/sum(w_lad)
  }

  w_nnls <- NA
  if ("NNLS" %in% search.method){
    message("Searching ENSEMBLE weight by NNLS -- Minimizing MSE of Y measurement")
    fitnnls <- nnls(A = as.matrix(xlist), b = c(getCPM0(yobs)))
    w_nnls <- fitnnls$x
    w_nnls <- w_nnls / sum(w_nnls)
  }

  wwfinal <- function(ww){
    temp <- matrix(0, ncol = length(ct.sub), nrow = dim(prop.list[[1]]$prop.est.mvw)[1])
    colnames(temp) <- ct.sub
    temp.y <- 0
    # ww=w_lad
    for (z in 1:length(prop.list)){
      wtemp <- prop.list[[z]]$prop.est.mvw * as.numeric(ww[z])
      mcol <- match(colnames(wtemp), colnames(temp))
      temp[,mcol] <- wtemp + temp[,mcol]
      ytemp <- as.numeric(ww[z])*prop.list[[z]]$yhat[gene.use,]
      temp.y <- ytemp + temp.y
    }
    rownames(temp) <- rownames(truep)
    w_eval <- list(evals.table = NULL)
    if (!is.null(truep)){
      w_eval <- SCDC_peval(ptrue = as.matrix(truep), pest = temp,
                         pest.names = c("SCDC"),
                         select.ct = ct.sub)
    }
    y_eval <- SCDC_yeval(y = yobs, yest = temp.y, yest.names = "SCDC")
    c(w_eval$evals.table,y_eval$yevals.table)
  }

  res_table <- rbind(wwfinal(w_pearson),
                     wwfinal(w_mad),
                     wwfinal(w_rmsd),
                     wwfinal(w_spearman),
                     wwfinal(w_lad),
                     wwfinal(w_nnls))
  w_vec <- rbind(w_pearson, w_mad, w_rmsd, w_spearman, w_lad, w_nnls)
  w_table <- cbind(w_vec[complete.cases(w_vec),],
                   res_table[complete.cases(res_table),])
  if ("Grid search" %in% search.method){
    if (!is.null(truep)){
      colnames(w_table)[(length(prop.list)+1):(length(prop.list)+6)] <- c("RMSD","mAD","Pearson","spearman_Y","RMSD_Y","mAD_Y")
      rownames(w_table) <- c('max.Pearson_prop','min.mAD_Y','min.RMSD_Y','max.Spearman_corr_Y',
                             'min.mAD_Y_LAD', 'min.RMSD_Y_NNLS')[complete.cases(res_table)]
      gridres <- cbind(t(testy), testp, gridmat)
      colnames(gridres) <- c("spearman_Y","RMSD_Y","mAD_Y","Pearson",colnames(gridmat))
    } else {
      w_table <- w_table[complete.cases(w_table),]
      colnames(w_table)[(length(prop.list)+1):(length(prop.list)+3)] <- c("spearman_Y","RMSD_Y","mAD_Y")
      rownames(w_table) <- c('min.mAD_Y','min.RMSD_Y','max.Spearman_corr_Y',
                             'min.mAD_Y_LAD', 'min.RMSD_Y_NNLS')
      gridres <- cbind(t(testy), gridmat)
      colnames(gridres) <- c("spearman_Y","RMSD_Y","mAD_Y",colnames(gridmat))
    }
  } else {
    if (!is.null(truep)){
      colnames(w_table)[(length(prop.list)+1):(length(prop.list)+6)] <- c("RMSD","mAD","Pearson","spearman_Y","RMSD_Y","mAD_Y")
      rownames(w_table) <- c('min.mAD_Y_LAD', 'min.RMSD_Y_NNLS')[complete.cases(res_table)]
      gridres <- NULL
    } else {
      w_table <- w_table[complete.cases(w_table),]
      colnames(w_table)[(length(prop.list)+1):(length(prop.list)+3)] <- c("spearman_Y","RMSD_Y","mAD_Y")
      rownames(w_table) <- c('min.mAD_Y_LAD', 'min.RMSD_Y_NNLS')
      gridres <- NULL
    }
  }
  return(list(w_table = w_table, prop.list = prop.list, gridres = gridres))
}


############################################################
#' ENSEMBL function for manually input deconvolution results
#' @description ENSEMBL function for manually input deconvolution results
#' @name SCDC_ENSEMBLE_subcl
#' @param bulk.eset ExpressionSet object for bulk samples
#' @param prop.list list of deconvolution results
#' @param ct.sub a subset of cell types that are selected to construct basis matrix
#' @param ct.varname variable name for 'cell types'
#' @param sample variable name for subject/samples
#' @param search.method method to derive the ENSEMBL weights, including "Grid search", "LAD", "NNLS".
#' @param search.length if using "Grid search", the step length used.
#' @param iter.max the maximum number of iteration in WNNLS
#' @param nu a small constant to facilitate the calculation of variance
#' @param epsilon a small constant number used for convergence criteria
#' @param truep true cell-type proportions for bulk samples if known
#' @export
SCDC_ENSEMBLE_subcl <- function(bulk.eset, prop.list, truep = NULL,
                                ct.sub, search.length = 0.05, search.method = c("Grid search", "LAD", "NNLS"),
                                iter.max = 2000, nu = 1e-04, epsilon = 0.001,
                                ...){

  row.list <- sapply(1:length(prop.list), function(x){
    rownames(prop.list[[x]]$yhat)
  })
  gene.use0 <- Reduce("intersect", row.list)
  gene.use <- intersect(gene.use0, rownames(bulk.eset))
  yobs <- exprs(bulk.eset)[gene.use,]
  w_pearson <- NA; w_mad <- NA; w_rmsd <- NA; w_spearman <- NA
  if ("Grid search" %in% search.method){
    message("Searching ENSEMBLE weight by Grid search:")
    # GRID SEARCH MATRIX:
    gridmat <- getSearchGrid(lengthby = search.length, nparam = length(prop.list))
    testp <- NA
    ptm <- proc.time()
    if (!is.null(truep)){
      message("Searching according to proportion--Pearson Correlation...")
      testp <- sapply(1:dim(gridmat)[1], function(gg){
        temp <- matrix(0, ncol = length(ct.sub), nrow = dim(prop.list[[1]]$prop.est)[1])
        colnames(temp) <- ct.sub
        wtemp <- NULL
        for (z in 1:length(prop.list)){
          wtemp <- gridmat[gg,z]*prop.list[[z]]$prop.est
          mcol <- match(colnames(wtemp), colnames(temp))
          temp[,mcol] <- wtemp + temp[,mcol]
        }
        rownames(temp) <- rownames(truep)
        w_eval <- SCDC_peval(ptrue = as.matrix(truep), pest = temp,
                           pest.names = c("SCDC"),
                           select.ct = ct.sub)
        w_eval$evals.table[3]
      })
    }

    message("Searching according to bulk expression measurement...")
    testy <- sapply(1:dim(gridmat)[1], function(gg){
      temp <- 0
      for (z in 1:length(prop.list)){
        # z=1
        # gg=1
        ytemp <- gridmat[gg,z]*prop.list[[z]]$yhat[gene.use,]
        temp <- ytemp + temp
      }
      head(temp)
      head(yobs)
      w_eval <- SCDC_yeval(y = yobs, yest = temp, yest.names = "SCDC")
      w_eval$yevals.table
    })
    ptm2 <- proc.time() - ptm
    message("Grid search used", ptm2[3]," seconds.")
    # summarize results:
    w_spearman <- gridmat[which.max(t(testy)[,1]),]
    w_rmsd <- gridmat[which.min(t(testy)[,2]),]
    w_mad <- gridmat[which.min(t(testy)[,3]),]
    w_pearson <- gridmat[which.max(testp),]
  }

  xlist <- sapply(1:length(prop.list), function(z){
    c(getCPM0(prop.list[[z]]$yhat[gene.use,]))
  })
  colnames(xlist) <- names(prop.list)

  w_lad <-NA
  if ("LAD" %in% search.method){
    message("Searching ENSEMBLE weight by LAD -- Minimizing mAD of Y measurement")
    dt <- data.frame(y = c(getCPM0(yobs)),xlist)
    fitlad <- lad(y~.-1, data = dt, method = c("BR", "EM"))
    w_lad <- fitlad$coefficients
    w_lad[w_lad <0] <- 0
    w_lad <- w_lad/sum(w_lad)
  }

  w_nnls <- NA
  if ("NNLS" %in% search.method){
    message("Searching ENSEMBLE weight by NNLS -- Minimizing MSE of Y measurement")
    fitnnls <- nnls(A = as.matrix(xlist), b = c(getCPM0(yobs)))
    w_nnls <- fitnnls$x
    w_nnls <- w_nnls / sum(w_nnls)
  }

  wwfinal <- function(ww){
    temp <- matrix(0, ncol = length(ct.sub), nrow = dim(prop.list[[1]]$prop.est)[1])
    colnames(temp) <- ct.sub
    temp.y <- 0
    # ww=w_lad
    for (z in 1:length(prop.list)){
      wtemp <- prop.list[[z]]$prop.est * as.numeric(ww[z])
      mcol <- match(colnames(wtemp), colnames(temp))
      temp[,mcol] <- wtemp + temp[,mcol]
      ytemp <- as.numeric(ww[z])*prop.list[[z]]$yhat[gene.use,]
      temp.y <- ytemp + temp.y
    }
    rownames(temp) <- rownames(truep)
    w_eval <- list(evals.table = NULL)
    if (!is.null(truep)){
      w_eval <- SCDC_peval(ptrue = as.matrix(truep), pest = temp,
                         pest.names = c("SCDC"),
                         select.ct = ct.sub)
    }
    y_eval <- SCDC_yeval(y = yobs, yest = temp.y, yest.names = "SCDC")
    c(w_eval$evals.table,y_eval$yevals.table)
  }

  res_table <- rbind(wwfinal(w_pearson),
                     wwfinal(w_mad),
                     wwfinal(w_rmsd),
                     wwfinal(w_spearman),
                     wwfinal(w_lad),
                     wwfinal(w_nnls))
  w_vec <- rbind(w_pearson, w_mad, w_rmsd, w_spearman, w_lad, w_nnls)
  w_table <- cbind(w_vec[complete.cases(w_vec),],
                   res_table[complete.cases(res_table),])
  if (!is.null(truep)){
    colnames(w_table)[(dim(gridmat)[2]+1):(dim(gridmat)[2]+6)] <- c("RMSD","mAD","Pearson","spearman_Y","RMSD_Y","mAD_Y")
    rownames(w_table) <- c('max.Pearson_prop','min.mAD_Y','min.RMSD_Y','max.Spearman_corr_Y',
                           'min.mAD_Y_LAD', 'min.RMSD_Y_NNLS')[complete.cases(res_table)]
    gridres <- cbind(t(testy), testp, gridmat)
    colnames(gridres) <- c("spearman_Y","RMSD_Y","mAD_Y","Pearson",colnames(gridmat))
  } else {
    w_table <- w_table[complete.cases(w_table),]
    colnames(w_table)[(dim(gridmat)[2]+1):(dim(gridmat)[2]+3)] <- c("spearman_Y","RMSD_Y","mAD_Y")
    rownames(w_table) <- c('min.mAD_Y','min.RMSD_Y','max.Spearman_corr_Y',
                           'min.mAD_Y_LAD', 'min.RMSD_Y_NNLS')
    gridres <- cbind(t(testy), gridmat)
    colnames(gridres) <- c("spearman_Y","RMSD_Y","mAD_Y",colnames(gridmat))
  }
  return(list(w_table = w_table, prop.list = prop.list, gridres = gridres))
}
