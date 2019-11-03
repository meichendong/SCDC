#######################################
####### DECONVOLUTION FUNCTIONS #######
#######################################

############################################
#' Basis Matrix
#' @description Basis matrix construction
#' @name SCDC_basis
#' @param x ExpressionSet object for single cells
#' @param ct.sub a subset of cell types that are selected to construct basis matrix
#' @param ct.varname variable name for 'cell types'
#' @param sample variable name for subject/samples
#' @return a list of basis matrix, sum of cell-type-specific library size, sample variance matrix, basis matrix by mvw, mvw matrix.
#' @export
SCDC_basis <- function(x, ct.sub = NULL, ct.varname, sample){
  # select only the subset of cell types of interest
  if (is.null(ct.sub)){
    ct.sub <- unique(x@phenoData@data[,ct.varname])
  }
  ct.sub <- ct.sub[!is.na(ct.sub)]
  x.sub <- x[,x@phenoData@data[,ct.varname] %in% ct.sub]
  # qc: remove non-zero genes
  x.sub <- x.sub[rowSums(exprs(x.sub)) > 0,]
  # calculate sample mean & sample variance matrix: genes by cell types
  countmat <- exprs(x.sub)
  ct.id <- droplevels(as.factor(x.sub@phenoData@data[,ct.varname]))
  sample.id <- as.character(x.sub@phenoData@data[,sample])
  ct_sample.id <- paste(ct.id,sample.id, sep = '%')

  mean.mat <- sapply(unique(ct_sample.id), function(id){
    y = as.matrix(countmat[, ct_sample.id %in% id])
    apply(y,1,sum, na.rm = TRUE)/sum(y)
  })
  mean.id <- do.call('rbind',strsplit(unique(ct_sample.id), split = '%'))

  sigma <- sapply(unique(mean.id[,1]), function(id){
    y = mean.mat[,mean.id[,1] %in% id]
    apply(y,1,var, na.rm = TRUE)
  })

  sum.mat2 <- sapply(unique(sample.id), function(sid){
    sapply(unique(ct.id), function(id){
      y = as.matrix(countmat[, ct.id %in% id & sample.id %in% sid])
      sum(y)/ncol(y)
    })
  })
  rownames(sum.mat2) <- unique(ct.id)
  colnames(sum.mat2) <- unique(sample.id)
  sum.mat <- rowMeans(sum.mat2, na.rm = T)

  basis <- sapply(unique(mean.id[,1]), function(id){
    z <- sum.mat[mean.id[,1]]
    mean.mat.z <- t(t(mean.mat)*z)
    y = as.matrix(mean.mat.z[,mean.id[,1] %in% id])
    apply(y,1,mean, na.rm = TRUE)
  })

  # weighted basis matrix
  my.max <- function(x,...){
    y <- apply(x,1,max, na.rm = TRUE)
    y / median(y, na.rm = T)
  }

  # MATCH DONOR, CELLTYPE, GENES!!!!!!!!!!!!!!!!
  var.adj <- sapply(unique(sample.id), function(sid) {
    my.max(sapply(unique(ct.id), function(id) {
      y = countmat[, ct.id %in% id & sample.id %in% sid,
                   drop = FALSE]
      apply(y,1,var, na.rm=T)
    }), na.rm = T)
  })
  colnames(var.adj) <- unique(sample.id)

  q15 <- apply(var.adj,2,quantile, probs = 0.15, na.rm =T)
  q85 <- apply(var.adj,2,quantile, probs = 0.85, na.rm =T)

  var.adj.q <- t(apply(var.adj, 1,
                       function(y){y[y<q15] <- q15[y<q15]
                       y[y>q85] <- q85[y>q85]
                       return(y)})) + 1e-4

  message("Creating Basis Matrix adjusted for maximal variance weight")
  mean.mat.mvw <- sapply(unique(ct_sample.id), function(id){
    sid = unlist(strsplit(id,'%'))[2]
    y = as.matrix(countmat[, ct_sample.id %in% id])
    yy = sweep(y, 1, sqrt(var.adj.q[,sid]), '/')
    apply(yy,1,sum, na.rm = TRUE)/sum(yy)
  })

  basis.mvw <- sapply(unique(mean.id[,1]), function(id){
    z <- sum.mat[mean.id[,1]]
    mean.mat.z <- t(t(mean.mat.mvw)*z)
    y = as.matrix(mean.mat.z[,mean.id[,1] %in% id])
    apply(y,1,mean, na.rm = TRUE)
  })

  # reorder columns
  basis.mvw <- basis.mvw[,ct.sub]
  sigma <- sigma[, ct.sub]
  basis <- basis[, ct.sub]
  sum.mat <- sum.mat[ct.sub]

  return(list(basis = basis, sum.mat = sum.mat,
              sigma = sigma, basis.mvw = basis.mvw, var.adj.q = var.adj.q))
}

#############################################
#' Basis matrix for single cells from one subject
#' @description Basis matrix construction for single cells from one subject
#' @name SCDC_basis_ONE
#' @param x ExpressionSet object for single cells
#' @param ct.sub a subset of cell types that are selected to construct basis matrix
#' @param ct.varname variable name for 'cell types'
#' @param sample variable name for subject/samples
#' @return a list of basis matrix, sum of cell-type-specific library size, sample variance matrix, basis matrix by mvw, mvw matrix.
#' @export
SCDC_basis_ONE <- function(x , ct.sub = NULL, ct.varname, sample){
  # select only the subset of cell types of interest
  if (is.null(ct.sub)){
    ct.sub <- unique(x@phenoData@data[,ct.varname])[!is.na(unique(x@phenoData@data[,ct.varname]))]
  }
  ct.sub <- ct.sub[!is.na(ct.sub)]
  x.sub <- x[,x@phenoData@data[,ct.varname] %in% ct.sub]
  # qc: remove non-zero genes
  x.sub <- x.sub[rowSums(exprs(x.sub)) > 0,]
  # calculate sample mean & sample variance matrix: genes by cell types
  countmat <- exprs(x.sub)
  # ct.id <- droplevels(as.factor(x.sub@phenoData@data[,ct.varname]))
  ct.id <- x.sub@phenoData@data[,ct.varname]
  sample.id <- x.sub@phenoData@data[,sample]
  ct_sample.id <- paste(ct.id,sample.id, sep = '%')

  mean.mat <- sapply(unique(ct_sample.id), function(id){
    y = as.matrix(countmat[, ct_sample.id %in% id])
    apply(y,1,sum, na.rm = TRUE)/sum(y)
  })
  mean.id <- do.call('rbind',strsplit(unique(ct_sample.id), split = '%'))

  # by subj, then take avg????
  sum.mat2 <- sapply(unique(sample.id), function(sid){
    sapply(unique(ct.id), function(id){
      y = as.matrix(countmat[, ct.id %in% id & sample.id %in% sid])
      sum(y)/ncol(y)
    })
  })
  rownames(sum.mat2) <- unique(ct.id)
  colnames(sum.mat2) <- unique(sample.id)
  sum.mat <- rowMeans(sum.mat2, na.rm = T)

  basis <- sapply(unique(mean.id[,1]), function(id){
    z <- sum.mat[mean.id[,1]]
    mean.mat.z <- t(t(mean.mat)*z)
    # id = unique(mean.id[,1])[1]
    y = as.matrix(mean.mat.z[,mean.id[,1] %in% id])
    apply(y,1,mean, na.rm = TRUE)
  })

  # weighted basis matrix
  my.max <- function(x,...){
    y <- apply(x,1,max, na.rm = TRUE)
    y / median(y, na.rm = T)
  }

  # MATCH DONOR, CELLTYPE, GENES!!!!!!!!!!!!!!!!
  var.adj <- sapply(unique(sample.id), function(sid) {
    my.max(sapply(unique(ct.id), function(id) {
      y = countmat[, ct.id %in% id & sample.id %in% sid,
                   drop = FALSE]
      apply(y,1,var, na.rm=T)
    }), na.rm = T)
  })
  colnames(var.adj) <- unique(sample.id)

  q15 <- apply(var.adj,2,quantile, probs = 0.15, na.rm =T)
  q85 <- apply(var.adj,2,quantile, probs = 0.85, na.rm =T)

  var.adj.q <- as.matrix(apply(var.adj, 1,
                               function(y){y[y<q15] <- q15[y<q15]
                               y[y>q85] <- q85[y>q85]
                               return(y)}) + 1e-4)

  message("Creating Basis Matrix adjusted for maximal variance weight")
  mean.mat.mvw <- sapply(unique(ct_sample.id), function(id){
    y = as.matrix(countmat[, ct_sample.id %in% id])
    yy = sweep(y, 1, sqrt(var.adj.q), '/')
    apply(yy,1,sum, na.rm = TRUE)/sum(yy)
  })

  basis.mvw <- sapply(unique(mean.id[,1]), function(id){
    z <- sum.mat[mean.id[,1]]
    mean.mat.z <- t(t(mean.mat.mvw)*z)
    y = as.matrix(mean.mat.z[,mean.id[,1] %in% id])
    apply(y,1,mean, na.rm = TRUE)
  })

  # reorder columns
  basis.mvw <- basis.mvw[,ct.sub]
  sigma <- NULL # in the one subject case, no variance is calculated.
  basis <- basis[, ct.sub]
  sum.mat <- sum.mat[ct.sub]

  return(list(basis = basis, sum.mat = sum.mat,
              sigma = sigma, basis.mvw = basis.mvw, var.adj.q = var.adj.q))
}

#################################
#' Clustering QC
#' @description Single cells Clustering QC
#' @name SCDC_qc
#' @import pheatmap
#' @param sc.eset ExpressionSet object for single cells
#' @param ct.varname variable name for 'cell type'
#' @param sample variable name for subject/sample
#' @param scsetname the name for the single cell dataset
#' @param ct.sub a subset of cell types that are selected to construct basis matrix
#' @param iter.max the maximum number of iteration in WNNLS
#' @param nu a small constant to facilitate the calculation of variance
#' @param epsilon a small constant number used for convergence criteria
#' @param arow annotation of rows for pheatmap
#' @param qcthreshold the probability threshold used to filter out questionable cells
#' @return a list including: 1) a probability matrix for each single cell input; 2) a clustering QCed ExpressionSet object; 3) a heatmap of QC result.
#' @export
SCDC_qc <- function (sc.eset, ct.varname, sample, scsetname = "Single Cell",
                   ct.sub, iter.max = 1000, nu = 1e-04, epsilon = 0.01, arow =NULL,
                   qcthreshold = 0.7, ...) {
  sc.basis = SCDC_basis(x = sc.eset, ct.sub = ct.sub, ct.varname = ct.varname, sample = sample)
  M.S <- sc.basis$sum.mat[ct.sub]
  xsc <- getCPM0(exprs(sc.eset)[rownames(sc.basis$basis.mvw),])
  N.sc <- ncol(xsc)

  m.basis <- sc.basis$basis.mvw[, ct.sub]
  sigma <- sc.basis$sigma[, ct.sub]
  valid.ct <- (colSums(is.na(sigma)) == 0) & (colSums(is.na(m.basis)) ==
                                                0) & (!is.na(M.S))
  if (sum(valid.ct) <= 1) {
    stop("Not enough valid cell type!")
  }
  message(paste("Used", sum(valid.ct), "cell types in deconvolution..."))
  m.basis <- m.basis[, valid.ct]
  M.S <- M.S[valid.ct]
  sigma <- sigma[, valid.ct]

  prop.qc <- NULL
  for (i in 1:N.sc) {
    message("Begin iterative weighted estimation...")
    basis.temp <- m.basis
    xsc.temp <- xsc[, i]
    sigma.temp <- sigma
    ### weighting scheme:
    lm.qc <- nnls::nnls(A=basis.temp,b=xsc.temp)
    delta <- lm.qc$residuals
    wt.gene <- 1/(nu + delta^2 + colSums((lm.qc$x)^2*t(sigma.temp)))
    x.wt <- xsc.temp*sqrt(wt.gene)
    b.wt <- sweep(basis.temp,1,sqrt(wt.gene),"*")
    lm.wt <- nnls::nnls(A=b.wt, b=x.wt)
    prop.wt <- lm.wt$x/sum(lm.wt$x)
    delta <- lm.wt$residuals
    for (iter in 1:iter.max){
      wt.gene <- 1/(nu + delta^2 + colSums((lm.wt$x)^2*t(sigma.temp)))
      x.wt <- xsc.temp*sqrt(wt.gene)
      b.wt <- sweep(basis.temp,1,sqrt(wt.gene),"*")
      lm.wt <- nnls::nnls(A=b.wt, b=x.wt)
      delta.new <- lm.wt$residuals
      prop.wt.new <- lm.wt$x/sum(lm.wt$x)
      if (sum(abs(prop.wt - prop.wt.new) < epsilon )){
        prop.wt <- prop.wt.new
        delta <- delta.new
        message("Converged at iteration ", iter)
        break
      }
      prop.wt <- prop.wt.new
      delta <- delta.new
    }

    prop.qc <- rbind(prop.qc, prop.wt)
  }
  # name col and row
  colnames(prop.qc) <- colnames(m.basis)
  rownames(prop.qc) <- colnames(xsc)
  heat.anno <- pheatmap(prop.qc, annotation_row = arow,
                        annotation_names_row=FALSE, show_rownames = F,
                        annotation_names_col=FALSE, cutree_rows = length(ct.sub),
                        color = cbPalette[2:4],
                        cluster_rows = T, cluster_cols = F)

  prop.qc.keep <- rowSums(prop.qc > qcthreshold) ==1 # truncated values -> F or T
  sc.eset.qc <- sc.eset[,prop.qc.keep]
  return(list(prop.qc = prop.qc, sc.eset.qc = sc.eset.qc, heatfig = heat.anno))
}

#################################
#' Clustering QC for single cells from one subject
#' @description Clustering QC for single cells from one subject
#' @name SCDC_qc_ONE
#' @param sc.eset ExpressionSet object for single cells
#' @param ct.varname variable name for 'cell type'
#' @param sample variable name for subject/sample
#' @param scsetname the name for the single cell dataset
#' @param ct.sub a subset of cell types that are selected to construct basis matrix
#' @param iter.max the maximum number of iteration in WNNLS
#' @param nu a small constant to facilitate the calculation of variance
#' @param epsilon a small constant number used for convergence criteria
#' @param arow annotation of rows for pheatmap
#' @param qcthreshold the probability threshold used to filter out questionable cells
#' @return a list including: 1) a probability matrix for each single cell input; 2) a clustering QCed ExpressionSet object; 3) a heatmap of QC result.
#' @export
SCDC_qc_ONE <- function(sc.eset, ct.varname, sample, scsetname = "Single Cell",
                  ct.sub, iter.max = 1000, nu = 1e-04, epsilon = 0.01,
                    arow = NULL, weight.basis = F, qcthreshold = 0.7, ...){
  sc.basis <- SCDC_basis_ONE(x = sc.eset, ct.sub = ct.sub, ct.varname = ct.varname, sample = sample)
  if (weight.basis){
    basis.mvw <- sc.basis$basis.mvw[, ct.sub]
  } else {
    basis.mvw <- sc.basis$basis[, ct.sub]
  }

  xsc <- getCPM0(exprs(sc.eset))
  N.sc <- ncol(xsc)

  ALS.S <- sc.basis$sum.mat[ct.sub]
  valid.ct <- (colSums(is.na(basis.mvw)) == 0) & (!is.na(ALS.S))
  if (sum(valid.ct) <= 1) {
    stop("Not enough valid cell type!")
  }
  message(paste("Used", sum(valid.ct), "cell types in deconvolution..."))
  basis.mvw <- basis.mvw[, valid.ct]
  ALS.S <- ALS.S[valid.ct]
  prop.est.mvw <- NULL

  # prop estimation for each sc sample:
  for (i in 1:N.sc) {
    xsc.i <- xsc[, i]*100 # why times 100 if not normalize???
    gene.use <- intersect(rownames(basis.mvw), names(xsc.i))
    basis.mvw.temp <- basis.mvw[gene.use,]
    xsc.temp <- xsc.i[gene.use]
    message(paste(colnames(xsc)[i], "has common genes", sum(xsc[, i] != 0), "..."))
    # first NNLS:

    lm <- nnls::nnls(A=basis.mvw.temp,b=xsc.temp)
    delta <- lm$residuals
    wt.gene <- 1/(nu + delta^2)
    x.wt <- xsc.temp*sqrt(wt.gene)
    b.wt <- sweep(basis.mvw.temp,1,sqrt(wt.gene),"*")

    lm.wt <- nnls::nnls(A=b.wt, b=x.wt)
    prop.wt <- lm.wt$x/sum(lm.wt$x)
    delta <- lm.wt$residuals

    for (iter in 1:iter.max){
      wt.gene <- 1/(nu + delta^2)
      x.wt <- xsc.temp * sqrt(wt.gene)
      b.wt <- sweep(basis.mvw.temp,1,sqrt(wt.gene),"*")
      lm.wt <- nnls::nnls(A=b.wt, b=x.wt)
      delta.new <- lm.wt$residuals
      prop.wt.new <- lm.wt$x/sum(lm.wt$x)

      if (sum(abs(prop.wt.new - prop.wt)) < epsilon){
        prop.wt <- prop.wt.new
        delta <- delta.new
        message("WNNLS Converged at iteration ", iter)
        break
      }
      prop.wt <- prop.wt.new
      delta <- delta.new
    }

    prop.est.mvw <- rbind(prop.est.mvw, prop.wt)
  }
  colnames(prop.est.mvw) <- colnames(basis.mvw)
  rownames(prop.est.mvw) <- colnames(xsc)

  ### plot steps:
  heat.anno <- pheatmap(prop.est.mvw, annotation_row = arow,
                        annotation_names_row=FALSE, show_rownames = F,
                        annotation_names_col=FALSE, cutree_rows = length(ct.sub),
                        color = cbPalette[2:4],
                        cluster_rows = T, cluster_cols = F) #, main = scsetname

  prop.qc.keep <- rowSums(prop.est.mvw > qcthreshold) ==1 # truncated values -> F or T
  sc.eset.qc <- sc.eset[,prop.qc.keep]
  return(list(prop.qc = prop.est.mvw, sc.eset.qc = sc.eset.qc, heatfig = heat.anno))
}

######################################
#' Proportion estimation
#' @description Proportion estimation function for multi-subject case
#' @name SCDC_prop
#' @param bulk.eset ExpressionSet object for bulk samples
#' @param sc.eset ExpressionSet object for single cell samples
#' @param ct.varname variable name for 'cell types'
#' @param sample variable name for subject/samples
#' @param ct.sub a subset of cell types that are selected to construct basis matrix
#' @param iter.max the maximum number of iteration in WNNLS
#' @param nu a small constant to facilitate the calculation of variance
#' @param epsilon a small constant number used for convergence criteria
#' @param truep true cell-type proportions for bulk samples if known
#' @param Transform_bisque The bulk sample transformation from bisqueRNA. Aiming to reduce the systematic difference between single cells and bulk samples.
#' @return Estimated proportion, basis matrix, predicted gene expression levels for bulk samples
#' @export
SCDC_prop <- function (bulk.eset, sc.eset, ct.varname, sample, ct.sub, iter.max = 1000,
                       nu = 1e-04, epsilon = 0.01, truep = NULL, weight.basis = T,
                       Transform_bisque = F, ...)
{
  bulk.eset <- bulk.eset[rowSums(exprs(bulk.eset)) > 0, , drop = FALSE]
  ct.sub <- intersect(ct.sub, unique(sc.eset@phenoData@data[,
                                                            ct.varname]))
  sc.basis <- SCDC_basis(x = sc.eset, ct.sub = ct.sub, ct.varname = ct.varname,
                         sample = sample)
  commongenes <- intersect(rownames(sc.basis$basis.mvw), rownames(bulk.eset))
  if (length(commongenes) < 0.2 * min(dim(sc.eset)[1], dim(bulk.eset)[1])) {
    stop("Too few common genes!")
  }
  message(paste("Used", length(commongenes), "common genes..."))
  if (weight.basis) {
    basis.mvw <- sc.basis$basis.mvw[commongenes, ct.sub]
  }
  else {
    basis.mvw <- sc.basis$basis[commongenes, ct.sub]
  }
  # link to bisqueRNA, bulk transformation method. https://github.com/cozygene/bisque
  if (Transform_bisque) {
    GenerateSCReference <- function(sc.eset, ct.sub) {
      cell.labels <- base::factor(sc.eset[[ct.sub]])
      all.cell.types <- base::levels(cell.labels)
      aggr.fn <- function(ct.sub) {
        base::rowMeans(Biobase::exprs(sc.eset)[,cell.labels == ct.sub, drop=F])
      }
      template <- base::numeric(base::nrow(sc.eset))
      sc.ref <- base::vapply(all.cell.types, aggr.fn, template)
      return(sc.ref)
    }
    sc.ref <- GenerateSCReference(sc.eset, cell.types)[genes, , drop = F]
    ncount <- table(sc.eset@phenoData@data[, sample], sc.eset@phenoData@data[, ct.varname])
    true.prop <- ncount/rowSums(ncount, na.rm = T)
    sc.props <- round(true.prop[complete.cases(true.prop), ], 2)
    Y.train <- sc.ref %*% t(sc.props[, colnames(sc.ref)])
    dim(Y.train)
    X.pred <- exprs(bulk.eset)[commongenes, ]
    sample.names <- base::colnames(Biobase::exprs(bulk.eset))
    template <- base::numeric(base::length(sample.names))
    base::names(template) <- sample.names
    SemisupervisedTransformBulk <- function(gene, Y.train, X.pred) {
      Y.train.scaled <- base::scale(Y.train[gene, , drop = T])
      Y.center <- base::attr(Y.train.scaled, "scaled:center")
      Y.scale <- base::attr(Y.train.scaled, "scaled:scale")
      n <- base::length(Y.train.scaled)
      shrink.scale <- base::sqrt(base::sum((Y.train[gene, , drop = T] - Y.center)^2)/n + 1)
      X.pred.scaled <- base::scale(X.pred[gene, , drop = T])
      Y.pred <- base::matrix((X.pred.scaled * shrink.scale) +
                               Y.center, dimnames = base::list(base::colnames(X.pred),
                                                               gene))
      return(Y.pred)
    }
    Y.pred <- base::matrix(base::vapply(X = commongenes,
                                        FUN = SemisupervisedTransformBulk, FUN.VALUE = template,
                                        Y.train, X.pred, USE.NAMES = TRUE), nrow = base::length(sample.names))
    indices <- base::apply(Y.pred, MARGIN = 2, FUN = function(column) {
      base::anyNA(column)
    })
    if (base::any(indices)) {
      if (sum(!indices) == 0) {
        base::stop("Zero genes left for decomposition.")
      }
      Y.pred <- Y.pred[, !indices, drop = F]
      sc.ref <- sc.ref[!indices, , drop = F]
    }

    results <- base::as.matrix(base::apply(Y.pred, 1, function(b) {
      sol <- lsei::pnnls(sc.ref, b, sum = 1)
      return(sol$x)
    }))

    prop.est.mvw <- t(results)
    colnames(prop.est.mvw) <- colnames(sc.ref)
    rownames(prop.est.mvw) <- colnames(bulk.eset)
    yhat <- sc.ref %*% results
    colnames(yhat) <- colnames(bulk.eset)
    yobs <- exprs(bulk.eset)
    yeval <- SCDC_yeval(y = yobs, yest = yhat, yest.names = c("SCDC"))
    peval <- NULL
    if (!is.null(truep)) {
      peval <- SCDC_peval(ptrue = truep, pest = prop.est.mvw,
                          pest.names = c("SCDC"), select.ct = ct.sub)
    }
  } else {
    xbulk <- getCPM0(exprs(bulk.eset)[commongenes, ])
    sigma <- sc.basis$sigma[commongenes, ct.sub]
    ALS.S <- sc.basis$sum.mat[ct.sub]
    N.bulk <- ncol(bulk.eset)
    valid.ct <- (colSums(is.na(sigma)) == 0) & (colSums(is.na(basis.mvw)) ==
                                                  0) & (!is.na(ALS.S))
    if (sum(valid.ct) <= 1) {
      stop("Not enough valid cell type!")
    }
    message(paste("Used", sum(valid.ct), "cell types in deconvolution..."))
    basis.mvw <- basis.mvw[, valid.ct]
    ALS.S <- ALS.S[valid.ct]
    sigma <- sigma[, valid.ct]
    prop.est.mvw <- NULL
    yhat <- NULL
    yhatgene.temp <- rownames(basis.mvw)
    for (i in 1:N.bulk) {
      basis.mvw.temp <- basis.mvw
      xbulk.temp <- xbulk[, i]*100
      sigma.temp <- sigma
      message(paste(colnames(xbulk)[i], "has common genes",
                    sum(xbulk[, i] != 0), "..."))
      lm <- nnls::nnls(A = basis.mvw.temp, b = xbulk.temp)
      delta <- lm$residuals
      wt.gene <- 1/(nu + delta^2 + colSums((lm$x * ALS.S)^2 *
                                             t(sigma.temp)))
      x.wt <- xbulk.temp * sqrt(wt.gene)
      b.wt <- sweep(basis.mvw.temp, 1, sqrt(wt.gene), "*")
      lm.wt <- nnls::nnls(A = b.wt, b = x.wt)
      prop.wt <- lm.wt$x/sum(lm.wt$x)
      delta <- lm.wt$residuals
      for (iter in 1:iter.max) {
        wt.gene <- 1/(nu + delta^2 + colSums((lm.wt$x * ALS.S)^2 *
                                               t(sigma.temp)))
        x.wt <- xbulk.temp * sqrt(wt.gene)
        b.wt <- sweep(basis.mvw.temp, 1, sqrt(wt.gene), "*")
        lm.wt <- nnls::nnls(A = b.wt, b = x.wt)
        delta.new <- lm.wt$residuals
        prop.wt.new <- lm.wt$x/sum(lm.wt$x)
        if (sum(abs(prop.wt.new - prop.wt)) < epsilon) {
          prop.wt <- prop.wt.new
          delta <- delta.new
          R2 <- 1 - var(xbulk.temp - basis.mvw.temp %*%
                          as.matrix(lm.wt$x))/var(xbulk.temp)
          message("WNNLS Converged at iteration ",
                  iter)
          break
        }
        prop.wt <- prop.wt.new
        delta <- delta.new
      }
      R2 <- 1 - var(xbulk.temp - basis.mvw.temp %*% as.matrix(lm.wt$x))/var(xbulk.temp)
      prop.est.mvw <- rbind(prop.est.mvw, prop.wt)
      yhat.temp <- basis.mvw.temp %*% as.matrix(lm.wt$x)
      yhatgene.temp <- intersect(rownames(yhat.temp), yhatgene.temp)
      yhat <- cbind(yhat[yhatgene.temp, ], yhat.temp[yhatgene.temp,
                                                     ])
    }
    colnames(prop.est.mvw) <- colnames(basis.mvw)
    rownames(prop.est.mvw) <- colnames(xbulk)
    colnames(yhat) <- colnames(xbulk)
    yobs <- exprs(bulk.eset)
    yeval <- SCDC_yeval(y = yobs, yest = yhat, yest.names = c("SCDC"))
    peval <- NULL
    if (!is.null(truep)) {
      peval <- SCDC_peval(ptrue = truep, pest = prop.est.mvw,
                          pest.names = c("SCDC"), select.ct = ct.sub)
    }
  }

  return(list(prop.est.mvw = prop.est.mvw, basis.mvw = basis.mvw,
              yhat = yhat, yeval = yeval, peval = peval))
}

############################################
#' Proportion estimation function for one-subject case
#' @description Proportion estimation function for one-subject case
#' @name SCDC_prop_ONE
#' @param bulk.eset ExpressionSet object for bulk samples
#' @param sc.eset ExpressionSet object for single cell samples
#' @param ct.varname variable name for 'cell types'
#' @param sample variable name for subject/samples
#' @param ct.sub a subset of cell types that are selected to construct basis matrix
#' @param iter.max the maximum number of iteration in WNNLS
#' @param nu a small constant to facilitate the calculation of variance
#' @param epsilon a small constant number used for convergence criteria
#' @param truep true cell-type proportions for bulk samples if known
#' @return Estimated proportion, basis matrix, predicted gene expression levels for bulk samples
#' @export
SCDC_prop_ONE <- function(bulk.eset, sc.eset, ct.varname, sample, truep = NULL,
                          ct.sub, iter.max = 2000, nu = 1e-10, epsilon = 0.01,
                          weight.basis = T, ...){
  bulk.eset <- bulk.eset[rowSums(exprs(bulk.eset))>0, , drop = FALSE]
  sc.basis <- SCDC_basis_ONE(x = sc.eset, ct.sub = ct.sub, ct.varname = ct.varname, sample = sample)
  # match genes / cells first
  if (weight.basis){
    basis <- sc.basis$basis.mvw
  } else {
    basis <- sc.basis$basis
  }
  commongenes <- intersect(rownames(basis), rownames(bulk.eset))
  # stop when few common genes exist...
  if (length(commongenes) < 0.2 * min(dim(sc.eset)[1], dim(bulk.eset)[1])){
    stop('Too few common genes!')
  }
  message(paste("Used", length(commongenes), "common genes..."))
  basis.mvw <- basis[commongenes, ct.sub]
  xbulk <- getCPM0(exprs(bulk.eset)[commongenes,])

  ALS.S <- sc.basis$sum.mat[ct.sub]
  N.bulk <- ncol(bulk.eset)
  valid.ct <- (colSums(is.na(basis.mvw)) == 0) & (!is.na(ALS.S))
  if (sum(valid.ct) <= 1) {
    stop("Not enough valid cell type!")
  }
  message(paste("Used", sum(valid.ct), "cell types in deconvolution..."))
  basis.mvw <- basis.mvw[, valid.ct]
  ALS.S <- ALS.S[valid.ct]
  prop.est.mvw <- NULL

  yhat <- NULL
  yhatgene.temp <- rownames(basis.mvw)
  # prop estimation for each bulk sample:
  for (i in 1:N.bulk) {
    xbulk.temp <- xbulk[, i]# *100 # why times 100 if not normalize???
    message(paste(colnames(xbulk)[i], "has common genes", sum(xbulk[, i] != 0), "..."))
    # first NNLS:
    lm <- nnls::nnls(A=basis.mvw,b=xbulk.temp)
    delta <- lm$residuals
    wt.gene <- 1/(nu + delta^2)
    x.wt <- xbulk.temp*sqrt(wt.gene)
    b.wt <- sweep(basis.mvw,1,sqrt(wt.gene),"*")

    lm.wt <- nnls::nnls(A=b.wt, b=x.wt)
    prop.wt <- lm.wt$x/sum(lm.wt$x)
    delta <- lm.wt$residuals

    for (iter in 1:iter.max){
      wt.gene <- 1/(nu + delta^2)
      x.wt <- xbulk.temp * sqrt(wt.gene)
      b.wt <- sweep(basis.mvw,1,sqrt(wt.gene),"*")
      lm.wt <- nnls::nnls(A=b.wt, b=x.wt)
      delta.new <- lm.wt$residuals
      prop.wt.new <- lm.wt$x/sum(lm.wt$x)

      if (sum(abs(prop.wt.new - prop.wt)) < epsilon){
        prop.wt <- prop.wt.new
        delta <- delta.new
        message("WNNLS Converged at iteration ", iter)
        break
      }
      prop.wt <- prop.wt.new
      delta <- delta.new
    }

    prop.est.mvw <- rbind(prop.est.mvw, prop.wt)
    yhat.temp <- basis.mvw %*% as.matrix(lm.wt$x)
    yhatgene.temp <- intersect(rownames(yhat.temp), yhatgene.temp)
    yhat <- cbind(yhat[yhatgene.temp,], yhat.temp[yhatgene.temp,])
  }
  colnames(prop.est.mvw) <- colnames(basis.mvw)
  rownames(prop.est.mvw) <- colnames(bulk.eset)
  colnames(yhat) <- colnames(bulk.eset)

  ### evaluation steps:
  yobs <- exprs(bulk.eset)
  yeval <- SCDC_yeval(y=yobs, yest = yhat,
                      yest.names=c("SCDC"))
  peval <- NULL
  if (!is.null(truep)){
    peval <- SCDC_peval(ptrue= truep, pest = prop.est.mvw, pest.names = c('SCDC'),
                       select.ct = ct.sub)
  }
  return(list(prop.est.mvw = prop.est.mvw, basis.mvw = basis.mvw, yhat = yhat,
              yeval = yeval, peval = peval))
}

############################################
#' Tree-guided proportion estimation
#' @description Proportion estimation function for multi-subject case, and apply tree-guided deconvolution
#' @name SCDC_prop_subcl_marker
#' @param bulk.eset ExpressionSet object for bulk samples
#' @param sc.eset ExpressionSet object for single cell samples
#' @param ct.varname variable name for 'cell types'
#' @param fl.varname variable name for first-level 'meta-clusters'
#' @param sample variable name for subject/samples
#' @param ct.sub a subset of cell types that are selected to construct basis matrix
#' @param ct.fl.sub 'cell types' for first-level 'meta-clusters'
#' @param iter.max the maximum number of iteration in WNNLS
#' @param nu a small constant to facilitate the calculation of variance
#' @param epsilon a small constant number used for convergence criteria
#' @param weight.basis logical, use basis matrix adjusted by MVW, default is T.
#' @param select.marker logical, select marker genes to perform deconvolution in tree-guided steps. Default is T.
#' @param markers A set of marker gene that input manully to be used in deconvolution. If NULL, then
#' @param marker.varname variable name of cluster groups when selecting marker genes. If NULL, then use ct.varname.
#' @param allgenes.fl logical, use all genes in the first-level deconvolution
#' @param pseudocount.use a constant number used when selecting marker genes, default is 1.
#' @param LFC.lim a threshold of log fold change when selecting genes as input to perform Wilcoxon's test.
#' @param truep true cell-type proportions for bulk samples if known
#' @return Estimated proportion, basis matrix, predicted gene expression levels for bulk samples
#' @export
SCDC_prop_subcl_marker <- function(bulk.eset, sc.eset, ct.varname, fl.varname, sample,
                                   ct.sub = NULL, ct.fl.sub, iter.max = 3000, nu = 1e-04, epsilon = 0.001,
                                   weight.basis = T, truep = NULL,
                                   select.marker = T, markers = NULL, marker.varname = NULL, allgenes.fl = F,
                                   pseudocount.use = 1, LFC.lim = 0.5, ...) {
  if (is.null(ct.sub)){
    ct.sub <- unique(sc.eset@phenoData@data[,ct.varname])[!is.na(unique(sc.eset@phenoData@data[,ct.varname]))]
  }
  ct.sub <- ct.sub[!is.na(ct.sub)]
  ct.fl.sub <- ct.fl.sub[!is.na(ct.fl.sub)]
  bulk.eset <- bulk.eset[rowSums(exprs(bulk.eset))>0, , drop = FALSE]
  sc.eset <- sc.eset[,sc.eset@phenoData@data[,ct.varname] %in% ct.sub]
  sc.basis <- SCDC_basis(x = sc.eset, ct.sub = ct.sub, ct.varname = ct.varname, sample = sample)
  sc.fl.basis <- SCDC_basis(x = sc.eset, ct.sub = ct.fl.sub[!is.na(ct.fl.sub)],
                            ct.varname = fl.varname, sample = sample)
  if (select.marker){
    if (is.null(marker.varname)){
      marker.varname <- ct.varname
    }
    # wilcox test on two groups of cells for marker gene selection... (refer to seurat::FindMarkers)
    countmat <- exprs(sc.eset)
    ct.group <- sc.eset@phenoData@data[,marker.varname]
    markers.wilcox <- NULL
    # u=1
    for(u in 1:length(unique(ct.group))){
      ct.group.temp <- ct.group == unique(ct.group)[u]
      group.1 <- apply(X = countmat[,ct.group.temp],
                       MARGIN = 1, FUN = function(x) log(x = mean(x = expm1(x = x)) +
                                                           pseudocount.use))
      group.2 <- apply(X = countmat[,! ct.group.temp],
                       MARGIN = 1, FUN = function(x) log(x = mean(x = expm1(x = x)) +
                                                           pseudocount.use))
      genes.diff <- rownames(sc.eset)[(group.1 - group.2) > LFC.lim]
      count.use <- countmat[rownames(sc.eset) %in% genes.diff,]

      ##
      p_val <- sapply(1:nrow(count.use), function(x){
        wilcox.test(count.use[x,] ~ ct.group.temp)$p.value
      })
      p_val_adj <- p.adjust(p = p_val, method = "bonferroni",
                            n = nrow(count.use))
      markers.temp <- rownames(count.use)[p_val_adj < 0.05]
      markers.wilcox <- c(markers.wilcox, markers.temp)
    }
    markers <- unique(markers.wilcox)
    message("Selected ",length(markers), " marker genes by Wilcoxon test...")
  } # else need input of marker genes for clustering

  # match genes / cells first
  if (weight.basis){
    basis <- sc.basis$basis.mvw
    basis.fl <- sc.fl.basis$basis.mvw
  } else {
    basis <- sc.basis$basis
    basis.fl <- sc.fl.basis$basis
  }
  if (!is.null(markers)){
    commongenes <- Reduce(intersect, list(rownames(basis), rownames(bulk.eset), markers))
    commongenes.fl <- Reduce(intersect, list(rownames(basis.fl), rownames(bulk.eset), markers))
  } else {
    commongenes <- intersect(rownames(basis), rownames(bulk.eset))
    commongenes.fl <- intersect(rownames(basis.fl), rownames(bulk.eset))
    # stop when few common genes exist...
    if (length(commongenes) < 0.2 * min(dim(sc.eset)[1], dim(bulk.eset)[1])){
      stop('Too few common genes!')
    }
  }

  message(paste("Used", length(commongenes), "common genes for all cell types, \n",
                "Used", length(commongenes.fl), "common genes for first level cell types..."))

  basis.mvw <- basis[commongenes, ct.sub]
  basis.mvw.fl <- basis.fl[commongenes.fl, ct.fl.sub]

  xbulk0 <- getCPM0(exprs(bulk.eset)[commongenes,])
  xbulk <- as.matrix(xbulk0) ## whether to normalize all /common genes
  colnames(xbulk) <- colnames(bulk.eset)
  xbulk1 <- getCPM0(exprs(bulk.eset)[commongenes.fl,])
  xbulk.fl <- as.matrix(xbulk1)

  ALS.S <- sc.basis$sum.mat[ct.sub]
  N.bulk <- ncol(bulk.eset)
  valid.ct <- (colSums(is.na(basis.mvw)) == 0) & (!is.na(ALS.S))
  ALS.S.fl <- sc.fl.basis$sum.mat[ct.fl.sub]
  valid.ct.fl <- (colSums(is.na(basis.mvw.fl)) == 0) & (!is.na(ALS.S.fl))

  if (sum(valid.ct) <= 1) {
    stop("Not enough valid cell type!")
  }
  message(paste("Used", sum(valid.ct), "cell types in deconvolution...\n",
                "Used", sum(valid.ct.fl),"first level cell types ..."))

  basis.mvw <- basis.mvw[, valid.ct]
  ALS.S <- ALS.S[valid.ct]
  basis.mvw.fl <- basis.mvw.fl[, valid.ct.fl]
  ALS.S.fl <- ALS.S[valid.ct.fl]

  prop.est <- NULL
  rsquared <- NULL

  # prop estimation for each bulk sample:
  for (i in 1:N.bulk) {
    # i=1
    xbulk.temp <- xbulk[, i] ## *1e3  will affect a little bit
    message(paste(colnames(xbulk)[i], "has common genes", sum(xbulk[, i] != 0), "..."))
    if (allgenes.fl){
      markers.fl <- names(xbulk.temp)
    } else {
      markers.fl <- Reduce(intersect, list(markers, names(xbulk.temp)))
    }

    # first level NNLS:
    lm <- nnls::nnls(A=basis.mvw.fl[markers.fl,],b=xbulk.temp[markers.fl])
    delta <- lm$residuals
    wt.gene <- 1/(nu + delta^2)
    x.wt <- xbulk.temp[markers.fl] *sqrt(wt.gene)
    b.wt <- sweep(basis.mvw.fl[markers.fl,],1,sqrt(wt.gene),"*")

    lm.wt <- nnls::nnls(A=b.wt, b=x.wt)
    prop.wt.fl <- lm.wt$x/sum(lm.wt$x)
    delta <- lm.wt$residuals

    for (iter in 1:iter.max){
      wt.gene <- 1/(nu + delta^2)
      x.wt <- xbulk.temp[markers.fl] * sqrt(wt.gene)
      b.wt <- sweep(basis.mvw.fl[markers.fl,],1,sqrt(wt.gene),"*")
      lm.wt <- nnls::nnls(A=b.wt, b=x.wt)
      delta.new <- lm.wt$residuals
      prop.wt.fl.new <- lm.wt$x/sum(lm.wt$x)

      if (sum(abs(prop.wt.fl.new - prop.wt.fl)) < epsilon){
        prop.wt.fl <- prop.wt.fl.new
        delta <- delta.new
        message("WNNLS for First level clusters Converged at iteration ", iter)
        break
      }
      prop.wt.fl <- prop.wt.fl.new
      delta <- delta.new
    }
    names(prop.wt.fl) <- colnames(basis.mvw.fl)

    # relationship between first level and overall
    rt <- table(sc.eset@phenoData@data[,ct.varname], sc.eset@phenoData@data[,fl.varname])
    rt <- rt[,ct.fl.sub]
    rt.list <- list()
    prop.wt <- NULL

    # prop.wt
    for (j in 1:ncol(rt)){ # for each first level cluster
      # j=1
      rt.list[[j]] <- rownames(rt)[rt[,j] >0]
      names(rt.list)[j] <- colnames(rt)[j]
      sub.cl <- rownames(rt)[rt[,j] >0]
      if (length(sub.cl) > 1 & prop.wt.fl[colnames(rt)[j]] > 0) {
        if (is.null(dim(prop.wt.fl))){
          # specify genes in xbulk.j ... first level genes ...
          xbulk.j <- basis.mvw.fl[,j]*prop.wt.fl[j] + (xbulk.temp - basis.mvw.fl %*% lm.wt$x)*prop.wt.fl[j]
        } else {
          xbulk.j <- basis.mvw.fl[,j]*prop.wt.fl[,j] + (xbulk.temp - basis.mvw.fl %*% lm.wt$x)*prop.wt.fl[,j]
        }

        markers.sl <- Reduce(intersect, list(markers, rownames(xbulk.j)))

        ##############################################################################
        # make markers.sub as a list, for each of the first-level intra clusters.
        ##############################################################################

        basis.sl <- basis.mvw[markers.sl,rownames(rt)[rt[,j] >0]]
        lm.sl <- nnls::nnls(A=basis.sl,b=xbulk.j[markers.sl,])
        delta.sl <- lm.sl$residuals
        wt.gene.sl <- 1/(nu + delta.sl^2)
        x.wt.sl <- xbulk.j[markers.sl,]*sqrt(wt.gene.sl)
        b.wt.sl <- sweep(basis.sl,1,sqrt(wt.gene.sl),"*")

        lm.wt.sl <- nnls::nnls(A=b.wt.sl, b=x.wt.sl)
        prop.wt.sl <- lm.wt.sl$x/sum(lm.wt.sl$x)
        delta.sl <- lm.wt.sl$residuals

        for (iter in 1:iter.max){
          wt.gene.sl <- 1/(nu + delta.sl^2)
          x.wt.sl <- xbulk.j[markers.sl,] * sqrt(wt.gene.sl)
          b.wt.sl <- sweep(basis.sl,1,sqrt(wt.gene.sl),"*")
          lm.wt.sl <- nnls::nnls(A=b.wt.sl, b=x.wt.sl)
          delta.sl.new <- lm.wt.sl$residuals
          prop.wt.sl.new <- lm.wt.sl$x/sum(lm.wt.sl$x)

          if (sum(abs(prop.wt.sl.new - prop.wt.sl)) < epsilon){
            prop.wt.sl <- prop.wt.sl.new
            delta.sl <- delta.sl.new
            cat("WNNLS for Second level clusters",rownames(rt)[rt[,j] >0],"Converged at iteration ", iter)
            break
          }
          prop.wt.sl <- prop.wt.sl.new
          delta.sl <- delta.sl.new
        }
        names(prop.wt.sl) <- sub.cl
        prop.wt <- c(prop.wt, prop.wt.sl*prop.wt.fl[colnames(rt)[j]])
      } else if (length(sub.cl) == 1){
        # j=2
        prop.wt <- c(prop.wt, prop.wt.fl[colnames(rt)[j]])
      } else if (length(sub.cl) > 1 & prop.wt.fl[colnames(rt)[j]] == 0){
        prop.wt.sl <- rep(0, length(sub.cl))
        names(prop.wt.sl) <- sub.cl
        prop.wt <- c(prop.wt, prop.wt.sl)
      }

    }
    prop.est <- rbind(prop.est, prop.wt)
  }
  rownames(prop.est) <- colnames(bulk.eset)

  peval <- NULL
  if (!is.null(truep)){
    peval <- SCDC_peval(ptrue= truep, pest = prop.est, pest.names = c('SCDC'),
                       select.ct = ct.sub)
  }

  return(list(prop.est = prop.est, prop.wt.fl = prop.wt.fl, basis.mvw = basis.mvw, peval = peval,
              sc.basis = sc.basis, sc.fl.basis = sc.fl.basis))
}


############################################
#' Tree-guided proportion estimation for ONE subject
#' @description Proportion estimation function for ONE-subject case, and apply tree-guided deconvolution
#' @name SCDC_prop_ONE_subcl_marker
#' @param bulk.eset ExpressionSet object for bulk samples
#' @param sc.eset ExpressionSet object for single cell samples
#' @param ct.varname variable name for 'cell types'
#' @param fl.varname variable name for first-level 'meta-clusters'
#' @param sample variable name for subject/samples
#' @param ct.sub a subset of cell types that are selected to construct basis matrix
#' @param ct.fl.sub 'cell types' for first-level 'meta-clusters'
#' @param iter.max the maximum number of iteration in WNNLS
#' @param nu a small constant to facilitate the calculation of variance
#' @param epsilon a small constant number used for convergence criteria
#' @param weight.basis logical, use basis matrix adjusted by MVW, default is T.
#' @param select.marker logical, select marker genes to perform deconvolution in tree-guided steps. Default is T.
#' @param markers A set of marker gene that input manully to be used in deconvolution. If NULL, then
#' @param marker.varname variable name of cluster groups when selecting marker genes. If NULL, then use ct.varname.
#' @param allgenes.fl logical, use all genes in the first-level deconvolution
#' @param pseudocount.use a constant number used when selecting marker genes, default is 1.
#' @param LFC.lim a threshold of log fold change when selecting genes as input to perform Wilcoxon's test.
#' @param truep true cell-type proportions for bulk samples if known
#' @return Estimated proportion, basis matrix, predicted gene expression levels for bulk samples
#' @export
SCDC_prop_ONE_subcl_marker <- function(bulk.eset, sc.eset, ct.varname, fl.varname, sample, truep = NULL,
                                       ct.sub = NULL, ct.fl.sub, iter.max = 3000, nu = 1e-04, epsilon = 0.001,
                                       weight.basis = F, bulk_disease = NULL, select.marker = T, markers = NULL, marker.varname = NULL,
                                       pseudocount.use = 1, LFC.lim = 0.5, allgenes.fl = F,
                                       ...)
{

  if (is.null(ct.sub)){
    ct.sub <- unique(sc.eset@phenoData@data[,ct.varname])[!is.na(unique(sc.eset@phenoData@data[,ct.varname]))]
  }
  ct.sub <- ct.sub[!is.na(ct.sub)]
  ct.fl.sub <- ct.fl.sub[!is.na(ct.fl.sub)]
  bulk.eset <- bulk.eset[rowSums(exprs(bulk.eset))>0, , drop = FALSE]
  sc.basis <- SCDC_basis_ONE(x = sc.eset, ct.sub = ct.sub, ct.varname = ct.varname, sample = sample)
  sc.fl.basis <- SCDC_basis_ONE(x = sc.eset, ct.sub = ct.fl.sub[!is.na(ct.fl.sub)],
                                ct.varname = fl.varname, sample = sample)
  if (select.marker){
    if (is.null(marker.varname)){
      marker.varname <- ct.varname
    }
    # wilcox test on two groups of cells for marker gene selection... (refer to seurat::FindMarkers)
    countmat <- exprs(sc.eset)
    ct.group <- sc.eset@phenoData@data[,marker.varname]
    markers.wilcox <- NULL
    # u=1
    for(u in 1:length(unique(ct.group))){
      ct.group.temp <- ct.group == unique(ct.group)[u]
      group.1 <- apply(X = countmat[,ct.group.temp],
                       MARGIN = 1, FUN = function(x) log(x = mean(x = expm1(x = x)) +
                                                           pseudocount.use))
      group.2 <- apply(X = countmat[,! ct.group.temp],
                       MARGIN = 1, FUN = function(x) log(x = mean(x = expm1(x = x)) +
                                                           pseudocount.use))
      genes.diff <- rownames(sc.eset)[(group.1 - group.2) > LFC.lim]
      count.use <- countmat[rownames(sc.eset) %in% genes.diff,]

      ##
      p_val <- sapply(1:nrow(count.use), function(x){
        wilcox.test(count.use[x,] ~ ct.group.temp)$p.value
      })
      p_val_adj <- p.adjust(p = p_val, method = "bonferroni",
                            n = nrow(count.use))
      markers.temp <- rownames(count.use)[p_val_adj < 0.05]
      markers.wilcox <- c(markers.wilcox, markers.temp)
    }
    markers <- unique(markers.wilcox)
    message("Selected ",length(markers), " marker genes by Wilcoxon test...")
  } # else need input of marker genes for clustering

  # match genes / cells first
  if (weight.basis){
    basis <- sc.basis$basis.mvw
    basis.fl <- sc.fl.basis$basis.mvw
  } else {
    basis <- sc.basis$basis
    basis.fl <- sc.fl.basis$basis
  }
  if (!is.null(markers)){
    commongenes <- Reduce(intersect, list(rownames(basis), rownames(bulk.eset), markers))
    commongenes.fl <- Reduce(intersect, list(rownames(basis.fl), rownames(bulk.eset), markers))
  } else {
    commongenes <- intersect(rownames(basis), rownames(bulk.eset))
    commongenes.fl <- intersect(rownames(basis.fl), rownames(bulk.eset))
    # stop when few common genes exist...
    if (length(commongenes) < 0.2 * min(dim(sc.eset)[1], dim(bulk.eset)[1])){
      stop('Too few common genes!')
    }
  }

  message(paste("Used", length(commongenes), "common genes for all cell types, \n",
                "Used", length(commongenes.fl), "common genes for first level cell types..."))

  basis.mvw <- basis[commongenes, ct.sub]
  basis.mvw.fl <- basis.fl[commongenes.fl, ct.fl.sub]

  xbulk0 <- getCPM0(exprs(bulk.eset)[commongenes,])
  xbulk <- as.matrix(xbulk0) ## whether to normalize all /common genes
  colnames(xbulk) <- colnames(bulk.eset)
  xbulk1 <- getCPM0(exprs(bulk.eset)[commongenes.fl,])
  xbulk.fl <- as.matrix(xbulk1)

  ALS.S <- sc.basis$sum.mat[ct.sub]
  N.bulk <- ncol(bulk.eset)
  valid.ct <- (colSums(is.na(basis.mvw)) == 0) & (!is.na(ALS.S))

  ALS.S.fl <- sc.fl.basis$sum.mat[ct.fl.sub]
  valid.ct.fl <- (colSums(is.na(basis.mvw.fl)) == 0) & (!is.na(ALS.S.fl))

  if (sum(valid.ct) <= 1) {
    stop("Not enough valid cell type!")
  }
  message(paste("Used", sum(valid.ct), "cell types in deconvolution...\n",
                "Used", sum(valid.ct.fl),"first level cell types ..."))

  basis.mvw <- basis.mvw[, valid.ct]
  ALS.S <- ALS.S[valid.ct]

  basis.mvw.fl <- basis.mvw.fl[, valid.ct.fl]
  ALS.S.fl <- ALS.S[valid.ct.fl]

  prop.est <- NULL
  rsquared <- NULL

  # prop estimation for each bulk sample:
  for (i in 1:N.bulk) {
    # i=1
    xbulk.temp <- xbulk[, i] *1e3 ## will affect a little bit
    message(paste(colnames(xbulk)[i], "has common genes", sum(xbulk[, i] != 0), "..."))
    if (allgenes.fl){
      markers.fl <- names(xbulk.temp)
    } else {
      markers.fl <- Reduce(intersect, list(markers, names(xbulk.temp)))
    }

    # first level NNLS:
    lm <- nnls::nnls(A=basis.mvw.fl[markers.fl,],b=xbulk.temp[markers.fl])
    delta <- lm$residuals
    wt.gene <- 1/(nu + delta^2)
    x.wt <- xbulk.temp[markers.fl] *sqrt(wt.gene)
    b.wt <- sweep(basis.mvw.fl[markers.fl,],1,sqrt(wt.gene),"*")

    lm.wt <- nnls::nnls(A=b.wt, b=x.wt)
    prop.wt.fl <- lm.wt$x/sum(lm.wt$x)
    delta <- lm.wt$residuals

    for (iter in 1:iter.max){
      wt.gene <- 1/(nu + delta^2)
      x.wt <- xbulk.temp[markers.fl] * sqrt(wt.gene)
      b.wt <- sweep(basis.mvw.fl[markers.fl,],1,sqrt(wt.gene),"*")
      lm.wt <- nnls::nnls(A=b.wt, b=x.wt)
      delta.new <- lm.wt$residuals
      prop.wt.fl.new <- lm.wt$x/sum(lm.wt$x)

      if (sum(abs(prop.wt.fl.new - prop.wt.fl)) < epsilon){
        prop.wt.fl <- prop.wt.fl.new
        delta <- delta.new
        message("WNNLS for First level clusters Converged at iteration ", iter)
        break
      }
      prop.wt.fl <- prop.wt.fl.new
      delta <- delta.new
    }
    names(prop.wt.fl) <- colnames(basis.mvw.fl)

    # relationship between first level and overall
    rt <- table(sc.eset@phenoData@data[,ct.varname], sc.eset@phenoData@data[,fl.varname])
    rt <- rt[,ct.fl.sub]
    rt.list <- list()
    prop.wt <- NULL

    # prop.wt
    for (j in 1:ncol(rt)){ # for each first level cluster
      # j=1
      rt.list[[j]] <- rownames(rt)[rt[,j] >0]
      names(rt.list)[j] <- colnames(rt)[j]
      sub.cl <- rownames(rt)[rt[,j] >0]
      if (length(sub.cl) > 1 & prop.wt.fl[colnames(rt)[j]] > 0) {
        if (is.null(dim(prop.wt.fl))){
          # specify genes in xbulk.j??? first level genes?
          xbulk.j <- basis.mvw.fl[,j]*prop.wt.fl[j] + (xbulk.temp - basis.mvw.fl %*% lm.wt$x)*prop.wt.fl[j]
        } else {
          xbulk.j <- basis.mvw.fl[,j]*prop.wt.fl[,j] + (xbulk.temp - basis.mvw.fl %*% lm.wt$x)*prop.wt.fl[,j]
        }

        markers.sl <- Reduce(intersect, list(markers, rownames(xbulk.j)))

        ##############################################################################
        # make markers.sub as a list, for each of the first-level intra clusters.
        ##############################################################################

        basis.sl <- basis.mvw[markers.sl,rownames(rt)[rt[,j] >0]]
        lm.sl <- nnls::nnls(A=basis.sl,b=xbulk.j[markers.sl,])
        delta.sl <- lm.sl$residuals
        wt.gene.sl <- 1/(nu + delta.sl^2)
        x.wt.sl <- xbulk.j[markers.sl,]*sqrt(wt.gene.sl)
        b.wt.sl <- sweep(basis.sl,1,sqrt(wt.gene.sl),"*")

        lm.wt.sl <- nnls::nnls(A=b.wt.sl, b=x.wt.sl)
        prop.wt.sl <- lm.wt.sl$x/sum(lm.wt.sl$x)
        delta.sl <- lm.wt.sl$residuals

        for (iter in 1:iter.max){
          wt.gene.sl <- 1/(nu + delta.sl^2)
          x.wt.sl <- xbulk.j[markers.sl,] * sqrt(wt.gene.sl)
          b.wt.sl <- sweep(basis.sl,1,sqrt(wt.gene.sl),"*")
          lm.wt.sl <- nnls::nnls(A=b.wt.sl, b=x.wt.sl)
          delta.sl.new <- lm.wt.sl$residuals
          prop.wt.sl.new <- lm.wt.sl$x/sum(lm.wt.sl$x)

          if (sum(abs(prop.wt.sl.new - prop.wt.sl)) < epsilon){
            prop.wt.sl <- prop.wt.sl.new
            delta.sl <- delta.sl.new
            cat("WNNLS for Second level clusters",rownames(rt)[rt[,j] >0],"Converged at iteration ", iter)
            break
          }
          prop.wt.sl <- prop.wt.sl.new
          delta.sl <- delta.sl.new
        }
        names(prop.wt.sl) <- sub.cl
        prop.wt <- c(prop.wt, prop.wt.sl*prop.wt.fl[colnames(rt)[j]])
      } else if (length(sub.cl) == 1){
        # j=2
        prop.wt <- c(prop.wt, prop.wt.fl[colnames(rt)[j]])
      } else if (length(sub.cl) > 1 & prop.wt.fl[colnames(rt)[j]] == 0){
        prop.wt.sl <- rep(0, length(sub.cl))
        names(prop.wt.sl) <- sub.cl
        prop.wt <- c(prop.wt, prop.wt.sl)
      }

    }
    prop.est <- rbind(prop.est, prop.wt)
  }
  rownames(prop.est) <- colnames(bulk.eset)

  peval <- NULL
  if (!is.null(truep)){
    peval <- SCDC_eval(ptrue= truep, pest = prop.est, pest.names = c('SCDC'),
                       dtname = 'Perou', select.ct = ct.sub, bulk_obj = bulk.eset,
                       bulk_disease = bulk_disease)
  }

  return(list(prop.est = prop.est, prop.wt.fl = prop.wt.fl, basis.mvw = basis.mvw, peval = peval,
              sc.basis = sc.basis, sc.fl.basis = sc.fl.basis))
}
