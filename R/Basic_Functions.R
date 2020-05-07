####################################
########## basic functions #########
####################################

#######################################################
#' Normalize matrix / vector
#' @description Normalize matrix by column or normalize a vector
#' @name getCPM0
#' @param x a matrix or a vector
#' @export
getCPM0 <- function(x, verbose = F){
  if (is.null(dim(x))){
    if (verbose){
      message("Normalizing a vector instead of a matrix")
    }
    vec = as.matrix(x/sum(x))
    vec
  } else {
    cpm <- t(t(x)/apply(x,2,sum))
    cpm
  }
}


#######################################################
#' Get ExpressionSet
#' @description Use Pdata, Fdata, and count matrix to derive ExpressionSet Object
#' @name getESET
#' @import Biobase
#' @param exprs raw count matrix
#' @param fdata feature data, for genes, usually it's gene name
#' @param pdata pheno data, for samples, usually it's the characteristics for each single cell/bulk sample, including name, gender, age, cluster, disease,...
#' @export
library(Biobase)
getESET <- function(exprs, fdata, pdata){
  pdata <- as.data.frame(pdata)
  fdata <- as.data.frame(fdata)
  exprs <- as.matrix(exprs)
  rownames(pdata) <- colnames(exprs)
  rownames(fdata) <- rownames(exprs)
  eset <- ExpressionSet(exprs,
                        AnnotatedDataFrame(pdata),
                        AnnotatedDataFrame(fdata))
}

######################################################
#' Construct Pseudo bulk samples
#' @description Construct Pseudo bulk samples by actual number of cells per subject
#' @name generateBulk_allcells
#' @import xbioc
#' @param eset ExpressionSet object for single cells
#' @param ct.varname variable name for 'cell types'
#' @param sample variable name for subject/samples
#' @param disease indicate the health condition of subjects
#' @param ct.sub a subset of cell types that are selected to construct pseudo bulk samples. If NULL, then all cell types are used.
#' @return pseudo bulk samples ExpressionSet, and actual cell-type proportions
#' @export
generateBulk_allcells <- function(eset, ct.varname, sample, disease = NULL, ct.sub = NULL){
  if (is.null(ct.sub)){
    ct.sub <- unique(eset@phenoData@data[,ct.varname])
  }
  eset <- eset[, eset@phenoData@data[,ct.varname] %in% ct.sub]
  cluster.id <- eset@phenoData@data[,ct.varname]
  sample.id <- eset@phenoData@data[,sample]
  condition.id <- eset@phenoData@data[,disease]

  ## expression
  pseudo.exprs <- sapply(unique(sample.id), function(sid){
    y <- exprs(eset)[, sample.id %in% sid]
    rowSums(y, na.rm = T)
  })
  colnames(pseudo.exprs) <- unique(sample.id)
  ## true proportion: sample by cell types
  ncount <- table(sample.id, cluster.id)
  true.prop <- ncount / rowSums(ncount, na.rm = T)
  true.prop <- true.prop[complete.cases(true.prop),]
  ## eset for pseudo bulk sample
  if (is.null(disease)){
    pseudo.disease <- NA
  } else {
    pseudo.disease <- sapply(unique(sample.id), function(sid){
      condition.id[sample.id == sid][1]
    })
  }
  pseudo.pdata <- data.frame(sample = colnames(pseudo.exprs),
                             disease = pseudo.disease)
  pseudo.fdata <- data.frame(genes = rownames(pseudo.exprs))
  rownames(pseudo.fdata) <- rownames(pseudo.exprs)
  pseudo_eset <- getESET(exprs = pseudo.exprs,
                         fdata = pseudo.fdata,
                         pdata = pseudo.pdata)
  return(list(truep = true.prop, pseudo_eset = pseudo_eset))
}



######################################################
#' Construct Pseudo bulk samples -- random
#' @description Construct Pseudo bulk samples by random sampled number of cells per subject, but based on the actual numbers.
#' @name generateBulk_norep
#' @param eset ExpressionSet object for single cells
#' @param ct.varname variable name for 'cell types'
#' @param sample variable name for subject/samples
#' @param disease indicate the health condition of subjects
#' @param ct.sub a subset of cell types that are selected to construct pseudo bulk samples. If NULL, then all cell types are used.
#' @param prop_mat manually input the cell-type proportion for pseudo bulk samples
#' @param nbulk number of pseudo bulk samples to be constructed
#' @param samplewithRep logical, randomly sample single cells with replacement. Default is F.
#' @return pseudo bulk samples ExpressionSet, and actual cell-type proportions
#' @export
generateBulk_norep <- function(eset, ct.varname, sample, disease = NULL, ct.sub,
                               prop_mat = NULL, nbulk=10, samplewithRep = F){
  x.sub <- eset[,eset@phenoData@data[,ct.varname] %in% ct.sub]
  # qc: remove non-zero genes
  x.sub <- x.sub[rowSums(exprs(x.sub)) > 0,]
  # calculate sample mean & sample variance matrix: genes by cell types
  ct.id <- droplevels(as.factor(x.sub@phenoData@data[,ct.varname]))
  sample.id <- x.sub@phenoData@data[,sample]

  pdatab <- x.sub@phenoData@data
  if (! is.null(disease)){
    disease <- x.sub@phenoData@data[,disease]
    names(disease) <- sample.id
  }

  # number of cell type of interest
  k <- length(unique(ct.id))
  message(paste('Using',k,'cell types to generate pseudo bulk samples...'))
  # select donors for each pseudo bulk sample
  pseudo_donors <- sample(sample.id, nbulk, replace = T)
  names(pseudo_donors) <- paste("bulk",1:nbulk, sep = "_")

  # generate random matrix for true proportions
  if (!is.null(prop_mat)){
    true.p1 <- prop_mat # manually input proportion matrix...
    colnames(true.p1) <- unique(ct.id)[order(unique(ct.id))]
    rownames(true.p1) <- names(pseudo_donors)
    message("Using input proportion matrix to create pseudo bulk samples...")
    true.ct <- matrix(data = 0,ncol = k, nrow = nbulk) # true number of cells per cell type for each sample
    colnames(true.ct) <- unique(ct.id)[order(unique(ct.id))]
    rownames(true.ct) <- paste(pseudo_donors,1:nbulk,sep = "_")
    # make sure if without replacement, number of cells matches the input prop mat...
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  } else {
    true.p1 <- matrix(data = 0,ncol = k, nrow = nbulk)
    colnames(true.p1) <- unique(ct.id)[order(unique(ct.id))]
    rownames(true.p1) <- paste(pseudo_donors,1:nbulk,sep = "_")
    true.ct <- matrix(data = 0,ncol = k, nrow = nbulk) # true number of cells per cell type for each sample
    colnames(true.ct) <- unique(ct.id)[order(unique(ct.id))]
    rownames(true.ct) <- paste(pseudo_donors,1:nbulk,sep = "_")
    message("Generating random cell type proportions...")
  }

  # create pseudo bulk sample.id according to total available number of cells:
  pseudo_bulk <- NULL
  for(xx in 1:length(pseudo_donors)){ #length(pseudo_donors)
    # xx = 1
    message('generating bulk ', xx, ' from donor ', pseudo_donors[xx], '...')
    idxd <- sample.id == pseudo_donors[xx] # for a selected donor
    temp <- exprs(x.sub)[,idxd] # his expression matrix
    temp.cluster <- ct.id[idxd] # cluster info for cells

    # match names!!!!!!!!!!!!!!!!!
    temp.ncellk <- table(factor(temp.cluster))

    if (is.null(prop_mat)){ #if using random proportions
      temp.nct <- ceiling(runif(length(temp.ncellk), min = 0.6, max = 1)*temp.ncellk) # take random number of available single cells
      true.p1[xx,names(temp.nct)] <- temp.nct/sum(temp.nct) # true proportions
      true.ct[xx,] <- temp.nct[colnames(true.ct)] # true number of cells in the pseudo bulk.
      true.p1[is.na(true.p1)] <- 0
      true.ct[is.na(true.ct)] <- 0
    } else { #if using the user-defined proportions
      temp.ntotal <- min(temp.ncellk / true.p1[xx,], na.rm = T)
      if (temp.ntotal <= k){
        message("Please check if your input prop_mat is reasonable. The number of cells of certain selected cell type might be too small.")
      }
      true.ct[xx,] <- round(temp.ntotal*true.p1[xx,names(temp.ncellk)] ) # true number of cells in the pseudo bulk.
      true.p1[is.na(true.p1)] <- 0
      true.ct[is.na(true.ct)] <- 0
    }

    temp.b1 <- sapply(ct.sub, function(ucluster){
      temp.vec <- temp[,temp.cluster %in% ucluster] # for a specific cell type
      if (is.null(dim(temp.vec))){
        temp.sum <- rep(0, length(temp.vec))
      } else if (! is.null(dim(temp.vec))) {
        if (dim(temp.vec)[2] == 0){
          temp.sum <- rep(0, dim(temp.vec)[1])
        } else {
          temp.sample <- sample(1:ncol(temp.vec), true.ct[xx, ucluster], replace = samplewithRep) # which cells in this cell type will be selected
          temp.mat <- temp.vec[,temp.sample] # select all those cells
          if (is.null(dim(temp.mat))){
            temp.sum <- temp.mat
          } else {
            temp.sum <- rowSums(temp.mat, na.rm = T) # expression sum for one cell type in this bulk, need to sum up all types.
          }
        }
      }
    })

    out = rowSums(temp.b1)
    pseudo_bulk <- cbind(pseudo_bulk, out)
    colnames(pseudo_bulk)[xx] <- paste(pseudo_donors[xx],xx,sep = "_")
  }
  # create pseudo eset for bulk sample.id
  if (! is.null(disease)){
    pseudo_pdata <- data.frame(sample.id = colnames(pseudo_bulk), donors = pseudo_donors, disease = disease[pseudo_donors])
  } else {
    pseudo_pdata <- data.frame(sample.id = colnames(pseudo_bulk), donors = pseudo_donors)
  }
  rownames(pseudo_pdata) <- colnames(pseudo_bulk)
  pseudo_fdata <- data.frame(labelDescription = rownames(pseudo_bulk),
                             row.names = rownames(pseudo_bulk))
  message("generating expression set object for pseudo bulk sample.id...")
  pseudo_eset <- ExpressionSet(pseudo_bulk,
                               AnnotatedDataFrame(pseudo_pdata),
                               AnnotatedDataFrame(pseudo_fdata))
  pseudo_eset0 <- pseudo_eset[,rowSums(!is.na(true.p1))>0] #non-zero eset
  true.p0 <- true.p1[rowSums(!is.na(true.p1))>0,]
  true.ct0 <- true.ct[rowSums(!is.na(true.p1))>0,]
  return(list(true_p = true.p1, pseudo_bulk = pseudo_bulk, pseudo_eset = pseudo_eset,
              num.real = true.ct, true_p0 = true.p0, true.ct0 = true.ct0, pseudo_eset0 = pseudo_eset0)) # , entropy = entropy
}




#' Grid search matrix
#' @description Generate Grid-search matrix
#' @name getSearchGrid
#' @param lengthby search length
#' @param nparam number of paramters/reference datasets
#' @export
getSearchGrid <- function(lengthby, nparam){
  wlist <- list()
  for(i in 1:nparam){
    wlist[[i]] <- seq(0,1+lengthby,by = lengthby)
  }
  w1grid <- round(expand.grid(wlist),digits = 2) # without the round function, there's missing combinations
  w1grid <- w1grid[round(rowSums(w1grid), digits = 2)==1.00,]
  colnames(w1grid) <- paste("w",1:nparam, sep = "")
  return(w1grid)
}

###############################################
#' Evaluate deconvolved proportions
#' @description Evaluation function, for deconvolved proportions and the actual/true proportions
#' @name SCDC_peval
#' @param ptrue a matrix of true/actual cell-type proportions for bulk samples
#' @param pest a list of estimated cell-type proportion matrix
#' @param pest.names method name for the estimated proportions in the pest list
#' @param select.ct selected cell types for deconvolution
#' @export
SCDC_peval <- function(ptrue, pest, pest.names, select.ct = NULL){
  if (!is.list(pest)){
    pest <- list(pest)
  }
  if (!is.data.frame(ptrue)){
  ptrue <- as.data.frame.matrix(ptrue)
  }
  n_est <- length(pest)
  sample_names <- lapply(pest, rownames)
  ctype_names <- lapply(pest, colnames)
  sample_common <- Reduce(intersect, sample_names)
  ctype_common <- Reduce(intersect, ctype_names)
  celltype <- intersect(colnames(ptrue), ctype_common)
  if (!is.null(select.ct)) {
    celltype <- intersect(celltype, select.ct)
  }
  sample <- intersect(rownames(ptrue), sample_common)
  N <- length(sample)
  K <- length(celltype)
  if (N < 1) {
    stop("No common Subjects! Check rowname!")
  }
  if (K <= 1) {
    stop("Not enough cell types!")
  }
  ptrue.use <- ptrue[intersect(rownames(ptrue), sample), intersect(colnames(ptrue), celltype)]
  ptrue.use <- as.data.frame.matrix(ptrue.use / apply(ptrue.use,1,sum))
  ptrue.use[is.na(ptrue.use)] <- 0

  # for each estimation method in the list
  evals <- lapply(pest, function(xx){
    pest.use <- xx[intersect(rownames(xx), sample), intersect(colnames(xx), celltype)]
    pest.use <- as.data.frame.matrix(pest.use / apply(pest.use,1,sum))
    pest.use <- pest.use[rownames(ptrue.use),colnames(ptrue.use)]
    RMSD_bysample <- round(sqrt(rowMeans((ptrue.use - pest.use)^2)), digits = 5)
    mAD_bysample <- round(rowMeans(abs(ptrue.use - pest.use)), digits = 5)
    Pearson_bysample <- sapply(1:nrow(ptrue.use), function(ss) {
	round(cor(c(as.matrix(ptrue.use[ss, ])), c(as.matrix(pest.use[ss, ]))), digits = 5)
	})
    RMSD <- round(sqrt(mean(as.matrix((ptrue.use - pest.use)^2), na.rm = T)), digits = 5)
    mAD <- round(mean(as.matrix(abs(ptrue.use - pest.use)), na.rm = T), digits = 5)
    Pearson <- round(cor(c(as.matrix(ptrue.use)), c(as.matrix(pest.use))), digits = 4)

    return(list(pest.use = pest.use, RMSD_bysample = RMSD_bysample, mAD_bysample = mAD_bysample, Pearson_bysample = Pearson_bysample,
                RMSD = RMSD, mAD = mAD, Pearson = Pearson))
  })
  evals.table <- NULL
  for (l in 1:n_est){
    evals.table <- rbind(evals.table, c(evals[[l]]$RMSD, evals[[l]]$mAD, evals[[l]]$Pearson))
  }
  colnames(evals.table) <- c("RMSD","mAD","R")
  rownames(evals.table) <- pest.names
  # evals per sample
  pearson.sample.table <- NULL
  for (l in 1:n_est){
    pearson.sample.table <- rbind(pearson.sample.table,evals[[l]]$Pearson_bysample)
  }
  rownames(pearson.sample.table) <- pest.names
  colnames(pearson.sample.table) <- rownames(ptrue.use)

  RMSD.sample.table <- NULL
  for (l in 1:n_est){
    RMSD.sample.table <- rbind(RMSD.sample.table,evals[[l]]$RMSD_bysample)
  }
  rownames(RMSD.sample.table) <- pest.names
  colnames(RMSD.sample.table) <- rownames(ptrue.use)

  mAD.sample.table <- NULL
  for (l in 1:n_est){
    mAD.sample.table <- rbind(mAD.sample.table,evals[[l]]$mAD_bysample)
  }
  rownames(mAD.sample.table) <- pest.names
  colnames(mAD.sample.table) <- rownames(ptrue.use)

  return(list(evals = evals, evals.table = evals.table, pearson.sample.table = pearson.sample.table,
              RMSD.sample.table = RMSD.sample.table, mAD.sample.table = mAD.sample.table))
}

###############################################
#' Evaluate predicted gene expression levels
#' @description Evaluation function, for observed and predicted gene expression levels
#' @name SCDC_yeval
#' @param y observed raw counts of a bulk sample
#' @param yest a list of predicted gene expression levels
#' @param yest.names method name for the predicted gene expression levels in yest list
#' @export
SCDC_yeval <- function(y, yest, yest.names=NULL){
  if (!is.list(yest)){
    yest <- list(yest)
  }
  n_est <- length(yest)
  # for each estimation method in the list
  evals <- lapply(yest, function(xx){
    if (!is.null(xx)){
      if (dim(xx)[2] >1){
        g.use = intersect(rownames(y), rownames(xx))
        y.norm <- getCPM0(y[g.use,])
        x = xx[g.use,colnames(y)]
        yest.norm = getCPM0(x)
        spearmany <-round(cor(c(yest.norm), c(y.norm), method = "spearman"), digits = 5)
        RMSDy = round(sqrt(mean((yest.norm - y.norm)^2)), digits = 7)
        mADy = round(mean(abs(yest.norm - y.norm)), digits = 8)

        spearmany_bysample <- sapply(1:ncol(x), function(ss) {round(cor(yest.norm[,ss],y.norm[,ss], method = "spearman"), digits = 4)})
        RMSDy_bysample = round(sqrt(colMeans((yest.norm[,] - y.norm[,])^2)), digits = 6)
        mADy_bysample = round(colMeans(abs(yest.norm - y.norm)), digits = 6)
      } else { # one sample case
        g.use = intersect(rownames(y), rownames(xx))
        y.norm <- getCPM0(y[g.use,])
        x = xx[g.use,colnames(y)]
        yest.norm = getCPM0(x)
        spearmany <-round(cor(c(yest.norm), c(y.norm), method = "spearman"), digits = 5)
        RMSDy = round(sqrt(mean((yest.norm - y.norm)^2)), digits = 7)
        mADy = round(mean(abs(yest.norm - y.norm)), digits = 8)

        spearmany_bysample <- spearmany
        RMSDy_bysample = RMSDy
        mADy_bysample = mADy
      }

    } else if (is.null(xx)){
      spearmany = NA
      RMSDy = NA
      mADy = NA
      spearmany_bysample = NA
      RMSDy_bysample = NA
      mADy_bysample =NA
    }
    return(list(spearmany = spearmany, RMSDy = RMSDy, mADy = mADy, spearmany_bysample = spearmany_bysample,
                RMSDy_bysample = RMSDy_bysample, mADy_bysample = mADy_bysample))
  })
  yevals.table <- NULL
  for (l in 1:n_est){
    yevals.table <- rbind(yevals.table, c(evals[[l]]$spearmany, evals[[l]]$RMSDy, evals[[l]]$mADy))
  }
  colnames(yevals.table) <- c("spearman_Y","RMSD_Y","mAD_Y")
  rownames(yevals.table) <- yest.names

  # evals per sample
  spearmany.sample.table <- NULL
  for (l in 1:n_est){
    spearmany.sample.table <- rbind(spearmany.sample.table,evals[[l]]$spearmany_bysample)
  }
  rownames(spearmany.sample.table) <- yest.names
  colnames(spearmany.sample.table) <- colnames(y)

  RMSDy.sample.table <- NULL
  for (l in 1:n_est){
    RMSDy.sample.table <- rbind(RMSDy.sample.table,evals[[l]]$RMSDy_bysample)
  }
  rownames(RMSDy.sample.table) <- yest.names
  colnames(RMSDy.sample.table) <- colnames(y)

  mADy.sample.table <- NULL
  for (l in 1:n_est){
    mADy.sample.table <- rbind(mADy.sample.table,evals[[l]]$mADy_bysample)
  }
  rownames(mADy.sample.table) <- yest.names
  colnames(mADy.sample.table) <- colnames(y)
  return(list(evals = evals, yevals.table = yevals.table, spearmany.sample.table =spearmany.sample.table,
              RMSDy.sample.table = RMSDy.sample.table, mADy.sample.table = mADy.sample.table))
}


###############################################
#' Calculate cell type proportions by suggested weights
#' @description Calculate proportions by linear combination of a list of proportions
#' @name wt_prop
#' @param wt a vector of weights for each reference dataset. should be of the same length as the proplist.
#' @param proplist a list of estimated proportions from each reference dataset
#' @export
wt_prop <- function(wt, proplist){
  wt <- as.numeric(wt)
  combo.list <- list()
  for (i in 1:length(proplist)){
    combo.list[[i]] <- proplist[[i]]*wt[i]
  }
  combo.prop <- Reduce("+", combo.list)
  return(combo.prop)
}

###############################################
#' Calculate cell type proportions by suggested weights
#' @description Calculate proportions by linear combination of a list of proportions
#' @name wt_y
#' @param wt a vector of weights for each reference dataset. Should be of the same length as the y.list.
#' @param y.list a matrix, with each column containing the vectorized estimated gene expressions from each reference dataset.
#' @export
wt_y <- function(wt, y.list = y.list){
  wt <- as.numeric(wt)
  combo.list <- list()
  for (i in 1:ncol(y.list)){
    combo.list[[i]] <- y.list[,i]*wt[i]
  }
  combo.y <- Reduce("+", combo.list)
  return(combo.y)
}


###############################################
#' Create user-defined SCDC_prop object, as the SCDC_ENSEMBLE input
#' @description Create user-defined SCDC_prop object, as the SCDC_ENSEMBLE input
#' @name CreateSCDCpropObj
#' @param p the estimated proportion matrix. sample by cell type.
#' @param b the estimated signature/feature/basis matrix. gene by cell type.
#' @export
CreateSCDCpropObj <- function(p,b){
  yhat <- b %*% t(p[,colnames(b)])
  obj <- list(prop.est.mvw = p,
              basis.mvw = b,
              yhat = yhat)
  return(obj)
}
