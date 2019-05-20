############### VISUALIZTION FUNCTIONS #######################
##############################################################
#' Demographic figure for single cells
#' @description Create a demographic figure for the input single cell dataset for each of the subject
#' @name DemoPlot
#' @param eset ExpressionSet object of the single cell data
#' @param cluster the variable name for "cell types" or "cluster name"
#' @param sample the variable name for subjects
#' @param select.ct the vector of cell types of interest, for example, c("alpha","beta")
#' @param Palette the color palette to be used in the demographic figure
#' @return A figure showing the number of cells and percentage of cell types for single cells after clustering
#' @export
DemoPlot <- function(eset, cluster, sample, select.ct, Palette = cbPalette){
  pdata <- eset@phenoData@data
  if (!is.factor(pdata[,sample])){
    pdata[,sample] <- as.factor(pdata[,sample])
  } else {
    pdata[,sample] <- droplevels(pdata[,sample])
  }
  levels(pdata[,cluster])[!levels(pdata[,cluster]) %in% select.ct] <- "Other"
  select.ct <- c(select.ct, "Other")
  count.tab <- as.data.frame(table(pdata[,cluster], pdata[, sample]))
  prop.tab <- as.data.frame(getCPM0(table(pdata[,cluster], pdata[, sample])))
  total.count <- as.matrix(table(pdata[,sample]))
  count.tab <- cbind(count.tab, total = total.count[match(count.tab$Var2, rownames(total.count)),1])
  prop.tab$type = "Percentage"
  prop.tab$total = 1
  count.tab$type = "Number of cells"

  dt <- rbind.data.frame(prop.tab, count.tab)
  p1<- ggplot(dt, aes(x=Var2, y=Freq, fill = factor(Var1, levels = select.ct))) +
    geom_bar(stat = "identity") +
    geom_text(data = count.tab ,aes(x= Var2, y= total+2, label = total)) +
    labs(title = "", fill = "Cell Type") + xlab("") + ylab("")+
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size=10),
          axis.text.y = element_text(size = 10),
          text = element_text(size = 10),
          plot.title = element_text(size=10, face = "bold"),
          plot.margin=unit(c(-5,1,-5,-3), "mm"),
          legend.position = "top", legend.title = element_blank(),
          legend.text = element_text(size=8),
          legend.box.spacing = unit(0, "mm")) +
    scale_fill_manual(values= Palette) +
    facet_grid(rows = vars(type), scales = "free")
  print(p1)
}


###############################################################
#' Heatmap for ENSEMBL weight selection, with four measures
#' @description Heatmap for ENSEMBL weight selection, with four measures
#' @name wheat_map
#' @param ensemble_res the EMSEMBL result object, derived from SCDC_ENSEMBL() or SCDC_ENSEMBL_subcl()
#' @param ref1 name for the first reference dataset
#' @param ref2 name for the second reference dataset
#' @return A figure of evaluated performance, varying the ENSEMBL weights for three reference datasets
wheat_map <- function(ensemble_res, ref1, ref2){
  pp <- ggplot( ensemble_res$gridres, aes(w1, w2, fill = Pearson)) +
    geom_tile(aes(fill = Pearson)) + #, colour = "white",size=0
    scale_fill_gradient(low = "#56B4E9", high = "#E69F00", name = "Pearson's R(prop)") +
    # xlab(paste("Weight for", ref1))+
    ylab(paste("Weight for", ref2)) +
    annotate(geom="text", x=mean(ensemble_res$gridres$w1[ensemble_res$gridres$Pearson == max(ensemble_res$gridres$Pearson)]),
             y=mean(ensemble_res$gridres$w2[ensemble_res$gridres$Pearson == max(ensemble_res$gridres$Pearson)]),
             label="*",
             color="#000000", size =5)+
    theme(text = element_text(size=10), plot.title = element_text(size=10),
          axis.title.x=element_blank(),
          axis.text.x = element_text(size=10, angle=45),
          axis.text.y = element_text(size=10), # axis.title.x = element_blank(),
          legend.position = c(0.56, 0.75)) #legend.position = "top"

  ps <- ggplot( ensemble_res$gridres, aes(w1, w2, fill = spearman_Y)) +
    geom_tile(aes(fill = spearman_Y) ) +
    scale_fill_gradient(low = "#56B4E9", high = "#E69F00", name = "Spearman's R(Y)") +
    # xlab(paste("Weight for", ref1))+
    # ylab(paste("Weight for", ref2)) +
    annotate(geom="text", x=mean(ensemble_res$gridres$w1[ensemble_res$gridres$spearman_Y==max(ensemble_res$gridres$spearman_Y)]),
             y=mean(ensemble_res$gridres$w2[ensemble_res$gridres$spearman_Y==max(ensemble_res$gridres$spearman_Y)]),
             label="*",
             color="#000000", size =5)+
    theme(text = element_text(size=10), plot.title = element_text(size=10),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x = element_text(size=10, angle=45),
          axis.text.y = element_text(size=10), # axis.title.x = element_blank(),
          legend.position = c(0.56, 0.75))

  pm <- ggplot( ensemble_res$gridres, aes(w1, w2, fill = mAD_Y)) +
    geom_tile(aes(fill = mAD_Y) ) +
    scale_fill_gradient(low = "#009E73", high = "#56B4E9", name = "mAD(Y)") +
    xlab(paste("Weight for", ref1))+
    ylab(paste("Weight for", ref2)) +
    annotate(geom="text", x=mean(ensemble_res$gridres$w1[ensemble_res$gridres$mAD_Y == min(ensemble_res$gridres$mAD_Y)]),
             y=mean(ensemble_res$gridres$w2[ensemble_res$gridres$mAD_Y == min(ensemble_res$gridres$mAD_Y)]),
             label="*",
             color="#000000", size =5)+
    theme(text = element_text(size=10), plot.title = element_text(size=10),
          axis.text.x = element_text(size=10, angle=45),
          axis.text.y = element_text(size=10), # axis.title.x = element_blank(),
          legend.position = c(0.65, 0.75))

  pr <- ggplot( ensemble_res$gridres, aes(w1, w2, fill = RMSD_Y)) +
    geom_tile(aes(fill = RMSD_Y) ) +
    scale_fill_gradient(low = "#009E73", high = "#56B4E9", name = "RMSD(Y)") +
    xlab(paste("Weight for", ref1))+
    # ylab(paste("Weight for", ref2)) +
    annotate(geom="text", x=mean(ensemble_res$gridres$w1[ensemble_res$gridres$RMSD_Y == min(ensemble_res$gridres$RMSD_Y)]),
             y=mean(ensemble_res$gridres$w2[ensemble_res$gridres$RMSD_Y == min(ensemble_res$gridres$RMSD_Y)]),
             label="*",
             color="#000000", size =5)+
    theme(text = element_text(size=10), plot.title = element_text(size=10),
          axis.title.y=element_blank(),
          axis.text.x = element_text(size=10, angle=45),
          axis.text.y = element_text(size=10), # axis.title.x = element_blank(),
          legend.position = c(0.65, 0.75))

  pout <- plot_grid(pp, ps, pm, pr, ncol = 2)
  return(pout)
}


################################################
