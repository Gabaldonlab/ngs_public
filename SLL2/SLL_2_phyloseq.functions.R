

# ****************************************************************************************************************** #
# ******************************************************* #
# Clustering and network analyses ####
# ******************************************************* #

# Based on the steps in the Bork Enterotype paper

# ************************************************************************************************ #
get_clusters <- function(tagl, tl, dist_meas, PS, glomTab, simulation=FALSE, clusSubset=FALSE, sampSize=NULL, 
                         iteration=NULL, saveClusPlots=FALSE, clusPlotsDir=NULL, nClustForce=NULL, CH.by.tagl=TRUE,
                         subPops=FALSE, subPop_Bdivs=NULL) {
  
  if (clusSubset) {
    clus_var_name <- sprintf("Stomatotype_subset_%s",iteration)
    s <- sample(sample_names(tagl), sampSize)
    ttc <- prune_samples(s, tagl)
  } else {
    clus_var_name <- "Stomatotype"
    s <- sample_names(tagl)
    ttc <- tagl
  }
  
  # tax.to.clust <- ttc@otu_table
  
  # Remove those genera for which the average abundance across all samples is below 0.01%
  noise.removal <- function(dataframe, percent=0.01, top=NULL){
    dataframe->Matrix
    bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent
    Matrix_1 <- Matrix[bigones,]
    # print(percent)
    return(Matrix_1)
  }
  
  # tax.to.clust <- noise.removal(tax.to.clust, percent=0.01)
  
  if (subPops==TRUE) {
    
    tax.to.clust <- ttc@otu_table
    dist.matrix <- subPop_Bdivs[[ dist_meas ]]
    
  } else {
    
    # Distance matrix object and abundance table (only different for JSD_alt)
    if (dist_meas == "JSD") dist.matrix = as.dist(jsd[ s, s ])
    else if (dist_meas == "Weighted_Unifrac") dist.matrix = as.dist(weighted_Unifrac[ s, s ])
    else if (dist_meas == "Unweighted_Unifrac") dist.matrix = as.dist(unweighted_Unifrac[ s, s ])
    else if (dist_meas == "VAW_GUnifrac") dist.matrix = as.dist(guni.VAW[ s, s ])
    else if (dist_meas == "a0_GUnifrac") dist.matrix = as.dist(guni.a0[ s, s ])
    else if (dist_meas == "a05_GUnifrac") dist.matrix = as.dist(guni.a05[ s, s ])
    else if (dist_meas == "Bray") dist.matrix = as.dist(bray[ s, s ])
    else if (dist_meas == "Jaccard") dist.matrix = as.dist(jaccard[ s, s ])
    else if (dist_meas == "Canberra") dist.matrix = as.dist(canberra[ s, s ])
    else if (dist_meas == "Aitchison") dist.matrix = as.dist(aitch[ s, s ])
    
    tax.to.clust <- ttc@otu_table
    
  }
  
  
  # Cluster using the Partitioning Around Medoids (PAM) algorithm
  pam.clustering = function(x,k) { # x is a distance matrix and k the number of clusters
    clust = as.vector(pam(as.dist(x), k, diss=TRUE))
    return(clust$clustering)
  }
  
  # Determine optimal number of clusters using the Calinski-Harabasz (CH) Index or silhouette
  # silhouette measures how similar one sample is to the others in the same cluster versus those in the neighbor cluster
  clus.meas <- matrix(NA, nrow = 10, ncol = 2)
  colnames(clus.meas) <- c("CH","obs")
  
  
  for (k in 1:10) { 
    if (k == 1) {
      clus.meas[k, "CH"] = NA
      clus.meas[k, "obs"] = NA
    } else {
      cluster_temp = pam.clustering(dist.matrix, k)
      if (CH.by.tagl==TRUE) {
        clus.meas[k, "CH"] = index.G1(t(tax.to.clust), cluster_temp,  d = dist.matrix, centrotypes = "medoids")
      } else {
        clus.meas[k, "CH"] = index.G1(t(noise.removal(glomTab[[ tl ]][ , sample_names(tagl) ], percent=0.01)), cluster_temp,  
                                      d = dist.matrix, centrotypes = "medoids")
      }
      clus.meas[k, "obs"] <- mean(silhouette(cluster_temp, dist.matrix)[,3])
    }
  }
  
  # incorporate prediction strength values as well
  clus.meas <- cbind(clus.meas, PS[[ dist_meas ]]$mean.pred)
  colnames(clus.meas) <- c("CH","obs","PS")
  clus.meas[1,"PS"] <- NA
  
  clus.meas.norm <- apply(clus.meas, 2, function(x) abs(x)/max(abs(x), na.rm = T))
  # clus.meas.norm <- apply(clus.meas, 2, function(x) x/max(x, na.rm = T))
  
  p.norm <- ggplot(reshape2::melt(clus.meas.norm), aes(x=Var1, y=value, fill=Var2)) +
    geom_bar(stat="identity", position = position_dodge(), width = 0.5) +
    scale_x_continuous(breaks = seq(0,10,1)) + theme(plot.title = element_text(hjust=0.5)) +
    ggtitle(sprintf('Optimal number of clusters by CH index, obs silhouette, Prediction Strength\ndist = %s', dist_meas)) +
    xlab("Number of clusters") + ylab('CH index / obs silhouette / Prediction Strength') + scale_fill_hue(name="score")
  
  p.orig <- ggplot(reshape2::melt(clus.meas[,c("obs","PS")]), aes(x=Var1, y=value, fill=Var2)) +
    geom_bar(stat="identity", position = position_dodge(), width = 0.5) +
    scale_x_continuous(breaks = seq(0,10,1)) + theme(plot.title = element_text(hjust=0.5)) +
    geom_hline(yintercept=0.5, linetype="dashed", color="indianred1") + 
    geom_hline(yintercept=0.8, linetype="dashed", color="turquoise3") + ylim(0,1) +
    annotate("text", 9, 0.8, vjust = -1, label = "Minimum suggested PS value", color="turquoise3") +
    annotate("text", 9, 0.5, vjust = -1, label = "Minimum suggested silhouette value", color="indianred1") +
    ggtitle(sprintf('Optimal number of clusters by obs silhouette, Prediction Strength\ndist = %s', dist_meas)) +
    xlab("Number of clusters") + ylab('obs silhouette / Prediction Strength') + scale_fill_hue(name="score")
  
  if (saveClusPlots) {
    ggsave(sprintf("%s/optimal.norm.%s.png", clusPlotsDir, dist_meas), plot = p.norm)
    ggsave(sprintf("%s/optimal.threshold.%s.png", clusPlotsDir, dist_meas), plot = p.orig)
  }
  
  
  # determine number of clusters by the max sum of normalized CH index and obs silhouette values
  nc <- sapply(1:nrow(clus.meas.norm), function(x) sum(clus.meas.norm[x,"CH"], clus.meas.norm[x,"obs"], clus.meas.norm[x,"PS"]))
  nclusters <- match(max(nc, na.rm = T), nc)
  
  # force use of indicated number of clusters if so desired
  if (! is.null(nClustForce)) {
    nclusters <- nClustForce
  }
  
  if (simulation) {
    return(nclusters)
    
  } else {
    # get cluster values
    cluster = pam.clustering( dist.matrix, k = nclusters )
    
    # must add values for those samples not included when appropriate
    cluster_full <- cluster
    cluster_full[ sample_names(tagl)[ ! sample_names(tagl) %in% s ] ] <- NA
    
    to.return <- list(clus_var_name, cluster_full, nclusters, tax.to.clust, dist.matrix, 
                      clus.meas, clus.meas.norm, p.norm, p.orig)
    names(to.return) <- c("clus_var_name","cluster_full","nclusters","tax.to.clust","dist.matrix",
                          "clus.meas","clus.meas.norm","p.norm","p.orig")
    return( to.return )
    # return( cluster_full )
  }
  
}
# ************************************************************************************************ #




# ******************************************************* #
# Get lists of taxa that are significantly greater in each Stomatotype of each type ####
# ******************************************************* #

get_sig_tax <- function(tagl, dist_meas, PS, glomTab, ncf=NULL, filt_anova=FALSE) {
  # list of tax ordered by overall abundance
  top.tax <- names(sort(taxa_sums(tagl), decreasing = T))
  
  # list of samples in each of the given Stomatotype
  clu <- get_clusters( tagl, tl, dist_meas = dist_meas, PS, glomTab, nClustForce = ncf, CH.by.tagl = FALSE ) 
  stoms <- sapply(unique(as.character( clu[[ "cluster_full" ]] )), function(x) 
    names(clu[[ "cluster_full" ]])[ clu[[ "cluster_full" ]] == x ])
  # stoms <- sapply(unique(as.character(as.matrix(SLL2@sam_data[,sprintf("Stomatotype_%s",dist_meas)]))), function(x) 
  #   rownames(SLL2@sam_data[SLL2@sam_data[,sprintf("Stomatotype_%s",dist_meas)] == x,]))
  
  # for each stomatotype, get a vector of taxa that are significantly greater in that stomatotype than in others
  st.dists <- sapply(names(stoms), function(st) {
    other_st <- names(stoms)[ st != names(stoms) ]
    
    sig <- sapply(top.tax, function(tax) {
      wilcox.test(as.numeric(as.matrix(otu_table(tagl)[ tax, stoms[[ st ]] ])), 
                  as.numeric(as.matrix(otu_table(tagl)[ tax, unlist(stoms[ other_st ]) ])),
                  alternative = "greater")$p.value
    })
    sig.adj <- p.adjust(sig, method = "fdr")
    names(sig.adj[sig.adj<0.05])
  })
  
  if (filt_anova==TRUE) return(filter_sig_tax_w_anova(st.dists, clu[[ "cluster_full" ]], tagl))
  else return(st.dists)
  
}
# ******************************************************* #
filter_sig_tax_w_anova <- function(st.dists, clusFull, tagl) {
  
  # In the case that there are more than 2 stoms, there are some taxa that appear 
  #    in multiple vectors in the list, since comparing each stom against all samples in all the different stoms
  #    - So now will try to use anova to determine which of those stoms, if any, has highest value of given taxon
  dups <- unname(unlist(st.dists))[ duplicated(unname(unlist(st.dists))) ]
  if (length(dups) > 0) {
    # create new list object keeping only those taxa that are not found in multiple stomatotypes
    new.std <- lapply(st.dists, function(x) x[ ! x %in% dups ])
    
    # here determine the Stomatotype in which each taxon is highest (or all Stoms if not one significantly higher)
    for (tax in unique(dups)) {
      # first get vector of Stomatotypes for which the given taxon is high
      dupIn <- sapply(names(st.dists), function(x) tax %in% st.dists[[ x ]])
      dupIn <- names(dupIn[dupIn])
      
      # run anova
      adf <- as.data.frame(cbind(as.matrix(clusFull), t(as.matrix(tagl@otu_table[tax,]))))
      colnames(adf) <- c("clust", tax)
      
      res.aov <- aov(formula = as.matrix(adf[, tax]) ~ as.factor(adf[, "clust"]), data = adf)
      thsd <- TukeyHSD(res.aov)$`as.factor(adf[, "clust"])`
      
      ### here determine which, if any, if the most significant Stomatotype for the given taxon
      if (sum(thsd[,"p adj"] < 0.05) == 0) {
        # if no signif diffs, just add tax to each Stomatotype that it already appears as high abundance
        print(tax)
        print(dupIn)
        print(thsd)
        for (st in dupIn) new.std[[ st ]] <- c(tax, new.std[[ st ]])
        
      } else {
        # then determine ultimate Stomatotype winner for given taxon
        #   initiate as first value of dupIn
        winner <- dupIn[1] 
        
        for (ro in rownames(thsd)) {
          s1 <- strsplit(ro,'-')[[1]][1]
          s2 <- strsplit(ro,'-')[[1]][2]
          
          if (s1 %in% dupIn & s2 %in% dupIn & thsd[ro,"p adj"] < 0.05) {
            # only check in the case that both are in dupIn
            
            if (thsd[ro,"diff"] > 0) higher <- s1
            else higher <- s2
            
            if (higher != winner) {
              if ("None" %in% c(higher, winner)) {
                # When filtering "consensus" stomatotypes, which includes a "None" category
                ro_check <- sprintf("None-%s", c(higher, winner)[c(higher, winner) != "None"])
                
                max.wi.hi <- "None"
                min.wi.hi <- c(higher, winner)[c(higher, winner) != "None"]
                
              } else {
                # Otherwise the table has rownames with higher number first
                max.wi.hi <- max( as.integer(winner), as.integer(higher) )
                min.wi.hi <- min( as.integer(winner), as.integer(higher) )
                
                ro_check <- sprintf("%s-%s", max.wi.hi, min.wi.hi)
              }
              
              # choose the better Stomatotype between the higher of the current row and the previously determined "winner"
              if (thsd[ ro_check, "diff"] > 0 & thsd[ ro_check, "p adj"] < 0.05) {
                winner <- as.character(max.wi.hi)
              } else if (thsd[ ro_check, "diff"] < 0 & thsd[ ro_check, "p adj"] < 0.05) {
                winner <- as.character(min.wi.hi)
              }
            }
          }
          
        }
        new.std[[ winner ]] <- c(tax, new.std[[ winner ]])
      }
    }
    
    # order the taxa in each Stomatotype by overall abundance
    top.tax <- names(sort(taxa_sums(tagl), decreasing = T))
    return( lapply(new.std, function(x) x[order(match(x, top.tax))]) )
    
  } else {
    return(st.dists)
  }
}
# ******************************************************* #






# ***************************************************************************************** #
# Check drivers and plot GRADIENTS ####

# ************************************* #
get_drivers <- function(tl, dist_or_var, glomTab, mTab, nTop, all_or_healthy="all", notStomatotype=T, strongest=T) {
  
  ah <- ifelse(all_or_healthy=="healthy","healthy.","")
  
  if (notStomatotype == T) {
    stom_col <- dist_or_var
  } else {
    stom_col <- sprintf("%sStomatotype_%s", ah, dist_or_var)
  }
  
  if (class(mTab[,stom_col]) == "factor" & sum(table(mTab[,stom_col])==0) > 0) {
    mTab[,stom_col] <- factor(mTab[,stom_col])
  }
  
  
  # ************************* #
  ## between-class analysis (BCA)
  k <- 10
  pca = dudi.pca(t(glomTab[[ tl ]][ , rownames(mTab) ]), scannf=F, nf=k)
  bet = bca(pca, fac=as.factor(mTab[ , stom_col ]), scannf=F, nf=k-1)
  # print(bet$tab)
  # ************************* #
  stomato.leaders <- list()
  
  for (sto in rownames(bet$tab)) {
    stomato.leaders[[ sto ]] <- sort(bet$tab[sto, ], decreasing = strongest)[1:nTop]
  }
  # for (i in 1:clu[[ "nclusters" ]]){
  #   stomato.leaders[i] <- colnames(bet$tab)[bet$tab[i,]==max(bet$tab[i,])]
  # }
  return(stomato.leaders)
  
}
# ************************************* #
library(gridExtra)
library(ggpubr)
library(grid)
# ************************************************************ #
plot_primary_drivers <- function(tl, mTab, glomTab, dist_or_var, nTop=1, all_or_healthy="all", 
                                 notStomatotype=T) {
  # ************************* #
  
  if (grepl("Diversity_group",dist_or_var))
    levs <- c("Low","Average","High")
  else
    levs <- sort(unique(mTab[, dist_or_var]))
  
  stomato.leaders <- get_drivers(tl, dist_or_var, glomTab, mTab, nTop = nTop,
                                 all_or_healthy = all_or_healthy, notStomatotype = notStomatotype)
  stomato.leaders <- stomato.leaders[levs]
  stomato.leaders <- lapply(stomato.leaders, names)
  
  # ************************* #
  # stomato.leaders <- c("Porphyromonas","Prevotella","Pseudomonas")
  stomato.melts <- list()
  # for (s in unlist(stomato.leaders)) {
  #   stomato.melts[[ s ]] <- cbind( reshape2::melt(t(clu[[ "tax.to.clust" ]][ s , ])), 
  #                                  reshape2::melt(as.matrix(clu[[ "cluster_full" ]])) )
  #   stomato.melts[[ s ]] <- stomato.melts[[ s ]][,c(1,2,3,6)]
  #   colnames(stomato.melts[[ s ]])[4] <- 'cluster_full'
  # }
  for (s in unlist(stomato.leaders)) {
    stomato.melts[[ s ]] <- cbind( reshape2::melt(t(glomTab[[ tl ]][s, rownames(mTab)])), 
                                   reshape2::melt(as.matrix(mTab[ , dist_or_var])) )
    stomato.melts[[ s ]] <- stomato.melts[[ s ]][,c(1,2,3,6)]
    colnames(stomato.melts[[ s ]])[4] <- 'cluster_full'
  }
  
  # ************************* #
  
  # to give the same scale in each boxplot
  # ymax <- max(clu[[ "tax.to.clust" ]][stomato.leaders,])
  ymax <- max(glomTab[[ tl ]][unlist(stomato.leaders),])
  
  # ************************* #
  stomato.box.leaders <- list()
  for (s in unlist(stomato.leaders)) {
    stomato.box.leaders[[ s ]] <- ggplot(stomato.melts[[ s ]], 
                                         aes(x=factor(cluster_full, levels=levs), y=value, 
                                             fill=factor(cluster_full, levels=levs))) +
      geom_boxplot(notch = T) +
      # geom_violin(scale = "count", trim = F, draw_quantiles = c(0.25, 0.5, 0.75)) +
      # geom_point() +
      ggtitle( s ) +# ylim(0, ymax) +
      xlab(NULL) + ylab('CLR') + theme(legend.position="none")
  }
  # ************************* #
  
  grid.arrange(grobs = stomato.box.leaders,
               top = textGrob(sprintf("Drivers at %s level for %s", tl, dist_or_var),
                              gp = gpar(fontsize=15)),
               # ncol = length( stomato.box.leaders ),
               nrow = ifelse(nTop > 1, length(levs), 1))
}
# ************************************************************ #

# ************************************************************ #
rank_drivers_for_subSamps <- function(variable, all_subs, mTab, glomTab, nTop, tl) {
  
  per_sub <- lapply(all_subs, function(i)
    lapply(get_drivers(tl, variable, glomTab, mTab[i,], nTop, notStomatotype = T), function(x) as.matrix(x)[1,])
  )
  
  sapply(names(per_sub$`1`), function(g) 
    sapply(rownames(glomTab[[ tl ]]), function(tax)
      mean(unlist(lapply(per_sub, function(i) i[[ g ]][ tax ] )))
    )
  )
  
}
# ************************************************************ #


# ****************************************************************************************************************** ####
# ****************************************************************************************************************** ####



# ****************************************************************************************************************** #
# PERMANOVA with adonis ####
# ****************************************************************************************************************** #


# ***************************************************************************************** #
# ***************************************************************************************** #
# get_adonis <- function(trait, covs, dist_meas, mTab, strataVar="School_ID.letter", customEffects=NULL) {
get_adonis <- function(trait, covs, dist_meas, mTab, strataVar="seqGroup", useStrata=T,
                       customEffects=NULL, printMods=F, ordObj=NULL) {
  
  # because adonis function wont work with the age bin variables that start with numbers:
  if (trait %in% c("13_20","20_30","30_40","40_50","50_60","60+")) {
    mTab[ , sprintf("AB_%s", gsub("\\+","_",trait)) ] <- mTab[ , trait ]
    trait <- sprintf("AB_%s", gsub("\\+","_",trait))
  }
  
  # ********************************** #
  b_dists <- c("JSD","Weighted_Unifrac","Unweighted_Unifrac",
               "VAW_GUnifrac","a0_GUnifrac","a05_GUnifrac",
               "Bray","Jaccard","Canberra", "Aitchison")
  bdist.codes <- structure(c("jsd","weighted_Unifrac","unweighted_Unifrac",
                             "guni.VAW","guni.a0","guni.a05",
                             "bray","jaccard","canberra", "aitch"), .Names = b_dists)
  # ********************************** #
  if ( ! is.null(ordObj) ) {
    dist_obj <- "as.matrix( ordObj )"
  } else {
    dist_obj <- bdist.codes[ dist_meas ]
  }
  
  # exclude trait from covs vector if already there to avoid duplicating
  if (trait %in% covs) covs <- covs[ covs != trait ]
  
  # get only those columns that are needed
  mTab <- mTab[ , c(trait, covs, strataVar) ]
  # change "No Sabe/No Contesta" to NA
  mTab[ mTab=="No Sabe/No Contesta"] <- NA
  # then keep only rows without NAs, or will get an error from adonis function
  mTab <- mTab[ sapply(rownames(mTab), function(x) sum(is.na(mTab[x, ]))==0), ]
  # for GM.index, only consider Cases, since control values are all far lower, cant really compare to cases
  if (trait == "GM.index") {
    mTab <- mTab[ mTab[,"Case.control"]=="Case", ]
    covs <- covs[ covs != "Case.control" ] # since only one value, cant use Case.control here
  }
  # sample names that are kept for this analysis
  sn <- rownames(mTab)
  
  # ********************************** #
  # each cov has to have at least 2 values...if not the case, exclude that cov
  covCounts <- sapply(colnames(mTab), function(x) length(unique(mTab[,x])))
  mTab <- mTab[ , names(covCounts[ covCounts > 1]) ]
  covs <- covs[ covs %in% colnames(mTab) ]
  
  # ********************************** #
  
  # if trait column is removed in previous step bc only 1 value, give non-signif values in table
  if ( ! trait %in% colnames(mTab)) {
    aov <- matrix(NA, nrow=length(covs)+3, ncol=6)
    rownames(aov) <- c(trait, covs, "Residuals","Total")
    colnames(aov) <- c("Df","SumOfSqs","MeanSqs","F.Model","R2","Pr(>F)")
    aov[c(trait,covs), "Pr(>F)"] <- rep(1, length(c(trait,covs)))
    return(aov)
  }
  
  # ********************************** #
  
  # get string of fixed effects
  if ( is.null(customEffects) )
    fixed_effs <- paste(c(trait, covs), collapse = " + ")
  else
    # if using a custom effects model (maybe includes interaction)
    fixed_effs <- customEffects
  # ********************************** #
  
  # adonis analysis
  if (trait == strataVar | ! strataVar %in% colnames(mTab))
    # dont use strata if that var is part of the fixed effects, or if the column was removed bc only 1 group used
    ado <- adonis(as.formula(sprintf('as.dist(%s[sn, sn]) ~ %s', dist_obj, fixed_effs)), 
                  data=mTab)
  else if (useStrata == T & ! is.null(strataVar))
    ado <- adonis(as.formula(sprintf('as.dist(%s[sn, sn]) ~ %s', dist_obj, fixed_effs)), 
                  data=mTab, 
                  strata = mTab[ , strataVar])
  
  # ********************************** #
  
  # check for covariates, and if signif, include an extra term in the formula for interaction
  aov <- as.data.frame(ado$aov.tab)
  sig.covs <- covs[ aov[covs, "Pr(>F)"] < 0.05 ]
  
  if (printMods) {
    print("first aov")
    print(aov)
    print(sig.covs)
  }
  
  
  # ********************************** #
  # RERUN MODEL after removing those fixed effects that are not
  fixeds.sigs <- paste(c(trait, sig.covs), collapse = " + ")
  # if all are significant, no need to rerun:
  if (length(sig.covs) == length(covs)) {
    aov.sigs <- aov
    sig.covs.RERUN <- sig.covs
    
  } else {
    if (trait == strataVar | ! strataVar %in% colnames(mTab))
      ado.sigs <- adonis(as.formula(sprintf('as.dist(%s[sn, sn]) ~ %s', dist_obj, fixeds.sigs)), 
                         data=mTab)
    else if (useStrata == T & ! is.null(strataVar))
      ado.sigs <- adonis(as.formula(sprintf('as.dist(%s[sn, sn]) ~ %s', dist_obj, fixeds.sigs)), 
                         data=mTab, strata = mTab[ , strataVar])
    
    aov.sigs <- as.data.frame(ado.sigs$aov.tab)
    sig.covs.RERUN <- sig.covs[ aov.sigs[sig.covs, "Pr(>F)"] < 0.05 ]
    
    if(printMods) {
      print("after removing non sigs")
      print(aov.sigs)
      print(sig.covs.RERUN)
    }
  }
  
  
  # will return this object if no covs were significant, which means will not enter into the next "if" section
  aov.to_return <- aov.sigs
  
  
  
  # ********************************** #
  # ********************************** #
  # rerun adonis with interactions included only if at least 1 covariate was also signif, 
  #   to check if interaction bt that cov/covs and trait is also significant
  if ( length(sig.covs.RERUN)>0 ) {
    
    fixeds.interact <- paste( paste(trait, sig.covs.RERUN, sep = ":"), collapse = " + " )
    
    if (trait == strataVar | ! strataVar %in% colnames(mTab))
      ado.interact <- adonis(as.formula(sprintf('as.dist(%s[sn, sn]) ~ %s + %s', 
                                                dist_obj, fixeds.sigs, fixeds.interact)), 
                             data=mTab)
    else if (useStrata == T & ! is.null(strataVar))
      ado.interact <- adonis(as.formula(sprintf('as.dist(%s[sn, sn]) ~ %s + %s', 
                                                dist_obj, fixeds.sigs, fixeds.interact)), 
                             data=mTab, strata = mTab[ , strataVar])
    
    aov.interact <- as.data.frame(ado.interact$aov.tab)
    sig_and_interact <- c(sig.covs.RERUN, paste(trait, sig.covs.RERUN, sep = ":"))
    
    sigs.interact <- sig_and_interact[ aov.interact[sig_and_interact, "Pr(>F)"] < 0.05 ]
    # since for some cases, presumably where there are multiple groups respresented by only 1 variable,
    #    such as the combo of Azythromycin and SampleStatus, or InhCalcineurin and SampleStatus,
    #    may end up with no row for the interaction term in the model, gives an NA here, must remove
    sigs.interact <- sigs.interact[ ! is.na(sigs.interact)]
    
    if (printMods) {
      print(ado.interact)
      print(sig_and_interact)
      
      print("after interact")
      print(aov.interact)
      print(sigs.interact)
    }
    
    
    # ********************************** #
    # AGAIN - RERUN MODEL after removing those fixed effects that are not
    fixeds.sigs.interact <- paste(c(trait, sigs.interact), collapse = " + ")
    
    if (printMods) {
      print(fixeds.sigs.interact)
    }
    
    # if all are significant, no need to rerun:
    if (length(sig.covs) == length(covs)) {
      aov.sigs.interact <- aov.interact
      sigs.interact.RERUN <- sigs.interact
      
    } else {
      if (trait == strataVar | ! strataVar %in% colnames(mTab))
        ado.sigs.interact <- adonis(as.formula(sprintf('as.dist(%s[sn, sn]) ~ %s', 
                                                       dist_obj, fixeds.sigs.interact)), 
                                    data=mTab)
      else if (useStrata == T & ! is.null(strataVar))
        ado.sigs.interact <- adonis(as.formula(sprintf('as.dist(%s[sn, sn]) ~ %s', 
                                                       dist_obj, fixeds.sigs.interact)), 
                                    data=mTab, strata = mTab[ , strataVar])
      
      aov.sigs.interact <- as.data.frame(ado.sigs.interact$aov.tab)
      sigs.interact.RERUN <- sigs.interact[ aov.sigs.interact[sigs.interact, "Pr(>F)"] < 0.05 ]
    }
    
    if (printMods) {
      print("after removing non sigs from interact")
      print(aov.sigs.interact)
      print(sigs.interact.RERUN)
    }
    
    aov.to_return <- aov.sigs.interact
    
  }
  # ********************************** #
  # ********************************** #
  
  # because adonis function wont work with the age bin variables that start with numbers:
  if (trait %in% c("AB_13_20","AB_20_30","AB_30_40","AB_40_50","AB_50_60","AB_60_")) {
    rownames(aov.to_return) <- gsub("60_", "60+", gsub("AB_","",rownames(aov.to_return)) )
  }
  
  # ********************************** #
  return( aov.to_return )
}
# ***************************************************************************************** #
# ***************************************************************************************** #





# ****************************************************************************************************************** #
# plot ordination with significant adonis variables ####
# *********************************************************** #

library(RColorBrewer)

plot_adonis_w_covars <- function(dist_meas, phy, metaTab, glomTab, shapeBy, colorBy, tl, 
                                 adoTab=NULL, adoCovList=NULL, anova.pTab=NULL, psByShape=F, 
                                 pcoaList=NULL, mce_tab=NULL, mce_centered="No",
                                 save.adoPlot=F, ellipses=T, log_of_color=F, swapShapeColor=F) {
  
  # *********************************************************************** #
  # *********************************************************************** #
  # *********************************************************************** #
  if ( ! is.null(mce_tab)) {
    # use MCE or ncMCE calculated coordinates instead of PCoA object to plot points
    if (colorBy %in% colnames(metaTab)) {
      plotVals <- cbind(mce_tab, metaTab[ , c(shapeBy, colorBy) ])
    } else if (tl %in% c("Phylum","Class","Order","Family","Genus","Species")) {
      plotVals <- cbind(mce_tab, metaTab[ , shapeBy ], as.numeric(glomTab[[ tl ]][ colorBy, ]))
    }
    
    xl <- "Dim1"
    yl <- "Dim2"
    # *********************************************************************** #
    # *********************************************************************** #
    # *********************************************************************** #
  } else {
    # otherwise, use the PCoA object as usual
    
    # ********************************* #
    # get PCoA object for the new distance matrix
    if ( ! is.null(pcoaList) ) {
      # use the precalculated pcoa coords if provided (useful if you want to the samples to stay in same coordinates)
      pcoaObj <- pcoaList[[ dist_meas ]]
      
      # ********************************* #
    } else {
      # otherwise caluclate the coordinates here using only non-NA samples for given trait
      
      # first keep only those samples for which colorBy is not NA
      if ( ! tl %in% c("Phylum","Class","Order","Family","Genus","Species")) {
        metaTab <- metaTab[ ! is.na(metaTab[,colorBy]), ]
        phy <- prune_samples(rownames(metaTab), phy)
      }
      
      # then get distance matrix using only samples without NAs for given variable
      if (dist_meas == "Weighted_Unifrac") dm <- UniFrac(phy, weighted = T, parallel = T)
      else if (dist_meas == "Unweighted_Unifrac") dm <- UniFrac(phy, weighted = F, parallel = T)
      else if (dist_meas == "VAW_GUnifrac") dm <- GUniFrac(t(phy@otu_table), phy@phy_tree)$unifracs[,,"d_VAW"]
      else if (dist_meas == "a05_GUnifrac") dm <- GUniFrac(t(phy@otu_table), phy@phy_tree)$unifracs[,,"d_0.5"]
      else if (dist_meas == "a0_GUnifrac") dm <- GUniFrac(t(phy@otu_table), phy@phy_tree)$unifracs[,,"d_0"]
      else if (dist_meas == "JSD") dm <- distance(phy, method = "jsd", parallel = T)
      else if (dist_meas == "Bray") dm <- vegdist(t(phy@otu_table), method = "bray")
      else if (dist_meas == "Jaccard") dm <- vegdist(decostand(as.data.frame(t(phy@otu_table)), method="pa"), method = "jaccard")
      else if (dist_meas == "Canberra") dm <- vegdist(t(phy@otu_table), method = "canberra")
      else if (dist_meas == "Aitchison") dm <- aDist( cmultRepl(t(phy@otu_table), method="CZM", label=0) )
      pcoaObj <- ape::pcoa(dm)
    }
    # ********************************* #
    
    # make matrix of plot values
    # coordTab <- as.data.frame(pcoaList[[ dist_meas ]]$vectors[ , 1:3 ])
    coordTab <- as.data.frame(pcoaObj$vectors[ , 1:3 ])
    
    if ("Rel_corr_eig" %in% names(pcoaObj$values)) {
      # varPers <- pcoaList[[ dist_meas ]]$values$Rel_corr_eig * 100
      varPers <- pcoaObj$values$Rel_corr_eig * 100
    } else if ("Relative_eig" %in% names(pcoaObj$values)) {
      # in case that this variable was not calculated (not totally sure why still, but for
      #   Jaccard it does not do this one) will just use relative eigenvalues without correction
      # varPers <- pcoaList[[ dist_meas ]]$values$Relative_eig * 100
      varPers <- pcoaObj$values$Relative_eig * 100
    } else {
      print("Rel_corr_eig and Relative_eig not calculated")
    }
    
    xl <- sprintf("Axis1 (%s%%)", round(varPers[1], 2))
    yl <- sprintf("Axis2 (%s%%)", round(varPers[2], 2))
    # ********************************* #
    
    if (colorBy %in% colnames(metaTab)) {
      plotVals <- cbind(coordTab,
                        metaTab[ rownames(coordTab), c(shapeBy, colorBy) ])
    } else if (tl %in% c("Phylum","Class","Order","Family","Genus","Species")) {
      plotVals <- cbind(coordTab,
                        metaTab[ rownames(coordTab), shapeBy ],
                        as.numeric(glomTab[[ tl ]][ colorBy, rownames(coordTab) ]))
    }
  }
  # *********************************************************************** #
  # *********************************************************************** #
  # *********************************************************************** #
  
  
  
  colnames(plotVals) <- c("Axis1","Axis2","Axis3", "shapeBy","colorBy")
  plotVals <- plotVals[ ! is.na(plotVals[,"shapeBy"]), ]
  plotVals <- plotVals[ ! is.na(plotVals[,"colorBy"]), ]
  
  plotVals$shapeBy <- as.factor(plotVals$shapeBy)
  
  if (log_of_color == T)
    plotVals$colorBy <- log(plotVals$colorBy + 1)
  
  # ********************************* #
  # for title 
  b_dists <- c("JSD","Weighted_Unifrac","Unweighted_Unifrac",
               "VAW_GUnifrac","a0_GUnifrac","a05_GUnifrac",
               "Bray","Jaccard","Canberra", "Aitchison")
  bdist.codes <- structure(c("jsd","weighted_Unifrac","unweighted_Unifrac",
                             "guni.VAW","guni.a0","guni.a05",
                             "bray","jaccard","canberra", "aitch"), .Names = b_dists)
  
  # ********************************* #
  # ********************************* #
  # pval of trait of interest
  #  depending on test of interest, will use one input or another
  # ********************************* #
  if ( ! is.null(adoTab) ) {
    if (colorBy == "Group") {
      ado.ps <- "\nAdonis pval = NA"
      ado.covs.ps <- ""
    } else if (tl %in% c("Phylum","Class","Order","Family","Genus","Species")) {
      ado.ps <- ""
      ado.covs.ps <- ""
    } else {
      # ado.pval <- round( as.data.frame(adoList[[ bdist.codes[dist_meas] ]]$aov.tab)[colorBy, "Pr(>F)"], 4)
      ado.pval <- adoTab[ dist_meas, colorBy]
      ado.ps <- sprintf("\nAdonis pval = %s", ado.pval)
      
      # pvals of covariates, if any signif
      ado.covs.pvals <- adoCovList[[ colorBy ]][ , dist_meas ]
      # take only those covariates which were significant
      acp <- ado.covs.pvals[ ! is.na(ado.covs.pvals) & ado.covs.pvals != "" ]
      # make string for title
      ado.covs.ps <- ifelse(length(acp) == 0,
                            "\nOther effects:   none",
                            sprintf("\nOther effects:   %s", 
                                    paste(paste(names(acp), acp, sep="="), collapse=" || ")))
    }
    ps.for.title <- sprintf("%s%s", ado.ps, ado.covs.ps)
    
    # ********************************* #
  } else if ( ! is.null(anova.pTab)) {
    # have to reverse the change of the characters that needed to be changed during lm tests so can check table
    if (colorBy == "Absconditabacteriales_(SR1)")
      colorBy.pTab <- "Absconditabacteriales_.SR1."
    else if ( ! colorBy %in% colnames(metaTab))
      # only change this for taxa variables
      colorBy.pTab <- gsub("-", "\\.", colorBy)
    else
      colorBy.pTab <- colorBy
    
    
    # pval of trait of interest
    if (colorBy.pTab %in% colnames(metaTab)) {
      # ano.P <- anova.pTab[ sprintf("Stomatotype_%s", dist_meas), colorBy.pTab ]
      ano.P <- anova.pTab[ shapeBy, colorBy.pTab ]
    } else {
      # ano.P <- anova.pTab[ colorBy.pTab, sprintf("Stomatotype_%s", dist_meas) ]
      ano.P <- anova.pTab[ colorBy.pTab, shapeBy ]
    }
    ano.P.title <- sprintf("\nAnova pval = %s", signif(as.numeric(ano.P), digits = 3))
    
    # pvals of covariates, if any signif (ignore seqGroup, as this was included just to account for its variation)
    if (colorBy.pTab %in% colnames(metaTab)) {
      covsToCheck <- colnames(anova.pTab)[ ! colnames(anova.pTab) %in% c(colorBy.pTab, "seqGroup") ]
      # ano.covs.P <- anova.pTab[ sprintf("Stomatotype_%s", dist_meas), covsToCheck ]
      ano.covs.P <- anova.pTab[ shapeBy, covsToCheck ]
    } else {
      # covsToCheck <- colnames(anova.pTab)[ ! colnames(anova.pTab) %in% c(sprintf("Stomatotype_%s", dist_meas), "seqGroup") ]
      covsToCheck <- colnames(anova.pTab)[ ! colnames(anova.pTab) %in% c(shapeBy, "seqGroup") ]
      ano.covs.P <- anova.pTab[ colorBy.pTab, covsToCheck ]
    }
    # in case only 1 cov, will lose the name, must ensure it keeps it
    names(ano.covs.P) <- covsToCheck
    
    # take only those covariates which were significant
    acp <- ano.covs.P[ ano.covs.P != "" ]
    # make string for title
    ano.covs.P.title <- ifelse(length(acp) == 0,
                               "\nOther effects:   none",
                               sprintf("\nOther effects:   %s", 
                                       paste(paste(names(acp), signif(as.numeric(acp), digits=3),
                                                   sep="="), collapse=" || ")))
    
    # finally, add these values to plot title
    ps.for.title <- sprintf("%s%s", ano.P.title, ano.covs.P.title)
    
    # ********************************* #
  } else {
    ps.for.title <- ""
  }
  # ********************************* #
  # ********************************* #
  
  if (swapShapeColor==T) {
    temp.shapeBy <- colorBy
    colorBy <- shapeBy
    shapeBy <- temp.shapeBy
    
    temp.pv.shapeBy <- plotVals$colorBy
    plotVals$colorBy <- plotVals$shapeBy
    plotVals$shapeBy <- temp.pv.shapeBy
  }
  
  # ********************************* #
  
  shape_vals <- c(16, 12, 8, 18, 2, 3, 4, 5)
  if (colorBy == shapeBy) shape_vals <- rep(16, 8)
  shape_vals <- shape_vals[ 1:length(levels(plotVals[,"shapeBy"])) ]
  # print((class(plotVals[,"colorBy"])))
  # ********************************* #
  
  # ********************************* #
  adoPlot <- ggplot(plotVals, aes(x=Axis1, y=Axis2, color=colorBy), col = colorBy) +
    geom_point(aes(col=colorBy, shape=shapeBy), size=4) +
    theme_minimal() +
    # geom_text(aes(x=Axis1, y=Axis2, label=labels)) +
    # scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(11,"RdYlBu")))(11)) +
    scale_shape_manual(values=shape_vals) +
    labs(shape=shapeBy, col=ifelse(log_of_color==T,paste0(colorBy,"_log"), colorBy)) +
    xlab(xl) + ylab(yl) +
    guides(shape=ifelse(colorBy==shapeBy | shapeBy=="Project",
                        "none", guide_legend(NULL))) + # remove shape legend if same as shapeVarName, or if Project
    theme(plot.title = element_text(size=17), legend.text = element_text(size=13),
          legend.title = element_text(size = 15, face = "bold")) +
    ggtitle(sprintf("%s || %s%s", dist_meas, colorBy, ps.for.title))
  # ********************************* #
  if (class(plotVals[,"colorBy"]) %in% c("numeric","integer"))
    adoPlot <- adoPlot + scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(11,"RdYlBu")))(11))
  # adoPlot <- adoPlot + scale_color_gradientn(colors = colorRampPalette((brewer.pal(11,"Blues")))(11))
  else if (ellipses == T) 
    adoPlot <- adoPlot + stat_ellipse()
  # ********************************* #
  
  # ********************************* #
  
  if (save.adoPlot == T) {
    
    # depending on which test was being used to print pvals
    if ( ! is.null(adoTab)) {
      testName <- "Adonis"
    } else if ( ! is.null(anova.pTab)) {
      testName <- "ANOVA"
    } else {
      testName <- "noTest"
    }
    
    # determine fname if shapeBy not Project
    shapeEnd <- ""
    ps.fname <- colorBy
    if (psByShape==F & shapeBy != "Project") {
      ps.fname <- colorBy
      shapeEnd <- sprintf("-shape.%s", shapeBy)
    }
    
    colorEnd <- ""
    if (psByShape==T & colorBy != "Project") {
      ps.fname <- shapeBy
      colorEnd <- sprintf("-color.%s", colorBy)
    }
    
    # get proper folder name based on which coordinates are being used
    if ( is.null(mce_tab)) {
      ggsave(sprintf("%s/figures/Adonis_signifs/%s-%s-PCoA-%s%s.png", 
                     p2_dir, testName, ps.fname, dist_meas, shapeEnd, colorEnd),
             width = 9.23, height = 7.04, device = "png", plot = adoPlot)
      
    } else {
      if (mce_centered == "No") {
        ggsave(sprintf("%s/figures/Adonis_signifs.ncMCE/%s-%s-ncMCE-%s%s.png", 
                       p2_dir, testName, ps.fname, dist_meas, shapeEnd, colorEnd),
               width = 9.23, height = 7.04, device = "png", plot = adoPlot)
      } else if (mce_centered == "Yes") {
        ggsave(sprintf("%s/figures/Adonis_signifs.MCE/%s-%s-MCE-%s%s.png", 
                       p2_dir, testName, ps.fname, dist_meas, shapeEnd, colorEnd),
               width = 9.23, height = 7.04, device = "png", plot = adoPlot)
      }
      
    }
    
    # ggsave(sprintf("%s/figures/%s/%s-%s-PCoA-%s%s.png", p2_dir, coords_dir, testName, ps.fname, dist_meas, shapeEnd, colorEnd),
    #        width = 9.23, height = 7.04, device = "png", plot = adoPlot)
    
  } else {
    adoPlot
  }
}
# *********************************************************** #







# ****************************************************************************************************************** #
# Run linear model ####


library(lme4)
library(car)
library(nnet)

# ************************************************************************ #
# ************************************************************************ #
get_lm <- function(fixeds, tl, mTab, glomTab, div_vars, noRemove=F, rerun.nonSig=T, silentAnova=F,
                   polyRegression=NULL, polyDegree="NULL", polyKnots="NULL", 
                   expModel=NULL, logModel=NULL, dv_toSkip=NULL, rescale_depVar=NULL) {
  
  # ******************************************** #
  if (tl %in% c("contVar","moreVars")) {
    
    pval.mat <- matrix(NA, nrow=length(div_vars), ncol = length(fixeds))
    rownames(pval.mat) <- div_vars
    colnames(pval.mat) <- fixeds
    # ******************************************** #
  } else {
    rownames(glomTab) <- gsub("-",".", gsub("\\(",".", gsub("\\)",".", gsub(" ", "___", rownames(glomTab)))))
    
    mTab <- cbind(mTab, t(glomTab))
    
    pval.mat <- matrix(NA, nrow=nrow(glomTab), ncol = length(fixeds))
    rownames(pval.mat) <- sort( rownames(glomTab) )
    colnames(pval.mat) <- fixeds
  }
  # ******************************************** #
  
  if ( ! is.null(polyRegression) ) {
    # update formula to run polynomial regression 
    # comes from this tutorial: https://rcompanion.org/rcompanion/e_03.html
    fixeds_for_poly <- fixeds
    fixeds_for_poly[ fixeds_for_poly == polyRegression ] <- sprintf("splines::bs(%s, knots = %s, degree = %s)", 
                                                                    polyRegression, polyKnots, polyDegree)
    fs.FULL <- paste(fixeds_for_poly, collapse = " + ")
    
  } else if ( ! is.null(expModel)) {
    fixeds_for_exp <- fixeds
    fixeds_for_exp[ fixeds_for_exp == expModel ] <- sprintf("exp(%s)", expModel)
    fs.FULL <- paste(fixeds_for_exp, collapse = " + ")
    
  } else if ( ! is.null(logModel)) {
    fixeds_for_log <- fixeds
    fixeds_for_log[ fixeds_for_log == logModel ] <- sprintf("log(%s)", logModel)
    fs.FULL <- paste(fixeds_for_log, collapse = " + ")
    
  } else {
    fs.FULL <- paste(fixeds, collapse = " + ")
  }
  
  
  # ******************************************** #
  for (dependVar in rownames(pval.mat)) {
    
    mTab.nonNA <- mTab[ , c(fixeds, dependVar)]
    mTab.nonNA <- mTab.nonNA[ sapply(rownames(mTab.nonNA), function(x) sum(is.na(mTab.nonNA[x, ]))==0), ]
    
    # ******************** #
    # if for some reason, the trait of interest has only 1 non-NA value here, just give 1 as pval and move on
    if (length(unique(mTab.nonNA[ , fixeds[1]])) == 1) {
      pval.mat[dependVar, fixeds] <- NA
      print(sprintf("**** %s only 1 value here ****", fixeds[1]))
      next
    }
    # ******************** #
    # if there are some dependVars that want to skip, give 
    if ( ! is.null( dv_toSkip ) & dependVar %in% dv_toSkip) {
      pval.mat[dependVar, fixeds] <- NA
      next
    }
    # ******************** #
    
    # in case there is fixed column with only 1 value, remove that fixed effect
    for (fi in fixeds) {
      if (length(unique(mTab.nonNA[ , fi])) == 1) {
        fixeds <- fixeds[ fixeds != fi ]
        
        if ( ! is.null(polyRegression) ) {
          # update formula to run polynomial regression 
          # comes from this tutorial: https://rcompanion.org/rcompanion/e_03.html
          fixeds_for_poly <- fixeds
          fixeds_for_poly[ fixeds_for_poly == polyRegression ] <- sprintf("splines::bs(%s, knots = %s, degree = %s)", 
                                                                          polyRegression, polyKnots, polyDegree)
          fs.FULL <- paste(fixeds_for_poly, collapse = " + ")
          
        } else if ( ! is.null(expModel)) {
          fixeds_for_exp <- fixeds
          fixeds_for_exp[ fixeds_for_exp == expModel ] <- sprintf("exp(%s)", expModel)
          fs.FULL <- paste(fixeds_for_exp, collapse = " + ")
          
        } else if ( ! is.null(logModel)) {
          fixeds_for_log <- fixeds
          fixeds_for_log[ fixeds_for_log == logModel ] <- sprintf("log(%s)", logModel)
          fs.FULL <- paste(fixeds_for_log, collapse = " + ")
          
        } else {
          fs.FULL <- paste(fixeds, collapse = " + ")
        }
        
      }
    }
    # skip this Family because it gives an error, and is not interesting anyway
    if (tl == "Family" & dependVar == "0319.6G20") next
    # 
    # if (length(table(mtab[,dependVar][!is.na(mtab[,"IL.8"])])) < 3) print(dependVar)
    
    # ************************ #
    # in case the values of the dependVar are very small and need to be rescaled, in order to avoid numeric issues
    #    due to floating point precision that can occur later when using the Anova() function
    # https://stats.stackexchange.com/questions/479377/error-reports-residual-sum-squares-is-0-in-anova-analysis-in-r
    if ( ! is.null(rescale_depVar) ) {
      mod_dependVar <- sprintf("I(%s * %s)", dependVar, rescale_depVar)
    } else {
      mod_dependVar <- dependVar
    }
    # ************************ #
    
    
    # ************************ #
    if (startsWith(dependVar,"Stomatotype_") | 
        dependVar %in% c("MALDI.Yeast_detected","MALDI.Mold_detected","Full_MALDI.Candida_albicans",
                         "Full_MALDI.Candida_dubliniensis","Full_MALDI.Candida_glabrata",
                         "Full_MALDI.Candida_parapsilosis","Full_MALDI.Candida_guillermondii",
                         "Full_MALDI.Candida_intermedia","Full_MALDI.Candida_krusei",
                         "Full_MALDI.Candida_lusitaniae","Full_MALDI.Cryptococcus_spp",
                         "Full_MALDI.Debaryomyces_hansenii","Full_MALDI.Rhodotorula_mucilaginosa","Allergy")) {
      
      # for these categorical variables, have to use nominal logistic regression, with multinom function
      
      mod <- try( multinom(as.formula(sprintf("%s ~ %s", mod_dependVar, fs.FULL)), 
                           data=mTab[ , c(fixeds, dependVar)], trace = F) )
      if (class(mod) == "try-error")
        next
      # name of anova test used depends on type of linear model: here uses likelihood-ratio chisquare
      anova.test.pName <- "Pr(>Chisq)"
      
    } else if (tl %in% c("Phylum","Class","Order","Family","Genus","Species")) {
      
      # for continuous variables, use standard linear model
      mod <- glm(as.formula(sprintf("%s ~ %s", mod_dependVar, fs.FULL)), 
                 data=mTab[ , c(fixeds, dependVar)])
      # name of anova test used depends on type of linear model: here uses F-test
      anova.test.pName <- "Pr(>Chisq)"
      
    } else {
      # for continuous variables, use standard linear model
      mod <- lm(as.formula(sprintf("%s ~ %s", mod_dependVar, fs.FULL)), 
                data=mTab[ , c(fixeds, dependVar)])
      # name of anova test used depends on type of linear model: here uses F-test
      anova.test.pName <- "Pr(>F)"
    }
    # ************************ #
    
    
    # ************************ #
    # ************************ #
    # ************************ #
    am <- try( as.data.frame(Anova(mod)), silent = silentAnova )
    # print(dependVar)
    # print(am)
    # if gives "Error ... residual sum of squares is 0 (within rounding error)", just give NAs
    if (class(am) == "try-error") {
      next
    } else {
      
      # ************************ #
      # ************************ #
      if (am[ 1, anova.test.pName] < 0.05 & rerun.nonSig==T) {
        # ************* #
        # RERUN MODEL if the trait of interest is significant after removing those fixed effects that are not
        f.sigs <- fixeds[ am[fixeds, anova.test.pName] < 0.05 ]
        
        if ( ! is.null(polyRegression) ) {
          # update formula to run polynomial regression 
          # comes from this tutorial: https://rcompanion.org/rcompanion/e_03.html
          f.sigs_for_poly <- f.sigs
          f.sigs_for_poly[ f.sigs_for_poly == polyRegression ] <- sprintf("splines::bs(%s, knots = %s, degree = %s)", 
                                                                          polyRegression, polyKnots, polyDegree)
          fs.FULL.sigs <- paste(f.sigs_for_poly, collapse = " + ")
          
        } else if ( ! is.null(expModel)) {
          f.sigs_for_exp <- f.sigs
          f.sigs_for_exp[ f.sigs_for_exp == expModel ] <- sprintf("exp(%s)", expModel)
          fs.FULL.sigs <- paste(f.sigs_for_exp, collapse = " + ")
          
        } else if ( ! is.null(logModel)) {
          f.sigs_for_log <- f.sigs
          f.sigs_for_log[ f.sigs_for_log == logModel ] <- sprintf("log(%s)", logModel)
          fs.FULL.sigs <- paste(f.sigs_for_log, collapse = " + ")
          
        } else {
          fs.FULL.sigs <- paste(f.sigs, collapse = " + ")
        }
        
        # print(sprintf("  %s", f.sigs))
        # print(sprintf("  %s", am[fixeds, anova.test.pName][am[fixeds, anova.test.pName]<0.05]))
        
        # ************* #
        if (startsWith(dependVar,"Stomatotype_") | 
            dependVar %in% c("MALDI.Yeast_detected","MALDI.Mold_detected","Full_MALDI.Candida_albicans",
                             "Full_MALDI.Candida_dubliniensis","Full_MALDI.Candida_glabrata",
                             "Full_MALDI.Candida_parapsilosis","Full_MALDI.Candida_guillermondii",
                             "Full_MALDI.Candida_intermedia","Full_MALDI.Candida_krusei",
                             "Full_MALDI.Candida_lusitaniae","Full_MALDI.Cryptococcus_spp",
                             "Full_MALDI.Debaryomyces_hansenii","Full_MALDI.Rhodotorula_mucilaginosa","Allergy")) {
          # for these categorical variables, have to use nominal logistic regression, with multinom function
          mod.sigs <- try( multinom(as.formula(sprintf("%s ~ %s", mod_dependVar, fs.FULL.sigs)), 
                                    data=mTab[ , c(f.sigs, dependVar)], trace = F) )
          if (class(mod.sigs) == "try-error")
            next
          # name of anova test used depends on type of linear model: here uses likelihood-ratio chisquare
          anova.test.pName <- "Pr(>Chisq)"
          
        } else if (tl %in% c("Phylum","Class","Order","Family","Genus","Species")) {
          # for continuous variables, use standard linear model
          mod.sigs <- glm(as.formula(sprintf("%s ~ %s", mod_dependVar, fs.FULL.sigs)), 
                          data=mTab[ , c(f.sigs, dependVar)])
          # name of anova test used depends on type of linear model: here uses F-test
          anova.test.pName <- "Pr(>Chisq)"
        } else {
          # for continuous variables, use standard linear model
          mod.sigs <- lm(as.formula(sprintf("%s ~ %s", mod_dependVar, fs.FULL.sigs)), 
                         data=mTab[ , c(f.sigs, dependVar)])
          # name of anova test used depends on type of linear model: here uses F-test
          anova.test.pName <- "Pr(>F)"
        }
        # ************* #
        am.sigs <- try( as.data.frame(Anova(mod.sigs)) )
        if (class(am.sigs) == "try-error") {
          next
        } else {
          if ( ! is.null(polyRegression) ) {
            # update formula to run polynomial regression 
            # comes from this tutorial: https://rcompanion.org/rcompanion/e_03.html
            f.sigs_for_poly <- f.sigs
            f.sigs_for_poly[ f.sigs_for_poly == polyRegression ] <- sprintf("splines::bs(%s, knots = %s, degree = %s)", 
                                                                            polyRegression, polyKnots, polyDegree)
            pval.mat[dependVar, f.sigs] <- signif( am.sigs[f.sigs_for_poly, anova.test.pName], 3 )
            
          } else if ( ! is.null(expModel)) {
            f.sigs_for_exp <- f.sigs
            f.sigs_for_exp[ f.sigs_for_exp == expModel ] <- sprintf("exp(%s)", expModel)
            pval.mat[dependVar, f.sigs] <- signif( am.sigs[f.sigs_for_exp, anova.test.pName], 3 )
            
          } else if ( ! is.null(logModel)) {
            f.sigs_for_log <- f.sigs
            f.sigs_for_log[ f.sigs_for_log == logModel ] <- sprintf("log(%s)", logModel)
            pval.mat[dependVar, f.sigs] <- signif( am.sigs[f.sigs_for_log, anova.test.pName], 3 )
            
          } else {
            pval.mat[dependVar, f.sigs] <- signif( am.sigs[f.sigs, anova.test.pName], 3 )
          }
        }
        # ************* #
        
      } else {
        # otherwise just write values to table, will get ignored later anyway
        if ( ! is.null(polyRegression) ) {
          # update formula to run polynomial regression 
          # comes from this tutorial: https://rcompanion.org/rcompanion/e_03.html
          fixeds_for_poly <- fixeds
          fixeds_for_poly[ fixeds_for_poly == polyRegression ] <- sprintf("splines::bs(%s, knots = %s, degree = %s)", 
                                                                          polyRegression, polyKnots, polyDegree)
          pval.mat[dependVar, fixeds] <- signif( am[fixeds_for_poly, anova.test.pName], 3 )
          
        } else if ( ! is.null(expModel)) {
          fixeds_for_exp <- fixeds
          fixeds_for_exp[ fixeds_for_exp == expModel ] <- sprintf("exp(%s)", expModel)
          pval.mat[dependVar, fixeds] <- signif( am[fixeds_for_exp, anova.test.pName], 3 )
          
        } else if ( ! is.null(logModel)) {
          fixeds_for_log <- fixeds
          fixeds_for_log[ fixeds_for_log == logModel ] <- sprintf("log(%s)", logModel)
          pval.mat[dependVar, fixeds] <- signif( am[fixeds_for_log, anova.test.pName], 3 )
          
        } else {
          pval.mat[dependVar, fixeds] <- signif( am[fixeds, anova.test.pName], 3 )
        }
        # ************* #
      }
      # print(pval.mat)
      
    }
    # ************************ #
    # ************************ #
    # ************************ #
  }
  
  # if tl is Species, had to alter white space to run lm(), now change names back to normal
  rownames(pval.mat) <- gsub("___", " ", rownames(pval.mat))
  
  # remove the rows for those dependVar that had all NAs
  dependVar_to_keep <- sapply(rownames(pval.mat), function(x) 
    sum(is.na(pval.mat[x, ])) < ncol(pval.mat))
  
  
  
  # adjust pvalues
  pval.mat <- apply(pval.mat, 2, p.adjust, method="fdr")
  
  # ************************ #
  # ************************ #
  # if want to keep all values from table, dont remove non-signif/NA rows
  if (noRemove == T) {
    return( as.matrix(pval.mat) )
  }
  # ************************ #
  # ************************ #
  
  
  # if only 1 fixed effect included, have to maintain format of pval.mat
  if (length(fixeds)==1) {
    pval.mat <- as.data.frame(pval.mat[ dependVar_to_keep, ])
    colnames(pval.mat) <- fixeds
  } else {
    pval.mat <- pval.mat[ dependVar_to_keep, ]
  }
  
  # print only significant values
  pval.mat[ pval.mat > 0.05 ] <- ""
  pval.mat[ is.na(pval.mat) ] <- ""
  
  # remove rows with dependVar values that had no significant effects
  if (length(fixeds)==1) {
    rpm <- rownames(pval.mat)[pval.mat[,1]!=""]
    pval.mat <- as.data.frame(pval.mat[ rowSums(pval.mat != "") > 0, ])
    rownames(pval.mat) <- rpm
    colnames(pval.mat) <- fixeds
  } else {
    # pval.mat <- as.data.frame(pval.mat)[ rowSums(pval.mat != "") > 0, ]
    pval.mat <- as.data.frame(pval.mat)[ pval.mat[ , 1] != "", ]
  }
  
  
  return( as.matrix(pval.mat) )
  
}
# ************************************************************************ #
# ************************************************************************ #












# ****************************************************************************************************************** ####
# ****************************************************************************************************************** ####





# ****************************************************************************************************************** #
# Control samples for comparison to disorders ####

# ************************************************ #
generate_control_samples <- function(disorder, phy, mTab) {
  
  if (disorder=="Celiac") {
    # prepare control samples
    conts <- sample_names(subset_samples(phy, Chronic_disorder=="No" & Celiac_family=="No"))
    # should only take controls from seqGroup 1 since all the disorder samples were from here anyway
    conts <- conts[ mTab[ conts, "seqGroup"] == "One" ]
    conts <- conts[ ! is.na(mTab[conts, "Age"])]
    # get proportionate numbers based on comunities
    conts.And <- conts[mTab[conts, "Community"]=="Andalucía"]
    conts.And <- sample(conts.And, length(conts.And)*(6/51))
    conts.Ara <- conts[mTab[conts, "Community"]=="Aragón"]
    conts.Ara <- sample(conts.Ara, length(conts.Ara)*(1/51))
    conts.Cat <- conts[mTab[conts, "Community"]=="Cataluña"]
    conts.Cat <- sample(conts.Cat, length(conts.Cat)*(35/51))
    conts.Mad <- conts[mTab[conts, "Community"]=="Comunidad de Madrid"]
    conts.Mad <- sample(conts.Mad, length(conts.Mad)*(3/51))
    conts.Val <- conts[mTab[conts, "Community"]=="Comunidad Valenciana"]
    conts.Val <- sample(conts.Val, length(conts.Val)*(1/51))
    conts.Gal <- conts[mTab[conts, "Community"]=="Galicia"]
    conts.Gal <- sample(conts.Gal, length(conts.Gal)*(1/51))
    conts.PV  <- conts[mTab[conts, "Community"]=="País Vasco"]
    conts.PV  <- sample(conts.PV, length(conts.PV)*(4/51))
    conts <- sort(c(conts.And, conts.Ara, conts.Cat, conts.Mad, conts.Val, conts.Gal, conts.PV))
    # get proportionate numbers based on ages
    conts <- sort(c(sample(conts[mTab[conts, "Age"]<19], size = 17),
                    sample(conts[19<=mTab[conts, "Age"] & mTab[conts, "Age"]<56], size = 27),
                    sample(conts[mTab[conts, "Age"]>=56], size = 7)))
    # ensure proportionate numbers based on gender (M~17 and F~34), will allow +- 3
    if (length(conts[mTab[conts, "Gender"]=="M"]) < 14 | 
        length(conts[mTab[conts, "Gender"]=="M"]) > 20) stop("Gender not balanced proportionately")
    
  } else if (disorder=="Cystic_fibrosis") {
    
    # how many are "Yes" for given variable, to know how many to match ("Two" for seqGroup)
    YesSamps <- rownames(mTab[ ! is.na(mTab[,"Age"]) & ! is.na(mTab[,disorder]) & mTab[,disorder] %in% c("Yes","Two"), ])
    
    # prepare controls 
    conts <- sample_names(subset_samples(phy, Chronic_disorder=="No" & CF_family=="No"))
    # should only take controls from seqGroup 1 since all the disorder samples were from here anyway
    conts <- conts[ mTab[ conts, "seqGroup"] == "One" ]
    # since only 1 cf sample has filtered tap water, will remove those from controls to avoid bias
    conts <- conts[ ! is.na(mTab[conts, "Age"]) & ! is.na(mTab[conts, "Water_type_home"]) &
                      ! mTab[conts, "Water_type_home"]=="Del Grifo (Filtrada)"]
    
    # get counts of the YesSamps from each community, will take that many 
    conts.communities <- sapply(sort(unique(mTab[ , "Community"])), function(x) length(YesSamps[mTab[YesSamps,"Community"]==x])) 
    
    
    # get proportionate numbers based on comunities
    conts.And <- conts[mTab[conts, "Community"]=="Andalucía"]
    conts.And <- sample(conts.And, length(conts.And)*(conts.communities["Andalucía"]/length(YesSamps)))
    conts.Ara <- conts[mTab[conts, "Community"]=="Aragón"]
    conts.Ara <- sample(conts.Ara, length(conts.Ara)*(conts.communities["Aragón"]/length(YesSamps)))
    conts.Can <- conts[mTab[conts, "Community"]=="Cantabria"]
    conts.Can <- sample(conts.Can, length(conts.Can)*(conts.communities["Cantabria"]/length(YesSamps)))
    conts.Cat <- conts[mTab[conts, "Community"]=="Cataluña"]
    conts.Cat <- sample(conts.Cat, length(conts.Cat)*(conts.communities["Cataluña"]/length(YesSamps)))
    conts.Mad <- conts[mTab[conts, "Community"]=="Comunidad de Madrid"]
    conts.Mad <- sample(conts.Mad, length(conts.Mad)*(conts.communities["Comunidad de Madrid"]/length(YesSamps)))
    conts.Gal <- conts[mTab[conts, "Community"]=="Galicia"]
    conts.Gal <- sample(conts.Gal, length(conts.Gal)*(conts.communities["Galicia"]/length(YesSamps)))
    conts.PV  <- conts[mTab[conts, "Community"]=="País Vasco"]
    conts.PV  <- sample(conts.PV, length(conts.PV)*(conts.communities["País Vasco"]/length(YesSamps)))
    conts <- sort(c(conts.And, conts.Ara, conts.Can, conts.Cat, conts.Mad, conts.Gal, conts.PV))
    
    # get proportionate numbers based on ages
    conts <- sort(c(sample(conts[mTab[conts, "Age"]<19], size = 10),
                    sample(conts[19<=mTab[conts, "Age"] & mTab[conts, "Age"]<56], size = 21)))
    # ensure proportionate numbers based on gender (M~17 and F~14), will allow +- 3
    if (length(conts[mTab[conts, "Gender"]=="M"]) < 14 | 
        length(conts[mTab[conts, "Gender"]=="M"]) > 20) stop("Gender not balanced proportionately")
    
  } else if (disorder=="Downs_Syndrome") {
    # prepare control samples
    conts <- sample_names(subset_samples(phy, Chronic_disorder=="No" & DS_family=="No"))
    # should only take controls from seqGroup 1 since all the disorder samples were from here anyway
    conts <- conts[ mTab[ conts, "seqGroup"] == "One" ]
    conts <- conts[ ! is.na(mTab[conts, "Age"])]
    # get proportionate numbers based on comunities
    conts.And <- conts[mTab[conts, "Community"]=="Andalucía"]
    conts.And <- sample(conts.And, length(conts.And)*(3/26))
    conts.Cat <- conts[mTab[conts, "Community"]=="Cataluña"]
    conts.Cat <- sample(conts.Cat, length(conts.Cat)*(13/26))
    conts.Gal <- conts[mTab[conts, "Community"]=="Galicia"]
    conts.Gal <- sample(conts.Gal, length(conts.Gal)*(1/26))
    conts.PV  <- conts[mTab[conts, "Community"]=="País Vasco"]
    conts.PV  <- sample(conts.PV, length(conts.PV)*(9/26))
    conts <- sort(c(conts.And, conts.Cat, conts.Gal, conts.PV))
    # get proportionate numbers based on ages
    conts <- sort(c(sample(conts[mTab[conts, "Age"]<19], size = 11),
                    sample(conts[19<=mTab[conts, "Age"] & mTab[conts, "Age"]<56], size = 15)))
    # ensure proportionate numbers based on gender (M~14 and F~12), will allow +- 2
    if (length(conts[mTab[conts, "Gender"]=="M"]) < 11 | 
        length(conts[mTab[conts, "Gender"]=="M"]) > 17) stop("Gender not balanced proportionately")
  }
  
  return(conts)
}

# ************************************************************************************************ #
# ************************************************************************************************ #
generate_control_samples.larger <- function(disorder, phy, mTab) {
  # differs from the original function in that, instead of ending with same number of controls
  #   as number of disorder.samples, will instead take max number of controls possible to get 
  #   same balance of location, age, and gender.
  #     - will use same proportion of samples from given locations as before, but instead of 
  #       fixed number of samples from each age group, and fixed range of numbers from each gender,
  #       will use proportions of age and gender that match disorder.samples
  
  # ******************************************************* #
  if (disorder=="Celiac") {
    # prepare control samples
    conts <- sample_names(subset_samples(phy, Chronic_disorder=="No" & Celiac_family=="No"))
    # should only take controls from seqGroup 1 since all the disorder samples were from here anyway
    conts <- conts[ mTab[ conts, "seqGroup"] == "One" ]
    conts <- conts[ ! is.na(mTab[conts, "Age"])]
    # get proportionate numbers based on comunities
    conts.And <- conts[mTab[conts, "Community"]=="Andalucía"]
    conts.And <- sample(conts.And, 1.25*length(conts.And)*(6/51))
    conts.Ara <- conts[mTab[conts, "Community"]=="Aragón"]
    conts.Ara <- sample(conts.Ara, 1.25*length(conts.Ara)*(1/51))
    conts.Cat <- conts[mTab[conts, "Community"]=="Cataluña"]
    conts.Cat <- sample(conts.Cat, 1.25*length(conts.Cat)*(35/51))
    conts.Mad <- conts[mTab[conts, "Community"]=="Comunidad de Madrid"]
    conts.Mad <- sample(conts.Mad, 1.25*length(conts.Mad)*(3/51))
    conts.Val <- conts[mTab[conts, "Community"]=="Comunidad Valenciana"]
    conts.Val <- sample(conts.Val, 1.25*length(conts.Val)*(1/51))
    conts.Gal <- conts[mTab[conts, "Community"]=="Galicia"]
    conts.Gal <- sample(conts.Gal, 1.25*length(conts.Gal)*(1/51))
    conts.PV  <- conts[mTab[conts, "Community"]=="País Vasco"]
    conts.PV  <- sample(conts.PV,  1.25*length(conts.PV)*(4/51))
    conts <- sort(c(conts.And, conts.Ara, conts.Cat, conts.Mad, conts.Val, conts.Gal, conts.PV))
    # get proportionate numbers based on ages
    conts.midAge <- conts[19<=mTab[conts, "Age"] & mTab[conts, "Age"]<56]
    # 51/27 is inverse of the proportion of CF samps that were in young range
    conts.new_total <- length(conts.midAge) * (51/27)
    
    if (length(conts[mTab[conts, "Age"]<19]) < floor(conts.new_total*(17/51))) 
      stop("Can't balance age proportionately with this subset")
    if (length(conts[mTab[conts, "Age"]>=56]) < floor(conts.new_total*(7/51)))
      stop("Can't balance age proportionately with this subset")
    
    # 17/51 is the proportion of CF samps that were young
    conts.young  <- sample(conts[mTab[conts, "Age"]<19], size = conts.new_total*(17/51))
    # 7/51 is the proportion of CF samps that were elderly
    conts.elder  <- sample(conts[mTab[conts, "Age"]>=56], size = conts.new_total*(7/51))
    conts <- c(conts.midAge, conts.young, conts.elder)
    # conts <- sort(c(sample(conts[mTab[conts, "Age"]<19], size = length(conts)*(17/51)),
    #                 sample(conts[19<=mTab[conts, "Age"] & mTab[conts, "Age"]<56], size = length(conts)*(27/51)),
    #                 sample(conts[mTab[conts, "Age"]>=56], size = length(conts)*(7/51))))
    # ensure proportionate numbers based on gender (M~17 and F~34), will allow +- 3
    if (length(conts[mTab[conts, "Gender"]=="M"]) < length(conts)*(14/51) | 
        length(conts[mTab[conts, "Gender"]=="M"]) > length(conts)*(20/51)) stop("Gender not balanced proportionately")
    # ******************************************************* #
    
  } else if (disorder=="Cystic_fibrosis") {
    # prepare controls 
    conts <- sample_names(subset_samples(phy, Chronic_disorder=="No" & CF_family=="No"))
    # should only take controls from seqGroup 1 since all the disorder samples were from here anyway
    conts <- conts[ mTab[ conts, "seqGroup"] == "One" ]
    # since only 1 cf sample has filtered tap water, will remove those from controls to avoid bias
    conts <- conts[ ! is.na(mTab[conts, "Age"]) & ! mTab[conts, "Water_type_home"]=="Del Grifo (Filtrada)"]
    # get proportionate numbers based on comunities
    conts.And <- conts[mTab[conts, "Community"]=="Andalucía"]
    conts.And <- sample(conts.And, 2*length(conts.And)*(4/31))
    conts.Ara <- conts[mTab[conts, "Community"]=="Aragón"]
    conts.Ara <- sample(conts.Ara, 2*length(conts.Ara)*(1/31))
    conts.Can <- conts[mTab[conts, "Community"]=="Cantabria"]
    conts.Can <- sample(conts.Can, 2*length(conts.Can)*(11/31))
    conts.Cat <- conts[mTab[conts, "Community"]=="Cataluña"]
    conts.Cat <- sample(conts.Cat, 2*length(conts.Cat)*(6/31))
    conts.Mad <- conts[mTab[conts, "Community"]=="Comunidad de Madrid"]
    conts.Mad <- sample(conts.Mad, 2*length(conts.Mad)*(6/31))
    conts.Gal <- conts[mTab[conts, "Community"]=="Galicia"]
    conts.Gal <- sample(conts.Gal, 2*length(conts.Gal)*(2/31))
    conts.PV  <- conts[mTab[conts, "Community"]=="País Vasco"]
    conts.PV  <- sample(conts.PV,  2*length(conts.PV)*(1/31))
    conts <- sort(c(conts.And, conts.Ara, conts.Can, conts.Cat, conts.Mad, conts.Gal, conts.PV))
    # get proportionate numbers based on ages
    conts.midAge <- conts[19<=mTab[conts, "Age"] & mTab[conts, "Age"]<56]
    # 31/21 is inverse of the proportion of CF samps that were in young range
    conts.new_total <- length(conts.midAge) * (31/21)
    
    if (length(conts[mTab[conts, "Age"]<19]) < floor(conts.new_total*(10/31))) 
      stop("Can't balance age proportionately with this subset")
    
    # 10/31 is the proportion of CF samps that were young
    conts.young  <- sample(conts[mTab[conts, "Age"]<19], size = conts.new_total*(10/31))
    conts <- c(conts.midAge, conts.young)
    # conts <- sort(c(sample(conts[mTab[conts, "Age"]<19], size = length(conts)*(10/31)),
    #                 sample(conts[19<=mTab[conts, "Age"] & mTab[conts, "Age"]<56], size = length(conts)*(21/31))))
    # ensure proportionate numbers based on gender (M~17 and F~14), will allow +- 2
    if (length(conts[mTab[conts, "Gender"]=="M"]) < length(conts)*(15/31) | 
        length(conts[mTab[conts, "Gender"]=="M"]) > length(conts)*(19/31)) stop("Gender not balanced proportionately")
    # ******************************************************* #
    
  } else if (disorder=="Downs_Syndrome") {
    # prepare control samples
    conts <- sample_names(subset_samples(phy, Chronic_disorder=="No" & DS_family=="No"))
    # should only take controls from seqGroup 1 since all the disorder samples were from here anyway
    conts <- conts[ mTab[ conts, "seqGroup"] == "One" ]
    conts <- conts[ ! is.na(mTab[conts, "Age"])]
    # get proportionate numbers based on comunities
    conts.And <- conts[mTab[conts, "Community"]=="Andalucía"]
    conts.And <- sample(conts.And, 2*length(conts.And)*(3/26))
    conts.Cat <- conts[mTab[conts, "Community"]=="Cataluña"]
    conts.Cat <- sample(conts.Cat, 2*length(conts.Cat)*(13/26))
    conts.Gal <- conts[mTab[conts, "Community"]=="Galicia"]
    conts.Gal <- sample(conts.Gal, 2*length(conts.Gal)*(1/26))
    conts.PV  <- conts[mTab[conts, "Community"]=="País Vasco"]
    conts.PV  <- sample(conts.PV,  2*length(conts.PV)*(9/26))
    conts <- sort(c(conts.And, conts.Cat, conts.Gal, conts.PV))
    # get proportionate numbers based on ages
    conts.midAge <- conts[19<=mTab[conts, "Age"] & mTab[conts, "Age"]<56]
    # 26/15 is inverse of the proportion of CF samps that were in young range
    conts.new_total <- length(conts.midAge) * (26/15)
    
    if (length(conts[mTab[conts, "Age"]<19]) < floor(conts.new_total*(11/26))) 
      stop("Can't balance age proportionately with this subset")
    
    # 11/26 is the proportion of CF samps that were young
    conts.young  <- sample(conts[mTab[conts, "Age"]<19], size = conts.new_total*(11/26))
    conts <- c(conts.midAge, conts.young)
    # conts <- sort(c(sample(conts[mTab[conts, "Age"]<19], size = length(conts)*(11/26)),
    #                 sample(conts[19<=mTab[conts, "Age"] & mTab[conts, "Age"]<56], size = length(conts)*(15/26))))
    # ensure proportionate numbers based on gender (M~14 and F~12), will allow +- 2
    if (length(conts[mTab[conts, "Gender"]=="M"]) < length(conts)*(12/26) | 
        length(conts[mTab[conts, "Gender"]=="M"]) > length(conts)*(16/26)) stop("Gender not balanced proportionately")
    # ******************************************************* #
    
  }
  
  return(conts)
}
# ************************************************************************************************ #
generate_control_samples.sameComs <- function(disorder, phy, mTab) {
  # differs from the original function in that, instead of ending with same number of controls
  #   as number of disorder.samples, will instead take max number of controls possible to get 
  #   same balance of location, age, and gender.
  #     - will use same proportion of samples from given locations as before, but instead of 
  #       fixed number of samples from each age group, and fixed range of numbers from each gender,
  #       will use proportions of age and gender that match disorder.samples
  
  # ******************************************************* #
  if (disorder=="Celiac") {
    # prepare control samples
    conts <- sample_names(subset_samples(phy, Chronic_disorder=="No" & Celiac_family=="No"))
    # should only take controls from seqGroup 1 since all the disorder samples were from here anyway
    conts <- conts[ mTab[ conts, "seqGroup"] == "One" ]
    conts <- conts[ ! is.na(mTab[conts, "Age"])]
    
    # get all samples from the communities in which Celiac samples were taken, unless only 1 or 2
    conts.takeFull <- conts[ mTab[conts, "Community"] %in% c("Andalucía","Cataluña","Comunidad de Madrid","País Vasco")]
    # in those cases with few samples in a given region, take a small sample of those
    conts.Ara <- sample(conts[mTab[conts, "Community"]=="Aragón"], 20)
    conts.Val <- sample(conts[mTab[conts, "Community"]=="Comunidad Valenciana"], 20)
    conts.Gal <- sample(conts[mTab[conts, "Community"]=="Galicia"], 20)
    
    conts <- sort(c(conts.takeFull, conts.Ara, conts.Val, conts.Gal))
    
    # will take a total of half of the samples from these locations
    conts.new_total <- length(conts) / 2
    
    
    if (length(conts[mTab[conts, "Age"]<19]) < floor(conts.new_total*(17/51))) 
      stop("Can't balance age proportionately with this subset")
    if (length(conts[mTab[conts, "Age"]>=56]) < floor(conts.new_total*(7/51)))
      stop("Can't balance age proportionately with this subset")
    
    # 27/51 is the proportion of Celiac samps that were midAge
    conts.midAge  <- sample(conts[19<=mTab[conts, "Age"] & mTab[conts, "Age"]<56], size = conts.new_total*(27/51))
    # 17/51 is the proportion of Celiac samps that were young
    conts.young  <- sample(conts[mTab[conts, "Age"]<19], size = conts.new_total*(17/51))
    # 7/51 is the proportion of Celiac samps that were elderly
    conts.elder  <- sample(conts[mTab[conts, "Age"]>=56], size = conts.new_total*(7/51))
    conts <- c(conts.midAge, conts.young, conts.elder)
    # conts <- sort(c(sample(conts[mTab[conts, "Age"]<19], size = length(conts)*(17/51)),
    #                 sample(conts[19<=mTab[conts, "Age"] & mTab[conts, "Age"]<56], size = length(conts)*(27/51)),
    #                 sample(conts[mTab[conts, "Age"]>=56], size = length(conts)*(7/51))))
    # ensure proportionate numbers based on gender (M~17 and F~34), will allow +- 3
    if (length(conts[mTab[conts, "Gender"]=="M"]) < length(conts)*(14/51) | 
        length(conts[mTab[conts, "Gender"]=="M"]) > length(conts)*(20/51)) stop("Gender not balanced proportionately")
    # ******************************************************* #
    
  } else if (disorder=="Cystic_fibrosis") {
    # prepare controls 
    conts <- sample_names(subset_samples(phy, Chronic_disorder=="No" & CF_family=="No"))
    # should only take controls from seqGroup 1 since all the disorder samples were from here anyway
    conts <- conts[ mTab[ conts, "seqGroup"] == "One" ]
    # since only 1 cf sample has filtered tap water, will remove those from controls to avoid bias
    conts <- conts[ ! is.na(mTab[conts, "Age"]) & ! mTab[conts, "Water_type_home"]=="Del Grifo (Filtrada)"]
    
    # get all samples from the communities in which Celiac samples were taken, unless only 1 or 2
    conts.takeFull <- conts[ mTab[conts, "Community"] %in% c("Andalucía","Cantabria","Cataluña","Comunidad de Madrid")]
    # in those cases with few samples in a given region, take a small sample of those
    conts.Ara <- sample(conts[mTab[conts, "Community"]=="Aragón"], 20)
    conts.Gal <- sample(conts[mTab[conts, "Community"]=="Galicia"], 20)
    conts.PV  <- sample(conts[mTab[conts, "Community"]=="País Vasco"], 20)
    
    conts <- sort(c(conts.takeFull, conts.Ara, conts.Gal, conts.PV))
    
    # will take a total of half of the samples from these locations
    conts.new_total <- length(conts) / 2
    
    
    if (length(conts[mTab[conts, "Age"]<19]) < floor(conts.new_total*(10/31))) 
      stop("Can't balance age proportionately with this subset")
    # midAge samples usually too few to get this proportion in this case, so will adjust the total size
    #   in this case, instead will take all the midAge samps that have been selected, then make new
    #   total the appropriate number from there
    if (length(conts[19<=mTab[conts, "Age"] & mTab[conts, "Age"]<56]) < floor(conts.new_total*(21/31))) {
      conts.new_total <- length(conts[19<=mTab[conts, "Age"] & mTab[conts, "Age"]<56]) * (31/21)
    }
    
    
    # 21/31 is the proportion of CF samps that were midAge
    conts.midAge  <- sample(conts[19<=mTab[conts, "Age"] & mTab[conts, "Age"]<56], size = conts.new_total*(21/31))
    # 10/31 is the proportion of CF samps that were young
    conts.young  <- sample(conts[mTab[conts, "Age"]<19], size = conts.new_total*(10/31))
    conts <- c(conts.midAge, conts.young)
    
    # conts <- sort(c(sample(conts[mTab[conts, "Age"]<19], size = length(conts)*(10/31)),
    #                 sample(conts[19<=mTab[conts, "Age"] & mTab[conts, "Age"]<56], size = length(conts)*(21/31))))
    # ensure proportionate numbers based on gender (M~17 and F~14), will allow +- 2
    if (length(conts[mTab[conts, "Gender"]=="M"]) < length(conts)*(15/31) | 
        length(conts[mTab[conts, "Gender"]=="M"]) > length(conts)*(19/31)) stop("Gender not balanced proportionately")
    # ******************************************************* #
    
  } else if (disorder=="Downs_Syndrome") {
    # prepare control samples
    conts <- sample_names(subset_samples(phy, Chronic_disorder=="No" & DS_family=="No"))
    # should only take controls from seqGroup 1 since all the disorder samples were from here anyway
    conts <- conts[ mTab[ conts, "seqGroup"] == "One" ]
    conts <- conts[ ! is.na(mTab[conts, "Age"])]
    
    # get all samples from the communities in which Celiac samples were taken, unless only 1 or 2
    conts.takeFull <- conts[ mTab[conts, "Community"] %in% c("Cataluña","País Vasco")]
    # in those cases with few samples in a given region, take a small sample of those
    conts.And <- sample(conts[mTab[conts, "Community"]=="Andalucía"], 20)
    conts.Gal <- sample(conts[mTab[conts, "Community"]=="Galicia"], 20)
    
    conts <- sort(c(conts.takeFull, conts.And, conts.Gal))
    
    # will take a total of half of the samples from these locations
    conts.new_total <- length(conts) / 2
    
    
    if (length(conts[mTab[conts, "Age"]<19]) < floor(conts.new_total*(11/26))) 
      stop("Can't balance age proportionately with this subset")
    # midAge samples usually too few to get this proportion in this case, so will adjust the total size
    #   in this case, instead will take all the midAge samps that have been selected, then make new
    #   total the appropriate number from there
    if (length(conts[19<=mTab[conts, "Age"] & mTab[conts, "Age"]<56]) < floor(conts.new_total*(15/26))) {
      conts.new_total <- length(conts[19<=mTab[conts, "Age"] & mTab[conts, "Age"]<56]) * (26/15)
    }
    
    
    # 15/26 is the proportion of DS samps that were midAge
    conts.midAge  <- sample(conts[19<=mTab[conts, "Age"] & mTab[conts, "Age"]<56], size = conts.new_total*(15/26))
    # 11/26 is the proportion of DS samps that were young
    conts.young  <- sample(conts[mTab[conts, "Age"]<19], size = conts.new_total*(11/26))
    conts <- c(conts.midAge, conts.young)
    # conts <- sort(c(sample(conts[mTab[conts, "Age"]<19], size = length(conts)*(11/26)),
    #                 sample(conts[19<=mTab[conts, "Age"] & mTab[conts, "Age"]<56], size = length(conts)*(15/26))))
    # ensure proportionate numbers based on gender (M~14 and F~12), will allow +- 2
    if (length(conts[mTab[conts, "Gender"]=="M"]) < length(conts)*(12/26) | 
        length(conts[mTab[conts, "Gender"]=="M"]) > length(conts)*(16/26)) stop("Gender not balanced proportionately")
    # ******************************************************* #
    
  }
  
  return(conts)
}
# ************************************************************************************************ #



# ************************************************ #
sub_to_fill_disorder_vectors <- function(disorder, disorder.samples, nsubs, statTest, comparison, tl, 
                                         rarefyPhyloseq=F, rarefyVegan=F, rarefyGUniFrac=F, moreControls=T) {
  
  dis.pTabsList <- list()
  dis.sampList <- list()
  
  # only need a list here when running correlations
  if (statTest == "Correlation") {
    dis.correlTab <- list()
  } else {
    dis.correlTab <- NULL
  }
  
  for (i in 1:nsubs) {
    # get appropriate subsampling of control samples for given disorder
    #  - sometimes it may be that the first set of subsampling in generate_control_samples() (by community)
    #    will not leave enough samples for a given age group in the second set of subsampling (by age)
    #    so have to use the error handling below to set control samples
    controls <- NULL
    attempt <- 0
    while( is.null(controls) ) {
      attempt <- attempt + 1
      try(
        if (moreControls) {
          controls <- generate_control_samples.larger(disorder)
        } else {
          controls <- generate_control_samples(disorder)
        },
        silent = TRUE
      )
    }
    if (attempt>1) print(sprintf("At i = %s, required %s attempts at subsampling", i,attempt))
    # controls <- generate_control_samples(disorder)
    
    # fill tables of p values taking comparable control samples multiple times for each disorder
    interest <- c(disorder.samples, controls)
    # keep track of controls used in each case, so can come back to check them if necessary
    dis.sampList[[ i ]] <- controls
    
    # make new SLL2 object if rarefying data
    if (rarefyPhyloseq==TRUE) {
      SLL2.dis <- prune_samples(interest, SLL2)
      SLL2.dis <- rarefy_even_depth(SLL2.dis, verbose = F)
      
      # include updated column for number of 16S gene counts per sample
      SLL2.dis@sam_data[,'Gene_counts'] <- colSums(SLL2.dis@otu_table)
      # include column for number of OTUs identified
      SLL2.dis@sam_data[,'Num_OTUs'] <- apply(SLL2.dis@otu_table, 2, function(x) sum( x != 0))
      
      # include updated column for Faith's Phylogenetic Diversity (alpha diversity that incorporates branch lengths in a tree)
      faith <- pd(t(as.data.frame(SLL2.dis@otu_table)), SLL2.dis@phy_tree, include.root = F)
      faith$PD[ is.na(faith$PD) ] <- 0 # for one sample that has only 1 species with non-0 count, which gives NA
      SLL2.dis@sam_data[,'Faiths.PD'] <- faith$PD
      SLL2.dis@sam_data[,'Species_Richness'] <- faith$SR
      
    } else if (rarefyVegan==TRUE) {
      SLL2.dis <- prune_samples(interest, SLL2)
      minReads <- min(colSums(SLL2.dis@otu_table))
      SLL2.dis.rar <- drarefy(SLL2.dis@otu_table, minReads)
      otu_table(SLL2.dis) <- otu_table(SLL2.dis.rar, taxa_are_rows = T)
      
      # include updated column for number of 16S gene counts per sample
      SLL2.dis@sam_data[,'Gene_counts'] <- colSums(SLL2.dis@otu_table)
      # include column for number of OTUs identified
      SLL2.dis@sam_data[,'Num_OTUs'] <- apply(SLL2.dis@otu_table, 2, function(x) sum( x != 0))
      
      # include updated column for Faith's Phylogenetic Diversity (alpha diversity that incorporates branch lengths in a tree)
      faith <- pd(t(as.data.frame(SLL2.dis@otu_table)), SLL2.dis@phy_tree, include.root = F)
      faith$PD[ is.na(faith$PD) ] <- 0 # for one sample that has only 1 species with non-0 count, which gives NA
      SLL2.dis@sam_data[,'Faiths.PD'] <- faith$PD
      SLL2.dis@sam_data[,'Species_Richness'] <- faith$SR
      
    } else if (rarefyGUniFrac==TRUE) {
      SLL2.dis <- prune_samples(interest, SLL2)
      otr <- Rarefy(t(SLL2.dis@otu_table))
      otu_table(SLL2.dis) <- otu_table(t(otr$otu.tab.rff), taxa_are_rows = T)
      
      # include updated column for number of 16S gene counts per sample
      SLL2.dis@sam_data[,'Gene_counts'] <- colSums(SLL2.dis@otu_table)
      # include column for number of OTUs identified
      SLL2.dis@sam_data[,'Num_OTUs'] <- apply(SLL2.dis@otu_table, 2, function(x) sum( x != 0))
      
      # include updated column for Faith's Phylogenetic Diversity (alpha diversity that incorporates branch lengths in a tree)
      faith <- pd(t(as.data.frame(SLL2.dis@otu_table)), SLL2.dis@phy_tree, include.root = F)
      faith$PD[ is.na(faith$PD) ] <- 0 # for one sample that has only 1 species with non-0 count, which gives NA
      SLL2.dis@sam_data[,'Faiths.PD'] <- faith$PD
      SLL2.dis@sam_data[,'Species_Richness'] <- faith$SR
      
    } else {
      SLL2.dis <- SLL2
    }
    
    # run given statistical test to get tables of p-values (and correlation coefficients for "Correlation") and
    # multiply table by nsubs as a multiple test correction for each p-value (since running subsamplings nsubs times)
    # ********************************* #
    if (statTest == "Correlation") {
      
      # will only use the controls samples here (NOT interest) in order to compare significant correlations
      #   found only within controls vs sig cors only within disorder samples
      if (comparison %in% c("taxa","questions")) {
        # if (comparison=="questions") {
        cors_and_ps <- fill_cor_tables(comparison, NA, controls, only_cont)
      } else {
        cors_and_ps <- fill_cor_tables(comparison, tl, controls, only_cont)
      }
      # then apply multiple testing correction to pvals
      cors_and_ps[[ "ps.adj" ]] <- apply(cors_and_ps[[ "ps" ]], 2, p.adjust, method='bonferroni')
      # must make diagonal of the p-vals table equal to 1, in order to ignore question with itself
      if (comparison != "TvsQ") diag( cors_and_ps[[ "ps.adj" ]] ) <- 1
      #For some questions that have all 0s (occurs when doing cities alone)
      cors_and_ps[[ "ps.adj" ]][ is.na(cors_and_ps[[ "ps.adj" ]]) ] <- 1
      cors_and_ps[[ "cor" ]][ is.na(cors_and_ps[[ "cor" ]]) ] <- 0
      
      dis.pTabsList[[ sprintf('p.%s', i) ]] <- cors_and_ps[[ "ps.adj" ]] # use only the adjusted pvals
      dis.correlTab[[ sprintf('p.%s', i) ]] <- cors_and_ps[[ "cor" ]]
      
      # ********************************* #  
    } else if (statTest == "Kruskal-Wallis") {
      
      rarefyData <- ifelse(TRUE %in% c(rarefyPhyloseq, rarefyVegan, rarefyGUniFrac), TRUE, FALSE)
      
      dis.pTabsList[[ sprintf('p.%s', i) ]] <- fill_kw_p_table(SLL2.dis, comparison, interest, 
                                                               disorder, only_cont, rarefyData) * nsubs
      
      # ********************************* #
    } else if (statTest == "Chi-squared") {
      
      if (comparison=="questions") {
        dis.pTabsList[[ sprintf('p.%s', i) ]] <- 
          t( fill_chi_table(comparison, NA, interest, disorder, keep_full = TRUE, gq_for_disorders = groupQs) * nsubs )
      } else {
        dis.pTabsList[[ sprintf('p.%s', i) ]] <- 
          t( fill_chi_table(comparison, tl, interest, disorder, keep_full = TRUE) * nsubs )
      }
      
      # ********************************* #
    } else if (statTest == "Ratios") {
      
      dcList <- list(disorder.samples, controls)
      names(dcList) <- c(disorder, "controls")
      
      dis.pTabsList[[ sprintf('p.%s', i) ]] <- get_taxa_ratio.diffs(comparison, dcList, "abunds")
      # ********************************* #
    }
    
    # in the case that all cont values were 0 (as may sometimes happen for rare species)
    #   kw test produces a p-value of NaN or NA => will change those to nsubs
    dis.pTabsList[[ sprintf('p.%s', i) ]][ is.nan(dis.pTabsList[[ sprintf('p.%s', i) ]]) ] <- 1
    dis.pTabsList[[ sprintf('p.%s', i) ]][ is.na(dis.pTabsList[[ sprintf('p.%s', i) ]]) ] <- 1
    # any p-values greater than 1 should be made as just 1
    # print(dis.pTabsList[[ sprintf('p.%s', i) ]][ "Div.Observed", ])
    dis.pTabsList[[ sprintf('p.%s', i) ]][ dis.pTabsList[[ sprintf('p.%s', i) ]] > 1 ] <- 1
    # print(dis.pTabsList[[ sprintf('p.%s', i) ]][ "Div.Observed", ])
  } 
  
  
  # ************ #
  # get proportion of significant p-values after correction for each KW test
  sigProps <- matrix(0, 
                     nrow = nrow(dis.pTabsList[[ sprintf('p.%s', i) ]]), 
                     ncol = ncol(dis.pTabsList[[ sprintf('p.%s', i) ]]))
  rownames(sigProps) <- rownames(dis.pTabsList[[ sprintf('p.%s', i) ]])
  colnames(sigProps) <- colnames(dis.pTabsList[[ sprintf('p.%s', i) ]])
  
  # also fill table with mean adjusted p-values for each KW test over the nsubs subsamplings
  meanPs <- sigProps
  
  for (ro in rownames(sigProps)) {
    for (co in colnames(sigProps)) {
      # get proportions
      sigProps[ ro, co ] <- 100 * sum(sapply(1:nsubs, function(i)
        dis.pTabsList[[ sprintf('p.%s', i) ]][ ro, co ] < 0.05 )) / nsubs
      # get mean p-values
      meanPs[ ro, co ] <- mean(sapply(1:nsubs, function(i) 
        dis.pTabsList[[ sprintf('p.%s', i) ]][ ro, co ] ))
    }
  }
  # ************ #
  return(list("pTabsList"=dis.pTabsList, "sigProps"=sigProps, "meanPs"=meanPs, 
              "sampList"=dis.sampList, "correlTab"=dis.correlTab))
}

# ************************************************ #

















healthySamps <- rownames(SLL2.meta[ SLL2.meta$Chronic_disorder == "No", ])
meta.healthy <- SLL2.meta[ healthySamps, ]




# ************************************************************************************************ #
generate_YesNo_subsamples <- function(variable, mTab) {
  
  # how many are "Yes" for given variable, to know how many to match ("Two" for seqGroup)
  YesSamps <- rownames(mTab[ ! is.na(mTab[,"Age"]) & ! is.na(mTab[,variable]) & mTab[,variable] %in% c("Yes","Two"), ])
  # have to use the "No" as main group here since there are many more "Yes" samples
  if (variable %in% c("Fluoride_toothpaste","Wash_hands_before_eat","Wash_hands_after_bathroom","Do_you_feel_well"))
    YesSamps <- rownames(mTab[ ! is.na(mTab[,"Age"]) & ! is.na(mTab[,variable]) & mTab[,variable]=="No", ])
  # get counts of the YesSamps from each community, will take that many 
  conts.communities <- sapply(sort(unique(mTab[ , "Community"])), function(x) length(YesSamps[mTab[YesSamps,"Community"]==x])) 
  
  # prepare control samples ("One" for seqGroup)
  conts <- rownames(mTab[ ! is.na(mTab[,variable]) & mTab[,variable] %in% c("No","One"), ])
  # have to use the "Yes" as conts group here since there are many more "Yes" samples
  if (variable %in% c("Fluoride_toothpaste","Wash_hands_before_eat","Wash_hands_after_bathroom","Do_you_feel_well"))
    conts <- rownames(mTab[ ! is.na(mTab[,variable]) & mTab[,variable]=="Yes", ])
  # must only use those samples that have a value for Age
  conts <- conts[ ! is.na(mTab[conts, "Age"])]
  # if getting controls for a disorder, must be samples without Chronic_disorder, 
  #    if not, still dont want to use Chronic_disorder samples in controls
  conts <- conts[ mTab[conts, "Chronic_disorder"] == "No" ]
  # print(conts.communities)
  
  # get representative number of samples from each community
  conts.And <- tryCatch(sample( conts[mTab[conts, "Community"]=="Andalucía"], conts.communities["Andalucía"]),
                        error = function(x) conts[mTab[conts, "Community"]=="Andalucía"])
  conts.Ara <- tryCatch(sample( conts[mTab[conts, "Community"]=="Aragón"], conts.communities["Aragón"]),
                        error = function(x) conts[mTab[conts, "Community"]=="Aragón"])
  conts.Can <- tryCatch(sample( conts[mTab[conts, "Community"]=="Cantabria"], conts.communities["Cantabria"]),
                        error = function(x) conts[mTab[conts, "Community"]=="Cantabria"])
  conts.Cat <- tryCatch(sample( conts[mTab[conts, "Community"]=="Cataluña"], conts.communities["Cataluña"]),
                        error = function(x) conts[mTab[conts, "Community"]=="Cataluña"])
  conts.Mad <- tryCatch(sample( conts[mTab[conts, "Community"]=="Comunidad de Madrid"], conts.communities["Comunidad de Madrid"]),
                        error = function(x) conts[mTab[conts, "Community"]=="Comunidad de Madrid"])
  conts.Val <- tryCatch(sample( conts[mTab[conts, "Community"]=="Comunidad Valenciana"], conts.communities["Comunidad Valenciana"]),
                        error = function(x) conts[mTab[conts, "Community"]=="Comunidad Valenciana"])
  conts.Gal <- tryCatch(sample( conts[mTab[conts, "Community"]=="Galicia"], conts.communities["Galicia"]),
                        error = function(x) conts[mTab[conts, "Community"]=="Galicia"])
  conts.Bal <- tryCatch(sample( conts[mTab[conts, "Community"]=="Islas Baleares"], conts.communities["Islas Baleares"]),
                        error = function(x) conts[mTab[conts, "Community"]=="Islas Baleares"])
  conts.Mur <- tryCatch(sample( conts[mTab[conts, "Community"]=="Murcia"], conts.communities["Murcia"]),
                        error = function(x) conts[mTab[conts, "Community"]=="Murcia"])
  conts.PV  <- tryCatch(sample( conts[mTab[conts, "Community"]=="País Vasco"], conts.communities["País Vasco"]),
                        error = function(x) conts[mTab[conts, "Community"]=="País Vasco"])
  conts <- sort(c(conts.And, conts.Ara, conts.Can, conts.Cat, conts.Mad, conts.Val, conts.Gal, conts.Bal, conts.Mur, conts.PV))
  
  
  # if (variable %in% c("Smoker","Pets")) {
  # conts.And <- sample( conts[mTab[conts, "Community"]=="Andalucía"], 24)
  # conts.Ara <- sample( conts[mTab[conts, "Community"]=="Aragón"], 15)
  # conts.Can <- sample( conts[mTab[conts, "Community"]=="Cantabria"], 1)
  # conts.Cat <- sample( conts[mTab[conts, "Community"]=="Cataluña"], 48)
  # conts.Mad <- sample( conts[mTab[conts, "Community"]=="Comunidad de Madrid"], 33)
  # conts.Val <- sample( conts[mTab[conts, "Community"]=="Comunidad Valenciana"], 34)
  # conts.Gal <- sample( conts[mTab[conts, "Community"]=="Galicia"], 21)
  # conts.Bal <- sample( conts[mTab[conts, "Community"]=="Islas Baleares"], 5)
  # conts.Mur <- sample( conts[mTab[conts, "Community"]=="Murcia"], 17)
  # conts.PV  <- sample( conts[mTab[conts, "Community"]=="País Vasco"], 8)
  # conts <- sort(c(conts.And, conts.Ara, conts.Can, conts.Cat, conts.Mad, conts.Val, conts.Gal, conts.Bal, conts.Mur, conts.PV))
  # get proportionate numbers based on ages
  var.age.props <- sapply(c("Child","Teen","Adult","Senior"), function(x) 
    nrow(mTab[YesSamps,][mTab[YesSamps, "Age_groups"]==x,])/length(YesSamps))
  # print(var.age.props)
  # print(sapply(c("Child","Teen","Adult","Senior"), function(x)
  #   nrow(mTab[conts,][mTab[conts, "Age_groups"]==x,])/length(conts)))
  
  
  
  # ranges of accepted percentages are somewhat arbitrary, but are estimated based on total meta.healthy Age_groups proportions:
  #   Child     Teen    Adult   Senior
  # 1.04712 81.00224 15.10845  3.29095
  
  # give these general ranges of variation for each age group (apv = age_proportion_variation)
  apv <- c("Child"=0.02, "Teen"=0.2, "Adult"=0.15, "Senior"=0.05)
  # in some cases have to be more lax in order to get samples, 
  #    since some communities may have had too many or too few of a particular age group
  if (variable %in% c("MALDI.Yeast_detected","Hair_in_mouth"))
    apv <- c("Child"=0.02, "Teen"=0.28, "Adult"=0.2, "Senior"=0.07)
  else if (variable %in% c("Cystic_fibrosis","Downs_Syndrome","Celiac","Hypothyroidism"))
    apv <- c("Child"=0.2, "Teen"=0.5, "Adult"=0.5, "Senior"=0.2)
  else if (variable %in% c("Hypertension"))
    apv <- c("Child"=0.2, "Teen"=0.5, "Adult"=0.3, "Senior"=0.4)
  else if (variable %in% c("Circulatory_issues","Kidney_issues"))
    apv <- c("Child"=0.05, "Teen"=0.35, "Adult"=0.2, "Senior"=0.2)
  else if (variable %in% c("Anemia"))
    apv <- c("Child"=0.05, "Teen"=0.3, "Adult"=0.3, "Senior"=0.06)
  
  
  if (length(conts[mTab[conts, "Age_groups"]=="Child"]) > (var.age.props["Child"]+apv["Child"])*length(YesSamps))
    stop("Can't balance age proportionately with this subset")
  if (length(conts[mTab[conts, "Age_groups"]=="Teen"]) < (var.age.props["Teen"]-apv["Teen"])*length(YesSamps) |
      length(conts[mTab[conts, "Age_groups"]=="Teen"]) > (var.age.props["Teen"]+apv["Teen"])*length(YesSamps))
    stop("Can't balance age proportionately with this subset")
  if (length(conts[mTab[conts, "Age_groups"]=="Adult"]) < (var.age.props["Adult"]-apv["Adult"])*length(YesSamps) |
      length(conts[mTab[conts, "Age_groups"]=="Adult"]) > (var.age.props["Adult"]+apv["Adult"])*length(YesSamps))
    stop("Can't balance age proportionately with this subset")
  if (length(conts[mTab[conts, "Age_groups"]=="Senior"]) < (var.age.props["Senior"]-apv["Senior"])*length(YesSamps) |
      length(conts[mTab[conts, "Age_groups"]=="Senior"]) > (var.age.props["Senior"]+apv["Senior"])*length(YesSamps))
    stop("Can't balance age proportionately with this subset")
  
  
  # ensure proportionate numbers based on gender (M~92 and F~114), will allow +- 3
  var.gender.props <- sapply(c("F","M"), function(x)
    nrow(mTab[YesSamps,][mTab[YesSamps, "Gender"]==x,])/length(YesSamps))
  
  if (length(conts[mTab[conts, "Gender"]=="M"]) < (var.gender.props["M"]-0.1)*length(YesSamps) |
      length(conts[mTab[conts, "Gender"]=="M"]) > (var.gender.props["M"]+0.1)*length(YesSamps))
    stop("Gender not balanced proportionately")
  # ******************************************************* #
  # }
  
  # ******************************************************* #
  
  return(conts)
}
# ************************************************************************************************ #
# ************************************************************************************************ #
generate_YesNo_subsamples.additional <- function(variable, mTab) {
  
  
  # "Gender"                    "Mouth_wounds.binary"       "Analgesics"                "Bite_nails"                "Chew_pens"                
  # [26] "Kissing_partner"           "seqGroup"
  
  # First get the "Yes/F/Two" samps
  YesSamps <- rownames(mTab[ ! is.na(mTab[,"Age"]) & ! is.na(mTab[,variable]) & mTab[,variable] %in% c("Yes","F","Two"), ])
  # take a sample of 100 of those
  YesSamps <- sample( YesSamps, 100 )
  
  # get counts of the YesSamps from each community, will take that many 
  counts.communities <- sapply(sort(unique(mTab[ , "Community"])), function(x) length(YesSamps[mTab[YesSamps,"Community"]==x])) 
  # print(counts.communities)
  
  # prepare the "No/M/One" samps
  NoSamps <- rownames(mTab[ ! is.na(mTab[,"Age"]) & ! is.na(mTab[,variable]) & mTab[,variable] %in% c("No","M","One"), ])
  if (variable %in% c("Kissing_partner","Pets.Cats")) {
    NoSamps <- NoSamps[ mTab[ NoSamps, "Age_groups"] != "Child" ]
  } else if (variable %in% c("seqGroup")) {
    NoSamps <- NoSamps[ ! mTab[ NoSamps, "Age_groups"] %in% c("Child","Senior") ]
  }
  # take a sample of 100 of those
  NoSamps <- sample( NoSamps, 100 )
  
  acceptDiffs <- c("Andalucía"=10, "Aragón"=8, "Cantabria"=5, "Cataluña"=10, "Comunidad de Madrid"=10,
                   "Comunidad Valenciana"=10, "Galicia"=10,"Islas Baleares"=7, "Murcia"=8, "País Vasco"=8)
  # print(acceptDiffs)
  for (com in names(counts.communities)) {
    
    # dont use samples if there is a community with 0 in one group and the other has some from that community
    # if (counts.communities[com] == 0 & length(NoSamps[mTab[NoSamps, "Community"]==com]) > 0) {
    #   print(com)
    #   stop("Can't balance communities proportionately with this subset")
    # }
    
    if (length(NoSamps[mTab[NoSamps, "Community"]==com]) < counts.communities[com]-5 |
        length(NoSamps[mTab[NoSamps, "Community"]==com]) > counts.communities[com]+5) {
    # if (length(NoSamps[mTab[NoSamps, "Community"]==com]) < counts.communities[com] - acceptDiffs[com] |
    #     length(NoSamps[mTab[NoSamps, "Community"]==com]) > counts.communities[com] + acceptDiffs[com]) {
      # print(com)
      # print(c(length(NoSamps[mTab[NoSamps, "Community"]==com]), counts.communities[com]))
      stop("Can't balance communities proportionately with this subset")
    }
  }
  # get representative number of samples from each community
  
  
  # get proportionate numbers based on ages
  var.age.props <- sapply(c("Child","Teen","Adult","Senior"), function(x) 
    nrow(mTab[YesSamps,][mTab[YesSamps, "Age_groups"]==x,])/length(YesSamps))
  # print(var.age.props)
  # print(sapply(c("Child","Teen","Adult","Senior"), function(x)
  #   nrow(mTab[NoSamps,][mTab[NoSamps, "Age_groups"]==x,])/length(NoSamps)))
  
  
  
  # ranges of accepted percentages are somewhat arbitrary, but are estimated based on total meta.healthy Age_groups proportions:
  #   Child     Teen    Adult   Senior
  # 1.04712 81.00224 15.10845  3.29095
  
  # give these general ranges of variation for each age group (apv = age_proportion_variation)
  apv <- c("Child"=0.02, "Teen"=0.2, "Adult"=0.15, "Senior"=0.05)
  # print(var.age.props)
  # print(sapply(c("Child","Teen","Adult","Senior"), function(x) 
  #   nrow(mTab[NoSamps,][mTab[NoSamps, "Age_groups"]==x,])/length(NoSamps)))
  # in some cases have to be more lax in order to get samples, 
  #    since some communities may have had too many or too few of a particular age group
  if (variable %in% c("MALDI.Yeast_detected","Hair_in_mouth"))
    apv <- c("Child"=0.02, "Teen"=0.28, "Adult"=0.2, "Senior"=0.07)
  else if (variable %in% c("Cystic_fibrosis","Downs_Syndrome","Celiac","Hypothyroidism"))
    apv <- c("Child"=0.2, "Teen"=0.5, "Adult"=0.5, "Senior"=0.2)
  else if (variable %in% c("Hypertension"))
    apv <- c("Child"=0.2, "Teen"=0.5, "Adult"=0.3, "Senior"=0.4)
  else if (variable %in% c("Circulatory_issues","Kidney_issues"))
    apv <- c("Child"=0.05, "Teen"=0.35, "Adult"=0.2, "Senior"=0.2)
  else if (variable %in% c("Anemia"))
    apv <- c("Child"=0.05, "Teen"=0.3, "Adult"=0.3, "Senior"=0.06)
  else if (variable %in% c("Gender"))
    apv <- c("Child"=0.05, "Teen"=0.3, "Adult"=0.3, "Senior"=0.15)
  else if (variable %in% c("Mouth_wounds.binary"))
    apv <- c("Child"=0.05, "Teen"=0.3, "Adult"=0.2, "Senior"=0.2)
  else if (variable %in% c("Bite_nails","Chew_pens"))
    apv <- c("Child"=0.05, "Teen"=0.3, "Adult"=0.2, "Senior"=0.25)
  else if (variable %in% c("Kissing_partner"))
    apv <- c("Child"=1, "Teen"=0.3, "Adult"=0.4, "Senior"=0.15)
  # print(apv)
  
  if (length(NoSamps[mTab[NoSamps, "Age_groups"]=="Child"]) > (var.age.props["Child"]+apv["Child"])*length(YesSamps)) {
    # print("Child")
    stop("Can't balance age proportionately with this subset")
  } else if (length(NoSamps[mTab[NoSamps, "Age_groups"]=="Teen"]) < (var.age.props["Teen"]-apv["Teen"])*length(YesSamps) |
      length(NoSamps[mTab[NoSamps, "Age_groups"]=="Teen"]) > (var.age.props["Teen"]+apv["Teen"])*length(YesSamps)) {
    # print("Teen")
    stop("Can't balance age proportionately with this subset")
  } else if (length(NoSamps[mTab[NoSamps, "Age_groups"]=="Adult"]) < (var.age.props["Adult"]-apv["Adult"])*length(YesSamps) |
      length(NoSamps[mTab[NoSamps, "Age_groups"]=="Adult"]) > (var.age.props["Adult"]+apv["Adult"])*length(YesSamps)) {
    # print("Adult")
    stop("Can't balance age proportionately with this subset")
  } else if (length(NoSamps[mTab[NoSamps, "Age_groups"]=="Senior"]) < (var.age.props["Senior"]-apv["Senior"])*length(YesSamps) |
      length(NoSamps[mTab[NoSamps, "Age_groups"]=="Senior"]) > (var.age.props["Senior"]+apv["Senior"])*length(YesSamps)) {
    # print("Senior")
    stop("Can't balance age proportionately with this subset")
  }
  
  
  if (variable != "Gender") {
    # ensure proportionate numbers based on gender (M~92 and F~114), will allow +- 3
    var.gender.props <- sapply(c("F","M"), function(x)
      nrow(mTab[YesSamps,][mTab[YesSamps, "Gender"]==x,])/length(YesSamps))
    
    if (length(NoSamps[mTab[NoSamps, "Gender"]=="M"]) < (var.gender.props["M"]-0.1)*length(YesSamps) |
        length(NoSamps[mTab[NoSamps, "Gender"]=="M"]) > (var.gender.props["M"]+0.1)*length(YesSamps)) {
      # print("Gender")
      stop("Gender not balanced proportionately")
    }
    # ******************************************************* #
    # }
  }
  
  
  # ******************************************************* #
  
  return(list("YFT"=YesSamps, "NMO"=NoSamps))
}
# ************************************************************************************************ #





sG <- list()
for ( i in as.character(1:100)) {
  sG[[ i ]] <- get_samps_for_test("Smoker", meta.healthy, i, YesNo = T)
}

table(meta.healthy$Smoker)
length(unique(unlist(lapply(sG, function(x) x))))
unique(unlist(lapply(sG, function(x) length(x))))
unique(unlist(lapply(sG, function(x) length(unique(x)))))


sG.2 <- list()
for ( i in as.character(1:100)) {
  sG.2[[ i ]] <- get_samps_for_test("Smoker", meta.healthy, i, YesNo = T)
}

table(meta.healthy$Smoker)
table(meta.healthy$Smoker, meta.healthy$Age_groups)
table(meta.healthy$Smoker, meta.healthy$Gender)
length(unique(unlist(lapply(sG.2, function(x) x))))
unique(unlist(lapply(sG.2, function(x) length(x))))
unique(unlist(lapply(sG.2, function(x) length(unique(x)))))



sG.3 <- list()
for ( i in as.character(1:100)) {
  sG.3[[ i ]] <- get_samps_for_test("Smoker", meta.healthy, i, YesNo = T)
}

table(meta.healthy$Smoker, meta.healthy$Community)
table(meta.healthy[sG.3$`4`, "Community"])



# *************** #

gG <- list()
for ( i in as.character(1:100)) {
  print(i)
  gG[[ i ]] <- get_samps_for_test("Gender", meta.healthy, i, YesNo.add = T)
}

table(meta.healthy$Gender)
length(unique(unlist(lapply(gG, function(x) x$YFT))))
length(unique(unlist(lapply(gG, function(x) x$NMO))))
unique(unlist(lapply(gG, function(x) c(length(x$YFT), length(x$NMO)))))
unique(unlist(lapply(gG, function(x) c(length(unique(x$YFT)), length(unique(x$NMO))))))

# *************** #





# ************************************************************************************************ #
generate_ageGroup_subsamples <- function(whichSamps, mTab) {
  
  # ******************************************************* #
  if (whichSamps=="healthy") {
    
    # healthy Child and Senior (12+42=54)
    
    # prepare teen samples
    teens <- rownames(mTab[ ! is.na(mTab[,"Age_groups"]) & mTab[,"Age_groups"]=="Teen", ])
    # get proportionate numbers based on comunities
    teens.And <- sample( teens[mTab[teens, "Community"]=="Andalucía"], 12) # (2+10 Child and Senior)
    # no healthySamps were teens from Cantabria, so will skip, though 1 Child and 4 Seniors
    # teens.Can <- teens[mTab[teens, "Community"]=="Cantabria"]
    # teens.Can <- sample(teens.Can, length(teens.Can)*(5/nSubSamps)) # (1+4 Child and Senior)
    teens.Cat <- sample( teens[mTab[teens, "Community"]=="Cataluña"], 13) # (4+9 Child and Senior)
    teens.Mad <- sample( teens[mTab[teens, "Community"]=="Comunidad de Madrid"], 11) # (1+10 Child and Senior)
    teens.Val <- sample( teens[mTab[teens, "Community"]=="Comunidad Valenciana"], 5) # (1+4 Child and Senior)
    teens.Gal <- sample( teens[mTab[teens, "Community"]=="Galicia"], 2) # (1+1 Child and Senior)
    teens.Mur <- sample( teens[mTab[teens, "Community"]=="Murcia"], 1) # (0+1 Child and Senior)
    teens.PV  <- sample( teens[mTab[teens, "Community"]=="País Vasco"], 5) # (2+3 Child and Senior)
    teens <- sort(c(teens.And, teens.Cat, teens.Mad, teens.Val, teens.Gal, teens.Mur, teens.PV))
    
    # prepare adult samples
    adults <- rownames(mTab[ ! is.na(mTab[,"Age_groups"]) & mTab[,"Age_groups"]=="Adult", ])
    # get proportionate numbers based on comunities
    adults.And <- sample( adults[mTab[adults, "Community"]=="Andalucía"], 12) # (2+10 Child and Senior)
    # no healthySamps were teens from Cantabria, so will skip, though 1 Child and 4 Seniors
    #   but there were 4 adults from here, so cannot reach total of Child+Senior, instead will take 2 adults
    adults.Can <- sample( adults[mTab[adults, "Community"]=="Cantabria"], 2) # (1+4 Child and Senior)
    adults.Cat <- sample( adults[mTab[adults, "Community"]=="Cataluña"], 13) # (4+9 Child and Senior)
    adults.Mad <- sample( adults[mTab[adults, "Community"]=="Comunidad de Madrid"], 11) # (1+10 Child and Senior)
    adults.Val <- sample( adults[mTab[adults, "Community"]=="Comunidad Valenciana"], 5) # (1+4 Child and Senior)
    adults.Gal <- sample( adults[mTab[adults, "Community"]=="Galicia"], 2) # (1+1 Child and Senior)
    adults.Mur <- sample( adults[mTab[adults, "Community"]=="Murcia"], 1) # (0+1 Child and Senior)
    adults.PV  <- sample( adults[mTab[adults, "Community"]=="País Vasco"], 5) # (2+3 Child and Senior)
    adults <- sort(c(adults.And, adults.Can, adults.Cat, adults.Mad, adults.Val, adults.Gal, adults.Mur, adults.PV))
    
    # ensure proportionate numbers based on gender (M~20(7+13 Child and Senior) and F~34(5+29 Child and Senior)), will allow +- 4
    if (length(teens[mTab[teens, "Gender"]=="M"]) < length(teens)*(16/54) | 
        length(teens[mTab[teens, "Gender"]=="M"]) > length(teens)*(24/54)) stop("Gender not balanced proportionately in teens")
    
    if (length(adults[mTab[adults, "Gender"]=="M"]) < length(adults)*(16/54) | 
        length(adults[mTab[adults, "Gender"]=="M"]) > length(adults)*(24/54)) stop("Gender not balanced proportionately in adults")
    # ******************************************************* #
    
  } else if (whichSamps=="healthy.nonBottle") {
    # ******************************************************* #
    
    # healthy Child and Senior (8+33=41)
    
    # prepare teen samples
    teens <- rownames(mTab[ ! is.na(mTab[,"Age_groups"]) & mTab[,"Age_groups"]=="Teen", ])
    # get proportionate numbers based on comunities
    teens.And <- sample( teens[mTab[teens, "Community"]=="Andalucía"], 7) # (1+6 Child and Senior)
    # no healthy.nonBottle were teens from Cantabria, so will skip, though 0 Child and 4 Seniors
    # teens.Can <- teens[mTab[teens, "Community"]=="Cantabria"]
    # teens.Can <- sample(teens.Can, length(teens.Can)*(5/nSubSamps)) # (0+4 Child and Senior)
    teens.Cat <- sample( teens[mTab[teens, "Community"]=="Cataluña"], 10) # (2+8 Child and Senior)
    teens.Mad <- sample( teens[mTab[teens, "Community"]=="Comunidad de Madrid"], 11) # (1+10 Child and Senior)
    teens.Val <- sample( teens[mTab[teens, "Community"]=="Comunidad Valenciana"], 2) # (1+1 Child and Senior)
    teens.Gal <- sample( teens[mTab[teens, "Community"]=="Galicia"], 1) # (1+0 Child and Senior)
    teens.Mur <- sample( teens[mTab[teens, "Community"]=="Murcia"], 1) # (0+1 Child and Senior)
    teens.PV  <- sample( teens[mTab[teens, "Community"]=="País Vasco"], 5) # (2+3 Child and Senior)
    teens <- sort(c(teens.And, teens.Cat, teens.Mad, teens.Val, teens.Gal, teens.Mur, teens.PV))
    
    # prepare adult samples
    adults <- rownames(mTab[ ! is.na(mTab[,"Age_groups"]) & mTab[,"Age_groups"]=="Adult", ])
    # get proportionate numbers based on comunities
    adults.And <- sample( adults[mTab[adults, "Community"]=="Andalucía"], 7) # (1+6 Child and Senior)
    # no healthy.nonBottle were teens from Cantabria, so will skip, though 0 Child and 4 Seniors
    #   but there were 3 adults from here, so cannot reach total of Child+Senior, instead will take 1 adult
    adults.Can <- sample( adults[mTab[adults, "Community"]=="Cantabria"], 1) # (0+4 Child and Senior)
    adults.Cat <- sample( adults[mTab[adults, "Community"]=="Cataluña"], 10) # (2+8 Child and Senior)
    adults.Mad <- sample( adults[mTab[adults, "Community"]=="Comunidad de Madrid"], 11) # (1+10 Child and Senior)
    adults.Val <- sample( adults[mTab[adults, "Community"]=="Comunidad Valenciana"], 2) # (1+1 Child and Senior)
    adults.Gal <- sample( adults[mTab[adults, "Community"]=="Galicia"], 1) # (1+0 Child and Senior)
    adults.Mur <- sample( adults[mTab[adults, "Community"]=="Murcia"], 1) # (0+1 Child and Senior)
    adults.PV  <- sample( adults[mTab[adults, "Community"]=="País Vasco"], 5) # (2+3 Child and Senior)
    adults <- sort(c(adults.And, adults.Can, adults.Cat, adults.Mad, adults.Val, adults.Gal, adults.Mur, adults.PV))
    
    # ensure proportionate numbers based on gender (M~15(5+10 Child and Senior) and F~26(3+23 Child and Senior)), will allow +- 3
    if (length(teens[mTab[teens, "Gender"]=="M"]) < length(teens)*(12/41) | 
        length(teens[mTab[teens, "Gender"]=="M"]) > length(teens)*(18/41)) stop("Gender not balanced proportionately in teens")
    
    if (length(adults[mTab[adults, "Gender"]=="M"]) < length(adults)*(12/41) | 
        length(adults[mTab[adults, "Gender"]=="M"]) > length(adults)*(18/41)) stop("Gender not balanced proportionately in adults")
    # ******************************************************* #
    
  }
  
  return(list("Teen"=teens, "Adult"=adults))
}
# ************************************************************************************************ #
# ************************************************************************************************ #
generate_ageGroup_subsamples.seniorSubs <- function(whichSamps, mTab) {
  
  # will aim to take 25 each of Teen, Adult, and Senior samples, 
  #   matching geography by taking double the amount of children in a given area, when possible
  #   when not, will compensate in other communities (** see note next to code below) - only needed for Teens
  
  
  # ******************************************************* #
  if (whichSamps=="healthy") {
    
    # prepare teen samples
    teens <- rownames(mTab[ ! is.na(mTab[,"Age_groups"]) & mTab[,"Age_groups"]=="Teen", ])
    # get proportionate numbers based on comunities
    teens.And <- sample( teens[mTab[teens, "Community"]=="Andalucía"], 4) # (2 Child)
    # no healthySamps were teens from Cantabria, so will skip, though 1 Child
    # teens.Can <- teens[mTab[teens, "Community"]=="Cantabria"]
    # teens.Can <- sample(teens.Can, length(teens.Can)*(5/nSubSamps)) # (1+4 Child and Senior)
    teens.Cat <- sample( teens[mTab[teens, "Community"]=="Cataluña"], 8) # (4 Child)
    teens.Mad <- sample( teens[mTab[teens, "Community"]=="Comunidad de Madrid"], 3) # (1 Child, **3 bc need extras, has many of each age)
    teens.Val <- sample( teens[mTab[teens, "Community"]=="Comunidad Valenciana"], 3) # (1 Child, **3 bc need extras, has many of each age)
    teens.Gal <- sample( teens[mTab[teens, "Community"]=="Galicia"], 2) # (1 Child)
    teens.Mur <- sample( teens[mTab[teens, "Community"]=="Murcia"], 1) # (0 Child, but 1 Senior)
    teens.PV  <- sample( teens[mTab[teens, "Community"]=="País Vasco"], 4) # (2 Child)
    teens <- sort(c(teens.And, teens.Cat, teens.Mad, teens.Val, teens.Gal, teens.Mur, teens.PV))
    
    
    # prepare adult samples
    adults <- rownames(mTab[ ! is.na(mTab[,"Age_groups"]) & mTab[,"Age_groups"]=="Adult", ])
    # get proportionate numbers based on comunities
    adults.And <- sample( adults[mTab[adults, "Community"]=="Andalucía"], 4) # (2 Child)
    # no healthySamps were teens from Cantabria, so will skip, though 1 Child and 4 Seniors
    #   but there were 4 adults from here, so cannot reach total of Child+Senior, instead will take 2 adults
    adults.Can <- sample( adults[mTab[adults, "Community"]=="Cantabria"], 2) # (1 Child)
    adults.Cat <- sample( adults[mTab[adults, "Community"]=="Cataluña"], 8) # (4 Child)
    adults.Mad <- sample( adults[mTab[adults, "Community"]=="Comunidad de Madrid"], 2) # (1 Child)
    adults.Val <- sample( adults[mTab[adults, "Community"]=="Comunidad Valenciana"], 2) # (1 Child)
    adults.Gal <- sample( adults[mTab[adults, "Community"]=="Galicia"], 2) # (1 Child)
    adults.Mur <- sample( adults[mTab[adults, "Community"]=="Murcia"], 1) # (0 Child, but 1 Senior)
    adults.PV  <- sample( adults[mTab[adults, "Community"]=="País Vasco"], 4) # (2 Child)
    adults <- sort(c(adults.And, adults.Can, adults.Cat, adults.Mad, adults.Val, adults.Gal, adults.Mur, adults.PV))
    
    
    # prepare senior samples
    seniors <- rownames(mTab[ ! is.na(mTab[,"Age_groups"]) & mTab[,"Age_groups"]=="Senior", ])
    # get proportionate numbers based on comunities
    seniors.And <- sample( seniors[mTab[seniors, "Community"]=="Andalucía"], 4) # (2 Child)
    # no healthySamps were teens from Cantabria, so will skip, though 1 Child and 4 Seniors
    #   but there were 4 adults from here, so cannot reach total of Child+Senior, instead will take 2 adults
    seniors.Can <- sample( seniors[mTab[seniors, "Community"]=="Cantabria"], 2) # (1 Child)
    seniors.Cat <- sample( seniors[mTab[seniors, "Community"]=="Cataluña"], 8) # (4 Child)
    seniors.Mad <- sample( seniors[mTab[seniors, "Community"]=="Comunidad de Madrid"], 2) # (1 Child)
    seniors.Val <- sample( seniors[mTab[seniors, "Community"]=="Comunidad Valenciana"], 2) # (1 Child)
    seniors.Gal <- sample( seniors[mTab[seniors, "Community"]=="Galicia"], 1) # (1 Child, but only 1 Senior)
    seniors.Mur <- sample( seniors[mTab[seniors, "Community"]=="Murcia"], 1) # (0 Child, but 1 Senior)
    seniors.PV  <- sample( seniors[mTab[seniors, "Community"]=="País Vasco"], 2) # (2 Child, but only 2 Senior)
    seniors <- sort(c(seniors.And, seniors.Can, seniors.Cat, seniors.Mad, seniors.Val, seniors.Gal, seniors.Mur, seniors.PV))
    
    # ensure proportionate numbers based on gender (M ~ 7 Child and F ~ 5 Child), will allow +- 2.5
    if (length(teens[mTab[teens, "Gender"]=="M"]) < length(teens)*(4.5/12) | 
        length(teens[mTab[teens, "Gender"]=="M"]) > length(teens)*(9.5/12)) stop("Gender not balanced proportionately in teens")
    
    if (length(adults[mTab[adults, "Gender"]=="M"]) < length(adults)*(4.5/12) | 
        length(adults[mTab[adults, "Gender"]=="M"]) > length(adults)*(9.5/12)) stop("Gender not balanced proportionately in adults")
    
    if (length(seniors[mTab[seniors, "Gender"]=="M"]) < length(seniors)*(4.5/12) | 
        length(seniors[mTab[seniors, "Gender"]=="M"]) > length(seniors)*(9.5/12)) stop("Gender not balanced proportionately in seniors")
    # ******************************************************* #
    
  }
  
  return(list("Teen"=teens, "Adult"=adults, "Senior"=seniors))
}
# ************************************************************************************************ #
# ************************************************************************************************ #
generate_ageGroup_subsamples.TAS <- function(whichSamps, mTab) {
  
  # will match Teen and Adult samples to the 42 Senior samples
  
  
  # ******************************************************* #
  if (whichSamps=="healthy") {
    
    # prepare teen samples
    teens <- rownames(mTab[ ! is.na(mTab[,"Age_groups"]) & mTab[,"Age_groups"]=="Teen", ])
    # get proportionate numbers based on comunities
    teens.And <- sample( teens[mTab[teens, "Community"]=="Andalucía"], 11) # (1 extra to match 42 Seniors total)
    # teens.Can <- sample( teens[mTab[teens, "Community"]=="Cantabria"], 1) # (only 1, though 4 Senior)
    teens.Cat <- sample( teens[mTab[teens, "Community"]=="Cataluña"], 10) # (1 extra to match 42 Seniors total)
    teens.Mad <- sample( teens[mTab[teens, "Community"]=="Comunidad de Madrid"], 11) # (1 extra to match 42 Seniors total)
    teens.Val <- sample( teens[mTab[teens, "Community"]=="Comunidad Valenciana"], 5) # (1 extra to match 42 Seniors total)
    teens.Gal <- sample( teens[mTab[teens, "Community"]=="Galicia"], 1) # 
    teens.Mur <- sample( teens[mTab[teens, "Community"]=="Murcia"], 1) # 
    teens.PV  <- sample( teens[mTab[teens, "Community"]=="País Vasco"], 3) # 
    teens <- sort(c(teens.And, teens.Cat, teens.Mad, teens.Val, teens.Gal, teens.Mur, teens.PV))
    
    
    # prepare adult samples
    adults <- rownames(mTab[ ! is.na(mTab[,"Age_groups"]) & mTab[,"Age_groups"]=="Adult", ])
    # get proportionate numbers based on comunities
    adults.And <- sample( adults[mTab[adults, "Community"]=="Andalucía"], 10) # 
    adults.Can <- sample( adults[mTab[adults, "Community"]=="Cantabria"], 4) #
    adults.Cat <- sample( adults[mTab[adults, "Community"]=="Cataluña"], 9) # 
    adults.Mad <- sample( adults[mTab[adults, "Community"]=="Comunidad de Madrid"], 10) 
    adults.Val <- sample( adults[mTab[adults, "Community"]=="Comunidad Valenciana"], 4) #
    adults.Gal <- sample( adults[mTab[adults, "Community"]=="Galicia"], 1) # 
    adults.Mur <- sample( adults[mTab[adults, "Community"]=="Murcia"], 1) # 
    adults.PV  <- sample( adults[mTab[adults, "Community"]=="País Vasco"], 3) #
    adults <- sort(c(adults.And, adults.Can, adults.Cat, adults.Mad, adults.Val, adults.Gal, adults.Mur, adults.PV))
    
    
    # ensure proportionate numbers based on gender (M ~ 13 Senior and F ~ 29 Senior), will allow +- 3
    if (length(teens[mTab[teens, "Gender"]=="M"]) < length(teens)*(10/42) | 
        length(teens[mTab[teens, "Gender"]=="M"]) > length(teens)*(16/42)) stop("Gender not balanced proportionately in teens")
    
    if (length(adults[mTab[adults, "Gender"]=="M"]) < length(adults)*(10/42) | 
        length(adults[mTab[adults, "Gender"]=="M"]) > length(adults)*(16/42)) stop("Gender not balanced proportionately in adults")
    
    # ******************************************************* #
    
  }
  
  return(list("Teen"=teens, "Adult"=adults))
}
# ************************************************************************************************ #
# ************************************************************************************************ #
generate_ageBins_subsamples.TAS <- function(whichSamps, mTab) {
  
  # will match Teen and Adult samples to the 42 Senior samples
  
  # only take non-Child samples
  mTab <- mTab[ ! is.na(mTab[,"Age_bins"]) & mTab[,"Age_bins"] != "0_13", ]
  
  bin_samplings <- list()
  # will set max number of samples to take from each community, some Age_bins have more, but will go for a more common proportion
  #   if a bin doesnt have that many from a given community, take all from that community
  nums <- c("Andalucía"=9, "Aragón"=3, "Cantabria"=1, "Cataluña"=7, "Comunidad de Madrid"=7,
            "Comunidad Valenciana"=7, "Galicia"=3, "Murcia"=3, "País Vasco"=2)
  # smaller values for the 13_20 group since they have many more samples from each region, have to cut them down a bit
  nums.13_20.40_50 <- c("Andalucía"=6, "Aragón"=3, "Cantabria"=1, "Cataluña"=6, "Comunidad de Madrid"=6,
                  "Comunidad Valenciana"=5, "Galicia"=3, "Murcia"=2, "País Vasco"=2)
  # ******************************************************* #
  
  for (bin_label in sort(unique(mTab$Age_bins))) {
    bin <- rownames(mTab[ mTab[,"Age_bins"]==bin_label, ])
    
    # prepare vector of samples to be selected from communities:
    bin.samp <- c()
    for (com in names(nums)) {
      bin_in_com <- bin[ mTab[bin, "Community"]==com ]
      # if fewer samples from bin in given community, use all those, otherwise use the indicated number of samples from that bin/community
      if (bin_label %in% c("13_20","40_50")) {
        sampSize.bin <- ifelse(length(bin_in_com) < nums.13_20.40_50[com], length(bin_in_com), nums.13_20.40_50[com])
      } else {
        sampSize.bin <- ifelse(length(bin_in_com) < nums[com], length(bin_in_com), nums[com])
      }
      bin.samp <- c(bin.samp, sample(bin_in_com, sampSize.bin))
    }
    bin.samp <- sort(bin.samp)
    
    # ensure proportionate numbers based on gender 
    # lowest proportion of males is in 60+: 13/42
    # highest proportion of males is in 30_40: 16/28
    # dont allow proportions to go wider than those
    if (length(bin.samp[mTab[bin.samp, "Gender"]=="M"]) < length(bin.samp)*(13/42) | 
        length(bin.samp[mTab[bin.samp, "Gender"]=="M"]) > length(bin.samp)*(16/28)) 
      stop(sprintf("Gender not balanced proportionately in %s", bin.samp))
    
    bin_samplings[[ bin_label ]] <- bin.samp
  }
  
  return(bin_samplings)
  
}
# ************************************************************************************************ #




# ************************************************************************************************ #
generate_waterType_subsamples <- function(whichSamps, mTab) {
  
  # will match filtered, unfiltered, and bottled samples to the 94 untreated samples
  
  
  # ******************************************************* #
  if (whichSamps=="healthy") {
    
    # prepare filtered samples
    filts <- rownames(mTab[ ! is.na(mTab[,"Age_groups"]) & ! is.na(mTab[,"Water_type_home"]) & 
                              mTab[,"Water_type_home"]=="Del Grifo (Filtrada)", ])
    # get proportionate numbers based on comunities
    filts.And <- sample( filts[mTab[filts, "Community"]=="Andalucía"], 4) #
    filts.Ara <- sample( filts[mTab[filts, "Community"]=="Aragón"], 1) #
    # filts.Can <- sample( filts[mTab[filts, "Community"]=="Cantabria"], 1) # 
    filts.Cat <- sample( filts[mTab[filts, "Community"]=="Cataluña"], 11) # 
    filts.Mad <- sample( filts[mTab[filts, "Community"]=="Comunidad de Madrid"], 3) # 
    filts.Val <- sample( filts[mTab[filts, "Community"]=="Comunidad Valenciana"], 39) # 
    filts.Gal <- sample( filts[mTab[filts, "Community"]=="Galicia"], 26) # 
    filts.IB  <- sample( filts[mTab[filts, "Community"]=="Islas Baleares"], 3) # 
    filts.Mur <- sample( filts[mTab[filts, "Community"]=="Murcia"], 4) # 
    filts.PV  <- sample( filts[mTab[filts, "Community"]=="País Vasco"], 3) # 
    filts <- sort(c(filts.And, filts.Ara, filts.Cat, filts.Mad, filts.Val, filts.Gal, filts.IB, filts.Mur, filts.PV))
    
    
    # prepare unfiltered samples
    unfilts <- rownames(mTab[ ! is.na(mTab[,"Age_groups"]) & ! is.na(mTab[,"Water_type_home"]) &
                                mTab[,"Water_type_home"]=="Del Grifo (No Filtrada)", ])
    # get proportionate numbers based on comunities
    unfilts.And <- sample( unfilts[mTab[unfilts, "Community"]=="Andalucía"], 7) #  (3 extra to match 94 untreated total)
    unfilts.Ara <- sample( unfilts[mTab[unfilts, "Community"]=="Aragón"], 3) # (2 extra to match 94 untreated total)
    # unfilts.Can <- sample( unfilts[mTab[unfilts, "Community"]=="Cantabria"], 4) #
    unfilts.Cat <- sample( unfilts[mTab[unfilts, "Community"]=="Cataluña"], 16) # (5 extra to match 94 untreated total)
    unfilts.Mad <- sample( unfilts[mTab[unfilts, "Community"]=="Comunidad de Madrid"], 4) # (1 extra to match 94 untreated total)
    unfilts.Val <- sample( unfilts[mTab[unfilts, "Community"]=="Comunidad Valenciana"], 26) # only 26 unfilts, have to make up for other 13
    unfilts.Gal <- sample( unfilts[mTab[unfilts, "Community"]=="Galicia"], 28) # (2 extra to match 94 untreated total)
    #   no Baleares unfilts, have to make up for those 3
    # unfilts.IB  <- sample( unfilts[mTab[unfilts, "Community"]=="Islas Baleares"], 3) # 
    unfilts.Mur <- sample( unfilts[mTab[unfilts, "Community"]=="Murcia"], 6) # (2 extra to match 94 untreated total)
    unfilts.PV  <- sample( unfilts[mTab[unfilts, "Community"]=="País Vasco"], 4) # (1 extra to match 94 untreated total)
    unfilts <- sort(c(unfilts.And, unfilts.Ara, unfilts.Cat, unfilts.Mad, unfilts.Val, unfilts.Gal, unfilts.Mur, unfilts.PV))
    
    
    # prepare bottled samples
    bottled <- rownames(mTab[ ! is.na(mTab[,"Age_groups"]) & ! is.na(mTab[,"Water_type_home"]) & 
                                mTab[,"Water_type_home"]=="Embotellada", ])
    # get proportionate numbers based on comunities
    bottled.And <- sample( bottled[mTab[bottled, "Community"]=="Andalucía"], 4) #
    bottled.Ara <- sample( bottled[mTab[bottled, "Community"]=="Aragón"], 1) #
    # bottled.Can <- sample( bottled[mTab[bottled, "Community"]=="Cantabria"], 1) # 
    bottled.Cat <- sample( bottled[mTab[bottled, "Community"]=="Cataluña"], 11) # 
    bottled.Mad <- sample( bottled[mTab[bottled, "Community"]=="Comunidad de Madrid"], 3) # 
    bottled.Val <- sample( bottled[mTab[bottled, "Community"]=="Comunidad Valenciana"], 39) # 
    bottled.Gal <- sample( bottled[mTab[bottled, "Community"]=="Galicia"], 26) # 
    bottled.IB  <- sample( bottled[mTab[bottled, "Community"]=="Islas Baleares"], 3) # 
    bottled.Mur <- sample( bottled[mTab[bottled, "Community"]=="Murcia"], 4) # 
    bottled.PV  <- sample( bottled[mTab[bottled, "Community"]=="País Vasco"], 3) # 
    bottled <- sort(c(bottled.And, bottled.Ara, bottled.Cat, bottled.Mad, bottled.Val, bottled.Gal, bottled.IB, bottled.Mur, bottled.PV))
    
    # ******************************************************* #
    # ensure proportionate numbers based on gender (M ~ 41 untreated and F ~ 53 untreated), will allow +- 5
    if (length(filts[mTab[filts, "Gender"]=="M"]) < length(filts)*(36/94) | 
        length(filts[mTab[filts, "Gender"]=="M"]) > length(filts)*(46/94)) stop("Gender not balanced proportionately in filts")
    
    if (length(unfilts[mTab[unfilts, "Gender"]=="M"]) < length(unfilts)*(36/94) | 
        length(unfilts[mTab[unfilts, "Gender"]=="M"]) > length(unfilts)*(46/94)) stop("Gender not balanced proportionately in unfilts")
    
    if (length(bottled[mTab[bottled, "Gender"]=="M"]) < length(bottled)*(36/94) | 
        length(bottled[mTab[bottled, "Gender"]=="M"]) > length(bottled)*(46/94)) stop("Gender not balanced proportionately in bottled")
    # ******************************************************* #
    
    untreat <- rownames(mTab[ ! is.na(mTab[,"Age_groups"]) & ! is.na(mTab[,"Water_type_home"]) &
                                mTab[,"Water_type_home"]=="No Tratada (Fuente, Pozo O Rio)", ])
    
    waterList <- list("Del Grifo (Filtrada)"=filts, "Del Grifo (No Filtrada)"=unfilts, "Embotellada"=bottled, 
                      "No Tratada (Fuente, Pozo O Rio)"=untreat)
    
    
    # ******************************************************* #
    
    # get proportionate numbers based on ages
    var.age.props <- sapply(c("Child","Teen","Adult","Senior"), function(x) 
      nrow(mTab[untreat,][mTab[untreat, "Age_groups"]==x,])/length(untreat))
    
    # ranges of accepted percentages are somewhat arbitrary, but are estimated based on total meta.healthy Age_groups proportions:
    #     Child      Teen     Adult    Senior 
    # 0.0212766 0.8191489 0.1489362 0.0106383 
    
    # give these general ranges of variation for each age group (apv = age_proportion_variation)
    apv <- c("Child"=0.05, "Teen"=0.3, "Adult"=0.3, "Senior"=0.05)
    
    for (wt in names(waterList)) {
      if (length(waterList[[ wt ]][mTab[waterList[[ wt ]], "Age_groups"]=="Child"]) < (var.age.props["Child"]-apv["Child"])*length(untreat) |
          length(waterList[[ wt ]][mTab[waterList[[ wt ]], "Age_groups"]=="Child"]) > (var.age.props["Child"]+apv["Child"])*length(untreat))
        stop("Can't balance age proportionately with this subset")
      if (length(waterList[[ wt ]][mTab[waterList[[ wt ]], "Age_groups"]=="Teen"]) < (var.age.props["Teen"]-apv["Teen"])*length(untreat) |
          length(waterList[[ wt ]][mTab[waterList[[ wt ]], "Age_groups"]=="Teen"]) > (var.age.props["Teen"]+apv["Teen"])*length(untreat))
        stop("Can't balance age proportionately with this subset")
      if (length(waterList[[ wt ]][mTab[waterList[[ wt ]], "Age_groups"]=="Adult"]) < (var.age.props["Adult"]-apv["Adult"])*length(untreat) |
          length(waterList[[ wt ]][mTab[waterList[[ wt ]], "Age_groups"]=="Adult"]) > (var.age.props["Adult"]+apv["Adult"])*length(untreat))
        stop("Can't balance age proportionately with this subset")
      if (length(waterList[[ wt ]][mTab[waterList[[ wt ]], "Age_groups"]=="Senior"]) < (var.age.props["Senior"]-apv["Senior"])*length(untreat) |
          length(waterList[[ wt ]][mTab[waterList[[ wt ]], "Age_groups"]=="Senior"]) > (var.age.props["Senior"]+apv["Senior"])*length(untreat))
        stop("Can't balance age proportionately with this subset")
    }
  }
  
  return(waterList)
}
# ************************************************************************************************ #




healthy.nonBottle <- healthySamps[ ! is.na(SLL2.meta[healthySamps, "Water_type_home"]) &
                                     SLL2.meta[healthySamps, "Water_type_home"] != "Embotellada" ]
meta.nonBotHeal <- SLL2.meta[ healthy.nonBottle, ]


# ************************************************************************************************ #
get_samps_for_test <- function(variable, mTab, i, phy=NULL, defaultDisorder=F, sameComs=F, equalSize=F,
                               printAttemptCount=T, ageGroups=F, ageGroups.seniorSubs=F, ageGroups.TAS=F,
                               ageBins.TAS=F, whichSamps=NULL, YesNo=F, YesNo.add=F, waterType=F) {
  # get appropriate subsampling of control samples for given variable
  #  - sometimes it may be that the first set of subsampling in generate_control_samples() (by community)
  #    will not leave enough samples for a given age group in the second set of subsampling (by age)
  #    so have to use the error handling below to set control samples
  controls <- NULL
  attempt <- 0
  while( is.null(controls) ) {
    attempt <- attempt + 1
    try( {
      if (sameComs == T)
        controls <- generate_control_samples.sameComs(variable, phy, mTab)
      else if (equalSize == T)
        controls <- generate_control_samples(variable, phy, mTab)
      else if (defaultDisorder == T)
        controls <- generate_control_samples.larger(variable, phy, mTab)
      else if (ageGroups == T)
        controls <- generate_ageGroup_subsamples(whichSamps, mTab)
      else if (ageGroups.seniorSubs == T)
        controls <- generate_ageGroup_subsamples.seniorSubs(whichSamps, mTab)
      else if (ageGroups.TAS == T)
        controls <- generate_ageGroup_subsamples.TAS(whichSamps, mTab)
      else if (ageBins.TAS == T)
        controls <- generate_ageBins_subsamples.TAS(whichSamps, mTab)
      else if (YesNo == T)
        controls <- generate_YesNo_subsamples(variable, mTab)
      else if (YesNo.add == T)
        controls <- generate_YesNo_subsamples.additional(variable, mTab)
      else if (waterType == T)
        controls <- generate_waterType_subsamples(whichSamps, mTab)
    },
    silent = TRUE
    )
  }
  # if (attempt>1) print(sprintf("At i = %s, required %s attempts at subsampling", i,attempt))
  if (printAttemptCount == T & sameComs ==T)
    print(sprintf("  At i = %s, required %s attempts at subsampling - sameComs", i,attempt))
  else if (printAttemptCount == T & equalSize ==T)
    print(sprintf("  At i = %s, required %s attempts at subsampling - equalSize", i,attempt))
  else if (printAttemptCount == T & defaultDisorder == T)
    print(sprintf("  At i = %s, required %s attempts at subsampling - larger", i,attempt))
  else if (printAttemptCount == T & ageGroups == T)
    print(sprintf("  At i = %s, required %s attempts at subsampling - ageGroups", i,attempt))
  else if (printAttemptCount == T & ageGroups.seniorSubs == T)
    print(sprintf("  At i = %s, required %s attempts at subsampling - ageGroups.seniorSubs", i,attempt))
  else if (printAttemptCount == T & ageGroups.TAS == T)
    print(sprintf("  At i = %s, required %s attempts at subsampling - ageGroups.TAS", i,attempt))
  else if (printAttemptCount == T & ageBins.TAS == T)
    print(sprintf("  At i = %s, required %s attempts at subsampling - ageBins.TAS", i,attempt))
  else if (printAttemptCount == T & YesNo == T)
    print(sprintf("  At i = %s, required %s attempts at subsampling - YesNo", i,attempt))
  else if (printAttemptCount == T & YesNo.add == T)
    print(sprintf("  At i = %s, required %s attempts at subsampling - YesNo.add", i,attempt))
  else if (printAttemptCount == T & waterType == T)
    print(sprintf("  At i = %s, required %s attempts at subsampling - waterType", i,attempt))
  
  return(controls)
  
}
# ************************************************************************************************ #



# test subsampling function for Age_groups:

aG <- list()
for ( i in as.character(1:100)) {
  aG[[ i ]] <- get_samps_for_test("Age_groups", meta.healthy, i, phy=NULL, ageGroups = T, whichSamps = "healthy")
}


length(unique(unlist(lapply(aG, function(x) x$Adult))))
length(unique(unlist(lapply(aG, function(x) x$Teen))))
table(meta.healthy$Age_groups)


nonAd <- rownames(meta.healthy[ ! is.na(meta.healthy$Age_groups) & meta.healthy$Age_groups == "Adult", ])
nonAd <- nonAd[ ! nonAd %in% unique(unlist(lapply(aG, function(x) x$Adult))) ]

nonTe <- rownames(meta.healthy[ ! is.na(meta.healthy$Age_groups) & meta.healthy$Age_groups == "Teen", ])
nonTe <- nonTe[ ! nonTe %in% unique(unlist(lapply(aG, function(x) x$Teen))) ]


table(meta.healthy[nonAd, "Community"])
table(meta.healthy[nonTe, "Community"])
table(meta.healthy$Community, meta.healthy$Age_groups)


# ************************************ #

aG.nb <- list()
for ( i in as.character(1:100)) {
  aG.nb[[ i ]] <- get_samps_for_test("Age_groups", meta.nonBotHeal, i, phy=NULL, ageGroups = T, whichSamps = "healthy.nonBottle")
}


length(unique(unlist(lapply(aG.nb, function(x) x$Adult))))
length(unique(unlist(lapply(aG.nb, function(x) x$Teen))))
table(meta.nonBotHeal$Age_groups)


nonAd.nb <- rownames(meta.nonBotHeal[ ! is.na(meta.nonBotHeal$Age_groups) & meta.nonBotHeal$Age_groups == "Adult", ])
nonAd.nb <- nonAd.nb[ ! nonAd.nb %in% unique(unlist(lapply(aG.nb, function(x) x$Adult))) ]

nonTe.nb <- rownames(meta.nonBotHeal[ ! is.na(meta.nonBotHeal$Age_groups) & meta.nonBotHeal$Age_groups == "Teen", ])
nonTe.nb <- nonTe.nb[ ! nonTe.nb %in% unique(unlist(lapply(aG.nb, function(x) x$Teen))) ]


table(meta.nonBotHeal[nonAd.nb, "Community"])
table(meta.nonBotHeal[nonTe.nb, "Community"])
table(meta.nonBotHeal$Community, meta.nonBotHeal$Age_groups)







# ************************************ #

aG.ss <- list()
for ( i in as.character(1:100)) {
  aG.ss[[ i ]] <- get_samps_for_test("Age_groups", meta.healthy, i, phy=NULL, ageGroups.seniorSubs = T, whichSamps = "healthy")
}


# ************************************ #

aG.tas <- list()
for ( i in as.character(1:100)) {
  aG.tas[[ i ]] <- get_samps_for_test("Age_groups", meta.healthy, i, phy=NULL, ageGroups.TAS = T, whichSamps = "healthy")
}





# ************************************ #
meta.healthy$Age.YAS <- factor(sapply(as.character(meta.healthy$Age_groups), function(x) ifelse(x %in% c("Child","Teen"),"Youth",x)),
                               levels = c("Youth","Adult","Senior"))





# ************************************ #
aB <- list()
for ( i in as.character(1:100)) {
  aB[[ i ]] <- get_samps_for_test("Age_bins", mTab.abTAS, i, phy=NULL, ageBins.TAS = T, whichSamps = "healthy")
}

length(unique(unlist(lapply(aB, function(x) x$`13_20`))))
length(unique(unlist(lapply(aB, function(x) x$`20_30`))))
length(unique(unlist(lapply(aB, function(x) x$`30_40`))))
length(unique(unlist(lapply(aB, function(x) x$`40_50`))))
length(unique(unlist(lapply(aB, function(x) x$`50_60`))))
length(unique(unlist(lapply(aB, function(x) x$`60+`))))
table(mTab.abTAS$Age_bins)


table(mTab.abTAS[unlist(aB$`1`),"Age_bins"])
table(mTab.abTAS[unlist(aB$`2`),"Age_bins"])
table(mTab.abTAS[unlist(aB$`3`),"Age_bins"])
table(mTab.abTAS[unlist(aB$`4`),"Age_bins"])
# ************************************ #






# ****************************************************************************************************************** ####
# Adonis and ANOVA for disorders (or other variables) with matched control subsamplings ####

# ********************************************************************** #
# function taken from here: 
#   https://sites.google.com/site/carlovittoriocannistraci/5-datasets-and-matlab-code/minimum-curvilinearity-ii-april-2012
library(igraph)
mce <- function(x, n, centring){
  
  # This code is performs just MCE and ncMCE (both SVD-based)
  
  #Given a distance or correlation matrix x, it performs Minimum Curvilinear 
  #Embedding (MCE) or non-centred MCE (ncMCE) (coded 24-March-2013 by 
  #Gregorio Alanis-Lobato and checked by Carlo Vittorio Cannistraci)
  
  #INPUT
  #   x => Distance (example: Euclidean) or distance-adjusted correlation matrix (example: x = 1 - Pearson_correlation)
  #   n => Dimension into which the data is to be embedded
  #   centring => 'yes' if x should be centred or 'no' if not
  #OUTPUT
  #   s => Coordinates of the samples being embedded into the reduced n-dimesional space
  
  # #Make sure the required library 'igraph' is installed and load it
  # if(require("igraph")){
  #   print("igraph has been loaded...");
  # } else{
  #   print("Trying to install igraph...");
  #   install.packages("igraph");
  #   if(require("igraph")){
  #     print("igraph has been installed and loaded...");
  #   } else{
  #     stop("Could not install igraph");
  #   }
  # }
  
  #Make sure the matrix is symmetric
  x <- pmax(x, t(x));
  
  #Create a graph object out of the adjacency matrix x
  g <- graph.adjacency(x, mode = "undirected", weighted = TRUE);
  
  #MC-kernel computation
  mst <- minimum.spanning.tree(g);
  kernel <- shortest.paths(mst);
  
  #Kernel centring
  if(centring == "yes"){
    N <- nrow(kernel);
    J <- diag(N) - (1/N)*matrix(1, N, N); #Form the centring matrix J
    kernel <- (-0.5)*(J %*% kernel^2 %*% J);
  }
  
  #SVD-based Embedding
  res <- svd(kernel);
  L <- diag(res$d);
  V <- res$v;
  
  sqrtL <- sqrt(L[1:n, 1:n]);
  V <- V[, 1:n];
  
  s <- t(sqrtL %*% t(V));
  
  return(s);
}
# ********************************************************************** #

library(fpc)
# ************************************************************ #
subsampling_ordination_objects <- function(samps, phy.subs, distsOnly=F, dists_and_ps_Only=F, pcoasOnly=F, dists_and_pcoas_only=F,
                                           print_current_dists=F, which_dists="all") {
  # function to produce distance matrices, pcoa object, and both centered and non-centered MCE coordinates
  #   for given subsamplings
  
  # remove taxa not seen at least once in at least 25 samples (2.3% of the 1085 samples).
  # This helps protect against an OTU with small mean & trivially large Coefficient of Variation.
  # from /users/tg/jwillis/SLL/Papers/Tools/phyloseq-article-source-files-figs-code-03/phyloseq_plos1_2012-source-doc.html
  phy.subs <- filter_taxa(phy.subs, function(x) sum(x > 0) > (0.025 * length(x)), prune = TRUE)
  
  # starting with distance matrices for given subsampling instance
  if (which_dists %in% c("all","Aitchison")) {
    if (print_current_dists==T) print("Aitchison")
    f.subs <- codaSeq.filter(phy.subs@otu_table, min.reads=1000, min.prop=0.001, min.occurrence=0.05, samples.by.row=FALSE)
    f.n0.subs <- cmultRepl(t(f.subs), method="CZM", label=0, suppress.print = T)
    aitch.subs <- aDist(f.n0.subs)
    # return only this distance matrix if desired
    if (which_dists == "Aitchison")
      return(aitch.subs)
  }
  
  if (which_dists %in% c("all","Weighted_Unifrac")) {
    if (print_current_dists==T) print("Weighted Unifrac")
    wu.subs <- UniFrac(phy.subs, weighted = T, parallel = T)
    # return only this distance matrix if desired
    if (which_dists == "Weighted_Unifrac")
      return(wu.subs)
  }
  
  if (which_dists %in% c("all","Unweighted_Unifrac")) {
    if (print_current_dists==T) print("Unweighted Unifrac")
    uu.subs <- UniFrac(phy.subs, weighted = F, parallel = T)
    # return only this distance matrix if desired
    if (which_dists == "Unweighted_Unifrac")
      return(uu.subs)
  }
  
  if (which_dists %in% c("all","Bray")) {
    if (print_current_dists==T) print("Bray")
    bray.subs <- vegdist(t(phy.subs@otu_table), method = "bray")
    # return only this distance matrix if desired
    if (which_dists == "Bray")
      return(bray.subs)
  }
  
  if (which_dists %in% c("all","Jaccard")) {
    if (print_current_dists==T) print("Jaccard")
    jaccard.subs <- vegdist(decostand(as.data.frame(t(phy.subs@otu_table)), method="pa"), method = "jaccard")
    # return only this distance matrix if desired
    if (which_dists == "Jaccard")
      return(jaccard.subs)
  }
  
  
  # ************************* #
  # ************************* #
  if (distsOnly == T) {
    return( list( "Aitchison"=aitch.subs, "Weighted_Unifrac"=wu.subs, "Unweighted_Unifrac"=uu.subs,
                  "Bray"=bray.subs, "Jaccard"=jaccard.subs) )
  }
  # ************************* #
  # ************************* #
  
  
  # also have to recalculate the prediction.strength values
  ps.dists.subs <- list()
  ps.dists.subs[[ "Aitchison" ]]          <- prediction.strength(aitch.subs, cutoff = 0.75)
  ps.dists.subs[[ "Weighted_Unifrac" ]]   <- prediction.strength(wu.subs, cutoff = 0.75)
  ps.dists.subs[[ "Unweighted_Unifrac" ]] <- prediction.strength(uu.subs, cutoff = 0.75)
  ps.dists.subs[[ "Bray" ]]               <- prediction.strength(bray.subs, cutoff = 0.75)
  ps.dists.subs[[ "Jaccard" ]]            <- prediction.strength(jaccard.subs, cutoff = 0.75)
  
  # ************************* #
  # ************************* #
  if (dists_and_ps_Only == T) {
    return( list( "Aitchison"=aitch.subs, "Weighted_Unifrac"=wu.subs, "Unweighted_Unifrac"=uu.subs,
                  "Bray"=bray.subs, "Jaccard"=jaccard.subs, "ps.dists"=ps.dists.subs) )
  }
  # ************************* #
  # ************************* #
  
  
  
  # then get pcoa objects
  pcoas.subs <- list()
  pcoas.subs[[ "Aitchison" ]]          <- ape::pcoa(aitch.subs)
  pcoas.subs[[ "Weighted_Unifrac" ]]   <- ape::pcoa(wu.subs)
  pcoas.subs[[ "Unweighted_Unifrac" ]] <- ape::pcoa(uu.subs)
  pcoas.subs[[ "Bray" ]]               <- ape::pcoa(bray.subs)
  pcoas.subs[[ "Jaccard" ]]            <- ape::pcoa(jaccard.subs)
  
  # ************************* #
  # ************************* #
  if (pcoasOnly == T) {
    return( list( "Aitchison"=aitch.subs, "Weighted_Unifrac"=wu.subs, "Unweighted_Unifrac"=uu.subs,
                  "Bray"=bray.subs, "Jaccard"=jaccard.subs, "pcoas"=pcoas.subs) )
  }
  # ************************* #
  # ************************* #
  if (dists_and_pcoas_only == T) {
    return( list("Aitchison"=aitch.subs, "Weighted_Unifrac"=wu.subs, "Unweighted_Unifrac"=uu.subs,
                 "Bray"=bray.subs, "Jaccard"=jaccard.subs, "ps.dists"=ps.dists.subs,
                 "pcoas"=pcoas.subs) )
  }
  # ************************* #
  # ************************* #
  
  # then non-centered MCE coordinates
  ncMCEs.subs <- list()
  ncMCEs.subs[[ "Aitchison" ]]          <- mce(as.matrix(aitch.subs), 3, "no")
  ncMCEs.subs[[ "Weighted_Unifrac" ]]   <- mce(as.matrix(wu.subs), 3, "no")
  ncMCEs.subs[[ "Unweighted_Unifrac" ]] <- mce(as.matrix(uu.subs), 3, "no")
  ncMCEs.subs[[ "Bray" ]]               <- mce(as.matrix(bray.subs), 3, "no")
  ncMCEs.subs[[ "Jaccard" ]]            <- mce(as.matrix(jaccard.subs), 3, "no")
  
  # then centered MCE coordinates
  MCEs.subs <- list()
  MCEs.subs[[ "Aitchison" ]]          <- mce(as.matrix(aitch.subs), 3, "yes")
  MCEs.subs[[ "Weighted_Unifrac" ]]   <- mce(as.matrix(wu.subs), 3, "yes")
  MCEs.subs[[ "Unweighted_Unifrac" ]] <- mce(as.matrix(uu.subs), 3, "yes")
  MCEs.subs[[ "Bray" ]]               <- mce(as.matrix(bray.subs), 3, "yes")
  MCEs.subs[[ "Jaccard" ]]            <- mce(as.matrix(jaccard.subs), 3, "yes")
  
  return( list("Aitchison"=aitch.subs, "Weighted_Unifrac"=wu.subs, "Unweighted_Unifrac"=uu.subs,
               "Bray"=bray.subs, "Jaccard"=jaccard.subs, "ps.dists"=ps.dists.subs,
               "pcoas"=pcoas.subs, "ncMCEs"=ncMCEs.subs, "MCEs"=MCEs.subs) )
}
# ************************************************************ #



# ************************************************************ #
run_full_subsampling_calcs.disorder <- function(disorder, phy, mTab, nsubs, sampMode="default", chosenControls=NULL, useAntibiotics=T) {
  
  dst <- list()
  
  disorder.samps <- rownames(mTab)[ mTab[ , disorder] == "Yes" ]
  
  # ******************************** #
  if (disorder == "Cystic_fibrosis") {
    if (useAntibiotics == T)
      cov_vars <- c("Antibiotics","Gender","Age","Population")#
    else
      cov_vars <- c("Gender","Age","Population")#
  } else {
    cov_vars <- c("Gender","Age","Population")#
  }
  # ******************************** #
  
  for ( i in as.character(1:nsubs)) {
    
    print(i)
    dst[[ i ]] <- list()
    
    # ************************************************************ #
    trait <- disorder
    tls <- c("contVar","Phylum","Class","Order","Family","Genus","Species")
    
    if ( ! is.null(chosenControls) ) {
      samps_of_interest <- chosenControls[[ i ]]
    } else {
      # first generate the control samples
      if (sampMode == "sameComs") {
        controls <- get_samps_for_test(disorder, mTab, i, phy=phy, sameComs = T)
      } else if (sampMode == "equalSize") {
        controls <- get_samps_for_test(disorder, mTab, i, phy=phy, equalSize = T)
      } else if (sampMode == "default") {
        controls <- get_samps_for_test(disorder, mTab, i, phy=phy, defaultDisorder = T)
        
      } else if (sampMode == "family") {
        controls <- NULL
        tls <- c("contVar","Phylum","Genus")
        
        if (disorder == "Downs_Syndrome")
          trait <- "DS_family"
        else if (disorder == "Cystic_fibrosis")
          trait <- "CF_family"
        
      } else if (sampMode == "useAntibiotics") {
        controls <- NULL
      }
      
      samps_of_interest <- c( disorder.samps, controls )
    }
    
    # ************************************************************ #
    dst[[ i ]][[ "samples" ]] <- samps_of_interest
    
    # make new phyloseq object with just the samples of interest
    phy.subs <- prune_samples( samps_of_interest, phy)
    mTab.subs <- mTab[ samps_of_interest, ]
    
    # get object containing the corresponding distance matrices, pcoa, and ncMCE/MCE objects
    ordObjects <- subsampling_ordination_objects( samps_of_interest, phy.subs )
    # ************************************************************ #
    
    # ************************************************************ #
    # update Stomatotype values with subsamplings
    for (dist_meas in c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard")) {
      
      clu <- get_clusters( phy.subs, "Species", dist_meas, ordObjects$ps.dists, gloms_clr, saveClusPlots = F,
                           subPops = T, subPop_Bdivs = ordObjects)
      
      #first must remove the current values of the column if they are there and being rewritten,
      # because can cause problems when adding factors with different NA values
      phy.subs@sam_data <- phy.subs@sam_data[ , colnames(phy.subs@sam_data) != sprintf("%s_%s",clu$clus_var_name, dist_meas)]
      # must specify row order in the case that not all samples were included for the clustering,
      # since the order in cluster_full.stool will be different
      phy.subs@sam_data[ names(clu$cluster_full), sprintf("%s_%s",clu$clus_var_name, dist_meas) ] <- as.factor(clu$cluster_full)
    }
    # then update mTab.subs
    mTab.subs <- meta(phy.subs)
    # and must put rows in same order as in other places
    mTab.subs <- mTab.subs[ samps_of_interest, ]
    
    # add this mTab with the new versions of Stomatotypes to the return object
    dst[[ i ]][[ "mTab" ]] <- mTab.subs
    # ************************************************************ #
    
    
    # ************************************************************ #
    dst[[ i ]][[ "Adonis" ]] <- list()
    # run adonis test for subsamples with each distance measure
    for (dist_meas in c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard")) {
      dst[[ i ]][[ "Adonis" ]][[ dist_meas ]] <- get_adonis(trait, cov_vars, dist_meas, mTab.subs)
    }
    # ************************************************************ #
    
    
    # ************************************************************ #
    # run linear model for subsamples
    dst[[ i ]][[ "Anova" ]] <- list()
    
    for (tl in tls) {
      
      # ******************** #
      if (tl == "contVar") {
        glomTab <- NULL
        dv <- c("Div.Shannon","Div.Simpson","Faiths.PD","Species_Richness","Gene_counts",
                "Stomatotype_Aitchison","Stomatotype_Weighted_Unifrac","Stomatotype_Unweighted_Unifrac",
                "Stomatotype_Bray","Stomatotype_Jaccard","pH","BMI")
        # ******************** #
      } else {
        glomTab <- gloms_clr
        dv <- NULL
      }
      # ******************** #
      
      dst[[ i ]][[ "Anova" ]][[ tl ]] <- get_lm( c(trait, cov_vars), 
                                                 tl, mTab.subs, glomTab[[ tl ]][ , samps_of_interest ], 
                                                 dv, noRemove = T, rerun.nonSig = F)
    }
    # ************************************************************ #
    
  }
  
  return(dst)
}
# ************************************************************ #



# ************************************************************ #
run_full_subsampling_calcs.ageGroups <- function(whichSamps, phy, mTab, nsubs, trait="Age_groups", chosenControls=NULL, seniorSubs=F, TAS=F, pH_BMI=F) {
  
  agst <- list()
  
  if (TAS == T)
    # in this instance, will ignore Child samples, subSamples of Teens and Adults will be matched only to Seniors
    child_senior.samps <- rownames(mTab[ ! is.na(mTab[ , "Age_groups"]) & mTab[ , "Age_groups"] %in% c("Senior"), ])
  else
    child_senior.samps <- rownames(mTab[ ! is.na(mTab[ , "Age_groups"]) & mTab[ , "Age_groups"] %in% c("Child","Senior"), ])
  
  # ******************************** #
  if (pH_BMI==T)
    cov_vars <- c("pH","BMI","Gender","Population")
  else
    cov_vars <- c("Gender","Population")
  # ******************************** #
  
  for ( i in as.character(1:nsubs)) {
    
    print(i)
    agst[[ i ]] <- list()
    
    # ************************************************************ #
    
    tls <- c("contVar","Phylum","Class","Order","Family","Genus","Species")
    
    if ( ! is.null(chosenControls) ) {
      samps_of_interest <- chosenControls[[ i ]]
    } else {
      # first generate the control samples
      if (seniorSubs == T)
        controls <- get_samps_for_test("Age_groups", mTab, i, phy=NULL, ageGroups.seniorSubs = T, whichSamps = whichSamps)
      else if (TAS == T)
        controls <- get_samps_for_test("Age_groups", mTab, i, phy=NULL, ageGroups.TAS = T, whichSamps = whichSamps)
      else
        controls <- get_samps_for_test("Age_groups", mTab, i, phy=NULL, ageGroups = T, whichSamps = whichSamps)
      
      samps_of_interest <- c( child_senior.samps, controls$Teen, controls$Adult )
    }
    
    # ************************************************************ #
    agst[[ i ]][[ "samples" ]] <- samps_of_interest
    
    # make new phyloseq object with just the samples of interest
    phy.subs <- prune_samples( samps_of_interest, phy)
    mTab.subs <- mTab[ samps_of_interest, ]
    
    # get object containing the corresponding distance matrices, pcoa, and ncMCE/MCE objects
    ordObjects <- subsampling_ordination_objects( samps_of_interest, phy.subs, dists_and_ps_Only = T )
    # ************************************************************ #
    
    # ************************************************************ #
    # update Stomatotype values with subsamplings
    for (dist_meas in c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard")) {
      
      clu <- get_clusters( phy.subs, "Species", dist_meas, ordObjects$ps.dists, gloms_clr, saveClusPlots = F,
                           subPops = T, subPop_Bdivs = ordObjects)
      
      #first must remove the current values of the column if they are there and being rewritten,
      # because can cause problems when adding factors with different NA values
      phy.subs@sam_data <- phy.subs@sam_data[ , colnames(phy.subs@sam_data) != sprintf("%s_%s",clu$clus_var_name, dist_meas)]
      # must specify row order in the case that not all samples were included for the clustering,
      # since the order in cluster_full.stool will be different
      phy.subs@sam_data[ names(clu$cluster_full), sprintf("%s_%s",clu$clus_var_name, dist_meas) ] <- as.factor(clu$cluster_full)
    }
    # then update mTab.subs
    mTab.subs <- meta(phy.subs)
    # and must put rows in same order as in other places
    mTab.subs <- mTab.subs[ samps_of_interest, ]
    
    # add this mTab with the new versions of Stomatotypes to the return object
    agst[[ i ]][[ "mTab" ]] <- mTab.subs
    # ************************************************************ #
    
    
    # ************************************************************ #
    agst[[ i ]][[ "Adonis" ]] <- list()
    # run adonis test for subsamples with each distance measure
    for (dist_meas in c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard")) {
      agst[[ i ]][[ "Adonis" ]][[ dist_meas ]] <- get_adonis(trait, cov_vars, dist_meas, mTab.subs)
    }
    # ************************************************************ #
    
    
    # ************************************************************ #
    # run linear model for subsamples
    agst[[ i ]][[ "Anova" ]] <- list()
    
    for (tl in tls) {
      
      # ******************** #
      if (tl == "contVar") {
        glomTab <- NULL
        dv <- c("Div.Shannon","Div.Simpson","Faiths.PD","Species_Richness","Gene_counts",
                "Stomatotype_Aitchison","Stomatotype_Weighted_Unifrac","Stomatotype_Unweighted_Unifrac",
                "Stomatotype_Bray","Stomatotype_Jaccard","pH","BMI")
        # ******************** #
      } else {
        glomTab <- gloms_clr
        dv <- NULL
      }
      # ******************** #
      
      agst[[ i ]][[ "Anova" ]][[ tl ]] <- get_lm( c(trait, cov_vars), 
                                                 tl, mTab.subs, glomTab[[ tl ]][ , samps_of_interest ], 
                                                 dv, noRemove = T, rerun.nonSig = F)
    }
    # ************************************************************ #
    
  }
  
  return(agst)
}
# ************************************************************ #



# ************************************************************ #
run_full_subsampling_calcs.ageBins <- function(whichSamps, phy, mTab, nsubs, trait="Age_bins", chosenControls=NULL, 
                                               seniorSubs=F, TAS=T, pH_BMI=F, adonis_only=F, anova_only=F,
                                               polyRegression=NULL, polyDegree="NULL", polyKnots="NULL",
                                               expModel=NULL, logModel=NULL, covVars=NULL, which_dists="all",
                                               dontReturnMTab=F, dontReturnSamps=F, subSampObjects=NULL) {
  
  abst <- list()
  
  # ******************************** #
  if (pH_BMI==T)
    cov_vars <- c("pH","BMI","Gender","Population")
  else if ( is.null(covVars))
    cov_vars <- c("Gender","Population")
  else
    cov_vars <- covVars
  # ******************************** #
  
  for ( i in as.character(1:nsubs)) {
    
    print(i)
    abst[[ i ]] <- list()
    
    # ************************************************************ #
    
    tls <- c("contVar","Phylum","Class","Order","Family","Genus","Species")
    
    if ( ! is.null(chosenControls) ) {
      samps_of_interest <- chosenControls[[ i ]]
    } else {
      # first generate the control samples
      controls <- get_samps_for_test("Age_bins", mTab, i, phy=NULL, ageBins.TAS = T, whichSamps = whichSamps)
      
      samps_of_interest <- unlist(controls)
    }
    
    if ( dontReturnSamps == F)
      abst[[ i ]][[ "samples" ]] <- samps_of_interest
    # ************************************************************ #
    
    
    if ( is.null(subSampObjects) ) {
      
      # make new phyloseq object with just the samples of interest
      phy.subs <- prune_samples( samps_of_interest, phy)
      mTab.subs <- mTab[ samps_of_interest, ]
      
      # get object containing the corresponding distance matrices, pcoa, and ncMCE/MCE objects
      if ( ! anova_only ) {
        ordObjects <- subsampling_ordination_objects( samps_of_interest, phy.subs, dists_and_ps_Only = T, which_dists = which_dists )
        # ************************************************************ #
        
        # ************************************************************ #
        # update Stomatotype values with subsamplings
        if (which_dists == "all")
          dists <- c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard")
        else
          dists <- which_dists
        
        for (dist_meas in dists) {
          
          clu <- get_clusters( phy.subs, "Species", dist_meas, ordObjects$ps.dists, gloms_clr, saveClusPlots = F,
                               subPops = T, subPop_Bdivs = ordObjects)
          
          #first must remove the current values of the column if they are there and being rewritten,
          # because can cause problems when adding factors with different NA values
          phy.subs@sam_data <- phy.subs@sam_data[ , colnames(phy.subs@sam_data) != sprintf("%s_%s",clu$clus_var_name, dist_meas)]
          # must specify row order in the case that not all samples were included for the clustering,
          # since the order in cluster_full.stool will be different
          phy.subs@sam_data[ names(clu$cluster_full), sprintf("%s_%s",clu$clus_var_name, dist_meas) ] <- as.factor(clu$cluster_full)
        }
        # then update mTab.subs
        mTab.subs <- meta(phy.subs)
        # and must put rows in same order as in other places
        mTab.subs <- mTab.subs[ samps_of_interest, ]
      }
      
      # add this mTab with the new versions of Stomatotypes to the return object
      if ( dontReturnMTab == F)
        abst[[ i ]][[ "mTab" ]] <- mTab.subs
      # ************************************************************ #
      
    } else {
      
      mTab.subs <- subSampObjects[[ i ]][[ "mTab.subs" ]]
      phy.subs  <- subSampObjects[[ i ]][[ "phy.subs" ]]
      ordObjects <- subSampObjects[[ i ]][[ "ordObjects" ]]
      
      if ( dontReturnMTab == F)
        abst[[ i ]][[ "mTab" ]] <- mTab.subs
      
    }
    
    # ************************************************************ #
    if ( ! anova_only ) {
      abst[[ i ]][[ "Adonis" ]] <- list()
      # run adonis test for subsamples with each distance measure
      if (which_dists == "all")
        dists <- c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard")
      else
        dists <- which_dists
      
      for (dist_meas in dists) {
        if (length(dists) == 1)
          oO <- ordObjects
        else
          oO <- ordObjects[[ dist_meas ]]
        
        abst[[ i ]][[ "Adonis" ]][[ dist_meas ]] <- get_adonis(trait, cov_vars, dist_meas, mTab.subs, ordObj = oO)
      }
    }
    # ************************************************************ #
    
    
    # ************************************************************ #
    if ( ! adonis_only ) {
      # run linear model for subsamples
      abst[[ i ]][[ "Anova" ]] <- list()
      
      for (tl in tls) {
        
        # ******************** #
        if (tl == "contVar") {
          glomTab <- NULL
          if (anova_only)
            dv <- c("Div.Shannon","Div.Simpson","Faiths.PD","Species_Richness","Gene_counts","pH","BMI")
          else 
            dv <- c("Div.Shannon","Div.Simpson","Faiths.PD","Species_Richness","Gene_counts",
                    "Stomatotype_Aitchison","Stomatotype_Weighted_Unifrac","Stomatotype_Unweighted_Unifrac",
                    "Stomatotype_Bray","Stomatotype_Jaccard","pH","BMI")
          
          if (trait %in% c("BMI","BMI_official","BMI_group"))
            dv <- dv[ dv != "BMI" ]
          if (trait %in% c("pH"))
            dv <- dv[ dv != "pH" ]
          
          # ******************** #
        } else {
          glomTab <- gloms_clr
          dv <- NULL
        }
        # ******************** #
        
        abst[[ i ]][[ "Anova" ]][[ tl ]] <- get_lm( c(trait, cov_vars), 
                                                    tl, mTab.subs, glomTab[[ tl ]][ , samps_of_interest ], 
                                                    dv, noRemove = T, rerun.nonSig = F,
                                                    polyRegression = polyRegression, polyKnots = polyKnots, polyDegree = polyDegree,
                                                    expModel = expModel, logModel = logModel)
      }
    }
    # ************************************************************ #
    
  }
  
  return(abst)
}
# ************************************************************ #

# ************************************************************ #
get_subsamp_objects <- function(chosenControls, phy, mTab, nsubs) {
  
  abst <- list()
  
  for ( i in as.character(1:nsubs)) {
    
    print(i)
    abst[[ i ]] <- list()
    
    # ************************************************************ #
    
    samps_of_interest <- chosenControls[[ i ]]
    
    # ************************************************************ #
    abst[[ i ]][[ "samples" ]] <- samps_of_interest
    
    # make new phyloseq object with just the samples of interest
    phy.subs <- prune_samples( samps_of_interest, phy)
    mTab.subs <- mTab[ samps_of_interest, ]
    
    # get object containing the corresponding distance matrices, pcoa, and ncMCE/MCE objects
    ordObjects <- subsampling_ordination_objects( samps_of_interest, phy.subs, dists_and_ps_Only = T )
    abst[[ i ]][[ "ordObjects" ]]  <- ordObjects
    
    # ************************************************************ #
    
    # ************************************************************ #
    # update Stomatotype values with subsamplings
    for (dist_meas in c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard")) {
      
      clu <- get_clusters( phy.subs, "Species", dist_meas, ordObjects$ps.dists, gloms_clr, saveClusPlots = F,
                           subPops = T, subPop_Bdivs = ordObjects)
      
      #first must remove the current values of the column if they are there and being rewritten,
      # because can cause problems when adding factors with different NA values
      phy.subs@sam_data <- phy.subs@sam_data[ , colnames(phy.subs@sam_data) != sprintf("%s_%s",clu$clus_var_name, dist_meas)]
      # must specify row order in the case that not all samples were included for the clustering,
      # since the order in cluster_full.stool will be different
      phy.subs@sam_data[ names(clu$cluster_full), sprintf("%s_%s",clu$clus_var_name, dist_meas) ] <- as.factor(clu$cluster_full)
    }
    # then update mTab.subs
    mTab.subs <- meta(phy.subs)
    # and must put rows in same order as in other places
    mTab.subs <- mTab.subs[ samps_of_interest, ]
    
    abst[[ i ]][[ "phy.subs" ]]  <- phy.subs
    abst[[ i ]][[ "mTab.subs" ]] <- mTab.subs
    # ************************************************************ #
  }
  
  return(abst)
}
# ************************************************************ #


# ************************************************************ #
get_subs_anova_per_taxon <- function(subsTests, variable, tl, calc) {
  if (calc == "mean") {
    sapply(rownames(subsTests$`1`$Anova[[ tl ]]), function(x)
      mean(p.adjust(unlist(lapply(subsTests, function(y)
        y$Anova[[ tl ]][ gsub("-","\\.",x), variable ]
      )), method = "fdr"))
    )
  } else if (calc == "num_sig") {
    sapply(rownames(subsTests$`1`$Anova[[ tl ]]), function(x)
      sum(p.adjust(unlist(lapply(subsTests, function(y)
        y$Anova[[ tl ]][ gsub("-","\\.",x), variable ]
      )), method = "fdr") < 0.05)
    )
  }
}
# ************************************************************ #



# ************************************************************ #
run_pairwise_subsampling_calcs.ageGroups <- function(whichSamps, phy, mTab, nsubs, test, chosenControls=NULL, TAS=F, which_dists="all") {
  
  agst <- list()
  
  if (TAS == T) {
    # in this instance, will ignore Child samples, subSamples of Teens and Adults will be matched only to Seniors
    child_senior.samps <- rownames(mTab[ ! is.na(mTab[ , "Age_groups"]) & mTab[ , "Age_groups"] %in% c("Senior"), ])
    ag1_vals <- c("Teen","Adult")
    ag2_vals <- c("Adult","Senior")
  } else {
    child_senior.samps <- rownames(mTab[ ! is.na(mTab[ , "Age_groups"]) & mTab[ , "Age_groups"] %in% c("Child","Senior"), ])
    ag1_vals <- c("Child","Teen","Adult")
    ag2_vals <- c("Teen","Adult","Senior")
  }
  
  # ******************************** #
  cov_vars <- c("Gender","Population")
  # ******************************** #
  
  
  for (ag1 in ag1_vals) {
    for (ag2 in ag2_vals) {
      
      agPair <- sprintf("%s-%s", ag1, ag2)
      
      if (ag1 == ag2 | 
          ( ! is.null(names(agst)) && 
            2 %in% unlist(lapply(strsplit(names(agst), "-"), function(x) sum(c(ag1,ag2) %in% x)))
            )) {
        # skip if same values or if combo already done
        next
      } 
        
      print(agPair)
      agst[[ agPair ]] <- list()
      
      # ******************************** #
      if (ag1=="Child" & ag2=="Senior") {
        # since these do not have subsamples, only run this combo once
        nsubs.agp <- 1
      } else {
        nsubs.agp <- nsubs
      }
      # ******************************** #
      
      for ( i in as.character(1:nsubs.agp)) {
        
        print(i)
        agst[[ agPair ]][[ i ]] <- list()
        
        # ************************************************************ #
        trait <- "Age_groups"
        
        if ( ! is.null(chosenControls) ) {
          samps_of_interest <- chosenControls[[ i ]]
        } else {
          # first generate the control samples
          if (TAS == T)
            controls <- get_samps_for_test("Age_groups", mTab, i, phy=NULL, ageGroups.TAS = T, whichSamps = whichSamps)
          else
            controls <- get_samps_for_test("Age_groups", mTab, i, phy=NULL, ageGroups = T, whichSamps = whichSamps)
          
          samps_of_interest <- c( child_senior.samps, controls$Teen, controls$Adult )
        }
        
        # use only those samples that are among the 2 values
        samps_of_interest <- samps_of_interest[ mTab[ samps_of_interest, "Age_groups"] %in% c(ag1, ag2) ]
        # ************************************************************ #
        agst[[ agPair ]][[ i ]][[ "samples" ]] <- samps_of_interest
        
        # make new phyloseq object with just the samples of interest
        phy.subs <- prune_samples( samps_of_interest, phy)
        mTab.subs <- mTab[ samps_of_interest, ]
        
        # get object containing the corresponding distance matrices, pcoa, and ncMCE/MCE objects
        ordObjects <- subsampling_ordination_objects( samps_of_interest, phy.subs, distsOnly = T, which_dists = which_dists )
        
        # add this mTab with the new versions of Stomatotypes to the return object
        agst[[ agPair ]][[ i ]][[ "mTab" ]] <- mTab.subs
        # ************************************************************ #
        
        
        # ************************************************************ #
        if (test == "Adonis") {
          
          agst[[ agPair ]][[ i ]][[ "Adonis" ]] <- list()
          # run adonis test for subsamples with each distance measure
          if (which_dists == "all")
            dists <- c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard")
          else
            dists <- which_dists
          
          for (dist_meas in dists) {
            if (length(dists) == 1)
              oO <- ordObjects
            else
              oO <- ordObjects[[ dist_meas ]]
            
            agst[[ agPair ]][[ i ]][[ "Adonis" ]][[ dist_meas ]] <- get_adonis(trait, cov_vars, dist_meas, mTab.subs,
                                                                               ordObj = oO)
          }
          # ******************************** #
        } else if (test == "ANOSIM") {
          
          if (trait %in% c("seqGroup")) { # and will have to add some disorder groups too, since controls only taken from seqGroup "One"
            agst[[ agPair ]][[ i ]][[ "ANOSIM" ]] <- anosim(as.dist( as.matrix(ordObjects$Aitchison)[samps_of_interest, samps_of_interest] ),
                                                            as.character(mTab.subs[ samps_of_interest, trait]))
          } else {
            agst[[ agPair ]][[ i ]][[ "ANOSIM" ]] <- anosim(as.dist( as.matrix(ordObjects$Aitchison)[samps_of_interest, samps_of_interest] ),
                                                            as.character(mTab.subs[ samps_of_interest, trait]),
                                                            strata = mTab.subs[ samps_of_interest, "seqGroup"])
          }
          
        }
        # ************************************************************ #
        
      }
    }
  }
  
  
  return(agst)
}
# ************************************************************ #


# ************************************************************ #
run_pairwise_subsampling_calcs.ageBins <- function(whichSamps, phy, mTab, nsubs, test, chosenControls=NULL, seniorSubs=F, 
                                                   TAS=T, pH_BMI=F, which_dists="all") {
  
  abst <- list()
  
  ag1_vals <- c("13_20","20_30","30_40","40_50","50_60")
  ag2_vals <- c("20_30","30_40","40_50","50_60","60+")
  
  # ******************************** #
  if (pH_BMI==T)
    cov_vars <- c("pH","BMI","Gender","Population")
  else
    cov_vars <- c("Gender","Population")
  # ******************************** #
  
  
  for (ag1 in ag1_vals) {
    for (ag2 in ag2_vals) {
      
      agPair <- sprintf("%s-%s", ag1, ag2)
      
      if (ag1 == ag2 | 
          ( ! is.null(names(abst)) && 
            2 %in% unlist(lapply(strsplit(names(abst), "-"), function(x) sum(c(ag1,ag2) %in% x)))
          )) {
        # skip if same values or if combo already done
        next
      } 
      
      print(agPair)
      abst[[ agPair ]] <- list()
      
      # ******************************** #
      if (ag1=="Child" & ag2=="Senior") {
        # since these do not have subsamples, only run this combo once
        nsubs.agp <- 1
      } else {
        nsubs.agp <- nsubs
      }
      # ******************************** #
      
      for ( i in as.character(1:nsubs.agp)) {
        
        print(i)
        abst[[ agPair ]][[ i ]] <- list()
        
        # ************************************************************ #
        trait <- "Age_bins"
        
        if ( ! is.null(chosenControls) ) {
          samps_of_interest <- chosenControls[[ i ]]
        } else {
          # first generate the control samples
          controls <- get_samps_for_test("Age_bins", mTab, i, phy=NULL, ageBins.TAS = T, whichSamps = whichSamps)
          
          samps_of_interest <- unlist(controls)
        }
        
        # use only those samples that are among the 2 values
        samps_of_interest <- samps_of_interest[ mTab[ samps_of_interest, "Age_bins"] %in% c(ag1, ag2) ]
        # ************************************************************ #
        abst[[ agPair ]][[ i ]][[ "samples" ]] <- samps_of_interest
        
        # make new phyloseq object with just the samples of interest
        phy.subs <- prune_samples( samps_of_interest, phy)
        mTab.subs <- mTab[ samps_of_interest, ]
        
        # get object containing the corresponding distance matrices, pcoa, and ncMCE/MCE objects
        ordObjects <- subsampling_ordination_objects( samps_of_interest, phy.subs, distsOnly = T, which_dists = which_dists )
        # ************************************************************ #
        
        # add this mTab with the new versions of Stomatotypes to the return object
        abst[[ agPair ]][[ i ]][[ "mTab" ]] <- mTab.subs
        # ************************************************************ #
        
        
        # ************************************************************ #
        if (test == "Adonis") {
          
          abst[[ agPair ]][[ i ]][[ "Adonis" ]] <- list()
          # run adonis test for subsamples with each distance measure
          if (which_dists == "all")
            dists <- c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard")
          else
            dists <- which_dists
          
          for (dist_meas in dists) {
            if (length(dists) == 1)
              oO <- ordObjects
            else
              oO <- ordObjects[[ dist_meas ]]
            
            abst[[ agPair ]][[ i ]][[ "Adonis" ]][[ dist_meas ]] <- get_adonis(trait, cov_vars, dist_meas, mTab.subs, 
                                                                               ordObj = oO)
          }
          # ******************************** #
        } else if (test == "ANOSIM") {
          
          if (trait %in% c("seqGroup")) { # and will have to add some disorder groups too, since controls only taken from seqGroup "One"
            abst[[ agPair ]][[ i ]][[ "ANOSIM" ]] <- anosim(as.dist( as.matrix(ordObjects$Aitchison)[samps_of_interest, samps_of_interest] ),
                                                            as.character(mTab.subs[ samps_of_interest, trait]))
          } else {
            abst[[ agPair ]][[ i ]][[ "ANOSIM" ]] <- anosim(as.dist( as.matrix(ordObjects$Aitchison)[samps_of_interest, samps_of_interest] ),
                                                            as.character(mTab.subs[ samps_of_interest, trait]),
                                                            strata = mTab.subs[ samps_of_interest, "seqGroup"])
          }
          
        }
        
        # ************************************************************ #
        
        
      }
      
    }
  }
  
  return(abst)
}
# ************************************************************ #









# ************************************************************ #
run_subsampling_anosims <- function(variable, phy, mTab, nsubs, chosenControls=NULL, YesNo.add=F, waterType=F) {
  
  ano_subs <- list()
  
  if (YesNo.add == F & waterType == F) {
    # get Yes samples first, then match controls with those
    YesSamps <- rownames(mTab[ ! is.na(mTab[,"Age"]) & ! is.na(mTab[,variable]) & mTab[,variable]=="Yes", ])
    # have to use the "No" as main group here since there are many more "Yes" samples
    if (variable %in% c("Fluoride_toothpaste","Wash_hands_before_eat","Wash_hands_after_bathroom","Do_you_feel_well"))
      YesSamps <- rownames(mTab[ ! is.na(mTab[,"Age"]) & ! is.na(mTab[,variable]) & mTab[,variable]=="No", ])
  }
  
  # ******************************** #
  
  for ( i in as.character(1:nsubs)) {
    
    # print(i)
    ano_subs[[ i ]] <- list()
    
    # ************************************************************ #
    tls <- c("contVar","Phylum","Class","Order","Family","Genus","Species")
    
    if ( ! is.null(chosenControls) ) {
      samps_of_interest <- chosenControls[[ i ]][[ "samples" ]]
    } else {
      # first generate the control samples
      if (YesNo.add == T) {
        controls <- get_samps_for_test(variable, mTab, i, YesNo.add = T)
        samps_of_interest <- c(controls$YFT, controls$NMO)
      } else if (waterType == T) {
        controls <- get_samps_for_test(variable, mTab, i, waterType = T, whichSamps = "healthy")
        samps_of_interest <- c(controls[[ "Del Grifo (Filtrada)" ]], 
                               controls[[ "Del Grifo (No Filtrada)" ]],
                               controls[[ "Embotellada" ]],
                               controls[[ "No Tratada (Fuente, Pozo O Rio)" ]])
      } else if (YesNo.add == F) {
        controls <- get_samps_for_test(variable, mTab, i, YesNo = T)
        samps_of_interest <- c( YesSamps, controls )
      }
      
    }
    
    # ************************************************************ #
    ano_subs[[ i ]][[ "samples" ]] <- samps_of_interest
    
    # make new phyloseq object with just the samples of interest
    phy.subs <- prune_samples( samps_of_interest, phy)
    mTab.subs <- mTab[ samps_of_interest, ]

    # get object containing the corresponding distance matrices, pcoa, and ncMCE/MCE objects
    ordObject <- subsampling_ordination_objects( samps_of_interest, phy.subs, which_dists = "Aitchison" )

    ano_subs[[ i ]][[ "ordObjects" ]] <- ordObject
    # ************************************************************ #



    # ************************************************************ #
    if (variable %in% c("seqGroup")) { # and will have to add some disorder groups too, since controls only taken from seqGroup "One"
      ano_subs[[ i ]][[ "ANOSIM" ]] <- anosim(as.dist( as.matrix(ordObject)[samps_of_interest, samps_of_interest] ),
                                              as.character(mTab.subs[ samps_of_interest, variable]))
    } else {
      ano_subs[[ i ]][[ "ANOSIM" ]] <- anosim(as.dist( as.matrix(ordObject)[samps_of_interest, samps_of_interest] ),
                                              as.character(mTab.subs[ samps_of_interest, variable]),
                                              strata = mTab.subs[ samps_of_interest, "seqGroup"])
    }
    
    # ************************************************************ #
    
  }
  
  return(ano_subs)
}
# ************************************************************ #


# ************************************************************ #
run_subsampling_adonis <- function(variable, phy, mTab, nsubs, chosenControls=NULL, chosenOrdObj=NULL, YesNo.add=F, which_dists="all") {
  
  ado_subs <- list()
  
  if (YesNo.add == F) {
    # get Yes samples first, then match controls with those
    YesSamps <- rownames(mTab[ ! is.na(mTab[,"Age"]) & ! is.na(mTab[,variable]) & mTab[,variable]=="Yes", ])
    # have to use the "No" as main group here since there are many more "Yes" samples
    if (variable %in% c("Fluoride_toothpaste","Wash_hands_before_eat","Wash_hands_after_bathroom","Do_you_feel_well"))
      YesSamps <- rownames(mTab[ ! is.na(mTab[,"Age"]) & ! is.na(mTab[,variable]) & mTab[,variable]=="No", ])
  }
  
  # ******************************** #
  if (variable %in% c("Age_groups","Age_groups.TAS","Age_bins",
                      "Children","Teens","Adults","Seniors",
                      "13_20","20_30","30_40","40_50","50_60","60+"))
    cov_vars <- c("Gender","Population")
  else if (variable=="Gender")
    cov_vars <- c("Age","Population")
  else
    cov_vars <- c("Age","Gender","Population")
  # ******************************** #
  
  # ******************************** #
  
  for ( i in as.character(1:nsubs)) {
    
    print(i)
    ado_subs[[ i ]] <- list()
    
    # ************************************************************ #
    tls <- c("contVar","Phylum","Class","Order","Family","Genus","Species")
    
    if ( ! is.null(chosenControls) ) {
      samps_of_interest <- chosenControls[[ i ]]
    } else {
      # first generate the control samples
      if (YesNo.add == F) {
        controls <- get_samps_for_test(variable, mTab, i, YesNo = T)
        samps_of_interest <- c( YesSamps, controls )
      } else if (YesNo.add == T) {
        controls <- get_samps_for_test(variable, mTab, i, YesNo.add = T)
        samps_of_interest <- c(controls$YFT, controls$NMO)
      }
      
    }
    
    # ************************************************************ #
    ado_subs[[ i ]][[ "samples" ]] <- samps_of_interest
    
    # make new phyloseq object with just the samples of interest
    phy.subs <- prune_samples( samps_of_interest, phy)
    mTab.subs <- mTab[ samps_of_interest, ]
    
    # get object containing the corresponding distance matrices, pcoa, and ncMCE/MCE objects
    if ( ! is.null(chosenOrdObj)) {
      
      # get distance matrices with just the samples of interest
      ordObject <- lapply(chosenOrdObj, function(x) as.dist(x[samps_of_interest, samps_of_interest]))
    } else {
      ordObject <- subsampling_ordination_objects( samps_of_interest, phy.subs, distsOnly = T, print_current_dists = F, which_dists = which_dists )
    }
    
    # ************************************************************ #
    
    
    
    # ************************************************************ #
    ado_subs[[ i ]][[ "Adonis" ]] <- list()
    # run adonis test for subsamples with each distance measure
    if (which_dists == "all")
      dists <- c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard")
    else
      dists <- which_dists
    
    for (dist_meas in dists) {
      # print(sprintf("%s - %s", dist_meas, Sys.time()))
      if (length(dists) == 1 & is.null(chosenOrdObj))
        oO <- ordObject
      else
        oO <- ordObject[[ dist_meas ]]
      
      ado_subs[[ i ]][[ "Adonis" ]][[ dist_meas ]] <- get_adonis(variable, cov_vars, dist_meas, mTab.subs, 
                                                                 ordObj = oO)
    }
    # ************************************************************ #
    
  }
  
  return(ado_subs)
}
# ************************************************************ #


# ************************************************************ #
library(fpc)
run_subsampling_Stomatotypes <- function(variable, phy, mTab, nsubs, ordsList, glomTab=gloms_clr, chosenControls=NULL, useAntibiotics=F) {
  
  stomato_subs <- list()
  
  # get Yes samples first, then match controls with those
  YesSamps <- rownames(mTab[ ! is.na(mTab[,"Age"]) & ! is.na(mTab[,variable]) & mTab[,variable]=="Yes", ])
  # have to use the "No" as main group here since there are many more "Yes" samples
  if (variable %in% c("Fluoride_toothpaste","Wash_hands_before_eat","Wash_hands_after_bathroom","Do_you_feel_well"))
    YesSamps <- rownames(mTab[ ! is.na(mTab[,"Age"]) & ! is.na(mTab[,variable]) & mTab[,variable]=="No", ])
  
  # ******************************** #
  if (variable == "Cystic_fibrosis") {
    if (useAntibiotics == T)
      cov_vars <- c("Antibiotics","Gender","Age","Population")#
    else
      cov_vars <- c("Gender","Age","Population")#
  } else if (variable == "Age_groups") {
    cov_vars <- c("Gender","Population")
  } else if (variable == "Gender") {
    cov_vars <- c("Age","Population")#
  } else {
    cov_vars <- c("Gender","Age","Population")#
  }
  # dont repeat cov_vars
  cov_vars <- cov_vars[ cov_vars != variable ]
  # ******************************** #
  
  
  
  # ******************************** #
  
  for ( i in as.character(1:nsubs)) {
    
    # print(i)
    stomato_subs[[ i ]] <- list()
    
    # ************************************************************ #
    
    if ( ! is.null(chosenControls) ) {
      samps_of_interest <- chosenControls[[ i ]]
    } else {
      # first generate the control samples
      controls <- get_samps_for_test(variable, mTab, i, YesNo = T)
      
      samps_of_interest <- c( YesSamps, controls )
    }
    
    stomato_subs[[ i ]][[ "samples" ]] <- samps_of_interest
    # ************************************************************ #
    
    # make new phyloseq object with just the samples of interest
    phy.subs <- prune_samples( samps_of_interest, phy)
    mTab.subs <- mTab[ samps_of_interest, ]
    
    # ************************************************************ #
    
    
    
    # ************************************************************ #
    # ************************************************************ #
    
    # get distance matrices with just the samples of interest
    ords.subs <- lapply(ordsList, function(x) as.dist(x[samps_of_interest, samps_of_interest]))
    
    # also have to recalculate the prediction.strength values
    ps.dists.subs <- list()
    ps.dists.subs[[ "Aitchison" ]]          <- prediction.strength(ords.subs$Aitchison, cutoff = 0.75, Gmax = 5, M=20)
    ps.dists.subs[[ "Weighted_Unifrac" ]]   <- prediction.strength(ords.subs$Weighted_Unifrac, cutoff = 0.75, Gmax = 5, M=20)
    ps.dists.subs[[ "Unweighted_Unifrac" ]] <- prediction.strength(ords.subs$Unweighted_Unifrac, cutoff = 0.75, Gmax = 5, M=20)
    ps.dists.subs[[ "Bray" ]]               <- prediction.strength(ords.subs$Bray, cutoff = 0.75, Gmax = 5, M=20)
    ps.dists.subs[[ "Jaccard" ]]            <- prediction.strength(ords.subs$Jaccard, cutoff = 0.75, Gmax = 5, M=20)
    
    # update Stomatotype values with subsamplings
    for (dist_meas in c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard")) {
      
      clu <- get_clusters( phy.subs, "Species", dist_meas, ps.dists.subs, glomTab, saveClusPlots = F,
                           subPops = T, subPop_Bdivs = ords.subs)
      #add stomatotype to sample_data
      mTab.subs[ names(clu$cluster_full), sprintf("%s_%s", clu$clus_var_name, dist_meas) ] <- clu$cluster_full#as.factor(cluster_full)
      
    }
    
    # add this mTab with the new versions of Stomatotypes to the return object
    stomato_subs[[ i ]][[ "mTab" ]] <- mTab.subs
    # ************************************************************ #
    
    # Will run 2 tests to check for association between Stomatotypes and the given variable
    # ****************** #
    # (1) Chi-square test
    stomato_subs[[ i ]][[ "Chi_Square" ]] <- list()
    for (dist_meas in c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard")) {
      
      stomato_subs[[ i ]][[ "Chi_Square" ]][[ dist_meas ]] <- chisq.test(mTab.subs[ , variable], 
                                                                         mTab.subs[ , sprintf("Stomatotype_%s", dist_meas)])
    }
    
    # ****************** #
    # (2) Multinomial Log-linear model
    dv <- c("Stomatotype_Aitchison","Stomatotype_Weighted_Unifrac","Stomatotype_Unweighted_Unifrac",
            "Stomatotype_Bray","Stomatotype_Jaccard")
    stomato_subs[[ i ]][[ "multinom" ]] <- get_lm( c(variable, cov_vars), 
                                                   "contVar", mTab.subs, NULL, 
                                                   dv, noRemove = T, rerun.nonSig = F)
    # ************************************************************ #
    
  }
  
  return(stomato_subs)
}
# ************************************************************ #




# ************************************************************ #
run_functional_lms <- function(variable, phy, meta_path, nsubs, pVars, chosenControls=NULL, useAntibiotics=F, YesNo.add=F, silentAnova=F,
                               polyRegression=NULL, polyDegree="NULL", polyKnots="NULL", expModel=NULL, logModel=NULL, useRescale=F) {
  
  if (YesNo.add == F) {
    # get Yes samples first, then match controls with those
    YesSamps <- rownames(meta_path[ ! is.na(meta_path[,"Age"]) & ! is.na(meta_path[,variable]) & meta_path[,variable]=="Yes", ])
    # have to use the "No" as main group here since there are many more "Yes" samples
    if (variable %in% c("Fluoride_toothpaste","Wash_hands_before_eat","Wash_hands_after_bathroom","Do_you_feel_well"))
      YesSamps <- rownames(meta_path[ ! is.na(meta_path[,"Age"]) & ! is.na(meta_path[,variable]) & meta_path[,variable]=="No", ])
  }
  
  # ******************************** #
  if (variable == "Cystic_fibrosis") {
    if (useAntibiotics == T)
      cov_vars <- c("Antibiotics","Gender","Age","Population")#
    else
      cov_vars <- c("Gender","Age","Population")#
  } else if (variable == "Age_groups") {
    cov_vars <- c("Gender","Population")
  } else if (variable == "Gender") {
    cov_vars <- c("Age","Population")#
  } else {
    cov_vars <- c("Gender","Age","Population")#
  }
  # dont repeat cov_vars
  cov_vars <- cov_vars[ cov_vars != variable ]
  # ******************************** #
  
  # ******************************** #
  
  path_subs <- list()
  
  for ( i in as.character(1:nsubs)) {
    
    print(i)
    path_subs[[ i ]] <- list()
    
    # ************************************************************ #
    
    if ( ! is.null(chosenControls) ) {
      samps_of_interest <- chosenControls[[ i ]]
    } else {
      # first generate the control samples
      if (YesNo.add == F) {
        controls <- get_samps_for_test(variable, meta_path, i, YesNo = T)
        samps_of_interest <- c( YesSamps, controls )
      } else if (YesNo.add == T) {
        controls <- get_samps_for_test(variable, meta_path, i, YesNo.add = T)
        samps_of_interest <- c(controls$YFT, controls$NMO)
      }
      
    }
    
    path_subs[[ i ]][[ "samples" ]] <- samps_of_interest
    # ************************************************************ #
    
    # make new phyloseq object with just the samples of interest
    phy.subs <- prune_samples( samps_of_interest, phy )
    meta_path.subs <- meta_path[ samps_of_interest, ]
    # ************************************************************ #
    
    # ************************************************************ #
    dv_toSkip <- NULL
    # check for need to rescale values (see description of rescale_depVar in get_lm())
    pathsOnly <- meta_path.subs[ , pVars ]
    smallest_non0 <- min(pathsOnly[ pathsOnly > 0 ])
    
    if (useRescale==T & smallest_non0 < 1e-4) {
      rescale_depVar <- 1*10^( -ceiling(log10( smallest_non0 )) ) # use ceiling to scale by 1 magnitude less than the inverse of the min
      print(sprintf("need to rescale, smallest_non0 = %s, rescale_depVar = %s", smallest_non0, rescale_depVar))
    } else {
      rescale_depVar <- NULL
    }
    # ****************** #
    path_subs[[ i ]][[ "Anova" ]] <- get_lm( c(variable, cov_vars), 
                                             "moreVars", meta_path.subs, NULL, pVars, 
                                             noRemove = T, rerun.nonSig = F, silentAnova = silentAnova, dv_toSkip=dv_toSkip,
                                             polyRegression = polyRegression, polyKnots = polyKnots, polyDegree = polyDegree,
                                             expModel = expModel, logModel = logModel, rescale_depVar=rescale_depVar)
    
    # ************************************************************ #
    
  }
  
  return(path_subs)
  
}
# ************************************************************ #


# ************************************************************ #
plot_freq_sigs <- function(variable, num_sigs, pMeans, minSig, tl, mTab, phy, subsTests, 
                           singleSubSamp=NULL, chosenSamps=NULL, adjustAlphas=F,
                           facetVar=NULL, facetRows=NULL) {
  
  if (tl == "Genus") {
    # only those for which at least 95/100 subsamplings were signif (ie p < 0.05)
    freqSig <- num_sigs[ num_sigs >= minSig]
    sort(pMeans[ names(freqSig) ])
    
    gensToCheck <- unique(names(freqSig))
    
    freq.mean.ps <- as.data.frame(cbind(num_sigs[ gensToCheck ], pMeans[ gensToCheck ]))
    colnames(freq.mean.ps) <- c("num_sig","meanP")
    freq.mean.ps <- freq.mean.ps[ rev(order(freq.mean.ps$num_sig, -freq.mean.ps$meanP)), ]
    freq.mean.ps
    print(freq.mean.ps)
    
    # subsampling.plots.box(variable, "Genus", variable, c("Lactobacillus","Fusobacterium","Porphyromonas","Prevotella","Leptotrichia","Veillonella","Neisseria","Rothia","Streptococcus","Alloprevotella"),
    subsampling.plots.box(variable, "Genus", variable, rownames(freq.mean.ps),
                          c(groupQs), only_cont, gloms_clr, mTab, phy, subsTests, dstStruc = "Age_groups", 
                          plotType = "box", plot_tukey = F, xAngle = 30, adjustAlphas = adjustAlphas,
                          chosenSamps = chosenSamps, singleSubSamp = singleSubSamp, 
                          facetVar = facetVar, facetRows = facetRows)
    
  } else if (tl == "Phylum") {
    
    phy.sig.order <- num_sigs[ rev(order(num_sigs, -pMeans)) ]
    phy.sig.order <- phy.sig.order[ pMeans[ names(phy.sig.order) ] < 0.05 ]
    print(phy.sig.order)
    subsampling.plots.box(variable,"Phylum",variable, names(phy.sig.order), c(groupQs), only_cont,
                          gloms_clr, mTab, phy, subsTests, dstStruc = "Age_groups", 
                          plotType = "box", plot_tukey = F, xAngle = 30, adjustAlphas = adjustAlphas,
                          chosenSamps = chosenSamps, singleSubSamp = singleSubSamp, 
                          facetVar = facetVar, facetRows = facetRows)
    
  } else if (tl == "contVar") {
    print(num_sigs)
    subsampling.plots.box(variable,"contVar",variable, 
                          c("Div.Shannon","Div.Simpson","Faiths.PD","Species_Richness","pH","BMI"), 
                          c(groupQs), only_cont,
                          gloms_clr, mTab, phy, subsTests, dstStruc = "Age_groups",
                          plotType = "box", plot_tukey = F, xAngle = 30, adjustAlphas = adjustAlphas,
                          chosenSamps = chosenSamps, singleSubSamp = singleSubSamp, 
                          facetVar = facetVar, facetRows = facetRows)
    
  }
  
}
# ************************************************************ #




# ****************************************************************************************************************** ####
# Correlate Age with diversity and taxa ####

# ******************************************************* #
get_age_correlation <- function(samps, mTab, dist_obj, varType, variable, glomTab) {
  
  if (varType == "B-div") {
    # dif.dis <- cbind(as.numeric(sapply(samps, function(x) abs(mTab[x, "Age"]-mTab[samps[samps!=x], "Age"] ))),
    #                  as.numeric(sapply(samps, function(x) as.numeric(as.matrix(dist_obj)[x, samps[samps!=x] ])) ))
    
    dif.dis <- cbind(unlist(sapply(1:(length(samps)-1), function(x) abs(mTab[samps[x], "Age"]-mTab[samps[(x+1):length(samps)], "Age"] ))),
                     unlist(sapply(1:(length(samps)-1), function(x) as.numeric(as.matrix(dist_obj)[samps[x], samps[(x+1):length(samps)] ])) ))
    print(class(dif.dis))
  } else if (varType == "contVar") {
    # dif.dis <- cbind(as.numeric(sapply(samps, function(x) (mTab[x, "Age"]-mTab[samps[samps!=x], "Age"] ))),
    #                  as.numeric(sapply(samps, function(x) (mTab[x, variable]-mTab[samps[samps!=x], variable] ) ) ))
    
    dif.dis <- cbind(unlist(sapply(1:(length(samps)-1), function(x) abs(mTab[samps[x], "Age"]-mTab[samps[(x+1):length(samps)], "Age"] ))),
                     unlist(sapply(1:(length(samps)-1), function(x) abs(mTab[samps[x], variable]-mTab[samps[(x+1):length(samps)], variable ])) ))
    # dif.dis <- dif.dis[ dif.dis[,1]>0, ]
  } else if (varType %in% c("Phylum","Class","Order","Family","Genus","Species")) {
    dif.dis <- cbind(unlist(sapply(1:(length(samps)-1), function(x) abs(mTab[samps[x], "Age"]-mTab[samps[(x+1):length(samps)], "Age"] ))),
                     unlist(sapply(1:(length(samps)-1), function(x) abs(glomTab[[ varType ]][variable, samps[x]] -
                                                                          glomTab[[ varType ]][variable, samps[(x+1):length(samps)] ])) ))
  }
  
  colnames(dif.dis) <- c("Age Difference", variable)
  
  ct <- cor.test(dif.dis[,"Age Difference"], dif.dis[,variable])
  
  return(list("dif.dis"=dif.dis, "cor"=ct))
}
# ******************************************************* #

dds <- get_age_correlation(TAS.samps, mTab.ageGroups, ordObj.agTAS$Aitchison, "B-div", "Aitchison")
dds <- get_age_correlation(TAS.samps, mTab.ageGroups, ordObj.agTAS$Aitchison, "contVar", "Div.Shannon")
dds <- get_age_correlation(TAS.samps, mTab.ageGroups, ordObj.agTAS$Aitchison, "contVar", "Faiths.PD")
dds <- get_age_correlation(TAS.samps, mTab.ageGroups, ordObj.agTAS$Aitchison, "contVar", "Species_Richness")
dds <- get_age_correlation(TAS.samps, mTab.ageGroups, ordObj.agTAS$Aitchison, "Genus", "Prevotella", glomTab = gloms_clr)
dds <- get_age_correlation(TAS.samps, mTab.ageGroups, ordObj.agTAS$Aitchison, "Genus", "Streptococcus", glomTab = gloms_clr)
dds <- get_age_correlation(TAS.samps, mTab.ageGroups, ordObj.agTAS$Aitchison, "Phylum", "Synergistetes", glomTab = gloms_clr)

plot(dds$dif.dis, main=sprintf("cor = %s, p = %s", round(dds$cor$estimate, 3), round(dds$cor$p.value, 3)))



# ***** TODO: ADAPT *******
# combine into a single table
dd33 <- as.data.frame( rbind(dif.dis.33, dif.dis.DS) )
dd33$samps <- c( rep("non-DS", nrow(dif.dis.33)), rep("DS", nrow(dif.dis.DS)) )

ggplot(dd33, aes(x=age_diffs, y=dists)) +
  geom_point(aes(shape = samps, color=samps), size=2.5) +
  scale_shape_manual(values=c(16,1), guide=F) +
  scale_color_manual(values=c("#CB6767","#73BDD3")) +
  labs(x="Age Difference", y="Aitchison Distance", color="") +
  new_scale("color") +
  geom_smooth(method = lm, aes(color=samps)) +
  scale_color_manual(values=c("#bf4342","#2f6690"), guide=F) +
  theme_classic() +
  theme(axis.title = element_text(size=17), axis.text = element_text(size=15),
        legend.title = element_text(size=17), legend.text = element_text(size=15))




# ****************************************************************************************************************** #




# ****************************************************************************************************************** ####
# ****************************************************************************************************************** ####

# ANOSIM for a list of disorders against each other ####

# plot differences in average distances for group members vs non group members

# ************************************************************************* #
get_anosim_plot_vals <- function(sam_tab, dist_obj, dist_name, comp.method, cols_to_check, var_groups,
                                 plotType="bar", anoList=NULL, adoList=NULL) {
  
  all_unit_labs <- list()
  unit_labs_pvals <- list()
  
  # ****************************** #
  for (col in cols_to_check) {
    
    if (var_groups == "Disorders") {
      # get vector of dists for each unit - members of unit against same unit
      disSamps         <- rownames(sam_tab[ sam_tab[ , col ] == "Yes", ] )
      other.disSamps   <- rownames(sam_tab[ sam_tab[ , "Chronic_disorder" ] == "Yes" & 
                                              sam_tab[ , col ] == "No", ] )
      all_non.disSamps <- rownames(sam_tab[ sam_tab[ , col ] == "No", ] )
      healthSamps      <- rownames(sam_tab[ sam_tab[ , "Chronic_disorder" ] == "No" &
                                              sam_tab[ , col ] == "No", ] )
      
      all_unit_labs[[ sprintf("%s.within", col) ]]    <- as.numeric( as.dist( dist_obj[ disSamps, disSamps ] ))
      all_unit_labs[[ sprintf("%s.other_dis", col) ]] <- as.numeric( as.matrix(dist_obj[ disSamps, other.disSamps ] ))
      all_unit_labs[[ sprintf("%s.non_dis", col) ]]   <- as.numeric( as.matrix(dist_obj[ disSamps, all_non.disSamps ] ))
      all_unit_labs[[ sprintf("%s.healthy", col) ]]   <- as.numeric( as.matrix(dist_obj[ disSamps, healthSamps ] ))
      
    } else if (var_groups == "Age_groups") {
      # get vector of dists for each unit - members of unit against same unit
      ageSamps       <- rownames(sam_tab[ sam_tab[ , col ] == "Yes", ] )
      other.ageSamps <- rownames(sam_tab[ sam_tab[ , col ] == "No", ] )
      
      all_unit_labs[[ sprintf("%s.within", col) ]]    <- as.numeric( as.dist( dist_obj[ ageSamps, ageSamps ] ))
      all_unit_labs[[ sprintf("%s.other_age", col) ]] <- as.numeric( as.matrix(dist_obj[ ageSamps, other.ageSamps ] ))
      
    } else if (var_groups == "Binary_variables") {
      
      # get all samples involved, if subsamplings were performed, get those samples, plus the mean pvals from anosim
      
    }
    
    
    if (comp.method == "ANOSIM") {
      unit_labs_pvals[[ sprintf("%s.pval", col) ]] <- anoList[[ dist_name ]][[ col ]]$signif
      unit_labs_pvals[[ sprintf("%s.R2", col) ]]   <- anoList[[ dist_name ]][[ col ]]$statistic
      
    } else if (comp.method == "PERMANOVA") {
      unit_labs_pvals[[ sprintf("%s.pval", col) ]] <- adoList[[ dist_name ]][[ col ]]$aov.tab$`Pr(>F)`[ 1 ]
      unit_labs_pvals[[ sprintf("%s.R2", col) ]]   <- sqrt(adoList[[ dist_name ]][[ col ]]$aov.tab$R2[ 1 ])
      unit_labs_pvals[[ sprintf("%s.F", col) ]]   <- adoList[[ dist_name ]][[ col ]]$aov.tab$F.Model[ 1 ]
    }
    
  }
  # ****************************** #
  # print(str(all_unit_labs))
  # print(str(unit_labs_pvals))
  
  
  if (var_groups == "Disorders") {
    plot_labs <- c("Same Disorder","Different Disorder","Healthy","All other samples")
    unit_reps <- 4
    plot_xLab <- sprintf("%s amongst disorders", comp.method)
    
  } else if (var_groups == "Age_groups") {
    plot_labs <- c("Same Age Group","Other Age Groups")
    unit_reps <- 2
    plot_xLab <- sprintf("%s amongst age groups", comp.method)
  }
  
  # ****************************** #
  if (plotType == "bar") {
    with.bet <- data.frame("vals" = unlist(lapply(all_unit_labs, mean)),
                           "labs" = rep(plot_labs, length(cols_to_check)),
                           "sds"  = unlist(lapply(all_unit_labs, sd)),
                           "unit" = rep(cols_to_check, each = unit_reps))
    
    # ****************************** #
  } else if (plotType == "box") {
    wb.labs <- sapply(names(all_unit_labs), function(x) {
      if (endsWith(x, "within")) rep(plot_labs[1], length(all_unit_labs[[ x ]]))
      else if (endsWith(x, "other_dis") | endsWith(x, "other_age")) rep(plot_labs[2], length(all_unit_labs[[ x ]]))
      else if (endsWith(x, "non_dis")) rep(plot_labs[3], length(all_unit_labs[[ x ]]))
      else if (endsWith(x, "healthy")) rep(plot_labs[4], length(all_unit_labs[[ x ]]))
    })
    
    wb.units <- sapply(cols_to_check, function(x) {
      aul.names <- names(all_unit_labs)[ startsWith(names(all_unit_labs), x) ]
      rep(x, sum(sapply(aul.names, function(y) length(all_unit_labs[[ y ]]))))
    })
    
    with.bet <- data.frame("vals" = unlist(all_unit_labs),
                           "labs" = unlist(wb.labs),
                           "unit" = unlist(wb.units))
  }
  
  with.bet$labs <- factor(with.bet$labs, levels = rev(plot_labs))
  with.bet$unit <- factor(with.bet$unit, levels = cols_to_check)
  # ****************************** #
  
  # print(unit_labs_pvals)
  star_labels <- sapply(c(cols_to_check), function(x) {
    if (unit_labs_pvals[[ sprintf("%s.pval", x) ]] < 0.001) {
      sprintf("R=%s **", round(unit_labs_pvals[[ sprintf("%s.R2", x) ]], digits = 3))
      # sprintf("F=%s**", round(unit_labs_pvals[[ sprintf("%s.F", x) ]], digits = 3))
    } else if (unit_labs_pvals[[ sprintf("%s.pval", x) ]] < 0.05) {
      sprintf("R=%s *", round(unit_labs_pvals[[ sprintf("%s.R2", x) ]], digits = 3))
      # sprintf("F=%s*", round(unit_labs_pvals[[ sprintf("%s.F", x) ]], digits = 3))
    } else {
      # sprintf("R=%s", round(unit_labs_pvals[[ sprintf("%s.R2", x) ]], digits = 3))
      ""
    }
  })
  # print(star_labels)
  # ****************************** #
  
  return(list("with.bet"=with.bet, "star_labels"=star_labels, "plot_xLab"=plot_xLab))
}
# ************************************************************************* #
plot_anosim_groups <- function(sam_tab, dist_obj, dist_name, comp.method, cols_to_check, var_groups,
                               plotType="bar", ast=NULL, anoList=NULL, adoList=NULL, whichGetFunc="original") {
  
  
  if (whichGetFunc == "original") {
    lab_vals <- get_anosim_plot_vals(sam_tab, dist_obj, dist_name, comp.method, cols_to_check, var_groups,
                                     plotType=plotType, anoList=anoList, adoList=adadoList)
  } else if (whichGetFunc == "subSamps") {
    lab_vals <- get_anosim_subSamps_plot_vals(ast, sam_tab, dist_obj, dist_name, comp.method, cols_to_check, var_groups,
                                              plotType=plotType, anoList=anoList, adoList=adadoList)
  }
  
  
  if (plotType == "bar") {
    
    if (whichGetFunc == "original") {
      annot.x <- 1:length(c(cols_to_check))
      annot.y <- lab_vals$with.bet$vals[lab_vals$with.bet$labs %in% c("Healthy","Other Age Groups","Other Age Bins")]+0.05
    } else if (whichGetFunc == "subSamps") {
      annot.x <- 1:length(c(cols_to_check))#+0.1
      annot.y <- lab_vals$with.bet$vals[lab_vals$with.bet$labs %in% c("No","No/M","No/M/One","Other Age Groups","Other Age Bins")]+0.05
    }
    
    ggplot(lab_vals$with.bet, aes(x=unit, y=vals, fill=labs)) +
      geom_bar(stat="identity", position = "dodge2") +
      geom_errorbar(aes(ymin=vals-sds, ymax=vals+sds), position = position_dodge2(width = 0.2, padding = 0.8)) +
      guides(fill=guide_legend(title=NULL, reverse = T)) +
      coord_flip() +
      annotate("text",
               x = annot.x,
               y = annot.y,
               label = (lab_vals$star_labels), size=5,
               col="red") +
      theme(axis.text = element_text(size=15), axis.title = element_text(size=17), 
            legend.text = element_text(size=15)) +#, legend.position = "bottom") +
      xlab(lab_vals$plot_xLab) + ylab(dist_name)
    
  } else if (plotType == "box") {
    
    if (whichGetFunc == "original") {
      annot.x <- 1:length(c(cols_to_check))
      annot.y <- max(lab_vals$with.bet$vals[ startsWith(as.character(lab_vals$with.bet$labs), "Same") ])
    } else if (whichGetFunc == "subSamps") {
      annot.x <- 1:length(c(cols_to_check))#+0.1
      annot.y <- max(lab_vals$with.bet$vals[ as.character(lab_vals$with.bet$labs) %in% c("No","No/M","No/M/One",
                                                                                         "Other Age Groups","Other Age Bins") ]) -
        ifelse(dist_name == "Aitchison", 5, 0.05)
    }
    
    ggplot(lab_vals$with.bet, aes(x=unit, y=vals, fill=labs)) +
      geom_boxplot(notch = T) +
      guides(fill=guide_legend(title=NULL, reverse = T)) +
      coord_flip() +
      annotate("text",
               x = annot.x,
               # y=max(lab_vals$with.bet$vals[lab_vals$with.bet$labs %in% c("Healthy","Other Age Groups")]),
               y = annot.y,
               label = (lab_vals$star_labels), size=6,
               col="red") +
      theme(axis.text = element_text(size=15), axis.title = element_text(size=17), 
            legend.text = element_text(size=15)) +#, legend.position = "bottom") +
      xlab(lab_vals$plot_xLab) + ylab(dist_name)
  }
  
}
# ************************************************************************* #




# first get blue bars (Yes values)
lv.Yes <- lab_vals$with.bet$labs[ lab_vals$with.bet$labs == "Yes/F" ]




# ************************************************************************* #
get_anosim_subSamps_plot_vals <- function(ast, sam_tab, dist_obj, dist_name, comp.method, cols_to_check, var_groups,
                                          plotType="bar", anoList=NULL, adoList=NULL) {
  
  all_unit_labs <- list()
  unit_labs_pvals <- list()
  
  # ****************************** #
  for (col in cols_to_check) {
    
    if (var_groups == "Disorders") {
      # get vector of dists for each unit - members of unit against same unit
      disSamps         <- rownames(sam_tab[ sam_tab[ , col ] == "Yes", ] )
      other.disSamps   <- rownames(sam_tab[ sam_tab[ , "Chronic_disorder" ] == "Yes" & 
                                              sam_tab[ , col ] == "No", ] )
      all_non.disSamps <- rownames(sam_tab[ sam_tab[ , col ] == "No", ] )
      healthSamps      <- rownames(sam_tab[ sam_tab[ , "Chronic_disorder" ] == "No" &
                                              sam_tab[ , col ] == "No", ] )
      
      all_unit_labs[[ sprintf("%s.within", col) ]]    <- as.numeric( as.dist( dist_obj[ disSamps, disSamps ] ))
      all_unit_labs[[ sprintf("%s.other_dis", col) ]] <- as.numeric( as.matrix(dist_obj[ disSamps, other.disSamps ] ))
      all_unit_labs[[ sprintf("%s.non_dis", col) ]]   <- as.numeric( as.matrix(dist_obj[ disSamps, all_non.disSamps ] ))
      all_unit_labs[[ sprintf("%s.healthy", col) ]]   <- as.numeric( as.matrix(dist_obj[ disSamps, healthSamps ] ))
      
    } else if (var_groups == "Age_groups") {
      # get vector of dists for each unit - members of unit against same unit
      
      teenSamps   <- rownames(sam_tab[ ! is.na(sam_tab[,"Age_groups"]) & sam_tab[,"Age_groups"] == "Teen", ])
      adultSamps  <- rownames(sam_tab[ ! is.na(sam_tab[,"Age_groups"]) & sam_tab[,"Age_groups"] == "Adult", ])
      seniorSamps <- rownames(sam_tab[ ! is.na(sam_tab[,"Age_groups"]) & sam_tab[,"Age_groups"] == "Senior", ])
      
      if (col == "Teens") {
        ag.samps <- teenSamps
        other.ag.samps <- c(adultSamps, seniorSamps)
        younger.samps <- NULL
        older.samps <- c(adultSamps, seniorSamps)
      } else if (col == "Adults") {
        ag.samps <- adultSamps
        other.ag.samps <- c(teenSamps, seniorSamps)
        younger.samps <- teenSamps
        older.samps <- seniorSamps
      } else if (col == "Seniors") {
        ag.samps <- seniorSamps
        other.ag.samps <- c(teenSamps, adultSamps)
        younger.samps <- c(teenSamps, adultSamps)
        older.samps <- NULL
      }
      # if ( TAS == F) {
      #   childSamps <- rownames(sam_tab[ sam_tab[,col] == "Child", ])
      #   all_unit_labs[[ sprintf("%s.Chlid", col) ]] <- as.numeric( as.matrix(dist_obj[ childSamps, c(teenSamps,adultSamps,seniorSamps) ] ))
      #   all_unit_labs[[ sprintf("%s.Teen", col) ]] <- as.numeric( as.matrix(dist_obj[ childSamps, c(teenSamps,adultSamps,seniorSamps) ] ))
      # } else {
      #   
      # }
      
      all_unit_labs[[ sprintf("%s.within", col) ]]    <- as.numeric( as.dist(dist_obj[ ag.samps, ag.samps ] ))
      all_unit_labs[[ sprintf("%s.other_age", col) ]] <- as.numeric( as.matrix(dist_obj[ ag.samps, other.ag.samps ] ))
      
      # if ( is.null(younger.samps))
      #   all_unit_labs[[ sprintf("%s.younger", col) ]]   <- 0
      # else 
      #   all_unit_labs[[ sprintf("%s.younger", col) ]] <- as.numeric( as.matrix(dist_obj[ ag.samps, younger.samps ] ))
      # 
      # if ( is.null(older.samps))
      #   all_unit_labs[[ sprintf("%s.older", col) ]]   <- 0
      # else 
      #   all_unit_labs[[ sprintf("%s.older", col) ]] <- as.numeric( as.matrix(dist_obj[ ag.samps, older.samps ] ))
      
      
      
    } else if (var_groups == "Age_bins") {
      # get vector of dists for each unit - members of unit against same unit
      
      samps.13_20 <- rownames(sam_tab[ ! is.na(sam_tab[,"Age_bins"]) & sam_tab[,"Age_bins"] == "13_20", ])
      samps.20_30 <- rownames(sam_tab[ ! is.na(sam_tab[,"Age_bins"]) & sam_tab[,"Age_bins"] == "20_30", ])
      samps.30_40 <- rownames(sam_tab[ ! is.na(sam_tab[,"Age_bins"]) & sam_tab[,"Age_bins"] == "30_40", ])
      samps.40_50 <- rownames(sam_tab[ ! is.na(sam_tab[,"Age_bins"]) & sam_tab[,"Age_bins"] == "40_50", ])
      samps.50_60 <- rownames(sam_tab[ ! is.na(sam_tab[,"Age_bins"]) & sam_tab[,"Age_bins"] == "50_60", ])
      samps.60_   <- rownames(sam_tab[ ! is.na(sam_tab[,"Age_bins"]) & sam_tab[,"Age_bins"] == "60+", ])
      
      if (col == "13_20") {
        ag.samps <- samps.13_20
        other.ag.samps <- c(samps.20_30, samps.30_40, samps.40_50, samps.50_60, samps.60_)
        younger.samps <- NULL
        older.samps <- c(samps.20_30, samps.30_40, samps.40_50, samps.50_60, samps.60_)
      } else if (col == "20_30") {
        ag.samps <- samps.20_30
        other.ag.samps <- c(samps.13_20, samps.30_40, samps.40_50, samps.50_60, samps.60_)
        younger.samps <- samps.13_20
        older.samps <- c(samps.30_40, samps.40_50, samps.50_60, samps.60_)
      } else if (col == "30_40") {
        ag.samps <- samps.30_40
        other.ag.samps <- c(samps.13_20, samps.20_30, samps.40_50, samps.50_60, samps.60_)
        younger.samps <- c(samps.13_20, samps.20_30)
        older.samps <- c(samps.40_50, samps.50_60, samps.60_)
      } else if (col == "40_50") {
        ag.samps <- samps.40_50
        other.ag.samps <- c(samps.13_20, samps.20_30, samps.30_40, samps.50_60, samps.60_)
        younger.samps <- c(samps.13_20, samps.20_30, samps.30_40)
        older.samps <- c(samps.50_60, samps.60_)
      } else if (col == "50_60") {
        ag.samps <- samps.50_60
        other.ag.samps <- c(samps.13_20, samps.20_30, samps.30_40, samps.40_50, samps.60_)
        younger.samps <- c(samps.13_20, samps.20_30, samps.30_40, samps.40_50)
        older.samps <- samps.60_
      } else if (col == "60+") {
        ag.samps <- samps.60_
        other.ag.samps <- c(samps.13_20, samps.20_30, samps.30_40, samps.40_50, samps.50_60)
        younger.samps <- c(samps.13_20, samps.20_30, samps.30_40, samps.40_50, samps.50_60)
        older.samps <- NULL
      }
      # if ( TAS == F) {
      #   childSamps <- rownames(sam_tab[ sam_tab[,col] == "Child", ])
      #   all_unit_labs[[ sprintf("%s.Chlid", col) ]] <- as.numeric( as.matrix(dist_obj[ childSamps, c(teenSamps,adultSamps,seniorSamps) ] ))
      #   all_unit_labs[[ sprintf("%s.Teen", col) ]] <- as.numeric( as.matrix(dist_obj[ childSamps, c(teenSamps,adultSamps,seniorSamps) ] ))
      # } else {
      #   
      # }
      
      all_unit_labs[[ sprintf("%s.within", col) ]]    <- as.numeric( as.dist(dist_obj[ ag.samps, ag.samps ] ))
      all_unit_labs[[ sprintf("%s.other_age", col) ]] <- as.numeric( as.matrix(dist_obj[ ag.samps, other.ag.samps ] ))
      
      # if ( is.null(younger.samps))
      #   all_unit_labs[[ sprintf("%s.younger", col) ]]   <- 0
      # else 
      #   all_unit_labs[[ sprintf("%s.younger", col) ]] <- as.numeric( as.matrix(dist_obj[ ag.samps, younger.samps ] ))
      # 
      # if ( is.null(older.samps))
      #   all_unit_labs[[ sprintf("%s.older", col) ]]   <- 0
      # else 
      #   all_unit_labs[[ sprintf("%s.older", col) ]] <- as.numeric( as.matrix(dist_obj[ ag.samps, older.samps ] ))
      
      
      
    } else if (var_groups == "Binary_variables") {
      
      # # get all samples involved, if subsamplings were performed, get those samples, plus the mean pvals from anosim
      # all_unit_labs <- lapply(ast[[ col ]], function(x) {
      #   samps <- x$samples
      #   sam_tab_sub <- sam_tab[ samps, ]
      #   YesSamps <- rownames(sam_tab_sub[ sam_tab_sub[ , col ] == "Yes", ])
      #   NoSamps  <- rownames(sam_tab_sub[ sam_tab_sub[ , col ] == "No", ])
      #   
      #   aul <- list()
      #   aul[[ sprintf("%s.Yes", col) ]]     <- as.numeric( as.dist(  as.matrix(x$ordObjects)[ YesSamps, YesSamps ] ))
      #   aul[[ sprintf("%s.No", col) ]]      <- as.numeric( as.dist(  as.matrix(x$ordObjects)[ NoSamps, NoSamps ] ))
      #   aul[[ sprintf("%s.between", col) ]] <- as.numeric( as.matrix(as.matrix(x$ordObjects)[ YesSamps, NoSamps ] ))
      #   
      #   return(aul)
      # })
      
      samps <- unique(unlist(lapply(ast[[ col ]], function(x) x$samples)))
      sam_tab_sub <- sam_tab[ samps, ]
      YesSamps <- rownames(sam_tab_sub[ sam_tab_sub[ , col ] %in% c("Yes","Two","F"), ])
      NoSamps  <- rownames(sam_tab_sub[ sam_tab_sub[ , col ] %in% c("No","One","M"), ])
      
      all_unit_labs[[ sprintf("%s.Yes", col) ]]     <- as.numeric( as.dist(  dist_obj[ YesSamps, YesSamps ] ))
      all_unit_labs[[ sprintf("%s.No", col) ]]      <- as.numeric( as.dist(  dist_obj[ NoSamps, NoSamps ] ))
      all_unit_labs[[ sprintf("%s.between", col) ]] <- as.numeric( as.matrix(dist_obj[ YesSamps, NoSamps ] ))
    }
    
    unit_labs_pvals[[ sprintf("%s.pval", col) ]] <- mean( p.adjust(unlist(lapply(ast[[ col ]], function(x) 
      x$ANOSIM$signif)), method = "fdr") )
    unit_labs_pvals[[ sprintf("%s.R2", col) ]] <- mean( unlist(lapply(ast[[ col ]], function(x) x$ANOSIM$statistic)) )
    
    
    # if (comp.method == "ANOSIM") {
    #   unit_labs_pvals[[ sprintf("%s.pval", col) ]] <- anoList[[ dist_name ]][[ col ]]$signif
    #   unit_labs_pvals[[ sprintf("%s.R2", col) ]]   <- anoList[[ dist_name ]][[ col ]]$statistic
    #   
    # } else if (comp.method == "PERMANOVA") {
    #   unit_labs_pvals[[ sprintf("%s.pval", col) ]] <- adoList[[ dist_name ]][[ col ]]$aov.tab$`Pr(>F)`[ 1 ]
    #   unit_labs_pvals[[ sprintf("%s.R2", col) ]]   <- sqrt(adoList[[ dist_name ]][[ col ]]$aov.tab$R2[ 1 ])
    #   unit_labs_pvals[[ sprintf("%s.F", col) ]]   <- adoList[[ dist_name ]][[ col ]]$aov.tab$F.Model[ 1 ]
    # }
    
  }
  # ****************************** #
  # print(str(all_unit_labs))
  # print(str(unit_labs_pvals))
  
  
  if (var_groups == "Disorders") {
    plot_labs <- c("Same Disorder","Different Disorder","Healthy","All other samples")
    unit_reps <- 4
    plot_xLab <- sprintf("%s amongst disorders", comp.method)
    
  } else if (var_groups == "Age_groups") {
    plot_labs <- c("Same Age Group","Other Age Groups")#,"Younger","Older")
    unit_reps <- 2 #4
    plot_xLab <- sprintf("%s amongst age groups", comp.method)
    
  } else if (var_groups == "Age_bins") {
    plot_labs <- c("Same Age Bin","Other Age Bins")#,"Younger","Older")
    unit_reps <- 2 #4
    plot_xLab <- sprintf("%s amongst age bins", comp.method)
    
  } else if (var_groups == "Binary_variables") {
    plot_labs <- c("Yes/F","No/M","Between")
    unit_reps <- 3
    plot_xLab <- sprintf("%s amongst variables", comp.method)
  }
  
  if (plotType == "bar") {
    with.bet <- data.frame("vals" = unlist(lapply(all_unit_labs, mean)),
                           "labs" = rep(plot_labs, length(cols_to_check)),
                           "sds"  = unlist(lapply(all_unit_labs, sd)),
                           "unit" = rep(cols_to_check, each = unit_reps))
    
  } else if (plotType == "box") {
    wb.labs <- sapply(names(all_unit_labs), function(x) {
      if (endsWith(x, "within") | endsWith(x, "Yes")) 
        rep(plot_labs[1], length(all_unit_labs[[ x ]]))
      else if (endsWith(x, "other_dis") | endsWith(x, "other_age") | endsWith(x, "No")) 
        rep(plot_labs[2], length(all_unit_labs[[ x ]]))
      else if (endsWith(x, "non_dis") | endsWith(x, "between") | endsWith(x, "younger")) 
        rep(plot_labs[3], length(all_unit_labs[[ x ]]))
      else if (endsWith(x, "healthy") | endsWith(x, "older")) 
        rep(plot_labs[4], length(all_unit_labs[[ x ]]))
    })
    
    wb.units <- sapply(cols_to_check, function(x) {
      aul.names <- names(all_unit_labs)[ startsWith(names(all_unit_labs), x) ]
      rep(x, sum(sapply(aul.names, function(y) length(all_unit_labs[[ y ]]))))
    })
    
    with.bet <- data.frame("vals" = unlist(all_unit_labs),
                           "labs" = unlist(wb.labs),
                           "unit" = as.character(unlist(wb.units)))
  }
  
  with.bet$labs <- factor(with.bet$labs, levels = rev(plot_labs))
  with.bet$unit <- factor(with.bet$unit, levels = cols_to_check)
  
  # ****************************** #
  
  # print(unit_labs_pvals)
  star_labels <- sapply(c(cols_to_check), function(x) {
    if (unit_labs_pvals[[ sprintf("%s.pval", x) ]] < 0.001) {
      sprintf("R = %s (**)", round(unit_labs_pvals[[ sprintf("%s.R2", x) ]], digits = 3))
      # sprintf("F=%s**", round(unit_labs_pvals[[ sprintf("%s.F", x) ]], digits = 3))
    } else if (unit_labs_pvals[[ sprintf("%s.pval", x) ]] < 0.05) {
      sprintf("R = %s (*)", round(unit_labs_pvals[[ sprintf("%s.R2", x) ]], digits = 3))
      # sprintf("F=%s*", round(unit_labs_pvals[[ sprintf("%s.F", x) ]], digits = 3))
    } else {
      sprintf("R = %s ( )", round(unit_labs_pvals[[ sprintf("%s.R2", x) ]], digits = 3))
      # ""
    }
  })
  # print(star_labels)
  # ****************************** #
  
  return(list("with.bet"=with.bet, "star_labels"=star_labels, "plot_xLab"=plot_xLab))
}

# ************************************************************************* #
# ************************************************************************* #
get_anosim_subSamps_plot_vals.alt <- function(ast, sam_tab, dist_obj, col_to_check, var_groups, plotType="box") {
  
  all_unit_labs <- list()
  
  # ****************************** #
  
  if (var_groups == "Binary_variables") {
    
    samps <- unique(unlist(lapply(ast[[ col_to_check ]], function(x) x$samples)))
    sam_tab_sub <- sam_tab[ samps, ]
    YesSamps <- rownames(sam_tab_sub[ sam_tab_sub[ , col_to_check ] %in% c("Yes","Two","F"), ])
    NoSamps  <- rownames(sam_tab_sub[ sam_tab_sub[ , col_to_check ] %in% c("No","One","M"), ])
    
    all_unit_labs[[ "Yes" ]]     <- as.numeric( as.dist(  dist_obj[ YesSamps, YesSamps ] ))
    all_unit_labs[[ "No" ]]      <- as.numeric( as.dist(  dist_obj[ NoSamps, NoSamps ] ))
    all_unit_labs[[ "Between" ]] <- as.numeric( as.matrix(dist_obj[ YesSamps, NoSamps ] ))
  }
  
  # ****************************** #
  # print(str(all_unit_labs))
  # print(str(unit_labs_pvals))
  
  
  if (var_groups == "Binary_variables") {
    plot_labs <- c("Yes/F","No/M","Between")
    unit_reps <- 3
    # plot_xLab <- sprintf("%s amongst variables", comp.method)
  }
  
  if (plotType == "box") {
    
    with.bet <- list("Yes/F"   = all_unit_labs$Yes,
                     "No/M"    = all_unit_labs$No,
                     "Between" = all_unit_labs$Between)
  }
  
  # ****************************** #
  
  return(with.bet)
}

# ************************************************************************* #


plot_anosim_groups(meta.healthy, as.matrix(ords.healthy$Aitchison), "Aitchison", "ANOSIM",
                   rev(c("Smoker","Braces.binary")), "Binary_variables", plotType="bar", whichGetFunc = "subSamps", ast=anosim_subSamps)





library(ggpubr)
# ************************************************************************* #
plot_within_dists <- function(ast, sam_tab, dist_obj, dist_name, cols_to_check, comps, print_tukey=F, plot_tukey=F) {
  
  within_dists <- list()
  unit_labs_pvals <- list()
  
  # ****************************** #
  for (col in cols_to_check) {
    
    samps <- unique(unlist(lapply(ast[[ col ]], function(x) x$samples)))
    sam_tab_sub <- sam_tab[ samps, ]
    YesSamps <- rownames(sam_tab_sub[ sam_tab_sub[ , col ] %in% c("Yes","Two","F"), ])
    NoSamps  <- rownames(sam_tab_sub[ sam_tab_sub[ , col ] %in% c("No","One","M"), ])
    
    if (comps == "Yes")
      within_dists[[ col ]] <- as.numeric( as.dist(  dist_obj[ YesSamps, YesSamps ] ))
    else if (comps == "No")
      within_dists[[ col ]] <- as.numeric( as.dist(  dist_obj[ NoSamps, NoSamps ] ))
    else if (comps == "Between")
      within_dists[[ col ]] <- as.numeric( as.matrix(dist_obj[ YesSamps, NoSamps ] ))
    
    
    unit_labs_pvals[[ sprintf("%s.pval", col) ]] <- mean( p.adjust(unlist(lapply(ast[[ col ]], function(x) 
      x$ANOSIM$signif)), method = "fdr") )
    unit_labs_pvals[[ sprintf("%s.R2", col) ]] <- mean( unlist(lapply(ast[[ col ]], function(x) x$ANOSIM$statistic)) )
    
  }
  # ****************************** #
  
  with.bet <- data.frame("vals"  = unlist(within_dists),
                         "unit"  = as.character(unlist(sapply(cols_to_check, function(x) rep(x, length(within_dists[[ x ]]))))),
                         "ano.p" = as.character(unlist(sapply(cols_to_check, function(x) 
                           rep(ifelse(unit_labs_pvals[[ sprintf("%s.pval", x) ]] > 0.05, "Not Signif", "Signif"), 
                               length(within_dists[[ x ]]))))))
  
  with.bet$unit <- factor(with.bet$unit, levels = cols_to_check)
  
  
  
  # ************************************************************* #
  if (plot_tukey == T | print_tukey == T) {
    res.aov <- aov(formula = as.numeric(as.matrix(with.bet$vals)) ~ as.factor(as.matrix(with.bet$unit)), data = with.bet)
    if (class(TukeyHSD(res.aov)[[1]][ TukeyHSD(res.aov)[[1]][ , "p adj"] < 0.05, ]) == "numeric") {
      # ******************** #
      # if only one row has signif dif, will need to print rowname separately
      rn <- rownames(TukeyHSD(res.aov)[[1]])[ TukeyHSD(res.aov)[[1]][ , "p adj"] < 0.05 ]
      # prepare table for plotting pvals with stat_pvalue_manual
      tukPval <- TukeyHSD(res.aov)[[1]][ TukeyHSD(res.aov)[[1]][ , "p adj"] < 0.05, "p adj"]
      cont.max <- max(with.bet$vals, na.rm = T)
      # arbitrarily make the space between lines based on the range of values
      rango <- (max(with.bet$vals, na.rm = T)-min(with.bet$vals, na.rm = T)) / 10
      cont.max <- cont.max + rango
      stat.test <- data.frame(".y."="cont", 
                              "group1"=strsplit(rn, "-")[[1]][1], 
                              "group2"=strsplit(rn, "-")[[1]][2],
                              "p"=NA,
                              "p.adj"=formatC(tukPval, format="e", digits=3),
                              "p.format"=formatC(tukPval, format="e", digits=3),
                              "p.signif"=symnum(tukPval, 
                                                cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                                symbols = c("****", "***", "**", "*", "ns")),
                              "method"="Tukey",
                              "y.position"=cont.max)
      if (print_tukey==T)
        print(c( rn, TukeyHSD(res.aov)[[1]][ TukeyHSD(res.aov)[[1]][ , "p adj"] < 0.05, ] ))
      
    } else {
      # ******************** #
      # prepare table for plotting pvals with stat_pvalue_manual
      TukeyTab <- TukeyHSD(res.aov)[[1]][ TukeyHSD(res.aov)[[1]][ , "p adj"] < 0.05, ]
      
      # # ******************** #
      # # for community and province, there are sometimes very strong differences between many regions
      # # so here will remove the pvalue comparison lines in plot when p=0 and its so obviously different
      # # or when there are many such lines and will have to filter some less signif ones
      # if (group_col %in% c("Community","Province") & nrow(TukeyTab[TukeyTab[,"p adj"]==0, ])>7) {
      #   TukeyTab <- TukeyTab[ TukeyTab[ , "p adj"] > 0, ]
      #   removed0s <- " (removed comparisons where p=0 to save space)"
      #   # also filter out some minimally significant lines if still many
      #   if (nrow(TukeyTab) > 10) {
      #     TukeyTab <- TukeyTab[ TukeyTab[ , "p adj"] < 0.005, ]
      #   }
      #   
      # } else if (nrow(TukeyTab) > 12) {
      if (nrow(TukeyTab) > 5) {
        # if many lines, but not because so many 0s, simply filter out some that are less significant
        TukeyTab <- TukeyTab[ 1e-15 < TukeyTab[ , "p adj"] & TukeyTab[ , "p adj"] < 0.005, ]
        # TukeyTab <- as.data.frame(TukeyTab)[ TukeyTab[ , "p adj"] > 0.005, ]
        print(TukeyTab)
      }
      
      # ******************** #
      # only works if there are sig diffs in Tukey
      if (nrow(TukeyTab) > 0) {
        group1s <- sapply(rownames(TukeyTab), function(x) strsplit(x, "-")[[1]][1])
        group2s <- sapply(rownames(TukeyTab), function(x) strsplit(x, "-")[[1]][2])
        cont.max <- max(with.bet$vals, na.rm = T)
        # arbitrarily make the space between lines based on the range of values
        rango <- (max(with.bet$vals, na.rm = T)-min(with.bet$vals, na.rm = T)) / 10
        # cont.max will determine the y.position for signif lines
        cont.max <- cont.max + rango
        cont.max <- seq( cont.max, cont.max+((nrow(TukeyTab)-1)*rango), by=rango)
        stat.test <- data.frame(".y."=rep("cont", nrow(TukeyTab)), 
                                "group1"=group1s, "group2"=group2s,
                                "p"=rep(NA, nrow(TukeyTab)), 
                                "p.adj"=formatC(TukeyTab[,"p adj"], format="e", digits=3),
                                "p.format"=formatC(TukeyTab[,"p adj"], format="e", digits=3),
                                "p.signif"=symnum(TukeyTab[,"p adj"], 
                                                  cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                                  symbols = c("****", "***", "**", "*", "ns")),
                                "method"=rep("Tukey", nrow(TukeyTab)),
                                "y.position"=cont.max)
        # print(stat.test)
        # print(str(stat.test))
        if (print_tukey==T)
          print(TukeyTab)
      } else {
        plot_tukey <- F
      }
      
    }
  }
  
  # ************************************************************* #
  
  p.x <- ifelse(dist_name=="Aitchison", 50, 
                ifelse(dist_name=="Weighted_Unifrac", 0.4,
                       ifelse(dist_name=="Unweighted_Unifrac", 0.6,
                              ifelse(dist_name %in% c("Bray","Jaccard"), 0.8, 0.5))))
  # ************************************************************* #
  boxes <- ggplot(with.bet, aes(x=unit, y=vals, fill=ano.p)) +
    geom_boxplot(notch = T) +
    guides(fill=guide_legend(title=NULL, reverse = T)) +
    # coord_flip() +
    # annotate("text",
    #          x = annot.x,
    #          # y=max(lab_vals$with.bet$vals[lab_vals$with.bet$labs %in% c("Healthy","Other Age Groups")]),
    #          y = annot.y,
    #          label = (lab_vals$star_labels), size=6,
    #          col="red") +
    # stat_compare_means(aes(group=unit, label=paste0("p = ",..p.format..)), color="red", hide.ns=T, method="kruskal",
    #                    label.x = c(3, 2, 1),
    #                    label.y =rep(p.x, length(cols_to_check))) +
    theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15, angle = 20),
          axis.title = element_text(size=17), 
          legend.text = element_text(size=15)) +#, legend.position = "bottom") +
    xlab(sprintf("\"%s\" samples for given variables", comps)) + ylab(dist_name)
  
  if (plot_tukey==T) {
    boxes <- boxes + stat_pvalue_manual(stat.test, label = "p = {p.adj}", inherit.aes=F)
  }
  
  print(boxes)
  
}
# ************************************************************************* #

plot_within_dists(anosim_subSamps, SLL2.meta, as.matrix(aitch), "Aitchison",
                   rev(names(anosim_subSamps)[1:5]), "Yes")

plot_within_dists(anosim_subSamps, SLL2.meta, as.matrix(aitch), "Aitchison",
                  rev(names(anosim_subSamps)[1:5]), "No", plot_tukey = T)

plot_within_dists(anosim_subSamps, SLL2.meta, as.matrix(aitch), "Aitchison",
                  rev(names(anosim_subSamps)[1:5]), "Between", plot_tukey = T)



library(ggpubr)
# ************************************************************************* #
get_YesVsNo_Kruskals <- function(allSubs, sam_tab, dist_obj, cols_to_check) {
  
  all_kruskals <- list()
  col_pvals    <- list()
  col_stats    <- list()
  
  # ****************************** #
  for (col in cols_to_check) {
    
    all_kruskals[[ col ]] <- list()
    col_pvals[[ col ]]    <- list()
    col_stats[[ col ]]    <- list()
    
    # ****************************** #
    if (col == "Age_groups") {
      
      for (i in names(allSubs[[ col ]])) {
        
        samps <- allSubs[[ col ]][[ i ]]
        sam_tab_sub <- sam_tab[ samps, ]
        Child_samps   <- rownames(sam_tab_sub[ sam_tab_sub[ , col ] == "Child", ])
        Teen_samps    <- rownames(sam_tab_sub[ sam_tab_sub[ , col ] == "Teen", ])
        Adult_samps   <- rownames(sam_tab_sub[ sam_tab_sub[ , col ] == "Adult", ])
        Senior_samps  <- rownames(sam_tab_sub[ sam_tab_sub[ , col ] == "Senior", ])
        
        dists <- list()
        dists[[ "Child_dists" ]]   <- as.numeric( as.dist(  dist_obj[ Child_samps, Child_samps ] ))
        dists[[ "Teen_dists" ]]    <- as.numeric( as.dist(  dist_obj[ Teen_samps, Teen_samps ] ))
        dists[[ "Adult_dists" ]]   <- as.numeric( as.dist(  dist_obj[ Adult_samps, Adult_samps ] ))
        dists[[ "Senior_dists" ]]  <- as.numeric( as.dist(  dist_obj[ Senior_samps, Senior_samps ] ))
        # dists[[ "Bet_dists" ]] <- as.numeric( as.matrix(dist_obj[ YesSamps, NoSamps ] ))
        
        all_kruskals[[ col ]][[ i ]] <- kruskal.test(dists)
      }
      
    } else if (col == "Age_groups.TAS") {
      # ****************************** #
      # this version does not include any child samples
      for (i in names(allSubs[[ col ]])) {
        
        samps <- allSubs[[ col ]][[ i ]]
        sam_tab_sub <- sam_tab[ samps, ]
        Teen_samps    <- rownames(sam_tab_sub[ sam_tab_sub[ , "Age_groups" ] == "Teen", ])
        Adult_samps   <- rownames(sam_tab_sub[ sam_tab_sub[ , "Age_groups" ] == "Adult", ])
        Senior_samps  <- rownames(sam_tab_sub[ sam_tab_sub[ , "Age_groups" ] == "Senior", ])
        
        dists <- list()
        dists[[ "Teen_dists" ]]    <- as.numeric( as.dist(  dist_obj[ Teen_samps, Teen_samps ] ))
        dists[[ "Adult_dists" ]]   <- as.numeric( as.dist(  dist_obj[ Adult_samps, Adult_samps ] ))
        dists[[ "Senior_dists" ]]  <- as.numeric( as.dist(  dist_obj[ Senior_samps, Senior_samps ] ))
        # dists[[ "Bet_dists" ]] <- as.numeric( as.matrix(dist_obj[ YesSamps, NoSamps ] ))
        
        all_kruskals[[ col ]][[ i ]] <- kruskal.test(dists)
      }
      
    } else if (col == "Age_bins") {
      # ****************************** #
      # also does not include any child samples
      for (i in names(allSubs[[ col ]])) {
        
        samps <- allSubs[[ col ]][[ i ]]
        sam_tab_sub <- sam_tab[ samps, ]
        samps.13_20 <- rownames(sam_tab_sub[ sam_tab_sub[ , col ] == "13_20", ])
        samps.20_30 <- rownames(sam_tab_sub[ sam_tab_sub[ , col ] == "20_30", ])
        samps.30_40 <- rownames(sam_tab_sub[ sam_tab_sub[ , col ] == "30_40", ])
        samps.40_50 <- rownames(sam_tab_sub[ sam_tab_sub[ , col ] == "40_50", ])
        samps.50_60 <- rownames(sam_tab_sub[ sam_tab_sub[ , col ] == "50_60", ])
        samps.60_   <- rownames(sam_tab_sub[ sam_tab_sub[ , col ] == "60+", ])
        
        dists <- list()
        dists[[ "dists.13_20" ]] <- as.numeric( as.dist(  dist_obj[ samps.13_20, samps.13_20 ] ))
        dists[[ "dists.20_30" ]] <- as.numeric( as.dist(  dist_obj[ samps.20_30, samps.20_30 ] ))
        dists[[ "dists.30_40" ]] <- as.numeric( as.dist(  dist_obj[ samps.30_40, samps.30_40 ] ))
        dists[[ "dists.40_50" ]] <- as.numeric( as.dist(  dist_obj[ samps.40_50, samps.40_50 ] ))
        dists[[ "dists.50_60" ]] <- as.numeric( as.dist(  dist_obj[ samps.50_60, samps.50_60 ] ))
        dists[[ "dists.60+" ]]   <- as.numeric( as.dist(  dist_obj[ samps.60_, samps.60_ ] ))
        # dists[[ "Bet_dists" ]] <- as.numeric( as.matrix(dist_obj[ YesSamps, NoSamps ] ))
        
        all_kruskals[[ col ]][[ i ]] <- kruskal.test(dists)
      }
      
    } else {#if (var_groups == "Binary_variables") {
      # ****************************** #
      for (i in names(allSubs[[ col ]])) {
        
        samps <- allSubs[[ col ]][[ i ]]
        sam_tab_sub <- sam_tab[ samps, ]
        YesSamps <- rownames(sam_tab_sub[ sam_tab_sub[ , col ] %in% c("Yes","Two","F"), ])
        NoSamps  <- rownames(sam_tab_sub[ sam_tab_sub[ , col ] %in% c("No","One","M"), ])
        
        dists <- list()
        dists[[ "Yes_dists" ]] <- as.numeric( as.dist(  dist_obj[ YesSamps, YesSamps ] ))
        dists[[ "No_dists" ]]  <- as.numeric( as.dist(  dist_obj[ NoSamps, NoSamps ] ))
        # dists[[ "Bet_dists" ]] <- as.numeric( as.matrix(dist_obj[ YesSamps, NoSamps ] ))
        
        all_kruskals[[ col ]][[ i ]] <- kruskal.test(dists)
      }
    }
    # ****************************** #
  }
  # ****************************** #
  return(all_kruskals)
}
# ************************************************************************* #
plot_YesVsNo_Kruskals <- function(yes_no_kruskals, dist_name, plotVal="pval") {
  
  dont_plot <- c("seqGroup","Do_you_feel_well","Lactose_intolerant","Anemia","Kidney_issues",
                 "Circulatory_issues","Lung_issues")
  vars_to_plot <- names(yes_no_kruskals)[ ! names(yes_no_kruskals) %in% dont_plot ]
  
  # ****************************** #
  # ****************************** #
  if (plotVal == "both") {
    
    # p-values
    krus.ps <- data.frame(vars = as.character(unlist(sapply(vars_to_plot, function(x) rep(x, length(yes_no_kruskals[[x]]))))),
                                pv   = unlist(lapply(yes_no_kruskals[ vars_to_plot ], 
                                                     function(v) p.adjust(unlist(lapply(v, function(n) n$p.value)), method = "fdr"))))
    # function(v) unlist(lapply(v, function(n) n$p.value)))))
    # for those p=0
    krus.ps$pv[ krus.ps$pv == 0 ] <- min(krus.ps$pv[ krus.ps$pv != 0 ])
    krus.ps$pv <- -log(krus.ps$pv)
    
    # statistic
    krus.stats <- data.frame(vars = as.character(unlist(sapply(vars_to_plot, function(x) rep(x, length(yes_no_kruskals[[x]]))))),
                                pv = unlist(lapply(yes_no_kruskals[ vars_to_plot ], function(v) unlist(lapply(v, function(n) n$statistic)))))
    
    kps <- rbind(krus.ps, krus.stats)
    kps$plot_val <- c(rep("-log(p.adj)", nrow(krus.ps)), rep("Kruskal-Wallis test statistic", nrow(krus.stats)))
    
    gg.krus <- ggplot(kps, aes(x=reorder(vars, pv, FUN=median), y=pv, fill=reorder(vars, pv, FUN=median))) +
      geom_boxplot(notch = T) +
      guides(fill=F) +#guide_legend(title=NULL, reverse = T)) +
      coord_flip() +
      ggtitle(sprintf("%s for difference between %s distances", "fdfs", dist_name)) + 
      theme(axis.text = element_text(size=15), axis.title = element_text(size=17),
            legend.text = element_text(size=15)) +#, legend.position = "bottom") +
      facet_wrap(~plot_val, scales = "free") +
      xlab("Variables") #+ ylab(yl)# + ylim(0,5000)
    
    # if (plotVal == "pval")
    #   gg.krus <- gg.krus + geom_hline(yintercept = -log(0.05), color="red")
    
    gg.krus
    
  } else {
    # ****************************** #
    # ****************************** #
    # make dataframe of pvals and statistics for each variable
    if (plotVal == "pval") {
      krus.ps.stats <- data.frame(vars = as.character(unlist(sapply(vars_to_plot, function(x) rep(x, length(yes_no_kruskals[[x]]))))),
                                  pv   = unlist(lapply(yes_no_kruskals[ vars_to_plot ], 
                                                       function(v) p.adjust(unlist(lapply(v, function(n) n$p.value)), method = "fdr"))))
      # function(v) unlist(lapply(v, function(n) n$p.value)))))
      # for those p=0
      krus.ps.stats$pv[ krus.ps.stats$pv == 0 ] <- min(krus.ps.stats$pv[ krus.ps.stats$pv != 0 ])
      krus.ps.stats$pv <- -log(krus.ps.stats$pv)
      
      yl <- "-log(p.adj)"
      
    } else if (plotVal == "stat") {
      
      krus.ps.stats <- data.frame(vars = as.character(unlist(sapply(vars_to_plot, function(x) rep(x, length(yes_no_kruskals[[x]]))))),
                                  pv = unlist(lapply(yes_no_kruskals[ vars_to_plot ], function(v) unlist(lapply(v, function(n) n$statistic)))))
      yl <- "Kruskal-Wallis test statistic"
      
    }
    
    # ****************************** #
    # print(table(krus.ps.stats$vars))
    # ggplot(krus.ps.stats, aes(x=vars, y=stat, fill=vars)) +
    gg.krus <- ggplot(krus.ps.stats, aes(x=reorder(vars, pv, FUN=median), y=pv, fill=reorder(vars, pv, FUN=median))) +
      geom_boxplot(notch = T) +
      guides(fill=F) +#guide_legend(title=NULL, reverse = T)) +
      coord_flip() +
      ggtitle(sprintf("%s for difference between %s distances", yl, dist_name)) + 
      # annotate("text",
      #          x = annot.x,
      #          # y=max(lab_vals$with.bet$vals[lab_vals$with.bet$labs %in% c("Healthy","Other Age Groups")]),
      #          y = annot.y,
      #          label = (lab_vals$star_labels), size=6,
      #          col="red") +
      theme(axis.text = element_text(size=15), axis.title = element_text(size=17),
            legend.text = element_text(size=15)) +#, legend.position = "bottom") +
      xlab("Variables") + ylab(yl)# + ylim(0,5000)
    
    if (plotVal == "pval")
      gg.krus <- gg.krus + geom_hline(yintercept = -log(0.05), color="red")
    
    gg.krus
    # ****************************** #
    # ****************************** #
  }
}
# ************************************************************************* #



# ************************************************************************* #
get_dists_per_group <- function(variable, group, mTab, dist_obj, chosenSamps=NULL, allSubs=NULL, mean_of_mean=F) {
  
  # ********* #
  if (mean_of_mean == T) {
    dists <- list()
    for (i in names(allSubs)) {
      samps <- allSubs[[ i ]][ mTab[allSubs[[ i ]], variable ] == group ]
      dists[[ i ]] <- as.numeric( as.dist(  as.matrix(dist_obj)[ samps, samps ] ))
    }
    return( dists )
    # ********* #
  } else {
    # ********* #
    if ( ! is.null(allSubs)) {
      samps <- unique(unlist(allSubs))
    } else if ( ! is.null(chosenSamps)) {
      samps <- chosenSamps
    }
    
    samps <- samps[ mTab[samps, variable] == group ]
    
    return( as.numeric( as.dist(  as.matrix(dist_obj)[ samps, samps ] )) )
    # ********* #
  }
  
}
# ************************************************************************* #


# ************************************************************************* #
library(reshape)
get_YesVsNo_anosim_plotVals <- function(varsOfInterest, stat, test) {
  
  yn_dist_diffs <- sapply(varsOfInterest, function(v) {
    lapply(all_subSamps[[ v ]], function(i) {
      # **** #
      yndd <- sapply(c("Yes/F","No/M"), function(g) {
        gVal <- ifelse(v=="Gender", strsplit(g, "/")[[1]][2], strsplit(g, "/")[[1]][1])
        if (stat == "mean") {
          mean(get_dists_per_group(v, gVal, SLL2.meta, aitch, chosenSamps = i))
        } else if (stat == "median") {
          median(get_dists_per_group(v, gVal, SLL2.meta, aitch, chosenSamps = i))
        }
        
      })
      yndd[[ "Yes/F" ]] - yndd[[ "No/M" ]]
      # **** #
    })
  })
  
  yn_dist_diffs <- as.list(as.data.frame(yn_dist_diffs))
  yn_dist_diffs <- lapply(yn_dist_diffs, unlist)
  
  # ************* #
  ynd <- melt.list(yn_dist_diffs)
  colnames(ynd) <- c("diff","Variable")
  
  if (stat == "mean") {
    ynd$KW.p.adj         <- sapply(ynd$Variable, function(x) yn_mean_diffs[x, "KW.p.adj"])
    ynd$`-log(KW.p.adj)` <- sapply(ynd$KW.p.adj, function(x) -log(x))
    ynd$KW.stat          <- sapply(ynd$Variable, function(x) yn_mean_diffs[x, "KW.stat"])
  } else if (stat == "median") {
    ynd$KW.p.adj         <- sapply(ynd$Variable, function(x) ifelse(yn_median_diffs[x, "KW.p.adj"]==0,
                                                                   min(yn_median_diffs[, "KW.p.adj"][yn_median_diffs[, "KW.p.adj"] > 0]),
                                                                   yn_median_diffs[x, "KW.p.adj"]))
    ynd$`-log(KW.p.adj)` <- sapply(ynd$KW.p.adj, function(x) -log(x))
    ynd$KW.stat          <- sapply(ynd$Variable, function(x) yn_median_diffs[x, "KW.stat"])
  }
  
  # ************* #
  if (test == "anosim") {
    ano_r_p  <- sapply(unique(ynd$Variable), function(x) {
      as <- readRDS(sprintf("%s/R_objects/anosim_subSamps.temp.%s.rds", p2_dir, x))
      as.p <- lapply(as, function(y) y$ANOSIM$signif)
      as.r <- lapply(as, function(y) y$ANOSIM$statistic)
      return(list("p"=as.p, "R"=as.r))
    })
    arp <- as.list(as.data.frame(ano_r_p))
    arp <- lapply(arp, function(x) lapply(x, unlist))
    
  } else if (test == "adonis") {
    ano_r_p  <- sapply(unique(ynd$Variable), function(x) {
      as <- readRDS(sprintf("%s/R_objects/adonis_subSamps.temp.%s.rds", p2_dir, x))
      as.p <- lapply(as, function(y) y$Adonis$Aitchison[x, "Pr(>F)"])
      as.r <- lapply(as, function(y) y$Adonis$Aitchison[x, "R2"])
      return(list("p"=as.p, "R"=as.r))
    })
    arp <- as.list(as.data.frame(ano_r_p))
    arp <- lapply(arp, function(x) lapply(x, unlist))
  }
  
  
  # ************* #
  if (stat == "mean") {
    ynd$anosim.p.adj         <- sapply(ynd$Variable, function(x) mean(p.adjust(arp[[ x ]]$p, method = "fdr")))
    ynd$`-log(anosim.p.adj)` <- sapply(ynd$anosim.p.adj, function(x) -log(x))
    ynd$anosim.R             <- sapply(ynd$Variable, function(x) mean(arp[[ x ]]$R))
  } else if (stat == "median") {
    ynd$anosim.p.adj         <- sapply(ynd$Variable, function(x) median(p.adjust(arp[[ x ]]$p, method = "fdr")))
    ynd$`-log(anosim.p.adj)` <- sapply(ynd$anosim.p.adj, function(x) -log(x))
    ynd$anosim.R             <- sapply(ynd$Variable, function(x) median(arp[[ x ]]$R))
  }
  
  # ************* #
  # update Variable name to include # of samples that are Yes/F
  ynd$Variable <- sapply(ynd$Variable, function(x) {
    if (x %in% c("Fluoride_toothpaste","Wash_hands_before_eat","Wash_hands_after_bathroom","Do_you_feel_well")) {
      mTab.sub <- meta.healthy[ all_subSamps[[ x ]]$`1`, ]
      var_num  <- nrow(mTab.sub[ ! is.na(mTab.sub[ , x]) & mTab.sub[ , x] == "No", ])
      
    } else if (x %in% c("Gender","Mouth_wounds.binary","Analgesics","Bite_nails","Chew_pens","Kissing_partner",
                      "Pets","Pets.Dogs","Pets.Cats","Pets.Both_dogs_cats")) {
      var_num <- 100
      
    } else if (x %in% c("Cystic_fibrosis","Downs_Syndrome","Celiac","Diabetes","Hypertension","Cholesterol",
                      "Depression","Anxiety","Headaches","Lactose_intolerant","Gastritis","Intestinal_issues",
                      "Anemia","Sinusitis","Fibrosis_carrier","Thyroid_issue","Hypothyroidism","Cancer","Transplant",
                      "Immune_issues","Skin_issues","Lung_issues","Circulatory_issues","Kidney_issues",
                      "Tonsil_issues","Samp_collect_issues","Birth_control","Antihistamines")) {
      mTab.sub <- SLL2.meta[ all_subSamps[[ x ]]$`1`, ]
      var_num  <- nrow(mTab.sub[ ! is.na(mTab.sub[ , x]) & mTab.sub[ , x] == "Yes", ])
      
    } else {
      mTab.sub <- meta.healthy[ all_subSamps[[ x ]]$`1`, ]
      var_num  <- nrow(mTab.sub[ ! is.na(mTab.sub[ , x]) & mTab.sub[ , x] == "Yes", ])
      
    }
    
    var_p    <- unique(ynd$anosim.p.adj[ ynd$Variable == x ])
    var_R    <- unique(ynd$anosim.R[ ynd$Variable == x ])
    
    # sprintf("%s (n=%s)\n%s p.adj = %s\n%s R = %s", x, var_num, stat, round(var_p, 3), stat, round(var_R, 3))
    sprintf("**%s** (n=%s)<br>%s p.adj = %s<br>%s R = %s", x, var_num, stat, round(var_p, 3), stat, round(var_R, 3))
    # check these links for properly using the markdown formatting to get bold text:
    #   https://stackoverflow.com/questions/32555531/how-to-italicize-part-one-or-two-words-of-an-axis-title
    #   https://github.com/wilkelab/ggtext
  })
  # ************* #
  
  return(ynd)
}
# ************************************************************************* #
plot_YesVsNo_anosim_plotVals <- function(ynd, stat, colorBy, dist_meas="Aitchison") {
  
  ynd$boxFill <- ynd[ , colorBy]
  
  ano_plot <- ggplot(ynd, aes(x=reorder(Variable, abs(diff), FUN=median), y=diff, fill=boxFill)) +
    geom_boxplot(notch = T) +
    geom_hline(yintercept = 0, color="red") +
    coord_flip() + labs(x=NULL, y=sprintf("%s %s( Yes/F - No/M )", stat, dist_meas)) +
    theme_minimal() +
    theme(axis.title = element_text(size=17), legend.title = element_text(size=12),
          # axis.text.x = element_text(size=17), axis.text.y = element_text(size=13))
          axis.text.x = element_text(size=17), axis.text.y = ggtext::element_markdown(size=13))
  
  if (colorBy == "anosim.R") {
    ano_plot +
      scale_fill_gradient2(low="darkred", mid="white", high="darkblue", midpoint=0, name=colorBy)
    # scale_fill_gradientn(colours = rev(colorspace::heat_hcl(10)), name=colorBy)
  } else {
    ano_plot +
      scale_fill_gradient2(low="white", high="darkblue", name=colorBy)
    # scale_fill_gradientn(colours = rev(colorspace::heat_hcl(10)), name=colorBy)
  }
  
}
# ************************************************************************* #





# ************************************************************************* #
library(reshape)
get_Age_anosim_plotVals <- function(varsOfInterest, stat, mTab, ano_sub, dist_obj) {
  
  age_dist_diffs <- sapply(varsOfInterest, function(v) {
    lapply(ano_sub[[ v ]], function(i) {
      # **** #
      a_d_d <- sapply(c("Yes/F","No/M"), function(g) {
        gVal <- ifelse(v=="Gender", strsplit(g, "/")[[1]][2], strsplit(g, "/")[[1]][1])
        if (stat == "mean") {
          mean(get_dists_per_group(v, gVal, mTab, dist_obj, chosenSamps = i$samples))
        } else if (stat == "median") {
          median(get_dists_per_group(v, gVal, mTab, dist_obj, chosenSamps = i$samples))
        }
        
      })
      a_d_d[[ "Yes/F" ]] - a_d_d[[ "No/M" ]]
      # **** #
    })
  })
  
  age_dist_diffs <- as.list(as.data.frame(age_dist_diffs))
  age_dist_diffs <- lapply(age_dist_diffs, unlist)
  
  # ************* #
  a_d <- melt.list(age_dist_diffs)
  colnames(a_d) <- c("diff","Variable")
  
  # ************* #
  
  ano_r_p  <- sapply(unique(a_d$Variable), function(x) {
    as.p <- lapply(ano_sub[[ x ]], function(x) x$ANOSIM$signif)
    as.r <- lapply(ano_sub[[ x ]], function(x) x$ANOSIM$statistic)
    return(list("p"=as.p, "R"=as.r))
  })
  arp <- as.list(as.data.frame(ano_r_p))
  arp <- lapply(arp, function(x) lapply(x, unlist))
  
  if (stat == "mean") {
    a_d$anosim.p.adj         <- sapply(a_d$Variable, function(x) mean(p.adjust(arp[[ x ]]$p, method = "fdr")))
    a_d$`-log(anosim.p.adj)` <- sapply(a_d$anosim.p.adj, function(x) -log(x))
    a_d$anosim.R             <- sapply(a_d$Variable, function(x) mean(arp[[ x ]]$R))
  } else if (stat == "median") {
    a_d$anosim.p.adj         <- sapply(a_d$Variable, function(x) median(p.adjust(arp[[ x ]]$p, method = "fdr")))
    a_d$`-log(anosim.p.adj)` <- sapply(a_d$anosim.p.adj, function(x) -log(x))
    a_d$anosim.R             <- sapply(a_d$Variable, function(x) median(arp[[ x ]]$R))
  }
  
  # ************* #
  # update Variable name to include # of samples that are Yes/F
  a_d$Variable <- sapply(a_d$Variable, function(x) {
    mTab.sub <- mTab[ ano_sub[[ x ]]$`1`$samples, ]
    
    var_num  <- nrow(mTab.sub[ ! is.na(mTab.sub[ , x]) & mTab.sub[ , x] == "Yes", ])
    var_p    <- unique(a_d$anosim.p.adj[ a_d$Variable == x ])
    var_R    <- unique(a_d$anosim.R[ a_d$Variable == x ])
    
    # sprintf("%s (n=%s)\n%s p.adj = %s\n%s R = %s", crayon::bold(x), var_num, stat, round(var_p, 3), stat, round(var_R, 3))
    sprintf("**%s** (n=%s)<br>%s p.adj = %s<br>%s R = %s", x, var_num, stat, round(var_p, 3), stat, round(var_R, 3))
    # check these links for properly using the markdown formatting to get bold text:
    #   https://stackoverflow.com/questions/32555531/how-to-italicize-part-one-or-two-words-of-an-axis-title
    #   https://github.com/wilkelab/ggtext
  })
  
  # reorder levels of factor
  a_d$Variable <- factor(a_d$Variable, levels=rev(unique(a_d$Variable)))
  return(a_d)
}
# ************************************************************************* #
plot_Age_anosim_plotVals <- function(a_d, stat, colorBy, dist_meas="Aitchison", paired=F) {
  
  a_d$boxFill <- a_d[ , colorBy]
  
  # ggplot(a_d, aes(x=reorder(Variable, abs(diff), FUN=median), y=diff, fill=boxFill)) +
  ano_plot <- ggplot(a_d, aes(x=Variable, y=diff, fill=boxFill)) +
    geom_boxplot(notch = T) +
    geom_hline(yintercept = 0, color="red") +
    coord_flip() + 
    labs(x=NULL, y=sprintf("%s %s( %s )", stat, dist_meas, ifelse(paired, "Younger - Older", "Yes - No"))) +
    theme_minimal() +
    theme(axis.title = element_text(size=17), legend.title = element_text(size=12),
          # axis.text.x = element_text(size=17), axis.text.y = element_text(size=13))
          axis.text.x = element_text(size=17), axis.text.y = ggtext::element_markdown(size=13))
  
  if (paired & length(unique(a_d$Variable)) > 4) {
    # for pairwise plots for Age_bins
    ano_plot <- ano_plot +
      geom_vline(xintercept = c(1.5, 3.5, 6.5, 10.5), color="darkgreen")
  }
  
  if (colorBy == "anosim.R") {
    ano_plot +
      scale_fill_gradient2(low="darkred", mid = "white", high="darkblue", midpoint=0, name=colorBy)
      # scale_fill_gradientn(colours = rev(colorspace::heat_hcl(10)), name=colorBy)
  } else {
    ano_plot +
      scale_fill_gradient2(low="white", high="darkblue", name=colorBy)
      # scale_fill_gradientn(colours = rev(colorspace::heat_hcl(10)), name=colorBy)
  }
  
}
# ************************************************************************* #




# ************************************************************************* #
library(reshape)
get_Age_pairwise_anosim_plotVals <- function(variable, stat, mTab, ano_sub, dist_obj) {
  
  age_dist_diffs <- sapply(names(ano_sub), function(agp) {
    lapply(ano_sub[[ agp ]], function(i) {
      ages <- strsplit(agp, "-")[[1]]
      a_d_d <- sapply(ages, function(g) {
        if (stat == "mean") {
          mean(get_dists_per_group(variable, g, mTab, dist_obj, chosenSamps = i$samples))
        } else if (stat == "median") {
          median(get_dists_per_group(variable, g, mTab, dist_obj, chosenSamps = i$samples))
        }
        
      })
      a_d_d[[ ages[1] ]] - a_d_d[[ ages[2] ]]
      # **** #
    })
  })
  
  
  # age_dist_diffs <- sapply(varsOfInterest, function(v) {
  #   lapply(ano_sub[[ v ]], function(i) {
  #     # **** #
  #     a_d_d <- sapply(c("Yes/F","No/M"), function(g) {
  #       gVal <- ifelse(v=="Gender", strsplit(g, "/")[[1]][2], strsplit(g, "/")[[1]][1])
  #       if (stat == "mean") {
  #         mean(get_dists_per_group(v, gVal, mTab, dist_obj, chosenSamps = i$samples))
  #       } else if (stat == "median") {
  #         median(get_dists_per_group(v, gVal, mTab, dist_obj, chosenSamps = i$samples))
  #       }
  #       
  #     })
  #     a_d_d[[ "Yes/F" ]] - a_d_d[[ "No/M" ]]
  #     # **** #
  #   })
  # })
  
  age_dist_diffs <- as.list(as.data.frame(age_dist_diffs))
  age_dist_diffs <- lapply(age_dist_diffs, unlist)
  
  # ************* #
  a_d <- melt.list(age_dist_diffs)
  colnames(a_d) <- c("diff","Variable")
  
  # ************* #
  
  ano_r_p  <- sapply(unique(a_d$Variable), function(x) {
    as.p <- lapply(ano_sub[[ x ]], function(x) x$ANOSIM$signif)
    as.r <- lapply(ano_sub[[ x ]], function(x) x$ANOSIM$statistic)
    return(list("p"=as.p, "R"=as.r))
  })
  arp <- as.list(as.data.frame(ano_r_p))
  arp <- lapply(arp, function(x) lapply(x, unlist))
  
  if (stat == "mean") {
    a_d$anosim.p.adj         <- sapply(a_d$Variable, function(x) mean(p.adjust(arp[[ x ]]$p, method = "fdr")))
    a_d$`-log(anosim.p.adj)` <- sapply(a_d$anosim.p.adj, function(x) -log(x))
    a_d$anosim.R             <- sapply(a_d$Variable, function(x) mean(arp[[ x ]]$R))
  } else if (stat == "median") {
    a_d$anosim.p.adj         <- sapply(a_d$Variable, function(x) median(p.adjust(arp[[ x ]]$p, method = "fdr")))
    a_d$`-log(anosim.p.adj)` <- sapply(a_d$anosim.p.adj, function(x) -log(x))
    a_d$anosim.R             <- sapply(a_d$Variable, function(x) median(arp[[ x ]]$R))
  }
  
  # ************* #
  # update Variable name to include # of samples that are Yes/F
  a_d$Variable <- sapply(a_d$Variable, function(agp) {
    mTab.sub <- mTab[ ano_sub[[ agp ]]$`1`$samples, ]
    
    ages <- strsplit(agp, "-")[[1]]
    
    var_num1  <- nrow(mTab.sub[ ! is.na(mTab.sub[ , variable]) & mTab.sub[ , variable] == ages[1], ])
    var_num2  <- nrow(mTab.sub[ ! is.na(mTab.sub[ , variable]) & mTab.sub[ , variable] == ages[2], ])
    var_p    <- unique(a_d$anosim.p.adj[ a_d$Variable == agp ])
    var_R    <- unique(a_d$anosim.R[ a_d$Variable == agp ])
    
    # sprintf("%s (n=%s)\n%s p.adj = %s\n%s R = %s", crayon::bold(x), var_num, stat, round(var_p, 3), stat, round(var_R, 3))
    sprintf("**%s** (n = %s, %s)<br>%s p.adj = %s<br>%s R = %s", gsub("-"," - ",agp), var_num1, var_num2, 
            stat, round(var_p, 3), 
            stat, round(var_R, 3))
    # check these links for properly using the markdown formatting to get bold text:
    #   https://stackoverflow.com/questions/32555531/how-to-italicize-part-one-or-two-words-of-an-axis-title
    #   https://github.com/wilkelab/ggtext
  })
  
  # reorder levels of factor
  a_d$Variable <- factor(a_d$Variable, levels=rev(unique(a_d$Variable)))
  return(a_d)
}
# ************************************************************************* #







# ****************************************************************************************************************** ####
# Networks for subSamps ####
library(nettools)
library(SpiecEasi)
library(seqtime)

# ************************************************************************* #
get_networks_for_subSamps <- function(variable, tagl, allSubList, mTab, glomTab, minCount=15, minOcc=20) {
  
  # how many are "Yes" for given variable, to know how many to match ("Two" for seqGroup)
  main_vals <- c("Yes","Two","F")
  subs_vals <- c("No","One","M")
  # have to use the "No" as main group here since there are many more "Yes" samples
  if (variable %in% c("Fluoride_toothpaste","Wash_hands_before_eat","Wash_hands_after_bathroom","Do_you_feel_well",
                      "MALDI.Yeast_detected","Full_MALDI.Candida")) {
    main_vals <- c("No","One","M")
    subs_vals <- c("Yes","Two","F")
  }
  
  
  # first get objects for main value
  samps.all <- unique(unlist(allSubList[[ variable ]]))
  samps.main <- samps.all[ mTab[samps.all, variable] %in% main_vals ]
  # then for all of the matched controls together
  samps.all_subs <- samps.all[ mTab[samps.all, variable] %in% subs_vals ]
  
  # ************************** #
  # from this tutorial: http://psbweb05.psb.ugent.be/conet/microbialnetworks/spieceasi.php
  # filter out taxa that dont have at least 15 reads in at least 20 of the CF + mC samples, but with this function,
  #    keeps a row of all the counts that were removed, so it can keep the total counts for each sample
  g10 <- as.data.frame(glomTab[ , c(samps.main, samps.all_subs) ])
  g10[ g10 <= minCount ] <- 0
  # filterobj <- filterTaxonMatrix(glomTab[ , c(samps.main, samps.all_subs) ], minocc=100, keepSum = TRUE, return.filtered.indices = TRUE)
  filterobj <- filterTaxonMatrix(g10, minocc=minOcc, keepSum = TRUE, return.filtered.indices = TRUE)
  otus.f <- filterobj$mat
  
  # replace the values that were made to be 0s for the filtering
  anti.g10 <- as.data.frame(glomTab[ , c(samps.main, samps.all_subs) ])
  anti.g10[ anti.g10 > minCount ] <- 0
  # replace them into g10
  g10[ anti.g10 != 0 ] <- anti.g10[ anti.g10 != 0 ]
  # first update the values that were made 0
  upd.otus_f <- otus.f
  upd.otus_f[ 1:(nrow(upd.otus_f)-1), ] <- as.data.frame(glomTab[ rownames(upd.otus_f)[ 1:(nrow(upd.otus_f)-1) ], c(samps.main, samps.all_subs) ])
  # then update the summed row from the removed values
  upd.otus_f[ nrow(upd.otus_f), ] <- as.numeric(upd.otus_f[ nrow(upd.otus_f), ] + (colSums(upd.otus_f) - colSums(otus.f)))
  
  # finally update all the names appropriately, including pseudo-name for the summed counts, which will be ignored later in the plots
  taxa.f <- tax_table(tagl)[setdiff(1:nrow(tax_table(tagl)),filterobj$filtered.indices),]
  dummyTaxonomy <- c("k__dummy","p__","c__","o__","f__","g__")
  taxa.f <- rbind(taxa.f, dummyTaxonomy)
  rownames(taxa.f)[nrow(taxa.f)] <- "0"
  rownames(upd.otus_f)[nrow(upd.otus_f)] <- "0"
  
  tagl.filtered <- phyloseq(otu_table(upd.otus_f, taxa_are_rows = T),
                            sample_data(tagl)[ c(samps.main, samps.all_subs), ],
                            tax_table(taxa.f))
  print("tagl.filtered:")
  print(tagl.filtered)
  # ************************** #
  
  # first get network objects for main value, which does not change across subSamps
  # tagl.main  <- subset_samples(tagl.filtered, sample_names(tagl.filtered) %in% samps.main)
  tagl.main  <- prune_samples(samps.main, tagl.filtered)
  se.mb.main <- spiec.easi(tagl.main, method="mb", lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=50))
  igraphs.main.allNodes <- adj2igraph(getRefit(se.mb.main),  vertex.attr=list(name=taxa_names(tagl.main)), rmEmptyNodes = F)
  igraphs.main.noEmpty <- adj2igraph(getRefit(se.mb.main),  vertex.attr=list(name=taxa_names(tagl.main)), rmEmptyNodes = T)
  # plot_network(igraphs.CF.Yes.noEmpty, tagl.CF.Yes, type='taxa', color="Phylum", title = "Co-occurrence network: CF samples")
  
  
  # then for all of the matched controls together (also unchanged since its all of them together)
  # tagl.all_subs  <- subset_samples(tagl.filtered, sample_names(tagl.filtered) %in% samps.all_subs)
  tagl.all_subs  <- prune_samples(samps.all_subs, tagl.filtered)
  se.mb.all_subs <- spiec.easi(tagl.all_subs, method="mb", lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=50))
  igraphs.all_subs.allNodes <- adj2igraph(getRefit(se.mb.all_subs),vertex.attr=list(name=taxa_names(tagl.all_subs)), rmEmptyNodes = F)
  igraphs.all_subs.noEmpty <- adj2igraph(getRefit(se.mb.all_subs),vertex.attr=list(name=taxa_names(tagl.all_subs)), rmEmptyNodes = T)
  
  
  
  # ************************** #
  # then get objects for each of the individual subSamp networks
  
  main_v_subs.hammings <- list()
  SE_objects <- list()
  igraphs.allNodes <- list()
  igraphs.noEmpty <- list()
  
  for (i in names(allSubList[[ variable ]])) {
    
    print(i)
    
    # get samples
    samps <- allSubList[[ variable ]][[ i ]]
    samps.subs <- samps[ mTab[samps, variable] %in% subs_vals ]
    
    # make phyloseq objects
    # tagl.subs   <- subset_samples(tagl.filtered, sample_names(tagl.filtered) %in% samps.subs)
    tagl.subs   <- prune_samples(samps.subs, tagl.filtered)
    
    # calculate spiec.easi object
    se.mb.subs  <- spiec.easi(tagl.subs, method="mb", lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=50))
    SE_objects[[ i ]] <- se.mb.subs
    
    # get igraphs with all nodes for hamming calculations
    ig2.mb.subs.allNodes <- adj2igraph(getRefit(se.mb.subs),  vertex.attr=list(name=taxa_names(tagl.subs)), rmEmptyNodes = F)
    igraphs.allNodes[[ i ]] <- ig2.mb.subs.allNodes
    
    # get igraphs without empty nodes for plotting
    ig2.mb.subs <- adj2igraph(getRefit(se.mb.subs),  vertex.attr=list(name=taxa_names(tagl.subs)), rmEmptyNodes = T)
    igraphs.noEmpty[[ i ]] <- ig2.mb.subs
    
    # calculating Hamming distances between CF and mC networks
    main_v_subs.hammings[[ i ]] <- netdist(ig2.mb.subs.allNodes, igraphs.main.allNodes, d = "Hamming")
    
  }
  
  # ************************** #
  # finally get pairwise distances between networks of all the matched control subsets
  asl_numbers <- names(allSubList[[ variable ]])
  
  subs_subs.hammings <- list()
  if (length(asl_numbers) > 1) {
    # ignore those variables for which subSamps were not calculated
    for ( i in asl_numbers[1]:(length(asl_numbers)-1) ) {
      for ( j in asl_numbers[2]:length(asl_numbers) ) {
        subs_subs.hammings[[ sprintf("%s_%s", i, j) ]] <- netdist(igraphs.allNodes[[ as.character(i) ]],
                                                                  igraphs.allNodes[[ as.character(j) ]],
                                                                  d = "Hamming")
      }
    }
  }
  
  # ************************** #
  return(list("se.mb.main"=se.mb.main,"igraphs.main.allNodes"=igraphs.main.allNodes,"igraphs.main.noEmpty"=igraphs.main.noEmpty,
              "se.mb.all_subs"=se.mb.all_subs,"igraphs.all_subs.allNodes"=igraphs.all_subs.allNodes,"igraphs.all_subs.noEmpty"=igraphs.all_subs.noEmpty,
              "SE_objects"=SE_objects,"igraphs.allNodes"=igraphs.allNodes,"igraphs.noEmpty"=igraphs.noEmpty,
              "main_v_subs.hammings"=main_v_subs.hammings,"subs_subs.hammings"=subs_subs.hammings))
}
# ************************************************************************* #

# ************************************************************************* #
get_networks_for_nonSubs <- function(variable, tagl, mTab, glomTab, minCount=15, minOcc=20) {
  
  var_vals <- as.character(unique(mTab[ , variable ]))
  var_vals <- var_vals[ ! is.na(var_vals) ]
  
  var_allSamps <- rownames(mTab[ ! is.na(mTab[ , variable]), ])
  
  # ************************** #
  # from this tutorial: http://psbweb05.psb.ugent.be/conet/microbialnetworks/spieceasi.php
  # filter out taxa that dont have at least 15 reads in at least 20 of the CF + mC samples, but with this function,
  #    keeps a row of all the counts that were removed, so it can keep the total counts for each sample
  g10 <- as.data.frame(glomTab[ , var_allSamps ])
  g10[ g10 <= minCount ] <- 0
  # filterobj <- filterTaxonMatrix(glomTab[ , var_allSamps ], minocc=100, keepSum = TRUE, return.filtered.indices = TRUE)
  filterobj <- filterTaxonMatrix(g10, minocc=minOcc, keepSum = TRUE, return.filtered.indices = TRUE)
  otus.f <- filterobj$mat
  
  # replace the values that were made to be 0s for the filtering
  anti.g10 <- as.data.frame(glomTab[ , var_allSamps ])
  anti.g10[ anti.g10 > minCount ] <- 0
  # replace them into g10
  g10[ anti.g10 != 0 ] <- anti.g10[ anti.g10 != 0 ]
  # first update the values that were made 0
  upd.otus_f <- otus.f
  upd.otus_f[ 1:(nrow(upd.otus_f)-1), ] <- as.data.frame(glomTab[ rownames(upd.otus_f)[ 1:(nrow(upd.otus_f)-1) ], var_allSamps ])
  # then update the summed row from the removed values
  upd.otus_f[ nrow(upd.otus_f), ] <- as.numeric(upd.otus_f[ nrow(upd.otus_f), ] + (colSums(upd.otus_f) - colSums(otus.f)))
  
  # finally update all the names appropriately, including pseudo-name for the summed counts, which will be ignored later in the plots
  taxa.f <- tax_table(tagl)[setdiff(1:nrow(tax_table(tagl)),filterobj$filtered.indices),]
  dummyTaxonomy <- c("k__dummy","p__","c__","o__","f__","g__")
  taxa.f <- rbind(taxa.f, dummyTaxonomy)
  rownames(taxa.f)[nrow(taxa.f)] <- "0"
  rownames(upd.otus_f)[nrow(upd.otus_f)] <- "0"
  
  tagl.filtered <- phyloseq(otu_table(upd.otus_f, taxa_are_rows = T),
                            sample_data(tagl)[ var_allSamps, ],
                            tax_table(taxa.f))
  
  # ************************** #
  
  vv_list <- list()
  
  for (val in var_vals) {
    print(val)
    vv_list[[ val ]] <- list()
    vv_samps <- rownames(mTab[ ! is.na(mTab[,variable]) & mTab[,variable] == val, ])
    vv_list[[ val ]][[ "samples" ]] <- vv_samps
    
    # first get network objects for main value, which does not change across subSamps
    vv_tagl  <- prune_samples(vv_samps, tagl.filtered)
    vv_se.mb <- spiec.easi(vv_tagl, method="mb", lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=50))
    vv_list[[ val ]][[ "se.mb" ]] <- vv_se.mb
    vv_list[[ val ]][[ "igraphs.allNodes" ]] <- adj2igraph(getRefit(vv_se.mb), vertex.attr=list(name=taxa_names(vv_tagl)), rmEmptyNodes = F)
    vv_list[[ val ]][[ "igraphs.noEmpty" ]] <- adj2igraph(getRefit(vv_se.mb), vertex.attr=list(name=taxa_names(vv_tagl)), rmEmptyNodes = T)
  }
  # ************************** #
  
  vv_list[[ "Hammings" ]] <- list()
  
  for (val1 in var_vals[1:(length(var_vals)-1)]) {
    for (val2 in var_vals[2:length(var_vals)]) {
      v1v2 <- sprintf("%s_%s", val1, val2)
      vv_list[[ "Hammings" ]][[ sprintf("%s_%s", val1, val2) ]] <- unname(netdist(vv_list[[ val1 ]]$igraphs.allNodes,
                                                                           vv_list[[ val2 ]]$igraphs.allNodes,
                                                                           d = "Hamming"))
        # c(vv_list[[ "Hammings" ]],
        #                            v1v2 = netdist(vv_list[[ val1 ]]$igraphs.allNodes,
        #                                           vv_list[[ val2 ]]$igraphs.allNodes,
        #                                           d = "Hamming"))
    }
  }
  
  vv_list[[ "Hammings" ]] <- unlist(vv_list[[ "Hammings" ]])
  
  # ************************** #
  return(vv_list)
}
# ************************************************************************* #

# ************************************************************************* #
get_networks_for_nonBinary_subSamps <- function(variable, tagl, allSubList, mTab, glomTab, minCount=15, minOcc=20) {
  # this function is a bit of a melding of the above functions for subSamps and nonSubs
  #   bc for Age_groups there were 100 subsamples taken, as usual, but for 2 groups together: Teen and Adult
  #   the function for subSamps only allowed for subsampling of binary variables, but accounted for the 100 subsamples
  #   the function for nonSubs accounted for potentially more than 2 values, but not for the 100 subsamples
  
  var_vals <- as.character(unique(mTab[ , variable ]))
  var_vals <- var_vals[ ! is.na(var_vals) ]
  
  # separate into the values that were or were not subsampled
  #   since there were relatively few Child and Senior samples, the same are always used
  #   while Teen and Adult samples were subSampled 100 times
  vv_nonSubs <- var_vals[ var_vals %in% c("Child","Senior") ]
  vv_subs <- var_vals[ var_vals %in% c("Teen","Adult") ]
  
  # first get objects for main value
  var_allSamps <- unique(unlist(allSubList[[ variable ]]))
  
  # ************************** #
  # from this tutorial: http://psbweb05.psb.ugent.be/conet/microbialnetworks/spieceasi.php
  # filter out taxa that dont have at least 15 reads in at least 20 of the CF + mC samples, but with this function,
  #    keeps a row of all the counts that were removed, so it can keep the total counts for each sample
  g10 <- as.data.frame(glomTab[ , var_allSamps ])
  g10[ g10 <= minCount ] <- 0
  # filterobj <- filterTaxonMatrix(glomTab[ , var_allSamps ], minocc=100, keepSum = TRUE, return.filtered.indices = TRUE)
  filterobj <- filterTaxonMatrix(g10, minocc=minOcc, keepSum = TRUE, return.filtered.indices = TRUE)
  otus.f <- filterobj$mat
  
  # replace the values that were made to be 0s for the filtering
  anti.g10 <- as.data.frame(glomTab[ , var_allSamps ])
  anti.g10[ anti.g10 > minCount ] <- 0
  # replace them into g10
  g10[ anti.g10 != 0 ] <- anti.g10[ anti.g10 != 0 ]
  # first update the values that were made 0
  upd.otus_f <- otus.f
  upd.otus_f[ 1:(nrow(upd.otus_f)-1), ] <- as.data.frame(glomTab[ rownames(upd.otus_f)[ 1:(nrow(upd.otus_f)-1) ], var_allSamps ])
  # then update the summed row from the removed values
  upd.otus_f[ nrow(upd.otus_f), ] <- as.numeric(upd.otus_f[ nrow(upd.otus_f), ] + (colSums(upd.otus_f) - colSums(otus.f)))
  
  # finally update all the names appropriately, including pseudo-name for the summed counts, which will be ignored later in the plots
  taxa.f <- tax_table(tagl)[setdiff(1:nrow(tax_table(tagl)),filterobj$filtered.indices),]
  dummyTaxonomy <- c("k__dummy","p__","c__","o__","f__","g__")
  taxa.f <- rbind(taxa.f, dummyTaxonomy)
  rownames(taxa.f)[nrow(taxa.f)] <- "0"
  rownames(upd.otus_f)[nrow(upd.otus_f)] <- "0"
  
  tagl.filtered <- phyloseq(otu_table(upd.otus_f, taxa_are_rows = T),
                            sample_data(tagl)[ var_allSamps, ],
                            tax_table(taxa.f))
  
  # ************************** #
  
  # first get networks for full sample sets of each value
  full_networks <- list()
  
  for (val in var_vals) {
    full_networks[[ val ]] <- list()
    vv_samps <- var_allSamps[ ! is.na(mTab[var_allSamps, variable]) & mTab[var_allSamps, variable] == val ]
    full_networks[[ val ]][[ "samples" ]] <- vv_samps
    
    vv_tagl  <- prune_samples(vv_samps, tagl.filtered)
    vv_se.mb <- spiec.easi(vv_tagl, method="mb", lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=50))
    full_networks[[ val ]][[ "se.mb" ]] <- vv_se.mb
    full_networks[[ val ]][[ "igraphs.allNodes" ]] <- adj2igraph(getRefit(vv_se.mb), vertex.attr=list(name=taxa_names(vv_tagl)), rmEmptyNodes = F)
    full_networks[[ val ]][[ "igraphs.noEmpty" ]] <- adj2igraph(getRefit(vv_se.mb), vertex.attr=list(name=taxa_names(vv_tagl)), rmEmptyNodes = T)
  }
  
  # then do so for each subSampling for the values that were subSampled
  subs_networks <- list()
  
  for (val in vv_subs) {
    print(val)
    subs_networks[[ val ]] <- list()
    
    # subs_networks[[ val ]][[ "main_v_subs.hammings" ]] <- list()
    subs_networks[[ val ]][[ "SE_objects" ]] <- list()
    subs_networks[[ val ]][[ "igraphs.allNodes" ]] <- list()
    subs_networks[[ val ]][[ "igraphs.noEmpty" ]] <- list()
    
    for (i in names(allSubList[[ variable ]])) {
      
      print(i)
      
      # get samples
      samps <- allSubList[[ variable ]][[ i ]]
      samps.subs <- samps[ mTab[samps, variable] == val ]
      
      # make phyloseq objects
      # tagl.subs   <- subset_samples(tagl.filtered, sample_names(tagl.filtered) %in% samps.subs)
      tagl.subs   <- prune_samples(samps.subs, tagl.filtered)
      
      # calculate spiec.easi object
      se.mb.subs  <- spiec.easi(tagl.subs, method="mb", lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=50))
      subs_networks[[ val ]][[ "SE_objects" ]][[ i ]] <- se.mb.subs
      
      # get igraphs with all nodes for hamming calculations
      ig2.mb.subs.allNodes <- adj2igraph(getRefit(se.mb.subs),  vertex.attr=list(name=taxa_names(tagl.subs)), rmEmptyNodes = F)
      subs_networks[[ val ]][[ "igraphs.allNodes" ]][[ i ]] <- ig2.mb.subs.allNodes
      
      # get igraphs without empty nodes for plotting
      ig2.mb.subs <- adj2igraph(getRefit(se.mb.subs),  vertex.attr=list(name=taxa_names(tagl.subs)), rmEmptyNodes = T)
      subs_networks[[ val ]][[ "igraphs.noEmpty" ]][[ i ]] <- ig2.mb.subs
      
      # # calculating Hamming distances between CF and mC networks
      # main_v_subs.hammings[[ i ]] <- netdist(ig2.mb.subs.allNodes, igraphs.main.allNodes, d = "Hamming")
    }
  }
  
  # ******************************* #
  return(list("full_networks"=full_networks, "subs_networks"=subs_networks))
}
# ************************************************************************* #




library(intergraph)
# ************************************************************************* #
connect_edges.subs <- function(variable, ssd) {
  
  # first get all connections for each genus in the "main" samples (non subsampled group - usually the "Yes" samples)
  ig.main.allNodes <- readRDS(sprintf("%s/%s/%s.igraphs.main.allNodes.rds", ssd, variable, variable))
  
  edge.mat.main <- as.matrix(asNetwork(ig.main.allNodes))
  # for each genus, get vector of those other genera with which it has an edge
  gen_edges.main <- list()
  for (genus in rownames(edge.mat.main)) {
    
    gen_edges.main[[ genus ]] <- colnames(edge.mat.main)[ edge.mat.main[ genus, ] != 0 ]
    
  }
  
  
  # ******************** #
  # do the same for each of the subsamplings of the other group (usually the "No" samples)
  ig.subs.allNodes <- readRDS(sprintf("%s/%s/%s.igraphs.allNodes.rds", ssd, variable, variable))
  
  edge.mats.subs <- list()
  gen_edges.subs <- list()
  for (i in names(ig.subs.allNodes)) {
    # get matrix of node connections
    edge.mats.subs[[ i ]] <- as.matrix(asNetwork(ig.subs.allNodes[[ i ]]))
    
    # get genus connection vectors
    gen_edges.subs[[ i ]] <- list()
    for (genus in rownames(edge.mats.subs[[ i ]])) {
      
      gen_edges.subs[[ i ]][[ genus ]] <- colnames(edge.mats.subs[[ i ]])[ edge.mats.subs[[ i ]][ genus, ] != 0 ]
      
    }
  }
  
  
  # ******************** #
  # then get the connections as follows:
  #   - for each genus, get vector of connected genera that are same in both main and subs
  #   - vector of genera only connected in main
  #   - vector of genera only connected in subs
  
  gen_connect <- list("all"=list(), "main"=list(), "subs"=list())
  
  for (i in names(gen_edges.subs)) {
    gen_connect$all[[ i ]] <- sapply(names(gen_edges.main), function(gen) 
      unique(c(gen_edges.main[[ gen ]][ gen_edges.main[[ gen ]] %in% gen_edges.subs[[ i ]][[ gen ]] ],
               gen_edges.subs[[ i ]][[ gen ]][ gen_edges.subs[[ i ]][[ gen ]] %in% gen_edges.main[[ gen ]] ])) )
    
    gen_connect$main[[ i ]] <- sapply(names(gen_edges.main), function(gen) 
      gen_edges.main[[ gen ]][ ! gen_edges.main[[ gen ]] %in% gen_edges.subs[[ i ]][[ gen ]] ])
    
    gen_connect$subs[[ i ]] <- sapply(names(gen_edges.main), function(gen) 
      gen_edges.subs[[ i ]][[ gen ]][ ! gen_edges.subs[[ i ]][[ gen ]] %in% gen_edges.main[[ gen ]] ])
  }
  
  
  
  
  # ******************** #
  # also get all connections for each genus in the full set of subs samples
  ig.all_subs.allNodes <- readRDS(sprintf("%s/%s/%s.igraphs.all_subs.allNodes.rds", ssd, variable, variable))
  
  edge.mat.subs.all_Subs <- as.matrix(asNetwork(ig.all_subs.allNodes))
  
  # for each genus, get vector of those other genera with which it has an edge
  gen_edges.subs.all_Subs <- list()
  for (genus in rownames(edge.mat.subs.all_Subs)) {
    
    gen_edges.subs.all_Subs[[ genus ]] <- colnames(edge.mat.subs.all_Subs)[ edge.mat.subs.all_Subs[ genus, ] != 0 ]
    
  }
  
  
  
  # ******************** #
  
  gen_connect.to_all_Subs <- list()
  
  gen_connect.to_all_Subs[[ "all" ]] <- sapply(names(gen_edges.main), function(gen) 
    gen_edges.main[[ gen ]][ gen_edges.main[[ gen ]] %in% gen_edges.subs.all_Subs[[ gen ]] ])
  
  gen_connect.to_all_Subs[[ "main" ]] <- sapply(names(gen_edges.main), function(gen) 
    gen_edges.main[[ gen ]][ ! gen_edges.main[[ gen ]] %in% gen_edges.subs.all_Subs[[ gen ]] ])
  
  gen_connect.to_all_Subs[[ "subs" ]] <- sapply(names(gen_edges.main), function(gen) 
    gen_edges.subs.all_Subs[[ gen ]][ ! gen_edges.subs.all_Subs[[ gen ]] %in% gen_edges.main[[ gen ]] ])
  
  
  # ************** #
  return(list("gen_connect"=gen_connect, "gen_connect.to_all_Subs"=gen_connect.to_all_Subs))
  
}
# ************************************************************************* #
# ************************************************************************* #
connect_edges.nonSubs <- function(nets.non, variable) {
  
  # first get the different values for the given variable for which network info is available
  vals <- names(nets.non[[ variable ]])[ names(nets.non[[ variable ]]) != "Hammings" ]
  
  igs.allNodes <- list()
  gen_edges <- list()
  
  for (val in vals) {
    igs.allNodes[[ val ]] <- nets.non[[ variable ]][[ val ]]$igraphs.allNodes
    
    edge.mat <- as.matrix(asNetwork(igs.allNodes[[ val ]]))
    # for each genus, get vector of those other genera with which it has an edge
    gen_edges[[ val ]] <- list()
    for (genus in rownames(edge.mat)) {
      
      gen_edges[[ val ]][[ genus ]] <- colnames(edge.mat)[ edge.mat[ genus, ] != 0 ]
      
    }
    
  }
  
  
  
  # ******************** #
  
  gen_connect <- list()
  
  # to get connections present in networks of all vals, have to do a few tryCatches depending on 
  #   first whether there are more than 1 genus connections for a given genus
  #   then whether or not there are any connections for that given genus
  gen_connect[[ "all" ]] <- sapply(names(gen_edges[[ vals[1] ]]), function(gen) {
    gen_edges[[ vals[1] ]][[ gen ]][ tryCatch(
      rowSums( sapply(vals[2:length(vals)], function(v) 
        gen_edges[[ vals[1] ]][[ gen ]] %in% gen_edges[[ v ]][[ gen ]]) ) == length(vals[2:length(vals)]),
      error = function(x) tryCatch(
        sum( sapply(vals[2:length(vals)], function(v) 
          gen_edges[[ vals[1] ]][[ gen ]] %in% gen_edges[[ v ]][[ gen ]]) ) == length(vals[2:length(vals)]),
        error = function(y) character(0))
      )]
  })
  
  # then for each val in vals, get those connections that are only present in that val
  for (val in vals) {
    gen_connect[[ val ]] <- sapply(names(gen_edges[[ val ]]), function(gen)
      gen_edges[[ val ]][[ gen ]][ ! gen_edges[[ val ]][[ gen ]] %in% 
                                     unique(unlist(lapply(gen_edges[ vals[ vals != val] ], function(x) x[[ gen ]])))
                                   ])
  }
  
  # finally, perhaps try combinations of vals if there are more than 2 (if only 2, would be the same as "all")
  
  
  # ******** #
  return(gen_connect)
  
}
# ************************************************************************* #

# *************************************** #
get_gen_con_freqs <- function(gc_list, gen, print_connections=FALSE) {
  
  gc_connections <- list()
  
  for (con in rev(names(gc_list))) {
    
    if (print_connections) {
      cat(sprintf("**** %s connections in %s %s ****\n", gen, con, ifelse(con=="all","","only")))
      print(sort(table(unlist(lapply(gc_list[[ con ]], function(x) x[[ gen ]])))))
    }
    
    gc_connections[[ con ]] <- sort(table(unlist(lapply(gc_list[[ con ]], function(x) x[[ gen ]]))))
  }
  
  return(gc_connections)
}
# *************************************** #
get_gen_con_freqs.all_No.or.nonSubs <- function(gc_list, gen, print_connections=FALSE) {
  
  gc_connections <- list()
  
  for (con in rev(names(gc_list))) {
    
    if (print_connections) {
      cat(sprintf("**** %s connections in %s %s ****\n", gen, con, ifelse(con=="all","","only")))
      print(sort(gc_list[[ con ]][[ gen ]]))
    }
    
    gc_connections[[ con ]] <- sort(gc_list[[ con ]][[ gen ]])
  }
  
  return(gc_connections)
}
# *************************************** #


# nc <- symBeta(getOptBeta(se.mb.MALDI))
# colnames(nc) <- rownames(nc) <- taxa_names(tagl.MALDI)
# inc <- graph.adjacency(nc, mode="undirected", add.rownames = T, weighted = T)
# edge.weights <- E(inc)$weight
# names(edge.weights) <- attr(E(inc),"vnames")
# names(edge.weights) <- sapply(names(edge.weights), function(x) paste( sort(strsplit(x, "\\|")[[1]]), collapse="."))
# 
# igm <- readRDS("/gpfs/projects/bsc40/current/jwillis/SLL/Part_2/R_objects/SpiecEasi/variable_subSamps/MALDI.Yeast_detected/MALDI.Yeast_detected.igraphs.main.allNodes.rds")
# as.matrix(asNetwork(igm))



# ************************************************************************* #
get_assoc_scores <- function(subVar, gcf, ss_sm_ms, anti_assoc=F) {
  
  if (ss_sm_ms == "main_main") {
    v1_sm <- "main"
    v2_sm <- "main"
  } else if (ss_sm_ms == "subs_subs") {
    v1_sm <- "subs"
    v2_sm <- "subs"
  } else if (ss_sm_ms == "subs_main") {
    v1_sm <- "subs"
    v2_sm <- "main"
  } else if (ss_sm_ms == "main_subs") {
    v1_sm <- "main"
    v2_sm <- "subs"
  }
  
  as_sco <- list()
  
  for (i in 1:(length(subVar)-1)) {
    for (j in (i+1):length(subVar)) {
      
      v1 <- subVar[ i ]
      v2 <- subVar[ j ]
      
      v1_v2 <- sprintf("%s-%s", v1, v2)
      print(v1_v2)
      as_sco[[ v1_v2 ]] <- list()
      
      # look at each genus as the first part of the association in each variable (except "0")
      gen_of_interest <- names(gcf[[ v1 ]])[ names(gcf[[ v1 ]]) != "0" ]
      
      for (gen_prime in gen_of_interest) {
        # check for each genus in assoc with the given gen_prime
        as_sco[[ v1_v2 ]][[ gen_prime ]] <- list()
        
        for (gen_con in gen_of_interest[ gen_of_interest != gen_prime ]) {
          # in each variable, check if gen_con is assoc at least once with gen_prime
          # if not, gets a count of 0
          v1_gg <- ifelse(gen_con %in% names(gcf[[ v1 ]][[ gen_prime ]][[ v1_sm ]]),
                          gcf[[ v1 ]][[ gen_prime ]][[ v1_sm ]][ gen_con ], 0)
          v2_gg <- ifelse(gen_con %in% names(gcf[[ v2 ]][[ gen_prime ]][[ v2_sm ]]),
                          gcf[[ v2 ]][[ gen_prime ]][[ v2_sm ]][ gen_con ], 0)
          
          # **************** #
          if (anti_assoc == F) {
            # then the score first calculates the difference between v1 and v2
            # subtracts from 100 (ideal value for an association)
            # finally weight by the mean percentage of subsamplings that showed that association in the 2 variables
            score <- (100 - abs(v1_gg - v2_gg)) * mean(c(v1_gg, v2_gg))/100
            
          } else if (anti_assoc == T) {
            # in this case, simply subtract v2 from v1
            score <- v1_gg - v2_gg
          }
          # **************** #
          
          as_sco[[ v1_v2 ]][[ gen_prime ]][[ gen_con ]] <- score
        }
        
        # then unlist the list of gen_cons, to have just a numeric vector
        as_sco[[ v1_v2 ]][[ gen_prime ]] <- unlist(as_sco[[ v1_v2 ]][[ gen_prime ]])
      }
      
    }
  }
  return(as_sco)
}
# ************************************************************************* #
# ************************************************************************* #
final_assoc_scores <- function(assoc_score, typeCombo, vars, tagl, subSamps, glomTab, tl, includeNegs=F) {
  
  sampType1 <- strsplit(typeCombo, "_")[[1]][1]
  sampType2 <- strsplit(typeCombo, "_")[[1]][2]

  if (sampType1 == "main") {
    # varVals <- c("Yes","One","F")
    network_file <- "main"
  } else {
    # varVals <- c("No","Two","M")
    network_file <- "all_subs"
  }
  
  var_scores_total <- list()
  var_to_var_scores <- list()
  
  for (v in vars) {
    print(v)
    var_to_var_scores[[ v ]] <- list()
    
    # ************************ #
    if (v %in% c("Cystic_fibrosis","Downs_Syndrome","Celiac","Hypertension"))
      mTab <- SLL2.meta
    else
      mTab <- meta.healthy
    
    # if (v == "MALDI.Yeast_detected") {
    #   # use no in this case, since according to the homogeneity test, the "No" samples are more consistent
    #   varVals <- c("No","Two","M")
    #   network_file <- "all_subs"
    # } else {
    #   varVals <- c("Yes","One","F")
    #   network_file <- "main"
    # }
    
    tagl.for_network <- get_tagl_for_networks(v, tagl, subSamps, mTab, glomTab[[ tl ]])
    
    # netSamps <- unique(unlist(lapply(subSamps[[ v ]], function(i) i)))
    # netSamps <- netSamps[ mTab[ netSamps, v] %in% varVals ]
    # net_tagl  <- subset_samples(tagl.for_network, sample_names(tagl.for_network) %in% netSamps)
    se.mb <- readRDS(sprintf("%s/%s/%s.se.mb.%s.rds", spiec_sub_dir, v, v, network_file))
    # ************************ #
    
    nc <- symBeta(getOptBeta(se.mb))
    colnames(nc) <- rownames(nc) <- taxa_names(tagl.for_network)
    inc <- graph.adjacency(nc, mode="undirected", add.rownames = T, weighted = T)
    edge.weights <- E(inc)$weight
    names(edge.weights) <- attr(E(inc),"vnames")
    # make names of edge weights match the names of the associations in assoc_score,
    #   and make sure taxa are alphabetized in the names as well
    names(edge.weights) <- sapply(names(edge.weights), function(x) paste( sort(strsplit(x, "\\|")[[1]]), collapse="."))
    # ************************ #
    
    for (varPair in names(assoc_score)) {
      
      pair_vs <- strsplit(varPair, "-")[[1]]
      
      if ( v %in% pair_vs ) {
        print(sprintf("   %s", varPair))
        
        # ********************************************* #
        if (includeNegs == T) {
          assocs <- assoc_score[[ varPair ]][ assoc_score[[ varPair ]] != 0 ]
          assocs <- assocs[ names(assocs) %in% names(edge.weights) ]
          if (startsWith(varPair, v)) {
            # no need to change sign
            assocs <- assocs * abs(edge.weights[ names(assocs) ])
          } else if (endsWith(varPair, v)) {
            # flip sign bc positive is for the other variable in this case
            assocs <- (-1) * assocs * abs(edge.weights[ names(assocs) ])
          }
          
          # ********************************************* #
        } else {
          # ************************ #
          if (startsWith(varPair, v)) {
            # take only the positive scores
            assocs <- assoc_score[[ varPair ]][ assoc_score[[ varPair ]] > 0 ]
            
          } else if (endsWith(varPair, v)) {
            # take only the negative scores
            assocs <- assoc_score[[ varPair ]][ assoc_score[[ varPair ]] < 0 ]
            
          }
          # ************************ #
          assocs <- abs(assocs * edge.weights[ names(assocs) ])
          # ************************ #
        }
        # ********************************************* #
        
        var_to_var_scores[[ v ]][[ pair_vs[ pair_vs != v ] ]] <- sum(assocs)
      }
    }
    # ************************ #
    var_scores_total[[ v ]] <- sum(unlist( var_to_var_scores[[ v ]] ))
  }
  return(list("var_scores_total"=var_scores_total, "var_to_var_scores"=var_to_var_scores))
}
# ************************************************************************* #

library(nettools)
library(SpiecEasi)
library(seqtime)
# **************************************************************************** #
get_tagl_for_networks <- function(variable, tagl, subSamps, mTab, glomTab, minCount=15, minOcc=20) {
  
  # how many are "Yes" for given variable, to know how many to match ("Two" for seqGroup)
  main_vals <- c("Yes","Two","F")
  subs_vals <- c("No","One","M")
  # have to use the "No" as main group here since there are many more "Yes" samples
  if (variable %in% c("Fluoride_toothpaste","Wash_hands_before_eat","Wash_hands_after_bathroom","Do_you_feel_well",
                      "MALDI.Yeast_detected","Full_MALDI.Candida")) {
    main_vals <- c("No","One","M")
    subs_vals <- c("Yes","Two","F")
  }
  
  # first get objects for main value
  samps.all <- unique(unlist(subSamps[[ variable ]]))
  samps.main <- samps.all[ mTab[samps.all, variable] %in% main_vals ]
  # then for all of the matched controls together
  samps.subs <- samps.all[ mTab[samps.all, variable] %in% subs_vals ]
  
  # ************************** #
  # from this tutorial: http://psbweb05.psb.ugent.be/conet/microbialnetworks/spieceasi.php
  # filter out taxa that dont have at least 15 reads in at least 20 of the Yes + mC samples, but with this function,
  #    keeps a row of all the counts that were removed, so it can keep the total counts for each sample
  g10 <- as.data.frame(glomTab[ , c(samps.main, samps.subs) ])
  g10[ g10 <= minCount ] <- 0
  # filterobj <- filterTaxonMatrix(glomTab[ , c(samps.main, samps.subs) ], minocc=100, keepSum = TRUE, return.filtered.indices = TRUE)
  filterobj <- filterTaxonMatrix(g10, minocc=minOcc, keepSum = TRUE, return.filtered.indices = TRUE)
  otus.f <- filterobj$mat
  
  # replace the values that were made to be 0s for the filtering
  anti.g10 <- as.data.frame(glomTab[ , c(samps.main, samps.subs) ])
  anti.g10[ anti.g10 > minCount ] <- 0
  # replace them into g10
  g10[ anti.g10 != 0 ] <- anti.g10[ anti.g10 != 0 ]
  # first update the values that were made 0
  upd.otus_f <- otus.f
  upd.otus_f[ 1:(nrow(upd.otus_f)-1), ] <- as.data.frame(glomTab[ rownames(upd.otus_f)[ 1:(nrow(upd.otus_f)-1) ], c(samps.main, samps.subs) ])
  # then update the summed row from the removed values
  upd.otus_f[ nrow(upd.otus_f), ] <- as.numeric(upd.otus_f[ nrow(upd.otus_f), ] + (colSums(upd.otus_f) - colSums(otus.f)))
  
  # finally update all the names appropriately, including pseudo-name for the summed counts, which will be ignored later in the plots
  taxa.f <- tax_table(tagl)[setdiff(1:nrow(tax_table(tagl)),filterobj$filtered.indices),]
  dummyTaxonomy <- c("k__dummy","p__","c__","o__","f__","g__")
  taxa.f <- rbind(taxa.f, dummyTaxonomy)
  rownames(taxa.f)[nrow(taxa.f)] <- "0"
  rownames(upd.otus_f)[nrow(upd.otus_f)] <- "0"
  
  tagl.network <- phyloseq(otu_table(upd.otus_f, taxa_are_rows = T),
                           sample_data(tagl)[ c(samps.main, samps.subs), ],
                           tax_table(taxa.f))
  # ************************** #
  return(tagl.network)
}
# **************************************************************************** #
get_network_objects <- function(sampType, spiecObj, phy, taxTab, printEdgeInfo=TRUE, i=NULL, 
                                plotNet=TRUE, addLabels=FALSE, Vprop=FALSE, Eprop=FALSE,
                                taxa_to_label=NULL, addLegend=TRUE, legendOnly=FALSE) {
  
  if (sampType == "mC_i") {
    # update objects as needed 
    spiecObj <- spiecObj[[ i ]]
  }
  
  n.c <- symBeta(getOptBeta( spiecObj ))
  betaMat <- as.matrix( n.c )
  print(dim(n.c))
  print(phy)
  print(n.c[1:5,1:5])
  print(summary(colnames(spiecObj$est$data) %in% taxa_names(phy)))
  phy <- prune_taxa(colnames(spiecObj$est$data), phy)
  print(phy)
  # to add abundance values to vertices
  colnames(n.c) <- rownames(n.c) <- taxa_names(phy)#rownames(gloms$Genus)
  # add log abundance as properties of vertex/nodes.
  # vsize <- log2(apply(gloms$Genus[ , samps.CF ], 1, mean)+1)
  vsize <- log2(apply(phy@otu_table, 1, mean)+1)
  
  # get number of positive and negative edges
  # We divide by two since an edge is represented by two entries in the matrix. 
  pos_edge <- length(betaMat[ betaMat > 0 ])/2
  neg_edge <- length(betaMat[ betaMat < 0 ])/2
  tot_edge <- length(betaMat[ betaMat != 0 ])/2
  
  
  # prepare data for plotting
  ig <- graph.adjacency(n.c, mode='undirected', add.rownames = TRUE, weighted = TRUE)
  # remove those nodes with no edges
  ig <- delete.vertices(ig, V(ig)[degree(ig) == 0])
  ig <- delete.vertices(ig, V(ig)[name == "0"])
  
  if (printEdgeInfo==TRUE){
    # # we can see all the attributes and weights
    print("**** EDGES ****")
    print(ig)
  }
  
  # one more object, to be able to check weights and genus names together
  ig.unchanged <- ig
  if (printEdgeInfo==TRUE){
    # # we can see all the attributes and weights
    cat("\n")
    print("**** POSITIVE EDGES ****")
    print(E(ig.unchanged)[weight > 0])
    cat("\n")
    print("**** NEGATIVE EDGES ****")
    print(E(ig.unchanged)[weight < 0])
    # print(E(ig.unchanged)[weight < 0]$weight)
    # print((log(abs(E(ig.unchanged)[weight < 0]$weight)+1, base = 50)))
  }
  
  # now color the edges based on their values positive is steelblue
  # E(igraphs.CF.Yes.noEmpty)[weight > 0]$color<-"steelblue"
  # # now color the edges based on their values
  # E(igraphs.CF.Yes.noEmpty)[weight < 0]$color<-"orange"
  E(ig)[weight > 0]$color <- "#73BDD3" #"steelblue"
  # now color the edges based on their values
  E(ig)[weight < 0]$color <- "#CB6767" # "orange"
  
  # before making coordinates, must remove weights
  ig.no_weight <- ig
  E(ig.no_weight)$weight <- 1
  
  # make coordinates
  legend.y <- 1.2
  legend.c <- 1.2
  if (sampType == "CF") {
    set.seed(172)
    net.title <- "(a) Co-occurrence network: CF samples"
    legend.x <- -1.83
  } else if (sampType == "mC") {
    set.seed(711)
    net.title <- "(b) Co-occurrence network: matched control samples"
    legend.x <- -1
  } else if (sampType == "mC_i") {
    net.title <- sprintf("Co-occurrence network: matched controls - subsampling #%s", i)
  } else if (sampType == "Child") {
    set.seed(111)
    net.title <- "() Co-occurrence network: Child samples"
  } else if (sampType == "Teen") {
    set.seed(111)
    net.title <- "() Co-occurrence network: Teen samples"
  } else if (sampType == "Adult") {
    set.seed(111)
    net.title <- "() Co-occurrence network: Adult samples"
  } else if (sampType == "Senior") {
    set.seed(111)
    net.title <- "() Co-occurrence network: Senior samples"
  } else if (sampType == "Younger") {
    set.seed(111)
    net.title <- "() Co-occurrence network: Younger samples"
  } else if (sampType == "Older") {
    set.seed(111)
    net.title <- "() Co-occurrence network: Older samples"
  } else if (sampType == "Smoker") {
    set.seed(172)
    net.title <- "(a) Co-occurrence network: Smokers"
    legend.x <- -1.83
  } else if (sampType == "MALDI") {
    set.seed(172)
    net.title <- "(a) Co-occurrence network: MALDI.Yeast_detected"
    legend.x <- -1.83
  }
  
  coords.fr <- layout_with_fr(ig.no_weight)#igraphs.CF.Yes.noEmpty)
  
  # color dots by phylum
  # prepare colors for plot
  colPal.cols <- c("#e38f87","#9abddf","#90C98A",
                   "#FFBC67","#a854a8","#ab8b67",
                   "darkblue","#bbbbbb","pink","#56423c",
                   "brown","black","gray","darkorange3")
  names(colPal.cols) <- c("Firmicutes","Bacteroidetes","Proteobacteria",
                          "Fusobacteria","Actinobacteria","Epsilonbacteraeota",
                          "Patescibacteria","Spirochaetes","Synergistetes","unclassified.P1",
                          "Verrucomicrobia","Tenericutes","Cyanobacteria","Chloroflexi")
  
  V(ig)$name <- unname(taxTables.both$Genus[V(ig)$name,"Phylum"])
  # V(ig)$color <- match(V(ig)$name, unique(V(ig)$name))
  V(ig)$color <- colPal.cols[ V(ig)$name ]
  
  # give same colors to no_weight graph, just in case
  E(ig.no_weight)$color <- E(ig)$color
  V(ig.no_weight)$color <- V(ig)$color
  
  
  # if indicated, print only the desired taxa names on nodes
  ig.no_weight.taxLabs <- ig.no_weight
  vertLabCex <- 0.75
  if ( ! is.null(taxa_to_label) ) {
    V(ig.no_weight.taxLabs)$name[ ! V(ig.no_weight.taxLabs)$name %in% taxa_to_label ] <- ""
    vertLabCex <- 1.8
    
    # update some unclassified names
    V(ig.no_weight.taxLabs)$name[ V(ig.no_weight.taxLabs)$name == "unclassified.G10" ] <- "F.Lentimicrobiaceae.UCG"
    # V(ig.no_weight.taxLabs)$name[ V(ig.no_weight.taxLabs)$name == "unclassified.G20" ] <- "F.Family_XIII.UCG"
    V(ig.no_weight.taxLabs)$name[ V(ig.no_weight.taxLabs)$name == "unclassified.G20" ] <- "Family_XIII"
    V(ig.no_weight.taxLabs)$name[ V(ig.no_weight.taxLabs)$name == "unclassified.G31" ] <- "C.Gracilibacteria.UCG"
    V(ig.no_weight.taxLabs)$name[ V(ig.no_weight.taxLabs)$name == "unclassified.G32" ] <- "F.Saccharimonadaceae.UCG"
    V(ig.no_weight.taxLabs)$name[ V(ig.no_weight.taxLabs)$name == "unclassified.G33" ] <- "O.Saccharimonadales.UCG"
  }
  
  
  # update some unclassified names
  # V(ig.no_weight)$name[ V(ig.no_weight)$name == "unclassified.G20" ] <- "Family_XIII"
  # V(ig.no_weight)$name[ V(ig.no_weight)$name == "unclassified.G31" ] <- "C.Gracilibacteria.UCG"
  # V(ig.no_weight)$name[ V(ig.no_weight)$name == "unclassified.G33" ] <- "O.Saccharimonadales.UCG"
  
  
  # ******************** #
  # ******************** #
  
  # make the negative edges a bit heavier since they tend to always be super weak
  E(ig)[weight < 0]$weight <- E(ig)[weight < 0]$weight * 3
  
  # netPlot <- netLegend <- NULL
  
  if (plotNet == TRUE) {
    
    # plot(ig, layout=coords.fr, vertex.size = vsize[ V(ig.no_weight)$name ], 
    #      vertex.label.cex = 0.75, edge.width = (abs(E(ig)$weight))*30, 
    #      main=net.title)
    vs.plot <- 6
    ew.plot <- 1
    
    if (Vprop == TRUE) vs.plot <- vsize[ V(ig.no_weight)$name ]
    
    if (Eprop == TRUE) ew.plot <- (abs(E(ig)$weight))*30
    # if (Eprop == TRUE) ew.plot <- log(abs(E(ig)$weight)+1, base = 100)*10
    
    if (addLabels == FALSE & legendOnly == FALSE) 
      plot(ig, layout=coords.fr, 
           vertex.size = vs.plot, edge.width = ew.plot, 
           vertex.label.cex = 0.001,
           main=net.title)
    
    else if (addLabels == TRUE & legendOnly == FALSE) 
      plot(ig.no_weight.taxLabs, layout=coords.fr, 
           vertex.size = vs.plot, edge.width = ew.plot, 
           vertex.label.cex = vertLabCex, vertex.label.color="darkblue",
           main=net.title)
    
    # ******************** #
    # if (addLabels == FALSE & Vprop == FALSE & Eprop == FALSE)
    #   # (1) plot with no labels or proportional anything
    #   plot(ig, layout=coords.fr, vertex.size = 6, 
    #        vertex.label.cex = 0.001, #edge.width = exp(abs(E(ig.CF_Yes)$weight))*3, 
    #        main=net.title)
    # 
    # if (addLabels == FALSE & Vprop == TRUE & Eprop == FALSE)
    #   # (2) plot with proportional vertex sizes
    #   plot(ig, layout=coords.fr, vertex.size = vsize[ V(ig.no_weight)$name ], 
    #        vertex.label.cex = 0.001, #edge.width = exp(E(ig.CF_Yes)$weight),
    #        main=net.title)
    # 
    # if (addLabels == FALSE & Vprop == FALSE & Eprop == TRUE)
    #   # (3) plot with proportional edge widths
    #   plot(ig, layout=coords.fr, vertex.size = 6, 
    #        vertex.label.cex = 0.001, edge.width = (abs(E(ig)$weight))*30, 
    #        main=net.title)
    # 
    # if (addLabels == FALSE & Vprop == TRUE & Eprop == TRUE)
    #   # (4) plot with proportional edge widths and vertex sizes
    #   plot(ig, layout=coords.fr, vertex.size = vsize[ V(ig.no_weight)$name ], 
    #        vertex.label.cex = 0.001, edge.width = (abs(E(ig)$weight))*30, 
    #        main=net.title)
    # 
    # if (addLabels == TRUE & Vprop == TRUE & Eprop == TRUE)
    #   # (5) plot with genus labels, proportional edge widths and vertex sizes
    #   plot(ig.no_weight.taxLabs, layout=coords.fr, vertex.size = vsize[ V(ig.no_weight)$name ], 
    #        vertex.label.cex = vertLabCex, edge.width = (abs(E(ig)$weight))*30, 
    #        main=net.title)
    # 
    # if (addLabels == TRUE & Vprop == FALSE & Eprop == FALSE)
    #   # (6) plot with genus labels
    #   plot(ig.no_weight.taxLabs, layout=coords.fr, vertex.size = 6, 
    #        vertex.label.cex = vertLabCex, #edge.width = exp(abs(E(ig)$weight))*3, 
    #        main=net.title)
    
    # ******************** #
    if (addLegend == TRUE | legendOnly == TRUE) {
      if (legendOnly == TRUE) {
        plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
        legend.x <- 0.25
        legend.y <- 0.75
        legend.c <- 1.75
      }
      legend(x=legend.x, y=legend.y, inset=0.2, title="Phylum", 
             gsub("unclassified.P1", "unclassified", sort(unique(V(ig)$name))),
             fill = colPal.cols[ sort(unique(V(ig)$name)) ], 
             cex=legend.c, bty="n")
    }
    # ******************** #
  }
  
  
  
  # ******************** #
  # ******************** #
  return(list("pos_edge"=pos_edge, "neg_edge"=neg_edge, "tot_edge"=tot_edge,
              "vsize"=vsize, "ig"=ig, "ig.no_weight"=ig.no_weight, 
              "coords.fr"=coords.fr, "colPal.cols"=colPal.cols, "net.title"=net.title))
  # ******************** #
  # ******************** #
}
# **************************************************************************** #



# ****************************************************************************************************************** #




























# ****************************************************************************************************************** ####
# Plots for subsamplings ####

# ************************************************************************************** #
subsampling.plots.box <- function(disorder, tl, trait, divOrTax, GQ, OC, glomTab, mTab, phy,
                                  dst, dstStruc="disorder.original", ignore_pvals=F, ignore_dstStruc=F, additional.dst=NULL,
                                  sampMode="default", singleSubSamp=NULL, chosenSamps=NULL, 
                                  print_tukey=F, plot_tukey=F, adjustAlphas=F,
                                  plotType=NULL, dist_meas=NULL, swapShapeColor=F, use_ggtitle=T, use_xlab=T,
                                  save_plot=F, plot_dir=NULL, xAngle=0, useLoess=F, specify_facet_labels=NULL,
                                  facetVar=NULL, facetRows=NULL, facetScales="free", forceFacet=F,
                                  choose_xl=NULL, choose_yl=NULL, choose_title=NULL, str.bck="grey90",
                                  atex.s=17, atit.s=18, pts=17, sts=13, flipCoords=F, include_legend_below=F) {
  
  # *************************** #
  if (startsWith(divOrTax, "Stomatotype")) {
    stop("ERROR: Stomatotypes should not be plotted with this function because they are different in each subsampling")
  }
  # *************************** #
  
  # *************************** #
  # *************************** #
  if (ignore_pvals == T & ignore_dstStruc == T) {
    mps <- NULL
    
    if (class(chosenSamps) == "list")
      samps_of_int <- unique(unlist(lapply(chosenSamps, function(x) x)))
    else
      samps_of_int <- chosenSamps
    
    mTab.sub.all <- mTab[ samps_of_int, ]
    phy.sub.all  <- prune_samples( samps_of_int, phy )
    # *************************** #
    # *************************** #
  } else {
    
    # *************************** #
    if (dstStruc == "Age_groups") {
      if ( ! is.null(chosenSamps)) {
        
        if (class(chosenSamps) == "list")
          samps_of_int <- unique(unlist(lapply(chosenSamps, function(x) x)))
        else
          samps_of_int <- chosenSamps
        
        mTab.sub.all <- mTab[ samps_of_int, ]
        phy.sub.all  <- prune_samples( samps_of_int, phy )
      } else if ( ! is.null(singleSubSamp) ) {
        samps_of_int <- dst[[ singleSubSamp ]]$samples
        mTab.sub.all <- mTab[ samps_of_int, ]
        phy.sub.all  <- prune_samples( samps_of_int, phy )
      } else {
        ag.all <- unique(unlist(lapply(dst, function(x) x$samples)))
        mTab.sub.all <- mTab[ ag.all, ]
        phy.sub.all <- prune_samples( ag.all, phy )
        samps_of_int <- ag.all
      }
      
    } else {
      # trying to get average values of all pvalues in anova tests, then plot this with all the controls used
      ds.conts.all <- unique(unlist(lapply(dst[[ sampMode ]][[ disorder ]], function(x) x$samples)))
      mTab.sub.all <- mTab[ds.conts.all, ]
      phy.sub.all <- prune_samples( ds.conts.all, phy)
      DS.samps <- ds.conts.all[ mTab.sub.all[ds.conts.all, disorder]=="Yes" ]
      cont.samps.all <- ds.conts.all[ mTab.sub.all[ds.conts.all, disorder]=="No" ]
      samps_of_int <- c(DS.samps, cont.samps.all)
      
      if (! is.null(singleSubSamp) ) {
        samps_of_int <- dst[[ sampMode ]][[ disorder ]][[ singleSubSamp ]]$samples
        mTab.sub.all <- mTab[ samps_of_int, ]
        phy.sub.all  <- prune_samples( samps_of_int, phy )
      }
      
    }
    
    # *************************** #
    
    
    # *************************** #
    # *************************** #
    # *************************** #
    # here are the mean pvalues for Downs_Syndrome and the other 3 covars 
    #  problem is that when using na.rm=T maybe only like 3 instances in which it worked
    if (length(divOrTax)==1) {
      
      if (dstStruc == "disorder.original") {
        mps <- as.data.frame(sapply(colnames(dst[[ sampMode ]][[ disorder ]]$`1`$Anova[[ tl ]]), function(y) 
          mean(p.adjust(unlist(lapply(dst[[ sampMode ]][[ disorder ]], function(x) {
            val <- x$Anova[[ tl ]][gsub("-","\\.",divOrTax), y]
            # replace pvals that are NA with 1
            ifelse(is.na(val), 1, val)
          })), method = "fdr" ) )))#, 
        # na.rm=F)))
      } else if (dstStruc == "Age_groups") {
        mps <- as.data.frame(sapply(colnames(dst$`1`$Anova[[ tl ]]), function(y) 
          mean(p.adjust(unlist(lapply(dst, function(x) {
            val <- x$Anova[[ tl ]][gsub("-","\\.",divOrTax), y]
            # replace pvals that are NA with 1
            ifelse(is.na(val), 1, val)
          })), method = "fdr" ) )))#, 
        # na.rm=F)))
      } else if (dstStruc == "additional") {
        mps <- as.data.frame(sapply(colnames(additional.dst[[ disorder ]]$`1`[[ tl ]]), function(y) 
          mean(p.adjust(unlist(lapply(additional.dst[[ disorder ]], function(x) {
            val <- x[[ tl ]][gsub("-","\\.",divOrTax), y]
            # replace pvals that are NA with 1
            ifelse(is.na(val), 1, val)
          })), method = "fdr" ) )))#, 
        # na.rm=F)))
      } else if (dstStruc == "water") {
        mps <- as.data.frame(sapply(colnames(additional.dst[[ disorder ]]$`1`[[ tl ]][[ trait ]]), function(y) 
          mean(p.adjust(unlist(lapply(additional.dst[[ disorder ]], function(x) {
            val <- x[[ tl ]][[ trait ]][gsub("-","\\.",divOrTax), y]
            # replace pvals that are NA with 1
            ifelse(is.na(val), 1, val)
          })), method = "fdr" ) )))#, 
        # na.rm=F)))
      }
      
      colnames(mps) <- divOrTax
      mps <- as.matrix(t(mps))
      print(mps)
      mps[2:length(mps)][ mps[2:length(mps)] > 0.05 ] <- ""
      # *************************** #
      
    } else {
      
      # *************************** #
      if (dstStruc == "disorder.original") {
        mps <- sapply(gsub("-","\\.",divOrTax), function(z) 
          as.data.frame(sapply(colnames(dst[[ sampMode ]][[ disorder ]]$`1`$Anova[[ tl ]]), function(y) 
            mean(p.adjust(unlist(lapply(dst[[ sampMode ]][[ disorder ]], function(x) {
              val <- x$Anova[[ tl ]][z, y]
              # replace pvals that are NA with 1
              ifelse(is.na(val), 1, val)
            })), method = "fdr" ) )))#, 
          # na.rm=F)))
        )
        mps_cols <- colnames(dst[[ sampMode ]][[ disorder ]]$`1`$Anova[[ tl ]])
        
      } else if (dstStruc == "Age_groups") {
        mps <- sapply(gsub("-","\\.",divOrTax), function(z) 
          as.data.frame(sapply(colnames(dst$`1`$Anova[[ tl ]]), function(y) 
            mean(p.adjust(unlist(lapply(dst, function(x) {
              val <- x$Anova[[ tl ]][z, y]
              # replace pvals that are NA with 1
              ifelse(is.na(val), 1, val)
            })), method = "fdr" ) )))#, 
          # na.rm=F)))
        )
        mps_cols <- colnames(dst$`1`$Anova[[ tl ]])
        
      } else if (dstStruc == "additional") {
        mps <- sapply(gsub("-","\\.",divOrTax), function(z) 
          as.data.frame(sapply(colnames(additional.dst[[ disorder ]]$`1`[[ tl ]]), function(y) 
            mean(p.adjust(unlist(lapply(additional.dst[[ disorder ]], function(x) {
              val <- x[[ tl ]][z, y]
              # replace pvals that are NA with 1
              ifelse(is.na(val), 1, val)
            })), method = "fdr" ) )))#, 
          # na.rm=F)))
        )
        mps_cols <- colnames(dst[[ sampMode ]][[ disorder ]]$`1`$Anova[[ tl ]])
        
      } else if (dstStruc == "combo") {
        # *************************** #
        mps <- sapply(gsub("-","\\.",divOrTax), function(z) {
          if (z %in% c("pH","BMI","MALDI.Num_Yeast_Colonies","MALDI.Num_Mold_Colonies")) {
            as.data.frame(sapply(colnames(additional.dst[[ disorder ]]$`1`$contVar), function(y) 
              mean(p.adjust(unlist(lapply(additional.dst[[ disorder ]], function(x) {
                val <- x$contVar[z, y]
                # replace pvals that are NA with 1
                ifelse(is.na(val), 1, val)
              })), method = "fdr" ) )))
          } else if (z %in% rownames(dst[[ sampMode ]][[ disorder ]]$`1`$Anova$contVar)) {
            as.data.frame(sapply(colnames(dst[[ sampMode ]][[ disorder ]]$`1`$Anova$contVar), function(y) 
              mean(p.adjust(unlist(lapply(dst[[ sampMode ]][[ disorder ]], function(x) {
                val <- x$Anova$contVar[z, y]
                # replace pvals that are NA with 1
                ifelse(is.na(val), 1, val)
              })), method = "fdr" ) )))#, 
            # na.rm=F)))
          } else {
            as.data.frame(sapply(colnames(dst[[ sampMode ]][[ disorder ]]$`1`$Anova[[ tl ]]), function(y) 
              mean(p.adjust(unlist(lapply(dst[[ sampMode ]][[ disorder ]], function(x) {
                val <- x$Anova[[ tl ]][z, y]
                # replace pvals that are NA with 1
                ifelse(is.na(val), 1, val)
              })), method = "fdr" ) )))#, 
            # na.rm=F)))
          }
        }
        )
        # *************************** #
        mps_cols <- colnames(dst[[ sampMode ]][[ disorder ]]$`1`$Anova[[ tl ]])
      }
      
      # *************************** #
      # *************************** #
      
      mps <- as.data.frame(mps)
      colnames(mps) <- divOrTax
      rownames(mps) <- mps_cols
      mps <- as.matrix(t(mps))
      print(mps)
      mps[ , 2:ncol(mps) ][ mps[ , 2:ncol(mps) ] > 0.05 ] <- ""
      print(mps)
      
      
    }
  }
  
  
  # *************************** #
  # *************************** #
  # *************************** #
  
  if (ignore_pvals == T)
    mps <- NULL
  
  
  # ******************** #
  if (tl == "contVar") {
    # arbitrarily make tlev Genus, this value wont be used in this case
    tlev <- "Genus"
    comparison <- "questions"
  } else {
    tlev <- tl
    comparison <- "TvsQ"
  }
  # ******************** #
  
  # # include the maldi values that were created later, in case making plots for those
  # GQ <- c(GQ, "Full_MALDI.Candida_albicans","Full_MALDI.Candida_dubliniensis","Full_MALDI.Candida_glabrata",
  #         "Full_MALDI.Candida_parapsilosis","Full_MALDI.Candida_guillermondii","Full_MALDI.Candida_intermedia",
  #         "Full_MALDI.Candida_krusei","Full_MALDI.Candida_lusitaniae","Full_MALDI.Cryptococcus_spp",
  #         "Full_MALDI.Debaryomyces_hansenii","Full_MALDI.Rhodotorula_mucilaginosa")
  
  
  # *************************** #
  if (plotType == "box") {
    
    group_vs_cont_box(mTab.sub.all, samps_of_int, tlev, glomTab[[ tlev ]], 
                      trait, divOrTax, GQ, OC, print_tukey = print_tukey, plot_tukey = plot_tukey,
                      save.boxes = F, plot_dirs = sprintf("LMM_signifs/%s/%s", nor, f),
                      anova.pTab = mps, xAngle = xAngle, flipCoords = flipCoords, include_legend_below=include_legend_below,
                      facet = facetVar, facetRows = facetRows, facetScales = facetScales, 
                      choose_yl=choose_yl, choose_xl=choose_xl, choose_title=choose_title, str.bck=str.bck,
                      atex.s=atex.s, atit.s=atit.s, pts=pts, sts=sts, adjustAlphas = adjustAlphas)
    
  } else if (plotType == "scatter") {
    
    plot_data.cont(trait, divOrTax, comparison, tl, samps_of_int, OC, glomTab, mTab.sub.all, 
                   save.scats=save_plot, plot_dirs = plot_dir, facetRows = facetRows, forceFacet=forceFacet,
                   anova.pTab = mps, facetVar = facetVar, facetScales=facetScales, useLoess = useLoess, choose_yl=choose_yl,
                   specify_facet_labels=specify_facet_labels, use_ggtitle=use_ggtitle, use_xlab=use_xlab )
    
  } else if (plotType == "gradient") {
    
    ordObj <- subsampling_ordination_objects(ds.conts.all, phy.sub.all, pcoasOnly=T)
    
    plot_adonis_w_covars(dist_meas, phy.sub.all, mTab.sub.all, glomTab, 
                         shapeBy = divOrTax, colorBy = trait, tl = ifelse(tl=="contVar",NA,tlev), 
                         anova.pTab = mps, pcoaList = ordObj$pcoas, 
                         # log_of_color = ifelse(f %in% colorLogs, T, F), 
                         save.adoPlot = F, swapShapeColor=swapShapeColor)
    
  } else if (plotType == "assoc") {
    
    print(table(mTab.sub.all[,trait], mTab.sub.all[, divOrTax], dnn = c(trait, divOrTax)))
    print(table(mTab.sub.all[, divOrTax]))
    print(chisq.test(table(mTab.sub.all[,trait], mTab.sub.all[, divOrTax], dnn = c(trait, divOrTax))))
    
    assoc_label <- gsub("MALDI.","",gsub("Full_MALDI.","",divOrTax))
    # assoc_label <- sprintf("%s (p = %s)", assoc_label, signif(as.numeric(mps[divOrTax, trait]), digits = 3))
    
    vcd::assoc(table(mTab.sub.all[,trait], mTab.sub.all[, divOrTax], 
                     dnn = c(trait, assoc_label)), shade=T)
    
  }
  
  
}
# ************************************************************************************** #











# ****************************************************************************************************************** ####
# ****************************************************************************************************************** ####




# ***************************************************************************************** #
# Plot GRADIENTS of abundances of particular taxa or groups of taxa ####
library(RColorBrewer)

# ************************************* #
plot_gradients <- function(tl, tagl, dist_meas, pcoaList, glomTab, metaTab, org.title, 
                           shapes.numbers=FALSE, overrideOrg=NULL, alt.col.to.plot=NULL,
                           save.grads=F, plot_dirs=NULL, anova.pTab=NULL, drivers=NULL) {
  
  
  # ********************************* #
  # get PCoA object for the new distance matrix
  if ( ! is.null(pcoaList) ) {
    # use the precalculated pcoa coords if provided (useful if you want to the samples to stay in same coordinates)
    pcoaObj <- pcoaList[[ dist_meas ]]
    
    # ********************************* #
  } else {
    # otherwise caluclate the coordinates here using only non-NA samples for given trait
    
    # first keep only those samples for which colorBy is not NA
    if ( ! tl %in% c("Phylum","Class","Order","Family","Genus","Species")) {
      metaTab <- metaTab[ ! is.na(metaTab[,colorBy]), ]
      phy <- prune_samples(rownames(metaTab), phy)
    }
    
    # then get distance matrix using only samples without NAs for given variable
    if (dist_meas == "Weighted_Unifrac") dm <- UniFrac(phy, weighted = T, parallel = T)
    else if (dist_meas == "Unweighted_Unifrac") dm <- UniFrac(phy, weighted = F, parallel = T)
    else if (dist_meas == "VAW_GUnifrac") dm <- GUniFrac(t(phy@otu_table), phy@phy_tree)$unifracs[,,"d_VAW"]
    else if (dist_meas == "a05_GUnifrac") dm <- GUniFrac(t(phy@otu_table), phy@phy_tree)$unifracs[,,"d_0.5"]
    else if (dist_meas == "a0_GUnifrac") dm <- GUniFrac(t(phy@otu_table), phy@phy_tree)$unifracs[,,"d_0"]
    else if (dist_meas == "JSD") dm <- distance(phy, method = "jsd", parallel = T)
    else if (dist_meas == "Bray") dm <- vegdist(t(phy@otu_table), method = "bray")
    else if (dist_meas == "Jaccard") dm <- vegdist(decostand(as.data.frame(t(phy@otu_table)), method="pa"), method = "jaccard")
    else if (dist_meas == "Canberra") dm <- vegdist(t(phy@otu_table), method = "canberra")
    else if (dist_meas == "Aitchison") dm <- aDist( cmultRepl(t(phy@otu_table), method="CZM", label=0) )
    pcoaObj <- ape::pcoa(dm)
  }
  # ********************************* #
  
  # make matrix of plot values
  # coordTab <- as.data.frame(pcoaList[[ dist_meas ]]$vectors[ , 1:3 ])
  coordTab <- as.data.frame(pcoaObj$vectors[ , 1:3 ])
  
  if ("Rel_corr_eig" %in% names(pcoaObj$values)) {
    # varPers <- pcoaList[[ dist_meas ]]$values$Rel_corr_eig * 100
    varPers <- pcoaObj$values$Rel_corr_eig * 100
  } else if ("Relative_eig" %in% names(pcoaObj$values)) {
    # in case that this variable was not calculated (not totally sure why still, but for
    #   Jaccard it does not do this one) will just use relative eigenvalues without correction
    # varPers <- pcoaList[[ dist_meas ]]$values$Relative_eig * 100
    varPers <- pcoaObj$values$Relative_eig * 100
  } else {
    print("Rel_corr_eig and Relative_eig not calculated")
  }
  
  xl <- sprintf("Axis1 (%s%%)", round(varPers[1], 2))
  yl <- sprintf("Axis2 (%s%%)", round(varPers[2], 2))
  
  # ********************************* #
  
  
  
  
  if (tl == "Phylum")
    tl.mult <- "phyla"
  else if (tl == "Class")
    tl.mult <- "classes"
  else if (tl == "Order")
    tl.mult <- "orders"
  else if (tl == "Family")
    tl.mult <- "families"
  else if (tl == "Genus")
    tl.mult <- "genera"
  else if (tl == "Species")
    tl.mult <- "species"
  
  
  # ******************** #
  # get proper values and labels
  if (org.title %in% colnames(metaTab)) {
    col.to.plot <- sprintf("Stomatotype_%s", dist_meas)
    plot.ti <- sprintf("Gradient of %s || Dist = %s", org.title, dist_meas)
  } else if (org.title %in% rownames(glomTab[[ tl ]])) {
    org <- org.title
    plot.ti <- sprintf("%s: %s", tl, org.title)
    col.to.plot <- sprintf("Stomatotype_%s", dist_meas)
  } else if (org.title %in% c("Stom 1 sig","Stom 2 sig","Stom 3 sig","Stom 4 sig")) {
    sn <- strsplit(org.title, ' ')[[1]][2]
    org <- sig_tax.anova.dists[[ dist_meas ]][[ sn ]]
    plot.ti <- sprintf("%s significantly greater in Stomatotype %s", tl.mult, sn)
    col.to.plot <- sprintf("Stomatotype_%s", dist_meas)
  } else if (org.title == "Insig") {
    org <- rownames(glomTab[[ tl ]])[! rownames(glomTab[[ tl ]]) %in% unname(unlist(sig_tax.anova.dists[[ dist_meas ]]))]
    plot.ti <- sprintf("%s not significantly greater in any Stomatotype", tl.mult)
    col.to.plot <- sprintf("Stomatotype_%s", dist_meas)
  } else if (org.title %in% c("Stom_1_drivers", "Stom_2_drivers", "Stom_3_drivers")) {
    sn <- strsplit(org.title, '_')[[1]][2]
    org <- names(drivers[[ sn ]])
    plot.ti <- sprintf("top %s driver %s in Stomatotype %s", length(org), tl.mult, sn)
    col.to.plot <- sprintf("Stomatotype_%s", dist_meas)
  }
  
  # ******************** #
  if ( ! is.null(alt.col.to.plot) )
    col.to.plot <- alt.col.to.plot
  
  # ******************** #
  # org.abund <- as.numeric(as.matrix(colSums(glomTab[[ tl ]][org,])))
  # pcoa.org <- cbind(pcoas[[ dist_meas ]]$li, 
  #                   as.numeric(as.matrix(colSums(glomTab[[ tl ]][org,]))), 
  #                   as.factor(clu.gradient[["cluster_full"]]))
  # colnames(pcoa.org) <- c(colnames(pcoas[[ dist_meas ]]$li), "org", "Stomatotype")
  
  # ******************** #
  if (org.title %in% colnames(metaTab)) {
    org.abund <- metaTab[ sample_names(tagl), org.title ]
    pcoa.org <- as.data.frame(cbind(pcoaList[[ dist_meas ]]$vectors[ , 1:3], 
                                    org.abund,
                                    as.character(metaTab[ rownames(pcoaList[[ dist_meas ]]$vectors[ , 1:3]), 
                                                          col.to.plot ])))
    
    full.plot.ti <- plot.ti
    # if org.title is not a taxon, just use org.title in file name when saving
    org.title.fname <- org.title
    
  } else {
    
    if (length(org) == 1)
      org.abund <- as.matrix(as.data.frame(glomTab[[ tl ]]))[org, sample_names(tagl)]
    else
      org.abund <- colSums(glomTab[[ tl ]][org, sample_names(tagl)])
    
    pcoa.org <- as.data.frame(cbind(pcoaList[[ dist_meas ]]$vectors[ , 1:3], 
                                    org.abund,
                                    as.character(metaTab[ rownames(pcoaList[[ dist_meas ]]$vectors[ , 1:3]), 
                                                          col.to.plot ])))
    full.plot.ti <- sprintf("Gradient of abundances - %s\nDist = %s", plot.ti, dist_meas)
    # if org.title is a taxon, include the tax level in file name when saving
    org.title.fname <- sprintf("%s.%s", tl, org.title)
  }
  # ******************** #
  colnames(pcoa.org) <- c("A1","A2","A3", "org", "Stomatotype")
  # make sure values are correct class
  pcoa.org$A1 <- as.numeric(as.character(pcoa.org$A1))
  pcoa.org$A2 <- as.numeric(as.character(pcoa.org$A2))
  pcoa.org$A3 <- as.numeric(as.character(pcoa.org$A3))
  pcoa.org$org <- as.numeric(as.character(pcoa.org$org))
  pcoa.org$Stomatotype <- as.factor(pcoa.org$Stomatotype)
  
  
  bins <- quantile(org.abund, probs=seq(0,1,0.1), na.rm = T)
  
  pcoa.org$bin <- as.factor(sapply(pcoa.org$org, function(x) round(max(bins[ x >= bins ]),2)))
  # when a value become "-Inf", it must be the lowest bin, change to that
  pcoa.org$bin[ pcoa.org$bin=="-Inf" ] <- round(min(bins), 2)
  pcoa.org$bin <- droplevels(pcoa.org$bin)
  
  
  if (length(table(pcoa.org$bin)) < 6) {
    # means there was a limited distribution of values - the value for most of the quantiles was 0
    # in this case, double the number of quantiles, so have greater range of colors to better catch spread of values
    bins <- quantile(org.abund, probs=seq(0, 1, 0.05), na.rm = T)
    pcoa.org$bin <- as.factor(sapply(pcoa.org$org, function(x) round(max(bins[ x >= bins ]),2)))
    
    # try once more if necessary
    if (length(table(pcoa.org$bin)) < 6) {
      bins <- quantile(org.abund, probs=seq(0, 1, 0.01), na.rm = T)
      pcoa.org$bin <- as.factor(sapply(pcoa.org$org, function(x) round(max(bins[ x >= bins ]),2)))
    }
  }
  
  
  # remove rows where org is NA
  pcoa.org <- pcoa.org[ ! is.na(pcoa.org$org), ]
  
  
  
  co <- colorRampPalette(rev(brewer.pal(11,"RdYlBu")))(length(unique(pcoa.org$bin)))
  # percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
  
  if (shapes.numbers == TRUE) shape_vals <- c(49, 50, 51, 52, 53, 54, 55, 56)
  else shape_vals <- c(16, 12, 8, 18, 1, 2, 3, 4)
  
  shape_vals <- shape_vals[ 1:length(unique(metaTab[ rownames(pcoaList[[ dist_meas ]]$vectors[ , 1:3]), col.to.plot ])) ]
  
  # just checking to see which is the sample that is Stomatotype 1, and was far along the first component grouped
  #    with other S1 samples, yet had quite low abundance of the S1 genera (in Weighted_Unifrac) - doesnt match gradient
  # print(rownames(pcoa.org)[sapply(rownames(pcoa.org), function(x) pcoa.org[x,"Stomatotype"]=="1" && as.numeric(as.matrix(pcoa.org)[x, "bin"])<56 && as.numeric(as.matrix(pcoa.org)[x,"A1"])<(-0.1))])
  # print(rownames(pcoa.org)[sapply(rownames(pcoa.org), function(x) as.numeric(as.matrix(pcoa.org)[x, "bin"])<56 && as.numeric(as.matrix(pcoa.org)[x,"A1"])<(-0.4))])
  
  
  
  # ******************** #
  if ( ! is.null(anova.pTab) ) {
    
    # have to reverse the change of the characters that needed to be changed during lm tests so can check table
    if (org.title == "Absconditabacteriales_(SR1)")
      org.title.pTab <- "Absconditabacteriales_.SR1."
    else if ( ! org.title %in% colnames(metaTab))
      # only change this for taxa variables
      org.title.pTab <- gsub("-", "\\.", org.title)
    else
      org.title.pTab <- org.title
    
    
    # pval of trait of interest
    if (org.title.pTab %in% colnames(metaTab)) {
      ano.P <- anova.pTab[ col.to.plot, org.title.pTab ]
    } else {
      ano.P <- anova.pTab[ org.title.pTab, col.to.plot ]
    }
    ano.P.title <- sprintf("\nAnova pval = %s", ano.P)
    
    # pvals of covariates, if any signif (ignore seqGroup, as this was included just to account for its variation)
    if (org.title.pTab %in% colnames(metaTab)) {
      covsToCheck <- colnames(anova.pTab)[ ! colnames(anova.pTab) %in% c(org.title.pTab, "seqGroup") ]
      ano.covs.P <- anova.pTab[ col.to.plot, covsToCheck ]
    } else {
      covsToCheck <- colnames(anova.pTab)[ ! colnames(anova.pTab) %in% c(col.to.plot, "seqGroup") ]
      ano.covs.P <- anova.pTab[ org.title.pTab, covsToCheck ]
    }
    # in case only 1 cov, will lose the name, must ensure it keeps it
    names(ano.covs.P) <- covsToCheck
    
    # take only those covariates which were significant
    acp <- ano.covs.P[ ano.covs.P != "" ]
    # make string for title
    ano.covs.P.title <- ifelse(length(acp) == 0,
                               "\nOther effects:   none",
                               sprintf("\nOther effects:   %s", 
                                       paste(paste(names(acp), acp, sep="="), collapse=" || ")))
    
    # finally, add these values to plot title
    full.plot.ti <- sprintf("%s%s%s", full.plot.ti, ano.P.title, ano.covs.P.title)
  }
  # ******************** #
  
  
  
  # ******************** #
  grads <- ggplot(pcoa.org, aes(x=A1, y=A2), col = bin) +
    geom_point(aes(col=bin, shape=Stomatotype), size=4) +
    theme_minimal() +
    scale_color_manual(values=co) +
    scale_shape_manual(values=shape_vals) +
    labs(shape=ifelse(startsWith(col.to.plot,"Stomatotype"),"Stomatotype",col.to.plot), col=org.title) +
    xlab(xl) + ylab(yl) +
    theme(plot.title = element_text(size=17), legend.text = element_text(size=13),
          legend.title = element_text(size = 15, face = "bold"),
          axis.title = element_text(size=15), axis.text = element_text(size=15)) +
    ggtitle(full.plot.ti)
  
  # ******************** #
  if (save.grads == T) {
    
    # make directory if necessary
    dir.create(sprintf("%s/figures/%s", p2_dir, plot_dirs), showWarnings = F)
    
    ggsave(sprintf("%s/figures/%s/%s-%s.Gradients.png", p2_dir, plot_dirs, dist_meas, org.title.fname),
           width = 9.23, height = 7.04, device = "png", plot = grads)
    
  } else {
    grads
  }
  # ******************** #
}
# ************************************* #













# ****************************************************************************************************************** #
# Heatmap of correlations for continuous/binary responses ####
# ****************************************************************************************************************** #


# determine which samples to look at
samps.of.int <- "All"

interest <- sample_names(SLL2)#names(cluster_vals)



# ********************************************************************************* #
fill_cor_tables <- function(comparison, tl, interest, contQs, cor.method="pearson") {
  
  # first get appropriate data
  if (comparison == "TvsQ") {
    # o.tmp <- gloms_rel[[ tl ]][ , gsub('\\.', '-', interest) ]
    # colnames(o.tmp) <- gsub('-', '\\.', colnames(o.tmp))
    # data.cont <- cbind(t(o.tmp), SLL2@sam_data[interest, contQs])
    data.cont <- cbind(t(gloms_rel[[ tl ]][ , interest ]), SLL2@sam_data[interest, contQs])
    cor.cols <- rownames(gloms_rel[[ tl ]])
    cor.rows <- contQs
    
  } else if (comparison == "taxa") {
    # data.cont <- t(gloms_rel[[ tl ]][ , interest ])
    # cor.rows <- cor.cols <- rownames(gloms_rel[[ tl ]])
    data.cont <- t(rbind(gloms_rel$Species, gloms_rel$Genus,
                         gloms_rel$Family, gloms_rel$Order,
                         gloms_rel$Class, gloms_rel$Phylum)[ , interest ])
    # include taxonomic level in column names for those cases where names are shared bt levels, and for clarity
    colnames(data.cont) <- paste0(unname(unlist(sapply(c("Species","Genus","Family","Order","Class","Phylum"), 
                                                       function(x) rep(sprintf("%s..",x), nrow(gloms_rel[[ x ]]))))),
                                  colnames(data.cont))
    cor.rows <- cor.cols <- colnames(data.cont)
    
  } else if (comparison == "questions") {
    data.cont <- SLL2@sam_data[interest, contQs]
    cor.rows <- cor.cols <- contQs
  }
  
  data.cont <- apply(data.cont, 2, as.numeric)
  
  # prepare matrices for values
  cor.matrix <- matrix(NA, nrow=length(cor.rows), ncol=length(cor.cols))
  rownames(cor.matrix) <- cor.rows
  colnames(cor.matrix) <- cor.cols
  ps.matrix <- matrix(NA, nrow=length(cor.rows), ncol=length(cor.cols))
  rownames(ps.matrix) <- cor.rows
  colnames(ps.matrix) <- cor.cols
  
  # fill matrices
  for (i in cor.rows) {
    for (j in cor.cols) {
      correl <- cor.test(data.cont[,j], data.cont[,i], na.rm=T, method=cor.method)
      cor.matrix[i,j] <- correl$estimate
      ps.matrix[i,j] <- correl$p.value
    }
  }
  
  return(list("cor"=cor.matrix, "ps"=ps.matrix))
}

# ********************************************************************************* #

adjust_cor_p_vals <- function(cor_table, p_table, comparison) {
  # get table of adjusted p-values
  p.adj <- apply(p_table, 2, p.adjust, method='bonferroni')
  
  # make p-values for diagonals = 1 when appropriate because we dont care about self-correlations,
  # and cause this will ruin the color scale of the heatmaps
  if (comparison != "TvsQ") diag( p.adj ) <- 1
  
  # for those values that are given as NaN because the given genus is not present in 
  # any samples with particular categories of the given
  cor_table[ is.nan(cor_table) ] <- 0
  p.adj[ is.nan(p.adj) ] <- 1
  
  #For some questions that have all 0s (occurs when doing cities alone)
  cor_table[ is.na(cor_table) ] <- 0
  p.adj[ is.na(p.adj) ] <- 1
  
  #at least 1 good p value within samples
  mins <- apply(p.adj, 2, min)
  goodPs <- mins[mins < 0.05] #-log10(0.05) == 1.3013 
  #at least 1 good p within taxa
  tmins <- apply(p.adj, 1, min)
  tgoodPs <- tmins[tmins < 0.05]
  
  if (length(tgoodPs) == 1) {
    cor_table <-  t(as.matrix(cor_table[ names(tgoodPs), names(goodPs) ]))
    p.adj <-  t(as.matrix(p.adj[ names(tgoodPs), names(goodPs) ]))
    
    rownames(cor_table) <- rownames(p.adj) <- names(tgoodPs)
    
  } else if (length(goodPs) == 1) {
    print(c(comparison, tgoodPs, goodPs))
    # cor_table <-  t(as.matrix(cor_table[ names(tgoodPs), names(goodPs) ]))
    # p.adj <-  t(as.matrix(p.adj[ names(tgoodPs), names(goodPs) ]))
    # 
    # rownames(cor_table) <- rownames(p.adj) <- tgoodPs
  } else {
    cor_table <- cor_table[ names(tgoodPs), names(goodPs) ]
    p.adj <- p.adj[ names(tgoodPs), names(goodPs) ]
  }
  
  # # get tables of only significant values
  # cor_table[ p.adj >= 0.05 ] <- ''
  # p.adj[ p.adj >= 0.05 ] <- ''
  
  return(list("cor.adj"=cor_table, "p.adj"=p.adj))
}
# ********************************************************************************* #


# ********************************************************************************* #
plot_cor_heatmap <- function(cor.tab, ps.tab, comparison, tl) {
  
  # get table of pluses for significant correlations so can print in heatmap later
  pluses <- matrix('', nrow=nrow(ps.tab), ncol=ncol(ps.tab))
  colnames(pluses) <- colnames(ps.tab)
  rownames(pluses) <- rownames(ps.tab)
  pluses[ ps.tab < 0.05 ] <- '+'
  
  # use hclust to cluster the rows and columns based on the values
  if (nrow(cor.tab) <= 2 | ncol(cor.tab) <= 2) {
    # because cannot cluster with 2 or fewer values
    cor.m <- cbind( reshape2::melt( cor.tab ), reshape2::melt( pluses ) )
    
  } else {
    ord <- hclust( dist( cor.tab, method="euclidean" ), method = "ward.D" )$order
    ord2 <- hclust( dist( t(cor.tab), method="euclidean" ), method = "ward.D" )$order
    # get melted object, 
    cor.m <- cbind( reshape2::melt( cor.tab[ ord, ord2 ] ), reshape2::melt( pluses[ ord, ord2 ] ) )
  }
  
  cor.m <- cor.m[,c(1,2,3,6)]
  colnames(cor.m)[4] <- 'signif'
  
  # then subset it to color only the lower right triangle of the heatmap 
  # (since both sides of the diagonal will be the same)
  if (comparison == "TvsQ") {
    cor.m.lower <- cor.m ## not using lower triangle here since cols dont match rows as in the other "vs" instances
    ti <- sprintf("Heatmap of correlations: questions vs %s", tl)
    file.to.save <- sprintf("%s/figures/Correlations/%s/cor_heatmap.%s.%s.png",
                            p2_dir, samps.of.int, comparison, tl)
  } else {
    cor.m.lower <- subset( cor.m[ upper.tri(cor.tab[ ord, ord2 ]), ] )
    if (comparison=="questions") {
      ti <- sprintf("Heatmap of correlations: questions vs questions")
      file.to.save <- sprintf("%s/figures/Correlations/%s/cor_heatmap.%s.png",
                              p2_dir, samps.of.int, comparison)
    }
    if (comparison=="taxa") {
      ti <- sprintf("Heatmap of correlations: %s vs %s", tl, tl)
      file.to.save <- sprintf("%s/figures/Correlations/%s/cor_heatmap.%s.%s.png",
                              p2_dir, samps.of.int, comparison, tl)
    }
  }
  
  # plot and save heatmap
  p <- ggplot(data = cor.m, aes(x=Var2, y=Var1, fill=value)) + 
    geom_tile(color="white", data=cor.m.lower) + 
    geom_text(aes(label=signif), size = 4.5, data=cor.m.lower) + 
    scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0, space="Lab",name="Pearson") +
    theme_minimal() + 
    # theme(axis.text = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 8, hjust = 1),
          axis.text.y = element_text(angle = 0, vjust = 1, size = 10, hjust = 1),
          plot.title = element_text(hjust=0.5)) +
    xlab('') + ylab('') + ggtitle(ti)
  
  ggsave(filename = file.to.save, plot = p, device = "png")
}
# ********************************************************************************* #















# ****************************************************************************************************************** #
# Scatterplot to observe correlations ####
# ************************** #
get_data.cont <- function(comparison, tl, interest, contQs, glomTab, mTab, facetVar=NULL) {
  
  # first get appropriate data
  if (comparison == "TvsQ") {
    data.cont <- cbind(t(glomTab[[ tl ]][ , interest ]), as.matrix(mTab)[interest, contQs])
    cor.cols <- rownames(glomTab[[ tl ]])
    cor.rows <- contQs[ sapply(contQs, function(x) sum( ! is.na(data.cont[, x])) != 0) ]
    
  } else if (comparison == "taxa") {
    data.cont <- t(glomTab[[ tl ]][ , interest ])
    cor.rows <- cor.cols <- rownames(glomTab[[ tl ]])
    
  } else if (comparison == "questions") {
    data.cont <- as.matrix(mTab)[interest, contQs]
    cor.rows <- cor.cols <- contQs[ sapply(contQs, function(x) sum( ! is.na(data.cont[, x])) != 0) ]
    
  }
  data.cont <- as.data.frame(data.cont[ , sapply(colnames(data.cont), function(x) sum( ! is.na(data.cont[, x])) != 0) ])
  data.cont <- apply(data.cont, 2, as.numeric)
  rownames(data.cont) <- interest
  
  # add column to facet variables if desired
  if ( ! is.null(facetVar)) {
    if (grepl("~", facetVar)) {
      # if faceting by combination of 2 variables
      fv1 <- strsplit(facetVar, "~")[[1]][1]
      fv2 <- strsplit(facetVar, "~")[[1]][2]
      data.cont <- as.data.frame(cbind(data.cont, as.matrix(mTab)[interest, c(fv1, fv2)]))
      colnames(data.cont)[ c(ncol(data.cont)-1, ncol(data.cont)) ] <- c(fv1, fv2)
      
    } else {
      # if faceting by just 1 variable
      data.cont <- as.data.frame(cbind(data.cont, as.matrix(mTab)[interest, facetVar]))
      colnames(data.cont)[ ncol(data.cont) ] <- facetVar
    }
    
  }
  
  return(data.cont)
}
# ************************************* #
library(ggnewscale)
# ************************************* #
plot_data.cont <- function(n1, n2, comparison, tl, interest, contQs, glomTab, mTab, useLoess=F, forceFacet=F,
                           facetVar=NULL, save.scats=F, plot_dirs=NULL, anova.pTab=NULL, facetRows=NULL,
                           specify_facet_labels=NULL, facetScales="free", use_ggtitle=T, use_xlab=T, choose_yl=NULL) {
  
  # first get table of values
  dc <- get_data.cont(comparison, tl, interest, contQs, glomTab, mTab, facetVar)
  
  
  # ******************** #
  if (length(n2) > 1)
    n2.ti <- ifelse(tl != "contVar" & sum(n2 %in% rownames(glomTab[[ tl ]]))==length(n2), 
                    sprintf("Indicated %s", 
                            gsub("Phylum","Phyla", 
                                 gsub("Genus","Genera", tl))), "Indicated variables")
  else if (n2 %in% contQs)
    n2.ti <- n2
  else
    # if n2 is a taxon, include the tax level in file name when saving
    n2.ti <- sprintf("%s: %s", tl, ifelse(n2=="unclassified.P1", "unclassified.Phylum", n2))
  # ******************** #
  
  # ******************** #
  if ( ! is.null(anova.pTab) & length(n2)==1 ) {
    
    # have to reverse the change of the characters that needed to be changed during lm tests so can check table
    if (n2 == "Absconditabacteriales_(SR1)")
      n2.pTab <- "Absconditabacteriales_.SR1."
    else
      n2.pTab <- gsub("-", "\\.", n2)
    
    # pval of trait of interest
    ano.P <- anova.pTab[ n2.pTab, n1 ]
    ano.P.title <- sprintf("\nAnova pval = %s", signif(as.numeric(ano.P), digits = 3))
    
    # pvals of covariates, if any signif (ignore seqGroup, as this was included just to account for its variation)
    covsToCheck <- colnames(anova.pTab)[ ! colnames(anova.pTab) %in% c(n1, "seqGroup") ]
    ano.covs.P <- anova.pTab[ n2.pTab, covsToCheck ]
    # in case only 1 cov, will lose the name, must ensure it keeps it
    names(ano.covs.P) <- covsToCheck
    
    # take only those covariates which were significant
    acp <- ano.covs.P[ ano.covs.P != "" ]
    # make string for title
    ano.covs.P.title <- ifelse(length(acp) == 0,
                               "\nOther effects:   none",
                               sprintf("\nOther effects:   %s", 
                                       paste(paste(names(acp), signif(as.numeric(acp), digits=3), 
                                                   sep="="), collapse=" || ")))
    
  } else {
    ano.P.title <- ""
    ano.covs.P.title <- ""
  }
  # ******************** #
  
  
  
  
  # plot figure
  # ******************************************************************* #
  # ******************************************************************* #
  # ******************************************************************* #
  
  if (length(n2) > 1 | forceFacet==T) {
    
    dat <- data.frame(dc[,n1], dc[,n2])
    
    colnames(dat) <- c("n1",n2)
    dat <- dat[ rowSums( is.na(dat) ) == 0, ]
    dat <- melt(dat, id.vars="n1")
    colnames(dat) <- c("n1","variable","value")
    
    # ******************** #
    if ( ! is.null(anova.pTab) ) {
      ano.P <- as.numeric(anova.pTab[ n2, n1 ])
      names(ano.P) <- n2
      
      facet_labels <- sapply(n2, function(x) {
        head.anoP <- sprintf("%s (pval %s)", x, 
                             # in case showing a covar, pval may not have been signif, will indicate that
                             ifelse(is.na(ano.P[x]), "> 0.05", 
                                    sprintf("= %s", signif(as.numeric(ano.P[x]), digits = 3))))
        
        # pvals of covariates, if any signif (ignore seqGroup, as this was included just to account for its variation)
        covsToCheck <- colnames(anova.pTab)[ ! colnames(anova.pTab) %in% c(n1, "seqGroup") ]
        ano.covs.P <- anova.pTab[ x, covsToCheck ]
        # in case only 1 cov, will lose the name, must ensure it keeps it
        names(ano.covs.P) <- covsToCheck
        # take only those covariates which were significant
        acp <- ano.covs.P[ ano.covs.P != "" ]
        # make string for title
        ano.covs.P.title <- ifelse(length(acp) == 0,
                                   "\nOther effects:   none",
                                   sprintf("\nOther effects:   %s", 
                                           paste(paste(names(acp), signif(as.numeric(acp), digits=3), 
                                                       sep="="), collapse=" || ")))
        return(sprintf("%s%s", head.anoP, ano.covs.P.title))
      })
    } else {
      facet_labels <- sapply(n2, function(x) x)
    }
    
    # rename unclassified.G43 as a more understandable name
    facet_labels <- gsub("unclassified.G43","unclass.Phy", 
                         gsub("unclassified.P1","unclassified.Phylum", facet_labels))
    # ******************** #
    if ( ! is.null(specify_facet_labels)) {
      facet_labels <- specify_facet_labels
      names(facet_labels) <- n2
    }
    # ******************** #
    
    scatPlot <- ggplot(dat, aes(x=n1, y=value)) +
      geom_point(shape=1)
      # geom_smooth() +
      # geom_smooth(method = lm) +
    
    if ( ! is.null(facetRows)) {
      scatPlot <- scatPlot +
        facet_wrap(~variable, scales=facetScales, labeller = as_labeller(facet_labels),
                   nrow = facetRows)
    } else {
      scatPlot <- scatPlot +
        facet_wrap(~variable, scales=facetScales, labeller = as_labeller(facet_labels))
    }
      
    scatPlot <- scatPlot +
      # ggtitle(sprintf("%s vs %s", n1, n2.ti)) +
      theme_minimal() +
      theme(axis.text = element_text(size=18), axis.title = element_text(size=20),
            plot.title = element_text(size=20), strip.text = element_text(size=14),
            strip.background = element_rect(fill="white")) +
      xlab(n1) + ylab(ifelse( ! is.null(choose_yl), choose_yl,
                              ifelse(tl=="contVar","","Abundance")))
    
    if (use_ggtitle)
      scatPlot <- scatPlot + ggtitle(sprintf("%s vs %s", n1, n2.ti))
    
    if (use_xlab)
      scatPlot <- scatPlot + xlab(n1)
    else
      scatPlot <- scatPlot + xlab(NULL)
    
    if (useLoess == T) {
      scatPlot <- scatPlot + geom_smooth(method = loess)
    } else {
      scatPlot <- scatPlot + geom_smooth(method = lm)
    }
    
    facet.fname <- ""
    
    # ******************************************************************* #
    # ******************************************************************* #
    
  } else {
    
    # ******************************************************************* #
    # ******************************************************************* #
    
    if ( ! is.null(facetVar)) {
      
      if (grepl("~", facetVar)) {
        # if faceting by combination of 2 variables
        fv1 <- strsplit(facetVar, "~")[[1]][1]
        fv2 <- strsplit(facetVar, "~")[[1]][2]
        
        dat <- data.frame(as.numeric(as.character(dc[,n1])), as.numeric(as.character(dc[,n2])), dc[,c(fv1, fv2)])
        colnames(dat) <- c("x","y","fv1","fv2")
        dat <- dat[!is.na(dat$x) & !is.na(dat$y) & !is.na(dat$fv1) & !is.na(dat$fv2),]
        
        scatPlot <- ggplot(dat, aes(x=x, y=y)) +
          geom_point(shape=1) + 
          # geom_smooth() +
          facet_wrap(fv1~fv2) +
          geom_smooth(method = lm) +
          ggtitle(sprintf("%s vs %s%s%s\nfacets: %s", n1, n2.ti, ano.P.title, ano.covs.P.title, facetVar)) +
          theme_minimal() +
          theme(axis.text = element_text(size=17), axis.title = element_text(size=18),
                plot.title = element_text(size=17), strip.text = element_text(size=16),
                strip.background = element_rect(fill="white")) +
          xlab(n1) + ylab(n2)
        
      } else {
        # if faceting by just 1 variable
        dat <- data.frame(as.numeric(as.character(dc[,n1])), as.numeric(as.character(dc[,n2])), dc[,facetVar])
        colnames(dat) <- c("x","y","facetVar")
        dat <- dat[!is.na(dat$x) & !is.na(dat$y) & !is.na(dat$facetVar),]
        
        if (facetVar %in% c("Downs_Syndrome","Cystic_fibrosis")) {
          # in order to use the same colors for disorders as in the age plot in the section:
          #   " # Scatterplot of healthyPop and DS age_diffs vs distances, combined "
          
          if (facetVar == "Downs_Syndrome") {
            dot_cols <- c("#CB6767","#73BDD3")
            lin_cols <- c("#bf4342","#2f6690")
          } else if (facetVar == "Cystic_fibrosis") {
            # dot_cols <- c("#e9c46a","#73BDD3")
            # lin_cols <- c("#B1871B","#2f6690")
            # dot_cols <- c("#f77f00","#73BDD3")
            # lin_cols <- c("#A35400","#2f6690")
            # dot_cols <- c("#fca311","#73BDD3")
            # lin_cols <- c("#A26402","#2f6690")
            dot_cols <- c("#E4BA4E","#73BDD3")
            lin_cols <- c("#A26402","#2f6690")
          }
          # scatPlot <- ggplot(dat, aes(x=x, y=y, color=facetVar)) +
          #   geom_point(shape=1) + 
          #   # geom_smooth() +
          #   facet_wrap(~facetVar) +
          #   geom_smooth(method = lm) +
          #   ggtitle(sprintf("%s vs %s%s%s\nfacets: %s", n1, n2.ti, ano.P.title, ano.covs.P.title, facetVar)) +
          #   theme_bw() +
          #   theme(axis.text = element_text(size=17), axis.title = element_text(size=18),
          #         plot.title = element_text(size=17), strip.text = element_text(size=16)) +
          #   xlab(n1) + ylab(ifelse(n2=="unclassified.P1", "unclassified.Phylum", n2)) +
          #   scale_fill_manual(values=c("#73BDD3","#CB6767"))
          
          scatPlot <- ggplot(dat, aes(x=x, y=y, color=facetVar)) +
            geom_point(aes(shape=facetVar, color=facetVar), size=2.5) + 
            scale_shape_manual(values=rev(c(16,1)), guide=F) +
            scale_color_manual(values=rev(dot_cols)) +
            labs(x="Age", y=ifelse(n2=="unclassified.P1", "unclassified.Phylum", n2), 
                 color=ifelse(facetVar=="Downs_Syndrome", "DS", 
                              ifelse(facetVar=="Cystic_fibrosis","CF", facetVar))) +
            new_scale("color") +
            geom_smooth(method = lm, aes(color=facetVar)) +
            scale_color_manual(values=rev(lin_cols), guide=F) +
            ggtitle(sprintf("%s vs %s%s%s\nfacets: %s", n1, n2.ti, ano.P.title, ano.covs.P.title, facetVar)) +
            theme_minimal() +
            theme(axis.text = element_text(size=17), axis.title = element_text(size=18),
                  plot.title = element_text(size=17), 
                  legend.title = element_text(size=17), legend.text = element_text(size=16))
          
          # ggplot(dda, aes(x=age_diffs, y=dists, color=samps)) +
          #   geom_point(aes(shape = samps, color=samps), size=2.5) +
          #   scale_shape_manual(values=c(16,1), guide=F) +
          #   scale_color_manual(values=c("#CB6767","#73BDD3")) +
          #   labs(x="Age Difference", y="Aitchison Distance", color="") +
          #   new_scale("color") +
          #   geom_smooth(method = lm, aes(color=samps)) +
          #   scale_color_manual(values=c("#bf4342","#2f6690"), guide=F) +
          #   theme_classic() +
          #   theme(axis.title = element_text(size=17), axis.text = element_text(size=15),
          #         legend.title = element_text(size=17), legend.text = element_text(size=15))
          
          
        } else {
          
          scatPlot <- ggplot(dat, aes(x=x, y=y)) +
            geom_point(shape=1) + 
            # geom_smooth() +
            facet_wrap(~facetVar) +
            geom_smooth(method = lm) +
            ggtitle(sprintf("%s vs %s%s%s\nfacets: %s", n1, n2.ti, ano.P.title, ano.covs.P.title, facetVar)) +
            theme_minimal() +
            theme(axis.text = element_text(size=17), axis.title = element_text(size=18),
                  plot.title = element_text(size=17), strip.text = element_text(size=16),
                  strip.background = element_rect(fill="white")) +
            xlab(n1) + ylab(n2)
          
        }
        
      }
      
      
      
      facet.fname <- sprintf("-facet.%s", facetVar)
      
    } else {
      
      dat <- data.frame(dc[,n1], dc[,n2])
      colnames(dat) <- c("x","y")
      dat <- dat[!is.na(dat$x) & !is.na(dat$y),]
      
      scatPlot <- ggplot(dat, aes(x=x, y=y)) +
        geom_point(shape=1) + 
        # geom_smooth() +
        geom_smooth(method = lm) +
        ggtitle(sprintf("%s vs %s%s%s", n1, n2.ti, ano.P.title, ano.covs.P.title)) +
        theme_minimal() +
        theme(axis.text = element_text(size=17), axis.title = element_text(size=18),
              plot.title = element_text(size=17)) +
        xlab(n1) + ylab(n2)
      
      facet.fname <- ""
    }
    
  }
  
  # ******************************************************************* #
  # ******************************************************************* #
  # ******************************************************************* #
  
  
  
  # ******************** #
  if (save.scats == T) {
    
    if (n2 %in% contQs)
      n2.fname <- n2
    else
      # if n2 is a taxon, include the tax level in file name when saving
      n2.fname <- sprintf("%s.%s", tl, n2)
    
    # make directory if necessary
    dir.create(sprintf("%s/figures/%s", p2_dir, plot_dirs), showWarnings = F)
    
    ggsave(sprintf("%s/figures/%s/%s-%s%s.png", p2_dir, plot_dirs, n1, n2.fname, facet.fname),
           width = 9.23, height = 7.04, device = "png", plot = scatPlot)
    
  } else {
    scatPlot
  }
  # ******************** #
  
}
# ************************************* #
# ************************************* #

















# ****************************************************************************************************************** #
# Kruskal-Wallis tests ####
# ****************************************************************************************************************** #

#kruskal-wallis to test for any separations within a given survey response

# ********************************************************************************* #
fill_kw_p_table <- function(phy, cont, interest, group_qs, contQs, rarefyData=F, externalGroupQ=NULL) {
  # first get appropriate data
  if (! is.null(externalGroupQ)) {
    gq <- as.matrix(externalGroupQ[ interest ])
    colnames(gq) <- "externalGroupQ"
  } else {
    gq <- as.matrix(phy@sam_data[ interest, group_qs ])
  }
  
  
  if (cont == "cont_vars") {
    data.cont <- phy@sam_data[ interest, contQs ]
  } else if (cont %in% c("Phylum","Class","Order","Family","Genus","Species")) {
    if (rarefyData==TRUE) {
      gr <- get_gloms(phy, cont, taxTables.both)
    } else {
      gr <- gloms_rel[[ cont ]]
    }
    data.cont <- t(gr[ , interest ])
    rownames(data.cont) <- rownames(data.cont)
  }
  data.cont[data.cont == "No Sabe/No Contesta"] <- NA
  data.cont <- as.data.frame(apply(data.cont, 2, factor))
  
  
  # prepare matrices for values, then fill them
  kw.p <- matrix(NA, nrow=ncol(data.cont), ncol=ncol(gq))
  colnames(kw.p) <- colnames(gq)
  rownames(kw.p) <- colnames(data.cont)
  
  #store kw$statistic, kw$p.value
  for (i in colnames(data.cont)) {
    # print(i)
    for (j in colnames(gq)) {
      # print(c(i,j))
      if (length(table(gq[, j], useNA = "no")) == 1) {
        kw.p[i,j] <- 1
      } else if (length(table(gq[, j], useNA = "no")) == 2 & min(table(gq[, j], useNA = "no")) == 1) {
        kw.p[i,j] <- 1
      } else if ( sum(rowSums(table(gq[,j], data.cont[,i], useNA = "no")) != 0) < 2) {
        # in the case that "all observations are in the same group"
        # length(table(gq[,j])) == 2 & 
        #          0 %in% rowSums(table(gq[,j], data.cont[,i]))) {
        # print(c(i,j))
        # where there is only a value for the cases, none for cont
        kw.p[i,j] <- 1
      } else {
        kw <- kruskal.test(as.numeric(as.matrix(data.cont)[,i]), as.factor(gq[,j]), na.rm=T )#na.action='na.exclude'
        kw.p[i,j] <- kw$p.value
      }
    }
  }
  
  return(kw.p)
  
}

# ********************************************************************************* #

adjust_kw_p_vals <- function(p_table) {
  # get table of adjusted p-values
  p.adj <- apply(p_table, 2, p.adjust, method='fdr')
  
  # for those values that are given as NaN because the given genus is not present in 
  # any samples with particular categories of the given
  p.adj[ is.nan(p.adj) ] <- 1
  
  #at least 1 good p value within samples
  mins <- apply(p.adj, 2, min)
  goodPs <- mins[mins < 0.05] #-log10(0.05) == 1.3013 
  #at least 1 good p within taxa
  tmins <- apply(p.adj, 1, min)
  tgoodPs <- tmins[tmins < 0.05]
  
  p.adj <- p.adj[ names(tgoodPs), names(goodPs) ]
  
  # get tables of only significant values
  p.adj[ p.adj >= 0.05 ] <- ''
  
  if (class(p.adj)=="character" & length(p.adj)==1) {
    p.adj <- as.data.frame(p.adj)
    colnames(p.adj) <- names(goodPs)
    rownames(p.adj) <- names(tgoodPs)
  }
  
  return(as.matrix(p.adj))
}

# ********************************************************************************* #




# ******************** #
# boxplots of particular group_q vs cont data/otu abundance ####

library(ggpubr)
# ******************** #
group_vs_cont_box <- function(mTab, interest, tlev, o_rel, group_col, cont_col, GQ, CONTS, 
                              plot_tukey=T, print_tukey=T, useNotch=T, xAngle=0, facet=NULL, adjustAlphas=F,
                              externalGroupQ=NULL, save.boxes=F, plot_dirs=NULL, anova.pTab=NULL, flipCoords=F,
                              atex.s=17, atit.s=18, pts=17, sts=13, facetRows=NULL, facetScales="free",
                              choose_xl=NULL, choose_yl=NULL, choose_title=NULL,
                              include_legend_below=F, str.bck="grey90") {
  
  # for communities, shorten some names
  if (group_col == "Community") {
    mTab[ , group_col] <- gsub("Comunidad de Madrid","Madrid",
                               gsub("Comunidad Valenciana","Valencia",
                                    gsub("Comunidad Foral de Navarra","Navarra",
                                         gsub("Región de Murcia","Murcia",
                                              gsub("Principado de Asturias","Asturias",
                                                   mTab[ , group_col])))))
  }
  
  
  if (! is.null(externalGroupQ)) {
    gqs <- as.matrix(externalGroupQ[ interest ])
    colnames(gqs) <- group_col
  } else {
    gqs <- as.matrix(mTab[ interest, GQ])
  }
  
  # in case treating certain numerical variables as the group variable here, may have to remove 
  #   leading white space created after converting to characters with "as.matrix" above
  gqs[ , group_col] <- trimws(gqs[ , group_col], which = "left")
  
  
  data.mix <- cbind(t(o_rel[ , interest]), mTab[ interest, CONTS])
  data.mix <- as.data.frame(apply(data.mix, 2, factor))
  
  kw.box <- as.data.frame( cbind(as.matrix(data.mix[, cont_col]), gqs[, group_col]) )
  
  if (length(cont_col) > 1) {
    
    colnames(kw.box) <- c(cont_col,"group")
    kw.box <- melt(kw.box, id.vars="group")
    colnames(kw.box) <- c("group","variable","value")
    kw.box[ kw.box == "No Sabe/No Contesta" ] <- NA
    kw.box <- kw.box[ ! is.na(kw.box$group), ]
  } else {
    colnames(kw.box) <- c("cont","group")
    kw.box[ kw.box == "No Sabe/No Contesta" ] <- NA
    kw.box <- kw.box[ ! is.na(kw.box$group), ]
  }
  
  
  # ******************** #
  removed0s <- ""
  if (plot_tukey == T | print_tukey == T) {
    # Check signif difs by anova 
    if (! is.null(externalGroupQ)) {
      adf <- cbind(t(o_rel[ , interest ]), mTab[ interest, ], gqs[interest,])
      colnames(adf) <- c(rownames(o_rel), colnames(mTab), group_col)
      adf <- as.data.frame(adf)
      adf[,group_col] <- as.factor(adf[,group_col])
    } else {
      adf <- cbind(t(o_rel[ , interest ]), mTab[ interest, ])
    }
    
    # for pH, since it is originally numeric, but here being treated as factor, 5.5, 6.5 and 7.5 are unchanged
    #   but 5.0, 6.0, 7.0 and 8.0 become 5, 6, 7, and 8, in the Tukeytable
    #   this causes problem in plot where ggplot is using the character values but the tukey table was not, 
    # so here I force the whole values to end in ".0"
    if (group_col == "pH") {
      adf[,group_col] <- as.character(adf[,group_col])
      adf[,group_col] <- sapply(adf[,group_col], function(x) ifelse(nchar(x)==1, sprintf("%s.0", x), x))
    }
    
    adf[ adf == "No Sabe/No Contesta" ] <- NA
    # adf <- adf[ ! is.na(adf[, cont_col]), ] # remove NAs from 
    res.aov <- aov(formula = as.numeric(as.matrix(adf[, cont_col])) ~ as.factor(as.matrix(adf[, group_col])), data = adf)
    # summary(res.aov)
    # print( TukeyHSD(res.aov) )
    if (print_tukey==T) {
      cat("Check signif difs by anova:\n\n")
    }
    if (class(TukeyHSD(res.aov)[[1]][ TukeyHSD(res.aov)[[1]][ , "p adj"] < 0.05, ]) == "numeric") {
      # ******************** #
      # if only one row has signif dif, will need to print rowname separately
      rn <- rownames(TukeyHSD(res.aov)[[1]])[ TukeyHSD(res.aov)[[1]][ , "p adj"] < 0.05 ]
      # prepare table for plotting pvals with stat_pvalue_manual
      tukPval <- TukeyHSD(res.aov)[[1]][ TukeyHSD(res.aov)[[1]][ , "p adj"] < 0.05, "p adj"]
      cont.max <- max(adf[,cont_col], na.rm = T)
      # arbitrarily make the space between lines based on the range of values
      rango <- (max(adf[,cont_col], na.rm = T)-min(adf[,cont_col], na.rm = T)) / 10
      cont.max <- cont.max + rango
      stat.test <- data.frame(".y."="cont", 
                              "group1"=strsplit(rn, "-")[[1]][1], 
                              "group2"=strsplit(rn, "-")[[1]][2],
                              "p"=NA,
                              "p.adj"=formatC(tukPval, format="e", digits=3),
                              "p.format"=formatC(tukPval, format="e", digits=3),
                              "p.signif"=symnum(tukPval, 
                                                cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                                symbols = c("****", "***", "**", "*", "ns")),
                              "method"="Tukey",
                              "y.position"=cont.max)
      if (print_tukey==T)
        print(c( rn, TukeyHSD(res.aov)[[1]][ TukeyHSD(res.aov)[[1]][ , "p adj"] < 0.05, ] ))
      
    } else {
      # ******************** #
      # prepare table for plotting pvals with stat_pvalue_manual
      TukeyTab <- TukeyHSD(res.aov)[[1]][ TukeyHSD(res.aov)[[1]][ , "p adj"] < 0.05, ]
      
      # # ******************** #
      # # for community and province, there are sometimes very strong differences between many regions
      # # so here will remove the pvalue comparison lines in plot when p=0 and its so obviously different
      # # or when there are many such lines and will have to filter some less signif ones
      # if (group_col %in% c("Community","Province") & nrow(TukeyTab[TukeyTab[,"p adj"]==0, ])>7) {
      #   TukeyTab <- TukeyTab[ TukeyTab[ , "p adj"] > 0, ]
      #   removed0s <- " (removed comparisons where p=0 to save space)"
      #   # also filter out some minimally significant lines if still many
      #   if (nrow(TukeyTab) > 10) {
      #     TukeyTab <- TukeyTab[ TukeyTab[ , "p adj"] < 0.005, ]
      #   }
      #   
      # } else if (nrow(TukeyTab) > 12) {
      if (nrow(TukeyTab) > 12) {
        # if many lines, but not because so many 0s, simply filter out some that are less significant
        TukeyTab <- TukeyTab[ TukeyTab[ , "p adj"] < 0.005, ]
      }
      
      # ******************** #
      # only works if there are sig diffs in Tukey
      if (nrow(TukeyTab) > 0) {
        group1s <- sapply(rownames(TukeyTab), function(x) strsplit(x, "-")[[1]][1])
        group2s <- sapply(rownames(TukeyTab), function(x) strsplit(x, "-")[[1]][2])
        cont.max <- max(adf[,cont_col], na.rm = T)
        # arbitrarily make the space between lines based on the range of values
        rango <- (max(adf[,cont_col], na.rm = T)-min(adf[,cont_col], na.rm = T)) / 10
        # cont.max will determine the y.position for signif lines
        cont.max <- cont.max + rango
        cont.max <- seq( cont.max, cont.max+((nrow(TukeyTab)-1)*rango), by=rango)
        stat.test <- data.frame(".y."=rep("cont", nrow(TukeyTab)), 
                                "group1"=group1s, "group2"=group2s,
                                "p"=rep(NA, nrow(TukeyTab)), 
                                "p.adj"=formatC(TukeyTab[,"p adj"], format="e", digits=3),
                                "p.format"=formatC(TukeyTab[,"p adj"], format="e", digits=3),
                                "p.signif"=symnum(TukeyTab[,"p adj"], 
                                                  cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                                  symbols = c("****", "***", "**", "*", "ns")),
                                "method"=rep("Tukey", nrow(TukeyTab)),
                                "y.position"=cont.max)
        # print(stat.test)
        # print(str(stat.test))
        if (print_tukey==T)
          print(TukeyTab)
      } else {
        plot_tukey <- F
      }
      
    }
  }
  # *********************************************************** #
  
  # ******************** #
  if (length(cont_col) > 1) {
    cont_col.ti <- ifelse(sum(cont_col %in% rownames(o_rel))==length(cont_col), 
                          sprintf("Indicated %s", 
                                  gsub("Phylum","Phyla", 
                                       gsub("Genus","Genera", tlev))), "Indicated variables")
    yl <- ifelse(sum(cont_col %in% rownames(o_rel))==length(cont_col), "Abundance", "")
  } else if (cont_col %in% CONTS) {
    cont_col.ti <- cont_col
  } else {
    # if cont_col is a taxon, include the tax level in file name when saving
    cont_col.ti <- sprintf("%s: %s", tlev, cont_col)
  }
  
  if ( ! is.null(choose_yl)) yl <- choose_yl
  if ( ! is.null(choose_xl)) xl <- choose_xl
  # ******************** #
  
  
  
  
  # ******************** #
  if ( ! is.null(anova.pTab) & length(cont_col)==1 ) {
    
    # have to reverse the change of the characters that needed to be changed during lm tests so can check table
    if (cont_col == "Absconditabacteriales_(SR1)")
      cont_col.pTab <- "Absconditabacteriales_.SR1."
    else if ( ! cont_col %in% CONTS)
      # only change this for taxa variables
      cont_col.pTab <- gsub("-", "\\.", cont_col)
    else
      cont_col.pTab <- cont_col
    
    
    # pval of trait of interest
    if (cont_col.pTab %in% rownames(anova.pTab) & group_col %in% colnames(anova.pTab))
      ano.P <- anova.pTab[ cont_col.pTab, group_col ]
    else
      ano.P <- anova.pTab[ group_col, cont_col.pTab ]
    
    ano.P.title <- sprintf("\nAnova pval = %s%s", signif(as.numeric(ano.P), digits = 3), removed0s )
    
    # pvals of covariates, if any signif (ignore seqGroup, as this was included just to account for its variation)
    if (cont_col.pTab %in% rownames(anova.pTab) & group_col %in% colnames(anova.pTab)) {
      covsToCheck <- colnames(anova.pTab)[ ! colnames(anova.pTab) %in% c(group_col, "seqGroup") ]
      ano.covs.P <- anova.pTab[ cont_col.pTab, covsToCheck ]
    } else {
      covsToCheck <- colnames(anova.pTab)[ ! colnames(anova.pTab) %in% c(cont_col.pTab, "seqGroup") ]
      ano.covs.P <- anova.pTab[ group_col, covsToCheck ]
    }
    # in case only 1 cov, will lose the name, must ensure it keeps it
    names(ano.covs.P) <- covsToCheck
    
    # take only those covariates which were significant
    acp <- ano.covs.P[ ano.covs.P != "" ]
    # make string for title
    ano.covs.P.title <- ifelse(length(acp) == 0,
                               "\nOther effects:   none",
                               sprintf("\nOther effects:   %s", 
                                       paste(paste(names(acp), signif(as.numeric(acp), digits=3), 
                                                   sep="="), collapse=" || ")))
    
  } else {
    ano.P.title <- ""
    ano.covs.P.title <- ""
  }
  # ******************** #
  
  
  
  # fix levels of Age_groups
  if (group_col == "Age_groups")
    kw.box$group <- factor(kw.box$group, levels = c("Child","Teen","Adult","Senior"))
  if (startsWith(group_col, "healthy.Diversity_group"))
    kw.box$group <- factor(kw.box$group, levels = c("Low","Average","High"))
  if (group_col == "pH_label.wide")
    kw.box$group <- factor(kw.box$group, levels = c("Highly_Acidic","Acidic","Neutral","Alkaline","Highly_Alkaline"))
  if (group_col == "pH_label")
    kw.box$group <- factor(kw.box$group, levels = c("Acidic","Neutral","Alkaline"))
  # if (startsWith(group_col, "Age.Senior_other"))
  #   kw.box$group <- factor(kw.box$group, levels = c("13-60",">60"))
  
  # ******************************************************************* #
  # ******************************************************************* #
  # ******************************************************************* #
  if (length(cont_col) > 1) {
    
    # make cont column numeric
    kw.box$value <- as.numeric(as.character(kw.box$value))
    
    # ******************** #
    if ( ! is.null(anova.pTab) ) {
      ano.P <- as.numeric(anova.pTab[ cont_col, group_col ])
      names(ano.P) <- cont_col
      
      facet_labels <- sapply(cont_col, function(x) {
        head.anoP <- sprintf("%s\npval = %s", x, signif(as.numeric(ano.P[x]), digits = 3))
        
        # pvals of covariates, if any signif (ignore seqGroup, as this was included just to account for its variation)
        covsToCheck <- colnames(anova.pTab)[ ! colnames(anova.pTab) %in% c(group_col, "seqGroup") ]
        ano.covs.P <- anova.pTab[ x, covsToCheck ]
        # in case only 1 cov, will lose the name, must ensure it keeps it
        names(ano.covs.P) <- covsToCheck
        # take only those covariates which were significant
        acp <- ano.covs.P[ ano.covs.P != "" ]
        # make string for title
        ano.covs.P.title <- ifelse(length(acp) == 0,
                                   "\nOther effects:   none",
                                   sprintf("\nOther effects:   %s", 
                                           paste(paste(names(acp), signif(as.numeric(acp), digits=3), 
                                                       sep="="), collapse=" || ")))
        return(sprintf("%s%s", head.anoP, ano.covs.P.title))
      })
    } else {
      facet_labels <- sapply(cont_col, function(x) x)
    }
    
    # rename unclassified.G43 as a more understandable name
    facet_labels <- gsub("unclassified.G43","unclass.Phy", 
                         gsub("unclassified.P1","unclassified.Phylum", 
                              gsub("unclassified.G20","F.Family_XIII.UCG", 
                                   gsub("unclassified.G19","F.Clostridiales_vadinBB60",
                                        gsub("unclassified.G33","O.Saccharimonadales.UCG", facet_labels)))))
    
    # reorder facets by the order that they are entered into the function
    kw.box$variable <- factor(kw.box$variable, levels = cont_col)
    # kw.box <- kw.box[ ! is.na(kw.box$value), ]
    
    if (adjustAlphas) {
      tax.alphas <- sapply(unique(kw.box$variable), function(x) {
        yesmed <- median(kw.box[kw.box$variable == x & kw.box$group %in% c("Yes","Female",">60","1"), "value"], na.rm = T)
        nomed  <- median(kw.box[kw.box$variable == x & kw.box$group %in% c("No","Male","13-60","2"), "value"], na.rm = T)
        if (yesmed > nomed)
          return(c("yes.alpha"=1, "no.alpha"=0.65))
        else if (nomed > yesmed)
          return(c("yes.alpha"=0.65, "no.alpha"=1))
        else
          return(c("yes.alpha"=1, "no.alpha"=1))
      })
      
      colnames(tax.alphas) <- unique(kw.box$variable)
      # print(tax.alphas)
      kw.box$Alphas <- sapply(rownames(kw.box), function(x) ifelse(kw.box[x, "group"] %in% c("Yes","Female",">60","1"), 
                                                                   tax.alphas["yes.alpha", kw.box[x, "variable"]],
                                                                   tax.alphas["no.alpha", kw.box[x, "variable"]]))
      # print(head(kw.box))
      # print(tail(kw.box))
    } else {
      kw.box$Alphas <- 1
    }
    
    if (group_col == "Age_bins")
      kw.box$group <- factor(gsub("60\\+", ">60", gsub("_","-",kw.box$group)), 
                             levels = c("13-20","20-30","30-40","40-50","50-60",">60"))
    
    # ******************** #
    boxes <- ggplot(kw.box, aes(x=group, y=value, fill=group)) +
      # geom_boxplot(notch = T) + 
      geom_boxplot(notch = useNotch, outlier.shape = NA, aes(alpha=Alphas)) + 
      geom_jitter(width = 0.03, height = 0.001, alpha=0.5) +
      theme_bw() +
      theme(legend.position = ifelse(include_legend_below==T, "bottom","none"), 
            legend.title = element_blank(), legend.text = element_text(size=15),
            axis.text = element_text(size=atex.s), axis.title = element_text(size=atit.s),
            axis.text.x = element_text(angle = xAngle, hjust = 0.5), plot.title = element_text(size=pts),
            strip.text = element_text(size=sts), strip.background = element_rect(fill=str.bck)) +
      xlab(ifelse(is.null(choose_xl), 
                  gsub(".Senior_other","",group_col),
                  xl)) +
      guides(alpha=F, fill=guide_legend(reverse=T)) +
      # scale_fill_manual(values=c("#CB6767","#73BDD3"), name=group_col) +
      # xlab(NULL) +
      facet_wrap(~variable, scales=facetScales, labeller = as_labeller(facet_labels), 
                 nrow = ifelse( ! is.null(facetRows), facetRows,
                                ifelse(length(cont_col)==2, 1, 
                                       ifelse(length(cont_col)>21, 4, 
                                              ifelse(length(cont_col)>12, 3, 2))))) +
      ggtitle(ifelse(is.null(choose_title),
                     sprintf("%s vs %s", gsub(".Senior_other","",group_col), cont_col.ti),
                     choose_title)) +
      ylab(yl) + scale_fill_hue(name=group_col) #+ ylim(0,5000)
    
    
    if (group_col %in% c("Downs_Syndrome")) {
      # in order to use the same colors for disorders as in the age plot in the section:
      #   " # Scatterplot of healthyPop and DS age_diffs vs distances, combined "
      boxes <- boxes + 
        scale_fill_manual(values=c("#73BDD3","#CB6767"), name=group_col)
      
    } else if (group_col %in% c("Cystic_fibrosis")) {
      # in order to use the same colors for disorders as in the age plot in the section:
      #   " # Scatterplot of healthyPop and DS age_diffs vs distances, combined "
      boxes <- boxes + 
        scale_fill_manual(values=c("#73BDD3","#E4BA4E"), name=group_col)
      # dot_cols <- c("#e9c46a","#73BDD3")
      # lin_cols <- c("#B1871B","#2f6690")
      # dot_cols <- c("#f77f00","#73BDD3")
      # lin_cols <- c("#A35400","#2f6690")
      # dot_cols <- c("#fca311","#73BDD3")
      # lin_cols <- c("#A26402","#2f6690")
    }
    
    facet.fname <- ""
    
    
    # ******************************************************************* #
    # ******************************************************************* #
    
  } else {
    
    # ******************************************************************* #
    # ******************************************************************* #
    
    # make cont column numeric
    kw.box$cont <- as.numeric(as.character(kw.box$cont))
    
    # ******************** #
    if (is.null(facet) ) {#| group_col=="seqGroup") {
      # ggplot(kw.box, aes(x=reorder(group,-as.numeric(as.character(cont)),median), y=as.numeric(as.character(cont)), 
      #                    fill=reorder(group,-as.numeric(as.character(cont)),median))) +
      boxes <- ggplot(kw.box, aes(x=group, y=cont, fill=group)) +
        # geom_boxplot(notch = T) + 
        geom_boxplot(notch = useNotch, outlier.shape = NA) + 
        geom_jitter(width = 0.03, height = 0.001, alpha=0.5) +
        theme(legend.position = "none", axis.text = element_text(size=atex.s), axis.title = element_text(size=atit.s),
              axis.text.x = element_text(angle = xAngle, hjust = 0.5), plot.title = element_text(size=pts)) +
        xlab(group_col) +
        # xlab(NULL) +
        ggtitle(sprintf("%s vs %s%s%s", group_col, cont_col.ti, ano.P.title, ano.covs.P.title)) +
        ylab(gsub("Div.Shannon","Shannon index", cont_col)) + scale_fill_hue(name=group_col) #+ ylim(0,5000)
      
      facet.fname <- ""
      
      # ******************** #
    } else {
      
      if (grepl("~", facet)) {
        # if faceting by combination of 2 variables
        fv1 <- strsplit(facet, "~")[[1]][1]
        fv2 <- strsplit(facet, "~")[[1]][2]
        
        kw.box$FV1 <- gqs[ , fv1]
        kw.box$FV2 <- gqs[ , fv2]
        
        # fix levels of Age_groups
        if (fv1 == "Age_groups")
          kw.box$FV1 <- factor(kw.box$FV1, levels = c("Child","Teen","Adult","Senior"))
        if (fv2 == "Age_groups")
          kw.box$FV2 <- factor(kw.box$FV2, levels = c("Child","Teen","Adult","Senior"))
        
        if (startsWith(fv1, "healthy.Diversity_group"))
          kw.box$FV1 <- factor(kw.box$FV1, levels = c("Low","Average","High"))
        if (startsWith(fv2, "healthy.Diversity_group"))
          kw.box$FV2 <- factor(kw.box$FV2, levels = c("Low","Average","High"))
        
        # ggplot(kw.box, aes(x=reorder(group,-as.numeric(as.character(cont)),median), y=as.numeric(as.character(cont)), 
        #                    fill=reorder(group,-as.numeric(as.character(cont)),median))) +
        boxes <- ggplot(kw.box, aes(x=group, y=cont, fill=group)) +
          # geom_boxplot(notch = T) + 
          geom_boxplot(notch = useNotch, outlier.shape = NA) + 
          geom_jitter(width = 0.03, height = 0.001, alpha=0.5) +
          theme(legend.position = "none", axis.text = element_text(size=atex.s), axis.title = element_text(size=atit.s),
                axis.text.x = element_text(angle = xAngle, hjust = 0.5), plot.title = element_text(size=pts)) +
          xlab(group_col) +
          # xlab(NULL) +
          facet_wrap(FV1~FV2) +
          ggtitle(sprintf("%s vs %s\nfacets: %s%s%s", group_col, cont_col.ti, facet, ano.P.title, ano.covs.P.title)) +
          ylab(gsub("Div.Shannon","Shannon index", cont_col)) + scale_fill_hue(name=group_col) #+ ylim(0,5000)
        
      } else {
        
        # kw.box$Group <- gqs[ , facet]
        kw.box$Group <- gqs[ ! is.na(gqs[ , group_col ]), facet]
        kw.box <- kw.box[ ! is.na(kw.box$Group), ]
        # fix levels of Age_groups
        if (facet == "Age_groups")
          kw.box$Group <- factor(kw.box$Group, levels = c("Child","Teen","Adult","Senior"))
        if (startsWith(facet, "healthy.Diversity_group"))
          kw.box$Group <- factor(kw.box$Group, levels = c("Low","Average","High"))
        
        # ggplot(kw.box, aes(x=reorder(group,-as.numeric(as.character(cont)),median), y=as.numeric(as.character(cont)), 
        #                    fill=reorder(group,-as.numeric(as.character(cont)),median))) +
        boxes <- ggplot(kw.box, aes(x=group, y=cont, fill=group)) +
          # geom_boxplot(notch = T) + 
          geom_boxplot(notch = useNotch, outlier.shape = NA) + 
          geom_jitter(width = 0.03, height = 0.001, alpha=0.5) +
          theme(legend.position = "none", axis.text = element_text(size=atex.s), axis.title = element_text(size=atit.s),
                axis.text.x = element_text(angle = xAngle, hjust = 0.5), plot.title = element_text(size=pts)) +
          xlab(group_col) +
          # xlab(NULL) +
          facet_wrap(~Group) +
          ggtitle(sprintf("%s vs %s\nfacets: %s%s%s", group_col, cont_col.ti, facet, ano.P.title, ano.covs.P.title)) +
          ylab(gsub("Div.Shannon","Shannon index", cont_col)) + scale_fill_hue(name=group_col) #+ ylim(0,5000)
      }
      
      
      facet.fname <- sprintf("-facet.%s", facet)
    }
    
  }
  
  # ******************************************************************* #
  # ******************************************************************* #
  # ******************************************************************* #
  
  
  
  # ******************** #
  if (plot_tukey==T & is.null(facet)) {
    if (length(table(kw.box$group)) >= 3)
      # no need to print pval if only 2 groups, since it is same as pval in title
      boxes <- boxes + stat_pvalue_manual(stat.test, label = "p = {p.adj}", inherit.aes=F)
  }
  
  # ******************** #
  if (flipCoords == T) {
    boxes <- boxes + coord_flip()
  }
  
  # ******************** #
  if (save.boxes == T) {
    
    if (cont_col %in% CONTS)
      cont_col.fname <- cont_col
    else
      # if cont_col is a taxon, include the tax level in file name when saving
      cont_col.fname <- sprintf("%s.%s", tlev, cont_col)
    
    # make directory if necessary
    dir.create(sprintf("%s/figures/%s", p2_dir, plot_dirs), showWarnings = F)
    
    suppressMessages(print( 
      boxes +
        ggsave(sprintf("%s/figures/%s/%s-%s%s.png", p2_dir, plot_dirs, group_col, cont_col.fname, facet.fname),
               width = 9.23, height = 7.04, device = "png"))) 
    # suppress the message that occurs because of notches going outsides hinges
    # plot = suppressMessages(print(boxes)) )
    
  } else {
    # suppress the message that occurs because of notches going outsides hinges
    suppressMessages(print(boxes))
  }
  # ******************** #
}
# ******************** #











# ****************************************************************************************************************** #
# *********************************************************** #
# Chi-squared test ####
# *********************************************************** #

library(graphics)
library(vcd)
library(corrplot)


# ********************************************************************************* #
fill_chi_table <- function(comparison, tl, interest, group_qs, keep_full=FALSE, 
                           gq_for_disorders=NULL, externalGroupQ=NULL, externalGroupName=NULL) {
  
  # first get appropriate data
  if (comparison == "TvsQ") {
    if (! is.null(externalGroupQ)) {
      data.group <- cbind(t( apply( gloms[[ tl ]][ , interest ], 2, function(x) ifelse(x==0, 0, 1) ) ),
                          as.matrix(externalGroupQ[ interest ]))
      colnames(data.group) <- c(rownames(gloms[[ tl ]]), externalGroupName)
      
      chi.cols <- externalGroupName
      chi.rows <- rownames(gloms[[ tl ]])
      
    } else {
      data.group <- cbind(t( apply( gloms[[ tl ]][ , interest ], 2, function(x) ifelse(x==0, 0, 1) ) ),
                          as.matrix( SLL2@sam_data[ interest, group_qs] ))
      chi.cols <- rownames(gloms[[ tl ]])
      chi.rows <- group_qs
    }
    
    
  } else if (comparison == "taxa") {
    data.group <- t( apply( gloms[[ tl ]][ , interest ], 2, function(x) ifelse(x==0, 0, 1) ) )
    chi.rows <- chi.cols <- rownames(gloms[[ tl ]])
    
  } else if (comparison == "questions") {
    
    if (! is.null(gq_for_disorders)) {
      # This is for the case where want to check only 1 group_q (eg one of the disorders) against all group_qs
      data.group <- as.matrix( SLL2@sam_data[ interest, gq_for_disorders] )
      chi.cols <- gq_for_disorders # will be the full list of group_qs
      chi.rows <- group_qs # will be the individual group_q to check (eg disorder)
      
    } else if (! is.null(externalGroupQ)) {
      # This is for the case where want to check only 1 group_q (eg clusters within one of the disorders) 
      #   against all group_qs
      data.group <- cbind(as.matrix( SLL2@sam_data[ interest, group_qs] ), as.matrix(externalGroupQ[ interest ]))
      colnames(data.group) <- c(group_qs, externalGroupName)
      
      chi.cols <- externalGroupName # will be the individual group_q to check (eg clusters within a disorder)
      chi.rows <- group_qs # will be the full list of group_qs
      
    } else {
      data.group <- as.matrix( SLL2@sam_data[ interest, group_qs] )
      chi.rows <- chi.cols <- group_qs
    }
    
  }
  
  # create table
  chi.matrix <- matrix(NA, nrow = length(chi.rows), ncol = length(chi.cols))
  rownames(chi.matrix) <- chi.rows
  colnames(chi.matrix) <- chi.cols
  
  ### fill table
  for (i in chi.rows) {
    for (j in chi.cols) {
      
      conting_tab <- table(data.group[ , i], data.group[ , j], dnn = c(i, j))
      conting_tab <- conting_tab[rownames(conting_tab) != "No Sabe/No Contesta", 
                                 colnames(conting_tab) != "No Sabe/No Contesta"]
      
      if (class(conting_tab)=="integer") {
        chi.matrix[i,j] <- 1
      } else if (nrow(conting_tab) < 2 | ncol(conting_tab) < 2) {
        chi.matrix[i,j] <- 1
      } else if (sum(conting_tab)==0) {
        chi.matrix[i,j] <- 1
      } else {
        chi <- chisq.test(conting_tab)
        
        # since shouldnt use Chi-square when any expected value is below 5, use Fishers test here instead
        # if (length( chi$expected[chi$expected<5]) > 0) {
        #   print(c(i,j, length( chi$expected[chi$expected<5])))
        #   chi <- fisher.test(conting_tab, simulate.p.value = T)
        # }
        
        chi.matrix[i,j] <- chi$p.value
      }
      
    }
  }
  
  if (keep_full) {
    return(chi.matrix)
    
  } else {
    
    # for those values that are given as NaN because at least one of the categories had all 0s 
    chi.matrix[ is.nan(chi.matrix) ] <- 1
    # adjust for multiple testing
    chi.matrix.adj <- apply(chi.matrix, 2, p.adjust, method='fdr')
    # make p-values for diagonals = 1 when appropriate because we dont care about self-correlations,
    # and cause this will ruin the color scale of the heatmaps
    if (comparison != "TvsQ") diag( chi.matrix.adj ) <- 1
    
    #at least 1 good p value within samples
    mins <- apply(chi.matrix.adj,2,min)
    goodPs <- mins[mins < 0.05] #-log10(0.05) == 1.3013
    #at least 1 good p within taxa
    tmins <- apply(chi.matrix.adj,1,min)
    tgoodPs <- tmins[tmins < 0.05]
    
    only.good.chi.matrix.adj <- chi.matrix.adj[names(tgoodPs),names(goodPs)]
    only.good.chi.matrix.adj[only.good.chi.matrix.adj >= 0.05] <- ''
    
    if (class(only.good.chi.matrix.adj)=="character" & length(only.good.chi.matrix.adj)==1) {
      only.good.chi.matrix.adj <- as.data.frame(only.good.chi.matrix.adj)
      colnames(only.good.chi.matrix.adj) <- names(goodPs)
      rownames(only.good.chi.matrix.adj) <- names(tgoodPs)
    }
    
    return(as.matrix(only.good.chi.matrix.adj))
  }
  
}
# ********************************************************************************* #  

# ********************************************************************************* #  
get_contingency_table <- function(ro, co, comparison, tl, interest, group_qs, phy, glomTab) {
  # first get appropriate data
  if (comparison == "TvsQ") {
    if (ro %in% group_qs) {
      t1 <- as.matrix( phy@sam_data )[ interest, ro ]
      t2 <- apply( glomTab[[ tl ]][ , interest ], 2, function(x) ifelse(x==0, "Absent", "Present") )[ co, ]
    } else {
      t1 <- apply( glomTab[[ tl ]][ , interest ], 2, function(x) ifelse(x==0, "Absent", "Present") )[ ro, ]
      t2 <- as.matrix( phy@sam_data )[ interest, co ]
    }
    
  } else if (comparison == "taxa") {
    t1 <- apply( glomTab[[ tl ]][ , interest ], 2, function(x) ifelse(x==0, "Absent", "Present") )[ ro, ]
    t2 <- apply( glomTab[[ tl ]][ , interest ], 2, function(x) ifelse(x==0, "Absent", "Present") )[ co, ]
    
  } else if (comparison == "questions") {
    t1 <- as.matrix( phy@sam_data )[ interest, ro ]
    t2 <- as.matrix( phy@sam_data )[ interest, co ]
    
  }
  
  contingency <- table(t1, t2, dnn = c(ro, co))
  
  # order diversity groups logically
  if (startsWith(ro, "Diversity_group")) {contingency <- contingency[ c("Low","Average","High"), ]
  } else if (startsWith(co, "Diversity_group")) {contingency <- contingency[ , c("Low","Average","High") ]}
  
  contingency <- contingency[rownames(contingency) != "No Sabe/No Contesta", 
                             colnames(contingency) != "No Sabe/No Contesta"]
  
  if (nrow(contingency) > ncol(contingency)) contingency <- t(contingency)
  
  return(contingency)
}
# ********************************************************************************* #  


























# ****************************************************************************************************************** ####
# ****************************************************************************************************************** ####


# ****************************************************************************************************************** #

# ******************************************************* #
# Pie charts / Donut charts of abundances ####
# ******************************************************* #

donuts <- function(x, inner, outer, group = 1, labels_out = NA, labels_in = NA, col = NULL, ti = NULL, radius = c(.7, 1, 1.05)) {
  group <- rep_len(group, length(x))
  ug  <- unique(group)
  tbl <- table(group)[order(ug)]
  
  col <- if (is.null(col))
    seq_along(ug) else rep_len(col, length(ug))
  col.main <- Map(rep, col[seq_along(tbl)], tbl)
  col.sub  <- lapply(col.main, function(x) {
    al <- head(seq(0, 1, length.out = length(x) + 2L)[-1L], -1L)
    # rev(Vectorize(adjustcolor)(x, alpha.f = al))
    Vectorize(adjustcolor)(x, alpha.f = al)
  })
  
  if (outer=="Genus") { ia <- -30; out.cex <- 1.35
  } else { ia <- 20; out.cex <- 1.35
  }
  
  plot.new()
  
  par(new = TRUE)
  pie(x, border = NA, radius = radius[2L],# init.angle = ia,
      col = unlist(col.sub), labels = labels_out, cex=out.cex,
      main = ti, cex.main = 2)
  
  par(new = TRUE)
  pie(x, border = NA, radius = radius[1L],# init.angle = ia,
      col = unlist(col.main), labels = NA)#labels_in, cex=0.6)
  
  # par(new = TRUE)
  # p <- floating.pie(0.5,0.5, x, radius = radius[2L]/3.5,
  #                   col = unlist(col.sub), main = ti, border = NA)
  # pie.labels(0.5,0.5, p, labels = labels_out, radius = radius[3L]/3.5, minangle = 0.08, cex = out.cex)
  # title(main = ti, cex.main = 2)
  # 
  # par(new = TRUE)
  # floating.pie(0.5,0.5, x, radius = radius[1L]/3.5, col = unlist(col.main), border = NA)
  
  legend("topright", legend = labels_in, col = names(col.main), pch = 15, bty = 'n', cex = 1.9, title = inner)
  # text(x = c(0.1, -.25, 0, .35, .5, 0.55), y = c(0.15, 0, -.3, -.25, -.1, -0.025), 
  #      labels = unique(phyla_for_pie$phylum), col = 'white', cex = 1.2)
}
# ******************************************************** #

# ******************************************************** #
get_donut_structures <- function(inner, outer) {
  
  if (inner=="Phylum" & outer=="Species") {
    in.num <- 3
    out.nums <- c(8, 4, 2)
    inner_struc <- c(1L,1L,1L,1L,1L,1L,1L,1L,1L,2L,2L,2L,2L,2L,3L,3L,3L,4L)
    outer_struc <- c(1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,11L,12L,13L,14L,15L,16L,17L,18L)
    rownames_val <- -18L
  } else if (inner=="Phylum" & outer=="Genus") {
    in.num <- 4
    out.nums <- c(5, 3, 2, 2)
    inner_struc <- c(1L,1L,1L,1L,1L,1L,2L,2L,2L,2L,3L,3L,3L,4L,4L,4L,5L)
    outer_struc <- c(1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,11L,12L,13L,14L,15L,16L,17L)
    rownames_val <- -17L
  } else if (inner=="Class" & outer=="Species") {
    in.num <- 4
    out.nums <- c(6, 4, 1, 1)
    inner_struc <- c(1L,1L,1L,1L,1L,1L,1L,2L,2L,2L,2L,2L,3L,3L,4L,4L,5L)
    outer_struc <- c(1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,11L,12L,13L,14L,15L,16L,17L)
    rownames_val <- -17L
  } else if (inner=="Class" & outer=="Genus") {
    in.num <- 4
    out.nums <- c(4, 3, 2, 1)
    inner_struc <- c(1L,1L,1L,1L,1L,2L,2L,2L,2L,3L,3L,3L,4L,4L,5L)
    outer_struc <- c(1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,11L,12L,13L,14L,15L)
    rownames_val <- -15L
  }
  
  return( list(in.num, out.nums, inner_struc, outer_struc, rownames_val) )
}

# ******************************************************** #

# ******************************************************** #
get_donut_objects <- function(inner, outer, samps, min_abund, glomTab) {
  
  samps <- sapply(samps, function(x) gsub('\\.', '-', x))
  
  # otus.in <- glomTab[[ inner ]][ , samps ]
  otus.out <- glomTab[[ outer ]][ , samps ]
  taxTs <- taxTables.both
  
  donut_strucs <- get_donut_structures(inner, outer)
  in.num <- donut_strucs[[ 1 ]]
  out.nums <- donut_strucs[[ 2 ]]
  inner_struc <- donut_strucs[[ 3 ]]
  outer_struc <- donut_strucs[[ 4 ]]
  rownames_val <- donut_strucs[[ 5 ]]
  
  # get vector of most abundant taxa at the level "inner" (for inner ring)
  # if (length(samps) == 1) {
  #   topInner <- c(names(sort(otus.in, decreasing = T))[1:in.num], "Other")
  # } else {
  #   topInner <- c(names(sort(rowSums(otus.in), decreasing = T))[1:in.num], "Other")
  # }
  # will use the top 5 "inner" from the full dataset, so all plots can be compared effectively
  topInner <- c(names(sort(rowSums(glomTab[[ inner ]]), decreasing = T))[1:in.num], "Other")
  
  topOuter <- list()
  for (i in topInner) {
    # take top 3 genera per phylum for the top 4 phyla, then only the top 2 for the 5th and "Other" 
    # since they are small and get cluttered on the pie chart
    if (i=="Other") {
      outer_in_inner <- rownames(taxTs[[ outer ]][ ! taxTs[[ outer ]][,inner] %in% topInner[1:in.num], ]);
      # outer_in_inner <- rownames(phy@tax_table[ ! phy@tax_table[,inner] %in% topInner[1:5],]);
      # topGenera[[i]] <- c(names(sort(rowSums(otu_table(SLL1)[gen_in_phy, which_samples]), decreasing = T))[1:2], "Other", NA)
      topOuter[[i]] <- c("Other", rep(NA, max(out.nums)))#c("Other", NA, NA, NA)
    } else {
      if ( class(taxTs[[ outer ]][taxTs[[ outer ]][,inner]==i, ])=="character" ) {
        # in the case that only one "outer" taxa within the "inner" taxa
        # outer_in_inner <- taxTs[[ outer ]][taxTs[[ outer ]][,inner]==i, outer]
        outer_in_inner <- rownames(taxTs[[ outer ]][taxTs[[ outer ]][,inner]==i, ])
        topOuter[[i]] <- c(outer_in_inner, rep(NA, max(out.nums)))
      } else {
        outer_in_inner <- rownames(taxTs[[ outer ]][taxTs[[ outer ]][,inner]==i, ])
        
        if (length(samps)==1) {
          # if only for 1 sample, otus.out will be a vector, not a table, rowSums not applicable
          topOuter[[i]] <- c(names(sort(otus.out[ outer_in_inner ], decreasing = T))[1:out.nums[match(i,topInner)]], "Other")
        } else {
          topOuter[[i]] <- c(names(sort(rowSums(otus.out[outer_in_inner, ]), decreasing = T))[1:out.nums[match(i,topInner)]], "Other")
        }
        
        topOuter[[i]] <- c(topOuter[[i]], rep(NA, max(out.nums)-out.nums[match(i,topInner)]))
      }
      
      # outer_in_inner <- rownames(phy@tax_table[phy@tax_table[,inner]==i,]);
      # ifelse(i %in% topInner[2:5],
      #        topOuter[[i]] <- c(names(sort(rowSums(otus.out[outer_in_inner, ]), decreasing = T))[1:2], "Other", NA), 
      #        topOuter[[i]] <- c(names(sort(rowSums(otus.out[outer_in_inner, ]), decreasing = T))[1:3], "Other"))
      
    }
  }
  
  
  # get list of abundances
  topOut.frame <- as.data.frame(topOuter)
  outerAbunds <- list()
  for (i in topInner) {
    outs <- list()
    add_to_other <- 0
    for (ou in topOut.frame[,i]) {
      if (is.na(ou)) {
        # do nothing 
      } else if (ou == "Other") {
        if (i == "Other") {
          # if genus and phylum are both "Other"
          outer_in_inner <- rownames(taxTs[[ outer ]][ ! taxTs[[ outer ]][,inner] %in% topInner[1:in.num], ])
          # outer_in_inner <- rownames(phy@tax_table[ ! phy@tax_table[,inner] %in% topInner[1:5],])
          o_in_i_other <- outer_in_inner[ ! outer_in_inner %in% topOuter[["Other"]] ]
          
          if (length(samps)==1) {
            # if only for 1 sample, cant take mean values
            outs[[ou]] <- sum(otus.out[ o_in_i_other ])
          } else {
            outs[[ou]] <- sum(rowMeans(otus.out[o_in_i_other, ]))
          }
          
        } else {
          # if genus is "Other", but not Phylum
          outer_in_inner <- rownames(taxTs[[ outer ]][ taxTs[[ outer ]][,inner] == i, ])
          # outer_in_inner <- rownames(phy@tax_table[ phy@tax_table[,inner] == i])
          o_in_i_other <- outer_in_inner[ ! outer_in_inner %in% topOuter[[i]] ]
          
          # if only 1 taxa in "Other", will only be one vector, not a table, cant use rowMeans
          if ( length(o_in_i_other) == 1 ) {
            
            if (length(samps)==1) { outs[[ou]] <- otus.out[ o_in_i_other ] # if only for 1 sample, cant take mean values
            } else { outs[[ou]] <- mean(otus.out[ o_in_i_other, ]) }
            
          } else if ( length(o_in_i_other) == 0) { 
            # if 0 taxa in "Other", just give value of 0
            outs[[ou]] <- 0
          } else {
            
            if (length(samps)==1) { outs[[ou]] <- sum(otus.out[ o_in_i_other ]) # if only for 1 sample, cant take mean values
            } else { outs[[ou]] <- sum(rowMeans(otus.out[o_in_i_other, ])) }
            
          }
          # in case any non-"Other" were less than min_abund, they were removed, must add values to "Other"
          outs[[ou]] <- outs[[ou]] + add_to_other
        }
      } else {
        # if neither genus nor phylum are "Other"
        if (length(samps)==1) {
          # if only for 1 sample, cant take mean values
          outs[[ou]] <- otus.out[ ou ]
        } else {
          outs[[ou]] <- mean(otus.out[ou, ])
        }
        
        # if the value here is less than min_abund, add its value to "Other", then remove it
        if (outs[[ou]] < min_abund) {
          # add_to_other <- add_to_other + outs[[ou]]
          outs[[ou]] <- NA
          topOuter[[ i ]][ topOuter[[ i ]] == ou ] <- NA
          # adjust structure of donut chart
          inner_struc <- inner_struc[ -match(match(i, names(topOuter)), inner_struc) ]
          outer_struc <- outer_struc[ -length(outer_struc) ]
          rownames_val <- rownames_val + 1
        }
      }
    }
    # then add the vector of genus abundances to each list element (a phylum name)
    outerAbunds[[i]] <- unname(unlist(outs))
    # in case some values were removed because less than min_abund, must adjust the order to put NAs at end
    outerAbunds[[i]] <- c( sort(outerAbunds[[i]][1:length(outerAbunds[[i]])-1], decreasing = T), 
                           outerAbunds[[i]][length(outerAbunds[[i]])] )
    # then adjust the orders of the names in topOuter in the same manner
    topOuter[[ i ]] <- c( topOuter[[ i ]][ ! is.na(topOuter[[ i ]]) ], topOuter[[ i ]][ is.na(topOuter[[ i ]]) ] )
  }
  
  pie_title <- sprintf("%s entre %s", outer, inner)
  # adjust "Other" names where appropriate - for clarification
  topInner[in.num+1] <- sprintf("Other-%s",inner)
  names(outerAbunds)[in.num+1] <- sprintf("Other-%s",inner)
  names(topOuter)[in.num+1] <- sprintf("Other-%s",inner)
  topOuter[[names(topOuter)[in.num+1]]][1] <- sprintf("Other-%s",inner)
  
  return( list(topInner, topOuter, outerAbunds, inner_struc, outer_struc, rownames_val, pie_title) )
}
# ******************************************************** #

# ******************************************************** #
plot_donut <- function(inner.tax, outer.tax, samps, min_abund, glomTab, to.save = F, group = NULL, fname = NULL) {
  
  donut_objects <- get_donut_objects(inner.tax, outer.tax, samps, min_abund, glomTab)
  topInner <- donut_objects[[ 1 ]]
  topOuter <- donut_objects[[ 2 ]]
  outerAbunds <- donut_objects[[ 3 ]]
  inner_struc <- donut_objects[[ 4 ]]
  outer_struc <- donut_objects[[ 5 ]]
  rownames_val <- donut_objects[[ 6 ]]
  pie_title <- donut_objects[[ 7 ]]
  
  fname <- gsub('\\.', '-', fname)
  if ( ! is.null(group) && group != "individuals") {
    pie_title <- sprintf("%s\n%s", pie_title, fname)
  }
  
  # make data.frame for plot
  taxa_for_pie <- structure(list(inner = structure(inner_struc,
                                                   .Label = topInner,#rep(topPhyla, each=4),
                                                   class = "factor"),
                                 outer = structure(outer_struc,
                                                   .Label = unname(unlist(topOuter))[ ! is.na(unname(unlist(topOuter))) ],
                                                   class = "factor"),
                                 abunds = unname(unlist(outerAbunds))),
                            .Names = c("inner", "outer", "abunds"),
                            row.names = c(NA, rownames_val), 
                            class = "data.frame")
  
  taxa_for_pie$total <- with(taxa_for_pie, ave(abunds, inner, FUN = sum))
  
  # plot donut
  if (to.save == T) {
    # first open the png device to prepare it to save an image at the given location
    png(filename = sprintf("%s/figures/Abundances/donut_charts/%s/%s.donut_chart.png", p2_dir, group, gsub(' ', '_', fname)), 
        # width = 1505, height = 951, pointsize = 16)
        width = 1600, height = 975, pointsize = 16)
    # then plot the image, which will get saved at the above location
    with(taxa_for_pie,
         donuts(abunds, inner.tax, outer.tax, inner, sprintf('%s: %s%%', outer, round(abunds,2)), unique(inner),
                col = c('cyan2','red','orange','green','dodgerblue2','grey'),
                ti = pie_title)
    )
    # then turn off png device so it can be repeated for next image
    dev.off()
  } else {
    # simply plot here to view figure without opening a separate device to save
    # print(taxa_for_pie$inner)
    with(taxa_for_pie,
         donuts(abunds, inner.tax, outer.tax, inner, sprintf('%s: %s%%', outer, round(abunds,2)), unique(inner),
                col = c('cyan2','red','orange','green','dodgerblue2','grey'),
                ti = pie_title)
    )
  }
  
}
# ******************************************************** #










# ****************************************************************************************************************** #
# STACKED BAR charts for comparison between groups ####
# ******************************************************************************************************************* #
library(RColorBrewer)

plot_stacked_bars <- function(samps.to.plot=NULL, lab=NULL, plotType="singleGroup", 
                              stp2=NULL, lab2=NULL, stp3=NULL, lab3=NULL, orientation=NULL, ang=0,
                              saveStacked=F, plotDir=NULL) {
  
  dis_cols <- c("Celíaco","Fibrosis Quística","Gingivitis-Periodontitis","Síndrome de Down",
                "Diabetes","Hipertensión","Migrañas","El resto de población analizada")
  
  # for most plots, will not show legend, but will just for the individual plots
  useGuide <- FALSE
  # and will have x tick labels for all plots except individuals, which will be blank
  if (toupper(orientation) %in% c("VERT","VERTICAL","V")) atx <- element_text(size=15, angle=ang, hjust = 1)
  else atx <- element_text(size=15)
  
  # *************************** #
  if (plotType=="singleGroup") {
    
    title.extra <- sprintf('\n%s', lab)
    aty <- element_blank()
    lel <- element_line()
    
    # if only an individual sample
    if (length(samps.to.plot)==1) {
      lab.file <- sprintf("individuals/%s",lab)
      useGuide <- guide_legend(reverse=T)
      aty <- element_text(size=14)
      atx <- element_blank()
      
      SLL2.top15.toplot <- as.matrix(SLL2.top15)[ , samps.to.plot]
      SLL2.top15.toplot <- c(100-sum(SLL2.top15.toplot), SLL2.top15.toplot)
      
    } else {
      if (lab %in% names(samps.regions)) {
        lab.file <- sprintf("regions/%s", gsub(' ','_',iconv(lab,to='ASCII//TRANSLIT')))
      } else if (endsWith(lab, "-healthy")) {
        lab.file <- sprintf("regions_healthy/%s", gsub(' ','_',iconv(lab,to='ASCII//TRANSLIT')))
        lel <- element_line(size = 4)
      } else if (lab %in% dis_cols) {
        lab.file <- sprintf("disorders/%s", gsub(' ','_',iconv(lab,to='ASCII//TRANSLIT')))
      } else {
        lab.file <- gsub(' ','_',iconv(lab,to='ASCII//TRANSLIT'))
      }
      
      SLL2.top15.toplot <- rowMeans(SLL2.top15[ , samps.to.plot])
      SLL2.top15.toplot <- c(100-sum(SLL2.top15.toplot), SLL2.top15.toplot)
    }
    
    names(SLL2.top15.toplot)[1] <- "Otros"
    tax.mps.m <- reshape2::melt(as.matrix(SLL2.top15.toplot))
    colnames(tax.mps.m) <- c("Genus","Samples","value")
    tax.mps.m$Samples <- lab
    
    # } else if (plotType=="idr") {
    #   
    #   dis_cols <- c("Celiac","Cystic_fibrosis","Gingivitis_periodontitis","Downs_Syndrome")
    #   disorders <- dis_cols[SLL2.meta[samps.to.plot, dis_cols]==1]
    #   
    #   aty <- element_text(size=17)
    
    # *************************** #
  } else if (plotType=="All_disorders") {
    
    title.extra <- ''
    aty <- element_text(size=17)
    lel <- element_line()
    lab.file <- "disorders/All_disorders"
    ang <- 90
    
    SLL2.top15.celiac <- SLL2.top15[ , samps.celiac ]
    SLL2.top15.cf     <- SLL2.top15[ , samps.cf ]
    SLL2.top15.gp     <- SLL2.top15[ , samps.gp ]
    SLL2.top15.down   <- SLL2.top15[ , samps.down ]
    SLL2.top15.diab   <- SLL2.top15[ , samps.diabetes ]
    SLL2.top15.hiper  <- SLL2.top15[ , samps.hiperten ]
    SLL2.top15.migr   <- SLL2.top15[ , samps.migr ]
    SLL2.top15.noDis  <- SLL2.top15[ , samps.noDisorder ]
    
    tax.mps <- cbind(rowMeans(SLL2.top15.noDis), rowMeans(SLL2.top15.migr), rowMeans(SLL2.top15.hiper), 
                     rowMeans(SLL2.top15.diab), rowMeans(SLL2.top15.down), rowMeans(SLL2.top15.gp), 
                     rowMeans(SLL2.top15.cf), rowMeans(SLL2.top15.celiac))
    colnames(tax.mps) <- rev(dis_cols)
    # tax.sd <- cbind(rowSds(as.matrix(SLL2.top15.region)), rowSds(as.matrix(SLL2.top15.disorder)), rep(0,15))
    tax.mps <- rbind(apply(tax.mps,2, function(x) 100-sum(x)), tax.mps)
    rownames(tax.mps)[1] <- "Otros"
    
    # to order bars appropriately
    if (toupper(orientation) %in% c("VERT","VERTICAL","V")) tax.mps <- tax.mps[ , rev(colnames(tax.mps))]
    
    tax.mps.m <- reshape2::melt(tax.mps)
    colnames(tax.mps.m) <- c("Genus","Samples","value")
    
    # *************************** #
  } else if (plotType=="All_regions") {
    
    title.extra <- ''
    aty <- element_text(size=17)
    lel <- element_line()
    lab.file <- "regions/All_regions"
    ang <- 90
    
    SLL2.top15.allRegions <- sapply(rev(sort(names(samps.regions))), function(x) SLL2.top15[ , samps.regions[[ x ]] ])
    
    tax.mps <- cbind(sapply(names(SLL2.top15.allRegions), function(x) rowMeans(SLL2.top15.allRegions[[ x ]]) ))
    colnames(tax.mps) <- sapply(colnames(tax.mps), function(x) sprintf("%s\n(n=%s)", x, length(samps.regions[[x]])))
    
    # tax.sd <- cbind(rowSds(as.matrix(SLL2.top15.region)), rowSds(as.matrix(SLL2.top15.disorder)), rep(0,15))
    tax.mps <- rbind(apply(tax.mps,2, function(x) 100-sum(x)), tax.mps)
    rownames(tax.mps)[1] <- "Otros"
    
    # to order bars appropriately
    if (toupper(orientation) %in% c("VERT","VERTICAL","V")) tax.mps <- tax.mps[ , rev(colnames(tax.mps))]
    
    tax.mps.m <- reshape2::melt(tax.mps)
    colnames(tax.mps.m) <- c("Genus","Samples","value")
    
    # *************************** #
  } else if (plotType=="All_regions_healthy") {
    
    title.extra <- ''
    aty <- element_text(size=17)
    lel <- element_line()
    lab.file <- "regions/All_regions_healthy"
    ang <- 90
    
    SLL2.top15.allRegions <- sapply(rev(sort(names(samps.regions.healthy))), function(x) SLL2.top15[ , samps.regions.healthy[[ x ]] ])
    
    tax.mps <- cbind(sapply(names(SLL2.top15.allRegions), function(x) rowMeans(SLL2.top15.allRegions[[ x ]]) ))
    colnames(tax.mps) <- sapply(colnames(tax.mps), function(x) sprintf("%s\n(n=%s)", x, length(samps.regions.healthy[[x]])))
    
    # tax.sd <- cbind(rowSds(as.matrix(SLL2.top15.region)), rowSds(as.matrix(SLL2.top15.disorder)), rep(0,15))
    tax.mps <- rbind(apply(tax.mps,2, function(x) 100-sum(x)), tax.mps)
    rownames(tax.mps)[1] <- "Otros"
    
    # to order bars appropriately
    if (toupper(orientation) %in% c("VERT","VERTICAL","V")) tax.mps <- tax.mps[ , rev(colnames(tax.mps))]
    
    tax.mps.m <- reshape2::melt(tax.mps)
    colnames(tax.mps.m) <- c("Genus","Samples","value")
    
    # *************************** #
  } else if (plotType=="Ages") {
    
    title.extra <- ''
    aty <- element_text(size=17)
    lel <- element_line()
    lab.file <- "Ages"
    
    SLL2.top15.teens  <- SLL2.top15[ , sample_names(subset_samples(SLL2, Age<19)) ]
    SLL2.top15.adult  <- SLL2.top15[ , sample_names(subset_samples(SLL2, Age>=19 & Age<56)) ]
    SLL2.top15.senior <- SLL2.top15[ , sample_names(subset_samples(SLL2, Age>=56)) ]
    
    if (! is.null(samps.to.plot)) { # to plot age groups for particular groups (celiacs, males, etc)
      SLL2.top15.teens <- as.matrix(SLL2.top15.teens[ , colnames(SLL2.top15.teens) %in% samps.to.plot ])
      SLL2.top15.adult <- as.matrix(SLL2.top15.adult[ , colnames(SLL2.top15.adult) %in% samps.to.plot ])
      SLL2.top15.senior <- as.matrix(SLL2.top15.senior[ , colnames(SLL2.top15.senior) %in% samps.to.plot ])
      
      title.extra <- sprintf('\n%s', lab)
      if (lab %in% dis_cols) {
        lab.file <- sprintf("disorders/Ages-%s", gsub(' ','_',iconv(lab,to='ASCII//TRANSLIT')))
      } else {
        lab.file <- sprintf("Ages-%s", gsub(' ','_',iconv(lab,to='ASCII//TRANSLIT')))
      }
    }
    
    tax.mps <- cbind(rowMeans(SLL2.top15.senior), rowMeans(SLL2.top15.adult), rowMeans(SLL2.top15.teens))
    colnames(tax.mps) <- c(sprintf("55+ años\n(n=%s)",ncol(SLL2.top15.senior)),
                           sprintf("19-55 años\n(n=%s)",ncol(SLL2.top15.adult)),
                           sprintf("7-18 años\n(n=%s)",ncol(SLL2.top15.teens)))
    
    if (!is.null(lab) && lab %in% c("Fibrosis Quística","Síndrome de Down")) tax.mps <- tax.mps[,2:3]
    # to order bars appropriately
    if (toupper(orientation) %in% c("VERT","VERTICAL","V")) tax.mps <- tax.mps[ , rev(colnames(tax.mps))]
    
    # tax.sd <- cbind(rowSds(as.matrix(SLL2.top15.region)), rowSds(as.matrix(SLL2.top15.disorder)), rep(0,15))
    tax.mps <- rbind(apply(tax.mps,2, function(x) 100-sum(x)), tax.mps)
    rownames(tax.mps)[1] <- "Otros"
    
    tax.mps.m <- reshape2::melt(tax.mps)
    colnames(tax.mps.m) <- c("Genus","Samples","value")
    
    
  }
  # *************************** #
  
  
  getPalette = colorRampPalette(brewer.pal(12,"Paired"))
  # colPal.set2=rev(getPalette(19))
  colPal=c("white",rev(getPalette(15)))
  
  if (toupper(orientation) %in% c("VERT","VERTICAL","V")) {
    
    stacked <- ggplot(tax.mps.m, aes(x=Samples, y=value, fill=Genus)) +
      geom_bar(aes(x=Samples, y=value, fill=Genus), stat="identity", color="black") +
      theme_minimal() + #ylim(0,100) + #scale_y_log10() +
      scale_fill_manual(values=colPal, name="Género") + 
      guides(fill = useGuide) +
      theme(plot.title = element_text(hjust=0.5, size=20), line = lel, 
            axis.title.x = element_text(size=17), axis.title.y = element_text(size=14),
            axis.text.y = aty, axis.text.x = atx, #element_text(size=15, angle=ang, hjust = 1),
            legend.title = element_text(size=17), legend.text = element_text(size=15)) +
      ggtitle(sprintf('Abundancia relativa de los géneros más comunes%s', title.extra)) +
      xlab('') + ylab('Abundancia media normalizada por muestra (%)') #scale_fill_hue(name=tlev)
    
    orient.file <- "vertical"
    
  } else if (toupper(orientation) %in% c("HOR","HORIZ","HORIZONTAL","H")) {
    
    stacked <- ggplot(tax.mps.m, aes(x=Samples, y=value, fill=Genus)) +
      geom_bar(aes(x=Samples, y=value, fill=Genus), stat="identity", color="black") +
      coord_flip() +
      theme_minimal() + #ylim(0,100) + #scale_y_log10() +
      scale_fill_manual(values=colPal, name="Género") + 
      guides(fill = useGuide) +
      theme(plot.title = element_text(hjust=0.5, size=20), line = lel, 
            axis.title.x = element_text(size=17), axis.title.y = element_text(size=14),
            axis.text.y = aty, axis.text.x = atx, #element_text(size=15),# angle=45),
            legend.title = element_text(size=17), legend.text = element_text(size=15)) +
      ggtitle(sprintf('Abundancia relativa de los géneros más comunes%s', title.extra)) +
      xlab('') + ylab('Abundancia media normalizada por muestra (%)') #scale_fill_hue(name=tlev)
    
    orient.file <- "horizontal"
    
  }
  
  
  # save plot to file
  if (saveStacked == T) {
    
    dir <- sprintf("%s/figures/Abundances/stacked_bars/%s", p2_dir, orient.file)
    if ( ! is.null(plotDir) )
      dir <- plotDir
    
    stacked +
      ggsave(sprintf("%s/%s.png", dir, lab.file), 
             width = 16.82, height = 8.41, device = "png")
  }
  
}
# ******************************************************************************************************************* #



# ******************************************************************************************* #
library(RColorBrewer)
plot_stacked_bars.disorders <- function(disorder, samps.dis, samps.conts, tl="Genus", glomTab=gloms_rel, 
                                        nTop=15, orientation="vertical", ang=0,
                                        saveStacked=F, plotDir=NULL) {
  
  # for most plots, will not show legend, but will just for the individual plots
  useGuide <- FALSE
  # and will have x tick labels for all plots except individuals, which will be blank
  if (toupper(orientation) %in% c("VERT","VERTICAL","V")) atx <- element_text(size=25, angle=ang)#, hjust = 1)
  else atx <- element_text(size=15)
  
  # *************************** #
  # top 10 most common OTUs
  top15 <- rev(names(sort(rowSums(glomTab[[ tl ]]), decreasing = T))[1:nTop])
  SLL2.top15 <- as.data.frame(glomTab[[ tl ]])[top15,]
  
  
  tax.mps <- cbind(rowMeans(SLL2.top15[ , samps.conts]), rowMeans(SLL2.top15[ , samps.dis ]) )
  colnames(tax.mps) <- c("Matched controls", ifelse(disorder=="Downs_Syndrome", "Down's Syndrome", 
                                                    ifelse(disorder=="Cystic_fibrosis","Cystic Fibrosis", "NA")))
  
  # tax.sd <- cbind(rowSds(as.matrix(SLL2.top15.region)), rowSds(as.matrix(SLL2.top15.disorder)), rep(0,15))
  tax.mps <- rbind(apply(tax.mps, 2, function(x) 100-sum(x)), tax.mps)
  rownames(tax.mps)[1] <- "Other"
  
  # to order bars appropriately
  if (toupper(orientation) %in% c("VERT","VERTICAL","V")) tax.mps <- tax.mps[ , rev(colnames(tax.mps))]
  
  tax.mps.m <- reshape2::melt(tax.mps)
  colnames(tax.mps.m) <- c("Genus","Samples","value")
  
  
  # *************************** #
  
  getPalette = colorRampPalette(brewer.pal(12,"Paired"))
  # colPal.set2=rev(getPalette(19))
  colPal=c("white", rev(getPalette(15)))
  
  if (toupper(orientation) %in% c("VERT","VERTICAL","V")) {
    
    stacked <- ggplot(tax.mps.m, aes(x=Samples, y=value, fill=Genus)) +
      geom_bar(aes(x=Samples, y=value, fill=Genus), stat="identity", color="black") +
      theme_minimal() + #ylim(0,100) + #scale_y_log10() +
      scale_fill_manual(values=colPal, name=tl) + 
      # guides(fill = useGuide) +
      theme(plot.title = element_text(hjust=0.5, size=20),# line = lel, 
            axis.title.x = element_text(size=17), axis.title.y = element_text(size=25),
            axis.text.y = element_text(size=20), axis.text.x = atx, #element_text(size=15, angle=ang, hjust = 1),
            legend.title = element_text(size=25), legend.text = element_text(size=20)) +
      # ggtitle(sprintf('Abundancia relativa de los géneros más comunes%s', title.extra)) +
      xlab('') + ylab('Mean relative abundance (%)') #scale_fill_hue(name=tlev)
    
    orient.file <- "vertical"
    
  } else if (toupper(orientation) %in% c("HOR","HORIZ","HORIZONTAL","H")) {
    
    stacked <- ggplot(tax.mps.m, aes(x=Samples, y=value, fill=Genus)) +
      geom_bar(aes(x=Samples, y=value, fill=Genus), stat="identity", color="black") +
      coord_flip() +
      theme_minimal() + #ylim(0,100) + #scale_y_log10() +
      scale_fill_manual(values=colPal, name=tl) + 
      guides(fill = useGuide) +
      theme(plot.title = element_text(hjust=0.5, size=20), line = lel, 
            axis.title.x = element_text(size=17), axis.title.y = element_text(size=15),
            axis.text.y = element_text(size=17), axis.text.x = atx, #element_text(size=15),# angle=45),
            legend.title = element_text(size=17), legend.text = element_text(size=15)) +
      # ggtitle(sprintf('Abundancia relativa de los géneros más comunes%s', title.extra)) +
      xlab('') + ylab('Mean relative abundance per sample (%)') #scale_fill_hue(name=tlev)
    
    orient.file <- "horizontal"
    
  }
  
  # *************************** #
  # save plot to file
  if (saveStacked == T) {
    
    ggsave(sprintf("%s/%s.png", plotDir, lab.file), 
           width = 16.82, height = 8.41, device = "png", plot = stacked)
    
  } else {
    
    stacked 
    
  }
  
  # *************************** #
  
}
# ******************************************************************************************* #





# ****************************************************************************************************************** ####
# ****************************************************************************************************************** ####














