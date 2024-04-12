#!/usr/bin/env Rscript

# Script to run HZAR with hybrid index data converted from gghybrid output

## Load the package
library(hzar)
library(scales)

#Set working directory
setwd("/path/to/hybrid_index_cline")

# List of locus ID of just the locus IDs and chromosome coordinates to tell hzar what to loop over
# These are the first 4 columns of the `populations.sumstats.tsv`from Stacks
snp_ID_file <- "snp_ids.tsv"

#This is the hzar file output by stacks
hzar_input_file <- "HI_Long_p8.hzar.csv"

# Name output files

name <- "HI_Long_p8"

##Read in data to R
## Also note that R is adding an 'X' to the beginning of each locus ID because they are numerical and it wants a letter there
manakinRAD <- read.csv(hzar_input_file, comment.char = '#', stringsAsFactors = FALSE)
#Modify the population names for sorting
manakinRAD$Population <- sapply(strsplit(manakinRAD$Population, '_'), function(manakinRAD){manakinRAD[2]})

# Sort by the numerical population ID
# Make sure that the population IDs are numerical so they get sorted in the right order
manakinRAD <- manakinRAD[order(manakinRAD$Population, decreasing = TRUE),]

## Population distances
cg_100 <- 00.00
#mr_095 <- 08.44
pr_090 <- 20.75
ru_080 <- 29.50
qp_060 <- 42.50
ro_050 <- 48.50
fc_040 <- 71.25
so_030 <- 79.00
ss_020 <- 92.50
dists <- c(cg_100, pr_090, ru_080, qp_060, ro_050, fc_040, so_030, ss_020)

# Add these distances to the manakinRAD data
manakinRAD$Distance <- dists

#Get a list of just the locus IDs and chromosome coordinates
# These are the first 4 columns of the `populations.sumstats.tsv`from Stacks
lociID <- read.delim(snp_ID_file, header=FALSE, stringsAsFactors = FALSE)
lociID$V5 <- paste(lociID$V1,'_',lociID$V4,sep='')
lociID$V6 <- paste('X',lociID$V5,sep='')
colnames(lociID) <- c('locus','chromosome','basepair','column','snp_id','hzar_id')

# Vectors for final per-snp dataframe
snps        <- c()
center.km   <- c()
width.km    <- c()
p.min       <- c()
p.max       <- c()
best.model  <- c()
scaffold    <- c()
basepair    <- c()
center.low  <- c()
center.med  <- c()
center.high <- c()
width.low   <- c()
width.med   <- c()
width.high  <- c()

## Distance Offset
offset=20
minDist=0
maxDist=95

# Dataframe with all observed clines
## The loop will add the values for the other clines
kms <- (minDist-offset):(maxDist+offset)
all_clines <- data.frame(kms)

## A typical chain length.  This value is the default setting in the package.
chainLength=1e5;

## Make each model run off a separate seed
mainSeed=
  list(A=c(596,528,124,978,544,99),
       B=c(528,124,978,544,99,596),
       C=c(124,978,544,99,596,528))

# Colors for the plots
# Re-add the observed data for the sites, with colors and legends
colors <- c('cornflowerblue', # CG
            #'coral1',         # MR
            'burlywood4',     # PR
            'deeppink2',      # RU
            'chartreuse3',    # QP
            'darkgoldenrod3', # RO
            'darkorchid3',    # FC
            'mediumseagreen', # SO
            'firebrick3')     # SS

#pop_ids <- c('10','9.5','9','8','6','5','4','3','2')
pop_ids <- c('10','9','8','6','5','4','3','2')

# Save plots in a PDF
#pdf('./individual_snp_clines.pdf',width=4,height=4)

# Loop over all the loci
for (i in 1:nrow(lociID)) {
	#for (i in 1:2) {
	# i=1
  
  # Select the current SNP
  curr.snp <- lociID[i,]
  snp_id   <- curr.snp$snp_id
  loc      <- curr.snp$hzar_id
  chrom    <- curr.snp$chromosome
  bp       <- curr.snp$basepair
  
  message('## ----- PROCESSING LOCUS ', i, ' - ', snp_id, ' ----- ##\n')

  ## Molecular Analysis
  
  ## Blank out space in memory to hold molecular analysis
  if(length(apropos("^mkn$",ignore.case=FALSE)) == 0 ||
     !is.list(mkn) ) mkn <- list()
  ## We are doing just the one allele at one locus, but it is
  ## good to stay organized.
  
  # ARC: Instead of using the literal string of the locus ID,  we  are using 
  # the variable `loc` and indexing with it
  mkn[[loc]] <- list()
  ## Space to hold the observed data
  mkn[[loc]][['obs']] <- list()
  ## Space to hold the models to fit
  mkn[[loc]][['models']] <- list()
  ## Space to hold the compiled fit request
  mkn[[loc]][['fitRs']] <- list()
  ## Space to hold the output data chains
  mkn[[loc]][['runs']] <- list()
  ## Space to hold the analyzed data
  mkn[[loc]][['analysis']] <- list()
  
  ## Load the locus of interest
  locA <- paste(loc,'.A',sep='')
  locN <- paste(loc,'.N',sep='')
  
  mkn[[loc]][['obs']] <- 
    hzar.doMolecularData1DPops(manakinRAD$Distance,
                               manakinRAD[,locA],
                               manakinRAD[,locN])
  
  ## Make a helper function
  mkn.loadLocAmodel <- function(scaling,tails,
                                id=paste(scaling,tails,sep="."))
    mkn[[loc]]$models[[id]] <<- hzar.makeCline1DFreq(mkn[[loc]]$obs, scaling, tails)
  
  mkn.loadLocAmodel("fixed","none","modelI");
  mkn.loadLocAmodel("free" ,"none","modelII");
  mkn.loadLocAmodel("free" ,"both","modelIII");
  
  ## Modify all models to focus on the region where the observed
  ## data were collected.
  ## Observations were between 0 and about 80 km.
  mkn[[loc]]$models <- sapply(mkn[[loc]]$models,
                              hzar.model.addBoxReq,
                              minDist-offset, maxDist+offset,
                              simplify=FALSE)
  
  ## Compile each of the models to prepare for fitting
  mkn[[loc]]$fitRs[['init']] <- sapply(mkn[[loc]]$models,
                                       hzar.first.fitRequest.old.ML,
                                       obsData=mkn[[loc]]$obs,
                                       verbose=FALSE,
                                       simplify=FALSE)
  ## Update the settings for the fitter if desired.
  ### Model I
  mkn[[loc]]$fitRs$init$modelI$mcmcParam$chainLength   <- chainLength        # 1e5 by default
  mkn[[loc]]$fitRs$init$modelI$mcmcParam$burnin        <- chainLength %/% 10 # 1e4 by default
  mkn[[loc]]$fitRs$init$modelI$mcmcParam$seed[[1]]     <- mainSeed$A
  ### Model II
  mkn[[loc]]$fitRs$init$modelII$mcmcParam$chainLength  <- chainLength        # 1e5 by default
  mkn[[loc]]$fitRs$init$modelII$mcmcParam$burnin       <- chainLength %/% 10 # 1e4 by default
  mkn[[loc]]$fitRs$init$modelII$mcmcParam$seed[[1]]    <- mainSeed$B
  ### Model III
  mkn[[loc]]$fitRs$init$modelIII$mcmcParam$chainLength <- chainLength        # 1e5 by default
  mkn[[loc]]$fitRs$init$modelIII$mcmcParam$burnin      <- chainLength %/% 10 # 1e4 by default
  mkn[[loc]]$fitRs$init$modelIII$mcmcParam$seed[[1]]   <- mainSeed$C
  
  # RUN MODEL I
  ## Run just one of the models for an initial chain
  mkn[[loc]]$runs$init <- list()
  mkn[[loc]]$runs$init$modelI <-
    hzar.doFit(mkn[[loc]]$fitRs$init$modelI)
  ## Plot the trace
  # plot(hzar.mcmc.bindLL(mkn[[loc]]$runs$init$modelI))
  
  # RUN MODEL II
  ## Run another model for an initial chain
  mkn[[loc]]$runs$init$modelII <-
    hzar.doFit(mkn[[loc]]$fitRs$init$modelII)
  ## Plot the trace
  # plot(hzar.mcmc.bindLL(mkn[[loc]]$runs$init$modelII))
  
  # RUN MODEL III
  ## Run another model for an initial chain
  mkn[[loc]]$runs$init$modelIII <-
    hzar.doFit(mkn[[loc]]$fitRs$init$modelIII)
  ## Plot the trace
  # plot(hzar.mcmc.bindLL(mkn[[loc]]$runs$init$modelIII))
  
  ## Compile a new set of fit requests using the initial chains 
  mkn[[loc]]$fitRs$chains <-
    lapply(mkn[[loc]]$runs$init,
           hzar.next.fitRequest)
  
  ## Replicate each fit request 3 times, keeping the original
  ## seeds while switching to a new seed channel.
  mkn[[loc]]$fitRs$chains <-
    hzar.multiFitRequest(mkn[[loc]]$fitRs$chains,
                         each=3,
                         baseSeed=NULL)
  
  ## Go ahead and run a chain of 3 runs for every fit requester
  ## ARC: This step can take quite a while to run because it looks for the models to
  ## converge. The code doesn't seem to behave properly when they do not. Might need
  ## to incorporate a check to stop the function after X amount of time when that happens.
  ## NOTE: This part is needed. 
  mkn[[loc]]$runs$chains <- hzar.doChain.multi(mkn[[loc]]$fitRs$chains,
                                               doPar=FALSE,
                                               inOrder=FALSE,
                                               count=3)
  
  ## Start aggregation of data for analysis
  
  ## Create a model data group for the null model (expected allele
  ## frequency independent of distance along cline) to include in
  ## analysis.
  
  mkn[[loc]]$analysis$initDGs <- list(
    nullModel =  hzar.dataGroup.null(mkn[[loc]]$obs))
  
  ## Create a model data group (hzar.dataGroup object) for each
  ## model from the initial runs.
  mkn[[loc]]$analysis$initDGs$modelI <-
    hzar.dataGroup.add(mkn[[loc]]$runs$init$modelI)
  mkn[[loc]]$analysis$initDGs$modelII <-
    hzar.dataGroup.add(mkn[[loc]]$runs$init$modelII)
  mkn[[loc]]$analysis$initDGs$modelIII <-
    hzar.dataGroup.add(mkn[[loc]]$runs$init$modelIII)
  
  ## Create a hzar.obsDataGroup object from the four hzar.dataGroup
  ## just created, copying the naming scheme (nullModel, modelI,
  ## modelII, modelIII).
  mkn[[loc]]$analysis$oDG <-
    hzar.make.obsDataGroup(mkn[[loc]]$analysis$initDGs)
  mkn[[loc]]$analysis$oDG <-
    hzar.copyModelLabels(mkn[[loc]]$analysis$initDGs,
                         mkn[[loc]]$analysis$oDG)
  
  ## Convert all 27 runs to hzar.dataGroup objects, adding them to
  ## the hzar.obsDataGroup object.
  mkn[[loc]]$analysis$oDG <-
    hzar.make.obsDataGroup(lapply(mkn[[loc]]$runs$chains,
                                  hzar.dataGroup.add),
                           mkn[[loc]]$analysis$oDG)
  
  ## Do model selection based on the AICc scores
  mkn[[loc]]$analysis$AICcTable <-
       hzar.AICc.hzar.obsDataGroup(mkn[[loc]]$analysis$oDG)
  
  ## Print out the model with the minimum AICc score
  best_model <- print(mkn[[loc]]$analysis$model.name <-
    rownames(mkn[[loc]]$analysis$AICcTable
    )[[ which.min(mkn[[loc]]$analysis$AICcTable$AICc )]])
  
  ## Extract the hzar.dataGroup object for the selected model
  mkn[[loc]]$analysis$model.selected <-
    mkn[[loc]]$analysis$oDG$data.groups[[mkn[[loc]]$analysis$model.name]]
  
  # If working with non-null models
  if (best_model != "nullModel"){
 
    ## Print the maximum likelihood cline for the selected model
    ## Need to capture this maximum likelihood output
    max_likelihood_cline <- hzar.get.ML.cline(mkn[[loc]]$analysis$model.selected)
    center <- max_likelihood_cline$param.all$center
    width  <- max_likelihood_cline$param.all$width
    pMin   <- max_likelihood_cline$param.all$pMin
    pMax   <- max_likelihood_cline$param.all$pMax
    # print(max_likelihood_cline)
    
    ## Obtain the confidence intervals
    conf        <- hzar.qScores.dataGroup(mkn[[loc]]$analysis$model.selected)
    center_low  <- conf[1,]$center
    center_med  <- conf[2,]$center
    center_high <- conf[3,]$center
    width_low   <- conf[1,]$width
    width_med   <- conf[2,]$width
    width_high  <- conf[3,]$width
    
    # @ARC: If you need to plot the individual clines, they are stored in
    #cline_line <- mkn[[loc]]$analysis$model.selected$obsData$frame
    cline_line <- max_likelihood_cline$clineFunc((minDist-offset):(maxDist+offset))
    ## Store current cline in all clines table
    all_clines[,snp_id] <- cline_line
  } else {
    center      <- NA
    width       <- NA
    pMin        <- NA
    pMax        <- NA
    center_low  <- NA
    center_med  <- NA
    center_high <- NA
    width_low   <- NA
    width_med   <- NA
    width_high  <- NA
  }
  
## This section makes graphs of individual snps. If you want one of these graphs, I suggest not running the loop and running a single iteration of the snp to get the graph
  # # Make the main cline plot
  # par(mar=c(5.1, 4.1, 4.1, 4.5), xpd=TRUE)
  # ## Plot the 95% credible cline region for the selected model
  # tryCatch(
  #   hzar.plot.fzCline(mkn[[loc]]$analysis$model.selected,
  #                     ylab='Allele Frequency',
  #                     xlab='Distance from M. vitellinus (Km)',
  #                     main=paste('SNP location:\n',chrom,bp),
  #                     ylim=c(0,1),
  #                     xlim=c(0,100),
  #                     las=1
  #   ),
  #   error = function(e) {
  #     message('ERROR: hzar.plot.fzCline (1st spot) ', i, ' - ', snp_id, ' - ', e)
  #   }
  # )
  # 
  # # Re-add the main cline line, but thicker
  # tryCatch(
  #   hzar.plot.cline(mkn[[loc]]$analysis$model.selected,
  #                   add=TRUE,
  #                   lwd=2),
  #   error = function(e) {
  #     message('ERROR: hzar.plot.cline (2nd spot) ', i, ' - ', snp_id, ' - ', e)
  #   }
  # )
  # ## Look at a graph of the observed data
  # tryCatch(
  #   hzar.plot.obsData(mkn[[loc]]$obs,
  #                     add=TRUE,
  #                     col='black',
  #                     bg=colors,
  #                     pch=21,
  #                     cex=1.25),
  #   error = function(e) {
  #     message('ERROR: hzar.plot.obsData (3rd spot) ', i, ' - ', snp_id, ' - ', e)
  #   }
  # )
  # 
  # # Add legend
  # legend('right',
  #        inset=c(-0.3,0),
  #        legend=pop_ids,
  #        pch=21,
  #        pt.bg=colors)
  
  pdf(paste("./", name ,".pdf", sep = ""), width=4, height=4)
  
  hzar.plot.fzCline(mkn[[loc]]$analysis$model.selected,
                    ylab='Hybrid Index',
                    xlab='Distance from M. vitellinus (Km)',
                    main='Hybrid Index',
                    ylim=c(0,1),
                    xlim=c(0,100),
                    las=1
  );
  dev.off()
  
  


  ## Add data to DF
  center.km   <- c(center.km, center)
  width.km    <- c(width.km, width)
  best.model  <- c(best.model, best_model)
  snps        <- c(snps, snp_id)
  p.min       <- c(p.min, pMin)
  p.max       <- c(p.max, pMax)
  scaffold    <- c(scaffold, chrom)
  basepair    <- c(basepair,bp)
  center.low  <- c(center.low, center_low)
  center.med  <- c(center.med, center_med)
  center.high <- c(center.high, center_high)
  width.low   <- c(width.low, width_low)
  width.med   <- c(width.med, width_med)
  width.high  <- c(width.high, width_high)
}

# Close plot
#.=dev.off()

# Create Final dataframe
cline_table <- data.frame(scaffold,
                          basepair,
                          snps,
                          center.km,
                          width.km,
                          p.min,
                          p.max,
                          best.model,
                          center.low,
                          center.med,
                          center.high,
                          width.low,
                          width.med,
                          width.high)

# Save the cline table dataframe into a file
write.table(cline_table,
            file=paste("./", name ,"_snp_clines_parameters.tsv", sep = ""),
            quote=FALSE,
            sep='\t',
            eol='\n',
            row.names = FALSE,
            col.names = colnames(cline_table))

# Save the all clines dataframe into a file
write.table(all_clines,
            file=paste("./", name ,"_all_clines.tsv", sep = ""),
            quote=FALSE,
            sep='\t',
            eol='\n',
            row.names = FALSE,
            col.names = colnames(all_clines))
