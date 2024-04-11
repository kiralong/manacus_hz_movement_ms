#!/usr/bin/env Rscript

# Script to run R package HZAR in parallel for geographic cline analysis

## Load the packages
library(optparse)
library(parallelly)
library(foreach)
library(doParallel)
library(hzar)
library(scales)

## Prase command line options
opt = parse_args(OptionParser(option_list=list(
  make_option(c("-w", "--workdir"), type="character", default="./",
              help="Set the R working directory", metavar="character"),
  make_option(c("-s", "--snpidfile"), type="character", default=NULL,
              help="List of snp IDs (tsv file)", metavar="character"),
  make_option(c("-z", "--hzarfile"), type="character", default=NULL,
              help="hzar .csv file from Stacks", metavar="character"),
  make_option(c("-c", "--chromosome"), type="character", default=NULL,
              help="Target chromosome or scaffold of interest", metavar="character"),
  make_option(c("-t", "--cores"), type="integer",  default= parallelly::availableCores(omit = 1) ,
              help="Number of cores to use for multicore parallelism")
)))

# Set working directory
setwd(opt$workdir)

# File with the list of the locus IDs and chromosome coordinates to tell hzar what to loop over
# These are the first 4 columns of the `populations.sumstats.tsv`from Stacks
snp_ID_file <- opt$snpidfile

# Hzar data file - output by stacks
hzar_input_file <- opt$hzarfile

# Target chromosome for this run
target_chromosome <- opt$chromosome

# TODO make the default in favor of explicit flag value
cl <- parallelly::makeClusterPSOCK(opt$cores, autoStop = TRUE)
registerDoParallel(cl)

## Read in data to R
load_hzar_input_data <- function(hzar_input_file) {
  manakinRAD <- read.csv(hzar_input_file, comment.char = '#', stringsAsFactors = FALSE)
  # Modify the population names for sorting
  manakinRAD$Population <- sapply(strsplit(manakinRAD$Population, '_'), function(manakinRAD){manakinRAD[2]})

  # Sort by the numerical population ID
  # Make sure that the population IDs are numerical so they get sorted in the right order
  manakinRAD <- manakinRAD[order(manakinRAD$Population, decreasing = TRUE),]

  ## Population distances
  # Add these distances to the manakinRAD data
  cg_100 <- 00.00
  mr_095 <- 08.44
  pr_090 <- 20.75
  ru_080 <- 29.50
  qp_060 <- 42.50
  ro_050 <- 48.50
  fc_040 <- 71.25
  so_030 <- 79.00
  ss_020 <- 92.50

  # Check the number of populations in the input file and adjust accordingly
  if (nrow(manakinRAD) == 9){
    dists <- c(cg_100, mr_095, pr_090, ru_080, qp_060, ro_050, fc_040, so_030, ss_020)
  } else if (nrow(manakinRAD) == 8){
    dists <- c(cg_100, pr_090, ru_080, qp_060, ro_050, fc_040, so_030, ss_020)
  }

  # Add these distances to the manakinRAD data
  manakinRAD$Distance <- dists

  return(manakinRAD)
}

manakinRAD <- load_hzar_input_data(hzar_input_file)

# Get a list of just the locus IDs and chromosome coordinates
# These are the first 4 columns of the `populations.sumstats.tsv`from Stacks
load_loci_table <- function(snp_ID_file,target_chromosome) {
  loci_table <- read.delim(snp_ID_file, header=FALSE, stringsAsFactors = FALSE)
  loci_table$V5 <- paste(loci_table$V1,'_',loci_table$V4,sep='')
  loci_table$V6 <- paste('X',loci_table$V5,sep='')
  colnames(loci_table) <- c('locus','chromosome','basepair','column','snp_id','hzar_id')

  # Filter the SNP ids to just include the target chromosomes
  loci_table <- subset(loci_table, chromosome == target_chromosome)
  message('Running HZAR for Chr ', target_chromosome, ', total loci: ', nrow(loci_table))

  return(loci_table)
}

loci_table <- load_loci_table(snp_ID_file, target_chromosome)

## Distance Offset
offset=20
minDist=0
maxDist=95
kms <- (minDist-offset):(maxDist+offset)

# Hzar Molecular Analysis function to be run on each loci.
# Returns a vector with the results of a single snp
doHzarMolecularAnalysis <- function(
  current_loci_table_row,
  progress_index,  # Simply for logging, progress, and readability Has no other function.
  manakinRAD,
  minDistWithOffset,
  maxDistWithOffset,
  kms
) {
  snp_id   <- current_loci_table_row[['snp_id']]
  # Reminder: hzar_id is ahzar-provided locus ID which is prefixed with an 'X'
  locus_id <- current_loci_table_row[['hzar_id']]
  chrom    <- current_loci_table_row[['chromosome']]
  bp       <- current_loci_table_row[['basepair']]

  # TODO Figure out parallel logging
  message(paste('## ----- PROCESSING LOCUS', progress_index, snp_id, ' ----- ##\n', sep = ' - '))

  ## A typical chain length. This value is the default setting in the package.
  chainLength=1e5;

  ## Make each model run off a separate seed
  mainSeed=list(
    A=c(596,528,124,978,544,99),
    B=c(528,124,978,544,99,596),
    C=c(124,978,544,99,596,528)
  )

  ## Molecular Analysis

  ## Space to hold the models to fit
  models <- list()
  ## Space to hold the compiled fit request
  fitRs <- list()
  ## Space to hold the output data chains
  runs <- list()
  ## Space to hold the analyzed data
  analysis <- list()

  ## Load the locus of interest (observed data)
  locA <- paste(locus_id,'.A',sep='')
  locN <- paste(locus_id,'.N',sep='')
  obs <- hzar.doMolecularData1DPops(
    manakinRAD$Distance,
    manakinRAD[,locA],
    manakinRAD[,locN]
  )

  # Helper function to set up the models
  loadLocAmodel <- function(scaling, tails, id=paste(scaling,tails,sep=".")) {
    models[[id]] <<- hzar.makeCline1DFreq(obs, scaling, tails)
  }
  loadLocAmodel("fixed","none","modelI");
  loadLocAmodel("free" ,"none","modelII");
  loadLocAmodel("free" ,"both","modelIII");

  ## Modify all models to focus on the region where the observed
  ## data were collected.
  ## Observations were between 0 and about 80 km.
  models <- sapply(
    models,
    hzar.model.addBoxReq,
    minDistWithOffset,
    maxDistWithOffset,
    simplify=FALSE
  )

  ## Compile each of the models to prepare for fitting
  fitRs[['init']] <- sapply(
    models,
    hzar.first.fitRequest.old.ML,
    obsData=obs,
    verbose=FALSE,
    simplify=FALSE
  )

  ## Update the settings for the fitter if desired.
  ### Model I
  fitRs$init$modelI$mcmcParam$chainLength   <- chainLength        # 1e5 by default
  fitRs$init$modelI$mcmcParam$burnin        <- chainLength %/% 10 # 1e4 by default
  fitRs$init$modelI$mcmcParam$seed[[1]]     <- mainSeed$A
  ### Model II
  fitRs$init$modelII$mcmcParam$chainLength  <- chainLength        # 1e5 by default
  fitRs$init$modelII$mcmcParam$burnin       <- chainLength %/% 10 # 1e4 by default
  fitRs$init$modelII$mcmcParam$seed[[1]]    <- mainSeed$B
  ### Model III
  fitRs$init$modelIII$mcmcParam$chainLength <- chainLength        # 1e5 by default
  fitRs$init$modelIII$mcmcParam$burnin      <- chainLength %/% 10 # 1e4 by default
  fitRs$init$modelIII$mcmcParam$seed[[1]]   <- mainSeed$C

  runs$init <- list()

  # RUN MODEL I
  ## Run just one of the models for an initial chain
  runs$init$modelI <- hzar.doFit(fitRs$init$modelI)

  # RUN MODEL II
  ## Run another model for an initial chain
  runs$init$modelII <- hzar.doFit(fitRs$init$modelII)

  # RUN MODEL III
  ## Run another model for an initial chain
  runs$init$modelIII <- hzar.doFit(fitRs$init$modelIII)

  ## Compile a new set of fit requests using the initial chains
  fitRs$chains <- lapply(runs$init, hzar.next.fitRequest)

  ## Replicate each fit request 3 times, keeping the original
  ## seeds while switching to a new seed channel.
  fitRs$chains <- hzar.multiFitRequest(fitRs$chains, each=3, baseSeed=NULL)

  ## Run a chain of 3 runs for every fit requester
  runs$chains <- hzar.doChain.multi(fitRs$chains, doPar=FALSE, inOrder=FALSE, count=3)

  ## Start aggregation of data for analysis

  ## Create a model data group for the null model (expected allele
  ## frequency independent of distance along cline) to include in
  ## analysis.

  analysis$initDGs <- list(nullModel = hzar.dataGroup.null(obs))

  ## Create a model data group (hzar.dataGroup object) for each
  ## model from the initial runs.
  analysis$initDGs$modelI <- hzar.dataGroup.add(runs$init$modelI)
  analysis$initDGs$modelII <- hzar.dataGroup.add(runs$init$modelII)
  analysis$initDGs$modelIII <- hzar.dataGroup.add(runs$init$modelIII)

  ## Create a hzar.obsDataGroup object from the four hzar.dataGroup
  ## just created, copying the naming scheme (nullModel, modelI,
  ## modelII, modelIII).
  analysis$oDG <- hzar.make.obsDataGroup(analysis$initDGs)
  analysis$oDG <- hzar.copyModelLabels(analysis$initDGs, analysis$oDG)

  ## Convert all 27 runs to hzar.dataGroup objects, adding them to
  ## the hzar.obsDataGroup object.
  analysis$oDG <- hzar.make.obsDataGroup(
    lapply(runs$chains, hzar.dataGroup.add),
    analysis$oDG
  )

  ## Do model selection based on the AICc scores
  analysis$AICcTable <- hzar.AICc.hzar.obsDataGroup(analysis$oDG)

  ## Assign and print out the model with the minimum AICc score
  best_model <- print(
    analysis$model.name <- rownames(analysis$AICcTable)[[which.min(analysis$AICcTable$AICc)]]
  )

  ## Extract the hzar.dataGroup object for the selected model
  analysis$model.selected <- analysis$oDG$data.groups[[analysis$model.name]]

  # If working with non-null models
  if (best_model != "nullModel"){

    ## The maximum likelihood cline for the selected model
    ## Need to capture this maximum likelihood output
    max_likelihood_cline <- hzar.get.ML.cline(analysis$model.selected)
    center <- max_likelihood_cline$param.all$center
    width  <- max_likelihood_cline$param.all$width
    pMin   <- max_likelihood_cline$param.all$pMin
    pMax   <- max_likelihood_cline$param.all$pMax

    ## Obtain the confidence intervals
    conf        <- hzar.qScores.dataGroup(analysis$model.selected)
    center_low  <- conf[1,]$center
    center_med  <- conf[2,]$center
    center_high <- conf[3,]$center
    width_low   <- conf[1,]$width
    width_med   <- conf[2,]$width
    width_high  <- conf[3,]$width

    cline_line <- max_likelihood_cline$clineFunc(kms)
  } else {
    # When null model
    # NOTICE NOT EVERY VALUE gets wiped out - some are preserved, such as best_model
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

    cline_line <- rep(NA, length(kms))
  }

  return(list(
    chromosome = chrom,
    snp_id = snp_id,
    base_pair = bp,
    best_model = best_model,
    center = center,
    width = width,
    pMin = pMin,
    pMax = pMax,
    center_low = center_low,
    center_med = center_med,
    center_high = center_high,
    width_low = width_low,
    width_med = width_med,
    width_high = width_high,
    cline_line = cline_line
  ))
}

# Process all the loci
hzar_results <- foreach(
    loci_table_row = split(loci_table, 1:nrow(loci_table)),
    progress_counter = 1:nrow(loci_table) # MUST be SAME length as loci_table
) %dopar% {
  # Provide packages to parallel sandbox environment
  library(hzar)
  library(scales)
  # Run analysis
  doHzarMolecularAnalysis(
    loci_table_row,
    progress_counter,
    # foreach/dopar can pull in expression values from the current scope
    manakinRAD,
    minDist - offset,
    maxDist + offset,
    kms
  )
}

# Create final dataframes

# Dataframe with all observed clines
all_clines_df <- data.frame(kms)

# Dataframe with snp clines parameters
snp_cline_parameters_column_names <- c(
  'chromosome',
  'snp_id',
  'base_pair',
  'best_model',
  'center',
  'width',
  'pMin',
  'pMax',
  'center_low',
  'center_med',
  'center_high',
  'width_low',
  'width_med',
  'width_high'
)
snp_cline_parameters_df <- setNames(
  as.data.frame(matrix(nrow = nrow(loci_table), ncol = length(snp_cline_parameters_column_names))),
  snp_cline_parameters_column_names
)

foreach(result = hzar_results, i = 1:length(hzar_results)) %do% {
  # Get out the cline data
  all_clines_df[,result$snp_id] <- result$cline_line

  # Get out the snp_cline_parameter data
  snp_cline_parameters_df[i,] <- c(
    result$chromosome,
    result$snp_id,
    result$base_pair,
    result$best_model,
    result$center,
    result$width,
    result$pMin,
    result$pMax,
    result$center_low,
    result$center_med,
    result$center_high,
    result$width_low,
    result$width_med,
    result$width_high
  )
}

# Save the cline table dataframe into a file
write.table(snp_cline_parameters_df,
            file='./snp_cline_parameters.tsv',
            quote=FALSE,
            sep='\t',
            eol='\n',
            row.names = FALSE,
            col.names = colnames(snp_cline_parameters_df))

# Save the all clines dataframe into a file
write.table(all_clines_df,
            file='./all_clines.tsv',
            quote=FALSE,
            sep='\t',
            eol='\n',
            row.names = FALSE,
            col.names = colnames(all_clines_df))
