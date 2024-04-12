## Hzar script for running on morphological data

## The package hzar was removed from CRAN in June 2022, so install from local 0.2-5 version tgz from archive
# Specify the filepath of the hzar tgz as the following environment variable
hzarLocalPackageLocation <- paste(Sys.getenv("HZAR_PACKAGE_LOCATION")[1], sep = '')
if (is.na(hzarLocalPackageLocation) || hzarLocalPackageLocation == '' ) {
  message("Did not find HZAR_PACKAGE_LOCATION environment variable set, so using hardcoded value.")
  hzarLocalPackageLocation = "/Users/hzar/hzar_0.2-5.tar.gz"
  
}

# Manually install dependencies because the "dependencies = TRUE" didn't seem to pull down the others automatically
install.packages("MCMCpack")
install.packages("foreach")
install.packages(hzarLocalPackageLocation,  repos = NULL, dependencies = TRUE, type = "source")


# Set working directory
#setwd("~/path/to/hzar/Morphological_data/2022_DC_measurement_redos/calibrated_spec_runs")

## Load packages
library(hzar)

# Hzar data file, a csv file
hzar_input_file <- "manacus_morphological_beard_epaulet_only.csv"

## Specify the trait you want to run

#trait <- "Tail"
#trait <- "epaulet"
trait <- "epaulet_width"
#trait <- "beard"
#trait <- "beard_length"
#trait <- "rel_belly_mean"
#trait <- "belly_mean"
#trait <- "throat_mean"
#trait <- "larger_hybrid_index"
#experiment <- "Brum"
experiment <- "Long"

## Read in data to R
manacus_morphology <- read.csv(hzar_input_file, comment.char = '#', stringsAsFactors = FALSE, na.strings = "None")
#manacus_morphology$rel_belly_mean <- manacus_morphology$belly_mean/max(manacus_morphology$belly_mean)

  
## Write out internal example trait data 
write.table(manacus_morphology,     # The data we just loaded
            file="mknExTrait.txt", # The file to overwrite
            col.names=TRUE,       # The columns are named
            row.names=FALSE,      # The rows are not named
            sep="\t",             # The file will be tab-delimited
            quote=TRUE)           # Use quotes as needed.

## As we no longer need the in-memory copy, drop the local reference
manakinMorphological <- NULL
hzar_distance_file <- "brum_manacusLocations.csv"
#hzar_distance_file <- "manacusLocations.csv"

## Write out internal example site data 
manacus_distances <- read.csv(hzar_distance_file, stringsAsFactors = TRUE)
#print(manacus_distances)
write.table(manacus_distances,     # The data we just loaded
            file="mknExSite.txt", # The file to overwrite
            col.names=TRUE,       # The columns are named
            row.names=FALSE,      # The rows are not named
            sep="\t",             # The file will be tab-delimited
            quote=TRUE)           # Use quotes as needed.

## As we no longer need the in-memory copy, drop the local reference
manakinLocations<- NULL


## Set chain length
chainLength=1e5;
#chainLength=1e6;

## Make each model run off a separate seed
mainSeed=
  list(
       A=c(978,544,99,596,528,124),
       B=c(544,99,596,528,124,978),
       C=c(99,596,528,124,978,544))

# Vectors for final per-trait dataframe
study        <- c()
morpho_trait <- c()
center_km    <- c()
width_km     <- c()
best_model   <- c()
center_low   <- c()
center_med   <- c()
center_high  <- c()
width_low    <- c()
width_med    <- c()
width_high   <- c()


## Trait Analysis

## Load example Quantitative trait data from the package
## Note that this is the individual data, with labels
## identifying the locality of each sample.
manakinMorphological <- read.table("mknExTrait.txt",header=TRUE)

## Load example locality data, matching each localility to
## a site ID and a transect distance.
manakinLocations <- read.table("mknExSite.txt",header=TRUE)

## Blank out space in memory to hold morphological analysis
if(length(apropos("^mkn$",ignore.case=FALSE)) == 0 ||
   !is.list(mkn) ) mkn <- list()
mkn[[trait]] <- list();
## Space to hold the observed data
mkn[[trait]]$obs <- list();
## Space to hold the models to fit
mkn[[trait]]$models <- list();
## Space to hold the compiled fit requests
mkn[[trait]]$fitRs <- list();
## Space to hold the output data chains
mkn[[trait]]$runs <- list();
## Space to hold the analysed data
mkn[[trait]]$analysis <- list();


## Setting up your Trait
mkn[[trait]]$obs <-
  hzar.doNormalData1DRaw(hzar.mapSiteDist(manakinLocations$LocalityID,
                                          manakinLocations$distance),
                         manakinMorphological$Locality,
                         manakinMorphological[,trait])

## Look at a graph of the observed data
 hzar.plot.obsData(mkn[[trait]]$obs);

## Make a helper function
mkn.loadTraitmodel <- function(scaling,tails,
                              id=paste(scaling,tails,sep=".")){
  mkn[[trait]]$models[[id]] <<-
    hzar.makeCline1DNormal(mkn[[trait]]$obs, tails)
  ## As there is no quick option for "fixed" scaling, and 
  ## some sites have a fair number of samples (> 20),
  ## fix the mean and variance of the left and right sides of
  ## the cline to values observed at sites.
  if (all(regexpr("fixed",scaling,ignore.case=TRUE) == 1 )){
    hzar.meta.fix(mkn[[trait]]$models[[id]])$muL <<- TRUE
    hzar.meta.fix(mkn[[trait]]$models[[id]])$muR <<- TRUE
    hzar.meta.fix(mkn[[trait]]$models[[id]])$varL <<- TRUE
    hzar.meta.fix(mkn[[trait]]$models[[id]])$varR <<- TRUE
  }
  ## cg_100 is the "left" side of the cline, so pull the
  ## initial values from there.
  hzar.meta.init(mkn[[trait]]$models[[id]])$muL <<-
    mkn[[trait]]$obs$frame["cg_100","mu"]
  hzar.meta.init(mkn[[trait]]$models[[id]])$varL <<-
    mkn[[trait]]$obs$frame["cg_100","var"]
  ## ss_020 is the "right" side of the cline, so pull the
  ## initial values from there.
  hzar.meta.init(mkn[[trait]]$models[[id]])$muR <<-
    mkn[[trait]]$obs$frame["ss_020","mu"]
  hzar.meta.init(mkn[[trait]]$models[[id]])$varR <<-
    mkn[[trait]]$obs$frame["ss_020","var"]

  
}
mkn.loadTraitmodel("fixed","none","modelI");
mkn.loadTraitmodel("free" ,"none","modelII");
mkn.loadTraitmodel("free" ,"both","modelIII");


## Modify all models to focus on the region where the observed
## data were collected.
## Observations were between 0 and 100 km.
mkn[[trait]]$models <- sapply(mkn[[trait]]$models,
                          hzar.model.addBoxReq,
                          -20 , 110,
                          simplify=FALSE)

## Due to the large number of free variables, it is prudent to
## reduce the tune setting of modelIII
hzar.meta.tune(mkn[[trait]]$models$modelIII)<-1.4

## Due to the large number of free variables, it is prudent to
## reduce the tune setting of modelIII
hzar.meta.tune(mkn[[trait]]$models$modelIII)<-1.2

## Compile each of the models to prepare for fitting
## Note that we are using hzar.first.fitRequest.gC for fitting
## guassian (aka "normal") clines.
mkn[[trait]]$fitRs$init <- sapply(mkn[[trait]]$models,
                         hzar.first.fitRequest.gC,
                         obsData=mkn[[trait]]$obs,
                         verbose=FALSE,
                         simplify=FALSE)
## Update the settings for the fitter if desired.
mkn[[trait]]$fitRs$init$modelI$mcmcParam$chainLength <-
  chainLength;                          #1e5
mkn[[trait]]$fitRs$init$modelI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
mkn[[trait]]$fitRs$init$modelI$mcmcParam$seed[[1]] <-
  mainSeed$A

mkn[[trait]]$fitRs$init$modelII$mcmcParam$chainLength <-
  chainLength;                          #1e5
mkn[[trait]]$fitRs$init$modelII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
mkn[[trait]]$fitRs$init$modelII$mcmcParam$seed[[1]] <-
  mainSeed$B 

mkn[[trait]]$fitRs$init$modelIII$mcmcParam$chainLength <-
  chainLength;                          #1e5
mkn[[trait]]$fitRs$init$modelIII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
mkn[[trait]]$fitRs$init$modelIII$mcmcParam$seed[[1]] <-
  mainSeed$C 

## Check fit request settings
# print(mkn[[trait]]$fitRs$init)


## Run just one of the models for an initial chain
mkn[[trait]]$runs$init <- list()
mkn[[trait]]$runs$init$modelI <-
  hzar.doFit(mkn[[trait]]$fitRs$init$modelI)

## Save the model I trace plot in a pdf
pdf(paste('./',experiment, '_', trait, 'p8_',chainLength,'_modelItrace.pdf', sep=''), width=12, height=12)

## Plot the trace
 plot(hzar.mcmc.bindLL(mkn[[trait]]$runs$init$modelI))

dev.off()
 
## Run another model for an initial chain
mkn[[trait]]$runs$init$modelII <-
  hzar.doFit(mkn[[trait]]$fitRs$init$modelII)

## Save the model II trace plot in a pdf
pdf(paste('./',experiment, '_', trait, 'p8_',chainLength,'_model2trace.pdf', sep=''), width=12, height=12)

## Plot the trace
 plot(hzar.mcmc.bindLL(mkn[[trait]]$runs$init$modelII))

 dev.off()
 

## Run another model for an initial chain
mkn[[trait]]$runs$init$modelIII <-
  hzar.doFit(mkn[[trait]]$fitRs$init$modelIII)

  
## Save the model III trace plot in a pdf
pdf(paste('./',experiment, '_', trait, 'p8_', chainLength, '_model3trace.pdf', sep=''), width=12, height=12)
  
## Plot the trace
 plot(hzar.mcmc.bindLL(mkn[[trait]]$runs$init$modelIII))
dev.off()

## Compile a new set of fit requests using the initial chains 
mkn[[trait]]$fitRs$chains <-
  lapply(mkn[[trait]]$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original
## seeds while switching to a new seed channel.
mkn[[trait]]$fitRs$chains <-
  hzar.multiFitRequest(mkn[[trait]]$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Just to be thorough, randomize the initial value for each fit

## runif(9,-30,600) center for modelI, modelII, modelIII
mkn[[trait]]$fitRs$chains[[1]]$modelParam$init["center"]=534.4945
mkn[[trait]]$fitRs$chains[[2]]$modelParam$init["center"]=231.9774 
mkn[[trait]]$fitRs$chains[[3]]$modelParam$init["center"]=422.3549
mkn[[trait]]$fitRs$chains[[4]]$modelParam$init["center"]=428.7615
mkn[[trait]]$fitRs$chains[[5]]$modelParam$init["center"]=357.8308
mkn[[trait]]$fitRs$chains[[6]]$modelParam$init["center"]=319.1885
mkn[[trait]]$fitRs$chains[[7]]$modelParam$init["center"]=214.0428
mkn[[trait]]$fitRs$chains[[8]]$modelParam$init["center"]=140.9664
mkn[[trait]]$fitRs$chains[[9]]$modelParam$init["center"]=537.0636

## runif(9,0,630) width for modelI, modelII, modelIII
mkn[[trait]]$fitRs$chains[[1]]$modelParam$init["width"]=120.29368 
mkn[[trait]]$fitRs$chains[[2]]$modelParam$init["width"]=498.97202
mkn[[trait]]$fitRs$chains[[3]]$modelParam$init["width"]=387.67701 
mkn[[trait]]$fitRs$chains[[4]]$modelParam$init["width"]= 13.07496
mkn[[trait]]$fitRs$chains[[5]]$modelParam$init["width"]=194.29939
mkn[[trait]]$fitRs$chains[[6]]$modelParam$init["width"]=330.13185
mkn[[trait]]$fitRs$chains[[7]]$modelParam$init["width"]=334.20075
mkn[[trait]]$fitRs$chains[[8]]$modelParam$init["width"]=396.23805
mkn[[trait]]$fitRs$chains[[9]]$modelParam$init["width"]=500.32902

## 10^runif(9,-1,1) varH for modelI, modelII, modelIII
mkn[[trait]]$fitRs$chains[[1]]$modelParam$init["varH"]=0.7017772
mkn[[trait]]$fitRs$chains[[2]]$modelParam$init["varH"]=0.8209124
mkn[[trait]]$fitRs$chains[[3]]$modelParam$init["varH"]=2.3528927
mkn[[trait]]$fitRs$chains[[4]]$modelParam$init["varH"]=4.2378902
mkn[[trait]]$fitRs$chains[[5]]$modelParam$init["varH"]=0.7912219
mkn[[trait]]$fitRs$chains[[6]]$modelParam$init["varH"]=0.1392787
mkn[[trait]]$fitRs$chains[[7]]$modelParam$init["varH"]=1.7187575
mkn[[trait]]$fitRs$chains[[8]]$modelParam$init["varH"]=0.7357516
mkn[[trait]]$fitRs$chains[[9]]$modelParam$init["varH"]=0.6686415

## runif(6,0,30) muL for modelII, modelIII
mkn[[trait]]$fitRs$chains[[4]]$modelParam$init["muL"]= 6.5385148
mkn[[trait]]$fitRs$chains[[5]]$modelParam$init["muL"]=14.5150005
mkn[[trait]]$fitRs$chains[[6]]$modelParam$init["muL"]=15.4981188
mkn[[trait]]$fitRs$chains[[7]]$modelParam$init["muL"]=23.4146040
mkn[[trait]]$fitRs$chains[[8]]$modelParam$init["muL"]=18.7001206
mkn[[trait]]$fitRs$chains[[9]]$modelParam$init["muL"]= 0.3482463

## runif(6,0,30) muR for modelII, modelIII 
mkn[[trait]]$fitRs$chains[[4]]$modelParam$init["muR"]= 2.633156
mkn[[trait]]$fitRs$chains[[5]]$modelParam$init["muR"]=25.877732
mkn[[trait]]$fitRs$chains[[6]]$modelParam$init["muR"]=18.239498
mkn[[trait]]$fitRs$chains[[7]]$modelParam$init["muR"]=23.209860
mkn[[trait]]$fitRs$chains[[8]]$modelParam$init["muR"]= 6.726390
mkn[[trait]]$fitRs$chains[[9]]$modelParam$init["muR"]=25.188848

## 10^runif(6,-1,1) varL for modelII, modelIII
mkn[[trait]]$fitRs$chains[[4]]$modelParam$init["varL"]=1.843783
mkn[[trait]]$fitRs$chains[[5]]$modelParam$init["varL"]=1.343432
mkn[[trait]]$fitRs$chains[[6]]$modelParam$init["varL"]=5.218808
mkn[[trait]]$fitRs$chains[[7]]$modelParam$init["varL"]=1.073516
mkn[[trait]]$fitRs$chains[[8]]$modelParam$init["varL"]=2.460740
mkn[[trait]]$fitRs$chains[[9]]$modelParam$init["varL"]=0.178316

## 10^runif(6,-1,1) varR for modelII, modelIII 
mkn[[trait]]$fitRs$chains[[4]]$modelParam$init["varR"]=0.1613381
mkn[[trait]]$fitRs$chains[[5]]$modelParam$init["varR"]=0.4282846
mkn[[trait]]$fitRs$chains[[6]]$modelParam$init["varR"]=2.9530271
mkn[[trait]]$fitRs$chains[[7]]$modelParam$init["varR"]=1.4397268
mkn[[trait]]$fitRs$chains[[8]]$modelParam$init["varR"]=0.3672374
mkn[[trait]]$fitRs$chains[[9]]$modelParam$init["varR"]=2.2332237

## runif(3,10,40) deltaL for modelIII
mkn[[trait]]$fitRs$chains[[7]]$modelParam$init["deltaL"]=19.34910
mkn[[trait]]$fitRs$chains[[8]]$modelParam$init["deltaL"]=19.03688
mkn[[trait]]$fitRs$chains[[9]]$modelParam$init["deltaL"]=12.02479

## runif(3,0,1) tauL for modelIII 
mkn[[trait]]$fitRs$chains[[7]]$modelParam$init["tauL"]=0.4973290
mkn[[trait]]$fitRs$chains[[8]]$modelParam$init["tauL"]=0.5396615
mkn[[trait]]$fitRs$chains[[9]]$modelParam$init["tauL"]=0.9699488

## runif(3,10,40) deltaR for modelIII
mkn[[trait]]$fitRs$chains[[7]]$modelParam$init["deltaR"]=24.87899
mkn[[trait]]$fitRs$chains[[8]]$modelParam$init["deltaR"]=33.96892
mkn[[trait]]$fitRs$chains[[9]]$modelParam$init["deltaR"]=14.27087

## runif(3,0,1) tauR for modelIII
mkn[[trait]]$fitRs$chains[[7]]$modelParam$init["tauR"]=0.1108214
mkn[[trait]]$fitRs$chains[[8]]$modelParam$init["tauR"]=0.5286597
mkn[[trait]]$fitRs$chains[[9]]$modelParam$init["tauR"]=0.7289642

## Run a chain of 3 runs for every fit request
mkn[[trait]]$runs$chains <-  hzar.doChain.multi(mkn[[trait]]$fitRs$chains,
                                            doPar=TRUE,
                                            inOrder=FALSE,
                                            count=3)


## Request an additional run for each chain of modelIII.  These
## three runs will be included in the analysis stage.

mkn[[trait]]$fitRs$extra.modelIII <-
  lapply(lapply(mkn[[trait]]$runs$chains[7:9], # Work from modelIII chains
                function(x) x[[3]]),        # Use third run from chain
         hzar.next.fitRequest)             

## Reduce the tune value as the parameters are near convergence.

mkn[[trait]]$fitRs$extra.modelIII[[1]]$modelParam$tune <- 
  as.list(as.numeric(mkn[[trait]]$fitRs$extra.modelIII[[1]]$modelParam$tune)*0.9)

mkn[[trait]]$fitRs$extra.modelIII[[2]]$modelParam$tune <- 
  as.list(as.numeric(mkn[[trait]]$fitRs$extra.modelIII[[1]]$modelParam$tune)*0.9)

mkn[[trait]]$fitRs$extra.modelIII[[3]]$modelParam$tune <- 
  as.list(as.numeric(mkn[[trait]]$fitRs$extra.modelIII[[1]]$modelParam$tune)*0.9)

## Only do a single run.
mkn[[trait]]$runs$extra.modelIII <-
  hzar.doFit.multi(mkn[[trait]]$fitRs$extra.modelIII,
                   doPar=TRUE,
                   inOrder=FALSE)

## Double check, are these three runs convergent?
# summary(do.call(mcmc.list,
#                 lapply(mkn[[trait]]$runs$extra.modelIII,
#                        hzar.mcmc.bindLL )))


## Start aggregation of data for analysis

## Clear out a spot to collect the data for analysis (note that
## there is currently no "null model" to compare against).
mkn[[trait]]$analysis$initDGs <- list()

## Create a model data group (hzar.dataGroup object) for each
## model from the initial runs.
mkn[[trait]]$analysis$initDGs$modelI <-
  hzar.dataGroup.add(mkn[[trait]]$runs$init$modelI)
mkn[[trait]]$analysis$initDGs$modelII <-
  hzar.dataGroup.add(mkn[[trait]]$runs$init$modelII)
mkn[[trait]]$analysis$initDGs$modelIII <-
  hzar.dataGroup.add(mkn[[trait]]$runs$init$modelIII)

## Create a hzar.obsDataGroup object from the three hzar.dataGroup
## just created, copying the naming scheme (modelI, modelII,
## modelIII).
mkn[[trait]]$analysis$oDG <-
  hzar.make.obsDataGroup(mkn[[trait]]$analysis$initDGs)
mkn[[trait]]$analysis$oDG <-
    hzar.copyModelLabels(mkn[[trait]]$analysis$initDGs,
                         mkn[[trait]]$analysis$oDG)

## Convert all 27 runs to hzar.dataGroup objects, adding them to
## the hzar.obsDataGroup object.
mkn[[trait]]$analysis$oDG <-
  hzar.make.obsDataGroup(lapply(mkn[[trait]]$runs$chains,
                                hzar.dataGroup.add),
                         mkn[[trait]]$analysis$oDG);

## Convert the additional 3 runs to hzar.dataGroup objects, adding
## them to the hzar.obsDataGroup object.
mkn[[trait]]$analysis$oDG <-
  hzar.make.obsDataGroup(lapply(mkn[[trait]]$runs$extra.modelIII,
                                hzar.dataGroup.add),
                         mkn[[trait]]$analysis$oDG);

## Check to make sure that there are only three hzar.dataGroup
## objects named modelI, modelII, and modelIII in the
## hzar.obsDataGroup object.
# print(summary(mkn[[trait]]$analysis$oDG$data.groups))


## Compare the 3 cline models to the null model graphically
# hzar.plot.cline(mkn[[trait]]$analysis$oDG, main="Beard Length");

### Do model selection based on the AICc scores
print(mkn[[trait]]$analysis$AICcTable <-
      hzar.AICc.hzar.obsDataGroup(mkn[[trait]]$analysis$oDG));

## Capture out the model with the minimum AICc score

mkn[[trait]]$analysis$model.name <-
  rownames(mkn[[trait]]$analysis$AICcTable
  )[[ which.min(mkn[[trait]]$analysis$AICcTable$AICc )]]

best.model <- mkn[[trait]]$analysis$model.name

## Extract the hzar.dataGroup object for the selected model
mkn[[trait]]$analysis$model.selected <-
  mkn[[trait]]$analysis$oDG$data.groups[[mkn[[trait]]$analysis$model.name]]

## Look at the variation in parameters for the selected model
# print(hzar.getLLCutParam(mkn[[trait]]$analysis$model.selected,
#                          names(mkn[[trait]]$analysis$model.selected$data.param)));

## Capture the maximum likelihood cline for the selected model
max_likelihood_cline <- hzar.get.ML.cline(mkn[[trait]]$analysis$model.selected)
center <- max_likelihood_cline$param.all$center
width  <- max_likelihood_cline$param.all$width

## Obtain the confidence intervals
conf        <- hzar.qScores.dataGroup(mkn[[trait]]$analysis$model.selected)
center.low  <- conf[1,]$center
center.med  <- conf[2,]$center
center.high <- conf[3,]$center
width.low   <- conf[1,]$width
width.med   <- conf[2,]$width
width.high  <- conf[3,]$width

## Plot the maximum likelihood cline for the selected model
# hzar.plot.cline(mkn[[trait]]$analysis$model.selected, main=paste(trait,experiment));

## Save the cline plot in a pdf
pdf(paste('./',experiment, '_', trait, 'p8_1e5_clineV2.pdf', sep=''), width=4, height=4)

## Plot the 95% credible cline region for the selected model
hzar.plot.fzCline(mkn[[trait]]$analysis$model.selected,
                  main=paste(trait,experiment),
                  ylab= 'Width (mm)',
                  xlab= 'Distance from M. vitellinus (km)',
                  ylim= c(0,14),
                  xlim= c(0,100),
                  las = 1);
dev.off()

# make final dataframe
study       <- c(study, experiment)
morpho_trait<- c(morpho_trait, trait)
center_km   <- c(center_km, center)
width_km    <- c(width_km, width)
best_model  <- c(best_model, best.model)
center_low  <- c(center_low, center.low)
center_med  <- c(center_med, center.med)
center_high <- c(center_high, center.high)
width_low   <- c(width_low, width.low)
width_med   <- c(width_med, width.med)
width_high  <- c(width_high, width.high)

## End Quantitative Trait Analysis


# Create Final dataframe
cline_table <- data.frame(study,
                          morpho_trait,
                          center_km,
                          width_km,
                          best_model,
                          center_low,
                          center_med,
                          center_high,
                          width_low,
                          width_med,
                          width_high)

# Save the cline table dataframe into a file
write.table(cline_table,
            file=paste('./', experiment, '_', trait, 'p8_1e5_cline_parametersV2.tsv', sep=''),
            quote=FALSE,
            sep='\t',
            eol='\n',
            row.names = FALSE,
            col.names = colnames(cline_table))


