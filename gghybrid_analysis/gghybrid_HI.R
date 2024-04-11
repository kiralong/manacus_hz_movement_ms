# Script for getting genomic hybrid index using R package gghybrid

#Original code written by Richard Ian Bailey 17 April 2018#
#Adapted code viewed here written by Kira Long and Angel Rivera-Colon

#Install the package from GitHub#
#install.packages("devtools")
#library(devtools)
#devtools::install_github("ribailey/gghybrid")

#Attach it#
library(gghybrid)
library(data.table)

#Take a look at the function help files#
#?read.data	#Read in a data file in structure format or similar#
#?data.prep	#Prepare the data for analysis#
#?esth		#Run hybrid index estimation#
#?plot_h	#Plot estimated hybrid indices and credible intervals#
#?ggcline	#Run genomic cline estimation#
#?plot_clinecurve	#Plot one or more fitted cline curves, with or without individual data points#
#?compare.models	#compare two models run on the same data set using the widely applicable information criterion#

#(Note: gghybrid relies on the data.table package for data manipulation)#

#Set working directory
#setwd("/path/to/gghybrid/Long_samples_maf01")

#Read in the data file (The formatted structure file output from stacks)#

dat=read.data("gghybrid_structure.csv",
              nprecol=2,MISSINGVAL=NA)

#Data preparation and filtering. Here I'm filtering out loci that have a minor allele
#frequency greater than 0.1 in both parental reference sets. There are also options for
#filtering by difference in parental allele frequencies, and for number of allele copies
#in each parental reference set (variable among loci due to missing data).

#The function uses objects produced by 'read.data'#

prepdata=data.prep(data=dat$data,
                   loci=dat$loci,
                   alleles=dat$alleles,
                   #S0=c("020SS"), #POPID names for the first parental reference set#
                   S0="P1", #POPID for parent 1 in simulated data
                   #S1="100CG", #POPID names for the second parental reference set#
                   S1="P2", #POPID for parent 2 in simulated data
                   precols=dat$precols,
                   max.S.MAF = 0.01,	#Filtering by parental minor allele frequency# @KML: Original cutoff was at 10%
                   return.genotype.table=T,
                   return.locus.table=T)

#'return.genotype.table=T' makes an optional table of genotypes, where for each locus
#an individual's genotype (assuming diploidy) will be 0 (two copies of the allele with
#relatively higher frequency in the 'S0' parental set), 1 (heterozygote), 2 (two copies
#of the designated 'S1' allele). This table isn't needed in downstream functions, but
#could be useful e.g. for estimating parental linkage disequilibria (associations of
#alleles from the same parent species).

#'return.locus.table=T' is also optional and not needed downstream. It's a table
#with one row per marker, giving some information on parental allele frequencies, sample
#sizes etc.

#Next, run hybrid index estimation

#This function uses objects produced by both the previous functions

hindlabel= esth(data.prep.object = prepdata$data.prep,
                read.data.precols = dat$precols,
                include.Source = TRUE,	#Set to TRUE if you want hybrid indices for the parental reference individuals
                nitt=10000,burnin=5000) #original
                #nitt=50000,burnin=10000) #high run
                #nitt=1000,burnin=100) #For testing that everything is running

write.table(hindlabel,file = "HI_table.txt",quote = FALSE, row.names = FALSE)
write.table(prepdata$locus.data,"locus_table.txt",quote = FALSE, row.names = FALSE)

quit("yes")


### Below is some optional plotting

#Plot a subset of the estimated hybrid indices (the resulting object 'abc' is useful for making a legend)#
pdf('./hybrid_index.pdf',8,8) #To make a pdf of the graph
#par(mar=c(5,5,4,2)) #To resize the margins of the graph on the pdf

setkey(hindlabel$hi,POPID)	#function from data.table, for rapid sorting and subsetting#

par(mar=c(5,5,4,5)) #To resize the margins of the graph for hybrid index

abc = plot_h(data=hindlabel$hi[c("020SS","030SO","040FC","050RO","060QP","080RU","090PR","100CG")],#Subset of POPIDs#
             test.subject=hindlabel$test.subject,
             mean.h.by="POPID",			#Calculate the mean hybrid index for each value of the "POPID" column#
             sort.by=c("mean_h","POPID","h_posterior_mode"),	#Order test subjects along the x axis by the mean hybrid
             #index calculated above and also by individual hybrid index#
             col.group="POPID",
             group.sep="POPID",
             fill.source=TRUE,
             basic.lines=FALSE,
             source.col=c("dodgerblue","firebrick"),
             source.limits=c("dodgerblue","firebrick"),
             cex=1,pch=16,
             cex.lab=1.5,cex.main=1.5,ylim=c(0,1),
             las=1)

#Reshape the plot window as you want#

#Add a legend using the 'plot_h' object 'abc'#

setkey(abc,rn)		#Order data by row number#
legend("bottomright",	#Place the legend in the top left of the figure# @KML: Overwrote legend position to outside the plot
       c("2","3","4","5","6","8","9","10"),   # @KML: Override the legend names to have the right IDs for CG and SS. Careful as it overrides the data itself.
       # abc[,POPID], 		#Name of the field by which data point colours are grouped# (@KML: These are the original legend labels)
       bg="white",			#Background colour#
       text.col=c("black"), #Text colour#
       pch=22, 				#Text size#
       col=abc[,col.Dark2], #Name of the field containing colour information#
       pt.bg=abc[,col.Dark2],	#Name of the field containing colour information#
       ncol=1,				#Number of columns for the group names#
       cex=1, pt.cex=1,
       xpd = TRUE,
       title = "Population")
dev.off() #to turn off the pdf
###


quit("yes")

