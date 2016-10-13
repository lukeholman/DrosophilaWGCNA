###############################
# PARTS YOU WILL NEED TO EDIT::
###############################
# First: set R's working directory to the location of this script!
setwd("~/Dropbox/Data and manuscripts - C 2 copy/Innocenti networks/scripts")

# Also, tell R where you put the raw microarray data (I store this on the hard drive rather than Dropbox, as the files are big)
# You can get the raw microarray data at: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17013
# I recommend you download the data from the link, and put it in a directory called "Innocenti microarray data" on the desktop. Remember to change the filepath's slashes etc if you are on Windows.
filepath <- "~/Desktop/Innocenti microarray data"
###############################

# Install these packages from Bioconductor if you haven't already:
# source("http://www.bioconductor.org/biocLite.R")
# biocLite("affy")
# biocLite("affyPLM")
# biocLite("drosophila2cdf")
# biocLite("org.Dm.eg.db")
# biocLite("drosophila2.db")
# biocLite("Rgraphviz")
# biocLite("topGO")
# biocLite("GOSemSim")
# biocLite("GOstats")

library(affy)
library(affyPLM)
library(drosophila2cdf)
library(org.Dm.eg.db)
library(drosophila2.db)
library(WGCNA)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(lme4)
library(topGO)
library(GOSemSim)
library(GOstats)
library(MCMCglmm)
source("functions.R")    # Load all the custom functions I wrote
options(stringsAsFactors = FALSE) 


# Now we load in the microarray data and normalise it. 
# SKIP THIS STEP if you already run it once and made the file with the normalised expression data - it takes ages, and my function will write the normalised data to disk to save doing it again.
source("Normalise_expression_data.R")
make.normalised.expression.data(filepath = "~/Desktop/Innocenti microarray data", make.sample.cluster.plot = F, write.normalised.data.to.disk = T)

# Load in the expression data, and the data describing each of the 120 samples (their sex, their hemiclone etc)
sampleID <- set.up.data("sampleID") # Now load the identifying data about each sample: its sex, its hemiclone, and the fitness of males and females in the line from which it originates
ned <- as.matrix(read.csv(file.path(filepath, "normalised expression data.csv"), row.names = 1)) # Load in the normalised expression data from the file saved in the above function

# Clean up the expression data: throw away some redundant probes (i.e. duplicated ones, and ones with no Entrez ID)
ned <- discard.redundant.probes(ned) 

# Make a dataset about average gene expression (overall, or by sex, and with and without logs)
# The column "sex.bias" is the log2 ratio of average female over average male expression (so positive numbers mean female biased)
expression <- data.frame(probe = gsub("X", "", colnames(ned)), mean = colMeans(ned), mean.log = colMeans(log(ned)), females = colMeans(ned[sampleID$sex == "female", ]), males = colMeans(ned[sampleID$sex == "male", ]), log.females = colMeans(log(ned)[sampleID$sex == "female", ]), log.males = colMeans(log(ned)[sampleID$sex == "male", ]), IQR = apply(ned,2,IQR), row.names=NULL)
expression$sex.bias <- log2(expression$females/expression$males)
write.csv(expression, file = "../outputs/expression levels.csv", row.names = F) # Write the genes' mean normalised expression levels to disk

ned <- log(ned) # We log-transform the normalised expression data for all subsequent analyses, because it varies over ~10 orders of magnitude

########################################################################################
# First let's make a consensus gene network, to look for genes that are co-expressed in both sexes

# Constructing a weighted gene network entails the choice of the soft thresholding power 'beta' to which co-expression similarity is raised to calculate adjacency (see WGCNA tutorials for info)

nSets <- 2 # We will make a consensus network with two sets: the male and female samples
powerTables <- vector(mode = "list", length = nSets) # Initialize a list to hold the results of scale-free analysis
multiExpr <- set.up.data("multiExpr", ned = ned) # Here is the expression data in the required format for building the consensus network (both sexes together)
powers <- 1:20 # Here are the powers we will assess with the 'pickSoftThreshold' function

# Call the network topology analysis function for each of the two sets in turn - this takes a while!
for (set in 1:nSets) powerTables[[set]] <- list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers, verbose = 2)[[2]])
collectGarbage()

# Make a nice plot to help choose the correct power. 
scale.free.topology.plot()

# After inspecting the figure, we chose the power 18, which is the lowest power for which the model assuming scale-free topology had a fit > 0.90 for both sexes
consensus.power <- 18

# Now let's make the consensus modules, using all 120 samples. The function makes an adjacency matrix, which is then converted to TOM dissimilarity. Note that the use a SIGNED network, meaning that the consensus network analysis looks for sets of transcripts which all go up or down together, in BOTH SEXES. 

both.net <- blockwiseConsensusModules(multiExpr, 
                                 networkType = "signed", 
                                 minModuleSize = 30,
                                 maxBlockSize = 9000, power = consensus.power, 
                                 deepSplit = 2, 
                                 pamRespectsDendro = FALSE, 
                                 mergeCutHeight = 0.25, numericLabels = TRUE,
                                 minKMEtoStay = 0,
                                 saveTOMs = FALSE, verbose = 5)

########################################################################################
# Now let's make two gene networks, one from the 60 female samples and one from the 60 male samples

# Choose a set of soft-thresholding powers to test
powers <- 1:16

# As before, male a nice plot to help choose the correct power. 
# power.picker.plot(ned[sampleID$sex == "female",], powers)
# power.picker.plot(ned[sampleID$sex == "male",], powers)

# These are the powers we get using the plot (those that give >0.9 fit for the scale-free topology model)
female.power <- 15
male.power <- 8

# Make the male and female networks (same parameters for all the network-building functions)
female.net <- blockwiseModules(ned[sampleID$sex == "female",], power = female.power,
                              networkType = "signed", 
                              minModuleSize = 30,
                              maxBlockSize = 9000,
                              deepSplit = 2, 
                              pamRespectsDendro = FALSE, 
                              mergeCutHeight = 0.25, numericLabels = TRUE, 
                              saveTOMs = FALSE, verbose = 5)

male.net <- blockwiseModules(ned[sampleID$sex == "male",], power = male.power,
                            networkType = "signed", 
                            minModuleSize = 30,
                            maxBlockSize = 9000,
                            deepSplit = 2, 
                            pamRespectsDendro = FALSE, 
                            mergeCutHeight = 0.25, numericLabels = TRUE, 
                            saveTOMs = FALSE, verbose = 5)


########################################################################################
# Now let's calculate connectedness for each gene, using the 3 different networks

colnames(ned) <- gsub("X", "", colnames(ned)) # Remove X's from the Affymetrix probe names
rownames(ned) <- gsub("X", "", rownames(ned))
nGenes <- dim(ned)[2]

# Calculate connectedness for each gene using the male and female data together
both.sexes.adjacency <- bigcor(ned, nblocks = 2, power = consensus.power) 
both.kTotal <- numeric(nGenes)
for(i in 1:nGenes) both.kTotal[i] <- sum(both.sexes.adjacency[,i]) - 1 
write.csv(both.kTotal, file = "../outputs/both.kTotal.csv", row.names = F)

rm(both.sexes.adjacency); rm(both.kTotal);  gc()

# Now calculate connectedness for each gene using only the female samples (quite a different set of associations)
female.adjacency <- bigcor(ned[sampleID$sex == "female",], nblocks = 2, power = female.power)
female.kTotal <- numeric(nGenes)
for(i in 1:nGenes) female.kTotal[i] <- sum(female.adjacency[,i])
write.csv(female.kTotal, file = "../outputs/female.kTotal.csv", row.names = F)
rm(female.adjacency);  rm(female.kTotal);  gc()

# Now calculate connectedness for each gene using only the male samples (quite a different set of associations again)
male.adjacency <- bigcor(ned[sampleID$sex == "male",], nblocks = 2, power = male.power)
male.kTotal <- numeric(nGenes)
for(i in 1:nGenes) male.kTotal[i] <- sum(male.adjacency[,i])
write.csv(male.kTotal, file = "../outputs/male.kTotal.csv", row.names = F)
rm(male.adjacency); rm(male.kTotal);  gc()


########################################################################################
# Now let's get the module eigengenes out of the network objects, and make some hemiclone-wide averages for each eigengene
# This is needed for a lot of the other analyses

# Here are the module eigengenes for the consensus network on all 120 samples
ME.consensus.females <- both.net$multiMEs[[1]][[1]] 
ME.consensus.males <- both.net$multiMEs[[2]][[1]]
ME.consensus.females <- ME.consensus.females[, match(paste("ME", 0:max(both.net$colors), sep=""), colnames(ME.consensus.females))] # Re-order the columns so that biggest modules come first
ME.consensus.males <- ME.consensus.males[, match(paste("ME", 0:max(both.net$colors), sep=""), colnames(ME.consensus.males))] # Re-order the columns so that biggest modules come first
MEconsensus.individual.data <- data.frame(rbind(ME.consensus.females, ME.consensus.males), sampleID[order(sampleID$sex),], row.names = NULL) # Link the module eigengene data to the sample ID data
# Take the clone-wide averages for the consensus, and link them to the clone data
MEclones.consensus.females <- apply(ME.consensus.females, 2, function(x) tapply(x, sampleID[sampleID$sex=="female", "hemiclone"], mean))
MEclones.consensus.males <- apply(ME.consensus.males, 2, function(x) tapply(x, sampleID[sampleID$sex=="male", "hemiclone"], mean))
colnames(MEclones.consensus.females) <- paste(colnames(MEclones.consensus.females), ".female", sep="")
colnames(MEclones.consensus.males) <- paste(colnames(MEclones.consensus.males), ".male", sep="")
MEconsensus.clone.data <- data.frame(clone = as.factor(c(rownames(MEclones.consensus.females))),
                          cbind(MEclones.consensus.females, MEclones.consensus.males), 
                          female.fitness = as.numeric(tapply(sampleID$female.fitness, sampleID$hemiclone, mean)), 
                          male.fitness = as.numeric(tapply(sampleID$male.fitness, sampleID$hemiclone, mean)),
                          row.names = NULL)
rm(list = c("ME.consensus.females", "ME.consensus.males", "MEclones.consensus.females", "MEclones.consensus.males"))

# And the module eigengenes for the female-only and male-only networks (60 samples each)
MEfemales.individual.data <- female.net$MEs
MEfemales.individual.data <- MEfemales.individual.data[, match(paste("ME", 0:max(female.net$colors), sep=""), colnames(MEfemales.individual.data))] # Re-order the columns so that biggest modules come first
MEmales.individual.data <- male.net$MEs
MEmales.individual.data <- MEmales.individual.data[, match(paste("ME", 0:max(male.net$colors), sep=""), colnames(MEmales.individual.data))] # Re-order the columns so that biggest modules come first

MEclones.females <- apply(MEfemales.individual.data, 2, function(x) tapply(x, sampleID[sampleID$sex=="female", "hemiclone"], mean))
MEclones.males <- apply(MEmales.individual.data, 2, function(x) tapply(x, sampleID[sampleID$sex=="male", "hemiclone"], mean))

MEfemales.individual.data <- data.frame(MEfemales.individual.data, sampleID[sampleID$sex=="female", -1], row.names = NULL)
MEmales.individual.data <- data.frame(MEmales.individual.data, sampleID[sampleID$sex=="male", -1], row.names = NULL)

MEfemale.clone.data <- data.frame(clone = as.factor(c(rownames(MEclones.females))),
                                  MEclones.females, 
                                  female.fitness = as.numeric(tapply(sampleID$female.fitness, sampleID$hemiclone, mean)), 
                                  male.fitness = as.numeric(tapply(sampleID$male.fitness, sampleID$hemiclone, mean)),
                                  row.names = NULL)
MEmale.clone.data <- data.frame(clone = as.factor(c(rownames(MEclones.males))),
                                  MEclones.males, 
                                  female.fitness = as.numeric(tapply(sampleID$female.fitness, sampleID$hemiclone, mean)), 
                                  male.fitness = as.numeric(tapply(sampleID$male.fitness, sampleID$hemiclone, mean)),
                                  row.names = NULL)
rm(list = c("MEclones.females", "MEclones.males"))

# We now have some nice clear datasets, giving eigengenes for individuals or clones, along with info about the fitness of each clone (and sex and microarry replicate, for the individual-level data)
str(MEconsensus.individual.data)
str(MEfemales.individual.data)
str(MEmales.individual.data)
str(MEconsensus.clone.data)
str(MEfemale.clone.data)
str(MEmale.clone.data)

########################################################################################
# Now let's estimate heritability for each module eigengene. Broad sense heritability is estimated as the proportion of variance explained by hemiclone, i.e. not including residual variance and variance due to batch effects ("replicate")
# This takes a long time, since we are running tons of MCMCglmm models

# For the consensus network, we have data from male and female samples. This means we can also measure the genetic correlation between male and female expression (measured from hemiclone-level variances and covariances)
mcmc.both.heritability <- run.MCMCglmm.heritability(MEconsensus.individual.data, single.sex=F)
mcmc.female.heritability <- run.MCMCglmm.heritability(MEfemales.individual.data, single.sex=T)
mcmc.male.heritability <- run.MCMCglmm.heritability(MEmales.individual.data, single.sex=T)
save(mcmc.both.heritability, file = "../outputs/mcmc.both.heritability.Rdata")
save(mcmc.female.heritability, file = "../outputs/mcmc.female.heritability.Rdata")
save(mcmc.male.heritability, file = "../outputs/mcmc.male.heritability.Rdata")

########################################################################################

# Make a table with information on each module, e.g. selection and genetic architecture (See 'functions' script)
module.data.consensus <- data.frame(structure.of.selection(MEconsensus.clone.data, round=F, consensus = T), 
                                    process.mcmc.models(mcmc.both.heritability)[,c(1,3,5,7,9,10)], 
                                    module.size = as.numeric(table(both.net$colors)), 
                                    mean.log.expression = tapply(read.csv("../outputs/expression levels.csv")$mean.log, both.net$colors, mean),
                                    mean.sex.bias = tapply(read.csv("../outputs/expression levels.csv")$sex.bias, both.net$colors, mean))

module.data.female <- data.frame(structure.of.selection(MEfemale.clone.data, round=F, consensus = F), 
                                    process.mcmc.models(mcmc.female.heritability, single.sex = T), 
                                    module.size = as.numeric(table(female.net$colors)), 
                                    mean.log.expression = tapply(read.csv("../outputs/expression levels.csv")$mean.log, female.net$colors, mean),
                                    mean.sex.bias = tapply(read.csv("../outputs/expression levels.csv")$sex.bias, female.net$colors, mean))

module.data.male <- data.frame(structure.of.selection(MEmale.clone.data, round=F, consensus = F), 
                                 process.mcmc.models(mcmc.male.heritability, single.sex = T), 
                                 module.size = as.numeric(table(male.net$colors)), 
                                 mean.log.expression = tapply(read.csv("../outputs/expression levels.csv")$mean.log, male.net$colors, mean),
                                 mean.sex.bias = tapply(read.csv("../outputs/expression levels.csv")$sex.bias, male.net$colors, mean))


########################################################################################
# Now make a rich dataset with information about every gene, and another one with information about every module

gene.data <- data.frame(affy.name = colnames(ned), consensus.module = both.net$colors, female.module = female.net$colors, male.module = male.net$colors)
chr.mappings <- melt(unlist(as.list(drosophila2CHR)))
gene.data$chromosome <- chr.mappings$value[match(gene.data$affy.name, rownames(chr.mappings))]; rm(chr.mappings)
gene.data$sex.linked <- NA
gene.data$sex.linked[gene.data$chromosome %in% c("2L", "2R", "3L", "3R", "4")] <- "Autosome"
gene.data$sex.linked[gene.data$chromosome  == "X"] <- "X"    # There are only 10 Y-linked markers, and all are in module 0, so let's ignore it
module.data.consensus$percent.sex.linked <- (t(as.matrix(table(gene.data$sex.linked, gene.data$consensus.module))) / rowSums(t(as.matrix(table(gene.data$sex.linked, gene.data$consensus.module)))))[,2]
module.data.female$percent.sex.linked <- (t(as.matrix(table(gene.data$sex.linked, gene.data$female.module))) / rowSums(t(as.matrix(table(gene.data$sex.linked, gene.data$female.module)))))[,2]
module.data.male$percent.sex.linked <- (t(as.matrix(table(gene.data$sex.linked, gene.data$male.module))) / rowSums(t(as.matrix(table(gene.data$sex.linked, gene.data$male.module)))))[,2]
gene.Entrez <- melt(unlist(as.list(drosophila2ENTREZID)))
gene.data$EntrezID <- gene.Entrez$value[match(gene.data$affy.name, rownames(gene.Entrez))]; rm(gene.Entrez)  # Add the Entrez IDs
gene.flybase <- melt(unlist(as.list(drosophila2FLYBASE)))
gene.data$flybase <- gene.flybase$value[match(gene.data$affy.name, rownames(gene.flybase))]; rm(gene.flybase)  # Add the Flybase IDs
gene.name <- melt(unlist(as.list(drosophila2GENENAME)))
gene.data$name <- gene.name$value[match(gene.data$affy.name, rownames(gene.name))]; rm(gene.name)
map.positions <- melt(unlist(lapply(as.list(drosophila2CHRLOC), mean))) # The start point for each gene on the focal chromo. Averaged, in the case of genes that have multiple start points. Chromosomal locations on the antisense strand have a leading "-" sign (e. g. -1234567)
gene.data$map.position <- map.positions$value[match(gene.data$affy.name, rownames(map.positions))]; rm(map.positions)

# Associate the gene info with the genes' effects on fitness, as recorded in Table S1 from Innocenti and Morrow:
# I don't actually use this info in the paper, but it lets you compare my estimates of SA selection with theirs
gene.data$female.fitness.associated <- 0
gene.data$male.fitness.associated <- 0
gene.data$overall.fitness.associated <- 0
gene.data$sexually.antagonistic <- 0
gene.data$female.fitness.associated[gene.data$affy.name %in% read.csv("../data/transcripts associated with female fitness.csv")$Probe] <- 1
gene.data$male.fitness.associated[gene.data$affy.name %in% read.csv("../data/transcripts associated with male fitness.csv")$Probe] <- 1
gene.data$overall.fitness.associated[gene.data$affy.name %in% read.csv("../data/transcripts associated with both sexes fitness.csv")$Probe] <- 1
gene.data$sexually.antagonistic[gene.data$affy.name %in% read.csv("../data/sexually antagonistic transcripts.csv")$Probe] <- 1
module.data.consensus$prop.SA <- as.numeric(tapply(gene.data$sexually.antagonistic, gene.data$consensus.module, sum) / tapply(gene.data$sexually.antagonistic, gene.data$consensus.module, length))
module.data.female$prop.SA <- as.numeric(tapply(gene.data$sexually.antagonistic, gene.data$female.module, sum) / tapply(gene.data$sexually.antagonistic, gene.data$female.module, length))
module.data.male$prop.SA <- as.numeric(tapply(gene.data$sexually.antagonistic, gene.data$male.module, sum) / tapply(gene.data$sexually.antagonistic, gene.data$male.module, length))

# Compute the SA score for every gene
gene.data <- data.frame(gene.data, compute.SA.score.per.gene(ned, sampleID))

# Associate the gene info with the genes' mean expression levels and variability across samples as measured by IQR:
expression <- read.csv("../outputs/expression levels.csv")
gene.data <- data.frame(gene.data, expression[,names(expression) %in% c("mean.log", "sex.bias")])
gene.data$IQR <- apply(ned,2,IQR)
# Associate the gene info with the genes' connectivities:
gene.data$both.kTotal <- read.csv("../outputs/both.kTotal.csv")[,1]
gene.data$female.kTotal <- read.csv("../outputs/female.kTotal.csv")[,1]
gene.data$male.kTotal <- read.csv("../outputs/male.kTotal.csv")[,1]

# Associate the gene info with the dN/dS data, which was calculated from a comparison of D. mel and D. simulans by Chuanzhu Fan in this paper: http://www.ncbi.nlm.nih.gov/pubmed/23314322
# Again I don't actually report this data (the paper was long enough already!) But you can see that I looked for dN/dS related results.
dNdS <- lapply(list.files("../data/dNdS data from Chuanzhu Fan/", full.names = T), clean.up.fan.files)
names(dNdS) <- sapply(strsplit(list.files("../data/dNdS data from Chuanzhu Fan/"), split="_"), function(x) x[2])
dNdS <- melt(dNdS)
dNdS$dNdS <- as.numeric(dNdS$dNdS)
gene.flybase <- melt(unlist(as.list(drosophila2FLYBASE)))
gene.data$dNdS <- dNdS$dNdS[match(gene.data$affy.name, rownames(gene.flybase))]
# sort(unique(gene.data$dNdS )) # There are a small number of genes with suspiciously massive dN/dS - let's delete them (cases where dS was 0, so dN/dS was arbitrarily defined as 99 by Fan?)
gene.data$dNdS[gene.data$dNdS > 90] <- NA
# First make a list giving the number of GO terms associated with every probe that was used in our network analysis 
x <- drosophila2GO
mapped_genes <- mappedkeys(x)
x <- sapply(as.list(x[mapped_genes]), length)
x <- x[which(names(x) %in%  colnames(ned))] # List of probes with >0 GO terms that are in the set of probes I used
gene.data$GOterms <- x[match(gene.data$affy.name, names(x))]
gene.data$GOterms[is.na(gene.data$GOterms)] <- 0
rm(list = c("x", "mapped_genes"))

# Calculate the GO terms that are enriched for each module (See 'functions' script)
GO.terms.per.module.consensus <- module.GOstats(gene.data, "consensus.module")
GO.terms.per.module.female <- module.GOstats(gene.data, "female.module")
GO.terms.per.module.male <- module.GOstats(gene.data, "male.module")
write.csv(GO.terms.per.module, file = "../outputs/GO.terms.per.module.csv")
# Add the number of enriched GO terms to the module data, and also the median number of GO terms that are listed per gene in that module
module.data.consensus$Enriched.GO.terms <- c(NA, as.numeric(sapply(GO.terms.per.module.consensus, nrow)))
module.data.female$Enriched.GO.terms <- c(NA, as.numeric(sapply(GO.terms.per.module.female, nrow)))
module.data.male$Enriched.GO.terms <- c(NA, as.numeric(sapply(GO.terms.per.module.male, nrow)))
module.data.consensus$Median.GO.terms.per.gene <- as.numeric(with(gene.data, tapply(GOterms, consensus.module, median)))
module.data.female$Median.GO.terms.per.gene <- as.numeric(with(gene.data, tapply(GOterms, female.module, median)))
module.data.male$Median.GO.terms.per.gene <- as.numeric(with(gene.data, tapply(GOterms, male.module, median)))

# Add the tissue specificity data from Flyatlas
flyatlas <- read.delim("http://flyatlas.org/20090519all.txt", sep = "\t")
probe <- flyatlas$Oligo
flyatlas <- with(flyatlas, data.frame(brain = Brain.fly,
                                      head = Head.fly,
                                      crop = ratio,
                                      midgut = Midgut.fly,
                                      hindgut = Hindgut.fly,
                                      tubule = tubule.fly,
                                      ovary = Ovary.fly,
                                      testis = Testis.fly,
                                      acc = Acc.fly, # male accessory glands
                                      fatbody = ratio.3,
                                      tag = tag.fly, # Thoracicoabdominal ganglion
                                      carcass = car.fly,
                                      saliv = ratio.1,
                                      virgin.sp = SptV.ratio,
                                      mated.sp = SptM.ratio,
                                      eye = eye.ratio,
                                      heart = heart.ratio,
                                      trachea = trachea.ratio  ))
flyatlas[flyatlas < 2] <- 0 # Define genes as enriched or not using a log2 threshold
flyatlas[flyatlas >= 2] <- 1
gene.data <- data.frame(gene.data, flyatlas[match(gene.data$affy.name, probe),])
rm(list = c("probe", "flyatlas"))

# Save all the data needed for the analyses as a list stored in an R object
all.data <- list(both.net=both.net, female.net=female.net, male.net=male.net, 
                 GO.terms.per.module.consensus=GO.terms.per.module.consensus, 
                 GO.terms.per.module.female=GO.terms.per.module.female,
                 GO.terms.per.module.male=GO.terms.per.module.male,
                 module.data.consensus=module.data.consensus,module.data.female=module.data.female,module.data.male=module.data.male,
                 MEconsensus.clone.data=MEconsensus.clone.data, MEfemale.clone.data=MEfemale.clone.data, MEmale.clone.data=MEmale.clone.data,
                 MEconsensus.individual.data=MEconsensus.individual.data,MEfemales.individual.data=MEfemales.individual.data,MEmales.individual.data=MEmales.individual.data,
                 gene.data=gene.data
                 )

save(all.data, file = "../outputs/allData.Rdata")




