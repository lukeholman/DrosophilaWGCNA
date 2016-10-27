# Discards probes with no Entrez ID. Also pick the most variable probe and discards the others, for Entrez IDs with multiple probes
discard.redundant.probes <- function(ned){
  entrez.id <- unlist(as.list(drosophila2ENTREZID)) # The IDs for each probe
  probes.before <- ncol(ned)
  ned <- ned[,which(!is.na(entrez.id))] # Discard all the probes that don't have an Entrez ID
  probes.after.entrez <- ncol(ned)
  entrez.id <- entrez.id[which(!is.na(entrez.id))]
  probes <-colnames(ned)
  IQR <- apply(ned,2,IQR) # The IQR for each of the remaining probes

  # Find duplicate probes that map to the same EntrezID (i.e. gene), and discard all but the most variable of the duplicated probes for each Entrez ID
  # This procedure copies the recommendations in the GOstats R package documentation
  probes.duplicated.Entrez <- names(table(entrez.id)[table(entrez.id) > 1])

  # probes.duplicated.Entrez <- probes.duplicated.Entrez[probes.duplicated.Entrez != "40940"] # Note: we explicitly keep the two probes for doublesex, to look for alternative splicing (one probe is in the female isoform, one in the male isoform)

  to.discard <- matrix("", ncol=100,nrow=length(probes.duplicated.Entrez))

  for(i in 1:length(probes.duplicated.Entrez)) {
    focal.probes <- probes[entrez.id %in% probes.duplicated.Entrez[i]]
    focal.IQRs <- IQR[entrez.id %in% probes.duplicated.Entrez[i]]
    focal.probes <- names(which(focal.IQRs != max(focal.IQRs)))
    to.discard[i, 1:length(focal.probes)] <- focal.probes
  }  
  to.discard <- c(to.discard)[c(to.discard) != ""]
  ned <- ned[,!(colnames(ned) %in% to.discard)]
  probes.after.duplicates <- ncol(ned)
  print(paste("Throwing away", probes.before-probes.after.entrez, "probes due missing Entrez IDs, plus a further", probes.after.entrez-probes.after.duplicates, "duplicated probes."))
  return(ned)
}



# Runs linear models on every gene, in order to compute a sexual antagonism score for each one (uses I+M 2010 Evolution formula for index 'I')
# Note that I don't need to correct for 'block', because I am dealing in line mean fitness and line mean transcript abundance. The 4 blocks are perfectly balanced with respect to line and sex.
compute.SA.score.per.gene <- function(ned, sampleID){
  # Find the hemiclone mean expression level for each probe
  female.clone.means <- apply(ned[sampleID$sex == "female", ], 2, function(x) as.numeric(tapply(x, sampleID$hemiclone[sampleID$sex == "female"], mean)))
  male.clone.means <- apply(ned[sampleID$sex == "male", ], 2, function(x) as.numeric(tapply(x, sampleID$hemiclone[sampleID$sex == "male"], mean)))
  
  # Scale these to have mean 0, variance 1
  female.clone.means <- apply(female.clone.means, 2, scale)
  male.clone.means <- apply(male.clone.means, 2, scale)
  
  # Get the mean fitness of each hemiclone
  female.fitness.by.clone <- as.numeric(tapply(sampleID$female.fitness, sampleID$hemiclone, mean))
  male.fitness.by.clone <- as.numeric(tapply(sampleID$male.fitness, sampleID$hemiclone, mean))
  
  # Make sure fitness is expressed as relative fitness, i.e. defined relative to the 15 clones in the study
  female.fitness.by.clone <- female.fitness.by.clone / mean(female.fitness.by.clone)
  male.fitness.by.clone <- male.fitness.by.clone / mean(male.fitness.by.clone)
  
  # Run "nGenes" linear models. Essentially we get the standardised selection gradient on every transcript. Sample size is 15 lines
  female.slopes <- sapply(1:ncol(female.clone.means), function(i) as.numeric(lm(female.fitness.by.clone ~ female.clone.means[,i])$coefficients[2]))
  male.slopes <- sapply(1:ncol(male.clone.means), function(i) as.numeric(lm(male.fitness.by.clone ~ male.clone.means[,i])$coefficients[2]))
  
  # Compute the Innocenti and Morrow's index for each transcript
  SA.score <- numeric(length(female.slopes))
  for(i in 1:length(female.slopes)) {
    SA.score[i] <- (female.slopes[i] * male.slopes[i]) / sqrt((female.slopes[i]^2 + male.slopes[i]^2) / 2)
  }
  
  # For the SA genes, record whether higher expression levels benefit females or males. "Female genes" are ones where females are selected for higher expression, males lower, and vice versa for "male genes"
  gene.gender <- rep(" ", length(SA.score))
  gene.gender[female.slopes > 0 & male.slopes < 0] <- "Female"
  gene.gender[male.slopes > 0 & female.slopes < 0] <- "Male"
  
  return(data.frame(SA.score = SA.score, gene.gender = gene.gender))
}


# Clean up the dN/dS data kindly provided by C. Fan into a more R-like format
# Note that I removed all mention of dN/dS from the paper for brevity
clean.up.fan.files <- function(fan.file){
  dd<-read.delim(fan.file, header=F, sep = "=")
  names(dd) <- c("flybaseID", "t", "S", "N", "dNdS", "dN", "dS")
  dd$flybaseID <- gsub("\tt", "", dd$flybaseID)
  dd <- as.data.frame(apply(dd, c(1,2), function(x) gsub(" ", "", x)))
  dd$t <- gsub("S", "", dd$t)
  dd$S <- gsub("N", "", dd$S)
  dd$N <- gsub("dN/dS", "", dd$N)
  dd$dNdS <- gsub("dN", "", dd$dNdS)
  dd$dN <- gsub("dS", "", dd$dN)
  split.ID <- strsplit(dd$flybaseID, split="_")
  dd$flybaseID <- sapply(split.ID, function(x) x[1])
  return(dd)
}


# Function to set up data in the required formats - this keeps the main scripts neat and readable
# "sampleID" - set up the sample ID data
# "multiExpr" - sets up the multi expression data for the soft thresholding and the consensus module building
set.up.data <- function(type, ned = NA){
  if(type=="sampleID")
  {
    sampleID <- read.csv("../data/SampleIDs.csv")
    sampleID$replicate <- as.character(sampleID$replicate)
    rownames(sampleID) <- sampleID[,1]
    sampleID <- sampleID[, 2:ncol(sampleID)]
    fem <-  melt(with(read.csv("../data/innocenti females.csv"), tapply(female.rel.fitness, H, mean)))  # Merge in the data on the fitness of these hemiclones
    male <- melt(with(read.csv("../data/innocenti males.csv"), tapply(male.rel.fitness, H, mean)))
    line.fitness <- fem
    line.fitness$male.fitness <- male$value[match(fem$Var1, male$Var1)]
    sampleID$female.fitness <- fem$value[match(sampleID$hemiclone, fem$Var1)]
    sampleID$male.fitness <- male$value[match(sampleID$hemiclone, male$Var1)]
    sampleID$type <- "balanced"
    sampleID$type[sampleID$male.fitness - sampleID$female.fitness > 0.17] <- "male beneficial"
    sampleID$type[sampleID$male.fitness - sampleID$female.fitness < -0.27] <- "female beneficial"
    return(sampleID)
  }
  
  if(type=="multiExpr")
  {
    # Form multi-set expression data: columns starting from 9 contain actual expression data. 
    multiExpr <- vector(mode = "list", length = nSets) 
    multiExpr[[1]] <- list(data = as.data.frame(ned[row.names(ned) %in% row.names(sampleID)[sampleID$sex == "female"], ]))
    names(multiExpr[[1]]$data) <- colnames(ned[row.names(ned) %in% row.names(sampleID)[sampleID$sex == "female"], ])
    rownames(multiExpr[[1]]$data) <- rownames(ned[row.names(ned) %in% row.names(sampleID)[sampleID$sex == "female"], ])
    multiExpr[[2]] <- list(data = as.data.frame(ned[row.names(ned) %in% row.names(sampleID)[sampleID$sex == "male"], ]))
    names(multiExpr[[2]]$data) <- colnames(ned[row.names(ned) %in% row.names(sampleID)[sampleID$sex == "male"], ])
    rownames(multiExpr[[2]]$data) <- rownames(ned[row.names(ned) %in% row.names(sampleID)[sampleID$sex == "male"], ])
    # Check that the data has the correct format for many functions operating on multiple sets: 
    print(checkSets(multiExpr))
    return(multiExpr)
  }
}


# Plots to help determine the best power to use to obtain scale-free topology - most code borrowed from WGCNA tutorials
scale.free.topology.plot <- function(){
  colors <- c("red", "blue")
  setLabels <- c("Females", "Males")  # For easier labeling of plots, create a vector holding descriptive names of the two sets. 
  # Will plot these columns of the returned scale free analysis tables
  plotCols <- c(2,5,6,7)
  colNames <- c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity", "Max connectivity")
  # Get the minima and maxima of the plotted points
  ylim <- matrix(NA, nrow = 2, ncol = 4)
  for (set in 1:nSets)
  {
    for (col in 1:length(plotCols))
    {
      ylim[1, col] <- min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE)
      ylim[2, col] <- max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE)
    }
  }
  # Plot the quantities in the chosen columns vs. the soft thresholding power
  sizeGrWindow(8, 6)
  par(mfcol = c(2,2));
  par(mar = c(4.2, 4.2 , 2.2, 0.5))
  cex1 <- 0.7
  for (col in 1:length(plotCols)) for (set in 1:nSets)
  {
    if (set==1)
    {
      plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2], xlab="Soft Threshold (power)", ylab=colNames[col], type="n", ylim = ylim[, col], main = colNames[col])
      addGrid()
    }
    if (col==1) text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2], labels=powers,cex=cex1,col=colors[set])
    else text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]], labels=powers,cex=cex1,col=colors[set])
    if (col==1) legend("bottomright", legend = setLabels, col = colors, pch = 20) 
    else legend("topright", legend = setLabels, col = colors, pch = 20) 
  }
}

# A similar function - only this one is for single-set expression data used to make male and female networks, not multiset data used to make consensus
power.picker.plot <- function(ned, powers){
  # Call the network topology analysis function
  sft = pickSoftThreshold(ned, powerVector = powers, verbose = 5)
  # Plot the results:
  sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.90,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
}


# Function to calculate a big correlation matrix in chunks (when it's too big to do in R's memory)
# Adapted from http://www.r-bloggers.com/bigcor-large-correlation-matrices-in-r/ - I corrected a typo involving MAT, and I added 2 lines (annotated)
# I use this because the built-in WGCNA functions for getting connectivity tend to fail for my dataset, since the correlation matrix is too big.
bigcor <- function(x, nblocks = 10, verbose = TRUE, make.dissimilarity=F, power)
{
  library(ff, quietly = TRUE)
  NCOL <- ncol(x)

  ## test if ncol(x) %% nblocks gives remainder 0
  if (NCOL %% nblocks != 0) stop("Choose different 'nblocks' so that ncol(x) %% nblocks = 0!")

  ## preallocate square matrix of dimension
  ## ncol(x) in 'ff' single format
  corMAT <- ff(vmode = "single", dim = c(NCOL, NCOL))

  ## split column numbers into 'nblocks' groups
  SPLIT <- split(1:NCOL, rep(1:nblocks, each = NCOL/nblocks))

  ## create all unique combinations of blocks
  COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
  COMBS <- t(apply(COMBS, 1, sort))
  COMBS <- unique(COMBS)

  ## iterate through each block combination, calculate correlation matrix
  ## between blocks and store them in the preallocated matrix on both
  ## symmetric sides of the diagonal
  for (i in 1:nrow(COMBS)) {
    COMB <- COMBS[i, ]
    G1 <- SPLIT[[COMB[1]]]
    G2 <- SPLIT[[COMB[2]]]
    if (verbose) cat("Block", COMB[1], "with Block", COMB[2], "\n")
    flush.console()
    COR <- cor(x[, G1], x[, G2])

    # Luke's addition: added abs() and raise by 'power'
    COR <- abs(COR)^power

    # Luke's addition: set this option to true if you'd like to subtract every correlation from 1, to get a dissimilarity matrix instead of a similarity matrix
    if(make.dissimilarity) COR <- 1 - COR

    corMAT[G1, G2] <- COR
    corMAT[G2, G1] <- t(COR)
    COR <- NULL
  }

  gc()
  return(corMAT)
}


# This function calculates the association between mean fitness and mean module eigengenes across hemiclones. Also calculates the I+M2010 index of sex-specific selection I for each module.
structure.of.selection <- function(MEclones, round=T, consensus=T)
{
  # First we express the fitness of each clone relative to the average of the 15 clones in this dataset (they are initially expressed relative to the 100 clones in Innocenti and Morrow)
  MEclones$female.fitness <- MEclones$female.fitness / mean(MEclones$female.fitness)
  MEclones$male.fitness <- MEclones$male.fitness / mean(MEclones$male.fitness)
  
  n.modules <- length(grep("ME", names(MEclones))) # Number of modules (includes module zero)
  if(consensus) n.modules <- n.modules / 2
  
  for(i in grep("ME", names(MEclones))) MEclones[,i] <- scale(MEclones[,i]) # Rescale the eigengenes to have mean 0, var 1 (So the fitness-trait slopes are comparable across eigengenes)
  output <- data.frame(Module = paste("ME", 0:(n.modules-1), sep=""), beta.female = 0, beta.male = 0, average.beta = 0,
                       female.p = 0, male.p = 0, SA.score = 0, advantaged.sex = " ")
  if(consensus) output <- data.frame(output, sex.specific.selection.pval = 0)
  
  module.names <- names(MEclones)[grep("ME", names(MEclones))]
  
  for(i in 1:n.modules){
    if(consensus){
      female.formula <- as.formula(paste("female.fitness ~", module.names[i]))
      male.formula <- as.formula(paste("male.fitness ~", module.names[i+n.modules]))
      beta.female <- as.data.frame(summary(lm(female.formula, data=MEclones))$coefficients) # Regress clone mean female fitness and clone mean female eigengene
      beta.male <- as.data.frame(summary(lm(male.formula, data=MEclones))$coefficients) # Same for males
    }
    else{
      beta.female <- as.data.frame(summary(lm(MEclones$female.fitness ~ MEclones[,names(MEclones) == module.names[i]]))$coefficients)
      beta.male <- as.data.frame(summary(lm(MEclones$male.fitness ~ MEclones[,names(MEclones) == module.names[i]]))$coefficients)
    }
    output$beta.female[i] <- beta.female[2,1]
    output$beta.male[i] <- beta.male[2,1]
    output$female.p[i] <- beta.female[2,4]
    output$male.p[i] <- beta.male[2,4]
    if(consensus) {
      fitness <- c(MEclones$female.fitness, MEclones$female.fitness)
      ME <- c(MEclones[,names(MEclones) == module.names[i]], MEclones[,names(MEclones) == module.names[i+n.modules]])
      sex <- rep(c("female", "male"), each=nrow(MEclones))
      compare.models <- anova(lm(fitness ~ ME * sex), lm(fitness ~ ME + sex))
      output$sex.specific.selection.pval[i] <- compare.models[6][2,] # test if the sex-by-eigengene term is significant using a partial F-test
    }
    
    # Calculate the SA score for the focal module, following Innocenti and Morrow 2010 Evolution 10.1111/j.1558-5646.2010.01021.x
    output$SA.score[i] <- with(output, (beta.female[i] * beta.male[i]) / sqrt((beta.female[i]^2 + beta.male[i]^2) / 2))
  }
  output$average.beta <- with(output, rowMeans(cbind(beta.female, beta.male)))
  
  if(round) output[,-1] <- round(output[,-1], 3)
  return(output)
}



# Run one MCMCglmm model per module, to calculate heritability of each module and the male-female genetic correlation
run.MCMCglmm.heritability <- function(MEdata, exclude.mod0 = F, single.sex = F, nitt = 550000, thin = 450, burnin = 100000){
  require(MCMCglmm)
  if(exclude.mod0) MEdata <- MEdata[, -grep("ME0", names(MEdata))]
  num.modules <- length(grep("ME", names(MEdata)))

  if(!single.sex) # for the concensus module analysis, make two separate traits for each module eigengene (ME) depending on focal indiv's sex
  {
    males <- MEdata[, grep("ME", names(MEdata))]
    names(males) <- paste(names(males), ".male", sep="")
    males[MEdata$sex == "female", ] <- NA
    names(MEdata)[grep("ME", names(MEdata))] <- paste(names(MEdata)[grep("ME", names(MEdata))], ".female", sep="")
    MEdata[MEdata$sex == "male", grep("ME", names(MEdata))] <- NA
    MEdata <- data.frame(males, MEdata)
  }
  
  MEdata[, grep("ME", names(MEdata))] <- apply(MEdata[, grep("ME", names(MEdata))], 2, scale) # Scale the eigengenes to have mean 0, var 1

  models <- vector(mode = "list", length = num.modules) 
  for(i in 1:num.modules){
    print(paste("Doing module", i, "of", num.modules))
    if(!single.sex) # For consensus analysis, fit the value of the eigengene as expressed in males and females as a two-column vector
    {
      prior <- list(G=list(G1=list(V=diag(2)*1e-3, nu=3, alpha.mu=rep(0,2), alpha.V=diag(2)*1000), # The random effect of hemiclone
                           G2=list(V=diag(2)*1e-3, nu=3, alpha.mu=rep(0,2), alpha.V=diag(2)*1000)), # The random effect of replicate (4 batches)
                    R=list(V=diag(2), nu=3))
      focal.formula <- as.formula(paste("cbind(", names(MEdata)[i], ", ", names(MEdata)[i+num.modules], ") ~ 1", sep=""))
      models[[i]] <-  MCMCglmm(focal.formula, 
                               random = ~us(trait):hemiclone + us(trait):replicate, 
                               rcov = ~idh(trait):units, # We can't estimate the residual covariance, since no individual is both male and female
                               prior = prior, family = rep("gaussian", 2), 
                               data = MEdata, nitt = nitt, thin = thin, burnin = burnin, verbose = F)
    }
    else 
    {
      prior <- list(G=list(G1=list(V=1e-3, nu=3, alpha.mu=0, alpha.V=1000),  # The random effect of hemiclone
                           G2=list(V=1e-3, nu=3, alpha.mu=0, alpha.V=1000)), # The random effect of replicate (4 batches)
                    R=list(V=1, nu=3))
      focal.formula <- as.formula(paste(names(MEdata)[i], "~ 1"))
      models[[i]] <-  MCMCglmm(focal.formula, 
                               random = ~ hemiclone + replicate, 
                               prior = prior, family = "gaussian", 
                               data = MEdata, nitt = nitt, thin = thin, burnin = burnin, verbose = F)
    }
  }
  return(models)
}

# Process the results of the MCMCglmm models into useful statistics
process.mcmc.models <- function(mcmc.models, single.sex=F){
  if(!single.sex)
  {
    mcmc.results <- as.data.frame(do.call("rbind", lapply(mcmc.models, function(x) posterior.mode(x$VCV))))
    CIs <- as.data.frame(do.call("rbind", lapply(lapply(mcmc.models, function(x) HPDinterval(x$VCV)[c(1,2,4), ]), function(x) apply(x, 1, function(xx) paste0(round(xx,2), collapse = " - ")))))
    mcmc.results <- data.frame(Male.var = mcmc.results[,1], Male.CI = CIs[,1], Intersex.cov = mcmc.results[,2], Intersex.CI = CIs[,2], Female.var = mcmc.results[,3], Female.CI = CIs[,3])
    mcmc.results$gen.corr <- sapply(mcmc.models, function(x) posterior.mode(x$VCV[,2] / sqrt(x$VCV[,1] * x$VCV[,4])))
    corr.CIs <- do.call("rbind", lapply(mcmc.models, function(x) HPDinterval(x$VCV[,2] / sqrt(x$VCV[,1] * x$VCV[,4]))))
    mcmc.results$corr.CIs <- apply(corr.CIs, 1, function(x) paste0(round(x,2), collapse = " - "))
    mcmc.results$male.heritability <- sapply(mcmc.models, function(x) posterior.mode(x$VCV[,1] / (x$VCV[,1]+x$VCV[,5]+x$VCV[,9])))
    mcmc.results$female.heritability <- sapply(mcmc.models, function(x) posterior.mode(x$VCV[,4] / (x$VCV[,4]+x$VCV[,8]+x$VCV[,10])))
  }
  else
  {
    mcmc.results <- data.frame(var = do.call("rbind", lapply(mcmc.models, function(x) posterior.mode(x$VCV)))[,1])
    mcmc.results$heritability <- sapply(mcmc.models, function(x) posterior.mode(x$VCV[,1] / (x$VCV[,1]+x$VCV[,2]+x$VCV[,3])))
  }
  return(mcmc.results)
}


# Do GO term enrichment on each module using GOstats package. "module.column" tells R where to look for the module info 
module.GOstats <- function(gene.data, module.column){
  
  gene.data$module <- gene.data[,names(gene.data) == module.column]
  
  # This whole analysis only uses probes that clustered into a module
  foc.gene.data <- gene.data[gene.data$module != 0, ]
  modules <- sort(unique(foc.gene.data$module))
  output <- vector('list', length(modules)) # Declare empty list to hold the output
  
  # Let's define the "gene universe" as all probes (in the duplicate-free probe set) that:
  # - are in a module other than 0, 
  # - that have Entrez IDs,
  # - have GO terms associated with them
  have.no.GO <- sapply(as.list(drosophila2GO), function(x) is.na(x)[1])
  have.no.GO <- names(have.no.GO)[have.no.GO]
  universe <- foc.gene.data$EntrezID[!is.na(foc.gene.data$EntrezID) & foc.gene.data$module != 0 & !(foc.gene.data$affy.name %in% have.no.GO)] 

  for(i in 1:length(modules)){   # Now for each of the modules...
    print(paste("Doing module", i, "of", length(modules)))
    # Let's define the gene set that is associated with the focal module. Again, this exclude genes without EntrezIDs, or which have no GO terms
    gene.set <- foc.gene.data$EntrezID[!is.na(foc.gene.data$EntrezID) & foc.gene.data$module == modules[i] & !(foc.gene.data$affy.name %in% have.no.GO)] 
    names(output)[i] <- paste("ME", modules[i], sep="")
    if(length(gene.set) == 0) output[[i]] <- "No GO terms known for this module"
    else
    {
      params <- new("GOHyperGParams",
                    geneIds = gene.set,
                    universeGeneIds = universe,
                    annotation="drosophila2.db",
                    ontology="BP",
                    pvalueCutoff=0.001,
                    conditional=FALSE,
                    testDirection="over")
      output[[i]] <- summary(hyperGTest(params))
    }
  }
  return(output)
}

# Plots the dendrogram of genes clustered by the WGCNA package. Takes a gene network object as its argument. Only plots the first 2 blocks for blockwise analyses
module.cluster.plot <- function(network){
  moduleColors <- labels2colors(c(0,1:max(network$colors)))
  sizeGrWindow(12,6)
  # Use the layout function for more involved screen sectioning
  layout(matrix(1:4, 2), heights = c(0.8, 0.2), widths = c(1,1))
  layout.show(4);
  nBlocks = length(both.net$dendrograms)
  # Plot the consensus network dendrogram and the module colors underneath for each block
  for (block in 1:2) # Note this only shows first 2 blocks
    plotDendroAndColors(network$dendrograms[[block]], moduleColors[1 + network$colors[network$blockGenes[[block]]]],
                        "Module colors",
                        main = paste("Gene dendrogram and module colors in block", block),
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05,
                        setLayout = FALSE)
}


# Function to plot each clone's mean module eigengene against its mean fitness, for males and female separately
fitness.vs.eigengenes.plot <- function(clone.data, module.data, consensus = F){
  
  clone.data <- clone.data[, !(names(clone.data) %in% "ME0")]   # Don't plot module zero
  module.data <- module.data[module.data$Module != "ME0", ]
  
  # First we express the fitness of each clone relative to the average of the 15 clones in this dataset (they are presently expressed relative to the 100 clones in Innocenti and Morrow)
  clone.data$female.fitness <- clone.data$female.fitness / mean(clone.data$female.fitness)
  clone.data$male.fitness <- clone.data$male.fitness / mean(clone.data$male.fitness)
  
  clone.data[, grep("ME", names(clone.data))] <- scale(clone.data[, grep("ME", names(clone.data))]) # SCALE THE MODULES so they have mean 0 and SD 1 (i.e. we're gonna find standardised selection gradients)
  plots <- vector(mode = "list", length = nrow(module.data)) 
  module.name <- module.data$Module
  female.colour <- rgb(228, 0, 27, maxColorValue = 255)
  male.colour <- rgb(3, 95, 254, maxColorValue = 255)
  
  xlim <- range(clone.data[, grep("ME", names(clone.data))])
  ylim <- range(clone.data[, grep("fitness", names(clone.data))])
  xlim <- xlim + c(-0.05*(xlim[2]-xlim[1]), 0.05*(xlim[2]-xlim[1]))
  ylim <- ylim + c(-0.05*(ylim[2]-ylim[1]), 0.05*(ylim[2]-ylim[1]))
  x.adjust <- (xlim[2] - xlim[1]) * 0.02
  y.adjust <- (ylim[2] - ylim[1]) * 0.02
  
  for(i in 1:length(plots))
  {
    if(!consensus) 
    {
      clone.data$focal <- clone.data[, names(clone.data) == module.name[i]]
      plots[[i]] <- ggplot(clone.data, aes(x = focal)) + geom_smooth(method="lm", aes(y = female.fitness), colour = female.colour) + geom_smooth(method="lm", aes(y = male.fitness), colour = male.colour) + geom_point(size=0.8, aes(y = female.fitness), colour = female.colour) + geom_point(size=0.8, aes(y = male.fitness), colour = male.colour)  + theme_bw()  + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), line = element_line(colour = "black"), panel.border = element_rect(colour="black"), plot.title = element_text(size = 8)) +xlab(NULL) + ylab(NULL) + coord_cartesian(xlim = xlim, ylim = ylim) + annotate("text", x = xlim[1]+x.adjust, y = ylim[2]-y.adjust, label = i, fontface="bold") + scale_x_continuous(breaks = c(-3,-1.5,0,1.5,3)) + scale_y_continuous(breaks = c(0.5, 1, 1.5))
    }
    else
    {
      clone.data$focal.f <- clone.data[, names(clone.data) == paste(module.name[i],".female",sep="")]
      clone.data$focal.m <- clone.data[, names(clone.data) == paste(module.name[i],".male",sep="")]
      plots[[i]] <- ggplot(clone.data) + geom_point(size=0.8, aes(x = focal.f, y = female.fitness), colour = female.colour) + geom_smooth(method="lm", aes(x = focal.f, y = female.fitness), colour = female.colour) + geom_smooth(method="lm", aes(x = focal.m, y = male.fitness), colour = male.colour) + geom_point(size=0.8, aes(x = focal.m, y = male.fitness), colour = male.colour)  + theme_bw()  + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), line = element_line(colour = "black"), panel.border = element_rect(colour="black"), plot.title = element_text(size = 10)) + xlab(NULL) + ylab(NULL) + coord_cartesian(xlim = xlim, ylim = ylim) + annotate("text", x = xlim[1]+x.adjust, y = ylim[2]-y.adjust, label = i, fontface="bold")+ scale_x_continuous(breaks = c(-3,-1.5,0,1.5,3)) + scale_y_continuous(breaks = c(0.5, 1, 1.5))
    }
  }
  x <- do.call(arrangeGrob, plots)
  return(grid.arrange(x, bottom = "Hemiclone mean eigengene", left = "Hemiclone mean fitness"))
}

# NOT IMPLEMENTED
# This function calculates the similarity between a specified set of modules in terms of the GO terms of the genes in them
# cluster.modules.by.GO.semantic.similarity <- function(gene.data, module.column, exclude0 = T){
#   module <- gene.data[, names(gene.data) == module.column]
#   
#   if(exclude0) {  # Get rid of module zero, if desired
#     gene.data <- gene.data[module != 0, ]
#     module <- module[module != 0]
#   }
#   
#   combos <- data.frame(t(combn(unique(module),2)))   # These are all the possible pairwise combinations of 2 different modules
#   names(combos) <- c("mod1", "mod2")
#   combos$GOsimilarity <- 0  # Add a column to hold the similarities
#   
#   nModules <- max(combos[,2])
#   output <- matrix(0, ncol=nModules, nrow=nModules)   # Make the output into a n*n matrix of distances, where n is number of modules
# 
#   do.one.combo <- function(mod1=NA, mod2=NA)
#   { 
#     cluster1 <- gene.data$EntrezID[module == mod1]
#     cluster2 <- gene.data$EntrezID[module == mod2]
#     return(clusterSim(cluster1, cluster2, ont="MF", organism="fly", measure="Wang")) # Note that we use MF ontology. The BP one seems to be broken for this function.
#   }
# 
#   for(i in 1:nrow(combos)){
#     print(paste("Doing combination", i, "out of", nrow(combos)))
#     combos$GOsimilarity[i] <- do.one.combo(mod1=combos[i,1], mod2=combos[i,2])
#   }
#   
#   for(i in 1:nrow(combos)) {
#       output[combos[i,1], combos[i,2]] <- combos$GOsimilarity[i]
#       output[combos[i,2], combos[i,1]] <- combos$GOsimilarity[i]
#   }
#   return(output)
# }


# Plots observed-expected figures by tissue and module
tissue.enrichment.plot <- function(module.data, module.column){
  
  focal.module <- gene.data[, names(gene.data) %in% module.column]
  nModules <- max(focal.module)
  tissue <- apply(gene.data[,(which(names(gene.data) == "brain")):(which(names(gene.data) == "trachea"))], 2, function(x) tapply(x, focal.module, sum))
  exp <- t(sapply(rowSums(tissue), function(x) x * colSums(tissue)/sum(colSums(tissue))))
  obs.exp <- melt((tissue - exp)/exp)    # The y-axis shows the number of enriched transcripts - expected, divided by expected.
  #obs.exp <- obs.exp[obs.exp$Var1 != 0,]
  ranks <- rank(-module.data$SA.score)
  ranks <- (1+max(ranks)) - ranks
  ranks <- match(1:length(ranks), ranks)
  ranks[is.na(ranks)] <- (1:nModules)[!((1:nModules) %in% ranks)]
  
  obs.exp$title <- paste("Module", obs.exp$Var1)
  obs.exp$title <- paste(obs.exp$title, ", I: ", format(round(module.data$SA.score[match(obs.exp$title, paste("Module", 0:nModules))],2), nsmall = 2), sep="")
  obs.exp$title <- factor(obs.exp$title, levels = unique(obs.exp$title)[ranks])
  obs.exp$Var2 <- as.character(obs.exp$Var2)
  obs.exp$Var2 <- paste0(toupper(substr(obs.exp$Var2, 1, 1)), substr(obs.exp$Var2, 2, nchar(obs.exp$Var2)))
  obs.exp$Var2[obs.exp$Var2 == "Tag"] <- "TAG"
  obs.exp$Var2[obs.exp$Var2 == "Saliv"] <- "Saliv. gland"
  obs.exp$Var2[obs.exp$Var2 == "Acc"] <- "Acc. gland"
  obs.exp$Var2[obs.exp$Var2 == "Virgin.sp"] <- "Virgin ST"
  obs.exp$Var2[obs.exp$Var2 == "Mated.sp"] <- "Mated ST"
  obs.exp$Var2 <- factor(obs.exp$Var2, levels = c("Head", "Brain", "Eye", "Crop", "Midgut", "Hindgut", "Tubule", 
                                                  "TAG", "Saliv. gland", "Fatbody", "Carcass",  "Testis", "Acc. gland", "Ovary", "Virgin ST", "Mated ST",
                                                  "Heart", "Trachea"))
  names(obs.exp) <- c("Module", "Tissue", "Fold_enrichment", "Panel")
  
  obs.exp$SA <- rep(module.data$SA.score, length(unique(obs.exp$Tissue)))
  
  ten.and.one <- levels(obs.exp$Panel) %in% c("Module 10, SA: 0", "Module 1, SA: 0") # R doesn't sort 10 and 1 very well, so manually correct it if needed
  if(sum(ten.and.one) == 2) 
  {
    new.levels <- levels(obs.exp$Panel)
    new.levels[ten.and.one] <- c("Module 1, SA: 0", "Module 10, SA: 0")
    obs.exp$Panel <- factor(obs.exp$Panel, levels = new.levels)
  }
  
  p1 <- ggplot(obs.exp, aes(x = Module, y= Fold_enrichment, fill=Tissue)) + geom_hline(aes(yintercept= 0))  + geom_hline(aes(yintercept= -0.5),linetype=2, colour = "grey") + geom_hline(aes(yintercept=2),linetype=2, colour = "grey") + geom_bar(stat="identity",position = "dodge") + facet_wrap(~Panel, scales="free_x", ncol=4) + xlab(NULL) +ylab("Tissue enrichment") + theme_bw() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(colour=NULL, fill = NULL), panel.border = element_rect(fill=NA, colour="black",size=1), legend.key = element_blank(), legend.title = element_blank())
  
  min.SA.across.all.networks <- min(c(module.data.consensus$SA.score, module.data.female$SA.score, module.data.male$SA.score))
  max.SA.across.all.networks <- max(c(module.data.consensus$SA.score, module.data.female$SA.score, module.data.male$SA.score))
  # Following code borrowed from Stack Overflow here: http://stackoverflow.com/questions/19440069/ggplot2-facet-wrap-strip-color-based-on-variable-in-data-set/21589891
  dummy <- ggplot(data = obs.exp, aes(x = Module, y = Fold_enrichment))+ facet_wrap(~Panel,ncol=4) + 
    geom_rect(aes(fill=SA), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, colour = "black",size=1) + scale_fill_gradient2(high = "#C7A9FD", mid = "white", low = "#F48B94", limits = c(min.SA.across.all.networks, max.SA.across.all.networks)) +
    theme_minimal()
  
  g1 <- ggplotGrob(p1)
  g2 <- ggplotGrob(dummy)
  
  gtable_select <- function (x, ...) 
  {
    matches <- c(...)
    x$layout <- x$layout[matches, , drop = FALSE]
    x$grobs <- x$grobs[matches]
    x
  }
  
  panels <- grepl(pattern="panel", g2$layout$name)
  strips <- grepl(pattern="strip_t", g2$layout$name)
  g2$layout$t[panels] <- g2$layout$t[panels] - 1
  g2$layout$b[panels] <- g2$layout$b[panels] - 1
  
  new_strips <- gtable_select(g2, panels | strips)
  grid.newpage()
  grid.draw(new_strips)
  
  gtable_stack <- function(g1, g2){
    g1$grobs <- c(g1$grobs, g2$grobs)
    g1$layout <- transform(g1$layout, z= z-max(z), name="g2")
    g1$layout <- rbind(g1$layout, g2$layout)
    g1
  }
  new_plot <- gtable_stack(g1, new_strips)
  grid.newpage()
  grid.draw(new_plot)
}



# Function to make the connectivity plot
# Note: this plot was removed from the finished paper
gene.connectivity.plot <- function(){
  plot.dat <- with(gene.data, data.frame(network = rep(c("Consensus network", "Female network", "Male network"), each = nrow(gene.data)),
                                         k =  c(both.kTotal / max(both.kTotal),
                                                female.kTotal / max(female.kTotal),
                                                male.kTotal / max(male.kTotal)),
                                         sex.bias = rep(sex.bias,3),
                                         SA.score = rep(SA.score,3)
  ))
  plot.dat$bias <- "Unbiased"
  plot.dat$bias[plot.dat$sex.bias > 1] <- "Female" # sex-biased genes are those with at least 2-fold difference in expression
  plot.dat$bias[plot.dat$sex.bias < -1] <- "Male"
  plot.dat$bias <- factor(plot.dat$bias, levels = c("Female", "Unbiased", "Male"))
  plot.dat$SA <- "Unselected"
  plot.dat$SA[plot.dat$SA.score <= -0.05] <- "SA"  # SA genes are those with a score of at least -0.05
  plot.dat$SA[plot.dat$SA.score >= 0.05] <- "Concordant"
  plot.dat$SA <- factor(plot.dat$SA, levels = c("SA", "Unselected", "Concordant"))
  #dat <- melt(tapply(plot.dat$k, list(plot.dat$network, plot.dat$bias, plot.dat$SA), median))
  print(melt(tapply(plot.dat$k, list(plot.dat$SA, plot.dat$bias, plot.dat$network), max)))
  dat <- melt(tapply(plot.dat$k, list(plot.dat$SA, plot.dat$bias, plot.dat$network), median)) %>%
    mutate(lower = melt(tapply(plot.dat$k, list(plot.dat$SA, plot.dat$bias, plot.dat$network), function(x) quantile(x, prob = 0.025)))[,4]) %>% mutate(higher = melt(tapply(plot.dat$k, list(plot.dat$SA, plot.dat$bias, plot.dat$network), function(x) quantile(x, prob = 0.975)))[,4])

  names(dat)[1:4] <- c("SA", "bias", "network", "k")

  ggplot(dat, aes(x = bias, y = k, group=SA)) + geom_errorbar(position=position_dodge(0.5), aes(ymin=lower, ymax=higher), width = 0)  + geom_point(position=position_dodge(0.5), size=4, aes(fill = SA), shape=21)  + theme_bw() + facet_wrap(~network)+ theme(legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour="black",size=1), strip.background = element_rect(fill="white", color="black",size=1), legend.key = element_blank()) + ylab("Relative connectivity \n\u00B1 95% quantiles") + xlab("Sex bias in expression") + scale_y_continuous(limits=c(0,1)) + scale_fill_manual(values=c("purple", "white", "yellow"))
}


# Function by Hadley Wickham that makes multiple ggplots that share a legend (used for the connectedness plot). Source: https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {

  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)

  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid.newpage()
  grid.draw(combined)
}


# Function to make plot showing that male-benefitting and female-benefitting transcripts are clustered, and these modules tend to be SA (Figure 3)
modularity.of.SA.plot <- function(){
  
  individual.pair.of.plots <- function(gene.data, module.data, mod.var, title, xlab = c(" ", " "), ylab = c(" ", " ")){
    gene.data <- gene.data[mod.var != 0, ] # Omit module zero
    mod.var <- droplevels(mod.var[mod.var != 0]) # Omit module zero
    module.data <- module.data[module.data$Module != "ME0", ] # Omit module zero
    dat <- melt(table(gene.data$gene.gender, mod.var))
    print(table(gene.data$gene.gender, mod.var)[2:3, ])
    print(chisq.test(table(gene.data$gene.gender, mod.var)[2:3, ]))
    dat[,2] <- factor(dat[,2], levels = levels(mod.var))
    dat <- dat[!is.na(dat[,2]), ]
    dat$value <- dat$value / rep(as.numeric(tapply(dat$value, dat$mod.var, sum)), each=3)
    dat <- dat[rev(order(dat$Var1)), ]
    dat2 <- data.frame(skew = tapply(dat$value, dat$mod.var, function(x) 2*(-0.5 + x[x == max(x[1:2])][1] / sum(x[1:2]))), 
                       SA.score = sort(module.data$SA.score))
    pink.col <- rgb(233, 197, 203, maxColorValue = 255) # rgb(237, 186, 198, maxColorValue = 255)
    blue.col <-  rgb(201, 237, 254, maxColorValue = 255)
    p1 <- ggplot(dat, aes(x = mod.var, y = 100*value, fill = Var1)) + geom_bar(stat="identity", colour = "black") + scale_fill_manual(values = c("white", pink.col,  blue.col)) + theme_bw() + scale_x_discrete(expand =c(0,0)) + scale_y_continuous(expand =c(0,0)) + xlab(NULL) + ylab(NULL) + theme(legend.position = "none", panel.border = element_blank())
    
    p2 <- ggplot(dat2, aes(x = skew, y = SA.score)) + stat_smooth(method="lm", colour = "black") + geom_hline(yintercept = 0, linetype = 2, colour = "grey") + geom_point(alpha=0.8)  +coord_cartesian(xlim =c(0, 1), ylim = c(-.18, .12)) + theme_bw() +theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour="black",size=1)) + xlab(NULL) + ylab(NULL) 
    return(list(p1,p2))
  } # end
  
  p1 <- individual.pair.of.plots(gene.data, module.data.consensus, gene.data$modules.in.order.of.SA.con, "Consensus network")
  p2 <- individual.pair.of.plots(gene.data, module.data.female, gene.data$modules.in.order.of.SA.fem, "Female network", ylab = c("% genes","Sexual antagonism score"))
  p3 <- individual.pair.of.plots(gene.data, module.data.male, gene.data$modules.in.order.of.SA.male, "Male network", xlab = c("Module", "Skew towards one sex"))
  
  return(grid.arrange(
    arrangeGrob(arrangeGrob(p1[[1]],p2[[1]],p3[[1]],ncol=1), bottom = "Module (from most to least sexually antagonistic)", left = "% genes benefitting males or females"), 
    arrangeGrob(arrangeGrob(p1[[2]],p2[[2]],p3[[2]],ncol=1), bottom = "Skew of SA genes towards one sex", left="Sex-specific selection index (I)"), 
    ncol=2, widths=c(1, 0.5)
  ))
}

# Function to make a neat LaTex table of the GO enrichment results
make.latex.table <- function(GO.data, suffix, nTerms){
  GO.data <- lapply(GO.data, function(x) {
    max <- nrow(x)
    if(max > nTerms) max <- nTerms
    x[1:max,]
  })
  lengths <- sapply(GO.data, nrow)
  nModules <- length(GO.data)
  GO.data <- do.call("rbind", GO.data)
  GO.data <- data.frame(Module = paste("M", unlist(mapply(rep,1:nModules, each=lengths)), suffix, sep = ""), 
                        GO.ID=GO.data$GOBPID, 
                        Term = GO.data$Term, 
                        OddsRatio = round(GO.data$OddsRatio,1), stringsAsFactors = F)
  GO.data$Term <-   paste0(toupper(substr(GO.data$Term, 1, 1)), substr(GO.data$Term, 2, nchar(GO.data$Term)))
  GO.data <- GO.data[apply(GO.data, 1, function(x) sum(is.na(x)))==0,]
  modules <- unique(GO.data$Module)
  last.row <- numeric(length(modules) - 1)
  for(i in 1:length(modules) - 1) last.row[i] <- tail(which(GO.data$Module == modules[i]), 1)
  
  names(GO.data)[2:4] <- c("GO ID", "GO Term", "Odds Ratio")
  print(xtable(GO.data), include.rownames = F, tabular.environment = "longtable", hline.after = c(-1,0, last.row), floating = F)
}


# Little function to make a table giving all the GO terms for the genes in 'supp.table', i.e. the top 1% most antagonistic ones
get.GO.for.SA.genes <- function(supp.table){
  x <- drosophila2GO
  mapped_genes <- mappedkeys(x)
  xx <- as.list(x[mapped_genes])
  xx <- melt(xx[names(xx) %in% supp.table[,1]])
  xx$value <- as.character(xx$value)
  
  SA.GOs <- vector(mode = "list", length = length(unique(xx$L1)))
  names(SA.GOs) <- unique(xx$L1)
  
  for(i in 1:length(unique(xx$L1))){
    foc <- xx[xx$L1 == unique(xx$L1)[i], ]
    GO.terms <- unique(foc$L2)
    keep <- foc$value[seq(3,nrow(foc),by=3)] == "BP"
    if(length(GO.terms[keep]) > 0) SA.GOs[[i]] <- GO.terms[keep]
  }
  
  output <- unlist(Term(GOTERM))
  output <- melt(lapply(SA.GOs, function(x) as.character(output[names(output) %in% x])))
  output$name <- gene.data$name[match(output[,2], gene.data$affy.name)]
  output$flybase <- gene.data$flybase[match(output[,2], gene.data$affy.name)]
  names(output) <- c("GO term", "Probe", "Gene name", "FlybaseID")
  return(output)
}


# Simulation to test if genes from the same module tend to be non-randomly positioned across the genome (answer: not really)
distance.simulation <- function(module.type, boots=1000){
  chrs <- c("2L", "2R", "3L", "3R", "X")
  gene.data <- gene.data[gene.data[, names(gene.data) %in% module.type] != 0 & gene.data$chromosome %in% chrs, ] # restrict to genes in a module which are on the main chromosomes, not chr 4 or Y
  gene.data$map.position <- abs(gene.data$map.position) # Get rid of the minus signs. This is because genes with positions e.g. -31850551 and 31819663 are right next to each other, just on opposite DNA strands (antisense and sense)
  module <- gene.data[, names(gene.data) %in% module.type]
  unique.modules <- sort(unique(module))
  n.modules <- length(unique.modules)
  chr.lengths <- as.numeric(with(gene.data, tapply(map.position, chromosome, function(x) max(x) - min(x)))) # distance between furthest genes on each chromosome
  
  output <- data.frame(expand.grid(unique.modules, chrs), mean.difference = 0, lowerCI=0, upperCI = 0)
  names(output)[1:2] <- c("Module", "Chromosome")
  output$Module <- factor(output$Module)
  counter <- 1
  
  # The test statistic is (median - mean) of the nearest neighbour distances - i.e. skew.  The distribution will be more right skewed under clustering, more left skewed under uniform distribution
  # We divide the skew statistic by the length of the chromosome, such that its maximum conceivable range is bounded by 1 and -1. Makes the results from different chromosomes more comparable
  test.statistic <- function(positions, chromosome) {
    x <- sort(positions) # sort the positions from smallest to largest
    n <- length(x)
    upward <- x[3:n] - x[2:(n-1)] # compare item 2 with item 3, item 3 with item 4...
    downward <- x[2:(n-1)] - x[1:(n-2)] # compare item 2 with item 1, item 3 with item 2... the nearest neighbour is the minimum value for upward[i] + downward[i]
    nearest.neighbour.dists <- c(x[2]-x[1], apply(rbind(upward, downward), 2, min), x[n]-x[n-1]) # We already know the distance for the first and last item, so append them here
    (mean(nearest.neighbour.dists) - median(nearest.neighbour.dists)) / chr.lengths[chrs == chromosome] # mean-median / chr length
  }
  
  for(i in 1:n.modules){
    for(j in 1:length(chrs)){
      print(paste("Doing combo ", counter, " of ", nrow(output), ".", sep=""))
      focal <- gene.data[module == unique.modules[i] & gene.data$chromosome == chrs[j], ]  # Focal dataset: just genes from this module, on this chromosome
      local.genes <- gene.data[gene.data$chromosome == chrs[j], ] # Local set of genes for comparison: just the genes on this chromosome (from any module, except the unassigned module 0 genes)
      module.size <- nrow(focal) # Number of genes in the module
      n.genes <- nrow(local.genes) # Number of genes on this chromosome
      
      observed <- test.statistic(focal$map.position, chromosome = chrs[j])
      
      # Now simulate the expected distribution of the test statistic, for this module size on this chromosome
      # Essentially we pick a random sample of 'module.size' genes (where each gene can only appear once per sample), 'boots' times. We will then calculate the test statistic on each random sample of genes
      rand <- matrix(ceiling(runif(module.size * boots) * n.genes), ncol = boots, nrow = module.size) # pick loads of random modules
      for(k in 1:ncol(rand)){
        duplicates <- which(duplicated(rand[,k])) # within each sample, check if there are any duplicated modules
        if(length(duplicates) > 0){ # if yes...
          rand[duplicates, k] <- sample((1:n.genes)[-duplicates], length(duplicates), replace = F) # draw some replacement values from the genes that haven't been used yet
        }
      }
      boot.data <- matrix(local.genes$map.position[c(rand)], ncol = boots, nrow = module.size)
      
      diff.with.null.simulation <-  apply(boot.data, 2, test.statistic, chromosome = chrs[j]) - observed # Subtract the observed value from the simualted ones. Positive differences = more clustered than expected
      output[output$Module == unique.modules[i] & output$Chromosome == chrs[j], 3:5] <- c(mean(diff.with.null.simulation), quantile(diff.with.null.simulation, probs = c(0.025, 0.975))) # Save the mean and 95% quantiles
      counter <- counter + 1
    }
  }
  return(output)
}


# Very similar simulation function that does a distance simulation on the SA and non-SA genes. By default, genes in the top 1% by SA.score are considered SA.
distance.simulation.SAgenes <- function(percentile = 0.01, boots=1000){
  chrs <- c("2L", "2R", "3L", "3R", "X")
  
  gene.data$map.position <- abs(gene.data$map.position) # Get rid of the minus signs. This is because genes with positions e.g. -31850551 and 31819663 are right next to each other, just on opposite DNA strands (antisense and sense)
  
  # restrict to proper chromosomes
  gene.data <- gene.data[gene.data$chromosome %in% chrs, ]
  gene.data$SA <- 0 # define which genes are SA
  gene.data$SA[gene.data$SA.score <= quantile(gene.data$SA.score, probs = percentile)] <- 1
  
  chr.lengths <- as.numeric(with(gene.data, tapply(map.position, chromosome, function(x) max(x) - min(x)))) # distance between furthest genes on each chromosome
  
  output <- data.frame(Chromosome = chrs, mean.difference = 0, lowerCI=0, upperCI = 0)
  counter <- 1
  
  # The test statistic is (median - mean) of the nearest neighbour distances - i.e. skew.  The distribution will be more right skewed under clustering, more left skewed under uniform distribution
  # We divide the skew statistic by the length of the chromosome, such that its maximum conceivable range is bounded by 1 and -1. Makes the results from different chromosomes more comparable
  test.statistic <- function(positions, chromosome) {
    x <- sort(positions) # sort the positions from smallest to largest
    n <- length(x)
    upward <- x[3:n] - x[2:(n-1)] # compare item 2 with item 3, item 3 with item 4...
    downward <- x[2:(n-1)] - x[1:(n-2)] # compare item 2 with item 1, item 3 with item 2... the nearest neighbour is the minimum value for upward[i] + downward[i]
    nearest.neighbour.dists <- c(x[2]-x[1], apply(rbind(upward, downward), 2, min), x[n]-x[n-1]) # We already know the distance for the first and last item, so append them here
    (mean(nearest.neighbour.dists) - median(nearest.neighbour.dists)) / chr.lengths[chrs == chromosome] # mean-median / chr length
  }
  
  for(j in 1:length(chrs)){
    print(paste("Doing combo ", counter, " of ", nrow(output), ".", sep=""))
    focal <- gene.data[gene.data$SA == 1 & gene.data$chromosome == chrs[j], ]  # Focal dataset: just the SA genes on this chromosome
    local.genes <- gene.data[gene.data$chromosome == chrs[j], ] # Local set of genes for comparison: all the genes on this chromosome (SA or not)
    n.SA <- nrow(focal) # Number of SA genes on focal chromosome
    n.genes <- nrow(local.genes) # Number of genes on this chromosome
    
    observed <- test.statistic(focal$map.position, chromosome = chrs[j])
    
    # Now simulate the expected distribution of the test statistic, for this chromosome
    # Essentially we pick a random sample of 'n.SA' genes (where each gene can only appear once per sample), 'boots' times. We will then calculate the test statistic on each random sample of genes
    rand <- matrix(ceiling(runif(n.SA * boots) * n.genes), ncol = boots, nrow = n.SA) # pick loads of random modules
    for(k in 1:ncol(rand)){
      duplicates <- which(duplicated(rand[,k])) # within each sample, check if there are any duplicated modules
      if(length(duplicates) > 0){ # if yes...
        rand[duplicates, k] <- sample((1:n.genes)[-duplicates], length(duplicates), replace = F) # draw some replacement values from the genes that haven't been used yet
      }
    }
    boot.data <- matrix(local.genes$map.position[c(rand)], ncol = boots, nrow = n.SA)
    
    diff.with.null.simulation <-  apply(boot.data, 2, test.statistic, chromosome = chrs[j]) - observed # Subtract the observed value from the simualted ones. Positive differences = more clustered than expected
    output[output$Chromosome == chrs[j], 2:4] <- c(mean(diff.with.null.simulation), quantile(diff.with.null.simulation, probs = c(0.025, 0.975))) # Save the mean and 95% quantiles
    counter <- counter + 1
  }
  return(output)
}


# Run GO term enrichment test. The gene universe is all probes that have Entrez IDs, and which have GO terms associated with them
# The test set of probes is either manually entered, or is equal to the set of probes with a given SA.score
# The arg restrict.to.sig.SA can be used to restrict the test to just genes that were found to be significantly SA by Innocenti and Morrow
# The arg SA.or.concordant tells the function to look below or above the SA.cutoff
GO.enrichment.test <- function(SA.cutoff=NULL, SA.or.concordant=NULL, manual.gene.set = NULL, restrict.to.sig.SA = F){
  
  # Let's define the gene universe as all probes that have Entrez IDs, and have GO terms associated with them
  have.no.GO <- sapply(as.list(drosophila2GO), function(x) is.na(x)[1])
  have.no.GO <- names(have.no.GO)[have.no.GO]
  universe <- gene.data$EntrezID[!is.na(gene.data$EntrezID) & !(gene.data$affy.name %in% have.no.GO)] 
  
  # Let's define the focal gene set. Again, this exclude genes without EntrezIDs, or which have no GO terms
  if(is.null(manual.gene.set) & !restrict.to.sig.SA & SA.or.concordant == "SA") gene.set <- gene.data$EntrezID[!is.na(gene.data$EntrezID) & !(gene.data$affy.name %in% have.no.GO) & gene.data$SA.score <= SA.cutoff] 
  else if(is.null(manual.gene.set) & restrict.to.sig.SA & SA.or.concordant == "SA") gene.set <- gene.data$EntrezID[!is.na(gene.data$EntrezID) & !(gene.data$affy.name %in% have.no.GO) & gene.data$SA.score <= SA.cutoff & gene.data$sexually.antagonistic==1]  
  if(is.null(manual.gene.set) & !restrict.to.sig.SA & SA.or.concordant == "concordant") gene.set <- gene.data$EntrezID[!is.na(gene.data$EntrezID) & !(gene.data$affy.name %in% have.no.GO) & gene.data$SA.score >= SA.cutoff] 
  else if(is.null(manual.gene.set) & restrict.to.sig.SA & SA.or.concordant == "concordant") gene.set <- gene.data$EntrezID[!is.na(gene.data$EntrezID) & !(gene.data$affy.name %in% have.no.GO) & gene.data$SA.score >= SA.cutoff & gene.data$sexually.antagonistic==1]  
  else if(!is.null(manual.gene.set) & !restrict.to.sig.SA) gene.set <- gene.data$EntrezID[gene.data$affy.name %in% manual.gene.set]  
  else if(!is.null(manual.gene.set) & restrict.to.sig.SA) gene.set <- gene.data$EntrezID[gene.data$affy.name %in% manual.gene.set & gene.data$sexually.antagonistic==1]  
  
  params <- new("GOHyperGParams",
                geneIds = gene.set,
                universeGeneIds = universe,
                annotation = "drosophila2.db",
                ontology = "BP",
                pvalueCutoff = 0.001,
                conditional = FALSE,
                testDirection = "over")
  summary(hyperGTest(params))
}




# Takes a list of focal genes, and plots their tissue enrichment scores. 
gene.tissue.plot <- function(focal.genes){
  
  in.or.not <- rep(0, nrow(gene.data))
  in.or.not[gene.data$affy.name %in% focal.genes] <- 1
  
  tissue <- apply(gene.data[,(which(names(gene.data) == "brain")):(which(names(gene.data) == "trachea"))], 2, function(x) tapply(x, in.or.not, sum))
  print(tissue)
  exp <- t(sapply(rowSums(tissue), function(x) x * colSums(tissue)/sum(colSums(tissue))))
  obs.exp <- melt((tissue - exp)/exp)    # The y-axis shows the number of enriched transcripts - expected, divided by expected.
  obs.exp <- obs.exp[obs.exp$Var1 != 0,2:3]
  
  obs.exp$Var2 <- as.character(obs.exp$Var2)
  obs.exp$Var2 <- paste0(toupper(substr(obs.exp$Var2, 1, 1)), substr(obs.exp$Var2, 2, nchar(obs.exp$Var2)))
  obs.exp$Var2[obs.exp$Var2 == "Tag"] <- "TAG"
  obs.exp$Var2[obs.exp$Var2 == "Saliv"] <- "Saliv. gland"
  obs.exp$Var2[obs.exp$Var2 == "Acc"] <- "Acc. gland"
  obs.exp$Var2[obs.exp$Var2 == "Virgin.sp"] <- "Virgin ST"
  obs.exp$Var2[obs.exp$Var2 == "Mated.sp"] <- "Mated ST"
  obs.exp$Var2 <- factor(obs.exp$Var2, levels = c("Head", "Brain", "Eye", "Crop", "Midgut", "Hindgut", "Tubule", 
                                                  "TAG", "Saliv. gland", "Fatbody", "Carcass",  "Testis", "Acc. gland", "Ovary", "Virgin ST", "Mated ST",
                                                  "Heart", "Trachea"))
  names(obs.exp) <- c("Tissue", "Fold_enrichment")
  
  ggplot(obs.exp, aes(x = 1, y= Fold_enrichment, fill=Tissue)) + geom_hline(aes(yintercept= 0))  + geom_hline(aes(yintercept= -0.5),linetype=2, colour = "grey")  + geom_bar(stat="identity",position = "dodge") + xlab(NULL) +ylab("Tissue enrichment") + theme_bw() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(colour=NULL, fill = NULL), panel.border = element_rect(fill=NA, colour="black",size=1), legend.key = element_blank(), legend.title = element_blank())
}
