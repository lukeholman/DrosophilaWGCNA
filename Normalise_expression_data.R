make.normalised.expression.data <- function(filepath, make.sample.cluster.plot = F, write.normalised.data.to.disk=T){
  
  original.filepath <- getwd()
  
  # Set to the directory where the microarray CEX files are located
  setwd(filepath)
  
  # read the data
  affydata <- ReadAffy()
  affydata
  
  # Normalise the microarray data
  # Using mas5 normalisation, which is supposedly better than RMA for making gene networks (See this paper:  http://bioinformatics.oxfordjournals.org/content/23/13/i282.full )
  nvals <- mas5(affydata)
  
  # Here are the normalised expression data
  ned <- exprs(nvals)
  rm(nvals)
  colnames(ned) <- gsub(".CEL.gz", "", colnames(ned)) # Remove the file extensions
  ned <- t(ned)
  
  print("Checking that all samples and genes are good, as classified by the WGCNA package")
  gsg = goodSamplesGenes(ned, verbose = 3) 
  print(gsg$allOK) # Yes they are!
  rm(gsg)
  
  # Write the data to the hard drive, in the directory in "filepath"
  if(write.normalised.data.to.disk) write.csv(ned, file = paste(filepath, "normalised expression data.csv", sep="/")
  
  if(make.sample.cluster.plot){
    # Now let's cluster the SAMPLES, to see if any samples are obvious outliers:
    sampleTree = hclust(dist(t(ned)), method = "average")
    # Plot the sample tree: Open a graphic output window of size 12 by 9 inches
    # The user should change the dimensions if the window is too large or too small.
    sizeGrWindow(12,9)
    #pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
    par(cex = 0.6);
    par(mar = c(0,4,2,0))
    plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
  }
  
  # Restore the working directory to the original one
  setwd(original.filepath) 
}