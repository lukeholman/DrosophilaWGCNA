# Load some libraries. I think but the last 2 can be installed from CRAN, otherwise try Bioconductor
library(ggplot2)
library(gridExtra)
library(WGCNA)
library(reshape2)
library(gtable)
library(grid)
library(car)
library(xtable)
library(dplyr)
library(drosophila2.db) 
library(GOstats)

# Set the working directory - change this as needed
setwd("~/Dropbox/Data and manuscripts - C 2 copy/Innocenti networks/scripts")

################################################
# FIRST, INPUT AND CLEAN UP THE DATA
################################################

load("../outputs/allData.Rdata")   # This is an R object with all the information needed to do the analyses and plot the figures
source("functions.R")
filepath <- "~/Desktop/Innocenti microarray data"
both.net <- all.data$both.net
female.net <- all.data$female.net
male.net <- all.data$male.net
GO.terms.per.module.consensus <- all.data$GO.terms.per.module.consensus
GO.terms.per.module.female <- all.data$GO.terms.per.module.female
GO.terms.per.module.male <- all.data$GO.terms.per.module.male
module.data.consensus <- all.data$module.data.consensus
module.data.female <- all.data$module.data.female
module.data.male <- all.data$module.data.male
MEconsensus.clone.data <- all.data$MEconsensus.clone.data
MEfemale.clone.data <- all.data$MEfemale.clone.data
MEmale.clone.data <- all.data$MEmale.clone.data
MEconsensus.individual.data <- all.data$MEconsensus.individual.data
MEfemales.individual.data <- all.data$MEfemales.individual.data
MEmales.individual.data <- all.data$MEmales.individual.data
gene.data <- all.data$gene.data
rm(all.data)
gene.data$sexually.antagonistic <- as.factor(gene.data$sexually.antagonistic)
gene.data$gene.gender <- factor(gene.data$gene.gender, levels = c("Female","Male", " ")) 
gene.data$modules.in.order.of.SA.con <- factor(gene.data$consensus.module, levels = (0:nrow(module.data.consensus))[rev(order(module.data.consensus$SA.score))]) 
gene.data$modules.in.order.of.SA.fem <- factor(gene.data$female.module, levels = (0:nrow(module.data.female))[rev(order(module.data.female$SA.score))]) 
gene.data$modules.in.order.of.SA.male <- factor(gene.data$male.module, levels = (0:nrow(module.data.male))[rev(order(module.data.male$SA.score))]) 

# Make a new dataframe to facilitate plotting:
plot.data <-rbind(module.data.consensus[-1,c("beta.female", "beta.male","mean.sex.bias", "module.size", "SA.score", "advantaged.sex")], module.data.female[-1,c("beta.female", "beta.male","mean.sex.bias", "module.size", "SA.score", "advantaged.sex")], module.data.male[-1,c("beta.female", "beta.male","mean.sex.bias", "module.size", "SA.score", "advantaged.sex")])
plot.data$network <- as.character(do.call("c", mapply(rep, c("Consensus network", "Female network", "Male network"), each = c(nrow(module.data.consensus[-1,]), nrow(module.data.female[-1,]), nrow(module.data.male[-1,])))))
plot.data$Female.herit <- c(module.data.consensus[-1,"female.heritability"], module.data.female[-1, "heritability"], rep(NA, nrow(module.data.male[-1,])))
plot.data$Male.herit <- c(module.data.consensus[-1,"male.heritability"], rep(NA, nrow(module.data.female[-1,])), module.data.male[-1, "heritability"])
plot.data$Herit.to.plot <- plot.data$Male.herit
plot.data$Herit.to.plot[is.na(plot.data$Herit.to.plot)] <- plot.data$Female.herit[is.na(plot.data$Herit.to.plot)]
plot.data$advantaged.sex <- factor(plot.data$advantaged.sex, levels = c(" ", "Female", "Male"))
plot.data$gen.corr <- c(module.data.consensus$gen.corr[-1], rep(NA, nrow(module.data.female) + nrow(module.data.male) - 2))



################################################
# SOME SIMPLE STATISTICS ABOUT THE MODULES
################################################

# number of probes in the microarray after trimming
length(both.net$colors)

# number of probes in modules (note that 'both' = the consensus network throughout)
length(both.net$colors[both.net$colors !=0])
length(female.net$colors[female.net$colors !=0])
length(male.net$colors[male.net$colors !=0])

# number of modules (excludes the unassigned genes, hence the minus 1)
length(unique(both.net$colors)) - 1
length(unique(female.net$colors)) - 1
length(unique(male.net$colors)) - 1

# Number of genes in each module, for the 3 networks 
table(both.net$colors)
table(female.net$colors)
table(male.net$colors)

# percent unassigned genes in each network
as.numeric(table(both.net$colors)[1]) / sum(table(both.net$colors)) * 100
as.numeric(table(female.net$colors)[1]) / sum(table(female.net$colors)) * 100
as.numeric(table(male.net$colors)[1]) / sum(table(male.net$colors)) * 100

# Gene dendrograms showing modules in each of the 3 networks - not in the paper
module.cluster.plot(both.net)
module.cluster.plot(female.net)
module.cluster.plot(male.net)

# This is not in the paper, but here is a histogram showing that the list of antagonistic genes identified by I+M are often not antagonistic, if one uses the 'I' measure of antagonism, and conducts the selection analysis on hemiclones means as I did:
# Remember, positive scores suggest sexually concordant selection, and scores closes to zero suggest sex-limited selection.
hist(gene.data$SA.score[gene.data$sexually.antagonistic == 1], xlab = "Index of sex-specific selection (I)") 
sum(gene.data$SA.score[gene.data$sexually.antagonistic == 1] > -0.05) / length(gene.data$SA.score[gene.data$sexually.antagonistic == 1])
sum(gene.data$SA.score[gene.data$sexually.antagonistic == 1] > 0) / length(gene.data$SA.score[gene.data$sexually.antagonistic == 1])
sum(gene.data$SA.score[gene.data$sexually.antagonistic == 1] > 0.05) / length(gene.data$SA.score[gene.data$sexually.antagonistic == 1])


#############################################################
# TABLE 1, PLUS ALL THE STATS GIVEN "LOOSE" IN THE TEXT
#############################################################

# Stats asking "Are genes from the same module randomly distributed across chromosomes?"

# First we find the proportion of genes that are on each chromosome
genes.per.chromosome <- with(gene.data[gene.data$chromosome %in% c("2L", "2R", "3L", "3R", "X"), ], table(chromosome))
genes.per.chromosome <- genes.per.chromosome / sum(genes.per.chromosome) # convert to frequencies

# Consensus module Chi-square test:
test.data <- gene.data[gene.data$chromosome %in% c("2L", "2R", "3L", "3R", "X") & gene.data$consensus.module != 0, ]
observed <- with(test.data, table(chromosome, consensus.module))
expected <- do.call("cbind", lapply(colSums(observed), function(x) x * genes.per.chromosome)) # expected number is number of genes in the module*freq they should appear on each chromosome under the null
chi <- sum((observed - expected)^2 / expected)
df <- (ncol(observed) - 1) * (nrow(observed) - 1) # The degrees of freedom for this test
paste("Chi = ", round(chi,2), ", df = ", df, ", p=", round(1-pchisq(chi,df),4), ".", sep="")

# Female module Chi-square test: (module 3 is strongly biased towards genes on the X - this is an ovary enriched module that is SA)
test.data <- gene.data[gene.data$chromosome %in% c("2L", "2R", "3L", "3R", "X") & gene.data$female.module != 0, ]
observed <- with(test.data, table(chromosome, female.module))
expected <- do.call("cbind", lapply(colSums(observed), function(x) x * genes.per.chromosome))
chi <- sum((observed - expected)^2 / expected)
paste("Chi = ", round(chi,2), ", df = ", df, ", p=", round(1-pchisq(chi,df),4), ".", sep="")

# Male module Chi-square test: (module 2 is strongly enriched for genes on chr 2L - this is a testis enriched module that is not SA)
test.data <- gene.data[gene.data$chromosome %in% c("2L", "2R", "3L", "3R", "X") & gene.data$male.module != 0, ]
observed <- with(test.data, table(chromosome, male.module))
expected <- do.call("cbind", lapply(colSums(observed), function(x) x * genes.per.chromosome))
chi <- sum((observed - expected)^2 / expected)
paste("Chi = ", round(chi,2), ", df = ", df, ", p=", round(1-pchisq(chi,df),4), ".", sep="")
rm(list=c("chi", "df", "test.data", "observed", "expected"))



# Table 1: The top 1% most SA genes are enriched for GO terms relating to alternative splicing
go.output.SA <- GO.enrichment.test(SA.cutoff = quantile(gene.data$SA.score, probs = 0.01), SA.or.concordant = "SA")
print(xtable(go.output.SA[-(1:2)][,c(5,1:4)]), include.rownames=F) # print table in LaTex format

# Similar table not in the text: The top 1% most sexually concordant genes are enriched for GO terms relating to programmed cell death
go.output.concordant <- GO.enrichment.test(SA.cutoff = quantile(gene.data$SA.score, probs = 0.99), SA.or.concordant = "concordant")
print(xtable(go.output.concordant[-(1:2)][,c(5,1:4)]), include.rownames=F) # print table in LaTex format


# Plot showing that genes in the top 1% most SA list are under-represented in the fly atlas list of testes genes (discussed but not shown in paper, since it just replicates conclusions already made by I+M)
gene.tissue.plot(gene.data$affy.name[gene.data$SA.score <= quantile(gene.data$SA.score, probs = 0.01)])

# Plot showing that genes in the top 1% most concordant list are under-represented in the fly atlas list of spermatheca genes (discussed but not shown in paper, since it just replicates conclusions already made by I+M)
gene.tissue.plot(gene.data$affy.name[gene.data$SA.score >= quantile(gene.data$SA.score, probs = 0.99)])

# Have a look at the data for some famous genes involved in the sex determination cascade - some are sexually antagonistic, some are not.
# Ovo = CG6824 gene product from transcript CG6824-RE
gene.data[gene.data$name %in% c("Sex lethal", "doublesex", "transformer", "transformer 2", "sans fille", "CG6824 gene product from transcript CG6824-RE", "ovarian tumor", "virilizer", "female lethal d", "fused"), c(1,9,15,16:18)]

# SA genes are not clustered on particular chromosomes
observed <- t(rbind(table(gene.data$chromosome[gene.data$chromosome %in% c("2L", "2R", "3L", "3R",  "X") & gene.data$SA.score <= quantile(gene.data$SA.score, probs = 0.01)]), table(gene.data$chromosome[gene.data$chromosome %in% c("2L", "2R", "3L", "3R",  "X") & gene.data$SA.score > quantile(gene.data$SA.score, probs = 0.01)])))
expected <- do.call("cbind", lapply(colSums(observed), function(x) x * genes.per.chromosome))
chi <- sum((observed - expected)^2 / expected)
df <- (ncol(expected)-1)*(nrow(expected)-1)
paste("Chi = ", round(chi,2), ", df = ", df, ", p=", round(1-pchisq(chi,df),4), ".", sep="")

# SA genes are not clustered within chromosomes
distance.simulation.SAgenes(boots = 10000, percentile = 0.99)



################################################
# THE MAIN PAPER'S FIGURES, AND ASSOCIATED STATS
################################################

# FIGURE 2: some modules are sexually antagonisitic, especially in the male and female networks
pdf("../figures/modules.are.SA.pdf", width = 8.5, height = 3.7)
ggplot(plot.data, aes(x = beta.female, y = beta.male)) + geom_hline(yintercept=0,linetype=2,colour="grey") + geom_vline(xintercept=0,linetype=2,colour="grey") + stat_smooth(method="lm", colour = "black") + geom_point(aes(size = log10(module.size), fill = mean.sex.bias),shape=21) + facet_wrap(~network)  + xlab("Selection estimate in females") +  ylab("Selection estimate in males") + theme_bw() + theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.background = element_rect(fill="white", color="black",size=1), panel.border = element_rect(colour="black",size=1), axis.title = element_text(size=12), strip.text = element_text(size=12)) + scale_fill_gradient2(low="blue", mid="white", high="red") + scale_x_continuous(limits = c(-0.2, 0.2))
dev.off()

# Stats for Figure 2:
summary(lm(beta.male ~ beta.female, data = plot.data[plot.data$network == "Consensus network", ]))
summary(lm(beta.male ~ beta.female, data = plot.data[plot.data$network == "Female network", ]))
summary(lm(beta.male ~ beta.female, data = plot.data[plot.data$network == "Male network", ]))



# FIGURE 3 - male-beneficial and female-beneficial transcripts are grouped into modules, and these modules are more SA.
pdf("../figures/modularity.of.SA.pdf", width = 9.5, height = 8)
modularity.of.SA.plot()
dev.off()



# FIGURE 4: sexually antagonistic modules tend to be more heritable
pdf("../figures/herit.pdf", width = 8.5, height = 3.7)
ggplot(plot.data, aes(x = SA.score, y = Herit.to.plot)) + stat_smooth(method="lm", colour = "black") + geom_vline(xintercept = 0, linetype = 2, colour = "grey") + geom_point(aes(size = module.size, fill = mean.sex.bias),shape=21)  + xlab("Sex-specific selection index (I)") + ylab("Heritability") +theme_bw() + theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour="black",size=1), axis.title = element_text(size=12), strip.text = element_text(size=12),strip.background = element_rect(fill="white", color="black",size=1))+ scale_fill_gradient2(low="blue", mid="white", high="red") + facet_wrap(~network) + coord_cartesian(ylim = c(0, 1.02))
dev.off()

# Stats for Figure 4:
summary(lm(Herit.to.plot ~ SA.score, data = plot.data[plot.data$network == "Consensus network", ]))
summary(lm(Herit.to.plot ~ SA.score, data = plot.data[plot.data$network == "Female network", ]))
summary(lm(Herit.to.plot ~ SA.score, data = plot.data[plot.data$network == "Male network", ]))

summary(lm(Herit.to.plot ~ mean.sex.bias + module.size, data = plot.data[plot.data$network == "Consensus network", ]))
summary(lm(Herit.to.plot ~ mean.sex.bias + module.size, data = plot.data[plot.data$network == "Female network", ]))
summary(lm(Herit.to.plot ~ mean.sex.bias + module.size, data = plot.data[plot.data$network == "Male network", ]))



# FIGURE 5: Modules with sex-specific genetic architecture tend to be less SA
pdf("../figures/inter.sex.corr.pdf", width = 3.7, height = 3.85)
ggplot(plot.data[plot.data$network=="Consensus network", ], aes(x = gen.corr, y = SA.score)) +geom_hline(yintercept = 0,linetype=2,colour="grey") +geom_vline(xintercept = 0,linetype=2,colour="grey")+ stat_smooth(method="lm", colour = "black") + geom_point(aes(size = module.size, fill = mean.sex.bias),shape=21)+  ylab("Sex-specific selection index (I)") +theme_bw() + theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour="black",size=1), axis.title = element_text(size=12), strip.text = element_text(size=12),strip.background = element_rect(fill="white", color="black",size=1))+ scale_fill_gradient2(low="blue", mid="white", high="red")  + xlab("Inter-sex genetic correlation\nfor module expression")+ scale_y_continuous(breaks = c(-0.1,-0.05,0,0.05,0.1)) + scale_x_continuous(breaks=c(-0.4,0,0.4,0.8))
dev.off()

# stats associated with Figure 5
summary(lm(SA.score ~ gen.corr, data = plot.data[plot.data$network=="Consensus network", ])) 
summary(lm(gen.corr ~  mean.sex.bias, data = plot.data[plot.data$network=="Consensus network", ])) 
summary(lm(SA.score ~  mean.sex.bias, data = plot.data[plot.data$network=="Consensus network", ])) 


# Figure 6 - export at 28 x 6.75 then finish in Inkscape
# Histogram of I values and associated GO enrichment tables
t1 <- tableGrob(go.output.SA[-(1:2)][,c(5,1:4)] %>% mutate(OddsRatio = round(OddsRatio,2), ExpCount = round(ExpCount, 2)), rows = NULL, base_size = 1)
t2 <- tableGrob(go.output.concordant[-(1:2)][,c(5,1:4)] %>% mutate(OddsRatio = round(OddsRatio,2), ExpCount = round(ExpCount, 2)), rows = NULL, base_size = 1)
p <- ggplot(gene.data, aes(x=SA.score)) +   geom_histogram(aes(y = ..count..), bins = 8, breaks = c(-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2), fill = c(cols[1], cols), colour = "black")+
  annotate("text", x = -0.175, y = 119, size = 4, label = paste(range.diff(cum.sum$percentile[cum.sum$I< -0.15]), "%", sep="")) + 
  annotate("text", x = -0.125, y = 525, size = 4, label = paste(range.diff(cum.sum$percentile[cum.sum$I> -0.15 & cum.sum$I< -0.1]), "%", sep="")) + 
  annotate("text", x = -0.075, y = 2210, size = 4, label = paste(range.diff(cum.sum$percentile[cum.sum$I> -0.1 & cum.sum$I< -0.05]), "%", sep="")) + 
  annotate("text", x = -0.025, y = 4499, size = 4, label = paste(range.diff(cum.sum$percentile[cum.sum$I> -0.05 & cum.sum$I< 0]), "%", sep="")) + 
  annotate("text", x = 0.025, y = 4160, size = 4, label = paste(range.diff(cum.sum$percentile[cum.sum$I> 0 & cum.sum$I< 0.05]), "%", sep="")) + 
  annotate("text", x = 0.075, y = 1470, size = 4, label = paste(range.diff(cum.sum$percentile[cum.sum$I> 0.05 & cum.sum$I< 0.1]), "%", sep="")) + 
  annotate("text", x = 0.125, y = 230, size = 4, label = paste(range.diff(cum.sum$percentile[cum.sum$I> 0.1 & cum.sum$I< 0.15]), "%", sep="")) + 
  annotate("text", x = 0.175, y = 83, size = 4, label = paste(range.diff(cum.sum$percentile[cum.sum$I> 0.15]), "%", sep="")) + 
  theme_bw(15) + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_rect(colour="black", fill = NA, size=1)) + 
  xlab("Index of sex-specific selection (I)") + ylab("Number of transcripts") + scale_y_continuous(expand = c(0,0), limits = c(0,4600)) + scale_x_continuous(breaks = c(-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2)) 
grid.newpage() 
grid.draw(arrangeGrob(t1,p,t2,ncol=3, as.table=T))
rm(list = c("t1", "t2", "p"))


# Removed from current draft of paper
# # FIGURE 6: Connectivity of SA and non-SA genes, split by gene type (male-biased, female-biased, or unbiased)
# pdf("../figures/connectivity.plot.pdf", width = 9.3, height = 3.7)
# gene.connectivity.plot()
# dev.off()
# Stats associated with Figure 6 (note that this treats the extent of sex bias, and the extent of SA, as continuous variables - they are binned into discrete categories in the figure for ease of presentation)
# xtable(Anova(lm(gene.data$both.kTotal / max(gene.data$both.kTotal) ~ SA.score * sex.bias, data = gene.data), type = "III"))
# xtable(Anova(lm(gene.data$female.kTotal / max(gene.data$female.kTotal)~ SA.score * sex.bias, data = gene.data), type = "III"))
# xtable(Anova(lm(gene.data$male.kTotal / max(gene.data$male.kTotal) ~ SA.score * sex.bias, data = gene.data), type = "III"))


################################################
# SUPPLEMENTARY TABLES
################################################

# GO term results
make.latex.table(GO.terms.per.module.consensus, "c", 10)
make.latex.table(GO.terms.per.module.female, "f", 10)
make.latex.table(GO.terms.per.module.male, "m", 10)

# Supplementary table listing the top 1% most SA genes
supp.table <- gene.data %>% arrange(SA.score)
supp.table <- supp.table[supp.table$SA.score <= quantile(supp.table$SA.score, probs = 0.01), c(1,9,7,16,18)] 
names(supp.table) <- c("Probe", "Gene name", "Entrez ID", "Benefiting sex", "Sex bias")
supp.table[,5] <- round(supp.table[,5], 2)
supp.table[grep("Glycosylphosphatidylinositol anchor attachment", supp.table[,2]),2] <- "Glycosylphosphatidylinositol anchor attachment 1 ortholog" # shorten this one's name
print(xtable(supp.table), include.rownames=F, tabular.environment = "longtable", floating = F)

# Supplementary table listing the top 1% most concordant genes
supp.table <- gene.data %>% arrange(-SA.score)
supp.table <- supp.table[supp.table$SA.score >= quantile(supp.table$SA.score, probs = 0.99), c(1,9,7,18)]
names(supp.table) <- c("Probe", "Gene name", "Entrez ID", "Benefiting sex", "Sex bias")
supp.table[,5] <- round(supp.table[,5], 2)
print(xtable(supp.table), include.rownames=F, tabular.environment = "longtable", floating = F)

# For fun, have a look at all the GO terms associated with the top 1% most SA genes 
get.GO.for.SA.genes(supp.table)



################################################
# SUPPLEMENTARY FIGURES
################################################

# Here is a histogram of I values for all the genes, as calculated by me using hemiclone means selection analysis
cols <- brewer.pal(7, "RdBu")
pdf("../figures/SA.genes.histogram.pdf", width = 7.5, height = 6.75)
ggplot(gene.data, aes(x=SA.score)) +   geom_histogram(aes(y = ..count..), bins = 8, breaks = c(-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2), fill = c(cols[1], cols), colour = "black")+
  annotate("text", x = -0.175, y = 119, size = 4, label = paste(range.diff(cum.sum$percentile[cum.sum$I< -0.15]), "%", sep="")) + 
  annotate("text", x = -0.125, y = 525, size = 4, label = paste(range.diff(cum.sum$percentile[cum.sum$I> -0.15 & cum.sum$I< -0.1]), "%", sep="")) + 
  annotate("text", x = -0.075, y = 2210, size = 4, label = paste(range.diff(cum.sum$percentile[cum.sum$I> -0.1 & cum.sum$I< -0.05]), "%", sep="")) + 
  annotate("text", x = -0.025, y = 4499, size = 4, label = paste(range.diff(cum.sum$percentile[cum.sum$I> -0.05 & cum.sum$I< 0]), "%", sep="")) + 
  annotate("text", x = 0.025, y = 4160, size = 4, label = paste(range.diff(cum.sum$percentile[cum.sum$I> 0 & cum.sum$I< 0.05]), "%", sep="")) + 
  annotate("text", x = 0.075, y = 1470, size = 4, label = paste(range.diff(cum.sum$percentile[cum.sum$I> 0.05 & cum.sum$I< 0.1]), "%", sep="")) + 
  annotate("text", x = 0.125, y = 230, size = 4, label = paste(range.diff(cum.sum$percentile[cum.sum$I> 0.1 & cum.sum$I< 0.15]), "%", sep="")) + 
  annotate("text", x = 0.175, y = 83, size = 4, label = paste(range.diff(cum.sum$percentile[cum.sum$I> 0.15]), "%", sep="")) + 
  theme_bw(15) + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_rect(colour="black", fill = NA, size=1)) + 
  xlab("Index of sex-specific selection (I)") + ylab("Number of transcripts") + scale_y_continuous(expand = c(0,0), limits = c(0,4600)) + scale_x_continuous(breaks = c(-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2))
dev.off()




# FIGURES S1-S3: relative fitness against standardized mean eigengene for each module, for males and females, at the level of hemiclones
fitness.vs.eigengenes.plot(MEconsensus.clone.data, module.data.consensus, consensus=T)
fitness.vs.eigengenes.plot(MEfemale.clone.data, module.data.female)
fitness.vs.eigengenes.plot(MEmale.clone.data, module.data.male)

pdf("../figures/fitness.consensus.pdf", width = 8, height = 8); fitness.vs.eigengenes.plot(MEconsensus.clone.data, module.data.consensus, consensus=T); dev.off()
pdf("../figures/fitness.female.pdf", width = 7, height = 8); fitness.vs.eigengenes.plot(MEfemale.clone.data, module.data.female); dev.off()
pdf("../figures/fitness.male.pdf", width = 8, height = 7); fitness.vs.eigengenes.plot(MEmale.clone.data, module.data.male); dev.off()

# FIGURES S4-S6: Tissue enrichment plots
tissue.enrichment.plot(module.data.consensus, "consensus.module")
tissue.enrichment.plot(module.data.female, "female.module")
tissue.enrichment.plot(module.data.male, "male.module")

# Write the tissue enrichment plots to disk (it's complicated to save them because I used the grid package)
pdf("../figures/consensus.tissue.plot%03d.pdf", height = 10, width = 10, onefile = F); tissue.enrichment.plot(module.data.consensus, "consensus.module"); dev.off()
file.rename("../figures/consensus.tissue.plot003.pdf", "../figures/consensus.tissue.plot.pdf")
unlink(c("../figures/consensus.tissue.plot001.pdf", "../figures/consensus.tissue.plot002.pdf"))

pdf("../figures/female.tissue.plot%03d.pdf", height = 10, width = 10, onefile = F); tissue.enrichment.plot(module.data.female, "female.module"); dev.off()
file.rename("../figures/female.tissue.plot003.pdf", "../figures/female.tissue.plot.pdf")
unlink(c("../figures/female.tissue.plot001.pdf", "../figures/female.tissue.plot002.pdf"))

pdf("../figures/male.tissue.plot%03d.pdf", height = 7, width = 10, onefile = F); tissue.enrichment.plot(module.data.male, "male.module"); dev.off()
file.rename("../figures/male.tissue.plot003.pdf", "../figures/male.tissue.plot.pdf")
unlink(c("../figures/male.tissue.plot001.pdf", "../figures/male.tissue.plot002.pdf"))


# FIGURES S7-S9: Simulation to determine if genes from the same module are clustered in the genome

# Run a simulation to check whether genes from the same module are more, or less, clustered on the chromosomes than expected by chance
# The simulation compares the observed distribution of nearest-neighbour distances to the one that expected if genes in modules are distributed randomly within chromosomes
consensus.boots <- distance.simulation("consensus.module", boots = 10000)
female.boots <- distance.simulation("female.module", boots = 10000)
male.boots <- distance.simulation("male.module", boots = 10000)

boot.distance.plot <- function(boot.data) ggplot(boot.data, aes(x = Module, y = mean.difference)) + geom_hline(yintercept = 0, colour = 'grey') + geom_point(size=1.8) + geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), width=0) + facet_wrap(~Chromosome) + ylab("Degree of clustering relative to the null expectation") + theme_bw(15) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour="black", size=1), strip.background = element_rect(colour="black", size=1))

ggsave("../figures/distance plot - consensus.eps", boot.distance.plot(consensus.boots), width = 13.1, height = 9.51)
ggsave("../figures/distance plot - female.eps", boot.distance.plot(female.boots), width = 13.1, height = 9.51)
ggsave("../figures/distance plot - male.eps", boot.distance.plot(male.boots), width = 13.1, height = 9.51)





# Snipped text from the paper about connectivity of SA and non-SA genes - I thought this was confusing, and the predictions are not clearly defined enough.

# %\vspace{4mm} 
# %\subsection*{Calculating connectivity for each transcript}
# %I calculated 'connectivity' for each transcript, separately for the consensus, male, and female networks. Connectivity is the sum of the absolute correlations in expression between the focal gene and all the other genes, after scaling the transcriptional network to have scale-independent topology. For the consensus network this was accomplished by making a $n\times n$ matrix of absolute correlation coefficients (where $n$ is the number of genes), raising each correlation to the power $p_{c}$ (where $p_{c}$ is the scaling power used when constructing the consensus network), and then summing the off-diagonal correlations in each column. The procedure was the same for the female and male networks, except that the adjacency matrix was constructed using only samples of the correct sex, and was scaled using the appropriate powers $p_{f}$ and $p_{m}$.

# \subsection*{Connectivity of sexually antagonistic and non-antagonistic genes}
# Next, I tested whether sexually antagonistic genes were more or less connected than non-antagonistic genes. Figure \ref{fig:connectedness_plot} shows connectivity (expressed as a proportion of the highest connectivity value in the network) for antagonistic and non-antagonistic genes, which have been further divided into genes with male-biased expression, female-biased expression, or no sex bias in expression.
# 
# For all three networks, sexual antagonism was significantly related to connectivity, though the strength of the relationship differed between between male-biased, female-biased or unbiased genes (Tables S4-S6). However, the amount of variation in connectivity explained by sexual antagonism was uniformly low (Tables S4-S6), as illustrated by the subtlety of the differences in Figure \ref{fig:connectedness_plot} and the large 95\% quantiles around the means. 
# 
# In the consensus network, antagonistic genes showing male- or female-biased expression were less connected. In the female network the negative effect of antagonism on connectivity was most apparent for female-biased genes, while in the male network the same was true of male-biased genes. 
# 
# \begin{figure*}[hb]
# \begin{center}
# \includegraphics[width=1.0\textwidth]{connectivity_plot.pdf}
# \end{center}
# \caption{The connectivity of a transcript depended on whether it was sexually antagonistic, whether it showed sex-biased expression, and the interaction between these two parameters. The genes have been binned into discrete categories for ease of interpretation: SA genes were arbitrarily defined as those with a sexual antagonism score $>0.05$, while sex-biased genes were defined as those with $>2$-fold differences in expression between sexes. The plots show the mean and 95\% quantiles of relative connectivity.}
# \label{fig:connectedness_plot}
# \end{figure*}
