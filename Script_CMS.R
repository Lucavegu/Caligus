###### Sea Lice Prediction Paper (CMS reproduction/analysis)

########################################################################################
#### Import appropriate R objects and/or data files
#### (also perform checks and any cleaning if necessary)
########################################################################################

##OTU TABLE
library(phyloseq)
Skin_ps_reordered <- readRDS("~/Skin_rare_SILVA.rds")
Skin_ps_reordered
#otu_table()   OTU Table:         [ 125 taxa and 107 samples ]
#sample_data() Sample Data:       [ 107 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 125 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 125 tips and 124 internal nodes ]
#refseq()      DNAStringSet:      [ 125 reference sequences ]

head(sample_data(Skin_ps_reordered))

## Arrange by sample_name

#install.packages("ggtext") # for rotated labels on ord_plot() 
#install.packages("ggraph") # for taxatree_plots()
#install.packages("DT") # for tax_fix_interactive()
#install.packages("corncob") # for beta binomial models in tax_model()
library(dplyr)
library(microViz)

Skin_ps_reordered <- Skin_ps_reordered %>% ps_arrange(Sample_name)
#otu_table()   OTU Table:         [ 125 taxa and 107 samples ]
#sample_data() Sample Data:       [ 107 samples by 17 sample variables ]
#tax_table()   Taxonomy Table:    [ 125 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 125 tips and 124 internal nodes ]
#refseq()      DNAStringSet:      [ 125 reference sequences ]

head(sample_data(Skin_ps_reordered))

## Try splitting the clearly bimodal distribution into two, then use
## means and std devs to discretize into "high" and "low" sea lice loads

library(VGAM)

# Will use Total_caligus (count) per Final_Length (cm) as a reasonable metric for "load"
tc_per_cm <- as.numeric(sample_data(Skin_ps_reordered)$Total_caligus)/as.numeric(sample_data(Skin_ps_reordered)$Final_Length)

# The mix2normal function will infer via ML the two normal distributions in the data
fit <- vglm(tc_per_cm ~ 1, mix2normal(eq.sd=FALSE))

fit

# Parameters from vglm and mix2normal fit
pars = as.vector(coef(fit))
w = logit(pars[1], inverse=TRUE)
m1 = pars[2]
sd1 = exp(pars[3])
m2 = pars[4]
sd2 = exp(pars[5])

hist(tc_per_cm, 30, col="grey90", freq=FALSE,
     main=NA,
     #main="Histogram of Lice load (total count per cm) with Gaussian Fit",
     xlab="Lice load", ylab="Density")

x <- seq(0, 2.5, 0.01)
lines(x, w*dnorm(x, m1, sd1) + (1-w)*dnorm(x, m2, sd2),
      col="red", lwd=2)

abline(v=c(m1, m2), col=c("aquamarine", "green"), lwd=2)
abline(v=c(m1+sd1, m1-sd1, m2+sd2, m2-sd2),
       col=c(rep("aquamarine",2), rep("green",2)),
       lty=2)

legend("topright",
       legend=c("Gaussian Fit", "Mean 1", "Mean 2", "±1 SD"),
       col=c("red", "aquamarine", "green", "grey40"),
       lty=c(1,1,1,2), lwd=c(2,1,1,1), bty="n")



# Superimpose the fitted distribution
x <- seq(0, 2.5, 0.01)
points(x, w*dnorm(x, m1, sd1)+(1-w)*dnorm(x,m2,sd2), "l", col="red", lwd=2)

abline(v=c(m1,m2), col=c("aquamarine","green"))
abline(v=c(m1+sd1,m1-sd1,
           m2+sd2,m2-sd2),
       col=c(rep("aquamarine",2),
             rep("green",2)), lty=2)

# proposed "cut-offs" are 1 sd above mean 1 and 1 sd below mean 2
m1+sd1 #0.7006043
m2-sd2 #1.066726


sum(tc_per_cm>=m2-sd2) #42 fish in the "High" group

sum(tc_per_cm<=m1+sd1) #48 fish in the "Low" group


## Update the original ps object to reflect the fixed categorization
metadata_fixed <- Skin_ps_reordered
sample_data(metadata_fixed)$liceLoad <- tc_per_cm

sample_data(metadata_fixed)$Burden <- ifelse(sample_data(metadata_fixed)$liceLoad<=m1+sd1, "Low",
                                             ifelse(sample_data(metadata_fixed)$liceLoad>=m2-sd2, "High", "Moderate"))


stripchart(sample_data(metadata_fixed)$liceLoad ~ sample_data(metadata_fixed)$Burden,
           vertical=TRUE, method="jitter", jitter=0.2, cex=0.4,
           ylab="Lice load")

summary(as.factor(sample_data(metadata_fixed)$Burden))
#High   Low   Moderate 
#42     48       17 


## High AND LOW BURDEN
High_Low_burden <- subset_samples(metadata_fixed, Burden!="Moderate")
#Arrange by sample_name
High_Low_microbiome <- High_Low_burden %>% ps_arrange(Sample_name)
#otu_table()   OTU Table:         [ 125 taxa and 90 samples ]
#sample_data() Sample Data:       [ 90 samples by 18 sample variables ]
#tax_table()   Taxonomy Table:    [ 125 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 125 tips and 124 internal nodes ]
#refseq()      DNAStringSet:      [ 125 reference sequences ]

##Figure 1. ABUNDANCE PLOT

metadata_skin <- sample_data(High_Low_microbiome)
# Convert 'Sex' from character to factor and Total_caligus as numeric
metadata_skin$Sex <- as.factor(metadata_skin$Sex)
metadata_skin$Total_caligus <- as.numeric(metadata_skin$Total_caligus)
# Reassign the modified metadata back to the phyloseq object
sample_data(High_Low_microbiome) <- metadata_skin

#SKIN samples
#from absolute to compositional data
ps.rel <- transform_sample_counts(High_Low_microbiome, function(x) x / sum(x))

#Aggregate at the Phylum level
physeq_phylum <- tax_glom(ps.rel, taxrank = "Phylum")

#Create a table for Phylum abundance
phylum_abund <- data.frame(tax_table(physeq_phylum), 
                           Abundance = taxa_sums(physeq_phylum))

#Calculate the total abundance (sum of all taxa)
total_abundance <- sum(phylum_abund$Abundance)

#Add a column with the percentage of each phylum
phylum_abund$Percentage <- (phylum_abund$Abundance / total_abundance) * 100

#Sort the table by abundance (descending order)
phylum_abund <- phylum_abund[order(-phylum_abund$Abundance), ]

#Save table

write.table(phylum_abund, "/home/lucas/Caligus_microbiome/First_sequencing/phylum_abund.tsv", sep = "\t", quote = FALSE, col.names = NA)

#Get the top 10 phyla
top_10_phyla <- head(phylum_abund, 10)

#Calculate the Most Abundant Taxon for Each Sample:
#Melt the phyloseq object
ps_melt <- psmelt(ps.rel)

#Summarize the total abundance of each phylum across all samples
phylum_total_abundance <- ps_melt %>%
  group_by(Phylum) %>%
  summarise(TotalAbundance = sum(Abundance)) %>%
  arrange(desc(TotalAbundance))

#Select the top 10 phyla
top_10_phyla <- phylum_total_abundance %>%
  top_n(10, wt = TotalAbundance) %>%
  pull(Phylum)

#Filter the ps_melt data to include only the top 10 phyla
ps_melt_top_10 <- ps_melt %>%
  filter(Phylum %in% top_10_phyla)

#Calculate the most abundant taxon for each sample (top taxon per sample)
top_taxa_per_sample <- ps_melt_top_10 %>%
  group_by(Sample, Phylum) %>%
  summarise(Abundance = sum(Abundance)) %>%
  top_n(1, wt = Abundance) %>%
  arrange(desc(Abundance))

#Order the samples by the most abundant taxon
sample_order <- top_taxa_per_sample %>%
  arrange(desc(Abundance)) %>%
  pull(Sample)

# Update sample order in ps_melt_top_10
ps_melt_top_10$Sample <- factor(ps_melt_top_10$Sample, levels = sample_order)

#Plot the relative abundance for the top 10 phyla
color_palette <- c("#8F6F00",  # Actinomycetota
                   "#01843a",  # Bacillota
                   "#fed976",  # Bacteroidota
                   "#FFA73B",  # Bdellovibrionota
                   "#D7301F",  # Patescibacteria
                   "#225EA8",  # Pseudomonadota
                   "#35978F")  # Verrucomicrobiota


p <- ggplot(ps_melt_top_10, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_x_discrete(limits = sample_order) +
  labs(
    y = "Relative Abundance",
    fill = "Phylum") +
  scale_fill_manual(values = color_palette) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right") # Adjust legend position as needed

# Print the plot
print(p)

#Aggregate at the Genus level
physeq_genus <- tax_glom(ps.rel, taxrank = "Genus")

#Create a table for Phylum abundance
genus_abund <- data.frame(tax_table(physeq_genus), 
                          Abundance = taxa_sums(physeq_genus))

#Calculate the total abundance (sum of all taxa)
total_abundance <- sum(genus_abund$Abundance)

#Add a column with the percentage of each phylum
genus_abund$Percentage <- (genus_abund$Abundance / total_abundance) * 100

#Sort the table by abundance (descending order)
genus_abund <- genus_abund[order(-genus_abund$Abundance), ]

#Save genus table
write.table(genus_abund, "/home/lucas/Caligus_microbiome/First_sequencing/genus_abund.tsv", sep = "\t", quote = FALSE, col.names = NA)

#Get the top 10 genera
top_10_genera <- head(genus_abund, 10)

#Calculate the Most Abundant Taxon for Each Sample:
#Melt the phyloseq object
ps_melt <- psmelt(ps.rel)

#Summarize the total abundance of each genus across all samples
genus_total_abundance <- ps_melt %>%
  group_by(Genus) %>%
  summarise(TotalAbundance = sum(Abundance)) %>%
  arrange(desc(TotalAbundance))

#Select the top 10 genera
top_10_genera <- genus_total_abundance %>%
  top_n(10, wt = TotalAbundance) %>%
  pull(Genus)

#Filter the ps_melt data to include only the top 10 genera
ps_melt_top_10_genus <- ps_melt %>%
  filter(Genus %in% top_10_genera)

#Calculate the most abundant taxon for each sample (top taxon per sample)
top_taxa_per_sample_genus <- ps_melt_top_10_genus %>%
  group_by(Sample, Genus) %>%
  summarise(Abundance = sum(Abundance)) %>%
  top_n(1, wt = Abundance) %>%
  arrange(desc(Abundance))

#Order the samples by the most abundant taxon
sample_order_genus <- top_taxa_per_sample_genus %>%
  arrange(desc(Abundance)) %>%
  pull(Sample)

# Update sample order in ps_melt_top_10
ps_melt_top_10_genus$Sample <- factor(ps_melt_top_10_genus$Sample, levels = sample_order)

#Plot the relative abundance for the top 10 genera
color_palette <- c("#1f77b4",  # Blue
                   "#ff7f0e",  # Orange
                   "#2ca02c",  # Green
                   "#d62728",  # Red
                   "#9467bd",  # Purple
                   "#8c564b",  # Brown
                   "#e377c2",  # Pink
                   "#7f7f7f",  # Gray
                   "#bcbd22",  # Olive
                   "#17becf")  # Cyan

q <- ggplot(ps_melt_top_10_genus, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_x_discrete(limits = sample_order) +
  labs(
    y = "Relative Abundance",
    fill = "Genus") +
  scale_fill_manual(values = color_palette) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right") # Adjust legend position as needed

# Print the plot
print(q)

#WATER
#from absolute to compositional data
Ps_water <- subset_samples(bacteria_physeq, Pittag =="CTRL")
ps.rel_water <- transform_sample_counts(Ps_water, function(x) x / sum(x))

#Aggregate at the Phylum level
physeq_phylum_water <- tax_glom(ps.rel_water, taxrank = "Phylum")

#Create a table for Phylum abundance
phylum_abund_water <- data.frame(tax_table(physeq_phylum_water), 
                                 Abundance = taxa_sums(physeq_phylum_water))

#Calculate the total abundance (sum of all taxa)
total_abundance_water <- sum(phylum_abund_water$Abundance)

#Add a column with the percentage of each phylum
phylum_abund_water$Percentage <- (phylum_abund_water$Abundance / total_abundance_water) * 100

#Sort the table by abundance (descending order)
phylum_abund_water <- phylum_abund_water[order(-phylum_abund_water$Abundance), ]

#Save table
write.table(phylum_abund_water, "/home/lucas/Caligus_microbiome/First_sequencing/phylum_abund_water.tsv", sep = "\t", quote = FALSE, col.names = NA)

#Get the top 10 phyla
top_10_phyla_water <- head(phylum_abund_water, 10)

#Calculate the Most Abundant Taxon for Each Sample:
#Melt the phyloseq object
ps_melt_water <- psmelt(ps.rel_water)

#Summarize the total abundance of each phylum across all samples
phylum_total_abundance_water <- ps_melt_water %>%
  group_by(Phylum) %>%
  summarise(TotalAbundance = sum(Abundance)) %>%
  arrange(desc(TotalAbundance))

#Select the top 10 phyla
top_10_phyla_water <- phylum_total_abundance_water %>%
  top_n(10, wt = TotalAbundance) %>%
  pull(Phylum)

#Filter the ps_melt data to include only the top 10 phyla
ps_melt_top_10_water <- ps_melt_water %>%
  filter(Phylum %in% top_10_phyla_water)

#Calculate the most abundant taxon for each sample (top taxon per sample)
top_taxa_per_sample_water <- ps_melt_top_10_water %>%
  group_by(Sample, Phylum) %>%
  summarise(Abundance = sum(Abundance)) %>%
  top_n(1, wt = Abundance) %>%
  arrange(desc(Abundance))

#Order the samples by the most abundant taxon
sample_order_water <- top_taxa_per_sample_water %>%
  arrange(desc(Abundance)) %>%
  pull(Sample)

# Update sample order in ps_melt_top_10
ps_melt_top_10_water$Sample <- factor(ps_melt_top_10_water$Sample, levels = sample_order_water)

#Plot the relative abundance for the top 10 phyla
color_palette <- c("#8F6F00",  # Actinomycetota
                   "#01843a",  # Bacillota
                   "#fed976",  # Bacteroidota
                   "#FFA73B",  # Bdellovibrionota
                   "#D7301F",  # Patescibacteria
                   "#225EA8",  # Pseudomonadota
                   "#35978F")  # Verrucomicrobiota

p1 <- ggplot(ps_melt_top_10_water, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_x_discrete(limits = sample_order_water) +
  labs(
    y = "Relative Abundance",
    fill = "Phylum") +
  scale_fill_manual(values = color_palette) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right") # Adjust legend position as needed

# Print the plot
print(p1)

par(mfrow = c(1, 2))
print(p1)
print(p)

# Combine the CTRL (unrarefied) and SG (rarefied) datasets

# ================================
# Load Required Libraries
# ================================
library(phyloseq)
library(biomformat)
library(ggplot2)
library(reshape2)
library(readr)
library(dplyr)
library(pheatmap)
library(ggpubr)

merged_physeq <- merge_phyloseq(Ps_water, High_Low_microbiome)


# ================================
# 1. Export OTU/ASV Table
# ================================
# Ensure your phyloseq object is called 'merged_physeq'
otu <- as(otu_table(merged_physeq), "matrix")
#if (taxa_are_rows(merged_physeq)) {
#  otu <- t(otu)
#}

# Save as tab-delimited text
write.table(otu, file = "~/otu_table.txt", sep = "\t", quote = FALSE)

# ================================
# 2. Export Sample Metadata
# ================================

metadata2 <- sample_data(merged_physeq)
Sample <- c(rep("Water",5),rep("Skin",90))
# Add to sample_data
sample_data(merged_physeq)$Type <- Sample


# ================================
# 3. Convert OTU Table to BIOM Format
# ================================
biom_obj <- make_biom(data = otu)
write_biom(biom_obj, biom_file = "~/otu_table.biom")

# Optional check
read_biom("~/otu_table.biom")

# ================================
# 4. Run SourceTracker2 (shell)
# ================================
# Run this line in the terminal, not in R
# Replace "~/Skin_water/" with your actual output folder
# system("sourcetracker2 gibbs -i ~/otu_table.biom -m ~/metadata.txt --sink_rarefaction_depth 0 --source_rarefaction_depth 0 -o ~/Skin_water/ --jobs 8")


## JUST ASV (Original)
#metadata
metadata_new <- High_Low_microbiome@sam_data

boxplot(metadata_new$liceLoad ~ metadata_new$Burden)

stripchart(metadata_new$liceLoad ~ metadata_new$Burden,
           vertical=TRUE, method="jitter", jitter=0.2)



#vector of discrete trait coded 1=EDD cases and 0=Healthy controls
y <- as.numeric(metadata_new$Burden=="High") #Low = 0; High = 1


## Save objects to load for future analysis
save(High_Low_microbiome, y,
     file = "~/Clay_analysis/CMS_90samples_new.RData")



########################################################################################
#### Prepare ASV table for ML classification models
#### 
########################################################################################

#Microbiome classification
#READ: https://microbiome.netlify.app/classification-using-microbiome
#READ: https://microbiome.github.io/course_2022_FindingPheno/supervised-learning.html

library(caret) ## Most popular R package for machine learning
library(randomForest)
library(pROC)
library(stringr)
library(xgboost)
library(e1071)
library(glmnet)
library(MASS)
library(neuralnet)
library(pROC)


#otu_table
otu_table <- as.data.frame(High_Low_microbiome@otu_table)
#check that order of samples in metadata and OTU table are identical
colnames(otu_table) <- gsub("[.]","-",colnames(otu_table))
otu_table <- otu_table[,rownames(metadata_new)]
#Selecting 500 more abundant ASVs
#Rank ASVs by total abundance
# Calculate the total abundance for each ASV
#total_abundance <- rowSums(otu_table)
# Select the 500 most abundant ASVs
# Order the ASVs by total abundance in descending order
#otu_table <- otu_table[order(total_abundance, decreasing = TRUE), ]
# Select the top 1,000 ASVs
#top_500_asvs <- otu_table[1:500, ]

x <- data.matrix(otu_table)
#calculate relative abundance of each OTU by dividing its value by the total reads per sample
x <- sweep(x,2,colSums(x),"/")
x <- t(x)
#Center the matrix
#x1<-scale(x,center=TRUE,scale=F)
x1<-x
id <- seq(1:length(y))
numsets <- 100



########################################################################################
#### Set up and run 5-fold cross validation for 8 ML classification models
#### (ASV count matrix transformation by total-sum scaling and NO centering)
########################################################################################

#actual y representing 100 random samples of sample IDs
y.sample.matrix <- matrix(nrow=length(y),ncol=numsets)
y.test.matrix <- matrix(nrow=length(y),ncol=numsets)
#predictions are probabilitities that y=1 for each sample
rf.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
xgboost.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
svm.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
lasso.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
ridge.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
enet.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
knn.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
#neural.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
lda.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
#hclust.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
#kmeans.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
set.seed(1000)


for (n in 1:numsets) {
  print(n)
  #randomly split samples into five groups for 5-fold cross-validation
  #recall that there are 90 (not 70) samples in this version of the data set
  subset <- matrix(0,nrow=18,ncol=5)
  subset[1:18,1] <- sort(sample(id,18))
  subset[1:18,2] <- sort(sample(which(!(id %in% subset)),18))
  subset[1:18,3] <- sort(sample(which(!(id %in% subset)),18))
  subset[1:18,4] <- sort(sample(which(!(id %in% subset)),18))
  subset[1:18,5] <- which(!(id %in% subset))
  y.test <- NULL
  rf.pred <- NULL
  xgboost.pred <- NULL
  svm.pred <- NULL
  lasso.pred <- NULL
  ridge.pred <- NULL
  enet.pred <- NULL
  knn.pred <- NULL
  #neural.pred <- NULL
  lda.pred <- NULL
  #use 4 groups as training set and 1 group as testing set for each supervised method, and repeat until each group is tested
  for (i in 1:5) {
    print(i)
    y.test <- c(y.test, y[subset[,i]])
    y.train <- y[-subset[,i]]
    x.train <- x1[-subset[,i],]
    #remove OTUs that have all zeroes in training set
    x.train <- subset(x.train,select=apply(x.train,2,var)>0)
    x.test <- x1[subset[,i],colnames(x.train)]
    #random forest
    rf.train <- randomForest(x.train,y=as.factor(y.train), importance = TRUE)
    rf.pred <- c(rf.pred,predict(rf.train,x.test,type="prob")[,2])
    #gradient boosting
    xgboost.train <- xgboost(x.train,label=y.train,nrounds=10,objective="binary:logistic",verbose=0)
    xgboost.pred <- c(xgboost.pred,predict(xgboost.train,x.test))
    #SVM
    svm.train <- svm(x.train,as.factor(y.train),probability=T)
    svm.model <- predict(svm.train,x.test,probability=T)
    svm.pred <- c(svm.pred,attr(svm.model,"probabilities")[,1])
    #lasso
    lasso.train <- train(x.train,as.factor(y.train),method="glmnet",tuneGrid=expand.grid(.alpha=1,.lambda=seq(0,1,by=0.1)))
    lasso.pred <- c(lasso.pred,predict(lasso.train,x.test,type="prob")[,2])
    #ridge
    ridge.train <- train(x.train,as.factor(y.train),method="glmnet",tuneGrid=expand.grid(.alpha=0,.lambda=seq(0,1,by=0.1)))
    ridge.pred <- c(ridge.pred,predict(ridge.train,x.test,type="prob")[,2])
    #elastic net
    enet.train <- train(x.train,as.factor(y.train),method="glmnet")
    enet.pred <- c(enet.pred,predict(enet.train,x.test,type="prob")[,2])
    #k-nearest neighbors
    knn.train <- knn3(x.train,as.factor(y.train))
    knn.pred <- c(knn.pred,predict(knn.train,x.test,type="prob")[,2])
    #neural networks
    #train <- data.frame(y.train,x.train)
    #names <- names(train)
    #formula <- as.formula(paste("y.train ~", paste(names[!names %in% "y.train"], collapse = " + ")))
    #neural.train <- neuralnet(formula,data=train,hidden=1,linear.output=F)
    #test <- data.frame(x.test)
    #neural.pred <- compute(neural.train, test)$net.result
    #LDA
    lda.train <- lda(x.train,y.train,tol=0)
    lda.pred <- c(lda.pred,predict(lda.train,x.test)$posterior[,2])
  }
  y.sample <- as.vector(subset)
  y.sample <- y.sample[which(y.sample>0)]
  y.sample.matrix[,n] <- y.sample
  y.test.matrix[,n] <- y.test
  rf.pred.matrix[,n] <- rf.pred
  xgboost.pred.matrix[,n] <- xgboost.pred
  svm.pred.matrix[,n] <- svm.pred
  lasso.pred.matrix[,n] <- lasso.pred
  ridge.pred.matrix[,n] <- ridge.pred
  enet.pred.matrix[,n] <- enet.pred
  knn.pred.matrix[,n] <- knn.pred
  #neural.pred.matrix[,n] <- neural.pred
  lda.pred.matrix[,n] <- lda.pred
  #x.sample <- x1[y.sample,]
  #hierarchial clustering
  #hclust.pred.matrix[,n] <- cutree(hclust(dist(x.sample)),k=2)
  #k-means clustering
  #kmeans.pred.matrix[,n] <- kmeans(x.sample,centers=2)$cluster
}


save(y.sample.matrix,y.test.matrix,rf.pred.matrix,xgboost.pred.matrix,
     svm.pred.matrix,lasso.pred.matrix,ridge.pred.matrix,enet.pred.matrix,
     knn.pred.matrix,#neural.pred.matrix,
     lda.pred.matrix,#hclust.pred.matrix,kmeans.pred.matrix,
     file="CMS_code/infiles/Caligus_microbiome_predict_ASV.Rdata")


load("CMS_code/infiles/Caligus_microbiome_predict_ASV.Rdata")
#calculate average predicted y by sample, then generate ROC between predicted and actual y
for (n in 1:numsets) {
  y.test.temp <- cbind(y.sample.matrix[,n],y.test.matrix[,n])
  y.test.temp <- y.test.temp[order(y.test.temp[,1]),]
  y.test.matrix[,n] <- y.test.temp[,2]
  rf.pred.temp <- cbind(y.sample.matrix[,n],rf.pred.matrix[,n])
  rf.pred.temp <- rf.pred.temp[order(rf.pred.temp[,1]),]
  rf.pred.matrix[,n] <- rf.pred.temp[,2]
  xgboost.pred.temp <- cbind(y.sample.matrix[,n],xgboost.pred.matrix[,n])
  xgboost.pred.temp <- xgboost.pred.temp[order(xgboost.pred.temp[,1]),]
  xgboost.pred.matrix[,n] <- xgboost.pred.temp[,2]
  svm.pred.temp <- cbind(y.sample.matrix[,n],svm.pred.matrix[,n])
  svm.pred.temp <- svm.pred.temp[order(svm.pred.temp[,1]),]
  svm.pred.matrix[,n] <- svm.pred.temp[,2]
  lasso.pred.temp <- cbind(y.sample.matrix[,n],lasso.pred.matrix[,n])
  lasso.pred.temp <- lasso.pred.temp[order(lasso.pred.temp[,1]),]
  lasso.pred.matrix[,n] <- lasso.pred.temp[,2]
  ridge.pred.temp <- cbind(y.sample.matrix[,n],ridge.pred.matrix[,n])
  ridge.pred.temp <- ridge.pred.temp[order(ridge.pred.temp[,1]),]
  ridge.pred.matrix[,n] <- ridge.pred.temp[,2]
  enet.pred.temp <- cbind(y.sample.matrix[,n],enet.pred.matrix[,n])
  enet.pred.temp <- enet.pred.temp[order(enet.pred.temp[,1]),]
  enet.pred.matrix[,n] <- enet.pred.temp[,2]
  knn.pred.temp <- cbind(y.sample.matrix[,n],knn.pred.matrix[,n])
  knn.pred.temp <- knn.pred.temp[order(knn.pred.temp[,1]),]
  knn.pred.matrix[,n] <- knn.pred.temp[,2]
  #neural.pred.temp <- cbind(y.sample.matrix[,n],neural.pred.matrix[,n])
  #neural.pred.temp <- neural.pred.temp[order(neural.pred.temp[,1]),]
  #neural.pred.matrix[,n] <- neural.pred.temp[,2]
  lda.pred.temp <- cbind(y.sample.matrix[,n],lda.pred.matrix[,n])
  lda.pred.temp <- lda.pred.temp[order(lda.pred.temp[,1]),]
  lda.pred.matrix[,n] <- lda.pred.temp[,2]
  #hclust.pred.temp <- cbind(y.sample.matrix[,n],hclust.pred.matrix[,n])
  #hclust.pred.temp <- hclust.pred.temp[order(hclust.pred.temp[,1]),]
  #hclust.pred.matrix[,n] <- hclust.pred.temp[,2]
  #kmeans.pred.temp <- cbind(y.sample.matrix[,n],kmeans.pred.matrix[,n])
  #kmeans.pred.temp <- kmeans.pred.temp[order(kmeans.pred.temp[,1]),]
  #kmeans.pred.matrix[,n] <- kmeans.pred.temp[,2]
}

identical(y,y.test.matrix[,1])
rf.pred.avg <- rowMeans(rf.pred.matrix)
rf <- roc(y,rf.pred.avg)
rf$auc
xgboost.pred.avg <- rowMeans(xgboost.pred.matrix)
xgboost <- roc(y,xgboost.pred.avg)
xgboost$auc
svm.pred.avg <- rowMeans(svm.pred.matrix)
svm <- roc(y,svm.pred.avg)
svm$auc
lasso.pred.avg <- rowMeans(lasso.pred.matrix)
lasso <- roc(y,lasso.pred.avg)
lasso$auc
ridge.pred.avg <- rowMeans(ridge.pred.matrix)
ridge <- roc(y,ridge.pred.avg)
ridge$auc
enet.pred.avg <- rowMeans(enet.pred.matrix)
enet <- roc(y,enet.pred.avg)
enet$auc
knn.pred.avg <- rowMeans(knn.pred.matrix)
knn <- roc(y,knn.pred.avg)
knn$auc
#neural.pred.avg <- rowMeans(neural.pred.matrix)
#neural <- roc(y,neural.pred.avg)
#neural$auc
lda.pred.avg <- rowMeans(lda.pred.matrix)
lda <- roc(y,lda.pred.avg)
lda$auc
#hclust.pred.avg <- rowMeans(hclust.pred.matrix)
#hclust <- roc(y,hclust.pred.avg)
#hclust$auc
#kmeans.pred.avg <- rowMeans(kmeans.pred.matrix)
#kmeans <- roc(y,kmeans.pred.avg)
#kmeans$auc

#combine AUCs into a vector for future plots, etc.
#(e.g. AUCs between original data and HFE features, 500 most abundant ASV data, etc.)
original <- c(rf$auc,xgboost$auc,svm$auc,lasso$auc,ridge$auc,enet$auc,knn$auc,
              #neural$auc,
              lda$auc)
original
#0.5218254 0.4503968 0.6755952 0.6944444 0.5758929 0.6299603 0.5687004 0.5528274


########################################################################################
#### Set up and run 5-fold cross validation for 8 ML classification models
#### (ASV count matrix transformation by total-sum scaling and mean-centering)
########################################################################################

#otu_table
otu_table <- as.data.frame(High_Low_microbiome@otu_table)
#check that order of samples in metadata and OTU table are identical
colnames(otu_table) <- gsub("[.]","-",colnames(otu_table))
otu_table <- otu_table[,rownames(metadata_new)]
#Selecting 500 more abundant ASVs
#Rank ASVs by total abundance
# Calculate the total abundance for each ASV
#total_abundance <- rowSums(otu_table)
# Select the 500 most abundant ASVs
# Order the ASVs by total abundance in descending order
#otu_table <- otu_table[order(total_abundance, decreasing = TRUE), ]
# Select the top 1,000 ASVs
#top_500_asvs <- otu_table[1:500, ]

x <- data.matrix(otu_table)
#calculate relative abundance of each OTU by dividing its value by the total reads per sample
x <- sweep(x,2,colSums(x),"/")
x <- t(x)
#Center the matrix
x1<-scale(x,center=TRUE,scale=F)
id <- seq(1:length(y))
numsets <- 100


#actual y representing 100 random samples of sample IDs
y.sample.matrix <- matrix(nrow=length(y),ncol=numsets)
y.test.matrix <- matrix(nrow=length(y),ncol=numsets)
#predictions are probabilitities that y=1 for each sample
rf.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
xgboost.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
svm.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
lasso.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
ridge.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
enet.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
knn.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
#neural.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
lda.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
#hclust.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
#kmeans.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
set.seed(1000)


for (n in 1:numsets) {
  print(n)
  #randomly split samples into five groups for 5-fold cross-validation
  #recall that there are 90 (not 70) samples in this version of the data set
  subset <- matrix(0,nrow=18,ncol=5)
  subset[1:18,1] <- sort(sample(id,18))
  subset[1:18,2] <- sort(sample(which(!(id %in% subset)),18))
  subset[1:18,3] <- sort(sample(which(!(id %in% subset)),18))
  subset[1:18,4] <- sort(sample(which(!(id %in% subset)),18))
  subset[1:18,5] <- which(!(id %in% subset))
  y.test <- NULL
  rf.pred <- NULL
  xgboost.pred <- NULL
  svm.pred <- NULL
  lasso.pred <- NULL
  ridge.pred <- NULL
  enet.pred <- NULL
  knn.pred <- NULL
  #neural.pred <- NULL
  lda.pred <- NULL
  #use 4 groups as training set and 1 group as testing set for each supervised method, and repeat until each group is tested
  for (i in 1:5) {
    print(i)
    y.test <- c(y.test, y[subset[,i]])
    y.train <- y[-subset[,i]]
    x.train <- x1[-subset[,i],]
    #remove OTUs that have all zeroes in training set
    x.train <- subset(x.train,select=apply(x.train,2,var)>0)
    x.test <- x1[subset[,i],colnames(x.train)]
    #random forest
    rf.train <- randomForest(x.train,y=as.factor(y.train), importance = TRUE)
    rf.pred <- c(rf.pred,predict(rf.train,x.test,type="prob")[,2])
    #gradient boosting
    xgboost.train <- xgboost(x.train,label=y.train,nrounds=10,objective="binary:logistic",verbose=0)
    xgboost.pred <- c(xgboost.pred,predict(xgboost.train,x.test))
    #SVM
    svm.train <- svm(x.train,as.factor(y.train),probability=T)
    svm.model <- predict(svm.train,x.test,probability=T)
    svm.pred <- c(svm.pred,attr(svm.model,"probabilities")[,1])
    #lasso
    lasso.train <- train(x.train,as.factor(y.train),method="glmnet",tuneGrid=expand.grid(.alpha=1,.lambda=seq(0,1,by=0.1)))
    lasso.pred <- c(lasso.pred,predict(lasso.train,x.test,type="prob")[,2])
    #ridge
    ridge.train <- train(x.train,as.factor(y.train),method="glmnet",tuneGrid=expand.grid(.alpha=0,.lambda=seq(0,1,by=0.1)))
    ridge.pred <- c(ridge.pred,predict(ridge.train,x.test,type="prob")[,2])
    #elastic net
    enet.train <- train(x.train,as.factor(y.train),method="glmnet")
    enet.pred <- c(enet.pred,predict(enet.train,x.test,type="prob")[,2])
    #k-nearest neighbors
    knn.train <- knn3(x.train,as.factor(y.train))
    knn.pred <- c(knn.pred,predict(knn.train,x.test,type="prob")[,2])
    #neural networks
    #train <- data.frame(y.train,x.train)
    #names <- names(train)
    #formula <- as.formula(paste("y.train ~", paste(names[!names %in% "y.train"], collapse = " + ")))
    #neural.train <- neuralnet(formula,data=train,hidden=1,linear.output=F)
    #test <- data.frame(x.test)
    #neural.pred <- compute(neural.train, test)$net.result
    #LDA
    lda.train <- lda(x.train,y.train,tol=0)
    lda.pred <- c(lda.pred,predict(lda.train,x.test)$posterior[,2])
  }
  y.sample <- as.vector(subset)
  y.sample <- y.sample[which(y.sample>0)]
  y.sample.matrix[,n] <- y.sample
  y.test.matrix[,n] <- y.test
  rf.pred.matrix[,n] <- rf.pred
  xgboost.pred.matrix[,n] <- xgboost.pred
  svm.pred.matrix[,n] <- svm.pred
  lasso.pred.matrix[,n] <- lasso.pred
  ridge.pred.matrix[,n] <- ridge.pred
  enet.pred.matrix[,n] <- enet.pred
  knn.pred.matrix[,n] <- knn.pred
  #neural.pred.matrix[,n] <- neural.pred
  lda.pred.matrix[,n] <- lda.pred
  #x.sample <- x1[y.sample,]
  #hierarchial clustering
  #hclust.pred.matrix[,n] <- cutree(hclust(dist(x.sample)),k=2)
  #k-means clustering
  #kmeans.pred.matrix[,n] <- kmeans(x.sample,centers=2)$cluster
}


save(y.sample.matrix,y.test.matrix,rf.pred.matrix,xgboost.pred.matrix,
     svm.pred.matrix,lasso.pred.matrix,ridge.pred.matrix,enet.pred.matrix,
     knn.pred.matrix,#neural.pred.matrix,
     lda.pred.matrix,#hclust.pred.matrix,kmeans.pred.matrix,
     file="CMS_code/infiles/Caligus_microbiome_predict_ASV_Cent.Rdata")


load("CMS_code/infiles/Caligus_microbiome_predict_ASV_Cent.Rdata")
#calculate average predicted y by sample, then generate ROC between predicted and actual y
for (n in 1:numsets) {
  y.test.temp <- cbind(y.sample.matrix[,n],y.test.matrix[,n])
  y.test.temp <- y.test.temp[order(y.test.temp[,1]),]
  y.test.matrix[,n] <- y.test.temp[,2]
  rf.pred.temp <- cbind(y.sample.matrix[,n],rf.pred.matrix[,n])
  rf.pred.temp <- rf.pred.temp[order(rf.pred.temp[,1]),]
  rf.pred.matrix[,n] <- rf.pred.temp[,2]
  xgboost.pred.temp <- cbind(y.sample.matrix[,n],xgboost.pred.matrix[,n])
  xgboost.pred.temp <- xgboost.pred.temp[order(xgboost.pred.temp[,1]),]
  xgboost.pred.matrix[,n] <- xgboost.pred.temp[,2]
  svm.pred.temp <- cbind(y.sample.matrix[,n],svm.pred.matrix[,n])
  svm.pred.temp <- svm.pred.temp[order(svm.pred.temp[,1]),]
  svm.pred.matrix[,n] <- svm.pred.temp[,2]
  lasso.pred.temp <- cbind(y.sample.matrix[,n],lasso.pred.matrix[,n])
  lasso.pred.temp <- lasso.pred.temp[order(lasso.pred.temp[,1]),]
  lasso.pred.matrix[,n] <- lasso.pred.temp[,2]
  ridge.pred.temp <- cbind(y.sample.matrix[,n],ridge.pred.matrix[,n])
  ridge.pred.temp <- ridge.pred.temp[order(ridge.pred.temp[,1]),]
  ridge.pred.matrix[,n] <- ridge.pred.temp[,2]
  enet.pred.temp <- cbind(y.sample.matrix[,n],enet.pred.matrix[,n])
  enet.pred.temp <- enet.pred.temp[order(enet.pred.temp[,1]),]
  enet.pred.matrix[,n] <- enet.pred.temp[,2]
  knn.pred.temp <- cbind(y.sample.matrix[,n],knn.pred.matrix[,n])
  knn.pred.temp <- knn.pred.temp[order(knn.pred.temp[,1]),]
  knn.pred.matrix[,n] <- knn.pred.temp[,2]
  #neural.pred.temp <- cbind(y.sample.matrix[,n],neural.pred.matrix[,n])
  #neural.pred.temp <- neural.pred.temp[order(neural.pred.temp[,1]),]
  #neural.pred.matrix[,n] <- neural.pred.temp[,2]
  lda.pred.temp <- cbind(y.sample.matrix[,n],lda.pred.matrix[,n])
  lda.pred.temp <- lda.pred.temp[order(lda.pred.temp[,1]),]
  lda.pred.matrix[,n] <- lda.pred.temp[,2]
  #hclust.pred.temp <- cbind(y.sample.matrix[,n],hclust.pred.matrix[,n])
  #hclust.pred.temp <- hclust.pred.temp[order(hclust.pred.temp[,1]),]
  #hclust.pred.matrix[,n] <- hclust.pred.temp[,2]
  #kmeans.pred.temp <- cbind(y.sample.matrix[,n],kmeans.pred.matrix[,n])
  #kmeans.pred.temp <- kmeans.pred.temp[order(kmeans.pred.temp[,1]),]
  #kmeans.pred.matrix[,n] <- kmeans.pred.temp[,2]
}

identical(y,y.test.matrix[,1])
rf.pred.avg <- rowMeans(rf.pred.matrix)
rf <- roc(y,rf.pred.avg)
rf$auc
xgboost.pred.avg <- rowMeans(xgboost.pred.matrix)
xgboost <- roc(y,xgboost.pred.avg)
xgboost$auc
svm.pred.avg <- rowMeans(svm.pred.matrix)
svm <- roc(y,svm.pred.avg)
svm$auc
lasso.pred.avg <- rowMeans(lasso.pred.matrix)
lasso <- roc(y,lasso.pred.avg)
lasso$auc
ridge.pred.avg <- rowMeans(ridge.pred.matrix)
ridge <- roc(y,ridge.pred.avg)
ridge$auc
enet.pred.avg <- rowMeans(enet.pred.matrix)
enet <- roc(y,enet.pred.avg)
enet$auc
knn.pred.avg <- rowMeans(knn.pred.matrix)
knn <- roc(y,knn.pred.avg)
knn$auc
#neural.pred.avg <- rowMeans(neural.pred.matrix)
#neural <- roc(y,neural.pred.avg)
#neural$auc
lda.pred.avg <- rowMeans(lda.pred.matrix)
lda <- roc(y,lda.pred.avg)
lda$auc
#hclust.pred.avg <- rowMeans(hclust.pred.matrix)
#hclust <- roc(y,hclust.pred.avg)
#hclust$auc
#kmeans.pred.avg <- rowMeans(kmeans.pred.matrix)
#kmeans <- roc(y,kmeans.pred.avg)
#kmeans$auc

#combine AUCs into a vector for future plots, etc.
#(e.g. AUCs between original data and HFE features, 500 most abundant ASV data, etc.)
original_Cent <- c(rf$auc,xgboost$auc,svm$auc,lasso$auc,ridge$auc,enet$auc,knn$auc,
                    #neural$auc,
                    lda$auc)
original_Cent
#0.5143849 0.4459325 0.6860119 0.6731151 0.5887897 0.6205357 0.5691964 0.5820933

plot(rf$specificities,rf$sensitivities,type="n",xlim=c(1,0),xlab="Specificity",ylab="Sensitivity")
lines(rf$specificities,rf$sensitivities,col="black")
par(new=T)
plot(xgboost$specificities,xgboost$sensitivities,type="n",xlim=c(1,0),xlab=NA,ylab=NA,axes=F)
lines(xgboost$specificities,xgboost$sensitivities,col="blue")
par(new=T)
plot(svm$specificities,svm$sensitivities,type="n",xlim=c(1,0),xlab=NA,ylab=NA,axes=F)
lines(svm$specificities,svm$sensitivities,col="red")
par(new=T)
plot(lasso$specificities,lasso$sensitivities,type="n",xlim=c(1,0),xlab=NA,ylab=NA,axes=F)
lines(lasso$specificities,lasso$sensitivities,col="green")
par(new=T)
plot(ridge$specificities,ridge$sensitivities,type="n",xlim=c(1,0),xlab=NA,ylab=NA,axes=F)
lines(ridge$specificities,ridge$sensitivities,col="orange")
par(new=T)
plot(enet$specificities,enet$sensitivities,type="n",xlim=c(1,0),xlab=NA,ylab=NA,axes=F)
lines(enet$specificities,enet$sensitivities,col="purple")
par(new=T)
plot(knn$specificities,knn$sensitivities,type="n",xlim=c(1,0),xlab=NA,ylab=NA,axes=F)
lines(knn$specificities,knn$sensitivities,col="brown")
par(new=T)
plot(lda$specificities,lda$sensitivities,type="n",xlim=c(1,0),xlab=NA,ylab=NA,axes=F)
lines(lda$specificities,lda$sensitivities,col="gold")
abline(1,-1,col="gray")
legend("bottomright",c("RF (0.51)","Gboost (0.45)","SVM (0.69)","Lasso (0.67)","Ridge (0.59)","Enet (0.62)","k-NN (0.57)","LDA (0.58)"),
       col=c("black","blue","red","green","orange","purple","brown","gold"),
       lty=1,
       cex=0.6)


#########################################################################################
#### Prepare ASV table for Hierarchical Feature Engineering (HFE) analysis
#### Then set up and run 5-fold cross validation for 8 ML classification models as before
#### (Count matrix transformation by total-sum scaling and mean-centering)
#########################################################################################

###
#HFE analysis
###

##OTU TABLE for HFE

# Transform the counts to relative abundances
#High_Low_microbiome_relative <- transform_sample_counts(High_Low_microbiome, function(x) x / sum(x))

#physeq_aggregated <- tax_glom(High_Low_microbiome, taxrank = "Genus")

tax_table <- as.data.frame(High_Low_microbiome@tax_table)
tax_table$Species <- row.names(tax_table)

# Function to combine taxonomic ranks, excluding NA values
combine_taxonomic_ranks <- function(row) {
  # Define the taxonomic ranks and their prefixes
  tax_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  prefixes <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")
  
  # Combine ranks and prefixes, excluding NA values
  combined <- paste0(prefixes, row[!is.na(row)])
  
  # Collapse the combined parts into a single string
  taxonomy <- paste(combined, collapse = "|")
  return(taxonomy)
}

# Apply the function to each row of the tax_table
tax_table$Taxonomy <- apply(tax_table, 1, combine_taxonomic_ranks)

#input OTU table file for HFE
#OTU TABLE
otu_table <- data.frame(tax_table[,8],
                        as.data.frame(High_Low_microbiome@otu_table))

#calculate relative abundance of each OTU by dividing its value by the total reads per sample
otu_table <- data.frame(otu_table[,1],
                        sweep(otu_table[,-1],2,colSums(otu_table[,-1]),"/"))
colnames(otu_table)[1] <- "clade_name"

write.table(otu_table, "CMS_code/infiles/HFE_otu_caligus.txt",
            quote=FALSE, row.names = F, col.names = T, sep = "\t")


#ps.relt_selected = subset_samples(ps_reordered, Status != "Neutral")

#Creating metadata
metadata_hfe <- High_Low_microbiome@sam_data
str(metadata_hfe)
metadata_hfe$Burden <- as.factor(metadata_hfe$Burden)

metadata_hfe <- metadata_hfe[,c(1,17)]
write.table(data.frame(metadata_hfe),
            "CMS_code/infiles/HFE_metadata_caligus.txt",
            quote=FALSE, row.names = F, sep = "\t")

###############################################################################################################
# Run HFE here locally from terminal:
# $ docker run --rm -it -v `pwd`:/home/docker -w /home/docker aoliver44/taxa_hfe:latest bash
# $ taxaHFE --subject_identifier Sample_name --label Burden --lowest_level 3 --ncores 8 --seed 123 \
#           ./HFE_metadata_caligus.txt ./HFE_otu_caligus.txt ./HFE_result_output.csv
###############################################################################################################

## Now, use HFE output for ML steps

library(dplyr)
library(microViz)
HFE_result <- read.csv("~/Clay_analysis/HFE_result_output.csv")
#Creating metadata
metadata <- High_Low_microbiome@sam_data
str(metadata)
metadata$Burden <- as.factor(metadata$Burden)

otu_HFE <- HFE_result[,-c(1,2)]
rownames(otu_HFE) <- metadata$Sample_name

y <- as.numeric(metadata$Burden =="High")
x <- data.matrix(otu_HFE)
x <- as.matrix(x)

#Center the matrix
x1<-scale(x,center=TRUE,scale=F)
id <- seq(1:length(y))
numsets <- 100


#actual y representing 100 random samples of sample IDs
y.sample.matrix <- matrix(nrow=length(y),ncol=numsets)
y.test.matrix <- matrix(nrow=length(y),ncol=numsets)
#predictions are probabilitities that y=1 for each sample
rf.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
xgboost.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
svm.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
lasso.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
ridge.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
enet.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
knn.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
#neural.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
lda.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
#hclust.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
#kmeans.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
set.seed(1000)
for (n in 1:numsets) {
  print(n)
  #randomly split samples into five groups for 5-fold cross-validation
  #(recall, there are 90 samples for this data set)
  subset <- matrix(0,nrow=18,ncol=5)
  subset[1:18,1] <- sort(sample(id,18))
  subset[1:18,2] <- sort(sample(which(!(id %in% subset)),18))
  subset[1:18,3] <- sort(sample(which(!(id %in% subset)),18))
  subset[1:18,4] <- sort(sample(which(!(id %in% subset)),18))
  subset[1:18,5] <- which(!(id %in% subset))
  y.test <- NULL
  rf.pred <- NULL
  xgboost.pred <- NULL
  svm.pred <- NULL
  lasso.pred <- NULL
  ridge.pred <- NULL
  enet.pred <- NULL
  knn.pred <- NULL
  #neural.pred <- NULL
  lda.pred <- NULL
  #use 4 groups as training set and 1 group as testing set for each supervised method, and repeat until each group is tested
  for (i in 1:5) {
    print(i)
    y.test <- c(y.test, y[subset[,i]])
    y.train <- y[-subset[,i]]
    #remove OTUs that have all zeroes in training set
    x.train <- x1[-subset[,i],]
    x.train <- subset(x.train,select=apply(x.train,2,var)>0)
    x.test <- x1[subset[,i],colnames(x.train)]
    #random forest
    rf.train <- randomForest(x.train,y=as.factor(y.train), importance = TRUE)
    rf.pred <- c(rf.pred,predict(rf.train,x.test,type="prob")[,2])
    #gradient boosting
    xgboost.train <- xgboost(x.train,label=y.train,nrounds=10,objective="binary:logistic",verbose=0)
    xgboost.pred <- c(xgboost.pred,predict(xgboost.train,x.test))
    #SVM
    svm.train <- svm(x.train,as.factor(y.train),probability=T)
    svm.model <- predict(svm.train,x.test,probability=T)
    svm.pred <- c(svm.pred,attr(svm.model,"probabilities")[,1])
    #lasso
    lasso.train <- train(x.train,as.factor(y.train),method="glmnet",tuneGrid=expand.grid(.alpha=1,.lambda=seq(0,1,by=0.1)))
    lasso.pred <- c(lasso.pred,predict(lasso.train,x.test,type="prob")[,2])
    #ridge
    ridge.train <- train(x.train,as.factor(y.train),method="glmnet",tuneGrid=expand.grid(.alpha=0,.lambda=seq(0,1,by=0.1)))
    ridge.pred <- c(ridge.pred,predict(ridge.train,x.test,type="prob")[,2])
    #elastic net
    enet.train <- train(x.train,as.factor(y.train),method="glmnet")
    enet.pred <- c(enet.pred,predict(enet.train,x.test,type="prob")[,2])
    #k-nearest neighbors
    knn.train <- knn3(x.train,as.factor(y.train))
    knn.pred <- c(knn.pred,predict(knn.train,x.test,type="prob")[,2])
    #neural networks
    #train <- data.frame(y.train,x.train)
    #names <- names(train)
    #formula <- as.formula(paste("y.train ~", paste(names[!names %in% "y.train"], collapse = " + ")))
    #neural.train <- neuralnet(formula,data=train,hidden=1,linear.output=F)
    #neural.pred <- c(neural.pred,compute(neural.train,x.test)$net.result)
    #LDA
    lda.train <- lda(x.train,y.train,tol=0)
    lda.pred <- c(lda.pred,predict(lda.train,x.test)$posterior[,2])
  }
  y.sample <- as.vector(subset)
  y.sample <- y.sample[which(y.sample>0)]
  y.sample.matrix[,n] <- y.sample
  y.test.matrix[,n] <- y.test
  rf.pred.matrix[,n] <- rf.pred
  xgboost.pred.matrix[,n] <- xgboost.pred
  svm.pred.matrix[,n] <- svm.pred
  lasso.pred.matrix[,n] <- lasso.pred
  ridge.pred.matrix[,n] <- ridge.pred
  enet.pred.matrix[,n] <- enet.pred
  knn.pred.matrix[,n] <- knn.pred
  #neural.pred.matrix[,n] <- neural.pred
  lda.pred.matrix[,n] <- lda.pred
  x.sample <- x[y.sample,]
  #hierarchial clustering
  #hclust.pred.matrix[,n] <- cutree(hclust(dist(x.sample)),k=2)
  #k-means clustering
  #kmeans.pred.matrix[,n] <- kmeans(x.sample,centers=2)$cluster
}


save(y.sample.matrix,
     y.test.matrix,
     rf.pred.matrix,
     xgboost.pred.matrix,
     svm.pred.matrix,
     lasso.pred.matrix,
     ridge.pred.matrix,
     enet.pred.matrix,
     knn.pred.matrix,
     #neural.pred.matrix,
     lda.pred.matrix,
     #hclust.pred.matrix,
     #kmeans.pred.matrix,
     file="CMS_code/infiles/Caligus_microbiome_predict_HFE.Rdata")


#calculate average predicted y by sample
load("CMS_code/infiles/Caligus_microbiome_predict_HFE.Rdata")
for (n in 1:numsets) {
  y.test.temp <- cbind(y.sample.matrix[,n],y.test.matrix[,n])
  y.test.temp <- y.test.temp[order(y.test.temp[,1]),]
  y.test.matrix[,n] <- y.test.temp[,2]
  rf.pred.temp <- cbind(y.sample.matrix[,n],rf.pred.matrix[,n])
  rf.pred.temp <- rf.pred.temp[order(rf.pred.temp[,1]),]
  rf.pred.matrix[,n] <- rf.pred.temp[,2]
  xgboost.pred.temp <- cbind(y.sample.matrix[,n],xgboost.pred.matrix[,n])
  xgboost.pred.temp <- xgboost.pred.temp[order(xgboost.pred.temp[,1]),]
  xgboost.pred.matrix[,n] <- xgboost.pred.temp[,2]
  svm.pred.temp <- cbind(y.sample.matrix[,n],svm.pred.matrix[,n])
  svm.pred.temp <- svm.pred.temp[order(svm.pred.temp[,1]),]
  svm.pred.matrix[,n] <- svm.pred.temp[,2]
  lasso.pred.temp <- cbind(y.sample.matrix[,n],lasso.pred.matrix[,n])
  lasso.pred.temp <- lasso.pred.temp[order(lasso.pred.temp[,1]),]
  lasso.pred.matrix[,n] <- lasso.pred.temp[,2]
  ridge.pred.temp <- cbind(y.sample.matrix[,n],ridge.pred.matrix[,n])
  ridge.pred.temp <- ridge.pred.temp[order(ridge.pred.temp[,1]),]
  ridge.pred.matrix[,n] <- ridge.pred.temp[,2]
  enet.pred.temp <- cbind(y.sample.matrix[,n],enet.pred.matrix[,n])
  enet.pred.temp <- enet.pred.temp[order(enet.pred.temp[,1]),]
  enet.pred.matrix[,n] <- enet.pred.temp[,2]
  knn.pred.temp <- cbind(y.sample.matrix[,n],knn.pred.matrix[,n])
  knn.pred.temp <- knn.pred.temp[order(knn.pred.temp[,1]),]
  knn.pred.matrix[,n] <- knn.pred.temp[,2]
  #neural.pred.temp <- cbind(y.sample.matrix[,n],neural.pred.matrix[,n])
  #neural.pred.temp <- neural.pred.temp[order(neural.pred.temp[,1]),]
  #neural.pred.matrix[,n] <- neural.pred.temp[,2]
  lda.pred.temp <- cbind(y.sample.matrix[,n],lda.pred.matrix[,n])
  lda.pred.temp <- lda.pred.temp[order(lda.pred.temp[,1]),]
  lda.pred.matrix[,n] <- lda.pred.temp[,2]
  #hclust.pred.temp <- cbind(y.sample.matrix[,n],hclust.pred.matrix[,n])
  #hclust.pred.temp <- hclust.pred.temp[order(hclust.pred.temp[,1]),]
  #hclust.pred.matrix[,n] <- hclust.pred.temp[,2]
  #kmeans.pred.temp <- cbind(y.sample.matrix[,n],kmeans.pred.matrix[,n])
  #kmeans.pred.temp <- kmeans.pred.temp[order(kmeans.pred.temp[,1]),]
  #kmeans.pred.matrix[,n] <- kmeans.pred.temp[,2]
}

identical(y,y.test.matrix[,1])
rf.pred.avg <- rowMeans(rf.pred.matrix)
rf <- roc(y,rf.pred.avg)
rf$auc
xgboost.pred.avg <- rowMeans(xgboost.pred.matrix)
xgboost <- roc(y,xgboost.pred.avg)
xgboost$auc
svm.pred.avg <- rowMeans(svm.pred.matrix)
svm <- roc(y,svm.pred.avg)
svm$auc
lasso.pred.avg <- rowMeans(lasso.pred.matrix)
lasso <- roc(y,lasso.pred.avg)
lasso$auc
ridge.pred.avg <- rowMeans(ridge.pred.matrix)
ridge <- roc(y,ridge.pred.avg)
ridge$auc
enet.pred.avg <- rowMeans(enet.pred.matrix)
enet <- roc(y,enet.pred.avg)
enet$auc
knn.pred.avg <- rowMeans(knn.pred.matrix)
knn <- roc(y,knn.pred.avg)
knn$auc
#neural.pred.avg <- rowMeans(neural.pred.matrix)
#neural <- roc(y,neural.pred.avg)
#neural$auc
lda.pred.avg <- rowMeans(lda.pred.matrix)
lda <- roc(y,lda.pred.avg)
lda$auc
#hclust.pred.avg <- rowMeans(hclust.pred.matrix)
#hclust <- roc(y,hclust.pred.avg)
#hclust$auc
#kmeans.pred.avg <- rowMeans(kmeans.pred.matrix)
#kmeans <- roc(y,kmeans.pred.avg)
#kmeans$auc

hfe <- c(rf$auc,xgboost$auc,svm$auc,lasso$auc,ridge$auc,enet$auc,knn$auc,
         #neural$auc,
         lda$auc)
hfe
#0.7750496 0.7430556 0.6413690 0.6364087 0.6463294 0.6473214 0.6168155 0.6378968

plot(rf$specificities,rf$sensitivities,type="n",xlim=c(1,0),xlab="Specificity",ylab="Sensitivity")
lines(rf$specificities,rf$sensitivities,col="black")
par(new=T)
plot(xgboost$specificities,xgboost$sensitivities,type="n",xlim=c(1,0),xlab=NA,ylab=NA,axes=F)
lines(xgboost$specificities,xgboost$sensitivities,col="blue")
par(new=T)
plot(svm$specificities,svm$sensitivities,type="n",xlim=c(1,0),xlab=NA,ylab=NA,axes=F)
lines(svm$specificities,svm$sensitivities,col="red")
par(new=T)
plot(lasso$specificities,lasso$sensitivities,type="n",xlim=c(1,0),xlab=NA,ylab=NA,axes=F)
lines(lasso$specificities,lasso$sensitivities,col="green")
par(new=T)
plot(ridge$specificities,ridge$sensitivities,type="n",xlim=c(1,0),xlab=NA,ylab=NA,axes=F)
lines(ridge$specificities,ridge$sensitivities,col="orange")
par(new=T)
plot(enet$specificities,enet$sensitivities,type="n",xlim=c(1,0),xlab=NA,ylab=NA,axes=F)
lines(enet$specificities,enet$sensitivities,col="purple")
par(new=T)
plot(knn$specificities,knn$sensitivities,type="n",xlim=c(1,0),xlab=NA,ylab=NA,axes=F)
lines(knn$specificities,knn$sensitivities,col="brown")
par(new=T)
plot(lda$specificities,lda$sensitivities,type="n",xlim=c(1,0),xlab=NA,ylab=NA,axes=F)
lines(lda$specificities,lda$sensitivities,col="gold")
abline(1,-1,col="gray")
legend("bottomright",c("RF (0.78)","Gboost (0.74)","SVM (0.64)","Lasso (0.64)","Ridge (0.65)","Enet (0.65)","k-NN (0.62)","LDA (0.64)"),
       col=c("black","blue","red","green","orange","purple","brown","gold"),
       lty=1,
       cex=0.6)

## plot AUC between original data and informative features from HFE
## (original data not mean-centered; HFE data mean-centered)

plot(original[1],hfe[1],xlim=c(0.4,1),ylim=c(0.4,1),
     xlab="125 ASVs (whole dataset)",ylab="10 informative features from HFE",
     main="Microbiome Prediction for Sea lice burden in Atlantic salmon",
     pch=19,cex=1.5, cex.main=.7, col="black")
par(new=T)
plot(original[2],hfe[2],xlim=c(0.4,1),ylim=c(0.4,1),xlab=NA,ylab=NA,main=NA,pch=19,cex=1.5,col="blue")
par(new=T)
plot(original[3],hfe[3],xlim=c(0.4,1),ylim=c(0.4,1),xlab=NA,ylab=NA,main=NA,pch=19,cex=1.5,col="red")
par(new=T)
plot(original[4],hfe[4],xlim=c(0.4,1),ylim=c(0.4,1),xlab=NA,ylab=NA,main=NA,pch=19,cex=1.5,col="green")
par(new=T)
plot(original[5],hfe[5],xlim=c(0.4,1),ylim=c(0.4,1),xlab=NA,ylab=NA,main=NA,pch=19,cex=1.5,col="orange")
par(new=T)
plot(original[6],hfe[6],xlim=c(0.4,1),ylim=c(0.4,1),xlab=NA,ylab=NA,main=NA,pch=19,cex=1.5,col="purple")
par(new=T)
plot(original[7],hfe[7],xlim=c(0.4,1),ylim=c(0.4,1),xlab=NA,ylab=NA,main=NA,pch=19,cex=1.5,col="brown")
par(new=T)
plot(original[8],hfe[8],xlim=c(0.4,1),ylim=c(0.4,1),xlab=NA,ylab=NA,main=NA,pch=19,cex=1.5,col="gold")
abline(0,1,col="gray")
legend("bottomright",c("RF","Gboost","SVM","Lasso","Ridge","Enet","k-NN","LDA"),
       col=c("black","blue","red","green","orange","purple","brown","gold"),pch=19)



## plot AUC between original data and informative features from HFE
## (both original AND HFE data mean-centered)

plot(original_Cent[1],hfe[1],xlim=c(0.4,1),ylim=c(0.4,1),
     xlab="125 ASVs\n(whole dataset, mean-centered)",ylab="10 informative features from HFE",
     main="Microbiome Prediction for Sea lice burden in Atlantic salmon",
     pch=19,cex=1.5, cex.main=.7, col="black")
par(new=T)
plot(original_Cent[2],hfe[2],xlim=c(0.4,1),ylim=c(0.4,1),xlab=NA,ylab=NA,main=NA,pch=19,cex=1.5,col="blue")
par(new=T)
plot(original_Cent[3],hfe[3],xlim=c(0.4,1),ylim=c(0.4,1),xlab=NA,ylab=NA,main=NA,pch=19,cex=1.5,col="red")
par(new=T)
plot(original_Cent[4],hfe[4],xlim=c(0.4,1),ylim=c(0.4,1),xlab=NA,ylab=NA,main=NA,pch=19,cex=1.5,col="green")
par(new=T)
plot(original_Cent[5],hfe[5],xlim=c(0.4,1),ylim=c(0.4,1),xlab=NA,ylab=NA,main=NA,pch=19,cex=1.5,col="orange")
par(new=T)
plot(original_Cent[6],hfe[6],xlim=c(0.4,1),ylim=c(0.4,1),xlab=NA,ylab=NA,main=NA,pch=19,cex=1.5,col="purple")
par(new=T)
plot(original_Cent[7],hfe[7],xlim=c(0.4,1),ylim=c(0.4,1),xlab=NA,ylab=NA,main=NA,pch=19,cex=1.5,col="brown")
par(new=T)
plot(original_Cent[8],hfe[8],xlim=c(0.4,1),ylim=c(0.4,1),xlab=NA,ylab=NA,main=NA,pch=19,cex=1.5,col="gold")
abline(0,1,col="gray")
legend("bottomright",c("RF","Gboost","SVM","Lasso","Ridge","Enet","k-NN","LDA"),
       col=c("black","blue","red","green","orange","purple","brown","gold"),pch=19)

