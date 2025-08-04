#Chapter 3 
#2024 data
#Freebayes pipeline#
#Load packages####

analysis_pkgs <- c("tidyverse", "openxlsx", "furrr", "progressr", "tictoc", "scales",
                   "glue",'plotly', 'readxl',# general packages
                   "paletteer", "viridis",
                   "vcfR","adegenet", "seqinr", "poppr", "mmod", "hierfstat",
                   "treemap", "ape", "DevonDeRaad/SNPfiltR", # pop gen packages
                   'geneHapR','pegas', 'dplyr')

#pak::pak(analysis_pkgs)

pacman::p_load(char = basename(analysis_pkgs), update = FALSE, install = FALSE)

# load Ido's magic functions####
devtools::source_gist("7f63547158ecdbacf31b54a58af0d1cc", 
                      filename = "util.R")

#Read vcf file
#Filtered: no heterozygotes; quality > 30; bi-allelic polymorphic; SNPs only
vcf <- read.vcfR('data/oyster_hybrids_hac.clair3.snps.Q30.DP7.CR0.8.poly.recode.vcf.gz')

head(vcf)
vcf_samples <- colnames(vcf@gt)[-1]
#Exclude data reps
# #There is something wrong with 24830_c
# samples_to_exclude <- c("24438", "24735", "24753", "24229", "24216", "24583", "24574", "24587", "24841", "24830", "24806", "24830_c")
# 
# exclude_samples <- vcf_samples %in% samples_to_exclude
# selected_columns <- c('FORMAT',vcf_samples[!exclude_samples]) #Step4 choosing replicated samples (singular)
# 
# # selected_columns <- c('FORMAT', grep(sample_name,vcf_samples, value = TRUE))
# vcf <- vcf[, selected_columns]

# filter vcf by sample
filter_vcfr_samples <- function(vcf, include_samples, exclude_samples) {
  require(vcfR)
  rlang::check_exclusive(include_samples, exclude_samples)
  sample_names <- colnames(vcf@gt)[-1]
  if (missing(include_samples)) {
    include_samples <- sample_names[!sample_names %in% exclude_samples]
  }
  return(vcf[,c('FORMAT', include_samples)])
}

vcf <- filter_vcfr_samples(vcf, exclude_samples = exclude_samples)

#SNPfiltR filtering####
vcf

dp <- extract.gt(vcf, element = "DP", as.numeric=TRUE)
heatmap.bp(dp, rlabels = FALSE)

popmap <- tibble(
  id = vcf_samples,
  pop = sub("_.+$", "", vcf_samples)
)

#Visualize distributions
# hard_filter(vcfR = vcf)
vcf_filt <- hard_filter(vcfR=vcf, depth = 7, gq = 20) %>% 
            missing_by_sample(., cutoff = 0.3) %>% 
            missing_by_snp(., cutoff = 0.7) # %>% 
            # filter_biallelic(.) %>% 
            # min_mac(., min.mac = 1) 
            
  
# 2.16% of genotypes fall below a read depth of 5 and were converted to NA
# 6.47% of genotypes fall below a genotype quality of 30 and were converted to NA
# max_depth(vcf)
# vcf <- max_depth(vcf, maxdepth = 90)
# maxdepth cutoff is specified, filtered vcfR object will be returned
# 0.09% of SNPs were above a mean depth of 90 and were removed from the vcf
missing_by_sample(vcfR = vcf, popmap = popmap )
vcf<-missing_by_sample(vcfR=vcf, cutoff = .2)
#10 samples are above a 0.2 missing data cutoff, and were removed from VCF
#filter biallelic SNPs

vcf_sample_ids <- colnames(vcf_filt@gt)[-1]
filtered_popmap <- popmap[popmap$id %in% vcf_sample_ids, ]
# miss.snp <-assess_missing_data_pca(vcfR=vcf, popmap = filtered_popmap, thresholds = .95, clustering = FALSE)
# filtered_vcf <- missing_by_snp(vcf, cutoff = .95)
#cutoff is specified, filtered vcfR object will be returned
#9.21% removed at 5% missingness
#19.98% of SNPs fell below a completeness cutoff of 1 and were removed from the VCF. 0.03% missiningness
# vcfR::write.vcf(filtered_vcf, "output/2024_vcf_filtered_5%.vcf")
# 
# #Informative loci####
# filtered_vcf <- read.vcfR("output/2024_vcf_filtered_5%.vcf")

genind_vcf <- vcfR2genind(vcf_filt, return.alleles = TRUE)
# ploidy(informative_genind_vcf) <- 1
gen.vcf.inform <- informloci(genind_vcf, quiet = FALSE)
#2159 uninformative loci at 5%
inform_loci <- locNames(gen.vcf.inform)
# filtered_vcf <- vcf_filt
# filtered_vcf@fix 

keep_loci <- vcf_filt@fix %>% as_tibble() %>% 
  mutate(lociname = paste(sub(".", "_", CHROM, fixed = TRUE),POS,sep = "_")) %>% 
  pull(lociname) %in% inform_loci
table(keep_loci)

filtered_vcf <- vcf_filt[keep_loci,]

# dp <- extract.gt(filtered_vcf, element = "DP", as.numeric=TRUE)
# heatmap.bp(dp, rlabels = FALSE)
# 
vcfR::write.vcf(filtered_vcf, "output/Oyster_hybrids_filtered_vcf_informloci.vcf.gz")

#PCA
#Convert to a genind object
filtered_vcf <- read.vcfR("output/Oyster_hybrids_filtered_vcf_informloci.vcf.gz")
gen <- vcfR2genind(filtered_vcf)
allelfreq <- scaleGen(gen, NA.method="mean", scale=TRUE)

pca <- dudi.pca(allelfreq, scale = FALSE) # , scannf = FALSE, nf = 6)

#need the coordinated that each isolate exist on
# genind_df <- data.frame(id = indNames(gen))

pcadata <- pca$li %>% rownames_to_column('id') %>% 
           left_join(popmap)
pcavar <- pca$li %>% summarise(across(starts_with('Axis'), var))

percent_var <- pcavar/sum(pcavar)
# library(scales)
# library(paletteer)
# library(viridis)
pop_levels <- length(unique(pcadata$pop))
# sample_levels <- 
# palettes_d_names %>% as_tibble() %>% 
#   dplyr::filter(length>=pop_levels,type =="qualitative")
ggplot(pcadata, aes(x= Axis1, y= Axis2, colour = pop)) + geom_point(size = 4) +  scale_colour_paletteer_d("awtools::mpalette") + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  labs(colour = "Isolate", x= glue("C1: {percent(percent_var$Axis1)}"),
       y= glue("C2: {percent(percent_var$Axis2)}")) +
  theme_bw(14)
ggsave("output/Oyster_samples_PCA.pdf", width = 7, height = 6)

#Plot points labeled instead
library(ggrepel)
ggplot(pcadata, aes(x = Axis1, y = Axis2, colour = genotype, label = genotype)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3, max.overlaps = Inf) +  # Set max.overlaps to Inf to allow all labels
  scale_colour_viridis_d(option = 'viridis') +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  labs(colour = "Isolate", x = glue("C1: {percent(percent_var$Axis1)}"),
       y = glue("C2: {percent(percent_var$Axis2)}")) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")
ggsave("output/vcf1_PCA_exp_pop_trainin_infomloci_wlabs.pdf", width = 7, height = 6)

#MLH
#Preferable to use this code in chapter_4 code file. For this purpose it checks what might be wrong with technical replicates and if any need to be removed from the dataset.
#Define MLH to determine DAPC priori clusters
gen <-  vcfR2genind(filtered_vcf)
gen.ml_grid <- mlg.id(gen) 
gen.ml_grid 
gen <- poppr::as.genclone(gen)
#X <- genclone2genind(gen_obj1)
X <- genind2loci(gen)
#Pairwise genetic distances
pairwise_genetic_distance <- dist.gene(X, method="pairwise", pairwise.deletion = FALSE, variance = FALSE)
hist(pairwise_genetic_distance)
# Convert the matrix to a data frame
pairdist_df <- as.data.frame(as.matrix(pairwise_genetic_distance))
write.csv(pairdist_df, file = "output/2024_pairdist_matrix_bitwise.csv", row.names = TRUE)
#define threshold if needed
#for 12 clusters 0.0311778291
#In this case is to match technical rep 753, all the others are fine
(thresholds <- mlg.filter(gen, distance = "bitwise.dist", stats = "THRESHOLDS",
                          threshold = 1))
(pcut <- cutoff_predictor(thresholds))
#Create the threshold so it is just as many clusters shown by controls clustering together
#threshold of 0.0341085271 to match-up the neccessary data of 12 cluster's
mlg.filter(gen, distance = "bitwise.dist") <- 0.0311778291#pcut
mlg.table(gen)
gen.ml_grid <- mlg.id(gen) 


#DAPC####
 filtered_vcf <- read.vcfR("output/2024_filtered_vcf_inform_5%.vcf.gz")
# #There is something wrong with 24830_a, filter out
# vcf_samples <- colnames(filtered_vcf@gt)[-1]
# samples_to_exclude <- c("24830_a", "24229_c")
# exclude_samples <- vcf_samples %in% samples_to_exclude
# selected_columns <- c('FORMAT',vcf_samples[!exclude_samples]) 
# # 
# filtered_vcf <- filtered_vcf[, selected_columns]
# filtered_vcf <- filter_vcfr_samples(filtered_vcf, exclude_samples = exclude_samples)
# vcfR::write.vcf(filtered_vcf, "output/2024_filtered_vcf_inform_5%.vcf.gz")
#As genlight

#Filter technical reps out 
#Do this after DAPC has been established. Until then,comment this section out.
#filter out the control samples
#filter vcf by sample
exclude_samples <- c("AR0128", "AR0022", "AR0023","24841_b", "24841_c", "24438_b", "24753_a", "24753_b", "24830_a", "AR0210", "AR0020", "24806_c", "24806_b", "24229_c", "24583_b", "24583_c", "24216_b", "Ar0128") 
vcf_rep_exclude <- filter_vcfr_samples(filtered_vcf, exclude_samples = exclude_samples)
gen <- vcfR2genlight(vcf_rep_exclude)
gen <- vcfR2genlight(filtered_vcf)
gl <- gen
ploidy(gl) <- 1



#Set metadata
#duplicated_genotypes <- joined_data$genotype[duplicated(joined_data$genotype)]
#joined_data %>% 
  #filter(genotype %in% c("24180", "24220", "24766"))

metadata <- read_excel("output/chp_4/dapc11_updated_metadata.xlsx")
colnames(metadata)
filtered_metadata <- metadata[metadata$sample_ID %in% indNames(gl), ]
gl_df <- data.frame(genotype = indNames(gl))
joined_data <- gl_df %>%
  left_join(filtered_metadata, by = c("genotype" = "sample_ID"))

strata_data <- joined_data %>%
  dplyr::select(genotype, everything()) %>%
  column_to_rownames(var = "genotype")

strata(gl) <- strata_data

popNames(gl)
setPop(gl) <- ~host
gl

#Clustering exploration
set.seed(999)
grp <- find.clusters(gl)
#50
#16
names(grp)
head(grp$Kstat, 20)
grp$stat
head(grp$grp, 366)

#Clusteres by a-score
#Method 1
dapc_a <- dapc(gl, n.da = 100, n.pca = 20)
temp <- optim.a.score(dapc_a)
#7
#Method 2
set.seed(999)
grp <- find.clusters(tab(gl, NA.method = "mean"), max.n.clust = 20)
#50
#12
head(grp$Kstat, 20)
dapc_init <- dapc(tab(gl, NA.method = "mean"), grp$grp, n.pca = 20)
#5
opt <- optim.a.score(dapc_init)

dapc_final <- dapc(tab(gl, NA.method = "mean"), grp$grp, n.pca = opt$best)
scatter(dapc_final)


scatter(dapc_final,
        posi.da = "topright", # legend position
        bg = "white", # background color
        pch = 20, # point shape
        cstar = 1, # draw star-like connections
        col = rainbow(length(unique(grp$grp))), # color by group
        scree.da = TRUE)# show scree plot of LDs

#Clusteres y K-score
set.seed(999)
grp <- find.clusters(gl, max.n.clust=20)
#7
#11
#&7 and 12
head(grp$Kstat, 20)
grp$stat
head(grp$grp, 344)
grp$size
pop(gl)

dapc_final <- dapc(tab(gl, NA.method = "mean"), grp$grp, n.pca =  , n.da = 2)

chp4_freq <- table(adegenet::pop(gl), grp$grp)  #It is this table that you need for the frequency calculations in chapter 4, or as just a reference point for the mlg.id clustering using mlg.filter (it will correctly account for any missing data)
chp4_MLG_freq <- table(adegenet::pop(gl), grp$grp)

write.xlsx(chp4_MLG_freq, file = 'output/chp_4/2024_dapc_MLG_freq.11clust.xlsx')

table(adegenet::pop(gl))

dapc <- dapc(gl, grp$grp, n.pca = 10, n.da = 5)
#PC=K-1
#10
#5
scatter(dapc)
myCol<-c("darkblue","purple","darkgreen","orange","red","brown" ,"magenta" ,"gold","green", "lightblue", "pink", "grey", "darkred", "black" )
scatter(dapc, scree.da=TRUE, bg="white", pch=20, cell=0, cstar=0, col=myCol, solid=.4,cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:11))
scatter(dapc, scree.da = TRUE, bg = "white", pch = 20, cell = 0, cstar = 0, col = myCol, solid = 0.4, cex = 3, 
        clab = 0, leg = TRUE, txt.leg = paste("Cluster", 1:14), posi.leg = "bottomleft")  


#Niloofar's method
gl
indNames(gl)
myCol<-c("darkblue","purple","darkgreen","orange","red","brown" ,"magenta" ,"gold","green", "lightblue", "pink", "grey", "darkred", "black" )

gl.clusters <- find.clusters(gl, max.n.clust = 20, n.start =1000)
#7
#12
head(gl.clusters$Kstat, 25)
#K = where the variance begins to even out. 
dapc <- dapc(gl, gl.clusters$grp, var.contrib = TRUE)
# Choose the number PCs to retain (>=1): PC=K-1
#11
#5
scatter(dapc, legend = TRUE, cleg = 1, col = myCol, clabel = 0, cex = 1.5, scree.da = TRUE, scree.pca = TRUE, posi.pca = "topleft", posi.da = "bottomleft", ratio.da = 0.15, ratio.pca = 0.15, inset.da = 0.01, inset.pca = 0.01, posi.leg = "bottomright")
dapc$var
xlim <- par("usr")[1:2]
mtext(paste0("Variance explained: ", round(dapc$var, 4)),
      side = 3, line = -2, at = xlim[2] - 0.3 * diff(xlim), cex = 0.9)

#Making the DAPC look prettier

# Extract coordinates and cluster assignments
coords <- dapc$ind.coord
# Match with group assignments from dapc_group
groups <- dapc$grp

# Check for all clusters
table(dapc_data$Cluster)

# Assuming dapc_data is already created as:
# dapc_data <- data.frame(ID = rownames(coords),
#                         X = coords[,1],
#                         Y = coords[,2],
#                         Cluster = paste("Cluster", groups))

# Make publishible figure
# Define custom colors
myCol <- c("darkblue", "purple", "darkgreen", "orange", "red", "brown", "magenta", "gold","lightblue", "pink", "lightgreen" )


pop_map <- c("Cluster 1" = "MLG.8",  "Cluster 2" = "MLG.3",  "Cluster 3" = "MLG.11",  "Cluster 4" = "MLG.7", "Cluster 5" = "MLG.2",  "Cluster 6" = "MLG.9",  "Cluster 7" = "MLG.1",  "Cluster 8" = "MLG.10",  "Cluster 9" = "MLG.5",  "Cluster 10" = "MLG.6",  "Cluster 11" = "MLG.4")

library(ggrepel)
# Step 1: Extract isolate names from VCF
inds <- colnames(vcf_rep_exclude@gt)[-1]  # Skip first column

# Step 2: Normalize isolate names

isolate_table <- tibble(
  id = inds,
  isolate = id,
  control = grepl("^AR", id),
  comment = case_when(
    grepl("^AR", id) ~ "Control",
    TRUE ~ NA_character_
  )
)


# Step 3: Define control isolates (excluding samples if needed)
controls <- isolate_table %>%
  filter(control == TRUE) %>%
  pull(isolate) %>%
  unique()

# Step 4: Prepare plot data from DAPC
plot_data <- dapc$ind.coord %>%
  as.data.frame() %>%
  rownames_to_column("isolate") %>%
  mutate(is_control = isolate %in% controls)

# Step 5: Prepare control labels
control_labels <- plot_data %>%
  filter(is_control) %>%
  select(ID = isolate)

# Step 6: Normalize dapc_data$ID and join population info
dapc_data <- data.frame(ID = rownames(coords),
                        X = coords[,1],
                        Y = coords[,2],
                        Cluster = paste("Cluster", groups)) %>%
  mutate(ID = toupper(trimws(ID))) %>%
  left_join(pop_map_df, by = c("Cluster" = "isolate"))

# Step 7: Plot
ggplot(dapc_data, aes(x = X, y = Y, color = population)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_hline(yintercept = 0, color = "black", size = 0.5) +
  geom_vline(xintercept = 0, color = "black", size = 0.5) +
  geom_label_repel(
    data = control_labels %>% left_join(dapc_data, by = "ID"),
    aes(x = X, y = Y, label = ID),
    color = "black",
    max.overlaps = Inf,
    box.padding = 1,
    point.padding = 1,
    segment.size = 0.4,
    show.legend = FALSE
  ) +
  scale_color_manual(values = myCol) +
  labs(x = "LD1", y = "LD2", color = "Sample Type") +
  theme_minimal(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 13),
    plot.title = element_blank()
  )

write.csv(dapc_data, "output/chp_4/dapc_coordinates_2024_11_clusters_by_group_8.csv", row.names = FALSE)

# check which isolates belong to each of the clusters
dapc
dapc$grp
# Assuming dapc$grp is a factor or named vector
grp_df <- data.frame(
  isolate = names(dapc$grp),
  group = as.character(dapc$grp),
  stringsAsFactors = FALSE
)
# Split isolates by group
grouped <- split(grp_df$isolate, grp_df$group)
# Find the maximum group size to align columns
max_len <- max(sapply(grouped, length))
# Pad each group with NA to make equal-length columns
grouped_padded <- lapply(grouped, function(x) {
  length(x) <- max_len
  return(x)
})
# Combine into a data frame
group_table <- as.data.frame(grouped_padded, stringsAsFactors = FALSE)
# View the result
View(group_table)
write.xlsx(group_table, file = "output/chp_4/dapclust_11_group_table.xlsx", rowNames = FALSE)

#To add the clustering to the metadata for further analysis

df_cluster <- data.frame(
  sample_ID = names(dapc$grp),
  cluster = as.vector(dapc$grp)
)


metadata_updated <- metadata %>%
  left_join(df_cluster, by = "sample_ID") #id needing to write out new metadata with cluster infomation

write.xlsx(metadata_updated, file = "output/chp_4/dapc11_updated_metadata.xlsx", rowNames = TRUE)


strata_data <- strata_data %>%
  rownames_to_column(var = "sample_ID") %>%
  left_join(df_cluster, by = "sample_ID") %>%
  column_to_rownames(var = "sample_ID")
strata(gl) <- strata_data
gl

#Determine the variance within each cluster
# Assign populations based on clusters
setPop(gl) <- ~cluster

# Calculate variance captured within each cluster
# Function to calculate variance within each cluster

cluster_variance <- function(dapc_obj, genlight_obj) {
  clusters <- levels(pop(genlight_obj))
  variance_list <- list()
  
  for (cl in clusters) {
    inds <- which(pop(genlight_obj) == cl)
    sub_dapc <- dapc_obj$ind.coord[inds, , drop = FALSE]
    var_total <- sum(apply(sub_dapc, 2, var))
    variance_list[[cl]] <- var_total
  }
  
  return(variance_list)
}

variance_by_cluster <- cluster_variance(dapc, gl)


# Run the function
variance_by_cluster <- cluster_variance(dapc_result, gl)

# Print results
print(variance_by_cluster)


#DAPC with outliers excluded
#Excluding outlier clusters
# Identify the clusters to exclude

strata_data_filtered <- strata_data %>%
  filter(!cluster %in% c(3, 10, 4, 3, 8, 7, 9))
# Step 4: Subset the genlight object
gl_filtered <- gl[indNames(gl) %in% rownames(strata_data_filtered), ]
strata_data_filtered <- strata_data_filtered[match(indNames(gl_filtered), rownames(strata_data_filtered)), ]
strata(gl_filtered) <- strata_data_filtered
#Make sure that the outlier clusters have been excluded
setPop(gl_filtered) <-  ~cluster
popNames(gl_filtered)
gl_filtered
#Check that the groups have more than one member
#Exclude the group that does not
pop_sizes <- table(pop(gl_filtered))
valid_pops <- names(pop_sizes[pop_sizes > 1])
# Subset the genlight object
gl_filtered <- gl_filtered[pop(gl_filtered) %in% valid_pops, ]

#Repeat DAPC
setPop(gl_filtered) <- ~cluster
popNames(gl_filtered)
cols3 <- c("darkblue","purple","darkgreen","orange","red" ,"magenta" ,"gold","green", "lightblue", "pink", "darkred")

dapc_cluster <- dapc(gl_filtered, var.contrib = TRUE, scale = FALSE, n.pca = 11, n.da = nPop(gl) - 1)

scatter(dapc_cluster, col = cols3, clabel=FALSE, xax = 1, yax = 2,inset.solid = 1, pch=20, solid = 0.6, cstar=0, cellipse = 2.5, legend = TRUE, cleg = 0.75, posi.leg = "bottomright", cex = 1.5, scree.da = TRUE, scree.pca = TRUE, posi.pca = "topleft", ratio.da = 0.15, ratio.pca = 0.15, inset.da = 0.03, inset.pca = 0.03, posi.da = "bottomleft")
dapc_cluster$var
xlim <- par("usr")[1:2]
mtext(paste0("Variance explained: ", round(dapc_cluster$var, 4)),
      side = 3, line = -2, at = xlim[2] - 0.25 * diff(xlim), cex = 0.9)
#Identify the controls in the mess
inds <- colnames(filtered_vcf@gt)[2:length(colnames(filtered_vcf@gt))]

# Create the isolate_table with the correct control identification
library(ggrepel)
isolate_table <- tibble(id = inds, 
                        isolate = sub('^2Ar', 'AR', inds), 
                        control = grepl("^AR", sub('^2Ar', 'AR', inds)), 
                        replicate = grepl("2Ar", inds))

isolate_table <- isolate_table %>%
  mutate(comment = case_when(
    grepl("Ar0", id) ~ "Control",
    grepl("AR", id) ~ "Control",
    TRUE ~ NA_character_
  ))

controls <- isolate_table %>%
  filter(comment == "Control") %>%
  pull(id)

plot_data <- dapc_cluster$ind.coord %>%
  as.data.frame() %>%
  rownames_to_column("isolate") %>%
  mutate(is_control = isolate %in% controls)


# ggplot(dapc_filt_cluster$ind.coord %>% as.data.frame() %>% rownames_to_column("isolate"), 
#        aes(x = LD1, y = LD2, label = isolate, colour = isolate %in% controls)) +
#   geom_point() + 
#   geom_label_repel(max.overlaps = Inf, box.padding = 0.2, point.padding = 0.2, 
#                    segment.size = 0.3) + 
#   scale_color_manual(name = "Key", 
#                      values = c("FALSE" = "blue", "TRUE" = "red"), 
#                      labels = c("FALSE" = "Progeny", "TRUE" = "Control")) 

ggplot(plot_data, aes(x = LD1, y = LD2, colour = is_control)) +
  geom_point() +
  geom_label_repel(
    data = subset(plot_data, is_control),
    aes(label = isolate),
    max.overlaps = Inf,
    box.padding = 0.2,
    point.padding = 0.2,
    segment.size = 0.3
  ) +
  scale_color_manual(
    name = "Key",
    values = c("FALSE" = "blue", "TRUE" = "red"),
    labels = c("FALSE" = "Disseminents", "TRUE" = "Control")
  )

#Pairwise of DAPC clusters based off Fst/compared with SNP contribution to clustering####
#must be done using a genind object. Re-import the data and line it up to 'updated_metadata"
library(hierfstat)
vcf_filtered <- read.vcfR("output/2024_filtered_vcf_inform_5%.vcf.gz")

gl_dist <- vcfR2genind(vcf_filtered)
metadata <- read_excel("output/dapc11_updated_metadata.xlsx")
colnames(metadata)
filtered_metadata <- metadata[metadata$sample_ID %in% indNames(gl_dist), ]
gl_df <- data.frame(genotype = indNames(gl_dist))
joined_data <- gl_df %>%
  left_join(filtered_metadata, by = c("genotype" = "sample_ID"))

strata_data <- joined_data %>%
  dplyr::select(genotype, everything()) %>%
  column_to_rownames(var = "genotype")

strata(gl_dist) <- strata_data

popNames(gl_dist)
setPop(gl_dist) <- ~cluster
gl_dist
ploidy(gl_dist) <- 1
#Get the eigenvalues or make clusters into factors
#eig.val <- dapc$eig
clusters <- dapc$assign
gl@pop <- as.factor(clusters)
hf_data <- genind2hierfstat(gl_dist)
#Pairwise test
dist_matrix <- hierfstat::genet.dist(hf_data, diploid = FALSE, method = "Dch")
#Plot as a heatmap
#Plot as a heatmap
dist_df <- as.data.frame(as.matrix(dist_matrix))
dist_df$Var1 <- rownames(dist_df)

library(reshape2)
library(ggplot2)

# Melt the distance matrix
dist_melt <- melt(dist_df, id.vars = "Var1")

pop_map <- c("1" = "MLG.8", "2" = "MLG.3", "3" = "MLG.11", "4" = "MLG.7", "5" = "MLG.2", "6" = "MLG.9", "7" = "MLG.1", "8" = "MLG.10", 
             "9" = "MLG.5", "10" = "MLG.6", "11" = "MLG.4")

# Apply mapping to Var1 and variable columns
dist_melt$Var1_mapped <- pop_map[as.character(dist_melt$Var1)]
dist_melt$variable_mapped <- pop_map[as.character(dist_melt$variable)]

# Define the desired MLG order
mlg_order <- paste0("MLG.", 1:11)

# Convert mapped columns to factors with specified levels
dist_melt$Var1_mapped <- factor(dist_melt$Var1_mapped, levels = mlg_order)
dist_melt$variable_mapped <- factor(dist_melt$variable_mapped)

# Plot the heatmap
ggplot(dist_melt, aes(x = Var1_mapped, y = variable_mapped, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(value, 3)), size = 3) +
  scale_fill_viridis_c(name = "Genetic Distance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Genetic Distance Heatmap", x = "Population", y = "Population")


#Loading plot comparison

loading_scores <- dapc$var.contr
all_load <- loadingplot(dapc$var.contr, axis = 1, thres = 0.005)
#Select DAPC clusters with little to no differentiation 
selected_clusters <- c("11", "8", "7", "5")

inds <- which(pop(gl) %in% selected_clusters)
genind_sub <- gl[inds, ]
pop(genind_sub) <- factor(pop(genind_sub)) 

dapc_sub <- dapc(genind_sub)
#3
#2
loadingplot(dapc_sub$var.contr, axis = 1, thres = 0.005)

#Select DAPC clusters with little to no differentiation 
selected_clusters <- c("11", "5")
inds <- which(pop(gl) %in% selected_clusters)
genind_sub <- gl[inds, ]

pop(genind_sub) <- factor(pop(genind_sub)) 
#check no missing data was introduced
dapc_sub <- dapc(genind_sub)# Only one discriminant axis for two groups
#20
#1
head(rownames(dapc_sub$var.contr))
loci_names <- locNames(gl)
rownames(dapc_sub$var.contr) <- loci_names
loadingplot(dapc_sub$var.contr, axis = 1, thres = 0.005)

threshold <- 0.005
loadings <- dapc_sub$var.contr[, 1]
loci_above_threshold <- names(loadings[abs(loadings) > threshold])
#[1] "ArME14_ctg_13_220743" "ArME14_ctg_18_639668" "ArME14_ctg_21_482082"

#Select DAPC clusters with little to no differentiation 
selected_clusters <- c("7", "8")
inds <- which(pop(gl) %in% selected_clusters)
genind_sub <- gl[inds, ]

pop(genind_sub) <- factor(pop(genind_sub)) 
#check no missing data was introduced
dapc_sub <- dapc(genind_sub)# Only one discriminant axis for two groups
#30
#1
head(rownames(dapc_sub$var.contr))
loci_names <- locNames(gl)
rownames(dapc_sub$var.contr) <- loci_names
loadingplot(dapc_sub$var.contr, axis = 1, thres = 0.006)

threshold <- 0.006
loadings <- dapc_sub$var.contr[, 1]
loci_above_threshold <- names(loadings[abs(loadings) > threshold])
#[1] "ArME14_ctg_05_135536" "ArME14_ctg_06_289485" "ArME14_ctg_08_590181" "ArME14_ctg_09_21335"  "ArME14_ctg_13_966205" "ArME14_ctg_13_966228"

#Select DAPC clusters with little to no differentiation 
selected_clusters <- c("7", "8", "11")
inds <- which(pop(gl) %in% selected_clusters)
genind_sub <- gl[inds, ]

pop(genind_sub) <- factor(pop(genind_sub)) 
dapc_sub <- dapc(genind_sub)# Only one discriminant axis for two groups
#20
#2
head(rownames(dapc_sub$var.contr))
loci_names <- locNames(gl)
rownames(dapc_sub$var.contr) <- loci_names
loadingplot(dapc_sub$var.contr, axis = 1, thres = 0.006)

threshold <- 0.0062
loadings <- dapc_sub$var.contr[, 1]
loci_above_threshold <- names(loadings[abs(loadings) > threshold])
#[1] "ArME14_ctg_11_314007" "ArME14_ctg_18_639668"

#Select DAPC clusters with little to no differentiation 
selected_clusters <- c("9", "3")
inds <- which(pop(gl) %in% selected_clusters)
genind_sub <- gl[inds, ]

pop(genind_sub) <- factor(pop(genind_sub)) 
#check no missing data was introduced
dapc_sub <- dapc(genind_sub)# Only one discriminant axis for two groups
#12
#1
head(rownames(dapc_sub$var.contr))
loci_names <- locNames(gl)
rownames(dapc_sub$var.contr) <- loci_names
loadingplot(dapc_sub$var.contr, axis = 1, thres = 0.05)

threshold <- 0.05
loadings <- dapc_sub$var.contr[, 1]
loci_above_threshold <- names(loadings[abs(loadings) > threshold])
#[1] "ArME14_ctg_03_397821"  "ArME14_ctg_04_78320"   "ArME14_ctg_05_191687"  "ArME14_ctg_07_163121"  "ArME14_ctg_07_1410954" "ArME14_ctg_10_1865349"
#[7] "ArME14_ctg_15_1103277" "ArME14_ctg_17_1286713" "ArME14_ctg_20_961346"

#Loading plot with only my control genotypes
selected_clusters <- c("3", "10", "4", "8", "7", "9")
inds <- which(pop(gl) %in% selected_clusters)
genind_sub <- gl[inds, ]

dapc_sub <- dapc(genind_sub)
#5
#4
loadingplot(dapc_sub$var.contr, axis = 1, thres = 0.0021)
#NULL

#DAPC with populations pre-defined based on host
setPop(gl) <- ~Host
popNames(gl)
cols <- c("darkorange1", "darkgreen", "red", "blue")
# need to define the no. of PCs retained for DAPC, which can have a substantial impact 
# on the results. Cross-validation (xvalDapc) is an objective way to identify no. of PCs to retain.
set.seed(999)
xvalgl <- xvalDapc(tab(gl, NA.method = "mean"), pop(gl), training.set = 0.9)
set.seed(999)
xvalgl <- xvalDapc(tab(gl, NA.method = "mean"), pop(gl), parallel = "multicore", ncpus = 4L, training.set = 0.9, n.rep = 1000, n.pca = 18:40)
xvalgl$`Number of PCs Achieving Highest Mean Success`
#20
popNames(gl)
dapc_host <- dapc(gl, var.contrib = TRUE, scale = FALSE, n.pca = 20, n.da = nPop(gl) - 1)
#final:
scatter(dapc_host, col = cols, clabel=FALSE, xax = 1, yax = 2,inset.solid = 1, pch=20, solid = 0.6, cstar=0, cellipse = 2.5, legend = FALSE, cleg = 0.75, posi.leg = "bottomright", cex = 1.5, scree.da = FALSE, scree.pca = FALSE, posi.pca = "topleft", ratio.da = 0.15, ratio.pca = 0.15, inset.da = 0.03, inset.pca = 0.03, posi.da = "bottomleft")
# The ellipses represent the maximum area spanned by 95% of the data in a population by year of collection
scatter(dapc_host, col = cols, clabel=FALSE, xax = 1, yax = 2,inset.solid = 1, pch=20, solid = 0.6, cstar=0, cellipse = 2.5, legend = TRUE, cleg = 0.75, posi.leg = "bottomright", cex = 1.5, scree.da = TRUE, scree.pca = TRUE, posi.pca = "topleft", ratio.da = 0.15, ratio.pca = 0.15, inset.da = 0.03, inset.pca = 0.03, posi.da = "bottomleft")
dapc_host$var
mtext(paste0("Variance explained: ", round(dapc_host$var, 4)),
      side = 3, line = -2, adj = 1, cex = 0.9)

contrib_Host_axis1 <- loadingplot(dapc_host$var.contr, axis = 1, thres = 0.05)
contrib_Host_axis1
contrib_Host_axis2 <- loadingplot(dapc_host$var.contr, axis = 2, thres = 0.05)
contrib_Host_axis2

#By priory DAPC clustering
#DAPC with populations pre-defined based on host
setPop(gl) <- ~cluster
popNames(gl)
cols2 <- c("darkblue","purple","darkgreen","orange","red","brown" ,"magenta" ,"gold","green", "lightblue", "pink", "grey", "darkred", "black" )
# need to define the no. of PCs retained for DAPC, which can have a substantial impact 
# on the results. Cross-validation (xvalDapc) is an objective way to identify no. of PCs to retain.
set.seed(999)
xvalgl <- xvalDapc(tab(gl, NA.method = "mean"), pop(gl), training.set = 0.9)
set.seed(999)
#xvalgl <- xvalDapc(tab(gl, NA.method = "mean"), pop(gl), parallel = "multicore", ncpus = 4L, training.set = 0.9, n.rep = 1000, n.pca = 18:70)
xvalgl$`Number of PCs Achieving Highest Mean Success`
#20
popNames(gl)
dapc_cluster <- dapc(gl, var.contrib = TRUE, scale = FALSE, n.pca = 11, n.da = nPop(gl) - 1)
#final:
scatter(dapc_cluster, col = cols2, clabel=FALSE, xax = 1, yax = 2,inset.solid = 1, pch=20, solid = 0.6, cstar=0, cellipse = 2.5, legend = FALSE, cleg = 0.75, posi.leg = "bottomright", cex = 1.5, scree.da = FALSE, scree.pca = FALSE, posi.pca = "topleft", ratio.da = 0.15, ratio.pca = 0.15, inset.da = 0.03, inset.pca = 0.03, posi.da = "bottomleft")
# The ellipses represent the maximum area spanned by 95% of the data in a population by year of collection
scatter(dapc_cluster, col = cols2, clabel=FALSE, xax = 1, yax = 2,inset.solid = 1, pch=20, solid = 0.6, cstar=0, cellipse = 2.5, legend = TRUE, cleg = 0.75, posi.leg = "bottomright", cex = 1.5, scree.da = TRUE, scree.pca = TRUE, posi.pca = "topleft", ratio.da = 0.15, ratio.pca = 0.15, inset.da = 0.03, inset.pca = 0.03, posi.da = "bottomleft")
dapc_cluster$var
xlim <- par("usr")[1:2]
mtext(paste0("Variance explained: ", round(dapc_cluster$var, 4)),
      side = 3, line = -2,at = xlim[2] - 0.3 * diff(xlim), cex = 0.9)

contrib_cluster_axis1 <- loadingplot(dapc_cluster$var.contr, axis = 1, thres = 0.05)
contrib_cluster_axis1
contrib_cluster_axis2 <- loadingplot(dapc_cluster$var.contr, axis = 2, thres = 0.05)
contrib_cluster_axis2

#DAPC by priory clustering 
#Excluding outlier clusters
# Identify the clusters to exclude

strata_data_filtered <- strata_data %>%
  filter(!cluster %in% c(1, 10, 12, 2, 4, 7))
# Step 4: Subset the genlight object
gl_filtered <- gl[indNames(gl) %in% rownames(strata_data_filtered), ]
strata_data_filtered <- strata_data_filtered[match(indNames(gl_filtered), rownames(strata_data_filtered)), ]
strata(gl_filtered) <- strata_data_filtered
#Make sure that the outlier clusters have been excluded
popNames(gl_filtered)
setPop(gl_filtered) <-  ~cluster
gl_filtered
#Check that the groups have more than one member
#Exclude the group that does not
pop_sizes <- table(pop(gl_filtered))
valid_pops <- names(pop_sizes[pop_sizes > 1])
# Subset the genlight object
gl_filtered <- gl_filtered[pop(gl_filtered) %in% valid_pops, ]

#Repeat DAPC
setPop(gl_filtered) <- ~cluster
popNames(gl_filtered)
cols3 <- c("darkblue","purple","darkgreen","orange","red" ,"magenta" ,"gold","green", "lightblue", "pink", "darkred", "brown")
set.seed(999)
xvalgl_filt <- xvalDapc(tab(gl_filtered, NA.method = "mean"), pop(gl_filtered), training.set = 0.9)
#xvalgl_filt <- xvalDapc(tab(gl_filtered, NA.method = "mean"), pop(gl_filtered), parallel = "multicore", ncpus = 4L, training.set = 0.9, n.rep = 100, n.pca = 40:60)
xvalgl_filt$`Number of PCs Achieving Highest Mean Success`
#20
popNames(gl_filtered)
dapc_filt_cluster <- dapc(gl_filtered, var.contrib = TRUE, scale = FALSE, n.pca = 11, n.da = nPop(gl_filtered) - 1)
#final:
scatter(dapc_filt_cluster, col = cols3, clabel=FALSE, xax = 1, yax = 2,inset.solid = 1, pch=20, solid = 0.6, cstar=0, cellipse = 2.5, legend = FALSE, cleg = 0.75, posi.leg = "bottomright", cex = 1.5, scree.da = FALSE, scree.pca = FALSE, posi.pca = "topleft", ratio.da = 0.15, ratio.pca = 0.15, inset.da = 0.03, inset.pca = 0.03, posi.da = "bottomleft")
# The ellipses represent the maximum area spanned by 95% of the data in a population by year of collection
scatter(dapc_filt_cluster, col = cols3, clabel=FALSE, xax = 1, yax = 2,inset.solid = 1, pch=20, solid = 0.6, cstar=0, cellipse = 2.5, legend = TRUE, cleg = 0.75, posi.leg = "bottomright", cex = 1.5, scree.da = TRUE, scree.pca = TRUE, posi.pca = "topleft", ratio.da = 0.15, ratio.pca = 0.15, inset.da = 0.03, inset.pca = 0.03, posi.da = "bottomleft")
dapc_filt_cluster$var
mtext(paste0("Variance explained: ", round(dapc_filt_cluster$var, 4)),
      side = 3, line = -2, adj = 1, cex = 0.9)

contrib_cluster_axis1 <- loadingplot(dapc_filt_cluster$var.contr, axis = 1, thres = 0.05)
contrib_cluster_axis1
contrib_cluster_axis2 <- loadingplot(dapc_filt_cluster$var.contr, axis = 2, thres = 0.05)
contrib_cluster_axis2

dapc_filt_cluster$loadings


#Not applicablez, target training set too small                    
# #DAPC using control set as training data
# #Use as a genlight object as recommended for genome wide snp data
# gl <- vcfR2genlight(filtered_vcf)
# 
# selected_inds <- c("AR0022", "AR0023", "AR0212", "AR0020", "AR0210", "AR0128")
# gl.train <- gl[which(indNames(gl) %in% selected_inds), ]
# gl.test <- gl[which(!indNames(gl) %in% selected_inds), ]
# indNames(gl.train)
# nInd(gl.test)
# #Assign the 6 control isolates as it's own single population otherwise there is no within group population
# #n.da must be the no. of groups-1 
# 
# 
# train_pops <- factor(rep("TrainGroup", nInd(gl.train)))
# 
# 
# dapc.train <- dapc(gl.train, pop = train_pops, n.pca = 5, n.da = 1)

# Distance tree ####
#Neighbor joining
#Check that it is a genlight object. Use gl first for all 355 samples. Then use gl_filtered.
#Make sure the population is assigned
#NJ tree
#gl
ind_dist <- bitwise.dist(gl)
indNames(gl)
# When differences_only = TRUE, the output will reflect the number of different loci. 
# The default setting, differences_only = FALSE, reflects the number of different alleles
# Note: this has no effect on haploid organisms since 1 locus = 1 allele
# Missing match determines whether two samples differing by missing data in a location should be counted as matching at that location. 
# Default set to TRUE, which forces missing data to match with anything. 
# FALSE forces missing data to not match with any other information, including other missing data
# percent = TRUE is default, that is, distance represented from 0 to 1
# FALSE will return the distance represented as integers from 1 to n where n is the number of loci.
# This option has no effect if euclidean = TRUE
ind_dist_no_of_loci <- bitwise.dist(gl, percent = FALSE)
library(MASS)
write.matrix (ind_dist_no_of_loci, file = "output/ind_dist_no_of_loci.gl")
cols3 <- c("darkblue","purple","darkgreen","orange","red" ,"magenta" ,"gold","green", "lightblue", "pink", "darkred", "brown")
#make tree
gl_genlight_tree <- aboot(gl, tree = "nj", distance = bitwise.dist, sample = 1000, showtree = T, cutoff = 70, quiet = T, missing = ignore)
tip_colors <- cols3[as.numeric(gl.clusters$grp)]
plot(gl_genlight_tree, type = "phylogram" , edge.w = 2, cex = 0.5, tip.color = tip_colors)
plot(gl_genlight_tree, type= "radial", edge.w=2, cex=0.5)

library(ape)
ape::write.tree(gl_genlight_tree, file='output/chp_4/individual_NJ_tree_gl.txt')

#gl_filtered
#Includes only clusters corresponding with controls
ind_dist_no_of_loci <- bitwise.dist(gl_filtered, percent = FALSE)
library(MASS)
write.matrix (ind_dist_no_of_loci, file = "output/chp_4/ind_dist_no_of_loci.gl_filtered")

gl_filtered_genlight_tree <- aboot(gl_filtered, tree = "nj", distance = bitwise.dist, sample = 1000, showtree = T, cutoff = 70, quiet = T, missing = ignore)
tip_colors <- cols3[as.numeric(gl.clusters$grp)]
plot(gl__filtered_genlight_tree, type = "phylogram" , edge.w = 2, cex = 0.5, tip.color = tip_colors)


library(ape)
ape::write.tree(gl__filtered_genlight_tree, file='output/chp_4/individual_NJ_tree_gl_filtered.txt')

#Write the vcf out as a fasta string 
library(dplyr)
library(tidyr)
library(seqinr)

# Extract the first allele from each genotype
genotype_table <- extract.gt(filtered_vcf, return.alleles = TRUE) %>% substr(., 1, 1) %>% as.data.frame() %>% t() %>% as.data.frame()

genotype_table[is.na(genotype_table)] <- "N"
# Create a FASTA-formatted string
genotype_table_string <- genotype_table %>% unite("fasta", everything(), sep = "") %>%
  rownames_to_column("isolate") %>% mutate(fasta = paste0(">", isolate, "\n", fasta))

# Write to FASTA file
writeLines(genotype_table_string$fasta, "output/filtered_vcf_2024.fasta")

# Determining significant groups of loci ####

# #Minimum spanning network####
# 
vcf <- read.vcfR("output/2024_filtered_vcf_inform_5%.vcf.gz")

gen_obj <- vcfR2genind(vcf)

ploidy(gen_obj) <- 1
tab(gen_obj)[1:5,1:5]
X <- genind2df(gen_obj)
head(X)
gen_obj <- df2genind(X, ploidy = 1, ncode=1)
#Need a genclone object for this
gen <- poppr::as.genclone(gen_obj)
#Set metadata
 metadata <- read_excel("output/chp_4/dapc11_updated_metadata.xlsx")
 colnames(metadata)
 filtered_metadata <- metadata[metadata$sample_ID %in% indNames(gen), ]
 genind_df <- data.frame(genotype = indNames(gen))
 
 joined_data <- genind_df %>%  left_join(filtered_metadata, by = c("genotype" = "sample_ID"))

strata(gen) <- joined_data
popNames(gen)
setPop(gen) <- ~ cluster
# #Check the genclone info
# # gen
# # This is a genclone object
# # -------------------------
# #   Genotype information:
# #   
# #   248 original multilocus genotypes 
# # 336 haploid individuals
# # 1290 codominant loci
# # 
# # Population information:
# #   
# #   15 strata - 
# #   genotype, No, species, ..., Disease_score, Sampling_no., cluster
# # 14 populations defined - 3, 2, 14, ..., 9, 6, 10
# 
# #interactive shiny msn
imsn()

gen_sub <- popsub(gen, exclude = character(0))
gen_nomiss <- missingno(gen, type = 'mean')
gen_dist <- provesti.dist(gen_nomiss)
min_span_net <- poppr.msn(gen_sub, gen_dist, showplot = FALSE, include.ties = TRUE)
#As a tree
set.seed(69)
plot_poppr_msn(gen,
               min_span_net,
               inds = "0",
               mlg = FALSE,
               gadj = 14,
               nodescale = 30,
               palette = funky,
               cutoff = NULL,
               quantiles = FALSE,
               beforecut = TRUE,
               pop.leg = FALSE,  # suppress population legend
               size.leg = TRUE, # suppress node size legend
               scale.leg = TRUE,# suppress scale bar
               layfun = igraph::layout_as_tree)

# Make it interactive?
g <- min_span_net$graph

pop_assignments <- pop(gen_sub)
ind_names <- indNames(gen_sub)
names(pop_assignments) <- indNames(gen_sub)

node_pops <- pop_assignments[V(g)$name]

pop_colors <- min_span_net$colors

node_colors <- pop_colors[as.character(node_pops)]
names(pop_colors)

nodes <- data.frame(
  id = V(g)$name,
  label = V(g)$name,
  title = paste("Sample:", V(g)$name, "<br>Population:", node_pops),
  color = node_colors
)


edges <- as_data_frame(g, what = "edges")

visNetwork(nodes, edges) %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)

# Generate the color palette
palette_colors <- funky(length(unique(pop(gen))))

# Then use it in the legend
legend("bottomleft",
       legend = sort(unique(pop(gen))),
       fill = palette_colors,
       title = "Populations",
       cex = 0.8)
#Layout nicely
set.seed(69)
plot_poppr_msn(gen,
               min_span_net,
               inds = "0",
               mlg = FALSE,
               gadj = 14,
               nodescale = 30,
               palette = funky,
               cutoff = NULL,
               quantiles = FALSE,
               beforecut = TRUE,
               pop.leg = FALSE,  # suppress population legend
               size.leg = TRUE, # suppress node size legend
               scale.leg = TRUE,# suppress scale bar
               layfun = igraph::layout_nicely)

# Generate the color palette
palette_colors <- funky(length(unique(pop(gen))))

# Then use it in the legend
legend("bottomleft",
       legend = sort(unique(pop(gen))),
       fill = palette_colors,
       title = "Populations",
       cex = 0.8)

#Loading plots per populations
#Make a function to extract the right chromosome
locNames(gl)

extract_chromosome <- function(snp_names) {
  sapply(snp_names, function(name) {
    match <- regmatches(name, regexpr("ctg_\\d+", name))
    return(match)
  })
}

#Extracti loading values
loci_names <- locNames(gl)

loadings <- dapc$var.contr
rownames(loadings) <- loci_names

loading_df <- as.data.frame(loadings)
loading_df$SNP <- rownames(loadings)

loading_df$Chromosome <- extract_chromosome(loading_df$SNP)

setPop(gl) <- ~host
popNames(gl)

ggplot(loading_df, aes(x=Chromosome, y=Axis1)) +
  geom_point(aes(color=Cluster, shape=Cultivar)) +
  geom_hline(yintercept=0.01, color="red", linetype="dashed") +
  labs(title="DAPC Loadings for SNPs", x="Chromosome", y="Loading Value") +
  theme_minimal()
# https://adegenet.r-forge.r-project.org/files/Glasgow2015/practical-GWAS_day4.pdf




