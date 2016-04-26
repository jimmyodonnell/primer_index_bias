# Code for analyses in O’Donnell et al. 2015. Indexed primers induce template-specific bias in large-scale DNA sequencing studies. PLoS One

# Compare the pairwise similarity of sequence data from samples amplified using two protocols

# See also function vegan::betadisper() and vegan::mrpp()
# though this assumes the same variables (species/OTUs) have been measured across samples (i.e. OTU1 is the same sequence in both sets)

# If RUN1 has 3 libraries and 2 tags per sample, and RUN2 has 2 libraries and 2 tags per sample, this is how many comparisons to expect:
# RUN, WITHIN_TAG, BETWEEN TAG
# RUN1, 6, 9
# RUN2, 1, 4

library(vegan)

setwd("~/Documents/GoogleDrive/Kelly_Lab/Projects/Lemonade/Analysis")

data_dir <- file.path("..", "Data")
fig_dir <- file.path("..", "Figures")


# Read in OTU table for each experiment
DAT_OTU <- list()
# DAT_OTU[[1]] <- t(read.csv("/Users/threeprime/temp_big_upload/16S/Analysis_20150614_1927/all_lib/OTU_table.csv", row.names = 1))
DAT_OTU[[1]] <- t(read.csv(file.path(data_dir, "OTU_table_1PCR.csv"), row.names = 1))
# DAT_OTU[[2]] <- t(read.csv("/Users/threeprime/temp_big_upload/16S/Analysis_20150522_0432/all_lib/OTU_table.csv", row.names = 1))
DAT_OTU[[2]] <- t(read.csv(file.path(data_dir, "OTU_table_2PCR.csv"), row.names = 1))

# OTU rownames should be unique sequencing sample identifiers.
# rownames(DAT_OTU[[1]])
# rownames(DAT_OTU[[2]])

# Set the first and last character that identifies the library and tag:
lib_start <- 1
lib_stop  <- 4
tag_start <- 10
tag_stop <- 15


# compute reads per sample before altering raw data
reads_per_sample <- lapply(DAT_OTU, rowSums)
# save(reads_per_sample, file = "reads_per_sample.rdata")

# remove chimeras
for(i in 1:length(DAT_OTU)){
  DAT_OTU[[i]] <- DAT_OTU[[i]][,-which(colnames(DAT_OTU[[i]]) == "chimera")]
}

# convert to proportions, and replace NaN with 0
for(i in 1:length(DAT_OTU)){
  DAT_OTU[[i]] <- DAT_OTU[[i]]/rowSums(DAT_OTU[[i]])
  DAT_OTU[[i]] <- replace(DAT_OTU[[i]], is.nan(DAT_OTU[[i]]), 0)
}

# calculate total number of OTUs per each dataset
N_otus_per_dataset <- lapply(DAT_OTU, ncol)



# Read in metadata for each sequencing run
DAT_pool <- list()
DAT_pool[[1]] <- read.csv("/Users/threeprime/temp_big/run_20141113_time_series/sample_data_update.csv")
DAT_pool[[2]] <- read.csv("/Users/threeprime/temp_big/run_20150401/20150317_sequencing_pool.csv")

# The sample pool data must include a column named "sample_type", and you must specify which type of sample to include
sample_type_include <- "environmental"

# which columns are in both of the sequencing metadata files?
# Reduce(intersect, lapply(DAT_pool, colnames))

# split(DAT_pool[[2]]$tag_sequence, DAT_pool[[2]]$sample_name)


# In order to make it easier to link metadata and sequencing data
for(i in 1:length(DAT_pool)){
  # make a vector that will link the OTU file and the sequencing pool
  lib_tag <- paste(DAT_pool[[i]]$library, "tag", DAT_pool[[i]]$tag_sequence, sep = "_")
  
  # ... and add it to the sequencing pool metadata
  DAT_pool[[i]] <- data.frame(DAT_pool[[i]], lib_tag)
}
  


################################################################################
# Define a function to reformat distance matrices
reformat_dists <- function(x){
  
  TEMP_labels <- t(combn(labels(x), 2))
  
  TEMP_mat <- matrix(
    data = cbind(
      substr(x = as.character(TEMP_labels), start = lib_start, stop = lib_stop),
      substr(x = as.character(TEMP_labels), start = tag_start, stop = tag_stop)
    ), 
    ncol = 4 , 
    dimnames = list(NULL, c("lib1", "lib2", "tag1", "tag2"))
  )
  
  TEMP <- data.frame(TEMP_mat, 
                     dist=as.numeric(x)
  )
  
  return(TEMP)
}

# define a function to get tag levels
# tag_levels <- function(X){
#   tag_df <- do.call(
#     rbind,
#       lapply(
#         lapply(
#           split(
#             X[,c("tag_sequence", "lib_tag")],
#             X$sample_name
#           ), 
#           droplevels
#         ), 
#         function(x) {x$tag_sequence <- as.numeric(x$tag_sequence); return(x) }
#       )
#     )
#   tag_lev <- tag_df$tag_sequence[match(tag_df$lib_tag, X$lib_tag)]
#   return(tag_lev)
# }
# 
# # bind tag levels to metadata
# for(i in 1:length(DAT_pool)){
#   taglev <- tag_levels(DAT_pool[[i]])
#   DAT_pool[[i]] <- 
#     data.frame(DAT_pool[[i]], tag_level = taglev)
# }


# Initiate a list to store distance matrices...
dist_samples <- list()
dist_reform <- list()
within_index <- list()
among_index <- list()

# ... and for each sequencing run...
for(i in 1:length(DAT_OTU)){
# for(i in 1){
  # link the pool file to the OTU file
  tag_to_sequencing_data <- match(rownames(DAT_OTU[[i]]), DAT_pool[[i]]$lib_tag)
  tag_to_samplename <- DAT_pool[[i]][tag_to_sequencing_data, "sample_name"]
  
  # DAT_pool[[i]][DAT_pool[[i]]$sample_type == sample_type_include,][tag_to_sequencing_data, "sample_name"]
  
  # boxplot(DAT_OTU[[i]][,5] ~ tag_to_samplename)
  
  # Compute pairwise dissimilarity for replicates of each sample
  dist_samples[[i]] <- lapply(split(as.data.frame(DAT_OTU[[i]]), tag_to_samplename), vegdist, method = "bray", upper = TRUE, diag = TRUE) # Add [,1:NUMBER] to include only first NUMBER OTUs
  
  # Exclude samples that don't match the intended sample type
  names_include <- unique(DAT_pool[[i]][DAT_pool[[i]]$sample_type == sample_type_include,]$sample_name)
  dist_samples[[i]] <- dist_samples[[i]][names_include]
  
  # reformat the distances
  dist_reform[[i]] <- lapply(dist_samples[[i]], reformat_dists)
  
  # subset relevant values
  # WITHIN index pairwise dissimilarity:
  within_index[[i]] <- lapply(dist_reform[[i]], function(x) subset(x$dist, x[,"tag1"] == x[,"tag2"]))
  
  # AMONG index pairwise dissimilarity:
  among_index[[i]] <- lapply(dist_reform[[i]], function(x) subset(x$dist, x[,"tag1"] != x[,"tag2"]))
}

# concatenate vectors of among_ and within_ index at first level
among <- lapply(among_index, function(x) do.call(what = c, args = x))
within <- lapply(within_index, function(x) do.call(what = c, args = x))

# put lists together by treatment
PCR_single <- list(within[[1]], among[[1]])
PCR_double <- list(within[[2]], among[[2]])
names(PCR_single) <- c("within", "among")
names(PCR_double) <- c("within", "among")

# conduct t-tests
t.test(x = PCR_single[[1]], y = PCR_single[[2]])
t.test(x = PCR_double[[1]], y = PCR_double[[2]])

# plot results
# boxplot(c(PCR_single, PCR_double), ylim = c(0,1))

# get means for each comparison (within, among index) for each treatment (single or double PCR)
mean_dist_within <- lapply(within_index, function(x) sapply(x, mean))
mean_dist_among <- lapply(among_index, function(x) sapply(x, mean))

# put lists together by treatment
PCR_single_mean <- list(mean_dist_within[[1]], mean_dist_among[[1]])
PCR_double_mean <- list(mean_dist_within[[2]], mean_dist_among[[2]])
names(PCR_single_mean) <- c("within", "among")
names(PCR_double_mean) <- c("within", "among")

# perform logit transformation:
logit_transform <- function(x) log(x/(1-x))
PCR_single_mean_logit <- lapply(PCR_single_mean, logit_transform)
PCR_double_mean_logit <- lapply(PCR_double_mean, logit_transform)

# For the single PCR treatment, there was little variation among OTU sequence counts generated with the same indexed primers (Bray-Curtis dissimilarity index:
mean(PCR_single_mean$within);
sd(PCR_single_mean$within)
# , FIGURE 1).
# By contrast—-and strikingly contrary to the assumption that multiplexing tags do not influence analytical outcomes—-there were large differences between OTU counts generated using different tagged primer sets on the same environmental sample (Bray-Curtis dissimilarity index; 
mean(PCR_single_mean$among);
sd(PCR_single_mean$among)
# , and these were significantly greater than comparisons within the same primer index (Welch Two Sample t(
t_means_single <- t.test(x = PCR_single_mean$among, y = PCR_single_mean$within)
t_means_single$statistic;
# p = 
t_means_single$p.value

# However, using a double-PCR approach kept dissimilarity low
# both within
mean(PCR_double_mean$within)
sd(PCR_double_mean$within)
# and among 
mean(PCR_double_mean$among)
sd(PCR_double_mean$among)
# primer index replicates. There was no difference in the mean dissimilarity between the within and among primer index comparisons using double PCR ((Welch Two Sample t(.
t_means_double <- t.test(x = PCR_double_mean$among, y = PCR_double_mean$within)
t_means_double$statistic
# p = 
t_means_double$p.value


# INCLUDE OR NOT?
# Finally, there was no difference in mean dissimilarity within tags for the single and double PCR approach.
t.test(x = PCR_double_mean$within, y = PCR_single_mean$within) 



########################################################################################################
# PLOTTING RESULTS
########################################################################################################
plot_input <- c(PCR_single_mean, PCR_double_mean)
# plot_input <- c(PCR_single_mean_logit, PCR_double_mean_logit)
# , FIGURE 1)
# plot results
### FOR PDF
# pdf("BC_dis_boxplot.pdf")

### FOR EPS
# setEPS()
# postscript("BC_dis_boxplot.eps")
par(mar = c(5,4,1,1))
boxplot(
  plot_input, 
  ylab = "Mean Bray-Curtis Dissimilarity", 
  ylim = c(0,1), 
  xlab = "Primer Index Comparison"
)
abline(v = 2.5, lty = 2)
legend("topleft", legend = "single PCR", bty = "n")
legend("topright", legend = "double PCR", bty = "n")
box(lwd = 2)
# dev.off()




########################################################################################################
# GENERATE A STRIPCHART
########################################################################################################
### FOR PDF
# pdf(file = "BC_dis_stripchart.pdf")

### FOR EPS
# setEPS()
# postscript("BC_dis_stripchart.eps")


# jittering is random; to control it, set a seed number
set.seed(1)
par(mar = c(5,4,1,1))
stripchart(
  x = plot_input, #c(PCR_single, PCR_double)
  vertical = TRUE,
  method = "jitter",
  jitter = 0.2,
  pch = 1, col = "black", #bg = "grey",
  cex = 0.8,
  # las = 2,
  ylab = "Mean Bray-Curtis Dissimilarity",
  ylim = c(0,1),
  xaxt = "n"
  # group.names = c("within tags,\nsingle PCR", "between tags,\nsingle PCR", "within tags,\ndouble PCR", "between tags,\ndouble PCR")
)
axis(
  side = 1,
  at = 1:4,
  tick = FALSE,
  line = 1,
  labels = rep(c("within primer\nindex", "among primer\nindex"), 2)
    # "within tags,\nsingle PCR", "between tags,\nsingle PCR", "within tags,\ndouble PCR", "between tags,\ndouble PCR"
    # paste("within tags,\namong libraries\nN = ", length(BC_within_tag), sep = ""),
    # paste("among tags,\namong libraries\nN = ", length(BC_between_tag_mean), sep = "")
  # )
)
abline(v = 2.5, lty = 2)
legend("topleft", legend = "single PCR", bty = "n")
# legend(x = 1, y = 1, legend = "single PCR", bty = "n")
legend("topright", legend = "double PCR", bty = "n")
# legend(x = 3, y = 1, legend = "double PCR", bty = "n")
box(lwd = 2)
# dev.off()
########################################################################################################






################################################################################
# SEQUENCE-SPECIFIC EFFECTS 
################################################################################
################################################################################
# Make some plots for a given OTU (referenced by number)
# You should check to make sure these are the same sequences (They are in the case of the current dataset)
target_OTU <- c(1,1)

# number of OTUs to isolate (e.g. the top 10, top 20, top 50)
n_OTUs <- 10

# referencing across lists of dataframes was driving me nuts. Instead...
# Create a list of data frames which combine metadata and OTU abundance data for top 10 OTUs
DAT_comb <- list()
for(i in 1:length(DAT_OTU)){
  DAT_comb[[i]] <- data.frame(DAT_pool[[i]], DAT_OTU[[i]][DAT_pool[[i]]$lib_tag,1:n_OTUs])
}
names(DAT_comb) <- c("Single PCR", "Double PCR")

# translate target OTU numbers into new column numbers
OTU_col <- numeric()
for(i in 1:length(DAT_OTU)){
  OTU_col[i] <- ncol(DAT_pool[[i]]) + target_OTU[i]
}

# isolate the proportion data
# OTU_props <- lapply(DAT_OTU, function(x) x[,target_OTU])
OTU_props <- numeric()
for(i in 1:length(DAT_OTU)){
  OTU_props <- c(OTU_props, DAT_comb[[i]][,OTU_col[i]])
}

# calculate the range in abundance
range_OTU <- range(OTU_props)

# a function to get midpoints
midPoints <- function(x){
  (x[-length(x)]+x[-1])/2
}

################################################################################
# GENERATE A PLOT
################################################################################

# FOR PDF
# pdf(file = "OTU_abundance.pdf", height = 5, width = 10)

### FOR EPS
# setEPS()
# postscript(file = "OTU_abundance.eps", height = 5, width = 10)

par(mar = c(5,4,1,1))
par(mfrow = c(1,2))
for(i in 1:length(DAT_comb)){
  
  # only include desired sample type
  plot_dat <- droplevels(DAT_comb[[i]][DAT_comb[[i]]$sample_type == sample_type_include,])
  
  # split the data frame by sample name, include only columns for target OTU abundance and tag sequence (reduced to numeric factor levels)
  TEMP <- lapply(split(plot_dat[,c(OTU_col[i], which(colnames(plot_dat)  == "tag_sequence"))], plot_dat$sample_name), function(x) {x$tag_sequence <- as.numeric(droplevels(x$tag_sequence)); return(x)})
  
  # bind those into a single data frame
  TEMP2 <- do.call(rbind, TEMP)
  
  # set vertical lines between each sample
  TEMP_L <- c(0,sapply(TEMP, nrow))
  v_lines <- cumsum(TEMP_L/2)+0.5
  
  # put ticks at the midpoints between vertical lines
  ticks <- (v_lines[-length(v_lines)] + v_lines[-1])/2
  
  plot(
    TEMP2[(TEMP2$tag_sequence == 1),1], 
    ylab = "proportional abundance", 
    xaxt = "n", 
    xlab = "Environmental Sample", 
    ylim = c(0,1)
  )
  points(TEMP2[(TEMP2$tag_sequence == 2),1], pch = 19) # pch = 19 for filled circles
  abline(v = v_lines[-1], lty = 2, col = "grey")
  axis(side = 1, at = ticks, labels = c(1:length(ticks)))
  legend("topright", legend = c("index 1", "index 2"), pch = c(1,19), bg = "white",  box.lty = 0)
  legend("topleft", legend = names(DAT_comb)[[i]], bg = "white", box.lty = 0)
  box(lwd = 2)
  
  # what is the range in sequence abundance among indexes for this OTU?
  print(TEMP)
  print(range_index <- sapply(TEMP, function(x) diff(range(x[,1]))))
  
  # what is the maximum difference in sequence abundance among indexes for this OTU
  print(max_diff <- max(range_index))
}
# dev.off()

# OTU3: 0.7188945 to 0.03476546

# OTU2: 0.7709074 to 0.01140888

# OTU1: 0.6496404 to 0.03382474








################################################################################
# GENERATE A STRIPCHART
################################################################################

# Initialize a list to store values by sample/index
OTU_prop_index <- list()

for(i in 1:length(DAT_pool)){
  # sort by sample name
  # TEMP <- DAT_pool[[1]][order(DAT_pool[[1]]$sample_name),]
  
  # only include target samples (e.g. environmental)
  TEMP <- droplevels(DAT_pool[[i]][DAT_pool[[i]]$sample_type == sample_type_include,])
  
  TEMP_lib_tag <- lapply(
    split(TEMP, TEMP$sample_name), 
    function(x) split(x$lib_tag, droplevels(x$tag_sequence))
  )
  
  index1 <- lapply(TEMP_lib_tag, function(x) x[[1]])
  index2 <- lapply(TEMP_lib_tag, function(x) x[[2]])
  
  OTU_prop_index[[i]] <- list(
    lapply(index1, function(x) OTU_props[[i]][x]),
    lapply(index2, function(x) OTU_props[[i]][x])
  )
  
}

### FOR PDF
# pdf(file = "OTU_abundance_stripchart.pdf")

### FOR EPS
# setEPS()
# postscript("OTU_abundance_stripchart.eps")

plot(OTU_prop_index[[1]][[2]])
length(OTU_prop_index)

# jittering is random; to control it, set a seed number
set.seed(1)
par(mar = c(5,4,1,1))
plot_OTU <- function(X, PTS = 1, ADD = FALSE){
  stripchart(
    x = X,
    vertical = TRUE,
    method = "jitter",
    jitter = 0.2,
    pch = PTS, 
    col = "black", #bg = "grey",
    cex = 1,
    # las = 2,
    ylab = "Proportional abundance in sample",
    ylim = c(0,1),
    # xaxt = "n",
    add = ADD, 
    group.names = names(OTU_prop_index[[2]][[1]])
  )
}
plot_OTU(X = DAT_comb[[1]][,OTU_col[1]])
plot_OTU(OTU_prop_index[[1]][[1]])
plot_OTU(OTU_prop_index[[1]][[2]], PTS = 16, ADD = TRUE)

abline(v = c(1:length(OTU_prop_index[[1]][[2]]))*2, lty = 2, col = "grey")

axis(
  side = 1,
  at = 1:4,
  tick = FALSE,
  line = 1,
  labels = rep(c("within primer\nindex", "among primer\nindex"), 2)
  # "within tags,\nsingle PCR", "between tags,\nsingle PCR", "within tags,\ndouble PCR", "between tags,\ndouble PCR"
  # paste("within tags,\namong libraries\nN = ", length(BC_within_tag), sep = ""),
  # paste("among tags,\namong libraries\nN = ", length(BC_between_tag_mean), sep = "")
  # )
)
abline(v = 2.5, lty = 2)
legend("topleft", legend = "single PCR", bty = "n")
# legend(x = 1, y = 1, legend = "single PCR", bty = "n")
legend("topright", legend = "double PCR", bty = "n")
# legend(x = 3, y = 1, legend = "double PCR", bty = "n")
box(lwd = 2)
dev.off()



########################################################################################################

# ALTERNATE PLOTTING
# plot(
#   do.call(c, OTU_prop_index[[2]][[1]]),
#   ylim = c(0,1)
# )
# points(
#   do.call(c, OTU_prop_index[[2]][[2]]),
#   pch = 16
# )

##############################################################################################
# GRAVEYARD

# ?unsplit(dist_reform, f = names(dist_reform))

# potentially useful:
within_tag <- function(x){
  which(x[,"tag1"] == x[,"tag2"])
}

between_tag <- function(x){
  which(x[,"tag1"] != x[,"tag2"])
}

which(TEMP_mat[,"lib1"] != TEMP_mat[,"lib2"])
which(TEMP_mat[,"tag1"] == TEMP_mat[,"tag2"])
which(TEMP_mat[,"lib1"] == TEMP_mat[,"lib2"])

################################################################################
stripchart(
  x = split(DAT_OTU[[i]][,5], tag_to_samplename),
  vertical = FALSE,
  method = "jitter",
  pch = 21, col = "black",
  las = 1,
  xlab = "Number of Clusters"
  # xlim = c(0,3000)
  # group.names=as.character(reps)
)


