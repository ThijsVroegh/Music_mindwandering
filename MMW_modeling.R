#    A temporal network approach to music-induced mind wandering   
#                      L.Taruffi & T.Vroegh 2024
#
# 1 Load libraries ----
library(tidyverse)
library(qgraph) 
library(psychonetrics) 
library(bootnet)      
library(NetworkComparisonTest)

no_cores <- parallel::detectCores() - 1

# remove all data in global environment
rm(list = ls())

# source functions
source("sources/plot_hist_facet.R")
source("sources/plot_bar_facet.R")
source("sources/generate_ri_clpm.R")

# random permutation seed
set.seed(1234) 

# read in cleaned data
mydata <- readRDS("rds/data_cleaned.rds")

## correlations
mydata_t0 <- mydata %>% select(t0_AF:t0_AE) %>% rename_at(vars(starts_with("t0_")),  ~ str_replace(., "t0_", ""))
mydata_t1 <- mydata %>% select(t1_AF:t1_AE) %>% rename_at(vars(starts_with("t1_")),~ str_replace(., "t1_", ""))
mydata_t2 <- mydata %>% select(t2_AF:t2_AE) %>% rename_at(vars(starts_with("t2_")),~ str_replace(., "t2_", ""))

# 2 GGM networks reading versus music ----

# Create a new dataframe for music with the mean values of the two music networks
df_mean <- data.frame(matrix(0, nrow = nrow(mydata_t1), ncol = ncol(mydata_t1)))

# Loop through each variable
for (col_index in seq_along(colnames(mydata_t1))) {
  df_mean[[col_index]] <- rowMeans(cbind(mydata_t1[[col_index]], mydata_t2[[col_index]]))
}

# Rename the columns based on original variable names
colnames(df_mean) <- colnames(mydata_t1)

reading <- mydata_t0 
music   <- df_mean   

network_reading <- EBICglasso(cor_auto(reading),n = nrow(reading)) 
network_music <- EBICglasso(cor_auto(music),n = nrow(music)) 

t_reading <- qgraph(network_reading, layout = 'spring', 
             details = FALSE, labels=colnames(network_reading),
             vsize=8, border.width=2, edge.labels=FALSE, DoNotPlot = TRUE)

t_music <- qgraph(network_music, layout = 'spring', 
             details = FALSE, labels=colnames(network_music),
             vsize=8, border.width=2, edge.labels=FALSE, DoNotPlot = TRUE)

r_m <- averageLayout(t_reading, t_music)

t_reading <- qgraph(network_reading, layout = L, details = FALSE, 
             labels=colnames(network_reading), vsize=12, border.width=1, 
             edge.labels=FALSE, title="Reading")

t_music <- qgraph(network_music, layout = L, details = FALSE, 
             labels=colnames(network_music), vsize=12, border.width=1, 
             edge.labels=FALSE, title="Music listening")

# Plot networks
tiff("cross sectional networks reading and music.tiff", width = 3000, height = 2000, units = "px", res = 300)
par(mfrow = c(1, 2))
plot(t_reading)
plot(t_music) 
dev.off()


## Three networks ----
network_t0 <- EBICglasso(cor_auto(mydata_t0),n = nrow(mydata_t0)) 
network_t1 <- EBICglasso(cor_auto(mydata_t1),n = nrow(mydata_t1)) 
network_t2 <- EBICglasso(cor_auto(mydata_t2),n = nrow(mydata_t2)) 

par(mfrow = c(1, 1))
t0 <- qgraph(network_t0, layout = 'spring', 
             details = FALSE, labels=colnames(network_t0),
             vsize=8, border.width=2, edge.labels=FALSE,DoNotPlot = TRUE)

t1 <- qgraph(network_t1, layout = 'spring', 
             details = FALSE, labels=colnames(network_t1),
             vsize=8, border.width=2, edge.labels=FALSE,DoNotPlot = TRUE)

t2 <- qgraph(network_t2, layout = 'spring', 
             details = FALSE, labels=colnames(network_t2),
             vsize=8, border.width=2, edge.labels=FALSE,DoNotPlot = TRUE)

# Create average layout to be used consistently across all networks
L <- averageLayout(t0,t1,t2)

t0 <- qgraph(network_t0, layout = L, details = FALSE, 
             labels=colnames(network_t0), vsize=12, border.width=1, 
             edge.labels=FALSE,title="Reading")

t1 <- qgraph(network_t1, layout = L, details = FALSE, 
             labels=colnames(network_t1), vsize=12, border.width=1, 
             edge.labels=FALSE,title="Halfway during music listening")

t2 <- qgraph(network_t2, layout = L, details = FALSE,
             labels=colnames(network_t2), vsize=12, border.width=1, 
             edge.labels=FALSE, title = "Directly after music listening")

# Plot networks
tiff("three network plots reading and music.tiff", width = 2200, height = 2200, units = "px", res = 300)
par(mfrow = c(2,2))
plot(t0)
plot(t1) 
plot(t2)
dev.off()

## Centralities of 3 networks ----

#expected influence as main measure of interest
cent0 <- centrality(t0)$InExpectedInfluence
cent1 <- centrality(t1)$InExpectedInfluence
cent2 <- centrality(t2)$InExpectedInfluence
met   <- cbind(cent0,cent1,cent2)

#plotting
means <- data.frame(met, stringsAsFactors = FALSE) %>%
  tibble::rownames_to_column() %>%
  dplyr::rename(Centrality = rowname) %>%
  reshape2::melt(id.vars = "Centrality", variable.name = "Time", value.name = "Mean")

means$Mean <- as.numeric(means$Mean)

gg_centralities <- means %>%
  ggplot(aes(Time, Mean, group = Centrality, color = Centrality)) +
  ylim(0, 1.5) +
  geom_point() + 
  geom_line() + 
  geom_text(aes(label = Centrality, hjust = 0.5, vjust = 0.5)) +
  labs(x = "Session", y = "z-score") + 
  scale_color_viridis_d(labels = colnames(network_t0)) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot centrality for three measurement points
tiff("centralities through time from reading to music.tiff", width = 2200, height = 1500, units = "px", res = 300)
par(mfrow = c(1,1))
plot(gg_centralities)
dev.off()

##  Estimate stability and accuracy

## Set seed
set.seed(1234)

# music
result_reading <- estimateNetwork(reading, default = "EBICglasso", corMethod = "cor_auto",tuning = 0.5)

## estimate the accuracy of the edge weights in the network
boot_reading <- bootnet(result_reading, nBoots = 5000, nCores = no_cores)

# Print edge weight CI
tiff("Fig S3 bootstrap_edge weight_accuracy_reading.tiff", width = 2200, 
     height = 2200, units = "px", res = 300)
par(mfrow = c(1, 1))
plot(boot_reading, plot = "interval", order = "sample", labels = TRUE)
dev.off()

tiff("Fig S4 bootstrap_edge weight_diff_reading.tiff", width = 2200, 
     height = 2200, units = "px", res = 300)
par(mfrow = c(1, 1))
plot(boot_reading,"edge", plot = "difference", onlyNonZero = TRUE, order = "sample")
dev.off()

boot_reading_case <- bootnet(result_reading, nBoots = 5000, type = "case", nCores = no_cores, 
                 statistics = c("strength"))

# CS-coefficients for strength in the person-dropping stability analysis
corStability(boot_reading_case, statistics = c("strength"))

tiff("Fig S5 bootstrap_centrality_stability_reading.tiff", width = 2200, 
     height = 2200, units = "px", res = 300)
par(mfrow = c(1, 1))
  plot(boot_reading_case, statistics = c("strength"))
dev.off()

## music
result_music <- estimateNetwork(music, default = "EBICglasso", corMethod = "cor_auto",tuning = 0.5)

## estimate the accuracy of the edge weights in the network
boot_music <- bootnet(result_music, nBoots = 5000, nCores = no_cores)

tiff("Fig S6 bootstrap_edge weight_accuracy_music.tiff", width = 2200, 
     height = 2200, units = "px", res = 300)
par(mfrow = c(1, 1))
plot(boot_music, plot = "interval", order = "sample", labels = TRUE)
dev.off()

tiff("Fig S7 bootstrap_edge weight_diff_music.tiff", width = 2200, 
     height = 2200, units = "px", res = 300)
par(mfrow = c(1, 1))
plot(boot_music,"edge", plot = "difference", onlyNonZero = TRUE, order = "sample")
dev.off()

boot_music_case <- bootnet(result_music, nBoots = 5000, type = "case", 
                           nCores = no_cores, statistics = c("strength"))

# CS-coefficients for strength in the person-dropping stability analysis
corStability(boot_music_case, statistics = c("strength"))

tiff("Fig S8 bootstrap_centrality_stability_music.tiff", width = 2200, 
     height = 2200, units = "px", res = 300)
par(mfrow = c(1, 1))
plot(boot_music_case, statistics = c("strength"))
dev.off()

## Network Comparison Test ----

#recalculate networks with bootnet; output is input for the NCT
network_reading <- estimateNetwork(reading, default = "EBICglasso", corMethod = "spearman")
network_music   <- estimateNetwork(music  , default = "EBICglasso", corMethod = "spearman")

set.seed(1234)
NCT_reading_music <- NCT(network_reading, network_music, 
                         it               = 1000,
                         paired           = TRUE, 
                         abs              = FALSE, 
                         test.edges       = TRUE,
                         edges            = "all", 
                         test.centrality  = TRUE, # investigate centrality differences
                         progressbar      = TRUE,
                         p.adjust.methods = c("BH"), # Holm-Bonferroni correction
                         centrality       = c("expectedInfluence"), 
                         nodes            = "all",
                         verbose          = FALSE)

summary(NCT_reading_music)

# Checking which edges are different: print only edge differences with p < .05

difference_value <- function(NCT, alpha = 0.05){
  
  diff_edges <- NCT$einv.pvals %>% dplyr::filter(`p-value` <= alpha)
  
  for (i in 1:nrow(diff_edges)) {
    var_1 <- as.character(diff_edges[i, 1])
    var_2 <- as.character(diff_edges[i, 2])
    
    value_net_1 <- NCT$nw1[var_1, var_2]
    value_net_2 <- NCT$nw2[var_1, var_2]
    
    abs_difference <- abs(value_net_1 - value_net_2)
    p_value <- diff_edges$`p-value`[i]
    
    cat("Test Edge", i, "\n----\n")
    cat(var_1, "and", var_2)
    cat("\nNetwork 1:", value_net_1,
        "\nNetwork 2:", value_net_2)
    cat("\nAbsolute difference:", abs_difference,
        "with p-value =", p_value, "\n----\n")
  }
}

difference_value(NCT_reading_music)

# 3 PanelGVAR modeling ----

## a) Checks ----
# check the variances across the measurement points
mydata %>% select(t0_AF:t2_AE) %>% 
  map_df(sd) %>% 
  pivot_longer(starts_with('t'), 
               names_to = c('time', '.value'), 
               names_sep = '\\_')

# check the means across the measurement points
mydata %>% select(t0_AF:t2_AE) %>% 
  map_df(mean) %>% 
  pivot_longer(starts_with('t'), 
               names_to = c('time', '.value'), 
               names_sep = '\\_')

## b) Scaling/ standardizing data ----
data <- mydata %>% select(t0_AF:t2_AE)

# standardize data across time points
data_scaled <- as.data.frame(scale(data))

# plot distributions of scaled variables
plot_hist_facet(data_scaled, bins = 8, ncol = 7)

saveRDS(data_scaled, file = "rds/data_scaled.rds")
# data_scaled <- readRDS("rds/rds/data_scaled.rds")

# Checking similarity of sd's
data_scaled %>% 
  map_df(sd) %>% 
  pivot_longer(starts_with('t'), 
               names_to = c('time', '.value'), 
               names_sep = '\\_')

# Checking similarity of means
data_scaled %>% 
  map_df(mean) %>% 
  pivot_longer(starts_with('t'), 
               names_to = c('time', '.value'), 
               names_sep = '\\_')

## c) Standardizing and detrending data ----
# Since GVAR models assume stationary relations across time, prior to fitting
# the models, data was detrended for linear time-related effects and was then
# standardized across time points
 vars <- c("AF","SA","IM", "VA","CALM","THO","AE")
# 
 variables        <- list()
 variables_lm     <- list()
 variables_scaled <- list()
# 
 mydata_detrend <- mydata %>%
   rename(id = ID) %>%
   select(id,t0_AF:t2_AE) %>%
   mutate(t0_VA   = log10(max(t0_VA   + 1) - t0_VA),
          t1_VA   = log10(max(t1_VA   + 1) - t1_VA),
          t2_VA   = log10(max(t2_VA   + 1) - t2_VA),
          t0_CALM = log10(max(t0_CALM + 1) - t0_CALM),
          t1_CALM = log10(max(t1_CALM + 1) - t1_CALM),
          t2_CALM = log10(max(t2_CALM + 1) - t2_CALM))
# 
# detrending loop 
 for (i in 1:7) {
#
   # reshape data
   variables[[i]] <- mydata_detrend %>% 
                     select(contains(vars[i]), id) %>%
                     gather(time, var, contains(vars[i])) %>% 
                     mutate(time = rep(c(1,2,3), each = 352),
                     # dummy variable for reading and music
                     read_music = case_when(
                     time == 1 ~ 1,
                     TRUE      ~ 0),
                     read_music = as.factor(read_music))
   
   # detrend linearly and quadratically
   variables_lm[[i]] <- lm(var ~ read_music + time + I(time^2), data = variables[[i]])
   
   # save detrended data
   variables[[i]]$var[!is.na(variables[[i]]$var)] <- residuals(variables_lm[[i]])
   
   # reshape data
   variables_scaled[[i]] <-   variables[[i]] %>% select(-read_music) %>% 
     spread(time, var) %>%
     select(-id) %>%
     as.matrix %>%
     as.vector() %>%
     scale() %>%
     matrix(nrow = 352, ncol = 3) %>% # ncol refers to number of waves
     as.data.frame()
   
 # save formatted and detrended data
    colnames(variables_scaled[[i]]) <- mydata_detrend %>% 
     select(contains(vars[i])) %>%
     colnames
 }
 
# # properly reorder variables
 data_detrended <- variables_scaled %>% 
   as.data.frame() %>% 
   select(contains("t0"),contains("t1"),contains("t2"))

# plot distributions of scaled and de-trended variables
 plot_hist_facet(data_detrended, bins = 8, ncol = 7)

## c) Define design matrix ----
design <- matrix(colnames(data_scaled), nrow = 7, ncol = 3)
colnames(design) <- c("t0", "t1", "t2")
rownames(design) <- c("AF","SA","IM","VA","CALM", "THO","AE")

layout(matrix(1:4,2,2))
qgraph(cor(data_scaled[,design[,1]]),maximum=1,title="cors t0", theme = "colorblind",vsize = 20)
qgraph(cor(data_scaled[,design[,2]]),maximum=1,title="cors t1", theme = "colorblind",vsize = 20)
qgraph(cor(data_scaled[,design[,3]]),maximum=1,title="cors t2", theme = "colorblind",vsize = 20)

ev1 <- eigen(cor(data_scaled[,design[,1]]))$values
ev2 <- eigen(cor(data_scaled[,design[,2]]))$values
ev3 <- eigen(cor(data_scaled[,design[,3]]))$values

# It appears that the second and third time points (music) are
# more unidimensional than the first time point (reading)
# Plot networks
tiff("Eigenvalues.tiff", width = 2200, height = 2200, units = "px", res = 300)
  par(mfrow = c(1, 1))
  matplot(cbind(ev1,ev2,ev3),type = "l", ylab = "Eigenvalue")
dev.off()

# Labels to be used in graphs
labels <- rownames(design)

## d) Estimate saturated model----
model1 <- panelgvar(data      = data_scaled, 
                    vars      = design,
                    estimator = "ML",
                    storedata = TRUE)

## ) Run the model
model1 <- model1 %>% runmodel

# The warning message (from package creator: https://github.com/SachaEpskamp/psychonetrics/issues/10)
# The warning “ The optimizer encountered at least one non-positive definite matrix
# and used a pseudoinverse in parameter estimation. Results may not be accurate.”
# can safely be ignored if the parameters look ok (e.g. no partial correlations of 1/-1).
# I would use the nlminb optimizer, it works the best.

# inspect model fit
model1 %>% fit 
               
# check parameters
model1 %>% parameters()

## Plot analytic confidence intervals (for the saturated model)
# contemporaneous
tiff(filename = "Final_CIplots_within_saturated_model.tiff", width = 4000, height = 4000, res = 450)
par(mfrow = c(1,1))
CIplot(model1, "omega_zeta_within")
dev.off()

# temporal
tiff(filename = "Final_CIplots_beta_saturated_model.tiff",width = 4000, height = 4000, res = 450)
par(mfrow = c(1,1))
CIplot(model1, "beta")
dev.off()


## e) Prune to find a sparse model ----
model2 <- model1 %>% prune(alpha = 0.05, recursive = FALSE)
fit(model2) 
            
## f) Stepup model ----
model3 <- model2 %>% stepup(criterion = "bic",alpha = 0.05)
fit(model3) 
            
model3 %>% parameters()

## g) Compare all models ----
compare(baseline     = model1, 
        pruned       = model2,
        stepup       = model3) # <- best AIC and BIC

# 4 Plotting output ----
# Note that the PDC and the beta matrix should align but are expressed in 
# opposite directions
temporal        <- getmatrix(model3,"PDC") # equals t(getmatrix(model2, "beta"))
contemporaneous <- getmatrix(model3,"omega_zeta_within")
between         <- getmatrix(model3,"omega_zeta_between")

graph1 <- qgraph(temporal,layout = "spring", DoNotPlot = TRUE )
graph2 <- qgraph(contemporaneous,layout = "spring", DoNotPlot = TRUE )
graph3 <- qgraph(between,layout = "spring", DoNotPlot = TRUE )

# Create average layout to be used consistently across all networks
L <- averageLayout(graph1,graph2, graph3)

## a) Community analysis ----
set.seed(1234)
network_t0 <- EBICglasso(cor_auto(mydata_t0),n = nrow(mydata_t0)) 
network_t1 <- EBICglasso(cor_auto(mydata_t1),n = nrow(mydata_t1)) 
network_t2 <- EBICglasso(cor_auto(mydata_t2),n = nrow(mydata_t2)) 

t0 <- qgraph(network_t0, layout = 'spring', 
             details = FALSE, labels = colnames(network_t0),
             vsize=8, border.width=2, edge.labels=FALSE)

t1 <- qgraph(network_t1, layout = 'spring', 
             details = FALSE, labels=colnames(network_t1),
             vsize=8, border.width=2, edge.labels=FALSE)

t2 <- qgraph(network_t2, layout = 'spring', 
             details = FALSE, labels=colnames(network_t2),
             vsize=8, border.width=2, edge.labels=FALSE)

net0 <- igraph::as.igraph(t0, attributes=TRUE)
net1 <- igraph::as.igraph(t1, attributes=TRUE)
net2 <- igraph::as.igraph(t2, attributes=TRUE)

matrix_spinglass <- matrix(NA, nrow = 1,ncol = 250)

for (i in 1:250) {
  set.seed(i)
  spinglass <- igraph::spinglass.community(net2, weights=NULL, vertex=NULL, spins=10,
                                           parupdate=FALSE, start.temp=1, stop.temp=0.01,
                                           cool.fact=0.99, update.rule="simple", gamma=0.5, 
                                           implementation="neg",gamma.minus=0)
  matrix_spinglass[1,i] <- max(spinglass$membership) 
}

mean(as.vector(matrix_spinglass)) 
max(as.vector(matrix_spinglass))  
min(as.vector(matrix_spinglass))    

net0_groups <- igraph::spinglass.community(net0) 
net1_groups <- igraph::spinglass.community(net1) 
net2_groups <- igraph::spinglass.community(net2)

# 3 communities identified in three measurements
gr_net0 <- list('Gr1'= c(4,5),   'Gr2'=c(1,6), 'Gr3'=c(2,3,7))
gr_net1 <- list('Gr1'= c(2,4,5), 'Gr2'=c(1,6), 'Gr3'=c(3,7))
gr_net2 <- list('Gr1'= c(4,5),   'Gr2'=c(1,6), 'Gr3'=c(2,3,7))

gr <- list('Emotional State' = c(4,5), 'Focus' = c(1,6), 'Altered Experience' = c(2,3,7))

names <- c("AF","SA","IM","VA","CALM", "THO","AE")

## b) Three networks in one plot ----
tiff(filename = "Contemporaneous and temporal networks with panelgvar.tiff",
     width  = 4500, 
     height = 2500, 
     res    = 600)

par(mfrow = c(1, 2))

cg <- qgraph(contemporaneous, 
             layout = L,
             title="Contemporaneous network",
             theme='colorblind', 
             negDashed=FALSE,
             legend=FALSE, 
             labels = names,
             nodeNames = names, 
             groups = gr,
             vsize=10,
             parallelEdge = TRUE,
             color=viridis::viridis_pal()(5)[3:5],
             details = F)

tg <- qgraph(temporal, 
             layout = L,
             title="Temporal network", 
             theme='colorblind', 
             negDashed=FALSE, 
             legend.cex=0.5, 
             legend=F, 
             labels = names,
             nodeNames = names,
             groups=gr,
             vsize=10, 
             parallelEdge = TRUE,
             color=viridis::viridis_pal()(5)[3:5], 
             asize=6,
             curve=0.75, 
             curveAll=F, 
             details = F)
dev.off()

## c) Strength centrality plots ----
par(mfrow=c(1,1))


tiff(filename="strength centrality contemporaneous network.tiff", width=1250, height=2500, res=450)
centralityPlot(cg,labels = vars, scale = "z-score", include = c("Strength"))
dev.off()

tiff(filename="strength centrality temporal network.tiff", width=1250, height=2500, res=450)
centralityPlot(tg,labels = vars,scale = "z-score", include = c("InStrength", "OutStrength"))
dev.off()

# 5. RI-CLPM -----                        

# This code was adapted from the example script provided by Freichel et al. (2023)
#Cross-Lagged Panel Models for Studying Psychopathology: A Comparative Overview of
#Structural Equation and Panel Network Approaches (tinyurl.com/4m7m78sm)

# remove thought variable because of non PD
no_thought <- data_scaled %>% select(-t0_THO,-t1_THO,-t2_THO)

design <- matrix(colnames(no_thought), nrow = 6, ncol = 3)
colnames(design) <- c("t0", "t1", "t2")
rownames(design) <- c("AF","SA","IM","VA","CALM","AE")


# use existing function for model estimation
riclpm_res <- generate_ri_clpm(thijs, design)

riclpm_summary <- summary(riclpm_res$lavres, 
                          standardized = TRUE, 
                          fit.measures = TRUE,
                          modindices   = TRUE,
                          rsquare      = TRUE)

lavaan::standardizedSolution(riclpm_res$lavres)

fitmeasures <- print(lavaan::fitMeasures(riclpm_res$lavres, 
                                         c("chisq", "df","pvalue", 
                                           "cfi","tli", "srmr","gfi",
                                           "rmsea","rmsea.ci.lower",
                                           "rmsea.ci.upper","rmsea.pvalue"), 
                                         output = "text"), add.h0 = TRUE)

# inspect model fit 
riclpm_summary$fit 

## a) Temporal network ----
temporal_thresholded_riclpm <- riclpm_res$matrices

## threshold temporal effects from RI-CLPM
temporal_thresholded_riclpm$PDC$est[temporal_thresholded_riclpm$PDC$p >= 0.001] <- 0

# visualized thresholded temporal network
temp_thresh_riclpm <- qgraph(temporal_thresholded_riclpm$PDC$est, 
                             #labels = labels_net, 
                             groups = gr, 
                             nodeNames = names,
                             theme = "colorblind", 
                             layout = L,
                             #vsize = 6,legend.cex = .50, legend.mode = "style1", # aesthetics,
                             legend = FALSE
)

tiff(filename = "temporal_thresh_riclpm.tiff", width = 3000, height = 2000, res = 300)
plot(temp_thresh_riclpm)
dev.off()

## b) Contemporaneous Network ---- 
contemp_thresholded_riclpm <- riclpm_res$matrices$contemporaneous_covariances

## threshold contemporaneous effects from RI-CLPM
contemp_thresholded_riclpm$est[contemp_thresholded_riclpm$p >= 0.001] <- 0
diag(contemp_thresholded_riclpm$est) <- 0

# visualized thresholded contemporaneous network
contemp_thresh_riclpm <- qgraph(contemp_thresholded_riclpm$est, 
                                groups = gr, 
                                nodeNames = names,
                                theme = "colorblind", 
                                layout = L,
                                #vsize = 6,legend.cex = .50, legend.mode = "style1" # aesthetics,
                                legend = FALSE
)

tiff(filename = "contemporaneous_thresh_riclpm.tiff", width = 3000, height = 2000, res=300)
plot(contemp_thresh_riclpm)
dev.off()


# 6 Stability Analysis with bootstrapping ----

# see Nur Hani Zainal & Michelle G. Newman" (2021) for a detailed description

## a) Bootstrapping ----
# keep number low (ca. 200), otherwise may take quite a long time
set.seed(1234)
nBoot <- 200

Bootstraps <- lapply(1:nBoot, function(x) {
  
  message("Simulation: ",x)
  
  # Sample from the data
  bootData <- data_scaled[sample(1:nrow(data_scaled), 
                                 round(0.75 * nrow(data_scaled))), ]
  
  # Form model
  bootstrapped_model <-
    panelgvar(
      bootData,
      vars      = design,
      estimator = "ML")  
  
  # Run first time
  bootstrapped_model <- bootstrapped_model %>% runmodel %>% 
  prune(alpha = 0.05, recursive = FALSE) %>% 
    stepup(criterion = "bic")
  
  return(bootstrapped_model)
})

# Extracting stability analyses to check the degree to which individual edges
# were included across all bootstrap samples 

## b) Bootstrapped results for temporal network ----
# check if individual edges are included (y/n)
resBoots_temp <-
  lapply(Bootstraps, function(x)
    ifelse(getmatrix(x, "PDC") > 0, 1, 0))

# how often (%) is edge included, over all bootstrap iterations
apply(simplify2array(resBoots_temp), 1:2, mean) 
Bootstraps_temporal_est_df            <- as.data.frame(apply(simplify2array(resBoots_temp), 1:2, sum))
row.names(Bootstraps_temporal_est_df) <- colnames(Bootstraps_temporal_est_df) <- colnames(mydata_t0)

write.csv(Bootstraps_temporal_est_df, "Bootstraps_temporal_est_df.csv")

## c) Bootstrapped results for contemporaneous network----
resBoots_cont <-
  lapply(Bootstraps, function(x)
    ifelse(getmatrix(x, "omega_zeta_within") > 0, 1, 0))

# how often (%) is edge included, over all bootstrap iterations
apply(simplify2array(resBoots_cont), 1:2, mean) 
Bootstraps_contemporaneous_est_df            <- as.data.frame(apply(simplify2array(resBoots_cont), 1:2, sum))
row.names(Bootstraps_contemporaneous_est_df) <- colnames(Bootstraps_contemporaneous_est_df) <- colnames(mydata_t0)

write.csv(Bootstraps_contemporaneous_est_df, "Bootstraps_contemporaneous_est_df.csv")

sessioninfo::session_info()
