#    A temporal network approach to music-induced mind wandering   
#                      L.Taruffi & T.Vroegh 2024
#
# 1 Load libraries ----
library(tidyverse)
library(qgraph) 
library(psychonetrics) 
library(bootnet)      
library(NetworkComparisonTest)

# remove all data in global environment
rm(list = ls())

# source functions
source("sources/plot_hist_facet.R")

# read in cleaned data
mydata <- readRDS("data_cleaned.rds")

# Modeling

## correlations ----
mydata_t0 <- mydata %>% select(t0_AF:t0_AE) %>% rename_at(vars(starts_with("t0_")),  ~ str_replace(., "t0_", ""))
mydata_t1 <- mydata %>% select(t1_AF:t1_AE) %>% rename_at(vars(starts_with("t1_")),~ str_replace(., "t1_", ""))
mydata_t2 <- mydata %>% select(t2_AF:t2_AE) %>% rename_at(vars(starts_with("t2_")),~ str_replace(., "t2_", ""))

# 2 GGM networks reading versus music ----

## Networks ----
network_t0 <- EBICglasso(cor_auto(mydata_t0),n = nrow(mydata_t0)) 
network_t1 <- EBICglasso(cor_auto(mydata_t1),n = nrow(mydata_t1)) 
network_t2 <- EBICglasso(cor_auto(mydata_t2),n = nrow(mydata_t2)) 

par(mfrow = c(1, 1))
t0 <- qgraph(network_t0, layout = 'spring', 
             details = FALSE, labels=colnames(network_t0),
             vsize=8, border.width=2, edge.labels=FALSE)

t1 <- qgraph(network_t1, layout = 'spring', 
             details = FALSE, labels=colnames(network_t1),
             vsize=8, border.width=2, edge.labels=FALSE)

t2 <- qgraph(network_t2, layout = 'spring', 
             details = FALSE, labels=colnames(network_t2),
             vsize=8, border.width=2, edge.labels=FALSE)

# Create average layout to be used consistently across all networks
L <- averageLayout(t0,t1,t2)

t0 <- qgraph(network_t0, layout = L, details = FALSE, 
             labels=colnames(network_t0), vsize=12, border.width=1, 
             edge.labels=FALSE,
             title="Reading")

t1 <- qgraph(network_t1, layout = L, details = FALSE, 
             labels=colnames(network_t1), vsize=12, border.width=1, 
             edge.labels=FALSE,
             title="Halfway during music listening")

t2 <- qgraph(network_t2, layout = L, details = FALSE,
             labels=colnames(network_t2), vsize=12, border.width=1, 
             edge.labels=FALSE,
             title = "Directly after music listening")

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
  melt(id.vars = "Centrality", variable.name = "Time", value.name = "Mean")

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

## Network Comparison Test ----

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

#recalculate networks with bootnet; output is input for the NCT
network_reading <- estimateNetwork(reading, default = "EBICglasso", corMethod = "spearman")
network_music   <- estimateNetwork(music  , default = "EBICglasso", corMethod = "spearman")

set.seed(1234)
NCT_reading_music <- NCT(network_reading, network_music, 
                         it               = 1000,
                         paired           = TRUE, # data is paired
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

## Scaling data ----
data <- mydata %>% select(t0_AF:t2_AE)
# standardize data across time points
data_scaled <- as.data.frame(scale(data))

# plot distributions of scaled variables
plot_hist_facet(data_scaled, bins = 8, ncol = 7)

# save file
saveRDS(data_scaled, file = "data_scaled.rds")
# read file
# data_scaled <- readRDS("data_scaled.rds")

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

# define design matrix (across 3 measurement points)
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

# it appears that the second and third time points (music) are more
# unidimensional than the first time point (reading)
par(mfrow = c(1, 1))
matplot(cbind(ev1,ev2,ev3),type = "l", ylab = "Eigenvalue")

# Labels to be used in graphs
labels <- rownames(design)

## Estimate saturated model----
model1 <- panelgvar(data      = data_scaled, 
                    vars      = design,
                    estimator = "ML",
                    storedata = TRUE)

## a) run the model ----
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

## Plot analytic confidence intervals ----
# The confidence intervals are only valid for the saturated model
tiff(filename = "Final_CIplots_within_saturated_model.tiff", width = 4000, height = 4000, res = 450)
par(mfrow = c(1,1))
CIplot(model1, "omega_zeta_within")
dev.off()

tiff(filename = "Final_CIplots_beta_saturated_model.tiff",width = 4000, height = 4000, res = 450)
par(mfrow = c(1,1))
CIplot(model1, "beta")
dev.off()

tiff("Final_CIplots_between_saturated_model.tiff",width = 4000, height = 4000, res = 450)
par(mfrow = c(1,1))
CIplot(model1, "omega_zeta_between")
dev.off()

## b) Prune to find a sparse model ----
model2 <- model1 %>% prune(alpha = 0.05, recursive = FALSE)
fit(model2) 
            
## c) stepup model ----
model3 <- model2 %>% stepup(criterion = "bic",alpha = 0.05)
fit(model3) 
            
model3 %>% parameters()

## e) compare all models ----
compare(baseline     = model1, 
        pruned       = model2,
        stepup       = model3) # <- best AIC and BIC

# 4 Plotting output ----
# Note that the PDC and the beta matrix should align but are expressed in 
# opposite directions
temporal        <- getmatrix(model3,"PDC") # equals t(getmatrix(model2, "beta"))
contemporaneous <- getmatrix(model3,"omega_zeta_within")
between         <- getmatrix(model3,"omega_zeta_between")

graph1 <- qgraph(temporal,layout = "spring")
graph2 <- qgraph(contemporaneous,layout = "spring")
graph3 <- qgraph(between,layout = "spring")

# Create average layout to be used consistently across all networks
L <- averageLayout(graph1,graph2, graph3)

names <- c("Focused Attention","Self Awareness","Visual Imagery",
           "Valence","Calmness", "One-topic centered thoughts","Dissociation")

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
             nodeNames = names, 
             groups = gr,
             vsize=12,
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
centralityPlot(cg,scale = "z-score", include = c("Strength"))
dev.off()

tiff(filename="strength centrality temporal network.tiff", width=1250, height=2500, res=450)
centralityPlot(tg,scale = "z-score", include = c("InStrength", "OutStrength"))
dev.off()

# 6 Bootstrap analysis ----
# keep number low (ca. 100), otherwise may take quite a long time
set.seed(1234)
nBoot <- 2

Bootstraps <- lapply(1:nBoot, function(x) {
  
  message("Simulation: ",x)
  
  # Sample from the data
  bootData <- mydata_t012[sample(1:nrow(mydata_t012), 
                                 round(0.75 * nrow(mydata_t012))), ]
  
  # Form model
  bootstrapped_model <-
    panelgvar(
      bootData,
      vars      = design,
      estimator = "ML")  
  
  # Run first time
  bootstrapped_model <- bootstrapped_model %>% runmodel %>% 
    
    # Fix the problematic parameters to zero:
    # diaglep <- diag(getmatrix(bootstrapped_model, "lowertri_epsilon_between"))
    # if (any(abs(diaglep) < 1e-6)){
    #   for (i in which(abs(diaglep) < 1e-6)){
    #     bootstrapped_model <- bootstrapped_model %>% 
    #       fixpar("lowertri_epsilon_between",i,i,value = 0)
    #   }
    # }
    
    
  prune(alpha = 0.05, recursive = FALSE) %>% 
    stepup(criterion = "bic")
  
  return(bootstrapped_model)
})

# 7 Stability Analysis ----
# 
# Extracting stability analyses to check the degree to which individual edges
# were included across all bootstrap samples 

##  bootstrapped results: temporal ----
# check if individual edges are included (y/n)
resBoots_temp <-
  lapply(Bootstraps, function(x)
    ifelse(getmatrix(x, "PDC") > 0, 1, 0))

# how often (%) is edge included, over all bootstrap iterations
apply(simplify2array(resBoots_temp), 1:2, mean) 
Bootstraps_temporal_est_df            <- as.data.frame(apply(simplify2array(resBoots_temp), 1:2, sum))
row.names(Bootstraps_temporal_est_df) <- colnames(Bootstraps_temporal_est_df) <- colnames(mydata_t0)

write.csv(Bootstraps_temporal_est_df, "Bootstraps_temporal_est_df.csv")

## bootstrapped results: contemporaneous ----
resBoots_cont <-
  lapply(Bootstraps, function(x)
    ifelse(getmatrix(x, "omega_zeta_within") > 0, 1, 0))

# how often (%) is edge included, over all bootstrap iterations
apply(simplify2array(resBoots_cont), 1:2, mean) 
Bootstraps_contemporaneous_est_df            <- as.data.frame(apply(simplify2array(resBoots_cont), 1:2, sum))
row.names(Bootstraps_contemporaneous_est_df) <- colnames(Bootstraps_contemporaneous_est_df) <- colnames(mydata_t0)

write.csv(Bootstraps_contemporaneous_est_df, "Bootstraps_contemporaneous_est_df.csv")

sessioninfo::session_info()
