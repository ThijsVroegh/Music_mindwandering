#    A temporal network approach to music-induced mind wandering   
#                      L.Taruffi & T.Vroegh 2025

# Exploratory data analysis ----

library(tidyverse)
library(reshape2)
library(corrplot)     
library(psy)          
library(networktools) 
library(psych)

# remove all data in global environment
rm(list = ls())

# source functions
source("sources/plot_hist_facet.R")

# read in cleaned data
mydata <- readRDS("rds/data_cleaned.rds")

## a) demographics Figure S0 ----

# sex
table(mydata$Sex) # 170 males, 180 females, 2 prefer not to say

# age
table(mydata$Age)

# Given data
categories <- c("Under 18", "18 - 24", "25 - 34", "35 - 44", "45 - 54", "55 - 64",
                "65 - 74", "75 - 84", "85 or older")
codes      <- 1:9
counts     <- c(0, 61, 105, 78, 52, 36, 17, 2, 1)
midpoints  <- c(8, 21, 29.5, 39.5, 49.5, 59.5, 69.5, 79.5, 90)

# Calculate the weighted mean
weighted_mean <- sum(midpoints * counts) / sum(counts)

# Calculate the weighted standard deviation
weighted_sd <- sqrt(sum(counts * (midpoints - weighted_mean)^2) / sum(counts - 1))

# Print the results
cat("Weighted Mean:", weighted_mean, "\n")
cat("Weighted Standard Deviation:", weighted_sd, "\n")

table(mydata$Musician)

# country
table(mydata$Country)

# summary table for demographics (counts and percentages)
age_labels <- c(
  "Under 18",
  "18 - 24",
  "25 - 34",
  "35 - 44",
  "45 - 54",
  "55 - 64",
  "65 - 74",
  "75 - 84",
  "85 or older",
  "Prefer not to say"
)

sex_labels <- c("Male", "Female", "Prefer not to say")

gender_labels <- c(
  "Male",
  "Female",
  "Non-binary / third gender",
  "Prefer not to say"
)

musician_labels <- c(
  "Non-musician",
  "Amateur musician",
  "Semi-professional musician",
  "Professional musician",
  "Prefer not to say"
)

# Country labels - only for countries with respondents
country_labels <- c(
  "8" = "Armenia",
  "9" = "Australia",
  "10" = "Austria",
  "14" = "Bangladesh",
  "17" = "Belgium",
  "24" = "Brazil",
  "26" = "Bulgaria",
  "31" = "Canada",
  "36" = "China",
  "42" = "Croatia",
  "44" = "Cyprus",
  "45" = "Czech Republic",
  "60" = "Finland",
  "61" = "France",
  "65" = "Germany",
  "67" = "Greece",
  "75" = "Hong Kong (S.A.R.)",
  "76" = "Hungary",
  "78" = "India",
  "84" = "Italy",
  "85" = "Jamaica",
  "87" = "Jordan",
  "104" = "Malaysia",
  "122" = "Netherlands",
  "123" = "New Zealand",
  "126" = "Nigeria",
  "128" = "Norway",
  "130" = "Pakistan",
  "136" = "Philippines",
  "137" = "Poland",
  "138" = "Portugal",
  "142" = "Romania",
  "143" = "Russian Federation",
  "156" = "Singapore",
  "158" = "Slovenia",
  "161" = "South Africa",
  "163" = "Spain",
  "169" = "Switzerland",
  "179" = "Turkey",
  "184" = "United Arab Emirates",
  "185" = "United Kingdom",
  "187" = "United States",
  "1357" = "Zimbabwe",
  "1358" = "Prefer not to say"
)

make_demo_table <- function(x, labels, var_name) {
  tibble::tibble(
    variable = var_name,
    category = factor(x, levels = seq_along(labels), labels = labels)
  ) %>%
    dplyr::count(variable, category, name = "n") %>%
    dplyr::mutate(percent = round(100 * n / sum(n), 1))
}

# Special function for Country (non-sequential codes)
make_country_table <- function(x, labels, var_name) {
  tibble::tibble(
    variable = var_name,
    category = labels[as.character(x)]
  ) %>%
    dplyr::count(variable, category, name = "n") %>%
    dplyr::mutate(percent = round(100 * n / sum(n), 1)) %>%
    dplyr::arrange(dplyr::desc(n))
}

demo_tables <- list(
  make_demo_table(mydata$Age, age_labels, "Age"),
  make_demo_table(mydata$Sex, sex_labels, "Sex"),
  make_demo_table(mydata$Musician, musician_labels, "Self-reported musician"),
  make_country_table(mydata$Country, country_labels, "Country")
)

if ("Gender" %in% names(mydata)) {
  demo_tables <- c(
    demo_tables[1:3],
    list(make_demo_table(mydata$Gender, gender_labels, "Gender")),
    demo_tables[4]
  )
}

demo_summary <- dplyr::bind_rows(demo_tables)

write.csv(demo_summary, "Demographics_summary.csv", row.names = FALSE)

# publication-ready demographics table (Word)
if (requireNamespace("flextable", quietly = TRUE) && requireNamespace("officer", quietly = TRUE)) {
  demo_summary_pub <- demo_summary %>%
    dplyr::mutate(n_percent = sprintf("%d (%.1f%%)", n, percent)) %>%
    dplyr::select(variable, category, n_percent)

  demo_ft <- flextable::flextable(demo_summary_pub)
  demo_ft <- flextable::set_header_labels(
    demo_ft,
    variable = "Variable",
    category = "Category",
    n_percent = "n (%)"
  )
  demo_ft <- flextable::autofit(demo_ft)
  demo_ft <- flextable::fontsize(demo_ft, size = 10, part = "all")
  demo_ft <- flextable::bold(demo_ft, part = "header")

  doc <- officer::read_docx()
  doc <- officer::body_add_par(doc, "S0. Participant demographics", style = "heading 1")
  doc <- flextable::body_add_flextable(doc, demo_ft)
  print(doc, target = "S0_demographics.docx")
} else {
  message("Packages 'flextable' and 'officer' are required for Word export.")
}


## b) contents of thought (time, subject, function) Figure 1 a b c ----

# 1) time orientation (Q9)
mydata_time_thought <- mydata %>% dplyr::select(Q9_1:Q9_7) %>% as.data.frame()

# Summarize the data to get the frequency of each option
summary_data <- mydata_time_thought %>% summarise_all(sum) 

# Convert the data to long format for plotting
long_data <- summary_data %>% gather(key = "tickbox_option", value = "count")

# Define labels for each tickbox option
tickbox_labels <- c("More than 1 year ago",
                    "Only recently",
                    "Present time",
                    "Near future",
                    "More than 1 year in the future",
                    "My thoughts had no clear time orientation",
                    "Not sure")

# Wrap the labels to the next line if they are too long
wrapped_labels <- str_wrap(tickbox_labels, width = 10)  # Adjust the width parameter as needed

## Export
tiff("Fig 1a Distribution of time orientation of reported thoughts.tiff", width = 2200, height = 1000, units = "px", res = 300)
ggplot(long_data, aes(x = tickbox_option, y = count, fill = tickbox_option)) +
  geom_bar(stat = "identity") +
  viridis::scale_fill_viridis(discrete = TRUE) +
  scale_x_discrete(labels = wrapped_labels) +  # Add wrapped labels to x-axis
  labs(title = "Distribution of time orientation of reported thoughts",
       subtitle = expression(italic("My thoughts involved events related to...")),
       y = "Count",
       x = "",
       fill = "Tickbox Option") +  # Label for legend
  theme_minimal() +
  theme(legend.position = "none",axis.text.x = element_text(size = 8))  # Remove legend
dev.off()

# 2) subject of thought (Q10) 
mydata_subject_thought <- mydata %>% dplyr::select(Q10_1:Q10_8) %>% as.data.frame()

# Summarize the data to get the frequency of each option
summary_data <- mydata_subject_thought %>% summarise_all(sum) 

# Convert the data to long format for plotting
long_data <- summary_data %>% gather(key = "tickbox_option", value = "count")

# Define labels for each tickbox option
tickbox_labels <- c("myself",
                    "other people",
                    "fictional situations and/or characters",
                    "work/study",
                    "the music",
                    "this study",
                    "other",
                    "Not sure")

# Wrap the labels to the next line if they are too long
wrapped_labels <- str_wrap(tickbox_labels, width = 10)  # Adjust the width parameter as needed

tiff("Fig 1b Distribution of topic of reported thoughts.tiff", width = 2200, height = 1000, units = "px", res = 300)
ggplot(long_data, aes(x = tickbox_option, y = count, fill = tickbox_option)) +
  geom_bar(stat = "identity") +
  viridis::scale_fill_viridis(discrete = TRUE) +
  scale_x_discrete(labels = wrapped_labels) +  # Add wrapped labels to x-axis
  labs(title = "Distribution of topic of reported thoughts",
       subtitle = expression(italic("My thoughts were about...")),
       y = "Count",
       x = "",
       fill = "Tickbox Option") + 
  theme_minimal() +
  theme(legend.position = "none",axis.text.x = element_text(size = 8))
dev.off()

# 3) function of thought (Q11) 
mydata_function_thought <- mydata %>% dplyr::select(Q11_1:Q11_10) %>% as.data.frame()

# Summarize the data to get the frequency of each option
summary_data <- mydata_function_thought %>% summarise_all(sum) 

long_data <- summary_data %>% gather(key = "tickbox_option", value = "count")

# Define labels for each tickbox option
tickbox_labels <- c("taking a decision",
                    "planning something",
                    "re-evaluating something",
                    "feeling better about myself",
                    "daydreaming or fantasizing",
                    "dealing with something personal",
                    "entertaining myself",
                    "solutions to problems",
                    "other",
                    "not sure")

# Wrap the labels to the next line if they are too long
wrapped_labels <- str_wrap(tickbox_labels, width = 10)  # Adjust the width parameter as needed

## Export
tiff("Fig 1c Distribution of functions of thought.tiff", width = 2200, height = 1000, units = "px", res = 300)
ggplot(long_data, aes(x = tickbox_option, y = count, fill = tickbox_option)) +
  geom_bar(stat = "identity") +
  viridis::scale_fill_viridis(discrete = TRUE) +
  scale_x_discrete(labels = wrapped_labels) + 
  labs(title = "Distribution of functions of thought",
       subtitle = expression(italic("My thoughts had to do with...")),
       y = "Count",
       x = "",
       fill = "Tickbox Option") +  
  theme_minimal() +
  theme(legend.position = "none",axis.text.x = element_text(size = 8))
dev.off()

## c) structure of thoughts correlationsFigure 2 ----
mydata_thoughts <- mydata %>% 
  dplyr::select(Amount_of_thoughts:thoughts_chaotic) %>% 
  as.data.frame() %>% 
  rename(number           = Amount_of_thoughts,
         valence          = thoughts_valence,
         duration         = thoughts_duration,
         importance       = thoughts_meaning,
         task_relatedness = thoughts_musicunrelated,
         intentionality   = thoughts_control,
         topics           = thoughts_amountoftopics,
         intrusiveness    = thoughts_not_intrusive,
         memory           = thoughts_memory,
         structure        = thoughts_chaotic)

#remove rows with missing values in any column
mydata_thoughts <- mydata_thoughts[complete.cases(mydata_thoughts), ]
cor_thoughts <- round(cor(mydata_thoughts, method = "spearman"),2)

# Set up the color mapping
col <- colorRampPalette(c("white", "red"))(100)

tiff(filename = "Fig 2 Correlations thoughts.tiff", width = 4000, height = 4000, res = 450)
par(mfrow = c(1,1))
corrplot(cor_thoughts, 
         method = "color", 
         addCoef.col = ifelse(cor_thoughts > 0.36 & row(cor_thoughts) != col(cor_thoughts), "red", "black"),
         number.cex = ifelse(cor_thoughts > 0.36 & row(cor_thoughts) != col(cor_thoughts), 0.9, 0.8),
         number.font = 1,
         cl.pos = "n",
         tl.col = "black",
         diag = TRUE, 
         insig = "blank")
dev.off()

## distributions of thought variables
tiff(filename = "Histograms of structure variables on thoughts.tiff", width = 6000, height = 2250, res = 450)
par(mfrow = c(1,1))
plot_hist_facet(mydata_thoughts, bins = 7)
dev.off()

tiff(filename = "Densities of thought items.tiff", width = 4500, height = 4500, res = 450)
par(mfrow = c(4, 3), mar = c(4, 2, 2, 2))
for (var in colnames(mydata_thoughts)) {
  x = mydata_thoughts[, var]
  hist(x, prob = TRUE, xlab = var, ylab = "", main = "", col = "ivory")
  lines(density(x), lwd = 2, col = "tomato")
  curve(dnorm(x, mean = mean(x), sd = sd(x)), from = min(x), to = max(x),
        add = TRUE, lwd = 2, col = "steelblue")
}
dev.off()


## d) means of experiential dimensions Figure 3 ----
temp <-  mydata %>% dplyr::select(t0_AF:t0_AE)
colnames(temp) <- c("AF", "SA","IM","VA","CALM","THO","AE")

# from wide to long format
datalong <- mydata %>% dplyr::select(t0_AF:t2_AE) %>% 
  pivot_longer(starts_with('t'), 
               names_to = c('time', '.value'), 
               names_sep = '\\_') 

freq <- as.data.frame(describeBy(datalong[,c(2:8)], 
                                 group = datalong$time, mat = T))
freq$err <- qnorm(.975) * freq$se
freq <- as.data.frame(freq[,c("group1","vars","mean", "n","se", "err")])

freq$vars <- factor(freq$vars, levels = c(1:7), labels = colnames(temp[,c(1:7)]))
freq$vars
rownames(freq) <- c()
freq

# Plot
freqbar <- ggplot(freq, aes(fill = group1, y = mean, x = vars)) + 
  geom_bar(position="dodge", stat="identity", color="black", size=0.3) +
  ylim(0,0.45) + labs(x="", y="") +
  geom_errorbar(aes(ymin=mean - err, ymax=mean + err), 
                width=.2, size=0.3, position=position_dodge(.9)) +
  theme(text = element_text(size=13), legend.position = c(0.88,0.83)) +
  theme_bw() +
  scale_fill_manual(values=viridis::viridis_pal()(5)[3:5])

## Export
tiff("Fig 3 Figure frequencies_ci_7dim.tiff", width = 2200, height = 1000, units = "px", res = 300)
print(freqbar + scale_y_continuous(name = "mean"))
dev.off()

## e) correlations between dimensions Figure S2 ----
mydata_t0 <- mydata %>% dplyr::select(t0_AF:t0_AE) %>% as.data.frame() 
mydata_t1 <- mydata %>% dplyr::select(t1_AF:t1_AE) %>% as.data.frame()
mydata_t2 <- mydata %>% dplyr::select(t2_AF:t2_AE) %>% as.data.frame()

cor_t0 <- round(cor(mydata_t0, method = "spearman"),2)
cor_t1 <- round(cor(mydata_t1, method = "spearman"),2)
cor_t2 <- round(cor(mydata_t2, method = "spearman"),2)

#check pd  
corpcor::is.positive.definite(cor_t0)
corpcor::is.positive.definite(cor_t1)
corpcor::is.positive.definite(cor_t2)

# Rename the column names of the dataframes so that they are equal and can be 
# used in the NCT (i.e. removing the t0, t1, and t2 in column names)
mydata_t0 <- mydata_t0 %>% 
  rename_at(vars(starts_with("t0_")),  ~ str_replace(., "t0_", ""))
mydata_t1 <- mydata_t1 %>% 
  rename_at(vars(starts_with("t1_")),~ str_replace(., "t1_", ""))
mydata_t2 <- mydata_t2 %>% 
  rename_at(vars(starts_with("t2_")),~ str_replace(., "t2_", ""))

# Plotting
tiff(filename = "Fig S2 Figure corr_pearson.tiff", width = 6000, height = 2250, res = 450)
par(mfrow = c(1,3))

corrplot(cor_t0, method = "color", addCoef.col = "black",
         number.cex = .9, cl.pos = "n", diag = T, insig = "blank")
title(main = "Reading task")

corrplot(cor_t1, method = "color", addCoef.col = "black",
         number.cex = .9, cl.pos = "n", diag = T, insig = "blank")
title(main = "Music- first measurement")

corrplot(cor_t2, method = "color", addCoef.col = "black",
         number.cex = .9, cl.pos = "n", diag = T, insig = "blank")
title(main = "Music - second measurement")
dev.off()

## f) boxplots ----
mycols <- c("#1a2a6c", "#b21f1f", "#fdbb2d")

par(mfrow = c(1, 1))
boxplot(mydata_t0,varwidth = FALSE,  col=(c(values = colorRampPalette(mycols)(10))), notch=TRUE, ylab = "mean intensity", xlab="Dimension", main ="T0")
boxplot(mydata_t1, varwidth = FALSE, col=(c(values = colorRampPalette(mycols)(10))), notch=TRUE, ylab = "mean intensity", xlab="Dimension", main ="T1")
boxplot(mydata_t2, varwidth = FALSE, col=(c(values = colorRampPalette(mycols)(10))), notch=TRUE, ylab = "mean intensity", xlab="Dimension", main ="T2")

## g) distributions ----
tiff(filename="densities_t0.tiff", width=4500, height=2500, res=450)
plot_hist_facet(mydata_t0)
dev.off()

tiff(filename="densities_t1.tiff", width=4500, height=2500, res=450)
plot_hist_facet(mydata_t1)
dev.off()

tiff(filename="densities_t2.tiff", width=4500, height=2500, res=450)
plot_hist_facet(mydata_t2)
dev.off()

## h) testing for collinearity with goldbricker function ----
## Use goldbricker function to search for potential "bad pairs"
## This function compares correlations in a psychometric network in order to identify
## nodes which most likely measure the same underlying construct (i.e., are colinear).
## A "bad pairs" includes two nodes that are highly correlated with each other and that 
## correlate with other nodes in the network in a highly similar pattern.
gb_t0 <- goldbricker(mydata_t0, p = 0.05, method = "hittner2003", threshold = 0.25, corMin = .50)
gb_t1 <- goldbricker(mydata_t1, p = 0.05, method = "hittner2003", threshold = 0.25, corMin = .50)
gb_t2 <- goldbricker(mydata_t2, p = 0.05, method = "hittner2003", threshold = 0.25, corMin = .50)

## Produce list of "bad pairs"
gb_t0 #  Less than 25 % of correlations are significantly different for the following pairs: CALM & VA 
gb_t1 #  "No suggested reductions"
gb_t2 #  "No suggested reductions"

