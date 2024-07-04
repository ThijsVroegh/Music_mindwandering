#    A temporal network approach to music-induced mind wandering   
#                      L.Taruffi & T.Vroegh 2024
#
# data cleaning
# input:  raw spss data (MW.sav)
# output: data_cleaned.rds 

# 1. Load libraries ----
library(tidyverse)
library(mice)   

# remove all data in global environment
rm(list = ls())

# random permutation seed
set.seed(1234) 

# 2. Reading in data ----
mydata <- haven::read_sav("../data/MW.sav")

# spss codebook
value_labels <- labelled::val_labels(mydata)

# Exclude participants (Prolific IDs)
remove_prolific <- c("5f9aa2ab79a9a81cf6f8937e","5e28456cb66c4d016fd05e57","5d7e3844fc5b4b0001bb0e1b","61671675737af5033ff16c47","61589d1197ce8247fff72c03","5a4d6973f4c9dc0001649fc7",
                     "5dcc845678cda09214282fb9","6111a5dc375335bf33003edd","587e7be80cfbeb0001a62ae4","5fce4b1b20baeb026e5da5d3","614907aa62e00009ab079c33")

# 3. Data preprocessing----
mydata_tbl <- mydata %>% 
  
  # filter data on: non-preview cases, agree to participate,
  # exclude the prolific cases emailed by Liila, 
  # must have selected a musical genre
  filter(Status == 0 & Q1 == 1,
         !Q150 %in% remove_prolific,
         !is.na(Q3)) %>% #280
  
  select(-StartDate:-Q1) %>% 
  mutate(ID = 1:n()) %>% 
  select(ID,everything()) %>%
  select(-(starts_with("time"))) %>% 
  
  #combine data of different musical genres 
  mutate(Q4_1 = select(., starts_with("Q4a_1")) %>% rowMeans( na.rm = TRUE),
         Q4_2 = select(., starts_with("Q4a_2")) %>% rowMeans( na.rm = TRUE),
         Q4_3 = select(., starts_with("Q4a_3")) %>% rowMeans( na.rm = TRUE),
         Q4_4 = select(., starts_with("Q4a_4")) %>% rowMeans( na.rm = TRUE),
         Q4_5 = select(., starts_with("Q4a_5")) %>% rowMeans( na.rm = TRUE),
         Q4_6 = select(., starts_with("Q4a_6")) %>% rowMeans( na.rm = TRUE),
         Q4_7 = select(., starts_with("Q4a_7")) %>% rowMeans( na.rm = TRUE),
         Q4_8 = select(., starts_with("Q4a_8")) %>% rowMeans( na.rm = TRUE),
         Q4_9 = select(., starts_with("Q4a_9")) %>% rowMeans( na.rm = TRUE),
         Q4_10 = select(., starts_with("Q4b_1")) %>% rowMeans( na.rm = TRUE),
         Q4_11 = select(., starts_with("Q4b_2")) %>% rowMeans( na.rm = TRUE),
         Q4_12 = select(., starts_with("Q4b_3")) %>% rowMeans( na.rm = TRUE),
         Q4_13 = select(., starts_with("Q4b_4")) %>% rowMeans( na.rm = TRUE),
         Q4_14 = select(., starts_with("Q4b_5")) %>% rowMeans( na.rm = TRUE),
         Q4_15 = select(., starts_with("Q4b_6")) %>% rowMeans( na.rm = TRUE),
         Q4_16 = select(., starts_with("Q4b_7")) %>% rowMeans( na.rm = TRUE),
         Q4_17 = select(., starts_with("Q4b_8")) %>% rowMeans( na.rm = TRUE),
         Q4_18 = select(., starts_with("Q4b_9")) %>% rowMeans( na.rm = TRUE)) %>% 
  
  rename_at(.vars = vars(Q2a_1:Q2b_9), 
            .funs = ~ str_c("t0_", 1:18)) %>% 
  rename_at(.vars = vars(Q4_1:Q4_18), 
            .funs = ~ str_c("t1_", 1:18)) %>% 
  rename_at(.vars = vars(Q5a_1:Q5b_9), 
            .funs = ~ str_c("t2_", 1:18)) %>%
  
  # replace the value '8' (not applicable) with 'NA' to avoid misunderstanding its meaning
  mutate(across(c(t0_1:t0_18), ~na_if(., 8))) %>%
  mutate(across(c(t1_1:t1_18), ~na_if(., 8))) %>% 
  mutate(across(c(t2_1:t2_18), ~na_if(., 8))) %>% 
  
  rename(M_genre = Q3) %>%
  
  sjlabelled::remove_all_labels() %>% 
  
  rename(Amount_of_thoughts     = Q6_1,
         
         thoughts_valence        = Q7_1,
         thoughts_duration       = Q7_2,
         thoughts_meaning        = Q7_3,
         thoughts_musicunrelated = Q7_4,
         thoughts_control        = Q7_5,
         thoughts_amountoftopics = Q7_6,
         thoughts_not_intrusive  = Q7_7, 
         thoughts_memory         = Q7_8,
         thoughts_chaotic        = Q7_9,
         
         thoughts_keyword1  = Q8_1,
         thoughts_keyword2  = Q8_2,
         thoughts_keyword3  = Q8_3,
         
         #Q9 thoughts and time aspect
         #Q10 thoughts and subject
         #Q11 function of thoughts
         #Q12 music-related thoughts 
         
         Age         = Q16,
         Sex         = Q17,
         Gender      = Q18,
         Country     = Q19,
         Musician    = Q20) %>% 
  
  # Remove unnecessary columns
  select(-Q24,-Q150,-Q4a_1:-Q4b_9.M) %>% 
  
  # Questions on thought categories -> feature engineering
  #Na's to 0 -> 1= yes, 0 = no
  mutate(across(c(Q9_1:Q12_9), ~replace_na(., 0))) %>% 
  
  #Na's to 0 -> 1= yes, 0 = no
  mutate(across(c(Q21_1:Q21_5), ~replace_na(., 0))) %>% 
  
  rename(Anxiety        = Q21_1,
         ADD_ADHD       = Q21_2,
         Compulsive     = Q21_3,
         Autism         = Q21_4,
         Depression     = Q21_5) %>% 
  
  # recode 2 to 0 -> 1= yes, 0 = no
  mutate(Q22_1 = if_else(Q22_1 == 2, 0, Q22_1),
         Q22_2 = if_else(Q22_2 == 2, 0, Q22_2),
         Q22_3 = if_else(Q22_3 == 2, 0, Q22_3),
         Q22_4 = if_else(Q22_4 == 2, 0, Q22_4),
         Q22_5 = if_else(Q22_5 == 2, 0, Q22_5),
         Q22_6 = if_else(Q22_6 == 2, 0, Q22_6)) %>% 
  
  rename(Headphones     = Q22_1,
         Speakers       = Q22_2,
         Eyesopen       = Q22_3,
         Interrupted    = Q22_4,
         Notfastforward = Q22_5,
         Quietplace     = Q22_6) %>% 
  
  # Coding error in Amount_of_thoughts 
  # recode 12 to 6
  mutate(Amount_of_thoughts = if_else(Amount_of_thoughts == 12,
                                      6, Amount_of_thoughts)) %>% 
  
  # properly reorder all variables
  select(ID:t0_18,M_genre,t1_1:t1_18,everything())

# count rowwise missing values per measurement (T1,T2,T3,Traits)
count_na_func <- function(x) sum(is.na(x)) 
sum_func      <- function(x) sum((x))

mydata_tbl2 <- mydata_tbl %>%
  mutate(count_na_Q2    = apply(.[2:19]   , 1, count_na_func), #col.135
         count_na_Q4    = apply(.[21:38]  , 1, count_na_func), #col.136
         count_na_Q5    = apply(.[39:56]  , 1, count_na_func), #col.137
         
         count_na_IRI   = apply(.[104:110], 1, count_na_func),
         count_na_spont = apply(.[111:114], 1, count_na_func),
         count_na_delib = apply(.[115:118], 1, count_na_func)) %>% 
  
  # sum of missing values in three waves of measurements       
  mutate(sum_missing    = apply(.[135:137], 1, sum_func)) %>%
  
  ungroup() %>% 
  
  # delete id's with too many missing values
  filter(sum_missing <= 5) 

## Missing values ----

# check the data for percentage of missing values, column and rowwise
pMiss <- function(x){round(sum(is.na(x))/length(x)*100,2)}

apply(mydata_tbl2,2,pMiss) #columns
mydata_tbl2$perc_miss <- apply(mydata_tbl2,1,pMiss) #rows

# note: 8 subjects reported not to have had a single thought; hence,11 missings

# missing completely at random (MCAR)? 

#This is the desirable scenario in case of missing data
# If the p value for Little's MCAR test is not significant, 
# then the data may be assumed to be MCAR
mydata_tbl2 %>% select(t0_1:t0_18, t1_1:t1_18,t2_1:t2_18) %>% 
  naniar::mcar_test()

# imputation of missing values for t0_1 to t2_18
init  <- mice(mydata_tbl2, maxit = 0)
meth  <- init$method
predM <- init$predictorMatrix

# exclude variables being used in imputation process 
predM[,c("ID",
         "thoughts_keyword1",
         "thoughts_keyword2",
         "thoughts_keyword3",
         "Country")] = 0

subset <- mydata_tbl2 %>% 
  select(t0_1:t0_18, t1_1:t1_18,t2_1:t2_18) %>% 
  names()

meth[!names(mydata_tbl2) %in% subset] <- ""

set.seed(1234)
imputed <- mice(mydata_tbl2, 
                method = meth, 
                predictorMatrix = predM,
                m = 5)

df_imputed <- complete(imputed)
sapply(df_imputed, function(x) sum(is.na(x)))

df_imputed <- df_imputed %>% 
  mutate(t1_8  = replace_na(t1_8 ,mean(t1_8, na.rm = TRUE)))
apply(df_imputed,2,pMiss)

## Summed constructs ----
mydata_tbl3 <- df_imputed %>% 
  
  # reverse items 3 and 4 of IRI fantasy (reverse-coded)
  mutate(Q13_IRI_Fantasy_3 = 6 - Q13_IRI_Fantasy_3,
         Q13_IRI_Fantasy_4 = 6 - Q13_IRI_Fantasy_4) %>% 
  
  rowwise() %>% 
  mutate(IRI_Fantasy    = sum(c_across(starts_with("Q13_IRI_Fantasy")), na.rm = TRUE)/7,
         Deliberate_MW  = sum(c_across(Q14_Del_Spont_MW_1:Q14_Del_Spont_MW_4), na.rm = TRUE)/4,
         Spontaneous_MW = sum(c_across(Q14_Del_Spont_MW_5:Q14_Del_Spont_MW_8), na.rm = TRUE)/4
  ) %>% 
  
  ungroup() %>% # undo rowwise()
  
  # Creating the nine dimensions T0 (1st measurement point)
  mutate(
    t0_AF   = (t0_11    + (8-t0_10))/2, # absorption
    t0_DIR  = (t0_1)  /1,               # inner directed attention
    t0_VC   = (8-t0_2)/1,               # self-control
    t0_SA   = (t0_3     + t0_12)/2,     # self-awareness
    t0_IM   = (t0_4     + t0_13)/2,     # visual imagery
    t0_ID   = (t0_5     + t0_14)/2,     # internal dialogue
    t0_VA   = (t0_6     + (8-t0_15))/2, # positive valence
    t0_CALM = (t0_7     + t0_16)/2,     # relaxation
    t0_THO  = (t0_8     + (8-t0_17))/2, # one-topic thoughts
    t0_AE   = (t0_9     + t0_18)/2      # altered awareness
  ) %>% 
  
  # Creating the nine dimensions T1 (2nd measurement point - halfway music)
  mutate(
    t1_AF   = (t1_11     + (8-t1_10))/2, 
    t1_DIR  = (t1_1) /1,     
    t1_VC   = (8-t1_2) /1,     
    t1_SA   = (t1_3     + t1_12)/2,     
    t1_IM   = (t1_4     + t1_13)/2,     
    t1_ID   = (t1_5     + t1_14)/2,     
    t1_VA   = (t1_6     + (8-t1_15))/2, 
    t1_CALM = (t1_7     + t1_16)/2,     
    t1_THO  = (t1_8     + (8-t1_17))/2, 
    t1_AE   = (t1_9     + t1_18)/2      
  ) %>% 
  
  # Creating the nine dimensions T2 (3rd measurement point)
  mutate(
    t2_AF   = (t2_11    + (8-t1_10))/2,
    t2_DIR  = (t2_1) /1,     
    t2_VC   = (8-t2_2)/1,     
    t2_SA   = (t2_3     + t2_12)/2,     
    t2_IM   = (t2_4     + t2_13)/2,     
    t2_ID   = (t2_5     + t2_14)/2,     
    t2_VA   = (t2_6     + (8-t2_15))/2, 
    t2_CALM = (t2_7     + t2_16)/2,     
    t2_THO  = (t2_8     + (8-t2_17))/2, 
    t2_AE   = (t2_9     + t2_18)/2
  ) 

## alphas and inter-item correlations ----

# alphas of measures
IRI_alpha <- alpha(mydata_tbl3[c('Q13_IRI_Fantasy_1',
                                 'Q13_IRI_Fantasy_2',
                                 'Q13_IRI_Fantasy_3',
                                 'Q13_IRI_Fantasy_4',
                                 'Q13_IRI_Fantasy_5',
                                 'Q13_IRI_Fantasy_6',
                                 'Q13_IRI_Fantasy_7')]) #.80

Deliberate_MW_alpha <- alpha(mydata_tbl3[c('Q14_Del_Spont_MW_1',
                                           'Q14_Del_Spont_MW_2',
                                           'Q14_Del_Spont_MW_3',
                                           'Q14_Del_Spont_MW_4')]) #.86

Spontaneous_MW <- alpha(mydata_tbl3[c('Q14_Del_Spont_MW_5',
                                      'Q14_Del_Spont_MW_6',
                                      'Q14_Del_Spont_MW_7',
                                      'Q14_Del_Spont_MW_8')])  #.81 

## inter-item correlations
correlations <- list(
  AF_alpha = c(),
  SA_alpha = c(),
  IM_alpha = c(),
  ID_alpha = c(),
  THO_alpha = c(),
  VAL_alpha = c(),
  CALM_alpha = c(),
  AE_alpha = c()
)

# Calculate correlations for each time point and store them in the list
time_points <- c(0, 1, 2)
for (t in time_points) {
  correlations$AF_alpha   <- c(correlations$AF_alpha, cor(mydata_tbl3[[paste0('t', t, '_11')]], 8 - mydata_tbl3[[paste0('t', t, '_10')]]))
  correlations$SA_alpha   <- c(correlations$SA_alpha, cor(mydata_tbl3[[paste0('t', t, '_3')]], mydata_tbl3[[paste0('t', t, '_12')]]))
  correlations$IM_alpha   <- c(correlations$IM_alpha, cor(mydata_tbl3[[paste0('t', t, '_4')]], mydata_tbl3[[paste0('t', t, '_13')]]))
  correlations$ID_alpha   <- c(correlations$ID_alpha, cor(mydata_tbl3[[paste0('t', t, '_5')]], mydata_tbl3[[paste0('t', t, '_14')]]))
  correlations$THO_alpha  <- c(correlations$THO_alpha, cor(mydata_tbl3[[paste0('t', t, '_8')]], 8 - mydata_tbl3[[paste0('t', t, '_17')]]))
  correlations$VAL_alpha  <- c(correlations$VAL_alpha, cor(mydata_tbl3[[paste0('t', t, '_6')]], 8 - mydata_tbl3[[paste0('t', t, '_15')]]))
  correlations$CALM_alpha <- c(correlations$CALM_alpha, cor(mydata_tbl3[[paste0('t', t, '_7')]], mydata_tbl3[[paste0('t', t, '_16')]]))
  correlations$AE_alpha   <- c(correlations$AE_alpha, cor(mydata_tbl3[[paste0('t', t, '_9')]], mydata_tbl3[[paste0('t', t, '_18')]]))
}

# Calculate the mean correlation for each variable
mean_correlations <- sapply(correlations, mean)

mydata_tbl4 <- mydata_tbl3 %>% 
  
  # reorder all variables in proper order
  select(ID, M_genre, 
         t0_AF:t2_AE, 
         Amount_of_thoughts: Q11_10,
         IRI_Fantasy: Spontaneous_MW,
         Age,Sex,Musician) %>% 
  
  # we don't include these variables in the network analysis
  select(-t0_DIR, -t1_DIR, -t2_DIR,
         -t0_VC, -t1_VC, -t2_VC,
         -t0_ID, -t1_ID, -t2_ID,
         -thoughts_keyword1, -thoughts_keyword2, -thoughts_keyword3,
         -M_genre) 

# save data for further exploratory data analysis and network modeling
saveRDS(mydata_tbl4,"data_cleaned.rds")
