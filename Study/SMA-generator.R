# this script is used to generate synthetic SMA data
set.seed(2023)

## loading data
load("~/rdborrow/Study/SMA_data_simulation.rda")

## target distribution
N1 = 100
N0 = 100
SMA_1_simulated = list()
SMA_0_simulated = list()

## fetch variables
A = SMA_Xover_data$A
S = SMA_Xover_data$S
Y0 = SMA_Xover_data$Y0
SMA_Type = SMA_Xover_data$SMA_Type
Scoliosis = SMA_Xover_data$Scoliosis
SMN2_Copy_Number = SMA_Xover_data$SMN2_Copy_Number
Age_Enrollment = SMA_Xover_data$Age_Enrollment

## remove missing values
SMA_Xover_data = na.omit(SMA_Xover_data)
N = nrow(SMA_Xover_data)

sAge = 0.5
sY0  = 1
sY1  = 1
sY2  = 1
sY3  = 1
sY4  = 1
snum1 = 1
snum0 = 1

## compute perturbed Age and Y0 using exponential mechanism


# SMA_Xover_data$Y0 = SMA_Xover_data$Y0 + rmutil::rlaplace(n = N, m = 0, s = sY0)
# SMA_Xover_data$Age_Enrollment = SMA_Xover_data$Age_Enrollment + rmutil::rlaplace(n = N, m = 0, s = sAge)
# SMA_Xover_data$Y1 = SMA_Xover_data$Y1 + rmutil::rlaplace(n = N, m = 0, s = sY1)
# SMA_Xover_data$Y2 = SMA_Xover_data$Y2 + rmutil::rlaplace(n = N, m = 0, s = sY2)
# SMA_Xover_data$Y3 = SMA_Xover_data$Y3 + rmutil::rlaplace(n = N, m = 0, s = sY3)
# SMA_Xover_data$Y4 = SMA_Xover_data$Y4 + rmutil::rlaplace(n = N, m = 0, s = sY4)


## Separate data according to trial status
SMA_0 = SMA_Xover_data %>% filter(S == 0) %>% 
  group_by(SMA_Type, SMN2_Copy_Number, Scoliosis) %>% 
  mutate(groupID = cur_group_id()) %>% ungroup()
SMA_1 = SMA_Xover_data %>% filter(S == 1)  %>% 
  group_by(SMA_Type, SMN2_Copy_Number, Scoliosis) %>% 
  mutate(groupID = cur_group_id()) %>% ungroup()
SMA_0_kde = list()
SMA_1_kde = list()


## build KDEs for S == 1
cat_levels = SMA_1 %>% select(SMA_Type, SMN2_Copy_Number, Scoliosis) %>%
  group_by(SMA_Type, SMN2_Copy_Number, Scoliosis) %>%
  summarize(num = n()) %>% ungroup()
cat_levels$perturb_num = sapply(cat_levels$num + rmutil::rlaplace(n = nrow(cat_levels), m = 0, s = snum1), function(x){max(x, 0)})
cat_levels$perturb_prob = cat_levels$perturb_num/sum(cat_levels$perturb_num)
SMA_1_simulated_num = rmultinom(1, N1, cat_levels$perturb_prob)
  
# for (i in 1:nrow(cat_levels)){
#   if (cat_levels$num[i] > 10){
#     SMA_1_simulated[[i]] = cbind(cat_levels %>%  
#                                    select(SMA_Type, SMN2_Copy_Number, Scoliosis) %>% slice(rep(i, SMA_1_simulated_num[i])), 
#                                  rkde(n = SMA_1_simulated_num[i], fhat = kde(SMA_1  %>% filter(groupID == i) %>% 
#                                                                                select(Age_Enrollment, Y0, Y1, Y2, Y3, Y4))))
#   }else if(cat_levels$num[i] > 0 && cat_levels$num[i] <= 10){
#     SMA_1_simulated[[i]] = cbind(cat_levels %>% 
#                                    select(SMA_Type, SMN2_Copy_Number, Scoliosis) %>% slice(rep(i, SMA_1_simulated_num[i])),
#                                  SMA_1 %>% filter(groupID == i) %>% sample_n(SMA_1_simulated_num[i], replace = TRUE) %>% 
#                                    select(Age_Enrollment, Y0, Y1, Y2, Y3, Y4))
#   }
# }

for (i in 1:nrow(cat_levels)){
  if(SMA_1_simulated_num[i] > 0){
    SMA_1_simulated[[i]] = cbind(cat_levels %>% 
                                   select(SMA_Type, SMN2_Copy_Number, Scoliosis) %>% slice(rep(i, SMA_1_simulated_num[i])),
                                 SMA_1 %>% filter(groupID == i) %>% sample_n(SMA_1_simulated_num[i], replace = TRUE) %>% 
                                   select(Age_Enrollment, Y0, Y1, Y2, Y3, Y4))
  }
    
}

SMA_1_simulated = bind_rows(SMA_1_simulated)
SMA_1_simulated = SMA_1_simulated %>% mutate(S = 1) %>% mutate(A = rbinom(N1, 1, prob = 2/3))
SMA_1_simulated$Age_Enrollment = SMA_1_simulated$Age_Enrollment + rmutil::rlaplace(n = N1, m = 0, s = sAge)
SMA_1_simulated$Y0 = SMA_1_simulated$Y0 + rmutil::rlaplace(n = N1, m = 0, s = sY0)
SMA_1_simulated$Y1 = SMA_1_simulated$Y1 + rmutil::rlaplace(n = N1, m = 0, s = sY1)
SMA_1_simulated$Y2 = SMA_1_simulated$Y2 + rmutil::rlaplace(n = N1, m = 0, s = sY2)
SMA_1_simulated$Y3 = SMA_1_simulated$Y3 + rmutil::rlaplace(n = N1, m = 0, s = sY3)
SMA_1_simulated$Y4 = SMA_1_simulated$Y4 + rmutil::rlaplace(n = N1, m = 0, s = sY4)
SMA_1_simulated$Age_Enrollment = ceiling(SMA_1_simulated$Age_Enrollment) + 1



## generalized additive mixed effects model







## build KDEs for S == 0
cat_levels = SMA_0 %>% select(SMA_Type, SMN2_Copy_Number, Scoliosis) %>%
  group_by(SMA_Type, SMN2_Copy_Number, Scoliosis) %>%
  summarize(num = n()) %>% ungroup()
cat_levels$perturb_num = sapply(cat_levels$num + rmutil::rlaplace(n = nrow(cat_levels), m = 0, s = snum0), function(x){max(x, 0)})
cat_levels$perturb_prob = cat_levels$perturb_num/sum(cat_levels$perturb_num)
SMA_0_simulated_num = rmultinom(1, N0, cat_levels$perturb_prob)

# for (i in 1:nrow(cat_levels)){
#   if (cat_levels$num[i] > 10){
#     SMA_0_simulated[[i]] = cbind(cat_levels %>%  
#                                    select(SMA_Type, SMN2_Copy_Number, Scoliosis) %>% slice(rep(i, SMA_0_simulated_num[i])), 
#                                  rkde(n = SMA_0_simulated_num[i], fhat = kde(SMA_0 %>% filter(groupID == i) %>% 
#                                                                                select(Age_Enrollment, Y0, Y1, Y2, Y3, Y4))))
#   }else if(cat_levels$num[i] > 0 && cat_levels$num[i] <= 10){
#     SMA_0_simulated[[i]] = cbind(cat_levels %>% 
#                                    select(SMA_Type, SMN2_Copy_Number, Scoliosis) %>% slice(rep(i, SMA_0_simulated_num[i])),
#                                  SMA_0 %>% filter(groupID == i) %>% sample_n(SMA_0_simulated_num[i], replace = TRUE) %>% 
#                                    select(Age_Enrollment, Y0, Y1, Y2, Y3, Y4))
#   }
# }

for (i in 1:nrow(cat_levels)){
  if(SMA_0_simulated_num[i]){
    SMA_0_simulated[[i]] = cbind(cat_levels %>% 
                                   select(SMA_Type, SMN2_Copy_Number, Scoliosis) %>% slice(rep(i, SMA_0_simulated_num[i])),
                                 SMA_0 %>% filter(groupID == i) %>% sample_n(SMA_0_simulated_num[i], replace = TRUE) %>% 
                                   select(Age_Enrollment, Y0, Y1, Y2, Y3, Y4))
  }
}

SMA_0_simulated = bind_rows(SMA_0_simulated)
SMA_0_simulated = SMA_0_simulated %>% mutate(S = 0) %>% mutate(A = 0)
SMA_0_simulated$Age_Enrollment = SMA_0_simulated$Age_Enrollment + rmutil::rlaplace(n = N0, m = 0, s = sAge)
SMA_0_simulated$Y0 = SMA_0_simulated$Y0 + rmutil::rlaplace(n = N0, m = 0, s = sY0)
SMA_0_simulated$Y1 = SMA_0_simulated$Y0 + rmutil::rlaplace(n = N0, m = 0, s = sY1)
SMA_0_simulated$Y2 = SMA_0_simulated$Y0 + rmutil::rlaplace(n = N0, m = 0, s = sY2)
SMA_0_simulated$Y3 = SMA_0_simulated$Y0 + rmutil::rlaplace(n = N0, m = 0, s = sY3)
SMA_0_simulated$Y4 = SMA_0_simulated$Y0 + rmutil::rlaplace(n = N0, m = 0, s = sY4)
SMA_0_simulated$Age_Enrollment = ceiling(SMA_0_simulated$Age_Enrollment)




