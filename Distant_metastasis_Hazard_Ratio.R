library(tidyverse)
library(lubridate)
library(survival)
library(survminer)
library(rethinking)
library(rstan)
library("writexl")

# import
dbapm <- read.csv('./dbapm.csv')
clin <- read.csv('./clin.csv')
dbapm$histo <- factor(dbapm$histo, 
                      levels = c("UPS","MFS",
                                 "DDLPS","MLPS",
                                 "LMS","GIST",
                                 "CHS","CH"))
dbapm$core <- factor(dbapm$core, 
                     levels =  c("High Grade",
                                 "Low Grade", 
                                 "High Inf.", 
                                 "Tum. Periphery", 
                                 "Tum. Border"))
dbapm$apm <- factor(dbapm$apm, 
                    levels = c("LMP10","TAP2", 
                               "B2M", "HC10", 
                               "HLAII") )

str(dbapm)

#########################

# we will model survival of patients
# we expect that loss of any apm component could affect Ag presentation
# since several areas were sampled we will limit the analysis on the HG

hg <- dbapm$core == 'High Grade'
dbsurv <- dbapm[hg,] #keep only HG cores

dbsurv <- dbsurv %>% 
  spread(apm, score_s) #spread the apm columns

dbsurv <- dbsurv %>% 
  group_by(pts_id,histo,Grade) %>%
  summarise(LMP10 = mean(LMP10, na.rm = TRUE), #take the mean for each APM value
            TAP2 = mean(TAP2, na.rm = TRUE),   # i.e. only pts with loss in all the HG cores will make 0
            B2M = mean(B2M, na.rm = TRUE),
            HC10 = mean(HC10, na.rm = TRUE),
            HLAII = mean(HLAII, na.rm = TRUE)
  )

apm <- ifelse(dbsurv[,4] < 0.6 | #create a vector value for apm
              dbsurv[,5] < 0.6 | 
              dbsurv[,6] < 0.6 | 
              dbsurv[,7] < 0.6 |
              dbsurv[,8] < 0.6 , 0,1)


apm[is.na(apm)] <- 1  # cases will have 'NA' if 'NaN & !=0'. makes them 1.
attr(apm, 'dimnames') <- NULL
dbsurv$APM <- as.factor(apm[,1])

dbsurv$Grade[is.na(dbsurv$Grade)] <- 2 # Chordoma has no grading system
str(dbsurv$Grade)

dbsurv <- left_join(dbsurv,clin, by = "pts_id")

dbsurv$time <- with(dbsurv,{ifelse(Distant_metastasis == 1,time_DM_m,last_fup_m)})


dat <- with(dbsurv, {list(
  DM = Distant_metastasis, 
  A = as.integer(APM),
  D = as.integer(histo), #   "UPS","MFS","DDLPS","MLPS","LMS","GIST","CHS","CH"
  G = as.integer(Grade),
  months_to_event = time)})

file <- file.path("DFS_grade.stan")
mod <- stan_model(file)

fit <- sampling(mod,data = dat,chains = 4)
# We show the summary of our experiment
summary(fit)
# We extract the sampling results
list_of_draws <- extract(fit)
# We need alphas to check the APm effect on histotypes
# med_s is 4000X2X8 means we have 4000 samples and each of them 2x8 matrices
# This is because we have 8 histotypes and for each of them we have
# 1 (apm exists) and 0 (apm does not exist)
med_s <- list_of_draws$alpha

# Since we have 4000 draws and 8 histotypes to control
med_s_diff <- matrix(NA, nrow = 4000, ncol = 8)

# Then we calculate the hazard ratio
# Information at the link below
# https://urldefense.proofpoint.com/v2/url?u=https-3A__bookdown.org_sestelo_sa-5Ffinancial_computing-2Dthe-2Dhazard-2Dratio.html&d=DwIGAg&c=5rLNXN0mp_7LMh3Fds96xpjyD06ZuE2RU7zikolS0lg&r=y76TMtVUYsFGbfDq9nJS0WbZeQ76UT6_9yhZbHeshcQ&m=M0ganPrQLcx2yUW9rUOmnXsCxjHbYMw-yfXRN0O7X8URCKqWP8lFhgZj0X77kQ93&s=2e_Xew36L0mu5TolTqHUnheOfdaLBPYDqSqB6u6Fk5g&e= 
# the med_s matrix was constructed as given below
#Example
#----                                                             -----
#|APM_zero_for_histo_1,APM_zero_for_histo_2,APM_zero_for_histo_3 .... |
#|APM_one_for_histo_1,APM_one_for_histo_2,APM_one_for_histo_3 .....   |
#----                                                             -----
# Since we want to know the effect of histo type which apm exists
# Over the histo type that apm does not exist
# We subtract first row from the second row of our result matrix
for(i in 1:8) med_s_diff[,i] <-  exp(med_s[,2,i]) / exp(med_s[,1,i])

# precis(as.data.frame(exp(med_s_diff)))
# And we save the differences in a dataframe

differences_df <- as.data.frame(med_s_diff)

m <- precis(as.data.frame(med_s_diff))@.Data[[1]]
l <- precis(as.data.frame(med_s_diff))@.Data[[3]]
h <- precis(as.data.frame(med_s_diff))@.Data[[4]]
coef_alpha <- data.frame(m,l,h, row.names = levels(dbapm$histo))
coef_alpha[ , 'Names'] <- c('UPS','MFS','DDLPS','MLPS','LMS','GIST','CHS','CH') 
write_xlsx(coef_alpha, "Values_for_CI_DFS.xlsx")
# plot(precis(as.data.frame(med_s_diff)))
# abline(v=1,lty=2)

# Since will graphs violin plots, first we change the column names into
# Original histo type names and save the dataframe into a 'DM_Graded.csv' file

colnames(differences_df) <- c('UPS','MFS','DDLPS','MLPS','LMS','GIST','CHS','CH') 
write.csv(differences_df, "DM_Graded.csv", row.names=FALSE)
med_s_diff <- read.csv("DM_Graded.csv")


colnames(med_s_diff) <- c('UPS','MFS','DDLPS','MLPS','LMS','GIST','CHS','CH') 
col_names = colnames(med_s_diff)
for (i in 1:8) cat(col_names[i],median(med_s_diff[,i]))
# png('Histogram_of_Graded_DM.pdf', width = 800, height = 4000)
# par(mfrow=c(8,1))
# for(i in 1:8) {
#   rethinking::dens(med_s_diff[,i] , xlim= c(0,44),bty = 'L', lwd = 2, col = i, adj = 1)
#   abline(v=1,lty = 2)
# }
# dev.off()
pdf('Histogram_of_Graded_DM.pdf', width = 10, height = 20)
par(mfrow=c(8,1))
for(i in 1:8) {
  rethinking::simplehist(med_s_diff[,i],xlab="",xlim= c(0,30), lwd = 4,col = i, main = col_names[i])
  #rethinking::dens(med_s_diff[,i] , xlim= c(0,44),bty = 'L', lwd = 2, col = i, adj = 1, add = TRUE)
  abline(v=1,lty = 2)
}
dev.off()

