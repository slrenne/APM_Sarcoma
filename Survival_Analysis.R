library(tidyverse)
library(lubridate)
library(survival)
library(survminer)
library(rethinking)

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

apm <- ifelse(  dbsurv[,4] < 0.6 | #create a vector value for apm
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

## This is distant metastasis
km_fit <- surv_fit(Surv(time,Distant_metastasis)~APM,data=dbsurv)
ggsurvplot_facet(fit = km_fit,
                 data = dbsurv,
                 xlab = "Months", 
                 ylab = "Disease free survival",
                 https://urldefense.proofpoint.com/v2/url?u=http-3A__facet.by&d=DwIGAg&c=5rLNXN0mp_7LMh3Fds96xpjyD06ZuE2RU7zikolS0lg&r=y76TMtVUYsFGbfDq9nJS0WbZeQ76UT6_9yhZbHeshcQ&m=M0ganPrQLcx2yUW9rUOmnXsCxjHbYMw-yfXRN0O7X8URCKqWP8lFhgZj0X77kQ93&s=UYFDX5hUbpYFOaLtwgcxcHmPKiC8X7dhO1e4SJaFEkY&e=  = "histo", 
                 legend.labs = c("APM-", "APM+"), 
                 legend.position= 'top',
                 short.panel.labs = TRUE,
                 https://urldefense.proofpoint.com/v2/url?u=http-3A__conf.int&d=DwIGAg&c=5rLNXN0mp_7LMh3Fds96xpjyD06ZuE2RU7zikolS0lg&r=y76TMtVUYsFGbfDq9nJS0WbZeQ76UT6_9yhZbHeshcQ&m=M0ganPrQLcx2yUW9rUOmnXsCxjHbYMw-yfXRN0O7X8URCKqWP8lFhgZj0X77kQ93&s=YTtjFDmjipKpPHGxi3VxqZz-EPwMdjiFdSICJ98j_-w&e=  = FALSE,
                 nrow = 2, 
                 palette = c(2,3))

ggsave('Disease_Free_Survival.png',width = 9, height = 6)

## plotting the kaplan meier
## This is overall survival
km_fit <- surv_fit(Surv(last_fup_m,Status) ~ APM,data=dbsurv)

ggsurvplot_facet(fit = km_fit,
                 data = dbsurv,
                 xlab = "Months",
                 ylab = "Overall survival",
                 https://urldefense.proofpoint.com/v2/url?u=http-3A__facet.by&d=DwIGAg&c=5rLNXN0mp_7LMh3Fds96xpjyD06ZuE2RU7zikolS0lg&r=y76TMtVUYsFGbfDq9nJS0WbZeQ76UT6_9yhZbHeshcQ&m=M0ganPrQLcx2yUW9rUOmnXsCxjHbYMw-yfXRN0O7X8URCKqWP8lFhgZj0X77kQ93&s=UYFDX5hUbpYFOaLtwgcxcHmPKiC8X7dhO1e4SJaFEkY&e=  = "histo", 
                 legend.labs = c("APM-", "APM+"), 
                 legend.position= 'top',
                 short.panel.labs = TRUE,
                 https://urldefense.proofpoint.com/v2/url?u=http-3A__conf.int&d=DwIGAg&c=5rLNXN0mp_7LMh3Fds96xpjyD06ZuE2RU7zikolS0lg&r=y76TMtVUYsFGbfDq9nJS0WbZeQ76UT6_9yhZbHeshcQ&m=M0ganPrQLcx2yUW9rUOmnXsCxjHbYMw-yfXRN0O7X8URCKqWP8lFhgZj0X77kQ93&s=YTtjFDmjipKpPHGxi3VxqZz-EPwMdjiFdSICJ98j_-w&e=  = FALSE,
                 nrow = 2, 
                 palette = c(2,3))

ggsave('Overall_Survival.png',width = 9, height = 6)