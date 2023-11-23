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


#list of arguments to pass to Stan
dat <- list(
  score=dbapm$score_s,
  core=as.integer(dbapm$core), # 1"High Grade"  2   "Low Grade"    3  "High Inf."   4   "Tum. Periphery" 5 "Tum. Border" 6   "nn"      
  apm=as.integer(dbapm$apm), # 1"LMP10" 2"TAP2"  3"B2M"   4"HC10"  5"HLAII"
  npl=as.integer(dbapm$histo), #   "UPS","MFS","DDLPS","MLPS","LMS","GIST","CHS","CH"
  id=as.integer(as.factor(dbapm$pts_id)))


#model
# We want to model the score of each patient by using Core, Id, APM, and NPL covariate
# Score has a binomial distribution that its rate depends on the our covariates
# We also model relationship of our covariates in themselves
# Each coefficient (Alpha, Beta and Delta) decided by these relationships
m_apm_intratum <- ulam(
  alist(
    score ~ dbinom( 1, p) ,
    logit(p) <- g[npl] + alpha[core,npl] + beta[id,npl] + delta[apm,npl],
    # adaptive priors non-centered
    transpars> matrix[core,8]:alpha <-
      compose_noncentered( sigma_core , L_Rho_core , z_core ),
    transpars> matrix[id,8]:beta <-
      compose_noncentered( sigma_id , L_Rho_id , z_id ),
    transpars> matrix[apm,8]:delta <-
      compose_noncentered( sigma_apm , L_Rho_apm , z_apm ),
    matrix[8,core]:z_core ~ normal( 0 , 1 ),
    matrix[8,id]:z_id ~ normal( 0 , 1 ),
    matrix[8,apm]:z_apm ~ normal( 0 , 1 ),
    # fixed priors
    g[npl] ~ normal(0,1),
    vector[8]:sigma_core ~ dexp(1),
    cholesky_factor_corr[8]:L_Rho_core ~ lkj_corr_cholesky( 2 ),
    vector[8]:sigma_id ~ dexp(1),
    cholesky_factor_corr[8]:L_Rho_id ~ lkj_corr_cholesky( 2 ),
    vector[8]:sigma_apm ~ dexp(1),
    cholesky_factor_corr[8]:L_Rho_apm ~ lkj_corr_cholesky( 2 ),
    # compute ordinary correlation matrixes from Cholesky factors
    gq> matrix[8,8]:Rho_core <<- multiply_lower_tri_self_transpose(L_Rho_core),
    gq> matrix[8,8]:Rho_id <<- multiply_lower_tri_self_transpose(L_Rho_id),
    gq> matrix[8,8]:Rho_apm <<- multiply_lower_tri_self_transpose(L_Rho_apm)
  ), data=dat , chains=4, cores = 4, iter = 4000)


# dashboard
dashboard(m_apm_intratum)

#posterior predictive simulations

post <- extract.samples(m_apm_intratum)

# First, we inspect APM effect from the posterior
# Draw it

pdf("post_pred_sim_APM.pdf",  
    width = 10, height = 16)
apm <- levels(dbapm$apm)
histo <- levels(dbapm$histo)
par(mfrow=c(length(histo),length(apm)),mar=c(4, 4.5, 2.1, 1.1)) 

for(j in 1:length(histo)){ 
  p_link <- function(a) {
    p <- with( post , g[,j] + delta[,a,j] ) 
    p <- inv_logit(p)
    return(p)
  }
  p_raw <- sapply( 1:5 , function(i) p_link( i ) ) 
  for(jj in 1:length(apm)){
    p_mu <- mean(p_raw[,jj])
    dens(p_raw[,jj], lwd = 3,   
         xlab="", ylab = 'Density', 
         main = paste(apm[jj],histo[j]), adj = 1)
    lines(c(p_mu, p_mu), c(0, 1000), lty = 2)
  }
}
dev.off()

# summary num from the model
pprob <- matrix(NA, nrow = length(apm) *length(histo), ncol = 3)
pprob <- data.frame(pprob)
for(a in 1:length(apm)){ 
  for(h in 1:length(histo)) {
    p <- with( post , g[,h] + delta[,a,h] ) 
    p <- inv_logit(p)
    idx <- ( (a-1) * length(histo) ) + h
    pprob[idx, 1] <- mean(p)
    pprob[idx, 2] <- HPDI(p)[1]
    pprob[idx, 3] <- HPDI(p)[2]
    rownames(pprob)[idx] <- paste(apm[a],histo[h],sep='_')
  }}

round(pprob,2)


# Similarly we inspect same effect, this time we add Core coefficient
# Draw it

pdf("Post_Pred_Sim_APM_core.pdf",
    width = 10, height = 16)
core <- levels(dbapm$core)[1:5] #no nn
apm <- levels(dbapm$apm)
histo <- levels(dbapm$histo)

par(mfrow=c(length(histo), 
            length(apm)),mar=c(4, 4.5, 2.1, 1.1)) 
for(jj in 1:length(histo)){
  for(j in 1:length(apm)){
    p_link_mu <- function(C) {
      mu <- with( post , g[,jj] + delta[,j,jj] + alpha[,C,jj]) 
      mu <- inv_logit(mu)
      return(mu)
    }
    p_raw <- sapply( 1:5 , function(i) p_link_mu( i ) ) 
    p_mu <- apply( p_raw , 2 , mean )
    p_ci <- apply( p_raw , 2 , HPDI )
    plot( NULL , xlab="Core" , ylab="Probability", 
          xaxt="n" , xlim=c(1,5), ylim=c(0,1),
          main = paste(apm[j],histo[jj]))
    axis( 1 , at=1:5 , 
          labels= c("HG","LG","Hi", "TP","TB"), 
          las = 2)
    lines( 1:5 , p_mu )
    shade( p_ci , 1:5 )}}
dev.off()

#contrast
p_link_mu <- function(C,A,H) {
  mu <- with( post , g[,H] + 
                delta[,A,H] + 
                alpha[,C,H]) 
  mu <- inv_logit(mu)
  return(mu)
} 
p_contrast_HGvsHI <- function(A,H){
  contr <- p_link_mu(1, A, H) - p_link_mu(3, A, H)
}

# contrasts 
# -	Contrast LMP10 DDLPS
L <- p_link_mu(1, 1, 3) - p_link_mu(3, 1, 3)
L_d <- p_contrast_HGvsHI(1,3)

# -	Contrast TAP2 ups gist
TAP_u <- p_contrast_HGvsHI(2,1)
TAP_g <- p_contrast_HGvsHI(2,6)
# -	Contrast b2 GIST
b2_g <- p_contrast_HGvsHI(3,6)

# -	Contrast  CH-chs
HC_chs <- p_contrast_HGvsHI(4,8)

x <- list(L_d, TAP_u, TAP_g, b2_g, HC_chs)
text_label <- c('LMP10 in DDLPS',
                'TAP2 in UPS', 
                'TAP2 in GIST', 
                'B2M in GIST',
                'HC10 in CHS')
pdf("contrast_HGvsHi.pdf", 
    width = 5, height = 5)
plot(NULL, xlim = c(-0.35,0.01), ylim = c(0,20), 
     main = 'Contrasts HG - Hi', 
     xlab = '', ylab = 'Density')
for(i in 1:5) dens(x[[i]], lwd = 3, 
                   col = i, add = TRUE)

legend('topleft', lwd = 3, col = 1:5, legend = text_label)
dev.off()

for(i in 1:5) {
  mu <- mean(x[[i]])
  print(round(mu,2))}

for(i in 1:5) {
  hpdi <- HPDI(x[[i]])
  print(round(hpdi,2))}



pdf("Post_Pred_Sim_APM_core_10ptv2.pdf", 
    width = 15, height = 20)
apm <- levels(dbapm$apm)
histo <- levels(dbapm$histo)

x <- rep(c(1,25, 25, 1), times = 63)
y <- rep(c(1, 1, 2, 2), times = 63) 
y <- c(1, 1, 2, 2) + rep(seq(0,124,2),each = 4)
plot(NULL, xlim = c(1,25), ylim = c(0,130),
     xlab="" , ylab="Patients", 
     yaxt="n", xaxt="n", bty = 'L', 
     main = 'Interpatient variability')
polygon(x,y, lty = 0, col = alpha(1, 0.05))
abline(h = 1:126 - 0.5 , lty = 2, lwd = 0.5)
abline(v = 1:4 * 5 + 0.5 , lty = 1, lwd = 1)
abline(v = 1:25, lty = 2, lwd = 0.5)
core_lab <- c("HG","LG","Hi", "TP","TB")
axis(1,1:25, rep(core_lab, times = 5))
text(1:5 * 5 - 2.5, rep(127, times = 5),
     labels = apm)
legend('top', legend = histo,
       col = 1:8, lwd = 2, 
       horiz = TRUE, bg = 'white')

for(jj in 1:length(histo)){
  ids <- dat$id[dat$npl == jj]
  for(j in 1:length(apm)){
    for (id in ids){
      p_link_mu <- function(C) {
        mu <- with( post , g[,jj] + 
                      delta[,j,jj] +
                      alpha[,C,jj] + 
                      beta[,id,C]) 
        mu <- inv_logit(mu)
        return(mu)
      }
      p_raw <- sapply( 1:5 , function(i) p_link_mu( i ) ) 
      for(i in 1:10) lines( 1:5 + (5 * (j-1)) , 
                            126 - id + p_raw[i,] , 
                            lwd = 0.2, col = jj)
    }}}

dev.off()