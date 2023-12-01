# APM_Sarcoma
Analysis of APM (antigen presenting and processing machinery) presence in Sarcoma

# Introduction

Sarcomas, a diverse range of malignant tumors, pose challenges due to limited treatment options and unfavorable outcomes. The antigen processing and presenting machinery (APM) is pivotal in creating tumor-specific antigens for immune recognition and targeting. Yet, the APM status in sarcomas remains poorly understood.

## Material and Methods

We approached the analysis of APM defects by modeling it as a binomial distribution represented by a score. This score is influenced by the interplay of NPL with Core, Id, and APM. The entire network was conceptualized as a hierarchical Bayesian model, wherein parameters were estimated for NPL, Core and NPL, patient ID and NPL, and finally APM and NPL. Samples for NPL followed a normal distribution with a mean of 0 and a variance of 1. Except for NPL, other samples in this model were drawn from a multivariate normal distribution. The mean for each relationship in this distribution was set to 0, and the covariance matrix was determined by the variance-covariance matrix for each relationship. Histotypes, denoted as 8, were considered in the analysis, resulting in a covariance matrix for these relationships with a shape of 8x8. The correlation factor (R value) for each matrix was sampled from the LKJcorr distribution, with a chosen value of 2 to prevent divergent transitions.


```math
\begin{align}
score &\sim \text{dbinom}(1, p)\\
Logit(score) &=\gamma_{NPL}+\alpha_{CORE,NPL}+\beta_{ID,NPL}+\delta_{APM,NPL}\\

\begin{bmatrix}
\alpha_{j,1}\\
\vdots\\
\alpha_{j,8}\\
\end{bmatrix}
&\sim
MVNormal\left(\begin{bmatrix}
0\\
\vdots\\
0\\
\end{bmatrix}, H_\alpha
\right)\nonumber\\
H_\alpha&=
\begin{pmatrix}
    \sigma_{\alpha1} & 0 & \dots & 0 \\
    0 & \sigma_{\alpha2} & \dots & 0 \\
    \vdots & \vdots & \ddots & \vdots \\
    0 & 0 & \dots & \sigma_{\alpha8}
  \end{pmatrix}
R_\alpha
\begin{pmatrix}
    \sigma_{\alpha1} & 0 & \dots & 0 \\
    0 & \sigma_{\alpha2} & \dots & 0 \\
    \vdots & \vdots & \ddots & \vdots \\
    0 & 0 & \dots & \sigma_{\alpha8}
  \end{pmatrix}\nonumber \\

\begin{bmatrix}
\beta_{j,1}\\
\vdots\\
\beta_{j,8}\\
\end{bmatrix}
&\sim
MVNormal\left(\begin{bmatrix}
0\\
\vdots\\
0\\
\end{bmatrix}, H_\beta
\right)\nonumber\\
H_\beta&=
\begin{pmatrix}
    \sigma_{\beta1} & 0 & \dots & 0 \\
    0 & \sigma_{\beta2} & \dots & 0 \\
    \vdots & \vdots & \ddots & \vdots \\
    0 & 0 & \dots & \sigma_{\beta8}
  \end{pmatrix}
R_\beta
\begin{pmatrix}
    \sigma_{\beta1} & 0 & \dots & 0 \\
    0 & \sigma_{\beta2} & \dots & 0 \\
    \vdots & \vdots & \ddots & \vdots \\
    0 & 0 & \dots & \sigma_{\beta8}
  \end{pmatrix}\nonumber \\

\begin{bmatrix}
\delta_{j,1}\\
\vdots\\
\delta_{j,8}\\
\end{bmatrix}
&\sim
MVNormal\left(\begin{bmatrix}
0\\
\vdots\\
0\\
\end{bmatrix}, H_\delta
\right)\nonumber\\
H_\delta&=
\begin{pmatrix}
    \sigma_{\delta1} & 0 & \dots & 0 \\
    0 & \sigma_{\delta2} & \dots & 0 \\
    \vdots & \vdots & \ddots & \vdots \\
    0 & 0 & \dots & \sigma_{\delta8}
  \end{pmatrix}
R_\delta
\begin{pmatrix}
    \sigma_{\delta1} & 0 & \dots & 0 \\
    0 & \sigma_{\delta2} & \dots & 0 \\
    \vdots & \vdots & \ddots & \vdots \\
    0 & 0 & \dots & \sigma_{\delta8}
  \end{pmatrix}\nonumber \\
\gamma_{NPL} &\sim Normal(0,1)\nonumber\\
R_\alpha &\sim LKJcorr (2)\\
R_\beta &\sim LKJcorr (2)\\
R_\delta &\sim LKJcorr (2)\\
\sigma_{\alpha j} &\sim Exponential(1)  \text{ \tiny{for $j=1..8$}}\nonumber\\
\sigma_{\beta j} &\sim Exponential(1)  \text{ \tiny{for $j=1..8$}}\nonumber\\
\sigma_{\delta j} &\sim Exponential(1)  \text{ \tiny{for $j=1..8$}}\nonumber\\
\end{align}
```

We conducted a survival analysis by modeling the time-to-event (T_i) as an exponential distribution and its log-complementary cumulative density function, akin to Cox survival analysis. The expected rate (λ_i) was transformed into the expected value μ_i, and a log link function was applied to the linear model, specifically an intercept-only model. The parameter α + β was modeled hierarchically using a Bayesian approach, estimating parameters for the APM status and histotype, and for the grading and histotype. These parameters were sampled from a multivariate normal distribution with a mean set to 0, and the covariance matrix was determined by the variance-covariance matrix of the presence of APM (and the grade) and the 8 histotypes denoted as H. Given a total of 8 sampled histotypes, the covariance matrix had a shape of 8x8. The correlation factor within the matrix was determined by R, which was sampled from the LKJcorr distribution with a chosen shape parameter of 2. To alleviate divergent transitions, we employed the non-centered version of the model.


```math
\begin{align}
T_{i} \mid Survival_{i} = 1 &\sim Exponential(\lambda_{i})\\
T_{i} \mid Survival_{i} = 0 &\sim ExponentialLCCDF(\lambda_{i})\\
\lambda_{i} &= 1/\mu_i\\
Exponential(\mu_i)&=\alpha_{APM,HIST}+\beta_{GRADE,HIST}\\
\begin{bmatrix}
\alpha_{j,1}\\
\vdots\\
\alpha_{j,8}\\
\end{bmatrix}
&\sim
MVNormal\left(\begin{bmatrix}
0\\
\vdots\\
0\\
\end{bmatrix}, H_\alpha
\right)\nonumber\\
H_\alpha&=
\begin{pmatrix}
    \sigma_{\alpha1} & 0 & \dots & 0 \\
    0 & \sigma_{\alpha2} & \dots & 0 \\
    \vdots & \vdots & \ddots & \vdots \\
    0 & 0 & \dots & \sigma_{\alpha8}
  \end{pmatrix}
R_\alpha
\begin{pmatrix}
    \sigma_{\alpha1} & 0 & \dots & 0 \\
    0 & \sigma_{\alpha2} & \dots & 0 \\
    \vdots & \vdots & \ddots & \vdots \\
    0 & 0 & \dots & \sigma_{\alpha8}
  \end{pmatrix}\nonumber \\
\begin{bmatrix}
\beta_{j,1}\\
\vdots\\
\beta_{j,8}\\
\end{bmatrix}
&\sim
MVNormal\left(\begin{bmatrix}
0\\
\vdots\\
0\\
\end{bmatrix}, H_\beta
\right)\nonumber\\
H_\beta &=
\begin{pmatrix}
    \sigma_{\beta1} & 0 & \dots & 0 \\
    0 & \sigma_{\beta2} & \dots & 0 \\
    \vdots & \vdots & \ddots & \vdots \\
    0 & 0 & \dots & \sigma_{\beta8}
  \end{pmatrix}
R_\beta
\begin{pmatrix}
    \sigma_{\beta1} & 0 & \dots & 0 \\
    0 & \sigma_{\beta2} & \dots & 0 \\
    \vdots & \vdots & \ddots & \vdots \\
    0 & 0 & \dots & \sigma_{\beta8}
  \end{pmatrix}\nonumber\\
R_\alpha &\sim LKJcorr(2)\\
R_\beta &\sim LKJcorr(2)\\
\sigma_{\alpha j} &\sim Exponential(1) \text{\tiny{for $j=1..8$}}\nonumber\\
\sigma_{\beta j} &\sim Exponential(1) \text{\tiny{for $j=1..8$}}\nonumber\\
\end{align}
```
## Results
Our study uncovered APM defects in all the investigated sarcomas, with HLA Class I β2-microglobulin and HLA Class II being the least affected. Noteworthy, histotype-specific defects include LMP10 in LMS, MLPS, and DDLPS, and TAP2 in UPS, GIST, and CH. Spatial variations in APM expression were observed, with high-grade areas exhibiting distinct patterns. Patient-level heterogeneity also was noted, with certain individuals prone to APM component deficiencies. Importantly, loss of any APM component predicted distant metastasis in LMS and DDLPS, and overall survival in LMS.

## Conclusion
Study reveals pronounced defects in Antigen Presenting Machinery (APM) components within sarcomas, showcasing notable variations across histotypes and tumoral regions. Key alterations involve HLA Class I subunit β2-microglobulin, HLA Class I subunit α (HC10), and MHC I transporting unit TAP2. Notably, the absence of APM components serves as a prognostic indicator for distant metastases and overall survival in LMS and DDLPS. This underscores the clinical significance of APM dysregulation in these sarcoma subtypes, providing crucial insights into molecular mechanisms and informing personalized therapeutic strategies.

## Repository Structure

`Survival Analysis.R` This file contains the Kaplan Meier Survival Analysis. It produces two Figures. They corresponds to Figure 5A and Figure 5B in the paper. 

`Overall Survival Hazard Ratio.R` R file to run the OS_Grade.stan model. It will produce the Supplementary Figure 4

`Distant Metastasis Hazard Ratio.R` R file to run the DFS_Grade.stan model. It will produce the Supplementary Figure 3

`APM_Presence_Spt_Differences_Interpatient_Variability.R` First analysis part of the paper. It will produce the Figure 1,2,3,4.

`OS_Grade.stan` This file contains the model for Overall Survival according to APM's presence on Histotypes

`DFS_Grade.stan` This file contains the model for Disease Free Survival according to APM's presence on Histotypes

## Dataset

The datasets used in this project can be downloaded from Zenodo at the following link: ZENODO

## Requirements

The following softwares were used for developing and running this project:

* R version 4.1.2
* Stan version 2.21.0


## License

This project is licensed under the Creative Commons Attribution 4.0 International.
