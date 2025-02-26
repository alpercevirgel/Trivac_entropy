
### libraries
packages = c("tidyverse", "ggpubr", "reshape2","gghalves", "ggalluvial","rstatix", "readxl","Hmisc","ordinal","writexl","mgcv","gghalves","gtsummary","sjPlot")

## Now load or install&load all
package.check <- lapply(packages,FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)}})

### %in& negate
`%!in%` = Negate(`%in%`)


dir.create(file.path(analysis.path, "results"), showWarnings = FALSE)
dir.create(file.path(analysis.path, "results/figures"), showWarnings = FALSE)
dir.create(file.path(analysis.path, "results/tables"), showWarnings = FALSE)
dir.create(file.path(analysis.path, "results/figures/validation"),showWarnings = FALSE)
dir.create(file.path(analysis.path, "results/figures/trivaccine"),showWarnings = FALSE)
dir.create(file.path(analysis.path, "results/figures/models"),showWarnings = FALSE)
dir.create(file.path(analysis.path, "results/figures/entropy"),showWarnings = FALSE)
dir.create(file.path(analysis.path, "results/figures/cytokine"),showWarnings = FALSE)

cytokine.dir <- file.path(analysis.path, "results/figures/cytokine/")
trivaccine.dir <- file.path(analysis.path, "results/figures/trivaccine/")
models.dir <- file.path(analysis.path, "results/figures/models/")
entropy.dir <- file.path(analysis.path, "results/figures/models/")


# colorblind friendly colors for clusters
cols_cluster <- c("1"= "#77AADD", "2"= "#99DDFF",
                  "3"= "#44BB99", "4"= "#BBCC33",
                  "5"= "#AAAA00", "6"= "#EEDD88",
                  "7"= "#EE8866", "8"= "#FFAABB", 
                  "9"= "#DDDDDD")

cols_model <- c("Immunotype 1"= "#77AADD", "Immunotype 2"= "#99DDFF","Immunotype 3"= "#44BB99", 
                "Immunotype 4"= "#BBCC33","Immunotype 5"= "#AAAA00", "Immunotype 6"= "#EEDD88",
                "Immunotype 7"= "#EE8866", "Immunotype 8"= "#FFAABB", "Immunotype 9"= "#DDDDDD", 
                "Day 0 HI Q1" = "#D2D2D2", "Day 0 HI Q2" = "#A0A0A0", "Day 0 HI Q3" = "#6E6E6E", "Day 0 HI Q4" = "#3C3C3C")

# colorblind friendly colors for age groups
cols_agegroup <- c(`25-49` = "#F0E442", `50-64` = "#85C0F9", `65-98` = "#F5793A")


# boxplot with log2 y axis
boxplot_cluster <- function()
{
  p <- ggplot(df_loop, aes(x=Xaxis, y=Yaxis)) +
    geom_boxplot(aes(fill=Groups), alpha=0.9,outlier.size=0,outlier.colour="white") +
    geom_jitter(aes(fill=Groups), alpha=0.4, width = 0.3, shape=21,size=1) +
    scale_fill_manual("Immunotypes", values=cols_cluster) +
    theme_classic()+
    stat_summary(aes(y=Yaxis, x=Xaxis),size=0.2)
}

# Function to fit models and store AIC, BIC, and gam.check() results
fit_and_compare_gams <- function(data, response, predictor, k_start, k_end) {
  results <- data.frame(k = integer(), AIC = numeric(), BIC = numeric())
  
  for (k in k_start:k_end) {
    # Fit the GAM model with the current value of k
    gam_model <- gam(formula = as.formula(paste(response, "~ s(", predictor, ", k=", k, ")")), 
                     data = data, method = "ML")
    
    # Extract the AIC and BIC
    model_aic <- AIC(gam_model)
    model_bic <- BIC(gam_model)
    
    # Store the results
    results <- rbind(results, data.frame(k = k, AIC = model_aic, BIC = model_bic))
  }
  
  return(results)
}
