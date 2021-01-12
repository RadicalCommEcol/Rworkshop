
# Adaptation of the code for chapter 4 in
# Beginner's Guide to GLM and GLMM with R by
# Alain Zuur, Joseph M Hilbe, and Elena N Ieno
# www.highstat.com

###################################################################
###################################################################
#Load the data

Spiders_raw <- read.table(file = "lmm-glmm/0 Raw_Data/Spiders.txt", 
                       header = TRUE, 
                       dec = ".")
names(Spiders_raw)
str(Spiders_raw)
###################################################################
###################################################################
#Load packages and library files
library(tidyverse)
library(lattice)
library(lme4) # LMMs
library(glmmTMB) # for GMMs
library(usdm) #for VIF
##################################################################
##################################################################


# We remove some values for demonstration purposes:
# Now there are plots with fewer observations

Spiders_raw <- Spiders_raw[Spiders_raw$Hlog10>0.62,]

# Create a factor for Plot: fPlot. It will be needed later on in the models
Spiders_raw$fPlot <- factor(Spiders_raw$Plot)


# We filter some plots, as Zuur et al. did.
Spiders <- Spiders_raw %>% filter(!fPlot %in% c("4","9","11","14","23"))

# Create a factor for Plot: fPlot. Now it has the correct number of levels (25)
Spiders$fPlot <- as.factor(as.numeric(Spiders$Plot))

##################################################################
# Visualizing the data

hist(Spiders$Hlog10)

ggplot(Spiders, aes(x = HerbLayer, y = Hlog10))+
  geom_point(alpha=0.5)+
  labs(x = "Percentage of herb layer cover",
       y = "Shannon index")

ggplot(Spiders, aes(x = HerbLayer, y = Hlog10))+
  geom_point(alpha=0.5)+
  facet_wrap(~fPlot)+
  labs(x = "Percentage of herb layer cover",
       y = "Shannon index")

##################################################################
# Linear regression

# We scale the covariable HerbLayer
Spiders$HerbLayerc <- scale(Spiders$HerbLayer)

M0 <- lm(Hlog10 ~ HerbLayerc, data = Spiders)

summary(M0)

ggplot(Spiders, aes(x = HerbLayerc, y = Hlog10))+
  geom_point(alpha=0.5)+
  geom_smooth(method = "lm")+
  labs(x = "(Scaled) Percentage of herb layer cover",
       y = "Shannon index")
  
# Checking model assumptions

# Homogeneity of variance

plot(M0, which = 1)  # Look alright
# ggplot way:
# ggplot(M0) + geom_point(aes(x = .fitted, y = .resid))+
#   geom_hline(aes(yintercept=0))+
#   geom_smooth(aes(x = .fitted, y = .resid), se = F)

# Normality of residuals
# qqplot - point should ideally fall onto the diagonal dashed line

plot(M0, which = 2)  # a bit off at the extremes, but looks OK
summary(M0)$sigma
# NOTE: Residuals deviate from the diagonal line in the upper tail.
# This plot indicated that the right tail is ‘lighter’ (have smaller values) than 
# what we would expect under the standard modeling assumptions.

# Compare the residuals with the normal distribution with its variance
library(ggResidpanel)
ggResidpanel::resid_panel(M0, plots = "hist")

# Check for influential observations
plot(M0, which = 5)

# Leverage values are the diagonal elements of hat matrix: y_hat = H %*% y
# H_ii determines the influence of y_i on y_hat_i. High H_ii indicates a potential
# influencial observation.

# One possible definition is that an outlier is any point that
# isn’t approximated well by the model (has a large residual) 
# and which significantly influences model fit (has large leverage).


# what about observation independence?
# Residuals VS covariate "Plot"

ggplot(Spiders %>% mutate(residuals=resid(M0)), 
       aes(x = fPlot, y = residuals)) + geom_boxplot()+
  geom_hline(aes(yintercept=0))+
  theme_bw()+
  labs(x="Plot",y="Residuals")

# Box plot suggest a Plot effect in the residuals,
# that can be visualize also as follows 

ggplot(Spiders,aes(x=HerbLayerc,y=Hlog10))+
  geom_point(show.legend = FALSE,alpha=0.5) +
  geom_abline(aes(intercept = coef(M0)[[1]], slope = coef(M0)[[2]],
                  color = 'black'),size = 1) + 
  facet_wrap(~fPlot)+
  scale_color_identity(labels=c("pooled"), guide="legend")+
  labs(x = "(Scaled) percentage of herb layer",y = "Shannon index",color="Model fit")

# To verify the Plot effect, we apply a linear
# regression in which we model the residuals as a function of Plot.

R0 <- rstandard(M0)
Spiders$R0 <- R0

Test0 <- lm(R0 ~ fPlot, data = Spiders)
drop1(Test0, test = "F")
# The drop1 function applies the full model and drops each term in turn, and, each
# time, an F test is applied to compare the full model and the nested model.
#PR is smaller than the significance level.
# We should abandon the M0 model.

####################################################
# Non-pooled model

# Plot is a categorical variable with 25 levels -> It requires 24 regression parameters
# Each of the latter is used as a correction of the intercept for a particular slope.

M1 <- lm(Hlog10 ~ HerbLayerc + fPlot, 
         data = Spiders)
summary(M1)
drop1(M1, test = "F") # Plot is significant at 5% level

length(coef(M1)) # There are 26 coefficients!!

Spiders$fitted_M1 <- fitted(M1)

ggplot(Spiders ,aes(x=HerbLayerc,y=Hlog10))+
  geom_point(show.legend = FALSE,alpha=0.5) +
  facet_wrap(~fPlot)+
  geom_abline(aes(intercept = coef(M0)[[1]], slope = coef(M0)[[2]],
                  color = 'black'),size = 1) + 
  geom_line(aes(x=HerbLayerc,y=fitted_M1, color = "red"),size=1)+
  scale_color_identity(labels=c("pooled","non-pooled"), guide="legend")+
  labs(x = "(Scaled) percentage of herb layer",y = "Shannon index",color="Model fit")


#####################################################
# Partially Pooled model
M2 <- lmer(Hlog10 ~ HerbLayerc + (1 | fPlot), 
           data = Spiders)
summary(M2)

Spiders$fitted_M2 <- fitted(M2)

ggplot(Spiders ,aes(x=HerbLayerc,y=Hlog10))+
  geom_point(show.legend = FALSE,alpha=0.5) +
  facet_wrap(~fPlot)+
  geom_abline(aes(intercept = coef(M0)[[1]], slope = coef(M0)[[2]],
                  color = 'black'),size = 1) + 
  geom_line(aes(x=HerbLayerc,y=fitted_M1, color = "red"),size=1)+
  geom_line(aes(x=HerbLayerc,y=fitted_M2, color = "blue"),size=1)+
  scale_color_identity(breaks = c("black", "red", "blue"),
                       labels=c("pooled","non-pooled",
                                "partially pooled"), guide="legend")+
  labs(x = "(Scaled) percentage of herb layer",y = "Shannon index",color="Model fit")


#library(sjPlot)

# Visualise random effects 
sjPlot::plot_model(M2, type = "re", show.values = TRUE)


# Obtaining a p-value based on a z-distribution, following Zuur et al.
Betas  <- fixef(M2)                  #Get the betas
SE     <-  sqrt(diag(vcov(M2)))      #Get the SEs
pval   <- 2*pnorm(-abs(Betas  / SE)) #Z distribution
Output <- cbind(Betas,SE, pval)
print(Output, digits = 3)

M2_TMB <- glmmTMB(Hlog10 ~ HerbLayerc + (1 | fPlot), 
               data = Spiders) # The estimates of p-values are quite close
# to the pevious ones.


#######################################################
# Full linear mixed model

# We add two more (scaled) variables as fixed effects
Spiders$GroundVegc <- scale(Spiders$GroundVeg)
Spiders$Litterc    <- scale(Spiders$Litter)

# We check that the variables are not collinear
usdm::vif(as.data.frame(dplyr::select(Spiders,HerbLayerc,GroundVegc,
                                      Litterc))) #OK!!

M3 <- lmer(Hlog10 ~ HerbLayerc + GroundVegc + Litterc + (1 | fPlot), 
           data = Spiders)

summary(M3)
# Note that correlation of Fixed Effects is about the expected correlation of 
# the regression coefficients (betas).
# So, it is not about the correlation of the variables!!

E3 <- resid(M3)
F3 <- fitted(M3) #provides y_hat


Betas     <- fixef(M3) # provides betas
X         <- model.matrix(M3) # provides model matrix
FitManual <- X %*% Betas # Fixed part of the y

RE <- ranef(M3)$fPlot$'(Intercept)' #contains the random intercepts of each plot (25)
AllRE <- RE[as.numeric(Spiders$fPlot)]# random intercepts of each observation (138)

FitManual + AllRE - fitted(M3) #Sanity check

ResRaw <- Spiders$Hlog10 - FitManual - AllRE
resid(M3) - ResRaw #Sanity check

summary(M3) 

# GroundVegc seems to be not significant

# Alternative method: LogLikelyhood ratio test
# REML assumes that the fixed effects structure is correct. 
# Use maximum likelihood when comparing models with different fixed effects, 
# as ML doesn’t rely on the coefficients of the fixed effects.


M4 <- lmer(Hlog10 ~ HerbLayerc + GroundVegc + Litterc + (1 | fPlot), 
           data = Spiders, REML = FALSE)

M4A <- update(M4, .~. - HerbLayerc)
M4B <- update(M4, .~. - GroundVegc)
M4C <- update(M4, .~. - Litterc)

anova(M4, M4A) # The larger the diff. between the log-Lik. values, the less likely
# it is that the null hyp. is true (loglik. values of the models are identical)
anova(M4, M4B)
anova(M4, M4C)

# Same analysis using the fucntion drop1
drop1(M4, test = "Chi") #diff~Xi^2_p where p is the diff. betw. the number of
# parameters in each model.

# These p-values are only reliable when the sample size - the number of parameters
# is larger than 50.

AIC(M4, M4A, M4B, M4C)

# Even though you use ML to compare models, you should report parameter estimates 
# from your final “best” REML model, as ML may underestimate variance of the 
# random effects.

#######################
# Final model

M5 <- lmer(Hlog10 ~ HerbLayerc + Litterc + (1 | fPlot), 
           data = Spiders, REML = TRUE)

summary(M5)

# Checking the model

res_M5 <- simulateResiduals(fittedModel = M5, n = 1000)
plot(res_M5)
plotResiduals(res_M5, Spiders$HerbLayerc)
plotResiduals(res_M5, Spiders$Litterc)
plotResiduals(res_M5, Spiders$GroundVegc)
plotResiduals(res_M5, Spiders$fPlot)

# Checking the variability of residuals 
Spiders$resid_M5 <- resid(M5)


# ggplot(Spiders %>% mutate(residuals=resid(M0)), 
#        aes(x = fPlot, y = residuals)) + geom_boxplot()
# ggplot(Spiders %>% mutate(residuals=resid(M2)), 
#        aes(x = fPlot, y = residuals)) + geom_boxplot()
ggplot(Spiders %>% mutate(residuals=resid(M5)), 
       aes(x = fPlot, y = residuals)) + geom_boxplot()

TestM5 <- lm(resid_M5 ~ fPlot, data = Spiders)
drop1(TestM5, test = "F") # fPlot is not significant -> OK!

r2(M5)

# VISUALIZATION M5
library(visreg)
library(RColorBrewer)

# Get discrete colors for our 25 fPlot from Brewer
# (https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r)

n <- length(Spiders$fPlot %>% unique())
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
color_vis=sample(col_vector, n)

# Add alpha values to colors
## Add an alpha value to a colour
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

# Labels for plot panels
plot_labs <- Spiders$fPlot %>% as.character() %>% unique()
plot_labs <- paste0("Plot ",plot_labs)
names(plot_labs) <- Spiders$fPlot %>% as.character() %>% unique()

M5_TMB <- glmmTMB(Hlog10 ~ HerbLayerc + Litterc + (1 | fPlot), 
           data = Spiders, REML = TRUE)

# Here we are plotting partial residuals

visreg(M5_TMB, "HerbLayerc",
       by="fPlot",
       overlay=F,
       line=list(col=color_vis, lwd=2,lty=rep(1,n)),
       fill=list(col=add.alpha(color_vis,0.5)),
       points=list(cex=0.5),
       ylab="Shannon ",
       xlab="% Herbal Layer cover")

visreg(M5_TMB, "HerbLayerc",
       by="fPlot",
       overlay=T,
       line=list(col=color_vis, lwd=2,lty=rep(1,n)),
       fill=list(col=add.alpha(color_vis,0.5)),
       points=list(col=add.alpha(color_vis,0.5),cex=1),
       ylab="Shannon ",
       xlab="% Herbal Layer cover")


visreg(M5_TMB, "Litterc",
       by="fPlot",
       overlay=F,
       line=list(col=color_vis, lwd=2,lty=rep(1,n)),
       fill=list(col=add.alpha(color_vis,0.5)),
       points=list(cex=0.5),
       ylab="Shannon ",
       xlab="% Litter content")

visreg2d(M5_TMB, "HerbLayerc", "Litterc", plot.type="rgl")



###########################################################
# Example of random intercept and random slope model

M6 <- lmer(Hlog10 ~ HerbLayerc + (HerbLayerc| fPlot), 
           data = Spiders, REML = TRUE)

summary(M6)

# Checking M6
res_M6 <- simulateResiduals(fittedModel = M6, n = 500)
plot(res_M6)
plotResiduals(res_M6, Spiders$HerbLayerc)
plotResiduals(res_M6, Spiders$Litterc)
plotResiduals(res_M6, Spiders$GroundVegc)
r2(M6)

Spiders <- Spiders %>% mutate(fitted_M6=fitted(M6))

ggplot(Spiders ,aes(x=HerbLayer,y=Hlog10))+
  geom_point(show.legend = FALSE,alpha=0.5) +
  facet_wrap(~fPlot)+
  geom_abline(aes(intercept = coef(M0)[[1]], slope = coef(M0)[[2]],
                  color = 'black'),size = 1) + 
  geom_line(aes(x=HerbLayer,y=fitted_M1, color = "red"),size=1)+
  geom_line(aes(x=HerbLayer,y=fitted_M2, color = "blue"),size=1)+
  geom_line(aes(x=HerbLayer,y=fitted_M6, color = "green"),size=1)+
  scale_color_identity(breaks = c("black", "red", "blue","green"),
                       labels=c("pooled","non-pooled",
                                "rd intercept","rd intercept + slope"),
                       guide="legend")+
  labs(x = "Percentage of herb layer",y = "Shannon index",color="Model fit")

ggplot(Spiders %>% filter(fPlot==15) ,aes(x=HerbLayer,y=Hlog10))+
  geom_point(show.legend = FALSE,alpha=0.5) +
  facet_wrap(~fPlot)+
  geom_abline(aes(intercept = coef(M0)[[1]], slope = coef(M0)[[2]],
                  color = 'black'),size = 1) + 
  geom_line(aes(x=HerbLayer,y=fitted_M1, color = "red"),size=1)+
  geom_line(aes(x=HerbLayer,y=fitted_M2, color = "blue"),size=1)+
  geom_line(aes(x=HerbLayer,y=fitted_M6, color = "green"),size=1)+
  scale_color_identity(breaks = c("black", "red", "blue","green"),
                       labels=c("pooled","non-pooled",
                                "rd intercept","rd intercept + slope"),
                       guide="legend")+
  labs(x = "Percentage of herb layer",y = "Shannon index",color="Model fit")

r2(M2)
r2(M6)

######################################
# Example of GLMM

M7 <- glmmTMB(Richness ~ HerbLayer + (1| fPlot),
              family = poisson(),
              data = Spiders)

summary(M7)

# Checking M7
res_M7 <- simulateResiduals(fittedModel = M7, n = 500)
plot(res_M7)
testDispersion(res_M7)
plotResiduals(res_M7, Spiders$HerbLayerc)
plotResiduals(res_M7, Spiders$Litterc) # It's not OK
plotResiduals(res_M7, Spiders$GroundVegc)
r2(M7)

# plotting model

visreg(M7, "HerbLayer",
       by="fPlot",
       overlay=F,
       line=list(col=color_vis, lwd=2,lty=rep(1,n)),
       fill=list(col=add.alpha(color_vis,0.5)),
       points=list(cex=0.5),
       ylab="Richness",
       xlab="% Herbal Layer cover")



visreg(M7, "HerbLayer", scale="response",
       by="fPlot",
       overlay=F,
       line=list(col=color_vis, lwd=2,lty=rep(1,n)),
       fill=list(col=add.alpha(color_vis,0.5)),
       points=list(cex=0.5),
       ylab="Richness",
       xlab="% Herbal Layer cover")
