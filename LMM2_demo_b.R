# aids data set

#library("JM")
#write.csv(aids, file="aids.csv", row.names = FALSE)

aids <- read.csv("Downloads/aids.csv", header = TRUE)
aids[aids$patient %in% c(1,2), c("patient", "CD4", "obstime")]




#### 2.4 Fitting random-intercepts and random-slopes model with REML
library(nlme)
# random-intercepts and random-slopes model for the AIDS dataset
lmeFit <- lme(CD4 ~ obstime, random = ~ obstime | patient, data = aids, method = "REML")

summary(lmeFit)

#The fixed effects coefficients: 
beta_hat=lmeFit$coefficients$fixed
beta_hat


#The variance parameter $\hat{\sigma}$. 
sigma_hat=lmeFit$sigma
sigma_hat


#The covariance matrix D.
D <- getVarCov(lmeFit)
D <- as.matrix(D)
dim(D)


#The random effects coefficients:
#The following code print $\hat{\bf b}_1,\ldots,\hat{\bf b}_{10}$.

lmeFit$coefficients$random$patient[1:10,]




#### 2.6 Subject specific covariance structure and prediction 
aids12 = aids[aids$patient==12,]
aids12

Xi <- Zi <- cbind(Int = 1, obstime = aids12$obstime)
Xi <- as.matrix(Xi) 
Zi <- as.matrix(Zi)

yi <- aids12$CD4

Sigma_i <- Zi%*%D%*%t(Zi)+(sigma_hat^2)*diag(1,5)
Sigma_i

### b_hat 
b_hat_i = D%*%t(Zi)%*%solve(Sigma_i)%*%(yi- Xi %*% beta_hat)
b_hat_i

lmeFit$coefficients$random$patient[12,]



### y_hat population
y_hat <- Xi%*%beta_hat 
y_hat


aids12$predict <- predict(lmeFit, aids12, level=0)
aids12[, c("patient","CD4","predict")]


### y_hat individual
beta_hat + b_hat_i 

y_hat_i <- Xi%*%beta_hat + Zi%*%b_hat_i 
y_hat_i


aids12$predict <- predict(lmeFit, aids12, level=1)
aids12[, c("patient","CD4","predict")]



































