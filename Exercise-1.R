# Shrikant Pawar, class exercise 1, 08/21/2018.

# Load the datasets, packages
require(datasets)
require(survival)

# Check first few lines of the dataset lung
head(lung)

# Transfer the lung data to object survData 
survData <- lung

# Surv function to determine the censored data
survObj <- Surv(time=survData$time, event=survData$status==2, type='right')

# survfit function to intercept
fit <- survfit(survObj~1)

# plot the fit
plot(fit, main="K-M Plot for Lung Data", xlab="Time in Days", ylab="Proportion Surviving")

# survfit function to sexes- males and females
fit <- survfit(survObj~survData$sex==1)

# re-plot the fit
plot(fit, main="K-M Plot for Lung Data", xlab="Time in Days", ylab="Proportion Surviving", col=c(1:2))
legend('topright', c("Male","Female"), lty=1, col=c(1:2))

# Cox regression (or proportional hazards regression) is method for investigating the effect of several variables upon the time a specified event takes to happen.
# It is a  nonparametric because it does assume that the effects of the predictor variables upon survival are constant over time and are additive in one scale.
coxFit <- coxph(survObj ~ survData$sex==1)
summary(coxFit)





