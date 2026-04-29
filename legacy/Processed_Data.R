rm(list=ls())
#install.packages("readxl")
library("readxl")           

data <- read_excel("All_data.xls", sheet = "Monthly")  
head(data)

monthly_FFR         <- ts(data$FFR, start = c(1975, 1), frequency = 12)
monthly_shadow_rate <- ts(data$shadow_rate, start = c(1975, 1), frequency = 12)
 
#Quarterly data
FFR         <- aggregate(round(monthly_FFR, 2), nfrequency = 4, FUN = mean)
SHADOW_RATE <- aggregate(monthly_shadow_rate,  nfrequency = 4, FUN = mean)
SHADOW_RATE[1:136] <- FFR[1:136] 


dates       <- time(SHADOW_RATE)
SHADOW_RATE <- data.frame(Date = dates, SHADOW_RATE)
SHADOW_RATE 
write.csv(SHADOW_RATE[1:184,], file="SHADOW_RATE.csv", row.names=FALSE)

