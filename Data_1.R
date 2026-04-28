rm(list=ls())
#install.packages("readxl")
library("readxl")           

house_data <- read_excel("All_data.xls", sheet = "Return")  
stock_data <- read.csv("stock_vol.csv", header=TRUE)

stock_vol    <- log(as.numeric(stock_data[4:555, 2]))
stock_vol[1] <- stock_vol[2] 

EI_DATA <- cbind(house_data, stock_vol)

write.csv(EI_DATA, file="EI_DATA.csv", row.names=FALSE)


