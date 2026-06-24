rm(list=ls())
os <- Sys.info()
if( os[1] == "Darwin" ){ # local machine
  mainDir="./real/"
  funcPath <- paste0("./Functions/")
}else{ # Univ cluster
  mainDir="/users/hchung13/mocl/real/"
  funcPath <- paste0("/users/hchung13/mocl/Functions/")
  .libPaths("/users/hchung13/R/x86_64-pc-linux-gnu-library/4.3/")
}
resultPath <- paste0(mainDir,"results/")
source(paste0(mainDir,"real_header/real_header.R"))

#######################################
mlb2021 <- read.csv( paste0(mainDir,"data/stats2021.csv") )
head(mlb2021)
mlb2021 <-  mlb2021[,-c(1,2,3,4)]
head(mlb2021)
colnames(mlb2021) <- c("strikeout","walk","batting","power","woba",
                       "velocity","angle","well-hit","hard","contact")

colnames(mlb2021)[5]
mlb_x <- scale(mlb2021[,-5])
# mlb_x <- mlb2021[,-5]
y <- mlb2021[,5]
#######################################


(file <- paste0(mainDir,"data/","baseball.rds"))
saveRDS(list(y = y, x = mlb_x), file = file)
















