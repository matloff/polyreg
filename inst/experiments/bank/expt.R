
# change 'unknown' to NA
f <- function()                                                                 {   
  for(i in 1:ncol(bank)) {
     b <- bank[,i]
     if (is.factor(b)) {
        unk <- which(b == 'unknown')
        bank[unk,i] <- NA
     }  
  } 
  bank  
}   
bank <- 
   read.table('~/Research/DataSets/Bank/bank-full.csv',header=T,sep=';')
bank$y <- as.integer(bank$y == 'yes')
bank$job <- NULL  # had data issues
bank$poutcome <- NULL  # almost all unknown
bnk <- f()  # see f below; changes 'unknown' to NA
bnk1 <- factorsToDummies(bnk)
doGenExpt(bnk1,regftn=glm)

