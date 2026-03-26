####################################################
## Data Pre-Processing for Stock Data Application ##
## Outputs stock.data.X.Rdata, and EDA plots from ##
## -data Obtained via Loggle R Package            ##
####################################################

### Load the Data (taken from loggle package, used in Yang and Peng '19)
rm(list=ls())
load("stockdata.rda")
date.index <- rownames(stockdata$data) 

start.date='2007-01-01'  ## great recession period 
end.date='2010-12-31'
date.index.sele=(date.index>=start.date)&(date.index<=end.date) 

stock.sector <- stockdata$info[, "sector"]
sp.it <- t(stockdata$data[date.index.sele, stock.sector == "Information Technology"])
sp.cd <- t(stockdata$data[date.index.sele, stock.sector == "Consumer Discretionary"])
sp.cs <- t(stockdata$data[date.index.sele, stock.sector == "Consumer Staples"])
sp.f <- t(stockdata$data[date.index.sele, stock.sector == "Financials"])
sp.i <- t(stockdata$data[date.index.sele, stock.sector == "Industrials"])
sp <- rbind(sp.it, sp.cd, sp.cs, sp.f, sp.i)
sp.num <- c(nrow(sp.it), nrow(sp.cd), nrow(sp.cs), nrow(sp.f), nrow(sp.i))
names(sp.num)=c("I. T.", "Cons. Disc.", "Cons. Staples", "Financials", "Industrials")
sector.num <- length(sp.num)
dim(sp) ##283 by 1008

## construct data matrix by taking log ratio of prices between adjacent time points
p <- dim(sp)[1]
n <- dim(sp)[2]-1
X <- matrix(0, n, p)
for(i in 1:p) {
  X[, i] <- scale(log(sp[i, -1] / sp[i, -(n+1)])) ##standardized 
}
dim(X)  ## dimension of data matrix: 283x 1007
colnames(X)=rownames(sp)
rownames(X)=colnames(sp)[-1]

## Save the data as a Matrix
save.image(file="stock.data.X.Rdata")


## examine and plot data matrix X
par(mfrow=c(2,3))
matplot(1:1007, X[,rownames(sp.it)], type='l', xlab="day", ylab="log return", main="Information Technology")
matplot(1:1007, X[,rownames(sp.cd)], type='l', xlab="day", ylab="log return", main="Consumer Discretionary")
matplot(1:1007, X[,rownames(sp.cs)], type='l', xlab="day", ylab="log return", main="Consumer Staples")
matplot(1:1007, X[,rownames(sp.f)], type='l', xlab="day", ylab="log return", main="Financials")
matplot(1:1007, X[,rownames(sp.i)], type='l', xlab="day", ylab="log return", main="Industrials")
matplot(1:1007, X, type='l', xlab="day", ylab="log return", main="All Sectors")
par(mfrow=c(1,1))



## examine autocorrelation
library(sarima)
maxlag=1006
autocorr=matrix(NA, ncol(X), maxlag) ## auto-correlation
pval.garch= matrix(NA, ncol(X), 4) ## pvalue of the whitenoise test based on autocorrelation; h0: garch
pval.iid= matrix(NA, ncol(X), 4) ## pvalue of the whitenoise test based on autocorrelation; h0:iid


for (i in 1: ncol(X)){
  auto=autocorrelations(X[,i], maxlag=maxlag)
  autocorr[i,]=(as.numeric(auto@data))[-1]
  
  test=whiteNoiseTest(auto, h0="garch", nlags = c(1,5,30, 365), x=X[,i])
  pval.garch[i,]=test$test[,"pval"]
  
  test=whiteNoiseTest(auto, h0="iid", nlags = c(1,5,30, 365), x=X[,i])
  pval.iid[i,]=test$test[,"pvalue"]
}

# autocorrelation decays with lag, max < 0.2; 
summary(pval.iid)
summary(pval.garch)
