\dontrun{
### Load data
load("data/stocks.Rdata")
n=nrow(stocks)
p=ncol(stocks)
model="LN"

### lambda and net.thre sequence
C.v=c(8,4,1)
lambda.v=C.v*sqrt(log(p)/n) 
rho.v=pmax(lambda.v, 0.01)

C.thre=exp(seq(log(1),log(0.1), length.out=12))
net.thre=C.thre*sqrt(log(p)/n) 

### Run GAR1_fit
S = var(stocks)*(n-1)/n
resList = GAR1_fit(S, n, lambda.v, net.thre, model) 
}
