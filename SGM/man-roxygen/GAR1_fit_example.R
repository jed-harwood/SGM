\dontrun{
### Load data
load("data/stocks.Rda")
n=nrow(stocks)
d=ncol(stocks)
model="LN"

### lambda and net.thre sequence
C.v=c(8,4,1)
lambda.v=C.v*sqrt(log(d)/n)
rho.v=pmax(lambda.v, 0.01)

C.thre=exp(seq(log(1),log(0.1), length.out=12))
net.thre=C.thre*sqrt(log(d)/n)

### Run GAR1_fit
S = var(stocks)*(n-1)/n
resList = GAR1_fit(S, n, lambda.v, net.thre, model, verbose = T,
                  max_iter_s1 = 100000, max_iter_s2 = 100000, max_iter_s3 = 100000) 
}
