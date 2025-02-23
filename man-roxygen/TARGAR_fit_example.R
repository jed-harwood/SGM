\dontrun{
  ### Load data
  load("data/stocks.Rda")
  n=nrow(stocks)
  p=ncol(stocks)
  model="LN"
  
  ### lambda and net.thre sequence
  C.v=c(8,4,1)
  lambda.v=C.v*sqrt(log(p)/(n-1)) 
  rho.v=pmax(lambda.v, 0.01)
  
  C.thre=exp(seq(log(1),log(0.1), length.out=12))
  net.thre=C.thre*sqrt(log(p)/(n-1)) 
  
  ### Run GAR1_fit
  resList = TARGAR_fit(stocks, lambda.v, net.thre, rho.v, num.pass = 3)

}
