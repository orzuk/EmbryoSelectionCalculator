# Sample the Poygenic Scores of set of embryos
# Input: 
# X.PP - polygenic scores of parental grandfather
# X.PM - polygenic scores of parental grandmother
# X.MP - polygenic scores of marental grandfather
# X.MM - polygenic scores of marental grandmother
# n - number of embryos
#
sample.embryos <- function(X.PP, X.PM, X.MP, X.MM, n)
{
  n.fam <- dim(X.PP)[1]
  M <- dim(X.PP)[2]  
  E = matrix(0, n.fam, n)
  for (i in c(1:n))
  {
    C.P <- matrix(rbinom(n.fam*M, 1,  0.5), nrow=n.fam)  
    C.M <- matrix(rbinom(n.fam*M, 1, 0.5), nrow=n.fam)  
    E[,i] = rowSums(X.PP*C.P + X.PM*(1-C.P) + X.MP*C.M + X.MM*(1-C.M)) 
  }
  return(E)
}  
