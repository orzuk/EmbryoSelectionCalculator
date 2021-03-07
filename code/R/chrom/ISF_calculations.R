# Below stuff is for grant

# From Shai: 

risk_lowest = function(r2,K,n)
{
  r = sqrt(r2)
  zk = qnorm(K) # , lower.tail=F)
  integrand_lowest = function(t)
  {
    arg = (zk-t*sqrt(1-r2/2)) / (r/sqrt(2))
    y = dnorm(t)*pnorm(arg)^n # , lower.tail=F
    return(y)
  }
  risk = integrate(integrand_lowest,-Inf,Inf)$value
  #  reduction = (K-risk)/K
  return(risk)
}




# Compute conditional risk 
risk_lowest_conditional = function(r2,K,n,qf,qm) # ,relative=T)
{
  r = sqrt(r2)
  zk = qnorm(K)
  zqf = qnorm(qf)
  zqm = qnorm(qm)
  c = (zqf+zqm)/2 * r
  baseline = pnorm((zk-c)/sqrt(1-r^2/2)) # without selection
  integrand_lowest_cond = function(t)
  {
    arg = (zk-c-t*sqrt(1-r^2)) / (r/sqrt(2))
    y = dnorm(t)*pnorm(arg)^n
    return(y)
  }
  risk = integrate(integrand_lowest_cond,-Inf,Inf)$value
  #  if (relative) {
  #    reduction = (baseline-risk)/baseline
  #  } else {
  #    reduction = baseline-risk
  #  }
  return(list(baseline=baseline, risk=risk))
}





# Give disease risk as function of number of embryos
n.vec <- 1:20
max.n <- max(n.vec)
K <- 0.01 # disease prevalence
r2 <- 0.3 # PRS explaining

#q.mf <- c(0.2, 0.25, 0.5, 0.75, 0.8)
q.mf <- c(0.1, 0.5, 0.9)
num.q <- length(q.mf)
risk.vec <- rep(0, max.n)
cond.risk.vec <- matrix(0, nrow=num.q, ncol=max.n) # rep(0, )

for(n in n.vec)
{
  risk.vec[n] <- risk_lowest(r2,K,n)
  for(j in 1:num.q)  
    cond.risk.vec[j,n] <- risk_lowest_conditional(r2, K, n, q.mf[j], q.mf[j])$risk # ,relative=T)
}
col.vec <- c("red", "green", "blue", "orange", "purple", "pink")

jpeg('risk_conditional.jpg')
plot(n.vec, risk.vec, col=col.vec[1], ylim = c(0,max(cond.risk.vec)), xlab = "n", ylab = "risk", pch=20, cex=1.5, cex.lab=1.5, cex.axis=1.5)
lines(n.vec, risk.vec, col=col.vec[1])
for(j in 1:num.q)
{
  points(n.vec, cond.risk.vec[j,], col=col.vec[j+1], pch=20, cex=1.5)
  lines(n.vec, cond.risk.vec[j,], col=col.vec[j+1])
}  
#legend(14, 0.02, legend=c("baseline", as.character(q.mf)),
#       col=col.vec[1:(num.q+1)], pch=rep(20, 4), cex=0.8)
dev.off()

jpeg('relative_risk_conditional.jpg')
# Next plot relative risk reduction 
plot(n.vec[2:max.n], 1 - risk.vec[2:max.n]/risk.vec[1], ylim = c(0.4,1), xlim=c(2,20), col=col.vec[1], pch=20, xlab = "n", ylab = "relative risk reduction", 
     cex=1.5, cex.lab=1.5, cex.axis=1.5)
lines(n.vec[2:max.n], 1 - risk.vec[2:max.n]/risk.vec[1],  col=col.vec[1])
for(j in 1:num.q)
{
  points(n.vec[2:max.n], 1 - cond.risk.vec[j,2:max.n]/cond.risk.vec[j,1], col=col.vec[j+1], pch=20, cex=1.5)
  lines(n.vec[2:max.n], 1 - cond.risk.vec[j,2:max.n]/cond.risk.vec[j,1], col=col.vec[j+1])
}
legend(12, 0.65, legend=c("baseline", as.character(q.mf)),
       col=col.vec[1:(num.q+1)], pch=rep(20, 4), cex=1.5)
dev.off()


