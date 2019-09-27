# Script for calculating variance decomposition for the gain
setwd('C:\\Users\\Or Zuk\\Dropbox\\Embryo Selection Or Zuk\\Code\\R')
source('gain_moments.R')
source('sample_embryos.R')

# Compute integrals for variance of gain (need to define f3,f4,f5)
n  <- 5
Var1 <- (n-1) * integrate(gauss.poly, -8, 8, n=c(n-2,2,1))$value # var1 = (n-1) .* quad( @(t) pnorm(t).^(n-2) .* dnorm(t).^2.*t, -8, 8)
Var2 <- -1/n + integrate(gauss.poly, -8, 8, n=c(n-1,1,2))$value # var2 = -1/(n-0) + quad( @(t) pnorm(t).^(n-1) .* dnorm(t).*t.^2, -8, 8)
var_const <- integrate(gauss.poly, -5, 5, n=c(n-1,1,0))$value # var_const =  quad( @(t) pnorm(t).^(n-1) .* dnorm(t), -5, 5) % should be 1/(n-1) 

x_vec <- seq(-5,10, by=0.01)
max_x <- max(x_vec) 

#figure;  plot(x_vec, n .* pnorm(x_vec).^(n-1) .* dnorm(x_vec).*x_vec.^2)
#figure;  plot(x_vec, (n-1) .* pnorm(x_vec).^(n-2) .* dnorm(x_vec).^2.*x_vec)
#hold on; plot(x_vec, pnorm(x_vec).^(n-1) .* dnorm(x_vec).*x_vec.^2, 'r');
max_x = 10

#pnorm(max_x).^(n-1) .* dnorm(max_x) * max_x
#n=10; X = randn(1000000, n); mu = mean(X, 2); M = max(X, [], 2); EE = mean(X(:,1).*M)
#E2 =  -(n-1) .* integral( @(t) pnorm(t).^(n-2) .* dnorm(t).^2.*t, -8, 8) + (1/n)*mean(M.^2)
#E3 = 1/n % (1-(n-1)*mean(M.^2))/n
#n=5
#var_G_approx = pi*pi / (6 * (n * dnorm( qnorm(1/n)))^2 ) - 1/n

#%%%%%%%%%%%%%%%%%%%%%%%
# New: Determine within and between family variance
gain_data_dir = 'C:\\Users\\Or Zuk\\Dropbox\\Embryo Selection Or Zuk\\Theory\\figs';

max_n = 50 # number of embryos
N_fam = 2000 # number of families 
iters_fam = 1 # iterations to sample families 
iters = 100 # iterations to sample embryos 
M = 70 # number of blocks s
subtract_mean = 0 #  0 - 'max', 1 - 'max_minus_mean'
sib.flag = 0 # 0 - just take independent Gaussians

Total_E_G <- Var_E_G <- E_Var_G <- Var_G <- matrix(0, 1, max_n)
E_G2 <- E_G <- matrix(0, N_fam, max_n)

#Var_G_Approx  = 0.5 * ( pi*pi ./ (6 .* ((1:max_n) .* dnorm( qnorm(1./(1:max_n)))).^2 ) - 1./(1:max_n) ); % analytic
Var_G_approx <- E_G_approx <- Var_G_exact <- E_G_exact <- matrix(0, 1, max_n) 
for (n in c(1:max_n))
{
  G_exact = gain.moments(n, 1, 'exact', subtract_mean, sib.flag) # divide to mean and variance
  E_G_exact[n] <- G_exact$E.G
  Var_G_exact[n] <- G_exact$Var.G  
  G_approx = gain.moments(n, 1, 'approx', subtract_mean, sib.flag)
  E_G_approx[n] <- G_approx$E.G
  Var_G_approx[n] <- G_approx$Var.G  
}

for (i in c(1:iters_fam))  # loop to sample families 
{
  X.PP <- matrix( rnorm(N_fam*M, sd=1/(2*M)), N_fam, M) # randn(N_fam, M) / sqrt(2*M);
  X.PM <- matrix( rnorm(N_fam*M, sd=1/(2*M)), N_fam, M)
  X.MP <- matrix( rnorm(N_fam*M, sd=1/(2*M)), N_fam, M)
  X.MM <- matrix( rnorm(N_fam*M, sd=1/(2*M)), N_fam, M)
  #    E_G = zeros(N_fam, max_n); E_G2 = zeros(N_fam, max_n);
  
  run_i = i
  for (j in c(1:iters))
  {
    if(sib.flag)
    {
      E = sample.embryos(X.PP, X.PM, X.MP, X.MM, max_n) 
    } else 
      E = randn(N_fam, max_n) # multivariate Gaussian  # Sample embryos
    
    
    G <- apply(E, 2, cummax) - subtract_mean * apply(E, 2, cummax) / 
      matrix(rep(c(1:max_n), N_fam), ncol=N_fam, byrow = TRUE) # TEMP! Take only MAX!!!
    
    Var_G = Var_G + var(G) # Overall variance 
    E_G = E_G + G 
    E_G2 = E_G2 + G^2
  } # end for loop on j 
  
  E_G <- E_G / iters # take average gain for each family 
  Total_E_G <- Total_E_G + mean(E_G)
  E_G2 = E_G2 / iters # take average square gain for each family 
  E_Var_G = E_Var_G + mean(E_G2 - E_G^2) # within family     
  Var_E_G = Var_E_G+var(E_G) # between family 
}  # loop on families 

# Normalize 
Var_G = Var_G / (iters_fam * iters)
E_Var_G = E_Var_G / iters_fam 
Var_E_G = Var_E_G / iters_fam 
E_G = Total_E_G / iters_fam 


# Another figure     
title_str = paste0(' subtract-mean= ', subtract_mean, ' sibs=', sib.flag, ' M=', M)
#figure; hold on; xlabel('n'); ylabel('Mean'); title(['E(G) ' title_str]);
plot(c(1:max_n), E_G, title(title_str), xlab('n'), ylab('Mean')) # Plot mean  'linewidth', 2,
lines(c(1:max_n), E_G_exact, col='k') # , 'linewidth', 2
#legend({'E(G)', 'E(G) Analytic'}, 'location', 'southeast'); legend('boxoff'); % , 'E(G) Approx'
lines(c(1:max_n), E_G_approx, col='green') # , 'linewidth', 2);
#legend({'E(G)', 'E(G) Analytic', 'E(G) Approx'}, 'location', 'southeast'); legend('boxoff'); % 


# Plot differnet variance components 
# figure; hold on; xlabel('n'); ylabel('Var'); title(['Var(G) ' title_str]);
plot(c(1:max_n), Var_G) #  'linewidth', 2); 
lines(c(1:max_n), Var_E_G, col='red')
lines(c(1:max_n), E_Var_G, col='green')
lines(c(1:max_n), E_Var_G+Var_E_G, col='magenta') # , 'linewidth', 1);
lines(c(1:max_n), Var_G_exact, col='black') # , 'linewidth', 2);
#legend({'Var(G)', 'Var[E(G)] (between)', 'E[Var(G)] (within)', 'between+within', 'Var(G) Analytic'}, 'location', 'southeast'); legend('boxoff'); % , 'Var(G) Approx'
lines(c(1:max_n), Var_G_approx, col='cyan') # , 'linewidth', 2); % pretty lousy approximation
#legend({'Var(G)', 'Var[E(G)] (between)', 'E[Var(G)] (within)', 'between+within', 'Var(G) Analytic', 'Var(G) Approx'}, 'location', 'northeast'); legend('boxoff'); % 
#my_saveas(gcf, fullfile(gain_data_dir, 'gain_var'), {'epsc', 'pdf', 'jpg'}); 



#%%%%%%%%%%%%%%%%%%
# Plot covariance histogram for one realization 
C = matrix(0, max_n, max_n)
for (i in c(1:N_fam))
  C = C + t(E[i,]) * E[i,]
C = C / N_fam; 
# How to do imagesc in R?
#figure; imagesc(C); colorbar; title(['Empirical Covariance for family, Var=' num2str(median(diag(C))) ', COV=' num2str(median(C(:)))]); 

I = 24; C = t(E[I,]) * E[I,] 
hist(C, 100) #  xlabel('X_i X_j'); ylabel('Freq.'); median(C(:))


