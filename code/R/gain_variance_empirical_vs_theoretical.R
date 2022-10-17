# Create Supp. Fig. 1 in the paper - variance of G
setwd('C:\\Code\\GitHub\\EmbryoSelectionCalculator\\code\\R')  # change to your path
source('gain_moments.R')
data = read.csv("../../Data/var_gain.txt",sep="\t")
ns = data$n_sibs # number of siblings
nm = max(ns)
ns = ns[2:nm]
var.emp = data$gain_var[2:nm] # empirical variance (simulations from longevity cohort)

sigma.z = 6 # variance of trait
r2ps = 0.243 # proportion of variance explained by PGS

var.th <- gain_moments(ns, 1, 'approx', 1, 1)$Var.G * sigma.z^2 * r2ps # asymptotic approximation
var.num <- gain_moments(ns, 1, 'exact', 1, 1)$Var.G * sigma.z^2 * r2ps # numeric integral

# png('../../Figures/GainVar.png', res=300) # save file 
png('../../Figures/GainVar.png', width = 350, height = 300, units='mm', res = 300)
par(mar=c(5,6,4,1)+.1) # adjust margins
plot(ns,var.emp,ylim=c(0,3), xlab="n", ylab="Var(G)", cex=3, cex.lab=2.5, cex.axis=2.5, lwd=4)
lines(ns, var.th, lwd=4, col="blue")
lines(ns, var.num, lwd=4, col="red")
dev.off()


