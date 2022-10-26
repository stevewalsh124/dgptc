library(gtools) #mixedsort()

# generate a design for SA main effects
p = 8

rm_avg <- T

# generate parameter settings that are uniform over [0,1]^8
# udesign = xx from the xstar.R script
xstar <- read.csv("/home/walsh124/CosmicEmu/2022-Mira-Titan-IV/P_tot/xstar.dat", sep="", header = F)

# generate the parameter variations for the SA
nv = 11
# prange = seq(0,1,length=nv) # already done in xstar.R

# Check to see if plot will be conditional or average effect
# That is, points fixed at 0.5 or LHS design
if(length(unique(xstar[,1]))==nv) condl <- T else condl <- F

# load k values
load("Mira-Titan-IV-data/precision_and_indexes.Rdata")
kvals = log10(k)

# make some emulator runs - each column is a model run
EMUfiles <- list.files("/home/walsh124/CosmicEmu/2022-Mira-Titan-IV/P_tot/", full.names = T, recursive = F, pattern = ".txt")
EMUfiles <- mixedsort(EMUfiles)
EMU_outs <- do.call(cbind, lapply(EMUfiles, function(x){read.csv(x, header = F, sep="")[,2]}))
modRuns = log10(EMU_outs)
m = ncol(modRuns)/(8*nv)  # number of MC runs to estimate mean effects

# do this for all 8 parameters
mainEffs = array(NA,c(length(kvals),nv,p))
for(k in 1:8){
  # collect sumulations
  sims1 = modRuns[,(k-1)*(m*nv)+1:(m*nv)]
  # reshape into a 3d array by kval,rep, paramval
  sims1a = array(sims1,c(length(kvals),m,nv))
  # compute main effect
  mainEffs[,,k] = apply(sims1a,c(1,3),mean)
}

if(rm_avg) avg_mean <- rowMeans(modRuns) else avg_mean <- rep(0, nrow(modRuns))


# make a plot - you could do this with ggplot if you'd rather
pnames = c('omega_m','omega_b','sigma_8','h','n_s','w_0','w_a','omega_nu','z')
PDF = TRUE  # switch to TRUE if you want R to make a pdf
if(PDF) pdf(paste0("pdf/mainEffs",if(rm_avg){"_rm_avg"},"_",m,if(condl){"_condl"},".pdf"),width=8,height=4)
par(mfrow=c(2,4),oma=c(4,4,1,1),mar=c(0,0,0,0))
yr = range(mainEffs - avg_mean)
grcolors = paste0("grey",round(seq(90,25,length=11)),sep='')
for(k in 1:8){
  matplot(kvals,mainEffs[,,k] - avg_mean,type='l',col=grcolors,lty=1,ylim=yr,axes=F); box()
  text(-3,yr[1]+0.01,pnames[k],adj=c(0,0))
  if(k %in% c(1,5)) axis(2)
  if(k %in% c(5,7)) axis(1)
}
mtext('k',side=1,line=2.5,outer=T)
mtext(paste0('log spectrum', if(rm_avg){" (average removed)"}),side=2,line=2.5,outer=T)
if(PDF) dev.off()

