# explore cosmology data

# Comparing the values between the different P(k) at different steps
#####
# aa <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-no_interp/STEP",163,"/pk_M",
#                         if(i<100){"0"},if(i<10){"0"},i,"_no_interp_test.dat"))
# bb <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-no_interp/STEP",189,"/pk_M",
#                         if(i<100){"0"},if(i<10){"0"},i,"_no_interp_test.dat"))
# cc <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-no_interp/STEP",247,"/pk_M",
#                         if(i<100){"0"},if(i<10){"0"},i,"_no_interp_test.dat"))
# dd <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-no_interp/STEP",300,"/pk_M",
#                         if(i<100){"0"},if(i<10){"0"},i,"_no_interp_test.dat"))
# ee <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-no_interp/STEP",347,"/pk_M",
#                         if(i<100){"0"},if(i<10){"0"},i,"_no_interp_test.dat"))
# ff <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-no_interp/STEP",401,"/pk_M",
#                         if(i<100){"0"},if(i<10){"0"},i,"_no_interp_test.dat"))
# gg <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-no_interp/STEP",453,"/pk_M",
#                         if(i<100){"0"},if(i<10){"0"},i,"_no_interp_test.dat"))
# hh <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-no_interp/STEP",499,"/pk_M",
#                         if(i<100){"0"},if(i<10){"0"},i,"_no_interp_test.dat"))
# all.equal(aa,bb)
# all.equal(aa,cc)
# all.equal(aa,dd)
# all.equal(aa,ee)
# all.equal(aa,ff)
# all.equal(aa,gg)
# all.equal(aa,hh)
#####

steps <- 499 #c(163, 189, 247, 300, 347, 401, 453, 499)
z <- 0 #c(2.02, 1.61, 1.006, 0.656, 0.434, 0.242, 0.101, 0)

##########################
# non-interp (3000 rows) #
##########################

for(step in steps){
  
  print(step)
  
  pdf(paste0("pdf/EDA_k_",step,".pdf"))
  for (i in 0:111) {
    k <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-no_interp/STEP",step,"/k_M",
                            if(i<100){"0"},if(i<10){"0"},i,"_no_interp_test.dat"))
    plot(k[(k[,18]!=0),18], type="l", log="xy", col="red", ylim = range(k[k!=0]), ylab="k",
         main = paste("k, model",i))
    for(j in 2:17) lines(k[,j], col="blue")
    for(j in 1) lines(k[,j], col="green")
    legend("topleft",c("lo","hi","pert"), col=c("blue","red","green"), lty=1)
  }
  dev.off()
  
  pdf(paste0("pdf/EDA_err_",step,".pdf"))
  for (i in 0:111) {
    k <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-no_interp/STEP",step,"/k_M",
                           if(i<100){"0"},if(i<10){"0"},i,"_no_interp_test.dat"))
    if(ncol(k)==18) k <- k[,-1]
    err <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-no_interp/STEP",step,"/err_M",
                            if(i<100){"0"},if(i<10){"0"},i,"_no_interp_test.dat"))
    plot(k[(k[,17]!=0)&(err[,17]!=0),17],err[(k[,17]!=0)&(err[,17]!=0),17], type="l",  
         ylim=range(err[err!=0], na.rm = T), log="xy", col="red", xlab="k", ylab="err",
         main = paste("err, model",i))
    for(j in 1:16) lines(k[(k[,j]!=0)&(err[,j]!=0),j],err[(k[,j]!=0)&(err[,j]!=0),j], col="blue")
    legend("bottomleft",c("lo","hi"), col=c("blue","red"), lty=1)
    
  }
  dev.off()
  
  pdf(paste0("pdf/EDA_pk_",step,".pdf"))
  for (i in 0:111) {
    k <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-no_interp/STEP",step,"/k_M",
                           if(i<100){"0"},if(i<10){"0"},i,"_no_interp_test.dat"))
    pk <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-no_interp/STEP",step,"/pk_M",
                            if(i<100){"0"},if(i<10){"0"},i,"_no_interp_test.dat"))
    plot(k[(k[,1]!=0)&(pk[,1]!=0),1],pk[(k[,1]!=0)&(pk[,1]!=0),1], type="l",  
         ylim=range(pk[pk!=0]), xlim=range(k[k!=0]), log="xy", col="green",
         xlab="k", ylab="P(k)", main = paste("P(k), model",i))
    for(j in 2:17) lines(k[(k[,j]!=0)&(pk[,j]!=0),j],pk[(k[,j]!=0)&(pk[,j]!=0),j], col="blue")
    for(j in 18) lines(k[(k[,j]!=0)&(pk[,j]!=0),j],pk[(k[,j]!=0)&(pk[,j]!=0),j], col="red", lwd=2)
    legend("bottomleft",c("lo","hi","pert"), col=c("blue","red","green"), lty=1)
    
  }
  dev.off()
  
}

###########################
# 2021 (interp; 351 rows) #
###########################

for(step in steps){
  
  print(step)
  
  # "k"s are the same for both 2021 and no_interp
  # pdf(paste0("pdf/EDA_k_",step,"_2021.pdf"))
  # for (i in 0:111) {
  #   k <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-2021/STEP",step,"/k_M",
  #                           if(i<100){"0"},if(i<10){"0"},i,"_test.dat"))
  #   plot(k[,1], type="l", ylim=c(0, max(k)))
  #   for(j in 2:17) lines(k[,j], col=j)
  # }
  # dev.off()
  
  pdf(paste0("pdf/EDA_err_",step,"_2021.pdf"))
  for (i in 0:111) {
    err2 <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-2021/STEP",step,"/err_M",
                            if(i<100){"0"},if(i<10){"0"},i,"_test.dat"))
    plot(err2[,1:2], type="l",  log="xy", xlab="k", ylab = "err", col="blue",
         ylim = range(err2[,-1][err2[,-1]!=0]), main = paste("err, model",i))
    for(j in 2:17) lines(err2[,c(1,j)], col="blue")
    for(j in 18) lines(err2[,c(1,j)], col="red")
    legend("bottomleft",c("lo","hi"), col=c("blue","red"), lty=1)
  }
  dev.off()
  
  pdf(paste0("pdf/EDA_pk_",step,"_2021.pdf"))
  for (i in 0:111) {
    pk2 <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-2021/STEP",step,"/pk_M",
                              if(i<100){"0"},if(i<10){"0"},i,"_test.dat"))
    plot(pk2[,1:2], type="l",  log="xy", xlab="k", ylab = "P(k)", col="green",
         ylim = range(pk2[,-1][pk2[,-1]!=0]), main = paste("P(k), model",i))
    for(j in 3:18) lines(pk2[,c(1,j)], col="blue")
    for(j in 19) lines(pk2[,c(1,j)], col="red")
    legend("bottomleft",c("lo","hi","pert"), col=c("blue","red","green"), lty=1)
  }
  dev.off()
  
}

####################
# EDA on prec data #
####################

pdf("pdf/prec_EDA.pdf")
step <- 499
i <- 1 # Model 1, choose from 000-111
pk2 <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-2021/STEP",step,"/pk_M",
                         if(i<100){"0"},if(i<10){"0"},i,"_test.dat"))

# load the precision data (k, prec_highres, prec_lowres, index_list)
load("Mira-Titan-IV-data/precision_and_indexes.Rdata")
plot(k, prec_lowres, type="l", main = "precision")
plot(k, 1/prec_lowres, type="l", main = "Var = 1/precision")
plot(k, sqrt(1/prec_lowres), type="l", main = "SD = sqrt(1/precision)")
plot(k, sqrt(1/prec_lowres), type="l", log="xy", main = "SD = sqrt(1/precision), log scale")
plot(k[index_list$lowres.ix], sqrt(1/prec_lowres[index_list$lowres.ix]), 
     type="l", log="xy", main = "SD = sqrt(1/precision), log scale, only lowres_index")

# Showing what the data originally look like
plot(pk2[,1:2], type="l",  log="xy", xlab="k", ylab = "P(k)", col="green",
     ylim = range(pk2[,-1][pk2[,-1]!=0]), main = paste("P(k), model",i))
for(j in 3:18) lines(pk2[,c(1,j)], col="blue")
for(j in 19) lines(pk2[,c(1,j)], col="red")
legend("bottomleft",c("lo","hi","pert"), col=c("blue","red","green"), lty=1)

# Show which low-res values have precisions
ind.plot <- ifelse(1:length(k) %in% index_list$lowres.ix, 1000, NA)
hi.ind.plot <- ifelse(1:length(k) %in% index_list$highres.ix, 1000, NA)
pt.ind.plot <- ifelse(1:length(k) %in% index_list$highres.ix, 1000, NA)

# Showing what the data originally look like
plot(cbind(pk2[,1], ifelse(pk2[,2] == 1, NA, pk2[,2])), type="l",  
     log="xy", xlab="k", ylab = "P(k)", col="green",
     ylim = range(pk2[,3:18][pk2[,3:18]!=0]), main = paste("P(k), model",i,"; line shows lowres_index"))
for(j in 3:18) lines(pk2[,c(1,j)], col="blue")
for(j in 19) lines(pk2[,c(1,j)], col="red")
lines(k,ind.plot)
lines(k,500+hi.ind.plot, lty=2)
lines(k,-500+pt.ind.plot, lty=2)
legend("bottomleft",c("lo","hi","pert","lowres index", "highres/pt index"), col=c("blue","red","green","black","black"), lty=c(1,1,1,1,2))

# all of the data, none removed
plot(pk2[,c(1,3)], col="blue", log="xy", type="l", main="seems too precise...")
legend("bottomleft",c("one lowres run"), col=c("blue"), lty=1)

plot(pk2[,c(1,3)], col="blue", log="xy", type="l", main="seems too precise...")
lines(pk2[,c(1,3)] - cbind(0, 2*1/sqrt(prec_lowres)))
lines(pk2[,c(1,3)] + cbind(0, 2*1/sqrt(prec_lowres)))
legend("bottomleft",c("one lowres run", "LB/UB"), col=c("blue","black"), lty=1)

# all of the nonzero data for low res
plot(pk2[pk2[,j]!=0,c(1,j)], col="blue", log="xy", type="l", main = "within vs between-run variation")
lines(pk2[pk2[,j]!=0,c(1,3)] - cbind(0, 2*sqrt(1/prec_lowres[pk2[,j]!=0])))
lines(pk2[pk2[,j]!=0,c(1,3)] + cbind(0, 2*sqrt(1/prec_lowres[pk2[,j]!=0])))


# the low res 
plot(pk2[index_list$lowres.ix,c(1,3)], col="blue", log="xy", type="l", 
     ylim = range(c(pk2[index_list$lowres.ix,c(3)] + 3*sqrt(1/prec_lowres[index_list$lowres.ix]),
                    pk2[index_list$lowres.ix,c(3)] - 3*sqrt(1/prec_lowres[index_list$lowres.ix]))),
     main = "all 18 low res, one set of 95% CIs from prec")
for(j in 3:18) lines(pk2[index_list$lowres.ix,c(1,j)], col="blue")
lines(cbind(pk2[index_list$lowres.ix,1],
            pk2[index_list$lowres.ix,j] - 2*sqrt(1/prec_lowres[index_list$lowres.ix])))
lines(cbind(pk2[index_list$lowres.ix,1],
            pk2[index_list$lowres.ix,j] + 2*sqrt(1/prec_lowres[index_list$lowres.ix])))
legend("bottomleft", c("lowres","CIs"), col=c("blue","black"), lty=1)

plot(pk2[index_list$lowres.ix,c(1,j)], col="blue", type="l", main = "orig scale, only chg P(k)")
for (j in 3:18) lines(pk2[index_list$lowres.ix,c(1,j)], col="blue")
lines(pk2[index_list$lowres.ix,c(1,j)] + cbind(0, 2*sqrt(1/prec_lowres[index_list$lowres.ix])))
lines(pk2[index_list$lowres.ix,c(1,j)] - cbind(0, 2*sqrt(1/prec_lowres[index_list$lowres.ix])))

plot(pk2[index_list$lowres.ix,c(1,j)], col="blue", type="l", log="xy", main = "log scale, only chg P(k)")
for (j in 3:18) lines(pk2[index_list$lowres.ix,c(1,j)], col="blue")
lines(pk2[index_list$lowres.ix,c(1,j)] + cbind(0, 2*sqrt(1/prec_lowres[index_list$lowres.ix])))
lines(pk2[index_list$lowres.ix,c(1,j)] - cbind(0, 2*sqrt(1/prec_lowres[index_list$lowres.ix])))

plot(pk2[index_list$lowres.ix,c(1,j)], col="blue", type="l",  main = "orig, only chg k & P(k) WRONG")
for (j in 3:18) lines(pk2[index_list$lowres.ix,c(1,j)], col="blue")
lines(pk2[index_list$lowres.ix,c(1,j)] + 2*sqrt(1/prec_lowres[index_list$lowres.ix]))
lines(pk2[index_list$lowres.ix,c(1,j)] - 2*sqrt(1/prec_lowres[index_list$lowres.ix]))

plot(pk2[index_list$lowres.ix,c(1,j)], col="blue", type="l", log="xy",  main = "log, chg k & P(k) WRONG")
for (j in 3:18) lines(pk2[index_list$lowres.ix,c(1,j)], col="blue")
lines(pk2[index_list$lowres.ix,c(1,j)] + 2*sqrt(1/prec_lowres[index_list$lowres.ix]))
lines(pk2[index_list$lowres.ix,c(1,j)] - 2*sqrt(1/prec_lowres[index_list$lowres.ix]))

dev.off()

# Boxplots
#####

# pdf("pdf/EDA_boxplots.pdf")
# i <- 0
# k <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-no_interp/STEP",step,"/k_M",
#                        if(i<100){"0"},if(i<10){"0"},i,"_no_interp_test.dat"))
# err <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-no_interp/STEP",step,"/err_M",
#                          if(i<100){"0"},if(i<10){"0"},i,"_no_interp_test.dat"))
# pk <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-no_interp/STEP",step,"/pk_M",
#                         if(i<100){"0"},if(i<10){"0"},i,"_no_interp_test.dat"))
# err2 <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-2021/STEP",step,"/err_M",
#                           if(i<100){"0"},if(i<10){"0"},i,"_test.dat"))
# pk2 <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-2021/STEP",step,"/pk_M",
#                          if(i<100){"0"},if(i<10){"0"},i,"_test.dat"))
# boxplot(k, main = paste0("k (no_interp), i=",i))
# boxplot(pk, main = paste0("pk (no_interp), i=",i))
# boxplot(pk2, main = paste0("pk (2021), i=",i))
# boxplot(err, main = paste0("err (no_interp), i=",i))
# boxplot(err2, main = paste0("err (2021), i=",i))
# 
# i <- 111
# k <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-no_interp/STEP",step,"/k_M",
#                        if(i<100){"0"},if(i<10){"0"},i,"_no_interp_test.dat"))
# err <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-no_interp/STEP",step,"/err_M",
#                          if(i<100){"0"},if(i<10){"0"},i,"_no_interp_test.dat"))
# pk <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-no_interp/STEP",step,"/pk_M",
#                         if(i<100){"0"},if(i<10){"0"},i,"_no_interp_test.dat"))
# err2 <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-2021/STEP",step,"/err_M",
#                           if(i<100){"0"},if(i<10){"0"},i,"_test.dat"))
# pk2 <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-2021/STEP",step,"/pk_M",
#                          if(i<100){"0"},if(i<10){"0"},i,"_test.dat"))
# boxplot(k, main = paste0("k (no_interp), i=",i))
# boxplot(pk, main = paste0("pk (no_interp), i=",i))
# boxplot(pk2, main = paste0("pk (2021), i=",i))
# boxplot(err, main = paste0("err (no_interp), i=",i))
# boxplot(err2, main = paste0("err (2021), i=",i))
# 
# dev.off()
# 
# pdf("pdf/box_k_no_interp.pdf")
# for(i in 0:111){
#   k <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-no_interp/STEP",step,"/k_M",
#                          if(i<100){"0"},if(i<10){"0"},i,"_no_interp_test.dat"))
#   boxplot(k, main = paste0("k (no_interp), i=",i))
# }
# dev.off()
# 
# pdf("pdf/box_pk_no_interp.pdf")
# for(i in 0:111){
#   pk <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-no_interp/STEP",step,"/pk_M",
#                           if(i<100){"0"},if(i<10){"0"},i,"_no_interp_test.dat"))
#   boxplot(pk, main = paste0("pk (no_interp), i=",i))
# }
# dev.off()
# 
# pdf("pdf/box_err_no_interp.pdf")
# for(i in 0:111){
#   err <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-no_interp/STEP",step,"/err_M",
#                            if(i<100){"0"},if(i<10){"0"},i,"_no_interp_test.dat"))
#   boxplot(err, main = paste0("err (no_interp), i=",i))
# }
# dev.off()
# 
# pdf("pdf/box_pk_2021.pdf")
# for(i in 0:111){
#   pk <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-2021/STEP",step,"/pk_M",
#                           if(i<100){"0"},if(i<10){"0"},i,"_test.dat"))
#   boxplot(pk, main = paste0("pk (2021), i=",i))
# }
# dev.off()
# 
# pdf("pdf/box_err_2021.pdf")
# for(i in 0:111){
#   err <- read.table(paste0("Mira-Titan-IV-data/Mira-Titan-2021/STEP",step,"/err_M",
#                            if(i<100){"0"},if(i<10){"0"},i,"_test.dat"))
#   boxplot(err, main = paste0("err (2021), i=",i))
# }
# dev.off()
