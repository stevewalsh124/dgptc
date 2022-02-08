#################################
# Inter-batch chain comparisons #
# Steve Walsh                   #
# February 8 2022               #
#################################

pdf("pdf/batch_compare.pdf")

# Load two different batches of summary stats (18 TCs per batch)
load("rda/FL_summaries/pmx_means_50000_thin100")
load("rda/FL_summaries/pmx_medians_50000_thin100")
pmx_means2 <- pmx_means
pmx_medians2 <- pmx_medians

load("rda/FL_summaries/pmx_means_50000_thin100orig50k")
load("rda/FL_summaries/pmx_medians_50000_thin100orig50k")

par(mfrow=c(2,3))
for (i in 1:6) {plot(pmx_means[,i], pmx_means2[,i],
                     main = paste(colnames(pmx_means)[i],
                                  round(cor(pmx_means[,i],pmx_means2[,i]),3))); abline(0,1)}
mtext("Means vs Means2", side = 3, line = -22, outer = TRUE)

for (i in 1:6) {plot(pmx_medians[,i], pmx_medians2[,i],
                     main = paste(colnames(pmx_means)[i],
                                  round(cor(pmx_medians[,i],pmx_medians2[,i]),3))); abline(0,1)}
mtext("Medians vs Medians2", side = 3, line = -22, outer = TRUE)

for (i in 1:6) {plot(pmx_means[,i], pmx_medians[,i],
                     main = paste(colnames(pmx_means)[i],
                                  round(cor(pmx_means[,i],pmx_medians[,i]),3))); abline(0,1)}
mtext("Means vs Medians", side = 3, line = -22, outer = TRUE)

for (i in 1:6) {plot(pmx_means2[,i], pmx_medians2[,i],
                     main = paste(colnames(pmx_means)[i],
                                  round(cor(pmx_means2[,i],pmx_medians2[,i]),3))); abline(0,1)}
mtext("Means2 vs Medians2", side = 3, line = -22, outer = TRUE)

dev.off()
