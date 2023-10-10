library(funFEM)
### Real Case study
library(BayesFMMM)

setwd("")

#################################################################
## Change relevant directories and make folders before running ##
#################################################################


### Peak alpha data
library(pracma)
library(gridExtra)
# Subject ID
subj_id <- sort(c(10,	11,	13,	14,	15,	23,	26,	30,	31,	35,	48,	49,	50,
                  53,	54,	55,	161,165,	184,	188,	189,	195,	201,
                  # 202,	excluded due to low counts
                  207,	210,	213,	214,	242,	255,	261,	282,	283,
                  284,	286,	287,	289,	290,	343,	351,	2,	3,	5,	6,
                  7,	8,	9,	12,	18,	19,	22,	24,	25,	27,	33,	34,	37,	38,
                  40,	41,	42,	43,	44,	47,	51,	401,	405,	406,	408,	411,
                  415,	416,	417,	418,	423,	426,	427,	430,
                  #431,	excluded due to low counts
                  433,	436,	438,	439,	440,	442,	444,	445,	446,	447,
                  448,	450,	451,	452,	453,	3019,	3024,	3026,	3029,	3032))
# Channel ID (order of chan_id corresponds to 1:25 labeling of regions)
chan_id <- c('FP1', 'FP2','F9','F7','F3','Fz','F4','F8','F10','T9','T7',
             'C3','CZ','C4','T8','T10','P9','P7','P3','PZ','P4','P8','P10','O1','O2')

# Demographic Data
demDat <- read.csv(file='demographic_data.csv', header = TRUE)
colnames(demDat) <- c("ID", "Gender", "Age", "Group", "VIQ", "NVIQ")
demDat <- demDat[which(demDat$ID %in% subj_id), ]

# Peak Alpha Data
load("pa.dat.Rdata")
# ID: subject ID
# group: TD(1) or ASD (2)
# func: frequency domain
# reg: electrode (order corresponds to chan_id above)
# Age: age in months
# y: alpha spectra density
out1 <- unique(pa.dat$func)
out3 <- unique(pa.dat$reg)
matplot(matrix(pa.dat$y, nrow = length(out1)), type = "l") # data
trapz(out1, pa.dat$y[1:33]) # all functional observations integrate to 1 (normalized across electordes, subjects)

### Convert to wide format
Y <- pa.dat
## paper used T8 electrode
Y <- Y[Y$reg == 15,]
Y$ID <- paste(Y$ID, Y$reg, sep = ".")
Y <- reshape(Y[,c(1,3,6)], idvar = "ID", timevar = "func", direction = "wide")
Y <- Y[,-1]
Y <- as.matrix(Y)
matplot(t(Y), type = "l")


######################
## Factor Analysis ###
######################
Y <- as.data.frame(Y)
fa_model <- factanal(Y[,1:32]- colMeans(Y[,1:32]),factors = 2)
plot(fa_model$loadings[,1])
plot(fa_model$loadings[,2])
# time <-(rep(seq(6,14, 0.25), nrow(Y)))
# x <- as.vector(t(Y))
#
# curve <- rep(0, length(time))
# for(i in 1:nrow(Y)){
#   curve[((i-1)*33 + 1):(i*33)] <- i
# }
#
# data <- list(x, time, curve)
# names(data) <- c("x", "time", "curve")
#
# m2 <- mocca(data, K = 2,q=8,h=1,lambda=1e-10,n.cores=2,EMstep.tol=1e-3)
# m2


##########################
## Clustering Analysis ###
##########################
basis <- create.bspline.basis(c(6,14), basis = 50, norder = 6)
fdobj <- smooth.basis(seq(6,14, 0.25), t(Y), basis)$fd
res<- funFEM(fdobj, K = 4)
par(mfrow=c(1,2))
#plot(fdobj, color=res$cls)
fdmeans <- fdobj
fdmeans$coefs <- t(res$prms$my)
plot(fdmeans, col=1:max(res$cls), lwd=2, xlab = "Frequency", ylab = " ")

matplot(seq(6,14,0.25), t(Y[1:10,]), type = "l", col = res$cls[1:10], xlab = "Frequency", ylab = " " )

res1<- funFEM(fdobj, K = 3)
res2<- funFEM(fdobj, K = 4)
res3<- funFEM(fdobj, K = 5)
res4<- funFEM(fdobj, K = 6)
res5<- funFEM(fdobj, K = 7)
par(mfrow=c(1,2))
#plot(fdobj, color=res$cls)
fdmeans <- fdobj
fdmeans$coefs <- t(res$prms$my)
plot(fdmeans, col=1:max(res$cls),  lty = 1:max(res$cls), lwd=2, xlab = "Frequency", ylab = "Relative Power", cex.lab = 1.25, cex.axis = 1.25)

matplot(seq(6,14,0.25), t(Y[1:10,]), type = "l", col = res$cls[1:10], lty = res$cls[1:10], xlab = "Frequency", ylab = "Relative Power", lwd = 1.5,  cex.lab = 1.25, cex.axis = 1.25)

fun_pca <- pca.fd(fdobj,nharm = 2)
plot(fun_pca$harmonics)

library(fdapace)
y_list <- split(Y, seq(nrow(Y)))
time <- seq(6, 14, 0.25)
time <- rep(list(time), 97)

FPCAobj <- FPCA(y_list, time)
plot(FPCAobj)

signal_ratio <- readRDS("./signal_ratio.rds")
demDat <- cbind(demDat, signal_ratio)

signal_ratio_ph <- signal_ratio$id
pa_id <- pa.dat$ID[pa.dat$reg == 15 & pa.dat$func == 6]
for(i in 1:length(signal_ratio)){
  signal_ratio_ph[i] <- signal_ratio[demDat$ID == pa_id[i]]
}

library(ggplot2)

df <- cbind(FPCAobj$xiEst[,1:2], signal_ratio)
df <- as.data.frame(df)
colnames(df) <- c("FPC 1", "FPC 2", "Sig_ratio")

ggplot(df, aes(`FPC 1`, `FPC 2`, colour = Sig_ratio)) + geom_point() + scale_color_gradient(low="blue", high="red") + theme_bw()

