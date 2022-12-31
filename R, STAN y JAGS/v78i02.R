library("lcmm")

#############        Section 4      #######################

hlme1 <- hlme(Y ~ Time * X1, random = ~ Time, subject = "ID", ng = 1,
              data = data_hlme)
hlme2 <- hlme(Y ~ Time * X1, random = ~ Time, subject = "ID", ng = 2,
              data = data_hlme, mixture = ~ Time, classmb = ~ X2 + X3, B = hlme1)

lcmm1 <- lcmm(Ydep2 ~ poly(Time, degree = 2, raw = TRUE), random = ~ Time,
              subject = "ID", data = data_lcmm)
lcmm2 <- lcmm(Ydep2 ~ poly(Time, degree = 2, raw = TRUE), random = ~ Time,
              subject = "ID", data = data_lcmm, link = "5-quant-splines")

mlcmm1 <- multlcmm(
    Ydep1 + Ydep2 + Ydep3 ~ X1 * poly(Time, degree = 2, raw = TRUE),
    random = ~ Time, subject = "ID", data = data_lcmm,
    link = c("linear","3-quant-splines","3-quant-splines"))

jlcmm1 <- Jointlcmm(
    Ydep1 ~ X1 * Time, random = ~ Time, subject = "ID",
    survival = Surv(Tevent, Event) ~ X1 + X2, hazard = "3-quant-splines",
    data = data_lcmm)

jlcmm2 <- Jointlcmm(
    Ydep1 ~ Time * X1, random = ~ Time, subject = "ID", 
    mixture = ~ Time, survival = Surv(Tevent, Event) ~ X1 + mixture(X2),
    hazard = "3-quant-splines", hazardtype = "PH", ng = 2, 
    data = data_lcmm, B = jlcmm1)

#############     Section 6.1    #######################

library("NormPsy")
paquid$normMMSE <- normMMSE(paquid$MMSE)
paquid$age65 <- (paquid$age - 65) / 10
head(paquid)

#############     Section 6.2    ####################### 

m1a <- hlme(normMMSE ~ poly(age65, degree = 2, raw = TRUE) * CEP, 
            random = ~ poly(age65, degree = 2, raw = TRUE),
            subject = "ID", data = paquid, ng = 1)
summary(m1a)

WaldMult(m1a, pos = c(5, 6),name = "CEP interaction with age65 & age65^2")

m1 <- hlme(normMMSE ~ poly(age65, degree = 2, raw = TRUE) + CEP, 
           random = ~ poly(age65, degree = 2, raw = TRUE),
           subject = "ID", ng = 1, data = paquid)
m2 <- hlme(normMMSE ~ poly(age65, degree = 2, raw = TRUE) + CEP, 
           random = ~ poly(age65, degree = 2, raw = TRUE),
           mixture = ~ poly(age65, degree = 2, raw = TRUE),
           subject = "ID", ng = 2, data = paquid, B = m1)
m3 <- hlme(normMMSE ~ poly(age65, degree = 2, raw = TRUE) + CEP, 
           random = ~ poly(age65, degree = 2, raw = TRUE),
           mixture = ~ poly(age65, degree = 2, raw = TRUE),
           subject = "ID", ng = 3, data = paquid, B = m1)
m2b <- hlme(normMMSE ~ poly(age65, degree = 2, raw = TRUE) + CEP, 
            random = ~ poly(age65, degree = 2, raw = TRUE),
            mixture = ~ poly(age65, degree = 2, raw = TRUE),
            subject = "ID", ng = 2, data = paquid, 
            B = c(0, 60, 40, 0, -4, 0, -10, 10, 212.869397, -216.421323, 
                  456.229910, 55.713775, -145.715516, 59.351000, 10.072221))

set.seed(1)
m2c <- hlme(normMMSE ~ poly(age65, degree = 2, raw = TRUE) + CEP, 
            random = ~ poly(age65, degree = 2, raw = TRUE),
            mixture = ~ poly(age65, degree = 2, raw = TRUE),
            subject = "ID", data = paquid, ng = 2, B = random(m1))

m2d <- gridsearch(hlme(normMMSE ~ poly(age65, degree = 2, raw = TRUE) + CEP, 
                       random = ~ poly(age65, degree = 2, raw = TRUE),
                       mixture = ~ poly(age65, degree = 2, raw = TRUE),
                       subject = "ID", data = paquid, ng = 2, verbose = FALSE),
                  rep = 30, maxiter = 15, minit = m1)
m3b <- gridsearch(hlme(normMMSE ~ poly(age65, degree = 2, raw = TRUE) + CEP, 
                       random = ~ poly(age65, degree = 2, raw = TRUE),
                       mixture = ~ poly(age65, degree = 2, raw = TRUE),
                       subject = "ID", data = paquid, ng = 3, verbose = FALSE),
                  rep = 30, maxiter = 15, minit = m1)
summarytable(m1, m2, m2b, m2c, m2d, m3, m3b)

postprob(m2)

plot(m2)  
plot(m2, which = "fit", var.time = "age65", bty = "l", ylab = "normMMSE",
     xlab = "(age-65)/10", lwd = 2)
plot(m2, which = "fit", var.time = "age65", bty = "l", ylab= "normMMSE",
     xlab = "(age-65)/10", lwd = 2, marg = FALSE)

datnew <- data.frame(age = seq(65, 95, length = 100))
datnew$age65 <- (datnew$age - 65) / 10
datnew$CEP <- 0
CEP0 <- predictY(m2, datnew, var.time = "age")
datnew$CEP <- 1
CEP1 <- predictY(m2, datnew, var.time = "age")
plot(CEP1, lty = 1,lwd = 2, type = "l", col = 1:2 , ylim = c(20, 100), 
     bty = "l", xlab = "age in year", ylab = "normalized MMSE", 
     legend = NULL)
plot(CEP0, lty = 2, lwd = 2, type = "l", col = 1 : 2, ylim = c(20, 100), 
     add = TRUE)
legend(x = "topright", bty = "n", ncol = 3, lty = c(NA, NA, 1, 1, 2, 2), 
       col = c(NA, NA, 1, 2, 1, 2),
       legend = c("G=1 (12.4%):", "G=2 (87.6%):", "EL+", "EL+", "EL-", "EL-"),
       lwd = 2)

#############      Section 6.3      #######################

mlin <- lcmm(CESD ~ age65 * male, random = ~ age65, subject = "ID", 
             data = paquid)
mbeta <- lcmm(CESD ~ age65 * male, random = ~ age65, subject = "ID",
              data = paquid, link = "beta")
mspl <- lcmm(CESD ~ age65 * male, random = ~ age65, subject = "ID", 
             data = paquid, link = "splines")
mspl5q <- lcmm(CESD ~ age65 * male, random = ~ age65, subject = "ID", 
               data = paquid, link = "5-quant-splines")

mord0 <- lcmm(CESD ~ age65 * male, random = ~ -1, subject = "ID", 
              data = paquid, link = "thresholds")
binit <- vector("numeric", length = 56)
binit[1:6] <- mspl$best[1:6]
binit[7:56] <- mord0$best[4:53]
mord <- lcmm(CESD ~ age65 * male, random = ~ age65, subject = "ID",
             data = paquid, link = "thresholds", B = binit)
summary(mspl5q)
summary(mord)

col <- rainbow(5)
plot(mlin, which = "linkfunction", bty = "l", ylab = "CES-D", lwd = 2,
     col = col[1], xlab = "underlying latent process")
plot(mbeta, which = "linkfunction", add = TRUE, col = col[2], lwd = 2)
plot(mspl, which = "linkfunction", add = TRUE, col = col[3], lwd = 2)
plot(mspl5q, which = "linkfunction", add = TRUE, col = col[4], lwd = 2)
plot(mord, which = "linkfunction", add = TRUE, col = col[5], lwd = 2)
legend(x = "topleft", legend = c("linear", "beta",
          "splines (5equidistant)", "splines (5 at quantiles)", "thresholds"), 
          lty = 1, col = col, bty = "n", lwd = 2)
linkspl5q <- predictlink(mspl5q, ndraws = 2000)
plot(linkspl5q, add = TRUE, col = col[4], lty = 2)
legend(legend = c("95% confidence bands", "for splines at quantiles"),
           x = "left", lty = c(2, NA), col = c(col[4], NA), bty = "n", lwd = 1)

plot(mspl5q)
plot(mspl5q, which = "fit", var.time = "age65", xlab = "(age - 65) / 10",
     bty = "l", break.times = 8, ylab = "latent process", lwd = 2, marg = FALSE,
     ylim = c(-1, 2))

datnew <- data.frame(age = seq(65, 95, length = 100))
datnew$age65 <- (datnew$age - 65) / 10
datnew$male <- 0
women <- predictY(mspl5q, newdata = datnew, var.time = "age", draws = TRUE)
datnew$male <- 1
men <- predictY(mspl5q, newdata = datnew, var.time = "age", draws = TRUE)

plot(women, lwd = c(2, 1), type = "l", col = 6, ylim = c(0, 20),
     xlab = "age in year", ylab = "CES-D", bty = "l", legend = NULL)
plot(men, add = TRUE, col = 4, lwd = c(2, 1))
legend(x = "topleft", bty = "n", ncol = 2, lty = c(1, 1, 2, 2), 
     col = c(6, 4, 6, 4), lwd = c(2, 2, 1, 1),
     legend = c("women", "men", "   95% CI", "   95% CI"))

#############      Section 6.4      #######################

paquid$time <- (paquid$age - paquid$age_init)
paquid$age0_centered <- paquid$age_init - 75
mult <- multlcmm(MMSE + IST + BVRT ~ age0_centered + male + contrast(male) +
                     time + I(time^2 / 10),
                 random = ~ time + I(time^2 / 10), subject = "ID",
                 data = paquid, randomY = TRUE, cor = BM(time),
                 link = c("beta", "beta", "beta"))
summary(mult)

plot(mult, which = "linkfunction", col = c(1, 4, 6), lwd = 2)
CI <- predictlink(mult)
plot(CI, col = c(1, 4, 6), lwd = 2)
head(CI$pred)

VarExpl(mult, values = data.frame(time = 0))
VarExpl(mult, values = data.frame(time = 5))

#############      Section 6.5      #######################

paquidS <- paquid[paquid$agedem > paquid$age_init, ]

mj1 <- Jointlcmm(normMMSE ~ poly(age65, degree = 2, raw = TRUE) + CEP,
                 random = ~ poly(age65, degree = 2, raw = TRUE), 
                 survival = Surv(age_init, agedem, dem) ~ CEP + male, 
                 hazard = "Weibull", subject = "ID", data = paquidS, ng = 1)
mj2 <- Jointlcmm(normMMSE ~ poly(age65, degree = 2, raw = TRUE) + CEP,
                 mixture = ~ poly(age65, degree = 2, raw = TRUE), 
                 random = ~ poly(age65, degree = 2, raw = TRUE), 
                 survival = Surv(age_init, agedem, dem) ~ CEP + male,
                 hazard = "Weibull", subject = "ID", data = paquidS, ng = 2, B = mj1)
mj3 <- Jointlcmm(normMMSE ~ poly(age65, degree = 2, raw = TRUE) + CEP,
                 mixture = ~ poly(age65, degree = 2, raw = TRUE), 
                 random = ~ poly(age65, degree = 2, raw = TRUE), 
                 survival = Surv(age_init, agedem, dem) ~ CEP + male,
                 hazard = "Weibull", subject = "ID", data = paquidS, ng = 3, B = mj1)
mj4 <- Jointlcmm(normMMSE ~ poly(age65, degree = 2, raw = TRUE) + CEP,
                 mixture = ~ poly(age65, degree = 2, raw = TRUE), 
                 random = ~ poly(age65, degree = 2, raw = TRUE), 
                 survival = Surv(age_init, agedem, dem) ~ CEP + male,
                 hazard = "Weibull", subject = "ID", data = paquidS, ng = 4, B = mj1)

summarytable(mj1, mj2, mj3, mj4)

Binit <- rep(0, length(mj2$best) + 6)
Binit[c(2, 5:10, 12, 13, 15, 16, 18, 19:(length(Binit)))] <- mj2$best
Binit[c(1, 3, 4, 11, 14, 17)] <- c(0, 0.11, 4, 70, 0, 0)
mj3b <- Jointlcmm(normMMSE ~ poly(age65, degree = 2, raw = TRUE) + CEP,
                  mixture = ~ poly(age65, degree = 2, raw = TRUE), 
                  random = ~ poly(age65, degree = 2, raw = TRUE), 
                  survival = Surv(age_init, agedem, dem) ~ CEP + male,
                  hazard = "Weibull", subject = "ID", data = paquidS, ng = 3, B = Binit)

Binit <- rep(0, length(mj3b$best) + 2 + 3 + 1)
Binit[c(1, 2, 4:7, 10:15, 17:19, 21:23, 25:length(Binit))] <- mj3b$best
Binit[c(3, 8, 9, 16, 20, 24)] <- c(0, 0.1, 10, 60, 5, -10)
mj4b <- Jointlcmm(normMMSE ~ poly(age65, degree = 2, raw = TRUE) + CEP,
                  mixture = ~ poly(age65, degree = 2, raw = TRUE), 
                  random = ~ poly(age65, degree = 2, raw = TRUE), 
                  survival = Surv(age_init, agedem, dem) ~ CEP + male,
                  hazard = "Weibull", subject = "ID", data = paquidS, ng = 4, B = Binit)
mj3c <- gridsearch(rep = 30, maxiter = 15, minit = mj1,
                   Jointlcmm(normMMSE ~ poly(age65, degree = 2, raw = TRUE) + CEP,
                             mixture = ~ poly(age65, degree = 2, raw = TRUE), 
                             random = ~ poly(age65, degree = 2, raw = TRUE), 
                             survival = Surv(age_init, agedem, dem) ~ CEP + male,
                             hazard = "Weibull", subject = "ID", data = paquidS, ng = 3,
                             verbose = FALSE)) 
mj4c <- gridsearch(rep = 30, maxiter = 15, minit = mj1,
                   Jointlcmm(normMMSE ~ poly(age65, degree = 2, raw = TRUE) + CEP,
                             mixture = ~ poly(age65, degree = 2, raw = TRUE), 
                             random = ~ poly(age65, degree = 2, raw = TRUE), 
                             survival = Surv(age_init, agedem, dem) ~ CEP + male,
                             hazard = "Weibull", subject = "ID", data = paquidS, ng = 4,
                             verbose = FALSE))

summarytable(mj1, mj2, mj3, mj3b, mj3c, mj4, mj4b, mj4c)

summary(mj4b)

plot(mj4b, which = "fit", var.time = "age65", marg = FALSE, break.times = 10,
     bty = "l", ylab = "normMMSE", xlab = "Age in decades from 65 years")

datnew <- data.frame(age65 = seq(0, 3, length = 100))
datnew$male <- 0
datnew$CEP <- 0
par(mfrow = c(1, 2))
mj4b.pred <- predictY(mj4b, newdata = datnew, var.time = "age65")
plot(mj4b.pred, bty = "l", ylim = c(0, 80), legend.loc = "bottomleft",
    ylab = "normMMSE", xlab = "age in decades from 65 years", lwd = 2)
plot(mj4b, which = "survival", lwd = 2, legend.loc = FALSE, bty = "l", 
    xlab = "age in years", ylab = "dementia-free probability")

postprob(mj4b)

landmark <- c(70, 72, 75, 77, 80, 82, 85, 87, 90)
epoce1 <- epoce(mj1, pred.times = landmark, var.time = "age65",
                fun.time = function(x) { 10 * x + 65 })
epoce2 <- epoce(mj2, pred.times = landmark, var.time = "age65",
                fun.time = function(x) { 10 * x + 65 })
epoce3 <- epoce(mj3b, pred.times = landmark, var.time = "age65",
                fun.time = function(x) { 10 * x + 65 })
epoce4 <- epoce(mj4b, pred.times = landmark, var.time = "age65",
                fun.time = function(x) { 10 * x + 65})
diff23 <- Diffepoce(epoce2, epoce3)
diff34 <- Diffepoce(epoce3, epoce4)
par(mfrow = c(1, 2))
plot(epoce1, ylim = c(0.5, 1.5), main = "cross-validated EPOCE estimates",
    bty = "l")
plot(epoce2, add = TRUE, col = 2, lty = 2)
plot(epoce3, add = TRUE, col = 3, lty = 3)
plot(epoce4, add = TRUE, col = 4, lty = 4)
legend("topright", legend = c("G=1", "G=2", "G=3", "G=4"), col = 1:4, 
       lty = 1:4, bty = "n")
plot(diff23, main = "Difference in EPOCE estimates", lty = c(1, 2, 2),
     pch = 20, ylim = c(-0.05, 0.30), bty = "l")
plot(diff34, add = TRUE, main = "Difference in EPOCE estimates", col = 4, 
     lty = c(1, 2, 2), pch = 18)
legend("topleft", legend = c("G=2/G=3", "G=3/G=4", "95%TI", "95%TI"), 
       ncol = 2, col = c(1, 4, 1, 4), lty = c(1, 1, 2, 2), 
       pch = c(20, 18, 20, 18), bty = "n")

paq72 <- paquid[which(paquid$ID == 72), ]
dynp <- dynpred(mj4b, paq72, landmark = c(80, 90), var.time = "age65", 
                horizon = c(1, 3, 5, 8, 9),
                fun.time = function(x) { 10 * x + 65 }, draws = TRUE)
plot(dynp, landmark = 80, ylim = c(55, 85, 0, 1), col = 1, pch = 20, 
     ylab = "normMMSE", main = "At landmark age 80", xlab = "age in years")
plot(dynp, landmark = 90, ylim = c(55, 85, 0, 1), col = 1, pch = 20, 
     ylab = "normMMSE", main = "At landmark age 90", xlab = "age in years")







