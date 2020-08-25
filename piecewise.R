library(piecewiseSEM)


modelA <- psem(lm(Moose.Density ~ LAI, df), lm(WolfDensit ~ Moose.Density, df), lm(caribou.lambda ~ disturb.p, df))
modelB <- psem(lm(Moose.Density ~ LAI, df), lm(WolfDensit ~ Moose.Density + disturb.p, df), lm(caribou.lambda ~ WolfDensit, df))
modelC <- psem(lm(Moose.Density ~ LAI, df), lm(WolfDensit ~ Moose.Density, df), lm(caribou.lambda ~ Moose.Density+LAI, df))
modelC <- update(modelC, disturb.p ~ 1)
modelD <- psem(lm(Moose.Density ~ LAI, df), lm(WolfDensit ~ Moose.Density, df), lm(caribou.lambda ~ WolfDensit, df))
modelD <- update(modelD, disturb.p ~ 1)
modelE <- psem(lm(Moose.Density ~ LAI, df), lm(WolfDensit ~ Moose.Density, df), lm(caribou.lambda ~ WolfDensit +disturb.p, df))
modelF <- psem(lm(Moose.Density ~ LAI, df), lm(WolfDensit ~ Moose.Density, df), lm(caribou.lambda ~ LAI, df))
modelF <- update(modelF, disturb.p ~ 1)


summary(modelA, .progressBar = F)
summary(modelB, .progressBar = F)
summary(modelC, .progressBar = F)
summary(modelD, .progressBar = F)
summary(modelE, .progressBar = F)
summary(modelF, .progressBar = F)
plot(modelE)

AIC(modelD, modelA)
AIC(modelD, modelB)
AIC(modelD, modelC)
AIC(modelD, modelE)
AIC(modelD, modelF)


modelA <- psem(lm(Moose.Density ~ LAI, df), lm(WolfDensit ~ Moose.Density, df), lm(caribou.lambda ~ disturb.p, df))
modelB <- psem(lm(Moose.Density ~ LAI, df), lm(WolfDensit ~ Moose.Density + disturb.p, df), lm(caribou.lambda ~ WolfDensit, df))
modelC <- psem(lm(Moose.Density ~ LAI, df), lm(WolfDensit ~ Moose.Density, df), lm(caribou.lambda ~ Moose.Density+LAI, df), lm(LAI ~ disturb.p, df))
modelD <- psem(lm(Moose.Density ~ LAI, df), lm(WolfDensit ~ Moose.Density, df), lm(caribou.lambda ~ WolfDensit, df), lm(LAI ~ disturb.p, df))
modelE <- psem(lm(Moose.Density ~ LAI, df), lm(WolfDensit ~ Moose.Density, df), lm(caribou.lambda ~ WolfDensit +LAI, df), lm(LAI ~ disturb.p, df))
modelF <- psem(lm(Moose.Density ~ LAI, df), lm(WolfDensit ~ Moose.Density, df), lm(caribou.lambda ~ LAI, df), lm(LAI ~ disturb.p, df))


plot(modelD)
plot(modelB)


summary(modelA, .progressBar = F)
summary(modelB, .progressBar = F)
summary(modelC, .progressBar = F)
summary(modelD, .progressBar = F)
summary(modelE, .progressBar = F)
summary(modelF, .progressBar = F)

AIC(modelD, modelA)
AIC(modelD, modelB)
AIC(modelD, modelC)
AIC(modelD, modelE)
AIC(modelD, modelF)


modelG <- psem(lm(Moose.Density ~ LAI +  disturb.p, df), lm(WolfDensit ~ Moose.Density, df), lm(caribou.lambda ~ WolfDensit, df))
AIC(modelD, modelG)