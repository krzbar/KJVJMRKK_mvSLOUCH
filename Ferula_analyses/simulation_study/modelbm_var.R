load("Ferula_Data.RData")
seed2 <- 1
set.seed(seed2)
BMmodel <- mvSLOUCH::BrownianMotionModel(ferulatreeape, feruladata_ME, M.error=M.error)
