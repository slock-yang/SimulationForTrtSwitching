library(stargazer)

dat = NULL
out = read_json("SimulationResults/Outcome_DR_1600.json", simplifyVector = T)
apply(out, 1, mean)[1]
a = apply(out - matrix(c(0.1, 0.075), 
                   nrow = 2, ncol = nrep), 1, mean)[1]
b = apply(out - matrix(c(0.1, 0.075), 
                   nrow = 2, ncol = nrep), 1, function(d)sqrt(sum(d^2)/(nrep-1)))[1]
c = apply(out, 1, sd)[1]

dat = cbind(dat, c(a, b, c))

out = read_json("SimulationResults/Outcome_DR_1600_ahaz_censoring.json", simplifyVector = T)
apply(out, 1, mean)[1]
a = apply(out - matrix(c(0.1, 0.075), 
                       nrow = 2, ncol = nrep), 1, mean)[1]
b = apply(out - matrix(c(0.1, 0.075), 
                       nrow = 2, ncol = nrep), 1, function(d)sqrt(sum(d^2)/(nrep-1)))[1]
c = apply(out, 1, sd)[1]

dat = cbind(dat, c(a, b, c))

out = read_json("SimulationResults/Outcome_DR_1600_ahaz_delete.json", simplifyVector = T)
apply(out, 1, mean)[1]
a = apply(out - matrix(c(0.1, 0.075), 
                       nrow = 2, ncol = nrep), 1, mean)[1]
b = apply(out - matrix(c(0.1, 0.075), 
                       nrow = 2, ncol = nrep), 1, function(d)sqrt(sum(d^2)/(nrep-1)))[1]
c = apply(out, 1, sd)[1]

dat = cbind(dat, c(a, b, c))

out = read_json("SimulationResults/Outcome_DR_1600_itt_tv.json", simplifyVector = T)
out1 = out$itt
out2 = out$tv
a = apply(out1 - matrix(c(0.1), 
                       nrow = 1, ncol = nrep), 1, mean)[1]
b = apply(out1 - matrix(c(0.1), 
                       nrow = 1, ncol = nrep), 1, function(d)sqrt(sum(d^2)/(nrep-1)))[1]
c = apply(out1, 1, sd)[1]

dat = cbind(dat, c(a, b, c))

a = apply(out2 - matrix(c(0.1, 0.075), 
                        nrow = 2, ncol = nrep), 1, mean)[1]
b = apply(out2 - matrix(c(0.1, 0.075), 
                        nrow = 2, ncol = nrep), 1, function(d)sqrt(sum(d^2)/(nrep-1)))[1]
c = apply(out2, 1, sd)[1]


dat = cbind(dat, c(a, b, c))

stargazer(dat, digits = 4)


out = read_json("SimulationResults/setting2/setting2outcome1600.json", simplifyVector = T)

dat = NULL
a = apply(out$dr[, out$Con]- rep(0.1, 7), 1, mean)[1]
b = apply(out$dr[, out$Con] - rep(0.1, 7), 1, function(d)sqrt(sum(d^2)/(nrep-1)))[1]
c = apply(out$dr[, out$Con], 1, sd)[1]
dat = cbind(dat, c(a, b, c))

a = apply(out$ahaz- rep(0.1, 7), 1, mean)[1]
b = apply(out$ahaz - rep(0.1, 7), 1, function(d)sqrt(sum(d^2)/(nrep-1)))[1]
c = apply(out$ahaz, 1, sd)[1]
dat = cbind(dat, c(a, b, c))

a = apply(out$delete- rep(0.1, 7), 1, mean)[1]
b = apply(out$delete - rep(0.1, 7), 1, function(d)sqrt(sum(d^2)/(nrep-1)))[1]
c = apply(out$delete, 1, sd)[1]
dat = cbind(dat, c(a, b, c))

a = apply(out$itt- rep(0.1, 7), 1, mean)[1]
b = apply(out$itt- rep(0.1, 7), 1, function(d)sqrt(sum(d^2)/(nrep-1)))[1]
c = apply(out$itt, 1, sd)[1]
dat = cbind(dat, c(a, b, c))

a = apply(out$timeVar- rep(0.1, 7), 1, mean)[1]
b = apply(out$timeVar- rep(0.1, 7), 1, function(d)sqrt(sum(d^2)/(nrep-1)))[1]
c = apply(out$timeVar, 1, sd)[1]
dat = cbind(dat, c(a, b, c))
stargazer(dat, digits = 4)

