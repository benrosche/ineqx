



dat1 <- crDat(2, 3, c("1000", "1000", "1000"), "10+10*(group==1)*x+250*(group==2)*x+300*(group==3)*x*year", "0.1", "within", X_const = T, seed=T, sav=F)

wibe(y="inc", groupvar = "group", timevar="year", dat = dat1)


test <- ineqx(x="c.x", y="c.inc", ystat = "CV2", groupvar = "c.group", timevar = "c.year", cf=list(n=c(1000,1000,100),
                                                                                                   mu=c(10,10,10),
                                                                                                   sigma=c(0.1,0.1,0.1),
                                                                                                   beta=c(0,0,0),
                                                                                                   lambda=c(0,0,0)), dat = dat1 %>% dplyr::mutate(group=1))

plot(test, type = "dPA")

plot(test, type = "dT")

