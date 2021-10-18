# Check of equations page 8
# In ineqx.R, go up until dat.f_cf
# Then, we see that the results from
# dW.out, dB.out, and dCD.out are the same as these equations from page 8

dat <-
  dat.f_cf %>%
  dplyr::mutate(across(c(n.f, n.cf), ~./sum(.)*10)) %>%
  dplyr::filter(time==2)

dB <- with(dat, sum(n.cf*((mu.cf+beta.f-sum(n.cf*(mu.cf+beta.f)))^2-(mu.cf+beta.cf-sum(n.cf*(mu.cf+beta.cf)))^2)))

dW <- with(dat, sum(n.cf*(lambda.f^2-lambda.cf^2+2*sigma.cf*(lambda.f-lambda.cf))))

dC <- with(dat, sum((n.f-n.cf)*((mu.f+beta.f-sum(n.f*(mu.f+beta.f)))^2+(sigma.f+lambda.f)^2)))

dD.B <- with(dat, sum(n.f*((mu.f+beta.f-sum(n.f*(mu.f+beta.f)))^2-(mu.cf+beta.f-sum(n.f*(mu.cf+beta.f)))^2+(mu.cf-sum(n.f*mu.cf))^2-(mu.f-sum(n.f*mu.f))^2)))
dD.W <- with(dat, sum(n.f*(2*lambda.f*(sigma.f-sigma.cf))))
dD.B+dD.W

dT <- with(dat,
           sum(n.f*((mu.f+beta.f-sum(n.f*(mu.f+beta.f)))^2-(mu.f-sum(n.f*(mu.f)))^2+(2*lambda.f*sigma.f+lambda.f^2)))-
             sum(n.cf*((mu.cf+beta.cf-sum(n.cf*(mu.cf+beta.cf)))^2-(mu.cf-sum(n.cf*(mu.cf)))^2+(2*lambda.cf*sigma.cf+lambda.cf^2))))

all.equal(dT, dB+dW+dC+dD.B+dD.W)
