library( gasplit) # of course!

library( offarray) # or die; needed for data org in this example
# ... but you tidyverse people don't generally need it for gasplitting

sim <- make_fake_ppns( 
    meanE= +2, # not-brilliant separation
    mean_nsamp=200, # reasonable samp sizes 
    interpow=0) # no interactions
hist( sim$DISCOSTAT, nc=40)

# Next blob of code, for fitting model & summarizing reslts
# is basically the same as:
# test <- test_gasplit( sim)
# but inlined here so you can see the steps

# Set up the conditional distros, based on sim@truth$meanE
e <- list2env( sim@truth, parent=.GlobalEnv)
d1 <- function( x) dt( x - meanE, df=df_t)
d2 <- function( x) dt( x + meanE, df=df_t) # meanW === -(meanE)
environment( d1) <- environment( d2) <- e

# Fit the model: could use gasplit() instead of gasplit2(), since no s() terms
gg1 <- gasplit2( DISCOSTAT ~ Y+Z-1, sim, d1=d1, d2=d2)

# Organize fitted ppn into an array like sim@truth$ppn1
# Could possibly use predict()... but deviousness would be needed

sim$fit_ppnE <- gg1$ppn
fitsim <- sim[ !duplicated( sim[ cq( Y, Z)]), cq( Y, Z, fit_ppnE)]
fit_ppnE <- sim@truth$ppnE * 0
fit_ppnE[ MATSUB=fitsim[ cq( Y, Z)]] <- fitsim$fit_ppnE

test <- returnList( 
    samp_ppnE= sim@truth$samp_ppnE, 
    tru_ppnE=sim@truth$ppnE, 
    fit_ppnE,
    gg1
  )
# ... end of fitting & summarizing

# Let's plot some stuff:
with( test, plot( tru_ppnE, fit_ppnE))
abline( 0,1)

# Posterior distro per point (also its derivs wrto beta params)
posty <- posterior( test$gg1, dbeta=TRUE) 

# Next plot is not very meaningful... NB that ratio of two t-distros with same var eventually *flattens* out at extreme values; basically saying "NFI but I don't believe either of them". That's why things aren't monotonic (I think!)

plot( sim$DISCOSTAT, logit( posty[,1]), pch='.')
abline( 0,1,col='red') # ... though there's no particular reason for y==x, I think

# Prior vs posterior: the willing eye can see an upwards trend...
plot( jitter( test$gg1$ppn), posty[,1])

# ... and a regression of posterior on prior confirms it 
# NB intercept & slope!
lm( posty[,1] ~ test$gg1$ppn)

# How does per-sample evidence affect the posterior?
DISCOSTAT_y <- with( test$gg1, d1( G$y) / d2( G$y))
plot( posty[,1], log( DISCOSTAT_y))

# Check that nothing changes if you change the sign of the stat, *and* the signs of the conditional means...
simneg <- sim
simneg$DISCOSTAT <- (-sim$DISCOSTAT)
simneg@truth$meanE <- (-sim@truth$meanE)
hist( simneg$DISCOSTAT, nc=40)
testneg <- test_gasplit( simneg) # for brevity
with( testneg, plot( tru_ppnE, fit_ppnE))
abline( 0,1)

plot( test$gg1$ppn, testneg$gg1$ppn)

# Posterior. 
