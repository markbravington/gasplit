# This is package gasplit 

"gasplit" <-
function( 
    formula, 
    data, 
    d1, d2, 
    start=0.001, 
    predict_from_previous=NULL
){
  if( !is.null( predict_from_previous)){
    # I think we have to actually "fit" the gam, so that we can use predict()
    # Result is garbage here!
    G <- gam( G=predict_from_previous$G)
    X <- predict( G, newdata=data, type='lpmatrix')
    pifodds <- rep( 0, nrow( X)) # just so nlglk can run (once only)
  } else { # normally... set up to fit gam, but don't do it
    G <- gam( formula=formula, data=data, fit=FALSE) # stuff needed for fit
    X <- G$X
    pifodds <- d1( G$y) / d2( G$y) - 1
  }
  
  ppn <- 0*pifodds - 1
  
  nlglk <- function( pars){
    ppn <<- plogis( X %*% pars)
    # prob = d1 * ppn + d2 * (1-ppn)
    # = ppn * (d1 - d2) + d2
    # = d2 * ( ppn * (d1/d2-1) + 1)
    # log( prob) = log( d2) + log( 1 + ppn * (d1/d2-1))
    lprob <- log( 1+pifodds*ppn)
  return( -sum( lprob))
  }

  if( !is.null( predict_from_previous)){
    nlglk( predict_from_previous$par)
return( ppn)
  }

  if( length( start)==1){
    start <- rep( start, ncol( X))
  }  # otherwise, user must ensure start is correct length
stopifnot( length( start)==ncol( X))

  fitto <- optim( start, nlglk, method='BFGS', 
      control=list( trace=5))
returnList( par=fitto$par, ppn, G)
}


"gasplit_gam" <-
function( formula, data, d1, d2){
## Smooth-ready version using RTMB

  G <- gam( formula=formula, data=data, fit=FALSE) # stuff needed for fit

  if( !length( G$S)){
return( gasplit_nogam( d1=d1, d2=d2, G=G))
  }
  
  X <- G$X
  pifodds <- d1( G$y) / d2( G$y) - 1
  S <- G$S
  
  ppn <- 0*pifodds - 1
  
  nlglk <- function( parmy){
    # Data:
    ppn <<- plogis( X %*% pars)
    lglk <- sum( log( 1+pifodds*ppn))
    
    # Priors (smooths):
    
  return( -sum( lprob))
  }
  
  fitto <- optim( rep( 1, ncol( X)), nlglk, method='BFGS', 
      control=list( trace=5))
returnList( par=fitto$par, ppn)
}


"gasplit_nogam" <-
function( formula, data, d1, d2, G){
  if( is.null( G)){
    G <- gam( formula=formula, data=data, fit=FALSE) 
    # ...stuff needed for fit
stopifnot( length( G$S)==0)
  }

  X <- G$X
  pifodds <- d1( G$y) / d2( G$y) - 1
  
  ppn <- 0*pifodds - 1
  
  nlglk <- function( pars){
    ppn <<- plogis( X %*% pars)
    # prob = d1 * ppn + d2 * (1-ppn)
    # = ppn * (d1 - d2) + d2
    # = d2 * ( ppn * (d1/d2-1) + 1)
    # log( prob) = log( d2) + log( 1 + ppn * (d1/d2-1))
    lprob <- log( 1+pifodds*ppn)
  return( -sum( lprob))
  }
  
  fitto <- optim( rep( 1, ncol( X)), nlglk, method='BFGS', 
      control=list( trace=5))
returnList( par=fitto$par, ppn)
}


"make_fake_ppns" <-
function(
  nyears= 4,
  nzones= 3,
  interpow= 0.1,
  df_t= 5,
  mean_nsamp= 100,
  meanE= 1,
  prange= 1,
  seed= 2
){
stopifnot( require( 'offarray'))

  YEARS <- 'Y' %&% (2000 + 1:nyears)
  ZONES <- LETTERS[ 1:nzones]
  
  rs <- .Random.seed
  on.exit( .Random.seed <<- rs)
  
  set.seed( 1)
  yeff <- offarray( prange*runif( nyears, -1, 1), 
      dimseq=list( YEARS))
  zeff <- 2 * offarray( prange*runif( nzones, -1, 1), 
      dimseq=list( ZONES)) # more oomph than year
  yzeff <- offarray( interpow * prange * runif( nyears*nzones, -1, 1), 
      dimseq=list( YEARS, ZONES)) # less oomph (presumably)
  
  # Could allow different numbers of total obs per stratum, but...
  extract.named( autoloop( Y=YEARS, Z=ZONES, {
      eff <- yeff[ Y] + zeff[ Z] + yzeff[ Y, Z];
      ppnE <- inv.logit( eff);
      nobsE <- rbinom( length( eff), size=mean_nsamp, prob=ppnE);
      nobsW <- mean_nsamp - nobsE
      samp_ppnE <- nobsE / mean_nsamp
    returnList( eff, ppnE, samp_ppnE, nobsE, nobsW)
    }))
  nstrat <- length( eff)
  
  # data.frame version combining the above:
  df <- as.data.frame( samp_ppnE, name_of_response='samp_ppnE')
  dfE <- df[ rep( seq_len( nstrat), c( nobsE)), ]
  dfE$whicho <- 'E'
  dfE$LGLR <- rt( nrow( dfE), df=df_t) + meanE
  
  meanW <- (-meanE)  # opposite mean
  dfW <- df[ rep( seq_len( nstrat), c( nobsW)),]
  dfW$whicho <- 'W'
  dfW$LGLR <- rt( nrow( dfW), df=df_t) + meanW
  
  dfall <- rbind( dfE, dfW)
  rownames( dfall) <- NULL # they are just annoying
  
  dfall@truth <- returnList( yeff, zeff, yzeff, df_t, 
      meanE, ppnE, samp_ppnE, prange, seed)
return( dfall)
}


"test_gasplit" <-
function( sim=NULL, ...){
  if( is.null( sim)){
    sim <- make_fake_ppns( ...)
  }
  simpure <- sim
  
  e <- list2env( sim@truth, parent=.GlobalEnv)
  d1 <- function( x) dt( x - meanE, df=df_t)
  d2 <- function( x) dt( x + meanE, df=df_t) # meanW === -(meanE)
  environment( d1) <- environment( d2) <- e
  gg1 <- gasplit( LGLR ~ Y+Z-1, sim, d1=d1, d2=d2)
  
  # Organize fitted ppn into an array like sim@truth$ppn1
  # Could possibly use predict()... but deviousness would be needed
  
  sim$fit_ppnE <- gg1$ppn
  fitsim <- sim[ !duplicated( sim[ cq( Y, Z)]), cq( Y, Z, fit_ppnE)]
  fit_ppnE <- sim@truth$ppnE * 0
  fit_ppnE[ MATSUB=fitsim[ cq( Y, Z)]] <- fitsim$fit_ppnE
  
  # fit_ppnE <- d2a( fitsim, data.col='fit_ppnE')
  
returnList( 
    samp_ppnE= sim@truth$samp_ppnE, 
    tru_ppnE=sim@truth$ppnE, 
    fit_ppnE,
    simpure,
    gg1
  )
}

