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
    ppn <- plogis( X %*% pars)
return( ppn)
  }

  # Else (normally) we want to fit. Set up to fit gam, but don't do it
  G <- gam( formula=formula, data=data, fit=FALSE) # stuff needed for fit
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

  if( length( start)==1){
    start <- rep( start, ncol( X))
  }  # otherwise, user must ensure start is correct length
stopifnot( length( start)==ncol( X))

  fitto <- optim( start, nlglk, method='BFGS', 
      control=list( trace=5))

return( c( returnList( 
    beta=fitto$par, ppn, G),
    fitto[ cq( convergence, message, evaluations)]
  ))
}


"gasplit2" <-
function( 
  formula, 
  data, 
  d1, d2, 
  predict_from_previous=NULL,
  ... # for gam()
){
## Allow REs (eg to do splines) and mgcv-style stuff
stopifnot( require( RTMB))

  if( !is.null( predict_from_previous)){
    # I think we have to actually "fit" the gam, so that we can use predict()
    # Result is garbage here!
    # Maybe this could be sped up by using existing estimates as startvals
    G <- gam( G=predict_from_previous$G)
    X <- predict( G, newdata=data, type='lpmatrix')
    ppn <- plogis( X %**% predict_from_previous$beta)
return( ppn)
  }
  
  # Else (ie normally), we fit
  # Stuff needed for fit, but don't actually fit:
  G <- gam( formula=formula, data=data, fit=FALSE, ...) 
  X <- G$X
  y <- G$y
  Slengths <- unlist( FOR( G$smooth, length( .$S)))
  Slist <- list()
  for( i in seq_along( G$smooth)){
    Si <- G$smooth[[i]]$S
    if( length( Si)){ # non-NULL, ie it's a smoother folks
      # IDNK if gam() returns S as a list if there's just one S...
      if( !is.list( Si)){ # ... so let's make sure it does
        Si <- list( Si)
      }
      Slist <- c( Slist, Si)
    }
  }
  
  # Separate rank for each penmat, if several...
  ranks <- unlist( FOR( G$smooth, .$rank))
  # but (?maybe?) the same subset of coefs for all
  # so, replicate 'first.para' if reqd
  repifreq <- function( smoo, name_of_thing){
    thing <- smoo[[ name_of_thing]]
    if( length( thing)==1){
      thing <- rep( thing, length( smoo$S))
    }
  }
    
  firsts <- unlist( FOR( G$smooth, repifreq( ., 'first.para')))
  # formo <- G$formula
  ncoef <- ncol( X)

  pifodds <- d1( G$y) / d2( G$y) - 1
  
  nlglk <- function( allpar){
    beta <- allpar[[1]]
    if( length( Slist)){
      log_lambda <- allpar[[2]]
      lambda <- exp( log_lambda)
    }

    ppn <- plogis( X %**% beta) # %**% to strip surplus dimension
    lprob <- log( 1+pifodds*ppn)
    lglk <- sum( lprob)

    REPORT( beta)
    REPORT( ppn)
    # Penalties
    for( i in seq_along( Slist)){
      m_i <- dim( Slist[[ i]])[1]
      irange <- seq( from=firsts[i], length=m_i)
      lglk <- lglk + 
          + 0.5 * ranks[ i] * log_lambda[i] + 
          - 0.5 * lambda[i] * beta[ irange] %**% Slist[[ i]] %**% beta[ irange]
    }

  return( -lglk)
  }
  
  if( length( Slist)){ # smooths
    allparz <- list( 
        beta= c( mean( pifodds>0), rep( 0, ncol( X)-1)),
        log_lambda= rep( 0.01, length( Slist)))

    obj <- MakeADFun( nlglk, allparz, random="beta")
  } else {
    # Fixed effects only: should match gasplit(), but use RTMB anyway
    allparz <- list( beta= c( mean( pifodds>0), rep( 0, ncol( X)-1)))
    obj <- MakeADFun( nlglk, allparz)
  }
  
  obj$fn( obj$par) # test here before nlminb()
  fitto <- nlminb( obj$par, obj$fn, obj$gr)
  obj <- c( 
      obj, 
      list( G=G),
      obj$report(), # beta and ppn
      fitto[ cq( convergence, message, evaluations)]
    )
return( obj)
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
function( 
    sim=NULL, 
    formulalala= LGLR ~ Y+Z-1, 
    use2= 's' %in% all.names( formula),
    ...
){
  if( is.null( sim)){
    sim <- make_fake_ppns( ...)
  }
  simpure <- sim
  
  e <- list2env( sim@truth, parent=.GlobalEnv)
  d1 <- function( x) dt( x - meanE, df=df_t)
  d2 <- function( x) dt( x + meanE, df=df_t) # meanW === -(meanE)
  environment( d1) <- environment( d2) <- e
  whichever_gasplit <- if( use2) gasplit2 else gasplit
  gg1 <- whichever_gasplit( formulalala, sim, d1=d1, d2=d2)
  
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
    formulalala,
    use2,
    gg1
  )
}

