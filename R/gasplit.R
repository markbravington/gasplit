# This is package gasplit 

"do_predict_from_previous" <-
function( nlocal=sys.parent()) mlocal({
  # I think we have to actually "fit" the gam, so that we can use predict()
  # Result is garbage here!
  # Maybe this could be sped up by using existing estimates as startvals
  
  PFP <- predict_from_previous # brevity
  G <- gam( G=PFP$G)
  # predict.gam() requires data to be either bona fide dataframe, or missing
  X <- if( is.null( data)) predict( G, type='lpmatrix') else
      predict( G, newdata=data, type='lpmatrix')
  linpred <- X %**% PFP$beta
  ppn <- plogis( linpred)

  if( dbeta){
    dbeta <- dlogis( linpred) * X
    # check_dbeta <- numvbderiv( function( b) plogis( X %**% b), PFP$beta)
    ppn <- cbind( ppn, dbeta)
  }
})


"example_prep4missy" <-
function( data){
## Make an extended version of 'data' where missing-Sex is replaced by multimps (multiple imputation). Use actual Sex if known. Otherwise, add extra rows with alternative sex and use PrFem_yl

  unk <- is.na( data$Sex)
  data <- within( data, {
    MULTIMP <- seq_along( Sex) # label
    PROBIMP <- ifelse( !is.na( Sex), 1, PrFem_yl)    
    Sex[ is.na( Sex)] <- 'F' # will be 'M' in the extra alternatives
  })
  
  extroid <- within( data[ unk,], {
    MULTIMP <- which( unk)
    Sex[] <- 'M'
    PROBIMP <- 1-PrFem_yl
  })
  
return( rbind( data, extroid))
}


"gasplit" <-
function( 
    formula, 
    data, 
    d1, d2, 
    start=0.001, 
    predict_from_previous=NULL,
    dbeta= FALSE
){
  if( !is.null( predict_from_previous)){
    do_predict_from_previous()
return( ppn)
  }

  # Else (normally) we want to fit. Set up to fit gam, but don't do it
  G <- gam( formula=formula, data=data, fit=FALSE) # stuff needed for fit
  X <- G$X
  pifodds <- d1( G$y) / d2( G$y) - 1
  
  ppn <- 0*pifodds - 1
  
  nlglk <- function( pars){
    ppn <<- plogis( X %**% pars)
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

  fitto <- optim( start, nlglk, method='BFGS', hessian=TRUE,
      control=list( trace=5))
      
  beta <- fitto$par    
  V_beta <- solve( fitto$hessian)
  SE_beta <- sqrt( diag( V_beta))

  if( dbeta){
    linpred <- X %**% fitto$par
    dbeta <- dlogis( linpred) * X
    ppn <- cbind( ppn, dbeta)
  }

  # Keep distros, eg for predictions/posteriors
  
return( c( returnList( 
    beta, SE_beta, V_beta,
    ppn, G, d1, d2),
    fitto[ cq( convergence, message, evaluations)]
  ))
}


"gasplit_missy" <-
function( 
  formula, 
  data, 
  d1, d2,
  link_field= 'MULTIMP',
  prob_field= 'PROBIMP',
  predict_from_previous=NULL,
  dbeta= FALSE,
  ... # for gam()
){
  if( !all( hasName( data, c( link_field, prob_field)))){
    warning( sprintf( 
      "Fields '%s' and/or '%s' not found; calling 'gasplit2' for ya instead",
      link_field, prob_field
    ))
    
    mc <- match.call( expand.dots=TRUE)
    mc[[1]] <- quote( 'gasplit2')
return( eval.parent( mc))
  }

  if( !is.null( predict_from_previous)){
    do_predict_from_previous()
    
    # DEAL WITH AGGREGATION...
return( ppn)    
  }

  # Prepare for FIML, integrating over all "imputations" of each case:
  if( any( c( link_field, prob_field)  %in% all.names( formula))) {
stop( "WHAAAAT are you thinking??? Link & prob fields don't belong in formula!")
  }
  linkid <- data[[ link_field]]
  linkid <- match( linkid, unique( linkid)) # integer: first will be 1
  probimp <- data[[ prob_field]]
stopifnot(
    all( probimp >= 0),
    all( probimp <= 1)
  )
  
  # Normalize probimp
  sump <- tapply( probimp, INDEX= list( linkid), sum)
  if( any( abs( log( sump)) > 0.01)){
    warning( "Multimp probs seem a bit... under-normalized. FT4U, but ...")
  }
  mm <- match( linkid, as.integer( names( sump)), 0)
  probimp <- probimp / sump[ mm]
  
  multab <- tabulate( linkid) # safe coz ints starting at 1
  primary <- which( linkid == seq_along( linkid)) # only these responses are used
  no_multi <- which( multab==1) # simples
  multi <- primary %except% no_multi # first multimp for each m-case

  # Organize data for efficient calcs: loop over number-of-multimps
  xlinkid <- linkid
  xlinkid[ primary] <- (-1L) # these are done
  max_n_multi <- max( multab)
  has_at_least <- nth_multi <- vector( 'list', max_n_multi)
  # has_at_least[[1]] <- primary
  for( im in 2 %upto% max_n_multi){
    this_or_more <- which( multab >= im)
    has_at_least[[ im]] <- this_or_more
    nth_instance <- match( this_or_more, xlinkid)
    nth_multi[[ im]] <- nth_instance
    xlinkid[ nth_instance] <- (-1L) # won't be found again
  }

  # Multimpized version   
  nlglk <- function( allpar){
    beta <- allpar[[1]]
    if( length( Slist)){
      log_lambda <- allpar[[2]]
      lambda <- exp( log_lambda)
    }

    ppn <- plogis( X %**% beta) # %**% to strip surplus dimension
    prob <- 1+pifodds*ppn # these are UNSCALED probs, already with const mult of d2() taken out. OK for logging.
    
    # lprob <- log( 1+pifodds*ppn)
    # lglk <- sum( lprob)
    
    # Non-multi ones are thereby done (but can safely be Xed by PROBIMP of 1 :). 
    # Accumulate multi ones into actual cases (only doing the first imps here)
    prob <- prob * probimp
    ppn <- ppn * probimp
    
    for( im in 2 %upto% max_n_multi){
      # Update only those cases with at least "im" multimps
      prob[ has_at_least[[ im]] ] <- prob[ has_at_least[[ im]] ] + prob[ nth_multi[[ im]] ]
      ppn[ has_at_least[[ im]] ] <- ppn[ has_at_least[[ im]] ] + ppn[ nth_multi[[ im]] ]
    }
    
    lglk <- sum( log( prob[ primary]))
    ppn <- ppn[ primary] # expected ppn just for each overall case (not imps); consistent with gasplit2

    REPORT( beta)
    REPORT( ppn)
    
    # Penalties
    for( i in seq_along( Slist)){
      m_i <- dim( Slist[[ i]])[1]
      irange <- seq( from=firsts[i], length=m_i)
      spi <- smoopar_map[ i]
      lglk <- lglk + 
          + 0.5 * ranks[ i] * log_lambda[ spi] + 
          - 0.5 * lambda[ spi] * beta[ irange] %**% Slist[[ i]] %**% beta[ irange]
    }

  return( -lglk)
  }

  guts_gasplit2() # all the GAM setup and fitting. Same as for 'gasplit2', but 'nlglk' is slightly different
return( retlist)
}


"gasplit2" <-
function( 
  formula, 
  data, 
  d1, d2, 
  predict_from_previous=NULL,
  dbeta= FALSE,
  ... # for gam()
){
## Allow REs (eg to do splines) and mgcv-style stuff
  if( !is.null( predict_from_previous)){
    do_predict_from_previous()
return( ppn)    
  }

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
      spi <- smoopar_map[ i]
      lglk <- lglk + 
          + 0.5 * ranks[ i] * log_lambda[ spi] + 
          - 0.5 * lambda[ spi] * beta[ irange] %**% Slist[[ i]] %**% beta[ irange]
    }

  return( -lglk)
  }
  
  guts_gasplit2() # all the GAM setup and fitting. Separate function coz also used for 'gasplit_multimp'--- only 'nlglk' is different
return( retlist)
}


"guts_gasplit2" <-
function( nlocal=sys.parent()) mlocal({
  # Else (ie normally), we fit
  # Stuff needed for fit, but don't actually fit:
  G <- gam( formula=formula, data=data, fit=FALSE, ...) 
  Gfake <- gam( formula=formula, data=data, G=G) # actually "fit" it, mainly for...
  coef_names <- names( coef( Gfake))
  nco <- length( coef_names)

  # Tidy up the ingredients. I am not too sure about this, when multiple smooths are used eg with "id=1" in "by"..
  
  X <- G$X
  y <- G$y
  n_smoopar <- length( G$sp) # several smooths may share same smoopar...
  smoopar_names <- names( G$sp)
  smoopar_map <- integer()
  
  Slengths <- unlist( FOR( G$smooth, length( .$S)))
  Slist <- list()
  
  for( i in seq_along( G$smooth)){
    Si <- G$smooth[[i]]$S
    if( length( Si)){ # non-NULL, ie it's a smoother folks
      # IDNK if gam() returns S as a list if there's just one S...
      if( !is.list( Si)){ # ... so let's make sure it does
        Si <- list( Si)
      }
      
      # Which smoopar(s) to use here? $first.sp and $last.sp are somehow related to expanded set of smoopars, not to the underlying ones that get optimzed over (eg if 'id' is used). I think matching on names is safest...
      
      spinds <- match( names( G$smooth[[i]]$sp), smoopar_names, 0)
      if( length( spinds) != length( Si)){
stop( sprintf(
        "Not sure how to handle this: more smoopars than smoomats for %i-th smooth term", i))      
      }
      smoopar_map <- c( smoopar_map, spinds)
      Slist <- c( Slist, Si)
    }
  }
  
stopifnot( 
    length( smoopar_map) == length( Slist),
    all( smoopar_map > 0),
    all( smoopar_map <= n_smoopar)
  )
  
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
  
  # "Smoothers" (ie models with random effects) need slightly different treatment to fixed-effect models, in order to get eg Hessian. Fixed-eff models can of course also be fitted with gasplit(), but should be consistent here.

  if( length( Slist)){ # smooths. Two-stage fit (see below):
    # NB: I think mgcv might allow constraining smoopars (ie several to take same value). NYI here; would need "map" 
  
    allparz <- list( 
        beta= c( mean( pifodds>0), rep( 0, ncol( X)-1)),
        log_lambda= rep( 0.01, n_smoopar))

    obj_outer <- RTMB::MakeADFun( nlglk, allparz, random="beta")
    outer_pars <- obj_outer$par
    obj_outer$fn( outer_pars) # test here before nlminb()
    
    # Can't reliably use nlminb coz often 1D (1 smoopar)! Only the RE var is "outer"
    opto <- with( obj_outer, optim( par, fn, gr, method='BFGS'))
    opto$evaluations <- opto$counts # similar to nlminb
    outer_pars <- opto$par

    # I'm happy with variances conditional on estimated outer pars (ie just log_lambda). RTMB does not yield that up easily, but we can re-fit with log_lambda fixed, and no "random effects":

    # Might as well start at the inner optimum...
    allparz$beta[] <- obj_outer$env$last.par.best[ 1:nco] 
    
    # Fix the "outer" par(s)
    allparz$log_lambda <- outer_pars
    fix_log_lambda <- list( log_lambda= factor( NA+allparz$log_lambda))
    obj <- RTMB::MakeADFun( nlglk, allparz, map=fix_log_lambda)
    fitto <- with( obj, nlminb( par, fn, gr))
  } else { # FIXED EFFECTS ONLY
    # Should match gasplit(), but use RTMB anyway
    allparz <- list( beta= c( mean( pifodds>0), rep( 0, ncol( X)-1)))
    obj <- RTMB::MakeADFun( nlglk, allparz)
    
    obj$fn( obj$par) # test here before nlminb()    
    fitto <- nlminb( obj$par, obj$fn, obj$gr)    
    outer_pars <- numeric(0)
  }

  rep <- obj$report()
  names( rep$beta) <- coef_names
  rep$ppn <- c( rep$ppn) # strip dim that makes it a 1D array! Grrr...
  
  H <- obj$he( fitto$par)
  dimnames( H) <- list( coef_names, coef_names)
  rep$V_beta <- solve( H)
  rep$SE_beta <- sqrt( diag( rep$V_beta))

  # Keep distros, eg for predictions/posteriors

  retlist <- c( 
      rep, # beta, SE, V, ppn
      returnList( G, d1, d2, outer_pars, obj),      
      fitto[ cq( convergence, message, evaluations)]
    )

  if( dbeta){
    #  ppn <- plogis( G$X %**% beta)
    linpred <- G$X %**% rep$beta
    dbeta <- dlogis( linpred) * G$X
    # check_dbeta <- numvbderiv( function( b) plogis( G$X %**% b), rep$beta)

    retlist$ppn <- cbind( retlist$ppn, dbeta)
  }
})


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
  if( !requireNamespace( 'offarray')){
stop( "This demo function requires 'offarray' package")
  }
  
  # Faaaaaaaaaaaaark the Craniacs are on the loose
  offarray <- offarray::offarray
  autoloop <- offarray::autoloop

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
  dfE$DISCOSTAT <- rt( nrow( dfE), df=df_t) + meanE
  
  meanW <- (-meanE)  # opposite mean
  dfW <- df[ rep( seq_len( nstrat), c( nobsW)),]
  dfW$whicho <- 'W'
  dfW$DISCOSTAT <- rt( nrow( dfW), df=df_t) + meanW
  
  dfall <- rbind( dfE, dfW)
  rownames( dfall) <- NULL # they are just annoying
  
  dfall@truth <- returnList( yeff, zeff, yzeff, df_t, 
      meanE, ppnE, samp_ppnE, prange, seed)
return( dfall)
}


"posterior" <-
function( 
  object, 
  newdata= NULL, 
  dbeta=FALSE
){
  d1 <- object$G$d1
  d2 <- object$G$d2  

  pE <- NULL # make it below  
  if( is.null( newdata)){
    y <- object$G$y
    if( !dbeta){
      pE <- object$ppn # fitted
    }
  } else {
    # Response variable (perhaps transformed, as per formula)
    y <- eval( object$G$formula[[2]], newdata)
  }
  
  LR21 <- d2( y) / d1( y)
  
  if( is.null( pE)){
    # Use gasplit2, even if original was gasplit()
    pE <- gasplit2( data=newdata, predict_from_previous=object, 
        dbeta=dbeta)
  }
  
  pEv <- if( dbeta) pE[,1] else pE
  posterior <- 1/( 1+LR21*(1/pEv-1))
  
  if( dbeta){
    # D( quote( 1/( 1+LR21*(1/pE-1))), 'pE')
    DpE <-  LR21 * (1/pEv^2)/(1 + LR21 * (1/pEv - 1))^2
    dpost_dbeta <- DpE * pE # pE is "really" cbind( pE, dpE/dbeta)
    dpost_dbeta[,1] <- posterior
    posterior <- dpost_dbeta
  }
  
return( posterior)
}


"test_gasplit" <-
function( 
    sim=NULL, 
    formulalala= DISCOSTAT ~ Y+Z-1, 
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

