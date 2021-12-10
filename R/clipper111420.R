
### Identifies one sided discoveries using Clipper
### score_exp: a vector or a matrix of measurements under experimental condition with rows being the features and the columns being the replicates
### score_back: a vector or a matrix of measurements under negative control or background with rows being the features and the columns being the replicates
### FDR: a vector of target FDR upper bounds
### ifuseknockoff: boolean.If set to TRUE, knockoffs will be constructed and used; if set to FALSE, no knockoffs will be constructed.
###                If not supplied (default), Clipper1sided will construct knockoffs if the numbers of replicates are different and skips knockoff construction if the number of replicates are the same
### nknockoff: the number of knockoffs to be constructed. only useful when ifuseknockoff is T or NULL.
### contrastScore_method: can be 'max'(default) or 'diff', only needed when knockoffs are constructed.
### importanceScore_method: can be 'max' or 'diff' (default)
### FDR_control_method: can be 'BH','BC','GZ'. Under knockoff construction, only 'GZ' will be used.
### ifpowerful: boolean. If set to TRUE (default), clipper1sided will use a heuristic approach when 'BC' or 'GZ' fail to identify any discovery.
### seed: random seed. Used for knockoff construction.

# remove(list = ls())
# FDR = c(0.0001,0.01, 0.05, 0.1)
# ifuseknockoff = NULL
# nknockoff = NULL
# contrastScore_method = 'max'
# importanceScore_method = 'diff'
# FDR_control_method = NULL
# ifpowerful = T
# seed = NULL
# temp1 = clipper1sided(score_exp,
#               score_back,
#               FDR ,
#               ifuseknockoff = T,
#               nknockoff = 2,
#               contrastScore_method = NULL,
#               importanceScore_method = 'diff',
#               FDR_control_method = NULL,
#               ifpowerful = T,
#               seed)
clipper1sided = function(score_exp,
                         score_back,
                         FDR = 0.05,
                         ifuseknockoff = NULL,
                         nknockoff = NULL,
                         contrastScore_method = NULL,
                         importanceScore_method = 'diff',
                         FDR_control_method = NULL,
                         ifpowerful = T,
                         seed = 12345){

  ### shift all measurements to be non-negative
  if(any(score_exp <0, na.rm = T) | any(score_back < 0,na.rm = T)){
    shift = min(min(score_exp[!is.na(score_exp)]),min(score_exp[!is.na(score_exp)]))
    score_exp = score_exp - shift
    score_back = score_back - shift
  }

  ### convert score_exp and score_back to matrices if they are numerical vectors
  if(is.null(dim(score_exp))){
    score_exp = matrix(score_exp, ncol = 1)
  }
  if(is.null(dim(score_back))){
    score_back = matrix(score_back, ncol = 1)
  }
  r1 = ncol(score_exp)
  r2 = ncol(score_back)

  ### check if score_exp and score_back have the same number of instances
  if(nrow(score_exp) != nrow(score_back)){
    stop('score_exp and score_back must have the same number of rows (features)')
  }

  ### default: use knockoffs when r1 neq r2.
  if( is.null(ifuseknockoff) ){
    if(r1 == r2){
      ifuseknockoff = F
    }else{
      ifuseknockoff = T
    }
  }

  ####################  use clipper_1sided_woknockoff  ####################
  if( !ifuseknockoff ){
    if( !r1 == r2 ){
      warnings( 'Caution: no knockoffs are constructed when the numbers of replicates are different; FDR control is not guaranteed' )
    }
    if( FDR_control_method == 'GZ'){
      warnings( 'FDR_control_method cannot be GZ when no knockoffs are constructed. Switching to BC.')
      FDR_control_method = 'BC'
    }
    re = suppressWarnings( clipper_1sided_woknockoff(score_exp = score_exp,
                                   score_back = score_back,
                                   r1 = r1,
                                   r2 = r2,
                                   FDR = FDR,
                                   importanceScore_method = importanceScore_method,
                                   FDR_control_method = FDR_control_method) )

    ###################### if FDR_control_method = BC or GZ but failt to identify any discovery at some FDR levels, switch to BH at those FDR levels.
    if( ifpowerful & FDR_control_method != 'BH'){
      FDR_nodisc = sapply(re$results, function(re_i){
        length(re_i$discovery) == 0
      })
      if( any(FDR_nodisc) ){
        warning(paste0('At FDR = ', paste0(FDR[FDR_nodisc], collapse = ', '), ', no discovery has been found using FDR control method ', FDR_control_method,' ; switching to BH...'))
        re_clipperbh = clipper_BH(contrastScore =  re$contrastScore, FDR = FDR[FDR_nodisc])
        re$results[FDR_nodisc] = re_clipperbh
      }
    }


  }


     #################### use clipper_1sided_wknockoff ####################
  if( ifuseknockoff ){

    if(r1 == 1 & r2 == 1){
      stop('Cannot generate knockoffs when both score_exp and score_back have one column. Please rerun clipper1sided by setting ifuseknockoff = F')
    }

    #### check if nknockoff is reasonable
    nknockoff_max = ifelse(r1 == r2, choose(r1+r2, r1)/2 - 1, choose(r1+r2, r1)-1)
    if(!is.null(nknockoff)){
      if(nknockoff > nknockoff_max | !is.integer(nknockoff) | nknockoff < 1){
        warnings('nknockoff must be a positive integer and must not exceed the maximum number of knockoff; using the exhaustive knockoffs instead.')
      }
      nknockoff = min(nknockoff, nknockoff_max)
    }else{
      warnings(paste0('nknockoff is not supplied; generate the maximum number of knockoffs: ', nknockoff_max))
      if (contrastScore_method == "diff"){
        nknockoff = nknockoff_max
      }else{
        nknockoff = 1
      }

    }

    re = suppressWarnings( clipper_1sided_wknockoff(score_exp = score_exp,
                                  score_back = score_back,
                                  r1 = r1,
                                  r2 = r2,
                                  FDR = FDR,
                                  importanceScore_method = importanceScore_method,
                                  contrastScore_method = contrastScore_method,
                                  nknockoff = nknockoff,
                                  nknockoff_max = nknockoff_max,
                                  seed = seed))


    ###################### if FDR_control_method = BC or GZ but failt to identify any discovery at some FDR levels, switch to BH at those FDR levels.
    if( ifpowerful & FDR_control_method != 'BH'){
      FDR_nodisc = sapply(re$results, function(re_i){
        length(re_i$discovery) == 0
      })
      if( any(FDR_nodisc) ){
        warning(paste0('At FDR = ', paste0(FDR[FDR_nodisc], collapse = ', '), ', no discovery has been found using FDR control method ',FDR_control_method,'; switching to BH...'))
        re_clipperbh = clipper_BH(contrastScore =  re$contrastScore, nknockoff = nknockoff, FDR = FDR[FDR_nodisc])
        re$results[FDR_nodisc] = re_clipperbh
      }
    }

  }


  return(re)
}


### Identifies two sided discoveries using Clipper
### score_exp: a vector or a matrix of measurements under experimental condition with rows being the features and the columns being the replicates
### score_back: a vector or a matrix of measurements under negative control or background with rows being the features and the columns being the replicates
### FDR: a vector of target FDR upper bounds
### nknockoff: the number of knockoffs to be constructed. Only positive integer is allowed.
### contrastScore_method: can be 'max'(default) or 'diff', only needed when knockoffs are constructed.
### importanceScore_method: can be 'max' or 'diff' (default)
### FDR_control_method: can be 'BH','BC','GZ'. Under knockoff construction, only 'GZ' will be used.
### ifpowerful: boolean. If set to TRUE (default), clipper1sided will use a heuristic approach when 'BC' or 'GZ' fail to identify any discovery.
### seed: random seed. Used for knockoff construction.
# score_exp = clipper_input$exp
# score_back = clipper_input$back
# FDR = 0.05
# importanceScore_method = 'diff'
# contrastScore_method = 'max'
# nknockoff = NULL
# ifpowerful = F
# FDR_control_method = 'GZ'
# seed = 12345
clipper2sided = function(score_exp,
                         score_back,
                         FDR = 0.05,
                         nknockoff = NULL,
                         contrastScore_method = 'max',
                         importanceScore_method = 'diff',
                         FDR_control_method = 'GZ',
                         ifpowerful = T,
                         seed = 12345){

  ### convert score_exp and score_back to matrices if they are numerical vectors
  if(is.null(dim(score_exp))){
    score_exp = matrix(score_exp, ncol = 1)
  }
  if(is.null(dim(score_back))){
    score_back = matrix(score_back, ncol = 1)
  }
  r1 = ncol(score_exp)
  r2 = ncol(score_back)

  if( r1 == 1 & r2 == 1){
    stop( 'clipper is not yet able to perform two sided identification when either condition has one replicate' )
  }

  #### check if nknockoff is reasonable
  nknockoff_max =  min(ifelse(r1 == r2, choose(r1+r2, r1)/2 - 1, choose(r1+r2, r1)-1), 200)
  if(!is.null(nknockoff)){
    if(nknockoff > nknockoff_max | !is.integer(nknockoff) | nknockoff < 1){
      warnings('nknockoff must be a positive integer and must not exceed the maximum number of knockoff; using the maximal number of knockoffs instead.')
    }
    nknockoff = min(nknockoff, nknockoff_max)
  }else{
    if (contrastScore_method == "diff"){
      nknockoff = min(nknockoff_max, 50)
    }else{
      nknockoff = 1
    }
  }

  ######################
  knockoffidx = generate_knockoffidx( r1 = r1, r2 = r2,
                                      nknockoff = nknockoff,
                                      nknockoff_max = nknockoff_max,
                                      seed = seed)
  kappatau_ls = compute_taukappa(score_exp = score_exp,
                                 score_back = score_back,
                                 r1, r2,
                                 if2sided = T,
                                 knockoffidx,
                                 importanceScore_method,
                                 contrastScore_method)

  re = clipper_GZ(tau = kappatau_ls$tau, kappa = kappatau_ls$kappa, nknockoff = nknockoff, FDR = FDR)
  re  = list(knockoffidx = knockoffidx,
             importanceScore_method = importanceScore_method,
             importanceScore = kappatau_ls$importanceScore,
             contrastScore_method = contrastScore_method,
             contrastScore = (2*kappatau_ls$kappa  - 1) * abs(kappatau_ls$tau),
             results = re)

  ###################### if FDR_control_method = BC or GZ but failt to identify any discovery at some FDR levels, switch to BH at those FDR levels.
  if( ifpowerful & FDR_control_method != 'BH'){
    FDR_nodisc = sapply(re$results, function(re_i){
      length(re_i$discovery) == 0
    })
    if( any(FDR_nodisc & contrastScore_method == 'max') ){
      warning(paste0('At FDR = ', paste0(FDR[FDR_nodisc], collapse = ', '), ', no discovery has been found using FDR control method ', FDR_control_method,' ; switching to BH...'))
      re_clipperbh = clipper_BH(contrastScore =  re$contrastScore, nknockoff = nknockoff, FDR = FDR[FDR_nodisc])
      re$results[FDR_nodisc] = re_clipperbh
    }
  }


  return(re)

}

clipper_1sided_woknockoff = function(score_exp,
                   score_back,
                   r1, r2,
                   FDR = 0.05,
                   aggregation_method = 'mean',
                   importanceScore_method,
                   FDR_control_method){


  ### aggregate multiple replicates into single replicate
  if(r1 > 1){
    score_exp = aggregate_clipper(score = score_exp, aggregation_method = aggregation_method)
  }
  if(r2 > 1){
    score_back = aggregate_clipper(score = score_back, aggregation_method = aggregation_method)
  }

  contrastscore = compute_importanceScore_wsinglerep(score_exp = score_exp,
                                        score_back = score_back,
                                        importanceScore_method = importanceScore_method)

  if(FDR_control_method == 'BC'){
    re = suppressWarnings( clipper_BC(contrastScore = contrastscore,
                                      FDR = FDR) )

  }

  if(FDR_control_method == 'BH'){
    re = clipper_BH(contrastScore = contrastscore, FDR = FDR)
  }

  re  = list(importanceScore = contrastscore,
             importanceScore_method = importanceScore_method,
             contrastScore = contrastscore,
             contrastScore_method = importanceScore_method,
             results = re)

  return(re)
}


clipper_1sided_wknockoff = function(score_exp,
                                    score_back,
                                    r1, r2,
                                    FDR = 0.05,
                                    importanceScore_method,
                                    contrastScore_method,
                                    nknockoff,
                                    nknockoff_max,
                                    seed){

  knockoffidx = generate_knockoffidx( r1, r2, nknockoff,nknockoff_max, seed = seed)
  kappatau_ls = compute_taukappa(score_exp, score_back, r1, r2, if2sided = F,
                              knockoffidx, importanceScore_method, contrastScore_method)

  re = clipper_GZ(tau = kappatau_ls$tau, kappa = kappatau_ls$kappa, nknockoff, FDR)
  re  = list(knockoffidx = knockoffidx,
              importanceScore_method = importanceScore_method,
              importanceScore = kappatau_ls$importanceScore,
              contrastScore_method = contrastScore_method,
             contrastScore = (2*kappatau_ls$kappa  - 1) * abs(kappatau_ls$tau),
             results = re)
  return(re)
}


aggregate_clipper = function(score, aggregation_method){

  if(aggregation_method == 'mean'){
    score_single = apply(score, 1, function(x){
      mean(x, na.rm = T)
    })
  }
  if(aggregation_method == 'median'){
    score_single = apply(score, 1, function(x){
      median(x, na.rm = T)
    })
  }
  return(score_single)

}

compute_importanceScore_wsinglerep = function(score_exp, score_back, importanceScore_method){
  if(importanceScore_method == 'diff'){
    contrastScore = score_exp - score_back
  }
  if(importanceScore_method == 'max'){
    contrastScore = pmax(score_exp, score_back)*sign(score_exp - score_back)
  }
  return(as.vector(contrastScore))
}

clipper_BC = function(contrastScore, FDR){

  contrastScore[is.na(contrastScore)] = 0 # impute missing contrast scores with 0
  c_abs = abs(contrastScore[contrastScore != 0])
  c_abs  = sort(unique(c_abs))

  i = 1
  emp_fdp = rep(NA, length(c_abs))
  emp_fdp[1] = 1
  while(i <= length(c_abs)){
    # print(i)
    t = c_abs[i]
    emp_fdp[i] = min((1 + sum(contrastScore <= -t))/ sum(contrastScore >= t),1)
    if (i >=2){emp_fdp[i] = min(emp_fdp[i], emp_fdp[i-1])}
    i = i + 1
  }

  c_abs = c_abs[!is.na(emp_fdp)]
  emp_fdp = emp_fdp[!is.na(emp_fdp)]
  q <- emp_fdp[match(contrastScore, c_abs)]
  q[which(is.na(q))] = 1

  re = lapply(FDR, function(FDR_i){
    thre = c_abs[min(which(emp_fdp <= FDR_i))]
    re_i = list(FDR = FDR_i,
                FDR_control = 'BC',
                thre = thre,
                q = q,
                discovery = which(contrastScore >= thre))
    return(re_i)
  })
  return(re)
}

clipper_BH = function(contrastScore, nknockoff = NULL, FDR){


  if( is.list(contrastScore) ){
    n = length(contrastScore[[1]])
    kappa = contrastScore$kappa
    tau = contrastScore$tau
    idx_na = is.na(tau) | is.na(kappa)
    tau = tau[!idx_na]
    kappa = kappa[!idx_na]
    pval = sapply(1:n, function(i){
      sum( !kappa & tau >= tau[i] , na.rm = T) / sum( !kappa, na.rm = T ) * nknockoff / (nknockoff+1)
    })

  }else{
    n = length(contrastScore)
    idx_na = is.na(contrastScore)
    contrastScore_nomiss = contrastScore[ !idx_na ] # exclude features with NA contrast scores

    cs_negative <- contrastScore_nomiss[ contrastScore_nomiss < 0 ]
    cs_null <- c( cs_negative, -cs_negative )
    pval <- sapply(contrastScore_nomiss, function(x){
      mean( x <= cs_null )
    })
  }


  qvalue <- p.adjust(pval, method = 'BH')


  re = lapply(FDR, function(FDR_i){
    re_i = list(FDR = FDR_i,
                FDR_control = 'BH',
                discovery = (1:n)[!idx_na][which(qvalue <= FDR_i)],
                q = qvalue)
    return(re_i)
  })

  return(re)
}

generate_knockoffidx = function( r1, r2, nknockoff,nknockoff_max, seed){
  set.seed(seed)
  if(nknockoff_max == 200){
    ### randomly permute nknockoff times
    knockoffidx = vector('list', length= nknockoff)
    i_knockoff = 1
    while(i_knockoff <= nknockoff){
      temp = sample(r1 + r2, r1, replace = F)
      if(any(!temp %in% 1:r1) & any(temp %in% 1:r1)){
        knockoffidx[[i_knockoff]] = temp
        i_knockoff = i_knockoff + 1
      }else{
        next
      }
    }


  }else{

    ### exhaust all possible combinations
    combination_all <- combn(r1 + r2, r1, simplify = F)
    if(r1 == r2){
      combination_all <- combination_all[1:(length(combination_all)*0.5)][-1]
    }else{
      combination_all <- combination_all[-1]
    }


    knockoffidx = combination_all[sample(1:nknockoff_max, size = nknockoff, replace = F)]

  }



  return(knockoffidx)
}

compute_taukappa = function(score_exp, score_back, r1, r2, if2sided,
                            knockoffidx, importanceScore_method, contrastScore_method){
  perm_idx = c(list(1:r1), knockoffidx)
  score_tot = cbind(score_exp, score_back)

  imp_ls = sapply(perm_idx, function(x){
    se = score_tot[, x,drop = F]
    sb = score_tot[, setdiff(1:(r1+r2), x),drop = F]
    se = aggregate_clipper(se, aggregation_method = 'mean')
    sb = aggregate_clipper(sb, aggregation_method = 'mean')
    imp = compute_importanceScore_wsinglerep(se,sb,importanceScore_method = importanceScore_method)
    if(if2sided){
      imp = abs(imp)
    }
    return(imp)
  })

  kappatau_ls = apply(imp_ls, 1, function(x){
    # idx_max = which(x == max(x))
    # kappa = sample(idx_max, 1)
    kappa = !any(x[-1] == max(x)) ### kappa should be true iff the maximum occurs only at the first index
    if(length(kappa) == 0){
      kappa = NA
    }
    x_sorted = sort(x, decreasing = T)
    if(contrastScore_method == 'diff'){
      tau = x_sorted[1] - x_sorted[2]
    }
    if(contrastScore_method == 'max'){
      tau = x_sorted[1]
    }
    return(list(kappa = kappa, tau = tau))
  })

  re = list(importanceScore = imp_ls,
    kappa = sapply(kappatau_ls,'[[','kappa'),
              tau = sapply(kappatau_ls,'[[','tau'))

  return(re)

}

# tau = kappatau_ls$tau
# kappa = kappatau_ls$kappa
clipper_GZ = function(tau, kappa, nknockoff, FDR){
  n = length(tau)
  contrastScore = (2*kappa  - 1) * abs(tau)
  contrastScore[is.na(contrastScore)] = 0 # impute missing contrast scores with 0
  c_abs = abs(contrastScore[contrastScore != 0])
  c_abs  = sort(unique(c_abs))

  i = 1
  emp_fdp = rep(NA, length(c_abs))
  emp_fdp[1] = 1
  while(i <= length(c_abs)){
    # print(i)
    t = c_abs[i]
    emp_fdp[i] = min((1/nknockoff + 1/nknockoff * sum(contrastScore <= -t))/ sum(contrastScore >= t),1)
    if (i >=2){emp_fdp[i] = min(emp_fdp[i], emp_fdp[i-1])}
    i = i + 1
  }

  c_abs = c_abs[!is.na(emp_fdp)]
  emp_fdp = emp_fdp[!is.na(emp_fdp)]
  q <- emp_fdp[match(contrastScore, c_abs)]
  q[which(is.na(q))] = 1


  re = lapply(FDR, function(FDR_i){
    thre = c_abs[min(which(emp_fdp <= FDR_i))]
    re_i = list(FDR = FDR_i,
                FDR_control = 'BC',
                thre = thre,
                discovery = which(contrastScore >= thre),
                q = q)
    return(re_i)
  })
  return(re)
}

# ######################### testing #########################
# ######################### single replicate with missing values #########################
# ## clipper deals with missing contrast scores differently under the BC and BH
# ## Under BC, missing contrast scores are imputed with 0
# ## Under BH, missing contrast scores are excluded from the construction of null distribution.
# ## In either case, features with missing contrast scores will never be called.
# set.seed(1)
# n = 1000
# trueidx = 1:(0.1*n)
# score_exp = rnorm(n)
# score_exp[1:(0.1*n)] = rnorm(0.1*n, mean = 5)
# score_exp[sample(1:(0.1*n), 10)] = NA
# score_back = rnorm(n)
# score_back[sample(1:(0.1*n), 10)] = NA
# FDR = c(0.001,0.01, 0.05, 0.1)
# aggregation_method = 'mean';
# contrastScore_method = 'diff';
# FDR_control_method = 'BC';
# ifpowerful = T
# re_clipper = clipper(score_exp,
#                      score_back,
#                      FDR = FDR,
#                      aggregation_method = aggregation_method,
#                      contrastScore_method = contrastScore_method,
#                      FDR_control_method = FDR_control_method,
#                      ifpowerful = ifpowerful )
#
#
# ######################### multiple replicates with missing values #########################
# set.seed(1)
# source('/Users/chenyiling/Dropbox/clipper/source042520/2sided/clipper2.R')
# n = 10000
# trueidx = 1:(0.1*n)
# r1 = 3
# r2 = 2
# prop = 0.05 ### proportion of NA
# mean_back = rnbinom(n, size = 20, prob = 0.5)
# mean_exp = mean_back
# mean_exp[1:length(trueidx)] = rpois(length(trueidx), lambda = 40)
#
# score_exp <- t(sapply(mean_exp, function(mean_exp_i){
#   rpois(r1, lambda = mean_exp_i)
# }))
# score_back <- t(sapply(mean_back, function(mean_back_i){
#   rpois(r2, lambda = mean_back_i)
# }))
#
# ### no knockoff needed
# temp1 = clipper1sided(score_exp,
#                       score_back,
#                       FDR ,
#                       ifuseknockoff = T,
#                       nknockoff = NULL,
#                       contrastScore_method = 'max',
#                       importanceScore_method = 'diff',
#                       FDR_control_method = 'BH',
#                       ifpowerful = T)
# temp2 = clipper(score_exp,
#                    score_back,
#                    FDR ,
#                    aggregation_method = 'mean',
#                    contrast_score_method = 'max',
#                    FDR_control_method = 'BC',
#                    ifpowerful = T,
#                    two.sides = T,
#                    n.knockoff = 2)
#
# all(unlist(temp1$results) == unlist(temp2$results))
#
# ##### with missing value
#
# if (prop > 0) {
#   score_back[sample(1:(r2*n) , (r2*n)*prop)] <- NA
#   score_exp[sample(1:(r1*n), (r1*n)*prop)] <- NA
# }
#
# FDR = c(0.001,0.01, 0.05, 0.1)
# score_exp = score_exp
# score_back = score_back
# FDR = FDR
# ifuseknockoff = T
# nknockoff = 1
# contrastScore_method = 'max'
# importanceScore_method = 'diff'
# FDR_control_method = 'BC'
# ifpowerful = T
#
# temp1 = clipper1sided(score_exp = score_exp,
#                       score_back = score_back,
#                       FDR = FDR,
#                       ifuseknockoff = ifuseknockoff,
#                       nknockoff = nknockoff,
#                       contrastScore_method = contrastScore_method,
#                       importanceScore_method = importanceScore_method,
#                       FDR_control_method = FDR_control_method,
#                       ifpowerful = ifpowerful)
#

# aggregation_method = 'mean';
# contrastScore_method = 'diff';
# FDR_control_method = 'BC';
# re_clipper = clipper(score_exp,
#                      score_back, FDR = FDR,
#                      aggregation_method = aggregation_method,
#                      contrastScore_method = contrastScore_method,
#                      FDR_control_method = FDR_control_method)
# lapply(re_clipper$results, function(re_i){
#   FDR = re_i$FDR
#   c(FDR, sum(!re_i$discovery %in% trueidx)/max(length(re_i$discovery),1))
# })


