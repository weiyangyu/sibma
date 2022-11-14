## SIB
##
## Find the global best and local best particles from a particle pool.
##
## @param X a list of particles.
## @param P_w a nested list. The largest list is composed of smaller lists in which
## each element is a matrix describing the orthogonal projection matrix onto
## the corresponding stratum variance.
## @param q_GB an integer describing how many columns of a particle
## should be mixed with the corresponding columns of \code{GB}.
## @param q_LB an integer describing how many columns of a particle
## should be mixed with the corresponding columns of \code{LB}.
## @param q_new an integer describing how many columns of a particle
## should be mixed with the corresponding columns of a random
## particle created by \code{sibma::create_particle()}.
## @param all_two_level logical. If all factors are two levels,
## \code{all_two_level} should be \code{TRUE}; otherwise it should be
## \code{FALSE}.
## @param factor_level a list. Each element is a numeric vector specifying
## levels of a factor.
## @param t an integer indicating iteration numbers.
## @param total_unit an integer representing the number of total run size in
## a full factorial design.
## @param multiply_len a numeric vector. Each element is used to modify
## the value of each column in a model matrix.
## @param incidence_matrix a list. Each element is a matrix specifying an
## incidence matrix.
## @param treatment_effect a list. Each element is a numeric vector specifying
## which factors should have the same treatment effect in terms of a certain
## incidence matrix. The order of elements in \code{treatment_effect} should be
## corresponding to that in \code{incidence_matrix}.
SIB <- function(X, P_w, q_GB, q_LB, q_new, all_two_level, factor_level, t,
                total_unit, multiply_len, incidence_matrix, treatment_effect){
  unit <- nrow(X[[1]])
  factor_number <- ncol(X[[1]])
  #######################

  # initial step
  initial_step <- initial(particle = X, all_two_level = all_two_level,
                          P_w = P_w, factor_level = factor_level,
                          total_unit = total_unit, unit = unit,
                          multiply_len = multiply_len)
  # mix and move step
  mixwGB <- lapply(initial_step$X, mix_with_GB, all_two_level = all_two_level,
                  GB = initial_step$GB,
                  P_w = P_w,
                  q_GB = q_GB, factor_level = factor_level,
                  total_unit = total_unit, unit = unit)

  X <- lapply(1:length(X), function(i){
    return(move(X = X[[i]], mixwGB = mixwGB[[i]], t_1 = TRUE, P_w = P_w,
                q_new = q_new, all_two_level = all_two_level,
                factor_level = factor_level, total_unit = total_unit,
                unit = unit, multiply_len = multiply_len,
                incidence_matrix = incidence_matrix,
                treatment_effect = treatment_effect))
  })
  # evaluate LB and GB
  LB <- lapply(1:length(X), function(i){
    return(evaluation(X = X[[i]], LB_0 = initial_step$LB[[i]], evaluate_GB = FALSE,
                      factor_number = factor_number, P_w = P_w,
                      all_two_level = all_two_level,
                      factor_level = factor_level,
                      total_unit = total_unit, unit = unit,
                      multiply_len = multiply_len))
  })

  eva <- evaluation(LB = LB, evaluate_GB = TRUE, factor_number = factor_number,
                     P_w = P_w, all_two_level = all_two_level,
                     factor_level = factor_level, total_unit = total_unit,
                     unit = unit, multiply_len = multiply_len)
  GB <- eva$GB
  ####################
  # t = 2 and upper
  for(a in 1:(t-1)){
    # mix and move step
    mixwGB <- lapply(X, mix_with_GB, all_two_level = all_two_level,
                     GB = GB, P_w = P_w,
                     q_GB = q_GB, factor_level = factor_level,
                     total_unit = total_unit, unit = unit)
    mixwLB <- lapply(1:length(X), function(i){
      return(mix_with_LB(X = X[[i]], all_two_level = all_two_level,
                         LB = LB[[i]], P_w = P_w, q_LB = q_LB,
                         factor_level = factor_level, total_unit = total_unit,
                         unit = unit))
    })
    X <- lapply(1:length(X), function(i){
      return(move(X = X[[i]], mixwGB = mixwGB[[i]], mixwLB = mixwLB[[i]],
                  t_1 = FALSE, P_w = P_w,
                  q_new = q_new, all_two_level = all_two_level,
                  factor_level = factor_level, total_unit = total_unit,
                  unit = unit, multiply_len = multiply_len,
                  incidence_matrix = incidence_matrix,
                  treatment_effect = treatment_effect))
    })
    # evaluate LB and GB
    LB <- lapply(1:length(X), function(i){
      return(evaluation(X = X[[i]], LB_0 = LB[[i]], evaluate_GB = FALSE,
                        factor_number = factor_number, P_w = P_w,
                        all_two_level = all_two_level,
                        factor_level = factor_level,
                        total_unit = total_unit, unit = unit,
                        multiply_len = multiply_len))
    })
    eva <- evaluation(LB = LB, evaluate_GB = TRUE, factor_number = factor_number,
                      P_w = P_w, all_two_level = all_two_level,
                      factor_level = factor_level, total_unit = total_unit,
                      unit = unit, multiply_len = multiply_len)
    GB <- eva$GB
  }

  GBMA <- eva$GB_MA

  if(all_two_level == FALSE){
    GB_orth <- particle_orth(particle = GB, factor_level = factor_level)
  }
  GB_U <- general_create_U(all_two_level = all_two_level,
                           particle = GB, orth_all = GB_orth,
                           multiply_len = multiply_len)

  result <- list("X" = X, "LB" = LB, "GB" = GB, "GBMA" = GBMA,
                 "U" = GB_U)
  return(result)
}

## Expand the Size of Particles
##
## Add more particles to the current particle pool, and then use \code{sibma::SIB()}
## to find global best and local best particles.
##
## @param X the \code{X} return from the previous
## \code{SIB} method.
## @param LB the \code{LB} return of the previous
## \code{SIB} method.
## @param additional_X a list representing additional particles.
## @param P_w a list. Each element is a matrix describing the orthogonal
## projection matrix onto the corresponding stratum variance.
## @param q_GB an integer describing how many columns of a particle
## should be mixed with the corresponding columns of \code{GB}.
## @param q_LB an integer describing how many columns of a particle
## should be mixed with the corresponding columns of \code{LB}.
## @param q_new an integer describing how many columns of a particle
## should be mixed with the corresponding columns of a random
## particle created by \code{sibma::create_particle()}.
## @param all_two_level logical. If all factors are two levels,
## \code{all_two_level} should be \code{TRUE}; otherwise it should be
## \code{FALSE}.
## @param factor_level a list. Each element is a numeric vector specifying
## levels of a factor.
## @param t an integer specifying additional iteration number.
## @param total_unit an integer representing the number of total run size in
## a full factorial design.
## @param multiply_len a numeric vector. Each element is used to modify
## the value of each column in a model matrix.
## @param incidence_matrix a list. Each element is a matrix specifying an
## incidence matrix.
## @param treatment_effect a list. Each element is a numeric vector specifying
## which factors should have the same treatment effect in terms of a certain
## incidence matrix. The order of elements in \code{treatment_effect} should be
## corresponding to that in \code{incidence_matrix}.
SIB_relay <- function(X, LB, additional_X, P_w, q_GB, q_LB, q_new,
                      all_two_level, factor_level, t,
                      total_unit, multiply_len, incidence_matrix,
                      treatment_effect){
  unit <- nrow(X[[1]][[1]])
  factor_number <- ncol(X[[1]][[1]])
  X <- append(X, additional_X)
  LB <- append(LB, additional_X)

  eva <- evaluation(LB = LB, evaluate_GB = TRUE, factor_number = factor_number,
                    P_w = P_w, all_two_level = all_two_level,
                    factor_level = factor_level, total_unit = total_unit,
                    unit = unit, multiply_len = multiply_len)
  GB <- eva$GB
  ##########################
  # t = 1 and upper
  for(a in 1:t){
    # mix and move step
    mixwGB <- lapply(X, mix_with_GB, all_two_level = all_two_level,
                     GB = GB, P_w = P_w,
                     q_GB = q_GB, factor_level = factor_level,
                     total_unit = total_unit, unit = unit)

    mixwLB <- lapply(1:length(X), function(i){
      return(mix_with_LB(X = X[[i]], all_two_level = all_two_level,
                         LB = LB[[i]], P_w = P_w, q_LB = q_LB,
                         factor_level = factor_level, total_unit = total_unit,
                         unit = unit))
    })

    X <- lapply(1:length(X), function(i){
      return(move(X = X[[i]], mixwGB = mixwGB[[i]], mixwLB = mixwLB[[i]],
                  t_1 = FALSE, P_w = P_w,
                  q_new = q_new, all_two_level = all_two_level,
                  factor_level = factor_level, total_unit = total_unit,
                  unit = unit, multiply_len = multiply_len,
                  incidence_matrix = incidence_matrix,
                  treatment_effect = treatment_effect))
    })

    # evaluate LB and GB
    LB <- lapply(1:length(X), function(i){
      return(evaluation(X = X[[i]], LB_0 = LB[[i]], evaluate_GB = FALSE,
                        factor_number = factor_number, P_w = P_w,
                        all_two_level = all_two_level,
                        factor_level = factor_level,
                        total_unit = total_unit, unit = unit,
                        multiply_len = multiply_len))
    })
    eva <- evaluation(LB = LB, evaluate_GB = TRUE, factor_number = factor_number,
                      P_w = P_w, all_two_level = all_two_level,
                      factor_level = factor_level, total_unit = total_unit,
                      unit = unit, multiply_len = multiply_len)
    GB <- eva$GB
  }

  GBMA <- eva$GB_MA

  if(all_two_level == FALSE){
    GB_orth <- particle_orth(particle = GB, factor_level = factor_level)
  }
  GB_U <- general_create_U(all_two_level = all_two_level,
                           particle = GB, orth_all = GB_orth,
                           multiply_len = multiply_len)

  result <- list("X" = X, "LB" = LB, "GB" = GB, "GBMA" = GBMA,
                 "U" = GB_U)
  return(result)
}

## Only Increase the number of Iteration
##
## @param X the \code{X} return of the previous
## \code{SIB} method.
## @param LB the \code{LB} return of the previous
## \code{SIB} method.
## @param GB the \code{GB} return of the previous
## \code{SIB} method.
## @param P_w a list. Each element is a matrix describing the orthogonal
## projection matrix onto the corresponding stratum variance.
## @param q_GB an integer describing how many columns of a particle
## should be mixed with the corresponding columns of \code{GB}.
## @param q_LB an integer describing how many columns of a particle
## should be mixed with the corresponding columns of \code{LB}.
## @param q_new an integer describing how many columns of a particle
## should be mixed with the corresponding columns of a random
## particle created by \code{sibma::create_particle()}.
## @param all_two_level logical. If all factors are two levels,
## \code{all_two_level} should be \code{TRUE}; otherwise it should be
## \code{FALSE}.
## @param factor_level a list. Each element is a numeric vector specifying
## levels of a factor.
## @param t an integer specifying additional iteration number.
## @param total_unit an integer representing the number of total run size in
## a full factorial design.
## @param multiply_len a numeric vector. Each element is used to modify
## the value of each column in a model matrix.
## @param incidence_matrix a list. Each element is a matrix specifying an
## incidence matrix.
## @param treatment_effect a list. Each element is a numeric vector specifying
## which factors should have the same treatment effect in terms of a certain
## incidence matrix. The order of elements in \code{treatment_effect} should be
## corresponding to that in \code{incidence_matrix}.
SIB_relay_only_t <- function(X, LB, GB, P_w, q_GB, q_LB, q_new, all_two_level,
                             factor_level, t, total_unit,
                             multiply_len, incidence_matrix,
                             treatment_effect){
  unit <- nrow(X[[1]][[1]])
  factor_number <- ncol(X[[1]][[1]])

  # t = 1 and upper
  for(a in 1:t){
    # mix and move step
    mixwGB <- lapply(X, mix_with_GB, all_two_level = all_two_level,
                     GB = GB, P_w = P_w,
                     q_GB = q_GB, factor_level = factor_level,
                     total_unit = total_unit, unit = unit)
    mixwLB <- lapply(1:length(X), function(i){
      return(mix_with_LB(X = X[[i]], all_two_level = all_two_level,
                         LB = LB[[i]], P_w = P_w, q_LB = q_LB,
                         factor_level = factor_level, total_unit = total_unit,
                         unit = unit))
    })
    X <- lapply(1:length(X), function(i){
      return(move(X = X[[i]], mixwGB = mixwGB[[i]], mixwLB = mixwLB[[i]],
                  t_1 = FALSE, P_w = P_w,
                  q_new = q_new, all_two_level = all_two_level,
                  factor_level = factor_level, total_unit = total_unit,
                  unit = unit, multiply_len = multiply_len,
                  incidence_matrix = incidence_matrix,
                  treatment_effect = treatment_effect))
    })
    # evaluate LB and GB
    LB <- lapply(1:length(X), function(i){
      return(evaluation(X = X[[i]], LB_0 = LB[[i]], evaluate_GB = FALSE,
                        factor_number = factor_number, P_w = P_w,
                        all_two_level = all_two_level,
                        factor_level = factor_level,
                        total_unit = total_unit, unit = unit,
                        multiply_len = multiply_len))
    })
    eva <- evaluation(LB = LB, evaluate_GB = TRUE, factor_number = factor_number,
                      P_w = P_w, all_two_level = all_two_level,
                      factor_level = factor_level, total_unit = total_unit,
                      unit = unit, multiply_len = multiply_len)
    GB <- eva$GB
  }
  GBMA <- eva$GB_MA

  if(all_two_level == FALSE){
    GB_orth <- particle_orth(particle = GB, factor_level = factor_level)
  }
  GB_U <- general_create_U(all_two_level = all_two_level,
                           particle = GB, orth_all = GB_orth,
                           multiply_len = multiply_len)

  result <- list("X" = X, "LB" = LB, "GB" = GB, "GBMA" = GBMA,
                 "U" = GB_U)
  return(result)
}
