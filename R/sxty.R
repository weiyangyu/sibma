## Conduct Multiple SIB
## @param factor_level a list. Each element is a numeric vector specifying
## levels of a factor.
## @param unit an integer describing the number of experimental units to be
## used in a particle.
## @param particle_number an integer describing the number of particles.
## @param SIB_time an integer indicating the number of SIB process you would
## like to run.
## @param all_two_level logical. If all factors are two levels,
## \code{all_two_level} should be \code{TRUE}; otherwise it should be
## \code{FALSE}.
## @param P_w a list. Each element is a matrix describing the orthogonal
## projection matrix onto the corresponding stratum variance.
## @param q_GB an integer describing how many columns of a particle
## should be mixed with the corresponding columns of \code{GB}.
## @param q_LB an integer describing how many columns of a particle
## should be mixed with the corresponding columns of \code{LB}.
## @param q_new an integer describing how many columns of a particle
## should be mixed with the corresponding columns of a random
## particle created by \code{sibma::create_particle()}.
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
sxty_SIB <- function(factor_level, unit, particle_number, SIB_time,
                     all_two_level, P_w, q_GB, q_LB, q_new, t, total_unit,
                     multiply_len, incidence_matrix, treatment_effect){

  test_particle <- lapply(1:SIB_time, function(i){
    return(create_particle(factor_level = factor_level,
                           unit = unit,
                           particle_number = particle_number,
                           incidence_matrix = incidence_matrix,
                           treatment_effect = treatment_effect))
  })

  testSIB <- lapply(1:SIB_time, function(i){
    return(SIB(all_two_level = all_two_level, factor_level = factor_level,
               X = test_particle[[i]], P_w = P_w,
               q_GB = q_GB, q_LB = q_LB, q_new = q_new, t = t,
               total_unit = total_unit, multiply_len = multiply_len,
               incidence_matrix = incidence_matrix,
               treatment_effect = treatment_effect))
  })

  MA_df <- lapply(1:SIB_time, function(i){
    lapply(testSIB[[i]]$GBMA, function(j){
      matrix(unlist(j), ncol = length(factor_level), byrow = TRUE)
    })
  })

  history <- data.frame(matrix(c(particle_number, t), ncol = 2))
  colnames(history) <- c("n_particle", "iteration")

  result <- list("testSIB" = testSIB, "history" = history, "MA_df" = MA_df)
}

## Conduct multiple SIB with Increased Size of Particles
## @param prior_SIB the \code{testSIB} return of the previous
## \code{sxty_SIB} method.
## @param history the \code{history} return of the previous
## \code{sxty_SIB} method.
## @param factor_level a list. Each element is a numeric vector specifying
## levels of a factor.
## @param unit an integer describing the number of experimental units to be
## used in a particle.
## @param particle_increase an integer describing the number of increased particles.
## @param SIB_time an integer indicating the number of SIB process you would
## like to run.
## @param all_two_level logical. If all factors are two levels,
## \code{all_two_level} should be \code{TRUE}; otherwise it should be
## \code{FALSE}.
## @param P_w a list. Each element is a matrix describing the orthogonal
## projection matrix onto the corresponding stratum variance.
## @param q_GB an integer describing how many columns of a particle
## should be mixed with the corresponding columns of \code{GB}.
## @param q_LB an integer describing how many columns of a particle
## should be mixed with the corresponding columns of \code{LB}.
## @param q_new an integer describing how many columns of a particle
## should be mixed with the corresponding columns of a random
## particle created by \code{sibma::create_particle()}.
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
sxty_SIB_relay <- function(prior_SIB, history, factor_level, unit,
                           particle_increase, SIB_time,
                           all_two_level, P_w, q_GB, q_LB, q_new, t,
                           total_unit, multiply_len,
                           incidence_matrix,
                           treatment_effect){

  additional_X <- lapply(1:SIB_time, function(i){
    lapply(1:particle_increase, function(j){
      create_particle(factor_level = factor_level, unit = unit,
                      particle_number = length(j),
                      incidence_matrix = incidence_matrix,
                      treatment_effect = treatment_effect)
    })
  })

  testSIB <- lapply(1:SIB_time, function(i){
    return(SIB_relay(X = prior_SIB[[i]]$X, LB = prior_SIB[[i]]$LB,
                     additional_X = additional_X[[i]],
                     P_w = P_w, q_GB = q_GB, q_LB = q_LB, q_new = q_new,
                     all_two_level = all_two_level,
                     factor_level = factor_level, t = t,
                     total_unit = total_unit, multiply_len = multiply_len,
                     incidence_matrix = incidence_matrix,
                     treatment_effect = treatment_effect))
  })

  MA_df <- lapply(1:SIB_time, function(i){
    lapply(testSIB[[i]]$GBMA, function(j){
      matrix(unlist(j), ncol = length(factor_level), byrow = TRUE)
    })
  })

  temp_history <- c(particle_increase + history[nrow(history),1],
                    t + history[nrow(history),2])
  history <- rbind(history, temp_history)

  result <- list("testSIB" = testSIB, "history" = history, "MA_df" = MA_df)
}

## Conduct multiple SIB with Increased Number of Iteration
## @param prior_SIB the \code{testSIB} return of the previous
## \code{sxty_SIB} method.
## @param history the \code{history} return of the previous
## \code{sxty_SIB} method.
## @param factor_level a list. Each element is a numeric vector specifying
## levels of a factor.
## @param SIB_time an integer indicating the number of SIB process you would
## like to run.
## @param all_two_level logical. If all factors are two levels,
## \code{all_two_level} should be \code{TRUE}; otherwise it should be
## \code{FALSE}.
## @param P_w a list. Each element is a matrix describing the orthogonal
## projection matrix onto the corresponding stratum variance.
## @param q_GB an integer describing how many columns of a particle
## should be mixed with the corresponding columns of \code{GB}.
## @param q_LB an integer describing how many columns of a particle
## should be mixed with the corresponding columns of \code{LB}.
## @param q_new an integer describing how many columns of a particle
## should be mixed with the corresponding columns of a random
## particle created by \code{sibma::create_particle()}.
## @param t an integer indicating the increased number of iteration.
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
sxty_SIB_relay_only_t <- function(prior_SIB, history, factor_level, SIB_time,
                                  all_two_level, P_w, q_GB, q_LB, q_new,
                                  t, total_unit, multiply_len,
                                  incidence_matrix, treatment_effect){

  testSIB <- lapply(1:SIB_time, function(i){
    return(SIB_relay_only_t(X = prior_SIB[[i]]$X, LB = prior_SIB[[i]]$LB,
                            GB = prior_SIB[[i]]$GB, P_w = P_w,
                            q_GB = q_GB, q_LB = q_LB, q_new = q_new,
                            all_two_level = all_two_level,
                            factor_level = factor_level, t = t,
                            total_unit = total_unit,
                            multiply_len = multiply_len,
                            incidence_matrix = incidence_matrix,
                            treatment_effect = treatment_effect))
  })

  MA_df <- lapply(1:SIB_time, function(i){
    lapply(testSIB[[i]]$GBMA, function(j){
      matrix(unlist(j), ncol = length(factor_level), byrow = TRUE)
    })
  })

  temp_history <- c(history[nrow(history),1],
                    t + history[nrow(history),2])
  history <- rbind(history, temp_history)

  result <- list("testSIB" = testSIB, "history" = history, "MA_df" = MA_df)
}
