## Update Current Particles
## @param X a list. Each element is a particle ready to be mixed.
## @param mixwGB a list. Each element is a particle after \code{X} mixed with
## \code{GB}.
## @param mixwLB a list. Each element is a particle after \code{X} mixed with
## \code{LB}.
## @param t_1 logical. If this is the first move step, \code{t_1} should be
## \code{TRUE}; otherwise, it should be \code{FALSE}.
## @param P_w a list. Each element is a matrix describing the orthogonal
## projection matrix onto the corresponding stratum variance.
## @param q_new an integer describing how many columns of each element of\code{X}
## should be mixed with the corresponding columns of a new particle created from
## \code{sibma::create_particle()}.
## @param all_two_level logical. If all factors are two levels,
## \code{all_two_level} should be \code{TRUE}; otherwise it should be
## \code{FALSE}.
## @param factor_level a list. Each element is a numeric vector specifying
## levels of a factor.
## @param total_unit an integer representing the number of total run size in
## a full factorial design.
## @param unit an integer describing the number of experimental units to be used in a
## particle.
## @param multiply_len a numeric vector. Each element is used to modify
## the value of each column in a model matrix.
## @param incidence_matrix a list. Each element is a matrix specifying an
## incidence matrix.
## @param treatment_effect a list. Each element is a numeric vector specifying
## which factors should have the same treatment effect in terms of a certain
## incidence matrix. The order of elements in \code{treatment_effect} should be
## corresponding to that in \code{incidence_matrix}.
move <- function(X, mixwGB, mixwLB, t_1, P_w, q_new, all_two_level,
                 factor_level, total_unit, unit, multiply_len,
                 incidence_matrix, treatment_effect){

  # transform the type of X to list
  if(is.data.frame(X) == TRUE){
    factor_number <- ncol(X)
    X <- list(X)
  }else{
    factor_number <- ncol(X[[1]])
  }

  # transform the type of mixwGB to list
  if(is.data.frame(mixwGB) == TRUE){
    mixwGB <- list(mixwGB)
  }

  if(t_1 == TRUE){
    candidate <- append(X, mixwGB)
    if(all_two_level == FALSE){
      orthParticle <- particle_orth(particle = candidate, factor_level = factor_level)
    }
    U <- general_create_U(all_two_level = all_two_level,
                          particle = candidate, orth_all = orthParticle,
                          multiply_len = multiply_len)

    MA <- lapply(P_w, function(i){
      for_each_P_w <- lapply(U, function(ii){
        wordlength_list <- lapply(1:length(i), function(j){
          return(worlen_pattern(factor_number = factor_number,
                                U = ii,
                                P_w = i[[j]],
                                total_unit = total_unit,
                                unit = unit))
        })
        wordlength <- Reduce("+", wordlength_list)
        return(wordlength)
      })
      return(for_each_P_w)
    })

    GB <- lapply(1:length(MA), function(i){
      comparison <- compare_MA_min(particle = candidate, MA = MA[[i]])
      remain_MA <- MA[-i]

      if(length(MA) > 1){
        if(length(comparison$bestset) > 1){
          for(j in 1:length(remain_MA)){
            comparison <- compare_MA_min(particle = comparison$GB,
                                         MA = remain_MA[[j]][comparison$bestset])
            if(length(comparison$bestset) == 1){
              break
            }
          }
        }
      }
      return(comparison$GB[[1]])
    })

    if(length(MA) > 1){
      if(all_equal_list(GB) == TRUE){
        GB <- list(GB[[1]])
      }
    }

    GB_equal_X <- lapply(GB, function(i){
      return(lapply(X, function(j){
        return(identical(i, j))
      }))
    })

    if(any(unlist(GB_equal_X)) == TRUE){
      X <- mix_with_new(X = X, all_two_level = all_two_level,
                        P_w = P_w, q_new = q_new,
                        factor_level = factor_level, unit = unit,
                        total_unit = total_unit,
                        incidence_matrix = incidence_matrix,
                        treatment_effect = treatment_effect)
    }else{
      X <- GB
    }
  }else{

    # transform the type of mixwLB to list
    if(is.data.frame(mixwLB) == TRUE){
      mixwLB <- list(mixwLB)
    }
    candidate <- Reduce(append, list(X, mixwGB, mixwLB))
    if(all_two_level == FALSE){
      orthParticle <- particle_orth(particle = candidate, factor_level = factor_level)
    }
    U <- general_create_U(all_two_level = all_two_level,
                          particle = candidate, orth_all = orthParticle,
                          multiply_len = multiply_len)
    MA <- lapply(P_w, function(i){
      for_each_P_w <- lapply(U, function(ii){
        wordlength_list <- lapply(1:length(i), function(j){
          return(worlen_pattern(factor_number = factor_number,
                                U = ii,
                                P_w = i[[j]],
                                total_unit = total_unit,
                                unit = unit))
        })
        wordlength <- Reduce("+", wordlength_list)
        return(wordlength)
      })
      return(for_each_P_w)
    })

    GB <- lapply(1:length(MA), function(i){
      comparison <- compare_MA_min(particle = candidate, MA = MA[[i]])
      remain_MA <- MA[-i]

      if(length(MA) > 1){
        if(length(comparison$bestset) > 1){
          for(j in 1:length(remain_MA)){
            comparison <- compare_MA_min(particle = comparison$GB,
                                         MA = remain_MA[[j]][comparison$bestset])
            if(length(comparison$bestset) == 1){
              break
            }
          }
        }
      }
      return(comparison$GB[[1]])
    })

    if(length(MA) > 1){
      if(all_equal_list(GB) == TRUE){
        GB <- list(GB[[1]])
      }
    }

    GB_equal_X <- lapply(GB, function(i){
      return(lapply(X, function(j){
        return(identical(i, j))
      }))
    })

    if(any(unlist(GB_equal_X)) == TRUE){
      X <- mix_with_new(X = X, all_two_level = all_two_level,
                        P_w = P_w, q_new = q_new,
                        factor_level = factor_level, unit = unit,
                        total_unit = total_unit,
                        incidence_matrix = incidence_matrix,
                        treatment_effect = treatment_effect)
    }else{
      X <- GB
    }

  }
  return(X)
}

## Update Local best or Global Best Particles
##
## @param X a list. Each element is a particle after \code{sibma::move()}.
## @param LB_0 representing the current local best particles of \code{X}.
## @param LB a list representing all \code{LB} particles.
## @param evaluate_GB logical. If you want to update the \code{GB} particle,
## \code{evaluate_GB} should be \code{TRUE}; otherwise, it should be
## \code{FALSE}.
## @param factor_number an integer describing the number of factors.
## @param P_w a list. Each element is a matrix describing the orthogonal
## projection matrix onto the corresponding stratum variance.
## @param all_two_level logical. If all factors are two levels,
## \code{all_two_level} should be \code{TRUE}; otherwise it should be
## \code{FALSE}.
## @param factor_level a list. Each element is a numeric vector specifying
## levels of a factor.
## @param total_unit an integer representing the number of total run size in
## a full factorial design.
## @param unit an integer describing the number of experimental units to be used in a
## particle.
## @param multiply_len a numeric vector. Each element is used to modify
## the value of each column in a model matrix.
evaluation <- function(X, LB_0, LB, evaluate_GB, factor_number, P_w, all_two_level,
                       factor_level, total_unit, unit, multiply_len){
  if(evaluate_GB == TRUE){

    # check if every element in LB is a data frame
    is_df <- lapply(LB, function(i){
      return(is.data.frame(i))
    })
    # transform those elements in LB from data frame to list
    LB <- lapply(1:length(is_df), function(i){
      if(is_df[[i]] == TRUE){
        return(list(LB[[i]]))
      }else{
        return(LB[[i]])
      }
    })

    LB <- Reduce(append, LB)

    if(all_two_level == FALSE){
      orthParticle <- particle_orth(particle = LB, factor_level = factor_level)
    }
    U <- general_create_U(all_two_level = all_two_level,
                          particle = LB, orth_all = orthParticle,
                          multiply_len = multiply_len)

    MA <- lapply(P_w, function(i){
      for_each_P_w <- lapply(U, function(ii){
        wordlength_list <- lapply(1:length(i), function(j){
          return(worlen_pattern(factor_number = factor_number,
                                U = ii,
                                P_w = i[[j]],
                                total_unit = total_unit,
                                unit = unit))
        })
        wordlength <- Reduce("+", wordlength_list)
        return(wordlength)
      })
      return(for_each_P_w)
    })

    GB <- lapply(1:length(MA), function(i){
      comparison <- compare_MA_min(particle = LB, MA = MA[[i]])
      remain_MA <- MA[-i]

      if(length(MA) > 1){
        if(length(comparison$bestset) > 1){
          for(j in 1:length(remain_MA)){
            comparison <- compare_MA_min(particle = comparison$GB,
                                         MA = remain_MA[[j]][comparison$bestset])
            if(length(comparison$bestset) == 1){
              break
            }
          }
        }
      }
      return(comparison$GB[[1]])
    })

    if(length(MA) > 1){
      if(all_equal_list(GB) == TRUE){
        GB <- list(GB[[1]])
      }
    }

    if(all_two_level == FALSE){
      orthParticle <- particle_orth(particle = GB, factor_level = factor_level)
    }
    U <- general_create_U(all_two_level = all_two_level,
                          particle = GB, orth_all = orthParticle,
                          multiply_len = multiply_len)

    GB_MA <- lapply(P_w, function(i){
      for_each_P_w <- lapply(U, function(ii){
        wordlength_list <- lapply(1:length(i), function(j){
          return(worlen_pattern(factor_number = factor_number,
                                U = ii,
                                P_w = i[[j]],
                                total_unit = total_unit,
                                unit = unit))
        })
        wordlength <- Reduce("+", wordlength_list)
        return(wordlength)
      })
      return(for_each_P_w)
    })

    ans <- list("GB" = GB, "GB_MA" = GB_MA)

  }else{

    if(is.data.frame(X) == TRUE){
      X <- list(X)
    }
    if(is.data.frame(LB_0) == TRUE){
      LB_0 <- list(LB_0)
    }

    LB_candidate <- append(X, LB_0)

    if(all_two_level == FALSE){
      orthParticle <- particle_orth(particle = LB_candidate, factor_level = factor_level)
    }

    U <- general_create_U(all_two_level = all_two_level,
                          particle = LB_candidate, orth_all = orthParticle,
                          multiply_len = multiply_len)

    MA <- lapply(P_w, function(i){
      for_each_P_w <- lapply(U, function(ii){
        wordlength_list <- lapply(1:length(i), function(j){
          return(worlen_pattern(factor_number = factor_number,
                                U = ii,
                                P_w = i[[j]],
                                total_unit = total_unit,
                                unit = unit))
        })
        wordlength <- Reduce("+", wordlength_list)
        return(wordlength)
      })
      return(for_each_P_w)
    })

    GB <- lapply(1:length(MA), function(i){
      comparison <- compare_MA_min(particle = LB_candidate, MA = MA[[i]])
      remain_MA <- MA[-i]

      if(length(MA) > 1){
        if(length(comparison$bestset) > 1){
          for(j in 1:length(remain_MA)){
            comparison <- compare_MA_min(particle = comparison$GB,
                                         MA = remain_MA[[j]][comparison$bestset])
            if(length(comparison$bestset) == 1){
              break
            }
          }
        }
      }
      return(comparison$GB[[1]])
    })

    if(length(MA) > 1){
      if(all_equal_list(GB) == TRUE){
        GB <- list(GB[[1]])
      }
    }
    ans <- GB
  }

  return(ans)
}



