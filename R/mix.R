## Mix with Global Best Particles
## @param X a list. Each element is a particle ready to be mixed.
## @param all_two_level logical. If all factors are two levels,
## \code{all_two_level} should be \code{TRUE}; otherwise it should be
## \code{FALSE}.
## @param GB a list. Each element is a global best particle.
## @param P_w a list. Each element is a matrix describing the orthogonal
## projection matrix onto the corresponding stratum variance.
## @param q_GB an integer describing how many columns of each element of\code{X}
## should be mixed with the corresponding columns of each element of \code{GB}.
## @param factor_level a list. Each element is a numeric vector specifying
## levels of a factor.
## @param total_unit an integer representing the number of total run size in a
## full factorial design.
## @param unit an integer describing the number of experimental units to be used
## in a particle.
mix_with_GB <- function(X, all_two_level, GB, P_w, q_GB,
                        factor_level, total_unit, unit){
  if(length(X) == 1){
    del_X_GB <- lapply(P_w, function(i){
      return(deletion(particle = X[[1]], all_two_level = all_two_level,
                      factor_level = factor_level,
                      P_w = i, q = q_GB, total_unit = total_unit,
                      unit = unit))
    })

    if(length(GB) == 1){
      mixw_GB <- lapply(1:length(P_w), function(i){
        return(addition(particle = del_X_GB[[i]]$temp_GB, GB = GB[[1]], left = del_X_GB[[i]]$left))
      })
    }else{
      mixw_GB <- lapply(1:length(P_w), function(i){
        return(addition(particle = del_X_GB[[i]]$temp_GB, GB = GB[[i]], left = del_X_GB[[i]]$left))
      })
    }
  }else{
    del_X_GB <- lapply(1:length(P_w), function(i){
      return(deletion(particle = X[[i]], all_two_level = all_two_level,
                      factor_level = factor_level,
                      P_w = P_w[[i]], q = q_GB, total_unit = total_unit,
                      unit = unit))
    })

    if(length(GB) == 1){
      mixw_GB <- lapply(1:length(P_w), function(i){
        return(addition(particle = del_X_GB[[i]]$temp_GB, GB = GB[[1]], left = del_X_GB[[i]]$left))
      })
    }else{
      mixw_GB <- lapply(1:length(P_w), function(i){
        return(addition(particle = del_X_GB[[i]]$temp_GB, GB = GB[[i]], left = del_X_GB[[i]]$left))
      })
    }
  }
  return(mixw_GB)
}

## Mix with Local Best Particles
## @param X a list. Each element is a particle ready to be mixed.
## @param all_two_level logical. If all factors are two levels,
## \code{all_two_level} should be \code{TRUE}; otherwise it should be
## \code{FALSE}.
## @param LB a list. Each element is a local best particle.
## @param P_w a list. Each element is a matrix describing the orthogonal
## projection matrix onto the corresponding stratum variance.
## @param q_LB an integer describing how many columns of each element of\code{X}
## should be mixed with the corresponding columns of each element of \code{LB}.
## @param factor_level a list. Each element is a numeric vector specifying
## levels of a factor.
## @param total_unit an integer representing the number of total run size in a
## full factorial design.
## @param unit an integer describing the number of experimental units to be used
## in a particle.
mix_with_LB <- function(X, all_two_level, LB, P_w, q_LB, factor_level,
                        total_unit, unit){
  if(length(X) == 1){
    del_X_LB <- lapply(P_w, function(i){
      return(deletion(particle = X[[1]], all_two_level = all_two_level,
                      factor_level = factor_level,
                      P_w = i, q = q_LB, total_unit = total_unit,
                      unit = unit))
    })

    if(length(LB) == 1){
      mixw_LB <- lapply(1:length(P_w), function(i){
        return(addition(particle = del_X_LB[[i]]$temp_GB, GB = LB[[1]], left = del_X_LB[[i]]$left))
      })
    }else{
      mixw_LB <- lapply(1:length(P_w), function(i){
        return(addition(particle = del_X_LB[[i]]$temp_GB, GB = LB[[i]], left = del_X_LB[[i]]$left))
      })
    }
  }else{
    del_X_LB <- lapply(1:length(P_w), function(i){
      return(deletion(particle = X[[i]], all_two_level = all_two_level,
                      factor_level = factor_level,
                      P_w = P_w[[i]], q = q_LB, total_unit = total_unit,
                      unit = unit))
    })

    if(length(LB) == 1){
      mixw_LB <- lapply(1:length(P_w), function(i){
        return(addition(particle = del_X_LB[[i]]$temp_GB, GB = LB[[1]], left = del_X_LB[[i]]$left))
      })
    }else{
      mixw_LB <- lapply(1:length(P_w), function(i){
        return(addition(particle = del_X_LB[[i]]$temp_GB, GB = LB[[i]], left = del_X_LB[[i]]$left))
      })
    }
  }
  return(mixw_LB)
}

## Mix with a New Particle
##
## @param X a list. Each element is a particle ready to be mixed.
## @param all_two_level logical. If all factors are two levels,
## \code{all_two_level} should be \code{TRUE}; otherwise it should be
## \code{FALSE}.
## @param P_w a list. Each element is a matrix describing the orthogonal
## projection matrix onto the corresponding stratum variance.
## @param q_new an integer describing how many columns of each element of\code{X}
## should be mixed with the corresponding columns of a new particle created from
## \code{sibma::create_particle()}.
## @param factor_level a list. Each element is a numeric vector specifying
## levels of a factor.
## @param total_unit an integer representing the number of total run size in a
## full factorial design.
## @param unit an integer describing the number of experimental units to be used
## in a particle.
## @param incidence_matrix a list. Each element is a matrix specifying an
## incidence matrix.
## @param treatment_effect a list. Each element is a numeric vector specifying
## which factors should have the same treatment effect in terms of a certain
## incidence matrix. The order of elements in \code{treatment_effect} should be
## corresponding to that in \code{incidence_matrix}.
mix_with_new <- function(X, all_two_level, P_w, q_new, factor_level, total_unit,
                         unit, incidence_matrix, treatment_effect){
  new <- create_particle(factor_level = factor_level, unit = unit,
                         particle_number = 1, incidence_matrix = incidence_matrix,
                         treatment_effect = treatment_effect)

  if(length(X) == 1){
    del_X_new <- lapply(P_w, function(i){
      return(deletion(particle = X[[1]], all_two_level = all_two_level,
                      factor_level = factor_level,
                      P_w = i, q = q_new, total_unit = total_unit,
                      unit = unit))
    })
    mixw_new <- lapply(1:length(P_w), function(i){
      return(addition(particle = del_X_new[[i]]$temp_GB, GB = new[[1]],
                      left = del_X_new[[i]]$left))
    })
  }else{
    del_X_new <- lapply(1:length(P_w), function(i){
      return(deletion(particle = X[[i]], all_two_level = all_two_level,
                      factor_level = factor_level,
                      P_w = P_w[[i]], q = q_new, total_unit = total_unit,
                      unit = unit))
    })
    mixw_new <- lapply(1:length(P_w), function(i){
      return(addition(particle = del_X_new[[i]]$temp_GB, GB = new[[1]],
                      left = del_X_new[[i]]$left))
    })
  }

  return(mixw_new)
}
