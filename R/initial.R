## Initial step
## @param particle a list containing all particles.
## @param all_two_level logical. If all factors are two levels,
## \code{all_two_level} should be \code{TRUE}; otherwise it should be
## \code{FALSE}.
## @param P_w a list. Each element is a matrix describing the orthogonal
## projection matrix onto the corresponding stratum variance.
## @param factor_level a list. Each element is a numeric vector specifying
## levels of a factor.
## @param total_unit an integer representing the number of total run size in a
## full factorial design.
## @param unit an integer describing the number of experimental units to be used in a
## particle.
## @param multiply_len a numeric vector. Each element is used to modify
## the value of each column in a model matrix.
initial <- function(particle, all_two_level, P_w, factor_level,
                    total_unit, unit, multiply_len){
  factor_number <- ncol(particle[[1]])

  if(all_two_level == TRUE){
    U <- general_create_U(all_two_level = all_two_level,
                          particle = particle,
                          multiply_len = multiply_len)
  }else{
    orthParticle <- particle_orth(particle = particle, factor_level = factor_level)
    U <- general_create_U(all_two_level = all_two_level,
                          orth_all =  orthParticle,
                          multiply_len = multiply_len)
  }

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
    comparison <- compare_MA_min(particle = particle, MA = MA[[i]])
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
    if(Reduce(identical, GB) == TRUE){
      GB <- list(GB[[1]])
    }
  }

  particle <- lapply(particle, function(i){
    return(list(i))
  })

  result <- list("X" = particle, "LB" = particle, "GB" = GB)
  return(result)
}
