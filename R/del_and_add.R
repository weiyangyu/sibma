## Delete some columns within a particle
## @param particle a matrix representing a particle.
## @param all_two_level logical. If all factors are two levels,
## \code{all_two_level} should be \code{TRUE}; otherwise it should be
## \code{FALSE}.
## @param factor_level a list. Each element is a numeric vector specifying
## levels of a factor.
## @param P_w a list. Each element is a matrix describing the orthogonal
## projection matrix onto the corresponding stratum variance.
## @param q an integer describing how many columns should be deleted.
## @param total_unit an integer representing the number of total run size in a
## full factorial design.
## @param unit an integer describing the number of experimental units to be used
## in a particle.
deletion <- function(particle, all_two_level, factor_level, P_w, q,
                     total_unit, unit){
  left <- 1:ncol(particle)
  for(i in 1:q){
    deleted_particle <- lapply(1:ncol(particle), function(j){
      return(particle[,-j])
    })
    deleted_factor_level <- lapply(1:ncol(particle), function(j){
      return(factor_level[-j])
    })
    if(all_two_level == FALSE){
      deleted_orth <- lapply(1:ncol(particle), function(j){
        orthParticle <- particle_orth(particle = list(deleted_particle[[j]]),
                                      factor_level = deleted_factor_level[[j]])
        return(orthParticle)
      })
    }
    multiply_len <- lapply(deleted_factor_level, function(j){
      return(rep(1, prod(lengths(j))))
    })
    temp_U <- lapply(1:length(multiply_len), function(j){
      return(general_create_U(all_two_level = all_two_level,
                              particle = list(deleted_particle[[j]]),
                              orth_all = deleted_orth[[j]],
                              multiply_len = multiply_len[[j]]))
    })
    factor_number <- ncol(particle) - 1
    MA <- lapply(temp_U, function(j){
      wordlength_list <- lapply(1:length(P_w), function(k){
        return(worlen_pattern(factor_number = factor_number,
                              U = j[[1]],
                              P_w = P_w[[k]],
                              total_unit = total_unit,
                              unit = unit))
      })
      wordlength <- Reduce("+", wordlength_list)
      return(wordlength)
    })
    comparison <- compare_MA_min(particle = deleted_particle, MA = MA)
    particle <- comparison$GB[[1]]
    factor_level <- factor_level[-comparison$bestset[1]]
    left <- left[-comparison$bestset[1]]
  }
  ans <- list("temp_GB" = particle, "left" = left)
  return(ans)
}

## Add some columns back to a particle
## @param particle a matrix representing a particle.
## @param GB a matrix representing a global best particle.
## @param left the \code{left} return from \code{sibma::deletion}.
addition <- function(particle, GB, left){
  for(i in 1:length(left)){
    GB[,left[i]] <- particle[,i]
  }
  return(GB)
}

