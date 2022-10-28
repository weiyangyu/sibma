## Create a list of particles
##
## Each particle is a kind of experimental design where the number of rows
## represent the number of experimental units and the number of columns represent
## the number of factors.
##
## We highly recommend that you use \code{c(-1,1)} as a code for a 2-level
## factor. If you have a factor with p-level, we recommend you use \code{0:(p-1)}
## as the code for this factor.
## @param factor_level a list. Each element is a numeric vector specifying
## levels of a factor.
## @param unit an integer describing the number of experimental units to be used
## in a particle.
## @param particle_number an integer describing the number of particles.
## @param incidence_matrix a list. Each element is a matrix specifying an
## incidence matrix.
## @param treatment_effect a list. Each element is a numeric vector specifying
## which factors should have the same treatment effect in terms of a certain
## incidence matrix. The order of elements in \code{treatment_effect} should be
## corresponding to that in \code{incidence_matrix}.
create_particle <- function(factor_level, unit, particle_number,
                            incidence_matrix, treatment_effect){
  # Return a list. Each element is an integer representing each level's occurrence
  # number of a factor.
  replicate <- lapply(1:length(factor_level), function(i){
    return(unit / length(factor_level[[i]]))
  })
  # Return a data.frame representing a design before randomized.
  standard_design <- sapply(1:length(factor_level), function(i){
    return(rep(factor_level[[i]], replicate[[i]]))
  })
  standard_design <- data.frame(standard_design)
  # Return a list of data.frame. Each represents a design, which is also a particle.
  # Each particle is the result of randomization of standard_design by factor.
  particle <- lapply(1:particle_number, function(i){
    standard_design[] <- lapply(standard_design, sample)
    return(standard_design)
  })
  # When the design is structured.
  if(is.null(incidence_matrix) == FALSE){
    # Return a list of matrices in which columns represents blocks and rows
    # represents units in that block.
    # Imply the size of each block should be the same under a certain stratum.
    idx_same_treatment <- lapply(incidence_matrix, function(i){
      return(apply(i, 2, function(j) which(j == 1)))
    })
    # Return a list.
    particle <- lapply(particle, function(i){
      unique_treatment <- lapply(1:length(idx_same_treatment), function(j){
        return(sapply(treatment_effect[[j]], function(k){
          return(sample(factor_level[[k]], ncol(idx_same_treatment[[j]]), replace = TRUE))
        }))
      })
      for(j in 1:length(idx_same_treatment)){
        for(k in 1:ncol(idx_same_treatment[[j]])){
          for(l in 1:length(treatment_effect[[j]])){
            i[idx_same_treatment[[j]][,k], treatment_effect[[j]][l]] <- unique_treatment[[j]][k,l]
          }
        }
      }
      return(i)
    })
  }
  return(particle)
}

## Transformation: a factor's design to corresponding effects
##
## Transform a particle's column, which represents a factor's design to its
## corresponding effects, such as a linear effect for a two-level factor or
## a linear effect and a quadratic effect for a three-level factor.
## @param col_particle an integer vector specifying a particle's column.
## @param factor_level an integer vector specifying col_particle's levels.
#' @importFrom stats poly
c_particle_orth <- function(col_particle, factor_level){

  # create orthogonal matrix
  main_effect_number <- length(factor_level) - 1
  ortho_matrix <- poly(factor_level, degree = main_effect_number, simple = T)

  # transformation
  output <- matrix(NA, nrow = length(col_particle), ncol = main_effect_number)

  for(i in 1:length(factor_level)){
    ind <- which(col_particle == factor_level[i])
    for(k in ind){
      output[k,] <- ortho_matrix[i,]
    }
  }
  return(output)
}

## Transformation: All particles to corresponding effects
##
## Transform every particle, namely, every experimental design to its corresponding
## effects.
## @param particle a list containing all particles.
## @param factor_level a list. Each element is a numeric vector specifying
## levels of a factor.
particle_orth <- function(particle, factor_level){
  orth_all <- lapply(particle, function(i){
    orth_each <- lapply(1:ncol(i), function(j){
      return(c_particle_orth(col_particle = i[,j], factor_level = factor_level[[j]]))
    })
    return(orth_each)
  })
  return(orth_all)
}
