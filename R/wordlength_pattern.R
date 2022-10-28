## Create model matrices for all particles
##
## Model matrix represents an experimental design's factorial effects. For
## a full factorial experiment, each effect's euclidean distance should be the
## same, i.e., each column of the full model matrix
## should have the same euclidean distance. However, we now are discussing in
## the realm of fractional factorial experiment. And it is not realistic to
## create a full model matrix first and pick a
## fraction of it when the factor number is large. Therefore, we first transform
## each factor's design to its corresponding orthogonal matrix using
## \code{sibma::particle_orth()}, and use \code{stats::model.matrix()} to create
## the model matrix. If all factors are 2-level, and you use \code{c(-1,1)} as
## your level coding, there is no need to use \code{sibma::particle_orth()}, you
## can just create the model matrix.
##
## But it is not over yet. After using \code{sibma::particle_orth()},
## every factor's orthogonal matrix is orthonormal, but
## \code{stats::model.matrix()} will not consider the effect of intercept.
## Therefore, we use \code{multiply_len} to specify the effect of intercept.
## For instance, if you have an experimental design with two 3-level factors.
## Each factor will have an orthogonal matrix specifying the linear and
## quadratic effects. Through \code{stats::model.matrix()}, you will receive a
## model matrix with a intercept, two linear effects, two quadratic effects, and
## four interaction between linear and quadratic effects. To specify the effect
## of intercept, \code{multiply_len = c(1,rep(1/sqrt(3),4), rep(1,4))}.
## It is fine to not correct the intercept, since it will not affect particles'
## ranking in terms of wordlength pattern.
## @param all_two_level logical. If all factors are two levels,
## \code{all_two_level} should be \code{TRUE}; otherwise it should be
## \code{FALSE}.
## @param particle a list containing all particles.
## @param orth_all a list containing all orthogonal matrices generating
## from all the particles.
## @param multiply_len a numeric vector. Each element is used to modify
## the value of each column in a model matrix.
#' @importFrom stats model.matrix
general_create_U <- function(all_two_level, particle, orth_all,
                             multiply_len){
  if(all_two_level == T){
    X <- as.matrix(particle[[1]])  # any particle
    unit <- nrow(X)
    ncol_particle <- ncol(X)
    each_particle <- lapply(1:length(particle), function(i){
      each_factor <- sapply(1:ncol_particle, function(j){
        return(paste0("particle[[", i, "]][,", j, "]"))
      })
      return(each_factor)
    })
    U <- lapply(each_particle, function(i){
      ans <- do.call(paste, c(as.list(i), sep = "+"))
      ans <- paste0("lm(rep(0,", unit, ") ~ ", "(", ans, ")^", unit, ")")
      result <- model.matrix(eval(parse(text=ans)))
      result <- sweep(result, 2, multiply_len, "*")
      return(result)
    })
  }else{
    unit <- nrow(orth_all[[1]][[1]])
    len_orth_each <- length(orth_all[[1]])
    each_orth <- lapply(1:length(orth_all), function(i){
      each_factor <- sapply(1:len_orth_each, function(j){
        return(paste0("orth_all[[", i, "]][[", j, "]]"))
      })
      return(each_factor)
    })
    U <- lapply(each_orth, function(i){
      ans <- do.call(paste, c(as.list(i), sep = "+"))
      ans <- paste0("lm(rep(0,", unit, ") ~ ", "(", ans, ")^", unit, ")")
      result <- model.matrix(eval(parse(text=ans)))
      result <- sweep(result, 2, multiply_len, "*")
      return(result)
    })
  }
  return(U)
}

## Create wordlength pattern
##
## Create an experimental design's wordlength pattern.
## @param factor_number an integer describing the number of factors.
## @param U a matrix describing a model matrix calculated from
## \code{sibma::general_create_U()}.
## @param P_w a matrix describing the orthogonal projection matrix onto the
## corresponding stratum variance.
## @param total_unit an integer representing the number of total run size in a
## full factorial design.
## @param unit an integer describing the number of experimental units to be used in a
## particle.
#' @importFrom psych tr
worlen_pattern <- function(factor_number, U, P_w,
                           total_unit, unit){
  result <- rep(0, factor_number)
  colnames_U <- colnames(U)
  ind_set <- stringr::str_count(colnames_U, pattern = ":")[-1]

  for(i in 1:factor_number){
    ind <- which(ind_set == (i-1)) + 1
    result[i] <- psych::tr(t(U[,ind]) %*% P_w %*% U[,ind])
  }
  result <- result * total_unit / unit
  result <- round(result, digits = 5)

  return(result)
}

## Compare multiple word-length pattern
##
## Compare a list of word-length pattern and yield the minimum aberration design.
##
## The function will return three values: \code{best}, \code{result}, \code{GB}.
## \code{best} returns a vector of integers indicating where  best particles are.
## \code{result} returns a matrix. Each row is a word-length pattern.
## \code{GB} returns a list of matrices specifying the best particles.
##
## @param particle a list containing all particles.
## @param MA a list containing all the word-length patterns of their
## corresponding particles.
compare_MA_min <- function(particle, MA){
  # create a MA dataframe
  MA_df <- t(as.data.frame(MA))
  rownames(MA_df) <- 1:nrow(MA_df)

  # compare by column
  bestset <- 1:nrow(MA_df)

  for(i in 1:ncol(MA_df)){
    mini <- min(MA_df[bestset, i])
    bestset <- as.numeric(names(which(MA_df[bestset, i] == mini)))
  }
  result <- MA_df[bestset, ]

  # output the GB
  GB <- particle[bestset]
  answer <- list("bestset" = bestset, "result" = result, "GB" = GB)
  return(answer)
}
