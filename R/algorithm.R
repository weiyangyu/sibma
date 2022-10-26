#' Check All Elements in a List are Identical
#' @param a_list a list.
#' @importFrom utils combn
all_equal_list <- function(a_list){
  total <- combn(length(a_list), 2)
  result <- lapply(1:ncol(total), function(i){
    ind <- total[,i]
    return(identical(a_list[[ind[1]]], a_list[[ind[2]]]))
  })
  ans <- all(result == TRUE)
  return(ans)
}

#' Return the Order of a Data Frame of Word-Length Patterns
#' @param MA_df a data frame containing word-length patterns.
MA_ordering <- function(MA_df){
  n_row <- nrow(MA_df)
  fullset <- 1:n_row
  bestset <- 1:n_row
  ordering <- c()

  repeat{
    if(length(bestset) == 1)
      break

    for(i in 1:ncol(MA_df)){
      mini <- min(MA_df[bestset, i])
      bestset <- bestset[which(MA_df[bestset, i] == mini)]
    }
    ordering <- append(ordering, bestset)
    bestset <- fullset[-match(ordering, fullset)]
  }
  ordering <- append(ordering, bestset)

  return(ordering)
}

#' Check Convergence
#' Check if the word-length patterns calculated from each SIB process are the same.
#' @param MA_df a list. Each element is a list describing the word-length patterns in
#' each SIB process.
check_convergence <- function(MA_df){
  # calculate the number of GB in each SIB_time
  GB_number <- lapply(MA_df, function(i){
    nrow(i[[1]])
  })

  if(all_equal_list(GB_number) == TRUE){
    # reorder every stratum's word-length patterns
    MA_df <- lapply(MA_df, function(i){
      lapply(i, function(j){
        ordering <- MA_ordering(j)
        return(j[ordering, ])
      })
    })

    P_w_number <- length(MA_df[[1]])
    checking <- lapply(1:P_w_number, function(i){
      temp <- lapply(1:length(MA_df), function(j){
        return(MA_df[[j]][[i]])
      })
      return(all_equal_list(temp))
    })
    result <- all(unlist(checking))

  }else{
    result <- FALSE
  }
  return(result)
}

#' SIB Algorithm
#' @param factor_level a list. Each element is a numeric vector specifying
#' levels of a factor.
#' @param unit an integer describing the number of experimental units to be
#' used in a particle.
#' @param particle_number an integer indicating the initial particle number. Default is 10.
#' @param particle_increase an integer indicating how many particle number
#' you would like to increase in a sequential way. Default is 10.
#' @param SIB_time an integer indicating the number of SIB process you would
#' like to run. Default is 3.
#' @param all_two_level logical. If all factors are two levels,
#' \code{all_two_level} should be \code{TRUE}; otherwise it should be
#' \code{FALSE}.
#' @param P_w a nested list. The larger list is composed of smaller list in which
#' each element is a matrix describing the orthogonal projection matrix onto
#' the corresponding stratum variance.
#' @param q_GB an integer describing how many columns of a particle
#' should be mixed with the corresponding columns of \code{GB}. Default is 1.
#' @param q_LB an integer describing how many columns of a particle
#' should be mixed with the corresponding columns of \code{LB}. Default is 1.
#' @param q_new an integer describing how many columns of a particle
#' should be mixed with the corresponding columns of a random
#' particle created by \code{sibma::create_particle()}. Default is 1.
#' @param t an integer indicating the initial iteration number. Default is 10.
#' @param total_unit an integer representing the number of total run size in
#' a full factorial design.
#' @param multiply_len a numeric vector. Each element is used to modify
#' the value of each column in a model matrix.
#' @param t_increase an integer indicating the increased iteration number
#' you would like to increase in a sequential way. Default is 10.
#' @param incidence_matrix a list. Each element is a matrix specifying an
#' incidence matrix.
#' @param treatment_effect a list. Each element is a numeric vector specifying
#' which factors should have the same treatment effect in terms of a certain
#' incidence matrix. The order of elements in \code{treatment_effect} should be
#' corresponding to that in \code{incidence_matrix}.
#' @examples
#' # For a design with 4 two-level factors, we only want 8 experimental units.
#' # Besides, we want the design has row-column structure, i.e. we want four experimental
#' # units have a kind of treatment generated from the first and third factors,
#' # so the other four experimental units have another setting.
#' # Moreover, we also want four experimental units have a kind of treatment
#' # generated from the second and fourth factors,
#' # so the other four experimental units have another setting.
#'
#' F1 <- matrix(rep(1,8), ncol = 1)  # unstructured incidence matrix
#' P_w_1 <- F1 %*% (solve(t(F1) %*% F1)) %*% t(F1) # unstructured orthogonal projection matrix
#'
#' F2 <- matrix(c(1,1,1,1,0,0,0,0, 0,0,0,0,1,1,1,1), ncol = 2)  # first stratum's incidence matrix
#' P_v_2 <- F2 %*% (solve(t(F2) %*% F2)) %*% t(F2)
#' P_w_2 <- P_v_2 - P_w_1  # first stratum's orthogonal projection matrix
#'
#' F3 <- matrix(c(1,1,0,0,1,1,0,0, 0,0,1,1,0,0,1,1), ncol = 2)  # second stratum's incidence matrix
#' P_v_3 <- F3 %*% (solve(t(F3) %*% F3)) %*% t(F3)
#' P_w_3 <- P_v_3 - P_w_1  # first stratum's orthogonal projection matrix
#'
#' P_w1 <- list(P_w_1)
#' P_w2 <- list(P_w_1, P_w_2)
#' P_w3 <- list(P_w_1, P_w_3)
#' P_w4 <- list(P_w_1, P_w_2, P_w_3)
#'
#' P_w <- list(P_w1, P_w2, P_w3, P_w4)
#'
#' SIB_algorithm(factor_level = list(c(-1,1), c(-1,1), c(-1,1), c(-1,1)),
#'               unit = 8, particle_number = 10,
#'               particle_increase = 10,
#'               SIB_time = 3, all_two_level = TRUE, P_w = P_w, q_GB = 1, q_LB = 1,
#'               q_new = 1, t = 10, total_unit = 16, multiply_len = rep(1/sqrt(2)^4, 16),
#'               t_increase = 10, incidence_matrix = list(F2, F3),
#'               treatment_effect = list(c(1,3), c(2,4)))
#'
#' @export
SIB_algorithm <- function(factor_level, unit, particle_number = 10,
                          particle_increase = 10,
                          SIB_time = 3, all_two_level, P_w, q_GB = 1, q_LB = 1,
                          q_new = 1, t = 10, total_unit, multiply_len,
                          t_increase = 10, incidence_matrix, treatment_effect){
  ##################
  a <- 0
  b <- 0 ## the number of convergence occurred
  ##################

  ## run the first SIB process
  SXTY_SIB <- sxty_SIB(factor_level = factor_level, unit = unit,
                       particle_number = particle_number, SIB_time = SIB_time,
                       all_two_level = all_two_level, P_w = P_w, q_GB = q_GB,
                       q_LB = q_LB, q_new = q_new, t = t, total_unit = total_unit,
                       multiply_len = multiply_len,
                       incidence_matrix = incidence_matrix,
                       treatment_effect = treatment_effect)
  current_SIB <- SXTY_SIB

  repeat{
    ########################## check the results are converged or not
    while(check_convergence(MA_df = current_SIB$MA_df) != TRUE){
      SXTY_SIB <- sxty_SIB_relay_only_t(prior_SIB = current_SIB$testSIB,
                                        history = current_SIB$history,
                                        factor_level = factor_level,
                                        SIB_time = SIB_time,
                                        all_two_level = all_two_level,
                                        P_w = P_w,
                                        q_GB = q_GB, q_LB = q_LB, q_new = q_new,
                                        t = t_increase,
                                        total_unit = total_unit,
                                        multiply_len = multiply_len,
                                        incidence_matrix = incidence_matrix,
                                        treatment_effect = treatment_effect)
      # reorder
      reorder_current <- lapply(current_SIB$MA_df, function(i){
        lapply(i, function(j){
          ordering <- MA_ordering(j)
          return(j[ordering, ])
        })
      })
      reorder_sxty <- lapply(SXTY_SIB$MA_df, function(i){
        lapply(i, function(j){
          ordering <- MA_ordering(j)
          return(j[ordering, ])
        })
      })

      check_improvement <- lapply(1:SIB_time, function(i){
        temp <- identical(reorder_current[[i]], reorder_sxty[[i]])
      })

      ## update current status
      current_SIB <- SXTY_SIB

      if(all(unlist(check_improvement) == TRUE)){
        SXTY_SIB <- sxty_SIB_relay(prior_SIB = current_SIB$testSIB,
                                   history = current_SIB$history,
                                   factor_level = factor_level,
                                   unit = unit,
                                   particle_increase = particle_increase,
                                   SIB_time = SIB_time,
                                   all_two_level = all_two_level,
                                   P_w = P_w,
                                   q_GB = q_GB, q_LB = q_LB, q_new = q_new,
                                   t = t_increase,
                                   total_unit = total_unit,
                                   multiply_len = multiply_len,
                                   incidence_matrix = incidence_matrix,
                                   treatment_effect = treatment_effect)
        ## update current status
        current_SIB <- SXTY_SIB
      }
    }

    ###################
    b <- b + 1
    if(b < 2){
      pre_MA <- current_SIB$MA_df
    }else{
      post_MA <- current_SIB$MA_df
    }
    ###################

    if(b == 2){
      # reorder
      reorder_pre <- lapply(pre_MA[[1]], function(i){
        ordering <- MA_ordering(i)
        return(i[ordering, ])
      })
      reorder_post <- lapply(post_MA[[1]], function(i){
        ordering <- MA_ordering(i)
        return(i[ordering, ])
      })

      check_improvement <- identical(reorder_pre[[1]], reorder_post[[1]])

      if(check_improvement == TRUE){
        a <- 1
      }else{
        b <- 1
        pre_MA <- post_MA
      }
    }

    if(a == 1)
      break

    SXTY_SIB <- sxty_SIB_relay(prior_SIB = current_SIB$testSIB,
                               history = current_SIB$history,
                               factor_level = factor_level,
                               unit = unit,
                               particle_increase = particle_increase,
                               SIB_time = SIB_time,
                               all_two_level = all_two_level,
                               P_w = P_w,
                               q_GB = q_GB, q_LB = q_LB, q_new = q_new,
                               t = t_increase,
                               total_unit = total_unit,
                               multiply_len = multiply_len,
                               incidence_matrix = incidence_matrix,
                               treatment_effect = treatment_effect)
    ## update current status
    current_SIB <- SXTY_SIB
  }
  return(current_SIB)
}
