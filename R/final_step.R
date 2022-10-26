#' Compare MA consistency
#'
#' @param MA_df MA_df
compare_MA_df <- function(MA_df){
  # compare by column
  bestset <- 1:nrow(MA_df)
  row_name <- row.names(MA_df)
  for(i in 1:ncol(MA_df)){
    mini <- min(MA_df[bestset, i])
    row_name <- row_name[which(MA_df[bestset, i] == mini)]
    bestset <- as.numeric(row_name)
  }
  return(bestset)
}



final_step <- function(SIB_result, SIB_time, P_w, total_unit, unit,
                       factor_number){
  # collect GB's U in multiple stratum
  stratum_GB <- lapply(SIB_result, function(i){
    ans <- sapply(1:SIB_time, function(j){
      return(i$testSIB[[j]]$GB)
    })
    return(unlist(ans, recursive = FALSE))
  })

  # collect GB's U in multiple stratum
  stratum_U <- lapply(SIB_result, function(i){
    ans <- sapply(1:SIB_time, function(j){
      return(i$testSIB[[j]]$U)
    })
    return(unlist(ans, recursive = FALSE))
  })
  stratum_MA <- lapply(stratum_U, function(i){
    MA_under_P_w <- lapply(P_w, function(j){
      MA <- lapply(i, function(k){
        wordlength_list <- lapply(1:length(j), function(z){
          return(worlen_pattern(factor_number = factor_number,
                                U = k,
                                P_w = j[[z]],
                                total_unit = total_unit,
                                unit = unit))
        })
        wordlength <- Reduce("+", wordlength_list)
        return(wordlength)
      })
      MA <- do.call(rbind.data.frame, MA)
      colnames(MA) <- paste0("X", 1:factor_number)
      return(MA)
    })
    return(MA_under_P_w)
  })
  intersection <- lapply(stratum_MA, function(i){
    comparison <- lapply(i, function(j){
      return(compare_MA_df(MA_df = j))
    })
    return(Reduce(intersect, comparison))
  })
  # random sample an index from intersection
  selected_index <- lapply(intersection, function(i){
    return(i[sample(1:length(i), 1)])
  })
  selected_GB <- lapply(1:length(selected_index), function(i){
    return(stratum_GB[[i]][selected_index[[i]]])
  })
  selected_MA <- lapply(1:length(selected_index), function(i){
    MA <- lapply(1:length(stratum_MA[[i]]), function(j){
      return(stratum_MA[[i]][[j]][selected_index[[i]], ])
    })
    MA <- do.call(rbind.data.frame, MA)
    rownames(MA) <- NULL
    return(MA)
  })
  consistency <- all_equal_list(selected_MA)
  result <- list("selected_GB" = selected_GB, "selected_MA" = selected_MA,
                 "consistency" = consistency)
}
