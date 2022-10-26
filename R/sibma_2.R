sibma <- function(factor_level, unit, particle_number, particle_increase,
                  SIB_time, all_two_level, P_w, q_GB, q_LB, q_new,
                  t, total_unit, multiply_len, t_increase, incidence_matrix,
                  treatment_effect){
  candidate <- lapply(1:length(P_w), function(i){
    return(SIB_algorithm(factor_level = factor_level, unit = unit,
                         particle_number = particle_number,
                         particle_increase = particle_increase,
                         SIB_time = SIB_time,
                         all_two_level = all_two_level,
                         P_w = P_w[[i]], q_GB = q_GB, q_LB = q_LB,
                         q_new = q_new, t = t, total_unit = total_unit,
                         multiply_len = multiply_len,
                         t_increase = t_increase,
                         incidence_matrix = incidence_matrix[[i]],
                         treatment_effect = treatment_effect[[i]]))
  })

  if(length(P_w) > 1){
    final_step <- final_step(SIB_result = candidate, SIB_time = SIB_time,
                              P_w = P_w, total_unit = total_unit, unit = unit,
                              factor_number = length(factor_level))
    result <- list("candidate" = candidate, "final_step" = final_step)
  }else{
    result <- list("candidate" = candidate)
  }

  return(result)
}
