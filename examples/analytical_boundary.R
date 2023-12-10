analytical_boundary_Gaussian = function(cutoff_out,
                                        cutoff_inner,
                                        eta,
                                        lift,
                                        arc_width){
  
  W_upper_bound = -log(cutoff_out)/eta
  W_lower_bound = -log(1-cutoff_inner)/eta
  
  out_starting_angle = asin(lift/W_upper_bound)
  out_angle_seq = seq(from = out_starting_angle, to = pi - out_starting_angle, by = arc_width/W_upper_bound)
  out_angle_seq = union(out_angle_seq, pi - out_starting_angle)
  
  inner_starting_angle = asin(lift/W_lower_bound)
  inner_angle_seq = seq(from = inner_starting_angle, to = pi - inner_starting_angle, by = arc_width/W_lower_bound)
  inner_angle_seq = union(inner_angle_seq, pi - inner_starting_angle)
  
  out_x_seq = W_upper_bound * cos(out_angle_seq)
  out_y_seq = W_upper_bound * sin(out_angle_seq)
  out_seq = rbind(out_x_seq, out_y_seq)
  out_seq = t(out_seq)
  
  inner_x_seq = W_lower_bound * cos(inner_angle_seq)
  inner_y_seq = W_lower_bound * sin(inner_angle_seq)
  inner_x_seq = rev(inner_x_seq)
  inner_y_seq = rev(inner_y_seq)
  inner_seq = rbind(inner_x_seq, inner_y_seq)
  inner_seq = t(inner_seq)
  
  left_seq = c( (out_x_seq[length(out_x_seq)] + inner_x_seq[1])/2, lift )
  #left_x_seq = seq(from = out_x_seq[length(out_x_seq)], to = inner_x_seq[length(inner_x_seq)], by = arc_width)
  #left_y_seq = 0*left_x_seq + lift
  #left_seq = rbind(left_x_seq, left_y_seq)
  #left_seq = t(left_seq)
  
  right_seq = c( (out_x_seq[1] + inner_x_seq[length(inner_x_seq)])/2, lift )
  
  # right_x_seq = seq(from = inner_x_seq[length(inner_x_seq)], to = out_x_seq[1], by = arc_width)
  # right_y_seq = 0*right_x_seq + lift
  # right_seq = rbind(right_x_seq, right_y_seq)
  # right_seq = t(right_seq)
  
  boundary_path = rbind(out_seq, left_seq,inner_seq, right_seq)
  
  return(boundary_path)
}

analytical_boundary_GP = function(cutoff,
                                  eta,
                                  bound_away_from_one
){
  W_upper_bound = -log(cutoff)/eta
  # W_lower_bound = -log(1-cutoff_inner)/eta
  boundary_path = c()
  boundary_path = rbind(boundary_path, c(0,1 - bound_away_from_one), c(0,0), c(W_upper_bound, 0), c( bound_away_from_one * W_upper_bound, 1 - bound_away_from_one ) )
  
  
  return(boundary_path)
}
