function metropolis(delta, temp)
  return exp(-1*delta/temp)
end

function clhs(
  x, # a DataFrame
  size, # Number of samples you want
  iter = 10000, # Number of max iterations
  temp = 1, # initial temperature
  tdecrease = 0.95, # temperature decrease rate
  obj.limit = -Inf, # Stopping criterion
  length.cycle = 10 # Number of cycles done at each constant temperature value
  )

  metropolis = metropolis(0, temp) # Initial Metropolis value
  n_data = nrow(x) # Number of individuals in the data set

  # Edge of the strata
  continuous_strata <- apply(
    data_continuous,
    2,
    function(x) {
      quantile(x, probs = seq(0, 1, length.out = size + 1), na.rm = TRUE)
    }
  )

  # Data correlation
  cor_mat <- cor(data_continuous, use = "complete.obs")

  # initialise, pick randomly
  n_remainings <- n_data - size # number of individuals remaining unsampled
  i_sampled <- sample(1:n_data, size = size, replace = FALSE) # individuals randomly chosen
  i_unsampled <- setdiff(1:n_data, i_sampled) # individuals remaining unsampled
  data_continuous_sampled <- data_continuous[i_sampled, , drop = FALSE] # sampled continuous data

  # objective function
  res <- .lhs_obj(size = size, data_continuous_sampled = data_continuous_sampled, data_factor_sampled = data_factor_sampled, continuous_strata = continuous_strata, cor_mat = cor_mat, factor_obj = factor_obj, weights = weights)

  obj <- res$obj # value of the objective function
  delta_obj_continuous <- res$delta_obj_continuous

  # vector storing the values of the objective function
  obj_values <- vector(mode = 'numeric', length = iter)

  for (i in 1:iter) {

    # storing previous values
    previous <- list()
    previous$obj <- obj
    previous$i_sampled <- i_sampled
    previous$i_unsampled <- i_unsampled
    previous$delta_obj_continuous <- delta_obj_continuous

    if (runif(1) < 0.5) {
      # pick a random sampled point and random unsampled point
      idx_unsampled <- sample(1:n_remainings, size = 1)
      idx_sampled <- sample(1:size, size = 1)
      # Swap these:
      i_sampled[idx_sampled] <- i_unsampled[idx_unsampled]
      i_unsampled[idx_unsampled] <- i_sampled[idx_sampled]

      # creating new data sampled
      data_continuous_sampled <- data_continuous[i_sampled, , drop = FALSE]
    }
    else {
      # remove the worse sampled & resample
      worse <- max(delta_obj_continuous)
      i_worse <- which(delta_obj_continuous == worse)
      n_worse <- length(i_worse)

      # swap with reservoir
      spl_removed <- i_sampled[i_worse] # will be removed from the sampled set
      idx_added <- sample(1:n_remainings, size = n_worse) # will take their place
      i_sampled[i_worse] <- i_unsampled[idx_added] # replacing worst sampled by new pick
      i_unsampled[1:n_worse] <- spl_removed # replacing the worst pick in the reservoir

      # creating new data sampled
      data_continuous_sampled <- data_continuous[i_sampled, , drop = FALSE]
    }

    # calc obj
    res <- .lhs_obj(size = size, data_continuous_sampled = data_continuous_sampled, data_factor_sampled = data_factor_sampled, continuous_strata = continuous_strata, cor_mat = cor_mat, factor_obj = factor_obj, weights = weights)

    obj <- res$obj
    delta_obj_continuous <- res$delta_obj_continuous
    # Compare with previous iterations
    delta_obj <- obj - previous$obj
    metropolis <- exp(-1*delta_obj/temp) #+ runif(1)*temp
    metropolis_cost <- Inf # runif(1) >= Inf is always FALSE

    # If the optimum has been reached
    if (obj <= obj.limit) {
      warning("\nThe objective function has reached its minimum value, as specified by the obj.limit option.")
      obj_values[i] <- obj
      break
    }

    # Revert change
    if (delta_obj > 0 & runif(1) >= metropolis | runif(1) >= metropolis_cost) {
      i_sampled <- previous$i_sampled
      i_unsampled <- previous$i_unsampled
      data_continuous_sampled <- data_continuous[i_sampled, , drop = FALSE]

      obj <- previous$obj
      delta_obj_continuous <- previous$delta_obj_continuous
    }

    # Storing the objective function value of the current iteration
    obj_values[i] <- obj

    # Temperature decrease
    if ((i %% length.cycle) == 0) temp <- temp*tdecrease
  }

  sampled_data <- data_continuous_sampled

  # Simple output - just the sampled object
  if (simple) res <- i_sampled
  else {
    # Making up the object to be returned
    res <- list(
      initial_object = x,
      index_samples = i_sampled,
      sampled_data = sampled_data,
      obj = obj_values,
      cost = op_cost_values
    )
    class(res) = c("cLHS_result","list")
  }

  return res

end
