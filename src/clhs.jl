using Iterators
using StatsBase
using DataFrames

function metropolis_value(delta, temp)
  return exp(-1*delta/temp)
end

function edge_col(x, size)
  # Initiate array storing result
  res = Float64[]
  # For each bin in the strata
  for q in linspace(0, 1, size + 1)
    push!(res, quantile(x, q))
  end
  # Return 1D-Array
  return res
end

function compute_edges(df, size)
  # Pre-allocating result Array
  res = Array(Float64, size + 1, ncol(df))
  # For each column
  for col in 1:ncol(df)
    res[:,col] = edge_col(df[col], size)
  end

  return res
end

function objective(
  size,
  data_continuous_sampled,
  continuous_strata,
  cor_mat
  )

  # Continuous variables
  #
  n_cont_variables <- ncol(data_continuous_sampled)

  cont_data_strata <- lapply(1:n_cont_variables, function(i) list(data_continuous_sampled[, i], continuous_strata[, i]) )
  cont_obj_sampled <- lapply(cont_data_strata, function(x) hist(x[[1]], breaks = x[[2]], plot = FALSE)$counts)
  cont_obj_sampled <- matrix(unlist(cont_obj_sampled), ncol = n_cont_variables, byrow = FALSE)

  delta_obj_continuous <- rowSums(abs(cont_obj_sampled - 1))

  # Correlation of continuous data
  #
  cor_sampled <- cor(data_continuous_sampled)
  cor_sampled[is.na(cor_sampled)] <- 1 # when there's only one observation, cor() throws NAs - we set these to 1

  delta_obj_cor <- sum(abs(cor_mat - cor_sampled))

  # Objective function
  #
  obj <- sum(delta_obj_continuous) + delta_obj_cor

  # Returning results
  #
  list(obj = obj, delta_obj_continuous = delta_obj_continuous, delta_obj_cor = delta_obj_cor)
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

  # Initiate Metropolis value
  metropolis = metropolis_value(0, temp)
  # Number of individuals in the data set
  n_data = nrow(x)

  # Edge of the strata
  continuous_strata = compute_edges(x, size)

  # Data correlation
  cor_mat = cor(x)

  # initialise, pick randomly
  n_remainings = n_data - size # number of individuals remaining unsampled

  i_sampled = sample(1:nrow(x), size, false) # individuals (rows) randomly chosen
  i_unsampled = find( [in(i, i_sampled) for i in 1:n_data] )

  data_continuous = x

  data_continuous_sampled = data_continuous[i_sampled,:] # sampled continuous data

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

    if rand() < 0.5
      # pick a random sampled point and random unsampled point
      idx_unsampled = sample(1:n_remainings)
      idx_sampled = sample(1:size)
      # Swap these:
      i_sampled[idx_sampled] <- i_unsampled[idx_unsampled]
      i_unsampled[idx_unsampled] <- i_sampled[idx_sampled]

      # creating new data sampled
      data_continuous_sampled <- data_continuous[i_sampled, , drop = FALSE]
    else
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
    end

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
