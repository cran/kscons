# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

estimate <- function(a, rho, iter, burn, msa = NULL) {
  #source('R/functions.R')
  if (nchar(a) != nchar(rho)) {
    stop("Amino acid sequence and secondary structure are not at the same length!");
  }
  if (burn > iter) {
    stop("Burn-in cannot exceed iteration")
  }
  a <- unlist(strsplit(a, ""));
  rho <- unlist(strsplit(rho, ""));
  alpha <- 1;
  omega <- 0.05;
  start_gamma_sum <- 10;

  # Load sampling model
  #load("data/lklh.RData");
  data("lklh", package = "kscons", envir = environment())
  model_local_true_free <- kscons::allmodel_local_true_free;
  model_local_true_knob <- kscons::allmodel_local_true_knob;
  model_local_true <- model_local_true_free + model_local_true_knob;
  model_local_false <- kscons::allmodel_local_false;

  model_nonlocal_true_free <- kscons::allmodel_nonlocal_true_free;
  model_nonlocal_true_knob <- kscons::allmodel_nonlocal_true_knob;
  model_nonlocal_true <- model_nonlocal_true_free + model_nonlocal_true_knob;
  model_nonlocal_false <- kscons::allmodel_nonlocal_false;

  model_local_true_free <- log((model_local_true_free + 1) / rowSums(model_local_true_free + 1));
  model_local_true_knob <- log((model_local_true_knob + 1) / rowSums(model_local_true_knob + 1));
  model_local_true <- log((model_local_true + 1) / rowSums(model_local_true + 1));
  model_local_false <- log((model_local_false + 1) / rowSums(model_local_false + 1));

  model_nonlocal_true_free <- log((model_nonlocal_true_free + 1) / sum(model_nonlocal_true_free + 1));
  model_nonlocal_true_knob <- log((model_nonlocal_true_knob + 1) / sum(model_nonlocal_true_knob + 1));
  model_nonlocal_true <- log((model_nonlocal_true + 1) / sum(model_nonlocal_true + 1));
  model_nonlocal_false <- log((model_nonlocal_false + 1) / sum(model_nonlocal_false + 1));

  # Load prior
  if (!is.null(msa)) {
    if(!isSymmetric.matrix(msa)) {
      stop("Multiple structure alignment matrix is not symmetric!")
    }
    n <- max(msa) + 0.001;
    phi <- msa + alpha;
    phi <- phi/(n + alpha);
  }

  # obtain the basic information
  L <- length(a);
  LL <- which(a != "_");
  rho_2 <- reparametrize_rho(rho);
  socket_all <- valid_socket_2(rho)$x;
  K <- dim(socket_all)[2];
  socket_hhh <- valid_socket_2(rho)$y;
  K_hhh <- dim(socket_hhh)[2];
  socket_allhhh <- cbind(socket_all, socket_hhh);
  local <- socket_allhhh[3,] - socket_allhhh[1,] <= 5;

  # Initialization
  E <- 5;
  gamma <- rep(0, K);
  if (!is.null(msa)) {
    gamma[sample.int(K, start_gamma_sum, replace = FALSE)] <- 1;
    if (K_hhh > 0) {
      gamma <- c(gamma, rep(1, K_hhh));
    }
    delta <- rep(0, length(gamma));
    delta[sample(which(gamma == 1), start_gamma_sum/2)] <- 1;
    z <- rep(0, length(gamma));
    for (i in which(delta == 1)) {
      z[i] <- sample(setdiff(1:L, which(rho_2 %in% unique(rho_2[socket_allhhh[, i]]))), 1);
    }
  } else {
    delta <- rep(0, length(gamma));
    z <- rep(0, length(gamma));
  }
  mpv <- rep(0, K);
  delta_mpv <- rep(0, length(delta));
  mpm <- matrix(0, nrow = L, ncol = L);
  z_map <- z;
  delta_map <- delta;

  start_time <- proc.time();
  # MCMC updates
  gamma_map <- gamma;
  for (i in 1:(burn + iter)) {

    # Between-model updates
    temp <- sample.int(3, 1);
    if (temp == 1) {
      # Add
      if (sum(gamma[1:K]) != K) {
        k <- sample(which(gamma[1:K] == 0), 1);
        if (!(a[socket_all[1, k]] == "Z" | a[socket_all[1, k]] == "X" | a[socket_all[2, k]] == "Z" | a[socket_all[2, k]] == "X" | a[socket_all[3, k]] == "Z" | a[socket_all[3, k]] == "X")) {
          posterior_ratio <- loglikelihood(a[socket_all[, k]], rho[socket_all[, k]], 1, local[k], model_local_true, model_local_false, model_nonlocal_true, model_nonlocal_false) - loglikelihood(a[socket_all[, k]], rho[socket_all[, k]], 0, local[k], model_local_true, model_local_false, model_nonlocal_true, model_nonlocal_false);
          if (is.null(msa)) {
            posterior_ratio <- posterior_ratio + log(omega) - log(1 - omega);
          } else {
            gamma_temp <- gamma[1:K];
            gamma_temp[k] <- 1;
            contact <- list2contact(socket_all, gamma[1:K], L);
            contact_temp <- list2contact(socket_all, gamma_temp, L);
            index_diff <- which(contact_temp != contact);
            posterior_ratio <- posterior_ratio + sum(log(phi[index_diff]))/2 - sum(log(1 - phi[index_diff]))/2;
          }
          hastings <- posterior_ratio + log(K - sum(gamma[1:K])) - log(sum(gamma[1:K]) + 1);
          if (hastings >= log(runif(1))) {
            gamma[k] <- 1;
            if (i > burn & posterior_ratio >= 0) {
              gamma_map <- gamma;
              if (!is.null(msa)) {
                delta_map <- delta;
                z_map <- z;
              }
            }
          }
        }
      }
    } else if (temp == 2) {
      # Delete
      if (sum(gamma[1:K]) > 1) {
        k <- sample(which(gamma[1:K] == 1), 1);
        if (!(a[socket_all[1, k]] == "Z" | a[socket_all[1, k]] == "X" | a[socket_all[2, k]] == "Z" | a[socket_all[2, k]] == "X" | a[socket_all[3, k]] == "Z" | a[socket_all[3, k]] == "X")) {
          posterior_ratio <- loglikelihood(a[socket_all[, k]], rho[socket_all[, k]], 0, local[k], model_local_true, model_local_false, model_nonlocal_true, model_nonlocal_false) - loglikelihood(a[socket_all[, k]], rho[socket_all[, k]], 1, local[k], model_local_true, model_local_false, model_nonlocal_true, model_nonlocal_false);
          if (is.null(msa)) {
            posterior_ratio <- posterior_ratio + log(1 - omega) - log(omega);
          } else {
            gamma_temp <- gamma[1:K];
            gamma_temp[k] <- 0;
            contact <- list2contact(socket_all, gamma[1:K], L);
            contact_temp <- list2contact(socket_all, gamma_temp, L);
            index_diff <- which(contact_temp != contact);
            posterior_ratio <- posterior_ratio + sum(log(1 - phi[index_diff]))/2 - sum(log(phi[index_diff]))/2;
          }
          hastings <- posterior_ratio + log(sum(gamma[1:K])) - log(K - sum(gamma[1:K]) + 1);
          if (hastings >= log(runif(1))) {
            gamma[k] <- 0;
            delta[k] <- 0;
            z[k] <- 0;
            if (i > burn & posterior_ratio >= 0) {
              gamma_map <- gamma;
              if (!is.null(msa)) {
                delta_map <- delta;
                z_map <- z;
              }
            }
          }
        }
      }
    } else {
      # Swap
      if (sum(gamma[1:K]) != 0 && sum(gamma[1:K]) != K) {
        k_0 <- sample(which(gamma[1:K] == 0), 1);
        k_1 <- sample(which(gamma[1:K] == 1), 1);
        if (!(a[socket_all[1, k_0]] == "Z" | a[socket_all[1, k_0]] == "X" | a[socket_all[2, k_0]] == "Z" | a[socket_all[2, k_0]] == "X" | a[socket_all[3, k_0]] == "Z" | a[socket_all[3, k_0]] == "X" | a[socket_all[1, k_1]] == "Z" | a[socket_all[1, k_1]] == "X" | a[socket_all[2, k_1]] == "Z" | a[socket_all[2, k_1]] == "X" | a[socket_all[3, k_1]] == "Z" | a[socket_all[3, k_1]] == "X")) {
          posterior_ratio <- loglikelihood(a[socket_all[, k_0]], rho[socket_all[, k_0]], 1, local[k_0], model_local_true, model_local_false, model_nonlocal_true, model_nonlocal_false) + loglikelihood(a[socket_all[, k_1]], rho[socket_all[, k_1]], 0, local[k_1], model_local_true, model_local_false, model_nonlocal_true, model_nonlocal_false) - loglikelihood(a[socket_all[, k_0]], rho[socket_all[, k_0]], 0, local[k_0], model_local_true, model_local_false, model_nonlocal_true, model_nonlocal_false) - loglikelihood(a[socket_all[, k_1]], rho[socket_all[, k_1]], 1, local[k_1], model_local_true, model_local_false, model_nonlocal_true, model_nonlocal_false);
          if (is.null(msa)) {

          } else {
            gamma_temp <- gamma[1:K];
            gamma_temp[k_0] <- 1;
            gamma_temp[k_1] <- 0;
            contact <- list2contact(socket_all, gamma[1:K], L);
            contact_temp <- list2contact(socket_all, gamma_temp, L);
            index_diff <- which(contact_temp != contact);
            posterior_ratio <- posterior_ratio + sum(contact_temp[index_diff]*log(phi[index_diff]) + (1 - contact_temp[index_diff])*(log(1 - phi[index_diff])))/2 - sum(contact[index_diff]*log(phi[index_diff]) + (1 - contact[index_diff])*(log(1 - phi[index_diff])))/2;
          }
          hastings <- posterior_ratio;
          if (hastings >= log(runif(1))) {
            gamma[k_0] <- 1;
            gamma[k_1] <- 0;
            delta[k_1] <- 0;
            z[k_1] <- 0;
            if (i > burn & posterior_ratio >= 0) {
              gamma_map <- gamma;
              if (!is.null(msa)) {
                delta_map <- delta;
                z_map <- z;
              }
            }
          }
        }
      }
    }

    if (!is.null(msa)) {
      # Within-model updates
      for (e in 1:E) {
        temp <- sample.int(2, 1);
        if (temp == 1) {
          # Add
          if (sum(delta) < sum(gamma)) {
            k <- sample(which(delta == 0 & gamma == 1), 1);
            if (length(setdiff(LL, which(rho_2 %in% unique(rho_2[socket_allhhh[, k]])))) > 0 & !(a[socket_allhhh[1, k]] == "Z" | a[socket_allhhh[1, k]] == "X" | a[socket_allhhh[2, k]] == "Z" | a[socket_allhhh[2, k]] == "X" | a[socket_allhhh[3, k]] == "Z" | a[socket_allhhh[3, k]] == "X")) {
              delta_temp <- delta;
              delta_temp[k] <- 1;
              posterior_ratio <- loglikelihood_2(a[socket_allhhh[, k]], rho[socket_allhhh[, k]], 1, local[k], model_local_true_free, model_local_true_knob, model_nonlocal_true_free, model_nonlocal_true_knob) - loglikelihood_2(a[socket_allhhh[, k]], rho[socket_allhhh[, k]], 0, local[k], model_local_true_free, model_local_true_knob, model_nonlocal_true_free, model_nonlocal_true_knob);
              z_temp <- z;
              z_temp[k] <- sample(setdiff(LL, which(rho_2 %in% unique(rho_2[socket_allhhh[, k]]))), 1);
              contact <- list2contact_2(socket_allhhh, z, gamma, delta, L);
              contact_temp <- list2contact_2(socket_allhhh, z_temp, gamma, delta_temp, L);
              index_diff <- which(contact_temp != contact);
              posterior_ratio <- posterior_ratio + sum(log(phi[index_diff]))/2 - sum(log(1 - phi[index_diff]))/2;
              hastings <- posterior_ratio + log(sum(gamma) - sum(delta)) - log(sum(delta) + 1);
              if (hastings >= log(runif(1))) {
                delta[k] <- 1;
                z[k] <- z_temp[k];
              }
            }
          }
        } else if (temp == 2) {
          # Delete
          if (sum(delta) > 1) {
            k <- sample(which(delta == 1), 1);
            if (!(a[socket_allhhh[1, k]] == "Z" | a[socket_allhhh[1, k]] == "X" | a[socket_allhhh[2, k]] == "Z" | a[socket_allhhh[2, k]] == "X" | a[socket_allhhh[3, k]] == "Z" | a[socket_allhhh[3, k]] == "X")) {
              delta_temp <- delta;
              delta_temp[k] <- 0;
              posterior_ratio <- loglikelihood_2(a[socket_allhhh[, k]], rho[socket_allhhh[, k]], 0, local[k], model_local_true_free, model_local_true_knob, model_nonlocal_true_free, model_nonlocal_true_knob) - loglikelihood_2(a[socket_allhhh[, k]], rho[socket_allhhh[, k]], 1, local[k], model_local_true_free, model_local_true_knob, model_nonlocal_true_free, model_nonlocal_true_knob);
              contact <- list2contact_2(socket_allhhh, z, gamma, delta, L);
              contact_temp <- list2contact_2(socket_allhhh, z, gamma, delta_temp, L);
              index_diff <- which(contact_temp != contact);
              posterior_ratio <- posterior_ratio + sum(log(1 - phi[index_diff]))/2 - sum(log(phi[index_diff]))/2;
              hastings <- posterior_ratio + log(sum(delta)) - log(sum(gamma) - sum(delta) + 1);
              if (hastings >= log(runif(1))) {
                delta[k] <- 0;
                z[k] <- 0;
              }
            }
          }
        }
      }
    }
    if (i > burn) {
      mpv <- mpv + gamma[1:K];
      delta_mpv <- delta_mpv + delta;
      contact <- list2contact_2(socket_allhhh, z, gamma, delta, L);
      mpm <- mpm + contact;
    }
  }

  mpv <- mpv/iter;
  delta_mpv <- delta_mpv/iter;
  mpm <- mpm/iter;
  if (!is.null(msa)) {
    contact_map <- list2contact_2(socket_allhhh, z_map, gamma_map, delta_map, L);
  } else {
    contact_map <- list2contact_2(socket_allhhh, z, gamma_map, delta, L);
  }
  temp <- valid_socket_2(rho)$y;
  if (length(temp) > 0) {
    for (j in 1:dim(temp)[2]) {
      mpm[temp[1, j], temp[2, j]] <- 1;
      mpm[temp[2, j], temp[1, j]] <- 1;
      mpm[temp[1, j], temp[3, j]] <- 1;
      mpm[temp[3, j], temp[1, j]] <- 1;
      mpm[temp[2, j], temp[3, j]] <- 1;
      mpm[temp[3, j], temp[2, j]] <- 1;
      contact_map[temp[1, j], temp[2, j]] <- 1;
      contact_map[temp[2, j], temp[1, j]] <- 1;
      contact_map[temp[1, j], temp[3, j]] <- 1;
      contact_map[temp[3, j], temp[1, j]] <- 1;
      contact_map[temp[2, j], temp[3, j]] <- 1;
      contact_map[temp[3, j], temp[2, j]] <- 1;
    }
  }
  if (length(mpv) != length(delta_mpv)) {
    mpv <- c(mpv, rep(1, length(delta_mpv) - length(mpv)))
  }
  end_time <- proc.time();
  time <- end_time - start_time;

  if (is.null(msa)) {
    return (list(time = time, a = a, rho = rho, socket = socket_allhhh, gamma_map = gamma_map, contact_map = contact_map, mpv = mpv, delta_mpv = delta_mpv, mpm = mpm));
    } else {
      return (list(time = time, a = a, rho = rho, socket = socket_allhhh, gamma_map = gamma_map, delta_map = delta_map, z_map = z_map, contact_map = contact_map, mpv = mpv, delta_mpv = delta_mpv, mpm = mpm));
    }
}
