# Load data
# load("data/validation.RData");
data("validation", package = "kscons", envir = environment())

# ====================================================================================================
blocklize = function(rho) {
  n <- length(rho);
  if (n > 1) {
    K <- 1;
    for (i in 2:n) {
      if (rho[i] != rho[i - 1]) {
        K <- K + 1;
      }
    }
    eta <- rep(NA, K);
    lambda <- rep(0L, K);
    eta[1] <- as.character(rho[1]);
    temp <- 1;
    k <- 2;
    for (i in 2:n) {
      if (rho[i] != rho[i - 1]) {
        lambda[k - 1] <- i - temp;
        eta[k] <- as.character(rho[i]);
        temp <- i;
        k <- k + 1;
      }
    }
    lambda[K] <- i - temp + 1;
  } else {
    K <- 1;
    eta <- rho;
    lambda <- 1;
  }

  return (list(K = K, eta = eta, lambda = lambda));
}

# ====================================================================================================
reparametrize_rho = function(rho) {
  lambda <- blocklize(rho)$lambda;
  K <- blocklize(rho)$K;
  eta <- blocklize(rho)$eta;
  rho_2 <- NULL;
  count <- 1;
  for (k in 1:K) {
    if (eta[k] == "_") {
      rho_2 <- c(rho_2, rep(eta[k], lambda[k]));
    } else {
      rho_2 <- c(rho_2, rep(paste0(eta[k], count), lambda[k]));
      count <- count + 1;
    }
  }
  return (rho_2);
}

# ====================================================================================================
valid_socket_1 = function(rho) {
  L <- length(rho);
  uvw <- combn(1:L, 3);
  K <- dim(uvw)[2];
  gamma <- rep(2, K);
  for (k in 1:K) {
    if (sum(rho[uvw[, k]] == "_") > 0) {
      gamma[k] <- 0;
      next;
    }
    cond_1 <- uvw[2, k] == uvw[1, k] + 1 & uvw[3, k] == uvw[1, k] + 4 & sum(rho[uvw[1, k]:uvw[3, k]] != "H") == 0;
    cond_2 <- uvw[2, k] == uvw[1, k] + 3 & uvw[3, k] == uvw[1, k] + 4 & sum(rho[uvw[1, k]:uvw[3, k]] != "H") == 0;
    if (cond_1 | cond_2) {
      gamma[k] <- 1;
      next;
    }
    if (kscons::validation_local[paste0(rho[uvw[1, k]], rho[uvw[2, k]], rho[uvw[3, k]]), paste0(uvw[2, k] - uvw[1, k], "_", uvw[3, k] - uvw[2, k])] != 1) {
      gamma[k] <- 0;
    }
  }
  return (list(x = uvw[, which(gamma == 2)], y = uvw[, which(gamma == 1)]));
}

# ====================================================================================================
valid_socket_2 = function(rho) {
  # rho <- unlist(strsplit(rho, split = ""));
  rho_2 <- rep(paste0(blocklize(rho)$eta, 1:blocklize(rho)$K), blocklize(rho)$lambda, each = 1);
  L <- length(rho);
  x <- NULL;
  y <- NULL;

  # Determine possible local sockets
  for (l in 1:(L - 5)) {
    rho_sub <- rho[l:(l + 5)];
    x <- cbind(x, valid_socket_1(rho_sub)$x + l - 1)
    y <- cbind(y, valid_socket_1(rho_sub)$y + l - 1)
  }

  # Determine possible non-local sockets (only consider EEE here)
  if (sum(rho == "E") != 0) {
    temp <- combn(which(rho == "E"), 3);
    gamma <- rep(2, dim(temp)[2])
    for (k in 1:dim(temp)[2]) {
      if (sum(rho[temp[, k]] == "_") > 0) {
        gamma[k] <- 0;
        next;
      }
      if (rho_2[temp[1, k]] == rho_2[temp[2, k]] & rho_2[temp[1, k]] == rho_2[temp[3, k]]) {
        gamma[k] <- 0;
      }
      if (!(temp[2, k] == temp[1, k] + 1 | temp[3, k] == temp[2, k] + 1 | temp[2, k] == temp[1, k] + 2 | temp[3, k] == temp[2, k] + 2)) {
        gamma[k] <- 0;
      }
    }
    temp <- temp[, which(gamma == 2)];
    x <- cbind(x, temp);
  }

  return(list(x = t(unique(t(x))), y = t(unique(t(y)))));
}

# ====================================================================================================
loglikelihood = function(a_sub, rho_sub, gamma, local, model_local_true, model_local_false, model_nonlocal_true, model_nonlocal_false) {
  if (local) {
    if (gamma == 1) {
      return(model_local_true[paste0(rho_sub, collapse = ""), paste0(a_sub, collapse = "")]);
    } else {
      return(model_local_false[paste0(rho_sub, collapse = ""), paste0(a_sub, collapse = "")]);
    }
  } else {
    if (gamma == 1) {
      return(model_nonlocal_true[paste0(a_sub, collapse = "")]);
    } else {
      return(model_nonlocal_false[paste0(a_sub, collapse = "")]);
    }
  }
}

# ====================================================================================================
loglikelihood_2 = function(a_sub, rho_sub, delta, local, model_local_true_free, model_local_true_knob, model_nonlocal_true_free, model_nonlocal_true_knob) {
  if (local) {
    if (delta == 0) {
      return(model_local_true_free[paste0(rho_sub, collapse = ""), paste0(a_sub, collapse = "")]);
    } else {
      return(model_local_true_knob[paste0(rho_sub, collapse = ""), paste0(a_sub, collapse = "")]);
    }
  } else {
    if (delta == 0) {
      return(model_nonlocal_true_free[paste0(a_sub, collapse = "")]);
    } else {
      return(model_nonlocal_true_knob[paste0(a_sub, collapse = "")]);
    }
  }
}

# ====================================================================================================
list2contact = function(socket_all, gamma, L) {
  C <- matrix(0, nrow = L, ncol = L);
  for (j in which(gamma == 1)) {
    C[socket_all[1, j], socket_all[2, j]] <- 1;
    C[socket_all[2, j], socket_all[1, j]] <- 1;
    C[socket_all[1, j], socket_all[3, j]] <- 1;
    C[socket_all[3, j], socket_all[1, j]] <- 1;
    C[socket_all[2, j], socket_all[3, j]] <- 1;
    C[socket_all[3, j], socket_all[2, j]] <- 1;
  }
  return (C)
}

# ====================================================================================================
list2contact_2 = function(socket_all, z, gamma, delta, L) {
  C <- matrix(0, nrow = L, ncol = L);
  for (j in which(gamma == 1)) {
    C[socket_all[1, j], socket_all[2, j]] <- 1;
    C[socket_all[2, j], socket_all[1, j]] <- 1;
    C[socket_all[1, j], socket_all[3, j]] <- 1;
    C[socket_all[3, j], socket_all[1, j]] <- 1;
    C[socket_all[2, j], socket_all[3, j]] <- 1;
    C[socket_all[3, j], socket_all[2, j]] <- 1;
  }
  for (j in which(delta == 1)) {
    C[socket_all[1, j], z[j]] <- 1;
    C[socket_all[2, j], z[j]] <- 1;
    C[socket_all[3, j], z[j]] <- 1;
    C[z[j], socket_all[1, j]] <- 1;
    C[z[j], socket_all[2, j]] <- 1;
    C[z[j], socket_all[3, j]] <- 1;
  }
  return (C)
}

# ====================================================================================================
tabulate_error = function(gamma_true, gamma) {
  table = matrix(0L, 2, 2);
  p <- length(gamma_true);
  for (i in 1:p) {
    table[gamma[i] + 1, gamma_true[i] + 1] <- table[gamma[i] + 1, gamma_true[i] + 1] + 1;
  }
  return (table);
}

mcc = function(tabulate) {
  tp <- tabulate[2, 2];
  tn <- tabulate[1, 1];
  fp <- tabulate[2, 1];
  fn <- tabulate[1, 2];
  mcc <- (tp*tn - fp*fn)/sqrt((tp + fp)*(tp + fn)*(tn + fp)*(tn + fn));
  return (mcc);
}

# ====================================================================================================

