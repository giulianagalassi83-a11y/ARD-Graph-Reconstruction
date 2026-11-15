
# ============================================================
# STEP 1 - SYNTHETIC GRAPH GENERATION FROM ARD
# Village51 dataset
# Latent Surface Model + MCMC + Sparse Graph Simulation
# ============================================================

# -----------------------------
# 0. Packages
# -----------------------------
library(movMF)

# -----------------------------
# 1. Load input data
# -----------------------------

setwd("C:/Users/Giuliana Galassi/OneDrive/Documentos/IMDEA/Grafos Angela/Indian villages/Indian villages")

ard_path <- "ard_matrix_village51.csv"
total_prop_path <- "total_prop_village51.txt"

ard_data <- read.csv(ard_path, row.names = 1)
ard_data <- as.matrix(ard_data)

# Remove all-zero columns (e.g. groupNA)
zero_cols <- which(colSums(ard_data, na.rm = TRUE) == 0)

if (length(zero_cols) > 0) {
  cat("Removing all-zero ARD columns:\n")
  print(colnames(ard_data)[zero_cols])
  
  ard_data <- ard_data[, -zero_cols, drop = FALSE]
}

total.prop <- as.numeric(
  gsub("Total Proportion:\\s*", "", readLines(total_prop_path))
)

cat("ARD dimensions:\n")
print(dim(ard_data))
cat("total.prop:", total.prop, "\n")

if (ncol(ard_data) <= 4) {
  stop("The model requires more than 4 ARD groups after removing zero columns.")
}


# -----------------------------
# 2. Helper functions
# -----------------------------

rot.mat <- function(r.v) {
  phi <- r.v[1]
  theta <- r.v[2]
  psi <- r.v[3]
  
  rot.tmp <- matrix(NA, 3, 3)
  
  rot.tmp[1, ] <- c(
    cos(theta) * cos(psi),
    -cos(phi) * sin(psi) + sin(phi) * sin(theta) * cos(psi),
    sin(phi) * sin(theta) + cos(phi) * sin(theta) * cos(psi)
  )
  
  rot.tmp[2, ] <- c(
    cos(theta) * sin(psi),
    cos(phi) * cos(psi) + sin(phi) * sin(theta) * sin(psi),
    -sin(phi) * cos(psi) + cos(phi) * sin(theta) * sin(psi)
  )
  
  rot.tmp[3, ] <- c(
    -sin(theta),
    sin(phi) * cos(theta),
    cos(phi) * cos(theta)
  )
  
  return(rot.tmp)
}


sphere.update <- function(u.mean, u.k, m = 3, rot = FALSE, rot.vec = NULL) {
  
  b <- ((-2 * u.k) + sqrt((4 * u.k^2) + (m - 1)^2)) / (m - 1)
  x0 <- (1 - b) / (1 + b)
  
  z.tmp.u <- rbeta(length(u.k), (m - 1) / 2, (m - 1) / 2)
  u.tmp.u <- runif(length(u.k), 0, 1)
  
  w.u <- (1 - ((1 + b) * z.tmp.u)) /
    (1 - ((1 - b) * z.tmp.u))
  
  iter <- 0
  
  rejump.ind <- which(
    log(u.tmp.u) >
      ((u.k * w.u) + ((m - 1) * log(1 - (x0 * w.u))))
  )
  
  while (length(rejump.ind) > 0 && iter < 200) {
    
    z.tmp.u.tmp <- rbeta(
      length(rejump.ind),
      (m - 1) / 2,
      (m - 1) / 2
    )
    
    u.tmp.u.tmp <- runif(length(rejump.ind), 0, 1)
    
    w.u.tmp <- (
      1 - ((1 + b[rejump.ind]) * z.tmp.u.tmp)
    ) / (
      1 - ((1 - b[rejump.ind]) * z.tmp.u.tmp)
    )
    
    w.u[rejump.ind] <- w.u.tmp
    u.tmp.u[rejump.ind] <- u.tmp.u.tmp
    z.tmp.u[rejump.ind] <- z.tmp.u.tmp
    
    iter <- iter + 1
    
    rejump.ind <- which(
      log(u.tmp.u) >
        ((u.k * w.u) + ((m - 1) * log(1 - (x0 * w.u))))
    )
  }
  
  v.tmp.u <- matrix(
    runif(length(u.k) * (m - 1), -100, 100),
    length(u.k),
    m - 1
  )
  
  v.tmp.u <- v.tmp.u / sqrt(rowSums(v.tmp.u^2))
  
  xt <- cbind(
    sqrt(1 - w.u^2) * v.tmp.u,
    w.u
  )
  
  if (rot) {
    
    if (!is.null(rot.vec)) {
      return(rot.mat(rot.vec) %*% t(xt))
    }
    
    den <- sqrt(sum(u.mean[1:2]^2))
    
    r.tmp <- c(
      0,
      acos(u.mean[3]),
      acos(u.mean[1] / den)
    )
    
    if (u.mean[2] < 0) {
      r.tmp <- c(
        0,
        acos(u.mean[3]),
        -acos(u.mean[1] / den)
      )
    }
    
    return(rot.mat(r.tmp) %*% t(xt))
  }
  
  return(xt)
}


cp.fcn <- function(kap, p = 3) {
  
  top <- kap^((p / 2) - 1)
  
  bot <- (
    ((2 * pi)^(p / 2)) *
      besselI(kap, ((p / 2) - 1))
  )
  
  return(top / bot)
}


f.datalik <- function(
    data,
    v.alpha,
    v.beta,
    v.eta,
    v.etak,
    v.z,
    v.muk
) {
  
  v.theta <- v.z %*% t(v.muk)
  
  v.lambda <- exp(
    outer(v.alpha, v.beta, "+")
  ) *
    cp.fcn(v.eta) *
    cp.fcn(v.etak) /
    (
      cp.fcn(0.000001) *
        cp.fcn(
          sqrt(
            v.eta^2 +
              v.etak^2 +
              2 * v.eta * v.etak *
              cos(acos(v.theta))
          )
        )
    )
  
  prob.y <- (
    (-v.lambda) +
      (data * log(v.lambda)) -
      log(factorial(data))
  )
  
  return(prob.y)
}


generateRandomInitial <- function(n, p) {
  
  z <- matrix(
    rnorm(n * p),
    nrow = n,
    ncol = p
  )
  
  z <- sweep(
    z,
    MARGIN = 1,
    1 / sqrt(rowSums(z^2)),
    `*`
  )
  
  return(z)
}


# -----------------------------
# 3. MCMC update
# -----------------------------

f.update <- function(
    data,
    v.alpha,
    v.beta,
    v.eta,
    v.etak,
    v.z,
    v.muk,
    a.mu,
    a.s,
    b.mu,
    b.s,
    ad.jump,
    bd.jump,
    eta.jump,
    etak.jump,
    muk.jump,
    z.jump,
    k.mu,
    k.s,
    vjind,
    total.prop,
    muk.fix,
    muk.fix.ind,
    ls.dim,
    is.muk.fix = TRUE
) {
  
  data.n <- nrow(data)
  data.j <- ncol(data)
  
  # alpha
  a.star <- v.alpha + rnorm(length(v.alpha), 0, sqrt(ad.jump))
  a.old <- v.alpha
  
  lik.old <- f.datalik(data, a.old, v.beta, v.eta, v.etak, v.z, v.muk)
  lik.star <- f.datalik(data, a.star, v.beta, v.eta, v.etak, v.z, v.muk)
  
  prob.alpha.diff <- (
    rowSums(lik.star, na.rm = TRUE) +
      log(dnorm(a.star, a.mu, sqrt(a.s))) -
      rowSums(lik.old, na.rm = TRUE) -
      log(dnorm(a.old, a.mu, sqrt(a.s)))
  )
  
  alpha.valid <- !is.infinite(prob.alpha.diff) & !is.na(prob.alpha.diff)
  
  jump.new <- rep(0, length(v.alpha))
  
  jump.new[alpha.valid] <- rbinom(
    sum(alpha.valid),
    1,
    exp(pmin(prob.alpha.diff[alpha.valid], 0))
  )
  
  v.alpha <- ifelse(jump.new, a.star, a.old)
  
  alpha.p <- rep(0, length(v.alpha))
  
  alpha.p[alpha.valid] <- exp(
    pmin(prob.alpha.diff[alpha.valid], 0)
  )
  
  a.mu <- rnorm(
    1,
    mean(v.alpha, na.rm = TRUE),
    sqrt(a.s / data.n)
  )
  
  a.s <- sum(
    (v.alpha - a.mu)^2,
    na.rm = TRUE
  ) / rchisq(1, data.n - 1)
  
  # beta
  b.new <- v.beta + rnorm(length(v.beta), 0, sqrt(bd.jump))
  b.old <- v.beta
  
  lik.old <- f.datalik(data, v.alpha, b.old, v.eta, v.etak, v.z, v.muk)
  lik.new <- f.datalik(data, v.alpha, b.new, v.eta, v.etak, v.z, v.muk)
  
  prob.beta.diff <- (
    colSums(lik.new, na.rm = TRUE) +
      log(dnorm(b.new, b.mu, sqrt(b.s))) -
      colSums(lik.old, na.rm = TRUE) -
      log(dnorm(b.old, b.mu, sqrt(b.s)))
  )
  
  beta.valid <- !is.infinite(prob.beta.diff) & !is.na(prob.beta.diff)
  
  jump.new <- rep(0, length(v.beta))
  beta.p <- rep(0, length(v.beta))
  
  if (sum(beta.valid) > 0) {
    
    jump.new[beta.valid] <- rbinom(
      sum(beta.valid),
      1,
      exp(pmin(prob.beta.diff[beta.valid], 0))
    )
    
    v.beta <- jump.new * b.new + (1 - jump.new) * b.old
    
    beta.p[beta.valid] <- exp(
      pmin(prob.beta.diff[beta.valid], 0)
    )
    
  } else {
    v.beta <- b.old
  }
  
  b.mu <- rnorm(1, mean(v.beta), sqrt(b.s / data.j))
  
  b.s <- sum((v.beta - b.mu)^2) /
    rchisq(1, data.j - 1)
  
  # normalization
  const <- log(sum(exp(v.beta)) / total.prop)
  
  v.alpha <- v.alpha + const
  a.mu <- a.mu + const
  
  v.beta <- v.beta - const
  b.mu <- b.mu - const
  
  # eta
  eta.star <- abs(v.eta + rnorm(1, 0, sqrt(eta.jump)))
  eta.old <- v.eta
  
  lik.old <- f.datalik(data, v.alpha, v.beta, eta.old, v.etak, v.z, v.muk)
  lik.star <- f.datalik(data, v.alpha, v.beta, eta.star, v.etak, v.z, v.muk)
  
  prob.eta.diff <- (
    sum(lik.star, na.rm = TRUE) +
      log(dgamma(eta.star, 0.5, 0.5)) -
      sum(lik.old, na.rm = TRUE) -
      log(dgamma(eta.old, 0.5, 0.5))
  )
  
  eta.valid <- !is.infinite(prob.eta.diff) & !is.na(prob.eta.diff)
  
  eta.p <- 0
  v.eta <- eta.old
  
  if (eta.valid) {
    
    jump <- rbinom(
      1,
      1,
      exp(pmin(prob.eta.diff, 0))
    )
    
    v.eta <- ifelse(jump, eta.star, eta.old)
    
    eta.p <- exp(min(prob.eta.diff, 0))
  }
  
  # etak
  etak.new <- abs(
    v.etak + rnorm(length(v.etak), 0, sqrt(20))
  )
  
  etak.new[etak.new > 250] <- abs(
    250 - abs(etak.new[etak.new > 250] - 250)
  )
  
  etak.old <- v.etak
  
  lik.old <- f.datalik(data, v.alpha, v.beta, v.eta, etak.old, v.z, v.muk)
  lik.new <- f.datalik(data, v.alpha, v.beta, v.eta, etak.new, v.z, v.muk)
  
  prob.etak.diff <- (
    colSums(lik.new, na.rm = TRUE) +
      log(dgamma(etak.new, 5, 0.1)) -
      colSums(lik.old, na.rm = TRUE) -
      log(dgamma(etak.old, 5, 0.1))
  )
  
  etak.valid <- (
    !is.infinite(prob.etak.diff) &
      !is.na(prob.etak.diff) &
      etak.new < 1450
  )
  
  jump.new <- rep(0, length(v.etak))
  etak.p <- rep(0, length(v.etak))
  
  if (sum(etak.valid) > 0) {
    
    jump.new[etak.valid] <- rbinom(
      sum(etak.valid),
      1,
      exp(pmin(prob.etak.diff[etak.valid], 0))
    )
    
    v.etak <- jump.new * etak.new +
      (1 - jump.new) * etak.old
    
    etak.p[etak.valid] <- exp(
      pmin(prob.etak.diff[etak.valid], 0)
    )
    
  } else {
    v.etak <- etak.old
  }
  
  k.mu <- rnorm(
    1,
    mean(log(v.etak)),
    sqrt(k.s / data.j)
  )
  
  k.s <- sum((log(v.etak) - k.mu)^2) /
    rchisq(1, data.j - 1)
  
  # z
  z.star <- array(dim = dim(v.z))
  
  for (mm in 1:nrow(z.star)) {
    
    z.star[mm, ] <- sphere.update(
      u.mean = v.z[mm, ],
      u.k = z.jump[mm],
      m = ls.dim,
      rot = TRUE
    )
  }
  
  z.old <- v.z
  
  lik.old <- f.datalik(data, v.alpha, v.beta, v.eta, v.etak, z.old, v.muk)
  lik.star <- f.datalik(data, v.alpha, v.beta, v.eta, v.etak, z.star, v.muk)
  
  prob.z.diff <- (
    rowSums(lik.star, na.rm = TRUE) -
      rowSums(lik.old, na.rm = TRUE)
  )
  
  z.valid <- !is.infinite(prob.z.diff) & !is.na(prob.z.diff)
  
  jump.new <- rep(0, nrow(v.z))
  
  jump.new[z.valid] <- rbinom(
    sum(z.valid),
    1,
    exp(pmin(prob.z.diff[z.valid], 0))
  )
  
  v.z[which(jump.new == 1), ] <- z.star[which(jump.new == 1), ]
  
  z.p <- rep(0, nrow(v.z))
  
  z.p[z.valid] <- exp(
    pmin(prob.z.diff[z.valid], 0)
  )
  
  # muk
  cm <- colMeans(muk.fix)
  cm <- cm / sqrt(sum(cm^2))
  
  muk.old <- v.muk
  muk.new <- array(dim = dim(v.muk))
  
  for (m in 1:nrow(muk.new)) {
    
    muk.new[m, ] <- rmovMF(
      n = 1,
      theta = (0.1 * cm) + (muk.old[m, ] * v.etak[m]),
      alpha = 1
    )
  }
  
  muk.p <- rep(0, nrow(v.muk))
  
  v.muk <- muk.new
  
  if (is.muk.fix) {
    v.muk[muk.fix.ind, ] <- muk.fix
  }
  
  return(list(
    vjind = vjind,
    v.alpha = v.alpha,
    v.beta = v.beta,
    a.mu = a.mu,
    a.s = a.s,
    b.mu = b.mu,
    b.s = b.s,
    k.mu = k.mu,
    k.s = k.s,
    v.eta = v.eta,
    v.etak = v.etak,
    v.z = v.z,
    v.muk = v.muk,
    jump.p = c(alpha.p, beta.p, eta.p, etak.p, z.p, muk.p)
  ))
}


# -----------------------------
# 4. MCMC wrapper
# -----------------------------

f.metro <- function(
    data,
    n.iter = 3000,
    m.iter = 3,
    n.thin = 10,
    ls.dim = 3,
    total.prop,
    z.pos.init,
    muk.fix,
    muk.fix.ind,
    is.muk.fix = TRUE
) {
  
  data.n <- nrow(data)
  data.j <- ncol(data)
  
  n.keep <- n.iter / n.thin
  
  cat("Data matrix dimensions:\n")
  print(dim(data))
  
  sims <- array(
    NA,
    c(
      n.keep,
      m.iter,
      data.n + data.j + 4 + data.j + 1
    )
  )
  
  sims.latent <- array(
    NA,
    c(
      n.keep,
      m.iter,
      data.n + data.j,
      ls.dim
    )
  )
  
  sims.p <- array(
    NA,
    c(
      n.keep,
      m.iter,
      data.n + data.j + data.n + data.j + data.j + 1
    )
  )
  
  for (m in 1:m.iter) {
    
    cat("MCMC chain:", m, "\n")
    
    mc.v.alpha <- rnorm(data.n, log(50), 2)
    
    mc.v.beta <- rnorm(data.j, log(0.001), 1)
    
    mc.v.beta[exp(mc.v.beta) > 1] <- log(0.03)
    
    mc.a.mu <- log(50)
    mc.a.s <- 1
    
    mc.b.mu <- log(0.001)
    mc.b.s <- 1
    
    mc.v.z <- z.pos.init
    
    mc.v.muk <- matrix(
      rnorm(data.j * ls.dim),
      data.j,
      ls.dim
    )
    
    mc.v.muk <- mc.v.muk /
      sqrt(rowSums(mc.v.muk^2))
    
    if (is.muk.fix) {
      mc.v.muk[muk.fix.ind, ] <- muk.fix
    }
    
    mc.v.eta <- 1
    mc.v.etak <- runif(data.j, 0.001, 50)
    
    mc.k.mu <- log(40)
    mc.k.s <- log(10)
    
    ad.jump <- rep(mc.a.s, data.n)
    bd.jump <- rep(mc.b.s, data.j)
    
    eta.jump <- 4
    etak.jump <- bd.jump
    
    muk.jump <- runif(length(bd.jump), 150, 300)
    z.jump <- runif(data.n, 1, 10)
    
    mc.j.ind <- rep(0, length(muk.jump))
    
    last.50p <- array(
      NA,
      c(
        (length(sims[1,1,]) +
           length(sims.latent[1,1,,1]) - 4),
        50
      )
    )
    
    p.ct <- 0
    
    for (t in 1:n.iter) {
      
      temp <- f.update(
        data = data,
        v.alpha = mc.v.alpha,
        v.beta = mc.v.beta,
        v.eta = mc.v.eta,
        v.etak = mc.v.etak,
        v.z = mc.v.z,
        v.muk = mc.v.muk,
        a.mu = mc.a.mu,
        a.s = mc.a.s,
        b.mu = mc.b.mu,
        b.s = mc.b.s,
        ad.jump = ad.jump,
        bd.jump = bd.jump,
        eta.jump = eta.jump,
        etak.jump = etak.jump,
        muk.jump = muk.jump,
        z.jump = z.jump,
        k.mu = mc.k.mu,
        k.s = mc.k.s,
        vjind = mc.j.ind,
        total.prop = total.prop,
        muk.fix = muk.fix,
        muk.fix.ind = muk.fix.ind,
        ls.dim = ls.dim,
        is.muk.fix = is.muk.fix
      )
      
      mc.v.alpha <- temp$v.alpha
      mc.v.beta <- temp$v.beta
      
      mc.a.mu <- temp$a.mu
      mc.a.s <- temp$a.s
      
      mc.b.mu <- temp$b.mu
      mc.b.s <- temp$b.s
      
      mc.v.z <- temp$v.z
      mc.v.muk <- temp$v.muk
      
      mc.v.eta <- temp$v.eta
      mc.v.etak <- temp$v.etak
      
      mc.k.mu <- temp$k.mu
      mc.k.s <- temp$k.s
      
      mc.j.ind <- temp$vjind
      
      p.ct <- p.ct + 1
      last.50p[, p.ct] <- temp$jump.p
      
      if (p.ct == 50) {
        
        jump.d <- c(
          ad.jump,
          bd.jump,
          eta.jump,
          etak.jump,
          z.jump,
          muk.jump
        )
        
        p.ct <- 0
        
        p.mean <- rowMeans(last.50p, na.rm = TRUE)
        
        p.mean[p.mean <= 0 | is.na(p.mean)] <- 0.001
        
        jump.d.keep <- jump.d
        
        jump.d <- pmax(
          pmin(log(0.4) * jump.d / log(p.mean), 10 * jump.d),
          0.0001
        )
        
        ad.jump <- jump.d[1:data.n]
        
        bd.jump <- jump.d[(data.n + 1):(data.n + data.j)]
        
        eta.jump <- jump.d[(1 + data.n + data.j)]
        
        etak.jump <- jump.d[
          (1 + data.n + data.j + 1):
            (1 + data.n + data.j + data.j)
        ]
        
        jump.d.keep <- pmax(
          pmin(
            1 / ((log(0.4) / log(p.mean)) * (1 / jump.d.keep)),
            750
          ),
          0.0001
        )
        
        z.jump <- jump.d.keep[
          (1 + data.n + data.j + data.j + 1):
            (1 + data.n + data.j + data.j + data.n)
        ]
        
        muk.jump <- jump.d.keep[
          (1 + data.n + data.j + data.j + data.n + 1):
            (1 + data.n + data.j + data.j + data.n + data.j)
        ]
      }
      
      if (t %% n.thin == 0) {
        
        if (t %% 100 == 0) {
          cat(
            "Iteration:",
            t,
            "| Chain:",
            m,
            "|",
            date(),
            "\n"
          )
        }
        
        sims[t / n.thin, m, ] <- c(
          mc.v.alpha,
          mc.v.beta,
          mc.a.mu,
          mc.a.s,
          mc.b.mu,
          mc.b.s,
          mc.v.etak,
          mc.v.eta
        )
        
        sims.latent[t / n.thin, m, , ] <- rbind(
          mc.v.z,
          mc.v.muk
        )
        
        sims.p[t / n.thin, m, ] <- temp$jump.p
      }
    }
  }
  
  return(list(
    sims = sims,
    sims.p = sims.p,
    sims.latent = sims.latent,
    data = data,
    data.n = data.n,
    data.j = data.j,
    n.iter = n.iter,
    m.iter = m.iter,
    n.thin = n.thin
  ))
}


# -----------------------------
# 5. Posterior extraction
# -----------------------------

getPosterior <- function(out, n.iter, m.iter, n.thin, n) {
  
  est.degrees <- NULL
  est.eta <- NULL
  est.latent.pos <- NULL
  
  start_ind <- (n.iter / n.thin / 2) + 1
  end_ind <- n.iter / n.thin
  
  for (n.ind in start_ind:end_ind) {
    
    for (m.ind in 1:m.iter) {
      
      est.degrees <- rbind(
        est.degrees,
        out$sims[n.ind, m.ind, 1:n]
      )
      
      est.eta <- c(
        est.eta,
        out$sims[n.ind, m.ind, length(out$sims[n.ind, m.ind, ])]
      )
      
      est.latent.pos <- rbind(
        est.latent.pos,
        c(out$sims.latent[n.ind, m.ind, 1:n, ])
      )
    }
  }
  
  return(list(
    est.degrees = est.degrees,
    est.eta = est.eta,
    est.latent.pos = est.latent.pos
  ))
}


getGi <- function(est.degrees, est.eta) {
  
  gi.m <- NULL
  
  for (ind in 1:length(est.eta)) {
    
    nexp.gi <- sqrt(
      sum(exp(est.degrees[ind, ])) *
        cp.fcn(est.eta[ind]) /
        cp.fcn(0.000001)
    )
    
    gi <- (
      exp(est.degrees[ind, ]) /
        nexp.gi *
        cp.fcn(est.eta[ind]) /
        cp.fcn(0.000001)
    )
    
    gi.m <- rbind(gi.m, gi)
  }
  
  return(log(gi.m))
}


# -----------------------------
# 6. Sparse graph simulation
# -----------------------------

simulate.graph.once <- function(
    z,
    g,
    eta,
    d.ARD,
    z.ARD,
    g.ARD
) {
  
  n.ARD <- length(g.ARD)
  
  adjexp <- matrix(
    0,
    nrow = n.ARD,
    ncol = n.ARD
  )
  
  for (i in 1:(n.ARD - 1)) {
    
    for (j in (i + 1):n.ARD) {
      
      val <- exp(
        g.ARD[i] +
          g.ARD[j] +
          eta * sum(z.ARD[i, ] * z.ARD[j, ])
      )
      
      adjexp[i, j] <- val
      adjexp[j, i] <- val
    }
  }
  
  const <- sum(exp(d.ARD)) / sum(adjexp)
  
  n <- length(g)
  
  ii <- integer(0)
  jj <- integer(0)
  
  for (i in 1:(n - 1)) {
    
    for (j in (i + 1):n) {
      
      p.ij <- exp(
        g[i] +
          g[j] +
          eta * sum(z[i, ] * z[j, ])
      ) * const
      
      if (runif(1) < min(p.ij, 1)) {
        
        ii <- c(ii, i)
        jj <- c(jj, j)
      }
    }
  }
  
  edges <- cbind(from = ii, to = jj)
  
  return(edges)
}


simulate.graph.all.sparse <- function(
    est.degrees.ARD,
    est.eta,
    est.latent.pos.ARD,
    est.gi.ARD,
    est.gi,
    est.latent.pos,
    ls.dim,
    out_dir = "village51_synthetic_graphs",
    prefix = "village51",
    compress = TRUE,
    progress_every = 10
) {
  
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  n.ARD <- dim(est.degrees.ARD)[2]
  n <- dim(est.gi)[2]
  n_samp <- length(est.eta)
  
  for (ind in seq_len(n_samp)) {
    
    z <- matrix(
      est.latent.pos[ind, ],
      byrow = FALSE,
      nrow = n,
      ncol = ls.dim
    )
    
    z.ARD <- matrix(
      est.latent.pos.ARD[ind, ],
      byrow = FALSE,
      nrow = n.ARD,
      ncol = ls.dim
    )
    
    edges <- simulate.graph.once(
      z = z,
      g = est.gi[ind, ],
      eta = est.eta[ind],
      d.ARD = est.degrees.ARD[ind, ],
      z.ARD = z.ARD,
      g.ARD = est.gi.ARD[ind, ]
    )
    
    file_name <- sprintf("%s_edges_%06d.csv", prefix, ind)
    
    file_out <- file.path(out_dir, file_name)
    
    if (compress) {
      
      con <- gzfile(
        paste0(file_out, ".gz"),
        open = "wt"
      )
      
      write.csv(edges, con, row.names = FALSE)
      
      close(con)
      
    } else {
      
      write.csv(edges, file_out)
    }
    
    if (
      progress_every > 0 &&
      (ind %% progress_every == 0 || ind == n_samp)
    ) {
      
      cat(
        sprintf(
          "[%s] Saved %d/%d graphs\n",
          format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
          ind,
          n_samp
        )
      )
    }
    
    rm(edges, z, z.ARD)
    gc()
  }
  
  meta <- list(
    out_dir = out_dir,
    n_graphs = n_samp,
    prefix = prefix,
    compressed = compress
  )
  
  saveRDS(
    meta,
    file = file.path(out_dir, "meta.rds")
  )
  
  return(meta)
}


# -----------------------------
# 7. Main pipeline
# -----------------------------

run_synthetic_generation <- function(
    ard_data,
    total.prop,
    n.iter = 1000,
    m.iter = 2,
    n.thin = 20,
    ls.dim = 3,
    out_dir = "village51_synthetic_graphs",
    prefix = "village51",
    seed = 123
) {
  
  set.seed(seed)
  
  n <- nrow(ard_data)
  K <- ncol(ard_data)
  
  if (K <= 4) {
    stop("The model requires more than 4 groups/circles.")
  }
  
  z.pos.init <- generateRandomInitial(n, ls.dim)
  
  muk.fix.ind <- sample(
    1:K,
    size = 4,
    replace = FALSE
  )
  
  muk.fix <- matrix(
    rnorm(4 * ls.dim),
    nrow = 4,
    ncol = ls.dim
  )
  
  muk.fix <- sweep(
    muk.fix,
    MARGIN = 1,
    1 / sqrt(rowSums(muk.fix^2)),
    `*`
  )
  
  out <- f.metro(
    data = ard_data,
    total.prop = total.prop,
    n.iter = n.iter,
    m.iter = m.iter,
    n.thin = n.thin,
    z.pos.init = z.pos.init,
    muk.fix = muk.fix,
    muk.fix.ind = muk.fix.ind,
    ls.dim = ls.dim,
    is.muk.fix = TRUE
  )
  
  cat("MCMC finished.\n")
  
  posterior <- getPosterior(
    out,
    n.iter,
    m.iter,
    n.thin,
    n
  )
  
  est.degrees <- posterior$est.degrees
  est.eta <- posterior$est.eta
  est.latent.pos <- posterior$est.latent.pos
  
  est.gi <- getGi(est.degrees, est.eta)
  
  cat("Starting sparse graph simulations.\n")
  
  meta <- simulate.graph.all.sparse(
    est.degrees.ARD = est.degrees,
    est.eta = est.eta,
    est.latent.pos.ARD = est.latent.pos,
    est.gi.ARD = est.gi,
    est.gi = est.gi,
    est.latent.pos = est.latent.pos,
    ls.dim = ls.dim,
    out_dir = out_dir,
    prefix = prefix,
    compress = TRUE
  )
  
  saveRDS(
    out,
    file = file.path(out_dir, "mcmc_output.rds")
  )
  
  saveRDS(
    posterior,
    file = file.path(out_dir, "posterior_samples.rds")
  )
  
  cat("Finished.\n")
  
  return(meta)
}


# -----------------------------
# 8. Run
# -----------------------------

result <- run_synthetic_generation(
  ard_data = ard_data,
  total.prop = total.prop,
  n.iter = 3000,
  m.iter = 3,
  n.thin = 10,
  ls.dim = 3,
  out_dir = "village51_synthetic_graphs",
  prefix = "village51",
  seed = 123
)

print(result)

