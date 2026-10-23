## =============================================================================
##  main_functions.R -- special functions, moments, simulation and plotting tools
##  (sourced by simulation.R and applications.R)
## =============================================================================

kummer <- function(a, b, z) hyperg_1F1(a, b, z)

gammaf <- function(x) exp(lngamma(x))
rgi    <- function(x) gammainv(x)

besselKf <- function(x, nu) bessel_Knu(nu, x)

hermite_imag <- function(p, y) {
  z <- -y^2
  2^p * sqrt(pi) * (kummer(-p / 2, 1 / 2, z) * rgi((1 - p) / 2) -
                    2i * y * kummer((1 - p) / 2, 3 / 2, z) * rgi(-p / 2))
}

M_herm <- function(p, mu, s) (-1i * s / sqrt(2))^p * hermite_imag(p, mu / (sqrt(2) * s))

pi_p <- function(p, theta) {
  f  <- function(x) exp(p * log(x) - (x - theta)^2 / 2 - 0.5 * log(2 * pi))
  hi <- max(theta, 0) + 12 + 3 * sqrt(max(p, 0))
  integrate(f, 0, hi, rel.tol = 1e-12, subdivisions = 2000L)$value
}

P_p <- function(p, mu, s) s^p * pi_p(p, mu / s)
N_p <- function(p, mu, s) s^p * pi_p(p, -mu / s)

parts_from_hermite <- function(p, mu, s) {
  M <- M_herm(p, mu, s)
  N <- -Im(M) / sin(p * pi)
  P <- Re(M) - cos(p * pi) * N
  c(P = P, N = N, A = P + N)
}

abs_kummer <- function(p, mu, s)
  s^p * 2^(p / 2) * gammaf((p + 1) / 2) / sqrt(pi) * kummer(-p / 2, 1 / 2, -mu^2 / (2 * s^2))

ladder_up <- function(p0, nmax, mu, s, side = c("pos", "neg")) {
  side <- match.arg(side); m <- if (side == "pos") mu else -mu
  out <- numeric(nmax + 2)
  out[1] <- s^(p0 - 1) * pi_p(p0 - 1, m / s)
  out[2] <- s^p0 * pi_p(p0, m / s)
  for (n in 1:nmax) {
    q <- p0 + n - 1
    out[n + 2] <- m * out[n + 1] + q * s^2 * out[n]
  }
  setNames(out, p0 + (-1):nmax)
}

mu_X    <- function(t, m) m$G(t)^m$alpha * (m$G(0)^(-m$alpha) * m$x0 + m$beta * t)
sigma_X <- function(t, m, H) m$sigma * m$G(t)^m$alpha * t^H
theta_t <- function(t, m, H) (m$G(0)^(-m$alpha) * m$x0 + m$beta * t) / (m$sigma * t^H)

mom_fsde <- function(p, t, m, H, side = 1)
  sigma_X(t, m, H)^p * pi_p(p, side * theta_t(t, m, H))

rho_t   <- function(t, m, H, dG) m$alpha * dG(t) / m$G(t) + H / t
kappa_t <- function(t, m, H) m$beta * m$G(t)^m$alpha - H * mu_X(t, m) / t

sde_coef <- function(m) list(
  drift     = sprintf("(%s) * (%s) / (%s) * x + (%s) * (%s)^(%s)",
                      m$alpha, m$dG_expr, m$G_expr, m$beta, m$G_expr, m$alpha),
  diffusion = sprintf("(%s) * (%s)^(%s)", m$sigma, m$G_expr, m$alpha))

yuima_paths <- function(m, H, T, n, M, keep = seq_len(n + 1), cores = n_cores) {
  cf    <- sde_coef(m)
  sizes <- diff(round(seq(0, M, length.out = cores + 1)))
  cl <- makeCluster(cores)
  on.exit(stopCluster(cl))
  registerDoParallel(cl)
  foreach(Mi = sizes, .combine = cbind, .packages = "yuima") %dorng% {
    mod  <- setModel(drift = cf$drift, diffusion = cf$diffusion, hurst = H,
                     solve.variable = "x", state.variable = "x", time.variable = "t")
    samp <- setSampling(Terminal = T, n = n)
    suppressWarnings(vapply(seq_len(Mi), function(j) {
      y <- simulate(mod, sampling = samp, xinit = m$x0)
      as.numeric(get.zoo.data(y)[[1]])[keep]
    }, numeric(length(keep))))
  }
}

fgn_acf <- function(n, H) {
  k <- 0:(n - 1)
  0.5 * (abs(k + 1)^(2 * H) - 2 * k^(2 * H) + abs(k - 1)^(2 * H))
}

fgn_profile <- function(H, y) {
  n <- length(y); L <- chol(toeplitz(fgn_acf(n, H)))
  one <- backsolve(L, rep(1, n), transpose = TRUE)
  yy  <- backsolve(L, y, transpose = TRUE)
  b   <- sum(one * yy) / sum(one^2)
  s2  <- sum((yy - b * one)^2) / n
  list(nll = 0.5 * (n * log(s2) + 2 * sum(log(diag(L))) + n), b = b, s = sqrt(s2))
}

fit_fgn <- function(y) {
  o  <- optimize(function(H) fgn_profile(H, y)$nll, c(0.01, 0.99))
  H  <- o$minimum; pr <- fgn_profile(H, y); h <- 1e-3
  d2 <- (fgn_profile(H + h, y)$nll - 2 * o$objective + fgn_profile(H - h, y)$nll) / h^2
  list(H = H, b = pr$b, s = pr$s, seH = 1 / sqrt(d2), n = length(y))
}

empirical_deficit <- function(y, d, t, p, x0 = 0) {
  cs <- c(0, cumsum(y - d)); n <- length(y)
  incr <- cs[(t + 1):(n + 1)] - cs[1:(n - t + 1)]
  if (p == 0) return(mean(x0 + incr < 0))
  mean(pmax(-(x0 + incr), 0)^p)
}

theme_fig <- function(legend_size = 10, strip_size = 11) {
  theme_minimal() +
    theme(
      panel.grid       = element_blank(),
      panel.grid.major = element_line(color = "grey80", linewidth = 0.5),
      panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.75),
      panel.spacing    = unit(0, "lines"),
      legend.position      = "top",
      legend.justification = "centre",
      legend.box           = "horizontal",
      legend.spacing.x     = unit(1, "cm"),
      legend.key.width     = unit(1.2, "cm"),
      legend.text          = element_text(size = legend_size, face = "bold.italic"),
      legend.margin        = margin(t = 0, r = 0, b = -10, l = 0),
      legend.box.spacing   = unit(0.5, "cm"),
      axis.title.x     = element_text(face = "bold", size = 13),
      axis.title.y     = element_text(face = "bold", size = 13),
      axis.text        = element_text(),
      plot.margin      = unit(c(0, 0.05, 0, 0), "cm"),
      strip.text       = element_text(size = strip_size, face = "bold.italic",
                                      margin = margin(2, 1, 2, 1)),
      strip.background = element_rect(fill = "#f0f0f0", color = "black", linewidth = 0.75)
    )
}

frac_label <- function(name, v) paste0(name, " = ", as.character(MASS::fractions(v)))
