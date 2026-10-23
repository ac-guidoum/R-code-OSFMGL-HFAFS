## =============================================================================
##  simulation.R -- numerical verification of Theorems 1-2 and Monte Carlo study
##  Output (in R): Tables 1-2, numerical checks and Figure 1
##  Monte Carlo: yuima, M = 1e5 paths (Figure 1)
## =============================================================================

source("packages.R")
source("main_functions.R")
set.seed(2026)

## =============================================================================
## Table 1
## =============================================================================

cfgs <- list(c(0.5, 1, 1), c(1.3, 1, 1), c(2.2, -0.5, 1), c(3.7, 0.8, 0.6),
             c(0.25, -2, 1.5), c(2.7, -1.2, 0.8))
tab1 <- t(sapply(cfgs, function(cf) {
  p <- cf[1]; mu <- cf[2]; s <- cf[3]; h <- parts_from_hermite(p, mu, s)
  c(p = p, mu = mu, s = s,
    P_hermite = h[["P"]], P_quad = P_p(p, mu, s),
    N_hermite = h[["N"]], N_quad = N_p(p, mu, s),
    A_hermite = h[["A"]], A_kummer = abs_kummer(p, mu, s))
}))
cat("\nTable 1 (Theorem 1 / Corollary 1):\n"); print(round(tab1, 5))

## =============================================================================
## Table 2
## =============================================================================

lad <- ladder_up(0.3, 8, mu = -0.7, s = 1.2)
dir <- sapply(as.numeric(names(lad)), function(q) P_p(q, -0.7, 1.2))
tab2 <- cbind(order = as.numeric(names(lad)), ladder = lad, direct = dir, rel_err = abs(lad / dir - 1))
cat("\nTable 2 (fractional ladder):\n"); print(tab2)

## =============================================================================
## Numerical checks quoted in the text
## =============================================================================

pp <- 1.7; mu <- 0.4; s <- 0.9
cat("\nComplex ladder (Eq. ladderM):", format(M_herm(pp + 1, mu, s)), " vs ",
    format(mu * M_herm(pp, mu, s) + pp * s^2 * M_herm(pp - 1, mu, s)), "\n")

th <- 1.3
cat("Bessel seed: pi_{-1/2}(-theta) =", pi_p(-0.5, -th), " Bessel form =",
    gammaf(0.5) / sqrt(2 * pi) * exp(-th^2 / 4) * sqrt(th / (2 * pi)) * besselKf(th^2 / 4, 0.25), "\n")

for (p in c(-0.6, 0.4, 2.3)) for (thv in c(-1.2, 0.7)) {
  e <- 1e-3
  d1 <- (pi_p(p, thv + e) - pi_p(p, thv - e)) / (2 * e)
  d2 <- (pi_p(p, thv + e) - 2 * pi_p(p, thv) + pi_p(p, thv - e)) / e^2
  cat(sprintf("Weber residual p = %.1f, theta = %.1f: %.2e\n", p, thv, d2 + thv * d1 - p * pi_p(p, thv)))
}

## =============================================================================
## Figure 1
## =============================================================================

mA <- list(G = function(t) t + 1, G_expr = "t + 1", dG_expr = "1",
           alpha = -0.5, beta = -1, sigma = 1, x0 = 2)
T <- 4; dt <- 1e-2; n <- round(T / dt); M <- 1e5; tt <- seq(0, T, length.out = n + 1)
n_cores <- detectCores() 
ps <- c(-0.25, 0.5, 1.5, 2.5)
tf  <- seq(0.02, T, length.out = 200)
idx <- round(seq(0.1, T, by = 0.2) / dt) + 1
fig1 <- NULL
for (H in c(0.6, 0.7, 0.8, 0.9)) {
  X <- yuima_paths(mA, H, T, n, M, keep = idx, cores = n_cores)
  for (p in ps) {
    e <- if (p > 0) 1 / p else 1
    fig1 <- rbind(fig1,
      data.frame(H = H, p = p, Part = "pos", kind = "line", x = tf,
                 value = sapply(tf, mom_fsde, p = p, m = mA, H = H, side = 1)^e),
      data.frame(H = H, p = p, Part = "neg", kind = "line", x = tf,
                 value = sapply(tf, mom_fsde, p = p, m = mA, H = H, side = -1)^e),
      data.frame(H = H, p = p, Part = "pos", kind = "point", x = tt[idx],
                 value = rowMeans(ifelse(X > 0, abs(X)^p, 0))^e),
      data.frame(H = H, p = p, Part = "neg", kind = "point", x = tt[idx],
                 value = rowMeans(ifelse(X < 0, abs(X)^p, 0))^e))
  }
}

df_all <- fig1 %>%
  mutate(source = ifelse(kind == "line", "Theorem 3", "Monte Carlo"),
         time   = x,
         p_fac  = factor(p),
         type   = Part,
         H_lab  = paste0("bolditalic('H = ", H, "')"))

p_unique    <- sort(unique(df_all$p))
labels_expr <- sapply(p_unique, function(v) frac_label("p", v))
Color <- ggsci::pal_lancet("lanonc", alpha = 0.5)(length(p_unique))
names(Color) <- as.character(p_unique)

ltype_map <- c("Theorem 3" = "dashed", "Monte Carlo" = "solid")

g1 <- ggplot() +
  geom_line(
    data = df_all %>% filter(source == "Monte Carlo"),
    aes(x = time, y = value, group = interaction(p_fac, type),
        color = factor(p), linetype = source),
    linewidth = 1
  ) +
  geom_line(
    data = df_all %>% filter(source == "Theorem 3"),
    aes(x = time, y = value, group = interaction(p_fac, type), linetype = source),
    color = "black", linewidth = 0.8
  ) +
  scale_linetype_manual(name = NULL, values = ltype_map) +
  scale_color_manual(name = NULL, values = Color, labels = labels_expr) +
  scale_x_continuous(breaks = pretty_breaks(n = 8), labels = number_format(accuracy = 0.1)) +
  scale_y_continuous(breaks = pretty_breaks(n = 10), labels = number_format(accuracy = 0.1)) +
  labs(x = "Time", y = bquote("One-sided fractional moments")) +
  facet_wrap(~ H_lab, labeller = label_parsed, nrow = 2, ncol = 2) +
  theme_fig(legend_size = 14) +
  theme(legend.key.width = unit(1.8, "cm"))+
  guides(
    linetype = guide_legend(order = 1, override.aes = list(linewidth = 1)),
    color    = guide_legend(order = 2, nrow = 1,
                            theme = theme(legend.justification = "right"),
                            keywidth = unit(0.2, "cm"),
                            override.aes = list(linewidth = 4, linetype = "solid"))
  )
  
g1
