## =============================================================================
##  applications.R -- applications to the linear fSDE (Theorem 3, Corollary 3)
##  Output (in R): numerical checks, Tables 3-4 and Figures 2-3
##  Monte Carlo (Figure 3): yuima, M = 1e5 paths
## =============================================================================

source("packages.R")
source("main_functions.R")

## =============================================================================
## Theorem 3(b): dynamics check
## =============================================================================

mA <- list(G = function(t) t + 1, alpha = -0.5, beta = -1, sigma = 1, x0 = 2)
dG <- function(t) 1 + 0 * t
H <- 0.7; p <- 1.5; t0 <- 1.4; eps <- 1e-5
Pf  <- function(t, q) mom_fsde(q, t, mA, H, side = 1)
lhs <- (Pf(t0 + eps, p) - Pf(t0 - eps, p)) / (2 * eps)
rhs <- p * rho_t(t0, mA, H, dG) * Pf(t0, p) + p * kappa_t(t0, mA, H) * Pf(t0, p - 1)
cat("Eq. (dyn):  finite difference =", lhs, "  right-hand side =", rhs, "\n")

kfd <- sigma_X(t0, mA, H) * (theta_t(t0 + eps, mA, H) - theta_t(t0 - eps, mA, H)) / (2 * eps)
cat("kappa: finite difference =", kfd, "  closed form =", kappa_t(t0, mA, H), "\n")

# Corollary 1(c): tail bounds
for (pv in c(0.5, 1.5, 3)) for (thv in c(3, 6)) {
  U <- gammaf(pv + 1) * dnorm(thv) * thv^(-pv - 1)
  cat(sprintf("p = %.1f, theta = %.0f: lower = %.4e  value = %.4e  upper = %.4e\n",
              pv, thv, U * (1 - (pv + 1) * (pv + 2) / (2 * thv^2)), pi_p(pv, -thv), U))
}

## =============================================================================
## Figure 2
## =============================================================================

mB <- list(G = function(t) 1 + 0 * t, alpha = 0, beta = 1, sigma = 1, x0 = 0.5)
tg <- seq(0.05, 30, length.out = 3000); Hs <- c(0.5, 0.6, 0.7, 0.8, 0.9)
fig2 <- do.call(rbind, lapply(Hs, function(H) {
  th <- theta_t(tg, mB, H)
  rbind(data.frame(H = H, Type = "exact", kind = "line", x = tg,
                   value = sapply(tg, mom_fsde, p = 1.5, m = mB, H = H, side = -1)),
        data.frame(H = H, Type = "tail", kind = "line", x = tg,
                   value = sigma_X(tg, mB, H)^1.5 * gammaf(2.5) * dnorm(th) * th^-2.5))
}))
fig2 <- fig2[fig2$value > 1e-12, ]
df2 <- fig2 %>%
  mutate(source = ifelse(Type == "exact", "Exact", "Tail approximation"),
         time   = x,
         H_fac  = factor(H))

H_unique  <- sort(unique(df2$H))
labels_H  <- sapply(H_unique, function(v) paste0("H = ", v))
Color2 <- ggsci::pal_lancet("lanonc", alpha = 0.7)(length(H_unique))
names(Color2) <- as.character(H_unique)

ltype_map2 <- c("Exact" = "solid", "Tail approximation" = "dashed")

g2 <- ggplot(df2, aes(x = time, y = value, group = interaction(H_fac, source),
                      color = H_fac, linetype = source)) +
  geom_line(linewidth = 0.9) +
  scale_linetype_manual(name = NULL, values = ltype_map2) +
  scale_color_manual(name = NULL, values = Color2, labels = labels_H) +
  scale_x_continuous(breaks = pretty_breaks(n = 12)) +
  scale_y_log10(breaks = scales::breaks_log(n = 10),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(x = "Time", y = expression(bold(E*"["*(X[t]^"-")^{3/2}*"]"))) +
  theme_fig(legend_size = 14) +
  theme(legend.key.width = unit(1.8, "cm"))+
  guides(
    linetype = guide_legend(order = 1, override.aes = list(linewidth = 1, color = "black")),
    color    = guide_legend(order = 2, nrow = 1, keywidth = unit(0.2, "cm"),
                            override.aes = list(linewidth = 4, linetype = "solid"))
  )
  
g2

## ggsave("figures/fig2_positivity.pdf", g2, width = 11.69, height = 6.5)

## =============================================================================
## Nile data
## =============================================================================

y <- as.numeric(Nile); years <- 1871:1970

## =============================================================================
## Table 3
## =============================================================================

fit_all  <- fit_fgn(y)
fit_post <- fit_fgn(y[years >= 1899])
tab3 <- data.frame(period = c("1871-1970", "1899-1970"),
                   n = c(fit_all$n, fit_post$n), H = c(fit_all$H, fit_post$H),
                   se_H = c(fit_all$seH, fit_post$seH), mu = c(fit_all$b, fit_post$b),
                   gamma = c(fit_all$s, fit_post$s))
cat("\nTable 3 (exact MLE):\n"); print(tab3, digits = 4)

## =============================================================================
## Table 4
## =============================================================================

nile_model <- function(f, d, H) list(G = function(t) 1 + 0 * t, G_expr = "1", dG_expr = "0",
                                     alpha = 0, beta = f$b - d, sigma = f$s, x0 = 0)
fit_bm <- list(b = mean(y), s = sd(y))
p_def  <- 1.5
hor    <- c(5, 10, 20, 30, 40)
b0 <- fit_all$b - mean(y)
cat(sprintf("\nd = ybar: beta = %.2f, beta/gamma = %.4f, theta_t at t = 5, 10, 20, 30, 40:", b0, b0 / fit_all$s),
    round(b0 * hor^(1 - fit_all$H) / fit_all$s, 4), "\n")
tab4 <- NULL
for (q in c(1, 0.9)) {
  d  <- q * mean(y)
  mF <- nile_model(fit_all, d); mB <- nile_model(fit_bm, d)
  for (t in hor) {
    tab4 <- rbind(tab4, data.frame(
      q = q, t = t,
      P_fsde = mom_fsde(0, t, mF, fit_all$H, side = -1),
      P_bm   = mom_fsde(0, t, mB, 0.5, side = -1),
      P_emp  = empirical_deficit(y, d, t, 0),
      D_fsde = mom_fsde(p_def, t, mF, fit_all$H, side = -1)^(1 / p_def),
      D_bm   = mom_fsde(p_def, t, mB, 0.5, side = -1)^(1 / p_def),
      D_emp  = empirical_deficit(y, d, t, p_def)^(1 / p_def)))
  }
}
cat("\nTable 4 (deficits):\n"); print(tab4, digits = 4)

## =============================================================================
## Figure 3
## =============================================================================

d  <- mean(y); tg <- 1:40
mF <- nile_model(fit_all, d); mB <- nile_model(fit_bm, d)
DF <- sapply(tg, function(t) mom_fsde(p_def, t, mF, fit_all$H, side = -1)^(1 / p_def))
DB <- sapply(tg, function(t) mom_fsde(p_def, t, mB, 0.5, side = -1)^(1 / p_def))
DE <- sapply(tg, function(t) empirical_deficit(y, d, t, p_def)^(1 / p_def))
set.seed(2026)
M <- 1e5; dt <- 1; tsim <- c(1, 2, 3, 5, 7, 10, 15, 20, 30, 40)
n_cores <- detectCores()
XF <- yuima_paths(mF, fit_all$H, T = 40, n = round(40 / dt), M = M, keep = round(tsim / dt) + 1, cores = n_cores)
XB <- yuima_paths(mB, 0.5,       T = 40, n = round(40 / dt), M = M, keep = round(tsim / dt) + 1, cores = n_cores)
SF <- rowMeans(pmax(-XF, 0)^p_def)^(1 / p_def)
SB <- rowMeans(pmax(-XB, 0)^p_def)^(1 / p_def)
lab_F <- "fractional"; lab_B <- "standard"
df3 <- rbind(data.frame(Model = lab_F, source = "Theorem 3",   time = tg,   value = DF),
             data.frame(Model = lab_B, source = "Theorem 3",   time = tg,   value = DB),
             data.frame(Model = lab_F, source = "Monte Carlo", time = tsim, value = SF),
             data.frame(Model = lab_B, source = "Monte Carlo", time = tsim, value = SB))
df3$Model <- factor(df3$Model, levels = c(lab_F, lab_B))
df3e <- data.frame(time = tg, value = DE, data = "Empirical")
Color3 <- setNames(ggsci::pal_lancet("lanonc", alpha = 0.6)(2), c(lab_F, lab_B))
labels_models <- c(bquote(bold("fractional, " * hat(H) == .(sprintf("%.3f", fit_all$H)))),
                   bquote(bold("standard, " * H == "1/2" * phantom(hat(H)))))
ltype_map3 <- c("Theorem 3" = "dashed", "Monte Carlo" = "solid")

g3a <- ggplot(data.frame(year = years, value = cumsum(y - mean(y)),
                         panel = "Cumulative departures"),
              aes(x = year, y = value)) +
  geom_hline(yintercept = 0, colour = "grey50", linetype = "dotted") +
  geom_vline(xintercept = 1898, colour = "#AD002AFF", linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = 1898, y = 300, label = "1898", hjust = -0.15,
           fontface = "bold.italic", size = 3.2) +
  geom_line(colour = "#00468BFF", linewidth = 1) +
  scale_x_continuous(breaks = pretty_breaks(n = 6)) +
  scale_y_continuous(breaks = pretty_breaks(n = 10)) +
  labs(x = "Year", y = "Cumulative departure from the mean") +
  facet_wrap(~ panel) +
  theme_fig()

g3b <- ggplot() +
  geom_line(
    data = df3 %>% filter(source == "Monte Carlo"),
    aes(x = time, y = value, group = Model, color = Model, linetype = source),
    linewidth = 1
  ) +
  geom_line(
    data = df3 %>% filter(source == "Theorem 3"),
    aes(x = time, y = value, group = Model, linetype = source),
    color = "black", linewidth = 0.8
  ) +
  geom_point(data = df3e, aes(x = time, y = value, shape = data),
             color = "#42B540FF", size = 3) +
  scale_linetype_manual(name = NULL, values = ltype_map3) +
  scale_color_manual(name = NULL, values = Color3, labels = labels_models) +
  scale_shape_manual(name = NULL, values = c(Empirical = 16)) +
  scale_x_log10(breaks = c(1, 2, 3, 5, 10, 20, 30, 40)) +
  scale_y_log10(breaks = scales::breaks_log(n = 12), labels = scales::label_number(big.mark = "")) +
  labs(x = "Horizon t (years)", y = expression(bold((E*"["*(X[t]^"-")^{3/2}*"]")^{2/3}))) +
  facet_wrap(~ "Deficit moment, d = mean flow") +
  theme_fig(legend_size = 14) +
  theme(legend.key.width            = unit(1.8, "cm"),
        legend.box                  = "vertical",
        legend.margin               = margin(0, 0, 0, 0),
        legend.box.spacing   = unit(0.2, "cm"),
        legend.spacing.y            = unit(0.1, "cm"),
        axis.title.y                = element_text(face = "bold", margin = margin(r = -4)),
        axis.text.y                 = element_text(margin = margin(r = 1)),
        legend.position.inside      = c(0.02, 0.98),
        legend.justification.inside = c(0, 1),
        legend.background           = element_blank()) +
  guides(
    linetype = guide_legend(order = 1, position = "inside",
                            override.aes = list(linewidth = 1)),
    shape    = guide_legend(order = 2, position = "inside"),
    color    = guide_legend(order = 3, nrow = 1, keywidth = unit(0.2, "cm"),
                            theme = theme(legend.text = element_text(size = 14),
                                          legend.key.spacing.x = unit(0.8, "cm")),
                            override.aes = list(linewidth = 4, linetype = "solid"))
  )

g3 <- g3a | g3b
g3

## ggsave("figures/fig3_nile.pdf", g3, width = 11.69, height = 6)