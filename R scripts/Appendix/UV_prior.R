# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : UV_prior.R                                                 ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.es                                  ║
# ║ Date          : 28-06-2025                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Description
# ─────────────────────────────────────────────────────────────────────────────
# Prior predictive checks for unit-vector factor loadings
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Load Packages
# ─────────────────────────────────────────────────────────────────────────────

library(ggplot2)
library(bayesplot)
library(ggforce)
library(patchwork)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: User-Defined Functions
# ─────────────────────────────────────────────────────────────────────────────

# Function to simulate standardized factor loadings based on communalities
sim_lambda_UV <- function(n_rows, beta_pars) {
  # Simulate communalities
  h2 <- rbeta(n_rows, beta_pars[1], beta_pars[2])
  # Independent gaussian values
  Z <- matrix(rnorm(n_rows * 2), ncol = 2)
  # Compute Standardized factor loadings
  Lambda <- matrix(NA, nrow = n_rows, ncol = 2)
  for(i in 1:n_rows) Lambda[i,] <- Z[i,]/sqrt(sum(Z[i,]^2)) * sqrt(h2[i])
  # Return final values
  return(data.frame(
    Lambda_1 = Lambda[,1],
    Lambda_2 = Lambda[,2]
  ))
}

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: Prior Predictive distribution: communalities and factor loadings
# ─────────────────────────────────────────────────────────────────────────────

# beta distribution parameters for communalities
beta_pars_list <- list(
  c1 = c(1, 1),
  c2 = c(1, 5),
  c3 = c(5, 5)
)

# Simulate all datasets
sim_list <- lapply(beta_pars_list, function(pars) {
  sim_lambda_UV(n_rows = 5e3, beta_pars = pars)
})

# Prior distribution: communalities (density plots)
dens_plots <- mapply(function(pars, nm) {
  ggplot(data.frame(x = seq(0,1,length.out = 500)), aes(x = x)) +
    stat_function(
      fun   = dbeta,
      args  = list(shape1 = pars[1], shape2 = pars[2]),
      geom  = "area",
      fill   = "royalblue1",
      alpha  = 0.3,
      colour = "black",
      size    = 0.8
    ) +
    bayesplot::theme_default(base_size = 16) +
    xlim(0, 1) +
    labs(
      title = bquote(italic(h)^2 ~ "~ Beta(" * .(pars[1]) * ", " * .(pars[2]) * ")"),
      x     = expression(italic(h)^2),
      y     = "Density"
    )
}, beta_pars_list, names(beta_pars_list), SIMPLIFY = FALSE)

# Remove redundant Y axis (and omit the density values)
dens_plots2 <- list(
  dens_plots[[1]] + theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  ),
  dens_plots[[2]] + theme(
    axis.line.y    = element_blank(),
    axis.title.y   = element_blank(),
    axis.text.y    = element_blank(),
    axis.ticks.y   = element_blank()
  ),
  dens_plots[[3]] + theme(
    axis.line.y    = element_blank(),
    axis.title.y   = element_blank(),
    axis.text.y    = element_blank(),
    axis.ticks.y   = element_blank()
  )
)

# Bivariate prior distribution: two factor loadings from two different factors
# loading in one item/experimental effect (scatter plot)
scatter_plots <- lapply(sim_list, function(df) {
  ggplot(df, aes(x = Lambda_1, y = Lambda_2)) +
    geom_point(
      shape  = 21,
      colour = "black",
      fill   = "royalblue1",
      stroke = 0.2,
      alpha  = 0.6
    ) +
    bayesplot::theme_default(base_size = 16) +
    xlim(-1, 1) + ylim(-1, 1) +
    labs(
      x = expression(lambda[11]),
      y = expression(lambda[12])
    )
})

# Remove redundant Y axis
scatter_plots2 <- list(
  scatter_plots[[1]],
  scatter_plots[[2]] + theme(
    axis.line.y    = element_blank(),
    axis.title.y   = element_blank(),
    axis.text.y    = element_blank(),
    axis.ticks.y   = element_blank()
  ),
  scatter_plots[[3]] + theme(
    axis.line.y    = element_blank(),
    axis.title.y   = element_blank(),
    axis.text.y    = element_blank(),
    axis.ticks.y   = element_blank()
  )
)

# Save combined plot
png(filename = "Figures/UV_prior.png", width = 4000, height = 2700, res = 300)
(dens_plots2[[1]] + dens_plots2[[2]] + dens_plots2[[3]]) /
  (scatter_plots2[[1]] + scatter_plots2[[2]] + scatter_plots2[[3]]) +
  plot_layout(heights = c(1, 1))
dev.off()

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: Square-to-circle parametric space
# ─────────────────────────────────────────────────────────────────────────────

# Simulate uniform bivariate distribution
set.seed(2025)
indep_lambdas <- data.frame(
  x = runif(n = 8e3, min = -1, max = 1),
  y = runif(n = 8e3, min = -1, max = 1)
)

# Plausible values assuming that h2 has a maximum value of 1
indep_lambdas$plausible <- rowSums(indep_lambdas^2) <= 1

# Prepare circle data
angle <- seq(0, 2*pi, length.out = 1e4)
circle_df <- data.frame(
  x = cos(angle), 
  y = sin(angle)
)

# Square-to-circle plot
p1 <- ggplot(indep_lambdas, aes(x = x, y = y, color = plausible)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_path(data = circle_df, aes(x = x, y = y), color = "black", linewidth = 1.3) +
  scale_color_manual(
    values = c("FALSE" = "firebrick4", "TRUE" = "royalblue"),
    labels = c(
      expression(italic(h^2) > 1),
      expression(italic(h^2) <= 1)
    ),
    name = expression(italic(h^2) * " values: ")
  ) +
  coord_fixed(xlim = c(-1.05, 1.05), ylim = c(-1.05, 1.05)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  labs(
    x = expression(italic(lambda[11])),
    y = expression(italic(lambda[12]))
  ) +
  bayesplot::theme_default(base_size = 30) + 
  theme(legend.position = "top", 
        legend.title = element_text(size = 30),
        legend.text  = element_text(size = 30)) + 
  guides(
    color = guide_legend(
      override.aes = list(size = 5)  # Adjust point size for legend
    ))
p1

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6: Same direction, different communality values
# ─────────────────────────────────────────────────────────────────────────────

# Define angles
angulos_deg <- seq(45, 45, length.out = 3)
radios      <- c(0.30, 0.60, 0.90)

mag_df <- data.frame(
  r     = radios,
  x     = radios * cos(pi * angulos_deg / 180),
  y     = radios * sin(pi * angulos_deg / 180),
  nivel = factor(sprintf("%.2f", radios), levels = sprintf("%.2f", radios))
)

ratio_arrow <- 0.97
arrow_df <- transform(mag_df,
                      xend = x * ratio_arrow,
                      yend = y * ratio_arrow)

circle_df <- data.frame(
  x = cos(seq(0, 2*pi, length.out = 4000)),
  y = sin(seq(0, 2*pi, length.out = 4000))
)


palette_rb <- c("0.30" = "#8FA8FF",
                "0.60" = "#5A7BEF",
                "0.90" = "#2643B7")
shapes_rb  <- c("0.30" = 16,
                "0.60" = 17,
                "0.90" = 15)
legend_title  <- expression("Magnitude (" * sqrt(italic(h^2)) * "): ")

p2 <- ggplot() +
  geom_path(data = circle_df, aes(x, y), linewidth = 1.3) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  
  # Arrows (hidden from legend)
  geom_segment(data = arrow_df,
               aes(x = 0, y = 0, xend = xend, yend = yend,
                   colour = "grey70"),
               arrow = arrow(length = unit(0.30, "cm")),
               linewidth = 1,
               show.legend = FALSE) +
  
  # Points (shown in legend)
  geom_point(data = mag_df,
             aes(x, y, colour = nivel),
             size = 3) +
  
  scale_colour_manual(values = palette_rb,
                      name   = legend_title,
                      labels = levels(mag_df$nivel)) +
  scale_shape_manual(values = shapes_rb,
                     name   = legend_title,
                     labels = levels(mag_df$nivel)) +
  
  guides(
    colour = guide_legend(override.aes = list(size = 4)),
    shape  = guide_legend(override.aes = list(size = 4))
  ) +
  
  coord_fixed(xlim = c(-1.05, 1.05), ylim = c(-1.05, 1.05)) +
  labs(
    x = expression(italic(lambda[11])),
    y = expression(italic(lambda[12]))
  ) +
  bayesplot::theme_default(base_size = 30) +
  theme(
    legend.position = "top",
    legend.title    = element_text(size = 30),
    legend.text     = element_text(size = 30)
  )

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 7: Different direction, same communality values
# ─────────────────────────────────────────────────────────────────────────────

r_val   <- 0.80
ang_deg <- c(20, 70)
ang_rad <- ang_deg * pi/180

dir_df <- data.frame(
  x      = r_val * cos(ang_rad),
  y      = r_val * sin(ang_rad),
  mag_id = c("m1", "m2")
)

dir_df$magn_val <- sprintf("%.2f", r_val)
ratio_arrow <- 0.98
arrow_df <- transform(dir_df,
                      xend = x * ratio_arrow,
                      yend = y * ratio_arrow)

circle_df <- data.frame(
  x = cos(seq(0, 2*pi, len = 4000)),
  y = sin(seq(0, 2*pi, len = 4000))
)

dir_df$angle_gg <- pi/2 - ang_rad
arc_df <- data.frame(
  start  = pi/2,
  end    = dir_df$angle_gg,
  r      = c(0.22, 0.28),
  mag_id = dir_df$mag_id
)

proj_df <- do.call(rbind, lapply(seq_len(nrow(dir_df)), function(i){
  xi <- dir_df$x[i]; yi <- dir_df$y[i]; mi <- dir_df$mag_id[i]
  data.frame(
    x = c(xi, xi), y = c(yi, yi),
    xend = c(xi, 0), yend = c(0, yi),
    mag_id = mi
  )
}))

off  <- 0.10
dx   <- 0.15

lab_df <- rbind(
  data.frame(
    x = dir_df$x,
    y = -off,
    label = sprintf("lambda[11]==%.2f", dir_df$x),
    mag_id = dir_df$mag_id,
    hjust = 0.5, vjust = 1.2
  ),
  data.frame(
    x = -off - dx,
    y = dir_df$y,
    label = sprintf("lambda[12]==%.2f", dir_df$y),
    mag_id = dir_df$mag_id,
    hjust = 1,  vjust = 0.5
  )
)

pal_colour <- c("m1" = "royalblue", "m2" = "seagreen")
pal_shape  <- c("m1" = 16,          "m2" = 17)
legend_labels <- c("m1" = "0.80", "m2" = "0.80")
legend_title  <- expression("Magnitude (" * sqrt(italic(h^2)) * "): ")

p3 <- ggplot() +
  geom_path(data = circle_df, aes(x, y), linewidth = 1.3) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_segment(data = arrow_df,
               aes(x = 0, y = 0, xend = xend, yend = yend,
                   colour = mag_id),
               arrow = arrow(length = unit(0.30, "cm")),
               linewidth = 1,
               show.legend = FALSE) +
  geom_point(data = dir_df,
             aes(x, y, colour = mag_id, shape = mag_id),
             size = 3) +
  geom_segment(data = proj_df,
               aes(x, y, xend = xend, yend = yend, colour = mag_id),
               linetype = "dotted", linewidth = 0.7,
               show.legend = FALSE) +
  geom_arc(data = arc_df,
           aes(x0 = 0, y0 = 0, r = r, start = start, end = end,
               colour = mag_id),
           linewidth = 1,
           show.legend = FALSE) +
  geom_text(data = lab_df,
            aes(x, y, label = label, colour = mag_id),
            parse = TRUE, size = 6, show.legend = FALSE) +
  scale_colour_manual(values = pal_colour,
                      breaks = names(legend_labels),
                      labels = legend_labels,
                      name   = legend_title) +
  scale_shape_manual(values = pal_shape,
                     breaks = names(legend_labels),
                     labels = legend_labels,
                     name   = legend_title) +
  guides(colour = guide_legend(override.aes = list(size = 4)),
         shape  = guide_legend(override.aes = list(size = 4))) +
  coord_fixed(xlim = c(-1.05, 1.05), ylim = c(-1.05, 1.05)) +
  labs(
    x = expression(italic(lambda[11])),
    y = expression(italic(lambda[12]))
  ) +
  bayesplot::theme_default(base_size = 30) +
  theme(
    legend.position = "top",
    legend.title    = element_text(size = 30),
    legend.text     = element_text(size = 30)
  )

p3

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 8: Final joint plot 
# ─────────────────────────────────────────────────────────────────────────────

# Set combined plot
combined_plot <- p1 + p2 + p3 + plot_layout(ncol = 3)
png(filename = "Figures/appendix_plot.png", width = 7000, height = 3000, res = 300)
print(combined_plot)
dev.off()

# ─────────────────────────────────────────────────────────────────────────────