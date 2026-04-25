# ============================================================
# Plot Figure 2 - Both R and EViews results
# ============================================================

# Load R results
r_gdp <- read.csv("est_1_R.csv")
r_pce <- read.csv("est_2_R.csv")
r_ffr <- read.csv("est_3_R.csv")

# Load original EViews results - fix column names
orig_gdp <- read.csv("~/Documents/Brown/Classes/DATA2020/Final Project/DATA_2020_Project/Figure 2/Fig_2/est_1.csv",
                     header = FALSE, skip = 1)
orig_pce <- read.csv("~/Documents/Brown/Classes/DATA2020/Final Project/DATA_2020_Project/Figure 2/Fig_2/est_2.csv",
                     header = FALSE, skip = 1)
orig_ffr <- read.csv("~/Documents/Brown/Classes/DATA2020/Final Project/DATA_2020_Project/Figure 2/Fig_2/est_3.csv",
                     header = FALSE, skip = 1)

# Assign correct column names
col_names <- c("point_linear","point_high","point_low",
               "linear1","linear2","linear3",
               "H1","H2","H3","L1","L2","L3")

colnames(orig_gdp) <- col_names
colnames(orig_pce) <- col_names
colnames(orig_ffr) <- col_names

# Keep only 20 rows (h=0 to h=19, last row has NA in original)
orig_gdp <- orig_gdp[1:20, ]
orig_pce <- orig_pce[1:20, ]
orig_ffr <- orig_ffr[1:20, ]

r_gdp <- r_gdp[1:20, ]
r_pce <- r_pce[1:20, ]
r_ffr <- r_ffr[1:20, ]

# Horizon axis
h <- 0:19
zero <- rep(0, 20)

# ============================================================
# Plot function
# ============================================================

plot_irf <- function(r_data, orig_data, title, ylim_range) {
  
  par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
  
  # --- Row 1: Original EViews results ---
  
  # Panel 1: Three Models (original)
  plot(h, orig_data$point_linear,
       type = "o", pch = 1, lwd = 1.5, col = "black",
       ylim = ylim_range, xlim = c(0, 19),
       xlab = "Horizon (Quarter)", ylab = "",
       main = "Three Models (EViews)", cex.main = 1)
  lines(h, orig_data$point_high, type = "o", pch = 2, lwd = 1.5, col = "blue")
  lines(h, orig_data$point_low,  type = "o", pch = 0, lwd = 1.5, col = "red")
  abline(h = 0, lwd = 0.5)
  legend("topleft", legend = c("Linear","High","Low"),
         col = c("black","blue","red"),
         pch = c(1,2,0), lwd = 1.5, bty = "n", cex = 0.8)
  grid()
  
  # Panel 2: Linear Model (original)
  plot(h, orig_data$linear2,
       type = "o", pch = 1, lwd = 1.5, col = "black",
       ylim = ylim_range, xlim = c(0, 19),
       xlab = "Horizon (Quarter)", ylab = "",
       main = "Linear Model (EViews)", cex.main = 1)
  polygon(c(h, rev(h)),
          c(orig_data$linear3, rev(orig_data$linear1)),
          col = rgb(0.85,0.85,0.85,0.5), border = NA)
  lines(h, orig_data$linear2, type = "o", pch = 1, lwd = 1.5)
  lines(h, orig_data$linear1, lty = 2, lwd = 1)
  lines(h, orig_data$linear3, lty = 2, lwd = 1)
  abline(h = 0, lwd = 0.5)
  grid()
  
  # Panel 3: State-Dependent (original)
  plot(h, orig_data$H2,
       type = "o", pch = 2, lwd = 1.5, col = "blue",
       ylim = ylim_range, xlim = c(0, 19),
       xlab = "Horizon (Quarter)", ylab = "",
       main = "State-Dependent (EViews)", cex.main = 1)
  polygon(c(h, rev(h)),
          c(orig_data$H3, rev(orig_data$H1)),
          col = rgb(0.85,0.85,1,0.5), border = NA)
  lines(h, orig_data$H2, type = "o", pch = 2, lwd = 1.5, col = "blue")
  lines(h, orig_data$H1, lty = 2, lwd = 1, col = "blue")
  lines(h, orig_data$H3, lty = 2, lwd = 1, col = "blue")
  lines(h, orig_data$L2, type = "o", pch = 0, lwd = 1.5, col = "red")
  lines(h, orig_data$L1, lty = 3, lwd = 2, col = "red")
  lines(h, orig_data$L3, lty = 3, lwd = 2, col = "red")
  abline(h = 0, lwd = 0.5)
  grid()
  
  # --- Row 2: R results ---
  
  # Panel 4: Three Models (R)
  plot(h, r_data$point_linear,
       type = "o", pch = 1, lwd = 1.5, col = "black",
       ylim = ylim_range, xlim = c(0, 19),
       xlab = "Horizon (Quarter)", ylab = "",
       main = "Three Models (R)", cex.main = 1)
  lines(h, r_data$point_high, type = "o", pch = 2, lwd = 1.5, col = "blue")
  lines(h, r_data$point_low,  type = "o", pch = 0, lwd = 1.5, col = "red")
  abline(h = 0, lwd = 0.5)
  legend("topleft", legend = c("Linear","High","Low"),
         col = c("black","blue","red"),
         pch = c(1,2,0), lwd = 1.5, bty = "n", cex = 0.8)
  grid()
  
  # Panel 5: Linear Model (R)
  plot(h, r_data$linear2,
       type = "o", pch = 1, lwd = 1.5, col = "black",
       ylim = ylim_range, xlim = c(0, 19),
       xlab = "Horizon (Quarter)", ylab = "",
       main = "Linear Model (R)", cex.main = 1)
  polygon(c(h, rev(h)),
          c(r_data$linear3, rev(r_data$linear1)),
          col = rgb(0.85,0.85,0.85,0.5), border = NA)
  lines(h, r_data$linear2, type = "o", pch = 1, lwd = 1.5)
  lines(h, r_data$linear1, lty = 2, lwd = 1)
  lines(h, r_data$linear3, lty = 2, lwd = 1)
  abline(h = 0, lwd = 0.5)
  grid()
  
  # Panel 6: State-Dependent (R)
  plot(h, r_data$H2,
       type = "o", pch = 2, lwd = 1.5, col = "blue",
       ylim = ylim_range, xlim = c(0, 19),
       xlab = "Horizon (Quarter)", ylab = "",
       main = "State-Dependent (R)", cex.main = 1)
  polygon(c(h, rev(h)),
          c(r_data$H3, rev(r_data$H1)),
          col = rgb(0.85,0.85,1,0.5), border = NA)
  lines(h, r_data$H2, type = "o", pch = 2, lwd = 1.5, col = "blue")
  lines(h, r_data$H1, lty = 2, lwd = 1, col = "blue")
  lines(h, r_data$H3, lty = 2, lwd = 1, col = "blue")
  lines(h, r_data$L2, type = "o", pch = 0, lwd = 1.5, col = "red")
  lines(h, r_data$L1, lty = 3, lwd = 2, col = "red")
  lines(h, r_data$L3, lty = 3, lwd = 2, col = "red")
  abline(h = 0, lwd = 0.5)
  grid()
  
  # Main title
  mtext(title, outer = TRUE, cex = 1.2, line = -1.5)
}

# ============================================================
# Generate the three figures
# ============================================================

# GDP
dev.new()
plot_irf(r_gdp, orig_gdp,
         title = "Figure 2 - GDP",
         ylim_range = c(-0.5, 1.5))

# PCE
dev.new()
plot_irf(r_pce, orig_pce,
         title = "Figure 2 - PCE Deflator",
         ylim_range = c(-0.5, 2.0))

# FFR
dev.new()
plot_irf(r_ffr, orig_ffr,
         title = "Figure 2 - FFR",
         ylim_range = c(-2, 1))