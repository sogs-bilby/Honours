# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# args = commandArgs(trailingOnly=TRUE)
# 
# #run_dir = args[1]
# 
# run_dir = "hd74156_12"
# 
# data = read.table(sprintf("runs/%s/data.txt", run_dir))
# post = read.table(sprintf("runs/%s/posterior_sample.txt", run_dir))
# 
# solve_kepler = function(M, ecc) {
#   M = M %% (2 * pi)
#   E = M
#   E[ecc > 0.8] = pi
#   
#   for (i in 1:50) {
#     f = E - ecc * sin(E) - M
#     f_prime = 1 - ecc * cos(E)
#     E_next = E - f / f_prime
#     
#     if (max(abs(E_next - E)) < 1e-10) {
#       return (E_next)
#     }
#     E = E_next
#   }
# }
# calc_rv_model = function(t, gamma, K, T_per, phi, ecc, omega) {
#   # mean anom
#   M_anom = (2 * pi * t / T_per) + phi
#   # eccentric anom
#   E_anom = solve_kepler(M_anom, ecc)
#   # true anom
#   f_true = 2 * atan(sqrt((1 + ecc) / (1 - ecc)) * tan(E_anom / 2))
#   # keplerian RV curve
#   mu = gamma + K * (cos(omega + f_true) + ecc * cos(omega))
#   
#   return(mu)
# }
# 
# t_min <- min(data$V1) - 10
# t_max <- max(data$V1) + 10
# t_grid <- seq(t_min, t_max, length.out = 1000)
# 
# png("test.png")
# 
# plot(data$V1, data$V2, 
#      pch = 16, col = "black", 
#      xlab = "time", ylab = "RV",
#      main = "HD 73526",
#      ylim = range(c(data$V2 - data$V3, data$V2 + data$V3)))
# 
# if ("V3" %in% colnames(data)) {
#   arrows(data$V1, data$V2 - data$V3, data$V1, data$V2 + data$V3, 
#          length = 0.05, angle = 90, code = 3, col = "gray50")
# }
# 
# num_samples_to_plot = 100
# sample_indices = sample(1:nrow(post), num_samples_to_plot)
# 
# for (i in sample_indices) {
#   # 1 posterior sample
#   p = post[i, ]
#   # get model line
#   rv_curve = calc_rv_model(t_grid, 
#                            p$V1, p$V2, p$V3, 
#                            p$V4, p$V5, p$V6)
#   # add line
#   lines(t_grid, rv_curve, col = rgb(0, 0.4, 0.8, alpha = 0.05), lwd = 1)
# }
# 
# dev.off()
# 







setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
args = commandArgs(trailingOnly=TRUE)

run_dir = "hd74156_12"

# Load data and posterior samples
data = read.table(sprintf("runs/%s/data.txt", run_dir))
post = read.table(sprintf("runs/%s/posterior_sample.txt", run_dir))

# Dynamic detection of instruments based on the 4th column of your data file
num_instruments <- max(data$V4) + 1 

# ==============================================================================
# POSTERIOR COLUMN MAP ASSUMPTION
# Verify this matches your C++ output exactly!
# Columns 1 to num_instruments          = gamma_0, gamma_1, ...
# Columns inst+1 to inst+5               = K1, P1, phi1, ecc1, omega1 (Planet 1)
# Columns inst+6 to inst+10              = K2, P2, phi2, ecc2, omega2 (Planet 2)
# ==============================================================================

# 1. Calculate median gamma offsets to shift the data points to a common frame
median_gammas <- c()
for (inst in 0:(num_instruments - 1)) {
  median_gammas[inst + 1] <- median(post[[paste0("V", inst + 1)]])
}

# 2. Create shifted RV data column centered around 0
data$V2_shifted <- data$V2
for (i in 1:nrow(data)) {
  inst_idx <- data$V4[i]
  data$V2_shifted[i] <- data$V2[i] - median_gammas[inst_idx + 1]
}

# Kepler Solver (unchanged)
solve_kepler = function(M, ecc) {
  M = M %% (2 * pi)
  E = M
  E[ecc > 0.8] = pi
  
  for (i in 1:50) {
    f = E - ecc * sin(E) - M
    f_prime = 1 - ecc * cos(E)
    E_next = E - f / f_prime
    
    if (max(abs(E_next - E)) < 1e-10) {
      return (E_next)
    }
    E = E_next
  }
}

# Zero-centered Keplerian RV curve function for an individual planet
calc_planet_rv = function(t, K, T_per, phi, ecc, omega) {
  M_anom = (2 * pi * t / T_per) + phi
  E_anom = solve_kepler(M_anom, ecc)
  f_true = 2 * atan(sqrt((1 + ecc) / (1 - ecc)) * tan(E_anom / 2))
  
  # Isolated Keplerian velocity contribution
  mu = K * (cos(omega + f_true) + ecc * cos(omega))
  return(mu)
}

# Set up time grid for smooth plotting
t_min <- min(data$V1) - 10
t_max <- max(data$V1) + 10
t_grid <- seq(t_min, t_max, length.out = 1000)

#png("test.png", width = 900, height = 600, res = 110)

# Define instrument colors
palette_colors <- c("firebrick", "dodgerblue3", "forestgreen", "darkorchid")
point_colors <- palette_colors[data$V4 + 1]

# Plot the shifted data points
# plot(data$V1, data$V2_shifted, 
#      pch = 16, col = point_colors, 
#      xlab = "Time (BJD)", ylab = "Offset-Corrected RV (m/s)",
#      main = "HD 74156 - Two-Planet Keplerian Fit",
#      ylim = range(c(data$V2_shifted - data$V3, data$V2_shifted + data$V3)))
plot(data$V1, data$V2_shifted, 
     pch = 16,
     xlab = "BJD", ylab = "RV",
     main = "HD74156 fit for all data",
     ylim = range(c(data$V2_shifted - data$V3, data$V2_shifted + data$V3)))

# Add error bars
if ("V3" %in% colnames(data)) {
  arrows(data$V1, data$V2_shifted - data$V3, data$V1, data$V2_shifted + data$V3, 
         length = 0.04, angle = 90, code = 3, col = "gray60")
}

# legend("topright", 
#        legend = paste("Instrument", 0:(num_instruments-1)), 
#        col = palette_colors[1:num_instruments], 
#        pch = 16)

# 3. Draw the combined two-planet model lines from posterior samples
num_samples_to_plot = 100
sample_indices = sample(1:nrow(post), num_samples_to_plot)

for (i in sample_indices) {
  p = post[i, ]
  
  # Extract Planet 1 Parameters
  K1     <- p[[paste0("V", num_instruments + 1)]]
  P1     <- p[[paste0("V", num_instruments + 2)]]
  phi1   <- p[[paste0("V", num_instruments + 3)]]
  ecc1   <- p[[paste0("V", num_instruments + 4)]]
  omega1 <- p[[paste0("V", num_instruments + 5)]]
  
  # Extract Planet 2 Parameters
  K2     <- p[[paste0("V", num_instruments + 6)]]
  P2     <- p[[paste0("V", num_instruments + 7)]]
  phi2   <- p[[paste0("V", num_instruments + 8)]]
  ecc2   <- p[[paste0("V", num_instruments + 9)]]
  omega2 <- p[[paste0("V", num_instruments + 10)]]
  
  # Calculate isolated signals for both planets
  rv_planet1 = calc_planet_rv(t_grid, K1, P1, phi1, ecc1, omega1)
  rv_planet2 = calc_planet_rv(t_grid, K2, P2, phi2, ecc2, omega2)
  
  # Superposition: Sum them up!
  rv_combined = rv_planet1 + rv_planet2
  
  # Draw the translucent combined model line
  lines(t_grid, rv_combined, col = rgb(0, 0.4, 0.8, alpha = 0.04), lwd = 1)
}

#dev.off()




















