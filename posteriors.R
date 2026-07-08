# args = commandArgs(trailingOnly=TRUE)
# 
# if (length(args) == 0) {
#   cat("Usage: posteriors.R <target_directory>")
#   quit(save="no", status=1)
# }
# 
# current = getwd()
# target_dir = args[1]
# target_dir = "hd74156_12"
# new_dir = file.path(current, "runs", target_dir)
# setwd(new_dir)

setwd("/Users/sogsbilby/Honours/code/pipeline/runs/hd74156_12")



data = read.table("data.txt")
num_instr = max(data$V4) + 1
post = read.table("posterior_sample.txt")
head(post)

period_cols = c(2, 7) + num_instr
ecc_cols = c(4, 9) + num_instr
omega_cols = c(5, 10) + num_instr

combined_periods = unlist(post[, period_cols])
combined_ecc = unlist(post[, ecc_cols])
combined_omega = unlist(post[, omega_cols])

if (!dir.exists("images")) dir.create("images")

hist(combined_periods, breaks=100, main="Periods")
hist(combined_ecc, breaks=100, main="Eccentricities")
hist(combined_omega, breaks=100, main="Arguments of Periastron")

# png("images/period_hist.png", width=800, height=600)
# hist(combined_periods, breaks=100, main="Periods")
# png("images/eccentricity_hist.png", width=800, height=600)
# hist(combined_ecc, breaks=100, main="Eccentricities")
# dev.off()
# png("images/periastron_hist.png", width=800, height=600)
# hist(combined_omega, breaks=100, main="Arguments of Periastron")
# dev.off()








