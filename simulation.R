setwd("/Users/arsene5240/R/")
source("Main_program.r")
load("ex3.health.rda")
## data compilation
ex3.health[, c(1, 5, 6, 7, 11, 12)] = lapply(c(1, 5, 6, 7, 11, 12), function(x) ex3.health[, x] = as.factor(ex3.health[, x]))
ex3.health = ex3.health[-ncol(ex3.health)]


# simulation
hl1 = coretree(prescrib + nonpresc ~ sex + age + income + levyplus + illness + actdays + hscore + chcond1 + doctorco + nondocco + hospadmi + hospdays, zero = "inflated", method = "constant", cdist = "negbin", data = ex3.health)

hl2 = coretree(prescrib + nonpresc ~ sex + age + income + levyplus + illness + actdays + hscore + chcond1 + doctorco + nondocco + hospadmi + hospdays, zero = "inflated", method = "multiple", cdist = "negbin", data = ex3.health)

hl3 = coretree(prescrib + nonpresc ~ sex + age + income + levyplus + illness + actdays + hscore + chcond1 + doctorco + nondocco + hospadmi + hospdays, zero = "inflated", method = "constant", cdist = "poisson", data = ex3.health)

hl4 = coretree(prescrib + nonpresc ~ sex + age + income + levyplus + illness + actdays + hscore + chcond1 + doctorco + nondocco + hospadmi + hospdays, zero = "inflated", method = "multiple", cdist = "poisson", data = ex3.health)

# graph
Names <- c("hl1", "hl2", "hl3", "hl4")
for (i in 1:length(u)) {
  pdf(paste(Names[i], ".pdf", sep = ""))
  print(plot(get(Names[i])))
  dev.off()
  print(i)
}
par(mfrow = c(2,2))
pdf("Merge.pdf", width = 16, height = 12)
print(plot(hl1), split = c(1, 1, 2, 2), more = TRUE)
print(plot(hl2), split = c(1, 2, 2, 2), more = TRUE)
print(plot(hl3), split = c(2, 2, 2, 2), more = TRUE)
print(plot(hl4), split = c(2, 1, 2, 2), more = TRUE)
dev.off()
