##############################################

# List input data files

## ATtRACT # Downloaded from https://attract.cnic.es/download
rbp_file = "ATtRACT/ATtRACT_db.txt"
pwm_file = "ATtRACT/pwm.txt"

################################################

# Import requirements
library(data.table)
library(TFBSTools)
source("covid-gene-expression/project/human-RBP-analysis/rbp_functions.R")

################################################

# Read input data

# RBPs
print("Reading RBPs")
rbp = fread(rbp_file, select=c(1:4, 7:9, 11:13))
rbp = unique(rbp)

# PWMs
print("Reading PWMs")
pwm = readPWMsFromFasta(pwm_file)

##############################################

# Count number of RBPs in database
initial_n = rbp[, length(unique(Gene_name))]

# Select only human RBPs
print("Filtering human RBPs")
rbp = rbp[Organism == "Homo_sapiens",]
print(paste0("Reduced number of RBPs from ", initial_n, " to ", rbp[, length(unique(Gene_name))]))

# Remove PWMs which consist of only a single motif
initial_n = rbp[, length(unique(Gene_name))]
rbp = rbp[Score!="1.000000**",]
rbp = unique(rbp[, Score:=NULL])
print(paste0("Reduced number of RBPs from ", initial_n, " to ", rbp[, length(unique(Gene_name))]))

# See experiment types
print("Types of experiments")
print(rbp[, .N, by=Experiment_description][order(N, decreasing = T),])

# See source databases
print("Source databases")
print(rbp[, .N, by=Database][order(N, decreasing = T),])

# How many proteins have multiple PWMs?
print("Number of PWMs per protein")
print(rbp[, length(unique(Matrix_id)), by=Gene_name][order(V1, decreasing = T),])

##############################################

# Select PWMs that match these RBPs
print("Selecting PWMs that match to human RBPs")
initial_n = length(pwm)
pwm = pwm[names(pwm) %in% rbp[, Matrix_id]]
print(paste0("Reduced number of PWMs from ", initial_n, " to ", length(pwm)))

# Calculate entropy of these PWMs
entropy = sapply(pwm, GetPWMEntropy)

# Filter out high-entropy PWMs
plot(density(entropy))
print("Removing high-entropy PWMs")
initial_n = length(pwm)
pwm = pwm[entropy < 8]
print(paste0("Reduced number of PWMs from ", initial_n, " to ", length(pwm)))

# Filter the RBPs corresponding to these PWMs
print("Filtering RBP-PWM matches after removing high-entropy PWMs")
initial_n = rbp[, length(unique(Gene_name))]
rbp = rbp[Matrix_id %in% names(pwm),]
print(paste0("Reduced number of RBPs in the search from ", initial_n, " to ", rbp[, length(unique(Gene_name))]))

# Save filtered data
print("Saving filtered RBP table")
save(rbp, file="output/filtered_rbp.RData")
print("Saving filtered PWMs")
save(pwm, file="output/filtered_pwm.RData")

