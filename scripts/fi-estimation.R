
library(data.table)
tbl <- read.csv("~/Desktop/COPD_MA.csv")

tbl <- as.data.table(tbl)

pop_prv <- tbl[Response == "No", c("Data_Value", "Stratification")]
pop_prv <- setnames(pop_prv, c("Data_Value", "Stratification"), c("prev", "age"))

pop_prv[, prev := 1 - as.numeric(prev) / 100]

# get the proportions of age groups in the data



# get the target prevalence
pop_prv <- cbind(
  pop_prv, 
  cnt = table(.bincode(data$age, c(-Inf, 24.5, 34.5, 44.5, 54.5, 64.5, Inf)))
)
pop_prv[, prev_miiv := cnt.N / sum(cnt.N)]
sum(pop_prv$prev * pop_prv$prev_miiv)

# get the true prevalence
mean((data$) > 0)

# get COPD prevalence
nrow(load_concepts("copd", "miiv", id_type = "patient",
                   patient_ids = data$subject_id)) / length(data$subject_id)
