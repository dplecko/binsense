
to_bits <- function(i, k) as.integer(intToBits(i)[seq_len(k)])

expit <- function(x) exp(x) / (1 + exp(x))

logit <- function(x) log(x / (1 - x))

unit_vec <- function(x) {x / sqrt(sum(x^2))}

vec_norm <- function(x) sqrt(sum(x^2))

mat_cor <- function(A, B) cor(as.vector(A), as.vector(B))
