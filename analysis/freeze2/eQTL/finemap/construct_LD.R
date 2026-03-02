






library(data.table)



args <- commandArgs(trailingOnly = TRUE)

variant_file=args[1]
ld_file=args[2]
output_file=args[3]


vars <- fread(variants_file, header = FALSE)[[1]]
N <- length(vars)

# 2. read PLINK LD output
ld <- fread(ld_file)

# 3. initialize matrix
ld_mat <- matrix(NA_real_, N, N,
                 dimnames = list(vars, vars))







# fill matrix
ld_mat[cbind(ld$SNP_A, ld$SNP_B)] <- ld$R2
ld_mat[cbind(ld$SNP_B, ld$SNP_A)] <- ld$R2

diag(ld_mat) <- 1



write.table(ld_mat,file=ld_file,quote = FALSE, sep = "\t", eol = "\n", na = NA, dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"))
