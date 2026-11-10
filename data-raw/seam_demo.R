# Run once from the package root to create data/seam_demo.rda.
set.seed(20261023)
n = 48
p = 12
ids = paste0('gene', sprintf('%02d', seq_len(2 * p)))
module_assignment = rep(c('blue', 'red'), each = p)
make_group = function(group){
  z1 = as.numeric(scale(rnorm(n)))
  z2 = as.numeric(scale(rnorm(n)))
  noise = matrix(rnorm(n * 2 * p, sd = 0.65), n, 2 * p)
  # Different dominant variables between groups; shared weaker effects.
  a = if(group == 'A') c(rep(1.6, 5), rep(0.35, 4), rep(0.05, 3)) else
    c(rep(0.15, 4), rep(1.7, 5), rep(0.05, 3))
  b = if(group == 'A') c(rep(1.5, 5), rep(0.35, 4), rep(0.05, 3)) else
    c(rep(0.15, 4), rep(1.6, 5), rep(0.05, 3))
  x = noise + tcrossprod(z1, c(a, rep(0, p))) +
    tcrossprod(z2, c(rep(0, p), b))
  colnames(x) = ids
  rownames(x) = paste0(group, sprintf('%02d', seq_len(n)))
  # Group B reverses the first trait's direction relative to blue module.
  traits = data.frame(
    trait_blue = (if(group == 'A') 1 else -1) * z1 + rnorm(n, sd = 0.3),
    trait_red = z2 + rnorm(n, sd = 0.3),
    row.names = rownames(x)
  )
  list(expression = x, traits = traits)
}
a = make_group('A')
b = make_group('B')
seam_demo = list(
  expression_data = list(A = a$expression, B = b$expression),
  trait_data = list(A = a$traits, B = b$traits),
  module_assignment = module_assignment,
  labels = c('A', 'B')
)
dir.create('data', showWarnings = FALSE)
save(seam_demo, file = 'data/seam_demo.rda', compress = 'xz')
