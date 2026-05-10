# GLM Freedman-Lane Null Generator

Builds an observed GLM contrast statistic map and deterministic
Freedman-Lane residual-permutation null draws.

## Usage

``` r
make_null_fun_glm(
  images,
  design,
  contrast,
  mask = NULL,
  data = NULL,
  permutation = permutation_glm(),
  statistic = c("t", "z"),
  template = NULL
)
```

## Arguments

- images:

  List of subject-level contrast images or matrix.

- design:

  Formula, model matrix, or data frame. If a formula is supplied, `data`
  must contain the variables.

- contrast:

  Named coefficient or numeric contrast vector.

- mask:

  Optional logical mask.

- data:

  Optional data for formula designs.

- permutation:

  A
  [`permutation_glm`](https://bbuchsbaum.github.io/neurothresh/reference/permutation_glm.md)
  object.

- statistic:

  "t" or "z".

- template:

  Optional template for matrix inputs.
