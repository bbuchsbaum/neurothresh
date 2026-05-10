# GLM Freedman-Lane null

Wraps
[`make_null_fun_glm()`](https://bbuchsbaum.github.io/neurothresh/reference/make_null_fun_glm.md)
(single-coefficient contrasts) and returns a
[`null_spec`](https://bbuchsbaum.github.io/neurothresh/reference/null_spec.md).

## Usage

``` r
null_glm(
  images,
  design,
  contrast,
  mask = NULL,
  data = NULL,
  n_perm = 5000,
  seed = NULL,
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

- n_perm:

  Number of permutations (default 5000).

- seed:

  Optional random seed.

- statistic:

  "t" or "z".

- template:

  Optional template for matrix inputs.

## Value

A
[`null_spec`](https://bbuchsbaum.github.io/neurothresh/reference/null_spec.md).
