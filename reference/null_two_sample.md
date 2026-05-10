# Two-sample label-permutation null

Builds an observed two-group t (or Z-equivalent) map and deterministic
label-permutation null draws. Default contrast direction is **level 2
minus level 1**; pin a different direction with `levels` or `contrast`.

## Usage

``` r
null_two_sample(
  images,
  group,
  mask = NULL,
  n_perm = 5000,
  seed = NULL,
  variance = c("welch", "pooled"),
  z_equiv = TRUE,
  template = NULL,
  levels = NULL,
  contrast = NULL
)
```

## Arguments

- images:

  Subject-level images (see
  [`as_subject_matrix()`](https://bbuchsbaum.github.io/neurothresh/reference/as_subject_matrix.md)).

- group:

  Vector of group labels (factor, character, or integer); coerced to a
  2-level factor.

- mask:

  Optional logical 3D mask.

- n_perm:

  Number of label permutations.

- seed:

  Optional random seed.

- variance:

  `"welch"` (default) or `"pooled"`.

- z_equiv:

  If `TRUE` (default), convert t to Z-equivalent.

- template:

  Optional template (only used for matrix input).

- levels:

  Optional length-2 character giving the group level order used for the
  contrast (e.g., `c("control", "case")` → contrast is
  `case - control`). Overrides `as.factor(group)` ordering.

- contrast:

  Optional length-2 character of the form `c("positive", "negative")`
  specifying that the contrast is `mean(positive) - mean(negative)`.
  Equivalent to `levels = rev(contrast)`.

## Value

A
[`null_spec`](https://bbuchsbaum.github.io/neurothresh/reference/null_spec.md).
