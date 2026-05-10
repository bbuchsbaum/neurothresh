# Coerce subject-level images to a canonical (subject by voxel) matrix

Converts heterogeneous subject-level image inputs into a shared
canonical representation used by every subject-level null factory.
Accepts a list of 3D arrays or `NeuroVol`s, a 4D array, a `NeuroVec`, a
`SparseNeuroVec`, or a numeric `matrix` already in subject-by-voxel
form.

## Usage

``` r
as_subject_matrix(x, mask = NULL, template = NULL, ...)
```

## Arguments

- x:

  Input images. One of: list of 3D arrays/`NeuroVol`s, 4D numeric array,
  `NeuroVec`, `SparseNeuroVec`, or a numeric `matrix` with one row per
  subject and one column per mask voxel.

- mask:

  Optional logical 3D mask. When `NULL`, the mask is inferred from the
  finite intersection across subjects (or, for `SparseNeuroVec`, from
  the intrinsic mask). When supplied, the mask is used exactly: any
  voxel inside it that is non-finite in the input is an error.

- template:

  Optional 3D template carrying spatial metadata. Required only for the
  `matrix` branch when `mask` does not already supply dimensions.

- ...:

  Additional arguments (currently unused; reserved for future per-branch
  options).

## Value

An object of S3 class `"neurothresh_image_matrix"` with components:

- X:

  `n_subjects` x `n_mask` numeric matrix.

- mask:

  3D logical array, or `LogicalNeuroVol` when input carried `NeuroSpace`
  metadata.

- mask_idx:

  Integer vector of mask voxel positions in column-major order.

- template:

  3D template volume (array or `NeuroVol`).

- dims:

  Integer length-3 dimensions.

- coords:

  `n_mask` x 3 numeric matrix of array indices.

- space:

  `NeuroSpace` if available, else `NULL`.

- source:

  Character: `"list"`, `"neurovec"`, `"sparseneurovec"`, `"array4d"`, or
  `"matrix"`.
