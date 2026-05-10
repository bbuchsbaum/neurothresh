# Construct Within-Subject Contrast Images

Builds one planned repeated-measures contrast image per subject. The
first implementation supports a complete 2 x 2 interaction contrast.

## Usage

``` r
make_within_subject_contrast_images(images, subject, within, contrast)
```

## Arguments

- images:

  List of 3D images, one per row of the design vectors.

- subject:

  Subject identifier for each image.

- within:

  Data frame or list containing within-subject factor columns.

- contrast:

  Contrast string. Currently supports "factor1:factor2" for two 2-level
  within-subject factors.

## Value

List of contrast images, one per subject.
