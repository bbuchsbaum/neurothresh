# Validate a neurothresh_null_spec

Checks that the object's slots are internally consistent: `draw(1)`
returns a finite numeric vector of length `sum(mask)`. Useful for
factories or tests; production code should not need to call this.

## Usage

``` r
validate_null_spec(ns)
```

## Arguments

- ns:

  A `neurothresh_null_spec`.

## Value

The object invisibly. Throws on failure.
