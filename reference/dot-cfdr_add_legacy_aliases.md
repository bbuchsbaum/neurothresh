# Attach legacy aliases to a cluster_fdr result

Allows older callers using `$sig_mask`, `$table`, `$n_clusters`, and
`$cluster_labels` to keep working.

## Usage

``` r
.cfdr_add_legacy_aliases(out, clusters, table)
```
