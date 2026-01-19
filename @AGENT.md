# neurothresh Build & Test Instructions

## Environment Setup

### Prerequisites

``` bash
# R >= 4.0 required
R --version

# Ensure Xcode CLI tools (macOS) or build-essential (Linux)
# macOS:
xcode-select --install

# Linux:
# sudo apt-get install build-essential
```

### Install Dependencies

``` bash
# Install R package dependencies
Rscript -e "install.packages(c('Rcpp', 'testthat', 'roxygen2', 'devtools', 'knitr', 'rmarkdown', 'remotes'))"

# Install neuroim2 from GitHub
Rscript -e "remotes::install_github('bbuchsbaum/neuroim2')"

# Optional: Install neuroatlas from GitHub
Rscript -e "remotes::install_github('bbuchsbaum/neuroatlas')"

# Install all package dependencies from DESCRIPTION (including Remotes)
Rscript -e "devtools::install_deps()"
```

------------------------------------------------------------------------

## Build Commands

### Generate Documentation

``` bash
Rscript -e "devtools::document()"
```

### Compile C++ Code

``` bash
Rscript -e "Rcpp::compileAttributes()"
Rscript -e "devtools::load_all()"
```

### Build Package

``` bash
Rscript -e "devtools::build()"
```

### Install Package

``` bash
R CMD INSTALL .
# Or:
Rscript -e "devtools::install()"
```

------------------------------------------------------------------------

## Test Commands

### Run All Tests

``` bash
Rscript -e "devtools::test()"
```

### Run Specific Test File

``` bash
Rscript -e "devtools::test(filter='canonicalize')"
```

### Check Test Coverage

``` bash
Rscript -e "covr::package_coverage()"
```

------------------------------------------------------------------------

## Quality Checks

### Full Package Check

``` bash
R CMD check .
```

### CRAN-Style Check

``` bash
R CMD check --as-cran .
```

### Linting (Optional)

``` bash
Rscript -e "lintr::lint_package()"
```

------------------------------------------------------------------------

## Development Workflow

### Load Package for Interactive Development

``` bash
Rscript -e "devtools::load_all()"
```

### Quick Iteration Cycle

``` bash
# After making changes:
Rscript -e "devtools::document(); devtools::load_all(); devtools::test()"
```

### Full Validation Before Commit

``` bash
Rscript -e "devtools::check()"
```

------------------------------------------------------------------------

## Rcpp-Specific Commands

### Regenerate Rcpp Exports

``` bash
Rscript -e "Rcpp::compileAttributes()"
```

This generates: - `src/RcppExports.cpp` - `R/RcppExports.R`

### Compile C++ with Verbose Output

``` bash
R CMD SHLIB src/*.cpp
```

------------------------------------------------------------------------

## Common Issues

### Issue: Rcpp not found

``` bash
Rscript -e "install.packages('Rcpp')"
```

### Issue: neuroim2 not found

``` bash
Rscript -e "devtools::install_local('~/code/neuroim2')"
```

### Issue: C++ compilation errors

1.  Check `src/Makevars` has correct flags
2.  Ensure C++11 or later: `CXX_STD = CXX11`
3.  Run
    [`Rcpp::compileAttributes()`](https://rdrr.io/pkg/Rcpp/man/compileAttributes.html)
    to regenerate exports

### Issue: NAMESPACE missing exports

``` bash
Rscript -e "devtools::document()"
```

------------------------------------------------------------------------

## Completion Checklist

Before marking any task complete:

Code compiles without errors

`devtools::test()` passes for relevant tests

Roxygen2 documentation is complete (@param, @return, @export)

`devtools::check()` passes with no errors

@fix_plan\.md is updated

Before marking project complete:

All items in @fix_plan.md are \[x\]

`R CMD check --as-cran .` passes with 0 errors, 0 warnings

`devtools::test()` shows \>80% coverage

Package installs cleanly via `devtools::install()`
