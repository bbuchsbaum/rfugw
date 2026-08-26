# Construct a Transport Plan Representation

Wraps a dense numeric matrix, a
[`Matrix::sparseMatrix`](https://rdrr.io/pkg/Matrix/man/sparseMatrix.html),
or a 1-based edge list without changing its mathematical coupling. Edge
lists are sorted by source then target. Duplicate edges are summed by
default (or rejected), and resulting zero-weight entries are removed.
Missing rows or columns are valid and have explicit zero transported
mass.

## Usage

``` r
as_transport_plan(
  x,
  n_source = NULL,
  n_target = NULL,
  duplicates = c("sum", "error")
)

# S3 method for class 'rfugw_transport_plan'
dim(x)

# S3 method for class 'rfugw_transport_plan'
as.matrix(x, ...)

# S3 method for class 'rfugw_transport_plan'
print(x, ...)
```

## Arguments

- x:

  Dense matrix,
  [`Matrix::sparseMatrix`](https://rdrr.io/pkg/Matrix/man/sparseMatrix.html),
  edge data frame/list, existing transport plan, or result containing a
  plan. Edge columns may be named `source,target,weight` or `i,j,x`.

- n_source, n_target:

  Required positive shape for an edge list; inferred from matrix inputs.

- duplicates:

  Edge-list duplicate policy: canonicalize by summing or reject.

- ...:

  Ignored by the matrix-conversion and print methods.

## Value

A serializable `rfugw_transport_plan`.

## Examples

``` r
edges <- data.frame(source = c(2, 1, 1), target = c(1, 2, 2),
                    weight = c(0.4, 0.2, 0.3))
plan <- as_transport_plan(edges, 2, 2)
transport_plan_mass(plan, "source")
#> [1] 0.5 0.4
```
