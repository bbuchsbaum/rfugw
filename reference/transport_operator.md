# Construct an Implicit Transport Operator

Creates a coupling-like linear operator when entries need not be stored.
Both forward and adjoint callbacks are mandatory; supplied row and
column masses define barycentric abstention. `adjoint` is the
mathematical transpose action, not an inverse or reverse conditional. A
materializer is optional and is called only by explicit
[`transport_plan_materialize()`](https://bbuchsbaum.github.io/rfugw/reference/transport_plan_materialize.md).

## Usage

``` r
transport_operator(
  n_source,
  n_target,
  apply,
  adjoint,
  source_mass,
  target_mass,
  materialize = NULL,
  metadata = list()
)
```

## Arguments

- n_source, n_target:

  Positive operator shape.

- apply, adjoint:

  Functions accepting a numeric destination/source matrix and returning
  the forward/adjoint product.

- source_mass, target_mass:

  Finite nonnegative row and column masses with equal totals.

- materialize:

  Optional zero-argument function returning a dense or sparse
  nonnegative plan of the declared shape.

- metadata:

  Serializable provenance list.

## Value

An implicit `rfugw_transport_plan`.
