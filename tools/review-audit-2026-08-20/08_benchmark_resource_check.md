# Benchmark-resource clean-package verification

Date: 2026-08-20  
Base commit: `5aef129cdb765a51f9b56b25b66fd2fd1a5f6c16`  
Patch under test: only `tests/testthat/test-bench-protocol.R` plus
`tests/testthat/helper-bench-resource.R`

The source-relative lookup reproduced by the audit was replaced with
`system.file("bench", ..., package = "rfugw")`; an explicit source-checkout
fallback is used only when no installed resource exists.

## Content identity

- `test-bench-protocol.R` SHA-256:
  `6c7dd5f52aa9a6549d193cf215773f7877b329e13ea335ad2213f8855fe3eef2`
- `helper-bench-resource.R` SHA-256:
  `28bc752c7a1955503884204477c2f852db8c645b1c442141b1f0956b25deb483`
- Built tarball SHA-256:
  `8a07c7b201ac82543197a0c3b41ea896201c9a9e165be15347e5fd65e8d9e4f9`

## Verification

Source-tree test:

```text
devtools::test(filter = "partial-asymmetric|bench-protocol")
bench-protocol: 34 passing expectations
```

Isolated build and installed-package check:

```text
LC_ALL=C LANG=C R CMD build src --no-manual
LC_ALL=C LANG=C R CMD check --as-cran --no-manual rfugw_0.0.1.9000.tar.gz
Status: 2 NOTEs, 0 WARNINGs, 0 ERRORs
tests/testthat.R: OK
vignette rebuild: OK
```

The check ran from `/tmp/rfugw-bench-check.hxOjDM`; the full check log was
`/tmp/rfugw-bench-check.hxOjDM/rfugw.Rcheck/00check.log`. The directory is an
ephemeral local artifact and is not publication or hosted-CI evidence.

Environment:

- R 4.5.1, arm64-apple-darwin20, macOS Sonoma 14.3
- package compiler: Homebrew clang 20.1.8
- SDK: MacOSX14.2.sdk
- locale: `LC_ALL=C`, `LANG=C`, session charset ASCII

The two reviewed notes are the development version/new-submission note and the
existing non-standard top-level `CONTRIBUTING.md`/`_pkgdown.yml` note. Neither
is caused by the benchmark-resource patch.
