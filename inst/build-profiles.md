# Compilation profiles

| Profile | Flags | Status |
|---|---|---|
| Conservative (default) | R's `-O2` / `$(SHLIB_OPENMP_CXXFLAGS)`, no `-march=native`, no `-ffast-math` | Certified by the package test suite. Use for CRAN-like checks, sanitizers, and release artifacts. Shipped `src/Makevars` is POSIX-portable. |
| Fast (opt-in) | `RFUGW_FAST_FLAGS="-march=native -ffast-math -fno-math-errno -fno-trapping-math"` | Machine-local. Can change IEEE numerics. Not used for 0.1 certification. |

Serial builds are the fallback when OpenMP flags are empty. OpenMP and serial
paths must remain equivalent within the solver-contract residuals.

On macOS Homebrew, set these if R's `SHLIB_OPENMP_*` flags are not enough:

```
RFUGW_OPENMP_FLAGS="-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include"
RFUGW_OPENMP_LIBS="-L/opt/homebrew/opt/libomp/lib -lomp"
```

Accuracy and wall-time deltas between profiles are recorded by
`inst/bench/PROTOCOL.md`. Only `profile=conservative` rows may update a
0.1 baseline. `RFUGW_FAST_FLAGS` runs must be labeled `profile=fast`.
