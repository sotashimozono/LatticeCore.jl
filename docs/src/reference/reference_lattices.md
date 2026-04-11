# Reference lattices API

LatticeCore ships two minimal concrete lattices so the interface
can be exercised end-to-end without any downstream package. These
are *reference implementations*, not production lattices — the full
2D catalogue (triangular, honeycomb, kagome, ...) lives in
`Lattice2D.jl`.

```@docs
LineLattice
SimpleSquareLattice
basis_vectors
```
