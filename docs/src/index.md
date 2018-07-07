# [Dendriform.jl](https://github.com/chakravala/Dendriform.jl)

*Symbolic parser generator for Julia language expressions using REDUCE algebra term rewrite system*


```@contents
```

## Setup
Installation of latest release version using Julia:
```Julia
julia> Pkg.add("Dendriform")
```
Provides the types `PBTree` for planar binary trees, `Grove` for tree collections of constant degree, and `GroveBin` to compress grove data. This package defines various essential operations on planar binary trees and groves like `∪` for `union`; `∨` for `graft`; `left` and `right` for branching; `<`, `>`, `≤`, `≥` for Tamari's partial ordering; `/` and `\` (i.e. `over` and `under`); and the `dashv` and `vdash` operations `⊣`, `⊢`, `+`, `*` for dendriform algebra.

## Background

Check the [README](https://github.com/chakravala/Dendriform.jl) for more background information.

## Usage
Basic usage examples:
```Julia
julia> using Dendriform

julia> Grove(3,7) ⊣ [1,2]∪[2,1]
[1,2,5,1,2] ↦ [3]∅∅[2,5][1,4] ↦ 20/42 or 21807
[1,2,5,2,1] ↦ [3]∅∅[2,4][1,5] ↦ 21/42 or 21906
[2,1,5,1,2] ↦ [3]∅∅[1,5][2,4] ↦ 22/42 or 22797
[2,1,5,2,1] ↦ [3]∅∅[1,4][2,5] ↦ 23/42 or 22896
[1,5,3,1,2] ↦ [2]∅[3][5][1,4] ↦ 27/42 or 30807
[1,5,2,1,3] ↦ [2]∅[5][3][1,4] ↦ 25/42 or 29007
[1,5,1,2,3] ↦ [2]∅[5][4][1,3] ↦ 24/42 or 28908
[1,5,3,2,1] ↦ [2]∅[3][4][1,5] ↦ 28/42 or 30906
[1,5,1,3,1] ↦ [2]∅[4]∅[1,3,5] ↦ 26/42 or 30186
267911168 Y5 #9/42 [0.006092%]

julia> Grove(2,3) * [1,2,3]∪[3,2,1] |> GroveBin
2981131286847743360614880957207748817969 Y6 #30/132 [54.75%]

julia> [2,1,7,4,1,3,1] < [2,1,7,4,3,2,1]
true
```
## References
* Dan Yasaki with Adriano Bruno, [The arithmetic of planar binary trees](http://libres.uncg.edu/ir/uncg/f/D_Yasaki_Arithmetic_2011.pdf), Involve 4 (2011), no. 1, 1-11. ([PDF](https://www.uncg.edu/mat/faculty/d_yasaki/publications/trees_for_print.pdf))
* Jean-Louis Loday, [Arithmetree](http://irma.math.unistra.fr/~loday/PAPERS/2002Loday(arithmetree).pdf), J. of Algebra (2002), no. 258, 275-309.
