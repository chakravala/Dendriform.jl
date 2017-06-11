# Dendriform.jl
 Dendriform dialgebra algorithms to compute using Loday's arithmetic on groves of planar binary trees (package extending the Julia language). Uses efficient planar binary tree indexing algorithm to compute partial grove indices.

[![Build Status](https://travis-ci.org/chakravala/Dendriform.jl.svg?branch=master)](https://travis-ci.org/chakravala/Dendriform.jl) [![Build status](https://ci.appveyor.com/api/projects/status/j7t3oc1doeot6i72?svg=true)](https://ci.appveyor.com/project/chakravala/grovealg-jl) [![Coverage Status](https://coveralls.io/repos/github/chakravala/Dendriform.jl/badge.svg?branch=master)](https://coveralls.io/github/chakravala/Dendriform.jl?branch=master) [![codecov.io](http://codecov.io/github/chakravala/Dendriform.jl/coverage.svg?branch=master)](http://codecov.io/github/chakravala/Dendriform.jl?branch=master)

Installation of latest release version using Julia:
```Julia
julia> Pkg.add("Dendriform")
```
Usage example:
```Julia
julia> using Dendriform

julia> Grove(3,7) ⊣ ∪([1,2],[2,1])
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

julia> Grove(2,3) * ∪([1,2,3],[3,2,1]) |> GroveBin
2981131286847743360614880957207748817969 Y6 #30/132 [54.75%]
```
Provides the types `PBTree` for planar binary trees, `Grove` for degree-invariant tree collections, and `GroveBin` to compress grove data. This package defines various essential operations on planar binary trees and groves like `∪` for `union`, `∨` for `graft`, `leftbranch` and `rightbranch` for branching, `↗` and `↖` (i.e. `over` and `under`), and the `leftsum` and `rightsum` operators `⊣`, `⊢`, `+`, `*` for dendriform algebra.
