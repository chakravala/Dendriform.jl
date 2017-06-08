# GroveAlg
 Dendriform dialgebra algorithms to compute using Loday's arithmetic on groves of planar binary trees (package extending the Julia language). Uses efficient planar binary tree indexing algorithm to compute partial grove indices.

[![Build Status](https://travis-ci.org/chakravala/GroveAlg.jl.svg?branch=master)](https://travis-ci.org/chakravala/GroveAlg.jl) [![Build status](https://ci.appveyor.com/api/projects/status/j7t3oc1doeot6i72?svg=true)](https://ci.appveyor.com/project/chakravala/grovealg-jl) [![Coverage Status](https://coveralls.io/repos/github/chakravala/GroveAlg.jl/badge.svg?branch=master)](https://coveralls.io/github/chakravala/GroveAlg.jl?branch=master) [![codecov.io](http://codecov.io/github/chakravala/GroveAlg.jl/coverage.svg?branch=master)](http://codecov.io/github/chakravala/GroveAlg.jl?branch=master)

Installation of latest release version using Julia:
```Julia
julia> Pkg.add("GroveAlg")
```
Usage example:
```Julia
julia> using GroveAlg

julia> Grove(3,7) ⊣ Grove(2,3)
[1,2,5,1,2] ↦ [3]∅∅[2,5][1,4] ↦ 20/42 or 510696
[1,2,5,2,1] ↦ [3]∅∅[2,4][1,5] ↦ 21/42 or 510795
[2,1,5,1,2] ↦ [3]∅∅[1,5][2,4] ↦ 22/42 or 511686
[2,1,5,2,1] ↦ [3]∅∅[1,4][2,5] ↦ 23/42 or 511785
[1,5,3,1,2] ↦ [2]∅[3][5][1,4] ↦ 27/42 or 519696
[1,5,2,1,3] ↦ [2]∅[5][3][1,4] ↦ 25/42 or 517896
[1,5,1,2,3] ↦ [2]∅[5][4][1,3] ↦ 24/42 or 517797
[1,5,3,2,1] ↦ [2]∅[3][4][1,5] ↦ 28/42 or 519795
[1,5,1,3,1] ↦ [2]∅[4]∅[1,3,5] ↦ 26/42 or 519075
267911168 Y5 #9/42 [0.006092%]

julia> Grove(2,3) * Grove(3,17) |> GroveBin
2981131286847743360614880957207748817969 Y6 #30/132 [54.75%]
```
Provides the types `PBTree` for planar binary trees, `Grove` for degree-invariant tree collections, and `GroveBin` to compress grove data. This package defines various essential operations on planar binary trees and groves like `∪` for union, `∨` for grafting, `leftbranch` and `rightbranch` for branching, `↗` and `↖` (i.e. `over` and `under`), and operators `⊣`, `⊢`, `+`, `*` for dendriform algebra.
