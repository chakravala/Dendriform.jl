# GroveAlg
Compute using Loday's dendriform dialgebra based arithmetic on groves of planar binary trees (using this package extending the Julia language). Uses efficient planar binary tree indexing algorithm to compute partial grove indices. 

```
julia> GrovePrint(Grove(3,7)+[1,4,2,1])
[1,2,7,1,4,2,1] ↦ [3]∅∅[5]∅[2,6][1,4,7] ↦ 4128174 or 243/429
[1,2,4,1,7,2,1] ↦ [5]∅∅[3]∅[2,6][1,4,7] ↦ 2328174 or 186/429
[1,2,3,4,7,2,1] ↦ [5]∅∅[4][3][2,6][1,7] ↦ 2221704 or 176/429
[2,1,7,1,4,2,1] ↦ [3]∅∅[5]∅[1,6][2,4,7] ↦ 4138074 or 245/429
[2,1,4,1,7,2,1] ↦ [5]∅∅[3]∅[1,6][2,4,7] ↦ 2338074 or 188/429
[2,1,3,4,7,2,1] ↦ [5]∅∅[4][3][1,6][2,7] ↦ 2222694 or 178/429
[1,7,5,1,4,2,1] ↦ [2]∅[3][5]∅[6][1,4,7] ↦ 5298174 or 292/429
[1,7,2,1,5,2,1] ↦ [2]∅[5]∅∅[3,6][1,4,7] ↦ 5118174 or 278/429
[1,7,1,2,5,2,1] ↦ [2]∅[5]∅∅[4,6][1,3,7] ↦ 5108184 or 276/429
[1,4,2,1,7,2,1] ↦ [5]∅∅[2]∅[3,6][1,4,7] ↦ 2418174 or 192/429
[1,4,1,2,7,2,1] ↦ [5]∅∅[2]∅[4,6][1,3,7] ↦ 2408184 or 190/429
[1,3,1,4,7,2,1] ↦ [5]∅∅[4][2][6][1,3,7] ↦ 2228184 or 180/429
Y7 #12/429, Grove [0.0%] 3978889433292738744351684718788268955154180851751998187580930331583304832228744501395456
```

[![Build Status](https://travis-ci.org/chakravala/GroveAlg.jl.svg?branch=master)](https://travis-ci.org/chakravala/GroveAlg.jl)

[![Coverage Status](https://coveralls.io/repos/chakravala/GroveAlg.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/chakravala/GroveAlg.jl?branch=master)

[![codecov.io](http://codecov.io/github/chakravala/GroveAlg.jl/coverage.svg?branch=master)](http://codecov.io/github/chakravala/GroveAlg.jl?branch=master)
