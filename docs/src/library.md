# Dendriform.jl Library

```@contents
```

## Index

```@index
```

## Dendriform Types

```@docs
PBTree
```

```@docs
Grove
```

```@docs
GroveBin
```

```@docs
Dendriform.BaseTree
```

## Dendriform Operators

```@docs
Base.∪
```

### Branching

```@docs
graft
```

```@docs
∨
```

```@docs
left
```

```@docs
right
```

### Partial Ordering

```@docs
Dendriform.posetnext
```

```@docs
Dendriform.posetprev
```

```@docs
Base.<
```

```@docs
Base.>
```

```@docs
Base.≤
```

```@docs
Base.≥
```

### Dialgebra Arithmetic

```@docs
over
```

```@docs
under
```

```@docs
↗
```

```@docs
↖
```

```@docs
dashv
```

```@docs
vdash
```

```@docs
⊣
```

```@docs
⊢
```

```@docs
Base.+
```

```@docs
Base.*
```

## Dendriform Morphisms

```@docs
σ
```

### Consistency Checks

```@docs
Dendriform.LeftInherited
```

```@docs
Dendriform.RightInherited
```

```@docs
Dendriform.PrimitiveTree
```

```@docs
treecheck
```

```@docs
grovecheck
```

```@docs
Dendriform.GroveError
```

### Index Data

```@docs
treeindex
```

```@docs
treeindexCn
```

```@docs
groveindex
```

```@docs
grovebit
```

### Transformations

```@docs
Dendriform.TreeBase
```

```@docs
Dendriform.TreeInteger
```

```@docs
Dendriform.TreeRational
```

## Tools & Options

```@docs
grovesort
```

```@docs
treeshift
```
