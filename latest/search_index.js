var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#[Dendriform.jl](https://github.com/chakravala/Dendriform.jl)-1",
    "page": "Home",
    "title": "Dendriform.jl",
    "category": "section",
    "text": "Symbolic parser generator for Julia language expressions using REDUCE algebra term rewrite system"
},

{
    "location": "index.html#Setup-1",
    "page": "Home",
    "title": "Setup",
    "category": "section",
    "text": "Installation of latest release version using Julia:julia> Pkg.add(\"Dendriform\")Provides the types PBTree for planar binary trees, Grove for tree collections of constant degree, and GroveBin to compress grove data. This package defines various essential operations on planar binary trees and groves like ∪ for union; ∨ for graft; left and right for branching; <, >, ≤, ≥ for Tamari's partial ordering; ↗ and ↖ (i.e. over and under); and the dashv and vdash operations ⊣, ⊢, +, * for dendriform algebra."
},

{
    "location": "index.html#Background-1",
    "page": "Home",
    "title": "Background",
    "category": "section",
    "text": "We call (Image: tree-symb) the name of a tree to represent it as a vector, where the sequence is made up of n integers. Collections of planar binary trees are encoded into an equivalence class of matrices:(Image: matrix-equivalence-class)where (Image: A~B) if there exists a permutation (Image: f  in Sk) so that (Image: condition). The binary tree grafting operation is computed<p align=\"center\"><img src=\"https://latex.codecogs.com/svg.latex?\\omega(\\alpha\\vee&space;\\beta)&space;=&space;\\omega(\\alpha)\\vee\\omega(\\beta)&space;:=&space;[\\omega(\\alpha),a&plus;1&plus;b,\\omega(\\beta)]\\in&space;\\Lambda_{a&plus;b&plus;1}^1\"/></p>The left and right addition are computed on the following recursive principle:<p align=\"center\"><img src=\"https://latex.codecogs.com/svg.latex?\\xi\\dashv&space;\\eta&space;&=&space;\\bigcup_{i}&space;\\bigcup_{\\tau&space;\\in&space;\\xi_i^r&space;&plus;&space;\\eta}&space;\\xi_i^l&space;\\vee&space;\\tau&space;\\qquad&space;&\\text{and}&space;\\qquad&space;\\qquad&space;\\xi\\vdash&space;\\eta&space;&=&space;\\bigcup_{j}&space;\\bigcup_{\\tau&space;\\in&space;\\xi&plus;\\eta_j^l}&space;\\tau\\vee&space;\\eta_j^r.\"/></p>Together these non-commutative binary operations satisfy the properties of an associative dendriform dialgebra. The structures induced by Loday's definition of the sum have the partial ordering of the associahedron known as Tamari lattice.<p align=\"center\">   <img src=\"https://raw.githubusercontent.com/wiki/chakravala/Fatou.jl/dendriform/grove-sum-1.png\" alt=\"tamari-grove-commutativity.png\"/> </p>Figure: Tamari associahedron, colored to visualize noncommutative sums of [1,2] and [2,1], code: gistHowever, in this computational package, a stricter total ordering is constructed using a function that transforms the set-vector isomorphism obtained from the descending greatest integer index position search method:<p align=\"center\"><img src=\"https://latex.codecogs.com/svg.latex?\\Theta(\\mu)&space;&=&space;\\sum_{j=n}^1&space;\\sum_{k=1}^{\\&hash;e_j}&space;(e_j)_k&space;\\cdot&space;10^{\\delta(j,k)},&space;\\qquad&space;&\\text{where}&space;\\qquad&space;\\delta(j,k)&space;&=&space;n&space;-&space;\\sum_{r=1}^{j-1}&space;\\sum_{s=1}^{\\&hash;e_r}&space;1&space;-&space;\\sum_{s=1}^{k}&space;1\"/></p>The structure obtained from this total ordering is used to construct a reliable binary groveindex representation that encodes the essential data of any grove, using the formula<p align=\"center\"><img src=\"https://latex.codecogs.com/svg.latex?\\zeta_\\gamma&space;:=&space;\\sum_{\\tau&space;\\in&space;\\gamma}&space;2^{\\theta_\\tau&space;-&space;1}\"/></p>These algorithms are used in order to facilitate computations that provide insight into the Loday arithmetic."
},

{
    "location": "index.html#Usage-1",
    "page": "Home",
    "title": "Usage",
    "category": "section",
    "text": "Basic usage examples:julia> using Dendriform\n\njulia> Grove(3,7) ⊣ [1,2]∪[2,1]\n[1,2,5,1,2] ↦ [3]∅∅[2,5][1,4] ↦ 20/42 or 21807\n[1,2,5,2,1] ↦ [3]∅∅[2,4][1,5] ↦ 21/42 or 21906\n[2,1,5,1,2] ↦ [3]∅∅[1,5][2,4] ↦ 22/42 or 22797\n[2,1,5,2,1] ↦ [3]∅∅[1,4][2,5] ↦ 23/42 or 22896\n[1,5,3,1,2] ↦ [2]∅[3][5][1,4] ↦ 27/42 or 30807\n[1,5,2,1,3] ↦ [2]∅[5][3][1,4] ↦ 25/42 or 29007\n[1,5,1,2,3] ↦ [2]∅[5][4][1,3] ↦ 24/42 or 28908\n[1,5,3,2,1] ↦ [2]∅[3][4][1,5] ↦ 28/42 or 30906\n[1,5,1,3,1] ↦ [2]∅[4]∅[1,3,5] ↦ 26/42 or 30186\n267911168 Y5 #9/42 [0.006092%]\n\njulia> Grove(2,3) * [1,2,3]∪[3,2,1] |> GroveBin\n2981131286847743360614880957207748817969 Y6 #30/132 [54.75%]\n\njulia> [2,1,7,4,1,3,1] < [2,1,7,4,3,2,1]\ntrue"
},

{
    "location": "index.html#References-1",
    "page": "Home",
    "title": "References",
    "category": "section",
    "text": "Dan Yasaki with Adriano Bruno, The arithmetic of planar binary trees, Involve 4 (2011), no. 1, 1-11. (PDF)\nJean-Louis Loday, Arithmetree, J. of Algebra (2002), no. 258, 275-309."
},

{
    "location": "library.html#",
    "page": "Library",
    "title": "Library",
    "category": "page",
    "text": ""
},

{
    "location": "library.html#Dendriform.jl-Library-1",
    "page": "Library",
    "title": "Dendriform.jl Library",
    "category": "section",
    "text": ""
},

{
    "location": "library.html#Index-1",
    "page": "Library",
    "title": "Index",
    "category": "section",
    "text": ""
},

{
    "location": "library.html#Dendriform.PBTree",
    "page": "Library",
    "title": "Dendriform.PBTree",
    "category": "Type",
    "text": "Planar Binary Tree with Loday's notation\n\nSummary\n\nmutable struct PBTree <: AbstractGrove\n\nFields\n\ndegr::UInt8 Y::Array{UInt8,1}\n\n\n\n"
},

{
    "location": "library.html#Dendriform.Grove",
    "page": "Library",
    "title": "Dendriform.Grove",
    "category": "Type",
    "text": "Grove of planar binary trees, matrix equivalence class\n\nSummary\n\nmutable struct Grove <: AbstractGrove\n\nFields\n\ndegr::UInt8 size::Int Y::Array{UInt8,2}\n\n\n\n"
},

{
    "location": "library.html#Dendriform.GroveBin",
    "page": "Library",
    "title": "Dendriform.GroveBin",
    "category": "Type",
    "text": "Compressed binary representation of grove\n\nSummary\n\nmutable struct GroveBin <: AbstractGrove\n\nFields\n\ndegr::UInt8 size::Int gbin::Integer ppos::Float16\n\n\n\n"
},

{
    "location": "library.html#Dendriform.BaseTree",
    "page": "Library",
    "title": "Dendriform.BaseTree",
    "category": "Type",
    "text": "Descending greatest integer search data for grove\n\nSummary\n\nmutable struct BaseTree <: Abstract Grove\n\nFields\n\nμ::Array{Array{UInt8,1},1}\n\n\n\n"
},

{
    "location": "library.html#Dendriform-Types-1",
    "page": "Library",
    "title": "Dendriform Types",
    "category": "section",
    "text": "PBTreeGroveGroveBinDendriform.BaseTree"
},

{
    "location": "library.html#Dendriform-Operators-1",
    "page": "Library",
    "title": "Dendriform Operators",
    "category": "section",
    "text": "Dendriform.∪"
},

{
    "location": "library.html#Dendriform.graft",
    "page": "Library",
    "title": "Dendriform.graft",
    "category": "Function",
    "text": "graft(left::AbstractPBTree, right::AbstractPBTree)\n\nGrafts the left and right PBTree with root vertex\n\n\n\n"
},

{
    "location": "library.html#Dendriform.:∨",
    "page": "Library",
    "title": "Dendriform.:∨",
    "category": "Function",
    "text": "∨(left::AbstractPBTree, right::AbstractPBTree)\n\nGrafts the left and right AbstractPBTree objects\n\n\n\n"
},

{
    "location": "library.html#Dendriform.left",
    "page": "Library",
    "title": "Dendriform.left",
    "category": "Function",
    "text": "left(::AbstractPBTree)\n\nReturns the left branch of an AbstractPBTree\n\n\n\n"
},

{
    "location": "library.html#Dendriform.right",
    "page": "Library",
    "title": "Dendriform.right",
    "category": "Function",
    "text": "right(::AbstractPBTree)\n\nReturns the right branch of an AbstractPBTree\n\n\n\n"
},

{
    "location": "library.html#Branching-1",
    "page": "Library",
    "title": "Branching",
    "category": "section",
    "text": "graft∨leftright"
},

{
    "location": "library.html#Dendriform.posetnext",
    "page": "Library",
    "title": "Dendriform.posetnext",
    "category": "Function",
    "text": "Dendriform.posetnex(::PBTree)\n\nReturns an Array{PBTree,1} of trees that are greater than it\n\n\n\n"
},

{
    "location": "library.html#Dendriform.posetprev",
    "page": "Library",
    "title": "Dendriform.posetprev",
    "category": "Function",
    "text": "Dendriform.posetprev(::PBTree)\n\nReturns an Array{PBTree,1} of trees that are less than it\n\n\n\n"
},

{
    "location": "library.html#Partial-Ordering-1",
    "page": "Library",
    "title": "Partial Ordering",
    "category": "section",
    "text": "Dendriform.posetnextDendriform.posetprevDendriform.<Dendriform.>Dendriform.≤Dendriform.≥"
},

{
    "location": "library.html#Dendriform.over",
    "page": "Library",
    "title": "Dendriform.over",
    "category": "Function",
    "text": "over(a::AbstractPBTree, b::AbstractPBTree)\n\nReturns PBTRee obtained from a over b operation\n\n\n\n"
},

{
    "location": "library.html#Dendriform.under",
    "page": "Library",
    "title": "Dendriform.under",
    "category": "Function",
    "text": "under(a::AbstractPBTree, b::AbstractPBTree)\n\nReturns PBTRee obtained from a under b operation\n\n\n\n"
},

{
    "location": "library.html#Dendriform.↗",
    "page": "Library",
    "title": "Dendriform.↗",
    "category": "Function",
    "text": "↗(a::AbstractPBTree, b::AbstractPBTree)\n\nReturns PBTRee obtained from a over b operation\n\n\n\n"
},

{
    "location": "library.html#Dendriform.↖",
    "page": "Library",
    "title": "Dendriform.↖",
    "category": "Function",
    "text": "↖(a::AbstractPBTree, b::AbstractPBTree)\n\nReturns PBTRee obtained from a under b operation\n\n\n\n"
},

{
    "location": "library.html#Dendriform.dashv",
    "page": "Library",
    "title": "Dendriform.dashv",
    "category": "Function",
    "text": "dashv(a::AbstractGrove, b::AbstractGrove)\n\nReturns Grove obtained from a ⊣ b operation\n\n\n\n"
},

{
    "location": "library.html#Dendriform.vdash",
    "page": "Library",
    "title": "Dendriform.vdash",
    "category": "Function",
    "text": "vdash(a::AbstractGrove, b::AbstractGrove)\n\nReturns Grove obtained from a ⊢ b operation\n\n\n\n"
},

{
    "location": "library.html#Dendriform.:⊣",
    "page": "Library",
    "title": "Dendriform.:⊣",
    "category": "Function",
    "text": "⊣(a::AbstractGrove, b::AbstractGrove)\n\nReturns Grove obtained from a ⊣ b operation\n\n\n\n"
},

{
    "location": "library.html#Dendriform.:⊢",
    "page": "Library",
    "title": "Dendriform.:⊢",
    "category": "Function",
    "text": "⊢(a::AbstractGrove, b::AbstractGrove)\n\nReturns Grove obtained from a ⊢ b operation\n\n\n\n"
},

{
    "location": "library.html#Dialgebra-Arithmetic-1",
    "page": "Library",
    "title": "Dialgebra Arithmetic",
    "category": "section",
    "text": "overunder↗↖dashvvdash⊣⊢Dendriform.+Dendriform.*"
},

{
    "location": "library.html#Dendriform.σ",
    "page": "Library",
    "title": "Dendriform.σ",
    "category": "Function",
    "text": "σ(g::AbstractGrove)\n\nApplies the involution to any PBTree or Grove object\n\n\n\n"
},

{
    "location": "library.html#Dendriform-Morphisms-1",
    "page": "Library",
    "title": "Dendriform Morphisms",
    "category": "section",
    "text": "σ"
},

{
    "location": "library.html#Dendriform.LeftInherited",
    "page": "Library",
    "title": "Dendriform.LeftInherited",
    "category": "Function",
    "text": "Dendriform.LeftInherited(::AbstractPBTree)\n\nReturns Bool that tells if PBTree is left inherited\n\n\n\n"
},

{
    "location": "library.html#Dendriform.RightInherited",
    "page": "Library",
    "title": "Dendriform.RightInherited",
    "category": "Function",
    "text": "Dendriform.RightInherited(::AbstractPBTree)\n\nReturns Bool that tells if PBTree is right inherited\n\n\n\n"
},

{
    "location": "library.html#Dendriform.PrimitiveTree",
    "page": "Library",
    "title": "Dendriform.PrimitiveTree",
    "category": "Function",
    "text": "Dendriform.PrimitiveTree(::AbstractPBTree)\n\nReturns Bool that tells if PBTree is primitive\n\n\n\n"
},

{
    "location": "library.html#Dendriform.treecheck",
    "page": "Library",
    "title": "Dendriform.treecheck",
    "category": "Function",
    "text": "treecheck(::AbstractPBTree)\n\nReturns Bool that tells if PBTree is valid\n\n\n\n"
},

{
    "location": "library.html#Dendriform.grovecheck",
    "page": "Library",
    "title": "Dendriform.grovecheck",
    "category": "Function",
    "text": "grovecheck(::AbstractGrove)\n\nReturns Bool that tells if Grove is valid\n\n\n\n"
},

{
    "location": "library.html#Dendriform.GroveError",
    "page": "Library",
    "title": "Dendriform.GroveError",
    "category": "Function",
    "text": "Dendriform.GroveError(::AbstractGrove)\n\nReturns Array with Grove sorting index error\n\n\n\n"
},

{
    "location": "library.html#Consistency-Checks-1",
    "page": "Library",
    "title": "Consistency Checks",
    "category": "section",
    "text": "Dendriform.LeftInheritedDendriform.RightInheritedDendriform.PrimitiveTreetreecheckgrovecheckDendriform.GroveError"
},

{
    "location": "library.html#Dendriform.treeindex",
    "page": "Library",
    "title": "Dendriform.treeindex",
    "category": "Function",
    "text": "treeindex(::AbstractGrove)\n\nReturns tree indices of any PBTree or Grove\n\n\n\n"
},

{
    "location": "library.html#Dendriform.groveindex",
    "page": "Library",
    "title": "Dendriform.groveindex",
    "category": "Function",
    "text": "groveindex(::AbstractGrove)\n\nReturns the grove index of any Grove\n\n\n\n"
},

{
    "location": "library.html#Dendriform.grovebit",
    "page": "Library",
    "title": "Dendriform.grovebit",
    "category": "Function",
    "text": "grovebit(::AbstractGrove)\n\nReturns a BitArray of tree indices\n\n\n\n"
},

{
    "location": "library.html#Index-Data-1",
    "page": "Library",
    "title": "Index Data",
    "category": "section",
    "text": "treeindextreeindexCngroveindexgrovebit"
},

{
    "location": "library.html#Dendriform.TreeBase",
    "page": "Library",
    "title": "Dendriform.TreeBase",
    "category": "Function",
    "text": "Dendriform.TreeBase(::AbstractGrove)\n\nReturns BaseTree objects for any AbstractGrove\n\n\n\n"
},

{
    "location": "library.html#Dendriform.TreeInteger",
    "page": "Library",
    "title": "Dendriform.TreeInteger",
    "category": "Function",
    "text": "Dendriform.TreeInteger(::AbstractGrove)\n\nReturns the tree integers of any AbstractGrove\n\n\n\n"
},

{
    "location": "library.html#Dendriform.TreeRational",
    "page": "Library",
    "title": "Dendriform.TreeRational",
    "category": "Function",
    "text": "Dendriform.TreeRational(::AbstractGrove)\n\nReturns the tree rationals of any AbstractGrove\n\n\n\n"
},

{
    "location": "library.html#Transformations-1",
    "page": "Library",
    "title": "Transformations",
    "category": "section",
    "text": "Dendriform.TreeBaseDendriform.TreeIntegerDendriform.TreeRational"
},

{
    "location": "library.html#Dendriform.grovesort",
    "page": "Library",
    "title": "Dendriform.grovesort",
    "category": "Function",
    "text": "grovesort(::Bool)\n\nToggles the grovesort algorithm (enabled by default - RECOMMENDED)\n\n\n\n"
},

{
    "location": "library.html#Dendriform.treeshift",
    "page": "Library",
    "title": "Dendriform.treeshift",
    "category": "Function",
    "text": "treeshift(::Bool)\n\nToggles the shift for the tree integers / rationals\n\n\n\n"
},

{
    "location": "library.html#Tools-and-Options-1",
    "page": "Library",
    "title": "Tools & Options",
    "category": "section",
    "text": "grovesorttreeshift"
},

]}
