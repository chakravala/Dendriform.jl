__precompile__()
module Dendriform

# This file is part of Dendriform.jl. It is licensed under the GPL license
# Dendriform Copyright (C) 2017 Michael Reed

using Combinatorics, AbstractLattices

export PBTree, Grove, GroveBin, ==, Cn, grovesort, grovesort!, σ, print, grovecomposition, grovedisplay

import Base: ==, +,*, /, ∪, <, >, ≤, ≥, promote_rule, convert, show, print

# definitions

abstract type AbstractGrove end

"""
Planar Binary Tree with Loday's notation

### Summary
mutable struct PBTree <: AbstractGrove

### Fields
degr::UInt8
Y::Array{UInt8,1}
"""
mutable struct PBTree <: AbstractGrove
    degr::UInt8
    Y::Array{UInt8,1}
end

"""
Grove of planar binary trees, matrix equivalence class

### Summary
mutable struct Grove <: AbstractGrove

### Fields
degr::UInt8
size::Int
Y::Array{UInt8,2}
"""
mutable struct Grove <: AbstractGrove
    degr::UInt8
    size::Int
    Y::Array{UInt8,2}
end

"""
Compressed binary representation of grove

### Summary
mutable struct GroveBin <: AbstractGrove

### Fields
degr::UInt8
size::Int
gbin::Integer
ppos::Float16
"""
mutable struct GroveBin <: AbstractGrove
    degr::UInt8
    size::Int
    gbin::Integer
    ppos::Float16
end

"""
Descending greatest integer search data for grove

### Summary
mutable struct BaseTree <: Abstract Grove

### Fields
μ::Array{Array{UInt8,1},1}
"""
mutable struct BaseTree <: AbstractGrove
    μ::Array{Array{UInt8,1},1}
end

Ar1UI8I=Union{Array{UInt8,1},Array{Int,1}}
Ar2UI8I=Union{Array{UInt8,2},Array{Int,2}}
AbstractPBTree = Union{PBTree,Ar1UI8I}
UI8I = Union{UInt8,Int}
NotGrove = Union{GroveBin,AbstractPBTree,Ar2UI8I,UI8I}
PureGrove = Union{Grove,GroveBin,Ar2UI8I}
Cn = catalannum

# constructors

"""
    PBTree(degree::Int, index::Int)

Looks up  PBTree using degree and tree index
"""
function PBTree(deg::UI8I,ind::Int)
    treecheck(deg,ind)
    deg == 0 && (return PBTree(UInt8(deg),Array{UInt8,1}(undef,0)))
    return PBTree(UInt8(deg),Υ(deg).Y[ind,:])
end

"""
    PBTree(t::Array{Int,1})

Creates PBTree using Loday label array
"""
PBTree(t::Ar1UI8I) = convert(PBTree,t)

"""
    Grove(g::Array{Int,n})

Creates Grove using Loday label array of dimension 1 or 2
"""
Grove(g::Ar2UI8I) = convert(Grove,g)
Grove(t::Ar1UI8I) = convert(Grove,t)
Grove(d::UI8I,g::Ar2UI8I) = Grove(UInt8(d),GroveSiz(g),convert(Array{UInt8,2},g))

"""
    Grove(g::AbstractGrove)

Creates Grove object from any NotGrove object
"""
Grove(t::PBTree) = convert(Grove,t)
Grove(d::UI8I) = convert(Grove,d)
Grove(d::UI8I,s::BitArray{1}) = TreeLoday(d,s)
Grove(g::GroveBin) = convert(Grove,g)
Grove(s::BitArray{1}) = Grove(CnInv(length(s)),s)
Grove(g::Grove) = g

"""
    Grove(degree::Int, index::Integer)

Looks up Grove using degree and grove index
"""
Grove(d::UI8I,s::Integer) = Grove(d,grovebit(d,s))

"""
    GroveBin(g::AbstractGrove)

Converts Grove g into a GroveBin
"""
GroveBin(g::Grove) = GroveBin(UInt8(g.degr),g.size,groveindex(g))
GroveBin(g::NotGrove) = GroveBin(convert(Grove,g))
GroveBin(d::UI8I,s::Int,i::Integer) = GroveBin(UInt8(d),s,i,Float16(100i//(2^Cn(d)-1)))

==(a::PBTree,b::PBTree)=(a.degr == b.degr && a.Y == b.Y)
==(a::Grove,b::Grove)=(a.degr==b.degr && a.size==b.size && grovesort!(a).Y==grovesort!(b).Y)
==(a::BaseTree,b::BaseTree)=(a.μ==b.μ)
==(a::GroveBin,b::GroveBin)=(a.degr==b.degr && a.size==b.size && a.gbin==b.gbin)

# conversions / promotions

function convert(::Type{PBTree},t::Ar1UI8I)
    return PBTree(isempty(t) ? 0 : UInt8(length(t)),convert(Array{UInt8,1},t))
end

function convert(::Type{Grove},t::PBTree)
    d=t.degr
    return Grove(d,1,(g=Array{Int,2}(undef,1,d); g[1,:]=t.Y[:]; g))
end

function convert(::Type{Grove},g::Ar2UI8I)
    return Grove(GroveDeg(g),GroveSiz(g),convert(Array{UInt8,2},g))
end

GroveDeg(g::Ar2UI8I) = isempty(g) ? 0 : UInt8(length(g[1,:]))
GroveSiz(g::Ar2UI8I) = isempty(g) ? 0 : length(g[:,1])
convert(::Type{Grove},t::Ar1UI8I) = Grove(PBTree(convert(Array{UInt8,1},t)))
convert(::Type{Grove},g::Array{Any,1}) = convert(Grove,convert(Array{UInt8,1},g))
convert(::Type{Grove},g::Array{Any,2}) = convert(Grove,convert(Array{UInt8,2},g))
convert(::Type{Grove},g::GroveBin) = Grove(g.degr,g.gbin)
convert(::Type{Grove},d::UI8I) = Υ(UInt8(d))
promote_rule(::Type{PBTree},::Type{T}) where T<:Ar1UI8I = PBTree
promote_rule(::Type{Grove},::Type{T}) where T<:Union{Ar1UI8I,Ar2UI8I,PBTree,UI8I} = Grove

# display

show(io::IO, ::MIME"text/plain", t::PBTree) = print(io,t)
show(io::IO, ::MIME"text/plain", g::Grove) = print(io,g)
show(io::IO, ::MIME"text/plain", k::GroveBin) = print(io,k)

# Sorting

grovesort!(g::Grove) = grovesort!(g.Y,TreeInteger(g))
grovesort!(g::NotGrove) = grovesort!(convert(Grove,g))
grovesort!(g::Ar2UI8I,Θ::Array{Int,1}) = grovesort!(Grove(g),Θ)
grovesort!(g::Grove,Θ::Array{Int,1}) = (g.Y[:,:] = g.Y[sortperm(Θ),:]; g)

"""
    grovesort(::Bool)

Toggles the grovesort algorithm (enabled by default - RECOMMENDED)
"""
grovesort = ( () -> begin
        gs=true
        return (tf=gs)->(gs≠tf && (gs=tf; ΥGS()); return gs)
    end)()

# Inverse catalannum & Involution

function CnInv(n::UI8I)
    d = 1
    k = Cn(d)
    while k < n
        d +=1
        k = Cn(d)
    end
    k == n ? (return d) : error("$n is not a Catalan number")
end

"""
    σ(g::AbstractGrove)

Applies the involution to any PBTree or Grove object
"""
σ(x::Union{Grove,GroveBin}) = σ(Grove(x).Y)
σ(x::Ar2UI8I) = Grove(x[:,end:-1:1])
σ(x::PBTree) = σ(x.Y)
σ(x::Ar1UI8I) = PBTree(x[end:-1:1])

include("morphism.jl")
include("arithmetic.jl")
include("poset.jl")

# Total Grove Repository

function GroveExtend() # initialize set of total groves
    Y = Array{Grove,1}(undef,1)
    Y[1]=Grove(0,0,Array{UInt8,2}(undef,0,1))
    R = Array{Array{Int,1},1}(undef,0)
    return (Y,R)
end

function GroveExtend!(Y::Array{Grove,1},R::Array{Array{Int,1},1},deg::UInt8)
    D = UInt8(length(Y)):deg
    !isempty(D) && println("Extend Grove Degree Level,")
    for n ∈ D # loop over all degree levels
        print(" $n")
        τ = 1
        cn = ceil(Int,n/2)
        fn = floor(Int,n/2)
        lYn = length(Y[n].Y[:,1])
        # initialize total grove
        Yn = Grove(Array{UInt8,2}(undef,Cn(n),n))
        # loop over left-branch grove
        n==1 ? (Yn.Y[τ,1] = n) : (Yn.Y[τ:lYn,n] .= n)
        for λ ∈ 1:lYn
            Yn.Y[τ,1:n-1] = Y[n].Y[λ,:]
            τ += 1
        end
        # loop over (right) root indices
        for ν ∈ n-1:-1:cn+1
            τ=GroveBuild!(Y[ν].Y,Y[n-ν+1].Y,Yn.Y,ν,n,τ,1,Y[ν].size,1,Y[n-ν+1].size)
        end
        # loop over center root next (if n is odd)
        cn ≠ fn && for ν = cn
            λs = length(Y[ν].Y[:,1])
            Λs = length(Y[n-ν+1].Y[:,1])
            fλs = floor(Int,λs/2)
            cλs = ceil(Int,λs/2)
            fΛs = floor(Int,Λs/2)
            cΛs = ceil(Int,Λs/2)
            τ = GroveBuild!(Y[ν].Y,Y[n-ν+1].Y,Yn.Y,ν,n,τ,1,fλs,1,fΛs)
            cΛs ≠ fΛs && (τ = GroveBuild!(Y[ν].Y,Y[n-ν+1].Y,Yn.Y,ν,n,τ,1,fλs,cΛs,cΛs))
            τ = GroveBuild!(Y[ν].Y,Y[n-ν+1].Y,Yn.Y,ν,n,τ,1,fλs,cΛs+1,Λs)
            cλs ≠ fλs && (τ = GroveBuild!(Y[ν].Y,Y[n-ν+1].Y,Yn.Y,ν,n,τ,cλs,cλs,1,Λs))
            τ = GroveBuild!(Y[ν].Y,Y[n-ν+1].Y,Yn.Y,ν,n,τ,cλs+1,λs,1,fΛs)
            cΛs ≠ fΛs && (τ = GroveBuild!(Y[ν].Y,Y[n-ν+1].Y,Yn.Y,ν,n,τ,cλs+1,λs,cΛs,cΛs))
            τ = GroveBuild!(Y[ν].Y,Y[n-ν+1].Y,Yn.Y,ν,n,τ,cλs+1,λs,cΛs+1,Λs)
        end
        for ν ∈ fn:-1:2 # loop over (right) root indices
            for Λ ∈ 1:Y[n-ν+1].size # loop over right grove
                for λ ∈ 1:Y[ν].size # loop left grove
                    Yn.Y[τ,:] = (Y[ν].Y[λ,:] ∨ Y[n-ν+1].Y[Λ,:]).Y
                    τ += 1
        end; end; end
        # loop over right-branch grove
        Yn.Y[τ:τ+lYn-1,1] .= n
        for Λ ∈ 1:lYn
            Yn.Y[τ,2:n] = Y[n].Y[Λ,:]
            τ += 1
        end
        grovesort() && (print("|Θ"); push!(R,TreeInteger(TreeBase(Yn))); grovesort!(Yn,R[n]); sort!(R[n]))
        push!(Y,Yn)
    end
    !isempty(D) && print("\n")
    return Y[deg+1]
end

function GroveBuild!(Yν::Array{UInt8,2},Ynν1::Array{UInt8,2},Yn::Array{UInt8,2},ν::UI8I,n::UInt8,τ::Int,λs1::Int,λs2::Int,Λs1::Int,Λs2::Int)
    #loop left/right grove, graft
    for λ ∈ λs1:λs2
        for Λ ∈ Λs1:Λs2
            Yn[τ,:] = (Yν[λ,:] ∨ Ynν1[Λ,:]).Y
            τ+=1
    end; end
    return τ
end

function GroveInteger!(Y::Array{Grove,1},R::Array{Array{Int,1},1},deg::UInt8)
    D = UInt8(length(R)+1):deg
    !isempty(D) && println("Extend Grove Integer Level,")
    for n ∈ D
        print(" $n")
        push!(R,TreeInteger(TreeBase(Y[n+1])))
    end
    !isempty(D) && print("\n")
    return R[deg]
end

(Υ,ΥI,ΥGS) = ( () -> begin
        (Y,R)=GroveExtend()
        return begin
                (d::UI8I)->(return GroveExtend!(Y,R,UInt8(d))),
                (d::UI8I)->(GroveExtend!(Y,R,UInt8(d)); !grovesort() && GroveInteger!(Y,R,UInt8(d)); return d==0 ? Array{Int,1}(undef,0) : R[d]),
                ()->((Y,R)=GroveExtend())
            end
    end)() # provides hidden total grove reference

# Grove Composition

function Compose(n::Int,η::Int=n)
    G = GroveSums(n)
    u=1
    !isempty(G) && (return G)
    n < η && (u = 2^Cn(n))
    n ≠ 0 && n < η && for i∈1:u-1
        push!(G,[GroveBin(Grove(n,i)),GroveBin(Grove(n,i))])
    end
    for s ∈ n-1:-1:1
        for i ∈ 1:2^Cn(s)-1
            g = Compose(n-s,η)
            gsi = Grove(s,i)
            gbsi = GroveBin(gsi)
            for r ∈ 1:length(g)
                sm = gsi + Grove(g[r][end].degr,g[r][end].gbin)
                push!(G,g[r])
                pushfirst!(G[u],gbsi)
                G[u][end] = GroveBin(sm)
                u += 1
    end; end; end
    c = Dict{Integer,Array{Integer,1}}()
    for k ∈ 1:length(G)
        ind = G[k][end].gbin
        try
            c[ind]
        catch
            push!(c,ind=>Array{Integer,1}(undef,0))
        end
        push!(c[ind],k)
    end
  GroveStore(n,G,c)
  return G
end

(GroveStore,GroveComp,GroveSums) = ( () -> begin
        GG=Array{Array{Array{GroveBin,1},1},1}(undef,0)
        CC=Array{Dict{Integer,Array{Integer,1}},1}(undef,0)
        return begin
                (n::Int,g::Array{Array{GroveBin,1},1},c::Dict{Integer,Array{Integer,1}})->(n==length(GG)+1 && (push!(GG,deepcopy(g)); push!(CC,deepcopy(c)))),
                (n::Int)->(Compose(n); (return deepcopy(CC[n]))),
                (n::Int)->(n<=length(GG) && n>0x00 ? (return deepcopy(GG[n])) : (return Array{Array{GroveBin,1},1}(undef,0) )  )
            end
    end)()

function grovecomposition(d::UI8I,ind::Integer)
    print(GroveBin(Grove(d,ind)))
    gi = GroveComp(d)
    com = Compose(d)
    try
        gi[ind]
    catch
        print(" has 1 composition (itself)\n")
        return 1
    end
    lg = length(gi[ind])
    print(" has $(lg+1) compositions\n")
    for k ∈ 1:lg
        print("(")
        print(com[gi[ind][k]][1])
        for t ∈ 2:length(com[gi[ind][k]])-1
            print(") + (")
            print(com[gi[ind][k]][t])
        end
        print(")\n")
    end
    return lg+1
end

# printing

"""
    grovedisplay(::Bool)

Toggles the display output of grove index data (disabled by default)
"""
grovedisplay = ( () -> begin
        gs=false
        return (tf=gs)->(gs≠tf && (gs=tf); return gs)
    end)()



function print(io::IO,υ::PBTree,μ::BaseTree)
    n = υ.degr
    if n == 0x00
        print(io,"∅"*(grovedisplay() ? " ↦ [∅] ↦ 0/1 or 0" : ""))
        return nothing
    end
    ti = TreeInteger(μ)
    tin = treeindex(n,ti)
    show(io,convert(Array{Int,1},υ.Y))
    if grovedisplay()
        print(io," ↦ ")
        for ω ∈ 1:length(υ.Y)
            μ.μ[ω]==[] ? print(io,'∅') : show(io,convert(Array{Int,1},μ.μ[ω]))
        end
        print(io," ↦ ")
        print(io,tin,"/",Cn(n)," or ")
        show(io,ti)
    end
    print(io,'\n')
end

print(io::IO,υ::PBTree) = print(io,υ,TreeBase(υ));
print(io::IO,μ::BaseTree) = print(io,TreeLoday(μ),μ)

function print(io::IO,Y::Grove) # given Loday label grove
    for η ∈ 1:Y.size
        print(io,PBTree(Y.Y[η,:]))
    end
    if grovedisplay()
        print(io,GroveBin(Y))
    else
        print(io,"Y$(Y.degr) "*'#'*"$(Y.size)/$(Cn(Y.degr))")
    end
end

function print(io::IO,Y::Array{BaseTree,1}) # given Index label grove
    for η ∈ 1:length(Y)
        print(io,Y[η])
    end
end

function print(io::IO,Y::Array{Grove,1})
    for n ∈ 1:length(Y)
        print(io,Y[n])
    end
end

function print(io::IO,k::GroveBin) 
    print(io,"$(k.gbin) Y$(k.degr) "*'#'*"$(k.size)/$(Cn(k.degr)) [$(k.ppos)"*'%'*"]")
end

end
