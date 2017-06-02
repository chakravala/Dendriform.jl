__precompile__()
module GroveAlg

# This file is part of GroveAlg.jl. It is licensed under the GPL license
# GroveAlg Copyright (C) 2017 Michael Reed

export PBTree, Grove, GroveBin, ==, GroveSort, GroveSort!, BranchLeft, BranchRight, Cn, Graft, GrovePrint

abstract AbstractGrove; importall Base; using Combinatorics
type PBTree <: AbstractGrove; degr::UInt8; Y::Array{UInt8,1}; end
type Grove <: AbstractGrove; degr::UInt8; size::Int; Y::Array{UInt8,2}; end
type GroveBin <: AbstractGrove; degr::UInt8; size::Int; gbin::Integer; ppos::Float16; end
type BaseTree <: AbstractGrove; μ::Array{Array{UInt8,1},1}; end
Ar1UI8I=Union{Array{UInt8,1},Array{Int,1}}; Ar2UI8I=Union{Array{UInt8,2},Array{Int,2}}
AbstractPBTree = Union{PBTree,Ar1UI8I}; UI8I = Union{UInt8,Int}; Cn = catalannum
NotGrove = Union{GroveBin,AbstractPBTree,Ar2UI8I,UI8I}

PBTree(t::Ar1UI8I) = convert(PBTree,t)
PBTree(deg::UI8I,ind::Int) = PBTree(UInt8(deg),Υ(deg).Y[ind,:])
Grove(t::Ar1UI8I) = convert(Grove,t); Grove(g::Ar2UI8I) = convert(Grove,g)
Grove(d::UI8I,g::Ar2UI8I) = Grove(UInt8(d),GroveSiz(g),convert(Array{UInt8,2},g))
Grove(t::PBTree) = convert(Grove,t); Grove(d::UI8I) = convert(Grove,d)
Grove(d::UI8I,s::BitArray{1}) = TreeLoday(d,s); Grove(g::GroveBin) = convert(Grove,g)
Grove(d::UI8I,s::Integer) = Grove(d,GroveBit(s))
GroveBin(g::Grove) = GroveBin(UInt8(g.degr),g.size,GroveIndex(g))
GroveBin(g::NotGrove) = GroveBin(convert(Grove,g))
GroveBin(d::UI8I,s::Int,i::Integer) = GroveBin(UInt8(d),s,i,Float16(100i//(2^Cn(d)-1)))
==(a::PBTree,b::PBTree)=(a.degr == b.degr && a.Y == b.Y)
==(a::Grove,b::Grove)=(a.degr==b.degr && a.size==b.size && GroveSort!(a).Y==GroveSort!(b).Y)
==(a::BaseTree,b::BaseTree)=(a.μ==b.μ)

# Conversions / Promotions

GroveDeg(g::Ar2UI8I) = isempty(g) ? 0 : UInt8(length(g[1,:]))
GroveSiz(g::Ar2UI8I) = isempty(g) ? 0 : length(g[:,1])
function convert(::Type{PBTree},t::Ar1UI8I)
  return PBTree(isempty(t) ? 0 : UInt8(length(t)),convert(Array{UInt8,1},t)); end
function convert(::Type{Grove},t::PBTree)
  d=t.degr; return Grove(d,1,(g=Array{Int,2}(1,d); g[1,:]=t.Y[:]; g)); end
function convert(::Type{Grove},g::Ar2UI8I)
  return Grove(GroveDeg(g),GroveSiz(g),convert(Array{UInt8,2},g)); end
convert(::Type{Grove},t::Ar1UI8I) = Grove(PBTree(convert(Array{UInt8,1},t)))
convert(::Type{Grove},g::GroveBin) = Grove(g.degr,g.gbin)
convert(::Type{Grove},d::UI8I) = Υ(UInt8(d))
#promote_rule(::Type{PBTree},::Type{Ar1UI8I}) = PBTree
#promote_rule(::Type{Grove},::Type{Union{Ar1UI8I,Ar2UI8I,PBTree,UI8I}})=Grove
show(io::IO,k::GroveBin) = print(io,"$(k.gbin) Y$(k.degr) \#$(k.size)/$(Cn(k.degr)) [$(k.ppos)\%]")

# Sorting

GroveSort!(g::Grove) = GroveSort!(g.Y,TreeInteger(g))
GroveSort!(g::NotGrove) = GroveSort!(convert(Grove,g))
GroveSort!(g::Ar2UI8I,Θ::Array{Int,1}) = GroveSort!(Grove(g),Θ)
GroveSort!(g::Grove,Θ::Array{Int,1}) = (g.Y[:,:] = g.Y[sortperm(Θ),:]; g)
GroveSort = (()->(gs=true; return (tf=gs)->(gs≠tf && (gs=tf; ΥGS()); return gs)))()

# Grafting

function ∨(L::PBTree,R::PBTree) # Graft()
  Ld = L.degr; Rd = R.degr; n = Ld + Rd; G = PBTree(n+1,Array{UInt8,1}(n+1));
  G.Y[Ld+1] = n+1; G.Y[1:Ld] = L.Y[:]; G.Y[Ld+2:Ld+Rd+1] = R.Y[:]; return G; end
∨(L::Ar1UI8I,R::PBTree) = PBTree(L) ∨ R; ∨(L::PBTree,R::Ar1UI8I) = L ∨ PBTree(R)
∨(L::Ar1UI8I,R::Ar1UI8I) = PBTree(L) ∨ PBTree(R)
Graft(x::AbstractPBTree,y::AbstractPBTree) = x ∨ y
function BranchLeft(t::PBTree); fx = findfirst(ξ->(ξ==t.degr),t.Y)
  fx>1 && (return PBTree(fx-1,t.Y[1:fx-1])); return PBTree(0x00,Array{UInt8,1}(0)); end
BranchLeft(t::Ar1UI8I) = BranchLeft(convert(PBTree,t))
function BranchRight(t::PBTree); fx = findfirst(ξ->(ξ==t.degr),t.Y)
  fx<t.degr && (return PBTree(t.Y[fx+1:end])); return PBTree(0x00,Array{UInt8,1}(0)); end
BranchRight(t::Ar1UI8I) = BranchRight(convert(PBTree,t))
LeftInherited(t::PBTree) = BranchRight(t).degr == 0
LeftInherited(t::Ar1UI8I) = LeftInherited(convert(PBTree,t))
RightInherited(t::PBTree) = BranchLeft(t).degr == 0
RightInherited(t::Ar1UI8I) = RightInherited(convert(PBTree,t))
PrimitiveTree(t::PBTree) = LeftInherited(t) || RightInherited(t)
PrimitiveTree(t::Ar1UI8I) = PrimitiveTree(convert(PBTree,t))

# Total Grove Repository

function GroveExtend() # initialize set of total groves
  Y = Array{Grove,1}(1); Y[1]=Grove(0,0,Array{UInt8,2}(0,1));
  R = Array{Array{Int,1},1}(0); return (Y,R); end
function GroveExtend!(Y::Array{Grove,1},R::Array{Array{Int,1},1},deg::UInt8)
  D = UInt8(length(Y)):deg; !isempty(D) && println("Extend Grove Degree Level,");
  for n ∈ D; print(" $n"); τ = 1 # loop over all degree levels
    cn = ceil(Int,n/2); fn = floor(Int,n/2); lYn = length(Y[n].Y[:,1])
    Yn = Grove(Array{UInt8,2}(Cn(n),n)) # initialize total grove
    n==1 ? Yn.Y[τ,1] = n : Yn.Y[τ:lYn,n] = n;  # loop over left-branch grove
    for λ ∈ 1:lYn; Yn.Y[τ,1:n-1] = Y[n].Y[λ,:]; τ += 1; end
    for ν ∈ n-1:-1:cn+1 # loop over (right) root indices
      τ=GroveBuild!(Y[ν].Y,Y[n-ν+1].Y,Yn.Y,ν,n,τ,1,Y[ν].size,1,Y[n-ν+1].size)
    end # loop over center root next (if n is odd)
    cn ≠ fn && for ν = cn
      λs = length(Y[ν].Y[:,1]); Λs = length(Y[n-ν+1].Y[:,1])
      fλs = floor(Int,λs/2); cλs = ceil(Int,λs/2)
      fΛs = floor(Int,Λs/2); cΛs = ceil(Int,Λs/2)
      τ = GroveBuild!(Y[ν].Y,Y[n-ν+1].Y,Yn.Y,ν,n,τ,1,fλs,1,fΛs)
      cΛs ≠ fΛs && (τ = GroveBuild!(Y[ν].Y,Y[n-ν+1].Y,Yn.Y,ν,n,τ,1,fλs,cΛs,cΛs))
      τ = GroveBuild!(Y[ν].Y,Y[n-ν+1].Y,Yn.Y,ν,n,τ,1,fλs,cΛs+1,Λs)
      cλs ≠ fλs && (τ = GroveBuild!(Y[ν].Y,Y[n-ν+1].Y,Yn.Y,ν,n,τ,cλs,cλs,1,Λs))
      τ = GroveBuild!(Y[ν].Y,Y[n-ν+1].Y,Yn.Y,ν,n,τ,cλs+1,λs,1,fΛs)
      cΛs ≠ fΛs && (τ = GroveBuild!(Y[ν].Y,Y[n-ν+1].Y,Yn.Y,ν,n,τ,cλs+1,λs,cΛs,cΛs))
      τ = GroveBuild!(Y[ν].Y,Y[n-ν+1].Y,Yn.Y,ν,n,τ,cλs+1,λs,cΛs+1,Λs); end
    for ν ∈ fn:-1:2 # loop over (right) root indices
      for Λ ∈ 1:Y[n-ν+1].size # loop over right grove
        for λ ∈ 1:Y[ν].size # loop left grove
          Yn.Y[τ,:] = (Y[ν].Y[λ,:] ∨ Y[n-ν+1].Y[Λ,:]).Y; τ += 1
        end; end; end
    Yn.Y[τ:τ+lYn-1,1] = n # loop over right-branch grove
    for Λ ∈ 1:lYn; Yn.Y[τ,2:n] = Y[n].Y[Λ,:]; τ += 1; end
    GroveSort() && (print("|Θ"); push!(R,TreeInteger(TreeBase(Yn))); GroveSort!(Yn,R[n]); sort!(R[n])); push!(Y,Yn); end; !isempty(D) && print("\n"); return Y[deg+1]; end
GroveBuild!(Yν::Array{UInt8,2},Ynν1::Array{UInt8,2},Yn::Array{UInt8,2},ν::UI8I,n::UInt8,τ::Int,λs1::Int,λs2::Int,Λs1::Int,Λs2::Int) = ( #loop left/right grove, Graft
  for λ ∈ λs1:λs2; for Λ ∈ Λs1:Λs2; Yn[τ,:] = (Yν[λ,:] ∨ Ynν1[Λ,:]).Y; τ+=1; end; end; τ)
function GroveInteger!(Y::Array{Grove,1},R::Array{Array{Int,1},1},deg::UInt8)
  D = UInt8(length(R)+1):deg; !isempty(D) && println("Extend Grove Integer Level,")
  for n ∈ D; print(" $n"); push!(R,TreeInteger(TreeBase(Y[n+1]))); end
  !isempty(D) && print("\n"); return R[deg]; end
(Υ,ΥI,ΥGS) = (()->((Y,R)=GroveExtend(); return ((d::UI8I)->(return GroveExtend!(Y,R,UInt8(d))), (d::UI8I)->(GroveExtend!(Y,R,UInt8(d)); !GroveSort() && GroveInteger!(Y,R,UInt8(d)); return d==0 ? Array{Int,1}(0) : R[d]), ()->((Y,R)=GroveExtend())) ))()
  # provides hidden total grove reference

include("morphism.jl")
include("arithmetic.jl")

# Involution

σ(x::Grove) = Grove(x.Y[:,end:-1:1]); σ(x::NotGrove) = σ(Grove(x));
#σ(x::PBTree) = PBTree(x.Y[end:-1:1]); function σ(Y::Array{Grove,1}); γ = length(Y);
#  r = Array{Grove,1}(γ); for n∈1:γ; r[n] = σ(Y[n]); end; return r; end

# Tree Label Print

function GrovePrint(υ::PBTree,μ::BaseTree)
  n=υ.degr; ti=TreeInteger(μ); tin=TreeIndex(n,ti)
  show(convert(Array{Int,1},υ.Y)); print(" ↦ "); for ω ∈ 1:length(υ.Y)
    μ.μ[ω]==[] ? print('∅') : show(convert(Array{Int,1},μ.μ[ω])); end; print(" ↦ ")
  show(ti); print(" or ",tin,"/",Cn(n)); print('\n'); end
GrovePrint(υ::PBTree) = GrovePrint(υ,TreeBase(υ));
GrovePrint(μ::BaseTree) = GrovePrint(TreeLoday(μ),μ)
function GrovePrint(Y::Grove) # given Loday label grove
  for η ∈ 1:Y.size; GrovePrint(PBTree(Y.Y[η,:])); end
  show(GroveBin(Y)); print("\n") end;
function GrovePrint(Y::Array{BaseTree,1}) # given Index label grove
  for η ∈ 1:length(Y); GrovePrint(Y[η]); end; end;
GrovePrint(Y::Array{Grove,1}) = for n ∈ 1:length(Y); GrovePrint(Y[n]); end
#GrovePrint(Y::Array{Array{BaseTree,1},1}) = for n ∈ 1:length(Y); GrovePrint(Y[n]); end
GrovePrint(deg::UI8I) = GrovePrint(Υ(deg)); # given deg

end
