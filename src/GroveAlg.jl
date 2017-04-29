__precompile__()
module GroveAlg

# GroveAlg Copyright (C) 2017 Michael Reed

export PBTree, Grove, GroveBin, ==, ===, GroveSort, GroveSort!, BranchLeft, BranchRight, TreeCheck, GroveCheck, GroveError, TreeIndex, TreeIndexCn, Cn, GroveIndex, GroveBit, TreeInteger, TreeRational, TreeShift, Graft, ⊣, ⊢, +, σ, GrovePart, GrovePrint

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
==(a::Grove,b::Grove) = (a.degr == b.degr && a.size == b.size && a.Y == b.Y)
===(a::Grove,b::Grove) = (GroveSort!(a) == GroveSort!(b))

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
promote_rule(::Type{PBTree},::Type{Ar1UI8I}) = PBTree
promote_rule(::Type{Grove},::Type{Union{Ar1UI8I,Ar2UI8I,PBTree,UI8I}})=Grove
show(io::IO,k::GroveBin) = print(io,"Y$(k.degr) \#$(k.size)/$(Cn(k.degr)), Grove [$(k.ppos)\%] $(k.gbin)")

# Sorting

GroveSort!(g::Grove) = GroveSort!(g.Y,TreeInteger(g))
GroveSort!(g::NotGrove) = GroveSort!(convert(Grove,g))
GroveSort!(g::Ar2UI8I,Θ::Array{Int,1}) = GroveSort!(Grove(g),Θ)
GroveSort!(g::Grove,Θ::Array{Int,1}) = (g.Y[:,:] = g.Y[sortperm(Θ),:]; g)
GroveSort = (()->(gs=true; return (tf=gs)->(gs!=tf && (gs=tf; ΥGS()); return gs)))()

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
LeftInherited(t::Ar1UI8I) = BranchRight(convert(PBTree,t))
RightInherited(t::PBTree) = BranchLeft(t).degr == 0
RightInherited(t::Ar1UI8I) = BranchLeft(convert(PBTree,t))
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
    cn != fn && for ν = cn
      λs = length(Y[ν].Y[:,1]); Λs = length(Y[n-ν+1].Y[:,1])
      fλs = floor(Int,λs/2); cλs = ceil(Int,λs/2)
      fΛs = floor(Int,Λs/2); cΛs = ceil(Int,Λs/2)
      τ = GroveBuild!(Y[ν].Y,Y[n-ν+1].Y,Yn.Y,ν,n,τ,1,fλs,1,fΛs)
      cΛs != fΛs && (τ = GroveBuild!(Y[ν].Y,Y[n-ν+1].Y,Yn.Y,ν,n,τ,1,fλs,cΛs,cΛs))
      τ = GroveBuild!(Y[ν].Y,Y[n-ν+1].Y,Yn.Y,ν,n,τ,1,fλs,cΛs+1,Λs)
      cλs != fλs && (τ = GroveBuild!(Y[ν].Y,Y[n-ν+1].Y,Yn.Y,ν,n,τ,cλs,cλs,1,Λs))
      τ = GroveBuild!(Y[ν].Y,Y[n-ν+1].Y,Yn.Y,ν,n,τ,cλs+1,λs,1,fΛs)
      cΛs != fΛs && (τ = GroveBuild!(Y[ν].Y,Y[n-ν+1].Y,Yn.Y,ν,n,τ,cλs+1,λs,cΛs,cΛs))
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


# Grove Error

TreeCheck(d::UI8I,t::Int) = t > 0 && t <= Cn(d)
TreeCheck(t::PBTree) = TreeCheck(t.degr,TreeIndex(t))
TreeCheck(t::Ar1UI8I) = TreeCheck(convert(PBTree,t))
TreeCheck(g::Grove) = findfirst(x->(x==0),TreeIndex(g)) == 0
TreeCheck(g::NotGrove) = TreeCheck(convert(Grove,g))
GroveCheck(d::UI8I,gi::Integer) = gi >= 0 && gi < 2^Cn(d)
GroveCheck(g::Grove) = GroveCheck(g.degr,GroveIndex(g))
GroveCheck(g::NotGrove) = GroveCheck(convert(Grove,g))
GroveError(n::UI8I) = TreeIndex(n)-sortperm(TreeInteger(n))
GroveError(g::Grove) = [1:g.size...]-sortperm(TreeInteger(g))
GroveError(g::NotGrove) = convert(Grove,g)

# TreeIndex()

function TreeIndex(d::UI8I,l::Int,g::Ar2UI8I); v=ΥI(d); ind=Array{Int,1}(l);
  i=TreeInteger(g); for c ∈ 1:l; ind[c] = findfirst(v.==i[c]); end; return ind; end
TreeIndex(g::Grove) = TreeIndex(g.degr,g.size,g.Y)
TreeIndex(g::Ar2UI8I) = TreeIndex(convert(Grove,g))
TreeIndex(t::Ar1UI8I) = TreeIndex(convert(PBTree,t))
TreeIndex(t::PBTree) = TreeIndex(t.degr,TreeInteger(t))
TreeIndex(d::UI8I,j::Int) = (findfirst(j .== TreeInteger(d))[1])
TreeIndex(deg::UI8I) = [1:Int(Cn(deg))...]
TreeIndexCn(deg::UI8I) = TreeIndex(deg).//Int(Cn(deg))

# GroveIndex() // BitArray

function GroveBit(d::UI8I,l::Int,g::Ar2UI8I)
  v=ΥI(d); s=falses(length(v)); gt = TreeInteger(TreeBase(d,l,g)) # TreeIndex()
  for c ∈ 1:l; s |= (v.==gt[c]); end; return s; end
GroveBit(g::Grove) = GroveBit(g.degr,g.size,g.Y)
GroveBit(g::Union{Ar1UI8I,Ar2UI8I}) = GroveBit(Grove(g))
GroveBit(s::Integer) = '1' .== flipdim(collect(bin(s)),1)
function GroveIndex(d::UI8I,l::Int,g::Ar2UI8I);
  s=BigInt(0); try; v=ΥI(d); gt = TreeInteger(TreeBase(d,l,g)) # TreeIndex()
  for c ∈ 1:l; s += BigInt(2)^(findfirst(v.==gt[c])-1); end;
catch TreeCheck(g) ? throw(DomainError()) : s=BigInt(-1); end; return s; end
GroveIndex(g::Grove) = GroveIndex(g.degr,g.size,g.Y)
GroveIndex(g::NotGrove) = GroveIndex(convert(Grove,g))
function GroveIndex(b::BitArray); s = BigInt(0); gi = find(b)
  for c ∈ length(gi):-1:1; s += BigInt(2)^(gi[c]-1); end; return s; end

# TreeLoday()

TreeLoday(d::UI8I,s::BitArray) = TreeLoday(d,find(s))
function TreeLoday(d::UI8I,ind::Array{Int,1}) # from TreeIndex
  y = Υ(d).Y; g = Grove(d,Array{UInt8,2}(length(ind),d)); c = 1;
  d==0 ? (g=Υ(0)) : (for i ∈ ind; g.Y[c,:] = y[i,:]; c +=1; end); return g; end
function TreeLoday(υ::BaseTree) # from TreeBase
  γ = UInt8(length(υ.μ)); L = PBTree(γ,Array{UInt8,1}(γ)); γ1 = γ+0x01
  for ω ∈ 0x01:γ; L.Y[υ.μ[ω]] = γ1-ω; end; return L; end;
function TreeLoday(Y::Array{BaseTree,1})
  γ = length(Y); L = Grove(Array{UInt8,2}(γ,length(Y[1].μ[:])))
  for η ∈ 1:γ; L.Y[η,:] = (TreeLoday(Y[η])).Y; end; return L; end;
function TreeLoday(Y::Array{Array{BaseTree,1},1})
  γ = length(Y); L = Array{Grove,1}(γ)
  for n ∈ 1:γ; L[n] = TreeLoday(Y[n]); end; return L; end;
TreeLoday(deg::UI8I,ind::Int) = TreeLoday(deg,[ind]) # TreeIndex
TreeLoday(deg::UI8I) = Υ(deg); # from degree
TreeLoday(d::UI8I,s::Integer) = Grove(d,s)

# TreeBase()

TreeBase(d::UI8I,s::Integer) = TreeBase(d,GroveBit(find(s)))
function TreeBase(d::UI8I,ind::BitArray) # from TreeIndex
  y = Υ(d).Y; g = Grove(Array{UInt8,2}(length(ind),d)); c = 1;
  for i ∈ ind; g.Y[c,:] = y[i,:]; c +=1; end; return TreeBase(g); end
function TreeBase(d::UI8I,υ::Ar1UI8I) # index label for tree
  μ = BaseTree(Array{Array{UInt8,1},1}(d)); γ1 = d +0x01
  for ω ∈ 0x01:d; μ.μ[ω] = find(ξ->(ξ==γ1-ω),υ); end; return μ; end;
TreeBase(υ::Ar1UI8I) = TreeBase(length(υ),υ)
TreeBase(t::PBTree) = TreeBase(t.degr,t.Y)
function TreeBase(d::UI8I,γ::Int,g::Ar2UI8I); μ = Array{BaseTree,1}(γ)
  for η ∈ 1:γ; μ[η] = TreeBase(d,g[η,:]); end; return μ; end;
TreeBase(g::Ar2UI8I) = TreeBase(GroveDeg(g),length(g[:,1]),g)
TreeBase(g::Grove) = TreeBase(g.degr,g.size,g.Y)
function TreeBase(Y::Union{Array{Ar2UI8I,1},Array{Grove,1}})
  γ = length(Y); μ = Array{Array{BaseTree,1},1}(γ)
  for n ∈ 1:γ; μ[n] = TreeBase(Y[n]); end; return μ; end;
TreeBase(deg::UI8I) = TreeBase(Υ(deg));

# TreeBase -> Max Integer

ΘMax = (()->(s = Array{Int,1}(1); s[1] = 1; return ((d::UI8I)->(
  for n ∈ length(s)+1:d; δ = d-0x01; σ = 0
    for i ∈ n:-1:1; σ += i*10^δ; δ -= 0x01; end; push!(s,σ); end; s[d])) ))()
function ΘInt(μ::Array{Array{UInt8,1},1}); d = UInt8(length(μ)-1); s = 0
  for n ∈ μ; for i ∈ n; s += i*10^d; d -= 0x01; end; end; return s; end

# TreeInteger()

TreeInteger(d::UI8I,s::Integer) = TreeInteger(d,GroveBit(s))
TreeInteger(d::UI8I,s::BitArray) = TreeInteger(d,find(s))
function TreeInteger(d::UI8I,ind::Array{Int,1}) # from TreeIndex
  y = ΥI(d); g = Array{Int,1}(length(ind)); c = 1;
  for i ∈ ind; g[c] = y[i]; c +=1; end; return g; end
TreeInteger(μ::BaseTree) = ΘMax(length(μ.μ))-ΘInt(μ.μ)
function TreeInteger(Y::Array{BaseTree,1})
  γ = length(Y); i = Array{Int,1}(γ)
  for n ∈ 1:γ; i[n] = TreeInteger(Y[n]); end; return i; end;
function TreeInteger(Y::Array{Array{BaseTree,1},1})
  γ = length(Y); i = Array{Array{Int,1},1}(γ)
  for n ∈ 1:γ; i[n] = TreeInteger(Y[n]); end; return i; end;
TreeInteger(g::Union{Grove,NotGrove,Array}) = TreeInteger(TreeBase(g))
TreeInteger(deg::UI8I) = ΥI(deg);

# TreeRational()

TreeRational(d::UI8I,s::Integer) = TreeRational(d,GroveBit(find(s)))
function TreeRational(d::UI8I,ind::BitArray) # from TreeIndex
  y = ΥI(d); g = Array{Int,1}(length(ind)); c = 1;
  for i ∈ ind; g[c] = y[i]; c +=1; end; return TreeRational(d,g); end
TreeRational(μ::BaseTree) = (s=TreeShift(); s+((-1)^s)*ΘInt(μ.μ)//ΘMax(length(μ.μ)))
TreeRational(deg::UI8I,Θ::Int) = (s=TreeShift(); s+((-1)^s)*Θ//ΘMax(deg))
TreeRational(deg::UI8I,Θ::Array{Int,1}) = (s=TreeShift(); s+((-1)^s)*Θ.//ΘMax(deg))
function TreeRational(Y::Array{BaseTree,1})
  γ = length(Y); r = Array{Rational,1}(γ)
  for n ∈ 1:γ; r[n] = TreeRational(Y[n]); end; return r; end;
function TreeRational(Y::Array{Array{BaseTree,1},1})
  γ = length(Y); r = Array{Array{Rational,1},1}(γ)
  for n ∈ 1:γ; r[n] = TreeRational(Y[n]); end; return r; end;
TreeRational(deg::UI8I) = TreeRational(deg,TreeInteger(deg));
TreeRational(υ::Any) = TreeRational(TreeBase(υ));
TreeShift = (()->(gs=true; return (tf=gs)->(gs!=tf && (gs=tf); return Int(gs))))()

# Arithmetic (Left)

function ⊣(x::PBTree,y::PBTree); sm = BranchRight(x)+y; blx = BranchLeft(x)
  isempty(sm.Y) && (return Grove(blx ∨ Array{UInt8,1}(0)))
  ls = sm.size; addl = Grove(Array{UInt8,2}(ls,sm.degr + blx.degr + 1))
  for i ∈ 1:ls; addl.Y[i,:] = (blx ∨ sm.Y[i,:]).Y; end; return addl; end
function ⊣(x::Grove,y::PBTree); γ = x.size; gr = Array{Array,1}(γ)
  for i ∈ 1:γ; gr[i] = (PBTree(x.Y[i,:]) ⊣ y).Y; end; return Grove(vcat(gr...)); end
function ⊣(x::PBTree,y::Grove); γ = y.size; gr = Array{Array,1}(γ)
  for i ∈ 1:γ; gr[i] = (x ⊣ PBTree(y.Y[i,:])).Y; end; return Grove(vcat(gr...)); end
function ⊣(x::Grove,y::Grove); γ = x.size; gr = Array{Array,1}(γ)
  for i ∈ 1:γ; gr[i] = (PBTree(x.Y[i,:]) ⊣ y).Y; end; return Grove(vcat(gr...)); end
⊣(x::NotGrove,y::NotGrove) = Grove(x) ⊣ Grove(y)
⊣(x::Ar1UI8I,y::Ar1UI8I) = PBTree(x) ⊣ PBTree(y)

# Arithmetic (Right)

function ⊢(x::PBTree,y::PBTree); sm = x+BranchLeft(y); bry = BranchRight(y)
  isempty(sm.Y) && (return Grove(Array{UInt8,1}(0) ∨ bry))
  ls = sm.size; addr = Grove(Array{UInt8,2}(ls,sm.degr + bry.degr + 1))
  for i ∈ 1:ls; addr.Y[i,:] = (sm.Y[i,:] ∨ bry).Y; end; return addr; end
function ⊢(x::PBTree,y::Grove); γ = y.size; gr = Array{Array,1}(γ)
  for i ∈ 1:γ; gr[i] = (x ⊢ PBTree(y.Y[i,:])).Y; end; return Grove(vcat(gr...)); end
function ⊢(x::Grove,y::PBTree); γ = x.size; gr = Array{Array,1}(γ)
  for i ∈ 1:γ; gr[i] = (PBTree(x.Y[i,:]) ⊢ y).Y; end; return Grove(vcat(gr...)); end
function ⊢(x::Grove,y::Grove); γ = x.size; gr = Array{Array,1}(γ);
  for i ∈ 1:γ; gr[i] = (PBTree(x.Y[i,:]) ⊢ y).Y; end; return Grove(vcat(gr...)); end
⊢(x::NotGrove,y::NotGrove) = Grove(x) ⊢ Grove(y)
⊢(x::Ar1UI8I,y::Ar1UI8I) = PBTree(x) ⊢ PBTree(y)

# Grove addition

function +(x::Grove,y::Grove)
  isempty(x.Y) && (return y); isempty(y.Y) && (return x)
  lx = x.size; ly = y.size; ij = Array{Array,2}(lx,ly)
  for i ∈ 1:lx
    for j ∈ 1:ly
      l = Grove(PBTree(x.Y[i,:]) ⊣ PBTree(y.Y[j,:]))
      r = Grove(PBTree(x.Y[i,:]) ⊢ PBTree(y.Y[j,:]))
      ij[i,j] = vcat(l.Y,r.Y); end; end; return Grove(vcat(ij...)); end
+(x::Grove,y::NotGrove) = x + convert(Grove,y)
+(x::NotGrove,y::Grove) = convert(Grove,x) + y
+(x::Union{GroveBin,PBTree},y::Union{GroveBin,PBTree}) = Grove(x) + Grove(y);

# Grove Composition

function GroveComposition(n::Int,η::Int=n)
  G = GroveSums(n); u=1
  !isempty(G) && (return G)
  n < η && (u = 2^Cn(n))
  n != 0 && n < η && for i∈1:u-1; push!(G,[(n,i),(n,i)]); end;
  for s ∈ n-1:-1:1
    for i ∈ 1:2^Cn(s)-1
      g = GrovePart(n-s,η); gsi = Grove(s,i);
      for r ∈ 1:length(g)
        deg = g[r][end][1]; gi = g[r][end][2]
        sm = gsi + Grove(deg,gi);
        push!(G,g[r]); unshift!(G[u],(s,i));
        G[u][end] = (UInt8(deg+s),GroveIndex(sm))
        u += 1; end; end; end;
  GroveStore(n,G); return G; end
(GroveStore,GL,GroveSums) = (()->(GG=Array{Array{Array{Tuple{UInt8,Integer},1},1},1}(0); return ((n::Int,g::Array{Array{Tuple{UInt8,Integer},1},1})->(n==length(GG)+1 && (g)), ()->(return length(GG)), (n::Int)->(n<=length(GG) && n>0x00 ? (return GG[n]) : (return Array{Array{Tuple{UInt8,Integer},1},1}(0) )  ) )))()

# Involution

σ(x::PBTree) = x.Y[end:-1:1]; σ(x::Grove) = x.Y[:,end:-1:1]; σ(x::Any) = σ(Grove(x))
function σ(Y::Array{Grove,1}); γ = length(Y);
  r = Array{Grove,1}(γ); for n∈1:γ; r[n] = σ(Y[n]); end; return r; end

# Tree Label Print

function GrovePrint(υ::PBTree,μ::BaseTree)
  n=υ.degr; ti=TreeInteger(μ); tin=TreeIndex(n,ti)
  show(convert(Array{Int,1},υ.Y)); print(" ↦ "); for ω ∈ 1:length(υ.Y)
    μ.μ[ω]==[] ? print('∅') : show(convert(Array{Int,1},μ.μ[ω])); end; print(" ↦ ")
  show(ti); print(" or ",tin,"/",Cn(n)); end
GrovePrint(υ::PBTree) = GrovePrint(υ,TreeBase(υ));
GrovePrint(μ::BaseTree) = GrovePrint(TreeLoday(μ),μ);
function GrovePrint(Y::Grove) # given Loday label grove
  for η ∈ 1:Y.size; GrovePrint(PBTree(Y.Y[η,:])); print('\n'); end
  show(GroveBin(Y)); print("\n") end;
function GrovePrint(Y::Array{BaseTree,1}) # given Index label grove
  for η ∈ 1:length(Y); GrovePrint(Y[η]); print('\n'); end; end;
GrovePrint(Y::Array{Grove,1}) = for n ∈ 1:length(Y); GrovePrint(Y[n]); end
GrovePrint(Y::Array{Array{BaseTree,1},1}) = for n ∈ 1:length(Y); GrovePrint(Y[n]); end
GrovePrint(deg::UI8I) = GrovePrint(Υ(deg)); # given deg

end
