# This file is part of GroveAlg.jl. It is licensed under the GPL license
# GroveAlg Copyright (C) 2017 Michael Reed

export TreeCheck, GroveCheck, GroveError, TreeIndex, TreeIndexCn, GroveIndex, GroveBit, TreeInteger, TreeRational, TreeShift

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
GroveError(g::NotGrove) = convert(Grove,g) |> GroveError

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
#function TreeLoday(Y::Array{Array{BaseTree,1},1})
#  γ = length(Y); L = Array{Grove,1}(γ)
#  for n ∈ 1:γ; L[n] = TreeLoday(Y[n]); end; return L; end;
#TreeLoday(deg::UI8I,ind::Int) = TreeLoday(deg,[ind]) # TreeIndex
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
#function TreeBase(Y::Union{Array{Ar2UI8I,1},Array{Grove,1}})
#  γ = length(Y); μ = Array{Array{BaseTree,1},1}(γ)
#  for n ∈ 1:γ; μ[n] = TreeBase(Y[n]); end; return μ; end;
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
#function TreeInteger(Y::Array{Array{BaseTree,1},1})
#  γ = length(Y); i = Array{Array{Int,1},1}(γ)
#  for n ∈ 1:γ; i[n] = TreeInteger(Y[n]); end; return i; end;
TreeInteger(g::Union{Grove,NotGrove,Array}) = TreeInteger(TreeBase(g))
TreeInteger(deg::UI8I) = ΥI(deg);

# TreeRational()

TreeRational(d::UI8I,s::Integer) = TreeRational(d,GroveBit(s))
function TreeRational(d::UI8I,indb::BitArray) # from TreeIndex
  y = ΥI(d); c = 1; ind = find(indb); g = Array{Int,1}(length(ind));
  for i ∈ ind; g[c] = y[i]; c +=1; end; return TreeRational(d,g); end
TreeRational(μ::BaseTree) = (s=TreeShift(); s+((-1)^s)*ΘInt(μ.μ)//ΘMax(length(μ.μ)))
#TreeRational(deg::UI8I,Θ::Int) = (s=TreeShift(); 1-s-((-1)^s)*Θ//ΘMax(deg))
TreeRational(deg::UI8I,Θ::Array{Int,1}) = (s=TreeShift(); 1-s-((-1)^s)*Θ.//ΘMax(deg))
function TreeRational(Y::Array{BaseTree,1}); γ = length(Y); r = Array{Rational,1}(γ)
  for n ∈ 1:γ; r[n] = TreeRational(Y[n]); end; return r; end;
#function TreeRational(Y::Array{Array{BaseTree,1},1})
#  γ = length(Y); r = Array{Array{Rational,1},1}(γ)
#  for n ∈ 1:γ; r[n] = TreeRational(Y[n]); end; return r; end;
TreeRational(deg::UI8I) = TreeRational(deg,TreeInteger(deg));
TreeRational(υ::Any) = TreeRational(TreeBase(υ));
TreeShift = (()->(gs=true; return (tf=gs)->(gs≠tf && (gs=tf); return Int(gs))))()
