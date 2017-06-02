# This file is part of GroveAlg.jl. It is licensed under the GPL license
# GroveAlg Copyright (C) 2017 Michael Reed

export ⊣, ⊢, +, σ, GroveComposition

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

function GroveCompos(n::Int,η::Int=n)
  G = GroveSums(n); u=1
  !isempty(G) && (return G)
  n < η && (u = 2^Cn(n))
  n ≠ 0 && n < η && for i∈1:u-1; push!(G,[GroveBin(Grove(n,i)),GroveBin(Grove(n,i))]); end;
  for s ∈ n-1:-1:1
    for i ∈ 1:2^Cn(s)-1
      g = GroveCompos(n-s,η); gsi = Grove(s,i); gbsi = GroveBin(gsi)
      for r ∈ 1:length(g)
        sm = gsi + Grove(g[r][end].degr,g[r][end].gbin);
        push!(G,g[r]); unshift!(G[u],gbsi);
        G[u][end] = GroveBin(sm); u += 1; end; end; end;
  c = Dict{Integer,Array{Integer,1}}(); for k ∈ 1:length(G); ind = G[k][end].gbin
    try; c[ind]; catch; push!(c,ind=>Array{Integer,1}(0)); end; push!(c[ind],k); end
  GroveStore(n,G,c); return G; end
(GroveStore,GroveComp,GroveSums) = (()->(GG=Array{Array{Array{GroveBin,1},1},1}(0); CC=Array{Dict{Integer,Array{Integer,1}},1}(0); return ((n::Int,g::Array{Array{GroveBin,1},1},c::Dict{Integer,Array{Integer,1}})->(n==length(GG)+1 && (push!(GG,deepcopy(g)); push!(CC,deepcopy(c)))), (n::Int)->(GroveCompos(n); (return deepcopy(CC[n]))), (n::Int)->(n<=length(GG) && n>0x00 ? (return deepcopy(GG[n])) : (return Array{Array{GroveBin,1},1}(0) )  ) )))()
function GroveComposition(d::UI8I,ind::Integer)
  show(GroveBin(Grove(d,ind))); gi = GroveComp(d); com = GroveCompos(d)
  try; gi[ind]; catch; print(" has 1 composition (iteslf)\n"); return 1; end
  lg = length(gi[ind]); print(" has $(lg+1) compositions\n")
  for k ∈ 1:lg; print("("); show(com[gi[ind][k]][1])
    for t ∈ 2:length(com[gi[ind][k]])-1; print(") + ("); show(com[gi[ind][k]][t]); end
    print(")\n"); end; return lg+1; end
