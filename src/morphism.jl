# This file is part of Dendriform.jl. It is licensed under the GPL license
# Dendriform Copyright (C) 2017 Michael Reed

export treecheck, grovecheck, treeindex, treeindexCn, groveindex, grovebit, treeshift

# Tree Structure

LeftInherited(t::PBTree) = rightbranch(t).degr == 0
LeftInherited(t::Ar1UI8I) = LeftInherited(convert(PBTree,t))
RightInherited(t::PBTree) = leftbranch(t).degr == 0
RightInherited(t::Ar1UI8I) = RightInherited(convert(PBTree,t))
PrimitiveTree(t::PBTree) = LeftInherited(t) || RightInherited(t)
PrimitiveTree(t::Ar1UI8I) = PrimitiveTree(convert(PBTree,t))

# Grove Error

treecheck(d::UI8I,t::Int) = t > 0 && t <= Cn(d)
treecheck(t::PBTree) = treecheck(t.degr,treeindex(t))
treecheck(t::Ar1UI8I) = treecheck(convert(PBTree,t))
treecheck(g::Grove) = findfirst(x->(x==0),treeindex(g)) == 0
treecheck(g::NotGrove) = treecheck(convert(Grove,g))
grovecheck(d::UI8I,gi::Integer) = gi >= 0 && gi < 2^Cn(d)
grovecheck(g::Grove) = grovecheck(g.degr,groveindex(g))
grovecheck(g::NotGrove) = grovecheck(convert(Grove,g))
GroveError(n::UI8I) = treeindex(n)-sortperm(TreeInteger(n))
GroveError(g::Grove) = [1:g.size...]-sortperm(TreeInteger(g))
GroveError(g::NotGrove) = convert(Grove,g) |> GroveError

# treeindex

function treeindex(d::UI8I,l::Int,g::Ar2UI8I)
    v=ΥI(d)
    ind=Array{Int,1}(l)
    i=TreeInteger(g)
    for c ∈ 1:l
        ind[c] = findfirst(v.==i[c])
    end
    return ind
end

treeindex(g::Grove) = treeindex(g.degr,g.size,g.Y)
treeindex(g::Ar2UI8I) = treeindex(convert(Grove,g))
treeindex(t::Ar1UI8I) = treeindex(convert(PBTree,t))
treeindex(t::PBTree) = treeindex(t.degr,TreeInteger(t))
treeindex(d::UI8I,j::Int) = (findfirst(j .== TreeInteger(d))[1])
treeindex(deg::UI8I) = [1:Int(Cn(deg))...]
treeindexCn(deg::UI8I) = treeindex(deg).//Int(Cn(deg))

# grovebit

function grovebit(d::UI8I,l::Int,g::Ar2UI8I)
    v=ΥI(d)
    s=falses(length(v))
    gt = TreeInteger(TreeBase(d,l,g)) # treeindex()
    for c ∈ 1:l
        s .|= (v.==gt[c])
    end
    return s
end

grovebit(g::Grove) = grovebit(g.degr,g.size,g.Y)
grovebit(g::Union{Ar1UI8I,Ar2UI8I}) = grovebit(Grove(g))
grovebit(g::GroveBin) = grovebit(g.degr,g.gbin)

function grovebit(d::UI8I,s::Integer)
    gb = '1' .== flipdim(collect(bin(s)),1)
    return vcat(gb,falses(Cn(d)-length(gb)))
end

# groveindex

function groveindex(d::UI8I,l::Int,g::Ar2UI8I)
    s=BigInt(0)
    try
        v=ΥI(d)
        gt = TreeInteger(TreeBase(d,l,g)) # treeindex()
        for c ∈ 1:l
            s += BigInt(2)^(findfirst(v.==gt[c])-1)
        end
    catch
        treecheck(g) ? throw(DomainError()) : s=BigInt(-1)
    end
    return s
end

groveindex(g::Grove) = groveindex(g.degr,g.size,g.Y)
groveindex(g::NotGrove) = groveindex(convert(Grove,g))

function groveindex(b::BitArray)
    s = BigInt(0)
    gi = find(b)
    for c ∈ length(gi):-1:1
        s += BigInt(2)^(gi[c]-1)
    end
    return s
end

# TreeLoday

TreeLoday(d::UI8I,s::BitArray) = TreeLoday(d,find(s))

function TreeLoday(d::UI8I,ind::Array{Int,1}) # from treeindex
    y = Υ(d).Y
    g = Grove(d,Array{UInt8,2}(length(ind),d))
    c = 1
    if d == 0
        g=Υ(0)
    else
        for i ∈ ind
            g.Y[c,:] = y[i,:]
            c +=1
        end
    end
    return g
end

function TreeLoday(υ::BaseTree) # from TreeBase
    γ = UInt8(length(υ.μ))
    L = PBTree(γ,Array{UInt8,1}(γ))
    γ1 = γ+0x01
    for ω ∈ 0x01:γ
        L.Y[υ.μ[ω]] = γ1-ω
    end
    return L
end

function TreeLoday(Y::Array{BaseTree,1})
    γ = length(Y)
    L = Grove(Array{UInt8,2}(γ,length(Y[1].μ[:])))
    for η ∈ 1:γ
        L.Y[η,:] = (TreeLoday(Y[η])).Y
    end
    return L
end

TreeLoday(deg::UI8I) = Υ(deg); # from degree
TreeLoday(d::UI8I,s::Integer) = Grove(d,s)

# TreeBase

TreeBase(d::UI8I,s::Integer) = TreeBase(d,grovebit(find(s)))

function TreeBase(d::UI8I,ind::BitArray) # from treeindex
    y = Υ(d).Y
    g = Grove(Array{UInt8,2}(length(ind),d))
    c = 1
    for i ∈ ind
        g.Y[c,:] = y[i,:]
        c +=1
    end
    return TreeBase(g)
end

function TreeBase(d::UI8I,υ::Ar1UI8I) # index label for tree
    μ = BaseTree(Array{Array{UInt8,1},1}(d))
    γ1 = d +0x01
    for ω ∈ 0x01:d
        μ.μ[ω] = find(ξ->(ξ==γ1-ω),υ)
    end
    return μ
end

TreeBase(υ::Ar1UI8I) = TreeBase(length(υ),υ)
TreeBase(t::PBTree) = TreeBase(t.degr,t.Y)

function TreeBase(d::UI8I,γ::Int,g::Ar2UI8I)
    μ = Array{BaseTree,1}(γ)
    for η ∈ 1:γ
        μ[η] = TreeBase(d,g[η,:])
    end
    return μ
end

TreeBase(g::Ar2UI8I) = TreeBase(GroveDeg(g),length(g[:,1]),g)
TreeBase(g::Grove) = TreeBase(g.degr,g.size,g.Y)
TreeBase(deg::UI8I) = TreeBase(Υ(deg));

# TreeBase -> Max Integer

ΘMax = ( () -> begin
        s = Array{Int,1}(1)
        s[1] = 1
        return ( (d::UI8I) -> begin
                for n ∈ length(s)+1:d
                    δ = d-0x01
                    σ = 0
                    for i ∈ n:-1:1
                        σ += i*10^δ
                        δ -= 0x01
                    end
                    push!(s,σ)
                end
                s[d]
            end)
    end)()

function ΘInt(μ::Array{Array{UInt8,1},1})
    d = UInt8(length(μ)-1)
    s = 0
    for n ∈ μ
        for i ∈ n
            s += i*10^d
            d -= 0x01
    end; end
    return s
end

# TreeInteger

TreeInteger(d::UI8I,s::Integer) = TreeInteger(d,grovebit(d,s))
TreeInteger(d::UI8I,s::BitArray) = TreeInteger(d,find(s))

function TreeInteger(d::UI8I,ind::Array{Int,1}) # from treeindex
    y = ΥI(d)
    g = Array{Int,1}(length(ind))
    c = 1
    for i ∈ ind
        g[c] = y[i]
        c +=1
    end
    return g
end

TreeInteger(μ::BaseTree) = ΘMax(length(μ.μ))-ΘInt(μ.μ)

function TreeInteger(Y::Array{BaseTree,1})
    γ = length(Y)
    i = Array{Int,1}(γ)
    for n ∈ 1:γ
        i[n] = TreeInteger(Y[n])
    end
    return i
end

TreeInteger(g::Union{Grove,NotGrove,Array}) = TreeInteger(TreeBase(g))
TreeInteger(deg::UI8I) = ΥI(deg);

# TreeRational

TreeRational(d::UI8I,s::Integer) = TreeRational(d,grovebit(d,s))

function TreeRational(d::UI8I,indb::BitArray) # from treeindex
    y = ΥI(d)
    c = 1
    ind = find(indb)
    g = Array{Int,1}(length(ind))
    for i ∈ ind
        g[c] = y[i]
        c +=1
    end
    return TreeRational(d,g)
end

TreeRational(μ::BaseTree) = (s=treeshift(); s+((-1)^s)*ΘInt(μ.μ)//ΘMax(length(μ.μ)))
TreeRational(deg::UI8I,Θ::Array{Int,1}) = (s=treeshift(); 1-s-((-1)^s)*Θ.//ΘMax(deg))

function TreeRational(Y::Array{BaseTree,1})
    γ = length(Y)
    r = Array{Rational,1}(γ)
    for n ∈ 1:γ
        r[n] = TreeRational(Y[n])
    end
    return r
end

TreeRational(deg::UI8I) = TreeRational(deg,TreeInteger(deg))
TreeRational(υ::Any) = TreeRational(TreeBase(υ))

treeshift = ( () -> begin
        gs=true
        return (tf=gs)->(gs≠tf && (gs=tf); return Int(gs))
    end)()

# Inequalities

<(x::AbstractPBTree,y::AbstractPBTree) = treeindex(x) < treeindex(y)
>(x::AbstractPBTree,y::AbstractPBTree) = treeindex(x) > treeindex(y)
≤(x::AbstractPBTree,y::AbstractPBTree) = treeindex(x) ≤ treeindex(y)
≥(x::AbstractPBTree,y::AbstractPBTree) = treeindex(x) ≥ treeindex(y)
<(x::PureGrove,y::PureGrove) = groveindex(x) < groveindex(y)
>(x::PureGrove,y::PureGrove) = groveindex(x) > groveindex(y)
≤(x::PureGrove,y::PureGrove) = groveindex(x) ≤ groveindex(y)
≥(x::PureGrove,y::PureGrove) = groveindex(x) ≥ groveindex(y)
