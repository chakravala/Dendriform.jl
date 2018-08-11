# This file is part of Dendriform.jl. It is licensed under the GPL license
# Dendriform Copyright (C) 2017 Michael Reed

export ⋖, ⋗, posetnext, posetprev, between, ⊴, over, under

# poset coverings

function posetnext_list(t::PBTree)
    g = Array{PBTree,1}()
    λ = left(t)
    ρ = right(t)
    x = left(λ) ∨ (right(λ) ∨ ρ)
    x.degr == t.degr && push!(g,x)
    if λ.degr ≠ 0x00
        gλ = posetnext_list(λ)
        for i ∈ 1:length(gλ)
            gλ[i] = gλ[i] ∨ ρ
        end
        push!(g,gλ...)
    end
    if ρ.degr ≠ 0x00
        gρ= posetnext_list(ρ)
        for i ∈ 1:length(gρ)
            gρ[i] = λ ∨ gρ[i]
        end
        push!(g,gρ...)
    end
    return g
end

"""
    posetnext(::AbstractPBTree)

Returns a Grove that covers the given tree
"""
posetnext(t::PBTree) = t |> posetnext_list |> Grove
posetnext(t::AbstractPBTree) = t |> PBTree |> posetnext

"""
    ⋖(a::AbstractPBTree, b::AbstractPBTree)

Returns Bool that tells if b covers a in Tamari partial order
"""
⋖(a::PBTree,b::PBTree) = b ∈ posetnext_list(a)
⋖(a::AbstractPBTree,b::AbstractPBTree) = PBTree(a) ⋖ PBTree(b)

function posetprev_list(t::PBTree)
    g = Array{PBTree,1}()
    λ = left(t)
    ρ = right(t)
    x = (λ ∨ left(ρ)) ∨ right(ρ)
    x.degr == t.degr && push!(g,x)
    if λ.degr ≠ 0x00
        gλ = posetprev_list(λ)
        for i ∈ 1:length(gλ)
            gλ[i] = gλ[i] ∨ ρ
        end
        push!(g,gλ...)
    end
    if ρ.degr ≠ 0x00
        gρ = posetprev_list(ρ)
        for i ∈ 1:length(gρ)
            gρ[i] = λ ∨ gρ[i]
        end
        push!(g,gρ...)
    end
    return g
end

"""
    posetprev(::AbstractPBTree)

Returns a Grove that covers the given tree
"""
posetprev(t::PBTree) = t |> posetprev_list |> Grove
posetprev(t::AbstractPBTree) = t |> PBTree |> posetprev

"""
    ⋗(a::AbstractPBTree, b::AbstractPBTree)

Returns Bool that tells if a covers b in Tamari partial order
"""
⋗(a::PBTree,b::PBTree) = b ∈ posetprev_list(a)
⋗(a::AbstractPBTree,b::AbstractPBTree) = PBTree(a) ⋗ PBTree(b)

# partial order

"""
    <(a::AbstractPBTree, b::AbstractPBTree)

Returns Bool that tells if a < b in Tamari partial order
"""
function <(a::PBTree,b::PBTree)
    h = posetnext_list(a)
    b ∈ h && return true
    less = false
    i = 0
    while !less & (i < length(h))
        i += 1
        less = h[i] < b
    end
    return less
end

<(a::AbstractPBTree,b::AbstractPBTree) = PBTree(a) < PBTree(b)

"""
    ≤(a::AbstractPBTree, b::AbstractPBTree)

Returns Bool that tells if a ≤ b in Tamari partial order
"""
≤(a::PBTree,b::PBTree) = (a == b) || (a < b)
≤(a::AbstractPBTree,b::AbstractPBTree) = PBTree(a) ≤ PBTree(b)

"""
    >(a::AbstractPBTree, b::AbstractPBTree)

Returns Bool that tells if a > b in Tamari partial order
"""
function >(a::PBTree,b::PBTree)
    h = posetprev_list(a)
    b ∈ h && return true
    gtr = false
    i = 0
    while !gtr & (i < length(h))
        i += 1
        gtr = h[i] > b
    end
    return gtr
end

>(a::AbstractPBTree,b::AbstractPBTree) = PBTree(a) > PBTree(b)

"""
    ≥(a::AbstractPBTree, b::AbstractPBTree)

Returns Bool that tells if a ≥ b in Tamari partial order
"""
≥(a::PBTree,b::PBTree) = (a == b) || (a > b)
≥(a::AbstractPBTree,b::AbstractPBTree) = PBTree(a) ≥ PBTree(b)

# between order

function between_list(a::PBTree,b::PBTree)
    g = Array{PBTree,1}()
    a == b && return push!(g,b)
    h = posetnext_list(a)
    push!(g,a)
    for i ∈ 1:length(h)
        if h[i] ≤ b
            p = between_list(h[i],b)
            for t ∈ p
                t ∉ g && push!(g,t)
            end
        end
    end
    return length(g) > 1 ? g : Array{PBTree,1}()
end

"""
    between(a::AbstractPBTree,b::AbstractPBTree)

Returns Grove of trees ordered between a and b
"""
between(a::PBTree,b::PBTree) = between_list(a,b) |> Grove
between(a::AbstractPBTree,b::AbstractPBTree) = between(PBTree(a),PBTree(b))

"""
    ⊴(a::AbstractPBTree,b::AbstractPBTree)

Returns Grove of trees ordered between a and b
"""
⊴(a::AbstractPBTree,b::AbstractPBTree) = between(a,b)

# over and under

import Base: /, \

"""
    over(a::AbstractPBTree, b::AbstractPBTree)

Returns PBTree obtained from a over b operation
"""
over(x::PBTree,y::PBTree) = y.degr > 0 ? over(x,left(y)) ∨ right(y) : x
over(x::AbstractPBTree,y::AbstractPBTree) = over(PBTree(x),PBTree(y))

"""
    /(a::PBTree, b::PBTree)

Returns PBTree obtained from a over b operation
"""
/(x::PBTree,y::PBTree) = over(x,y)

"""
    under(a::AbstractPBTree, b::AbstractPBTree)

Returns PBTree obtained from a under b operation
"""
under(x::PBTree,y::PBTree) = x.degr > 0 ? left(x) ∨ under(right(x),y) : y
under(x::AbstractPBTree,y::AbstractPBTree) = under(PBTree(x),PBTree(y))

"""
    \\(a::PBTree, b::PBTree)

Returns PBTree obtained from a under b operation
"""
\(x::PBTree,y::PBTree) = under(x,y)

# interval tools

function intervals(d::Int)
    g = Array{BigInt,1}()
    for i∈1:Int(Cn(d)),j∈1:Int(Cn(d))
        t = PBTree(d,i) ⊴ PBTree(d,j)
        t ≠ Grove(0) && push!(g,GroveBin(t).gbin)
    end
    return sort!(g)
end

function intcomp(d::Int)
    ins = intervals(d)
    lins = length(ins)
    cc = zeros(Int,lins)
    cn = 0
    for q ∈ 1:d-1
        for i ∈ 1:Int(groveindex(Grove(q))), j ∈ 1:Int(groveindex(Grove(d-q)))
            z = GroveBin(Grove(q,i) + Grove(d-q,j))
            f = findall(s->s==z.gbin, ins)
            length(f) == 0 ? (cn += 1) : (cc[f] .+= 1)
        end
    end
    @info "Non-intervals: $cn"
    return cc
end

function intcompt(d::Int)
    ins = intervals(d)
    lins = length(ins)
    cc = zeros(Int,lins)
    cn = 0
    for q ∈ 1:d-1
        for i ∈ 1:Int(Cn(q)), j ∈ 1:Int(Cn(d-q))
            z = GroveBin(PBTree(q,i) + PBTree(d-q,j))
            f = findall(s->s==z.gbin, ins)
            length(f) == 0 ? (cn += 1) : (cc[f] .+= 1)
        end
    end
    #@info "Non-intervals: $cn"
    return cc
end

function between_list_full(a::PBTree,b::PBTree)
    a == b && return true
    h = posetnext_list(a)
    not = true
    for i ∈ 1:length(h)
        if h[i] ≤ b
            not = not & between_list_full(h[i],b)
        else
            not = false
        end
    end
    return not
end

function intervals_full(d::Int)
    g = Array{BigInt,1}()
    f = BitArray{1}()
    for i∈1:Int(Cn(d)),j∈1:Int(Cn(d))
        t = PBTree(d,i) ⊴ PBTree(d,j)
        t ≠ Grove(0) && (push!(g,GroveBin(t).gbin);push!(f,between_list_full(PBTree(d,i),PBTree(d,j))))
    end
    return f[sortperm(g)]
end

function print_interval_bin(d::Int)
    ins = intervals(d)
    for i ∈ ins
        lpad(string(i,base=2),Cn(d),"0") |> println
    end
end

function print_intcomp_bin(d::Int)
    ins = intervals(d)[findall(s->s>0,intcomp(d))]
    inst = intervals(d)[findall(s->s>0,intcompt(d))]
    for i ∈ ins
        i ∉ inst && (lpad(string(i,base=2),Cn(d),"0") |> println)
    end
end

function print_intcompt_bin(d::Int)
    ins = intervals(d)[findall(s->s>0,intcompt(d))]
    for i ∈ ins
        lpad(string(i,base=2),Cn(d),"0") |> println
    end
end
