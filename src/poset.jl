# This file is part of Dendriform.jl. It is licensed under the GPL license
# Dendriform Copyright (C) 2017 Michael Reed

export ⋖, ⋗, posetnext, posetprev, between, ⊴, over, under, ↗, ↖

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

Returns PBTRee obtained from a over b operation
"""
over(x::PBTree,y::PBTree) = y.degr > 0 ? over(x,left(y)) ∨ right(y) : x
over(x::AbstractPBTree,y::AbstractPBTree) = over(PBTree(x),PBTree(y))

"""
    ↗(a::AbstractPBTree, b::AbstractPBTree)

Returns PBTRee obtained from a over b operation
"""
↗(x::AbstractPBTree,y::AbstractPBTree) = over(x,y)

"""
    /(a::PBTree, b::PBTree)

Returns PBTRee obtained from a over b operation
"""
/(x::PBTree,y::PBTree) = over(x,y)

"""
    under(a::AbstractPBTree, b::AbstractPBTree)

Returns PBTRee obtained from a under b operation
"""
under(x::PBTree,y::PBTree) = x.degr > 0 ? left(x) ∨ under(right(x),y) : y
under(x::AbstractPBTree,y::AbstractPBTree) = under(PBTree(x),PBTree(y))

"""
    ↖(a::AbstractPBTree, b::AbstractPBTree)

Returns PBTRee obtained from a under b operation
"""
↖(x::AbstractPBTree,y::AbstractPBTree) = under(x,y)

"""
    \(a::PBTree, b::PBTree)

Returns PBTRee obtained from a under b operation
"""
\(x::PBTree,y::PBTree) = under(x,y)
