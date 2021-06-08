"""Orbit defined by its minimum (L2-norm) and """
mutable struct Orbit{T<:AbstractFloat}
    N::Int              # number of variables
    length::Int         # length of orbit
    min::Vector{T}      # identify orbit by the point with minimum L2 norm
    norm::T             # that L2 norm
    basin::Int          # basin of attractions, number of random initial conditions that yield that orbit
    basin_norm::Real    # normalised to a fraction
end

# genereator function
Orbit(l::Int,min::Vector{T}) where T = Orbit{T}(length(min),l,min,norm(min),1,0.0)

# two orbits are equivalent when their βs, length and min agree
function Base.:(==)(o1::Orbit{T},o2::Orbit{T}) where T
    # if the norms don't match the orbits can't be identical
    # if the lengths match then the mins still might disagree
    condition1 = (o1.norm == o2.norm) && (o1.length == o2.length)

    if condition1
        o2min = copy(o2.min)
        not_equiv = o1.min != o2min
        nshifts = length(o2min)-1
        ishift = 0
        while not_equiv && ishift < nshifts
            o2min = circshift(o2min,1)
            not_equiv = o1.min != o2min ? true : false
            ishift += 1
        end
        return ~not_equiv
    else
        return false
    end
end

# two equivalent orbits added adds their basins
function Base.:(+)(o1::Orbit{T},o2::Orbit{T}) where T
    @assert o1 == o2 "$o1 != $o2."
    o1.basin += o2.basin
    return o1
end

# addition without uniqueness check
function unsafe_add(o1::Orbit{T},o2::Orbit{T}) where T
    o1.basin += o2.basin
    return o1
end

# return index of first element in x equal to y
function findin(x::Array{T,1},y::T) where T
    @inbounds for i in 1:length(x)
        x[i] == y && return i
    end
end

# orbit array reduction methods for distributed calculation
# check for uniqueness and extend orbit array if so
reduce_orbits(o1::Orbit{T},o2::Orbit{T}) where T = o1 == o2 ? [unsafe_add(o1,o2)] : [o1,o2]

function reduce_orbits(orbit::Orbit{T},orbits::Array{Orbit{T},1}) where T
    if orbit in orbits  # check whether orbit already exists in orbits
        # if so add their basins of attraction
        orbits[findin(orbits,orbit)] += orbit
        return orbits
    else  # append the orbit to the orbits array
        return push!(orbits,orbit)
    end
end

function reduce_orbits(orbits::Array{Orbit{T},1},orbit::Orbit{T}) where T
    if orbit in orbits  # check whether orbit already exists in orbits
        # if so add their basins of attraction
        orbits[findin(orbits,orbit)] += orbit
        return orbits
    else  # append the orbit to the orbits array
        return push!(orbits,orbit)
    end
end

function reduce_orbits(orbits1::Array{Orbit{T},1},orbits2::Array{Orbit{T},1}) where T
    for orbit in orbits1
        orbits2 = reduce_orbits(orbit,orbits2)
    end
    return orbits2
end

function normalise_basins!(orbits::Array{Orbit{T},1}) where T
    s = 0
    for o in orbits
        s += o.basin
    end
    for o in orbits
        o.basin_norm = o.basin/s
    end
end
     
# define < for sorting 
function Base.isless(o1::Orbit{T},o2::Orbit{T}) where T
    o1.length < o2.length
end

function Base.show(io::IO, o::Orbit{T}) where T

    N = o.N
    l = @sprintf("%10d",o.length)
    m = @sprintf("%s",repr(o.min))
    b = @sprintf("%s",repr(o.basin_norm))

    print(io,"Orbit{$T,N=$N}(length=$l, min=$m, basin=$b)")
end