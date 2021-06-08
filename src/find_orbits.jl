import Distributed: @distributed
import LinearAlgebra: norm

"""Find the Lorenz96 orbit length that is reached starting from X0.
Recursive function that starts testing for orbits over period length l,
but enlarges that testing period if no orbit was found."""
function orbit_length(  X0::Vector{T},      # initial condition
                        l::Int) where T     # max period length
    
    X = L96(T,X=X0,n=length(X0),N=l-1,output=true)
    Xend = X[:,end]

    n = 0                           # orbit length (0 = not found yet)
    j = l                           # respective index of X

    while n == 0 && j > 1               # stop when orbit is found (n>0) or when max period is reached
        j -= 1                          # walk backwards through the array to minimize spin-up confliction
        n = Xend == X[:,j] ? l-j : 0    # check for periodicity
    end

    # recursive call on orbit_length
    if n == 0         # in case no orbit was found in l steps, repeat with 5l steps
        n,Xend = orbit_length(Xend,5l)
    else
        return n,Xend    # return orbit length and an X that's on the orbit
    end
end

"""Returns the minimum X on the orbit of length l, starting from X0 on the orbit."""
function orbit_minimum( X0::Vector{T},      # initial condition on the orbit
                        l::Int) where T     # the period length of the orbit

    X = L96(T,X=X0,n=length(X0),N=l-1,output=true)
    Xmin = copy(X0)
    Xnorm = norm(Xmin)

    # find the L2-norm minimum by running through X
    for i in 2:l
        normXi = norm(X[:,i])
        if normXi < Xnorm
            Xmin = X[:,i]
            Xnorm = normXi
        end
    end
    return Xmin
end

"""Returns both the length and minimum of the orbit reached from initial condition X0."""
function orbit_length_minimum(  X0::Vector{T},      # initial condition on the orbit
                                l::Int=100) where T # initial period length to test for
    o,x = orbit_length(X0,l)    
    xmin = orbit_minimum(x,o)
    return Orbit(o,xmin)
end

"""For a given number format T, number of variables N and n initial conditions, find orbits of
Lorenz96 and test for their uniqueness."""
function find_orbits(   ::Type{T},                  # Number format
                        N::Int,                     # number of variables in L96
                        n::Int=100) where T         # n initial conditions

    # pre-allocate empty array of orbits
    orbits = Orbit[]

    # create initial conditions, discard spinup
    spinup = 100
    ini = L96(Float64,n=N,N=10*n+spinup,output=true,η=0.005+0.01*rand())
    ini = ini[:,unique(rand(spinup:10*n+spinup,10*n))[1:n]]

    tic = time()
    orbits = @distributed (reduce_orbits) for i in 1:n            # for n ICs calculate orbit lengths & x
        X = convert(Vector{T},ini[:,i])
        orbit_length_minimum(X)    
    end

    # in case n=1 pack into a Vector
    orbits = isa(orbits,Orbit) ? [orbits] : orbits   
    
    sort!(orbits)
    normalise_basins!(orbits)
    
    toc = time()
    time_elapsed = readable_secs(toc-tic)
    println("Found $(length(orbits)) orbits in $time_elapsed.")
    return orbits
end