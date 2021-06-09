import Distributed: @distributed
import LinearAlgebra: norm

"""Find the Lorenz96 orbit length that is reached starting from X0.
Recursive function that starts testing for orbits over period length l,
but enlarges that testing period if no orbit was found."""
function orbit_length(  X0::Vector{T};          # initial condition
                        lspinup::Int=100_000,   # spinup
                        lchunk::Int=100_000,    # chunk size
                        verbose::Bool=true
                        ) where T      
    
    X0 = L96(T,X=X0,n=length(X0),N=lspinup)
    
    Xtest1 = copy(X0)   # use X after spinup to test for periodicity
    Xtest2 = copy(X0)   # second X for testing will slowly move forwards in time
                        # in case Xtest1 isn't part of the orbit yet
    ltest2 = 0          # time stamp for Xtest2

    l = 0               # orbit length (0 = not found yet)
    nchunk = 0          # count the chunks

    while l == 0
        # compute next chunk of L96 trajectories
        X = L96(T,X=X0,n=length(X0),N=lchunk,output=true)
        X0 = X[:,end]   # new initial condition for next chunk

        j = 2   # index within chunk, skip initial condition
        @inbounds while l == 0 && j < lchunk  # stop when orbit is found (l>0)
            # check for periodicity   
            Xj = X[:,j] 
            l = Xtest1 == Xj ? nchunk*lchunk+j-1 : 0
            l = Xtest2 == Xj ? nchunk*lchunk+j-ltest2-1 : 0
            j += 1
        end

        nchunk += 1

        # move Xtest2 slower than chunks forward in time 
        if log2(nchunk) % 1 == 0.0
            Xtest2 = copy(X0)
            ltest2 = nchunk*lchunk
        end
    end

    if verbose
        println("Orbit of length=$l found after $nchunk chunks of size $lchunk.")
    end

    return l,X0    # return orbit length and an X that's on the orbit
end

"""Returns the minimum X on the orbit of length l, starting from X0 on the orbit."""
function orbit_minimum( X0::Vector{T},               # initial condition on the orbit
                        l::Int;                      # the period length of the orbit
                        lchunk::Int=100_000) where T     # chunk size to reduce memory stress
    
    Xmin = copy(X0)         # start with X0 as minimum
    Xnorm = norm(Xmin)

    ll = 0                  # counter to count up till l for each chunk
    while ll < l 
        X = L96(T,X=X0,n=length(X0),N=lchunk,output=true)
        X0 = X[:,end]       # last time step are the new initial conditions

        # find the L2-norm minimum by running through X
        @inbounds for i in 2:lchunk
            normXi = norm(X[:,i])
            if normXi < Xnorm
                Xmin = X[:,i]
                Xnorm = normXi
            end
        end

        ll += lchunk
    end

    return Xmin
end

"""Returns both the length and minimum of the orbit reached from initial condition X0."""
function orbit_length_minimum(  X0::Vector{T};              # initial condition on the orbit
                                lspinup::Int=100_000,       # spin-up length
                                lchunk::Int=100_000,
                                verbose::Bool=true) where T    # chunk size
    len,X = orbit_length(X0;lspinup,lchunk,verbose)    
    Xmin = orbit_minimum(X,len;lchunk)
    return Orbit(len,Xmin)
end

"""For a given number format T, number of variables N and n initial conditions, find orbits of
Lorenz96 and test for their uniqueness."""
function find_orbits(   ::Type{T},                  # Number format
                        N::Int,                     # number of variables in L96
                        n::Int=100;                 # number of initial conditions
                        lini::Int=100_000,          # initial spin-up length in Float64
                        lspinup::Int=100_000,       # spinup with format T
                        lchunk::Int=100_000,        # chunk size to reduce memory stress
                        verbose::Bool=true          # report every orbit found?
                        ) where T         # n initial conditions

    # pre-allocate empty array of orbits
    orbits = Orbit[]

    # create initial conditions, discard spinup
    lini = max(lini,10*n)       # increase lini in case of many initial conditions
    spinup = 100
    ini = L96(Float64,n=N,N=lini+spinup,output=true,η=0.005+0.01*rand())
    ini = ini[:,unique(rand(spinup:lini+spinup,10*n))[1:n]]

    tic = time()
    orbits = @distributed (reduce_orbits) for i in 1:n            # for n ICs calculate orbit lengths & x
        X = convert(Vector{T},ini[:,i])
        orbit = orbit_length_minimum(X;lspinup,lchunk,verbose)
        if verbose
            println(orbit)
        end
        orbit
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