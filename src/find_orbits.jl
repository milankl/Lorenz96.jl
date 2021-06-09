import Distributed: @distributed
import LinearAlgebra: norm

"""Find the Lorenz96 orbit length that is reached starting from X0.
Recursive function that starts testing for orbits over period length l,
but enlarges that testing period if no orbit was found."""
function orbit_length(  X0::Vector{T};          # initial condition
                        lspinup::Int=100_000,   # spinup
                        lchunk::Int=100_000,    # chunk size
                        verbose::Bool=true,     # feedback?
                        nchunkmax::Int=10_000,  # abort after this many chunks
                        ) where T      
    
    X0 = L96(T,X=X0,n=length(X0),N=lspinup)

    # use several X for testing
    #   1: initial conditions, static
    # dynamic X that move forward in time (in case 1 isn't part of the orbit yet)
    #   2: sqrt(nchunk)
    #   3: log1.5(nchunk)
    #   4: cbrt(nchunk)
    #   5: log2(nchunk)
    log1p5(x) = log(x)/log(1.5)
    tests = [sqrt,log1p5,cbrt,log2]
    ntests = length(tests)+1

    Xtest = Array{T,2}(undef,length(X0),ntests)
    for i in 1:ntests
        Xtest[:,i] = X0     # initialise all with X0
    end

    # respective time stamps
    ltest = zeros(Int,ntests)

    l = 0               # orbit length (0 = not found yet)
    nchunk = 0          # count the chunks

    while l == 0 && nchunk < nchunkmax
        # compute next chunk of L96 trajectories
        X = L96(T,X=X0,n=length(X0),N=lchunk,output=true)
        X0 = X[:,end]   # new initial condition for next chunk

        j = 2   # index within chunk, skip initial condition
        @inbounds while l == 0 && j < lchunk  # stop when orbit is found (l>0)
            # check for periodicity   
            Xj = X[:,j] 
            for i in 1:ntests
                l = Xtest[:,i] == Xj ? nchunk*lchunk+j-ltest[i]-1 : 0
            end
            j += 1
        end

        nchunk += 1

        # move Xtest[2:end] slowly forwards in time
        for (i,func) in enumerate(tests)
            if func(nchunk) % 1 == 0.0
                Xtest[:,i+1] = copy(X0)
                ltest[i+1] = nchunk*lchunk
            end
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
                                nchunkmax::Int=10_000,
                                verbose::Bool=true) where T    # chunk size
    len,X = orbit_length(X0;lspinup,lchunk,nchunkmax,verbose)    
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
                        nchunkmax::Int=10_000,      # maximum number of chunks to abort
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
        orbit = orbit_length_minimum(X;lspinup,lchunk,nchunkmax,verbose)
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