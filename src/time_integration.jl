function RK4(   ::Type{T},
                ::Type{Tprog},
                N::Int,
                X::Array{Float64,1},
                F::Float64,
                s::Float64,
                Δt::Float64) where {T<:AbstractFloat,Tprog<:AbstractFloat}

    # number of variables
    n = length(X)

    # preallocate for storing results - store without scaling in Float64
    Xout = Array{Float64,2}(undef,n,N+1)
    Xout[:,1] = X

    # Runge Kutta 4th order coefficients including time step and sigma for x
    RKα = [1/6.,1/3.,1/3.,1/6.]*Δt
    RKβ = [1/2.,1/2.,1.]*Δt

    # convert everything to the desired number system determined by T and scale
    X = Tprog.(X*s)
    F = T.(F*s)
    s_inv = T(1.0 / s)
    RKα = Tprog.(RKα)
    RKβ = Tprog.(RKβ)

    # preallocate memory for intermediate results
    X0 = deepcopy(X)
    X1 = deepcopy(X)                # prognostic variables of type Tprog
    X1rhs = zeros(T,size(X1))       # prognostic variables when converted to type T
    dXrhs = zeros(T,size(X))        # tendencies of type T
    dX = zeros(Tprog,size(X))       # tendencies when converted to type Tprog

    for i = 1:N
        @simd for j in 1:n
            @inbounds X1[j] = X[j]
        end

        for rki = 1:4
            X1rhs = convert(X1rhs,X1)
            rhs!(dXrhs,X1rhs,F,s_inv)
            dX = convert(dX,dXrhs)          # change from T to Tprog

            if rki < 4
                @simd for j in 1:n
                    @inbounds X1[j] = X[j] + dX[j] * RKβ[rki]
                end
            end

            # sum the RK steps on the go
            @simd for j in 1:n
                @inbounds X0[j] += dX[j] * RKα[rki]
            end
        end

        @simd for j in 1:n
            @inbounds X[j] = X0[j]
        end

        # store as 64bit, undo scaling
        Xout[:,i+1] = Float64.(X)/s
    end

    return Xout
end

"""Convert function for two 1-dim arrays, X1, X2, in case their eltypes differ. Convert every element from X1 and store it in X2."""
function Base.convert(X2::Array{T2,1},X1::Array{T1,1}) where {T1<:AbstractFloat,T2<:AbstractFloat}
    m = length(X2)
    @boundscheck m == length(X1) || throw(BoundsError())

    @inbounds for i in 1:m
            X2[i] = T2(X1[i])
    end

    return X2
end

"""Convert function for two 1-dim arrays, X1, X2, in case their eltypes are identical.
Just pass X1, such that X2 is pointed to the same place in memory."""
function Base.convert(X2::Array{T,1},X1::Array{T,1}) where {T<:AbstractFloat}
    m = length(X2)
    @boundscheck m == length(X1) || throw(BoundsError())
    return X1
end
