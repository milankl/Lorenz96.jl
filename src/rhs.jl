function rhs!(dX::Vector{T},X::Vector{T},α::T,F::T,s_inv::T) where {T<:AbstractFloat}

    n = length(dX)
    @boundscheck n == length(X) || throw(BoundsError())

    # periodic boundary conditions
    @inbounds dX[1] = (X[2]-X[n-1])*X[n]*s_inv - α*X[1] + F
    @inbounds dX[2] = (X[3]-X[n])*X[1]*s_inv - α*X[2] + F

    # in the middle of the domain
    for i in 3:n-1
        @inbounds dX[i] = (X[i+1] - X[i-2])*X[i-1]*s_inv - α*X[i] + F
    end

    # periodic boundary conditions
    @inbounds dX[n] = (X[1] - X[n-2])*X[n-1]*s_inv - α*X[n] + F
end
