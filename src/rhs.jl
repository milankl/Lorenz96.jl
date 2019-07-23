function rhs!(dX::Array{T,1},X::Array{T,1},F::T,s_inv::T) where {T<:AbstractFloat}

    n = length(dX)
    @boundscheck n == length(X) || throw(BoundsError())

    @inbounds dX[1] = (X[2]-X[n-1])*X[n]*s_inv - X[1] + F
    @inbounds dX[2] = (X[3]-X[n])*X[1]*s_inv - X[2] + F

    @simd for i in 3:n-1
        @inbounds dX[i] = (X[i+1] - X[i-2])*X[i-1]*s_inv - X[i] + F
    end

    @inbounds dX[n] = (X[1] - X[n-2])*X[n-1]*s_inv - X[n] + F
end
