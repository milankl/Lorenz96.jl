function RK4(T::Type,N::Int,X::Array{Float64,1},F::Float64,s::Float64,Δt::Float64)

    # number of variables
    n = length(X)

    # preallocate for storing results - store without scaling in Float64
    Xout = Array{Float64,2}(undef,n,N+1)
    Xout[:,1] = X

    # scale the initial conditions
    for i = 1:n
        X[i] = s*X[i]
    end

    # Runge Kutta 4th order coefficients including time step and sigma for x
    RKα = [1/6.,1/3.,1/3.,1/6.]*Δt
    RKβ = [1/2.,1/2.,1.]*Δt

    # convert everything to the desired number system determined by T
    X = T.(X)
    F = T.(F*s)
    s_inv = T(1.0 / s)
    RKα = T.(RKα)
    RKβ = T.(RKβ)

    # preallocate memory for intermediate results
    X0 = deepcopy(X)
    X1 = deepcopy(X)
    dX = zero(X)       # tendencies

    for i = 1:N
        @simd for j in 1:n
            @inbounds X1[j] = X[j]
        end

        for rki = 1:4
            rhs!(dX,X1,F,s_inv)

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
