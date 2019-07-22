function L96(::Type{T}=Float64;
            N::Int=10_000,
            n::Int=36,
            X::Array{Float64,1}=zeros(36),
            F::Float64=8.0,
            s::Float64=1.0,
            η::Float64=0.01,
            Δt::Float64=0.01,
            scheme::String="RK4") where {T<:AbstractFloat}

            # if both n and X are specified then they should match!
            if n != 36 && X != zeros(36) && length(X) != n
                throw(error("n=$n and length(X)=$(length(X)) do not match."))
            end

            # add the forcing for equilibrium initial conditions
            if X == zeros(36)   # the default
                X = zeros(n) .+ F
                X[1] += η # add small perturbation on the first variable
            end

            if scheme == "RK4"
                return RK4(T,N,X,F,s,Δt)
            else
                throw(error("Other schemes than RK4 not implemented yet."))
            end
end

X₁
