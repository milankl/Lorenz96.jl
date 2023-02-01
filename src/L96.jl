"""

    X = L96(Float32)

Integrates the Lorenz 1996 model. Standard parameters

    L96(::Type{T}=Float64,              # number format for RHS
    ::Type{Tprog}=T;                    # number format for prognostic variables
    N::Int=10_000,                      # number of time steps
    n::Int=36,                          # number of variables
    X::Array{Float64,1}=zeros(36),      # initial conditions
    α::Real=1.0,                        # friction parameter
    F::Float64=8.0,                     # forcing constant
    s::Float64=1.0,                     # scaling
    η::Float64=0.01,                    # strength of initial perturbation at X[1] if no X provided
    Δt::Float64=0.1,                    # time step
    scheme::String=:RK4                 # time integration scheme, :RK2, :RK3, :RK4 supported
    output::Bool=false                  # return every time step?

# Examples
```jldoc
julia> X1 = L96(Float32);
julia> X2 = L96(Float32,n=6,Δt=0.005);
julia> X3 = L96(Float16,Float32);
```
"""
function L96(::Type{T},                         # number format for RHS
            ::Type{Tprog};                      # number format for prognostic variables
            N::Int = 10_000,                    # number of time steps
            n::Int = 36,                        # number of variables
            X::Vector = [],                     # initial conditions
            α::Real = 1.0,                      # friction parameter
            F::Real = 8.0,                      # forcing constant
            s::Real = 1.0,                      # scaling
            η::Real = 0.01,                     # strength of initial perturbation at X[1] if no X provided
            Δt::Real = 0.1,                     # time step
            scheme::Symbol = :RK4,              # time integration scheme
            output::Bool = false                # return every time step?
            ) where {T<:AbstractFloat,Tprog<:AbstractFloat}

            # if both n and X are specified then they should match!
            if n != 36 && X != [] && length(X) != n
                throw(error("n=$n and length(X)=$(length(X)) do not match."))
            end

            # add the forcing for equilibrium initial conditions
            if X == []      # the default
                X = zeros(n) .+ F
                X[1] += η   # add small perturbation on the first variable
            end

            # check time integration scheme and order and run simulation
            scheme in [:RK2, :RK3, :RK4] || throw(error("Scheme $scheme not implemented."))
            order = parse(Int,string(scheme)[end])
            return RKn(T,Tprog,N,X,α,F,s,Δt,output,RKo=order)
end

"""Calls L96 if no type is specified. Float64 for T and Tprog."""
function L96(;kwargs...)                 # if no type specified
    L96(Float64,Float64;kwargs...)       # use Float64 for everything
end

"""Calls L96 if only one type T is specified. T for RHS and prognostic variables."""
function L96(::Type{T};kwargs...) where {T<:AbstractFloat}
    L96(T,T;kwargs...)
end
