# Lorenz96.jl - A type-flexible Lorenz 1996 model

![attractor](figs/hovmoeller.png?raw=true "L96 Hovmoeller diagram")

Lorenz96.jl simulates the [Lorenz 96 system](https://en.wikipedia.org/wiki/Lorenz_96_model) with one level (two and three level version planned) for any given number type, as long as conversions (to and from Float64) and arithmetics (+,-,*) are defined - the scaled equations are written division-free. Output always in Float64.

# Usage
```julia
using Lorenz96
X = L96()
```
simulates the Lorenz system with `Float64` and standard parameters. Providing the type (which has to be a subtype of `AbstractFloat`), returns the simulation calculated in that type (though output in `Float64`)
```julia
using SoftPosit
X = L96(Posit32)
```
Change parameters by specifying optional arguments
```julia
X = L96(Float32,N=100_000,n=36,X=zeros(36),F=8.0,s=1.0,η=0.01,Δt=0.01,scheme="RK4")
```
with `N` the number of time steps, `n` number of variables, `X` the initial conditions, `F` the forcing constant, `s` a scaling factor of the equations (that will be undone for storage), `η` the perturbation that is added on `X₁` for the default `X`, `Δt` the time step, and `scheme` the time integration scheme.

# Equations

The Lorenz system is scaled with `s` and therefore the prognostic variables are actually  `sX -> X`. The RHS then reads with `s_inv = 1/s`
```
dX_i/dt = (X_i+1 - X_i-2)*X_i-1*s_inv - X_i + F
```

# Installation

In the package manager do
```julia
add https://github.com/milankl/Lorenz96.jl
```
