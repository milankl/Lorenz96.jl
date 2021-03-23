# Lorenz96.jl - A type-flexible Lorenz 1996 model
[![CI](https://github.com/milankl/Lorenz96.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/milankl/Lorenz96.jl/actions/workflows/CI.yml)
[![DOI](https://zenodo.org/badge/198242642.svg)](https://zenodo.org/badge/latestdoi/198242642)

![attractor](figs/hovmoeller.png?raw=true "L96 Hovmoeller diagram")

Lorenz96.jl simulates the [Lorenz 96 system](https://en.wikipedia.org/wiki/Lorenz_96_model) with one level (two and three level version planned) for any given number type, as long as conversions (to and from Float64) and arithmetics (+,-,*) are defined - the scaled equations are written division-free. Output always in Float64.

Also supports mixed precision: Different number types can be defined for prognostic variables and calculations on the right-hand side, with automatic conversion on every time step.

## Usage
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
with `N` the number of time steps, `n` number of variables, `X` the initial conditions, `F` the forcing constant, `s` a scaling factor of the equations (that will be undone for storage), `η` the perturbation that is added on `X₁` for the default `X`, `Δt` the time step, and `scheme` the time integration scheme. `α` is additional parameter that increases the friction for `α>1`.

For mixed precision you also specify a type `Tprog` as the second argument, which is the type used for the prognostic variables
```julia
X = L96(Float16,Float32)
```
Here, the prognostic variables are kept at single-precision (Float32), but calculations on the right-hand side are performed in half-precision (Float16).

## Equations

The Lorenz system is scaled with `s` and therefore the prognostic variables are actually  `sX -> X`. The RHS then reads with `s_inv = 1/s`
```
dX_i/dt = (X_i+1 - X_i-2)*X_i-1*s_inv - α*X_i + F
```

## Installation

Lorenz96.jl is registered so simply do
```julia
] add Lorenz96
```
where `]` opens the package manager
