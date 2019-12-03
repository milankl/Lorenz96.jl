[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://img.shields.io/badge/repo_status-active-brightgreen)](https://www.repostatus.org/#active)
[![Travis](https://img.shields.io/travis/com/milankl/Lorenz96.jl?label=Linux%20%26%20osx&logo=travis)](https://travis-ci.com/milankl/Lorenz96.jl)
[![AppVeyor](https://img.shields.io/appveyor/ci/milankl/Lorenz96-jl?label=Windows&logo=appveyor&logoColor=white)](https://ci.appveyor.com/project/milankl/Lorenz96-jl)
[![Cirrus CI](https://img.shields.io/cirrus/github/milankl/Lorenz96.jl?label=FreeBSD&logo=cirrus-ci&logoColor=white)](https://cirrus-ci.com/github/milankl/Lorenz96.jl)
[![DOI](https://zenodo.org/badge/198242642.svg)](https://zenodo.org/badge/latestdoi/198242642)

# Lorenz96.jl - A type-flexible Lorenz 1996 model

![attractor](figs/hovmoeller.png?raw=true "L96 Hovmoeller diagram")

Lorenz96.jl simulates the [Lorenz 96 system](https://en.wikipedia.org/wiki/Lorenz_96_model) with one level (two and three level version planned) for any given number type, as long as conversions (to and from Float64) and arithmetics (+,-,*) are defined - the scaled equations are written division-free. Output always in Float64.

Also supports mixed precision: Different number types can be defined for prognostic variables and calculations on the right-hand side, with automatic conversion on every time step.

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

For mixed precision you also specify a type `Tprog` as the second argument, which is the type used for the prognostic variables
```julia
X = L96(Float16,Float32)
```
Here, the prognostic variables are kept at single-precision (Float32), but calculations on the right-hand side are performed in half-precision (Float16).

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
