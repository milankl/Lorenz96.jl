using Lorenz96
using LogFixPoint16s
using PyPlot
using Statistics, StatsBase
using BFloat16s
# using ColorSchemes
# cmap = ColorMap(ColorSchemes.oslo.colors)

LogFixPoint16s.set_nfrac(7)

# parameters
Nens = 100
N = 10
n = 200
Δt = 0.1
F = 8.0

# coordinates
t = collect(0:N)*Δt
nvec = collect(1:n)

# friction parameter space
αs = collect(1:1e-3:2.0)
nαs = length(αs)

# # preallocate
# f16 = fill(0.0,nαs,Nens)
# lfx16 = fill(0.0,nαs,Nens)

# for (iα,α) = enumerate(αs)

#     # control run
#     X0 = L96(;n,Δt,F,α)[:,1000:end]

#     for ie in 1:Nens
#         # pick random start date from control run
#         X = X0[:,rand(1:size(X0)[2])]

#         # experiments
#         XF64 = L96(Float64;X,N,n,Δt,α,F)[:,end]
#         XF16 = L96(BFloat16;X,N,n,Δt,α,F)[:,end]
#         XLFX16 = L96(LogFixPoint16;X,N,n,Δt,α,F)[:,end]

#         f16[iα,ie] = sqrt(mean((XF64-XF16).^2))
#         lfx16[iα,ie] = sqrt(mean((XF64-XLFX16).^2))
#     end
# end

# E = 5.13    # climatological error 

# # monte carlo sampling for mean
# m = Nens
# f16m = fill(0.0,nαs,m)
# lfx16m = fill(0.0,nαs,m)

# for i in 1:m
#     # calculate median from a random half of the ensemble members
#     f16m[:,i] = mean(f16[:,rand(1:Nens,Nens÷2)],dims=2)
#     lfx16m[:,i] = mean(lfx16[:,rand(1:Nens,Nens÷2)],dims=2)
# end

# confidence interval from those
f16lo = [percentile(f16m[i,:],5) for i in 1:nαs]
f16hi = [percentile(f16m[i,:],95) for i in 1:nαs]

lfx16lo = [percentile(lfx16m[i,:],5) for i in 1:nαs]
lfx16hi = [percentile(lfx16m[i,:],95) for i in 1:nαs]

fig,ax = subplots(figsize=(10,3))

ax.semilogy(αs,mean(f16,dims=2)/E,label="Float16")
ax.plot(αs,mean(lfx16,dims=2)/E,label="LogFixPoint16")

ax.fill_between(αs,f16lo/E,f16hi/E,alpha=0.3)
ax.fill_between(αs,lfx16lo/E,lfx16hi/E,alpha=0.3)

ax.legend()

ax.set_ylabel("Error after $(N÷10) mtu [%]")
ax.set_xlabel("Friction parameter α")
ax.set_xlim(1.0,2.0)
tight_layout()