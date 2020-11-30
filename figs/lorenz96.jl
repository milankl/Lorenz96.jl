using Lorenz96
using LogFixPoint16s
using PyPlot
using JLD2, FileIO
using ColorSchemes
using Statistics, StatsBase
cmap = ColorMap(ColorSchemes.oslo.colors)

# LogFixPoint16s.set_nfrac(10)

N = 300
n = 200
Δt = 0.1
F = 8.0

t = collect(0:N)*Δt
nvec = collect(1:n)

αs = [1.0,1.3,2.0]
# Xf1 = L96(Float32;N,n,Δt,α=αs[1])
# Xf2 = L96(Float32;N,n,Δt,α=αs[2])
# Xf3 = L96(Float32;N,n,Δt,α=αs[3])

X  = fill(F,n)
X[100] = F+0.1

Xa1 = L96(LogFixPoint16,X=X/αs[1];N,n,Δt,α=αs[1])
Xa2 = L96(LogFixPoint16,X=X/αs[2];N,n,Δt,α=αs[2])
Xa3 = L96(LogFixPoint16,X=X/αs[3];N,n,Δt,α=αs[3])

## ENSEMBLE

# parameters
Nens = 100
N2 = 10
α2 = collect(1:1e-3:2.0)
nαs = length(α2)
E = 5.13    # climatological error 

f16,lfx16 = load("N10n200.jld2","f16","lfx16")
bf16,blfx16 = load("BF_N10n200.jld2","f16","lfx16")

# monte carlo sampling for mean
m = Nens
f16m = fill(0.0,nαs,m)
lfx16m = fill(0.0,nαs,m)
bf16m = fill(0.0,nαs,m)
blfx16m = fill(0.0,nαs,m)

for i in 1:m
    # calculate median from a random half of the ensemble members
    f16m[:,i] = mean(f16[:,rand(1:Nens,Nens÷2)],dims=2)
    lfx16m[:,i] = mean(lfx16[:,rand(1:Nens,Nens÷2)],dims=2)
    bf16m[:,i] = mean(bf16[:,rand(1:Nens,Nens÷2)],dims=2)
    blfx16m[:,i] = mean(blfx16[:,rand(1:Nens,Nens÷2)],dims=2)
end

# confidence interval from those
f16lo = [percentile(f16m[i,:],2.5) for i in 1:nαs]
f16hi = [percentile(f16m[i,:],97.5) for i in 1:nαs]

lfx16lo = [percentile(lfx16m[i,:],2.5) for i in 1:nαs]
lfx16hi = [percentile(lfx16m[i,:],97.5) for i in 1:nαs]

bf16lo = [percentile(bf16m[i,:],2.5) for i in 1:nαs]
bf16hi = [percentile(bf16m[i,:],97.5) for i in 1:nαs]

blfx16lo = [percentile(blfx16m[i,:],2.5) for i in 1:nαs]
blfx16hi = [percentile(blfx16m[i,:],97.5) for i in 1:nαs]

## PLOT

ioff()
fig,axs = subplots(2,3,sharey=true,figsize=(10,6))
tight_layout(rect=[0.03,0.04,0.95,0.96],w_pad=0.01)
pos = axs[1,end].get_position()
pos1 = axs[1,1].get_position()
cax = fig.add_axes([pos.x1+0.02,pos.y0,0.015,pos.y1-pos.y0])
bax = fig.add_axes([pos1.x0,0.1,pos.x1-pos1.x0,0.95*(pos.y1-pos.y0)])
axs[2,1].remove()
axs[2,2].remove()
axs[2,3].remove()

vmin,vmax = -8,8
axs[1,1].pcolormesh(t,nvec,Xa1;vmin,vmax,cmap)
axs[1,2].pcolormesh(t,nvec,Xa2;vmin,vmax,cmap)
q = axs[1,3].pcolormesh(t,nvec,Xa3;vmin,vmax,cmap)
cbar = colorbar(q;cax)
cbar.set_ticks([-8,-4,0,4,8])

axs[1,1].set_ylabel("X")
axs[1,1].set_yticks([])

l = ["a","b","c"]
for (ia,ax) in enumerate(axs[1,:])
    ax.set_title("LogFixPoint16, α=$(αs[ia])",loc="left")
    ax.set_title(l[ia],loc="right",fontweight="bold")
    ax.set_xticks([0,10,20,30])
    ax.set_xticklabels([0,10,20,30])
end

axs[1,2].set_xlabel("time [mtu]")
axs[1,3].set_xlabel("time [mtu]")

bax.semilogy(α2,mean(f16,dims=2)/E,"k",lw=2.5,label="Float16")
bax.plot(α2,mean(lfx16,dims=2)/E,"C1",lw=2.5,label="LogFixPoint16")
bax.plot(α2,mean(bf16,dims=2)/E,"grey",lw=2.5,label="BFloat16")
bax.plot(α2,mean(blfx16,dims=2)/E,"C3",lw=2.5,label="BLogFixPoint16")

alfa=0.5
bax.fill_between(α2,f16lo/E,f16hi/E,alpha=alfa,color="k",label="95% confidence")
bax.fill_between(α2,lfx16lo/E,lfx16hi/E,alpha=alfa,color="C1",label="95% confidence")
bax.fill_between(α2,bf16lo/E,bf16hi/E,alpha=alfa,color="grey",label="95% confidence")
bax.fill_between(α2,blfx16lo/E,blfx16hi/E,alpha=alfa,color="C3",label="95% confidence")

bax.legend(loc=1,ncol=2,fontsize=9)
bax.set_ylim(5e-4,6e-2)

bax.set_ylabel("RMSE (normalised)")
bax.set_xlabel("Friction parameter α")
bax.set_xlim(1.0,2.0)

bax.set_title("Error to Float64 after 1 mtu",loc="left")
bax.set_title("d",loc="right",fontweight="bold")

savefig("lorenz96_lfx16.png",dpi=200)
close(fig)