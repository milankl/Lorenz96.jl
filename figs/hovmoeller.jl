using Lorenz96
using SoftPosit
using PyPlot

Xp = L96(Posit16,F=8.0,N=1000,n=256,Δt=0.05)
Xf = L96(Float64,F=8.0,N=1000,n=256,Δt=0.05)

##
ioff()

cmin = -10
cmax = 10

fig,(ax1,ax2) = subplots(1,2,sharex=true,sharey=true,figsize=(8,3))

ax1.imshow(Xf,vmin=cmin,vmax=cmax,interpolation="bilinear",aspect="auto",origin="lower")
ax2.pcolormesh(Xp,vmin=cmin,vmax=cmax)

ax1.set_xticks([])
ax1.set_yticks([])

xlim(0,size(Xf)[2])
ylim(0,size(Xf)[1])

ax1.set_title("Float64",loc="left")
ax2.set_title("Posit16",loc="left")

ax1.set_ylabel(L"$X_i$")
ax1.set_xlabel("time")
ax2.set_xlabel("time")

tight_layout()
savefig("figs/hovmoeller.png",dpi=200)
