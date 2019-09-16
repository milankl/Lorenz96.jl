using Lorenz96
using Test

X = L96(N=1000)

@testset "Bounded" begin
    @test maximum(X) < 100.0
    @test minimum(X) > -100.0
end

@testset "Chaotic" begin
    @test all(X[:,end-1] .!= X[:,end])
end

@testset "NoForcing" begin
    X = L96(N=1000,F=0.0)
    @test all(X[2:end,:] .== 0.0)
end

@testset "InitialCond" begin
    ini = randn(36)
    Xf64 = L96(Float64,X=ini,N=10)
    Xf32 = L96(Float32,X=ini,N=10)
    Xf16 = L96(Float16,X=ini,N=10)

    @test Xf64[:,1] == ini
    @test Xf32[:,1] == ini
    @test Xf16[:,1] == ini
end
