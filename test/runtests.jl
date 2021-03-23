using Lorenz96
using Test

@testset "Bounded" begin
    X = L96(N=1000)
    @test maximum(X) < 100.0
    @test minimum(X) > -100.0
end

@testset "Chaotic" begin
    X = L96(N=1000,output=true)
    @test all(X[:,end-1] .!= X[:,end])
end

@testset "NoForcing" begin
    X = L96(N=1000,F=0.0,output=true)
    @test all(X[2:end,:] .== 0.0)
end

@testset "InitialCond" begin
    ini = randn(36)
    Xf64 = L96(Float64,X=ini,N=10,output=true)
    Xf32 = L96(Float32,X=ini,N=10,output=true)
    Xf16 = L96(Float16,X=ini,N=10,output=true)

    @test Xf64[:,1] == ini
    @test Xf32[:,1] == ini
    @test Xf16[:,1] == ini
end

@testset "Mixed Precision" begin
    X1 = L96(Float32,Float32,N=1000)
    X2 = L96(Float32,N=1000)

    @test X1 == X2

    X3 = L96(Float16,Float32,N=1000)
    @test X1 != X3
end
