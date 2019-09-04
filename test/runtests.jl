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
