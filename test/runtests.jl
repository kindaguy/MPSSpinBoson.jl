using MPSSpinBoson
using Test

@testset "MPSSpinBoson.jl" begin
    # Write your tests here.
    @test MPSSpinBoson.test_exsistence() == "MPSSpinBoson package does exist!"
    @test MPSSpinBoson.test_exsistence() != "Hello world!"
end
