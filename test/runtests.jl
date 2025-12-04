using GCATCodes
using Test

@testset "GCATCodes.jl" begin
    @test GCATCodes.demo_function("hello") == "HELLO"
    @test GCATCodes.demo_function("abc") == "ABC"
end
