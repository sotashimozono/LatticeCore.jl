using LatticeCore
using Test
using Aqua

@testset "Aqua" begin
    # Run all Aqua checks. We start with `test_all` and only carve out
    # checks individually if a false positive shows up.
    Aqua.test_all(LatticeCore; ambiguities = false, persistent_tasks = false)

    # `test_all` skips ambiguity checks above because Base/stdlib can
    # introduce method ambiguities outside our control. Restrict the
    # ambiguity check to our own package so regressions inside
    # LatticeCore are still caught.
    Aqua.test_ambiguities(LatticeCore)
end
