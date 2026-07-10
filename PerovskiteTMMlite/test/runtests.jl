# test/runtests.jl
# Numerical smoke tests for the TMM core. No plotting deps — runs under
#   julia --project=.
#   julia> ] test
# or:  julia --project=. -e 'using Pkg; Pkg.test()'

using PerovskiteTMMlite
using Test

@testset "PerovskiteTMMlite TMM core" begin

    λ    = collect(300.0:10.0:1000.0)
    x_nm = collect(0.0:10.0:1000.0)

    @testset "single-wavelength coh_tmm" begin
        stack  = build_stack1()
        n_list = layer_nk_at(stack, 500.0)
        coh    = coh_tmm("s", n_list, stack.d_nm, 0.0, 500.0)
        @test 0.0 <= coh["R"] <= 1.0
        @test 0.0 <= coh["T"] <= 1.0
        A = absorp_in_each_layer(coh)
        @test length(A) == length(stack.d_nm)
        @test all(a -> a > -1e-9, A)          # no (meaningfully) negative absorption
        @test isapprox(sum(A), 1.0; atol = 1e-3)  # energy conservation
    end

    @testset "Stack 1 spectra" begin
        stack = build_stack1()
        res   = compute_G_matrix(stack, λ, x_nm)
        A_tot = vec(sum(res.A_per_layer[2:end-1, :], dims = 1))
        closure = res.R .+ A_tot .+ res.T
        @test all(c -> isapprox(c, 1.0; atol = 5e-3), closure)   # R+ΣA+T ≈ 1
        A_pvk = res.A_per_layer[stack.active_layer, :]
        tol = 1e-6
        @test all(-tol .<= A_pvk .<= 1.0+tol)
        @test maximum(A_pvk) > 0.5              # perovskite absorbs strongly
    end

    @testset "generation decays from the front" begin
        stack     = build_stack1()
        starts    = layer_starts_nm(stack.d_nm)
        ℓ         = stack.active_layer
        xg        = collect(starts[ℓ]:5.0:(starts[ℓ] + stack.d_nm[ℓ]))
        res       = compute_G_matrix(stack, [500.0], xg)
        g         = res.G[:, 1]
        @test g[1] > g[end]                     # blue absorbed near the front
        @test all(g .>= 0.0)
    end

    @testset "Ag back reflector is opaque" begin
        stack = build_stack1_Ag()
        res   = compute_G_matrix(stack, λ, x_nm)
        @test maximum(res.T) < 0.02             # ~no light passes through Ag
    end

    @testset "all three stacks build & run" begin
        for build in (build_stack1, build_stack2, build_stack1_Ag)
            stack = build()
            @test stack.active_layer == 6
            res = compute_G_matrix(stack, [400.0, 600.0], x_nm)
            @test size(res.A_per_layer, 2) == 2
        end
    end
end
