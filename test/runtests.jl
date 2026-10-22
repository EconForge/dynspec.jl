using Test
using Dynspec
using Dynspec.Dyno
using LinearAlgebra

@testset "Dynspec.jl Tests" begin

    @testset "Types & Distributions" begin
        # Normal 1D
        n1 = Normal(0.016)
        @test n1.Μ == [0.0]
        @test n1.Σ ≈ fill(0.016^2, 1, 1)

        # Normal with mean and std
        n2 = Normal(0.1, 0.05)
        @test n2.Μ == [0.1]
        @test n2.Σ ≈ fill(0.05^2, 1, 1)

        # N() helper
        @test N(0.02).Σ ≈ fill(0.02^2, 1, 1)
        @test N(0.5, 0.02).Μ == [0.5]

        # TimedVar
        v1 = TimedVar(:k, -1)
        @test string(v1) == "k[t-1]"
        v2 = TimedVar(:c, 1)
        @test string(v2) == "c[t+1]"
        v3 = TimedVar(:y, 0)
        @test string(v3) == "y[t]"
        v4 = TimedVar(:z, :~)
        @test string(v4) == "z[~]"
    end

    @testset "Context Structure & Access" begin
        m = parse_model("""
        α <- 0.387
        β <- 0.96
        k[~] <- 5.15
        """)
        ctx = m.context
        @test ctx isa Context
        @test ctx.constants isa Dict{Symbol, Any}
        @test ctx.steady_states isa Dict{Symbol, Any}

        # Explicit dictionary access
        @test ctx.constants[:α] == 0.387
        @test ctx.constants[:β] == 0.96
        @test ctx.steady_states[:k] == 5.15

        # Function taking Context
        f(c) = c.constants[:α] * c.steady_states[:k] + c.constants[:β]
        @test f(ctx) ≈ 0.387 * 5.15 + 0.96
    end

    @testset "Complementarity Conditions" begin
        m1 = parse_model("x[t] + y[t] ⟂ 0 <= x[t] <= inf\n")
        eq1 = m1.equations[1]
        @test eq1.complementarity !== nothing
        @test eq1.complementarity.variable.name == :x
        @test eq1.complementarity.variable.shift == 0

        # One-sided bound: upper bound
        m2 = parse_model("x[t] + y[t] ⟂ x[t] <= h\n")
        eq2 = m2.equations[1]
        @test eq2.complementarity !== nothing
        @test eq2.complementarity.variable.name == :x
        @test eq2.complementarity.lb == -Inf

        # One-sided bound: lower bound
        m3 = parse_model("x[t] + y[t] ⟂ 0 <= x[t]\n")
        eq3 = m3.equations[1]
        @test eq3.complementarity !== nothing
        @test eq3.complementarity.variable.name == :x
        @test eq3.complementarity.ub == Inf

        # Legacy separator |
        m4 = parse_model("x[t] + y[t] | -inf <= x[t]\n")
        eq4 = m4.equations[1]
        @test eq4.complementarity !== nothing
        @test eq4.complementarity.variable.name == :x

        # Ambiguous variable bound should error
        @test_throws ErrorException parse_model("a[t] ⟂ x[t] <= y[t]\n")
    end

    @testset "Neoclassical Model (neo.dyno)" begin
        model_path = joinpath(@__DIR__, "models", "neo.dyno")
        model = read_model(model_path)

        # Context: Metadata
        @test model.context.metadata[:name] == "Neoclassical"
        @test haskey(model.context.metadata, :run)
        runs = model.context.metadata[:run]
        @test any(r -> r[:command] == :steady, runs)
        @test any(r -> r[:command] == :solve, runs)

        # Context: Constants
        @test model.context.constants[:α] ≈ 0.387
        @test model.context.constants[:β] ≈ 0.96
        @test model.context.constants[:γ] ≈ 4.0
        @test model.context.constants[:δ] ≈ 0.1
        @test model.context.constants[:ρ] ≈ 0.9

        # Context: Steady states
        @test model.context.steady_states[:z] ≈ 0.0
        @test model.context.steady_states[:k] ≈ ((1/0.96 - (1 - 0.1))/0.387)^(1/(0.387 - 1))
        @test model.context.steady_states[:y] ≈ model.context.steady_states[:k]^0.387
        @test model.context.steady_states[:i] ≈ 0.1 * model.context.steady_states[:k]
        @test model.context.steady_states[:c] ≈ model.context.steady_states[:y] - model.context.steady_states[:i]

        # Context: Processes
        @test model.context.processes[:e_z] isa Normal
        @test model.context.processes[:e_y] isa Normal

        # Equations
        @test length(model.equations) == 5
        labels = [get(eq.metadata, :label, "") for eq in model.equations]
        @test "Productivity" in labels
        @test "Production" in labels
        @test "Euler Equation" in labels

        # Residuals at steady state
        res = residuals(model)
        @test length(res) == 5
        @test maximum(abs.(res)) < 1e-12

        # Variable listing
        vars = extract_variables(model)
        @test vars[:k] == Set([-1, 0])
        @test vars[:c] == Set([0, 1])
        @test vars[:y] == Set([0, 1])
        @test vars[:z] == Set([-1, 0])

        var_list = list_variables(model)
        @test :k in var_list
        @test :c in var_list
        @test :y in var_list
        @test :z in var_list

        # List constants
        const_list = list_constants(model)
        @test :α in const_list
        @test :β in const_list
        @test :γ in const_list

        # First-order DSGE model NamedTuple
        fo = Dyno.first_order_model(model)
        @test :k in fo.symbols.endogenous
        @test :e_z in fo.symbols.exogenous
        @test :α in fo.symbols.constants
        @test length(fo.residuals) == 5
        @test size(fo.jacobians.A) == (5, 5)
        @test size(fo.jacobians.B) == (5, 5)
        @test size(fo.jacobians.C) == (5, 5)
        @test size(fo.jacobians.D) == (5, 2)

        # Endogenous & Exogenous variable extraction
        vars_ee = endogenous_exogenous(model)
        @test :k in vars_ee.endogenous
        @test :c in vars_ee.endogenous
        @test :e_z in vars_ee.exogenous
        @test :e_y in vars_ee.exogenous

        # Tuple destructuring and alias
        endo_l, exo_l = extract_endogenous_exogenous(model)
        @test endo_l == vars_ee.endogenous
        @test exo_l == vars_ee.exogenous
    end

    @testset "RBC Model with Quantified Assignments (rbc.dyno)" begin
        model_path = joinpath(@__DIR__, "models", "rbc.dyno")
        model = read_model(model_path)

        # Context: Constants
        @test model.context.constants[:alpha] ≈ 0.36
        @test model.context.constants[:beta] ≈ 1 / 1.01

        # Context: Values & Quantified loop
        @test model.context.values[:k][0] ≈ model.context.steady_states[:k] * 1.01
        @test length(model.context.values[:e]) == 10
        @test model.context.values[:e][0] ≈ 0.99
        @test model.context.values[:e][9] ≈ 0.99 / 10
        @test length(model.context.values[:u]) == 9
        @test model.context.values[:u][1] ≈ model.context.values[:e][0] * 0.5

        # Equations
        @test length(model.equations) == 6
        @test maximum(abs.(residuals(model))) < 0.002
    end

    @testset "RBC Dolo Model with Complementarities (rbc_dolo_comp.dyno)" begin
        model_path = joinpath(@__DIR__, "models", "rbc_dolo_comp.dyno")
        model = read_model(model_path)

        # Check block tags
        def_eqs = filter_equations(model, :definitions)
        @test length(def_eqs) == 4

        trans_eqs = filter_equations(model, :transition)
        @test length(trans_eqs) == 2

        arb_eqs = filter_equations(model, :arbitrage)
        @test length(arb_eqs) == 2

        # Check complementarities
        comps = [eq.complementarity for eq in arb_eqs]
        @test all(c !== nothing for c in comps)
        @test [c.variable.name for c in comps] == [:n, :i]

        # Steady state residuals
        @test maximum(abs.(residuals(model))) < 1e-12
    end

    @testset "Dynamic Compilation & Jacobians via ForwardDiff (neo.dyno)" begin
        model_path = joinpath(@__DIR__, "models", "neo.dyno")
        model = read_model(model_path)

        endo = [:z, :k, :y, :i, :c]
        exo = [:e_z, :e_y]
        dm = compile_dynamic(model; endogenous=endo, exogenous=exo)

        @test dm isa DynamicModel
        @test dm.endogenous == endo
        @test dm.exogenous == exo

        # Test steady-state residuals
        res = residuals(dm)
        @test maximum(abs.(res)) < 1e-12

        # Test compute_jacobians at steady state
        r, A, B, C, D = compute_jacobians(dm)
        @test maximum(abs.(r)) < 1e-12
        @test size(A) == (5, 5)
        @test size(B) == (5, 5)
        @test size(C) == (5, 5)
        @test size(D) == (5, 2)

        # Exact match with dyno derivatives:
        # Eq 1: z[t] = ρ*z[t-1] + e_z[t] + e_y[t]
        @test B[1, 1] ≈ -1.0 # ∂/∂z_t
        @test C[1, 1] ≈ 0.9  # ∂/∂z_{t-1} = ρ
        @test D[1, 1] ≈ 1.0  # ∂/∂e_z
        @test D[1, 2] ≈ 1.0  # ∂/∂e_y

        # Eq 2: y[t] = exp(z[t])*k[t-1]^α
        @test B[2, 1] ≈ 1.8859751917924765 # ∂/∂z_t = y_ss
        @test B[2, 3] ≈ -1.0                # ∂/∂y_t
        @test C[2, 2] ≈ 0.14166666666666672 # ∂/∂k_{t-1} = α*y_ss/k_ss

        # Eq 3: k[t] = k[t-1]*(1-δ) + i[t]
        @test B[3, 2] ≈ -1.0 # ∂/∂k_t
        @test B[3, 4] ≈ 1.0  # ∂/∂i_t
        @test C[3, 2] ≈ 0.9  # ∂/∂k_{t-1} = 1 - δ

        # Eq 5 Euler equation derivatives
        @test A[5, 3] ≈ -0.07211123486239622 # ∂/∂y_{t+1}
        @test A[5, 5] ≈ 2.918065509161387   # ∂/∂c_{t+1}
        @test B[5, 2] ≈ 0.026397308196828936 # ∂/∂k_t
        @test B[5, 5] ≈ -2.9180655091613867 # ∂/∂c_t

        # Test direct call from ModelSpec using compute_dsge_jacobians
        r2, A2, B2, C2, D2 = compute_dsge_jacobians(model; endogenous=endo, exogenous=exo)
        @test A2 ≈ A
        @test B2 ≈ B
        @test C2 ≈ C
        @test D2 ≈ D

        # Test compute_jacobians backwards compatibility alias
        r_alias, A_alias, B_alias, C_alias, D_alias = compute_jacobians(model; endogenous=endo, exogenous=exo)
        @test A_alias ≈ A

        # Test in-place zero-allocation views
        r_view, A_view, B_view, C_view, D_view = compute_jacobians!(dm)
        @test A_view ≈ A
        @test B_view ≈ B
        @test C_view ≈ C
        @test D_view ≈ D

        # Test type inference stability
        @test @inferred(residuals(dm)) isa Vector{Float64}
        @test @inferred(compute_dsge_jacobians(dm)) isa Tuple{Vector{Float64}, Vararg{Matrix{Float64}, 4}}
        @test dm.model.context.constants[:α] == 0.387
    end

    @testset "AbstractModelBlock Architecture, AST Retention & Tag Filtering" begin
        # 1. Test AST retention and tag filtering on RBC with dolo complementarities
        rbc_comp_path = joinpath(@__DIR__, "models", "rbc_dolo_comp.dyno")
        block = read_model(rbc_comp_path)

        @test block isa AbstractModelBlock
        @test block isa ModelBlock
        @test block isa ModelSpec
        @test length(block) == 8

        # Verify parsed AST trees are preserved
        @test block.tree !== nothing
        @test block.equations[1].tree !== nothing

        # Iteration
        eq_list = [eq for eq in block]
        @test length(eq_list) == 8
        @test eq_list[1] === block.equations[1]

        # Tag-filtered residuals
        r_trans = residuals(block; filter_by=:transition)
        @test length(r_trans) == 2
        @test maximum(abs.(r_trans)) < 1e-10

        r_defs = residuals(block; filter_by=:definitions)
        @test length(r_defs) == 4
        @test maximum(abs.(r_defs)) < 1e-10

        r_arb = residuals(block; filter_by=:arbitrage)
        @test length(r_arb) == 2

        # Multi-tag and string-based filtering
        r_multi = residuals(block; filter_by=[:transition, :arbitrage])
        @test length(r_multi) == 4

        r_str = residuals(block; filter_by="transition")
        @test length(r_str) == 2

        # Function predicate filtering
        r_pred = residuals(block; filter_by=eq -> :transition in get(eq.metadata, :tags, Symbol[]))
        @test length(r_pred) == 2

        # Residuals with evaluation point overrides
        r_override = residuals(block; filter_by=:transition, at=Dict(:z => 0.5))
        @test length(r_override) == 2

        # equations accessor and collection indexing
        @test length(equations(block)) == 8
        @test block[1] === equations(block)[1]

        # Tag-filtered jacobian on RBC model with variable targets
        J_trans = jacobian(block; wrt=[(:z, 0), (:k, 0)], filter_by=:transition)
        @test size(J_trans) == (2, 2)

        # 2. Generic Differentiation on Neo model
        neo_path = joinpath(@__DIR__, "models", "neo.dyno")
        neo = read_model(neo_path)
        @test neo isa AbstractModelBlock
        @test neo.tree !== nothing

        # Full block jacobian with respect to endogenous variables at t=0
        endo_vars = [(:z, 0), (:k, 0), (:y, 0), (:i, 0), (:c, 0)]
        J_t = jacobian(neo; wrt=endo_vars)
        @test size(J_t) == (5, 5)

        # Parameter sensitivity w.r.t bare symbols (constants)
        J_params = jacobian(neo; wrt=[:α, :β, :δ])
        @test size(J_params) == (5, 3)

        # Jacobian with mixed targets: Symbols (constants) and Tuples / TimedVars (variables)
        wrt_mixed = [:α, (:k, -1), TimedVar(:c, 1)]
        J_mixed = jacobian(neo; wrt=wrt_mixed)
        @test size(J_mixed) == (5, 3)

        # Grouped wrt: returns tuple of matrices
        wrt_grouped = [ [(:z, 0), (:k, 0)], [(:c, 1), (:y, 1)] ]
        J_grp1, J_grp2 = jacobian(neo; wrt=wrt_grouped)
        @test size(J_grp1) == (5, 2)
        @test size(J_grp2) == (5, 2)

        r_grp, (J_g1, J_g2) = residuals_and_jacobian(neo; wrt=wrt_grouped)
        @test length(r_grp) == 5
        @test J_g1 ≈ J_grp1
        @test J_g2 ≈ J_grp2

        # Verify taylor() produces identical order-1 results
        r_taylor, (J_t1, J_t2) = taylor(neo; wrt=wrt_grouped, order=1)
        @test r_taylor == r_grp
        @test J_t1 ≈ J_g1
        @test J_t2 ≈ J_g2
        @test_throws ArgumentError taylor(neo; wrt=wrt_grouped, order=2)

        # Filtered differentiation: differentiate only first 2 equations via predicate
        J_sub = jacobian(neo; wrt=[(:z, 0), (:k, 0)], filter_by=eq -> eq.raw in (neo.equations[1].raw, neo.equations[2].raw))
        @test size(J_sub) == (2, 2)
        @test J_sub[1, 1] ≈ -1.0
        @test J_sub[2, 1] ≈ 1.8859751917924765

        # Simultaneous residuals and jacobian via residuals_and_jacobian() and derivatives alias
        r_all, J_all = residuals_and_jacobian(neo; wrt=endo_vars)
        @test length(r_all) == 5
        @test maximum(abs.(r_all)) < 1e-12
        @test J_all ≈ J_t
        @test derivatives(neo; wrt=endo_vars) == (r_all, J_all)

        # compute_dsge_jacobians on ModelBlock directly
        r_dsge, A_dsge, B_dsge, C_dsge, D_dsge = compute_dsge_jacobians(neo; endogenous=[:z, :k, :y, :i, :c], exogenous=[:e_z, :e_y])
        @test maximum(abs.(r_dsge)) < 1e-12
        @test size(A_dsge) == (5, 5)
        @test size(B_dsge) == (5, 5)
        @test size(C_dsge) == (5, 5)
        @test size(D_dsge) == (5, 2)
        @test B_dsge[1, 1] ≈ -1.0
        @test C_dsge[1, 1] ≈ 0.9
        @test D_dsge[1, 1] ≈ 1.0

        # Dyno.solve_dsge and Dyno.simulate_linear
        sol = Dyno.solve_dsge(neo; endogenous=[:z, :k, :y, :i, :c], exogenous=[:e_z, :e_y])
        @test size(sol.P) == (5, 5)
        @test size(sol.Q) == (5, 2)
        @test maximum(abs.(eigvals(sol.P))) < 1.0 # stable
        @test norm(sol.A * sol.P^2 + sol.B * sol.P + sol.C) < 1e-10

        # Simulate IRF
        Y_irf = Dyno.simulate_linear(sol.P, sol.Q; T=30, shock_idx=1, shock_size=0.01)
        @test size(Y_irf) == (5, 30)
        @test Y_irf[1, 1] ≈ sol.Q[1, 1] * 0.01

        # Type stability checks
        @test @inferred(residuals(neo)) isa Vector{Float64}
        @test @inferred(residuals(neo; filter_by=:transition)) isa Vector{Float64}
        @test @inferred(jacobian(neo; wrt=endo_vars)) isa Matrix{Float64}
        @test @inferred(residuals_and_jacobian(neo; wrt=endo_vars)) isa Tuple{Vector{Float64}, Matrix{Float64}}
        @test @inferred(taylor(neo; wrt=endo_vars)) isa Tuple{Vector{Float64}, Matrix{Float64}}
        @test @inferred(derivatives(neo; wrt=endo_vars)) isa Tuple{Vector{Float64}, Matrix{Float64}}
    end

end
