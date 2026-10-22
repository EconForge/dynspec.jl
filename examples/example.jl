using Dynspec

println("="^65)
println(" Dynspec.jl — Clean Feature Showcase")
println("="^65)

# -------------------------------------------------------------
# 1. Parsing a Model (String DSL) & Inspecting Context
# -------------------------------------------------------------
println("\n--- 1. Model Parsing & Context ---")
model_str = """
@name: NeoclassicalGrowth

# Parameters
α <- 0.36
β <- 0.96
δ <- 0.08
ρ <- 0.90

# Calibrated steady states
z[~] <- 0.0
k[~] <- ((1/β - (1-δ))/α)^(1/(α-1))
y[~] <- k[~]^α
i[~] <- δ*k[~]
c[~] <- y[~] - i[~]

# Exogenous shock
e_z[t] <- N(0, 0.01^2)

# Equations organized by blocks
[transition] :: {
    z[t] = ρ*z[t-1] + e_z[t]           :: "Productivity Process"
    k[t] = (1-δ)*k[t-1] + i[t]         :: "Capital Accumulation"
}

[definitions] :: {
    y[t] = exp(z[t])*k[t-1]^α          :: "Production"
    c[t] = y[t] - i[t]                 :: "Resource Constraint"
}

[arbitrage] :: {
    1 - β*(c[t+1]/c[t])^(-1)*(1-δ + α*y[t+1]/k[t]) ⟂ 0 <= c[t] <= inf :: "Euler Equation"
}
"""

block = parse_model(model_str)
println("Model:            ", block.context.metadata[:name])
println("Parameters:       ", block.context.constants)
println("Steady states:    ", block.context.steady_states)
println("Preserved AST:    ", typeof(block.tree), " (complete grammar tree retained)")

# Dictionary access on context:
println("Parameters:       α = ", block.context.constants[:α], ", steady state k = ", round(block.context.steady_states[:k]; digits=4))

# -------------------------------------------------------------
# 2. Variable Listing
# -------------------------------------------------------------
println("\n--- 2. Variable Occurrences & Time Shifts ---")
println("Variables:        ", list_variables(block))
vars = list_variables(block; shifts=true)
for (v, shifts) in sort(collect(vars), by=x->string(x[1]))
    println("  ", rpad(string(v), 6), " => time shifts: ", sort(collect(shifts)))
end

# List constants
println("Constants:        ", list_constants(block))

# Extract endogenous and exogenous variables via Dyno submodule
endo, exo = Dyno.endogenous_exogenous(block)
println("Endogenous:       ", endo)
println("Exogenous:        ", exo)

# First order model representation
fo = Dyno.first_order_model(block)
println("First order model symbols: ", fo.symbols)

# -------------------------------------------------------------
# 3. Residuals & Tag Filtering (Direct Tree Evaluation)
# -------------------------------------------------------------
println("\n--- 3. Residuals & Tag Filtering ---")
# Residuals for all equations evaluated directly on AST:
r_all = residuals(block)
println("All residuals max abs: ", maximum(abs.(r_all)))

# Residuals filtered by tag:
r_trans = residuals(block; filter_by=:transition)
println("Transition equations residuals (", length(r_trans), " eqs): ", r_trans)

r_defs  = residuals(block; filter_by=:definitions)
println("Definitions equations residuals (", length(r_defs), " eqs):  ", r_defs)

# Residuals evaluated at an off-steady-state point:
r_perturbed = residuals(block; filter_by=:transition, at=Dict(:z => 0.1))
println("Transition residuals at z=0.1: ", r_perturbed)

# -------------------------------------------------------------
# 4. Analytical Jacobians via ForwardDiff
# -------------------------------------------------------------
println("\n--- 4. Analytical Jacobians via ForwardDiff ---")

# Differentiate with respect to dynamic variables (time-dated tuples):
endo_vars = [(:z, 0), (:k, 0), (:y, 0), (:i, 0), (:c, 0)]
J_all = jacobian(block; wrt=endo_vars)
println("Full Jacobian (all eqs wrt $(endo_vars)):")
display(round.(J_all; digits=4))

# Differentiate with respect to parameters (bare symbols):
J_params = jacobian(block; wrt=[:α, :β, :δ])
println("\nParameter sensitivity Jacobian (wrt [:α, :β, :δ]):")
display(round.(J_params; digits=4))

# Differentiate only a specific block:
wrt_trans = [(:z, 0), (:k, -1), (:i, 0)]
J_trans = jacobian(block; wrt=wrt_trans, filter_by=:transition)
println("\nTransition block Jacobian wrt [z[t], k[t-1], i[t]]:")
display(round.(J_trans; digits=4))

# Grouped arguments: returns multiple Jacobians in one pass
wrt_grouped = [ [(:z, 0), (:k, 0)], [(:c, 1), (:y, 1)] ]
J_curr, J_lead = jacobian(block; wrt=wrt_grouped)
println("\nGrouped Jacobian 1 (current [z[t], k[t]]): size ", size(J_curr))
println("Grouped Jacobian 2 (leads [c[t+1], y[t+1]]): size ", size(J_lead))

# Simultaneous residuals and Jacobian evaluation:
r_arb, J_arb = derivatives(block; wrt=[(:c, 0), (:c, 1)], filter_by=:arbitrage)
println("\nArbitrage block residual: ", r_arb[1])
println("Arbitrage Jacobian wrt [c[t], c[t+1]]: ", round.(J_arb; digits=4))

# -------------------------------------------------------------
# 5. Perturbation Matrices (A, B, C, D) for DSGE Solvers
# -------------------------------------------------------------
println("\n--- 5. DSGE Perturbation Jacobians (A, B, C, D) ---")
# A: ∂f / ∂y_{t+1}
# B: ∂f / ∂y_t
# C: ∂f / ∂y_{t-1}
# D: ∂f / ∂e_t
res, A, B, C, D = Dyno.compute_dsge_jacobians(block; endogenous=[:z, :k, :y, :i, :c], exogenous=[:e_z])
println("A (∂f/∂y_{t+1}) size: ", size(A))
println("B (∂f/∂y_t)     size: ", size(B))
println("C (∂f/∂y_{t-1}) size: ", size(C))
println("D (∂f/∂e_t)     size: ", size(D))
println("Residuals max abs:    ", maximum(abs.(res)))
println("\nDone!")
