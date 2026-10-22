### A Pluto.jl notebook ###
# v1.0.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ bdded1f8-1550-427f-8a1b-79cd80cb3b81
begin
    import Pkg
    Pkg.activate(@__DIR__)
end

# ╔═╡ e688d134-14bc-4c66-a2f6-5aebe2b5c33d
using Dynspec

# ╔═╡ f75db5ae-d2ea-4e0d-b64a-7b2a5abb9a55
using Dynspec.Dyno

# ╔═╡ 9caf2ee2-d43d-40c8-853f-da1c9fc483f1
using PlutoUI

# ╔═╡ 4e2e2f91-8740-45d2-9f6e-dcffcae30c36
using PlotlyLight

# ╔═╡ cd236097-9068-4bcf-a85e-fce0f03ee769
using LinearAlgebra

# ╔═╡ 6df13a7f-7efa-44bf-9097-32004b84e9a4
using Random

# ╔═╡ 954c5513-dbe5-469e-a3a5-2484d04bf1e7
md"""
# 📈 Real Business Cycle (RBC) Model
### Model Specification, Analytical Derivatives, First-Order Solution & Interactive Simulation

This notebook demonstrates an end-to-end DSGE workflow in **Dynspec.jl**:
1. **Define** the nonlinear RBC model in Dynspec syntax.
2. **Compute** exact analytical perturbation Jacobians ($\mathbf{A}, \mathbf{B}, \mathbf{C}, \mathbf{D}$) via direct AST ForwardDiff.
3. **Solve** the linear rational expectations system for the recursive decision rule $\hat{y}_t = \mathbf{P}\hat{y}_{t-1} + \mathbf{Q}e_t$.
4. **Simulate** and visualize Impulse Response Functions (IRFs) and stochastic sample paths using **PlotlyLight** (< 0.2s TTFX).
"""

# ╔═╡ 63195518-d459-4b84-a7fa-7d9cd8c8a89f
md"""
---
## 1. Interactive Structural Calibration
Adjust structural parameters in real time to observe their immediate effect on the model solution and dynamics:
"""

# ╔═╡ d379eac9-8539-41a7-82dc-32aef78fd228
@bind alpha Slider(0.20:0.01:0.50; default=0.36, show_value=true)

# ╔═╡ 24cb0f1d-6a7b-431c-b479-25d6a4bc0938
@bind rho Slider(0.50:0.01:0.99; default=0.95, show_value=true)

# ╔═╡ a34c9204-5ad3-4bad-9821-6869e41ff452
@bind delta Slider(0.01:0.005:0.10; default=0.025, show_value=true)

# ╔═╡ 4b54f12a-8f4a-4684-9b99-b7db8a1967a2
@bind r_rate Slider(1.005:0.002:1.05; default=1.01, show_value=true)

# ╔═╡ 7cc45984-c42b-47d6-afd1-b2f4634fc818
@bind theta Slider(1.0:0.05:5.0; default=2.95, show_value=true)

# ╔═╡ 2ddbb302-5aa6-4a3c-a192-4a3b5d89eda0
@bind tau Slider(0.00:0.005:0.10; default=0.025, show_value=true)

# ╔═╡ c262c7af-3b94-4115-96e0-17a1290ac8bb
begin
    beta_val = 1 / r_rate

    model_str = """
@name: RBC_Interactive

# Parameters
alpha <- $alpha
rho   <- $rho
tau   <- $tau
r     <- $r_rate
beta  <- $beta_val
delta <- $delta
psi   <- 0.0
theta <- $theta
phi   <- 0.1

# Initial Steady State
y[~] <- 1.08068253095672
c[~] <- 0.80359242014163
h[~] <- 0.29175631001732
k[~] <- 11.08360443260358
a[~] <- 0.0
b[~] <- 0.0

# Define Exogenous Shocks
e[t] <- N(0, 0.002)
u[t] <- N(0, 0.002)

# Dynamic Equations
c[t]*theta*h[t]^(1+psi) = (1-alpha)*y[t]                                                                         :: "Labor Supply"
k[t] = beta*(((exp(b[t])*c[t])/(exp(b[t+1])*c[t+1]))*(exp(b[t+1])*alpha*y[t+1]+(1-delta)*k[t]))                  :: "Euler Equation"
y[t] = exp(a[t])*(k[t-1]^alpha)*(h[t]^(1-alpha))                                                                  :: "Production"
k[t] = exp(b[t])*(y[t]-c[t])+(1-delta)*k[t-1]                                                                    :: "Capital Accumulation"
a[t] = rho*a[t-1] + tau*b[t-1] + e[t]                                                                            :: "Productivity Shock"
b[t] = tau*a[t-1] + rho*b[t-1] + u[t]                                                                            :: "Preference Shock"
"""
    block = parse_model(model_str)
end

# ╔═╡ 7462972f-5c07-4b3c-ac57-396734388553
begin
    endo_vars = [:a, :b, :y, :k, :c, :h]
    exo_vars  = [:e, :u]
    res, A, B, C, D = Dyno.compute_dsge_jacobians(block; endogenous=endo_vars, exogenous=exo_vars)
    
    md"""
---
## 2. Analytical Derivatives (DSGE Jacobians)
Dynspec differentiates the equations directly on the syntax tree using `ForwardDiff.Dual`:
$$\mathbf{A} \mathbb{E}_t[\hat{y}_{t+1}] + \mathbf{B} \hat{y}_t + \mathbf{C} \hat{y}_{t-1} + \mathbf{D} e_t = 0$$
- **Endogenous Variables**: `$(endo_vars)`
- **Exogenous Shocks**: `$(exo_vars)`
- **$\mathbf{A} = \partial f / \partial y_{t+1}$**: size `$(size(A))`
- **$\mathbf{B} = \partial f / \partial y_t$**: size `$(size(B))`
- **$\mathbf{C} = \partial f / \partial y_{t-1}$**: size `$(size(C))`
- **$\mathbf{D} = \partial f / \partial e_t$**: size `$(size(D))`
- **Steady State Residuals Max Error**: `$(round(maximum(abs.(res)); digits=6))`
"""
end

# ╔═╡ 61b8f5dc-cacf-4ea4-abe0-424650dbb16e
begin
    sol = Dyno.solve_dsge(A, B, C, D)
    P_mat, Q_mat = sol.P, sol.Q
    spec_rad = maximum(abs.(eigvals(P_mat)))
    quad_err = norm(A * P_mat^2 + B * P_mat + C)

    md"""
---
## 3. First-Order DSGE Solution
The stable recursive law of motion is solved via cyclic reduction:
$$\hat{y}_t = \mathbf{P} \hat{y}_{t-1} + \mathbf{Q} e_t$$
- **Spectral Radius $\rho(\mathbf{P})$**: **`$(round(spec_rad; digits=4))`** $(spec_rad < 1.0 ? "✅ Stationary & Stable" : "❌ Explosive")
- **Solution Residual $\|\mathbf{A}\mathbf{P}^2 + \mathbf{B}\mathbf{P} + \mathbf{C}\|$**: `$(round(quad_err; digits=14))`
"""
end

# ╔═╡ 5a621e4d-808e-497c-a1ed-c8d3f4070a8a
md"""
---
## 4. Impulse Response Functions (IRFs)
Choose the shock and response horizon to inspect the transition dynamics:
"""

# ╔═╡ ca056f72-f286-4033-8fcb-538015a3c559
@bind selected_shock Select(["Technology Shock (e)", "Preference Shock (u)"])

# ╔═╡ 0d645113-4c81-420c-b5f3-421a98d7e7ad
@bind horizon Slider(10:5:80; default=40, show_value=true)

# ╔═╡ dc44c81d-5354-443e-aaea-d3877f0bcfce
begin
    s_idx = startswith(selected_shock, "Technology") ? 1 : 2
    Y_irf = Dyno.simulate_linear(P_mat, Q_mat; T=horizon, shock_idx=s_idx, shock_size=1.0)
    
    idx_y = findfirst(==( :y ), endo_vars)
    idx_c = findfirst(==( :c ), endo_vars)
    idx_k = findfirst(==( :k ), endo_vars)
    idx_h = findfirst(==( :h ), endo_vars)
    
    t_axis = 1:horizon
    
    p_irf = PlotlyLight.Plot()
    p_irf.data = [
        Config(x=collect(t_axis), y=Y_irf[idx_y, :], name="Output (y)", mode="lines+markers", line=Config(color="#2563eb", width=3)),
        Config(x=collect(t_axis), y=Y_irf[idx_c, :], name="Consumption (c)", mode="lines+markers", line=Config(color="#10b981", width=3)),
        Config(x=collect(t_axis), y=Y_irf[idx_k, :], name="Capital (k)", mode="lines+markers", line=Config(color="#f59e0b", width=3)),
        Config(x=collect(t_axis), y=Y_irf[idx_h, :], name="Hours (h)", mode="lines+markers", line=Config(color="#8b5cf6", width=3))
    ]
    p_irf.layout.title.text = "Impulse Response to 1-StdDev $(selected_shock) (Horizon = $(horizon) Periods)"
    p_irf.layout.xaxis.title.text = "Quarters (t)"
    p_irf.layout.yaxis.title.text = "Percentage Deviation from Steady State"
    p_irf.layout.hovermode = "x unified"
    p_irf
end

# ╔═╡ aaf13705-2219-4608-9eb1-5d1de4024852
md"""
---
## 5. Stochastic Sample Path Simulation
Generate a synthetic time series of macro aggregates with random Gaussian innovations:
"""

# ╔═╡ 457e3b5c-0e53-426a-8ee3-2258181e80cd
@bind sim_T Slider(50:25:300; default=150, show_value=true)

# ╔═╡ bb64243e-af5c-422e-abc7-d58ab0a35823
@bind sim_seed Slider(1:1:100; default=42, show_value=true)

# ╔═╡ 62194cc6-29ba-4abc-b554-d7b9a7a30e65
begin
    rng = MersenneTwister(sim_seed)
    shocks = randn(rng, length(exo_vars), sim_T) .* 0.01
    
    Y_sim = Dyno.simulate_linear(P_mat, Q_mat; T=sim_T, shocks=shocks)
    
    p_sim = PlotlyLight.Plot()
    p_sim.data = [
        Config(x=collect(1:sim_T), y=Y_sim[idx_y, :], name="Output (y)", mode="lines", line=Config(color="#2563eb", width=2)),
        Config(x=collect(1:sim_T), y=Y_sim[idx_c, :], name="Consumption (c)", mode="lines", line=Config(color="#10b981", width=2)),
        Config(x=collect(1:sim_T), y=Y_sim[idx_k, :], name="Capital (k)", mode="lines", line=Config(color="#f59e0b", width=2))
    ]
    p_sim.layout.title.text = "Stochastic Simulation (T = $(sim_T) Quarters, Seed = $(sim_seed))"
    p_sim.layout.xaxis.title.text = "Time (t)"
    p_sim.layout.yaxis.title.text = "Deviation"
    p_sim.layout.hovermode = "x unified"
    p_sim
end

# ╔═╡ Cell order:
# ╠═954c5513-dbe5-469e-a3a5-2484d04bf1e7
# ╠═bdded1f8-1550-427f-8a1b-79cd80cb3b81
# ╠═e688d134-14bc-4c66-a2f6-5aebe2b5c33d
# ╠═f75db5ae-d2ea-4e0d-b64a-7b2a5abb9a55
# ╠═9caf2ee2-d43d-40c8-853f-da1c9fc483f1
# ╠═4e2e2f91-8740-45d2-9f6e-dcffcae30c36
# ╠═cd236097-9068-4bcf-a85e-fce0f03ee769
# ╠═6df13a7f-7efa-44bf-9097-32004b84e9a4
# ╠═63195518-d459-4b84-a7fa-7d9cd8c8a89f
# ╠═24cb0f1d-6a7b-431c-b479-25d6a4bc0938
# ╠═a34c9204-5ad3-4bad-9821-6869e41ff452
# ╠═4b54f12a-8f4a-4684-9b99-b7db8a1967a2
# ╠═7cc45984-c42b-47d6-afd1-b2f4634fc818
# ╠═2ddbb302-5aa6-4a3c-a192-4a3b5d89eda0
# ╠═c262c7af-3b94-4115-96e0-17a1290ac8bb
# ╠═7462972f-5c07-4b3c-ac57-396734388553
# ╠═61b8f5dc-cacf-4ea4-abe0-424650dbb16e
# ╠═5a621e4d-808e-497c-a1ed-c8d3f4070a8a
# ╠═ca056f72-f286-4033-8fcb-538015a3c559
# ╠═0d645113-4c81-420c-b5f3-421a98d7e7ad
# ╠═d379eac9-8539-41a7-82dc-32aef78fd228
# ╠═dc44c81d-5354-443e-aaea-d3877f0bcfce
# ╠═aaf13705-2219-4608-9eb1-5d1de4024852
# ╠═bb64243e-af5c-422e-abc7-d58ab0a35823
# ╠═457e3b5c-0e53-426a-8ee3-2258181e80cd
# ╠═62194cc6-29ba-4abc-b554-d7b9a7a30e65
