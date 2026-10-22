# DSGE Models & Perturbation (`Dyno` Submodule)

The `Dyno` submodule provides domain-specific tools for solving and simulating dynamic stochastic general equilibrium (DSGE) models.

While `Dynspec` core is a general dynamic model specification and differentiation engine, `Dyno` specializes in first-order perturbation approximations and recursive linear state-space representations.

```julia
using Dynspec
using Dynspec.Dyno
```

---

## 1. Endogenous and Exogenous Variables

In DSGE models, variables are traditionally partitioned into endogenous variables ($y_t$) and exogenous stochastic disturbances ($e_t$).

### `Dyno.endogenous_exogenous`
Extract the ordered lists of endogenous and exogenous symbols:

```julia
vars = Dyno.endogenous_exogenous(block)
vars.endogenous   # [:z, :k, :y, :i, :c]
vars.exogenous    # [:e_z]

# Or via tuple destructuring:
endo, exo = Dyno.endogenous_exogenous(block)
```

By default, exogenous variables are identified from stochastic processes defined in the context (`block.context.processes`), and endogenous variables are identified from the remaining equation variables and calibrated steady states. You can also explicitly specify custom orderings:

```julia
vars = Dyno.endogenous_exogenous(block;
    endogenous = [:k, :c, :y, :i, :z],
    exogenous  = [:e_z]
)
```

---

## 2. First-Order Model Representation (`first_order_model`)

The `first_order_model` function packages the complete first-order perturbation description of a DSGE model into a unified, transparent `NamedTuple`.

### System Definition
A first-order DSGE model is represented by the equilibrium condition:
```math
\mathbb{E}_t \left[ f(y_{t+1}, y_t, y_{t-1}, e_t; \theta) \right] = 0
```
Linearized around the steady state, this becomes:
```math
A \cdot \hat{y}_{t+1} + B \cdot \hat{y}_t + C \cdot \hat{y}_{t-1} + D \cdot e_t = 0
```
where:
- $A = \frac{\partial f}{\partial y_{t+1}}$ (lead Jacobian)
- $B = \frac{\partial f}{\partial y_t}$ (current Jacobian)
- $C = \frac{\partial f}{\partial y_{t-1}}$ (lag Jacobian)
- $D = \frac{\partial f}{\partial e_t}$ (shock Jacobian)

### Calling `first_order_model`
```julia
model = read_model("test/models/neo.dyno")
fo = Dyno.first_order_model(model)
```

### Returned Structure
`first_order_model` returns a structured `NamedTuple`:
```julia
(
    symbols = (
        endogenous = [:z, :k, :y, :i, :c],
        exogenous  = [:e_z, :e_y],
        constants  = [:α, :β, :γ, :δ, :ρ]
    ),
    residuals = [0.0, 0.0, 0.0, 0.0, 0.0],
    jacobians = (
        A = Matrix{Float64}(...),  # ∂f / ∂y_{t+1}
        B = Matrix{Float64}(...),  # ∂f / ∂y_t
        C = Matrix{Float64}(...),  # ∂f / ∂y_{t-1}
        D = Matrix{Float64}(...)   # ∂f / ∂e_t
    )
)
```

Access components directly with dot syntax:
```julia
fo.symbols.endogenous   # [:z, :k, :y, :i, :c]
fo.jacobians.A          # 5×5 Matrix{Float64}
fo.residuals            # Vector{Float64} of steady-state residuals
```

---

## 3. Direct Perturbation Jacobians (`compute_dsge_jacobians`)

If you only need the matrices $(r, A, B, C, D)$ as a flat tuple:

```julia
r, A, B, C, D = Dyno.compute_dsge_jacobians(block)
```

---

## 4. Solving the Model (`solve_dsge`)

`Dyno.solve_dsge` solves for the stable recursive law of motion using cyclic reduction:
```math
y_t = P \cdot y_{t-1} + Q \cdot e_t
```

You can pass either the `first_order_model` NamedTuple or the `ModelBlock` directly:

```julia
# Solve from first_order_model:
fo = Dyno.first_order_model(model)
sol = Dyno.solve_dsge(fo)

# Or solve directly from block:
sol = Dyno.solve_dsge(model)

sol.P   # 5×5 state transition matrix
sol.Q   # 5×2 shock transition matrix
```

---

## 5. Linear Simulation and Impulse Responses (`simulate_linear`)

Compute impulse response functions (IRFs) or simulate trajectories under shock sequences:

```julia
# 1. One-time unit shock on the first shock (IRF):
irf = Dyno.simulate_linear(sol.P, sol.Q; T=40, shock_idx=1, shock_size=1.0)
# Returns a 5×40 matrix of simulated variable deviations

# 2. Simulation with an explicit shock matrix (n_shocks × T):
shocks = randn(length(sol.exogenous), 100) * 0.01
sim = Dyno.simulate_linear(sol.P, sol.Q; T=100, shocks=shocks)
```
