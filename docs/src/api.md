# API Reference

```@meta
CurrentModule = Dynspec
```

## Model Construction & Types

```@docs
AbstractModelBlock
ModelBlock
Equation
Context
Normal
TimedVar
Complementarity
equations
```

## Model Ingestion

```@docs
parse_model
read_model
```

## Analysis & Evaluation

```@docs
list_variables
residuals
jacobian
derivatives
filter_equations
```

## Dyno Submodule

```@docs
Dyno
Dyno.model_properties
Dyno.endogenous_exogenous
Dyno.compute_dsge_jacobians
Dyno.compile_dynamic
```
