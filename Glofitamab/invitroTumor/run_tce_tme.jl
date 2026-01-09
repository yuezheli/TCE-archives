# purpose of this code: run in vitro TCE model (ModelTrace version) for glofit example

using ModelTrace # This is the julia package that is being developed
using Unitful
using Plots
using ModelingToolkit
using DifferentialEquations
using DataFrames
using ProjectRoot
using MRGEvents
using CSV
using Optim
using OptimizationOptimJL

include("tce_tme.jl");

pobj = ModelTrace.load_key(@projectroot("keys.yml"),variant = "human");

#### THIS IS IMPORTANT
### Model Name must match file name.jl
ModelTrace.@register tce_tme = tceTMEModel(pobj;name=:tce_tme);
#### THIS IS IMPORTANT

tce_tme = structural_simplify(tce_tme);

Tumor_cells_init = 1000 # placeholder 
ETratio = 5
T_cells_init = Tumor_cells_init * ETratio
test_dose = 1.0 # [nM]

# Define simulation time span
const s_per_hr = 60
final_time_hr = 24.0; # [hr]
tspan_s = (0.0, final_time_hr*s_per_hr);  # [s]
prob0 = ODEProblem(tce_tme, [], tspan_s, [tce_tme.TCE => test_dose])

sol0 = solve(prob0, alg=AutoTsit5(Rosenbrock23()), saveat = s_per_hr, abstol=1e-6, reltol=1e-3);

#prob_glo = ODEProblem(tce_tme, [], tspan_s, [tce_tme.TCE => test_dose, tce_tme.Tumor_cells_init => Tumor_cells_init, tce_tme.total_T_cells_init => T_cells_init])
prob_glo = ModelTrace.update(prob0, p = [tce_tme.Tumor_cells_init => Tumor_cells_init, tce_tme.total_T_cells_init => T_cells_init]); 

solglo = solve(prob_glo, alg=AutoTsit5(Rosenbrock23()), saveat = s_per_hr, abstol=1e-6, reltol=1e-3);
