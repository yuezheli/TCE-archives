# date: 10/24/24
# author: Yuezhe Li 
# purpose of this code: run in vitro cytotoxicity for teclistamab

using Pkg; Pkg.activate("."); 
# Pkg.instantiate()  # run this line when running script in this repo for the first time 

using DifferentialEquations, ModelingToolkit 
using Plots
using DataFrames, CSV
using JLD2
using ModelingToolkit: getdefault
using Optimization, OptimizationNLopt

# observed data, from Pillarisetti 2020 H929 cells (Figure 1A); https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7509877/
obs = DataFrame(
    tcb_conc = [0.0005, 0.002, 0.01, 0.03, 0.1, 0.5, 2, 8.5, 33.5, 131., 512], 
    cytotox_median = [28.45, 23.04, 26.24, 24.32, 53.79, 84.21, 87.67, 88.71, 90.26, 89.72, 89.69], 
    cytotox_lowerrorbar = [23.35, 17.4, 20.6, 20.6, 44.45, 80.71, 84.71, 86.29, 87.61, 87.84, 87.54],
    cytotox_highererrorbar = [34.64, 27.08, 30.81, 27.82, 61.39, 87.16, 89.82, 91.13, 92.18, 91.87, 91.57]
); 

# load model 
include("tce_invitro_cleanup.jl");

Tumor_cells_init = 2E4
T_cells_init = Tumor_cells_init * 5
Tumor_cell_doubling_time = 50.
tspan = (0., 48*s_per_hr); 

@named tce = TCE_invitro(T_cells_init, Tumor_cells_init, Tumor_cell_doubling_time);
tce = structural_simplify(tce);
prob0 = ODEProblem(tce, [], tspan, []);

function cytotox(p_update, obs = obs; prob = prob0, runOpt = true)
    # update parameters 
    newparam = Dict([tce.lambda => exp(p_update[1]), 
                     tce.Thalf_TAA => exp(p_update[2]), 
                     tce.k_apop => exp(p_update[3]), 
                     tce.Thalf_CD3 => exp(p_update[4]), 
                     tce.ec50_kill => exp(p_update[5]) ]); 
    prob1 = remake(prob, p = newparam);
    # solve default case
    sol0 = solve(prob, saveat = s_per_hr);
    surv0 = [];
    for c_basb in obs.tcb_conc
        u0_2 = Dict([tce.TCE => c_basb]); 
        prob2 = remake(prob1, u0 = u0_2); 
        sol = solve(prob2, saveat = s_per_hr, alg=QNDF());
        append!(surv0, sol[:Tumor_cell][end]/sol0[:Tumor_cell][end])
        # append!(surv0, sol[:Tumor_cell][end]/Tumor_cells_init)
    end
    cytotox = (1 .- surv0)*100
    if runOpt
        return sum( (obs.cytotox_median .- cytotox).^2 );
    else
        return cytotox
    end
end

f_loss = OptimizationFunction(cytotox, Optimization.AutoForwardDiff());
sol = solve(OptimizationProblem(f_loss, log.([3., getdefault(tce.Thalf_TAA), 1.7E-5, getdefault(tce.Thalf_CD3), getdefault(tce.ec50_kill)]), obs ), NLopt.LN_NELDERMEAD() );

# visualization
cytotox0 = cytotox(sol, obs; runOpt = false);

p_opt = plot(xlabel = "Teclistamab concentration (nM)", ylabel = "Cytotoxicity (%)", xaxis = :log10, seriestype = :scatter, dpi = 300);
scatter!(obs.tcb_conc, obs.cytotox_median, ms = 6, yerror = (obs.cytotox_median .- obs.cytotox_lowerrorbar, obs.cytotox_highererrorbar .- obs.cytotox_median), label = "Pillarisetti 2020");
plot!(obs.tcb_conc, cytotox0, label = "Sim");
plot!(legend = :outerright);
display(p_opt)

savefig(p_opt, "figure/teclistamab-in-vitro.png");

# save outcome
using JLD2
jldsave("data/param/optim_teclistamab.jld2"; lambda=exp(sol.u[1]), thalf_taa=exp(sol.u[2]), k_apop=exp(sol.u[3]), Thalf_CD3 = exp(sol.u[4]), ec50_kill = exp(sol.u[5]))
