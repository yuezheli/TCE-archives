# date: 9/18/24
# author: Yuezhe Li 
# purpose of this code: to optimize a hypothetical TCE to generate parameters for Hosseini style in vitro model 

using Pkg; Pkg.activate("")

using DifferentialEquations, ModelingToolkit 
using Plots

include("model/tce_invitro_cleanup.jl");

#--------------- base model ---------------#
# parameters obtained from Utsey et al., 2023; # https://www.metrumrg.com/wp-content/uploads/2023/11/poster-final.pdf
# parameters for teclistamab
Tumor_cells_init = 1.5E3
T_cells_init = Tumor_cells_init * 5
Tumor_cell_doubling_time = 50.
tspan = (0., 48*s_per_hr); 

@named tce = TCE_invitro(T_cells_init, Tumor_cells_init, Tumor_cell_doubling_time);
tce = structural_simplify(tce);

param_teclistamab = Dict([tce.KD_CD3 => 28.03, tce.KD_TAA => 0.18, tce.CD3_per_cell => 3E4, tce.TAA_per_cell => 13173.0, tce.k_apop => 1.7E-5, tce.Thalf_CD3 => 0.1, tce.Thalf_TAA => 24.])

prob0 = ODEProblem(tce, [], tspan, param_teclistamab);

function cytotox(tcb_conc; prob = prob0)
    # solve default case
    sol0 = solve(prob, saveat = s_per_hr);
    surv0 = [];
    for c_basb in tcb_conc
        u0_2 = Dict([tce.TCE => c_basb]); 
        prob2 = remake(prob, u0 = u0_2); 
        sol = solve(prob2, saveat = s_per_hr, alg=QNDF());
        append!(surv0, max(sol[:Tumor_cell][end], 0.)/sol0[:Tumor_cell][end])
    end
    return (1 .- surv0)*100
end

tcb_conc0 = [1E-4, 1E-3, 2E-3, 1E-2, 2E-2, 0.1, 0.2, 0.5, 1., 2., 5., 1E1, 1E2, 1E3]; # [nM]
cytotox0 = cytotox(tcb_conc0);

# plot(tcb_conc0, cytotox0, label = false, xlabel = "TCE conc (nM)", ylabel = "Cytotoxicity (%)", xaxis = :log10, seriestype = :scatter)

#--------------- optimization towards Hosseini model ---------------#
using Optimization, OptimizationOptimJL

include("model/tdb-in-vitro-hosseini.jl");

p_base = ComponentArray(
    kBtumorprolif = p_homo_3.kBtumorprolif, KmBT_act = p_homo_3.KmBT_act, S = p_homo_3.S, KdrugactT = p_homo_3.KdrugactT, ndrugactT = p_homo_3.ndrugactT, VmT = p_homo_3.VmT, 
    KmTB_kill = p_homo_3.KmTB_kill, nkill = p_homo_3.nkill, KdrugB = p_homo_3.KdrugB, VmB = p_homo_3.VmB, kTact = p_homo_3.kTact, fTadeact = p_homo_3.fTadeact, fTa0deact = p_homo_3.fTa0deact); 

const MW_EDG = 150000.0; # Da
u0_0 = ComponentArray(Btumor = 1.5E3, bs_ugperml = 0., restT = 7.5E3, actT = 0., act0T = 0.); 

# observed data, from Pillarisetti 2020 H929 cells (Figure 1A); https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7509877/
using DataFrames
obs = DataFrame(
    tcb_conc = [0.0005, 0.002, 0.01, 0.03, 0.1, 0.5, 2, 8.5, 33.5, 131., 512], 
    cytotox_median = [28.45, 23.04, 26.24, 24.32, 53.79, 84.21, 87.67, 88.71, 90.26, 89.72, 89.69], 
    cytotox_lowerrorbar = [23.35, 17.4, 20.6, 20.6, 44.45, 80.71, 84.71, 86.29, 87.61, 87.84, 87.54],
    cytotox_highererrorbar = [34.64, 27.08, 30.81, 27.82, 61.39, 87.16, 89.82, 91.13, 92.18, 91.87, 91.57]
); 

function loss2(p, obs = obs; runOpt = true)
    p_new = deepcopy(p_base);
    p_new.KmBT_act = exp(p[1]);
    p_new.KdrugactT = exp(p[2]);
    p_new.KdrugB = exp(p[3]);
    test_new = [];
    for conc in obs.tcb_conc
        u0_tmp = deepcopy(u0_0);
        u0_tmp.bs_ugperml = conc/1E9 * MW_EDG * 1E3; # [ug/mL]
        sol_tmp = solve(ODEProblem(bs_in_vitro_hosseini!, u0_tmp, (0., 2.), p_new), reltol = 1E-18); # time in days
        tmp_kill_ratio = 1 - sol_tmp.u[end].Btumor/u0_tmp.Btumor
        append!(test_new, tmp_kill_ratio)
    end
    diff = sum((test_new * 100 .- obs.cytotox_median).^2);
    if runOpt
        return diff
    else
        return test_new * 100
    end
end

p0 = log.([p_homo_3.KmBT_act, p_homo_3.KdrugactT, p_homo_3.KdrugB]);
f2 = OptimizationFunction(loss2, Optimization.AutoForwardDiff());
oprob2 = OptimizationProblem(f2, p0, obs);
psol2 = solve(oprob2, NelderMead()) 

using JLD2
jldsave("data/params/optim_teclistamab_semimechanistic.jld2"; KmBT_act = exp(psol2.u[1]), KdrugactT = exp(psol2.u[2]), KdrugB = exp(psol2.u[3]));

cytotox_opt = loss2(psol2.u, runOpt = false);

p_opt = plot(xlabel = "Teclistamab concentration (nM)", ylabel = "Cytotoxicity (%)", xaxis = :log10, seriestype = :scatter);
scatter!(obs.tcb_conc, obs.cytotox_median, ms = 6, yerror = (obs.cytotox_median .- obs.cytotox_lowerrorbar, obs.cytotox_highererrorbar .- obs.cytotox_median), label = "Pillarisetti 2020");
plot!(obs.tcb_conc, cytotox_opt, label = "Optimized");
plot!(legend = :outerright);
display(p_opt)

savefig(p_opt, "deliv/figure/invitro/teclistamab-optimization.png");
