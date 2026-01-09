# date: 10/16/24 
# author: Yuezhe Li 
# purpose of this code: to clean up mosun in vitro

using Pkg; Pkg.activate("")

using DifferentialEquations, ModelingToolkit 
using DataFrames
using ModelingToolkit: getdefault
using DataFrames
using JLD2
using Plots
using ComponentArrays

include("tce_invitro_cleanup.jl");

#--------------- obs ---------------#
# in vitro data from Bacac et al., 2018, Fig 1B; # https://pubmed.ncbi.nlm.nih.gov/29716920/
tcb_conc = [0.001, 0.01, 0.1, 1., 10., 100., 1000.]; # [pmol/L]
toledo_monocd20 = DataFrame(TCBconc = tcb_conc, tumor_lysis = [-3.7, -3.7, -3.7, -3.7, -0.85, 32.32, 46.65]);

#--------------- base model ---------------#
# initial condition obtained from Bacac et al., 2018; # https://pubmed.ncbi.nlm.nih.gov/29716920/
# Toledo cell line doubling time obtained from https://www.cellosaurus.org/CVCL_3611
# in vitro binding affinities obtained from BLA filing; https://www.accessdata.fda.gov/drugsatfda_docs/nda/2023/761263Orig1s000MultidisciplineR.pdf
# CD20 expression level obtaiend from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC500695/pdf/jclinpath00266-0020.pdf, though different from https://pubmed.ncbi.nlm.nih.gov/8817090/
# CD3 expression level was obtained from https://pubmed.ncbi.nlm.nih.gov/8817090/ and didn't change here 
Tumor_cells_init = 1.5E5
T_cells_init = Tumor_cells_init * 5
Tumor_cell_doubling_time = 24.
tspan = (0., 20*s_per_hr); 

@named tce = TCE_invitro(T_cells_init, Tumor_cells_init, Tumor_cell_doubling_time);
tce = structural_simplify(tce);
prob0 = ODEProblem(tce, [], tspan, [tce.KD_CD3 => 40., tce.KD_TAA => 63., tce.TAA_per_cell => 9.4E4]);

function cytotox(tcb_conc, p_update; prob = prob0)
    # update parameters 
    newparam = Dict([tce.lambda => p_update.lambda, 
                     tce.Thalf_TAA => p_update.Thalf_TAA, 
                     tce.k_apop => p_update.k_apop, 
                     tce.ec50_kill => p_update.ec50_kill]); 
    prob1 = remake(prob, p = newparam);
    # solve default case
    sol0 = solve(prob1, saveat = s_per_hr);
    surv0 = [];
    for c_basb in tcb_conc
        u0_2 = Dict([tce.TCE => c_basb/1E3]); 
        prob2 = remake(prob1, u0 = u0_2); 
        sol = solve(prob2, saveat = s_per_hr, alg=QNDF());
        # display(plot(sol.t/s_per_hr, sol[:Tumor_cell]./sol0[:Tumor_cell], ylims = [0, 1], title = string(c_basb), label = false));
        append!(surv0, max(sol[:Tumor_cell][end], 0.)/sol0[:Tumor_cell][end])
    end
    return (1 .- surv0)*100
end

p_update = ComponentVector(lambda = getdefault(tce.lambda), Thalf_TAA = getdefault(tce.Thalf_TAA), k_apop = getdefault(tce.k_apop), ec50_kill = 1 );
cytotox0 = cytotox(tcb_conc, p_update);

p_preopt = plot(xlabel = "TCE conc (nM)", ylabel = "Cytotoxicity (%)", xaxis = :log10, legend = :topleft, dpi = 300);
plot!(tcb_conc/1E3, cytotox0, label = "sims");
plot!(tcb_conc/1E3, toledo_monocd20.tumor_lysis, label = "Bacac 2018", seriestype=:scatter);
display(p_preopt);


#--------------- optimization ---------------#
function obj(p, obs; tcb_conc = tcb_conc, prob = prob0)
    newsims = cytotox(tcb_conc, p, prob = prob);
    return sum( (obs .- newsims).^2 );
end

using Optimization, OptimizationOptimJL
opt_prob = OptimizationProblem(obj, p_update, toledo_monocd20.tumor_lysis);
solopt = solve(opt_prob, NelderMead(), maxiters=1E4);

cytotox_opt = cytotox(tcb_conc, solopt.u);

p_postopt = plot(xlabel = "TCE conc (nM)", ylabel = "Cytotoxicity (%)", xaxis = :log10, legend = :topleft, dpi = 300);
plot!(tcb_conc/1E3, cytotox_opt, label = "sims");
plot!(tcb_conc/1E3, toledo_monocd20.tumor_lysis, label = "Bacac 2018", seriestype=:scatter);
display(p_postopt);

# aim to generate the bell curve
tcb_conc2 = [0.001, 0.01, 0.1, 1., 10., 100., 1E3, 1E4, 1E5, 1E6, 1E7, 1E8, 1E9, 1E10]; # [pmol/L]

cytotox2 = cytotox(tcb_conc2, solopt.u);

p_2 = plot(xlabel = "Mosunetuzumab concentration (nM)", ylabel = "Cytotoxicity (%)", xaxis = :log10, legend = :topleft, dpi = 300);
plot!(tcb_conc2/1E3, cytotox2, label = "sims");
plot!(tcb_conc/1E3, toledo_monocd20.tumor_lysis, label = "Bacac 2018", seriestype=:scatter);
display(p_2);

savefig(p_2, "figure/mosunetuzumab-in-vitro.png")

using JLD2
jldsave("data/param/optim_mosunetuzumab.jld2"; lambda=solopt.u.lambda, thalf_taa=solopt.u.Thalf_TAA, k_apop=solopt.u.k_apop, ec50_kill = solopt.u.ec50_kill)
