# date: 9/24/24
# purpose of this code: fit for in vitro model of epcoritamab 

using Pkg; Pkg.activate("")

using DifferentialEquations, ModelingToolkit 
using ModelingToolkit: getdefault
using DataFrames
using JLD2
using Plots
using ComponentArrays

# in vitro data from Engelberts et al., 2020, Fig 1F; # https://pubmed.ncbi.nlm.nih.gov/31981978/
const MW_EDG = 150000.0; # Da
epco_conc_RI1 = [2.88306E-06, 3.24301E-05, 3.29798E-04, 2.88306E-03, 3.41073E-02, 3.13581E-01, 3.70974, 34.1073, 298.162]; # [ng/mL]
RI1_epco = DataFrame(TCBconc = epco_conc_RI1, tumor_lysis = [3.12269, 9.23908E-01, 3.25167, 5.17748, 18.1812, 50.1040, 60.6419, 67.4963, 67.7782]);

epco_conc_daudi = [3.3015E-06, 3.31284E-05, 0.000337174, 0.003052746, 0.028595947, 0.262800248, 2.578298474, 23.79369186, 249.6927377]; # [ng/mL]
daudi_epco = DataFrame(TCBconc = epco_conc_daudi, tumor_lysis = [2.957686583, -8.456304486, -5.728258558, 5.277482458, 32.24431214, 67.08909207, 84.96704231, 88.09483309, 88.19795875]);

include("model/tce_invitro_cleanup.jl");

# initial condition obtained from Engelbert et al., 2020; # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6992935/
# Daudi cell line doubling time obtained from https://www.cellosaurus.org/CVCL_0008
# in vitro binding affinities obtained from BLA filing; https://www.accessdata.fda.gov/drugsatfda_docs/nda/2023/761324Orig1s000MultidisciplineR.pdf
# CD20 expression level range from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC500695/pdf/jclinpath00266-0020.pdf and https://pubmed.ncbi.nlm.nih.gov/8817090/
# CD3 expression level was obtained from https://pubmed.ncbi.nlm.nih.gov/8817090/ and didn't change here 
Tumor_cells_init = 1.5E4
T_cells_init = Tumor_cells_init * 2
Tumor_cell_doubling_time = 23.7
tspan = (0., 48*s_per_hr); 

@named tce = TCE_invitro(T_cells_init, Tumor_cells_init, Tumor_cell_doubling_time);
tce = structural_simplify(tce);
prob0 = ODEProblem(tce, [], tspan, [tce.KD_CD3 => 4.73, tce.KD_TAA => 2.47]);

function cytotox(p_update, obs = RI1_epco; prob = prob0, runOpt = true)
    # update parameters 
    newparam = Dict([tce.lambda => p_update[1], 
                     tce.Thalf_TAA => p_update[2], 
                     tce.k_apop => p_update[3], 
                     tce.TAA_per_cell => p_update[4]]); 
    prob1 = remake(prob, p = newparam);
    # solve default case
    sol0 = solve(prob, saveat = s_per_hr);
    surv0 = [];
    for c_basb in obs.TCBconc
        u0_2 = Dict([tce.TCE => c_basb*1E-6/MW_EDG*1E9]); 
        prob2 = remake(prob1, u0 = u0_2); 
        sol = solve(prob2, saveat = s_per_hr, alg=QNDF());
        append!(surv0, sol[:Tumor_cell][end]/sol0[:Tumor_cell][end])
        # append!(surv0, sol[:Tumor_cell][end]/Tumor_cells_init)
    end
    cytotox = (1 .- surv0)*100
    if runOpt
        return sum( (obs.tumor_lysis .- cytotox).^2 );
    else
        return cytotox
    end
end

p0 = [getdefault(tce.lambda), getdefault(tce.Thalf_TAA), getdefault(tce.k_apop), getdefault(tce.TAA_per_cell)]
cytotox0 = cytotox(p0, daudi_epco, runOpt = false);

p_preopt = plot(title = "pre-optimization", titlefontsize = 8, xlabel = "bsAb conc (ng/mL)", ylabel = "Cytotoxicity (%)", xaxis = :log10, legend = :outerright);
plot!(RI1_epco.TCBconc, cytotox0, label = "sims");
plot!(RI1_epco.TCBconc, RI1_epco.tumor_lysis, label = "Engelbert 2020", seriestype=:scatter);
display(p_preopt);

# Optimization 
using Optimization, OptimizationBBO

opt_prob = OptimizationProblem(cytotox, p0, daudi_epco, lb = [1., 0.1, 1E-10, 1E3], ub = [50., 10., 1E-4, 150E3]);
solopt = solve(opt_prob, BBO_adaptive_de_rand_1_bin_radiuslimited(), show_progress = true); # [33.94, 9.99, 3.4E-6, 37387]
jldsave("data/params/optim_epcoritamab.jld2"; lambda=solopt.u[1], thalf_taa=solopt.u[2], k_apop=solopt.u[3], TAA_per_cell = solopt.u[4]);

cytotox_opt = cytotox(solopt.u, daudi_epco, runOpt = false);

p_postopt = plot(title = "post-optimization", titlefontsize = 8, xlabel = "bsAb conc (ng/mL)", ylabel = "Cytotoxicity (%)", xaxis = :log10, legend = :outerright);
plot!(daudi_epco.TCBconc, cytotox_opt, label = "sims");
plot!(daudi_epco.TCBconc, daudi_epco.tumor_lysis, label = "Engelbert 2020", seriestype=:scatter);
display(p_postopt);

# expand bsAb conc for the bell curve
TCBconc2 = [1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1, 1E0, 1E1, 1E2, 1E3, 1E4, 1E5, 1E6, 1E7]
sims_tcb = DataFrame(TCBconc = TCBconc2, tumor_lysis = zeros(size(TCBconc2))); # placeholder
cytotox2 = cytotox(solopt.u, sims_tcb, runOpt = false);

p_bell_curve = plot(xlabel = "epcoritamab conc (ng/mL)", ylabel = "Cytotoxicity (%)", xaxis = :log10, legend = :bottomright);
plot!(TCBconc2, cytotox2, label = "sims");
plot!(daudi_epco.TCBconc, daudi_epco.tumor_lysis, label = "Engelbert 2020", seriestype=:scatter);
plot!(xticks = TCBconc2);
display(p_bell_curve);

# save figure 
savefig(p_bell_curve, "deliv/figure/invitro/epcoritamab_bell_curve.png");
