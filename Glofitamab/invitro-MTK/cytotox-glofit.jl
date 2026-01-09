# date: 2/28/2025
# author: Jimena Davis
# purpose of this code: test in vitro TCE model against public glofitamab cytotox data

using ProjectRoot
using Pkg; Pkg.activate(@projectroot)

using DifferentialEquations, ModelingToolkit 
using ModelingToolkit: getdefault
using SymbolicIndexingInterface
using DataFrames
using Statistics
using Plots

#--------------- obs ---------------#
# in vitro data from Bacac et al., 2018, Fig 1B (Toledo cells); # https://pubmed.ncbi.nlm.nih.gov/29716920/
tcb_conc = [0.01, 0.1, 1., 10., 100., 1000.]; # [pmol/L]
toledo_glofit = DataFrame(TCBconc = tcb_conc, tumor_lysis = [3.49, 14.89, 29.96, 41.51, 47.16, 49.45]);

#--------------- base model ---------------#
include(@projectroot("model_constants.jl"));
include(@projectroot("invitro/invitro_tce.jl"));

Tumor_cells_init = 1000 # placeholder 
ETratio = 5
T_cells_init = Tumor_cells_init * ETratio
Tumor_cell_doubling_time = 24. # [hr]; Toledo cell line doubling time obtained from https://www.cellosaurus.org/CVCL_3611
tspan = (0., 24*s_per_hr); 

@mtkbuild tce = invitro(Tumor_cell_doubling_time);
prob0 = ODEProblem(tce, [], tspan, [tce.TAA_per_cell => 9.4E4]); # CD20 expression level obtained from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC500695/pdf/jclinpath00266-0020.pdf, though different from https://pubmed.ncbi.nlm.nih.gov/8817090/

#--------------- output function ---------------#
function varies_readout(bsab, prob = prob0)
    # solve system when TCE conc = 0
    sol0 = solve(remake(prob, u0 =  Dict([tce.TCE => 0]) ), saveat = s_per_hr, alg=AutoTsit5(Rosenbrock23()), abstol=1e-6, reltol=1e-3);
    tumor_cell_0 = sol0[tce.Tumor_cell][end]

    # save outcomes 
    trimer_per_T = []
    trimer_per_tumor = []
    taa_occupancy = []
    cd3_occupancy = []
    tumor_cell = []
    trimer_1_per_tumor_cell =[]
    trimer_2_per_tumor_cell = []
    time = []
    for init_bsab_conc in bsab
        tmp = solve(remake(prob, u0 =  Dict([tce.TCE => init_bsab_conc]) ), saveat = s_per_hr, alg=AutoTsit5(Rosenbrock23()), abstol=1e-6, reltol=1e-3);
        append!(trimer_per_T, tmp[tce.trimer_per_T][end]);
        append!(trimer_per_tumor, tmp[tce.trimer_per_Tumor][end]);
        append!(taa_occupancy, tmp[tce.R0_TAA][end]);
        append!(cd3_occupancy, tmp[tce.R0_CD3][end]);
        append!(tumor_cell, tmp[tce.Tumor_cell][end]);
        append!(trimer_1_per_tumor_cell, tmp[tce.trimer_1_per_tumor_cell][end]);
        append!(trimer_2_per_tumor_cell, tmp[tce.trimer_2_per_tumor_cell][end]);
        append!(time, tmp.t[end]); # saved to check for potential solver issue
    end
    tmp_Dict = Dict("trimer_per_T_cell" => trimer_per_T, "trimer_per_tumor_cell" => trimer_per_tumor, 
                    "TAA_RO" => taa_occupancy * 100, "CD3_RO" => cd3_occupancy * 100, 
                    "cytotox" => (1 .- tumor_cell/tumor_cell_0) * 100,
                    "CD3_TCE_TAA" => trimer_1_per_tumor_cell, "CD3_TCE_TAA_TAA" => trimer_2_per_tumor_cell, 
                    "Endtime" => time) 
    return tmp_Dict
end


#--------------- glofitamab param scan ---------------#
# expanded TCE concentration; aim to generate the bell curve
tcb_conc2 = [0.001, 0.01, 0.1, 1., 10., 100., 1E3, 1E4, 1E5, 1E6, 1E7, 1E8, 1E9, 1E10]; # [pmol/L]

# glofitamab KD for CD20 was obtained from obinutuzumab BLA; https://www.accessdata.fda.gov/drugsatfda_docs/nda/2013/125486Orig1s000PharmR.pdf
# glofitmab was based on obinutuzumab; Bacac 2018; https://pubmed.ncbi.nlm.nih.gov/29716920/
# Binding affinity from cell binding data

update_KD_CD3 = 21.1
update_KD_TAA = 4.3
update_Thalf_TAA = 200
try_tt50_kill = 0.01
try_emax_apop = 8e-6
try_tt50_t_activate = 0.01
try_emax_t_activate = 8e-6
try_n_act = 0.6 

prob_glo = remake(prob0, p = Dict([tce.KD_CD3 => update_KD_CD3, tce.KD_TAA => update_KD_TAA, tce.Thalf_TAA => update_Thalf_TAA, tce.tt50_kill => try_tt50_kill, tce.emax_apop => try_emax_apop, tce.Tumor_cells_init => Tumor_cells_init, tce.total_T_cells_init => T_cells_init, tce.tt50_t_activate => try_tt50_t_activate, tce.emax_t_activate => try_emax_t_activate, tce.n_act => try_n_act])); 

output_0 = varies_readout(tcb_conc2*1E-3, prob_glo);
output_2 = varies_readout(tcb_conc2*1E-3, remake(prob_glo, p = Dict([tce.lambda_taa_avidity => 1E2])));
output_4 = varies_readout(tcb_conc2*1E-3, remake(prob_glo, p = Dict([tce.lambda_taa_avidity => 1E4])));
output_6 = varies_readout(tcb_conc2*1E-3, remake(prob_glo, p = Dict([tce.lambda_taa_avidity => 1E6])));

plot_cis_lambda_cytotox = plot(xaxis = :log10, xlabel = "Glofitamab concentration (nM)", ylabel = "Cytotoxicity (%)", title = "Cytotoxicity in Toledo", dpi = 300);
plot!(tcb_conc2/1E3, output_0["cytotox"], label = "cis avidity = 1", lw = 2); 
plot!(tcb_conc2/1E3, output_2["cytotox"], label = "cis avidity = 1E2", lw = 2); 
plot!(tcb_conc2/1E3, output_4["cytotox"], label = "cis avidity = 1E4", lw = 2); 
plot!(tcb_conc2/1E3, output_6["cytotox"], label = "cis avidity = 1E6", lw = 2); 
plot!(tcb_conc/1E3, toledo_glofit.tumor_lysis, label = "Bacac et al., 2018 (glofit)", seriestype = :scatter); 
plot!(legend = :bottomright, background_color_legend = nothing);
plot!(xticks = tcb_conc2/1E3);
display(plot_cis_lambda_cytotox)

savefig(plot_cis_lambda_cytotox, @projectroot("figures/invitro/glofitamab-cytotox.png")) 

plot_cis_lambda_ro = plot(xaxis = :log10, xlabel = "Glofit concentration (nM)", ylabel = "CD20 occupancy (%)", title ="Receptor Occupancy", dpi = 300);
plot!(tcb_conc2/1E3, output_0["TAA_RO"], label = "cis avidity = 1", lw = 2); 
plot!(tcb_conc2/1E3, output_2["TAA_RO"], label = "cis avidity = 1E2", lw = 2); 
plot!(tcb_conc2/1E3, output_4["TAA_RO"], label = "cis avidity = 1E4", lw = 2); 
plot!(tcb_conc2/1E3, output_6["TAA_RO"], label = "cis avidity = 1E6", lw = 2); 
plot!(legend = :bottomright, background_color_legend = nothing);
display(plot_cis_lambda_ro)

savefig(plot_cis_lambda_ro, @projectroot("figures/invitro/glofitamab-cytotox-RO.png")) 
