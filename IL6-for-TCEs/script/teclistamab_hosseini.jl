# date: 9/20/24
# author: Yuezhe Li 
# purpose of this code: IL6 simulation for teclistamab

using Pkg; Pkg.activate("")

using DifferentialEquations, ModelingToolkit 
using ModelingToolkit: getdefault
using Plots
using DataFrames
using DataFramesMeta
using CSV
using JLD2

# observed PK (IV), from Figure 6, BLA filing
# https://www.accessdata.fda.gov/drugsatfda_docs/nda/2022/761291Orig1s000MultidisciplineR.pdf
obs_iv_point3 = DataFrame(t = [0, 0.78, 2.34, 6.25, 13.28], conc = [0.00444, 0.00338, 0.00225, 0.00100, 0.00076]);
obs_iv_point6 = DataFrame(t = [0.77, 2.32, 3.87, 6.96, 14.7], conc = [0.00694, 0.00528, 0.00350, 0.00203, 0.00089]); 
obs_iv_2point4 = DataFrame(t = [0., 0.78, 3.1, 0, 0.78, 2.3, 0, 2.3, 3.1, 7, 7, 14, 14, 14], conc = [0.0706, 0.0468, 0.0355, 0.0236, 0.018, 0.0136, 0.0091, 0.0104, 0.0136, 0.0136, 0.004, 0.0059, 0.0015, 0.001]); 
obs_iv_9point6 = DataFrame(t = [1.7017, 2.4606, 1.648, 3.2, 4.8, 4.8, 8.6, 8.6, 15.6, 15.4], conc = [0.1851, 0.1408, 0.071, 0.0934, 0.0934, 0.063, 0.0471, 0.0313, 0.1228, 0.0091]);
obs_iv_19point2 = DataFrame(t = [0, 0.78, 0, 2.3, 3.1, 6.2, 13.2], conc = [0.4539, 0.2995, 0.1981, 0.1976, 0.1498, 0.0749, 0.0124] );

# observed PK, from Figure 9, BLA filing 
# dosing: 0.06mg/kg, 0.3mg/kg step up, then 1.5mg/kg QW
obs_sc_ = CSV.read("data/teclistamab_sc_JandJ.csv", DataFrame);
@select!(obs_sc_, :time_day = :time_weeks*7 .+ 7, :conc_ugmL);

# observed IL6, from Figure 1 of Willenmin 2024; # https://ascpt.onlinelibrary.wiley.com/doi/full/10.1002/psp4.13144
obs_il6 = CSV.read("data/teclistamab-IL6.csv", DataFrame);

# original hosseini model 
include("model/tdb_homo.jl");
@mtkbuild tdb = TDB_homo();

const BW = 74. # [kg]

# update PK parameters, obtained from teclistamab FDA BLA document
pk_teclistamab = Dict([tdb.Vm_tdb => 0, tdb.kabs_TDB => 0.14, tdb.CL_TDB => 999.8, tdb.Q_TDB => 47.3, tdb.V1_TDB => 4090., tdb.V2_TDB => 1290.]); 

# optimized in vitro parameters 
tec_in_vitro = load("data/params/optim_teclistamab_semimechanistic.jld2");
pd_teclistamab = Dict([tdb.KdrugactT => tec_in_vitro["KdrugactT"], tdb.KmBT_act => tec_in_vitro["KmBT_act"], tdb.KdrugB => tec_in_vitro["KdrugB"]]); 

pkpd_teclistamab = merge(pk_teclistamab, pd_teclistamab); 

# create base model 
prob0 = ODEProblem(tdb, [], (0., 14.), pkpd_teclistamab);

# simulation of different IV doses 
sol_point3 = solve(remake(prob0, u0 = Dict([tdb.TDBc_ugperml => 0.3*BW/pk_teclistamab[tdb.V1_TDB], tdb.TCEinjection_effect => 10.])), reltol=1e-18, saveat = 0.1); 
sol_point6 = solve(remake(prob0, u0 = Dict([tdb.TDBc_ugperml => 0.6*BW/pk_teclistamab[tdb.V1_TDB], tdb.TCEinjection_effect => 10.])), reltol=1e-18, saveat = 0.1); 
sol_2point4 = solve(remake(prob0, u0 = Dict([tdb.TDBc_ugperml => 2.4*BW/pk_teclistamab[tdb.V1_TDB], tdb.TCEinjection_effect => 10.])), reltol=1e-18, saveat = 0.1); 
sol_9point6 = solve(remake(prob0, u0 = Dict([tdb.TDBc_ugperml => 9.6*BW/pk_teclistamab[tdb.V1_TDB], tdb.TCEinjection_effect => 10.])), reltol=1e-18, saveat = 0.1); 
sol_19point2 = solve(remake(prob0, u0 = Dict([tdb.TDBc_ugperml => 19.2*BW/pk_teclistamab[tdb.V1_TDB], tdb.TCEinjection_effect => 10.])), reltol=1e-18, saveat = 0.1); 

p_pk_iv = plot(yaxis = :log10, ylims = [1E-4, 1E2], xlabel = "Time (days)", ylabel = "Serum conc (ug/mL)", legend = :outerright, dpi = 300);
plot!(sol_point3.t, sol_point3[:TDBc_ugperml], label = "sims, 0.3ug/kg", color = :darkgreen, alpha = 0.5);
plot!(obs_iv_point3.t, obs_iv_point3.conc, label = "obs, 0.3ug/kg", seriestype = :scatter, color = :darkgreen, alpha = 0.5);
plot!(sol_point6.t, sol_point6[:TDBc_ugperml], label = "sims, 0.6ug/kg", color = :turquoise3, alpha = 0.5);
plot!(obs_iv_point6.t, obs_iv_point6.conc, label = "obs, 0.6ug/kg", seriestype = :scatter, color = :turquoise3, alpha = 0.5);
plot!(sol_2point4.t, sol_2point4[:TDBc_ugperml], label = "sims, 1.2ug/kg", color = :gray16, alpha = 0.5);
plot!(obs_iv_2point4.t, obs_iv_2point4.conc, label = "obs, 2.4ug/kg", seriestype = :scatter, color = :gray16, alpha = 0.5);
plot!(sol_9point6.t, sol_9point6[:TDBc_ugperml], label = "sims, 9.6ug/kg", color = :slateblue2, alpha = 0.5);
plot!(obs_iv_9point6.t, obs_iv_9point6.conc, label = "obs, 9.6ug/kg", seriestype = :scatter, color = :slateblue2, alpha = 0.5);
plot!(sol_19point2.t, sol_19point2[:TDBc_ugperml], label = "sims, 19.2ug/kg", color = :firebrick4, alpha = 0.5);
plot!(obs_iv_19point2.t, obs_iv_19point2.conc, label = "obs, 9.6ug/kg", seriestype = :scatter, color = :firebrick4, alpha = 0.5);
display(p_pk_iv); 

savefig(p_pk_iv, "deliv/figure/invivo/teclistamab-pk-iv.png")

# simulation for sc dosing 
prob_sc = ODEProblem(tdb, [], (-1E-12, 35.), pkpd_teclistamab);

function tdb_sc_dosing(dose_amt, dose_time; injection_effect_init = 10.)
    global cbs_human = [];
    if length(dose_time) > 1
        for i in 1:length(dose_time)
            function affect!(integrator)
                integrator[:TDBdepot_ug] += dose_amt[i];
                integrator[:TCEinjection_effect] += injection_effect_init;
            end
            cb = PresetTimeCallback(dose_time[i],affect!);
            global cbs_human = push!(cbs_human, cb);
        end
    end
    cbset_human = CallbackSet(cbs_human...);
    return cbset_human
end

# dosing event 
dose_amt = [0.06, 0.3, 1.5, 1.5, 1.5, 1.5] * BW * 1e3;  # [ug]
dose_time = [0, 3, 7, 14, 21, 28];                 # [day] 
cbset_sc = tdb_sc_dosing(dose_amt, dose_time); 

sol_sc = solve(prob_sc, callback = cbset_sc, reltol = 1E-12, saveat = 0.1);

p_pk_sc = plot(yaxis = :log10, ylims = [1E-2, 1E2], xlabel = "Time (days)", ylabel = "Serum conc (ug/mL)", legend = :outerright, dpi = 300);
plot!(sol_sc.t, sol_sc[:TDBc_ugperml], label = "sims");
plot!(obs_sc_.time_day, obs_sc_.conc_ugmL, label = "obs", seriestype = :scatter);
display(p_pk_sc);

savefig(p_pk_sc, "deliv/figure/invivo/teclistamab-pk-sc.png")

# T cell & B cell check
p_tcell = plot(yaxis = :log10, xlabel = "Time (days)", ylabel = "peripheral blood T cell (#/uL)", legend = :outerright, dpi = 300);
plot!(sol_sc.t, (sol_sc[:restTpb] .+ sol_sc[:actTpb] .+ sol_sc[:act0Tpb])/(getdefault(tdb.Vpb)*1E3), label = "sims");
display(p_tcell);

p_bcell = plot(yaxis = :log10, xlabel = "Time (days)", ylabel = "peripheral blood B cell (#/uL)", legend = :outerright, dpi = 300);
plot!(sol_sc.t, sol_sc[:Bpb]/(getdefault(tdb.Vpb)*1E3), label = "sims");
display(p_bcell);

# IL6 prediction 
# protocol obtained from https://www.nejm.org/doi/full/10.1056/NEJMoa2203478

p_il6 = plot(xlabel = "Time (hours)", ylabel = "plasma IL6 conc (pg/mL)", legend = :outerright, dpi = 300);
plot!(sol_sc.t * 24, sol_sc[:IL6pb], label = "sims");
plot!(obs_il6.time_hr, obs_il6.IL6_ngL, label = "mean IL6 from MajesTEC-1 trial", seriestype = :scatter, alpha = 0.8);
plot!(xlims = [0, 500], ylims = [0.1, 100], yaxis = :log10);
display(p_il6);

savefig(p_il6, "deliv/figure/invivo/teclistamab-il6-sc.png")
