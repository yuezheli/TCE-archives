# date: 9/23/24
# author: Yuezhe Li 
# purpose of this code: simulation of Hosseini model on epcoritamab IL6 response 

using Pkg; Pkg.activate("")

using DifferentialEquations, ModelingToolkit 
using ModelingToolkit: getdefault
using Plots
using CSV, DataFrames

# observed peripheral blood T cell count; Figure 4 of Hutchings et al., 2021; # https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(21)00889-8/abstract
obs_Tcount = DataFrame(
    time_day = [0, 0.25, 1, 3, 7, 7.25, 8, 13, 14, 14.25, 15, 21, 21.25, 22], 
    Tcell_peruL_median = [552.48, 243.25, 332.21, 454.44, 467.95, 183.95, 184.01, 315.69, 315.98, 135.23, 110.96, 335.62, 291.51, 365.72], 
    Tcell_peruL_lowIQR = [428.23, 138.08, 211.21, 250.77, 352.56, 113.68, 78.70, 200.69, 184.52, 66.63, 66.66, 201.63, 160.86, 219.71], 
    Tcell_peruL_highIQR = [776.00, 393.60, 522.58, 620.47, 676.11, 289.35, 324.15, 418.97, 456.55, 266.79, 190.00, 608.22, 471.71, 626.26]
);

# observed IL6; Figure 4 of Hutchings et al., 2021; # https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(21)00889-8/abstract
obs_IL6 = DataFrame(
    time_day = [0, 0.25, 1, 3, 7, 7.25, 8, 10, 14, 14.25, 15, 18, 21, 21.25, 22, 25, 28, 28.05, 35, 35.05], 
    IL6_pgmL_median = [1.78, 0.7, 1.99, 1.49, 2.47, 1.72, 1.72, 1.99, 2.76, 0.72, 17.47, 9.1, 2.86, 0.7, 2.06, 2.22, 2.47, 1.78, 4.1, 2.22], 
    IL6_pgmL_lowIQR = [0.67, 0.7, 0.72, 0.7, 0.72, 0.7, 0.7, 0.7, 1.6, 0.67, 4.74, 3.19, 1.6, 0.67, 0.7, 0.72, 1.49, 0.72, 1.08, 0.7], 
    IL6_pgmL_highIQR = [3.19, 2.3, 4.74, 2.3, 4.92, 5.29, 7.88, 9.44, 9.79, 3.82, 183.85, 59.84, 44.79, 11.31, 53.68, 38.75, 18.78, 9.44, 24.2, 16.25]
);

# epco PK were obtained from popPK sims used in the BLA filing: 
# https://www.accessdata.fda.gov/drugsatfda_docs/nda/2023/761324Orig1s000MultidisciplineR.pdf
pk_median = CSV.read("data/epco_pk/epco-pk-med.csv", DataFrame);
pk_low = CSV.read("data/epco_pk/epco-pk-lo.csv", DataFrame);
pk_high = CSV.read("data/epco_pk/epco-pk-hi.csv", DataFrame);

pk_48_bla = vcat(pk_median, pk_low, pk_high);

# parameter updates for epcoritamab
include("model/ParamUpdate.jl");
p_epco = TDB_param_update("epcoritamab");

# original hosseini model 
include("model/tdb_homo.jl");
@mtkbuild tdb = TDB_homo();

# update Hosseini model to epcoritamab
prob0 = ODEProblem(tdb, [], (-1E-12, 84.), [tdb.tumor_on => 1., tdb.CL_TDB => p_epco.CL_TDB, tdb.Q_TDB => p_epco.Q_TDB, tdb.V1_TDB => p_epco.V1_TDB, tdb.V2_TDB => p_epco.V2_TDB, tdb.Vm_tdb => 0., tdb.KmBT_act => p_epco.KmBT_act, tdb.KdrugactT => p_epco.KdrugactT, tdb.KdrugB => p_epco.KdrugB]);

# dosing function 
function tdb_sc_dosing(dose_amt, dose_time; injection_effect_init = 10., bioavailability = 0.95)
    global cbs_human = [];
    if length(dose_time) > 1
        for i in 1:length(dose_time)
            function affect!(integrator)
                integrator[:TDBdepot_ug] += dose_amt[i] * bioavailability;
                integrator[:TCEinjection_effect] += injection_effect_init;
            end
            cb = PresetTimeCallback(dose_time[i],affect!);
            global cbs_human = push!(cbs_human, cb);
        end
    end
    cbset_human = CallbackSet(cbs_human...);
    return cbset_human
end

dose_time = 7 * [0, 1, 2, 3, 4, 5, 6, 7, 8, 10];
dose_amount_12 = [0.16, 0.8, 12, 12, 12, 12, 12, 12, 12, 12] * 1E3; # [ug]
dose_amount_24 = [0.16, 0.8, 24, 24, 24, 24, 24, 24, 24, 24] * 1E3; # [ug]
dose_amount_48 = [0.16, 0.8, 48, 48, 48, 48, 48, 48, 48, 48] * 1E3; # [ug]
dose_amount_60 = [0.16, 0.8, 60, 60, 60, 60, 60, 60, 60, 60] * 1E3; # [ug]

cbset_12 = tdb_sc_dosing(dose_amount_12, dose_time);
cbset_24 = tdb_sc_dosing(dose_amount_24, dose_time);
cbset_48 = tdb_sc_dosing(dose_amount_48, dose_time);
cbset_60 = tdb_sc_dosing(dose_amount_60, dose_time);

sol_0_12 = solve(prob0, reltol=1e-18, saveat = 0.2, callback = cbset_12); 
sol_0_24 = solve(prob0, reltol=1e-18, saveat = 0.2, callback = cbset_24); 
sol_0_48 = solve(prob0, reltol=1e-18, saveat = 0.2, callback = cbset_48); 
sol_0_60 = solve(prob0, reltol=1e-18, saveat = 0.2, callback = cbset_60);

## plot T cell 
p_tcell = plot(ylims = [40, 1E4], yscale = :log10, xlabel = "Time (day)", ylabel = "T cell count (/uL peripheral blood)", legend = :bottomright, dpi = 300);
plot!(sol_0_12.t, (sol_0_12[:restTpb] .+ sol_0_12[:actTpb] .+ sol_0_12[:act0Tpb])/(getdefault(tdb.Vpb)*1E3), label = "12 mg", alpha = 0.6, lw = 4, linestyle = :dashdot);
plot!(sol_0_24.t, (sol_0_24[:restTpb] .+ sol_0_24[:actTpb] .+ sol_0_24[:act0Tpb])/(getdefault(tdb.Vpb)*1E3), label = "24 mg", alpha = 0.7, lw = 3, linestyle = :dashdotdot);
plot!(sol_0_48.t, (sol_0_48[:restTpb] .+ sol_0_48[:actTpb] .+ sol_0_48[:act0Tpb])/(getdefault(tdb.Vpb)*1E3), label = "48 mg", alpha = 0.8, lw = 2, linestyle = :dash);
plot!(sol_0_60.t, (sol_0_60[:restTpb] .+ sol_0_60[:actTpb] .+ sol_0_60[:act0Tpb])/(getdefault(tdb.Vpb)*1E3), label = "60 mg", alpha = 0.9, lw = 1, linestyle = :solid);
scatter!(obs_Tcount.time_day, obs_Tcount.Tcell_peruL_median, ms = 6, yerror = (obs_Tcount.Tcell_peruL_median .- obs_Tcount.Tcell_peruL_lowIQR, obs_Tcount.Tcell_peruL_highIQR .- obs_Tcount.Tcell_peruL_median), label = "Hutchings 2021");
plot!(xticks = ([0, 7, 14, 21, 28, 35], string.([0, 7, 14, 21, 28, 35])), xlims = [0, 35]);
display(p_tcell);

## plot B cell 
p_bcell = plot(xlabel = "Time (day)", ylabel = "B cell count (/uL peripheral blood)", legend = :outerright, dpi = 300);
plot!(sol_0_12.t, sol_0_12[:Bpb]/(getdefault(tdb.Vpb)*1E3), label = "12 mg", alpha = 0.6, lw = 4, linestyle = :dashdot);
plot!(sol_0_24.t, sol_0_24[:Bpb]/(getdefault(tdb.Vpb)*1E3), label = "24 mg", alpha = 0.7, lw = 3, linestyle = :dashdotdot);
plot!(sol_0_48.t, sol_0_48[:Bpb]/(getdefault(tdb.Vpb)*1E3), label = "48 mg", alpha = 0.8, lw = 2, linestyle = :dash);
plot!(sol_0_60.t, sol_0_60[:Bpb]/(getdefault(tdb.Vpb)*1E3), label = "60 mg", alpha = 0.9, lw = 1, linestyle = :solid);
plot!(xticks = ([0, 7, 14, 21, 28, 35, 42], string.([0, 7, 14, 21, 28, 35, 42])), xlims = [0, 42]);
display(p_bcell);

## plot IL6
p_il6 = plot(ylims = [0.25, 252], yscale = :log10, xlabel = "Time (day)", ylabel = "plasma IL6 concentration (pg/mL)", legend = :topright, dpi = 300);
plot!(sol_0_12.t, sol_0_12[:IL6pb], label = "12 mg", alpha = 0.6, lw = 4, linestyle = :dashdot);
plot!(sol_0_24.t, sol_0_24[:IL6pb], label = "24 mg", alpha = 0.7, lw = 3, linestyle = :dashdotdot);
plot!(sol_0_48.t, sol_0_48[:IL6pb], label = "48 mg", alpha = 0.8, lw = 2, linestyle = :dash);
plot!(sol_0_60.t, sol_0_60[:IL6pb], label = "60 mg", alpha = 0.9, lw = 1, linestyle = :solid);
scatter!(obs_IL6.time_day, obs_IL6.IL6_pgmL_median, ms = 6, yerror = (obs_IL6.IL6_pgmL_median .- obs_IL6.IL6_pgmL_lowIQR, obs_IL6.IL6_pgmL_highIQR .- obs_IL6.IL6_pgmL_median), label = "Hutchings 2021");
plot!(xticks = ([0, 7, 14, 21, 28, 35, 42], string.([0, 7, 14, 21, 28, 35, 42])), xlims = [0, 42]);
display(p_il6);

## plot PK 
p_epco_plasma = plot(ylims = [0.0002, 20], yscale = :log10, xlabel = "Time (day)", ylabel = "plasma epcoritamab concentration (ug/mL)", legend = :bottomright, dpi = 300);
plot!(sol_0_12.t, sol_0_12[:TDBc_ugperml], label = "12 mg", alpha = 0.6, lw = 4, linestyle = :dashdot);
plot!(sol_0_24.t, sol_0_24[:TDBc_ugperml], label = "24 mg", alpha = 0.7, lw = 3, linestyle = :dashdotdot);
plot!(sol_0_48.t, sol_0_48[:TDBc_ugperml], label = "48 mg", alpha = 0.8, lw = 2, linestyle = :dash);
plot!(sol_0_60.t, sol_0_60[:TDBc_ugperml], label = "60 mg", alpha = 0.9, lw = 1, linestyle = :solid);
plot!(pk_48_bla.time_week * 7, pk_48_bla.conc_ugmL, seriestype = :scatter, label = "48mg, simulated in popPK for BLA filing", alpha = 0.5); 
plot!(xticks = ([0, 7, 14, 21, 28, 35, 42], string.([0, 7, 14, 21, 28, 35, 42])), xlims = [0, 42]);
plot!(yticks = [0.01, 0.1, 1, 10]);
display(p_epco_plasma);

## save figures
savefig(p_epco_plasma, "deliv/figure/invivo/epco-hosseini-plasma-pk.png");
savefig(p_bcell, "deliv/figure/invivo/epco-hosseini-pb-b.png");
savefig(p_tcell, "deliv/figure/invivo/epco-hosseini-pb-t.png");
savefig(p_il6, "deliv/figure/invivo/epco-hosseini-plasma-il6.png");
