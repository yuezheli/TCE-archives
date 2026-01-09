# date: 9/23/24
# author: Yuezhe Li 
# purpose of this code: test IL6 production in glofitamab, using Hosseini model 

using Pkg; Pkg.activate("")

using DifferentialEquations, ModelingToolkit 
using ModelingToolkit: getdefault
using Plots
using CSV, DataFrames, DataFramesMeta

# read in observed data 
pk = CSV.read("data/Frances2021_Fig2A.csv", DataFrame, header=true);
# pk1 = @rsubset(pk, :Dose_ug_kg == 1);
pk3 = @rsubset(pk, :Dose_ug_kg == 3);
pk10 = @rsubset(pk, :Dose_ug_kg == 10);
pk30 = @rsubset(pk, :Dose_ug_kg == 30);
il6 = CSV.read("data/Frances2021_Fig2B.csv", DataFrame, header=true);
il6_3 = @rsubset(il6, :Dose_ug_kg == 3);
il6_10 = @rsubset(il6, :Dose_ug_kg == 10);
il6_30 = @rsubset(il6, :Dose_ug_kg == 30);

# parameter updates for glofitamab
include("model/ParamUpdate.jl");
p_glofit = TDB_param_update("glofitamab");

# original hosseini model 
include("model/tdb_homo.jl");
@mtkbuild tdb = TDB_homo();

# update Hosseini model to glofitamab and cyno PK
const bw_cyno = 3.; # [kg];  Frances et al., 2022; https://jpharmsci.org/article/S0022-3549(21)00697-3/fulltext
prob0 = ODEProblem(tdb, [tdb.TCEinjection_effect => 10., tdb.Bpb => 0., tdb.Btiss => 0, tdb.Btiss2 => 0, tdb.B1920tiss3 => 0], (0., 3.), 
                    [tdb.KmBT_act => p_glofit.KmBT_act, tdb.KdrugactT => p_glofit.KdrugactT, tdb.KdrugB => p_glofit.KdrugB, # glofit in vitro parameters 
                    # glofit cyno PK 
                    tdb.CL_TDB => 500*bw_cyno, tdb.Q_TDB => 100*bw_cyno, tdb.V1_TDB => 34.7*bw_cyno, tdb.V2_TDB => 40*bw_cyno, tdb.Vm_tdb => 0.,
                    # cyno params
                    tdb.kIL6prod => 4., tdb.Vpb => 380, tdb.Vtissue => 7., tdb.Vtissue2 => 25, tdb.Vtissue3 => 50, 
                    tdb.KTrp => 500, tdb.KTrp2 => 500, tdb.KTrp3 => 50, tdb.KBp => 900, tdb.KBp2 => 600, tdb.KBp3 => 60]);

# update initial dosing 
prob_3 = remake(prob0, u0 = Dict([tdb.TDBc_ugperml => 3*bw_cyno/(34.7*bw_cyno) ])); 
prob_10 = remake(prob0, u0 = Dict([tdb.TDBc_ugperml => 10*bw_cyno/(34.7*bw_cyno) ]));
prob_30 = remake(prob0, u0 = Dict([tdb.TDBc_ugperml => 30*bw_cyno/(34.7*bw_cyno) ]));

# simulation 
sol_3 = solve(prob_3, reltol=1e-18, saveat = 0.1); 
sol_10 = solve(prob_10, reltol=1e-18, saveat = 0.1); 
sol_30 = solve(prob_30, reltol=1e-18, saveat = 0.1); 

# pk plot
p_pk = plot(xlabel = "Time (h)", ylabel = "Glofitamab concentration (ug/mL)", legendtitle = "Dose", legendtitlefontsize = 8, legend = :outerright, dpi = 300); 
plot!(sol_3.t * DayToHour, sol_3[:TDBc_ugperml], lw = 3, alpha = 0.5, label = "3 ug/kg", color = :blue);
plot!(sol_10.t * DayToHour, sol_10[:TDBc_ugperml], lw = 3, alpha = 0.5, label = "10 ug/kg", color = :green);
plot!(sol_30.t * DayToHour, sol_30[:TDBc_ugperml], lw = 3, alpha = 0.5, label = "30 ug/kg", color = :grey10);
plot!(pk3.Hour, pk3.glofit_ug_mL, lw = 3, alpha = 0.5, seriestype=:scatter, label = "3 ug/kg; Frances et al., 2022", color = :blue);
plot!(pk10.Hour, pk10.glofit_ug_mL, lw = 3, alpha = 0.5, seriestype=:scatter, label = "10 ug/kg; Frances et al., 2022", color = :green);
plot!(pk30.Hour, pk30.glofit_ug_mL, lw = 3, alpha = 0.5, seriestype=:scatter, label = "30 ug/kg; Frances et al., 2022", color = :grey10);
ylims!(1E-3, 1E1); plot!(yaxis=:log);
display(p_pk);

# IL6 
p_il6 = plot(xlabel = "Time (h)", ylabel = "IL6 concentration (pg/mL)", legendtitle = "Dose", legendtitlefontsize = 8, legend = :outerright, dpi = 300); 
plot!(sol_3.t * DayToHour, sol_3[:IL6pb], lw = 3, alpha = 0.5, label = "3 ug/kg", color = :blue);
plot!(sol_10.t * DayToHour, sol_10[:IL6pb], lw = 3, alpha = 0.5, label = "10 ug/kg", color = :green);
plot!(sol_30.t * DayToHour, sol_30[:IL6pb], lw = 3, alpha = 0.5, label = "30 ug/kg", color = :grey10);
plot!(il6_3.Hour, il6_3.IL6_pg_mL, lw = 3, alpha = 0.5, seriestype=:scatter, label = "3 ug/kg; Frances et al., 2022", color = :blue);
plot!(il6_10.Hour, il6_10.IL6_pg_mL, lw = 3, alpha = 0.5, seriestype=:scatter, label = "10 ug/kg; Frances et al., 2022", color = :green);
plot!(il6_30.Hour, il6_30.IL6_pg_mL, lw = 3, alpha = 0.5, seriestype=:scatter, label = "30 ug/kg; Frances et al., 2022", color = :grey10);
plot!(xlims = (0, 60), ylims = (0.01, 1E4), yaxis = :log10);
display(p_il6); 

# save figure
savefig(p_pk, "deliv/figure/invivo/glofit-cyno-hosseini-plasma-pk.png");
savefig(p_il6, "deliv/figure/invivo/glofit-cyno-hosseini-plasma-il6.png");
