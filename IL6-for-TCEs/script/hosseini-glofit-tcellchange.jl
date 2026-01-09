# date: 9/23/24
# author: Yuezhe Li 
# purpose of this code: simulate T cell change in human after glofit dosing in the Hosseini model 

using Pkg; Pkg.activate("")

using DifferentialEquations, ModelingToolkit 
using ModelingToolkit: getdefault
using Plots
using CSV, DataFrames, DataFramesMeta

# data obtained from Hutchings et al., 2021, Fig S7; https://pubmed.ncbi.nlm.nih.gov/33739857/
obs_tpb = data = CSV.read("data/Hutchings2021_JClinOncol_Tchange.csv", DataFrame, header=true);
obs_16 = @rsubset(obs_tpb, :Dose_mg == "2.5_10_16");
obs_30 = @rsubset(obs_tpb, :Dose_mg == "2.5_10_30");

# parameter updates for glofitamab
include("model/ParamUpdate.jl");
p_glofit = TDB_param_update("glofitamab");

# original hosseini model 
include("model/tdb_homo.jl");
@mtkbuild tdb = TDB_homo();

# base Hosseini model for glofitamab (tumor left off)
prob0 = ODEProblem(tdb, [tdb.TCEinjection_effect => 10., tdb.Bpb => 0., tdb.Btiss => 0, tdb.Btiss2 => 0, tdb.B1920tiss3 => 0], (0., 28.), 
                [tdb.CL_TDB => p_glofit.CL_TDB, tdb.Q_TDB => p_glofit.Q_TDB, tdb.V1_TDB => p_glofit.V1_TDB, tdb.V2_TDB => p_glofit.V2_TDB, tdb.Vm_tdb => 0.,
                tdb.KmBT_act => p_glofit.KmBT_act, tdb.KdrugactT => p_glofit.KdrugactT, tdb.KdrugB => p_glofit.KdrugB]);

# dosing function 
function tdb_iv_dosing(dose_amt, dose_time, V1_TDB; injection_effect_init = 10.)
    global cbs_human = [];
    if length(dose_time) > 1
        for i in 1:length(dose_time)
            function affect!(integrator)
                integrator[:TDBc_ugperml] += dose_amt[i]/ V1_TDB;
                integrator[:TCEinjection_effect] += injection_effect_init;
            end
            cb = PresetTimeCallback(dose_time[i],affect!);
            global cbs_human = push!(cbs_human, cb);
        end
    end
    cbset_human = CallbackSet(cbs_human...);
    return cbset_human
end

dose_amt_30 = [2.5, 10., 30.] .* 1e3;  # [ug]
dose_amt_16 = [2.5, 10., 16.] .* 1e3;  # [ug]
dose_time = [7, 14, 21];                 # [day] 
cbset_30 = tdb_iv_dosing(dose_amt_30, dose_time, p_glofit.V1_TDB); 
cbset_16 = tdb_iv_dosing(dose_amt_16, dose_time, p_glofit.V1_TDB); 

# simulation 
sol_30 = solve(prob0, reltol=1e-18, saveat = 0.1, callback = cbset_30); 
sol_16 = solve(prob0, reltol=1e-18, saveat = 0.1, callback = cbset_16); 

function pbTchange(sol_glofit)
    sdf_glofit = DataFrame(sol_glofit); 
    @select!(sdf_glofit, :Time = :timestamp, :actTpb, :restTpb, :act0Tpb);
    @transform!(sdf_glofit, :totalTpb = :actTpb .+ :restTpb .+ :act0Tpb)
    @transform!(sdf_glofit, :log2_Tpb_foc = log2.( :totalTpb / first(:totalTpb) ) );
    return sdf_glofit;
end

df16 = pbTchange(sol_16);
df30 = pbTchange(sol_30);

p_Tpb_change = plot(xlabel = "Time (day)", ylabel = "log2 Fold-change T cells in PB", yticks = (-6:2:6, string.(-6:2:6)), dpi = 300); 
plot!(df16.Time, df16.log2_Tpb_foc, lw = 3, alpha = 0.5, label = "2.5/10/16mg", color = :red);
scatter!(obs_16.Time, obs_16.log2_Tpb_FoC, ms=5, ma=0.5, markerstrokewidth=0, label = "2.5/10/16mg; Hutchings et al., 2021", color = :red);
plot!(df30.Time, df30.log2_Tpb_foc, lw = 2, linestyle = :dash, alpha = 0.5, label = "2.5/10/30mg", color = :blue);
scatter!(obs_30.Time, obs_30.log2_Tpb_FoC, ms=5, ma=0.3, markerstrokewidth=0, label = "2.5/10/30mg; Hutchings et al., 2021", color = :blue);
xlims!(0, 21);
plot!(legend = :topleft);

savefig(p_Tpb_change, "deliv/figure/invivo/glofit-homo-hosseini-pb-t-change.png");

