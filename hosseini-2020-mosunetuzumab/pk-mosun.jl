# date: 8/30/24
# author: Yuezhe Li 
# purpose of this code: to test the TDB model in MTK

using Pkg; Pkg.activate("")

using DifferentialEquations, ModelingToolkit 
using DataFrames, CSV
using Plots
using ModelingToolkit: getdefault

# read in clincial PK
mosu2 = CSV.read("lunsumio-pk-genetech.csv",DataFrame);

# build model 
include("TDB_homo.jl");
@mtkbuild tdb = TDB_homo();
prob = ODEProblem(tdb, [], (-1E-12, 80.));

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

# dosing event 
dose_amt = [1., 2., 60., 60., 30.] .* 1e3;  # [ug]
dose_time = [0, 7, 14, 21, 42];                 # [day] 
cbset_human = tdb_iv_dosing(dose_amt, dose_time, getdefault(tdb.V1_TDB)); 

# simulation 
sdf2 = solve(prob, reltol=1e-18, saveat = 0.2, callback = cbset_human); 

# visualization
p_pk = plot(title = "Mosunetuzumab PK (dose = 1, 2, 60, 60, 30mg)", titlefont = 8, legend = :outerright, xticks = ([0, 7, 14, 21, 42, 63], string.([0, 7, 14, 21, 42, 63])), dpi = 300);
scatter!(mosu2.day, mosu2.mosun_ug_mL, label = "Lunsumio data", alpha = 0.5, color = :blue, markerstrokewidth=0, markersize=5);
plot!(sdf2.t, sdf2[:TDBc_ugperml], label = "sims", color = :blue);
xlabel!("time (day)", xguidefontsize = 8); 
ylabel!("Mosunetuzumab plasma conc (ug/mL)", yguidefontsize = 8);
plot!(yscale = :log10); 
ylims!(0.01, 100);

savefig(p_pk, "pk-lunsumio.png")

# plot(sdf2.t, sdf2[:B19no20tiss3], label = false)

# additional digging around steady state 
sol0 = solve(prob, reltol=1e-18, saveat = 0.2); 

# plot(sol0.t, sol0[:B1920tiss3])
# B cell output from BM 
getdefault(tdb.tissue3on)*getdefault(tdb.kBtiss3exit)*getdefault(tdb.Vpb)*max.(0,sol0[:B1920tiss3]/getdefault(tdb.Vtissue3)-sol0[:Bpb]/getdefault(tdb.Vpb)*getdefault(tdb.KBp3))

