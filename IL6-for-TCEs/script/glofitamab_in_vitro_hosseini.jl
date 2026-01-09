# date: 9/24/24
# author: Yuezhe Li 
# purpose of this code: to produce glofit in vitro fit

using Pkg; Pkg.activate("")

using DifferentialEquations, ModelingToolkit 
using ModelingToolkit: getdefault
using Plots
using DataFrames

# parameter updates for epcoritamab
include("model/ParamUpdate.jl");
p_glofit = TDB_param_update("glofitamab");

# in vitro params
p_glofitamab = ComponentArray(
    kBtumorprolif = p_glofit.kBtumorprolif, KmBT_act = p_glofit.KmBT_act, S = p_glofit.S, KdrugactT = p_glofit.KdrugactT, ndrugactT = p_glofit.ndrugactT, VmT = p_glofit.VmT, 
    KmTB_kill = p_glofit.KmTB_kill, nkill = p_glofit.nkill, KdrugB = p_glofit.KdrugB, VmB = p_glofit.VmB, kTact = p_glofit.kTact, fTadeact = p_glofit.fTadeact, fTa0deact = p_glofit.fTa0deact); 

# in vitro data from Bacac et al., 2018, Fig 1B; # https://pubmed.ncbi.nlm.nih.gov/29716920/
tcb_conc = [0.001, 0.01, 0.1, 1., 10., 100., 1000.]; # [pmol/L]
toledo_monocd20 = DataFrame(TCBconc = tcb_conc, tumor_lysis = [-3.7, -3.7, -3.7, -3.7, -0.85, 32.32, 46.65]); # this is for mosun
toledo_glofit = DataFrame(TCBconc = tcb_conc, tumor_lysis = [-1.14, 3.49, 14.89, 29.96, 41.51, 47.16, 49.45]);

# load hosseini in vitro model 
const MW_EDG = 150000.0; # Da
include("model/tdb-in-vitro-hosseini.jl");

# init condition 
u0_0 = ComponentArray(Btumor = 1.0E3, bs_ugperml = 0., restT = 5.0E3, actT = 0., act0T = 0.); 

# simulation 
test_glofit = [];
for conc in tcb_conc
    u0_tmp = deepcopy(u0_0);
    u0_tmp.bs_ugperml = conc/1E9 * MW_EDG; # [ug/mL]
    sol_tmp = solve(ODEProblem(bs_in_vitro_hosseini!, u0_tmp, (0., 1.), p_glofitamab), alg=AutoTsit5(Rosenbrock23()), reltol=1e-18, saveat = 0.2);
    tmp_kill_ratio = 1- sol_tmp.u[end].Btumor/u0_tmp.Btumor
    append!(test_glofit, tmp_kill_ratio)
end

# visualization
p_invitro_glofit = plot(xlabel = "Glofitamab concentration (pmol/L)", ylabel = "% Tumor cell lysis",  legend = :outerright, dpi = 300);
plot!(tcb_conc, test_glofit * 100, label = "Simulation", alpha = 0.6, lw = 3);
plot!(tcb_conc, toledo_glofit.tumor_lysis, label = "Bacac et al., 2018", seriestype = :scatter, markersize = 5, alpha = 0.6);
plot!(xaxis=:log10);
display(p_invitro_glofit);

savefig(p_invitro_glofit, "deliv/figure/invitro/glofit-hosseini-in-vitro.png");
