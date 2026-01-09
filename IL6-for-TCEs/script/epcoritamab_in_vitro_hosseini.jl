# date: 9/24/24
# author: Yuezhe Li 
# purpose of this code: to produce epcoritamab in vitro fit

using Pkg; Pkg.activate("")

using DifferentialEquations 
using Plots
using DataFrames

# parameter updates for epcoritamab
include("model/ParamUpdate.jl");
p_epco = TDB_param_update("epcoritamab");

# in vitro data from Engelberts et al., 2020, Fig 1F; # https://pubmed.ncbi.nlm.nih.gov/31981978/
epco_conc_daudi = [3.3015E-06, 3.31284E-05, 0.000337174, 0.003052746, 0.028595947, 0.262800248, 2.578298474, 23.79369186, 249.6927377]; # [ng/mL]
daudi_epco = DataFrame(TCBconc = epco_conc_daudi, tumor_lysis = [2.957686583, -8.456304486, -5.728258558, 5.277482458, 32.24431214, 67.08909207, 84.96704231, 88.09483309, 88.19795875]);

epco_conc_RI1 = [2.88306E-06, 3.24301E-05, 3.29798E-04, 2.88306E-03, 3.41073E-02, 3.13581E-01, 3.70974, 34.1073, 298.162]; # [ng/mL]
RI1_epco = DataFrame(TCBconc = epco_conc_RI1, tumor_lysis = [3.12269, 9.23908E-01, 3.25167, 5.17748, 18.1812, 50.1040, 60.6419, 67.4963, 67.7782]);

# load hosseini in vitro model 
const MW_EDG = 150000.0; # Da
include("model/tdb-in-vitro-hosseini.jl");

# define initial condition 
u0_0 = ComponentArray(Btumor = 1.0E3, bs_ugperml = 0., restT = 5.0E3, actT = 0., act0T = 0.); 

# define epcoritamab in vitro parameters 
p_epcoritamab = ComponentArray(
    kBtumorprolif = p_epco.kBtumorprolif, KmBT_act = p_epco.KmBT_act, S = p_epco.S, KdrugactT = p_epco.KdrugactT, ndrugactT = p_epco.ndrugactT, VmT = p_epco.VmT, 
    KmTB_kill = p_epco.KmTB_kill, nkill = p_epco.nkill, KdrugB = p_epco.KdrugB, VmB = p_epco.VmB, kTact = p_epco.kTact, fTadeact = p_epco.fTadeact, fTa0deact = p_epco.fTa0deact); 

# simulation 
test_epco = [];
for conc in epco_conc_RI1
    sol0 = solve(ODEProblem(bs_in_vitro_hosseini!, u0_0, (0., 2.), p_epcoritamab), alg=AutoTsit5(Rosenbrock23()), reltol=1e-18, saveat = 0.2);
    u0_tmp = deepcopy(u0_0);
    u0_tmp.Btumor = 2.5E3;
    u0_tmp.restT = 5E3;
    u0_tmp.bs_ugperml = conc/1E3; # [ug/mL]
    sol_tmp = solve(ODEProblem(bs_in_vitro_hosseini!, u0_tmp, (0., 2.), p_epcoritamab), alg=AutoTsit5(Rosenbrock23()), reltol=1e-18, saveat = 0.2);
    #tmp_kill_ratio = 1- sol_tmp.u[end].Btumor/sol0.u[end].Btumor
    tmp_kill_ratio = 1- sol_tmp.u[end].Btumor/u0_tmp.Btumor
    append!(test_epco, tmp_kill_ratio)
end

p_epco_in_vitro = plot(xlabel = "Epcoritamab conc (pmol/L)", ylabel = "% Tumor cell lysis",  legend = :outerright, dpi = 300);
plot!(epco_conc_RI1*1E6/MW_EDG, test_epco * 100, label = "Simulation", alpha = 0.8);
plot!(epco_conc_RI1*1E6/MW_EDG, RI1_epco.tumor_lysis, label = "Engelberts et al., 2020", seriestype = :scatter, alpha = 0.7, markersize = 6);
plot!(xaxis=:log10);
display(p_epco_in_vitro)

savefig(p_epco_in_vitro, "deliv/figure/invitro/epco-hosseini-in-vitro.png");
