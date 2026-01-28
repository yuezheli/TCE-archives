# date: 5/1/2025
# author: Yuezhe Li 
# purpose of this code: parameter optimization for mouse data 

using ProjectRoot
using Pkg; Pkg.activate(@projectroot)
using DifferentialEquations, ModelingToolkit 
using ModelingToolkit: getdefault
using SymbolicIndexingInterface
using SymbolicIndexingInterface: parameter_values, state_values
using DataFrames
using DataFramesMeta
using Chain
using Statistics
using Plots
using Plots.PlotMeasures
using Printf
using CSV
using Optimization, OptimizationOptimJL
using ForwardDiff
using Random
using JLD2

#--------------- base model ---------------#
include(@projectroot("model/model_constants.jl"));
include(@projectroot("model/invivo/mouse_pkpd_v8_2.jl"));
include(@projectroot("script/invivo/helper-mouse-pkpd-visualization.jl"));
include(@projectroot("script/invivo/helper-mouse-pkpd-simulation.jl"));

@mtkbuild pkpd = Mouse_PKPD(); #

#--------------- obs ---------------#
obs_mousePKPD = CSV.read(@projectroot("data/derived/mouse.csv"), DataFrame, header=true, missingstring=".");
replace!(obs_mousePKPD.DOSE, -99=>0); #-99 stands for 'vehicle' or 'control' or dose = 0

# Tumor volume, 2.2
obs_mouseTV_2p2 = @rsubset(obs_mousePKPD, :DVID == 0, :SOURCE == "2.2"); # DVID == 0 for tumor volume 

# Tumor volume, 2.1
obs_mouseTV_2p1_singledose = @rsubset(obs_mousePKPD, :DVID == 0, :SOURCE == "2.1", :TRT in [0, 0.5, 1,2,3]);  # filter for data in 2.1 of single dose
obs_mouseTV_2p1_multidose = @rsubset(obs_mousePKPD, :DVID == 0, :SOURCE == "2.1", :TRT in [0, 0.5, 4,5]);  # filter for data in 2.1 of multiple dose

# CD8 (from 2.2)
obs_mus_cd8_obs = @rsubset(obs_mousePKPD, :SOURCE == "2.2", :DVID == 3);

function ComputeE2T(obs)
    obs_mus_cd8 = @select(obs, :DOSE, :DAY, :DV, :TRT, :DVID); 
    # conversion based on Mi et al., 2020; https://www.frontiersin.org/journals/physiology/articles/10.3389/fphys.2020.583333/full
    @transform!(obs_mus_cd8, :cd8_per_mm3 = :DV / 0.011); # convert cell/mm2 to cell/mm3, t=5um, D=7um => d = 11um (Acrit=10um2) 
    @transform!(obs_mus_cd8, :DV = :cd8_per_mm3*mm3_per_L*cell_vol/getdefault(pkpd.f_tum_cell) / 5); # fudge factor of 5 was applied to make the E2T ratio close to BI data 
    @select!(obs_mus_cd8, :DOSE, :DAY, :TRT, :DV, :DVID)
    return obs_mus_cd8
end

obs_mus_cd8 = ComputeE2T(obs_mus_cd8_obs)

# compute median
function dv_median(obs)
    df_summary = @chain obs begin
        @groupby :DOSE :DAY :TRT :DVID
        @combine :med_dv = quantile(skipmissing(:DV), 0.5)
        @select :DOSE :DAY :DV_med = :med_dv :TRT :DVID
    end;
    return df_summary
end

# simulation
tspan_hr = (0, 25*hr_per_day);  # [hr]

prob0 = ODEProblem(pkpd, [], tspan_hr); 

# read in parameters from in vitro 
param_hpac_v8_il10 = CSV.read(@projectroot("data/param/invitro-hpac-v8-il10.csv"), DataFrame, header=true, missingstring=".");
param_hpac_v8_tnfa = CSV.read(@projectroot("data/param/invitro-hpac-v8-tnfa-1.csv"), DataFrame, header=true, missingstring=".");
param_hpac_v8_ifng = CSV.read(@projectroot("data/param/invitro-hpac-v8-ifng-1.csv"), DataFrame, header=true, missingstring=".");

# read in slow activating donor
param_donor2_1p6 = CSV.read(@projectroot("data/param/invitro-hpac-v8-1point6-donor-2.csv"), DataFrame, header=true, missingstring=".");
param_donor1_1p6 = CSV.read(@projectroot("data/param/invitro-hpac-v8-1point6-donor-1.csv"), DataFrame, header=true, missingstring=".");

# update base parameters 
p_00 = Dict([pkpd.tumor_T_cells_init => 0, 
            pkpd.Trpbref_perml => 1,  # assumed
            pkpd.f_T_cell_survival => 0.2, # assumed
            pkpd.kTgen => 0,  # turn off T cell proliferation outside tumor 
            pkpd.k_t_n_prolif => 0, # turn off naive T cell proliferation in tumor
            pkpd.TAA_per_cell => TAA_HPAC, 
            pkpd.lambda_taa_avidity => @rsubset(param_hpac_v8_il10, :Name == "lambda_taa_avidity").Value[1], 
            # update params 
            pkpd.emax_t_activate => @rsubset(param_hpac_v8_il10, :Name == "emax_t_activate").Value[1]*s_per_hr, 
            pkpd.emax_apop => @rsubset(param_hpac_v8_il10, :Name == "emax_apop").Value[1]*s_per_hr, 
            pkpd.k_syn_cytokine1 => @rsubset(param_hpac_v8_il10, :Name == "k_syn_cytokine").Value[1]*s_per_hr,  # set cytokine 1 to be IL10
            pkpd.k_syn_cytokine2 => @rsubset(param_hpac_v8_tnfa, :Name == "k_syn_cytokine").Value[1]*s_per_hr,  # set cytokine 2 to be TNFa
            pkpd.k_syn_cytokine3 => @rsubset(param_hpac_v8_ifng, :Name == "k_syn_cytokine").Value[1]*s_per_hr,  # set cytokine 3 to be IFNg
            pkpd.tt50_t_activate => @rsubset(param_hpac_v8_il10, :Name == "tt50_t_activate").Value[1],
            pkpd.tt50_kill => @rsubset(param_hpac_v8_il10, :Name == "tt50_kill").Value[1],
            pkpd.tt50_cytokine => @rsubset(param_hpac_v8_il10, :Name == "tt50_cytokine").Value[1],
            # cytokine degredation rate calculated from their serum half life
            pkpd.k_deg_cytokine1 => 0.25, # https://pmc.ncbi.nlm.nih.gov/articles/PMC4454631/
            pkpd.k_deg_cytokine2 => 0.25, # https://pmc.ncbi.nlm.nih.gov/articles/PMC7887858/, https://advanced.onlinelibrary.wiley.com/doi/full/10.1002/advs.202004433
            pkpd.k_deg_cytokine3 => 0.25, # https://pmc.ncbi.nlm.nih.gov/articles/PMC7887858/
            # PK 
            pkpd.CL_TCE => 0.12/hr_per_day/mL_per_L, 
            pkpd.Q_TCE => 0.5/hr_per_day/mL_per_L, 
            # turn on TCE driven T cell infiltration into tumor
            pkpd.FLG_tce_driven_T_cell_enter_tumor => 1, 
            ]); 

prob1 = remake(prob0, p = p_00, u0 = Dict());

function sim_2p1_2p2(; prob_ = prob1, newparam = nothing)
    mouse_2p2_2 = SimulateMultiDose(prob_, dose_mgkg_list = [0, 0.1], parameters = merge(Dict(pkpd.TV_mm3_0 => 200.), newparam) );

    mouse_2p1_singldose = simulateDose(prob_, [0, 0.01, 0.1, 1], tspan_hr, doseday = 4, parameters = newparam); 

    mouse_2p1_multdose_control = simulateDose(prob_, [0], tspan_hr, doseday = 4, parameters = newparam); 
    mouse_2p1_multdose_point01 = SimulateMultiDose(prob_, dose_mgkg_list = [0.01], parameters = newparam, doseday = [4, 18] );
    mouse_2p1_multdose_point1 = SimulateMultiDose(prob_, dose_mgkg_list = [0.1], parameters = newparam, doseday = [4,7,11,14,18,21,25] );
    mouse_2p1_multdose = merge(mouse_2p1_multdose_control, mouse_2p1_multdose_point01, mouse_2p1_multdose_point1); 

    return mouse_2p1_singldose, mouse_2p1_multdose, mouse_2p2_2
end

function ProcessSims(mus = mouse_2p2_2; pkpd=pkpd, dvid = 3)
    df = DataFrame(DAY = [], DV = [], DVID = [], DOSE = [])
    for ks in keys(mus)
        tmpsims = mus[ks]
        tmpdays = tmpsims.t ./ hr_per_day
        tmptv_mm3 = tmpsims[pkpd.TV_mm3]
        tmpe2t = (tmpsims[pkpd.actTtumor] .+ tmpsims[pkpd.restTtumor]) ./ (tmpsims[pkpd.Ntot])
        tmpdf1 = DataFrame(DAY = tmpdays, DV = tmptv_mm3, DVID = 0, DOSE = ks); # tumor volume; dvid = 0
        tmpdf2 = DataFrame(DAY = tmpdays, DV = tmpe2t, DVID = dvid, DOSE = ks);
        df = vcat(df, tmpdf1, tmpdf2)
    end
    return df
end

function loss(p, obs = nothing; prob_ = prob1, pkpd = pkpd, 
              obs_mus_cd8 = obs_mus_cd8, obs_mouseTV_2p2 = obs_mouseTV_2p2, 
              obs_mouseTV_2p1_singledose = obs_mouseTV_2p1_singledose, obs_mouseTV_2p1_multidose = obs_mouseTV_2p1_multidose, 
              cd8weight = 1, tvweight = 1)
              println(p)

    p_2 = Dict([pkpd.q_t_in => p[1],  pkpd.kTrexit => p[2], pkpd.KTrptumor => p[3], 
                pkpd.fTaprolif => p[4], pkpd.frac_act_T_death => p[5], pkpd.emax_apop => p[6]])

    mouse_2p1_singldose, mouse_2p1_multdose, mouse_2p2_2 = sim_2p1_2p2(prob_ = prob_, newparam = p_2);
    
    comb_2p2_cd8 = innerjoin( dv_median(obs_mus_cd8), ProcessSims(mouse_2p2_2), on = [:DAY, :DVID, :DOSE])
    comb_2p2_tv = innerjoin( dv_median(obs_mouseTV_2p2), ProcessSims(mouse_2p2_2), on = [:DAY, :DVID, :DOSE])

    comb_2p1_singledose_tv = innerjoin( dv_median(obs_mouseTV_2p1_singledose), ProcessSims(mouse_2p1_singldose), on = [:DAY, :DVID, :DOSE])
    comb_2p1_multidose_tv = innerjoin( dv_median(obs_mouseTV_2p1_multidose), ProcessSims(mouse_2p1_multdose), on = [:DAY, :DVID, :DOSE])

    tmpdiff_2p2_cd8 = sum( (comb_2p2_cd8.DV .- comb_2p2_cd8.DV_med).^2 )
    tmpdiff_2p2_tv = sum( (comb_2p2_tv.DV .- comb_2p2_tv.DV_med).^2 )

    tmpdiff_2p1_tv = sum( (comb_2p1_singledose_tv.DV./median(comb_2p1_singledose_tv.DV) .- comb_2p1_singledose_tv.DV_med./median(comb_2p1_singledose_tv.DV_med)).^2 ) + sum( (comb_2p1_multidose_tv.DV./median(comb_2p1_multidose_tv.DV) .- comb_2p1_multidose_tv.DV_med./median(comb_2p1_multidose_tv.DV_med)).^2 )

    diff = cd8weight*tmpdiff_2p2_cd8 + tvweight*(tmpdiff_2p2_tv + tmpdiff_2p1_tv)
    return diff

end

# note the optimization was set up as such bc of the PSO limitation
p_pre = [0.8, 4e-4, 1, 0.02, 0.01, @rsubset(param_hpac_v8_il10, :Name == "emax_apop").Value[1]*s_per_hr * 0.5]

loss(p_pre)

function PlotOpt(p; prob_ = prob1, obs_mouseTV_2p1_singledose = obs_mouseTV_2p1_singledose, obs_mouseTV_2p1_multidose = obs_mouseTV_2p1_multidose, 
                 obs_mouseTV_2p2 = obs_mouseTV_2p2, obs_mus_cd8 = obs_mus_cd8, pkpd = pkpd)
    
    p_dict = Dict([pkpd.q_t_in => p[1],  pkpd.kTrexit => p[2], pkpd.KTrptumor => p[3], 
                pkpd.fTaprolif => p[4], pkpd.frac_act_T_death => p[5], pkpd.emax_apop => p[6]])

    mus_2p1_single, mus_2p1_multi, mus_2p2 = sim_2p1_2p2(; prob_ = prob_, newparam = p_dict)

    p_comb = plot(
    plot_TV_sim([0, 0.1], mus_2p2, obs_mouseTV = obs_mouseTV_2p2, yaxis_log = true, legendpos = :bottomleft, yrange = [10, 1E4]),  
    Plot_E2T_calculatedE2T([0., 0.1], mus_2p2, @transform!(obs_mus_cd8, :e2t = :DV); xrange = missing, legendpos = :topleft),
    plot_TV_sim([0, 0.01, 0.1, 1], mus_2p1_single, obs_mouseTV = obs_mouseTV_2p1_singledose, yaxis_log = true, legendpos = :topleft, yrange = [10, 1E4]), 
    plot_TV_sim([0, 0.01, 0.1], mus_2p1_multi, obs_mouseTV = obs_mouseTV_2p1_multidose, yaxis_log = true, legendpos = :topleft, yrange = [10, 1E4]), 
    layout = (2,2), size = (800, 700), margin = 30px
    );
    return p_comb
end

# q_t_in, kTrexit, KTrptumor, fTaprolif, frac_act_T_death, emax_apop
p_lower_bound = [0, 1E-4, 0.1, 1E-4, 0, 0];
p_upper_bound = [5, 1, 10, 0.1, 1, @rsubset(param_hpac_v8_il10, :Name == "emax_apop").Value[1]*s_per_hr];

p_lower_bound = p_pre / 100
p_upper_bound = p_pre * 3
p_upper_bound[end] = @rsubset(param_hpac_v8_il10, :Name == "emax_apop").Value[1]*s_per_hr

using StaticArrays
p_lower_bound = SVector{length(p_lower_bound)}(p_lower_bound)
p_upper_bound = SVector{length(p_upper_bound)}(p_upper_bound)
p_pre = SVector{length(p_pre)}(p_pre)


using ParallelParticleSwarms 
using KernelAbstractions
optfn = OptimizationFunction(loss)
optprob = OptimizationProblem(optfn,p_pre, lb = p_lower_bound, ub = p_upper_bound)

p_optim_gl = Optimization.solve(optprob, ParallelSyncPSOKernel(100, backend = CPU(static=false)),maxiters=100)

p_optim_gl.stats


function SaveOutcome(p_optim_gl)
    p_tosave = @chain DataFrame(q_t_in = p_optim_gl.u[1], kTrexit = p_optim_gl.u[2], KTrptumor = p_optim_gl.u[3], fTaprolif = p_optim_gl.u[4], frac_act_T_death = p_optim_gl.u[5], emax_apop = p_optim_gl.u[6]) begin
            stack(_)
            rename!(_, :variable => :Name, :value => :Value)
    end
    return p_tosave
end

CSV.write(@projectroot("data/param/mus-pso-trial-2.csv"), SaveOutcome(p_optim_gl))

# save all solver outcome
# JLD2.save(@projectroot("data/param/mouse-parallel-pso-outcome.jld2"),Dict("sol"=>p_optim_gl.u, "stats" =>p_optim_gl.stats, "retcode"=>p_optim_gl.retcode, "objective"=>p_optim_gl.objective, "alg" =>p_optim_gl.alg))