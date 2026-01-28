# date: 4/8/2025
# author: Yuezhe Li 
# purpose of this code: to run simulation for mouse pkpd model 

using OrdinaryDiffEq
using DifferentialEquations

function getVariableIndex(integrator, sym::Symbol)
    return ModelingToolkit.variable_index(integrator.f.sys, sym)
end

# use this function for single dose simulation 
function simulateDose(prob, dose_mgkg_list, tspan_hr; saveat_time = 1.0, parameters = nothing, doseday = 0)

    sol_PKPD = Dict();

    for dose_i in dose_mgkg_list

        if doseday == 0.
            # update parameter values and initial conditions that depend on updated parameter values
            if !isnothing(parameters)
                prob_new = remake(prob; p = vcat(parameters..., [pkpd.dose_mgkg => dose_i]), u0 = Dict(), tspan = tspan_hr)
            else
                prob_new = remake(prob; p = [pkpd.dose_mgkg => dose_i], u0 = Dict(), tspan = tspan_hr)
            end

            # Simulate model with updated dose and initial conditions
            sol_new = solve(prob_new, alg=Rodas5P(), abstol=1e-6, reltol=1e-3, saveat=saveat_time)

            sol_PKPD[dose_i] = sol_new
        
        elseif doseday > 0 
            if !isnothing(parameters)
                prob_new = remake(prob; p = vcat(parameters...), u0 = Dict(), tspan = tspan_hr)
            else
                prob_new = remake(prob; u0 = Dict(), tspan = tspan_hr)
            end

            affect!(integrator) = integrator.u[getVariableIndex(integrator, :TCEc_nM)] += dose_i * integrator.ps[:BW] / (mg_per_g) / (MW_TCE) * (nmol_per_mol) / integrator.ps[:V1_TCE]
            cb = PresetTimeCallback(doseday*hr_per_day, affect!)

            sol_new = solve(prob_new, alg=Rodas5P(), abstol=1e-6, reltol=1e-3, saveat=saveat_time, callback = cb)

            sol_PKPD[dose_i] = sol_new

        else
            println("Dosing time has issue")
        end
    end

    return sol_PKPD;
   
end

# use this function for multi-dose simulation; default dosing days were set to be doing days in dataset 2.2
function SimulateMultiDose(prob; dose_mgkg_list = [0.], saveat_time = 1.0, parameters = nothing, doseday = [4,7,11,18,21,25])

    if length(doseday)< 2
        println("DO NOT USE THIS FUNCTION")
        return nothing 
    else
        if !isnothing(parameters)
            prob__0 = remake(prob, p = parameters, u0 = Dict())
        else 
            prob__0 = prob
        end

        sol_PKPD = Dict();

        for dose_mgkg in dose_mgkg_list
            affect(integrator) = integrator.u[getVariableIndex(integrator, :TCEc_nM)] += dose_mgkg * integrator.ps[:BW] / (mg_per_g) / (MW_TCE) * (nmol_per_mol) / integrator.ps[:V1_TCE]

            if doseday[1] > 0
                cb = PresetTimeCallback(doseday .* hr_per_day, affect)
                sol__ = solve(prob__0, alg=Rodas5P(), abstol=1e-6, reltol=1e-3, saveat=saveat_time, callback = cb)
                sol_PKPD[dose_mgkg] = sol__
            elseif doseday[1] == 0
                prob__1 = remake(prob__0; p = [pkpd.dose_mgkg => dose_mgkg], u0 = Dict())
                # affect(integrator) = integrator.u[getVariableIndex(integrator, :TCEc_nM)] += dose_mgkg * integrator.ps[:BW] / (mg_per_g) / (MW_TCE) * (nmol_per_mol) / integrator.ps[:V1_TCE]
                cb = PresetTimeCallback(doseday[2:end]*hr_per_day, affect)
                sol__ = solve(prob__1, alg=Rodas5P(), abstol=1e-6, reltol=1e-3, saveat=saveat_time, callback = cb)
                sol_PKPD[dose_mgkg] = sol__
            else
                println("Dosing time has issue")
            end
        end
        return sol_PKPD
    end
end

