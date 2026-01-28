# date: 4/8/2025
# author: Yuezhe Li 
# purpose of this code: to host helper visualization functions

#--------------- visualization, PK ---------------#
function plot_PK_sim(dose_mgkg_list, pksim, obs_mousePK)
    
    usepalette = cgrad(:glasgowS)

    plot_TCE_plasma = plot(size = (400, 400), dpi = 300, background_color_legend = nothing);
    
    for dose_i in dose_mgkg_list
        plot_TCE_plasma = plot!(pksim[dose_i].t./hr_per_day, pksim[dose_i][pkpd.TCE_plasma_ngpermL], yaxis=:log, label="sim $(round(dose_i, digits=2)) mgkg", width=2, linestyle=:dashdot, color = usepalette[dose_i]);
        plot!((@rsubset(obs_mousePK, :DOSE == dose_i)).DAY, (@rsubset(obs_mousePK, :DOSE == dose_i)).DV, seriestype = :scatter, markershape=:circle, markerstrokewidth=0, alpha = 0.5, color = usepalette[dose_i], label="data $(round(dose_i, digits=2)) mgkg");
    end

    Plots.xlabel!("Time (d)");
    Plots.ylabel!("plasma TCE concentration (ng/mL)");
    ylims!(1e0, 1e6);
    yticks!([1e1, 1e2, 1e3, 1e4, 1e5])

    return plot_TCE_plasma

end

# PK diagnostic plot (plot TCE in central and in peripheral)
function plot_TCE_(dose_mgkg_list, sims_outcome; legendpos = :topright, showPeripheral = false, showCentral = false, showTumor = false, 
    yrange = [1E-2, 2E3], legendFontSiz = 6, maxLegendRow = 2)
    usepalette = cgrad(:glasgowS, length(dose_mgkg_list))
    yticks_ticks = 10 .^ (floor(log10.(yrange)[1]):1:floor(log10.(yrange)[2]))
    plot_TCEt= plot(ylims = yrange, yticks = yticks_ticks, yaxis=:log, background_color_legend = nothing, xlabel = "Time (Day)", ylabel = "TCE (nM)", dpi = 300, size = (400, 400));
    if showPeripheral
        for dose_i in dose_mgkg_list
            plot!(sims_outcome[dose_i].t/hr_per_day, sims_outcome[dose_i][pkpd.TCEp_nM], label= "Peripheral $(round(dose_i, digits=2)) mgkg", lw=2, linestyle=:dash, color = usepalette[dose_i]);
        end
    end
    if showCentral
        for dose_i in dose_mgkg_list
            plot!(sims_outcome[dose_i].t/hr_per_day, sims_outcome[dose_i][pkpd.TCEc_nM], label= "Plasma $(round(dose_i, digits=2)) mgkg", lw=2, linestyle=:dashdot, color = usepalette[dose_i]);
        end
    end
    if showTumor
        for dose_i in dose_mgkg_list
            plot!(sims_outcome[dose_i].t/hr_per_day, sims_outcome[dose_i][pkpd.TCEt_nM], label= "TME $(round(dose_i, digits=2)) mgkg", lw=2, linestyle=:solid, color = usepalette[dose_i]);
        end
    end
    plot!(xlims = (0., tspan_hr[2]/hr_per_day + 1), legend = legendpos, legendfontsize = legendFontSiz, legend_columns = maxLegendRow); 
    return plot_TCEt
end

#--------------- visualization function TV  ---------------#
function plot_TV_sim(dose_mgkg_list, pksim; obs_mouseTV = missing, legendpos = :outerright, yaxis_log = false, yrange = missing)
    usepalette = cgrad(:glasgowS, length(dose_mgkg_list))
    plot_TV = plot(background_color_legend = nothing, xlabel = "Time (day)", ylabel = "Tumor Volume (mm^3)", legend = legendpos, legend_columns=2, dpi = 300, size = (400, 400));
    for dose_i in dose_mgkg_list
        plot_TV = plot!(pksim[dose_i].t./hr_per_day, pksim[dose_i][pkpd.TV_mm3], label="sim $(round(dose_i, digits=2)) mgkg", width=2, linestyle=:dashdot, color = usepalette[dose_i]);
        if !ismissing(obs_mouseTV)
            plot!((@rsubset(obs_mouseTV, :DOSE == dose_i)).DAY, (@rsubset(obs_mouseTV, :DOSE == dose_i)).DV, seriestype = :scatter, markershape=:circle, markerstrokewidth=0, alpha = 0.5, color = usepalette[dose_i], label="data $(round(dose_i, digits=2)) mgkg");
        end
    end
    if !ismissing(yrange)
        plot!(ylims = yrange); 
    end
    if yaxis_log
        plot!(yaxis = :log);
    end
    return plot_TV
end


#--------------- visualization, T cells  ---------------#
function PlotTcells(sol; showBM = false, showLN = false, showSpleen=false, showTumor = true, legendpos = :bottomright, xInDays = false, xtickInDays = false, showNaive = true, showAct = true, showTotal = false, titlestr = missing)
    if xInDays
        plot__ = plot(xlabel = "Time (Day)", ylabel = "Cell count (#)", size = (400, 400), dpi = 300, legend = legendpos, background_color_legend = nothing);
        time_plots_ = sol.t / hr_per_day
        if xtickInDays
            plot!(xticks = (0:hr_per_day:maximum(sol.t)/hr_per_day));
        end
    else
        plot__ = plot(xlabel = "Time (hr)", ylabel = "Cell count (#)", size = (400, 400), dpi = 300, legend = legendpos, background_color_legend = nothing);
        time_plots_ = sol.t 
        if xtickInDays
            plot!(xticks = (0:hr_per_day:maximum(sol.t)));
        end
    end
    if showNaive
        plot!(time_plots_, sol[:restTpb], label = "naive T, PB", color = :cyan3, linestyle=:dash, lw = 2); 
    end
    if showAct
        plot!(time_plots_, sol[:actTpb], label = "active T, PB", color = :cyan3, linestyle=:solid); 
    end
    if showTotal
        plot!(time_plots_, sol[:restTpb] .+ sol[:actTpb], label = "total T, PB", color = :cyan3, linestyle=:dashdot); 
    end
    if showSpleen
        if showNaive
            plot!(time_plots_, sol[:restTtiss], label = "naive T, spleen", color = :deepskyblue, linestyle=:dash, lw = 2); 
        end
        if showAct
            plot!(time_plots_, sol[:actTtiss], label = "active T, spleen", color = :deepskyblue, linestyle=:solid); 
        end
        if showTotal
            plot!(time_plots_, sol[:restTtiss] .+ sol[:actTtiss], label = "total T, spleen", color = :deepskyblue, linestyle=:dashdot); 
        end
    end
    if showLN
        if showNaive
            plot!(time_plots_, sol[:restTtiss2], label = "naive T, LN", color = :violet, linestyle=:dash, lw = 2); 
        end
        if showAct
            plot!(time_plots_, sol[:actTtiss2], label = "active T, LN", color = :violet, linestyle=:solid); 
        end
        if showTotal
            plot!(time_plots_, sol[:restTtiss2] .+ sol[:actTtiss2], label = "total T, LN", color = :violet, linestyle=:dashdot); 
        end
    end
    if showBM
        if showNaive
            plot!(time_plots_, sol[:restTtiss3], label = "naive T, BM", color = :limegreen, linestyle=:dash, lw = 2); 
        end
        if showAct
            plot!(time_plots_, sol[:actTtiss3], label = "active T, BM", color = :limegreen, linestyle=:solid); 
        end
        if showTotal
            plot!(time_plots_, sol[:restTtiss3] .+ sol[:actTtiss3], label = "total T, BM", color = :limegreen, linestyle=:dashdot); 
        end
    end
    if showTumor
        if showNaive
            plot!(time_plots_, sol[:restTtumor], label = "naive T, tumor", color = :purple3, linestyle=:dash, lw = 2); 
        end
        if showAct
            plot!(time_plots_, sol[:actTtumor], label = "active T, tumor", color = :purple3, linestyle=:solid); 
        end
        if showTotal
            plot!(time_plots_, sol[:restTtumor] .+ sol[:actTtumor], label = "total T, tumor", color = :purple3, linestyle=:dashdot); 
        end
    end
    if !ismissing(titlestr)
        plot!(title = titlestr, titlefontsize = 12);
    end
    return plot__
end

function Plot_T_doubleAxis(sol; showPB = true, legendpos1 = :bottomright, legendpos2 = :topright, xInDays = false, xtickInDays = false, showNaive = false, showAct = true)
    if xInDays
        plot__ = plot(xlabel = "Time (Day)", ylabel = "Cell count, tumor (#)", size = (500, 400), dpi = 300, legend = legendpos1, background_color_legend = nothing);
        time_plots_ = sol.t / hr_per_day
        if xtickInDays
            plot!(xticks = (0:hr_per_day:maximum(sol.t)/hr_per_day));
        end
    else
        plot__ = plot(xlabel = "Time (hr)", ylabel = "Cell count, tumor (#)", size = (500, 400), dpi = 300, legend = legendpos1, background_color_legend = nothing);
        time_plots_ = sol.t 
        if xtickInDays
            plot!(xticks = (0:hr_per_day:maximum(sol.t)));
        end
    end
    axis2 = twinx();
    plot!(axis2, ylabel = "Cell count, not tumor (#)", legend = legendpos2, background_color_legend = nothing);
    if showNaive
        plot!(time_plots_, sol[:restTtumor], label = "naive T, tumor", color = :purple3, linestyle=:dash, lw = 2); 
        if showPB
            plot!(axis2, time_plots_, sol[:restTpb], label = "naive T, PB", color = :cyan3, linestyle=:dash, lw = 2); 
        end
    end
    if showAct
        plot!(time_plots_, sol[:actTtumor], label = "active T, tumor", color = :purple3, linestyle=:solid, lw = 1.5); 
        if showPB
            plot!(axis2, time_plots_, sol[:actTpb], label = "active T, PB", color = :cyan3, linestyle=:solid, lw = 1.5); 
        end
    end
    plot!(right_margin = 15Plots.mm, left_margin = 15Plots.mm); 
    return plot__
end

function Plot_Tumor_T_multidose(dose_mgkg_list, sims_outcome; showNaive = false, showAct = false, showTotal = false, xInDays = false, legendpos = :right, yrange = missing)
    usepalette = cgrad(:glasgowS, length(dose_mgkg_list))
    if xInDays
        plot_TuT = plot(background_color_legend = nothing, xlabel = "Time (Day)", ylabel = "Tumor T cell count (#)", dpi = 300, size = (400, 400), legend = legendpos);
        for (i, dose_i) in enumerate(dose_mgkg_list)
            if showNaive
                plot!(sims_outcome[dose_i].t/hr_per_day, sims_outcome[dose_i][pkpd. restTtumor], label="naive T, $(round(dose_i, digits=2)) mgkg", width=2, linestyle=:dashdot, color = usepalette[dose_i]);
            end
            if showAct
                plot!(sims_outcome[dose_i].t/hr_per_day, sims_outcome[dose_i][pkpd. actTtumor], label="active T, $(round(dose_i, digits=2)) mgkg", width=2, linestyle=:solid, color = usepalette[dose_i]);
            end
            if showTotal
                plot!(sims_outcome[dose_i].t/hr_per_day, sims_outcome[dose_i][pkpd. actTtumor] .+ sims_outcome[dose_i][pkpd. restTtumor], label="total T, $(round(dose_i, digits=2)) mgkg", width=2, linestyle=:dash, color = usepalette[dose_i]);
            end
        end
    else
        plot_TuT = plot(background_color_legend = nothing, xlabel = "Time (hr)", ylabel = "Tumor T cell count (#)", dpi = 300, size = (400, 400), legend = legendpos);
        for (i, dose_i) in enumerate(dose_mgkg_list)
            if showNaive
                plot!(sims_outcome[dose_i].t, sims_outcome[dose_i][pkpd. restTtumor], label="naive T, $(round(dose_i, digits=2)) mgkg", width=2, linestyle=:dashdot, color = usepalette[dose_i]);
            end
            if showAct
                plot!(sims_outcome[dose_i].t, sims_outcome[dose_i][pkpd. actTtumor], label="active T, $(round(dose_i, digits=2)) mgkg", width=2, linestyle=:solid, color = usepalette[dose_i]);
            end
            if showTotal
                plot!(sims_outcome[dose_i].t, sims_outcome[dose_i][pkpd. actTtumor] .+ sims_outcome[dose_i][pkpd. restTtumor], label="total T, $(round(dose_i, digits=2)) mgkg", width=2, linestyle=:dash, color = usepalette[dose_i]);
            end
        end
    end
    if !ismissing(yrange)
        plot!(ylims = yrange)
    end
    
    return plot_TuT
end
#--------------- visualization, E:T ratio  ---------------#
 
function Plot_E2T_multidose(dose_mgkg_list, sims_outcome; xrange = missing, showTotalE2T = false, showActE2T = false, showNaiveE2T = false, titlestr = "", legendpos = :topleft, yrange = [0, 1], ylabel_set = "Tumor E:T ratio")
    usepalette = cgrad(:glasgowS, length(dose_mgkg_list))
    plot_ETratio= plot(background_color_legend = nothing, xlabel = "Time (Day)", ylabel = ylabel_set, dpi = 300, size = (400, 400), legendtitle = "Dose", legendtitlefontsize = 8, ylims = yrange);
    if showTotalE2T
        for dose_i in dose_mgkg_list
            plot!(sims_outcome[dose_i].t/hr_per_day, (sims_outcome[dose_i][pkpd.actTtumor].+sims_outcome[dose_i][pkpd.restTtumor])./sims_outcome[dose_i][pkpd.Ntot], label= "$(round(dose_i, digits=2)) mgkg", width=2, linestyle=:solid, color = usepalette[dose_i], alpha = 0.8);
        end
    end
    if showActE2T
        for dose_i in dose_mgkg_list
            plot!(sims_outcome[dose_i].t/hr_per_day, sims_outcome[dose_i][pkpd.actTtumor]./sims_outcome[dose_i][pkpd.Ntot], label= "$(round(dose_i, digits=2)) mgkg", width=2, linestyle=:dash, color = usepalette[dose_i], alpha = 0.8);
        end
    end
    if showNaiveE2T
        for dose_i in dose_mgkg_list
            plot!(sims_outcome[dose_i].t/hr_per_day, sims_outcome[dose_i][pkpd.restTtumor]./sims_outcome[dose_i][pkpd.Ntot], label="$(round(dose_i, digits=2)) mgkg", width=2, linestyle=:dashdotdot, color = usepalette[dose_i], alpha = 0.8);
        end
    end
    if !ismissing(xrange)
        plot!(xlims = xrange);
    end
    plot!(title = titlestr, titlefontsize = 8, legend = legendpos);
    
    return plot_ETratio
end

function Plot_totalE2T_wObs(dose_mgkg_list, sims_outcome, obs; obs_doubleaxis = false, xrange = missing, showActE2T = false, titlestr = "", legendpos = :topleft, ylabelstr2 = "CD3 (cells/mm2)", yrange = [0, 1])
    usepalette = cgrad(:glasgowS, length(dose_mgkg_list));
    plot_ETratio = plot(background_color_legend = nothing, xlabel = "Time (Day)", ylabel = "Tumor E:T ratio", dpi = 300, size = (400, 400), legendtitle = "Dose", legendtitlefontsize = 8, legend = legendpos, ylims = yrange);
    axis2 = twinx(); 
    for dose_i in dose_mgkg_list
        plot!(sims_outcome[dose_i].t/hr_per_day, (sims_outcome[dose_i][pkpd.actTtumor].+sims_outcome[dose_i][pkpd.restTtumor])./sims_outcome[dose_i][pkpd.Ntot], label= "Total E:T, $(round(dose_i, digits=2)) mgkg", width=2, linestyle=:solid, color = usepalette[dose_i], alpha = 0.8);
        if obs_doubleaxis
            obs_tmp = @rsubset(obs, :DOSE == dose_i);
            scatter!(axis2, obs_tmp.DAY, obs_tmp.DV, label= "$(round(dose_i, digits=2)) mgkg", color = usepalette[dose_i], ma = 0.6, ms = 6, markerstrokewidth = 0, ylabel = ylabelstr2, background_color_legend = nothing, legend = :topleft);
        else
            obs_tmp = @rsubset(obs, :DOSE == dose_i);
            scatter!(obs_tmp.DAY, obs_tmp.DV, label= "$(round(dose_i, digits=2)) mgkg", mc = usepalette[dose_i], ma = 0.6, ms = 6, markerstrokewidth = 0, background_color_legend = nothing, legend = :topleft);
        end
    end
    if showActE2T
        for dose_i in dose_mgkg_list
            plot!(sims_outcome[dose_i].t/hr_per_day, sims_outcome[dose_i][pkpd.actTtumor]./sims_outcome[dose_i][pkpd.Ntot], label= "Active E:T $(round(dose_i, digits=2)) mgkg", width=2, linestyle=:dash, color = usepalette[dose_i], alpha = 0.8);
        end
    end
    if !ismissing(xrange)
        plot!(xlims = xrange);
    end
    plot!(title = titlestr, titlefontsize = 8);
    plot!(right_margin = 15Plots.mm); 
    return plot_ETratio
end

function Plot_E2T_calculatedE2T(dose_mgkg_list, sims_outcome, obs; xrange = missing, showActE2T = false, titlestr = "", legendpos = :topleft)
    usepalette = cgrad(:glasgowS, length(dose_mgkg_list));
    plot_ETratio = plot(background_color_legend = nothing, xlabel = "Time (Day)", ylabel = "Tumor E:T ratio", dpi = 300, size = (400, 400), legendtitle = "Dose", legendtitlefontsize = 8, legend = legendpos);
    for dose_i in dose_mgkg_list
        plot!(sims_outcome[dose_i].t/hr_per_day, (sims_outcome[dose_i][pkpd.actTtumor].+sims_outcome[dose_i][pkpd.restTtumor])./(sims_outcome[dose_i][pkpd.Ntot]), label= "Total E:T, $(round(dose_i, digits=2)) mgkg", width=2, linestyle=:solid, color = usepalette[dose_i], alpha = 0.8);
        obs_tmp = @rsubset(obs, :DOSE == dose_i);
        scatter!(obs_tmp.DAY, obs_tmp.e2t, label= "$(round(dose_i, digits=2)) mgkg", mc = usepalette[dose_i], ma = 0.6, ms = 6, markerstrokewidth = 0, background_color_legend = nothing, legend = :topleft);
    end
    if showActE2T
        for dose_i in dose_mgkg_list
            plot!(sims_outcome[dose_i].t/hr_per_day, sims_outcome[dose_i][pkpd.actTtumor]./sims_outcome[dose_i][pkpd.Ntot], label= "Active E:T $(round(dose_i, digits=2)) mgkg", width=2, linestyle=:dash, color = usepalette[dose_i], alpha = 0.8);
        end
    end
    if !ismissing(xrange)
        plot!(xlims = xrange);
    end
    plot!(title = titlestr, titlefontsize = 8);
    return plot_ETratio
end

#--------------- visualization, trimer  ---------------#
function Plot_trimer_multidose(dose_mgkg_list, sims_outcome; xrange = missing, legendpos = :right)
    usepalette = cgrad(:glasgowS, length(dose_mgkg_list))
    plot_abs_Trim = plot(background_color_legend = nothing, xlabel = "Time (Day)", ylabel = "Trimer in tumor (#)", legend = legendpos, dpi = 300);
    for dose_i in dose_mgkg_list
        plot!(sims_outcome[dose_i].t/hr_per_day, sims_outcome[dose_i][pkpd.trimer_num], label="Trimer#  $(round(dose_i, digits=2)) mgkg", width=2, linestyle=:dash, color = usepalette[dose_i]);
    end
    if !ismissing(xrange)
        plot!(xlims = xrange);
    end
    return plot_abs_Trim
end

function Plot_trimer_conc_multidose(dose_mgkg_list, sims_outcome; xrange = missing, legendpos = :right, logy = false, yrange = missing)
    usepalette = cgrad(:glasgowS, length(dose_mgkg_list))
    plot_trimer_conc = plot(background_color_legend = nothing, xlabel = "Time (Day)", ylabel = "Trimer in TME (nM)", legend = legendpos, dpi = 300);
    for dose_i in dose_mgkg_list
        plot!(sims_outcome[dose_i].t/hr_per_day, sims_outcome[dose_i][pkpd.CD3_TCE_TAA_t].+ sims_outcome[dose_i][pkpd.CD3_TCE_TAA_TAA_t], label="$(round(dose_i, digits=2)) mgkg", width=2, linestyle=:dash, color = usepalette[dose_i]);
    end
    if !ismissing(xrange)
        plot!(xlims = xrange);
    end
    if !ismissing(yrange)
        plot!(ylims = yrange);
    end
    if logy
        plot!(yaxis = :log);
    end
    return plot_trimer_conc
end

function Plot_trimer_per_cell(dose_mgkg_list, sims_outcome; xrange = missing, legendpos = :right, showTri_Tu = false, showTri_Tu_actT = false, showTri_T = false, 
                                legendtitlestr = "", yaxis_log = false, yrange = missing, ylabelstr = "Trimer per cell (#)", tt50_kill = missing)
    usepalette = cgrad(:glasgowS, length(dose_mgkg_list))
    plot_Trim = plot(background_color_legend = nothing, xlabel = "Time (Day)", ylabel = ylabelstr, yguidefontsize=10, legend = legendpos, dpi = 300, size = (400, 400));
    if showTri_Tu
        for dose_i in dose_mgkg_list
            plot!(sims_outcome[dose_i].t/hr_per_day, sims_outcome[dose_i][pkpd.Tri_per_Tu], label="$(round(dose_i, digits=2)) mg/kg", width=2, linestyle=:dashdot, color = usepalette[dose_i]);
        end
    end
    if showTri_Tu_actT
        for dose_i in dose_mgkg_list
            plot!(sims_outcome[dose_i].t/hr_per_day, sims_outcome[dose_i][pkpd.Trimer_per_Tu_cell_on_act_T], label="$(round(dose_i, digits=2)) mg/kg", width=2, linestyle=:dash, color = usepalette[dose_i]);
        end
    end
    if showTri_T
        for dose_i in dose_mgkg_list
            plot!(sims_outcome[dose_i].t, sims_outcome[dose_i][pkpd.Tri_per_T], label="$(round(dose_i, digits=2)) mg/kg", width=2, linestyle=:solid, color = usepalette[dose_i]);
        end
    end
    if !ismissing(yrange)
        plot!(ylims = yrange);
    end
    if yaxis_log
        plot!(yaxis = :log);
    end
    if !ismissing(xrange)
        plot!(xlims = xrange);
    end
    if !ismissing(tt50_kill)
        hline!([tt50_kill], linestyle = :dot, lw = 2, alpha = 0.3, color = :black, label = "TT50_kill");
    end
    plot!(legendtitle = legendtitlestr, legendtitlefontsize = 6);
    return plot_Trim
end

#--------------- visualization, cytokine  ---------------#
function Plot_cytokines_(dose_mgkg_list, sims_outcome, obs = missing; showIL10 = false, showTNFa = false, showIFNg = false, MW_IL10 = MW_IL10, MW_TNFa = MW_TNFa, MW_IFNg = MW_IFNg, 
                         ylabelstr = "Tumor cytokine (nM)", ylabelstr2 = "PB cytokine (nM)", legendpos = :topleft, legendpos2 = :right)

    usepalette = cgrad(:glasgowS, length(dose_mgkg_list))
    plot_cytokine = plot(background_color_legend = nothing, xlabel = "Time (Day)", ylabel = ylabelstr, yguidefontsize=10, legend = legendpos, dpi = 300, size = (400, 300));
    axis2 = twinx(); 
    if showIL10
        for dose_i in dose_mgkg_list
            plot!(sims_outcome[dose_i].t/hr_per_day, sims_outcome[dose_i][pkpd.cytokine1_t], label="IL10 $(round(dose_i, digits=2)) mg/kg", width=2, linestyle=:solid, color = usepalette[dose_i]);
        end
    end
    if showTNFa
        for dose_i in dose_mgkg_list
            plot!(sims_outcome[dose_i].t/hr_per_day, sims_outcome[dose_i][pkpd.cytokine2_t], label="TNFa $(round(dose_i, digits=2)) mg/kg", width=2, linestyle=:dash, color = usepalette[dose_i]);
        end
    end
    if showIFNg
        for dose_i in dose_mgkg_list
            plot!(sims_outcome[dose_i].t/hr_per_day, sims_outcome[dose_i][pkpd.cytokine3_t], label="IFNg $(round(dose_i, digits=2)) mg/kg", width=2, linestyle=:dashdotdot, color = usepalette[dose_i]);
        end
    end
    if !ismissing(obs)
        if showIL10
            for dose_i in dose_mgkg_list
                obs_tmp = @rsubset(obs, :DOSE == dose_i);
                scatter!(axis2, obs_tmp.DAY, obs_tmp.DV / MW_IL10, label="$(round(dose_i, digits=2)) mg/kg", mc = usepalette[dose_i], ma = 0.4, ms = 6, markerstrokewidth = 0, ylabel = ylabelstr2, yguidefontsize=10, background_color_legend = nothing, legend = legendpos2);
            end
        end
        if showTNFa
            for dose_i in dose_mgkg_list
                obs_tmp = @rsubset(obs, :DOSE == dose_i);
                scatter!(axis2, obs_tmp.DAY, obs_tmp.DV / MW_TNFa, label="$(round(dose_i, digits=2)) mg/kg", mc = usepalette[dose_i], ma = 0.4, ms = 6, markerstrokewidth = 0, ylabel = ylabelstr2, yguidefontsize=10, background_color_legend = nothing, legend = legendpos2);
            end
        end
        if showIFNg
            for dose_i in dose_mgkg_list
                obs_tmp = @rsubset(obs, :DOSE == dose_i);
                scatter!(axis2, obs_tmp.DAY, obs_tmp.DV / MW_IFNg, label="$(round(dose_i, digits=2)) mg/kg", mc = usepalette[dose_i], ma = 0.4, ms = 6, markerstrokewidth = 0, ylabel = ylabelstr2, yguidefontsize=10, background_color_legend = nothing, legend = legendpos2);
            end
        end
    end
    plot!(right_margin = 10Plots.mm, left_margin = 10Plots.mm); 
    return plot_cytokine
end


function Plot_cytokines_pb(dose_mgkg_list, sims_outcome, obs = missing; showIL10 = false, showTNFa = false, showIFNg = false, MW_IL10 = MW_IL10, MW_TNFa = MW_TNFa, MW_IFNg = MW_IFNg, 
                         ylabelstr = "PB cytokine (nM)", ylabelstr2 = "PB cytokine (nM)", legendpos = :topleft, legendpos2 = :right)

    usepalette = cgrad(:glasgowS, length(dose_mgkg_list))
    plot_cytokine = plot(background_color_legend = nothing, xlabel = "Time (Day)", ylabel = ylabelstr, yguidefontsize=10, legend = legendpos, dpi = 300, size = (400, 300));
    axis2 = twinx(); 
    if showIL10
        for dose_i in dose_mgkg_list
            plot!(sims_outcome[dose_i].t/hr_per_day, sims_outcome[dose_i][pkpd.cytokine1_pb], label="IL10 $(round(dose_i, digits=2)) mg/kg", width=2, linestyle=:solid, color = usepalette[dose_i]);
        end
    end
    if showTNFa
        for dose_i in dose_mgkg_list
            plot!(sims_outcome[dose_i].t/hr_per_day, sims_outcome[dose_i][pkpd.cytokine2_pb], label="TNFa $(round(dose_i, digits=2)) mg/kg", width=2, linestyle=:dash, color = usepalette[dose_i]);
        end
    end
    if showIFNg
        for dose_i in dose_mgkg_list
            plot!(sims_outcome[dose_i].t/hr_per_day, sims_outcome[dose_i][pkpd.cytokine3_pb], label="IFNg $(round(dose_i, digits=2)) mg/kg", width=2, linestyle=:dashdotdot, color = usepalette[dose_i]);
        end
    end
    if !ismissing(obs)
        if showIL10
            for dose_i in dose_mgkg_list
                obs_tmp = @rsubset(obs, :DOSE == dose_i);
                scatter!(axis2, obs_tmp.DAY, obs_tmp.DV / MW_IL10, label="$(round(dose_i, digits=2)) mg/kg", mc = usepalette[dose_i], ma = 0.4, ms = 6, markerstrokewidth = 0, ylabel = ylabelstr2, yguidefontsize=10, background_color_legend = nothing, legend = legendpos2);
            end
        end
        if showTNFa
            for dose_i in dose_mgkg_list
                obs_tmp = @rsubset(obs, :DOSE == dose_i);
                scatter!(axis2, obs_tmp.DAY, obs_tmp.DV / MW_TNFa, label="$(round(dose_i, digits=2)) mg/kg", mc = usepalette[dose_i], ma = 0.4, ms = 6, markerstrokewidth = 0, ylabel = ylabelstr2, yguidefontsize=10, background_color_legend = nothing, legend = legendpos2);
            end
        end
        if showIFNg
            for dose_i in dose_mgkg_list
                obs_tmp = @rsubset(obs, :DOSE == dose_i);
                scatter!(axis2, obs_tmp.DAY, obs_tmp.DV / MW_IFNg, label="$(round(dose_i, digits=2)) mg/kg", mc = usepalette[dose_i], ma = 0.4, ms = 6, markerstrokewidth = 0, ylabel = ylabelstr2, yguidefontsize=10, background_color_legend = nothing, legend = legendpos2);
            end
        end
    end
    plot!(right_margin = 10Plots.mm, left_margin = 10Plots.mm); 
    return plot_cytokine
end


function Plot_cytokines_pb_2(dose_mgkg_list, sims_outcome, obs = missing; showIL10 = false, showTNFa = false, showIFNg = false, 
                             MW_IL10 = MW_IL10, MW_TNFa = MW_TNFa, MW_IFNg = MW_IFNg, ylabelstr = "PB cytokine (pg/mL)", 
                             legendpos = :topleft, ytype = :log, yrange = [1E-3, 1E5])

    usepalette = cgrad(:glasgowS, length(dose_mgkg_list))
    plot_cytokine = plot(background_color_legend = nothing, xlabel = "Time (Day)", ylabel = ylabelstr, yguidefontsize=10, legend = legendpos, dpi = 300, size = (400, 300));
    if showIL10
        for dose_i in dose_mgkg_list
            plot!(sims_outcome[dose_i].t/hr_per_day, sims_outcome[dose_i][pkpd.cytokine1_pb] * MW_IL10, label="IL10 $(round(dose_i, digits=2)) mg/kg", width=2, linestyle=:solid, color = usepalette[dose_i]);
        end
    end
    if showTNFa
        for dose_i in dose_mgkg_list
            plot!(sims_outcome[dose_i].t/hr_per_day, sims_outcome[dose_i][pkpd.cytokine2_pb] * MW_TNFa, label="TNFa $(round(dose_i, digits=2)) mg/kg", width=2, linestyle=:dash, color = usepalette[dose_i]);
        end
    end
    if showIFNg
        for dose_i in dose_mgkg_list
            plot!(sims_outcome[dose_i].t/hr_per_day, sims_outcome[dose_i][pkpd.cytokine3_pb] * MW_IFNg, label="IFNg $(round(dose_i, digits=2)) mg/kg", width=2, linestyle=:dashdotdot, color = usepalette[dose_i]);
        end
    end
    if !ismissing(obs)
        if showIL10
            for dose_i in dose_mgkg_list
                obs_tmp = @rsubset(obs, :DOSE == dose_i);
                scatter!(obs_tmp.DAY, obs_tmp.DV, label="$(round(dose_i, digits=2)) mg/kg", mc = usepalette[dose_i], ma = 0.4, ms = 6, markerstrokewidth = 0);
            end
        end
        if showTNFa
            for dose_i in dose_mgkg_list
                obs_tmp = @rsubset(obs, :DOSE == dose_i);
                scatter!(obs_tmp.DAY, obs_tmp.DV, label="$(round(dose_i, digits=2)) mg/kg", mc = usepalette[dose_i], ma = 0.4, ms = 6, markerstrokewidth = 0);
            end
        end
        if showIFNg
            for dose_i in dose_mgkg_list
                obs_tmp = @rsubset(obs, :DOSE == dose_i);
                scatter!(obs_tmp.DAY, obs_tmp.DV, label="$(round(dose_i, digits=2)) mg/kg", mc = usepalette[dose_i], ma = 0.4, ms = 6, markerstrokewidth = 0);
            end
        end
    end
    plot!(yaxis = ytype, ylims = yrange); 
    plot!(right_margin = 10Plots.mm, left_margin = 10Plots.mm); 
    return plot_cytokine
end