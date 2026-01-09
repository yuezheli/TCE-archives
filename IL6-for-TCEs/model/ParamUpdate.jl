# date: 8/30/24
# author: Yuezhe Li 
# purpose of this code: to update parameters for different TDBs (mosunetuzumab, glofitamab, epcoritamab)

include("param.jl");

function TDB_param_update(TDB = "mosunetuzumab"; p_base = p_homo_3)
    if TDB == "mosunetuzumab"
        p_new = deepcopy(p_base);
    elseif TDB == "glofitamab"
        p_new = deepcopy(p_base);
        # glofit PK update; note this PK only work for patients treated with 1000mg obinutuzumab 7 days prior to the first dose of glofitamab
        p_new.CL_TDB = 240.        # [mL/d]
        p_new.Q_TDB  = 432.        # [mL/d]
        p_new.V1_TDB = 3126.       # [mL]
        p_new.V2_TDB = 12.1E3      # [mL]
        # glofit PD; optimized based on literature data
        p_new.KmBT_act = 0.98;   # p_mosun.KmBT_act = 0.72
        p_new.KdrugactT = 0.41;  # p_mosun.KdrugactT = 130.08
        p_new.KdrugB = 1E-4;     # p_mosun.KdrugB = 1.3
        p_new.Vm_tdb = 0.        # remove binding partner in glofit
    elseif TDB == "epcoritamab"
        p_new = deepcopy(p_base);
        # epco PK update; allometric scaling of cyno epco PK parameters
        BW_cyno = 3.4; # [kg]; # Hosseini et al., 2020
        # scaling factor between cyno and human 
        f_scale_CL = (p_new.BW / BW_cyno)^(0.7); 
        f_scale_V = (p_new.BW / BW_cyno)^(0.9);
        p_cyno = ComponentArray(CL = 81.6, Q = 82.187, V1 = 125.052, V2 = 306., ka = 0.26); # cyno param from previuos fit
        p_new.CL_TDB =  p_cyno.CL*f_scale_CL       # [mL/d]
        p_new.Q_TDB  = p_cyno.Q*f_scale_CL       # [mL/d]
        p_new.V1_TDB = p_cyno.V1*f_scale_V      # [mL]
        p_new.V2_TDB = p_cyno.V2*f_scale_V      # [mL]
        p_new.fbio_TDB = 0.95             # unitless
        # parameter fit with RI-1 (Figure 1F Engelberts et al 2020)
        p_new.KmBT_act =  1.0;  # p_mosun.KmBT_act = 0.72
        p_new.KdrugactT = 0.1;  # p_mosun.KdrugactT = 130.08
        p_new.KdrugB = 0.19;     # p_mosun.KdrugB = 1.3
        p_new.Vm_tdb = 0.        # remove binding partner in epco
    else
        println("Need TDB choice; otherwise return default values for mosunetuzumab")
        p_new = deepcopy(p_base);
    end
    return p_new; 
end

