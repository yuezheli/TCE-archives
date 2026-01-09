function tceTMEModel(pobject ;name)
    # Define Independent Variable
    @independent_variables t, [unit = u"s"]
    # Load some common constants
    ModelTrace.@common_constants
    pars = @parameters_from_spec pobject begin
        Tumor_cells_init
        total_T_cells_init
        Vwell
        kon_taa_tce
        kon_cd3_tce
        KD_TAA
        KD_CD3
        lambda_taa_avidity
        lambda_trimer_T
        lambda_trimer_tumor
        lambda_int
        Thalf_TAA
        Thalf_CD3
        CD3_per_cell
        TAA_per_cell
        Tumor_cell_doubling_time
        k_growth
        emax_apop
        tt50_kill
        n_kill
        FLG_trimer_per_tumor_kill
        k_t_prolif
        frac_act_T_death
        tt50_t_activate
        emax_t_activate
        n_act
        k_syn_cytokine
        tt50_cytokine
        n_syn_cytokine
        k_deg_cytokine    
    end 

    vars = @variables_from_spec pobject begin
        TCE(t)
        CD3(t)
        TAA(t)
        CD3_TCE(t)
        TAA_TCE(t)
        TAA_TAA_TCE(t)
        CD3_TCE_TAA(t)
        CD3_TCE_TAA_TAA(t)
        Tumor_cell(t)
        T_cell(t)
        T_act(t)
        cytokine(t)
        deg_TAA_TCE(t)
        dead_T(t)
        dead_Tumor(t)
    end
    D = Differential(t)

    koff_CD3 = kon_cd3_tce * KD_CD3                # [s-1]
    kdeg_CD3 = log(2)/(Thalf_CD3 * s_per_hr)       # CD3 degredation rate; [s-1]
    ksyn_CD3 = kdeg_CD3 * (CD3_per_cell * T_cell / N_Av * nmol_per_mol)/Vwell # CD3 synthesis rate per well; [nM.s-1]
    # TAA related parameters 
    koff_TAA = kon_taa_tce * KD_TAA                              # [s-1]
    kdeg_TAA = log(2)/(Thalf_TAA * s_per_hr)              # TAA degredation rate; [s-1]
    ksyn_TAA = kdeg_TAA * (TAA_per_cell * Tumor_cell / N_Av * nmol_per_mol)/Vwell # TAA synthesis rate per well; [nM.s-1]

    # binding dynamics
    flux_TAA_TCE_binding = kon_taa_tce * TAA * TCE - koff_TAA * TAA_TCE  # TAA + TCE <-> TAA:TCE
    flux_TAA_TAA_TCE_binding = lambda_taa_avidity * kon_taa_tce * TAA_TCE * TAA - koff_TAA * TAA_TAA_TCE  # TAA + TAA:TCE <-> TAA:TAA:TCE
    flux_CD3_TCE_binding = kon_cd3_tce * CD3 * TCE - koff_CD3 * CD3_TCE  # CD3 + TCE <-> CD3:TCE
    flux_TAA_trimer_binding_1 = lambda_trimer_tumor * kon_taa_tce * TAA * CD3_TCE - koff_TAA * CD3_TCE_TAA  # TAA + CD3:TCE <-> CD3:TCE:TAA 
    flux_TAA_trimer_binding_2 = lambda_trimer_T * kon_cd3_tce * CD3 * TAA_TCE - koff_CD3 * CD3_TCE_TAA  # CD3 + TAA:TCE <-> CD3:TCE:TAA 
    flux_TAA_tetramer_binding_1 = lambda_taa_avidity * kon_taa_tce * CD3_TCE_TAA * TAA - koff_TAA * CD3_TCE_TAA_TAA # CD3:TCE:TAA + TAA <-> CD3:TCE:TAA:TAA ## AKW Comment: it isn't truly a tetramer since the trimer is between the 2 cells and the molecule. The equations are fine but I'd consider a different name. YL: tetramer in the sense of the protein, not cells
    flux_TAA_tetramer_binding_2 = lambda_trimer_T * kon_cd3_tce * CD3 * TAA_TAA_TCE - koff_CD3 * CD3_TCE_TAA_TAA # TAA:TAA:TCE + CD3 <-> CD3:TCE:TAA:TAA
    # trimer & cell dynamics 
    trimer = (CD3_TCE_TAA + CD3_TCE_TAA_TAA)/nmol_per_mol*N_Av*Vwell 
    trimer_per_Tumor_cell = max(trimer / max(Tumor_cell, 1), 0)
    trimer_per_Tumor_cell_on_act_T = trimer_per_Tumor_cell * max(T_act, 0) / max(T_cell, 1)
    trimer_per_T = max(trimer / max(T_cell, 1), 0)
    # tumor killing (driven by either trimer on tumor cells or trimer on T cell)
    flux_tce_apop_tumor = emax_apop * (trimer_per_Tumor_cell_on_act_T^n_kill) / (tt50_kill^n_kill + trimer_per_Tumor_cell_on_act_T^n_kill)
    flux_tce_apop_T = emax_apop * (trimer_per_T^n_kill) / (tt50_kill^n_kill + trimer_per_T^n_kill) * max(T_act, 0) / max(T_cell, 1) # QUESTIONABLE IMPLEMENTATION
    flux_tce_apop = flux_tce_apop_tumor * FLG_trimer_per_tumor_kill + flux_tce_apop_T * (1 - FLG_trimer_per_tumor_kill)
    # T cell activation 
    flux_t_act = emax_t_activate * (trimer_per_T^n_act)/(tt50_t_activate^n_act + trimer_per_T^n_act)

    eqs = [
        # TCE concentration in medium 
        D(TCE) ~ ( - flux_TAA_TCE_binding - flux_CD3_TCE_binding ), 
        # CD3 concentration in medium 
        D(CD3) ~ ( ksyn_CD3 - kdeg_CD3*CD3  # synthesis and degredation of CD3 on T cells 
                    - flux_CD3_TCE_binding - flux_TAA_trimer_binding_2 - flux_TAA_tetramer_binding_2  # binding activity 
                    + flux_tce_apop*(CD3_TCE_TAA + CD3_TCE_TAA_TAA) * (1-frac_act_T_death)  # CD3 release from trimer/tetramer after tumorlysis
                    + k_t_prolif*T_act*(CD3_per_cell*nmol_per_mol/N_Av/Vwell) ),  # CD3 from T cell proliferation
        # TAA concentration in medium 
        D(TAA) ~ ( ksyn_TAA - kdeg_TAA*TAA  # synthesis and degredation of TAA on tumor cells 
                    - flux_TAA_TCE_binding - flux_TAA_TAA_TCE_binding - flux_TAA_trimer_binding_1 - flux_TAA_tetramer_binding_1  # binding activity 
                    - flux_tce_apop*TAA # TAA loss from tumor cell death 
                    + k_growth*Tumor_cell*(TAA_per_cell*nmol_per_mol)/N_Av/Vwell ), # TAA from tumor cell proliferation 
        # CD3:TCE concentration in medium 
        D(CD3_TCE) ~ (flux_CD3_TCE_binding - flux_TAA_trimer_binding_1  # binding activity 
                        - kdeg_CD3*CD3_TCE), # CD3:TCE loss due to CD3 turnover 
        # TAA:TCE in medium 
        D(TAA_TCE) ~ (flux_TAA_TCE_binding - flux_TAA_TAA_TCE_binding - flux_TAA_trimer_binding_2  # binding activity 
                        - kdeg_TAA*TAA_TCE  # TAA:TCE loss due to TAA turnover
                        - flux_tce_apop*TAA_TCE ), # TAA:TCE loss tumorlysis
        # TAA:TAA:TCE in medium 
        D(TAA_TAA_TCE) ~ (flux_TAA_TAA_TCE_binding - flux_TAA_tetramer_binding_2  # binding activity 
                            - lambda_int*kdeg_TAA*TAA_TAA_TCE  # TAA:TAA:TCE loss due to TAA turnover 
                            - flux_tce_apop*TAA_TAA_TCE ),  # TAA:TAA:TCE loss tumorlysis
        # CD3:TCE:TAA in medium 
        D(CD3_TCE_TAA) ~ (flux_TAA_trimer_binding_1 + flux_TAA_trimer_binding_2 - flux_TAA_tetramer_binding_1  # binding activity 
                            - flux_tce_apop*CD3_TCE_TAA ),  # CD3:TCE:TAA loss tumorlysis
        # CD3:TCE:TAA:TAA in medium 
        D(CD3_TCE_TAA_TAA) ~ ( flux_TAA_tetramer_binding_1 + flux_TAA_tetramer_binding_2  # binding activity 
                                - flux_tce_apop*CD3_TCE_TAA_TAA ),    # CD3:TCE:TAA:TAA loss tumorlysis
        # Tumor cells in medium 
        D(Tumor_cell) ~ ( k_growth*Tumor_cell - flux_tce_apop*Tumor_cell  ), 
        # T cell in medium 
        D(T_cell) ~ (k_t_prolif*T_act - frac_act_T_death*flux_tce_apop*T_cell ), 
        D(T_act) ~ (k_t_prolif*T_act + flux_t_act*(T_cell - T_act)), 
        # cytokine 
        D(cytokine) ~ (k_syn_cytokine*(trimer_per_T^n_syn_cytokine)/(tt50_cytokine^n_syn_cytokine + trimer_per_T^n_syn_cytokine) * max(T_act, 0) - k_deg_cytokine*cytokine), 
        # dummy variable 
        D(deg_TAA_TCE) ~ ( kdeg_TAA*TAA_TCE + lambda_int*kdeg_TAA*TAA_TAA_TCE ), 
        D(dead_T) ~ ( frac_act_T_death*flux_tce_apop*T_cell ),
        D(dead_Tumor) ~ (flux_tce_apop*Tumor_cell),
    ]    

    R0_CD3 = (CD3_TCE + CD3_TCE_TAA + CD3_TCE_TAA_TAA) / (CD3 + CD3_TCE + CD3_TCE_TAA + CD3_TCE_TAA_TAA)
    R0_TAA = (TAA_TCE + TAA_TAA_TCE + CD3_TCE_TAA_TAA) / (TAA + TAA_TCE + TAA_TAA_TCE + CD3_TCE_TAA_TAA)
    trimer_per_T = trimer / max(T_cell, 1)
    trimer_per_Tumor = trimer / max(Tumor_cell, 1)
    trimer_act_T_per_Tumor = trimer / max(Tumor_cell, 1) * max(T_act, 0) / max(T_cell, 1)
    trimer_1_per_tumor_cell = CD3_TCE_TAA / nmol_per_mol * N_Av * Vwell / max(Tumor_cell, 1)
    trimer_2_per_tumor_cell = CD3_TCE_TAA_TAA / nmol_per_mol * N_Av * Vwell / max(Tumor_cell, 1)

     return MRGODESystem(eqs, t; name=name, observed = [@observed(R0_CD3), @observed(R0_TAA), @observed(trimer_per_T), @observed(trimer_per_Tumor), @observed(trimer_act_T_per_Tumor), @observed(trimer_1_per_tumor_cell), @observed(trimer_2_per_tumor_cell)])
end
