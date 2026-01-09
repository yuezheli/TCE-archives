# date: 3/12/2025
# author: Yuezhe Li
# purpose of this code: see readme.txt

function invitro(;name)
    pars = @parameters begin 
        # initialization of the model 
        Tumor_cells_init = 2E4          
        total_T_cells_init = 2E5
        # Binding parameters
        kon_taa_tce = 1.2E-4            # protein binding rate constant; [nM-1.s-1]
        kon_cd3_tce = 3.1E-4            # protein binding rate constant; [nM-1.s-1]
        KD_TAA = 20                     # TCE:TAA Binding Affinity; [nM]
        KD_CD3 = 230                    # TCE:CD3 binding affinity; [nM] human
        lambda_taa_avidity = 1          # scaling factor for cis-avidity of TAA:TCE binding; should be between 1E2 - 1E5; Harms et al., 2014; https://pubmed.ncbi.nlm.nih.gov/23872324/
        lambda_trimer_T = 1             # placeholder trans avidity for trimer formation on T cells 
        lambda_trimer_tumor = 1         # placeholder trans avidity for trimer formation on tumor cells 
        lambda_int = 0                  # placeholder parameter for internalization rates ratio for TAA:TAA:TCE than TAA:TCE
        # Receptor kinetics
        Thalf_TAA = 0.1                 # Internalization half-life of TAA; [hr]; tuned based on internalization data
        Thalf_CD3 = 3.5                 # Internalization half-life of CD3; [hr]; https://journals.aai.org/jimmunol/article/173/1/384/72940/Constitutive-and-Ligand-Induced-TCR-Degradation-1
        # Cell line parameters
        CD3_per_cell = 57000            # Expression of CD3 per T cell (molecule/cell); https://pubmed.ncbi.nlm.nih.gov/8817090/; could range between 5000 - 70000 # https://pubmed.ncbi.nlm.nih.gov/19834968/
        TAA_per_cell = 9.4E4            # Expression of TAA per tumor cell (molecule/cell); # CD20 expression level obtained from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC500695/pdf/jclinpath00266-0020.pdf, though different from https://pubmed.ncbi.nlm.nih.gov/8817090/
        # Tumor cell kinetics
        Tumor_cell_doubling_time = 24   # Tumor doubling time [hr] Toledo cell line doubling time obtained from https://www.cellosaurus.org/CVCL_3611
        k_growth = log(2)/(Tumor_cell_doubling_time * s_per_hr)  # Tumor growth rate [s-1]
        emax_apop = 6.0e-6              # Rate of tumor cell apoptosis due to trimers (CD3:TCE:TAA) on tumor cell; [s-1]
        tt50_kill = 1                   # Number of trimers per cell resulting in half-max cell killing; [#/ cell]
        n_kill = 1                      # Hill coefficient on trimer-induced tumor cell killing
        FLG_trimer_per_tumor_kill = 1   # FLG for triggering tumor cell death based on trimers on tumor cell (if turned off, then trigger death based on trimer per T cell)
        # T cell kinetics 
        k_t_prolif = 0                  # placeholder for active T cell proliferation 
        frac_act_T_death = 0            # placeholder for the fraction of active T cell death after trimer formation 
        tt50_t_activate = 1             # EC50 of trimer to induce T cell activation 
        emax_t_activate = 0             # maximum rate for trimer-induced T cell activation 
        n_act = 1                       # Hill coefficient on trimer-induced T cell activation 
        # cytokine dynamics
        k_syn_cytokine = 0              # maximum rate for trimer-induced (per T cell) cytokine production per active T cell
        tt50_cytokine = 1               # EC50 of trimer-induced (per T cell) cytokine production 
        n_syn_cytokine = 1              # Hill coefficient trimer-induced (per T cell) cytokine production 
        k_deg_cytokine = 0   
    end
     # define variables
     @independent_variables t; 
     Dt = Differential(t); 
     vars = @variables begin 
        TCE(t) = 0.
        CD3(t) = CD3_per_cell*total_T_cells_init/NAv/Vwell * nmol_per_mol          # init conc of CD3; [nM]; total CD3 (on both naive T and active T)
        TAA(t) = TAA_per_cell*Tumor_cells_init/NAv/Vwell * nmol_per_mol      # init conc of TAA; [nM] 
        CD3_TCE(t) = 0.0
        TAA_TCE(t) = 0.0
        TAA_TAA_TCE(t) = 0.0
        CD3_TCE_TAA(t) = 0.0
        CD3_TCE_TAA_TAA(t) = 0.0
        Tumor_cell(t) = Tumor_cells_init
        T_cell(t) = total_T_cells_init            # total T cells in the well (naive T + act T)
        T_act(t) = 0.
        cytokine(t) = 0
        # dummy variable 
        deg_TAA_TCE(t) = 0.
        dead_T(t) = 0. 
        dead_Tumor(t) = 0.
    end
    # CD3 related parameters 
    koff_CD3 = kon_cd3_tce * KD_CD3                # [s-1]
    kdeg_CD3 = log(2)/(Thalf_CD3 * s_per_hr)       # CD3 degredation rate; [s-1]
    ksyn_CD3 = kdeg_CD3 * (CD3_per_cell * T_cell / NAv * nmol_per_mol)/Vwell # CD3 synthesis rate per well; [nM.s-1]
    # TAA related parameters 
    koff_TAA = kon_taa_tce * KD_TAA                              # [s-1]
    kdeg_TAA = log(2)/(Thalf_TAA * s_per_hr)              # TAA degredation rate; [s-1]
    ksyn_TAA = kdeg_TAA * (TAA_per_cell * Tumor_cell / NAv * nmol_per_mol)/Vwell # TAA synthesis rate per well; [nM.s-1]
    # binding dynamics
    flux_TAA_TCE_binding = kon_taa_tce * TAA * TCE - koff_TAA * TAA_TCE  # TAA + TCE <-> TAA:TCE
    flux_TAA_TAA_TCE_binding = lambda_taa_avidity * kon_taa_tce * TAA_TCE * TAA - koff_TAA * TAA_TAA_TCE  # TAA + TAA:TCE <-> TAA:TAA:TCE
    flux_CD3_TCE_binding = kon_cd3_tce * CD3 * TCE - koff_CD3 * CD3_TCE  # CD3 + TCE <-> CD3:TCE
    flux_TAA_trimer_binding_1 = lambda_trimer_tumor * kon_taa_tce * TAA * CD3_TCE - koff_TAA * CD3_TCE_TAA  # TAA + CD3:TCE <-> CD3:TCE:TAA 
    flux_TAA_trimer_binding_2 = lambda_trimer_T * kon_cd3_tce * CD3 * TAA_TCE - koff_CD3 * CD3_TCE_TAA  # CD3 + TAA:TCE <-> CD3:TCE:TAA 
    flux_TAA_tetramer_binding_1 = lambda_taa_avidity * kon_taa_tce * CD3_TCE_TAA * TAA - koff_TAA * CD3_TCE_TAA_TAA # CD3:TCE:TAA + TAA <-> CD3:TCE:TAA:TAA ## AKW Comment: it isn't truly a tetramer since the trimer is between the 2 cells and the molecule. The equations are fine but I'd consider a different name. YL: tetramer in the sense of the protein, not cells
    flux_TAA_tetramer_binding_2 = lambda_trimer_T * kon_cd3_tce * CD3 * TAA_TAA_TCE - koff_CD3 * CD3_TCE_TAA_TAA # TAA:TAA:TCE + CD3 <-> CD3:TCE:TAA:TAA
    # trimer & cell dynamics 
    trimer = (CD3_TCE_TAA + CD3_TCE_TAA_TAA)/nmol_per_mol*NAv*Vwell 
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
        Dt(TCE) ~ ( - flux_TAA_TCE_binding - flux_CD3_TCE_binding ), 
        # CD3 concentration in medium 
        Dt(CD3) ~ ( ksyn_CD3 - kdeg_CD3*CD3  # synthesis and degredation of CD3 on T cells 
                    - flux_CD3_TCE_binding - flux_TAA_trimer_binding_2 - flux_TAA_tetramer_binding_2  # binding activity 
                    + flux_tce_apop*(CD3_TCE_TAA + CD3_TCE_TAA_TAA) * (1-frac_act_T_death)  # CD3 release from trimer/tetramer after tumorlysis
                    + k_t_prolif*T_act*(CD3_per_cell*nmol_per_mol/NAv/Vwell) ),  # CD3 from T cell proliferation
        # TAA concentration in medium 
        Dt(TAA) ~ ( ksyn_TAA - kdeg_TAA*TAA  # synthesis and degredation of TAA on tumor cells 
                    - flux_TAA_TCE_binding - flux_TAA_TAA_TCE_binding - flux_TAA_trimer_binding_1 - flux_TAA_tetramer_binding_1  # binding activity 
                    - flux_tce_apop*TAA # TAA loss from tumor cell death 
                    + k_growth*Tumor_cell*(TAA_per_cell*nmol_per_mol)/NAv/Vwell ), # TAA from tumor cell proliferation 
        # CD3:TCE concentration in medium 
        Dt(CD3_TCE) ~ (flux_CD3_TCE_binding - flux_TAA_trimer_binding_1  # binding activity 
                        - kdeg_CD3*CD3_TCE), # CD3:TCE loss due to CD3 turnover 
        # TAA:TCE in medium 
        Dt(TAA_TCE) ~ (flux_TAA_TCE_binding - flux_TAA_TAA_TCE_binding - flux_TAA_trimer_binding_2  # binding activity 
                        - kdeg_TAA*TAA_TCE  # TAA:TCE loss due to TAA turnover
                        - flux_tce_apop*TAA_TCE ), # TAA:TCE loss tumorlysis
        # TAA:TAA:TCE in medium 
        Dt(TAA_TAA_TCE) ~ (flux_TAA_TAA_TCE_binding - flux_TAA_tetramer_binding_2  # binding activity 
                            - lambda_int*kdeg_TAA*TAA_TAA_TCE  # TAA:TAA:TCE loss due to TAA turnover 
                            - flux_tce_apop*TAA_TAA_TCE ),  # TAA:TAA:TCE loss tumorlysis
        # CD3:TCE:TAA in medium 
        Dt(CD3_TCE_TAA) ~ (flux_TAA_trimer_binding_1 + flux_TAA_trimer_binding_2 - flux_TAA_tetramer_binding_1  # binding activity 
                            - flux_tce_apop*CD3_TCE_TAA ),  # CD3:TCE:TAA loss tumorlysis
        # CD3:TCE:TAA:TAA in medium 
        Dt(CD3_TCE_TAA_TAA) ~ ( flux_TAA_tetramer_binding_1 + flux_TAA_tetramer_binding_2  # binding activity 
                                - flux_tce_apop*CD3_TCE_TAA_TAA ),    # CD3:TCE:TAA:TAA loss tumorlysis
        # Tumor cells in medium 
        Dt(Tumor_cell) ~ ( k_growth*Tumor_cell - flux_tce_apop*Tumor_cell  ), 
        # T cell in medium 
        Dt(T_cell) ~ (k_t_prolif*T_act - frac_act_T_death*flux_tce_apop*T_cell ), 
        Dt(T_act) ~ (k_t_prolif*T_act + flux_t_act*(T_cell - T_act)), 
        # cytokine 
        Dt(cytokine) ~ (k_syn_cytokine*(trimer_per_T^n_syn_cytokine)/(tt50_cytokine^n_syn_cytokine + trimer_per_T^n_syn_cytokine) * max(T_act, 0) - k_deg_cytokine*cytokine), 
        # dummy variable 
        Dt(deg_TAA_TCE) ~ ( kdeg_TAA*TAA_TCE + lambda_int*kdeg_TAA*TAA_TAA_TCE ), 
        Dt(dead_T) ~ ( frac_act_T_death*flux_tce_apop*T_cell ),
        Dt(dead_Tumor) ~ (flux_tce_apop*Tumor_cell),
    ]    

    ODESystem(eqs, t, vars, pars; name=name, observed = [
        (@variables R0_CD3(t))... ~ ((CD3_TCE + CD3_TCE_TAA + CD3_TCE_TAA_TAA)/(CD3 + CD3_TCE + CD3_TCE_TAA + CD3_TCE_TAA_TAA)), 
        (@variables R0_TAA(t))... ~ ((TAA_TCE + TAA_TAA_TCE + CD3_TCE_TAA_TAA)/(TAA + TAA_TCE + TAA_TAA_TCE + CD3_TCE_TAA_TAA)), 
        (@variables trimer_per_T(t))... ~ (trimer / max(T_cell, 1)), 
        (@variables trimer_per_Tumor(t))... ~ (trimer / max(Tumor_cell, 1)), 
        (@variables trimer_act_T_per_Tumor(t))... ~ (trimer / max(Tumor_cell, 1) * max(T_act, 0) / max(T_cell, 1) ), 
        (@variables trimer_1_per_tumor_cell(t))... ~ (CD3_TCE_TAA/nmol_per_mol*NAv*Vwell/max(Tumor_cell, 1)), 
        (@variables trimer_2_per_tumor_cell(t))... ~ (CD3_TCE_TAA_TAA/nmol_per_mol*NAv*Vwell/max(Tumor_cell, 1)), 
        ]);
end