# date: 4/8/2025 
# author: Yuezhe Li 
# purpose of this code: to incorporate v8 of in vitro model to Hosseini backbone 

function Mouse_PKPD(; name)
    @independent_variables t
    Dt = Differential(t)

    # default parameters for mouse 
    pars = @parameters begin

        frac_targpos = 1 # [description = "fraction of initially target positive cells] baseline value

        # initialization, dosing 
        BW = mouse_WT                                                           # mouse body weight [kg]
        dose_mgkg = 0.0                                                         # TCE dose administered [mg/kg]
        dose_nmol = dose_mgkg * (BW) / (mg_per_g) / (MW_TCE) * (nmol_per_mol)   # convert TCE dose to nmol

        # initialization, tumor 
        TV_mm3_0= 100.0                                                         # [unit = u"mm^3", description = "initial tumor volume"]
        TV_L0 = TV_mm3_0/mm3_per_L # [unit = u"L", description = "initial tumor volume"]
        f_tum_isf = 0.5 # [description = "interstitial fraction of tumor volume"]
        f_tum_cell = 0.375 # [description = "cellular fraction of tumor volume"]
        Nc1_0 = TV_L0*f_tum_cell/cell_vol*frac_targpos # [description = "initial total number of TAA+ tumor cells"]
        Ncneg_0 = (1 - frac_targpos)*TV_L0*f_tum_cell/cell_vol # [description = "initial total number of TAA- tumor cells"]
        V_TME_isf0 = f_tum_isf*TV_L0 # [unit = u"L", description = "volume of interstitial fluid in tumor"]

        # initialization, PBMC infusion 
        frac_Tcell_PBMC = 0.5 #unitless fraction of naive T cells in PBMC
        PBMC_dose = 1E7  #Protocol note: "Injected 10M tumor cells, leave until it can be measured by day 4/5 (tumor size of 100mm3), then inoculate 10M PBMC (about 50% T cells) and wait for 4 days to start TCE dosing. Unclear how many T cells survive."
        f_T_cell_survival = 0.5 # unitless fraction of surviving naive T cells post inoculation
        #Protocol note: "Injected 10M tumor cells, leave until it can be measured by day 4/5 (tumor size of 100mm3), then inoculate 10M PBMC (about 50% T cells) and wait for 4 days to start TCE dosing. Unclear how many T cells survive."
        T_cells_dose = f_T_cell_survival  * PBMC_dose * frac_Tcell_PBMC
        # Trpbref_perml = T_cells_dose/Vpb_mL #PBMC dosing replacing hosseini baseline conc of resting CD8+ T-cells in PB [cells/mL] 
        Trpbref_perml = 1 # ASSUMED 

        # tumor microenvironment & penetration 
        R_cap = 8 # [unit = u"um", description = "capillary radius"]
        R_Krogh = 75 # [unit = u"um", description = "Krogh cylinder radius"]
        perm_TCE = 50.34/hr_per_day # [unit = u"um/hr", description = "vascular permeability of TCE"] from Scheuher et al 2023 (https://link.springer.com/article/10.1007/s10928-023-09884-6) 
        diff_TCE = 1.74E5/hr_per_day # [unit = u"um^2/hr", description = "diffusivity of an TCE"] from Scheuher et al 2023 (https://link.springer.com/article/10.1007/s10928-023-09884-6)
        epsilon = 0.24 # partition coeffficient of TCE between plasma:tumor; Singh et al., 2016 # https://pmc.ncbi.nlm.nih.gov/articles/PMC5234806/

        # tumor physiology
        Tumor_doubling_time_d = 7       # [unit = u"d", description = "doubling time of untreated tumor"]
        TAA_per_cell = 13173            # Expression of TAA per tumor cell (molecule/cell) (placeholder) 
        kon_taa_tce = 1.2E-4 *s_per_hr  # protein binding rate constant; [nM-1.h-1]; (BI 1.1)
        KD_TAA = 20                     # TCE:TAA Binding Affinity; [nM] (BI 1.1)
        lambda_taa_avidity = 1          # scaling factor for cis-avidity of TAA:TCE binding; should be between 1E2 - 1E5; Harms et al., 2014; https://pubmed.ncbi.nlm.nih.gov/23872324/
        lambda_trimer_tumor = 1         # placeholder trans avidity for trimer formation on tumor cells 
        lambda_int = 0.002              # scaling factor for internalization rates between TAA:TAA:TCE and TAA:TCE (value set based on in vitro tuning)
        Thalf_TAA = 0.1                 # Internalization half-life of TAA; [hr]; tuned based on in vitro experiment 
        k_growth = log(2)/(Tumor_doubling_time_d*hr_per_day)  # Tumor growth rate [h-1]

        # binding parameters for TCE to CD3 in central compartment (keep for now - will be useful when PBMC are dosed)
        CD3_per_cell = 57000                                # Expression of CD3 per T cell (molecule/cell); https://pubmed.ncbi.nlm.nih.gov/8817090/; could range between 5000 - 70000 # https://pubmed.ncbi.nlm.nih.gov/19834968/
        kon_cd3_tce =3.1E-4 *s_per_hr                       # protein binding rate constant; [nM-1.h-1]; cyno (BI 1.1) 
        KD_CD3 = 230                                        # TCE:CD3 binding affinity; [nM] cyno (BI 1.1) 
        koff_CD3 = kon_cd3_tce * KD_CD3                     # [h-1]
        lambda_trimer_T = 1                                 # placeholder trans avidity for trimer formation on T cells 
        Thalf_CD3 = 3.5                                     # Internalization half-life of CD3; [hr]; https://journals.aai.org/jimmunol/article/173/1/384/72940/Constitutive-and-Ligand-Induced-TCR-Degradation-1
        kdeg_CD3 = log(2) / (Thalf_CD3)                     # CD3 degredation rate; [h-1]

        # tumor killing
        emax_apop = 3e-6*s_per_hr       # Max rate of tumor cell apoptosis due to trimers (CD3:TCE:TAA) on tumor cell; [h-1]
        tau = 25.8                      # [unit = u"hr", description = "delay in cell killing"] value from prior experience;
        tt50_kill = 0.01                # Number of trimers per tumor cell resulting in half-max cell killing; [#/tumor cell]
        n_kill = 1 
        FLG_trimer_per_tumor_kill = 1   # FLG for triggering tumor cell death based on trimers on tumor cell (if turned off, then trigger death based on trimer per T cell)
        
        # T cell activation & cytokine release
        emax_t_activate = 8.0e-6*s_per_hr  # maximum rate for trimer-induced T cell activation ; [h-1]
        tt50_t_activate = 0.01             # Number of trimers per T cell resulting in half-max cell killing; [#/naive T cell]
        n_act = 0.6                        # Hill coefficient, T cell activation 

        # cytokine production 
        k_syn_cytokine1 = 0.              # maximum rate for trimer-induced (per T cell) cytokine1 production [nmol/hr/cell] 
        k_syn_cytokine2 = 0.              # maximum rate for trimer-induced (per T cell) cytokine2 production [nmol/hr/cell] 
        k_syn_cytokine3 = 0.              # maximum rate for trimer-induced (per T cell) cytokine3 production [nmol/hr/cell]
        tt50_cytokine = 0.01              # EC50 of trimer-induced (per T cell) cytokine production 
        n_syn_cytokine = 1.               # Hill coefficient trimer-induced (per T cell) cytokine production 
        k_deg_cytokine1 = 0.              # degredation rate for cytokine [1/hr]
        k_deg_cytokine2 = 0.              # degredation rate for cytokine [1/hr]
        k_deg_cytokine3 = 0.              # degredation rate for cytokine [1/hr]
        k_tr_cytokine1 = 0.1              # distribution rate for cytokine1 from tumor to PB [1/hr]; # https://ascpt.onlinelibrary.wiley.com/doi/full/10.1002/psp4.12574
        k_tr_cytokine2 = 0.1              # distribution rate for cytokine2 from tumor to PB [1/hr]; # https://ascpt.onlinelibrary.wiley.com/doi/full/10.1002/psp4.12574
        k_tr_cytokine3 = 0.1              # distribution rate for cytokine3 from tumor to PB [1/hr]; # https://ascpt.onlinelibrary.wiley.com/doi/full/10.1002/psp4.12574
        k_cl_cytokine1 = log(2)/0.33      # clearance of cytokine1 in PB [1/hr] # assumed half life ~ 20 mins
        k_cl_cytokine2 = log(2)/0.33      # clearance of cytokine2 in PB [1/hr]
        k_cl_cytokine3 = log(2)/0.33      # clearance of cytokine3 in PB [1/hr]
   
        ## T cells - not relevant in mouse. 
        kTprolif = 0.7                             # fraction of post-activated cells that proliferate 
        fTaprolif = 2.11/hr_per_day                # proliferation rate of post activated CD8+ T-cells [1/h] -- See Hosseini appendix page 3; # https://static-content.springer.com/esm/art%3A10.1038%2Fs41540-020-00145-7/MediaObjects/41540_2020_145_MOESM1_ESM.pdf
        kTaexit = 0.12/hr_per_day                  # trafficking rate of activated CD8+ T-cells from PB to tissue [1/h] 
        kTact = 0.0                                # T cell activation rate outside tumor -- was 9.83 day-1 in Hosseini 
        fTadeact = 0.01                            # fraction of activated CD8+ T-cells that deactivate to resting or "exhausted" CD8+ T-cells
        fTap = 3.77                                # Ratio of the partition coeff for activated CD8+ T-cells to that of resting CD8+ T-cells
        kTaapop = 0.06/hr_per_day                  # apoptosis rate of CD69+ CD8+ T-cells (activated CD8+ T-cells) [1/h]
        fTrapop = 0.2                              # Ratio of the apoptosis rate of resting CD8+ T-cells to that of CD69+CD8+ T-cells
        kTrexit = 0.05/hr_per_day                  # trafficking rate of resting/ post-activated CD8+ T-cells from PB to tissue [1/h]
        fAICD = 1.5                                # Effect of activation-induced cell death (AICD) on CD69+CD8+ T-cells
        fTa0apop = 2.02                            # Ratio of the apoptosis rate of post-activated CD8+ T-cells to that of CD69+CD8+ T-cells
        fTa0deact = 0.0/hr_per_day                 # Conversion rate of post-activated CD8+ T-cells to resting CD8+ T-cells [1/h]
        kTgen = 1.0/hr_per_day                     # for generation of resting (naive) T-cells in PB [1/h]   
        # VmT = 0.9                                # maximum rate of T-cell activation [unitless] in Hosseini model; replaced by in vitro mechanistic model
        fTact = 0.25                               # VmT is multiplied by this factor (<1) to capture the less efficient activation of tissue T-cells

        # additional T cell dynamics parameters (in tumor only)
        k_t_n_prolif = log(2)/(48)      # naive T cell proliferation rate [h-1]; 
        k_t_a_prolif = fTaprolif        # active T cell proliferate rate [hr-1]; 
        frac_act_T_death = 0            # fraction of active T cell death after trimer formation (turned off in vitro) 
        kTaapop_tumor = kTaapop
        
        # organ volumes 
        Vpb_mL = 0.944                     # physiological volume of peripheral blood [mL] (Shah& Betts, 2012)
        Vtissue_mL = 0.127                 # physiological volume of spleen [mL](Shah& Betts, 2012)
        Vtissue2_mL = 0.113                # total physiological volume of LNs [mL](Shah& Betts, 2012)
        Vtissue3_mL = 0.35                 # Physiology volume of BM [mL]; https://pmc.ncbi.nlm.nih.gov/articles/PMC5738992/ 
        
        ## other 
        tinjhalf = 1. * hr_per_day         # placeholder from Hosseini model describing the T cell mobilization after injection; not used in mouse

        # FLGs (internal decision on 4/9/25 -- not to include spleen, LN, BM in mouse)
        act0on = 0.
        tissue1on = 0.                  
        tissue2on = 0.
        tissue3on = 0.
        tumor_on = 1. 
        FLG_T_death_PB = 1
        FLG_tce_driven_T_cell_enter_tumor = 0

        # T cell entering tumor through TCE -driven manner (Ma et al., 2020) # https://pubmed.ncbi.nlm.nih.gov/33615174/
        q_t_in = 5.58E-5*s_per_hr # [hr-1.L-1.nM-1]

        # Partition coeff: ratio of the baseline conc. of resting CD8+ T-cells in spleen to that in PB
        KTrp = 0.5 # human default should not apply here ....500.                     
        KTrp2 = 0.5 # human default should not apply here ....500.
        KTrp3 = 0.5 # human default should not apply here ....50
        KTrptumor = 0.5 # human default should not apply here ....125. 
        
        pb_T_cells_init = T_cells_dose  
        tissue_T_cells_init = 0
        tissue2_T_cells_init = 0 
        tissue3_T_cells_init = 0
        tumor_T_cells_init = 0

        ## PK   
        CL_TCE = 0.5e-5             # [L/h] Cao et al, 2014 Table 1
        Q_TCE =  0.5/hr_per_day/mL_per_L            # [L/h] fit to mouse PK
        V1_TCE =  Vpb_mL/mL_per_L            # [L]   (Shah&Betts, 2012)(value in Cao et al, 2014 Table 1 is 0.85mL), 1.2 in Singh
        V2_TCE = 0.0948*BW            # [L] 2.6544mL from Singh et al, 2018 173.47*BW/mL_per_L [L] assumed to scale as in Hosseini = 4.85716mL;   Ma et al. 2021 uses 36mL 
    
    end 

    vars = @variables begin 
        ## T cells in PB
        restTpb(t) = pb_T_cells_init #Trpbref_perml * Vpb_mL
        actTpb(t) = 0. 
        act0Tpb(t) = 0.
        ## T cells in spleen
        restTtiss(t) = tissue_T_cells_init #tissue1on * KTrp * Trpbref_perml * Vtissue_mL
        actTtiss(t) = 0.
        act0Ttiss(t) = 0.
        ## T cells in LN
        restTtiss2(t) = tissue2_T_cells_init #tissue2on * KTrp2 * Trpbref_perml * Vtissue2_mL
        actTtiss2(t) = 0.
        act0Ttiss2(t) = 0.
        ### T cells in bone marrow
        restTtiss3(t) = tissue3_T_cells_init
        actTtiss3(t) = 0.
        act0Ttiss3(t) = 0.
        # T cells in tumor
        restTtumor(t) = tumor_T_cells_init 
        actTtumor(t) = 0.
        act0Ttumor(t) = 0.

        ## PK 
        TCEc_nM(t) = dose_nmol/V1_TCE
        TCEp_nM(t) = 0.
   
        # injection effect
        TCEinjection_effect(t) = 0.
        
        #plasma compartment
        CD3c_nM(t) = restTpb*CD3_per_cell/NAv*nmol_per_mol/V1_TCE            # init conc of CD3 in central cmpt; [nM] 
        CD3c_TCEc_nM(t) = 0.0

        #tumor compartment
        # TCEt_nM(t) = tumor_on*Kptumor*max(dose_nmol/V1_TCE , 0.)  #  this was used in Hosseini model; updated with Krogh cylinder model for TCE penetration into solid tumor
        TCEt_nM(t) = 0.
        CD3_t(t) = CD3_per_cell*tumor_T_cells_init /NAv/V_TME_isf0 * nmol_per_mol             # init conc of total CD3; [nM]  
        TAA_t(t) = TAA_per_cell*Nc1_0/NAv/V_TME_isf0 * nmol_per_mol                           # init conc of TAA; [nM] 
        CD3_TCE_t(t) = 0.0
        TAA_TCE_t(t) = 0.0
        TAA_TAA_TCE_t(t) = 0.0
        CD3_TCE_TAA_t(t) = 0.0
        CD3_TCE_TAA_TAA_t(t) = 0.0

        cytokine1_t(t) = 0.     # cytokine concentration in tumor interstitium [nM] 
        cytokine2_t(t) = 0.     # 
        cytokine3_t(t) = 0.     # 

        cytokine1_pb(t) = 0.  # cytokine concenrtation in PB [nM]
        cytokine2_pb(t) = 0.
        cytokine3_pb(t) = 0.
        
        Nc1(t) = Nc1_0 # [description = "number of TAA+ cells in stage 1"] (source: BI CDX mouse protocol)
        Nc2(t) = 0.0 # [description = "number of TAA+ cells in stage 2"] 
        Nc3(t) = 0.0 # [description = "number of TAA+ cells in stage 3"] 
        Nc4(t) = 0.0 # [description = "number of TAA+ cells in stage 4"] 
        Nc_neg(t) = Ncneg_0 # [description = "number of TAA- cells"] 
    end

    # tumor volume (L)
    Nctot = Nc1 + Nc2 + Nc3 + Nc4 + Nc_neg
    TV_L = (Nctot/f_tum_cell) * cell_vol
    Vtumor_mL = TV_L*mL_per_L
    V_TME_isf = f_tum_isf*TV_L

    # CD3
    ksyn_CD3 = kdeg_CD3 * (CD3_per_cell *  restTpb / NAv * nmol_per_mol)/V1_TCE # CD3 synthesis rate in central cmpt; [nM.h-1]
    ksyn_CD3t = kdeg_CD3 * (CD3_per_cell * (restTtumor+actTtumor) / NAv * nmol_per_mol)/V_TME_isf # CD3 synthesis rate in tumor compartment on T cells; [nM.h-1] 
    
    

    # TAA related parameters 
    koff_TAA = kon_taa_tce * KD_TAA                # [h-1]
    kdeg_TAA = log(2)/(Thalf_TAA )       # TAA degredation rate; [h-1]
    ksyn_TAA = kdeg_TAA * (TAA_per_cell * Nc1 / NAv * nmol_per_mol)/V_TME_isf # TAA synthesis rate per well; [nM.s-1]    
   
    # tumor PK as in Hosseini
    #TCEt_nM = tumor_on*Kptumor*max(TCEc_nM , 0.) #assume partition coeff same for molar vs weight units --> converted to tumor diff equation below

    # tumor radius (mm) 
    R_tumor_mm = (3.0*(maximum([TV_L,0.0])*mm3_per_L)/(4.0*pi))^(1/3)
    # tumor radius (um)
    R_tumor_um = R_tumor_mm*um_per_mm
    # lymph flow for tumor
    tumor_penetration = (2*perm_TCE*R_cap/(R_Krogh^2)) + (6*diff_TCE/(R_tumor_um^2))
    flux_tce_plasma2tumor = tumor_penetration * (TCEc_nM*epsilon - TCEt_nM)

    # T cells
    restTpb_perml = restTpb/Vpb_mL
    act0Tpb_perml = act0Tpb/Vpb_mL
    actTpb_perml = actTpb/Vpb_mL
    totTpb_perml = restTpb_perml + act0Tpb_perml + actTpb_perml

    # binding dynamics of TAA, TCE, CD3 in tumor (copied from in vitro model)
    flux_TAA_TCE_binding = kon_taa_tce * TAA_t * TCEt_nM - koff_TAA * TAA_TCE_t  # TAA + TCE <-> TAA:TCE
    flux_TAA_TAA_TCE_binding = lambda_taa_avidity * kon_taa_tce * TAA_TCE_t * TAA_t - koff_TAA * TAA_TAA_TCE_t  # TAA + TAA:TCE <-> TAA:TAA:TCE
    flux_CD3_TCE_binding = kon_cd3_tce * CD3_t * TCEt_nM - koff_CD3 * CD3_TCE_t  # CD3 + TCE <-> CD3:TCE
    flux_TAA_trimer_binding_1 = lambda_trimer_tumor * kon_taa_tce * TAA_t * CD3_TCE_t - koff_TAA * CD3_TCE_TAA_t  # TAA + CD3:TCE <-> CD3:TCE:TAA 
    flux_TAA_trimer_binding_2 = lambda_trimer_T * kon_cd3_tce * CD3_t * TAA_TCE_t - koff_CD3 * CD3_TCE_TAA_t  # CD3 + TAA:TCE <-> CD3:TCE:TAA 
    flux_TAA_tetramer_binding_1 = lambda_taa_avidity * kon_taa_tce * CD3_TCE_TAA_t * TAA_t - koff_TAA * CD3_TCE_TAA_TAA_t # CD3:TCE:TAA + TAA <-> CD3:TCE:TAA:TAA 
    flux_TAA_tetramer_binding_2 = lambda_trimer_T * kon_cd3_tce * CD3_t * TAA_TAA_TCE_t - koff_CD3 * CD3_TCE_TAA_TAA_t # TAA:TAA:TCE + CD3 <-> CD3:TCE:TAA:TAA

    # trimer & cell dynamics (adapted for in vivo model)
    trimer = (CD3_TCE_TAA_t + CD3_TCE_TAA_TAA_t)/nmol_per_mol*NAv*V_TME_isf # number
    trimer_per_Tu_cell = max(trimer / max(Nc1, 1), 0) # number
    trimer_per_Tu_cell_on_act_T = trimer_per_Tu_cell * max(actTtumor, 0) / max(restTtumor + actTtumor + act0Ttumor, 1)
    trimer_per_T = max(trimer / max((actTtumor+restTtumor+act0Ttumor), 1), 0)  # average trimer per T; this should be the same for restT, actT, act0T
    trimer_per_Ta = trimer_per_T
    trimer_per_Tr = trimer_per_T

    # the hill coefficient part related to cytokine 
    hill_cytokines_tumor = (trimer_per_Ta^n_syn_cytokine)/(tt50_cytokine^n_syn_cytokine + trimer_per_Ta^n_syn_cytokine)

    # tumor killing (trimer driven)
    flux_tce_apop_tumor = emax_apop * (trimer_per_Tu_cell_on_act_T^n_kill) / (tt50_kill^n_kill + trimer_per_Tu_cell_on_act_T^n_kill)
    flux_tce_apop_T = emax_apop * (trimer_per_T^n_kill) / (tt50_kill^n_kill + trimer_per_T^n_kill) * max(actTtumor, 0) / max(actTtumor+restTtumor+act0Ttumor, 1) 
    flux_tce_apop = flux_tce_apop_tumor * FLG_trimer_per_tumor_kill + flux_tce_apop_T * (1 - FLG_trimer_per_tumor_kill)
    
    # T cell activation (trimer driven)
    flux_t_act = emax_t_activate * (trimer_per_Tr^n_act)/(tt50_t_activate^n_act + trimer_per_Tr^n_act)

    # TCE-induced T cell activation (set to 0 because there is no TCE-induced T cell activation in wt mouse_pk_tce)
    drugTpbact = 0
    drugTtissact = 0
    drugTtissact2 = 0
    drugTtissact3 = 0
    #drugTtumoract = fTact*VmT*(max(BTrratio_tumor, 0)^S/(KmBT_act^S+max(BTrratio_tumor, 0)^S))*((TDBtumor_ugperml*1000)^ndrugactT/(KdrugactT^ndrugactT+(TDBtumor_ugperml*1000)^ndrugactT))
    drugTtumoract = flux_t_act

    # T cell activation (all activation were turned off through kTact =  0)
    T_act_rest_PB = (kTact*(drugTpbact*restTpb-fTadeact*actTpb)); 
    T_act_rest_spleen = kTact*(drugTtissact*restTtiss-fTadeact*actTtiss); 
    T_act_rest_LN = (kTact*(drugTtissact2*restTtiss2-fTadeact*actTtiss2)); 
    T_act_rest_BM = (kTact*(drugTtissact3*restTtiss3-fTadeact*actTtiss3)); 
    # T_act_rest_Tumor = (kTact*(drugTtumoract*restTtumor-fTadeact*actTtumor)); #Hosseini
    T_act_rest_Tumor = flux_t_act * restTtumor - fTadeact*actTtumor;

    T_act_postact_PB = (act0on*kTact*(drugTpbact*act0Tpb-fTadeact*actTpb)); 
    T_act_postact_spleen = (act0on*kTact*(drugTtissact*act0Ttiss-fTadeact*actTtiss));
    T_act_postact_LN = (act0on*kTact*(drugTtissact2*act0Ttiss2-fTadeact*actTtiss2)); 
    T_act_postact_BM = (act0on*kTact*(drugTtissact3*act0Ttiss3-fTadeact*actTtiss3)); 
    T_act_postact_tumor = (act0on*kTact*(drugTtumoract*act0Ttumor-fTadeact*actTtumor));

    # T cell conversion from post-activated to resting state - we decided to omit
    T_postact_rest_conv_PB = (act0on*fTa0deact*act0Tpb)
    T_postact_rest_conv_spleen = (act0on*fTa0deact*act0Ttiss); 
    T_postact_rest_conv_LN = (act0on*fTa0deact*act0Ttiss2); 
    T_postact_rest_conv_BM = (act0on*fTa0deact*act0Ttiss3);
    T_postact_rest_conv_Tumor = (act0on*fTa0deact*act0Ttumor);

    # T cell trafficking 
    rest_T_enter_spleen = (tissue1on*kTrexit*Vpb_mL*((1+TCEinjection_effect)*restTpb/Vpb_mL*KTrp-restTtiss/Vtissue_mL));
    rest_T_enter_LN = (tissue2on*kTrexit*Vpb_mL*((1+TCEinjection_effect)*restTpb/Vpb_mL*KTrp2-restTtiss2/Vtissue2_mL));
    rest_T_enter_BM = (tissue3on*kTrexit*Vpb_mL*((1+TCEinjection_effect)*restTpb/Vpb_mL*KTrp3-restTtiss3/Vtissue3_mL));
    rest_T_enter_Tumor = tumor_on*(kTrexit*Vpb_mL*((1+TCEinjection_effect)*restTpb/Vpb_mL*KTrptumor-restTtumor/Vtumor_mL) + FLG_tce_driven_T_cell_enter_tumor * TCEc_nM * q_t_in * TV_L * restTpb) ;
    T_act_traffick_spleen = (tissue1on*kTaexit*Vpb_mL*(actTpb/Vpb_mL*KTrp*fTap-actTtiss/Vtissue_mL)); 
    T_act_traffick_LN = (tissue2on*kTaexit*Vpb_mL*(actTpb/Vpb_mL*KTrp2*fTap-actTtiss2/Vtissue2_mL));
    T_act_traffick_BM = (tissue3on*kTaexit*Vpb_mL*(actTpb/Vpb_mL*KTrp3*fTap-actTtiss3/Vtissue3_mL));
    T_act_traffick_Tumor = (tumor_on*kTaexit*Vpb_mL*(actTpb/Vpb_mL*KTrptumor*fTap-actTtumor/Vtumor_mL));
    postact_T_traffick_spleen = (tissue1on*act0on*kTrexit*Vpb_mL*((1+TCEinjection_effect)*act0Tpb/Vpb_mL*KTrp-act0Ttiss/Vtissue_mL)); 
    postact_T_traffick_LN = (tissue2on*act0on*kTrexit*Vpb_mL*((1+TCEinjection_effect)*act0Tpb/Vpb_mL*KTrp2-act0Ttiss2/Vtissue2_mL));
    postact_T_traffick_BM = (tissue3on*act0on*kTrexit*Vpb_mL*((1+TCEinjection_effect)*act0Tpb/Vpb_mL*KTrp3-act0Ttiss3/Vtissue3_mL)); 
    postact_T_traffick_Tumor = (tumor_on*act0on*kTrexit*Vpb_mL*((1+TCEinjection_effect)*act0Tpb/Vpb_mL*KTrptumor-act0Ttumor/Vtumor_mL));

    dTV_control = 1/V_TME_isf*( ((k_growth*(Nc1+Nc_neg)  - flux_tce_apop*Nc1  - (1/tau)*(Nc2 + Nc3+ Nc4)) /f_tum_cell ) * cell_vol)*f_tum_isf   # volume adjustment (YL correction)

    dTCEc_dt = (- 1/V1_TCE * Q_TCE * (TCEc_nM - TCEp_nM) - CL_TCE/V1_TCE * TCEc_nM - kon_cd3_tce * CD3c_nM * TCEc_nM + koff_CD3 * CD3c_TCEc_nM ) # this is plasma TCE concentration change 

    # track CD3 across naive T, active T, and exhausted T in tumor algebracially 
    CD3_t_naive = CD3_t * restTtumor / (restTtumor + actTtumor + act0Ttumor)
    CD3_t_active = CD3_t * actTtumor / (restTtumor + actTtumor + act0Ttumor)
    CD3_t_postactive = CD3_t * act0Ttumor / (restTtumor + actTtumor + act0Ttumor)

    # CD3 change in tumor due to cell proliferation
    cd3_per_T_in_tumor = ( CD3_per_cell*nmol_per_mol/NAv/V_TME_isf )
    cd3_naive_t_prolif = k_t_n_prolif * restTtumor * cd3_per_T_in_tumor
    cd3_active_t_prolif = k_t_a_prolif * actTtumor * cd3_per_T_in_tumor
    cd3_postactive_t_prolif = act0on*fTaprolif*kTprolif*actTtumor*cd3_per_T_in_tumor 
    cd3_active_t_death = kTaapop_tumor*actTtumor*cd3_per_T_in_tumor
    # CD3 change in tumor due to cell traffickng 
    cd3_naive_T_traffic = rest_T_enter_Tumor * cd3_per_T_in_tumor
    cd3_active_T_traffic = T_act_traffick_Tumor * cd3_per_T_in_tumor
    cd3_postactive_T_traffic = postact_T_traffick_Tumor * cd3_per_T_in_tumor
    # CD3 change due to T cell death 
    cd3_postactive_T_death = act0on*fTa0apop*kTaapop*act0Ttumor*cd3_per_T_in_tumor

    eqs = [
        # T cells dynamics in tumor
        Dt(restTtumor) ~ (T_postact_rest_conv_Tumor - T_act_rest_Tumor + rest_T_enter_Tumor), 
        # Dt(actTtumor) ~ (T_act_postact_tumor + T_act_rest_Tumor  + T_act_traffick_Tumor), # Hosseini 
        Dt(actTtumor) ~ (T_act_postact_tumor + T_act_rest_Tumor  + T_act_traffick_Tumor + k_t_a_prolif * actTtumor - frac_act_T_death*flux_tce_apop*actTtumor - kTaapop_tumor*actTtumor),  # incorporate placeholder active T cell proliferation & death
        Dt(act0Ttumor) ~ ( 0 - T_postact_rest_conv_Tumor - T_act_postact_tumor + (act0on*fTaprolif*kTprolif*actTtumor) + postact_T_traffick_Tumor - (fTa0apop*kTaapop*act0Ttumor) ), 
        
        # 2 compartment TCE PK (copied from cyno)
        Dt(TCEc_nM) ~ dTCEc_dt - flux_tce_plasma2tumor * V_TME_isf/V1_TCE, 
        Dt(TCEp_nM) ~ ( 1/V2_TCE * Q_TCE * (TCEc_nM - TCEp_nM) ), 
        
        # 2 compartment TCE PK TMDD via T cell CD3 (same as cyno)
        Dt(CD3c_nM) ~ ( ksyn_CD3 - kdeg_CD3*CD3c_nM - kon_cd3_tce * CD3c_nM * TCEc_nM + koff_CD3 * CD3c_TCEc_nM ), 
        Dt(CD3c_TCEc_nM) ~ (kon_cd3_tce * CD3c_nM * TCEc_nM - koff_CD3 * CD3c_TCEc_nM - kdeg_CD3*CD3c_TCEc_nM), 

        # tumor PK  - TCE concentration in tumor interstitium
        Dt(TCEt_nM) ~ ( - flux_TAA_TCE_binding - flux_CD3_TCE_binding + flux_tce_plasma2tumor - TCEt_nM*dTV_control ) * tumor_on, 
        
        # CD3 concentration in tumor interstitium
        Dt(CD3_t) ~ ( ksyn_CD3t  - kdeg_CD3*CD3_t  # synthesis and degredation of CD3 on T cells 
                    - flux_CD3_TCE_binding - flux_TAA_trimer_binding_2 - flux_TAA_tetramer_binding_2  # binding activity 
                    + flux_tce_apop * actTtumor * cd3_per_T_in_tumor * (1-frac_act_T_death)  # CD3 release from trimer/tetramer after tumorlysis
                    # + (k_t_n_prolif*restTtumor + fTaprolif*actTtumor)*(CD3_per_cell*nmol_per_mol/NAv/V_TME_isf )  # CD3 from T cell proliferation in vitro (this is replaced in mouse)
                    + cd3_naive_t_prolif + cd3_active_t_prolif + cd3_postactive_t_prolif - cd3_active_t_death
                    + cd3_naive_T_traffic + cd3_active_T_traffic + cd3_postactive_T_traffic 
                    - cd3_postactive_T_death 
                    - CD3_t*dTV_control ),
        # TAA concentration in tumor interstitium
        Dt(TAA_t) ~ ( ksyn_TAA - kdeg_TAA*TAA_t  # synthesis and degredation of TAA on tumor cells 
                    - flux_TAA_TCE_binding - flux_TAA_TAA_TCE_binding - flux_TAA_trimer_binding_1 - flux_TAA_tetramer_binding_1  # binding activity 
                    - flux_tce_apop*TAA_t # TAA loss from tumor cell death 
                    + k_growth*Nc1*(TAA_per_cell*nmol_per_mol)/NAv/V_TME_isf # TAA from tumor cell proliferation 
                    - TAA_t*dTV_control ), 
        # CD3:TCE concentration in tumor interstitium 
        Dt(CD3_TCE_t) ~ (flux_CD3_TCE_binding - flux_TAA_trimer_binding_1  # binding activity 
                        - kdeg_CD3*CD3_TCE_t # CD3:TCE loss due to CD3 turnover 
                        - kTaapop_tumor * CD3_TCE_t # CD3:TCE loss due to active T cell turnover 
                        - CD3_TCE_t*dTV_control ),
        # TAA:TCE in tumor interstitium
        Dt(TAA_TCE_t) ~ (flux_TAA_TCE_binding - flux_TAA_TAA_TCE_binding - flux_TAA_trimer_binding_2  # binding activity 
                        - kdeg_TAA*TAA_TCE_t  # TAA:TCE loss due to TAA turnover
                        - flux_tce_apop*TAA_TCE_t # TAA:TCE loss tumorlysis
                        - TAA_TCE_t*dTV_control ),
        # TAA:TAA:TCE in tumor interstitium
        Dt(TAA_TAA_TCE_t) ~ (flux_TAA_TAA_TCE_binding - flux_TAA_tetramer_binding_2  # binding activity 
                            - lambda_int*kdeg_TAA*TAA_TAA_TCE_t  # TAA:TAA:TCE loss due to TAA turnover 
                            - flux_tce_apop*TAA_TAA_TCE_t   # TAA:TAA:TCE loss tumorlysis
                            - TAA_TAA_TCE_t*dTV_control ),
        # CD3:TCE:TAA in tumor interstitium 
        Dt(CD3_TCE_TAA_t) ~ (flux_TAA_trimer_binding_1 + flux_TAA_trimer_binding_2 - flux_TAA_tetramer_binding_1  # binding activity 
                            - flux_tce_apop*CD3_TCE_TAA_t  # CD3:TCE:TAA loss tumorlysis
                            - kTaapop_tumor * CD3_TCE_TAA_t # CD3:TCE:TAA loss due to active T cell turnover 
                            - CD3_TCE_TAA_t*dTV_control ),
        # CD3:TCE:TAA:TAA in tumor interstitium
        Dt(CD3_TCE_TAA_TAA_t) ~ ( flux_TAA_tetramer_binding_1 + flux_TAA_tetramer_binding_2  # binding activity 
                                - flux_tce_apop*CD3_TCE_TAA_TAA_t     # CD3:TCE:TAA:TAA loss tumorlysis
                                - kTaapop_tumor * CD3_TCE_TAA_TAA_t   # CD3:TCE:TAA:TAA loss due to active T cell turnover 
                                - CD3_TCE_TAA_TAA_t*dTV_control ),

        # cytokine in tumor interstitium
        Dt(cytokine1_t) ~ (k_syn_cytokine1*hill_cytokines_tumor*actTtumor/V_TME_isf - k_deg_cytokine1*cytokine1_t - k_tr_cytokine1*cytokine1_t) - cytokine1_t*dTV_control ,  
        Dt(cytokine2_t) ~ (k_syn_cytokine2*hill_cytokines_tumor*actTtumor/V_TME_isf - k_deg_cytokine2*cytokine2_t - k_tr_cytokine2*cytokine2_t) - cytokine2_t*dTV_control ,
        Dt(cytokine3_t) ~ (k_syn_cytokine3*hill_cytokines_tumor*actTtumor/V_TME_isf - k_deg_cytokine3*cytokine3_t - k_tr_cytokine3*cytokine3_t) - cytokine3_t*dTV_control ,

        # cytokine in pb
        Dt(cytokine1_pb) ~ ( k_tr_cytokine1*cytokine1_t * TV_L * mL_per_L / Vpb_mL ) - k_cl_cytokine1*cytokine1_pb, 
        Dt(cytokine2_pb) ~ ( k_cl_cytokine2*cytokine2_t * TV_L * mL_per_L / Vpb_mL ) - k_cl_cytokine2*cytokine2_pb, 
        Dt(cytokine3_pb) ~ ( k_cl_cytokine3*cytokine3_t * TV_L * mL_per_L / Vpb_mL ) - k_cl_cytokine3*cytokine3_pb, 
 
        # other 
        Dt(TCEinjection_effect) ~ (-(log(2)/tinjhalf*TCEinjection_effect)),

        # Tumor cells
        Dt(Nc1) ~ ( k_growth*Nc1 - flux_tce_apop*Nc1 ), 
        Dt(Nc2) ~ (flux_tce_apop*Nc1 - 1/tau*Nc2),
        Dt(Nc3) ~ 1/tau*(Nc2 - Nc3),
        Dt(Nc4) ~ 1/tau*(Nc3 - Nc4),

        Dt(Nc_neg) ~ ( k_growth*Nc_neg),

        # T cell dynamics in PB (copied from Hosseini model)
        Dt(restTpb) ~ ( T_postact_rest_conv_PB - rest_T_enter_spleen - rest_T_enter_LN - rest_T_enter_BM - rest_T_enter_Tumor - T_act_rest_PB + 
                        (kTgen*Vpb_mL*Trpbref_perml*(max(0,1-totTpb_perml/Trpbref_perml)^2)) - FLG_T_death_PB*fTrapop*kTaapop*max(0,restTpb-Trpbref_perml*Vpb_mL)  ), 
        Dt(actTpb) ~ ( T_act_rest_PB + T_act_postact_PB - T_act_traffick_spleen - T_act_traffick_LN - T_act_traffick_BM - T_act_traffick_Tumor - (FLG_T_death_PB*kTaapop*(actTpb+fAICD*actTpb^2/(Vpb_mL*Trpbref_perml))) ), 
        
        # T cells (copied from Hosseini model -- NOT USED IN THE MOUSE MODEL)
        Dt(act0Tpb) ~ ( (act0on*fTaprolif*kTprolif*actTpb) - (FLG_T_death_PB*fTa0apop*kTaapop*act0Tpb) - T_postact_rest_conv_PB - T_act_postact_PB
                         - postact_T_traffick_spleen - postact_T_traffick_LN - postact_T_traffick_BM - postact_T_traffick_Tumor),
        Dt(restTtiss) ~ ( 0 - T_act_rest_spleen + rest_T_enter_spleen + T_postact_rest_conv_spleen ), 
        Dt(actTtiss) ~ ( T_act_rest_spleen + T_act_traffick_spleen + T_act_postact_spleen ), 
        Dt(act0Ttiss) ~ ( -(fTa0apop*kTaapop*act0Ttiss) - T_act_postact_spleen - T_postact_rest_conv_spleen + postact_T_traffick_spleen + (act0on*fTaprolif*kTprolif*actTtiss)), 
        Dt(act0Ttiss2) ~ ((act0on*fTaprolif*kTprolif*actTtiss2) + postact_T_traffick_LN - T_postact_rest_conv_LN - T_act_postact_LN - (fTa0apop*kTaapop*act0Ttiss2)), 
        Dt(restTtiss2) ~ (T_postact_rest_conv_LN + rest_T_enter_LN - T_act_rest_LN ), 
        Dt(actTtiss2) ~ ( T_act_postact_LN + T_act_traffick_LN + T_act_rest_LN ), 
        Dt(restTtiss3) ~ (rest_T_enter_BM + T_postact_rest_conv_BM - T_act_rest_BM ), 
        Dt(act0Ttiss3) ~ (postact_T_traffick_BM - (fTa0apop*kTaapop*act0Ttiss3) - T_act_postact_BM - T_postact_rest_conv_BM + (act0on*fTaprolif*kTprolif*actTtiss3)), 
        Dt(actTtiss3) ~ (T_act_traffick_BM + T_act_postact_BM + T_act_rest_BM ), 
    ]; 

    ODESystem(eqs, t, vars, pars; name=name, observed=[
        (@variables TCE_plasma_ngpermL(t))... ~ (TCEc_nM / (nmol_per_mol) * (MW_TCE) * (ng_per_g) / (mL_per_L))
        (@variables TV_mm3(t))... ~ (TV_L * (mm3_per_L))
        (@variables Tri_per_Tu(t))... ~ (trimer_per_Tu_cell)
        (@variables Trimer_per_Tu_cell_on_act_T(t))... ~ (trimer_per_Tu_cell_on_act_T)  
        (@variables Tri_per_T(t))... ~ (trimer_per_T)
        (@variables trimer_num(t))... ~ (trimer)
        (@variables Ntot(t))... ~ (Nctot)
        (@variables total_e2t(t))... ~ ((restTtumor + actTtumor) / Nctot)
    ])
     
end