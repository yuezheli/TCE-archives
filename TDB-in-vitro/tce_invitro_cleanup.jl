# date: 9/12/24
# author: Yuezhe Li 
# purpose of this code: a cleanup version of TCE model, used in ACoP 2023 workshop
# https://docs.google.com/presentation/d/1qgINuEx62biMEYuzXSd4ZVbo9W1hClWt9TEkEqh4tYU/edit?usp=sharing
# modification since the workshop: add a hill equation on killing function 

using DifferentialEquations, ModelingToolkit 

# define constants 
const nmol_per_mol = 1.0e9 
const NAv = 6.022e23 #  Avogadro's number (1/mole)
const s_per_hr = 3600.0 
const s_per_d = 3600.0*24.0 
const mm3_per_L = 1.0e6
const g_cell_per_l = 1000 #  density of tumor tissue; [g/L]
const Vwell = 2E-4 # Volume of well [L]
const Vc = 4/3 * pi * (4E-6)^3 * 1E3; # tumor cell volume; [L]

function TCE_invitro(T_cells_init, Tumor_cells_init, Tumor_cell_doubling_time; name)
    pars = @parameters begin 
        # Binding parameters
        kon = 7E-4                      # protein binding rate constant; [nM-1.s-1]; typical value; https://www.bio-rad.com/webroot/web/pdf/lsr/literature/Bulletin_6044A.pdf
        KD_CD3 = 28.03                  # CD3 Binding Affinity; [nM]  
        KD_TAA = 0.18                   # Tumor Associated Antigen (TAA) Binding Affinity; [nM] 
        lambda = 1                      # scaling factor for trimer formation
        # Receptor kinetics
        Thalf_CD3 = 0.5                 # Internalization half-life of CD3; [hr]
        Thalf_TAA = 5.8                # Internalization half-life of TAA; [hr]
        # Cell line parameters
        CD3_per_cell = 57000            # Expression of CD3 per T- cell (molecule/cell) 
        TAA_per_cell = 13173           # Expression of TAA per tumor cell (molecule/cell)
        # Tumor cell kinetics
        k_growth = log(2)/(Tumor_cell_doubling_time * s_per_hr)  # Tumor growth rate [s-1]
        k_apop = 6.0e-6                 # Rate of tumor cell apoptosis due to trimers (CD3:TCE:TAA) on tumor cell; [s-1]
        ec50_kill = 1E-4                # trimer per cell to achive 50% of maximum killing 
        n_kill = 1                      
    end
    # define variables
    @independent_variables t; 
    Dt = Differential(t); 
    vars = @variables begin 
        TCE(t) = 0.
        CD3(t) = CD3_per_cell*T_cells_init/NAv/Vwell * nmol_per_mol                     # init conc of CD3; [nM]  
        TAA(t) = TAA_per_cell*Tumor_cells_init/NAv/Vwell * nmol_per_mol               # init conc of TAA; [nM] 
        CD3_TCE(t) = 0.0
        TAA_TCE(t) = 0.0
        CD3_TCE_TAA(t) = 0.0
        Tumor_cell(t) = Tumor_cells_init
    end
    # calculate other parameters 
    koff_CD3 = kon * KD_CD3                                # [s-1]
    koff_TAA = kon * KD_TAA                              # [s-1]
    kdeg_CD3 = log(2)/(Thalf_CD3 * s_per_hr)                # CD3 degredation rate; [s-1]
    kdeg_TAA = log(2)/(Thalf_TAA * s_per_hr)              # TAA degredation rate; [s-1]
    ksyn_CD3 = kdeg_CD3 * (CD3_per_cell * T_cells_init / NAv * nmol_per_mol)/Vwell # CD3 synthesis rate per well; [nM.s-1]
    ksyn_TAA = kdeg_TAA * (TAA_per_cell * Tumor_cell / NAv * nmol_per_mol)/Vwell # TAA synthesis rate per well; [nM.s-1]
    # binding reactions
    flux_cd3_TCE_binding = kon * CD3 * TCE - koff_CD3 * CD3_TCE
    flux_TAA_TCE_binding = kon * TAA * TCE - koff_TAA * TAA_TCE
    flux_cd3_TAA_TCE_binding = lambda * kon * TAA_TCE * CD3 - koff_CD3 * CD3_TCE_TAA
    flux_TAA_cd3_TCE_binding = lambda * kon * CD3_TCE * TAA - koff_TAA * CD3_TCE_TAA
    # tumor killing
    trimer_per_tumor_cell = (CD3_TCE_TAA/nmol_per_mol*NAv*Vwell)/max(Tumor_cell, 1)
    flux_TCE_apop = k_apop * ( trimer_per_tumor_cell^n_kill / (trimer_per_tumor_cell^n_kill + ec50_kill^n_kill) ) * Tumor_cell
    # define diff eqs
    eqs = [
        Dt(TCE) ~ -flux_cd3_TCE_binding - flux_TAA_TCE_binding,
        Dt(CD3) ~ (ksyn_CD3 - kdeg_CD3*CD3) - flux_cd3_TCE_binding - flux_cd3_TAA_TCE_binding + flux_TCE_apop*(CD3_TCE_TAA/max(Tumor_cell, 1)), 
        Dt(TAA) ~ (ksyn_TAA - kdeg_TAA*TAA) -flux_TAA_TCE_binding - flux_TAA_cd3_TCE_binding - flux_TCE_apop*(TAA/max(Tumor_cell, 1)) + k_growth*Tumor_cell*(TAA_per_cell*nmol_per_mol)/(NAv*Vwell),
        Dt(CD3_TCE) ~ flux_cd3_TCE_binding - flux_TAA_cd3_TCE_binding - kdeg_CD3*CD3_TCE,
        Dt(TAA_TCE) ~ flux_TAA_TCE_binding - flux_cd3_TAA_TCE_binding - kdeg_TAA*TAA_TCE - flux_TCE_apop*(TAA_TCE/max(Tumor_cell, 1)),
        Dt(CD3_TCE_TAA) ~ flux_cd3_TAA_TCE_binding + flux_TAA_cd3_TCE_binding - flux_TCE_apop*(CD3_TCE_TAA/max(Tumor_cell, 1)),
        Dt(Tumor_cell) ~ k_growth*Tumor_cell - flux_TCE_apop, 
    ];
    ODESystem(eqs, t, vars, pars; name=name);
end


