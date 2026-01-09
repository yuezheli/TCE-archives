# author: Yuezhe LI
# date: 10/3/2023
# purpose of this code: to define parameters for bs-lonca-3.jl
# this parameter was based on the cleaned up model from bs-lonca-2.jl

using ComponentArrays
using Parameters: @unpack

const DayToHour = 24.
const HourToMinute = 60.

# define tumor cell volume
const Vc = 4/3 * pi * (4E-6)^3 * 1E6; # [mL]

# define plasma volume 
const V_Plasma = 3.126; # [L]

# define default parameters 
p_homo_3 = ComponentArray(
    # TDB : T-cell act
    VmT = 0.9,                      # maximum rate of T-cell activation [unitless]
    fTact = 0.25,                   # VmT is multiplied by this factor (<1) to capture the less efficient activation of tissue T-cells
    KmBT_act = 0.72,                # The B/T rario, at which the rate is half of VmT (at high concentrations of TDB)
    ndrugactT = 0.8,                # Hill coefficient for the TDB portion of the Michaelis-Menten equation
    S = 1.4,                        # Hill coefficient for the B:T portion of the Michaelis-Menten equation
    KdrugactT = 130.08,             # the TDB conc at which T-cell activation was half of VmT [ng/mL]
    # TDB : B-cell kill
    VmB = 0.95,                     # maximum rate of B-cell killing 
    KmTB_kill = 0.75,               # the Tact/B ratio that the rate if half of VmB
    fKmTB_kill = 10.0,              # KmTB_kill is multiplied by this factor (>1) to compensate for the fact the tissues are not well-stirred environments like PB
    nkill = 1.03,                   # Hill coefficient for the Tact/B portion of the Michaelis-Menten equation
    KdrugB = 1.3,                   # the TDB conc at which B-cell killing rate if half of VmB [ng/mL]
    #
    # B-cells
    kBapop = 0.02,                  # Apoptosis rate of B-cells [1/d]
    kBkill = 275.19,                # killing rate of B-cells [1/d]
    fBkill = 0.1,                   # kBkill is multiplied by this factor (<1) to capture the less efficient killing of tissue B-cells
    kBprolif = 0.05,                # proliferation rate of B cells [1/d]
    # 
    # T-cells
    kTprolif = 0.7 ,                # fraction of activated cells that proliferate
    kTaexit = 0.12,                 # trafficking rate of activated CD8+ T-cells from PB to tissue [1/d]
    kTact = 9.83,                   # conversion rate of resting or post-activated CD8+ T cells (and vice versa) [1/d]
    fTadeact = 0.01,                # fraction of activated CD8+ T-cells that deactivate to resting or "exhausted" CD8+ T-cells
    fTap = 3.77,                    # Ratio of the partition coeff for activated CD8+ T-cells to that of resting CD8+ T-cells
    kTaapop = 0.06,                 # apoptosis rate of CD69+ CD8+ T-cells (activated CD8+ T-cells) [1/d]
    fTaprolif = 2.11,               # proliferation rate of activated CD8+ T-cells [1/d]
    fTrapop = 0.2,                  # Ratio of the apoptosis rate of resting CD8+ T-cells to that of CD69+CD8+ T-cells
    kTrexit = 0.05,                 # trafficking rate of resting/ post-activated CD8+ T-cells from PB to tissue [1/d]
    fa0 = 0.,                       # Ratio of the placebo/TDB injection effect on traficking of activated CD8+ T-cells from PB to tissues to that of resting CD8+ T-cells
    fAICD = 1.5,                    # Effect of activation-induced cell death (AICD) on CD69+CD8+ T-cells
    fTa0deact = 0,                  # Conversion rate of post-activated CD8+ T-cells to resting CD8+ T-cells [1/d]
    fTa0apop = 2.02,                # Ratio of the apoptosis rate of post-activated CD8+ T-cells to that of CD69+CD8+ T-cells
    finj = 1.,                      # The placebo injection effect on traficking of resting or "exhausted" CD8+ T-cells from PB to tissues
    #
    # PB Physiology
    Trpbref_perml = 2.0e6,          # baseline conc of resting CD8+ T-cells in PB [cells/mL]
    Bpbref_perml = 0.5e6,           # baseline conc of CD19+CD20+ B-cells in PB [cells/mL]
    Vpb = 5000.,                    # physiological volume of peripheral blood [mL]
    # 
    # Spleen Physiology 
    # Trtiss_perml = 1e9            # Baseline conc. of resting CD8+ T-cells in spleen [cells/ml]
    # Btiss_perml = 9e8             # Baseline conc. of CD19+CD20+ B-cells in spleen [cells/ml]
    Vtissue = 210.,                  # physiological volume of spleen [mL]
    KTrp = 200.0,                   # Partition coeff: ratio of the baseline conc. of resting CD8+ T-cells in spleen to that in PB          
    KBp = 333. ,                    # Partition coeff: ratio of the baseline conc. of B-cells in spleen to that in PB
    Kp = 0.14,                      # Partition coeff: ratio of the drug conc. in spleen to that in PB at any given time
    # 
    # LNs Physiology
    # Trtiss2_perml = 1e9           # Baseline conc. of resting CD8+ T-cells in LNs [cells/ml]
    # Btiss2_perml = 6e8            # Baseline conc. of CD19+CD20+ B-cells in LNs [cells/ml]
    Vtissue2 = 400. ,               # total physiological volume of LNs [mL]
    KTrp2 = 190.,                   # Partition coeff: ratio of the baseline conc. of B-cells in LNs to that in PB
    KBp2 = 190.,                    # Partition coeff: ratio of the baseline conc. of resting CD8+ T-cells in LNs to that in PB
    Kp2 = 0.07,                     # Partition coeff: ratio of the drug conc. in LNs to that in PB at any given time
    # 
    # BM Physiology
    # Trtiss3_perml = 1e8           # Baseline conc. of resting CD8+ T-cells in BM [cells/ml]
    # B1920tiss3_perml = 6e7        # Baseline conc. of CD19+CD20+ B-cells in BM [cells/ml]
    # B19no20tiss3_perml = 1.5e7    # Baseline conc. of CD19+CD20- B-cells in BM [cells/ml]
    # B19tiss3_perml = 7.5e7        # Baseline conc. of CD19+ B cells in BM [cells/ml]
    Vtissue3 = 500.0,               # Physiology volume of BM [mL]
    KTrp3 = 60.0,                   # Partition coeff: ratio of the baseline conc. of resting CD8+ T-cells in BM to that in PB
    KBp3 = 80.0,                    # Partition coeff: ratio of the baseline conc. of CD19+CD20+ B-cells in BM to that in PB
    Kp3 = 0.07,                     # Partition coeff: ratio of the drug conc. in BM to that in PB at any given time
    kBmat_kBapop_ratio = 0.25,      # The ratio of maturation rate of CD19+CD20- B-cells in BM to B-cell apoptosis rate
    B19no20_B1920_ratio = 0.25,     # CD19+CD20- B-cell:CD19+CD20+ B-cell ratio in BM
    # 
    # Tumor Physiology (THIS PART IS TURNED OFF IN MONKEY)
    tumor_on = 1.0,                 # this turns on/ off the tumor compartment 
    # Trtumor_perml                 # Baseline conc. of resting CD8+ T-cells in tumor [cells/ml]
    # Btumor_perml                  # Baseline conc. of CD19+CD20+ B-cells in tumor [cells/ml]
    Vtumor = 49.6,                  # Physiological volume of tumor [ml]
    KTrptumor = 125.,               # Partition coeff: ratio of the baseline conc. of resting CD8+ T-cells in tumor to that in PB
    KBptumor = 6.5E3,               # Partition coeff: ratio of the baseline conc. of B-cells in tumor to that in PB
    Kptumor = 0.05,                 # Partition coeff: ratio of the drug conc. in tumor to that in PB at any given time
    kBtumorprolif = 0.025,          # proliferation rate of B-cells in tumor [1/d]
    # 
    # Cytokine 
    kIL6prod = 14.0,                # production rate of IL6 [pg/d]
    IL6_tiss_contribution = 0.,   # Contribution of tissue-produced IL6 to the systematic IL6 conc
    thalfIL6 = 20.0,                # [minute]
    #
    # Init-related values
    Trpbo_perml = 2.0e6,    # [1/mL]
    Bpbo_perml = 1.0e6,     # [1/mL]
    # 
    # not documented in Supp 2 
    Vm_tdb = 212 ,          # [μg/d]
    Km_tdb = 4.74 ,         # [μg/d]
    act0on = 1.0,
    tissue1on = 1.0,
    tissue2on = 1.0,
    tissue3on = 1.0,
    Bcell_tumor_trafficking_on = 0.,
    BW = 75. ,  # [kg]
    fBexit = 1.0,
    kTgen = 1.0 ,  # [1/d] 
    tinjhalf = 1.0,   # [d]
    kBtiss3exit = 0.03,
    fBtissue3_v1 = 1.0,
    kBapop_cll = 1.0,
    kabs_TDB = 0.26,
    fbio_TDB = 0.7, 
    # mosun PK; note mosun parameters here are already computed out, and outdated parameters were removed 
    # note the clerance of the mosun was based on cyno parameter; this adjustment was based on PK simulation comparison to the clinical data
    CL_TDB = 637.5,       # [mL/day]
    Q_TDB = 1812.75,      # [mL/day]
    V1_TDB = 2758.5,      # [mL]
    V2_TDB = 13010.25,    # [mL]
    infusion_TDB = 0.0,
    # Lonca PK 
    CL_lonca = 0.015e3*DayToHour, 
    V1_lonca = 3.1e3, 
    V2_lonca=3.13e3, 
    Q_lonca= 0.13e3*DayToHour, 
    infusion_lonca = 0.0, 
    # lonca PD
    lonca_kill = 0.05 * 0.15 * DayToHour, 
    lonca_EC50 = 1.76E-5, 
    n_kill_lonca = 0.708,
    k_trans = 0.05, 
    lonca_kill_neg = 0., 
    k_trans_neg = 0.02, 
    lonca_EC50_neg = 1.76E-3
);
