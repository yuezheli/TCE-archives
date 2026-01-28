# date: 2/12/2025
# purpose of this code: to host constants in the model 

const nmol_per_mol = 1.0e9 
const NAv = 6.022e23 #  Avogadro's number (1/mole)
const s_per_hr = 3600.0 
const s_per_d = 3600.0*24.0 
const mm3_per_L = 1.0e6
const mm3_per_mL = 1E3
const mm_per_cm = 10
const g_cell_per_l = 1000 #  density of tumor tissue; [g/L]
const Vwell = 1E-4 # Volume of well [L]
const Vc = 4/3 * pi * (4E-6)^3 * 1E3; # tumor cell volume; [L]
const MW_TCE = 15E4
const cyno_WT = 2.6 # assume a 2.6 kg cyno monkey [kg]
const mouse_WT = 0.028 # assume a 28 g male mouse [kg] (Shah& Betts, 2012)
const human_WT = 70 # assume a 70 kg human [kg] 
const hr_per_day = 24
const ng_per_ug = 1E3
const ug_per_mg = 1E3
const um_per_mm = 1E3
const mg_per_g = 1E3
const ng_per_g = 1E9
const mL_per_L = 1E3
const cell_vol = 4E-12 # [unit = u"L", description = "volume of a cell"] (source: https://pubmed.ncbi.nlm.nih.gov/15736406/) 

const MW_IFNg = 17E3   # https://pubmed.ncbi.nlm.nih.gov/3109913/
const MW_IL10 = 18E3   # https://www.rndsystems.com/resources/articles/interleukin-10-il-10-family
const MW_TNFa = 17.4E3 # https://www.novusbio.com/products/tnf-alpha-antibody-52b83_nb600-1422
const MW_IL6 = 23E3 # https://pubmed.ncbi.nlm.nih.gov/2033081/

# TAA membrane count per cell (QiFI number) (obtained from Table 1.4)
TAA_AspC1 = 17635
TAA_HPAC = 66322
TAA_RMUGS = 1863 # BI communication, 2/18/2025
