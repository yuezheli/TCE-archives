# date: 9/18/24
# author: Yuezhe Li 
# purpose of this code: in vitro model isolated from Hosseini model 
# this was copied from `glofit-in-vitro.jl` from Project ADC0501; the original script was fully QCed; 

include("param.jl");

#=
p_mosun = ComponentArray(
    kBtumorprolif = p_homo_3.kBtumorprolif, KmBT_act = p_homo_3.KmBT_act, S = p_homo_3.S, KdrugactT = p_homo_3.KdrugactT, ndrugactT = p_homo_3.ndrugactT, VmT = p_homo_3.VmT, 
    KmTB_kill = p_homo_3.KmTB_kill, nkill = p_homo_3.nkill, KdrugB = p_homo_3.KdrugB, VmB = p_homo_3.VmB, kTact = p_homo_3.kTact, fTadeact = p_homo_3.fTadeact, fTa0deact = p_homo_3.fTa0deact); 


p_glofit = deepcopy(p_mosun);
p_glofit.KmBT_act = 0.98;   # p_mosun.KmBT_act = 0.72
p_glofit.KdrugactT = 0.41;  # p_mosun.KdrugactT = 130.08
p_glofit.KdrugB = 1E-4;     # p_mosun.KdrugB = 1.3

=#

function bs_in_vitro_hosseini!(du, u, p, t)
    @unpack Btumor, bs_ugperml, restT, actT, act0T = u
    @unpack kBtumorprolif, KmBT_act, S, KdrugactT, ndrugactT, VmT, KmTB_kill, nkill, KdrugB, VmB, kTact, fTadeact, fTa0deact = p

    BTrratio = max(Btumor,1)/max(restT+act0T,1)
    TaBratio = max(actT, 1.)/max(Btumor,1)

    drugTact = VmT*(BTrratio^S/(KmBT_act^S+BTrratio^S))*((bs_ugperml*1000)^ndrugactT/(KdrugactT^ndrugactT+(bs_ugperml*1000)^ndrugactT))
    drugBkill = VmB*max(0, TaBratio)^nkill/(max(0, KmTB_kill)^nkill+max(0, TaBratio)^nkill)*((bs_ugperml*1000)/(KdrugB+(bs_ugperml*1000)))


    du.Btumor = kBtumorprolif*Btumor - drugBkill*Btumor
        
    du.restT = ( -kTact*drugTact*restT + fTa0deact*act0T  )
    du.actT = kTact*drugTact*restT - kTact*fTadeact*actT
    du.act0T = ( -fTa0deact*act0T + kTact*fTadeact*actT )
end
