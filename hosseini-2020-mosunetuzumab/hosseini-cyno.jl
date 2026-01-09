# date: 9/26/24
# author: Yuezhe Li 
# purpose of this code: to use the cyno version of the Hosseini model 

using Pkg; Pkg.activate("")

using DifferentialEquations, ModelingToolkit 
using ModelingToolkit: getdefault

# original hosseini model 
include("tdb_homo.jl");
@mtkbuild tdb = TDB_homo();

# update Hosseini model to glofitamab and cyno PK
const bw_cyno = 3.; # [kg];  Frances et al., 2022; https://jpharmsci.org/article/S0022-3549(21)00697-3/fulltext
prob0 = ODEProblem(tdb, [tdb.TCEinjection_effect => 10., tdb.Bpb => 0., tdb.Btiss => 0, tdb.Btiss2 => 0, tdb.B1920tiss3 => 0], (0., 3.), 
                    [
                    # glofit cyno PK 
                    tdb.CL_TDB => 500*bw_cyno, tdb.Q_TDB => 100*bw_cyno, tdb.V1_TDB => 34.7*bw_cyno, tdb.V2_TDB => 40*bw_cyno, tdb.Vm_tdb => 0.,
                    # cyno params
                    tdb.kIL6prod => 4., tdb.Vpb => 380, tdb.Vtissue => 7., tdb.Vtissue2 => 25, tdb.Vtissue3 => 50, 
                    tdb.KTrp => 500, tdb.KTrp2 => 500, tdb.KTrp3 => 50, tdb.KBp => 900, tdb.KBp2 => 600, tdb.KBp3 => 60]);

