#pragma once
#include <iostream>
#include <string>
#include <unordered_map>
#include "bsim4v82.hpp"


namespace bsim4{

// Parameter map for BSIM4
const std::unordered_map<std::string, int> bsim4ParamMap = {
{"l", BSIM4_L}, // Length
{"w", BSIM4_W}, // Width
{"m", BSIM4_M}, // Separate Parallel multiplier
{"nf", BSIM4_NF}, // Number of fingers
{"sa", BSIM4_SA}, // distance between  OD edge to poly of one side 
{"sb", BSIM4_SB}, // distance between  OD edge to poly of the other side
{"sd", BSIM4_SD}, // distance between neighbour fingers
{"sca", BSIM4_SCA}, // Integral of the first distribution function for scattered well dopant
{"scb", BSIM4_SCB}, // Integral of the second distribution function for scattered well dopant
{"scc", BSIM4_SCC}, // Integral of the third distribution function for scattered well dopant
{"sc", BSIM4_SC}, // Distance to a single well edge 
{"min", BSIM4_MIN}, // Minimize either D or S
{"ad", BSIM4_AD}, // Drain area
{"as", BSIM4_AS}, // Source area
{"pd", BSIM4_PD}, // Drain perimeter
{"ps", BSIM4_PS}, // Source perimeter
{"nrd", BSIM4_NRD}, // Number of squares in drain
{"nrs", BSIM4_NRS}, // Number of squares in source
{"off", BSIM4_OFF}, // Device is initially off
{"rbdb", BSIM4_RBDB}, // Body resistance
{"rbsb", BSIM4_RBSB}, // Body resistance
{"rbpb", BSIM4_RBPB}, // Body resistance
{"rbps", BSIM4_RBPS}, // Body resistance
{"rbpd", BSIM4_RBPD}, // Body resistance
{"delvto", BSIM4_DELVTO}, // Zero bias threshold voltage variation
{"delvt0", BSIM4_DELVTO}, // Zero bias threshold voltage variation
// {"mulu0", BSIM4_MULU0}, // Low field mobility multiplier
{"xgw", BSIM4_XGW}, // Distance from gate contact center to device edge
{"ngcon", BSIM4_NGCON}, // Number of gate contacts
// {"wnflag", BSIM4_WNFLAG}, // W/NF device flag for bin selection

{"trnqsmod", BSIM4_TRNQSMOD}, // Transient NQS model selector
{"acnqsmod", BSIM4_ACNQSMOD}, // AC NQS model selector
{"rbodymod", BSIM4_RBODYMOD}, // Distributed body R model selector
{"rgatemod", BSIM4_RGATEMOD}, // Gate resistance model selector
{"geomod", BSIM4_GEOMOD}, // Geometry dependent parasitics model selector
{"rgeomod", BSIM4_RGEOMOD}, // S/D resistance and contact model selector
{"ic", BSIM4_IC}, // Vector of DS,GS,BS initial voltages
{"gmbs", BSIM4_GMBS}, // Gmb
{"gm", BSIM4_GM}, // Gm
{"gds", BSIM4_GDS}, // Gds
{"vdsat", BSIM4_VDSAT}, // Vdsat
{"vth", BSIM4_VON}, // Vth
{"id", BSIM4_CD}, // Ids
{"ibd", BSIM4_CBD}, // Ibd
{"ibs", BSIM4_CBS}, // Ibs
{"gbd", BSIM4_GBD}, // gbd
{"gbs", BSIM4_GBS}, // gbs
{"isub", BSIM4_CSUB}, // Isub
{"igidl", BSIM4_IGIDL}, // Igidl
{"igisl", BSIM4_IGISL}, // Igisl
{"igs", BSIM4_IGS}, // Igs
{"igd", BSIM4_IGD}, // Igd
{"igb", BSIM4_IGB}, // Igb
{"igcs", BSIM4_IGCS}, // Igcs
{"igcd", BSIM4_IGCD}, // Igcd
{"vbs", BSIM4_VBS}, // Vbs
{"vgs", BSIM4_VGS}, // Vgs
{"vds", BSIM4_VDS}, // Vds
{"cgg", BSIM4_CGGB}, // Cggb
{"cgs", BSIM4_CGSB}, // Cgsb
{"cgd", BSIM4_CGDB}, // Cgdb
{"cbg", BSIM4_CBGB}, // Cbgb
{"cbd", BSIM4_CBDB}, // Cbdb
{"cbs", BSIM4_CBSB}, // Cbsb
{"cdg", BSIM4_CDGB}, // Cdgb
{"cdd", BSIM4_CDDB}, // Cddb
{"cds", BSIM4_CDSB}, // Cdsb
{"csg", BSIM4_CSGB}, // Csgb
{"csd", BSIM4_CSDB}, // Csdb
{"css", BSIM4_CSSB}, // Cssb
{"cgb", BSIM4_CGBB}, // Cgbb
{"cdb", BSIM4_CDBB}, // Cdbb
{"csb", BSIM4_CSBB}, // Csbb
{"cbb", BSIM4_CBBB}, // Cbbb
{"capbd", BSIM4_CAPBD}, // Capbd
{"capbs", BSIM4_CAPBS}, // Capbs
{"qg", BSIM4_QG}, // Qgate
{"qb", BSIM4_QB}, // Qbulk
{"qd", BSIM4_QD}, // Qdrain
{"qs", BSIM4_QS}, // Qsource
{"qinv", BSIM4_QINV}, // Qinversion
{"qdef", BSIM4_QDEF}, // Qdef
{"gcrg", BSIM4_GCRG}, // Gcrg
{"gtau", BSIM4_GTAU}, // Gtau
// {"vgsteff", BSIM4_VGSTEFF}, // Vgsteff
// {"vdseff", BSIM4_VDSEFF}, // Vdseff
// {"cgso", BSIM4_CGSO}, // Cgso
// {"cgdo", BSIM4_CGDO}, // Cgdo
// {"cgbo", BSIM4_CGBO}, // Cgbo
// {"weff", BSIM4_WEFF}, // Weff
// {"leff", BSIM4_LEFF}, // Leff
};

// Model parameter map for BSIM4
const std::unordered_map<std::string, int> bsim4ModelParamMap = {
{"cvchargemod", BSIM4_MOD_CVCHARGEMOD}, // Capacitance Charge model selector
{"capmod", BSIM4_MOD_CAPMOD}, // Capacitance model selector
{"diomod", BSIM4_MOD_DIOMOD}, // Diode IV model selector
{"rdsmod", BSIM4_MOD_RDSMOD}, // Bias-dependent S/D resistance model selector
{"trnqsmod", BSIM4_MOD_TRNQSMOD}, // Transient NQS model selector
{"acnqsmod", BSIM4_MOD_ACNQSMOD}, // AC NQS model selector
{"mobmod", BSIM4_MOD_MOBMOD}, // Mobility model selector
{"rbodymod", BSIM4_MOD_RBODYMOD}, // Distributed body R model selector
{"rgatemod", BSIM4_MOD_RGATEMOD}, // Gate R model selector
{"permod", BSIM4_MOD_PERMOD}, // Pd and Ps model selector
{"geomod", BSIM4_MOD_GEOMOD}, // Geometry dependent parasitics model selector
{"rgeomod", BSIM4_MOD_GEOMOD}, // S/D resistance and contact model selector (BSIM4_MOD_RGEOMOD)
{"fnoimod", BSIM4_MOD_FNOIMOD}, // Flicker noise model selector
{"tnoimod", BSIM4_MOD_TNOIMOD}, // Thermal noise model selector
{"mtrlmod", BSIM4_MOD_MTRLMOD}, // parameter for non-silicon substrate or metal gate selector
{"mtrlcompatmod", BSIM4_MOD_MTRLCOMPATMOD}, // New Material Mod backward compatibility selector
{"igcmod", BSIM4_MOD_IGCMOD}, // Gate-to-channel Ig model selector
{"igbmod", BSIM4_MOD_IGBMOD}, // Gate-to-body Ig model selector
{"tempmod", BSIM4_MOD_TEMPMOD}, // Temperature model selector
{"gidlmod", BSIM4_MOD_GIDLMOD}, // parameter for GIDL selector
{"paramchk", BSIM4_MOD_PARAMCHK}, // Model parameter checking selector
{"binunit", BSIM4_MOD_BINUNIT}, // Bin  unit  selector
{"version", BSIM4_MOD_VERSION}, // parameter for model version
{"eot", BSIM4_MOD_EOT}, // Equivalent gate oxide thickness in meters
{"vddeot", BSIM4_MOD_VDDEOT}, // Voltage for extraction of Equivalent gate oxide thickness
{"tempeot", BSIM4_MOD_TEMPEOT}, //  Temperature for extraction of EOT
{"leffeot", BSIM4_MOD_LEFFEOT}, //  Effective length for extraction of EOT
{"weffeot", BSIM4_MOD_WEFFEOT}, // Effective width for extraction of EOT
{"ados", BSIM4_MOD_ADOS}, // Charge centroid parameter
{"bdos", BSIM4_MOD_BDOS}, // Charge centroid parameter
{"toxe", BSIM4_MOD_TOXE}, // Electrical gate oxide thickness in meters
{"toxp", BSIM4_MOD_TOXP}, // Physical gate oxide thickness in meters
{"toxm", BSIM4_MOD_TOXM}, // Gate oxide thickness at which parameters are extracted
{"toxref", BSIM4_MOD_TOXREF}, // Target tox value
{"dtox", BSIM4_MOD_DTOX}, // Defined as (toxe - toxp) 
{"epsrox", BSIM4_MOD_EPSROX}, // Dielectric constant of the gate oxide relative to vacuum
{"cdsc", BSIM4_MOD_CDSC}, // Drain/Source and channel coupling capacitance
{"cdscb", BSIM4_MOD_CDSCB}, // Body-bias dependence of cdsc
{"cdscd", BSIM4_MOD_CDSCD}, // Drain-bias dependence of cdsc
{"cit", BSIM4_MOD_CIT}, // Interface state capacitance
{"nfactor", BSIM4_MOD_NFACTOR}, // Subthreshold swing Coefficient
{"xj", BSIM4_MOD_XJ}, // Junction depth in meters
{"vsat", BSIM4_MOD_VSAT}, // Saturation velocity at tnom
{"at", BSIM4_MOD_AT}, // Temperature coefficient of vsat
{"a0", BSIM4_MOD_A0}, // Non-uniform depletion width effect coefficient.
{"ags", BSIM4_MOD_AGS}, // Gate bias  coefficient of Abulk.
{"a1", BSIM4_MOD_A1}, // Non-saturation effect coefficient
{"a2", BSIM4_MOD_A2}, // Non-saturation effect coefficient
{"keta", BSIM4_MOD_KETA}, // Body-bias coefficient of non-uniform depletion width effect.
{"phig", BSIM4_MOD_PHIG}, // Work function of gate
{"epsrgate", BSIM4_MOD_EPSRGATE}, // Dielectric constant of gate relative to vacuum
{"easub", BSIM4_MOD_EASUB}, // Electron affinity of substrate
{"epsrsub", BSIM4_MOD_EPSRSUB}, // Dielectric constant of substrate relative to vacuum
{"ni0sub", BSIM4_MOD_NI0SUB}, // Intrinsic carrier concentration of substrate at 300.15K
{"bg0sub", BSIM4_MOD_BG0SUB}, // Band-gap of substrate at T=0K
{"tbgasub", BSIM4_MOD_TBGASUB}, // First parameter of band-gap change due to temperature
{"tbgbsub", BSIM4_MOD_TBGBSUB}, // Second parameter of band-gap change due to temperature
{"nsub", BSIM4_MOD_NSUB}, // Substrate doping concentration
{"ndep", BSIM4_MOD_NDEP}, // Channel doping concentration at the depletion edge
{"nsd", BSIM4_MOD_NSD}, // S/D doping concentration
{"phin", BSIM4_MOD_PHIN}, // Adjusting parameter for surface potential due to non-uniform vertical doping
{"ngate", BSIM4_MOD_NGATE}, // Poly-gate doping concentration
{"gamma1", BSIM4_MOD_GAMMA1}, // Vth body coefficient
{"gamma2", BSIM4_MOD_GAMMA2}, // Vth body coefficient
{"vbx", BSIM4_MOD_VBX}, // Vth transition body Voltage
{"vbm", BSIM4_MOD_VBM}, // Maximum body voltage

{"xt", BSIM4_MOD_XT}, // Doping depth
{"k1", BSIM4_MOD_K1}, // Bulk effect coefficient 1
{"kt1", BSIM4_MOD_KT1}, // Temperature coefficient of Vth
{"kt1l", BSIM4_MOD_KT1L}, // Temperature coefficient of Vth
{"kt2", BSIM4_MOD_KT2}, // Body-coefficient of kt1
{"k2", BSIM4_MOD_K2}, // Bulk effect coefficient 2
{"k3", BSIM4_MOD_K3}, // Narrow width effect coefficient
{"k3b", BSIM4_MOD_K3B}, // Body effect coefficient of k3
{"w0", BSIM4_MOD_W0}, // Narrow width effect parameter
{"dvtp0", BSIM4_MOD_DVTP0}, // First parameter for Vth shift due to pocket
{"dvtp1", BSIM4_MOD_DVTP1}, // Second parameter for Vth shift due to pocket
{"dvtp2", BSIM4_MOD_DVTP2}, // 3rd parameter for Vth shift due to pocket
{"dvtp3", BSIM4_MOD_DVTP3}, // 4th parameter for Vth shift due to pocket
{"dvtp4", BSIM4_MOD_DVTP4}, // 5th parameter for Vth shift due to pocket
{"dvtp5", BSIM4_MOD_DVTP5}, // 6th parameter for Vth shift due to pocket
{"lpe0", BSIM4_MOD_LPE0}, // Equivalent length of pocket region at zero bias
{"lpeb", BSIM4_MOD_LPEB}, // Equivalent length of pocket region accounting for body bias
{"dvt0", BSIM4_MOD_DVT0}, // Short channel effect coeff. 0
{"dvt1", BSIM4_MOD_DVT1}, // Short channel effect coeff. 1
{"dvt2", BSIM4_MOD_DVT2}, // Short channel effect coeff. 2
{"dvt0w", BSIM4_MOD_DVT0W}, // Narrow Width coeff. 0
{"dvt1w", BSIM4_MOD_DVT1W}, // Narrow Width effect coeff. 1
{"dvt2w", BSIM4_MOD_DVT2W}, // Narrow Width effect coeff. 2
{"drout", BSIM4_MOD_DROUT}, // DIBL coefficient of output resistance
{"dsub", BSIM4_MOD_DSUB}, // DIBL coefficient in the subthreshold region
{"vth0", BSIM4_MOD_VTH0}, // Threshold voltage
{"vtho", BSIM4_MOD_VTH0}, // Threshold voltage
{"ua", BSIM4_MOD_UA}, // Linear gate dependence of mobility
{"ua1", BSIM4_MOD_UA1}, // Temperature coefficient of ua
{"ub", BSIM4_MOD_UB}, // Quadratic gate dependence of mobility
{"ub1", BSIM4_MOD_UB1}, // Temperature coefficient of ub
{"uc", BSIM4_MOD_UC}, // Body-bias dependence of mobility
{"uc1", BSIM4_MOD_UC1}, // Temperature coefficient of uc
{"ud", BSIM4_MOD_UD}, // Coulomb scattering factor of mobility
{"ud1", BSIM4_MOD_UD1}, // Temperature coefficient of ud
{"up", BSIM4_MOD_UP}, // Channel length linear factor of mobility
{"lp", BSIM4_MOD_LP}, // Channel length exponential factor of mobility
{"u0", BSIM4_MOD_U0}, // Low-field mobility at Tnom
{"eu", BSIM4_MOD_EU}, // Mobility exponent
{"ucs", BSIM4_MOD_UCS}, // Colombic scattering exponent
{"ute", BSIM4_MOD_UTE}, // Temperature coefficient of mobility
{"ucste", BSIM4_MOD_UCSTE}, // Temperature coefficient of colombic mobility
{"voff", BSIM4_MOD_VOFF}, // Threshold voltage offset
{"minv", BSIM4_MOD_MINV}, // Fitting parameter for moderate inversion in Vgsteff
{"minvcv", BSIM4_MOD_MINVCV}, // Fitting parameter for moderate inversion in Vgsteffcv
{"voffl", BSIM4_MOD_VOFFL}, // Length dependence parameter for Vth offset
{"voffcvl", BSIM4_MOD_VOFFCVL}, // Length dependence parameter for Vth offset in CV
{"tnom", BSIM4_MOD_TNOM}, // Parameter measurement temperature
{"cgso", BSIM4_MOD_CGSO}, // Gate-source overlap capacitance per width
{"cgdo", BSIM4_MOD_CGDO}, // Gate-drain overlap capacitance per width
{"cgbo", BSIM4_MOD_CGBO}, // Gate-bulk overlap capacitance per length
{"xpart", BSIM4_MOD_XPART}, // Channel charge partitioning
{"delta", BSIM4_MOD_DELTA}, // Effective Vds parameter
{"rsh", BSIM4_MOD_RSH}, // Source-drain sheet resistance
{"rdsw", BSIM4_MOD_RDSW}, // Source-drain resistance per width
{"rdswmin", BSIM4_MOD_RDSWMIN}, // Source-drain resistance per width at high Vg
{"rsw", BSIM4_MOD_RSW}, // Source resistance per width
{"rdw", BSIM4_MOD_RDW}, // Drain resistance per width
{"rdwmin", BSIM4_MOD_RDWMIN}, // Drain resistance per width at high Vg
{"rswmin", BSIM4_MOD_RSWMIN}, // Source resistance per width at high Vg

{"prwg", BSIM4_MOD_PRWG}, // Gate-bias effect on parasitic resistance 
{"prwb", BSIM4_MOD_PRWB}, // Body-effect on parasitic resistance 

{"prt", BSIM4_MOD_PRT}, // Temperature coefficient of parasitic resistance 
{"eta0", BSIM4_MOD_ETA0}, // Subthreshold region DIBL coefficient
{"etab", BSIM4_MOD_ETAB}, // Subthreshold region DIBL coefficient
{"pclm", BSIM4_MOD_PCLM}, // Channel length modulation Coefficient
{"pdiblc1", BSIM4_MOD_PDIBL1}, // Drain-induced barrier lowering coefficient
{"pdiblc2", BSIM4_MOD_PDIBL2}, // Drain-induced barrier lowering coefficient
{"pdiblcb", BSIM4_MOD_PDIBLB}, // Body-effect on drain-induced barrier lowering
{"fprout", BSIM4_MOD_FPROUT}, // Rout degradation coefficient for pocket devices
{"pdits", BSIM4_MOD_PDITS}, // Coefficient for drain-induced Vth shifts
{"pditsl", BSIM4_MOD_PDITSL}, // Length dependence of drain-induced Vth shifts
{"pditsd", BSIM4_MOD_PDITSD}, // Vds dependence of drain-induced Vth shifts
{"pscbe1", BSIM4_MOD_PSCBE1}, // Substrate current body-effect coefficient
{"pscbe2", BSIM4_MOD_PSCBE2}, // Substrate current body-effect coefficient
{"pvag", BSIM4_MOD_PVAG}, // Gate dependence of output resistance parameter

{"jss", BSIM4_MOD_JSS}, // Bottom source junction reverse saturation current density
{"jsws", BSIM4_MOD_JSWS}, // Isolation edge sidewall source junction reverse saturation current density
{"jswgs", BSIM4_MOD_JSWGS}, // Gate edge source junction reverse saturation current density
{"pbs", BSIM4_MOD_PBS}, // Source junction built-in potential
{"njs", BSIM4_MOD_NJS}, // Source junction emission coefficient
{"xtis", BSIM4_MOD_XTIS}, // Source junction current temperature exponent
{"mjs", BSIM4_MOD_MJS}, // Source bottom junction capacitance grading coefficient
{"pbsws", BSIM4_MOD_PBSWS}, // Source sidewall junction capacitance built in potential
{"mjsws", BSIM4_MOD_MJSWS}, // Source sidewall junction capacitance grading coefficient
{"pbswgs", BSIM4_MOD_PBSWGS}, // Source (gate side) sidewall junction capacitance built in potential
{"mjswgs", BSIM4_MOD_MJSWGS}, // Source (gate side) sidewall junction capacitance grading coefficient
{"cjs", BSIM4_MOD_CJS}, // Source bottom junction capacitance per unit area
{"cjsws", BSIM4_MOD_CJSWS}, // Source sidewall junction capacitance per unit periphery
{"cjswgs", BSIM4_MOD_CJSWGS}, // Source (gate side) sidewall junction capacitance per unit width

{"jsd", BSIM4_MOD_JSD}, // Bottom drain junction reverse saturation current density
{"jswd", BSIM4_MOD_JSWD}, // Isolation edge sidewall drain junction reverse saturation current density
{"jswgd", BSIM4_MOD_JSWGD}, // Gate edge drain junction reverse saturation current density
{"pbd", BSIM4_MOD_PBD}, // Drain junction built-in potential
{"njd", BSIM4_MOD_NJD}, // Drain junction emission coefficient
{"xtid", BSIM4_MOD_XTID}, // Drainjunction current temperature exponent
{"mjd", BSIM4_MOD_MJD}, // Drain bottom junction capacitance grading coefficient
{"pbswd", BSIM4_MOD_PBSWD}, // Drain sidewall junction capacitance built in potential
{"mjswd", BSIM4_MOD_MJSWD}, // Drain sidewall junction capacitance grading coefficient
{"pbswgd", BSIM4_MOD_PBSWGD}, // Drain (gate side) sidewall junction capacitance built in potential
{"mjswgd", BSIM4_MOD_MJSWGD}, // Drain (gate side) sidewall junction capacitance grading coefficient
{"cjd", BSIM4_MOD_CJD}, // Drain bottom junction capacitance per unit area
{"cjswd", BSIM4_MOD_CJSWD}, // Drain sidewall junction capacitance per unit periphery
{"cjswgd", BSIM4_MOD_CJSWGD}, // Drain (gate side) sidewall junction capacitance per unit width

{"vfbcv", BSIM4_MOD_VFBCV}, // Flat Band Voltage parameter for capmod=0 only
{"vfb", BSIM4_MOD_VFB}, // Flat Band Voltage
{"tpb", BSIM4_MOD_TPB}, // Temperature coefficient of pb
{"tcj", BSIM4_MOD_TCJ}, // Temperature coefficient of cj
{"tpbsw", BSIM4_MOD_TPBSW}, // Temperature coefficient of pbsw
{"tcjsw", BSIM4_MOD_TCJSW}, // Temperature coefficient of cjsw
{"tpbswg", BSIM4_MOD_TPBSWG}, // Temperature coefficient of pbswg
{"tcjswg", BSIM4_MOD_TCJSWG}, // Temperature coefficient of cjswg
{"acde", BSIM4_MOD_ACDE}, // Exponential coefficient for finite charge thickness
{"moin", BSIM4_MOD_MOIN}, // Coefficient for gate-bias dependent surface potential
{"noff", BSIM4_MOD_NOFF}, // C-V turn-on/off parameter
{"voffcv", BSIM4_MOD_VOFFCV}, // C-V lateral-shift parameter
{"dmcg", BSIM4_MOD_DMCG}, // Distance of Mid-Contact to Gate edge
{"dmci", BSIM4_MOD_DMCI}, // Distance of Mid-Contact to Isolation
{"dmdg", BSIM4_MOD_DMDG}, // Distance of Mid-Diffusion to Gate edge
{"dmcgt", BSIM4_MOD_DMCGT}, // Distance of Mid-Contact to Gate edge in Test structures
{"xgw", BSIM4_MOD_XGW}, // Distance from gate contact center to device edge
{"xgl", BSIM4_MOD_XGL}, // Variation in Ldrawn
{"rshg", BSIM4_MOD_RSHG}, // Gate sheet resistance
{"ngcon", BSIM4_MOD_NGCON}, // Number of gate contacts
{"xrcrg1", BSIM4_MOD_XRCRG1}, // First fitting parameter the bias-dependent Rg
{"xrcrg2", BSIM4_MOD_XRCRG2}, // Second fitting parameter the bias-dependent Rg
{"lambda", BSIM4_MOD_LAMBDA}, //  Velocity overshoot parameter
{"vtl", BSIM4_MOD_VTL}, //  thermal velocity
{"lc", BSIM4_MOD_LC}, //  back scattering parameter
{"xn", BSIM4_MOD_XN}, //  back scattering parameter
{"vfbsdoff", BSIM4_MOD_VFBSDOFF}, // S/D flatband voltage offset
{"tvfbsdoff", BSIM4_MOD_TVFBSDOFF}, // Temperature parameter for vfbsdoff
{"tvoff", BSIM4_MOD_TVOFF}, // Temperature parameter for voff
{"tnfactor", BSIM4_MOD_TNFACTOR}, // Temperature parameter for nfactor
{"teta0", BSIM4_MOD_TETA0}, // Temperature parameter for eta0
{"tvoffcv", BSIM4_MOD_TVOFFCV}, // Temperature parameter for tvoffcv

{"lintnoi", BSIM4_MOD_LINTNOI}, // lint offset for noise calculation
{"lint", BSIM4_MOD_LINT}, // Length reduction parameter
{"ll", BSIM4_MOD_LL}, // Length reduction parameter
{"llc", BSIM4_MOD_LLC}, // Length reduction parameter for CV
{"lln", BSIM4_MOD_LLN}, // Length reduction parameter
{"lw", BSIM4_MOD_LW}, // Length reduction parameter
{"lwc", BSIM4_MOD_LWC}, // Length reduction parameter for CV
{"lwn", BSIM4_MOD_LWN}, // Length reduction parameter
{"lwl", BSIM4_MOD_LWL}, // Length reduction parameter
{"lwlc", BSIM4_MOD_LWLC}, // Length reduction parameter for CV
{"lmin", BSIM4_MOD_LMIN}, // Minimum length for the model
{"lmax", BSIM4_MOD_LMAX}, // Maximum length for the model

{"wr", BSIM4_MOD_WR}, // Width dependence of rds
{"wint", BSIM4_MOD_WINT}, // Width reduction parameter
{"dwg", BSIM4_MOD_DWG}, // Width reduction parameter
{"dwb", BSIM4_MOD_DWB}, // Width reduction parameter

{"wl", BSIM4_MOD_WL}, // Width reduction parameter
{"wlc", BSIM4_MOD_WLC}, // Width reduction parameter for CV
{"wln", BSIM4_MOD_WLN}, // Width reduction parameter
{"ww", BSIM4_MOD_WW}, // Width reduction parameter
{"wwc", BSIM4_MOD_WWC}, // Width reduction parameter for CV
{"wwn", BSIM4_MOD_WWN}, // Width reduction parameter
{"wwl", BSIM4_MOD_WWL}, // Width reduction parameter
{"wwlc", BSIM4_MOD_WWLC}, // Width reduction parameter for CV
{"wmin", BSIM4_MOD_WMIN}, // Minimum width for the model
{"wmax", BSIM4_MOD_WMAX}, // Maximum width for the model

{"b0", BSIM4_MOD_B0}, // Abulk narrow width parameter
{"b1", BSIM4_MOD_B1}, // Abulk narrow width parameter

{"cgsl", BSIM4_MOD_CGSL}, // New C-V model parameter
{"cgdl", BSIM4_MOD_CGDL}, // New C-V model parameter
{"ckappas", BSIM4_MOD_CKAPPAS}, // S/G overlap C-V parameter 
{"ckappad", BSIM4_MOD_CKAPPAD}, // D/G overlap C-V parameter
{"cf", BSIM4_MOD_CF}, // Fringe capacitance parameter
{"clc", BSIM4_MOD_CLC}, // Vdsat parameter for C-V model
{"cle", BSIM4_MOD_CLE}, // Vdsat parameter for C-V model
{"dwc", BSIM4_MOD_DWC}, // Delta W for C-V model
{"dlc", BSIM4_MOD_DLC}, // Delta L for C-V model
{"xw", BSIM4_MOD_XW}, // W offset for channel width due to mask/etch effect
{"xl", BSIM4_MOD_XL}, // L offset for channel length due to mask/etch effect
{"dlcig", BSIM4_MOD_DLCIG}, // Delta L for Ig model
{"dlcigd", BSIM4_MOD_DLCIGD}, // Delta L for Ig model drain side
{"dwj", BSIM4_MOD_DWJ}, // Delta W for S/D junctions

{"alpha0", BSIM4_MOD_ALPHA0}, // substrate current model parameter
{"alpha1", BSIM4_MOD_ALPHA1}, // substrate current model parameter
{"beta0", BSIM4_MOD_BETA0}, // substrate current model parameter

{"agidl", BSIM4_MOD_AGIDL}, // Pre-exponential constant for GIDL
{"bgidl", BSIM4_MOD_BGIDL}, // Exponential constant for GIDL
{"cgidl", BSIM4_MOD_CGIDL}, // Parameter for body-bias dependence of GIDL
{"rgidl", BSIM4_MOD_RGIDL}, // GIDL vg parameter
{"kgidl", BSIM4_MOD_KGIDL}, // GIDL vb parameter
{"fgidl", BSIM4_MOD_FGIDL}, // GIDL vb parameter
{"egidl", BSIM4_MOD_EGIDL}, // Fitting parameter for Bandbending
{"agisl", BSIM4_MOD_AGISL}, // Pre-exponential constant for GISL
{"bgisl", BSIM4_MOD_BGISL}, // Exponential constant for GISL
{"cgisl", BSIM4_MOD_CGISL}, // Parameter for body-bias dependence of GISL
{"rgisl", BSIM4_MOD_RGISL}, // GISL vg parameter
{"kgisl", BSIM4_MOD_KGISL}, // GISL vb parameter
{"fgisl", BSIM4_MOD_FGISL}, // GISL vb parameter
{"egisl", BSIM4_MOD_EGISL}, // Fitting parameter for Bandbending
{"aigc", BSIM4_MOD_AIGC}, // Parameter for Igc
{"bigc", BSIM4_MOD_BIGC}, // Parameter for Igc
{"cigc", BSIM4_MOD_CIGC}, // Parameter for Igc
{"aigsd", BSIM4_MOD_AIGSD}, // Parameter for Igs,d
{"bigsd", BSIM4_MOD_BIGSD}, // Parameter for Igs,d
{"cigsd", BSIM4_MOD_CIGSD}, // Parameter for Igs,d
{"aigs", BSIM4_MOD_AIGS}, // Parameter for Igs
{"bigs", BSIM4_MOD_BIGS}, // Parameter for Igs
{"cigs", BSIM4_MOD_CIGS}, // Parameter for Igs
{"aigd", BSIM4_MOD_AIGD}, // Parameter for Igd
{"bigd", BSIM4_MOD_BIGD}, // Parameter for Igd
{"cigd", BSIM4_MOD_CIGD}, // Parameter for Igd
{"aigbacc", BSIM4_MOD_AIGBACC}, // Parameter for Igb
{"bigbacc", BSIM4_MOD_BIGBACC}, // Parameter for Igb
{"cigbacc", BSIM4_MOD_CIGBACC}, // Parameter for Igb
{"aigbinv", BSIM4_MOD_AIGBINV}, // Parameter for Igb
{"bigbinv", BSIM4_MOD_BIGBINV}, // Parameter for Igb
{"cigbinv", BSIM4_MOD_CIGBINV}, // Parameter for Igb
{"nigc", BSIM4_MOD_NIGC}, // Parameter for Igc slope
{"nigbinv", BSIM4_MOD_NIGBINV}, // Parameter for Igbinv slope
{"nigbacc", BSIM4_MOD_NIGBACC}, // Parameter for Igbacc slope
{"ntox", BSIM4_MOD_NTOX}, // Exponent for Tox ratio
{"eigbinv", BSIM4_MOD_EIGBINV}, // Parameter for the Si bandgap for Igbinv
{"pigcd", BSIM4_MOD_PIGCD}, // Parameter for Igc partition
{"poxedge", BSIM4_MOD_POXEDGE}, // Factor for the gate edge Tox

{"ijthdfwd", BSIM4_MOD_IJTHDFWD}, // Forward drain diode forward limiting current
{"ijthsfwd", BSIM4_MOD_IJTHSFWD}, // Forward source diode forward limiting current
{"ijthdrev", BSIM4_MOD_IJTHDREV}, // Reverse drain diode forward limiting current
{"ijthsrev", BSIM4_MOD_IJTHSREV}, // Reverse source diode forward limiting current
{"xjbvd", BSIM4_MOD_XJBVD}, // Fitting parameter for drain diode breakdown current
{"xjbvs", BSIM4_MOD_XJBVS}, // Fitting parameter for source diode breakdown current
{"bvd", BSIM4_MOD_BVD}, // Drain diode breakdown voltage
{"bvs", BSIM4_MOD_BVS}, // Source diode breakdown voltage

{"jtss", BSIM4_MOD_JTSS}, // Source bottom trap-assisted saturation current density
{"jtsd", BSIM4_MOD_JTSD}, // Drain bottom trap-assisted saturation current density
{"jtssws", BSIM4_MOD_JTSSWS}, // Source STI sidewall trap-assisted saturation current density
{"jtsswd", BSIM4_MOD_JTSSWD}, // Drain STI sidewall trap-assisted saturation current density
{"jtsswgs", BSIM4_MOD_JTSSWGS}, // Source gate-edge sidewall trap-assisted saturation current density
{"jtsswgd", BSIM4_MOD_JTSSWGD}, // Drain gate-edge sidewall trap-assisted saturation current density
{"jtweff", BSIM4_MOD_JTWEFF}, // TAT current width dependence
{"njts", BSIM4_MOD_NJTS}, // Non-ideality factor for bottom junction
{"njtssw", BSIM4_MOD_NJTSSW}, // Non-ideality factor for STI sidewall junction
{"njtsswg", BSIM4_MOD_NJTSSWG}, // Non-ideality factor for gate-edge sidewall junction
{"njtsd", BSIM4_MOD_NJTSD}, // Non-ideality factor for bottom junction drain side
{"njtsswd", BSIM4_MOD_NJTSSWD}, // Non-ideality factor for STI sidewall junction drain side
{"njtsswgd", BSIM4_MOD_NJTSSWGD}, // Non-ideality factor for gate-edge sidewall junction drain side
{"xtss", BSIM4_MOD_XTSS}, // Power dependence of JTSS on temperature
{"xtsd", BSIM4_MOD_XTSD}, // Power dependence of JTSD on temperature
{"xtssws", BSIM4_MOD_XTSSWS}, // Power dependence of JTSSWS on temperature
{"xtsswd", BSIM4_MOD_XTSSWD}, // Power dependence of JTSSWD on temperature
{"xtsswgs", BSIM4_MOD_XTSSWGS}, // Power dependence of JTSSWGS on temperature
{"xtsswgd", BSIM4_MOD_XTSSWGD}, // Power dependence of JTSSWGD on temperature
{"tnjts", BSIM4_MOD_TNJTS}, // Temperature coefficient for NJTS
{"tnjtssw", BSIM4_MOD_TNJTSSW}, // Temperature coefficient for NJTSSW
{"tnjtsswg", BSIM4_MOD_TNJTSSWG}, // Temperature coefficient for NJTSSWG
{"tnjtsd", BSIM4_MOD_TNJTSD}, // Temperature coefficient for NJTSD
{"tnjtsswd", BSIM4_MOD_TNJTSSWD}, // Temperature coefficient for NJTSSWD
{"tnjtsswgd", BSIM4_MOD_TNJTSSWGD}, // Temperature coefficient for NJTSSWGD
{"vtss", BSIM4_MOD_VTSS}, // Source bottom trap-assisted voltage dependent parameter
{"vtsd", BSIM4_MOD_VTSD}, // Drain bottom trap-assisted voltage dependent parameter
{"vtssws", BSIM4_MOD_VTSSWS}, // Source STI sidewall trap-assisted voltage dependent parameter
{"vtsswd", BSIM4_MOD_VTSSWD}, // Drain STI sidewall trap-assisted voltage dependent parameter
{"vtsswgs", BSIM4_MOD_VTSSWGS}, // Source gate-edge sidewall trap-assisted voltage dependent parameter
{"vtsswgd", BSIM4_MOD_VTSSWGD}, // Drain gate-edge sidewall trap-assisted voltage dependent parameter

{"gbmin", BSIM4_MOD_GBMIN}, // Minimum body conductance
{"rbdb", BSIM4_MOD_RBDB}, // Resistance between bNode and dbNode
{"rbpb", BSIM4_MOD_RBPB}, // Resistance between bNodePrime and bNode
{"rbsb", BSIM4_MOD_RBSB}, // Resistance between bNode and sbNode
{"rbps", BSIM4_MOD_RBPS}, // Resistance between bNodePrime and sbNode
{"rbpd", BSIM4_MOD_RBPD}, // Resistance between bNodePrime and bNode

{"rbps0", BSIM4_MOD_RBPS0}, // Body resistance RBPS scaling
{"rbpsl", BSIM4_MOD_RBPSL}, // Body resistance RBPS L scaling
{"rbpsw", BSIM4_MOD_RBPSW}, // Body resistance RBPS W scaling
{"rbpsnf", BSIM4_MOD_RBPSNF}, // Body resistance RBPS NF scaling

{"rbpd0", BSIM4_MOD_RBPD0}, // Body resistance RBPD scaling
{"rbpdl", BSIM4_MOD_RBPDL}, // Body resistance RBPD L scaling
{"rbpdw", BSIM4_MOD_RBPDW}, // Body resistance RBPD W scaling
{"rbpdnf", BSIM4_MOD_RBPDNF}, // Body resistance RBPD NF scaling

{"rbpbx0", BSIM4_MOD_RBPBX0}, // Body resistance RBPBX  scaling
{"rbpbxl", BSIM4_MOD_RBPBXL}, // Body resistance RBPBX L scaling
{"rbpbxw", BSIM4_MOD_RBPBXW}, // Body resistance RBPBX W scaling
{"rbpbxnf", BSIM4_MOD_RBPBXNF}, // Body resistance RBPBX NF scaling
{"rbpby0", BSIM4_MOD_RBPBY0}, // Body resistance RBPBY  scaling
{"rbpbyl", BSIM4_MOD_RBPBYL}, // Body resistance RBPBY L scaling
{"rbpbyw", BSIM4_MOD_RBPBYW}, // Body resistance RBPBY W scaling
{"rbpbynf", BSIM4_MOD_RBPBYNF}, // Body resistance RBPBY NF scaling

{"rbsbx0", BSIM4_MOD_RBSBX0}, // Body resistance RBSBX  scaling
{"rbsby0", BSIM4_MOD_RBSBY0}, // Body resistance RBSBY  scaling
{"rbdbx0", BSIM4_MOD_RBDBX0}, // Body resistance RBDBX  scaling
{"rbdby0", BSIM4_MOD_RBDBY0}, // Body resistance RBDBY  scaling

{"rbsdbxl", BSIM4_MOD_RBSDBXL}, // Body resistance RBSDBX L scaling
{"rbsdbxw", BSIM4_MOD_RBSDBXW}, // Body resistance RBSDBX W scaling
{"rbsdbxnf", BSIM4_MOD_RBSDBXNF}, // Body resistance RBSDBX NF scaling
{"rbsdbyl", BSIM4_MOD_RBSDBYL}, // Body resistance RBSDBY L scaling
{"rbsdbyw", BSIM4_MOD_RBSDBYW}, // Body resistance RBSDBY W scaling
{"rbsdbynf", BSIM4_MOD_RBSDBYNF}, // Body resistance RBSDBY NF scaling

{"lcdsc", BSIM4_MOD_LCDSC}, // Length dependence of cdsc
{"lcdscb", BSIM4_MOD_LCDSCB}, // Length dependence of cdscb
{"lcdscd", BSIM4_MOD_LCDSCD}, // Length dependence of cdscd
{"lcit", BSIM4_MOD_LCIT}, // Length dependence of cit
{"lnfactor", BSIM4_MOD_LNFACTOR}, // Length dependence of nfactor
{"lxj", BSIM4_MOD_LXJ}, // Length dependence of xj
{"lvsat", BSIM4_MOD_LVSAT}, // Length dependence of vsat
{"lat", BSIM4_MOD_LAT}, // Length dependence of at
{"la0", BSIM4_MOD_LA0}, // Length dependence of a0
{"lags", BSIM4_MOD_LAGS}, // Length dependence of ags
{"la1", BSIM4_MOD_LA1}, // Length dependence of a1
{"la2", BSIM4_MOD_LA2}, // Length dependence of a2
{"lketa", BSIM4_MOD_LKETA}, // Length dependence of keta
{"lnsub", BSIM4_MOD_LNSUB}, // Length dependence of nsub
{"lndep", BSIM4_MOD_LNDEP}, // Length dependence of ndep
{"lnsd", BSIM4_MOD_LNSD}, // Length dependence of nsd
{"lphin", BSIM4_MOD_LPHIN}, // Length dependence of phin
{"lngate", BSIM4_MOD_LNGATE}, // Length dependence of ngate
{"lgamma1", BSIM4_MOD_LGAMMA1}, // Length dependence of gamma1
{"lgamma2", BSIM4_MOD_LGAMMA2}, // Length dependence of gamma2
{"lvbx", BSIM4_MOD_LVBX}, // Length dependence of vbx
{"lvbm", BSIM4_MOD_LVBM}, // Length dependence of vbm
{"lxt", BSIM4_MOD_LXT}, // Length dependence of xt
{"lk1", BSIM4_MOD_LK1}, // Length dependence of k1
{"lkt1", BSIM4_MOD_LKT1}, // Length dependence of kt1
{"lkt1l", BSIM4_MOD_LKT1L}, // Length dependence of kt1l
{"lkt2", BSIM4_MOD_LKT2}, // Length dependence of kt2
{"lk2", BSIM4_MOD_LK2}, // Length dependence of k2
{"lk3", BSIM4_MOD_LK3}, // Length dependence of k3
{"lk3b", BSIM4_MOD_LK3B}, // Length dependence of k3b
{"lw0", BSIM4_MOD_LW0}, // Length dependence of w0
{"ldvtp0", BSIM4_MOD_LDVTP0}, // Length dependence of dvtp0
{"ldvtp1", BSIM4_MOD_LDVTP1}, // Length dependence of dvtp1
{"ldvtp2", BSIM4_MOD_LDVTP2}, // Length dependence of dvtp2
{"ldvtp3", BSIM4_MOD_LDVTP3}, // Length dependence of dvtp3
{"ldvtp4", BSIM4_MOD_LDVTP4}, // Length dependence of dvtp4
{"ldvtp5", BSIM4_MOD_LDVTP5}, // Length dependence of dvtp5
{"llpe0", BSIM4_MOD_LLPE0}, // Length dependence of lpe0
{"llpeb", BSIM4_MOD_LLPEB}, // Length dependence of lpeb
{"ldvt0", BSIM4_MOD_LDVT0}, // Length dependence of dvt0
{"ldvt1", BSIM4_MOD_LDVT1}, // Length dependence of dvt1
{"ldvt2", BSIM4_MOD_LDVT2}, // Length dependence of dvt2
{"ldvt0w", BSIM4_MOD_LDVT0W}, // Length dependence of dvt0w
{"ldvt1w", BSIM4_MOD_LDVT1W}, // Length dependence of dvt1w
{"ldvt2w", BSIM4_MOD_LDVT2W}, // Length dependence of dvt2w
{"ldrout", BSIM4_MOD_LDROUT}, // Length dependence of drout
{"ldsub", BSIM4_MOD_LDSUB}, // Length dependence of dsub
{"lvth0", BSIM4_MOD_LVTH0}, // Length dependence of vth0
IOPR("lvtho", BSIM4_MOD_LVTH0, IF_REAL,"Length dependence of vtho"),
{"lua", BSIM4_MOD_LUA}, // Length dependence of ua
{"lua1", BSIM4_MOD_LUA1}, // Length dependence of ua1
{"lub", BSIM4_MOD_LUB}, // Length dependence of ub
{"lub1", BSIM4_MOD_LUB1}, // Length dependence of ub1
{"luc", BSIM4_MOD_LUC}, // Length dependence of uc
{"luc1", BSIM4_MOD_LUC1}, // Length dependence of uc1
{"lud", BSIM4_MOD_LUD}, // Length dependence of ud
{"lud1", BSIM4_MOD_LUD1}, // Length dependence of ud1
{"lup", BSIM4_MOD_LUP}, // Length dependence of up
{"llp", BSIM4_MOD_LLP}, // Length dependence of lp
{"lu0", BSIM4_MOD_LU0}, // Length dependence of u0
{"lute", BSIM4_MOD_LUTE}, // Length dependence of ute
{"lucste", BSIM4_MOD_LUCSTE}, // Length dependence of ucste
{"lvoff", BSIM4_MOD_LVOFF}, // Length dependence of voff
{"lminv", BSIM4_MOD_LMINV}, // Length dependence of minv
{"lminvcv", BSIM4_MOD_LMINVCV}, // Length dependence of minvcv
{"ldelta", BSIM4_MOD_LDELTA}, // Length dependence of delta
{"lrdsw", BSIM4_MOD_LRDSW}, // Length dependence of rdsw 
{"lrsw", BSIM4_MOD_LRSW}, // Length dependence of rsw
{"lrdw", BSIM4_MOD_LRDW}, // Length dependence of rdw

{"lprwg", BSIM4_MOD_LPRWG}, // Length dependence of prwg 
{"lprwb", BSIM4_MOD_LPRWB}, // Length dependence of prwb 

{"lprt", BSIM4_MOD_LPRT}, // Length dependence of prt 
{"leta0", BSIM4_MOD_LETA0}, // Length dependence of eta0
{"letab", BSIM4_MOD_LETAB}, // Length dependence of etab
{"lpclm", BSIM4_MOD_LPCLM}, // Length dependence of pclm
{"lpdiblc1", BSIM4_MOD_LPDIBL1}, // Length dependence of pdiblc1
{"lpdiblc2", BSIM4_MOD_LPDIBL2}, // Length dependence of pdiblc2
{"lpdiblcb", BSIM4_MOD_LPDIBLB}, // Length dependence of pdiblcb
{"lfprout", BSIM4_MOD_LFPROUT}, // Length dependence of pdiblcb
{"lpdits", BSIM4_MOD_LPDITS}, // Length dependence of pdits
{"lpditsd", BSIM4_MOD_LPDITSD}, // Length dependence of pditsd
{"lpscbe1", BSIM4_MOD_LPSCBE1}, // Length dependence of pscbe1
{"lpscbe2", BSIM4_MOD_LPSCBE2}, // Length dependence of pscbe2
{"lpvag", BSIM4_MOD_LPVAG}, // Length dependence of pvag
{"lwr", BSIM4_MOD_LWR}, // Length dependence of wr
{"ldwg", BSIM4_MOD_LDWG}, // Length dependence of dwg
{"ldwb", BSIM4_MOD_LDWB}, // Length dependence of dwb
{"lb0", BSIM4_MOD_LB0}, // Length dependence of b0
{"lb1", BSIM4_MOD_LB1}, // Length dependence of b1
{"lcgsl", BSIM4_MOD_LCGSL}, // Length dependence of cgsl
{"lcgdl", BSIM4_MOD_LCGDL}, // Length dependence of cgdl
{"lckappas", BSIM4_MOD_LCKAPPAS}, // Length dependence of ckappas
{"lckappad", BSIM4_MOD_LCKAPPAD}, // Length dependence of ckappad
{"lcf", BSIM4_MOD_LCF}, // Length dependence of cf
{"lclc", BSIM4_MOD_LCLC}, // Length dependence of clc
{"lcle", BSIM4_MOD_LCLE}, // Length dependence of cle
{"lalpha0", BSIM4_MOD_LALPHA0}, // Length dependence of alpha0
{"lalpha1", BSIM4_MOD_LALPHA1}, // Length dependence of alpha1
{"lbeta0", BSIM4_MOD_LBETA0}, // Length dependence of beta0

{"lagidl", BSIM4_MOD_LAGIDL}, // Length dependence of agidl
{"lbgidl", BSIM4_MOD_LBGIDL}, // Length dependence of bgidl
{"lcgidl", BSIM4_MOD_LCGIDL}, // Length dependence of cgidl
{"lrgidl", BSIM4_MOD_LRGIDL}, // Length dependence of rgidl
{"lkgidl", BSIM4_MOD_LKGIDL}, // Length dependence of kgidl
{"lfgidl", BSIM4_MOD_LFGIDL}, // Length dependence of fgidl
{"legidl", BSIM4_MOD_LEGIDL}, // Length dependence of egidl
{"lagisl", BSIM4_MOD_LAGISL}, // Length dependence of agisl
{"lbgisl", BSIM4_MOD_LBGISL}, // Length dependence of bgisl
{"lcgisl", BSIM4_MOD_LCGISL}, // Length dependence of cgisl
{"lrgisl", BSIM4_MOD_LRGISL}, // Length dependence of rgisl
{"lkgisl", BSIM4_MOD_LKGISL}, // Length dependence of kgisl
{"lfgisl", BSIM4_MOD_LFGISL}, // Length dependence of fgisl
{"legisl", BSIM4_MOD_LEGISL}, // Length dependence of egisl
{"laigc", BSIM4_MOD_LAIGC}, // Length dependence of aigc
{"lbigc", BSIM4_MOD_LBIGC}, // Length dependence of bigc
{"lcigc", BSIM4_MOD_LCIGC}, // Length dependence of cigc
{"laigsd", BSIM4_MOD_LAIGSD}, // Length dependence of aigsd
{"lbigsd", BSIM4_MOD_LBIGSD}, // Length dependence of bigsd
{"lcigsd", BSIM4_MOD_LCIGSD}, // Length dependence of cigsd
{"laigs", BSIM4_MOD_LAIGS}, // Length dependence of aigs
{"lbigs", BSIM4_MOD_LBIGS}, // Length dependence of bigs
{"lcigs", BSIM4_MOD_LCIGS}, // Length dependence of cigs
{"laigd", BSIM4_MOD_LAIGD}, // Length dependence of aigd
{"lbigd", BSIM4_MOD_LBIGD}, // Length dependence of bigd
{"lcigd", BSIM4_MOD_LCIGD}, // Length dependence of cigd
{"laigbacc", BSIM4_MOD_LAIGBACC}, // Length dependence of aigbacc
{"lbigbacc", BSIM4_MOD_LBIGBACC}, // Length dependence of bigbacc
{"lcigbacc", BSIM4_MOD_LCIGBACC}, // Length dependence of cigbacc
{"laigbinv", BSIM4_MOD_LAIGBINV}, // Length dependence of aigbinv
{"lbigbinv", BSIM4_MOD_LBIGBINV}, // Length dependence of bigbinv
{"lcigbinv", BSIM4_MOD_LCIGBINV}, // Length dependence of cigbinv
{"lnigc", BSIM4_MOD_LNIGC}, // Length dependence of nigc
{"lnigbinv", BSIM4_MOD_LNIGBINV}, // Length dependence of nigbinv
{"lnigbacc", BSIM4_MOD_LNIGBACC}, // Length dependence of nigbacc
{"lntox", BSIM4_MOD_LNTOX}, // Length dependence of ntox
{"leigbinv", BSIM4_MOD_LEIGBINV}, // Length dependence for eigbinv
{"lpigcd", BSIM4_MOD_LPIGCD}, // Length dependence for pigcd
{"lpoxedge", BSIM4_MOD_LPOXEDGE}, // Length dependence for poxedge

{"lvfbcv", BSIM4_MOD_LVFBCV}, // Length dependence of vfbcv
{"lvfb", BSIM4_MOD_LVFB}, // Length dependence of vfb
{"lacde", BSIM4_MOD_LACDE}, // Length dependence of acde
{"lmoin", BSIM4_MOD_LMOIN}, // Length dependence of moin
{"lnoff", BSIM4_MOD_LNOFF}, // Length dependence of noff
{"lvoffcv", BSIM4_MOD_LVOFFCV}, // Length dependence of voffcv
{"lxrcrg1", BSIM4_MOD_LXRCRG1}, // Length dependence of xrcrg1
{"lxrcrg2", BSIM4_MOD_LXRCRG2}, // Length dependence of xrcrg2
{"llambda", BSIM4_MOD_LLAMBDA}, // Length dependence of lambda
{"lvtl", BSIM4_MOD_LVTL}, //  Length dependence of vtl
{"lxn", BSIM4_MOD_LXN}, //  Length dependence of xn
{"leu", BSIM4_MOD_LEU}, //  Length dependence of eu
{"lucs", BSIM4_MOD_LUCS}, // Length dependence of lucs
{"lvfbsdoff", BSIM4_MOD_LVFBSDOFF}, // Length dependence of vfbsdoff
{"ltvfbsdoff", BSIM4_MOD_LTVFBSDOFF}, // Length dependence of tvfbsdoff
{"ltvoff", BSIM4_MOD_LTVOFF}, // Length dependence of tvoff
{"ltnfactor", BSIM4_MOD_LTNFACTOR}, // Length dependence of tnfactor
{"lteta0", BSIM4_MOD_LTETA0}, // Length dependence of teta0
{"ltvoffcv", BSIM4_MOD_LTVOFFCV}, // Length dependence of tvoffcv

{"wcdsc", BSIM4_MOD_WCDSC}, // Width dependence of cdsc
{"wcdscb", BSIM4_MOD_WCDSCB}, // Width dependence of cdscb
{"wcdscd", BSIM4_MOD_WCDSCD}, // Width dependence of cdscd
{"wcit", BSIM4_MOD_WCIT}, // Width dependence of cit
{"wnfactor", BSIM4_MOD_WNFACTOR}, // Width dependence of nfactor
{"wxj", BSIM4_MOD_WXJ}, // Width dependence of xj
{"wvsat", BSIM4_MOD_WVSAT}, // Width dependence of vsat
{"wat", BSIM4_MOD_WAT}, // Width dependence of at
{"wa0", BSIM4_MOD_WA0}, // Width dependence of a0
{"wags", BSIM4_MOD_WAGS}, // Width dependence of ags
{"wa1", BSIM4_MOD_WA1}, // Width dependence of a1
{"wa2", BSIM4_MOD_WA2}, // Width dependence of a2
{"wketa", BSIM4_MOD_WKETA}, // Width dependence of keta
{"wnsub", BSIM4_MOD_WNSUB}, // Width dependence of nsub
{"wndep", BSIM4_MOD_WNDEP}, // Width dependence of ndep
{"wnsd", BSIM4_MOD_WNSD}, // Width dependence of nsd
{"wphin", BSIM4_MOD_WPHIN}, // Width dependence of phin
{"wngate", BSIM4_MOD_WNGATE}, // Width dependence of ngate
{"wgamma1", BSIM4_MOD_WGAMMA1}, // Width dependence of gamma1
{"wgamma2", BSIM4_MOD_WGAMMA2}, // Width dependence of gamma2
{"wvbx", BSIM4_MOD_WVBX}, // Width dependence of vbx
{"wvbm", BSIM4_MOD_WVBM}, // Width dependence of vbm
{"wxt", BSIM4_MOD_WXT}, // Width dependence of xt
{"wk1", BSIM4_MOD_WK1}, // Width dependence of k1
{"wkt1", BSIM4_MOD_WKT1}, // Width dependence of kt1
{"wkt1l", BSIM4_MOD_WKT1L}, // Width dependence of kt1l
{"wkt2", BSIM4_MOD_WKT2}, // Width dependence of kt2
{"wk2", BSIM4_MOD_WK2}, // Width dependence of k2
{"wk3", BSIM4_MOD_WK3}, // Width dependence of k3
{"wk3b", BSIM4_MOD_WK3B}, // Width dependence of k3b
{"ww0", BSIM4_MOD_WW0}, // Width dependence of w0
{"wdvtp0", BSIM4_MOD_WDVTP0}, // Width dependence of dvtp0
{"wdvtp1", BSIM4_MOD_WDVTP1}, // Width dependence of dvtp1
{"wdvtp2", BSIM4_MOD_WDVTP2}, // Width dependence of dvtp2
{"wdvtp3", BSIM4_MOD_WDVTP3}, // Width dependence of dvtp3
{"wdvtp4", BSIM4_MOD_WDVTP4}, // Width dependence of dvtp4
{"wdvtp5", BSIM4_MOD_WDVTP5}, // Width dependence of dvtp5
{"wlpe0", BSIM4_MOD_WLPE0}, // Width dependence of lpe0
{"wlpeb", BSIM4_MOD_WLPEB}, // Width dependence of lpeb
{"wdvt0", BSIM4_MOD_WDVT0}, // Width dependence of dvt0
{"wdvt1", BSIM4_MOD_WDVT1}, // Width dependence of dvt1
{"wdvt2", BSIM4_MOD_WDVT2}, // Width dependence of dvt2
{"wdvt0w", BSIM4_MOD_WDVT0W}, // Width dependence of dvt0w
{"wdvt1w", BSIM4_MOD_WDVT1W}, // Width dependence of dvt1w
{"wdvt2w", BSIM4_MOD_WDVT2W}, // Width dependence of dvt2w
{"wdrout", BSIM4_MOD_WDROUT}, // Width dependence of drout
{"wdsub", BSIM4_MOD_WDSUB}, // Width dependence of dsub
{"wvth0", BSIM4_MOD_WVTH0}, // Width dependence of vth0
IOPR("wvtho", BSIM4_MOD_WVTH0, IF_REAL,"Width dependence of vtho"),
{"wua", BSIM4_MOD_WUA}, // Width dependence of ua
{"wua1", BSIM4_MOD_WUA1}, // Width dependence of ua1
{"wub", BSIM4_MOD_WUB}, // Width dependence of ub
{"wub1", BSIM4_MOD_WUB1}, // Width dependence of ub1
{"wuc", BSIM4_MOD_WUC}, // Width dependence of uc
{"wuc1", BSIM4_MOD_WUC1}, // Width dependence of uc1
{"wud", BSIM4_MOD_WUD}, // Width dependence of ud
{"wud1", BSIM4_MOD_WUD1}, // Width dependence of ud1
{"wup", BSIM4_MOD_WUP}, // Width dependence of up
{"wlp", BSIM4_MOD_WLP}, // Width dependence of lp
{"wu0", BSIM4_MOD_WU0}, // Width dependence of u0
{"wute", BSIM4_MOD_WUTE}, // Width dependence of ute
{"wucste", BSIM4_MOD_WUCSTE}, // Width dependence of ucste
{"wvoff", BSIM4_MOD_WVOFF}, // Width dependence of voff
{"wminv", BSIM4_MOD_WMINV}, // Width dependence of minv
{"wminvcv", BSIM4_MOD_WMINVCV}, // Width dependence of minvcv
{"wdelta", BSIM4_MOD_WDELTA}, // Width dependence of delta
{"wrdsw", BSIM4_MOD_WRDSW}, // Width dependence of rdsw 
{"wrsw", BSIM4_MOD_WRSW}, // Width dependence of rsw
{"wrdw", BSIM4_MOD_WRDW}, // Width dependence of rdw

{"wprwg", BSIM4_MOD_WPRWG}, // Width dependence of prwg 
{"wprwb", BSIM4_MOD_WPRWB}, // Width dependence of prwb 

{"wprt", BSIM4_MOD_WPRT}, // Width dependence of prt
{"weta0", BSIM4_MOD_WETA0}, // Width dependence of eta0
{"wetab", BSIM4_MOD_WETAB}, // Width dependence of etab
{"wpclm", BSIM4_MOD_WPCLM}, // Width dependence of pclm
{"wpdiblc1", BSIM4_MOD_WPDIBL1}, // Width dependence of pdiblc1
{"wpdiblc2", BSIM4_MOD_WPDIBL2}, // Width dependence of pdiblc2
{"wpdiblcb", BSIM4_MOD_WPDIBLB}, // Width dependence of pdiblcb
{"wfprout", BSIM4_MOD_WFPROUT}, // Width dependence of pdiblcb
{"wpdits", BSIM4_MOD_WPDITS}, // Width dependence of pdits
{"wpditsd", BSIM4_MOD_WPDITSD}, // Width dependence of pditsd
{"wpscbe1", BSIM4_MOD_WPSCBE1}, // Width dependence of pscbe1
{"wpscbe2", BSIM4_MOD_WPSCBE2}, // Width dependence of pscbe2
{"wpvag", BSIM4_MOD_WPVAG}, // Width dependence of pvag
{"wwr", BSIM4_MOD_WWR}, // Width dependence of wr
{"wdwg", BSIM4_MOD_WDWG}, // Width dependence of dwg
{"wdwb", BSIM4_MOD_WDWB}, // Width dependence of dwb
{"wb0", BSIM4_MOD_WB0}, // Width dependence of b0
{"wb1", BSIM4_MOD_WB1}, // Width dependence of b1
{"wcgsl", BSIM4_MOD_WCGSL}, // Width dependence of cgsl
{"wcgdl", BSIM4_MOD_WCGDL}, // Width dependence of cgdl
{"wckappas", BSIM4_MOD_WCKAPPAS}, // Width dependence of ckappas
{"wckappad", BSIM4_MOD_WCKAPPAD}, // Width dependence of ckappad
{"wcf", BSIM4_MOD_WCF}, // Width dependence of cf
{"wclc", BSIM4_MOD_WCLC}, // Width dependence of clc
{"wcle", BSIM4_MOD_WCLE}, // Width dependence of cle
{"walpha0", BSIM4_MOD_WALPHA0}, // Width dependence of alpha0
{"walpha1", BSIM4_MOD_WALPHA1}, // Width dependence of alpha1
{"wbeta0", BSIM4_MOD_WBETA0}, // Width dependence of beta0

{"wagidl", BSIM4_MOD_WAGIDL}, // Width dependence of agidl
{"wbgidl", BSIM4_MOD_WBGIDL}, // Width dependence of bgidl
{"wcgidl", BSIM4_MOD_WCGIDL}, // Width dependence of cgidl
{"wrgidl", BSIM4_MOD_WRGIDL}, // Width dependence of rgidl
{"wkgidl", BSIM4_MOD_WKGIDL}, // Width dependence of kgidl
{"wfgidl", BSIM4_MOD_WFGIDL}, // Width dependence of fgidl
{"wegidl", BSIM4_MOD_WEGIDL}, // Width dependence of egidl
{"wagisl", BSIM4_MOD_WAGISL}, // Width dependence of agisl
{"wbgisl", BSIM4_MOD_WBGISL}, // Width dependence of bgisl
{"wcgisl", BSIM4_MOD_WCGISL}, // Width dependence of cgisl
{"wrgisl", BSIM4_MOD_WRGISL}, // Width dependence of rgisl
{"wkgisl", BSIM4_MOD_WKGISL}, // Width dependence of kgisl
{"wfgisl", BSIM4_MOD_WFGISL}, // Width dependence of fgisl
{"wegisl", BSIM4_MOD_WEGISL}, // Width dependence of egisl
{"waigc", BSIM4_MOD_WAIGC}, // Width dependence of aigc
{"wbigc", BSIM4_MOD_WBIGC}, // Width dependence of bigc
{"wcigc", BSIM4_MOD_WCIGC}, // Width dependence of cigc
{"waigsd", BSIM4_MOD_WAIGSD}, // Width dependence of aigsd
{"wbigsd", BSIM4_MOD_WBIGSD}, // Width dependence of bigsd
{"wcigsd", BSIM4_MOD_WCIGSD}, // Width dependence of cigsd
{"waigs", BSIM4_MOD_WAIGS}, // Width dependence of aigs
{"wbigs", BSIM4_MOD_WBIGS}, // Width dependence of bigs
{"wcigs", BSIM4_MOD_WCIGS}, // Width dependence of cigs
{"waigd", BSIM4_MOD_WAIGD}, // Width dependence of aigd
{"wbigd", BSIM4_MOD_WBIGD}, // Width dependence of bigd
{"wcigd", BSIM4_MOD_WCIGD}, // Width dependence of cigd
{"waigbacc", BSIM4_MOD_WAIGBACC}, // Width dependence of aigbacc
{"wbigbacc", BSIM4_MOD_WBIGBACC}, // Width dependence of bigbacc
{"wcigbacc", BSIM4_MOD_WCIGBACC}, // Width dependence of cigbacc
{"waigbinv", BSIM4_MOD_WAIGBINV}, // Width dependence of aigbinv
{"wbigbinv", BSIM4_MOD_WBIGBINV}, // Width dependence of bigbinv
{"wcigbinv", BSIM4_MOD_WCIGBINV}, // Width dependence of cigbinv
{"wnigc", BSIM4_MOD_WNIGC}, // Width dependence of nigc
{"wnigbinv", BSIM4_MOD_WNIGBINV}, // Width dependence of nigbinv
{"wnigbacc", BSIM4_MOD_WNIGBACC}, // Width dependence of nigbacc
{"wntox", BSIM4_MOD_WNTOX}, // Width dependence of ntox
{"weigbinv", BSIM4_MOD_WEIGBINV}, // Width dependence for eigbinv
{"wpigcd", BSIM4_MOD_WPIGCD}, // Width dependence for pigcd
{"wpoxedge", BSIM4_MOD_WPOXEDGE}, // Width dependence for poxedge
{"wvfbcv", BSIM4_MOD_WVFBCV}, // Width dependence of vfbcv
{"wvfb", BSIM4_MOD_WVFB}, // Width dependence of vfb
{"wacde", BSIM4_MOD_WACDE}, // Width dependence of acde
{"wmoin", BSIM4_MOD_WMOIN}, // Width dependence of moin
{"wnoff", BSIM4_MOD_WNOFF}, // Width dependence of noff
{"wvoffcv", BSIM4_MOD_WVOFFCV}, // Width dependence of voffcv
{"wxrcrg1", BSIM4_MOD_WXRCRG1}, // Width dependence of xrcrg1
{"wxrcrg2", BSIM4_MOD_WXRCRG2}, // Width dependence of xrcrg2
{"wlambda", BSIM4_MOD_WLAMBDA}, // Width dependence of lambda
{"wvtl", BSIM4_MOD_WVTL}, // Width dependence of vtl
{"wxn", BSIM4_MOD_WXN}, // Width dependence of xn
{"weu", BSIM4_MOD_WEU}, // Width dependence of eu
{"wucs", BSIM4_MOD_WUCS}, // Width dependence of ucs
{"wvfbsdoff", BSIM4_MOD_WVFBSDOFF}, // Width dependence of vfbsdoff
{"wtvfbsdoff", BSIM4_MOD_WTVFBSDOFF}, // Width dependence of tvfbsdoff
{"wtvoff", BSIM4_MOD_WTVOFF}, // Width dependence of tvoff
{"wtnfactor", BSIM4_MOD_WTNFACTOR}, // Width dependence of tnfactor
{"wteta0", BSIM4_MOD_WTETA0}, // Width dependence of teta0
{"wtvoffcv", BSIM4_MOD_WTVOFFCV}, // Width dependence of tvoffcv

{"pcdsc", BSIM4_MOD_PCDSC}, // Cross-term dependence of cdsc
{"pcdscb", BSIM4_MOD_PCDSCB}, // Cross-term dependence of cdscb
{"pcdscd", BSIM4_MOD_PCDSCD}, // Cross-term dependence of cdscd
{"pcit", BSIM4_MOD_PCIT}, // Cross-term dependence of cit
{"pnfactor", BSIM4_MOD_PNFACTOR}, // Cross-term dependence of nfactor
{"pxj", BSIM4_MOD_PXJ}, // Cross-term dependence of xj
{"pvsat", BSIM4_MOD_PVSAT}, // Cross-term dependence of vsat
{"pat", BSIM4_MOD_PAT}, // Cross-term dependence of at
{"pa0", BSIM4_MOD_PA0}, // Cross-term dependence of a0
{"pags", BSIM4_MOD_PAGS}, // Cross-term dependence of ags
{"pa1", BSIM4_MOD_PA1}, // Cross-term dependence of a1
{"pa2", BSIM4_MOD_PA2}, // Cross-term dependence of a2
{"pketa", BSIM4_MOD_PKETA}, // Cross-term dependence of keta
{"pnsub", BSIM4_MOD_PNSUB}, // Cross-term dependence of nsub
{"pndep", BSIM4_MOD_PNDEP}, // Cross-term dependence of ndep
{"pnsd", BSIM4_MOD_PNSD}, // Cross-term dependence of nsd
{"pphin", BSIM4_MOD_PPHIN}, // Cross-term dependence of phin
{"pngate", BSIM4_MOD_PNGATE}, // Cross-term dependence of ngate
{"pgamma1", BSIM4_MOD_PGAMMA1}, // Cross-term dependence of gamma1
{"pgamma2", BSIM4_MOD_PGAMMA2}, // Cross-term dependence of gamma2
{"pvbx", BSIM4_MOD_PVBX}, // Cross-term dependence of vbx
{"pvbm", BSIM4_MOD_PVBM}, // Cross-term dependence of vbm
{"pxt", BSIM4_MOD_PXT}, // Cross-term dependence of xt
{"pk1", BSIM4_MOD_PK1}, // Cross-term dependence of k1
{"pkt1", BSIM4_MOD_PKT1}, // Cross-term dependence of kt1
{"pkt1l", BSIM4_MOD_PKT1L}, // Cross-term dependence of kt1l
{"pkt2", BSIM4_MOD_PKT2}, // Cross-term dependence of kt2
{"pk2", BSIM4_MOD_PK2}, // Cross-term dependence of k2
{"pk3", BSIM4_MOD_PK3}, // Cross-term dependence of k3
{"pk3b", BSIM4_MOD_PK3B}, // Cross-term dependence of k3b
{"pw0", BSIM4_MOD_PW0}, // Cross-term dependence of w0
{"pdvtp0", BSIM4_MOD_PDVTP0}, // Cross-term dependence of dvtp0
{"pdvtp1", BSIM4_MOD_PDVTP1}, // Cross-term dependence of dvtp1
{"pdvtp2", BSIM4_MOD_PDVTP2}, // Cross-term dependence of dvtp2
{"pdvtp3", BSIM4_MOD_PDVTP3}, // Cross-term dependence of dvtp3
{"pdvtp4", BSIM4_MOD_PDVTP4}, // Cross-term dependence of dvtp4
{"pdvtp5", BSIM4_MOD_PDVTP5}, // Cross-term dependence of dvtp5
{"plpe0", BSIM4_MOD_PLPE0}, // Cross-term dependence of lpe0
{"plpeb", BSIM4_MOD_PLPEB}, // Cross-term dependence of lpeb
{"pdvt0", BSIM4_MOD_PDVT0}, // Cross-term dependence of dvt0
{"pdvt1", BSIM4_MOD_PDVT1}, // Cross-term dependence of dvt1
{"pdvt2", BSIM4_MOD_PDVT2}, // Cross-term dependence of dvt2
{"pdvt0w", BSIM4_MOD_PDVT0W}, // Cross-term dependence of dvt0w
{"pdvt1w", BSIM4_MOD_PDVT1W}, // Cross-term dependence of dvt1w
{"pdvt2w", BSIM4_MOD_PDVT2W}, // Cross-term dependence of dvt2w
{"pdrout", BSIM4_MOD_PDROUT}, // Cross-term dependence of drout
{"pdsub", BSIM4_MOD_PDSUB}, // Cross-term dependence of dsub
{"pvth0", BSIM4_MOD_PVTH0}, // Cross-term dependence of vth0
IOPR("pvtho", BSIM4_MOD_PVTH0, IF_REAL,"Cross-term dependence of vtho"),
{"pua", BSIM4_MOD_PUA}, // Cross-term dependence of ua
{"pua1", BSIM4_MOD_PUA1}, // Cross-term dependence of ua1
{"pub", BSIM4_MOD_PUB}, // Cross-term dependence of ub
{"pub1", BSIM4_MOD_PUB1}, // Cross-term dependence of ub1
{"puc", BSIM4_MOD_PUC}, // Cross-term dependence of uc
{"puc1", BSIM4_MOD_PUC1}, // Cross-term dependence of uc1
{"pud", BSIM4_MOD_PUD}, // Cross-term dependence of ud
{"pud1", BSIM4_MOD_PUD1}, // Cross-term dependence of ud1
{"pup", BSIM4_MOD_PUP}, // Cross-term dependence of up
{"plp", BSIM4_MOD_PLP}, // Cross-term dependence of lp
{"pu0", BSIM4_MOD_PU0}, // Cross-term dependence of u0
{"pute", BSIM4_MOD_PUTE}, // Cross-term dependence of ute
{"pucste", BSIM4_MOD_PUCSTE}, // Cross-term dependence of ucste
{"pvoff", BSIM4_MOD_PVOFF}, // Cross-term dependence of voff
{"pminv", BSIM4_MOD_PMINV}, // Cross-term dependence of minv
{"pminvcv", BSIM4_MOD_PMINVCV}, // Cross-term dependence of minvcv
{"pdelta", BSIM4_MOD_PDELTA}, // Cross-term dependence of delta
{"prdsw", BSIM4_MOD_PRDSW}, // Cross-term dependence of rdsw 
{"prsw", BSIM4_MOD_PRSW}, // Cross-term dependence of rsw
{"prdw", BSIM4_MOD_PRDW}, // Cross-term dependence of rdw

{"pprwg", BSIM4_MOD_PPRWG}, // Cross-term dependence of prwg 
{"pprwb", BSIM4_MOD_PPRWB}, // Cross-term dependence of prwb 

{"pprt", BSIM4_MOD_PPRT}, // Cross-term dependence of prt 
{"peta0", BSIM4_MOD_PETA0}, // Cross-term dependence of eta0
{"petab", BSIM4_MOD_PETAB}, // Cross-term dependence of etab
{"ppclm", BSIM4_MOD_PPCLM}, // Cross-term dependence of pclm
{"ppdiblc1", BSIM4_MOD_PPDIBL1}, // Cross-term dependence of pdiblc1
{"ppdiblc2", BSIM4_MOD_PPDIBL2}, // Cross-term dependence of pdiblc2
{"ppdiblcb", BSIM4_MOD_PPDIBLB}, // Cross-term dependence of pdiblcb
{"pfprout", BSIM4_MOD_PFPROUT}, // Cross-term dependence of pdiblcb
{"ppdits", BSIM4_MOD_PPDITS}, // Cross-term dependence of pdits
{"ppditsd", BSIM4_MOD_PPDITSD}, // Cross-term dependence of pditsd
{"ppscbe1", BSIM4_MOD_PPSCBE1}, // Cross-term dependence of pscbe1
{"ppscbe2", BSIM4_MOD_PPSCBE2}, // Cross-term dependence of pscbe2
{"ppvag", BSIM4_MOD_PPVAG}, // Cross-term dependence of pvag
{"pwr", BSIM4_MOD_PWR}, // Cross-term dependence of wr
{"pdwg", BSIM4_MOD_PDWG}, // Cross-term dependence of dwg
{"pdwb", BSIM4_MOD_PDWB}, // Cross-term dependence of dwb
{"pb0", BSIM4_MOD_PB0}, // Cross-term dependence of b0
{"pb1", BSIM4_MOD_PB1}, // Cross-term dependence of b1
{"pcgsl", BSIM4_MOD_PCGSL}, // Cross-term dependence of cgsl
{"pcgdl", BSIM4_MOD_PCGDL}, // Cross-term dependence of cgdl
{"pckappas", BSIM4_MOD_PCKAPPAS}, // Cross-term dependence of ckappas
{"pckappad", BSIM4_MOD_PCKAPPAD}, // Cross-term dependence of ckappad
{"pcf", BSIM4_MOD_PCF}, // Cross-term dependence of cf
{"pclc", BSIM4_MOD_PCLC}, // Cross-term dependence of clc
{"pcle", BSIM4_MOD_PCLE}, // Cross-term dependence of cle
{"palpha0", BSIM4_MOD_PALPHA0}, // Cross-term dependence of alpha0
{"palpha1", BSIM4_MOD_PALPHA1}, // Cross-term dependence of alpha1
{"pbeta0", BSIM4_MOD_PBETA0}, // Cross-term dependence of beta0

{"pagidl", BSIM4_MOD_PAGIDL}, // Cross-term dependence of agidl
{"pbgidl", BSIM4_MOD_PBGIDL}, // Cross-term dependence of bgidl
{"pcgidl", BSIM4_MOD_PCGIDL}, // Cross-term dependence of cgidl
{"prgidl", BSIM4_MOD_PRGIDL}, // Cross-term dependence of rgidl
{"pkgidl", BSIM4_MOD_PKGIDL}, // Cross-term dependence of kgidl
{"pfgidl", BSIM4_MOD_PFGIDL}, // Cross-term dependence of fgidl
{"pegidl", BSIM4_MOD_PEGIDL}, // Cross-term dependence of egidl
{"pagisl", BSIM4_MOD_PAGISL}, // Cross-term dependence of agisl
{"pbgisl", BSIM4_MOD_PBGISL}, // Cross-term dependence of bgisl
{"pcgisl", BSIM4_MOD_PCGISL}, // Cross-term dependence of cgisl
{"pegisl", BSIM4_MOD_PEGISL}, // Cross-term dependence of egisl
{"prgisl", BSIM4_MOD_PRGISL}, // Cross-term dependence of rgisl
{"pkgisl", BSIM4_MOD_PKGISL}, // Cross-term dependence of kgisl
{"pfgisl", BSIM4_MOD_PFGISL}, // Cross-term dependence of fgisl
{"paigc", BSIM4_MOD_PAIGC}, // Cross-term dependence of aigc
{"pbigc", BSIM4_MOD_PBIGC}, // Cross-term dependence of bigc
{"pcigc", BSIM4_MOD_PCIGC}, // Cross-term dependence of cigc
{"paigsd", BSIM4_MOD_PAIGSD}, // Cross-term dependence of aigsd
{"pbigsd", BSIM4_MOD_PBIGSD}, // Cross-term dependence of bigsd
{"pcigsd", BSIM4_MOD_PCIGSD}, // Cross-term dependence of cigsd
{"paigs", BSIM4_MOD_PAIGS}, // Cross-term dependence of aigs
{"pbigs", BSIM4_MOD_PBIGS}, // Cross-term dependence of bigs
{"pcigs", BSIM4_MOD_PCIGS}, // Cross-term dependence of cigs
{"paigd", BSIM4_MOD_PAIGD}, // Cross-term dependence of aigd
{"pbigd", BSIM4_MOD_PBIGD}, // Cross-term dependence of bigd
{"pcigd", BSIM4_MOD_PCIGD}, // Cross-term dependence of cigd
{"paigbacc", BSIM4_MOD_PAIGBACC}, // Cross-term dependence of aigbacc
{"pbigbacc", BSIM4_MOD_PBIGBACC}, // Cross-term dependence of bigbacc
{"pcigbacc", BSIM4_MOD_PCIGBACC}, // Cross-term dependence of cigbacc
{"paigbinv", BSIM4_MOD_PAIGBINV}, // Cross-term dependence of aigbinv
{"pbigbinv", BSIM4_MOD_PBIGBINV}, // Cross-term dependence of bigbinv
{"pcigbinv", BSIM4_MOD_PCIGBINV}, // Cross-term dependence of cigbinv
{"pnigc", BSIM4_MOD_PNIGC}, // Cross-term dependence of nigc
{"pnigbinv", BSIM4_MOD_PNIGBINV}, // Cross-term dependence of nigbinv
{"pnigbacc", BSIM4_MOD_PNIGBACC}, // Cross-term dependence of nigbacc
{"pntox", BSIM4_MOD_PNTOX}, // Cross-term dependence of ntox
{"peigbinv", BSIM4_MOD_PEIGBINV}, // Cross-term dependence for eigbinv
{"ppigcd", BSIM4_MOD_PPIGCD}, // Cross-term dependence for pigcd
{"ppoxedge", BSIM4_MOD_PPOXEDGE}, // Cross-term dependence for poxedge
{"pvfbcv", BSIM4_MOD_PVFBCV}, // Cross-term dependence of vfbcv
{"pvfb", BSIM4_MOD_PVFB}, // Cross-term dependence of vfb
{"pacde", BSIM4_MOD_PACDE}, // Cross-term dependence of acde
{"pmoin", BSIM4_MOD_PMOIN}, // Cross-term dependence of moin
{"pnoff", BSIM4_MOD_PNOFF}, // Cross-term dependence of noff
{"pvoffcv", BSIM4_MOD_PVOFFCV}, // Cross-term dependence of voffcv
{"pxrcrg1", BSIM4_MOD_PXRCRG1}, // Cross-term dependence of xrcrg1
{"pxrcrg2", BSIM4_MOD_PXRCRG2}, // Cross-term dependence of xrcrg2
{"plambda", BSIM4_MOD_PLAMBDA}, // Cross-term dependence of lambda
{"pvtl", BSIM4_MOD_PVTL}, // Cross-term dependence of vtl
{"pxn", BSIM4_MOD_PXN}, // Cross-term dependence of xn
{"peu", BSIM4_MOD_PEU}, // Cross-term dependence of eu
{"pucs", BSIM4_MOD_PUCS}, // Cross-term dependence of ucs
{"pvfbsdoff", BSIM4_MOD_PVFBSDOFF}, // Cross-term dependence of vfbsdoff
{"ptvfbsdoff", BSIM4_MOD_PTVFBSDOFF}, // Cross-term dependence of tvfbsdoff
{"ptvoff", BSIM4_MOD_PTVOFF}, // Cross-term dependence of tvoff
{"ptnfactor", BSIM4_MOD_PTNFACTOR}, // Cross-term dependence of tnfactor
{"pteta0", BSIM4_MOD_PTETA0}, // Cross-term dependence of teta0
{"ptvoffcv", BSIM4_MOD_PTVOFFCV}, // Cross-term dependence of tvoffcv

/* stress effect*/
{"saref", BSIM4_MOD_SAREF}, // Reference distance between OD edge to poly of one side
{"sbref", BSIM4_MOD_SBREF}, // Reference distance between OD edge to poly of the other side
{"wlod", BSIM4_MOD_WLOD}, // Width parameter for stress effect
{"ku0", BSIM4_MOD_KU0}, // Mobility degradation/enhancement coefficient for LOD
{"kvsat", BSIM4_MOD_KVSAT}, // Saturation velocity degradation/enhancement parameter for LOD
{"kvth0", BSIM4_MOD_KVTH0}, // Threshold degradation/enhancement parameter for LOD
{"tku0", BSIM4_MOD_TKU0}, // Temperature coefficient of KU0
{"llodku0", BSIM4_MOD_LLODKU0}, // Length parameter for u0 LOD effect
{"wlodku0", BSIM4_MOD_WLODKU0}, // Width parameter for u0 LOD effect
{"llodvth", BSIM4_MOD_LLODVTH}, // Length parameter for vth LOD effect
{"wlodvth", BSIM4_MOD_WLODVTH}, // Width parameter for vth LOD effect
{"lku0", BSIM4_MOD_LKU0}, // Length dependence of ku0
{"wku0", BSIM4_MOD_WKU0}, // Width dependence of ku0
{"pku0", BSIM4_MOD_PKU0}, // Cross-term dependence of ku0
{"lkvth0", BSIM4_MOD_LKVTH0}, // Length dependence of kvth0
{"wkvth0", BSIM4_MOD_WKVTH0}, // Width dependence of kvth0
{"pkvth0", BSIM4_MOD_PKVTH0}, // Cross-term dependence of kvth0
{"stk2", BSIM4_MOD_STK2}, // K2 shift factor related to stress effect on vth
{"lodk2", BSIM4_MOD_LODK2}, // K2 shift modification factor for stress effect
{"steta0", BSIM4_MOD_STETA0}, // eta0 shift factor related to stress effect on vth
{"lodeta0", BSIM4_MOD_LODETA0}, // eta0 shift modification factor for stress effect
/* Well Proximity Effect */
{"web", BSIM4_MOD_WEB}, // Coefficient for SCB
{"wec", BSIM4_MOD_WEC}, // Coefficient for SCC
{"kvth0we", BSIM4_MOD_KVTH0WE}, // Threshold shift factor for well proximity effect
{"k2we", BSIM4_MOD_K2WE}, //  K2 shift factor for well proximity effect 
{"ku0we", BSIM4_MOD_KU0WE}, //  Mobility degradation factor for well proximity effect 
{"scref", BSIM4_MOD_SCREF}, //  Reference distance to calculate SCA, SCB and SCC
{"wpemod", BSIM4_MOD_WPEMOD}, //  Flag for WPE model (WPEMOD=1 to activate this model) 
{"lkvth0we", BSIM4_MOD_LKVTH0WE}, // Length dependence of kvth0we
{"lk2we", BSIM4_MOD_LK2WE}, //  Length dependence of k2we 
{"lku0we", BSIM4_MOD_LKU0WE}, //  Length dependence of ku0we 
{"wkvth0we", BSIM4_MOD_WKVTH0WE}, // Width dependence of kvth0we
{"wk2we", BSIM4_MOD_WK2WE}, //  Width dependence of k2we 
{"wku0we", BSIM4_MOD_WKU0WE}, //  Width dependence of ku0we 
{"pkvth0we", BSIM4_MOD_PKVTH0WE}, // Cross-term dependence of kvth0we
{"pk2we", BSIM4_MOD_PK2WE}, //  Cross-term dependence of k2we 
{"pku0we", BSIM4_MOD_PKU0WE}, //  Cross-term dependence of ku0we 

{"noia", BSIM4_MOD_NOIA}, // Flicker noise parameter
{"noib", BSIM4_MOD_NOIB}, // Flicker noise parameter
{"noic", BSIM4_MOD_NOIC}, // Flicker noise parameter
{"tnoia", BSIM4_MOD_TNOIA}, // Thermal noise parameter
{"tnoib", BSIM4_MOD_TNOIB}, // Thermal noise parameter
{"tnoic", BSIM4_MOD_TNOIC}, // Thermal noise parameter
{"rnoia", BSIM4_MOD_RNOIA}, // Thermal noise coefficient
{"rnoib", BSIM4_MOD_RNOIB}, // Thermal noise coefficient
{"rnoic", BSIM4_MOD_RNOIC}, // Thermal noise coefficient
{"gidlclamp", BSIM4_MOD_GIDLCLAMP}, // gidl clamp value
{"idovvds", BSIM4_MOD_IDOVVDSC}, // noise clamping limit parameter
{"ntnoi", BSIM4_MOD_NTNOI}, // Thermal noise parameter
{"em", BSIM4_MOD_EM}, // Flicker noise parameter
{"ef", BSIM4_MOD_EF}, // Flicker noise frequency exponent
{"af", BSIM4_MOD_AF}, // Flicker noise exponent
{"kf", BSIM4_MOD_KF}, // Flicker noise coefficient

// {"vgs_max", BSIM4_MOD_VGS_MAX}, // maximum voltage G-S branch
// {"vgd_max", BSIM4_MOD_VGD_MAX}, // maximum voltage G-D branch
// {"vgb_max", BSIM4_MOD_VGB_MAX}, // maximum voltage G-B branch
// {"vds_max", BSIM4_MOD_VDS_MAX}, // maximum voltage D-S branch
// {"vbs_max", BSIM4_MOD_VBS_MAX}, // maximum voltage B-S branch
// {"vbd_max", BSIM4_MOD_VBD_MAX}, // maximum voltage B-D branch
// {"vgsr_max", BSIM4_MOD_VGSR_MAX}, // maximum voltage G-S branch
// {"vgdr_max", BSIM4_MOD_VGDR_MAX}, // maximum voltage G-D branch
// {"vgbr_max", BSIM4_MOD_VGBR_MAX}, // maximum voltage G-B branch
// {"vbsr_max", BSIM4_MOD_VBSR_MAX}, // maximum voltage B-S branch
// {"vbdr_max", BSIM4_MOD_VBDR_MAX}, // maximum voltage B-D branch

{"nmos", BSIM4_MOD_NMOS}, // Flag to indicate NMOS
{"pmos", BSIM4_MOD_PMOS}, // Flag to indicate PMOS
};

} // namespace bsim4