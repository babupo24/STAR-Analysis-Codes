#include <iostream>
#include "TH1.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include <cmath> 
#include "TCanvas.h"
#include "TLatex.h"

using namespace std;
void setPad(int i); 
double calErr(double a, double ea, double b, double eb);
void  asymPt_AfterPIDRelativeDiff(){
 gStyle -> SetOptStat(0);
 //gStyle -> SetOptFit(0);
 //gStyle -> SetLegendBorderSize(0);
 //-1<nSigmaPion<2
double  pT1[9]= {3.26865, 3.57259, 3.84804, 4.14584, 4.48818, 4.90199, 5.45597, 6.32023, 8.59682};
double  pT2[9]= {3.26189, 3.57787, 3.86589, 4.17821, 4.54056, 4.98886, 5.59061, 6.52354, 8.99652};
double  pT3[9]= {3.25046, 3.58372, 3.89162, 4.22767, 4.6186, 5.10109, 5.74469, 6.73754, 9.3112};
double  pT4[9]= {3.19279, 3.50986, 3.80316, 4.12717, 4.51068, 4.99119, 5.64314, 6.65054, 9.27631};
double  pT5[9]= {3.49378, 3.95072, 4.34082, 4.75126, 5.21713, 5.78142, 6.52177, 7.64473, 10.3651};


 double avgM[5]={0.39,0.58,0.76,0.96,1.42};


 // with start less tof 
 // polynomial TOF cut applied (F-inal result)  
 //<<<<<<<<<<<       Eta < 0, Weighted Average A_{UT}   >>>>>>>>>>>>>>>>>
 double  A_F1Lt[9]    	= {0.00348548, -0.00320122, 0.00412137, 0.00234069, 0.000779033, 0.000201483, -0.00105705, 0.00130495, 0.000728564};
 double  AWerr_F1Lt[9]	= {0.00375107, 0.00362694, 0.00358834, 0.00346801, 0.00347846, 0.00344461, 0.00341071, 0.00338781, 0.00342116}     ;
 double  A_F2Lt[9]	= {0.00442192, -0.00179046, -0.00107806, 0.00371675, -0.00159781, 0.00180338, -0.00177111, 0.0100875, 0.00481067}  ;
 double  AWerr_F2Lt[9] 	= {0.00378167, 0.00366882, 0.00358811, 0.00351899, 0.00346813, 0.00342893, 0.00338408, 0.00335888, 0.00331322}     ;
 double A_F3Lt[9] 	= {-0.00505565, 0.0080467, -0.00171752, 0.0051388, 0.0116025, 0.00814332, -0.00218784, 0.00502388, 0.00809537}     ;
 double AWerr_F3Lt[9] 	= {0.00380225, 0.00368214, 0.00357991, 0.00350593, 0.00345122, 0.00339904, 0.00335227, 0.00331955, 0.00327086}     ;
 double  A_F4Lt[9] 	= {0.00341239, -0.00412455, -0.000653273, -0.00261282, 0.00695626, 0.00565052, 0.00556933, 0.00887667, 0.00747906} ;
 double AWerr_F4Lt[9] 	= {0.00400277, 0.00384886, 0.00372137, 0.00362233, 0.00353502, 0.00347057, 0.00340738, 0.00335546, 0.00326837}     ;
 double  A_F5Lt[9] 	= {0.00454411, 0.00680305, -0.000636729, 0.00281044, -0.00182767, -0.000549067, 0.00608078, 0.00593821, 0.0102812} ;
 double AWerr_F5Lt[9] 	= {0.00398988, 0.00381732, 0.00371753, 0.00362689, 0.00356792, 0.00352274, 0.00347634, 0.00345971, 0.00341648}     ;
 //<<<<<<<<<<<       Eta >0, Average A_{UT}   >>>>>>>>>>>>>>>>>
 double  A_F1Gt[9] 	= {0.00229973, 0.00463731, 0.00013667, 0.00804334, 0.00453314, 0.00915381, 0.00352464, 0.00667502, 0.0100966}      ;
 double  AWerr_F1Gt[9] 	= {0.0037528, 0.00363439, 0.00357339, 0.00347235, 0.003483, 0.0034433, 0.00340907, 0.00338246, 0.00342133}         ;
 double  A_F2Gt[9] 	= {0.00451102, 0.00866882, 0.0111155, 0.017502, 0.014146, 0.0148322, 0.0156611, 0.0153059, 0.0249271}              ;
 double  AWerr_F2Gt[9]  = {0.00378481, 0.00367466, 0.00360961, 0.00354, 0.0034761, 0.00343937, 0.00340406, 0.00336604, 0.00335045}         ;
 double  A_F3Gt[9]  	= {0.00895856, 0.00588134, 0.0102752, 0.0181503, 0.0146603, 0.0201062, 0.029061, 0.0266236, 0.053982}              ;
 double  AWerr_F3Gt[9]  = {0.00380873, 0.00368305, 0.00358678, 0.00352513, 0.00345702, 0.00341665, 0.00340377, 0.00335804, 0.00343369}     ;
 double  A_F4Gt[9]  	= {0.00302912, 0.0149648, 0.00498253, 0.00442801, 0.0106955, 0.0226177, 0.0225974, 0.0289299, 0.056158}            ;
 double  AWerr_F4Gt[9] 	= {0.00401012, 0.00386289, 0.00373546, 0.00363212, 0.00354491, 0.00349476, 0.00343496, 0.00340427, 0.00344471}     ;
 double  A_F5Gt[9] 	= {-0.00329592, 0.00257373, -0.00230884, 0.00582794, -0.000161914, 0.00175628, 0.0160709, 0.0206936, 0.0369948}    ;
 double  AWerr_F5Gt[9] 	= {0.00399831, 0.00382364, 0.00371902, 0.00363197, 0.00357699, 0.00353266, 0.00350116, 0.00347386, 0.00343944}     ;
 
//------Final result<<<


 //<<<<<<<<<<<  TPC only     Eta < 0, Weighted Average A_{UT}   >>>>>>>>>>>>>>>>>
 double  A_Tpc1Lt[9]     = {0.00148781, -0.0041467, 0.00239031, 0.00167791, 0.000242881, 0.00105581, -0.000699193, 0.000816474, 0.000658892};
 double  AWerr_Tpc1Lt[9] = {0.00353169, 0.00346298, 0.00345464, 0.00336567, 0.00339766, 0.00338129, 0.00336233, 0.00335055, 0.00339171}     ;
 double  A_Tpc2Lt[9]     = {0.00321058, -0.00264875, -0.000661371, 0.00404281, -0.00122838, 0.00268216, -0.00148145, 0.00931278, 0.00537388};
 double  AWerr_Tpc2Lt[9] = {0.00353949, 0.00348286, 0.00344461, 0.00340933, 0.0033847, 0.00336514, 0.00333537, 0.0033211, 0.00328545}       ;
 double A_Tpc3Lt[9]      = {-0.00399108, 0.00680065, -0.00256975, 0.00491628, 0.0111401, 0.00844275, -0.0014852, 0.00585025, 0.00756285}    ;
 double AWerr_Tpc3Lt[9]  = {0.00354232, 0.00348939, 0.00343825, 0.00339949, 0.0033721, 0.0033413, 0.00330662, 0.00328497, 0.00324258}       ;
 double  A_Tpc4Lt[9]     = {0.00140875, -0.00433457, -2.09806e-05, -0.00203192, 0.00662898, 0.00507674, 0.00514478, 0.00849403, 0.00745116} ;
 double AWerr_Tpc4Lt[9]  = {0.00365862, 0.00358002, 0.00351588, 0.00346592, 0.00342095, 0.00338507, 0.003345, 0.00331019, 0.0032363}        ;
 double  A_Tpc5Lt[9]     = {0.00402429, 0.00640172, -0.00222258, 0.00339737, -0.000664207, 0.000214951, 0.00593108, 0.00658557, 0.00900641} ;
 double AWerr_Tpc5Lt[9]  = {0.00367489, 0.00359106, 0.0035469, 0.00349777, 0.00346187, 0.00343715, 0.00340584, 0.00339704, 0.00336598} ;
 //<<<<<<<<<<<  TPConly     Eta > 0, Average A_{UT}   >>>>>>>>>>>>>>>>>
 double  A_Tpc1Gt[9]     = {0.0025708, 0.0065227, -0.000246177, 0.00896932, 0.00508233, 0.00884426, 0.00357167, 0.00658349, 0.00986304};
 double  AWerr_Tpc1Gt[9] = {0.00353532, 0.00346848, 0.00344198, 0.00337176, 0.00340507, 0.00338151, 0.00336303, 0.00334573, 0.00339209};
 double  A_Tpc2Gt[9]     = {0.0054225, 0.00699914, 0.00914777, 0.0161481, 0.0126812, 0.014782, 0.0150711, 0.0152443, 0.0253258}        ;
 double  AWerr_Tpc2Gt[9] = {0.00354172, 0.00348835, 0.00346254, 0.00343079, 0.00339092, 0.00337701, 0.00335519, 0.00332823, 0.00332341};
 double  A_Tpc3Gt[9]     = {0.00855265, 0.0044207, 0.00986776, 0.0176423, 0.0142387, 0.01956, 0.0271261, 0.0270546, 0.0536495}         ;
 double  AWerr_Tpc3Gt[9] = {0.00355233, 0.00348871, 0.00344308, 0.00341857, 0.00337911, 0.00335786, 0.00335327, 0.00332279, 0.00340313};
 double  A_Tpc4Gt[9]     = {0.00501529, 0.010025, 0.00351953, 0.00512257, 0.0102015, 0.0223209, 0.0233967, 0.0295134, 0.0558711}       ;
 double  AWerr_Tpc4Gt[9] = {0.00366139, 0.00358495, 0.00352766, 0.00347483, 0.00343158, 0.00341303, 0.00337754, 0.00336049, 0.00341318};
 double  A_Tpc5Gt[9]     = {-0.00424291, 0.000871567, -0.00244794, 0.00568318, 0.000821109, 0.00226176, 0.016264, 0.0207675, 0.036248} ;
 double  AWerr_Tpc5Gt[9] = {0.00367703, 0.00359753, 0.00355175, 0.00350256, 0.00347305, 0.00344637, 0.00343129, 0.00341345, 0.00339174};

 // TOF only eta>0
 //<<<<<<<<<<<       Eta < 0, Weighted Average A_{UT}   >>>>>>>>>>>>>>>>>
 double  A_Tof1Lt[9]     = {-0.00426733, -0.00522767, -0.000345471, 0.00251373, -0.000633467, 0.0028124, -0.00413704, 0.00225719, -0.000921868};
 double  AWerr_Tof1Lt[9] = {0.0057057, 0.0055578, 0.00551057, 0.00536831, 0.00538514, 0.00534688, 0.00532028, 0.00529064, 0.00535289}          ;
 double  A_Tof2Lt[9]     = {-0.00203451, 0.000396127, 0.00556187, 0.00442653, -0.00489747, -0.00128452, -0.0019166, 0.00895087, -0.000680901}  ;
 double  AWerr_Tof2Lt[9] = {0.0057332, 0.00560913, 0.0055089, 0.00542678, 0.0053439, 0.00529229, 0.00524884, 0.00520203, 0.0051848}            ;
 double A_Tof3Lt[9]      = {-0.00626754, 0.00500823, -0.00602793, 0.00586334, 0.0132168, -0.000250217, -0.00169531, 0.00469488, 0.00918192}    ;
 double AWerr_Tof3Lt[9]  = {0.00573771, 0.00560565, 0.00550786, 0.00540325, 0.00533959, 0.00526553, 0.00518804, 0.00513517, 0.00507258}        ;
 double  A_Tof4Lt[9]     = {-0.00344903, -0.00215106, -0.0109362, -0.00152686, 0.00314727, 0.00405867, 0.00307198, 0.00683286, 0.00681364}     ;
 double AWerr_Tof4Lt[9]  = {0.00597989, 0.00579122, 0.00565579, 0.00553056, 0.00544224, 0.00533743, 0.00525782, 0.00519273, 0.00509074}        ;
 double  A_Tof5Lt[9]     = {0.00305685, 0.00624364, -0.00406964, -0.00252273, -0.00428285, 0.00115815, 0.00674366, -0.000278485, 0.00958226}   ;
 double AWerr_Tof5Lt[9]  = {0.00598696, 0.0057939, 0.00567837, 0.00556203, 0.00550167, 0.00544786, 0.00538952, 0.00525469, 0.00519169}         ;
 //<<<<<<<<<<<       Eta >00, Average A_{UT}   >>>>>>>>>>>>>>>>>
 double  A_Tof1Gt[9] 	 = {0.003594, 0.00493318, -0.00498515, 0.009226, 0.00278908, 0.00664479, 0.00232062, 0.00353825, 0.00893677}           ;
 double  AWerr_Tof1Gt[9] = {0.00570282, 0.00556635, 0.00550422, 0.00535864, 0.0053899, 0.00535337, 0.00531021, 0.00527827, 0.00535512}         ;
 double  A_Tof2Gt[9]	 = {0.0058868, 0.00817059, 0.00704804, 0.0134224, 0.0142114, 0.0159216, 0.0127315, 0.00847549, 0.027675}               ;
 double  AWerr_Tof2Gt[9] = {0.00573343, 0.00560965, 0.00551577, 0.0054289, 0.00535289, 0.00530828, 0.00524881, 0.00521671, 0.00522481}         ;
 double  A_Tof3Gt[9] 	 = {0.00887719, 0.0029311, 0.0101039, 0.0268273, 0.0131828, 0.0165602, 0.0270581, 0.0240076, 0.0545419}                ;
 double  AWerr_Tof3Gt[9] = {0.00574342, 0.00561553, 0.00550466, 0.00543578, 0.00532502, 0.00526682, 0.00521588, 0.00514811, 0.00518769}        ;
 double  A_Tof4Gt[9] 	 = {0.00443251, 0.00938992, 0.000808333, 0.00378168, 0.00268466, 0.0158177, 0.0200894, 0.0258779, 0.0538489}           ;
 double  AWerr_Tof4Gt[9] = {0.00598855, 0.00580144, 0.00565604, 0.00555135, 0.00543907, 0.00535029, 0.00527424, 0.00522606, 0.00518815}        ;
 double  A_Tof5Gt[9] 	 = {-0.0154572, 0.00343071, 0.00111607, 0.00800457, -0.000412353, 0.00223463, 0.0168569, 0.0228139, 0.0374865}         ;
 double  AWerr_Tof5Gt[9] = {0.00601495, 0.00581361, 0.00570355, 0.00558529, 0.00552049, 0.00545801, 0.00541073, 0.00535327, 0.00529594}        ;
 //pair fractions for TOF&&TPC only----------->>
 //pion
 double pionPairPtGtT[5]={0.995685,0.99413,0.992819,0.981137,0.951615};
 double pionPairPtLtT[5]={0.997217,0.995865,0.994613,0.982959,0.954391};
 double pionPairMGtT[5]={0.978901,0.98193,0.981563,0.980261,0.971465};
 double pionPairMLtT[5]={0.981487,0.984854,0.98391,0.982889,0.973322};
 //kaon+proton
 double kpPairPtGtT[5]={0,0,0,2.08895e-05,0.000131916};
 double kpPairPtLtT[5]={0,0,0,1.90933e-05,0.000122819};
 double kpPairMGtT[5]={1.88864e-05,3.60047e-05,3.40103e-05,4.00047e-05,7.8545e-05};
 double kpPairMLtT[5]={1.49303e-05,2.6513e-05,2.64161e-05,3.2977e-05,7.32954e-05};
 // pair purity in eta bins
 double pionPairEtaT[9]={0.97593,0.98328,0.982627,0.978905,0.975671,0.978282,0.983113,0.983608,0.976228};
 double kpPairEtaT[9]={3.55615e-05,2.98475e-05,3.96205e-05,5.86822e-05,7.06504e-05,5.12718e-05,2.85773e-05,2.04411e-05,2.54495e-05};
 double electronPairEtaT[9]={3.69162e-05,8.40101e-06,5.82411e-06,8.56687e-06,1.40317e-05,1.3517e-05,9.81749e-06,1.35601e-05,4.67372e-05};
 // pair purity in integrated pT bins
 double pionPairIntMGtT[9]={0.975651,0.987897,0.985057,0.981433,0.984766,0.983844,0.98232,0.979519,0.967469};
 double pionPairIntMLtT[9]={0,0,0,0,0,0,0,0,0};
 double kpPairIntMGtT[9]={1.03398e-05,1.02371e-05,1.9807e-05,3.59803e-05,1.7999e-05,2.05885e-05,2.80107e-05,4.46095e-05,9.39217e-05};
 double kpPairIntMLtT[9]={0,0,0,0,0,0,0,0,0};
 double electronPairIntMGtT[9]={8.00411e-05,8.2071e-06,9.2425e-06,1.09599e-05,1.15676e-05,1.26905e-05,1.28323e-05,1.30567e-05,4.49852e-05};
 double electronPairIntMLtT[9]={0,0,0,0,0,0,0,0,0};
 //---------------------------------------->>
 //pair fractions for TOF Or TPC ---------------->>
 //pion
 double pionPairPtGt[5]={0.928937,0.9147,0.899586,0.879687,0.840236};
 double pionPairPtLt[5]={0.935433,0.928957,0.912599,0.892914,0.853353};
 double pionPairMGt[5]={0.878543,0.884971,0.895744,0.878525,0.849068};
 double pionPairMLt[5]={0.89301,0.899083,0.907726,0.891743,0.863924};
 //kaon+proton
 double kpPairPtGt[5]={0.00124171,0.00175025,0.00232503,0.00307886,0.00464377};
 double kpPairPtLt[5]={0.000981375,0.00119706,0.00174034,0.00239529,0.00377215};
 double kpPairMGt[5]={0.0029743,0.00298951,0.00236606,0.00333788,0.0051196};
 double kpPairMLt[5]={0.00222713,0.00227888,0.0018332,0.00262824,0.0040893};
 // pair purity in eta bins
 double pionPairEta[9]={0.861208,0.901232,0.902267,0.894404,0.885842,0.889964,0.895621,0.890979,0.843378};
 double kpPairEta[9]={0.00370107,0.00222934,0.00222452,0.00252551,0.00280255,0.00266003,0.00248695,0.00269267,0.00505025};
 double electronPairEta[9]={0.000119078,1.14592e-05,8.60324e-06,1.57037e-05,3.27612e-05,2.43846e-05,1.40358e-05,1.73235e-05,0.000108841};
 // pair purity in integrated pT bins
 double pionPairIntMGt[9]={0.866901,0.896771,0.888716,0.886221,0.903343,0.892411,0.879951,0.86373,0.837623};
 double pionPairIntMLt[9]={0.876213,0.90534,0.898047,0.895626,0.90986,0.899965,0.888411,0.873801,0.849522};
 double kpPairIntMGt[9]={0.00332073,0.00236653,0.0027781,0.00289667,0.0019887,0.00253699,0.00324735,0.00429827,0.00587125};
 double kpPairIntMLt[9]={0.00279147,0.00200249,0.00233454,0.00244492,0.00173725,0.00219747,0.0028164,0.00366644,0.00498812};
 double electronPairIntMGt[9]={0.000126754,1.87604e-05,2.07219e-05,2.25915e-05,2.4395e-05,2.4383e-05,2.42277e-05,2.55197e-05,6.62633e-05};
 double electronPairIntMLt[9]={0.000122624,1.35428e-05,1.56263e-05,1.65356e-05,1.91616e-05,1.92729e-05,1.85872e-05,2.14145e-05,5.82379e-05};
 //------------------------------------------>>
 //trigger bias
double  ratiopTGt1[9]={0.9543,0.9543,0.9543,0.9543,0.9543,0.9543,0.9543,0.9543,0.9543};
double  ratiopTGt2[9]={1.001,1.001,1.001,1.001,1.001,1.001,1.001,1.001,1.001};
double  ratiopTGt3[9]={0.9893,0.9893,0.9893,0.9893,0.9893,0.9893,0.9893,0.9893,0.9893};
double  ratiopTGt4[9]={1.064,1.064,1.064,1.064,1.064,1.064,1.064,1.064,1.064};
double  ratiopTGt5[9]={1.063,1.063,1.063,1.063,1.063,1.063,1.063,1.063,1.063};

double  ratiopTLt1[9]={1.056,1.056,1.056,1.056,1.056,1.056,1.056,1.056,1.056};
double  ratiopTLt2[9]={1.044,1.044,1.044,1.044,1.044,1.044,1.044,1.044,1.044};
double  ratiopTLt3[9]={0.9292,0.9292,0.9292,0.9292,0.9292,0.9292,0.9292,0.9292,0.9292};
double  ratiopTLt4[9]={1.049,1.049,1.049,1.049,1.049,1.049,1.049,1.049,1.049};
double  ratiopTLt5[9]={1.134,1.134,1.134,1.134,1.134,1.134,1.134,1.134,1.134};
 
double biaspTGt1[9]={0};double biaspTGt2[9]={0};double biaspTGt3[9]={0};double biaspTGt4[9]={0};double biaspTGt5[9]={0};
 double biaspTLt1[9]={0};double biaspTLt2[9]={0};double biaspTLt3[9]={0};double biaspTLt4[9]={0};double biaspTLt5[9]={0};


 //relataive difference in asymmetry between only TPC and only TOF
double rrdiff1Gt[9]={0}; double rrdiff2Gt[9]={0}; double rrdiff3Gt[9]={0}; double  rrdiff4Gt[9]={0}; double  rrdiff5Gt[9]={0};
double rrdiff1errGt[9]={0}; double rrdiff2errGt[9]={0}; double  rrdiff3errGt[9]={0}; double  rrdiff4errGt[9]={0}; double  rrdiff5errGt[9]={0};
double rrdiff1Lt[9]={0}; double rrdiff2Lt[9]={0}; double rrdiff3Lt[9]={0}; double  rrdiff4Lt[9]={0}; double  rrdiff5Lt[9]={0};
double rrdiff1errLt[9]={0}; double rrdiff2errLt[9]={0}; double  rrdiff3errLt[9]={0}; double  rrdiff4errLt[9]={0}; double  rrdiff5errLt[9]={0};

 double rdiff1Gt[9]={0}; double rdiff2Gt[9]={0}; double rdiff3Gt[9]={0}; double  rdiff4Gt[9]={0}; double  rdiff5Gt[9]={0};
 double rdiff1errGt[9]={0}; double rdiff2errGt[9]={0}; double  rdiff3errGt[9]={0}; double  rdiff4errGt[9]={0}; double  rdiff5errGt[9]={0};
 double trdiff1Gt[9]={0}; double trdiff2Gt[9]={0}; double trdiff3Gt[9]={0}; double  trdiff4Gt[9]={0}; double  trdiff5Gt[9]={0};
 double trdiff1errGt[9]={0}; double trdiff2errGt[9]={0}; double  trdiff3errGt[9]={0}; double  trdiff4errGt[9]={0}; double  trdiff5errGt[9]={0};
 double rdiff1Lt[9]={0}; double rdiff2Lt[9]={0}; double rdiff3Lt[9]={0}; double  rdiff4Lt[9]={0}; double  rdiff5Lt[9]={0};
 double rdiff1errLt[9]={0}; double rdiff2errLt[9]={0}; double  rdiff3errLt[9]={0}; double  rdiff4errLt[9]={0}; double  rdiff5errLt[9]={0};
 double trdiff1Lt[9]={0}; double trdiff2Lt[9]={0}; double trdiff3Lt[9]={0}; double  trdiff4Lt[9]={0}; double  trdiff5Lt[9]={0};
 double trdiff1errLt[9]={0}; double trdiff2errLt[9]={0}; double  trdiff3errLt[9]={0}; double  trdiff4errLt[9]={0}; double  trdiff5errLt[9]={0};
 //pid systematic uncertainty
 double pidSysGt1[9]={0};double pidSysGt2[9]={0};double pidSysGt3[9]={0};double pidSysGt4[9]={0};double pidSysGt5[9]={0};
 double fpidSysGt1[9]={0};double fpidSysGt2[9]={0};double fpidSysGt3[9]={0};double fpidSysGt4[9]={0};double fpidSysGt5[9]={0};
 double fp0Gt[5]={0};
 double pidSysLt1[9]={0};double pidSysLt2[9]={0};double pidSysLt3[9]={0};double pidSysLt4[9]={0};double pidSysLt5[9]={0};
 double fpidSysLt1[9]={0};double fpidSysLt2[9]={0};double fpidSysLt3[9]={0};double fpidSysLt4[9]={0};double fpidSysLt5[9]={0};
 double fp0Lt[5]={0};

 double combSyspTGt1[9]={0};double combSyspTGt2[9]={0};double combSyspTGt3[9]={0};double combSyspTGt4[9]={0};double combSyspTGt5[9]={0};
 double combSyspTLt1[9]={0};double combSyspTLt2[9]={0};double combSyspTLt3[9]={0};double combSyspTLt4[9]={0};double combSyspTLt5[9]={0};

 for(int i=0; i<9; i++){
  //trigger bias eta_pair>0
  biaspTGt1[i]=(1-ratiopTGt1[i])*TMath::Max(A_Tpc1Gt[i],AWerr_Tpc1Gt[i]);
  biaspTGt2[i]=(1-ratiopTGt2[i])*TMath::Max(A_Tpc2Gt[i],AWerr_Tpc2Gt[i]);
  biaspTGt3[i]=(1-ratiopTGt3[i])*TMath::Max(A_Tpc3Gt[i],AWerr_Tpc3Gt[i]);
  biaspTGt4[i]=(1-ratiopTGt4[i])*TMath::Max(A_Tpc4Gt[i],AWerr_Tpc4Gt[i]);
  biaspTGt5[i]=(1-ratiopTGt5[i])*TMath::Max(A_Tpc5Gt[i],AWerr_Tpc5Gt[i]);
  //trigger bias eta_pair<0      
  biaspTLt1[i]=(1-ratiopTLt1[i])*TMath::Max(A_Tpc1Lt[i],AWerr_Tpc1Lt[i]);
  biaspTLt2[i]=(1-ratiopTLt2[i])*TMath::Max(A_Tpc2Lt[i],AWerr_Tpc2Lt[i]);
  biaspTLt3[i]=(1-ratiopTLt3[i])*TMath::Max(A_Tpc3Lt[i],AWerr_Tpc3Lt[i]);
  biaspTLt4[i]=(1-ratiopTLt4[i])*TMath::Max(A_Tpc4Lt[i],AWerr_Tpc4Lt[i]);
  biaspTLt5[i]=(1-ratiopTLt5[i])*TMath::Max(A_Tpc5Lt[i],AWerr_Tpc5Lt[i]);

 rrdiff1Gt[i]=(double)(A_Tpc1Gt[i]-A_Tof1Gt[i])/A_Tpc1Gt[i];//relative difference
 rrdiff1errGt[i]=calErr(A_Tpc1Gt[i],AWerr_Tpc1Gt[i],A_Tof1Gt[i],AWerr_Tof1Gt[i]);//relative difference
 rrdiff2Gt[i]=(double)(A_Tpc2Gt[i]-A_Tof2Gt[i])/A_Tpc2Gt[i];
 rrdiff2errGt[i]=calErr(A_Tpc2Gt[i],AWerr_Tpc2Gt[i],A_Tof2Gt[i],AWerr_Tof2Gt[i]);//relative difference
 rrdiff3Gt[i]=(double)(A_Tpc3Gt[i]-A_Tof3Gt[i])/A_Tpc3Gt[i];
 rrdiff3errGt[i]=calErr(A_Tpc3Gt[i],AWerr_Tpc3Gt[i],A_Tof3Gt[i],AWerr_Tof3Gt[i]);//relative difference
 rrdiff4Gt[i]=(double)(A_Tpc4Gt[i]-A_Tof4Gt[i])/A_Tpc4Gt[i];
 rrdiff4errGt[i]=calErr(A_Tpc4Gt[i],AWerr_Tpc4Gt[i],A_Tof4Gt[i],AWerr_Tof4Gt[i]);//relative difference
 rrdiff5Gt[i]=(double)(A_Tpc5Gt[i]-A_Tof5Gt[i])/A_Tpc5Gt[i];
 rrdiff5errGt[i]=calErr(A_Tpc5Gt[i],AWerr_Tpc5Gt[i],A_Tof5Gt[i],AWerr_Tof5Gt[i]);//relative difference

 rrdiff1Lt[i]=(double)(A_Tpc1Lt[i]-A_Tof1Lt[i])/A_Tpc1Lt[i];//relative difference
 rrdiff1errLt[i]=calErr(A_Tpc1Lt[i],AWerr_Tpc1Lt[i],A_Tof1Lt[i],AWerr_Tof1Lt[i]);//relative difference
 rrdiff2Lt[i]=(double)(A_Tpc2Lt[i]-A_Tof2Lt[i])/A_Tpc2Lt[i];
 rrdiff2errLt[i]=calErr(A_Tpc2Lt[i],AWerr_Tpc2Lt[i],A_Tof2Lt[i],AWerr_Tof2Lt[i]);//relative difference
 rrdiff3Lt[i]=(double)(A_Tpc3Lt[i]-A_Tof3Lt[i])/A_Tpc3Lt[i];
 rrdiff3errLt[i]=calErr(A_Tpc3Lt[i],AWerr_Tpc3Lt[i],A_Tof3Lt[i],AWerr_Tof3Lt[i]);//relative difference
 rrdiff4Lt[i]=(double)(A_Tpc4Lt[i]-A_Tof4Lt[i])/A_Tpc4Lt[i];
 rrdiff4errLt[i]=calErr(A_Tpc4Lt[i],AWerr_Tpc4Lt[i],A_Tof4Lt[i],AWerr_Tof4Lt[i]);//relative difference
 rrdiff5Lt[i]=(double)(A_Tpc5Lt[i]-A_Tof5Lt[i])/A_Tpc5Lt[i];
 rrdiff5errLt[i]=calErr(A_Tpc5Lt[i],AWerr_Tpc5Lt[i],A_Tof5Lt[i],AWerr_Tof5Lt[i]);//relative difference






  //difference between A_UT with TPC only and TOF only
  rdiff1Gt[i]=(double)abs(A_Tpc1Gt[i]-A_Tof1Gt[i]);
  rdiff1errGt[i]=(double)sqrt(pow(AWerr_Tpc1Gt[i],2)+pow(AWerr_Tof1Gt[i],2));
  rdiff2Gt[i]=(double)abs(A_Tpc2Gt[i]-A_Tof2Gt[i]);
  rdiff2errGt[i]=(double)sqrt(pow(AWerr_Tpc2Gt[i],2)+pow(AWerr_Tof2Gt[i],2));
  rdiff3Gt[i]=(double)abs(A_Tpc3Gt[i]-A_Tof3Gt[i]);
  rdiff3errGt[i]=(double)sqrt(pow(AWerr_Tpc3Gt[i],2)+pow(AWerr_Tof3Gt[i],2));
  rdiff4Gt[i]=(double)abs(A_Tpc4Gt[i]-A_Tof4Gt[i]);
  rdiff4errGt[i]=(double)sqrt(pow(AWerr_Tpc4Gt[i],2)+pow(AWerr_Tof4Gt[i],2));
  rdiff5Gt[i]=(double)abs(A_Tpc5Gt[i]-A_Tof5Gt[i]);
  rdiff5errGt[i]=(double)sqrt(pow(AWerr_Tpc5Gt[i],2)+pow(AWerr_Tof5Gt[i],2));

  rdiff1Lt[i]=(double)abs(A_Tpc1Lt[i]-A_Tof1Lt[i]);
  rdiff1errLt[i]=(double)sqrt(pow(AWerr_Tpc1Lt[i],2)+pow(AWerr_Tof1Lt[i],2));
  rdiff2Lt[i]=(double)abs(A_Tpc2Lt[i]-A_Tof2Lt[i]);
  rdiff2errLt[i]=(double)sqrt(pow(AWerr_Tpc2Lt[i],2)+pow(AWerr_Tof2Lt[i],2));
  rdiff3Lt[i]=(double)abs(A_Tpc3Lt[i]-A_Tof3Lt[i]);
  rdiff3errLt[i]=(double)sqrt(pow(AWerr_Tpc3Lt[i],2)+pow(AWerr_Tof3Lt[i],2));
  rdiff4Lt[i]=(double)abs(A_Tpc4Lt[i]-A_Tof4Lt[i]);
  rdiff4errLt[i]=(double)sqrt(pow(AWerr_Tpc4Lt[i],2)+pow(AWerr_Tof4Lt[i],2));
  rdiff5Lt[i]=(double)abs(A_Tpc5Lt[i]-A_Tof5Lt[i]);
  rdiff5errLt[i]=(double)sqrt(pow(AWerr_Tpc5Lt[i],2)+pow(AWerr_Tof5Lt[i],2));

  pidSysGt1[i]=(1.0-pionPairMGtT[0])*abs(max(A_F1Gt[i],AWerr_F1Gt[i]));
  pidSysGt2[i]=(1.0-pionPairMGtT[1])*abs(max(A_F2Gt[i],AWerr_F2Gt[i]));
  pidSysGt3[i]=(1.0-pionPairMGtT[2])*abs(max(A_F3Gt[i],AWerr_F3Gt[i]));
  pidSysGt4[i]=(1.0-pionPairMGtT[3])*abs(max(A_F4Gt[i],AWerr_F4Gt[i]));
  pidSysGt5[i]=(1.0-pionPairMGtT[4])*abs(max(A_F5Gt[i],AWerr_F5Gt[i]));

  pidSysLt1[i]=(1.0-pionPairMLtT[0])*abs(max(A_F1Lt[i],AWerr_F1Lt[i]));
  pidSysLt2[i]=(1.0-pionPairMLtT[1])*abs(max(A_F2Lt[i],AWerr_F2Lt[i]));
  pidSysLt3[i]=(1.0-pionPairMLtT[2])*abs(max(A_F3Lt[i],AWerr_F3Lt[i]));
  pidSysLt4[i]=(1.0-pionPairMLtT[3])*abs(max(A_F4Lt[i],AWerr_F4Lt[i]));
  pidSysLt5[i]=(1.0-pionPairMLtT[4])*abs(max(A_F5Lt[i],AWerr_F5Lt[i]));


  // cout<<A_F1[i]<<" , "<<Aerr_F1[i]<<", frac: "<<pionPairPtGtT[0]<<", imp: "<< double (1-pionPairPtGtT[0])<<" A: "<<double (abs(TpTath::Max(A_F1[i],Aerr_F1[i])))<< ", Sys: "<< pidSysGt1[i]<<endl;
 }

 TGraphErrors *gr_AvspTFGt[5]; 
 TGraphErrors *gr_AvspTFLt[5]; 
 TGraphErrors *gr_AvspTgtS[5]; 
 TGraphErrors *gr_AvspTgtSc[5]; 
 TGraphErrors *gr_AvspTgtTOF[5]; 
 TGraphErrors *gr_rdiffGt[5]; 
 TGraphErrors *gr_rdiffLt[5]; 
 TGraphErrors *gr_trdiff[5]; 
 TGraphErrors *gr_rrdiffGt[5];
 TGraphErrors *gr_rrdiffLt[5];
 double pT[9]={0};

 double rdiffGt[9]={0};
 double rdiffLt[9]={0};
 double rdifferrGt[9]={0};
 double rdifferrLt[9]={0};

 double rrdiffGt[9]={0};
 double rrdiffLt[9]={0};
 double rrdifferrGt[9]={0};
 double rrdifferrLt[9]={0};

 double biaspTGt[9]={0};
 double biaspTLt[9]={0};
 TLegend *leg[5];
 TLine *line[6];
 TLatex tex[6];

 double errx[9]={0};
 double errxps[9]={0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15};
 double errxts[9]={0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15};

 double A_FGt[9]={0};
 double Aerr_FGt[9]={0};
 double A_FLt[9]={0};
 double Aerr_FLt[9]={0};

 TGraphErrors *gr_FGtS[5]={0};
 TGraphErrors *gr_FLtS[5]={0};
double frp0Gt[5]={0};
 double frp0Lt[5]={0};




 for(int pad=0;pad<5; pad++ ){ 
  //if(pad==5) break;
  if(pad==0)for(int mbin=0; mbin<9; mbin++){pT[mbin]=pT1[mbin];  A_FGt[mbin]=A_Tpc1Gt[mbin]; Aerr_FGt[mbin]=AWerr_Tpc1Gt[mbin]; A_FLt[mbin]=A_Tpc1Lt[mbin]; Aerr_FLt[mbin]=AWerr_Tpc1Lt[mbin]; rdiffGt[mbin]=rdiff1Gt[mbin]; rdifferrGt[mbin]=rdiff1errGt[mbin]; rdiffLt[mbin]=rdiff1Lt[mbin]; rdifferrLt[mbin]=rdiff1errLt[mbin];rrdiffGt[mbin]=rrdiff1Gt[mbin]; rrdifferrGt[mbin]=rrdiff1errGt[mbin]; rrdiffLt[mbin]=rrdiff1Lt[mbin]; rrdifferrLt[mbin]=rrdiff1errLt[mbin];} 
  if(pad==1)for(int mbin=0; mbin<9; mbin++){pT[mbin]=pT2[mbin];  A_FGt[mbin]=A_Tpc2Gt[mbin]; Aerr_FGt[mbin]=AWerr_Tpc2Gt[mbin]; A_FLt[mbin]=A_Tpc2Lt[mbin]; Aerr_FLt[mbin]=AWerr_Tpc2Lt[mbin]; rdiffGt[mbin]=rdiff2Gt[mbin]; rdifferrGt[mbin]=rdiff2errGt[mbin]; rdiffLt[mbin]=rdiff2Lt[mbin]; rdifferrLt[mbin]=rdiff2errLt[mbin];rrdiffGt[mbin]=rrdiff2Gt[mbin]; rrdifferrGt[mbin]=rrdiff2errGt[mbin]; rrdiffLt[mbin]=rrdiff2Lt[mbin]; rrdifferrLt[mbin]=rrdiff2errLt[mbin];} 
  if(pad==2)for(int mbin=0; mbin<9; mbin++){pT[mbin]=pT3[mbin];  A_FGt[mbin]=A_Tpc3Gt[mbin]; Aerr_FGt[mbin]=AWerr_Tpc3Gt[mbin]; A_FLt[mbin]=A_Tpc3Lt[mbin]; Aerr_FLt[mbin]=AWerr_Tpc3Lt[mbin]; rdiffGt[mbin]=rdiff3Gt[mbin]; rdifferrGt[mbin]=rdiff3errGt[mbin]; rdiffLt[mbin]=rdiff3Lt[mbin]; rdifferrLt[mbin]=rdiff3errLt[mbin];rrdiffGt[mbin]=rrdiff3Gt[mbin]; rrdifferrGt[mbin]=rrdiff3errGt[mbin]; rrdiffLt[mbin]=rrdiff3Lt[mbin]; rrdifferrLt[mbin]=rrdiff3errLt[mbin];} 
  if(pad==3)for(int mbin=0; mbin<9; mbin++){pT[mbin]=pT4[mbin];  A_FGt[mbin]=A_Tpc4Gt[mbin]; Aerr_FGt[mbin]=AWerr_Tpc4Gt[mbin]; A_FLt[mbin]=A_Tpc4Lt[mbin]; Aerr_FLt[mbin]=AWerr_Tpc4Lt[mbin]; rdiffGt[mbin]=rdiff4Gt[mbin]; rdifferrGt[mbin]=rdiff4errGt[mbin]; rdiffLt[mbin]=rdiff4Lt[mbin]; rdifferrLt[mbin]=rdiff4errLt[mbin];rrdiffGt[mbin]=rrdiff4Gt[mbin]; rrdifferrGt[mbin]=rrdiff4errGt[mbin]; rrdiffLt[mbin]=rrdiff4Lt[mbin]; rrdifferrLt[mbin]=rrdiff4errLt[mbin];} 
  if(pad==4)for(int mbin=0; mbin<9; mbin++){pT[mbin]=pT5[mbin];  A_FGt[mbin]=A_Tpc5Gt[mbin]; Aerr_FGt[mbin]=AWerr_Tpc5Gt[mbin]; A_FLt[mbin]=A_Tpc5Lt[mbin]; Aerr_FLt[mbin]=AWerr_Tpc5Lt[mbin]; rdiffGt[mbin]=rdiff5Gt[mbin]; rdifferrGt[mbin]=rdiff5errGt[mbin]; rdiffLt[mbin]=rdiff5Lt[mbin]; rdifferrLt[mbin]=rdiff5errLt[mbin];rrdiffGt[mbin]=rrdiff5Gt[mbin]; rrdifferrGt[mbin]=rrdiff5errGt[mbin]; rrdiffLt[mbin]=rrdiff5Lt[mbin]; rrdifferrLt[mbin]=rrdiff5errLt[mbin];} 
  //draw final asymmetries
  gr_AvspTFGt[pad] = new TGraphErrors(9,pT, A_FGt,errx,Aerr_FGt);//this draws asymmetry after polynomial TOF cut with all other tracks included
  gr_AvspTFGt[pad]->GetYaxis()-> SetTitle("#font[22]{A_{UT}^{Sin(#Phi_{S}-#Phi_{R})}}");
  gr_AvspTFGt[pad]->SetTitle("");
  gr_AvspTFGt[pad]->GetYaxis()->SetTitleOffset(1.2);
  gr_AvspTFGt[pad]->GetYaxis()->SetTitleSize(0.06);
  gr_AvspTFGt[pad]->GetYaxis()->SetLabelSize(0.06);
  gr_AvspTFGt[pad]->GetYaxis()->SetLabelFont(22);
  gr_AvspTFGt[pad]->GetXaxis()->SetTitle("#font[22]{p^{#pi^{+}#pi^{-}}_{T}(GeV/c)}");
  gr_AvspTFGt[pad]->GetXaxis()->SetTitleSize(0.06);
  gr_AvspTFGt[pad]->GetXaxis()->SetLabelSize(0.06);
  gr_AvspTFGt[pad]->GetXaxis()-> SetLabelFont(22);
  gr_AvspTFGt[pad]->GetXaxis()->SetLimits(2.5,11.5);
  gr_AvspTFGt[pad]->GetYaxis()->SetRangeUser(-0.015, 0.068);
  gr_AvspTFGt[pad]-> SetMarkerStyle(20);
  gr_AvspTFGt[pad]-> SetMarkerColor(2);
  gr_AvspTFGt[pad]-> SetLineColor(2);
  gr_AvspTFGt[pad]-> SetLineWidth(1);

  gr_AvspTFLt[pad] = new TGraphErrors(9,pT, A_FLt,errx,Aerr_FLt);//this draws asymmetry after polynomial TOF cut with all other tracks included
  gr_AvspTFLt[pad]-> SetMarkerStyle(20);
  gr_AvspTFLt[pad]-> SetMarkerColor(4);
  gr_AvspTFLt[pad]-> SetLineColor(4);
  gr_AvspTFLt[pad]-> SetLineWidth(1);

  gr_rdiffGt[pad] = new TGraphErrors(9,pT, rdiffGt, 0, rdifferrGt);//this draws asymmetry after polynomial TOF cut with all other tracks included
  gr_rdiffGt[pad]->GetYaxis()-> SetTitle("#font[22]{|A_{Tpc} - A_{Tof}|}");
  gr_rdiffGt[pad]->SetTitle("");
  gr_rdiffGt[pad]->GetYaxis()->SetTitleOffset(1.2);
  gr_rdiffGt[pad]->GetYaxis()->SetTitleSize(0.06);
  gr_rdiffGt[pad]->GetYaxis()->SetLabelSize(0.06);
  gr_rdiffGt[pad]->GetYaxis()->SetLabelFont(22);
  gr_rdiffGt[pad]->GetXaxis()->SetTitle("#font[22]{p^{#pi^{+}#pi^{-}}_{T}(GeV/c)}");
  gr_rdiffGt[pad]->GetXaxis()->SetTitleSize(0.06);
  gr_rdiffGt[pad]->GetXaxis()->SetLabelSize(0.06);
  gr_rdiffGt[pad]->GetXaxis()-> SetLabelFont(22);
  gr_rdiffGt[pad]->GetXaxis()->SetLimits(2.5,11.5);
  gr_rdiffGt[pad]->GetYaxis()->SetRangeUser(-0.015, 0.068);
  gr_rdiffGt[pad]-> SetMarkerStyle(20);
  gr_rdiffGt[pad]-> SetMarkerColor(2);
  gr_rdiffGt[pad]-> SetLineColor(2);
  gr_rdiffGt[pad]-> SetLineWidth(1);

  gr_rdiffLt[pad] = new TGraphErrors(9,pT, rdiffLt, 0, rdifferrLt);//this draws asymmetry after polynomial TOF cut with all other tracks included
  gr_rdiffLt[pad]-> SetMarkerStyle(20);
  gr_rdiffLt[pad]-> SetMarkerColor(4);
  gr_rdiffLt[pad]-> SetLineColor(4);
  gr_rdiffLt[pad]-> SetLineWidth(1);

gr_rrdiffGt[pad] = new TGraphErrors(9,pT, rrdiffGt, 0, rrdifferrGt);//this draws asymmetry after polynomial TOF cut with all other tracks included
gr_rrdiffGt[pad]->GetYaxis()-> SetTitle("#font[22]{A_{tpc} - A_{tof}/A_{tpc}}");
gr_rrdiffGt[pad]->SetTitle("");
gr_rrdiffGt[pad]->GetYaxis()->SetTitleOffset(1.2);
gr_rrdiffGt[pad]->GetYaxis()->SetTitleSize(0.06);
gr_rrdiffGt[pad]->GetYaxis()->SetLabelSize(0.06);
gr_rrdiffGt[pad]->GetYaxis()->SetLabelFont(22);
gr_rrdiffGt[pad]->GetXaxis()->SetTitle("#font[22]{p^{#pi^{+}#pi^{-}_{T}}(GeV/c)}");
gr_rrdiffGt[pad]->GetXaxis()->SetTitleSize(0.06);
gr_rrdiffGt[pad]->GetXaxis()->SetLabelSize(0.06);
gr_rrdiffGt[pad]->GetXaxis()-> SetLabelFont(22);
gr_rrdiffGt[pad]->GetXaxis()->SetLimits(2.5,11.5);
gr_rrdiffGt[pad]->GetYaxis()->SetRangeUser(-35., 35);
gr_rrdiffGt[pad]-> SetMarkerStyle(20);
gr_rrdiffGt[pad]-> SetMarkerColor(2);
gr_rrdiffGt[pad]-> SetLineColor(2);
gr_rrdiffGt[pad]-> SetLineWidth(1);

gr_rrdiffLt[pad] = new TGraphErrors(9,pT, rrdiffLt, 0, rrdifferrLt);//this draws asymmetry after polynomial TOF cut with all other tracks included
gr_rrdiffLt[pad]-> SetMarkerStyle(20);
gr_rrdiffLt[pad]-> SetMarkerColor(4);
gr_rrdiffLt[pad]-> SetLineColor(4);
gr_rrdiffLt[pad]-> SetLineWidth(1);
 }
TCanvas *crdiff=new TCanvas("crdiff","",150,10,1100,700);
crdiff->Divide(3,2);
for(int pad=0; pad<6; pad++){
 crdiff->cd(pad+1);
 setPad(pad);
 if(pad==5){
  gPad->SetGrid(0,0);
  tex[pad].DrawLatex(0.1, 0.7, "#font[22]{#color[2]{#eta^{#pi^{+}#pi^{-}} #GT 0 }}");
  tex[pad].DrawLatex(0.1, 0.6, "#font[22]{#color[4]{#eta^{#pi^{+}#pi^{-}} #LT 0 }}");
  tex[pad].DrawLatex(0.1, 0.5, "#font[22]{A_{1} = A_{UT}(-1<n#sigma_{#pi}<2) }");
  tex[pad].DrawLatex(0.1, 0.4, "#font[22]{A_{2} = A_{UT}(-1<n#sigma_{#pi}<2) && Pol. TOF cut}");

 }else if(pad<5){
  gPad->SetGrid(0,0);
  gr_rrdiffGt[pad]->Draw("AP");
  if(pad==1 || pad==2 || pad==4)gr_rrdiffGt[pad]->GetYaxis()->SetLabelSize(0);
  gr_rrdiffGt[pad]->GetYaxis()->SetNdivisions(505);
  gr_rrdiffGt[pad]->GetXaxis()->SetNdivisions(505);

  gr_rrdiffGt[pad]->Fit("pol0");
  frp0Gt[pad]=(double)gr_rrdiffGt[pad]->GetFunction("pol0")->GetParameter(0);
  gr_rrdiffLt[pad]->Draw("SAME P");
  gr_rrdiffLt[pad]->Fit("pol0");
  frp0Lt[pad]=(double)gr_rrdiffLt[pad]->GetFunction("pol0")->GetParameter(0);
  gr_rrdiffLt[pad]->GetFunction("pol0")->SetLineColor(4);
  cout<<frp0Gt[pad]<<", "<<frp0Lt[pad]<<endl;

  gPad->Update();

  line[pad]=  new TLine(crdiff->cd(pad+1)->GetUxmin(),0.,crdiff->cd(pad+1)->GetUxmax(),0.);
  line[pad]->SetLineStyle(2);
  line[pad]->SetLineWidth(1);
  line[pad]->Draw();


  tex[pad].SetTextSize(0.04);
  tex[pad].SetTextAlign(13);
  tex[pad].DrawLatex(4,34,Form("#font[22]{#color[2]{chi2/ndf  = %g / %i}}",gr_rrdiffGt[pad]->GetFunction("pol0")->GetChisquare(),gr_rrdiffGt[pad]->GetFunction("pol0")->GetNDF()));
  tex[pad].DrawLatex(4,30,Form("#font[22]{#color[2]{p0  = %g #pm %g}}",gr_rrdiffGt[pad]->GetFunction("pol0")->GetParameter(0),gr_rrdiffGt[pad]->GetFunction("pol0")->GetParError(0)));

  tex[pad].DrawLatex(4,25,Form("#font[22]{#color[4]{chi2/ndf  = %g / %i}}",gr_rrdiffLt[pad]->GetFunction("pol0")->GetChisquare(),gr_rrdiffLt[pad]->GetFunction("pol0")->GetNDF()));
  tex[pad].DrawLatex(4,21,Form("#font[22]{#color[4]{p0  = %g #pm %g}}",gr_rrdiffLt[pad]->GetFunction("pol0")->GetParameter(0),gr_rrdiffLt[pad]->GetFunction("pol0")->GetParError(0)));

  tex[pad].DrawLatex(8,34,Form("#font[22]{#color[1]{< M_{inv} > = %g GeV/c}}",avgM[pad]));

  gPad->Update();
 }
}
crdiff->Update();
crdiff->SaveAs("./Plots/pidsysRelativeDiffWithTpc.pdf");






 TCanvas *cdiff=new TCanvas("cdiff","",150,10,1100,700);
 cdiff->Divide(3,2);
 for(int pad=0; pad<6; pad++){
  cdiff->cd(pad+1);
  setPad(pad);
  if(pad==5){
   gPad->SetGrid(0,0);
   tex[pad].DrawLatex(0.1, 0.7, "#font[22]{#color[2]{#eta^{#pi^{+}#pi^{-}} #GT 0 }}");
   tex[pad].DrawLatex(0.1, 0.6, "#font[22]{#color[4]{#eta^{#pi^{+}#pi^{-}} #LT 0 }}");
   tex[pad].DrawLatex(0.1, 0.5, "#font[22]{A_{1} = A_{UT}(-1<n#sigma_{#pi}<2) }");
   tex[pad].DrawLatex(0.1, 0.4, "#font[22]{A_{2} = A_{UT}(-1<n#sigma_{#pi}<2) && Pol. TOF cut}");

  }else if(pad<5){
   gPad->SetGrid(0,0);
   gr_rrdiffGt[pad]->Draw("AP");
   if(pad==1 || pad==2 || pad==4)gr_rrdiffGt[pad]->GetYaxis()->SetLabelSize(0);
   gr_rrdiffGt[pad]->GetYaxis()->SetNdivisions(505);
   gr_rrdiffGt[pad]->GetXaxis()->SetNdivisions(505);

   gr_rrdiffGt[pad]->Fit("pol0");
   fp0Gt[pad]=(double)gr_rrdiffGt[pad]->GetFunction("pol0")->GetParameter(0);
   gr_rrdiffLt[pad]->Draw("SAME P");
   gr_rrdiffLt[pad]->Fit("pol0");
   fp0Lt[pad]=(double)gr_rrdiffLt[pad]->GetFunction("pol0")->GetParameter(0);
   gr_rrdiffLt[pad]->GetFunction("pol0")->SetLineColor(4);
   cout<<fp0Gt[pad]<<", "<<fp0Lt[pad]<<endl;

   gPad->Update();

   line[pad]=  new TLine(cdiff->cd(pad+1)->GetUxmin(),0.,cdiff->cd(pad+1)->GetUxmax(),0.);
   line[pad]->SetLineStyle(2);
   line[pad]->SetLineWidth(1);
   line[pad]->Draw();


   tex[pad].SetTextSize(0.04);
   tex[pad].SetTextAlign(13);
   tex[pad].DrawLatex(3,33,Form("#font[22]{#color[2]{chi2/ndf  = %g #pm %i}}",gr_rrdiffGt[pad]->GetFunction("pol0")->GetChisquare(),gr_rrdiffGt[pad]->GetFunction("pol0")->GetNDF()));
   tex[pad].DrawLatex(3,29,Form("#font[22]{#color[2]{p0  = %g #pm %g}}",gr_rrdiffGt[pad]->GetFunction("pol0")->GetParameter(0),gr_rrdiffGt[pad]->GetFunction("pol0")->GetParError(0)));

   tex[pad].DrawLatex(3.0,25,Form("#font[22]{#color[4]{chi2/ndf  = %g #pm %i}}",gr_rrdiffLt[pad]->GetFunction("pol0")->GetChisquare(),gr_rrdiffLt[pad]->GetFunction("pol0")->GetNDF()));
   tex[pad].DrawLatex(3.0,21,Form("#font[22]{#color[4]{p0  = %g #pm %g}}",gr_rrdiffLt[pad]->GetFunction("pol0")->GetParameter(0),gr_rrdiffLt[pad]->GetFunction("pol0")->GetParError(0)));

   tex[pad].DrawLatex(8.0,33,Form("#font[22]{#color[1]{< M_{inv} > = %g GeV/c^{2}}}",avgM[pad]));

   gPad->Update();
  }
 }
 cdiff->Update();
 cdiff->SaveAs("./Plots/pidsysRatioPlot_pT.pdf");
 //final systematic error
 for(int i=0; i<9; i++){
  fpidSysGt1[i]=sqrt(pow(fp0Gt[0],2)+pow(pidSysGt1[0],2)); 
  fpidSysGt2[i]=sqrt(pow(fp0Gt[1],2)+pow(pidSysGt2[1],2)); 
  fpidSysGt3[i]=sqrt(pow(fp0Gt[2],2)+pow(pidSysGt3[2],2)); 
  fpidSysGt4[i]=sqrt(pow(fp0Gt[3],2)+pow(pidSysGt4[3],2)); 
  fpidSysGt5[i]=sqrt(pow(fp0Gt[4],2)+pow(pidSysGt5[4],2)); 

  fpidSysLt1[i]=sqrt(pow(fp0Lt[0],2)+pow(pidSysLt1[0],2)); 
  fpidSysLt2[i]=sqrt(pow(fp0Lt[1],2)+pow(pidSysLt2[1],2)); 
  fpidSysLt3[i]=sqrt(pow(fp0Lt[2],2)+pow(pidSysLt3[2],2)); 
  fpidSysLt4[i]=sqrt(pow(fp0Lt[3],2)+pow(pidSysLt4[3],2)); 
  fpidSysLt5[i]=sqrt(pow(fp0Lt[4],2)+pow(pidSysLt5[4],2)); 

//final systematics. Combined trigger bias and PID as quadrature sum
 combSyspTGt1[i]=sqrt(pow(biaspTGt1[i],2) + pow(frp0Gt[0]*A_Tpc1Gt[i],2));
 combSyspTGt2[i]=sqrt(pow(biaspTGt2[i],2) + pow(frp0Gt[1]*A_Tpc2Gt[i],2));
 combSyspTGt3[i]=sqrt(pow(biaspTGt3[i],2) + pow(frp0Gt[2]*A_Tpc3Gt[i],2));
 combSyspTGt4[i]=sqrt(pow(biaspTGt4[i],2) + pow(frp0Gt[3]*A_Tpc4Gt[i],2));
 combSyspTGt5[i]=sqrt(pow(biaspTGt5[i],2) + pow(frp0Gt[4]*A_Tpc5Gt[i],2));
 //cout<<"trig bias: "<<biaspTGt1[i]<<", pid: "<<Gt1[i]<<",  combined: "<<combSyspTGt1[i]<<endl;
 combSyspTLt1[i]=sqrt(pow(biaspTLt1[i],2) + pow(frp0Lt[0]*A_Tpc1Lt[i],2));
 combSyspTLt2[i]=sqrt(pow(biaspTLt2[i],2) + pow(frp0Lt[1]*A_Tpc2Lt[i],2));
 combSyspTLt3[i]=sqrt(pow(biaspTLt3[i],2) + pow(frp0Lt[2]*A_Tpc3Lt[i],2));
 combSyspTLt4[i]=sqrt(pow(biaspTLt4[i],2) + pow(frp0Lt[3]*A_Tpc4Lt[i],2));
 combSyspTLt5[i]=sqrt(pow(biaspTLt5[i],2) + pow(frp0Lt[4]*A_Tpc5Lt[i],2));

/*  //final systematics. Combined trigger bias and PID as quadrature sum
  combSyspTGt1[i]=sqrt(pow(biaspTGt1[i],2) + pow(fpidSysGt1[i],2));
  combSyspTGt2[i]=sqrt(pow(biaspTGt2[i],2) + pow(fpidSysGt2[i],2));
  combSyspTGt3[i]=sqrt(pow(biaspTGt3[i],2) + pow(fpidSysGt3[i],2));
  combSyspTGt4[i]=sqrt(pow(biaspTGt4[i],2) + pow(fpidSysGt4[i],2));
  combSyspTGt5[i]=sqrt(pow(biaspTGt5[i],2) + pow(fpidSysGt5[i],2));
  //cout<<"trig bias: "<<biaspTGt1[i]<<", pid: "<<fpidSysGt1[i]<<",  combined: "<<combSyspTGt1[i]<<endl;
  combSyspTLt1[i]=sqrt(pow(biaspTLt1[i],2) + pow(fpidSysLt1[i],2));
  combSyspTLt2[i]=sqrt(pow(biaspTLt2[i],2) + pow(fpidSysLt2[i],2));
  combSyspTLt3[i]=sqrt(pow(biaspTLt3[i],2) + pow(fpidSysLt3[i],2));
  combSyspTLt4[i]=sqrt(pow(biaspTLt4[i],2) + pow(fpidSysLt4[i],2));
  combSyspTLt5[i]=sqrt(pow(biaspTLt5[i],2) + pow(fpidSysLt5[i],2));
*/
 }

 double finalSysGt[9]={0};
 double finalSysLt[9]={0};

 TCanvas *cfinal=new TCanvas("cfinal","",150,10,1100,700);
 cfinal->Divide(3,2,0,0);

 for(int pad=0;pad<6; pad++ ){
  if(pad==5){
   cfinal->cd(pad+1);
   setPad(pad);

  }else{
   if(pad==0)for(int mbin=0; mbin<9; mbin++){ pT[mbin]=pT1[mbin];  A_FGt[mbin]=A_Tpc1Gt[mbin]; Aerr_FGt[mbin]=AWerr_Tpc1Gt[mbin]; A_FLt[mbin]=A_Tpc1Lt[mbin]; Aerr_FLt[mbin]=AWerr_Tpc1Lt[mbin]; finalSysGt[mbin]=combSyspTGt1[mbin]; finalSysLt[mbin]=combSyspTLt1[mbin]; }
   if(pad==1)for(int mbin=0; mbin<9; mbin++){pT[mbin]=pT2[mbin];  A_FGt[mbin]=A_Tpc2Gt[mbin]; Aerr_FGt[mbin]=AWerr_Tpc2Gt[mbin]; A_FLt[mbin]=A_Tpc2Lt[mbin]; Aerr_FLt[mbin]=AWerr_Tpc2Lt[mbin];  finalSysGt[mbin]=combSyspTGt2[mbin]; finalSysLt[mbin]=combSyspTLt2[mbin]; }
   if(pad==2)for(int mbin=0; mbin<9; mbin++){pT[mbin]=pT3[mbin];  A_FGt[mbin]=A_Tpc3Gt[mbin]; Aerr_FGt[mbin]=AWerr_Tpc3Gt[mbin]; A_FLt[mbin]=A_Tpc3Lt[mbin]; Aerr_FLt[mbin]=AWerr_Tpc3Lt[mbin];  finalSysGt[mbin]=combSyspTGt3[mbin]; finalSysLt[mbin]=combSyspTLt3[mbin]; }
   if(pad==3)for(int mbin=0; mbin<9; mbin++){pT[mbin]=pT4[mbin];  A_FGt[mbin]=A_Tpc4Gt[mbin]; Aerr_FGt[mbin]=AWerr_Tpc4Gt[mbin]; A_FLt[mbin]=A_Tpc4Lt[mbin]; Aerr_FLt[mbin]=AWerr_Tpc4Lt[mbin];  finalSysGt[mbin]=combSyspTGt4[mbin]; finalSysLt[mbin]=combSyspTLt4[mbin];  }
   if(pad==4)for(int mbin=0; mbin<9; mbin++){pT[mbin]=pT5[mbin];  A_FGt[mbin]=A_Tpc5Gt[mbin]; Aerr_FGt[mbin]=AWerr_Tpc5Gt[mbin]; A_FLt[mbin]=A_Tpc5Lt[mbin]; Aerr_FLt[mbin]=AWerr_Tpc5Lt[mbin];  finalSysGt[mbin]=combSyspTGt5[mbin]; finalSysLt[mbin]=combSyspTLt5[mbin];  }

   gr_FGtS[pad] = new TGraphErrors(9,pT,A_FGt,errxts,finalSysGt);
   gr_FGtS[pad]->GetYaxis()-> SetTitle("#font[22]{A_{UT}^{Sin(#Phi_{S}-#Phi_{R})}} ");
   gr_FGtS[pad]->SetTitle("");
   gr_FGtS[pad]->GetYaxis()->SetTitleOffset(1.2);
   gr_FGtS[pad]->GetYaxis()->SetTitleSize(0.06);
   gr_FGtS[pad]->GetYaxis()->SetLabelSize(0.06);
   gr_FGtS[pad]->GetYaxis()->SetLabelFont(22);
   gr_FGtS[pad]->GetXaxis()->SetTitle("#font[22]{p^{#pi^{+}#pi^{-}}_{T} (GeV/c)}");
   gr_FGtS[pad]->GetXaxis()->SetTitleSize(0.06);
   gr_FGtS[pad]->GetXaxis()->SetLabelSize(0.06);
   gr_FGtS[pad]->GetXaxis()-> SetLabelFont(22);
   gr_FGtS[pad]->GetXaxis()->SetLimits(2.5,11.5);
   gr_FGtS[pad]->GetYaxis()->SetRangeUser(-0.015, 0.068);
   gr_FGtS[pad]->SetFillStyle(0000);
   gr_FGtS[pad]->SetLineColor(1);
   gr_FGtS[pad]->SetLineWidth(2);


   gr_FLtS[pad] = new TGraphErrors(9,pT,A_FLt,errxts,finalSysLt);
   gr_FLtS[pad]->SetFillStyle(0000);
   gr_FLtS[pad]->SetLineColor(1);
   gr_FLtS[pad]->SetLineWidth(2);

   cfinal->cd(pad+1);
   setPad(pad);


   gr_FGtS[pad]->Draw("AP2");
   if(pad==1 ||pad==2 || pad==4){
    gr_FGtS[pad]->GetYaxis()->SetLabelSize(0);
    gr_FGtS[pad]->GetYaxis()->SetTitle("");
   }
   gr_FGtS[pad]->GetYaxis()->SetNdivisions(505);
   gr_FGtS[pad]->GetXaxis()->SetNdivisions(505);
   gr_FLtS[pad]->Draw("SAME 2");
   gr_AvspTFGt[pad]->Draw("SAME P");
   gr_AvspTFLt[pad]->Draw("SAME P");

   gPad->Update();
   line[pad]=  new TLine(cfinal->cd(pad+1)->GetUxmin(),0.,cfinal->cd(pad+1)->GetUxmax(),0.);
   line[pad]->SetLineStyle(2);
   line[pad]->SetLineWidth(1);
   line[pad]->Draw();


   tex[pad].SetTextSize(0.06);
   tex[pad].SetTextAlign(13);
   tex[pad].DrawLatex(3.0,0.062,Form("#font[22]{#color[1]{< M_{inv} > = %g GeV/c^{2}}}",avgM[pad])); 

   if(pad==0){
    leg[pad] = new TLegend(0.7,.61, 0.97, 0.81);
    leg[pad]->AddEntry(gr_AvspTFGt[pad], " #font[22]{ #eta^{#pi^{+}#pi^{-}} > 0}", "lp");
    leg[pad]->AddEntry(gr_AvspTFLt[pad], " #font[22]{ #eta^{#pi^{+}#pi^{-}} < 0}", "lp");
    leg[pad]->AddEntry(gr_FGtS[pad], " #font[22]{ Syst. Error}", "f");
    leg[pad]->SetTextSize(0.05);
    leg[pad]->Draw();
   }



  }



 }
 cfinal->Update();
 cfinal->SaveAs("./Plots/finalAut_pT.pdf");

}

//error calculation
double calErr(double a, double ea, double b, double eb){
 double err=0;
 err=sqrt(pow(((b/a)*ea),2)+pow(eb/a,2));
 return err;
}


 void setPad(int i){

   if(i==2){
   gPad->SetTopMargin(0.137);
   gPad->SetLeftMargin(0.02);
   gPad->SetBottomMargin(0.15);
   gPad->SetRightMargin(0.02);
   gPad->SetPad(.6634,0.4197,0.97,.9897);
   gPad->SetFillStyle(4000);
   gPad->Update();
   }else if(i==5){
   gPad->SetTopMargin(0.4);
   gPad->SetLeftMargin(0.01);
   gPad->SetBottomMargin(0.15);
   gPad->SetRightMargin(0.02);
   gPad->SetPad(.6634,0.01,0.97,.49);
   gPad->SetFillStyle(4000);
   gPad->SetFrameFillColor(0);
   gPad->SetFrameFillStyle(0);
   gPad->SetFrameBorderMode(0);

   gPad->Update();
   }else if(i==1){
  gPad->SetTopMargin(0.16);
  gPad->SetLeftMargin(0.0);
  gPad->SetBottomMargin(0.01);
  gPad->SetRightMargin(0.0);
  gPad->SetPad(.37,0.5,0.67,.99);
  gPad->Update();
  }else if(i==0){
  gPad->SetTopMargin(0.16);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.01);
  gPad->SetRightMargin(0.0);
  gPad->SetPad(0.,0.5,0.37,.99);
  }else if(i==3){
  gPad->SetTopMargin(0.0);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.15);
  gPad->SetRightMargin(0.0);
  gPad->SetPad(0.0,0.02052,0.37,.5052);
  }else if(i==4){
  gPad->SetTopMargin(0.0);
  gPad->SetLeftMargin(0.0);
  gPad->SetBottomMargin(0.15);
  gPad->SetRightMargin(0.0035);
  gPad->SetPad(.37,0.02052,0.67,.5052);
  }
}


