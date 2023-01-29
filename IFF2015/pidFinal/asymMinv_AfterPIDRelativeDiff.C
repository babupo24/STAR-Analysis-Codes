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
void  asymMinv_AfterPIDRelativeDiff(){
 gStyle -> SetOptStat(0);
 //gStyle -> SetOptFit(0);
 //gStyle -> SetLegendBorderSize(0);
 //-1<nSigmaPion<2
 double  M1[9]= {0.338101, 0.434923, 0.540022, 0.62909, 0.718894, 0.807297, 0.911696, 1.03117, 1.18508};
 double  M2[9]= {0.339518, 0.435199, 0.539328, 0.631967, 0.726091, 0.822506, 0.946478, 1.10889, 1.32425};
 double  M3[9]= {0.345042, 0.443618, 0.553353, 0.649368, 0.745406, 0.84824, 0.987234, 1.18438, 1.47069};
 double  M4[9]= {0.353874, 0.461777, 0.573937, 0.675656, 0.771731, 0.883841, 1.038, 1.26821, 1.66068};
 double  M5[9]= {0.373609, 0.503572, 0.621629, 0.7279, 0.8248, 0.95165, 1.14034, 1.41618, 2.03475};

 double avgPt[5]={3.4, 4.0, 4.7, 5.7, 8.3};
 // with start less tof 
 // polynomial TOF cut applied (F-inal result)  
 //<<<<<<<<<<<       Eta < 0, Weighted Average A_{UT}   >>>>>>>>>>>>>>>>>
 double  A_F1Lt[9] ={0.00489823, -0.00430312, -0.00356112, 0.00528712, 0.000337676, 0.000384846, 0.00558977, -0.0048005, 0.000888098};
 double  AWerr_F1Lt[9]={0.00368883, 0.00368594, 0.00371496, 0.00374052, 0.00372239, 0.00377178, 0.00386506, 0.00394903, 0.00399618};
 double  A_F2Lt[9]={0.00099076, 0.00537243, 0.00405629, -0.000573236, 0.00550309, -0.00118831, -0.00630739, 0.00419717, 0.00815201};
 double  AWerr_F2Lt[9] ={0.00355515, 0.00355438, 0.0035561, 0.00356225, 0.00354478, 0.0035875, 0.0036438, 0.00374183, 0.0038106};
 double A_F3Lt[9] ={-0.0024375, 0.00448947, -0.00286706, 0.00513057, 0.0082413, 0.0138714, 0.00347054, -0.00161353, 0.00212526};
 double AWerr_F3Lt[9] ={0.00347165, 0.0034547, 0.00346721, 0.00344753, 0.00342828, 0.00347836, 0.00351123, 0.00358481, 0.00368574};
 double  A_F4Lt[9] ={-0.00141132, 0.00242315, 0.00507164, -0.00191126, 0.000668047, 0.00796082, 0.00161702, 0.00287624, 0.00150298};
 double AWerr_F4Lt[9] ={0.00341729, 0.00340919, 0.00336966, 0.00337142, 0.00332125, 0.00338701, 0.00342922, 0.00349393, 0.00358027};
 double  A_F5Lt[9] ={-0.000545041, 0.00816163, 0.00270653, 0.00255125, 0.0128721, 0.00817112, 0.00817139, 0.00651971, 0.00572501};
 double AWerr_F5Lt[9] ={0.00338253, 0.00336336, 0.00329832, 0.0032906, 0.00330443, 0.00330024, 0.00335064, 0.00341485, 0.00350231};
 //<<<<<<<<<<<       Eta > 0, Average A_{UT}   >>>>>>>>>>>>>>>>>
 double  A_F1Gt[9] ={0.00518558, 0.00146925, 0.00825373, 0.00420507, 0.00561308, 0.0095551, 0.0157253, 0.00398878, -0.00107999};
 double  AWerr_F1Gt[9] ={0.00368517, 0.00369548, 0.00371671, 0.00375544, 0.00372459, 0.00378217, 0.0038757, 0.00395504, 0.00400719};
 double  A_F2Gt[9] ={0.00675689, 0.000471164, 0.0127479, 0.0153323, 0.0157607, 0.011665, 0.00209254, 0.00177399, 0.00163639};
 double  AWerr_F2Gt[9] ={0.00355905, 0.00355154, 0.00357443, 0.00357947, 0.00355633, 0.00359372, 0.00365525, 0.00374649, 0.00381887};
 double  A_F3Gt[9] ={0.00148311, 0.0147317, 0.0132459, 0.0161542, 0.0174046, 0.0184481, 0.0156865, 0.00798066, -0.00436376};
 double  AWerr_F3Gt[9] ={0.00347808, 0.00346723, 0.00346675, 0.00346084, 0.00344748, 0.0034848, 0.00353194, 0.0035916, 0.00370114};
 double  A_F4Gt[9] ={0.0105824, -0.000321865, 0.0143367, 0.0200754, 0.0249354, 0.0286684, 0.0165614, 0.00659794, 0.000170155};
 double  AWerr_F4Gt[9] ={0.00342492, 0.00340405, 0.00338378, 0.00339427, 0.00335975, 0.00344295, 0.00344854, 0.00349915, 0.00358906};
 double  A_F5Gt[9] ={0.0061358, 0.0155081, 0.0266134, 0.0384272, 0.0491812, 0.0486391, 0.0399287, 0.022047, 0.0217442};
 double  AWerr_F5Gt[9] ={0.00337906, 0.00337159, 0.00334915, 0.00338252, 0.00341918, 0.00343078, 0.00345388, 0.00342654, 0.00352053};
 //------Final result<<<


 //<<<<<<<<<<r  TPC only     Eta < 0, Weighted Average A_{UT}   >>>>>>>>>>>>>>>>>
 double  A_Tpc1Lt[9] ={0.00226636, -0.00459173, -0.00268284, 0.00206945, 0.000474561, -0.000288791, 0.00484177, -0.00576341, 0.00157324};
 double  AWerr_Tpc1Lt[9]={0.00349819, 0.00349526, 0.00350758, 0.00351971, 0.00350203, 0.00353291, 0.00358688, 0.00363153, 0.00367385};
 double  A_Tpc2Lt[9]={-0.00107499, 0.00520563, 0.00466177, -0.000681774, 0.00419886, -0.00132757, -0.00572289, 0.00485792, 0.00584329};
 double  AWerr_Tpc2Lt[9] ={0.0034348, 0.00343865, 0.0034307, 0.00342935, 0.00341854, 0.0034452, 0.0034745, 0.00353498, 0.00359222};
 double A_Tpc3Lt[9] ={-0.00160088, 0.00368117, -0.00188615, 0.00550266, 0.00832736, 0.0128054, 0.00348407, -0.00215124, 0.00235355};
 double AWerr_Tpc3Lt[9] ={0.00339838, 0.00338268, 0.00338999, 0.00336865, 0.0033555, 0.00339131, 0.00340384, 0.00346038, 0.00353963};
 double  A_Tpc4Lt[9] ={-0.00164831, 0.00262459, 0.00497436, -0.00121943, 0.00124444, 0.00816545, 0.000841496, 0.00398923, 0.0024252};
 double AWerr_Tpc4Lt[9] ={0.00337095, 0.00336458, 0.0033222, 0.0033232, 0.00327572, 0.00333165, 0.00336227, 0.00341328, 0.00347884};
 double  A_Tpc5Lt[9] ={-0.00081091, 0.00867389, 0.00260866, 0.00256634, 0.0130571, 0.00804663, 0.00771658, 0.00663896, 0.00497128};
 double AWerr_Tpc5Lt[9] ={0.00335129, 0.00333338, 0.00326706, 0.00326008, 0.00327395, 0.00326402, 0.00330672, 0.00336186, 0.00342742};
 //<<<<<<<<<<<  TPConly     Eta > 0, Average A_{UT}   >>>>>>>>>>>>>>>>>
 double  A_Tpc1Gt[9] ={0.00610489, 0.00265511, 0.00703812, 0.00419068, 0.00575031, 0.00722378, 0.0132474, 0.00456309, -0.00315424};
 double  AWerr_Tpc1Gt[9] ={0.00349595, 0.00350191, 0.00350922, 0.00353207, 0.0035053, 0.00354418, 0.00359277, 0.00363364, 0.00367802};
 double  A_Tpc2Gt[9] ={0.00678043, 0.00120099, 0.0101624, 0.0139553, 0.0152412, 0.00974799, 0.00227968, 0.000874331, 0.00109361};
 double  AWerr_Tpc2Gt[9] ={0.00343829, 0.003435, 0.0034443, 0.00344742, 0.00343057, 0.00344979, 0.00348286, 0.00353989, 0.00360025};
 double  A_Tpc3Gt[9] ={0.0011622, 0.0152358, 0.011955, 0.0157361, 0.0176004, 0.0176008, 0.0153222, 0.0063462, -0.00326211};
 double  AWerr_Tpc3Gt[9] ={0.0034039, 0.00339881, 0.00338927, 0.00338237, 0.00337619, 0.00339817, 0.00342631, 0.00346572, 0.0035534};
 double  A_Tpc4Gt[9] ={0.0100371, -0.000112963, 0.0141281, 0.0189976, 0.0241479, 0.0286815, 0.0176132, 0.00764366, 0.00175196};
 double  AWerr_Tpc4Gt[9] ={0.00338025, 0.00335924, 0.00333704, 0.00334438, 0.00331288, 0.00338736, 0.0033854, 0.00341885, 0.00348958};
 double  A_Tpc5Gt[9] ={0.00630835, 0.015704, 0.0264532, 0.038504, 0.0496949, 0.0483116, 0.039685, 0.0224035, 0.0200129};
 double  AWerr_Tpc5Gt[9] ={0.00334791, 0.00334174, 0.00331658, 0.00335106, 0.0033905, 0.00339471, 0.00340928, 0.00337489, 0.00344997};

 //Preliminary result eta_pair>0 
 double A_Mgt1[9]   = {0.0086015, 0.00383147, 0.00528821, 0.00161433, 0.00545937, 0.00776295, 0.0134748, 0.00365541, 0.000953803};
 double Aerr_Mgt1[9]={0.00345409, 0.00346363, 0.00346483, 0.00349142, 0.00347018, 0.00350349, 0.00355958, 0.00359207, 0.00363856};
 double A_Mgt2[9]   ={0.00723479, -0.000339341, 0.0115569, 0.0125207, 0.0144476, 0.0135844, 0.00466871, 0.000840637, 0.00140434};
 double Aerr_Mgt2[9]={0.00339975, 0.00339972, 0.00341085, 0.00340838, 0.00339196, 0.00341634, 0.00344911, 0.00350431, 0.00356671};
 double A_Mgt3[9]   ={0.00265034, 0.0138686, 0.0130309, 0.0139354, 0.0153492, 0.0156167, 0.0171356, 0.0050152, -0.00253808};
 double Aerr_Mgt3[9]={0.00336237, 0.00336156, 0.00335115, 0.00334333, 0.00333155, 0.0033593, 0.00339086, 0.00342809, 0.00351591};
 double A_Mgt4[9]   ={0.00740953, -0.000793648, 0.0119921, 0.0167145, 0.0251308, 0.0289349, 0.0192252, 0.00679528, 0.00354698};
 double Aerr_Mgt4[9]={0.00334129, 0.00331971, 0.00329577, 0.0033025, 0.00328059, 0.00335543, 0.00335095, 0.00338351, 0.00345398};
 double A_Mgt5[9]   ={0.00763535, 0.0157211, 0.0254242, 0.0412055, 0.0473758, 0.0490571, 0.0380578, 0.0230834, 0.0210097};
 double Aerr_Mgt5[9]={0.00331411, 0.00330394, 0.00327693, 0.00332278, 0.00334376, 0.00336544, 0.00336635, 0.00333871, 0.00341542};

 // TOF only eta>0
 //<<<<<<<<<<<       Eta < 0, Weighted Average A_{UT}   >>>>>>>>>>>>>>>>>
 double  A_Tof1Lt[9] ={-0.000919145, -0.00830345, -0.00411713, 0.00103787, -0.00130816, 0.000410677, 0.00110655, -0.00819202, 0.00268105};
 double  AWerr_Tof1Lt[9]={0.00560996, 0.00565905, 0.00564939, 0.00569758, 0.00565342, 0.00570627, 0.00580905, 0.00592825, 0.00599257};
 double  A_Tof2Lt[9]={-0.00528504, 0.00761524, 0.0130568, -0.000926411, 0.000465328, -0.00669712, -0.0116902, 0.00429111, 0.00502006};
 double  AWerr_Tof2Lt[9] ={0.00547525, 0.00549645, 0.00547845, 0.00547476, 0.00546357, 0.00549235, 0.00555401, 0.00569012, 0.00580166};
 double A_Tof3Lt[9] ={-0.00258572, 0.00344853, -0.00204231, 0.00454596, 0.00445146, 0.0111637, 0.00172159, -0.00224111, -0.00488648};
 double AWerr_Tof3Lt[9] ={0.00535524, 0.00538836, 0.0053373, 0.00531031, 0.00530723, 0.00535804, 0.00538749, 0.00550718, 0.00564995};
 double  A_Tof4Lt[9] ={-0.00505455, -0.00214218, 0.00622171, -0.0025953, 0.000423149, 0.00117565, 0.00564727, 0.00271607, -0.0025862};
 double AWerr_Tof4Lt[9] ={0.0053118, 0.0053268, 0.00521997, 0.00521013, 0.0051532, 0.00524236, 0.00529411, 0.00539485, 0.00553587};
 double  A_Tof5Lt[9] ={0.00322494, 0.0109131, -0.0105239, 0.000960933, 0.0164255, 0.00770273, 0.00786609, 0.00142949, 0.00351508};
 double AWerr_Tof5Lt[9] ={0.00527003, 0.00527966, 0.00513624, 0.00510263, 0.00511923, 0.00513083, 0.00520193, 0.0052091, 0.0053391};
 //<<<<<<<<<<<       Eta > 0, Average A_{UT}   >>>>>>>>>>>>>>>>>
 double  A_Tof1Gt[9] ={0.0060482, 0.00203195, 0.00460549, 0.0100715, 0.00796831, 0.00116065, 0.0134777, 0.00286132, -0.0130274};
 double  AWerr_Tof1Gt[9] ={0.00560291, 0.00565826, 0.00563217, 0.00570939, 0.00565277, 0.00573308, 0.00580973, 0.00593438, 0.00601517};
 double  A_Tof2Gt[9] ={0.000469872, 0.00390591, 0.00708019, 0.0124322, 0.0171216, 0.0163397, -0.000881484, 0.000987981, 0.00379593};
 double  AWerr_Tof2Gt[9] ={0.0054615, 0.00549126, 0.00548408, 0.00548028, 0.00547549, 0.00550828, 0.0055617, 0.00570737, 0.00580659};
 double  A_Tof3Gt[9] ={-0.00723279, 0.0202064, 0.0088815, 0.0168209, 0.0176925, 0.00849573, 0.0109812, 0.00315995, 0.00171479};
 double  AWerr_Tof3Gt[9] ={0.00535547, 0.00539899, 0.00533798, 0.00532407, 0.00531889, 0.005345, 0.0054068, 0.00553291, 0.00567632};
 double  A_Tof4Gt[9] ={0.0100551, -0.00352617, 0.0145973, 0.0170954, 0.0210753, 0.0266076, 0.0140168, 0.0112816, -0.000348832};
 double  AWerr_Tof4Gt[9] ={0.00530421, 0.0053188, 0.00523719, 0.00521606, 0.00515647, 0.00527679, 0.00530765, 0.00540792, 0.00554633};
 double  A_Tof5Gt[9] ={0.00790273, 0.0139246, 0.0238977, 0.0385186, 0.0485609, 0.0419761, 0.0424984, 0.0245977, 0.0186763};
 double  AWerr_Tof5Gt[9] ={0.00526907, 0.00528719, 0.00516858, 0.00515849, 0.00519463, 0.00518767, 0.00526863, 0.0053092, 0.00544972};


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
 double  ratioMGt1[9]={0.9006,0.9006,0.9006,0.9006,0.9006,0.9006,0.9006,0.9006,0.9006};
 double  ratioMGt2[9]={0.9156,0.9156,0.9156,0.9156,0.9156,0.9156,0.9156,0.9156,0.9156};
 double  ratioMGt3[9]={1.108,1.108,1.108,1.108,1.108,1.108,1.108,1.108,1.108};
 double  ratioMGt4[9]={1.083,1.083,1.083,1.083,1.083,1.083,1.083,1.083,1.083};
 double  ratioMGt5[9]={1.032,1.032,1.032,1.032,1.032,1.032,1.032,1.032,1.032};

 //bias in eta<0
 double  ratioMLt1[9]={0.9051,0.9051,0.9051,0.9051,0.9051,0.9051,0.9051,0.9051,0.9051};
 double  ratioMLt2[9]={1.103,1.103,1.103,1.103,1.103,1.103,1.103,1.103,1.103};
 double  ratioMLt3[9]={1.085,1.085,1.085,1.085,1.085,1.085,1.085,1.085,1.085};
 double  ratioMLt4[9]={1.042,1.042,1.042,1.042,1.042,1.042,1.042,1.042,1.042};
 double  ratioMLt5[9]={1.058,1.058,1.058,1.058,1.058,1.058,1.058,1.058,1.058};

 double biasMGt1[9]={0};double biasMGt2[9]={0};double biasMGt3[9]={0};double biasMGt4[9]={0};double biasMGt5[9]={0};
 double biasMLt1[9]={0};double biasMLt2[9]={0};double biasMLt3[9]={0};double biasMLt4[9]={0};double biasMLt5[9]={0};


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

 double combSysMGt1[9]={0};double combSysMGt2[9]={0};double combSysMGt3[9]={0};double combSysMGt4[9]={0};double combSysMGt5[9]={0};
 double combSysMLt1[9]={0};double combSysMLt2[9]={0};double combSysMLt3[9]={0};double combSysMLt4[9]={0};double combSysMLt5[9]={0};

 for(int i=0; i<9; i++){
  //trigger bias eta_pair>0
  biasMGt1[i]=(1-ratioMGt1[i])*TMath::Max(A_Tpc1Gt[i],AWerr_Tpc1Gt[i]);
  biasMGt2[i]=(1-ratioMGt2[i])*TMath::Max(A_Tpc2Gt[i],AWerr_Tpc2Gt[i]);
  biasMGt3[i]=(1-ratioMGt3[i])*TMath::Max(A_Tpc3Gt[i],AWerr_Tpc3Gt[i]);
  biasMGt4[i]=(1-ratioMGt4[i])*TMath::Max(A_Tpc4Gt[i],AWerr_Tpc4Gt[i]);
  biasMGt5[i]=(1-ratioMGt5[i])*TMath::Max(A_Tpc5Gt[i],AWerr_Tpc5Gt[i]);
  //trigger bias eta_pair<0      
  biasMLt1[i]=(1-ratioMLt1[i])*TMath::Max(A_Tpc1Lt[i],AWerr_Tpc1Lt[i]);
  biasMLt2[i]=(1-ratioMLt2[i])*TMath::Max(A_Tpc2Lt[i],AWerr_Tpc2Lt[i]);
  biasMLt3[i]=(1-ratioMLt3[i])*TMath::Max(A_Tpc3Lt[i],AWerr_Tpc3Lt[i]);
  biasMLt4[i]=(1-ratioMLt4[i])*TMath::Max(A_Tpc4Lt[i],AWerr_Tpc4Lt[i]);
  biasMLt5[i]=(1-ratioMLt5[i])*TMath::Max(A_Tpc5Lt[i],AWerr_Tpc5Lt[i]);


  //difference between A_UT with TPC only and TOF only

  rdiff1Gt[i]=(double)abs(A_Tpc1Gt[i]-A_Tof1Gt[i]);//absolute difference
  rdiff1errGt[i]=(double)sqrt(pow(AWerr_Tpc1Gt[i],2)+pow(AWerr_Tof1Gt[i],2));
  rdiff2Gt[i]=(double)abs(A_Tpc2Gt[i]-A_Tof2Gt[i]);
  rdiff2errGt[i]=(double)sqrt(pow(AWerr_Tpc2Gt[i],2)+pow(AWerr_Tof2Gt[i],2));
  rdiff3Gt[i]=(double)abs(A_Tpc3Gt[i]-A_Tof3Gt[i]);
  rdiff3errGt[i]=(double)sqrt(pow(AWerr_Tpc3Gt[i],2)+pow(AWerr_Tof3Gt[i],2));
  rdiff4Gt[i]=(double)abs(A_Tpc4Gt[i]-A_Tof4Gt[i]);
  rdiff4errGt[i]=(double)sqrt(pow(AWerr_Tpc4Gt[i],2)+pow(AWerr_Tof4Gt[i],2));
  rdiff5Gt[i]=(double)abs(A_Tpc5Gt[i]-A_Tof5Gt[i]);
  rdiff5errGt[i]=(double)sqrt(pow(AWerr_Tpc5Gt[i],2)+pow(AWerr_Tof5Gt[i],2));


  rrdiff1Gt[i]=(double)(A_Tpc1Gt[i]-A_Tof1Gt[i])/A_Tpc1Gt[i];//relative difference
  rrdiff1errGt[i]=(double)calErr(A_Tpc1Gt[i],AWerr_Tpc1Gt[i],A_Tof1Gt[i],AWerr_Tof1Gt[i]);//relative difference
  rrdiff2Gt[i]=(double)(A_Tpc2Gt[i]-A_Tof2Gt[i])/A_Tpc2Gt[i];
  rrdiff2errGt[i]=(double)calErr(A_Tpc2Gt[i],AWerr_Tpc2Gt[i],A_Tof2Gt[i],AWerr_Tof2Gt[i]);//relative difference
  rrdiff3Gt[i]=(double)(A_Tpc3Gt[i]-A_Tof3Gt[i])/A_Tpc3Gt[i];
  rrdiff3errGt[i]=(double)calErr(A_Tpc3Gt[i],AWerr_Tpc3Gt[i],A_Tof3Gt[i],AWerr_Tof3Gt[i]);//relative difference
  rrdiff4Gt[i]=(double)(A_Tpc4Gt[i]-A_Tof4Gt[i])/A_Tpc4Gt[i];
  rrdiff4errGt[i]=(double)calErr(A_Tpc4Gt[i],AWerr_Tpc4Gt[i],A_Tof4Gt[i],AWerr_Tof4Gt[i]);//relative difference
  rrdiff5Gt[i]=(double)(A_Tpc5Gt[i]-A_Tof5Gt[i])/A_Tpc5Gt[i];
  rrdiff5errGt[i]=(double)calErr(A_Tpc5Gt[i],AWerr_Tpc5Gt[i],A_Tof5Gt[i],AWerr_Tof5Gt[i]);//relative difference

  rrdiff1Lt[i]=(double)(A_Tpc1Lt[i]-A_Tof1Lt[i])/A_Tpc1Lt[i];//relative difference
  rrdiff1errLt[i]=(double)calErr(A_Tpc1Lt[i],AWerr_Tpc1Lt[i],A_Tof1Lt[i],AWerr_Tof1Lt[i]);//relative difference
  rrdiff2Lt[i]=(double)(A_Tpc2Lt[i]-A_Tof2Lt[i])/A_Tpc2Lt[i];
  rrdiff2errLt[i]=(double)calErr(A_Tpc2Lt[i],AWerr_Tpc2Lt[i],A_Tof2Lt[i],AWerr_Tof2Lt[i]);//relative difference
  rrdiff3Lt[i]=(double)(A_Tpc3Lt[i]-A_Tof3Lt[i])/A_Tpc3Lt[i];
  rrdiff3errLt[i]=(double)calErr(A_Tpc3Lt[i],AWerr_Tpc3Lt[i],A_Tof3Lt[i],AWerr_Tof3Lt[i]);//relative difference
  rrdiff4Lt[i]=(double)(A_Tpc4Lt[i]-A_Tof4Lt[i])/A_Tpc4Lt[i];
  rrdiff4errLt[i]=(double)calErr(A_Tpc4Lt[i],AWerr_Tpc4Lt[i],A_Tof4Lt[i],AWerr_Tof4Lt[i]);//relative difference
  rrdiff5Lt[i]=(double)(A_Tpc5Lt[i]-A_Tof5Lt[i])/A_Tpc5Lt[i];
  rrdiff5errLt[i]=(double)calErr(A_Tpc5Lt[i],AWerr_Tpc5Lt[i],A_Tof5Lt[i],AWerr_Tof5Lt[i]);//relative difference
cout<<"Bin:  "<<i<<" a: "<<A_Tpc1Lt[i]<<", ae: "<<AWerr_Tpc1Lt[i]<<", b: "<<A_Tof1Lt[i]<<", eb: "<<AWerr_Tof1Lt[i]<<", rdiff: "<<rrdiff1Lt[i]<<" ,rdiff err.: "<<rrdiff1errLt[i]<<endl;

  pidSysGt1[i]=(1.0-pionPairPtGtT[0])*abs(TMath::Max(A_F1Gt[i],AWerr_F1Gt[i]));
  pidSysGt2[i]=(1.0-pionPairPtGtT[1])*abs(TMath::Max(A_F2Gt[i],AWerr_F2Gt[i]));
  pidSysGt3[i]=(1.0-pionPairPtGtT[2])*abs(TMath::Max(A_F3Gt[i],AWerr_F3Gt[i]));
  pidSysGt4[i]=(1.0-pionPairPtGtT[3])*abs(TMath::Max(A_F4Gt[i],AWerr_F4Gt[i]));
  pidSysGt5[i]=(1.0-pionPairPtGtT[4])*abs(TMath::Max(A_F5Gt[i],AWerr_F5Gt[i]));

  pidSysLt1[i]=(1.0-pionPairPtLtT[0])*abs(TMath::Max(A_F1Lt[i],AWerr_F1Lt[i]));
  pidSysLt2[i]=(1.0-pionPairPtLtT[1])*abs(TMath::Max(A_F2Lt[i],AWerr_F2Lt[i]));
  pidSysLt3[i]=(1.0-pionPairPtLtT[2])*abs(TMath::Max(A_F3Lt[i],AWerr_F3Lt[i]));
  pidSysLt4[i]=(1.0-pionPairPtLtT[3])*abs(TMath::Max(A_F4Lt[i],AWerr_F4Lt[i]));
  pidSysLt5[i]=(1.0-pionPairPtLtT[4])*abs(TMath::Max(A_F5Lt[i],AWerr_F5Lt[i]));


  // cout<<A_F1[i]<<" , "<<Aerr_F1[i]<<", frac: "<<pionPairPtGtT[0]<<", imp: "<< double (1-pionPairPtGtT[0])<<" A: "<<double (abs(TMath::Max(A_F1[i],Aerr_F1[i])))<< ", Sys: "<< pidSysGt1[i]<<endl;
 }

 TGraphErrors *gr_AvsMFGt[5]; 
 TGraphErrors *gr_AvsMFLt[5]; 
 TGraphErrors *gr_AvsMgtS[5]; 
 TGraphErrors *gr_AvsMgtSc[5]; 
 TGraphErrors *gr_AvsMgtTOF[5]; 
 TGraphErrors *gr_rdiffGt[5]; 
 TGraphErrors *gr_rdiffLt[5]; 
 TGraphErrors *gr_trdiff[5]; 
 TGraphErrors *gr_rrdiffGt[5]; 
 TGraphErrors *gr_rrdiffLt[5]; 

 double M[9]={0};

 double rdiffGt[9]={0};
 double rdiffLt[9]={0};
 double rdifferrGt[9]={0};
 double rdifferrLt[9]={0};

 double rrdiffGt[9]={0};
 double rrdiffLt[9]={0};
 double rrdifferrGt[9]={0};
 double rrdifferrLt[9]={0};

 double biasMGt[9]={0};
 double biasMLt[9]={0};
 TLegend *leg[5];
 TLine *line[6];
 TLatex tex[6];

 double errx[9]={0};
 double errxps[9]={0.042,0.042,0.042,0.042,0.042,0.042,0.042,0.042,0.042};
 double errxts[9]={0.042,0.042,0.042,0.042,0.042,0.042,0.042,0.042,0.042};
 //double errxts[9]={0.038,0.038,0.038,0.038,0.038,0.038,0.038,0.038,0.038};

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
  if(pad==0)for(int mbin=0; mbin<9; mbin++){M[mbin]=M1[mbin];  A_FGt[mbin]=A_Tpc1Gt[mbin]; Aerr_FGt[mbin]=AWerr_Tpc1Gt[mbin]; A_FLt[mbin]=A_Tpc1Lt[mbin]; Aerr_FLt[mbin]=AWerr_Tpc1Lt[mbin]; rdiffGt[mbin]=rdiff1Gt[mbin]; rdifferrGt[mbin]=rdiff1errGt[mbin]; rdiffLt[mbin]=rdiff1Lt[mbin]; rdifferrLt[mbin]=rdiff1errLt[mbin];rrdiffGt[mbin]=rrdiff1Gt[mbin]; rrdifferrGt[mbin]=rrdiff1errGt[mbin]; rrdiffLt[mbin]=rrdiff1Lt[mbin]; rrdifferrLt[mbin]=rrdiff1errLt[mbin];} 
  if(pad==1)for(int mbin=0; mbin<9; mbin++){M[mbin]=M2[mbin];  A_FGt[mbin]=A_Tpc2Gt[mbin]; Aerr_FGt[mbin]=AWerr_Tpc2Gt[mbin]; A_FLt[mbin]=A_Tpc2Lt[mbin]; Aerr_FLt[mbin]=AWerr_Tpc2Lt[mbin]; rdiffGt[mbin]=rdiff2Gt[mbin]; rdifferrGt[mbin]=rdiff2errGt[mbin]; rdiffLt[mbin]=rdiff2Lt[mbin]; rdifferrLt[mbin]=rdiff2errLt[mbin];rrdiffGt[mbin]=rrdiff2Gt[mbin]; rrdifferrGt[mbin]=rrdiff2errGt[mbin]; rrdiffLt[mbin]=rrdiff2Lt[mbin]; rrdifferrLt[mbin]=rrdiff2errLt[mbin];} 
  if(pad==2)for(int mbin=0; mbin<9; mbin++){M[mbin]=M3[mbin];  A_FGt[mbin]=A_Tpc3Gt[mbin]; Aerr_FGt[mbin]=AWerr_Tpc3Gt[mbin]; A_FLt[mbin]=A_Tpc3Lt[mbin]; Aerr_FLt[mbin]=AWerr_Tpc3Lt[mbin]; rdiffGt[mbin]=rdiff3Gt[mbin]; rdifferrGt[mbin]=rdiff3errGt[mbin]; rdiffLt[mbin]=rdiff3Lt[mbin]; rdifferrLt[mbin]=rdiff3errLt[mbin];rrdiffGt[mbin]=rrdiff3Gt[mbin]; rrdifferrGt[mbin]=rrdiff3errGt[mbin]; rrdiffLt[mbin]=rrdiff3Lt[mbin]; rrdifferrLt[mbin]=rrdiff3errLt[mbin];} 
  if(pad==3)for(int mbin=0; mbin<9; mbin++){M[mbin]=M4[mbin];  A_FGt[mbin]=A_Tpc4Gt[mbin]; Aerr_FGt[mbin]=AWerr_Tpc4Gt[mbin]; A_FLt[mbin]=A_Tpc4Lt[mbin]; Aerr_FLt[mbin]=AWerr_Tpc4Lt[mbin]; rdiffGt[mbin]=rdiff4Gt[mbin]; rdifferrGt[mbin]=rdiff4errGt[mbin]; rdiffLt[mbin]=rdiff4Lt[mbin]; rdifferrLt[mbin]=rdiff4errLt[mbin];rrdiffGt[mbin]=rrdiff4Gt[mbin]; rrdifferrGt[mbin]=rrdiff4errGt[mbin]; rrdiffLt[mbin]=rrdiff4Lt[mbin]; rrdifferrLt[mbin]=rrdiff4errLt[mbin];} 
  if(pad==4)for(int mbin=0; mbin<9; mbin++){M[mbin]=M5[mbin];  A_FGt[mbin]=A_Tpc5Gt[mbin]; Aerr_FGt[mbin]=AWerr_Tpc5Gt[mbin]; A_FLt[mbin]=A_Tpc5Lt[mbin]; Aerr_FLt[mbin]=AWerr_Tpc5Lt[mbin]; rdiffGt[mbin]=rdiff5Gt[mbin]; rdifferrGt[mbin]=rdiff5errGt[mbin]; rdiffLt[mbin]=rdiff5Lt[mbin]; rdifferrLt[mbin]=rdiff5errLt[mbin];rrdiffGt[mbin]=rrdiff5Gt[mbin]; rrdifferrGt[mbin]=rrdiff5errGt[mbin]; rrdiffLt[mbin]=rdiff5Lt[mbin]; rrdifferrLt[mbin]=rrdiff5errLt[mbin];} 
  //draw final asymmetries
  gr_AvsMFGt[pad] = new TGraphErrors(9,M, A_FGt,errx,Aerr_FGt);//this draws asymmetry after polynomial TOF cut with all other tracks included
  gr_AvsMFGt[pad]->GetYaxis()-> SetTitle("#font[22]{A_{UT}^{Sin(#Phi_{S}-#Phi_{R})}}");
  gr_AvsMFGt[pad]->SetTitle("");
  gr_AvsMFGt[pad]->GetYaxis()->SetTitleOffset(1.2);
  gr_AvsMFGt[pad]->GetYaxis()->SetTitleSize(0.06);
  gr_AvsMFGt[pad]->GetYaxis()->SetLabelSize(0.06);
  gr_AvsMFGt[pad]->GetYaxis()->SetLabelFont(22);
  gr_AvsMFGt[pad]->GetXaxis()->SetTitle("#font[22]{M_{#pi^{+}#pi^{-}}(GeV/c^{2})}");
  gr_AvsMFGt[pad]->GetXaxis()->SetTitleSize(0.06);
  gr_AvsMFGt[pad]->GetXaxis()->SetLabelSize(0.06);
  gr_AvsMFGt[pad]->GetXaxis()-> SetLabelFont(22);
  gr_AvsMFGt[pad]->GetXaxis()->SetLimits(0.2,2.3);
  gr_AvsMFGt[pad]->GetYaxis()->SetRangeUser(-0.015, 0.06);
  gr_AvsMFGt[pad]-> SetMarkerStyle(20);
  gr_AvsMFGt[pad]-> SetMarkerColor(2);
  gr_AvsMFGt[pad]-> SetLineColor(2);
  gr_AvsMFGt[pad]-> SetLineWidth(1);

  gr_AvsMFLt[pad] = new TGraphErrors(9,M, A_FLt,errx,Aerr_FLt);//this draws asymmetry after polynomial TOF cut with all other tracks included
  gr_AvsMFLt[pad]-> SetMarkerStyle(20);
  gr_AvsMFLt[pad]-> SetMarkerColor(4);
  gr_AvsMFLt[pad]-> SetLineColor(4);
  gr_AvsMFLt[pad]-> SetLineWidth(1);

  gr_rdiffGt[pad] = new TGraphErrors(9,M, rdiffGt, 0, rdifferrGt);//this draws asymmetry after polynomial TOF cut with all other tracks included
  gr_rdiffGt[pad]->GetYaxis()-> SetTitle("#font[22]{|A_{1} - A_{2}|}");
  gr_rdiffGt[pad]->SetTitle("");
  gr_rdiffGt[pad]->GetYaxis()->SetTitleOffset(1.2);
  gr_rdiffGt[pad]->GetYaxis()->SetTitleSize(0.06);
  gr_rdiffGt[pad]->GetYaxis()->SetLabelSize(0.06);
  gr_rdiffGt[pad]->GetYaxis()->SetLabelFont(22);
  gr_rdiffGt[pad]->GetXaxis()->SetTitle("#font[22]{M_{#pi^{+}#pi^{-}}(GeV/c^{2})}");
  gr_rdiffGt[pad]->GetXaxis()->SetTitleSize(0.06);
  gr_rdiffGt[pad]->GetXaxis()->SetLabelSize(0.06);
  gr_rdiffGt[pad]->GetXaxis()-> SetLabelFont(22);
  gr_rdiffGt[pad]->GetXaxis()->SetLimits(0.2,2.3);
  gr_rdiffGt[pad]->GetYaxis()->SetRangeUser(-0.015, 0.06);
  gr_rdiffGt[pad]-> SetMarkerStyle(20);
  gr_rdiffGt[pad]-> SetMarkerColor(2);
  gr_rdiffGt[pad]-> SetLineColor(2);
  gr_rdiffGt[pad]-> SetLineWidth(1);

  gr_rdiffLt[pad] = new TGraphErrors(9,M, rdiffLt, 0, rdifferrLt);//this draws asymmetry after polynomial TOF cut with all other tracks included
  gr_rdiffLt[pad]-> SetMarkerStyle(20);
  gr_rdiffLt[pad]-> SetMarkerColor(4);
  gr_rdiffLt[pad]-> SetLineColor(4);
  gr_rdiffLt[pad]-> SetLineWidth(1);

  gr_rrdiffGt[pad] = new TGraphErrors(9,M, rrdiffGt, 0, rrdifferrGt);//this draws asymmetry after polynomial TOF cut with all other tracks included
  gr_rrdiffGt[pad]->GetYaxis()-> SetTitle("#font[22]{A_{tpc} - A_{tof}/A_{tpc}}");
  gr_rrdiffGt[pad]->SetTitle("");
  gr_rrdiffGt[pad]->GetYaxis()->SetTitleOffset(1.2);
  gr_rrdiffGt[pad]->GetYaxis()->SetTitleSize(0.06);
  gr_rrdiffGt[pad]->GetYaxis()->SetLabelSize(0.06);
  gr_rrdiffGt[pad]->GetYaxis()->SetLabelFont(22);
  gr_rrdiffGt[pad]->GetXaxis()->SetTitle("#font[22]{M_{#pi^{+}#pi^{-}}(GeV/c^{2})}");
  gr_rrdiffGt[pad]->GetXaxis()->SetTitleSize(0.06);
  gr_rrdiffGt[pad]->GetXaxis()->SetLabelSize(0.06);
  gr_rrdiffGt[pad]->GetXaxis()-> SetLabelFont(22);
  gr_rrdiffGt[pad]->GetXaxis()->SetLimits(0.2,2.3);
  gr_rrdiffGt[pad]->GetYaxis()->SetRangeUser(-20., 20);
  gr_rrdiffGt[pad]-> SetMarkerStyle(20);
  gr_rrdiffGt[pad]-> SetMarkerColor(2);
  gr_rrdiffGt[pad]-> SetLineColor(2);
  gr_rrdiffGt[pad]-> SetLineWidth(1);

  gr_rrdiffLt[pad] = new TGraphErrors(9,M, rrdiffLt, 0, rrdifferrLt);//this draws asymmetry after polynomial TOF cut with all other tracks included
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
   tex[pad].DrawLatex(0.35,34,Form("#font[22]{#color[2]{chi2/ndf  = %g / %i}}",gr_rrdiffGt[pad]->GetFunction("pol0")->GetChisquare(),gr_rrdiffGt[pad]->GetFunction("pol0")->GetNDF()));
   tex[pad].DrawLatex(0.35,30,Form("#font[22]{#color[2]{p0  = %g #pm %g}}",gr_rrdiffGt[pad]->GetFunction("pol0")->GetParameter(0),gr_rrdiffGt[pad]->GetFunction("pol0")->GetParError(0)));

   tex[pad].DrawLatex(0.35,25,Form("#font[22]{#color[4]{chi2/ndf  = %g / %i}}",gr_rrdiffLt[pad]->GetFunction("pol0")->GetChisquare(),gr_rrdiffLt[pad]->GetFunction("pol0")->GetNDF()));
   tex[pad].DrawLatex(0.35,21,Form("#font[22]{#color[4]{p0  = %g #pm %g}}",gr_rrdiffLt[pad]->GetFunction("pol0")->GetParameter(0),gr_rrdiffLt[pad]->GetFunction("pol0")->GetParError(0)));

   tex[pad].DrawLatex(1.3,34,Form("#font[22]{#color[1]{< p_{T} > = %g GeV/c}}",avgPt[pad]));

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
   gr_rdiffGt[pad]->Draw("AP");
   if(pad==1 || pad==2 || pad==4)gr_rdiffGt[pad]->GetYaxis()->SetLabelSize(0);
   gr_rdiffGt[pad]->GetYaxis()->SetNdivisions(505);
   gr_rdiffGt[pad]->GetXaxis()->SetNdivisions(505);

   gr_rdiffGt[pad]->Fit("pol0");
   fp0Gt[pad]=(double)gr_rdiffGt[pad]->GetFunction("pol0")->GetParameter(0);
   gr_rdiffLt[pad]->Draw("SAME P");
   gr_rdiffLt[pad]->Fit("pol0");
   fp0Lt[pad]=(double)gr_rdiffLt[pad]->GetFunction("pol0")->GetParameter(0);
   gr_rdiffLt[pad]->GetFunction("pol0")->SetLineColor(4);
   cout<<fp0Gt[pad]<<", "<<fp0Lt[pad]<<endl;

   gPad->Update();

   line[pad]=  new TLine(cdiff->cd(pad+1)->GetUxmin(),0.,cdiff->cd(pad+1)->GetUxmax(),0.);
   line[pad]->SetLineStyle(2);
   line[pad]->SetLineWidth(1);
   line[pad]->Draw();


   tex[pad].SetTextSize(0.04);
   tex[pad].SetTextAlign(13);
   tex[pad].DrawLatex(0.35,0.055,Form("#font[22]{#color[2]{chi2/ndf  = %g #pm %i}}",gr_rdiffGt[pad]->GetFunction("pol0")->GetChisquare(),gr_rdiffGt[pad]->GetFunction("pol0")->GetNDF()));
   tex[pad].DrawLatex(0.35,0.05,Form("#font[22]{#color[2]{p0  = %g #pm %g}}",gr_rdiffGt[pad]->GetFunction("pol0")->GetParameter(0),gr_rdiffGt[pad]->GetFunction("pol0")->GetParError(0)));

   tex[pad].DrawLatex(0.35,0.045,Form("#font[22]{#color[4]{chi2/ndf  = %g #pm %i}}",gr_rdiffLt[pad]->GetFunction("pol0")->GetChisquare(),gr_rdiffLt[pad]->GetFunction("pol0")->GetNDF()));
   tex[pad].DrawLatex(0.35,0.040,Form("#font[22]{#color[4]{p0  = %g #pm %g}}",gr_rdiffLt[pad]->GetFunction("pol0")->GetParameter(0),gr_rdiffLt[pad]->GetFunction("pol0")->GetParError(0)));

   tex[pad].DrawLatex(1.25,0.055,Form("#font[22]{#color[1]{< p_{T} > = %g GeV/c}}",avgPt[pad]));

   gPad->Update();
  }
 }
 cdiff->Update();
 cdiff->SaveAs("./Plots/pidsysRatioPlot_Minv.pdf");
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

 /* //final systematics. Combined trigger bias and PID as quadrature sum
  combSysMGt1[i]=sqrt(pow(biasMGt1[i],2) + pow(fpidSysGt1[i],2));
  combSysMGt2[i]=sqrt(pow(biasMGt2[i],2) + pow(fpidSysGt2[i],2));
  combSysMGt3[i]=sqrt(pow(biasMGt3[i],2) + pow(fpidSysGt3[i],2));
  combSysMGt4[i]=sqrt(pow(biasMGt4[i],2) + pow(fpidSysGt4[i],2));
  combSysMGt5[i]=sqrt(pow(biasMGt5[i],2) + pow(fpidSysGt5[i],2));
  //cout<<"trig bias: "<<biasMGt1[i]<<", pid: "<<fpidSysGt1[i]<<",  combined: "<<combSysMGt1[i]<<endl;
  combSysMLt1[i]=sqrt(pow(biasMLt1[i],2) + pow(fpidSysLt1[i],2));
  combSysMLt2[i]=sqrt(pow(biasMLt2[i],2) + pow(fpidSysLt2[i],2));
  combSysMLt3[i]=sqrt(pow(biasMLt3[i],2) + pow(fpidSysLt3[i],2));
  combSysMLt4[i]=sqrt(pow(biasMLt4[i],2) + pow(fpidSysLt4[i],2));
  combSysMLt5[i]=sqrt(pow(biasMLt5[i],2) + pow(fpidSysLt5[i],2));
*/
  //final systematics. Combined trigger bias and PID as quadrature sum
  combSysMGt1[i]=sqrt(pow(biasMGt1[i],2) + pow(frp0Gt[0]*A_Tpc1Gt[i],2));
  combSysMGt2[i]=sqrt(pow(biasMGt2[i],2) + pow(frp0Gt[1]*A_Tpc2Gt[i],2));
  combSysMGt3[i]=sqrt(pow(biasMGt3[i],2) + pow(frp0Gt[2]*A_Tpc3Gt[i],2));
  combSysMGt4[i]=sqrt(pow(biasMGt4[i],2) + pow(frp0Gt[3]*A_Tpc4Gt[i],2));
  combSysMGt5[i]=sqrt(pow(biasMGt5[i],2) + pow(frp0Gt[4]*A_Tpc5Gt[i],2));
  //cout<<"trig bias: "<<biasMGt1[i]<<", pid: "<<Gt1[i]<<",  combined: "<<combSysMGt1[i]<<endl;
  combSysMLt1[i]=sqrt(pow(biasMLt1[i],2) + pow(frp0Lt[0]*A_Tpc1Lt[i],2));
  combSysMLt2[i]=sqrt(pow(biasMLt2[i],2) + pow(frp0Lt[1]*A_Tpc2Lt[i],2));
  combSysMLt3[i]=sqrt(pow(biasMLt3[i],2) + pow(frp0Lt[2]*A_Tpc3Lt[i],2));
  combSysMLt4[i]=sqrt(pow(biasMLt4[i],2) + pow(frp0Lt[3]*A_Tpc4Lt[i],2));
  combSysMLt5[i]=sqrt(pow(biasMLt5[i],2) + pow(frp0Lt[4]*A_Tpc5Lt[i],2));

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
   if(pad==0)for(int mbin=0; mbin<9; mbin++){ M[mbin]=M1[mbin];  A_FGt[mbin]=A_Tpc1Gt[mbin]; Aerr_FGt[mbin]=AWerr_Tpc1Gt[mbin]; A_FLt[mbin]=A_Tpc1Lt[mbin]; Aerr_FLt[mbin]=AWerr_Tpc1Lt[mbin]; finalSysGt[mbin]=combSysMGt1[mbin]; finalSysLt[mbin]=combSysMLt1[mbin]; }
   if(pad==1)for(int mbin=0; mbin<9; mbin++){M[mbin]=M2[mbin];  A_FGt[mbin]=A_Tpc2Gt[mbin]; Aerr_FGt[mbin]=AWerr_Tpc2Gt[mbin]; A_FLt[mbin]=A_Tpc2Lt[mbin]; Aerr_FLt[mbin]=AWerr_Tpc2Lt[mbin];  finalSysGt[mbin]=combSysMGt2[mbin]; finalSysLt[mbin]=combSysMLt2[mbin]; }
   if(pad==2)for(int mbin=0; mbin<9; mbin++){M[mbin]=M3[mbin];  A_FGt[mbin]=A_Tpc3Gt[mbin]; Aerr_FGt[mbin]=AWerr_Tpc3Gt[mbin]; A_FLt[mbin]=A_Tpc3Lt[mbin]; Aerr_FLt[mbin]=AWerr_Tpc3Lt[mbin];  finalSysGt[mbin]=combSysMGt3[mbin]; finalSysLt[mbin]=combSysMLt3[mbin]; }
   if(pad==3)for(int mbin=0; mbin<9; mbin++){M[mbin]=M4[mbin];  A_FGt[mbin]=A_Tpc4Gt[mbin]; Aerr_FGt[mbin]=AWerr_Tpc4Gt[mbin]; A_FLt[mbin]=A_Tpc4Lt[mbin]; Aerr_FLt[mbin]=AWerr_Tpc4Lt[mbin];  finalSysGt[mbin]=combSysMGt4[mbin]; finalSysLt[mbin]=combSysMLt4[mbin];  }
   if(pad==4)for(int mbin=0; mbin<9; mbin++){M[mbin]=M5[mbin];  A_FGt[mbin]=A_Tpc5Gt[mbin]; Aerr_FGt[mbin]=AWerr_Tpc5Gt[mbin]; A_FLt[mbin]=A_Tpc5Lt[mbin]; Aerr_FLt[mbin]=AWerr_Tpc5Lt[mbin];  finalSysGt[mbin]=combSysMGt5[mbin]; finalSysLt[mbin]=combSysMLt5[mbin];  }

   gr_FGtS[pad] = new TGraphErrors(9,M,A_FGt,errxts,finalSysGt);
   gr_FGtS[pad]->GetYaxis()-> SetTitle("#font[22]{A_{UT}^{Sin(#Phi_{S}-#Phi_{R})}} ");
   gr_FGtS[pad]->SetTitle("");
   gr_FGtS[pad]->GetYaxis()->SetTitleOffset(1.2);
   gr_FGtS[pad]->GetYaxis()->SetTitleSize(0.06);
   gr_FGtS[pad]->GetYaxis()->SetLabelSize(0.06);
   gr_FGtS[pad]->GetYaxis()->SetLabelFont(22);
   gr_FGtS[pad]->GetXaxis()->SetTitle("#font[22]{M_{#pi^{+}#pi^{-}} (GeV/c^{2})}   ");
   gr_FGtS[pad]->GetXaxis()->SetTitleSize(0.06);
   gr_FGtS[pad]->GetXaxis()->SetLabelSize(0.06);
   gr_FGtS[pad]->GetXaxis()-> SetLabelFont(22);
   gr_FGtS[pad]->GetXaxis()->SetLimits(0.2,2.3);
   gr_FGtS[pad]->GetYaxis()->SetRangeUser(-0.015, 0.063);
   gr_FGtS[pad]->SetFillStyle(0000);
   gr_FGtS[pad]->SetLineColor(1);
   gr_FGtS[pad]->SetLineWidth(2);


   gr_FLtS[pad] = new TGraphErrors(9,M,A_FLt,errxts,finalSysLt);
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
   gr_AvsMFGt[pad]->Draw("SAME P");
   gr_AvsMFLt[pad]->Draw("SAME P");

   gPad->Update();
   line[pad]=  new TLine(cfinal->cd(pad+1)->GetUxmin(),0.,cfinal->cd(pad+1)->GetUxmax(),0.);
   line[pad]->SetLineStyle(2);
   line[pad]->SetLineWidth(1);
   line[pad]->Draw();


   tex[pad].SetTextSize(0.06);
   tex[pad].SetTextAlign(13);
   tex[pad].DrawLatex(1.25,0.055,Form("#font[22]{#color[1]{< p_{T} > = %g GeV/c}}",avgPt[pad])); 

   if(pad==0){
    leg[pad] = new TLegend(0.24,.6, 0.55, 0.8);
    leg[pad]->AddEntry(gr_AvsMFGt[pad], " #font[22]{ #eta^{#pi^{+}#pi^{-}} > 0}", "lp");
    leg[pad]->AddEntry(gr_AvsMFLt[pad], " #font[22]{ #eta^{#pi^{+}#pi^{-}} < 0}", "lp");
    leg[pad]->AddEntry(gr_FGtS[pad], " #font[22]{ Syst. Error}", "f");
    leg[pad]->SetTextSize(0.05);
    leg[pad]->Draw();
   }



  }



 }
 cfinal->Update();
 cfinal->SaveAs("./Plots/finalAut_Minv.pdf");

}

//error calculation
double calErr(double a, double ea, double b, double eb){
 double err=-999;
 //err=(b/a)*sqrt(pow(ea/a,2)+pow(eb/b,2));
 err=(1/(a*a))*sqrt(pow(b*ea,2)+pow(a*eb,2));
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


