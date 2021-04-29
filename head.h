#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <vector>
#include <cmath>
#include <complex>
#include <quadmath.h>
#include <gsl/gsl_sf_dilog.h>
#include <qcdloop/qcdloop.h>


/*  Numeric values for the kinematic parameters in GeV.                         */

#define Etot 10.583                     /* total energy in center of mom. frame */
#define ME 0.0005109989461              /*     Electron Rest mass from pdg      */
#define MMU 0.1056583745                /*       Muon Rest mass from pdg        */
#define MT 1.77686                      /*      Tauon Rest mass from pdg        */

#define MU 0.002160                     /*     up-quark Rest mass from pdg      */
#define MC 1.27                         /*   charm-quark Rest mass from pdg     */
#define MTO 172.9                       /*    top-quark Rest mass from pdg      */

#define MD 0.00467                      /*    down-quark Rest mass from pdg     */
#define MS 0.093                        /*  strange-quark Rest mass from pdg    */
#define MB 4.18                         /*   bottom-quark Rest mass from pdg    */

#define MW 80.379                       /*      W Boson rest mass from pdg      */
#define MZ 91.1876                      /*      Z Boson rest mass from pdg      */
#define MH 125                          /*    Higgs Boson rest mass from pdg    */

#define Eem 7.0                         /*   Electron energy in the lab frame   */
#define Eep 4.0                         /*  Positron energy in the lab system   */

const double shat = Etot*Etot,          /*             Mandelstam s             */
             ME2 = ME*ME,               /*        squared Electron Mass         */
             MMU2 = MMU*MMU,            /*          squared Muon Mass           */
             MT2 = MT*MT,               /*          squared Tauon Mass          */
             MTO2 = MTO*MTO,            /*           squared Top Mass           */
             MW2 = MW*MW,               /*            squared W Mass            */
             MZ2 = MZ*MZ,               /*            squared Z Mass            */
             MUT = MT/Etot;             /*          reduced Tauon mass          */
double absp3,                           /*        final state 3-momentum        */
       ymax;                            /*           maximum rapidity           */
//       propg = 1/(Etot*Etot),           /*         Photon-Propagator/I          */
//       propZ = 1/(Etot*Etot - MZ*MZ);   /*            Z-Propagator/I            */


/*  Numeric values for the coupling parameters.
    Do the gfm and gfp actually depend on fermion charge sign?                  */

#define sw 0.47183314385900446          /*      sine of the Weinberg angle      */
#define cw 0.88168786106882973          /*     cosine of the Weinberg angle     */
#define sw2 0.22262651564387198         /*  squared sine of the Weinberg angle  */
//#define alpha 0.0072973525664           /*       fine structure constant        */
#define alpha 0.007555915561462338      /*       fine structure constant        */
#define gfermi 0.000011663787           /*           fermi constant             */
#define hc2oGeV2 0.3893793721           /*    conversion rate to mb from pdg    */

double cw2 = 1.0 - sw2,                 /*  squared sine of the Weinberg angle  */
       af = -1.0/4/sw/cw,               /*     axial vector coupling of the Z   */
       vf =(4*sw2-1)/4/sw/cw,           /*      vector coupling of the Z        */
       gp = sw/cw,                      /*      right hand LEPTON coupling      */
       gm = (sw2-0.5)/sw/cw,            /*      left hand LEPTON coupling       */
       e2 = alpha*4*M_PI,               /*   second power of electron charge    */
       e4 = alpha*alpha*16*M_PI*M_PI,   /*   fourth power of electron charge    */
       e6 = e4*e2,                      /*   fourth power of electron charge    */
       CF = 1.0;                        /*   casimir F  ??? Check the value pls */
       

/*  Phase space constants.                                                      */
//#define ABST 4.97734804292406           /*    3-momentum of final state taus    */
//           = sqrt(shat/4-MT2)

// const double CPS2 = ABST/8/M_PI         /* constant for 2 particle phase space  */
//                     /Etot/Etot/Etot,
//              CPS3 = (Etot*Etot-4*MT*MT) /* constant for 3 particle phase space  */
//                     /128/M_PI/M_PI/M_PI/Etot/Etot;

       
/*  Control sequence for the integration.                                       */

#define NSMP 2000000000                 /*          number of samples           */
#define NIT 15                          /*         number iterations            */
#define NBINS 50                        /*           number of bins             */

const double scut = 1e-7,               /*        Soft cutoff in E5^2/s         */
             ccut = 1e-7,               /*     Collinear cutoff in pi.p5/s      */
             // 1e-8 tested and stable at NSMP=2000000000
             ecut = 1e-14;               /*      upper cut for x integration     */
double acut_II,acut_IF,acut_FI,acut_FF, /*    Dipole upper cutoffs in E5/Etot   */
       acut_FFp,yp,                     /*   parameters of final-final Dipole   */
       Ieik[2],                         /*      eikonal integral plus part      */
       IgQcol[2],                       /*          collinear integral          */
       Del_Ieik,                        /*        final-final alpha term        */
       Del_IgQcol,                      /*        final-final alpha term        */
       wmax,                            /*           maximum weight             */
       itBin[NIT][3][NBINS+2],          /*     0: costh, 1: p_T, 2: rapidity    */
       fbin[3][NBINS+2],                /*     0: costh, 1: p_T, 2: rapidity    */
       dbin[3][NBINS+2];                /* errors 0: costh, 1: p_T, 2: rapidity */
int bin = 1,                            /*          binning momentum            */
    cEps = 0,                           /*      calculate poles in 1/eps        */
    colQED = 1,                         /*    calculate electron mass terms     */
    countIt;                            /*           iteration counter          */
unsigned long countN;                   /*             event counter            */

/*  Control sequence for the virtual correction.                                 */

// const double _Complex I = sqrt(-1);
typedef std::complex<double> complex;   /*       convenient use of type         */
const complex I(0.0,1.0);               /*             this is i                */
double muR,                             /*       renormalization scale          */   
       muR2,                            /*    squared renormalization scale     */                    
       muF2,                            /*        factorization scale           */
       logMu,                           /*              log(mu2)                */
       logMuMT,                         /*            log(mu2/MT2)              */
       logMuMMU,                        /*            log(mu2/MMU2)             */
       logMuME,                         /*            log(mu2/ME2)              */
       logMuShat,                       /*            log(mu2/shat)             */
       logShat,                         /*             log(shat)                */
       logMe,                           /*              log(ME2)                */
       logMeShat;                       /*            log(ME2/shat)             */
       

/*  Global variables for storage of PS independent parts of virtual Correction. */
std::vector<complex> ElCT(4),           /*    Electron triangle counter term    */
                     TauCT(4),          /*     Tauon triangle counter term      */
                     ElBCT(4),          /*     Electron bubble counter term     */
                     mElBCT(4),         /* massive Electron bubble counter term */
                     MuBCT(4),          /*       Muon bubble counter term       */
                     TauBCT(4),         /*      Tauon bubble counter term       */
                     A0MT(4),           /*              A0(MT)                  */
                     A0MMU(4),          /*              A0(MMU)                 */
                     A0ME(4),           /*              A0(ME)                  */
                     B00(4),            /*             B0(0,0,0)                */
                     B0s(4),            /*           B0(p1+p2,0,0)              */
                     B0MT(4),           /*            B0(MT,0,MT)               */
                     B00MT(4),          /*            B0(0,MT,MT)               */
                     B0MT0(4),          /*            B0(0,MT,0)                */
                     B00MMU(4),         /*           B0(0,MMU,MMU)              */
                     B00ME(4),          /*            B0(0,ME,ME)               */
                     B0sMT(4),          /*          B0(p3+p4,MT,MT)             */
                     B0sMMU(4),         /*         B0(p3+p4,MMU,MMU)            */
                     B0sME(4),          /*          B0(p3+p4,ME,ME)             */
                     C0s(4),            /*         C0(p1,p1+p2,0,0,0)           */
                     C0MTs(4),          /*         C0(MT,p1+p2,0,MT,0)          */
                     C0sMT(4);          /*        C0(p3,p3+p4,MT,0,MT)          */
             
             
const double twoPi2 = 4*M_PI*M_PI,
             twoPi3 = 8*M_PI*M_PI*M_PI,
             twoPi4 = 16*M_PI*M_PI*M_PI*M_PI,
             twoPi5 = 32*M_PI*M_PI*M_PI*M_PI*M_PI;

             
