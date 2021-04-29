C         Last modification on December 20th 1999 by M.S.
C ==================================================================
C ================= PROGRAM HDECAY: COMMENTS =======================
C ==================================================================
C
C                       ***************
C                       * VERSION 2.0 *
C                       ***************
C
C
C  This program calculates the total decay widths and the branching 
C  ratios of the C Standard Model Higgs boson (HSM) as well as those 
C  of the neutral (HL= the light CP-even, HH= the heavy CP-even, HA= 
C  the pseudoscalar) and the charged (HC) Higgs bosons of the Minimal
C  Supersymmetric extension of the Standard Model (MSSM). It includes:
C
C - All the decay channels which are kinematically allowed and which
C   have branching ratios larger than 10**(-4). 
C
C - All QCD corrections to the fermionic and gluonic decay modes.
C   Most of these corrections are mapped into running masses in a
C   consistent way with some freedom for including high order terms. 
C
C - Below--threshold three--body decays with off--shell top quarks
C   or ONE off-shell gauge boson, as well as some decays with one
C   off-shell Higgs boson in the MSSM. 
C
C - Double off-shell decays: HSM,HL,HH --> W*W*,Z*Z* -->4 fermions,
C   which could be important for Higgs masses close to MW or MZ.
C
C - In the MSSM, the radiative corrections with full squark mixing and 
C   uses the RG improved values of Higgs masses and couplings with the 
C   main NLO corrections implemented (based on M.Carena, M. Quiros and
C   C.E.M. Wagner, Nucl. Phys. B461 (1996) 407, hep-ph/9508343). 
C
C - In the MSSM, all the decays into CHARGINOS, NEUTRALINOS, SLEPTONS 
C   and SQUARKS (with mixing in the stop and sbottom sectors). 
C
C - Chargino, slepton and squark loops in the 2 photon decays and squark
C   loops in the gluonic decays (including QCD corrections). 
C
C  ===================================================================
C  This program has been written by A.Djouadi, J.Kalinowski and M.Spira.
C  For details on how to use the program see: Comp. Phys. Commun. 108
C  (1998) 56, hep-ph/9704448. For any question, comment, suggestion or
C  complaint, please contact us at:
C          djouadi@lpm.univ-montp2.fr
C          kalino@fuw.edu.pl
C          Michael.Spira@cern.ch


C ================ IT USES AS INPUT PARAMETERS:
C
C   IHIGGS: =0: CALCULATE BRANCHING RATIOS OF SM HIGGS BOSON
C           =1: CALCULATE BRANCHING RATIOS OF MSSM h BOSON
C           =2: CALCULATE BRANCHING RATIOS OF MSSM H BOSON
C           =3: CALCULATE BRANCHING RATIOS OF MSSM A BOSON
C           =4: CALCULATE BRANCHING RATIOS OF MSSM H+ BOSON
C           =5: CALCULATE BRANCHING RATIOS OF ALL MSSM HIGGS BOSONS
C
C TGBET:    TAN(BETA) FOR MSSM
C MABEG:    START VALUE OF M_A FOR MSSM AND M_H FOR SM
C MAEND:    END VALUE OF M_A FOR MSSM AND M_H FOR SM
C NMA:      NUMBER OF ITERATIONS FOR M_A
C ALS(MZ):  VALUE FOR ALPHA_S(M_Z)
C MSBAR(1): MSBAR MASS OF STRANGE QUARK AT SCALE Q=1 GEV
C MC:       CHARM POLE MASS
C MB:       BOTTOM POLE MASS
C MT:       TOP POLE MASS
C MTAU:     TAU MASS
C MMUON:    MUON MASS
C ALPH:     INVERSE QED COUPLING
C GF:       FERMI CONSTANT
C GAMW:     W WIDTH
C GAMZ:     Z WIDTH
C MZ:       Z MASS
C MW:       W MASS
C VUS:      CKM PARAMETER V_US
C VCB:      CKM PARAMETER V_CB
C VUB/VCB:  RATIO V_UB/V_CB
C 1ST AND 2ND GENERATION:
C MSL1:      SUSY BREAKING MASS PARAMETERS OF LEFT HANDED SLEPTONS 
C MER1:      SUSY BREAKING MASS PARAMETERS OF RIGHT HANDED SLEPTONS 
C MQL1:      SUSY BREAKING MASS PARAMETERS OF LEFT HANDED SUPS
C MUR1:      SUSY BREAKING MASS PARAMETERS OF RIGHT HANDED SUPS
C MDR1:      SUSY BREAKING MASS PARAMETERS OF RIGHT HANDED SDOWNS 
C 3RD GENERATION:
C MSL:      SUSY BREAKING MASS PARAMETERS OF LEFT HANDED STAUS 
C MER:      SUSY BREAKING MASS PARAMETERS OF RIGHT HANDED STAUS 
C MSQ:      SUSY BREAKING MASS PARAMETERS OF LEFT HANDED STOPS
C MUR:      SUSY BREAKING MASS PARAMETERS OF RIGHT HANDED STOPS
C MDR:      SUSY BREAKING MASS PARAMETERS OF RIGHT HANDED SBOTTOMS 
C AL:       STAU TRILINEAR SOFT BREAKING TERMS 
C AU:       STOP TRILINEAR SOFT BREAKING TERMS.
C AD:       SBOTTOM TRILINEAR SOFT BREAKING TERMS.
C MU:       SUSY HIGGS MASS PARAMETER
C M2:       GAUGINO MASS PARAMETER. 
C
C NNLO (M): =1: USE O(ALPHA_S) FORMULA FOR POLE MASS --> MSBAR MASS
C           =2: USE O(ALPHA_S**2) FORMULA FOR POLE MASS --> MSBAR MASS
C
C ON-SHELL: =0: INCLUDE OFF_SHELL DECAYS H,A --> T*T*, A --> Z*H,
C               H --> W*H+,Z*A, H+ --> W*A, W*H, T*B
C           =1: EXCLUDE THE OFF-SHELL DECAYS ABOVE
C
C ON-SH-WZ: =0: INCLUDE DOUBLE OFF-SHELL PAIR DECAYS PHI --> W*W*,Z*Z*
C           =1: INCLUDE ONLY SINGLE OFF-SHELL DECAYS PHI --> W*W,Z*Z
C
C IPOLE:    =0 COMPUTES RUNNING HIGGS MASSES (FASTER) 
C           =1 COMPUTES POLE HIGGS MASSES 
C
C OFF-SUSY: =0: INCLUDE DECAYS (AND LOOPS) INTO SUPERSYMMETRIC PARTICLES
C           =1: EXCLUDE DECAYS (AND LOOPS) INTO SUPERSYMMETRIC PARTICLES
C
C INIDEC:   =0: PRINT OUT SUMS OF CHARGINO/NEUTRALINO/SFERMION DECAYS
C           =1: PRINT OUT INDIVIDUAL CHARGINO/NEUTRALINO/SFERMION DECAYS
C
C NF-GG:    NUMBER OF LIGHT FLAVORS INCLUDED IN THE GLUONIC DECAYS 
C            PHI --> GG* --> GQQ (3,4 OR 5)
C           
C =======================================================================
C ============== BEGINNING OF THE MAIN PROGRAM ==========================
C =======================================================================
C
      subroutine HDECAY
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(K=6,NI=87,NSA=85,NSB=86,NLA=88,NLB=89,NHA=90,NHB=91,
     .          NHC=92,NAA=93,NAB=94,NCA=95,NCB=96,NRA=97,NRB=98,
     .          NSUSYL=81,NSUSYA=82,NSUSYH=83,NSUSYC=84,NPAR=80,
     .          NSUSYLA=79,NSUSYLB=78,NSUSYLC=77,NSUSYLD=76,NSUSYLE=75,
     .          NSUSYLF=59,NSUSYHF=58,
     .          NSUSYHA=74,NSUSYHB=73,NSUSYHC=72,NSUSYHD=71,NSUSYHE=70,
     .          NSUSYAA=69,NSUSYAB=68,NSUSYAC=67,NSUSYAD=66,NSUSYAE=65,
     .          NSUSYCA=64,NSUSYCB=63,NSUSYCC=62,NSUSYCD=61,NSUSYCE=60)
      DIMENSION GMN(4),XMN(4),GMC(2),GMST(2),GMSB(2),GMSL(2),
     .          GMSU(2),GMSD(2),GMSE(2),GMSN(2)
      DIMENSION HLBRSC(2,2),HLBRSN(4,4),HHBRSC(2,2),HHBRSN(4,4),
     .          HABRSC(2,2),HABRSN(4,4),HCBRSU(2,4),
     .          HHBRST(2,2),HHBRSB(2,2),HCBRSTB(2,2) 
      DIMENSION AC1(2,2),AC2(2,2),AC3(2,2),
     .          AN1(4,4),AN2(4,4),AN3(4,4),
     .          ACNL(2,4),ACNR(2,4)
      DIMENSION GLTT(2,2),GLBB(2,2),GHTT(2,2),GHBB(2,2),GCTB(2,2),
     .          GLEE(2,2),GHEE(2,2),GCEN(2,2)
      COMMON/MASSES/AMS,AMC,AMB,AMT
      COMMON/STRANGE/AMSB
      COMMON/PARAM/GF,ALPH,AMTAU,AMMUON,AMZ,AMW
      COMMON/CKMPAR/VUS,VCB,VUB
      COMMON/HMASS/AMSM,AMA,AML,AMH,AMCH,AMAR
      COMMON/BREAK/AMEL,AMER,AMSQ,AMUR,AMDR,AL,AU,AD,AMU,AM2
      COMMON/SFER1ST/AMQL1,AMUR1,AMDR1,AMEL1,AMER1
      COMMON/WZWDTH/GAMC0,GAMT0,GAMT1,GAMW,GAMZ
      COMMON/COUP/GAT,GAB,GLT,GLB,GHT,GHB,GZAH,GZAL,
     .            GHHH,GLLL,GHLL,GLHH,GHAA,GLAA,GLVV,GHVV,
     .            GLPM,GHPM,B,A
      COMMON/ALS/XLAMBDA,AMC0,AMB0,AMT0,N0
      COMMON/FLAG/IHIGGS,NNLO,IPOLE
      COMMON/ONSHELL/IONSH,IONWZ,IOFSUSY
      COMMON/OLDFASH/NFGG
      COMMON/WIDTHSM/SMBRB,SMBRL,SMBRM,SMBRS,SMBRC,SMBRT,SMBRG,SMBRGA,
     .               SMBRZGA,SMBRW,SMBRZ,SMWDTH
      COMMON/WIDTHA/ABRB,ABRL,ABRM,ABRS,ABRC,ABRT,ABRG,ABRGA,ABRZGA,
     .              ABRZ,AWDTH
      COMMON/WIDTHHL/HLBRB,HLBRL,HLBRM,HLBRS,HLBRC,HLBRT,HLBRG,HLBRGA,
     .               HLBRZGA,HLBRW,HLBRZ,HLBRA,HLBRAZ,HLBRHW,HLWDTH
      COMMON/WIDTHHH/HHBRB,HHBRL,HHBRM,HHBRS,HHBRC,HHBRT,HHBRG,HHBRGA,
     .               HHBRZGA,HHBRW,HHBRZ,HHBRH,HHBRA,HHBRAZ,HHBRHW,
     .               HHWDTH
      COMMON/WIDTHHC/HCBRB,HCBRL,HCBRM,HCBRBU,HCBRS,HCBRC,HCBRT,HCBRW,
     .               HCBRA,HCWDTH
      COMMON/WISUSY/HLBRSC,HLBRSN,HHBRSC,HHBRSN,HABRSC,HABRSN,HCBRSU,
     .              HLBRCHT,HHBRCHT,HABRCHT,HLBRNET,HHBRNET,HABRNET,
     .              HCBRCNT,HLBRSL,HHBRSL,HCBRSL,HABRSL,HABRST,HABRSB,
     .              HHBRSQ,HHBRST,HHBRSB,HHBRSQT,HCBRSQ,HCBRSTB,
     .              HCBRSQT,HLBRSQ,HLBRSQT
      COMMON/WISFER/BHLSLNL,BHLSLEL,BHLSLER,BHLSQUL,BHLSQUR,BHLSQDL,
     .              BHLSQDR,BHLST(2,2),BHLSB(2,2),BHLSTAU(2,2),
     .              BHHSLNL,BHHSLEL,BHHSLER,BHHSQUL,BHHSQUR,BHHSQDL,
     .              BHHSQDR,BHHST(2,2),BHHSB(2,2),BHHSTAU(2,2),
     .              BHASTAU,BHASB,BHAST,
     .              BHCSL00,BHCSL11,BHCSL21,BHCSQ,BHCSTB(2,2)
      COMMON/SMASS/GMN,XMN,GMC,GMST,GMSB,GMSL,GMSU,GMSD,GMSE,GMSN 

      common /hdparms/ alsmz, amabeg, amaend, tgbet

      parameter( PI = 3.14159265358979323846d0 )
C      PI = 4*DATAN(1D0)

C      OPEN(NI,FILE='hdecay.in')
C      OPEN(NPAR,FILE='br.in')

C      READ(NI,101)IHIGGS
C      READ(NI,100)TGBET
C      READ(NI,100)AMABEG
C      READ(NI,100)AMAEND
C      READ(NI,101)NMA
      nma = 1
C      READ(NI,100)ALSMZ
      alsmz = 0.118d0
C      READ(NI,100)AMS
      ams = 0.190d0
C      READ(NI,100)AMC
C      READ(NI,100)AMB
C      READ(NI,100)AMT
      amc = 1.42d0
      amb = 4.62d0
      amt = 175.d0
C      READ(NI,100)AMTAU
C      READ(NI,100)AMMUON
      amtau  = 1.7771d0
      ammuon = 0.105658389d0
C      READ(NI,100)ALPH
      alph = 137.0359895d0
C      READ(NI,100)GF
      gf = 1.16639d-5
C      READ(NI,100)GAMW
C      READ(NI,100)GAMZ
C      READ(NI,100)AMZ
C      READ(NI,100)AMW
      gamw = 2.062d0
      gamz = 2.485d0
      amz  = 91.190d0
      amw  = 79.941d0
C      READ(NI,100)VUS
C      READ(NI,100)VCB
C      READ(NI,100)RVUB
      vus  = 0.2205d0
      vcb  = 0.04d0
      rvub = 0.08d0
C      READ(NI,100)AMU
C      READ(NI,100)AM2
C      READ(NI,100)AMEL1
C      READ(NI,100)AMER1
C      READ(NI,100)AMQL1
C      READ(NI,100)AMUR1
C      READ(NI,100)AMDR1
C      READ(NI,100)AMEL
C      READ(NI,100)AMER
C      READ(NI,100)AMSQ
C      READ(NI,100)AMUR
C      READ(NI,100)AMDR
C      READ(NI,100)AL
C      READ(NI,100)AU
C      READ(NI,100)AD
C      READ(NI,101)NNLO
C      READ(NI,101)IONSH
C      READ(NI,101)IONWZ
C      READ(NI,101)IPOLE
C      READ(NI,101)IOFSUSY
C      READ(NI,101)INDIDEC
      indidec = 0
C      READ(NI,101)NFGG
      nfgg = 5

      VUB=RVUB*VCB
      ALPH=1.D0/ALPH
      AMSB = AMS

C--CHECK NFGG
      IF(NFGG.GT.5.OR.NFGG.LT.3)THEN
       WRITE(6,*)'NF-GG NOT VALID. TAKING THE DEFAULT NF-GG = 3....'
       NFGG = 3
      ENDIF

100   FORMAT(10X,G30.20)
101   FORMAT(10X,I30)

      goto 125

      IF(IHIGGS.EQ.0) THEN
       OPEN(NSA,FILE='br.sm1')
       OPEN(NSB,FILE='br.sm2')
      ENDIF
      IF(IHIGGS.EQ.1.OR.IHIGGS.EQ.5) THEN
       OPEN(NLA,FILE='br.l1')
       OPEN(NLB,FILE='br.l2')
      IF(IOFSUSY.EQ.0)THEN 
       OPEN(NSUSYL,FILE='br.ls')
       IF(INDIDEC.NE.0)THEN 
        OPEN(NSUSYLA,FILE='br.ls1')
        OPEN(NSUSYLB,FILE='br.ls2')
        OPEN(NSUSYLC,FILE='br.ls3')
        OPEN(NSUSYLD,FILE='br.ls4')
        OPEN(NSUSYLE,FILE='br.ls5')
        OPEN(NSUSYLF,FILE='br.ls6')
       ENDIF
      ENDIF
      ENDIF
      IF(IHIGGS.EQ.2.OR.IHIGGS.EQ.5) THEN
       OPEN(NHA,FILE='br.h1')
       OPEN(NHB,FILE='br.h2')
       OPEN(NHC,FILE='br.h3')
      IF(IOFSUSY.EQ.0)THEN 
       OPEN(NSUSYH,FILE='br.hs')
       IF(INDIDEC.NE.0)THEN 
        OPEN(NSUSYHA,FILE='br.hs1')
        OPEN(NSUSYHB,FILE='br.hs2')
        OPEN(NSUSYHC,FILE='br.hs3')
        OPEN(NSUSYHD,FILE='br.hs4')
        OPEN(NSUSYHE,FILE='br.hs5')
        OPEN(NSUSYHF,FILE='br.hs6')
       ENDIF
      ENDIF
      ENDIF
      IF(IHIGGS.EQ.3.OR.IHIGGS.EQ.5) THEN
       OPEN(NAA,FILE='br.a1')
       OPEN(NAB,FILE='br.a2')
      IF(IOFSUSY.EQ.0)THEN 
       OPEN(NSUSYA,FILE='br.as')
       IF(INDIDEC.NE.0)THEN 
        OPEN(NSUSYAA,FILE='br.as1')
        OPEN(NSUSYAB,FILE='br.as2')
        OPEN(NSUSYAC,FILE='br.as3')
        OPEN(NSUSYAD,FILE='br.as4')
       ENDIF
      ENDIF
      ENDIF
      IF(IHIGGS.EQ.4.OR.IHIGGS.EQ.5) THEN
       OPEN(NCA,FILE='br.c1')
       OPEN(NCB,FILE='br.c2')
      IF(IOFSUSY.EQ.0)THEN 
       OPEN(NSUSYC,FILE='br.cs')
       IF(INDIDEC.NE.0)THEN 
        OPEN(NSUSYCA,FILE='br.cs1')
        OPEN(NSUSYCB,FILE='br.cs2')
        OPEN(NSUSYCC,FILE='br.cs3')
        OPEN(NSUSYCD,FILE='br.cs4')
       ENDIF
      ENDIF
      ENDIF

 125  continue

      AMC0=AMC
      AMB0=AMB
      AMT0=AMT
      ACC=1.D-8
      NLOOP=2
      XLAMBDA=XITLA(NLOOP,ALSMZ,ACC)
      N0=5
      CALL ALSINI(ACC)
  
C--INITIALIZE COEFFICIENTS FOR POLYLOGARITHMS
      NBER = 18
      CALL BERNINI(NBER)

C--WRITE THE INPUT PARAMTERS TO A DATA-FILE

      goto 126

      WRITE(NPAR,8)'HIGGS    = ',IHIGGS
      WRITE(NPAR,9)'TGBET    = ',TGBET
      WRITE(NPAR,9)'MABEG    = ',AMABEG
      WRITE(NPAR,9)'MAEND    = ',AMAEND
      WRITE(NPAR,7)'NMA      = ',NMA
      WRITE(NPAR,9)'ALS(MZ)  = ',ALSMZ
      WRITE(NPAR,9)'MSBAR(1) = ',AMS
      WRITE(NPAR,9)'MC       = ',AMC
      WRITE(NPAR,9)'MB       = ',AMB
      WRITE(NPAR,9)'MT       = ',AMT
      WRITE(NPAR,9)'MTAU     = ',AMTAU
      WRITE(NPAR,9)'MMUON    = ',AMMUON
      WRITE(NPAR,9)'ALPH     = ',1.D0/ALPH
      WRITE(NPAR,9)'GF       = ',GF
      WRITE(NPAR,9)'GAMW     = ',GAMW
      WRITE(NPAR,9)'GAMZ     = ',GAMZ
      WRITE(NPAR,9)'MZ       = ',AMZ
      WRITE(NPAR,9)'MW       = ',AMW
      WRITE(NPAR,9)'VUS      = ',VUS
      WRITE(NPAR,9)'VCB      = ',VCB
      WRITE(NPAR,9)'VUB/VCB  = ',RVUB
      WRITE(NPAR,9)'MU       = ',AMU
      WRITE(NPAR,9)'M2       = ',AM2
      WRITE(NPAR,9)'MEL1      = ',AMEL1
      WRITE(NPAR,9)'MER1      = ',AMER1
      WRITE(NPAR,9)'MQL1      = ',AMQL1
      WRITE(NPAR,9)'MUR1      = ',AMUR1
      WRITE(NPAR,9)'MDR1      = ',AMDR1
      WRITE(NPAR,9)'MEL      = ',AMEL
      WRITE(NPAR,9)'MER      = ',AMER
      WRITE(NPAR,9)'MSQ      = ',AMSQ
      WRITE(NPAR,9)'MUR      = ',AMUR
      WRITE(NPAR,9)'MDR      = ',AMDR
      WRITE(NPAR,9)'AL       = ',AL
      WRITE(NPAR,9)'AU       = ',AU
      WRITE(NPAR,9)'AD       = ',AD
      WRITE(NPAR,8)'NNLO (M) = ',NNLO
      WRITE(NPAR,8)'ON-SHELL = ',IONSH
      WRITE(NPAR,8)'ON-SH-WZ = ',IONWZ
      WRITE(NPAR,8)'OFF-SUSY = ',IOFSUSY
      WRITE(NPAR,8)'IPOLE    = ',IPOLE 
      WRITE(NPAR,8)'NF-GG    = ',NFGG
      WRITE(NPAR,9)'LAMBDA_5 = ',XLAMBDA

      CLOSE(NPAR)

 126  continue

7     FORMAT(A11,I7)
8     FORMAT(A11,I4)
9     FORMAT(A11,G15.6)

C--SETUP THE HEADS OF THE TABLES IN THE DATA-FILES

      goto 127

      IF(IHIGGS.EQ.0) THEN
      WRITE(NSA,70)'MHSM  ','BB   ','TAU TAU','MU MU ','SS ','CC ','TT '
      WRITE(NSA,69)
      WRITE(NSA,*)
      WRITE(NSB,70)'MHSM  ','GG ','GAM GAM','Z GAM ','WW ','ZZ ','WIDTH'
      WRITE(NSB,69)
      WRITE(NSB,*)
      ENDIF

      IF(IHIGGS.EQ.1.OR.IHIGGS.EQ.5) THEN
      WRITE(NLA,70)'MHL   ','BB   ','TAU TAU','MU MU ','SS ','CC ','TT '
      WRITE(NLA,69)
      WRITE(NLA,*)
      WRITE(NLB,70)'MHL   ','GG ','GAM GAM','Z GAM ','WW ','ZZ ','WIDTH'
      WRITE(NLB,69)
      WRITE(NLB,*)
      ENDIF

      IF(IHIGGS.EQ.2.OR.IHIGGS.EQ.5) THEN
      WRITE(NHA,70)'MHH   ','BB   ','TAU TAU','MU MU ','SS ','CC ','TT '
      WRITE(NHA,69)
      WRITE(NHA,*)
      WRITE(NHB,72)'MHH   ','GG ','GAM GAM','Z GAM ','WW ','ZZ '
      WRITE(NHB,69)
      WRITE(NHB,*)
      WRITE(NHC,72)'MHH   ','hh ','AA ','Z A ','W+- H-+','WIDTH '
      WRITE(NHC,69)
      WRITE(NHC,*)
      ENDIF

      IF(IHIGGS.EQ.3.OR.IHIGGS.EQ.5) THEN
      WRITE(NAA,70)'MHA   ','BB   ','TAU TAU','MU MU ','SS ','CC ','TT '
      WRITE(NAA,69)
      WRITE(NAA,*)
      WRITE(NAB,72)'MHA   ','GG ','GAM GAM','Z GAM ','Z HL ','WIDTH '
      WRITE(NAB,69)
      WRITE(NAB,*)
      ENDIF

      IF(IHIGGS.EQ.4.OR.IHIGGS.EQ.5) THEN
      WRITE(NCA,70)'MHC   ','BC   ','TAU NU ','MU NU ','SU ','CS ','TB '
      WRITE(NCA,69)
      WRITE(NCA,*)
      WRITE(NCB,70)'MHC   ','hW ','AW ','WIDTH '
      WRITE(NCB,69)
      WRITE(NCB,*)
      ENDIF

 127  continue

69    FORMAT(79('_'))
70    FORMAT(A9,6(1X,A10))
71    FORMAT(A9,4(1X,A10))
72    FORMAT(A9,5(1X,A10))
73    FORMAT(A9,3(1X,A10))

      IWRISU = 1

C      DO 9999 II=1,NMA
C       IF(NMA.NE.1)THEN
C        AMAR = AMABEG + (AMAEND-AMABEG)/(NMA-1D0)*(II-1D0)
C       ELSE
        AMAR = AMABEG
C       ENDIF
       AMSM = AMAR
       AMA = AMAR

      IF(IHIGGS.NE.0)THEN 
C *******************************  SUSY OUTPUT 
      CALL SUSYCP(TGBET)

      IF(IOFSUSY.EQ.0.AND.IWRISU.NE.0)THEN
C--WRITE THE GAUGINO MASSES/ TB, MU AND M2 IN THE SUSY DATA-FILE
C--WRITE THE SFERMION MASSES/ SUSY MASSES AND COUPLINGS IN SUSY DATA-FILE
       CALL GAUGINO(AMU,AM2,B,A,GMC,GMN,XMN,AC1,AC2,AC3,
     .              AN1,AN2,AN3,ACNL,ACNR)         
       CALL SFERMION(AMSQ,AMUR,AMDR,AMEL,AMER,AL,AU,AD,AMU,
     .               GMST,GMSB,GMSL,GMSU,GMSD,GMSE,GMSN, 
     .               GLEE,GLTT,GLBB,GHEE,GHTT,GHBB,
     .               GAEE,GATT,GABB,GCEN,GCTB)

C 
       IF(IHIGGS.EQ.1.OR.IHIGGS.EQ.5) THEN
       WRITE(NSUSYL,347) TGBET,AM2,AMU,AMSQ
       WRITE(NSUSYL,348) GMC(1),GMC(2),GMN(1),GMN(2),GMN(3),GMN(4)
       WRITE(NSUSYL,349) GMST(1),GMST(2),GMSU(1),GMSU(2)
       WRITE(NSUSYL,350) GMSB(1),GMSB(2),GMSD(1),GMSD(2)
       WRITE(NSUSYL,351) GMSL(1),GMSL(2),GMSE(1),GMSE(2),GMSN(1)
       WRITE(NSUSYL,*)
       WRITE(NSUSYL,*)'   MHL        CHARGINOS  NEUTRALS   '//
     . 'SLEPTONS   SQUARKS'
       WRITE(NSUSYL,69)
       WRITE(NSUSYL,*)
        IF(INDIDEC.NE.0)THEN
         WRITE(NSUSYLA,73)'MHL   ','C1 C1 ','C2 C2 ','C1 C2 '
         WRITE(NSUSYLA,69)
         WRITE(NSUSYLA,*)
         WRITE(NSUSYLB,71)'MHL   ','N1 N1 ','N2 N2 ','N3 N3 ','N4 N4 '
         WRITE(NSUSYLB,69)
         WRITE(NSUSYLB,*)
         WRITE(NSUSYLC,70)'MHL   ','N1 N2 ','N1 N3 ','N1 N4 ','N2 N3 ',
     .                    'N2 N4 ','N3 N4 '
         WRITE(NSUSYLC,69)
         WRITE(NSUSYLC,*)
         WRITE(NSUSYLD,*)'   MHL        SNL SNL    SEL SEL    '//
     .   'SER SER    STA1 STA1  STA1 STA2  STA2 STA2' 
         WRITE(NSUSYLD,69)
         WRITE(NSUSYLD,*)
         WRITE(NSUSYLE,*)'   MHL        SUL SUL    SUR SUR    '//
     .   'SDL SDL    SDR SDR'
         WRITE(NSUSYLE,69)
         WRITE(NSUSYLE,*)
         WRITE(NSUSYLF,*)'   MHL        SB1 SB1    SB1 SB2    '//
     .   'SB2 SB2    ST1 ST1    ST1 ST2    ST2 ST2'
         WRITE(NSUSYLF,69)
         WRITE(NSUSYLF,*)
        ENDIF
       ENDIF

       IF(IHIGGS.EQ.2.OR.IHIGGS.EQ.5) THEN
       WRITE(NSUSYH,347) TGBET,AM2,AMU,AMSQ
       WRITE(NSUSYH,348) GMC(1),GMC(2),GMN(1),GMN(2),GMN(3),GMN(4)
       WRITE(NSUSYH,349) GMST(1),GMST(2),GMSU(1),GMSU(2)
       WRITE(NSUSYH,350) GMSB(1),GMSB(2),GMSD(1),GMSD(2)
       WRITE(NSUSYH,351) GMSL(1),GMSL(2),GMSE(1),GMSE(2),GMSN(1)
       WRITE(NSUSYH,*)
       WRITE(NSUSYH,*)'   MHH        CHARGINOS  NEUTRALS   '//
     . 'SLEPTONS   SQUARKS'
       WRITE(NSUSYH,69)
       WRITE(NSUSYH,*)
        IF(INDIDEC.NE.0)THEN
         WRITE(NSUSYHA,73)'MHH   ','C1 C1 ','C2 C2 ','C1 C2 '
         WRITE(NSUSYHA,69)
         WRITE(NSUSYHA,*)
         WRITE(NSUSYHB,71)'MHH   ','N1 N1 ','N2 N2 ','N3 N3 ','N4 N4 '
         WRITE(NSUSYHB,69)
         WRITE(NSUSYHB,*)
         WRITE(NSUSYHC,70)'MHH   ','N1 N2 ','N1 N3 ','N1 N4 ','N2 N3 ',
     .                    'N2 N4 ','N3 N4 '
         WRITE(NSUSYHC,69)
         WRITE(NSUSYHC,*)
         WRITE(NSUSYHD,*)'   MHH        SNL SNL    SEL SEL    '//
     .   'SER SER    STA1 STA1  STA1 STA2  STA2 STA2' 
         WRITE(NSUSYHD,69)
         WRITE(NSUSYHD,*)
         WRITE(NSUSYHE,*)'   MHH        SUL SUL    SUR SUR    '//
     .   'SDL SDL    SDR SDR'
         WRITE(NSUSYHE,69)
         WRITE(NSUSYHE,*)
         WRITE(NSUSYHF,*)'   MHH        SB1 SB1    SB1 SB2    '//
     .   'SB2 SB2    ST1 ST1    ST1 ST2    ST2 ST2'
         WRITE(NSUSYHF,69)
         WRITE(NSUSYHF,*)
        ENDIF
       ENDIF

       IF(IHIGGS.EQ.3.OR.IHIGGS.EQ.5) THEN
       WRITE(NSUSYA,347) TGBET,AM2,AMU,AMSQ
       WRITE(NSUSYA,348) GMC(1),GMC(2),GMN(1),GMN(2),GMN(3),GMN(4)
       WRITE(NSUSYA,349) GMST(1),GMST(2),GMSU(1),GMSU(2)
       WRITE(NSUSYA,350) GMSB(1),GMSB(2),GMSD(1),GMSD(2)
       WRITE(NSUSYA,351) GMSL(1),GMSL(2),GMSE(1),GMSE(2),GMSN(1)
       WRITE(NSUSYA,*)
       WRITE(NSUSYA,*)'   MHH        CHARGINOS  NEUTRALS   '//
     . 'SLEPTONS   SQUARKS'
       WRITE(NSUSYA,69)
       WRITE(NSUSYA,*)
        IF(INDIDEC.NE.0)THEN
         WRITE(NSUSYAA,73)'MHA   ','C1 C1 ','C2 C2 ','C1 C2 '
         WRITE(NSUSYAA,69)
         WRITE(NSUSYAA,*)
         WRITE(NSUSYAB,71)'MHA   ','N1 N1 ','N2 N2 ','N3 N3 ','N4 N4 '
         WRITE(NSUSYAB,69)
         WRITE(NSUSYAB,*)
         WRITE(NSUSYAC,70)'MHA   ','N1 N2 ','N1 N3 ','N1 N4 ','N2 N3 ',
     .                    'N2 N4 ','N3 N4 '
         WRITE(NSUSYAC,69)
         WRITE(NSUSYAC,*)
         WRITE(NSUSYAD,*)
         WRITE(NSUSYAD,*)'   MHA        STA1 STA2  SB1 SB2    ST1 ST2'
         WRITE(NSUSYAD,69)
         WRITE(NSUSYAD,*)
        ENDIF
       ENDIF

       IF(IHIGGS.EQ.4.OR.IHIGGS.EQ.5) THEN
       WRITE(NSUSYC,347) TGBET,AM2,AMU,AMSQ
       WRITE(NSUSYC,348) GMC(1),GMC(2),GMN(1),GMN(2),GMN(3),GMN(4)
       WRITE(NSUSYC,349) GMST(1),GMST(2),GMSU(1),GMSU(2)
       WRITE(NSUSYC,350) GMSB(1),GMSB(2),GMSD(1),GMSD(2)
       WRITE(NSUSYC,351) GMSL(1),GMSL(2),GMSE(1),GMSE(2),GMSN(1)
       WRITE(NSUSYC,*)
       WRITE(NSUSYC,*)'   MHC        CHARG/NEU  SLEPTONS   SQUARKS'
       WRITE(NSUSYC,69)
       WRITE(NSUSYC,*)
        IF(INDIDEC.NE.0)THEN
         WRITE(NSUSYCA,70)'MHC   ','C1 N1 ','C1 N2 ','C1 N3 ','C1 N4 '
         WRITE(NSUSYCA,69)
         WRITE(NSUSYCA,*)
         WRITE(NSUSYCB,70)'MHC   ','C2 N1 ','C2 N2 ','C2 N3 ','C2 N4 '
         WRITE(NSUSYCB,69)
         WRITE(NSUSYCB,*)
         WRITE(NSUSYCC,*)'   MHC        SEL SNL    STAU1 SNL  STAU2 SNL'
         WRITE(NSUSYCC,69)
         WRITE(NSUSYCC,*)
         WRITE(NSUSYCD,*)'   MHC        SUL SDL    ST1 SB1    '//
     .   'ST1 SB2    ST2 SB1    ST2 SB2'
         WRITE(NSUSYCD,69)
         WRITE(NSUSYCD,*)
        ENDIF
       ENDIF

347    FORMAT('TB=',G12.6,1X,'M2=',G12.6,1X,'MU=',G12.6,1X,
     .        'MSQ=',G12.6)
348    FORMAT('C1=',F7.3,1X,'C2=',F8.3,1X,'N1=',F7.3,1X,'N2=',F7.3,1X,
     .        'N3=',F8.3,1X,'N4=',F8.3)
349    FORMAT('MST1=',G12.6,1X,'MST2=',G12.6,1X,
     .        'MSUL=',G12.6,1X,'MSUR=',G12.6) 
350    FORMAT('MSB1=',G12.6,1X,'MSB2=',G12.6,1X,
     .        'MSDL=',G12.6,1X,'MSDR=',G12.6) 
351    FORMAT('TAU1=',F8.3,1X,'TAU2=',F8.3,1X,'EL=',F8.3,1X,
     .        'ER=',F8.3,1X,'NL=',F8.3)
C
C
C **************************************************************
      IWRISU = 0
      ENDIF
      ENDIF

      CALL HDEC(TGBET)

C      IF(IHIGGS.EQ.0)THEN
C      WRITE(NSA,20)AMSM,SMBRB,SMBRL,SMBRM,SMBRS,SMBRC,SMBRT
C      WRITE(NSB,20)AMSM,SMBRG,SMBRGA,SMBRZGA,SMBRW,SMBRZ,SMWDTH
C      ENDIF

      IF(IHIGGS.EQ.1.OR.IHIGGS.EQ.5)THEN
      WRITE(NLA,20)AML,HLBRB,HLBRL,HLBRM,HLBRS,HLBRC,HLBRT
      WRITE(NLB,20)AML,HLBRG,HLBRGA,HLBRZGA,HLBRW,HLBRZ,HLWDTH
c      IF(IOFSUSY.EQ.0)THEN 
c       WRITE(NSUSYL,21)AML,HLBRCHT,HLBRNET,HLBRSL,HLBRSQT
c       IF(INDIDEC.NE.0)THEN
c        WRITE(NSUSYLA,23)AML,HLBRSC(1,1),HLBRSC(2,2),
c     .                   HLBRSC(1,2)+HLBRSC(2,1)
c        WRITE(NSUSYLB,21)AML,HLBRSN(1,1),HLBRSN(2,2),HLBRSN(3,3),
c     .                   HLBRSN(4,4)
c        WRITE(NSUSYLC,20)AML,HLBRSN(1,2)+HLBRSN(2,1),
c     .                   HLBRSN(1,3)+HLBRSN(3,1),
c     .                   HLBRSN(1,4)+HLBRSN(4,1),
c     .                   HLBRSN(2,3)+HLBRSN(3,2),
c     .                   HLBRSN(2,4)+HLBRSN(4,2),
c     .                   HLBRSN(3,4)+HLBRSN(4,3)
c        WRITE(NSUSYLD,20)AML,BHLSLNL,BHLSLEL,BHLSLER,BHLSTAU(1,1),
c     .                   BHLSTAU(1,2)+BHLSTAU(2,1),BHLSTAU(2,2)
c        WRITE(NSUSYLE,21)AML,BHLSQUL,BHLSQUR,BHLSQDL,BHLSQDR
c      WRITE(NSUSYLF,20)AML,BHLSB(1,1),BHLSB(1,2)+BHLSB(2,1),BHLSB(2,2),
c     .                   BHLST(1,1),BHLST(1,2)+BHLST(2,1),BHLST(2,2)
c       ENDIF 
c      ENDIF
      ENDIF

      IF(IHIGGS.EQ.2.OR.IHIGGS.EQ.5)THEN
      WRITE(NHA,20)AMH,HHBRB,HHBRL,HHBRM,HHBRS,HHBRC,HHBRT
      WRITE(NHB,20)AMH,HHBRG,HHBRGA,HHBRZGA,HHBRW,HHBRZ
      WRITE(NHC,20)AMH,HHBRH,HHBRA,HHBRAZ,HHBRHW,HHWDTH
c      IF(IOFSUSY.EQ.0)THEN 
c       WRITE(NSUSYH,21)AMH,HHBRCHT,HHBRNET,HHBRSL,HHBRSQT
c       IF(INDIDEC.NE.0)THEN
c        WRITE(NSUSYHA,23)AMH,HHBRSC(1,1),HHBRSC(2,2),
c     .                  HHBRSC(1,2)+HHBRSC(2,1)
c        WRITE(NSUSYHB,21)AMH,HHBRSN(1,1),HHBRSN(2,2),HHBRSN(3,3),
c     .                   HHBRSN(4,4)
c        WRITE(NSUSYHC,20)AMH,HHBRSN(1,2)+HHBRSN(2,1),
c     .                   HHBRSN(1,3)+HHBRSN(3,1),
c     .                   HHBRSN(1,4)+HHBRSN(4,1),
c     .                   HHBRSN(2,3)+HHBRSN(3,2),
c     .                   HHBRSN(2,4)+HHBRSN(4,2),
c     .                   HHBRSN(3,4)+HHBRSN(4,3)
c        WRITE(NSUSYHD,20)AMH,BHHSLNL,BHHSLEL,BHHSLER,BHHSTAU(1,1),
c     .                   BHHSTAU(1,2)+BHHSTAU(2,1),BHHSTAU(2,2)
c        WRITE(NSUSYHE,21)AMH,BHHSQUL,BHHSQUR,BHHSQDL,BHHSQDR
c      WRITE(NSUSYHF,20)AMH,BHHSB(1,1),BHHSB(1,2)+BHHSB(2,1),BHHSB(2,2),
c     .                   BHHST(1,1),BHHST(1,2)+BHHST(2,1),BHHST(2,2)
c       ENDIF
c      ENDIF
      ENDIF

      IF(IHIGGS.EQ.3.OR.IHIGGS.EQ.5)THEN
      WRITE(NAA,20)AMA,ABRB,ABRL,ABRM,ABRS,ABRC,ABRT
      WRITE(NAB,22)AMA,ABRG,ABRGA,ABRZGA,ABRZ,AWDTH
c      IF(IOFSUSY.EQ.0)THEN 
c       WRITE(NSUSYA,21)AMA,HABRCHT,HABRNET,HABRSL,HABRST+HABRSB 
c       IF(INDIDEC.NE.0)THEN
c        WRITE(NSUSYAA,23)AMA,HABRSC(1,1),HABRSC(2,2),
c     .                   HABRSC(1,2)+HABRSC(2,1)
c        WRITE(NSUSYAB,21)AMA,HABRSN(1,1),HABRSN(2,2),HABRSN(3,3),
c     .                   HABRSN(4,4)
c        WRITE(NSUSYAC,20)AMA,HABRSN(1,2)+HABRSN(2,1),
c     .                   HABRSN(1,3)+HABRSN(3,1),
c     .                   HABRSN(1,4)+HABRSN(4,1),
c     .                   HABRSN(2,3)+HABRSN(3,2),
c     .                   HABRSN(2,4)+HABRSN(4,2),
c     .                   HABRSN(3,4)+HABRSN(4,3)
c        WRITE(NSUSYAD,23)AMA,BHASTAU,BHASB,BHAST
c       ENDIF
c      ENDIF
      ENDIF

      IF(IHIGGS.EQ.4.OR.IHIGGS.EQ.5)THEN
      WRITE(NCA,20)AMCH,HCBRB,HCBRL,HCBRM,HCBRS,HCBRC,HCBRT
      WRITE(NCB,23)AMCH,HCBRW,HCBRA,HCWDTH
c      IF(IOFSUSY.EQ.0)THEN 
c       WRITE(NSUSYC,23)AMCH,HCBRCNT,HCBRSL,HCBRSQT 
c       IF(INDIDEC.NE.0)THEN
c        WRITE(NSUSYCA,21)AMCH,HCBRSU(1,1),HCBRSU(1,2),
c     .                   HCBRSU(1,3),HCBRSU(1,4)
c        WRITE(NSUSYCB,21)AMCH,HCBRSU(2,1),HCBRSU(2,2),
c     .                   HCBRSU(2,3),HCBRSU(2,4)
c        WRITE(NSUSYCC,23)AMCH,BHCSL00,BHCSL11,BHCSL21
c        WRITE(NSUSYCD,22)AMCH,BHCSQ,BHCSTB(1,1),BHCSTB(1,2),
c     .                   BHCSTB(2,1),BHCSTB(2,2)
c       ENDIF
c      ENDIF
      ENDIF
  
20    FORMAT(G12.6,6(1X,G10.4))
21    FORMAT(G12.6,4(1X,G10.4))
22    FORMAT(G12.6,5(1X,G10.4))
23    FORMAT(G12.6,3(1X,G10.4))

C9999  CONTINUE
c      CLOSE(NI)

      IF(IHIGGS.EQ.0) THEN
c       CLOSE(NSA)
c       CLOSE(NSB)
      ENDIF

      IF(IHIGGS.EQ.1.OR.IHIGGS.EQ.5) THEN
c       CLOSE(NLA)
c       CLOSE(NLB) 
c       CLOSE(NSUSYL)
      ENDIF

      IF(IHIGGS.EQ.2.OR.IHIGGS.EQ.5) THEN
       CLOSE(NHA)
       CLOSE(NHB) 
       CLOSE(NHC)
       CLOSE(NSUSYH)
      ENDIF

      IF(IHIGGS.EQ.3.OR.IHIGGS.EQ.5) THEN
c       CLOSE(NAA)
c       CLOSE(NAB) 
c       CLOSE(NSUSYA)
      ENDIF

      IF(IHIGGS.EQ.4.OR.IHIGGS.EQ.5) THEN
c       CLOSE(NCA)
c       CLOSE(NCB) 
c       CLOSE(NSUSYC)
      ENDIF

      return
      end

C =====================================================================
C =========== BEGINNING OF THE SUBROUTINE FOR THE DECAYS ==============
C !!!!!!!!!!!!!! Any change below this line is at your own risk!!!!!!!!
C =====================================================================

      SUBROUTINE HDEC(TGBET)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION LAMB
      DIMENSION XX(4),YY(4)
      DIMENSION AMCHAR(2),AMNEUT(4),XMNEUT(4),
     .          AC1(2,2),AC2(2,2),AC3(2,2),
     .          AN1(4,4),AN2(4,4),AN3(4,4),
     .          ACNL(2,4),ACNR(2,4),
     .          AMST(2),AMSB(2),AMSL(2),
     .          AMSU(2),AMSD(2),AMSE(2),AMSN(2),
     .          GLTT(2,2),GLBB(2,2),GLEE(2,2),
     .          GHTT(2,2),GHBB(2,2),GHEE(2,2),
     .          GCTB(2,2),GCEN(2,2)
      DIMENSION GMST(2),GMSB(2),GMSL(2),GMSU(2),GMSD(2),GMSE(2),
     .          GMSN(2)
      DIMENSION HLBRSC(2,2),HLBRSN(4,4),HHBRSC(2,2),
     .          HHBRSN(4,4),HABRSC(2,2),HABRSN(4,4),HCBRSU(2,4),
     .          HHBRST(2,2),HHBRSB(2,2),HCBRSTB(2,2) 
      DIMENSION WHLCH(2,2),WHLNE(4,4),WHHCH(2,2),WHHNE(4,4),
     .          WHACH(2,2),WHANE(4,4),WHCCN(2,4),
     .          WHHST(2,2),WHHSB(2,2),WHHSTAU(2,2),WHCSTB(2,2), 
     .          WHLST(2,2),WHLSB(2,2),WHLSTAU(2,2)
      COMPLEX*16 CF,CG,CI1,CI2,CA,CB,CTT,CTB,CTC,CTW,CLT,CLB,CLW,
     .           CAT,CAB,CAC,CAW,CAH,CTH,CLH,CX1,CX2,CAX1,CAX2,CTL,CAL,
     .           CSL,CSQ,CSB1,CSB2,CST1,CST2,CSL1,CSL2,
     .           CXL,CXQ,CXB1,CXB2,CXT1,CXT2,CXL1,CXL2
      COMPLEX*16 CSEL,CSER,CSUL,CSUR,CSDL,CSDR,
     .           CXEL,CXER,CXUL,CXUR,CXDL,CXDR
      COMMON/HMASS/AMSM,AMA,AML,AMH,AMCH,AMAR
      COMMON/CHIMASS/AMCHI
      COMMON/MASSES/AMS,AMC,AMB,AMT
      COMMON/ALS/XLAMBDA,AMC0,AMB0,AMT0,N0
      COMMON/PARAM/GF,ALPH,AMTAU,AMMUON,AMZ,AMW
      COMMON/CKMPAR/VUS,VCB,VUB
      COMMON/BREAK/AMEL,AMER,AMSQ,AMUR,AMDR,AL,AU,AD,AMU,AM2
      COMMON/WZWDTH/GAMC0,GAMT0,GAMT1,GAMW,GAMZ
      COMMON/ONSHELL/IONSH,IONWZ,IOFSUSY
      COMMON/OLDFASH/NFGG
      COMMON/FLAG/IHIGGS,NNLO,IPOLE
      COMMON/WIDTHSM/SMBRB,SMBRL,SMBRM,SMBRS,SMBRC,SMBRT,SMBRG,SMBRGA,
     .               SMBRZGA,SMBRW,SMBRZ,SMWDTH
      COMMON/WIDTHA/ABRB,ABRL,ABRM,ABRS,ABRC,ABRT,ABRG,ABRGA,ABRZGA,
     .              ABRZ,AWDTH
      COMMON/WIDTHHL/HLBRB,HLBRL,HLBRM,HLBRS,HLBRC,HLBRT,HLBRG,HLBRGA,
     .               HLBRZGA,HLBRW,HLBRZ,HLBRA,HLBRAZ,HLBRHW,HLWDTH
      COMMON/WIDTHHH/HHBRB,HHBRL,HHBRM,HHBRS,HHBRC,HHBRT,HHBRG,HHBRGA,
     .               HHBRZGA,HHBRW,HHBRZ,HHBRH,HHBRA,HHBRAZ,HHBRHW,
     .               HHWDTH
      COMMON/WIDTHHC/HCBRB,HCBRL,HCBRM,HCBRBU,HCBRS,HCBRC,HCBRT,HCBRW,
     .               HCBRA,HCWDTH
      COMMON/WISUSY/HLBRSC,HLBRSN,HHBRSC,HHBRSN,HABRSC,HABRSN,HCBRSU,
     .              HLBRCHT,HHBRCHT,HABRCHT,HLBRNET,HHBRNET,HABRNET,
     .              HCBRCNT,HLBRSL,HHBRSL,HCBRSL,HABRSL,HABRST,HABRSB,
     .              HHBRSQ,HHBRST,HHBRSB,HHBRSQT,HCBRSQ,HCBRSTB,
     .              HCBRSQT,HLBRSQ,HLBRSQT
      COMMON/WISFER/BHLSLNL,BHLSLEL,BHLSLER,BHLSQUL,BHLSQUR,BHLSQDL,
     .              BHLSQDR,BHLST(2,2),BHLSB(2,2),BHLSTAU(2,2),
     .              BHHSLNL,BHHSLEL,BHHSLER,BHHSQUL,BHHSQUR,BHHSQDL,
     .              BHHSQDR,BHHST(2,2),BHHSB(2,2),BHHSTAU(2,2),
     .              BHASTAU,BHASB,BHAST,
     .              BHCSL00,BHCSL11,BHCSL21,BHCSQ,BHCSTB(2,2)
      COMMON/SMASS/AMNEUT,XMNEUT,AMCHAR,AMST,AMSB,AMSL,
     .              AMSU,AMSD,AMSE,AMSN 
      COMMON/COUP/GAT,GAB,GLT,GLB,GHT,GHB,GZAH,GZAL,
     .            GHHH,GLLL,GHLL,GLHH,GHAA,GLAA,GLVV,GHVV,
     .            GLPM,GHPM,B,A
C      write(6,*) 'entering HDEC'
      HVV(X,Y)= GF/(4.D0*PI*DSQRT(2.D0))*X**3/2.D0*BETA(Y)
     .            *(1.D0-4.D0*Y+12.D0*Y**2)
C      write(6,*) 'calculated HVV'
      AFF(X,Y)= GF/(4*PI*DSQRT(2.D0))*X**3*Y*(BETA(Y))
      HFF(X,Y)= GF/(4*PI*DSQRT(2.D0))*X**3*Y*(BETA(Y))**3
      CFF(Z,TB,X,Y)= GF/(4*PI*DSQRT(2.D0))*Z**3*LAMB(X,Y)
     .              *((1.D0-X-Y)*(X*TB**2+Y/TB**2)-4.D0*X*Y)
C      write(6,*) 'calculated xFF'
      HV(V)=3.D0*(1.D0-8.D0*V+20.D0*V**2)/DSQRT((4.D0*V-1.D0))
     .      *DACOS((3.D0*V-1.D0)/2.D0/DSQRT(V**3))
     .      -(1.D0-V)*(47.D0/2.D0*V-13.D0/2.D0+1.D0/V)
     .      -3.D0/2.D0*(1.D0-6.D0*V+4.D0*V**2)*DLOG(V)
C      write(6,*) 'calculated HV'
      HVH(X,Y)=0.25D0*( (1-X)*(-2+4*X-2*X**2+9*Y+9*X*Y-6*Y**2)
     .        /(3*Y)-2*(1-X-X**2+X**3-3*Y-2*X*Y-3*X**2*Y+3*Y**2
     .        +3*X*Y**2-Y**3)*(-PI/2- DATAN((1-2*X+X**2-Y-X*Y)/
     .         ((1-X)*DSQRT(-1.D0+2*X+2*Y-(X-Y)**2))))/DSQRT(-1.D0
     .         +2*X-(X-Y)**2+2*Y)-(1+X**2-2*Y-2*X*Y+Y**2)*DLOG(X))
C      write(6,*) 'calculated HVH'
      QCD0(X) = (1+X**2)*(4*SP((1-X)/(1+X)) + 2*SP((X-1)/(X+1))
     .        - 3*DLOG((1+X)/(1-X))*DLOG(2/(1+X))
     .        - 2*DLOG((1+X)/(1-X))*DLOG(X))
     .        - 3*X*DLOG(4/(1-X**2)) - 4*X*DLOG(X)
C      write(6,*) 'calculated QCD0'
      HQCDM(X)=QCD0(X)/X+(3+34*X**2-13*X**4)/16/X**3*DLOG((1+X)/(1-X))
     .        + 3.D0/8/X**2*(7*X**2-1)
C      write(6,*) 'calculated HQCDM'
      AQCDM(X)=QCD0(X)/X + (19+2*X**2+3*X**4)/16/X*DLOG((1+X)/(1-X))
     .        + 3.D0/8*(7-X**2)
C      write(6,*) 'calculated AQCDM'
      HQCD(X)=(4.D0/3*HQCDM(BETA(X))
     .        +2*(4.D0/3-DLOG(X))*(1-10*X)/(1-4*X))*ASH/PI
     .       + (29.14671D0 + RATCOUP*(1.570D0 - 2*DLOG(HIGTOP)/3
     .                                     + DLOG(X)**2/9))*(ASH/PI)**2
     .       + (164.14D0 - 25.77D0*5 + 0.259D0*5**2)*(ASH/PI)**3
C      write(6,*) 'calculated HQCD'
      AQCD(X)=(4.D0/3*AQCDM(BETA(X))
     .        +2*(4.D0/3-DLOG(X))*(1-6*X)/(1-4*X))*ASH/PI
     .       + (29.14671D0 + RATCOUP*(23/6.D0 - DLOG(HIGTOP)
     .                                     + DLOG(X)**2/6))*(ASH/PI)**2
     .       + (164.14D0 - 25.77D0*5 + 0.259D0*5**2)*(ASH/PI)**3
C      write(6,*) 'calculated AQCD'
      QCDH(X)=1.D0+HQCD(X)
      QCDA(X)=1.D0+AQCD(X)
      TQCDH(X)=1.D0+4.D0/3*HQCDM(BETA(X))*ASH/PI
      TQCDA(X)=1.D0+4.D0/3*AQCDM(BETA(X))*ASH/PI
C      write(6,*) 'calculated TQCDx'
      QCDC(X,Y)=1.D0+4/3.D0*ASH/PI*(9/4.D0 + (3-2*X+2*Y)/4*DLOG(X/Y)
     .         +((1.5D0-X-Y)*LAMB(X,Y)**2+5*X*Y)/2/LAMB(X,Y)
     .         /(1-X-Y)*DLOG(XI(X,Y)*XI(Y,X))
     .         + BIJ(X,Y))
     .         + ASH/PI*(2*(4/3.D0-DLOG(X))
     .         - (X*2*(4/3.D0-DLOG(X)) + Y*2*(4/3.D0-DLOG(Y)))/(1-X-Y)
     .         - (X*2*(4/3.D0-DLOG(X))*(1-X+Y)
     .           +Y*2*(4/3.D0-DLOG(Y))*(1+X-Y))/LAMB(X,Y)**2)
      QCDCI(X,Y)=1.D0+4/3.D0*ASH/PI*(3 + (Y-X)/2*DLOG(X/Y)
     .         +(2*(1-X-Y)+LAMB(X,Y)**2)/2/LAMB(X,Y)
     .         *DLOG(XI(X,Y)*XI(Y,X))
     .         + BIJ(X,Y))
     .         + ASH/PI*(2*(4/3.D0-DLOG(X)) + 2*(4/3.D0-DLOG(Y))
     .         - (X*2*(4/3.D0-DLOG(X))*(1-X+Y)
     .           +Y*2*(4/3.D0-DLOG(Y))*(1+X-Y))/LAMB(X,Y)**2)
      QCDCM(X,Y)=1.D0+4/3.D0*ASH/PI*(9/4.D0 + (3-2*X+2*Y)/4*DLOG(X/Y)
     .         +((1.5D0-X-Y)*LAMB(X,Y)**2+5*X*Y)/2/LAMB(X,Y)
     .         /(1-X-Y)*DLOG(4*X*Y/(1-X-Y+LAMB(X,Y))**2)
     .         + BIJ(X,Y))
      QCDCMI(X,Y)=1.D0+4/3.D0*ASH/PI*(3 + (Y-X)/2*DLOG(X/Y)
     .         +(2*(1-X-Y)*LAMB(X,Y)**2)/2/LAMB(X,Y)
     .         *DLOG(4*X*Y/(1-X-Y+LAMB(X,Y))**2)
     .         + BIJ(X,Y))
      CQCD(Z,TB,X,Y)= GF/(4*PI*DSQRT(2.D0))*Z**3*LAMB(X,Y)
     .              *((1.D0-X-Y)*(X*TB**2*QCDC(X,Y)
     .                           +Y/TB**2*QCDC(Y,X))
     .               -4.D0*X*Y*QCDCI(X,Y))
      CQCDM(Z,TB,X,Y)= GF/(4*PI*DSQRT(2.D0))*Z**3*LAMB(X,Y)
     .              *((1.D0-X-Y)*(X*TB**2*QCDCM(X,Y)
     .                           +Y/TB**2*QCDCM(Y,X))
     .               -4.D0*X*Y*QCDCMI(X,Y))
      ELW(AMH,AMF,QF,ACF)=ALPH/PI*3.D0/2*QF**2
     .                              *(3.D0/2-DLOG(AMH**2/AMF**2))
     .      +GF/8/DSQRT(2.D0)/PI**2*(ACF*AMT**2
     .        +AMW**2*(3*DLOG(CS)/SS-5)+AMZ**2*(0.5D0
     .          -3*(1-4*SS*DABS(QF))**2))
      CF(CA) = -CDLOG(-(1+CDSQRT(1-CA))/(1-CDSQRT(1-CA)))**2/4
      CG(CA) = CDSQRT(1-CA)/2*CDLOG(-(1+CDSQRT(1-CA))/(1-CDSQRT(1-CA)))
      CI1(CA,CB) = CA*CB/2/(CA-CB)
     .           + CA**2*CB**2/2/(CA-CB)**2*(CF(CA)-CF(CB))
     .           + CA**2*CB/(CA-CB)**2*(CG(CA)-CG(CB))
      CI2(CA,CB) = -CA*CB/2/(CA-CB)*(CF(CA)-CF(CB))
      HGGQCD(ASG,NF)=1.D0+ASG/PI*(95.D0/4.D0-NF*7.D0/6.D0)
      SGGQCD(ASG)=ASG/PI*17.D0/6.D0
      AGGQCD(ASG,NF)=1.D0+ASG/PI*(97.D0/4.D0-NF*7.D0/6.D0)
      HFFSELF(AMH)=1.D0+GF*AMH**2/16.D0/PI**2/DSQRT(2.D0)*2.117203D0
     .            -(GF*AMH**2/16.D0/PI**2/DSQRT(2.D0))**2*32.6567D0
      HVVSELF(AMH)=1.D0+GF*AMH**2/16.D0/PI**2/DSQRT(2.D0)*2.800952D0
     .            +(GF*AMH**2/16.D0/PI**2/DSQRT(2.D0))**2*62.0308D0

      PI=4D0*DATAN(1D0)
      SS=1.D0-(AMW/AMZ)**2
      CS=1.D0-SS

C--DECOUPLING THE TOP QUARK FROM ALPHAS
      AMT0=3.D8

C--TOP QUARK DECAY WIDTH
      GAMT0 = GF*AMT**3/8/DSQRT(2D0)/PI*(1-AMW**2/AMT**2)**2
     .                                 *(1+2*AMW**2/AMT**2)
      IF(IHIGGS.NE.0.AND.AMT.GT.AMCH+AMB)THEN
       GAMT1 = GF*AMT**3/8/DSQRT(2D0)/PI*(1-AMCH**2/AMT**2)**2
     .        *((AMB/AMT)**2*TGBET**2 + 1/TGBET**2)
      ELSE
       GAMT1 = 0
      ENDIF
      GAMT1 = GAMT0+GAMT1

      IF(IHIGGS.EQ.0)THEN

C        =========================================================
C                              SM HIGGS DECAYS
C        =========================================================
      AMXX=AMH
      AMH=AMSM
C     =============  RUNNING MASSES 
      RMS = RUNM(AMH,3)
      RMC = RUNM(AMH,4)
      RMB = RUNM(AMH,5)
      RMT = RUNM(AMH,6)
      RATCOUP = 1
      HIGTOP = AMH**2/AMT**2

      ASH=ALPHAS(AMH,2)
      AMC0=1.D8
      AMB0=2.D8
      AS3=ALPHAS(AMH,2)
      AMC0=AMC
      AS4=ALPHAS(AMH,2)
      AMB0=AMB
C     AMT0=AMT
C     =============== PARTIAL WIDTHS 
C  H ---> G G
C
       EPS=1.D-8
       NFEXT = 3
       ASG = AS3
       CTT = 4*AMT**2/AMH**2*DCMPLX(1D0,-EPS)
       CTB = 4*AMB**2/AMH**2*DCMPLX(1D0,-EPS)
       CAT = 2*CTT*(1+(1-CTT)*CF(CTT))
       CAB = 2*CTB*(1+(1-CTB)*CF(CTB))
       FQCD=HGGQCD(ASG,NFEXT)
       XFAC = CDABS(CAT+CAB)**2*FQCD
       HGG=HVV(AMH,0.D0)*(ASG/PI)**2*XFAC/8

C  H ---> G G* ---> G CC   TO BE ADDED TO H ---> CC
       NFEXT = 4
       ASG = AS4
       FQCD=HGGQCD(ASG,NFEXT)
       XFAC = CDABS(CAT+CAB)**2*FQCD
       DCC=HVV(AMH,0.D0)*(ASG/PI)**2*XFAC/8 - HGG
C  H ---> G G* ---> G BB   TO BE ADDED TO H ---> BB
       NFEXT = 5
       ASG = ASH
       FQCD=HGGQCD(ASG,NFEXT)
       XFAC = CDABS(CAT+CAB)**2*FQCD
       DBB=HVV(AMH,0.D0)*(ASG/PI)**2*XFAC/8 - HGG - DCC

      IF(NFGG.EQ.5)THEN
       HGG = HGG + DBB + DCC
       DBB = 0
       DCC = 0
      ELSEIF(NFGG.EQ.4)THEN
       HGG = HGG + DCC
       DCC = 0
      ENDIF

C  H ---> MU MU
      IF(AMH.LE.2*AMMUON) THEN
       HMM = 0
      ELSE
      HMM=HFF(AMH,(AMMUON/AMH)**2)
     .    *(1+ELW(AMH,AMMUON,-1.D0,7.D0))
     .    *HFFSELF(AMH)
      ENDIF
C  H ---> TAU TAU
      IF(AMH.LE.2*AMTAU) THEN
       HLL = 0
      ELSE
      HLL=HFF(AMH,(AMTAU/AMH)**2)
     .    *(1+ELW(AMH,AMTAU,-1.D0,7.D0))
     .    *HFFSELF(AMH)
      ENDIF
C  H --> SS
      IF(AMH.LE.2*AMS) THEN
       HSS = 0
      ELSE
       HS2=3.D0*HFF(AMH,(RMS/AMH)**2)
     .    *QCDH(RMS**2/AMH**2)
     .    *(1+ELW(AMH,RMS,-1.D0/3.D0,7.D0))
     .    *HFFSELF(AMH)
       IF(HS2.LT.0.D0) HS2 = 0
       HS1=3.D0*HFF(AMH,(AMS/AMH)**2)
     .    *TQCDH(AMS**2/AMH**2)
     .    *HFFSELF(AMH)
       RAT = 2*AMS/AMH
       HSS = QQINT(RAT,HS1,HS2)
      ENDIF
C  H --> CC
      IF(AMH.LE.2*AMC) THEN
       HCC = 0
      ELSE
       HC2=3.D0*HFF(AMH,(RMC/AMH)**2)
     .    *QCDH(RMC**2/AMH**2)
     .    *(1+ELW(AMH,RMC,2.D0/3.D0,7.D0))
     .    *HFFSELF(AMH)
     .   + DCC
       IF(HC2.LT.0.D0) HC2 = 0
       HC1=3.D0*HFF(AMH,(AMC/AMH)**2)
     .    *TQCDH(AMC**2/AMH**2)
     .    *HFFSELF(AMH)
       RAT = 2*AMC/AMH
       HCC = QQINT(RAT,HC1,HC2)
      ENDIF
C  H --> BB :
      IF(AMH.LE.2*AMB) THEN
       HBB = 0
      ELSE
       HB2=3.D0*HFF(AMH,(RMB/AMH)**2)
     .    *QCDH(RMB**2/AMH**2)
     .    *(1+ELW(AMH,RMB,-1.D0/3.D0,1.D0))
     .    *HFFSELF(AMH)
     .   + DBB
       IF(HB2.LT.0.D0) HB2 = 0
       HB1=3.D0*HFF(AMH,(AMB/AMH)**2)
     .    *TQCDH(AMB**2/AMH**2)
     .    *HFFSELF(AMH)
       RAT = 2*AMB/AMH
       HBB = QQINT(RAT,HB1,HB2)
      ENDIF
C  H ---> TT
      RATCOUP = 0
      IF(IONSH.EQ.0)THEN
       DLD=3D0
       DLU=5D0
       XM1 = 2D0*AMT-DLD
       XM2 = 2D0*AMT+DLU
       IF (AMH.LE.AMT+AMW+AMB) THEN
       HTT=0.D0
       ELSEIF (AMH.LE.XM1) THEN
        FACTT=6.D0*GF**2*AMH**3*AMT**2/2.D0/128.D0/PI**3
        CALL HTOTTS(AMH,AMT,AMB,AMW,HTTS)
        HTT=FACTT*HTTS
       ELSEIF (AMH.LE.XM2) THEN
        XX(1) = XM1-1D0
        XX(2) = XM1
        XX(3) = XM2
        XX(4) = XM2+1D0
        FACTT=6.D0*GF**2*XX(1)**3*AMT**2/2.D0/128.D0/PI**3
        CALL HTOTTS(XX(1),AMT,AMB,AMW,HTTS)
        YY(1)=FACTT*HTTS
        FACTT=6.D0*GF**2*XX(2)**3*AMT**2/2.D0/128.D0/PI**3
        CALL HTOTTS(XX(2),AMT,AMB,AMW,HTTS)
        YY(2)=FACTT*HTTS
        XMT = RUNM(XX(3),6)
        XY2=3.D0*HFF(XX(3),(XMT/XX(3))**2)
     .    *QCDH(XMT**2/XX(3)**2)
     .    *HFFSELF(XX(3))
        IF(XY2.LT.0.D0) XY2 = 0
        XY1=3.D0*HFF(XX(3),(AMT/XX(3))**2)
     .    *TQCDH(AMT**2/XX(3)**2)
     .    *HFFSELF(XX(3))
        RAT = 2*AMT/XX(3)
        YY(3) = QQINT(RAT,XY1,XY2)
        XMT = RUNM(XX(4),6)
        XY2=3.D0*HFF(XX(4),(XMT/XX(4))**2)
     .    *QCDH(XMT**2/XX(4)**2)
     .    *HFFSELF(XX(4))
        IF(XY2.LT.0.D0) XY2 = 0
        XY1=3.D0*HFF(XX(4),(AMT/XX(4))**2)
     .    *TQCDH(AMT**2/XX(4)**2)
     .    *HFFSELF(XX(4))
        RAT = 2*AMT/XX(4)
        YY(4) = QQINT(RAT,XY1,XY2)
        HTT = FINT(AMH,XX,YY)
       ELSE
        HT2=3.D0*HFF(AMH,(RMT/AMH)**2)
     .    *QCDH(RMT**2/AMH**2)
     .    *HFFSELF(AMH)
        IF(HT2.LT.0.D0) HT2 = 0
        HT1=3.D0*HFF(AMH,(AMT/AMH)**2)
     .    *TQCDH(AMT**2/AMH**2)
     .    *HFFSELF(AMH)
        RAT = 2*AMT/AMH
        HTT = QQINT(RAT,HT1,HT2)
       ENDIF
      ELSE
       IF (AMH.LE.2.D0*AMT) THEN
        HTT=0.D0
       ELSE
        HT2=3.D0*HFF(AMH,(RMT/AMH)**2)
     .    *QCDH(RMT**2/AMH**2)
     .    *HFFSELF(AMH)
        IF(HT2.LT.0.D0) HT2 = 0
        HT1=3.D0*HFF(AMH,(AMT/AMH)**2)
     .    *TQCDH(AMT**2/AMH**2)
     .    *HFFSELF(AMH)
        RAT = 2*AMT/AMH
        HTT = QQINT(RAT,HT1,HT2)
       ENDIF
      ENDIF
C  H ---> GAMMA GAMMA
       EPS=1.D-8
       CTT = 4*AMT**2/AMH**2*DCMPLX(1D0,-EPS)
       CTB = 4*AMB**2/AMH**2*DCMPLX(1D0,-EPS)
       CTC = 4*AMC**2/AMH**2*DCMPLX(1D0,-EPS)
       CTL = 4*AMTAU**2/AMH**2*DCMPLX(1D0,-EPS)
       CTW = 4*AMW**2/AMH**2*DCMPLX(1D0,-EPS)
       CAW = -(2+3*CTW+3*CTW*(2-CTW)*CF(CTW))
       CAT = 4/3D0 * 2*CTT*(1+(1-CTT)*CF(CTT))
       CAB = 1/3D0 * 2*CTB*(1+(1-CTB)*CF(CTB))
       CAC = 4/3D0 * 2*CTC*(1+(1-CTC)*CF(CTC))
       CAL =         2*CTL*(1+(1-CTL)*CF(CTL))
       XFAC = CDABS(CAT+CAB+CAC+CAL+CAW)**2
       HGA=HVV(AMH,0.D0)*(ALPH/PI)**2/16.D0*XFAC
C  H ---> Z GAMMA
      IF(AMH.LE.AMZ)THEN
       HZGA=0
      ELSE
       EPS=1.D-8
       TS = SS/CS
       FT = -3*2D0/3*(1-4*2D0/3*SS)/DSQRT(SS*CS)
       FB = 3*1D0/3*(-1+4*1D0/3*SS)/DSQRT(SS*CS)
       CTT = 4*AMT**2/AMH**2*DCMPLX(1D0,-EPS)
       CTB = 4*AMB**2/AMH**2*DCMPLX(1D0,-EPS)
       CTW = 4*AMW**2/AMH**2*DCMPLX(1D0,-EPS)
       CLT = 4*AMT**2/AMZ**2*DCMPLX(1D0,-EPS)
       CLB = 4*AMB**2/AMZ**2*DCMPLX(1D0,-EPS)
       CLW = 4*AMW**2/AMZ**2*DCMPLX(1D0,-EPS)
       CAT = FT*(CI1(CTT,CLT) - CI2(CTT,CLT))
       CAB = FB*(CI1(CTB,CLB) - CI2(CTB,CLB))
       CAW = -1/DSQRT(TS)*(4*(3-TS)*CI2(CTW,CLW)
     .     + ((1+2/CTW)*TS - (5+2/CTW))*CI1(CTW,CLW))
       XFAC = CDABS(CAT+CAB+CAW)**2
       ACOUP = DSQRT(2D0)*GF*AMZ**2*SS*CS/PI**2
       HZGA = GF/(4.D0*PI*DSQRT(2.D0))*AMH**3*(ALPH/PI)*ACOUP/16.D0
     .        *XFAC*(1-AMZ**2/AMH**2)**3
      ENDIF
C  H ---> W W
      IF(IONWZ.EQ.0)THEN
       DLD=2D0
       DLU=2D0
       XM1 = 2D0*AMW-DLD
       XM2 = 2D0*AMW+DLU
       IF (AMH.LE.XM1) THEN
        CALL HTOVV(AMH,AMW,GAMW,HTWW)
        HWW = 3D0/2D0*GF*AMW**4/DSQRT(2D0)/PI/AMH**3*HTWW
       ELSEIF (AMH.LE.XM2) THEN
        XX(1) = XM1-1D0
        XX(2) = XM1
        XX(3) = XM2
        XX(4) = XM2+1D0
        CALL HTOVV(XX(1),AMW,GAMW,HTWW)
        YY(1)=3D0/2D0*GF*AMW**4/DSQRT(2D0)/PI/XX(1)**3*HTWW
        CALL HTOVV(XX(2),AMW,GAMW,HTWW)
        YY(2)=3D0/2D0*GF*AMW**4/DSQRT(2D0)/PI/XX(2)**3*HTWW
        YY(3)=HVV(XX(3),AMW**2/XX(3)**2)
     .       *HVVSELF(XX(3))
        YY(4)=HVV(XX(4),AMW**2/XX(4)**2)
     .       *HVVSELF(XX(4))
        HWW = FINT(AMH,XX,YY)
       ELSE
        HWW=HVV(AMH,AMW**2/AMH**2)
     .     *HVVSELF(AMH)
       ENDIF
      ELSE
       DLD=2D0
       DLU=2D0
       XM1 = 2D0*AMW-DLD
       XM2 = 2D0*AMW+DLU
      IF (AMH.LE.AMW) THEN
       HWW=0
      ELSE IF (AMH.LE.XM1) THEN
       CWW=3.D0*GF**2*AMW**4/16.D0/PI**3
       HWW=HV(AMW**2/AMH**2)*CWW*AMH
      ELSE IF (AMH.LT.XM2) THEN
       CWW=3.D0*GF**2*AMW**4/16.D0/PI**3
       XX(1) = XM1-1D0
       XX(2) = XM1
       XX(3) = XM2
       XX(4) = XM2+1D0
       YY(1)=HV(AMW**2/XX(1)**2)*CWW*XX(1)
       YY(2)=HV(AMW**2/XX(2)**2)*CWW*XX(2)
       YY(3)=HVV(XX(3),AMW**2/XX(3)**2)
     .      *HVVSELF(XX(3))
       YY(4)=HVV(XX(4),AMW**2/XX(4)**2)
     .      *HVVSELF(XX(4))
       HWW = FINT(AMH,XX,YY)
      ELSE
       HWW=HVV(AMH,AMW**2/AMH**2)
     .     *HVVSELF(AMH)
      ENDIF
      ENDIF
C  H ---> Z Z
      IF(IONWZ.EQ.0)THEN
       DLD=2D0
       DLU=2D0
       XM1 = 2D0*AMZ-DLD
       XM2 = 2D0*AMZ+DLU
       IF (AMH.LE.XM1) THEN
        CALL HTOVV(AMH,AMZ,GAMZ,HTZZ)
        HZZ = 3D0/4D0*GF*AMZ**4/DSQRT(2D0)/PI/AMH**3*HTZZ
       ELSEIF (AMH.LE.XM2) THEN
        XX(1) = XM1-1D0
        XX(2) = XM1
        XX(3) = XM2
        XX(4) = XM2+1D0
        CALL HTOVV(XX(1),AMZ,GAMZ,HTZZ)
        YY(1)=3D0/4D0*GF*AMZ**4/DSQRT(2D0)/PI/XX(1)**3*HTZZ
        CALL HTOVV(XX(2),AMZ,GAMZ,HTZZ)
        YY(2)=3D0/4D0*GF*AMZ**4/DSQRT(2D0)/PI/XX(2)**3*HTZZ
        YY(3)=HVV(XX(3),AMZ**2/XX(3)**2)/2
     .       *HVVSELF(XX(3))
        YY(4)=HVV(XX(4),AMZ**2/XX(4)**2)/2
     .       *HVVSELF(XX(4))
        HZZ = FINT(AMH,XX,YY)
       ELSE
        HZZ=HVV(AMH,AMZ**2/AMH**2)/2.D0
     .      *HVVSELF(AMH)
       ENDIF
      ELSE
      DLD=2D0
      DLU=2D0
      XM1 = 2D0*AMZ-DLD
      XM2 = 2D0*AMZ+DLU
      IF (AMH.LE.AMZ) THEN
      HZZ=0
      ELSE IF (AMH.LE.XM1) THEN
      CZZ=3.D0*GF**2*AMZ**4/192.D0/PI**3*(7-40/3.D0*SS+160/9.D0*SS**2)
      HZZ=HV(AMZ**2/AMH**2)*CZZ*AMH
      ELSE IF (AMH.LT.XM2) THEN
      CZZ=3.D0*GF**2*AMZ**4/192.D0/PI**3*(7-40/3.D0*SS+160/9.D0*SS**2)
      XX(1) = XM1-1D0
      XX(2) = XM1
      XX(3) = XM2
      XX(4) = XM2+1D0
      YY(1)=HV(AMZ**2/XX(1)**2)*CZZ*XX(1)
      YY(2)=HV(AMZ**2/XX(2)**2)*CZZ*XX(2)
      YY(3)=HVV(XX(3),AMZ**2/XX(3)**2)/2
     .      *HVVSELF(XX(3))
      YY(4)=HVV(XX(4),AMZ**2/XX(4)**2)/2
     .      *HVVSELF(XX(4))
      HZZ = FINT(AMH,XX,YY)
      ELSE
      HZZ=HVV(AMH,AMZ**2/AMH**2)/2.D0
     .   *HVVSELF(AMH)
      ENDIF
      ENDIF
C
C    ==========  TOTAL WIDTH AND BRANCHING RATIOS 
C
      WTOT=HLL+HMM+HSS+HCC+HBB+HTT+HGG+HGA+HZGA+HWW+HZZ
      SMBRT=HTT/WTOT
      SMBRB=HBB/WTOT
      SMBRL=HLL/WTOT
      SMBRM=HMM/WTOT
      SMBRC=HCC/WTOT
      SMBRS=HSS/WTOT
      SMBRG=HGG/WTOT
      SMBRGA=HGA/WTOT
      SMBRZGA=HZGA/WTOT
      SMBRW=HWW/WTOT
      SMBRZ=HZZ/WTOT
      SMWDTH=WTOT

      AMH=AMXX

      endif

      IF(IHIGGS.GT.0)THEN

C +++++++++++++++++++++++  SUSY HIGGSSES +++++++++++++++++++++++
C
      CALL GAUGINO(AMU,AM2,B,A,AMCHAR,AMNEUT,XMNEUT,AC1,AC2,AC3,
     .           AN1,AN2,AN3,ACNL,ACNR)         
C
      CALL SFERMION(AMSQ,AMUR,AMDR,AMEL,AMER,AL,AU,AD,AMU,
     .               GMST,GMSB,GMSL,GMSU,GMSD,GMSE,GMSN, 
     .               GLEE,GLTT,GLBB,GHEE,GHTT,GHBB,
     .               GAEE,GATT,GABB,GCEN,GCTB)
C
      ENDIF

      IF(IHIGGS.EQ.1.OR.IHIGGS.EQ.5)THEN
C        =========================================================
C                           LIGHT CP EVEN HIGGS DECAYS
C        =========================================================
C     =============  RUNNING MASSES 
      RMS = RUNM(AML,3)
      RMC = RUNM(AML,4)
      RMB = RUNM(AML,5)
      RMT = RUNM(AML,6)
      RATCOUP = GLT/GLB
      HIGTOP = AML**2/AMT**2

      ASH=ALPHAS(AML,2)
      AMC0=1.D8
      AMB0=2.D8
C     AMT0=3.D8
      AS3=ALPHAS(AML,2)
      AMC0=AMC
      AS4=ALPHAS(AML,2)
      AMB0=AMB
C     AMT0=AMT

C     =============== PARTIAL WIDTHS 
C  H ---> G G
       EPS=1.D-8
       NFEXT = 3
       ASG = AS3
       CTT = 4*AMT**2/AML**2*DCMPLX(1D0,-EPS)
       CTB = 4*AMB**2/AML**2*DCMPLX(1D0,-EPS)
       CAT = 2*CTT*(1+(1-CTT)*CF(CTT))*GLT
       CAB = 2*CTB*(1+(1-CTB)*CF(CTB))*GLB
       CTC = 4*AMC**2/AML**2*DCMPLX(1D0,-EPS)
       CAC = 2*CTC*(1+(1-CTC)*CF(CTC))*GLT
C
       IF(IOFSUSY.EQ.0) THEN 
       CSB1= 4*AMSB(1)**2/AML**2*DCMPLX(1D0,-EPS)
       CSB2= 4*AMSB(2)**2/AML**2*DCMPLX(1D0,-EPS)
       CST1= 4*AMST(1)**2/AML**2*DCMPLX(1D0,-EPS)
       CST2= 4*AMST(2)**2/AML**2*DCMPLX(1D0,-EPS)
       CXB1=-AMZ**2/AMSB(1)**2*CSB1*(1-CSB1*CF(CSB1))*GLBB(1,1)
       CXB2=-AMZ**2/AMSB(2)**2*CSB2*(1-CSB2*CF(CSB2))*GLBB(2,2)
       CXT1=-AMZ**2/AMST(1)**2*CST1*(1-CST1*CF(CST1))*GLTT(1,1)
       CXT2=-AMZ**2/AMST(2)**2*CST2*(1-CST2*CF(CST2))*GLTT(2,2)

       CSUL = 4*AMSU(1)**2/AML**2*DCMPLX(1D0,-EPS)
       CSUR = 4*AMSU(2)**2/AML**2*DCMPLX(1D0,-EPS)
       CSDL = 4*AMSD(1)**2/AML**2*DCMPLX(1D0,-EPS)
       CSDR = 4*AMSD(2)**2/AML**2*DCMPLX(1D0,-EPS)
       CXUL=2*(1.D0/2.D0-2.D0/3.D0*SS)*AMZ**2/AMSU(1)**2*DSIN(A+B)
     .      *CSUL*(1-CSUL*CF(CSUL))
       CXUR=2*(2.D0/3.D0*SS)*AMZ**2/AMSU(2)**2*DSIN(A+B)
     .      *CSUR*(1-CSUR*CF(CSUR))
       CXDL=2*(-1.D0/2.D0+1.D0/3.D0*SS)*AMZ**2/AMSD(1)**2*DSIN(A+B)
     .      *CSDL*(1-CSDL*CF(CSDL))
       CXDR=2*(-1.D0/3.D0*SS)*AMZ**2/AMSD(2)**2*DSIN(A+B)
     .      *CSDR*(1-CSDR*CF(CSUR))

       ELSE
       CXB1=0.D0 
       CXB2=0.D0 
       CXT1=0.D0 
       CXT2=0.D0 

       CXUL=0.D0 
       CXUR=0.D0 
       CXDL=0.D0 
       CXDR=0.D0 
       ENDIF

       FQCD=HGGQCD(ASG,NFEXT)
       SQCD=SGGQCD(ASG)
       XFAC = CDABS(CAT+CAB+CAC+CXB1+CXB2+CXT1+CXT2
     .             +CXUL+CXUR+CXDL+CXDR)**2*FQCD
     .      + DREAL(DCONJG(CAT+CAB+CAC)*(CXB1+CXB2+CXT1+CXT2
     .             +CXUL+CXUR+CXDL+CXDR))*SQCD
       HGG=HVV(AML,0.D0)*(ASG/PI)**2*XFAC/8

C  H ---> G G* ---> G CC   TO BE ADDED TO H ---> CC
       NFEXT = 4
       ASG = AS4
       FQCD=HGGQCD(ASG,NFEXT)
       SQCD=SGGQCD(ASG)
       XFAC = CDABS(CAT+CAB+CAC+CXB1+CXB2+CXT1+CXT2
     .             +CXUL+CXUR+CXDL+CXDR)**2*FQCD
     .      + DREAL(DCONJG(CAT+CAB+CAC)*(CXB1+CXB2+CXT1+CXT2
     .             +CXUL+CXUR+CXDL+CXDR))*SQCD
       DCC=HVV(AML,0.D0)*(ASG/PI)**2*XFAC/8 - HGG

C  H ---> G G* ---> G BB   TO BE ADDED TO H ---> BB
       NFEXT = 5
       ASG = ASH
       FQCD=HGGQCD(ASG,NFEXT)
       SQCD=SGGQCD(ASG)
       XFAC = CDABS(CAT+CAB+CAC+CXB1+CXB2+CXT1+CXT2
     .             +CXUL+CXUR+CXDL+CXDR)**2*FQCD
     .      + DREAL(DCONJG(CAT+CAB+CAC)*(CXB1+CXB2+CXT1+CXT2
     .             +CXUL+CXUR+CXDL+CXDR))*SQCD
       DBB=HVV(AML,0.D0)*(ASG/PI)**2*XFAC/8 - HGG - DCC

      IF(NFGG.EQ.5)THEN
       HGG = HGG + DBB + DCC
       DBB = 0
       DCC = 0
      ELSEIF(NFGG.EQ.4)THEN
       HGG = HGG + DCC
       DCC = 0
      ENDIF

C  H ---> MU MU
      IF(AML.LE.2*AMMUON) THEN
       HMM = 0
      ELSE
      HMM=HFF(AML,(AMMUON/AML)**2)*GLB**2
      ENDIF
C  H ---> TAU TAU
      IF(AML.LE.2*AMTAU) THEN
       HLL = 0
      ELSE
      HLL=HFF(AML,(AMTAU/AML)**2)*GLB**2
      ENDIF
C  H --> SS
      IF(AML.LE.2*AMS) THEN
       HSS = 0
      ELSE
       HS1=3.D0*HFF(AML,(AMS/AML)**2)
     .    *GLB**2
     .    *TQCDH(AMS**2/AML**2)
       HS2=3.D0*HFF(AML,(RMS/AML)**2)*GLB**2
     .    *QCDH(RMS**2/AML**2)
       IF(HS2.LT.0.D0) HS2 = 0
       RAT = 2*AMS/AML
       HSS = QQINT(RAT,HS1,HS2)
      ENDIF
C  H --> CC
      RATCOUP = 1
      IF(AML.LE.2*AMC) THEN
       HCC = 0
      ELSE
       HC1=3.D0*HFF(AML,(AMC/AML)**2)
     .    *GLT**2
     .    *TQCDH(AMC**2/AML**2)
       HC2=3.D0*HFF(AML,(RMC/AML)**2)*GLT**2
     .    *QCDH(RMC**2/AML**2)
     .   + DCC
       IF(HC2.LT.0.D0) HC2 = 0
       RAT = 2*AMC/AML
       HCC = QQINT(RAT,HC1,HC2)
      ENDIF
C  H --> BB :
      RATCOUP = GLT/GLB
      IF(AML.LE.2*AMB) THEN
       HBB = 0
      ELSE
       HB1=3.D0*HFF(AML,(AMB/AML)**2)
     .    *GLB**2
     .    *TQCDH(AMB**2/AML**2)
       HB2=3.D0*HFF(AML,(RMB/AML)**2)*GLB**2
     .    *QCDH(RMB**2/AML**2)
     .   + DBB
       IF(HB2.LT.0.D0) HB2 = 0
       RAT = 2*AMB/AML
       HBB = QQINT(RAT,HB1,HB2)
      ENDIF
C  H ---> TT
      RATCOUP = 0
      IF (AML.LE.2*AMT) THEN
       HTT=0.D0
      ELSE
       HT1=3.D0*HFF(AML,(AMT/AML)**2)*GLT**2
     .    *TQCDH(AMT**2/AML**2)
       HT2=3.D0*HFF(AML,(RMT/AML)**2)*GLT**2
     .    *QCDH(RMT**2/AML**2)
       IF(HT2.LT.0.D0) HT2 = 0
       RAT = 2*AMT/AML
       HTT = QQINT(RAT,HT1,HT2)
      ENDIF
C  H ---> GAMMA GAMMA
       EPS=1.D-8
       CTT = 4*AMT**2/AML**2*DCMPLX(1D0,-EPS)
       CTB = 4*AMB**2/AML**2*DCMPLX(1D0,-EPS)
       CTC = 4*AMC**2/AML**2*DCMPLX(1D0,-EPS)
       CTL = 4*AMTAU**2/AML**2*DCMPLX(1D0,-EPS)
       CTW = 4*AMW**2/AML**2*DCMPLX(1D0,-EPS)
       CTH = 4*AMCH**2/AML**2*DCMPLX(1D0,-EPS)
       CAT = 4/3D0 * 2*CTT*(1+(1-CTT)*CF(CTT))*GLT
       CAB = 1/3D0 * 2*CTB*(1+(1-CTB)*CF(CTB))*GLB
       CAC = 4/3D0 * 2*CTC*(1+(1-CTC)*CF(CTC))*GLT
       CAL = 1.D0  * 2*CTL*(1+(1-CTL)*CF(CTL))*GLB
       CAW = -(2+3*CTW+3*CTW*(2-CTW)*CF(CTW))*GLVV
       CAH = -AMZ**2/2/AMCH**2*CTH*(1-CTH*CF(CTH))*GLPM
       IF(IOFSUSY.EQ.0) THEN 
       CX1 = 4*AMCHAR(1)**2/AML**2*DCMPLX(1D0,-EPS)
       CX2 = 4*AMCHAR(2)**2/AML**2*DCMPLX(1D0,-EPS)
       CSB1= 4*AMSB(1)**2/AML**2*DCMPLX(1D0,-EPS)
       CSB2= 4*AMSB(2)**2/AML**2*DCMPLX(1D0,-EPS)
       CST1= 4*AMST(1)**2/AML**2*DCMPLX(1D0,-EPS)
       CST2= 4*AMST(2)**2/AML**2*DCMPLX(1D0,-EPS)
       CSL1= 4*AMSL(1)**2/AML**2*DCMPLX(1D0,-EPS)
       CSL2= 4*AMSL(2)**2/AML**2*DCMPLX(1D0,-EPS)
       CAX1= AMW/AMCHAR(1) * 2*CX1*(1+(1-CX1)*CF(CX1))*2*AC2(1,1) 
       CAX2= AMW/AMCHAR(2) * 2*CX2*(1+(1-CX2)*CF(CX2))*2*AC2(2,2) 

       CSEL = 4*AMSE(1)**2/AML**2*DCMPLX(1D0,-EPS)
       CSER = 4*AMSE(2)**2/AML**2*DCMPLX(1D0,-EPS)
       CSUL = 4*AMSU(1)**2/AML**2*DCMPLX(1D0,-EPS)
       CSUR = 4*AMSU(2)**2/AML**2*DCMPLX(1D0,-EPS)
       CSDL = 4*AMSD(1)**2/AML**2*DCMPLX(1D0,-EPS)
       CSDR = 4*AMSD(2)**2/AML**2*DCMPLX(1D0,-EPS)
       CXEL=2*(-1/2D0+SS)*AMZ**2/AMSE(1)**2*DSIN(A+B)
     .      *CSEL*(1-CSEL*CF(CSEL))
       CXER=-2*(SS)*AMZ**2/AMSE(2)**2*DSIN(A+B)
     .      *CSER*(1-CSER*CF(CSER))
       CXUL=2*4.D0/3.D0*(1.D0/2.D0-2.D0/3.D0*SS)
     .      *AMZ**2/AMSU(1)**2*DSIN(A+B)*CSUL*(1-CSUL*CF(CSUL))
       CXUR=2*4.D0/3.D0*(2.D0/3.D0*SS)
     .      *AMZ**2/AMSU(2)**2*DSIN(A+B)*CSUR*(1-CSUR*CF(CSUR))
       CXDL=2/3.D0*(-1.D0/2.D0+1.D0/3.D0*SS)
     .      *AMZ**2/AMSD(1)**2*DSIN(A+B)*CSDL*(1-CSDL*CF(CSDL))
       CXDR=2/3.D0*(-1.D0/3.D0*SS)
     .      *AMZ**2/AMSD(2)**2*DSIN(A+B)*CSDR*(1-CSDR*CF(CSUR))

       CXB1=-1/3D0*AMZ**2/AMSB(1)**2*CSB1*(1-CSB1*CF(CSB1))*GLBB(1,1)
       CXB2=-1/3D0*AMZ**2/AMSB(2)**2*CSB2*(1-CSB2*CF(CSB2))*GLBB(2,2)
       CXT1=-4/3D0*AMZ**2/AMST(1)**2*CST1*(1-CST1*CF(CST1))*GLTT(1,1)
       CXT2=-4/3D0*AMZ**2/AMST(2)**2*CST2*(1-CST2*CF(CST2))*GLTT(2,2)
       CSL1= 4*AMSL(1)**2/AML**2*DCMPLX(1D0,-EPS)
       CSL2= 4*AMSL(2)**2/AML**2*DCMPLX(1D0,-EPS)
       CXL1=      -AMZ**2/AMSL(1)**2*CSL1*(1-CSL1*CF(CSL1))*GLEE(1,1)
       CXL2=      -AMZ**2/AMSL(2)**2*CSL2*(1-CSL2*CF(CSL2))*GLEE(2,2)
       XFAC = CDABS(CAT+CAB+CAC+CAL+CAW+CAH+CAX1+CAX2
     .     +  CXEL+CXER+CXUL+CXUR+CXDL+CXDR
     .     +  CXB1+CXB2+CXT1+CXT2+CXL1+CXL2)**2
       ELSE 
       XFAC = CDABS(CAT+CAB+CAC+CAL+CAW+CAH)**2
       ENDIF
       HGA=HVV(AML,0.D0)*(ALPH/PI)**2/16.D0*XFAC
C  H ---> Z GAMMA
      IF(AML.LE.AMZ)THEN
       HZGA=0
      ELSE
       TS = SS/CS
       FT = -3*2D0/3*(1-4*2D0/3*SS)/DSQRT(SS*CS)*GLT
       FB = 3*1D0/3*(-1+4*1D0/3*SS)/DSQRT(SS*CS)*GLB
       EPS=1.D-8
       CTT = 4*AMT**2/AML**2*DCMPLX(1D0,-EPS)
       CTB = 4*AMB**2/AML**2*DCMPLX(1D0,-EPS)
       CTW = 4*AMW**2/AML**2*DCMPLX(1D0,-EPS)
       CTH = 4*AMCH**2/AML**2*DCMPLX(1D0,-EPS)
       CLT = 4*AMT**2/AMZ**2*DCMPLX(1D0,-EPS)
       CLB = 4*AMB**2/AMZ**2*DCMPLX(1D0,-EPS)
       CLW = 4*AMW**2/AMZ**2*DCMPLX(1D0,-EPS)
       CLH = 4*AMCH**2/AMZ**2*DCMPLX(1D0,-EPS)
       CAT = FT*(CI1(CTT,CLT) - CI2(CTT,CLT))
       CAB = FB*(CI1(CTB,CLB) - CI2(CTB,CLB))
       CAW = -1/DSQRT(TS)*(4*(3-TS)*CI2(CTW,CLW)
     .     + ((1+2/CTW)*TS - (5+2/CTW))*CI1(CTW,CLW))*GLVV
       CAH = (1-2*SS)/DSQRT(SS*CS)*AMZ**2/2/AMCH**2*CI1(CTH,CLH)*GLPM
       XFAC = CDABS(CAT+CAB+CAW+CAH)**2
       ACOUP = DSQRT(2D0)*GF*AMZ**2*SS*CS/PI**2
       HZGA = GF/(4.D0*PI*DSQRT(2.D0))*AML**3*(ALPH/PI)*ACOUP/16.D0
     .        *XFAC*(1-AMZ**2/AML**2)**3
      ENDIF
C  H ---> W W
      IF(IONWZ.EQ.0)THEN
       DLD=2D0
       DLU=2D0
       XM1 = 2D0*AMW-DLD
       XM2 = 2D0*AMW+DLU
       IF (AML.LE.XM1) THEN
        CALL HTOVV(AML,AMW,GAMW,HTWW)
        HWW = 3D0/2D0*GF*AMW**4/DSQRT(2D0)/PI/AML**3*HTWW*GLVV**2
       ELSEIF (AML.LE.XM2) THEN
        XX(1) = XM1-1D0
        XX(2) = XM1
        XX(3) = XM2
        XX(4) = XM2+1D0
        CALL HTOVV(XX(1),AMW,GAMW,HTWW)
        YY(1)=3D0/2D0*GF*AMW**4/DSQRT(2D0)/PI/XX(1)**3*HTWW
        CALL HTOVV(XX(2),AMW,GAMW,HTWW)
        YY(2)=3D0/2D0*GF*AMW**4/DSQRT(2D0)/PI/XX(2)**3*HTWW
        YY(3)=HVV(XX(3),AMW**2/XX(3)**2)
        YY(4)=HVV(XX(4),AMW**2/XX(4)**2)
        HWW = FINT(AML,XX,YY)*GLVV**2
       ELSE
        HWW=HVV(AML,AMW**2/AML**2)*GLVV**2
       ENDIF
      ELSE
      DLD=2D0
      DLU=2D0
      XM1 = 2D0*AMW-DLD
      XM2 = 2D0*AMW+DLU
      IF (AML.LE.AMW) THEN
       HWW=0
      ELSE IF (AML.LE.XM1) THEN
       CWW=3.D0*GF**2*AMW**4/16.D0/PI**3
       HWW=HV(AMW**2/AML**2)*CWW*AML*GLVV**2
      ELSE IF (AML.LT.XM2) THEN
       CWW=3.D0*GF**2*AMW**4/16.D0/PI**3
       XX(1) = XM1-1D0
       XX(2) = XM1
       XX(3) = XM2
       XX(4) = XM2+1D0
       YY(1)=HV(AMW**2/XX(1)**2)*CWW*XX(1)
       YY(2)=HV(AMW**2/XX(2)**2)*CWW*XX(2)
       YY(3)=HVV(XX(3),AMW**2/XX(3)**2)
       YY(4)=HVV(XX(4),AMW**2/XX(4)**2)
       HWW = FINT(AML,XX,YY)*GLVV**2
      ELSE
       HWW=HVV(AML,AMW**2/AML**2)*GLVV**2
      ENDIF
      ENDIF
C  H ---> Z Z
      IF(IONWZ.EQ.0)THEN
       DLD=2D0
       DLU=2D0
       XM1 = 2D0*AMZ-DLD
       XM2 = 2D0*AMZ+DLU
       IF (AML.LE.XM1) THEN
        CALL HTOVV(AML,AMZ,GAMZ,HTZZ)
        HZZ = 3D0/4D0*GF*AMZ**4/DSQRT(2D0)/PI/AML**3*HTZZ*GLVV**2
       ELSEIF (AML.LE.XM2) THEN
        XX(1) = XM1-1D0
        XX(2) = XM1
        XX(3) = XM2
        XX(4) = XM2+1D0
        CALL HTOVV(XX(1),AMZ,GAMZ,HTZZ)
        YY(1)=3D0/4D0*GF*AMZ**4/DSQRT(2D0)/PI/XX(1)**3*HTZZ
        CALL HTOVV(XX(2),AMZ,GAMZ,HTZZ)
        YY(2)=3D0/4D0*GF*AMZ**4/DSQRT(2D0)/PI/XX(2)**3*HTZZ
        YY(3)=HVV(XX(3),AMZ**2/XX(3)**2)/2
        YY(4)=HVV(XX(4),AMZ**2/XX(4)**2)/2
        HZZ = FINT(AML,XX,YY)*GLVV**2
       ELSE
        HZZ=HVV(AML,AMZ**2/AML**2)/2.D0*GLVV**2
       ENDIF
      ELSE
      DLD=2D0
      DLU=2D0
      XM1 = 2D0*AMZ-DLD
      XM2 = 2D0*AMZ+DLU
      IF (AML.LE.AMZ) THEN
       HZZ=0
      ELSE IF (AML.LE.XM1) THEN
       CZZ=3.D0*GF**2*AMZ**4/192.D0/PI**3*(7-40/3.D0*SS+160/9.D0*SS**2)
       HZZ=HV(AMZ**2/AML**2)*CZZ*AML*GLVV**2
      ELSE IF (AML.LT.XM2) THEN
       CZZ=3.D0*GF**2*AMZ**4/192.D0/PI**3*(7-40/3.D0*SS+160/9.D0*SS**2)
       XX(1) = XM1-1D0
       XX(2) = XM1
       XX(3) = XM2
       XX(4) = XM2+1D0
       YY(1)=HV(AMZ**2/XX(1)**2)*CZZ*XX(1)
       YY(2)=HV(AMZ**2/XX(2)**2)*CZZ*XX(2)
       YY(3)=HVV(XX(3),AMZ**2/XX(3)**2)/2D0
       YY(4)=HVV(XX(4),AMZ**2/XX(4)**2)/2D0
       HZZ = FINT(AML,XX,YY)*GLVV**2
      ELSE
       HZZ=HVV(AML,AMZ**2/AML**2)/2.D0*GLVV**2
      ENDIF
      ENDIF
C  H ---> A A
      IF (AML.LE.2.D0*AMA) THEN
      HAA=0
      ELSE
      HAA=GF/16.D0/DSQRT(2D0)/PI*AMZ**4/AML*BETA(AMA**2/AML**2)*GLAA**2
      ENDIF
C  H ---> A Z
      IF (AML.LE.AMZ+AMA) THEN
      HAZ=0
      ELSE
      CAZ=LAMB(AMA**2/AML**2,AMZ**2/AML**2)
     .   *LAMB(AML**2/AMZ**2,AMA**2/AMZ**2)**2
      HAZ=GF/8.D0/DSQRT(2D0)/PI*AMZ**4/AML*CAZ*GZAL**2
      ENDIF
C  H ---> H+ W+
      IF (AML.LE.AMW+AMCH) THEN
      HHW=0
      ELSE
      CHW=LAMB(AMCH**2/AML**2,AMW**2/AML**2)
     .   *LAMB(AML**2/AMW**2,AMCH**2/AMW**2)**2
      HHW=GF/8.D0/DSQRT(2D0)/PI*AMZ**2*AMW**2/AML*CHW*GHVV**2
      ENDIF

C  ============================ SUSY DECAYS 
      IF(IOFSUSY.EQ.0) THEN
C
C  HL ----> CHARGINOS
C
      DO 711 I=1,2
      DO 711 J=1,2
      IF (AML.GT.AMCHAR(I)+AMCHAR(J)) THEN
      WHLCH(I,J)=GF*AMW**2/(2*PI*DSQRT(2.D0))/AML 
     .     *LAMB(AMCHAR(I)**2/AML**2,AMCHAR(J)**2/AML**2)
     .     *( (AC2(I,J)**2+AC2(J,I)**2)*(AML**2-AMCHAR(I)
     .         **2-AMCHAR(J)**2)-4.D0*AC2(I,J)*AC2(J,I)* 
     .         AMCHAR(I)*AMCHAR(J) ) 
      ELSE
      WHLCH(I,J)=0.D0
      ENDIF
      WHLCHT=WHLCH(1,1)+WHLCH(1,2)+WHLCH(2,1)+WHLCH(2,2)
 711  CONTINUE
C
C  HL ----> NEUTRALINOS 
C
      DO 712 I=1,4
      DO 712 J=1,4
      IF (AML.GT.AMNEUT(I)+AMNEUT(J)) THEN
      WHLNE(I,J)=GF*AMW**2/(2*PI*DSQRT(2.D0))/AML 
     .         *AN2(I,J)**2*(AML**2-(XMNEUT(I)+XMNEUT(J))**2)
     .         *LAMB(AMNEUT(I)**2/AML**2,AMNEUT(J)**2/AML**2)
      ELSE 
      WHLNE(I,J)=0.D0
      ENDIF
 712  CONTINUE
      WHLNET= WHLNE(1,1)+WHLNE(1,2)+WHLNE(1,3)+WHLNE(1,4)
     .       +WHLNE(2,1)+WHLNE(2,2)+WHLNE(2,3)+WHLNE(2,4)
     .       +WHLNE(3,1)+WHLNE(3,2)+WHLNE(3,3)+WHLNE(3,4)
     .       +WHLNE(4,1)+WHLNE(4,2)+WHLNE(4,3)+WHLNE(4,4)
CCC
C  HL ----> SLEPTONS 
C
      IF (AML.GT.2.D0*AMSE(1)) THEN
      WHLSLEL=2*GF/2.D0/DSQRT(2D0)/PI*AMZ**4/AML*DSIN(B+A)**2
     .      *BETA(AMSE(1)**2/AML**2)*(-0.5D0+SS)**2
      ELSE
      WHLSLEL=0.D0
      ENDIF

      IF (AML.GT.2.D0*AMSE(2)) THEN
      WHLSLER=2*GF/2.D0/DSQRT(2D0)/PI*AMZ**4/AML*DSIN(B+A)**2
     .      *BETA(AMSE(2)**2/AML**2)*SS**2
      ELSE
      WHLSLER=0.D0
      ENDIF

      IF (AML.GT.2.D0*AMSN(1)) THEN
      WHLSLNL=3*GF/2.D0/DSQRT(2D0)/PI*AMZ**4/AML*DSIN(B+A)**2
     .      *BETA(AMSN(1)**2/AML**2)*0.5D0**2
      ELSE
      WHLSLNL=0.D0
      ENDIF

      DO 718 I=1,2
      DO 718 J=1,2
      IF(AML.GT.AMSL(I)+AMSL(J)) THEN
      WHLSTAU(I,J)=GF*AMZ**4/2.D0/DSQRT(2.D0)/PI*GLEE(I,J)**2*
     .      LAMB(AMSL(I)**2/AML**2,AMSL(J)**2/AML**2)/AML
      ELSE
      WHLSTAU(I,J)=0.D0
      ENDIF
 718  CONTINUE

      WHLSLT=WHLSTAU(1,1)+WHLSTAU(2,1)+WHLSTAU(1,2)+WHLSTAU(2,2) 
     .       +WHLSLEL+WHLSLER+WHLSLNL
C
C  HL ----> SQUARKS 
C
      IF (AML.GT.2.D0*AMSU(1)) THEN
      WHLSQUL=6*GF/2.D0/DSQRT(2D0)/PI*AMZ**4/AML*DSIN(B+A)**2
     .      *BETA(AMSU(1)**2/AML**2)*(0.5D0-2.D0/3.D0*SS)**2
      ELSE
      WHLSQUL=0.D0
      ENDIF

      IF (AML.GT.2.D0*AMSU(2)) THEN
      WHLSQUR=6*GF/2.D0/DSQRT(2D0)/PI*AMZ**4/AML*DSIN(B+A)**2
     .      *BETA(AMSU(2)**2/AML**2)*(-2.D0/3.D0*SS)**2
      ELSE
      WHLSQUR=0.D0
      ENDIF

      IF (AML.GT.2.D0*AMSD(1)) THEN
      WHLSQDL=6*GF/2.D0/DSQRT(2D0)/PI*AMZ**4/AML*DSIN(B+A)**2
     .      *BETA(AMSD(1)**2/AML**2)*(-0.5D0+1.D0/3.D0*SS)**2
      ELSE
      WHLSQDL=0.D0
      ENDIF

      IF (AML.GT.2.D0*AMSD(2)) THEN
      WHLSQDR=6*GF/2.D0/DSQRT(2D0)/PI*AMZ**4/AML*DSIN(B+A)**2
     .      *BETA(AMSD(2)**2/AML**2)*(+1.D0/3.D0*SS)**2
      ELSE
      WHLSQDR=0.D0
      ENDIF

      WHLSQ=WHLSQUL+WHLSQUR+WHLSQDL+WHLSQDR
      
C
C  HL ----> STOPS 
      DO 713 I=1,2
      DO 713 J=1,2
      IF(AML.GT.AMST(I)+AMST(J)) THEN
      WHLST(I,J)=3*GF*AMZ**4/2.D0/DSQRT(2.D0)/PI*GLTT(I,J)**2*
     .      LAMB(AMST(I)**2/AML**2,AMST(J)**2/AML**2)/AML
      ELSE
      WHLST(I,J)=0.D0
      ENDIF
 713  CONTINUE
C
C  HL ----> SBOTTOMS 
      DO 714 I=1,2
      DO 714 J=1,2
      IF(AML.GT.AMSB(I)+AMSB(J)) THEN
      WHLSB(I,J)=3*GF*AMZ**4/2.D0/DSQRT(2.D0)/PI*GLBB(I,J)**2*
     .      LAMB(AMSB(I)**2/AML**2,AMSB(J)**2/AML**2)/AML
      ELSE
      WHLSB(I,J)=0.D0
      ENDIF
 714  CONTINUE
C
      WHLSTT=WHLST(1,1)+WHLST(1,2)+WHLST(2,1)+WHLST(2,2) 
      WHLSBB=WHLSB(1,1)+WHLSB(1,2)+WHLSB(2,1)+WHLSB(2,2) 
      WHLSQT=WHLSTT+WHLSBB+WHLSQ

      ELSE 
      WHLCHT=0.D0
      WHLNET=0.D0
      WHLSLT=0.D0
      WHLSQT=0.D0
C--Change thanks to Elzbieta Richter-Was
      DO I=1,2
       DO J=1,2
        WHLCH(I,J)=0.D0
        WHLST(I,J)=0.D0
        WHLSB(I,J)=0.D0
        WHLSTAU(I,J)=0.D0
       ENDDO
      ENDDO
      DO I=1,4
       DO J=1,4
        WHLNE(I,J)=0.D0
       ENDDO
      ENDDO
      ENDIF

C    ==========  TOTAL WIDTH AND BRANCHING RATIOS 
      WTOT=HLL+HMM+HSS+HCC+HBB+HTT+HGG+HGA+HZGA+HWW+HZZ+HAA+HAZ+HHW
     .    +WHLCHT+WHLNET+WHLSLT+WHLSQT
      HLBRT=HTT/WTOT
      HLBRB=HBB/WTOT
      HLBRL=HLL/WTOT
      HLBRM=HMM/WTOT
      HLBRS=HSS/WTOT
      HLBRC=HCC/WTOT
      HLBRG=HGG/WTOT
      HLBRGA=HGA/WTOT
      HLBRZGA=HZGA/WTOT
      HLBRW=HWW/WTOT
      HLBRZ=HZZ/WTOT
      HLBRA=HAA/WTOT
      HLBRAZ=HAZ/WTOT
      HLBRHW=HHW/WTOT
      DO 811 I=1,2
      DO 811 J=1,2
      HLBRSC(I,J)=WHLCH(I,J)/WTOT
811   CONTINUE
      DO 812 I=1,4
      DO 812 J=1,4
      HLBRSN(I,J)=WHLNE(I,J)/WTOT
812   CONTINUE
      HLBRCHT=WHLCHT/WTOT 
      HLBRNET=WHLNET/WTOT 
      HLBRSL=WHLSLT/WTOT 
      HLBRSQ=WHLSQ/WTOT 
      HLBRSQT=WHLSQT/WTOT 
      HLWDTH=WTOT

      BHLSLNL = WHLSLNL/WTOT
      BHLSLEL = WHLSLEL/WTOT
      BHLSLER = WHLSLER/WTOT
      BHLSQUL = WHLSQUL/WTOT
      BHLSQUR = WHLSQUR/WTOT
      BHLSQDL = WHLSQDL/WTOT
      BHLSQDR = WHLSQDR/WTOT
      DO I = 1,2
       DO J = 1,2
        BHLST(I,J) = WHLST(I,J)/WTOT
        BHLSB(I,J) = WHLSB(I,J)/WTOT
        BHLSTAU(I,J) = WHLSTAU( I,J)/WTOT
       ENDDO
      ENDDO

      ENDIF

      IF(IHIGGS.GT.1)THEN
      

C        =========================================================
C                       CHARGED HIGGS DECAYS
C        =========================================================
      TB=TGBET
C     =============  RUNNING MASSES 
      RMS = RUNM(AMCH,3)
      RMC = RUNM(AMCH,4)
      RMB = RUNM(AMCH,5)
      RMT = RUNM(AMCH,6)
      ASH=ALPHAS(AMCH,2)
C     =============== PARTIAL WIDTHS 
C  H+ ---> MU NMU
      IF(AMCH.LE.AMMUON) THEN
       HMN = 0
      ELSE
      HMN=CFF(AMCH,TB,(AMMUON/AMCH)**2,0.D0)
      ENDIF
C  H+ ---> TAU NTAU
      IF(AMCH.LE.AMTAU) THEN
       HLN = 0
      ELSE
      HLN=CFF(AMCH,TB,(AMTAU/AMCH)**2,0.D0)
      ENDIF
C  H+ --> SU
      EPS = 1.D-12
      IF(AMCH.LE.AMS+EPS) THEN
       HSU = 0
      ELSE
       HSU1=3.D0*VUS**2*CQCDM(AMCH,TB,(AMS/AMCH)**2,EPS)
       HSU2=3.D0*VUS**2*CQCD(AMCH,TB,(RMS/AMCH)**2,EPS)
       IF(HSU2.LT.0.D0) HSU2 = 0
       RAT = AMS/AMCH
       HSU = QQINT(RAT,HSU1,HSU2)
      ENDIF
C  H+ --> CS
      IF(AMCH.LE.AMS+AMC) THEN
       HSC = 0
      ELSE
       HSC1=3.D0*CQCDM(AMCH,TB,(AMS/AMCH)**2,(AMC/AMCH)**2)
       HSC2=3.D0*CQCD(AMCH,TB,(RMS/AMCH)**2,(RMC/AMCH)**2)
       IF(HSC2.LT.0.D0) HSC2 = 0
       RAT = (AMS+AMC)/AMCH
       HSC = QQINT(RAT,HSC1,HSC2)
      ENDIF
C  H+ --> CB
      IF(AMCH.LE.AMB+AMC) THEN
       HBC = 0
      ELSE
       HBC1=3.D0*VCB**2*CQCDM(AMCH,TB,(AMB/AMCH)**2,(AMC/AMCH)**2)
       HBC2=3.D0*VCB**2*CQCD(AMCH,TB,(RMB/AMCH)**2,(RMC/AMCH)**2)
       IF(HBC2.LT.0.D0) HBC2 = 0
       RAT = (AMB+AMC)/AMCH
       HBC = QQINT(RAT,HBC1,HBC2)
      ENDIF
C  H+ --> BU
      EPS = 1.D-12
      IF(AMCH.LE.AMB+EPS) THEN
       HBU = 0
      ELSE
       HBU1=3.D0*VUB**2*CQCDM(AMCH,TB,(AMB/AMCH)**2,EPS)
       HBU2=3.D0*VUB**2*CQCD(AMCH,TB,(RMB/AMCH)**2,EPS)
       IF(HBU2.LT.0.D0) HBU2 = 0
       RAT = AMB/AMCH
       HBU = QQINT(RAT,HBU1,HBU2)
      ENDIF
C  H+ --> TB :
      IF(IONSH.EQ.0)THEN
       DLD=2D0
       DLU=2D0
       XM1 = AMT+AMB-DLD
       XM2 = AMT+AMB+DLU
       IF (AMCH.LE.AMW+2*AMB) THEN
        HBT=0.D0
       ELSEIF (AMCH.LE.XM1) THEN
        FACTB=3.D0*GF**2*AMCH*AMT**4/32.D0/PI**3/TB**2
        CALL CTOTT(AMCH,AMT,AMB,AMW,CTT0)
        HBT=FACTB*CTT0
       ELSEIF (AMCH.LE.XM2) THEN
        XX(1) = XM1-1D0
        XX(2) = XM1
        XX(3) = XM2
        XX(4) = XM2+1D0
        FACTB=3.D0*GF**2*XX(1)*AMT**4/32.D0/PI**3/TB**2
        CALL CTOTT(XX(1),AMT,AMB,AMW,CTT0)
        YY(1)=FACTB*CTT0
        FACTB=3.D0*GF**2*XX(2)*AMT**4/32.D0/PI**3/TB**2
        CALL CTOTT(XX(2),AMT,AMB,AMW,CTT0)
        YY(2)=FACTB*CTT0
        XMB = RUNM(XX(3),5)
        XMT = RUNM(XX(3),6)
        XYZ2 = 3.D0*CQCD(XX(3),TB,(XMB/XX(3))**2,(XMT/XX(3))**2)
        IF(XYZ2.LT.0.D0) XYZ2 = 0
        XYZ1 = 3.D0*CQCDM(XX(3),TB,(AMB/XX(3))**2,(AMT/XX(3))**2)
        RAT = (AMB+AMT)/XX(3)
        YY(3) = QQINT(RAT,XYZ1,XYZ2)
        XMB = RUNM(XX(4),5)
        XMT = RUNM(XX(4),6)
        XYZ2 = 3.D0*CQCD(XX(4),TB,(XMB/XX(4))**2,(XMT/XX(4))**2)
        IF(XYZ2.LT.0.D0) XYZ2 = 0
        XYZ1 = 3.D0*CQCDM(XX(4),TB,(AMB/XX(4))**2,(AMT/XX(4))**2)
        RAT = (AMB+AMT)/XX(4)
        YY(4) = QQINT(RAT,XYZ1,XYZ2)
        HBT = FINT(AMCH,XX,YY)
       ELSE
        HBT2=3.D0*CQCD(AMCH,TB,(RMB/AMCH)**2,(RMT/AMCH)**2)
        IF(HBT2.LT.0.D0) HBT2 = 0
        HBT1=3.D0*CQCDM(AMCH,TB,(AMB/AMCH)**2,(AMT/AMCH)**2)
        RAT = (AMB+AMT)/AMCH
        HBT = QQINT(RAT,HBT1,HBT2)
       ENDIF
      ELSE
       IF (AMCH.LE.AMT+AMB) THEN
        HBT=0.D0
       ELSE
        HBT2=3.D0*CQCD(AMCH,TB,(RMB/AMCH)**2,(RMT/AMCH)**2)
        IF(HBT2.LT.0.D0) HBT2 = 0
        HBT1=3.D0*CQCDM(AMCH,TB,(AMB/AMCH)**2,(AMT/AMCH)**2)
        RAT = (AMB+AMT)/AMCH
        HBT = QQINT(RAT,HBT1,HBT2)
       ENDIF
      ENDIF
C  H+ ---> W H
      IF(IONSH.EQ.0)THEN
       DLD=3D0
       DLU=5D0
       XM1 = AMW+AML-DLD
       XM2 = AMW+AML+DLU
       IF (AMCH.LT.AML) THEN
        HWH=0
       ELSEIF (AMCH.LE.XM1) THEN
        IF(AMCH.LE.DABS(AMW-AML))THEN
         HWH=0
        ELSE
         HWH=9.D0*GF**2/16.D0/PI**3*AMW**4*AMCH*GHVV**2	
     .      *HVH((AML/AMCH)**2,(AMW/AMCH)**2)
        ENDIF
       ELSEIF (AMCH.LT.XM2) THEN
        XX(1) = XM1-1D0
        XX(2) = XM1
        XX(3) = XM2
        XX(4) = XM2+1D0
        YY(1) = 9.D0*GF**2/16.D0/PI**3*AMW**4*XX(1)
     .         *HVH((AML/XX(1))**2,(AMW/XX(1))**2)
        YY(2) = 9.D0*GF**2/16.D0/PI**3*AMW**4*XX(2)
     .         *HVH((AML/XX(2))**2,(AMW/XX(2))**2)
        CWH=LAMB(AML**2/XX(3)**2,AMW**2/XX(3)**2)
     .     *LAMB(XX(3)**2/AMW**2,AML**2/AMW**2)**2
        YY(3)=GF/8.D0/DSQRT(2D0)/PI*AMW**4/XX(3)*CWH
        CWH=LAMB(AML**2/XX(4)**2,AMW**2/XX(4)**2)
     .     *LAMB(XX(4)**2/AMW**2,AML**2/AMW**2)**2
        YY(4)=GF/8.D0/DSQRT(2D0)/PI*AMW**4/XX(4)*CWH
        HWH = FINT(AMCH,XX,YY)*GHVV**2
       ELSE
        CWH=LAMB(AML**2/AMCH**2,AMW**2/AMCH**2)
     .     *LAMB(AMCH**2/AMW**2,AML**2/AMW**2)**2
        HWH=GF/8.D0/DSQRT(2D0)/PI*AMW**4/AMCH*GHVV**2*CWH
       ENDIF
      ELSE
       IF (AMCH.LT.AMW+AML) THEN
        HWH=0
       ELSE
        CWH=LAMB(AML**2/AMCH**2,AMW**2/AMCH**2)
     .     *LAMB(AMCH**2/AMW**2,AML**2/AMW**2)**2
        HWH=GF/8.D0/DSQRT(2D0)/PI*AMW**4/AMCH*GHVV**2*CWH
       ENDIF
      ENDIF
C  H+ ---> W A
      IF(IONSH.EQ.0)THEN
       IF (AMCH.LT.AMA) THEN
        HWA=0
       ELSEIF (AMCH.LT.AMW+AMA) THEN
        IF(AMCH.LE.DABS(AMW-AMA))THEN
         HWA=0
        ELSE
         HWA=9.D0*GF**2/16.D0/PI**3*AMW**4*AMCH	
     .      *HVH((AMA/AMCH)**2,(AMW/AMCH)**2)
        ENDIF
       ELSE
        HWA=0.D0
       ENDIF
      ELSE
       IF (AMCH.LT.AMW+AMA) THEN
        HWA=0
       ELSE
        HWA=0.D0
       ENDIF
      ENDIF

C  ======================= SUSY DECAYS 
      IF(IOFSUSY.EQ.0) THEN
C
C  H+ ----> CHARGINOS+NEUTRALINOS
C
      DO 751 I=1,2
      DO 751 J=1,4
      IF (AMCH.GT.AMCHAR(I)+AMNEUT(J)) THEN
      WHCCN(I,J)=GF*AMW**2/(2*PI*DSQRT(2.D0))/AMCH
     .   *LAMB(AMCHAR(I)**2/AMCH**2,AMNEUT(J)**2/AMCH**2)*(
     .   (ACNL(I,J)**2+ACNR(I,J)**2)*(AMCH**2-AMCHAR(I)**2-XMNEUT(J)
     .   **2)-4.D0*ACNL(I,J)*ACNR(I,J)*AMCHAR(I)*XMNEUT(J) )
      ELSE
      WHCCN(I,J)=0.D0
      ENDIF
 751  CONTINUE

      WHCCNT=WHCCN(1,1)+WHCCN(1,2)+WHCCN(1,3)+WHCCN(1,4)
     .      +WHCCN(2,1)+WHCCN(2,2)+WHCCN(2,3)+WHCCN(2,4)
C
C  H+ ----> SLEPTONS 
C
      IF (AMCH.GT.AMSE(1)+AMSN(1)) THEN
      WHCSL00=2*GF/4.D0/DSQRT(2D0)/PI*AMW**4/AMCH*DSIN(2.D0*B)**2
     .     *LAMB(AMSE(1)**2/AMCH**2,AMSN(1)**2/AMCH**2)
      ELSE 
      WHCSL00=0.D0
      ENDIF

      IF (AMCH.GT.AMSL(1)+AMSN(1)) THEN
      WHCSL11=GF/4.D0/DSQRT(2D0)/PI*AMW**4/AMCH*GCEN(1,1)**2
     .     *LAMB(AMSL(1)**2/AMCH**2,AMSN(1)**2/AMCH**2)
      ELSE 
      WHCSL11=0.D0
      ENDIF

      IF (AMCH.GT.AMSL(2)+AMSN(1)) THEN
      WHCSL21=GF/4.D0/DSQRT(2D0)/PI*AMW**4/AMCH*GCEN(1,2)**2
     .     *LAMB(AMSL(2)**2/AMCH**2,AMSN(1)**2/AMCH**2)
      ELSE 
      WHCSL21=0.D0
      ENDIF

      WHCSLT=WHCSL00+WHCSL11+WHCSL21

C
C  H+ ----> SQUARKS 
C
      IF (AMCH.GT.AMSU(1)+AMSD(1)) THEN
      WHCSQ=6*GF/4.D0/DSQRT(2D0)/PI*AMW**4/AMCH*DSIN(2.D0*B)**2
     .     *LAMB(AMSU(1)**2/AMCH**2,AMSD(1)**2/AMCH**2)
      ELSE 
      WHCSQ=0.D0
      ENDIF
C
      DO 753 I=1,2
      DO 753 J=1,2
      IF(AMCH.GT.AMST(I)+AMSB(J)) THEN
      WHCSTB(I,J)=3*GF*AMW**4/4.D0/DSQRT(2.D0)/PI*GCTB(I,J)**2
     .      *LAMB(AMST(I)**2/AMCH**2,AMSB(J)**2/AMCH**2)/AMCH
      ELSE
      WHCSTB(I,J)=0.D0
      ENDIF

 753  CONTINUE
C
      WHCSQT=WHCSQ+WHCSTB(1,1)+WHCSTB(1,2)+WHCSTB(2,1)+WHCSTB(2,2) 

      ELSE 
      WHCCNT=0.D0
      WHCSLT=0.D0
      WHCSQT=0.D0
C--Change thanks to Elzbieta Richter-Was
      DO I=1,2
       DO J=1,2
        WHCSTB(I,J)=0.D0
       ENDDO
      ENDDO
      DO I=1,2
       DO J=1,4
        WHCCN(I,J)=0.D0
       ENDDO
      ENDDO
      ENDIF
C
C    ==========  TOTAL WIDTH AND BRANCHING RATIOS 
C
      WTOT=HLN+HMN+HSU+HBU+HSC+HBC+HBT+HWH+HWA+WHCCNT+WHCSLT+WHCSQT

      HCBRL=HLN/WTOT
      HCBRM=HMN/WTOT
      HCBRS=HSU/WTOT
      HCBRBU=HBU/WTOT
      HCBRC=HSC/WTOT
      HCBRB=HBC/WTOT
      HCBRT=HBT/WTOT
      HCBRW=HWH/WTOT
      HCBRA=HWA/WTOT
      DO 851 I=1,2
      DO 851 J=1,4
      HCBRSU(I,J)=WHCCN(I,J)/WTOT
851   CONTINUE
      HCBRCNT=WHCCNT/WTOT
      HCBRSL=WHCSLT/WTOT 
      HCBRSQ=WHCSQ/WTOT 
      HCBRSQT=WHCSQT/WTOT 
      DO 853 I=1,2
      DO 853 J=1,2
      HCBRSTB(I,J)=WHCSTB(I,J)/WTOT
853   CONTINUE
      HCWDTH=WTOT

      BHCSL00 = WHCSL00/WTOT
      BHCSL11 = WHCSL11/WTOT
      BHCSL21 = WHCSL21/WTOT
      BHCSQ = WHCSQ/WTOT
      DO I = 1,2
       DO J = 1,2
        BHCSTB(I,J) = WHCSTB(I,J)/WTOT
       ENDDO
      ENDDO

      GAMC0 = WTOT

      ENDIF

      IF(IHIGGS.EQ.2.OR.IHIGGS.EQ.5)THEN
     
C        =========================================================
C                       HEAVY CP EVEN HIGGS DECAYS
C        =========================================================
C     =============  RUNNING MASSES 
      RMS = RUNM(AMH,3)
      RMC = RUNM(AMH,4)
      RMB = RUNM(AMH,5)
      RMT = RUNM(AMH,6)
      RATCOUP = GHT/GHB
      HIGTOP = AMH**2/AMT**2

      ASH=ALPHAS(AMH,2)
      AMC0=1.D8
      AMB0=2.D8
C     AMT0=3.D8
      AS3=ALPHAS(AMH,2)
      AMC0=AMC
      AS4=ALPHAS(AMH,2)
      AMB0=AMB
C     AMT0=AMT

C     =============== PARTIAL WIDTHS 
C  H ---> G G
       EPS=1.D-8
       NFEXT = 3
       ASG = AS3
       CTT = 4*AMT**2/AMH**2*DCMPLX(1D0,-EPS)
       CTB = 4*AMB**2/AMH**2*DCMPLX(1D0,-EPS)
       CAT = 2*CTT*(1+(1-CTT)*CF(CTT))*GHT
       CAB = 2*CTB*(1+(1-CTB)*CF(CTB))*GHB
       CTC = 4*AMC**2/AMH**2*DCMPLX(1D0,-EPS)
       CAC = 2*CTC*(1+(1-CTC)*CF(CTC))*GHT
C
       IF(IOFSUSY.EQ.0) THEN 
       CSB1= 4*AMSB(1)**2/AMH**2*DCMPLX(1D0,-EPS)
       CSB2= 4*AMSB(2)**2/AMH**2*DCMPLX(1D0,-EPS)
       CST1= 4*AMST(1)**2/AMH**2*DCMPLX(1D0,-EPS)
       CST2= 4*AMST(2)**2/AMH**2*DCMPLX(1D0,-EPS)
C
       CXB1=-AMZ**2/AMSB(1)**2*CSB1*(1-CSB1*CF(CSB1))*GHBB(1,1)
       CXB2=-AMZ**2/AMSB(2)**2*CSB2*(1-CSB2*CF(CSB2))*GHBB(2,2)
       CXT1=-AMZ**2/AMST(1)**2*CST1*(1-CST1*CF(CST1))*GHTT(1,1)
       CXT2=-AMZ**2/AMST(2)**2*CST2*(1-CST2*CF(CST2))*GHTT(2,2)
C
       CSUL = 4*AMSU(1)**2/AMH**2*DCMPLX(1D0,-EPS)
       CSUR = 4*AMSU(2)**2/AMH**2*DCMPLX(1D0,-EPS)
       CSDL = 4*AMSD(1)**2/AMH**2*DCMPLX(1D0,-EPS)
       CSDR = 4*AMSD(2)**2/AMH**2*DCMPLX(1D0,-EPS)
       CXUL=-2*(1.D0/2.D0-2.D0/3.D0*SS)*AMZ**2/AMSU(1)**2*DCOS(A+B)
     .      *CSUL*(1-CSUL*CF(CSUL))
       CXUR=-2*(2.D0/3.D0*SS)*AMZ**2/AMSU(2)**2*DCOS(A+B)
     .      *CSUR*(1-CSUR*CF(CSUR))
       CXDL=-2*(-1.D0/2.D0+1.D0/3.D0*SS)*AMZ**2/AMSD(1)**2*DCOS(A+B)
     .      *CSDL*(1-CSDL*CF(CSDL))
       CXDR=-2*(-1.D0/3.D0*SS)*AMZ**2/AMSD(2)**2*DCOS(A+B)
     .      *CSDR*(1-CSDR*CF(CSUR))
       ELSE
       CXB1=0.D0 
       CXB2=0.D0 
       CXT1=0.D0 
       CXT2=0.D0 
       ENDIF

       FQCD=HGGQCD(ASG,NFEXT)
       SQCD=SGGQCD(ASG)
       XFAC = CDABS(CAT+CAB+CAC+CXB1+CXB2+CXT1+CXT2
     .             +CXUL+CXUR+CXDL+CXDR)**2*FQCD
     .      + DREAL(DCONJG(CAT+CAB+CAC)*(CXB1+CXB2+CXT1+CXT2
     .             +CXUL+CXUR+CXDL+CXDR))*SQCD
       HGG=HVV(AMH,0.D0)*(ASG/PI)**2*XFAC/8

C  H ---> G G* ---> G CC   TO BE ADDED TO H ---> CC
       NFEXT = 4
       ASG = AS4
       FQCD=HGGQCD(ASG,NFEXT)
       SQCD=SGGQCD(ASG)
       XFAC = CDABS(CAT+CAB+CAC+CXB1+CXB2+CXT1+CXT2
     .             +CXUL+CXUR+CXDL+CXDR)**2*FQCD
     .      + DREAL(DCONJG(CAT+CAB+CAC)*(CXB1+CXB2+CXT1+CXT2
     .             +CXUL+CXUR+CXDL+CXDR))*SQCD
       DCC=HVV(AMH,0.D0)*(ASG/PI)**2*XFAC/8 - HGG

C  H ---> G G* ---> G BB   TO BE ADDED TO H ---> BB
       NFEXT = 5
       ASG = ASH
       FQCD=HGGQCD(ASG,NFEXT)
       SQCD=SGGQCD(ASG)
       XFAC = CDABS(CAT+CAB+CAC+CXB1+CXB2+CXT1+CXT2
     .             +CXUL+CXUR+CXDL+CXDR)**2*FQCD
     .      + DREAL(DCONJG(CAT+CAB+CAC)*(CXB1+CXB2+CXT1+CXT2
     .             +CXUL+CXUR+CXDL+CXDR))*SQCD
       DBB=HVV(AMH,0.D0)*(ASG/PI)**2*XFAC/8 - HGG - DCC

      IF(NFGG.EQ.5)THEN
       HGG = HGG + DBB + DCC
       DBB = 0
       DCC = 0
      ELSEIF(NFGG.EQ.4)THEN
       HGG = HGG + DCC
       DCC = 0
      ENDIF

C  H ---> MU MU
      IF(AMH.LE.2*AMMUON) THEN
       HMM = 0
      ELSE
      HMM=HFF(AMH,(AMMUON/AMH)**2)*GHB**2
      ENDIF
C  H ---> LL
      IF(AMH.LE.2*AMTAU) THEN
       HLL = 0
      ELSE
      HLL=HFF(AMH,(AMTAU/AMH)**2)*GHB**2
      ENDIF
C  H --> SS
      IF(AMH.LE.2*AMS) THEN
       HSS = 0
      ELSE
       HS1=3.D0*HFF(AMH,(AMS/AMH)**2)
     .    *GHB**2
     .    *TQCDH(AMS**2/AMH**2)
       HS2=3.D0*HFF(AMH,(RMS/AMH)**2)*GHB**2
     .    *QCDH(RMS**2/AMH**2)
       IF(HS2.LT.0.D0) HS2 = 0
       RAT = 2*AMS/AMH
       HSS = QQINT(RAT,HS1,HS2)
      ENDIF
C  H --> CC
      RATCOUP = 1
      IF(AMH.LE.2*AMC) THEN
       HCC = 0
      ELSE
       HC1=3.D0*HFF(AMH,(AMC/AMH)**2)
     .    *GHT**2
     .    *TQCDH(AMC**2/AMH**2)
       HC2=3.D0*HFF(AMH,(RMC/AMH)**2)*GHT**2
     .    *QCDH(RMC**2/AMH**2)
     .   + DCC
       IF(HC2.LT.0.D0) HC2 = 0
       RAT = 2*AMC/AMH
       HCC = QQINT(RAT,HC1,HC2)
      ENDIF
C  H --> BB :
      RATCOUP = GHT/GHB
      IF(AMH.LE.2*AMB) THEN
       HBB = 0
      ELSE
       HB1=3.D0*HFF(AMH,(AMB/AMH)**2)
     .    *GHB**2
     .    *TQCDH(AMB**2/AMH**2)
       HB2=3.D0*HFF(AMH,(RMB/AMH)**2)*GHB**2
     .    *QCDH(RMB**2/AMH**2)
     .   + DBB
       IF(HB2.LT.0.D0) HB2 = 0
       RAT = 2*AMB/AMH
       HBB = QQINT(RAT,HB1,HB2)
      ENDIF
C  H ---> TT
      RATCOUP = 0
      IF(IONSH.EQ.0)THEN
       DLD=3D0
       DLU=5D0
       XM1 = 2D0*AMT-DLD
       XM2 = 2D0*AMT+DLU
       IF (AMH.LE.AMT+AMW+AMB) THEN
        HTT=0.D0
       ELSEIF (AMH.LE.XM1) THEN
        FACTT=6.D0*GF**2*AMH**3*AMT**2/2.D0/128.D0/PI**3
        CALL HTOTT(AMH,AMT,AMB,AMW,AMCH,TB,GHT,GAT,GHVV,HTT0)
        HTT=FACTT*HTT0
       ELSEIF (AMH.LE.XM2) THEN
        ZZMA=AMAR
        XX(1) = XM1-1D0
        XX(2) = XM1
        XX(3) = XM2
        XX(4) = XM2+1D0
        CALL AMHAMA(2,XX(1),TGBET)
        FACTT=6.D0*GF**2*XX(1)**3*AMT**2/2.D0/128.D0/PI**3
        CALL HTOTT(XX(1),AMT,AMB,AMW,AMCH,TB,GHT,GAT,GHVV,HTT0)
        YY(1)=FACTT*HTT0
        CALL AMHAMA(2,XX(2),TGBET)
        FACTT=6.D0*GF**2*XX(2)**3*AMT**2/2.D0/128.D0/PI**3
        CALL HTOTT(XX(2),AMT,AMB,AMW,AMCH,TB,GHT,GAT,GHVV,HTT0)
        YY(2)=FACTT*HTT0
        CALL AMHAMA(2,XX(3),TGBET)
        XMT = RUNM(XX(3),6)
        HT1=3.D0*HFF(XX(3),(AMT/XX(3))**2)*GHT**2
     .    *TQCDH(AMT**2/XX(3)**2)
        HT2=3.D0*HFF(XX(3),(XMT/XX(3))**2)*GHT**2
     .    *QCDH(XMT**2/XX(3)**2)
        IF(HT2.LT.0.D0) HT2 = 0
        RAT = 2*AMT/XX(3)
        YY(3) = QQINT(RAT,HT1,HT2)
        CALL AMHAMA(2,XX(4),TGBET)
        XMT = RUNM(XX(4),6)
        HT1=3.D0*HFF(XX(4),(AMT/XX(4))**2)*GHT**2
     .    *TQCDH(AMT**2/XX(4)**2)
        HT2=3.D0*HFF(XX(4),(XMT/XX(4))**2)*GHT**2
     .    *QCDH(XMT**2/XX(4)**2)
        IF(HT2.LT.0.D0) HT2 = 0
        RAT = 2*AMT/XX(4)
        YY(4) = QQINT(RAT,HT1,HT2)
        AMA = ZZMA
        CALL SUSYCP(TGBET)
        HTT=FINT(AMH,XX,YY)
       ELSE
        HT1=3.D0*HFF(AMH,(AMT/AMH)**2)*GHT**2
     .    *TQCDH(AMT**2/AMH**2)
        HT2=3.D0*HFF(AMH,(RMT/AMH)**2)*GHT**2
     .    *QCDH(RMT**2/AMH**2)
        IF(HT2.LT.0.D0) HT2 = 0
        RAT = 2*AMT/AMH
        HTT = QQINT(RAT,HT1,HT2)
       ENDIF
      ELSE
       IF (AMH.LE.2.D0*AMT) THEN
        HTT=0.D0
       ELSE
        HT1=3.D0*HFF(AMH,(AMT/AMH)**2)*GHT**2
     .    *TQCDH(AMT**2/AMH**2)
        HT2=3.D0*HFF(AMH,(RMT/AMH)**2)*GHT**2
     .    *QCDH(RMT**2/AMH**2)
        IF(HT2.LT.0.D0) HT2 = 0
        RAT = 2*AMT/AMH
        HTT = QQINT(RAT,HT1,HT2)
       ENDIF
      ENDIF
C  H ---> GAMMA GAMMA
       EPS=1.D-8
       CTT = 4*AMT**2/AMH**2*DCMPLX(1D0,-EPS)
       CTB = 4*AMB**2/AMH**2*DCMPLX(1D0,-EPS)
       CTL = 4*AMTAU**2/AMH**2*DCMPLX(1D0,-EPS)
       CTW = 4*AMW**2/AMH**2*DCMPLX(1D0,-EPS)
       CTH = 4*AMCH**2/AMH**2*DCMPLX(1D0,-EPS)
       CTC = 4*AMC**2/AMH**2*DCMPLX(1D0,-EPS)
       CAC = 4/3D0 * 2*CTC*(1+(1-CTC)*CF(CTC))*GHT
       CAT = 4/3D0 * 2*CTT*(1+(1-CTT)*CF(CTT))*GHT
       CAB = 1/3D0 * 2*CTB*(1+(1-CTB)*CF(CTB))*GHB
       CAL = 1.D0  * 2*CTL*(1+(1-CTL)*CF(CTL))*GHB
       CAW = -(2+3*CTW+3*CTW*(2-CTW)*CF(CTW))*GHVV
       CAH = -AMZ**2/2/AMCH**2*CTH*(1-CTH*CF(CTH))*GHPM
       IF(IOFSUSY.EQ.0) THEN 
       CX1 = 4*AMCHAR(1)**2/AMH**2*DCMPLX(1D0,-EPS)
       CX2 = 4*AMCHAR(2)**2/AMH**2*DCMPLX(1D0,-EPS)
       CAX1= AMW/AMCHAR(1) * 2*CX1*(1+(1-CX1)*CF(CX1))*2*AC1(1,1) 
       CAX2= AMW/AMCHAR(2) * 2*CX2*(1+(1-CX2)*CF(CX2))*2*AC1(2,2) 
       CSL1= 4*AMSL(1)**2/AMH**2*DCMPLX(1D0,-EPS)
       CSL2= 4*AMSL(2)**2/AMH**2*DCMPLX(1D0,-EPS)
       CSB1= 4*AMSB(1)**2/AMH**2*DCMPLX(1D0,-EPS)
       CSB2= 4*AMSB(2)**2/AMH**2*DCMPLX(1D0,-EPS)
       CST1= 4*AMST(1)**2/AMH**2*DCMPLX(1D0,-EPS)
       CST2= 4*AMST(2)**2/AMH**2*DCMPLX(1D0,-EPS)

       CSEL = 4*AMSE(1)**2/AMH**2*DCMPLX(1D0,-EPS)
       CSER = 4*AMSE(2)**2/AMH**2*DCMPLX(1D0,-EPS)
       CSUL = 4*AMSU(1)**2/AMH**2*DCMPLX(1D0,-EPS)
       CSUR = 4*AMSU(2)**2/AMH**2*DCMPLX(1D0,-EPS)
       CSDL = 4*AMSD(1)**2/AMH**2*DCMPLX(1D0,-EPS)
       CSDR = 4*AMSD(2)**2/AMH**2*DCMPLX(1D0,-EPS)
       CXEL=-2*(-1/2D0+SS)*AMZ**2/AMSE(1)**2*DCOS(A+B)
     .      *CSEL*(1-CSEL*CF(CSEL))
       CXER=2*(SS)*AMZ**2/AMSE(2)**2*DCOS(A+B)
     .      *CSER*(1-CSER*CF(CSER))
       CXUL=-2*4.D0/3.D0*(1.D0/2.D0-2.D0/3.D0*SS)
     .      *AMZ**2/AMSU(1)**2*DCOS(A+B)*CSUL*(1-CSUL*CF(CSUL))
       CXUR=-2*4.D0/3.D0*(2.D0/3.D0*SS)
     .      *AMZ**2/AMSU(2)**2*DCOS(A+B)*CSUR*(1-CSUR*CF(CSUR))
       CXDL=-2/3.D0*(-1.D0/2.D0+1.D0/3.D0*SS)
     .      *AMZ**2/AMSD(1)**2*DCOS(A+B)*CSDL*(1-CSDL*CF(CSDL))
       CXDR=-2/3.D0*(-1.D0/3.D0*SS)
     .      *AMZ**2/AMSD(2)**2*DCOS(A+B)*CSDR*(1-CSDR*CF(CSUR))

       CXB1= -1/3D0*AMZ**2/AMSB(1)**2*CSB1*(1-CSB1*CF(CSB1))*GHBB(1,1)
       CXB2= -1/3D0*AMZ**2/AMSB(2)**2*CSB2*(1-CSB2*CF(CSB2))*GHBB(2,2)
       CXT1= -4/3D0*AMZ**2/AMST(1)**2*CST1*(1-CST1*CF(CST1))*GHTT(1,1)
       CXT2= -4/3D0*AMZ**2/AMST(2)**2*CST2*(1-CST2*CF(CST2))*GHTT(2,2)
       CXL1=       -AMZ**2/AMSL(1)**2*CSL1*(1-CSL1*CF(CSL1))*GHEE(1,1)
       CXL2=       -AMZ**2/AMSL(2)**2*CSL2*(1-CSL2*CF(CSL2))*GHEE(2,2)
       XFAC = CDABS(CAT+CAB+CAC+CAL+CAW+CAH+CAX1+CAX2
     .     +  CXEL+CXER+CXUL+CXUR+CXDL+CXDR
     .     +  CXB1+CXB2+CXT1+CXT2+CXL1+CXL2)**2
       ELSE 
       XFAC = CDABS(CAT+CAB+CAC+CAL+CAW+CAH)**2
       ENDIF
       HGA=HVV(AMH,0.D0)*(ALPH/PI)**2/16.D0*XFAC
C  H ---> Z GAMMA
      IF(AMH.LE.AMZ)THEN
       HZGA=0
      ELSE
       TS = SS/CS
       FT = -3*2D0/3*(1-4*2D0/3*SS)/DSQRT(SS*CS)*GHT
       FB = 3*1D0/3*(-1+4*1D0/3*SS)/DSQRT(SS*CS)*GHB
       EPS=1.D-8
       CTT = 4*AMT**2/AMH**2*DCMPLX(1D0,-EPS)
       CTB = 4*AMB**2/AMH**2*DCMPLX(1D0,-EPS)
       CTW = 4*AMW**2/AMH**2*DCMPLX(1D0,-EPS)
       CTH = 4*AMCH**2/AMH**2*DCMPLX(1D0,-EPS)
       CLT = 4*AMT**2/AMZ**2*DCMPLX(1D0,-EPS)
       CLB = 4*AMB**2/AMZ**2*DCMPLX(1D0,-EPS)
       CLW = 4*AMW**2/AMZ**2*DCMPLX(1D0,-EPS)
       CLH = 4*AMCH**2/AMZ**2*DCMPLX(1D0,-EPS)
       CAT = FT*(CI1(CTT,CLT) - CI2(CTT,CLT))
       CAB = FB*(CI1(CTB,CLB) - CI2(CTB,CLB))
       CAW = -1/DSQRT(TS)*(4*(3-TS)*CI2(CTW,CLW)
     .     + ((1+2/CTW)*TS - (5+2/CTW))*CI1(CTW,CLW))*GHVV
       CAH = (1-2*SS)/DSQRT(SS*CS)*AMZ**2/2/AMCH**2*CI1(CTH,CLH)*GHPM
       XFAC = CDABS(CAT+CAB+CAW+CAH)**2
       ACOUP = DSQRT(2D0)*GF*AMZ**2*SS*CS/PI**2
       HZGA = GF/(4.D0*PI*DSQRT(2.D0))*AMH**3*(ALPH/PI)*ACOUP/16.D0
     .        *XFAC*(1-AMZ**2/AMH**2)**3
      ENDIF
C  H ---> W W
      IF(IONWZ.EQ.0)THEN
       DLD=2D0
       DLU=2D0
       XM1 = 2D0*AMW-DLD
       XM2 = 2D0*AMW+DLU
       IF (AMH.LE.XM1) THEN
        CALL HTOVV(AMH,AMW,GAMW,HTWW)
        HWW = 3D0/2D0*GF*AMW**4/DSQRT(2D0)/PI/AMH**3*HTWW*GHVV**2
       ELSEIF (AMH.LE.XM2) THEN
        XX(1) = XM1-1D0
        XX(2) = XM1
        XX(3) = XM2
        XX(4) = XM2+1D0
        CALL HTOVV(XX(1),AMW,GAMW,HTWW)
        YY(1)=3D0/2D0*GF*AMW**4/DSQRT(2D0)/PI/XX(1)**3*HTWW
        CALL HTOVV(XX(2),AMW,GAMW,HTWW)
        YY(2)=3D0/2D0*GF*AMW**4/DSQRT(2D0)/PI/XX(2)**3*HTWW
        YY(3)=HVV(XX(3),AMW**2/XX(3)**2)
        YY(4)=HVV(XX(4),AMW**2/XX(4)**2)
        HWW = FINT(AMH,XX,YY)*GHVV**2
       ELSE
        HWW=HVV(AMH,AMW**2/AMH**2)*GHVV**2
       ENDIF
      ELSE
      DLD=2D0
      DLU=2D0
      XM1 = 2D0*AMW-DLD
      XM2 = 2D0*AMW+DLU
      IF (AMH.LE.AMW) THEN
       HWW=0
      ELSE IF (AMH.LE.XM1) THEN
       CWW=3.D0*GF**2*AMW**4/16.D0/PI**3
       HWW=HV(AMW**2/AMH**2)*CWW*AMH*GHVV**2
      ELSE IF (AMH.LT.XM2) THEN
       CWW=3.D0*GF**2*AMW**4/16.D0/PI**3
       XX(1) = XM1-1D0
       XX(2) = XM1
       XX(3) = XM2
       XX(4) = XM2+1D0
       YY(1)=HV(AMW**2/XX(1)**2)*CWW*XX(1)
       YY(2)=HV(AMW**2/XX(2)**2)*CWW*XX(2)
       YY(3)=HVV(XX(3),AMW**2/XX(3)**2)
       YY(4)=HVV(XX(4),AMW**2/XX(4)**2)
       HWW = FINT(AMH,XX,YY)*GHVV**2
      ELSE
       HWW=HVV(AMH,AMW**2/AMH**2)*GHVV**2
      ENDIF
      ENDIF
C  H ---> Z Z
      IF(IONWZ.EQ.0)THEN
       DLD=2D0
       DLU=2D0
       XM1 = 2D0*AMZ-DLD
       XM2 = 2D0*AMZ+DLU
       IF (AMH.LE.XM1) THEN
        CALL HTOVV(AMH,AMZ,GAMZ,HTZZ)
        HZZ = 3D0/4D0*GF*AMZ**4/DSQRT(2D0)/PI/AMH**3*HTZZ*GHVV**2
       ELSEIF (AMH.LE.XM2) THEN
        XX(1) = XM1-1D0
        XX(2) = XM1
        XX(3) = XM2
        XX(4) = XM2+1D0
        CALL HTOVV(XX(1),AMZ,GAMZ,HTZZ)
        YY(1)=3D0/4D0*GF*AMZ**4/DSQRT(2D0)/PI/XX(1)**3*HTZZ
        CALL HTOVV(XX(2),AMZ,GAMZ,HTZZ)
        YY(2)=3D0/4D0*GF*AMZ**4/DSQRT(2D0)/PI/XX(2)**3*HTZZ
        YY(3)=HVV(XX(3),AMZ**2/XX(3)**2)/2
        YY(4)=HVV(XX(4),AMZ**2/XX(4)**2)/2
        HZZ = FINT(AMH,XX,YY)*GHVV**2
       ELSE
        HZZ=HVV(AMH,AMZ**2/AMH**2)/2.D0*GHVV**2
       ENDIF
      ELSE
      DLD=2D0
      DLU=2D0
      XM1 = 2D0*AMZ-DLD
      XM2 = 2D0*AMZ+DLU
      IF (AMH.LE.AMZ) THEN
       HZZ=0
      ELSE IF (AMH.LE.XM1) THEN
       CZZ=3.D0*GF**2*AMZ**4/192.D0/PI**3*(7-40/3.D0*SS+160/9.D0*SS**2)
       HZZ=HV(AMZ**2/AMH**2)*CZZ*AMH*GHVV**2
      ELSE IF (AMH.LT.XM2) THEN
       CZZ=3.D0*GF**2*AMZ**4/192.D0/PI**3*(7-40/3.D0*SS+160/9.D0*SS**2)
       XX(1) = XM1-1D0
       XX(2) = XM1
       XX(3) = XM2
       XX(4) = XM2+1D0
       YY(1)=HV(AMZ**2/XX(1)**2)*CZZ*XX(1)
       YY(2)=HV(AMZ**2/XX(2)**2)*CZZ*XX(2)
       YY(3)=HVV(XX(3),AMZ**2/XX(3)**2)/2D0
       YY(4)=HVV(XX(4),AMZ**2/XX(4)**2)/2D0
       HZZ = FINT(AMH,XX,YY)*GHVV**2
      ELSE
       HZZ=HVV(AMH,AMZ**2/AMH**2)/2.D0*GHVV**2
      ENDIF
      ENDIF
C  H ---> h h
      IF(IONSH.EQ.0)THEN
       ZZMA = AMAR
       AMREAL = AMH
       AMA = 20.D0
       CALL SUSYCP(TGBET)
       AMLOW = AMH
       AMDEL = AMREAL - AMLOW
       DLD = 0.3D0*(TGBET-1.3D0)
       DLD = DMAX1(0.1D0,DLD)
       DLU=DLD
       AMA = ZZMA
       CALL SUSYCP(TGBET)
       XM1 = 2*AML-DLD
       XM2 = 2*AML+DLU
       IF (AMH.LE.AML) THEN
        HHH=0
       ELSEIF (AMH.LT.XM1) THEN
        XH=AML**2/AMH**2
        XH1=(XH-1.D0)*(2.D0-.5D0*DLOG(XH))+(1.D0-5.D0*XH)
     .    *(DATAN((2.D0*XH-1.D0)/DSQRT(4.D0*XH-1.D0))
     .     -DATAN(1.D0/DSQRT(4.D0*XH-1.D0)))/DSQRT(4.D0*XH-1.D0)
        XH2=3*GF**2/16.D0/PI**3*AMZ**4/AMH*GHLL**2*GLB**2*AMB**2
        HHH=XH1*XH2
       ELSEIF (AMH.LT.XM2) THEN
        IFLON0 = 0
        IFLON1 = 0
        ZZMA=AMAR
        AMACRIT = AMAR
        AMA0 = AMAR
        AMA1 = AMAR
510     AMA0 = AMA0 - 1
        AMA1 = AMA1 + 1
        AMA = AMA0
        CALL SUSYCP(TGBET)
        IF(AMH.LT.2*AML) THEN
         IFLON0 = -1
        ELSE
         IFLON0 = 1
        ENDIF
        AMA = AMA1
        CALL SUSYCP(TGBET)
        IF(AMH.LT.2*AML) THEN
         IFLON1 = -1
        ELSE
         IFLON1 = 1
        ENDIF
        IF(IFLON0*IFLON1.NE.-1) GOTO 510
501     AMA = (AMA0+AMA1)/2
        CALL SUSYCP(TGBET)
        IF(AMH.LT.2*AML) THEN
         IF(IFLON0.EQ.-1) THEN
          AMA0 = AMAR
         ELSE
          AMA1 = AMAR
         ENDIF
        ELSE
         IF(IFLON0.EQ.-1) THEN
          AMA1 = AMAR
         ELSE
          AMA0 = AMAR
         ENDIF
        ENDIF
        AMACRIT = (AMA0+AMA1)/2
        DEL = 1.D-8
        AMDEL = 2*DABS(AMA1-AMA0)/(AMA1+AMA0)
        IF(AMDEL.GT.DEL) GOTO 501
       AMA = AMACRIT
       CALL SUSYCP(TGBET)
       YM1 = AMACRIT
       YM2 = AMACRIT
       AMA0 = AMACRIT
       AMA1 = AMACRIT
       DELSTEP = 1.D0
511    AMA0 = AMA0 - DELSTEP
       AMA1 = AMA1 + DELSTEP
       AMA = AMACRIT
       CALL SUSYCP(TGBET)
       IF(AMH.LT.2*AML-DLD) THEN
        IFLONC = -1
       ELSE
        IFLONC = 1
       ENDIF
       AMA = AMA0
       CALL SUSYCP(TGBET)
       IF(AMH.LT.2*AML-DLD) THEN
        IFLON0 = -1
       ELSE
        IFLON0 = 1
       ENDIF
       AMA = AMA1
       CALL SUSYCP(TGBET)
       IF(AMH.LT.2*AML-DLD) THEN
        IFLON1 = -1
       ELSE
        IFLON1 = 1
       ENDIF
       IF(IFLON0*IFLONC.NE.-1.AND.IFLONC*IFLON1.NE.-1) GOTO 511
       IF(IFLON0*IFLONC.EQ.-1) THEN
         AMA1 = AMACRIT
         IFLON1 = IFLONC
       ELSE
         AMA0 = AMACRIT
         IFLON0 = IFLONC
       ENDIF
512    AMA = (AMA0+AMA1)/2
       CALL SUSYCP(TGBET)
       IF(AMH.LT.2*AML-DLD) THEN
        IF(IFLON0.EQ.-1) THEN
         AMA0 = AMAR
        ELSE
         AMA1 = AMAR
        ENDIF
       ELSE
        IF(IFLON0.EQ.-1) THEN
         AMA1 = AMAR
        ELSE
         AMA0 = AMAR
        ENDIF
       ENDIF
       YM1 = (AMA0+AMA1)/2
       DEL = 1.D-8
       AMDEL = 2*DABS(AMA1-AMA0)/(AMA1+AMA0)
       IF(AMDEL.GT.DEL) GOTO 512
       AMA = YM1
       CALL SUSYCP(TGBET)
       AMA0 = AMACRIT
       AMA1 = AMACRIT
       DELSTEP = 1.D0
513    AMA0 = AMA0 - DELSTEP
       AMA1 = AMA1 + DELSTEP
       AMA = AMACRIT
       CALL SUSYCP(TGBET)
       IF(AMH.LT.2*AML+DLU) THEN
        IFLONC = -1
       ELSE
        IFLONC = 1
       ENDIF
       AMA = AMA0
       CALL SUSYCP(TGBET)
       IF(AMH.LT.2*AML+DLU) THEN
        IFLON0 = -1
       ELSE
        IFLON0 = 1
       ENDIF
       AMA = AMA1
       CALL SUSYCP(TGBET)
       IF(AMH.LT.2*AML+DLU) THEN
        IFLON1 = -1
       ELSE
        IFLON1 = 1
       ENDIF
       IF(IFLON0*IFLONC.NE.-1.AND.IFLONC*IFLON1.NE.-1) GOTO 513
       IF(IFLON0*IFLONC.EQ.-1) THEN
         AMA1 = AMACRIT
         IFLON1 = IFLONC
       ELSE
         AMA0 = AMACRIT
         IFLON0 = IFLONC
       ENDIF
514    AMA = (AMA0+AMA1)/2
       CALL SUSYCP(TGBET)
       IF(AMH.LT.2*AML+DLU) THEN
        IF(IFLON0.EQ.-1) THEN
         AMA0 = AMAR
        ELSE
         AMA1 = AMAR
        ENDIF
       ELSE
        IF(IFLON0.EQ.-1) THEN
         AMA1 = AMAR
        ELSE
         AMA0 = AMAR
        ENDIF
       ENDIF
       YM2 = (AMA0+AMA1)/2
       DEL = 1.D-8
       AMDEL = 2*DABS(AMA1-AMA0)/(AMA1+AMA0)
       IF(AMDEL.GT.DEL) GOTO 514
       AMA = YM2
       CALL SUSYCP(TGBET)
       DEL = 1.D-4
        XX(1) = YM1 - DEL
        XX(2) = YM1
        XX(3) = YM2
        XX(4) = YM2 + DEL
        AMAR = ZZMA
        DO J=1,4
         AMA = XX(J)
         CALL SUSYCP(TGBET)
         XX(J) = AMH
         IF(AMH.GE.2*AML)THEN
          YY(J)=GF/16D0/DSQRT(2D0)/PI*AMZ**4/XX(J)
     .          *BETA(AML**2/XX(J)**2)
         ELSEIF(AMH.LE.AML)THEN
          YY(J) = 0
         ELSE
          XH=AML**2/XX(J)**2
          XH1=(XH-1.D0)*(2.D0-.5D0*DLOG(XH))+(1.D0-5.D0*XH)
     .    *(DATAN((2.D0*XH-1.D0)/DSQRT(4.D0*XH-1.D0))
     .     -DATAN(1.D0/DSQRT(4.D0*XH-1.D0)))/DSQRT(4.D0*XH-1.D0)
          XH2=3*GF**2/16.D0/PI**3*AMZ**4/XX(J)*GLB**2*AMB**2
          YY(J)=XH1*XH2
         ENDIF
        ENDDO
        AMA = ZZMA
        CALL SUSYCP(TGBET)
        HHH = FINT(AMH,XX,YY)*GHLL**2
       ELSE
        HHH=GF/16D0/DSQRT(2D0)/PI*AMZ**4/AMH*BETA(AML**2/AMH**2)
     .      *GHLL**2
       ENDIF
      ELSE
       IF (AMH.LE.2*AML) THEN
        HHH=0
       ELSE
        HHH=GF/16D0/DSQRT(2D0)/PI*AMZ**4/AMH*BETA(AML**2/AMH**2)
     .      *GHLL**2
       ENDIF
      ENDIF
C  H ---> A A
      IF(IONSH.EQ.0)THEN
       DLD = 0.3D0*(TGBET-1.3D0)
       DLD = DMAX1(0.1D0,DLD)
       DLU=DLD
       ALD = DLD/2
       ALU = DLU/2
       XM1 = 2*AMA-DLD
       XM2 = 2*AMA+DLU
       IF (AMH.LE.AMA) THEN
        HAA=0
       ELSEIF (AMH.LT.XM1) THEN
        XA=AMA**2/AMH**2
        XA1=(XA-1.D0)*(2.D0-.5D0*DLOG(XA))+(1.D0-5.D0*XA)
     .    *(DATAN((2.D0*XA-1.D0)/DSQRT(4.D0*XA-1.D0))
     .     -DATAN(1.D0/DSQRT(4.D0*XA-1.D0)))/DSQRT(4.D0*XA-1.D0)
        XA2=3*GF**2/16.D0/PI**3*AMZ**4/AMH*GHAA**2*GAB**2*AMB**2
        HAA=XA1*XA2
       ELSEIF (AMH.LT.XM2) THEN
        ZZMA=AMAR
        AMACRIT = AMAR
        AMA0 = 10.D0
        AMA1 = AMAR + 50.D0
        AMA = AMA0
        CALL SUSYCP(TGBET)
        IF(AMH.LT.2*AMA) THEN
         IFLON0 = -1
        ELSEIF(AMH.EQ.2*AMA) THEN
         IFLON0 = 0
         AMACRIT = AMAR
        ELSE
         IFLON0 = 1
        ENDIF
        AMA = AMA1
        CALL SUSYCP(TGBET)
        IF(AMH.LT.2*AMA) THEN
         IFLON1 = -1
        ELSEIF(AMH.EQ.2*AMA) THEN
         IFLON1 = 0
         AMACRIT = AMAR
        ELSE
         IFLON1 = 1
        ENDIF
        IF(IFLON0*IFLON1.EQ.0)THEN
         IFLON0 = 0
         IFLON1 = 0
        ENDIF
        IF(IFLON0.NE.IFLON1)THEN
502      AMA = (AMA0+AMA1)/2
         CALL SUSYCP(TGBET)
         IF(AMH.LT.2*AMA) THEN
          IF(IFLON0.EQ.-1) THEN
           AMA0 = AMAR
          ELSE
           AMA1 = AMAR
          ENDIF
         ELSEIF(AMH.EQ.2*AMA) THEN
          IFLON0 = 0
          IFLON1 = 0
          AMACRIT = AMAR
         ELSE
          IF(IFLON0.EQ.-1) THEN
           AMA1 = AMAR
          ELSE
           AMA0 = AMAR
          ENDIF
         ENDIF
         IF(IFLON0.NE.0)THEN
          AMACRIT = (AMA0+AMA1)/2
          DEL = 1.D-8
          AMDEL = 2*DABS(AMA1-AMA0)/(AMA1+AMA0)
          IF(AMDEL.GT.DEL) GOTO 502
         ENDIF
        ENDIF
        DEL = 1.D-4
        XX(1) = AMACRIT - ALD - DEL
        XX(2) = AMACRIT - ALD
        XX(3) = AMACRIT + ALU
        XX(4) = AMACRIT + ALU + DEL
        DO J=1,4
         AMA = XX(J)
         CALL SUSYCP(TGBET)
         XX(J) = AMH
         IF(AMH.GE.2*AMA)THEN
          YY(J)=GF/16D0/DSQRT(2D0)/PI*AMZ**4/XX(J)
     .          *BETA(AMA**2/XX(J)**2)
         ELSEIF(AMH.LE.AMA)THEN
          YY(J) = 0
         ELSE
          XA=AMA**2/XX(J)**2
          XA1=(XA-1.D0)*(2.D0-.5D0*DLOG(XA))+(1.D0-5.D0*XA)
     .    *(DATAN((2.D0*XA-1.D0)/DSQRT(4.D0*XA-1.D0))
     .     -DATAN(1.D0/DSQRT(4.D0*XA-1.D0)))/DSQRT(4.D0*XA-1.D0)
          XA2=3*GF**2/16.D0/PI**3*AMZ**4/XX(J)*GAB**2*AMB**2
          YY(J)=XA1*XA2
         ENDIF
        ENDDO
        AMA = ZZMA
        CALL SUSYCP(TGBET)
        HAA = FINT(AMH,XX,YY)*GHAA**2
       ELSE
        HAA=GF/16D0/DSQRT(2D0)/PI*AMZ**4/AMH*BETA(AMA**2/AMH**2)
     .       *GHAA**2
       ENDIF
      ELSE
       IF (AMH.LE.2*AMA) THEN
        HAA=0
       ELSE
        HAA=GF/16D0/DSQRT(2D0)/PI*AMZ**4/AMH*BETA(AMA**2/AMH**2)
     .       *GHAA**2
       ENDIF
      ENDIF
C  H ---> A Z
      IF(IONSH.EQ.0)THEN
       DLD=1D0
       DLU=8D0
       XM1 = AMA+AMZ-DLD
       XM2 = AMA+AMZ+DLU
       IF (AMH.LT.AMA) THEN
        HAZ=0
       ELSEIF (AMH.LT.XM1) THEN
        IF(AMH.LE.DABS(AMZ-AMA))THEN
         HAZ=0
        ELSE
        HAZ=9.D0*GF**2/8.D0/PI**3*AMZ**4*AMH*GZAH**2*
     .      (7.D0/12.D0-10.D0/9.D0*SS+40.D0/27.D0*SS**2)
     .      *HVH((AMA/AMH)**2,(AMZ/AMH)**2)
        ENDIF
       ELSEIF (AMH.LT.XM2) THEN
        ZZMA=AMAR
165     AMA = AMAR - 1.D0
        CALL SUSYCP(TGBET)
        IF(AMH.LT.AMA+AMZ+DLU.AND.AMH.GT.AMA+AMZ-DLD) GOTO 165
        XX(1) = AMAR-1D0
        XX(2) = AMAR
        AMA = ZZMA
        CALL SUSYCP(TGBET)
166     AMA = AMAR + 1.D0
        CALL SUSYCP(TGBET)
        IF(AMH.LT.AMA+AMZ+DLU.AND.AMH.GT.AMA+AMZ-DLD) GOTO 166
        XX(3) = AMAR
        XX(4) = AMAR+1D0
        DO IJ=1,4
         AMA = XX(IJ)
         CALL SUSYCP(TGBET)
         XX(IJ) = AMH
         IF(AMH.LE.AMA+AMZ) THEN
          YY(IJ)=9.D0*GF**2/8.D0/PI**3*AMZ**4*XX(IJ)*
     .          (7.D0/12.D0-10.D0/9.D0*SS+40.D0/27.D0*SS**2)
     .          *HVH((AMA/XX(IJ))**2,(AMZ/XX(IJ))**2)
         ELSE
          CAZ=LAMB(AMA**2/XX(IJ)**2,AMZ**2/XX(IJ)**2)
     .       *LAMB(XX(IJ)**2/AMZ**2,AMA**2/AMZ**2)**2
          YY(IJ)=GF/8.D0/DSQRT(2D0)/PI*AMZ**4/XX(IJ)*CAZ
         ENDIF
        ENDDO
        AMA = ZZMA
        CALL SUSYCP(TGBET)
        HAZ = FINT(AMH,XX,YY)*GZAH**2
       ELSE
        CAZ=LAMB(AMA**2/AMH**2,AMZ**2/AMH**2)
     .     *LAMB(AMH**2/AMZ**2,AMA**2/AMZ**2)**2
        HAZ=GF/8.D0/DSQRT(2D0)/PI*AMZ**4/AMH*CAZ*GZAH**2
       ENDIF
      ELSE
       IF (AMH.LT.AMZ+AMA) THEN
        HAZ=0
       ELSE
        CAZ=LAMB(AMA**2/AMH**2,AMZ**2/AMH**2)
     .     *LAMB(AMH**2/AMZ**2,AMA**2/AMZ**2)**2
        HAZ=GF/8.D0/DSQRT(2D0)/PI*AMZ**4/AMH*CAZ*GZAH**2
       ENDIF
      ENDIF
C  H ---> H+ W+
      IF(IONSH.EQ.0)THEN
       DLD=3D0
       DLU=9D0
       XM1 = AMCH+AMW-DLD
       XM2 = AMCH+AMW+DLU
       IF (AMH.LT.AMCH) THEN
        HHW=0.D0
       ELSEIF (AMH.LT.XM1) THEN
        IF(AMH.LE.DABS(AMW-AMCH))THEN
         HHW=0
        ELSE
        HHW=9.D0*GF**2/8.D0/PI**3*AMW**4*AMH*GLVV**2*2
     .      *HVH((AMCH/AMH)**2,(AMW/AMH)**2)
        ENDIF
       ELSEIF (AMH.LT.XM2) THEN
        ZZMA=AMAR
167     AMA = AMAR - 1.D0
        CALL SUSYCP(TGBET)
        IF(AMH.LT.AMCH+AMW+DLU) GOTO 167
        XX(1) = AMAR-1D0
        XX(2) = AMAR
        AMA = ZZMA
        CALL SUSYCP(TGBET)
168     AMA = AMAR + 1.D0
        CALL SUSYCP(TGBET)
        IF(AMH.GT.AMCH+AMW-DLD) GOTO 168
        XX(3) = AMAR
        XX(4) = AMAR+1D0
        AMA = XX(1)
        CALL SUSYCP(TGBET)
        XX(1) = AMH
        CHW=LAMB(AMCH**2/XX(1)**2,AMW**2/XX(1)**2)
     .     *LAMB(XX(1)**2/AMW**2,AMCH**2/AMW**2)**2
        YY(1)=2*GF/8.D0/DSQRT(2D0)/PI*AMW**4/XX(1)*CHW
        AMA = XX(2)
        CALL SUSYCP(TGBET)
        XX(2) = AMH
        CHW=LAMB(AMCH**2/XX(2)**2,AMW**2/XX(2)**2)
     .     *LAMB(XX(2)**2/AMW**2,AMCH**2/AMW**2)**2
        YY(2)=2*GF/8.D0/DSQRT(2D0)/PI*AMW**4/XX(2)*CHW
        AMA = XX(3)
        CALL SUSYCP(TGBET)
        XX(3) = AMH
        YY(3)=9.D0*GF**2/8.D0/PI**3*AMW**4*XX(3)*2
     .       *HVH((AMCH/XX(3))**2,(AMW/XX(3))**2)
        AMA = XX(4)
        CALL SUSYCP(TGBET)
        XX(4) = AMH
        YY(4)=9.D0*GF**2/8.D0/PI**3*AMW**4*XX(4)*2
     .       *HVH((AMCH/XX(4))**2,(AMW/XX(4))**2)
        AMA = ZZMA
        CALL SUSYCP(TGBET)
        HHW=FINT(AMH,XX,YY)*GLVV**2
       ELSE
        CHW=LAMB(AMCH**2/AMH**2,AMW**2/AMH**2)
     .     *LAMB(AMH**2/AMW**2,AMCH**2/AMW**2)**2
        HHW=2*GF/8.D0/DSQRT(2D0)/PI*AMW**4/AMH*CHW*GLVV**2
       ENDIF
      ELSE
       IF (AMH.LT.AMW+AMCH) THEN
        HHW=0.D0
       ELSE
        CHW=LAMB(AMCH**2/AMH**2,AMW**2/AMH**2)
     .     *LAMB(AMH**2/AMW**2,AMCH**2/AMW**2)**2
        HHW=2*GF/8.D0/DSQRT(2D0)/PI*AMW**4/AMH*CHW*GLVV**2
       ENDIF
      ENDIF

C ========================== SUSY DECAYS 
C
      IF(IOFSUSY.EQ.0) THEN
C  HH ----> CHARGINOS
      DO 741 I=1,2
      DO 741 J=1,2
      IF (AMH.GT.AMCHAR(I)+AMCHAR(J)) THEN
      WHHCH(I,J)=GF*AMW**2/(2*PI*DSQRT(2.D0))/AMH 
     .     *LAMB(AMCHAR(I)**2/AMH**2,AMCHAR(J)**2/AMH**2)
     .     *( (AC1(I,J)**2+AC1(J,I)**2)*(AMH**2-AMCHAR(I)
     .         **2-AMCHAR(J)**2)-4.D0*AC1(I,J)*AC1(J,I)* 
     .         AMCHAR(I)*AMCHAR(J) ) 
      ELSE
      WHHCH(I,J)=0.D0
      ENDIF
 741  CONTINUE
      WHHCHT=WHHCH(1,1)+WHHCH(1,2)+WHHCH(2,1)+WHHCH(2,2)
C
C  HH ----> NEUTRALINOS 
      DO 742 I=1,4
      DO 742 J=1,4
      IF (AMH.GT.AMNEUT(I)+AMNEUT(J)) THEN
      WHHNE(I,J)=GF*AMW**2/(2*PI*DSQRT(2.D0))/AMH
     .         *AN1(I,J)**2*(AMH**2-(XMNEUT(I)+XMNEUT(J))**2)
     .         *LAMB(AMNEUT(I)**2/AMH**2,AMNEUT(J)**2/AMH**2)
      ELSE 
      WHHNE(I,J)=0.D0
      ENDIF
 742  CONTINUE
      WHHNET= WHHNE(1,1)+WHHNE(1,2)+WHHNE(1,3)+WHHNE(1,4)
     .       +WHHNE(2,1)+WHHNE(2,2)+WHHNE(2,3)+WHHNE(2,4)
     .       +WHHNE(3,1)+WHHNE(3,2)+WHHNE(3,3)+WHHNE(3,4)
     .       +WHHNE(4,1)+WHHNE(4,2)+WHHNE(4,3)+WHHNE(4,4)
C
C  HH ----> SLEPTONS 
C
      IF (AMH.GT.2.D0*AMSE(1)) THEN
      WHHSLEL=2*GF/2.D0/DSQRT(2D0)/PI*AMZ**4/AMH*DCOS(B+A)**2
     .      *BETA(AMSE(1)**2/AMH**2)*(-0.5D0+SS)**2
      ELSE
      WHHSLEL=0.D0
      ENDIF

      IF (AMH.GT.2.D0*AMSE(2)) THEN
      WHHSLER=2*GF/2.D0/DSQRT(2D0)/PI*AMZ**4/AMH*DCOS(B+A)**2
     .      *BETA(AMSE(2)**2/AMH**2)*SS**2
      ELSE
      WHHSLER=0.D0
      ENDIF

      IF (AMH.GT.2.D0*AMSN(1)) THEN
      WHHSLNL=3*GF/2.D0/DSQRT(2D0)/PI*AMZ**4/AMH*DCOS(B+A)**2
     .      *BETA(AMSN(1)**2/AMH**2)*0.5D0**2
      ELSE
      WHHSLNL=0.D0
      ENDIF

      DO 748 I=1,2
      DO 748 J=1,2
      IF(AMH.GT.AMSL(I)+AMSL(J)) THEN
      WHHSTAU(I,J)=GF*AMZ**4/2.D0/DSQRT(2.D0)/PI*GHEE(I,J)**2*
     .      LAMB(AMSL(I)**2/AMH**2,AMSL(J)**2/AMH**2)/AMH
      ELSE
      WHHSTAU(I,J)=0.D0
      ENDIF
 748  CONTINUE

      WHHSLT=WHHSTAU(1,1)+WHHSTAU(1,2)+WHHSTAU(2,1)+WHHSTAU(2,2) 
     .       +WHHSLEL+WHHSLER+WHHSLNL
C
C  HH ----> SQUARKS 
C
      IF (AMH.GT.2.D0*AMSU(1)) THEN
      WHHSQUL=6*GF/2.D0/DSQRT(2D0)/PI*AMZ**4/AMH*DCOS(B+A)**2
     .      *BETA(AMSU(1)**2/AMH**2)*(0.5D0-2.D0/3.D0*SS)**2
      ELSE
      WHHSQUL=0.D0
      ENDIF

      IF (AMH.GT.2.D0*AMSU(2)) THEN
      WHHSQUR=6*GF/2.D0/DSQRT(2D0)/PI*AMZ**4/AMH*DCOS(B+A)**2
     .      *BETA(AMSU(2)**2/AMH**2)*(-2.D0/3.D0*SS)**2
      ELSE
      WHHSQUR=0.D0
      ENDIF

      IF (AMH.GT.2.D0*AMSD(1)) THEN
      WHHSQDL=6*GF/2.D0/DSQRT(2D0)/PI*AMZ**4/AMH*DCOS(B+A)**2
     .      *BETA(AMSD(1)**2/AMH**2)*(-0.5D0+1.D0/3.D0*SS)**2
      ELSE
      WHHSQDL=0.D0
      ENDIF

      IF (AMH.GT.2.D0*AMSD(2)) THEN
      WHHSQDR=6*GF/2.D0/DSQRT(2D0)/PI*AMZ**4/AMH*DCOS(B+A)**2
     .      *BETA(AMSD(2)**2/AMH**2)*(+1.D0/3.D0*SS)**2
      ELSE
      WHHSQDR=0.D0
      ENDIF

      WHHSQ=WHHSQUL+WHHSQUR+WHHSQDL+WHHSQDR
C
C  HH ----> STOPS 
      DO 743 I=1,2
      DO 743 J=1,2
      IF(AMH.GT.AMST(I)+AMST(J)) THEN
      WHHST(I,J)=3*GF*AMZ**4/2.D0/DSQRT(2.D0)/PI*GHTT(I,J)**2*
     .      LAMB(AMST(I)**2/AMH**2,AMST(J)**2/AMH**2)/AMH
      ELSE
      WHHST(I,J)=0.D0
      ENDIF
 743  CONTINUE
C
C  HH ----> SBOTTOMS 
      DO 744 I=1,2
      DO 744 J=1,2
      IF(AMH.GT.AMSB(I)+AMSB(J)) THEN
      WHHSB(I,J)=3*GF*AMZ**4/2.D0/DSQRT(2.D0)/PI*GHBB(I,J)**2*
     .      LAMB(AMSB(I)**2/AMH**2,AMSB(J)**2/AMH**2)/AMH
      ELSE
      WHHSB(I,J)=0.D0
      ENDIF
 744  CONTINUE
C
      WHHSTT=WHHST(1,1)+WHHST(1,2)+WHHST(2,1)+WHHST(2,2) 
      WHHSBB=WHHSB(1,1)+WHHSB(1,2)+WHHSB(2,1)+WHHSB(2,2) 
      WHHSQT=WHHSTT+WHHSBB+WHHSQ
C
      ELSE 
      WHHCHT=0.D0
      WHHNET=0.D0
      WHHSLT=0.D0
      WHHSQT=0.D0
C--Change thanks to Elzbieta Richter-Was
      DO I=1,2
       DO J=1,2
        WHHCH(I,J)=0.D0
        WHHST(I,J)=0.D0
        WHHSB(I,J)=0.D0
        WHHSTAU(I,J)=0.D0
       ENDDO
      ENDDO
      DO I=1,4
       DO J=1,4
        WHHNE(I,J)=0.D0
       ENDDO
      ENDDO
      ENDIF
C
C    ==========  TOTAL WIDTH AND BRANCHING RATIOS 
      WTOT=HLL+HMM+HSS+HCC+HBB+HTT+HGG+HGA+HZGA+HWW+HZZ+HHH+HAA+HAZ
     .    +HHW+WHHCHT+WHHNET+WHHSLT+WHHSQT
      HHBRT=HTT/WTOT
      HHBRB=HBB/WTOT
      HHBRL=HLL/WTOT
      HHBRM=HMM/WTOT
      HHBRS=HSS/WTOT
      HHBRC=HCC/WTOT
      HHBRG=HGG/WTOT
      HHBRGA=HGA/WTOT
      HHBRZGA=HZGA/WTOT
      HHBRW=HWW/WTOT
      HHBRZ=HZZ/WTOT
      HHBRH=HHH/WTOT
      HHBRA=HAA/WTOT
      HHBRAZ=HAZ/WTOT
      HHBRHW=HHW/WTOT
      DO 841 I=1,2
      DO 841 J=1,2
      HHBRSC(I,J)=WHHCH(I,J)/WTOT
841   CONTINUE
      DO 842 I=1,4
      DO 842 J=1,4
      HHBRSN(I,J)=WHHNE(I,J)/WTOT
842   CONTINUE
      HHBRCHT=WHHCHT/WTOT 
      HHBRNET=WHHNET/WTOT 
      HHBRSL=WHHSLT/WTOT
      HHBRSQ=WHHSQ/WTOT
      HHBRSQT=WHHSQT/WTOT
      DO 843 I=1,2
      DO 843 J=1,2
      HHBRST(I,J)=WHHST(I,J)/WTOT
843   CONTINUE
      DO 844 I=1,2
      DO 844 J=1,2
      HHBRSB(I,J)=WHHSB(I,J)/WTOT
844   CONTINUE
      HHWDTH=WTOT

      BHHSLNL = WHHSLNL/WTOT
      BHHSLEL = WHHSLEL/WTOT
      BHHSLER = WHHSLER/WTOT
      BHHSQUL = WHHSQUL/WTOT
      BHHSQUR = WHHSQUR/WTOT
      BHHSQDL = WHHSQDL/WTOT
      BHHSQDR = WHHSQDR/WTOT
      DO I = 1,2
       DO J = 1,2
        BHHST(I,J) = WHHST(I,J)/WTOT
        BHHSB(I,J) = WHHSB(I,J)/WTOT
        BHHSTAU(I,J) = WHHSTAU( I,J)/WTOT
       ENDDO
      ENDDO

      ENDIF

      IF(IHIGGS.EQ.3.OR.IHIGGS.EQ.5)THEN 
C
C        =========================================================
C                       CP ODD  HIGGS DECAYS
C        =========================================================
C     =============  RUNNING MASSES 
      RMS = RUNM(AMA,3)
      RMC = RUNM(AMA,4)
      RMB = RUNM(AMA,5)
      RMT = RUNM(AMA,6)
      RATCOUP = GAT/GAB
      HIGTOP = AMA**2/AMT**2

      ASH=ALPHAS(AMA,2)
      AMC0=1.D8
      AMB0=2.D8
C     AMT0=3.D8
      AS3=ALPHAS(AMA,2)
      AMC0=AMC
      AS4=ALPHAS(AMA,2)
      AMB0=AMB
C     AMT0=AMT

C     =============== PARTIAL WIDTHS
C  A ---> G G
       EPS=1.D-8
       NFEXT = 3
       ASG = AS3
       CTT = 4*AMT**2/AMA**2*DCMPLX(1D0,-EPS)
       CTB = 4*AMB**2/AMA**2*DCMPLX(1D0,-EPS)
       CAT = CTT*CF(CTT)*GAT
       CAB = CTB*CF(CTB)*GAB
       CTC = 4*AMC**2/AMA**2*DCMPLX(1D0,-EPS)
       CAC = CTC*CF(CTC)*GAT
       FQCD=AGGQCD(ASG,NFEXT)
       XFAC = CDABS(CAT+CAB+CAC)**2*FQCD
       HGG=GF/(16.D0*PI*DSQRT(2.D0))*AMA**3*(ASG/PI)**2*XFAC

C  A ---> G G* ---> G CC   TO BE ADDED TO A ---> CC
       NFEXT = 4
       ASG = AS4
       FQCD=AGGQCD(ASG,NFEXT)
       XFAC = CDABS(CAT+CAB+CAC)**2*FQCD
       DCC=GF/(16.D0*PI*DSQRT(2.D0))*AMA**3*(ASG/PI)**2*XFAC
     .     - HGG

C  A ---> G G* ---> G BB   TO BE ADDED TO A ---> BB
       NFEXT = 5
       ASG = ASH
       FQCD=AGGQCD(ASG,NFEXT)
       XFAC = CDABS(CAT+CAB+CAC)**2*FQCD
       DBB=GF/(16.D0*PI*DSQRT(2.D0))*AMA**3*(ASG/PI)**2*XFAC
     .     - HGG - DCC

      IF(NFGG.EQ.5)THEN
       HGG = HGG + DBB + DCC
       DBB = 0
       DCC = 0
      ELSEIF(NFGG.EQ.4)THEN
       HGG = HGG + DCC
       DCC = 0
      ENDIF

C  A ---> MU MU
      IF(AMA.LE.2*AMMUON) THEN
       HMM = 0
      ELSE
      HMM=AFF(AMA,(AMMUON/AMA)**2)*GAB**2
      ENDIF
C  A ---> LL
      IF(AMA.LE.2*AMTAU) THEN
       HLL = 0
      ELSE
      HLL=AFF(AMA,(AMTAU/AMA)**2)*GAB**2
      ENDIF
C  A --> SS
      IF(AMA.LE.2*AMS) THEN
       HSS = 0
      ELSE
       HS1=3.D0*HFF(AMA,(AMS/AMA)**2)
     .    *GAB**2
     .    *TQCDA(AMS**2/AMA**2)
       HS2=3.D0*AFF(AMA,(RMS/AMA)**2)
     .    *GAB**2
     .    *QCDA(RMS**2/AMA**2)
       IF(HS2.LT.0.D0) HS2 = 0
       RAT = 2*AMS/AMA
       HSS = QQINT(RAT,HS1,HS2)
      ENDIF
C  A --> CC
      RATCOUP = 1
      IF(AMA.LE.2*AMC) THEN
       HCC = 0
      ELSE
       HC1=3.D0*HFF(AMA,(AMC/AMA)**2)
     .    *GAT**2
     .    *TQCDA(AMC**2/AMA**2)
       HC2=3.D0*HFF(AMA,(RMC/AMA)**2)
     .    *GAT**2
     .    *QCDA(RMC**2/AMA**2)
     .   + DCC
       IF(HC2.LT.0.D0) HC2 = 0
       RAT = 2*AMC/AMA
       HCC = QQINT(RAT,HC1,HC2)
      ENDIF
C  A --> BB :
      RATCOUP = GAT/GAB
      IF(AMA.LE.2*AMB) THEN
       HBB = 0
      ELSE
       HB1=3.D0*HFF(AMA,(AMB/AMA)**2)
     .    *GAB**2
     .    *TQCDA(AMB**2/AMA**2)
       HB2=3.D0*HFF(AMA,(RMB/AMA)**2)
     .    *GAB**2
     .    *QCDA(RMB**2/AMA**2)
     .   + DBB
       IF(HB2.LT.0.D0) HB2 = 0
       RAT = 2*AMB/AMA
       HBB = QQINT(RAT,HB1,HB2)
      ENDIF
C  A --> TT :
      RATCOUP = 0
      IF(IONSH.EQ.0)THEN
       DLD=3D0
       DLU=4D0
       XM1 = 2D0*AMT-DLD
       XM2 = 2D0*AMT+DLU
       IF (AMA.LE.AMT+AMW+AMB) THEN
        HTT=0.D0
       ELSEIF (AMA.LE.XM1) THEN
        FACTT=6.D0*GF**2*AMA**3*AMT**2/2.D0/128.D0/PI**3*GAT**2
        CALL ATOTT(AMA,AMT,AMB,AMW,AMCH,ATT0)
        HTT=FACTT*ATT0
       ELSEIF (AMA.LE.XM2) THEN
        XX(1) = XM1-1D0
        XX(2) = XM1
        XX(3) = XM2
        XX(4) = XM2+1D0
        FACTT=6.D0*GF**2*XX(1)**3*AMT**2/2.D0/128.D0/PI**3
        CALL ATOTT(XX(1),AMT,AMB,AMW,AMCH,ATT0)
        YY(1)=FACTT*ATT0
        FACTT=6.D0*GF**2*XX(2)**3*AMT**2/2.D0/128.D0/PI**3
        CALL ATOTT(XX(2),AMT,AMB,AMW,AMCH,ATT0)
        YY(2)=FACTT*ATT0
        XMT = RUNM(XX(3),6)
        XYZ1 =3.D0*AFF(XX(3),(AMT/XX(3))**2)
     .    *TQCDA(AMT**2/XX(3)**2)
        XYZ2 =3.D0*AFF(XX(3),(XMT/XX(3))**2)
     .    *QCDA(XMT**2/XX(3)**2)
        IF(XYZ2.LT.0.D0) XYZ2 = 0
        RAT = 2*AMT/XX(3)
        YY(3) = QQINT(RAT,XYZ1,XYZ2)
        XMT = RUNM(XX(4),6)
        XYZ1 =3.D0*AFF(XX(4),(AMT/XX(4))**2)
     .    *TQCDA(AMT**2/XX(4)**2)
        XYZ2 =3.D0*AFF(XX(4),(XMT/XX(4))**2)
     .    *QCDA(XMT**2/XX(4)**2)
        IF(XYZ2.LT.0.D0) XYZ2 = 0
        RAT = 2*AMT/XX(4)
        YY(4) = QQINT(RAT,XYZ1,XYZ2)
        HTT = FINT(AMA,XX,YY)*GAT**2
       ELSE
        HT1=3.D0*AFF(AMA,(AMT/AMA)**2)*GAT**2
     .    *TQCDA(AMT**2/AMA**2)
        HT2=3.D0*AFF(AMA,(RMT/AMA)**2)*GAT**2
     .    *QCDA(RMT**2/AMA**2)
        IF(HT2.LT.0.D0) HT2 = 0
        RAT = 2*AMT/AMA
        HTT = QQINT(RAT,HT1,HT2)
       ENDIF
      ELSE
       IF (AMA.LE.2.D0*AMT) THEN
        HTT=0.D0
       ELSE
        HT1=3.D0*AFF(AMA,(AMT/AMA)**2)*GAT**2
     .    *TQCDA(AMT**2/AMA**2)
        HT2=3.D0*AFF(AMA,(RMT/AMA)**2)*GAT**2
     .    *QCDA(RMT**2/AMA**2)
        IF(HT2.LT.0.D0) HT2 = 0
        RAT = 2*AMT/AMA
        HTT = QQINT(RAT,HT1,HT2)
       ENDIF
      ENDIF
C  A ---> GAMMA GAMMA
       EPS=1.D-8
       CTT = 4*AMT**2/AMA**2*DCMPLX(1D0,-EPS)
       CTB = 4*AMB**2/AMA**2*DCMPLX(1D0,-EPS)
       CAT = 4/3D0 * CTT*CF(CTT)*GAT
       CAB = 1/3D0 * CTB*CF(CTB)*GAB
       CTC = 4*AMC**2/AMA**2*DCMPLX(1D0,-EPS)
       CAC = 4/3D0 * CTC*CF(CTC)*GAT
       CTL = 4*AMTAU**2/AMA**2*DCMPLX(1D0,-EPS)
       CAL = 1.D0  * CTL*CF(CTL)*GAB
       IF(IOFSUSY.EQ.0) THEN 
       CX1 = 4*AMCHAR(1)**2/AMA**2*DCMPLX(1D0,-EPS)
       CX2 = 4*AMCHAR(2)**2/AMA**2*DCMPLX(1D0,-EPS)
       CAX1= AMW/AMCHAR(1) * CX1*CF(CX1) * 2*AC3(1,1) 
       CAX2= AMW/AMCHAR(2) * CX2*CF(CX2) * 2*AC3(2,2) 
       XFAC = CDABS(CAT+CAB+CAC+CAL+CAX1+CAX2)**2
       ELSE 
       XFAC = CDABS(CAT+CAB+CAC+CAL)**2
       ENDIF
       HGA=GF/(32.D0*PI*DSQRT(2.D0))*AMA**3*(ALPH/PI)**2*XFAC
C  A ---> Z GAMMA
      IF(AMA.LE.AMZ)THEN
       HZGA=0
      ELSE
       TS = SS/CS
       FT = -3*2D0/3*(1-4*2D0/3*SS)/DSQRT(SS*CS)*GAT
       FB = 3*1D0/3*(-1+4*1D0/3*SS)/DSQRT(SS*CS)*GAB
       EPS=1.D-8
       CTT = 4*AMT**2/AMA**2*DCMPLX(1D0,-EPS)
       CTB = 4*AMB**2/AMA**2*DCMPLX(1D0,-EPS)
       CLT = 4*AMT**2/AMZ**2*DCMPLX(1D0,-EPS)
       CLB = 4*AMB**2/AMZ**2*DCMPLX(1D0,-EPS)
       CAT = FT*(- CI2(CTT,CLT))
       CAB = FB*(- CI2(CTB,CLB))
       XFAC = CDABS(CAT+CAB)**2
       ACOUP = DSQRT(2D0)*GF*AMZ**2*SS*CS/PI**2
       HZGA = GF/(4.D0*PI*DSQRT(2.D0))*AMA**3*(ALPH/PI)*ACOUP/16.D0
     .        *XFAC*(1-AMZ**2/AMA**2)**3
      ENDIF
C  A ---> H Z* ---> HFF
      IF(IONSH.EQ.0)THEN
       DLD=3D0
       DLU=5D0
       XM1 = AML+AMZ-DLD
       XM2 = AML+AMZ+DLU
       IF (AMA.LE.AML) THEN
        HAZ=0
       ELSEIF (AMA.LE.XM1) THEN
        IF (AMA.LE.DABS(AMZ-AML)) THEN
         HAZ=0
        ELSE
         HAZ=9.D0*GF**2/8.D0/PI**3*AMZ**4*AMA*GZAL**2*
     .      (7.D0/12.D0-10.D0/9.D0*SS+40.D0/27.D0*SS**2)
     .      *HVH((AML/AMA)**2,(AMZ/AMA)**2)
        ENDIF
       ELSEIF (AMA.LE.XM2) THEN
        XX(1) = XM1-1D0
        XX(2) = XM1
        XX(3) = XM2
        XX(4) = XM2+1D0
        YY(1)=9.D0*GF**2/8.D0/PI**3*AMZ**4*XX(1)*
     .      (7.D0/12.D0-10.D0/9.D0*SS+40.D0/27.D0*SS**2)
     .      *HVH((AML/XX(1))**2,(AMZ/XX(1))**2)
        YY(2)=9.D0*GF**2/8.D0/PI**3*AMZ**4*XX(2)*
     .      (7.D0/12.D0-10.D0/9.D0*SS+40.D0/27.D0*SS**2)
     .      *HVH((AML/XX(2))**2,(AMZ/XX(2))**2)
        CAZ=LAMB(AML**2/XX(3)**2,AMZ**2/XX(3)**2)
     .     *LAMB(XX(3)**2/AMZ**2,AML**2/AMZ**2)**2
        YY(3)=GF/8D0/DSQRT(2D0)/PI*AMZ**4/XX(3)*CAZ
        CAZ=LAMB(AML**2/XX(4)**2,AMZ**2/XX(4)**2)
     .     *LAMB(XX(4)**2/AMZ**2,AML**2/AMZ**2)**2
        YY(4)=GF/8D0/DSQRT(2D0)/PI*AMZ**4/XX(4)*CAZ
        HAZ = FINT(AMA,XX,YY)*GZAL**2
       ELSE
        CAZ=LAMB(AML**2/AMA**2,AMZ**2/AMA**2)
     .     *LAMB(AMA**2/AMZ**2,AML**2/AMZ**2)**2
        HAZ=GF/8D0/DSQRT(2D0)/PI*AMZ**4/AMA*GZAL**2*CAZ
       ENDIF
      ELSE
       IF (AMA.LE.AMZ+AML) THEN
        HAZ=0
       ELSE
        CAZ=LAMB(AML**2/AMA**2,AMZ**2/AMA**2)
     .     *LAMB(AMA**2/AMZ**2,AML**2/AMZ**2)**2
        HAZ=GF/8D0/DSQRT(2D0)/PI*AMZ**4/AMA*GZAL**2*CAZ
       ENDIF
      ENDIF
C
C ========================== SUSY DECAYS  
C
      IF(IOFSUSY.EQ.0) THEN 
C  A ----> CHARGINOS
      DO 731 I=1,2
      DO 731 J=1,2
      IF (AMA.GT.AMCHAR(I)+AMCHAR(J)) THEN
      WHACH(I,J)=GF*AMW**2/(2*PI*DSQRT(2.D0))/AMA
     .     *LAMB(AMCHAR(I)**2/AMA**2,AMCHAR(J)**2/AMA**2)
     .     *( (AC3(I,J)**2+AC3(J,I)**2)*(AMA**2-AMCHAR(I)
     .         **2-AMCHAR(J)**2)+4.D0*AC3(I,J)*AC3(J,I)* 
     .         AMCHAR(I)*AMCHAR(J) ) 
      ELSE 
      WHACH(I,J)=0.D0
      ENDIF
 731  CONTINUE
      WHACHT=WHACH(1,1)+WHACH(1,2)+WHACH(2,1)+WHACH(2,2)
C  A ----> NEUTRALINOS 
      DO 732 I=1,4
      DO 732 J=1,4
      IF (AMA.GT.AMNEUT(I)+AMNEUT(J)) THEN
      WHANE(I,J)=GF*AMW**2/(2*PI*DSQRT(2.D0))/AMA
     .         *AN3(I,J)**2*(AMA**2-(XMNEUT(I)-XMNEUT(J))**2)
     .         *LAMB(AMNEUT(I)**2/AMA**2,AMNEUT(J)**2/AMA**2)
      ELSE 
      WHANE(I,J)=0.D0
      ENDIF
 732  CONTINUE
      WHANET= WHANE(1,1)+WHANE(1,2)+WHANE(1,3)+WHANE(1,4)
     .       +WHANE(2,1)+WHANE(2,2)+WHANE(2,3)+WHANE(2,4)
     .       +WHANE(3,1)+WHANE(3,2)+WHANE(3,3)+WHANE(3,4)
     .       +WHANE(4,1)+WHANE(4,2)+WHANE(4,3)+WHANE(4,4)

C  A ----> STAU'S 
C
      IF(AMA.GT.AMSL(1)+AMSL(2)) THEN
      WHASL=GF*AMZ**4/DSQRT(2.D0)/PI*GAEE**2*
     .      LAMB(AMSL(1)**2/AMA**2,AMSL(2)**2/AMA**2)/AMA
      ELSE
      WHASL=0.D0
      ENDIF
C
C  A ----> STOPS 
C
      IF(AMA.GT.AMST(1)+AMST(2)) THEN
      WHAST=3*GF*AMZ**4/DSQRT(2.D0)/PI*GATT**2*
     .      LAMB(AMST(1)**2/AMA**2,AMST(2)**2/AMA**2)/AMA
      ELSE
      WHAST=0.D0
      ENDIF
C
C  A ----> SBOTTOMS 
C
      IF(AMA.GT.AMSB(1)+AMSB(2)) THEN
      WHASB=3*GF*AMZ**4/DSQRT(2.D0)/PI*GABB**2*
     .      LAMB(AMSB(1)**2/AMA**2,AMSB(2)**2/AMA**2)/AMA
      ELSE
      WHASB=0.D0
      ENDIF
C
      ELSE 
      WHACHT=0.D0
      WHANET=0.D0
      WHASL=0.D0
      WHAST=0.D0
      WHASB=0.D0
C--Change thanks to Elzbieta Richter-Was
      DO I=1,2
       DO J=1,2
        WHACH(I,J)=0.D0
       ENDDO
      ENDDO
      DO I=1,4
       DO J=1,4
        WHANE(I,J)=0.D0
       ENDDO
      ENDDO
      ENDIF
C
C    ==========  TOTAL WIDTH AND BRANCHING RATIOS 
      WTOT=HLL+HMM+HSS+HCC+HBB+HGG+HGA+HZGA+HAZ+HTT
     .    +WHACHT+WHANET+WHASL+WHAST+WHASB 
      ABRT=HTT/WTOT
      ABRB=HBB/WTOT
      ABRL=HLL/WTOT
      ABRM=HMM/WTOT
      ABRS=HSS/WTOT
      ABRC=HCC/WTOT
      ABRG=HGG/WTOT
      ABRGA=HGA/WTOT
      ABRZGA=HZGA/WTOT
      ABRZ=HAZ/WTOT
      DO 831 I=1,2
      DO 831 J=1,2
      HABRSC(I,J)=WHACH(I,J)/WTOT
831   CONTINUE
      DO 832 I=1,4
      DO 832 J=1,4
      HABRSN(I,J)=WHANE(I,J)/WTOT
832   CONTINUE
      HABRCHT=WHACHT/WTOT      
      HABRNET=WHANET/WTOT      
      HABRSL=WHASL/WTOT 
      HABRST=WHAST/WTOT 
      HABRSB=WHASB/WTOT 
C 
      AWDTH=WTOT

      BHASTAU = WHASL/WTOT
      BHASB = WHASB/WTOT
      BHAST = WHAST/WTOT

C    ==============================================================
      ENDIF

      RETURN
      END
 
      DOUBLE PRECISION FUNCTION BIJ(X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION LAMB
      BIJ = (1-X-Y)/LAMB(X,Y)*(
     .          4*SP(XI(X,Y)*XI(Y,X))
     .        - 2*SP(-XI(X,Y)) - 2*SP(-XI(Y,X))
     .        + 2*DLOG(XI(X,Y)*XI(Y,X))*DLOG(1-XI(X,Y)*XI(Y,X))
     .        - DLOG(XI(X,Y))*DLOG(1+XI(X,Y))
     .        - DLOG(XI(Y,X))*DLOG(1+XI(Y,X))
     .          )
     .        -4*(DLOG(1-XI(X,Y)*XI(Y,X))
     .      +XI(X,Y)*XI(Y,X)/(1-XI(X,Y)*XI(Y,X))*DLOG(XI(X,Y)*XI(Y,X)))
     .        + (LAMB(X,Y)+X-Y)/LAMB(X,Y)*(DLOG(1+XI(X,Y))
     .                 - XI(X,Y)/(1+XI(X,Y))*DLOG(XI(X,Y)))
     .        + (LAMB(X,Y)-X+Y)/LAMB(X,Y)*(DLOG(1+XI(Y,X))
     .                 - XI(Y,X)/(1+XI(Y,X))*DLOG(XI(Y,X)))
      RETURN
      END

      DOUBLE PRECISION FUNCTION BETA(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      BETA=DSQRT(1.D0-4.D0*X)
      RETURN
      END

      DOUBLE PRECISION FUNCTION LAMB(X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LAMB=DSQRT((1.D0-X-Y)**2-4.D0*X*Y)
      RETURN
      END

      DOUBLE PRECISION FUNCTION XI(X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION LAMB
      XI = 2*X/(1-X-Y+LAMB(X,Y))
      RETURN
      END

C *****************************************************************
C ************* SUBROUTINE FOR THE SUSY COUPLINGS *****************
C *****************************************************************
      SUBROUTINE SUSYCP(TGBET)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      DOUBLE PRECISION LA1,LA2,LA3,LA4,LA5,LA6,LA7,LA3T
      COMPLEX*16 F0
      DIMENSION MST(2),GLTT(2,2),GHTT(2,2),
     .          MSB(2),GLBB(2,2),GHBB(2,2)
      COMMON/FLAG/IHIGGS,NNLO,IPOLE
      COMMON/MASSES/AMS,AMC,AMB,AMT
      COMMON/HMASS/AMSM,AMA,AML,AMH,AMCH,AMAR
      COMMON/CHIMASS/AMCHI
      COMMON/HSELF/LA1,LA2,LA3,LA4,LA5,LA6,LA7
      COMMON/BREAK/AMEL,AMER,AMSQ,AMUR,AMDR,AL,AU,AD,AMU,AM2
      COMMON/PARAM/GF,ALPH,AMTAU,AMMUON,AMZ,AMW
      COMMON/COUP/GAT,GAB,GLT,GLB,GHT,GHB,GZAH,GZAL,
     .            GHHH,GLLL,GHLL,GLHH,GHAA,GLAA,GLVV,GHVV,
     .            GLPM,GHPM,B,A
      PI=4*DATAN(1D0)
      V=1.D0/DSQRT(DSQRT(2.D0)*GF)
      BET=DATAN(TGBET)
      SB = DSIN(BET)
      CB = DCOS(BET)
      AMAR = AMA
C  ============ HEAVIEST CHARGINO MASS NEEDED FOR SUBH ========== 
      AMCHI2=AM2**2+AMU**2+2.D0*AMW**2+DSQRT((AM2**2-AMU**2)**2
     .      +4.D0*AMW**4*DCOS(2.D0*BET)**2+4.D0*AMW**2*
     .      (AM2**2+AMU**2+2.D0*AMU*AM2*DSIN(2.D0*BET) ) ) 
      AMCHI=DSQRT(0.5D0*AMCHI2)
C ===============================================================
C ========== RUNNING MASSES
***
      CALL SUBH(AMA,TGBET,AMSQ,AMUR,AMDR,AMT,AU,AD,AMU,AMCHI,
     .          AMLR,AMHR,AMCH,SA,CA,TANBA)
      write(6,*) 'SUBH called OK'
      call pole(3,amchi,ama,tgbet,amsq,amur,amdr,amt,au,ad,amu,
     &          dum1,amlr,dum2,amhr,duma,sa,ca,dumsa,dumca,
     &          dumst1,dumst2,dumsb1,dumsb2,dumtba,500d0)
C      call pole(3,amchi,ama,tgbet,amsq,amur,amdr,amt,au,ad,amu,
C     &          amlr,amhr,duma,sa,ca,dumtba,500d0)
      write(6,*) 'test: ',ama,duma,tgbet,dumtba
****
      LA3T=LA3+LA4+LA5
      AMA2=AMAR**2
      AML2=AMLR**2
      AMH2=AMHR**2
      AMP2=AMCH**2
C ========== HIGGS COUPLINGS 
      SBMA = SB*CA-CB*SA
      CBMA = CB*CA+SB*SA
      SBPA = SB*CA+CB*SA
      CBPA = CB*CA-SB*SA
      S2A = 2*SA*CA
      C2A = CA**2-SA**2
      S2B = 2*SB*CB
      C2B = CB**2-SB**2
      GLZZ = 1/V/2*AML2*SBMA
      GHZZ = 1/V/2*AMH2*CBMA
      GLWW = 2*GLZZ
      GHWW = 2*GHZZ
      GLAZ = 1/V*(AML2-AMA2)*CBMA
      GHAZ = -1/V*(AMH2-AMA2)*SBMA
      GLPW = -1/V*(AMP2-AML2)*CBMA
      GLMW = GLPW
      GHPW = 1/V*(AMP2-AMH2)*SBMA
      GHMW = GHPW
      GAPW = 1/V*(AMP2-AMA2)
      GAMW = -GAPW
      GHHH = V/2*(LA1*CA**3*CB + LA2*SA**3*SB + LA3T*SA*CA*SBPA
     .     + LA6*CA**2*(3*SA*CB+CA*SB) + LA7*SA**2*(3*CA*SB+SA*CB))
      GLLL = -V/2*(LA1*SA**3*CB - LA2*CA**3*SB + LA3T*SA*CA*CBPA
     .     - LA6*SA**2*(3*CA*CB-SA*SB) + LA7*CA**2*(3*SA*SB-CA*CB))
      GLHH = -3*V/2*(LA1*CA**2*CB*SA - LA2*SA**2*SB*CA
     .     + LA3T*(SA**3*CB-CA**3*SB+2*SBMA/3)
     .     - LA6*CA*(CB*C2A-SA*SBPA) - LA7*SA*(C2A*SB+CA*SBPA))
      GHLL = 3*V/2*(LA1*SA**2*CB*CA + LA2*CA**2*SB*SA
     .     + LA3T*(SA**3*SB+CA**3*CB-2*CBMA/3)
     .     - LA6*SA*(CB*C2A+CA*CBPA) + LA7*CA*(C2A*SB+SA*CBPA))
      GLAA = -V/2*(LA1*SB**2*CB*SA - LA2*CB**2*SB*CA
     .     - LA3T*(SB**3*CA-CB**3*SA) + 2*LA5*SBMA
     .     - LA6*SB*(CB*SBPA+SA*C2B) - LA7*CB*(C2B*CA-SB*SBPA))
      GHAA = V/2*(LA1*SB**2*CB*CA + LA2*CB**2*SB*SA
     .     + LA3T*(SB**3*SA+CB**3*CA) - 2*LA5*CBMA
     .     - LA6*SB*(CB*CBPA+CA*C2B) + LA7*CB*(SB*CBPA+SA*C2B))
      GLPM = 2*GLAA + V*(LA5 - LA4)*SBMA
      GHPM = 2*GHAA + V*(LA5 - LA4)*CBMA
      GLZZ = 2*GLZZ
      GHZZ = 2*GHZZ
      GLLL = 6*GLLL
      GHHH = 6*GHHH
      GLHH = 2*GLHH
      GHLL = 2*GHLL
      GLAA = 2*GLAA
      GHAA = 2*GHAA
      XNORM = AMZ**2/V
      GLLL = GLLL/XNORM
      GHLL = GHLL/XNORM
      GLHH = GLHH/XNORM
      GHHH = GHHH/XNORM
      GHAA = GHAA/XNORM
      GLAA = GLAA/XNORM
      GLPM = GLPM/XNORM
      GHPM = GHPM/XNORM
      GAT=1.D0/TGBET
      GAB=TGBET
      GLT=CA/SB
      GLB=-SA/CB
      GHT=SA/SB
      GHB=CA/CB
      GZAL=-CBMA
      GZAH=SBMA
      GLVV=SBMA
      GHVV=CBMA
      B=BET
      A=DATAN(SA/CA)
      IF(CA.LT.0D0)THEN
       IF(SA.LT.0D0)THEN
        A = A-PI
       ELSE
        A = A+PI
       ENDIF
      ENDIF
C ===============================================================
C ========== POLE MASSES 
      IF(IPOLE.EQ.1) THEN 
      MT=RUNM(AMT,6)
      MB=RUNM(AMT,5)
      SW2=1.D0-AMW**2/AMZ**2
C===== STOP MASSES
      MSTL2=AMSQ**2+(0.5D0-2.D0/3.D0*SW2)*AMZ**2*DCOS(2.D0*B)
      MSTR2=AMUR**2+2.D0/3.D0*SW2*AMZ**2*DCOS(2.D0*B)
      MLRT=AU-AMU/TGBET
      DELT=(MSTL2-MSTR2)**2+4*MT**2*MLRT**2
      MST12=MT**2+0.5D0*(MSTL2+MSTR2-DSQRT(DELT))
      MST22=MT**2+0.5D0*(MSTL2+MSTR2+DSQRT(DELT))
        IF(MST12.LT.0.D0)GOTO 111
      MST(1)=DSQRT(MST12)
      MST(2)=DSQRT(MST22)
      THET=0.5D0*DATAN(2.D0*MT*MLRT / (MSTL2-MSTR2) )
      IF(MSTL2.GT.MSTR2) THET = THET + PI/2
      CST= DCOS(THET)
      SST= DSIN(THET)
C===== SBOTTOM MASSES
      MSBL2=AMSQ**2+(-0.5D0+1.D0/3.D0*SW2)*AMZ**2*DCOS(2.D0*B)
      MSBR2=AMDR**2-1.D0/3.D0*SW2*AMZ**2*DCOS(2.D0*B)
      MLRB=AD-AMU*TGBET
      DELB=(MSBL2-MSBR2)**2+4*MB**2*MLRB**2
      MSB12=MB**2+0.5D0*(MSBL2+MSBR2-DSQRT(DELB))
      MSB22=MB**2+0.5D0*(MSBL2+MSBR2+DSQRT(DELB))
        IF(MSB12.LT.0.D0)GOTO 111
      MSB(1)=DSQRT(MSB12)
      MSB(2)=DSQRT(MSB22)
      THEB=0.5D0*DATAN(2.D0*MB*MLRB / (MSBL2-MSBR2) )
      IF(MSBL2.GT.MSBR2) THEB = THEB + PI/2
      CSB= DCOS(THEB)
      SSB= DSIN(THEB)
C===== LIGHT HIGGS COUPLINGS 
      GLTT(1,1)=-SBPA*(0.5D0*CST**2-2.D0/3.D0*SW2*DCOS(2*THET) )
     .    +MT**2/AMZ**2*GLT + MT*SST*CST/AMZ**2*(AU*GLT+AMU*GHT)
      GLTT(2,2)=-SBPA*(0.5D0*SST**2+2.D0/3.D0*SW2*DCOS(2*THET) )
     .    +MT**2/AMZ**2*GLT - MT*SST*CST/AMZ**2*(AU*GLT+AMU*GHT)
      GLTT(1,2)=-2*SBPA*SST*CST*(2.D0/3.D0*SW2-0.25D0)
     .    + MT*DCOS(2*THET)/2.D0/AMZ**2*(AU*GLT+AMU*GHT)
      GLTT(2,1)=-2*SBPA*SST*CST*(2.D0/3.D0*SW2-0.25D0)
     .    + MT*DCOS(2*THET)/2.D0/AMZ**2*(AU*GLT+AMU*GHT)
      GLBB(1,1)=-SBPA*(-0.5D0*CSB**2+1.D0/3.D0*SW2*DCOS(2*THEB))
     .    +MB**2/AMZ**2*GLB + MB*SSB*CSB/AMZ**2*(AD*GLB-AMU*GHB)
      GLBB(2,2)=-SBPA*(-0.5D0*SSB**2-1.D0/3.D0*SW2*DCOS(2*THEB))
     .    +MB**2/AMZ**2*GLB - MB*SSB*CSB/AMZ**2*(AD*GLB-AMU*GHB)
      GLBB(1,2)=-2*SBPA*SSB*CSB*(-1.D0/3.D0*SW2+0.25D0)
     .    + MB*DCOS(2*THEB)/2.D0/AMZ**2*(AD*GLB-AMU*GHB)
      GLBB(2,1)=-2*SBPA*SSB*CSB*(-1.D0/3.D0*SW2+0.25D0)
     .    + MB*DCOS(2*THEB)/2.D0/AMZ**2*(AD*GLB-AMU*GHB)
C===== HEAVY HIGGS COUPLINGS 
      GHTT(1,1)=CBPA*(0.5D0*CST**2-2.D0/3.D0*SW2*DCOS(2*THET))
     .    +MT**2/AMZ**2*GHT + MT*SST*CST/AMZ**2*(AU*GHT-AMU*GLT)
      GHTT(2,2)=CBPA*(0.5D0*SST**2+2.D0/3.D0*SW2*DCOS(2*THET))
     .    +MT**2/AMZ**2*GHT - MT*SST*CST/AMZ**2*(AU*GHT-AMU*GLT)
      GHTT(1,2)=2*CBPA*SST*CST*(2.D0/3.D0*SW2-0.25D0)
     .    +MT*DCOS(2*THET)/2.D0/AMZ**2*(AU*GHT-AMU*GLT)
      GHTT(2,1)=2*CBPA*SST*CST*(2.D0/3.D0*SW2-0.25D0)
     .    + MT*DCOS(2*THET)/2.D0/AMZ**2*(AU*GHT-AMU*GLT)
      GHBB(1,1)=CBPA*(-0.5D0*CSB**2+1.D0/3.D0*SW2*DCOS(2*THEB))
     .    +MB**2/AMZ**2*GHB + MB*SSB*CSB/AMZ**2*(AD*GHB+AMU*GLB)
      GHBB(2,2)=CBPA*(-0.5D0*SSB**2-1.D0/3.D0*SW2*DCOS(2*THEB))
     .    + MB**2/AMZ**2*GHB - MB*SSB*CSB/AMZ**2*(AD*GHB+AMU*GLB)
      GHBB(1,2)=2*CBPA*SSB*CSB*(-1.D0/3.D0*SW2+0.25D0)
     .    + MB*DCOS(2*THEB)/2.D0/AMZ**2*(AD*GHB+AMU*GLB)
      GHBB(2,1)=2*CBPA*SSB*CSB*(-1.D0/3.D0*SW2+0.25D0)
     .    + MB*DCOS(2*THEB)/2.D0/AMZ**2*(AD*GHB+AMU*GLB)
C===== PSEUDOSCALAR HIGGS COUPLINGS 
      GATT=-MT/2.D0/AMZ**2*(AMU+AU*GAT) 
      GABB=-MB/2.D0/AMZ**2*(AMU+AD*GAB) 
C======= LOOP CORRECTIONS  
      XDLT=GF/(2.D0*DSQRT(2.D0)*PI**2)*GLT**2*(-2.D0*MT**2+0.5D0*AML2)
     .    *DREAL(F0(MT,MT,AML2))
     .    *3*MT**2
      XDLB=GF/(2.D0*DSQRT(2.D0)*PI**2)*GLB**2*(-2.D0*MB**2+0.5D0*AML2)
     .    *DREAL(F0(MB,MB,AML2))
     .    *3*MB**2
C--BUG IN CARENA ET AL. FIXED
     .    +GF/(2.D0*DSQRT(2.D0)*PI**2)*GLB**2*(0.5D0*AML2)
     .    *DLOG(MB**2/MT**2)
     .    *3*MB**2
      XDHT=GF/(2.D0*DSQRT(2.D0)*PI**2)*GHT**2*(-2.D0*MT**2+0.5D0*AMH2)
     .    *DREAL(F0(MT,MT,AMH2))
     .    *3*MT**2
      XDHB=GF/(2.D0*DSQRT(2.D0)*PI**2)*GHB**2*(-2.D0*MB**2+0.5D0*AMH2)
     .    *DREAL(F0(MB,MB,AMH2))
     .    *3*MB**2
C--BUG IN CARENA ET AL. FIXED
     .    +GF/(2.D0*DSQRT(2.D0)*PI**2)*GHB**2*(0.5D0*AMH2)
     .    *DLOG(MB**2/MT**2)
     .    *3*MB**2
      XDAT=GF/(2.D0*DSQRT(2.D0)*PI**2)*GAT**2*(-0.5D0*AMA2)
     .    *DREAL(F0(MT,MT,AMA2))
     .    *3*MT**2
      XDAB=GF/(2.D0*DSQRT(2.D0)*PI**2)*GAB**2*(-0.5D0*AMA2)
     .    *DREAL(F0(MB,MB,AMA2))
     .    *3*MB**2
C--BUG IN CARENA ET AL. FIXED
     .    +GF/(2.D0*DSQRT(2.D0)*PI**2)*GAB**2*(-0.5D0*AMA2)
     .    *DLOG(MB**2/MT**2)
     .    *3*MB**2
      XDLST=0.D0
      XDLSB=0.D0
      XDHST=0.D0
      XDHSB=0.D0
         DO 311 I=1,2
         DO 311 J=1,2
      XDLST=XDLST+GF/(2.D0*DSQRT(2.D0)*PI**2)*GLTT(I,J)**2*
     .      DREAL(F0(MST(I),MST(J),AML2))
     .    *3*AMZ**4
      XDLSB=XDLSB+GF/(2.D0*DSQRT(2.D0)*PI**2)*GLBB(I,J)**2*
     .      DREAL(F0(MSB(I),MSB(J),AML2))
     .    *3*AMZ**4
      XDHST=XDHST+GF/(2.D0*DSQRT(2.D0)*PI**2)*GHTT(I,J)**2*
     .      DREAL(F0(MST(I),MST(J),AMH2))
     .    *3*AMZ**4
      XDHSB=XDHSB+GF/(2.D0*DSQRT(2.D0)*PI**2)*GHBB(I,J)**2*
     .      DREAL(F0(MSB(I),MSB(J),AMH2))
     .    *3*AMZ**4
311    CONTINUE
      XDAST=GF/(1.D0*DSQRT(2.D0)*PI**2)*GATT**2*
     .      DREAL(F0(MST(1),MST(2),AMA2))
     .    *3*AMZ**4
      XDASB=GF/(1.D0*DSQRT(2.D0)*PI**2)*GABB**2*
     .      DREAL(F0(MSB(1),MSB(2),AMA2))
     .    *3*AMZ**4
      
       AML=DSQRT(AML2+XDLT+XDLB+XDLST+XDLSB)
       AMH=DSQRT(AMH2+XDHT+XDHB+XDHST+XDHSB)  
       AMA=DSQRT(AMA2+XDAT+XDAB+XDAST+XDASB)  
      ELSE
       AML=AMLR
       AMH=AMHR     
       AMA=AMAR     
      ENDIF 
      RETURN
111   STOP
      END

C ===================== THE FUNCTION F0 ===============
      COMPLEX FUNCTION F0*16(M1,M2,QSQ)
      IMPLICIT REAL*8 (A-H,M,O-Z)
      COMPLEX*16 CD,CR,CQ2,IEPS,CBET,CXX
      M1SQ = M1*M1
      M2SQ = M2*M2
      AQSQ = DABS(QSQ)
      IEPS = DCMPLX(1.D0,1.D-12)
      CQ2 = QSQ*IEPS
      CD = (M1SQ-M2SQ)/CQ2
      CR = CDSQRT((1+CD)**2 - 4*M1SQ/CQ2)
      IF(QSQ.EQ.0.D0) THEN
       F0 = 0.D0
      ELSE
       IF(M1.EQ.M2) THEN
        F0 = -2.D0 + CR*CDLOG(-(1+CR)/(1-CR))
       ELSE
        CBET = CDSQRT(1-4*M1*M2/(CQ2 - (M1-M2)**2))
        CXX = (CBET-1)/(CBET+1)
        F0 = -1 + ((QSQ+M2SQ-M1SQ)/2/QSQ - M2SQ/(M2SQ-M1SQ))
     .                                           *DLOG(M2SQ/M1SQ)
     .     - (QSQ-(M1-M2)**2)/QSQ*CBET*CDLOG(CXX)
       ENDIF
      ENDIF
      RETURN
      END

C     ************************************************************
C     SUBROUTINE FOR HSM ---> V*V* ---> 4F
C     ************************************************************
      SUBROUTINE HTOVV(AMH,AMV,GAMV,HTVV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/VVOFF/AMH1,AMV1,GAMV1
      COMMON/PREC/IP
      EXTERNAL FTOVV1
      IP=20
      AMH1=AMH
      AMV1=AMV
      GAMV1=GAMV
      DLT=1D0/IP
      SUM=0D0
      DO 1 I=1,IP
       UU=DLT*I
       DD=UU-DLT
       CALL QGAUS1(FTOVV1,DD,UU,RES)
       SUM=SUM+RES
1     CONTINUE
      HTVV=SUM
      RETURN
      END

      DOUBLE PRECISION FUNCTION FTOVV1(XX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/FIRST/X1
      COMMON/PREC/IP
      EXTERNAL FTOVV2
      X1=XX
      DLT=1D0/IP
      SUM=0D0
      DO 1 I=1,IP
       UU=DLT*I
       DD=UU-DLT
       CALL QGAUS2(FTOVV2,DD,UU,RES)
       SUM=SUM+RES
1     CONTINUE
      FTOVV1=SUM
      RETURN
      END

      DOUBLE PRECISION FUNCTION FTOVV2(XX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION YY(2)
      COMMON/FIRST/X1
      YY(1)=X1
      YY(2)=XX
      FTOVV2=FTOVV(YY)
      RETURN
      END

      DOUBLE PRECISION FUNCTION FTOVV(XX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION LAMB
      DIMENSION XX(2)
      COMMON/VVOFF/AMH,AMW,GAMW
      LAMB(X,Y)=DSQRT((1.D0-X-Y)**2-4.D0*X*Y)
      PI=4D0*DATAN(1D0)
      ICASE = 1
      IF(ICASE.EQ.0)THEN
       YY = AMH**2
       Y1 = DATAN((YY-AMW**2)/AMW/GAMW)
       Y2 = -DATAN((AMW**2)/AMW/GAMW)
       DJAC = Y1-Y2
       T1 = TAN(Y1*XX(1)+Y2*(1.D0-XX(1)))
       SP = AMW**2 + AMW*GAMW*T1
       YY = (AMH-DSQRT(SP))**2
       Y1 = DATAN((YY-AMW**2)/AMW/GAMW)
       Y2 = -DATAN((AMW**2)/AMW/GAMW)
       DJAC = DJAC*(Y1-Y2)
       T2 = TAN(Y1*XX(2)+Y2*(1.D0-XX(2)))
       SM = AMW**2 + AMW*GAMW*T2
       AM2=AMH**2
       GAM = AM2*LAMB(SP/AM2,SM/AM2)*(1+LAMB(SP/AM2,SM/AM2)**2*AMH**4
     .                               /SP/SM/12)
       PRO1 = SP/AMW**2
       PRO2 = SM/AMW**2
       FTOVV = PRO1*PRO2*GAM*DJAC/PI**2
      ELSE
       SP = AMH**2*XX(1)
       SM = (AMH-DSQRT(SP))**2*XX(2)
       DJAC = AMH**2*(AMH-DSQRT(SP))**2/PI**2
       AM2=AMH**2
       GAM = AM2*LAMB(SP/AM2,SM/AM2)*(1+LAMB(SP/AM2,SM/AM2)**2*AMH**4
     .                               /SP/SM/12)
       PRO1 = SP*GAMW/AMW/((SP-AMW**2)**2+AMW**2*GAMW**2)
       PRO2 = SM*GAMW/AMW/((SM-AMW**2)**2+AMW**2*GAMW**2)
       FTOVV = PRO1*PRO2*GAM*DJAC
      ENDIF
      RETURN
      END

C     ************************************************************
C     SUBROUTINE FOR HSM ---> TT* ---> TBW
C     ************************************************************
      SUBROUTINE HTOTTS(AMH,AMT,AMB,AMW,HTTS)
      IMPLICIT REAL*8(A-Z)
      INTEGER IP,K
      COMMON/PREC1/IP
      EXTERNAL FUNSTT1
      COMMON/IKSY0/X1,X2,M1,M2,M3,ECM,S
      COMMON/TOP0/AMH0,AMT0,AMB0,AMW0
      AMH0=AMH
      AMT0=AMT
      AMB0=AMB
      AMW0=AMW
      IP=5
      M1=AMB
      M2=AMT
      M3=AMW
C     FIRST INTEGRATE OVER X2, i.e. (1+3) SYSTEM
C        CHECK WHETHER ENOUGH PHASE SPACE
      MASTOT=M1+M2+M3
      IF(MASTOT.GE.AMH) GOTO 12
      ECM=AMH
      S=ECM**2
      U1=(ECM-M2)**2
      D1=(M1+M3)**2
      U=(S-D1+M2**2)/s
      D=(S-U1+M2**2)/s
      DEL=(U-D)/IP
      U=D+DEL
      XSEC=0.D0
      DO K=1,IP
      CALL QGAUS1(FUNSTT1,D,U,SS)
      D=U
      U=D+DEL
      XSEC=XSEC+SS
      ENDDO
      HTTS=XSEC
12    CONTINUE
      RETURN
      END

      DOUBLE PRECISION FUNCTION FUNSTT1(XL)
      IMPLICIT REAL*8(A-Z)
      INTEGER IP,I
      COMMON/IKSY0/X1,X2,M1,M2,M3,ECM,S
      COMMON/PREC1/IP
      EXTERNAL FUNSTT2
      X2=XL
      S13=S-S*X2+M2**2
      TEM=2.D0*DSQRT(S13)
      E2S=(S-S13-M2**2)/TEM
      E3S=(S13+M3**2-M1**2)/TEM
C     SECOND INTEGRAL OVER X1, i.e. (2+3) SYSTEM
      U1=(E2S+E3S)**2-(DSQRT(E2S**2-M2**2)-DSQRT(E3S**2-M3**2))**2
      D1=(E2S+E3S)**2-(DSQRT(E2S**2-M2**2)+DSQRT(E3S**2-M3**2))**2
      U=(S-D1+M1**2)/s
      D=(S-U1+M1**2)/s
      DEL=(U-D)/IP
      FUNSTT1=0.d0
      U=D+DEL
      DO I=1,IP
      CALL QGAUS2(FUNSTT2,D,U,SS)
      FUNSTT1=FUNSTT1+SS
      D=U
      U=D+DEL
      ENDDO
      RETURN
      END

      DOUBLE PRECISION FUNCTION FUNSTT2(XK)
      IMPLICIT REAL*8(A-Z)
      COMMON/IKSY0/X1,X2,M1,M2,M3,ECM,S
      X1=XK
      CALL ELEMSTT(SS)
      FUNSTT2=SS
      RETURN
      END

      SUBROUTINE ELEMSTT(RES)
      IMPLICIT REAL*8(A-Z)
      COMMON/IKSY0/X1,X2,M1,M2,M3,ECM,S
      COMMON/TOP0/AMH,AMT,AMB,AMW
      COMMON/WZWDTH/GAMC0,GAMT0,GAMT1,GAMW0,GAMZ0
      GAMT=GAMT0**2*AMT**2/AMH**4
      GAMW=GAMW0**2*AMW**2/AMH**4
      W=AMW**2/AMH**2
      T=AMT**2/AMH**2
      Y1=1-X2
      Y2=1-X1
      X0=2.D0-X1-X2
      W1=(1.D0-X2)
      W3=(1.-X1-X2)
      W11=1.D0/((1.D0-X2)**2+GAMT)
      W33=1.D0/(W3**2+GAMW**2)
      W13=W1*W3*W11*W33

      R11=4*T*W-16.*T*W*Y1-4.*T*Y2*Y1+8.*T*Y1+32.*T*W**2-20
     . .*T*Y1**2+8.*W*Y2*Y1+4.*W*Y1**2-4.*Y2*Y1**2-16.*T**2*W-
     .  32.*T**2*Y1+4.*T**2-16.*T**3-8.*W**2+4.*Y1**2-4.*Y1**3
      R33=-4.*T*W+4.*T*W*Y2-2.*T*W*Y2*Y1+4.*T*W*Y1+T*W*Y2**2-
     .  3.*T*W*Y1**2+2.*T*Y2*Y1-3.*T*Y2*Y1**2+4.*T*W**2-4.*T*W**3
     .  +T*Y2**2-3.*T*Y2**2*Y1-T*Y2**3+T*Y1**2-T*Y1**3+4.*T**2
     .  *W-4.*T**2*W*Y2-4.*T**2*W*Y1-2.*T**2*Y2*Y1-4.*T**2*W**2-
     .  T**2*Y2**2-T**2*Y1**2+4.*W**2*Y2*Y1-8.*W**3*Y2-8.*W**3*Y1
     .  +4.*W**3+8.*W**4
      R13=8.*W-24.*T*W+16.*T*W*Y1 -4.*T*Y2+16.*T*Y2*Y1-4.*T*
     .  Y1+16.*T*W**2+4.*T*Y2**2+12.*T*Y1**2-8.*W*Y2-12.*W*Y2*Y1
     .  -8.*W*Y1+4.*W*Y1**2-4.*Y2*Y1+8.*Y2*Y1**2+16.*T**2*W+8.
     .  *T**2*Y2+8.*T**2*Y1+16.*W**2*Y2+24.*W**2*Y1+4.*Y2**2*Y1-
     .  32.*W**3-4.*Y1**2+4.*Y1**3
      RES=R11*W11+4.D0*R33*W33/T-2.D0*R13*W13
      RETURN
      END

C     **************************************************
C     SUBROUTINE FOR A -> TT* -> TBW
C     **************************************************

      SUBROUTINE ATOTT(AMA,AMT,AMB,AMW,AMCH,ATT0)
      IMPLICIT REAL*8(A-Z)
      INTEGER IP,K
      COMMON/PREC1/IP
      EXTERNAL FUNATT1
      COMMON/IKSY1/X1,X2,M1,M2,M3,ECM,S
      COMMON/TOP1/AMA1,AMT1,AMB1,AMW1,AMCH1
      AMA1=AMA
      AMT1=AMT
      AMB1=AMB
      AMW1=AMW
      AMCH1=AMCH
      IP=5
      M1=AMB
      M2=AMT
      M3=AMW
C        FIRST INTEGRATE OVER X2, i.e. (1+3) SYSTEM
C        CHECK WHETHER ENOUGH PHASE SPACE
      MASTOT=M1+M2+M3
      IF(MASTOT.GE.AMA) GOTO 12
      ECM=AMA
      S=ECM**2
      U1=(ECM-M2)**2
      D1=(M1+M3)**2
      U=(S-D1+M2**2)/s
      D=(S-U1+M2**2)/s
      DEL=(U-D)/IP
      U=D+DEL
      XSEC=0.D0
      DO K=1,IP
      CALL QGAUS1(FUNATT1,D,U,SS)
      D=U
      U=D+DEL
      XSEC=XSEC+SS
      ENDDO
      ATT0=XSEC
 12   CONTINUE
      RETURN
      END

      DOUBLE PRECISION FUNCTION FUNATT1(XL)
      IMPLICIT REAL*8(A-Z)
      INTEGER IP,I
      COMMON/IKSY1/X1,X2,M1,M2,M3,ECM,S
      COMMON/PREC1/IP
      EXTERNAL FUNATT2
      X2=XL
      S13=S-S*X2+M2**2
      TEM=2.D0*DSQRT(S13)
      E2S=(S-S13-M2**2)/TEM
      E3S=(S13+M3**2-M1**2)/TEM
C     SECOND INTEGRAL OVER X1, i.e. (2+3) SYSTEM
      U1=(E2S+E3S)**2-(DSQRT(E2S**2-M2**2)-DSQRT(E3S**2-M3**2))**2
      D1=(E2S+E3S)**2-(DSQRT(E2S**2-M2**2)+DSQRT(E3S**2-M3**2))**2
      U=(S-D1+M1**2)/s
      D=(S-U1+M1**2)/s
      DEL=(U-D)/IP
      FUNATT1=0.d0
      U=D+DEL
      DO I=1,IP
      CALL QGAUS2(FUNATT2,D,U,SS)
      FUNATT1=FUNATT1+SS
      D=U
      U=D+DEL
      ENDDO
      RETURN
      END

      DOUBLE PRECISION FUNCTION FUNATT2(XK)
      IMPLICIT REAL*8(A-Z)
      COMMON/IKSY1/X1,X2,M1,M2,M3,ECM,S
      X1=XK
      CALL ELEMATT(SS)
      FUNATT2=SS
      RETURN
      END

      SUBROUTINE ELEMATT(RES)
      IMPLICIT REAL*8(A-Z)
      COMMON/IKSY1/X1,X2,M1,M2,M3,ECM,S
      COMMON/TOP1/AMA,AMT,AMB,AMW,AMCH
      COMMON/WZWDTH/GAMC0,GAMT0,GAMT1,GAMW,GAMZ
      GAMT=GAMT1**2*AMT**2/AMA**4
      GAMC=GAMC0**2*AMCH**2/AMA**4
      CH=AMCH**2/AMA**2
      W=AMW**2/AMA**2
      T=AMT**2/AMA**2
      Y1=1-X1
      Y2=1-X2
      X0=2.D0-X1-X2
      W1=(1.D0-x2)
      W2=(1.D0-X0+W-CH)
      W22=1.D0/ ((1.D0-X0+W-CH)**2+GAMC)
      W11=1.D0/((1.D0-X2)**2+GAMT)
      W12=W1*W2*W11*W22
      R11=4.D0*T*W-4.D0*T*Y1*Y2+8.D0*T*Y2-4.D0*T*Y2**2+8.D0*W*Y1*Y2+4.D0
     .  *W*Y2**2-4.D0*Y1*Y2**2+4.D0*T**2-8.D0*W**2+4.D0*Y2**2-4.D0*Y2**3
      R22=-16.D0*W+16.D0*T*W-8.D0*T*Y1*Y2-4.D0*T*Y1**2-4.D0*T*Y2**2+16.
     .D0*W*Y1+8.D0*W*Y1*Y2+16.D0*W*Y2+4.D0*W*Y1**2+4.D0*W*Y2**2+8.D0*Y1*
     . Y2-12.D0*Y1*Y2**2-12.D0*Y1**2*Y2-16.D0*W**2+4.D0*Y1**2-4.D0*Y1**3
     . +4.D0*Y2**2-4.D0*Y2**3
      R12=16.D0*W-16.D0*T*W-8.D0*T*Y1+16.D0*T*Y1*Y2-8.D0*T*Y2+8.D0*T*Y1
     . **2+8.D0*T*Y2**2-16.D0*W*Y1-8.D0*W*Y1*Y2-16.D0*W*Y2-8.D0*W*Y2**2-
     . 8.D0*Y1*Y2+16.D0*Y1*Y2**2+8.D0*Y1**2*Y2+16.D0*W**2-8.D0*Y2**2
     . +8.D0*Y2**3
      RES=R11*W11+R22*W22+R12*W12
      RETURN
      END

C     ************************************************************
C     SUBROUTINE FOR H ---> TT* ---> TBW
C     ************************************************************
      SUBROUTINE HTOTT(AMH,AMT,AMB,AMW,AMCH,TB,GHT,GAT,GHVV,HTT0)
      IMPLICIT REAL*8(A-Z)
      INTEGER IP,K
      COMMON/PREC1/IP
      EXTERNAL FUNHTT1
      COMMON/IKSY2/X1,X2,M1,M2,M3,ECM,S
      COMMON/TOP2/AMH2,AMT2,AMB2,AMW2,AMCH2,TB2,GHT2,GAT2,GHVV2
      AMH2=AMH
      AMT2=AMT
      AMB2=AMB
      AMW2=AMW
      AMCH2=AMCH
      TB2=TB
      GHT2=GHT
      GAT2=GAT
      GHVV2=GHVV
      IP=5
      M1=AMB
      M2=AMT
      M3=AMW
C     FIRST INTEGRATE OVER X2, i.e. (1+3) SYSTEM
C        CHECK WHETHER ENOUGH PHASE SPACE
      MASTOT=M1+M2+M3
      IF(MASTOT.GE.AMH) GOTO 12
      ECM=AMH
      S=ECM**2
      U1=(ECM-M2)**2
      D1=(M1+M3)**2
      U=(S-D1+M2**2)/s
      D=(S-U1+M2**2)/s
      DEL=(U-D)/IP
      U=D+DEL
      XSEC=0.D0
      DO K=1,IP
      CALL QGAUS1(FUNHTT1,D,U,SS)
      D=U
      U=D+DEL
      XSEC=XSEC+SS
      ENDDO
      HTT0=XSEC
 12   CONTINUE
      RETURN
      END

      DOUBLE PRECISION FUNCTION FUNHTT1(XL)
      IMPLICIT REAL*8(A-Z)
      INTEGER IP,I
      COMMON/IKSY2/X1,X2,M1,M2,M3,ECM,S
      COMMON/PREC1/IP
      EXTERNAL FUNHTT2
      X2=XL
      S13=S-S*X2+M2**2
      TEM=2.D0*DSQRT(S13)
      E2S=(S-S13-M2**2)/TEM
      E3S=(S13+M3**2-M1**2)/TEM
C     SECOND INTEGRAL OVER X1, i.e. (2+3) SYSTEM
      U1=(E2S+E3S)**2-(DSQRT(E2S**2-M2**2)-DSQRT(E3S**2-M3**2))**2
      D1=(E2S+E3S)**2-(DSQRT(E2S**2-M2**2)+DSQRT(E3S**2-M3**2))**2
      U=(S-D1+M1**2)/s
      D=(S-U1+M1**2)/s
      DEL=(U-D)/IP
      FUNHTT1=0.d0
      U=D+DEL
      DO I=1,IP
      CALL QGAUS2(FUNHTT2,D,U,SS)
      FUNHTT1=FUNHTT1+SS
      D=U
      U=D+DEL
      ENDDO
      RETURN
      END

      DOUBLE PRECISION FUNCTION FUNHTT2(XK)
      IMPLICIT REAL*8(A-Z)
      COMMON/IKSY2/X1,X2,M1,M2,M3,ECM,S
      X1=XK
      CALL ELEMHTT(SS)
      FUNHTT2=SS
      RETURN
      END

      SUBROUTINE ELEMHTT(RES)
      IMPLICIT REAL*8(A-Z)
      COMMON/IKSY2/X1,X2,M1,M2,M3,ECM,S
      COMMON/TOP2/AMH,AMT,AMB,AMW,AMCH,TB,GHT,GAT,GHVV
      COMMON/WZWDTH/GAMC0,GAMT0,GAMT1,GAMW0,GAMZ0
      GAMT=GAMT1**2*AMT**2/AMH**4
      GAMC=GAMC0**2*AMCH**2/AMH**4
      GAMW=GAMW0**2*AMW**2/AMH**4
      CH=AMCH**2/AMH**2
      W=AMW**2/AMH**2
      T=AMT**2/AMH**2
      Y1=1-X2
      Y2=1-X1
      X0=2.D0-X1-X2
      W1=(1.D0-X2)
      W2=(1.D0-X0+W-CH)
      W3=-(1.-X1-X2)
      W22=1.D0/ ((1.D0-X0+W-CH)**2+GAMC)
      W11=1.D0/((1.D0-X2)**2+GAMT)
      W33=1.D0/(W3**2+GAMW**2)
      W12=W1*W2*W11*W22
      W13=W1*W3*W11*W33
      W23=W2*W3*W22*W33

      R11=4*T*W-16.*T*W*Y1-4.*T*Y2*Y1+8.*T*Y1+32.*T*W**2-20
     . .*T*Y1**2+8.*W*Y2*Y1+4.*W*Y1**2-4.*Y2*Y1**2-16.*T**2*W-
     .  32.*T**2*Y1+4.*T**2-16.*T**3-8.*W**2+4.*Y1**2-4.*Y1**3
      R22=-16.*W+16.*T*W-8.*T*Y2*Y1-4.*T*Y2**2-4.*T*Y1**2+16
     .  .*W*Y2 + 8.*W*Y2*Y1 + 16.*W*Y1 + 4.*W*Y2**2 + 4.*W*Y1**2+8.*Y2*
     .  Y1-12.*Y2*Y1**2-12.*Y2**2*Y1-16.*W**2+4.*Y2**2-4.*Y2**3
     .  +4.*Y1**2-4.*Y1**3
      R33=-4.*T*W+4.*T*W*Y2-2.*T*W*Y2*Y1+4.*T*W*Y1+T*W*Y2**2-
     .  3.*T*W*Y1**2+2.*T*Y2*Y1-3.*T*Y2*Y1**2+4.*T*W**2-4.*T*W**3
     .  +T*Y2**2-3.*T*Y2**2*Y1-T*Y2**3+T*Y1**2-T*Y1**3+4.*T**2
     .  *W-4.*T**2*W*Y2-4.*T**2*W*Y1-2.*T**2*Y2*Y1-4.*T**2*W**2-
     .  T**2*Y2**2-T**2*Y1**2+4.*W**2*Y2*Y1-8.*W**3*Y2-8.*W**3*Y1
     .  +4.*W**3+8.*W**4
      R12=-16.*W+48.*T*W-16.*T*W*Y2+16.*T*W*Y1+8.*T*Y2-32.*T
     .  *Y2*Y1+8.*T*Y1-8.*T*Y2**2 - 24.*T*Y1**2+16.*W*Y2+8.*W*Y2*
     .  Y1+16.*W*Y1+8.*W*Y1**2+8.*Y2*Y1-16.*Y2*Y1**2-16.*T**2*Y2
     .  -16.*T**2*Y1-8.*Y2**2*Y1-16.*W**2+8.*Y1**2-8.*Y1**3
      R13=8.*W-24.*T*W+16.*T*W*Y1 -4.*T*Y2+16.*T*Y2*Y1-4.*T*
     .  Y1+16.*T*W**2+4.*T*Y2**2+12.*T*Y1**2-8.*W*Y2-12.*W*Y2*Y1
     .  -8.*W*Y1+4.*W*Y1**2-4.*Y2*Y1+8.*Y2*Y1**2+16.*T**2*W+8.
     .  *T**2*Y2+8.*T**2*Y1+16.*W**2*Y2+24.*W**2*Y1+4.*Y2**2*Y1-
     .  32.*W**3-4.*Y1**2+4.*Y1**3
      R23=16.*W-16.*T*W+8.*T*W*Y2+8.*T*W*Y1+8.*T*Y2*Y1+4.*T*
     .  Y2**2+4.*T*Y1**2-16.*W*Y2-16.*W*Y1-4.*W*Y2**2+4.*W*Y1**2
     .  -8.*Y2*Y1+12.*Y2*Y1**2+8.*W**2*Y2-8.*W**2*Y1+12.*Y2**2*
     .  Y1-4.*Y2**2+4.*Y2**3-4.*Y1**2+4.*Y1**3
      GLVV=DSQRT(1.D0-GHVV**2)
      RES=GHT**2*R11*W11+GLVV**2*GAT**2*R22*W22+
     .    4.D0*GHVV**2*R33*W33/T+2.D0*GHT*GLVV*GAT*R12*W12+
     .    2.D0*GHT*GHVV*R13*W13+2.D0*GHVV*GLVV*GAT*R23*W23
      RETURN
      END

C     ************************************************************
C     SUBROUTINE FOR H+ ---> BT* ---> BBW
C     ************************************************************
      SUBROUTINE CTOTT(AMCH,AMT,AMB,AMW,CTT0)
      IMPLICIT REAL*8(A-Z)
      INTEGER IP,K
      COMMON/PREC1/IP
      EXTERNAL FUNCTT1
      COMMON/IKSY3/X1,X2,M1,M2,M3,ECM,S
      COMMON/TOP3/AMH3,AMT3,AMB3,AMW3
      AMH3=AMCH
      AMT3=AMT
      AMB3=AMB
      AMW3=AMW
      IP=5
      M1=AMB
      M2=AMB
      M3=AMW
C     FIRST INTEGRATE OVER X2, i.e. (1+3) SYSTEM
C        CHECK WHETHER ENOUGH PHASE SPACE
      MASTOT=M1+M2+M3
      IF(MASTOT.GE.AMCH) GOTO 12
      ECM=AMCH
      S=ECM**2
      U1=(ECM-M2)**2
      D1=(M1+M3)**2
      U=(S-D1+M2**2)/s
      D=(S-U1+M2**2)/s
      DEL=(U-D)/IP
      U=D+DEL
      XSEC=0.D0
      DO K=1,IP
      CALL QGAUS1(FUNCTT1,D,U,SS)
      D=U
      U=D+DEL
      XSEC=XSEC+SS
      ENDDO
      CTT0=XSEC
12    CONTINUE
      RETURN
      END

      DOUBLE PRECISION FUNCTION FUNCTT1(XL)
      IMPLICIT REAL*8(A-Z)
      INTEGER IP,I
      COMMON/IKSY3/X1,X2,M1,M2,M3,ECM,S
      COMMON/PREC1/IP
      EXTERNAL FUNCTT2
      X2=XL
      S13=S-S*X2+M2**2
      TEM=2.D0*DSQRT(S13)
      E2S=(S-S13-M2**2)/TEM
      E3S=(S13+M3**2-M1**2)/TEM
C     SECOND INTEGRAL OVER X1, i.e. (2+3) SYSTEM
      U1=(E2S+E3S)**2-(DSQRT(E2S**2-M2**2)-DSQRT(E3S**2-M3**2))**2
      D1=(E2S+E3S)**2-(DSQRT(E2S**2-M2**2)+DSQRT(E3S**2-M3**2))**2
      U=(S-D1+M1**2)/s
      D=(S-U1+M1**2)/s
      DEL=(U-D)/IP
      FUNCTT1=0.d0
      U=D+DEL
      DO I=1,IP
      CALL QGAUS2(FUNCTT2,D,U,SS)
      FUNCTT1=FUNCTT1+SS
      D=U
      U=D+DEL
      ENDDO
      RETURN
      END

      DOUBLE PRECISION FUNCTION FUNCTT2(XK)
      IMPLICIT REAL*8(A-Z)
      COMMON/IKSY3/X1,X2,M1,M2,M3,ECM,S
      X1=XK
      CALL ELEMCTT(SS)
      FUNCTT2=SS
      RETURN
      END

      SUBROUTINE ELEMCTT(RES)
      IMPLICIT REAL*8(A-Z)
      COMMON/IKSY3/X1,X2,M1,M2,M3,ECM,S
      COMMON/TOP3/AMCH,AMT,AMB,AMW
      COMMON/WZWDTH/GAMC0,GAMT0,GAMT1,GAMW,GAMZ
      GAMT=GAMT1**2*AMT**2/AMCH**4
      W=AMW**2/AMCH**2
      T=AMT**2/AMCH**2
      B=AMB**2/AMCH**2
      RES=((1.D0-X1-W)*(1.D0-X2-W)+W*(X1+X2-1.D0+W))/
     .   ((1.D0-X2+B-T)**2+GAMT)
      RETURN
      END

C   *****************  INTEGRATION ROUTINE ***********************
C    Returns SS as integral of FUNC from A to B, by 10-point Gauss-
C    Legendre integration
      SUBROUTINE QGAUS1(FUNC,A,B,SS)
      IMPLICIT REAL*8(A-Z)
      INTEGER J
      DIMENSION X(5),W(5)
      EXTERNAL FUNC
      DATA X/.1488743389D0,.4333953941D0,.6794095682D0
     .  ,.8650633666D0,.9739065285D0/
      DATA W/.2955242247D0,.2692667193D0,.2190863625D0
     .  ,.1494513491D0,.0666713443D0/
      XM=0.5D0*(B+A)
      XR=0.5D0*(B-A)
      SS=0.D0
      DO 11 J=1,5
        DX=XR*X(J)
        SS=SS+W(J)*(FUNC(XM+DX)+FUNC(XM-DX))
11    CONTINUE
      SS=XR*SS
      RETURN
      END

C     Returns SS as integral of FUNC from A to B, by 10-point Gauss-
C      Legendre integration
      SUBROUTINE QGAUS2(FUNC,A,B,SS)
      IMPLICIT REAL*8(A-Z)
      INTEGER J
      DIMENSION X(5),W(5)
      EXTERNAL FUNC
      DATA X/.1488743389D0,.4333953941D0,.6794095682D0
     .  ,.8650633666D0,.9739065285D0/
      DATA W/.2955242247D0,.2692667193D0,.2190863625D0
     .  ,.1494513491D0,.0666713443D0/
      XM=0.5D0*(B+A)
      XR=0.5D0*(B-A)
      SS=0.D0
      DO 11 J=1,5
        DX=XR*X(J)
        SS=SS+W(J)*(FUNC(XM+DX)+FUNC(XM-DX))
11    CONTINUE
      SS=XR*SS
      RETURN
      END
C
C=====================================================================
C
C  New h,H,A mass program by Carena, Quiros & Wagner
C
C=====================================================================
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   
C    This subroutine computes the CP-even Higgs and CP-odd pole
C    Higgs masses and mixing angles. 
C
C    Inputs: ihiggs(explained below),mchi,ma,tanb,mq,mur,mdr,mtop,
C    at,ab,mu
C
C    where mchi is the largest chargino mass, ma is the running
C    CP-odd Higgs mass, tanb is the value of the ratio of vacuum
C    expectaion values at the scale mtop,mq is the third generation
C    left handed squark mass parameter, mur is the third generation
C    right handed stop mass parameter, mdr is the third generation
C    right handed sbottom mass parameter, mtop is the pole top quark
C    mass; at,ab are the soft supersymmetry breaking trilinear 
C    couplings of the stop and sbottoms, respectively, and mu is the
C    supersymmetric mass parameter. 
C
C
C    Output: mh and mhp which are the lightest CP-even Higgs running
C    and pole masses, respectively; hm and hmp are the heaviest CP-even
C    Higgs running and pole masses, repectively; sa and ca are the
C    sin(alpha) and cos(alpha) where alpha is the Higgs mixing angle.   
C    amp is the CP-odd Higgs pole mass. stop1,stop2,sbot1 and sbot2
C    are the stop and sbottom mass eigenvalues. Finally, tanbA is
C    the value of tanb at the CP-odd Higgs mass scale.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCcC
ccccccccccccccccccccccccccccccccccccccccccc
ccccc  The parameter ihiggs=0,1,2,3 corresponds to the
ccccc  number of Higgses whose pole mass is computed 
ccccc   by the subroutine vac(...). If ihiggs=0 only running
ccccc   masses are given, what makes the running of the program
ccccc   much faster and it is quite generally a good approximation
ccccc   (for a theoretical discussion see Ref. below). 
ccccc    If ihiggs=1, only the pole
ccccc   mass for h is computed. If ihiggs=2, then h and H, and
ccccc   if ihiggs=3, then h,H,A polarizations are computed
ccccccccccccccccccccccccccccccccccccccccccccccc 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Program based on the work by M. Carena, M. Quiros
c       and C.E.M. Wagner, "Effective potential methods and
c       the Higgs mass spectrum in the MSSM", CERN-TH/95-157,
c       
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


 	subroutine pole(ihiggs,mchi,ma,tanb,mq,mur,mdr,mtop,at,ab,mu,
     *        mh,mhp,hm,hmp,amp,sa,ca,sab,cab,stop1,stop2,sbot1,sbot2,
     *  tanbA,
     *  mglu)
C 	subroutine pole(ihiggs,mchi,ma,tanb,mq,mur,mdr,mtop,at,ab,mu,
C     *        mhp,hmp,amp,sa,ca,
C     *  tanbA, mglu)

	implicit real*8(a-h,m,o-z)
	call vac(ihiggs,mchi,ma,tanb,mq,mur,mdr,mtop,at,ab,mu,
     *        mh,mhp,hm,hmp,amp,stop1,stop2,sbot1,sbot2,
     *         sa,ca,sab,cab,stop1w,stop2w,tanbA,mglu)
	sinb = tanb/(tanb**2+1.)**.5
	cosb = 1./(tanb**2+1.)**.5
	sinbma = sinb*ca - cosb*sa
	return
	end

       subroutine vac(ihiggs,mchi,ma,tanb,mq,mur,mdr,
     *  mtop,at,ab,mu,mh,mhp,hm,hmp,amp,stop1,stop2,
     *  sbot1,sbot2,sa,ca,sab,cab,stop1w,stop2w,tanbA,mglu)
       implicit real*8(a-h,m,o-z)

       dimension delta(2,2),coupt(2,2),T(2,2),sstop2(2),
     *ssbot2(2),B(2,2),coupb(2,2),
     *hcoupt(2,2),hcoupb(2,2),
     *acoupt(2,2),acoupb(2,2),pr(3), polar(3)

      delta(1,1) = 1.
      delta(2,2) = 1.
      delta(1,2) = 0.
      delta(2,1) = 0.
       v = 174.1
  	mz=91.18
	pi=3.14159
	alpha3z=.12
     	alpha3=1./(1./alpha3z+23./6./pi*log(mtop/mz))	

         rmtop = mtop/(1.+4*alpha3/3./pi)

	
       ht = rmtop /v

       call rghm(mchi,ma,tanb,mq,mur,mdr,mtop,at,ab,
     *   mu,mh,hm,sa,ca,sab,cab,tanbA,mglu,deltamt,deltamb)

       sinb = tanb/(tanb**2+1.)**.5
       cosb = 1./(tanb**2+1.)**.5
       cos2b = sinb**2 - cosb**2
       sinbpa = sinb*ca + cosb*sa
       cosbpa = cosb*ca - sinb*sa	
       rmbot = 3./(1.+deltamb)
       mq2 = mq**2
	if(abs(mur).lt.0.1) mur  = 0.1
        mur2 = abs(mur)/mur*mur**2
        mdr2 = mdr**2
        mst11 = Rmtop**2 + mq2  - 0.35*MZ**2*cos2b
        mst22 = Rmtop**2 + mur2 - 0.15*MZ**2*cos2b
        if(mst11.lt.0.) goto 3333
        if(mst22.lt.0.) goto 3333
        msb11 = Rmbot**2 + mq2  + 0.42*MZ**2*cos2b
        msb22 = Rmbot**2 + mdr2 + 0.08*MZ**2*cos2b
        if(msb11.lt.0.) goto 3333
        if(msb22.lt.0.) goto 3333
        wmst11 = Rmtop**2 + mq2
        wmst22 = Rmtop**2 + mur2
        mst12 = Rmtop*(At - mu/tanb)
        msb12 = Rmbot*(Ab - mu*tanb)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     
C
C            Stop Eigenvalues calculation
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                                                      
 
        Stop12 = 0.5*(mst11+mst22) +
     * 0.5*((mst11+mst22)**2 - 
     * 4.*(mst11*mst22 - mst12**2))**.5     
        Stop22 = 0.5*(mst11+mst22) -
     * 0.5*((mst11+mst22)**2 - 4.*(mst11*mst22 - mst12**2))**.5     
        if(Stop22.lt.0.) goto 3333
        sstop2(1) = stop12
        sstop2(2) = stop22
        stop1 = Stop12**.5
        stop2 = Stop22**.5 
	stop1w = stop1
	stop2w = stop2
c        if(stop2.lt.45.) goto 3333
c       write(6,*) stop22,stop2

	if(mst12.eq.0.) xst11 = 1.
	if(mst12.eq.0.) xst12 = 0.
	if(mst12.eq.0.) xst21 = 0.
	if(mst12.eq.0.) xst22 = 1.

	if(mst12.eq.0.) goto 223


 22     xst11 = mst12/(mst12**2+(mst11-stop12)**2)**.5
        xst12 = - (mst11-stop12)/(mst12**2+(mst11-stop12)**2)**.5
        xst21 = mst12/(mst12**2+(mst11-stop22)**2)**.5
        xst22 = - (mst11-stop22)/(mst12**2+(mst11-stop22)**2)**.5



 223	T(1,1) = xst11
        T(2,2) = xst22
        T(1,2) = xst12
        T(2,1) = xst21
c	write(*,*)'T=',T
 
c       T(1,1) = xst11
c       T(2,2) = xst22
c      T(1,2) = xst12
c     T(2,1) = xst12



        Sbot12 = 0.5*(msb11+msb22) +
     * 0.5*((msb11+msb22)**2 - 
     * 4.*(msb11*msb22 - msb12**2))**.5     
        Sbot22 = 0.5*(msb11+msb22) -
     * 0.5*((msb11+msb22)**2 - 4.*(msb11*msb22 - msb12**2))**.5     
        if(Sbot22.lt.0.) goto 3333
        sbot1 = Sbot12**.5
        sbot2 = Sbot22**.5
c       if(sbot2.lt.45.) goto 3333
        ssbot2(1) = sbot12
        ssbot2(2) = sbot22

	if(msb12.eq.0.) xsb11 = 1.
	if(msb12.eq.0.) xsb12 = 0.
	if(msb12.eq.0.) xsb21 = 0.
	if(msb12.eq.0.) xsb22 = 1.

	if(msb12.eq.0.) goto 2233





 23     xsb11 = msb12/(msb12**2+(msb11-sbot12)**2)**.5
        xsb12 = - (msb11-sbot12)/(msb12**2+(msb11-sbot12)**2)**.5
        xsb21 = msb12/(msb12**2+(msb11-sbot22)**2)**.5
        xsb22 = - (msb11-sbot22)/(msb12**2+(msb11-sbot22)**2)**.5
 
 2233	B(1,1) = xsb11
        B(2,2) = xsb22
        B(1,2) = xsb12
        B(2,1) = xsb21


        sint = 0.2320
        sqr = 2.**.5
        vp = 174.1*sqr

cccccccccccccccccccccccccccccccccccc
ccc    starting of light higgs
cccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccc      
	if(ihiggs.eq.0)goto 3524
ccccccccccccccccccccccccccccccccccccccccc
        do 4646 i = 1,2
        do 4576 j = 1,2
        coupt(i,j) = 
     * sint*mz**2*2.*sqr/174.1/3.*sinbpa*(delta(i,j) +
     * (3. - 8.*sint)/4./sint*T(1,i)*T(1,j)) 
     * -rmtop**2/174.1**2*vp/sinb*ca*delta(i,j)
     * -rmtop/vp/sinb*(At*ca + mu*sa)*(T(1,i)*T(2,j) +
     *T(1,j)*T(2,i))
 4576   continue
 4646   continue


        do 1646 i = 1,2
        do 1576 j = 1,2
        coupb(i,j) = 
     * -sint*mz**2*2.*sqr/174.1/6.*sinbpa*(delta(i,j) +
     * (3. - 4.*sint)/2./sint*B(1,i)*B(1,j)) 
     * +rmbot**2/174.1**2*vp/cosb*sa*delta(i,j)
     * +rmbot/vp/cosb*(Ab*sa + mu*ca)*(B(1,i)*B(2,j) +
     * B(1,j)*B(2,i))
 1576   continue
 1646   continue


      prun = mh                                                                 
      epss = 1D-4*prun                                                           
      iter = 0                                                                  
7007  iter = iter + 1                                                           
C     WRITE(6,*) 'iter = ',iter                                                 
      DO 7980 i3 = 1,3                                                          
        pr(i3)=prun+(i3-2)*epss/2                                                
        p2=pr(i3)**2                                                            
        polt = 0.                                                               
        DO 7979 i = 1,2                                                         
        DO 7978 j = 1,2                                                         
         polt = polt + coupt(i,j)**2*3.*                                        
     *   fint2(p2,sstop2(i),sstop2(j))/16./pi**2                                 
 7978   CONTINUE                                                                
 7979   CONTINUE                                                                
        polb = 0.                                                               


        DO 9979 i = 1,2                                                         
        DO 9978 j = 1,2                                                         
          polb = polb + coupb(i,j)**2*3.*                                       
     *    fint2(p2,ssbot2(i),ssbot2(j))/16./pi**2                                
 9978   CONTINUE                                                                
 9979   CONTINUE                                                                


        rmtop2 = rmtop**2                                                       
        mtop2=mtop**2                                                           
C                                                                               
        poltt =                                                                 
Cpaj * 3.*rmtop**2/8./pi**2/174.1**2*                                           
     * 3.*rmtop**2/8./pi**2/  v  **2*                                           
     * ca**2/sinb**2 *                                                          
     *   (-2.*mtop**2+.5*p2)*                                                   
     *  fint2(p2,mtop2,mtop2)                                                    
C                                                                               
        pol = polt + polb + poltt                                               
        polar(i3) = p2 - mh**2 - pol                                            
 7980 CONTINUE                                                                  


C     WRITE(6,*) 'Pr    = ',pr                                                  
C     WRITE(6,*) 'Polar = ',polar                                               
      deriv = (polar(3)-polar(1))/epss                                           
      drun = - polar(2)/deriv                                                   
C     WRITE(6,*) 'drun  = ',drun                                                
      prun = prun + drun                                                        
      p2 = prun**2                                                              
      IF ( ABS(drun) .LT. 1D-4 ) GOTO 7777                                      
      GOTO 7007                                                                 
 7777 CONTINUE                                                                  
c	write(*,*) ' HERE'

 
         mhp = p2**.5

cccccccccccccccccccccccccccccccccccccccc
ccc   end of light higgs
cccccccccccccccccccccccccccccccccccccccc      
 3340	 if(ihiggs.eq.1)goto 3524
ccccccccccccccccccccccccccccccccccccccccc
ccc starting of heavy higgs
cccccccccccccccccccccccccccccccccccccccccc



        do 1446 i = 1,2
        do 1476 j = 1,2
        hcoupt(i,j) = 
     * -sint*mz**2*2.*sqr/174.1/3.*cosbpa*(delta(i,j) +
     * (3. - 8.*sint)/4./sint*T(1,i)*T(1,j)) 
     * -rmtop**2/174.1**2*vp/sinb*sa*delta(i,j)
     * -rmtop/vp/sinb*(At*sa - mu*ca)*(T(1,i)*T(2,j) +
     *T(1,j)*T(2,i))
 1476   continue
 1446   continue
 
       do 1146 i = 1,2
        do 1176 j = 1,2
        hcoupb(i,j) = 
     * sint*mz**2*2.*sqr/174.1/6.*cosbpa*(delta(i,j) +
     * (3. - 4.*sint)/2./sint*B(1,i)*B(1,j)) 
     * -rmbot**2/174.1**2*vp/cosb*ca*delta(i,j)
     * -rmbot/vp/cosb*(Ab*ca - mu*sa)*(B(1,i)*B(2,j) +
     * B(1,j)*B(2,i))
	hcoupb(i,j)=0.
 1176   continue
 1146   continue

      prun = hm                                                                 
      epss = 1D-4*prun                                                           
      iter = 0                                                                  
 1001 iter = iter + 1                                                           
C     WRITE(6,*) 'iter = ',iter                                                 
      DO 1780 i3 = 1,3                                                          
        pr(i3)=prun+(i3-2)*epss/2                                                
        hp2=pr(i3)**2                                                           
C                                                                               
        hpolt = 0.                                                              
        do 1779 i = 1,2                                                         
        do 1778 j = 1,2                                                         
        hpolt = hpolt + hcoupt(i,j)**2*3.*                                      
     *  fint2(hp2,sstop2(i),sstop2(j))/16./pi**2                                 
 1778 CONTINUE                                                                  
 1779 CONTINUE                                                                  
C                                                                               
      hpolb = 0.                                                                
      DO 1979 I = 1,2                                                           
      DO 1978 J = 1,2                                                           
        hpolb = hpolb + hcoupb(i,j)**2*3.*                                      
     *  fint2(hp2,ssbot2(i),ssbot2(j))/16./pi**2                                 
 1978 CONTINUE                                                                  
 1979 CONTINUE                                                                  
C                                                                               
      rmtop2 = rmtop**2                                                         
      mtop2  = mtop**2                                                          
C                                                                               
      hpoltt =                                                                  
Cpaj * 3.*rmtop**2/8./pi**2/174.1**2*                                           
     * 3.*rmtop**2/8./pi**2/  v  **2*                                           
     *  sa**2/sinb**2 *                                                         
     *   (-2.*mtop**2+.5*hp2)*                                                  
     *  fint2(hp2,mtop2,mtop2)                                                   
C                                                                               
      hpol = hpolt + hpolb + hpoltt                                             
      polar(i3) =hp2-hm**2-hpol                                                 
 1780 CONTINUE                                                                  
C     WRITE(6,*) 'Pr    = ',pr                                                  
C     WRITE(6,*) 'Polar = ',polar                                               
      deriv = (polar(3)-polar(1))/epss                                           
      drun = - polar(2)/deriv                                                   
C     WRITE(6,*) 'drun  = ',drun                                                
      prun = prun + drun                                                        
      hp2 = prun**2                                                             
      IF ( ABS(drun) .LT. 1D-4 ) GOTO 1111                                      
      GOTO 1001                                                                 
 1111 CONTINUE                                                                  


 2222	continue
           hmp = hp2**.5
ccccccccccccccccccccccccccccccccccccccccccc
ccc  end of heavy higgs
cccccccccccccccccccccccccccccccccccccccccccc
	if(ihiggs.eq.2)goto 3524
cccccccccccccccccccccccccccccccccccccccccccc
ccc  beginning of pseudoscalar higgs
cccccccccccccccccccccccccccccccccccccccccccc


        do 3446 i = 1,2
        do 3476 j = 1,2
        acoupt(i,j) = 
     * -rmtop/vp/sinb*(At*cosb + mu*sinb)*
     *  (T(1,i)*T(2,j) -T(1,j)*T(2,i))
 3476   continue
 3446   continue
        do 3146 i = 1,2
        do 3176 j = 1,2
        acoupb(i,j) = 
     * rmbot/vp/cosb*(Ab*sinb + mu*cosb)*
     *  (B(1,i)*B(2,j) -B(1,j)*B(2,i))
 3176   continue
 3146   continue

      prun = ma                                                                 
      epss = 1D-4*prun                                                           
      iter = 0                                                                  
 6006 iter = iter + 1                                                           
C     WRITE(6,*) 'iter = ',iter                                                 
      DO 3780 i3 = 1,3                                                          
        pr(i3)=prun+(i3-2)*epss/2                                                
        ap2=pr(i3)**2                                                           
        apolt = 0.                                                              
        DO 3779 I = 1,2                                                         
        DO 3778 J = 1,2                                                         
          apolt = apolt + acoupt(i,j)**2*3.*                                    
     *    fint2(ap2,sstop2(i),sstop2(j))/16./pi**2                               
 3778   CONTINUE                                                                
 3779   CONTINUE                                                                
        apolb = 0.                                                              
        DO 3979 I = 1,2                                                         
        DO 3978 J = 1,2                                                         
          apolb = apolb + acoupb(i,j)**2*3.*                                    
     *    fint2(ap2,ssbot2(i),ssbot2(j))/16./pi**2                               
 3978   CONTINUE                                                                
 3979   CONTINUE                                                                
        rmtop2 = rmtop**2                                                       
        mtop2=mtop**2                                                           
        apoltt =                                                                
Cpaj *  3.*rmtop**2/8./pi**2/174.1**2*                                          
     *  3.*rmtop**2/8./pi**2/  v  **2*                                          
     *  cosb**2/sinb**2 *                                                       
     *   (-.5*ap2)*                                                             
     *  fint2(ap2,mtop2,mtop2)                                                   
        apol = apolt + apolb + apoltt                                           
        polar(i3) = ap2 - ma**2 -apol                                           
 3780 CONTINUE                                                                  
C     WRITE(6,*) 'Pr    = ',pr                                                  
C     WRITE(6,*) 'Polar = ',polar                                               
      deriv = (polar(3)-polar(1))/epss                                           
      drun = - polar(2)/deriv                                                   
C     WRITE(6,*) 'drun  = ',drun                                                
      prun = prun + drun                                                        
      ap2 = prun**2                                                             
      IF ( ABS(drun) .LT. 1D-4 ) GOTO 6666                                      
      GOTO 6006                                                                 
 6666 CONTINUE                                                                  




        amp = ap2**.5

ccccccccccccccccccccccccccccccccccccccccccc
ccc end of pseudoscalar higgs
cccccccccccccccccccccccccccccccccccccccccccc
	if(ihiggs.eq.3)goto 3524
cccccccccccccccccccccccccccccccccccccccccccc 

3524	return
 	
      
3333   stop 
       end

       FUNCTION Z(x)
       implicit real*8(a-h,m,o-z)
       SA = 1. - 4.*X
c       write(6,*) x,sa
       if(Sa.lt.0.) sa1 = abs(sa)
       if(SA.lt.0.) Z = 2.*SA1**.5*atan(1./SA1**.5)
       if(SA.gt.0.) Z = SA**.5*log((1.+SA**.5)/(1.-SA**.5))
c       write(6,*) Sa,Z
       return
       end

      FUNCTION FINTAN(p2,y1,y2)
      implicit real*8(a-h,m,o-z)
	delta=(y1-y2)/p2
	erre=(abs((1.+delta)**2-4.*y1/p2))**.5
	fintan=-1.+.5*((y1+y2)/(y1-y2)-delta)*log(y2/y1)	
     *    +.5*erre*log(abs((delta**2-(1.+erre)**2)) / 
     *      abs((delta**2-(1.-erre)**2)) )
      return
      end




      subroutine rghm(mchi,ma,tanb,mq,mur,md,mtop,au,ad,mu,
     *    mhp,hmp,sa,ca,sab,cab,tanbA,mglu,deltamt,deltamb)
      implicit real*8(a-h,l,m,o-z)
      dimension vh(2,2),m2(2,2),m2p(2,2)
c      read(5,*) ma,tanb,mq,mur,md,mtop,Au,Ad,mu

      mz = 91.18
      alpha1 = 0.0101
      alpha2 = 0.0337
      alpha3Z = 0.12    
      v = 174.1
      pi = 3.14159
      tanbA = tanb
      tanbt = tanb	

C     mbottom(mtop) = 3. GeV
      mb = 3.
      alpha3 = alpha3Z/(1. +(11. - 10./3.)/4./pi*alpha3Z*
     *log(mtop**2/mz**2))

C     rmtop= running top quark mass     
      rmtop = mtop/(1.+4.*alpha3/3./pi)
      tq = log((mq**2+mtop**2)/mtop**2)
      tu = log((mur**2 + mtop**2)/mtop**2)
      td = log((md**2 + mtop**2)/mtop**2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    New definition, tglu.
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
      tglu = log(mglu**2/mtop**2)
      sinb = tanb/((1. + tanb**2)**.5)
      cosb = sinb/tanb
      if(ma.gt.mtop)
     *tanbA = tanb*(1.-3./32./pi**2*
     *(rmtop**2/v**2/sinb**2-mb**2/v**2/cosb**2)*
     *log(ma**2/mtop**2))
      if(ma.lt.mtop.or.ma.eq.mtop) tanbt = tanbA
      sinb = tanbt/((1. + tanbt**2)**.5)
      cosb = 1./((1. + tanbt**2)**.5)
      cos2b = (tanbt**2 - 1.)/(tanbt**2 + 1.)
      g1 = (alpha1*4.*pi)**.5 
      g2 = (alpha2*4.*pi)**.5 
      g3 = (alpha3*4.*pi)**.5 
      hu = rmtop/v/sinb
      hd =  mb/v/cosb

      call Gfun2(ma,tanbA,mq,mur,md,mtop,Au,Ad,mu,mglu,vh,stop1,stop2,
     *sbot1,sbot2,deltamt,deltamb)

      if(mq.gt.mur) tp = tq - tu
      if(mq.lt.mur.or.mq.eq.mur) tp = tu - tq
      if(mq.gt.mur) tdp = tu
      if(mq.lt.mur.or.mq.eq.mur) tdp = tq
      if(mq.gt.md) tpd = tq - td
      if(mq.lt.md.or.mq.eq.md) tpd = td - tq
      if(mq.gt.md) tdpd = td
      if(mq.lt.md.or.mq.eq.md) tdpd = tq

      if(mq.gt.md) dlambda1 = 6./96./pi**2*g1**2*hd**2*tpd
      if(mq.lt.md.or.mq.eq.md) dlambda1 = 3./32./pi**2*
     * hd**2*(g1**2/3.+g2**2)*tpd

      if(mq.gt.mur) dlambda2 =12./96./pi**2*g1**2*hu**2*tp	
      if(mq.lt.mur.or.mq.eq.mur) dlambda2 = 3./32./pi**2*
     * hu**2*(-g1**2/3.+g2**2)*tp

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  dlambdap1 and dlambdap2 are the new log corrections due to
c  the presence of the gluino mass. They are in general very small, 
c  and only present if there is a hierarchy of masses between the
c  two stops.
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    

	if(mglu.lt.mur.or.mglu.lt.mq) then
	if(mq.gt.mur.and.mglu.gt.mur) then
	dlambdap2 = -4./(16.*pi**2)**2*hu**4*(tq**2-tglu**2)
	endif

	if(mq.gt.mur.and.mglu.lt.mur) then
	dlambdap2 = -4./(16.*pi**2)**2*hu**4*(tq**2-tu**2)
	endif

	if(mq.gt.mur.and.mglu.eq.mur) then
	dlambdap2 = -4./(16.*pi**2)**2*hu**4*(tq**2-tu**2)
	endif

	if(mur.gt.mq.and.mglu.gt.mq) then
	dlambdap2 = -4./(16.*pi**2)**2*hu**4*(tu**2-tglu**2)
	endif

	if(mur.gt.mq.and.mglu.lt.mq) then
	dlambdap2 = -4./(16.*pi**2)**2*hu**4*(tu**2-tq**2)
	endif

	if(mur.gt.mq.and.mglu.eq.mq) then
	dlambdap2 = -4./(16.*pi**2)**2*hu**4*(tu**2-tq**2)
	endif
	endif

	


      dlambda3 = 0.
      dlambda4 = 0.

      if(mq.gt.md) dlambda3 = -1./32./pi**2*g1**2*hd**2*tpd
      if(mq.lt.md.or.mq.eq.md) dlambda3 = 3./64./pi**2*hd**2*
     *(g2**2-g1**2/3.)*tpd
      

      if(mq.gt.mur) dlambda3 = dlambda3 - 
     *1./16./pi**2*g1**2*hu**2*tp
      if(mq.lt.mur.or.mq.eq.mur) dlambda3 = dlambda3 + 
     * 3./64./pi**2*hu**2*(g2**2+g1**2/3.)*tp

      if(mq.lt.mur) dlambda4 = -3./32./pi**2*g2**2*hu**2*tp
      if(mq.lt.md) dlambda4 = dlambda4 - 3./32./pi**2*g2**2*
     *hd**2*tpd






      lambda1 = ((g1**2 + g2**2)/4.)*
     * (1.-3.*hd**2*(tpd + tdpd)/8./pi**2)
     *+(3.*hd**4./16./pi**2) *tpd*(1.   
     *+ (3.*hd**2/2. + hu**2/2.       
     *- 8.*g3**2) * (tpd + 2.*tdpd)/16./pi**2) 
     *+(3.*hd**4./8./pi**2) *tdpd*(1.  + (3.*hd**2/2. + hu**2/2.       
     *- 8.*g3**2) * tdpd/16./pi**2) + dlambda1 
      lambda2 = ((g1**2 + g2**2)/4.)*(1.-3.*hu**2*
     *(tp + tdp)/8./pi**2)
     *+(3.*hu**4./16./pi**2) *tp*(1.   
     *+ (3.*hu**2/2. + hd**2/2.       
     *- 8.*g3**2) * (tp + 2.*tdp)/16./pi**2) 
     *+(3.*hu**4./8./pi**2) *tdp*(1. + (3.*hu**2/2. + hd**2/2.       
     *- 8.*g3**2) * tdp/16./pi**2) + dlambda2 + dlambdap2
      lambda3 = ((g2**2 - g1**2)/4.)*(1.-3.*
     *(hu**2)*(tp + tdp)/16./pi**2 -3.*
     *(hd**2)*(tpd + tdpd)/16./pi**2) +dlambda3 
      lambda4 = (- g2**2/2.)*(1.
     *-3.*(hu**2)*(tp + tdp)/16./pi**2
     *-3.*(hd**2)*(tpd + tdpd)/16./pi**2) +dlambda4
      
	lambda5 = 0.
	lambda6 = 0.
	lambda7 = 0.

      m2(1,1) = 2.*v**2*(lambda1*cosb**2+2.*lambda6*
     *cosb*sinb + lambda5*sinb**2) + ma**2*sinb**2

      m2(2,2) = 2.*v**2*(lambda5*cosb**2+2.*lambda7*
     *cosb*sinb + lambda2*sinb**2) + ma**2*cosb**2
      m2(1,2) = 2.*v**2*(lambda6*cosb**2+(lambda3+lambda4)*
     *cosb*sinb + lambda7*sinb**2) - ma**2*sinb*cosb

      m2(2,1) = m2(1,2)
ccccccccccccccccccccccccccccccccccccccccccccccccc
ccc  this is the contribution from light charginos/neutralinos
ccccccccccccccccccccccccccccccccccccccccccccccccc
	
	mssusy=(.5*(mq**2+mur**2)+mtop**2)**.5

	if(mchi.gt.mssusy)goto 3790
	if(mchi.lt.mtop) mchi=mtop

	tchar=log(mssusy**2/mchi**2)

	deltal12=(9./64./pi**2*g2**4+5./192./pi**2*g1**4)*tchar
	deltal3p4=(3./64./pi**2*g2**4+7./192./pi**2*g1**4
     *       +4./32/pi**2*g1**2*g2**2)*tchar

	deltam112=2.*deltal12*v**2*cosb**2
	deltam222=2.*deltal12*v**2*sinb**2
	deltam122=2.*deltal3p4*v**2*sinb*cosb

	m2(1,1)=m2(1,1)+deltam112
	m2(2,2)=m2(2,2)+deltam222
	m2(1,2)=m2(1,2)+deltam122
	m2(2,1)=m2(2,1)+deltam122
		
 3790	continue

ccccccccccccccccccccccccccccccccccccccccccc
ccc  end of charginos/neutralinos
ccccccccccccccccccccccccccccccccccccccccccc

      do 9800 i = 1,2
      do 9801 j = 1,2
      m2p(i,j) = m2(i,j) + vh(i,j)
c      write(886,*) m2p(i,j),m2(i,j)
 9801 continue
 9800 continue

      Trm2p = m2p(1,1) + m2p(2,2)
      detm2p = m2p(1,1)*m2p(2,2) - m2p(1,2)*m2p(2,1)


c      write(6,*) Trm2,Trm2p,detm2,detm2p

      mh2p = (Trm2p - (Trm2p**2 - 4.* detm2p)**.5)/2.
      HM2p = (Trm2p + (Trm2p**2 - 4.* detm2p)**.5)/2.
      HMp = Hm2p**.5 
c      write(886,*) mh2
c      if(mh2.lt.0.) goto 5555
      if(mh2p.lt.0.) goto 5555
c      mh = mh2**.5 
      mhp = mh2p**.5
      sin2alpha = 2.*m2p(1,2)/(Trm2p**2-4.*detm2p)**.5
      cos2alpha = (m2p(1,1)-m2p(2,2))/(Trm2p**2-4.*detm2p)**.5
      if(cos2alpha.gt.0.) alpha = asin(sin2alpha)/2.
      if(cos2alpha.lt.0.) alpha = -pi/2.-asin(sin2alpha)/2.
      sa = sin(alpha)
      ca = cos(alpha)

      write(*,*) 'deltamb = ',deltamb,'deltamt= ',deltamt 

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        Here the values of sab and cab are defined, in order
c        to define the new couplings of the lightest and 
c        heavy CP-even Higgs to the bottom quark.
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	

      sab = sa*(1.-deltamb/(1.+deltamb)*(1.+ca/sa/tanb))
      cab = ca*(1.-deltamb/(1.+deltamb)*(1.-sa/ca/tanb)) 	
      sqbma = (sinb*ca - cosb*sa)**2
c     write(6,*) mhp,Hmp,sa,ca,tanb,stop1,stop2,sqbma
 5555  xin = 1.
 2242 return
      end





        SUBROUTINE GFUN2(ma,tanb,mq,mur,md,mtop,At,Ab,mu,mglu,vh,
     *  stop1,stop2,sbot1,sbot2,deltamt,deltamb)
        implicit real*8 (a-h,l,m,o-z)
        dimension diah(2),vh(2,2),vh1(2,2),vh2(2,2),
     *  vh3t(2,2),vh3b(2,2),
     *  hmix(2,2),al(2,2),m2(2,2)

        if(dabs(mu).lt.0.000001) mu = 0.000001
        mq2 = mq**2
        mur2 = mur**2
        md2 = md**2
        tanbA = tanb
        sinbA = tanbA/(tanbA**2+1.)**.5
        cosbA = sinbA/tanbA        

        sinb = tanb/(tanb**2+1.)**.5
        cosb = sinb/tanb
        pi = 3.14159
        g2 = (0.0336*4.*pi)**.5
        g12 = (0.0101*4.*pi)
        g1 = g12**.5
        mz = 91.18
        v = 174.1
        mw = (g2**2*v**2/2.)**.5
        alpha3 = 0.12/(1.+23/12./pi*0.12*log(mtop**2/mz**2))

        mb = 3.
        if(mq.gt.mur) mst = mq
        if(mur.gt.mq.or.mur.eq.mq) mst = mur

        msusyt = (mst**2  + mtop**2)**.5

	if(mq.gt.md) msb = mq
	if(md.gt.mq.or.md.eq.mq) msb = md
	
	msusyb = (msb**2 + mb**2)**.5
c	if(mq.gt.md.and.mq.gt.mur) msusyb = msusyt

	tt = log(msusyt**2/mtop**2)
	tb = log(msusyb**2/mtop**2)

 645    format(2x,'mu=',2x,f13.4)
c        if(mq.gt.mur) t = log((mq2 + mtop**2)/mtop**2)
c        if(mur.gt.mq.or.mq.eq.mur) t = log((mur2+mtop**2)/mtop**2)
        rmtop = mtop/(1.+4.*alpha3/3./pi)
        ht = rmtop/(174.1*sinb)
        htst = rmtop/174.1
        hb = mb/174.1/cosb
        g32 = alpha3*4.*pi
        bt2 = -(8.*g32 - 9.*ht**2/2. - hb**2/2.)/(4.*pi)**2
	bb2 = -(8.*g32 - 9.*hb**2/2. - ht**2/2.)/(4.*pi)**2
        al2 = 3./8./pi**2*ht**2
        bt2st = -(8.*g32 - 9.*htst**2/2.)/(4.*pi)**2
        alst = 3./8./pi**2*htst**2
        al1 = 3./8./pi**2*hb**2

        al(1,1) = al1
        al(1,2) = (al2+al1)/2.
        al(2,1) = (al2+al1)/2.
        al(2,2) = al2

	


c	write(6,*) mbot4,mbot2,mbot2**.5

	if(ma.gt.mtop) then
        vi = 174.1*(1. + 3./32./pi**2*htst**2*
     *  log(mtop**2/ma**2))
        h1i = vi* cosbA
        h2i = vi*sinbA
        h1t = h1i*(1.+3./8./pi**2*hb**2*log(ma**2/msusyt**2))**.25
        h2t = h2i*(1.+3./8./pi**2*ht**2*log(ma**2/msusyt**2))**.25
        h1b = h1i*(1.+3./8./pi**2*hb**2*log(ma**2/msusyb**2))**.25
        h2b = h2i*(1.+3./8./pi**2*ht**2*log(ma**2/msusyb**2))**.25
	else
	vi = 174.1
	h1i = vi*cosb
	h2i = vi*sinb
        h1t = h1i*(1.+3./8./pi**2*hb**2*log(mtop**2/msusyt**2))**.25
        h2t = h2i*(1.+3./8./pi**2*ht**2*log(mtop**2/msusyt**2))**.25
        h1b = h1i*(1.+3./8./pi**2*hb**2*log(mtop**2/msusyb**2))**.25
        h2b = h2i*(1.+3./8./pi**2*ht**2*log(mtop**2/msusyb**2))**.25
	end if





        tanbst = h2t/h1t
        sinbt = tanbst/(1.+tanbst**2)**.5
        cosbt = sinbt/tanbst


        tanbsb = h2b/h1b
        sinbb = tanbsb/(1.+tanbsb**2)**.5
        cosbb = sinbb/tanbsb

	deltamt = 0.
	deltamb = 0.

        mtop4 = rmtop**4.*(1.+2.*bt2*tt- al2*tt - 4.*deltamt)
c     * /(1.+deltamt)**4.
        mtop2 = mtop4**.5
	mbot4 = mb**4.*(1.+2.*bb2*tb - al1*tb)
     * /(1.+deltamb)**4.
	mbot2 = mbot4**.5


        stop12 = (mq2 + mur2)*.5 + mtop2 
     *   +1./8.*(g2**2+g1**2)*(h1t**2-h2t**2)
     *   +(((g2**2-5.*g1**2/3.)/4.*(h1t**2-h2t**2) +
     *   mq2 - mur2)**2*0.25 + mtop2*(At-mu/tanbst)**2)**.5
        stop22 = (mq2 + mur2)*.5 + mtop2 
     *  +1./8.*(g2**2+g1**2)*(h1t**2-h2t**2) 
     *   - (((g2**2-5.*g1**2/3.)/4.*(h1t**2-h2t**2) +
     *  mq2 - mur2)**2*0.25 
     *  + mtop2*(At-mu/tanbst)**2)**.5
c	write(6,*) stop22**.5, mtop2**.5, rmtop
        if(stop22.lt.0.) goto 4237
        sbot12 = (mq2 + md2)*.5  
     *   - 1./8.*(g2**2+g1**2)*(h1b**2-h2b**2)
     *  + (((g1**2/3.-g2**2)/4.*(h1b**2-h2b**2) +
     *  mq2 - md2)**2*0.25 + mbot2*(Ab-mu*tanbsb)**2)**.5
        sbot22 = (mq2 + md2)*.5  
     *   - 1./8.*(g2**2+g1**2)*(h1b**2-h2b**2)
     *   - (((g1**2/3.-g2**2)/4.*(h1b**2-h2b**2) +
     *   mq2 - md2)**2*0.25 + mbot2*(Ab-mu*tanbsb)**2)**.5
        if(sbot22.lt.0.) sbot22 = 10000.


        stop1 = stop12**.5
        stop2 = stop22**.5
        sbot1 = sbot12**.5
        sbot2 = sbot22**.5
        write(*,*) 'TEST',At,mtop2,stop1,stop2
c	write(6,*) sbot1,sbot2,stop1,stop2

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Here is the definition of deltamb and deltamt, which
c     are the vertex corrections to the bottom and top quark
c     mass, keeping the dominant QCD and top Yukawa coupling
c     induced corrections.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	deltamb = -2*alpha3/3./pi*mglu*(ab-mu*tanb)*
     *  T(sbot1,sbot2,mglu)
     *  + ht**2/(4.*pi)**2*(at-mu/tanb)*mu*tanb*
     *  T(stop1,stop2,mu)

  
	deltamt = -2.*alpha3/3./pi*(at-mu/tanb)*mglu*
     *  T(stop1,stop2,mglu)

c	write(*,*) deltamb,deltamt,mglu,T(stop1,stop2,mglu),tanb,
c     *  alpha3,at,mu,stop1,stop2,T(sbot1,sbot2,mglu)


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Here the new values of the top and bottom quark masses at
c   the scale MS are defined, to be used in the effective 
c   potential approximation. They are just the old ones, but
c   including the finite corrections deltamt and deltamb.
c   The deltamb corrections can become large and are resummed
c   to all orders, as suggested in the two recent works by M. Carena, 
c   S. Mrenna and C.E.M. Wagner, as well as in the work by M. Carena,
c   D. Garcia, U. Nierste and C.E.M. Wagner, to appear. The top
c   quark mass corrections are small and are kept in the perturbative
c   formulation.  The function T(X,Y,Z) is necessary for the calculation.
c   the entries are masses and NOT their squares !
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        mtop4 = rmtop**4.*(1.+2.*bt2*tt- al2*tt - 4.*deltamt)
c     * /(1.+deltamt)**4.
        mtop2 = mtop4**.5
	mbot4 = mb**4.*(1.+2.*bb2*tb - al1*tb)
     * /(1.+deltamb)**4.
	mbot2 = mbot4**.5


        stop12 = (mq2 + mur2)*.5 + mtop2 
     *   +1./8.*(g2**2+g1**2)*(h1t**2-h2t**2)
     *   +(((g2**2-5.*g1**2/3.)/4.*(h1t**2-h2t**2) +
     *   mq2 - mur2)**2*0.25 + mtop2*(At-mu/tanbst)**2)**.5
        stop22 = (mq2 + mur2)*.5 + mtop2 
     *  +1./8.*(g2**2+g1**2)*(h1t**2-h2t**2) 
     *   - (((g2**2-5.*g1**2/3.)/4.*(h1t**2-h2t**2) +
     *  mq2 - mur2)**2*0.25 
     *  + mtop2*(At-mu/tanbst)**2)**.5
c	write(6,*) stop22**.5, mtop2**.5, rmtop
        if(stop22.lt.0.) goto 4237
        sbot12 = (mq2 + md2)*.5  
     *   - 1./8.*(g2**2+g1**2)*(h1b**2-h2b**2)
     *  + (((g1**2/3.-g2**2)/4.*(h1b**2-h2b**2) +
     *  mq2 - md2)**2*0.25 + mbot2*(Ab-mu*tanbsb)**2)**.5
        sbot22 = (mq2 + md2)*.5  
     *   - 1./8.*(g2**2+g1**2)*(h1b**2-h2b**2)
     *   - (((g1**2/3.-g2**2)/4.*(h1b**2-h2b**2) +
     *   mq2 - md2)**2*0.25 + mbot2*(Ab-mu*tanbsb)**2)**.5
        if(sbot22.lt.0.) goto 4237


        stop1 = stop12**.5
        stop2 = stop22**.5
        sbot1 = sbot12**.5
        sbot2 = sbot22**.5




        vh1(1,1) = 1./tanbst
        vh1(2,1) = -1.
        vh1(1,2) = -1.
        vh1(2,2) = tanbst
        vh2(1,1) = tanbst
        vh2(1,2) = -1.
        vh2(2,1) = -1.
        vh2(2,2) = 1./tanbst
cccccccccccccccccccccccccccccccc
ccc   D-terms
cccccccccccccccccccccccccccccccc
	stw=.2320

	f1t=(mq2-mur2)/(stop12-stop22)*(.5-4./3.*stw)*
     *         log(stop1/stop2)
     *        +(.5-2./3.*stw)*log(stop1*stop2/(mq2+mtop2))
     *        + 2./3.*stw*log(stop1*stop2/(mur2+mtop2))

	f1b=(mq2-md2)/(sbot12-sbot22)*(-.5+2./3.*stw)*
     *        log(sbot1/sbot2)
     *        +(-.5+1./3.*stw)*log(sbot1*sbot2/(mq2+mbot2))
     *        - 1./3.*stw*log(sbot1*sbot2/(md2+mbot2))

	f2t=mtop2**.5*(at-mu/tanbst)/(stop12-stop22)*
     *         (-.5*log(stop12/stop22)
     *        +(4./3.*stw-.5)*(mq2-mur2)/(stop12-stop22)*
     *         g(stop12,stop22))

	f2b=mbot2**.5*(ab-mu*tanbsb)/(sbot12-sbot22)*
     *         (.5*log(sbot12/sbot22)
     *        +(-2./3.*stw+.5)*(mq2-md2)/(sbot12-sbot22)*
     *        g(sbot12,sbot22))

	

        vh3b(1,1) = mbot4/(cosbb**2)*(log(sbot1**2*sbot2**2/
     *  (mq2+mbot2)/(md2+mbot2)) 
     *  + 2.*(Ab*(Ab-mu*tanbsb)/(sbot1**2-sbot2**2))*
     *  log(sbot1**2/sbot2**2)) + 
     *  Mbot4/(cosbb**2)*(Ab*(Ab-mu*tanbsb)/
     *  (sbot1**2-sbot2**2))**2*g(sbot12,sbot22) 

	vh3t(1,1) =
     *  mtop4/(sinbt**2)*(mu*(-At+mu/tanbst)/(stop1**2
     * -stop2**2))**2*g(stop12,stop22)

	vh3b(1,1)=vh3b(1,1)+
     *    mz**2*(2*mbot2*f1b-mbot2**.5*ab*f2b)

	vh3t(1,1) = vh3t(1,1) + 
     *  mz**2*(mtop2**.5*mu/tanbst*f2t)  

        vh3t(2,2) = mtop4/(sinbt**2)*(log(stop1**2*stop2**2/
     *  (mq2+mtop2)/(mur2+mtop2)) 
     *  + 2.*(At*(At-mu/tanbst)/(stop1**2-stop2**2))*
     *  log(stop1**2/stop2**2)) + 
     *  mtop4/(sinbt**2)*(At*(At-mu/tanbst)/
     *  (stop1**2-stop2**2))**2*g(stop12,stop22) 

	vh3b(2,2) =
     *  Mbot4/(cosbb**2)*(mu*(-Ab+mu*tanbsb)/(sbot1**2
     * -sbot2**2))**2*g(sbot12,sbot22)




	vh3t(2,2)=vh3t(2,2)+
     *    mz**2*(-2*mtop2*f1t+mtop2**.5*at*f2t)

	vh3b(2,2) = vh3b(2,2) -mz**2*mbot2**.5*mu*tanbsb*f2b  
 

        vh3t(1,2) = -
     *   mtop4/(sinbt**2)*mu*(At-mu/tanbst)/
     * (stop1**2-stop2**2)*(log(stop1**2/stop2**2) + At*
     * (At - mu/tanbst)/(stop1**2-stop2**2)*g(stop12,stop22))

	vh3b(1,2) =
     * - mbot4/(cosbb**2)*mu*(At-mu*tanbsb)/
     * (sbot1**2-sbot2**2)*(log(sbot1**2/sbot2**2) + Ab*
     * (Ab - mu*tanbsb)/(sbot1**2-sbot2**2)*g(sbot12,sbot22))


	vh3t(1,2)=vh3t(1,2) +
     *      mz**2*(mtop2/tanbst*f1t-mtop2**.5*(at/tanbst+mu)/2.*f2t)

	vh3b(1,2)=vh3b(1,2)
     *  +mz**2*(-mbot2*tanbsb*f1b+mbot2**.5*(ab*tanbsb+mu)/2.*f2b)

        vh3t(2,1) = vh3t(1,2)
        vh3b(2,1) = vh3b(1,2)

       tq = log((mq2 + mtop2)/mtop2)
       tu = log((mur2+mtop2)/mtop2)
       tqd = log((mq2 + mb**2)/mb**2)
       td = log((md2+mb**2)/mb**2)


        do 8910 i = 1,2
        do 8911 j = 1,2

        vh(i,j) = 
c     ((
c     *(g2**2+g12)* (h1t**2+h2t**2)*
c     *  (sinbt*cosbt)*vh1(i,j)+
c     *  vh2(i,j)*ma**2*(2.*sinb*cosb) + 
     *  6./(8.*pi**2*(h1t**2+h2t**2))
     *  *vh3t(i,j)*0.5*(1.-al(i,j)*tt/2.) +
     *  6./(8.*pi**2*(h1b**2+h2b**2))
     *  *vh3b(i,j)*0.5*(1.-al(i,j)*tb/2.)

c	write(6,*) vh3t(i,j),vh3b(i,j),vh(i,j),f1b,f2b,f1t,f2t
 8911   continue
 8910   continue

        goto 4236
 4237   do 6868 i =1,2
        do 6867 j = 1,2
        vh(i,j) = -1.d+15
 6867   continue
 6868   continue

        
 4236   return
        end

      FUNCTION g(X,Y)
      implicit real*8(a-h,l,m,o-z)
      g = 2. - (X+Y)/(X-Y)*log(X/Y)
      return
      end


      FUNCTION T(X,Y,Z) 
      implicit real*8(a-h,l,m,o-z)
      if(x.eq.y) x = x - 0.00001
      if(x.eq.z) x = x - 0.00002
      if(y.eq.z) y = y - 0.00003		
c	write(*,*) 'xyz',x,y,z
      T = (X**2*Y**2*log(X**2/Y**2) + X**2*Z**2*log(Z**2/X**2) 
     * + Y**2*Z**2*log(Y**2/Z**2))/((X**2-Y**2)*(Y**2-Z**2)*(X**2-Z**2))
      return
      end	

      DOUBLE PRECISION FUNCTION FINT2(a,b,c)                                     
      IMPLICIT REAL*8(A-H,O-Z)                                                  
      COMMON /CEM/ y1,y2,p2,eps                                                     
      EXTERNAL fintsub,dgauss   
	eps=1.d-02                                                      
      p2 = a                                                                    
      y1 = b                                                                    
      y2 = c                                                                    
      xlo = 0.D0                                                                 
      xhi = 1.D0                                                                 
      fint2  = DGAUSS(fintsub,xlo,xhi,eps)                                     
      RETURN                                                                    
      END                                                                       
                                                                                
      FUNCTION fintsub(x)                                             
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /CEM/ y1,y2,p2,eps                                                     
        fintsub = LOG(ABS(x*y1+(1-x)*y2-x*(1-x)*p2)                   
     .                /(x*(y1-y2)+y2))                                       

      RETURN                                                                    
      END                                                                       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C    DGAUSS stuff
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      FUNCTION DGAUSS(F,A,B,EPS)
C.----------------------------------------------------------------------
C.
C.    GAUSS INTEGRAL OF THE FUNCTION F IN INTERVAL A,B
C.    LAST UPDATE: 12/03/87
C.
C.----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION W(12),X(12)
      EXTERNAL F
      DATA CONST/1.E-12/
      DATA W
     &/0.101228536290376, 0.222381034453374, 0.313706645877887,
     & 0.362683783378362, 0.027152459411754, 0.062253523938648,
     & 0.095158511682493, 0.124628971255534, 0.149595988816577,
     & 0.169156519395003, 0.182603415044924, 0.189450610455069/
      DATA X
     &/0.960289856497536, 0.796666477413627, 0.525532409916329,
     & 0.183434642495650, 0.989400934991650, 0.944575023073233,
     & 0.865631202387832, 0.755404408355003, 0.617876244402644,
     & 0.458016777657227, 0.281603550779259, 0.095012509837637/
C--
C--   INITIALISE
      DELTA=CONST*ABS(A-B)
      DGAUSS=0.
      AA=A
C--
C--   ITERATION LOOP
   10 Y=B-AA
C--
C--   EPSILON REACHED ??
      IF (ABS(Y).LE.DELTA) RETURN
   20 BB=AA+Y
      C1=0.5*(AA+BB)
      C2=C1-AA
      S8=0.
      S16=0.
      DO 30 I=1,4
         U=X(I)*C2
   30 S8=S8+W(I)*(F(C1+U)+F(C1-U))
      DO 40 I=5,12
         U=X(I)*C2
   40 S16=S16+W(I)*(F(C1+U)+F(C1-U))
      S8=S8*C2
      S16=S16*C2
      IF (ABS(S16-S8).GT.EPS*(1.0+ABS(S16))) GOTO 50
      DGAUSS=DGAUSS+S16
      AA=BB
      GOTO 10
   50 Y=0.5*Y
      IF (ABS(Y).GT.DELTA) GOTO 20
C      WRITE (6,9000)
      DGAUSS=0.
      RETURN
 9000 FORMAT(1H ,'****** DGAUSS... TOO HIGH ACCURACY REQUIRED ******')
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     THIS PROGRAM COMPUTES THE RENORMALIZATION GROUP IMPROVED
C     VALUES OF HIGGS MASSES AND COUPLINGS IN THE MSSM.
C
C     INPUT: MA,TANB = TAN(BETA),MQ,MUR,MDR,MTOP,AU,AD,MU,MCHI
C
C     ALL MASSES IN GEV UNITS. MA IS THE CP-ODD HIGGS MASS,
C     MTOP IS THE PHYSICAL TOP MASS, MQ AND MUR/MDR ARE THE SOFT
C     SUPERSYMMETRY BREAKING MASS PARAMETERS OF LEFT HANDED
C     AND RIGHT HANDED STOPS RESPECTIVELY, AU AND AD ARE THE
C     STOP AND SBOTTOM TRILINEAR SOFT BREAKING TERMS,
C     RESPECTIVELY,  AND MU IS THE SUPERSYMMETRIC
C     HIGGS MASS PARAMETER. WE USE THE  CONVENTIONS FROM
C     THE PHYSICS REPORT OF HABER AND KANE: LEFT RIGHT
C     STOP MIXING TERM PROPORTIONAL TO (AU - MU/TANB).
C     MCHI IS THE HEAVIEST CHARGINO MASS. 
C     WE USE AS INPUT TANB DEFINED AT THE SCALE MTOP.

C     OUTPUT: MH,HM,MCH, SA = SIN(ALPHA), CA= COS(ALPHA), TANBA
C     WHERE MHP AND HPM ARE THE LIGHTEST AND HEAVIEST CP-EVEN
C     HIGGS MASSES, MHCH IS THE CHARGED HIGGS MASS AND
C     ALPHA IS THE HIGGS MIXING ANGLE.
C     TANBA IS THE ANGLE TANB AT THE CP-ODD HIGGS MASS SCALE.

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Program based on the work by M. Carena, M. Quiros
c       and C.E.M. Wagner, "Effective potential methods and
c       the Higgs mass spectrum in the MSSM", Nucl. Phys.
c       B461 (1996) 407. 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE SUBH(MA,TANB,MQ,MUR,MD,MTOP,AU,AD,MU,MCHI0,
     *                MHP,HMP,MCH,SA,CA,TANBA)

      IMPLICIT REAL*8(A-H,L,M,O-Z)
      DIMENSION VH(2,2),M2(2,2),M2P(2,2)
      COMMON/PARAM/GF,ALPH,AMTAU,AMMUON,AMZ,AMW
      COMMON/HSELF/LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,LAMBDA5,
     .             LAMBDA6,LAMBDA7

      MCHI = MCHI0
      TANBA = TANB
      TANBT = TANB
      
      PI = 4*DATAN(1D0)
      MZ = AMZ
      MW = AMW
      V  = 1/DSQRT(2*DSQRT(2D0)*GF)
      CW = AMW**2/AMZ**2
      SW = 1-CW
      ALPHA2  = (2*AMW/V/DSQRT(2D0))**2/4/PI
      ALPHA1  = ALPHA2*SW/CW
      ALPHA3Z = ALPHAS(AMZ,2)
      ALPHA3  = ALPHAS(MTOP,2)
      MB      = RUNM(MTOP,5)
      RMTOP   = RUNM(MTOP,6)

      TQ = LOG((MQ**2+MTOP**2)/MTOP**2)
      TU = LOG((MUR**2 + MTOP**2)/MTOP**2)
      TD = LOG((MD**2 + MTOP**2)/MTOP**2)
      SINB = TANB/DSQRT(1.D0 + TANB**2)
      COSB = SINB/TANB

      IF(MA.GT.MTOP)
     *       TANBA = TANB*(1.D0-3.D0/32.D0/PI**2*
     *       (RMTOP**2/V**2/SINB**2-MB**2/V**2/COSB**2)*
     *       DLOG(MA**2/MTOP**2))
      IF(MA.LT.MTOP.OR.MA.EQ.MTOP) TANBT = TANBA

      SINB = TANBT/DSQRT(1.D0 + TANBT**2)
      COSB = 1.D0/DSQRT(1.D0 + TANBT**2)
      COS2B = (TANBT**2 - 1.D0)/(TANBT**2 + 1.D0)
      G1 = DSQRT(ALPHA1*4.D0*PI)
      G2 = DSQRT(ALPHA2*4.D0*PI)
      G3 = DSQRT(ALPHA3*4.D0*PI)
      HU = RMTOP/V/SINB
      HD =  MB/V/COSB
C

      IF(MQ.GT.MUR) TP = TQ - TU
      IF(MQ.LT.MUR.OR.MQ.EQ.MUR) TP = TU - TQ
      IF(MQ.GT.MUR) TDP = TU
      IF(MQ.LT.MUR.OR.MQ.EQ.MUR) TDP = TQ
      IF(MQ.GT.MD) TPD = TQ - TD
      IF(MQ.LT.MD.OR.MQ.EQ.MD) TPD = TD - TQ
      IF(MQ.GT.MD) TDPD = TD
      IF(MQ.LT.MD.OR.MQ.EQ.MD) TDPD = TQ

      IF(MQ.GT.MD) DLAMBDA1 = 6./96./PI**2*G1**2*HD**2*TPD
      IF(MQ.LT.MD.OR.MQ.EQ.MD) DLAMBDA1 = 3./32./PI**2*
     * HD**2*(G1**2/3.+G2**2)*TPD

      IF(MQ.GT.MUR) DLAMBDA2 =12./96./PI**2*G1**2*HU**2*TP
      IF(MQ.LT.MUR.OR.MQ.EQ.MUR) DLAMBDA2 = 3./32./PI**2*
     * HU**2*(-G1**2/3.+G2**2)*TP

      DLAMBDA3 = 0.
      DLAMBDA4 = 0.

      IF(MQ.GT.MD) DLAMBDA3 = -1./32./PI**2*G1**2*HD**2*TPD
      IF(MQ.LT.MD.OR.MQ.EQ.MD) DLAMBDA3 = 3./64./PI**2*HD**2*
     *(G2**2-G1**2/3.)*TPD
      
      IF(MQ.GT.MUR) DLAMBDA3 = DLAMBDA3 - 
     *1./16./PI**2*G1**2*HU**2*TP
      IF(MQ.LT.MUR.OR.MQ.EQ.MUR) DLAMBDA3 = DLAMBDA3 + 
     * 3./64./PI**2*HU**2*(G2**2+G1**2/3.)*TP

      IF(MQ.LT.MUR) DLAMBDA4 = -3./32./PI**2*G2**2*HU**2*TP
      IF(MQ.LT.MD) DLAMBDA4 = DLAMBDA4 - 3./32./PI**2*G2**2*
     *                        HD**2*TPD
C
      LAMBDA1 = ((G1**2 + G2**2)/4.)*
     *(1.-3.*HD**2*(TPD + TDPD)/8./PI**2)
     *+(3.*HD**4./16./PI**2) *TPD*(1.   
     *+ (3.*HD**2/2. + HU**2/2.       
     *- 8.*G3**2) * (TPD + 2.*TDPD)/16./PI**2) 
     *+(3.*HD**4./8./PI**2) *TDPD*(1.  + (3.*HD**2/2. + HU**2/2.       
     *- 8.*G3**2) * TDPD/16./PI**2) + DLAMBDA1 
C
      LAMBDA2 = ((G1**2 + G2**2)/4.)*(1.-3.*HU**2*
     *(TP + TDP)/8./PI**2)
     *+(3.*HU**4./16./PI**2) *TP*(1.   
     *+ (3.*HU**2/2. + HD**2/2.       
     *- 8.*G3**2) * (TP + 2.*TDP)/16./PI**2) 
     *+(3.*HU**4./8./PI**2) *TDP*(1. + (3.*HU**2/2. + HD**2/2.       
     *- 8.*G3**2) * TDP/16./PI**2) + DLAMBDA2 
C
      LAMBDA3 = ((G2**2 - G1**2)/4.)*(1.-3.*
     *(HU**2)*(TP + TDP)/16./PI**2 -3.*
     *(HD**2)*(TPD + TDPD)/16./PI**2) +DLAMBDA3 
C
      LAMBDA4 = (- G2**2/2.)*(1.
     *-3.*(HU**2)*(TP + TDP)/16./PI**2
     *-3.*(HD**2)*(TPD + TDPD)/16./PI**2) +DLAMBDA4
C     
	LAMBDA5 = 0.
	LAMBDA6 = 0.
	LAMBDA7 = 0.

      M2(1,1) = 2.*V**2*(LAMBDA1*COSB**2+2.*LAMBDA6*
     *COSB*SINB + LAMBDA5*SINB**2) + MA**2*SINB**2
      M2(2,2) = 2.*V**2*(LAMBDA5*COSB**2+2.*LAMBDA7*
     *COSB*SINB + LAMBDA2*SINB**2) + MA**2*COSB**2
      M2(1,2) = 2.*V**2*(LAMBDA6*COSB**2+(LAMBDA3+LAMBDA4)*
     *COSB*SINB + LAMBDA7*SINB**2) - MA**2*SINB*COSB
      M2(2,1) = M2(1,2)

C
C     THIS IS THE CONTRIBUTION FROM LIGHT CHARGINOS/NEUTRALINOS
C     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  	 MSSUSY=DSQRT(0.5D0*(MQ**2+MUR**2)+MTOP**2)
	IF(MCHI.GT.MSSUSY)GOTO 3790
	IF(MCHI.LT.MTOP) MCHI=MTOP
	TCHAR=LOG(MSSUSY**2/MCHI**2)
	DELTAL12=(9./64./PI**2*G2**4+5./192./PI**2*G1**4)*TCHAR
	DELTAL3P4=(3./64./PI**2*G2**4+7./192./PI**2*G1**4
     *       +4./32/PI**2*G1**2*G2**2)*TCHAR
	DELTAM112=2.*DELTAL12*V**2*COSB**2
	DELTAM222=2.*DELTAL12*V**2*SINB**2
	DELTAM122=2.*DELTAL3P4*V**2*SINB*COSB
	M2(1,1)=M2(1,1)+DELTAM112
	M2(2,2)=M2(2,2)+DELTAM222
	M2(1,2)=M2(1,2)+DELTAM122
	M2(2,1)=M2(2,1)+DELTAM122
 3790	CONTINUE
CCCCCCCCCCCCCCC    END OF CHARGINOS AND NEUTRALINOS  CCCCCCCCCCCC 


      CALL GFUN(MA,TANBA,MQ,MUR,MD,MTOP,AU,AD,MU,VH)
      
      DO 9800 I = 1,2
      DO 9801 J = 1,2
      M2P(I,J) = M2(I,J) + VH(I,J)
 9801 CONTINUE
 9800 CONTINUE

      TRM2P  = M2P(1,1) + M2P(2,2)
      DETM2P = M2P(1,1)*M2P(2,2) - M2P(1,2)*M2P(2,1)

      MH2P = (TRM2P - DSQRT(TRM2P**2 - 4.D0* DETM2P))/2.D0
      HM2P = (TRM2P + DSQRT(TRM2P**2 - 4.D0* DETM2P))/2.D0
C !!!!!!!!!!!!!!!!!!!
      MCH2=MA**2+(LAMBDA5-LAMBDA4)*V**2
C !!!!!!!!!!!!!!!!!!!
      MCH=DSQRT(MCH2)
      HMP = DSQRT(HM2P) 
      IF(MH2P.LT.0.)GOTO 5555
      MHP = DSQRT(MH2P) 
C
      SIN2ALPHA = 2.*M2P(1,2)/DSQRT(TRM2P**2-4.D0*DETM2P)
      COS2ALPHA = (M2P(1,1)-M2P(2,2))/DSQRT(TRM2P**2-4.D0*DETM2P)
      IF(COS2ALPHA.GT.0.) ALPHA = DASIN(SIN2ALPHA)/2.D0
      IF(COS2ALPHA.LT.0.) ALPHA = -PI/2.D0-DASIN(SIN2ALPHA)/2.D0
      SA = DSIN(ALPHA)
      CA = DCOS(ALPHA)  
      SQBMA = (SINB*CA - COSB*SA)**2

5555  RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE AMHAMA(ICASE,MH,TANB)
C--CALCULATION OF PSEUDOSCALAR HIGGS MASS FROM HIGGS MASS MH
C--ICASE=0: MH=PSEUDOSCALAR MASS
C--ICASE=1: MH=LIGHT SCALAR MASS
C--ICASE=2: MH=HEAVY SCALAR MASS
C--ICASE=3: MH=CHARGED HIGGS MASS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT REAL*8(A-H,L,M,O-Z)
      DIMENSION VH(2,2),M2(2,2),M2P(2,2)
      COMMON/HMASS/AMSM,AMA,AML,AMH,AMCH,AMAR
      IF(ICASE.EQ.0)THEN
       MA = MH
      ELSE
       DEL0 = 1.D-4
       MA0 = 1.D0
       MA1 = 1.D4
1      MA = (MA0+MA1)/2
C      CALL SUBH(MA,TANB,MQ,MUR,MD,MTOP,AU,AD,MU,MCHI0,
C    *                 MHP,HMP,MCH,SA,CA,TANBA)
       AMA = MA
       CALL SUSYCP(TANB)
       IF(ICASE.EQ.1)THEN
        MX = AML
       ELSEIF(ICASE.EQ.2)THEN
        MX = AMH
       ELSEIF(ICASE.EQ.3)THEN
        MX = AMCH
       ENDIF
       DEL = DABS(MA1 - MA0)/MA
       IF(DEL.GT.DEL0) THEN
        IF(MX.GT.MH) MA1 = MA
        IF(MX.LT.MH) MA0 = MA
        GOTO 1
       ENDIF
       FAC = 1
       MAX = DINT(FAC*MA+0.5D0)/FAC
C      CALL SUBH(MAX,TANB,MQ,MUR,MD,MTOP,AU,AD,MU,MCHI0,
C    *                 MHP,HMP,MCH,SA,CA,TANBA)
       AMA = MAX
       CALL SUSYCP(TANB)
       IF(ICASE.EQ.1)THEN
        MX = AML
       ELSEIF(ICASE.EQ.2)THEN
        MX = AMH
       ELSEIF(ICASE.EQ.3)THEN
        MX = AMCH
       ENDIF
       IF(MX.EQ.MH)THEN
        MA = MAX
       ELSE
        DEL0 = 1.D-8
2       MA = (MA0+MA1)/2
C       CALL SUBH(MA,TANB,MQ,MUR,MD,MTOP,AU,AD,MU,MCHI0,
C    *                  MHP,HMP,MCH,SA,CA,TANBA)
        AMA = MA
        CALL SUSYCP(TANB)
        IF(ICASE.EQ.1)THEN
         MX = AML
        ELSEIF(ICASE.EQ.2)THEN
         MX = AMH
        ELSEIF(ICASE.EQ.3)THEN
         MX = AMCH
        ENDIF
        DEL = DABS(MA1 - MA0)/MA
        IF(DEL.GT.DEL0) THEN
         IF(MX.GT.MH) MA1 = MA
         IF(MX.LT.MH) MA0 = MA
         GOTO 2
        ENDIF
       ENDIF
      ENDIF
      AMA = MA
      CALL SUSYCP(TANB)
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCC NON DEGENERATE STOP/SBOTTOM EFFECTS CCCCCCCCC
C
        SUBROUTINE GFUN(MA,TANB,MQ,MUR,MD,MTOP,AT,AB,MU,VH)
        IMPLICIT REAL*8 (A-H,L,M,O-Z)
        DIMENSION VH(2,2),VH1(2,2),VH2(2,2),
     *            VH3T(2,2),VH3B(2,2),AL(2,2)
        COMMON/PARAM/GF,ALPH,AMTAU,AMMUON,AMZ,AMW
        G(X,Y) = 2.D0 - (X+Y)/(X-Y)*DLOG(X/Y)

        IF(DABS(MU).LT.0.000001) MU = 0.000001
        MQ2   = MQ**2
        MUR2  = MUR**2
        MD2   = MD**2
        TANBA = TANB
        SINBA = TANBA/DSQRT(TANBA**2+1.D0)
        COSBA = SINBA/TANBA        
        SINB = TANB/DSQRT(TANB**2+1.D0)
        COSB = SINB/TANB

      MB = RUNM(MTOP,5)
      PI = 4*DATAN(1D0)
      MZ = AMZ
      MW = AMW
      V  = 1/DSQRT(2*DSQRT(2D0)*GF)
      CW = AMW**2/AMZ**2
      SW = 1-CW
      ALPHA2  = (2*AMW/V/DSQRT(2D0))**2/4/PI
      ALPHA1  = ALPHA2*SW/CW
      ALPHA3Z = ALPHAS(AMZ,2)
      ALPHA3  = ALPHAS(MTOP,2)

      G1 = DSQRT(ALPHA1*4.*PI)
      G2 = DSQRT(ALPHA2*4.*PI)
      G3 = DSQRT(ALPHA3*4.*PI)
      
        IF(MQ.GT.MUR) MST = MQ
        IF(MUR.GT.MQ.OR.MUR.EQ.MQ) MST = MUR
        MSUSYT = DSQRT(MST**2  + MTOP**2)

	IF(MQ.GT.MD) MSB = MQ
	IF(MD.GT.MQ.OR.MD.EQ.MQ) MSB = MD
	MSUSYB = DSQRT(MSB**2 + MB**2)

	TT = LOG(MSUSYT**2/MTOP**2)
	TB = LOG(MSUSYB**2/MTOP**2)

        RMTOP   = RUNM(MTOP,6)

        HT = RMTOP/V/SINB
        HTST = RMTOP/V
        HB =  MB/V/COSB
        G32 = ALPHA3*4.*PI

        BT2 = -(8.*G32 - 9.*HT**2/2. - HB**2/2.)/(4.*PI)**2
	BB2 = -(8.*G32 - 9.*HB**2/2. - HT**2/2.)/(4.*PI)**2
        AL2 = 3./8./PI**2*HT**2
        BT2ST = -(8.*G32 - 9.*HTST**2/2.)/(4.*PI)**2
        ALST = 3./8./PI**2*HTST**2
        AL1 = 3./8./PI**2*HB**2

        AL(1,1) = AL1
        AL(1,2) = (AL2+AL1)/2.
        AL(2,1) = (AL2+AL1)/2.
        AL(2,2) = AL2

        MTOP4 = RMTOP**4.*(1.+2.*BT2*TT- AL2*TT)
        MTOP2 = DSQRT(MTOP4)
	MBOT4 = MB**4.*(1.+2.*BB2*TB - AL1*TB)
	MBOT2 = DSQRT(MBOT4)

	IF(MA.GT.MTOP) THEN
        VI = V*(1. + 3./32./PI**2*HTST**2*LOG(MTOP**2/MA**2))
        H1I = VI*COSBA
        H2I = VI*SINBA
        H1T = H1I*(1.+3./8./PI**2*HB**2*LOG(MA**2/MSUSYT**2))**.25
        H2T = H2I*(1.+3./8./PI**2*HT**2*LOG(MA**2/MSUSYT**2))**.25
        H1B = H1I*(1.+3./8./PI**2*HB**2*LOG(MA**2/MSUSYB**2))**.25
        H2B = H2I*(1.+3./8./PI**2*HT**2*LOG(MA**2/MSUSYB**2))**.25
	ELSE
	VI =  V
	H1I = VI*COSB
	H2I = VI*SINB
        H1T = H1I*(1.+3./8./PI**2*HB**2*LOG(MTOP**2/MSUSYT**2))**.25
        H2T = H2I*(1.+3./8./PI**2*HT**2*LOG(MTOP**2/MSUSYT**2))**.25
        H1B = H1I*(1.+3./8./PI**2*HB**2*LOG(MTOP**2/MSUSYB**2))**.25
        H2B = H2I*(1.+3./8./PI**2*HT**2*LOG(MTOP**2/MSUSYB**2))**.25
	END IF

        TANBST = H2T/H1T
        SINBT = TANBST/(1.+TANBST**2)**.5
        COSBT = SINBT/TANBST

        TANBSB = H2B/H1B
        SINBB = TANBSB/(1.+TANBSB**2)**.5
        COSBB = SINBB/TANBSB

        STOP12 = (MQ2 + MUR2)*.5 + MTOP2 
     *   +1./8.*(G2**2+G1**2)*(H1T**2-H2T**2)
     *   +(((G2**2-5.*G1**2/3.)/4.*(H1T**2-H2T**2) +
     *   MQ2 - MUR2)**2*0.25 + MTOP2*(AT-MU/TANBST)**2)**.5

        STOP22 = (MQ2 + MUR2)*.5 + MTOP2 
     *  +1./8.*(G2**2+G1**2)*(H1T**2-H2T**2) 
     *   - (((G2**2-5.*G1**2/3.)/4.*(H1T**2-H2T**2) +
     *  MQ2 - MUR2)**2*0.25 
     *  + MTOP2*(AT-MU/TANBST)**2)**.5

        IF(STOP22.LT.0.) GOTO 4237

        SBOT12 = (MQ2 + MD2)*.5  
     *   - 1./8.*(G2**2+G1**2)*(H1B**2-H2B**2)
     *  + (((G1**2/3.-G2**2)/4.*(H1B**2-H2B**2) +
     *  MQ2 - MD2)**2*0.25 + MBOT2*(AB-MU*TANBSB)**2)**.5

        SBOT22 = (MQ2 + MD2)*.5  
     *   - 1./8.*(G2**2+G1**2)*(H1B**2-H2B**2)
     *   - (((G1**2/3.-G2**2)/4.*(H1B**2-H2B**2) +
     *   MQ2 - MD2)**2*0.25 + MBOT2*(AB-MU*TANBSB)**2)**.5

        IF(SBOT22.LT.0.) GOTO 4237

        STOP1 = STOP12**.5
        STOP2 = STOP22**.5
        SBOT1 = SBOT12**.5
        SBOT2 = SBOT22**.5

        VH1(1,1) = 1./TANBST
        VH1(2,1) = -1.
        VH1(1,2) = -1.
        VH1(2,2) = TANBST
        VH2(1,1) = TANBST
        VH2(1,2) = -1.
        VH2(2,1) = -1.
        VH2(2,2) = 1./TANBST

C CCCCCCCCCCCCCCCCCCCCCCCCCCC  D-terms CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	STW=SW

	F1T=(MQ2-MUR2)/(STOP12-STOP22)*(.5-4./3.*STW)*
     *         LOG(STOP1/STOP2)
     *        +(.5-2./3.*STW)*LOG(STOP1*STOP2/(MQ2+MTOP2))
     *        + 2./3.*STW*LOG(STOP1*STOP2/(MUR2+MTOP2))

	F1B=(MQ2-MD2)/(SBOT12-SBOT22)*(-.5+2./3.*STW)*
     *        LOG(SBOT1/SBOT2)
     *        +(-.5+1./3.*STW)*LOG(SBOT1*SBOT2/(MQ2+MBOT2))
     *        - 1./3.*STW*LOG(SBOT1*SBOT2/(MD2+MBOT2))

	F2T=MTOP2**.5*(AT-MU/TANBST)/(STOP12-STOP22)*
     *         (-.5*LOG(STOP12/STOP22)
     *        +(4./3.*STW-.5)*(MQ2-MUR2)/(STOP12-STOP22)*
     *         G(STOP12,STOP22))

	F2B=MBOT2**.5*(AB-MU*TANBSB)/(SBOT12-SBOT22)*
     *         (.5*LOG(SBOT12/SBOT22)
     *        +(-2./3.*STW+.5)*(MQ2-MD2)/(SBOT12-SBOT22)*
     *        G(SBOT12,SBOT22))

        VH3B(1,1) = MBOT4/(COSBB**2)*(LOG(SBOT1**2*SBOT2**2/
     *  (MQ2+MBOT2)/(MD2+MBOT2)) 
     *  + 2.*(aB*(aB-MU*TANBSB)/(SBOT1**2-SBOT2**2))*
     *  LOG(SBOT1**2/SBOT2**2)) + 
     *  mBOT4/(COSBB**2)*(aB*(aB-MU*TANBSB)/
     *  (SBOT1**2-SBOT2**2))**2*G(SBOT12,SBOT22) 

	VH3T(1,1) =
     *  MTOP4/(SINBT**2)*(MU*(-aT+MU/TANBST)/(STOP1**2
     * -STOP2**2))**2*G(STOP12,STOP22)

	VH3B(1,1)=VH3B(1,1)+
     *    MZ**2*(2*MBOT2*F1B-MBOT2**.5*AB*F2B)

	VH3T(1,1) = VH3T(1,1) + 
     *  MZ**2*(MTOP2**.5*MU/TANBST*F2T)  

        VH3T(2,2) = MTOP4/(SINBT**2)*(LOG(STOP1**2*STOP2**2/
     *  (MQ2+MTOP2)/(MUR2+MTOP2)) 
     *  + 2.*(aT*(aT-MU/TANBST)/(STOP1**2-STOP2**2))*
     *  LOG(STOP1**2/STOP2**2)) + 
     *  MTOP4/(SINBT**2)*(aT*(aT-MU/TANBST)/
     *  (STOP1**2-STOP2**2))**2*G(STOP12,STOP22) 

	VH3B(2,2) =
     *  mBOT4/(COSBB**2)*(MU*(-aB+MU*TANBSB)/(SBOT1**2
     * -SBOT2**2))**2*G(SBOT12,SBOT22)

	VH3T(2,2)=VH3T(2,2)+
     *    MZ**2*(-2*MTOP2*F1T+MTOP2**.5*AT*F2T)

	VH3B(2,2) = VH3B(2,2) -MZ**2*MBOT2**.5*MU*TANBSB*F2B  
 
        VH3T(1,2) = -
     *   MTOP4/(SINBT**2)*MU*(aT-MU/TANBST)/
     * (STOP1**2-STOP2**2)*(LOG(STOP1**2/STOP2**2) + aT*
     * (AT - MU/TANBST)/(STOP1**2-STOP2**2)*G(STOP12,STOP22))

	VH3B(1,2) =
     * - MBOT4/(COSBB**2)*MU*(aT-MU*TANBSB)/
     * (SBOT1**2-SBOT2**2)*(LOG(SBOT1**2/SBOT2**2) + AB*
     * (AB - MU*TANBSB)/(SBOT1**2-SBOT2**2)*G(SBOT12,SBOT22))

	VH3T(1,2)=VH3T(1,2) +
     *      MZ**2*(MTOP2/TANBST*F1T-MTOP2**.5*(AT/TANBST+MU)/2.*F2T)

	VH3B(1,2)=VH3B(1,2)
     *  +MZ**2*(-MBOT2*TANBSB*F1B+MBOT2**.5*(AB*TANBSB+MU)/2.*F2B)

        VH3T(2,1) = VH3T(1,2)
        VH3B(2,1) = VH3B(1,2)

       TQ = LOG((MQ2 + MTOP2)/MTOP2)
       TU = LOG((MUR2+MTOP2)/MTOP2)
       TQD = LOG((MQ2 + MB**2)/MB**2)
       TD = LOG((MD2+MB**2)/MB**2)

        DO 8910 I = 1,2
        DO 8911 J = 1,2
        VH(I,J) = 
     *  6./(8.*PI**2*(H1T**2+H2T**2))
     *  *VH3T(I,J)*0.5*(1.-AL(I,J)*TT/2.) +
     *  6./(8.*PI**2*(H1B**2+H2B**2))
     *  *VH3B(I,J)*0.5*(1.-AL(I,J)*TB/2.)
 8911   CONTINUE
 8910   CONTINUE

        GOTO 4236
 4237   DO 6868 I =1,2
        DO 6867 J = 1,2
        VH(I,J) = -1.D+15
 6867   CONTINUE
 6868   CONTINUE

4236    RETURN
        END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       End of program from M. Carena, M. Quiros and C.E.M. Wagner.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     ****************************************************************
C	  CHARGINO AND NEUTRALINO MASS MATRICES AND COUPLINGS
C     ****************************************************************
      SUBROUTINE GAUGINO(MU,M2,B,A,MC,MN,XMN,AC1,AC2,AC3,AN1,AN2,AN3
     .                 ,ACNL,ACNR)            
      IMPLICIT REAL*8(A-H,K-Z)
      COMPLEX*16 CXA,CXB,CXC,CXD,CX1,CX2,CX3
      DIMENSION MC(2),MN(4),XMN(4),Z(4,4),ZX(4,4),U(2,2),V(2,2),
     .          QQ(4,4),SS(4,4),S(2,2),Q(2,2),AC1(2,2),AC2(2,2),
     .          AC3(2,2),AN1(4,4),AN2(4,4),AN3(4,4),ACNL(2,4),
     .          ACNR(2,4),IORD(4),IREM(2)
      DIMENSION X(2,2)
      DIMENSION YMN(4),YZ(4,4),XMC(2),BU(2),BV(2)
      COMMON/PARAM/GF,ALPH,AMTAU,AMMUON,MZ,MW
      CW=MW/MZ
      SW=DSQRT(1-CW**2)
      PI=4.D0*DATAN(1.D0)
      SB=DSIN(B)
      CB=DCOS(B)
      TW=SW/CW
      M1=5.D0/3.D0*TW**2*M2
C     ************  NEUTRALINO MASSES AND MATRIX ELEMENTS ***********
      EPS=-1.D-10
      XC2=(M1*M2-MZ**2-MU**2)-3.D0/8.D0*(M1+M2)**2
      XC3=-1.D0/8.D0*(M1+M2)**3+1.D0/2.D0*(M1+M2)*(M1*M2-MZ**2
     .    -MU**2)+(M1+M2)*MU**2+(M1*CW**2+M2*SW**2)*MZ**2
     .    -MU*MZ**2*DSIN(2.D0*B)
      XC4=+(M1*CW**2+M2*SW**2)*MU*MZ**2*DSIN(2.D0*B)-M1*M2*MU**2
     .    +1.D0/4.D0*(M1+M2)*( (M1+M2)*MU**2+(M1*CW**2+M2*SW**2)
     .    *MZ**2-MU*MZ**2*DSIN(2.D0*B) )+1.D0/16.D0*(M1+M2)**2*
     .    (M1*M2-MZ**2-MU**2)-3.D0/256.D0*(M1+M2)**4
      XS=-XC3**2-2.D0/27.D0*XC2**3+8.D0/3.D0*XC2*XC4
      XU=-1.D0/3.D0*XC2**2-4.D0*XC4
      CXD=(-4*XU**3-27*XS**2)*DCMPLX(1.D0,EPS)
      CXC=1.D0/2.D0*(-XS+DCMPLX(0.D0,1.D0)*CDSQRT(CXD/27.D0))
      CXA=DREAL(CXC**(1.D0/3.D0))*DCMPLX(1.D0,-EPS)
      CXB=8.D0*CXA-8.D0/3.D0*XC2*DCMPLX(1.D0,-EPS)
C     *********** MASSES AND COUPLINGS:
      X0=(M1+M2)/4.D0
      CX1= CXA/2.D0-XC2/6.D0*DCMPLX(1.D0,-EPS)
      CX2=-CXA/2.D0-XC2/3.D0*DCMPLX(1.D0,-EPS)
      CX3=XC3*DCMPLX(1.D0,-EPS)/CDSQRT(CXB)
      XMN(1)=X0-CDABS(CDSQRT(CX1))+CDABS(CDSQRT(CX2+CX3))
      XMN(2)=X0+CDABS(CDSQRT(CX1))-CDABS(CDSQRT(CX2-CX3))
      XMN(3)=X0-CDABS(CDSQRT(CX1))-CDABS(CDSQRT(CX2+CX3))
      XMN(4)=X0+CDABS(CDSQRT(CX1))+CDABS(CDSQRT(CX2-CX3))
      DO 10 I=1,4
      MN(I)=DABS(XMN(I))
      YMN(I)=XMN(I)
      ZX(I,2)=-CW/SW*(M1-XMN(I))/(M2-XMN(I))
      ZX(I,3)=(MU*(M2-XMN(I))*(M1-XMN(I))-MZ**2*SB*CB*((M1-M2)*CW**2
     .       +M2-XMN(I)))/MZ/(M2-XMN(I))/SW/(MU*CB+XMN(I)*SB)
      ZX(I,4)=(-XMN(I)*(M2-XMN(I))*(M1-XMN(I))-MZ**2*CB*CB*((M1-M2)
     .       *CW**2+M2-XMN(I)))/MZ/(M2-XMN(I))/SW/(MU*CB+XMN(I)*SB)
      ZX(I,1)=1.D0/DSQRT(1.D0+ZX(I,2)**2+ZX(I,3)**2+ZX(I,4)**2) 
      YZ(I,1)=ZX(I,1)
      YZ(I,2)=ZX(I,2)*ZX(I,1)
      YZ(I,3)=ZX(I,3)*ZX(I,1)
      YZ(I,4)=ZX(I,4)*ZX(I,1)
 10   CONTINUE
C     *************  ORDERING THE DISORDER ******************
      XX0 = DMIN1(MN(1),MN(2),MN(3),MN(4))
      XX1 = DMAX1(MN(1),MN(2),MN(3),MN(4))
      IDUMMY = 1
      DO I = 1,4
       IF(MN(I).EQ.XX0)THEN
        IORD(1) = I
       ELSEIF(MN(I).EQ.XX1)THEN
        IORD(4) = I
       ELSE
        IREM(IDUMMY) = I
        IDUMMY = IDUMMY+1
       ENDIF
      ENDDO
      IF(MN(IREM(1)).LE.MN(IREM(2)))THEN
       IORD(2) = IREM(1)
       IORD(3) = IREM(2)
      ELSE
       IORD(2) = IREM(2)
       IORD(3) = IREM(1)
      ENDIF
C 
      DO 98 J=1,4
      I=IORD(J)
      XMN(J)=YMN(I)
      MN(J) =DABS(YMN(I))
        DO I1=1,4
        Z(J,I1)=YZ(I,I1)
        ENDDO
 98   CONTINUE
C     ************  NEUTRALINO COUPLINGS TO HIGGS BOSONS ***********
	DO 11 I=1,4
	DO 11 J=1,4
	QQ(I,J)=1.D0/2.D0*(Z(I,3)*(Z(J,2)-TW*Z(J,1))+Z(J,3)*
     .		(Z(I,2)-TW*Z(I,1)))
	SS(I,J)=1.D0/2.D0*(Z(I,4)*(Z(J,2)-TW*Z(J,1))+Z(J,4)*
     .		(Z(I,2)-TW*Z(I,1)))
 11	CONTINUE
	DO 21 I=1,4
	DO 21 J=1,4
	AN1(I,J)= QQ(I,J)*DCOS(A)-SS(I,J)*DSIN(A)
	AN2(I,J)=-QQ(I,J)*DSIN(A)-SS(I,J)*DCOS(A)
	AN3(I,J)= QQ(I,J)*DSIN(B)-SS(I,J)*DCOS(B)
 21	CONTINUE

C       ************* CHARGINO MASSES AND MATRIX ELEMENTS ***********
	DELTA=DABS(B-.25*PI)
	DDD=MU*DCOS(B)+M2*DSIN(B)
	CCC=MU*DSIN(B)+M2*DCOS(B)
	IF(DELTA.LT.0.01D0) THEN
	PHIM=PI/4.D0-.5D0*DATAN((M2-MU)/(2.D0*MW))
	PHIP=PHIM
	ELSE IF	(DABS(CCC).LT.1.D-5) THEN
	PHIM=0.D0
	PHIP=DATAN(DSQRT(2.D0)*MW*DSIN(B)/(M2+1.D-5))
	ELSE IF	(DABS(DDD).LT.1.D-5) THEN
	PHIP=0.D0
	PHIM=DATAN(DSQRT(2.D0)*MW*DCOS(B)/(M2+1.D-5))
	ELSE
	RAD=DSQRT((M2**2-MU**2)**2+4.D0*MW**4*DCOS(2.D0*B)**2
     +	+4.D0*MW**2*(M2**2+MU**2+2.D0*M2*MU*DSIN(2.D0*B)))
	PHIP=DATAN((RAD-(M2**2-MU**2+2.D0*MW**2*DCOS(2.D0*B)))
     +	/(2.D0*DSQRT(2.D0)*MW*(MU*DCOS(B)+M2*DSIN(B))))
	PHIM=DATAN((RAD-(M2**2-MU**2-2.D0*MW**2*DCOS(2.D0*B)))
     +	/(2.D0*DSQRT(2.D0)*MW*(MU*DSIN(B)+M2*DCOS(B))))
	ENDIF
	CP=DCOS(PHIP)
	SP=DSIN(PHIP)
	CM=DCOS(PHIM)
	SM=DSIN(PHIM)
C  MY CONVENTION
	U(2,2)=CM
	U(2,1)=-SM
	U(1,2)=SM
	U(1,1)=CM
	V(1,1)=CP
	V(1,2)=SP
	V(2,1)=-SP
	V(2,2)=CP
        X(1,1)=M2
        X(1,2)=DSQRT(2.D0)*MW*DSIN(B)
        X(2,1)=DSQRT(2.D0)*MW*DCOS(B)
        X(2,2)=MU
 555    CONTINUE
       XMC(1)=(U(1,1)*X(1,1)+U(1,2)*X(2,1))*V(1,1)
     .       +(U(1,1)*X(1,2)+U(1,2)*X(2,2))*V(1,2)
       XMC(2)=(U(2,1)*X(1,1)+U(2,2)*X(2,1))*V(2,1)
     .       +(U(2,1)*X(1,2)+U(2,2)*X(2,2))*V(2,2)
        IF(XMC(1).LT.0.D0) THEN
	V(1,1)=-CP
	V(1,2)=-SP
	V(2,1)=-SP
	V(2,2)=CP
        GOTO 555
        ENDIF
        IF(XMC(2).LT.0.D0) THEN
	V(1,1)=CP
	V(1,2)=SP
	V(2,1)=SP
	V(2,2)=-CP
        GOTO 555
        ENDIF
        IF(XMC(1).GT.XMC(2)) THEN
        MTEMP=XMC(1)
        XMC(1)=XMC(2)
        XMC(2)=MTEMP
        DO J=1,2
        BU(J)=U(1,J)
        U(1,J)=U(2,J)
        U(2,J)=BU(J)
        BV(J)=V(1,J)
        V(1,J)=V(2,J)
        V(2,J)=BV(J)
        ENDDO
        ENDIF        
        MC(1)=DABS(XMC(1))
        MC(2)=DABS(XMC(2))

C     ************  CHARGINO COUPLINGS TO HIGGS BOSONS ***********
	DO 12 I=1,2
	DO 12 J=1,2
	Q(I,J)=DSQRT(1.D0/2.D0)*U(J,2)*V(I,1)
	S(I,J)=DSQRT(1.D0/2.D0)*U(J,1)*V(I,2)
 12	CONTINUE
	DO 22 I=1,2
	DO 22 J=1,2	
	AC1(I,J)= Q(I,J)*DCOS(A)+S(I,J)*DSIN(A)
	AC2(I,J)=-Q(I,J)*DSIN(A)+S(I,J)*DCOS(A)
	AC3(I,J)= Q(I,J)*DSIN(B)+S(I,J)*DCOS(B)
 22	CONTINUE
C     **** CHARGINO-NEUTRALINO COUPLINGS TO CHARGED HIGGS BOSONS 
	DO 13 I=1,2
	DO 13 J=1,4
        ACNL(I,J)=DCOS(B)*(Z(J,4)*V(I,1)+(Z(J,2)+Z(J,1)*TW)
     .       *V(I,2)/DSQRT(2.D0)) 
        ACNR(I,J)=DSIN(B)*(Z(J,3)*U(I,1)-(Z(J,2)+Z(J,1)*TW)
     .       *U(I,2)/DSQRT(2.D0)) 
 13     CONTINUE

       RETURN
       END

C   ****************************************************************
C     SUBROUTINE FOR SFERMION MASSES, MIXING AND COUPLINGS 
C   ****************************************************************

       SUBROUTINE SFERMION(MQL,MUR,MDR,MEL,MER,AL,AT,AB,MU,
     .                    MST,MSB,MSL,MSU,MSD,MSE,MSN, 
     .                    GLEE,GLTT,GLBB,GHEE,GHTT,GHBB,
     .                    GAEE,GATT,GABB,GCEN,GCTB)

      IMPLICIT REAL*8(A-H,K-Z)
      DIMENSION MST(2),MSB(2),MSL(2),MSU(2),MSD(2),MSE(2),MSN(2),
     .          GCEN(2,2),GCTB(2,2),GLEE(2,2),GLTT(2,2),GLBB(2,2),
     .          GHEE(2,2),GHTT(2,2),GHBB(2,2)
      COMMON/MASSES/AMS,AMC,AMB,AMT
      COMMON/PARAM/GF,ALPH,AMTAU,AMMUON,MZ,MW
      COMMON/COUP/GAT,GAB,GLT,GLB,GHT,GHB,GZAH,GZAL,GHHH,GLLL,GHLL,
     .            GLHH,GHAA,GLAA,GLVV,GHVV,GLPM,GHPM,B,A
      COMMON/SFER1ST/MQL1,MUR1,MDR1,MEL1,MER1
C
      PI = 4*DATAN(1.D0)
      SW2=1.D0-MW**2/MZ**2
      TB=DTAN(B)
      MB = RUNM(AMT,5)
      MT = RUNM(AMT,6)
      ML = AMTAU
C FIRST TWO GENERATIONS:  NO MIXING INCLUDED 
C UP SQUARKS: 
      MSTL2=MQL1**2+(0.5D0-2.D0/3.D0*SW2)*MZ**2*DCOS(2.D0*B)
      MSTR2=MUR1**2+2.D0/3.D0*SW2*MZ**2*DCOS(2.D0*B) 
      MSU(1)=DSQRT(MSTL2)
      MSU(2)=DSQRT(MSTR2)
C DOWN SQUARKS
      MSBL2=MQL1**2+(-0.5D0+1.D0/3.D0*SW2)*MZ**2*DCOS(2.D0*B)
      MSBR2=MDR1**2-1.D0/3.D0*SW2*MZ**2*DCOS(2.D0*B) 
      MSD(1)=DSQRT(MSBL2)
      MSD(2)=DSQRT(MSBR2)
C SLEPTONS
      MSEL2=MEL1**2+(-0.5D0+SW2)*MZ**2*DCOS(2.D0*B)
      MSER2=MER1**2- SW2*MZ**2*DCOS(2.D0*B) 
      MSNL2=MEL1**2+0.5D0*MZ**2*DCOS(2.D0*B)
      MSE(1)=DSQRT(MSEL2)
      MSE(2)=DSQRT(MSER2)
      MSN(1)=DSQRT(MSNL2)
      MSN(2)=1.D+15

C NOW THE THIRD GENERATION
C
C STOP MASSES/MIXING
C
      MSTL2=MQL**2+(0.5D0-2.D0/3.D0*SW2)*MZ**2*DCOS(2.D0*B)
      MSTR2=MUR**2+2.D0/3.D0*SW2*MZ**2*DCOS(2.D0*B) 
      MLRT=AT-MU/TB
      DELT=(MSTL2-MSTR2)**2+4*MT**2*MLRT**2
      MST12=MT**2+0.5D0*(MSTL2+MSTR2-DSQRT(DELT))
      MST22=MT**2+0.5D0*(MSTL2+MSTR2+DSQRT(DELT))
        IF(MST12.LT.0.D0)THEN 
      PRINT *, 'MSTOP**2 is negative!!!!'
      GOTO 111 
      ELSE 
      MST(1)=DSQRT(MST12)
      MST(2)=DSQRT(MST22)
      THET=0.5D0*DATAN(2.D0*MT*MLRT / (MSTL2-MSTR2) )
      IF(MSTL2.GT.MSTR2) THET = THET + PI/2
        ENDIF 
      CT= DCOS(THET)
      ST= DSIN(THET) 
C
C SBOTTOM MASSES/MIXING
C
      MSBL2=MQL**2+(-0.5D0+1.D0/3.D0*SW2)*MZ**2*DCOS(2.D0*B)
      MSBR2=MDR**2-1.D0/3.D0*SW2*MZ**2*DCOS(2.D0*B) 
      MLRB=AB-MU*TB
      DELB=(MSBL2-MSBR2)**2+4*MB**2*MLRB**2
      MSB12=MB**2+0.5D0*(MSBL2+MSBR2-DSQRT(DELB))
      MSB22=MB**2+0.5D0*(MSBL2+MSBR2+DSQRT(DELB))
        IF(MSB12.LT.0.D0)THEN
      PRINT *, 'MSBOT**2 is negative!!!!'
      GOTO 111
        ELSE
      MSB(1)=DSQRT(MSB12)
      MSB(2)=DSQRT(MSB22)
      THEB=0.5D0*DATAN(2.D0*MB*MLRB / (MSBL2-MSBR2) )
      IF(MSBL2.GT.MSBR2) THEB = THEB + PI/2
        ENDIF  
      CB= DCOS(THEB)
      SB= DSIN(THEB) 
C
C  STAU MASSES/MIXING
C
      MSEL2=MEL**2+(-0.5D0+SW2)*MZ**2*DCOS(2.D0*B)
      MSER2=MER**2- SW2*MZ**2*DCOS(2.D0*B) 
      MSNL2=MEL**2+0.5D0*MZ**2*DCOS(2.D0*B)
      MLRE=AL-MU*TB
      DELE=(MSEL2-MSER2)**2+4*ML**2*MLRE**2
      MSE12=ML**2+0.5D0*(MSEL2+MSER2-DSQRT(DELE))
      MSE22=ML**2+0.5D0*(MSEL2+MSER2+DSQRT(DELE))
        IF(MSE12.LT.0.D0)THEN
      PRINT *, 'MSTAU**2 is negative!!!!'
      GOTO 111
        ELSE
      MSL(1)=DSQRT(MSE12)
      MSL(2)=DSQRT(MSE22)
      THEL=0.5D0*DATAN(2.D0*ML*MLRE / (MSEL2-MSER2) )
      IF(MSEL2.GT.MSER2) THEL = THEL + PI/2
        ENDIF  
      CL= DCOS(THEL)
      SL= DSIN(THEL) 
C
C LIGHT CP--EVEN HIGGS COUPLINGS TO STOPS
C 
      GLTT(1,1)=-DSIN(B+A)*(0.5D0*CT**2-2.D0/3.D0*SW2*DCOS(2*THET)) 
     .    + MT**2/MZ**2*GLT + MT*ST*CT/MZ**2*(AT*GLT+MU*GHT)
      GLTT(2,2)=-DSIN(B+A)*(0.5D0*ST**2+2.D0/3.D0*SW2*DCOS(2*THET))
     .    + MT**2/MZ**2*GLT - MT*ST*CT/MZ**2*(AT*GLT+MU*GHT)
      GLTT(1,2)=-2*DSIN(B+A)*ST*CT*(2.D0/3.D0*SW2-0.25D0)
     .    + MT*DCOS(2*THET)/2.D0/MZ**2*(AT*GLT+MU*GHT) 
      GLTT(2,1)=-2*DSIN(B+A)*ST*CT*(2.D0/3.D0*SW2-0.25D0)
     .    + MT*DCOS(2*THET)/2.D0/MZ**2*(AT*GLT+MU*GHT) 
C
C LIGHT CP--EVEN HIGGS COUPLINGS TO SBOTTOMS
C
      GLBB(1,1)=-DSIN(B+A)*(-0.5D0*CB**2+1.D0/3.D0*SW2*DCOS(2*THEB)) 
     .    + MB**2/MZ**2*GLB + MB*SB*CB/MZ**2*(AB*GLB-MU*GHB)
      GLBB(2,2)=-DSIN(B+A)*(-0.5D0*SB**2-1.D0/3.D0*SW2*DCOS(2*THEB)) 
     .    + MB**2/MZ**2*GLB - MB*SB*CB/MZ**2*(AB*GLB-MU*GHB)
      GLBB(1,2)=-2*DSIN(B+A)*SB*CB*(-1.D0/3.D0*SW2+0.25D0)
     .    + MB*DCOS(2*THEB)/2.D0/MZ**2*(AB*GLB-MU*GHB) 
      GLBB(2,1)=-2*DSIN(B+A)*SB*CB*(-1.D0/3.D0*SW2+0.25D0)
     .    + MB*DCOS(2*THEB)/2.D0/MZ**2*(AB*GLB-MU*GHB) 

C
C LIGHT CP--EVEN HIGGS COUPLINGS TO STAU'S 
C
      GLEE(1,1)=-DSIN(B+A)*(-0.5D0*CL**2+SW2*DCOS(2*THEL)) 
     .    + ML**2/MZ**2*GLB + ML*SL*CL/MZ**2*(AL*GLB-MU*GHB)
      GLEE(2,2)=-DSIN(B+A)*(-0.5D0*SL**2-SW2*DCOS(2*THEL)) 
     .    + ML**2/MZ**2*GLB - ML*SL*CL/MZ**2*(AL*GLB-MU*GHB)
      GLEE(1,2)=-2*DSIN(B+A)*SL*CL*(-SW2+0.25D0)
     .    + ML*DCOS(2*THEL)/2.D0/MZ**2*(AL*GLB-MU*GHB) 
      GLEE(2,1)=-2*DSIN(B+A)*SL*CL*(-SW2+0.25D0)
     .    + ML*DCOS(2*THEL)/2.D0/MZ**2*(AL*GLB-MU*GHB) 
C
C HEAVY CP--EVEN HIGGS COUPLINGS TO STOPS
C
      GHTT(1,1)=DCOS(B+A)*(0.5D0*CT**2-2.D0/3.D0*SW2*DCOS(2*THET)) 
     .    + MT**2/MZ**2*GHT + MT*ST*CT/MZ**2*(AT*GHT-MU*GLT)
      GHTT(2,2)= DCOS(B+A)*(0.5D0*ST**2+2.D0/3.D0*SW2*DCOS(2*THET))
     .    + MT**2/MZ**2*GHT - MT*ST*CT/MZ**2*(AT*GHT-MU*GLT)
      GHTT(1,2)=2*DCOS(B+A)*ST*CT*(2.D0/3.D0*SW2-0.25D0)
     .    + MT*DCOS(2*THET)/2.D0/MZ**2*(AT*GHT-MU*GLT) 
      GHTT(2,1)=2*DCOS(B+A)*ST*CT*(2.D0/3.D0*SW2-0.25D0)
     .    + MT*DCOS(2*THET)/2.D0/MZ**2*(AT*GHT-MU*GLT) 
C
C HEAVY CP--EVEN HIGGS COUPLINGS TO SBOTTOMS
C
      GHBB(1,1)= DCOS(B+A)*(-0.5D0*CB**2+1.D0/3.D0*SW2*DCOS(2*THEB)) 
     .    + MB**2/MZ**2*GHB + MB*SB*CB/MZ**2*(AB*GHB+MU*GLB)
      GHBB(2,2)= DCOS(B+A)*(-0.5D0*SB**2-1.D0/3.D0*SW2*DCOS(2*THEB)) 
     .    + MB**2/MZ**2*GHB - MB*SB*CB/MZ**2*(AB*GHB+MU*GLB)
      GHBB(1,2)=2*DCOS(B+A)*SB*CB*(-1.D0/3.D0*SW2+0.25D0)
     .    + MB*DCOS(2*THEB)/2.D0/MZ**2*(AB*GHB+MU*GLB) 
      GHBB(2,1)=2*DCOS(B+A)*SB*CB*(-1.D0/3.D0*SW2+0.25D0)
     .    + MB*DCOS(2*THEB)/2.D0/MZ**2*(AB*GHB+MU*GLB) 
C
C HEAVY CP--EVEN HIGGS COUPLINGS TO STAU'S 
C
      GHEE(1,1)= DCOS(B+A)*(-0.5D0*CL**2+SW2*DCOS(2*THEL)) 
     .    + ML**2/MZ**2*GHB + ML*SL*CL/MZ**2*(AL*GHB+MU*GLB)
      GHEE(2,2)= DCOS(B+A)*(-0.5D0*SL**2-SW2*DCOS(2*THEL)) 
     .    + ML**2/MZ**2*GHB - ML*SL*CL/MZ**2*(AL*GHB+MU*GLB)
      GHEE(1,2)=2*DCOS(B+A)*SL*CL*(-SW2+0.25D0)
     .    + ML*DCOS(2*THEL)/2.D0/MZ**2*(AL*GHB+MU*GLB) 
      GHEE(2,1)=2*DCOS(B+A)*SL*CL*(-SW2+0.25D0)
     .    + ML*DCOS(2*THEL)/2.D0/MZ**2*(AL*GHB+MU*GLB) 

C
C PSEUDOSCALAR COUPLINGS 
C
      GATT=-MT/2.D0/MZ**2*(MU+AT*GAT) 
      GABB=-MB/2.D0/MZ**2*(MU+AB*GAB) 
      GAEE=-ML/2.D0/MZ**2*(MU+AL*GAB) 
C
C CHARGED HIGGS COUPLINGS STOPS/SBOTTOMS 
C
      CLL=(MW**2*DSIN(2*B)-MT**2*GAT-MB**2*GAB)/DSQRT(2.D0)/MW**2
      CRR=-MT*MB*(GAT+GAB)/DSQRT(2.D0)/MW**2
      CLR=-MB*(MU+AB*GAB)/DSQRT(2.D0)/MW**2
      CRL=-MT*(MU+AT*GAT)/DSQRT(2.D0)/MW**2
      GCTB(1,1)=+CT*CB*CLL+ST*SB*CRR+CT*SB*CLR+ST*CB*CRL
      GCTB(1,2)=-CT*SB*CLL+ST*CB*CRR+CT*CB*CLR-ST*SB*CRL
      GCTB(2,1)=-ST*CB*CLL+CT*SB*CRR-ST*SB*CLR+CT*CB*CRL
      GCTB(2,2)=+ST*SB*CLL+CT*CB*CRR-ST*CB*CLR-CT*SB*CRL

C
C CHARGED HIGGS COUPLINGS TAU'S AND NEUTRINOS 
C
      CLL=(MW**2*DSIN(2*B)-ML**2*GAB)/DSQRT(2.D0)/MW**2
      CLR=-ML*(MU+AL*GAB)/DSQRT(2.D0)/MW**2
      GCEN(1,1)=CL*CLL+SL*CLR
      GCEN(1,2)=-SL*CLL+CL*CLR
      GCEN(2,1)=0.D0
      GCEN(2,2)=0.D0 

      RETURN 
111   STOP
      END 

C ******************************************************************

C      DOUBLE PRECISION FUNCTION RUNP(Q,NF)
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      COMMON/RUN/XMSB,XMHAT,XKFAC
C      RUNP = RUNM(Q,NF)
C      RUNP = RUNM(Q/2.D0,NF)*XKFAC
C      RETURN
C      END

      DOUBLE PRECISION FUNCTION RUNM(Q,NF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NN=6)
      PARAMETER (ZETA3 = 1.202056903159594D0)
      DIMENSION AM(NN),YMSB(NN)
      COMMON/ALS/XLAMBDA,AMCA,AMBA,AMTA,N0A
      COMMON/MASSES/AMS,AMC,AMB,AMT
      COMMON/STRANGE/AMSB
      COMMON/RUN/XMSB,XMHAT,XKFAC
      COMMON/FLAG/IHIGGS,NNLO,IPOLE
      SAVE ISTRANGE
      B0(NF)=(33.D0-2.D0*NF)/12D0
      B1(NF) = (102D0-38D0/3D0*NF)/16D0
      B2(NF) = (2857D0/2D0-5033D0/18D0*NF+325D0/54D0*NF**2)/64D0
      G0(NF) = 1D0
      G1(NF) = (202D0/3D0-20D0/9D0*NF)/16D0
      G2(NF) = (1249D0-(2216D0/27D0+160D0/3D0*ZETA3)*NF
     .       - 140D0/81D0*NF**2)/64D0
      C1(NF) = G1(NF)/B0(NF) - B1(NF)*G0(NF)/B0(NF)**2
      C2(NF) = ((G1(NF)/B0(NF) - B1(NF)*G0(NF)/B0(NF)**2)**2
     .       + G2(NF)/B0(NF) + B1(NF)**2*G0(NF)/B0(NF)**3
     .       - B1(NF)*G1(NF)/B0(NF)**2 - B2(NF)*G0(NF)/B0(NF)**2)/2D0
      TRAN(X,XK)=1D0+4D0/3D0*ALPHAS(X,2)/PI+XK*(ALPHAS(X,2)/PI)**2
      CQ(X,NF)=(2D0*B0(NF)*X)**(G0(NF)/B0(NF))
     .            *(1D0+C1(NF)*X+C2(NF)*X**2)
      DATA ISTRANGE/0/
      PI=4D0*DATAN(1D0)
      ACC = 1.D-8
      AM(1) = 0
      AM(2) = 0
C--------------------------------------------
      IMSBAR = 0
      IF(IMSBAR.EQ.1)THEN
       IF(ISTRANGE.EQ.0)THEN
C--STRANGE POLE MASS FROM MSBAR-MASS AT 1 GEV
        AMSD = XLAMBDA
        AMSU = 1.D8
123     AMS  = (AMSU+AMSD)/2
        AM(3) = AMS
        XMSB = AMS/CQ(ALPHAS(AMS,2)/PI,3)
     .            *CQ(ALPHAS(1.D0,2)/PI,3)/TRAN(AMS,0D0)
        DD = (XMSB-AMSB)/AMSB
        IF(DABS(DD).GE.ACC)THEN
         IF(DD.LE.0.D0)THEN
          AMSD = AM(3)
         ELSE
          AMSU = AM(3)
         ENDIF
         GOTO 123
        ENDIF
        ISTRANGE=1
       ENDIF
       AM(3) = AMSB
      ELSE
       AMS=AMSB
       AM(3) = AMS
      ENDIF
C--------------------------------------------
      AM(3) = AMSB
      AM(4) = AMC
      AM(5) = AMB
      AM(6) = AMT
      XK = 16.11D0
      DO 1 I=1,NF-1
       XK = XK - 1.04D0*(1.D0-AM(I)/AM(NF))
1     CONTINUE
      IF(NF.GE.4)THEN
       XMSB = AM(NF)/TRAN(AM(NF),0D0)
       XMHAT = XMSB/CQ(ALPHAS(AM(NF),2)/PI,NF)
      ELSE
       XMSB = 0
       XMHAT = 0
      ENDIF
      YMSB(3) = AMSB
      IF(NF.EQ.3)THEN
       YMSB(4) = YMSB(3)*CQ(ALPHAS(AM(4),2)/PI,3)/
     .                   CQ(ALPHAS(1.D0,2)/PI,3)
       YMSB(5) = YMSB(4)*CQ(ALPHAS(AM(5),2)/PI,4)/
     .                   CQ(ALPHAS(AM(4),2)/PI,4)
       YMSB(6) = YMSB(5)*CQ(ALPHAS(AM(6),2)/PI,5)/
     .                   CQ(ALPHAS(AM(5),2)/PI,5)
      ELSEIF(NF.EQ.4)THEN
       YMSB(4) = XMSB
       YMSB(5) = YMSB(4)*CQ(ALPHAS(AM(5),2)/PI,4)/
     .                   CQ(ALPHAS(AM(4),2)/PI,4)
       YMSB(6) = YMSB(5)*CQ(ALPHAS(AM(6),2)/PI,5)/
     .                   CQ(ALPHAS(AM(5),2)/PI,5)
      ELSEIF(NF.EQ.5)THEN
       YMSB(5) = XMSB
       YMSB(4) = YMSB(5)*CQ(ALPHAS(AM(4),2)/PI,4)/
     .                   CQ(ALPHAS(AM(5),2)/PI,4)
       YMSB(6) = YMSB(5)*CQ(ALPHAS(AM(6),2)/PI,5)/
     .                   CQ(ALPHAS(AM(5),2)/PI,5)
      ELSEIF(NF.EQ.6)THEN
       YMSB(6) = XMSB
       YMSB(5) = YMSB(6)*CQ(ALPHAS(AM(5),2)/PI,5)/
     .                   CQ(ALPHAS(AM(6),2)/PI,5)
       YMSB(4) = YMSB(5)*CQ(ALPHAS(AM(4),2)/PI,4)/
     .                   CQ(ALPHAS(AM(5),2)/PI,4)
      ENDIF
      IF(Q.LT.AMC)THEN
       N0=3
       Q0 = 1.D0
      ELSEIF(Q.LE.AMB)THEN
       N0=4
       Q0 = AMC
      ELSEIF(Q.LE.AMT)THEN
       N0=5
       Q0 = AMB
      ELSE
       N0=6
       Q0 = AMT
      ENDIF
      IF(NNLO.EQ.1.AND.NF.GT.3)THEN
       XKFAC = TRAN(AM(NF),0D0)/TRAN(AM(NF),XK)
      ELSE
       XKFAC = 1D0
      ENDIF
      RUNM = YMSB(N0)*CQ(ALPHAS(Q,2)/PI,N0)/
     .               CQ(ALPHAS(Q0,2)/PI,N0)
     .       * XKFAC
      RETURN
      END

      DOUBLE PRECISION FUNCTION ALPHAS(Q,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XLB(6)
      COMMON/ALSLAM/XLB1(6),XLB2(6)
      COMMON/ALS/XLAMBDA,AMC,AMB,AMT,N0
      B0(NF)=33.D0-2.D0*NF
      B1(NF)=6.D0*(153.D0-19.D0*NF)/B0(NF)**2
      ALS1(NF,X)=12.D0*PI/(B0(NF)*DLOG(X**2/XLB(NF)**2))
      ALS2(NF,X)=12.D0*PI/(B0(NF)*DLOG(X**2/XLB(NF)**2))
     .          *(1.D0-B1(NF)*DLOG(DLOG(X**2/XLB(NF)**2))
     .           /DLOG(X**2/XLB(NF)**2))
      PI=4.D0*DATAN(1.D0)
      IF(N.EQ.1)THEN
       DO 1 I=1,6
        XLB(I)=XLB1(I)
1      CONTINUE
      ELSE
       DO 2 I=1,6
        XLB(I)=XLB2(I)
2      CONTINUE
      ENDIF
      IF(Q.LT.AMC)THEN
       NF=3
      ELSEIF(Q.LE.AMB)THEN
       NF=4
      ELSEIF(Q.LE.AMT)THEN
       NF=5
      ELSE
       NF=6
      ENDIF
      IF(N.EQ.1)THEN
        ALPHAS=ALS1(NF,Q)
      ELSE
        ALPHAS=ALS2(NF,Q)
      ENDIF
      RETURN
      END

      SUBROUTINE ALSINI(ACC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XLB(6)
      COMMON/ALSLAM/XLB1(6),XLB2(6)
      COMMON/ALS/XLAMBDA,AMC,AMB,AMT,N0
      PI=4.D0*DATAN(1.D0)
      XLB1(1)=0D0
      XLB1(2)=0D0
      XLB2(1)=0D0
      XLB2(2)=0D0
      IF(N0.EQ.3)THEN
       XLB(3)=XLAMBDA
       XLB(4)=XLB(3)*(XLB(3)/AMC)**(2.D0/25.D0)
       XLB(5)=XLB(4)*(XLB(4)/AMB)**(2.D0/23.D0)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2.D0/21.D0)
      ELSEIF(N0.EQ.4)THEN
       XLB(4)=XLAMBDA
       XLB(5)=XLB(4)*(XLB(4)/AMB)**(2.D0/23.D0)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2.D0/27.D0)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2.D0/21.D0)
      ELSEIF(N0.EQ.5)THEN
       XLB(5)=XLAMBDA
       XLB(4)=XLB(5)*(XLB(5)/AMB)**(-2.D0/25.D0)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2.D0/27.D0)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2.D0/21.D0)
      ELSEIF(N0.EQ.6)THEN
       XLB(6)=XLAMBDA
       XLB(5)=XLB(6)*(XLB(6)/AMT)**(-2.D0/23.D0)
       XLB(4)=XLB(5)*(XLB(5)/AMB)**(-2.D0/25.D0)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2.D0/27.D0)
      ENDIF
      DO 1 I=1,6
       XLB1(I)=XLB(I)
1     CONTINUE
      IF(N0.EQ.3)THEN
       XLB(3)=XLAMBDA
       XLB(4)=XLB(3)*(XLB(3)/AMC)**(2.D0/25.D0)
     .             *(2.D0*DLOG(AMC/XLB(3)))**(-107.D0/1875.D0)
       XLB(4)=XITER(AMC,XLB(3),3,XLB(4),4,ACC)
       XLB(5)=XLB(4)*(XLB(4)/AMB)**(2.D0/23.D0)
     .             *(2.D0*DLOG(AMB/XLB(4)))**(-963.D0/13225.D0)
       XLB(5)=XITER(AMB,XLB(4),4,XLB(5),5,ACC)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2.D0/21.D0)
     .            *(2.D0*DLOG(AMT/XLB(5)))**(-321.D0/3381.D0)
       XLB(6)=XITER(AMT,XLB(5),5,XLB(6),6,ACC)
      ELSEIF(N0.EQ.4)THEN
       XLB(4)=XLAMBDA
       XLB(5)=XLB(4)*(XLB(4)/AMB)**(2.D0/23.D0)
     .             *(2.D0*DLOG(AMB/XLB(4)))**(-963.D0/13225.D0)
       XLB(5)=XITER(AMB,XLB(4),4,XLB(5),5,ACC)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2.D0/27.D0)
     .             *(2.D0*DLOG(AMC/XLB(4)))**(107.D0/2025.D0)
       XLB(3)=XITER(AMC,XLB(4),4,XLB(3),3,ACC)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2.D0/21.D0)
     .            *(2.D0*DLOG(AMT/XLB(5)))**(-321.D0/3381.D0)
       XLB(6)=XITER(AMT,XLB(5),5,XLB(6),6,ACC)
      ELSEIF(N0.EQ.5)THEN
       XLB(5)=XLAMBDA
       XLB(4)=XLB(5)*(XLB(5)/AMB)**(-2.D0/25.D0)
     .             *(2.D0*DLOG(AMB/XLB(5)))**(963.D0/14375.D0)
       XLB(4)=XITER(AMB,XLB(5),5,XLB(4),4,ACC)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2.D0/27.D0)
     .             *(2.D0*DLOG(AMC/XLB(4)))**(107.D0/2025.D0)
       XLB(3)=XITER(AMC,XLB(4),4,XLB(3),3,ACC)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2.D0/21.D0)
     .            *(2.D0*DLOG(AMT/XLB(5)))**(-321.D0/3381.D0)
       XLB(6)=XITER(AMT,XLB(5),5,XLB(6),6,ACC)
      ELSEIF(N0.EQ.6)THEN
       XLB(6)=XLAMBDA
       XLB(5)=XLB(6)*(XLB(6)/AMT)**(-2.D0/23.D0)
     .            *(2.D0*DLOG(AMT/XLB(6)))**(321.D0/3703.D0)
       XLB(5)=XITER(AMT,XLB(6),6,XLB(5),5,ACC)
       XLB(4)=XLB(5)*(XLB(5)/AMB)**(-2.D0/25.D0)
     .             *(2.D0*DLOG(AMB/XLB(5)))**(963.D0/14375.D0)
       XLB(4)=XITER(AMB,XLB(5),5,XLB(4),4,ACC)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2.D0/27.D0)
     .             *(2.D0*DLOG(AMC/XLB(4)))**(107.D0/2025.D0)
       XLB(3)=XITER(AMC,XLB(4),4,XLB(3),3,ACC)
      ENDIF
      DO 2 I=1,6
       XLB2(I)=XLB(I)
2     CONTINUE
      RETURN
      END

      DOUBLE PRECISION FUNCTION XITER(Q,XLB1,NF1,XLB,NF2,ACC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      B0(NF)=33.D0-2.D0*NF
      B1(NF)=6.D0*(153.D0-19.D0*NF)/B0(NF)**2
      ALS2(NF,X,XLB)=12.D0*PI/(B0(NF)*DLOG(X**2/XLB**2))
     .              *(1.D0-B1(NF)*DLOG(DLOG(X**2/XLB**2))
     .              /DLOG(X**2/XLB**2))
      AA(NF)=12D0*PI/B0(NF)
      BB(NF)=B1(NF)/AA(NF)
      XIT(A,B,X)=A/2.D0*(1D0+DSQRT(1D0-4D0*B*DLOG(X)))
      PI=4.D0*DATAN(1.D0)
      XLB2=XLB
      II=0
1     II=II+1
      X=DLOG(Q**2/XLB2**2)
      ALP=ALS2(NF1,Q,XLB1)
      A=AA(NF2)/ALP
      B=BB(NF2)*ALP
      XX=XIT(A,B,X)
      XLB2=Q*DEXP(-XX/2.D0)
      Y1=ALS2(NF1,Q,XLB1)
      Y2=ALS2(NF2,Q,XLB2)
      DY=DABS(Y2-Y1)/Y1
      IF(DY.GE.ACC) GOTO 1
      XITER=XLB2
      RETURN
      END

      DOUBLE PRECISION FUNCTION FINT(Z,XX,YY)
C--ONE-DIMENSIONAL CUBIC INTERPOLATION
C--Z  = WANTED POINT
C--XX = ARRAY OF 4 DISCRETE X-VALUES AROUND Z
C--YY = ARRAY OF 4 DISCRETE FUNCTION-VALUES AROUND Z
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XX(4),YY(4)
      X = DLOG(Z)
      X0=DLOG(XX(1))
      X1=DLOG(XX(2))
      X2=DLOG(XX(3))
      X3=DLOG(XX(4))
      Y0=DLOG(YY(1))
      Y1=DLOG(YY(2))
      Y2=DLOG(YY(3))
      Y3=DLOG(YY(4))
      A0=(X-X1)*(X-X2)*(X-X3)/(X0-X1)/(X0-X2)/(X0-X3)
      A1=(X-X0)*(X-X2)*(X-X3)/(X1-X0)/(X1-X2)/(X1-X3)
      A2=(X-X0)*(X-X1)*(X-X3)/(X2-X0)/(X2-X1)/(X2-X3)
      A3=(X-X0)*(X-X1)*(X-X2)/(X3-X0)/(X3-X1)/(X3-X2)
      FINT=DEXP(A0*Y0+A1*Y1+A2*Y2+A3*Y3)
      RETURN
      END

      DOUBLE PRECISION FUNCTION SP(X)
C--REAL DILOGARITHM (SPENCE-FUNCTION)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMPLEX*16 CX,LI2
      CX = DCMPLX(X,0.D0)
      SP = DREAL(LI2(CX))
      RETURN
      END
 
      COMPLEX FUNCTION LI2*16(X)
C--COMPLEX DILOGARITHM (SPENCE-FUNCTION)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMPLEX*16 X,Y,CLI2
      COMMON/CONST/ZETA2,ZETA3
      ZERO=1.D-16
      XR=DREAL(X)
      XI=DIMAG(X)
      R2=XR*XR+XI*XI
      LI2=0
      IF(R2.LE.ZERO)THEN
        LI2=X
        RETURN
      ENDIF
      RR=XR/R2
      IF(R2.EQ.1.D0.AND.XI.EQ.0.D0)THEN
        IF(XR.EQ.1.D0)THEN
          LI2=DCMPLX(ZETA2)
        ELSE
          LI2=-DCMPLX(ZETA2/2.D0)
        ENDIF
        RETURN
      ELSEIF(R2.GT.1.D0.AND.RR.GT.0.5D0)THEN
        Y=(X-1.D0)/X
        LI2=CLI2(Y)+ZETA2-CDLOG(X)*CDLOG(1.D0-X)+0.5D0*CDLOG(X)**2
        RETURN
      ELSEIF(R2.GT.1.D0.AND.RR.LE.0.5D0)THEN
        Y=1.D0/X
        LI2=-CLI2(Y)-ZETA2-0.5D0*CDLOG(-X)**2
        RETURN
      ELSEIF(R2.LE.1.D0.AND.XR.GT.0.5D0)THEN
        Y=1.D0-X
        LI2=-CLI2(Y)+ZETA2-CDLOG(X)*CDLOG(1.D0-X)
       RETURN
      ELSEIF(R2.LE.1.D0.AND.XR.LE.0.5D0)THEN
        Y=X
        LI2=CLI2(Y)
        RETURN
      ENDIF
      END
 
      COMPLEX FUNCTION CLI2*16(X)
C--TAYLOR-EXPANSION FOR COMPLEX DILOGARITHM (SPENCE-FUNCTION)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMPLEX*16 X,Z
      COMMON/BERNOULLI/B2(18),B12(18),B3(18)
      COMMON/POLY/NBER
      N=NBER-1
      Z=-CDLOG(1.D0-X)
      CLI2=B2(NBER)
      DO 111 I=N,1,-1
        CLI2=Z*CLI2+B2(I)
111   CONTINUE
      CLI2=Z**2*CLI2+Z
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION FACULT(N)
C--DOUBLE PRECISION VERSION OF FACULTY
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      FACULT=1.D0
      IF(N.EQ.0)RETURN
      DO 999 I=1,N
        FACULT=FACULT*DFLOAT(I)
999   CONTINUE
      RETURN
      END
 
      SUBROUTINE BERNINI(N)
C--INITIALIZATION OF COEFFICIENTS FOR POLYLOGARITHMS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION B(18),PB(19)
      COMMON/BERNOULLI/B2(18),B12(18),B3(18)
      COMMON/CONST/ZETA2,ZETA3
      COMMON/POLY/NBER
 
      NBER=N
      PI=4.D0*DATAN(1.D0)
 
      B(1)=-1.D0/2.D0
      B(2)=1.D0/6.D0
      B(3)=0.D0
      B(4)=-1.D0/30.D0
      B(5)=0.D0
      B(6)=1.D0/42.D0
      B(7)=0.D0
      B(8)=-1.D0/30.D0
      B(9)=0.D0
      B(10)=5.D0/66.D0
      B(11)=0.D0
      B(12)=-691.D0/2730.D0
      B(13)=0.D0
      B(14)=7.D0/6.D0
      B(15)=0.D0
      B(16)=-3617.D0/510.D0
      B(17)=0.D0
      B(18)=43867.D0/798.D0
      ZETA2=PI**2/6.D0
      ZETA3=1.202056903159594D0
 
      DO 995 I=1,18
        B2(I)=B(I)/FACULT(I+1)
        B12(I)=DFLOAT(I+1)/FACULT(I+2)*B(I)/2.D0
        PB(I+1)=B(I)
        B3(I)=0.D0
995   CONTINUE
      PB(1)=1.D0
      DO 996 I=1,18
      DO 996 J=0,I
        B3(I)=B3(I)+PB(J+1)*PB(I-J+1)/FACULT(I-J)/FACULT(J+1)
     .                                            /DFLOAT(I+1)
996   CONTINUE
 
      RETURN
      END

      DOUBLE PRECISION FUNCTION QQINT(RAT,H1,H2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      N = 2
      QQINT = RAT**N * H1 + (1-RAT**N) * H2
      RETURN
      END

      DOUBLE PRECISION FUNCTION XITLA(NO,ALP,ACC)
C--ITERATION ROUTINE TO DETERMINE IMPROVED LAMBDA'S
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/PARAM/GF,ALPH,AMTAU,AMMUON,AMZ,AMW
      B0(NF)=33.D0-2.D0*NF
      B1(NF)=6.D0*(153.D0-19.D0*NF)/B0(NF)**2
      ALS2(NF,X,XLB)=12.D0*PI/(B0(NF)*DLOG(X**2/XLB**2))
     .              *(1.D0-B1(NF)*DLOG(DLOG(X**2/XLB**2))
     .              /DLOG(X**2/XLB**2))
      AA(NF)=12D0*PI/B0(NF)
      BB(NF)=B1(NF)/AA(NF)
      XIT(A,B,X)=A/2.D0*(1D0+DSQRT(1D0-4D0*B*DLOG(X)))
      PI=4.D0*DATAN(1.D0)
      NF=5
      Q=AMZ
      XLB=Q*DEXP(-AA(NF)/ALP/2.D0)
      IF(NO.EQ.1)GOTO 111
      II=0
1     II=II+1
      X=DLOG(Q**2/XLB**2)
      A=AA(NF)/ALP
      B=BB(NF)*ALP
      XX=XIT(A,B,X)
      XLB=Q*DEXP(-XX/2.D0)
      Y1=ALP
      Y2=ALS2(NF,Q,XLB)
      DY=DABS(Y2-Y1)/Y1
      IF(DY.GE.ACC) GOTO 1
111   XITLA=XLB
      RETURN
      END
