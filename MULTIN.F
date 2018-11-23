
C       PROGRAM I-MULTI-PHONON, MODIFIED MAY 1987 FROM P-MOLE OF K.R. RAO
C       THE SUBROUTINE TOW IS REWRITTEN AND CHANGES MADE IN ALL ROUTINES.
C       S.L. CHAPLOT, MAY 1987. FEB. 1991.
C
C       GIVES S(Q,W) AND DENSITY OF STATES INCLUDING MULTI-PHONON
C       AND WEIGHTED FOR NEUTRON SCATTERING.
      IMPLICIT REAL*8(A-H,O-Z)
        DIMENSION TITLE(20),G0(990),G1(990),G2(990), G11(3,990),
     1  G10(5,990),G1N(990),GPN(990),GP(990),G16(990),
     1  GTOT(990),GBROD1(990),GBROD2(990),GTNEW(990),GNEW(9,990)
        DIMENSION FWH(5),NNI(5),FNORMS(99),U2XS(99),ANGLES(99),QQS(99)
        DIMENSION JSUM(100)
c       DIMENSION JSUM(59)
        COMMON/B/ANGLE,QQ,EIFIX,EFFIX,FNORM,U2X
        COMMON/B1/S1,E3(990),COT(990),DOW2(990),IMODE,NORM,TEMP,DELE
c        COMMON/DFR/DFREQ(281,59),ISUM(59),AMI(59),B2I(59),FSTEP,FFWHM
c       COMMON/DFR/DFREQ(380,100),ISUM(100),AMI(100),B2I(100),FSTEP,FFWHM
        COMMON/DFR/DFREQ(700,100),ISUM(100),AMI(100),B2I(100),FSTEP,FFWHM
        open (5,FILE='LUMU')
        open (6,FILE='o-lumun')
        open (7,FILE='MULTI-IN')
c        open (8,FILE='z-wzfe5')

  1     FORMAT (20A4)
  2     FORMAT (8F10.4)
  3     FORMAT (10F8.4)
  4     FORMAT (20I4)
  5     FORMAT (20F4.0)
  6     FORMAT (10F8.2)
  78    FORMAT ('--------------------------------------------------')
 178    FORMAT ('DOS-----------------------------------------------')   
  79    FORMAT ('      E        G1        G2        G3        G4
     *  G5        G6')
  778   FORMAT (F10.2,7F10.4,2F12.6)
  777   CONTINUE
 7777   FORMAT(I5)

 100    READ  (5,  1,END=1000)TITLE
        PRINT 1, TITLE

C       GIVE FIXED VALUES OF EI OR EF FOR CALCULATION OF SCATTERING VECTOR.
C       GIVE ZERO FOR VARIABLE EI OR EF. ANGLE IS SCATTERING ANGLE.
C       DOS=0 OUTPUT IS S(Q,W).  DOS=1 OUTPUT IS WEIGHTED DENSITY OF STATES.
        READ (5,*)EIFIX,EFFIX,DOS
        IF(EIFIX.EQ.0.0.AND.EFFIX.EQ.0.0)GO TO 1000
        PRINT 2, EIFIX,EFFIX,DOS
        READ (5,*)IAA,IQQ
        PRINT 4,IAA,IQQ
        IF(IAA.GT.0)IQQ=0
        IF(IAA.GT.0)READ (5,*)(ANGLES(I),I=1,IAA)
        IF(IQQ.GT.0)READ (5,*)(QQS(I),I=1,IQQ)
c        IF(IAA.GT.0)PRINT 2,(ANGLES(I),I=1,IAA)
        IF(IQQ.GT.0)PRINT 2,(QQS(I),I=1,IQQ)
C       GIVE EMIN=0, AND DELE AS POSITIVE QUANTITY.
C       NUMBER OF POINTS= EN RESTRICTED TO 124.
C       EMIN AND EMAX ARE ENERGY TRANSFERS AND NOT NEUTRON ENERGIES
C       AMODE=-1.0  IF ENERGY LOSS,=1.0 IF ENERGY GAIN
C       SIGNAL ZERO MEANS USE GIVEN TWOW
C       SIGNAL NONZERO MEANS USE COMPUTED TWOW
        READ  (5,*)EMIN,EN,DELE,AMODE,SIGNAL,TWOW
        EMIN=0.0
        PRINT 2, EMIN,EN,DELE,AMODE,SIGNAL,TWOW
        IMODE=-1
        IF(AMODE.GT.0.0)IMODE=1
          N=EN+0.5
C       EMAX=EMIN+(N-1)*DELE
        IF(N.GT.199)N=199
        NN2=2*N-1

C       ATOM= NUMBER OF PARTIAL DENSITIES.
C       TEMP= TEMPERATURE.
        READ (5,*) NATOM,TEMP
        PRINT 4, NATOM
        PRINT 2, TEMP

C       B2I= B**2 FOR I-TH PARTIAL DENSITY.
        READ(5,3)(B2I(I),I=1,NATOM)
        PRINT 3,(B2I(I),I=1,NATOM)

C       AMI= MASS  FOR I-TH PARTIAL DENSITY.
        READ(5,3)(AMI(I),I=1,NATOM)
        PRINT 3,(AMI(I),I=1,NATOM)

C       ISUM= WEIGHTAGE  FOR I-TH PARTIAL DENSITY,
C       e.g., NUMBER OF ATOMS OF I-TH TYPE IN UNIT CELL.
C       JSUM= ATOM NUMBER OF I-TH PARTIAL DENSITY.
        READ (5,*)(ISUM(I),I=1,NATOM)
        READ (5,*)(JSUM(I),I=1,NATOM)
        PRINT 4,(ISUM(I),I=1,NATOM)
        PRINT 4,(JSUM(I),I=1,NATOM)
        READ (5,5)(FNORMS(I),I=1,NATOM)
        READ (5,5)(U2XS(I),I=1,NATOM)
        PRINT 2,(FNORMS(I),I=1,NATOM)
        PRINT 2,(U2XS(I),I=1,NATOM)
        DO 2224 I=1,NATOM
 2224   IF(FNORMS(I).EQ.0.0)FNORMS(I)=1.0
C       IMAX= NUMBER OF GAUSSIANS FOR BROADENINGS  OF SPECTRA.
C       FWH= FWHM OF THE I-TH GAUSSIAN IN MILLI-EV.
        READ (5,*) IMAX
c        PRINT 4, IMAX
        IF(IMAX.EQ.0)IMAX=1
        READ(5,*)(FWH(I),I=1,IMAX)
        IF(IAA.GT.1.OR.IQQ.GT.1)IMAX=1
c        PRINT 2,(FWH(I),I=1,IMAX)
        WRITE(6,   1)TITLE
        WRITE (6,2) EIFIX,EFFIX,DOS
        WRITE(6,4)IAA,IQQ
        IF(IAA.GT.0)WRITE(6,2)(ANGLES(I),I=1,IAA)
        IF(IQQ.GT.0)WRITE(6,2)(QQS(I),I=1,IQQ)
        WRITE(6,2)EMIN,EN,DELE,AMODE,SIGNAL,TWOW
        WRITE(6,4) NATOM
        WRITE(6,2) TEMP
        WRITE(6,3)(B2I(I),I=1,NATOM)
        WRITE(6,3)(AMI(I),I=1,NATOM)
        WRITE(6,4)(ISUM(I),I=1,NATOM)
        WRITE(6,4)(JSUM(I),I=1,NATOM)
        WRITE(6,3)(FNORMS(I),I=1,NATOM)
        WRITE(6,3)(U2XS(I),I=1,NATOM)
        WRITE(6,4) IMAX
        WRITE(6,2)(FWH(I),I=1,IMAX)
c        REWIND 10
!        READ(7,*)IQN,ITWT,IWT,IH,IK,IL
        WRITE(*,7)IQN,ITWT,IWT,IH,IK,IL
c        WRITE(8,7)IQN,ITWT,IWT,IH,IK,IL
7       FORMAT(16I5)
C       WRITE(7,*)IQN,ITWT,IWT,IH,IK,IL
!        READ(7,8)VOL,VOL0,PS,PMIN
c        WRITE(8,8)VOL,VOL0,PS,PMIN
8       FORMAT(6f12.6)
C       WRITE(7,*)VOL,VOL0,PS,PMIN
        READ(7,9)DFREQ
c        WRITE(8,9)DFREQ
c         9       FORMAT(8F10.6)
9       FORMAT(8F10.6)
C       WRITE(7,17)DFREQ
  17  FORMAT(4F7.4)

C       READ(5,*)((DFREQ(I,IJJ),I=1,N),IJJ=1,NATOM)
C       PRINT *,((DFREQ(I,IJJ),I=1,N),IJJ=1,NATOM)
C       WRITE(6,*)((DFREQ(I,IJJ),I=1,N),IJJ=1,NATOM)
       DO 21 IJJ=1,NATOM
       IF(ISUM(IJJ).EQ.0)GO TO 21
       J=JSUM(IJJ)
       DO 22 IKK=IJJ+1,NATOM
       IF(ISUM(IKK).EQ.0)GO TO 22
       IF(J.NE.JSUM(IKK))GO TO 22
       ISUM(IJJ)=ISUM(IJJ)+ISUM(IKK)
       ISUM(IKK)=0
       DO 23 I=1,N
 23    DFREQ(I,IJJ)=DFREQ(I,IJJ)+DFREQ(I,IKK)
 22    CONTINUE
 21    CONTINUE

        ET=2.0*TEMP/11.605
C       ET IN MILLI-EV
        DO   11   I=1,990
        E1=EMIN+DELE*(I-1)
        E1=E1*IMODE
        E3(I)=E1
        COT(I)=0.0
        IF(ABS(E1).GT.1.0D-6) COT(I)=COTH(E1/ET)
  11    CONTINUE

        DO 9001 I=1,990
        DO 9001 J=1,3
 9001   G11(J,I)=0.0
        IF(IAA.GT.0)IRR=IAA
        IF(IQQ.GT.0)IRR=IQQ
        DO 9000 IRRR=1,IRR
c        PRINT  78
c        PRINT  78
        ANGLE=0.0
        QQ=0.0
        IF(IAA.GT.0)ANGLE=ANGLES(IRRR)
        IF(IQQ.GT.0)QQ=QQS(IRRR)
        DO 2222 I=1,990
        GTOT(I)=0.0
        G16(I)=0.0
        GP(I)=0.0
        GPN(I)=0.0
        G0(I)=0.0
 2222   CONTINUE

        NNI(1)=2*(N-1)+1
        NNI(2)=3*(N-1)+1
        NNI(3)=4*(N-1)+1
        NNI(4)=3*(N-1)+1
        NNI(5)=2*(N-1)+1

        IKK=0
        DO 1222 IJJ=1,NATOM
       IF(ISUM(IJJ).EQ.0)GO TO 1222
        IKK=IKK+1
        IF(IKK.GT.12)STOP 12
        BATMSQ=B2I(IJJ)*ISUM(IJJ)
        AMATOM=AMI(IJJ)
c        PRINT 2,BATMSQ,AMATOM
C       READ G0(W)= F(W),  ONE-PHONON D.O.S.
C       G(W)=G1(W) IS  NORMALISED TO UNITY BET. -INF TO INF
C        READ (5,*)(G0(I),I=1,N)
        DO 2223 I=1,N
 2223   G0(I)= DFREQ(I,IJJ)
        FNORM=FNORMS(IJJ)
        U2X=U2XS(IJJ)
        WRITE(6,2)FNORM,U2X
        G0(1)=0.0
        WRITE(6,2)(G0(I),I=1,N)
        SUM=0.0
        DO 226 I=1,N
        SUM=SUM+G0(I)
 226    CONTINUE
        SUM=SUM*DELE
        P=FNORM/SUM
        DO 227 I=1,N
        G0(I)=G0(I)*P
 227    CONTINUE

        DO   10   I=2,N
        E1=E3(I)
        POP1=COT(I)-1.0
        G1(I)=G0(I)*POP1/E1
  10    CONTINUE
        G1(1)=0.5*(G1(2)+G1(3))

        WRITE(6, 78)
        NORM=0
        CALL    GW(G1,G1N,N)
        SG1=S1
        IF(SIGNAL.NE.0.0) CALL    TOW(G0,N,AMATOM)
        DO   125   I=1,N
        GP(I)=G1(I)
        GPN(I)=G1N(I)
  125   CONTINUE

        NN=N
        DO   122  IJ=1,5
        NNP=NN
        NN=NNI(IJ)
        IF(NN.GT.990)STOP 990
        DO   12   I=1,NN
        GSUM=0.0

        DO   13   J=1,NN2
C          W1=-(N-1) TO N-1;   G2(W)= INTEGRAL[ G1(W1)GP(W-W1) DW1 ]
           N1=J-N
        N1=N1*IMODE
        IF(N1.LT.0)FF=G1N(-N1+1)
        IF(N1.GE.0)FF=G1(N1+1)
C                                M=N
C       FOR IMODE=-1,  W=0, W1=(-M, M); W=-4M, W1=(-M, M)
C        W2=    W - W1;     W2=( M,-M);           W2=(-3M,-5M)
C       FOR IMODE= 1,  W=0, W1=(-M, M); W= 4M, W1=(-M, M)
C        W2=    W - W1;     W2=( M,-M);           W2=( 5M, 3M)
        N2=IMODE*(I-1)-(J-N)
        N2=N2*IMODE
        IF(N2)669,670,670
  669   SF=GPN(-N2+1)
        GO    TO   671
  670   IF(N2.GE.NNP)GO TO 13
        SF=GP(N2+1)
  671   GSUM=GSUM+FF*SF
  13    CONTINUE
        G2(I)=GSUM *DELE
  12    CONTINUE

        DO   20   I=1,NN
        GP(I)=G2(I)
        G10(IJ,I)=G2(I)
  20    CONTINUE
        CALL    GW(GP,GPN,NN)
  122   CONTINUE

        WRITE(6,   79)
        DO   107    I=1,NN2
        IF(ABS(SIGNAL).GT.1.0D-6) TWOW=DOW2(I)
        IF(TWOW.EQ.0.0)GO TO 107
        NM=I

        X1=   TWOW*(G1(I)+         (TWOW/2.0)*(G10(1,I)+
     1      (TWOW/3.0)*(G10(2,I)+  (TWOW/4.0)*(G10(3,I)+
     1      (TWOW/5.0)*(G10(4,I)+  (TWOW/6.0)* G10(5,I)  )))))

        X=EXP(-TWOW)*BATMSQ
        GNEW(IKK,I)= X1*X
      IF(I-1.EQ.(I-1)/5*5)WRITE(6,778)E3(I),G1(I),(G10(J,I),J=1,5)
     * ,TWOW,X1,X
        GTOT(I)=GTOT(I)+G1(I)*TWOW*X
        IF(I.GT.1)G16(I)=G16(I)+ X*TWOW*(COT(I)-1.0)/(E3(I)*SG1)
  107   CONTINUE
        WRITE(6,   78)
 1222   CONTINUE
        MATOM=IKK

        DO 1223 I=1,NM
        GTNEW(I)=0.0
        DO 1224 J=1,MATOM
        GTNEW(I)=GNEW(J,I)+GTNEW(I)
 1224   CONTINUE
        GTOT(I)=GTOT(I)*1000.0
        GTNEW(I)=GTNEW(I)*1000.0
 1223   CONTINUE

        DO 2226 K=1,IMAX
        FWHM=FWH(K)
        CALL BROAD(GTNEW,GBROD2,FWHM,NM)
        CALL BROAD(GTOT,GBROD1,FWHM,NM)
        IF(ABS(DOS).LT.1.0D-6)WRITE(6,15)
        IF(ABS(DOS).GE.1.0D-6)WRITE(6,16)
 15     FORMAT(20X,' S(Q,W)')
 16     FORMAT(20X,' F(W) - DENSITY OF STATES')
        WRITE(6, 1778)
 1778   FORMAT(5X,'THZ',5X,'MEV ONE-PHONON  MULTI-PH   DIFF',
     * 11X,'BROADENED         S(Q,W)/F(W)')
        DO 225 I=1,NM
        X=-E3(I)/4.136
        XT=GTOT(I)
        XN=GTNEW(I)
        IF(ABS(DOS).LT.1.0D-6)GO TO 1228
        IF(I.EQ.1)GO TO 225
        XT=XT/G16(I)
        XN=XN/G16(I)
        GBROD1(I)=GBROD1(I)/G16(I)
        GBROD2(I)=GBROD2(I)/G16(I)
 1228   CONTINUE
        X1=XN-XT
        X2=GBROD2(I)-GBROD1(I)
        G11(1,I)=G11(1,I)+GBROD1(I)
        G11(2,I)=G11(2,I)+GBROD2(I)
        G11(3,I)=G11(3,I)+X2
        WRITE (6,14)X,E3(I),XT,XN,X1,GBROD1(I),GBROD2(I),X2,G16(I)
  14    FORMAT(2F8.2,2(2F9.3,F10.4),F10.5)
 225   CONTINUE
c        PRINT  78
c        PRINT 6,(GBROD1(I),I=1,NM)
c        PRINT  78
c        PRINT 6,(GBROD2(I),I=1,NM)
c        PRINT  78
c        PRINT 6,(GBROD2(I)-GBROD1(I),I=1,NM)
 2226   CONTINUE
 9000   CONTINUE
        X=0.0
        DO 9003 I=1,NM
 9003   IF(X.LT.G11(2,I)) X=G11(2,I)
        X=100.0/X
        DO 9004 I=1,NM
        DO 9004 J=1,3
 9004   G11(J,I)=G11(J,I)*X
        J=1.0/DELE +0.1
c        PRINT  78
c        PRINT 4,J
c        PRINT 3,X
c        PRINT 6,(G11(1,I),I=1,NM,J)
c        PRINT  78
c        PRINT 6,(G11(2,I),I=1,NM,J)
c        PRINT  78
c        PRINT 6,(G11(3,I),I=1,NM,J)
c        PRINT  78
        WRITE (6,178)
        WRITE (6,3)X
        WRITE (6,6)(G11(1,I),I=1,NM)
        WRITE (6,178)
        WRITE (6,6)(G11(2,I),I=1,NM)
        WRITE (6,178)
        WRITE (6,6)(G11(3,I),I=1,NM)
        WRITE (6,78)
        GO TO 100
 1000   STOP
        END
C
C------------------------------------------------------------------
        SUBROUTINE    GW (A,B,N)
C       CALCULATION  OF G(-W) GIVEN G(W) AND
C       NORMALISING G(W)  BET.-INF.AND. INF.
      IMPLICIT REAL*8(A-H,O-Z)
        COMMON/B1/S1,E3(990),COT(990),DOW2(990),IMODE,NORM,TEMP,DELE
        DIMENSION  A(990),B(990)
 6666  FORMAT(I5,4F10.4)
 6667  FORMAT ('    I    E     G(E)    E-     G(-E)    ')
C       ET=2.0*TEMP/11.605
C       ET IN MILLI-EV
        B(1)=A(1)
        SUMT=A(1)+B(1)
        DO   10   I=2,N
C       E1=DELE*(I-1)
        E1=E3(I)
C       POP1=COTH(E1/2KT)-1;  POP2=COTH(-E1/2KT)-1.
        POP1=COT(I)-1.0
        POP2=-COT(I)-1.0
        B(I)=-A(I)*POP2/POP1
        SUMT=SUMT+A(I)+B(I)
  10    CONTINUE
        SUMT=SUMT*DELE
        S1=SUMT
        IF(NORM.NE.0)GO TO 35
        DO 25 I=1,N
        A(I)=A(I)/SUMT
        B(I)=B(I)/SUMT
  25    CONTINUE
  35    CONTINUE
        NORM=1
        RETURN
        END
C
C------------------------------------------------------------------
        SUBROUTINE TOW(G0,N,AMATOM)
      IMPLICIT REAL*8(A-H,O-Z)
        DIMENSION G0(990)
        COMMON/B/ANGLE,QQ,EIFIX,EFFIX,FNORM,U2X
        COMMON/B1/S1,E3(990),COT(990),DOW2(990),IMODE,NORM,TEMP,DELE

        U2=U2X
        IF(U2X.NE.0.0)GO TO 11
C       ADD TO SUM FOR I=1 TERM.
         SUM=0.0
C        SUM=0.25*(G0(2)*COT(2)/E3(2) + G0(3)*COT(3)/E3(3))
        SUMD=G0(1)
        DO 10 I=2,N
        SUM=SUM+ G0(I)*COT(I)/E3(I)
        SUMD=SUMD+ G0(I)
 10     CONTINUE
        U2=FNORM*SUM*0.5*4.177/(SUMD*AMATOM)
 11     CONTINUE

  2     FORMAT(8F10.4)
  3     FORMAT(5X,'THZ',7X,'MEV',7X,'Q2',8X,'U2   EXP(-Q2*U2)')
       WRITE(6,3)

      IF(ANGLE.NE.0.0)COSTH=COS(ANGLE*3.141593/180.0)

C     UNITS MILLIEV.
      NN2=2*N-1
      EMIN=0.0
      DO 100 I=1,NN2
        E1=EMIN+DELE*(I-1)
C     E1, DE ARE NEGATIVE FOR NEUTRON ENERGY LOSS. DE=EF-EI
      DE=E1*IMODE
      DET=DE/4.136
      EI=EIFIX
      EF=EFFIX
C     IF EF IS FIXED.
      IF(ABS(EIFIX).LT.1.0D-6) EI=EFFIX-DE
C     IF EI IS FIXED.
      IF(ABS(EFFIX).LT.1.0D-6) EF=EIFIX+DE
      DOW2(I)=0.0
      IF(EF*EI.LE.0.0)GO TO 100
      Q2=0.0
      IF(ANGLE.NE.0.0)Q2=EF+EI-2.0*SQRT(EF*EI)*COSTH
      Q2=Q2/2.0718
      IF(ANGLE.EQ.0.0)Q2=QQ**2
C     Q2 IS IN ANG**-2
      DOW2(I)= Q2*U2
      IF(ANGLE.EQ.0.0.AND.I.GT.1)GO TO 100
      DW=EXP(-Q2*U2)
c      IF(I.EQ.1)PRINT 2,DET,DE,Q2,U2,DW
      IF(I-1.EQ.(I-1)/5*5)WRITE (6,2)DET,DE,Q2,U2,DW
100    CONTINUE
        RETURN
        END
C
C----------------------------------------------------------------------
        SUBROUTINE BROAD(A,B,FWHM,N)
      IMPLICIT REAL*8(A-H,O-Z)
        COMMON/B1/S1,E3(990),COT(990),DOW2(990),IMODE,NORM,TEMP,DELE
C       FWHM IN MEV
        DIMENSION A(990),B(990)
 3      FORMAT('--I AM IN BROAD-----------')
c         PRINT 3
        FWHM1=FWHM/DELE
        SIG=FWHM1/(2*1.1774)
        ISIG=2*SIG+1
        DO 12 I=1,N
        B(I)=A(I)
        IF(FWHM.EQ.0.0)GO TO 12
        B(I)=0.0
        DO 10 J=-ISIG,ISIG
        K=I+J
        IF(K.LT.0)K=-K+1
        IF(K.EQ.0)K=1
        IF(K.GT.N)K=N
        B(I)=B(I)+A(K)*EXP(-(J**2)/(2*SIG**2))
  10    CONTINUE
        B(I)=B(I)/(SQRT(2.0*3.141592)*SIG)
  12    CONTINUE
   2    FORMAT(10F8.4)
C        WRITE(6,2)(B(I),I=1,200)
        RETURN
        END
C
C-------------------------------------------------------------------
      FUNCTION COTH(X)
      IMPLICIT REAL*8(A-H,O-Z)
      Y=ABS(X)/X
      IF(ABS(X).GT.30.0) GO TO 1
      Y=EXP(2.0*X)
      Y=1.0+ 2.0/(Y-1.0)
   1  COTH=Y
      RETURN
      END
C
C--------------------------------------------------------------------
      SUBROUTINE ISQW(GW)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION GW(100),SQW(100),TITLE(10)

      ANGLE=133.0 *3.141593/180.0
      COSTH=COS(ANGLE)

      TEMP=20.0/11.605
C     UNITS MILLIEV.
      EF=4.0
C10    READ (5,1,END=200)TITLE
 1    FORMAT(10A8)
      WRITE(6,1)TITLE
C      READ (5,*)GW
      WRITE(6,7)
 7    FORMAT('    THZ    MEV    GW    SQW     Q2     N     DW')
      DO 100 I=1,100
      DET=(I-1)/10.0
      IF(I.EQ.1)DET=0.025
      DE=DET*4.136

      EI=EF+DE
      Q2=EF+EI-2.0*SQRT(EF*EI)*COSTH
      Q2=Q2/2.0718
C     Q2 IS IN ANG**-2
      DW=EXP(-Q2*0.002)

      PN=EXP(DE/TEMP)-1.0
      PN=1.0+1.0/PN

      SQW(I)=GW(I)*DW*Q2*PN/DE
      WRITE(6,2)DET,DE,GW(I),SQW(I),Q2,PN,DW
 2    FORMAT(2F8.2,2F8.4,F8.2,2F8.3)
 100  CONTINUE

      WRITE(6,3)
      WRITE(6,3)SQW
 3    FORMAT(10F8.4)
C      GO TO 10
C 200  STOP
      RETURN
      END
C------------------------------------------------------------------

CC TL-2122 FEB. 1991
CC  00.0    30.0        1.0
CC  0 3
CC    5.0       6.0       7.0
CC     .0000  190.0000    0.5000   -1.0000    1.0000    0.0000
CC   23  298.0000
CC  9.7093  9.7093  9.7093  3.4636  3.4636  3.4636  7.4855  7.4855  7.4855  4.2346
CC  4.2346  4.2346  4.2346  4.2346  4.2346  4.2346  4.2346  4.2346  4.2346  4.2346
CC  4.2346  3.0172  3.0172
CC204.3800204.3800204.3800137.3300137.3300137.3300 63.5460 63.5460 63.5460 16.0000
CC 16.0000 16.0000 16.0000 16.0000 16.0000 16.0000 16.0000 16.0000 16.0000 16.0000
CC 16.0000 40.0800 40.0800
CC  4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 2
CC  1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 7 7 7 8 8
CC 0.0 0.0
CC 0.0 0.0
CC 0.0 0.0
CC 0.0 0.0
CC    1
CC    6.0000     .0000     .0000