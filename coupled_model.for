CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      MODULE TEMPARRAY_1
      IMPLICIT NONE
      DOUBLE PRECISION MODEL_COOR(10000,5),ULOCAL_TEMP(100000,3),
     1            ALOCAL_TEMP(100000,9),HH(100000,2),BASEMENT(100000),
     2            H_HISTORY(100000)
      INTEGER ELE(10000,14),TEMP_KINC(100000),TEMP_TIME(100000),
     1            TEMP_NODE(100000),COR(4),TEMP_EAR(100000,8),TEK(1),
     2            TEMP_C(100000,8)
      INTEGER MINI_NODE,MAX_NODE,NUM_NODE,LIE,OPEN_EAR
      CHARACTER*256 BASE_PATH,PATH,PATH1,PATH2,PATH3,
     1              PATH4,PATH5,PATH6,PATH7,PATH8,PATH9
      
      
      SAVE
      CONTAINS
      
      SUBROUTINE SetPath()
      BASE_PATH = 'C:/Users/Administrator/Desktop/work12/'      
ccc      BASE_PATH ='/data3/sunxinyao/CAGS/qwe2/'
      PATH = TRIM(adjustl(BASE_PATH)) // 'mo_pre.txt'
      PATH1 = TRIM(adjustl(BASE_PATH)) // 'model_param.txt'
      PATH2 = TRIM(adjustl(BASE_PATH)) // 'ada_ele.txt'
      
      PATH3 = TRIM(adjustl(BASE_PATH)) // 'u_slip.txt'
      PATH4 = TRIM(adjustl(BASE_PATH)) // 'hh.txt'
      PATH5 = TRIM(adjustl(BASE_PATH)) // 'drainage_area.txt'
      PATH6 = TRIM(adjustl(BASE_PATH)) // 'total_erosion.txt'
      PATH7 = TRIM(adjustl(BASE_PATH)) // 'erosion_rate.txt'
      PATH8 = TRIM(adjustl(BASE_PATH)) // 'catchment.txt'
      PATH9 = TRIM(adjustl(BASE_PATH)) // 'basement.txt'
      RETURN
      END SUBROUTINE SetPath
      
      
      END MODULE
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      SUBROUTINE SDVINI(STATEV,COORDS,NSTATV,NCRDS,NOEL,NPT,
     1 LAYER,KSPT)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION STATEV(NSTATV),COORDS(NCRDS)
      STATEV(1:NSTATV)=0
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  
      SUBROUTINE uexternaldb(lop,lrestart,time,dtime,kstep,kinc)
      include 'aba_param.inc' 
      dimension time(2)
      
      if (lop.eq.1) then
      call mutexinit(1)
      call mutexinit(2)
      call mutexinit(3)
      call mutexinit(4)
      call mutexinit(5)
      call mutexinit(6)
      endif
      
      RETURN
      END
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C     
      USE TEMPARRAY_1
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4)
      
      INTEGER ZRO,ONE,N_Newton,MAX_Newton,NI,NJ,IN1,C_RANGE
      DOUBLE PRECISION E,BSB,ND,C,AF,D,Q,D2,D3,D4,XGM,DSGM,
     1                 I1,J2,ALPHA,BETA,F2,DP1,DP2,DP,DP3,
     2                 DP4,DP5,DSTRESS,DLP1,DLP2,DLP3,DL,
     3                 HMO,G,AMP2,TIME_EAR,DH,I1_EXC,K2,ROU,
     4                 C_DRO,STRESS2,temp2,STRESS_K,STRESS_TRI,
     5                 DSTRAN2,TIME_TEMP,FF1,DEVI,FF2,MIN_DL,
     6                 TIME_NOEAR(2),R,AMP,DL_TEMP,STRESS_K2,FF3
      PARAMETER (ONTHIRD=0.33333,ONSIX=0.16666)
      DIMENSION D(NTENS,NTENS),Q(NTENS,NTENS),D2(NTENS,NTENS),
     2          D3(NTENS,NTENS),D4(NTENS,NTENS),XGM(NTENS,NTENS),
     3          DSGM(NTENS),F2(NTENS),G(NTENS),DP1(NTENS),
     4          DP2(NTENS,NTENS),DP(NTENS,NTENS),DP3(NTENS),DP5(NTENS),
     5          DSTRESS(NTENS),DLP1(NTENS),DLP2(NTENS),DSTRAN2(NTENS),
     6          FE(NTENS),STRESS2(NTENS),STRESS_K2(NTENS),
     7          STRESS_K(NTENS),STRESS_TRI(NTENS)
      PARAMETER (ZRO=0,ONE=1)
      
C ----------------------------------------------------------------------------------------------------------------      
C                UMAT FOR MODIFIED BURGERS MODEL
C ----------------------------------------------------------------------------------------------------------------
C    PARAMETERS OF MODIFIED BURGERS MODEL AT TRANSIENT TEMPERATURE 
C    材料参数说明：后四项输入0则为纯黏弹性
C    PROPS(1)- E    Young's modulus
C    PROPS(2)- μ   Poisson's ratio
C    PROPS(3)- R   Density    
C    PROPS(4)- η   Viscosity
C    PROPS(5)- C    Cohesion
C    PROPS(6)- Cd   Cohesion reduction  
C    PROPS(7)- β   Angle of internal friction
C    PROPS(8)- Hp  Plastic modulus
C    Three-dimensional stress and strain components (S11, S22, S33, S12, S13, S23) (E11, E22, E33, E12, E13, E23)
C ---------------------------------------------------------------------------------------------------------------
C     Initial material setup
      E=PROPS(1)      ! Young's modulus
      BSB=PROPS(2)    ! Poisson's ratio
      ROU=PROPS(3)    ! Density
      ND=PROPS(4)     ! Viscosity (ND>1E+30: pure elasticity)
      C=PROPS(5)      ! Cohesion
      C_DRO=PROPS(6)  ! Cohesion reduction
      AF=PROPS(7)     ! Angle of internal friction
      HMO=PROPS(8)    ! Plastic modulus (Stain hardening: HMO>0; Stain softening: HMO<0)
      ALPHA=2.0*SIND(AF)/(SQRT(3.0)*(3+SIND(AF)))
      BETA=6.0*C*COSD(AF)/(SQRT(3.0)*(3+SIND(AF)))
      I1_EXC=COORDS(3)*9.8*ROU*3    ! The default top coordinate is 0
      IF (I1_EXC.GT.0) I1_EXC=0
      I1_EXC=0
      
C     Initialization set
      DSTRESS=0
      DSGM=0
      DP=0
      DL=0
      DP4=0
      FF1=0
      K2=0
      FF2=0
      DH=0
      DP5=0
      STRESS2=0
      DSTRAN2=DSTRAN
      DEVI=1E-2
      MAX_Newton=100
      C_RANGE=2E+6
      
      TIME_EAR=1                 ! Earthquake occurrence time
      TIME_NOEAR(1)=(3.15E+7)*10000         ! NO Earthquake occurrence time
      TIME_NOEAR(2)=3.15E+7         ! High erosion phase after an earthquake
      
	IF ((KINC.LE.1).AND.(JSTEP(1).EQ.1)) THEN
      STATEV(1:NSTATV)=ZRO
      
      call MutexLock(3)
      OPEN_EAR=0
      TEMP_EAR(:,:)=0
      IF (SUM(ABS(TEMP_C(:,:))).LT.1) THEN
      CALL random_seed()
      DO NI = 1,100000
      DO NJ = 1,8
      CALL random_number(R)
      TEMP_C(NI,NJ)=R*(2*C_RANGE+1)-C_RANGE
      END DO
      END DO
      END IF
      call MutexUnlock(3)
      END IF
      
      IF (C_DRO.GT.1E-10) C=C+TEMP_C(NOEL,NPT)
      
C     Parameter changes at coseismic moments (e.g., decreased cohesion and shortened time steps)        
      IF (OPEN_EAR.GT.0.5) PNEWDT=TIME_NOEAR(1)/DTIME 
      
      IF ((TEMP_EAR(NOEL,NPT).GT.0.5).AND.(DTIME.LT.TIME_EAR+0.1)) THEN
      C=C-C_DRO
      BETA=6*C*COSD(AF)/(SQRT(3.0)*(3+SIND(AF)))
      E=E
      END IF
      
      IF (KINC.GT.TEK(1)) THEN
      call MutexLock(6)
      TEMP_EAR(:,:)=0
      call MutexUnlock(6)
      END IF
      
      IF (NOEL.EQ.439) THEN
      WRITE(6,*) '*******************************'
      WRITE(6,*) 'UMAT_NOEL,NPT:',NOEL,NPT
      WRITE(6,*) 'DTIME:',DTIME
      WRITE(6,*) 'UMAT_STEP,KINC:',JSTEP(1),KINC
      WRITE(6,*) 'TOTLE TIME:',TIME(2)
      WRITE(6,*) 'STRESS1:',STRESS
      WRITE(6,*) 'DSTRAN:',DSTRAN
      WRITE(6,*) 'STATEV(1:NSTATV):',STATEV(1:NSTATV)
      WRITE(6,*) 'TIME_NOEAR1,TIME_NOEAR2:',TIME_NOEAR(1:2)
      WRITE(6,*) 'TEMP_C(NOEL,NPT):',TEMP_C(NOEL,NPT)
      WRITE(6,*) 'OPEN_EAR',OPEN_EAR
      WRITE(6,*) 'COORDS(3),ROU,I1_EXC',COORDS(3),ROU,I1_EXC
      WRITE(6,*) 'ALPHA,BETA,E',ALPHA,BETA,E
      WRITE(6,*) 'PROPS(1:8)',E,BSB,ROU,ND,C,C_DRO,AF,HMO
      END IF
      
C     Viscoelastic material property matrix            
C     D=[D], Q=[Q]^(-1), D2=[D]^(-1), D4=[D~]=([D]^(-1)+[Q]^(-1)*DTIME)^(-1), DSGM={dσ~}
      CALL ELAMATRIX(D(:NTENS,:NTENS),E,BSB,NTENS,NDI)
      CALL VISMATRIX(Q(:NTENS,:NTENS),ND,NTENS,NDI)
      CALL INVERT(D(:NTENS,:NTENS),D2(:NTENS,:NTENS),NTENS)
      CALL MATRIX_ADD(D2(:NTENS,:NTENS),Q(:NTENS,:NTENS)*DTIME,
     1                                D3(:NTENS,:NTENS),NTENS)
      CALL INVERT(D3(:NTENS,:NTENS),D4(:NTENS,:NTENS),NTENS)
      CALL MATRIX_MUL(-1.0*D4(:NTENS,:NTENS),Q(:NTENS,:NTENS)*DTIME,
     1                                XGM(:NTENS,:NTENS),NTENS)
      CALL MATRIX_MUL2(XGM(:NTENS,:NTENS),STRESS(:NTENS),
     1                                    DSGM(:NTENS),NTENS)
      
C     Elastic predictor      
      IF (ND.GT.1E+30) D4(:NTENS,:NTENS)=D(:NTENS,:NTENS)
      IF (ND.GT.1E+30) DSGM(:NTENS)=0
      
      CALL STRESS_IJ(STRESS(:NTENS),NTENS,I1,J2)
      CALL FMATRIX(STRESS(:NTENS),ALPHA,SQRT(J2),F2(:NTENS),NTENS)
      CALL GMATRIX(STRESS(:NTENS),SQRT(J2),G(:NTENS),NTENS) 
      IF (ABS(J2).LT.1E-10) G=0
      IF (ABS(J2).LT.1E-10) F2=0
      CALL MATRIX_MUL4(G(:NTENS),G(:NTENS),DH,NTENS)
      
C     DP2=[D~]{dG/dσ}{dF/dσ}; DP=[D~]{dG/dσ}{dF/dσ}[D~]
      CALL MATRIX_MUL2(D4(:NTENS,:NTENS),G(:NTENS),DP1(:NTENS),NTENS)
      CALL MATRIX_MUL3(DP1(:NTENS),F2(:NTENS),DP2(:NTENS,:NTENS),NTENS)
      CALL MATRIX_MUL(DP2(:NTENS,:NTENS),D4(:NTENS,:NTENS),
     1                                        DP(:NTENS,:NTENS),NTENS)  

C     DP5={dσp~}=[D~]{dG/dσ}{dF/dσ}{dσ~}
      CALL MATRIX_MUL2(DP2(:NTENS,:NTENS),DSGM(:NTENS),
     1                                        DP5(:NTENS),NTENS)
      
C     DP4={dF/dσ}[D~]{dG/dσ}
      CALL MATRIX_MUL22(F2(:NTENS),D4(:NTENS,:NTENS),DP3(:NTENS),NTENS)
      CALL MATRIX_MUL4(DP3(:NTENS),G(:NTENS),DP4,NTENS)

C     Combination Matrix
      DP(:NTENS,:NTENS)=DP(:NTENS,:NTENS)/(DP4-HMO*DH)
      DP5(:NTENS)=DP5(:NTENS)/(DP4-HMO*DH)
      IF (ABS(J2).LT.1E-10) DP=0
      IF (ABS(J2).LT.1E-10) DP5=0
      IF (ABS(C).LT.1E-10) DP=0
      
C     Jacobian matrix      
      DDSDDE(:NTENS,:NTENS)=D4(:NTENS,:NTENS)-DP(:NTENS,:NTENS)
      
C     Test stress calculation
      CALL MATRIX_MUL2(D4(:NTENS,:NTENS),DSTRAN(:NTENS),
     1                                        DSTRESS(:NTENS),NTENS)
      STRESS2(:NTENS)=STRESS(:NTENS)+DSTRESS(:NTENS)+DSGM(:NTENS)
      CALL STRESS_IJ(STRESS2(:NTENS),NTENS,I1,J2)
CCCC      FF1=ALPHA*(I1+I1_EXC)+SQRT(J2)-(BETA+HMO*STATEV(6))
      FF1=ALPHA*(I1+I1_EXC)+SQRT(J2)-BETA
      
      IF (NOEL.EQ.439) THEN
      WRITE(6,*) 'ALPHA*I1,SQRT(J2),BETA,FF:',
     1    ALPHA*(I1+I1_EXC),SQRT(J2),BETA,FF1
      END IF
      
      IF (TEMP_EAR(NOEL,NPT).GT.0.5) GOTO 214
      IF (FF1.LE.(DEVI*100)) STRESS(:NTENS)=STRESS2(:NTENS) ! This step is elastic and does not require plastic correction
      IF (FF1.LE.(DEVI*100)) DDSDDE(:NTENS,:NTENS)=D4(:NTENS,:NTENS)
      IF (C.LT.1E-10) STRESS(:NTENS)=STRESS2(:NTENS) ! pure viso-elastic
      IF (C.LT.1E-10) DDSDDE(:NTENS,:NTENS)=D4(:NTENS,:NTENS) 
      IF ((C.LT.1E-10).OR.(FF1.LE.(DEVI*100))) GOTO 112
      
214   CONTINUE       
      
C     Drucker-Prager yield function (Newton)
      FF2=FF1
      FF3=FF1
      STRESS_TRI(:NTENS)=STRESS2(:NTENS)
      N_Newton=0
      DL=0
      MIN_DL=0
      DH=0
      AMP=1.0
      
C     Initialization of plastic correction
      CALL STRESS_IJ(STRESS_TRI(:NTENS),NTENS,I1,J2)
      CALL FMATRIX(STRESS_TRI(:NTENS),ALPHA,SQRT(J2),F2(:NTENS),NTENS)
      CALL GMATRIX(STRESS_TRI(:NTENS),SQRT(J2),G(:NTENS),NTENS)  
      
C     σK=σtri-dλ(k)*D4*G
123   N_Newton=N_Newton+1
      CALL MATRIX_MUL2(D4(:NTENS,:NTENS),G(:NTENS),
     1                    DSTRESS(:NTENS),NTENS)
      STRESS_K(:NTENS)=STRESS_TRI(:NTENS)-DL*DSTRESS(:NTENS)+DSGM
      
      CALL STRESS_IJ(STRESS_K(:NTENS),NTENS,I1,J2)
      CALL FMATRIX(STRESS_K(:NTENS),ALPHA,SQRT(J2),F2(:NTENS),NTENS)
      CALL GMATRIX(STRESS_K(:NTENS),SQRT(J2),G(:NTENS),NTENS)  
      CALL MATRIX_MUL4(G(:NTENS),G(:NTENS),DH,NTENS) ! DH=2/3*{dG/dσ}{dG/dσ})^(1/2)
      DH=SQRT(2.0/3.0*DH)
CCCC      K2=STATEV(6)+DH*DL
      FF1=ALPHA*(I1+I1_EXC)+SQRT(J2)-(BETA+HMO*K2)
      
      IF (FF1.LE.DEVI) GOTO 124
      IF (N_Newton.EQ.MAX_Newton) GOTO 124
      
C     R(dλ(k))=F(σ(k))=0
C     dR/d(dλ)=-dF/dσ*[D4]*dG/dσ
      CALL MATRIX_MUL22(F2(:NTENS),D4(:NTENS,:NTENS),DP3(:NTENS),NTENS)
      CALL MATRIX_MUL4(DP3(:NTENS),G(:NTENS),DP4,NTENS)

C     dλ(k+1)=dλ(k)-R(dλ(k))/(dR/d(dλ))
      DL_TEMP=FF1/(-DP4-HMO*DH)*AMP
      DL=DL-DL_TEMP
      
      FF2=FF1
      STRESS_K2(:NTENS)=STRESS_K(:NTENS)
      MIN_DL=DL
      
      GOTO 123
      
124   CONTINUE
      
      IF (N_Newton.EQ.MAX_Newton) then
      WRITE(6,*) 'ERROR,ERROR,ERROR'
      WRITE(6,*) 'NAN_NOEL,N,DL,FF1,FF2',NOEL,N_Newton,DL,FF1,FF2
      WRITE(6,*) 'STRESS_K(:NTENS):',STRESS_K(:NTENS)
      WRITE(6,*) 'STRESS_TRI(:NTENS):',STRESS_TRI(:NTENS)
      WRITE(6,*) 'F2(:NTENS):',F2(:NTENS)
      WRITE(6,*) 'G(:NTENS):',G(:NTENS)
      WRITE(6,*) 'D4*G:',DSTRESS(:NTENS)
      WRITE(6,*) 'DP4:',DP4
cccc      CALL XIT
      STRESS(:NTENS)=STRESS2(:NTENS)
      END IF
      
C     The iterative stress and the historical total effective plastic strain are stored
      STRESS(:NTENS)=STRESS_K(:NTENS)
CCCC      STATEV(6)=K2
      
      CALL STRESS_IJ(STRESS(:NTENS),NTENS,I1,J2)
      CALL FMATRIX(STRESS(:NTENS),ALPHA,SQRT(J2),F2(:NTENS),NTENS)
      CALL GMATRIX(STRESS(:NTENS),SQRT(J2),G(:NTENS),NTENS)  
      CALL MATRIX_MUL4(G(:NTENS),G(:NTENS),DH,NTENS) 
      DH=SQRT(2.0/3.0*DH)
CCCC      K2=STATEV(6)+DH*DL
      FF1=ALPHA*(I1+I1_EXC)+SQRT(J2)-(BETA+HMO*K2)
      
      IF (NOEL.EQ.439) THEN
      WRITE(6,*) 'N_Newton,DL,HMO*K2,FF1,FF2:',N_Newton,
     1            DL,HMO*K2,FF1,FF2
      WRITE(6,*) 'STRESS_K2(:NTENS):',STRESS_K(:NTENS)
      WRITE(6,*) 'STRESS_TRI(:NTENS):',STRESS_TRI(:NTENS)
      WRITE(6,*) 'DSGM:',DSGM
      WRITE(6,*) 'F2(:NTENS):',F2(:NTENS)
      WRITE(6,*) 'G(:NTENS):',G(:NTENS)
      WRITE(6,*) 'DP4:',DP4
      WRITE(6,*) '0+DEVI:',0+DEVI
      WRITE(6,*) 'D4*G:',DSTRESS(:NTENS)
      END IF
      
C     inter-seismic information
      IF ((OPEN_EAR.LT.0.5).AND.(DTIME.LT.TIME_EAR+0.1)
     1                            .AND.(JSTEP(1).GT.2)) then
      call MutexLock(4)
      OPEN_EAR=ONE
      call MutexUnlock(4)
      END IF
      
C     Information marking the time when the earthquake occurred
      IF ((C_DRO.GT.1E-10).AND.(JSTEP(1).GT.2).AND.(KINC.GE.2)) THEN
      call MutexLock(5)
      TEMP_EAR(NOEL,NPT)=1
      call MutexUnlock(5)
c      DDSDDE(:NTENS,:NTENS)=D4(:NTENS,:NTENS)
      STATEV(1)=ONE
      PNEWDT=TIME_EAR/DTIME
      WRITE(6,*) 'NOEL,NPT,N_Newton,FF1,FF3',NOEL,NPT,N_Newton,FF1,FF3
      GOTO 121
      END IF
      
112   CONTINUE
      STATEV(1)=ZRO

121   CONTINUE    
      
      STATEV(2)=SQRT(J2)
      STATEV(3)=-(I1+I1_EXC)/3
c      STATEV(4)=FF1
      
      TEK(1)=KINC
      
      IF (NOEL.EQ.439) THEN
      WRITE(6,*) 'ALPHA*I1,SQRT(J2),BETA,HMO*K2:',ALPHA*(I1+I1_EXC),
     1            SQRT(J2),BETA,HMO*K2
      WRITE(6,*) 'FF1:',FF1
      WRITE(6,*) 'STRESS2:',STRESS
      WRITE(6,*) 'DSTRESS(:NTENS):',DSTRESS(:NTENS)
      WRITE(6,*) 'C,C_DRO,SUM(DSTRAN:',C,C_DRO,ABS(SUM(DSTRAN2(:NTENS)))
      WRITE(6,*) 'PNEWDT:',PNEWDT
      WRITE(6,*) 'DP4,DH,DP4-DH:',DP4,DH,DP4-DH
      WRITE(6,*) 'STATEV(1:NSTATV):',STATEV(1:NSTATV)
      WRITE(6,*) 'TEMP_EAR(NOEL,NPT):',TEMP_EAR(NOEL,NPT)
      WRITE(6,*) 'DSGM:',DSGM
      WRITE(6,*) 'DP:',DP
      WRITE(6,*) 'DL:',DL
      WRITE(6,*) 'DDSDDE:',DDSDDE
      END IF
      
      RETURN
      END
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C     Stress invariants  
      SUBROUTINE STRESS_IJ(STRESS,NTENS,I1,J2)
      INTEGER NTENS
      DOUBLE PRECISION, INTENT(IN):: STRESS(NTENS)
      DOUBLE PRECISION, INTENT(OUT):: I1,J2
      
      I1=SUM(STRESS(1:3))    ! Positive compression and negative tension in rock and soil
      J2=((STRESS(1)-STRESS(2))**2+(STRESS(2)-STRESS(3))**2+
     1    (STRESS(3)-STRESS(1))**2+
     2    6*(STRESS(4)**2+STRESS(5)**2+STRESS(6)**2))/6.0
      
      RETURN
      END
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C     Elastic Matrix 
      SUBROUTINE ELAMATRIX(D,E,BSB,NTENS,NDI)
      INTEGER NTENS,I,J,NDI
      DOUBLE PRECISION, INTENT(IN):: E,BSB
      DOUBLE PRECISION, INTENT(OUT):: D(NTENS,NTENS)
      
      D=0
      DO I=1,NTENS-NDI
      DO J=1,NTENS-NDI
      D(I,J)=BSB
      END DO
      D(I,I)=1-BSB
      D(I+NDI,I+NDI)=0.5-BSB
      END DO
      D=D*(E/((1.0+BSB)*(1.0-2.0*BSB)))
      
      RETURN
      END
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C     Viscosity Matrix
      SUBROUTINE VISMATRIX(Q,ND,NTENS,NDI)
      INTEGER NTENS,I,NDI
      DOUBLE PRECISION, INTENT(IN):: ND
      DOUBLE PRECISION, INTENT(OUT):: Q(NTENS,NTENS)
      
      Q=0
      DO I=1,NTENS-NDI
      DO J=1,NTENS-NDI
      Q(I,J)=-1.0/6.0
      END DO
      Q(I,I)=1.0/3
      Q(I+3,I+3)=1.0
      END DO
      Q=Q/ND
      
      RETURN
      END 
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C     The partial derivative matrix of FF
      SUBROUTINE FMATRIX(AMAT,ALPHA,FF,BMAT,NTENS)
      INTEGER NTENS
      DOUBLE PRECISION, INTENT(IN):: AMAT(NTENS),ALPHA,FF
      DOUBLE PRECISION, INTENT(OUT):: BMAT(NTENS)
      
      BMAT(1)=(2*AMAT(1)-AMAT(2)-AMAT(3))/(FF*6.0)+ALPHA*1.0
      BMAT(2)=(2*AMAT(2)-AMAT(3)-AMAT(1))/(FF*6.0)+ALPHA*1.0
      BMAT(3)=(2*AMAT(3)-AMAT(1)-AMAT(2))/(FF*6.0)+ALPHA*1.0
      BMAT(4)=AMAT(4)/(FF*1.0)
      BMAT(5)=AMAT(5)/(FF*1.0)
      BMAT(6)=AMAT(6)/(FF*1.0)
      
      RETURN
      END
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C     The partial derivative matrix of G
      SUBROUTINE GMATRIX(AMAT,FF,BMAT,NTENS)
      INTEGER NTENS
      DOUBLE PRECISION, INTENT(IN)::  AMAT(NTENS),FF
      DOUBLE PRECISION, INTENT(OUT)::  BMAT(NTENS)
      
      BMAT(1)=(2*AMAT(1)-AMAT(2)-AMAT(3))/(FF*6.0)
      BMAT(2)=(2*AMAT(2)-AMAT(3)-AMAT(1))/(FF*6.0)
      BMAT(3)=(2*AMAT(3)-AMAT(1)-AMAT(2))/(FF*6.0)
      BMAT(4)=AMAT(4)/(FF*1.0)
      BMAT(5)=AMAT(5)/(FF*1.0)
      BMAT(6)=AMAT(6)/(FF*1.0)
      
      RETURN
      END 
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C     Solving for the inverse of a matrix
C     Gauss-Jordan Elimination Method to Compute Inverse Matrix
      SUBROUTINE INVERT(A2,AINV,N)
      INTEGER N,I,J,K
      DOUBLE PRECISION, INTENT(IN):: A2(N,N)
      DOUBLE PRECISION, INTENT(OUT):: AINV(N,N)
      DOUBLE PRECISION AMAT(N,N),FACTOR
      
C     Initialize the inverse matrix
      AINV=0.
      FACTOR=0
      AMAT=A2
      
      DO I=1,N
      AINV(I,I)=1.
      END DO
      
      DO K=1,N
      FACTOR=AMAT(K,K)
      
      DO J=1,N
      AMAT(K,J)=AMAT(K,J)/FACTOR
      AINV(K,J)=AINV(K,J)/FACTOR
      END DO
      
      DO I=1,N
      IF (I.NE.K) THEN
      FACTOR = AMAT(I,K)
      DO J=1,N
      AMAT(I,J) = AMAT(I,J) - FACTOR * AMAT(K,J)
      AINV(I,J) = AINV(I,J) - FACTOR * AINV(K,J)
      END DO
      END IF
      END DO
      END DO
      
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C     Square matrix addition
      SUBROUTINE MATRIX_ADD(AMAT,BMAT,C,N)
      INTEGER N,I,J
      DOUBLE PRECISION, INTENT(IN):: AMAT(N,N),BMAT(N,N)
      DOUBLE PRECISION, INTENT(OUT):: C(N,N)
      
C     Initialize the matrix
      C=0.
      DO I=1,N
      DO J=1,N
      C(I,J)=AMAT(I,J)+BMAT(I,J)
      ENDDO
      ENDDO

      RETURN
      END    

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C     Addition of rows or columns of matrices
      SUBROUTINE MATRIX_ADD2(AMAT,BMAT,C,N)
      INTEGER N,I
      DOUBLE PRECISION, INTENT(IN):: AMAT(N),BMAT(N)
      DOUBLE PRECISION, INTENT(OUT):: C(N)
      
C     Initialize the matrix
      C=0.
      DO I=1,N
      C(I)=AMAT(I)+BMAT(I)
      ENDDO

      RETURN
      END    
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C     Square matrix multiplication
      SUBROUTINE MATRIX_MUL(AMAT,BMAT,C,N)
      INTEGER N,I,J,K
      DOUBLE PRECISION, INTENT(IN):: AMAT(N,N),BMAT(N,N)
      DOUBLE PRECISION, INTENT(OUT):: C(N,N)
      
C     Initialize the matrix
      C=0.
      DO I=1,N
      DO J=1,N
      DO K=1,N
      C(I,J)=C(I,J)+AMAT(I,K)*BMAT(K,J)
      ENDDO
      ENDDO
      ENDDO

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C     Multiplying a square matrix with a row-column matrix 
C     C(N)=A(N,M)*B(M)
      SUBROUTINE MATRIX_MUL2(AMAT,BMAT,C,N)
      INTEGER N,I,J
      DOUBLE PRECISION, INTENT(IN):: AMAT(N,N),BMAT(N)
      DOUBLE PRECISION, INTENT(OUT):: C(N)
      
C     Initialize the matrix
      C=0.
      DO I=1,N
      DO J=1,N
      C(I)=C(I)+AMAT(I,J)*BMAT(J)
      ENDDO
      ENDDO

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C     Multiplying a row-column matrix with a square matrix
C     C(N)=B(N)*A(N,M)
      SUBROUTINE MATRIX_MUL22(BMAT,AMAT,C,N)
      INTEGER N,I,J
      DOUBLE PRECISION, INTENT(IN):: BMAT(N),AMAT(N,N)
      DOUBLE PRECISION, INTENT(OUT):: C(N)
      
C     Initialize the matrix
      C=0.
      DO I=1,N
      DO J=1,N
      C(I)=C(I)+BMAT(J)*AMAT(J,I)
      ENDDO
      ENDDO

      RETURN
      END
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C     Multiplying row and column matrices to get a square matrix
C     C(N,N)=A(N)*B(N)
      SUBROUTINE MATRIX_MUL3(AMAT,BMAT,C,N)
      INTEGER N,I,J
      DOUBLE PRECISION, INTENT(IN):: AMAT(N),BMAT(N)
      DOUBLE PRECISION, INTENT(OUT):: C(N,N)
      
      C=0.
      DO I=1,N
      DO J=1,N
      C(I,J)=AMAT(I)*BMAT(J)
      END DO
      END DO
      
      RETURN
      END
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C     Multiplying rows and columns together to get a single number       
C     C=A(N)*B(N)
      SUBROUTINE MATRIX_MUL4(AMAT,BMAT,C,N)
      INTEGER N,I
      DOUBLE PRECISION, INTENT(IN):: AMAT(N),BMAT(N)
      DOUBLE PRECISION, INTENT(OUT):: C
      
      C=0
      DO I=1,N
      C=C+AMAT(I)*BMAT(I)
      END DO
      
      RETURN
      END
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE  DISP(U,KSTEP,KINC,TIME,NODE,NOEL,JDOF,COORDS)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION U(3),TIME(2),COORDS(3)
      double precision :: point_x(10), point_z(10)
      logical :: inside
      INTEGER NUM_POINT
      
c     Boundary node information of one side loading range (in order)
      NUM_POINT=4
      point_x(1)=0
      point_z(1)=0
      point_x(2)=50E+03
      point_z(2)=0
      point_x(3)=50E+03
      point_z(3)=-50E+03
      point_x(4)=0
      point_z(4)=-50E+03
      
      IF (JDOF.EQ.2) then
c     The number of border nodes needs to be modified
c      call polygon(NUM_POINT,point_x(:NUM_POINT),point_z(:NUM_POINT),
c     1            COORDS(1),COORDS(3),inside)
c      if (inside) then
c      U(1)=(9.51e-11)*TIME(1)
c      U(2)=9.51e-11
c      else
c      U(1)=(-9.51e-11)*TIME(1)
c      U(2)=-9.51e-11
c      end if
      IF (COORDS(1).LT.50E+3) THEN
      U(1)=(9.51e-11)*TIME(1)
      U(2)=9.51e-11
      ELSE IF (COORDS(1).GT.50E+3) THEN
      U(1)=(-9.51e-11)*TIME(1)
      U(2)=-9.51e-11
      END IF
      END IF
      
      RETURN
      END
      
      
      
C############################################################################## 
      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME, 
     1                 NUVARM,NOEL,NPT,NLAYER,NSPT,KSTEP,KINC, 
     2                 NDI,NSHR,COORD,JMAC,JMATYP,MATLAYO, LACCFLG) 
CC 
      INCLUDE 'ABA_PARAM.INC' 
CC 
      CHARACTER*80 CMNAME,ORNAME 
      DIMENSION UVAR(*),TIME(2),DIRECT(3,3),T(3,3),COORD(*), 
     1     JMAC(*),JMATYP(*)  
CC     USER DEFINED DIMENSION STATEMENTS 
      CHARACTER*3 FLGRAY(15) 
      DIMENSION ARRAY(15),JARRAY(15) 
C
      DOUBLE PRECISION DSIG(6),DSIGPS(3),J2,I1
	  
C      WRITE(6,*)'IP:',IP
C---------非断层单元          
      CALL GETVRM('S',ARRAY,JARRAY,FLGRAY,JRCD,
     1     JMAC,JMATYP,MATLAYO, LACCFLG)
      DSIG(1)=ARRAY(1)
      DSIG(2)=ARRAY(2)
      DSIG(3)=ARRAY(3)
      DSIG(4)=ARRAY(4)
      DSIG(5)=ARRAY(5)
      DSIG(6)=ARRAY(6)

C---------主应力      
      CALL GETVRM('SP',ARRAY,JARRAY,FLGRAY,JRCD,
     1     JMAC,JMATYP,MATLAYO, LACCFLG)
      DSIGPS(1)=ARRAY(3)
      DSIGPS(2)=ARRAY(2)
      DSIGPS(3)=ARRAY(1)
      
      I1=DSIG(1)+DSIG(2)+DSIG(3)
      J2=1.0/6.0*((DSIG(1)-DSIG(2))**2+(DSIG(2)-DSIG(3))**2+
     1    (DSIG(3)-DSIG(1))**2+6*(DSIG(4)**2+DSIG(5)**2+DSIG(6)**2))
      UVAR(1)=SQRT(J2)
      UVAR(2)=I1*0.0325+SQRT(J2)
      
      return 
      END

C###################################################################
C####################################################################
      subroutine polygon(pn,point_xx,point_yy,pxx,pyy,inside)
      implicit none
      integer, intent(in) :: pn
      double precision, intent(in) :: pxx,pyy
      double precision, intent(in), dimension(*) :: point_xx,point_yy
      logical, intent(out) :: inside
      integer :: i, j
      logical :: on_edge
      inside = .false.
      on_edge = .false.
      
      j = pn
      do i = 1, pn
      ! 判断是否在边上
      call segment(point_xx(i),point_yy(i),point_xx(j),
     1                            point_yy(j),pxx,pyy,on_edge)
      if (on_edge) then
      inside = .true.
      exit
      end if
      
      ! 射线法判断是否在内部
      if (((point_yy(i)>pyy).neqv.(point_yy(j)>pyy)).and.
     1    (pxx<(point_xx(j)-point_xx(i))*(pyy-point_yy(i))/
     2                (point_yy(j)-point_yy(i))+point_xx(i))) then
      inside = .not. inside
      end if
      j = i
      end do
      
      return
      end subroutine polygon
      
      
C###################################################################
C####################################################################      
      subroutine segment(x1, y1, x2, y2, pxx, pyy,on_segment)
      implicit none
      double precision, intent(in) :: x1, y1, x2, y2, pxx, pyy
      logical :: on_segment
      double precision :: cross, dot, lendot
      
      cross = (pxx - x1) * (y2 - y1) - (pyy - y1) * (x2 - x1)
      if (abs(cross) > 1.0D-6) then
         on_segment = .false.
         goto 10
      end if
      
      dot = (pxx - x1) * (x2 - x1) + (pyy - y1) * (y2 - y1)
      if (dot < 0.0D0) then
         on_segment = .false.
         goto 10
      end if
      
      lendot = (x2 - x1) ** 2 + (y2 - y1) ** 2
      if (dot > lendot) then
         on_segment = .false.
         goto 10
      end if
      
      on_segment = .true.
10    return
      end subroutine segment

C###################################################################
C####################################################################
C###################################################################
      SUBROUTINE UMESHMOTION(UREF,ULOCAL,NODE,NNDOF,
     *    LNODETYPE,ALOCAL,NDIM,TIME,DTIME,PNEWDT,
     *    KSTEP,KINC,KMESHSWEEP,JMATYP,JGVBLOCK,LSMOOTH)
      
      USE TEMPARRAY_1
      INCLUDE 'ABA_PARAM.INC'

      DIMENSION ULOCAL(NDIM)
      DIMENSION ALOCAL(NDIM,*),TIME(2)
      DIMENSION JMATYP(*),JGVBLOCK(*)
      DIMENSION ARRAY(15)
      DIMENSION JELEMLIST(12),JELEMTYPE(12)
      
      INTEGER ISTEP,NSTEP,NODE_OR,I,LEV1,LEV2,LEV3,LEV4,
     1        K1,K2,N_MESHSWEEP,POSI,BONC,BONC2,I1,I2,
     2        J,NCOR,DEP_IN
      DOUBLE PRECISION KFSED,M,N,KDSED,GGB,D_TIME,H_LIMIT,
     1                 INCE_H,DIST_N1,RAP_ERO,SED_Z,RAP_KRO,
     2                 Z_LOW,COR_DEP,RAND_Z,Ocean_h
      DOUBLE PRECISION U_ARRAY(3),X_COR(3),WVGLOBAL(3),WVLOCAL(3),
     1                 ALOCAL2(3,3),ALOCAL3(3,3),ALOCAL4(3,3),UL2(3)
      DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE:: HIGHT,KFF,KDD,UU,
     1    CHIT,HHT2,ACDA,HIGHT2,PRE,ERO,ERO_RA,CATC,H3,U_BASEMENT

      PARAMETER (MONTH=2592000,YEAR=31500000)
      
      IF (NODE.EQ.855) THEN
      WRITE(6,*) '****************************************'
      WRITE(6,*) 'UMESHMOTION_KSTEP:',KSTEP
      WRITE(6,*) 'KINC:',KINC
      WRITE(6,*) 'KMESHSWEEP:', KMESHSWEEP
      WRITE(6,*) 'NODE,LNODETYPE:', NODE,LNODETYPE
      WRITE(6,*) 'TIME(2),DTIME:',TIME(2),DTIME
      WRITE(6,*) 'TEMP_NODE,TEMP_KINC',TEMP_NODE(NODE),TEMP_KINC(NODE)
      WRITE(6,*) 'TEMP_TIME(NODE)',TEMP_TIME(NODE)
      WRITE(6,*) 'MINI_NODE,MAX_NODE,NUM_NODE,LIE)',MINI_NODE,MAX_NODE,
     1                                                    NUM_NODE,LIE
      END IF
      
      IF (KINC.LE.1) TEMP_KINC(NODE)=1
      IF (KINC.LE.1) TEMP_TIME(NODE)=0
      IF (KINC.GT.TEMP_KINC(NODE)) TEMP_NODE(NODE)=0
      
C     Read information
      IF (KINC.EQ.1) THEN
      CALL SetPath()
      call MutexLock(2)
      MODEL_COOR=0
      ELE=0
      H_HISTORY=0
      COR=0
      CALL READ_FILE3(PATH1,MINI_NODE,MAX_NODE,NUM_NODE,LIE) ! MINI_NODE/MAX_NODE/NUM_NODE/ELE_NM
      CALL READ_FILE1(PATH,NUM_NODE,5,MODEL_COOR(:NUM_NODE,:5)) ! MODEL_COOR(NN,5)=NODE/X/Y/Z/precipitation
      CALL READ_FILE2(PATH2,ELE(:NUM_NODE,:LIE),NUM_NODE,LIE)  ! ELE(NN,<=14)=NODE/AREAS/BOUNDARY/NODE_NUM/NODE1/NODE2/NODE3....(Left, right, top, and bottom node join order)
                                                                      !   BOUNDARY-(Corner points belong to the left and right boundary numbers)
      CALL CORNER_POINT(MODEL_COOR(:NUM_NODE,3),
     1    ELE(:NUM_NODE,3),ELE(:NUM_NODE,4),COR(:4),NUM_NODE)
      call MutexUnlock(2)
      END IF
      
C     Initial definition of array
      ALLOCATE(HIGHT(NUM_NODE),KFF(NUM_NODE),KDD(NUM_NODE),UU(NUM_NODE),
     1        CHIT(NUM_NODE),ACDA(NUM_NODE),HIGHT2(NUM_NODE),
     3        ERO(NUM_NODE),ERO_RA(NUM_NODE),CATC(NUM_NODE),
     4        HHT2(NUM_NODE),H3(NUM_NODE),PRE(NUM_NODE),
     5        U_BASEMENT(NUM_NODE))
      
      HIGHT=0
      KFF=0
      KDD=0
      UU=0
      CHIT=0
      H3=0
      ACDA=0
      HIGHT2=0
      ERO=0
      ERO_RA=0
      CATC=0
      U_BASEMENT=0
      HHT2=0
      PRE=0
      
      DO I1=1,NUM_NODE
      IF (INT(MODEL_COOR(I1,1)).EQ.NODE) THEN
      POSI=I1  ! Determine the location of the point
      EXIT
      END IF
      END DO
      
C      N_MESHSWEEP=3       ! Maximum number of grid scans +1
      OPEN_KINC=3         ! Kinc to start the surface process(advices >=2 )
      
      IF (DTIME.LT.MONTH) GOTO 97 ! Skip Earthquake step
      IF (KSTEP.LE.2) GOTO 97 ! Step 1 is Geostress Balance; Step 2 is to detect the stress balance.
      IF ((KINC.LT.OPEN_KINC).AND.(KSTEP.LE.3)) GOTO 97  ! The initial terrain is generated by the previous kinc
      IF (TEMP_NODE(NODE).GE.1) GOTO 96
      IF (KMESHSWEEP.GT.0) GOTO 96
      
      
      LSMOOTH=0       ! Surface smoothing (Set LSMOOTH to 1)
      LTRN=0          ! Global coordinate system=0, Local transformed system=1
      
      KFF = 1E-6     ! Bedrock river incision rate(nx*ny)--m*year^(-1)
      KFSED = 1E-6   ! Sediment river incision rate(nx*ny)--m*year^(-1)
      M = 0.4        ! Drainage area exponent
      N = 1.0        ! Slope exponent
      KDD = 1E-2     ! Bedrock diffusivity coefficient for hillslope processes(nx*ny)--m^(2)/year
      KDSED = 1E-2   ! Sediment diffusivity coefficient for hillslope processes(nx*ny)--m^(2)/year
      GGB = 5        ! Bedrock dimensionless deposition coefficient(must >= 0)
      
      UU = 0          ! uplift velocity(set to 0 in Abaqus)(nx*ny)--m/year
      NSTEP = 1       ! Built-in time step number (set to 1 in Abaqus)
      PRE=2         ! precipitation rate(nx*ny)--m/year
cccc      PRE(:NUM_NODE)=MODEL_COOR(:NUM_NODE,5)
      
      BONC=0100       ! Boundary conditions of the model (1: fixed boundary - ocean; 0: variable boundary)(Order: Down/Right/Up/Left)
      BONC2=0001      ! For variable boundary, when the option 0 is selected in bonc, you can also choose to construct only -1 or construct plus erosion -0.
      
      RAP_ERO=1      ! KFSED amplification factor (Effective within 10 years after the earthquake) 50
      RAP_KRO=1     ! KDSED amplification factor (Effective within 10 years after the earthquake) 100
      Z_LOW=0         ! The lowest z coordinate of the model--m
      SED_Z=5000      ! Thickness of Weathered Cover--m
      Ocean_h=99990       ! Ocean constant height (when BONC=1)--m
      COR_DEP=0       ! Initial terrain drop at the border--m
      RAND_Z=1        ! Initial random terrain elevation--m
      
      D_TIME=(DTIME*1.0)/(YEAR*1.0)
      IF (TEMP_TIME(NODE).EQ.1) KFSED = KFSED*RAP_ERO ! The rapid erosion stage after an earthquake
      IF (TEMP_TIME(NODE).EQ.1) KDSED = KDSED*RAP_KRO ! The rapid erosion stage after an earthquake
      
C     FastScape Program
      IF ((NODE.EQ.MINI_NODE).AND.(KMESHSWEEP.EQ.0)) THEN
C-----initialize
      call MutexLock(1)
      call FastScape_Init ()
      
C-----set grid size
      call FastScape_Set_NX_NY (NUM_NODE,LIE-2)
      
C-----allocate memory
      call FastScape_Setup ()
      
C-----set time step
      call FastScape_Set_DT (D_TIME)
      
C-----set random initial topography
      DO I=1,NUM_NODE
      CALL GETVRN(INT(MODEL_COOR(I,1)),'COORD',ARRAY,JRCD,JGVBLOCK,LTRN)
      MODEL_COOR(I,2:4)=ARRAY(1:3)    ! Update horizontal coordinates
      HIGHT(I)=ARRAY(3)
      END DO
      
      call random_number (HIGHT2)
      IF ((KINC.EQ.OPEN_KINC).AND.(KSTEP.EQ.3)) 
     1    HIGHT(:NUM_NODE)=HIGHT(:NUM_NODE)+HIGHT2(:NUM_NODE)*RAND_Z
      
C     Calculate the base elevation at this step
      IF ((KINC.EQ.OPEN_KINC).AND.(KSTEP.EQ.3)) 
     1    H_HISTORY(:NUM_NODE)=HIGHT(:NUM_NODE)-SED_Z
      IF ((KINC.GT.OPEN_KINC).OR.(KSTEP.GT.3)) THEN
      H_HISTORY(1:NUM_NODE)=H_HISTORY(:NUM_NODE)
     1                +HIGHT(:NUM_NODE)-HH(:NUM_NODE,1)
      H_HISTORY(1:NUM_NODE)=MIN(HIGHT(:NUM_NODE),H_HISTORY(:NUM_NODE))
      END IF
      BASEMENT(1:NUM_NODE)=H_HISTORY(1:NUM_NODE)
      U_BASEMENT(1:NUM_NODE)=(HIGHT(1:NUM_NODE)-HH(1:NUM_NODE,1))/D_TIME
      
      call FastScape_Init_H (HIGHT(1:NUM_NODE),BASEMENT(1:NUM_NODE))
C-----Set the node's rainfall amount
      CALL FastScape_Set_Precip (PRE(:NUM_NODE))
C-----set erosional parameters
      call FastScape_Set_Erosional_Parameters (KFF(1:NUM_NODE),KFSED,
     1        M,N,KDD(1:NUM_NODE),KDSED,GGB,GGB,-2.d0)

C-----set uplift rate (uniform while keeping boundaries at base level)
CCCC      call FastScape_Set_U (UU(1:NUM_NODE),ELE(1:NUM_NODE,3:LIE),
CCCC     1                        NUM_NODE,LIE-2)
      
C-----set boundary conditions
      call FastScape_Set_BC (BONC,ELE(1:NUM_NODE,3:LIE),  
     1                NUM_NODE,LIE-2,MODEL_COOR(1:NUM_NODE,2:3),2)
C-----set number of time steps and initialize counter istep (istep=0)
      
      call FastScape_Get_Step (ISTEP)
      
C-----set model dimensions
      call FastScape_Set_XL_YL (ELE(1:NUM_NODE,2))
C-----loop on time stepping
      do while (istep<nstep)
C     execute step
      call FastScape_Execute_Step()
C     get value of time step counter
      call FastScape_Get_Step (istep)
C     extract solution
      call FastScape_Copy_Chi (CHIT(1:NUM_NODE))
C     outputs values
      call FastScape_Copy_H(H3(1:NUM_NODE))
      call FastScape_Copy_HH2(HHT2(1:NUM_NODE))
      call FastScape_Copy_Drainage_Area(ACDA(1:NUM_NODE))
      call FastScape_Copy_Total_Erosion(ERO(1:NUM_NODE))
      call FastScape_Copy_Erosion_Rate(ERO_RA(1:NUM_NODE))
      call FastScape_Copy_Catchment(CATC(1:NUM_NODE))
      enddo

C     output timing
      call FastScape_Debug()
      
C     Stores procedure variables
      HH(1:NUM_NODE,1)=H3(1:NUM_NODE)
      HH(1:NUM_NODE,2)=HHT2(1:NUM_NODE)
      
C     Fixed elevation of ocean terrain edge (When BONC=1)
      IF ((KINC.LE.OPEN_KINC).AND.(KSTEP.EQ.3)) THEN
      CALL FIX_OCEAN(ELE(:NUM_NODE,3:4),MODEL_COOR(:NUM_NODE,3),
     1                HH(:NUM_NODE,1),BONC,NUM_NODE,Ocean_h)
      END IF
      
      CALL Uneroded_boundary(ELE(:NUM_NODE,3:4),MODEL_COOR(:NUM_NODE,3),
     1                HH(:NUM_NODE,1),HIGHT(:NUM_NODE),BONC2,NUM_NODE)
      
C     Initial stage: Terrain edge elevation setting
      IF ((KINC.EQ.OPEN_KINC).AND.(KSTEP.EQ.3)) THEN
      DO I=1,NUM_NODE
      IF ((ELE(I,3).GT.0).OR.(I.EQ.COR(1)).OR.(I.EQ.COR(2)).OR.
     1        (I.EQ.COR(3)).OR.(I.EQ.COR(4))) HH(I,1)=HH(I,1)+COR_DEP
c      IF (MODEL_COOR(I,2).LE.150E+3) HH(I,1)=HH(I,1)+1000
      END DO
      END IF
      
C     Print step information
      CALL Print_ste(KINC,KSTEP,TIME(2),HH(:NUM_NODE,1),ACDA(:NUM_NODE),
     1                ERO(:NUM_NODE),ERO_RA(:NUM_NODE),
     2                CATC(:NUM_NODE),H_HISTORY(:NUM_NODE),
     3                U_BASEMENT(:NUM_NODE))
      
C     end FastScape run
      WRITE(6,*) 'TEMP_TIME,KFSED,KDSED', TEMP_TIME(NODE),KFSED,KDSED
      WRITE(6,*) 'success_surface process', NODE
      WRITE(6,*) 'COR(I)',(COR(I),I=1,4)
      WRITE(6,*) 'HH4',(HH(COR(I),1),I=1,4)
      WRITE(6,*) 'DEP_IN,HH(DEP_IN,1)',DEP_IN,HH(DEP_IN,1)
      
      call FastScape_Destroy ()
      call MutexUnlock(1)
      END IF
      
C     Update elevation for first scan
      IF (NODE.EQ.855) THEN
      WRITE(6,*) '1111'
      WRITE(6,*) 'ALOCAL',(ALOCAL(I,:3),I=1,3)
      WRITE(6,*) 'ULOCAL1',ULOCAL(:3)
      END IF
      
      WVLOCAL=0
      CALL GETVRN(NODE,'U',U_ARRAY,JRCD,JGVBLOCK,LTRN)
      CALL GETVRN(NODE,'COORD',ARRAY,JRCD,JGVBLOCK,LTRN)
      INCE_H=HH(POSI,1)-ARRAY(3)
C      H_LIMIT=KFF(1)*(DTIME/(YEAR*1.0))
C      IF (ABS(INCE_H).GE.H_LIMIT) INCE_H = SIGN(H_LIMIT,INCE_H) 
      IF (LNODETYPE.LE.5) THEN    ! Temporarily assume that all node types require a conversion matrix
      X_COR(1:2)=ARRAY(1:2)
      X_COR(3)=Z_LOW
      DIST_N1=SQRT((X_COR(1)-ARRAY(1))**2+(X_COR(2)-ARRAY(2))**2+
     1        (X_COR(3)-ARRAY(3))**2)
      WVGLOBAL(1)=INCE_H*(-X_COR(1)+ARRAY(1))/DIST_N1
      WVGLOBAL(2)=INCE_H*(-X_COR(2)+ARRAY(2))/DIST_N1
      WVGLOBAL(3)=INCE_H*(-X_COR(3)+ARRAY(3))/DIST_N1
      CALL MATRIX_MUL22(WVGLOBAL(:NDIM),ALOCAL(:NDIM,:NDIM),
     1                                        WVLOCAL(:NDIM),NDIM)
      ULOCAL(:NDIM)=ULOCAL(:NDIM)+WVLOCAL(:NDIM)
      ELSE
      ULOCAL(NDIM)=ULOCAL(NDIM)+INCE_H
      END IF
      
      ULOCAL_TEMP(POSI,:3)=ULOCAL(:3)
      ALOCAL_TEMP(POSI,:3)=ALOCAL(1,:3)
      ALOCAL_TEMP(POSI,4:6)=ALOCAL(2,:3)
      ALOCAL_TEMP(POSI,7:9)=ALOCAL(3,:3)
      
      IF (NODE.EQ.855) THEN
      WRITE(6,*) 'POSI',POSI
      WRITE(6,*) 'U_ARRAY',U_ARRAY(:3)
      WRITE(6,*) 'COORD',ARRAY(:3)
      WRITE(6,*) 'ULOCAL2',ULOCAL(:3)
      WRITE(6,*) 'ULOCAL_TEMP',ULOCAL_TEMP(POSI,:3)
      WRITE(6,*) 'H_HIS,hjght,hh',H_HISTORY(POSI),ARRAY(3),HH(POSI,1)
      WRITE(6,*) 'INCE_H,HH2,H_LIMIT',INCE_H,HH(POSI,2),H_LIMIT
      WRITE(6,*) 'WVGLOBAL(:3)',WVGLOBAL(:3)
      WRITE(6,*) 'WVLOCAL(:3)',WVLOCAL(:3)
      END IF
      
      IF (NODE.EQ.855) THEN
      WRITE(6,*) 'finish1'
      END IF
      
      GOTO 97
      
96    IF (NODE.EQ.855) THEN
      WRITE(6,*) '2222'
      WRITE(6,*) 'ALOCAL',(ALOCAL(I,:3),I=1,3)
      WRITE(6,*) 'ULOCAL1',ULOCAL(:3)
      END IF
      
C      IF (TEMP_NODE(NODE).GE.N_MESHSWEEP) GOTO 499
      
C     Update elevation for subsequent scans   
      ALOCAL2(1,:3)=ALOCAL_TEMP(POSI,:3)
      ALOCAL2(2,:3)=ALOCAL_TEMP(POSI,4:6)
      ALOCAL2(3,:3)=ALOCAL_TEMP(POSI,7:9)
      CALL INVERT(ALOCAL2(:3,:3),ALOCAL3(:3,:3),3)
      CALL MATRIX_MUL(ALOCAL3(:3,:3),ALOCAL(:3,:3),ALOCAL4(:3,:3),3)
      CALL MATRIX_MUL22(ULOCAL_TEMP(POSI,:3),ALOCAL4(:3,:3),UL2(:3),3)
      ULOCAL(:3)=UL2(:3)
      
499   CONTINUE
      
      CALL GETVRN(NODE,'U',U_ARRAY,JRCD,JGVBLOCK,LTRN)
      CALL GETVRN(NODE,'COORD',ARRAY,JRCD,JGVBLOCK,LTRN)
      INCE_H=HH(POSI,1)-ARRAY(3)
      
      IF (NODE.EQ.855) THEN
      WRITE(6,*) 'U_ARRAY',U_ARRAY(:3)
      WRITE(6,*) 'COORD',ARRAY(:3)
      WRITE(6,*) 'ULOCAL_TEMP',ULOCAL_TEMP(POSI,:3)
      WRITE(6,*) 'ULOCAL2',ULOCAL(:3)
      WRITE(6,*) 'POSI',POSI
      WRITE(6,*) 'INCE_H,HH2',INCE_H,HH(POSI,2)
      END IF
      
97    CONTINUE
      
      TEMP_NODE(NODE)=TEMP_NODE(NODE)+1
      TEMP_KINC(NODE)=KINC
      IF (DTIME.LT.MONTH) TEMP_TIME(NODE)=1
      IF (DTIME.GT.MONTH) TEMP_TIME(NODE)=0
      
c      deallocate (H3,HHT2,UU,KFF,KDD,CHIT,HIGHT,
c     1            HIGHT2,ACDA,PRE,ERO,ERO_RA,CATC)
      if (allocated(H3)) deallocate(H3)
      if (allocated(HHT2)) deallocate(HHT2)
      if (allocated(UU)) deallocate(UU)
      if (allocated(KFF)) deallocate(KFF)
      if (allocated(KDD)) deallocate(KDD)
      if (allocated(CHIT)) deallocate(CHIT)
      if (allocated(HIGHT)) deallocate(HIGHT)
      if (allocated(HIGHT2)) deallocate(HIGHT2)
      if (allocated(ACDA)) deallocate(ACDA)
      if (allocated(PRE)) deallocate(PRE)
      if (allocated(ERO)) deallocate(ERO)
      if (allocated(ERO_RA)) deallocate(ERO_RA)
      if (allocated(CATC)) deallocate(CATC)
      if (allocated(U_BASEMENT)) deallocate(U_BASEMENT)
      
      IF (NODE.EQ.855) THEN
      WRITE(6,*) 'finish2'
      END IF
      
      RETURN
      END
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCC---FastScapeContext---CCCCCCCCCCCCCCCCCCCCCCCCCCCC
      module FastScapeContext
      
      ! Context module for FastScape api
      ! should not be accessed or changed
      ! see API for name of routines and externally accessible variables
      
      implicit none
      
      integer :: nn, nstack, nn1
      integer :: bounds_ibc
      integer :: bounds_i1, bounds_i2, bounds_j1, bounds_j2
      logical :: bounds_xcyclic, bounds_ycyclic
      logical, dimension(:), allocatable :: bounds_bc
      integer :: step
      integer :: nGSStreamPowerLaw, nGSMarine
      logical :: setup_has_been_run
      double precision :: tol_rel, tol_abs
      integer :: nGSStreamPowerLawMax
      double precision, target, dimension(:), allocatable :: h,u,vx,
     1    vy,a,erate,etot,catch,catch0,b,precip,kf,kd,length,hh2
      double precision, dimension(:,:), allocatable :: length2
      double precision, target, dimension(:), allocatable :: Sedflux,
     1    Fmix
      double precision, target, dimension(:), allocatable :: gob,sear
      double precision, target, dimension(:), allocatable :: p_mfd_exp
      double precision, dimension(:), pointer, contiguous :: h2, vx2,
     1    vy2, etot2, b2
      double precision :: dt,kfsed,m,n,kdsed,g1,g2,p
      double precision :: sealevel, poro1, poro2, zporo1, zporo2, ratio,
     1    layer, kdsea1, kdsea2
      integer, dimension(:), allocatable :: stack, ndon, rec
      integer, dimension(:,:), allocatable :: don
      integer, dimension(:,:), allocatable :: bond_num
      logical :: runSPL, runAdvect, runDiffusion, runUplift,
     1    runMarine
      real :: timeSPL, timeAdvect, timeDiffusion,timeUplift, timeMarine
      double precision, dimension(:,:), allocatable :: reflector
      double precision, dimension(:,:,:), allocatable :: fields
      integer nfield, nfreq, nreflector, nfreqref, ireflector
      double precision :: vexref
      logical :: SingleFlowDirection
      double precision, dimension(:), allocatable :: lake_depth, hwater
      integer, dimension(:), allocatable :: mnrec,mstack
      integer, dimension(:,:), allocatable :: mrec
      double precision, dimension(:,:), allocatable :: mwrec,mlrec
      contains
c      
C---------------------------------------------------------------      
      subroutine Init()
      
      step=0
      setup_has_been_run = .false.
      timeSPL = 0.
      timeAdvect = 0.
      timeDiffusion = 0.
      timeMarine = 0.
      timeUplift = 0.
      
      end subroutine Init
      
C---------------------------------------------------------------
      subroutine SetUp()
      implicit none
      
      if (nn.eq.0) THEN
      WRITE(6,*) 'FastScapeSetup - You need to set nn first'
      END IF
      
      call Destroy()
      
      allocate (h(nn),u(nn),vx(nn),vy(nn),stack(nn),ndon(nn),rec(nn),
     1        don(10,nn),catch0(nn),catch(nn),precip(nn),
     2        bond_num(nn,nn1),hh2(nn))
      allocate (gob(nn),sear(nn))
      allocate (bounds_bc(nn))
      allocate (p_mfd_exp(nn))
      allocate (length(nn),a(nn),erate(nn),etot(nn),b(nn),Sedflux(nn),
     1        Fmix(nn),kf(nn),kd(nn),length2(nn,2))
      allocate (lake_depth(nn),hwater(nn),mrec(10,nn),mnrec(nn),
     1        mwrec(10,nn),mlrec(10,nn),mstack(nn))
      
      h2(1:nn) => h
      b2(1:nn) => b
      vx2(1:nn) => vx
      vy2(1:nn) => vy
      etot2(1:nn) => etot
      
C      call SetBC (1111)
      h = 0.d0
      hh2 = 0.d0
      u = 0.d0
      vx = 0.d0
      vy = 0.d0
      etot = 0.d0
      b = 0
      precip = 2.d0
      p_mfd_exp(1:nn) = 1.d0
      catch0=0
      sealevel = 0.d0
      Fmix = 0.5d0
      lake_depth = 0.d0
      
      runSPL = .false.
      runAdvect = .false.
      runDiffusion = .false.
      runMarine = .false.
      runUplift = .false.
      
      nGSStreamPowerLaw = 0
      nGSMarine = 0
      
      setup_has_been_run = .true.
      
      tol_rel = 1.d-4
      tol_abs = 1.d-4
      nGSStreamPowerLawMax = 100
      return
      
      end subroutine SetUp
      
C---------------------------------------------------------------
      subroutine Destroy()
      
      if (allocated(h)) deallocate(h)
      if (allocated(hh2)) deallocate(hh2)
      if (allocated(u)) deallocate(u)
      if (allocated(vx)) deallocate(vx)
      if (allocated(vy)) deallocate(vy)
      if (allocated(stack)) deallocate(stack)
      if (allocated(ndon)) deallocate(ndon)
      if (allocated(rec)) deallocate(rec)
      if (allocated(don)) deallocate(don)
      if (allocated(catch0)) deallocate(catch0)
      if (allocated(catch)) deallocate(catch)
      if (allocated(length)) deallocate (length)
      if (allocated(length2)) deallocate (length2)
      if (allocated(a)) deallocate (a)
      if (allocated(b)) deallocate (b)
      if (allocated(sedflux)) deallocate (sedflux)
      if (allocated(Fmix)) deallocate (Fmix)
      if (allocated(erate)) deallocate(erate)
      if (allocated(etot)) deallocate(etot)
      if (allocated(precip)) deallocate(precip)
      if (allocated(kd)) deallocate(kd)
      if (allocated(kf)) deallocate(kf)
      if (allocated(reflector)) deallocate(reflector)
      if (allocated(fields)) deallocate(fields)
      if (allocated(lake_depth)) deallocate(lake_depth)
      if (allocated(hwater)) deallocate(hwater)
      if (allocated(mrec)) deallocate(mrec)
      if (allocated(mnrec)) deallocate(mnrec)
      if (allocated(mwrec)) deallocate(mwrec)
      if (allocated(mlrec)) deallocate(mlrec)
      if (allocated(mstack)) deallocate(mstack)
      if (allocated(gob)) deallocate(gob)
      if (allocated(p_mfd_exp)) deallocate(p_mfd_exp)
      if (allocated(bounds_bc)) deallocate(bounds_bc)
      if (allocated(bond_num)) deallocate(bond_num)
      if (allocated(sear)) deallocate(sear)
      return
      
      end subroutine Destroy
      
C---------------------------------------------------------------
      subroutine CopyH (hp)
      double precision, intent(out), dimension(*) :: hp
      
      if (.not.setup_has_been_run) then
      write(6,*) 'CopyH - You need to run SetUp first'
      end if
      hp(1:nn)=h(1:nn)
      return
      
      end subroutine CopyH

C---------------------------------------------------------------
      subroutine CopyHH2 (hp2)
      double precision, intent(out), dimension(*) :: hp2
      
      if (.not.setup_has_been_run) then
      write(6,*) 'CopyH - You need to run SetUp first'
      end if
      hp2(1:nn)=hh2(1:nn)
      return
      
      end subroutine CopyHH2
      
C---------------------------------------------------------------
      subroutine CopyBasement (bp)
      double precision, intent(out), dimension(*) :: bp
      
      if (.not.setup_has_been_run) then
      write(6,*) 'CopyB - You need to run SetUp first'
      end if
      bp(1:nn)=b(1:nn)
      return
      
      end subroutine CopyBasement
      
C---------------------------------------------------------------
      subroutine CopyEtot (etotp)
      double precision, intent(inout), dimension(*) :: etotp
      
      if (.not.setup_has_been_run) 
     1       write(6,*) 'CopyEtot - You need to run SetUp first'
      etotp(1:nn)=etot(1:nn)
      return
      
      end subroutine CopyEtot
      
C---------------------------------------------------------------
      subroutine CopyArea (ap)
      double precision, intent(inout), dimension(*) :: ap
      
      if (.not.setup_has_been_run) then
      write(6,*) 'CopyArea - You need to run SetUp first'
      end if
      ap(1:nn)=a(1:nn)
      return
      
      end subroutine CopyArea
      
C---------------------------------------------------------------
      subroutine CopyErate (eratep)
      double precision, intent(inout), dimension(*) :: eratep
      
      if (.not.setup_has_been_run)
     1        write(6,*) 'CopyErate - You need to run SetUp first'
      eratep(1:nn)=erate(1:nn)
      return
      
      end subroutine CopyErate
      
C---------------------------------------------------------------
      subroutine Copychi (chip)
      double precision, intent(inout), dimension(*) :: chip
      double precision, dimension(:), allocatable :: chi
      integer ij,ijk
      double precision dx,dy,a0
      
      if (.not.setup_has_been_run) then
      write(6,*) 'CopyChi - You need to run SetUp first'
      end if
      
      allocate (chi(nn))
      chi=0.d0
      
      do ij=1,nn
        a0=sear(ij)*10.d0
        ijk=stack(ij)
        if (a(ijk).gt.a0) chi(ijk)=chi(rec(ijk))+
     1            (a0/a(ijk))**(m/n)*length(ijk)
      enddo
      chip(1:nn)=chi(1:nn)
      if (allocated(chi)) deallocate(chi)
      return
      
      end subroutine CopyChi
      
C---------------------------------------------------------------
      subroutine CopyCatchment (catchp)
      double precision, intent(inout), dimension(*) :: catchp
      
      if (.not.setup_has_been_run)
     1    write(6,*) 'CopyCatchment - You need to run SetUp first'
      catchp(1:nn)=catch(1:nn)
      return
      
      end subroutine CopyCatchment
      
cC---------------------------------------------------------------
c      subroutine CopyF (Fmixp)
c      double precision, intent(out), dimension(*) :: Fmixp
c      
c      if (.not.setup_has_been_run)
c     1    stop 'CopyF - You need to run SetUp first'
c      Fmixp(1:nn) = Fmix
c      return
c      
c      end subroutine CopyF
c      
cC---------------------------------------------------------------
c      subroutine CopyLakeDepth (Lp)
c      double precision, intent(out), dimension(*) :: Lp
c      
c      if (.not.setup_has_been_run)
c     1        stop 'CopyLakeDepth - You need to run SetUp first'
c      
c      Lp(1:nn) = lake_depth
c      return
c      
c      end subroutine CopyLakeDepth
c      
C---------------------------------------------------------------
      subroutine InitH (hp,hpp2)
      double precision, intent(in), dimension(*) :: hp,hpp2
      
      IF (.not.setup_has_been_run) THEN
      WRITE(6,*) 'InitH - You need to run SetUp first'
      END IF
      
      h(1:nn) = hp(1:nn)
      b(1:nn) = hpp2(1:nn)
      continue
      return
      
      end subroutine InitH
      
C---------------------------------------------------------------
c      subroutine InitF (Fmixp)
c      double precision, intent(in), dimension(*) :: Fmixp
c      
c      if (.not.setup_has_been_run)
c     1    stop 'InitF - You need to run SetUp first'
c      Fmix = Fmixp(1:nn)
c      return
c      
c      end subroutine InitF
c      
cC---------------------------------------------------------------
c      subroutine ResetCumulativeErosion ()
c      
c      etot = 0.d0
c      return
c      
c      end subroutine ResetCumulativeErosion
c      
cC---------------------------------------------------------------
c      subroutine View()
c      
c      write (*,*) 'FastScapeContext:'
c      write (*,*) 'nx,ny',nx,ny
c      write (*,*) 'nn',nn
c      write (*,*) 'step',step
c      write (*,*) 'xl,yl',xl,yl
c      write (*,*) 'dt',dt
c      write (*,*) 'Kf,Kfsed,,m,n,Kd,Kdsed,G1,G2',sum(kf)/nn,kfsed,m,n,
c     1            sum(kd)/nn,kdsed,g1,g2
c      write (*,*) 'ibc',bounds_ibc
c      write (*,*) 'h',minval(h),sum(h)/nn,maxval(h)
c      write (*,*) 'u',minval(u),sum(u)/nn,maxval(u)
c      return
c      
c      end subroutine View
c      
C---------------------------------------------------------------
      subroutine SetNXNY (nn_1,nn_2)
      integer, intent(in) :: nn_1,nn_2
      
      nn = nn_1
      nn1 = nn_2
      return
      
      end subroutine SetNXNY
     
C---------------------------------------------------------------
      subroutine SetXLYL (xlyl)
      integer, intent(in), dimension(*) :: xlyl

      sear(1:nn)=xlyl(1:nn)*precip(1:nn)   ! 采用Voronoi面积计算降雨量，precip是降雨速率m/yr
      return
      
      end subroutine SetXLYL
      
C---------------------------------------------------------------
      subroutine SetErosionalParam (kkf,kkfsed,mm,nnn,kkd,kkdsed,
     1                                gg1,gg2,pp)
      double precision, intent(in), dimension(*) :: kkf,kkd
      double precision, intent(in) :: kkfsed,mm,nnn,kkdsed,gg1,gg2,pp
      
      runSPL = .true.
      kf(1:nn) = kkf(1:nn)	! 基岩的河流切割参数（m^(1-2m)*yr^(-1)）
      kfsed = kkfsed			! 沉积物的河流切割参数（m^(1-2m)*yr^(-1)）；kfsed<0则全部采用基岩参数
      m = mm
      n = nnn
      kd(1:nn) = kkd(1:nn)	! 基岩的扩散系数在山坡过程中m2/yr
      kdsed = kkdsed			! 沉积物的扩散系数m2/yr （kdsed < 0时全部使用基岩扩散系数）
      g1 = gg1	! 基岩的沉积/运输系数（无量纲系数），g1是SPL模型中的G项；当G>0时,执行Gauss-Seidel迭代；G会影响迭代次数（建议值为1）
      g2 = gg2	! 沉积物的沉积/运输系数；当g2<0时，设置为g1
      p = pp		! 多向流分配系数
      p_mfd_exp(1:nn) = pp
      SingleFlowDirection = .false.
      if (pp.lt.-1.5d0) then
        SingleFlowDirection = .true.
        p = 1.d0
      endif
      
c      WRITE(6,*) 'ErosionalParam',kf(1),kfsed,m,n,kd(1),kdsed,g1,g2,p
      
      if (maxval(kd).gt.tiny(kd).or.kdsed.gt.tiny(kdsed)) 
     1        runDiffusion = .true.
      return
      
      end subroutine SetErosionalParam
      
C---------------------------------------------------------------
CCc      subroutine SetMarineParam (sl, p1, p2, z1, z2, r, l, kds1, kds2)
CCc      double precision, intent(in) :: sl,p1,p2,z1,z2,r,l,kds1,kds2
CCc      
CCc      runMarine = .true.
CCc      sealevel = sl
CCc      poro1 = p1
CCc      poro2 = p2
CCc      zporo1 = z1
CCc      zporo2 = z2
CCc      ratio = r
CCc      layer = l
CCc      kdsea1 = kds1
CCc      kdsea2 = kds2
CCc      return
CCc      
CCc      end subroutine SetMarineParam
CCc      
C---------------------------------------------------------------
      subroutine SetDT (dtt)
      double precision, intent(in) :: dtt
      
      dt = dtt
      return
      
      end subroutine SetDT
      
C---------------------------------------------------------------
CCc      subroutine GetSizes (nnx,nny)
CCc      integer, intent(out) :: nnx,nny
CCc      
CCc      nnx = nx
CCc      nny = ny
CCc      return
CCc      
CCc      end subroutine GetSizes
CCc      
C---------------------------------------------------------------
      subroutine GetStep (nstep)
      integer, intent(out) :: nstep
      
      nstep = step
      return
      
      end subroutine GetStep
      
C---------------------------------------------------------------
      subroutine Debug ()
      implicit none
      integer i,j,ij,counter
      
      write (6,*) '----------------------------------------------------'
      write (6,*) 'Time step', step
      write (6,*) 'Debug information'
      write (6,*) 'Total number of nodes (nx*ny)',nn
      write (6,*) 'Stack size',nstack
      write (6,*) 'Number of nil elements in stack',count(stack==0)
      write (6,*) 'Total number of donors',sum(ndon)
      
      counter=0
      do ij=1,nn
        if (rec(ij)==ij) counter=counter+1
      enddo
      
      write (6,*) 'Total number of self donors',counter
      
      counter=0
      do ij=1,nn
          if (bounds_bc(ij)) cycle
          if (rec(ij)==ij) counter=counter+1
      enddo

      
      write (6,*) 'Total number of local minima',counter
      write (6,*) 'Number of Gauss-Siedel iterations (SPL)',
     1                                        nGSStreamPowerLaw
      write (6,*) 'Number of Crank-Nicholson iterations (Marine)',
     1                                        nGSMarine
      
      write (6,*) 'Timing:'
      if (runSPL) write (6,*) 'SPL:',timeSPL
      if (runDiffusion) write (6,*) 'Diffusion:',timeDiffusion
      if (runMarine) write (6,*) 'Marine:',timeMarine
      if (runAdvect) write (6,*) 'Advection:',timeAdvect
      if (runUplift) write (6,*) 'Uplift:',timeUplift
      
      return
      end subroutine Debug
      
C---------------------------------------------------------------
      subroutine SetBC (jbc,bon,bon_1,bon_2,len_nodes,len_1)
      implicit none
      
      INTEGER :: I,bon_1,bon_2,len_1
      integer, intent(in) :: jbc
      integer bon(bon_1,bon_2)
      character :: cbc*4
      double precision len_nodes(bon_1,len_1)
      
      bounds_ibc = jbc
      bond_num(1:nn,1:nn1)=bon(1:bon_1,1:bon_2)
      length2(1:nn,1:2)=len_nodes(1:nn,1:2)
      
      write (cbc,'(i4)') jbc
      bounds_bc=.FALSE.
c      bounds_i1=1
c      bounds_i2=nx
c      bounds_j1=1
c      bounds_j2=ny
c      if (cbc(4:4).eq.'1') bounds_i1=2
c      if (cbc(2:2).eq.'1') bounds_i2=nx-1
c      if (cbc(1:1).eq.'1') bounds_j1=2
c      if (cbc(3:3).eq.'1') bounds_j2=ny-1
      DO I=1,NN
      if ((cbc(4:4).EQ.'1').AND.(BON(I,1).EQ.4)) bounds_bc(I)=.TRUE.
      if ((cbc(3:3).EQ.'1').AND.(BON(I,1).EQ.3)) bounds_bc(I)=.TRUE.
      if ((cbc(2:2).EQ.'1').AND.(BON(I,1).EQ.2)) bounds_bc(I)=.TRUE.
      if ((cbc(1:1).EQ.'1').AND.(BON(I,1).EQ.1)) bounds_bc(I)=.TRUE.
      END DO
C      if (cbc(4:4).eq.'1') bounds_bc(1:nn:nx)=.TRUE.
C      if (cbc(2:2).eq.'1') bounds_bc(nx:nn:nx)=.TRUE.
C      if (cbc(1:1).eq.'1') bounds_bc(1:nx)=.TRUE.
C      if (cbc(3:3).eq.'1') bounds_bc(nx*(ny-1)+1:nn)=.TRUE.
      bounds_xcyclic=.FALSE.
      bounds_ycyclic=.FALSE.
      if (cbc(4:4).ne.'1'.and.cbc(2:2).ne.'1') bounds_xcyclic=.TRUE.
      if (cbc(1:1).ne.'1'.and.cbc(3:3).ne.'1') bounds_ycyclic=.TRUE.
      return
      
      end subroutine SetBC
      
C---------------------------------------------------------------
      subroutine SetU (up,bon,bon_1,bon_2)
      double precision :: up(*)
      integer bon_1,bon_2
      integer bon(bon_1,bon_2)
      integer i
      
      runUplift = .true.
      u(:nn) = up(:nn)
cc      DO I=1,nn
cc      IF (bon(I,1).EQ.4) u(I)=-20/1000
cc      IF (bon(I,1).EQ.3) u(I)=-0.2/1000
cc      IF (bon(I,1).EQ.2) u(I)=-20/1000
cc      IF (bon(I,1).EQ.1) u(I)=-0.2/1000
cc      END DO
      
      return
      
      end subroutine SetU
      
C---------------------------------------------------------------
c      subroutine SetV (ux,uy)
c      implicit none
c      double precision, intent(in) :: ux(*),uy(*)
c      integer i
c      
c      runAdvect = .true.
c      do i=1,nn
c        vx(i) = ux(i)
c        vy(i) = uy(i)
c      enddo
c      return
c      
c      end subroutine SetV
c      
cC---------------------------------------------------------------
c      subroutine SetH (hp)
c      double precision, intent(in), dimension(*) :: hp
c      
c      h = hp(1:nn)
c      return
c      
c      end subroutine SetH
c      
C---------------------------------------------------------------
      subroutine SetPrecip (precipp)
      double precision, intent(in), dimension(*) :: precipp
      
      precip(1:nn) = precipp(1:nn)
      return
      
      end subroutine SetPrecip
      
cC---------------------------------------------------------------
c      subroutine SetBasement (bp)
c      double precision, intent(in), dimension(*) :: bp
c      
c      b = bp(1:nn)
c      return
c      
c      end subroutine SetBasement
c      
cC---------------------------------------------------------------
c      subroutine compute_fluxes (tectonic_flux, erosion_flux,
c     1                            boundary_flux)
c      implicit none
c      double precision, intent(out) :: tectonic_flux, erosion_flux,
c     1                                boundary_flux
c      double precision :: surf
c      double precision, dimension(:), allocatable :: flux
c      integer ij,ijk,k
c      
c      surf = xl*yl/(nx - 1)/(ny - 1)
c      tectonic_flux = sum(u)*surf
c      erosion_flux = sum(erate)*surf
c      
c      ! computes receiver and stack information for multi-direction flow
c      !allocate (mrec(10,nn), mnrec(nn), mwrec(10,nn), mlrec(10,nn), mstack(nn), hwater(nn)
c      !call find_mult_rec (h, rec, stack, hwater, mrec, mnrec, mwrec, mlrec, mstack, nx, ny, xl/(nx-1), yl/(ny-1), p, bounds, p_mfd_exp)
c      ! computes sediment flux
c      allocate (flux(nn))
c      flux = erate
c      
c      if (SingleFlowDirection) then
c        do ij = nn ,1, -1
c          ijk = stack(ij)
c          flux(ijk)=max(0.d0,flux(ijk))
c          flux(rec(ijk)) = flux(rec(ijk)) + flux(ijk)
c        enddo
c      else
c        do ij = 1, nn
c          ijk = mstack(ij)
c          flux(ijk)=max(0.d0,flux(ijk))
c          do k = 1, mnrec(ijk)
c            flux(mrec(k,ijk)) = flux(mrec(k,ijk)) + 
c     1                            flux(ijk)*mwrec(k,ijk)
c          enddo
c        enddo
c      endif
c      
c      ! compute boundary flux
c      boundary_flux = sum(flux,bounds_bc)*surf
c      deallocate (flux)
c      return
c      
c      end subroutine compute_fluxes
c      
cC---------------------------------------------------------------
c      subroutine set_tolerance (rel, abs, nGSSmax)
c      ! internal routine to set the relative and aboslute tolerance and maximum 
c      ! number of GSS iterations for the streamPowerLaw that includes sediment
c      implicit none
c      double precision :: rel, abs
c      integer :: nGSSmax
c      
c      tol_rel = rel		! GAUSS-SEIDEL迭代的相对容差（应用于当前最低地形高程）（建议：1e-4）
c      tol_abs = abs		! GAUSS-SEIDEL迭代的绝对容差（建议：1e-4）
c      nGSStreamPowerLawMax = nGSSmax		! GAUSS-SEIDEL迭代最大次数（建议：100）
c      return
c      
c      end subroutine set_tolerance
c        
cC---------------------------------------------------------------
c      subroutine get_nGSSiterations (nGSS)
c      ! internal routine to get the number of GSS iterations for the
c      ! StreamPowerLaw computations (with sediment)
c      implicit none
c      
c      integer :: nGSS
c      nGSS = nGSStreamPowerLaw
c      return
c      
c      end subroutine get_nGSSiterations
c      
cC---------------------------------------------------------------        
      end module FastScapeContext
 
      
      
      
      
      
      
      
      
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCC---FastScape_API---CCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine FastScape_Init()
      use FastScapeContext
      implicit none
      
      call Init()
      return
      
      end subroutine FastScape_Init
      
C-------------------------------------------------------------------------        
      subroutine FastScape_Set_NX_NY (nn2,liee)
      use FastScapeContext
      implicit none
      integer, intent(in) :: nn2,liee
      
      call SetNXNY (nn2,liee)
      return
      
      end subroutine FastScape_Set_NX_NY
      
C-------------------------------------------------------------------------    
      subroutine FastScape_Set_XL_YL (xlyl)
      use FastScapeContext
      implicit none
      integer, intent(in), dimension(*) :: xlyl
      
      call SetXLYL (xlyl(1:nn))
      return
    
      end subroutine FastScape_Set_XL_YL

C-------------------------------------------------------------------------     
      subroutine FastScape_Set_DT (dtt)
      use FastScapeContext
      implicit none
      double precision, intent(in) :: dtt
      
      call SetDT (dtt)
      return
      
      end subroutine FastScape_Set_DT
      
C-------------------------------------------------------------------------           
      subroutine FastScape_Init_H(HP,HPP2)
      use FastScapeContext
      implicit none
      double precision, intent(in), dimension(*) :: hp,hpp2
      
      call InitH(hp(1:nn),hpp2(1:nn))

      return
      
      end subroutine FastScape_Init_H

C-------------------------------------------------------------------------             
      subroutine FastScape_Set_Precip (precipp)
      use FastScapeContext
      implicit none
      double precision, intent(in), dimension(*) :: precipp
      
      call SetPrecip (precipp(1:nn))
      return
      
      end subroutine FastScape_Set_Precip
      
      
C-------------------------------------------------------------------------           
      subroutine FastScape_Setup()
      use FastScapeContext
      implicit none
      
      call SetUp()
      return
      
      end subroutine FastScape_Setup
      
C-------------------------------------------------------------------------          
      subroutine FastScape_Set_Erosional_Parameters (kkf,kkfsed,mm,nnn,
     1            kkd,kkdsed,gg1,gg2,pp)
      use FastScapeContext
      implicit none
      double precision:: kkf(10201),kkd(10201)
      double precision, intent(in) :: kkfsed,mm,nnn,kkdsed,gg1,gg2,pp
      
      call SetErosionalParam (kkf(1:nn),kkfsed,mm,nnn,kkd(1:nn),
     1                            kkdsed,gg1,gg2,pp)
      return
      
      end subroutine FastScape_Set_Erosional_Parameters
      
C-------------------------------------------------------------------------      
      subroutine FastScape_Set_U (up,bon,bon_11,bon_22)
      use FastScapeContext
      implicit none
      double precision, intent(in):: up(*)
      integer bon_11,bon_22
      integer bon(bon_11,bon_22)
      
      call SetU(up(1:nn),bon(1:bon_11,1:bon_22),bon_11,bon_22)
      return
      
      end subroutine FastScape_Set_U
      
C-------------------------------------------------------------------------        
      subroutine FastScape_Set_BC(jbc,bon,bon_1,bon_2,len_nodes,len_1)
      use FastScapeContext
      implicit none
      integer bon_1,bon_2,len_1
      integer, intent(in) :: jbc
      integer bon(bon_1,bon_2)
      double precision len_nodes(bon_1,len_1)
      
      call SetBC (jbc,bon(1:nn,1:bon_2),bon_1,bon_2,
     1            len_nodes(1:nn,1:len_1),len_1)
      return
      
      end subroutine FastScape_Set_BC
      
C-------------------------------------------------------------------------   
      subroutine FastScape_Get_Step (sstep)
      use FastScapeContext
      implicit none
      integer, intent(out) :: sstep
      
      call GetStep (sstep)
      return
      
      end subroutine FastScape_Get_Step
      
C-------------------------------------------------------------------------        
      subroutine FastScape_Execute_Step()
      use FastScapeContext
      implicit none
      real :: time_in, time_out
      
C      if (runAdvect) then
C        call cpu_time (time_in)
C        call Advect ()
C        call cpu_time (time_out)
C        timeAdvect = timeAdvect + time_out-time_in
C      endif
      
      if (runUplift) then
        call cpu_time (time_in)
        call Uplift()
        call cpu_time (time_out)
        timeUplift = timeUplift + time_out-time_in
      endif
      
      if (runSPL) then
        call cpu_time (time_in)
        if (SingleFlowDirection) then
          call FlowRoutingSingleFlowDirection ()
          call FlowAccumulationSingleFlowDirection ()
          call StreamPowerLawSingleFlowDirection ()
        else
          call FlowRouting ()
          call FlowAccumulation ()
          call StreamPowerLaw ()
        endif
        call cpu_time (time_out)
        timeSPL = timeSPL + time_out-time_in
      endif

      if (runDiffusion) then
        call cpu_time (time_in)
        call Diffusion ()
        call cpu_time (time_out)
        timeDiffusion = timeDiffusion + time_out-time_in
      endif
      
c      if (runMarine) then
c         call cpu_time (time_in)
c         call Marine ()
c         call cpu_time (time_out)
c         timeMarine = timeMarine + time_out-time_in
c      endif
c     
c      if (runStrati) then
c         call cpu_time (time_in)
c         call Run_Strati ()
c         call cpu_time (time_out)
c         timeStrati = timeStrati + time_out-time_in
c      endif
      
      step=step+1
      return
      
      end subroutine FastScape_Execute_Step
      
C-------------------------------------------------------------------------              
      subroutine FastScape_Copy_Chi (chip)
      use FastScapeContext
      implicit none
      double precision, intent(inout), dimension(*) :: chip
      
      call CopyChi(chip(1:nn))
      return
      
      end subroutine FastScape_Copy_Chi
      
C-------------------------------------------------------------------------       
      subroutine FastScape_Copy_H(hp)
      use FastScapeContext
      implicit none
      double precision, intent(inout), dimension(*) :: hp
      
      call CopyH(hp(1:nn))
      return
      
      end subroutine FastScape_Copy_H

C-------------------------------------------------------------------------       
      subroutine FastScape_Copy_HH2(hp2)
      use FastScapeContext
      implicit none
      double precision, intent(inout), dimension(*) :: hp2
      
      call CopyHH2(hp2(1:nn))
      return
      
      end subroutine FastScape_Copy_HH2
      
C-------------------------------------------------------------------------         
      subroutine FastScape_Debug()
      use FastScapeContext
      implicit none

      call Debug()
      return

      end subroutine FastScape_Debug
      
C-------------------------------------------------------------------------           
      subroutine FastScape_Destroy()
      use FastScapeContext
      implicit none
      
      call Destroy()
      return
      
      end subroutine FastScape_Destroy
      
C-------------------------------------------------------------------------        
      subroutine FastScape_Copy_Drainage_Area (ap)
      use FastScapeContext
      implicit none
      
      double precision, intent(inout), dimension(*) :: ap
      call CopyArea(ap)
      return
      
      end subroutine FastScape_Copy_Drainage_Area
      
C-------------------------------------------------------------------------          
      subroutine FastScape_Copy_Total_Erosion (etotp)
      use FastScapeContext
      implicit none
      double precision, intent(inout), dimension(*) :: etotp
      
      call CopyEtot(etotp)
      return

      end subroutine FastScape_Copy_Total_Erosion
      
C-------------------------------------------------------------------------       
      subroutine FastScape_Copy_Erosion_Rate (eratep)
      use FastScapeContext
      implicit none
      double precision, intent(inout), dimension(*) :: eratep

      call CopyERate(eratep)
      return

      end subroutine FastScape_Copy_Erosion_Rate
      
C-------------------------------------------------------------------------        
      subroutine FastScape_Copy_Catchment (catchp)
      use FastScapeContext
      implicit none
      double precision, intent(inout), dimension(*) :: catchp
      
      call CopyCatchment (catchp)
      return

      end subroutine FastScape_Copy_Catchment
      
C-------------------------------------------------------------------------      
      subroutine FastScape_Copy_Basement(bp)
      use FastScapeContext
      implicit none
      double precision, intent(inout), dimension(*) :: bp

      call CopyBasement(bp)
      return

      end subroutine FastScape_Copy_Basement
   
C-------------Read file1 information  
      SUBROUTINE READ_FILE1(PATH,N1,N2,ATA2)
      INTEGER N1,N2,I
      DOUBLE PRECISION ATA2(N1,N2)
      CHARACTER*256 PATH
      
      ATA2=0
      OPEN(101,FILE=TRIM(PATH),STATUS='OLD')
      DO I=1,N1
      READ(101,*)  ATA2(I,1:N2)
      END DO
      CLOSE(101)
      
      RETURN
      END SUBROUTINE

C-------------Read file2 information
      SUBROUTINE READ_FILE2(PATH,ATA,NUM_NODE1,NUM_NODE2)
      INTEGER NUM_NODE1,NUM_NODE2,I
      INTEGER ATA(NUM_NODE1,NUM_NODE2)
      CHARACTER*256 PATH
      
      ATA=0
      OPEN(102,FILE=TRIM(PATH),STATUS='OLD')
      DO I=1,NUM_NODE1
      READ(102,*)  ATA(I,1:NUM_NODE2)
      END DO
      CLOSE(102)
      
      RETURN
      END SUBROUTINE
      
C-------------Read file3 information
      SUBROUTINE READ_FILE3(PATH,TA1,TA2,TA3,TA4)
      INTEGER TA1,TA2,TA3,TA4
      DOUBLE PRECISION TEMP(3)
      CHARACTER*256 PATH
      
      OPEN(103,FILE=TRIM(PATH),STATUS='OLD')
      READ(103,*)  TA1,TA2,TA3,TA4
      CLOSE(103)
      
      RETURN
      END SUBROUTINE
      
C-------------读取文件3信息      
      SUBROUTINE Print_ste(KI,KI2,TIM,TEMP1,TEMP2,TEMP3,TEMP4,
     1                                            TEMP5,TEMP6,TEMP7)
      USE TEMPARRAY_1
      
      INTEGER I2,KI,KI2
      double precision TIM
      double precision, intent(in), dimension(*) :: TEMP1,TEMP2,TEMP3,
     1                                        TEMP4,TEMP5,TEMP6,TEMP7
      
      OPEN(104,FILE=PATH3,STATUS='OLD',POSITION='APPEND',
     1                                            ACTION='READWRITE')
      OPEN(105,FILE=PATH4,STATUS='OLD',POSITION='APPEND',
     1                                            ACTION='READWRITE')
      OPEN(106,FILE=PATH5,STATUS='OLD',POSITION='APPEND',
     1                                            ACTION='READWRITE')
      OPEN(107,FILE=PATH6,STATUS='OLD',POSITION='APPEND',
     1                                            ACTION='READWRITE')
      OPEN(108,FILE=PATH7,STATUS='OLD',POSITION='APPEND',
     1                                            ACTION='READWRITE')
      OPEN(109,FILE=PATH8,STATUS='OLD',POSITION='APPEND',
     1                                            ACTION='READWRITE')
      OPEN(110,FILE=PATH9,STATUS='OLD',POSITION='APPEND',
     1                                            ACTION='READWRITE')
      WRITE(104,*) 'KINC:', KI,TIM,KI2
      WRITE(105,*) 'KINC:', KI,TIM,KI2
      WRITE(106,*) 'KINC:', KI,TIM,KI2
      WRITE(107,*) 'KINC:', KI,TIM,KI2
      WRITE(108,*) 'KINC:', KI,TIM,KI2
      WRITE(109,*) 'KINC:', KI,TIM,KI2
      WRITE(110,*) 'KINC:', KI,TIM,KI2
      DO I2=1,NUM_NODE
      WRITE(104,*)  MODEL_COOR(I2,2:3),TEMP7(I2)
      WRITE(105,*)  MODEL_COOR(I2,2:3),TEMP1(I2)
      WRITE(106,*)  MODEL_COOR(I2,2:3),TEMP2(I2)
      WRITE(107,*)  MODEL_COOR(I2,2:3),TEMP3(I2)
      WRITE(108,*)  MODEL_COOR(I2,2:3),TEMP4(I2)
      WRITE(109,*)  MODEL_COOR(I2,2:3),TEMP5(I2)
      WRITE(110,*)  MODEL_COOR(I2,2:3),TEMP6(I2)
      END DO
      CLOSE(104)
      CLOSE(105)
      CLOSE(106)
      CLOSE(107)
      CLOSE(108)
      CLOSE(109)
      CLOSE(110)
      
      RETURN
      END SUBROUTINE
      
C-------------寻找模型角点      
      SUBROUTINE CORNER_POINT(TEMP1,TEMP2,TEMP22,TEMP3,NUM)
      double precision, intent(in), dimension(*) :: TEMP1
      integer, intent(in), dimension(*) :: TEMP2,TEMP22
      integer, intent(inout), dimension(*) :: TEMP3
      integer NUM,I,I1,max_pos(1),min_pos(1),index_temp(4,2),
     1        TEMPN,TEMPN2
      double precision TEMP(4),MAX_TEMP,MAX_TEMP2
      
      I1=0
      DO I=1,NUM
      IF (TEMP22(I).EQ.3) THEN
      I1=I1+1
      TEMP(I1)=TEMP1(I)
      index_temp(I1,1)=I
      index_temp(I1,2)=TEMP2(I)
      END IF
      END DO
      
      MAX_TEMP=-9999999999
      MAX_TEMP2=-9999999999
      DO I=1,4
      IF (index_temp(I,2).EQ.4) THEN
      IF (TEMP(I).GT.MAX_TEMP) THEN
      TEMPN=I
      MAX_TEMP=TEMP(I)
      END IF
      ELSEIF (index_temp(I,2).EQ.2) THEN
      IF (TEMP(I).GT.MAX_TEMP2) THEN
      TEMPN2=I
      MAX_TEMP2=TEMP(I)
      END IF
      END IF 
      END DO
      
      TEMP3(1)=index_temp(TEMPN,1)
      TEMP3(3)=index_temp(TEMPN2,1)
      
      MAX_TEMP=9999999999
      MAX_TEMP2=9999999999
      DO I=1,4
      IF (index_temp(I,2).EQ.4) THEN
      IF (TEMP(I).LT.MAX_TEMP) THEN
      TEMPN=I
      MAX_TEMP=TEMP(I)
      END IF
      ELSEIF (index_temp(I,2).EQ.2) THEN
      IF (TEMP(I).LT.MAX_TEMP2) THEN
      TEMPN2=I
      MAX_TEMP2=TEMP(I)
      END IF
      END IF 
      END DO
      
      TEMP3(2)=index_temp(TEMPN,1)
      TEMP3(4)=index_temp(TEMPN2,1)
      
      RETURN
      END SUBROUTINE
      
C-------------
      SUBROUTINE FIX_OCEAN(TEMP1,TEMP2,TEMP3,BONC,NUM,H_OCE)
      INTEGER TEMP1(NUM,2)
      DOUBLE PRECISION TEMP2(NUM),TEMP3(NUM)
      INTEGER :: I,NUM,BONC,LEV1,LEV2,LEV3,LEV4
      DOUBLE PRECISION MIN_value,MAX_value,H_OCE

      LEV1 = MOD(BONC/1000,10)    ! Down
      LEV2 = MOD(BONC/100,10)     ! Right
      LEV3 = MOD(BONC/10,10)      ! Up 
      LEV4 = MOD(BONC,10)         ! Left
      
      MIN_value = MINVAL(TEMP2)
      MAX_value = MAXVAL(TEMP2)
      
      DO I=1,NUM
      IF ((TEMP1(I,1).EQ.1).AND.(LEV1.EQ.1)) TEMP3(I)=H_OCE
      IF ((TEMP1(I,1).EQ.2).AND.(LEV2.EQ.1)) TEMP3(I)=H_OCE
      IF ((TEMP1(I,1).EQ.3).AND.(LEV3.EQ.1)) TEMP3(I)=H_OCE
      IF ((TEMP1(I,1).EQ.4).AND.(LEV4.EQ.1)) TEMP3(I)=H_OCE
      IF ((LEV1.EQ.1).AND.(ABS(TEMP2(I)-MIN_value).LE.10))TEMP3(I)=H_OCE
      IF ((LEV3.EQ.1).AND.(ABS(TEMP2(I)-MAX_value).LE.10))TEMP3(I)=H_OCE
      END DO
      
      RETURN
      END SUBROUTINE
      
C-------------
      SUBROUTINE Uneroded_boundary(TEMP1,TEMP2,TEMP3,TEMP4,BONC,NUM)
      INTEGER TEMP1(NUM,2)
      DOUBLE PRECISION TEMP2(NUM),TEMP3(NUM),TEMP4(NUM)
      INTEGER :: I,NUM,BONC,LEV1,LEV2,LEV3,LEV4
      DOUBLE PRECISION MIN_value,MAX_value

      LEV1 = MOD(BONC/1000,10)    ! Down
      LEV2 = MOD(BONC/100,10)     ! Right
      LEV3 = MOD(BONC/10,10)      ! Up 
      LEV4 = MOD(BONC,10)         ! Left
      
      MIN_value = MINVAL(TEMP2)
      MAX_value = MAXVAL(TEMP2)
      
      DO I=1,NUM
      IF ((TEMP1(I,1).EQ.1).AND.(LEV1.EQ.1)) TEMP3(I)=TEMP4(I)
      IF ((TEMP1(I,1).EQ.2).AND.(LEV2.EQ.1)) TEMP3(I)=TEMP4(I)
      IF ((TEMP1(I,1).EQ.3).AND.(LEV3.EQ.1)) TEMP3(I)=TEMP4(I)
      IF ((TEMP1(I,1).EQ.4).AND.(LEV4.EQ.1)) TEMP3(I)=TEMP4(I)
      IF ((LEV1.EQ.1).AND.(ABS(TEMP2(I)-MIN_value).LE.10))
     1    TEMP3(I)=TEMP4(I)
      IF ((LEV3.EQ.1).AND.(ABS(TEMP2(I)-MAX_value).LE.10))
     1    TEMP3(I)=TEMP4(I)
      END DO
      
      RETURN
      END SUBROUTINE
      
C-------------计算节点的高程抬升
      subroutine Uplift ()
C      subroutine to apply an uplift step using the uplift function/array u
C      设置高程抬升速率(nx*ny)（m/yr）--固定边界也可以设置抬升速率
      use FastScapeContext
      implicit none
      
      h(1:nn) = h(1:nn) + u(1:nn)*dt	! u是单节点的抬升率
      b(1:nn) = b(1:nn) + u(1:nn)*dt
      return
      
      end subroutine Uplift

C---------------------------------------------------------------------------------------------------------
      subroutine FlowRouting ()
      use FastScapeContext
      implicit none
      
      ! finds receiver
      call find_receiver (h,nn,nn1,rec,length,bounds_xcyclic,
     2                    bounds_ycyclic,bounds_bc,bond_num(1:nn,1:nn1),
     3                    length2(:nn,:2))
      
      ! finds donors
      call find_donor (rec, nn, ndon, don)
      
      ! find stack
      call find_stack (rec, don, ndon, nn, catch0, stack, catch)
      
      ! removes local minima
      call LocalMinima (stack,rec,bounds_bc,ndon,don,h,
     2                    length,nn,nn1,bond_num,length2(:nn,:2))
      
      ! computes receiver and stack information for mult-direction flow （计算后的mstack是由顶点到低点）
      call find_mult_rec (h,rec,stack,hwater,mrec,mnrec,mwrec,mlrec,
     1                mstack,nn,nn1,p,p_mfd_exp,bounds_xcyclic,
     3                bounds_ycyclic,bounds_bc,bond_num,length2(:nn,:2))
      
      ! compute lake depth
      lake_depth(1:nn) = hwater(1:nn) - h(1:nn)
      return
      
      end subroutine FlowRouting
      
C--------------------------------------------------------------------------------------------
      subroutine FlowRoutingSingleFlowDirection ()
      use FastScapeContext
      implicit none
      integer :: i, ijk, ijr
      double precision :: deltah
      integer, dimension(nn,nn1):: bond_num1
      
      do i=1,nn
      bond_num1(i,1:nn1)=bond_num(i,1:nn1)
      end do
      ! finds receiver
      call find_receiver (h,nn,nn1,rec,length,bounds_xcyclic,
     2                    bounds_ycyclic,bounds_bc,
     3                    bond_num1(1:nn,1:nn1),length2(:nn,:2))
      
      ! finds donors
      call find_donor (rec, nn, ndon, don)
      
      ! find stack
      call find_stack (rec, don, ndon, nn, catch0, stack, catch)
      
      ! removes local minima
      call LocalMinima (stack,rec,bounds_bc,ndon,don,h,
     2                    length,nn,nn1,bond_num,length2(:nn,:2))
      
      ! find hwater
      hwater(1:nn) = h(1:nn)
      
      ! fill the local minima with a nearly planar surface
      deltah = 1.d-8
      do i=1,nn
        ijk = stack(i)
        ijr = rec(ijk)
        if (ijr.ne.0) then
          if (hwater(ijr).gt.hwater(ijk)) then
            hwater(ijk) = hwater(ijr) + deltah
          endif
        endif
      enddo
      
      ! compute lake depth
      lake_depth(1:nn) = hwater(1:nn) - h(1:nn)
      return
      
      end subroutine FlowRoutingSingleFlowDirection
      
C--------------------------------------------------------------------------------------------
      subroutine FlowAccumulation () ! 多流向计算河流流量
      use FastScapeContext
      implicit none
      integer :: ij, ijk, k
      
      a(1:nn)=sear(1:nn)*precip(1:nn) ! 采用方格面积计算降雨量，precip是降雨速率m/yr
      do ij=1,nn
        ijk=mstack(ij)
        do k =1,mnrec(ijk)
          a(mrec(k,ijk))=a(mrec(k,ijk))+a(ijk)*mwrec(k,ijk)
        enddo
      enddo
      return
      
      end subroutine FlowAccumulation
      
C--------------------------------------------------------------------------------------------
      subroutine FlowAccumulationSingleFlowDirection ()   ! 单流向计算河流流量
      use FastScapeContext
      implicit none
      integer :: ij,ijk,k
      
      a(1:nn)=sear(1:nn)*precip(1:nn) ! 采用方格面积计算降雨量，precip是降雨速率m/yr
      do ij=nn,1,-1
        ijk=stack(ij)
        a(rec(ijk))=a(rec(ijk))+a(ijk)
      enddo
      
      return
      
      end subroutine FlowAccumulationSingleFlowDirection
      
C--------------------------------------------------------------------------------------------
      subroutine find_mult_rec (h,rec0,stack0,water,rec,nrec,wrec,lrec,
     1                   stack,nn,nn1,p,p_mfd_exp,bounds_xcyclic,
     3                   bounds_ycyclic,bounds_bc,bond_num,length22)
      
      ! subroutine to find multiple receiver information
      ! in input:
      ! h is topography
      ! rec0 is single receiver information
      ! stack0 is stack (from bottom to top) obtained by using single receiver information
      ! water is the surface of the lakes (or topography where there is no lake)
      ! nx, ny resolution in x- and y-directions
      ! dx, dy grid spacing in x- and y-directions
      ! p is exponent to which the slope is put to share the water/sediment among receivers
      ! bounds: boundaries type (boundary conditions)
      ! in output:
      ! rec: multiple receiver information
      ! nrec: number of receivers for each node
      ! wrec: weight for each receiver
      ! lrec: distance to each receiver
      ! stack: stoack order for multiple receivers (from top to bottom)
      
      logical, intent(in) :: bounds_xcyclic, bounds_ycyclic
      integer nn,ij2,nn1
      double precision h(nn),wrec(10,nn),lrec(10,nn),p,
     1                    water(nn),p_mfd_exp(nn)
      integer rec(10,nn),nrec(nn),stack(nn),rec0(nn),
     1                stack0(nn)
      integer :: i,j,ii,jj,iii,jjj,ijk,k,ijr,nparse,nstack,ijn,ij
      double  precision :: slopemax,sumweight,deltah,slope
      integer, dimension(:), allocatable :: ndon,vis,parse
      integer, dimension(:,:), allocatable :: don
      double precision, dimension(:), allocatable :: h0
      double precision, dimension(nn,2), intent(in):: length22
      integer, dimension(nn,nn1), intent(in) :: bond_num
      logical, dimension(nn), intent(in):: bounds_bc
      allocate (h0(nn))
      
      h0(1:nn)=h(1:nn)
      
      ! fill the local minima with a nearly planar surface
      deltah = 1.d-8	! 填充湖水的量
      do i=1,nn
        ijk=stack0(i)
        ijr=rec0(ijk)
        if (ijr.ne.0) then
          if (h0(ijr).gt.h0(ijk)) then	! 前面利用假高判断的可能的内盆地水流向其他盆地
            h0(ijk)=h0(ijr)+deltah
          endif
        endif
      enddo
      
      water(1:nn) = h0(1:nn)
      nrec=0
      wrec=0.d0
      ! loop on all nodes 寻找节点接收者（abaqus修改）
      do ij=1,nn
      if (bounds_bc(ij)) cycle
      slopemax = 0.
      bi=bond_num(ij,2)
      do jj=1,bi
      ijk = bond_num(ij,2+jj)
      if (h0(ij).gt.h0(ijk)) then	! 与之前的判别不同，不是判别最低点，而是判别所有比节点低的接收点
        nrec(ij)=nrec(ij)+1	! 记录有几个可能的低接收点
        rec(nrec(ij),ij) = ijk
        lrec(nrec(ij),ij) = sqrt((length22(ij,1)-length22(ijk,1))**2+
     1            (length22(ij,2)-length22(ijk,2))**2)
        wrec(nrec(ij),ij) = (h0(ij) - h0(ijk))/lrec(nrec(ij),ij)
      endif
      enddo
      enddo
      
      do ij =1,nn
        if (p<0.d0) then	! 多向流分配系数；当p<0时，p = 0.5 + 0.6*slope
          slope = 0.d0
          if (nrec(ij).ne.0) 
     1        slope = real(sum(wrec(1:nrec(ij),ij))/nrec(ij))
          p_mfd_exp(ij) = 0.5 + 0.6*slope
        endif
        do k=1,nrec(ij)
          wrec(k,ij) = wrec(k,ij)**p_mfd_exp(ij)
        enddo
        sumweight = sum(wrec(1:nrec(ij),ij))
        do ij2=1,nrec(ij)     !计算分配系数
        wrec(ij2,ij) = wrec(ij2,ij)/sumweight
        end do
      enddo
      
      allocate (ndon(nn),don(10,nn))
      
      ndon=0
      
      ! 记录分支节点接收来自几个节点的分配
      do ij=1,nn
        do k=1,nrec(ij)
          ijk = rec(k,ij)
          ndon(ijk)=ndon(ijk)+1
          don(ndon(ijk),ijk) = ij	
        enddo
      enddo
      
      allocate (vis(nn),parse(nn))
      
      nparse=0
      nstack=0
      stack=0
      vis=0
      parse=0
      
      ! we go through the nodes
      do ij=1,nn
        ! when we find a "summit" (ie a node that has no donors)
        ! we parse it (put it in a stack called parse)
        if (ndon(ij).eq.0) then	! 记录顶峰节点（即无贡献者的节点），并放入parse堆栈中
          nparse=nparse+1
          parse(nparse)=ij
        endif
        ! we go through the parsing stack
        do while (nparse.gt.0)
          ijn=parse(nparse)
          nparse=nparse-1
          ! we add the node to the stack
          nstack=nstack+1
          stack(nstack)=ijn		! stack堆栈从顶点开始记录
          ! for each of its receivers we increment a counter called vis
          do  ijk=1,nrec(ijn)
            ijr=rec(ijk,ijn)
            vis(ijr)=vis(ijr)+1
            ! if the counter is equal to the number of donors for that node we add it to the parsing stack
            if (vis(ijr).eq.ndon(ijr)) then
              nparse=nparse+1
              parse(nparse)=ijr
            endif
          enddo
        enddo
        enddo
        
      if (nstack.ne.nn) then
      write(6,*) 'error in stack'
      end if
      
c      deallocate (ndon,don,vis,parse,h0)
      if (allocated(ndon)) deallocate(ndon)
      if (allocated(don)) deallocate(don)
      if (allocated(vis)) deallocate(vis)
      if (allocated(parse)) deallocate(parse)
      if (allocated(h0)) deallocate(h0)
      
      return

      end subroutine find_mult_rec
      
C--------------------------------------------------------------------------------------------
      subroutine find_receiver (h,nn,nn1,rec,length,bounds_xcyclic,
     1                        bounds_ycyclic,bounds_bc,bond_num1,
     2                        length22)
      
      implicit none
      logical, intent(in) :: bounds_xcyclic, bounds_ycyclic
      integer, intent(in) :: nn,nn1
      double precision, dimension(nn), intent(in) :: h
      integer, dimension(nn), intent(out) :: rec
      double precision, dimension(nn), intent(out):: length
      double precision, dimension(nn,2), intent(in):: length22
      integer, dimension(nn,nn1):: bond_num1
      logical, dimension(nn):: bounds_bc
      integer :: ij, i, j, ii, jj, ijk,bi
      double precision :: l, smax, slope
      
      ! resets receiver and distance between node and its receiver
      do ij=1,nn
        rec(ij)=ij
      enddo
      
      ! finds receiver using steepest descent/neighbour method
      do ij=1,nn
      if (bounds_bc(ij)) cycle
          smax=tiny(smax)
          bi=bond_num1(ij,2)   !bond/num_node/node1/node2/node3...
          if ((bounds_xcyclic.and.
     1        ((bond_num1(ij,1).eq.2).or.
     2        (bond_num1(ij,1).eq.4)))
     3        .or.(bounds_ycyclic.and.
     4        ((bond_num1(ij,1).eq.1).or.
     5        (bond_num1(ij,1).eq.3)))) then
              jj=0
              do while(bond_num1(ij,jj+3).ne.0)
              jj=jj+1
              if (((jj.eq.10).and.(nn1.eq.12)).or.
     1                    ((jj.eq.8).and.(nn1.eq.10))) exit
              end do
              bi=jj
          end if
          do ii=1,bi
              ijk=bond_num1(ij,2+ii)
              l=sqrt((length22(ij,1)-length22(ijk,1))**2+
     1                    (length22(ij,2)-length22(ijk,2))**2)
              slope=(h(ij)-h(ijk))/l
              if (slope.gt.smax) then
                smax=slope
                rec(ij)=ijk
                length(ij)=l
              endif
            enddo
      end do
      return
      
      end subroutine find_receiver
      
C--------------------------------------------------------------------------------------------
      subroutine find_donor (rec, nn, ndon, don)
      implicit none
      integer, intent(in) :: nn
      integer, dimension(nn), intent(in) :: rec
      integer, dimension(nn), intent(out) :: ndon
      integer, dimension(10,nn), intent(out) :: don
      integer ij, ijk
      
      ! inverts receiver array to compute donor arrays
      ndon=0
      do ij=1,nn
        if (rec(ij).ne.ij) then
          ijk=rec(ij)
          ndon(ijk)=ndon(ijk)+1
          don(ndon(ijk),ijk)=ij
        endif
      enddo
      return
      
      end subroutine find_donor
      
C--------------------------------------------------------------------------------------------
      subroutine find_stack (rec, don, ndon, nn, catch0, stack, catch)
      implicit none
      integer, intent(in) :: nn
      integer, intent(in), dimension(nn) :: rec, ndon
      integer, intent(in), dimension(10,nn) :: don
      double precision, intent(in), dimension(nn) :: catch0
      integer, intent(out), dimension(nn) :: stack
      double precision, dimension(nn), intent(out) :: catch
      
      integer :: ij, nstack
      
      ! computes stack by recursion
      nstack=0
      catch(:nn)=catch0(:nn)
      do ij=1,nn
      if (rec(ij).eq.ij) then
      nstack=nstack+1
      stack(nstack)=ij
      call find_stack_recursively (ij,don,ndon,nn,stack,nstack,catch)
      endif
      enddo
      return
      
      end subroutine find_stack
      
C--------------------------------------------------------------------------------------------
      recursive subroutine find_stack_recursively (ij,don,ndon,nn,stack,
     1                                                nstack,catch)
      ! recursive routine to go through all nodes following donor information
      implicit none
      integer k,ij,ijk,nn,nstack
      integer don(10,nn),ndon(nn),stack(nn)
      double precision catch(nn)
      
      do k=1,ndon(ij)
        ijk=don(k,ij)
        nstack=nstack+1
        stack(nstack)=ijk
        catch(ijk)=catch(ij)
        call find_stack_recursively (ijk,don,ndon,nn,stack,nstack,catch)
      enddo
      return
      
      end subroutine find_stack_recursively
      
C--------------------------------------------------------------------------------------------
      subroutine LocalMinima (stack,rec,bc,ndon,donor,h,length,
     1                        nn,nn1,bond_num,length22)
      
      ! subroutine to compute and remove local inima by recomputing the receiver connectivty
      ! using Guillaume Cordonnier s algorithm as helped by Benoit Bovy and debuged with Jean Braun
      ! This is (at first) a fortran translation of a python routine provided to me by Gullaume
      ! on October 24 2017
      
      ! in input stack(nn) : stack order as defined by Braun and Willett (2013)
      !          rec(nn) : receiver node information as defined by Braun and Willet (2013)
      !          ndon(nn) : number of donors per node
      !          donor(10,nn) : list of donors per node
      !          h(nn) : nodal heigh field
      !          length(nn) : distance between node and receiver (per node)
      !          nx-ny : dimension of rectangular mesh in x- and y-directions
      !          dx-dy : grid spacing in the x- and y-directions
      
      ! in output stack : updated stack
      !           rec : updated receiver information
      !           ndon : updated number of donor information
      !           donor : updated donor list
      !           length : updated distance between donot and receiver
      
      ! The updated receiver list takes into account the existence of local minima in the original receiver list,
      ! i.e. nodes taht are their own receivers. The algorithm finds the geometry of lakes centered on these local
      ! minima as well as their sills. By connecting the lakes through their sills, the algorithm finds the steepest
      ! path (for all nodes) to a node on the boundary. A fake height arrays is created such that lakes
      ! are drained in the order given by the new receiver information. The algorithm is O(n) + (NlogN) ou N is the
      ! number of lakes
      
      ! Note that the user can use the routines loc_min_3__find_receivers, loc_min_3_find_donors and loc_min_3_find_stack
      ! provided below to compute the receiver, donor information and the stack before calling this routine.
      
      implicit none
      integer nn,k,nlocmin,nn1
      integer stack(nn),rec(nn),ndon(nn),donor(10,nn),maybe_basin(nn)
      double precision h(nn),length(nn)
      logical bc(nn)
      integer, dimension(nn,nn1), intent(in) :: bond_num
      integer, dimension(:), allocatable :: basins,outlets,
     1                                    tree,basin_stack
      integer, dimension(:,:), allocatable :: conn_basins,conn_nodes,
     1                                    sills,junki
      double precision, dimension(:), allocatable :: conn_weights,junkd
      logical, dimension(:), allocatable :: active_nodes
      logical continuous_flow,continuous_flow_v2
      integer nbasins,basin0,nconn,nconn_max,tree_size
      double precision, dimension(nn,2), intent(in):: length22
      
      
      !continuous_flow = .true.
      continuous_flow = .false.
      continuous_flow_v2 = .true.
      
      allocate(basins(nn),outlets(nn))
      call compute_basins (stack,rec,basins,outlets,nn,nbasins,
     1                    h,maybe_basin)
      
      nconn_max=nbasins*6
      allocate (conn_basins(nconn_max,2),conn_nodes(nconn_max,2),
     1            conn_weights(nconn_max))
      allocate(active_nodes(nn))
      active_nodes(1:nn)=.not.bc !流动节点是true
      
      nlocmin=0
      do k=1,nn
        if (rec(k).eq.k .and. active_nodes(k)) nlocmin=nlocmin+1 !内盆地数量
      enddo
      !if (nlocmin.gt.0) print*,'nlocmin before',nlocmin
      
      call connect_basins (nbasins,basins,outlets,rec,stack,
     1                    active_nodes,h,nn,nn1,basin0,nconn,
     2                    conn_basins,conn_nodes,conn_weights,
     3                    nconn_max,bond_num,maybe_basin)
      
      ! 将盆地连接的信息数组替换至刚好大小的数组内
      allocate (junki(nconn_max,2))
      junki(1:nconn_max,1:2)=conn_basins(1:nconn_max,1:2)
      if (allocated(conn_basins)) deallocate(conn_basins)
      allocate(conn_basins(nconn,2))
      conn_basins(1:nconn,:)=junki(1:nconn,:)
      junki(1:nconn_max,1:2)=conn_nodes(1:nconn_max,1:2)
      if (allocated(conn_nodes)) deallocate(conn_nodes)
      allocate(conn_nodes(nconn,2))
      conn_nodes(1:nconn,:)=junki(1:nconn,:)
      if (allocated(junki)) deallocate(junki)
      allocate (junkd(nconn_max))
      junkd(1:nconn_max)=conn_weights(1:nconn_max)
      if (allocated(conn_weights)) deallocate(conn_weights)
      allocate(conn_weights(nconn))
      conn_weights(1:nconn)=junkd(1:nconn)
      if (allocated(junkd)) deallocate(junkd)
      allocate (tree(nbasins-1))
      tree_size=0
      
      ! 建立最小生成树
      call mst_kruskal(conn_weights,conn_basins,nbasins,
     1                    nconn,tree,tree_size)
      
      allocate (sills(nbasins,2))
      allocate (basin_stack(nbasins))
      
      call order_tree (conn_basins(1:nconn,1),conn_basins(1:nconn,2),
     1                conn_nodes(1:nconn,1),conn_nodes(1:nconn,2),
     2                nconn,nbasins,basin0,tree,tree_size,sills,
     3                basin_stack,continuous_flow)
      
      if (.not.continuous_flow) then
        if (continuous_flow_v2) then
          call correct_receivers_v2 (rec,length,outlets,conn_basins,
     1            conn_nodes,tree,nn,nbasins,nconn,tree_size,
     2            length22(:nn,:2))
        else
          call correct_receivers (rec,length,outlets,conn_basins,
     1     conn_nodes,tree,h,nn,nbasins,nconn,tree_size,
     2    length22(:nn,:2))
        endif
      else
        allocate (junkd(nn))
        junkd(1:nn) = h(1:nn)	! 储存真实高程
        call update_fake_topography (sills,basin_stack,basins,h,nn,nn1,
     1                                nbasins,bond_num,length22(:nn,:2))
        call loc_min_3_find_receivers (h,rec,length,bc,nn,nn1,
     1                                    bond_num,length22(:nn,:2))	!利用假高程计算节点的接受者信息
        h(1:nn) = junkd(1:nn)	! 返回真实高程
        if (allocated(junkd)) deallocate(junkd)
      endif
      
      nlocmin=0
      do k=1,nn
        if (rec(k).eq.k.and.active_nodes(k)) nlocmin=nlocmin+1
      enddo
      if (nlocmin.gt.0) then
      write(6,*) 'nlocmin after:',nlocmin
      end if
      
cc      deallocate (basins,outlets,conn_basins,conn_nodes,conn_weights,
cc     1            active_nodes,tree,sills)
      if (allocated(basins)) deallocate(basins)
      if (allocated(outlets)) deallocate(outlets)
      if (allocated(conn_basins)) deallocate(conn_basins)
      if (allocated(conn_nodes)) deallocate(conn_nodes)
      if (allocated(conn_weights)) deallocate(conn_weights)
      if (allocated(active_nodes)) deallocate(active_nodes)
      if (allocated(junki)) deallocate(junki)
      if (allocated(junkd)) deallocate(junkd)
      if (allocated(tree)) deallocate(tree)
      if (allocated(sills)) deallocate(sills)
      if (allocated(basin_stack)) deallocate(basin_stack)
      
      call loc_min_3_find_donors (rec,ndon,donor,nn)
      call loc_min_3_find_stack (rec,ndon,donor,stack,nn)
      return
      
      end subroutine LocalMinima
      
C----------------------
      subroutine compute_basins (stack,rec,basins,outlets,n,nbasins,
     1                        h,maybe_basin)
      !Input:
      !stack: parse order from lower to upper nodes
      !receivers: array of recievers
      
      !Outputs:
      !basins: id of the basin for each node
      !outlets: minimal node of each basin

      !return nbasins: the number of basins
      
      implicit none
      integer n,nbasins
      integer stack(n),rec(n),basins(n),outlets(n),maybe_basin(n)
      integer ibasins,i,istack,irec,min_pos(1)
      double precision min1, max1, min2, max2
      double precision h(n),temp_h(n)
      
      ibasins=0
      maybe_basin=0
      temp_h=1e+20
      ! 循环寻找盆地（即接收者为自己的节点）
      do i=1,n
        istack=stack(i)
        irec=rec(istack)
        if (irec.eq.istack) then
          ibasins=ibasins+1
          outlets(ibasins)=istack
          temp_h(i)=h(i)
        endif
        basins(istack)=ibasins
        nbasins=ibasins
      enddo
      
      min_pos = minloc(temp_h(:n))
      maybe_basin(min_pos(1))=1
      return
      
      end subroutine compute_basins
      
C----------------------
      subroutine connect_basins (nbasins,basins,outlets,rec,stack,
     1                            active_nodes,h,n,nn1,basin0,nconn,
     2                            conn_basins,conn_nodes,conn_weights,
     3                            nconn_max,bond_num,maybe_basin)
      !Input:
      !nbasins: the number of basins
      !basins: id of the basin for each node
      !outlets: minimal node of each basin
      !receivers: array of recievers
      !stack: parse order from lower to upper nodes
      !active_nodes: nodes not considered as glabal minima (not base nodes)
      !elevation: elevation of each node
      !nx, ny: size of the dem
      
      !Output:
      !basin0: id of the root basin (the one representing the sea)
      !nconn: number of connections between basins
      !conn_basins: array of pairs of basins (b0, b1) that share a pass
      !conn_nodes: array of pairs (p0, p1) : nodes of the passes for each pair of basin
      !conn_weights : height of the passes (max of elevation[p0], elevation[p1])
      
      !Connect adjacent basins together through their lowest pass.
      
      !Creates an (undirected) graph of basins and their connections.
      
      !The following information is stored for each edge of the graph:
      
      !- the two grid nodes that togheter form the pass of lowest
      !  elevation found between the two basins ;
      !- a weight that corresponds to the elevation of the pass, i.e.,
      !  the highest elevation among the two nodes forming the pass.
      
      !Notes
      !-----
      !Connections between open basins are handled differently:
      
      !Instead of finding connections between adjacent basins,
      !virtual connections are added between one given basin
      !and all other basins.
      !This may save a lot of uneccessary computation, while it
      !ensures a connected graph (i.e., every node has at least
      !an edge), as required for applying minimum spanning tree
      !optimization.
      
      implicit none
      
      integer n,basin0,nconn,nconn_max,nn1,ace
      integer nbasins,basins(n),outlets(n),rec(n),stack(n),
     1        maybe_basin(n)
      integer conn_basins(nconn_max,2),conn_nodes(nconn_max,2)
      double precision conn_weights(nconn_max)
      logical active_nodes(n)
      double precision h(n)
      integer, dimension(n,nn1), intent(in) :: bond_num
      integer i,j,istack,iused,irec,iistack,iiused,ii,ki,kj,k
      integer ineighbor,ineighbor_basin,ineighbor_outlet
      integer iconn,ibasin,conn_pos_used_size,conn_idx
      integer, dimension(:), allocatable :: conn_pos,conn_pos_used
      logical active,all_true
      double precision weight
      
      all_true = all(active_nodes)    !判断是否边界全为0000
      if (all_true) then
        ace=1
      else
        ace=0
      end if
      
      ! theory of planar graph -> max nb. of connections known
      iconn=1
      basin0=-1
      ibasin=1
      
      allocate (conn_pos(nbasins),conn_pos_used(nbasins))

      conn_pos=-1
      conn_pos_used_size=0
      
      active=.False.
      
      ! king (D8) neighbor lookup
      
      do iistack=1,n
        istack=stack(iistack)
        irec=rec(istack)
        ! new basin
        if (irec.eq.istack) then
          ibasin=basins(istack) ! 盆地编号1、2....
          active=active_nodes(istack) ! 流动点是true/大海流入点是fault
          do iiused=1,conn_pos_used_size
            iused=conn_pos_used(iiused)
            conn_pos(iused)=-1	! 重置盆地连接许可
          enddo
          conn_pos_used_size=0  ! 重置连接许可的暂时储存数组
          if ((.not.active).or.((ace.eq.1).and.
     1                (maybe_basin(iistack).eq.1))) then	! 开放盆地间的连接
            !      print*,ibasin,basin0
            if (basin0.eq.-1) then
              basin0=ibasin
            else
              conn_basins(iconn,1)=basin0
              conn_basins(iconn,2)=ibasin
              conn_nodes(iconn,1)=-1
              conn_nodes(iconn,2)=-1
              conn_weights(iconn)=-1.d10	! 给予足够小的权重，使得计算最小生成树时两个开放盆地一定相连
              iconn=iconn+1
              !        print*,iconn
              if (ace.eq.1) cycle
            endif
          endif
        endif
      
        if (active) then
          ii=bond_num(istack,2)
          do k=1,ii	!识别节点的8或者10个方向关联节点--模型侧边和拐角的节点仅5/3个关联的节点
      
            ineighbor = bond_num(istack,k+2) ! 临近点编号
            ineighbor_basin = basins(ineighbor) ! 临近点所属的盆地编号
            ineighbor_outlet = outlets(ineighbor_basin) ! 临近点所属的盆地最低点
      
            ! skip same basin or already connected adjacent basin
            ! don't skip adjacent basin if it's an open basin
            if (ibasin.ge.ineighbor_basin.and.
     1            active_nodes(ineighbor_outlet)) goto 111	! 开放盆地（入海口）可以链接所有盆地
      
            weight = max(h(istack), h(ineighbor))
            conn_idx = conn_pos(ineighbor_basin)
      
            ! add new connection
            if (conn_idx .eq. -1) then	! 检查盆地连接许可，已被连接的盆地不可添加新的通道，而是更新已有连接（两个盆地间只能有一个通道）
              conn_basins(iconn,1) = ibasin
              conn_basins(iconn,2) = ineighbor_basin
              conn_nodes(iconn,1) = istack
              conn_nodes(iconn,2) = ineighbor
              conn_weights(iconn) = weight
              conn_pos(ineighbor_basin) = iconn
              iconn = iconn+1
              conn_pos_used_size = conn_pos_used_size+1
              conn_pos_used(conn_pos_used_size) = ineighbor_basin
              !update existing connection
            elseif (weight .lt. conn_weights(conn_idx)) then ! 更新已有连接（两个盆地间只能有一个通道），最低点设为连接通道
              conn_nodes(conn_idx,1) = istack
              conn_nodes(conn_idx,2) = ineighbor
              conn_weights(conn_idx) = weight
            endif
      
111           continue
          enddo
        endif
      enddo
      
      if (allocated(conn_pos)) deallocate(conn_pos)
      if (allocated(conn_pos_used)) deallocate(conn_pos_used)
      
      nconn = iconn - 1
      return
      
      end subroutine connect_basins
      
C----------------------
      subroutine mst_kruskal(conn_weights, conn_basins, nbasins, nconn,
     1                        mstree, mstree_size)
      
      !kruskal algorithm to compute a Minimal Spanning Tree
      
      !Input:
      !conn_basins: array of pairs of basins (b0, b1) that share a pass
      !conn_weights : height of the passes (max of elevation[p0], elevation[p1])
      !nbasins: number of basins
      
      !Output:
      !mstree: id of the connections in the minimal spanning tree
      
      implicit none
      integer nbasins,nconn
      double precision conn_weights(nconn)
      integer conn_basins(nconn,2)
      integer mstree(nbasins-1)
      integer, dimension(:), allocatable :: sort_id,parent,rank
      integer mstree_size,eid,eeid,f0,f1,b0,b1
      allocate (sort_id(nconn))
      
      mstree_size = 0
      
      ! sort edges
      call loc_min_3_indexx (nconn,conn_weights,sort_id)
      !print*,'weights',conn_weights(sort_id)
      
      allocate (parent(nbasins),rank(nbasins))
      call UnionFindInit (parent,rank,nbasins)
      
      do eeid=1,nconn
        eid = sort_id(eeid)
      
        b0 = conn_basins(eid, 1)
        b1 = conn_basins(eid, 2)
      
        call Unionfind (b0,parent,nbasins,f0)
        call Unionfind (b1,parent,nbasins,f1)
      
        if (f0 .ne. f1) then
          mstree_size = mstree_size+1
          mstree(mstree_size) = eid
          call DoUnion (b0,b1,parent,rank,nbasins)
        endif
      enddo
      
c      deallocate (parent,rank)
      if (allocated(parent)) deallocate(parent)
      if (allocated(rank)) deallocate(rank)
      if (allocated(sort_id)) deallocate(sort_id)
      return
      
      end subroutine mst_kruskal
      
C----------------------
      subroutine order_tree(edges_n1, edges_n2, edges_p1, edges_p2,
     1                        nconn, num_nodes, root,tree, ntree, 
     2                        sills, basin_stack, continuous_flow)
      
      !Swap node order in each edge to follow the inverse of the flow
      !Input:
      !edges_n1, edges_n2:  pairs of basins (b0, b1) that share a pass
      !edges_p1, edges_p2: nodes of the passes for each pair of basin
      !num_nodes: number of basins
      !root: root basin (the sea)
      !tree: tree that need ordering
      
      !Output:
      !None: edges_n1, edges_n2, edges_p1, edges_p2 are updated in place
      
      implicit none
      integer num_nodes,root,ntree,nconn
      integer edges_n1(nconn),edges_n2(nconn)
      integer edges_p1(nconn),edges_p2(nconn)
      integer tree(ntree), sills(num_nodes, 2)
      integer basin_stack(num_nodes)
      logical continuous_flow
      
      integer, dimension(:), allocatable :: nodes_connects_size,
     1                            nodes_connects_ptr,nodes_adjacency
      integer, dimension(:,:), allocatable :: stack
      
      integer i,ii,n1,n2,stack_size,nodes_adjacency_size,node,parent,
     1        edge_id,edge_swap, basin_stack_size
      
      !nodes connections
      allocate (nodes_connects_size(num_nodes),
     1            nodes_connects_ptr(num_nodes))
      nodes_connects_size=0
      
      !print*,'num_nodes,ntree,root',num_nodes,ntree,root
      
      !parse the edges to compute the number of edges per node
      do ii = 1,ntree
        i=tree(ii)
        nodes_connects_size(edges_n1(i)) = 
     1            nodes_connects_size(edges_n1(i))+1
        nodes_connects_size(edges_n2(i)) = 
     1            nodes_connects_size(edges_n2(i))+1
      enddo
      !write(*,'(a,1024i4)')'nodes_connect_size',nodes_connects_size(1:num_nodes)
      
      !compute the id of first edge in adjacency table
      nodes_connects_ptr(1) = 1
      do i=2,num_nodes
        nodes_connects_ptr(i) = nodes_connects_ptr(i-1) + 
     1                                nodes_connects_size(i-1)
        nodes_connects_size(i-1) = 0
      enddo
      !write(*,'(a,1024i4)')'nodes_connect_ptr',nodes_connects_ptr(1:num_nodes)
      !write(*,'(a,1024i4)')'nodes_connect_size',nodes_connects_size(1:num_nodes)
      
      !create the adjacency table
      nodes_adjacency_size = nodes_connects_ptr(num_nodes) + 
     1                                nodes_connects_size(num_nodes)-1
      nodes_connects_size(num_nodes) = 0
      !print*,'node_adjacency_size',nodes_adjacency_size
      allocate (nodes_adjacency(nodes_adjacency_size))
      nodes_adjacency(1:nodes_adjacency_size)=0
      
      !parse the edges to update the adjacency
      do ii = 1,ntree
        i=tree(ii)
        n1 = edges_n1(i)
        n2 = edges_n2(i)
        nodes_adjacency(nodes_connects_ptr(n1) + 
     1                        nodes_connects_size(n1)) = i
        nodes_adjacency(nodes_connects_ptr(n2) + 
     1                        nodes_connects_size(n2)) = i
        nodes_connects_size(n1) = nodes_connects_size(n1)+1
        nodes_connects_size(n2) = nodes_connects_size(n2)+1
      enddo
      !write(*,'(a,1024i4)')'nodes_connect_ptr',nodes_connects_ptr(1:num_nodes)
      !write(*,'(a,1024i4)')'nodes_connect_size',nodes_connects_size(1:num_nodes)
      
      !depth-first parse of the tree, starting from root
      !stack of node, parent
      allocate(stack(num_nodes,2))
      stack_size = 1
      stack(1,1) = root
      stack(1,2) = root
      
      ! sea has no sill
      sills(root,1) = -1
      sills(root,2) = -1
      
      !print*,'nodes_adjacency',nodes_adjacency
      !print*,'node_ptr',nodes_connects_ptr
      !print*,'nodes_size',nodes_connects_size
      
      if (continuous_flow) then
        sills(root, 1) = -1
        sills(root, 2) = -1
        basin_stack(1) = root
        basin_stack_size = 1
      endif
      
      do while (stack_size.gt.0)
        !get parsed node
        node = stack(stack_size, 1)
        parent = stack(stack_size, 2)
        stack_size = stack_size-1
        !  print*,'node',node
        !print*,'range',nodes_connects_ptr(node), nodes_connects_size(node)
      
        !for each edge of the graph
        !print*,'node',node
        do i=nodes_connects_ptr(node), nodes_connects_ptr(node) + 
     1                                    nodes_connects_size(node)-1
          !    print*,i
          edge_id = nodes_adjacency(i)
          !the edge comming from the parent node has already been updated.
        !in this case, the edge is (parent, node)
        if (edges_n1(edge_id).eq.parent.and.node.ne.parent) then
      
          if (continuous_flow) then
            basin_stack_size = basin_stack_size+1
            basin_stack(basin_stack_size) = node
      
            sills(node, 1) = edges_p1(edge_id)
            sills(node, 2) = edges_p2(edge_id)
          endif
      
        else
      
          !check if the node is not in n1
          if (node.ne.edges_n1(edge_id)) then
            !swap n1 and n2
            edge_swap = edges_n1(edge_id)
            edges_n1(edge_id) = edges_n2(edge_id)
            edges_n2(edge_id) = edge_swap
            !swap p1 and p2
            edge_swap = edges_p1(edge_id)
            edges_p1(edge_id) = edges_p2(edge_id)
            edges_p2(edge_id) = edge_swap
          endif
          !add the opposite node to the stack
          stack_size = stack_size+1
          stack(stack_size,1) =  edges_n2(edge_id)
          stack(stack_size,2) = node
      
        endif
      
      enddo
      enddo
      
      if (allocated(nodes_connects_size))deallocate(nodes_connects_size)
      if (allocated(nodes_connects_ptr)) deallocate(nodes_connects_ptr)
      if (allocated(nodes_adjacency)) deallocate(nodes_adjacency)
      if (allocated(stack)) deallocate(stack)
      return
      
      end subroutine order_tree
      
C----------------------
      subroutine correct_receivers(receivers,dist2receivers,outlets,
     1                        conn_basins,conn_nodes,tree,elevation,
     2                        nn,nbasins,nconn,ntree,length22)
      
      !Correct receivers: correct the receivers according the the tree order
      
      !Input:
      !receivers : array of receviers
      !dist2receivers : distance between bode and its receiver
      !outlets: outlets: minimal node of each basin
      !conn_basins: array of pairs of basins (b0, b1) that share a pass
      !conn_nodes: array of pairs (p0, p1) : nodes of the passes for each pair of basin
      !tree: id of the connections in the minimal spanning tree
      !elevation: elevation of each node
      !nx size of the demin x direction
      !dx, dy: size of a cell
      
      !Output: none: receivers and dist2receivers are updated in place
      
      implicit none
      
      integer nbasins,ntree,nconn,nn
      integer receivers(nn),outlets(nbasins),tree(ntree)
      integer conn_basins(nconn,2),conn_nodes(nconn,2)
      double precision dist2receivers(nn),elevation(nn)
      double precision, dimension(nn,2), intent(in):: length22
      
      integer i,ii,node_from,node_to,outlet_from
      
      do ii=1,ntree
        i=tree(ii)
        ! tree order: inverse of water flow
        node_to = conn_nodes(i, 1)
        node_from = conn_nodes(i, 2)
        if (node_from.eq.-1) goto 111
      
        outlet_from = outlets(conn_basins(i, 2))
      
        ! -> no river erosion on sinks
        dist2receivers(outlet_from) = 1.d10
      
        if (elevation(node_from).lt.elevation(node_to)) then
          receivers(outlet_from) = node_to
        else
          !receivers[outlet_from] = node_to
          !continue
          receivers(outlet_from) = node_from
          receivers(node_from) = node_to
          !distance based on previous king (4D) neighbor lookup
C          if (mod(node_from,nx).eq.mod(node_to,nx)) then
C            dist2receivers(node_from) = dx
C          else
C            dist2receivers(node_from) = dy
C          endif
          dist2receivers(node_from)=SQRT((length22(node_from,1)-
     1            length22(node_to,1))**2+(length22(node_from,2)-
     2            length22(node_to,2))**2)
          
        endif
      
111     continue
      enddo
      return
      
      end subroutine correct_receivers
      
C----------------------
      subroutine correct_receivers_v2(receivers,dist2receivers,outlets,
     1    conn_basins,conn_nodes,tree,nn,nbasins,nconn,ntree,length22)
      
      !Correct receivers: correct the receivers according the the tree order
      
      !Input:
      !receivers : array of receviers
      !dist2receivers : distance between bode and its receiver
      !outlets: outlets: minimal node of each basin
      !conn_basins: array of pairs of basins (b0, b1) that share a pass
      !conn_nodes: array of pairs (p0, p1) : nodes of the passes for each pair of basin
      !tree: id of the connections in the minimal spanning tree
      !nx size of the demin x direction
      !dx, dy: size of a cell
      
      !Output: none: receivers and dist2receivers are updated in place
      
      implicit none
      
      integer nn,nbasins,ntree,nconn
      double precision, dimension(nn,2), intent(in):: length22
      integer receivers(nn),outlets(nbasins),tree(ntree)
      integer conn_basins(nconn,2),conn_nodes(nconn,2)
      integer next_node,cur_node,rcv_next_node
      double precision dist2receivers(nn)
      double precision ddx,ddy,previous_dist,tmp
      
      integer i,ii,node_from,node_to,outlet_from
      
      do ii=1,ntree
        i=tree(ii)
        ! tree order: inverse of water flow
        node_to=conn_nodes(i,1)
        node_from=conn_nodes(i,2)
        if (node_from.eq.-1) goto 111
      
        outlet_from=outlets(conn_basins(i,2))
      
        next_node=receivers(node_from)
        cur_node=node_from
        previous_dist=dist2receivers(node_from)
      
        if ((next_node.eq.593).or.(node_from.eq.593)) then
        continue
        end if
        receivers(node_from)=node_to
        
        ddx=length22(node_from,1)-length22(node_to,1)
        ddy=length22(node_from,2)-length22(node_to,2)
      
        dist2receivers(node_from)=sqrt(ddx*ddx+ddy*ddy)
      
        do while (cur_node.ne.outlet_from)
          rcv_next_node=receivers(next_node)
          receivers(next_node)=cur_node
          tmp=previous_dist
          previous_dist=dist2receivers(next_node)
          dist2receivers(next_node)=tmp
          cur_node=next_node
          next_node=rcv_next_node
        end do
      
111   continue
      enddo
      
      return
      end subroutine correct_receivers_v2
      
      
C----------------------
      subroutine update_fake_topography (sills,basin_stack,basins,
     1                elevation,nn,nn1,nbasins,bond_num,length22)
      !Input:         
      
      ! sills: pairs of possible sills for each basin
      ! basin_stack: order of basin parsing from sea to top
      ! basins: basin id for each node
      ! nx, ny dx, dy: size and scale of dem
      
      implicit none
      
      integer nn,nn1,nbasins,i_sill, parse_begin, parse_end, i_node,bi
      integer sills(nbasins, 2), basin_stack(nbasins)
      integer basins(nn)
      double precision elevation(nn)
      integer i_b, b, nb_dir_i, i_nb, j_nb, i_neighbor
      double precision slope, dist_x, dist_y, new_height
      integer, dimension(:), allocatable :: parse, parsed
      integer, dimension(nn,nn1), intent(in) :: bond_num
      double precision, dimension(nn,2), intent(in):: length22
      
      slope = 1.d-10	! 假设内盆地流向其他盆地时导致该盆地高程变高
      
      ! parse queue
      ! sills can be parsed twice, so allocate nn*2 for safety
      allocate (parse(2*nn))
      
      !state of the nodes: "basins" means belong to basin and is not parsed ; -1 means parsed
      allocate (parsed(nn))
      parsed(1:nn) = basins(1:nn)
      
      !a breadth first parse will be performed, a queue structure is then necessary. "end" means last+1
      parse_begin = 1
      parse_end = 1
      
      !parse basins in flow order
      do i_b = 1, nbasins
        b = basin_stack(i_b)
      
        !boudary basins have sills -1.
        if (sills(b, 1) .eq. -1) goto 111
      
        !find sill
        if (elevation(sills(b, 1)) .ge. elevation(sills(b, 2))) then
          i_sill = sills(b, 1)
        else
          i_sill = sills(b, 2)
        endif
      
        ! start parsing at sill
        parse(parse_begin) = i_sill
        parse_end = parse_end+1
      
        parsed(i_sill) = -1
      
        do while (parse_begin.ne.parse_end)
      
          i_node = parse(parse_begin) ! 当前节点
          parse_begin = parse_begin+1
      
          ! parse neighbors 识别节点的8或者10个方向关联节点（abaqus修改）
          bi=bond_num(i_node,2)
          do nb_dir_i=1,bi	
            i_neighbor = bond_num(i_node,2+nb_dir_i)	! 当前循环的关联节点号（abaqus修改）
            if (parsed(i_neighbor) .eq. b) then
              parsed(i_neighbor) = -1
              dist_x = length22(i_node,1)-length22(i_neighbor,1)
              dist_y = length22(i_node,2)-length22(i_neighbor,2)
              new_height = elevation(i_node) + slope * 
     1                          sqrt(dist_x*dist_x + dist_y*dist_y)	! 假设内盆地的水流向开放盆地溢出导致周边节点高程变化
              if (elevation(i_neighbor) .le. new_height) then
                elevation(i_neighbor) = new_height
                parse(parse_end) = i_neighbor
                parse_end = parse_end+1
              endif
            endif
          enddo
        enddo
      
111     continue
      enddo
      
      if (allocated(parse)) deallocate(parse)
      if (allocated(parsed)) deallocate(parsed)
      return
        
      end subroutine update_fake_topography
      
C----------------------
      subroutine UnionFindInit (parent,rank,n)
      implicit none
      integer n,i
      integer parent(n),rank(n)
      
      do i=1,n
        parent(i)=i
        rank(i)=0
      enddo
      return
      
      end subroutine UnionFindInit
      
C----------------------
      subroutine DoUnion (x,y,parent,rank,n)
      implicit none
      integer n,x,y
      integer parent(n),rank(n)
      integer xroot,yroot
      
      call UnionFind(x,parent,n,xroot)
      call UnionFind(y,parent,n,yroot)
      
      if (xroot.ne.yroot) then
        if (rank(xroot).lt.rank(yroot)) then
          parent(xroot)=yroot
        else
          parent(yroot)=xroot
          if (rank(xroot).eq.rank(yroot)) then
            rank(xroot)=rank(xroot)+1
          endif
        endif
      endif
      return
      
      end subroutine DoUnion
      
C----------------------
      subroutine UnionFind (x,parent,n,found)
      implicit none
      integer x,xp,xc,n,found
      integer parent(n)
      
      xp=x
      xc=-1000
      do while (xc.ne.xp)
        xc=xp
        xp=parent(xc)
      enddo
      parent(x)=xc
      found=xc
      return
      
      end subroutine UnionFind
      
C----------------------
      subroutine loc_min_3_find_receivers (h,rec,length,bc,nn,nn1,
     1                                        bond_num,length22)
      implicit none
      integer ij,ii,ijk,bi,nn,nn1
      double precision h(nn),length(nn)
      integer rec(nn)
      logical bc(nn)
      double precision smax,l,slope
      integer, dimension(nn,nn1), intent(in) :: bond_num
      double precision, dimension(nn,2), intent(in):: length22
      
      do ij=1,nn
         smax=tiny(smax)
         rec(ij)=ij
         if (bc(ij)) goto 111
         bi=bond_num(ij,2)
         do ii=1,bi
            ijk=bond_num(ij,2+ii)
            if (ijk.ne.ij) then
              l=sqrt((length22(ij,1)-length22(ijk,1))**2+
     1                    (length22(ij,2)-length22(ijk,2))**2)
              slope=(h(ij)-h(ijk))/l
              if (slope.gt.smax) then
                smax=slope
                rec(ij)=ijk
                length(ij)=l
              endif
            endif
         enddo
111      continue
      enddo
      return
      
      end subroutine loc_min_3_find_receivers
      
C----------------------
      subroutine loc_min_3_find_donors (rec,ndon,donor,nn)
      implicit none
      integer nn
      integer rec(nn),ndon(nn),donor(10,nn)
      integer ij,ijk
      
      ndon=0
      do ij=1,nn
        if (rec(ij).ne.ij) then
          ijk=rec(ij)
          ndon(ijk)=ndon(ijk)+1
          donor(ndon(ijk),ijk)=ij
        endif
      enddo
      return
      
      end subroutine loc_min_3_find_donors
      
C----------------------
      subroutine loc_min_3_find_stack (rec,ndon,donor,stack,nn)
      implicit none
      integer nn
      integer rec(nn),stack(nn),ndon(nn),donor(10,nn)
      integer nstack,ij
      
      nstack=0
      do ij=1,nn
        if (rec(ij).eq.ij) then
          nstack=nstack+1
          stack(nstack)=ij
          call find_stack_recursively_locmin (ij,donor,ndon,nn,
     1                                        stack,nstack)
        endif
      enddo
      return
      
      end subroutine loc_min_3_find_stack
      
C-----------------------
      subroutine loc_min_3_indexx(n,arr,indx) ! 快速排序和插入排序混合算法
      implicit none
      integer n,M,NSTACK
      integer indx(n)
      double precision arr(n)
      parameter (M=7,NSTACK=50)	! M是插入排序的区间大小，NSTACK是快速排序区间大小
      integer i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      double precision a
      
      do j=1,n
        indx(j)=j
      enddo
      
      jstack=0
      l=1
      ir=n
      
1     if(ir-l.lt.M)then
        do j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
          enddo
          i=0
2         indx(i+1)=indxt
        enddo
        if(jstack.eq.0) return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
        i=i+1
        if(arr(indx(i)).lt.a) goto 3
4       continue
        j=j-1
        if(arr(indx(j)).gt.a) goto 4
        if(j.lt.i) goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
      
        if(jstack.gt.NSTACK) then
        write(6,*) 'NSTACK too small in loc_min_3_indexx'
        end if
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      
      end subroutine loc_min_3_indexx
      
C----------------------
      recursive subroutine find_stack_recursively_locmin (ij,donor,ndon,
     1                                            nn,stack,nstack)
      implicit none
      integer ij,nn,nstack
      integer donor(10,nn),ndon(nn),stack(nn)
      integer k,ijk
      
      do k=1,ndon(ij)
      ijk=donor(k,ij)
      nstack=nstack+1
      stack(nstack)=ijk
      call find_stack_recursively_locmin (ijk,donor,ndon,nn,
     1                                        stack,nstack)
      enddo
      return
      
      end subroutine find_stack_recursively_locmin

C-------------------------------------------------------------------------------------------
      subroutine StreamPowerLaw ()
        ! subroutine to solve the stream power law equation following the FastScape method described
        ! in Braun and Willett, Geomorphology, 2015
        use FastScapeContext
        implicit none
        integer :: ij,ijk,ijr,k,ijr1
        double precision :: fact,tol,err
        double precision :: f,df,errp,h0,hn,omega,tolp,w_rcv
        double precision, dimension(:), allocatable :: htt,kfint,dhh,hp
        double precision, dimension(:), allocatable :: elev
        double precision, dimension(:), allocatable :: water,
     1                                lake_water_volume,lake_sediment
        integer, dimension(:), allocatable :: lake_sill
        allocate (htt(nn),kfint(nn),dhh(nn),hp(nn))
        allocate (elev(nn))
        allocate (water(nn),lake_water_volume(nn),lake_sediment(nn),
     1                lake_sill(nn))
      
        ! set g, dimensionless parameter for sediment transport and deposition
        ! if g1<0, skip and use g values directly from FastScapeContext (not in API!!!)
        if (g1.ge.0.d0) then
          gob(1:nn)=g1	! 默认为基岩的沉积/运输系数（g）
          if (g2.gt.0.d0) then    ! 如果该点的厚度超过基岩厚度，则使用沉积物g值
          where ((h-b).gt.1e-5) 
              gob(1:nn)=g2
          end where
          end if
        endif
      
        ! set kf / kfsed
        kfint(1:nn)=kf(1:nn)	! 默认为基岩的河流切割参数（kf）
        if (kfsed.gt.0.d0) then   ! 如果该点的厚度超过基岩厚度，则使用沉积物kf值
        where ((h-b).gt.0.01) ! 大于1米才认为是沉积物
        kfint(1:nn)=kfsed
c        write(6,*) 'Sediment was found1'
        end where
        end if
        
        if (count(mstack==0).ne.0) then
        write(6,*) 'incomplete stack',count(mstack==0),nn
        end if
        ! modified by Jean Braun (20/11/2022) to allow for relative versus
        ! absolute tolerance 计算迭代时的容差
        tol = tol_rel*maxval(abs(h)) + tol_abs
        err=2.d0*tol
      
        ! store the elevation at t
        htt(1:nn)=h(1:nn)
      
        ! Gauss-Seidel iteration 迭代步计数器
        nGSStreamPowerLaw=0
        lake_sediment=0.d0
        lake_sill=0.d0
        dhh=0.d0
        hp(1:nn)=h(1:nn)
        
        ! htt是t时刻的高程，hp是迭代更新后（t+dt）的高程
        do while (err.gt.tol.and.nGSStreamPowerLaw.lt.
     1                nGSStreamPowerLawMax-1)
          nGSStreamPowerLaw=nGSStreamPowerLaw+1
      
          where (bounds_bc)	! true是固定约束（部分边界），fault是可变约束
            elev(1:nn)=htt(1:nn)
          elsewhere
            elev(1:nn)=htt(1:nn)+(dhh(1:nn)-(htt(1:nn)-hp(1:nn)))*
     1                    gob(1:nn)*sear(1:nn)/a(1:nn)	! 高程+沉积量=htt+G/pA*(∫(U-dhh)dA)
          endwhere
      
            ! apply modified stream power law using lake surface (hwater)
      	
            if (abs(n-1.d0).lt.tiny(n)) then	! n系数等于1不迭代
      
              do ij=nn,1,-1
                ijk=mstack(ij)
                ijr1=rec(ijk)
                if (ijr1.eq.ijk) then
                  water(ijk)=htt(ijk)
                  lake_sill(ijk)=ijk
                  lake_water_volume(ijk)=0.d0
                else
                  w_rcv=water(ijr1)
                  if (elev(ijk).gt.w_rcv) then
                    if (mnrec(ijk).gt.0) then
                      if (h(ijk).ge.sealevel.or..not.runMarine) then
                        f = elev(ijk)
                        df = 1.d0
                        do k=1,mnrec(ijk)
                          if (htt(ijk).ge.htt(mrec(k,ijk))) then
                            fact = kfint(ijk)*dt*(a(ijk)*
     1                                mwrec(k,ijk))**m/mlrec(k,ijk)
                            f = f + fact*h(mrec(k,ijk))
                            df = df + fact
                          endif
                        enddo
                        h(ijk)=f/df
                      endif
                    endif
                    lake_sill(ijk)=ijk
                    lake_water_volume(ijk)=0.d0
                    if (h(ijk).lt.w_rcv) h(ijk)=w_rcv
                  else
                    h(ijk)=elev(ijk)
                    lake_sill(ijk)=lake_sill(ijr1)
                    if (lake_sill(ijk).ne.0) 
     1                lake_water_volume(lake_sill(ijk)) = 
     2                lake_water_volume(lake_sill(ijk))+(w_rcv-h(ijk))
                  endif
                  water(ijk)=max(w_rcv,h(ijk))
                endif
              enddo
      
            else
      
              do ij=nn,1,-1	! 倒叙循环
                ijk=mstack(ij)	! stack堆栈是从顶点开始记录
                ijr1=rec(ijk)
                if (ijr1.eq.ijk) then	! 判断节点是否是开放盆地的最低点
                  water(ijk)=htt(ijk)	! 记录湖水深
                  lake_sill(ijk)=ijk	! 记录湖水的节点编号
                  lake_water_volume(ijk)=0.d0	! 记录湖水量为0（开放盆地，即入海口）
                else
                  w_rcv=water(ijr1)
                  if (elev(ijk).gt.w_rcv) then	! 上节点是否高于下节点
                    if (mnrec(ijk).gt.0) then	! mnrec是多流向节点的接收者数
                      if (htt(ijk).ge.sealevel.or..not.runMarine) then	! sealevel海平面高度，runMarine默认是true
                        omega=0.875d0/n	! SOR松弛因子
                        tolp=1.d-3		! 迭代收敛阈值
                        errp=2.d0*tolp	! 初始化误差
                        h0=elev(ijk)		! 初始节点的高程
                        do while (errp.gt.tolp)	! Newton-Raphson-SOR迭代法
                          f=h(ijk)-h0		! 节点的高程增量
                          df=1.d0			! 初始导数1
                          do k=1,mnrec(ijk)	! 多流向分支计算剥蚀
                            if (htt(ijk).gt.htt(mrec(k,ijk))) then	! 上节点是否高于下节点
                              fact = kfint(ijk)*dt*(a(ijk)*
     1                                mwrec(k,ijk))**m/mlrec(k,ijk)**n	! 剥蚀导致的高程变化=Kf*(P^m)*(A^m)*(S^n)
                              f=f+fact*max(0.d0,h(ijk)-
     1                                        h(mrec(k,ijk)))**n							! S^n=((h(ijk)-h(mrec(k,ijk)))/mlrec)**n	剥蚀量
                              df=df+fact*n*max(0.d0,h(ijk)-
     1                                    h(mrec(k,ijk)))**(n-1.d0) ! (S^n)'=n*((h(ijk)-h(mrec(k,ijk)))/mlrec)**(n-1)	剥蚀量的导数
                            endif
                          enddo
                          hn=h(ijk)-f/df	! 牛顿迭代法求得的添加剥蚀量后的高程xk+1=xk-f(xk)/f'(xk)
                          errp=abs(hn-h(ijk))	! 迭代后的高程（hn）与原高程（h）的差值
                          h(ijk)=h(ijk)*(1.d0-omega)+hn*omega	! 使用松弛因子moega更新h，加速收敛
                        enddo
                      endif
                    endif
                    lake_sill(ijk)=ijk
                    lake_water_volume(ijk)=0.d0	! 因为上节点一定高于下节点，所有上节点无蓄水
                    if (h(ijk).lt.w_rcv) h(ijk)=w_rcv	! 更新后的节点低于其下节点高度（当前最低剥蚀高度），则更新该节点与下节点高度一致
                  else	! 若该节点为内盆地最低点（但是该盆地有比高的节点接收者，即该盆地连接了其他开放盆地）
                    h(ijk)=elev(ijk)	! 未被剥蚀，高程不变
                    lake_sill(ijk)=lake_sill(ijr1)	! 同属于其连接的湖
                    if (lake_sill(ijk).ne.0) 
     1                lake_water_volume(lake_sill(ijk))=
     2                lake_water_volume(lake_sill(ijk))+(w_rcv-h(ijk))	! 默认该盆地已填满水（湖量为该点到上点的距离）
                  endif
                  water(ijk)=max(w_rcv,h(ijk))	! 选择接收点和该点的最大高程为水深（若该点为内盆地最低点且向外流动，则默认该盆地已蓄满水）
                endif
              enddo
      
            endif
      
            err=sqrt(sum((h-hp)**2)/nn)	! 更新迭代误差
          ! Jean Braun modification 18/11/2022: moved the computation of redistribution of sediment in lakes
          ! following Sebastian Wolf's suggestion; this ensures mass conservation in multi-minima cases
          ! guess/update the elevation at t+Dt (k)
          hp(1:nn)=h(1:nn)
      
          ! calculate erosion/deposition at each node
          dhh(1:nn)=htt(1:nn)-hp(1:nn)
      
          ! sum the erosion in stack order
          do ij=1,nn
            ijk=mstack(ij)
            ijr1=rec(ijk)
            if (ijr1.ne.ijk) then
              dhh(ijk)=dhh(ijk)-(htt(ijk)-hp(ijk))	! 来自上方点的沉积量=该点总沉积量-该点剥蚀量
              if (lake_sill(ijk).eq.ijk) then		! 判断是否属于正常流动的盆地点（有比起低的接受点）
                if (dhh(ijk).le.0.d0) then			! 判断有无上方点的沉积量
                  lake_sediment(ijk)=0.d0			! 该点为顶点无沉积量
                else
                  lake_sediment(ijk)=dhh(ijk)		! 该点的沉积量=来自上方点的沉积量
                endif
              endif
              dhh(ijk)=dhh(ijk)+(htt(ijk)-hp(ijk))	! 恢复该点总沉积量
              do k=1,mnrec(ijk)					! 计算该点对下方点的沉积量贡献
                ijr=mrec(k,ijk)
                dhh(ijr)=dhh(ijr)+dhh(ijk)*mwrec(k,ijk)	! 将上方点的剥蚀量作为下方点的沉积量进行多流向的分配
              enddo
            else
              lake_sediment(ijk)=dhh(ijk)	! 盆地最低点的沉积量=该点剥蚀量+上方点的剥蚀量分配
            endif
          enddo
          
          if (maxval(gob(1:nn)).lt.tiny(gob(1:nn))) err=0.d0	! 若沉积参数为0，则默认收敛；因为剥蚀量已经迭代过
      enddo
      
          b(1:nn)=min(h(1:nn),b(1:nn))	! 设置最低的高程为基岩
          do ij=1,nn
            if (lake_sill(ij).ne.0) then
      		! 判断内盆地点（即该点有蓄水）
              if (lake_water_volume(lake_sill(ij)).gt.0.d0) 
     1            h(ij)=h(ij)+max(0.d0,min(lake_sediment(lake_sill(ij)),
     2            lake_water_volume(lake_sill(ij))))/
     3            lake_water_volume(lake_sill(ij))*(water(ij)-h(ij))	! 内盆地点高度=原高程+(上方沉积量/蓄水量)*蓄水量  ------若上方沉积量大于蓄水量，则认为该点被填平与蓄水量一致
            endif
          enddo
      
      	! 更新剥蚀量/剥蚀速率/流量
          ! stores total erosion, erosion rate and flux for output
          etot(1:nn)=etot(1:nn)+htt(1:nn)-h(1:nn)
          erate(1:nn)=(htt(1:nn)-h(1:nn))/dt
          Sedflux(1:nn)=htt(1:nn)-h(1:nn)
          !if (runMarine) where (h.lt.sealevel) Sedflux=0.d0
      
cc          deallocate (htt,kfint,dhh,hp,elev,water,lake_water_volume,
cc     1                lake_sediment,lake_sill)
          if (allocated(htt)) deallocate(htt)
          if (allocated(kfint)) deallocate(kfint)
          if (allocated(dhh)) deallocate(dhh)
          if (allocated(hp)) deallocate(hp)
          if (allocated(elev)) deallocate(elev)
          if (allocated(water)) deallocate(water)
          if (allocated(lake_water_volume))deallocate(lake_water_volume)
          if (allocated(lake_sediment)) deallocate(lake_sediment)
          if (allocated(lake_sill)) deallocate(lake_sill)
          
          return
      
        end subroutine StreamPowerLaw

C--------------------------------------------------------------------------------------------
      subroutine StreamPowerLawSingleFlowDirection ()
        ! subroutine to solve the stream power law equation following the FastScape method described
        ! in Braun and Willett, Geomorphology, 2015
        use FastScapeContext
        implicit none
        integer :: ij,ijk,ijr
        double precision :: fact,tol,err
        double precision :: f,df,errp,h0,hn,omega,tolp,w_rcv
        double precision, dimension(:), allocatable :: htt,kfint,dhh,hp
        double precision, dimension(:), allocatable :: elev
        double precision, dimension(:), allocatable :: water,
     1                            lake_water_volume,lake_sediment
        integer, dimension(:), allocatable :: lake_sill
      
        allocate (htt(nn),kfint(nn),dhh(nn),hp(nn))
        allocate (elev(nn))
        allocate (water(nn),lake_water_volume(nn),lake_sediment(nn),
     1            lake_sill(nn))
        
        ! set g, dimensionless parameter for sediment transport and deposition
        ! if g1<0, skip and use g values directly from FastScapeContext (not in API!!!)
        if (g1.ge.0.d0) then
          gob(1:nn)=g1
          if (g2.gt.0.d0) then
          where ((h-b).gt.1e-5) 
              gob(1:nn)=g2
          end where
          end if
        endif
      
        ! set kf / kfsed
        kfint(1:nn)=kf(1:nn)
        if (kfsed.gt.0.d0) then
        where ((h-b).gt.0.01)
            kfint(1:nn)=kfsed
c            write(6,*) 'Sediment was found2',kfsed
        end where
        end if
        
        ! modified by Jean Braun (20/11/2022) to allow for relative versus
        ! absolute tolerance
        tol = tol_rel*maxval(abs(h)) + tol_abs 
        err=2.d0*tol
      
        ! store the elevation at t
        htt(1:nn)=h(1:nn)
      
        ! Gauss-Seidel iteration
        nGSStreamPowerLaw=0
      
        lake_sediment=0.d0
        dhh=0.d0
        hp(1:nn)=h(1:nn)
      
        do while (err.gt.tol.and.nGSStreamPowerLaw.lt.
     1                nGSStreamPowerLawMax-1)
          nGSStreamPowerLaw=nGSStreamPowerLaw+1
          where (bounds_bc)
            elev(1:nn)=htt(1:nn)
          elsewhere
            elev(1:nn)=htt(1:nn)+(dhh(1:nn)-(htt(1:nn)-hp(1:nn)))*
     1                gob(1:nn)*sear(1:nn)/a(1:nn)
          endwhere
      
            ! apply modified stream power law using lake surface (hwater)
            if (abs(n-1.d0).lt.tiny(n)) then
              do ij=1,nn
                ijk=stack(ij)
                ijr=rec(ijk)
                if (ijr.eq.ijk) then
                  water(ijk)=htt(ijk)
                  lake_sill(ijk)=ijk
                  lake_water_volume(ijk)=0.d0
                else
                  w_rcv=water(ijr)
                  if (elev(ijk).gt.w_rcv) then
                    if (h(ijk).ge.sealevel.or..not.runMarine) then
                    f = elev(ijk)
                    df = 1.d0
                    ! todo: check if we don't need those checks for single flow
C                    if (htt(ijk).ge.htt(ijr)) then
                      fact = kfint(ijk)*dt*a(ijk)**m/length(ijk)
                      f = f + fact*h(ijr)
                      df = df + fact
C                    endif
                    h(ijk)=f/df
C                    h(ijk)=min(f/df,minval(h(don(1:ndon(ijk),ijk))))
                      endif
                    lake_sill(ijk)=ijk
                    lake_water_volume(ijk)=0.d0
                    if (h(ijk).lt.w_rcv) h(ijk)=w_rcv
                  else
                    h(ijk)=elev(ijk)
                    lake_sill(ijk)=lake_sill(ijr)
                    if (lake_sill(ijk).ne.0) 
     1                lake_water_volume(lake_sill(ijk))=
     2                lake_water_volume(lake_sill(ijk))+(w_rcv-h(ijk))
                  endif
                  water(ijk)=max(w_rcv,h(ijk))
                endif
              enddo
      
            else
      
              do ij=1,nn
                ijk=stack(ij)
                ijr=rec(ijk)
                if (ijr.eq.ijk) then
                  water(ijk)=htt(ijk)
                  lake_sill(ijk)=ijk
                  lake_water_volume(ijk)=0.d0
                else
                  w_rcv=water(ijr)
                  if (elev(ijk).gt.w_rcv) then
                    if (htt(ijk).ge.sealevel.or..not.runMarine) then	
                    omega=0.875d0/n
                    tolp=1.d-3
                    errp=2.d0*tolp
                    h0=elev(ijk)
                    do while (errp.gt.tolp)
                      f=h(ijk)-h0
                      df=1.d0
                      if (htt(ijk).gt.htt(ijr)) then
                        fact = kfint(ijk)*dt*a(ijk)**m/length(ijk)**n
                        f=f+fact*max(0.d0,h(ijk)-h(ijr))**n
                        df=df+fact*n*max(0.d0,h(ijk)-h(ijr))**(n-1.d0)
                      endif
                      hn=h(ijk)-f/df
                      errp=abs(hn-h(ijk))
                      h(ijk)=h(ijk)*(1.d0-omega)+hn*omega
                    enddo
                    endif
                    lake_sill(ijk)=ijk
                    lake_water_volume(ijk)=0.d0
                    if (h(ijk).lt.w_rcv) h(ijk)=w_rcv
                  else
                    h(ijk)=elev(ijk)
                    lake_sill(ijk)=lake_sill(ijr)
                    if (lake_sill(ijk).ne.0) 
     1                lake_water_volume(lake_sill(ijk)) = 
     2                lake_water_volume(lake_sill(ijk))+(w_rcv-h(ijk))
                  endif
                  water(ijk)=max(w_rcv,h(ijk))
                endif
              enddo
      
            endif
      
            err=sqrt(sum((h(1:nn)-hp(1:nn))**2)/nn)
      
          ! Jean Braun modification 18/11/2022: moved the computation of redistribution of sediment in lakes
          ! following Sebastian Wolf's suggestion; this ensures mass conservation in multi-minima cases
          ! guess/update the elevation at t+Dt (k)
          hp(1:nn)=h(1:nn)
      
           ! calculate erosion/deposition at each node
          dhh(1:nn)=htt(1:nn)-hp(1:nn)
      
          ! sum the erosion in stack order
          do ij=nn,1,-1
            ijk=stack(ij)
            ijr=rec(ijk)
            if (ijr.ne.ijk) then
              dhh(ijk)=dhh(ijk)-(htt(ijk)-hp(ijk))
              if (lake_sill(ijk).eq.ijk) then
                if (dhh(ijk).le.0.d0) then
                  lake_sediment(ijk)=0.d0
                else
              lake_sediment(ijk)=min(dhh(ijk),lake_water_volume(ijk))
                  dhh(ijk) = dhh(ijk)-lake_water_volume(ijk) !remove the sediment that is going to be deposited in the lake from the dhh stack
                  if (dhh(ijk)<0.d0) dhh(ijk)=0.d0
                endif
              endif
              dhh(ijk)=dhh(ijk)+(htt(ijk)-hp(ijk))
              dhh(ijr)=dhh(ijr)+dhh(ijk)
            else
              lake_sediment(ijk)=dhh(ijk)
            endif
          enddo
      
          if (maxval(gob(1:nn)).lt.tiny(gob(1:nn))) err=0.d0
          enddo
      
          b(1:nn)=min(h(1:nn),b(1:nn))
          do ij=1,nn
            if (lake_sill(ij).ne.0) then
              if (lake_water_volume(lake_sill(ij)).gt.0.d0) h(ij)=h(ij)
     1                +max(0.d0,min(lake_sediment(lake_sill(ij)),
     2                lake_water_volume(lake_sill(ij))))/
     3                lake_water_volume(lake_sill(ij))*(water(ij)-h(ij))
            endif
          enddo
      
          ! stores total erosion, erosion rate and flux for output
          etot(1:nn)=etot(1:nn)+htt(1:nn)-h(1:nn)
          erate(1:nn)=(htt(1:nn)-h(1:nn))/dt
          Sedflux(1:nn)=htt(1:nn)-h(1:nn)
          !if (runMarine) where (h.lt.sealevel) Sedflux=0.d0
      
cc          deallocate (htt,kfint,dhh,hp,elev,water,lake_water_volume,
cc     1                lake_sediment,lake_sill)
        if (allocated(htt)) deallocate(htt)
        if (allocated(kfint)) deallocate(kfint)
        if (allocated(dhh)) deallocate(dhh)
        if (allocated(hp)) deallocate(hp)
        if (allocated(elev)) deallocate(elev)
        if (allocated(water)) deallocate(water)
        if (allocated(lake_water_volume)) deallocate(lake_water_volume)
        if (allocated(lake_sediment)) deallocate(lake_sediment)
        if (allocated(lake_sill)) deallocate(lake_sill)
      
          return

        end subroutine StreamPowerLawSingleFlowDirection

C--------------------------------------------------------------------------------------------
      subroutine Diffusion ()
      ! subroutine to solve the diffusion equation by ADI
      ! 嵌入的变量名要重新赋予在abaqus子程序内
      use FastScapeContext
      implicit none
      double precision, dimension(:), allocatable :: f,diag,sup,
     1                                    zintp,zint,kdint,inf,res
                                                    
      integer i,j,ij,bi
      double precision factxp,factxm,factyp,factym,dx1,dx2,dy1,dy2
      character cbc*4
      
      !print*,'Diffusion'
      write (cbc,'(i4)') bounds_ibc
      
      ! creates 2D internal arrays to store topo and kd
      
      allocate (zint(nn),kdint(nn),zintp(nn))
      do ij=1,nn
      zint(ij)=h(ij)	! 高程
      kdint(ij)=kd(ij)	! 基岩的扩散系数
      if (kdsed.gt.0.d0.and.(h(ij)-b(ij)).gt.1.d-6) kdint(ij)=kdsed   ! 沉积物上的扩散系数（kdsed < 0时全部使用基岩扩散系数）
      enddo
      
      zintp(1:nn) = zint(1:nn)
      
      ! 交替方向隐式法
      ! first pass along the x-axis
      ! 第一步先对x方向隐式
      allocate (f(nn),diag(nn),sup(nn),inf(nn),res(nn))
      f=0.d0
      diag=0.d0
      sup=0.d0
      inf=0.d0
      res=0.d0
      ! 矩阵组装
      do ij=1,nn
      bi=bond_num(ij,2)
      if (bi.ge.6) then
      dx1=abs(length2(bond_num(ij,3),1)-length2(ij,1))
      dx2=abs(length2(bond_num(ij,4),1)-length2(ij,1))
      dy1=abs(length2(bond_num(ij,5),2)-length2(ij,2))
      dy2=abs(length2(bond_num(ij,6),2)-length2(ij,2))
      factxm=(kdint(bond_num(ij,3))+kdint(ij))/2.d0*(dt/2.)/dx1**2
      factxp=(kdint(bond_num(ij,4))+kdint(ij))/2.d0*(dt/2.)/dx2**2
      factyp=(kdint(bond_num(ij,5))+kdint(ij))/2.d0*(dt/2.)/dy1**2
      factym=(kdint(bond_num(ij,6))+kdint(ij))/2.d0*(dt/2.)/dy2**2
      diag(ij)=1.d0+factxp+factxm	! 中心对角系数
      sup(ij)=-factxp				! 上对角系数
      inf(ij)=-factxm				! 下对角系数
      f(ij)=zintp(ij)+factyp*zintp(bond_num(ij,5))
     1    +factym*zintp(bond_num(ij,6))-(factyp+factym)*zintp(ij)
      
cccc      elseif ((bi.eq.5).and.(bond_num(ij,1).eq.1)
cccc     1                        .and.(cbc(1:1).ne.'1')) then
cccc      dx1=abs(length2(bond_num(ij,3),1)-length2(ij,1))
cccc      dx2=abs(length2(bond_num(ij,4),1)-length2(ij,1))
cccc      dy1=abs(length2(bond_num(ij,5),2)-length2(ij,2))
cccc      factxm=(kdint(bond_num(ij,3))+kdint(ij))/2.d0*(dt/2.)/dx1**2
cccc      factxp=(kdint(bond_num(ij,4))+kdint(ij))/2.d0*(dt/2.)/dx2**2
cccc      factyp=(kdint(bond_num(ij,5))+kdint(ij))/2.d0*(dt/2.)/dy1**2
cccc      diag(ij)=1.d0+factxp+factxm	! 中心对角系数
cccc      sup(ij)=-factxp				! 上对角系数
cccc      inf(ij)=-factxm				! 下对角系数
cccc      f(ij)=zintp(ij)+factyp*zintp(bond_num(ij,5))-factyp*zintp(ij)
      
      elseif ((bi.eq.5).and.(bond_num(ij,1).eq.2)
     1                        .and.(cbc(2:2).ne.'1')) then
      dx2=abs(length2(bond_num(ij,3),1)-length2(ij,1))
      dy1=abs(length2(bond_num(ij,4),2)-length2(ij,2))
      dy2=abs(length2(bond_num(ij,5),2)-length2(ij,2))
      factxm=(kdint(bond_num(ij,3))+kdint(ij))/2.d0*(dt/2.)/dx2**2
      factyp=(kdint(bond_num(ij,4))+kdint(ij))/2.d0*(dt/2.)/dy1**2
      factym=(kdint(bond_num(ij,5))+kdint(ij))/2.d0*(dt/2.)/dy2**2
      diag(ij)=1.d0+factxm	! 中心对角系数
      sup(ij)=0			! 上对角系数
      inf(ij)=-factxm				! 下对角系数
      f(ij)=zintp(ij)+factyp*zintp(bond_num(ij,4))
     1    +factym*zintp(bond_num(ij,5))-(factyp+factym)*zintp(ij)
      
cccc      elseif ((bi.eq.5).and.(bond_num(ij,1).eq.3)
cccc     1                        .and.(cbc(3:3).ne.'1')) then
cccc      dx1=abs(length2(bond_num(ij,3),1)-length2(ij,1))
cccc      dx2=abs(length2(bond_num(ij,4),1)-length2(ij,1))
cccc      dy2=abs(length2(bond_num(ij,5),2)-length2(ij,2))
cccc      factxm=(kdint(bond_num(ij,3))+kdint(ij))/2.d0*(dt/2.)/dx1**2
cccc      factxp=(kdint(bond_num(ij,4))+kdint(ij))/2.d0*(dt/2.)/dx2**2
cccc      factym=(kdint(bond_num(ij,5))+kdint(ij))/2.d0*(dt/2.)/dy2**2
cccc      diag(ij)=1.d0+factxp+factxm	! 中心对角系数
cccc      sup(ij)=-factxp				! 上对角系数
cccc      inf(ij)=-factxm				! 下对角系数
cccc      f(ij)=zintp(ij)+factym*zintp(bond_num(ij,5))-factym*zintp(ij)
      
      elseif ((bi.eq.5).and.(bond_num(ij,1).eq.4)
     1                        .and.(cbc(4:4).ne.'1')) then
      dx1=abs(length2(bond_num(ij,3),1)-length2(ij,1))
      dy1=abs(length2(bond_num(ij,4),2)-length2(ij,2))
      dy2=abs(length2(bond_num(ij,5),2)-length2(ij,2))
      factxp=(kdint(bond_num(ij,3))+kdint(ij))/2.d0*(dt/2.)/dx1**2
      factyp=(kdint(bond_num(ij,4))+kdint(ij))/2.d0*(dt/2.)/dy1**2
      factym=(kdint(bond_num(ij,5))+kdint(ij))/2.d0*(dt/2.)/dy2**2
      diag(ij)=1.d0+factxp	! 中心对角系数
      sup(ij)=-factxp				! 上对角系数
      inf(ij)=0				! 下对角系数
      f(ij)=zintp(ij)+factyp*zintp(bond_num(ij,4))
     1    +factym*zintp(bond_num(ij,5))-(factyp+factym)*zintp(ij)
      
      else    ! 判断是固定边界不进行扩散计算
      diag(ij)=1.
      sup(ij)=0.
      inf(ij)=0.
      f(ij)=zintp(ij)
      end if
      end do
      call tridag (inf(1:nn),diag(1:nn),sup(1:nn),f(1:nn),res(1:nn),nn)
      do ij=1,nn
      zint(ij)=res(ij)
      enddo
      
C     second pass along y-axis
C     第二步对y方向隐式      
      f=0.d0
      diag=0.d0
      sup=0.d0
      inf=0.d0
      res=0.d0
      ! 矩阵组装
      do ij=1,nn
      bi=bond_num(ij,2)
      if (bi.eq.8) then   ! 中心边界
      dx1=abs(length2(bond_num(ij,3),1)-length2(ij,1))
      dx2=abs(length2(bond_num(ij,4),1)-length2(ij,1))
      dy1=abs(length2(bond_num(ij,5),2)-length2(ij,2))
      dy2=abs(length2(bond_num(ij,6),2)-length2(ij,2))
      factxm=(kdint(bond_num(ij,3))+kdint(ij))/2.d0*(dt/2.)/dx1**2
      factxp=(kdint(bond_num(ij,4))+kdint(ij))/2.d0*(dt/2.)/dx2**2
      factyp=(kdint(bond_num(ij,5))+kdint(ij))/2.d0*(dt/2.)/dy1**2
      factym=(kdint(bond_num(ij,6))+kdint(ij))/2.d0*(dt/2.)/dy2**2
      diag(ij)=1.d0+factyp+factym	! 中心对角系数
      sup(ij)=-factyp				! 上对角系数
      inf(ij)=-factym				! 下对角系数
      f(ij)=zint(ij)+factxp*zint(bond_num(ij,4))
     1    +factxm*zint(bond_num(ij,3))-(factxp+factxm)*zint(ij)
c      write(6,*)'ij,bi,dx1,dx2,dy1,dy2',ij,bi,dx1,dx2,dy1,dy2
c      write(6,*)'bond_num(ij,3:6)',bond_num(ij,3:6)
c      write(6,*)'factxp,factxm,fatyp,factym',factxp,factxm,factyp,factym
c      write(6,*)'diag,sup,inf,f',diag(ij),sup(ij),inf(ij),f(ij)
      
      elseif ((bi.eq.5).and.(bond_num(ij,1).eq.1)     ! 下边界
     1                        .and.(cbc(1:1).ne.'1')) then
      dx1=abs(length2(bond_num(ij,3),1)-length2(ij,1))
      dx2=abs(length2(bond_num(ij,4),1)-length2(ij,1))
      dy1=abs(length2(bond_num(ij,5),2)-length2(ij,2))
      factxm=(kdint(bond_num(ij,3))+kdint(ij))/2.d0*(dt/2.)/dx1**2
      factxp=(kdint(bond_num(ij,4))+kdint(ij))/2.d0*(dt/2.)/dx2**2
      factyp=(kdint(bond_num(ij,5))+kdint(ij))/2.d0*(dt/2.)/dy1**2
      diag(ij)=1.d0+factyp	! 中心对角系数
      sup(ij)=-factyp				! 上对角系数
      inf(ij)=0			! 下对角系数
      f(ij)=zint(ij)+factxp*zint(bond_num(ij,4))
     1    +factxm*zint(bond_num(ij,3))-(factxp+factxm)*zint(ij)
      
cccc      elseif ((bi.eq.5).and.(bond_num(ij,1).eq.2)     ! 右边界
cccc     1                        .and.(cbc(2:2).ne.'1')) then
cccc      dx2=abs(length2(bond_num(ij,3),1)-length2(ij,1))
cccc      dy1=abs(length2(bond_num(ij,4),2)-length2(ij,2))
cccc      dy2=abs(length2(bond_num(ij,5),2)-length2(ij,2))
cccc      factxm=(kdint(bond_num(ij,3))+kdint(ij))/2.d0*(dt/2.)/dx2**2
cccc      factyp=(kdint(bond_num(ij,4))+kdint(ij))/2.d0*(dt/2.)/dy1**2
cccc      factym=(kdint(bond_num(ij,5))+kdint(ij))/2.d0*(dt/2.)/dy2**2
cccc      diag(ij)=1.d0+factyp+factym	! 中心对角系数
cccc      sup(ij)=-factyp				! 上对角系数
cccc      inf(ij)=-factym				! 下对角系数
cccc      f(ij)=zint(ij)+factxm*zint(bond_num(ij,3))-factxm*zint(ij)
      
      elseif ((bi.eq.5).and.(bond_num(ij,1).eq.3)     ! 上边界
     1                        .and.(cbc(3:3).ne.'1')) then
      dx1=abs(length2(bond_num(ij,3),1)-length2(ij,1))
      dx2=abs(length2(bond_num(ij,4),1)-length2(ij,1))
      dy2=abs(length2(bond_num(ij,5),2)-length2(ij,2))
      factxm=(kdint(bond_num(ij,3))+kdint(ij))/2.d0*(dt/2.)/dx1**2
      factxp=(kdint(bond_num(ij,4))+kdint(ij))/2.d0*(dt/2.)/dx2**2
      factym=(kdint(bond_num(ij,5))+kdint(ij))/2.d0*(dt/2.)/dy2**2
      diag(ij)=1.d0+factym	! 中心对角系数
      sup(ij)=0			! 上对角系数
      inf(ij)=-factym				! 下对角系数
      f(ij)=zint(ij)+factxp*zint(bond_num(ij,4))
     1    +factxm*zint(bond_num(ij,3))-(factxp+factxm)*zint(ij)
      
cccc      elseif ((bi.eq.5).and.(bond_num(ij,1).eq.4)     ! 左边界
cccc     1                        .and.(cbc(4:4).ne.'1')) then
cccc      dx1=abs(length2(bond_num(ij,3),1)-length2(ij,1))
cccc      dy1=abs(length2(bond_num(ij,4),2)-length2(ij,2))
cccc      dy2=abs(length2(bond_num(ij,5),2)-length2(ij,2))
cccc      factxp=(kdint(bond_num(ij,3))+kdint(ij))/2.d0*(dt/2.)/dx1**2
cccc      factyp=(kdint(bond_num(ij,4))+kdint(ij))/2.d0*(dt/2.)/dy1**2
cccc      factym=(kdint(bond_num(ij,5))+kdint(ij))/2.d0*(dt/2.)/dy2**2
cccc      diag(ij)=1.d0+factyp+factym	! 中心对角系数
cccc      sup(ij)=-factyp				! 上对角系数
cccc      inf(ij)=-factym				! 下对角系数
cccc      f(ij)=zint(ij)+factxp*zint(bond_num(ij,3))-factxp*zint(ij)
      
      else    ! 判断是固定边界不进行扩散计算
      diag(ij)=1.
      sup(ij)=0.
      inf(ij)=0.
      f(ij)=zint(ij)
c      write(6,*)'ij,bi,bond_num(ij,1)',ij,bi,bond_num(ij,1)
      end if
      end do
      
      call tridag (inf(1:nn),diag(1:nn),sup(1:nn),f(1:nn),res(1:nn),nn)
      do ij=1,nn
      zintp(ij)=res(ij)
      enddo
      

      
      ! stores result in 1D array
      ! 更新节点高程，剥蚀总量，剥蚀速率
      hh2(:nn)=zintp(:nn)-h(:nn)
      do ij=1,nn
      etot(ij)=etot(ij)+h(ij)-zintp(ij)
      erate(ij)=erate(ij)+(h(ij)-zintp(ij))/dt
      h(ij)=zintp(ij)
      enddo
      
c      write(6,*)'zintp(1:20)',zintp(1:20)
c      write(6,*)'h(1:20)',h(1:20)
      
      b(1:nn)=min(h(1:nn),b(1:nn))
cc      deallocate (f,diag,sup,inf,res)
cc      deallocate (zint,kdint,zintp)
      if (allocated(f)) deallocate(f)
      if (allocated(diag)) deallocate(diag)
      if (allocated(sup)) deallocate(sup)
      if (allocated(inf)) deallocate(inf)
      if (allocated(res)) deallocate(res)
      if (allocated(zint)) deallocate(zint)
      if (allocated(kdint)) deallocate(kdint)
      if (allocated(zintp)) deallocate(zintp)
      
      
      return
      
      end subroutine Diffusion
      
C----------
      ! 解三对角线性方程的子程序
      ! subroutine to solve a tri-diagonal system of equations (from Numerical Recipes)
      SUBROUTINE tridag(a,b,c,r,u,n)
      implicit none
      INTEGER n,j
      double precision a(n),b(n),c(n),r(n),u(n)
      double precision bet
      double precision,dimension(:),allocatable::gam
      allocate (gam(n))
      
      if(b(1).eq.0.d0) THEN
      WRITE(6,*) 'in tridag'
      END IF
      
      ! first pass
      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet.eq.0.) then
          WRITE(6,*) 'tridag failed'
        endif
        u(j)=(r(j)-a(j)*u(j-1))/bet
11    continue
      
        ! second pass
       do 12 j=n-1,1,-1
         u(j)=u(j)-gam(j+1)*u(j+1)
12     continue

      if (allocated(gam)) deallocate(gam)
      return
      
      END
