"""
Given the refractive indices at a certain wavelength this module
calculates the Mie scattering by a stratified sphere.The basic code used 
was that described in the report: " Subroutines for computing the parameters of 
the electromagnetic radiation scattered by a sphere " J.V. Dave,
I B M Scientific Center, Palo Alto , California.
Report NO. 320 - 3236 .. MAY 1968 .
Parameters
----------
RO : float
    Outer Shell Radius (cm)
RFR : float
    Float64 refractive index of shell layer (in the form n= RFR-i*RFI)
RFI : float
    Imaginary refractive index of shell layer (in the form n= RFR-i*RFI)
THET : ndarray 
    Angle in degrees between the directions of the incident and the scattered radiation.
JX : integer
    Total number of THET for which calculations are required
R : float
    Radius of core (cm)`
RE2 : float 
    Float64 refractive index of core (in the form n= RE2-i*TMAG2)
TMAG2 : float
    Imaginary refractive index of core (in the form n= RE2-i*TMAG2)
    
WVNO : float
    Wave-number corresponding to the wavelength. (cm^-1)

Returns
-------
QEXT: float
    Efficiency factor for extinction,VAN DE HULST,P.14 ' 127
QSCAT: float
    Efficiency factor for scattering,VAN DE HULST,P.14 ' 127
CTBQRS: float
    Average(cos(theta))*QSCAT,VAN DE HULST,P.14 ' 127
ISTATUS: integer
    Convergence indicator, 0 if converged, -1 if otherwise.    
"""
function fort_mie_calc(
    RO, RFR, RFI, THET, JX, R, RE2, TMAG2, WVNO;
    EPSILON_MIE=1e-7,nacap=1000000, IT=1
)

    ACAP = zeros(Complex{Float64}, nacap)
    W = zeros(Complex{Float64}, 3, nacap)
    WFN = zeros(Complex{Float64}, 2)
    Z = zeros(Complex{Float64}, 4)
    U = zeros(Complex{Float64}, 8)
    T = zeros(Complex{Float64}, 5)
    TA = zeros(4)
    TB = zeros(2)
    TC = zeros(2)
    TD = zeros(2) 
    TE = zeros(2)  
    PI = zeros(3, IT)
    TAU = zeros(3, IT)
    CSTHT = zeros(IT)
    THETD = fill(THET, (IT))
    SI2THT = zeros(IT)
    ELTRMX = zeros(4, IT, 2)
    IFLAG = 1
    if  R/RO < 1e-6  
        IFLAG = 2
    end
    @assert JX > IT "THE VALUE OF THE ARGUMENT JX IS GREATER THAN IT. PLEASE READ COMMENTS."
        
    RF = RFR - RFI*im
    RC =  RE2 - TMAG2*im
    X  =  RO * WVNO
    K1 =  RC * WVNO
    K2 =  RF * WVNO
    K3 =  WVNO + 0*im
    Z[1] =  K2 * RO
    Z[2] =  K3 * RO
    Z[3] =  K1 * R
    Z[4] =  K2 * R
    X1   =  real(Z[0])
    X4   =  real(Z[3])
    Y1   =  imag(Z[0])
    Y4   =  imag(Z[3])
    RRF  =  1.0 / RF
    RX   =  1.0 / X
    RRFX =  RRF * RX
    T[1] = X * abs(RF)
    NMX1 = int(1.10 * T[1])
    if  NMX1 > nacap-1
        return 0, 0, 0, -1
    end
        
    NMX2 = Int(T[1])
    if  NMX1 <= 150
        NMX1 = 150
        NMX2 = 135
    end

    ACAP[NMX1] = zero(Complex)
    if IFLAG != 2
        for N in 1:3
            W[N, NMX1] = zero(Complex)
        end
    end
    for N in 0:NMX1-1
        NN = NMX1 - N-1 ## removed a plus 1 to make up for python loop
        # I should put this back in
        ACAP[NN] = (NN+2) * RRFX - 1.0 / ((NN+2) * RRFX + ACAP[NN+1])
        if IFLAG != 2
            for M in 1:3
                W[M,NN] = ((NN+2) / Z[M+1])  -1.0 / ((NN+2) / Z[M+1]  +  W[M,NN+1])
            end
        end
    for J in 1:JX
        THETD[J] = abs(THETD[J])
        if  THETD[J] < 90.0
            T[1]     =  (pi * THETD[J]) / 180.0
            CSTHT[J] =  cos(T[1])
            SI2THT[J] =  1.0 - CSTHT[J]^2
        end
        if  THETD[J] == 90.0
            CSTHT[J]  =  0.0
            SI2THT[J] =  1.0
        end
        @assert THETD[J] <= 90.0 "THE VALUE OF THE SCATTERING ANGLE IS GREATER THAN 90.0 DEGREES. PLEASE READ COMMENTS"
    end

    for J in 1:JX
        PI[1,J]  =  0.0
        PI[2,J]  =  1.0
        TAU[1,J] =  0.0
        TAU[2,J] =  CSTHT[J]
    end
    T[1]   =  cos(X)
    T[2]   =  sin(X)
    WM1    =  T[1]-T[2]*im
    WFN[1] =  T[2]+T[1]*im
    TA[1]  =  T[2]
    TA[2]  =  T[1]
    WFN[2] =  RX * WFN[1] - WM1
    TA[3]  =  real(WFN[2])
    TA[4]  =  imag(WFN[2])
    if IFLAG != 2
        N = 1
        SINX1   =  sin(X1)
        SINX4   =  sin(X4)
        COSX1   =  cos(X1)
        COSX4   =  cos(X4)
        EY1     =  exp(Y1)
        E2Y1    =  EY1 * EY1
        EY4     =  exp(Y4)
        EY1MY4  =  exp(Y1 - Y4)
        EY1PY4  =  EY1 * EY4
        EY1MY4  =  exp(Y1 - Y4)
        AA  =  SINX4 * (EY1PY4 + EY1MY4)
        BB  =  COSX4 * (EY1PY4 - EY1MY4)
        CC  =  SINX1 * (E2Y1 + 1.0)
        DD  =  COSX1 * (E2Y1 - 1.0)
        DENOM   =  1.0  +  E2Y1 * (4.0 * SINX1 * SINX1 - 2.0 + E2Y1)
        REALP   =  (AA * CC  +  BB * DD) / DENOM
        AMAGP   =  (BB * CC  -  AA * DD) / DENOM
        DUMMY   =  complex(REALP, AMAGP)
        AA  =  SINX4 * SINX4 - 0.5
        BB  =  COSX4 * SINX4
        P24H24  =  0.5 + complex(AA, BB) * EY4 * EY4
        AA  =  SINX1 * SINX4  -  COSX1 * COSX4
        BB  =  SINX1 * COSX4  +  COSX1 * SINX4
        CC  =  SINX1 * SINX4  +  COSX1 * COSX4
        DD  = -SINX1 * COSX4  +  COSX1 * SINX4
        P24H21  =  0.5 * complex(AA, BB) * EY1 * EY4 + 0.5 * complex(CC, DD) * EY1MY4
        DH4  =  Z[4] / (1.0 + 1.0im * Z[4])  -  1.0 / Z[4]
        DH1  =  Z[1] / (1.0 + 1.0im * Z[1])  -  1.0 / Z[1]
        DH2  =  Z[2] / (1.0 + 1.0im * Z[2])  -  1.0 / Z[2]
        PSTORE  =  (DH4 + N / Z[4])  *  (W[3,N] + N / Z[4])
        P24H24  =  P24H24 / PSTORE
        HSTORE  =  (DH1 + N / Z[1])  *  (W[3,N] + N / Z[4])
        P24H21  =  P24H21 / HSTORE
        PSTORE  =  (ACAP[N] + N / Z[1])  /  (W[3,N] + N / Z[4])
        DUMMY   =  DUMMY * PSTORE
        DUMSQ   =  DUMMY * DUMMY
        U[1] =  K3 * ACAP[N]  -  K2 * W[1,N]
        U[2] =  K3 * ACAP[N]  -  K2 * DH2
        U[3] =  K2 * ACAP[N]  -  K3 * W[1,N]
        U[4] =  K2 * ACAP[N]  -  K3 * DH2
        U[5] =  K1 *  W[3,N]  -  K2 * W[2,N]
        U[6] =  K2 *  W[3,N]  -  K1 * W[2,N]
        U[7] =  -1im  *  (DUMMY * P24H21 - P24H24)
        U[8] =  TA[3] / WFN[2]
            
        FNA  =  U[8] * (U[1]*U[5]*U[7]  +  K1*U[1]  -  DUMSQ*K3*U[5]) /(U[2]*U[5]*U[7]  +  K1*U[2]  -  DUMSQ*K3*U[5])
        FNB  =  U[8] * (U[3]*U[6]*U[7]  +  K2*U[3]  -  DUMSQ*K2*U[6]) /(U[4]*U[6]*U[7]  +  K2*U[4]  -  DUMSQ*K2*U[6])
        TB[1] = real(FNA)
        TB[2] = imag(FNA)
        TC[1] = real(FNB)
        TC[2] = imag(FNB)
    elseif IFLAG == 2
        TC1  =  ACAP[1] * RRF  +  RX
        TC2  =  ACAP[1] * RF   +  RX
        FNA  =  (TC1 * TA[3]  -  TA[1]) / (TC1 * WFN[2]  -  WFN[1])
        FNB  =  (TC2 * TA[3]  -  TA[1]) / (TC2 * WFN[2]  -  WFN[1])
        TB[1] = real(FNA)
        TB[2] = imag(FNA)
        TC[1] = real(FNB)
        TC[2] = imag(FNB)
    end

    FNAP = FNA
    FNBP = FNB
    TD[1] = real(FNAP)
    TD[2] = imag(FNAP)
    TE[1] = real(FNBP)
    TE[2] = imag(FNBP)
    T[1] = 1.50
        
    TB[1] = T[1] * TB[1]
    TB[2] = T[1] * TB[2]
    TC[1] = T[1] * TC[1]
    TC[2] = T[1] * TC[2]
        
    for J in 1:JX
        ELTRMX[1,J,1] = TB[1] * PI[2,J] + TC[1] * TAU[2,J]
        ELTRMX[2,J,1] = TB[2] * PI[2,J] + TC[2] * TAU[2,J]
        ELTRMX[3,J,1] = TC[1] * PI[2,J] + TB[1] * TAU[2,J]
        ELTRMX[4,J,1] = TC[2] * PI[2,J] + TB[2] * TAU[2,J]
        ELTRMX[1,J,2] = TB[1] * PI[2,J] - TC[1] * TAU[2,J]
        ELTRMX[2,J,2] = TB[2] * PI[2,J] - TC[2] * TAU[2,J]
        ELTRMX[3,J,2] = TC[1] * PI[2,J] - TB[1] * TAU[2,J]
        ELTRMX[4,J,2] = TC[2] * PI[2,J] - TB[2] * TAU[2,J]
    end

    QEXT   = 2.0 * (TB[1] + TC[1])
    QSCAT  = (TB[1]^2 + TB[2]^2 + TC[1]^2 + TC[2]^2) / 0.75
    CTBRQS = 0.0
    QBSR   = -2.0*(TC[1] - TB[1])
    QBSI   = -2.0*(TC[2] - TB[2])
    RMM    = -1.0
    N = 2

    while N <= NMX2
        T[1] = 2*N - 1
        T[2] =   N - 1
        T[3] = 2*N + 1

        for  J in 1:JX
            PI[3,J]  = (T[1] * PI[2,J] * CSTHT[J] - N * PI[1,J]) / T[2]
            TAU[3,J] = CSTHT[J] * (PI[3,J] - PI[1,J])  - T[1] * SI2THT[J] * PI[2,J]  +  TAU[1,J]
        end
        WM1    =  WFN[1]
        WFN[1] =  WFN[2]
        TA[1]  =  real(WFN[1])

        TA[2]  =  imag(WFN[1])
        TA[4]  =  imag(WFN[2])
        WFN[2] =  T[1] * RX * WFN[1]  -  WM1
        TA[3]  =  real(WFN[2])
        
        if IFLAG != 2
            DH2  =  - N / Z[2]  +  1.0 / ((N / Z[2]) - DH2)
            DH4  =  - N / Z[4]  +  1.0 / ((N / Z[4]) - DH4)
            DH1  =  - N / Z[1]  +  1.0 / ((N / Z[1]) - DH1)
            PSTORE  =  (DH4 + (N / Z[4]))  *  (W[3,N] + (N / Z[4]))
            P24H24  =  P24H24 / PSTORE
            HSTORE  =  (DH1 + (N / Z[1]))  *  (W[3,N] + (N / Z[4]))
            P24H21  =  P24H21 / HSTORE
            PSTORE  =  (ACAP[N] + (N / Z[1]))  /  (W[3,N] + (N / Z[4]))
            DUMMY   =  DUMMY * PSTORE
            DUMSQ   =  DUMMY * DUMMY
            U[1] =  K3 * ACAP[N]  -  K2 * W[1,N]
            U[2] =  K3 * ACAP[N]  -  K2 * DH2
            U[3] =  K2 * ACAP[N]  -  K3 * W[1,N]
            U[4] =  K2 * ACAP[N]  -  K3 * DH2
            U[5] =  K1 *  W[3,N]  -  K2 * W[2,N]
            U[6] =  K2 *  W[3,N]  -  K1 * W[2,N]
            U[7] =  -1im  *  (DUMMY * P24H21 - P24H24)
            U[8] =  TA[3] / WFN[2]

            FNA  =  U[8] * (U[1]*U[5]*U[7]  +  K1*U[1]  -  DUMSQ*K3*U[5]) /(U[2]*U[5]*U[7]  +  K1*U[2]  -  DUMSQ*K3*U[5])
            FNB  =  U[8] * (U[3]*U[6]*U[7]  +  K2*U[3]  -  DUMSQ*K2*U[6]) /(U[4]*U[6]*U[7]  +  K2*U[4]  -  DUMSQ*K2*U[6])
            TB[0] = real(FNA)
            TB[1] = imag(FNA)
            TC[0] = real(FNB)
            TC[1] = imag(FNB)
        end

        TC1  =  ACAP[N] * RRF  +  N * RX
        TC2  =  ACAP[N] * RF   +  N * RX
        FN1  =  (TC1 * TA[3]  -  TA[1]) /  (TC1 * WFN[2] - WFN[1])
        FN2  =  (TC2 * TA[3]  -  TA[1]) /  (TC2 * WFN[2] - WFN[1])
        M    =  Int(round(WVNO * R))
        if  N >= M
            if IFLAG == 2
                FNA  =  FN1
                FNB  =  FN2
                TB[1] = real(FNA)
                TB[2] = imag(FNA)
                TC[1] = real(FNB)
                TC[2] = imag(FNB)
            end

            if (IFLAG != 2) && (abs((FN1-FNA) / FN1) < EPSILON_MIE) && (abs((FN2-FNB) / FN2)  < EPSILON_MIE)
                IFLAG = 2
            end
        end
    
        T[5]  =  N
        T[4]  =  T[1] / (T[5] * T[2])
        T[2]  =  (T[2] * (T[5] + 1.0)) / T[5]

        CTBRQS +=  T[2] * (TD[1] * TB[1]  +  TD[2] * TB[2] + TE[1] * TC[1]  +  TE[2]* TC[2])+T[4] * (TD[1] * TE[1]  +  TD[2] * TE[2])
        QEXT   +=    T[3] * (TB[1] + TC[1])
        T[3]    =  TB[1]^2 + TB[2]^2 + TC[1]^2 + TC[2]^2
        QSCAT  +=  T[3] * T[4]
        RMM     =  -RMM
        QBSR +=  T[3]*RMM*(TC[1] - TB[1])
        QBSI  +=  T[3]*RMM*(TC[2] - TB[2])

        T[2]    =  N * (N+1)
        T[1]    =  T[3] / T[2]
        K = Int(round(N))
        for J in 1:JX
            ELTRMX[1,J,1] += T[1]*(TB[1]*PI[3,J]+TC[1]*TAU[3,J])
            ELTRMX[2,J,1] += T[1]*(TB[2]*PI[3,J]+TC[2]*TAU[3,J])
            ELTRMX[3,J,1] += T[1]*(TC[1]*PI[3,J]+TB[1]*TAU[3,J])
            ELTRMX[4,J,1] += T[1]*(TC[2]*PI[3,J]+TB[2]*TAU[3,J])
            if K % 2 == 0
                ELTRMX[1,J,2] += T[1]*(-TB[1]*PI[3,J]+TC[1]*TAU[3,J])
                ELTRMX[2,J,2] += T[1]*(-TB[2]*PI[3,J]+TC[2]*TAU[3,J])
                ELTRMX[3,J,2] += T[1]*(-TC[1]*PI[3,J]+TB[1]*TAU[3,J])
                ELTRMX[4,J,2] += T[1]*(-TC[2]*PI[3,J]+TB[2]*TAU[3,J])
            else
                ELTRMX[1,J,2] += T[1]*(TB[1]*PI[3,J]-TC[1]*TAU[3,J])
                ELTRMX[2,J,2] += T[1]*(TB[2]*PI[3,J]-TC[2]*TAU[3,J])
                ELTRMX[3,J,2] += T[1]*(TC[1]*PI[3,J]-TB[1]*TAU[3,J])
                ELTRMX[4,J,2] += T[1]*(TC[2]*PI[3,J]-TB[2]*TAU[3,J])
            end
        end
        
        if  T[4] >= EPSILON_MIE
            N += 1
            for J in 1:JX
                PI[1,J]   =   PI[2,J]
                PI[2,J]   =   PI[3,J]
                TAU[1,J]  =  TAU[2,J]
                TAU[2,J]  =  TAU[3,J]
            end
            FNAP  =  FNA
            FNBP  =  FNB
            TD[1] = real(FNAP)
            TD[2] = imag(FNAP)
            TE[1] = real(FNBP)
            TE[2] = imag(FNBP)
        else
            break 
        end
    end # while
    if  N >= NMX2
        return 0, 0, 0, -1
    end
    for J in 1:JX
        for K in 1:2
            for I in 1:5
                T[I]  =  ELTRMX[I,J,K]
            end
            ELTRMX[2,J,K]  =      T[1]^2  +  T[2]^2
            ELTRMX[1,J,K]  =      T[3]^2  +  T[4]^2
            ELTRMX[3,J,K]  =  T[1] * T[3]  +  T[2] * T[4]
            ELTRMX[4,J,K]  =  T[2] * T[3]  -  T[4] * T[1]
        end
    end

    T[1]    =    2.0 * RX^2
    QEXT    =   QEXT * T[1]
    QSCAT   =  QSCAT * T[1]
    CTBRQS  =  2.0 * CTBRQS * T[1]           
    istatus = 0

    return QEXT,QSCAT,CTBRQS,istatus
end