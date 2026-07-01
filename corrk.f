      MODULE corrkmodule
        IMPLICIT NONE
        include 'nwno.inc'
        ! COMMON BLOCK SHAPES HERE:
        REAL :: OPAC_CORRK(NTGRID, NPGRID, NWNO, 8)
        REAL :: TS_CORRK(NTGRID), PS_CORRK(NPGRID), TS_LOG_CORRK(NTGRID), PS_LOG_CORRK(NPGRID), WGTS_CORRK(8) 
        REAL :: WNO_EDGES(NWNO+1), WNO_CTRS(NWNO), STEL_SPEC(NWNO), INT_SPEC(NWNO)
        REAL :: OPAC_CIA(NTGRID, NPGRID, NWNO), TAURAY_PER_DPG(NTGRID, NPGRID, NWNO)
        REAL :: PLANCK_INTS(NWNO, 3925), PLANCK_TS(3925)
        INTEGER :: MINWNOSTEL
        real :: CLOUD_KEXT(13, 100, NWNO), CLOUD_A(13, 100, NWNO), CLOUD_G(13, 100, NWNO)
        

        CONTAINS
        SUBROUTINE corrk_setup(METALLICITY, C_TO_O, FBASEFLUX, GASCON, with_TiO_and_VO,TS_CORRK, PS_CORRK, 
     &           TS_LOG_CORRK, PS_LOG_CORRK, WGTS_CORRK, WNO_EDGES, WNO_CTRS, STEL_SPEC, INT_SPEC, TAURAY_PER_DPG,
     &           OPAC_CORRK, PLANCK_INTS, PLANCK_TS, NWNO, MINWNOSTEL, CLOUD_KEXT, CLOUD_A, CLOUD_G)
            implicit none
            integer :: NWNO
            INTEGER :: I, J, K, L
            REAL :: METALLICITY, C_TO_O, FBASEFLUX
            REAL :: TINT
            INTEGER :: TINT_INDEX
            REAL :: GASCON, MMW
            REAL :: PLANCK_INTS(NWNO, 3925), PLANCK_TS(3925)
            REAL :: with_TiO_and_VO

            CHARACTER(LEN=32) :: dummyMETALLICITY_str, dummyC_TO_O_str, tiovo_str
            CHARACTER(len=30) :: FOLDER
            CHARACTER(len=80) :: FNM, CIAFNM, RAYFNM
            REAL :: TEMP_OPAC_CORRK(8,NTGRID*NPGRID*NWNO) ! 8 gauss pts, NTGRID T, NPGRID P, NWNO wavelengths
            REAL :: TEMP_OPAC_CIA(NWNO,NTGRID*NPGRID) ! NTGRID T, NPGRID P, NWNO wavelengths
            REAL :: TEMP_TAURAY_PER_DPG_PERMU(NWNO, NTGRID*NPGRID)

            ! COMMON BLOCK SHAPES HERE:
            REAL :: OPAC_CORRK(NTGRID, NPGRID, NWNO, 8)
            REAL :: TS_CORRK(NTGRID), PS_CORRK(NPGRID), TS_LOG_CORRK(NTGRID), PS_LOG_CORRK(NPGRID), WGTS_CORRK(8) 
            REAL :: WNO_EDGES(NWNO+1), WNO_CTRS(NWNO), STEL_SPEC(NWNO), INT_SPEC(NWNO)
            REAL :: OPAC_CIA(NTGRID, NPGRID, NWNO), TAURAY_PER_DPG(NTGRID, NPGRID, NWNO)
            INTEGER :: MINWNOSTEL

            ! Declare cloud opacity stuff:
            real :: CLOUD_KEXT(13, 100, NWNO), CLOUD_A(13, 100, NWNO), CLOUD_G(13, 100, NWNO)

            MMW = 8314.462618 / GASCON ! Mean molecular weight in g/mol (or H masses per molecule)

            ! Convert METALLICITY and C_TO_O to strings
            WRITE(dummyMETALLICITY_str, '(SP,I4.3)') NINT(METALLICITY * 100)
            WRITE(dummyC_TO_O_str, '(I3.3)') NINT(C_TO_O * 100)
            write(*,*) 'C to O ratio: ', C_TO_O
            write(*,*) 'C to O string: ', dummyC_TO_O_str
            ! write(*,*) 'with_TiO_and_VO: ', with_TiO_and_VO

            if (with_TiO_and_VO .eq. 1.) then
                tiovo_str = '_witiovo'
            else if (with_TiO_and_VO .eq. 2.) then
                tiovo_str = '_notiovo'
            else
                write(*,*) "tiovo flag not recognized, setting with_TiO_and_VO = 1 for the first loop"
                tiovo_str = '_witiovo'
            end if
            ! write(*,*) 'METALLICITY: ', dummyMETALLICITY_str
            ! write(*,*) 'C_TO_O: ', dummyC_TO_O_str
            IF (NWNO.eq.11) THEN
                FOLDER = '../11bands/'
            ELSE IF (NWNO.eq.30) THEN
                FOLDER = '../30bands/'
            ELSE
                WRITE(*,*) 'NWNO is wrong in nwno.inc, options are 11 and 30'
            STOP
                END IF
            ! FOLDER = '../11bands/'
            FNM= trim(FOLDER)//"kcoeffs/feh"//trim(dummyMETALLICITY_str)//"_co_"//trim(dummyC_TO_O_str)//trim(tiovo_str)//".txt"
            CIAFNM= trim(FOLDER)//"continuum/cia+hmin_feh"//trim(dummyMETALLICITY_str)//"_co_"//trim(dummyC_TO_O_str)
     &            //trim(tiovo_str)//".txt"
            RAYFNM= trim(FOLDER)//"rayleigh/rayleigh_feh"//trim(dummyMETALLICITY_str)//"_co_"//trim(dummyC_TO_O_str)
     &              //trim(tiovo_str)//".txt"
            ! Read in k-coefficients from file
            write(*,*) 'Reading in k-coefficients from file: ', FNM
            write(*,*) 'Reading in CIA opacities from file: ', CIAFNM
            write(*,*) 'Reading in Rayleigh opacities from file: ', RAYFNM
            OPEN(UNIT=1, FILE=FNM)
            READ(1,*) TEMP_OPAC_CORRK
            CLOSE(1)
            OPEN(UNIT=1, FILE=CIAFNM)
            READ(1,*) TEMP_OPAC_CIA
            CLOSE(1)
            OPEN(UNIT=1, FILE=RAYFNM)
            READ(1,*) TEMP_TAURAY_PER_DPG_PERMU
            CLOSE(1)

            ! Reshape into 3D array
            ! Index order of new array: T, P, wavelength (no gauss pt because gray in each band)
            DO I = 1, NTGRID
                DO J = 1, NPGRID
                    DO K = 1, NWNO
                    ! after this line we have m^2/kg
                    ! file is in units of m^2/molecule
                    OPAC_CIA(I,J,K) = TEMP_OPAC_CIA(K,(I-1)*NPGRID+J) * 6.022e23 / (MMW/1000)
                    ! after this line we should have log10(tau / (deltap / gravity), all in SI units)
                    TAURAY_PER_DPG(I,J,K) = log10(TEMP_TAURAY_PER_DPG_PERMU(K,(I-1)*NPGRID+J) / MMW)
                    END DO
                END DO
            END DO

            ! Reshape into 4D array
            ! Index order of new array: T, P, wavelength, gauss pt
            DO I = 1, NTGRID
                DO J = 1, NPGRID
                    DO K = 1, NWNO
                        DO L = 1, 8
                            ! After exponential, we're in cm^2/molecule units which we convert to m^2/molecule since RT uses SI units
                            ! avagadro + MMW/1000 conversion factor get us to m^2/kg
                            OPAC_CORRK(I,J,K,L) = EXP(TEMP_OPAC_CORRK(L,(I-1)*NPGRID*NWNO+(J-1)*NWNO+K)) * 1e-4 * 6.022e23 / 
     &                                              (MMW/1000)
                            ! Add in CIA and H- opacities
                            OPAC_CORRK(I,J,K,L) = OPAC_CORRK(I,J,K,L) + OPAC_CIA(I,J,K)
                            ! Convert to log space
                            OPAC_CORRK(I,J,K,L) = LOG10(OPAC_CORRK(I,J,K,L))
                        END DO
                    END DO
                END DO
            END DO

            ! Set up the gauss weights of the 8 gauss points
            WGTS_CORRK = (/0.16523105, 0.30976895, 0.30976895, 0.16523105, 0.00869637,
     &              0.01630363, 0.01630363, 0.00869637/)

            ! Read in the temperature and pressure points of our grid
            OPEN(UNIT=1, FILE=trim(FOLDER)//"kcoeffs/Tpts.txt")
            READ(1,*) TS_CORRK
            CLOSE(1)

            OPEN(UNIT=1, FILE=trim(FOLDER)//"kcoeffs/Ppts.txt")
            READ(1,*) PS_CORRK
            CLOSE(1)
            PS_CORRK = PS_CORRK * 1.0E5 ! Convert bars to Pascals, which the RT uses

            ! Invert T and log P for interpolation
            TS_LOG_CORRK = LOG10(TS_CORRK)
            PS_LOG_CORRK = LOG10(PS_CORRK)

            ! Read in the wavenumber points of our grid
            OPEN(UNIT=1, FILE=trim(FOLDER)//"kcoeffs/wno_edges.txt")
            READ(1,*) WNO_EDGES
            CLOSE(1)

            OPEN(UNIT=1, FILE=trim(FOLDER)//"kcoeffs/wno_ctrs.txt")
            READ(1,*) WNO_CTRS
            CLOSE(1)

            ! Read cloud opacities:
            CALL CLOUD_SETUP(NWNO, FOLDER, CLOUD_KEXT, CLOUD_A, CLOUD_G)

            ! Read in the stellar spectrum
            CALL BIN_STELLAR_SPECTRUM(trim(FOLDER)//"kcoeffs/stellar_spectrum.txt", STEL_SPEC, WNO_EDGES, NWNO)
            ! WRITE(*,*) 'STELLAR SPECTRUM: ', STEL_SPEC
            ! WRITE(*,*) 'total stellar flux: ', SUM(STEL_SPEC)

            ! Grab Planck function integrals for each tmeperature and bin:
            OPEN(UNIT=1, FILE=trim(FOLDER)//"kcoeffs/planck_integrals.txt")
            READ(1,*) PLANCK_INTS
            CLOSE(1)
            OPEN(UNIT=1, FILE=trim(FOLDER)//"kcoeffs/planck_temps.txt")
            READ(1,*) PLANCK_TS
            CLOSE(1)
            ! window, temp
            ! write(*,*) 'Planck integrals: ', PLANCK_INTS(1,1), PLANCK_INTS(1,2), PLANCK_INTS(2,1)
            ! write(*,*) 'Planck temps: ', PLANCK_TS(1), PLANCK_TS(2), PLANCK_TS(3)

            ! Now split up internal flux by channel:
            TINT = (FBASEFLUX / 5.670367E-8) ** 0.25 ! Convert flux to temperature
            TINT_INDEX = NEAREST_INDEX(PLANCK_TS, 3925, TINT)
            INT_SPEC = PLANCK_INTS(:, TINT_INDEX)
            INT_SPEC = INT_SPEC / SUM(INT_SPEC)
            write(*,*) 'TINT: ', TINT, 'TINT_INDEX: ', TINT_INDEX
            write(*,*) 'BINNED INTERIOR SPECTRUM: ', INT_SPEC
            write(*,*) 'BINNED STELLAR SPECTRUM: ', STEL_SPEC
            MINWNOSTEL = 0
            ! the number in the if statement is the fraction of starlight you want to capture
            DO I=1, NWNO
                IF (SUM(STEL_SPEC(I:)) > 1.) then
                    MINWNOSTEL = I-1
                END IF
            END DO
            ! MINWNOSTEL = 1
            WRITE(*,*) 'WNO MIN: ', MINWNOSTEL
        END SUBROUTINE corrk_setup

        SUBROUTINE BIN_STELLAR_SPECTRUM(file_name, STEL_SPEC, WNO_EDGES, NWNO)
            implicit none
            ! include 'nwno.inc'
            integer :: NWNO
            CHARACTER(len=*) :: file_name
            REAL :: STEL_SPEC(NWNO), STEL_SPEC_IN(2, 5000)
            REAL :: WNO_EDGES(NWNO+1)
            INTEGER :: I, J
            REAL :: WNO(5000), FLUX(5000)
            ! Read in the stellar spectrum
            OPEN(UNIT=1, FILE=file_name)
            READ(1,*) STEL_SPEC_IN
            CLOSE(1)
            DO I = 1, 5000
                WNO(I) = STEL_SPEC_IN(1,I)
                FLUX(I) = STEL_SPEC_IN(2,I)
            END DO
            ! Trapezoidal integration in each spectral band
            DO J = 1, NWNO
                STEL_SPEC(J) = 0.0
                DO I = 1, 4999
                    IF (WNO(I) >= WNO_EDGES(J) .AND. WNO(I) < WNO_EDGES(J+1)) THEN
                        STEL_SPEC(J) = STEL_SPEC(J) + ((FLUX(I) + FLUX(I+1))/2 *((WNO(I+1) - WNO(I)) * 3e10))
                    END IF
                END DO
            END DO

            ! Normalize the spectrum so the fort.7's SOLC_IN will define the total flux:
            STEL_SPEC = STEL_SPEC / SUM(STEL_SPEC)

        END SUBROUTINE BIN_STELLAR_SPECTRUM

        SUBROUTINE CLOUD_SETUP(NWNO, FOLDER, CLOUD_KEXT, CLOUD_A, CLOUD_G)
            implicit none
            integer :: NWNO
            real :: CLOUD_KEXT(13, 100, NWNO) ! species, radius, wavenumber bin 
            real :: CLOUD_A(13, 100, NWNO)
            real :: CLOUD_G(13, 100, NWNO)
            real :: DUMMY(NWNO, 100) ! to hold each array as it's read in
            character(len=30) :: FOLDER
            character(len=15) :: CLOUD_SPECIES(13)
            integer :: I, J, K
            
            
            CLOUD_SPECIES(1)  = 'KCl'
            CLOUD_SPECIES(2)  = 'ZnS'
            CLOUD_SPECIES(3)  = 'Na2S'
            CLOUD_SPECIES(4)  = 'MnS'
            CLOUD_SPECIES(5)  = 'Cr'
            CLOUD_SPECIES(6)  = 'SiO2_amorph'
            CLOUD_SPECIES(7)  = 'Mg2SiO4_amorph'
            CLOUD_SPECIES(8)  = 'VO'
            CLOUD_SPECIES(9)  = 'Ni'
            CLOUD_SPECIES(10) = 'Fe'
            CLOUD_SPECIES(11) = 'CaSiO4'
            CLOUD_SPECIES(12) = 'CaTiO3'
            CLOUD_SPECIES(13) = 'Al2O3'

            ! Read in the cloud kext values
            DO I = 1, 13
                OPEN(UNIT=1, FILE=trim(FOLDER)//"clouds/"//trim(CLOUD_SPECIES(I))//"_corrk_kext.txt")
                READ(1,*) DUMMY
                CLOSE(1)
                DO J = 1, 100
                    DO K = 1, NWNO
                        CLOUD_KEXT(I,J,K) = DUMMY(K,J) !TODO: figure out what units this should be in
                        ! tested, indices are what you expect (species, radius, wavenumber)
                    END DO
                END DO
            END DO
            ! Read in the cloud albedos
            DO I = 1, 13
                OPEN(UNIT=1, FILE=trim(FOLDER)//"clouds/"//trim(CLOUD_SPECIES(I))//"_corrk_pi0.txt")
                READ(1,*) DUMMY
                CLOSE(1)
                DO J = 1, 100
                    DO K = 1, NWNO
                        CLOUD_A(I,J,K) = DUMMY(K,J) !TODO: figure out what units this should be in
                        ! tested, indices are what you expect (species, radius, wavenumber)
                    END DO
                END DO
            END DO
            ! Read in the cloud asymmetry factors
            DO I = 1, 13
                OPEN(UNIT=1, FILE=trim(FOLDER)//"clouds/"//trim(CLOUD_SPECIES(I))//"_corrk_gg.txt")
                READ(1,*) DUMMY
                CLOSE(1)
                DO J = 1, 100
                    DO K = 1, NWNO
                        CLOUD_G(I,J,K) = DUMMY(K,J) !TODO: figure out what units this should be in
                        ! tested, indices are what you expect (species, radius, wavenumber)
                    END DO
                END DO
            END DO
        END SUBROUTINE CLOUD_SETUP

        INTEGER FUNCTION NEAREST_INDEX(ARR, N, VAL)
            ! Explicit-loop nearest-value lookup. Replaces MINLOC(ABS(ARR-VAL),1),
            ! which forces the compiler to materialize a temporary array; under
            ! -recursive these temporaries can land in shared scratch storage and
            ! be clobbered by other OpenMP threads, corrupting the result. This
            ! matters most for large lookup tables (e.g. PLANCK_TS, size 3925)
            ! that are queried every layer/column inside the parallel region.
            INTEGER, INTENT(IN) :: N
            REAL, INTENT(IN) :: ARR(N), VAL
            INTEGER, AUTOMATIC :: KK
            REAL, AUTOMATIC :: BESTDIFF, DIFFVAL
            NEAREST_INDEX = 1
            BESTDIFF = ABS(ARR(1) - VAL)
            DO KK = 2, N
                DIFFVAL = ABS(ARR(KK) - VAL)
                IF (DIFFVAL .LT. BESTDIFF) THEN
                    BESTDIFF = DIFFVAL
                    NEAREST_INDEX = KK
                ENDIF
            END DO
        END FUNCTION NEAREST_INDEX

      END MODULE corrkmodule





    !     SUBROUTINE get_gas_opacity_corrk_wrapper
    !         include 'rcommons.h'
    !         ! include 'nwno.inc'
    !         WRITE(*,*) "In get_gas_opacity_corrk_wrapper"
    !         call get_gas_opacity_corrk(METALLICITY,0.0, FBASEFLUX, GASCON, with_TiO_and_VO,NWNO)
    !     END SUBROUTINE get_gas_opacity_corrk_wrapper

    !     SUBROUTINE get_gas_opacity_corrk(METALLICITY,C_TO_O, FBASEFLUX, GASCON, with_TiO_and_VO,NWNO)
    !         implicit none
    !         COMMON/CORRKGAS/OPAC_CORRK, TS_CORRK, PS_CORRK, TS_LOG_CORRK, PS_LOG_CORRK, WGTS_CORRK, WNO_EDGES, WNO_CTRS, STEL_SPEC,
    !  &                    INT_SPEC, TAURAY_PER_DPG
    !         COMMON/PLANCK_INT/PLANCK_INTS, PLANCK_TS
    !         ! include 'nwno.inc'
    !         ! integer :: NWNO
    !         ! parameter (NWNO=11)
    !         integer :: NWNO
    !         INTEGER :: I, J, K, L
    !         REAL, INTENT(IN) :: METALLICITY, C_TO_O, FBASEFLUX
    !         REAL :: TINT
    !         INTEGER :: TINT_INDEX
    !         REAL :: GASCON, MMW
    !         REAL :: PLANCK_INTS(NWNO, 3925), PLANCK_TS(3925)
    !         REAL :: with_TiO_and_VO

    !         CHARACTER(LEN=32) :: dummyMETALLICITY_str, dummyC_TO_O_str, tiovo_str
    !         CHARACTER(len=30) :: FOLDER
    !         CHARACTER(len=80) :: FNM, CIAFNM, RAYFNM
    !         REAL :: TEMP_OPAC_CORRK(8,73*20*NWNO) ! 8 gauss pts, 73 T, 20 P, NWNO wavelengths
    !         REAL :: TEMP_OPAC_CIA(NWNO,73*20) ! 73 T, 20 P, NWNO wavelengths
    !         REAL :: TEMP_TAURAY_PER_DPG_PERMU(NWNO, 73*20)

    !         ! COMMON BLOCK SHAPES HERE:
    !         REAL :: OPAC_CORRK(73, 20, NWNO, 8)
    !         REAL :: TS_CORRK(73), PS_CORRK(20), TS_LOG_CORRK(73), PS_LOG_CORRK(20), WGTS_CORRK(8) 
    !         REAL :: WNO_EDGES(NWNO+1), WNO_CTRS(NWNO), STEL_SPEC(NWNO), INT_SPEC(NWNO)
    !         REAL :: OPAC_CIA(73, 20, NWNO), TAURAY_PER_DPG(73, 20, NWNO)

    !         MMW = 8314.462618 / GASCON ! Mean molecular weight in g/mol (or H masses per molecule)

    !         ! Convert METALLICITY and C_TO_O to strings
    !         if (METALLICITY .eq. 0.0) then
    !         dummyMETALLICITY_str = '+000'
    !         else
    !         write(*,*) "Thomas hasn't coded non-solar metallicity k-tables yet"
    !         stop
    !         end if
    !         if (C_TO_O .eq. 0.0) then
    !         dummyC_TO_O_str = '100'
    !         else
    !         write(*,*) "Thomas hasn't coded non-solar C/O k-tables yet"
    !         stop
    !         end if
    !         ! write(*,*) 'with_TiO_and_VO: ', with_TiO_and_VO
    !         if (with_TiO_and_VO .eq. 1.) then
    !         tiovo_str = '_witiovo'
    !         else if (with_TiO_and_VO .eq. 2.) then
    !         tiovo_str = '_notiovo'
    !         else
    !         write(*,*) "tiovo flag not recognized, setting with_TiO_and_VO = 1 for the first loop"
    !         tiovo_str = '_witiovo'
    !         end if
    !         ! write(*,*) 'METALLICITY: ', dummyMETALLICITY_str
    !         ! write(*,*) 'C_TO_O: ', dummyC_TO_O_str
    !         IF (NWNO.eq.11) THEN
    !         FOLDER = '../11bands/'
    !         ELSE IF (NWNO.eq.30) THEN
    !         FOLDER = '../30bands/'
    !         ELSE
    !         WRITE(*,*) 'NWNO is wrong in params.i, options are 11 and 30'
    !         STOP
    !         END IF
    !         FOLDER = '../11bands/'
    !         FNM= trim(FOLDER)//"kcoeffs/feh"//trim(dummyMETALLICITY_str)//"_co_"//trim(dummyC_TO_O_str)//trim(tiovo_str)//".txt"
    !         CIAFNM= trim(FOLDER)//"continuum/cia+hmin_feh"//trim(dummyMETALLICITY_str)//"_co_"//trim(dummyC_TO_O_str)
    !  &              //trim(tiovo_str)//".txt"
    !         RAYFNM= trim(FOLDER)//"rayleigh/rayleigh_feh"//trim(dummyMETALLICITY_str)//"_co_"//trim(dummyC_TO_O_str)
    !  &              //trim(tiovo_str)//".txt"
    !         ! Read in k-coefficients from file
    !         write(*,*) 'Reading in k-coefficients from file: ', FNM
    !         write(*,*) 'Reading in CIA opacities from file: ', CIAFNM
    !         write(*,*) 'Reading in Rayleigh opacities from file: ', RAYFNM
    !         OPEN(UNIT=1, FILE=FNM)
    !         READ(1,*) TEMP_OPAC_CORRK
    !         CLOSE(1)
    !         OPEN(UNIT=1, FILE=CIAFNM)
    !         READ(1,*) TEMP_OPAC_CIA
    !         CLOSE(1)
    !         OPEN(UNIT=1, FILE=RAYFNM)
    !         READ(1,*) TEMP_TAURAY_PER_DPG_PERMU
    !         CLOSE(1)

    !         ! Reshape into 3D array
    !         ! Index order of new array: T, P, wavelength (no gauss pt because gray in each band)
    !         DO I = 1, 73
    !         DO J = 1, 20
    !             DO K = 1, NWNO
    !             ! after this line we have m^2/kg
    !             ! file is in units of m^2/molecule
    !             OPAC_CIA(I,J,K) = TEMP_OPAC_CIA(K,(I-1)*20+J) * 6.022e23 / (MMW/1000)
    !             ! after this line we should have log(tau / (deltap / gravity), all in SI units)
    !             TAURAY_PER_DPG(I,J,K) = log10(TEMP_TAURAY_PER_DPG_PERMU(K,(I-1)*20+J) / MMW)
    !             END DO
    !         END DO
    !         END DO

    !         ! Reshape into 4D array
    !         ! Index order of new array: T, P, wavelength, gauss pt
    !         DO I = 1, 73
    !         DO J = 1, 20
    !             DO K = 1, NWNO
    !             DO L = 1, 8
    !                 ! After exponential, we're in cm^2/molecule units which we convert to m^2/molecule since RT uses SI units
    !                 ! avagadro + MMW/1000 conversion factor get us to m^2/kg
    !                 OPAC_CORRK(I,J,K,L) = EXP(TEMP_OPAC_CORRK(L,(I-1)*20*NWNO+(J-1)*NWNO+K)) * 1e-4 * 6.022e23 / (MMW/1000)
    !                 ! Add in CIA and H- opacities
    !                 OPAC_CORRK(I,J,K,L) = OPAC_CORRK(I,J,K,L) + OPAC_CIA(I,J,K)
    !                 ! Convert to log space
    !                 OPAC_CORRK(I,J,K,L) = LOG10(OPAC_CORRK(I,J,K,L))
    !             END DO
    !             END DO
    !         END DO
    !         END DO





    !         ! Set up the gauss weights of the 8 gauss points
    !         WGTS_CORRK = (/0.16523105, 0.30976895, 0.30976895, 0.16523105, 0.00869637,
    !  &              0.01630363, 0.01630363, 0.00869637/)

    !         ! Read in the temperature and pressure points of our grid
    !         OPEN(UNIT=1, FILE=trim(FOLDER)//"kcoeffs/Tpts.txt")
    !         READ(1,*) TS_CORRK
    !         CLOSE(1)

    !         OPEN(UNIT=1, FILE=trim(FOLDER)//"kcoeffs/Ppts.txt")
    !         READ(1,*) PS_CORRK
    !         CLOSE(1)
    !         PS_CORRK = PS_CORRK * 1.0E5 ! Convert bars to Pascals, which the RT uses

    !         ! Invert T and log P for interpolation
    !         TS_LOG_CORRK = LOG10(TS_CORRK)
    !         PS_LOG_CORRK = LOG10(PS_CORRK)

    !         ! Read in the wavenumber points of our grid
    !         OPEN(UNIT=1, FILE=trim(FOLDER)//"kcoeffs/wno_edges.txt")
    !         READ(1,*) WNO_EDGES
    !         CLOSE(1)

    !         OPEN(UNIT=1, FILE=trim(FOLDER)//"kcoeffs/wno_ctrs.txt")
    !         READ(1,*) WNO_CTRS
    !         CLOSE(1)

    !         ! Read in the stellar spectrum
    !         CALL BIN_STELLAR_SPECTRUM(trim(FOLDER)//"kcoeffs/stellar_spectrum.txt", STEL_SPEC, WNO_EDGES, NWNO)
    !         ! WRITE(*,*) 'STELLAR SPECTRUM: ', STEL_SPEC
    !         ! WRITE(*,*) 'total stellar flux: ', SUM(STEL_SPEC)

    !         ! Grab Planck function integrals for each tmeperature and bin:
    !         OPEN(UNIT=1, FILE=trim(FOLDER)//"kcoeffs/planck_integrals.txt")
    !         READ(1,*) PLANCK_INTS
    !         CLOSE(1)
    !         OPEN(UNIT=1, FILE=trim(FOLDER)//"kcoeffs/planck_temps.txt")
    !         READ(1,*) PLANCK_TS
    !         CLOSE(1)
    !         ! window, temp
    !         ! write(*,*) 'Planck integrals: ', PLANCK_INTS(1,1), PLANCK_INTS(1,2), PLANCK_INTS(2,1)
    !         ! write(*,*) 'Planck temps: ', PLANCK_TS(1), PLANCK_TS(2), PLANCK_TS(3)

    !         ! Now split up internal flux by channel:
    !         TINT = (FBASEFLUX / 5.670367E-8) ** 0.25 ! Convert flux to temperature
    !         TINT_INDEX = MINLOC(ABS(PLANCK_TS - TINT), 1)
    !         INT_SPEC = PLANCK_INTS(:, TINT_INDEX)
    !         INT_SPEC = INT_SPEC / SUM(INT_SPEC)
    !         write(*,*) 'TINT: ', TINT, 'TINT_INDEX: ', TINT_INDEX
    !         write(*,*) 'INT_SPEC: ', INT_SPEC
    !     END SUBROUTINE get_gas_opacity_corrk

    !     SUBROUTINE BIN_STELLAR_SPECTRUM(file_name, STEL_SPEC, WNO_EDGES, NWNO)
    !         implicit none
    !         ! include 'nwno.inc'
    !         integer :: NWNO
    !         CHARACTER(len=80) :: file_name
    !         REAL :: STEL_SPEC(NWNO), STEL_SPEC_IN(2, 5000)
    !         REAL :: WNO_EDGES(NWNO+1)
    !         INTEGER :: I, J
    !         REAL :: WNO(5000), FLUX(5000)
    !         ! Read in the stellar spectrum
    !         OPEN(UNIT=1, FILE=file_name)
    !         READ(1,*) STEL_SPEC_IN
    !         CLOSE(1)
    !         DO I = 1, 5000
    !         WNO(I) = STEL_SPEC_IN(1,I)
    !         FLUX(I) = STEL_SPEC_IN(2,I)
    !         END DO
    !         ! Trapezoidal integration in each spectral band
    !         DO J = 1, NWNO
    !         STEL_SPEC(J) = 0.0
    !         DO I = 1, 4999
    !             IF (WNO(I) >= WNO_EDGES(J) .AND. WNO(I) < WNO_EDGES(J+1)) THEN
    !             STEL_SPEC(J) = STEL_SPEC(J) + ((FLUX(I) + FLUX(I+1))/2 *((WNO(I+1) - WNO(I)) * 3e10))
    !             END IF
    !         END DO
    !         END DO

    !         ! Normalize the spectrum so the fort.7's SOLC_IN will define the total flux:
    !         STEL_SPEC = STEL_SPEC / SUM(STEL_SPEC)

    !     END SUBROUTINE BIN_STELLAR_SPECTRUM

    !     END MODULE corrkmodules