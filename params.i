! Picket-fence or double-gray
!C     T2 L50  full sphere
!      PARAMETER(NN=2,MM=2,NHEM=2, NL=50,MOCT=1,MG=8,JG=2,NWJ2=2
!     +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1,NIR=2,NSOL=3,NTOTAL=5,
!     +NKGAUSS=3,NBATCH=5,NDYE=2)

!C     T31 L50  full sphere
!     NTRAC is the number of SPECTRAL tracers: 1, water vapour only.
!     NDYE is the number of cloud dye tracers, on the Gaussian grid
!     (dyeadv.f), numbered 2..NDYE+1 with a radius ADYE(KK) each in
!     fort.7. KDYEMODE=2 needs NDYE even.
!     KOLOUR in fort.7 takes NTRAC (=1) entries.
      PARAMETER(NN=31,MM=31,NHEM=2, NL=50,MOCT=1,MG=96,JG=24,NWJ2=256
     +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1,NIR=2,NSOL=3,NTOTAL=5,
     +NKGAUSS=3,NBATCH=5,NDYE=2)

! 11 bin correlated-k
!C     T2 L50  full sphere
!      PARAMETER(NN=2,MM=2,NHEM=2, NL=50,MOCT=1,MG=8,JG=2,NWJ2=2
!     +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1,NIR=88,NSOL=88,NTOTAL=176,
!     +NKGAUSS=8,NBATCH=16,NDYE=2)

!C     T31 L50  full sphere
!      PARAMETER(NN=31,MM=31,NHEM=2, NL=50,MOCT=1,MG=96,JG=24,NWJ2=256
!     +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1,NIR=88,NSOL=88,NTOTAL=176,
!     +NKGAUSS=8,NBATCH=16,NDYE=2)

! 30 bin correlated-k
!C     T2 L50  full sphere
!      PARAMETER(NN=2,MM=2,NHEM=2, NL=50,MOCT=1,MG=8,JG=2,NWJ2=2
!     +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1,NIR=240,NSOL=240,NTOTAL=480,
!     +NKGAUSS=8,NBATCH=16,NDYE=2)

!C     T31 L50  full sphere
!      PARAMETER(NN=31,MM=31,NHEM=2, NL=50,MOCT=1,MG=96,JG=24,NWJ2=256
!     +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1,NIR=240,NSOL=240,NTOTAL=480,
!     +NKGAUSS=8,NBATCH=16,NDYE=2)
