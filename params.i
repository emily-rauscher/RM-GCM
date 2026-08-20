! Picket-fence or double-gray
!C     T2 L50  full sphere
!      PARAMETER(NN=2,MM=2,NHEM=2, NL=50,MOCT=1,MG=8,JG=2,NWJ2=2
!     +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1,NIR=2,NSOL=3,NTOTAL=5,
!     +NKGAUSS=3,NBATCH=5)

!C     T31 L50  full sphere
!     NTRAC = 1 (water vapour) + one slot per settling dye tracer, so
!     NTRAC=4 gives three dyes. The dyes are identical in every respect
!     except their grain radius ADYE(KK), set per tracer in fort.7.
!     Raising NTRAC here is all that is needed to add another dye:
!     every dye loop below runs KK=2,NTRAC. Remember to add matching
!     ADYE and KOLOUR entries in fort.7.
      PARAMETER(NN=31,MM=31,NHEM=2, NL=50,MOCT=1,MG=96,JG=24,NWJ2=256
     +,NCRAY=8,JGL=JG,NTRAC=4,NLEVRF=1,NIR=2,NSOL=3,NTOTAL=5,
     +NKGAUSS=3,NBATCH=5)

! 11 bin correlated-k
!C     T2 L50  full sphere
!      PARAMETER(NN=2,MM=2,NHEM=2, NL=50,MOCT=1,MG=8,JG=2,NWJ2=2
!     +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1,NIR=88,NSOL=88,NTOTAL=176,
!     +NKGAUSS=8,NBATCH=16)

!C     T31 L50  full sphere
!      PARAMETER(NN=31,MM=31,NHEM=2, NL=50,MOCT=1,MG=96,JG=24,NWJ2=256
!     +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1,NIR=88,NSOL=88,NTOTAL=176,
!     +NKGAUSS=8,NBATCH=16)

! 30 bin correlated-k
!C     T2 L50  full sphere
!      PARAMETER(NN=2,MM=2,NHEM=2, NL=50,MOCT=1,MG=8,JG=2,NWJ2=2
!     +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1,NIR=240,NSOL=240,NTOTAL=480,
!     +NKGAUSS=8,NBATCH=16)

!C     T31 L50  full sphere
!      PARAMETER(NN=31,MM=31,NHEM=2, NL=50,MOCT=1,MG=96,JG=24,NWJ2=256
!     +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1,NIR=240,NSOL=240,NTOTAL=480,
!     +NKGAUSS=8,NBATCH=16)
