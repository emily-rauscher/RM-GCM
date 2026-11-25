C**********************************************************
C             MODULE DATAPATHS
C**********************************************************
      MODULE datapaths
C     Where the model's input data live: DATA_DIR holds CLOUD_DATA/ and
C     the corr-k band sets (11bands/, 30bands/, ...). Set by DATA_DIR in
C     fort.7 (&INVARPARAM), absolute or relative to the run directory;
C     the default, '../', is the repository when the run directory is
C     one level below it. INIVARPARAM adds a trailing '/' if missing.
      SAVE
      CHARACTER(len=256) :: DATA_DIR = '../'
      END MODULE datapaths
