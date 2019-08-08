!  MATTHEW BILSKIE
!  SUMMER 2013
!  COPYRIGHT 2013
!
!  GRD2DEM.F90
!
!  VERSION 5.6 12/11/2015
!
!  VERSION 5.1 10/17/2013
!     - ACCEPTS MAXELE.63 INPUTS AFTER THE GRID FILE
!     - IF MAXELE.63 IS PRESENT THE CODE WILL USE THOSE WSE FOR INTERPOLATION
!     - PROGRAM NOW DETERMINES THE SIZE OF ALLOCATION NEEDED FOR THE
!       RASTER AND PROMPTS THE USER IF OVER 1 GB.
!     - NOW WRITES OUT .prj FILE BASED ON THE UTM_ZONE SPECIFIED
!     - NOW WRITES OUT .stx FILE WITH RASTER STATISTICS (MIN,MAX,MEAN,SD)
!  VERSION 5.2
!     - MINOR BUG FIXES
!  VERSION 5.3
!     - FIXED SOME BUGS IN READING THE INPUT FILE WHEN NOT INCLUDING A
!       MAXELE FILE - CREATED LOGICAL MAXFILE
!     - X,Y MIN MAX IN INPUT FILE ARE NOW SPECIFIED AS GEOGRAPHIC LAT/LONG
!     - 01-28-2014: FIXED BUG IN READING FORT.14,MAXELE.63
!  VERSION 5.4
!     - SWITCHED OUT THE CONVERSION FROM LAT/LONG TO UTM SUBROUTINE
!  VERSION 5.5
!     - ONLY INTERPOLATE CELLS THAT ARE WITHIN A FULLY WETTED ELEMENT.
!  VERSION 5.6 12/11/15
!     - ADDED OPTION TO COMPUTE INUNDATION DEPTH
!     - THIS IS INCORPORATED INTO THE INPUT FILE
!     - REMOVE LEGACY LATLON2UTM SUBROUTINE THAT WAS NO LONGER USED
!
!  THIS CODE READS IN AN ADCIRC MESH AND INPUT FILE AND CONVERTS THE
!  ADCIRC MESH TO A FLOAT/GRID FORMAT DEM
!
!  THIS CODE WAS CREATED TO RE-ENACT ADCDEM BY DAVID COGGIN.
!  THANKS TO DAVID COGGIN FOR HELP IN OPTIMIZING THE ALGORITHM
!
!  VERSION 5 REMOVES THE KDTREE2 SEARCH ALGORITHM
!  VERSION 5 SEARCHES THE ELEMENT TABLE RATHER THAN GRID CELL BY GRID
!  CELL, DRASTICALLY REDUCING THE NUMBER OF ITERATIONS
!
!  FOR MAXIMIMUM EFFICIENTY, COMPILE IN 64-BIT
!  COMPATABILITY w/ -O3 optimization, USING CMD:
!  x86_64-w64-mingw32-gfortran -O3 -m64 GRD2DEMv05.F90
!
!  ----------------------------------------------------------------------------
!
!  COMPILES WITH GFORTRAN
!
!  ----------------------------------------------------------------------------
!
!  THE INPUT FILE SHOULD CONTAIN THE FOLLOWING
!
!  ADCIRC MESH FILE (RENUMBERED & GEOGRAPHIC)
!  Inundation Depth Flag (0 = off, 1 = compute inundation depth)
!  OUTPUT DEM FILE NAME (NO FILE EXTENSION)
!  GRID CELL RESOLUTION
!  href (0 = mesh is in geo, 1 = mesh is in UTM)
!  Z MULTIPLICATION FACTOR
!  XMIN
!  YMIN
!  XMAX
!  YMAX
!
!  ----------------------------------------------------------------------------

      PROGRAM GRD2DEM

      IMPLICIT NONE

      INTEGER                             :: i,j,k
      INTEGER                             :: NE, NP
      INTEGER                             :: NHY = 3
      INTEGER                             :: NUM_ROWS, NUM_COLS
      INTEGER                             :: ROW, COL
      INTEGER                             :: NUM_BB_ROWS, NUM_BB_COLS
      INTEGER                             :: min_row, max_row
      INTEGER                             :: min_col, max_col
      INTEGER                             :: CURR_ROW, CURR_COL
      INTEGER                             :: counter
      INTEGER                             :: UTM_Zone
      INTEGER                             :: HDATUM
      INTEGER,PARAMETER                   :: BAND = 1
      INTEGER                             :: href
      INTEGER                             :: dumI,tempI
      INTEGER                             :: ID_Flag
      INTEGER                             :: RasterFormat

      INTEGER, ALLOCATABLE                :: NID(:),EID(:)
      INTEGER, ALLOCATABLE                :: NM(:,:)

      CHARACTER(60)                       :: gridfile
      CHARACTER(60)                       :: demfile
      CHARACTER(60)                       :: infile,ofile
      CHARACTER(24)                       :: AGRID
      CHARACTER(100)                      :: dumC

      REAL(8)                              :: gridsize
      REAL(8)                              :: mult_fac
      REAL(8)                              :: xMin, xMax
      REAL(8)                              :: yMin, yMax
      REAL(8)                              :: xMin2, xMax2
      REAL(8)                              :: yMin2, yMax2
      REAL(8)                              :: x_diff, y_diff
      REAL(8)                              :: elem_xMin, elem_xMax
      REAL(8)                              :: elem_yMin, elem_yMax
      REAL(8)                              :: elem_x(3), elem_y(3)
      REAL(8)                              :: x_centroid, y_centroid
      REAL(8)                              :: CURR_x, CURR_y
      REAL(8)                              :: x1,x2,x3
      REAL(8)                              :: y1,y2,y3
      REAL(8)                              :: z1,z2,z3
      REAL(8)                              :: a1,a2,a3
      REAL(8)                              :: b1,b2,b3
      REAL(8)                              :: elem_area
      REAL(8)                              :: phi(3)
      REAL(8)                              :: xNodeMin, yNodeMin
      REAL(8)                              :: xNodeMax, yNodeMax
      REAL(8)                              :: finalVal
      REAL(8)                              :: fileSize
      REAL(8)                              :: rstx_min,rstx_max
      REAL(8)                              :: rstx_mean,rstx_sd
      REAL(8),PARAMETER                    :: DRY = -99999.0
      REAL                                 :: version 
      REAL(8)                              :: tempR
      
      REAL(8), ALLOCATABLE                 :: x(:)
      REAL(8), ALLOCATABLE                 :: y(:)
      REAL(8), ALLOCATABLE                 :: z(:)
      REAL(8),ALLOCATABLE                 :: UTMX(:),UTMY(:)
      REAL(8),ALLOCATABLE                 :: WSE(:)

      REAL(4), ALLOCATABLE                :: Raster_Z(:,:)

      LOGICAL                             :: found,maxfile

      version = 5.5

      WRITE(*,*)""
      WRITE(*,'(A)',ADVANCE='YES')"-------------------------------------------"
      WRITE(*,'(A)',ADVANCE='YES')"          GRD2DEM VERSION 5.5              "
      WRITE(*,'(A)',ADVANCE='YES')"            MATTHEW BILSKIE                "
      WRITE(*,'(A)',ADVANCE='YES')"            CHAMPS Lab UCF                 "
      WRITE(*,'(A)',ADVANCE='YES')"            COPYRIGHT 2013                 "
      WRITE(*,'(A)',ADVANCE='YES')"-------------------------------------------"
      WRITE(*,*)""

      IF (IARGC().NE.1) THEN
          WRITE(*,'(A)',ADVANCE='YES')"NAME OF INPUT FILE: "
          WRITE(*,'(A)',ADVANCE='NO')"==> "
          READ(*,*)infile
          WRITE(*,*)""
          !WRITE(*,'(A)',ADVANCE='YES')"USAGE: *.exe <stdin>"
          !WRITE(*,'(A)',ADVANCE='NO')"PRESS ANY KEY TO EXIT"
          !READ(*,*)
          !STOP
      ELSE
          CALL GETARG(1,infile)
      ENDIF

      !CALL GETARG(1,infile)
      INQUIRE(FILE = TRIM(infile),EXIST=found) ! CHECKS TO SEE IF INPUT FILE EXISTS
      IF (found) THEN
          WRITE(*,1010)infile
          OPEN(10,FILE=infile,STATUS='OLD',ACTION='READ')
      ELSE
          WRITE(*,1011)infile
          WRITE(*,'(A)',ADVANCE='YES')"PRESS ANY KEY TO EXIT"
          READ(*,*)
          STOP
      ENDIF

      maxfile = .TRUE.
      READ(10,'(A70)')dumC
      loopA: DO i=1,LEN_TRIM(dumC)
          IF(dumC(i:i).EQ.",")THEN
              gridfile=(dumC(1:i-1))
              ofile = TRIM(dumC(i+1:LEN_TRIM(dumC)))
          ENDIF
      ENDDO loopA
      !IF ((i-1).eq.LEN_TRIM(dumC)) THEN
      IF ((i-1).LE.LEN_TRIM(gridfile)) THEN
          gridfile = dumC
          ofile = ""
          maxfile = .FALSE.
      ENDIF
      
      READ(10,*)ID_Flag
      READ(10,*)demfile
      READ(10,*)gridsize
      REAd(10,*)href
      READ(10,*)mult_fac
      READ(10,*)xMin2
      READ(10,*)yMin2
      READ(10,*)xMax2
      READ(10,*)yMax2
      READ(10,*)UTM_Zone
      READ(10,*)HDATUM
      READ(10,*)RasterFormat
      CLOSE(10)
      IF(xMin2.GT.xMax2)THEN
          TempR = xMin2
          xMin2 = xMax2
          xMax2 = TempR
      ENDIF
      IF(yMin2.GT.yMax2)THEN
          TempR = yMin2
          yMin2 = yMax2
          yMax2 = TempR
      ENDIF
      IF(href.EQ.1)THEN
          CALL LatLon2UTM(xMin2,yMin2,UTM_Zone,HDATUM,xMin,yMin)
          CALL LatLon2UTM(xMax2,yMax2,UTM_Zone,HDATUM,xMax,yMax)
      ELSE
          xMin = xMin2
          yMin = yMin2
          xMax = xMax2
          yMax = yMax2
      ENDIF

      INQUIRE(FILE = TRIM(gridfile),EXIST=found) ! CHECKS TO SEE IF INPUT FILE EXISTS
      IF (found) THEN
          WRITE(*,1010)gridfile
          OPEN(14,FILE=gridfile,STATUS='OLD',ACTION='READ')
      ELSE
          WRITE(*,1011)gridfile
          WRITE(*,'(A)',ADVANCE='YES')"PRESS ANY KEY TO EXIT"
          READ(*,*)
          STOP
      ENDIF

      IF (maxfile) THEN
          INQUIRE(FILE = TRIM(ofile),EXIST=found) ! CHECKS TO SEE IF INPUT FILE EXISTS
          IF (found) THEN
              WRITE(*,1010)ofile
              OPEN(63,FILE=ofile,STATUS='OLD',ACTION='READ')
          ELSE
              WRITE(*,1011)ofile
              WRITE(*,'(A)',ADVANCE='YES')"PRESS ANY KEY TO EXIT"
              READ(*,*)
              STOP
          ENDIF
      ENDIF

      READ(14,*)AGRID
      READ(14,*)NE,NP
      IF (ofile.NE."") THEN
          READ(63,*)
          READ(63,*)tempI,dumI
          IF (dumI.ne.NP) THEN
              WRITE(*,*)'GEOMETRY IN ',TRIM(ofile),' DOES NOT MATCH GRID FILE.'
              WRITE(*,*)
              STOP
          ENDIF
          READ(63,*)
      ENDIF

      ALLOCATE ( x(NP) ) 
      ALLOCATE ( y(NP) ) 
      ALLOCATE ( z(NP) ) 
      ALLOCATE ( NID(NP) )
      ALLOCATE ( EID(NE) )
      ALLOCATE ( NM(NE,3) )
      ALLOCATE ( WSE(NP) )
      ALLOCATE (UTMX(NP) )
      ALLOCATE (UTMY(NP) )
      
      WRITE(*,'(A)',ADVANCE='YES')"READING NODE TABLE..."
      DO i = 1, NP

          READ(14,*)NID(i),x(i),y(i),z(i)
          IF(ID_Flag.EQ.1) z(i) = -1*z(i)
          IF(maxfile) READ(63,*)dumI,WSE(i)

          IF(href.EQ.1) THEN
            CALL LatLon2UTM(x(i),y(i),UTM_Zone,HDATUM,UTMX(i),UTMY(i))
          ELSE
            UTMX(i) = x(i)
            UTMY(i) = y(i)
          ENDIF

          IF (ID_Flag.EQ.1) THEN
            IF((WSE(i).NE.-99999).AND.(z(i).GT.0)) THEN
                z(i) = WSE(i) - z(i)
            ELSE
                z(i) = -99999.0D0
            ENDIF
            IF(z(i).LE.0) z(i) = -99999.0D0
          ELSEIF (maxfile) THEN
            IF(WSE(i).GT.(-1*z(i)))THEN
                z(i) = WSE(i)
            ELSE
                z(i) = -99999.0D0
            ENDIF
          ENDIF

          IF (NID(i).ne.i) THEN
              WRITE(*,'(A)',ADVANCE='YES')"NODE SHOULD BE RENUMBERED"
              WRITE(*,'(A)',ADVANCE='NO')"PRESS ANY KEY TO EXIT"
              READ(*,*)
          ENDIF

          IF (UTMx(i).LT.xNodeMin) xNodeMin = UTMX(i)
          IF (UTMY(i).LT.yNodeMin) yNodeMin = UTMY(i)
          IF (UTMX(i).GT.xNodeMax) xNodeMax = UTMX(i)
          IF (UTMY(i).GT.yNodeMax) yNodeMax = UTMY(i)

      ENDDO

      WRITE(*,'(A)',ADVANCE='YES')"READING ELEMENT TABLE..."
      DO i = 1, NE

          READ(14,*)EID(i),dumC,NM(i,1),NM(i,2),NM(i,3)

      ENDDO
      WRITE(*,'(A)')""
      CLOSE(14)

      ! DETERMINE IF MESH NODE XNODEMIN AND YNODEMIN ARE INSIDE OR
      ! OUTSIDE THE DESIRTED RASTER BBOX. THIS WILL REMOVE UNNECESSARY
      ! COLUMNS AND ROWS OF NODATA VALUES.  IT WILL ALSO START THE
      ! LLCORNER AT THE LLC OF THE MESH.
      IF (xNodeMin.GE.xMin) xMin = xNodeMin
      IF (yNodeMin.GE.yMin) yMin = yNodeMin
      IF (xNodeMax.LE.xMax) xMax = xNodeMax
      IF (yNodeMax.LE.yMax) yMax = yNodeMax

      ! COMPUTE INFO ABOUT OUTPUT RASTER (NUM. ROWS & COLS)
      NUM_ROWS = CEILING((yMax - yMin) / gridsize)
      NUM_COLS = CEILING((xMax - xMin) / gridsize)
      ! REFERENCE GRID FROM LOWER LEFT CORNER
      ! RECOMPUTE yMax and xMax BASED ON gridsize
      xMax = xMin + (NUM_COLS*gridsize)
      yMax = yMin + (NUM_ROWS*gridsize)

      fileSize = NUM_COLS*NUM_ROWS
      fileSize = fileSize * 4 ! BYTES (EACH CELL IS 32 BITES = 4 BYTES)
      fileSize = fileSize * (9.5367E-7) ! bytes to megabytes
      IF (fileSize.GT.1024.0D0) THEN
          WRITE(*,301)fileSize
          READ(*,*)dumC
      ENDIF
      IF ((dumC=='N').OR.(dumC=='n')) STOP
      ALLOCATE ( Raster_Z(NUM_COLS,NUM_ROWS) )
      DO i = 1, NUM_ROWS
          DO j = 1, NUM_COLS
              Raster_Z(j,i) = -99999
          ENDDO
      ENDDO

      WRITE(*,302)fileSize,fileSize/1024.0D0
      
      !WRITE(*,'(I12,A)',ADVANCE='YES')NUM_ROWS," NUMBER OF ROWS CREATED."
      !WRITE(*,'(I12,A)',ADVANCE='YES')NUM_COLS," NUMBER OF COLUMNS CREATED."
      !WRITE(*,*)

      WRITE(*,'(A)',ADVANCE='YES')"INTERPOLATING..."

      rstx_min = 99999
      rstx_max = -99999
      rstx_mean = 0.0d0
      rstx_sd = 0.0d0
      counter = 0
      DO i = 1, NE

          !WRITE(*,*)"ELEMENT ",i, NM(i,1),NM(i,2),NM(i,3)

          ! DETERMINE IF ELEMENT i IS INSIDE OUTPUT RASTER BOUNDS
          ! FIRST, FIND THE CENTROID OF THE ELEMENT
          x_centroid = ( (UTMX(NM(i,1)) + UTMX(NM(i,2)) + UTMX(NM(i,3)) ) / 3 )
          y_centroid = ( (UTMY(NM(i,1)) + UTMY(NM(i,2)) + UTMY(NM(i,3)) ) / 3 )

          IF ( (x_centroid.GE.xMin).and.(x_centroid.LE.xMax).and. &
             (y_centroid.GE.yMin).and.(y_centroid.LE.yMax) ) THEN
              
              ! ELEMENT i IS INSIDE THE RASTER BOUNDS

              ! COMPUTE THE BOUNDING BOX OF ELEMENT i
              DO j = 1, 3
                  elem_x(j) = UTMX( NM(i,j) ) 
                  elem_y(j) = UTMY( NM(i,j) ) 
              ENDDO
              elem_xMin = MINVAL(elem_x)
              elem_xMax = MAXVAL(elem_x)
              elem_yMin = MINVAL(elem_y)
              elem_yMax = MAXVAL(elem_y)

              ! FIND THE RASTER CELLS THAT ARE INSIDE THIS BBOX
              ! (REFERENCE UPPER LEFT CORNER)
              ! FIRST, FIND ROW,COL OF CENTER OF ELEMENT
              ROW = CEILING((yMax - y_centroid) / gridsize)
              COL = CEILING((x_centroid - xMin) / gridsize)

              ! FIND THE CELL OF ELEMENT i's BBOX (xMin,yMin & xMax,yMax)
              min_row = CEILING((yMax - elem_yMin) / gridsize)
              max_row = CEILING((yMax - elem_yMax) / gridsize)
              min_col = CEILING((elem_xMin - xMin) / gridsize)
              max_col = CEILING((elem_xMax - xMin) / gridsize)

              ! COMPUTER THE NUMBER OF GRID CELLS W/IN ELEMENT i's BBOX
              NUM_BB_ROWS = min_row - max_row + 1
              NUM_BB_COLS = max_col - min_col + 1

              ! SEARCH THESE GRID CELLS AND DETERMINE IF THE CELLS
              ! CENTER IS IN ELEMENT i
              DO j = 1, NUM_BB_ROWS

                  DO k = 1, NUM_BB_COLS

                      ! DETERMINE ROW & COLUMN
                      CURR_ROW = max_row + (j-1)
                      CURR_COL = min_col + (k-1)
                      !WRITE(*,*)CURR_ROW,CURR_COL
                      IF (CURR_ROW.le.0) CURR_ROW = 1
                      IF (CURR_COL.le.0) CURR_COL = 1

                      ! DETERMINE x,y of CELL CENTER
                      CURR_y = yMax - (CURR_ROW * gridsize) + (gridsize/2)
                      CURR_x = xMin + (CURR_COL * gridsize) - (gridsize/2)

                      ! COMPUTE THE LINEAR BASIS FUNCTIONS
                      x1 = UTMX(NM(i,1))
                      x2 = UTMX(NM(i,2))
                      x3 = UTMX(NM(i,3))
                      y1 = UTMY(NM(i,1))
                      y2 = UTMY(NM(i,2))
                      y3 = UTMY(NM(i,3))
                      z1 = z(NM(i,1))
                      z2 = z(NM(i,2))
                      z3 = z(NM(i,3))

                      a1 = x3 - x2
                      a2 = x1 - x3
                      a3 = x2 - x1
                      b1 = y2 - y3
                      b2 = y3 - y1
                      b3 = y1 - y2

                      elem_area = (b1*a2 - b2*a1) / 2

                      phi(1) = ((x2*y3)-(x3*y2)+(b1*CURR_x)+(a1*CURR_y)) / (2*elem_area)
                      phi(2) = ((x3*y1)-(x1*y3)+(b2*CURR_x)+(a2*CURR_y)) / (2*elem_area)
                      phi(3) = ((x1*y2)-(x2*y1)+(b3*CURR_x)+(a3*CURR_y)) / (2*elem_area)
                      !WRITE(*,*) phi(1),phi(2),phi(3)

                      IF ( (MINVAL(phi).GE.0d0).AND.(MAXVAL(phi).LE.1d0) ) THEN
                          ! RASTER CELL CENTER IS INSIDE THE ELEMENT
                          IF ( (CURR_COL.LE.NUM_COLS).AND.(CURR_ROW.LE.NUM_ROWS) ) THEN
                              IF( (z1.NE.-99999).AND.(z2.NE.-99999).AND.(z3.NE.-99999) )THEN
                                  finalVal = mult_fac * ((z1*phi(1)) + (z2*phi(2)) + (z3*phi(3)))
                              ELSE
                                  finalVal = -99999.0D0
                              ENDIF
                              !IF ((LEN_TRIM(ofile).GT.0).AND.(finalVal.LE.-10.0D0)) finalVal = -99999
                              Raster_Z(CURR_COL,CURR_ROW) = finalVal
                              IF (finalVal.NE.-99999) THEN
                                  IF (finalVal.lt.rstx_min)rstx_min=finalVal
                                  IF (finalVal.gt.rstx_max)rstx_max=finalVal
                                  rstx_mean = rstx_mean + finalVal
                                  counter = counter + 1
                              ENDIF
                          ENDIF
                          !WRITE(100,*)CURR_x,CURR_y
                      ENDIF

                  ENDDO

              ENDDO

          ENDIF

      ENDDO

      WRITE(*,'(A)',ADVANCE='YES')""
      WRITE(*,'(A)',ADVANCE='YES')"FINISHED INTERPOLATING!"
      WRITE(*,'(A)',ADVANCE='YES')""

      WRITE(*,'(A)',ADVANCE='YES')"WRITING GRID..."
      OPEN(UNIT=15,FILE=TRIM(demfile)//'.hdr',ACTION='WRITE',STATUS='UNKNOWN')
      OPEN(UNIT=16,FILE=TRIM(demfile)//'.flt',ACTION='WRITE',  &
            STATUS='UNKNOWN',FORM='UNFORMATTED',ACCESS='STREAM')

      WRITE(15,200)NUM_COLS
      WRITE(15,201)NUM_ROWS
      WRITE(15,202)xMin
      WRITE(15,203)yMin
      WRITE(15,204)gridsize
      WRITE(15,205)-99999
      WRITE(15,206)
      WRITE(15,207)
      WRITE(15,208)
200   format('NCOLS',2x,I12)
201   format('NROWS',2x,I12)
202   format('XLLCORNER',2x,f20.5)
203   format('YLLCORNER',2x,f20.5)
204   format('CELLSIZE',2x,f15.5)
205   format('NODATA_VALUE',2x,I6)
206   format('BYTEORDER LSBFIRST')
207   format('NBITS 32')
208   format('PIXELTYPE FLOAT')
300   format('RASTER FILE SIZE:',x,f10.6,x,'mb.')
301   format('YOU ARE ABOUT TO CREATE A RASTER THAT IS',x,f15.6,x,'MB. DO YOU WANT TO CONTINUE (Y/N)?')
302   format('ALLOCATING RASTER ARRAY ->',x,f15.6,x,'MB in file size (',x,f15.6,x,'GB).')

      DO i = 1, NUM_ROWS
          WRITE(16)(Raster_Z(j,i),j=1,NUM_COLS)
      ENDDO

      CLOSE(15)
      CLOSE(16)
    
      ! WRITE THE .PRJ FILE
      tempR = UTM_Zone*6.0D0-183.0D0
      OPEN(UNIT=15,FILE=TRIM(demfile)//'.prj',ACTION='WRITE',STATUS='UNKNOWN')
      WRITE(15,'(A,I2,A)')'PROJCS["NAD_1983_UTM_Zone_',UTM_Zone,'N",'
      WRITE(15,*)'GEOGCS["GCS_North_American_1983",'
      WRITE(15,*)'DATUM["D_North_American_1983",'
      WRITE(15,*)'SPHEROID["GRS_1980",6378137,298.257222101]],'
      WRITE(15,*)'PRIMEM["Greenwich",0],'
      WRITE(15,*)'UNIT["Degree",0.0174532925199433]],'
      WRITE(15,*)'PROJECTION["Transverse_Mercator"],'
      WRITE(15,*)'PARAMETER["False_Easting",500000.0],'
      WRITE(15,*)'PARAMETER["False_Northing",0.0],'
      WRITE(15,*)'PARAMETER["Central_Meridian",',tempR,'],'
      WRITE(15,*)'PARAMETER["Scale_Factor",0.9996],'
      WRITE(15,*)'PARAMETER["Latitude_of_Origin",0.0],'
      WRITE(15,*)'UNIT["Meter",1.0]]'
      CLOSE(15)

      ! WRITE THE .STX (STATISTICS FILE)
      ! FIRST, LETS SETUP/COMPUTE SOME OF THE STATS
      rstx_min = rstx_min
      rstx_max = rstx_max
      rstx_mean = rstx_mean / counter
      DO i = 1, NUM_ROWS
          DO j = 1, NUM_COLS
              IF (Raster_Z(j,i).NE.-99999) THEN
                  rstx_sd = rstx_sd + ((Raster_Z(j,i) - rstx_mean)**2.0d0)
              ENDIF
          ENDDO
      ENDDO
      DEALLOCATE(Raster_Z)
      rstx_sd = SQRT(rstx_sd/(counter-1))
      OPEN(UNIT=15,FILE=TRIM(demfile)//'.stx',ACTION='WRITE', &
         STATUS='UNKNOWN')
      WRITE(15,'(I1,4(x,1pE20.10E3))')BAND,rstx_min,rstx_max,rstx_mean,rstx_sd
      CLOSE(15)

      IF(RasterFormat.EQ.1)THEN !... Convert to *.img using gdal tools
          dumC = 'gdal_translate -of HFA ' // TRIM(demfile) // '.flt ' &
                  // TRIM(demfile) // '.img'
          CALL SYSTEM(TRIM(dumC))
          dumC = 'rm '//TRIM(demfile)
          CALL SYSTEM(TRIM(dumC)//'.flt')
          CALL SYSTEM(TRIM(dumC)//'.hdr')
          CALL SYSTEM(TRIM(dumC)//'.prj')
          CALL SYSTEM(TRIM(dumC)//'.stx')
      ENDIF


      WRITE(*,*)''
      WRITE(*,'(A)',ADVANCE='YES')"FINISHED!"
      WRITE(*,*)

1010  format(' FILE ',a60,/,'WAS FOUND!   OPENING FILE',/)      
1011  format(' FILE ',a60,/,'WAS NOT FOUND!  PROGRAM TERMINATING',/)

      DEALLOCATE (NID)
      DEALLOCATE (EID)
      DEALLOCATE (NM)
      DEALLOCATE (x)
      DEALLOCATE (y)
      DEALLOCATE (z)

      END PROGRAM


    SUBROUTINE LatLon2UTM(MyX,MyY,MyZone,MyHorizontal,UTM_X,UTM_Y)

        IMPLICIT NONE

        INTRINSIC :: DACOS

        !..IN/OUT Variables
        INTEGER,INTENT(IN)  :: MyHorizontal
        INTEGER,INTENT(IN)  :: MyZone
        REAL(8),INTENT(IN)  :: MyX
        REAL(8),INTENT(IN)  :: MyY
        REAL(8),INTENT(OUT) :: UTM_X
        REAL(8),INTENT(OUT) :: UTM_Y

        !...Local Variables
        REAL(8) :: DLat
        REAL(8) :: DLon
        REAL(8) :: RLat
        REAL(8) :: RLon
        REAL(8) :: X
        REAL(8) :: Y
        REAL(8),PARAMETER :: UTMScaleFactor = 0.9996d0
        REAL(8) :: A
        REAL(8) :: B
        REAL(8) :: PI
        REAL(8) :: DN
        REAL(8) :: SALPHA
        REAL(8) :: SBETA
        REAL(8) :: SGAMMA
        REAL(8) :: SDELTA
        REAL(8) :: SEPSILON
        REAL(8) :: SLENGTH
        REAL(8) :: CMERIDIAN
        REAL(8) :: SEP2
        REAL(8) :: SNU2
        REAL(8) :: SN
        REAL(8) :: T
        REAL(8) :: T2
        REAL(8) :: TMP
        REAL(8) :: S1
        REAL(8) :: SL
        REAL(8) :: SL3COEF
        REAL(8) :: SL4COEF
        REAL(8) :: SL5COEF
        REAL(8) :: SL6COEF
        REAL(8) :: SL7COEF
        REAL(8) :: SL8COEF
        INTEGER :: I

        SELECT CASE(MyHorizontal)
            CASE(1)
               A = 6378137.d0      ! Equatorial Radius
               B = 6356752.3141d0  ! Polar Radius
            CASE(2)
               A = 6378137.d0      ! Equatorial Radius
               B = 6356752.3142d0  ! Polar Radius
            CASE(3)
               A = 6378135.d0      ! Equatorial Radius
               B = 6356750.5000d0  ! Polar Radius
           CASE DEFAULT
               WRITE(*,'(A)') "Invalid Horizontal System."
               STOP
        END SELECT

        PI = DACOS(-1.0D0)

        DLON = MyX
        DLAT = MyY
        RLAT = DLAT * PI / 180.0D0
        RLON = DLON * PI / 180.0D0
        DN = (A-B) / (A+B)
        SALPHA = ((A+B)/2.D0) * ( 1.D0 + DN**(2)/4.D0 + DN**(4)/ 64.D0 )
        SBETA = ( -3.D0*DN/2.D0 ) + ( 9.D0*DN**(3)/16.D0 ) + (-3.D0*DN**(5)/32.D0)
        SGAMMA = ( 15.D0*DN**2)/16.D0 - (15.D0*DN**(4))/32.D0
        SDELTA = (-35.D0*DN**(3))/48.D0 + (105.D0*DN**(5))/256.D0
        SEPSILON = 315.D0 * DN ** (4) / 512.D0
        SLENGTH = SALPHA * ( RLAT + SBETA  * DSIN(2.D0*RLAT) + SGAMMA * DSIN(4.D0*RLAT) &
                           + SDELTA * DSIN(6.D0*RLAT) + SEPSILON * DSIN(8.D0*RLAT) )
        CMERIDIAN = ( -183.D0 + MyZone*6.D0 ) * PI / 180.D0

        SEP2 = (A**(2) - B **(2)) / (B**(2))
        SNU2 = SEP2 * (DCOS(RLAT) ** (2))
        SN   = A*A / ( B*DSQRT(1.D0+SNU2) )
        T    = DTAN(RLAT)
        T2   = T * T
        TMP  = ( T2*T2*T2 ) - T**(6)
        SL   = RLON - CMERIDIAN
        SL3COEF = 1.D0 - T2 + SNU2
        SL4COEF = 5.D0 - T2 + 9.D0*SNU2 + 4.D0*SNU2*SNU2
        SL5COEF = 5.D0 - 18.D0*T2 + T2*T2 + 14.D0*SNU2 - 58.D0*  T2*SNU2
        SL6COEF = 61.D0 - 58.D0*T2 + T2*T2 + 270.D0*SNU2 - 330.D0*  T2*SNU2
        SL7COEF = 61.D0 - 479.D0*T2 + 179.D0*T2*T2 - T2*T2*T2
        SL8COEF = 1385.D0 - 3311.D0*T2 + 543.D0*T2*T2 - T2*T2*T2
        X = SN * DCOS(RLAT)                * SL                &
            + SN * DCOS(RLAT)**(3) * SL3COEF * SL**(3) /    6.D0 &
            + SN * DCOS(RLAT)**(5) * SL5COEF * SL**(5) /  120.D0 &
            + SN * DCOS(RLAT)**(7) * SL7COEF * SL**(7) / 5040.D0
        Y = SLENGTH &
            + T * SN * DCOS(RLAT)**(2)           * SL**(2) /     2.D0 &
            + T * SN * DCOS(RLAT)**(4) * SL4COEF * SL**(4) /    24.D0 &
            + T * SN * DCOS(RLAT)**(6) * SL6COEF * SL**(6) /   720.D0 &
            + T * SN * DCOS(RLAT)**(8) * SL8COEF * SL**(8) / 40320.D0

        x = x * UTMScaleFactor + 500000.d0
        y = y * UTMScaleFactor

        IF(Y.LT.0.0D0)THEN
            Y = Y + 10000000.D0
        ENDIF

        UTM_X = X
        UTM_Y = Y

        RETURN

    END SUBROUTINE
