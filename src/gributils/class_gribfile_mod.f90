MODULE class_gribfile

    !*****************************************************************
    !                                                                *
    !  This is a class-like structure for providing access to the    *
    !  metadata and data within a GRIB met file                      *
    !                                                                *
    !                                                                *
    !*****************************************************************

    IMPLICIT NONE
    PRIVATE   ! The default is that everything is private, unless
              ! specified otherwise

    ! Note that all public interfaces and variables should have a
    ! GRIBFILE_ prefix

    PUBLIC :: gribfile_object,                &
              gribfile_object_create,         &
              gribfile_printobj,              &
              gribfile_centre,                &
              gribfile_num_xlon,              &
              gribfile_num_ylat,              &
              gribfile_num_zlevel


    PUBLIC :: GRIBFILE_TYPE_ECMWF_GRIB1,      &
              GRIBFILE_TYPE_ECMWF_GRIB2,      &
              GRIBFILE_TYPE_ECMWF_GRIB1_2,    &
              GRIBFILE_TYPE_NCEP_GRIB1,       &
              GRIBFILE_TYPE_NCEP_GRIB2,       &
              GRIBFILE_TYPE_UNKNOWN,          &
              GRIBFILE_CENTRE_NCEP,           &
              GRIBFILE_CENTRE_ECMWF,          &
              GRIBFILE_CENTRE_UNKNOWN

    ! These are codes for designating the type of GRIB file  
    ! being looked at
    INTEGER, PARAMETER :: GRIBFILE_TYPE_ECMWF_GRIB1 = 1,     &
                          GRIBFILE_TYPE_ECMWF_GRIB2 = 2,     &
                          GRIBFILE_TYPE_ECMWF_GRIB1_2 = 3,   &
                          GRIBFILE_TYPE_NCEP_GRIB1 = 4,      &
                          GRIBFILE_TYPE_NCEP_GRIB2 = 5,      &
                          GRIBFILE_TYPE_UNKNOWN = -9999,     &
                          GRIBFILE_CENTRE_NCEP = 1,          &
                          GRIBFILE_CENTRE_ECMWF = 2,         &
                          GRIBFILE_CENTRE_UNKNOWN = -9999

    ! These are the official centre codes for NCEP and ECMWF in grib files.
    INTEGER, PARAMETER :: CENTRE_NCEP = 7, CENTRE_ECMWF = 98

    TYPE gribfile_object
        PRIVATE    ! Make everything in here private so it's not directly manipulated outside
        LOGICAL :: is_instantiated = .FALSE.
        CHARACTER(LEN=256) :: file_path = ''
        INTEGER :: grib_edition = 0  ! Not sure we want this, since it can vary on hybrid files
        INTEGER :: grib_centre = GRIBFILE_CENTRE_UNKNOWN
        INTEGER :: gribfile_type = GRIBFILE_TYPE_UNKNOWN
        INTEGER :: num_xlon = -9999
        INTEGER :: num_ylat = -9999
        INTEGER :: num_zlevel = -9999
    END TYPE gribfile_object



CONTAINS

    SUBROUTINE gribfile_testhello()
        PRINT *, 'Hello gribfile'
    END SUBROUTINE gribfile_testhello



    TYPE(gribfile_object) FUNCTION gribfile_object_create(filepath)

        ! This is the "constructor" for the pseudo gribfile object.  Given the path to a gribfile,
        ! fills in attributes that can be accessed through methods in this
        ! module

        USE grib_api

        IMPLICIT NONE

        CHARACTER(LEN=*), INTENT(IN) :: filepath  ! full path to GRIB file

        TYPE(gribfile_object) :: returned_object

        INTEGER :: ifile, iret, igrib, grib_centre, gribfile_type
        INTEGER :: xlon, ylat, zlev


        CALL get_centre_and_type(filepath, grib_centre, gribfile_type)
        returned_object%grib_centre = grib_centre
        returned_object%gribfile_type = gribfile_type

        ! Get dimensions of 3d u field
        CALL get_3d_u_dims(filepath, gribfile_type, xlon, ylat, zlev)
        returned_object%num_xlon = xlon
        returned_object%num_ylat = ylat
        returned_object%num_zlevel = zlev


        returned_object%is_instantiated = .TRUE.
        returned_object%file_path = TRIM(filepath)

        gribfile_object_create = returned_object

    END FUNCTION gribfile_object_create


    SUBROUTINE gribfile_printobj(gribfile_obj)

        ! Pretty prints the attributes of the gribfile pseudo-object
        TYPE(gribfile_object), INTENT(IN) :: gribfile_obj

        PRINT *, 'is_instantiated: ', gribfile_obj%is_instantiated
        PRINT *, 'filepath: ', TRIM(gribfile_obj%file_path)
        PRINT *, 'grib_centre: ', gribfile_obj%grib_centre
        PRINT *, 'gribfile_type: ', gribfile_obj%gribfile_type
        PRINT *, 'num_xlon: ', gribfile_obj%num_xlon
        PRINT *, 'num_ylat: ', gribfile_obj%num_ylat
        PRINT *, 'num_zlevel: ', gribfile_obj%num_zlevel

    END SUBROUTINE gribfile_printobj

    INTEGER FUNCTION gribfile_centre(filepath)

        ! Returns an integer constant denoting the grib centre (currently either ECMWF, NCEP or UNKNOWN)
        ! for the specified filepath.  Returns one of the GRIBFILE_CENTRE_ constants defined at top of this
        ! module.


        USE grib_api

        IMPLICIT NONE

        CHARACTER(LEN=*), INTENT(IN) :: filepath ! full path to GRIB file



        INTEGER :: ifile, iret, igrib, grib_centre

        CALL grib_open_file(ifile, filepath, 'r', iret)
        IF (iret == 0) THEN
            ! Use first record to detect centre, which is assumed constant
            ! amongst all messages
            CALL grib_new_from_file(ifile, igrib, iret)
            CALL grib_get(igrib, 'centre', grib_centre)
            CALL grib_close_file(ifile)
        ELSE
            PRINT *, "WARNING: problem opening GRIB file: ", filepath
            grib_centre = -999
        END IF





        IF (grib_centre == CENTRE_NCEP) THEN
            gribfile_centre = GRIBFILE_CENTRE_NCEP
        ELSE IF (grib_centre == CENTRE_ECMWF) THEN
            gribfile_centre = GRIBFILE_CENTRE_ECMWF
        ELSE
            gribfile_centre = GRIBFILE_CENTRE_UNKNOWN
        END IF

    END FUNCTION gribfile_centre


    ! This is currently a PRIVATE subroutine
    SUBROUTINE get_centre_and_type(filepath, grib_centre, gribfile_type)
        ! Given specified grib file, passes back centre and gribfile
        ! type to the calling program.  Numeric codes are defined as integer parameters
        ! in this module

        ! To get this information, we have to iterate through the entire file in order to
        ! determine if it is hybrid or not
        !

        USE grib_api
        IMPLICIT NONE
        CHARACTER(LEN=*), INTENT(IN) :: filepath  ! full path to GRIB file
        INTEGER, INTENT(OUT) :: grib_centre, gribfile_type

        INTEGER :: fileptr, iret, igrib, centre, grib_edition
        LOGICAL :: end_of_file
        LOGICAL :: grib1_detected, grib2_detected


        grib1_detected = .FALSE.
        grib2_detected = .FALSE.

        CALL grib_open_file(fileptr, filepath, 'r', iret)
        IF (iret /= 0) THEN
            PRINT *, 'class_gributils:get_centre_and_type()...'
            PRINT *, '     unable to open filepath: ', filepath
            STOP
        END IF


        ! Use first record to detect centre and and grib version of first messages.  We will
        ! then assume that all following messages have same centre, but not necessarily same
        ! GRIB version
        CALL grib_new_from_file(fileptr, igrib, iret)
        CALL grib_get(igrib, 'centre', grib_centre)
        CALL grib_get(igrib, 'edition', grib_edition)

        IF (grib_edition == 1) grib1_detected = .TRUE.
        IF (grib_edition == 2) grib2_detected = .TRUE.

        ! Now, iterate through the rest of records to determine if this is a mixed edition file
        end_of_file = .FALSE.
        DO WHILE (.NOT. end_of_file)
            CALL grib_new_from_file(fileptr, igrib, iret)
            IF (iret .eq. GRIB_END_OF_FILE) THEN
                end_of_file = .TRUE.
            ELSE

                ! Get edition from file
                CALL grib_get(igrib, 'edition', grib_edition)
                IF (grib_edition .eq. 1) grib1_detected = .TRUE.
                IF (grib_edition .eq. 2) grib2_detected = .TRUE.
            END IF
        END DO

    CALL grib_close_file(fileptr)

    ! Determine the gribfile type depending on centre and edition(s)
    IF (grib_centre == CENTRE_ECMWF) THEN
        IF (grib1_detected .AND. grib2_detected) THEN
            gribfile_type = GRIBFILE_TYPE_ECMWF_GRIB1_2
        ELSE IF (grib1_detected .AND. .NOT. grib2_detected) THEN
            gribfile_type = GRIBFILE_TYPE_ECMWF_GRIB1
        ELSE IF (.NOT. grib1_detected .AND. grib2_detected) THEN
            gribfile_type = GRIBFILE_TYPE_ECMWF_GRIB2
        ELSE
            gribfile_type = GRIBFILE_TYPE_UNKNOWN
        ENDIF
    ELSE IF (grib_centre == CENTRE_NCEP) THEN
        IF (grib1_detected .AND. .NOT. grib2_detected) THEN
             gribfile_type = GRIBFILE_TYPE_NCEP_GRIB1
        ELSE IF (.NOT. grib1_detected .AND. grib2_detected) THEN
            gribfile_type = GRIBFILE_TYPE_NCEP_GRIB2
        ELSE
            gribfile_type = GRIBFILE_TYPE_UNKNOWN
        ENDIF
    ELSE
        gribfile_type = GRIBFILE_TYPE_UNKNOWN
    ENDIF

    END SUBROUTINE get_centre_and_type


    SUBROUTINE get_3d_u_dims(filepath, gribfile_type, xlon, ylat, zlev)

        ! Looks at the 3d u fields in the GRIBFILE to get x and y dims, as well as number of levels
        USE grib_api

        IMPLICIT NONE

        CHARACTER(LEN=*), INTENT(IN) :: filepath  ! full path to GRIB file
        INTEGER, INTENT(IN) :: gribfile_type
        INTEGER, INTENT(OUT) :: xlon, ylat, zlev

        INTEGER :: ifile, iret, igrib, grib_centre
        LOGICAL :: end_of_file

        ! These will be assigned according to type of grib file, then used to filter
        ! for the 3d messages
        ! Name of the key being sought
        CHARACTER(LEN=64) :: keyname_leveltype, keyname_shortname, keyname_level, &
                             keyname_xlon, keyname_ylat

        ! The key value being filtered for
        CHARACTER(LEN=64) :: keyvalue_leveltype, keyvalue_shortname

        ! Actual values read in from the grib file
        CHARACTER(LEN=64) :: value_leveltype, value_shortname
        INTEGER :: value_level

        INTEGER :: num_levels

        ! Get the field names to read, based on the type of grib file

        !!! Note that ALL of these have the same key names, except that
        !!! leveltype is 'hybrid' in ECMWF and 'isobaricInhPa' in NCEP.
        !!! This could probably be consolidated, but because these are
        !!! files that go through some preprocessing, and aren't 
        !!! necessarily standard (at least for ECMWF), I'm going to be
        !!! safe and leave it as is for now, so it would be easier to
        !!! modify for one type, if necessary.
        IF (gribfile_type == GRIBFILE_TYPE_ECMWF_GRIB1) THEN
            keyname_leveltype = 'typeOfLevel'
            keyname_shortname = 'shortName'
            keyname_level = 'level'
            keyvalue_leveltype = 'hybrid'
            keyvalue_shortname = 'u'
            keyname_xlon = 'Ni'
            keyname_ylat = 'Nj'
        ELSE IF (gribfile_type == GRIBFILE_TYPE_ECMWF_GRIB1_2) THEN
            keyname_leveltype = 'typeOfLevel'
            keyname_shortname = 'shortName'
            keyname_level = 'level'
            keyvalue_leveltype = 'hybrid'
            keyvalue_shortname = 'u'
            keyname_xlon = 'Ni'
            keyname_ylat = 'Nj'
        ELSE IF (gribfile_type == GRIBFILE_TYPE_ECMWF_GRIB2) THEN
            keyname_leveltype = 'typeOfLevel'
            keyname_shortname = 'shortName'
            keyname_level = 'level'
            keyvalue_leveltype = 'hybrid'
            keyvalue_shortname = 'u'
            keyname_xlon = 'Ni'
            keyname_ylat = 'Nj'
        ELSE IF (gribfile_type == GRIBFILE_TYPE_NCEP_GRIB1) THEN
            keyname_leveltype = 'typeOfLevel'
            keyname_shortname = 'shortName'
            keyname_level = 'level'
            keyvalue_leveltype = 'isobaricInhPa'
            keyvalue_shortname = 'u'
            keyname_xlon = 'Ni'
            keyname_ylat = 'Nj'
        ELSE IF (gribfile_type == GRIBFILE_TYPE_NCEP_GRIB2) THEN
            keyname_leveltype = 'typeOfLevel'
            keyname_shortname = 'shortName'
            keyname_level = 'level'
            keyvalue_leveltype = 'isobaricInhPa'
            keyvalue_shortname = 'u'
            keyname_xlon = 'Ni'
            keyname_ylat = 'Nj'
        ELSE
            PRINT *, 'class_gribfile:get_3d_u_dims(): Unsupported gribfile type: ', gribfile_type
            STOP
        ENDIF

        CALL grib_open_file(ifile, filepath, 'r', iret)
        IF (iret == 0) THEN

            ! Iterate through all messages to count 3d u messages (levels) and get x,y dimensions
            end_of_file = .FALSE.
            num_levels = 0
            DO WHILE (.NOT. end_of_file)
                CALL grib_new_from_file(ifile, igrib, iret)
                IF (iret .eq. GRIB_END_OF_FILE) THEN
                    end_of_file = .TRUE.
                ELSE

                    ! Get relevant keys and filter for the 3d U wind
                    CALL grib_get(igrib, keyname_shortname, value_shortname)
                    CALL grib_get(igrib, keyname_leveltype, value_leveltype)
                    CALL grib_get(igrib, keyname_level, value_level)
                    IF ( TRIM(value_leveltype) == TRIM(keyvalue_leveltype) .AND. &
                         TRIM(value_shortname) == TRIM(keyvalue_shortname) ) THEN

                        ! If this is first 3d U wind message, get dimensions
                        IF (num_levels == 0) THEN
                            CONTINUE
                            CALL grib_get(igrib, keyname_xlon, xlon)
                            CALL grib_get(igrib, keyname_ylat, ylat)
                        ENDIF
                        !PRINT *, TRIM(value_shortname), '  ', TRIM(value_leveltype), '  ', value_level
                        num_levels = num_levels + 1
                    END IF
                END IF
            END DO


        ELSE
            PRINT *, "ERROR: class_gribfile::get_3d_u_dims(): problem opening GRIB file: ", filepath
            STOP
        END IF

        CALL grib_close_file(ifile)

        zlev = num_levels

    END SUBROUTINE get_3d_u_dims



    !!! Getter methods
    INTEGER FUNCTION gribfile_num_xlon(gribfile_obj)

        ! Returns x (lon) dimension of met data
        TYPE(gribfile_object), INTENT(IN) :: gribfile_obj

        IF (.NOT. gribfile_obj%is_instantiated) THEN
            PRINT *, 'ERROR: class_gribfile: gribfile_obj not instantiated'
        ENDIF
        gribfile_num_xlon = gribfile_obj%num_xlon

    END FUNCTION gribfile_num_xlon

    INTEGER FUNCTION gribfile_num_ylat(gribfile_obj)

        ! Returns y (lat) dimension of met data
        TYPE(gribfile_object), INTENT(IN) :: gribfile_obj

        IF (.NOT. gribfile_obj%is_instantiated) THEN
            PRINT *, 'ERROR: class_gribfile: gribfile_obj not instantiated'
        ENDIF
        gribfile_num_ylat = gribfile_obj%num_ylat

    END FUNCTION gribfile_num_ylat

    INTEGER FUNCTION gribfile_num_zlevel(gribfile_obj)

        ! Returns z (level) dimension of met data
        TYPE(gribfile_object), INTENT(IN) :: gribfile_obj

        IF (.NOT. gribfile_obj%is_instantiated) THEN
            PRINT *, 'ERROR: class_gribfile: gribfile_obj not instantiated'
        ENDIF
        gribfile_num_zlevel = gribfile_obj%num_zlevel

    END FUNCTION gribfile_num_zlevel

END MODULE class_gribfile
