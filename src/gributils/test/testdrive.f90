PROGRAM testdrive

    USE class_gribfile

    IMPLICIT NONE

    CHARACTER(LEN=256) :: the_file_path

    INTEGER, PARAMETER :: NCASES = 6
    CHARACTER(LEN=256), DIMENSION(NCASES) :: file_paths
    CHARACTER(LEN=256), DIMENSION(NCASES) :: descriptions

    TYPE(gribfile_object) :: my_gribfile

    INTEGER :: case_number

    !!!!!!!!!!!!!!!!!!!
    ! These define the header and the grib file used for each test
    descriptions(1) = "ECMWF GRIB1 on global 1.0 degree domain"
    file_paths(1) = "../../../devtest/case_data/met_data/ecmwf/t1_03h_ec1p0d/EN13062503"

    descriptions(2) = "ECMWF GRIB1/2 on tiny domain"
    file_paths(2) = "../../../devtest/case_data/met_data/ecmwf/t1_33h_ec1p0d/EL14091909"

    descriptions(3) = "ECMWF GRIB1/2 on global 1.0 degree domain"
    file_paths(3) = "../../../devtest/case_data/met_data/ecmwf/t1_03h_ec1p0d_grib1-2/EE13110700"

    descriptions(4) = "ECMWF GRIB2 on global 1.0 degree domain"
    file_paths(4) = "../../../devtest/case_data/met_data/ecmwf/t1_03h_ec1p0d_grib2/EN13110700"

    descriptions(5) = "NCEP GRIB1 on global 1.0 degree domain"
    file_paths(5) = "../../../devtest/case_data/met_data/ncep/t1_06h_nc1p0d_grib1/GD05051406"

    descriptions(6) = "NCEP GRIB2 on global 1.0 degree domain"
    file_paths(6) = "../../../devtest/case_data/met_data/ncep/t1_03h_nc1p0d/GF15021603"
    !!!!!!!!!!!!!!!!!!!

    DO case_number = 1,NCASES
        PRINT *,
        PRINT *, TRIM( descriptions( case_number) )

        my_gribfile = gribfile_object_create( file_paths(case_number) )
        PRINT *,
        PRINT *, 'Output from calling gribfile_printobj()...'
        PRINT *, '++++++++++'
        CALL gribfile_printobj(my_gribfile)
        PRINT *, '++++++++++'
        PRINT *,
        
        PRINT *, 'Output from the getter methods...'
        PRINT *, 'gribfile_center()   : ', gribfile_centre( file_paths(case_number) )

        PRINT *, 'gribfile_num_xlon   : ', gribfile_num_xlon(my_gribfile)
        PRINT *, 'gribfile_num_ylat   : ', gribfile_num_ylat(my_gribfile)
        PRINT *, 'gribfile_num_zlevel : ', gribfile_num_zlevel(my_gribfile)
        PRINT *,
        PRINT *, '----------------------------'
    ENDDO


END PROGRAM testdrive
