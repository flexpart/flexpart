
module class_vtable

    implicit none
    private

    ! Note that all public interfaces and variables should have a 
    ! "vtable" prefix
    public :: Vtable,                    &
              Vtable_record,             &
              vtable_load_by_name,       &
              vtable_dump_records,       &
              vtable_get_fpname,         &
              vtable_gribfile_inventory

    integer, parameter :: VTABLE_MISSING_ENTRY = -9999
                          
    ! These codes depict the origin/model of the GRIB File
    integer, parameter :: GRIB_CENTRE_NCEP = 7, &
                          GRIB_CENTRE_ECMWF = 98

    type Vtable_record
        character(len=15) :: fp_name
        character(len=15) :: fp_units
        character(len=25) :: fp_description        
        integer :: grib1_param
        integer :: grib1_leveltype
        integer :: grib2_discipline
        integer :: grib2_category
        integer :: grib2_param
        integer :: grib2_leveltype
    end type Vtable_record

    type Vtable
        logical :: initialized=.FALSE. 
        character(len=255) :: path_to_vtable_file
        integer :: num_entries = 0
        type(Vtable_record), allocatable :: the_entries(:)
    end type Vtable

contains




    subroutine vtable_load_by_name(vtable_name, the_vtable_data)
        implicit none

        character(len=*), intent(in) :: vtable_name  ! Full path to vtable file

        logical :: lexist
        integer :: ierr
        integer :: num_vrecs = 0 
        integer :: vrec_idx
        character(len=255) :: file_line = ' ' 

        type(Vtable), intent(out) :: the_vtable_data    ! Data structure holding the vtable 

        type(Vtable_record) :: vrecord

        ! Make sure the file exists
        inquire(file=trim(vtable_name), exist=lexist)
        if (.not. lexist) then
            print *, 'file: ', trim(vtable_name), ' does not exist...'
            stop
        endif

        ! Open file
        open(10, file=trim(vtable_name), status='old', form='formatted', iostat=ierr)
        if (ierr .ne. 0) then
            print *, 'file: ', trim(vtable_name), ' open failed...'
            stop
        endif

        ! Go through the file once and count the vtable_records
        ! Read past headers
        file_line = ' '
        do while(file_line(1:1) .ne. '-')
            read(10, '(A255)', iostat=ierr) file_line
        enddo 

        ! Now we are at the '----------' line - process everything between
        ! here and the next '----------' line.  In this case, we just want to
        ! count
        file_line = ' '
        num_vrecs = 0
        do while(file_line(1:1) .ne. '-')
            read(10, '(A255)', iostat=ierr) file_line
            !print *, file_line
            num_vrecs = num_vrecs + 1
        enddo 
        num_vrecs = num_vrecs - 1

        ! print *, 'num_vrecs: ', num_vrecs

        ! Rewind
        rewind(unit=10)

        ! Allocate array for storing the vtable records, and store
        ! num_entries
        !print *, 'Ready to allocate the_vtable_data'
        allocate(the_vtable_data%the_entries(num_vrecs))
        !print *, 'Allocated the_vtable_data'
        the_vtable_data%num_entries = num_vrecs

        ! Read, parse and store the vtable records
        ! Read past headers
        file_line = ' '
        do while(file_line(1:1) .ne. '-')
            read(10, '(A255)', iostat=ierr) file_line
            !print *, file_line
        enddo 

        ! Now we are at the '----------' line - process everything between
        ! here and the next '----------' line.  In this case, we just want to
        ! count
        file_line = ' '
        vrec_idx = 0
        do while(file_line(1:1) .ne. '-')
            read(10, '(A255)', iostat=ierr) file_line
            if (file_line(1:1) .ne. '-') then
                ! PROCESS THE LINE
                vrec_idx = vrec_idx + 1

                ! Parse the line and put it in the vtable structure
                the_vtable_data%the_entries(vrec_idx) = vtable_parse_record(file_line)
              
                !print *, the_vtable_data%the_entries(vrec_idx) 
                !print *, file_line
                !print *, 'hello'
            endif
        enddo 
        num_vrecs = num_vrecs - 1

        ! Close the file
        close(unit=10)

        the_vtable_data%initialized = .TRUE. 
        
        !print *, the_vtable_data%the_entries(1)
    end subroutine vtable_load_by_name



    type(Vtable_record) function vtable_parse_record(vtable_line)

    !!! Using a vtable line as input argument, parses into a Vtable_record, and returns that
    !!! record
        implicit none
        character(LEN=255), intent(in) :: vtable_line

        !!! This is a sample of what a Vtable header and first two lines
        !!! will look like
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GRIB1| Level| flexpart | flexpart | flexpart      |GRIB2|GRIB2|GRIB2|GRIB2|
! Param| Type | Name     | Units    | Description   |Discp|Catgy|Param|Level|
! -----+------+----------+----------+---------------+-----------------------+
!  130 | 109  | TT       | K        | Temperature   |  0  |  0  |  0  | 105 |
!  131 | 109  | UU       | m s-1    | U             |  0  |  2  |  2  | 105 |
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        ! Storage for Vtable record tokens
        character(25) :: token_fp_name='',&
                         token_fp_units='', &
                         token_fp_description='', &
                         token_grib1_param='', &
                         token_grib1_leveltype='', &
                         token_grib2_discipline='', &
                         token_grib2_category='', &
                         token_grib2_param='', &
                         token_grib2_leveltype=''

        ! These indices mark the locations of the '|' delimiter in a Vtable record
        integer :: delim_fp_name_idx, &
                   delim_fp_units_idx, &
                   delim_fp_description_idx, &
                   delim_grib1_param_idx, &
                   delim_grib1_leveltype_idx, &
                   delim_grib2_discipline_idx, &
                   delim_grib2_category_idx, &
                   delim_grib2_param_idx, &
                   delim_grib2_leveltype_idx

        type(Vtable_record) :: vrecord

        integer :: istat    ! Error indicator for some I/O routines

        ! Calculate the indices of each field so we can extract later
        delim_grib1_param_idx = index(vtable_line, '|')
        delim_grib1_leveltype_idx = index(vtable_line(delim_grib1_param_idx+1:), '|') &
                               + delim_grib1_param_idx
        delim_fp_name_idx = index(vtable_line(delim_grib1_leveltype_idx+1:), '|') &
                               + delim_grib1_leveltype_idx
        delim_fp_units_idx = index(vtable_line(delim_fp_name_idx+1:), '|') &
                            + delim_fp_name_idx
        delim_fp_description_idx = index(vtable_line(delim_fp_units_idx+1:), '|') &
                           + delim_fp_units_idx
        delim_grib2_discipline_idx = index(vtable_line(delim_fp_description_idx+1:), '|') &
                            + delim_fp_description_idx
        delim_grib2_category_idx = index(vtable_line(delim_grib2_discipline_idx+1:), '|') &
                            + delim_grib2_discipline_idx
        delim_grib2_param_idx = index(vtable_line(delim_grib2_category_idx+1:), '|') &
                            + delim_grib2_category_idx
        delim_grib2_leveltype_idx = index(vtable_line(delim_grib2_param_idx+1:), '|') &
                            + delim_grib2_param_idx

        ! Extract the tokens
        token_grib1_param = vtable_line(:delim_grib1_param_idx-1)
        token_grib1_leveltype = vtable_line(delim_grib1_param_idx+1:delim_grib1_leveltype_idx-1)
        token_fp_name = ADJUSTL(vtable_line(delim_grib1_leveltype_idx+1:delim_fp_name_idx-1))
        token_fp_units = ADJUSTL(vtable_line(delim_fp_name_idx+1:delim_fp_units_idx-1))
        token_fp_description = ADJUSTL(vtable_line(delim_fp_units_idx+1:delim_fp_description_idx-1))
        token_grib2_discipline = vtable_line(delim_fp_description_idx+1:delim_grib2_discipline_idx-1)
        token_grib2_category = vtable_line(delim_grib2_discipline_idx+1:delim_grib2_category_idx-1)
        token_grib2_param = vtable_line(delim_grib2_category_idx+1:delim_grib2_param_idx-1)
        token_grib2_leveltype = vtable_line(delim_grib2_param_idx+1:delim_grib2_leveltype_idx-1)

        ! Jam the data in the record for return
        ! I "think" I have the read statement as a way to cleanly
        ! convert the strings to integers
        read(token_grib1_param, *, iostat=istat) vrecord%grib1_param
        if (istat .ne. 0) vrecord%grib1_param = VTABLE_MISSING_ENTRY

        read(token_grib1_leveltype, *, iostat=istat) vrecord%grib1_leveltype
        if (istat .ne. 0) vrecord%grib1_leveltype = VTABLE_MISSING_ENTRY

        vrecord%fp_name = token_fp_name
        vrecord%fp_units = token_fp_units
        vrecord%fp_description = token_fp_description

        read(token_grib2_discipline, *, iostat=istat) vrecord%grib2_discipline
        if (istat .ne. 0) vrecord%grib2_discipline = VTABLE_MISSING_ENTRY

        read(token_grib2_category, *, iostat=istat) vrecord%grib2_category
        if (istat .ne. 0) vrecord%grib2_category = VTABLE_MISSING_ENTRY

        read(token_grib2_param, *, iostat=istat) vrecord%grib2_param
        if (istat .ne. 0) vrecord%grib2_param = VTABLE_MISSING_ENTRY

        read(token_grib2_leveltype, *, iostat=istat) vrecord%grib2_leveltype
        if (istat .ne. 0) vrecord%grib2_leveltype = VTABLE_MISSING_ENTRY

        vtable_parse_record = vrecord

        !print *, "Hello vtable_parse_record()"
        !print *, vrecord
    end function vtable_parse_record


    SUBROUTINE vtable_dump_records(the_vtable_data)
    
        IMPLICIT NONE
        
        !  Prints the contents of the Vtable that's been loaded into
        !  class/module.  Original intent is for debugging dumps

        type(Vtable), intent(in) :: the_vtable_data    ! Data structure holding the vtable 

        type(Vtable_record) :: vrecord        

        INTEGER :: i
                
        IF (the_vtable_data%initialized) THEN
            DO i=1,the_vtable_data%num_entries
                WRITE(6,1000) the_vtable_data%the_entries(i)%grib1_param, &
                &          the_vtable_data%the_entries(i)%grib1_leveltype, &
                &          the_vtable_data%the_entries(i)%fp_name, &
                &          the_vtable_data%the_entries(i)%fp_units, &
                &          the_vtable_data%the_entries(i)%fp_description, &
                &          the_vtable_data%the_entries(i)%grib2_discipline, &
                &          the_vtable_data%the_entries(i)%grib2_category, &
                &          the_vtable_data%the_entries(i)%grib2_param, &
                &          the_vtable_data%the_entries(i)%grib2_leveltype
            ENDDO
        ELSE
            WRITE(6,*) 'Unitialized Vtable - nothing to dump...'
        END IF

1000 FORMAT(1X, I6, 1X, I3, 1X, A10, 1X, A10, 1X, A25, 1X, I3, 1X, I3, 1X, &
     &      I3, 1X, I3)
        
    END SUBROUTINE vtable_dump_records



    CHARACTER(LEN=15) FUNCTION vtable_get_fpname(igrib, vtable_object)
    
        !!! Assumes that a calling routine has opened up a GRIB file and obtained the
        !!! grib id for a specific message.
        !!! Given a grib message and a Vtable, looks up the message parameters in the Vtable
        !!! and, if found, returns the fpname
        
        use grib_api
        implicit none
        
        integer, intent(in) :: igrib
        type(Vtable), intent(in) :: vtable_object
        
        integer :: parameter_id, category, number, discipline, edition, surface_type, &
                   level, indicator_of_parameter
        character(len=10) :: center
        
        integer :: idx
        logical :: record_match
        
        call grib_get(igrib, 'editionNumber', edition)
        call grib_get(igrib, 'level', level)
                
        if (edition .eq. 1) then
            call grib_get(igrib, 'indicatorOfParameter', indicator_of_parameter)
            call grib_get(igrib, 'indicatorOfTypeOfLevel', surface_type)
            !print *, '(edition, indicator_of_parameter, surftype, level): ', edition, indicator_of_parameter, surface_type,&
            !          level
        else if (edition .eq. 2) then
            call grib_get(igrib, 'parameterNumber', number)
            call grib_get(igrib, 'parameterCategory', category)
            call grib_get(igrib, 'discipline', discipline)
            call grib_get(igrib, 'typeOfFirstFixedSurface', surface_type)
            !print *, '(edition, number, cat, disc, surftype, level): ', edition, number, &
            !          category, discipline, surface_type, level
        else
            print *, 'vtable_get_fpname(): Illegal edition: ', edition
            stop
        endif       

        ! Iterate through Vtable and look for a match
        vtable_get_fpname = 'None'
        record_match = .FALSE.
        idx = 1
        do while (.NOT. record_match .AND. idx .LE. vtable_object%num_entries) 

            if (edition .eq. 1) then
                if (vtable_object%the_entries(idx)%grib1_param .eq. indicator_of_parameter .and. &
                    vtable_object%the_entries(idx)%grib1_leveltype .eq. surface_type) then
                    vtable_get_fpname = vtable_object%the_entries(idx)%fp_name
                    record_match = .TRUE.
                end if                               
            else if (edition .eq. 2) then
                if (vtable_object%the_entries(idx)%grib2_discipline .eq. discipline .and.    &
                    vtable_object%the_entries(idx)%grib2_param .eq. number .and.   &
                    vtable_object%the_entries(idx)%grib2_category .eq. category .and.   &
                    vtable_object%the_entries(idx)%grib2_leveltype .eq. surface_type) then
                
                    vtable_get_fpname = vtable_object%the_entries(idx)%fp_name
                    record_match = .TRUE.
                end if
            else
                print *, 'Illegal edition: ', edition
                stop
            endif                  
            idx = idx + 1    
        end do
        
        
    END FUNCTION vtable_get_fpname

    SUBROUTINE vtable_gribfile_inventory(gribfile_path, my_vtable)


        ! Given full path to a gribfile, and an initialized Vtable, prints an
        ! inventory of the gribfile relative to the vtable.
        

        USE grib_api

        IMPLICIT NONE

        CHARACTER(LEN=*), INTENT(IN) :: gribfile_path  ! Full path to grib file
        TYPE(Vtable), INTENT(IN) :: my_vtable            ! Data structure holding the vtable 

        ! Parallel arrays to keep track of information in gribfile
        CHARACTER(LEN=15), ALLOCATABLE :: fpnames(:)
        CHARACTER(LEN=15) :: fpname
        INTEGER, ALLOCATABLE :: num_messages(:)
        
        INTEGER :: i, fpname_idx

        INTEGER :: grib_fh, iret, grib_msg_id
        LOGICAL :: file_open_success, end_of_grib_file


        file_open_success = .FALSE.
        CALL grib_open_file(grib_fh, gribfile_path, 'r', iret)
        IF (iret .EQ. GRIB_SUCCESS) THEN
            file_open_success = .TRUE.
            PRINT *, 'Opened grib file...'
        ELSE
            WRITE(6,*) 'vtable_gribfile_inventory(): Failed to open grib file...', &
        &        TRIM(gribfile_path)
        ENDIF

        IF (.NOT. my_vtable%initialized) THEN
            PRINT *, 'vtable_gribfile_inventory(): Vtable was not initialized...'
        END IF

        IF (file_open_success .AND. my_vtable%initialized) THEN
        
            !PRINT *, "gribfile Inventory: ", TRIM(gribfile_path)
            !PRINT *, ' '
            !PRINT *, ' '
            !CALL vtable_dump_records(my_vtable)
            !PRINT *, ' '
            !PRINT *, ' '
            
            ! Create the "counting" arrays - add one for "non FP"
            ALLOCATE(fpnames(my_vtable%num_entries+1))
            ALLOCATE(num_messages(my_vtable%num_entries+1))
            DO i=1,my_vtable%num_entries
                fpnames(i) = my_vtable%the_entries(i)%fp_name
            END DO
            fpnames(my_vtable%num_entries+1) = 'NOFP'
            num_messages = 0  ! This is a full array that's getting initialized
            
            !PRINT *, 'fpnames: ', fpnames
            !PRINT *, 'num_messages: ', num_messages
            

            ! Iterate through the gribfile and count occurences of each variable            
            end_of_grib_file = .FALSE.
            DO WHILE (.NOT. end_of_grib_file) 
                CALL grib_new_from_file(grib_fh, grib_msg_id, iret)
                IF (iret .EQ. GRIB_END_OF_FILE) THEN
                    end_of_grib_file = .TRUE.
                ELSE
                    !PRINT *, 'Reading grib_msg_id: ', grib_msg_id
                    fpname = vtable_get_fpname(grib_msg_id, my_vtable)         
                    !PRINT *, '    fpname: ', TRIM(fpname)
                    
                    ! First thing is change fpname to the NOFP label if
                    ! it didn't match a Vtable entry
                    IF (TRIM(fpname) .eq. "None" ) fpname = "NOFP"
                    !PRINT *, '    fpname: ', TRIM(fpname)
                    
                    ! Find the index in the counting array
                    fpname_idx = findidx(fpnames, fpname)
                    
                    ! Increment appropriate index
                    IF (fpname_idx .GT. 0) THEN
                        num_messages(fpname_idx) = num_messages(fpname_idx) + 1
                    END IF                   
                ENDIF
            END DO 

            WRITE(6,*) "GRIBFILE Inventory: ", TRIM(gribfile_path)
            WRITE(6,*)
            WRITE(6,*) "FPNAME    OCCURRENCES"
            DO i=1,SIZE(fpnames)
                WRITE(6,2000) ADJUSTL(fpnames(i)), num_messages(i)
            END DO      
            
        END IF


2000 FORMAT(1X, A10, 3x, I3)
  
    END SUBROUTINE vtable_gribfile_inventory    



    INTEGER FUNCTION findidx(myarr, key)
    
        ! Returns index of key in myarr
        ! Returns 0 if not found
    
        IMPLICIT NONE
    
        CHARACTER(LEN=15), DIMENSION (:), INTENT(IN) :: myarr
        CHARACTER(LEN=15), INTENT(IN) :: key
    
        LOGICAL :: found
        INTEGER :: i, myarr_size
    
        !PRINT *, 'myarr in findidx: ', myarr
        !PRINT *, 'size of myarr in findidx: ', SIZE(myarr)
    
        myarr_size = SIZE(myarr)
        findidx = 0
        found = .FALSE.
        i = 1
        DO WHILE(.NOT. found .AND. i .LE. myarr_size)
            IF (TRIM(myarr(i)) .EQ. TRIM(key)) THEN
                findidx = i
                found = .TRUE.
            ENDIF
            i = i + 1
    
        END DO
    
    
        RETURN
    
    END FUNCTION findidx





end module class_vtable


