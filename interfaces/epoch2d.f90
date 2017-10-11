MODULE tight_focusing

  USE shared_data
  USE, INTRINSIC :: iso_c_binding
  IMPLICIT NONE
 
INTERFACE
  
  SUBROUTINE calculate_2d(rank, nproc, laser_start, laser_end, &
    fwhm_time, t_0, omega, amplitude, horizontal_pos, w_0, direction, id, boundary, focus, horizontal_min, horizontal_max, &
    cells_horizontal, cpml_thickness, t_end, cell_size_horizontal, dt, output_path) bind(c)
    
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    
    INTEGER(c_int), INTENT(IN) :: rank, nproc, direction, id, cells_horizontal, cpml_thickness
    CHARACTER(kind=c_char), DIMENSION(*), INTENT(IN) :: output_path
    REAL(c_double), INTENT(IN) :: laser_start, laser_end, fwhm_time, t_0, omega, amplitude, & 
    horizontal_pos, w_0, boundary, focus, horizontal_min, horizontal_max, t_end, cell_size_horizontal, dt
  
  END SUBROUTINE calculate_2d
  
  SUBROUTINE retrieve_2d(buffer, laser_id, output_path, field, timestep, size_global, first, last) bind(c)
    
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    
    INTEGER(c_int), INTENT(IN) :: laser_id, timestep, size_global, first, last
    CHARACTER(kind=c_char), DIMENSION(*), INTENT(IN) :: output_path, field
    REAL(c_double), DIMENSION(*), INTENT(OUT) :: buffer
  
  END SUBROUTINE retrieve_2d

END INTERFACE

CONTAINS

  SUBROUTINE Maxwell_consistent_computation_of_EM_fields
    
    TYPE(laser_block), POINTER :: current

    current => laser_x_min
    DO WHILE(ASSOCIATED(current))
      CALL calculate_2d(rank, nproc, current%t_start, current%t_end, &
        current%fwhm_time, current%t_0, current%omega, current%amp, current%horizontal_pos, current%w_0, 1, &
        current%id, x_min, current%focus, y_min, y_max, ny_global, cpml_thickness, t_end, &
        dy, dt, TRIM(data_dir)//C_NULL_CHAR)
      current => current%next
    ENDDO

    current => laser_x_max
    DO WHILE(ASSOCIATED(current))
      CALL calculate_2d(rank, nproc, current%t_start, current%t_end, &
        current%fwhm_time, current%t_0, current%omega, current%amp, current%horizontal_pos, current%w_0, -1, &
        current%id, x_max, current%focus, y_min, y_max, ny_global, cpml_thickness, t_end, &
        dy, dt, TRIM(data_dir)//C_NULL_CHAR)
      current => current%next
    ENDDO

    current => laser_y_min
    DO WHILE(ASSOCIATED(current))
      CALL calculate_2d(rank, nproc, current%t_start, current%t_end, &
        current%fwhm_time, current%t_0, current%omega, current%amp, current%horizontal_pos, current%w_0, 1, &
        current%id, y_min, current%focus, x_min, x_max, nx_global, cpml_thickness, t_end, &
        dx, dt, TRIM(data_dir)//C_NULL_CHAR)
      current => current%next
    ENDDO
      
    current => laser_y_max
    DO WHILE(ASSOCIATED(current))
      CALL calculate_2d(rank, nproc, current%t_start, current%t_end, &
        current%fwhm_time, current%t_0, current%omega, current%amp, current%horizontal_pos, current%w_0, -1, &
        current%id, y_max, current%focus, x_min, x_max, nx_global, cpml_thickness, t_end, &
        dx, dt, TRIM(data_dir)//C_NULL_CHAR)
      current => current%next
    ENDDO

  END SUBROUTINE Maxwell_consistent_computation_of_EM_fields

  SUBROUTINE get_source_x_boundary(buffer, laser_id)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: buffer
    INTEGER, INTENT(IN) :: laser_id

    CALL retrieve_2d(buffer, laser_id, TRIM(data_dir)//C_NULL_CHAR, "e_x"//C_NULL_CHAR, &
      step, ny_global, ny_global_min, ny_global_max)

  END SUBROUTINE get_source_x_boundary
  
  SUBROUTINE get_source_y_boundary(buffer, laser_id)
    
    REAL(num), DIMENSION(:), INTENT(INOUT) :: buffer
    INTEGER, INTENT(IN) :: laser_id

    CALL retrieve_2d(buffer, laser_id, TRIM(data_dir)//C_NULL_CHAR, "e_x"//C_NULL_CHAR, &
      step, nx_global, nx_global_min, nx_global_max)

  END SUBROUTINE get_source_y_boundary


END MODULE tight_focusing
