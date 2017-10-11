MODULE tight_focusing

  USE shared_data
  USE, INTRINSIC :: iso_c_binding
  IMPLICIT NONE
 
INTERFACE
  
  SUBROUTINE calculate_3d(rank, nproc, laser_start, laser_end, &
      fwhm_time, t_0, omega, amplitude, horizontal_pos, vertical_pos, w_0, direction, id, &
      boundary, focus, horizontal_min, horizontal_max, vertical_min, &
      vertical_max, cells_horizontal, cells_vertical, cpml_thickness, t_end, &
      cell_size_horizontal, cell_size_vertical, dt, output_path) bind(c)

    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    
    INTEGER(c_int), INTENT(IN) :: rank, nproc, direction, id, cells_horizontal, cells_vertical, cpml_thickness
    CHARACTER(kind=c_char), DIMENSION(*), INTENT(IN) :: output_path
    REAL(c_double), INTENT(IN) :: laser_start, laser_end, fwhm_time, t_0, omega, amplitude, &
      horizontal_pos, vertical_pos, w_0, boundary, focus, horizontal_min, horizontal_max, &
      vertical_min, vertical_max, t_end, cell_size_horizontal, cell_size_vertical, dt
  
  END SUBROUTINE calculate_3d
  
  SUBROUTINE retrieve_3d(buffer, laser_id, output_path, field, timestep, horizontal_global, &
      vertical_global, horizontal_first, horizontal_last, vertical_first, vertical_last) bind(c)
    
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE
    
    INTEGER(c_int), INTENT(IN) :: laser_id, timestep, horizontal_global, vertical_global, &
      horizontal_first, horizontal_last, vertical_first, vertical_last
    CHARACTER(kind=c_char), DIMENSION(*), INTENT(IN) :: output_path, field
    REAL(c_double), DIMENSION(*), INTENT(OUT) :: buffer
  
  END SUBROUTINE retrieve_3d

END INTERFACE

CONTAINS

  SUBROUTINE Maxwell_consistent_computation_of_EM_fields
    
    TYPE(laser_block), POINTER :: current

    current => laser_x_min
    DO WHILE(ASSOCIATED(current))
      CALL calculate_3d(rank, nproc, current%t_start, current%t_end, &
        current%fwhm_time, current%t_0, current%omega, current%amp, current%horizontal_pos, &
        current%vertical_pos, current%w_0, 1, current%id, x_min, current%focus, y_min, y_max, & 
        z_min, z_max, ny_global, nz_global, cpml_thickness, t_end, dy, dz, dt, TRIM(data_dir)//C_NULL_CHAR)
      current => current%next
    ENDDO

    current => laser_x_max
    DO WHILE(ASSOCIATED(current))
      CALL calculate_3d(rank, nproc, current%t_start, current%t_end, &
        current%fwhm_time, current%t_0, current%omega, current%amp, current%horizontal_pos, &
        current%vertical_pos, current%w_0, -1, current%id, x_max, current%focus, y_min, y_max, & 
        z_min, z_max, ny_global, nz_global, cpml_thickness, t_end, dy, dz, dt, TRIM(data_dir)//C_NULL_CHAR)
      current => current%next
    ENDDO

    current => laser_y_min
    DO WHILE(ASSOCIATED(current))
      CALL calculate_3d(rank, nproc, current%t_start, current%t_end, &
        current%fwhm_time, current%t_0, current%omega, current%amp, current%horizontal_pos, &
        current%vertical_pos, current%w_0, 1, current%id, y_min, current%focus, x_min, x_max, & 
        z_min, z_max, nx_global, nz_global, cpml_thickness, t_end, dx, dz, dt, TRIM(data_dir)//C_NULL_CHAR)
      current => current%next
    ENDDO
      
    current => laser_y_max
    DO WHILE(ASSOCIATED(current))
      CALL calculate_3d(rank, nproc, current%t_start, current%t_end, &
        current%fwhm_time, current%t_0, current%omega, current%amp, current%horizontal_pos, &
        current%vertical_pos, current%w_0, -1, current%id, y_max, current%focus, x_min, x_max, & 
        z_min, z_max, nx_global, nz_global, cpml_thickness, t_end, dx, dz, dt, TRIM(data_dir)//C_NULL_CHAR)
      current => current%next
    ENDDO

    current => laser_z_min
    DO WHILE(ASSOCIATED(current))
      CALL calculate_3d(rank, nproc, current%t_start, current%t_end, &
        current%fwhm_time, current%t_0, current%omega, current%amp, current%horizontal_pos, &
        current%vertical_pos, current%w_0, 1, current%id, z_min, current%focus, x_min, x_max, & 
        y_min, y_max, nx_global, ny_global, cpml_thickness, t_end, dx, dy, dt, TRIM(data_dir)//C_NULL_CHAR)
      current => current%next
    ENDDO

    current => laser_z_max
    DO WHILE(ASSOCIATED(current))
      CALL calculate_3d(rank, nproc, current%t_start, current%t_end, &
        current%fwhm_time, current%t_0, current%omega, current%amp, current%horizontal_pos, &
        current%vertical_pos, current%w_0, -1, current%id, z_max, current%focus, x_min, x_max, & 
        y_min, y_max, nx_global, ny_global, cpml_thickness, t_end, dx, dy, dt, TRIM(data_dir)//C_NULL_CHAR)
      current => current%next
    ENDDO

  END SUBROUTINE Maxwell_consistent_computation_of_EM_fields

  SUBROUTINE get_source_x_boundary(buffer, laser_id)

    REAL(num), DIMENSION(:,:), INTENT(INOUT) :: buffer
    INTEGER, INTENT(IN) :: laser_id

    CALL retrieve_3d(buffer, laser_id, TRIM(data_dir)//C_NULL_CHAR, "e_x"//C_NULL_CHAR, &
      step, ny_global, nz_global, ny_global_min, ny_global_max, nz_global_min, nz_global_max)

  END SUBROUTINE get_source_x_boundary
  
  SUBROUTINE get_source_y_boundary(buffer, laser_id)
    
    REAL(num), DIMENSION(:,:), INTENT(INOUT) :: buffer
    INTEGER, INTENT(IN) :: laser_id

    CALL retrieve_3d(buffer, laser_id, TRIM(data_dir)//C_NULL_CHAR, "e_x"//C_NULL_CHAR, &
      step, nx_global, nz_global, nx_global_min, nx_global_max, nz_global_min, nz_global_max)

  END SUBROUTINE get_source_y_boundary

  SUBROUTINE get_source_z_boundary(buffer, laser_id)
    
    REAL(num), DIMENSION(:,:), INTENT(INOUT) :: buffer
    INTEGER, INTENT(IN) :: laser_id

    CALL retrieve_3d(buffer, laser_id, TRIM(data_dir)//C_NULL_CHAR, "e_x"//C_NULL_CHAR, &
      step, nx_global, ny_global, nx_global_min, nx_global_max, ny_global_min, ny_global_max)

  END SUBROUTINE get_source_z_boundary

END MODULE tight_focusing
