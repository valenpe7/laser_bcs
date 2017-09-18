#include <iostream>
#include <fstream>
#include <string>
#include <complex>
#include <cmath>
#include <mpi.h>

#include "inc/laser_bcs.hpp"

laser_bcs::laser_bcs(const domain_param* dp, const laser_param* lp) {
	this->set_param(dp, lp);
	this->fields_computed = false;
}

void laser_bcs::set_param(const domain_param* dp, const laser_param* lp) {
	this->domain = dp;
	this->laser = lp;
	this->e_x.resize(boost::extents[this->domain->nx][this->domain->ny][this->domain->nt]);
	this->e_y.resize(boost::extents[this->domain->nx][this->domain->ny][this->domain->nt]);
	this->e_z.resize(boost::extents[this->domain->nx][this->domain->ny][this->domain->nt]);
	this->b_x.resize(boost::extents[this->domain->nx][this->domain->ny][this->domain->nt]);
	this->b_y.resize(boost::extents[this->domain->nx][this->domain->ny][this->domain->nt]);
	this->b_z.resize(boost::extents[this->domain->nx][this->domain->ny][this->domain->nt]);
	int half_nt = static_cast<int>(ceil(this->domain->nt_global / 2.0));
	int half_nx = static_cast<int>(ceil(this->domain->nx / 2.0));
	int half_ny = static_cast<int>(ceil(this->domain->ny / 2.0));
  std::vector<int> tmp1(2 * half_nt), tmp2(2* half_nx), tmp3(2 * half_ny);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for(size_t i = 0; i < tmp1.size(); i++) {
    tmp1[i] = (i < half_nt) ? i : 0;
  }
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
 for(size_t i = 0; i < tmp2.size(); i++) {
    tmp2[i] = (i < half_nx) ? i : i - 2 * half_nx;
  }
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
 for(size_t i = 0; i < tmp3.size(); i++) {
    tmp3[i] = (i < half_ny) ? i : i - 2 * half_ny;
  }
	this->omega.resize(this->domain->nt);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
	for(auto i = 0; i < this->domain->nt; i++) {
		this->omega[i] = 2.0 * constants::pi * static_cast<double>(tmp1[i + this->domain->nt_start]) / (this->domain->dt * this->domain->nt_global);
	}
	this->k_x.resize(this->domain->nx);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
	for(auto i = 0; i < this->domain->nx; i++) {
		this->k_x[i] = 2.0 * constants::pi * static_cast<double>(tmp2[i]) / (this->domain->dx * this->domain->nx);
	}
  this->k_y.resize(this->domain->ny);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
	for(auto i = 0; i < this->domain->ny; i++) {
		this->k_y[i] = 2.0 * constants::pi * static_cast<double>(tmp3[i]) / (this->domain->dy * this->domain->ny);
	}
	this->k_z.resize(boost::extents[this->domain->nx][this->domain->ny][this->domain->nt]);
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(static)
#endif
	for(auto i = 0; i < this->domain->nx; i++) {
		for(auto j = 0; j < this->domain->ny; j++) {
			for(auto k = 0; k < this->domain->nt; k++) {
				k_z[i][j][k] = std::real(sqrt(static_cast<std::complex<double>>(pow(this->omega[k] / constants::c, 2) - pow(k_x[i], 2) - pow(k_y[j], 2))));
			}
		}
	}
}

void laser_bcs::prescribe_field_at_focus(array_3d<complex>& field) const {
	if(this->laser->t_start < this->domain->t_lim[0] || (this->laser->t_start + this->domain->time_shift) > this->domain->t_lim[1]) {
		std::cout << "Warning: pulse not captured at focus" << std::endl;
	}
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(static)
#endif
	for(auto i = 0; i < this->domain->nx; i++) {
		for(auto j = 0; j < this->domain->ny; j++) {
			for(auto k = 0; k < this->domain->nt; k++) {
				if((this->domain->t_coord[k] - this->domain->time_shift) >= this->laser->t_start && (this->domain->t_coord[k] - this->domain->time_shift) <= this->laser->t_end) {
					field[i][j][k] = {this->laser->amp * exp(
						- pow((this->domain->x_coord[i] - this->laser->x_0) / this->laser->w_0, 2)
						- pow((this->domain->y_coord[j] - this->laser->y_0) / this->laser->w_0, 2)
						- pow((this->domain->t_coord[k] - this->laser->t_0 - this->domain->time_shift) * (2.0 * sqrt(log(2.0)))
							/ this->laser->fwhm_time, 2)) * cos(this->laser->omega * this->domain->t_coord[k]), 0.0};

				} else {
					field[i][j][k] = {0.0, 0.0};
				}
			}
		}
	}
}

void laser_bcs::dft_time(array_3d<complex>& field, int sign) const {
	std::vector<complex> global_extent(this->domain->nt_global);
	fft::create_plan_1d(this->domain->nt_global, sign);
	for(auto i = 0; i < this->domain->nx; i++) {
		for(auto j = 0; j < this->domain->ny; j++) {
			MPI_Gatherv(field[boost::indices[i][j][range()]].origin(), this->domain->nt, MPI_DOUBLE_COMPLEX, global_extent.data(), this->domain->t_counts.data(), this->domain->t_displs.data(), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
			if(this->domain->rank == 0) {
				global_extent = fft::execute_plan(global_extent);
			}
			MPI_Scatterv(global_extent.data(), this->domain->t_counts.data(), this->domain->t_displs.data(), MPI_DOUBLE_COMPLEX, field[boost::indices[i][j][range()]].origin(), this->domain->nt, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
		}
	}
	fft::destroy_plan();
}

void laser_bcs::dft_space(array_3d<complex>& field, int sign) const {
	fft::create_plan_2d(this->domain->nx, this->domain->ny, sign);
	for(auto k = 0; k < this->domain->nt; k++) {
		view_2d slice = field[boost::indices[range()][range()][k]];
		tools::vec_to_slice(slice, fft::execute_plan(tools::slice_to_vec(slice)));
	}
	fft::destroy_plan();
}

void laser_bcs::calculate_transverse_electric_field() {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(static)
#endif
	for(auto i = 0; i < this->domain->nx; i++) {
		for(auto j = 0; j < this->domain->ny; j++) {
			for(auto k = 0; k < this->domain->nt; k++) {
				if(this->k_z[i][j][k] > 0) {
					this->e_x[i][j][k] *= exp(constants::imag_unit * this->k_z[i][j][k] * (this->domain->z_boundary - this->domain->z_focus));
					this->e_y[i][j][k] *= exp(constants::imag_unit * this->k_z[i][j][k] * (this->domain->z_boundary - this->domain->z_focus));
				} else {
					this->e_x[i][j][k] = {0.0, 0.0};
					this->e_y[i][j][k] = {0.0, 0.0};
				}
			}
		}
	}
}

void laser_bcs::calculate_longitudinal_electric_field() {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(static)
#endif
	for(auto i = 0; i < this->domain->nx; i++) {
		for(auto j = 0; j < this->domain->ny; j++) {
			for(auto k = 0; k < this->domain->nt; k++) {
				if(this->k_z[i][j][k] > 0) {
					this->e_z[i][j][k] = -(this->k_x[i] * this->e_x[i][j][k] + this->k_y[j] * this->e_y[i][j][k]) / this->k_z[i][j][k];
				} else {
					this->e_z[i][j][k] = {0.0, 0.0};
				}
			}
		}
	}
}

void laser_bcs::calculate_magnetic_field() {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(static)
#endif
	for(auto i = 0; i < this->domain->nx; i++) {
		for(auto j = 0; j < this->domain->ny; j++) {
			for(auto k = 0; k < this->domain->nt; k++) {
				if(this->k_z[i][j][k] > 0) {
					this->b_x[i][j][k] = (-this->k_x[i] * this->k_y[j] * this->e_x[i][j][k] + (pow(this->k_x[i], 2) - pow(this->omega[k] / constants::c, 2) * this->e_y[i][j][k]))
					/ (this->omega[k] * this->k_z[i][j][k]);
					this->b_y[i][j][k] = ((pow(this->omega[k] / constants::c, 2) - pow(this->k_y[j], 2)) * this->e_x[i][j][k] + this->k_x[i] * this->k_y[j] * this->e_y[i][j][k])
					/ (this->omega[k] * this->k_z[i][j][k]);
					this->b_z[i][j][k] = (-this->k_y[j] * this->e_x[i][j][k] + this->k_x[i] * this->e_y[i][j][k]) / this->omega[k];
				} else {
					this->b_x[i][j][k] = {0.0, 0.0};
					this->b_y[i][j][k] = {0.0, 0.0};
					this->b_z[i][j][k] = {0.0, 0.0};
				}
			}
		}
	}
}

void laser_bcs::dump_fields(std::string output_path) const {
	if(!this->fields_computed) {
		std::cout << "warning: cannot dump fiedls - fields are not computed" << std::endl;
		return;
	}
	std::array<int, 6> local_extent = {0, this->domain->nx - 1, 0, this->domain->ny - 1, this->domain->nt_start, this->domain->nt_start + this->domain->nt - 1};
	std::array<int, 6> global_extent = {0, this->domain->nx - 1, 0, this->domain->ny - 1, 0, this->domain->nt_global - this->domain->ghost_cells - 1};
	this->dump_to_shared_file(this->e_x, local_extent, global_extent, output_path + "/e_x_" + std::to_string(this->laser->id) + ".raw");
	this->dump_to_shared_file(this->e_y, local_extent, global_extent, output_path + "/e_y_" + std::to_string(this->laser->id) + ".raw");
	this->dump_to_shared_file(this->e_z, local_extent, global_extent, output_path + "/e_z_" + std::to_string(this->laser->id) + ".raw");
	this->dump_to_shared_file(this->b_x, local_extent, global_extent, output_path + "/b_x_" + std::to_string(this->laser->id) + ".raw");
	this->dump_to_shared_file(this->b_y, local_extent, global_extent, output_path + "/b_y_" + std::to_string(this->laser->id) + ".raw");
	this->dump_to_shared_file(this->b_z, local_extent, global_extent, output_path + "/b_z_" + std::to_string(this->laser->id) + ".raw");
}

void laser_bcs::dump_to_shared_file(array_3d<complex> field, std::array<int, 6> local_extent, std::array<int, 6> global_extent, std::string filename) const {
	MPI_File file;
	MPI_Offset offset = 0;
	MPI_Status status;
	MPI_Datatype local_array;
	if(local_extent[4] > global_extent[5]) {
		return;
	}
	if(local_extent[5] > global_extent[5]) {
		local_extent[5] = global_extent[5];
	}
	std::array<int, 3> size_local = {local_extent[1] - local_extent[0] + 1, local_extent[3] - local_extent[2] + 1, local_extent[5] - local_extent[4] + 1};
	std::array<int, 3> size_global = {global_extent[1] - global_extent[0] + 1, global_extent[3] - global_extent[2] + 1, global_extent[5] - global_extent[4] + 1};
	std::array<int, 3> start_coords = {local_extent[0], local_extent[2], local_extent[4]};
	MPI_Type_create_subarray(3, size_global.data(), size_local.data(), start_coords.data(), MPI_ORDER_FORTRAN, MPI_DOUBLE, &local_array);
	MPI_Type_commit(&local_array);
	MPI_File_open(MPI_COMM_WORLD, filename.data(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
	MPI_File_set_view(file, offset, MPI_DOUBLE, local_array, "native", MPI_INFO_NULL);
	MPI_File_write_all(file, tools::get_real(field, size_local).data(), size_local[0] * size_local[1] * size_local[2], MPI_DOUBLE, &status);
	MPI_File_close(&file);
	MPI_Type_free(&local_array);
}

void laser_bcs::run_computation() {
	this->prescribe_field_at_focus(this->e_x);
	this->prescribe_field_at_focus(this->e_y);
	this->dft_time(this->e_x, 1);
	this->dft_time(this->e_y, 1);
	this->dft_space(this->e_x, -1);
	this->dft_space(this->e_y, -1);
	this->calculate_transverse_electric_field();
	this->calculate_longitudinal_electric_field();
	this->calculate_magnetic_field();
	this->dft_space(this->e_x, 1);
	this->dft_space(this->e_y, 1);
	this->dft_space(this->e_z, 1);
	this->dft_space(this->b_x, 1);
	this->dft_space(this->b_y, 1);
	this->dft_space(this->b_z, 1);
	this->dft_time(this->e_x, -1);
	this->dft_time(this->e_y, -1);
	this->dft_time(this->e_z, -1);
	this->dft_time(this->b_x, -1);
	this->dft_time(this->b_y, -1);
	this->dft_time(this->b_z, -1);
	tools::multiply_array(this->e_x, 1.0 / (this->domain->nx * this->domain->ny * this->domain->nt_global));
	tools::multiply_array(this->e_y, 1.0 / (this->domain->nx * this->domain->ny * this->domain->nt_global));
	tools::multiply_array(this->e_z, 1.0 / (this->domain->nx * this->domain->ny * this->domain->nt_global));
	tools::multiply_array(this->b_x, 1.0 / (this->domain->nx * this->domain->ny * this->domain->nt_global));
	tools::multiply_array(this->b_y, 1.0 / (this->domain->nx * this->domain->ny * this->domain->nt_global));
	tools::multiply_array(this->b_z, 1.0 / (this->domain->nx * this->domain->ny * this->domain->nt_global));
	this->fields_computed = true;
}
