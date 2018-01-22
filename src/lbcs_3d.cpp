#include <iostream>
#include <fstream>
#include <string>
#include <complex>
#include <cmath>
#include <mpi.h>

#include "inc/lbcs_3d.hpp"

lbcs_3d::lbcs_3d(param_3d param) {
	this->fields_computed = false;
	this->param = std::make_unique<param_3d>(param);
	this->e_x.resize(boost::extents[this->param->nx][this->param->ny][this->param->nt]);
	this->e_y.resize(boost::extents[this->param->nx][this->param->ny][this->param->nt]);
	this->e_z.resize(boost::extents[this->param->nx][this->param->ny][this->param->nt]);
	this->b_x.resize(boost::extents[this->param->nx][this->param->ny][this->param->nt]);
	this->b_y.resize(boost::extents[this->param->nx][this->param->ny][this->param->nt]);
	this->b_z.resize(boost::extents[this->param->nx][this->param->ny][this->param->nt]);
	int half_nt = static_cast<int>(ceil(this->param->nt_global / 2.0));
	int half_nx = static_cast<int>(ceil(this->param->nx / 2.0));
	int half_ny = static_cast<int>(ceil(this->param->ny / 2.0));
	std::vector<int> tmp1(2 * half_nt), tmp2(2 * half_nx), tmp3(2 * half_ny);
	if (this->param->direction == 1) {
		for (auto i = 0; i < static_cast<int>(tmp1.size()); i++) {
			tmp1[i] = (i < half_nt) ? i : 0;
		}
	}
	else {
		for (auto i = 0; i < static_cast<int>(tmp1.size()); i++) {
			tmp1[i] = (i < half_nt) ? 0 : i - 2 * half_nt;
		}
	}
	for (auto i = 0; i < static_cast<int>(tmp2.size()); i++) {
		tmp2[i] = (i < half_nx) ? i : i - 2 * half_nx;
	}
	for (auto i = 0; i < static_cast<int>(tmp3.size()); i++) {
		tmp3[i] = (i < half_ny) ? i : i - 2 * half_ny;
	}
	this->omega.resize(this->param->nt);
	for (auto i = 0; i < this->param->nt; i++) {
		this->omega[i] = 2.0 * constants::pi * static_cast<double>(tmp1[i + this->param->nt_start]) / (this->param->dt * this->param->nt_global);
	}
	this->k_x.resize(this->param->nx);
	for (auto j = 0; j < this->param->nx; j++) {
		this->k_x[j] = 2.0 * constants::pi * static_cast<double>(tmp2[j]) / (this->param->dx * this->param->nx);
	}
	this->k_y.resize(this->param->ny);
	for (auto j = 0; j < this->param->ny; j++) {
		this->k_y[j] = 2.0 * constants::pi * static_cast<double>(tmp3[j]) / (this->param->dy * this->param->ny);
	}
	this->k_z.resize(boost::extents[this->param->nx][this->param->ny][this->param->nt]);
	for (auto i = 0; i < this->param->nx; i++) {
		for (auto j = 0; j < this->param->ny; j++) {
			for (auto k = 0; k < this->param->nt; k++) {
				k_z[i][j][k] = std::real(sqrt(static_cast<complex>(pow(this->omega[k] / constants::c, 2) - pow(k_x[i], 2) - pow(k_y[j], 2))));
			}
		}
	}
}

void lbcs_3d::prescribe_field_at_focus(m_array<complex, 3>& field) const {
	if (this->param->t_start < this->param->t_lim[0] || (this->param->t_start + this->param->time_shift) > this->param->t_lim[1]) {
		if (this->param->rank == 0) std::cout << "Warning: pulse not captured at focus" << std::endl;
	}
	for (auto i = 0; i < this->param->nx; i++) {
		for (auto j = 0; j < this->param->ny; j++) {
			for (auto k = 0; k < this->param->nt; k++) {
				if ((this->param->t_coord[k] - this->param->time_shift) >= this->param->t_start && (this->param->t_coord[k] - this->param->time_shift) <= this->param->t_end) {
					field[i][j][k] = { this->param->amp * exp(
						-pow((this->param->x_coord[i] - this->param->x_0) / this->param->w_0, 2)
						- pow((this->param->y_coord[j] - this->param->y_0) / this->param->w_0, 2)
						- pow((this->param->t_coord[k] - this->param->t_0 - this->param->time_shift) * (2.0 * sqrt(log(2.0)))
							/ this->param->fwhm_time, 2)) * cos(this->param->omega * (this->param->t_coord[k] - this->param->t_0 - this->param->time_shift) + this->param->phase), 0.0 };
				}
				else {
					field[i][j][k] = { 0.0, 0.0 };
				}
			}
		}
	}
}

void lbcs_3d::dft_time(m_array<complex, 3>& field, int sign) const {
	std::vector<complex> global_extent(this->param->nt_global);
	fft::create_plan_1d(this->param->nt_global, sign);
	for (auto i = 0; i < this->param->nx; i++) {
		for (auto j = 0; j < this->param->ny; j++) {
			MPI_Gatherv(field[boost::indices[i][j][range()]].origin(), this->param->nt, MPI_DOUBLE_COMPLEX, global_extent.data(), const_cast<int*>(this->param->t_counts.data()), const_cast<int*>(this->param->t_displs.data()), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
			if (this->param->rank == 0) global_extent = fft::execute_plan(global_extent);
			MPI_Scatterv(global_extent.data(), const_cast<int*>(this->param->t_counts.data()), const_cast<int*>(this->param->t_displs.data()), MPI_DOUBLE_COMPLEX, field[boost::indices[i][j][range()]].origin(), this->param->nt, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
		}
	}
	fft::destroy_plan();
}

void lbcs_3d::dft_space(m_array<complex, 3>& field, int sign) const {
	fft::create_plan_2d(this->param->nx, this->param->ny, sign);
	for (auto k = 0; k < this->param->nt; k++) {
		auto slice = field[boost::indices[range()][range()][k]];
		tools::vec_to_array(slice, fft::execute_plan(tools::array_to_vec(slice)));
	}
	fft::destroy_plan();
}

void lbcs_3d::calculate_transverse_electric_field() {
	for (auto i = 0; i < this->param->nx; i++) {
		for (auto j = 0; j < this->param->ny; j++) {
			for (auto k = 0; k < this->param->nt; k++) {
				if (this->k_z[i][j][k] > 0) {
					this->e_x[i][j][k] *= exp(constants::imag_unit * this->k_z[i][j][k] * (this->param->z_boundary - this->param->z_focus));
					this->e_y[i][j][k] *= exp(constants::imag_unit * this->k_z[i][j][k] * (this->param->z_boundary - this->param->z_focus));
				}
				else {
					this->e_x[i][j][k] = { 0.0, 0.0 };
					this->e_y[i][j][k] = { 0.0, 0.0 };
				}
			}
		}
	}
}

void lbcs_3d::calculate_longitudinal_electric_field() {
	for (auto i = 0; i < this->param->nx; i++) {
		for (auto j = 0; j < this->param->ny; j++) {
			for (auto k = 0; k < this->param->nt; k++) {
				if (this->k_z[i][j][k] > 0) {
					this->e_z[i][j][k] = -(this->k_x[i] * this->e_x[i][j][k] + this->k_y[j] * this->e_y[i][j][k]) / this->k_z[i][j][k];
				}
				else {
					this->e_z[i][j][k] = { 0.0, 0.0 };
				}
			}
		}
	}
}

void lbcs_3d::calculate_magnetic_field() {
	for (auto i = 0; i < this->param->nx; i++) {
		for (auto j = 0; j < this->param->ny; j++) {
			for (auto k = 0; k < this->param->nt; k++) {
				if (this->k_z[i][j][k] > 0) {
					this->b_x[i][j][k] = (-this->k_x[i] * this->k_y[j] * this->e_x[i][j][k] + (pow(this->k_x[i], 2) - pow(this->omega[k] / constants::c, 2) * this->e_y[i][j][k]))
						/ (this->omega[k] * this->k_z[i][j][k]);
					this->b_y[i][j][k] = ((pow(this->omega[k] / constants::c, 2) - pow(this->k_y[j], 2)) * this->e_x[i][j][k] + this->k_x[i] * this->k_y[j] * this->e_y[i][j][k])
						/ (this->omega[k] * this->k_z[i][j][k]);
					this->b_z[i][j][k] = (-this->k_y[j] * this->e_x[i][j][k] + this->k_x[i] * this->e_y[i][j][k]) / this->omega[k];
				}
				else {
					this->b_x[i][j][k] = { 0.0, 0.0 };
					this->b_y[i][j][k] = { 0.0, 0.0 };
					this->b_z[i][j][k] = { 0.0, 0.0 };
				}
			}
		}
	}
}

void lbcs_3d::normalize(m_array<complex, 3>& field) const {
	tools::multiply_array<complex, 3>(field, 2.0 / (this->param->nx * this->param->ny * this->param->nt_global));
}

void lbcs_3d::dump_field(m_array<complex, 3> field, std::string name, std::string output_path) const {
	if (!this->fields_computed) {
		if (this->param->rank == 0) std::cout << "Warning: cannot dump field " << name << " - fields are not computed" << std::endl;
		return;
	}
	std::array<int, 6> local_extent = { 0, this->param->nx - 1, 0, this->param->ny - 1, this->param->nt_start, this->param->nt_start + this->param->nt - 1 };
	std::array<int, 6> global_extent = { 0, this->param->nx - 1, 0, this->param->ny - 1, 0, this->param->nt_global - this->param->ghost_cells - 1 };
	this->dump_to_shared_file(field, local_extent, global_extent, output_path + "/" + name + "_" + std::to_string(this->param->id) + ".raw");
}

void lbcs_3d::dump_to_shared_file(m_array<complex, 3> field, std::array<int, 6> local_extent, std::array<int, 6> global_extent, std::string filename) const {
	MPI_File file;
	MPI_Offset offset = 0;
	MPI_Status status;
	MPI_Datatype local_array;
	MPI_Group group, group_world;
	MPI_Comm comm;
	std::vector<int> members;
	if (local_extent[4] <= global_extent[5]) members.push_back(this->param->rank);
	MPI_Comm_group(MPI_COMM_WORLD, &group_world);
	MPI_Group_incl(group_world, static_cast<int>(members.size()), members.data(), &group);
	MPI_Comm_create(MPI_COMM_WORLD, group, &comm);
	if (local_extent[4] > global_extent[5]) return;
	if (local_extent[5] > global_extent[5]) local_extent[5] = global_extent[5];
	std::array<int, 3> size_local = { local_extent[1] - local_extent[0] + 1, local_extent[3] - local_extent[2] + 1, local_extent[5] - local_extent[4] + 1 };
	std::array<int, 3> size_global = { global_extent[1] - global_extent[0] + 1, global_extent[3] - global_extent[2] + 1, global_extent[5] - global_extent[4] + 1 };
	std::array<int, 3> start_coords = { local_extent[0], local_extent[2], local_extent[4] };
	MPI_Type_create_subarray(3, size_global.data(), size_local.data(), start_coords.data(), MPI_ORDER_FORTRAN, MPI_DOUBLE, &local_array);
	MPI_Type_commit(&local_array);
	MPI_File_open(comm, const_cast<char*>(filename.data()), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
	MPI_File_set_view(file, offset, MPI_DOUBLE, local_array, "native", MPI_INFO_NULL);
	MPI_File_write_all(file, tools::get_real(field, size_local).data(), size_local[0] * size_local[1] * size_local[2], MPI_DOUBLE, &status);
	MPI_File_close(&file);
	MPI_Type_free(&local_array);
}

void lbcs_3d::calculate_fields(std::string output_path) {
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
	this->normalize(this->e_x);
	this->normalize(this->e_y);
	this->normalize(this->e_z);
	this->normalize(this->b_x);
	this->normalize(this->b_y);
	this->normalize(this->b_z);
	this->fields_computed = true;
	this->dump_field(this->e_x, "e_x", output_path);
	this->dump_field(this->e_y, "e_y", output_path);
	this->dump_field(this->e_z, "e_z", output_path);
	this->dump_field(this->b_x, "b_x", output_path);
	this->dump_field(this->b_y, "b_y", output_path);
	this->dump_field(this->b_z, "b_z", output_path);
}
