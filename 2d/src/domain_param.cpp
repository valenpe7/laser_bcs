#include <iostream>
#include <cmath>
#include <mpi.h>

#include "inc/domain_param.hpp"
#include "inc/global.hpp"

domain_param::domain_param(int rank, int size, double z_boundary, double z_focus, double x_min, double x_max, int nx,
	int cpml, double t_max, double dx, double dt) {
	this->set_values(rank, size, z_boundary, z_focus, x_min, x_max, nx, cpml, t_max, dx, dt);
}

void domain_param::set_values(int rank, int size, double z_boundary, double z_focus, double x_min, double x_max, int nx,
	int cpml, double t_max, double dx, double dt) {
	if(x_min > x_max || t_max < 0.0) {
		std::cerr << "error: bad value" << std::endl;
	}
	this->rank = rank;
	this->size = size;
	this->x_lim = {x_min, x_max};
	this->z_boundary = z_boundary;
	this->z_focus = z_focus;
	this->time_shift = std::abs((this->z_boundary - this->z_focus) / constants::c);
	this->t_lim = {0.0, t_max + time_shift};
	this->dx = dx;
	this->dt = dt;
	this->ghost_cells = static_cast<int>(ceil(this->time_shift / this->dt));
	this->cpml = cpml;
	this->nx = nx;
	this->nt_global = static_cast<int>(ceil((this->t_lim[1] - this->t_lim[0]) / this->dt));
	this->nt = static_cast<int>(round(this->nt_global / this->size));
	this->nt_start = this->rank * this->nt;
	if(this->rank == this->size - 1) {
		if((this->nt_global - this->size * this->nt) != 0) {
			this->nt += (this->nt_global - this->size * this->nt);
		}
	}
	this->t_counts.resize(this->size);
	this->t_displs.resize(this->size);
	MPI_Gather(&this->nt, 1, MPI_INT, t_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Gather(&this->nt_start, 1, MPI_INT, t_displs.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
	this->x_coord.resize(this->nx);
	this->t_coord.resize(this->nt);
	for(auto i = 0; i < this->nx; i++) {
		this->x_coord[i] = this->x_lim[0] + (i - this->cpml) * this->dx;
	}
	for(auto k = 0; k < this->nt; k++) {
		this->t_coord[k] = this->t_lim[0] + (this->nt_start + k) * this->dt;
	}
}
