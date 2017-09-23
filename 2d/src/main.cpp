#include <iostream>
#include <fstream>
#include <chrono>
#include <mpi.h>

#include "inc/laser_bcs.hpp"

#ifdef __cplusplus
extern "C" {
#endif
	void compute_laser_at_boundary(int* rank, int* size, double* t_start, double* t_end, double* fwhm_time, double* t_0, double* omega,
		double* amplitude, double* x_0, double* w_0, int* direction, int* id, double* z_boundary, double* z_focus, double* x_min, double* x_max,
		int* nx, int* cpml, double* t_max, double* dx, double* dt, const char* data_dir);
	void populate_laser_at_boundary(double* buffer, int* id, const char* data_dir, const char* field, int* timestep, int* size_global, int* first, int* last);
#ifdef __cplusplus
}
#endif

void compute_laser_at_boundary(int* rank, int* size, double* t_start, double* t_end, double* fwhm_time, double* t_0, double* omega,
	double* amplitude, double* x_0, double* w_0, int* direction, int* id, double* z_boundary, double* z_focus, double* x_min, double* x_max,
	int* nx, int* cpml, double* t_max, double* dx, double* dt, const char* data_dir) {

	std::string output_path(data_dir);

	domain_param domain(*rank, *size, *z_boundary, *z_focus, *x_min, *x_max, *nx, *cpml, *t_max, *dx, *dt);
	laser_param laser(*t_start, *t_end, *fwhm_time, *t_0, *x_0, *omega, *amplitude, *w_0, *direction, *id);

	MPI_Barrier(MPI_COMM_WORLD);
	auto t1 = std::chrono::high_resolution_clock::now();

	laser_bcs bcs(&domain, &laser);
	bcs.run_computation();
	bcs.dump_fields(output_path);

	MPI_Barrier(MPI_COMM_WORLD);
	auto t2 = std::chrono::high_resolution_clock::now();

	if(*rank == 0) {
		std::cout << " Fields of laser " << *id << " dumped in: " << output_path << " (runtime: " << std::chrono::duration_cast<std::chrono::duration<float>>(t2 - t1).count() << " s)" << std::endl;
	}
}

void populate_laser_at_boundary(double* buffer, int* id, const char* data_dir, const char* field, int* timestep, int* size_global, int* first, int* last) {
	double num = 0.0;
	std::string laser_id = std::to_string(*id);
	std::string path(data_dir);
	std::string name(field);
	std::ifstream in;
	in.open(path + name + laser_id + ".raw", std::ios::binary);
	if(in.is_open()) {
		in.seekg(((*timestep) * (*size_global) + (*first) - 1) * sizeof(num));
		for(auto i = 0; i < *last - *first + 1; i++) {
			in.read(reinterpret_cast<char*>(&num), sizeof(num));
			buffer[i] = num;
		}
		in.close();
	} else {
		std::cout << "error: cannot read file " << path + name + laser_id + ".raw" << std::endl;
	}
}
