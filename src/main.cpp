#include <iostream>
#include <fstream>
#ifdef _DEBUG
	#include <chrono>
	#include <mpi.h>
#endif

#include "inc/lbcs_2d.hpp"
#include "inc/lbcs_3d.hpp"

#ifdef __cplusplus
extern "C" {
#endif
	void calculate_2d(int* rank, int* size, double* t_start, double* t_end, double* fwhm_time, double* t_0, double* omega,
		double* amplitude, double* x_0, double* w_0, int* order, int* direction, int* id, double* phase, double* z_boundary, double* z_focus, double* x_min, double* x_max,
		int* nx, int* cpml, double* t_max, double* dx, double* dt, const char* data_dir);
	void calculate_3d(int* rank, int* size, double* t_start, double* t_end, double* fwhm_time, double* t_0, double* omega,
		double* amplitude, double* x_0, double* y_0, double* w_0, int* order, int * direction, int* id, double* phase, double* z_boundary, double* z_focus, double* x_min, double* x_max,
		double* y_min, double* y_max, int* nx, int* ny, int* cpml, double* t_max, double* dx, double* dy, double* dt, const char* data_dir);
	void retrieve_2d(double* buffer, int* id, const char* data_dir, const char* field, int* timestep,
		int* size_global, int* first, int* last);
	void retrieve_3d(double* buffer, int* id, const char* data_dir, const char* field, int* timestep,
		int* horizontal_global, int* vertical_global, int* horizontal_first, int* horizontal_last, int* vertical_first, int* vertical_last);
#ifdef __cplusplus
}
#endif

void calculate_2d(int* rank, int* size, double* t_start, double* t_end, double* fwhm_time, double* t_0, double* omega,
	double* amplitude, double* x_0, double* w_0, int* order, int* direction, int* id, double* phase, double* z_boundary, double* z_focus, double* x_min, double* x_max,
	int* nx, int* cpml, double* t_max, double* dx, double* dt, const char* data_dir) {
	param_2d param;
	param.set_domain(*rank, *size, *z_boundary, *z_focus, *x_min, *x_max, *nx, *cpml, *t_max, *dx, *dt);
	param.set_laser(*t_start, *t_end, *fwhm_time, *t_0, *x_0, *omega, *amplitude, *w_0, *order, *direction, *id, *phase);
#ifdef _DEBUG
	MPI_Barrier(MPI_COMM_WORLD);
	auto t1 = std::chrono::high_resolution_clock::now();
	if (*rank == 0) {
		std::cout << "Calculating fields of laser " << *id << " at boundary..." << std::endl;
	}
#endif
	lbcs_2d lbcs(&param);
	lbcs.calculate_fields(std::string(data_dir));
#ifdef _DEBUG
	MPI_Barrier(MPI_COMM_WORLD);
	auto t2 = std::chrono::high_resolution_clock::now();
	if (*rank == 0) {
		std::cout << "Fields of laser " << *id << " dumped in: " << data_dir << " (runtime: " << std::chrono::duration_cast<std::chrono::duration<float>>(t2 - t1).count() << " s)" << std::endl;
	}
#endif
}

void calculate_3d(int* rank, int* size, double* t_start, double* t_end, double* fwhm_time, double* t_0, double* omega,
	double* amplitude, double* x_0, double* y_0, double* w_0, int* order, int * direction, int* id, double* phase, double* z_boundary, double* z_focus, double* x_min, double* x_max,
	double* y_min, double* y_max, int* nx, int* ny, int* cpml, double* t_max, double* dx, double* dy, double* dt, const char* data_dir) {
	param_3d param;
	param.set_domain(*rank, *size, *z_boundary, *z_focus, *x_min, *x_max, *nx, *y_min, *y_max, *ny, *cpml, *t_max, *dx, *dy, *dt);
	param.set_laser(*t_start, *t_end, *fwhm_time, *t_0, *x_0, *y_0, *omega, *amplitude, *w_0, *order, *direction, *id, *phase);
#ifdef _DEBUG
	MPI_Barrier(MPI_COMM_WORLD);
	auto t1 = std::chrono::high_resolution_clock::now();
	if (*rank == 0) {
		std::cout << "Calculating fields of laser " << *id << " at boundary..." << std::endl;
	}
#endif
	lbcs_3d lbcs(&param);
	lbcs.calculate_fields(std::string(data_dir));
#ifdef _DEBUG
	MPI_Barrier(MPI_COMM_WORLD);
	auto t2 = std::chrono::high_resolution_clock::now();
	if (*rank == 0) {
		std::cout << "Fields of laser " << *id << " dumped in: " << data_dir << " (runtime: " << std::chrono::duration_cast<std::chrono::duration<float>>(t2 - t1).count() << " s)" << std::endl;
	}
#endif
}

void retrieve_2d(double* buffer, int* id, const char* data_dir, const char* field, int* timestep, int* size_global, int* first, int* last) {
	double num = 0.0;
	std::string laser_id = std::to_string(*id);
	std::string path(data_dir);
	std::string name(field);
	std::ifstream in;
	in.open(path + "/" + name + "_" + laser_id + ".raw", std::ios::binary);
	if (in.is_open()) {
		in.seekg(((*timestep) * (*size_global) + (*first) - 1) * sizeof(num));
		for (auto i = 0; i < *last - *first + 1; i++) {
			in.read(reinterpret_cast<char*>(&num), sizeof(num));
			buffer[i] = num;
		}
		in.close();
	}
	else {
		std::cout << "Error: cannot read file " << path + "/" + name + "_" + laser_id + ".raw" << std::endl;
	}
}

void retrieve_3d(double* buffer, int* id, const char* data_dir, const char* field, int* timestep, int* horizontal_global, int* vertical_global,
	int* horizontal_first, int* horizontal_last, int* vertical_first, int* vertical_last) {
	double num = 0.0;
	std::string laser_id = std::to_string(*id);
	std::string path(data_dir);
	std::string name(field);
	std::ifstream in;
	in.open(path + "/" + name + "_" + laser_id + ".raw", std::ios::binary);
	if (in.is_open()) {
		for (auto j = 0; j < *vertical_last - *vertical_first + 1; j++) {
			in.seekg((((*timestep) * (*vertical_global) + (*vertical_first + j - 1)) * (*horizontal_global) + (*horizontal_first - 1)) * sizeof(num));
			for (auto i = 0; i < *horizontal_last - *horizontal_first + 1; i++) {
				in.read(reinterpret_cast<char*>(&num), sizeof(num));
				buffer[i + j * (*horizontal_last - *horizontal_first + 1)] = num;
			}
		}
		in.close();
	}
	else {
		std::cout << "Error: cannot read file " << path + "/" + name + "_" + laser_id + ".raw" << std::endl;
	}
}
