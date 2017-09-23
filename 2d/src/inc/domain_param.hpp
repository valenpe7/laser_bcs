#pragma once

#include <vector>
#include <array>

class domain_param {
	friend class laser_bcs;
public:
	domain_param(int rank, int size, double z_boundary, double z_focus, double x_min, double x_max, int nx,
		int cpml, double t_max, double dx, double dt);
	void set_values(int rank, int size, double z_boundary, double z_focus, double x_min, double x_max, int nx,
		int cpml, double t_max, double dx, double dt);
	~domain_param() = default;
protected:
	int rank;
	int size;
	std::array<double, 2> x_lim;
	std::array<double, 2> t_lim;
	double z_boundary;
	double z_focus;
	double time_shift;
	int ghost_cells;
	int cpml;
	double dx;
	double dt;
	int nx;
	int nt;
	int nt_global;
	int nt_start;
	std::vector<int> t_counts;
	std::vector<int> t_displs;
	std::vector<double> x_coord;
	std::vector<double> t_coord;
};
