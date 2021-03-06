#pragma once

#include <vector>
#include <array>

class param_2d {
	friend class lbcs_2d;
public:
	param_2d() = default;
	void set_domain(int rank, int size, double z_boundary, double z_focus, double x_min, double x_max, int nx, int cpml, double t_max, double dx, double dt);
	void set_laser(double t_start, double t_end, double fwhm_time, double t_0, double x_0, double omega, double amplitude, double w_0, int order, int direction, int id, double phase);
	~param_2d() = default;
private:
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
	double omega;
	double t_start;
	double t_end;
	double fwhm_time;
	double t_0;
	double x_0;
	double amplitude;
	double w_0;
	int order;
    int direction;
	int id;
	double phase;
};
