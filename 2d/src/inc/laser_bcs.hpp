#pragma once

#include "domain_param.hpp"
#include "laser_param.hpp"
#include "global.hpp"

class laser_bcs {
public:
	laser_bcs(const domain_param* dp, const laser_param* lp);
	void set_param(const domain_param* dp, const laser_param* lp);
	void run_computation();
	void dump_fields(std::string output_path) const;
	~laser_bcs() = default;
private:
	void prescribe_field_at_focus(array_2d<complex>& field) const;
	void calculate_transverse_electric_field();
	void calculate_longitudinal_electric_field();
	void calculate_magnetic_field();
	void dft_time(array_2d<complex>& field, int sign) const;
	void dft_space(array_2d<complex>& field, int sign) const;
	void dump_to_shared_file(array_2d<complex> field, std::array<int, 4> local_extent, std::array<int, 4> global_extent, std::string filename) const;
	const domain_param* domain;
	const laser_param* laser;
	array_2d<complex> e_x;
	array_2d<complex> e_y;
	array_2d<complex> e_z;
	array_2d<complex> b_x;
	array_2d<complex> b_y;
	array_2d<complex> b_z;
	array_2d<double> k_z;
	std::vector<double> k_x;
	std::vector<double> omega;
	bool fields_computed;
};
