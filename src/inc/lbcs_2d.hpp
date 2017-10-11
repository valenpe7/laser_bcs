#pragma once

#include "param_2d.hpp"
#include "global.hpp"

class lbcs_2d {
public:
	lbcs_2d(const param_2d* param);
	void set_param(const param_2d* param);
	void calculate_fields(std::string output_path);
	~lbcs_2d() = default;
private:
	void prescribe_field_at_focus(array_2d<complex>& field) const;
	void calculate_transverse_electric_field();
	void calculate_longitudinal_electric_field();
	void calculate_magnetic_field();
	void dft_time(array_2d<complex>& field, int sign) const;
	void dft_space(array_2d<complex>& field, int sign) const;
	void normalize(array_2d<complex>& field) const;
	void dump_field(array_2d<complex> field, std::string name, std::string output_path) const;
	void dump_to_shared_file(array_2d<complex> field, std::array<int, 4> local_extent, std::array<int, 4> global_extent, std::string filename) const;
	const param_2d* param;
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