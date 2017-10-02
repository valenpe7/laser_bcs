#pragma once

#include "param_3d.hpp"
#include "global.hpp"

class lbcs_3d {
public:
	lbcs_3d(const param_3d* param);
	void set_param(const param_3d* param);
	void calculate_fields();
	void dump_fields(std::string output_path) const;
	~lbcs_3d() = default;
private:
	void prescribe_field_at_focus(array_3d<complex>& field) const;
	void calculate_transverse_electric_field();
	void calculate_longitudinal_electric_field();
	void calculate_magnetic_field();
	void dft_time(array_3d<complex>& field, int sign) const;
	void dft_space(array_3d<complex>& field, int sign) const;
	void dump_to_shared_file(array_3d<complex> field, std::array<int, 6> local_extent, std::array<int, 6> global_extent, std::string filename) const;
	const param_3d* param;
	array_3d<complex> e_x;
	array_3d<complex> e_y;
	array_3d<complex> e_z;
	array_3d<complex> b_x;
	array_3d<complex> b_y;
	array_3d<complex> b_z;
	array_3d<double> k_z;
	std::vector<double> k_x;
	std::vector<double> k_y;
	std::vector<double> omega;
	bool fields_computed;
};