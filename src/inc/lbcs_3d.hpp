#pragma once

#include "param_3d.hpp"
#include "global.hpp"

class lbcs_3d {
public:
	lbcs_3d(param_3d param);
	void calculate_fields(std::string output_path);
	~lbcs_3d() = default;
private:
	void prescribe_field_at_focus(m_array<complex, 3>& field) const;
	void calculate_transverse_electric_field();
	void calculate_longitudinal_electric_field();
	void calculate_magnetic_field();
	void dft_time(m_array<complex, 3>& field, int sign) const;
	void dft_space(m_array<complex, 3>& field, int sign) const;
	void normalize(m_array<complex, 3>& field) const;
	void dump_field(m_array<complex, 3> field, std::string output_path) const;
	void dump_to_shared_file(m_array<complex, 3> field, std::array<int, 6> local_extent, std::array<int, 6> global_extent, std::string filename) const;
	std::unique_ptr<param_3d> param;
	m_array<complex, 3> e_x;
	m_array<complex, 3> e_y;
	m_array<complex, 3> e_z;
	m_array<complex, 3> b_x;
	m_array<complex, 3> b_y;
	m_array<complex, 3> b_z;
	m_array<double, 3> k_z;
	std::vector<double> k_x;
	std::vector<double> k_y;
	std::vector<double> omega;
};
