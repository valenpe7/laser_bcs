#pragma once

#include "param_2d.hpp"
#include "global.hpp"

class lbcs_2d {
public:
	lbcs_2d(param_2d param);
	void calculate_fields(std::string output_path);
	~lbcs_2d() = default;
private:
	void prescribe_field_at_focus(m_array<complex, 2>& field) const;
	void calculate_transverse_electric_field();
	void calculate_longitudinal_electric_field();
	void calculate_magnetic_field();
	void dft_time(m_array<complex, 2>& field, int sign) const;
	void dft_space(m_array<complex, 2>& field, int sign) const;
	void normalize(m_array<complex, 2>& field) const;
	void dump_field(m_array<complex, 2> field, std::string name, std::string output_path) const;
	void dump_to_shared_file(m_array<complex, 2> field, std::array<int, 4> local_extent, std::array<int, 4> global_extent, std::string filename) const;
	std::unique_ptr<param_2d> param;
	m_array<complex, 2> e_x;
	m_array<complex, 2> e_y;
	m_array<complex, 2> e_z;
	m_array<complex, 2> b_x;
	m_array<complex, 2> b_y;
	m_array<complex, 2> b_z;
	m_array<double, 2> k_z;
	std::vector<double> k_x;
	std::vector<double> omega;
	bool fields_computed;
};
