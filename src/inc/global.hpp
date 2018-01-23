#pragma once

#include <vector>
#include <array>
#include <complex>
#include <fftw3.h>
#include <boost/multi_array.hpp>

using complex = std::complex<double>;
using view_1d = boost::multi_array<complex, 2>::array_view<1>::type;
using view_2d = boost::multi_array<complex, 3>::array_view<2>::type;
using range = boost::multi_array_types::index_range;
template <typename T, int N> struct m_array {
	boost::multi_array<T, N> data;
	std::string name;
	bool calculated;
};

namespace constants {
	extern const double pi;
	extern const double c;
	extern const double epsilon_0;
	extern const complex imag_unit;
}

namespace fft {
	extern fftw_plan plan;
	extern std::vector<complex> data;
	void create_plan_1d(int dim_1, int sign);
	void create_plan_2d(int dim_1, int dim_2, int sign);
	std::vector<complex> execute_plan(std::vector<complex> in);
	void destroy_plan();
}

namespace tools {
	template <typename T, int N> void multiply_array(m_array<T, N>& field, const double scalar);
	std::vector<complex> array_to_vec(view_1d view);
	std::vector<complex> array_to_vec(view_2d view);
	void vec_to_array(view_1d& view, std::vector<complex> vec);
	void vec_to_array(view_2d& view, std::vector<complex> vec);
	std::vector<double> get_real(m_array<complex, 2> field, std::array<int, 2> size);
	std::vector<double> get_real(m_array<complex, 3> field, std::array<int, 3> size);
}
