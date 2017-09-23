#include <iostream>

#include "inc/laser_param.hpp"
#include "inc/global.hpp"

laser_param::laser_param(double t_start, double t_end, double fwhm_time, double t_0, double x_0, double y_0, double omega, double amplitude, double w_0, int direction, int id) {
	this->set_values(t_start, t_end, fwhm_time, t_0, x_0, y_0, omega, amplitude, w_0, direction, id);
}

void laser_param::set_values(double t_start, double t_end, double fwhm_time, double t_0, double x_0, double y_0, double omega, double amplitude, double w_0, int direction, int id) {
	if(t_start > t_end || fwhm_time < 0 || omega < 0 || amplitude < 0 || w_0 < 2.0 * constants::c / omega) {
		std::cerr << "error: bad value" << std::endl;
		return;
	}
	this->t_start = t_start;
	this->t_end = t_end;
	this->fwhm_time = fwhm_time;
	this->t_0 = t_0;
	this->x_0 = x_0;
	this->y_0 = y_0;
	this->omega = omega;
	this->amplitude = amplitude;
	this->w_0 = w_0;
  this->direction = direction;
	this->id = id;
}
