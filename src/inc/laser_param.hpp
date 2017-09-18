#pragma once

class laser_param {
	friend class laser_bcs;
public:
	laser_param(double t_start, double t_end, double fwhm_time, double t_0, double x_0, double y_0, double omega, double amp, double w_0, int id);
	void set_values(double t_start, double t_end, double fwhm_time, double t_0, double x_0, double y_0, double omega, double amp, double w_0, int id);
	~laser_param() = default;
protected:
	double omega;
	double t_start;
	double t_end;
	double fwhm_time;
	double t_0;
	double x_0;
	double y_0;
	double amp;
	double w_0;
	int id;
};
