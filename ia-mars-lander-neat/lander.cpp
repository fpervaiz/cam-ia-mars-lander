// Mars lander simulator
// Version 1.10
// Mechanical simulation functions
// Gabor Csanyi and Andrew Gee, August 2017

// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation, to make use of it
// for non-commercial purposes, provided that (a) its original authorship
// is acknowledged and (b) no modified versions of the source code are
// published. Restriction (b) is designed to protect the integrity of the
// exercise for future generations of students. The authors would be happy
// to receive any suggested modifications by private correspondence to
// ahg@eng.cam.ac.uk and gc121@eng.cam.ac.uk.

#include <vector>
#include <time.h>
#include "lander.h"

void autopilot(void)
// Autopilot to adjust the engine throttle, parachute and attitude control
{
	// INSERT YOUR CODE HERE

	double altitude = position.abs() - MARS_RADIUS;

	// Deploy parachute when appropriate
	if (parachute_status == NOT_DEPLOYED && safe_to_deploy_parachute() && altitude < 50000.0) {
		parachute_status = DEPLOYED;
	}

	// De-orbit burn if periapsis outside atmosphere and attitude stabilisation enabled
	// Burn at or near apoapsis for better fuel usage
	bool near_apoapsis;
	if (abs(altitude - apoapsis) < apoapsis / 20) {
		near_apoapsis = true;
	}
	else {
		near_apoapsis = false;
	}
	if (periapsis > 80000.0 && near_apoapsis) {
		stabilized_attitude = true;
		throttle = 1;
		return;
	}

	/*
	switch (scenario) {

	case 0: // Working - 68.1 litres consumption
		Kh = 0.02;
		Kp = 0.5;
		break;

	case 1: // Working - 60.6 litres consumption
		Kh = 0.02;
		Kp = 0.5;
		break;

	case 2: // Working - 67.3 litres consumption
		Kh = 0.02;
		Kp = 0.5;
		break;

	case 3: // Working - 84.1 litres consumption
		Kh = 0.02;
		Kp = 0.05;
		break;

	case 4: // Working - 66.0 litres consumption
		Kh = 0.02;
		Kp = 0.5;
		break;

	case 5: // Working - 81.9 litres consumption
		Kh = 0.02;
		Kp = 0.5;
		break;

	case 6: // Working - 94.8 litres consumption
		Kh = 0.02;
		Kp = 0.5;
		break;
	}*/

	static vector<double> a_list, v_list;
	static bool data_written = false;

	double sink_rate = velocity * position.norm();
	double error = -(0.5 + Kh * altitude + sink_rate);
	double Pout = Kp * error;

	double lander_mass = UNLOADED_LANDER_MASS + FUEL_DENSITY * fuel * FUEL_CAPACITY;
	double lander_weight = GRAVITY * MARS_MASS * lander_mass / position.abs2();
	double density = atmospheric_density(position);
	double lander_drag_force = (0.5) * density * DRAG_COEF_LANDER * 3.14159*LANDER_SIZE*LANDER_SIZE * velocity.abs2() * (velocity.norm() * position.norm());
	double delta = (lander_weight - lander_drag_force) / MAX_THRUST; // delta equals resultant downward force on lander as proportion of max thrust

	if (delta > 1) {
		delta = 1;
	}
	else if (delta < 0) {
		delta = 0;
	}

	stabilized_attitude = true;

	if (Pout <= -delta) {
		throttle = 0;
	}
	else if (Pout > -delta && Pout < (1 - delta)) {
		throttle = delta + Pout;
	}
	else {
		throttle = 1;
	}

	a_list.push_back(altitude);
	v_list.push_back(sink_rate);

	/* Uncomment to enable log write
	if ((altitude < 1) && !(data_written)) {
		ofstream fout;
		fout.open("profile.txt");
		for (int i = 0; i < a_list.size(); i = i + 1) {
			fout << a_list[i] << ' ' << v_list[i] << endl;
		}
		data_written = true;
	}*/

}

void numerical_dynamics(void)
// This is the function that performs the numerical integration to update the
// lander's pose. The time step is delta_t (global variable).
{
	// INSERT YOUR CODE HERE
	// Typical Martian wind speeds: https://sciencing.com/average-wind-speed-mars-3805.html
	static vector3d surface_wind, wind;
	static bool gusting = false, steady = false;
	static double steady_duration, gust_duration;

	if (wind_enabled) {
		if (!gusting && !steady) {
			if (rand() % 6 == 0) {
				// Gusting 1 in 6 times (arbitrary)
				surface_wind = vector3d(rand() % 40 - 20, rand() % 40 - 20, rand() % 40 - 20);
				gust_duration = rand() % 4 + 5;
				gust_start_time = simulation_time;
				gusting = true;
			}
			else {
				// Steady
				surface_wind = vector3d(rand() % 12 - 6, rand() % 12 - 6, rand() % 12 - 6);
				steady_duration = rand() % 120 + 15;
				steady_start_time = simulation_time;
				steady = true;
			}
		}

		// Reset steady/gusting when duration elapsed
		if (gusting && simulation_time > gust_start_time + gust_duration) {
			gusting = false;
		}
		else if (steady && simulation_time > steady_start_time + steady_duration) {
			steady = false;
		}

		// Wind strengh variation with altitude (peak 3x)
		// Roughly based on https://www.researchgate.net/profile/Alexey_Pankine/publication/252307634/figure/fig3/AS:298062665797634@1448075085348/Mars-atmospheric-wind-profile.png
		double altitude = position.abs() - MARS_RADIUS;
		if (altitude < 16000) {
			wind = surface_wind * (1 + altitude / 6400);
		}
		else if (altitude >= 16000 && altitude < 36000) {
			wind = surface_wind * (12 - altitude / 6000);
		}
		else {
			wind = surface_wind * 0;
		}

		wind_at_pos = wind.abs();
		wind_at_surface = surface_wind.abs();

		velocity = velocity - wind;
	}

	vector3d rocket_thrust_force, lander_drag_force, chute_drag_force, gravity_force, resultant_force, acceleration, new_position;
	static vector3d previous_position;
	double density, lander_mass;

	rocket_thrust_force = thrust_wrt_world();
	density = atmospheric_density(position);
	lander_mass = UNLOADED_LANDER_MASS + fuel * FUEL_CAPACITY;

	lander_drag_force = (-0.5*density*DRAG_COEF_LANDER*3.14159*LANDER_SIZE*LANDER_SIZE*velocity.abs2())*velocity.norm();

	if (parachute_status == DEPLOYED) {
		double parachute_area = (2.0*LANDER_SIZE)*(2.0*LANDER_SIZE) * 5;
		chute_drag_force = (-0.5*density*DRAG_COEF_CHUTE*parachute_area*velocity.abs2())*velocity.norm();
	}
	else {
		chute_drag_force = vector3d(0.0, 0.0, 0.0);
	}
	gravity_force = ((-GRAVITY * MARS_MASS*lander_mass) / position.abs2())*position.norm();

	resultant_force = rocket_thrust_force + lander_drag_force + chute_drag_force + gravity_force;

	// Numerical integration
	acceleration = resultant_force / lander_mass;
	if (simulation_time == 0.0) {
		// Euler method
		new_position = position + delta_t * velocity;
		velocity = velocity + delta_t * acceleration;
	}
	else {
		// Verlet method
		new_position = 2 * position - previous_position + delta_t * delta_t*acceleration;
		velocity = (new_position - position) / delta_t;
	}

	previous_position = position;
	position = new_position;

	// Periapsis calculation
	vector3d h = position ^ velocity;
	vector3d node = vector3d(0, 0, 1) ^ h;
	double mu = GRAVITY * MARS_MASS;
	static int counter = 0;

	vector3d eccentricity = (1 / mu)*((velocity.abs2() - (mu / position.abs()))*position - (position*velocity)*velocity);
	double e = eccentricity.abs();
	double spec_mech_energy = velocity.abs2() / 2 - mu / position.abs();
	double a = -mu / (2 * spec_mech_energy);
	periapsis = ((1 - e)*a) - MARS_RADIUS;
	apoapsis = ((1 + e)*a) - MARS_RADIUS;

	// Here we can apply an autopilot to adjust the thrust, parachute and attitude
	if (autopilot_enabled) autopilot();

	// Here we can apply 3-axis stabilization to ensure the base is always pointing downwards
	if (stabilized_attitude) attitude_stabilization();
}

void initialize_simulation(void)
// Lander pose initialization - selects one of 10 possible scenarios
{
	srand(time(NULL));
	
	// The parameters to set are:
	// position - in Cartesian planetary coordinate system (m)
	// velocity - in Cartesian planetary coordinate system (m/s)
	// orientation - in lander coordinate system (xyz Euler angles, degrees)
	// delta_t - the simulation time step
	// boolean state variables - parachute_status, stabilized_attitude, autopilot_enabled
	// scenario_description - a descriptive string for the help screen

	scenario_description[0] = "circular orbit";
	scenario_description[1] = "descent from 10km";
	scenario_description[2] = "elliptical orbit, thrust changes orbital plane";
	scenario_description[3] = "polar launch at escape velocity (but drag prevents escape)";
	scenario_description[4] = "elliptical orbit that clips the atmosphere and decays";
	scenario_description[5] = "descent from 200km";
	scenario_description[6] = "areostationary circular orbit";
	scenario_description[7] = "";
	scenario_description[8] = "";
	scenario_description[9] = "";

	switch (scenario) {

	case 0:
		// a circular equatorial orbit
		position = vector3d(1.2*MARS_RADIUS, 0.0, 0.0);
		velocity = vector3d(0.0, -3247.087385863725, 0.0);
		orientation = vector3d(0.0, 90.0, 0.0);
		delta_t = 0.1;
		parachute_status = NOT_DEPLOYED;
		stabilized_attitude = false;
		autopilot_enabled = true;
		break;

	case 1:
		// a descent from rest at 10km altitude
		position = vector3d(0.0, -(MARS_RADIUS + 10000.0), 0.0);
		velocity = vector3d(0.0, 0.0, 0.0);
		orientation = vector3d(0.0, 0.0, 90.0);
		delta_t = 0.1;
		parachute_status = NOT_DEPLOYED;
		stabilized_attitude = false;
		autopilot_enabled = true;
		break;

	case 2:
		// an elliptical polar orbit
		position = vector3d(0.0, 0.0, 1.2*MARS_RADIUS);
		velocity = vector3d(3500.0, 0.0, 0.0);
		orientation = vector3d(0.0, 0.0, 90.0);
		delta_t = 0.1;
		parachute_status = NOT_DEPLOYED;
		stabilized_attitude = false;
		autopilot_enabled = true;
		break;

	case 3:
		// polar surface launch at escape velocity (but drag prevents escape)
		position = vector3d(0.0, 0.0, MARS_RADIUS + LANDER_SIZE / 2.0);
		velocity = vector3d(0.0, 0.0, 5027.0);
		orientation = vector3d(0.0, 0.0, 0.0);
		delta_t = 0.1;
		parachute_status = NOT_DEPLOYED;
		stabilized_attitude = false;
		autopilot_enabled = true;
		break;

	case 4:
		// an elliptical orbit that clips the atmosphere each time round, losing energy
		position = vector3d(0.0, 0.0, MARS_RADIUS + 100000.0);
		velocity = vector3d(4000.0, 0.0, 0.0);
		orientation = vector3d(0.0, 90.0, 0.0);
		delta_t = 0.1;
		parachute_status = NOT_DEPLOYED;
		stabilized_attitude = false;
		autopilot_enabled = true;
		break;

	case 5:
		// a descent from rest at the edge of the exosphere
		position = vector3d(0.0, -(MARS_RADIUS + EXOSPHERE), 0.0);
		velocity = vector3d(0.0, 0.0, 0.0);
		orientation = vector3d(0.0, 0.0, 90.0);
		delta_t = 0.1;
		parachute_status = NOT_DEPLOYED;
		stabilized_attitude = true;
		autopilot_enabled = true;
		break;

	case 6:
		// areostationary circular orbit
		position = vector3d((MARS_RADIUS + AREOSTATIONARY_ALTITUDE), 0.0, 0.0);
		velocity = vector3d(0.0, sqrt((GRAVITY*MARS_MASS) / (MARS_RADIUS + AREOSTATIONARY_ALTITUDE)), 0.0);
		orientation = vector3d(0.0, 90.0, 0.0);
		delta_t = 0.1;
		parachute_status = NOT_DEPLOYED;
		stabilized_attitude = false;
		autopilot_enabled = true;
		break;

	case 7:
		break;

	case 8:
		break;

	case 9:
		break;

	}
}
