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

#include "lander.h"

void autopilot (void)
  // Autopilot to adjust the engine throttle, parachute and attitude control
{
  // INSERT YOUR CODE HERE
	// Deploy parachute when appropriate
	if (parachute_status == NOT_DEPLOYED && safe_to_deploy_parachute() && (position.abs() - MARS_RADIUS) < 50000.0) {
		parachute_status = DEPLOYED;
	}

	// Orientation in opposite direction to velocity
	// De-orbit burn if periapsis outside atmosphere and attitude stabilisation enabled
	if (periapsis > 80000.0) {
		//double angle = acos((position*velocity) / (position.abs()*velocity.abs())) * (180 / M_PI);
		//orientation = vector3d(0.0, 0.0, -180 + angle);
		stabilized_attitude = true;
		throttle = 1;
		//if ((int)simulation_time % 10 == 0) { cout << angle << endl; }
		return;
	}

	double Kh = 0.01;
	double Kp = 0.5;
	double delta = 0.5;

	static vector3d radial = vector3d(1.0, 1.0, 1.0);
	radial = radial.norm();
	double sink_rate = -velocity * radial;
	double error = -(0.5 + Kh * (position.abs() - MARS_RADIUS) + sink_rate);
	double Pout = Kp * error;

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

}

void numerical_dynamics (void)
  // This is the function that performs the numerical integration to update the
  // lander's pose. The time step is delta_t (global variable).
{
  // INSERT YOUR CODE HERE
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
	gravity_force = ((-GRAVITY*MARS_MASS*lander_mass) / position.abs2())*position.norm();

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

void initialize_simulation (void)
  // Lander pose initialization - selects one of 10 possible scenarios
{
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
    autopilot_enabled = false;
    break;

  case 1:
    // a descent from rest at 10km altitude
    position = vector3d(0.0, -(MARS_RADIUS + 10000.0), 0.0);
    velocity = vector3d(0.0, 0.0, 0.0);
    orientation = vector3d(0.0, 0.0, 90.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = true;
    autopilot_enabled = false;
    break;

  case 2:
    // an elliptical polar orbit
    position = vector3d(0.0, 0.0, 1.2*MARS_RADIUS);
    velocity = vector3d(3500.0, 0.0, 0.0);
    orientation = vector3d(0.0, 0.0, 90.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;
    break;

  case 3:
    // polar surface launch at escape velocity (but drag prevents escape)
    position = vector3d(0.0, 0.0, MARS_RADIUS + LANDER_SIZE/2.0);
    velocity = vector3d(0.0, 0.0, 5027.0);
    orientation = vector3d(0.0, 0.0, 0.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;
    break;

  case 4:
    // an elliptical orbit that clips the atmosphere each time round, losing energy
    position = vector3d(0.0, 0.0, MARS_RADIUS + 100000.0);
    velocity = vector3d(4000.0, 0.0, 0.0);
    orientation = vector3d(0.0, 90.0, 0.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;
    break;

  case 5:
    // a descent from rest at the edge of the exosphere
    position = vector3d(0.0, -(MARS_RADIUS + EXOSPHERE), 0.0);
    velocity = vector3d(0.0, 0.0, 0.0);
    orientation = vector3d(0.0, 0.0, 90.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = true;
    autopilot_enabled = false;
    break;

  case 6:
	// areostationary circular orbit
	position = vector3d((MARS_RADIUS + AREOSTATIONARY_ALTITUDE), 0.0, 0.0);
	velocity = vector3d(0.0, sqrt((GRAVITY*MARS_MASS) / (MARS_RADIUS + AREOSTATIONARY_ALTITUDE)), 0.0);
	orientation = vector3d(0.0, 90.0, 0.0);
	delta_t = 0.1;
	parachute_status = NOT_DEPLOYED;
	stabilized_attitude = false;
	autopilot_enabled = false;
    break;

  case 7:
    break;

  case 8:
    break;

  case 9:
    break;

  }
}
