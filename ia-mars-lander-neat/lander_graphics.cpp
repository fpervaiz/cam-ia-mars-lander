// Mars lander simulator
// Version 1.10
// Graphics functions
// Gabor Csanyi and Andrew Gee, August 2017

// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation, to make use of it
// for non-commercial purposes, provided that (a) its original authorship
// is acknowledged and (b) no modified versions of the source code are
// published. Restriction (b) is designed to protect the integrity of the
// exercise for future generations of students. The authors would be happy
// to receive any suggested modifications by private correspondence to
// ahg@eng.cam.ac.uk and gc121@eng.cam.ac.uk.

// Some functions adapted from freeglut_geometry.c, which is covered by the
// following license:
//
// Copyright (c) 1999-2000 Pawel W. Olszta. All Rights Reserved.
// Written by Pawel W. Olszta, <olszta@sourceforge.net>
// Creation date: Fri Dec 3 1999
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
// PAWEL W. OLSZTA BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
// IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
// CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// Some functions adapted from trackball.cpp by Gavin Bell, which is covered by
// the following license:
//
// (c) Copyright 1993, 1994, Silicon Graphics, Inc.
// ALL RIGHTS RESERVED
// Permission to use, copy, modify, and distribute this software for
// any purpose and without fee is hereby granted, provided that the above
// copyright notice appear in all copies and that both the copyright notice
// and this permission notice appear in supporting documentation, and that
// the name of Silicon Graphics, Inc. not be used in advertising
// or publicity pertaining to distribution of the software without specific,
// written prior permission.
//
// THE MATERIAL EMBODIED ON THIS SOFTWARE IS PROVIDED TO YOU "AS-IS"
// AND WITHOUT WARRANTY OF ANY KIND, EXPRESS, IMPLIED OR OTHERWISE,
// INCLUDING WITHOUT LIMITATION, ANY WARRANTY OF MERCHANTABILITY OR
// FITNESS FOR A PARTICULAR PURPOSE.  IN NO EVENT SHALL SILICON
// GRAPHICS, INC.  BE LIABLE TO YOU OR ANYONE ELSE FOR ANY DIRECT,
// SPECIAL, INCIDENTAL, INDIRECT OR CONSEQUENTIAL DAMAGES OF ANY
// KIND, OR ANY DAMAGES WHATSOEVER, INCLUDING WITHOUT LIMITATION,
// LOSS OF PROFIT, LOSS OF USE, SAVINGS OR REVENUE, OR THE CLAIMS OF
// THIRD PARTIES, WHETHER OR NOT SILICON GRAPHICS, INC.  HAS BEEN
// ADVISED OF THE POSSIBILITY OF SUCH LOSS, HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, ARISING OUT OF OR IN CONNECTION WITH THE
// POSSESSION, USE OR PERFORMANCE OF THIS SOFTWARE.
//
// US Government Users Restricted Rights
// Use, duplication, or disclosure by the Government is subject to
// restrictions set forth in FAR 52.227.19(c)(2) or subparagraph
// (c)(1)(ii) of the Rights in Technical Data and Computer Software
// clause at DFARS 252.227-7013 and/or in similar or successor
// clauses in the FAR or the DOD or NASA FAR Supplement.
// Unpublished-- rights reserved under the copyright laws of the
// United States.  Contractor/manufacturer is Silicon Graphics,
// Inc., 2011 N.  Shoreline Blvd., Mountain View, CA 94039-7311.
//
// OpenGL(TM) is a trademark of Silicon Graphics, Inc.

#define DECLARE_GLOBAL_VARIABLES
#include "lander.h"

#include <stdlib.h>
#include <string>

void invert(double m[], double mout[])
// Inverts a 4x4 OpenGL rotation matrix
{
	double zero_three, one_three, two_three;
	zero_three = -m[12] * m[0] - m[13] * m[1] - m[14] * m[2];
	one_three = -m[12] * m[4] - m[13] * m[5] - m[14] * m[6];
	two_three = -m[12] * m[8] - m[13] * m[9] - m[14] * m[10];
	mout[1] = m[4]; mout[4] = m[1]; mout[2] = m[8]; mout[8] = m[2];
	mout[6] = m[9]; mout[9] = m[6]; mout[12] = zero_three; mout[13] = one_three;
	mout[14] = two_three; mout[0] = m[0]; mout[5] = m[5]; mout[10] = m[10];
	mout[15] = 1.0; mout[3] = 0.0; mout[7] = 0.0; mout[11] = 0.0;
}

void xyz_euler_to_matrix(vector3d ang, double m[])
// Constructs a 4x4 OpenGL rotation matrix from xyz Euler angles
{
	double sin_a, sin_b, sin_g, cos_a, cos_b, cos_g;
	double ra, rb, rg;

	// Pre-calculate radian angles
	ra = ang.x*M_PI / (double)180;
	rb = ang.y*M_PI / (double)180;
	rg = ang.z*M_PI / (double)180;

	// Pre-calculate sines and cosines
	cos_a = cos(ra);
	cos_b = cos(rb);
	cos_g = cos(rg);
	sin_a = sin(ra);
	sin_b = sin(rb);
	sin_g = sin(rg);

	// Create the correct matrix coefficients
	m[0] = cos_a * cos_b;
	m[1] = sin_a * cos_b;
	m[2] = -sin_b;
	m[3] = 0.0;
	m[4] = cos_a * sin_b * sin_g - sin_a * cos_g;
	m[5] = sin_a * sin_b * sin_g + cos_a * cos_g;
	m[6] = cos_b * sin_g;
	m[7] = 0.0;
	m[8] = cos_a * sin_b * cos_g + sin_a * sin_g;
	m[9] = sin_a * sin_b * cos_g - cos_a * sin_g;
	m[10] = cos_b * cos_g;
	m[11] = 0.0;
	m[12] = 0.0;
	m[13] = 0.0;
	m[14] = 0.0;
	m[15] = 1.0;
}

vector3d matrix_to_xyz_euler(double m[])
// Decomposes a 4x4 OpenGL rotation matrix into xyz Euler angles
{
	double tmp;
	vector3d ang;

	// Catch degenerate elevation cases
	if (m[2] < -0.99999999) {
		ang.y = 90.0;
		ang.x = 0.0;
		ang.z = acos(m[8]);
		if ((sin(ang.z) > 0.0) ^ (m[4] > 0.0)) ang.z = -ang.z;
		ang.z *= 180.0 / M_PI;
		return ang;
	}
	if (m[2] > 0.99999999) {
		ang.y = -90.0;
		ang.x = 0.0;
		ang.z = acos(m[5]);
		if ((sin(ang.z) < 0.0) ^ (m[4] > 0.0)) ang.z = -ang.z;
		ang.z *= 180.0 / M_PI;
		return ang;
	}

	// Non-degenerate elevation - between -90 and +90
	ang.y = asin(-m[2]);

	// Now work out azimuth - between -180 and +180
	tmp = m[0] / cos(ang.y); // the denominator will not be zero
	if (tmp <= -1.0) ang.x = M_PI;
	else if (tmp >= 1.0) ang.x = 0.0;
	else ang.x = acos(tmp);
	if (((sin(ang.x) * cos(ang.y)) > 0.0) ^ ((m[1]) > 0.0)) ang.x = -ang.x;

	// Now work out roll - between -180 and +180
	tmp = m[10] / cos(ang.y); // the denominator will not be zero
	if (tmp <= -1.0) ang.z = M_PI;
	else if (tmp >= 1.0) ang.z = 0.0;
	else ang.z = acos(tmp);
	if (((sin(ang.z) * cos(ang.y)) > 0.0) ^ ((m[6]) > 0.0)) ang.z = -ang.z;

	// Convert to degrees
	ang.y *= 180.0 / M_PI;
	ang.x *= 180.0 / M_PI;
	ang.z *= 180.0 / M_PI;

	return ang;
}

void normalize_quat(quat_t &q)
// Normalizes a quaternion
{
	double mag;
	mag = (q.v.x*q.v.x + q.v.y*q.v.y + q.v.z*q.v.z + q.s*q.s);
	if (mag > 0.0) {
		q.v.x /= mag; q.v.y /= mag;
		q.v.z /= mag; q.s /= mag;
	}
}

quat_t axis_to_quat(vector3d a, const double phi)
// Given an axis and angle, compute quaternion
{
	quat_t q;
	q.v = a.norm() * sin(phi / 2.0);
	q.s = cos(phi / 2.0);
	return q;
}

double project_to_sphere(const double r, const double x, const double y)
// Project an x,y pair onto a sphere of radius r or a hyperbolic sheet if
// we are away from the centre of the sphere
{
	double d, t, z;

	d = sqrt(x*x + y * y);
	if (d < (r * 0.70710678118654752440)) z = sqrt(r*r - d * d);
	else { // on hyperbola
		t = r / 1.41421356237309504880;
		z = t * t / d;
	}
	return z;
}

quat_t add_quats(quat_t q1, quat_t q2)
// Given two rotations q1 and q2, calculate single equivalent quaternion
{
	quat_t s;
	vector3d t1 = q1.v * q2.s;
	vector3d t2 = q2.v * q1.s;
	vector3d t3 = q2.v ^ q1.v;
	vector3d tf = t1 + t2;
	s.v = tf + t3;
	s.s = q1.s * q2.s - q1.v * q2.v;
	normalize_quat(s);
	return s;
}

void quat_to_matrix(double m[], const quat_t q)
// Convert quaternion into a rotation matrix
{
	m[0] = 1.0 - 2.0 * (q.v.y * q.v.y + q.v.z * q.v.z);
	m[1] = 2.0 * (q.v.x * q.v.y - q.v.z * q.s);
	m[2] = 2.0 * (q.v.z * q.v.x + q.v.y * q.s);
	m[3] = 0.0;

	m[4] = 2.0 * (q.v.x * q.v.y + q.v.z * q.s);
	m[5] = 1.0 - 2.0 * (q.v.z * q.v.z + q.v.x * q.v.x);
	m[6] = 2.0 * (q.v.y * q.v.z - q.v.x * q.s);
	m[7] = 0.0;

	m[8] = 2.0 * (q.v.z * q.v.x - q.v.y * q.s);
	m[9] = 2.0 * (q.v.y * q.v.z + q.v.x * q.s);
	m[10] = 1.0 - 2.0 * (q.v.y * q.v.y + q.v.x * q.v.x);
	m[11] = 0.0;

	m[12] = 0.0; m[13] = 0.0; m[14] = 0.0; m[15] = 1.0;
}

quat_t track_quats(const double p1x, const double p1y, const double p2x, const double p2y)
// Derive quaternion from x and y mouse displacements
{
	double t, phi;
	vector3d a, p1, p2, d;
	quat_t q;

	if ((p1x == p2x) && (p1y == p2y)) {
		q.v.x = 0.0; q.v.y = 0.0; q.v.z = 0.0; q.s = 1.0;
		return q;
	}

	p1.x = p1x; p1.y = p1y;
	p1.z = project_to_sphere(0.5, p1x, p1y);
	p2.x = p2x; p2.y = p2y;
	p2.z = project_to_sphere(0.5, p2x, p2y);
	a = p2 ^ p1; d = p1 - p2; t = d.abs();
	if (t > 1.0) t = 1.0;
	if (t < -1.0) t = -1.0;
	phi = 2.0 * asin(t);
	return axis_to_quat(a, phi);
}

void microsecond_time(unsigned long long &t)
// Returns system time in microseconds
{
#ifdef _WIN32
	LARGE_INTEGER counter, frequency;
	QueryPerformanceFrequency(&frequency);
	QueryPerformanceCounter(&counter);
	counter.QuadPart *= 1000000;
	t = (unsigned long long)(counter.QuadPart / frequency.QuadPart);
#else
	struct timeval tv;
	gettimeofday(&tv, NULL);
	t = (unsigned long long)tv.tv_usec + 1000000 * (unsigned long long)tv.tv_sec;
#endif
}

void fghCircleTable(double **sint, double **cost, const int n)
// Borrowed from freeglut source code, used to draw hemispheres and open cones
{
	int i;
	const int size = abs(n);
	const double angle = 2 * M_PI / (double)((n == 0) ? 1 : n);

	*sint = (double*)calloc(sizeof(double), size + 1);
	*cost = (double*)calloc(sizeof(double), size + 1);
	if (!(*sint) || !(*cost)) exit(1);

	(*sint)[0] = 0.0;
	(*cost)[0] = 1.0;

	for (i = 1; i < size; i++) {
		(*sint)[i] = sin(angle*i);
		(*cost)[i] = cos(angle*i);
	}

	(*sint)[size] = (*sint)[0];
	(*cost)[size] = (*cost)[0];
}

double atmospheric_density(vector3d pos)
// Simple exponential model between surface and exosphere (around 200km), surface density is approximately 0.017 kg/m^3,
// scale height is approximately 11km
{
	double alt;

	alt = pos.abs() - MARS_RADIUS;
	if ((alt > EXOSPHERE) || (alt < 0.0)) return 0.0;
	else return (0.017 * exp(-alt / 11000.0));
}

void update_closeup_coords(void)
// Updates the close-up view's coordinate frame, based on the lander's current position and velocity.
// This needs to be called every time step, even if the view is not being rendered, since any-angle
// attitude stabilizers reference closeup_coords.right
{
	vector3d s, tv, t;
	double tmp;

	// Direction from surface to lander (radial) - this must map to the world y-axis
	s = position.norm();

	// Direction of tangential velocity - this must map to the world x-axis
	tv = velocity_from_positions - (velocity_from_positions*s)*s;
	if (tv.abs() < SMALL_NUM) // vertical motion only, use last recorded tangential velocity
		tv = closeup_coords.backwards ? (closeup_coords.right*s)*s - closeup_coords.right : closeup_coords.right - (closeup_coords.right*s)*s;
	if (tv.abs() > SMALL_NUM) t = tv.norm();

	// Check these two vectors are non-zero and perpendicular (they should be, unless s and closeup_coords.right happen to be parallel)
	if ((tv.abs() <= SMALL_NUM) || (fabs(s*t) > SMALL_NUM)) {
		// Set t to something perpendicular to s
		t.x = -s.y; t.y = s.x; t.z = 0.0;
		if (t.abs() < SMALL_NUM) { t.x = -s.z; t.y = 0.0; t.z = s.x; }
		t = t.norm();
	}

	// Adjust the terrain texture angle if the lander has changed direction. The motion will still be along
	// the x-axis, so we need to rotate the texture to compensate.
	if (closeup_coords.initialized) {
		if (closeup_coords.backwards) {
			tmp = -closeup_coords.right*t;
			if (tmp > 1.0) tmp = 1.0; if (tmp < -1.0) tmp = -1.0;
			if ((-closeup_coords.right^t)*position.norm() < 0.0) terrain_angle += (180.0 / M_PI)*acos(tmp);
			else terrain_angle -= (180.0 / M_PI)*acos(tmp);
		}
		else {
			tmp = closeup_coords.right*t;
			if (tmp > 1.0) tmp = 1.0; if (tmp < -1.0) tmp = -1.0;
			if ((closeup_coords.right^t)*position.norm() < 0.0) terrain_angle += (180.0 / M_PI)*acos(tmp);
			else terrain_angle -= (180.0 / M_PI)*acos(tmp);
		}
		while (terrain_angle < 0.0) terrain_angle += 360.0;
		while (terrain_angle >= 360.0) terrain_angle -= 360.0;
	}

	// Normally we maintain motion to the right, the one exception being when the ground speed passes
	// through zero and changes sign. A sudden 180 degree change of viewpoint would be confusing, so
	// in this instance we allow the lander to fly to the left.
	if (closeup_coords.initialized && (closeup_coords.right*t < 0.0)) {
		closeup_coords.backwards = true;
		closeup_coords.right = -1.0*t;
	}
	else {
		closeup_coords.backwards = false;
		closeup_coords.right = t;
		closeup_coords.initialized = true;
	}
}

bool safe_to_deploy_parachute(void)
// Checks whether the parachute is safe to deploy at the current position and velocity
{
	double drag;

	// Assume high Reynolds number, quadratic drag = -0.5 * rho * v^2 * A * C_d
	drag = 0.5*DRAG_COEF_CHUTE*atmospheric_density(position)*5.0*2.0*LANDER_SIZE*2.0*LANDER_SIZE*velocity_from_positions.abs2();
	// Do not use the global variable "altitude" here, in case this function is called from within the
	// numerical_dynamics function, before altitude is updated in the update_visualization function
	if ((drag > MAX_PARACHUTE_DRAG) || ((velocity_from_positions.abs() > MAX_PARACHUTE_SPEED) && ((position.abs() - MARS_RADIUS) < EXOSPHERE))) return false;
	else return true;
}

void update_visualization(void)
// The visualization part of the idle function. Re-estimates altitude, velocity, climb speed and ground
// speed from current and previous positions. Updates throttle and fuel levels, then redraws all subwindows.
{
	static vector3d last_track_position;
	vector3d av_p, d;
	double a, b, c, mu;

	simulation_time += delta_t;
	altitude = position.abs() - MARS_RADIUS;

	// Use average of current and previous positions when calculating climb and ground speeds
	av_p = (position + last_position).norm();
	if (delta_t != 0.0) velocity_from_positions = (position - last_position) / delta_t;
	else velocity_from_positions = vector3d(0.0, 0.0, 0.0);
	climb_speed = velocity_from_positions * av_p;
	ground_speed = (velocity_from_positions - climb_speed * av_p).abs();

	// Check to see whether the lander has landed
	if (altitude < LANDER_SIZE / 2.0) {
		// Estimate position and time of impact
		d = position - last_position;
		a = d.abs2();
		b = 2.0*last_position*d;
		c = last_position.abs2() - (MARS_RADIUS + LANDER_SIZE / 2.0) * (MARS_RADIUS + LANDER_SIZE / 2.0);
		mu = (-b - sqrt(b*b - 4.0*a*c)) / (2.0*a);
		position = last_position + mu * d;
		simulation_time -= (1.0 - mu)*delta_t;
		altitude = LANDER_SIZE / 2.0;
		landed = true;
		if ((fabs(climb_speed) > MAX_IMPACT_DESCENT_RATE) || (fabs(ground_speed) > MAX_IMPACT_GROUND_SPEED)) crashed = true;
		velocity_from_positions = vector3d(0.0, 0.0, 0.0);
	}

	// Update throttle and fuel (throttle might have been adjusted by the autopilot)
	if (throttle < 0.0) throttle = 0.0;
	if (throttle > 1.0) throttle = 1.0;
	fuel -= delta_t * (FUEL_RATE_AT_MAX_THRUST*throttle) / FUEL_CAPACITY;
	if (fuel <= 0.0) fuel = 0.0;
	if (landed || (fuel == 0.0)) throttle = 0.0;
	throttle_control = (short)(throttle*THROTTLE_GRANULARITY + 0.5);

	// Check to see whether the parachute has vaporized or the tethers have snapped
	if (parachute_status == DEPLOYED) {
		if (!safe_to_deploy_parachute() || parachute_lost) {
			parachute_lost = true; // to guard against the autopilot reinstating the parachute!
			parachute_status = LOST;
		}
	}

	// Update record of lander's previous positions, but only if the position or the velocity has 
	// changed significantly since the last update
	if (!track.n || (position - last_track_position).norm() * velocity_from_positions.norm() < TRACK_ANGLE_DELTA
		|| (position - last_track_position).abs() > TRACK_DISTANCE_DELTA) {
		track.pos[track.p] = position;
		track.n++; if (track.n > N_TRACK) track.n = N_TRACK;
		track.p++; if (track.p == N_TRACK) track.p = 0;
		last_track_position = position;
	}
}

void attitude_stabilization(void)
// Three-axis stabilization to ensure the lander's base is always pointing downwards 
{
	vector3d up, left, out;
	double m[16];

	up = velocity.norm() * -1; // this is the direction we want the lander's nose to point in
	stabilized_attitude_angle = 2;

	// !!!!!!!!!!!!! HINT TO STUDENTS ATTEMPTING THE EXTENSION EXERCISES !!!!!!!!!!!!!!
	// For any-angle attitude control, we just need to set "up" to something different,
	// and leave the remainder of this function unchanged. For example, suppose we want
	// the attitude to be stabilized at stabilized_attitude_angle to the vertical in the
	// close-up view. So we need to rotate "up" by stabilized_attitude_angle degrees around
	// an axis perpendicular to the plane of the close-up view. This axis is given by the
	// vector product of "up"and "closeup_coords.right". To calculate the result of the
	// rotation, search the internet for information on the axis-angle rotation formula.

	// Set left to something perpendicular to up
	left.x = -up.y; left.y = up.x; left.z = 0.0;
	if (left.abs() < SMALL_NUM) { left.x = -up.z; left.y = 0.0; left.z = up.x; }
	left = left.norm();
	out = left ^ up;
	// Construct modelling matrix (rotation only) from these three vectors
	m[0] = out.x; m[1] = out.y; m[2] = out.z; m[3] = 0.0;
	m[4] = left.x; m[5] = left.y; m[6] = left.z; m[7] = 0.0;
	m[8] = up.x; m[9] = up.y; m[10] = up.z; m[11] = 0.0;
	m[12] = 0.0; m[13] = 0.0; m[14] = 0.0; m[15] = 1.0;
	// Decomponse into xyz Euler angles
	orientation = matrix_to_xyz_euler(m);
}

vector3d thrust_wrt_world(void)
// Works out thrust vector in the world reference frame, given the lander's orientation
{
	double m[16], k, delayed_throttle, lag = ENGINE_LAG;
	vector3d a, b;
	static double lagged_throttle = 0.0;
	static double last_time_lag_updated = -1.0;

	if (simulation_time < last_time_lag_updated) lagged_throttle = 0.0; // simulation restarted
	if (throttle < 0.0) throttle = 0.0;
	if (throttle > 1.0) throttle = 1.0;
	if (landed || (fuel == 0.0)) throttle = 0.0;

	if (simulation_time != last_time_lag_updated) {

		// Delayed throttle value from the throttle history buffer
		if (throttle_buffer_length > 0) {
			delayed_throttle = throttle_buffer[throttle_buffer_pointer];
			throttle_buffer[throttle_buffer_pointer] = throttle;
			throttle_buffer_pointer = (throttle_buffer_pointer + 1) % throttle_buffer_length;
		}
		else delayed_throttle = throttle;

		// Lag, with time constant ENGINE_LAG
		if (lag <= 0.0) k = 0.0;
		else k = pow(exp(-1.0), delta_t / lag);
		lagged_throttle = k * lagged_throttle + (1.0 - k)*delayed_throttle;

		last_time_lag_updated = simulation_time;
	}

	if (stabilized_attitude && (stabilized_attitude_angle == 0)) { // specific solution, avoids rounding errors in the more general calculation below
		b = lagged_throttle * MAX_THRUST*position.norm();
	}
	else {
		a.x = 0.0; a.y = 0.0; a.z = lagged_throttle * MAX_THRUST;
		xyz_euler_to_matrix(orientation, m);
		b.x = m[0] * a.x + m[4] * a.y + m[8] * a.z;
		b.y = m[1] * a.x + m[5] * a.y + m[9] * a.z;
		b.z = m[2] * a.x + m[6] * a.y + m[10] * a.z;
	}
	return b;
}

void update_lander_state(void)
// The GLUT idle function, called every time round the event loop
{

	// This needs to be called every time step, even if the close-up view is not being rendered,
	// since any-angle attitude stabilizers reference closeup_coords.right
	update_closeup_coords();

	// Update historical record
	last_position = position;

	// Mechanical dynamics
	numerical_dynamics();

	// Refresh the visualization
	update_visualization();
}

void reset_simulation(void)
// Resets the simulation to the initial state
{
	vector3d p, tv;
	unsigned long i;

	// Reset these three lander parameters here, so they can be overwritten in initialize_simulation() if so desired
	stabilized_attitude_angle = 0;
	throttle = 0.0;
	fuel = 1.0;

	// Restore initial lander state
	initialize_simulation();

	// Check whether the lander is underground - if so, make sure it doesn't move anywhere
	landed = false;
	crashed = false;
	altitude = position.abs() - MARS_RADIUS;
	if (altitude < LANDER_SIZE / 2.0) {
		landed = true;
		velocity = vector3d(0.0, 0.0, 0.0);
	}

	if (wind_enabled) steady_start_time = 0, gust_start_time = 0;

	// Visualisation routine's record of various speeds and velocities
	velocity_from_positions = velocity;
	last_position = position - delta_t * velocity_from_positions;
	p = position.norm();
	climb_speed = velocity_from_positions * p;
	tv = velocity_from_positions - climb_speed * p;
	ground_speed = tv.abs();

	// Miscellaneous state variables
	throttle_control = (short)(throttle*THROTTLE_GRANULARITY + 0.5);
	simulation_time = 0.0;
	track.n = 0;
	parachute_lost = false;
	closeup_coords.initialized = false;
	closeup_coords.backwards = false;
	closeup_coords.right = vector3d(1.0, 0.0, 0.0);
	update_closeup_coords();

	// Initialize the throttle history buffer
	if (delta_t > 0.0) throttle_buffer_length = (unsigned long)(ENGINE_DELAY / delta_t + 0.5);
	else throttle_buffer_length = 0;
	if (throttle_buffer_length > 0) {
		if (throttle_buffer != NULL) delete[] throttle_buffer;
		throttle_buffer = new double[throttle_buffer_length];
		for (i = 0; i < throttle_buffer_length; i++) throttle_buffer[i] = throttle;
		throttle_buffer_pointer = 0;
	}
}

int main(int argc, char* argv[])
// Initializes GLUT windows and lander state, then enters GLUT main loop
{
	if (argc < 4) {
		cout << "Usage: lander <scenario number> <Kh value> <Kp value> <parachute deploy altitude> <-r readable output>";
		return 0;
	}

	bool readable;

	if (argc == 5) {
		readable = true;
	}
	else {
		readable = false;
	}
	wind_enabled = true;
	// Initialise the simulation state
	scenario = atoi(argv[1]);
	Kh = atof(argv[2]);
	Kp = atof(argv[3]);
	if (readable) {
		cout << "\nRunning scenario " << scenario << " with Kh " << Kh << " and Kp " << Kp;
	}
	reset_simulation();
	microsecond_time(time_program_started);

	double simulation_time_limit = 100000;

	while (!landed && !crashed && simulation_time < simulation_time_limit) {
		update_lander_state();
	}


	if (readable) {

		cout << "\n\n";

		if (crashed) {
			cout << "Crashed";
		}
		else if (landed) {
			cout << "Landed";
		}
		else {
			cout << "Timed out";
		}

		cout << "\nTime " << simulation_time << " s";
		cout << "\nFuel used: " << FUEL_CAPACITY * (1.0 - fuel) << " l";
		cout << "\nLanding vertical speed " << -climb_speed << " m/s";
		cout << "\nLanding ground speed " << ground_speed << " m/s";
		cout << "\n ";
	}
	else {
		if (crashed) {
			cout << "0,";
		}
		else if (landed) {
			cout << "1,";
		}
		else {
			cout << "-1,";
		}
		cout << simulation_time << "," << FUEL_CAPACITY * (1.0 - fuel) << "," << -climb_speed << "," << ground_speed;
	}
	
	
}
