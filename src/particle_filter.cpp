/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h>
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	// Set number of particles
	num_particles = 100;
	// Create a default random engine
	default_random_engine gen;
	// (Gaussian) distribution for x
	normal_distribution<double> dist_x(x, std[0]);
	// (Gaussian) distribution for y
	normal_distribution<double> dist_y(y, std[1]);
	// (Gaussian) distribution for theta
	normal_distribution<double> dist_theta(theta, std[2]);

	//particles = (num_particles);
	Particle *t_particle;
	//cout << "size " << particles.size() << endl;
	for (int i = 0; i < num_particles; ++i) {
		t_particle = new Particle;
		//Sample  x, y and theta from these normal distrubtions
		t_particle->id = i;
		t_particle->x = dist_x(gen);
		t_particle->y = dist_y(gen);
		t_particle->theta = dist_theta(gen);
		t_particle->weight = 1.0;

		particles.push_back(*t_particle);
	}
	is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	// Create a default random engine
	default_random_engine gen;
	double x_f, y_f, theta_f;
	for (int i = 0; i < num_particles; ++i) {
		if (fabs(yaw_rate) > 0.00001) {
			//cout << "yaw " << yaw_rate << endl;
			x_f = particles[i].x + (velocity/yaw_rate) * (sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
			y_f = particles[i].y + (velocity/yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
			theta_f = particles[i].theta + yaw_rate*delta_t;
		} else {
			x_f = particles[i].x + velocity * delta_t * cos(particles[i].theta);
			y_f = particles[i].y + velocity * delta_t * sin(particles[i].theta);
			theta_f = particles[i].theta;
		}

		// (Gaussian) distribution for x
		normal_distribution<double> dist_x(0, std_pos[0]);
		// (Gaussian) distribution for y
		normal_distribution<double> dist_y(0, std_pos[1]);
		// (Gaussian) distribution for theta
		normal_distribution<double> dist_theta(0, std_pos[2]);

		particles[i].x = x_f + dist_x(gen);
		particles[i].y = y_f + dist_y(gen);
		particles[i].theta = theta_f + dist_theta(gen);

		//angle normalization
    //while (particles[i].theta > M_PI) particles[i].theta -= 2.*M_PI;
    //while (particles[i].theta < -M_PI) particles[i].theta += 2.*M_PI;

	}
}

void ParticleFilter::dataAssociation(Particle& particle, std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.
	//cout << "obs size " << observations.size() << endl;
	//cout << "pred size " << predicted.size() << endl;

	std::vector<int> associations;
	std::vector<double> sense_x;
	std::vector<double> sense_y;
	if (predicted.size() > 0) {
		for (int n=0; n<observations.size(); n++) {
			double dist_min = 100000;
			double test_dist;
			int id = -1;
			for (int m=0; m<predicted.size(); m++) {
				test_dist = dist(observations[n].x, observations[n].y, predicted[m].x, predicted[m].y);
				if (test_dist < dist_min) {
					id = m;
					dist_min = test_dist;
				}
			}
			associations.push_back(predicted[id].id);
			sense_x.push_back(observations[n].x);
			sense_y.push_back(observations[n].y);
		}
	}

	particle = SetAssociations(particle, associations, sense_x, sense_y);

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	weights.clear();
	for (int i=0; i<num_particles; i++) {
		std::vector<LandmarkObs> trans_observations;
		// Transform observations from car coord sys to map coord sys at possition of particle
		for (int j=0; j<observations.size(); j++) {
			LandmarkObs obs;
			obs.x = particles[i].x + cos(particles[i].theta)*observations[j].x - sin(particles[i].theta)*observations[j].y;
			obs.y = particles[i].y + sin(particles[i].theta)*observations[j].x + cos(particles[i].theta)*observations[j].y;
			trans_observations.push_back(obs);
		}

		// Get a lst of ladmarks from the map that are within sensor_rage of particle
		std::vector<LandmarkObs> predicted;
		//int count = 0;
		for (int j=0; j<map_landmarks.landmark_list.size(); j++) {
			if (dist(particles[i].x, particles[i].y, map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f) <= sensor_range) {
				LandmarkObs obs;
				obs.id = map_landmarks.landmark_list[j].id_i;
				obs.x = map_landmarks.landmark_list[j].x_f;
				obs.y = map_landmarks.landmark_list[j].y_f;
				predicted.push_back(obs);
				//count++;
			}
		}
		//cout << "node " <<  i << " has " << count << " landmarks\n";

		// Call dataAssociation() to determin the best fit landmark for each sensor observation
		//if (count > 0)
		/*  DIBUG CODE
		if (predicted.size() == 0) {
			cout << "obs\n";
			for (int o=0; o<observations.size(); o++) {
				cout << "ob " << o << " x = " << observations[o].x << " y = " << observations[o].y << endl;
			}
		}
		if (predicted.size() == 0) {
			cout << "trans_obs\n";
			for (int o=0; o<trans_observations.size(); o++) {
				cout << "trans_ob " << o << " x = " << trans_observations[o].x << " y = " << trans_observations[o].y << endl;
			}
		}
		if (predicted.size() == 0) {
			cout << "theta = " << particles[i].theta << endl;
			cout << "obs.x = " << particles[i].x + cos(particles[i].theta)*observations[0].x - sin(particles[i].theta)*observations[0].y << endl;
			cout << "particle x = " << particles[i].x << " y = " << particles[i].y << endl;
		}
		*/
		dataAssociation(particles[i], predicted, trans_observations);

		// Calculate weights using multivariate Gaussian
		double p;
		if (particles[i].associations.size() > 0) {
			p = 1.0;
			for (int j=0; j<particles[i].associations.size(); j++) {
				int id = particles[i].associations[j] - 1;
				double x = particles[i].sense_x[j];
				double y = particles[i].sense_y[j];
				double u_x = map_landmarks.landmark_list[id].x_f;
				double u_y = map_landmarks.landmark_list[id].y_f;
				double p_exp = exp(-(((u_x-x)*(u_x-x)/(2*std_landmark[0]*std_landmark[0]))+((u_y-y)*(u_y-y)/(2*std_landmark[1]*std_landmark[1]))));
				double p_f = 1/(2*M_PI*std_landmark[0]*std_landmark[1]);
				p = p * p_f * p_exp;
			}
		} else {
				p = 1.0;
		}
		particles[i].weight = p;
		weights.push_back(p);

	} // End of particles for loop

	// Normalize weights
	double sum_of_elems = std::accumulate(weights.begin(), weights.end(), 0.0);
	for (int i=0; i<weights.size(); i++) {
		weights[i] = weights[i]/sum_of_elems;
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	// Create a default random engine
	default_random_engine gen;
	// Set descrete distribution with weights of particles
	discrete_distribution<> dist_weights(weights.begin(), weights.end());

	// Set of new list of particles
	std::vector<Particle> new_particles;
	Particle new_particle;

	// Sample new particles
	int new_id = -1;
	for (int i=0; i<num_particles; i++) {
		new_id = dist_weights(gen);
		new_particle.x = particles[new_id].x;
		new_particle.y = particles[new_id].y;
		new_particle.theta = particles[new_id].theta;
		new_particle.weight = particles[new_id].weight;
		new_particle.associations = particles[new_id].associations;
		new_particle.sense_x = particles[new_id].sense_x;
		new_particle.sense_y = particles[new_id].sense_y;
		new_particles.push_back(new_particle);
		//cout << "new particle id " << new_id << endl;
	}
	particles = new_particles;

}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations = associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
