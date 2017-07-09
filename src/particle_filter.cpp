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
//  std::cout << "Init: " << x << " " << y << " " << theta << " " << std[0] << " " << std[1] << " " << std[2] << std::endl;


  if(num_particles == 0) num_particles = 20;

  std::default_random_engine generator;
  std::normal_distribution<double> dist_x(x, std[0]);
  std::normal_distribution<double> dist_y(y, std[1]);
  std::normal_distribution<double> dist_theta(theta, std[2]);


  for(int i = 0; i < num_particles; i++) {
      Particle p;
      p.x = dist_x(generator);
      p.y = dist_y(generator);
      p.theta = dist_theta(generator);
      p.weight = 1;
      particles.push_back(p);
  }

  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
  // construct a trivial random generator engine from a time-based seed:
//  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator;
  std::normal_distribution<double> dist_x(0, std_pos[0]);
  std::normal_distribution<double> dist_y(0, std_pos[1]);
  std::normal_distribution<double> dist_theta(0, std_pos[2]);

  for (auto &p : particles){

      if(fabs(yaw_rate) > 0.001) {
	p.x = p.x + velocity/yaw_rate * (sin(p.theta + yaw_rate*delta_t) - sin(p.theta)) + dist_x(generator);
	p.y = p.y + velocity/yaw_rate * (cos(p.theta) - cos(p.theta + yaw_rate*delta_t)) + dist_y(generator);
      } else {
	p.x = p.x + velocity*cos(p.theta)*delta_t + dist_x(generator);
	p.y = p.y + velocity*sin(p.theta)*delta_t + dist_y(generator);
      }
      p.theta = p.theta + yaw_rate*delta_t + dist_theta(generator);
  }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

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
//  for(int i = 0; i < observations.size(); i++) {
//      LandmarkObs obs = observations[i];
//      std::cout << "observation x:" << obs.x << " y:" << obs.y << std::endl;
//  }
  for(int i = 0; i < num_particles; i++) {
      Particle &p = particles[i];
      std::vector<LandmarkObs> obs_map;

      p.weight = 1.0;

      //loop through and convert each observation into map coordinates
//      for(int j = 0; j < observations.size(); j++){
//	  LandmarkObs obs = observations[j];
      for(auto &obs : observations) {
//	  std::cout << "p:" << p.x << "," << p.y << std::endl;
//	  std::cout << "obs prior:" << obs.x << "," << obs.y << std::endl;

	  double newX = p.x + obs.x*cos(p.theta) - obs.y*sin(p.theta);
	  double newY = p.y + obs.x*sin(p.theta) + obs.y*cos(p.theta);

//	  std::cout << "obs:" << obs.x << "," << obs.y << std::endl;
//	  std::cout << "obs original:" << observations[j].x << "," << observations[j].y << std::endl;

	  double shortest_dist = -1;
	  Map::single_landmark_s *closest;
//	  for(int k = 0; k < map_landmarks.landmark_list.size(); k++){
//	      Map::single_landmark_s &landmark = map_landmarks.landmark_list[k];
	  for(auto &landmark : map_landmarks.landmark_list) {
	      double distance = dist(newX, newY, landmark.x_f, landmark.y_f);

	      if(distance > sensor_range) continue;
	      if(shortest_dist == -1 || distance < shortest_dist){
		  closest = &landmark;
		  shortest_dist = distance;
//		  std::cout << "assigned closest" << std::endl;
	      }

//	      std::cout << "landmark:" << landmark.x_f << ", " << landmark.y_f << " distance:" << distance << std::endl;
	  }

//	  std::cout << "closest landmark:" << closest->x_f << ", " << closest->y_f << " distance:" << shortest_dist << std::endl;
	  //what happens when std_landmark[0] is 0?
	  double left = 1 / (2*M_PI*std_landmark[0]*std_landmark[1]);
	  double right1 = pow((newX - closest->x_f), 2)/(2.0*pow(std_landmark[0], 2));
	  double right2 = pow((newY - closest->y_f), 2)/(2.0*pow(std_landmark[1], 2));
	  double multivariate_gaussian = left * exp(-(right1 + right2));
//	  std::cout << "test: " << (2*pow(std_landmark[0], 2)) << std::endl;
//	  std::cout << "left: " << left << std::endl;
//	  std::cout << "right1: " << right1 << std::endl;
//	  std::cout << "right2: " << right2 << std::endl;
//	  std::cout << "multivariate_gaussian: " << multivariate_gaussian << std::endl;
	  p.weight *= multivariate_gaussian;
      }

  }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
  std::vector<double> weights;
  for(int i = 0; i < particles.size(); i++) {
      weights.push_back(particles[i].weight);
  }

  std::random_device rd;
  std::mt19937 gen(rd());
  std::discrete_distribution<> dist(std::begin(weights), std::end(weights));

  std::vector<Particle> new_particles;
  for(int i = 0; i < num_particles; i++) {
      new_particles.push_back(particles[dist(gen)]);
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

	particle.associations= associations;
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
