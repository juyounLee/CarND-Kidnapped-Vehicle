  
/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using namespace std;

std::default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  if(is_initialized) {
		return;
	} 
  num_particles = 100;  // TODO: Set the number of particles
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);
  
  for (int i=0;i<num_particles;i++){
    Particle particle;
    particle.id = i;
    particle.x = dist_x(gen);
    particle.y = dist_y(gen);
    particle.theta = dist_theta(gen);
    particle.weight = 1.0;
    
    particles.push_back(particle);
    weights.push_back(particle.weight);
    
  }
  
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  
  // random number generator
  std::default_random_engine gen;
  
  //use 0 as mean values
  normal_distribution<double> dist_x(0, std_pos[0]);
  normal_distribution<double> dist_y(0, std_pos[1]);
  normal_distribution<double> dist_theta(0, std_pos[2]);
  
  
  // Predict the vehicle state
  for(auto& particle : particles){
  	// If yaw rate is not 0:
  	if(yaw_rate!=0.0){
  	particle.x += (velocity/yaw_rate)*(sin(particle.theta+(yaw_rate * delta_t)) - sin(particle.theta));
  	particle.y += (velocity/yaw_rate)*(cos(particle.theta)-cos(particle.theta+yaw_rate * delta_t));
    particle.theta += yaw_rate * delta_t;
  }
  //If yaw rate is 0
  else {
    particle.x += velocity * delta_t * cos(particle.theta);
    particle.y += velocity * delta_t * sin(particle.theta);
    particle.theta += 0;  
  }
  //add random gaussian noise
   particle.x += dist_x(gen);
   particle.y += dist_y(gen);
   particle.theta += dist_theta(gen);
 }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  for (unsigned int i=0;i<observations.size();i++){
    	//initialize the minimum distance
    	double min_dist = numeric_limits<double>::max();
    	
    	int closest_landmark;
    	//while looping through the obserations vector, loop through the predicted vector
    	for (unsigned int j=0;j<predicted.size();j++){

    		double distance = dist(observations[i].x,observations[i].y,predicted[j].x,predicted[j].y);
    		if(distance<min_dist){
      			min_dist = distance;
      			closest_landmark = predicted[j].id;
    		}
    	observations[i].id=closest_landmark;
    	}
  	}
}
void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  
  // Homogenous Transformation

  double weight_sum = 0;
  
  for(int i = 0; i < num_particles; i++ ){
    vector<LandmarkObs> predictions;     
    
    for(unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++) {
           
      if( dist(map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f,
               particles[i].x, particles[i].y) <= sensor_range){
        predictions.push_back(LandmarkObs{ map_landmarks.landmark_list[j].id_i, 
                                           map_landmarks.landmark_list[j].x_f,
                                           map_landmarks.landmark_list[j].y_f });
      }
    }
    
    vector<LandmarkObs> transform_obs;
    double transform_x;
    double transform_y;
    
    for(unsigned int j = 0; j < observations.size(); j++){
      transform_x = (cos(particles[i].theta)*observations[j].x) - (sin(particles[i].theta)*observations[j].y) + particles[i].x;
      transform_y = (sin(particles[i].theta)*observations[j].x) + (cos(particles[i].theta)*observations[j].y) + particles[i].y;
      transform_obs.push_back(LandmarkObs{ -1, transform_x, transform_y });
    }
    
    dataAssociation(predictions, transform_obs);
    
    particles[i].weight = 1.0;
    
    double pred_x, pred_y;
    double power;
    double weight;
    
    for(unsigned int j = 0; j < transform_obs.size(); j++){
      for(unsigned int k = 0; k < predictions.size(); k++){
        if(predictions[k].id == transform_obs[j].id){
          pred_x = predictions[k].x;
          pred_y = predictions[k].y;
          break;
        }
    }
      
      power = (pow(transform_obs[j].x - pred_x, 2) / (2 * pow(std_landmark[0], 2)))
               + (pow(transform_obs[j].y - pred_y, 2) / (2 * pow(std_landmark[1], 2)));
      weight = 1 / (2 * M_PI * std_landmark[0] * std_landmark[1]) * exp(-power);
      particles[i].weight *= weight;
  }
      weights[i] = particles[i].weight;
      weight_sum += weights[i];
  }
  if(fabs(weight_sum) > 0.0){
    for(unsigned int i=0;i<weights.size();i++){
      weights[i] = weights[i] / weight_sum;
    }
  }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  std::default_random_engine gen;
  
  vector<Particle> resample;
   discrete_distribution<size_t> randomIndex(weights.begin(), weights.end());
  
    for (unsigned int i = 0; i < particles.size(); i++) {
        resample.push_back(particles[randomIndex(gen)]);
    }  
    particles = resample;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}