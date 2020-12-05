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

void ParticleFilter::init(double x, double y, double theta, double std[])
{
  /**
   * TODO: Set the number of particles. Initialize all particles to
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1.
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method
   *   (and others in this file).
   */
  // TODO: Set the number of particles. I will begin with 25. and Add 25 until I see the best result
  num_particles = 25;

  // According to lesson 5: Program Gaussian Sampling: Code
  // First, I will create the normal (Gaussian) distribution for X, Y and theta
  using std::normal_distribution;
  std::default_random_engine gen;

  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

  // in the particle_filter.h is a Particle structure with variables
  // id, x, y, theta and weight.
  // Set of current particles is a vector called particles
  // Set of weights is a vector called weights. I need to change the vector size
  // to the same size of the number of particles and initialize them with 1.0
  //https://www.cplusplus.com/reference/vector/vector/resize/
  weights.resize(num_particles, 1.0);

  // Loop over the particles
  for(int p = 0; p < num_particles; ++p)
  {
    Particle particle; // initialize the particle structure
    particle.id = p;
    particle.x = dist_x(gen);
    particle.y = dist_y(gen);
    particle.theta = dist_theta(gen);
    particle.weight = 1.0;
    particles.push_back(particle);

  }is_initialized = true; // Let the program know that the initialization is complete
}

void ParticleFilter::prediction(double delta_t, double std_pos[],
                                double velocity, double yaw_rate)
{
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  using std::normal_distribution;
  std::default_random_engine gen;

  // First, loop over all the particles, to predict one by one the new position
  for(int cont = 0; cont < num_particles;++cont)
  {
    Particle particle = particles[cont];
    // In the prediction equations, there are two cases
    // First, when yaw rate is zero, and the second when yaw rate is not zero
    if(yaw_rate < 0.0001)
    {
      // I will use the given equations for yaw rate equal to 0
      //X_f = X_0 + V(dt)*cos(theta_0);
      //Y_f = Y_0 + V(dt)*sin(theta_0);
      //theta_f = theta_0;
      particle.x = particle.x + velocity*delta_t*cos(particle.theta);
      particle.y = particle.y + velocity* delta_t* sin(particle.theta);
      particle.theta = particle.theta; // Maybe unnecessary, just for the sake of representing the OG equations
    }
    else
    {
      // I will use the given equations for yaw rate different than zero
      // X_f = X_0 + (V/theta_dot)*(sin(theta_0 + theta_dot*delta_t)- sin(theta_0))
      // Y_f = Y_0 + (V/theta_dot)*(cos(theta_0)- cos(theta_0) +theta_dot*delta_t)
      // Theta_f = theta_0 + theta_dot*delta_t
      particle.x = particle.x + (velocity / yaw_rate) * (sin(particle.theta + yaw_rate * delta_t) - sin(particle.theta));
      particle.y = particle.y + (velocity / yaw_rate) * (cos(particle.theta) - cos(particle.theta + yaw_rate * delta_t));
      particle.theta = particle.theta + yaw_rate * delta_t;
    }

    // Like in the question of a previous student, I will add noise to each predicted value
    // https://knowledge.udacity.com/questions/148019

    normal_distribution<double> dist_x(particle.x, std_pos[1]);
    normal_distribution<double> dist_y(particle.y, std_pos[2]);
    normal_distribution<double> dist_theta(particle.theta, std_pos[3]);

    // Update the particle with the gaussian noise
    particle.x = dist_x(gen);
    particle.y = dist_y(gen);
    particle.theta = dist_theta(gen);
    particles[cont] = particle;
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted,
                                     vector<LandmarkObs>& observations)
{
  /**
   * TODO: Find the predicted measurement that is closest to each
   *   observed measurement and assign the observed measurement to this
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will
   *   probably find it useful to implement this method and use it as a helper
   *   during the updateWeights phase.
   */

  // At first, I didn't know what were each landmark vector, but this question help me a lot
  // https://knowledge.udacity.com/questions/360973

  // First, loop over each observation
  for (int i = 0; i < observations.size(); i++)
  {
    // As a Landmark object[id, x, y], I assign the observation to it
    LandmarkObs obs_landmark = observations[i];

    // Since I want a minimum distance, I will set a standard distance as minimum
    // The max range of a landmark should be 50 (sensor range) so I use a value above that value
    double min_dist = 1000.;

    // I need to "Remember" which particle was the closest to the landmark.
    int closest_id = -1;

    // Now, I will loop over each predicted particle and check which one is closer to the landmark
    for (int j = 0; j < predicted.size(); j++)
    {
      // As a Landmark object, I assign the predicted particle  to it
      LandmarkObs pred_landmark = predicted[j];

      // I will use a helper_functions.h method to calculate the distance between
      // the actual landmark and the predicted landmark
      // calculate distance between landmark observation and prediciton
      double obs_pred_dist = dist(obs_landmark.x, obs_landmark.y, pred_landmark.x, pred_landmark.y);

      // Here is where I compare the calculated distance, with the previous minimum distance
      // If it is less, I assign it to the minimum distance variable

      if (obs_pred_dist < min_dist)
      {
        min_dist = obs_pred_dist;
        closest_id = pred_landmark.id;
      }
    }
    // Finally, I assign the ID to the observation object ID
    observations[i].id = closest_id;
  }
}

double ParticleFilter::transFunct(double x_part, double y_part, double x_obs, double y_obs, double theta)
{
  // I will create a new method for the observation transformation
  // and will do as in the Landmarks class of Implementation of a Particle Filter
  // transform to map x coordinate
  double x_map;
  x_map = x_part + (cos(theta) * x_obs) - (sin(theta) * y_obs);

  // transform to map y coordinate
  double y_map;
  y_map = y_part + (sin(theta) * x_obs) + (cos(theta) * y_obs);

  return x_map, y_map;
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   const vector<LandmarkObs> &observations,
                                   const Map &map_landmarks)
{
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian
   *   distribution. You can read more about this distribution here:
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system.
   *   Your particles are located according to the MAP'S coordinate system.
   *   You will need to transform between the two systems. Keep in mind that   <----------
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

  // I need the landmark standard deviations
  double std_x = std_landmark[0];
  double std_y = std_landmark[1];

  // For the dataAssociation function, I need TWO LandmarkObs vectors.
  // I was somewhat confused for what I needed for the dataAssociation method.
  // But according to this question, I need ACTUAL landmarks and predicted observations
  // https://knowledge.udacity.com/questions/360973
  //vector<LandmarksObs> close_landmarks;
  //vector<LandmarksObs> predicted_obs;
  // I need a vector with ALL the landmarks in the map for the Association
  // https://knowledge.udacity.com/questions/56469
  vector<Map::single_landmark_s> landmarks = map_landmarks.landmark_list;

  // I need two things to calculate the particles weights.
  // First, I need the landmarks positions that are in the range of the sensor.
  // Secondly, I need the predicted observations from the car IN MAP COORDINATES

  // First, loop over all the particles
  for(int p = 0; p <num_particles;p++)
  {
    //assign each predicted position to a variable
    double x_particle = particles[p].x;
    double y_particle = particles[p].y;
    double theta_particle = particles[p].theta;
    // For each particle, I have to reset to 1 the weight value
    double weight = 1.;

    // my second loop will be for the observations
    for(int o = 0; o < observations.size(); o++)
    {
      double vobs_x = observations[o].x;
      double vobs_y = observations[o].y;
      // I will use the Transformation method to get the sensor observations in MAP COORDINATES
      double mobs_x;
      double mobs_y;
      mobs_x, mobs_y = ParticleFilter::transFunct(x_particle, y_particle, vobs_x, vobs_y, theta_particle);

      //predicted_obs.push_back(LandmarkObs{observations[o].id, mobs_x, mobs_y });
      // As recomended, I will use a nested loop for the landmarks inside the Observations loop
      // If I am not wrong, I can check the actual observation to all the Landmarks
      // This way I can check which observation is closer to which landmarks (data association)
      float min_dist = 1000.;
      double closest_lm_x;
      double closest_lm_y;

      for(int l = 0; l < landmarks.size(); l++)
      {
        // get x and y coordinates
        double actual_lm_x = landmarks[l].x_f;
        double actual_lm_y = landmarks[l].y_f;
        //int actual_lm_id = landmarks[l].id_i;

        double landmark_dist = dist(actual_lm_x, actual_lm_y, mobs_x, mobs_y);
        // Now, I just need to check which landmark within the sensor range is closest to the actual observation

        if(landmark_dist < min_dist && landmark_dist < sensor_range)
        {
          min_dist = landmark_dist;
          closest_lm_x = actual_lm_x;
          closest_lm_y = actual_lm_y;

          // I will add the landmarks that are close to the particles
          //close_landmarks.push_back(LandmarkObs{ actual_lm_id, actual_lm_x, actual_lm_y });
        }
      }
      //For the data association, I will take all the close Landmarks
      //dataAssociation(predicted_obs[o], close_landmarks);

      // Now I have the closest landmark to my actual observation.
      // I will calculate the weights for the particle

      // First, I calculate normalization term
      double gauss_norm;
      gauss_norm = 1 / (2 * M_PI * std_x * std_y);
      // calculate exponent
      double exponent;
      exponent = (pow(mobs_x - closest_lm_x, 2) / (2 * pow(std_x, 2))) +
        		(pow(mobs_y - closest_lm_y, 2) / (2 * pow(std_y, 2)));
      // calculate weight using normalization terms and exponent
      weight *= (gauss_norm * exp(-exponent));
    }
    particles[p].weight = weight;

  }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	std::random_device rd;
    std::mt19937 gen(rd());
    std::discrete_distribution<> d(weights.begin(), weights.end());
    std::vector<Particle> resampled_particles;

    for(int n = 0; n < num_particles; ++n)
    {
        Particle particle = particles[d(gen)];
        resampled_particles.push_back(particle);
    }
    particles = resampled_particles;
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
