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

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  
  // Already initialized 
  if (is_initialized)
    return;
  
  // Set number of particles
  num_particles = 1000;  // TODO: Set the number of particles
  
  std::default_random_engine gen;
  double std_x, std_y, std_theta;
  
  // Set std params for x, y and theta
  std_x=std[0];
  std_y=std[1];
  std_theta=std[2];
  
  // This line creates a normal (Gaussian) distribution for x, y and theta
  normal_distribution<double> dist_x(x, std_x);
  normal_distribution<double> dist_y(y, std_y);
  normal_distribution<double> dist_theta(theta, std_theta);

  // Iterate over each particle to assign x,y, theta randomly
  for (int i = 0; i < num_particles; i++)
  {
  	Particle curr_p;
    curr_p.id = i;
    curr_p.x = dist_x(gen);
    curr_p.y = dist_y(gen);;
    curr_p.theta = dist_theta(gen);
    curr_p.weight = 1;
    particles.push_back(curr_p);
  }
  // Set is_initialized flag to 1
  is_initialized = 1;

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
   std::default_random_engine gen;
  
  // This line creates a normal (Gaussian) distribution for sensor noise
   normal_distribution<double> dist_x(0, std_pos[0]);
   normal_distribution<double> dist_y(0, std_pos[1]);
   normal_distribution<double> dist_theta(0, std_pos[2]);
  
   for (int i = 0; i < num_particles; i++)
   {
     
     double particle_theta = particles[i].theta;
     
     // Define new states
     if (fabs(yaw_rate) < 0.0001)
     {
     	particles[i].x += velocity*delta_t*cos(particle_theta);
        particles[i].y += velocity*delta_t*sin(particle_theta);
     } else {
        particles[i].x += velocity / yaw_rate * (sin(particle_theta + yaw_rate * delta_t) - sin(particle_theta));
        particles[i].y += velocity/ yaw_rate * (cos(particle_theta)- cos(particle_theta + yaw_rate * delta_t));
        particles[i].theta += yaw_rate *delta_t;
     }
     // Add noise
     particles[i].x += dist_x(gen);
     particles[i].y += dist_y(gen);
     particles[i].theta += dist_theta(gen);
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
   // Get observations and predictions size
   unsigned int landmark_observations = observations.size();
   unsigned int landmark_predictions = predicted.size();
  
  // Iterate over each observations
   for (unsigned int i = 0; i < landmark_observations; i++)
   {
     //Initialize to very big number
     double min_distance = numeric_limits<double>::max();
     
     //Initialize map id
     int map_id = -1;
     
     // iterate over each predictions
     for (unsigned int j = 0; j < landmark_predictions; j++)
     {
       // Calculate distance
       double x_dist = observations[i].x - predicted[j].x;
       double y_dist = observations[i].y - predicted[j].y;
       double distance = x_dist*x_dist + y_dist *y_dist;
       
       // Set observation close to predicted using nearest neighbour 
       // Update min_distance
       if (distance < min_distance)
       {
         min_distance = distance;
         map_id = predicted[j].id;
       }
     }
     
     // set the observation id to the nearest predicted landmark id
     observations[i].id = map_id;
      
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
  
  
   double std_landmark1 = std_landmark[0]; 
   double std_landmark2 = std_landmark[1]; 
  
  //Iterate over each particles
  for (int i = 0; i< num_particles; i++)
  {
    //Get x,y and theta for each particle
    double x = particles[i].x;
    double y = particles[i].y;
    double theta = particles[i].theta;
    
    double sensor_range2 = sensor_range * sensor_range;
    
    //Create vector to hold transformed data
    vector<LandmarkObs> transformed;
    
    //Iterate over each landmark
    //Find landmark in particle range
    for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++)
    {
      // Set x,y and id of landmark
      float landmarkx = map_landmarks.landmark_list[j].x_f;
      float landmarky = map_landmarks.landmark_list[j].y_f;
      int id = map_landmarks.landmark_list[j].id_i;
      
      // Calculate distance
      double distance_x = x - landmarkx;
      double distance_y = y - landmarky;
      double distance = distance_x * distance_x + distance_y * distance_y ;
      
      // Update distance and transformed vector
      if (distance <= sensor_range2)
      {
        transformed.push_back(LandmarkObs{id, landmarkx, landmarky});
      }
      
    }
    //Transform observation cordinates
    vector<LandmarkObs> mapped_observations;
    for (unsigned int j = 0; j < observations.size(); j++)
    {
      double dx = x + cos(theta)*observations[j].x - sin(theta)*observations[j].y;
      double dy = y + sin(theta)*observations[j].x + cos(theta)*observations[j].y;
      mapped_observations.push_back(LandmarkObs{observations[j].id, dx, dy});
    }
    // Association of transformed cordinates with landmarks
    dataAssociation(transformed, mapped_observations);
    
    //Reset weights
    particles[i].weight = 1.0;
    
    //Calculate weights
    for (unsigned int j = 0; j <transformed.size(); j++)
    {
      double tx = transformed[j].x;
      double ty = transformed[j].y;
      int ti = transformed[j].id;
      double px,py;
      for (unsigned int k = 0; k < mapped_observations.size(); k++)
      {
        if (mapped_observations[k].id == ti)
        {
          px = mapped_observations[k].x;
          py = mapped_observations[k].y;
        }
      }
      //Update weights for particles
      double obs_w = (1/(2*M_PI*std_landmark1*std_landmark2))* exp(-(pow(px-tx,2)/(2*pow(std_landmark1,2))+(pow(py-ty,2)/(2*pow(std_landmark2,2))) ) );
      if (obs_w == 0)
      {
      	particles[i].weight *= 0.00001;
      }
      else
      {
      	particles[i].weight *= obs_w;
      }
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
   vector<double> w;
  default_random_engine gen;
   //Get weight and maximum weight
   double max_w = numeric_limits<double>::min();
   for (int i = 0; i < num_particles; i++)
   {
     w.push_back(particles[i].weight);
     if (particles[i].weight > max_w)
     {
       max_w = particles[i].weight;
     }
   }
  // Create uniform distributions for weights and particles
   uniform_real_distribution<double> distDouble(0.0, max_w);
   uniform_int_distribution<int> distInt(0, num_particles-1);
  
  //Update beta and index using Wheel approach
  int index = distInt(gen);
  double beta = 0.0;
  
  vector<Particle> resampled;
  for (int i = 0; i< num_particles; i++)
  {
    beta += distDouble(gen)*2.0;
    while (beta > w[index]) 
    {
      beta -= w[index];
      index = (index +1)%num_particles;
    }
    resampled.push_back(particles[index]);
    
  }
  particles = resampled;
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