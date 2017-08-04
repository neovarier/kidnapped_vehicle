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
	default_random_engine gen;
    Particle particle;

    num_particles = 1000;
	// This line creates a normal (Gaussian) distribution for x.
	normal_distribution<double> dist_x(x, std[0]);
	// This line creates a normal (Gaussian) distribution for y.
	normal_distribution<double> dist_y(y, std[1]);
	// This line creates a normal (Gaussian) distribution for psi.
	normal_distribution<double> dist_theta(theta, std[2]);

    
    weights.resize(num_particles,1.0);

    for(int i = 0;i < num_particles; i++)
    {
      particle.id = i;
	  particle.x = dist_x(gen);
	  particle.y = dist_y(gen);
	  particle.theta = dist_theta(gen);
	  particle.weight = 1.0;
      particles.push_back(particle);
    }
    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
    double x_update;
    double y_update;
    double theta_update;
    default_random_engine gen;

    for (int i = 0;i < num_particles; i++)
    {
      if (yaw_rate == 0)
      {
        x_update = particles[i].x + velocity*delta_t*cos(particles[i].theta);
        y_update = particles[i].y + velocity*delta_t*sin(particles[i].theta);

        // This line creates a normal (Gaussian) distribution for x.
        normal_distribution<double> dist_ctrl_x(x_update, std_pos[0]);
        // This line creates a normal (Gaussian) distribution for y.
        normal_distribution<double> dist_ctrl_y(y_update, std_pos[1]);
        particles[i].x = dist_ctrl_x(gen);
        particles[i].y = dist_ctrl_y(gen);
      }
      else
      {
        theta_update = particles[i].theta + yaw_rate*delta_t;
        x_update = particles[i].x + (velocity/yaw_rate)*(sin(theta_update) - sin(particles[i].theta));
        y_update = particles[i].y + (velocity/yaw_rate)*(cos(particles[i].theta) - cos(theta_update));

        // This line creates a normal (Gaussian) distribution for control x.
        normal_distribution<double> dist_ctrl_x(x_update, std_pos[0]);
        // This line creates a normal (Gaussian) distribution for control y.
        normal_distribution<double> dist_ctrl_y(y_update, std_pos[1]);
        // This line creates a normal (Gaussian) distribution for control theta.
        normal_distribution<double> dist_ctrl_theta(theta_update, std_pos[2]);

        particles[i].x = dist_ctrl_x(gen);
        particles[i].y = dist_ctrl_y(gen);
        particles[i].theta = dist_ctrl_theta(gen);
      }
    }

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> close_landmarks, std::vector<LandmarkObs> observations, int particle_index) {
	// TODO: Find the observed  measurement that is closest to each close landmark and associate the 
	//   the landmark with observed measurement.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

    /*Find the closest observation  w.r.t a landmark*/
    /*Remove the observations once it is associated*/
    /*There might be cases where number of close landmarks are greater than the number of observations*/
    /*If all the observations are associated then stop the loop*/
    
    std::vector<double> c_sense_x;
    std::vector<double> c_sense_y;
    std::vector<int> p_associations;
    double distance;
    double min_dist;
    int closest_obs;

    for (int k =0; ((k < close_landmarks.size() && (observations.size() > 0))); k++)
    {
      for (int j =0; j < observations.size(); j++)
      {
        distance = dist(close_landmarks[k].x, close_landmarks[k].y, observations[j].x, observations[j].y);
        if (j == 0)
        {
          min_dist = distance;
          closest_obs = j;
        }
        else
        {
          if (min_dist > distance)
          {
            min_dist=distance;
            closest_obs=j;
          }
        }
      }
      /*Set the id of the associated landmark for the correxponding observation*/
      p_associations.push_back(close_landmarks[k].id);
      c_sense_x.push_back(observations[closest_obs].x);
      c_sense_y.push_back(observations[closest_obs].y);

      /*Remove the landmark that is already associated from the list*/
      observations.erase(observations.begin()+closest_obs);
    }
    particles[particle_index] = SetAssociations(particles[particle_index], p_associations, c_sense_x, c_sense_y);
}

int ParticleFilter::getLMIndex(int landmark_id, Map map_landmarks)
{
  int lm_index;
  /*Find the index of the landmark from the landmark list based on the landmark id*/
  for(int i=0; i<map_landmarks.landmark_list.size(); i++)
  {
    if(map_landmarks.landmark_list[i].id_i == landmark_id)
    {
      lm_index = i;
      break;
    }
  }
  return lm_index;
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
    double p_sense_x;
    double p_sense_y;
    double distance;
    std::vector<LandmarkObs> close_landmarks;
    std::vector<LandmarkObs> trans_observations;

    for(int i =0; i< num_particles; i++)
    {
      for (int k=0; k< map_landmarks.landmark_list.size(); k++)
      {
        distance = dist((double)map_landmarks.landmark_list[k].x_f, (double)map_landmarks.landmark_list[k].y_f, particles[i].x, particles[i].y);
        if (distance < sensor_range)
        {
          LandmarkObs landmark;
          landmark.x = (double)map_landmarks.landmark_list[k].x_f;
          landmark.y = (double)map_landmarks.landmark_list[k].y_f;
		  landmark.id = map_landmarks.landmark_list[k].id_i;
          close_landmarks.push_back(landmark);
        }
      }
      LandmarkObs obs;
      obs.x = 0.0;
      obs.y = 0.0;
      obs.id = 0;
      trans_observations.resize(observations.size(), obs);
      for (int j =0; j<observations.size(); j++)
      {
        /*Transform the observation from car coordinates to map coordinates*/
        /*Rotation based on theta of the particle*/
        /*Translation based on the (x,y) position of the particle*/
        /*Do this for each observation*/

        trans_observations[j].x  =  observations[j].x*cos(particles[i].theta) - observations[j].y*sin(particles[i].theta) + particles[i].x;
        trans_observations[j].y  =  observations[j].x*sin(particles[i].theta) + observations[j].y*cos(particles[i].theta) + particles[i].y;
        trans_observations[j].id =  observations[j].id;
      }
       
      /*Find the closest map landmark w.r.t each observation and set in the particle structure*/
      dataAssociation(close_landmarks, trans_observations, i);
      close_landmarks.clear();
      trans_observations.clear();
    }

    /*Calcualte weights for each particle*/
    for(int i =0; i<num_particles; i++)
    {
      particles[i].weight = 1.0;
      for(int j =0; j<particles[i].associations.size(); j++)
      {
        /*Get index of the associated landmark based on its id*/
        int lm_index = getLMIndex(particles[i].associations[j], map_landmarks);
        /*Calculated the weight based on multi-variate normal distribution*/
        /*weight will be product of probabilities of each observation w.r.t the associated landmark */
        double exp_num_x = squared(particles[i].sense_x[j] - map_landmarks.landmark_list[lm_index].x_f)/(2*squared(std_landmark[0]));
        double exp_num_y = squared(particles[i].sense_y[j] - map_landmarks.landmark_list[lm_index].y_f)/(2*squared(std_landmark[1]));
        particles[i].weight *=  exp(-(exp_num_x+exp_num_y))/(2*M_PI*std_landmark[0]*std_landmark[1]);
      }
    }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    std::vector<double> int_weights;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::vector<Particle> rsmpl_particles;
    double sum = 0.0;
    /*Calcualte the sum of the weights for normalization*/
    for(int i=0;i< num_particles; i++)
    {
      sum += particles[i].weight;
    }
    /*Calcualte the integer weights based on total number of particles*/
    for(int i=0; i< num_particles; i++)
    {
      int_weights.push_back(int(num_particles*particles[i].weight/sum));
    }
    /*Create a discrete distrbution based on the integer weights*/
    discrete_distribution<> weight_dist(int_weights.begin(), int_weights.end());
    /*pick particles based on the distribution*/
    for(int i=0; i< num_particles; i++)
    {
      int index = weight_dist(gen);
      rsmpl_particles.push_back(particles[index]);
      
    }
    particles.clear();
    particles = rsmpl_particles;
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
