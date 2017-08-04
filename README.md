# Kidnapped Vehicle Localization

Steps followed for complete the Particle Filter to localize the Kidnapped Vehicle

## Initialization
* Initialized the number of particles to 1000.
* Initialized the x,y,theta parameter with initial GPS coordinates considering the sensor noise for each particle.
* Initialized the weights of each particle to 1.0.

## Prediction
* Predicted each particle's next location (x,y,theta) based on the measured yaw rate and velocity
* Added measurement noise
* Considered cases where yaw rate is zero and non-zero

## Update Weights
Following steps are done for each particle
* Found the closest landmarks within the sensor range for a particle
* Transform each observation from car coordinate system to map coordinate system
* Find the observation closest to each landmark by nearest neighbour method
* It could be possible that the number of landmarks are more than observations.
In such a case, do the above step only for the available observations.
* Maintain the association of landmarks with the closest observation.

Find the weight of each particle based on the multi-variate Gaussian distribution
The multi-variate Gaussian distribution is the product of a multi-variate Gaussian distribution for each measurement
For each measurment, the multi-variate gaussian distributions a mean (x,y) coordinates of landmark and std deviation of landmark measurement.

## Resample particles
Based on the weights, particles are resamples.
Particles in the new set have instances as per their weights.
Using discrete_distribution, the particles are resampled

## Results
The simulator gives the following message "Success! Your particle filter passed!"
Some of the recorded error in positions:
* x =.111, y =.101, yaw=0.004
* x =.109, y =.102, yaw=0.004
* x =.109, y =.103, yaw=0.004

As we reduce the number of particles, the total time for the program decreases. But the error in the predicted position increases.
For num_particles=3, the particle filter does not succeed.
For num_particles=4, the particle filter succeeds.

Beyond 200 particles, the filter tends to achieve the above good error range.

The code from lessons were helpful.
Understood the usage of discrete_distribution from http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution.
