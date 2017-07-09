#include "gtest/gtest.h"
#include "particle_filter.h"

double sigma_pos[3] =
  { 0.3, 0.3, 0.01 }; // GPS measurement uncertainty [x [m], y [m], theta [rad]]
double sigma_landmark[2] =
  { 0.3, 0.3 }; // Landmark measurement uncertainty [x [m], y [m]]

void printParticle(Particle p) {
  std::cout << "X:" << p.x << " Y:" << p.y << " Theta:" << p.theta << " Weight:" << p.weight << std::endl;
}

std::vector<LandmarkObs> getObservations() {
  std::vector<LandmarkObs> observations;
  LandmarkObs obs;
  obs.x = 2;
  obs.y = 2;
  observations.push_back(obs);

  obs = LandmarkObs();
  obs.x = 3;
  obs.y = -2;
  observations.push_back(obs);

  obs = LandmarkObs();
  obs.x = 0;
  obs.y = -4;
  observations.push_back(obs);
  return observations;
}

Map getMap() {
  Map map;
  Map::single_landmark_s landmark;

  landmark.id_i = 1;
  landmark.x_f = 5;
  landmark.y_f = 3;
  map.landmark_list.push_back(landmark);

  landmark = Map::single_landmark_s();
  landmark.id_i = 2;
  landmark.x_f = 2;
  landmark.y_f = 1;
  map.landmark_list.push_back(landmark);

  landmark = Map::single_landmark_s();
  landmark.id_i = 3;
  landmark.x_f = 6;
  landmark.y_f = 1;
  map.landmark_list.push_back(landmark);

  landmark = Map::single_landmark_s();
  landmark.id_i = 4;
  landmark.x_f = 7;
  landmark.y_f = 4;
  map.landmark_list.push_back(landmark);

  landmark = Map::single_landmark_s();
  landmark.id_i = 5;
  landmark.x_f = 4;
  landmark.y_f = 7;
  map.landmark_list.push_back(landmark);

  return map;
}



TEST(ParticleFilter, test_init) {
  ParticleFilter pf(2);

  pf.init(6.315, 1.9479, 0.0111, sigma_pos);
  EXPECT_TRUE(pf.particles.size()==2);

//  printParticle(pf.particles[0]);
}

TEST(ParticleFilter, test_prediction_no_noise) {
  ParticleFilter pf(1);

  double zero_noise[3] = {0.0, 0.0, 0.0};
  pf.init(102.0, 65.0, 5*M_PI/8, zero_noise);
  pf.prediction(0.1, zero_noise, 110.0, M_PI/8);

  Particle p = pf.particles[0];
  EXPECT_DOUBLE_EQ(p.x, 97.592046082729098);
  EXPECT_DOUBLE_EQ(p.y, 75.07741997215382);
  EXPECT_DOUBLE_EQ(p.theta, 51*M_PI/80);

//  printParticle(p);
}

TEST(ParticleFilter, test_prediction_no_noise_zero_yawrate) {
  ParticleFilter pf(1);

  double zero_noise[3] = {0.0, 0.0, 0.0};

//  std::cout << "test:" << atan2(4, 3) << " " << atan2(4, 3)/M_PI * 180 << std::endl;

  pf.init(3, 4, atan2(4, 3), zero_noise);
  pf.prediction(1, zero_noise, 5, 0);

  Particle p = pf.particles[0];
  EXPECT_DOUBLE_EQ(p.x, 6);
  EXPECT_DOUBLE_EQ(p.y, 8);
  EXPECT_DOUBLE_EQ(p.theta, atan2(4,3));

//  printParticle(p);
}

TEST(ParticleFilter, test_prediction_with_noise) {
  ParticleFilter pf(1);

  double zero_noise[3] = {0.0, 0.0, 0.0};
  pf.init(102.0, 65.0, 5*M_PI/8, zero_noise);

  pf.prediction(0.1, sigma_pos, 110.0, M_PI/8);

  Particle p = pf.particles[0];
//  EXPECT_DOUBLE_EQ(p.x, 97.592046082729098);
//  EXPECT_DOUBLE_EQ(p.y, 75.07741997215382);
//  EXPECT_DOUBLE_EQ(p.theta, 51*M_PI/80);
//  std::cout << std::chrono::high_resolution_clock::now() << std::endl;
//  printParticle(p);
}


TEST(ParticleFilter, test_update_weight_no_noise) {
  //using the simple scenario in the lectures
  ParticleFilter pf(1);

  double zero[3] = {0.0, 0.0, 0.0};
  pf.init(4, 5, -M_PI/2, zero);

  double sensor_range = 50;
  double zero_noise[2] = {0.3, 0.3};
  std::vector<LandmarkObs> observations = getObservations();
  Map map = getMap();

  pf.updateWeights(sensor_range, zero_noise, observations, map);

  Particle p = pf.particles[0];
//  std::cout << "Map:" << map.landmark_list[0].x_f << " " << map.landmark_list[0].y_f << std::endl;
//  printParticle(p);
}


TEST(ParticleFilter, test_resample) {
  //using the simple scenario in the lectures
  ParticleFilter pf(3);

  double zero[3] = {0.0, 0.0, 0.0};
  pf.init(4, 5, -M_PI/2, zero);

  pf.particles[1].x = 3.98;
  pf.particles[1].y = 5.02;
  pf.particles[2].x = 3.95;
  pf.particles[2].y = 5.05;

  double sensor_range = 50;
  double zero_noise[2] = {0.3, 0.3};
  std::vector<LandmarkObs> observations = getObservations();
  Map map = getMap();

  pf.updateWeights(sensor_range, zero_noise, observations, map);

  pf.resample();
  std::cout << "Map:" << map.landmark_list[0].x_f << " " << map.landmark_list[0].y_f << std::endl;
  printParticle(pf.particles[0]);
  printParticle(pf.particles[1]);
  printParticle(pf.particles[2]);




//  std::vector<double> weights;
//  for(int i = 0; i < pf.particles.size(); i++) {
//      weights.push_back(pf.particles[i].weight);
//  }
//
//  std::random_device rd;
//  std::mt19937 gen(rd());
//  std::discrete_distribution<> dist(std::begin(weights), std::end(weights));
//  std::map<int, int> m;
//  for(int n=0; n<10000; ++n) {
//      ++m[dist(gen)];
//  }
//  for(auto p : m) {
//      std::cout << p.first << " generated " << p.second << " times\n";
//  }

}






