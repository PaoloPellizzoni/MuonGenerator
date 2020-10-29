#ifndef UTILS_H_
#define UTILS_H_

#include <bits/stdc++.h>
#include "geometry.h"
#include "optimizer.h"

using namespace std;

Line trajectory_deviation(Line trj, double z_start, double L, double B, double q);
double momentum_from_deviation(Line in, Line out, double z_start, double L, double B, double q);

vector<vector<vector<Vec3D>>> measures_from_parameters( Line& positron_trj,
                                                        Line& mu_pos_trj,
                                                        Line& mu_neg_trj,
                                                        vector<Plane>& positron_det,
                                                        vector<Plane>& muon_det,
                                                        vector<Plane>& muon_det_2,
                                                        double det_std,
                                                        double poisson_noise,
                                                        double z_start, double L, double B, double q,
                                                        mt19937 gen
                                                    );

vector<Line> parameters_from_measures(vector<vector<vector<Vec3D>>>& measures,
                                                vector<Plane>& positron_det,
                                                vector<Plane>& muon_det,
                                                vector<Plane>& muon_det_2,
                                                double det_std,
                                                double z_start, double L, double B, double q
                                            );

vector<Line> compensate_offset_divergence(vector<Line> event);


#endif
