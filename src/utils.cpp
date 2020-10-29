#include "utils.h"

using namespace std;

void subset_creator_util(int i, vector<vector<Vec3D>>& measures, vector<vector<Vec3D>>& ans, vector<Vec3D>& tmp){
    if(i == measures.size()){
        vector<Vec3D> v = vector<Vec3D>(tmp.size());
        for(int i=0; i<v.size(); i++)
            v[i] = tmp[i];
        ans.push_back(v);
        return;
    }
    for(Vec3D p : measures[i]){
        tmp.push_back(p);
        subset_creator_util(i+1, measures, ans, tmp);
        tmp.pop_back();
    }
}

vector<vector<Vec3D>> subset_creator(vector<vector<Vec3D>>& measures){
    vector<vector<Vec3D>> ans = vector<vector<Vec3D>>();
    vector<Vec3D> tmp = vector<Vec3D>();
    subset_creator_util(0, measures, ans, tmp);
    return ans;
}




Line trajectory_deviation(Line trj, double z_start, double L, double B, double q){
    double px = trj.dir.x;  //momentum in x
    double py = trj.dir.y;  //momentum in y
    double pz = trj.dir.z;  //momentum in z
    double p_ort = sqrt(px*px + pz*pz); //momentum orthogonal to B
    double a = atan2(px, pz);           //angle of momentum with z axis
    double r = p_ort*5.36e-19/(q*B);             //curvature radius
    Plane start_plane = Plane(z_start); //the plane between the no-field and the field areas
    Vec3D start_point = start_plane.line_intersection(trj);         //the point where the particle enters the field
    double tmp = atan2( (q > 0 ? 1 : -1)*(L + r*sin(a)) , sqrt(r*r - (L+r*sin(a))*(L+r*sin(a)) ) );
    Vec3D dirf = Vec3D(p_ort*sin(tmp), py, p_ort*cos(tmp));         //new momentum
    double theta = acos((px*dirf.x + pz*dirf.z)/(sqrt(px*px + pz*pz)*sqrt(dirf.x*dirf.x + dirf.z*dirf.z)));

    double dx = r*cos(a) - (q > 0 ? 1 : -1)*sqrt(r*r - (L + r*sin(a))*(L + r*sin(a)) ); //delta x from entry to exit point
    //double dy = r*(asin(L/r + sin(a)) -a)*(py/p_ort);                                   //delta y from entry to exit point
    double dy = fabs(r)*theta*(py/p_ort);
    double dz = L;                                                                      //delta z from entry to exit point
    Vec3D dp = Vec3D(dx, dy, dz);                                                       //delta position

    return Line(start_point + dp, dirf);                            // Line(exit point, new momentum)
}

double momentum_from_deviation(Line in, Line out, double z_start, double L, double B, double q){
    double px_in = in.dir.x; //entering momentum in x
    double pz_in = in.dir.z; //enterng momentum in z
    double px_out = out.dir.x; //exiting momentum in x
    double pz_out = out.dir.z; //exiting momentum in z
    Plane start_plane = Plane(z_start);
    Plane end_plane = Plane(z_start+L);
    Vec3D start_point = start_plane.line_intersection(in);
    Vec3D end_point = end_plane.line_intersection(out);
    double theta = acos((px_in*px_out + pz_in*pz_out)/(sqrt(px_in*px_in + pz_in*pz_in)*sqrt(px_out*px_out + pz_out*pz_out)));
    double dx = end_point.x - start_point.x;
    double dy = end_point.y - start_point.y;
    double dz = end_point.z - start_point.z;
    double d = sqrt(dx*dx + dz*dz);
    double a = atan2(px_in, pz_in);
    double r = d/(2*sin(theta/2));
    double p_ort = r*q*B/5.36e-19;
    double py = (p_ort*dy)/(r*(asin(L/r + sin(a)) -a));
    double p = sqrt(py*py + p_ort*p_ort);
    return p;
}

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
                                                        ){
    // this function creates the synthetic dataset of measures
    vector<vector<vector<Vec3D>>> ret = vector<vector<vector<Vec3D>>>(3);
    poisson_distribution<int> poiss(poisson_noise);
    normal_distribution<double> gauss(0, det_std);
    // compute positron measures
    vector<vector<Vec3D>> positron_measures = vector<vector<Vec3D>>();
    for(Plane& det : positron_det){
        uniform_real_distribution<double> unif_h(-det.height/2, det.height/2);
        uniform_real_distribution<double> unif_w(-det.width/2, det.width/2);
        Vec3D vertical = Vec3D(0,1,0);
        Vec3D horizontal = Vec3D(1,0,0);

        vector<Vec3D> det_measures = vector<Vec3D>();
        Vec3D measure = det.line_intersection(positron_trj);
        measure = measure + vertical*gauss(gen) + horizontal*gauss(gen);   //add smearing
        det_measures.push_back(measure);

        int n = poiss(gen);
        for(int i=0; i<n; i++){
            measure = Vec3D(0,0,det.z) + vertical*unif_h(gen) + horizontal*unif_w(gen); // uniform
            det_measures.push_back(measure);
        }
        positron_measures.push_back(det_measures);
    }
    // compute muon measures before the field
    vector<vector<Vec3D>> muon_measures = vector<vector<Vec3D>>();

    for(Plane& det : muon_det){
        uniform_real_distribution<double> unif_h(-det.height/2, det.height/2);
        uniform_real_distribution<double> unif_w(-det.width/2, det.width/2);
        Vec3D vertical = Vec3D(0,1,0);
        Vec3D horizontal = Vec3D(1,0,0);

        vector<Vec3D> det_measures = vector<Vec3D>();
        Vec3D measure = det.line_intersection(mu_pos_trj);
        measure = measure + vertical*gauss(gen) + horizontal*gauss(gen);   //add smearing
        det_measures.push_back(measure);
        measure = det.line_intersection(mu_neg_trj);
        measure = measure + vertical*gauss(gen) + horizontal*gauss(gen);   //add smearing
        det_measures.push_back(measure);

        int n = poiss(gen);
        for(int i=0; i<n; i++){
            measure = Vec3D(0,0,det.z) + vertical*unif_h(gen) + horizontal*unif_w(gen); // uniform
            det_measures.push_back(measure);
        }
        muon_measures.push_back(det_measures);
    }

    // compute muon measures after the field
    vector<vector<Vec3D>> muon_measures_2 = vector<vector<Vec3D>>();
    Line mu_pos_trj_2 = trajectory_deviation(mu_pos_trj, z_start, L, B, q);
    Line mu_neg_trj_2 = trajectory_deviation(mu_neg_trj, z_start, L, B, -q);

    for(Plane& det : muon_det_2){
        uniform_real_distribution<double> unif_h(-det.height/2, det.height/2);
        uniform_real_distribution<double> unif_w(-det.width/2, det.width/2);
        Vec3D vertical = Vec3D(0,1,0);
        Vec3D horizontal = Vec3D(1,0,0);

        vector<Vec3D> det_measures = vector<Vec3D>();
        Vec3D measure = det.line_intersection(mu_pos_trj_2);
        measure = measure + vertical*gauss(gen) + horizontal*gauss(gen);   //add smearing
        det_measures.push_back(measure);
        measure = det.line_intersection(mu_neg_trj_2);
        measure = measure + vertical*gauss(gen) + horizontal*gauss(gen);   //add smearing
        det_measures.push_back(measure);

        int n = poiss(gen);
        for(int i=0; i<n; i++){
            measure = Vec3D(0,0,det.z) + vertical*unif_h(gen) + horizontal*unif_w(gen); // uniform
            det_measures.push_back(measure);
        }
        muon_measures_2.push_back(det_measures);
    }


    ret[0] = positron_measures;
    ret[1] = muon_measures;
    ret[2] = muon_measures_2;
    return ret;
}


vector<Line> parameters_from_measures(vector<vector<vector<Vec3D>>>& measures,
                                                vector<Plane>& positron_det,
                                                vector<Plane>& muon_det,
                                                vector<Plane>& muon_det_2,
                                                double det_std,
                                                double z_start, double L, double B, double q
                                                ){
    // this funtction computes the parameters of the particles from the measures
    vector<Line> ret = vector<Line>(3);
    vector<double> best_param(5);
    if(positron_det.size()>0){
        //optimization for positron  trj
        vector<vector<Vec3D>> possible_pos_meas = subset_creator(measures[0]);
        best_param[4] = 1e30;
        for(vector<Vec3D> pos_meas : possible_pos_meas){
            OptimizerGM opt = OptimizerGM(pos_meas, positron_det, det_std);
            vector<double> ans = opt.optimize();
            if(ans[4] < best_param[4])
                best_param = ans;
        }
    }
    vector<double> positron_param = vector<double>(4);
    for(int i=0; i<4; i++)
        positron_param[i] = best_param[i];

    //optimization for muon 1 trjs
    vector<vector<Vec3D>> possible_mu_meas = subset_creator(measures[1]);
    best_param = vector<double>(6);
    best_param[4] = 1e30;
    vector<double> sec_param = vector<double>(6);
    sec_param[4] = 1e30;
    for(vector<Vec3D> mu_meas : possible_mu_meas){
        OptimizerGM opt = OptimizerGM(mu_meas, muon_det, det_std);
        vector<double> ans = opt.optimize();
        double sumY = 0;
        for(Vec3D hitt : mu_meas)
            sumY += hitt.y;
        if(ans[4] < best_param[4]){
            sec_param = best_param;
            best_param = ans;
            best_param.push_back(sumY);
        }
        else if(ans[4] >= best_param[4] && ans[4] < sec_param[4]){
            sec_param = ans;
            sec_param.push_back(sumY);
        }
    }

    double mu_y_1 = best_param[5];
    vector<double> mu_param1 = vector<double>(4);
    for(int i=0; i<4; i++)
        mu_param1[i] = best_param[i];
    double mu_y_2 = sec_param[5];
    vector<double> mu_param2 = vector<double>(4);
    for(int i=0; i<4; i++)
        mu_param2[i] = sec_param[i];


    // optimization for muon 2 trjs
    possible_mu_meas = subset_creator(measures[2]);
    best_param = vector<double>(5);
    best_param[4] = 1e30;
    sec_param = vector<double>(5);
    sec_param[4] = 1e30;
    for(vector<Vec3D> mu_meas : possible_mu_meas){
        OptimizerGM opt = OptimizerGM(mu_meas, muon_det_2, det_std);
        vector<double> ans = opt.optimize();
        double sumY = 0;
        for(Vec3D hitt : mu_meas)
            sumY += hitt.y;
        if(ans[4] < best_param[4]){
            sec_param = best_param;
            best_param = ans;
            best_param.push_back(sumY);
        }
        else if(ans[4] >= best_param[4] && ans[4] < sec_param[4]){
            sec_param = ans;
            sec_param.push_back(sumY);
        }
    }

    double mu2_y_1 = best_param[5];
    vector<double> mu2_param1 = vector<double>(4);
    for(int i=0; i<4; i++)
        mu2_param1[i] = best_param[i];
    double mu2_y_2 = sec_param[5];
    vector<double> mu2_param2 = vector<double>(4);
    for(int i=0; i<4; i++)
        mu2_param2[i] = sec_param[i];



    Line mu_up_trj_1 = Line(mu_y_1 > mu_y_2 ? mu_param1 : mu_param2);
    Line mu_down_trj_1 = Line(mu_y_1 > mu_y_2 ? mu_param2 : mu_param1);
    Line mu_up_trj_2 = Line(mu2_y_1 > mu2_y_2 ? mu2_param1 : mu2_param2);
    Line mu_down_trj_2 = Line(mu2_y_1 > mu2_y_2 ? mu2_param2 : mu2_param1);

    double momentum_up = momentum_from_deviation(mu_up_trj_1, mu_up_trj_2, z_start, L, B, q);
    double momentum_down = momentum_from_deviation(mu_down_trj_1, mu_down_trj_2, z_start, L, B, q);
    //double tot_momentum = momentum_up + momentum_down;

    double mult_up = momentum_up/sqrt(mu_up_trj_1.dir.norm2());
    double mult_down = momentum_down/sqrt(mu_down_trj_1.dir.norm2());
    //double mult_tot = tot_momentum/sqrt(positron_param[2]*positron_param[2] + positron_param[3]*positron_param[3] + 1);


    Vec3D dir_up = Vec3D(mu_up_trj_1.dir.x*mult_up, mu_up_trj_1.dir.y*mult_up, mu_up_trj_1.dir.z*mult_up);
    Vec3D dir_down = Vec3D(mu_down_trj_1.dir.x*mult_down, mu_down_trj_1.dir.y*mult_down, mu_down_trj_1.dir.z*mult_down);
    // ignore first detc
    Vec3D dir_pos = dir_up + dir_down;

    // this used the positron detectors
    //Vec3D dir_pos_dec = Vec3D(positron_param[2], positron_param[3], 1)* (sqrt(dir_pos.norm2())/sqrt((Vec3D(positron_param[2], positron_param[3], 1)).norm2()));


    Vec3D point_up = Vec3D(mu_up_trj_1.point.x, mu_up_trj_1.point.y, 0);
    Vec3D point_down = Vec3D(mu_down_trj_1.point.x, mu_down_trj_1.point.y, 0);
    // this used the positron detectors
    //Vec3D point_pos = Vec3D(positron_param[0], positron_param[1], 0);


    Line mu_up_trj_tmp = Line(point_up, dir_up);
    Line mu_down_trj_tmp = Line(point_down, dir_down);

    Vec3D point_pos = mu_up_trj_tmp.line_quasi_intersection(mu_down_trj_tmp);
    Line pos_trj = Line(point_pos, dir_pos);
    Line mu_up_trj = Line(point_pos, dir_up);
    Line mu_down_trj = Line(point_pos, dir_down);

    ret[0] = pos_trj;
    ret[1] = mu_up_trj;
    ret[2] = mu_down_trj;

    return ret;
}


vector<Line> compensate_offset_divergence(vector<Line> event){
    // needs (positron, muon1, muon2) as Lines
    double alpha = atan2(event[0].dir.y, event[0].dir.z); // around x
    double beta = atan2(event[0].dir.x, hypot(event[0].dir.z, event[0].dir.x) );
    double cosA = cos(alpha);
    double sinA = sin(alpha);
    double cosB = cos(beta);
    double sinB = sin(beta);
    double z = event[0].point.z;
    vector<Line> ret;
    for(Line l : event){
        double px = l.dir.x;
        double py = l.dir.y;
        double pz = l.dir.z;
        double px_n = px*cosB - py*sinA*sinB - pz*cosA*sinB;
        double py_n = py*cosA - pz*sinA;
        double pz_n = px*sinB + py*sinA*cosB + pz*cosA*cosB;
        ret.push_back(Line(Vec3D(0,0,z), Vec3D(px_n, py_n, pz_n)));
    }
    return ret;
}
