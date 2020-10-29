#ifndef EV_G_H
#define EV_G_H

#include <bits/stdc++.h>

using namespace std;

class Generator{
private:
    static constexpr double a = 1/137.0;
    static constexpr double m_e = 0.511e-3;
    static constexpr double m_u = 0.105683;
    static constexpr double pi = acos(-1.0);
    static constexpr double x0 = 5.06e15 * 0.3528;

    mt19937 gen;
    uniform_real_distribution<double> unif;
    normal_distribution<double> gauss;
    vector<double> theta_pdf;
    vector<double> depth_pdf;
    discrete_distribution<int> theta_distribution;
    discrete_distribution<int> depth_distribution;
    double d_theta, dz, genertorDist;
    double nominal_energy, energy_std, elec_density;
    int depth_subd, theta_subd;

public:
    Generator(double _nominal_energy, double _energy_std, double _elec_density, double thickness, int _depth_subd, int _theta_subd, double _genertorDist, mt19937 _gen = mt19937(time(0))){
        gen = _gen;
        unif = uniform_real_distribution<double>(0, 1);
        gauss = normal_distribution<double>(0, 1);
        nominal_energy = _nominal_energy;
        energy_std = _energy_std;
        depth_subd = _depth_subd;
        theta_subd = _theta_subd;
        elec_density = _elec_density;
        genertorDist = _genertorDist;
        theta_pdf = vector<double>(theta_subd);
        depth_pdf = vector<double>(depth_subd);
        d_theta = pi/theta_subd;
        dz = thickness/depth_subd;

    }

    vector<double> generate_event(bool rot = false){
        double E = nominal_energy + energy_std*gauss(gen);
        double s0 = 2*E*m_e;

        // generate the depth of the event
        double s;
        double sum = 0, pos = 0;
        for(int i=0; i<depth_subd; i++){
            s = s0 * exp(-pos/x0);
            double sigma = a*a*pi/s * sqrt( (1 - m_u*m_u/s) / (1 - m_e*m_e/s)) * (1 + 4/s*(m_e*m_e + m_u*m_u) + (1- 4*m_e*m_e/s)*(1- 4*m_u*m_u/s)/3 );
            depth_pdf[i] = (1-sum) * (elec_density * dz * sigma);
            sum += depth_pdf[i];
            pos += dz;
        }

        depth_distribution = discrete_distribution<int>(depth_pdf.begin(), depth_pdf.end());
        int depth_id = depth_distribution(gen);
        double depth = (depth_id + unif(gen))*dz;

        // now with the depth generate the angle
        s = s0 * exp(-depth/x0);
        double sigma = a*a*pi/s * sqrt( (1 - m_u*m_u/s) / (1 - m_e*m_e/s)) * (1 + 4/s*(m_e*m_e + m_u*m_u) + (1- 4*m_e*m_e/s)*(1- 4*m_u*m_u/s)/3 );
        for(int j=0; j<theta_subd; j++){
            double theta = d_theta*(j+0.5);
            theta_pdf[j] = ( 0.5*pi*a*a/s * sqrt( (1 - m_u*m_u/s) / (1 - m_e*m_e/s)) * (1 + 4/s*(m_e*m_e + m_u*m_u) + (1- 4*m_e*m_e/s)*(1- 4*m_u*m_u/s)*cos(theta)*cos(theta)) ) * sin(theta) * d_theta;
        }
        theta_distribution = discrete_distribution<int>(theta_pdf.begin(), theta_pdf.end());
        int theta_id = theta_distribution(gen);
        double theta = (theta_id + unif(gen))*d_theta;
        double phi = unif(gen)*2*pi; // uniform on [0, 2pi]

        // generate offeset position uniform on disc
        double x = 1, y = 1;
        while(x*x + y*y > 0.0001){
            x = -0.01 +0.02*unif(gen);
            y = -0.01 +0.02*unif(gen);
        }

        double gamma = sqrt(s)/(2*m_e);
        double beta = sqrt(1 - 4*m_e*m_e/s);
        // in mass center
        double pz_cm = 0.5*cos(theta)*sqrt(s - 4 * m_u * m_u);
        double px_cm = 0.5*sin(theta)*cos(phi)*sqrt(s - 4 * m_u * m_u);
        double py_cm = 0.5*sin(theta)*sin(phi)*sqrt(s - 4 * m_u * m_u);

        // in lab frame
        double pz_1 = gamma*(pz_cm + beta*(sqrt(s)/2));
        double px_1 = px_cm;
        double py_1 = py_cm;

        double pz_2 = gamma*(-pz_cm + beta*(sqrt(s)/2));
        double px_2 = -px_cm;
        double py_2 = -py_cm;

        vector<double> ans(12);
        if(s < 4 * m_u * m_u){ // failed to generate event
            //cout << "Error, energy "<< s <<" too low! " << 4*m_u*m_u << "\n";
            return ans; //return 0s vector
        }


        ans[0] = x;
        ans[1] = y;
        ans[2] = depth/5.06e15; //performed calcs in GeV-1, transform in meters

        //rotation
        if(rot){
            double alpha = unif(gen)*0.001; //around x
            double beta = unif(gen)*0.001;  //around y
            double cosA = cos(alpha);
            double sinA = sin(alpha);
            double cosB = cos(beta);
            double sinB = sin(beta);

            ans[0] += genertorDist*tan(beta); //add offset on x due to divergence
            ans[1] += genertorDist*tan(alpha); //add offset on y due to divergence
            //muon 1
            ans[3] = px_1*cosB + pz_1*sinB;
            ans[4] = px_1*sinA*sinB - pz_1*sinA*cosB + py_1*cosA;
            ans[5] = -px_1*cosA*sinB + pz_1*cosA*cosB + py_1*sinA;
            //muon 2
            ans[6] = px_2*cosB + pz_2*sinB;
            ans[7] = px_2*sinA*sinB - pz_2*sinA*cosB + py_2*cosA;;
            ans[8] = -px_2*cosA*sinB + pz_2*cosA*cosB + py_2*sinA;
            //positron
            double pz_3 = pz_1 + pz_2;
            ans[9] = pz_3*sinB;
            ans[10] = - pz_3*sinA*cosB;
            ans[11] = pz_3*cosA*cosB;

        }
        else{
            //muon 1
            ans[3] = px_1;
            ans[4] = py_1;
            ans[5] = pz_1;
            //muon 2
            ans[6] = px_2;
            ans[7] = py_2;
            ans[8] = pz_2;
            //positron
            ans[9] = 0;
            ans[10] = 0;
            ans[11] = pz_1 + pz_2;
        }
        return ans;

    }

};

#endif
