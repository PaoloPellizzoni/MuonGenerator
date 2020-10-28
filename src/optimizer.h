#ifndef OPT_H
#define OPT_H
using namespace std;
class OptimizerGM{
private:
    vector<Vec3D>* _y;
    vector<Plane>* _detectors;
    double _detc_std;
    int n;

    vector<vector<double>> multiply(vector<vector<double>> mat1, vector<vector<double>> mat2){
        int m1 = mat1.size();
        int m2 = mat1[0].size();
        int n2 = mat2[0].size();
        vector<vector<double>> ret = vector<vector<double>>(m1);
        for(int i=0; i<m1; i++){
            ret[i] = vector<double>(n2);
            for(int j=0; j<n2; j++){
                for(int k=0; k<m2; k++)
                    ret[i][j] += mat1[i][k]*mat2[k][j];
            }
        }
        return ret;
    }

    vector<vector<double>> transpose(vector<vector<double>> mat1){
        int m1 = mat1.size();
        int m2 = mat1[0].size();
        vector<vector<double>> ret = vector<vector<double>>(m2);
        for(int i=0; i<m2; i++){
            ret[i] = vector<double>(m1);
            for(int j=0; j<m1; j++){
                ret[i][j] = mat1[j][i];
            }
        }
        return ret;
    }

    double log_likelihood(vector<double> &theta, vector<Vec3D> &y, vector<Plane> &detectors, double detc_std){
        vector<Vec3D> y0;
        Line trj(theta);
        for(int i=0; i<detectors.size(); i++){
            y0.push_back(detectors[i].line_intersection(trj));
        }

        double ans = 0;
        //for each detector compute its likelihood
        for(int i=0; i< y.size(); i++){
            double distance2 = (y[i] - y0[i]).norm2();
            double log_p = log_normal_pdf(detc_std, distance2);
            ans += log_p;
        }

        return ans;
    }

    double log_prior(vector<double> &theta, int cas){
        return 0;
    }

    double log_normal_pdf(double std, double x){
        return (- x / (std*std));
    }

    double minus_log_posterior(vector<double> &theta, vector<Vec3D> &y, vector<Plane> &detectors, double detc_std){
        return - log_likelihood(theta, y, detectors, detc_std) - log_prior(theta, 0);
    }

public:
    OptimizerGM(vector<Vec3D> &y, vector<Plane> &detectors, double detc_std){
        _y = &y;
        _detectors = &detectors;
        _detc_std = detc_std;
        n = _y->size();
    }

    vector<double> optimize(){
        vector<vector<double>> phi = vector<vector<double>>(n); // linear model
        for(int i=0; i<n; i++){
            phi[i] = vector<double>(2);
            phi[i][0] = 1;
            phi[i][1] = (*_detectors)[i].z;
        }
        vector<vector<double>> sigma = vector<vector<double>>(n); //covariance matrix
        for(int i=0; i<n; i++){
            sigma[i] = vector<double>(n);
            sigma[i][i] = 1/(_detc_std*_detc_std);
        }

        vector<vector<double>> tmp = multiply(multiply(transpose(phi), sigma), phi);
        vector<vector<double>> eta(2);
        eta[0] = vector<double>(2);
        eta[1] = vector<double>(2);
        // invert tmp
        double d = tmp[0][0]*tmp[1][1] - tmp[0][1]*tmp[1][0];
        eta[0][0] = tmp[1][1]/d;
        eta[0][1] = -tmp[1][0]/d;
        eta[1][0] = -tmp[0][1]/d;
        eta[1][1] = tmp[0][0]/d;

        tmp = multiply(multiply(eta, transpose(phi)), sigma);

        vector<vector<double>> y_x(n);
        vector<vector<double>> y_y(n);
        for(int i=0; i<n; i++){
            y_x[i] = vector<double>(1);
            y_x[i][0] = (*_y)[i].x;
            y_y[i] = vector<double>(1);
            y_y[i][0] = (*_y)[i].y;
        }

        vector<vector<double>> theta_hat_x = multiply(tmp, y_x);
        vector<vector<double>> theta_hat_y = multiply(tmp, y_y);

        vector<double> ans(4);
        ans[0] = theta_hat_x[0][0];
        ans[1] = theta_hat_y[0][0];
        ans[2] = theta_hat_x[1][0];
        ans[3] = theta_hat_y[1][0];
        ans.push_back(minus_log_posterior(ans, *_y, *_detectors, _detc_std));

        return ans;
    }

};

#endif
