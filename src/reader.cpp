#include "geometry.h"

using namespace std;

vector<vector<Line>> read_parameters(string filename){
    ifstream fin(filename);
    //freopen(filename.c_str(), "r", stdin);
    int n;
    double x,y,z;
    double px,py,pz;
    fin >> n;
    vector<vector<Line>> ret(n);
    for(int i=0; i<n; i++){
        fin >> x >> y >> z;
        Vec3D p = Vec3D(x, y, z);
        fin >> px >> py >> pz;
        Vec3D l = Vec3D(px, py, pz);
        Line mu1 = Line(p, l);
        fin >> px >> py >> pz;
        l = Vec3D(px, py, pz);
        Line mu2 = Line(p, l);
        fin >> px >> py >> pz;
        l = Vec3D(px, py, pz);
        Line po = Line(p, l);
        vector<Line> tmp(3);
        tmp[0] = po;
        tmp[1] = mu1;
        tmp[2] = mu2;
        ret[i] = tmp;
    }
    //fin.close();
    return ret;
}

vector<vector<vector<vector<Vec3D>>>> read_measures(string filename){
    ifstream fin(filename);
    int n, d, k;
    int d1, d2, d3;
    double x,y,z;
    fin >> n;
    vector<vector<vector<vector<Vec3D>>>> ret(n);
    for(int i=0; i<n; i++){
        vector<vector<vector<Vec3D>>> ret_i(3);
        fin >> d1;
        vector<vector<Vec3D>> pos(d1);
        for(int j=0; j<d1; j++){
            fin >> k;
            for(int jj=0; jj<k; jj++){
                fin >> x >> y >> z;
                pos[j].push_back(Vec3D(x,y,z));
            }
        }
        fin >> d2;
        vector<vector<Vec3D>> mu1(d2);
        for(int j=0; j<d2; j++){
            fin >> k;
            for(int jj=0; jj<k; jj++){
                fin >> x >> y >> z;
                mu1[j].push_back(Vec3D(x,y,z));
            }
        }
        fin >> d3;
        vector<vector<Vec3D>> mu2(d3);
        for(int j=0; j<d3; j++){
            fin >> k;
            for(int jj=0; jj<k; jj++){
                fin >> x >> y >> z;
                mu2[j].push_back(Vec3D(x,y,z));
            }
        }
        ret_i[0] = pos;
        ret_i[1] = mu1;
        ret_i[2] = mu2;
        ret[i] = ret_i;
    }
    //fin.close();
    return ret;
}
