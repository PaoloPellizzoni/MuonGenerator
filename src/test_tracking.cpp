#include "utils.h"
#include "reader.h"
#include "event_generator.h"

using namespace std;

int main(){
    mt19937 gen(time(0));

    vector<Plane> p_dects; //it's empty!

    Plane dect4(2, 0.1, 0.1);
    Plane dect5(4, 0.1, 0.1);
    Plane dect6(6, 0.1, 0.1);
    Plane dect6b(8, 0.1, 0.1);
    vector<Plane> m_dects;
    m_dects.push_back(dect4);
    m_dects.push_back(dect5);
    m_dects.push_back(dect6);
    m_dects.push_back(dect6b);

    Plane dect7(12, 0.1, 0.1);
    Plane dect8(14, 0.1, 0.1);
    Plane dect9(16, 0.1, 0.1);
    Plane dect9b(18, 0.1, 0.1);
    vector<Plane> m2_dects;
    m2_dects.push_back(dect7);
    m2_dects.push_back(dect8);
    m2_dects.push_back(dect9);
    m2_dects.push_back(dect9b);

    uniform_real_distribution<double> u(0, 1000000000);
    vector<vector<Line>> lines = read_parameters("dataset_param_c.txt");
    int cnt = 0;

    for(vector<Line> event : lines){
        if(cnt++>= 2000) break;
        vector<vector<vector<Vec3D>>> meas = measures_from_parameters(event[0], event[1], event[2], p_dects, m_dects, m2_dects, 0.0001, 1, 9, 2, 1.7, 1.602e-19, mt19937((int)u(gen)));
        vector<Line> lines_r = parameters_from_measures(meas, p_dects, m_dects, m2_dects, 0.0001, 9, 2, 1.7, 1.602e-19);
        Vec3D p_error_1 = event[0].dir - lines_r[0].dir;
        cout <<sqrt(p_error_1.norm2()) << endl;
    }

}
