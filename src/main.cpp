//#pragma GCC optimize ("O3")

#include "utils.h"
#include "reader.h"
#include "event_generator.h"

using namespace std;

int main(){

    /* Test tracking
    {
    mt19937 gen(time(0));

    Plane dect1(-8, 4, 4);
    Plane dect2(-5, 4, 4);
    Plane dect3(-2, 4, 4);
    Plane dect3b(0, 4, 4);
    vector<Plane> p_dects;
    p_dects.push_back(dect1);
    p_dects.push_back(dect2);
    p_dects.push_back(dect3);
    p_dects.push_back(dect3b);

    Plane dect4(2, 4, 4);
    Plane dect5(4, 4, 4);
    Plane dect6(8, 4, 4);
    Plane dect6b(10, 4, 4);
    vector<Plane> m_dects;
    m_dects.push_back(dect4);
    m_dects.push_back(dect5);
    m_dects.push_back(dect6);
    m_dects.push_back(dect6b);

    Plane dect7(12, 4, 4);
    Plane dect8(14, 4, 4);
    Plane dect9(18, 4, 4);
    Plane dect9b(20, 4, 4);
    vector<Plane> m2_dects;
    m2_dects.push_back(dect7);
    m2_dects.push_back(dect8);
    m2_dects.push_back(dect9);
    m2_dects.push_back(dect9b);

    Line lp = Line(Vec3D(-0.0158,0.0604,0.0074), Vec3D(0.0,0.0,1.0));
    Line lm1 = Line(Vec3D(-0.0158,0.0604,0.0074), Vec3D(-0.0284,0.0181,27.186));
    Line lm2 = Line(Vec3D(-0.0158,0.0604,0.0074), Vec3D(0.0284,-0.0181,21.617));
    vector<vector<vector<Vec3D>>> meas = measures_from_parameters(lp, lm1, lm2, p_dects, m_dects, m2_dects, 0.0001, 1,
                                                                10, 1, 1.7, 1.602e-19);




    vector<Line> param = parameters_from_measures(meas, p_dects, m_dects, m2_dects, 0.0001,
                                                            10, 1, 1.7, 1.602e-19);
    cout << "Pos\n";
    cout << param[0].str() << endl;
    cout << "Mu+\n";
    cout << param[1].str() << endl;
    cout << "Mu-\n";
    cout << param[2].str() << endl;
    }
    //*/


    /* Test Generator
    {
    double thickness = 5.06e15 * 0.03; //GeV^-1
    double eta = 4* 1e-18; // GeV^3
    double E = 49; // GeV
    Generator g = Generator(E, 0.5, eta, thickness, 1000, 360);
    int n = 10000;
    //cout << n << endl;
    for(int i=0; i<n; i++){
        vector<double> tmp = g.generate_event();


        //for(double dd : tmp)
        //    cout << dd << " ";
        //cout << endl;


        Line mu_pos = Line(Vec3D(tmp[0], tmp[1], tmp[2]), Vec3D(tmp[3], tmp[4], tmp[5]));
        Line mu_neg = Line(Vec3D(tmp[0], tmp[1], tmp[2]), Vec3D(tmp[6], tmp[7], tmp[8]));
        Line mu_pos_2 = trajectory_deviation(mu_pos, 9, 2, 1.7, 1.602e-19);
        Line mu_neg_2 = trajectory_deviation(mu_neg, 9, 2, 1.7, -1.602e-19);
        Plane p = Plane(2);
        Plane p2 = Plane(10);
        Vec3D h1 = p.line_intersection(mu_pos);
        Vec3D h2 = p.line_intersection(mu_neg);
        Vec3D h12 = p2.line_intersection(mu_pos);
        Vec3D h22 = p2.line_intersection(mu_neg);

        cout << sqrt((h1 - h2).norm2()) << " " << sqrt((h12 - h22).norm2()) <<endl;

        //double z = 2;
        //cout << z << " "<< mu_pos_2.point.x + (z-11)*(mu_pos_2.dir.x/mu_pos_2.dir.z) << " " << z << " " << mu_neg_2.point.x + (z-11)*(mu_neg_2.dir.x/mu_neg_2.dir.z) << endl;
        //cout << mu_pos_2.point.x + (z-11)*(mu_pos_2.dir.x/mu_pos_2.dir.z) << " " <<mu_pos_2.point.y + z*(mu_pos_2.dir.y/mu_pos_2.dir.z) << " " << mu_neg_2.point.x + (z-11)*(mu_neg_2.dir.x/mu_neg_2.dir.z) << " "<< mu_neg_2.point.y + z*(mu_neg_2.dir.y/mu_neg_2.dir.z) << endl;


    }


    }
    //*/


    /* Test hit generation */
    {
    mt19937 gen(time(0));

    Plane dect1(-6, 0.1, 0.1);
    Plane dect2(-4, 0.1, 0.1);
    Plane dect3(-2, 0.1, 0.1);
    Plane dect3b(0, 0.1, 0.1);
    vector<Plane> p_dects;
    p_dects.push_back(dect1);
    p_dects.push_back(dect2);
    p_dects.push_back(dect3);
    p_dects.push_back(dect3b);

    Plane dect4(2, 0.1, 0.1);
    Plane dect5(4, 0.1, 0.1);
    Plane dect6(6, 0.1, 0.1);
    Plane dect6b(8, 0.1, 0.1);
    vector<Plane> m_dects;
    m_dects.push_back(dect4);
    m_dects.push_back(dect5);
    m_dects.push_back(dect6);
    m_dects.push_back(dect6b);

    vector<Plane> m_dects2;
    m_dects2.push_back(dect4);
    m_dects2.push_back(dect5);
    m_dects2.push_back(dect6b);

    vector<Plane> m_dects3;
    m_dects3.push_back(dect4);
    m_dects3.push_back(dect5);
    m_dects3.push_back(dect6);
    m_dects3.push_back(Plane(7));
    m_dects3.push_back(Plane(8));

    Plane dect7(12, 0.1, 0.1);
    Plane dect8(14, 0.1, 0.1);
    Plane dect9(16, 0.1, 0.1);
    Plane dect9b(18, 0.1, 0.1);
    vector<Plane> m2_dects;
    m2_dects.push_back(dect7);
    m2_dects.push_back(dect8);
    m2_dects.push_back(dect9);
    m2_dects.push_back(dect9b);

    vector<Plane> m2_dects2;
    m2_dects2.push_back(dect7);
    m2_dects2.push_back(dect8);
    m2_dects2.push_back(dect9b);

    vector<Plane> m2_dects3;
    m2_dects3.push_back(dect7);
    m2_dects3.push_back(dect8);
    m2_dects3.push_back(dect9);
    m2_dects3.push_back(Plane(17));
    m2_dects3.push_back(Plane(18));

    vector<vector<Line>> lines = read_parameters("dataset_param_c.txt");
    int cnt = 0;
    for(vector<Line> event : lines){
        if(cnt++ >= 0) break;
        vector<vector<vector<Vec3D>>> meas = measures_from_parameters(event[0], event[1], event[2], p_dects, m_dects, m2_dects, 0.0002, 1, 9, 2, 1.7, 1.602e-19);
        vector<Line> lines_r = parameters_from_measures(meas, p_dects, m_dects, m2_dects, 0.0002, 9, 2, 1.7, 1.602e-19);
        /*this is error of point
        Vec3D p0 = (Plane(0)).line_intersection(event[0]);
        Vec3D mp0 = (Plane(0)).line_intersection(event[1].dir.z > event[2].dir.z ? event[1] : event[2]);
        Vec3D mn0 = (Plane(0)).line_intersection(event[1].dir.z > event[2].dir.z ? event[2] : event[1]);
        Vec3D p0_r = lines_r[0].point;
        Vec3D mp0_r = ((lines_r[1].dir.z > lines_r[2].dir.z) ? lines_r[1].point : lines_r[2].point);
        Vec3D mn0_r = ((lines_r[1].dir.z > lines_r[2].dir.z) ? lines_r[2].point : lines_r[1].point);
        double dist_error_1 = (p0 - p0_r).norm2();
        double dist_error_2 = (mp0 - mp0_r).norm2();
        double dist_error_3 = (mn0 - mn0_r).norm2();
        cout << dist_error_1 << " " << dist_error_2 << " " << dist_error_3<< endl;
        */

        // this is error of momentum
        Vec3D p_error_0 = event[0].dir - lines_r[3].dir;
        Vec3D p_error_1 = event[0].dir - lines_r[0].dir;
        Vec3D p_error_2 = (event[1].dir.z > event[2].dir.z ? event[1].dir : event[2].dir) - ((lines_r[1].dir.z > lines_r[2].dir.z) ? lines_r[1] : lines_r[2]).dir;
        Vec3D p_error_3 = (event[1].dir.z > event[2].dir.z ? event[2].dir : event[1].dir) - ((lines_r[1].dir.z > lines_r[2].dir.z) ? lines_r[2] : lines_r[1]).dir;

        //cout << p_error_1.norm2() << " " << p_error_0.norm2()<< endl;
        //if(p_error_1.norm2()+p_error_0.norm2() > 1){
        //    cout << p_error_1.norm2() << " " << p_error_0.norm2()<< endl;
        //    cout << event[0].dir.str() << ":  " << lines_r[0].dir.str() << " " << lines_r[3].dir.str()<< endl;
        //    cout << event[1].dir.str() << " " << event[2].dir.str() << ":  " << lines_r[1].dir.str()<< " " << lines_r[2].dir.str()<< endl;
        //}
        //cout << (p_error_0.x*p_error_0.x + p_error_0.y*p_error_0.y) << " " << ((p_error_1.x*p_error_1.x + p_error_1.y*p_error_1.y)) << endl;
        cout <<p_error_1.norm2() << " ";
    }

    for(vector<Line> event : lines){
        if(cnt++ >= 2000) break;
        vector<vector<vector<Vec3D>>> meas = measures_from_parameters(event[0], event[1], event[2], p_dects, m_dects2, m2_dects2, 0.0002, 1, 9, 2, 1.7, 1.602e-19);
        vector<Line> lines_r = parameters_from_measures(meas, p_dects, m_dects2, m2_dects2, 0.0001, 9, 2, 1.7, 1.602e-19);
        Vec3D p_error_1 = event[0].dir - lines_r[0].dir;
        cout <<p_error_1.norm2() << " ";

        meas = measures_from_parameters(event[0], event[1], event[2], p_dects, m_dects, m2_dects, 0.0001, 1, 9, 2, 1.7, 1.602e-19);
        lines_r = parameters_from_measures(meas, p_dects, m_dects, m2_dects, 0.0001, 9, 2, 1.7, 1.602e-19);
        p_error_1 = event[0].dir - lines_r[0].dir;
        cout <<p_error_1.norm2() << " ";

        meas = measures_from_parameters(event[0], event[1], event[2], p_dects, m_dects3, m2_dects3, 0.0001, 1, 9, 2, 1.7, 1.602e-19);
        lines_r = parameters_from_measures(meas, p_dects, m_dects3, m2_dects3, 0.0001, 9, 2, 1.7, 1.602e-19);
        p_error_1 = event[0].dir - lines_r[0].dir;
        cout <<p_error_1.norm2() << " ";

        cout << endl;
    }


    }
    //*/

}
