#include "utils.h"
#include "reader.h"
#include "event_generator.h"

using namespace std;

int main(){
    double thickness = 5.06e15 * 0.03; //GeV^-1
    double eta = 4* 1e-18; // GeV^3
    double E = 49; // GeV
    Generator g = Generator(E, 0.5, eta, thickness, 1000, 360, 10.0);

    int n = 10000;
    for(int i=0; i<n; i++){
        vector<double> tmp = g.generate_event(false);
        //print generated event data
        for(double dd : tmp)
            cout << dd << " ";
        cout << endl;

        /* test hit generation
        Line mu_pos = Line(Vec3D(tmp[0], tmp[1], tmp[2]), Vec3D(tmp[3], tmp[4], tmp[5]));
        Line mu_neg = Line(Vec3D(tmp[0], tmp[1], tmp[2]), Vec3D(tmp[6], tmp[7], tmp[8]));
        Line postr = Line(Vec3D(tmp[0], tmp[1], tmp[2]), Vec3D(tmp[9], tmp[10], tmp[11]));
        Line mu_pos_2 = trajectory_deviation(mu_pos, 9, 2, 1.7, 1.602e-19);
        Line mu_neg_2 = trajectory_deviation(mu_neg, 9, 2, 1.7, -1.602e-19);
        Plane p = Plane(10); //a detector at z = 10
        Vec3D hit1 = p.line_intersection(mu_pos);
        Vec3D hit2 = p.line_intersection(mu_neg);

        cout << sqrt((hit1 - hit2).norm2()) <<endl; //print distance between hits
        //*/
    }
    /* test offset and divergence compensation
    vector<double> tmp = g.generate_event(true);
    Line mu_pos = Line(Vec3D(tmp[0], tmp[1], tmp[2]), Vec3D(tmp[3], tmp[4], tmp[5]));
    Line mu_neg = Line(Vec3D(tmp[0], tmp[1], tmp[2]), Vec3D(tmp[6], tmp[6], tmp[8]));
    Line postr = Line(Vec3D(tmp[0], tmp[1], tmp[2]), Vec3D(tmp[9], tmp[10], tmp[11]));
    vector<Line> event = {postr, mu_pos, mu_neg};
    for(Line l : event) cout << l.str() << " " << l.dir.norm2() << endl;
    cout << endl;
    vector<Line> event_comp = compensate_offset_divergence(event);
    for(Line l : event_comp) cout << l.str() << " " << l.dir.norm2() << endl;
    cout << endl;
    */
}
