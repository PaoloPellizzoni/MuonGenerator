#ifndef GEOM_H_
#define GEOM_H_

#include <bits/stdc++.h>

struct Vec3D{
    double x, y, z;
    Vec3D(){
    }

    Vec3D(double _x, double _y, double _z){
        x = _x;
        y = _y;
        z = _z;
    }

    Vec3D operator + (const Vec3D &o){
        return Vec3D(x + o.x, y + o.y, z + o.z);
    }
    Vec3D operator - (const Vec3D &o){
        return Vec3D(x - o.x, y - o.y, z - o.z);
    }
    Vec3D operator * (const double k){
        return Vec3D(x * k, y * k, z * k);
    }

    Vec3D cross(const Vec3D &o){
        return Vec3D( y*o.z- z*o.y,  z*o.x- x*o.z,  x*o.y- y*o.x);
    }
    double dot(const Vec3D &o){
        return x*o.x + y*o.y + z*o.z;
    }
    double norm2(){
        return x*x + y*y + z*z;
    }
    Vec3D normalize(){
        return (*this) * (1/sqrt(this->norm2()));
    }

    std::string str(){
        return "("+std::to_string(x)+" "+std::to_string(y)+" "+std::to_string(z)+")";
    }
};

struct Line{
    Vec3D point, dir;

    Line(){}

    Line(Vec3D _p, Vec3D _d){
        point = _p;
        dir = _d;
    }

    Line(std::vector<double> theta){
        point = Vec3D(theta[0], theta[1], 0);
        dir = Vec3D(theta[2], theta[3], 1);
    }

    Vec3D line_quasi_intersection(Line o){
        Vec3D u = dir.normalize();
        Vec3D v = o.dir.normalize();
        Vec3D w = o.point - point;
        double a = u.dot(u);
        double b = u.dot(v);
        double c = v.dot(v);
        double d = u.dot(w);
        double e = v.dot(w);
        double D = a*c - b*b;
        if(D < 1e-20) // parallel lines
            return point;
        double sc = (b*e - c*d)/D;
        double tc = (a*e - b*d)/D;
        Vec3D s = point - u * sc;
        Vec3D t = o.point - v * tc;
        Vec3D mid = (s + t) * 0.5;
        //std::cout << s.str() << " " << t.str() << " " <<mid.str() << " " << std::endl;
        return mid;
    }

    std::string str(){
        return "["+point.str()+" + k*" + dir.str()+"]";
    }
};

struct Plane{
    double z, height, width;

    Plane(double _z){
        z = _z;
        height = 0;
        width = 0;
    }

    Plane(double _z, double _h, double _w){
        z = _z;
        height = _h;
        width = _w;
    }

    Vec3D line_intersection(Line &l){
        double d = (z - l.point.z)/(l.dir.z);
		double x = l.point.x + l.dir.x*d;
		double y = l.point.y + l.dir.y*d;
		return Vec3D(x,y,z);
    }
};

struct Plane3D{
    Vec3D point, normal;
    double height, width;

    Plane3D(Vec3D _p, Vec3D _n){
        point = _p;
        normal = _n.normalize();
    }

    Plane3D(Vec3D _p, Vec3D _n, double _h, double _w){
        point = _p;
        normal = _n.normalize();
        height = _h;
        width = _w;
    }

    Vec3D line_intersection(Line &l){
        double d = ((point - l.point).dot(normal))/(l.dir.dot(normal));
		double x = l.point.x + l.dir.x*d;
		double y = l.point.y + l.dir.y*d;
		double z = l.point.z + l.dir.z*d;
		return Vec3D(x,y,z);
    }
};


#endif
