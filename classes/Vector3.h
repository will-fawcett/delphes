#ifndef Vector3_h
#define Vector3_h

#include<iostream>
#include<cmath>
#include<string>
#include<vector>

class vector3
{
  friend ostream & operator<<(ostream &os, const vector3 &vec3){
    os<<"("<<vec3.x<<","<<vec3.y<<","<<vec3.z<<")"<<endl;
    return os;
  }

private:
	double x,y,z;

public:

  // Default constructor
	vector3(){x=y=z=0;}

  // Paramitarised constructor
	vector3(double xin, double yin, double zin){x=xin; y=yin= z=zin;}

  // Destructor
	~vector3(){}

	double getx() const {return x;}
	double gety() const {return y;}
	double getz() const {return z;}

  // Function to print out vector (redundant with friend)
	void show() const {cout<<"("<<x<<","<<y<<","<<z<<")"<<endl;}

  // Function to add scalar to each component
	void addscalar(double s) {x+=s; y=+s; z+=s;}

  // Function to overload + for vector sum
	vector3 operator+(const vector3 &vec3)
	{
		vector3 temp(x+vec3.getx(), y+vec3.gety(), z+vec3.getz());
		return temp;
	}

  // Function to overload * for dot product
	double operator*(vector3 &vec3)
	{
		double dot(x*vec3.getx() + y*vec3.gety() + z*vec3.getz());
		return dot;
	}

};

#endif /* Vector3_h */
