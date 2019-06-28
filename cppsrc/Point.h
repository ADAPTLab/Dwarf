#ifndef POINT_H
#define POINT_H

#include <vector>
#include <string>
#include <stdlib.h>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <cstring>
#include <cmath>

#include "cereal/types/unordered_map.hpp"
#include "cereal/types/memory.hpp"
#include "cereal/archives/binary.hpp"
#include "cereal/types/vector.hpp"

#define MAXSTRINGLEN 8
typedef double PointEle;
static const int MaxEles = 3;
#define distanceEuclidean(A,B)	(sqrt(	\
		((A).point_arr[0] - (B).point_arr[0])	\
	*	((A).point_arr[0] - (B).point_arr[0])	\
	+		\
		((A).point_arr[1] - (B).point_arr[1])	\
	*	((A).point_arr[1] - (B).point_arr[1])	\
	+		\
		((A).point_arr[2] - (B).point_arr[2])	\
	*	((A).point_arr[2] - (B).point_arr[2])	\
))
class Point
{
//private:
public:
PointEle point_arr[MaxEles];
int NumEles;
//public:
	Point()
	{
		NumEles = 0;
	}
	template <class Archive>
	void serialize(Archive &ar)
	{
		ar(CEREAL_NVP(point_arr),NumEles);
	}
	int GetSize();
	int Size();
	int InitializePoint();
	int InitializePoint(double val);
	PointEle GetPointEleAtIndex(int index);
	//double GetEleAtIndex(int index);
	inline double GetEleAtIndex(int index)
	{
		if(index >= 0 && index < NumEles)
			return this->point_arr[index];
		throw std::out_of_range("\n**index out of bounds**\n");
	}
	//double euclidDist(Point &p);
	void InsertPointEleAt(PointEle &pe, int index);
	void InsertInteger(int val);
	void InsertBoolean(bool val);
	void InsertString(std::string val);
	void InsertFloat(double val);
	void SetIntegerAt(int index, int val);
	void SetBooleanAt(int index, bool val);
	void SetStringAt(int index, std::string val);
	void SetFloatAt(int index, double val);
	void SetDoubleAt(int index, double val);
	PointEle operator[] (int pos);
	Point& operator+= (Point &p);
	Point operator+ (Point &p);
	Point operator- (Point p);
	Point& operator-= (Point &p);
	Point operator/ (int p);
	Point operator* (double val);
	Point operator/ (double val);
	bool operator== (Point p);
};
std::ostream& operator<< (std::ostream& out, Point p);
#endif
