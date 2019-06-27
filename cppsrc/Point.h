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
#define GENERAL_POINT
static const int MaxEles = 0;
struct PointEle
{
	enum {INT, FLOAT, BOOLEAN, STRING, ERROR} tag;
	union
	{
		int a;
		double b;
		bool c;
		char d[MAXSTRINGLEN];
	};
	template <class Archive>
	void serialize( Archive &ar)
	{
		ar(tag);
		switch(tag)
		{
			case INT:	ar(a);
					break;
			case FLOAT:	ar(b);
					break;
			case BOOLEAN:	ar(c);
					break;
			case STRING:	ar(d);
					break;
			case ERROR:	break;
		}
	}
	PointEle();
	PointEle(int val);
	PointEle(bool val);
	PointEle(std::string val);
	PointEle(double val);
	PointEle operator= (int val);
	operator double() const
	{
		switch(tag)
		{
			case INT:	return a;
			case FLOAT:	return b;
			case BOOLEAN:	return c;
			default:	throw std::domain_error("\n**function not defined for current Point Element**\n");
		}
	}
	PointEle operator+= (PointEle &p);
	PointEle operator-= (PointEle &p);
	PointEle operator- (PointEle &p);
	PointEle operator+ (PointEle &p);
	PointEle operator/ (int p);
};
class Point
{
//private:
public:
std::vector <PointEle> point_arr;
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
