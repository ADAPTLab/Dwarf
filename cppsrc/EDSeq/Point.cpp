#include "Point.h"
int Point::GetSize()
{
	return this->NumEles;
}

int Point::Size()
{
	return this->NumEles;
}

int Point::InitializePoint()
{
	for (int i = 0; i < NumEles; ++i)
	{
		this->point_arr[i] = 0;
	}
	return 0;
}
int Point::InitializePoint(double val)
{
	for (int i = 0; i < NumEles; ++i)
	{
		this->point_arr[i] = val;
	}
	return 0;
}

PointEle Point::GetPointEleAtIndex(int index)
{
	if (index >= 0 && index < NumEles)
		return this->point_arr[index];
	throw std::out_of_range("\n**index out of bounds**\n");
}

/*double Point::GetEleAtIndex(int index)
{
	if(index >= 0 && index < NumEles)
		return this->point_arr[index];
	throw std::out_of_range("\n**index out of bounds**\n");
}

double Point::euclidDist(Point &p)
{
	double ret = 0.0;
	for (int i = 0; i < NumEles; ++i)
	{
		ret+=pow(this->point_arr[i] - p.point_arr[i], 2);
	}
	return sqrt(ret);
}*/

PointEle Point::operator[](int index)
{
	return this->GetPointEleAtIndex(index);
}

Point &Point::operator+=(Point &p)
{
	for (int i = 0; i < NumEles; ++i)
	{
		this->point_arr[i] += p.point_arr[i];
	}
	return *this;
}

Point Point::operator+(Point &p)
{
	Point ret;
	for (int i = 0; i < NumEles; ++i)
	{
		PointEle temp = this->GetEleAtIndex(i) + p.GetEleAtIndex(i);
		ret.InsertPointEleAt(temp, i);
	}
	return ret;
}
Point &Point::operator-=(Point &p)
{
	for (int i = 0; i < NumEles; ++i)
	{
		this->point_arr[i] -= p.point_arr[i];
	}
	return *this;
}

Point Point::operator-(Point p)
{
	Point ret;
	for (int i = 0; i < NumEles; ++i)
	{
		PointEle temp = this->GetEleAtIndex(i) - p.GetEleAtIndex(i);
		ret.InsertPointEleAt(temp, i);
	}
	return ret;
}

Point Point::operator/(int p)
{
	Point ret;
	for (int i = 0; i < NumEles; ++i)
	{
		PointEle temp = this->GetEleAtIndex(i) / p;
		ret.InsertPointEleAt(temp, i);
	}
	return ret;
}

Point Point::operator/(double p)
{
	Point ret;
	for (int i = 0; i < NumEles; ++i)
	{
		PointEle temp = this->GetEleAtIndex(i) / p;
		ret.InsertPointEleAt(temp, i);
	}
	return ret;
}

Point Point::operator*(double p)
{
	Point ret;
	for (int i = 0; i < NumEles; ++i)
	{
		PointEle temp = this->GetEleAtIndex(i) * p;
		ret.InsertPointEleAt(temp, i);
	}
	return ret;
}
bool Point::operator==(Point p)
{
	for (int i = 0; i < NumEles; ++i)
	{
		if (!(this->GetEleAtIndex(i) == p.GetEleAtIndex(i)))
			return false;
	}
	return true;
	;
}
void Point::InsertPointEleAt(PointEle &pe, int index)
{
	if (index >= MaxEles)
	{
		throw std::out_of_range("\n**index out of bounds**\n");
	}
	if (0 <= index && index < NumEles)
	{
		this->point_arr[index] = pe;
	}
	else if (index == NumEles)
	{
		this->point_arr[index] = pe;
		NumEles++;
	}
	else
	{
		throw std::out_of_range("\n**index out of bounds**\n");
	}
}

void Point::InsertInteger(int val)
{
	PointEle temp(val);
	this->point_arr[NumEles] = temp;
	NumEles++;
}

void Point::InsertBoolean(bool val)
{
	throw std::domain_error("\n**function not defined for current Point Element**\n");
}

void Point::InsertString(std::string val)
{
	throw std::domain_error("\n**function not defined for current Point Element**\n");
}

void Point::InsertFloat(double val)
{
	PointEle temp(val);
	this->point_arr[NumEles] = temp;
	NumEles++;
}

void Point::SetIntegerAt(int index, int val)
{
	if (index >= MaxEles)
	{
		throw std::out_of_range("\n**index out of bounds**\n");
	}
	if (0 <= index && index < NumEles)
	{
		this->point_arr[index] = val;
	}
	else if (index == NumEles)
	{
		PointEle temp(val);
		this->point_arr[index] = temp;
		NumEles++;
	}
	else
	{
		throw std::out_of_range("\n**index out of bounds**\n");
	}
}

void Point::SetBooleanAt(int index, bool val)
{
	throw std::domain_error("\n**function not defined for current Point Element**\n");
}

void Point::SetStringAt(int index, std::string val)
{
	throw std::domain_error("\n**function not defined for current Point Element**\n");
}

void Point::SetFloatAt(int index, double val)
{
	if (index >= MaxEles)
	{
		throw std::out_of_range("\n**index out of bounds**\n");
	}
	if (0 <= index && index < NumEles)
	{
		this->point_arr[index] = val;
	}
	else if (index == NumEles)
	{
		PointEle temp(val);
		this->point_arr[index] = temp;
		NumEles++;
	}
	else
	{
		throw std::out_of_range("\n**index out of bounds**\n");
	}
}
void Point::SetDoubleAt(int index, double val)
{
	if (index >= MaxEles)
	{
		throw std::out_of_range("\n**index out of bounds**\n");
	}
	if (0 <= index && index < NumEles)
	{
		this->point_arr[index] = val;
	}
	else if (index == NumEles)
	{
		PointEle temp(val);
		this->point_arr[index] = temp;
		NumEles++;
	}
	else
	{
		throw std::out_of_range("\n**index out of bounds**\n");
	}
}
std::ostream &operator<<(std::ostream &out, Point p)
{
	int NumEles = p.GetSize();
	out << "<";
	for (int i = 0; i < NumEles; ++i)
	{
		out << p.GetPointEleAtIndex(i) << ", ";
	}
	out << ">";
	return out;
}
