#include "Point.h"
PointEle::PointEle()
{
	// do nothing
}

PointEle::PointEle(int val)
{
	this->tag = PointEle::INT;
	this->a = val;
}

PointEle::PointEle(bool val)
{
	this->tag = PointEle::BOOLEAN;
	this->c = val;
}

PointEle::PointEle(std::string val)
{
	this->tag = PointEle::STRING;
	strncpy(this->d, val.c_str(), MAXSTRINGLEN);
	if(strlen(val.c_str()) >= MAXSTRINGLEN)
		this->d[MAXSTRINGLEN-1] = '\0';
}

PointEle::PointEle(double val)
{
	this->tag = PointEle::FLOAT;
	this->b = val;
}

PointEle PointEle::operator= (int val)
{
	switch(tag)
	{
		case INT:	this->a = val;
				break;
		case FLOAT:	this->b = val;
				break;
		case BOOLEAN:	this->c = val;
				break;
		default:	throw std::domain_error("\n**function not defined for current Point Element**\n");
	}
	return *this;
}

PointEle PointEle::operator+= (PointEle &p)
{
	if(this->tag == p.tag)
	{
		switch(tag)
		{
			case INT:	this->a+=p.a;
					break;
			case FLOAT:	this->b+=p.b;
					break;
			case BOOLEAN:	this->c+=p.c;
					break;
			default:	throw std::runtime_error("\n**Runtime error in PointEle**\n");
		}
	}
	else
	{
		// neither can be a string or bool or error now
		// otherwise at least 1 of them is a double
		// so we will typecast result to double
		if(this->tag == PointEle::STRING || p.tag == PointEle::STRING)
		{
			this->tag = ERROR;
		}
		else if(this->tag == PointEle::BOOLEAN || p.tag == PointEle::BOOLEAN)
		{
			this->tag = ERROR;
		}
		else if(this->tag == PointEle::ERROR || p.tag == PointEle::ERROR)
		{
			this->tag = ERROR;
		}
		else
		{
			double val = 0.0;
			if(this->tag == PointEle::INT)
			{
				val+=this->a;
				val+=p.b;
			}
			else
			{
				val+=this->b;
				val+=p.a;
			}
			this->tag = PointEle::FLOAT;
			this->b = val;
		}
	}
	return *this;
}

PointEle PointEle::operator- (PointEle &p)
{
	PointEle ans;
	if(this->tag != PointEle::INT && this->tag != PointEle::FLOAT)
	{
		ans.tag = PointEle::ERROR;
		return ans;
	}
	if(p.tag != PointEle::INT && p.tag != PointEle::FLOAT)
	{
		ans.tag = PointEle::ERROR;
		return ans;
	}
	if(this->tag == p.tag)
	{
		switch(tag)
		{
			case INT:	ans.tag = PointEle::INT;
					ans.a = this->a - p.a;
					break;
			case FLOAT:	ans.tag = PointEle::FLOAT;
					ans.b = this->b - p.b;
					break;
		}
	}
	else
	{
		// 1 is int the other float
		double val = 0.0;
		if(this->tag == PointEle::INT)
		{
			val = this->a - p.b;
		}
		else
		{
			val = this->b - p.a;
		}
		ans.tag = PointEle::FLOAT;
		ans.b = val;
	}
	return ans;
}
PointEle PointEle::operator-= (PointEle &p)
{
	PointEle ans;
	if(this->tag != PointEle::INT && this->tag != PointEle::FLOAT)
	{
		ans.tag = PointEle::ERROR;
		return ans;
	}
	if(p.tag != PointEle::INT && p.tag != PointEle::FLOAT)
	{
		ans.tag = PointEle::ERROR;
		return ans;
	}
	if(this->tag == p.tag)
	{
		switch(tag)
		{
			case INT:	ans.tag = PointEle::INT;
					ans.a -= p.a;
					break;
			case FLOAT:	ans.tag = PointEle::FLOAT;
					ans.b -= p.b;
					break;
		}
	}
	else
	{
		// 1 is int the other float
		double val = 0.0;
		if(this->tag == PointEle::INT)
		{
			val -= p.b;
		}
		else
		{
			val -= p.a;
		}
		ans.tag = PointEle::FLOAT;
		ans.b = val;
	}
	return ans;
}
PointEle PointEle::operator+ (PointEle &p)
{
	PointEle ans = *this;
	ans+=p;
	return ans;
}

PointEle PointEle::operator/ (int p)
{
	PointEle ans;
	if(p == 0)
	{
		ans.tag = ERROR;
		return ans;
	}
	switch(tag)
	{
		case INT:	ans.tag = PointEle::INT;
				ans.a = this->a/p;
				break;
		case FLOAT:	ans.tag = PointEle::FLOAT;
				ans.b = this->b/p;
				break;
		default:	ans.tag = ERROR;
	}
	return ans;
}
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

PointEle Point::GetPointEleAtIndex (int index)
{
	if(index >= 0 && index < NumEles)
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

PointEle Point::operator[] (int index)
{
	return this->GetPointEleAtIndex(index);
}

Point& Point::operator+= (Point &p)
{
	for (int i = 0; i < NumEles; ++i)
	{
		this->point_arr[i]+=p.point_arr[i];
	}
	return *this;
}

Point Point::operator+ (Point &p)
{
	Point ret;
	for (int i = 0; i < NumEles; ++i)
	{
		PointEle temp = this->GetEleAtIndex(i) + p.GetEleAtIndex(i);
		ret.InsertPointEleAt(temp,i);
	}
	return ret;
}
Point& Point::operator-= (Point &p)
{
	for (int i = 0; i < NumEles; ++i)
	{
		this->point_arr[i]-=p.point_arr[i];
	}
	return *this;
}

Point Point::operator- (Point p)
{
	Point ret;
	for (int i = 0; i < NumEles; ++i)
	{
		PointEle temp = this->GetEleAtIndex(i) - p.GetEleAtIndex(i);
		ret.InsertPointEleAt(temp,i);
	}
	return ret;
}

Point Point::operator/ (int p)
{
	Point ret;
	for (int i = 0; i < NumEles; ++i)
	{
		PointEle temp = this->GetEleAtIndex(i)/p;
		ret.InsertPointEleAt(temp, i);
	}
	return ret;
}

Point Point::operator/ (double p)
{
	Point ret;
	for (int i = 0; i < NumEles; ++i)
	{
		PointEle temp = this->GetEleAtIndex(i)/p;
		ret.InsertPointEleAt(temp, i);
	}
	return ret;
}

Point Point::operator* (double p)
{
	Point ret;
	for (int i = 0; i < NumEles; ++i)
	{
		PointEle temp = this->GetEleAtIndex(i)*p;
		ret.InsertPointEleAt(temp, i);
	}
	return ret;
}
bool Point::operator== (Point p)
{
	for (int i = 0; i < NumEles; ++i)
	{
		if(!(this->GetEleAtIndex(i) == p.GetEleAtIndex(i)))
			return false;
	}
	return true;;
}
void Point::InsertPointEleAt(PointEle &pe, int index)
{
	if(0 <= index && index < NumEles)
	{
		this->point_arr[index] = pe;
	}
	else if(index == NumEles)
	{
		this->point_arr.push_back(pe);
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
	this->point_arr.push_back(temp);
	NumEles++;
}

void Point::InsertBoolean(bool val)
{
	PointEle temp(val);
	this->point_arr.push_back(temp);
	NumEles++;
}

void Point::InsertString(std::string val)
{
	PointEle temp(val);
	this->point_arr.push_back(temp);
	NumEles++;
}

void Point::InsertFloat(double val)
{
	PointEle temp(val);
	this->point_arr.push_back(temp);
	NumEles++;
}

void Point::SetIntegerAt(int index, int val)
{
	if(0 <= index && index < NumEles)
	{
		this->point_arr[index] = val;
	}
	else if(index == NumEles)
	{
		PointEle temp(val);
		this->point_arr.push_back(temp);
		NumEles++;
	}
	else
	{
		throw std::out_of_range("\n**index out of bounds**\n");
	}
}

void Point::SetBooleanAt(int index, bool val)
{
	if(0 <= index && index < NumEles)
	{
		this->point_arr[index] = val;
	}
	else if(index == NumEles)
	{
		PointEle temp(val);
		this->point_arr.push_back(temp);
		NumEles++;
	}
	else
	{
		throw std::out_of_range("\n**index out of bounds**\n");
	}
}

void Point::SetStringAt(int index, std::string val)
{
	if(0 <= index && index < NumEles)
	{
		this->point_arr[index] = val;
	}
	else if(index == NumEles)
	{
		PointEle temp(val);
		this->point_arr.push_back(temp);
		NumEles++;
	}
	else
	{
		throw std::out_of_range("\n**index out of bounds**\n");
	}
}

void Point::SetFloatAt(int index, double val)
{
	if(0 <= index && index < NumEles)
	{
		this->point_arr[index] = val;
	}
	else if(index == NumEles)
	{
		PointEle temp(val);
		this->point_arr.push_back(temp);
		NumEles++;
	}
	else
	{
		throw std::out_of_range("\n**index out of bounds**\n");
	}
}

void Point::SetDoubleAt(int index, double val)
{
	if(0 <= index && index < NumEles)
	{
		this->point_arr[index] = val;
	}
	else if(index == NumEles)
	{
		PointEle temp(val);
		this->point_arr.push_back(temp);
		NumEles++;
	}
	else
	{
		throw std::out_of_range("\n**index out of bounds**\n");
	}
}

std::ostream& operator<< (std::ostream &out, Point p)
{
	int NumEles = p.GetSize();
	out<<"<";
	for (int i = 0; i < NumEles; ++i)
	{
		PointEle pl = p.GetPointEleAtIndex(i);
		switch(pl.tag)
		{
			case PointEle::INT:	out<<pl.a;
						break;
			case PointEle::FLOAT:	out<<pl.b;
						break;
			case PointEle::BOOLEAN:	out<<pl.c;
						break;
			case PointEle::STRING:	out<<pl.d;
						break;
		}
		out<<", ";
	}
	out<<">";
	return out;
}
