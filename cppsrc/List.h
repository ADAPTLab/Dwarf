// defines the interface to the Type List<T>
#ifndef LIST_H
#define LIST_H
#include <memory>
#include <vector>
#include <climits>
#include <stdexcept>
#include <cxxabi.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include "mpi.h"

#include "Point.h"

#include "cereal/types/unordered_map.hpp"
#include "cereal/types/memory.hpp"
#include "cereal/archives/binary.hpp"
#include "cereal/types/vector.hpp"
using namespace std;
typedef unsigned long ulong;

typedef char myChar[MAXSTRINGLEN];

template <typename T>
class List
{
private :
    std::shared_ptr<std::vector<T>> Elems;
public :
	List<T>& operator+=(List<T> &other) {
		for (int i = 0; i < other.Size(); i++) {
			this->AddEle(other[i]);
		}
		return *this;
	}
	void SetElementsPtr(std::shared_ptr<std::vector<T>> ptr)
	{
		this->Elems = ptr;
	}
    void AddEle(T const& Ele); // adds an ele to the end of the list
    void SetEleAt(int index, T const& Ele); // adds an ele at the specified index in the vector, does not shift other eles
    int FindEle(T const& Ele); // returns the index of the ele if present else -1
    bool IsPresent(T const& Ele);
    void RemoveEle(T const& Ele); // removes the first occurence of the ele in the list
    void RemoveAllOcc(T const& Ele); // removes all occurrences of the ele in the list
    T& GetEleAtIndex(int index); // Returns an immutable pointer to the ele at the index if valid else NULL pointer
    int Size(); // returns the size of the list
    int Max();  //just for testing
    List<T>* Sublist(int beg_range, int end_range) const; // returns a pointer to the sublist of the org list from beg_range to end_range inclusive
    std::vector<T>& Elements();
    //std::vector<T>& RElements();
    T& operator[] (int index);
    //List<T>& operator= (const List<T>& l);
    //List<T>& operator= (List<T>& l);
    T Value (int index);
    std::string Type(int index);
    std::string Type();
    bool Free();
    // rohan
    std::shared_ptr<std::vector<T>> ElementsPtr() const;
    List<T> getLocalList();
    int getDistributionCount(int procNo, int numprocs);
    void distributeList(int numprocs);
    int getLowerLimit(int lower, int range, int processor, int numprocs);
    void distributeRangeList(int lower, int upper, int numprocs);
    void sendBackNonSerialized();
    void gatherListNonSerialized(int numprocs);
    void makeRangeList(int lower, int upper);
    void serializeSend(int processor, int lower_list_limit, int upper_list_limit);
    void nonSerializeSend(int processor, int lower_list_limit, int upper_list_limit);
    List<T> deserializeRecv(int recvFromProcessor);
    List<T> nonDeserializeRecv(int recvFromProcessor);
    void sendBack();
    void sendBackSpawn(MPI_Comm comm);
    void gatherList(int numprocs);
    T reduceList1(char* op);
    T reduceList2(char* op);
    T reduceList3(char* op);
    void broadcastList(int numprocs);
    void SetElements(List<T>& v);
    List<T>* Sort(int beg_range, int end_range);
    List<T>* Sort(int beg_range, int end_range,bool (*comparator)(T a, T b));
    template <class Archive>
    void serialize( Archive & ar )
    {
        ar( CEREAL_NVP(Elems) );
    }

    List() {
	//cout << "List creation" << endl;
        Elems = make_shared<std::vector<T>> ();
    }

    List(int num_elements)
    {
        Elems = make_shared<std::vector<T>> (num_elements);
    }

    //~List() {
        //delete Elems;
    //}
};

template<typename T>
std::vector<T>& List<T>::Elements()
{
    return *Elems;
}

// rohan
template<typename T>
std::shared_ptr<std::vector<T>> List<T>::ElementsPtr() const
{
    return Elems;
}

/*
template<typename T>
std::vector<T>& List<T>::RElements()
{
    return Elems;
}*/
template<typename T>
void List<T>::SetElements(List<T>& v)
{
    std::vector<T> vv = v.Elements();
    std::cout<< vv.size()<<"before element set\n";
    *Elems = vv;
    std::cout<< Elems->size()<<"after element set\n";
}

template <typename T>
void List<T>::AddEle(T const& Ele)
{
    try
    {
        //int i = Elems.capacity();
        Elems->push_back(Ele);
        /*if(Elems.capacity() > i)
                    std::cout << "capacity changed from " << i << " to " << Elems.capacity() << endl;
        */
    }
    catch (std::exception &e)
    {
        std::cout << std::string(e.what()) + " during List push_back" << std::endl;
    }
}

template <typename T>
void List<T>::SetEleAt(int index, T const& Ele)
{
    if (index >= 0 && index < Elems->size())
    {
        (*Elems)[index] = Ele;
    }
    else
        throw std::out_of_range("\n**index out of bounds**\n");
}

template <typename T>
int List<T>::FindEle(T const& Ele)
{
    bool flag = false;
    int i = 0;
    for (i = 0 ; i  < Elems->size(); i++)
    {
        if ((*Elems)[i] == Ele)
        {
            flag = true;
            break;
        }
    }
    if (flag)
        return i;
    else
        return -1;
}

template<typename T>
bool List<T>::IsPresent(T const& Ele)
{
    int x = FindEle(Ele);
    if (x < 0) return false;
    return true;

}

template <typename T>
T& List<T>::GetEleAtIndex(int index)
{
    if (index >= 0 && index < Elems->size())
        return ((*Elems)[index]);
    else
        throw std::out_of_range("\n**index out of bounds**\n");
}

template <typename T>
void List<T>::RemoveEle(T const& Ele)
{
    int index = FindEle(Ele);
    if (index == -1) return;
    for (int i = index; i < Elems->size()-1; i++)
    {
        (*Elems)[i] = (*Elems)[i+1];
    }
    Elems->pop_back();
}



template <typename T>
void List<T>::RemoveAllOcc(T const& Ele)
{
    int j = 0, count = 0;
    for (int i = 0; i < Elems->size(); i++)
    {
        if ((*Elems)[i] == Ele)
        {
            count++;
        }
        else
        {
            (*Elems)[j] = (*Elems)[i];
            j++;
        }
    }
    for (int i = 0; i < count; i++)
        Elems->pop_back();
//   Elems.erase(j, Elems.end());
}

template <typename T>
int List<T>::Size()
{
    return Elems->size();
}

template <typename T>
List<T>* List<T>::Sublist(int beg_range, int end_range) const
{
    if (beg_range < 0 || beg_range > end_range || end_range >= Elems->size())
        return NULL;
    List<T> *sublist = new List();
    for (int i = beg_range; i <= end_range; i++)
    {
        sublist->AddEle((*Elems)[i]);
    }
    return sublist;
}

template <typename T>
List<T> operator+ (List<T> &l1, List<T> &l2)
{
    List<T> l3;
    for (int i = 0; i < l1.Size(); i++) {
	l3.AddEle(l1.GetEleAtIndex(i));
    }
    for (int i = 0; i < l2.Size(); i++)
    {
        l3.AddEle(l2.GetEleAtIndex(i));
    }

    return l3;
    /*for (int i = 0; i < l1.Size(); i++)
    {
        l1.AddEle(l2.GetEleAtIndex(i));
    }
    return l1;*/
}

/*template <typename T>
List<T>& List<T>::operator= (const List<T>& l)
{
    Elems = l.ElementsPtr();
    return *this;
}*/

template <typename T>
List<T> operator+ (List<T> &l1, T l2)
{
    //cout << l2 << endl;
    l1.AddEle(l2);
    return l1;
}

template <typename T>
List<T> operator- (List<T> l1, List<T> l2)
{
    List<T>l3;
    for (int i = 0; i < l1.Size(); i++)
    {
        l3.AddEle(l1.GetEleAtIndex(i));
    }
    for (int i = 0; i < l2.Size(); i++)
    {
        l3.RemoveEle(l2.GetEleAtIndex(i));
    }
    return l3;
}

// removes all occurences of eles in the list
template <typename T>
List<T> operator/ (List<T> l1, List<T> l2)
{
    for (int i = 0; i < l2.Size(); i++)
    {
        l1.RemoveAllOcc(l2.GetEleAtIndex(i));
    }
    return l1;
}

// == operator for equality testing
template <typename T>
bool operator== (List<T> l1, List<T> l2)
{
    if (l1.Size() != l2.Size()) return false;
    for (int i = 0; i < l1.Size(); i++)
    {
        if ( (l1.GetEleAtIndex(i)) != (l2.GetEleAtIndex(i)))
            return false;
    }
    return true;
}

template <typename T>
std::ostream& operator << (std::ostream& out, List<T> l1)
{
    out << "[ ";
    for (auto &i : *(l1.ElementsPtr()))
        out << i << ", ";
    out << "]";
    return out;
}

template <typename T>
T& List<T>::operator[] (int index)
{
    if (index >= 0 && index < this->Elems->size())
        return (*(this->Elems))[index];
    else
        throw std::out_of_range("\n**index out of bounds**\n");
}

template <typename T>
T List<T>::Value(int index)
{
    if (index >= 0 && index < Elems->size())
        return (*Elems)[index];
    else
        throw std::out_of_range("\n**index out of bounds**\n");
}

template <typename T>
std::string List<T>::Type()
{
    int status;
//    this->
    std::string tname = typeid(T).name();
    char *demangled_name = abi::__cxa_demangle(tname.c_str(), NULL, NULL, &status);
    if(status == 0)
    {
        tname = demangled_name;
        std::free(demangled_name);
    }
    return tname;
}
template<typename T>
std::string List<T>::Type(int index)
{
    int status;
    std::string tname = typeid(T).name();
    /*char *demangled_name = abi::__cxa_demangle(tname.c_str(), NULL, NULL, &status);
    if(status == 0)
    {
        tname = demangled_name;
        std::free(demangled_name);
    }*/
    return tname;
}

template<typename T>
bool List<T>::Free()
{
    int oldSize = Size(), newSize;
    // std::cout <<"Before: " << oldSize;
    Elems->clear();
    newSize = Size();
    // std::cout <<"\tAfter: " << newSize<<endl;
    if (newSize == 0)
        return true;
    else
        return false;
}

// ROHAN
template <typename T>
int List<T>::getDistributionCount(int procNo, int numprocs)
{
    int size_b = Elems->size();
    if (size_b % (numprocs - 1) == 0)
    {
        return size_b / (numprocs - 1);
    }
    else
    {
        int remaining_b = size_b % (numprocs - 1);
        if ((remaining_b - procNo + 1) > 0)
        {
            return 1 + (size_b / (numprocs - 1));
        }
        else
        {
            return (size_b / (numprocs - 1));
        }
    }
}

#ifdef GENERAL_POINT
template <> inline
void List<Point>::serializeSend(int processor, int lower_list_limit, int upper_list_limit)
{
    int count;
    if ((upper_list_limit - lower_list_limit) < 0)
    {
        //if the size is greater than 0 only then send otherwise just send count.
        count = 0;
        MPI_Send(&count, 1, MPI_INT, processor, 0, MPI_COMM_WORLD);
        return;
    }
    List<Point> *subList = Sublist(lower_list_limit, upper_list_limit);
    count = subList->Size();
    MPI_Send(&count, 1, MPI_INT, processor, 0, MPI_COMM_WORLD);
    //cout << "serializeSend<Point>(): Sending Back Count = " << count << endl;
    Point p = (*subList)[0];
    //cout << p << endl;
    int NumEles = p.GetSize();
	//cout << "NumEles : " << NumEles << endl;
    MPI_Send(&NumEles, 1, MPI_INT, processor, 0, MPI_COMM_WORLD);
    //cout << "serializeSend<Point>(): Sending Back NumEles = " << NumEles << endl;
    for (int i = 0; i < NumEles; i++) {
        //Point p1 = p;
        PointEle pe = p[i];
        //cout << "there2" << endl;
        int tag = pe.tag;
        MPI_Send(&tag, 1, MPI_INT, processor, i, MPI_COMM_WORLD);
        //cout << "serializeSend<Point>(): Sending Back tag = " << tag << endl;
            if (tag == PointEle::INT) {
                //cout << "oh no 1" << endl;
                int *arri = new int[count];
                for (int j = 0; j < count; j++) {
                    arri[j] = (*subList)[j][i].a;
                }
                MPI_Send(arri, count, MPI_INT, processor, 0, MPI_COMM_WORLD);
                delete[] arri;
            } else if (tag == PointEle::FLOAT) {
                double *arrd = new double[count];
                for (int j = 0; j < count; j++) {
                    arrd[j] = (*subList)[j][i].b;
                }
                /*cin >> x;
                cout << "here" << endl;
                for (int j = 0; j < 10; j++) {
                    cout << "at j: " << arrd[j] << "; at opp: " << arrd[count-1-j] << endl;
                }*/
                MPI_Send(arrd, count, MPI_DOUBLE, processor, 0, MPI_COMM_WORLD);
                delete[] arrd;
                /*for (int j = 0; j < 10; j++) {
                    cout << "at j: " << arrd[j] << "; at opp: " << arrd[count-1-j] << endl;
                }
                cin >> x;
                cout << "here2" << endl;
                cin >> x;*/
            } else if (PointEle::BOOLEAN) {
      //          cout << "oh no 2" << endl;
                bool *arrb = new bool[count];
                for (int j = 0; j < count; j++) {
                    arrb[j] = (*subList)[j][i].c;
                }
                MPI_Send(arrb, count, MPI::BOOL, processor, 0, MPI_COMM_WORLD);
                delete[] arrb;
            } else if (PointEle::STRING) {
    //            cout << "oh no 3" << endl;
                myChar *arrc = new char[count][MAXSTRINGLEN];
                for (int j = 0; j < count; j++) {
                    for (int k = 0; k < MAXSTRINGLEN; k++) {
                        arrc[j][k] = (*subList)[j][i].d[k];
                    }
                }
                MPI_Send(arrc, count * MAXSTRINGLEN, MPI_BYTE, processor, 0, MPI_COMM_WORLD);
                delete[] arrc;
            }
    }
//    cout << "what?!" << endl;
    delete subList;
}
#endif // GENERAL_POINT

template <typename T>
void List<T>::serializeSend(int processor, int lower_list_limit, int upper_list_limit)
{
    ulong count;
    if ((upper_list_limit - lower_list_limit) < 0)
    {
        //if the size is greater than 0 only then send otherwise just send count.
        count = 0;
        MPI_Send(&count, 1, MPI_UNSIGNED_LONG, processor, 0, MPI_COMM_WORLD);
        return;
    }
    List<T> *subList = Sublist(lower_list_limit, upper_list_limit);
    std::stringstream ss;
    cereal::BinaryOutputArchive archive { ss };
    archive( *subList );
    delete subList;
    count = ss.str().size();
    string mystr = ss.str();
    //cout << "serializeSend(): Sending Back Count = " << count << endl;
    //cout<<ss.str()<<endl<<endl;
    // char * send_bytestream = (char *) malloc(count);
    MPI_Send(&count, 1, MPI_UNSIGNED_LONG, processor, 0, MPI_COMM_WORLD);
    ulong times = count / INT_MAX + 1;
    ulong start = 0;
    char * send_bytestream;
    if (times > 1) {
        send_bytestream = new char[INT_MAX]();
    } else {
        send_bytestream = new char[count]();
    }
    while (times > 1) {
        for (ulong i = start; i < start + INT_MAX; i++)
        {
            send_bytestream[i - start] = mystr[i];
        }
        MPI_Send(send_bytestream, INT_MAX, MPI_BYTE, processor, 0, MPI_COMM_WORLD);
        times--;
        start += INT_MAX;
    }
    count = count % INT_MAX;
    //cout << "count % INT_MAX = " << count << endl;
    if (count > 0)
    {
        for (ulong i = start; i < start + count; i++)
        {
            send_bytestream[i - start] = mystr[i];
        }
        //cout << "serializeSend(): Sending Back Count = " << count << endl;
        MPI_Send(send_bytestream, count, MPI_BYTE, processor, 0, MPI_COMM_WORLD);
    }
    // free(send_bytestream);
    delete[] send_bytestream;
}

template <typename T>
void List<T>::nonSerializeSend(int processor, int lower_list_limit, int upper_list_limit)
{
    int count;
    if ((upper_list_limit - lower_list_limit) < 0)
    {
        //if the size is greater than 0 only then send otherwise just send count.
        count = 0;
        MPI_Send(&count, 1, MPI_INT, processor, 0, MPI_COMM_WORLD);
        return;
    }
    List<T> *subList = Sublist(lower_list_limit, upper_list_limit);
    count = subList->Size();
    MPI_Send(&count, 1, MPI_INT, processor, 0, MPI_COMM_WORLD);
    cout << "nonSerializeSend(): Sending Back Count = " << count << endl;
    if (count > 0) {
        MPI_Send(&subList->ElementsPtr()->front(), count * sizeof(T), MPI_BYTE, processor, 0, MPI_COMM_WORLD);
    }
    cout << "here after condition" << endl;
    delete subList;
}

template <typename T>
void List<T>::distributeList(int numprocs)
{
    int sum = 0;
    for (int i = 1; i < numprocs; i++)
    {
        //type
        int count = getDistributionCount(i, numprocs);
        //MPI_Send(&count, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        // for(int j = 0; j<count; j++){
        //   MPI_Send(&Elems.front()+sum, sizeof(T), MPI_BYTE, i, i, MPI_COMM_WORLD);
        //   sum += 1;
        // }
        //MPI_Send(&Elems.front()+sum, count*sizeof(T), MPI_BYTE, i, 0, MPI_COMM_WORLD);
        serializeSend(i, sum, sum + count - 1);
        sum += count;
    }
}

template <typename T>
int List<T>::getLowerLimit(int lower, int range, int processor, int numprocs)
{
    if (range <= 0)
    {
        return -1;
    }
    else
    {
        if (range % (numprocs - 1) == 0)
        {
            return lower + (processor - 1) * (range / (numprocs - 1));
        }
        else
        {
            int remaining_b = range % (numprocs - 1);
            if ((remaining_b - (processor - 1)) > 0)
            {
                return lower + (processor - 1) * (range / (numprocs - 1)) + 1 * (processor - 1);
            }
            else
            {
                return lower + (processor - 1) * (range / (numprocs - 1)) + (1 * remaining_b);
            }
        }
    }
}


template <typename T>
void List<T>::distributeRangeList(int lower, int upper, int numprocs)
{
    int lower_limit_i;
    int upper_limit_prev_i;
    for (int i = 1; i < numprocs; i++)
    {
        //type
        lower_limit_i = getLowerLimit(lower, upper - lower + 1, i, numprocs);
        MPI_Send(&lower_limit_i, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        //int upper_limit_i = getRightLimit(lower, upper - lower + 1, i+1, numprocs) - 1;
        //MPI_Send(&upper_limit_i, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        if (i > 1)
        {
            upper_limit_prev_i = lower_limit_i - 1;
            MPI_Send(&upper_limit_prev_i, 1, MPI_INT, i - 1, 0, MPI_COMM_WORLD);
        }
    }
    MPI_Send(&upper, 1, MPI_INT, numprocs - 1, 0, MPI_COMM_WORLD);
}

#ifdef GENERAL_POINT
template <> inline
List<Point> List<Point>::deserializeRecv(int recvFromProcessor)
{
    int count = 0;
    MPI_Status status;
    MPI_Recv(&count, 1, MPI_INT, recvFromProcessor, 0, MPI_COMM_WORLD, &status);
    //cout << "<Point> Count recd = " << count << ", from proc = " << recvFromProcessor << endl;
    List<Point> localList(count);
    if (count > 0) {
        int NumEles = 0;
        MPI_Recv(&NumEles, 1, MPI_INT, recvFromProcessor, 0, MPI_COMM_WORLD, &status);
        //cout << "<Point> NumEles recd = " << NumEles << endl;
        for (int i = 0; i < NumEles; i++) {
            int tag = -1;
            MPI_Recv(&tag, 1, MPI_INT, recvFromProcessor, i, MPI_COMM_WORLD, &status);
            //cout << "<Point> tag recd = " << tag << endl;
                if (tag == PointEle::INT) {
                    int *arri = new int[count];
                    MPI_Recv(arri, count, MPI_INT, recvFromProcessor, 0, MPI_COMM_WORLD, &status);
                    for (int j = 0; j < count; j++) {
                        localList[j].InsertInteger(arri[j]);
                    }
                    delete[] arri;
                } else if (tag == PointEle::FLOAT) {
                    double *arrd = new double[count];
                    //cout << "here3" << endl;
                    MPI_Recv(arrd, count, MPI_DOUBLE, recvFromProcessor, 0, MPI_COMM_WORLD, &status);
                    //cout << "here4" << endl;
                    for (int j = 0; j < count; j++) {
                        localList[j].InsertFloat(arrd[j]);
                    }
                    delete[] arrd;
                    //cout << "here5" << endl;
                } else if (tag == PointEle::BOOLEAN) {
                    bool *arrb = new bool[count];
                    MPI_Recv(arrb, count, MPI::BOOL, recvFromProcessor, 0, MPI_COMM_WORLD, &status);
                    for (int j = 0; j < count; j++) {
                        localList[j].InsertBoolean(arrb[j]);
                    }
                    delete[] arrb;
                } else if (tag == PointEle::STRING) {
                    myChar *arrc = new char[count][MAXSTRINGLEN];
                    MPI_Recv(arrc, count * MAXSTRINGLEN, MPI_BYTE, recvFromProcessor, 0, MPI_COMM_WORLD, &status);
                    for (int j = 0; j < count; j++) {
                        localList[j].InsertString(arrc[j]);
                    }
                    delete[] arrc;
                }
        }
    }
    return localList;
}
#endif // GENERAL_POINT

template <typename T>
List<T> List<T>::deserializeRecv(int recvFromProcessor)
{
    List<T> localList;
    std::stringstream ss;
    ulong count = 0;
    MPI_Status status;
    MPI_Recv(&count, 1, MPI_UNSIGNED_LONG, recvFromProcessor, 0, MPI_COMM_WORLD, &status);
    //cout << "Count recd = " << count << ", from proc = " << recvFromProcessor << endl;
    if (count > 0)
    {
        // char * recv_bytestream = (char *) malloc(count);
        char * recv_bytestream = new char[count]();
        memset(recv_bytestream, 0, count);
        ulong times = count / INT_MAX + 1;
        ulong start = 0;
        while (times > 1) {
            MPI_Recv(recv_bytestream + start, INT_MAX, MPI_BYTE, recvFromProcessor, 0, MPI_COMM_WORLD, &status);
            times--;
            start += INT_MAX;
        }
        ulong count1 = count % INT_MAX;
        //cout << "count % INT_MAX = " << count1 << endl;
        if (count1 > 0)
        {
            MPI_Recv(recv_bytestream + start, count1, MPI_BYTE, recvFromProcessor, 0, MPI_COMM_WORLD, &status);
        }
        //cout << "here" << endl;
        ss.write(&recv_bytestream[0], count);
        delete[] recv_bytestream;
        //cout << "here2" << endl;

        cereal::BinaryInputArchive iarchive { ss };
        //cout << "here3" << endl;
        iarchive( localList );
        //cout << "here4" << endl;
        //DEBUG
        //for (int i = 0; i<localList.Size(); i++) {
          //T temp;
        //temp = localList.GetEleAtIndex(i);
        //cout<<"Temp.var["<<i<<"] = "<<temp.var<<", "<<temp.ab<<endl;
        //}

        // free(recv_bytestream);
    }
    return localList;
}

template <typename T>
List<T> List<T>::nonDeserializeRecv(int recvFromProcessor)
{
    List<T> localList;
    int count = 0;
    MPI_Status status;
    MPI_Recv(&count, 1, MPI_INT, recvFromProcessor, 0, MPI_COMM_WORLD, &status);
    cout << "Count recd = " << count << ", from proc = " << recvFromProcessor << endl;
    if (count > 0) {
        localList.ElementsPtr()->resize(count);
        MPI_Recv(&localList.ElementsPtr()->front(), count * sizeof(T), MPI_BYTE, recvFromProcessor, 0, MPI_COMM_WORLD, &status);
    }
    cout << "here after condition 2" << endl;
    return localList;
}

template <typename T>
List<T> List<T>::getLocalList()
{
    return deserializeRecv(0);
}

template <typename T>
void List<T>::sendBackNonSerialized()
{
    //send counts
    int count = Size();
    MPI_Send(&count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    //send data
    if (count > 0)
    {
        MPI_Send(&Elems->front(), count * sizeof(T), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
    }
}

template <typename T>
void List<T>::sendBack()
{
    //default: serialized send. Use in combination with gatherList. The unserailized version is written as sendBackNonSerialized (use in Combination with gatherListNonSerialized).
    int count = Size();
    // cout << "sendBack(): Count = " << count << endl;
    serializeSend(0, 0, count - 1);
}

template <typename T>
void List<T>::gatherListNonSerialized(int numprocs)
{
    //send counts
    int counts[numprocs - 1];
    MPI_Status status;
    int sum = 0;
    for (int i = 0; i < numprocs - 1; i++)
    {
        MPI_Recv((counts + i), 1, MPI_INT, i + 1, 0, MPI_COMM_WORLD, &status);
        //std::cout<<"Counts["<<i+1<<"] = "<<counts[i]<<std::endl;
        std::vector<T> *localElems = ElementsPtr();
        localElems->resize(sum + counts[i]);
        if (counts[i] > 0)
        {
            MPI_Recv(&Elems->front() + sum, counts[i] * sizeof(T), MPI_BYTE, i + 1, 0, MPI_COMM_WORLD, &status);
            sum += counts[i];
        }
    }
}

template <typename T>
void List<T>::gatherList(int numprocs)
{
    int sum = 0;
    this->Free();
    for (int i = 1; i < numprocs; i++)
    {
        List<T> temp = deserializeRecv(i);
        if (temp.Size() > 0)
        {
            for (auto &j : *(temp.ElementsPtr()))
            {
                this->AddEle(j);
            }
            sum += temp.Size();
        }
        temp.Free();
    }
}

template <typename T>
T List<T>::reduceList1(char* op)
{
    T res = (*Elems)[0] ;
    for (int i = 1; i < Elems->size(); i++)
    {
        if(strcmp(op,"+") == 0)
            res += (*Elems)[i];
        else if(strcmp(op,"*") == 0)
            res *= (*Elems)[i];
        else if(strcmp(op,"&&") == 0)
            res = res && (*Elems)[i];
        else if(strcmp(op,"||") == 0)
            res = res || (*Elems)[i];

    }
    return res;
}
template <typename T>
T List<T>::reduceList2(char* op)
{
    T res = (*Elems)[0] ;
    for (int i = 1; i < Elems->size(); i++)
    {
        if(strcmp(op,"+") == 0)
            res += (*Elems)[i];
        else if(strcmp(op,"*") == 0)
            res *= (*Elems)[i];
        else if(strcmp(op,"&") == 0)
            res = res & (*Elems)[i];
        else if(strcmp(op,"|") == 0)
            res = res | (*Elems)[i];
        else if(strcmp(op,"&&") == 0)
            res = res && (*Elems)[i];
        else if(strcmp(op,"||") == 0)
            res = res || (*Elems)[i];

    }
    return res;
}
template <typename T>
T List<T>::reduceList3(char *op)
{
    T res = (*Elems)[0];
    for (int i = 1; i < Elems->size(); i++) {
        if (strcmp(op, "+") == 0) {
            res = res + (*Elems)[i];
        } else if (strcmp(op, "*") == 0) {
            res = res * (*Elems)[i];
        }
    }
    return res;
}

template <typename T>
void List<T>::makeRangeList(int lower, int upper)
{
    for (int i = lower; i <= upper; i++)
    {
        this->AddEle(i);
    }
}

template <typename T>
void List<T>::broadcastList(int numprocs)
{
    int sum = 0;
    for (int i = 1; i < numprocs; i++)
    {
        //type
        int count = Size();
        serializeSend(i, 0, count - 1);
    }
}

template <typename T>
List<T>* List<T>::Sort(int beg_range, int end_range)
{
    sort(Elems->begin()+beg_range,Elems->begin()+end_range);
    return this;
}

template <typename T>
List<T>* List<T>::Sort(int beg_range, int end_range,bool (*comparator)(T a, T b))
{
    sort(Elems->begin()+beg_range,Elems->begin()+end_range,comparator);
    return this;
}

#endif
