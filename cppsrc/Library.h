#ifndef REDUCE_H
#define REDUCE_H
#include "File.h"
#include <algorithm>
#include <vector>

template <typename T>
void* Reduce(T f(T, T), std::vector<T> &l, bool b)
{
    if (l.empty())
        return NULL;
    //std::cout<< "size of l:" << l.size() << std::endl;
    void* result;

    T Res = l[0];
    for (int i = 1; i < l.size(); i++)
    {
        Res = f(Res, l[i]);
        //if(b)   std::cout << "Point processed: " << i << " out of " << l.size() << std::endl;
    }
    if (b)
    {
        //std::cout << "before point is: " << Res << std::endl;
        //result = (T *)malloc(sizeof(T));
        result = new T();
        if (result != NULL)
            *((T *)result) = Res;
        else
            std::cout << "EXIT" << std::endl;
        //std::cout << "after point is: " << result << std::endl;
    }
    else
    {

        static int pos;
	pos = distance(l.begin(), std::find(l.begin(), l.end(), Res));
        //result = (int *)malloc(sizeof(int));
        //*((int *)result) = __pos__;
	result = &pos;
    }
    return result;
}
#ifdef GENERAL_POINT
double distanceEuclidean(Point &A, Point &B) { 
	double temp; 
	int dim; 
	double d; 
	temp = 0.0; 
	d = 0.0; 
	dim = A.Size(); 
	for (int i = 0; i < dim; i++) { 
		temp = A.GetEleAtIndex(i) - B.GetEleAtIndex(i); 
		d = d + (temp * temp); 
	} 
	d = sqrt(d); 
	return d; 
} 
#endif

/*inline double distanceEuclidean(Point &A, Point &B) {
    double temp;
    int dim;
    double d;
    temp = 0.0;
    d = 0.0;
    dim = A.Size();
    for (int i = 0; i < dim; i++) {
        temp = A.GetEleAtIndex(i) - B.GetEleAtIndex(i);
        d = d + (temp * temp);
    }
    d = sqrt(d);
    return d;
}

inline float FdistanceEuclidean(Point &A, Point &B) {
    float temp;
    int dim;
    float d;
    temp = 0.0;
    d = 0.0;
    dim = A.Size();
    for (int i = 0; i < dim; i++) {
        temp = A.GetEleAtIndex(i) - B.GetEleAtIndex(i);
        d = d + (temp * temp);
    }
    d = sqrt(d);
    return d;
}*/

double MIN(double a, double b) {
    if (a <= b) {
        return a;
    } else {
        return b;
    }
}

double MEAN(double a, double b) {
    return (a + b) / 2;
}
/*
Point MEAN(Point A, Point B) {
    return (A + B) / 2;
}

Point SUM(Point A, Point B) {
    return (A + B);
}
*/
double SUM(double A, double B) {
    return (A + B);
}

/*
template <typename T>
void* ReduceMean(std::vector<T> l, bool b)
{
    if (l.empty())
        return NULL;
    //std::cout<< "size of l:" << l.size() << std::endl;
    void* result;
    int n, i , j, listSize, pointSize;

    listSize = l.size();
    T Res = l[0];
    //std::cout << "\nBefore: " << l[0];
    double listSizef = (double)listSize;
    pointSize = Res.GetSize();
    for (i = 1; i < listSize; i++)
    {
        //if(b)   std::cout << "Point processed: " << i << " out of " << l.size() << std::endl;
        for (j = 0; j < pointSize; j++)
        {
            PointEle pe1, pe2;
            pe1 = Res.GetPointEleAtIndex(j);
            pe2 = (l[i]).GetPointEleAtIndex(j);
            if (pe1.tag == pe2.tag)
            {
                //pele.tag = pe1.tag;
                if (i < listSize )
                {
                    if (pe1.tag == PointEle::INT)
                    {
                        //pele.a = pe1.a + pe2.a;
                        Res.SetIntegerAt(j, pe1.a + pe2.a);
                    }
                    else
                    {
                        //pele.b = pe1.b + pe2.b;
                        Res.SetFloatAt(j, pe1.b + pe2.b);
                    }
                }
//                else if(i == listSize -1)
//                {
//                    if(pe1.tag == PointEle::INT)
//                    {
//                         //pele.a = (pe1.a + pe2.a) / listSize;
//                         Res.SetIntegerAt(j, (pe1.a + pe2.a) / listSize);
//                    }
//                    else
//                    {
//                        //pele.b = (pe1.b + pe2.b) / listSizef;
//                        Res.SetFloatAt(j, (pe1.b + pe2.b) / listSizef);
//                    }
//                }



                //Res.InsertPointEleAt(pele, j);
            }
        }

    }


    // std::cout << " After: " << Res;
    //std::cout << "\nlist size: " << listSize << " list size: "<< listSizef << " POint size: " << Res.Size();
    if (listSize != 0)
    {
        double listSizef = (double)listSize;
        for (i = 0; i < pointSize; i++)
        {
            PointEle pele, pe;
            pe = Res.GetPointEleAtIndex(i);
            pele.tag = pe.tag;
            if (pe.tag == PointEle::INT)
            {
                pele.a = pe.a / listSize;
                Res.SetIntegerAt(i, pele.a);
            }
            else if (pe.tag == PointEle::FLOAT)
            {
                pele.b = pe.b / listSizef;
                Res.SetFloatAt(i, pele.b);
            }


            //Res.InsertPointEleAt(pele, i);

        }
    }

    if (b)
    {
        result = new T();
        if (result != NULL)
            *((T *)result) = Res;
        else
            std::cout << "EXIT" << std::endl;

    }
    else
    {

        int pos = distance(l.begin(), std::find(l.begin(), l.end(), Res));
        result = (int *)malloc(sizeof(int));
        *((int *)result) = pos;
    }
    return result;
}*/

// function that broadcasts an object of type T
template <typename T>
void broadcastSend(T p)
{
    std::stringstream ss;
    cereal::BinaryOutputArchive archive { ss };
    archive( p );
    int count = ss.str().size();
    string mystr = ss.str();
    //char * send_bytestream = (char *) malloc(count);
    char * send_bytestream = new char[count];
    for (int i = 0; i < count; i++)
    {
        send_bytestream[i] = mystr[i];
    }
    MPI_Bcast(&count, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (count > 0)
    {
        MPI_Bcast(send_bytestream, count, MPI_BYTE, 0, MPI_COMM_WORLD);
    }
    delete[] send_bytestream;
}

// function that receives and returns an object of type T
template <typename T>
T broadcastReceive()
{
    std::stringstream ss;
    int count = 0;
    MPI_Status status;
    MPI_Bcast(&count, 1, MPI_INT, 0, MPI_COMM_WORLD);
    T p;
    if (count > 0)
    {
        char * recv_bytestream = new char[count];
        memset(recv_bytestream, 0, count);
        MPI_Bcast(recv_bytestream, count, MPI_BYTE, 0, MPI_COMM_WORLD);
        ss.write(&recv_bytestream[0], count);

        cereal::BinaryInputArchive iarchive { ss };
        iarchive( p );
        delete[] recv_bytestream;
    }
    return p;
}

// function that sends object
template <typename T>
void ObjectSend(T t, int processor)
{

    std::stringstream ss;
    cereal::BinaryOutputArchive archive { ss };
    archive( t );
    int count = ss.str().size();
    string mystr = ss.str();
    //char * send_bytestream = (char *) malloc(count);
    char * send_bytestream = new char[count];
    for (int i = 0; i < count; i++)
    {
        send_bytestream[i] = mystr[i];
    }
    MPI_Send(&count, 1, MPI_INT, processor, 0, MPI_COMM_WORLD);
    if (count > 0)
    {
        // cout << "broadcastPointSend(): Send Count = " << count << endl;
        MPI_Send(send_bytestream, count, MPI_BYTE, processor, 0, MPI_COMM_WORLD);
    }
    // free(send_bytestream);
    delete[] send_bytestream;
}

// function to receive objects
template <typename T>
T ObjectReceive(int recvFromProcessor)
{
    std::stringstream ss;
    int count = 0;
    MPI_Status status;
    MPI_Recv(&count, 1, MPI_INT, recvFromProcessor, 0, MPI_COMM_WORLD, &status);
    // cout << "broadcastPointReceive(): Count recd = " << count << ", from proc = " << 0 << endl;
    T p;
    if (count > 0)
    {
        char * recv_bytestream = new char[count]();
        memset(recv_bytestream, 0, count);
        MPI_Recv(recv_bytestream, count, MPI_BYTE, recvFromProcessor, 0, MPI_COMM_WORLD, &status);
        ss.write(&recv_bytestream[0], count);

        cereal::BinaryInputArchive iarchive { ss };
        iarchive( p );
        //cout << p;
        // free(recv_bytestream);
        delete[] recv_bytestream;
    }
    return p;
}

int int_ceil(int thread_num, int num_threads, int list_size)
{
    if (num_threads >= list_size) return thread_num;
    // Fast ceiling of an integer division x / y == (x + y - 1) / y
    int my_ceil = (list_size + num_threads - 1) / num_threads;
    int my_quot = list_size / num_threads;
    int my_rem = list_size % num_threads;
    return (thread_num <= my_rem) ? thread_num * my_ceil : my_rem * my_ceil + (thread_num - my_rem) * my_quot;
}

#endif
