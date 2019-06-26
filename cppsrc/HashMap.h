#ifndef HASHMAP_H
#define HASHMAP_H
#include <memory>
#include <vector>
#include <stdexcept>
#include <cxxabi.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <mpi.h>

#include "cereal/types/unordered_map.hpp"
#include "cereal/types/memory.hpp"
#include "cereal/archives/binary.hpp"
#include "cereal/types/vector.hpp"
#include "cereal/types/polymorphic.hpp"
using namespace std;

template <typename T>
class HashMap
{    
public :
    std::unordered_map<std::string, List<T>> htmap;
    void serializeSend(int processor);
    HashMap<T> deserializeRecv(int recvFromProcessor);
    void broadcastMap(int numprocs);
    std::unordered_map<std::string, List<T>>& Map();
    HashMap<T>* SubList() const;
    template <class Archive>
    void serialize( Archive & ar )
    {
        ar( htmap );
    }

    int Size()
    {
        return htmap.size();
    }

};

template<typename T>
std::unordered_map<std::string, List<T> >& HashMap<T>::Map()
{   
    return htmap;
}

template <typename T>
void HashMap<T>::serializeSend(int processor)
{
    int count;
    
    std::stringstream ss;
    cereal::BinaryOutputArchive archive { ss };
    

    //HashMap<T> pointer;// = new HashMap();
    //pointer.Map() = Map();
    archive( *this );
    
    count = ss.str().size();
    //cout<<"c = "<<count<<endl;
    string mystr = ss.str();
    //cout<<ss.str()<<endl<<endl;
    // char * send_bytestream = (char *) malloc(count);
    char * send_bytestream = new char[count]();
    for (int i = 0; i < count; i++)
    {
        send_bytestream[i] = mystr[i];
    }
    MPI_Send(&count, 1, MPI_INT, processor, 0, MPI_COMM_WORLD);
    if (count > 0)
    {
        // cout << "serializeSend(): Sending Back Count = " << count << endl;
        MPI_Send(send_bytestream, count, MPI_BYTE, processor, 0, MPI_COMM_WORLD);
    }
    // free(send_bytestream);
    //delete subList;
    delete[] send_bytestream;
}

template <typename T>
HashMap<T>* HashMap<T>::SubList() const
{
    HashMap<T> *sublist = new HashMap();
    //sublist->this;
    return sublist;
}

template <typename T>
HashMap<T> HashMap<T>::deserializeRecv(int recvFromProcessor)
{
    HashMap<T> localList;
    
    std::stringstream ss;
    int count = 0;
    MPI_Status status;
    MPI_Recv(&count, 1, MPI_INT, recvFromProcessor, 0, MPI_COMM_WORLD, &status);
    // cout << "Count recd = " << count << ", from proc = " << recvFromProcessor << endl;
    if (count > 0)
    {
        // char * recv_bytestream = (char *) malloc(count);
        char * recv_bytestream = new char[count]();
        memset(recv_bytestream, 0, count);
        MPI_Recv(recv_bytestream, count, MPI_BYTE, recvFromProcessor, 0, MPI_COMM_WORLD, &status);
        ss.write(&recv_bytestream[0], count);

        cereal::BinaryInputArchive iarchive { ss };
        iarchive( localList );
        /*DEBUG
        for (int i = 0; i<localList.Size(); i++) {
          T temp;
        temp = localList.GetEleAtIndex(i);
        //cout<<"Temp.var["<<i<<"] = "<<temp.var<<", "<<temp.ab<<endl;
        }
        */
        // free(recv_bytestream);
        delete[] recv_bytestream;
    }
    return localList;
}

template <typename T>
void HashMap<T>::broadcastMap(int numprocs)
{
    int sum = 0;
    for (int i = 1; i < numprocs; i++)
    {
        //type
        int count = Size();
        serializeSend(i);
    }
}
#endif
