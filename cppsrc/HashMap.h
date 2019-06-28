#ifndef HASHMAP_H
#define HASHMAP_H

//#define COMPRESS
#define INT_HT 1
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

#ifndef INT_HT
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
#else
template <typename T>
class HashMap
{    
public :
    std::unordered_map<uint, T> htmap;
    void serializeSend(int processor);
    HashMap<T> deserializeRecv(int recvFromProcessor);
    void broadcastMap(int numprocs);
    std::unordered_map<uint, T>& Map();
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
std::unordered_map<uint, T>& HashMap<T>::Map()
{   
    return htmap;
}

#endif
#ifndef COMPRESS
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

#else
template <typename T>
void HashMap<T>::serializeSend(int processor)
{
    int count;

    std::stringstream ss;
    cereal::BinaryOutputArchive archive { ss };
    archive( *this );
    count = ss.str().size();
    string mystr = ss.str();

       ulong c_size = mystr.size() + 100, size_taken; // Set size of compressed buffer
        char *compressed = (char*)malloc(sizeof(char) * c_size);
        // auto cbeg = MPI_Wtime();
        size_taken = Z_Compress(mystr.c_str(), mystr.size(), compressed, c_size);
        // auto cend = MPI_Wtime();
        cerr << "Uncompressed = " << count << " and compressed = " << size_taken << endl;
        MPI_Send(&size_taken, 1, MPI_UNSIGNED_LONG, processor, 0, MPI_COMM_WORLD); // Send compressed buffer size
	MPI_Send(&count, 1, MPI_UNSIGNED_LONG, processor, 0, MPI_COMM_WORLD);
        ulong times = size_taken / INT_MAX + 1;
        ulong start = 0;
        char * send_bytestream;
       	if (times > 1) {
               	send_bytestream = new char[INT_MAX]();
        } else {
                send_bytestream = new char[size_taken]();
        }
        // double send_sum = 0;
        while (times > 1) {
                for (ulong i = start; i < start + INT_MAX; i++) {
                        send_bytestream[i - start] = compressed[i];
                }
                // auto sbeg = MPI_Wtime();
                MPI_Send(send_bytestream, INT_MAX, MPI_BYTE, processor, 0,
                                MPI_COMM_WORLD);
                // auto send = MPI_Wtime();
                // send_sum += send - sbeg;
                times--;
                start += INT_MAX;
        }
	size_taken = size_taken % INT_MAX;
        if (size_taken > 0) {
                for (ulong i = start; i < start + size_taken; i++) {
                        send_bytestream[i - start] = compressed[i];
                }
                // auto sbeg = MPI_Wtime();
                MPI_Send(send_bytestream, size_taken, MPI_BYTE, processor, 0, MPI_COMM_WORLD);
                // auto send = MPI_Wtime();
                // send_sum += send - sbeg;
        }
	// cerr << "Compress Time: " << cend - cbeg << ", send time: " << send_sum <<
                //" c_size: " << size_taken << " uc_size: " << count << endl;
        delete[] send_bytestream;
        delete[] compressed;

/*
    char * send_bytestream = new char[count]();
    for (int i = 0; i < count; i++)
    {
        send_bytestream[i] = mystr[i];
    }
    MPI_Send(&count, 1, MPI_INT, processor, 0, MPI_COMM_WORLD);
    if (count > 0)
    {
        MPI_Send(send_bytestream, count, MPI_BYTE, processor, 0, MPI_COMM_WORLD);
    }
    delete[] send_bytestream;
*/
}
#endif
template <typename T>
HashMap<T>* HashMap<T>::SubList() const
{
    HashMap<T> *sublist = new HashMap();
    //sublist->this;
    return sublist;
}
#ifndef COMPRESS
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
#else
template <typename T>
HashMap<T> HashMap<T>::deserializeRecv(int recvFromProcessor)
{
	HashMap<T> localList;
       	std::stringstream ss;
        ulong count = 0, c_size = 0;
        MPI_Status status;
        MPI_Recv(&c_size, 1, MPI_UNSIGNED_LONG, recvFromProcessor, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&count, 1, MPI_UNSIGNED_LONG, recvFromProcessor, 0, MPI_COMM_WORLD, &status);
        // cerr << "Recv: " << c_size << " - " << count << endl;
        if(c_size > 0) {
                char *recv_bytestream = new char[c_size + 3]();
                memset(recv_bytestream, 0, c_size + 3);
                ulong times = c_size / INT_MAX + 1;
                ulong start = 0;
                while(times > 1) {
                        MPI_Recv(recv_bytestream + start, INT_MAX, MPI_BYTE,
                                        recvFromProcessor, 0, MPI_COMM_WORLD, &status);
                        times--;
                        start += INT_MAX;
                }
               	ulong c_size1 = c_size % INT_MAX;
                if(c_size1 > 0) {
                        MPI_Recv(recv_bytestream + start, c_size1, MPI_BYTE,
                                recvFromProcessor, 0, MPI_COMM_WORLD, &status);
                }
                char *decomp = new char[count + 3];
                Z_Decompress(recv_bytestream, c_size, decomp, count + 3);
                ss.write(&decomp[0], count);
                delete[] decomp;
                delete[] recv_bytestream;
                cereal::BinaryInputArchive iarchive { ss };
                iarchive(localList);
        }
        return localList;
}
#endif
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
