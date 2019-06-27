#ifndef FILE_H
#define FILE_H
#include "arff_parser.h"
#include <exception>
#include <iostream>
#include <fstream>
#include <queue>
using namespace std;

//double distanceEuclidean(Point &A, Point &B);

std::string GetFileExtension(const std::string& FileName)
{
    if(FileName.find_last_of(".") != std::string::npos)
        return FileName.substr(FileName.find_last_of(".")+1);
    return "";
}
/*
bool compareVertices(Edge<Point,double>p, Edge<Point,double>q)
{
    return p.getData()>q.getData();
}
*/
class File
{
public:
    std::string filename;
    ios_base::openmode filemode;
    fstream f;

    File()
    {
        /*filename = NULL;
        filemode = NULL;
        f = NULL;*/
        filemode = std::fstream::in | std::fstream::out;
    }
    File(std::string name)
    {
        filename = name;
        this->Open(filename);
    }
    File(std::string name, ios_base::openmode mode)
    {
        filename = name;
        filemode = mode;
    }
    void Open(std::string name)
    {
        filename = name;
        filemode = std::fstream::in | std::fstream::out;
        try
        {
            //this.f = new f(filename, filemode);
            f.open(filename, filemode);
        }
        catch(exception &e)
        {
            std::cout << std::string(e.what()) << std::endl;
        }
    }
    void Open(std::string name, ios_base::openmode mode)
    {
        filename = name;
        filemode = mode;
        try
        {
            //this.f = new f(filename, filemode);
            f.open(filename, filemode);
        }
        catch(exception &e)
        {
            std::cout << std::string(e.what()) << std::endl;
        }
    }

    //Function for base case of variadic function
    template <typename T> void Write(T t)
    {
        f << t;
    }

    //recurcise function defiition for template type
    template <typename T, typename... Args> void Write(T t, Args... args)
    {
        f << t;
        Write(args...) ;
    }

    template <typename T> void Read(T *t)
    {
        f >> *t;
    }

    template <typename T, typename... Args> void Read(T *t, Args... args)
    {
        f >> *t;
        Read(args...) ;
    }
    void Close()
    {
        f.close();
    }

    List<Point> convertPointFormat(ArffData* aData)
    {
        //std::cout << "number of attributes " << aData->num_attributes()<<endl;
        //std::cout << "number of instances " << aData->num_instances()<<endl;
        int i, j;
        int32 val;
        double f;
        const char *ccc;
        std::string str;
        List<Point> pData;
	long instanceCount = aData->num_instances();
        long attribCount = aData->num_attributes();
	for(i = 0; i < instanceCount; i++)
        {
            Point p;
            for(j = 0; j < attribCount; j++)
            {

                switch(aData->get_instance(i)->get(j)->type())
                {
                case (INTEGER):
                    val = *(aData->get_instance(i)->get(j));
                    p.InsertInteger(val);
                    break;
                case (FLOAT):
                    f = *(aData->get_instance(i)->get(j));
                    p.InsertFloat(f);
                    break;
                case (NUMERIC):
                    f = *(aData->get_instance(i)->get(j));
                    p.InsertFloat(f);
                    break;

                case (STRING):
                    str= (std::string)*(aData->get_instance(i)->get(j));
                    p.InsertString(str.c_str());
                    break;
                case (NOMINAL):
                    str= (std::string)*(aData->get_instance(i)->get(j));
                    p.InsertString(str.c_str());
                    break;
                case (DATE):
                case (UNKNOWN_VAL):

                default:
                    std::cerr << "Value of this type not supported. ";
                    break;
                }
            }
            pData.AddEle(p);
            // std::cout << endl;

        }


        return pData;
    }
    void readARFFDataset(string name, List<Point>* datalist)
    {
        ArffParser parser(name);
        parser.parse(datalist);
    }
    void ReadDataset(std::string name, List<Point>* datalist)
    {
        filename = name;
        filemode = std::fstream::in;
        try
        {
            string extn = GetFileExtension(filename);
            if(extn.compare("ARFF") == 0 || extn.compare("arff") == 0)
                return (readARFFDataset(filename, datalist));
            else
                throw std::runtime_error("\nOnly ARFF File Format Supported for Dataset\n");
        }
        catch(exception &e)
        {
            std::cout << std::string(e.what()) << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
};
#endif
