// frotranRest.C
// Class to create the rest array for fortran
// Author: Sam Hewitt
// Date : 4th September 2016

#include<iostream>
#include<string.h>

namespace std
{
class fortranRest
{
    private:
    	int* Node_;	// Node Number
	int* x_;	// X restraint
	int* y_;	// Y restraint
	int* z_;	// Z restraint
	int counter_;	// Counter
    	int size_;	// Number of Restrained Nodes
    public:
	fortranRest(int size);  // constructor
	~fortranRest();		// Destructor
	void addNode(int i,int x, int y, int z);
	void editNode(int i, int x, int y, int z);
	void print(int* array);
	void printENSI(int* array);
	void printFromENSI(int* array,int arraySize);
	int getCounter();
	int getNumRestNodes();
};


//- Constructor
fortranRest::fortranRest(int size) 
{
    Node_ = new int [size];
    x_ = new int [size];
    y_ = new int [size];
    z_ = new int [size];
    counter_ = 0;
    size_ = size;
};

//-Destructor
fortranRest::~fortranRest()
{
    delete[] Node_; 
    delete[] x_;
    delete[] y_;
    delete[] z_;
};
 	
//- Add a Node to the object
void fortranRest::addNode(int i, int x, int y ,int z)
{
    Node_[counter_] = i;
    x_[counter_] = x;
    y_[counter_] = y;
    z_[counter_] = z;
    counter_++;
};

//- Edit Node (Can Only restrain in 1 Direction)
void fortranRest::editNode(int i, int x, int y, int z)
{
    if((x==2) and (y==2))
    {
	z_[i] = z;
    }
    else if((x==2) and (z==2))
    {
	y_[i] = y;       
    }
    else if((y==2) and (z==2))
    {
	x_[i] = x;
    }
    else
    {
	x_[i] = x;
	y_[i] = y;
	z_[i] = z;
    } 
};

//- If Object in REST format print directly
void fortranRest::print(int* array)
{
    memcpy(&array[size_*0],Node_, size_*sizeof(float));
    memcpy(&array[size_*1],x_,	size_*sizeof(float));
    memcpy(&array[size_*2],y_,	size_*sizeof(float));
    memcpy(&array[size_*3],z_,	size_*sizeof(float));
};

//- Print Object to ENSI Format
void fortranRest::printENSI(int* array)
{
    for(int i=0; i<(counter_); i++)
    {
	array[i*4+0] = Node_[i];
	array[i*4+1] = x_[i];
	array[i*4+2] = y_[i];
	array[i*4+3] = z_[i];
    }
};

//- Print the ENSI Format to REST array
void fortranRest::printFromENSI(int* array,int arraySize)
{
    int k=0;
    for(int j=0; j<counter_; j++)
    {
	if( (x_[j]==1) or (y_[j]==1) or (z_[j]==1) )
	{
	array[k+(arraySize*0)] = Node_[j]+1;
	array[k+(arraySize*1)] = !x_[j];
	array[k+(arraySize*2)] = !y_[j];
	array[k+(arraySize*3)] = !z_[j];
	k++;
	}
    }
};

//- Get total number of Nodes
int fortranRest::getCounter()
{	
    return counter_;
};

//- Return Total number of restrained Nodes
int fortranRest::getNumRestNodes()
{
    int numRestNodes = 0;
    for(int j=0; j<counter_; j++)
    {
	if( (x_[j]==1) | (y_[j]==1) | (z_[j]==1) )
	{
	    numRestNodes++;
	}
    }
    return numRestNodes;
};

}

