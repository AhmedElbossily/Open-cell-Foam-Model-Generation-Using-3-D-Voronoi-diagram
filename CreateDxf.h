#pragma once
#include <iostream>
#include <string.h>
#include <sstream>
#include <fstream> 
using namespace std;

class CreateDxf
{
public:
	CreateDxf();
	~CreateDxf();
	void DxfBegin (string &str);
	void DxfEnd   (string &str); 
	void line ( float x1, float y1,float z1, float x2, float y2,float z2,string &str);
};

