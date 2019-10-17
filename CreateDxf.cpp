#include "CreateDxf.h"

CreateDxf::CreateDxf()
{
	cout <<"constructor called..." << endl;

}

CreateDxf::~CreateDxf()
{
	cout <<"Destructor called..." << endl;

}

void  CreateDxf::DxfBegin (string &str)
{
	// Creation of an output stream object in text mode. Header section of every dxf file. 
	ofstream To_Dxf (str, ios::out);

	To_Dxf << 0		  << endl;
	To_Dxf << "SECTION"  << endl;
	To_Dxf << 2		  << endl;
	To_Dxf << "ENTITIES" << endl;

	To_Dxf.close();
}

void  CreateDxf::DxfEnd (string &str)
{
	// Creation of an output stream objet in text mode. End of sequence objects of dxf file.
	ofstream To_Dxf (str, ios::app);

	To_Dxf << 0		  << endl;
	To_Dxf << "ENDSEC"   << endl;
	To_Dxf << 0		  << endl;
	To_Dxf << "EOF"<<endl;

	To_Dxf.close();
}

void CreateDxf::line ( float x1, float y1,float z1, float x2, float y2,float z2,string &str)
{	
	// Creation of an output stream objet in text mode.	
	ofstream To_Dxf (str, ios::app);

	// Draw the line
	To_Dxf << 0<< endl;
	To_Dxf << "LINE" << endl;
	To_Dxf << 8 << endl;	
	To_Dxf << 0 << endl;	// Layer number (default layer in autocad)
	To_Dxf << 10 << endl;	// XYZ is the Center point of circle
	To_Dxf << x1 << endl;	// X in UCS (User Coordinate System)coordinates
	To_Dxf << 20 << endl;
	To_Dxf << y1 << endl;	// Y in UCS (User Coordinate System)coordinates
	To_Dxf << 30 << endl;
	To_Dxf << z1 << endl;	// Z in UCS (User Coordinate System)coordinates
	To_Dxf << 11 << endl;	// XYZ is the Center point of circle
	To_Dxf << x2 << endl;	// X in UCS (User Coordinate System)coordinates
	To_Dxf << 21 << endl;
	To_Dxf << y2 << endl;	// Y in UCS (User Coordinate System)coordinates
	To_Dxf << 31 << endl;
	To_Dxf << z2 << endl;

	//close ios stream 
	To_Dxf.close();
}