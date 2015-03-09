/** 
	Plot 3D spins using C++
	Shuhua Liang
	September 6th, 2009

COMPILE: g++ -O3 -g main.cpp -std=c++0x  `pkg-config --cflags eigen3`   -I/usr/include/eigen3
 **/


#include <iostream>
#include <fstream>      // std::ifstream
#include <math.h>
#include <Eigen/Dense>
using namespace Eigen;
using namespace std;


     const double PI=3.14159265358979 ;
     const int POINTS=40;   //Resolution for spin
     const double ANG1=PI*0.3, ANG2=PI*0.1;   //ANG1 is around x, ANG2 is around z
     const double ARROW_LENGTH=0.2,ARROW_SHAPE=0.12;
     const double COLUMN_LENGTH=1.0,COLUMN_SHAPE=0.032;
     const double SHIFT=0.0;


  void ellipse(double ratio, int points, double shape, double r[], double list_x[],double list_y[]);
  void draw_spin(ofstream &outputFile,double ratio, int points, double arrow_shape,  double column_shape, double r[],double rotate[]);
int main()
{
  vector<vector<vector<double> > > array3D_theta, array3D_phi;
  int ns_x=4;
  int ns_y=4;
  int ns_z=2;
  double r[3];
  double r_old[3];
  double ratio;  
  double rotate[3];


  // Set up sizes. 
  array3D_theta.resize(ns_x);
  array3D_phi.resize(ns_x);
  for (int i = 0; i < ns_x; i++) {
    array3D_theta[i].resize(ns_y);
    array3D_phi[i].resize(ns_y);
    for (int j = 0; j < ns_y; j++){
      array3D_theta[i][j].resize(ns_z);
      array3D_phi[i][j].resize(ns_z);
    }
  }

  ifstream myfile ("input.dat", ios::in);    
  ofstream outputFile;
  outputFile.open("outfile.txt");

  int a,b,c;
  for (int k = 0; k < ns_y; k++){
  for (int j = 0; j < ns_y; j++){
  for (int i = 0; i < ns_x; i++){
     myfile >> a >> b >> c>>array3D_theta[i][j][k]>>array3D_phi[i][j][k];
//     cout <<i<<" "<<j<<" "<<k<<" "<<ARROW_LENGTH*cos(array3D_theta[i][j][k])<<" "<<array3D_theta[i][j][k]<<" "<<array3D_phi[i][j][k]<<endl;
  }
  }
  }
  myfile.close();

  for (int k = 0; k < ns_z; k++){
  for (int j = 0; j < ns_y; j++){
  for (int i = 0; i < ns_x; i++){
    r_old[0]=ARROW_LENGTH*sin(array3D_theta[i][j][k])*cos(array3D_phi[i][j][k]);
    r_old[1]=ARROW_LENGTH*sin(array3D_theta[i][j][k])*sin(array3D_phi[i][j][k]);
    r_old[2]=ARROW_LENGTH*cos(array3D_theta[i][j][k]);

//     cout <<i<<" "<<j<<" "<<k<<" "<<r_old[2]<<endl;
    r[0]= cos(ANG2)*r_old[0] +cos(ANG1)*sin(ANG2)*r_old[1] +sin(ANG1)*sin(ANG2)*r_old[2];
    r[1]=-sin(ANG2)*r_old[0] +cos(ANG1)*cos(ANG2)*r_old[1] +sin(ANG1)*cos(ANG2)*r_old[2];
    r[2]=      0.0         -sin(ANG1)          *r_old[1] +cos(ANG1)          *r_old[2];

    ratio=abs(r[1]/ARROW_LENGTH);

    rotate[0]= cos(ANG2)*double(i) +cos(ANG1)*sin(ANG2)*double(j) +sin(ANG1)*sin(ANG2)*double(k);
    rotate[1]=-sin(ANG2)*double(i) +cos(ANG1)*cos(ANG2)*double(j) +sin(ANG1)*cos(ANG2)*double(k);
    rotate[2]=      0.0            -sin(ANG1)          *double(j) +cos(ANG1)          *double(k);

//draw spin:
   draw_spin(outputFile, ratio,  POINTS, ARROW_SHAPE, COLUMN_SHAPE, r, rotate);
  }
  }
  }

  outputFile.close();
}

//*******************************************************************************************
  void draw_spin(ofstream &outputFile,double ratio, int points, double arrow_shape,  double column_shape, double r[], double rotate[])
{
     double q1, q2,q3,q4,p1,p2,p3,p4,p5;
  double list_whole[POINTS][POINTS][2];
  double list_x[4*POINTS],list_y[4*POINTS];

     ellipse(ratio, points, arrow_shape, r, list_x,list_y);

  for( int l1=0; l1<points; l1++){ //axis direction
    q1=1.0/pow(double(points),2.0)*pow(double(l1),2.0);
    q2=1.0-q1;

  for( int l=0; l<points; l++){ //circle direction
    list_whole[l][l1][0]=q1*r[0]+q2*list_x[l]; 
    list_whole[l][l1][1]=q1*r[2]+q2*list_y[l]; 
  }
  }

  for( int l1=0; l1<points; l1++){ //axis direction
  for( int l=0; l<points; l++){ //circle direction
    if (r[1]>=0) {
      p4=rotate[0]-SHIFT*r[0]+list_whole[l][l1][0];
      p5=rotate[2]-SHIFT*r[2]+list_whole[l][l1][1];
      outputFile << p4<<" "<<p5 << endl;
    }
    else{
    q4=atan2(r[0], r[2]);
    p1= list_whole[l][l1][0]*cos(q4)-list_whole[l][l1][1]*sin(q4);
    p2= list_whole[l][l1][0]*sin(q4)+list_whole[l][l1][1]*cos(q4);
    p3= pow(arrow_shape, 2.0);

      if (pow(p1, 2.0)+pow(p2/ratio, 2.0)>=0.97*p3){
        p4=rotate[0]-SHIFT*r[0]+list_whole[l][l1][0];
        p5=rotate[2]-SHIFT*r[2]+list_whole[l][l1][1];
        outputFile << p4<<" "<<p5 << endl;
      }
    }

  }
  }

 for(int l=1; l<4*points; l++){ 
   outputFile << rotate[0]+(-SHIFT*r[0]+list_x[l])<<" "<<rotate[2]+(-SHIFT*r[2]+list_y[l]) << endl;
 }




//----draw column
         ellipse(ratio, points, column_shape, r, list_x,list_y);


  for( int l1=0; l1<points; l1++){ //axis direction
      q1=COLUMN_LENGTH/double(points)*double(l1);
      q2=COLUMN_LENGTH/double(points)*double(points-l1);
  for( int l=0; l<points; l++){ //circle direction
       list_whole[l][l1][0]=q1*r[0]+list_x[l]; 
       list_whole[l][l1][1]=q1*r[2]+list_y[l]; 
   }
   }

  for( int l1=0; l1<points; l1++){ //axis direction
  for( int l=0; l<points; l++){ //circle direction
    if (r[1]>=0) {
    p4=rotate[0]-SHIFT*r[0]-(q1)*r[0]+list_whole[l][l1][0];
    p5=rotate[2]-SHIFT*r[2]-(q1)*r[2]+list_whole[l][l1][1];
    outputFile << p4<<" "<<p5 << endl;
    }
    else{
    q4=atan2(r[0], r[2]);
    p1= list_whole[l][l1][0]*cos(q4)-list_whole[l][l1][1]*sin(q4);
    p2= list_whole[l][l1][0]*sin(q4)+list_whole[l][l1][1]*cos(q4);
    p3= pow(column_shape, 2.0);
      if (pow(p1, 2.0)+pow(p2/ratio, 2.0)>=0.95*p3){
        p4=rotate[0]-SHIFT*r[0]-(q1)*r[0]+list_whole[l][l1][0];
        p5=rotate[2]-SHIFT*r[2]-(q1)*r[2]+list_whole[l][l1][1];
        outputFile << p4<<" "<<p5 << endl;
      }
    }
  }
  }


     for(int l=1; l<4*points; l++){ 
       outputFile << rotate[0]+(-SHIFT*r[0]-(q1)*r[0]+list_x[l])<<" "<<rotate[2]+(-SHIFT*r[2]-(q1)*r[2]+list_y[l])<< endl;
     }
}


//*******************************************************************************************
  void ellipse(double ratio, int points, double shape, double r[], double list_x[],double list_y[])
  {
    double phase;
    double list_x_0, list_y_0;

   phase= atan(r[0]/r[2]);  /* atan2 (y, x) */

   for (int i=0;i<4*points; i++){
     list_x_0=shape*cos(2.0*PI*double(i)/double(4*points));
     list_y_0=shape*ratio*sin(2.0*PI*double(i)/double(4*points));

     list_x[i]= cos(phase)*list_x_0+sin(phase)* list_y_0  ;
     list_y[i]=-sin(phase)*list_x_0+cos(phase)* list_y_0  ;
   }
  }

