#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <ctime>

using namespace std;

void TriSolve(int M,
	      vector<double> a, 
	      vector<double> b, 
	      vector<double> c, 
	      vector<double> &d)
{
	for (int i = 1; i < M; i++){
		b[i] = b[i] - a[i]/b[i-1]*c[i-1];
		d[i] = d[i] - a[i]/b[i-1]*d[i-1];
	}
	d[M-1] = d[M-1]/b[M-1];
	b[M-1] = 1;
	for (int i = M-2; i > -1; i--){
		d[i] = (d[i] - c[i]/b[i+1]*d[i+1])/b[i];
		b[i] = 1;
	}
}
void TriVectorInitialize(int Mx,
						 int My,
						 double B,
						 double dx,
						 double dy,
						 vector<double> &a1,
						 vector<double> &b1,
						 vector<double> &c1,
						 vector<double> &a2,
						 vector<double> &b2,
						 vector<double> &c2,
						 vector<double> &Dx,
						 vector<double> &Dy)
{
	// Initalize values for the xdirection
	a1.resize(Mx,-B/(dx*dx));
	c1=a1;
	b1.resize(Mx,1+2*B/(dx*dx));

	// Initialize values for the y direction
	a2.resize(My,-B/(dy*dy));
	c2=a2;
	b2.resize(My,1+2*B/(dy*dy));

	// for boundary condition dT/dx = 0
	b1[0] = b1[0] + a1[0];
	b1[Mx-1] = b1[Mx-1] + c1[Mx-1];

	// for boundary condition dT/dy = 0
	b2[0] = b2[0] + a2[0];
	b2[My-1] = b2[My-1] + c2[My-1];

	//Initialize Solution vectors
	Dx.resize(Mx);
	Dy.resize(My);
}
//Define external boundary conditions 
void BC(int Mx,
		int My,
		vector<vector<double> > &T)
{
// x-z surfaces (front & back) dT/dy = 0
	for (int i = 1; i < Mx+1; i++){
			T[i][0] = T[i][1];
			T[i][My+1] = T[i][My];
	}

// y-z Surfaces (left & right) dT/dx = 0;
	for (int j = 1; j < My+1; j++){
			T[0][j] = T[1][j];
			T[Mx+1][j] = T[Mx][j];
	}
}



void LumpMass(int Mx,
			  int My,
			  double dt,
			  double dx,
			  double dy,
			  double dz,
			  double &TMold, //old mass
			  double &TMnew, //new mass
			  vector<vector<double> > &Tnew,
			  vector<vector<double> > &Told,
			  double TopCons,
			  double t)

{
	TMold = TMnew;
	double Qs;
	double Qt = 0;
	double R = 240;
	double mass = .0040;
	double Cp = 896;
	if (t == 0)
		cout << "The total resistance between the bottom copper layer and lump mass is " << R/(Mx*My) << endl;
//	double k = 204; // conduction of sink material
//	double l = dz/2 + ; //half thickness of sink
	double TopCons2 = dt/(dx*dy*dz*TopCons);
	for (int j = 1; j < My+1; j++){
		for (int i = 1; i < Mx+1; i++){
			Qs = 1/R*(TMold - Told[i][j]);
			Tnew[i][j] = Tnew[i][j] + Qs*TopCons2;
			Qt = Qt+Qs;
		}
	}
	//mas temp change from conduction
	TMnew = TMold - Qt*dt/(mass*Cp);

	//convection off lump mass
	double R2 = .375;
	double Tf = 41;
	TMnew = TMnew+dt/(R2*mass*Cp)*(Tf-TMold);
}

//Conduction between DBC and die
void ConductionBoundary(int lx,
			int ly,
			int r,
			vector<vector<double> > &TTopold,
			vector<vector<double> > &TBotold,
			vector<vector<double> > &TTopnew,
			vector<vector<double> > &TBotnew,
			double dxtop,
			double dytop,
			double dztop,
			double dxbot,
			double dybot,
			double dzbot,
			double diex,
			double diey,
			double &Qt,
			double TopCons,
			double BotCons,
			double dtR,
			double dt)
{
	Qt = 0;
	double dTs;
	double dTt=0;
	int lbx;
	int ubx;
	int lby;
	int uby;
	double TopCons1 = TopCons*dxtop*dytop*dztop;
	double BotCons1 = BotCons*dxbot*dybot*dzbot;
	double TopCons2 = dtR/TopCons1;
	double BotCons2 = dtR/BotCons1;

	//conduction through bottom
	for (int n = diey; n < diey+ly; n++){
		lby = (n-diey)*r+1;
		uby = (n-diey)*r+r+1;
		for (int m = diex; m < diex+lx; m++){
			lbx = (m-diex)*r+1;
			ubx = (m-diex)*r+r+1;
			for (int j = lby; j < uby; j++){
				for (int i = lbx; i < ubx; i++){
					dTs = TBotold[m][n] - TTopold[i][j];
					TTopnew[i][j] = TTopnew[i][j] + dTs*TopCons2;
					dTt = dTt+dTs;
				}
			}
			TBotnew[m][n] = TBotnew[m][n] - dTt*BotCons2;
			Qt = Qt + dTt;	//Running value of sum of temperature differences between surfaces
			dTt = 0;	//reset for next DBC element
		}
	}
	//Total power between DBC and die
	Qt = Qt*dtR/dt;;	//Convert sum of temperature differences into a total power between surfaces
}

void InternalHeatGeneration(double qt,
			int Mx,
			int My,
			double dx,
			double dy,
			double dz,
			double dt,
			double t,
			double Cons,
			vector<vector<double> > &Tnew)
{
	double q = qt/(Mx*My*dx*dy*dz*.7225); //.85*.85 = 72.25% of crosssectional area is active
	int lbx = 3*Mx / 40 + 1;
	int	ubx = Mx - lbx + 2;
	int lby = 3*My / 40 + 1;
	int uby = My - lby + 2;
	for (int i = lbx; i < ubx; i++){
		for(int j = lby; j < uby; j++){
			Tnew[i][j] = Tnew[i][j] + q*dt/(Cons);
		}
	}
}


void ADI2d(int Mx,
		   int My,
		   vector<double> &a, //lowercase for x direction
		   vector<double> &b,
		   vector<double> &c,
		   vector<double> &A, //uppercase for y direction
		   vector<double> &B,
		   vector<double> &C,
		   vector<double> &Dx, //solution vector for x
		   vector<double> &Dy, //solution vector for y
		   double betax, //beta divided by dx^2
		   double betay, //beta divided by dy^2
		   vector<vector<double> > &Tnew,
		   vector<vector<double> > &Told)
{
	Told = Tnew; //update old vector with the newest
	//x direction solving
	for (int j = 1; j < My+1; j++){
		for (int i = 1; i < Mx+1; i++){
			Dx[i-1] = betay*(Told[i][j+1]-2*Told[i][j]+Told[i][j-1])+Told[i][j];
		}
		TriSolve(Mx,a,b,c,Dx);
		for (int i = 0; i < Mx; i++)
			Tnew[i+1][j] = Dx[i];
		Tnew[0][j] = Tnew[1][j];	//apply BCs to intermediate step
		Tnew[Mx+1][j] = Tnew[Mx][j];
	}
	//y direction
	for (int i = 1; i < Mx+1; i++){
		for (int j = 1; j < My+1; j++){
			Dy[j-1] = betax*(Tnew[i+1][j]-2*Tnew[i][j]+Tnew[i-1][j])+Tnew[i][j];
		}
		TriSolve(My,A,B,C,Dy);
		for (int j = 0; j < My; j++)
			Tnew[i][j+1] = Dy[j];
	}
}
void TempAverages(int Mx,
		  int My,
		  double &T1,
		  double &T2,
		  vector<vector<double> > &TNew,
		  int diex1,
		  int diey1)
{
	T1 = 0;
	T2 = 0;
	//Active region average
	for (int j = 4; j < My-2; j++){
		for (int i = 4; i < Mx-2; i++){
			T1 = T1 + TNew[i][j];
		}
	}
	T1 = T1/((Mx-6)*(My-6));

	//Hottest values
	for (int j = My/2; j < My/2+2; j++){
		for (int i = Mx/2; i < Mx/2+2; i++){
			T2 = T2 + TNew[i][j];
		}
	}
	T2 = T2/4;
}
void fileoutput(int Mx1,
		int My1,
		int Mx2,
		int My2,
		int Mx3,
		int My3,
		int Mx4,
		int My4,
		int Mx5,
		int My5,
		int Mx6,
		int My6,
		vector<double> x1,
		vector<double> y1,
		double z1,
		vector<double> x2,
		vector<double> y2,
		double z2,
		vector<double> x3,
		vector<double> y3,
		double z3,
		vector<double> x4,
		vector<double> y4,
		double z4,
		vector<double> x5,
		vector<double> y5,
		double z5,
		vector<double> x6,
		vector<double> y6,
		double z6,
		vector<vector<double> > &T1,
		vector<vector<double> > &T2,
		vector<vector<double> > &T3,
		vector<vector<double> > &T4,
		vector<vector<double> > &T5,
		vector<vector<double> > &T6,
		double t)
{
	string stringbase;
	string string1 = "t = ";
	string string2 = ".txt";

	//Convert double class time to string
	string time;
	ostringstream conversion;
	conversion << t;
	time = conversion.str();

	//Construct filename string
	stringbase.append(string1);
	stringbase.append(time);
	stringbase.append(string2);

	ofstream Temps;
	Temps << scientific;
	Temps.open(stringbase.c_str(), ofstream::out);
	Temps << setw(8) << "x1" << "	" << setw(8) << "y1" << "	" << setw(8) << "z" << "	" 
			<< setw(12) << "Temp (C)" << endl;
	
	//for (int j = 1; j < My1+1; j++){
	//	for (int i = 1; i < Mx1+1; i++){
	//		Temps.precision(3);
	//		Temps << setw(8) << x1[i-1] << "	" << setw(8) << y1[j-1] << "	" << setw(8) << z1 << "	";
	//		Temps.precision(6);
	//		Temps << setw(16) << T1[i][j] << endl;
	//	}
	//}
	for (int j = 1; j < My2+1; j++){
		for (int i = 1; i < Mx2+1; i++){
			Temps.precision(3);
			Temps << setw(8) << x2[i-1] << "	" << setw(8) << y2[j-1] << "	" << setw(8) << z2 << "	";
			Temps.precision(6);
			Temps << setw(16) << T2[i][j] << endl;
		}
	}
	for (int j = 1; j < My3+1; j++){
		for (int i = 1; i < Mx3+1; i++){
			Temps.precision(3);
			Temps << setw(8) << x3[i-1] << "	" << setw(8) << y3[j-1] << "	" << setw(8) << z3 << "	";
			Temps.precision(6);
			Temps << setw(16) << T3[i][j] << endl;
		}
	}
	for (int j = 1; j < My4+1; j++){
		for (int i = 1; i < Mx4+1; i++){
			Temps.precision(3);
			Temps << setw(8) << x4[i-1] << "	" << setw(8) << y4[j-1] << "	" << setw(8) << z4 << "	";
			Temps.precision(6);
			Temps << setw(16) << T4[i][j] << endl;
		}
	}
	for (int j = 1; j < My5+1; j++){
		for (int i = 1; i < Mx5+1; i++){
			Temps.precision(3);
			Temps << setw(8) << x5[i-1] << "	" << setw(8) << y5[j-1] << "	" << setw(8) << z5 << "	";
			Temps.precision(6);
			Temps << setw(16) << T5[i][j] << endl;
		}
	}
	for (int j = 1; j < My6+1; j++){
		for (int i = 1; i < Mx6+1; i++){
			Temps.precision(3);
			Temps << setw(8) << x6[i-1] << "	" << setw(8) << y6[j-1] << "	" << setw(8) << z6 << "	";
			Temps.precision(6);
			Temps << setw(16) << T6[i][j] << endl;
		}
	}
	Temps.close();
}
		
	//Main function call
int main()
{
//Define variables, units in standard metric (m,kg,J,s)

	double tstart = clock();
	double runtime;

	//silicon properties
	const double SiDens = 2400;
	const double SiCond = 100;
	const double SiCp = 790;
	const double SiCons = SiDens*SiCp;
	const double alpha1 = SiCond/(SiDens*SiCp); //Thermal diffusivity Silicon (die)
	cout << "Silicon Alpha	" << alpha1 << endl;

	//copper properties
	const double CuDens = 8920;
	const double CuCond = 380;
	const double CuCp = 383;
	const double CuCons = CuDens*CuCp;
	const double alpha2 = CuCond/(CuDens*CuCp); //Thermal diffusivitiy Copper (dbc)
	cout << "Copper Alpha	" << alpha2 << endl;

	//Alumina properties (using 94.5%)
	const double AlDens = 3690;
	const double AlCond = 18;
	const double AlCp = 880;
	const double AlCons = AlDens*AlCp;
	const double alpha3 = AlCond/(AlDens*AlCp); //Thermal diffusivity alumina (dbc)
	cout << "Alumina Alpha	" << alpha3 << endl;

	// High Temp Solder ()
	const double SoDens = 7420;
	const double SoCond = 81.7;
	const double SoCp = 228;
	const double SoCons = SoDens*SoCp;
	const double alpha4 = SoCond/(SoDens*SoCp);
	cout << "Solder Alpha	" << alpha4 << endl;

	//mesh expansion ratios
	double r1 = 1; //between active and inactive die layers
	double r2 = 2;	//between inactive die and solder layers
	double r3 = 2;	//between solder layer and copper pad
	double r4 = 2;	//between copper pad and alumina
	double r5 = 1;	//between alumin and bottom copper layer

	//Die information
	//Active region, top layer (5 microns deep)
	const double xdim = .00625;
	const double ydim = .0055;
	const double zdim = .00013335;
	double Mx = 40;
	double My = 40;
	int Mx1 = Mx;
	int My1 = My;
	double dx1 = xdim/Mx; //Mesh spacing 
	double dy1 = ydim/My;
	double dz1 = 0e-6;
	vector<double> a1;
	vector<double> b1;
	vector<double> c1;
	vector<double> A1;
	vector<double> B1;
	vector<double> C1;
	vector<double> Dx1;
	vector<double> Dy1;
	double DieArea = Mx*My*dx1*dy1;
	cout << "Total area of the die is " << DieArea << endl;

	//Nonactive die region (the remainder)
	double dx2 = dx1*r1;
	double dy2 = dy1*r1;
	double dz2 = zdim - dz1;
	int Mx2 = Mx1/r1;
	int My2 = My1/r1;
	vector<double> a2;
	vector<double> b2;
	vector<double> c2;
	vector<double> A2;
	vector<double> B2;
	vector<double> C2;
	vector<double> Dx2;
	vector<double> Dy2;

	//Solder layer (first layer of mesh expansion)
	double dx3 = dx2*r2;
	double dy3 = dy2*r2;
	double dz3 = 75e-6;
	int Mx3 = Mx1/(r1*r2);
	int My3 = My1/(r1*r2);
	vector<double> a3;
	vector<double> b3;
	vector<double> c3;
	vector<double> A3;
	vector<double> B3;
	vector<double> C3;
	vector<double> Dx3;
	vector<double> Dy3;

	//DBC information
	//Top layer	(reduced size)
	double dx4 = dx3*r3;	//at 4x the size, total DBC1 dimensions are 3x those of die
	double dy4 = dy3*r3;
	double dz4 = 300e-6;	//Copper (DBC) thickness
	int Mx4 = Mx1/(r1*r2*r3)+2;	//1 element of overlap in each direction
	int My4 = My1/(r1*r2*r3)*2.6;	//1 extra copper pad element as buffer past die ylength
	vector<double> a4;
	vector<double> b4;
	vector<double> c4;
	vector<double> A4;
	vector<double> B4;
	vector<double> C4;
	vector<double> Dx4;
	vector<double> Dy4;
	double DBC4Area = Mx4*My4*dx4*dy4;
	cout << "Total area of the top DBC layer is " << DBC4Area << endl;

	//Middle layer (full size)
	double dx5 = dx4*r4;
	double dy5 = dy4*r4;
	double dz5 = 380e-6;	//Alumina layer thickness total
	int Mx5 = Mx1/(r1*r2*r3*r4)*5;	//5x length of die
	int My5 = Mx1/(r1*r2*r3*r4)*7;	//7x height of die 
	vector<double> a5;
	vector<double> b5;
	vector<double> c5;
	vector<double> A5;
	vector<double> B5;
	vector<double> C5;
	vector<double> Dx5;
	vector<double> Dy5;
	double DBC5Area = Mx5*My5*dx5*dy5;
	cout << "Total area of the middle DBC layer is " << DBC5Area << endl;

	//Bottom layer (full size)
	double dx6 = dx5*r5;
	double dy6 = dy5*r5;
	double dz6 = 200e-6; //bottom copper layer thickness
	int Mx6 = Mx5/r5;
	int My6 = My5/r5;
	vector<double> a6;
	vector<double> b6;
	vector<double> c6;
	vector<double> A6;
	vector<double> B6;
	vector<double> C6;
	vector<double> Dx6;
	vector<double> Dy6;
	double DBC6Area = Mx6*My6*dx6*dy6;
	cout << "Total area of the bottom DBC layer is " << DBC6Area << endl;

//Temperature/Spatial domain structure
	//Initialize Temperature arrays
	vector<vector<double> >  TDie1new;	//Die1 is active layer of silicon
	vector<vector<double> >  TDie1old;
	vector<vector<double> >  TDie2new;	//Die2 is inactive layer of silicon
	vector<vector<double> >  TDie2old;
	vector<vector<double> >  TSold3new;	//Sold3 is solder layer
	vector<vector<double> >  TSold3old;
	vector<vector<double> >  TDBC4new;	//DBC4 is copper pad
	vector<vector<double> >  TDBC4old;
	vector<vector<double> >  TDBC5new;	//DBC5 is DBC layer of alumina
	vector<vector<double> >  TDBC5old;
	vector<vector<double> >  TDBC6new;	//DBC6 is DBC copper bottom layer 
	vector<vector<double> >  TDBC6old;

//Sizing & Initializing 2D temperature arrays (IC = 42.17C)
	//active region temperatures
	TDie1new.resize(Mx1+2);
	for(int i = 0; i < Mx1+2; i++)
		TDie1new[i].resize(My1+2,40);
	TDie1old = TDie1new;

	//inactive region temperatures
	TDie2new.resize(Mx2+2);
	for (int i = 0; i < Mx2+2; i++)
		TDie2new[i].resize(My2+2,40);
	TDie2old = TDie2new;

	//solder temperatures
	TSold3new.resize(Mx3+2);
	for (int i = 0; i < Mx3+2; i++)
		TSold3new[i].resize(My3+2,40);
	TSold3old = TSold3new;

	//copper pad temperatures
	TDBC4new.resize(Mx4+2);
	for  (int i = 0; i < Mx4+2; i++)
		TDBC4new[i].resize(My4+2,40);
	TDBC4old = TDBC4new;

	//Alumina temperatures
	TDBC5new.resize(Mx5+2);
	for  (int i = 0; i < Mx5+2; i++)
		TDBC5new[i].resize(My5+2,40);
	TDBC5old = TDBC5new;
	//Bottom copper layer temperatures
	TDBC6new.resize(Mx6+2);
	for (int i = 0; i < Mx6+2; i++)
		TDBC6new[i].resize(My6+2,40);
	TDBC6old = TDBC6new;

	double TActive;
	double TMax;

	double TMold = 41;
	double TMnew = TMold;

	//Initialize coordinate vectors
	vector<double> x1; vector<double> y1;
	vector<double> x2; vector<double> y2;
	vector<double> x3; vector<double> y3;
	vector<double> x4; vector<double> y4;
	vector<double> x5; vector<double> y5;
	vector<double> x6; vector<double> y6;

	//Initializing Coordinate lengths
	x1.resize(Mx1,0); y1.resize(My1,0);	//die1
	x2.resize(Mx2,0); y2.resize(My2,0);	//die2
	x3.resize(Mx3,0); y3.resize(My3,0);	//solder3
	x4.resize(Mx4,0); y4.resize(My4,0);	//dbc4
	x5.resize(Mx5,0); y5.resize(My5,0);	//dbc5
	x6.resize(Mx6,0); y6.resize(My6,0);	//dbc6

//DBC1 element index for die location
	double diex1 = 1; double diey1 = 1;	//location of die1 on die2
	double diex2 = 1; double diey2 = 1; 	//location of die2 on solder3
	double diex3 = 2; double diey3 = My4-My3/(r3); //location of solder3 on dbc4 (copper pad), offset by 1 length from bottom

	double diex4 = Mx1/(r1*r2*r3*r4)*1.8; //1.8 die lengths from left edge
	double diey4 = My5-My1/(r1*r2*r3*r3*r4)*1.6-My4/r4; //1.6 die lengths from top

	double diex5 = 1; //location of alumina on copper pad
	double diey5 = 1;

//Spatial Coordinates, origin is bottom left of DBC (looking down)
	for (int i = 0; i < Mx1; i++)
		x1[i] = (0.5+i)*dx1+(diex1-1)*dx2+(diex2-1)*dx3+(diex3-1)*dx4+(diex4-1)*dx5+(diex5-1)*dx6;
	for (int i = 0; i < My1; i++)
		y1[i] = (0.5+i)*dy1+(diey1-1)*dy2+(diey2-1)*dy3+(diey3-1)*dy4+(diey4-1)*dy5+(diey5-1)*dy6;
	for (int i = 0; i < Mx2; i++)
		x2[i] = (0.5+i)*dx2+(diex2-1)*dx3+(diex3-1)*dx4+(diex4-1)*dx5+(diex5-1)*dx6;
	for (int i = 0; i < My2; i++)
		y2[i] = (0.5+i)*dy2+(diey2-1)*dy3+(diey3-1)*dy4+(diey4-1)*dy5+(diey5-1)*dy6;
	for (int i = 0; i < Mx3; i++)
		x3[i] = (0.5+i)*dx3+(diex3-1)*dx4+(diex4-1)*dx5+(diex5-1)*dx6;
	for (int i = 0; i < My3; i++)
		y3[i] = (0.5+i)*dy3+(diey3-1)*dy4+(diey4-1)*dy5+(diey5-1)*dy6;
	for (int i = 0; i < Mx4; i++)
		x4[i] = (0.5+i)*dx4+(diex4-1)*dx5+(diex5-1)*dx6;
	for (int i = 0; i < My4; i++)
		y4[i] = (0.5+i)*dy4+(diey4-1)*dy5+(diey5-1)*dy6;
	for (int i = 0; i < Mx5; i++)
		x5[i] = (0.5+i)*dx5+(diex5-1)*dx6;
	for (int i = 0; i < My5; i++)
		y5[i] = (0.5+i)*dy5+(diey5-1)*dy6;
	for (int i = 0; i < Mx6; i++)
		x6[i] = (0.5+i)*dx6;
	for (int i = 0; i < My6; i++)
		y6[i] = (0.5+i)*dy6;

	double z1 = dz6+dz5+dz4+dz3+dz2+dz1*.5;
	double z2 = dz6+dz5+dz4+dz3+dz2*.5;
	double z3 = dz6+dz5+dz4+dz3*.5;
	double z4 = dz6+dz5+dz4*.5;
	double z5 = dz6+dz5*.5;
	double z6 = dz6*.5;

	//Time vs. Temperature logging
	double t = 0; //initialize time at 0
	double dt = 75e-6;
	cout << "The chosen timestep is " << dt*1e6 << " microseconds." << endl;

//Dependent constants
	double BetaSi = alpha1*dt/2;	//beta for silicon
	double BetaCu = alpha2*dt/2;	//beta for copper
	double BetaAl = alpha3*dt/2;	//beta for alumina
	double BetaSo = alpha4*dt/2;	//beta for solder
	double Beta1x = BetaSi/(dx1*dx1);	//active die layer
	double Beta1y = BetaSi/(dy1*dy1);
	double Beta2x = BetaSi/(dx2*dx2);	//inactive die layer
	double Beta2y = BetaSi/(dy2*dy2);
	double Beta3x = BetaSo/(dx3*dx3);	//solder layer
	double Beta3y = BetaSo/(dy3*dy3);
	double Beta4x = BetaCu/(dx4*dx4);	//copper layer, top (pad)
	double Beta4y = BetaCu/(dy4*dy4);	
	double Beta5x = BetaAl/(dx5*dx5);	//alumina layer
	double Beta5y = BetaAl/(dy5*dy5);
	double Beta6x = BetaCu/(dx6*dx6);	//copper layer, bottom
	double Beta6y = BetaCu/(dy6*dy6);

//Initialize tri-vectors for each layer
	TriVectorInitialize(Mx1,My1,BetaSi,dx1,dy1,a1,b1,c1,A1,B1,C1,Dx1,Dy1);	//Active die layer
	TriVectorInitialize(Mx2,My2,BetaSi,dx2,dy2,a2,b2,c2,A2,B2,C2,Dx2,Dy2);	//Inactive die layer
	TriVectorInitialize(Mx3,My3,BetaSo,dx3,dy3,a3,b3,c3,A3,B3,C3,Dx3,Dy3);	//Solder layer
	TriVectorInitialize(Mx4,My4,BetaCu,dx4,dy4,a4,b4,c4,A4,B4,C4,Dx4,Dy4);	//Copper pad layer
	TriVectorInitialize(Mx5,My5,BetaAl,dx5,dy5,a5,b5,c5,A5,B5,C5,Dx5,Dy5);	//Alumina layer
	TriVectorInitialize(Mx6,My6,BetaCu,dx6,dy6,a6,b6,c6,A6,B6,C6,Dx6,Dy6);	//Bottom Copper layer
	
//Interface information
	//Between Active and Inactive die layers
	double R1 = dz1/(2*SiCond*dx1*dy1) + dz2/(2*SiCond*dx1*dy1);
	double Qt1 = 0;	//net power between solder layers
	
	//Between Inactive die and Solder
	double R2 = dz2/(2*SiCond*dx2*dy2)+dz3/(2*SoCond*dx2*dy2)+50;
	double Qt2 = 0;
	cout << "Total resistance between Inactive Die and Solder layers	" << R2/(Mx2*My2) << endl;

	//Between Solder and Copper pad
	double R3 = dz3/(2*SoCond*dx3*dy3) + dz4/(2*CuCond*dx3*dy3)+25;
	double Qt3 = 0;
	cout << "Total resistance between Solder layer and Copper pad	" << R3/(Mx3*My3) << endl;
	
	//Between Copper pad and Alumina
	double R4 = dz4/(2*CuCond*dx4*dy4) + dz5/(2*AlCond*dx4*dy4)+43;
	double Qt4 = 0;
	cout << "Total resistance between Copper pad and Alumina layers	" << R4/(Mx4*My4) << endl;

	//Between Alumina and bottom Copper layer
	double R5 = dz5/(2*AlCond*dx5*dy5) + dz6/(2*CuCond*dx5*dy5);
	double Qt5 = 0;
	cout << "Total resistance between Alumina layer and bottom Copper layer	" << R5/(Mx5*My5) << endl;

	//Times for full temperature array output
	double outtime[19] = {0, .001, 0.01, 0.04, 0.2, 1, 5, 10, 14.9, 
		15, 15.001, 15.01, 15.04, 15.2, 16, 20, 25, 29.9, 100}; //100 to avoid bad indexing

	TempAverages(Mx1,My1,TActive,TMax,TDie1new,diex1,diey1);
	ofstream Masstemp;
	Masstemp.open("TimevTemps.txt");
	Masstemp << setw(8) << "Time (s)"<< "	" << setw(11) << "Mass Temp (C)" << "	" << setw(14) << "Die Top Temp (C)" << "	" << setw(10) << "Contact Power" << endl;
	Masstemp.precision(8);
	Masstemp << scientific <<  setw(15) << t;
	Masstemp.precision(3);
	Masstemp << "	" << fixed << setw(8) << TMnew << "	" << setw(8) << TActive << "	" << setw(8) << -Qt5 << "	" << setw(8) << TMax << endl;

	//Modulus values for die logging
	int dielog;
	double logtimes [7] = {0, .1, 1, 15, 15.01, 16, 100};
	int logsteps [7] = {10, 50, 500, 10, 50, 500, 1};

	//Dependent Constants set 2
	double dtR1 = dt/R1;
	double dtR2 = dt/R2;
	double dtR3 = dt/R3;
	double dtR4 = dt/R4;
	double dtR5 = dt/R5;

	int iteration = 0;
	int outcount = 0;
	int dtcount = 0;
	int logcount = 0;
	//Begin solution loop 

	while (t<=30){

		//BC(Mx1,My1,TDie1new);
		BC(Mx2,My2,TDie2new);
		BC(Mx3,My3,TSold3new);
		BC(Mx4,My4,TDBC4new);
		BC(Mx5,My5,TDBC5new);
		BC(Mx6,My6,TDBC6new);

		//ADI2d(Mx1,My1,a1,b1,c1,A1,B1,C1,Dx1,Dy1,Beta1x,Beta1y,TDie1new,TDie1old);
		ADI2d(Mx2,My2,a2,b2,c2,A2,B2,C2,Dx2,Dy2,Beta2x,Beta2y,TDie2new,TDie2old);
		ADI2d(Mx3,My3,a3,b3,c3,A3,B3,C3,Dx3,Dy3,Beta3x,Beta3y,TSold3new,TSold3old);
		ADI2d(Mx4,My4,a4,b4,c4,A4,B4,C4,Dx4,Dy4,Beta4x,Beta4y,TDBC4new,TDBC4old);
		ADI2d(Mx5,My5,a5,b5,c5,A5,B5,C5,Dx5,Dy5,Beta5x,Beta5y,TDBC5new,TDBC5old);
		ADI2d(Mx6,My6,a6,b6,c6,A6,B6,C6,Dx6,Dy6,Beta6x,Beta6y,TDBC6new,TDBC6old);
		
		//ConductionBoundary(Mx1/r1,My1/r1,r1,TDie1old,TDie2old,TDie1new,TDie2new,dx1,dy1,dz1,dx2,dy2,dz2,diex1,diey1,Qt1,SiCons,SiCons,dtR1,dt);
		ConductionBoundary(Mx2/r2,My2/r2,r2,TDie2old,TSold3old,TDie2new,TSold3new,dx2,dy2,dz2,dx3,dy3,dz3,diex2,diey2,Qt2,SiCons,SoCons,dtR2,dt);
		ConductionBoundary(Mx3/r3,My3/r3,r3,TSold3old,TDBC4old,TSold3new,TDBC4new,dx3,dy3,dz3,dx4,dy4,dz4,diex3,diey3,Qt3,SoCons,CuCons,dtR3,dt);
		ConductionBoundary(Mx4/r4,My4/r4,r4,TDBC4old,TDBC5old,TDBC4new,TDBC5new,dx4,dy4,dz4,dx5,dy5,dz5,diex4,diey4,Qt4,CuCons,AlCons,dtR4,dt);
		ConductionBoundary(Mx5/r5,My5/r5,r5,TDBC5old,TDBC6old,TDBC5new,TDBC6new,dx5,dy5,dz5,dx6,dy6,dz6,diex5,diey5,Qt5,AlCons,CuCons,dtR5,dt);

		//lump mass on bottom surface, with convection to fluid sink
		LumpMass(Mx6,My6,dt,dx6,dy6,dz6,TMold,TMnew,TDBC6new,TDBC6old,CuCons,t);
		if (t < 15) //Power on for 15s
			InternalHeatGeneration(42,Mx2,My2,dx2,dy2,dz2,dt,t,SiCons,TDie2new);

		t = t+dt;
		iteration++;

		if (t >= logtimes[logcount])
		{
			dielog = logsteps[logcount];
			logcount++;
		}

		if (iteration % dielog == 0)
		{	
			TempAverages(Mx1,My1,TActive,TMax,TDie2new,diex1,diey1);
			Masstemp.precision(8);
			Masstemp << scientific <<  setw(15) << t;
			Masstemp.precision(3);
			Masstemp << "	" << fixed << setw(8) << TMold << "	" << setw(8) << TActive << "	" << setw(8) << -Qt5 << "	" << setw(8) << TMax <<  endl;
		}


		if (iteration%10000 == 0){
			runtime = (clock()-tstart)/CLOCKS_PER_SEC;
			cout << "The current time is " << t << " after a runtime of " << runtime << " seconds." << endl;
		}
		if (t > outtime[outcount]){
			fileoutput(Mx1,My1,Mx2,My2,Mx3,My3,Mx4,My4,Mx5,My5,Mx6,My6,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,
				TDie1new,TDie2new,TSold3new,TDBC4new,TDBC5new,TDBC6new,outtime[outcount]);
			outcount++;
		}
	}

//Output results to command window
		runtime = (clock()-tstart)/CLOCKS_PER_SEC;
		cout << "A final simulated time of: " << endl << t << "s was attained after "<< iteration << " timesteps." << endl;
		cout << "The total time spent running is calculated as " << runtime << " seconds." << endl;
}
