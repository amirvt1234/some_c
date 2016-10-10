//    An input file example to run a case based on GPUSPH code.

#include <cmath>
#include <iostream>
#include <stdexcept>

#include "NEES.h"
#include "particledefine.h"
#include "GlobalData.h"

#define MK_par 2

NEES::NEES(const GlobalData *_gdata) : Problem(_gdata)
{
	// Size and origin of the simulation domain
	lx = 42.0; //42
	ly = 2.65;  // 26.5
	lz = 2.0;
	m_size = make_double3(lx, ly, lz);
	m_origin = make_double3(0.0, 0.0, 0.0);

	// Data for problem setup
	h_length = 22.0f; //22
	slope_length = lx-h_length;
	height = .73f;
	beta = 4.2364*M_PI/180.0;

	// We have at least 1 moving boundary, the paddle
	m_mbnumber = 1;
	m_simparams.mbcallback = true;

	// Add objects to the tank
	icyl = 0;	// icyl = 0 means no cylinders
	icone = 0;	// icone = 0 means no cone
	// If presents, cylinders and cone are moving alltogether with
	// the same velocity
	if (icyl || icone)
		m_mbnumber++;

	i_use_bottom_plane = 1; // 1 for real plane instead of boundary parts

	// SPH parameters
	set_deltap(0.022f);  //0.005f;
	m_simparams.dt = 0.00013f;
	m_simparams.xsph = false;
	m_simparams.dtadapt = true;
	m_simparams.dtadaptfactor = 0.3;
	m_simparams.buildneibsfreq = 10;
	m_simparams.shepardfreq = 20;
	m_simparams.mlsfreq = 0;
	//m_simparams.visctype = ARTVISC;
	//m_simparams.visctype = KINEMATICVISC;
	m_simparams.visctype = SPSVISC;
	m_simparams.usedem = false;
	m_simparams.tend = 30.0;

	m_simparams.vorticity = true;
	m_simparams.boundarytype = LJ_BOUNDARY;  //LJ_BOUNDARY or MK_BOUNDARY

	// Physical parameters
	H = 0.73f;
	m_physparams.gravity = make_float3(0.0f, 0.0f, -9.81f);
	float g = length(m_physparams.gravity);

	m_physparams.set_density(0, 1000.0f, 7.0f, 20.f);
	m_physparams.numFluids = 1;
	float r0 = m_deltap;
	m_physparams.r0 = r0;

	m_physparams.artvisccoeff = 0.3f;
	m_physparams.kinematicvisc =  1.0e-6f;
	m_physparams.smagfactor = 0.12*0.12*m_deltap*m_deltap;
	m_physparams.kspsfactor = (2.0/3.0)*0.0066*m_deltap*m_deltap;
	m_physparams.epsartvisc = 0.01*m_simparams.slength*m_simparams.slength;

	// BC when using LJ
	m_physparams.dcoeff = 5.0f*g*H;
	//set p1coeff,p2coeff, epsxsph here if different from 12.,6., 0.5

	// BC when using MK
	m_physparams.MK_K = g*H;
	m_physparams.MK_d = 1.1*m_deltap/MK_par;
	m_physparams.MK_beta = MK_par;

	//Wave paddle definition:  location, start & stop times
	MbCallBack& mbpistondata = m_mbcallbackdata[0];
	mbpistondata.type = PISTONPART;
	mbpistondata.origin = make_float3(0.0, 0.0, 0.0);
	float amplitude = 0.2f;
	m_Hoh = amplitude/H;
	float kappa = sqrt((3*m_Hoh)/(4.0*H*H));
	float cel = sqrt(g*(H + amplitude));
	m_S = sqrt(16.0*amplitude*H/3.0);
	m_tau = 2.0*(3.8 + m_Hoh)/(kappa*cel);
//	std::cout << "m_tau: " << m_tau << "\n";
	mbpistondata.tstart = 0.0f;
	mbpistondata.tend = 4.0; //m_tau
	// Call mb_callback for piston a first time to initialise
	// values set by the call back function
	mb_callback(0.0, 0.0, 0);

	// Moving boundary initialisation data for cylinders and cone
	// used only if needed(cyl = 1 or cone = 1)
	MbCallBack& mbcyldata = m_mbcallbackdata[1];
	mbcyldata.type = GATEPART;
	mbcyldata.tstart = 0.0f;
	mbcyldata.tend =  1.0f;
	// Call mb_callback  for cylindres and cone a first time to initialise
	// values set by the call back function
	mb_callback(0.0, 0.0, 0);




	// Free surface detection
	m_simparams.surfaceparticle = true;
	m_simparams.savenormals = true;


	/*
	wavegages
	1-5: Resistance Wave Gages
	6-12: Sonic Run1
	13-19: Sonic Run2
	20-26: Sonic Run3
	Please note that all wave gages are located in l_y/2 despite the real position in experiment!
	*/
	add_gage(3.39, ly/2.);
	add_gage(8.27, ly/2.);
	add_gage(14.37, ly/2.);
	add_gage(19.25, ly/2.);
	add_gage(26.35, ly/2.);

	add_gage(36.54, ly/2.);
	add_gage(35.00, ly/2.);
	add_gage(34.18, ly/2.);
	add_gage(38.30, ly/2.);
	add_gage(35.44, ly/2.);
	add_gage(39.99, ly/2.);
	add_gage(37.20, ly/2.);

	add_gage(37.22, ly/2.);
	add_gage(35.00, ly/2.);
	add_gage(34.00, ly/2.);
	add_gage(38.29, ly/2.);
	add_gage(36.10, ly/2.);
	add_gage(39.45, ly/2.);
	add_gage(36.77, ly/2.);

	add_gage(32.79, ly/2.);
	add_gage(32.78, ly/2.);
	add_gage(35.00, ly/2.);
	add_gage(37.67, ly/2.);
	add_gage(36.11, ly/2.);
	add_gage(38.95, ly/2.);
	add_gage(37.20, ly/2.);

	// Drawing and saving times
	set_timer_tick(0.001f);
	add_writer(VTKWRITER, 10);

	// Name of problem used for directory creation
	m_name = "NEES";
}


NEES::~NEES(void)
{
	release_memory();
}


void NEES::release_memory(void)
{
	parts.clear();
	boundary_parts.clear();
	gate_parts.clear();
	piston_parts.clear();
}


MbCallBack& NEES::mb_callback(const float t, const float dt, const int i)
{
	switch (i) {
		// Piston
		case 0:
			{
			MbCallBack& mbpistondata = m_mbcallbackdata[0];
			mbpistondata.type = PISTONPART;
			const float posx = mbpistondata.origin.x;
			if (t >= mbpistondata.tstart && t < mbpistondata.tend) {
				//float arg = 2.0*((3.8 + m_Hoh)*((t - mbpistondata.tstart)/m_tau - 0.5)
				//			- 2.0*m_Hoh*((posx/m_S) - 0.5));
				//mbpistondata.disp.x = m_S*(1.0 + tanh(arg))/2.0;
				//mbpistondata.vel.x = (3.8 + m_Hoh)*m_S/(m_tau*cosh(arg)*cosh(arg));
				//float arg = 2.0*((3.8 + m_Hoh)*((t - mbpistondata.tstart)/m_tau - 0.5)
				//			- 2.0*m_Hoh*((posx/m_S) - 0.5));
				mbpistondata.disp.x = 1.0 + erf((t-2.0)/0.8);
				mbpistondata.vel.x = 2.0/(0.8*sqrt(3.14))*exp(-(t-2.0)/0.8*(t-2.0)/0.8   );
			}
			else {
				mbpistondata.vel.x = 0;
				}
			}
			break;

		// Cylinders and cone
		case 1:
			{
			MbCallBack& mbcyldata = m_mbcallbackdata[1];
			if (t >= mbcyldata.tstart && t < mbcyldata.tend) {
				mbcyldata.vel = make_float3(0.0f, 0.0f, 0.5f);
				mbcyldata.disp += mbcyldata.vel*dt;
				}
			else
				mbcyldata.vel = make_float3(0.0f, 0.0f, 0.0f);
			break;
			}

		default:
			throw runtime_error("Incorrect moving boundary object number");
			break;
		}

	return m_mbcallbackdata[i];
}


int NEES::fill_parts()
{
	const float r0 = m_physparams.r0;
	const float width = ly;

	const float br = (m_simparams.boundarytype == MK_BOUNDARY ? m_deltap/MK_par : r0);

	experiment_box = Cube(Point(0, 0, 0), Vector(h_length + slope_length, 0, 0),
                     Vector(0, width, 0), Vector(0, 0, height));

	boundary_parts.reserve(100);
	parts.reserve(34000);
	gate_parts.reserve(2000);
	piston_parts.reserve(500);

	MbCallBack& mbpistondata = m_mbcallbackdata[0];
	Rect piston = Rect(Point(mbpistondata.origin),
						Vector(0, width, 0), Vector(0, 0, lz));
	piston.SetPartMass(m_deltap, m_physparams.rho0[0]);
	piston.Fill(piston_parts, br, true);

	experiment_box1 = Rect(Point(0,0,0  ), Vector(0, ly, 0),
                          Vector(h_length, 0.0, 0.0));
        experiment_box1.SetPartMass(br, m_physparams.rho0[0]);
        experiment_box1.Fill(boundary_parts,br,true);
	experiment_box2 = Rect(Point(h_length,0,0  ), Vector(0, ly, 0),
                          Vector(slope_length, 0.0, slope_length*tan(beta)));
        experiment_box2.SetPartMass(br, m_physparams.rho0[0]);
        experiment_box2.Fill(boundary_parts,br,true);

	if (i_use_bottom_plane == 0) {
	   experiment_box1 = Rect(Point(h_length, 0, 0), Vector(0, width, 0),
			Vector(slope_length/cos(beta), 0.0, slope_length*tan(beta)));
	   experiment_box1.SetPartMass(m_deltap, m_physparams.rho0[0]);
	   experiment_box1.Fill(boundary_parts,br,true);
	   std::cout << "bottom rectangle defined" <<"\n";
	   }



	Rect fluid;
	float z = 0;
	int n = 0;
	while (z < H) {
		z = n*m_deltap + 1.5*r0;    //z = n*m_deltap + 1.5*r0;
		float x = mbpistondata.origin.x + r0;
		float l = h_length + z/tan(beta) - 1.5*r0/sin(beta) - x;
		fluid = Rect(Point(x,  r0, z),
				Vector(0, width-2.0*r0, 0), Vector(l, 0, 0));
		fluid.SetPartMass(m_deltap, m_physparams.rho0[0]);
		fluid.Fill(parts, m_deltap, true);
		n++;
	 }

    return parts.size() + boundary_parts.size() + gate_parts.size() + piston_parts.size();
}


uint NEES::fill_planes()
{

    if (i_use_bottom_plane == 0) {
		return 5;
		}
	else {
		return 6;
		} //corresponds to number of planes
}


void NEES::copy_planes(float4 *planes, float *planediv)
{
	const float w = m_size.y;
	const float l = h_length + slope_length;

	//  plane is defined as a x + by +c z + d= 0
	planes[0] = make_float4(0, 0, 1.0, 0);   //bottom, where the first three numbers are the normal, and the last is d.
	planediv[0] = 1.0;
	planes[1] = make_float4(0, 1.0, 0, 0);   //wall
	planediv[1] = 1.0;
	planes[2] = make_float4(0, -1.0, 0, w); //far wall
	planediv[2] = 1.0;
	planes[3] = make_float4(1.0, 0, 0, 0);  //end
	planediv[3] = 1.0;
	planes[4] = make_float4(-1.0, 0, 0, l);  //one end
	planediv[4] = 1.0;
	if (i_use_bottom_plane == 1)  {
		planes[5] = make_float4(-sin(beta),0,cos(beta), h_length*sin(beta));  //sloping bottom starting at x=h_length
		planediv[5] = 1.0;
	}
}


void NEES::copy_to_array(BufferList &buffers)
{
	float4 *pos = buffers.getData<BUFFER_POS>();
	hashKey *hash = buffers.getData<BUFFER_HASH>();
	float4 *vel = buffers.getData<BUFFER_VEL>();
	particleinfo *info = buffers.getData<BUFFER_INFO>();

	std::cout << "\nBoundary parts: " << boundary_parts.size() << "\n";
		std::cout << "      "<< 0  <<"--"<< boundary_parts.size() << "\n";
	for (uint i = 0; i < boundary_parts.size(); i++) {
		vel[i] = make_float4(0, 0, 0, m_physparams.rho0[0]);
		info[i]= make_particleinfo(BOUNDPART, 0, i);  // first is type, object, 3rd id
		calc_localpos_and_hash(boundary_parts[i], info[i], pos[i], hash[i]);
	}
	int j = boundary_parts.size();
	std::cout << "Boundary part mass:" << pos[j-1].w << "\n";

	std::cout << "\nPiston parts: " << piston_parts.size() << "\n";
	std::cout << "     " << j << "--" << j + piston_parts.size() << "\n";
	for (uint i = j; i < j + piston_parts.size(); i++) {
		vel[i] = make_float4(0, 0, 0, m_physparams.rho0[0]);
		info[i] = make_particleinfo(PISTONPART, 0, i);
		calc_localpos_and_hash(piston_parts[i - j], info[i], pos[i], hash[i]);
	}
	j += piston_parts.size();
	std::cout << "Piston part mass:" << pos[j-1].w << "\n";

	std::cout << "\nGate parts: " << gate_parts.size() << "\n";
	std::cout << "       " << j << "--" << j+gate_parts.size() <<"\n";
	for (uint i = j; i < j + gate_parts.size(); i++) {
		vel[i] = make_float4(0, 0, 0, m_physparams.rho0[0]);
		info[i] = make_particleinfo(GATEPART, 1, i);
		calc_localpos_and_hash(gate_parts[i - j], info[i], pos[i], hash[i]);
	}
	j += gate_parts.size();
	std::cout << "Gate part mass:" << pos[j-1].w << "\n";


	std::cout << "\nFluid parts: " << parts.size() << "\n";
	std::cout << "      "<< j  <<"--"<< j+ parts.size() << "\n";
	for (uint i = j; i < j + parts.size(); i++) {
		vel[i] = make_float4(0, 0, 0, m_physparams.rho0[0]);
	    info[i]= make_particleinfo(FLUIDPART,0,i);
		calc_localpos_and_hash(parts[i - j], info[i], pos[i], hash[i]);
		// initializing density
		//       float rho = m_physparams.rho0*pow(1.+g*(H-pos[i].z)/m_physparams.bcoeff,1/m_physparams.gammacoeff);
		//        vel[i] = make_float4(0, 0, 0, rho);
	}
	j += parts.size();
	std::cout << "Fluid part mass:" << pos[j-1].w << "\n";

	std::cout << " Everything uploaded" <<"\n";
}
#undef MK_par
