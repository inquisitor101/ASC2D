#include "initial_structure.hpp"





CInitial::CInitial
(
 CConfig    	 *config_container,
 CGeometry  	 *geometry_container,
 CElement   	**element_container,
 unsigned short iZone
)
 /*
	* Constructor, used to initialize CInitial.
	*/
{
	// Assign zone ID.
	zoneID = iZone;
}


CInitial::~CInitial
(
 void
)
 /*
	* Destructor for CInitial class, frees allocated memory.
	*/
{

}


CGaussianInitial::CGaussianInitial
(
 CConfig    	 *config_container,
 CGeometry  	 *geometry_container,
 CElement   	**element_container,
 unsigned short iZone
)
	:
		CInitial
		(
		 config_container,
		 geometry_container,
		 element_container,
		 iZone
		)
 /*
	* Constructor, used to initialize CGaussianInitial.
	*/
{
  // Initialize Gaussian pulse center.
  x0 = config_container->GetCenterX0()[0];
  y0 = config_container->GetCenterX0()[1];
  // Initialize pulse strength, perctantage of background condition.
  A0 = config_container->GetDisturbanceRatio();
  // Initialize pulse radial distance.
  b  = config_container->GetDisturbanceWidth();
  // Pulse attenuation factor.
  kappa = log(2.0)/(b*b);

  // Background flow properties.
  Mach   = config_container->GetMachInf();
  Tinf   = 300.0;
  pInf   = 101325.0;
  // Deduce density.
  rhoInf = pInf/(GAS_CONSTANT*Tinf);

  // Reference speed of sound.
  aInf = sqrt(pInf*GAMMA/rhoInf);

  // Flow direction [degrees].
  theta  = config_container->GetFlowAngle();
  // Flow velocity.
  VelInf = Mach*aInf;
  uInf   = VelInf*cos(theta*PI_CONSTANT/180.0);
  vInf   = VelInf*sin(theta*PI_CONSTANT/180.0);

  // Initialize dispersion-relaxation correction.
  betaPML = uInf/( aInf*aInf - uInf*uInf );

	// Set the tuning constants to unity, by default. Implying a 2D pulse.
	sigmax = 1.0;
	sigmay = 1.0;

	// Check if this is a 1D pulse or not.
	if( config_container->GetTypeIC(iZone) == IC_GAUSSIAN_PRESSURE_1D_X ) sigmay = 0.0;
	if( config_container->GetTypeIC(iZone) == IC_GAUSSIAN_PRESSURE_1D_Y ) sigmax = 0.0;

	// Select whether or not to use a linear, 1-way acoustic pulse.
	// Note, this is only for a 1D pulse in the x-direction and moving
	// towards the right, i.e. w(+) = p' + rho*a*u'.
	if( config_container->GetTypeIC(iZone) == IC_GAUSSIAN_PRESSURE_1D_X ) 
		OneWayPulse = true;
	else
		OneWayPulse = false;
}


CGaussianInitial::~CGaussianInitial
(
 void
)
 /*
	* Destructor for CGaussianInitial class, frees allocated memory.
	*/
{

}


void CGaussianInitial::ComputeTargetStatePrimitive
(
  const as3data1d<as3double> &Coords,
  as3data1d<as3double>       &TargetState,
  unsigned short              nData
)
 /*
  * Function that computes the target state in a Gaussian pulse.
  */
{
  // Loop across all data points.
  for(unsigned short l=0; l<nData; l++){

    // Prescribe the primitive data.
    TargetState[0][l] = rhoInf;
    TargetState[1][l] = uInf;
    TargetState[2][l] = vInf;
    TargetState[3][l] = pInf;
  }
}


void CGaussianInitial::ComputeTargetStatePrimitivePerDOF
(
  const as3vector1d<as3double> &xy,
  as3vector1d<as3double>       &Qinf
)
 /*
  * Function that computes the target state in a Gaussian pulse.
	* This only does the calculation on a single DOF.
  */
{    
	// Prescribe the primitive data.
	Qinf[0] = rhoInf;
	Qinf[1] = uInf;
	Qinf[2] = vInf;
	Qinf[3] = pInf;
}


void CGaussianInitial::SetInitialCondition
(
  const as3data1d<as3double> &grid_nodes,
  as3data1d<as3double>       &data_nodes,
  unsigned short              nNode,
  as3double                   time
)
 /*
	* Function that sets a Gaussian initial condition in this zone.
	*/
{
  // Abbreviation.
  const as3double ovgm1 = 1.0/GAMMA_MINUS_ONE;
	const as3double ovra  = 1.0/(rhoInf*aInf);

	// Select whether to use a 1-way acoustic pulse or a generic pressure pulse.
	if( OneWayPulse )
	{	
		// Loop over every node and compute the solution.
		for(unsigned short l=0; l<nNode; l++){

			// Coordinates x and y.
			const as3double x = grid_nodes[0][l];
			const as3double y = grid_nodes[1][l];

			// Compute relative position w.r.t. pulse center.
			as3double rxPos2 = (x-x0); rxPos2 *= rxPos2;
			as3double ryPos2 = (y-y0); ryPos2 *= ryPos2;

			// Radial distance squared.
			const as3double r2 = sigmax*rxPos2 + sigmay*ryPos2;

			// Gaussian profile perturbation, based on background pressure.
			const as3double f = pInf*A0*exp( -kappa*r2 );

			// Compute initial condition in primitive form.
			const as3double rho = rhoInf;
			const as3double u   = uInf + f*ovra;
			const as3double v   = vInf;
			const as3double p   = pInf + f;
			// Compute total energy.
			const as3double rhoE = p*ovgm1 + 0.5*rho*( u*u + v*v );

			// Assemble conservative form.
			data_nodes[0][l] = rho;
			data_nodes[1][l] = rho*u;
			data_nodes[2][l] = rho*v;
			data_nodes[3][l] = rhoE;
		}
	}
	else
	{
		// Loop over every node and compute the solution.
		for(unsigned short l=0; l<nNode; l++){

			// Coordinates x and y.
			const as3double x = grid_nodes[0][l];
			const as3double y = grid_nodes[1][l];

			// Compute relative position w.r.t. pulse center.
			as3double rxPos2 = (x-x0); rxPos2 *= rxPos2;
			as3double ryPos2 = (y-y0); ryPos2 *= ryPos2;

			// Radial distance squared.
			const as3double r2 = sigmax*rxPos2 + sigmay*ryPos2;

			// Compute initial condition in primitive form.
			const as3double rho = rhoInf;
			const as3double u   = uInf;
			const as3double v   = vInf;
			const as3double p   = pInf*( 1.0 + A0*exp(-kappa*r2) );
			// Compute total energy.
			const as3double rhoE = p*ovgm1 + 0.5*rho*( u*u + v*v );

			// Assemble conservative form.
			data_nodes[0][l] = rho;
			data_nodes[1][l] = rho*u;
			data_nodes[2][l] = rho*v;
			data_nodes[3][l] = rhoE;
		}
	}
}


CIsentropicVortexInitial::CIsentropicVortexInitial
(
 CConfig    	 *config_container,
 CGeometry  	 *geometry_container,
 CElement   	**element_container,
 unsigned short iZone
)
	:
		CInitial
		(
		 config_container,
		 geometry_container,
		 element_container,
		 iZone
		)
 /*
	* Constructor, used to initialize CIsentropicVortexInitial.
	*/
{
  // Initialize vortex center.
  x0 = config_container->GetCenterX0()[0];
  y0 = config_container->GetCenterX0()[1];
  // Initialize vortex.
  A0 = config_container->GetDisturbanceRatio();
  // Initialize vortex radial distance.
  Rv = config_container->GetDisturbanceWidth();

  // Background flow properties.
  Mach   = config_container->GetMachInf();
  Tinf   = 300.0;
  pInf   = 101325.0;
  // Deduce density.
  rhoInf = pInf/(GAS_CONSTANT*Tinf);

  // Reference speed of sound.
  aInf = sqrt(pInf*GAMMA/rhoInf);

  // Flow direction [degrees].
  theta  = config_container->GetFlowAngle();
  // Flow velocity.
  VelInf = Mach*aInf;
  uInf   = VelInf*cos(theta*PI_CONSTANT/180.0);
  vInf   = VelInf*sin(theta*PI_CONSTANT/180.0);

  // Initialize dispersion-relaxation correction.
  betaPML = uInf/( aInf*aInf - uInf*uInf );
}


CIsentropicVortexInitial::~CIsentropicVortexInitial
(
 void
)
 /*
	* Destructor for CIsentropicVortexInitial class, frees allocated memory.
	*/
{

}


void CIsentropicVortexInitial::ComputeTargetStatePrimitive
(
  const as3data1d<as3double> &Coords,
  as3data1d<as3double>       &TargetState,
  unsigned short              nData
)
 /*
  * Function that computes the target state in an isentropic vortex.
  */
{
  // Loop across all data points.
  for(unsigned short l=0; l<nData; l++){

    // Prescribe the primitive data.
    TargetState[0][l] = rhoInf;
    TargetState[1][l] = uInf;
    TargetState[2][l] = vInf;
    TargetState[3][l] = pInf;
  }
}


void CIsentropicVortexInitial::ComputeTargetStatePrimitivePerDOF
(
  const as3vector1d<as3double> &xy,
  as3vector1d<as3double>       &Qinf
)
 /*
  * Function that computes the target state in a Gaussian pulse.
	* This only does the calculation on a single DOF.
  */
{
	// Prescribe the primitive data.
	Qinf[0] = rhoInf;
	Qinf[1] = uInf;
	Qinf[2] = vInf;
	Qinf[3] = pInf;
}


void CIsentropicVortexInitial::SetInitialCondition
(
  const as3data1d<as3double> &grid_nodes,
  as3data1d<as3double>       &data_nodes,
  unsigned short              nNode,
  as3double                   time
)
 /*
	* Function that sets an isentropic vortex initial condition in this zone.
	*/
{
  // Abbreviations.
  const as3double ovgm1 =  1.0/GAMMA_MINUS_ONE;
  const as3double ovRv2 =  1.0/(Rv*Rv);
  const as3double alpha =  A0/(aInf*Rv);
  const as3double omega = -0.5*GAMMA*alpha*alpha;

	// Loop over every node and set the solution.
	for(unsigned short l=0; l<nNode; l++){

		// Coordinates x and y.
		const as3double x = grid_nodes[0][l];
		const as3double y = grid_nodes[1][l];

		// Compute relative position w.r.t. vortex center.
		const as3double dx  = x-x0; const as3double dx2 = dx*dx;
		const as3double dy  = y-y0; const as3double dy2 = dy*dy;

		// Radial distance squared.
		const as3double r2 = dx2 + dy2;

    // Compute stream function.
    const as3double psi = A0*exp(-0.5*r2*ovRv2);

		// Compute initial condition in primitive form.
		const as3double u   = -ovRv2*dy*psi + uInf;
		const as3double v   =  ovRv2*dx*psi + vInf;
		const as3double p   =  pInf*exp( omega*exp(-r2*ovRv2) );
    const as3double rho =  p/(GAS_CONSTANT*Tinf);
		// Compute total energy.
		const as3double rhoE = p*ovgm1 + 0.5*rho*( u*u + v*v );

		// Assemble conservative form.
		data_nodes[0][l] = rho;
		data_nodes[1][l] = rho*u;
		data_nodes[2][l] = rho*v;
		data_nodes[3][l] = rhoE;
	}
}


CEntropyWave::CEntropyWave
(
 CConfig    	 *config_container,
 CGeometry  	 *geometry_container,
 CElement   	**element_container,
 unsigned short iZone
)
	:
		CInitial
		(
		 config_container,
		 geometry_container,
		 element_container,
		 iZone
		)
 /*
	* Constructor, used to initialize CEntropyWave.
	*/
{
  // Initialize Gaussian center.
  x0 = config_container->GetCenterX0()[0];
  y0 = config_container->GetCenterX0()[1];
  // Initialize pulse strength, perctantage of background condition.
  A0 = config_container->GetDisturbanceRatio();
  // Initialize pulse radial distance.
  b  = config_container->GetDisturbanceWidth();
  // Pulse attenuation factor.
  kappa = log(2.0)/(b*b);

  // Background flow properties.
  Mach   = config_container->GetMachInf();
  Tinf   = 300.0;
  pInf   = 101325.0;
  // Deduce density.
  rhoInf = pInf/(GAS_CONSTANT*Tinf);

  // Reference speed of sound.
  aInf = sqrt(pInf*GAMMA/rhoInf);

  // Flow direction [degrees].
  theta  = config_container->GetFlowAngle();
  // Flow velocity.
  VelInf = Mach*aInf;
  uInf   = VelInf*cos(theta*PI_CONSTANT/180.0);
  vInf   = VelInf*sin(theta*PI_CONSTANT/180.0);

  // Initialize dispersion-relaxation correction.
  betaPML = uInf/( aInf*aInf - uInf*uInf );
}


CEntropyWave::~CEntropyWave
(
 void
)
 /*
	* Destructor for CEntropyWave class, frees allocated memory.
	*/
{

}


void CEntropyWave::ComputeTargetStatePrimitive
(
  const as3data1d<as3double> &Coords,
  as3data1d<as3double>       &TargetState,
  unsigned short              nData
)
 /*
  * Function that computes the target state in an entropy wave.
  */
{
  // Loop across all data points.
  for(unsigned short l=0; l<nData; l++){

    // Prescribe the primitive data.
    TargetState[0][l] = rhoInf;
    TargetState[1][l] = uInf;
    TargetState[2][l] = vInf;
    TargetState[3][l] = pInf;
  }
}


void CEntropyWave::ComputeTargetStatePrimitivePerDOF
(
  const as3vector1d<as3double> &xy,
  as3vector1d<as3double>       &Qinf
)
 /*
  * Function that computes the target state in a Gaussian pulse.
	* This only does the calculation on a single DOF.
  */
{
	// Prescribe the primitive data.
	Qinf[0] = rhoInf;
	Qinf[1] = uInf;
	Qinf[2] = vInf;
	Qinf[3] = pInf;
}


void CEntropyWave::SetInitialCondition
(
  const as3data1d<as3double> &grid_nodes,
  as3data1d<as3double>       &data_nodes,
  unsigned short              nNode,
  as3double                   time
)
 /*
	* Function that sets an entropy wave initial condition in this zone.
	*/
{
  // Abbreviation.
  const as3double ovgm1 = 1.0/GAMMA_MINUS_ONE;
  const as3double ovrg  = 1.0/GAS_CONSTANT;

	// Loop over every node and compute the solution.
	for(unsigned short l=0; l<nNode; l++){

		// Coordinates x and y.
		const as3double x = grid_nodes[0][l];
		const as3double y = grid_nodes[1][l];

		// Compute relative position w.r.t. pulse center.
		as3double rxPos2 = (x-x0); rxPos2 *= rxPos2;
		as3double ryPos2 = (y-y0); ryPos2 *= ryPos2;

		// Radial distance squared.
		const as3double r2 = rxPos2 + ryPos2;

    // Set the mean primitive data.
    const as3double u = uInf;
		const as3double v = vInf;
		const as3double p = pInf;
    // Compute the temperature by superimposive a Gaussian pulse.
    const as3double T = Tinf*( 1.0 + A0*exp(-kappa*r2) );
    // Deduce the density.
    const as3double rho = ovrg*p/T;
		// Compute total energy.
		const as3double rhoE = p*ovgm1 + 0.5*rho*( u*u + v*v );

		// Assemble conservative form.
		data_nodes[0][l] = rho;
		data_nodes[1][l] = rho*u;
		data_nodes[2][l] = rho*v;
		data_nodes[3][l] = rhoE;
	}
}


CVortexRollup::CVortexRollup
(
 CConfig    	 *config_container,
 CGeometry  	 *geometry_container,
 CElement   	**element_container,
 unsigned short iZone
)
	:
		CInitial
		(
		 config_container,
		 geometry_container,
		 element_container,
		 iZone
		)
 /*
	* Constructor, used to initialize CVortexRollup.
  * NOTE, this is a non-dimensional problem, as defined in the paper by Hu.
  * If there need be dimensionalization, then the dispersion-relation needs
  * to be taken care of and not use the 1/4 ratio defined by the paper of Hu.
	*/
{
  // Background flow properties.
  Mach   = config_container->GetMachInf();
  Tinf   = 1.0;
  pInf   = 1.0/GAMMA;
  // Deduce density.
  rhoInf = 1.0;

  // Reference speed of sound.
  aInf   = 1.0;

  // Flow direction [degrees].
  theta  = 0.0;
  // Flow velocity.
  VelInf = 1.0;
  uInf   = 1.0;
  vInf   = 0.0;

  // Temperature on top.
  T1 = 1.0;
  // Temperature in bottom.
  T2 = 0.8;
  // Velocity on top.
  U1 = 0.8;
  // Velocity in bottom.
  U2 = 0.2;

  // Difference velocity.
  Ujmp  =     ( U1 - U2 );
  // Average velocity.
  Uavg  = 0.5*( U1 + U2 );
  // Initialize pulse strength, perctantage of background condition.
  A0    = config_container->GetDisturbanceRatio();
  // Flow width coefficient, around the shear layer center.
  delta = 0.4;

  // Initialize dispersion-relaxation correction.
  betaPML = 1.0/1.4;

  // Abbreviations.
  tovd   = 2.0/delta;
  ovujmp = 1.0/Ujmp;
  ovgm1  = 1.0/GAMMA_MINUS_ONE;
  gm1ov2 = 0.5*GAMMA_MINUS_ONE;
}


CVortexRollup::~CVortexRollup
(
 void
)
 /*
	* Destructor for CVortexRollup class, frees allocated memory.
	*/
{

}


void CVortexRollup::ComputeTargetStatePrimitive
(
  const as3data1d<as3double> &Coords,
  as3data1d<as3double>       &TargetState,
  unsigned short              nData
)
 /*
  * Function that computes the target state in a vortex roll-up.
  */
{
  // Loop across all data points.
#pragma omp simd
  for(unsigned short l=0; l<nData; l++){

    // Extract the y-coordinate.
    const as3double y = Coords[YDIM][l];

    // Compute the u-velocity.
    const as3double u   = Uavg + 0.5*Ujmp*tanh(y*tovd);
    // Compute temperature.
    const as3double T   = T1*ovujmp*(u-U2) + T2*ovujmp*(U1-u) + gm1ov2*(U1-u)*(u-U2);
    // Compute density.
    const as3double rho = 1.0/T;

    // Prescribe the primitive data.
    TargetState[0][l] = rho;
    TargetState[1][l] = u;
    TargetState[2][l] = vInf;
    TargetState[3][l] = pInf;
  }
}


void CVortexRollup::ComputeTargetStatePrimitivePerDOF
(
  const as3vector1d<as3double> &xy,
  as3vector1d<as3double>       &Qinf
)
 /*
  * Function that computes the target state in a Gaussian pulse.
	* This only does the calculation on a single DOF.
  */
{  
	// Loop across all data points.
    
	// Extract the y-coordinate.
  const as3double y = xy[1];

  // Compute the u-velocity.
  const as3double u   = Uavg + 0.5*Ujmp*tanh(y*tovd);
  // Compute temperature.
  const as3double T   = T1*ovujmp*(u-U2) + T2*ovujmp*(U1-u) + gm1ov2*(U1-u)*(u-U2);
  // Compute density.
  const as3double rho = 1.0/T;

  // Prescribe the primitive data.
  Qinf[0] = rho;
  Qinf[1] = u;
  Qinf[2] = vInf;
  Qinf[3] = pInf;
}


void CVortexRollup::SetInitialCondition
(
  const as3data1d<as3double> &grid_nodes,
  as3data1d<as3double>       &data_nodes,
  unsigned short              nNode,
  as3double                   time
)
 /*
	* Function that sets a vortex roll-up initial condition in this zone.
	*/
{
	// Loop over every node and compute the solution.
	for(unsigned short l=0; l<nNode; l++){

		// Coordinates x and y.
		const as3double y    = grid_nodes[1][l];
    // Compute the u-velocity.
    const as3double u    = Uavg + 0.5*Ujmp*tanh(y*tovd);
    // Compute temperature.
    const as3double T    = T1*ovujmp*(u-U2) + T2*ovujmp*(U1-u) + gm1ov2*(U1-u)*(u-U2);
    // Compute density.
    const as3double rho  = 1.0/T;
    // Compute pressure.
    const as3double p    = pInf;
    // Use background v-velocity.
    const as3double v    = vInf;
    // Compute total energy.
		const as3double rhoE = p*ovgm1 + 0.5*rho*( u*u + v*v );

		// Assemble conservative form.
		data_nodes[0][l] = rho;
		data_nodes[1][l] = rho*u;
		data_nodes[2][l] = rho*v;
		data_nodes[3][l] = rhoE;
	}
}


CAcousticPlane::CAcousticPlane
(
 CConfig    	 *config_container,
 CGeometry  	 *geometry_container,
 CElement   	**element_container,
 unsigned short iZone
)
	:
		CInitial
		(
		 config_container,
		 geometry_container,
		 element_container,
		 iZone
		)
 /*
	* Constructor, used to initialize CAcousticPlane.
	*/
{
  // Initialize Gaussian pulse center.
  x0 = config_container->GetCenterX0()[0];
  y0 = config_container->GetCenterX0()[1];
  // Time lag [sec].
  t0 = 20.0e-6;
  // Pressure disturbance coefficient.
  pa = 10.0;
  // Gaussian temporal-width [sec].
  st = 5.0e-6;
  // Initialize pulse radial distance.
  sx = config_container->GetDisturbanceWidth();
  // Pulse temporal attenuation factor.
  kappat = 0.5/(st*st);
  // Pulse spatial attenuation factor.
  kappax = 0.5/(sx*sx);

  // Background flow properties.
  Mach   = config_container->GetMachInf();
  Tinf   = 293.17;
  pInf   = 101325.0;
  // Deduce density.
  rhoInf = pInf/(GAS_CONSTANT*Tinf);

  // Reference speed of sound.
  aInf = sqrt(pInf*GAMMA/rhoInf);

  // Flow direction [degrees].
  theta  = config_container->GetFlowAngle();
  // Flow velocity.
  VelInf = Mach*aInf;
  uInf   = VelInf*cos(theta*PI_CONSTANT/180.0);
  vInf   = VelInf*sin(theta*PI_CONSTANT/180.0);

  // Initialize pulse strength amplification coefficient.
  A0 = pInf*config_container->GetDisturbanceRatio();

  // Initialize dispersion-relaxation correction.
  betaPML = uInf/( aInf*aInf - uInf*uInf );
}


CAcousticPlane::~CAcousticPlane
(
 void
)
 /*
	* Destructor for CAcousticPlane class, frees allocated memory.
	*/
{

}


void CAcousticPlane::ComputeTargetStatePrimitive
(
  const as3data1d<as3double> &Coords,
  as3data1d<as3double>       &TargetState,
  unsigned short              nData
)
 /*
  * Function that computes the target state in an acoustic Gaussian plane wave.
  */
{
  // Loop across all data points.
  for(unsigned short l=0; l<nData; l++){

    // Prescribe the primitive data.
    TargetState[0][l] = rhoInf;
    TargetState[1][l] = uInf;
    TargetState[2][l] = vInf;
    TargetState[3][l] = pInf;
  }
}


void CAcousticPlane::ComputeTargetStatePrimitivePerDOF
(
  const as3vector1d<as3double> &xy,
  as3vector1d<as3double>       &Qinf
)
 /*
  * Function that computes the target state in a Gaussian pulse.
	* This only does the calculation on a single DOF.
  */
{
	// Prescribe the primitive data.
	Qinf[0] = rhoInf;
	Qinf[1] = uInf;
	Qinf[2] = vInf;
	Qinf[3] = pInf;
}


void CAcousticPlane::SetInitialCondition
(
  const as3data1d<as3double> &grid_nodes,
  as3data1d<as3double>       &data_nodes,
  unsigned short              nNode,
  as3double                   time
)
 /*
	* Function that sets an acoustic Gaussian plane wave initial condition in this zone.
	*/
{
  // Abbreviation.
  const as3double ovgm1 = 1.0/GAMMA_MINUS_ONE;

	// Loop over every node and compute the solution.
	for(unsigned short l=0; l<nNode; l++){

    // Set the primitive variables.
    const as3double rho = rhoInf;
    const as3double u   = uInf;
    const as3double v   = vInf;
    const as3double p   = pInf;

		// Compute total energy.
		const as3double rhoE = p*ovgm1 + 0.5*rho*( u*u + v*v );

		// Assemble conservative form.
		data_nodes[0][l] = rho;
		data_nodes[1][l] = rho*u;
		data_nodes[2][l] = rho*v;
		data_nodes[3][l] = rhoE;
	}
}


