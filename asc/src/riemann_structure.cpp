#include "riemann_structure.hpp"




CRiemann::CRiemann
(
  CConfig *config_container
)
 /*
  * Constructor, used to initialize CRiemann.
  */
{
  // Abbreviation involving specific heat ratio.
  gm1 = GAMMA_MINUS_ONE;
}


CRiemann::~CRiemann
(
  void
)
 /*
  * Destructor for CRiemann class, frees allocated memory
  */
{

}


CRoeRiemann::CRoeRiemann
(
  CConfig *config_container
)
	:
		CRiemann
		(
		 config_container
		)
 /*
  * Constructor, used to initialize CRoeRiemann.
  */
{
  // Entropy fix coefficient.
  Delta = 0.001;
}


CRoeRiemann::~CRoeRiemann
(
  void
)
 /*
  * Destructor for CRoeRiemann class, frees allocated memory
  */
{

}


void CRiemann::ComputeVariableStateUpwinding
(
  const as3vector1d<as3double> &UnitNormal,
  const as3vector1d<as3double> &weights,
  const as3vector1d<as3double> &hElem,
  as3double                     Velocity,
  as3double                   **VarI,
  as3double                   **VarJ,
  as3double                   **VarF
)
 /*
  * Function that computes the variable at a given boundary face using an upwinded scheme.
  */
{
  // Deduce number of integration points on face.
  unsigned short nDOFsInt1D = weights.size();

  // Extract unit-normal.
  const as3double nx = UnitNormal[XDIM];
  const as3double ny = UnitNormal[YDIM];

  // Some abbreviations.
  const as3double unovh = Velocity*( nx/hElem[0] + ny/hElem[1] );
  // Determine the sign of the normal-velocity component.
  const as3double sign  = ( fabs(unovh) > 1.0e-15 ) ? unovh/fabs(unovh) : 0.0;

  // Use full-upwinding to determine the boundary state.
  for(unsigned short i=0; i<nVar; i++){
#pragma omp simd
    for(unsigned short l=0; l<nDOFsInt1D; l++){

      // Compute the weights and Jacobian. Note, the factor 2.0 in the
      // Jacobian is ignored, since it is being multiplied by a 0.5 in the
      // upwinding/central-differencing step.
      const as3double scale = -unovh*weights[l];

      // Compute upwinding. Note, the half is removed because it is
      // implicitly taken into account by the 2.0 in the Jacobian.
      VarF[i][l] = scale*(       ( VarJ[i][l] + VarI[i][l] )
                 -          sign*( VarJ[i][l] - VarI[i][l] ) );
    }
  }
}


void CRoeRiemann::ComputeFluxState
(
  const as3vector1d<as3double> &UnitNormal,
  const as3vector1d<as3double> &weights,
  const as3vector1d<as3double> &hElem,
  as3double                   **VarI,
  as3double                   **VarJ,
  as3double                   **Flux
)
 /*
  * Function that computes the flux at a given boundary face using the Roe scheme.
  */
{
  // Deduce number of integration points on face.
  unsigned short nDOFsInt1D = weights.size();

  // Extract unit-normal.
  const as3double nx = UnitNormal[XDIM];
  const as3double ny = UnitNormal[YDIM];

  // Length of the face, included is the inverse matrix jacobian.
  const as3double lenFace = -fabs( nx/hElem[XDIM] + ny/hElem[YDIM] );

  // Loop over all integration points on the face and compute the flux.
#pragma omp simd
  for(unsigned short l=0; l<nDOFsInt1D; l++){

    // Multiplication factor, used in flux computations.
    const as3double halfArea = weights[l]*lenFace;

    // Compute primitive data for internal face, i.e. left state.
  	as3double tmp 	  	= 1.0/VarI[0][l];
  	const as3double	vxL = tmp*VarI[1][l];
  	const as3double vyL = tmp*VarI[2][l];
  	const as3double pL  = gm1*(VarI[3][l]
                        - 0.5*(vxL*VarI[1][l] + vyL*VarI[2][l]) );

  	// Compute primitive data for external face, i.e. right state.
  	tmp                 = 1.0/VarJ[0][l];
  	const as3double vxR = tmp*VarJ[1][l];
  	const as3double vyR = tmp*VarJ[2][l];
  	const as3double pR  = gm1*(VarJ[3][l]
                        - 0.5*(vxR*VarJ[1][l] + vyR*VarJ[2][l]) );

  	// Compute the difference of the conservative mean flow variables.
  	const as3double dr  = VarJ[0][l] - VarI[0][l];
  	const as3double dru = VarJ[1][l] - VarI[1][l];
  	const as3double drv = VarJ[2][l] - VarI[2][l];
  	const as3double drE = VarJ[3][l] - VarI[3][l];

  	// Compute the Roe average state.
  	const as3double zL = sqrt(VarI[0][l]);
  	const as3double zR = sqrt(VarJ[0][l]);
  	tmp                = 1.0/(zL + zR);

  	const as3double rHL = VarI[3][l] + pL;
  	const as3double rHR = VarJ[3][l] + pR;

  	const as3double uAvg = tmp*(zL*vxL + zR*vxR);
  	const as3double vAvg = tmp*(zL*vyL + zR*vyR);
  	const as3double HAvg = tmp*(rHL/zL + rHR/zR);

    // Compute some abbreviations.
  	const as3double alphaAvg = 0.5*(uAvg*uAvg + vAvg*vAvg);
  	tmp                      = gm1*(HAvg - alphaAvg);
  	const as3double a2Avg    = fabs(tmp);
  	const as3double aAvg     = sqrt(a2Avg);
  	const as3double vnAvg  	 = uAvg*nx + vAvg*ny;
  	const as3double ovaAvg   = 1.0/aAvg;
  	const as3double ova2Avg  = 1.0/a2Avg;

  	// Compute absolute values of eigenvalues.
  	as3double lam1 = fabs(vnAvg + aAvg);
  	as3double lam2 = fabs(vnAvg - aAvg);
  	as3double lam3 = fabs(vnAvg);

    // Apply entropy fix.
    tmp  = Delta*std::max(lam1, lam2);
  	lam1 = std::max(lam1, tmp);
  	lam2 = std::max(lam2, tmp);
  	lam3 = std::max(lam3, tmp);

  	// Some abbreviations.
  	const as3double abv1 = 0.5*(lam1 + lam2);
  	const as3double abv2 = 0.5*(lam1 - lam2);
  	const as3double abv3 = abv1 - lam3;

  	const as3double abv4 = gm1*(alphaAvg*dr  - uAvg*dru - vAvg*drv + drE);
  	const as3double abv5 = nx*dru + ny*drv   - vnAvg*dr;
  	const as3double abv6 = abv3*abv4*ova2Avg + abv2*abv5*ovaAvg;
  	const as3double abv7 = abv2*abv4*ovaAvg  + abv3*abv5;

  	// Compute Roe flux vector: 0.5*( FR + FL - |A|*(UR - UL) ).
  	const as3double vnL = vxL*nx + vyL*ny;
  	const as3double vnR = vxR*nx + vyR*ny;
  	const as3double pa  = pL + pR;

    // There should be a 0.5 here, but it was cancelled by the 2.0 in the jacobian.
    // That is why lenFace is not O( 2.0/h ).
    Flux[0][l] = halfArea*(VarI[0][l]*vnL + VarJ[0][l]*vnR
               - 				  (lam3*dr        + abv6));
    Flux[1][l] = halfArea*(VarI[1][l]*vnL + VarJ[1][l]*vnR + pa*nx
               -          (lam3*dru       + uAvg*abv6      + nx*abv7));
    Flux[2][l] = halfArea*(VarI[2][l]*vnL + VarJ[2][l]*vnR + pa*ny
               -          (lam3*drv       + vAvg*abv6      + ny*abv7));
    Flux[3][l] = halfArea*(VarI[3][l]*vnL + VarJ[3][l]*vnR + pL*vnL + pR*vnR
               -          (lam3*drE       + HAvg*abv6      + vnAvg*abv7));
  }
}


CRusanovRiemann::CRusanovRiemann
(
  CConfig *config_container
)
	:
		CRiemann
		(
		 config_container
		)
 /*
  * Constructor, used to initialize CRusanovRiemann.
  */
{

}


CRusanovRiemann::~CRusanovRiemann
(
  void
)
 /*
  * Destructor for CRusanovRiemann class, frees allocated memory
  */
{

}


void CRusanovRiemann::ComputeFluxState
(
  const as3vector1d<as3double> &UnitNormal,
  const as3vector1d<as3double> &weights,
  const as3vector1d<as3double> &hElem,
  as3double                   **VarI,
  as3double                   **VarJ,
  as3double                   **Flux
)
 /*
  * Function that computes the flux at a given boundary face using the Rusanov scheme.
  */
{
  // Deduce number of integration points on face.
  unsigned short nDOFsInt1D = weights.size();

  // Extract unit-normal.
  const as3double nx = UnitNormal[XDIM];
  const as3double ny = UnitNormal[YDIM];

  // Length of the face, included is the inverse matrix jacobian.
  const as3double lenFace = -fabs( nx/hElem[XDIM] + ny/hElem[YDIM] );

  // Loop over all integration points on the face and compute the flux.
#pragma omp simd
  for(unsigned short l=0; l<nDOFsInt1D; l++){

    // Multiplication factor, used in flux computations.
    const as3double halfArea = weights[l]*lenFace;

    // Compute primitive data for internal face, i.e. left state.
    as3double tmp 	  	= 1.0/VarI[0][l];
    const as3double	vxL = tmp*VarI[1][l];
    const as3double vyL = tmp*VarI[2][l];
    const as3double pL  = gm1*(VarI[3][l]
                        - 0.5*(vxL*VarI[1][l] + vyL*VarI[2][l]) );
    // Compute the speed-of-sound on the left-state.
    const as3double aL  = sqrt(GAMMA*pL*tmp);

    // Compute primitive data for external face, i.e. right state.
    tmp                 = 1.0/VarJ[0][l];
    const as3double vxR = tmp*VarJ[1][l];
    const as3double vyR = tmp*VarJ[2][l];
    const as3double pR  = gm1*(VarJ[3][l]
                        - 0.5*(vxR*VarJ[1][l] + vyR*VarJ[2][l]) );
    // Compute the speed-of-sound on the right-state.
    const as3double aR  = sqrt(GAMMA*pR*tmp);

    // Compute the difference of the conservative mean flow variables.
    const as3double dr  = VarJ[0][l] - VarI[0][l];
    const as3double dru = VarJ[1][l] - VarI[1][l];
    const as3double drv = VarJ[2][l] - VarI[2][l];
    const as3double drE = VarJ[3][l] - VarI[3][l];

  	// Compute Rusanov flux vector: 0.5*( FR + FL - max(lambda)*(UR - UL) ).
  	const as3double vnL = vxL*nx + vyL*ny;
  	const as3double vnR = vxR*nx + vyR*ny;
  	const as3double pa  = pL + pR;

    // Maximum eigenvalue, in between the two states.
    const as3double maxlam = std::max( fabs(vnL)+aL, fabs(vnR)+aR );

    // There should be a 0.5 here, but it was cancelled by the 2.0 in the jacobian.
    // That is why lenFace is not O( 2.0/h ).
    Flux[0][l] = halfArea*(VarI[0][l]*vnL + VarJ[0][l]*vnR
               - 				   maxlam*dr);
    Flux[1][l] = halfArea*(VarI[1][l]*vnL + VarJ[1][l]*vnR + pa*nx
               -           maxlam*dru);
    Flux[2][l] = halfArea*(VarI[2][l]*vnL + VarJ[2][l]*vnR + pa*ny
               -           maxlam*drv);
    Flux[3][l] = halfArea*(VarI[3][l]*vnL + VarJ[3][l]*vnR + pL*vnL + pR*vnR
               -           maxlam*drE);
  }
}


CRoeIsmailRiemann::CRoeIsmailRiemann
(
  CConfig *config_container
)
	:
		CRiemann
		(
		 config_container
		)
 /*
  * Constructor, used to initialize CRoeIsmailRiemann.
  */
{
  // Abbreviations involving gamma.
  ovgm1  =  1.0/(GAMMA-1.0);
  gp1Ovg = (GAMMA+1.0)/GAMMA;
  gm1Ovg =  gm1/GAMMA;

  // Values to scale the acoustic eigenvalues to obtain an adequate amount
  // of dissipation to be entropy satisfying in the Ismail_Roe flux.
  // Note that only alphaMax is taken, assuming that the jump in Mach number
  // over the interface is less than 0.5. For alphaMax = 0 the EC1 flux is obtained.
  beta     = 1.0/6.0;
  // alphaMax = 2.0;
  alphaMax = 0.0;
}


CRoeIsmailRiemann::~CRoeIsmailRiemann
(
  void
)
 /*
  * Destructor for CRoeIsmailRiemann class, frees allocated memory
  */
{

}


void CRoeIsmailRiemann::ComputeFluxState
(
  const as3vector1d<as3double> &UnitNormal,
  const as3vector1d<as3double> &weights,
  const as3vector1d<as3double> &hElem,
  as3double                   **VarI,
  as3double                   **VarJ,
  as3double                   **Flux
)
 /*
  * Function that computes the flux at a given boundary face using the Roe-Ismail scheme.
  */
{
  // Deduce number of integration points on face.
  unsigned short nDOFsInt1D = weights.size();

  // Extract unit-normal.
  const as3double nx = UnitNormal[XDIM];
  const as3double ny = UnitNormal[YDIM];

  // Length of the face, included is the inverse matrix jacobian.
  const as3double lenFace = -2.0*fabs( nx/hElem[XDIM] + ny/hElem[YDIM] );

  // Loop over all integration points on the face and compute the flux.
#pragma omp simd
  for(unsigned short l=0; l<nDOFsInt1D; l++){

    // Multiplication factor, used in flux computations.
    const as3double Area = weights[l]*lenFace;

    // Compute the primitive variables of the left and right state.
    as3double tmp 	  	= 1.0/VarI[0][l];
    const as3double	vxL = tmp*VarI[1][l];
    const as3double vyL = tmp*VarI[2][l];
    const as3double pL  = gm1*(VarI[3][l]
                        - 0.5*(vxL*VarI[1][l] + vyL*VarI[2][l]) );

    // Compute primitive data for external face, i.e. right state.
    tmp                 = 1.0/VarJ[0][l];
    const as3double vxR = tmp*VarJ[1][l];
    const as3double vyR = tmp*VarJ[2][l];
    const as3double pR  = gm1*(VarJ[3][l]
                        - 0.5*(vxR*VarJ[1][l] + vyR*VarJ[2][l]) );

    const as3double rhoPInvL = VarI[0][l]/pL;
    const as3double rhoPInvR = VarJ[0][l]/pR;

    // Compute the entropy variables of the left and right state.
    tmp = log(pL/pow(VarI[0][l], GAMMA));

    const as3double V0L =  (GAMMA-tmp)*ovgm1
                        -  0.5*rhoPInvL*(vxL*vxL + vyL*vyL);
    const as3double V1L =  rhoPInvL*vxL;
    const as3double V2L =  rhoPInvL*vyL;
    const as3double V3L = -rhoPInvL;

    tmp = log(pR/pow(VarJ[0][l], GAMMA));

    const as3double V0R =  (GAMMA-tmp)*ovgm1
                        -  0.5*rhoPInvR*(vxR*vxR + vyR*vyR);
    const as3double V1R =  rhoPInvR*vxR;
    const as3double V2R =  rhoPInvR*vyR;
    const as3double V3R = -rhoPInvR;

    // Compute the difference in entropy variables.
    const as3double dV0 = V0R - V0L, dV1 = V1R - V1L,
                    dV2 = V2R - V2L, dV3 = V3R - V3L;

    // Compute the z-variables of the left and right states.
    const as3double z0L = sqrt(rhoPInvL),      z0R = sqrt(rhoPInvR);
    const as3double z1L = vxL*z0L,             z1R = vxR*z0R;
    const as3double z2L = vyL*z0L,             z2R = vyR*z0R;
    const as3double z3L = sqrt(VarI[0][l]*pL), z3R = sqrt(VarJ[0][l]*pR);

    // Compute the arithmetic average of the z-variables.
    const as3double z0Avg = 0.5*(z0L + z0R);
    const as3double z1Avg = 0.5*(z1L + z1R);
    const as3double z2Avg = 0.5*(z2L + z2R);
    const as3double z3Avg = 0.5*(z3L + z3R);

    // Compute the logarithmic mean of z0.
    as3double zeta, f, u, F;

    zeta = z0L/z0R;
    f    = (zeta-1.0)/(zeta+1.0);
    u    = f*f;
    if(u < (as3double) 0.01)
      F = 1.0 + u/3.0 + u*u/5.0 + u*u*u/7.0 + u*u*u*u/9.0;
    else
      F = log(zeta)/(2.0*f);

    const as3double z0LogAvg = z0Avg/F;

    // Compute the logarithmic mean of z4.
    zeta = z3L/z3R;
    f    = (zeta-1.0)/(zeta+1.0);
    u    = f*f;
    if(u < (as3double) 0.01)
      F = 1.0 + u/3.0 + u*u/5.0 + u*u*u/7.0 + u*u*u*u/9.0;
    else
      F = log(zeta)/(2.0*f);

    const as3double z3LogAvg = z3Avg/F;

    // Compute the other averaged quantities that are necessary.
    const as3double oneOvz0Avg = 1.0/z0Avg;
    const as3double rhoAvg = z0Avg*z3LogAvg;
    const as3double p1Avg  = oneOvz0Avg*z3Avg;
    const as3double p2Avg  = 0.5*(gp1Ovg*z3LogAvg/z0LogAvg + gm1Ovg*p1Avg);
    const as3double vxAvg  = oneOvz0Avg*z1Avg;
    const as3double vyAvg  = oneOvz0Avg*z2Avg;

    const as3double vnAvg  = vxAvg*nx + vyAvg*ny;
    const as3double kinAvg = 0.5*(vxAvg*vxAvg + vyAvg*vyAvg);
    const as3double a2Avg  = GAMMA*p2Avg/rhoAvg;
    const as3double aAvg   = sqrt(a2Avg);
    const as3double HAvg   = a2Avg*ovgm1 + kinAvg;
    const as3double EAvg   = HAvg - p2Avg/rhoAvg;

    const as3double ovaAvg  = 1.0/aAvg;
    const as3double ova2Avg = 1.0/a2Avg;

    // Define the difference in conservative variables as dU/dV deltaV, where
    // the transformation matrix dU/dV must be evaluated at the averaged state.
    as3double abv1 = rhoAvg*(vxAvg*dV1 + vyAvg*dV2);
    as3double abv2 = abv1 + rhoAvg*(dV0 + HAvg*dV3);

    const as3double dr  = abv1 + rhoAvg*(dV0 + EAvg*dV3);
    const as3double dru = vxAvg*abv2 + p2Avg*dV1;
    const as3double drv = vyAvg*abv2 + p2Avg*dV2;
    const as3double drE = HAvg*abv1  + rhoAvg*EAvg*(dV0 + HAvg*dV3)
                        + p2Avg*kinAvg*dV3;

    // Compute the absolute values of the eigenvalues of the flux Jacobian.
    as3double lam1 = fabs(vnAvg + aAvg);
    as3double lam2 = fabs(vnAvg - aAvg);
    as3double lam3 = fabs(vnAvg);

    // Scale the acoustic eigenvalue, such that the EC2 (or EC1) flux of Ismail
    // and Roe is obtained. Also multiply all eigenvalues by half to obtain
    // the correct scaling.
    lam1 *= 0.5*(1.0 + beta + alphaMax);
    lam2 *= 0.5*(1.0 + beta + alphaMax);
    lam3 *= 0.5;

    // Some abbreviations, which occur quite often in the dissipation terms.
    abv1 = 0.5*(lam1 + lam2);
    abv2 = 0.5*(lam1 - lam2);

    const as3double abv3 = abv1 - lam3;
    const as3double abv4 = gm1*(kinAvg*dr - vxAvg*dru - vyAvg*drv + drE);
    const as3double abv5 = nx*dru + ny*drv - vnAvg*dr;
    const as3double abv6 = abv3*abv4*ova2Avg + abv2*abv5*ovaAvg;
    const as3double abv7 = abv2*abv4*ovaAvg  + abv3*abv5;

    // Compute the central flux.
    Flux[0][l] = rhoAvg*vnAvg;
    Flux[1][l] = rhoAvg*vnAvg*vxAvg + p1Avg*nx;
    Flux[2][l] = rhoAvg*vnAvg*vyAvg + p1Avg*ny;
    Flux[3][l] = rhoAvg*vnAvg*HAvg;

    // Subtract the dissipation terms.
    Flux[0][l] -= (lam3*dr + abv6);
    Flux[1][l] -= (lam3*dru + vxAvg*abv6 + nx*abv7);
    Flux[2][l] -= (lam3*drv + vyAvg*abv6 + ny*abv7);
    Flux[3][l] -= (lam3*drE + HAvg*abv6 + vnAvg*abv7);

    // Scale the fluxes with Area.
    Flux[0][l] *= Area;
    Flux[1][l] *= Area;
    Flux[2][l] *= Area;
    Flux[3][l] *= Area;
  }
}

