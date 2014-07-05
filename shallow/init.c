
/************************************************************************
*									*
* Commonwealth Scientific and Industrial Research Organisation (CSIRO)	*
*	- Division of Information Technology	(DIT)			*
*	- Division of Atmospheric Research	(DAR)			*
*									*
* Shallow water weather model - Distributed Memory Version		*
*									*
* Finite difference model of shallow water equations based on :-	*
* "The dynamics of finite difference models of the shallow water	*
* equations" by R. Sadourney, JAS, 32, 1975.				*
* Code from:-								*
* "An introduction to three-dimensional climate modelling"		*
* by Washington and Parkinson						*
*									*
* Programmers	= David Abramson	(DIT) rcoda@koel.co.rmit.oz	*
*		= Paul Whiting		(DIT) rcopw@koel.co.rmit.oz	*
*		= Martin Dix		(DAR) mrd@koel.co.rmit.oz	*
* Language	= BSD c using Argonne NL macros				*
* O/S		= Unix System V						*
* H/W		= Encore Multimax 320					*
*									*
************************************************************************/
#include <stdio.h>
#include <math.h>
#include "decs.h"

void initialise(p, u, v, psi, pold, uold, vold, di, dj, z)
float p[n][m];
float u[n][m];
float v[n][m];
float psi[n][m];
float pold[n][m];
float uold[n][m];
float vold[n][m];
float di, dj;
float z[n][m];
{
  int	i,j,ip,jp;

  /* initialise values of the streamfunction */
  for (j = 0; j < n; j++){
    for (i = 0; i < m; i++){
      float sin1 = sin((double)((i+0.5)*di));
      float sin2 = sin((double)((j+0.5)*dj));
      psi[j][i] = a*sin1*sin2;
    }
  }

  /* initialise velocities */
  for (j = 0; j < n; j++){
    jp = (j+1) % n;
    for (i = 0; i < m; i++){
	ip = (i+1) % m;
	u[j][ip] = -(psi[jp][ip]-psi[j][ip])/dy;
	v[jp][i] = -(psi[jp][ip]-psi[jp][i])/dx;
    }
  }

  for (j = 0; j < n; j++) {
    for (i = 0; i < m; i++) {
      uold[j][i] = u[j][i];
      vold[j][i] = v[j][i];
      /* free surface height * gravitational acceleration */
      pold[j][i] = 50000.; 
      p[j][i] = 50000.;
    }
  }

  for (j = 0; j < n; j++) {
    for (i = 0; i < m; i++) {
      z[j][i] = 0.;
    }
  }


  printf("\n");
  printf("Shallow water weather model - Distributed Memory Version 0.6\n\n");
  printf("Number of points in the X direction%8d\n", n);
  printf("Number of points in the Y direction%8d\n", m);
  printf("Grid spacing in the X direction      %8.2f\n", dx);
  printf("Grid spacing in the Y direction      %8.2f\n", dy);
  printf("Time step                             %8.3f\n", dt);
  printf("Time filter parameter                 %8.3f\n", alpha);
} 
