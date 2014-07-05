
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

#include "decs.h"

void tstep(
int m_,
int n_,
float alpha_,
int jstart,
int jend,
float pold[n][m],
float uold[n][m],
float vold[n][m],
float p[n][m],
float u[n][m],
float v[n][m],
float pnew[n][m],
float unew[n][m],
float vnew[n][m],
float dpdt[n][m],
float dudt[n][m],
float dvdt[n][m],
int firststep,
float tdt)
{
  int i,j;

  for (j = jstart; j <= jend; j++){
    for (i = 0; i < m; i++){
      pnew[j][i] = pold[j][i] + tdt*dpdt[j][i]; 
      unew[j][i] = uold[j][i] + tdt*dudt[j][i];
      vnew[j][i] = vold[j][i] + tdt*dvdt[j][i];
    }
  }

  /* Don't apply time filter on first step */
  if ( !firststep ) {
    for (j = jstart; j <= jend; j++) {
      for (i = 0; i < m; i++) {
	pold[j][i] = p[j][i]+alpha*(pnew[j][i]-2.*p[j][i]+pold[j][i]);
	uold[j][i] = u[j][i]+alpha*(unew[j][i]-2.*u[j][i]+uold[j][i]);
	vold[j][i] = v[j][i]+alpha*(vnew[j][i]-2.*v[j][i]+vold[j][i]);
      }
    }
  }

  for (j = jstart; j <= jend; j++) {
    for (i = 0; i < m; i++) {
      p[j][i] = pnew[j][i];
      u[j][i] = unew[j][i];
      v[j][i] = vnew[j][i];
    }
  }
}
