
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

#define	m	512	/* 18 Number of points in x direction */
#define	n	512	/* 18 Number of points in y direction */
#define	a	1.e6	/* Nominally the radius of the earth but here just
			a length scale*/
#define	dt	90.0	/* Time step in seconds */
#define secs_pd	86400.	/* Seconds per day */
#define dx	1.e5	/* Grid spacing in x direction */
#define dy	1.e5	/* Grid spacing in y direction */
#define	alpha	0.001	/* Asselin time filter parameter */
#define itmax	1000	/* Number of time steps in run */
#define mprint	50	/* Print diagnostics every mprint steps */

#define lower	1	/* low bound of range of processors */
#define upper	20	/* hi bound of range of processors */
#define version	0.6	/* version number of program */

/* constants used by dump() to determine what data structure to print */
#define one_dim		1
#define two_dim		2
#define p_label		0
#define u_label		1
#define v_label		2
#define pold_label	3
#define uold_label	4
#define vold_label	5
#define psi_label	6
#define cu_label	7
#define cv_label	8
#define h_label		9
#define z_label		10
#define dudt_label	11
#define dvdt_label	12

/* message types */
#define START_SIGNAL	0
#define END_SIGNAL	4
#define CALC1a		20
#define CALC1b		21
#define CALC1c		22
#define CALC2a		23
#define CALC2b		24
#define	TIME1a		30
#define	TIME1b		31
#define	TIME1c		32
#define	TIME1d		33
#define TIME2		34
#define P_ROW		50
#define U_ROW		51
#define V_ROW		52
#define PSI_ROW		53
#define POLD_ROW	54
#define UOLD_ROW	55
#define VOLD_ROW	56
#define H_ROW		57
#define Z_ROW		58

#define	PREV	0
#define	NEXT	1
#define	JSTART	2
#define	JEND	3

#define debug_data      0x1
#define debug_master    0x2
#define debug_worker    0x4
#define debug_call      0x8
#define debug           0x0

struct res
{
	float	row[m];
	int	indx;
};
