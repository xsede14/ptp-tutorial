
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

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include "decs.h"

MPI_Datatype *	setup_res();

extern void tstep(
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
float tdt);

void
worker()
{
	int	firststep, ncycle;
	float	tdt, time;
	int	i,j,ip,jp,jstart,jend;
	int	prv;
	int	nxt;
	
	int	msg_type;
	int	master_id, my_id;
	int	nprocs;
	int	nbytes;

	float 	p[n][m];	/* Pressure (or free surface height) */
	float 	u[n][m];	/* Zonal wind */
	float 	v[n][m];	/* Meridional wind */
	float 	psi[n][m];	/* Velocity streamfunction */
	float 	pold[n][m];
	float 	uold[n][m];
	float 	vold[n][m];
	float	pnew[n][m];
	float	unew[n][m];
	float	vnew[n][m];
	float	dpdt[n][m];
	float	dudt[n][m];	/* Time tendency of u */
	float	dvdt[n][m];
	float	cu[n][m];	/* Mass weighted u */
	float	cv[n][m];	/* Mass weighted v */
	float	h[n][m];       
	float	z[n][m];	/* Potential enstrophy */
	float	dummy1[m];
	float	dummy2[n];
	float	fsdx = 4./dx;
	float	fsdy = 4./dy;

	int		worker[4];
	float		p_start[m];
	float		u_start[m];
	float		v_start[m];
	float		psi_start[m];
	float		pold_start[m];
	float		uold_start[m];
	float		vold_start[m];
	MPI_Datatype *	res_type;
	MPI_Status	status;

	/*
	initialise control variables
	*/

	firststep = 1;
	ncycle = 0;
	tdt = dt;
	time = 0.;

	/*
	set up environment for worker
	*/

	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_id);

	MPI_Recv(&worker, 4, MPI_INT, MPI_ANY_SOURCE, START_SIGNAL,
		MPI_COMM_WORLD, &status);

	prv = worker[PREV];
	nxt = worker[NEXT];
	jstart = worker[JSTART];
	jend = worker[JEND];

	/*
	receive initialisation packets from master
	*/

	for (i = 0; i < n; i++){
		MPI_Recv(&p_start, m, MPI_FLOAT, 0, P_ROW, 
			MPI_COMM_WORLD, &status);
		acopy_one_to_two(p_start, p, i);

		MPI_Recv(&u_start, m, MPI_FLOAT, 0, U_ROW, 
			MPI_COMM_WORLD, &status);
		acopy_one_to_two(u_start, u, i);

		MPI_Recv(&v_start, m, MPI_FLOAT, 0, V_ROW, 
			MPI_COMM_WORLD, &status);
		acopy_one_to_two(v_start, v, i);

		MPI_Recv(&psi_start, m, MPI_FLOAT, 0, PSI_ROW, 
			MPI_COMM_WORLD, &status);
		acopy_one_to_two(psi_start, psi, i);

		MPI_Recv(&pold_start, m, MPI_FLOAT, 0, POLD_ROW, 
			MPI_COMM_WORLD, &status);
		acopy_one_to_two(pold_start, pold, i);

		MPI_Recv(&uold_start, m, MPI_FLOAT, 0, UOLD_ROW, 
			MPI_COMM_WORLD, &status);
		acopy_one_to_two(uold_start, uold, i);

		MPI_Recv(&vold_start, m, MPI_FLOAT, 0, VOLD_ROW, 
			MPI_COMM_WORLD, &status);
		acopy_one_to_two(vold_start, vold, i);
	}

	while (ncycle < itmax) {
		/*
		loop over latitudes calculating U, V, z and h
		do the block of latitudes from jstart to jend inclusive
		*/

		calc_load(prv, nxt, my_id, jstart, jend, p, u,v);
		calcuvzh(jstart, jend, p, u, v, cu, cv, h, z, fsdx, fsdy);
		calc_unload(prv, nxt, my_id, jstart, jend, cv, z);


		/*
		Calculate time tendencies of p, u and v
		*/

		time_load(prv, nxt, my_id, jstart, jend, cu, cv, h, z);
		timetend(jstart, jend, dpdt, dudt, dvdt, cu, cv, h, z);
		time_unload(prv, nxt, my_id, jstart, jend, dvdt);

		if ((my_id == 1) && (ncycle%mprint==0)) {
			diag(ncycle, time, p, u, v, h, z);
		}
		
		time += dt;

		tstep(m, n, alpha,
			jstart, jend, pold, uold, vold, p, u, v, pnew,
			unew, vnew, dpdt, dudt, dvdt, firststep, tdt);

		if ( firststep ) {
			/* Double tdt because all future steps are leapfrog */
			firststep = 0;
			tdt = tdt+tdt;
		}
		
		ncycle++;
	}  /* End of time step loop */

	/*
	send local data structures (results) back to master
	*/

	res_type = setup_res();

	send_updated_ds(res_type, jstart, jend, p, P_ROW, 0);
	send_updated_ds(res_type, jstart, jend, u, U_ROW, 0);
	send_updated_ds(res_type, jstart, jend, v, V_ROW, 0);
	send_updated_ds(res_type, jstart, jend, h, H_ROW, 0);
	send_updated_ds(res_type, jstart, jend, z, Z_ROW, 0);

	MPI_Recv(&worker, 4, MPI_INT, 0, END_SIGNAL,
		MPI_COMM_WORLD, &status);

	if (debug & debug_call) {
		printf("worker %d sent TIDY_UP to master\n", my_id);
		printf("worker %d got END_SIGNAL from master\n", my_id);
	}
}

send_updated_ds(res_type, jstart, jend, ds, indx, master_id)
	MPI_Datatype *	res_type;
	int		jstart;
	int		jend;
	float		ds[n][m];
	int		indx;
	int		master_id;
{
	int		j;
	struct res	res;
	MPI_Request	rq[2];
	MPI_Status	stat[2];

	for (j = jstart; j <= jend; j++) {
		acopy_two_to_one(ds, res.row, j);
		res.indx = j;

		MPI_Send(&res, 1, *res_type, master_id, indx,
			MPI_COMM_WORLD);
	}
}

/*
this procedure does all the message passing before the call to _calcuvzh_
*/
calc_load(prv,nxt,my_id,jstart,jend,p,u,v)
	int	prv;
	int	nxt;
	int	my_id;
	int	jstart;
	int	jend;
	float	p[n][m];
	float	u[n][m];
	float	v[n][m];
{
	neighbour_send(prv, my_id, CALC1a, p, jstart);
	neighbour_send(prv, my_id, CALC1b, u, jstart);
	neighbour_send(prv, my_id, CALC1c, v, jstart);
	neighbour_receive(nxt, my_id, CALC1a, p, (jend+1) % n);
	neighbour_receive(nxt, my_id, CALC1b, u, (jend+1) % n);
	neighbour_receive(nxt, my_id, CALC1c, v, (jend+1) % n);
}

/*
this procedure does all the message passing after the call to _calcuvzh_
*/
calc_unload(prv,nxt,my_id,jstart,jend,cv,z)
	int	prv;
	int	nxt;
	int	my_id;
	int	jstart;
	int	jend;
	float	cv[n][m];
	float	z[n][m];
{
	neighbour_send(nxt, my_id, CALC2a, cv, (jend+1) % n);
	neighbour_send(nxt, my_id, CALC2b, z, (jend+1) % n);
	neighbour_receive(prv, my_id, CALC2a, cv, jstart);
	neighbour_receive(prv, my_id, CALC2b, z, jstart);
}

/*
this procedure does all the message passing before the call to _timetend_
*/
time_load(prv,nxt,my_id,jstart,jend,cu,cv,h,z)
	int	prv;
	int	nxt;
	int	my_id;
	int	jstart;
	int	jend;
	float	cu[n][m];
	float	cv[n][m];
	float	h[n][m];
	float	z[n][m];
{
	neighbour_send(prv, my_id, TIME1a, cu, jstart);
	neighbour_send(prv, my_id, TIME1b, cv, jstart);
	neighbour_send(prv, my_id, TIME1c, h, jstart);
	neighbour_send(prv, my_id, TIME1d, z, jstart);
	neighbour_receive(nxt, my_id, TIME1a, cu, (jend+1) % n);
	neighbour_receive(nxt, my_id, TIME1b, cv, (jend+1) % n);
	neighbour_receive(nxt, my_id, TIME1c, h, (jend+1) % n);
	neighbour_receive(nxt, my_id, TIME1d, z, (jend+1) % n);
}

/*
this procedure does all the message passing after the call to _timetend_
*/
time_unload(prv,nxt,tu_my_id,jstart,jend,dvdt)
	int	prv;
	int	nxt;
	int	tu_my_id;
	int	jstart;
	int	jend;
	float	dvdt[n][m];
{
	neighbour_send(nxt, tu_my_id, TIME2, dvdt, (jend+1) % n);
	neighbour_receive(prv, tu_my_id, TIME2, dvdt, jstart);
}

/*
this is a general purpose function for sending packets b/w workers
*/
neighbour_send(ns_neighbour,ns_my_id,ns_rec_id,ns_ds,ns_edge)
	int	ns_neighbour;
	int ns_my_id;
	int	ns_rec_id;
	float	ns_ds[n][m];
	int	ns_edge;
{
	float		ns_rec[m];
	MPI_Request	rq;

	acopy_two_to_one(ns_ds, ns_rec, ns_edge);

	MPI_Isend(&ns_rec, m, MPI_FLOAT, ns_neighbour, ns_rec_id,
		MPI_COMM_WORLD, &rq);

	if (debug & debug_worker)
		printf("worker %d sent packet %d to worker %d\n", ns_my_id,
			ns_rec_id, ns_neighbour);
}

/*
this is a general purpose function for receiving packets b/w workers
*/
neighbour_receive(nr_neighbour,nr_my_id,nr_rec_id,nr_ds,nr_edge)
	int	nr_neighbour;
	int	nr_my_id;
	int	nr_rec_id;
	float	nr_ds[n][m];
	int	nr_edge;
{
	float		nr_rec[m];
	MPI_Status	status;

	MPI_Recv(&nr_rec, m, MPI_FLOAT, nr_neighbour, nr_rec_id,
		MPI_COMM_WORLD, &status);

	acopy_one_to_two(nr_rec, nr_ds, nr_edge);
}
