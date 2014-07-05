
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

extern void worker();
MPI_Datatype *	setup_res();

main (argc, argv)
	int	argc;
	char *	argv[];
{
	float	pi=4.*(float)atan((double)1.);
	float 	p[n][m];	/* Pressure (or free surface height) */
	float 	u[n][m];	/* Zonal wind */
	float 	v[n][m];	/* Meridional wind */
	float 	psi[n][m];	/* Velocity streamfunction */
	float 	pold[n][m];
	float 	uold[n][m];
	float 	vold[n][m];
	float	h[n][m];
	float	z[n][m];
	float	dummy1[m];
	float	dummy2[n][m];
	float	tpi=pi+pi;
	float	di=tpi/(float)m;
	float	dj=tpi/(float)n;
	int	i, j, chunk_size, nxt, prv;

	int	master_packet[4];
	float	p_start[m];
	float	u_start[m];
	float	v_start[m];
	float	psi_start[m];
	float	pold_start[m];
	float	uold_start[m];
	float	vold_start[m];
	int	proc_cnt;
	int	tid;
	MPI_Datatype *	res_type;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &proc_cnt);
	MPI_Comm_rank(MPI_COMM_WORLD, &tid);

	if ( proc_cnt < 2 )
	{
		fprintf(stderr, "must have at least 2 processes, not %d\n", proc_cnt);
		MPI_Finalize();
		return 1;
	}

	if ( (n % (proc_cnt - 1)) != 0 )
	{
		if ( tid == 0 )
			fprintf(stderr, "(number of processes - 1) must be a multiple of %d\n", n);

		MPI_Finalize();
		return 1;
	}

	if (tid != 0) {
		worker();
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
	} else {

	/* master process */

	chunk_size = n / (proc_cnt - 1);

	for (i = 1; i < proc_cnt; i++) {
		/* calculate each worker's boundary */
		master_packet[JSTART] = (i - 1) * chunk_size;

		if (i == proc_cnt - 1)
			master_packet[JEND] = n - 1;
		else
			master_packet[JEND] = i * chunk_size - 1;

		if (i == 1)
			prv = proc_cnt-1;
		else
			prv = i-1;

		master_packet[PREV] = prv;

		if (i == proc_cnt - 1)
			nxt = 1;
		else
			nxt = i+1;

		master_packet[NEXT] = nxt;

		MPI_Send(&master_packet, 4, MPI_INT, i, START_SIGNAL,
			MPI_COMM_WORLD);

	printf("jstart=%d, jend=%d, next=%d, prev=%d\n", 
		master_packet[JSTART],
		master_packet[JEND],
		master_packet[NEXT],
		master_packet[PREV]);
	}


	/*
	initialise data structures and construct packets to be sent to workers
	*/

	initialise(p, u, v, psi, pold, uold, vold, di, dj, z);
	diag(1, 0., p, u, v, h, z);

	for (i = 1; i < proc_cnt; i++) {
		for (j = 0; j < n; j++) {
			acopy_two_to_one(p, p_start, j);
			MPI_Send(&p_start, m, MPI_FLOAT, i, P_ROW, 
				MPI_COMM_WORLD);

			acopy_two_to_one(u, u_start, j);
			MPI_Send(&u_start, m, MPI_FLOAT, i, U_ROW, 
				MPI_COMM_WORLD);

			acopy_two_to_one(v, v_start, j);
			MPI_Send(&v_start, m, MPI_FLOAT, i, V_ROW, 
				MPI_COMM_WORLD);

			acopy_two_to_one(psi, psi_start, j);
			MPI_Send(&psi_start, m, MPI_FLOAT, i, PSI_ROW, 
				MPI_COMM_WORLD);

			acopy_two_to_one(pold, pold_start, j);
			MPI_Send(&pold_start, m, MPI_FLOAT, i, POLD_ROW, 
				MPI_COMM_WORLD);

			acopy_two_to_one(uold, uold_start, j);
			MPI_Send(&uold_start, m, MPI_FLOAT, i, UOLD_ROW, 
				MPI_COMM_WORLD);

			acopy_two_to_one(vold, vold_start, j);
			MPI_Send(&vold_start, m, MPI_FLOAT, i, VOLD_ROW, 
				MPI_COMM_WORLD);
		}
	}

	/*
	receive packets back from the workers
	*/
	res_type = setup_res();

	if ( debug & debug_master )
		printf("receiving P\n");

	update_global_ds(res_type, P_ROW, p);

	if ( debug & debug_master )
		printf("receiving U\n");

	update_global_ds(res_type, U_ROW, u);

	if ( debug & debug_master )
		printf("receiving V\n");

	update_global_ds(res_type, V_ROW, v);

	if ( debug & debug_master )
		printf("receiving H\n");

	update_global_ds(res_type, H_ROW, h);

	if ( debug & debug_master )
		printf("receiving Z\n");

	update_global_ds(res_type, Z_ROW, z);

	for (i = 1; i < proc_cnt; i++){
		MPI_Send(&master_packet, 4, MPI_INT, i, END_SIGNAL,
			MPI_COMM_WORLD);
	}

	/* wait for all workers to end */
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();

	diag(itmax, itmax*dt, p, u, v, h, z);
	}
	
	return 0;
}

MPI_Datatype *
setup_res()
{
	struct res		res;
	MPI_Aint		res_disp[2];
	static int		res_done = 0;
	static int		res_len[2] = { m, 1 };
	static MPI_Datatype	res_old[2] = { MPI_FLOAT, MPI_INT };
	static MPI_Datatype	res_type;

	if ( res_done )
		return &res_type;

	res_done++;
	MPI_Address(&res.row[0], &res_disp[0]);
	MPI_Address(&res.indx, &res_disp[1]);
	res_disp[1] -= res_disp[0];
	res_disp[0] = 0;
	MPI_Type_struct(2, res_len, res_disp, res_old, &res_type);
	MPI_Type_commit(&res_type);

	return &res_type;
}

/*
this function waits for all the workers to return the packets of
a particular type and then updates the master's copy of the same type
*/
update_global_ds(res_type, indx, ds)
	MPI_Datatype *	res_type;
	int		indx;
	float		ds[n][m];

{
	int		i;
	int		row;
	struct res	res;
	MPI_Status	status;

	for (i = 0; i < n; i++) {
		MPI_Recv(&res, 1, *res_type, MPI_ANY_SOURCE, indx,
			MPI_COMM_WORLD, &status);

		acopy_one_to_two(res.row, ds, res.indx);
	}

}
