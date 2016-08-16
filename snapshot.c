
#include "snapshot.h"

#define COPY(type) int copy_ ## type (type * a, type * b) {    \
	if (a==NULL || b==NULL) return 255; \
	memcpy(a,b,sizeof(type));           \
	return 0;                           \
}
COPY(atom);
COPY(Box3);
const char s_timestep[] = "ITEM: TIMESTEP";
const char s_n_atoms[] = "ITEM: NUMBER OF ATOMS";
const char s_box_bounds[] = "ITEM: BOX BOUNDS pp pp pp";
const char s_box_bounds_z[] = "ITEM: BOX BOUNDS pp pp ff";
const char s_box_bounds_xy[] = "ITEM: BOX BOUNDS ff ff pp";
const char s_box_bounds_xyz[] = "ITEM: BOX BOUNDS ff ff ff";
const char s_atoms[]    = "ITEM: ATOMS id type xu yu zu mux muy muz";
const char s_atoms_pos[]    = "ITEM: ATOMS id type xu yu zu";
const char delimeter[] = " ";
int read_dump_OnlyCheck( FILE* fp) {
	const int i_timestep = strlen(s_timestep);
	const int i_n_atoms = strlen(s_n_atoms);
	const int i_box_bounds = strlen(s_box_bounds);
	const int i_atoms = strlen(s_atoms);
	const int i_atoms_pos = strlen(s_atoms_pos);

	bigint timestep;
	int n_atoms;
	real xlow,xhigh,ylow,yhigh,zlow,zhigh;
	int id,type;
	real xu,yu,zu,mux,muy,muz;
	int i;
	atom* p_atom;
#define FAIL 0
#define SUCCESS n_atoms
	read_lines(1,fp);
	if( strncmp(s_timestep,line,i_timestep) !=0) {
		return FAIL;
	}

	read_lines(1,fp);
	timestep = atol(line);
	fprintf(stderr,"atol(line) = %ld\n"
			"atoi(line) = %d\n"
			, timestep,atoi(line));

	read_lines(1,fp);
	if( strncmp(s_n_atoms,line,i_n_atoms) !=0)
		return FAIL;

	read_lines(1,fp);
	n_atoms = atoi(line);
	int  pbc[3]={0,0,0};
	read_lines(1,fp);
	if( strncmp(s_box_bounds,line,i_box_bounds) ==0) {
		pbc[0] = 1; pbc[1]=1; pbc[2]=1;
	}
	else if( strncmp(s_box_bounds_z,line,i_box_bounds) ==0) {
		pbc[0] = 1; pbc[1]=1; pbc[2]=0;
	}
	else if( strncmp(s_box_bounds_xy,line,i_box_bounds) ==0) {
		pbc[0] = 0; pbc[1]=0; pbc[2]=1;
	}
	else if( strncmp(s_box_bounds_xyz,line,i_box_bounds) ==0) {
		pbc[0] = 0; pbc[1]=0; pbc[2]=0;
	}
	else
		return FAIL;

	read_lines(1,fp); 
	xlow = atof(strtok(line,delimeter));
	xhigh = atof(strtok(NULL,delimeter));
	read_lines(1,fp); 
	ylow = atof(strtok(line,delimeter));
	yhigh = atof(strtok(NULL,delimeter));
	read_lines(1,fp); 
	zlow = atof(strtok(line,delimeter));
	zhigh = atof(strtok(NULL,delimeter));

	
	Snapshot *snap   = new_Snapshot(timestep,n_atoms);
	Box3 box = {xlow,xhigh,ylow,yhigh,zlow,zhigh, {pbc[0],pbc[1],pbc[2]}};
	copy_Box3( &snap->box,  &box);

	error( (char*)s_timestep);
//	fprintf(stderr,"%ld\n", timestep);
	error((char*)s_n_atoms);
//	fprintf(stderr,"%d\n", snap->n_atoms);
	error((char*)s_box_bounds);
//	fprintf(stderr,"%f %f\n", snap->box.xlow,snap->box.xhigh);
//	fprintf(stderr,"%f %f\n", snap->box.ylow,snap->box.yhigh);
//	fprintf(stderr,"%f %f\n", snap->box.zlow,snap->box.zhigh);

	read_lines(1,fp);

	if( strncmp(s_atoms,line,i_atoms) ==0) {
		for (i=0; i<n_atoms; i++){
			read_lines(1,fp);
			/*	fprintf(stderr,"%d %d %f %f %f %f %f %f\n", 
					id,type,
					xu,yu,zu,
					mux,muy,muz);*/
		}
	}
	else if (strncmp(s_atoms_pos,line,i_atoms_pos) ==0 ) {
		for (i=0; i<n_atoms; i++){
			read_lines(1,fp);
			/*	fprintf(stderr,"%d %d %f %f %f %f %f %f\n", 
					id,type,
					xu,yu,zu,
					mux,muy,muz);*/
		}

	}
	else  {
		return FAIL;
	}

	error_code =0;
	return SUCCESS;
#undef FAIL
#undef SUCCESS
}
Snapshot* read_dump( FILE* fp) {
	const char s_timestep[] = "ITEM: TIMESTEP";
	const char s_n_atoms[] = "ITEM: NUMBER OF ATOMS";
	const char s_box_bounds[] = "ITEM: BOX BOUNDS pp pp pp";
	const char s_box_bounds_z[] = "ITEM: BOX BOUNDS pp pp ff";
	const char s_box_bounds_xy[] = "ITEM: BOX BOUNDS ff ff pp";
	const char s_box_bounds_xyz[] = "ITEM: BOX BOUNDS ff ff ff";
	const char s_atoms[]    = "ITEM: ATOMS id type xu yu zu mux muy muz";
	const char s_atoms_pos[]    = "ITEM: ATOMS id type xu yu zu";
	const int i_timestep = strlen(s_timestep);
	const int i_n_atoms = strlen(s_n_atoms);
	const int i_box_bounds = strlen(s_box_bounds);
	const int i_atoms = strlen(s_atoms);
	const int i_atoms_pos = strlen(s_atoms_pos);
	const char delimeter[] = " ";

	bigint timestep;
	int n_atoms;
	real xlow,xhigh,ylow,yhigh,zlow,zhigh;
	int id,type;
	real xu,yu,zu,mux,muy,muz;
	int i;
	atom* p_atom;

	read_lines(1,fp);
	if( strncmp(s_timestep,line,i_timestep) !=0) {
		return (Snapshot*)(error("not ITEM: TIMESTEP"));
	}

	read_lines(1,fp);
	timestep = atol(line);
	fprintf(stderr,"atol(line) = %ld\n"
			"atoi(line) = %d\n"
			, timestep,atoi(line));

	read_lines(1,fp);
	if( strncmp(s_n_atoms,line,i_n_atoms) !=0)
		return (Snapshot*)(error("not ITEM: NUMBER OF ATOMS"));

	read_lines(1,fp);
	n_atoms = atoi(line);
	int  pbc[3]={0,0,0};
	read_lines(1,fp);
	if( strncmp(s_box_bounds,line,i_box_bounds) ==0) {
		pbc[0] = 1; pbc[1]=1; pbc[2]=1;
	}
	else if( strncmp(s_box_bounds_z,line,i_box_bounds) ==0) {
		pbc[0] = 1; pbc[1]=1; pbc[2]=0;
	}
	else if( strncmp(s_box_bounds_xy,line,i_box_bounds) ==0) {
		pbc[0] = 0; pbc[1]=0; pbc[2]=1;
	}
	else if( strncmp(s_box_bounds_xyz,line,i_box_bounds) ==0) {
		pbc[0] = 0; pbc[1]=0; pbc[2]=0;
	}
	else
		return (Snapshot*)(error("not ITEM: BOX BOUNDS pp pp pp(ff)"));

	read_lines(1,fp); 
	xlow = atof(strtok(line,delimeter));
	xhigh = atof(strtok(NULL,delimeter));
	read_lines(1,fp); 
	ylow = atof(strtok(line,delimeter));
	yhigh = atof(strtok(NULL,delimeter));
	read_lines(1,fp); 
	zlow = atof(strtok(line,delimeter));
	zhigh = atof(strtok(NULL,delimeter));

	
	Snapshot *snap   = new_Snapshot(timestep,n_atoms);
	Box3 box = {xlow,xhigh,ylow,yhigh,zlow,zhigh, {pbc[0],pbc[1],pbc[2]}};
	copy_Box3( &snap->box,  &box);

	error( (char*)s_timestep);
//	fprintf(stderr,"%ld\n", timestep);
	error((char*)s_n_atoms);
//	fprintf(stderr,"%d\n", snap->n_atoms);
	error((char*)s_box_bounds);
//	fprintf(stderr,"%f %f\n", snap->box.xlow,snap->box.xhigh);
//	fprintf(stderr,"%f %f\n", snap->box.ylow,snap->box.yhigh);
//	fprintf(stderr,"%f %f\n", snap->box.zlow,snap->box.zhigh);

	read_lines(1,fp);

	if( strncmp(s_atoms,line,i_atoms) ==0) {
		for (i=0; i<n_atoms; i++){
			read_lines(1,fp);
			id =atoi(strtok(line, delimeter));
			type = atoi(strtok(NULL,delimeter));
			xu   = atof(strtok(NULL,delimeter));
			yu   = atof(strtok(NULL,delimeter));
			zu   = atof(strtok(NULL,delimeter));
			mux   = atof(strtok(NULL,delimeter));
			muy   = atof(strtok(NULL,delimeter));
			muz   = atof(strtok(NULL,delimeter));

			p_atom = &(snap->atoms[i]);
			make_atom( p_atom,id,type,xu,yu,zu,mux,muy,muz);
			/*	fprintf(stderr,"%d %d %f %f %f %f %f %f\n", 
					id,type,
					xu,yu,zu,
					mux,muy,muz);*/
		}
	}
	else if (strncmp(s_atoms_pos,line,i_atoms_pos) ==0 ) {
		for (i=0; i<n_atoms; i++){
			read_lines(1,fp);
			id =atoi(strtok(line, delimeter));
			type = atoi(strtok(NULL,delimeter));
			xu   = atof(strtok(NULL,delimeter));
			yu   = atof(strtok(NULL,delimeter));
			zu   = atof(strtok(NULL,delimeter));
			mux   = 0.;
			muy   = 0.;
			muz   = 0.;

			p_atom = &(snap->atoms[i]);
			make_atom( p_atom,id,type,xu,yu,zu,mux,muy,muz);
			/*	fprintf(stderr,"%d %d %f %f %f %f %f %f\n", 
					id,type,
					xu,yu,zu,
					mux,muy,muz);*/
		}

	}
	else  {
		return (Snapshot*)(error("not ITEM: ATOMS id type xu yu zu mux muy muz"));
	}

	error_code =0;
	return snap;
}
void* error( char string[MAXLINE] ) {
	fputs( string, stderr );
	fputs( "\n", stderr );
	error_code =1;
	return NULL;
	//	exit(1);
}
void read_lines(int n,FILE* fp)  // from lammps reader_native.cpp
{
	char *eof;
	for (int i = 0; i < n; i++) eof = fgets(line,MAXLINE,fp);
	if (eof == NULL) error("Unexpected end of dump file");
}

Snapshot* new_Snapshot(bigint timestep, int n) {
	Snapshot* snap = (Snapshot*) malloc(sizeof(Snapshot));
	snap->timestep = timestep;
	snap-> n_atoms = n;
	snap->atoms = (atom*) malloc(sizeof(atom)*n);
	return snap;
}


void free_Snapshot(Snapshot* snap) {
 if(snap->atoms !=NULL)
	 free(snap->atoms);
 if(snap !=NULL)
	 free(snap);
}

int make_atom(atom* col,int id, int type, 
		real x, real y, real z,
		real mux,real muy, real muz) {
	real mu1;
	if (col ==NULL)
		return 255;
	col->id=id; col->type=type;
	col->x=x;col->y=y;col->z=z;
	mu1 = sqrt(mux*mux+muy*muy+muz*muz);
	col->mu1 = mu1;

	if (mu1>0.001) { 
		col->mux=mux/mu1;col->muy=muy/mu1;col->muz=muz/mu1;
	}
	return 0;
}
int dump_stream(atomstream* stream, FILE* fp, int nTime, int n_atoms,int s_id, int s_type) 
{
	int id,type;
	real xu,yu,zu,mux,muy,muz;
	atom* p_atom;
	const int i_atoms = strlen(s_atoms);
	const int i_atoms_pos = strlen(s_atoms_pos);
	fseek(fp, 0, SEEK_SET); // 

	
	for( int i=0; i<nTime; i++) {
		read_lines(9,fp);

		if( strncmp(s_atoms,line,i_atoms) ==0) {
			for (i=0; i<n_atoms; i++){
				read_lines(1,fp);
				id =atoi(strtok(line, delimeter));
				type = atoi(strtok(NULL,delimeter));
				if ( s_id==id && s_type == type) {
					xu   = atof(strtok(NULL,delimeter));
					yu   = atof(strtok(NULL,delimeter));
					zu   = atof(strtok(NULL,delimeter));
					mux   = atof(strtok(NULL,delimeter));
					muy   = atof(strtok(NULL,delimeter));
					muz   = atof(strtok(NULL,delimeter));

					stream->x[i] = xu;
					stream->y[i] = yu;
					stream->z[i] = zu;
					stream->mux[i] = mux;
					stream->muy[i] = muy;
					stream->muz[i] = muz;
				}

			}
		}
		else if (strncmp(s_atoms_pos,line,i_atoms_pos) ==0 ) {
			for (i=0; i<n_atoms; i++){
				read_lines(1,fp);
				id =atoi(strtok(line, delimeter));
				type = atoi(strtok(NULL,delimeter));
				if ( s_id==id && s_type == type) {
					xu   = atof(strtok(NULL,delimeter));
					yu   = atof(strtok(NULL,delimeter));
					zu   = atof(strtok(NULL,delimeter));
					stream->x[i] = xu;
					stream->y[i] = yu;
					stream->z[i] = zu;
				}

			}

		}

	}
	return 0;

}
int malloc_stream( atomstream* stream, int nTime) 
{
	stream->x = (real*) malloc( nTime* sizeof(real));
	stream->y = (real*) malloc( nTime* sizeof(real));
	stream->z = (real*) malloc( nTime* sizeof(real));
	stream->mux = (real*) malloc( nTime* sizeof(real));
	stream->muy = (real*) malloc( nTime* sizeof(real));
	stream->muz = (real*) malloc( nTime* sizeof(real));
}
int free_stream( atomstream* stream)
{
	free(&(stream->x));
	free(&(stream->y));
	free(&(stream->z));
	free(&(stream->mux));
	free(&(stream->muy));
	free(&(stream->muz));
}
