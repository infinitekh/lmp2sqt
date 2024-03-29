
#include "snapshot.h"
#include <stdbool.h>

#define COPY(type) int copy_ ## type (type * a, type * b) {    \
	if (a==NULL || b==NULL) return 255; \
	memcpy(a,b,sizeof(type));           \
	return 0;                           \
}
COPY(Atom);
COPY(Box3);
const char s_timestep[] = "ITEM: TIMESTEP";
const char s_n_atoms[] = "ITEM: NUMBER OF ATOMS";
const char s_box_bounds[] = "ITEM: BOX BOUNDS pp pp pp";
const char s_box_bounds_z[] = "ITEM: BOX BOUNDS pp pp ff";
const char s_box_bounds_xy[] = "ITEM: BOX BOUNDS ff ff pp";
const char s_box_bounds_xyz[] = "ITEM: BOX BOUNDS ff ff ff";
const char s_atoms_dipole[]    = "ITEM: ATOMS id type xu yu zu mux muy muz";
const char s_atoms_vel[]    = "ITEM: ATOMS id type xu yu zu vx vy vz";
const char s_atoms_all[]    = "ITEM: ATOMS id type xu yu zu mux muy muz vx vy vz";
const char s_atoms[]    = "ITEM: ATOMS id type xu yu zu";
const char delimeter[] = " ";
char * strtok_null_safe ( char * str, const char * delimiters );

int ReadDumpForOnlyCheck(FILE *fp) {
  const int i_timestep = strlen(s_timestep);
  const int i_n_atoms = strlen(s_n_atoms);
  const int i_box_bounds = strlen(s_box_bounds);
  const int i_atoms = strlen(s_atoms);
  const int i_atoms_vel = strlen(s_atoms_vel);
  const int i_atoms_dipole = strlen(s_atoms_dipole);
  const int i_atoms_all = strlen(s_atoms_all);

  bigint timestep;
  int n_atoms;
  real xlow, xhigh, ylow, yhigh, zlow, zhigh;
  error_code = 0;
  /* 	atom_column atom_ids = {
   * 		.id = -1 ,.type=-1,
   * 	 .x= -1,.y=-1,.z=-1,
   * 	 .mux=-1,.muy=-1,.muz=-1,
   * 	 .vx = -1,.vy=-1,.vz =-1
   * 	};
   */
  int id, type;
  real xu, yu, zu;
  real mux, muy, muz;
  real vx, vy, vz;
  int i;
  Atom *p_atom;
#define FAIL false
#define SUCCESS n_atoms
	read_lines(1,fp);
	if( strncmp(s_timestep,line,i_timestep) !=0) {
		return FAIL;
	}

	read_lines(1,fp);
	timestep = atol(line);
#ifndef NDEBUG
	fprintf(stderr,"atol(line) = %ld\n"
			"atoi(line) = %d\n"
			, timestep,atoi(line));
#endif

	read_lines(1,fp);
	if( strncmp(s_n_atoms,line,i_n_atoms) !=0)
		return FAIL;

	read_lines(1,fp);
	n_atoms = atoi(line);
//	int  pbc[3]={0,0,0};
	int  pbcTYPE=0;
	read_lines(1,fp);
	if( strncmp(s_box_bounds,line,i_box_bounds) ==0) {
		pbcTYPE = PBC_X|PBC_Y|PBC_Z;
//		pbc[0] = 1; pbc[1]=1; pbc[2]=1;
	}
	else if( strncmp(s_box_bounds_z,line,i_box_bounds) ==0) {
		pbcTYPE = PBC_X|PBC_Y;
//		pbc[0] = 1; pbc[1]=1; pbc[2]=0;
	}
	else if( strncmp(s_box_bounds_xy,line,i_box_bounds) ==0) {
		pbcTYPE = PBC_Z;
//		pbc[0] = 0; pbc[1]=0; pbc[2]=1;
	}
	else if( strncmp(s_box_bounds_xyz,line,i_box_bounds) ==0) {
		 pbcTYPE =0;
	}
	else
		return FAIL;

	read_lines(1,fp); 
	xlow = atof(strtok_null_safe(line,delimeter));
	xhigh = atof(strtok_null_safe(NULL,delimeter));
	read_lines(1,fp); 
	ylow = atof(strtok_null_safe(line,delimeter));
	yhigh = atof(strtok_null_safe(NULL,delimeter));
	read_lines(1,fp); 
	zlow = atof(strtok_null_safe(line,delimeter));
	zhigh = atof(strtok_null_safe(NULL,delimeter));

        Snapshot *snap = NewSnapshot(timestep, n_atoms);
        Box3 box = {xlow,xhigh,ylow,yhigh,zlow,zhigh, pbcTYPE};
//	Box3 box = {xlow,xhigh,ylow,yhigh,zlow,zhigh, {pbc[0],pbc[1],pbc[2]}};
        copy_Box3(&snap->Box, &box);

        //	error( (char*)s_timestep);
        //	fprintf(stderr,"%ld\n", timestep);
        //	error((char*)s_n_atoms);
        //	fprintf(stderr,"%d\n", snap->n_atoms);
        //	error((char*)s_box_bounds);
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
			if (error_code == 1) return FAIL;
		}
	}
	else if (strncmp(s_atoms_vel,line,i_atoms_vel) ==0 ) {
		for (i=0; i<n_atoms; i++){
			read_lines(1,fp);
			if (error_code == 1) return FAIL;
			/*	fprintf(stderr,"%d %d %f %f %f %f %f %f\n", 
					id,type,
					xu,yu,zu,
					mux,muy,muz);*/
		}
	}
	else if (strncmp(s_atoms_dipole,line,i_atoms_dipole) ==0 ) {
		for (i=0; i<n_atoms; i++){
			read_lines(1,fp);
			if (error_code == 1) return FAIL;
			/*	fprintf(stderr,"%d %d %f %f %f %f %f %f\n", 
					id,type,
					xu,yu,zu,
					mux,muy,muz);*/
		}
	}
	else if (strncmp(s_atoms_all,line,i_atoms_all) ==0 ) {
		for (i=0; i<n_atoms; i++){
			read_lines(1,fp);
			if (error_code == 1) return FAIL;
			/*	fprintf(stderr,"%d %d %f %f %f %f %f %f\n", 
					id,type,
					xu,yu,zu,
					mux,muy,muz);*/
		}
	}
	else  {
		return FAIL;
	}

	return SUCCESS;
}
Snapshot *ReadDump(FILE *fp) {
  const char s_timestep[] = "ITEM: TIMESTEP";
  const char s_n_atoms[] = "ITEM: NUMBER OF ATOMS";
  const char s_box_bounds[] = "ITEM: BOX BOUNDS pp pp pp";
  const char s_box_bounds_z[] = "ITEM: BOX BOUNDS pp pp ff";
  const char s_box_bounds_xy[] = "ITEM: BOX BOUNDS ff ff pp";
  const char s_box_bounds_xyz[] = "ITEM: BOX BOUNDS ff ff ff";
  const char s_atoms_dipole[] = "ITEM: ATOMS id type xu yu zu mux muy muz";
  const char s_atoms_vel[] = "ITEM: ATOMS id type xu yu zu vx vy vz";
  const char s_atoms_all[] =
      "ITEM: ATOMS id type xu yu zu mux muy muz vx vy vz";
  const char s_atoms[] = "ITEM: ATOMS id type xu yu zu";
  const int i_timestep = strlen(s_timestep);
  const int i_n_atoms = strlen(s_n_atoms);
  const int i_box_bounds = strlen(s_box_bounds);
  const int i_atoms = strlen(s_atoms);
  const int i_atoms_vel = strlen(s_atoms_vel);
  const int i_atoms_dipole = strlen(s_atoms_dipole);
  const int i_atoms_all = strlen(s_atoms_all);
  const char delimeter[] = " ";

  bigint timestep;
  int n_atoms;
  real xlow, xhigh, ylow, yhigh, zlow, zhigh;
  int id, type;
  real xu, yu, zu, mux, muy, muz;
  real vx, vy, vz;
  int i;
  Atom *p_atom;
  error_code = 0;

  read_lines(1, fp);
  if (strncmp(s_timestep, line, i_timestep) != 0) {
    return (Snapshot *)(error("not ITEM: TIMESTEP"));
  }

  read_lines(1, fp);
  timestep = atol(line);
#ifndef NDEBUG
	fprintf(stderr,"atol(line) = %ld\n"
			"atoi(line) = %d\n"
			, timestep,atoi(line));
#endif

	read_lines(1,fp);
	if( strncmp(s_n_atoms,line,i_n_atoms) !=0)
		return (Snapshot*)(error("not ITEM: NUMBER OF ATOMS"));

	read_lines(1,fp);
	n_atoms = atoi(line);
	int  pbcTYPE=0;
	read_lines(1,fp);
	if( strncmp(s_box_bounds,line,i_box_bounds) ==0) {
		pbcTYPE = PBC_X|PBC_Y|PBC_Z;
//		pbc[0] = 1; pbc[1]=1; pbc[2]=1;
	}
	else if( strncmp(s_box_bounds_z,line,i_box_bounds) ==0) {
		pbcTYPE = PBC_X|PBC_Y;
//		pbc[0] = 1; pbc[1]=1; pbc[2]=0;
	}
	else if( strncmp(s_box_bounds_xy,line,i_box_bounds) ==0) {
		pbcTYPE = PBC_Z;
//		pbc[0] = 0; pbc[1]=0; pbc[2]=1;
	}
	else if( strncmp(s_box_bounds_xyz,line,i_box_bounds) ==0) {
		 pbcTYPE =0;
	}
	else
		return (Snapshot*)(error("not ITEM: BOX BOUNDS pp pp pp(ff)"));

	read_lines(1,fp); 
	xlow = atof(strtok_null_safe(line,delimeter));
	xhigh = atof(strtok_null_safe(NULL,delimeter));
	read_lines(1,fp); 
	ylow = atof(strtok_null_safe(line,delimeter));
	yhigh = atof(strtok_null_safe(NULL,delimeter));
	read_lines(1,fp); 
	zlow = atof(strtok_null_safe(line,delimeter));
	zhigh = atof(strtok_null_safe(NULL,delimeter));

        Snapshot *snap = NewSnapshot(timestep, n_atoms);
        //	Box3 box = {xlow,xhigh,ylow,yhigh,zlow,zhigh,
        //{pbc[0],pbc[1],pbc[2]}};
	Box3 box = {xlow,xhigh,ylow,yhigh,zlow,zhigh, pbcTYPE};
        copy_Box3(&snap->Box, &box);

        error( (char*)s_timestep);
//	fprintf(stderr,"%ld\n", timestep);
	error((char*)s_n_atoms);
//	fprintf(stderr,"%d\n", snap->n_atoms);
	error((char*)s_box_bounds);
//	fprintf(stderr,"%f %f\n", snap->box.xlow,snap->box.xhigh);
//	fprintf(stderr,"%f %f\n", snap->box.ylow,snap->box.yhigh);
//	fprintf(stderr,"%f %f\n", snap->box.zlow,snap->box.zhigh);

	read_lines(1,fp);
	
	if( strncmp(s_atoms_dipole,line,i_atoms_dipole) ==0) {
#ifndef NDEBUG
		fprintf(stderr,"s_atoms_dipole style snapshot\n");
#endif
		for (i=0; i<n_atoms; i++){
			read_lines(1,fp);
			if (error_code == 1) return FAIL;
			id =atoi(strtok_null_safe(line, delimeter));
			type = atoi(strtok_null_safe(NULL,delimeter));
			xu   = atof(strtok_null_safe(NULL,delimeter));
			yu   = atof(strtok_null_safe(NULL,delimeter));
			zu   = atof(strtok_null_safe(NULL,delimeter));
			mux   = atof(strtok_null_safe(NULL,delimeter));
			muy   = atof(strtok_null_safe(NULL,delimeter));
			muz   = atof(strtok_null_safe(NULL,delimeter));

			p_atom = &(snap->atoms[i]);
			make_atom_dipole( p_atom,id,type,xu,yu,zu,mux,muy,muz);
		}
	}
	else if (strncmp(s_atoms_vel,line,i_atoms_vel) ==0 ) {
#ifndef NDEBUG
		fprintf(stderr,"s_atoms_vel style snapshot\n");
#endif
		for (i=0; i<n_atoms; i++){
			read_lines(1,fp);
			if (error_code == 1) return FAIL;
			id =atoi(strtok_null_safe(line, delimeter));
			type = atoi(strtok_null_safe(NULL,delimeter));
			xu   = atof(strtok_null_safe(NULL,delimeter));
			yu   = atof(strtok_null_safe(NULL,delimeter));
			zu   = atof(strtok_null_safe(NULL,delimeter));
			vx   = atof(strtok_null_safe(NULL,delimeter));
			vy   = atof(strtok_null_safe(NULL,delimeter));
			vz   = atof(strtok_null_safe(NULL,delimeter));

			p_atom = &(snap->atoms[i]);
			make_atom_vel( p_atom,id,type,xu,yu,zu,vx,vy,vz);
			/*	fprintf(stderr,"%d %d %f %f %f %f %f %f\n", 
					id,type,
					xu,yu,zu,
					mux,muy,muz);*/
		}

	}
	else if (strncmp(s_atoms_all,line,i_atoms_all) ==0 ) {
#ifndef NDEBUG
		fprintf(stderr,"s_atoms_all style snapshot\n");
#endif
		for (i=0; i<n_atoms; i++){
			read_lines(1,fp);
			if (error_code == 1) return FAIL;
			id =atoi(strtok_null_safe(line, delimeter));
			type = atoi(strtok_null_safe(NULL,delimeter));
			xu   = atof(strtok_null_safe(NULL,delimeter));
			yu   = atof(strtok_null_safe(NULL,delimeter));
			zu   = atof(strtok_null_safe(NULL,delimeter));
			mux   = atof(strtok_null_safe(NULL,delimeter));
			muy   = atof(strtok_null_safe(NULL,delimeter));
			muz   = atof(strtok_null_safe(NULL,delimeter));
			vx   = atof(strtok_null_safe(NULL,delimeter));
			vy   = atof(strtok_null_safe(NULL,delimeter));
			vz   = atof(strtok_null_safe(NULL,delimeter));

			p_atom = &(snap->atoms[i]);
			make_atom_all( p_atom,id,type,xu,yu,zu,mux,muy,muz,vx,vy,vz);
		}

	}
	else if (strncmp(s_atoms_vel,line,i_atoms_vel) ==0 ) {
#ifndef NDEBUG
		fprintf(stderr,"s_atoms_vel style snapshot\n");
#endif
		for (i=0; i<n_atoms; i++){
			read_lines(1,fp);
			if (error_code == 1) return FAIL;
			id =atoi(strtok_null_safe(line, delimeter));
			type = atoi(strtok_null_safe(NULL,delimeter));
			xu   = atof(strtok_null_safe(NULL,delimeter));
			yu   = atof(strtok_null_safe(NULL,delimeter));
			zu   = atof(strtok_null_safe(NULL,delimeter));
			vx   = atof(strtok_null_safe(NULL,delimeter));
			vy   = atof(strtok_null_safe(NULL,delimeter));
			vz   = atof(strtok_null_safe(NULL,delimeter));

			p_atom = &(snap->atoms[i]);
			make_atom_vel( p_atom,id,type,xu,yu,zu,vx,vy,vz);
			/*	fprintf(stderr,"%d %d %f %f %f %f %f %f\n", 
					id,type,
					xu,yu,zu,
					mux,muy,muz);*/
		}

	}
	else if (strncmp(s_atoms,line,i_atoms) ==0 ) {
		for (i=0; i<n_atoms; i++){
			read_lines(1,fp);
			if (error_code == 1) return FAIL;
			id =atoi(strtok_null_safe(line, delimeter));
			type = atoi(strtok_null_safe(NULL,delimeter));
			xu   = atof(strtok_null_safe(NULL,delimeter));
			yu   = atof(strtok_null_safe(NULL,delimeter));
			zu   = atof(strtok_null_safe(NULL,delimeter));

			p_atom = &(snap->atoms[i]);
			make_atom( p_atom,id,type,xu,yu,zu);
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
#ifndef NDEBUG
	fputs( string, stderr );
	fputs( "\n", stderr );
#endif
	return NULL;
	//	exit(1);
}
char * strtok_null_safe ( char * str, const char * delimiters ) {

	char* ret_str = strtok(str,delimiters);
	if (ret_str == NULL ) {
		error_code = 1 ; 
		return "0";
	}
	return  ret_str;
}

void read_lines(int n,FILE* fp)  // from lammps reader_native.cpp
{
	char *eof;int i;
	memset(line, 0, sizeof(line));                                           

	for (i = 0; i < n; i++) eof = fgets(line,MAXLINE,fp);
	if (eof == NULL) {
		error_code = 1;
		error("Unexpected end of dump file");
	}
}

Snapshot *NewSnapshot(bigint timestep, int n) {
  Snapshot *snap = (Snapshot *)malloc(sizeof(Snapshot));
  snap->timestep = timestep;
  snap->NumAtoms = n;
  snap->atoms = (Atom *)malloc(sizeof(Atom) * n);
  return snap;
}

void FreeSnapshot(Snapshot *snap) {
  if (snap->atoms != NULL)
    free(snap->atoms);
  if (snap != NULL)
    free(snap);
}

int make_atom(Atom *col, int id, int type, real x, real y, real z) {
  real mu1;
  if (col == NULL)
    return 255;
  col->id = id;
  col->type = type;
  col->x = x;
  col->y = y;
  col->z = z;
  col->atomType = 0;
  return 0;
}
int make_atom_dipole(Atom *col, int id, int type, real x, real y, real z,
                     real mux, real muy, real muz) {
  real mu1;
  if (col == NULL)
    return 255;
  col->id = id;
  col->type = type;
  col->x = x;
  col->y = y;
  col->z = z;
  mu1 = sqrt(mux * mux + muy * muy + muz * muz);
  col->mu1 = mu1;

  if (mu1 > 0.001) {
    col->mux = mux / mu1;
    col->muy = muy / mu1;
    col->muz = muz / mu1;
  }
  col->atomType = ATOM_DIPOLE;
  return 0;
}
int make_atom_vel(Atom *col, int id, int type, real x, real y, real z, real vx,
                  real vy, real vz) {
  real mu1;
  if (col == NULL)
    return 255;
  col->id = id;
  col->type = type;
  col->x = x;
  col->y = y;
  col->z = z;
  col->vx = vx;
  col->vy = vy;
  col->vz = vz;
  col->atomType = ATOM_VEL;
  return 0;
}
int make_atom_all(Atom *col, int id, int type, real x, real y, real z, real mux,
                  real muy, real muz, real vx, real vy, real vz) {
  real mu1;
  if (col == NULL)
    return 255;
  col->id = id;
  col->type = type;
  col->x = x;
  col->y = y;
  col->z = z;
  mu1 = sqrt(mux * mux + muy * muy + muz * muz);
  col->mu1 = mu1;
  col->vx = vx;
  col->vy = vy;
  col->vz = vz;
  col->atomType = ATOM_ALL;

  if (mu1 > 0.001) {
    col->mux = mux / mu1;
    col->muy = muy / mu1;
    col->muz = muz / mu1;
  }
  return 0;
}
int DumpAtomStream(AtomStream *stream, FILE *fp, int nTime, int n_atoms,
                   int s_id, int s_type) {
  int id, type;
  real xu, yu, zu, mux, muy, muz;
  real vx, vy, vz;
  Atom *p_atom;
  const int i_atoms = strlen(s_atoms);
  const int i_atoms_vel = strlen(s_atoms_vel);
  const int i_atoms_dipole = strlen(s_atoms_dipole);
  const int i_atoms_all = strlen(s_atoms_all);
  fseek(fp, 0, SEEK_SET); //
  int i;

  for (i = 0; i < nTime; i++) {
    read_lines(9, fp);

    if (strncmp(s_atoms_dipole, line, i_atoms_dipole) == 0) {
      for (i = 0; i < n_atoms; i++) {
        read_lines(1, fp);
        if (error_code == 1)
          return FAIL;
        id = atoi(strtok_null_safe(line, delimeter));
        type = atoi(strtok_null_safe(NULL, delimeter));
        if (s_id == id && s_type == type) {
          xu = atof(strtok_null_safe(NULL, delimeter));
          yu = atof(strtok_null_safe(NULL, delimeter));
          zu = atof(strtok_null_safe(NULL, delimeter));
          mux = atof(strtok_null_safe(NULL, delimeter));
          muy = atof(strtok_null_safe(NULL, delimeter));
          muz = atof(strtok_null_safe(NULL, delimeter));

          stream->x[i] = xu;
          stream->y[i] = yu;
          stream->z[i] = zu;
          stream->mux[i] = mux;
          stream->muy[i] = muy;
          stream->muz[i] = muz;
          stream->atomType = ATOM_DIPOLE;
        }
      }
    } else if (strncmp(s_atoms_vel, line, i_atoms_vel) == 0) {
      for (i = 0; i < n_atoms; i++) {
        read_lines(1, fp);
        if (error_code == 1)
          return FAIL;
        id = atoi(strtok_null_safe(line, delimeter));
        type = atoi(strtok_null_safe(NULL, delimeter));
        if (s_id == id && s_type == type) {
          xu = atof(strtok_null_safe(NULL, delimeter));
          yu = atof(strtok_null_safe(NULL, delimeter));
          zu = atof(strtok_null_safe(NULL, delimeter));
          vx = atof(strtok_null_safe(NULL, delimeter));
          vy = atof(strtok_null_safe(NULL, delimeter));
          vz = atof(strtok_null_safe(NULL, delimeter));

          stream->x[i] = xu;
          stream->y[i] = yu;
          stream->z[i] = zu;
          stream->vx[i] = vx;
          stream->vy[i] = vy;
          stream->vz[i] = vz;
          stream->atomType = ATOM_VEL;
        }
      }
    } else if (strncmp(s_atoms_all, line, i_atoms_all) == 0) {
      for (i = 0; i < n_atoms; i++) {
        read_lines(1, fp);
        if (error_code == 1)
          return FAIL;
        id = atoi(strtok_null_safe(line, delimeter));
        type = atoi(strtok_null_safe(NULL, delimeter));
        if (s_id == id && s_type == type) {
          xu = atof(strtok_null_safe(NULL, delimeter));
          yu = atof(strtok_null_safe(NULL, delimeter));
          zu = atof(strtok_null_safe(NULL, delimeter));
          mux = atof(strtok_null_safe(NULL, delimeter));
          muy = atof(strtok_null_safe(NULL, delimeter));
          muz = atof(strtok_null_safe(NULL, delimeter));
          vx = atof(strtok_null_safe(NULL, delimeter));
          vy = atof(strtok_null_safe(NULL, delimeter));
          vz = atof(strtok_null_safe(NULL, delimeter));

          stream->x[i] = xu;
          stream->y[i] = yu;
          stream->z[i] = zu;
          stream->vx[i] = vx;
          stream->vy[i] = vy;
          stream->vz[i] = vz;
          stream->mux[i] = mux;
          stream->muy[i] = muy;
          stream->muz[i] = muz;
          stream->atomType = ATOM_ALL;
        }
      }
    } else if (strncmp(s_atoms, line, i_atoms) == 0) {
      for (i = 0; i < n_atoms; i++) {
        read_lines(1, fp);
        if (error_code == 1)
          return FAIL;
        id = atoi(strtok_null_safe(line, delimeter));
        type = atoi(strtok_null_safe(NULL, delimeter));
        if (s_id == id && s_type == type) {
          xu = atof(strtok_null_safe(NULL, delimeter));
          yu = atof(strtok_null_safe(NULL, delimeter));
          zu = atof(strtok_null_safe(NULL, delimeter));
          stream->x[i] = xu;
          stream->y[i] = yu;
          stream->z[i] = zu;
          stream->atomType = 0;
        }
      }
    }
  }
  return 0;
}
int MallocAtomStream(AtomStream *stream, int nTime) {
  stream->x = (real *)malloc(nTime * sizeof(real));
  stream->y = (real *)malloc(nTime * sizeof(real));
  stream->z = (real *)malloc(nTime * sizeof(real));
  stream->mux = (real *)malloc(nTime * sizeof(real));
  stream->muy = (real *)malloc(nTime * sizeof(real));
  stream->muz = (real *)malloc(nTime * sizeof(real));
}
int FreeAtomStream(AtomStream *stream) {
  free(stream->x);
  free(stream->y);
  free(stream->z);
  free(stream->mux);
  free(stream->muy);
  free(stream->muz);
}
#undef FAIL
#undef SUCCESS
