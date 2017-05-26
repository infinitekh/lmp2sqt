#ifndef __snapshot_h__ 
#define __snapshot_h__ 

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#define MAXLINE 256
#include "common.h"
typedef enum { SUCCESS, FAIL } checktype;
typedef long int bigint;
enum{ 
	ERR_NONE
};
int error_code;

static char line[MAXLINE];
#define ATOM_VEL 0x1
#define ATOM_DIPOLE 0x2
typedef struct {
	real *x,*y,*z;
	real *mux,*muy,*muz;
	real *vx,*vy,*vz;
	int atomType;
	int nTime;
} atomstream;
typedef struct {
	int id,type;
	int x,y,z;
	int mux,muy,muz;
	int vx,vy,vz;
} atom_column;
typedef struct {
	int id,type;
	real x,y,z;
	real mux,muy,muz;
	real vx,vy,vz;
	real mu1;
	int atomType;
} atom;


#define PBC_X 0x1
#define PBC_Y 0x2
#define PBC_Z 0x4

typedef struct { 
	real xlow,xhigh;
	real ylow,yhigh;
	real zlow,zhigh;
	int pbcTYPE;
} Box3;


typedef struct  {
	bigint timestep;
	int    n_atoms;
	Box3 box;
	atom* atoms;
} Snapshot;


int make_atom(atom* col,int id, int type, 
		real x, real y, real z) ;
int make_atom_vel(atom* col,int id, int type, 
		real x, real y, real z,
		real vx,real vy, real vz) ;
int make_atom_dipole(atom* col,int id, int type, 
		real x, real y, real z,
		real mux,real muy, real muz) ;
int make_atom_all(atom* col,int id, int type, 
		real x, real y, real z,
		real mux,real muy, real muz,
		real vx,real vy, real vz) ;
Snapshot* read_dump(FILE*);
int read_dump_OnlyCheck(FILE*);
void read_lines(int n,FILE*);
void* error(char[]);
void free_Snapshot(Snapshot *snap) ;
Snapshot* new_Snapshot(bigint timestep, int n) ;


/*-----------------------------------------------------------------------------
 *  edit 160811
 *-----------------------------------------------------------------------------*/
int dump_stream(atomstream *,FILE*, int nTime,int nAtom, int id, int type);
int free_stream(atomstream *);
int malloc_stream(atomstream *, int nTime);

#endif
