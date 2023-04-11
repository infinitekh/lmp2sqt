/*!
 *    \file  snapshot.h
 *   \brief  
 *
 *  Snaptshot struct definition
 *
 *  \author  KIM Hyeok (kh), ekh0324@gmail.com
 *
 *  \internal
 *       Created:  2017년 05월 29일
 *      Revision:  none
 *      Compiler:  gcc
 *  Organization:  Konkuk University
 *     Copyright:  Copyright (c) 2017, KIM Hyeok
 *
 *  This source code is released for free distribution under the terms of the
 *  GNU General Public License as published by the Free Software Foundation.
 */

#ifndef __snapshot_h__ 
#define __snapshot_h__ 
/* #####   HEADER FILE INCLUDES   ################################################### */
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#define MAXLINE 1024
#include "common.h"
typedef enum { SUCCESS, FAIL } checktype;
typedef long int bigint;
enum{ 
	ERR_NONE
};
static int error_code;

static char line[MAXLINE];
#define ATOM_VEL 0x1
#define ATOM_DIPOLE 0x2
#define ATOM_ALL (ATOM_VEL|ATOM_DIPOLE)
/* #####   EXPORTED TYPE DEFINITIONS   ############################################## */
typedef struct {
	real *x,*y,*z;   
	real *mux,*muy,*muz; 
	real *vx,*vy,*vz;    
	int atomType;
	int nTime;
} AtomStream;
typedef struct {
	int id,type;
	int x,y,z;
	int mux,muy,muz;
	int vx,vy,vz;
} atom_column;
typedef struct {
	int id,                                       ///< unique atom id
			type;                                     ///< atom type
	real x,                                       ///< atom position
			 y,                                       ///< atom position z
			 z;                                       /*!< \brief atom position z */
	real vx,                                      /*!< \brief velocity of x-axis */
			 vy,                                      /*!< \brief velocity of x-axis */
			 vz;                                      /*!< \brief velocity of x-axis */
	int atomType;
} AtomVel;
typedef struct {
	int id,                                       ///< unique atom id
			type;                                     ///< atom type
	real x,                                       ///< atom position
			 y,                                       ///< atom position z
			 z;                                       /*!< \brief atom position z */
	real mux,                                     /*!< \brief magnetic moment of x-axis */
			 muy,                                     /*!< \brief magnetic moment of x-axis */
			 muz;                                     /*!< \brief magnetic moment of x-axis */
	real mu1;                                     /*!< \brief \f$ |\mu| \f$ */
	int atomType;
} AtomMu;
typedef struct {
	int id,                                       ///< unique atom id
			type;                                     ///< atom type
	real x,                                       ///< atom position
			 y,                                       ///< atom position z
			 z;                                       /*!< \brief atom position z */
	real mux,                                     /*!< \brief magnetic moment of x-axis */
			 muy,                                     /*!< \brief magnetic moment of x-axis */
			 muz;                                     /*!< \brief magnetic moment of x-axis */
	real vx,                                      /*!< \brief velocity of x-axis */
			 vy,                                      /*!< \brief velocity of x-axis */
			 vz;                                      /*!< \brief velocity of x-axis */
	real mu1;                                     /*!< \brief \f$ |\mu| \f$ */
	int atomType;
} AtomAll;
typedef struct {
	int id,                                       ///< unique atom id
			type;                                     ///< atom type
	real x,                                       ///< atom position
			 y,                                       ///< atom position z
			 z;                                       /*!< \brief atom position z */
	real mux,                                     /*!< \brief magnetic moment of x-axis */
			 muy,                                     /*!< \brief magnetic moment of x-axis */
			 muz;                                     /*!< \brief magnetic moment of x-axis */
	real vx,                                      /*!< \brief velocity of x-axis */
			 vy,                                      /*!< \brief velocity of x-axis */
			 vz;                                      /*!< \brief velocity of x-axis */
	real mu1;                                     /*!< \brief \f$ |\mu| \f$ */
	int atomType;
} Atom;
/*!
 * \struct atom
 *  \brief atom degree of freedom 3(dimension) x 3 (pos,velo,magn)
 *
 *  \f$ \vec{r}=\{x,y,z\}, \; \vec{\mu} =\{ mux,muy,muz\}, \f$...
 */


#define PBC_X 0x1                               /*!< \brief periodic boundary condition flag of x axis */
#define PBC_Y 0x2                               /*!< \brief periodic boundary condition flag of y axis */
#define PBC_Z 0x4                               /*!< \brief periodic boundary condition flag of z axis */

typedef struct { 
	real xlow                                     /*!< \brief Lower limit of x dimension */
		,xhigh;                                     /*!< \brief Upper limit of x dimension */
	real ylow,                                    /*!< \brief Lower limit of x dimension */
			 yhigh                                    /*!< \brief Upper limit of x dimension */;
	real zlow,                                    /*!< \brief Lower limit of x dimension */
			 zhigh                                    /*!< \brief Upper limit of x dimension */;
	int pbcTYPE;                                  /*!< \brief pbc Type flag. Usable by  bit And operation */
} Box3;
/*!
 *  \brief  system box properties
 */


typedef struct  {
	bigint timestep;                              /*!< \brief t : time  */
        int NumAtoms; /*!< \brief n : the number of atom */
        Box3 Box;     /*!< \brief Box property */
        Atom *atoms;  /*!< \brief atom[0..n-1] */
        int snapFlag;
} Snapshot;
/*!
 *  \brief  struct  for system snapshot
 */

int make_atom(Atom *col, int id, int type, real x, real y, real z);
int make_atom_vel(Atom *col, int id, int type, real x, real y, real z, real vx,
                  real vy, real vz);
int make_atom_dipole(Atom *col, int id, int type, real x, real y, real z,
                     real mux, real muy, real muz);
int make_atom_all(Atom *col, int id, int type, real x, real y, real z, real mux,
                  real muy, real muz, real vx, real vy, real vz);
Snapshot *ReadDump(FILE *);
int ReadDumpForOnlyCheck(FILE *);
void read_lines(int n,FILE*);
void *error(char[MAXLINE]);
void FreeSnapshot(Snapshot *snap);
Snapshot *NewSnapshot(bigint timestep, int n);

/*-----------------------------------------------------------------------------
 *  edit 160811
 *-----------------------------------------------------------------------------*/
int DumpAtomStream(AtomStream *, FILE *, int nTime, int nAtom, int id,
                   int type);
int FreeAtomStream(AtomStream *);
int MallocAtomStream(AtomStream *, int nTime);

#endif
