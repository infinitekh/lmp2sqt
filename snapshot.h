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
/* #####   EXPORTED TYPE DEFINITIONS   ############################################## */
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
} atom;
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
	int    n_atoms;                               /*!< \brief n : the number of atom */
	Box3 box;                                     /*!< \brief Box property */
	atom* atoms;                                  /*!< \brief atom[0..n-1] */
} Snapshot;
/*!
 *  \brief  struct  for system snapshot
 */


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
