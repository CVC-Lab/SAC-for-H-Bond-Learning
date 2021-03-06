/**-----------------------------------------------------------------------------                                                                                                  
**                                                                                                                                                                               
**  Copyright (C) : Structural Bioinformatics Laboratory, Boston University.                                                                                                                        
**                                                                                                                                                                               
**  This software was developed at the Boston University 2006-2011, by                                                                                                      
**  Structural Bioinformatics Laboratory, as part of NIH funded research.                                                                                                                      
**                                                                                                                                                                               
**  Explicit permission is hereby granted to US Universities and US                                                                                                     
**  Government supported Research Institutions to copy and modify this                                                                                                           
**  software for educational and research purposes, provided copies include                                                                                                      
**  this notice. This software (or modified copies thereof) may not be                                                                                                           
**  distributed to any other institution without express permission from the                                                                                                     
**  Structural Bioinformatics Laboratory and  Boston University. Requests to use this software (or                                                                                 **  modified copies therof) in any other way should be sent to Dima Kozakov,                                                                                                     
**  Department of Biomedical Engineering: "midas@bu.edu".                                                                                                                  
**                                                                                                                                                                               
**---------------------------------------------------------------------------*/
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include _MOL_INCLUDE_

struct zmat* cart2zmat (struct protein* protein)
{
	struct zmat* zmat = (struct zmat*) _mol_malloc (sizeof (struct zmat));
	zmat->natoms = p->natoms;
	zmat->zdats = (struct zdat*) _mol_malloc (p->natoms * sizeof (struct zdat));

	// sums of atom positions
	float sumX = 0;
	float sumY = 0;
	float sumZ = 0;

	int i;
	for (i = 0; i < protein->natoms; i++) // sum the atoms
	{
		sumX += protein->atoms[i].X;
		sumY += protein->atoms[i].Y;
		sumZ += protein->atoms[i].Z;
	}

	struct tvector* tvec = (struct tvector*) _mol_malloc (sizeof (struct tvector));

	tvec->X = -(sumX / protein->natoms);
	tvec->Y = -(sumY / protein->natoms);
	tvec->Z = -(sumZ / protein->natoms);

	return zmat;
}
