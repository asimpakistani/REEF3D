/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2022 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Arun Kamath
--------------------------------------------------------------------*/

#include"force_fit.h"
#include"gradient.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"ioflow.h"
#include<sys/stat.h>
#include<sys/types.h>

force_fit::force_fit(lexer* p, fdm_fnpf *c, ghostcell *pgc, int qn) : ID(qn){}

force_fit::~force_fit(){}

void force_fit::ini(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    force_fitprintcount=0;

    // Read cylinder force input - xc,yc,zc,rc,lc
    xc = p->P86_x[ID];
	yc = p->P86_y[ID];
	zc = p->P86_z[ID];
	rc = p->P86_r[ID];
	lc = p->P86_l[ID];
	
    // Open files
    print_ini(p,c,pgc);

    // Ini arrays
//	p->Darray(un, p->knoz);
	//p->Darray(u2n, p->knoz);
//	p->Darray(vn, p->knoz);

    // Ini eta
//	etan=p->wd;
	//eta2n=p->wd;

    // Ini time
    //dtn=0;

    // Ini processor boundaries
	xstart = p->originx;
	ystart = p->originy;
	zstart = p->originz;
	xend = p->endx;
	yend = p->endy;
	zend = p->endz;

}

void force_fit::start(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    if (xc >= xstart && xc < xend && yc >= ystart && yc < yend && zc >= zstart && zc < zend) // cylinder in processor
    {
        i = p->posc_i(xc);
        j = p->posc_j(yc);
        k = p->posc_k(zc);
        

        // Calculate force
        force_fit_force(p,c,pgc);
    }
    else
    {
        Fx = Fy = Fz= 0.0;
    }

    // Sum up to distribute forces
    Fx = pgc->globalsum(Fx);
    Fy = pgc->globalsum(Fy);
    Fz = pgc->globalsum(Fz);

    // Print
    if(p->mpirank==0)
    {
        print_force_fit(p,c,pgc);
    }
}
