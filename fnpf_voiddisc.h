/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"fnpf_convection.h"
#include"increment.h"

#ifndef FNPF_VOIDDISC_H_
#define FNPF_VOIDDISC_H_

using namespace std;

class fnpf_voiddisc : public fnpf_convection, public increment
{
public:
	fnpf_voiddisc(lexer*);
	virtual ~fnpf_voiddisc();

    virtual double fx(lexer*, field&, double, double);
	virtual double fy(lexer*, field&, double, double);
	virtual double fz(lexer*, field&, double, double);
    
    virtual double sx(lexer*, slice&, double);
	virtual double sy(lexer*, slice&, double);
    virtual double sz(lexer*, double*);

private:
   

	double L,grad;
	
	double gradx, grady, gradz;
	double fu1,fv1,fw1,fu2,fv2,fw2;

};

#endif
