/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"sflow_convection.h"
#include"increment.h"

class sflow_flux;
class fdm2D;

using namespace std;

#ifndef SFLOW_CFOU_H_
#define SFLOW_CFOU_H_

class sflow_cfou : public sflow_convection, public increment
{
public:
	sflow_cfou(lexer*,fdm2D*);
	virtual ~sflow_cfou();

	virtual void start(lexer*,fdm2D*,slice&,int,slice&,slice&);

private:
    double aij(lexer*, fdm2D*, slice&, int, slice&, slice&);
	
	double dx,dy,L;
	double ul,ur,vl,vr;
    
    sflow_flux *pflux;
    
    double ivel1,ivel2,jvel1,jvel2;

};

#endif
