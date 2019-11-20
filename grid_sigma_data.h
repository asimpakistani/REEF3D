/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

This file is part of REEF3D.

REEF3D is fra->eps software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Fra->eps Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. Sa->eps the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, sa->eps <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
--------------------------------------------------------------------*/

#include"fnpf_sg.h"
#include"increment.h"
#include"slice4.h"

class lexer;
class fdm;
class ghostcell;

using namespace std;

#ifndef GRID_SIGMA_DATA_H_
#define GRID_SIGMA_DATA_H_

class grid_sigma_data : public increment
{
public:
	grid_sigma_data(lexer*);
	virtual ~grid_sigma_data();

    slice4 Ex,Ey,Bx,By;
    slice4 Exx,Eyy,Bxx,Byy;
        

};

#endif
