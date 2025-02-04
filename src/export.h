/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"increment.h"


class lexer;
class fdm;
class ghostcell;

using namespace std;

#ifndef EXPORTFILE_H_
#define EXPORTFILE_H_

class exportfile : public increment
{

public:
	exportfile(lexer*,fdm*,ghostcell*);
	virtual ~exportfile();
	virtual void start(lexer*,fdm*,ghostcell*);
	
private:
    virtual void filename(lexer*,fdm*,ghostcell*,int);
    virtual void preproc(lexer*,fdm*,ghostcell*);

    char name[200];
    float ffn;
	int iin;
	double ddn;
	int printcount;
    double **eta;
    
    //header
    double xs,xe,ys,ye,zs,ze;
    
    
};

#endif



