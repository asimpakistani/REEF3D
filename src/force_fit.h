/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2022 Hans Bihs

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
Author: Arun Kamath
--------------------------------------------------------------------*/

#include"nodefill.h"
#include"fieldint5.h"
#include"field5.h"
#include<iostream>
#include<fstream>

class lexer;
class fdm_fnpf;
class ghostcell;

#ifndef FORCE_FIT_H_
#define FORCE_FIT_H_

using namespace std;

class force_fit :  public increment
{

public:
	force_fit(lexer*,fdm_fnpf*,ghostcell*,int);
	virtual ~force_fit();
	virtual void start(lexer*,fdm_fnpf*,ghostcell*);
    virtual void ini(lexer*,fdm_fnpf*,ghostcell*);

private:	
	
    void force_fit_force(lexer*,fdm_fnpf*,ghostcell*);
	void print_force_fit(lexer*,fdm_fnpf*,ghostcell*);
    void print_ini(lexer*,fdm_fnpf*,ghostcell*);
	/*double dndt(lexer*, fdm_fnpf*, ghostcell*);
	double dudsig(lexer*, fdm_fnpf*, ghostcell*);
	double dvdsig(lexer*, fdm_fnpf*, ghostcell*);
	double dudxi(lexer*, fdm_fnpf*, ghostcell*);
	double dvdxi(lexer*, fdm_fnpf*, ghostcell*);*/
	
    // force fit variabes
    double Fx,Fy,Fz, xc,yc,zc,rc,lc;
	//double *un, *u2n, *vn;
    const int ID;
    
    // printing
    char name[100],pname[100],epsvar[100];
    int iin,offset[100];
    float ffn;
    int force_fitprintcount;
    ofstream fout;

    // parallelisation
    double xstart,ystart,zstart,xend,yend,zend;
};

#endif
