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

#include"data_void.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

data_void::data_void(lexer* p, fdm *a, ghostcell* pgc)
{
}

data_void::~data_void()
{
}

void data_void::start(lexer* p, fdm* a, ghostcell* pgc)
{
}

void data_void::print_3D(lexer* p, fdm *a, ghostcell *pgc, ofstream &result)
{
}

void data_void::name_pvtu(lexer *p, fdm *a, ghostcell *pgc, ofstream &result)
{
}

void data_void::name_vtu(lexer *p, fdm *a, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
}

void data_void::offset_vtu(lexer *p, fdm *a, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
}


