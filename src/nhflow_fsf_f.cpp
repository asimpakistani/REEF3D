/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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
--------------------------------------------------------------------*/

#include"nhflow_fsf_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"fluid_update_void.h"
#include"heat.h"
#include"concentration.h"
#include"momentum.h"

nhflow_fsf_f::nhflow_fsf_f(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow) : epsi(p->A440*p->DXM), etark1(p),etark2(p)
{
	pupdate = new fluid_update_void();
}

nhflow_fsf_f::~nhflow_fsf_f()
{
}

void nhflow_fsf_f::start(lexer* p, fdm* a, ghostcell* pgc, ioflow* pflow)
{
    /*pgc->start4(p,a->ro,1);
    pgc->start4(p,a->visc,1);
    
    // fill eta_n
    SLICELOOP4
    a->eta_n(i,j) = a->eta(i,j);

    pgc->gcsl_start4(p,a->eta_n,gcval_phi);
    
    
    // Calculate Eta
    SLICELOOP1
    a->P(i,j)=0.0;
    
    SLICELOOP2
    a->Q(i,j)=0.0;

    ULOOP
    a->P(i,j) += a->u(i,j,k)*p->DZN[KP]*0.5*(p->sigz[IJ]+p->sigz[Ip1J]);

    VLOOP
	a->Q(i,j) += a->v(i,j,k)*p->DZN[KP]*0.5*(p->sigz[IJ]+p->sigz[IJp1]);

	pgc->gcsl_start1(p,a->P,10);
    pgc->gcsl_start2(p,a->Q,11);

    // fsf equation
    SLICELOOP4
    a->eta(i,j) -= p->dt*((a->P(i,j)-a->P(i-1,j))/p->DXN[IP] + (a->Q(i,j)-a->Q(i,j-1))/p->DYN[JP]);	  
    
    pflow->eta_relax(p,pgc,a->eta);
    pgc->gcsl_start4(p,a->eta,1);
    p->sigma_update(p,a,pgc,a->eta);
    */
}

void nhflow_fsf_f::ltimesave(lexer* p, fdm *a, slice &ls)
{
}

void nhflow_fsf_f::update(lexer *p, fdm *a, ghostcell *pgc, slice &f)
{
    pupdate->start(p,a,pgc);
}

void nhflow_fsf_f::ini(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow)
{
}

void nhflow_fsf_f::step1(lexer* p, fdm* a, ghostcell* pgc, ioflow* pflow, field &u, field&v, field&w, double alpha)
{
    SLICELOOP1
    a->P(i,j)=0.0;
    
    SLICELOOP2
    a->Q(i,j)=0.0;

    ULOOP
    a->P(i,j) += u(i,j,k)*p->DZN[KP]*0.5*(p->sigz[IJ]+p->sigz[Ip1J]);

    VLOOP
	a->Q(i,j) += v(i,j,k)*p->DZN[KP]*0.5*(p->sigz[IJ]+p->sigz[IJp1]);

	pgc->gcsl_start1(p,a->P,10);
    pgc->gcsl_start2(p,a->Q,11);

    SLICELOOP4
    etark1(i,j) =      a->eta(i,j)

                -      p->dt*((a->P(i,j)-a->P(i-1,j))/p->DXN[IP] + (a->Q(i,j)-a->Q(i,j-1))/p->DYN[JP]);

    pflow->eta_relax(p,pgc,etark1);
    pgc->gcsl_start4(p,a->eta,1);
    
    p->sigma_update(p,a,pgc,etark1);
    p->omega_update(p,a,pgc,u,v,w);
    
    ULOOP
	a->F(i,j,k) = -PORVAL1*fabs(p->W22)*(etark1(i+1,j)-etark1(i,j))/p->DXP[IP];
    
    VLOOP
	a->G(i,j,k) = -PORVAL2*fabs(p->W22)*(etark1(i,j+1)-etark1(i,j))/p->DYP[JP];
}

void nhflow_fsf_f::step2(lexer* p, fdm* a, ghostcell* pgc, ioflow* pflow, field &u, field&v, field&w, double alpha)
{
    SLICELOOP1
    a->P(i,j)=0.0;
    
    SLICELOOP2
    a->Q(i,j)=0.0;

    ULOOP
    a->P(i,j) += u(i,j,k)*p->DZN[KP]*0.5*(p->sigz[IJ]+p->sigz[Ip1J]);

    VLOOP
	a->Q(i,j) += v(i,j,k)*p->DZN[KP]*0.5*(p->sigz[IJ]+p->sigz[IJp1]);

	pgc->gcsl_start1(p,a->P,10);
    pgc->gcsl_start2(p,a->Q,11);

    SLICELOOP4
    etark2(i,j) = 0.75*a->eta(i,j) + 0.25*etark1(i,j)

                - 0.25*p->dt*((a->P(i,j)-a->P(i-1,j))/p->DXN[IP] + (a->Q(i,j)-a->Q(i,j-1))/p->DYN[JP]);

    pflow->eta_relax(p,pgc,etark2);
    pgc->gcsl_start4(p,a->eta,1);
    
    p->sigma_update(p,a,pgc,etark2);
    p->omega_update(p,a,pgc,u,v,w);
    
    ULOOP
	a->F(i,j,k) = -PORVAL1*fabs(p->W22)*(etark2(i+1,j)-etark2(i,j))/p->DXP[IP];
    
    VLOOP
	a->G(i,j,k) = -PORVAL2*fabs(p->W22)*(etark2(i,j+1)-etark2(i,j))/p->DYP[JP];
}

void nhflow_fsf_f::step3(lexer* p, fdm* a, ghostcell* pgc, ioflow* pflow, field &u, field&v, field&w, double alpha)
{
    SLICELOOP1
    a->P(i,j)=0.0;
    
    SLICELOOP2
    a->Q(i,j)=0.0;

    ULOOP
    a->P(i,j) += u(i,j,k)*p->DZN[KP]*0.5*(p->sigz[IJ]+p->sigz[Ip1J]);

    VLOOP
	a->Q(i,j) += v(i,j,k)*p->DZN[KP]*0.5*(p->sigz[IJ]+p->sigz[IJp1]);

	pgc->gcsl_start1(p,a->P,10);
    pgc->gcsl_start2(p,a->Q,11);

    SLICELOOP4
    a->eta(i,j) = (1.0/3.0)*a->eta(i,j) + (2.0/3.0)*etark2(i,j)

                - (2.0/3.0)*p->dt*((a->P(i,j)-a->P(i-1,j))/p->DXN[IP] + (a->Q(i,j)-a->Q(i,j-1))/p->DYN[JP]);

    pflow->eta_relax(p,pgc,a->eta);
    pgc->gcsl_start4(p,a->eta,1);
    
    p->sigma_update(p,a,pgc,a->eta);
    p->omega_update(p,a,pgc,u,v,w);
    
    ULOOP
	a->F(i,j,k) = -PORVAL1*fabs(p->W22)*(a->eta(i+1,j)-a->eta(i,j))/p->DXP[IP];
    
    VLOOP
	a->G(i,j,k) = -PORVAL2*fabs(p->W22)*(a->eta(i,j+1)-a->eta(i,j))/p->DYP[JP];
}


