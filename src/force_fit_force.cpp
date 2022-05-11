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
Author: Arun Kamath
--------------------------------------------------------------------*/

#include"force_fit.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include <math.h>

/*double force_fit::fbouy(lexer *p, fdm_fnpf *c, ghostcell *pgc) // to calculate buoyancy force
{
    double fbouy = p->W1*p->W22*DZN[KP];
    //double dndt= (3*c->eta(i,j) - 4*etan + eta2n)/(p->dt+dtn);
		 
    return fbouy;
}

double force_fit::dudsig(lexer *p, fdm_fnpf *c, ghostcell *pgc) 	// getting dudsig for ax2 and 3
{
    double dudsig_ = 0;
	double dudsig2_ = 0;
    
    if(k<p->knoz)
    {
        dudsig_ =  (c->U[FIJKp1] - c->U[FIJKm1])/(p->DZN[KP1] + p->DZN[KM1]);
	// dudsig2_ =  (c->U[FIp1JK] - c->U[FIm1JK])/(p->DZN[KP] + p->DZN[KM1]);
    }

    if(k==p->knoz)
    { 
        dudsig_ = (c->U[FIJK] - c->U[FIJKm1])/(p->DZN[KM1]);
	//  dudsig2_ = (c->U[FIJK] - c->U[FIm1JK])/(p->DZN[KP]);
    }

    // cout<<"ddsig: "<<dudsig_<<" ddsig2: "<<dudsig2_<<endl;
	return dudsig_;        
}

double force_fit::dvdsig(lexer *p, fdm_fnpf *c, ghostcell *pgc) 	// getting dvdsig for ax2 and 3
{
	double dvdsig_ = 0;

    if(k<p->knoz)
    { 
        dvdsig_ =  (c->V[FIJKp1] - c->V[FIJKm1])/(p->DZN[KP1] + p->DZN[KM1]);
    }

    if(k==p->knoz)
    {
        dvdsig_ = (c->V[FIJK] - c->V[FIJKm1])/(p->DZN[KM1]);
    }

	return dvdsig_;        
}


double force_fit::dudxi(lexer *p, fdm_fnpf *c, ghostcell *pgc) 	// getting dudxi
{
    return (c->U[FIp1JK] - c->U[FIm1JK])/(p->DXN[IP1] + p->DXN[IM1]); 
}

double force_fit::dvdxi(lexer *p, fdm_fnpf *c, ghostcell *pgc) 	// getting dvdxi
{
    return (c->V[FIJp1K] - c->V[FIJm1K])/(p->DYN[JP1] + p->DYN[JM1]); 
}

*/
/*----------------------------------------------------------------------*/

// bouyancy force

/*
for(k= p->posc_k(zc); k<p->knoz; ++k)// loop for body only
 {
  Fbouy= p->W1*p->W22*Fi[FIJK]*DZN[KP]; // rho*g*phi*dsig
 


// Froude-Krylov force

	Ifkry= p->W1*Fi[FIJK]*DZN[KP]; // rho*g*phi*dsig
	Ffkry= (Ifkry-Ifkryn)/p->dt;

// RD fixed body

	Irdfix= p->W1*Fi[FIJK]*DZN[KP];
	Frdfix= (Irdfix-Irdfixn)/p->dt;
	
// RD FS
	
// Total force
	Ffimx= Fbouy + Ffkry + Frdfix;
	
	Ffimx + = Ffimx;
 }
 */
 
 /* ----------------------------------------------------------------------------------- */
 
 
 
 void force_fit::force_fit_force(lexer* p, fdm_fnpf *c, ghostcell *pgc)
{	


	Fx= Fy= Fz= 0;
	
    for(k=3; k<p->knoz; ++k)
	{
      //  dudsig_= dudsig(p, c, pgc); 
     // double dudsig2_= dudsig(p, c, pgc); // cleanup alt ddsig. diff values, no change to force
      //  dvdsig_= dvdsig(p, c, pgc); 
        
     /*   // Term 1 from eqn (9) of Pakozdi et al (2021) MS
        ax1= (c->U[FIJK] - un[k])/(p->dt);  
     // ax1= (3*c->U[FIJK] - 4*un[k] + u2n[k])/(p->dt+dtn); // 2nd order backward diff
        ay1= (c->V[FIJK] - vn[k])/ (p->dt);
        
        // Term 2
        ax2 = c->U[FIJK]*(dudxi(p,c,pgc) + (dudsig_*p->sigx[FIJK]));
        ay2 = c->V[FIJK]*(dvdxi(p,c,pgc) + (dvdsig_*p->sigy[FIJK]));
        
        // Term 3
        ax3 = (c->W[FIJK] - (p->sig[FIJK]*dndt(p, c, pgc)))* dudsig_*p->sigz[IJ];
        ay3 = (c->W[FIJK] - (p->sig[FIJK]*dndt(p, c, pgc)))* dvdsig_*p->sigz[IJ];
      
        // Sum up acceleration
        ax = ax1 + ax2 + ax3;
        ay = ay1 + ay2 + ay3;
    */
        
        // Force on current strip
			// Fx1 = (p->wd + c->eta(i,j))*((cm*ax*p->W1*PI*rc*rc*p->DZN[KP])+ (cd*c->U[FIJK]*fabs(c->U[FIJK])*0.5*p->W1*2*rc*p->DZN[KP]));
			 // Fy1 = (p->wd + c->eta(i,j))*((cm*ay*p->W1*PI*rc*rc*p->DZN[KP])+ (cd*c->V[FIJK]*fabs(c->V[FIJK])*0.5*p->W1*2*rc*p->DZN[KP]));
        	Fx1=0;
        	Fy1=0;
        	Fz1 = p->W1*p->W22*p->wd*p->sigz[IJ]*PI*rc*rc*p->DZN[KP];
        		
        // Sum up forces
        Fx += Fx1;
        Fy += Fy1;
        Fz += Fz1;
        cout<< "Fz1: "<<Fz1<<" Fz: " <<Fz<< " sigz: " <<p->sigz[IJ]<<" dzn: "<<p->DZN[KP] <<endl;
        // Storing current time step information for next time step gradient calculation
        //dtn=p->dt;
        //u2n[k]= un[k];
      //  un[k] = c->U[FIJK]; 
       // vn[k] = c->V[FIJK];
	 
        // cout<< "ax1: "<<ax1<<" ax2: " <<ax2<<" ax3: " <<ax3<< endl;
	    // cout<<"km1: "<<p->ZN[KM1]<<" Dkm1: "<<p->DZN[KM1]<<" kp: "<<p->ZN[KP]<<" sig: "<<p->sig[FIJK]<<"  uvel: "<<c->U[FIJK]<<"  ax: "<<ax<<endl;
	}
	
	//cout<<"ztot: "<<ztot<<endl;
	
	// store current eta value for gradient in next step
	//eta2n=etan;
	//etan=c->eta(i,j);
	
}