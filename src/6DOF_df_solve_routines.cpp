/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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
--------------------------------------------------------------------*/

#include"6DOF_df_object.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void sixdof_df_object::rk3(lexer *p, fdm *a, ghostcell *pgc, double alpha)
{
    if (alpha == 1.0)
    {
        pk_ = p_;
        ck_ = c_;
        hk_ = h_;
        ek_ = e_;

        get_trans(p,a,pgc, dp_, dc_, p_, c_);    
        get_rot(dh_, de_, h_, e_);

        p_ = p_ + alpha*p->dt*dp_;
        c_ = c_ + alpha*p->dt*dc_;
        h_ = h_ + alpha*p->dt*dh_;
        e_ = e_ + alpha*p->dt*de_;
    }
    else if (alpha == 0.25)
    {
        get_trans(p,a,pgc, dp_, dc_, p_, c_);    
        get_rot(dh_, de_, h_, e_);        
        
        p_ = 0.75*pk_ + 0.25*p_ + alpha*p->dt*dp_;
        c_ = 0.75*ck_ + 0.25*c_ + alpha*p->dt*dc_;
        h_ = 0.75*hk_ + 0.25*h_ + alpha*p->dt*dh_;
        e_ = 0.75*ek_ + 0.25*e_ + alpha*p->dt*de_;         
    }  
    else
    {
        get_trans(p,a,pgc, dp_, dc_, p_, c_);    
        get_rot(dh_, de_, h_, e_);        
        
        p_ = (1.0/3.0)*pk_ + (2.0/3.0)*p_ + alpha*p->dt*dp_;
        c_ = (1.0/3.0)*ck_ + (2.0/3.0)*c_ + alpha*p->dt*dc_;
        h_ = (1.0/3.0)*hk_ + (2.0/3.0)*h_ + alpha*p->dt*dh_;
        e_ = (1.0/3.0)*ek_ + (2.0/3.0)*e_ + alpha*p->dt*de_;         
    }
}


void sixdof_df_object::rk4(lexer *p, fdm *a, ghostcell *pgc, double alpha)
{
    get_trans(p,a,pgc, dp_, dc_, p_, c_);    
    get_rot(dh_, de_, h_, e_);
        
    Eigen::Vector3d p1 = alpha*p->dt*dp_;
    Eigen::Vector3d c1 = alpha*p->dt*dc_;
    Eigen::Vector3d h1 = alpha*p->dt*dh_;
    Eigen::Vector4d e1 = alpha*p->dt*de_;

    pk_ = p_ + 0.5*p1;
    ck_ = c_ + 0.5*c1;
    hk_ = h_ + 0.5*h1;
    ek_ = e_ + 0.5*e1;
        
    get_trans(p,a,pgc, dp_, dc_, pk_, ck_);    
    get_rot(dh_, de_, hk_, ek_);        
        
    Eigen::Vector3d p2 = alpha*p->dt*dp_;
    Eigen::Vector3d c2 = alpha*p->dt*dc_;
    Eigen::Vector3d h2 = alpha*p->dt*dh_;
    Eigen::Vector4d e2 = alpha*p->dt*de_;
        
    pk_ = p_ + 0.5*p2;
    ck_ = c_ + 0.5*c2;
    hk_ = h_ + 0.5*h2;
    ek_ = e_ + 0.5*e2;        
        
    get_trans(p,a,pgc, dp_, dc_, pk_, ck_);    
    get_rot(dh_, de_, hk_, ek_);        
        
    Eigen::Vector3d p3 = alpha*p->dt*dp_;
    Eigen::Vector3d c3 = alpha*p->dt*dc_;
    Eigen::Vector3d h3 = alpha*p->dt*dh_;
    Eigen::Vector4d e3 = alpha*p->dt*de_;
        
    pk_ = p_ + 0.5*p3;
    ck_ = c_ + 0.5*c3;
    hk_ = h_ + 0.5*h3;
    ek_ = e_ + 0.5*e3;  

    get_trans(p,a,pgc, dp_, dc_, pk_, ck_);    
    get_rot(dh_, de_, hk_, ek_);        
        
    p_ += 1.0/6.0*(p1 + 2.0*p2 + 2.0*p3 + alpha*p->dt*dp_);
    c_ += 1.0/6.0*(c1 + 2.0*c2 + 2.0*c3 + alpha*p->dt*dc_);
    h_ += 1.0/6.0*(h1 + 2.0*h2 + 2.0*h3 + alpha*p->dt*dh_);
    e_ += 1.0/6.0*(e1 + 2.0*e2 + 2.0*e3 + alpha*p->dt*de_);     
}


Eigen::Vector3d sixdof_df_object::ab4_3
(
    lexer *p,
    const Eigen::Vector3d& fn, 
    const Eigen::Vector3d& fn1, 
    const Eigen::Vector3d& fn2, 
    const Eigen::Vector3d& fn3,
    double alpha
)
{
    Eigen::Vector3d ab4_3 = 
        (
            alpha*p->dt * (((((fn2 - fn3) * dtn1 + dtn3 * (fn - fn1)) * pow(dtn2,2) 
            + ((fn2 - fn3) * pow(dtn1,2) - 2 * dtn3 * (fn1 - fn2) * dtn1 + pow(dtn3,2) 
            * (fn - fn1)) * dtn2 - dtn1 * dtn3 * (dtn1 + dtn3) * (fn1 - fn2)) * pow(alpha*p->dt,3)) 
            + (((0.4e1 / 0.3e1 * fn2 - 0.4e1 / 0.3e1 * fn3) * dtn1 + 0.8e1 / 0.3e1 * dtn3 
            * (fn - fn1)) * (pow(dtn2,3)) + (((4 * fn2 - 4 * fn3) * pow(dtn1,2) + 4 * dtn3 
            * (fn - 2 * fn1 + fn2) * dtn1 + 4 * pow(dtn3,2) * (fn - fn1)) * pow(dtn2,2)) 
            + ((-0.8e1 / 0.3e1 * fn3 + 0.8e1 / 0.3e1 * fn2) * (pow(dtn1,3)) - (8 * dtn3 * (fn1 - fn2) * pow(dtn1,2)) 
            + (4 * pow(dtn3,2) * (fn - 2 * fn1 + fn2) * dtn1) + 0.4e1 / 0.3e1 * (pow(dtn3,3)) * (fn - fn1)) * dtn2 
            - 0.8e1 / 0.3e1 * dtn3 * dtn1 * (dtn1 + dtn3 / 0.2e1) * (dtn1 + dtn3) * (fn1 - fn2)) * (pow(alpha*p->dt,2)) 
            + ((2 * dtn3 * (fn - fn1) * pow(dtn2,4)) + (((2 * fn2 - 2 * fn3) * pow(dtn1,2) + 8 * dtn3 * (fn - fn1) 
            * dtn1 + 4 * pow(dtn3,2) * (fn - fn1)) * pow(dtn2,3)) + (((4 * fn2 - 4 * fn3) * pow(dtn1,3) + 6 * dtn3 
            * (fn - 2 * fn1 + fn2) * pow(dtn1,2) + 12 * pow(dtn3,2) * (fn - fn1) * dtn1 + 2 * pow(dtn3,3) * (fn - fn1)) * pow(dtn2,2)) 
            + 0.6e1 * dtn1 * ((fn2 / 0.3e1 - fn3 / 0.3e1) * (pow(dtn1,3)) - 0.4e1 / 0.3e1 * dtn3 * (fn1 - fn2) * (pow(dtn1,2)) 
            + (pow(dtn3,2) * (fn - 2 * fn1 + fn2) * dtn1) + 0.2e1 / 0.3e1 * (pow(dtn3,3)) * (fn - fn1)) * dtn2 - (2 * pow(dtn1,2) * dtn3
            * pow(dtn1 + dtn3,2) * (fn1 - fn2))) * alpha*p->dt + (4 * fn * dtn1 * dtn2 * dtn3 * (dtn2 + dtn3) * (dtn1 + dtn2) 
            * (dtn1 + dtn2 + dtn3))) / dtn1 / (dtn1 + dtn2) / (dtn1 + dtn2 + dtn3) / dtn2 / (dtn2 + dtn3) / dtn3 / 0.4e1
        ); 
  
    return ab4_3;
}


Eigen::MatrixXd sixdof_df_object::ab4_matrix
(
    lexer *p,
    const Eigen::MatrixXd& fn, 
    const Eigen::MatrixXd& fn1, 
    const Eigen::MatrixXd& fn2, 
    const Eigen::MatrixXd& fn3,
    double alpha
)
{
    Eigen::MatrixXd ab4_3 = 
        (
            alpha*p->dt * (((((fn2 - fn3) * dtn1 + dtn3 * (fn - fn1)) * pow(dtn2,2) 
            + ((fn2 - fn3) * pow(dtn1,2) - 2 * dtn3 * (fn1 - fn2) * dtn1 + pow(dtn3,2) 
            * (fn - fn1)) * dtn2 - dtn1 * dtn3 * (dtn1 + dtn3) * (fn1 - fn2)) * pow(alpha*p->dt,3)) 
            + (((0.4e1 / 0.3e1 * fn2 - 0.4e1 / 0.3e1 * fn3) * dtn1 + 0.8e1 / 0.3e1 * dtn3 
            * (fn - fn1)) * (pow(dtn2,3)) + (((4 * fn2 - 4 * fn3) * pow(dtn1,2) + 4 * dtn3 
            * (fn - 2 * fn1 + fn2) * dtn1 + 4 * pow(dtn3,2) * (fn - fn1)) * pow(dtn2,2)) 
            + ((-0.8e1 / 0.3e1 * fn3 + 0.8e1 / 0.3e1 * fn2) * (pow(dtn1,3)) - (8 * dtn3 * (fn1 - fn2) * pow(dtn1,2)) 
            + (4 * pow(dtn3,2) * (fn - 2 * fn1 + fn2) * dtn1) + 0.4e1 / 0.3e1 * (pow(dtn3,3)) * (fn - fn1)) * dtn2 
            - 0.8e1 / 0.3e1 * dtn3 * dtn1 * (dtn1 + dtn3 / 0.2e1) * (dtn1 + dtn3) * (fn1 - fn2)) * (pow(alpha*p->dt,2)) 
            + ((2 * dtn3 * (fn - fn1) * pow(dtn2,4)) + (((2 * fn2 - 2 * fn3) * pow(dtn1,2) + 8 * dtn3 * (fn - fn1) 
            * dtn1 + 4 * pow(dtn3,2) * (fn - fn1)) * pow(dtn2,3)) + (((4 * fn2 - 4 * fn3) * pow(dtn1,3) + 6 * dtn3 
            * (fn - 2 * fn1 + fn2) * pow(dtn1,2) + 12 * pow(dtn3,2) * (fn - fn1) * dtn1 + 2 * pow(dtn3,3) * (fn - fn1)) * pow(dtn2,2)) 
            + 0.6e1 * dtn1 * ((fn2 / 0.3e1 - fn3 / 0.3e1) * (pow(dtn1,3)) - 0.4e1 / 0.3e1 * dtn3 * (fn1 - fn2) * (pow(dtn1,2)) 
            + (pow(dtn3,2) * (fn - 2 * fn1 + fn2) * dtn1) + 0.2e1 / 0.3e1 * (pow(dtn3,3)) * (fn - fn1)) * dtn2 - (2 * pow(dtn1,2) * dtn3
            * pow(dtn1 + dtn3,2) * (fn1 - fn2))) *alpha* p->dt + (4 * fn * dtn1 * dtn2 * dtn3 * (dtn2 + dtn3) * (dtn1 + dtn2) 
            * (dtn1 + dtn2 + dtn3))) / dtn1 / (dtn1 + dtn2) / (dtn1 + dtn2 + dtn3) / dtn2 / (dtn2 + dtn3) / dtn3 / 0.4e1
        ); 
  
    return ab4_3;
}


Eigen::Vector4d sixdof_df_object::ab4_4
(
    lexer *p,
    const Eigen::Vector4d& fn, 
    const Eigen::Vector4d& fn1, 
    const Eigen::Vector4d& fn2, 
    const Eigen::Vector4d& fn3,
    double alpha
)
{
    Eigen::Vector4d ab4_4 =   
    (
        alpha*p->dt * (((((fn2 - fn3) * dtn1 + dtn3 * (fn - fn1)) * pow(dtn2,2) 
        + ((fn2 - fn3) * pow(dtn1,2) - 2 * dtn3 * (fn1 - fn2) * dtn1 + pow(dtn3,2) 
        * (fn - fn1)) * dtn2 - dtn1 * dtn3 * (dtn1 + dtn3) * (fn1 - fn2)) * pow(alpha*p->dt,3)) 
        + (((0.4e1 / 0.3e1 * fn2 - 0.4e1 / 0.3e1 * fn3) * dtn1 + 0.8e1 / 0.3e1 * dtn3 
        * (fn - fn1)) * (pow(dtn2,3)) + (((4 * fn2 - 4 * fn3) * pow(dtn1,2) + 4 * dtn3 
        * (fn - 2 * fn1 + fn2) * dtn1 + 4 * pow(dtn3,2) * (fn - fn1)) * pow(dtn2,2)) 
        + ((-0.8e1 / 0.3e1 * fn3 + 0.8e1 / 0.3e1 * fn2) * (pow(dtn1,3)) - (8 * dtn3 * (fn1 - fn2) * pow(dtn1,2)) 
        + (4 * pow(dtn3,2) * (fn - 2 * fn1 + fn2) * dtn1) + 0.4e1 / 0.3e1 * (pow(dtn3,3)) * (fn - fn1)) * dtn2 
        - 0.8e1 / 0.3e1 * dtn3 * dtn1 * (dtn1 + dtn3 / 0.2e1) * (dtn1 + dtn3) * (fn1 - fn2)) * (pow(alpha*p->dt,2)) 
        + ((2 * dtn3 * (fn - fn1) * pow(dtn2,4)) + (((2 * fn2 - 2 * fn3) * pow(dtn1,2) + 8 * dtn3 * (fn - fn1) 
        * dtn1 + 4 * pow(dtn3,2) * (fn - fn1)) * pow(dtn2,3)) + (((4 * fn2 - 4 * fn3) * pow(dtn1,3) + 6 * dtn3 
        * (fn - 2 * fn1 + fn2) * pow(dtn1,2) + 12 * pow(dtn3,2) * (fn - fn1) * dtn1 + 2 * pow(dtn3,3) * (fn - fn1)) * pow(dtn2,2)) 
        + 0.6e1 * dtn1 * ((fn2 / 0.3e1 - fn3 / 0.3e1) * (pow(dtn1,3)) - 0.4e1 / 0.3e1 * dtn3 * (fn1 - fn2) * (pow(dtn1,2)) 
        + (pow(dtn3,2) * (fn - 2 * fn1 + fn2) * dtn1) + 0.2e1 / 0.3e1 * (pow(dtn3,3)) * (fn - fn1)) * dtn2 - (2 * pow(dtn1,2) * dtn3
        * pow(dtn1 + dtn3,2) * (fn1 - fn2))) * alpha*p->dt + (4 * fn * dtn1 * dtn2 * dtn3 * (dtn2 + dtn3) * (dtn1 + dtn2) 
        * (dtn1 + dtn2 + dtn3))) / dtn1 / (dtn1 + dtn2) / (dtn1 + dtn2 + dtn3) / dtn2 / (dtn2 + dtn3) / dtn3 / 0.4e1
    );  
    
    return ab4_4;
}


Eigen::Vector3d sixdof_df_object::am4_3
(
    lexer *p,
    const Eigen::Vector3d& f,
    const Eigen::Vector3d& fn, 
    const Eigen::Vector3d& fn1, 
    const Eigen::Vector3d& fn2,
    double alpha
)
{
    Eigen::Vector3d am4_3 = 
    (
        ((2 * dtn2 * (f + fn) * pow(dtn1,4)) + (((4 * f + 4 * fn) * pow(dtn2,2)) + 0.8e1 / 0.3e1 *alpha* p->dt * (f + 2 * fn) * 
        dtn2 - 0.2e1 / 0.3e1 * (pow(alpha*p->dt,2)) * (fn1 - fn2)) * (pow(dtn1,3)) + (((2 * f + 2 * fn) * pow(dtn2,3) + 4 * alpha*p->dt * (f + 2 * fn) * 
        pow(dtn2,2) + pow(alpha*p->dt,2) * (f + 5 * fn - 2 * fn1) * dtn2 - pow(alpha*p->dt,3) * (fn1 - fn2)) * pow(dtn1,2)) + alpha*p->dt * ((0.4e1 / 0.3e1 * 
        f + 0.8e1 / 0.3e1 * fn) * (pow(dtn2,3)) + (alpha*p->dt * (f + 5 * fn - 2 * fn1) * pow(dtn2,2)) + (2 * pow(alpha*p->dt,2) * (fn - fn1) * 
        dtn2) - (pow(alpha*p->dt,3) * (fn1 - fn2)) / 0.3e1) * dtn1 + (pow(alpha*p->dt,2) * dtn2 * (alpha*p->dt + 2 * dtn2) * (alpha*p->dt + dtn2) * 
        (fn - fn1)) / 0.3e1) * alpha*p->dt / (alpha*p->dt + dtn1) / (alpha*p->dt + dtn1 + dtn2) / dtn1 / (dtn1 + dtn2) / dtn2 / 0.4e1
    );
    
    return am4_3;
}


Eigen::Vector4d sixdof_df_object::am4_4
(
    lexer *p,
    const Eigen::Vector4d& f,
    const Eigen::Vector4d& fn, 
    const Eigen::Vector4d& fn1, 
    const Eigen::Vector4d& fn2,
    double alpha
)
{
    Eigen::Vector4d am4_4 = 
    (
        ((2 * dtn2 * (f + fn) * pow(dtn1,4)) + (((4 * f + 4 * fn) * pow(dtn2,2)) + 0.8e1 / 0.3e1 * alpha*p->dt * (f + 2 * fn) * 
        dtn2 - 0.2e1 / 0.3e1 * (pow(alpha*p->dt,2)) * (fn1 - fn2)) * (pow(dtn1,3)) + (((2 * f + 2 * fn) * pow(dtn2,3) + 4 * alpha*p->dt * (f + 2 * fn) * 
        pow(dtn2,2) + pow(alpha*p->dt,2) * (f + 5 * fn - 2 * fn1) * dtn2 - pow(alpha*p->dt,3) * (fn1 - fn2)) * pow(dtn1,2)) + alpha*p->dt * ((0.4e1 / 0.3e1 * 
        f + 0.8e1 / 0.3e1 * fn) * (pow(dtn2,3)) + (alpha*p->dt * (f + 5 * fn - 2 * fn1) * pow(dtn2,2)) + (2 * pow(alpha*p->dt,2) * (fn - fn1) * 
        dtn2) - (pow(alpha*p->dt,3) * (fn1 - fn2)) / 0.3e1) * dtn1 + (pow(alpha*p->dt,2) * dtn2 * (alpha*p->dt + 2 * dtn2) * (alpha*p->dt + dtn2) * 
        (fn - fn1)) / 0.3e1) * alpha*p->dt / (alpha*p->dt + dtn1) / (alpha*p->dt + dtn1 + dtn2) / dtn1 / (dtn1 + dtn2) / dtn2 / 0.4e1
    );
    
    return am4_4;
}
