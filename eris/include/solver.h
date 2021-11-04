//
//  solver.h
//  eris
//
//  Created by Nuno de Sousa on 20/10/15.
//  Copyright © 2015 Nuno de Sousa. All rights reserved.
//
//*****************************************************************
//
// goal: Solve AX = B
//
// Direct Method
// Use the standard method defined in the armadillo library to
// solve the linear system of equations.
// It is a direct allocation of the result
//
// BiCStab - Biconjugate Stablized method dynamic loader
// Solve the system of equations by using the method described in
// Wikipedia. We start with an initial guess (X0) of the solution
// as .zeros. This method doesn't load the A Matrix.
//
//
// BiCStab - Biconjugate Stablized method static loader
// Solve the system of equations by using the method described in
// Wikipedia. We start with an initial guess (X0) of the solution
// as .zeros.
//
//

#ifndef solver_h
#define solver_h

#include <armadillo>
#include <stdlib.h>
#include "particle.h"
//#include "green_tensor.h"

using namespace std;
using namespace arma;


void direct_solver(Equations *eq, cx_mat const*A, cx_mat const*B)
{
    cout << "Direct Solver Initialized" << endl;
    (*eq).set_X(solve((*A),(*B)));
    cout << "Solver terminated" << endl;
}

void solver_BiCStab_dynamic(Equations *eq, cx_mat const*A, cx_mat const*B, Particle *part, Parameters param, int iterations, double error)
{
    
    cout << "Dynamic BiCStab Solver Initialized not working for the moment.\nProgram stoped." << endl;
    exit(1);
    
    
    
    cout << "Dynamic BiCStab Solver Initialized" << endl;
    int Nx = 3*param.N_part;
    int Ny = 3;
    
    (*eq).zeros_X(Nx, Ny);
    
    //step0 - A guess of x0
    cx_mat X0(Nx,Ny);
    cx_mat prod_mat(Ny,Nx);
    X0.eye();
    prod_mat.eye();
    
    cx_mat r0temp(3,Ny);
    r0temp.zeros();
    
    cx_mat Bshit(3,Nx);
    Bshit.zeros();
 
    Bshit = *(B);
    
    
    cx_mat Btemp(3,3);
    Btemp.zeros();
    
    //step1
    cx_mat r0(Nx,Ny);
    
    
    #pragma omp parallel for
    for(int i = 0; i < param.N_part; i = i + 1)
    {
        for(int j = 0; j < param.N_part; j = j + 1)
        {
            
            vec vector_dist (3);
            vec vector_dist_norm = zeros<vec>(3);
            double dist;
            cx_mat output_tensor(3,3);
            
            if(i != j)
            {
                vector_dist(0)=part[i].x - part[j].x;
                vector_dist(1)=part[i].y - part[j].y;
                vector_dist(2)=part[i].z - part[j].z;
                
                dist = norm(vector_dist,2);
                vector_dist_norm = vector_dist/dist;
                
                output_tensor = calc_green_normv3(vector_dist_norm, output_tensor, param, dist);
                
                //I delete the imaginary part
                
                output_tensor = -(pow(param.wavenumber,2)/(param.eps_0*param.eps_m))*output_tensor*part[j].polarizability_tensor;
                
                for(int k = 0; k < 3; k = k + 1)
                {
                    for(int l = 0; l < 3; l = l + 1)
                    {
                        prod_mat(k,3*j+l) = output_tensor(k,l);
                    }
                }
                
            }
            
        }
        
                cout << "chega aqui1" << endl;
        
        for(int k = 0; k < 3; k = k + 1)
        {
            for(int l = 0; l < 3; l = l + 1)
            {
                Btemp(k,l) = Bshit(3*i+k,l);
            }
        }
        
        cout << "chega aqui2" << endl;
        //A falha está aqui.
        
        cout << conj(prod_mat.t())*X0 << endl;
        
        r0temp = Btemp - conj(prod_mat.t())*X0;
        
        for(int k = 0; k < 3; k = k + 1)
        {
            for(int l = 0; l < 3; l = l + 1)
            {
                r0(3*i+k,l) = r0temp(k,l);
            }
        }
        
    }
    
    cout << "chega aqui3" << endl;
    
    
    
    //-----------------------
    
    
    
    //step2
    cx_mat vr0(Nx,Ny);
    vr0 = r0;
    
    //step3
    dcmplx rho0(1,0);
    dcmplx alpha(1,0);
    dcmplx omega0(1,0);
    dcmplx rhoi, rhoim1, omegai, omegaim1;
    dcmplx beta;
    
    rhoim1 = rho0;
    omegaim1 = omega0;
    
    //step4
    cx_mat v0(Nx,Ny);
    cx_mat p0(Nx,Ny);
    cx_mat vi(Nx,Ny);
    cx_mat pi(Nx,Ny);
    cx_mat vim1(Nx,Ny);
    cx_mat pim1(Nx,Ny);
    cx_mat s(Nx,Ny);
    cx_mat t(Nx,Ny);
    cx_mat xi(Nx,Ny);
    cx_mat xim1(Nx,Ny);
    cx_mat ri(Nx,Ny);
    cx_mat rim1(Nx,Ny);
    
    v0.zeros();
    p0.zeros();
    
    rim1 = r0;
    pim1 = p0;
    vim1 = v0;
    
    
    for(int i = 0; i < iterations; i ++)
    {
        //1
        rhoi = dot(vr0,rim1);
        //cout << "rhoi = " << rhoi << endl;
        //2
        beta = (rhoi/rhoim1)*(alpha/omegaim1);
        //3
        pi = rim1 + beta*(pim1-omegaim1*vim1);
        //4
        
#pragma omp parallel for
        for(int i = 0; i < param.N_part; i = i + 1)
        {
            for(int j = 0; j < param.N_part; j = j + 1)
            {
                
                vec vector_dist (3);
                vec vector_dist_norm = zeros<vec>(3);
                double dist;
                cx_mat output_tensor(3,3);
                
                if(i != j)
                {
                    vector_dist(0)=part[i].x - part[j].x;
                    vector_dist(1)=part[i].y - part[j].y;
                    vector_dist(2)=part[i].z - part[j].z;
                    
                    dist = norm(vector_dist,2);
                    vector_dist_norm = vector_dist/dist;
                    
                    output_tensor = calc_green_normv3(vector_dist_norm, output_tensor, param, dist);
                    
                    //I delete the imaginary part
                    
                    output_tensor = -(pow(param.wavenumber,2)/(param.eps_0*param.eps_m))*output_tensor*part[j].polarizability_tensor;
                    
                    for(int k = 0; k < 3; k = k + 1)
                    {
                        for(int l = 0; l < 3; l = l + 1)
                        {
                            prod_mat(k,3*j+l) = output_tensor(k,l);
                        }
                    }
                    
                }
            }
            
            r0temp = prod_mat*pi;
            
            for(int k = 0; k < 3; k = k + 1)
            {
                for(int l = 0; l < 3; l = l + 1)
                {
                    vi(3*i+k,l) = r0temp(k,l);
                }
            }
            
        }
        
        //vi=(*A)*pi;
        
        //5
        alpha = rhoi/dot(vr0,vi);
        //6
        s = rim1 - alpha*vi;
        //7
        
#pragma omp parallel for
        for(int i = 0; i < param.N_part; i = i + 1)
        {
            for(int j = 0; j < param.N_part; j = j + 1)
            {
                
                vec vector_dist (3);
                vec vector_dist_norm = zeros<vec>(3);
                double dist;
                cx_mat output_tensor(3,3);
                
                if(i != j)
                {
                    vector_dist(0)=part[i].x - part[j].x;
                    vector_dist(1)=part[i].y - part[j].y;
                    vector_dist(2)=part[i].z - part[j].z;
                    
                    dist = norm(vector_dist,2);
                    vector_dist_norm = vector_dist/dist;
                    
                    output_tensor = calc_green_normv3(vector_dist_norm, output_tensor, param, dist);
                    
                    //I delete the imaginary part
                    
                    output_tensor = -(pow(param.wavenumber,2)/(param.eps_0*param.eps_m))*output_tensor*part[j].polarizability_tensor;
                    
                    for(int k = 0; k < 3; k = k + 1)
                    {
                        for(int l = 0; l < 3; l = l + 1)
                        {
                            prod_mat(k,3*j+l) = output_tensor(k,l);
                        }
                    }
                    
                }
            }
            
            r0temp = prod_mat*s;
            
            for(int k = 0; k < 3; k = k + 1)
            {
                for(int l = 0; l < 3; l = l + 1)
                {
                    t(3*i+k,l) = r0temp(k,l);
                }
            }
            
        }
        
        //t = (*A)*s;
        
        //8
        omegai = dot(t,s)/dot(t,t);
        //9
        xi = xim1+ alpha*pi + omegai*s;
        //10 if x is accurate quitedirect2_mex
        
        //11
        ri = s - omegai*t;
        
        rim1 = ri;
        pim1 = pi;
        omegaim1 = omegai;
        vim1 = vi;
        xim1 = xi;
        rhoim1 = rhoi;
        
        if(norm(ri)< error)
        {
            (*eq).set_X(xi);
            //(*X)=xi;
            break;
        }
        
    }
    
    (*eq).set_X(xi);
    //(*X)=xi;
    
    cout << "tolerance = " << norm(ri) << endl;
    cout << "Solver terminated" << endl;
    
    
    
    
}

void solver_BiCStab_static(Equations *eq, cx_mat const*A, cx_mat const*B, Parameters param, int iterations, double error)
{
    cout << "Static BiCStab Solver Initialized" << endl;
    int Nx = 3*param.N_part;
    int Ny = 3;
    
    (*eq).zeros_X(Nx, Ny);
    
    //step0 - A guess of x0
    cx_mat X0(Nx,Ny);
    X0.eye();
    
    //step1
    cx_mat r0(Nx,Ny);
    r0 = (*B) - (*A)*X0;
    
    
    //step2
    cx_mat vr0(Nx,Ny);
    vr0 = r0;
    
    //step3
    dcmplx rho0(1,0);
    dcmplx alpha(1,0);
    dcmplx omega0(1,0);
    dcmplx rhoi, rhoim1, omegai, omegaim1;
    dcmplx beta;
    
    rhoim1 = rho0;
    omegaim1 = omega0;
    
    //step4
    cx_mat v0(Nx,Ny);
    cx_mat p0(Nx,Ny);
    cx_mat vi(Nx,Ny);
    cx_mat pi(Nx,Ny);
    cx_mat vim1(Nx,Ny);
    cx_mat pim1(Nx,Ny);
    cx_mat s(Nx,Ny);
    cx_mat t(Nx,Ny);
    cx_mat xi(Nx,Ny);
    cx_mat xim1(Nx,Ny);
    cx_mat ri(Nx,Ny);
    cx_mat rim1(Nx,Ny);
    
    v0.zeros();
    p0.zeros();
    
    rim1 = r0;
    pim1 = p0;
    vim1 = v0;
    
    
    for(int i = 0; i < iterations; i ++)
    {
        //1
        rhoi = dot(vr0,rim1);
        //cout << "rhoi = " << rhoi << endl;
        //2
        beta = (rhoi/rhoim1)*(alpha/omegaim1);
        //3
        pi = rim1 + beta*(pim1-omegaim1*vim1);
        //4
        vi=(*A)*pi;
        //5
        alpha = rhoi/dot(vr0,vi);
        //6
        s = rim1 - alpha*vi;
        //7
        t = (*A)*s;
        //8
        omegai = dot(t,s)/dot(t,t);
        //9
        xi = xim1+ alpha*pi + omegai*s;
        //10 if x is accurate quitedirect2_mex
        
        //11
        ri = s - omegai*t;
        
        rim1 = ri;
        pim1 = pi;
        omegaim1 = omegai;
        vim1 = vi;
        xim1 = xi;
        rhoim1 = rhoi;
        
        if(norm(ri)< error)
        {
            (*eq).set_X(xi);
            //(*X)=xi;
            break;
        }
        
    }
    
    (*eq).set_X(xi);
    //(*X)=xi;
    
    cout << "tolerance = " << norm(ri) << endl;
    cout << "Solver terminated" << endl;
}

#endif /* solver_h */
