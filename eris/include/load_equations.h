//
//  Equations.h
//  eris
//
//  Created by Nuno de Sousa on 10/01/14.
//  Copyright (c) 2014 Nuno de Sousa. All rights reserved.
//
//Problems to solve:
// In the future, call the calc_green_tensor externally

#ifndef eris_load_equations_h
#define eris_load_equations_h

typedef std::complex<double> dcmplx;
using namespace arma;

class Equations{
    
private:
    cx_mat _A, _B, _X;
    mat _im_parts;
    dcmplx prefactor;
    
public:
    
    Equations():A(_A), B(_B), im_parts(_im_parts), X(_X)
    {
    }
    
    void set_equations(int n)
    {
        _A.eye(3*n,3*n);
        
        this-> _A = A;
    }

    void zeros_X(int n, int m)
    {
        _X.eye(3*n,3*m);
        
        this-> _X = X;
    }
    
    void set_X(cx_mat result_X)
    {
        this-> _X = result_X;
    }
  
    void load_equations(Parameters param, Particle *part)
    {

        _A.zeros(3*param.N_part,3*param.N_part);
        _im_parts.zeros(3*param.N_part,3*param.N_part);
        
        _A.eye(3*param.N_part,3*param.N_part);
        
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
                    
                    output_tensor = calc_green_norm(vector_dist_norm, output_tensor, param, dist);
                    
                    
                    for(int k = 0; k < 3; k = k + 1)
                    {
                        for(int l = 0; l < 3; l = l + 1)
                        {
                            {
                            _im_parts(3*i+k,3*j+l) = imag(output_tensor(k,l));
                            _im_parts(3*j+k,3*i+l) = imag(output_tensor(k,l));
                            }
                            
                            }
                    }
                    output_tensor = -(pow(param.wavenumber,2)/(param.eps_0*param.eps_m))*output_tensor*part[j].polarizability_tensor;
                    
                    for(int k = 0; k < 3; k = k + 1)
                    {
                        for(int l = 0; l < 3; l = l + 1)
                        {
                            _A(3*i+k,3*j+l) = output_tensor(k,l);
                        }
                    }
                    
                }
            }
        }
        
        this-> _A = A;
        this-> _im_parts = im_parts;
    }
    
    
    //Load independent terms.
    //Planar wave
    void load_ind_term_plw(Parameters param, Particle *part)
    {
        
        cx_mat ind_tensor;
        _B.zeros(3*param.N_part,3);
        ind_tensor.zeros(3,3);
        dcmplx I (0,1);
        
        //This must be modified because there is no physical meaning for ind_tensor(2,2) != 0 when the k is in the z direction
        //#pragma omp parallel for
        for(int i = 0; i < param.N_part; i = i + 1)
        {
            ind_tensor(0,0) = param.E0*exp(I*param.u(1)*param.wavenumber*part[i].y+I*param.u(2)*param.wavenumber*part[i].z);
            ind_tensor(1,1) = param.E0*exp(I*param.u(0)*param.wavenumber*part[i].x+I*param.u(2)*param.wavenumber*part[i].z);
            ind_tensor(2,2) = param.E0*exp(I*param.u(0)*param.wavenumber*part[i].x+I*param.u(1)*param.wavenumber*part[i].y);
            
            for(int k = 0; k < 3; k = k + 1)
            {
                for(int l = 0; l < 3; l = l + 1)
                {
                    _B(3*i+k,l)=ind_tensor(k,l);
                }
            }
        }
        
        this-> _B = B;
    }
    
    //load the independent terms.
    //Point source
    void load_ind_term_source(Parameters param, Particle *part)
    {
        
        //This is a example
        //Single source placed at (0,0,0)
        
        vec vector_dist(3);
        vec vector_dist_norm = zeros<vec>(3);
        cx_mat output_tensor(3,3);
        _B.zeros(3*param.N_part,3);
        double dist;
        
        for(int i = 0; i < param.N_part; i = i + 1)
        {
            vector_dist(0)=part[i].x - param.xsource_pos;
            vector_dist(1)=part[i].y - param.ysource_pos;
            vector_dist(2)=part[i].z - param.zsource_pos;
            
            dist = norm(vector_dist,2);
            vector_dist_norm = vector_dist/dist;
            
            output_tensor = calc_green_norm(vector_dist_norm, output_tensor, param, dist);
            
            for(int k = 0; k < 3; k = k + 1)
            {
                for(int l = 0; l < 3; l = l + 1)
                {
                    _B(3*i+k,l)=output_tensor(k,l);
                }
            }
        }
        
        _B=_B*(pow(param.wavenumber,2)/(param.eps_0*param.eps_m));
        
        
        //Copy of the independent term
        _X = _B;
        
        this-> _X = X;
        this-> _B = B;
        
        //Free memory of output tensor
        
    }
    
/*
    void solver(void)
    {
        _X = solve(_A, _B);
        this-> _X = X;
    }
    
    
    void solver_BiCStab(Parameters param, int iterations, double error)
    {
        _X.zeros();
        int Nx = 3*param.N_part;
        int Ny = 3;
        
        //step0 - A guess of x0
        cx_mat X0(Nx,Ny);
        X0.eye();
        
        //step1
        cx_mat r0(Nx,Ny);
        r0 = _B - _A*X0;
        
        
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
            vi=_A*pi;
            //5
            alpha = rhoi/dot(vr0,vi);
            //6
            s = rim1 - alpha*vi;
            //7
            t = _A*s;
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
                
                _X=xi;
                break;
            }
            
        }
        
        _X=xi;
        
        cout << "tolerance = " << norm(ri) << endl;
        this-> _X = X;
    }
*/    
    
    cx_mat calc_green_norm(vec vector_dist_norm, cx_mat output_tensor, Parameters param, double dist)
    {
        dcmplx exponential, factor;
        double kr;
        dcmplx Im(0,1);
        
        for(int m = 0; m < 3; m = m + 1)
        {
            for(int n = m; n < 3; n = n + 1)
            {
                output_tensor(m,n)=vector_dist_norm(m)*vector_dist_norm(n);
                output_tensor(n,m)=output_tensor(m,n);
                
            }
        }
        
        kr = param.wavenumber*dist;
        exponential = exp(Im*kr)/(4*M_PI*dist);
        factor = (1. - Im*kr)/(pow(kr,2));
        
        for(int m = 0; m < 3; m = m + 1)
        {
            for(int n = 0; n < 3; n = n + 1)
            {
                output_tensor(m,n) = (3.*factor-1.)*output_tensor(m,n);
            }
        }
        
        for(int m = 0; m < 3; m = m +1)
        {
            output_tensor(m,m) = output_tensor(m,m) - (factor - 1.);
        }
        
        output_tensor = output_tensor*exponential;
        
        return output_tensor;
        
    }
    
    const cx_mat &A;
    const cx_mat &B;
    const cx_mat &X;
    const mat &im_parts;
};
#endif
