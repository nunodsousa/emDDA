//
//  forces.h
//  eris
//
//  Created by Nuno de Sousa on 28/09/15.
//  Copyright © 2015 Nuno de Sousa. All rights reserved.
//

#ifndef forces_h
#define forces_h

cx_mat calc_green_norm(vec vector_dist_norm, cx_mat output_tensor, Parameters param, double dist);

mat forces_eval(cx_mat const*X, Parameters param, Particle *part, double lambda, bool numerical_verification)
{
    
    double delta = 1E-7;
    mat force;
    force.zeros(3*param.N_part,3);
    mat force_N;
    force_N.zeros(3*param.N_part,3);
    
    cx_mat X_temp;
    X_temp = *X;
    cx_mat pi(3,3), pj(3,3);
    cx_mat Ei(3,3), Ej(3,3);
    //solution.zeros(3,3);
    double r, k;
    dcmplx Im(0,1);
    dcmplx a, b ,c;
    cx_mat diffGx, diffGy, diffGz;
    diffGx.zeros(3,3);
    diffGy.zeros(3,3);
    diffGz.zeros(3,3);
    cx_mat E_temp;
    
    cx_mat E_temp_N;
    
    cx_mat temp_matrix(3,3);
    cx_mat temp_matrix_N(3,3);
    
    //Analytical evaluation
    //X component
    //i is the particle where we evaluate
    
    //X component
    
    for(int i = 0; i < param.N_part; i = i +1) //exerted force over the particle i
    {
        
        E_temp.zeros(3,3);
        E_temp_N.zeros(3,3);
        pi.zeros();
        
        for(int l = 0; l < 3; l = l +1)
        {
            for(int m = 0; m < 3; m = m + 1)
            {
                Ei(l,m) = X_temp(3*i+l,m);
            }
        }
        
        pi = part[i].polarizability_tensor* Ei;
        
        
        for(int j = 0; j < param.N_part; j = j + 1)
        {
            
            for(int l = 0; l < 3; l = l +1)
            {
                for(int m = 0; m < 3; m = m + 1)
                {
                    Ej(l,m) = X_temp(3*j+l,m);
                }
            }
            
            
            pj.zeros();
            pj = part[j].polarizability_tensor*Ej;
            
            
            
            //Analytical diff Green Tensor
            if(i != j)
            {
                vec vector_r(3);
                vec vector_r_norm = zeros<vec>(3);
                cx_mat RoutR(3,3);
                cx_mat dG0t;
                dG0t.eye(3,3);
                
                vector_r(0)=part[i].x - part[j].x;
                vector_r(1)=part[i].y - part[j].y;
                vector_r(2)=part[i].z - part[j].z;
                
                r = norm(vector_r,2);
                vector_r_norm = vector_r/r;
                
                k = param.wavenumber;
                
                //Diagonal elements
                dcmplx dG0 =  (Im - 2./(k*r) - Im*3./pow(k*r,2)+3./pow(k*r,3))*((part[i].x-part[j].x)/r)*k*exp(Im*k*r)/(4.*M_PI*r);
                
                dG0t = dG0*dG0t;
                
                //Out of diagonal elements
                
                dcmplx dG1 = -k*exp(Im*k*r)/(4.*M_PI*r)*(part[i].x-part[j].x)/pow(r,3)*(Im-6./(k*r)-Im*15./pow(k*r,2)+15./pow(k*r,3));
                
                for(int m = 0; m < 3; m = m + 1)
                {
                    for(int n = m; n < 3; n = n + 1)
                    {
                        RoutR(m,n)=vector_r(m)*vector_r(n);
                        RoutR(n,m)=RoutR(m,n);
                        
                    }
                }
                
                RoutR = dG1*RoutR;
                
                cx_mat Rx;
                Rx.zeros(3,3);
                Rx(0,0)=2*vector_r(0);
                Rx(0,1)=vector_r(1);
                Rx(0,2)=vector_r(2);
                Rx(1,0)=vector_r(1);
                Rx(2,0)=vector_r(2);
                
                
                dcmplx G1 = -exp(Im*k*r)/(4.*M_PI*r)*(1.+Im*3./(k*r)-3./pow(k*r,2))*(1./pow(r,2));
                
                Rx = G1*Rx;
                
                diffGx = dG0t + RoutR + Rx;
                
                E_temp = E_temp + diffGx*pj;
                
            }
            
        }
        
        
        temp_matrix = (pow(k,2)/(2.*param.eps_0*param.eps_m))*conj(pi.t())*conj(E_temp);
 
        
        force(3*i,0) = real(temp_matrix(0,0));
        force(3*i,1) = real(temp_matrix(1,1));
        force(3*i,2) = real(temp_matrix(2,2));
        
    }
    
    //Y component
    for(int i = 0; i < param.N_part; i = i +1) //exerted force over the particle i
    {
        
        E_temp.zeros(3,3);
        E_temp_N.zeros(3,3);
        pi.zeros();
        
        for(int l = 0; l < 3; l = l +1)
        {
            for(int m = 0; m < 3; m = m + 1)
            {
                Ei(l,m) = X_temp(3*i+l,m);
            }
        }
        
        pi = part[i].polarizability_tensor* Ei;
        
        
        for(int j = 0; j < param.N_part; j = j + 1)
        {
            
            for(int l = 0; l < 3; l = l +1)
            {
                for(int m = 0; m < 3; m = m + 1)
                {
                    Ej(l,m) = X_temp(3*j+l,m);
                }
            }
            
            
            pj.zeros();
            pj = part[j].polarizability_tensor*Ej;
            
            
            
            //Analytical diff Green Tensor
            if(i != j)
            {
                vec vector_r(3);
                vec vector_r_norm = zeros<vec>(3);
                cx_mat RoutR(3,3);
                cx_mat dG0t;
                dG0t.eye(3,3);
                
                vector_r(0)=part[i].x - part[j].x;
                vector_r(1)=part[i].y - part[j].y;
                vector_r(2)=part[i].z - part[j].z;
                
                r = norm(vector_r,2);
                vector_r_norm = vector_r/r;
                
                k = param.wavenumber;
                
                //Diagonal elements
                dcmplx dG0 =  (Im - 2./(k*r) - Im*3./pow(k*r,2)+3./pow(k*r,3))*((part[i].y-part[j].y)/r)*k*exp(Im*k*r)/(4.*M_PI*r);
                
                dG0t = dG0*dG0t;
                
                //Out of diagonal elements
                
                dcmplx dG1 = -k*exp(Im*k*r)/(4.*M_PI*r)*(part[i].y-part[j].y)/pow(r,3)*(Im-6./(k*r)-Im*15./pow(k*r,2)+15./pow(k*r,3));
                
                for(int m = 0; m < 3; m = m + 1)
                {
                    for(int n = m; n < 3; n = n + 1)
                    {
                        RoutR(m,n)=vector_r(m)*vector_r(n);
                        RoutR(n,m)=RoutR(m,n);
                        
                    }
                }
                
                RoutR = dG1*RoutR;
                
                cx_mat Ry;
                Ry.zeros(3,3);
                Ry(1,1)=2*vector_r(1);
                Ry(0,1)=vector_r(0);
                Ry(1,0)=vector_r(0);
                Ry(1,2)=vector_r(2);
                Ry(2,1)=vector_r(2);
                
                
                dcmplx G1 = -exp(Im*k*r)/(4.*M_PI*r)*(1.+Im*3./(k*r)-3./pow(k*r,2))*(1./pow(r,2));
                
                Ry = G1*Ry;
                
                diffGy = dG0t + RoutR + Ry;
                
                E_temp = E_temp + diffGy*pj;
                
            }
            
            
        }
        
        
        temp_matrix = (pow(k,2)/(2.*param.eps_0*param.eps_m))*conj(pi.t())*conj(E_temp);
        
        force(3*i+1,0) = real(temp_matrix(0,0));
        force(3*i+1,1) = real(temp_matrix(1,1));
        force(3*i+1,2) = real(temp_matrix(2,2));
        
    }
    
    
    //Z component
    for(int i = 0; i < param.N_part; i = i +1) //exerted force over the particle i
    {
        
        E_temp.zeros(3,3);
        E_temp_N.zeros(3,3);
        pi.zeros();
        
        for(int l = 0; l < 3; l = l +1)
        {
            for(int m = 0; m < 3; m = m + 1)
            {
                Ei(l,m) = X_temp(3*i+l,m);
            }
        }
        
        pi = part[i].polarizability_tensor* Ei;
        
        
        for(int j = 0; j < param.N_part; j = j + 1)
        {
            
            for(int l = 0; l < 3; l = l +1)
            {
                for(int m = 0; m < 3; m = m + 1)
                {
                    Ej(l,m) = X_temp(3*j+l,m);
                }
            }
            
            
            pj.zeros();
            pj = part[j].polarizability_tensor*Ej;
            
            
            
            //Analytical diff Green Tensor
            if(i != j)
            {
                vec vector_r(3);
                vec vector_r_norm = zeros<vec>(3);
                cx_mat RoutR(3,3);
                cx_mat dG0t;
                dG0t.eye(3,3);
                
                vector_r(0)=part[i].x - part[j].x;
                vector_r(1)=part[i].y - part[j].y;
                vector_r(2)=part[i].z - part[j].z;
                
                r = norm(vector_r,2);
                vector_r_norm = vector_r/r;
                
                k = param.wavenumber;
                
                //Diagonal elements
                dcmplx dG0 =  (Im - 2./(k*r) - Im*3./pow(k*r,2)+3./pow(k*r,3))*((part[i].z-part[j].z)/r)*k*exp(Im*k*r)/(4.*M_PI*r);
                
                dG0t = dG0*dG0t;
                
                //Out of diagonal elements
                
                dcmplx dG1 = -k*exp(Im*k*r)/(4.*M_PI*r)*(part[i].z-part[j].z)/pow(r,3)*(Im-6./(k*r)-Im*15./pow(k*r,2)+15./pow(k*r,3));
                
                for(int m = 0; m < 3; m = m + 1)
                {
                    for(int n = m; n < 3; n = n + 1)
                    {
                        RoutR(m,n)=vector_r(m)*vector_r(n);
                        RoutR(n,m)=RoutR(m,n);
                        
                    }
                }
                
                RoutR = dG1*RoutR;
                
                cx_mat Rz;
                Rz.zeros(3,3);
                Rz(2,2)=2*vector_r(2);
                Rz(0,2)=vector_r(0);
                Rz(2,0)=vector_r(0);
                Rz(1,2)=vector_r(1);
                Rz(2,1)=vector_r(1);
                
                
                dcmplx G1 = -exp(Im*k*r)/(4.*M_PI*r)*(1.+Im*3./(k*r)-3./pow(k*r,2))*(1./pow(r,2));
                
                Rz = G1*Rz;
                
                diffGz = dG0t + RoutR + Rz;
                
                E_temp = E_temp + diffGz*pj;
                
            }
            
        }
        
        
        temp_matrix = (pow(k,2)/(2.*param.eps_0*param.eps_m))*conj(pi.t())*conj(E_temp);
        
        force(3*i+2,0) = real(temp_matrix(0,0));
        force(3*i+2,1) = real(temp_matrix(1,1));
        force(3*i+2,2) = real(temp_matrix(2,2));
        
    }
     
    
    //cout << "force diff = " <<force << endl;
    
    //planar wave component
    //If the incidence if not in the z direction, this must be changed
    //
    
    for(int i = 0; i < param.N_part; i = i + 1)
    {
        
        E_temp.zeros(3,3);
        E_temp_N.zeros(3,3);
        pi.zeros();
        
        for(int l = 0; l < 3; l = l +1)
        {
            for(int m = 0; m < 3; m = m + 1)
            {
                Ei(l,m) = X_temp(3*i+l,m);
            }
        }
        
        pi = part[i].polarizability_tensor* Ei;
        
        //because the incident wave cames from the z-direction
        
        //param.E0*exp(I*param.u(1)*param.wavenumber*part[i].y+I*param.u(2)*param.wavenumber*part[i].z);
        
        force(3*i+2,0) = force(3*i+2,0) - 0.5*param.wavenumber*param.E0*real(Im*pi(0,0)*exp(-Im*param.wavenumber*part[i].z));
        force(3*i+2,1) = force(3*i+2,1) - 0.5*param.wavenumber*param.E0*real(Im*pi(1,1)*exp(-Im*param.wavenumber*part[i].z));;
        //force(3*i,0) = force(3*i,0) + real(pi(0,0) * conj(Ei(0,0)*Im*param.wavenumber*(part[i].x)*exp(Im*param.wavenumber*r)/(4*M_PI*r)));
        //force(3*i+1,1) = force(3*i+1,1) + real(pi(1,1) * conj(Ei(1,1)*Im*param.wavenumber*(part[i].y)*exp(Im*param.wavenumber*r)/(4*M_PI*r)));
        force(3*i+2,2) = 0-force(3*i+2,2);
        
    }
    
    
    
    
    //------------------------------------------------------------------------------------------------------------
    //Numberical evaluation of the Green Tensor
    //Evaluation of the force
    
    //X component
    if(numerical_verification == true)
    {
        for(int i = 0; i < param.N_part; i = i +1) //exerted force over the particle i
        {
            
            E_temp.zeros(3,3);
            E_temp_N.zeros(3,3);
            pi.zeros();
            k = param.wavenumber;
            
            for(int l = 0; l < 3; l = l +1)
            {
                for(int m = 0; m < 3; m = m + 1)
                {
                    Ei(l,m) = X_temp(3*i+l,m);
                }
            }
            
            pi = part[i].polarizability_tensor* Ei;
            
            
            for(int j = 0; j < param.N_part; j = j + 1)
            {
                
                for(int l = 0; l < 3; l = l +1)
                {
                    for(int m = 0; m < 3; m = m + 1)
                    {
                        Ej(l,m) = X_temp(3*j+l,m);
                    }
                }
                
                
                pj.zeros();
                pj = part[j].polarizability_tensor*Ej;
                
                //Numerical diff Green Tensor
                vec vector_dist_dt (3), vector_dist (3);
                vec vector_dist_norm_dt = zeros<vec>(3);
                vec vector_dist_norm = zeros<vec>(3);
                double dist_dt, dist;
                cx_mat output_tensor_dt(3,3);
                cx_mat output_tensor(3,3);
                cx_mat diffGx_N(3,3);
                
                
                if(i != j)
                {
                    vector_dist_dt(0)=part[i].x+ delta - part[j].x;
                    vector_dist_dt(1)=part[i].y - part[j].y;
                    vector_dist_dt(2)=part[i].z - part[j].z;
                    
                    vector_dist(0)=part[i].x - part[j].x;
                    vector_dist(1)=part[i].y - part[j].y;
                    vector_dist(2)=part[i].z - part[j].z;
                    
                    dist_dt = norm(vector_dist_dt,2);
                    vector_dist_norm_dt = vector_dist_dt/dist_dt;
                    
                    dist = norm(vector_dist,2);
                    vector_dist_norm = vector_dist/dist;
                    
                    output_tensor_dt = calc_green_norm(vector_dist_norm_dt, output_tensor_dt, param, dist_dt);
                    output_tensor = calc_green_norm(vector_dist_norm, output_tensor, param, dist);
                    
                    diffGx_N = (output_tensor_dt - output_tensor)/delta;
                    
                    E_temp_N = E_temp_N + diffGx_N*pj;
                    
                }
                
            }
            
            temp_matrix_N = (pow(k,2)/(2.*param.eps_0*param.eps_m))*conj(pi.t())*conj(E_temp_N);
            
            force_N(3*i,0) = real(temp_matrix_N(0,0));
            force_N(3*i,1) = real(temp_matrix_N(1,1));
            force_N(3*i,2) = real(temp_matrix_N(2,2));

            //return forces
        }
        
        //Y component
        for(int i = 0; i < param.N_part; i = i +1) //exerted force over the particle i
        {
            
            E_temp.zeros(3,3);
            E_temp_N.zeros(3,3);
            pi.zeros();
            
            for(int l = 0; l < 3; l = l +1)
            {
                for(int m = 0; m < 3; m = m + 1)
                {
                    Ei(l,m) = X_temp(3*i+l,m);
                }
            }
            
            pi = part[i].polarizability_tensor* Ei;
            
            
            for(int j = 0; j < param.N_part; j = j + 1)
            {
                
                for(int l = 0; l < 3; l = l +1)
                {
                    for(int m = 0; m < 3; m = m + 1)
                    {
                        Ej(l,m) = X_temp(3*j+l,m);
                    }
                }
                
                
                pj.zeros();
                pj = part[j].polarizability_tensor*Ej;
                
                //Numerical diff Green Tensor
                vec vector_dist_dt (3), vector_dist (3);
                vec vector_dist_norm_dt = zeros<vec>(3);
                vec vector_dist_norm = zeros<vec>(3);
                double dist_dt, dist;
                cx_mat output_tensor_dt(3,3);
                cx_mat output_tensor(3,3);
                cx_mat diffGy_N(3,3);
                
                
                if(i != j)
                {
                    vector_dist_dt(0)=part[i].x - part[j].x;
                    vector_dist_dt(1)=part[i].y+ delta - part[j].y;
                    vector_dist_dt(2)=part[i].z - part[j].z;
                    
                    vector_dist(0)=part[i].x - part[j].x;
                    vector_dist(1)=part[i].y - part[j].y;
                    vector_dist(2)=part[i].z - part[j].z;
                    
                    dist_dt = norm(vector_dist_dt,2);
                    vector_dist_norm_dt = vector_dist_dt/dist_dt;
                    
                    dist = norm(vector_dist,2);
                    vector_dist_norm = vector_dist/dist;
                    
                    output_tensor_dt = calc_green_norm(vector_dist_norm_dt, output_tensor_dt, param, dist_dt);
                    output_tensor = calc_green_norm(vector_dist_norm, output_tensor, param, dist);
                    
                    diffGy_N = (output_tensor_dt - output_tensor)/delta;
                    
                    E_temp_N = E_temp_N + diffGy_N*pj;
                    
                }
                
            }
            
            temp_matrix_N = (pow(k,2)/(2.*param.eps_0*param.eps_m))*conj(pi.t())*conj(E_temp_N);
            
            force_N(3*i+1,0) = real(temp_matrix_N(0,0));
            force_N(3*i+1,1) = real(temp_matrix_N(1,1));
            force_N(3*i+1,2) = real(temp_matrix_N(2,2));
            //return forces
        }
        
        
        //Z component
        for(int i = 0; i < param.N_part; i = i +1) //exerted force over the particle i
        {
            
            E_temp.zeros(3,3);
            E_temp_N.zeros(3,3);
            pi.zeros();
            
            for(int l = 0; l < 3; l = l +1)
            {
                for(int m = 0; m < 3; m = m + 1)
                {
                    Ei(l,m) = X_temp(3*i+l,m);
                }
            }
            
            pi = part[i].polarizability_tensor* Ei;
            
            
            for(int j = 0; j < param.N_part; j = j + 1)
            {
                
                for(int l = 0; l < 3; l = l +1)
                {
                    for(int m = 0; m < 3; m = m + 1)
                    {
                        Ej(l,m) = X_temp(3*j+l,m);
                    }
                }
                
                
                pj.zeros();
                pj = part[j].polarizability_tensor*Ej;
                
                //Numerical diff Green Tensor
                vec vector_dist_dt (3), vector_dist (3);
                vec vector_dist_norm_dt = zeros<vec>(3);
                vec vector_dist_norm = zeros<vec>(3);
                double dist_dt, dist;
                cx_mat output_tensor_dt(3,3);
                cx_mat output_tensor(3,3);
                cx_mat diffGz_N(3,3);
                
                
                if(i != j)
                {
                    vector_dist_dt(0)=part[i].x - part[j].x;
                    vector_dist_dt(1)=part[i].y - part[j].y;
                    vector_dist_dt(2)=part[i].z + delta - part[j].z;
                    
                    vector_dist(0)=part[i].x - part[j].x;
                    vector_dist(1)=part[i].y - part[j].y;
                    vector_dist(2)=part[i].z - part[j].z;
                    
                    dist_dt = norm(vector_dist_dt,2);
                    vector_dist_norm_dt = vector_dist_dt/dist_dt;
                    
                    dist = norm(vector_dist,2);
                    vector_dist_norm = vector_dist/dist;
                    
                    output_tensor_dt = calc_green_norm(vector_dist_norm_dt, output_tensor_dt, param, dist_dt);
                    output_tensor = calc_green_norm(vector_dist_norm, output_tensor, param, dist);
                    
                    diffGz_N = (output_tensor_dt - output_tensor)/delta;
                    
                    E_temp_N = E_temp_N + diffGz_N*pj;
                    
                }
                
            }
            
            temp_matrix_N = (pow(k,2)/(2.*param.eps_0*param.eps_m))*conj(pi.t())*conj(E_temp_N);
            
            force_N(3*i+2,0) = real(temp_matrix_N(0,0));
            force_N(3*i+2,1) = real(temp_matrix_N(1,1));
            force_N(3*i+2,2) = real(temp_matrix_N(2,2));
            
            //return forces
            +
        }
    }
    
    if(numerical_verification == true)
    {
        cout << "force diff Num = " <<force_N << endl;
        cout << "force diff = " <<(force_N - force)/force_N << endl;
    }
    
    
    
    return force;

}


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



#endif /* forces_h */
