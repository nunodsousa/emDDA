//
//  sections.h
//  eris
//
//  Created by Nuno de Sousa on 18/11/14.
//  Copyright (c) 2014 Nuno de Sousa. All rights reserved.
//

#ifndef eris_sections_h
#define eris_sections_h
/*
 
It Works, but there is a problem with the signs
I think I solve it already
 
vec abs_cross_section_old_version(Equations eq, Parameters param, Particle *part)
{
    cx_mat solution(3,3);
    solution.zeros(3,3);
    cx_mat p1(3,3);
    vec abs_sec(3);
    abs_sec.zeros(3);
    mat eta(3,3), phi(3,3);
    
    
    for(int i = 0; i < param.N_part; i = i +1)
    {
        p1.zeros();
        
        for(int l = 0; l < 3; l = l +1)
        {
            for(int m = 0; m < 3; m = m + 1)
            {
                solution(l,m) = eq.X(3*i+l,m);
            }
        }
        
        cx_mat identity(3,3);
        identity.eye();
        
        p1 = part[i].polarizability_tensor* solution;
        
        eta = -imag(p1.t()*solution);
        
        cout << "eta = " << eta << endl;
        phi = (pow(param.wavenumber,3))/(6*M_PI*param.eps_0*param.eps_m)*real(p1.t()*p1);
        cout << "phi = " << phi << endl;
        
        for(int j = 0; j < 3; j = j + 1)
        {
            abs_sec(j) = abs_sec(j) + (eta(j,j) - phi(j,j));
        }
    }
    
    //É PRECISO CORRIGIR O |E|^2
    
    abs_sec = (param.wavenumber)/(param.eps_0*param.eps_m*pow(param.E0,2))*abs_sec;
    
    return abs_sec;
}

*/



vec abs_cross_section(Equations eq, Parameters param, Particle *part)
{
    cx_mat solution(3,3);
    solution.zeros(3,3);
    cx_mat p1(3,3);
    vec abs_sec(3);
    abs_sec.zeros(3);
    mat eta(3,3), phi(3,3);
    
    
    for(int i = 0; i < param.N_part; i = i +1)
    {
        p1.zeros();
        
        for(int l = 0; l < 3; l = l +1)
        {
            for(int m = 0; m < 3; m = m + 1)
            {
                solution(l,m) = eq.X(3*i+l,m);
            }
        }
        
        cx_mat identity(3,3);
        identity.eye();
        
        p1 = part[i].polarizability_tensor* solution;
        
        
        //As diagonais não estão a ser sumadas
        //eta = imag(p1*solution.t());
        eta = imag(solution.t()*p1);
        
        //cout << "eta = " << eta << endl;
        phi = ((pow(param.wavenumber,3))/(6*M_PI))*real(p1.t()*p1);
        //cout << "phi = " << phi << endl;
        
        for(int j = 0; j < 3; j = j + 1)
        {
            abs_sec(j) = abs_sec(j) + (eta(j,j) - phi(j,j));
        }
    }
    
    //É PRECISO CORRIGIR O |E|^2
    
    abs_sec = (param.wavenumber)/(pow(param.E0,2))*abs_sec;
    
    return abs_sec;
}

vec sca_cross_section(Equations eq, Parameters param, Particle *part)
{
    cx_mat p1(3,3), p2(3,3), solution(3,3), solution2(3,3);
    cx_mat output_tensor(3,3), temporary(3,3);
    vec scat_sec(3);
    
    scat_sec.zeros();
    
    
    
    for(int i = 0; i < param.N_part; i = i + 1)
    {
        for(int k = 0; k < param.N_part; k = k + 1)
        {
            if(i != k) //possibilidade de exitirem erros por aqui
            {
                for(int l = 0; l < 3; l = l +1)
                {
                    for(int m = 0; m < 3; m = m +1)
                    {
                        solution(l,m) = eq.X(3*k+l,m);
                        solution2(l,m) = eq.X(3*i+l,m);
                        
                        output_tensor(l,m) = eq.im_parts(3*i+l,3*k+m);
                        
                    }
                }
                
                //cout << "Imag_Green = " << output_tensor << endl;
                
                p1 = part[k].polarizability_tensor*solution;
                p2 = part[i].polarizability_tensor*solution2;
                temporary = p2.t()*(output_tensor*p1);
                
                //cout << "temporary = " << temporary << endl;
                
                //scat_sec(0) = scat_sec(0) + real(temporary(0,0));
                
                
                for(int j = 0; j < 3; j++)
                {
                    scat_sec(j) = scat_sec(j) + real(temporary(j,j));
                }
                
                //cout << "scat_sec = " << scat_sec << endl;
                
                
            }
            
            //cout << "temp = " << scat_sec << endl;
            
            if( i == k)
            {
                for(int l = 0; l < 3; l = l +1)
                {
                    for(int m = 0; m < 3; m = m +1)
                    {
                        solution(l,m) = eq.X(3*i+l,m);
                        //output_tensor(l,m) = eq.im_parts(3*i+l,3*i+m);
                        
                    }
                }
                
                p1 = part[i].polarizability_tensor*solution;
                temporary = p1.t()*p1;
                
                //cout << "produto p = " << temporary << endl;
                
                for(int j = 0; j < 3; j++)
                {
                    scat_sec(j) = scat_sec(j) + (param.wavenumber/(6*M_PI))*real(temporary(j,j));
                }
                
                //cout << " a scat_sec = " << scat_sec << endl;
            }
        }
    }
    
    //cout << "temporary scat_sec = " << scat_sec << endl;
    
    /*for(int i = 0; i < param.N_part; i = i +1)
     {
     
     p1.zeros();
     
     for(int l = 0; l < 3; l = l +1)
     {
     for(int m = 0; m < 3; m = m + 1)
     {
     solution(l,m) = eq.X(3*i+l,m);
     }
     }
     
     p1 = part[i].polarizability_tensor* solution;
     
     temporary = p1.t()*p1;
     }*/
    
    scat_sec = pow(param.wavenumber,3)/(pow(param.E0,2))*scat_sec;
    
    return scat_sec;
}


vec ext_cross_section(Equations eq, Parameters param, Particle *part)
{
    cx_mat solution(3,3);
    solution.zeros(3,3);
    cx_mat p1(3,3);
    vec ext_sec(3);
    ext_sec.zeros(3);
    mat temp(3,3);
    cx_mat ind_tensor;
    ind_tensor.zeros(3,3);
    dcmplx I (0,1);
    
    for(int i = 0; i < param.N_part; i = i +1)
    {
        ind_tensor(0,0) = param.E0*exp(I*param.u(1)*param.wavenumber*part[i].y+I*param.u(2)*param.wavenumber*part[i].z);
        ind_tensor(1,1) = param.E0*exp(I*param.u(0)*param.wavenumber*part[i].x+I*param.u(2)*param.wavenumber*part[i].z);
        ind_tensor(2,2) = param.E0*exp(I*param.u(0)*param.wavenumber*part[i].x+I*param.u(1)*param.wavenumber*part[i].y);
        
        p1.zeros();
        
        for(int l = 0; l < 3; l = l +1)
        {
            for(int m = 0; m < 3; m = m + 1)
            {
                solution(l,m) = eq.X(3*i+l,m);
            }
        }
        
        p1 = part[i].polarizability_tensor* solution;
        
        temp =imag(ind_tensor.t()*p1);
        
        for(int j = 0; j < 3; j = j + 1)
        {
            ext_sec(j) = ext_sec(j) + temp(j,j);
        }
    }
    
    ext_sec = (pow(param.wavenumber,1)/pow(param.E0,2))*ext_sec;
    
    return ext_sec;
    
    
}

#endif

