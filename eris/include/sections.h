//
//  sections.h
//  eris
//
//  Created by Nuno de Sousa on 18/11/14.
//  Copyright (c) 2014 Nuno de Sousa. All rights reserved.
//

#ifndef eris_sections_h
#define eris_sections_h

vec abs_cross_section(cx_mat const*X, Parameters param, Particle *part)
{
    vec abs_sec;
    abs_sec.zeros(3);
    double abs_secx = 0;
    double abs_secy = 0;
    double abs_secz = 0;

    #pragma omp parallel for reduction(+:abs_secx,abs_secy,abs_secz)
    for(int i = 0; i < param.N_part; i = i +1)
    {
        cx_mat solution(3,3);
        solution.zeros(3,3);
        cx_mat p1(3,3);
        mat eta(3,3), phi(3,3);
        cx_mat identity(3,3);
        
        p1.zeros();
        
        for(int l = 0; l < 3; l = l +1)
        {
            for(int m = 0; m < 3; m = m + 1)
            {
                solution(l,m) = (*X)(3*i+l,m);
            }
        }
        
        identity.eye();
        
        p1 = part[i].polarizability_tensor*solution;
        
        eta = imag(solution.t()*p1);
        phi = ((pow(param.wavenumber,3))/(6*M_PI*param.eps_0*param.eps_m))*real(p1.t()*p1);
        
        abs_secx = abs_secx + (eta(0,0) - phi(0,0));
        abs_secy = abs_secy + (eta(1,1) - phi(1,1));
        abs_secz = abs_secz + (eta(2,2) - phi(2,2));
        
    }
    
    abs_sec(0) = abs_secx;
    abs_sec(1) = abs_secy;
    abs_sec(2) = abs_secz;
    
    //É PRECISO CORRIGIR O |E|^2
    
    abs_sec = (param.wavenumber)/(param.eps_0*param.eps_m*pow(param.E0,2))*abs_sec;
    
    return abs_sec;
}


vec sca_cross_section(cx_mat const*X, mat const*im_parts, Parameters param, Particle *part)
{
    
    vec scat_sec(3);
    scat_sec.zeros();

    mat im_parts_temp;
    im_parts_temp = *im_parts;
    
    for(int i = 0; i < param.N_part; i = i + 1)
    {
        double scat_secx = 0;
        double scat_secy = 0;
        double scat_secz = 0;
        
        #pragma omp parallel for reduction(+:scat_secx,scat_secy,scat_secz)
        for(int k = 0; k < param.N_part; k = k + 1)
        {
            cx_mat p1(3,3), p2(3,3), solution(3,3), solution2(3,3);
            cx_mat output_tensor(3,3), temporary(3,3);
            
            if(i != k)
            {
                for(int l = 0; l < 3; l = l +1)
                {
                    for(int m = 0; m < 3; m = m +1)
                    {
                        solution(l,m) = (*X)(3*k+l,m);
                        solution2(l,m) = (*X)(3*i+l,m);
                        
                        output_tensor(l,m) = im_parts_temp(3*i+l,3*k+m);
                        
                    }
                }

                p1 = part[k].polarizability_tensor*solution;
                p2 = part[i].polarizability_tensor*solution2;
                temporary = p2.t()*(output_tensor*p1);
                                
                scat_secx = scat_secx + real(temporary(0,0));
                scat_secy = scat_secy + real(temporary(1,1));
                scat_secz = scat_secz + real(temporary(2,2));
                
            }
            
            
            
            if( i == k)
            {
                for(int l = 0; l < 3; l = l +1)
                {
                    for(int m = 0; m < 3; m = m +1)
                    {
                        solution(l,m) = (*X)(3*i+l,m);
                    }
                }
                
                p1 = part[i].polarizability_tensor*solution;
                temporary = p1.t()*p1;
                
                scat_secx = scat_secx + (param.wavenumber/(6*M_PI))*real(temporary(0,0));
                scat_secy = scat_secy + (param.wavenumber/(6*M_PI))*real(temporary(1,1));
                scat_secz = scat_secz + (param.wavenumber/(6*M_PI))*real(temporary(2,2));
            }
        }
        
        scat_sec(0) = scat_sec(0) + scat_secx;
        scat_sec(1) = scat_sec(1) + scat_secy;
        scat_sec(2) = scat_sec(2) + scat_secz;
        //fecha aqui
        
        
    }
    
    scat_sec = pow(param.wavenumber,3)/(pow(param.eps_0*param.eps_m,2)*pow(param.E0,2))*scat_sec;
    
    return scat_sec;
}


vec ext_cross_section(cx_mat const*X, Parameters param, Particle *part)
{
    vec ext_sec(3);
    ext_sec.zeros(3);
    
    double ext_secx = 0;
    double ext_secy = 0;
    double ext_secz = 0;
    
    #pragma omp parallel for reduction(+:ext_secx,ext_secy,ext_secz)
    for(int i = 0; i < param.N_part; i = i +1)
    {
        cx_mat solution(3,3);
        solution.zeros(3,3);
        cx_mat p1(3,3);

        mat temp(3,3);
        cx_mat ind_tensor;
        ind_tensor.zeros(3,3);
        dcmplx I (0,1);
        
        ind_tensor(0,0) = param.E0*exp(I*param.u(1)*param.wavenumber*part[i].y+I*param.u(2)*param.wavenumber*part[i].z);
        ind_tensor(1,1) = param.E0*exp(I*param.u(0)*param.wavenumber*part[i].x+I*param.u(2)*param.wavenumber*part[i].z);
        ind_tensor(2,2) = param.E0*exp(I*param.u(0)*param.wavenumber*part[i].x+I*param.u(1)*param.wavenumber*part[i].y);
        
        p1.zeros();
        
        for(int l = 0; l < 3; l = l +1)
        {
            for(int m = 0; m < 3; m = m + 1)
            {
                solution(l,m) = (*X)(3*i+l,m);
            }
        }
        
        p1 = part[i].polarizability_tensor* solution;
        
        temp =imag(ind_tensor.t()*p1);
        
        ext_secx = ext_secx + temp(0,0);
        ext_secy = ext_secy + temp(1,1);
        ext_secz = ext_secz + temp(2,2);
    }
    
    ext_sec(0) = ext_secx;
    ext_sec(1) = ext_secy;
    ext_sec(2) = ext_secz;
    
    ext_sec = param.wavenumber/(param.eps_0*param.eps_m*pow(param.E0,2))*ext_sec;
    
    return ext_sec;

    
}

#endif

