//
//  power.h
//  eris
//
//  Created by Nuno de Sousa on 09/02/14.
//  Copyright (c) 2014 Nuno de Sousa. All rights reserved.
//

#include "green_tensor.h"

#ifndef eris_power_h
#define eris_power_h

cx_mat calc_green_normv2(vec vector_dist_norm, Parameters param);

vec power_method_1(cx_mat const*X, Parameters param, Particle *part)
{

    vec LDOS1(3);
    vec vector_dist(3);
    cx_mat p(3,3);
    vec vector_dist_norm = zeros<vec>(3);
    cx_mat green_tensor(3,3);
    double dist;
    cx_mat solution(3,3);
    cx_mat G_p(3,3);
    cx_mat temporary(3,3);
    mat temporaryr(3,3);
    
    cx_mat X_temp;
    X_temp = *X;

    for(int i = 0; i < param.N_part; i = i + 1)
    {
        vector_dist(0)=part[i].x - param.xsource_pos;
        vector_dist(1)=part[i].y - param.ysource_pos;
        vector_dist(2)=part[i].z - param.zsource_pos;
        
        dist = norm(vector_dist,2);
        
        vector_dist_norm = vector_dist/dist;
        
        green_tensor = calc_green_normv2(vector_dist_norm, param, dist);
        
        for(int m = 0; m < 3; m = m + 1)
        {
            for(int n = 0; n < 3; n = n + 1)
            {
                solution(m,n) = X_temp(3*i+m,n);
            }
       }
        
        p = part[i].polarizability_tensor* solution;
        G_p = G_p + green_tensor*p;
    }
    
    for(int m = 0; m < 3; m = m + 1)
    {
        LDOS1(m)=imag(G_p(m,m));
    }
    
    LDOS1= 1+(6*M_PI)/(param.wavenumber)*LDOS1;
    
    return LDOS1;
}

vec power_method_2(cx_mat const*X, mat const*im_parts, Parameters param, Particle *part,vec *absorption)
{
    vec LDOS2(3);
    vec vector_dist(3);
    cx_mat p1(3,3);
    cx_mat p2(3,3);
    vec vector_dist_norm = zeros<vec>(3);
    cx_mat green_tensor(3,3);
    double dist;
    cx_mat solution(3,3);
    cx_mat solution2(3,3);
    cx_mat G_p(3,3);
    cx_mat temporary(3,3);
    mat temporaryr(3,3);
    mat eta(3,3);
    mat phi(3,3);
    vec real_vector(3);
    mat output_tensor(3,3);
    vec absorved_power(3);
    vec abs_power(3);
    cx_mat identity(3,3);
    identity.eye();
    vec test(3);
    
    test.zeros();
    LDOS2.ones();
    real_vector.zeros();

    cx_mat X_temp;
    mat im_parts_temp;
    
    X_temp = *X;
    im_parts_temp = *im_parts;
    
    // primeira parte
    for(int i = 0; i <  param.N_part; i = i + 1)
    {
        vector_dist(0)=part[i].x - param.xsource_pos;
        vector_dist(1)=part[i].y - param.ysource_pos;
        vector_dist(2)=part[i].z - param.zsource_pos;
        
        dist = norm(vector_dist,2);
        
        vector_dist_norm = vector_dist/dist;
        
        green_tensor = calc_green_normv2(vector_dist_norm, param, dist);
        
        for(int k = 0; k < 3; k = k + 1)
        {
            for(int l = 0; l < 3; l = l + 1)
            {
                solution(k,l) = X_temp(3*i+k,l);
                output_tensor(k,l) = imag(green_tensor(k,l));
            }
            
        }
        
        p1 = part[i].polarizability_tensor* solution;
        temporary = p1.t()*output_tensor;
        
        for(int k = 0; k < 3; k ++)
        {
            real_vector(k) = real_vector(k) + real(temporary(k,k));
        }
    }
    
    for(int k = 0; k < 3; k ++)
    {
        real_vector(k) = real_vector(k)*6*M_PI/param.wavenumber;
    }
    
    for(int k = 0; k < 3; k ++)
    {
        LDOS2(k) = LDOS2(k) + real_vector(k);
    }
    
    real_vector.zeros();
    
    // segunda parte
    
    for(int i = 0; i <  param.N_part; i = i + 1)
    {
        vector_dist(0)= param.xsource_pos - part[i].x;
        vector_dist(1)= param.ysource_pos - part[i].y;
        vector_dist(2)= param.zsource_pos - part[i].z;

        dist = norm(vector_dist,2);
        
        vector_dist_norm = vector_dist/dist;
        
        green_tensor = calc_green_normv2(vector_dist_norm, param, dist);
        
        for(int k = 0; k < 3; k = k + 1)
        {
            for(int l = 0; l < 3; l = l + 1)
            {
                solution(k,l) = X_temp(3*i+k,l);
                output_tensor(k,l) = imag(green_tensor(k,l));
            }
            
        }
        
        //Aqui está diferente do código original
        
        p1 = part[i].polarizability_tensor* solution;
        temporary = output_tensor*p1;
        
        
        
        for(int k = 0; k < 3; k ++)
        {
            real_vector(k) = real_vector(k) + real(temporary(k,k));
        }
    }
    
    for(int k = 0; k < 3; k ++)
    {
        real_vector(k) = real_vector(k)*6*M_PI/param.wavenumber;
    }
    
    for(int k = 0; k < 3; k ++)
    {
        LDOS2(k) = LDOS2(k) + real_vector(k);
    }
    
    real_vector.zeros();
    
    //terceira parte
    
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
                        solution(l,m) = X_temp(3*k+l,m);
                        solution2(l,m) = X_temp(3*i+l,m);

                        output_tensor(l,m) = im_parts_temp(3*i+l,3*k+m);

					}
				}
                
                p1 = part[k].polarizability_tensor*solution;
                p2 = part[i].polarizability_tensor*solution2;
                temporary = p2.t()*(output_tensor*p1);
                
				for(int l = 0; l < 3; l ++)
				{
					real_vector(l) = real_vector(l) + real(temporary(l,l));
				}
			}
		}					
	}
    
    for(int l = 0; l < 3; l ++)
	{
		real_vector(l) = real_vector(l)*6*M_PI/param.wavenumber;
	}
    
    for(int l = 0; l < 3; l ++)
	{
		LDOS2(l) = LDOS2(l) + real_vector(l);
	}
    
	for(int k = 0; k < 3; k ++)
	{
		real_vector(k) = 0.;
	}
    
    // Quarta parte
    
    for(int i = 0; i < param.N_part; i = i +1)
    {
        
        p1.zeros();
        
        for(int l = 0; l < 3; l = l +1)
        {
            for(int m = 0; m < 3; m = m + 1)
            {
                solution(l,m) = X_temp(3*i+l,m);
            }
        }
        
        p1 = part[i].polarizability_tensor* solution;
        
        temporary = p1.t()*p1;
        
        
        for(int l = 0; l < 3; l ++)
        {
            real_vector[l] = real_vector[l] + real(temporary(l,l));
        }
    }
    
    for(int i = 0; i < 3; i++)
    {
        LDOS2(i)=LDOS2(i)+ real_vector(i);
    }

    //Absorption
    abs_power.zeros();
    
    for(int i = 0; i < param.N_part; i = i +1)
	{
        
        p1.zeros();
        
		for(int l = 0; l < 3; l = l +1)
		{
			for(int m = 0; m < 3; m = m + 1)
			{
				solution(l,m) = X_temp(3*i+l,m);
			}
		}
        
        cx_mat identity(3,3);
        identity.eye();
        
        // This must be corrected
        // the polarization have a k^2 and it shouldn't have it.
        // For that reason, the temporaryr have a k^7. I think it should be a k^3
        //pow(param.wavenumber,2)/(param.eps_0*param.eps_m)*
        p1 = part[i].polarizability_tensor* solution;
        temporaryr=-((6*M_PI)/(pow(param.wavenumber,3)*(part[i].volume)))*imag(p1.t()*inv(part[i].dielectric_tensor - identity)*p1);
        
        //cout << "volume = " << part[i].volume << endl;
        
		for(int l = 0; l < 3; l ++)
		{
            abs_power(l) = abs_power(l) + temporaryr(l,l);
		}
	}
    
    for(int k = 0; k < 3; k ++)
	{
		real_vector(k) = 0.;
	}
    
	for(int i = 0; i < 3; i++)
	{
		LDOS2(i)=LDOS2(i) + abs_power(i);
        //cout << "Absorption power["<< i << "] = " << abs_power(i) << endl;
	}
    
    *absorption = abs_power;

    return LDOS2;
}

#endif
