//
//  green_tensor.h
//  eris
//
//  Created by Nuno de Sousa on 13/04/14.
//  Copyright (c) 2014 Nuno de Sousa. All rights reserved.
//

typedef std::complex<double> dcmplx;
using namespace arma;
using namespace std;


cx_mat calc_green_normv2(vec vector_dist_norm, Parameters param, double dist)
{
    dcmplx exponential, factor;
    cx_mat output_tensor(3,3);
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
    
    //cout << "vector = " << vector_dist_norm << endl;
    
    
    
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
    
    //cout << "output_tensor = " << output_tensor << endl;
    
    return output_tensor;
    
}

cx_mat calc_green_normv3(vec vector_dist_norm, cx_mat output_tensor, Parameters param, double dist)
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
