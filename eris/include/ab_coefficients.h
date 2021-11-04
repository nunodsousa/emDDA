//
//  Header.h
//  eris
//
//  Created by Nuno de Sousa on 03/07/15.
//  Copyright (c) 2015 Nuno de Sousa. All rights reserved.
//

#ifndef eris_Header_h
#define eris_Header_h

double besselj_0(double rho);
double besselj_1(double rho);
double besselj_2(double rho);

dcmplx print_ae01(cx_mat const*X, Parameters param, Particle *part, double lambda, int pol)
{
    cx_mat p1(3,3);
    cx_mat solution(3,3);
    solution.zeros(3,3);
    double dist;
    double kr;
    dcmplx Im(0,1);
    dcmplx b(0,0);
    
    cx_mat X_temp;
    X_temp = *X;
    
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
        
        
        dist = sqrt(pow(part[i].x,2) + pow(part[i].y,2) + pow(part[i].z,2));
        kr = param.wavenumber*dist;
        dcmplx alpha;
        
        if(kr == 0.)
        {
            dcmplx temp(-Im*pow(param.wavenumber,3.)/(6.*M_PI*param.eps_0)*p1(0,pol));
            alpha = temp;
        }
        
        else
        {
            dcmplx prefactor =-(Im*pow(param.wavenumber,3.)/(4.*M_PI*param.eps_0));
            dcmplx z_contrib = ((besselj_1(kr)/kr)*(3.*((pow(part[i].z,2))/pow(dist,2.))-1.)+besselj_0(kr)*(1.-((pow(part[i].z,2.))/pow(dist,2.))))*p1(2,pol);
            dcmplx x_contrib = ((part[i].x*part[i].z)/pow(dist,2.))*besselj_2(kr)*p1(0,pol);
            dcmplx y_contrib = ((part[i].y*part[i].z)/pow(dist,2.))*besselj_2(kr)*p1(1,pol);
            
            dcmplx temp(prefactor*(z_contrib + x_contrib + y_contrib));
            alpha = temp;
        }
        
        
        b = b + alpha;
        
    }
    
    //cout << "ae01 = " << b << endl;
    
    return b;}

dcmplx print_be01(cx_mat const*X, Parameters param, Particle *part, double lambda, int pol)
{
    cx_mat p1(3,3);
    cx_mat solution(3,3);
    solution.zeros(3,3);
    double dist;
    double kr;
    dcmplx b(0,0);
    dcmplx Im(0,1);
    
    cx_mat X_temp;
    X_temp = *X;
    
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
        
        
        dist = sqrt(pow(part[i].x,2) + pow(part[i].y,2) + pow(part[i].z,2));
        kr = param.wavenumber*dist;
        dcmplx alpha;
        
        if(kr == 0.)
        {
            dcmplx temp(0,0);
            alpha = temp;
        }
        
        else
        {
            dcmplx prefactor (-pow(param.wavenumber,3.)/(4.*M_PI*param.eps_0));
            dcmplx x_contrib ((part[i].y)/dist*p1(0,pol)*besselj_1(kr));
            dcmplx y_contrib (-(part[i].x)/dist*p1(1,pol)*besselj_1(kr));
            dcmplx temp(prefactor*(x_contrib + y_contrib));
            alpha = temp;
        }
        
        b = b + alpha;
        
    }
    
    //cout << "be01 = " << b << endl;
    
    return b;
}


dcmplx print_ae11(cx_mat const*X, Parameters param, Particle *part, double lambda, int pol)
{
    cx_mat p1(3,3);
    cx_mat solution(3,3);
    solution.zeros(3,3);
    double dist;
    double kr;
    dcmplx Im(0,1);
    dcmplx b(0,0);
    
    cx_mat X_temp;
    X_temp = *X;
    
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
        
        
        dist = sqrt(pow(part[i].x,2) + pow(part[i].y,2) + pow(part[i].z,2));
        kr = param.wavenumber*dist;
        dcmplx alpha;
        
        if(kr == 0.)
        {
            dcmplx temp(-Im*pow(param.wavenumber,3.)/(6.*M_PI*param.eps_0)*p1(0,pol));
            alpha = temp;
        }
        
        else
        {
            dcmplx prefactor =-(Im*pow(param.wavenumber,3.)/(4.*M_PI*param.eps_0));
            dcmplx x_contrib = ((besselj_1(kr)/kr)*(3.*((pow(part[i].x,2))/pow(dist,2.))-1.)+besselj_0(kr)*(1.-((pow(part[i].x,2.))/pow(dist,2.))))*p1(0,pol);
            dcmplx y_contrib = ((part[i].x*part[i].y)/pow(dist,2.))*besselj_2(kr)*p1(1,pol);
            dcmplx z_contrib = ((part[i].x*part[i].z)/pow(dist,2.))*besselj_2(kr)*p1(2,pol);
            
            dcmplx temp(prefactor*(x_contrib + y_contrib + z_contrib));
            alpha = temp;
        }
        
        
        b = b + alpha;
        
    }
    
    //cout << "ae11 = " << b << endl;
    
    return b;
}


dcmplx print_be11(cx_mat const*X, Parameters param, Particle *part, double lambda, int pol)
{
    cx_mat p1(3,3);
    cx_mat solution(3,3);
    solution.zeros(3,3);
    double dist;
    double kr;
    dcmplx b(0,0);
    dcmplx Im(0,1);
    
    cx_mat X_temp;
    X_temp = *X;
    
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
        
        
        dist = sqrt(pow(part[i].x,2) + pow(part[i].y,2) + pow(part[i].z,2));
        kr = param.wavenumber*dist;
        dcmplx alpha;
        
        if(kr == 0.)
        {
            dcmplx temp(0,0);
            alpha = temp;
        }
        
        else
        {
            dcmplx prefactor (-pow(param.wavenumber,3.)/(4.*M_PI*param.eps_0));
            dcmplx x_contrib ((part[i].z)/dist*p1(1,pol)*besselj_1(kr));
            dcmplx y_contrib (-(part[i].y)/dist*p1(2,pol)*besselj_1(kr));
            dcmplx temp(prefactor*(x_contrib + y_contrib));
            alpha = temp;
        }
        
        b = b + alpha;
        
    }
    
    //cout << "be11 = " << b << endl;
    
    return b;
}


dcmplx print_ao11(cx_mat const*X, Parameters param, Particle *part, double lambda, int pol)
{
    cx_mat p1(3,3);
    cx_mat solution(3,3);
    solution.zeros(3,3);
    double dist;
    double kr;
    dcmplx Im(0,1);
    dcmplx b(0,0);
    
    cx_mat X_temp;
    X_temp = *X;
    
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
        
        
        dist = sqrt(pow(part[i].x,2) + pow(part[i].y,2) + pow(part[i].z,2));
        kr = param.wavenumber*dist;
        dcmplx alpha;
        
        if(kr == 0.)
        {
            dcmplx temp(-Im*pow(param.wavenumber,3.)/(6.*M_PI*param.eps_0)*p1(0,pol));
            alpha = temp;
        }
        
        else
        {
            dcmplx prefactor =-(Im*pow(param.wavenumber,3.)/(4.*M_PI*param.eps_0));
            dcmplx y_contrib = ((besselj_1(kr)/kr)*(3.*((pow(part[i].y,2))/pow(dist,2.))-1.)+besselj_0(kr)*(1.-((pow(part[i].y,2.))/pow(dist,2.))))*p1(1,pol);
            dcmplx x_contrib = ((part[i].x*part[i].y)/pow(dist,2.))*besselj_2(kr)*p1(0,pol);
            dcmplx z_contrib = ((part[i].y*part[i].z)/pow(dist,2.))*besselj_2(kr)*p1(2,pol);
            
            dcmplx temp(prefactor*(x_contrib + y_contrib + z_contrib));
            alpha = temp;
        }
        
        
        b = b + alpha;
        
    }
    
    //cout << "ao11 = " << b << endl;
    
    return b;
}


dcmplx print_bo11(cx_mat const*X, Parameters param, Particle *part, double lambda, int pol)
{
    cx_mat p1(3,3);
    cx_mat solution(3,3);
    solution.zeros(3,3);
    double dist;
    double kr;
    dcmplx b(0,0);
    dcmplx Im(0,1);
    
    cx_mat X_temp;
    X_temp = *X;
    
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
        
        
        dist = sqrt(pow(part[i].x,2) + pow(part[i].y,2) + pow(part[i].z,2));
        kr = param.wavenumber*dist;
        dcmplx alpha;
        
        if(kr == 0.)
        {
            dcmplx temp(0,0);
            alpha = temp;
        }
        
        else
        {
            dcmplx prefactor (-pow(param.wavenumber,3.)/(4.*M_PI*param.eps_0));
            dcmplx x_contrib ((part[i].z)/dist*p1(0,pol)*besselj_1(kr));
            dcmplx z_contrib (-(part[i].x)/dist*p1(2,pol)*besselj_1(kr));
            dcmplx temp(prefactor*(x_contrib + z_contrib));
            alpha = temp;
        }
        
        b = b + alpha;
        
    }
    
    //cout << "bo11 = " << b << endl;
    
    return b;
}





void print_aeo(double lambda, dcmplx ae01, dcmplx ae11, dcmplx ao11, string name)
{
    FILE *file;
    
    //cout << "name = " << name <<endl;
    
    char name1[50];
    for(int i=0;i<50;i++)
        name1[i]=name[i];
    
    
    if(file_exists(name1) == 0) //if doesn't exit
    {
        file = fopen(name1, "w");
        fprintf(file,"#lambda\tre(ae01)\tim(ae01)\tre(ae11)\tim(ae11)\tre(ao11)\tim(ao11)\n");
        fclose(file);
    }
    
    if(file_exists(name1) != 0)
    {
        //cout << "entra em ficheiro existente" << endl;
        file = fopen(name1, "a");
        fprintf(file, "%lf\t%6.6e\t%6.6e\t%6.6e\t%6.6e\t%6.6e\t%6.6e\n", lambda, real(ae01),imag(ae01),real(ae11),imag(ae11),real(ao11),imag(ao11));
        fclose(file);
    }
}

void print_beo(double lambda, dcmplx be01, dcmplx be11, dcmplx bo11, string name)
{
    FILE *file;
    
    //cout << "name = " << name <<endl;
    
    char name1[50];
    for(int i=0;i<50;i++)
        name1[i]=name[i];
    
    if(file_exists(name1) == 0) //if doesn't exit
    {
        file = fopen(name1, "w");
        fprintf(file,"#lambda\tre(be01)\tim(be01)\tre(be11)\tim(be11)\tre(bo11)\tim(bo11)\n");
        fclose(file);
    }
    
    if(file_exists(name1) != 0)
    {
        //cout << "entra em ficheiro existente" << endl;
        file = fopen(name1, "a");
        fprintf(file, "%lf\t%6.6e\t%6.6e\t%6.6e\t%6.6e\t%6.6e\t%6.6e\n", lambda, real(be01),imag(be01),real(be11),imag(be11),real(bo11),imag(bo11));
        fclose(file);
    }
}

double besselj_0(double rho)
{
    return sin(rho)/rho;
}

double besselj_1(double rho)
{
    return sin(rho)/pow(rho,2) - cos(rho)/rho;
    
}

double besselj_2(double rho)
{
    return ((3./pow(rho,2))-1.)*sin(rho)/rho - 3.*cos(rho)/pow(rho,2);
    
}


#endif
