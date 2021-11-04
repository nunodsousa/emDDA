//
//  mie_polarizability.h
//  eris
//
//  Created by Nuno de Sousa on 03/03/15.
//  Copyright (c) 2015 Nuno de Sousa. All rights reserved.
//
//Complex bessel functions from here
//https://github.com/valandil/complex_bessel

#ifndef eris_mie_polarizability_h
#define eris_mie_polarizability_h

#include "bessel_functions/besselFunctions.h"

using namespace sp_bessel; //necessary for complex_bessel

dcmplx RicBesf(double l, dcmplx rho);
dcmplx RicBesfprime(double l, dcmplx rho);
dcmplx RicHanf(double l, dcmplx rho);
dcmplx RicHanfprime(double l, dcmplx rho);
dcmplx an(double l, dcmplx refmed, dcmplx refsca, dcmplx k_1, dcmplx k_2, double radius);
dcmplx qfunction(dcmplx k, double radius);
cx_mat Mie_Polarizability(double lambda, dcmplx refmed, dcmplx refsca, double radius);
/*
int main(void)
{
    cx_mat pol (3,3);
    double lambda = 1000;
    dcmplx refmed (1,0); //refraction index of the medium
    dcmplx eps (12.25,0.0); // refraction index of the scatterer
    double radius = 63; //radius of the particle in nanometers
    pol = Mie_Polarizability(lambda, refmed, eps, radius);
    
    cout << "polarization = " <<pol << endl;
}
*/

cx_mat Mie_Polarizability(double lambda, dcmplx refmed, dcmplx eps, double radius)
{
    cx_mat polarizability;
    polarizability.eye(3,3);
    dcmplx I (0,1);
    
    dcmplx refsca (sqrt(real(eps)),0);
    
    dcmplx lambda_med = lambda/(refmed); //wavelength in external medium in micrometers
    dcmplx lambda_sca = lambda/(refsca); //wavelength inside the sphere in micrometers
    
    //    dcmplx k_0 = 2*M_PI/lambda; //wavenumber in vaccum in micrometer^-1
    dcmplx k_1 = 2*M_PI/lambda_med; //wavenumber in external medium in micrometer^-1
    dcmplx k_2 = 2*M_PI/lambda_sca; //wavenumber inside the sphere in micrometer^-1
    
    polarizability = I*(6*M_PI)/(pow(k_1,3))*an(1, refmed, refsca, k_1, k_2, radius)* polarizability;
    
    return polarizability;
    
}

dcmplx an(double l, dcmplx refmed, dcmplx refsca, dcmplx k_1, dcmplx k_2, double radius)
{
    dcmplx temp;
    temp =(refsca*RicBesfprime(l,qfunction(k_1, radius))*RicBesf(l,qfunction(k_2, radius)) - refmed*RicBesf(l,qfunction(k_1, radius))*RicBesfprime(l,qfunction(k_2, radius)))/
    (refsca*RicHanfprime(l,qfunction(k_1, radius))*RicBesf(l,qfunction(k_2, radius)) - refmed*RicHanf(l,qfunction(k_1, radius))*RicBesfprime(l,qfunction(k_2, radius)));
    return temp;
    
}


dcmplx qfunction(dcmplx k, double radius)
{
    return k*radius;
}


//Ricatti-Bessel Function of the first kind
dcmplx RicBesf(double l, dcmplx rho)
{
    return rho*sqrt(M_PI/(2.*rho))*besselJ(l+0.5,rho);
}

//Partial derivative of Riccati-Bessel functions
dcmplx RicBesfprime(double l, dcmplx rho)
{
    dcmplx temp = sqrt(M_PI/(8.))*sqrt(1./rho)*besselJ(l+0.5,rho) + (sqrt(M_PI/2.)*(besselJ(l-0.5,rho)-besselJ(l+1.5,rho)))/(2.*sqrt(1./rho));
    return temp;
}

//Ricatti-Hankel function
dcmplx RicHanf(double l, dcmplx rho)
{
    dcmplx result;
    dcmplx I (0,1);
    dcmplx temp1;
    dcmplx temp2;
    temp1 = RicBesf(l, rho);
    temp2 = rho*sqrt(M_PI/(2.*rho))*besselY(l+0.5, rho);
    result = temp1+ I*temp2;
    return result;
}

//Partial Derivative of Ricatti-Hankel Function
dcmplx RicHanfprime(double l, dcmplx rho)
{
    dcmplx temp1, temp2;
    dcmplx result;
    dcmplx I (0,1);
    temp1 = sqrt(M_PI/(8.*rho))*besselJ(l+0.5,rho) + (sqrt(M_PI/2)*(besselJ(l-0.5,rho)-besselJ(l+1.5,rho)))/(2.*sqrt(1./rho));
    temp2 = I*sqrt(M_PI/(8.*rho))*besselY(l+0.5,rho) + I*(sqrt(M_PI/2)*(besselY(l-0.5,rho)-besselY(l+1.5,rho)))/(2.*sqrt(1./rho));
    result = temp1 + temp2;
    return result;
}



//End of the program

#endif
