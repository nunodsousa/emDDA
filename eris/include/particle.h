//
//  particle.h
//  eris
//
//  Created by Nuno de Sousa on 28/12/13.
//  Copyright (c) 2013 Nuno de Sousa. All rights reserved.
//

#ifndef eris_particle_h
#define eris_particle_h

#include "parameters.h"

typedef std::complex<double> dcmplx;
using namespace arma;


class Particle{

private:
    bool dielectric_controller;

//read only
private:
    
    double _x, _y, _z; // Note: different name than public, read-only interface

    double _volume;
    
    string _material;
    
    cx_mat _dielectric_tensor;
    
    cx_mat _polarizability_tensor;
    
public:
    
    Particle():x(_x), y(_y), z(_z), volume(_volume), dielectric_tensor(_dielectric_tensor), material(_material), polarizability_tensor(_polarizability_tensor)
    {
        dielectric_controller = false;
        
        _x = 0;
        _y = 0;
        _z = 0;
        
        _volume = 0;
        
        _dielectric_tensor.zeros(3,3);
    
        _polarizability_tensor.zeros(3,3);
    }
    
    void set_positions(double x, double y, double z)
    {
        this->_x = x;
        this->_y = y;
        this->_z = z;
    }
    
    void set_volume(double volume)
    {
        this->_volume = volume;
        
    }
    
    void set_material(string material)
    {
        
        this ->_material = material;
    }
    
    void set_dielectric_const(dcmplx eps){
        
        dielectric_controller = true;
        cx_mat temporary;
        temporary.zeros(3,3);
        temporary = eps*temporary.eye();
        
        _dielectric_tensor=temporary;
        
        this->_dielectric_tensor = dielectric_tensor;
    }
    
    void set_dielectric_tensor(cx_mat dielectric_tensor){
        
        dielectric_controller = true;
        this->_dielectric_tensor=dielectric_tensor;
        
    }
    
    void set_cube_polarizability(double edge, Parameters param)
    {
        if(dielectric_controller == false)
        {
            cout << "Dielectric tensor not loaded. Program stoped." << endl;
            exit(0);
        }
        cx_mat alpha_0, alpha;
        cx_mat identity(3,3);
        double volume;
        identity.eye();
        dcmplx radiative_factor (0,(pow(param.wavenumber,3)/(6*M_PI*param.eps_0*param.eps_m)));
        
        volume = pow(edge,3);
        
        alpha_0 = (3*volume*param.eps_0*param.eps_m)*(_dielectric_tensor-param.eps_m*identity)*inv(_dielectric_tensor+param.eps_m*2.*identity);
        
        alpha = alpha_0*inv(identity-radiative_factor*alpha_0);
        
        _polarizability_tensor = alpha;
        
        this-> _polarizability_tensor = polarizability_tensor;
        this-> _volume = volume;
        
    }
    
    void set_cLDR_polarizability(double edge, Parameters param)
    {
        if(dielectric_controller == false)
        {
            cout << "Dielectric tensor not loaded. Program stoped." << endl;
            exit(0);
        }
        
        double b1 = -1.8915316;
        double b2 = 0.1648469;
        double b3 = -1.7700004;
        double S = 0.25;
        
        cx_mat alpha_0, alpha;
        cx_mat identity(3,3);
        double volume;
        identity.eye();
        dcmplx radiative_factor (0,(pow(param.wavenumber,3)/(6*M_PI*param.eps_0*param.eps_m)));
        
        dcmplx temp_eps = _dielectric_tensor(0,0);
        
        volume = pow(edge,3);
        
        alpha_0 = (3*volume*param.eps_0*param.eps_m)*(_dielectric_tensor-param.eps_m*identity)*inv(_dielectric_tensor+param.eps_m*2.*identity);
        
        dcmplx eta = (alpha_0(0,0)/(4*M_PI*edge))*pow(param.wavenumber,2)*(b1+temp_eps*b2+temp_eps*b3*S);
        
        alpha = alpha_0*inv(identity + eta*identity - radiative_factor*alpha_0);
        
        _polarizability_tensor = alpha;
        
        this-> _polarizability_tensor = polarizability_tensor;
        this-> _volume = volume;
        
    }
    
    
    void set_CFD_polarizability(double edge, Parameters param)
    {
        if(dielectric_controller == false)
        {
            cout << "Dielectric tensor not loaded. Program stoped." << endl;
            exit(0);
        }
        

        cx_mat alpha_0, alpha;
        cx_mat identity(3,3);
        double volume;
        dcmplx cp (0,1);
        identity.eye();
        //dcmplx radiative_factor (0,(pow(param.wavenumber,3)/(6*M_PI*param.eps_0*param.eps_m)));
        
        volume = pow(edge,3);
        
        alpha_0 = (3*volume*param.eps_0*param.eps_m)*(_dielectric_tensor-param.eps_m*identity)*inv(_dielectric_tensor+param.eps_m*2.*identity);
        
        dcmplx eta = (1/(3*M_PI*edge))*pow(param.wavenumber,2);
        dcmplx phi = (pow(param.wavenumber,3)/(6*M_PI))*(cp + log((M_PI-param.wavenumber*edge)/(M_PI+param.wavenumber*edge))/M_PI);
        
        alpha = alpha_0*inv(identity - eta*identity - phi*alpha_0);
        
        _polarizability_tensor = alpha;
        
        this-> _polarizability_tensor = polarizability_tensor;
        this-> _volume = volume;
        
    }
    
    void set_Mie_polarizability()
    {
        
        
        
    }

    
    // Under development. This is the polarizability defined by Vincenzo in his paper.
    // Not verified for the moment
    // It must be called in the following way:
    // part[i].set_cubic_VG_polarizability(pow(part[i].volume,(1./3.)), param);
    void set_cubic_VG_polarizability(double edge, Parameters param)
    {
        
        if(dielectric_controller == false)
        {
            cout << "Dielectric tensor not loaded. Program stoped." << endl;
            exit(0);
        }
        
        dcmplx Im(0,1.);
        cx_mat iden(3,3);
        iden.eye();
        
        dcmplx epsilon = _dielectric_tensor(0,0);
        
        double volume = pow(edge,3);
        double beta = 12.6937*pow(edge,2.);
        double omega = 2.*M_PI/3.;
        dcmplx delta = (8.*volume/(pow(3.,(3./2.))*pow(volume,1./2.)))*(1./epsilon);
        
        
        dcmplx alpha = 8.*volume/(1./(epsilon-1.)-(1./(4.*M_PI))*(-2.*omega-delta+pow(param.wavenumber,2)*beta/2.+(16./3.)*Im*pow(param.wavenumber,3)*volume));
        
        _polarizability_tensor = alpha * iden;
        
        this-> _polarizability_tensor = polarizability_tensor;
        this-> _volume = volume;
    }
    
    
    void set_sphere_polarizability(double radius, Parameters param)
    {
        if(dielectric_controller == false)
        {
            cout << "Dielectric tensor not loaded. Program stoped." << endl;
            exit(0);
        }
        cx_mat alpha_0, alpha;
        cx_mat identity(3,3);
        double volume;
        identity.eye();
        dcmplx radiative_factor (0,(pow(param.wavenumber,3)/(6*M_PI*param.eps_0*param.eps_m)));
        
        volume = (4./3.)*M_PI*pow(radius,3);
        
        alpha_0 = (3*volume*param.eps_0*param.eps_m)*(_dielectric_tensor-param.eps_m*identity)*inv(_dielectric_tensor+param.eps_m*2.*identity);
                
        alpha = alpha_0*inv(identity-radiative_factor*alpha_0);
        
        _polarizability_tensor = alpha;
        
        
        
        this-> _polarizability_tensor = polarizability_tensor;
        this-> _volume = volume;
        
        
    }
    
    void set_oblate_polarizability(double a, double c, Parameters param)
    {
        if(dielectric_controller == false)
        {
            cout << "Dielectric tensor not loaded. Program stoped." << endl;
            exit(0);
        }
        //To be implemented
        
        double e,g,L1;
        double volume;
        cx_mat alpha_0, alpha;
        cx_mat identity(3,3);
        identity.eye();
        dcmplx radiative_factor (0,(pow(param.wavenumber,3)/(6*M_PI*param.eps_0*param.eps_m)));
        
        mat L(3,3);
        L.zeros(3,3);
        
        
        e = sqrt(1-pow(c,2)/pow(a,2));
        g=sqrt((1-pow(e,2))/pow(e,2));
        L1 = g/(2*pow(e,2))*((M_PI)/2 - atan(g))-pow(g,2)/2;
        
        L(0,0) = L1;
        L(1,1) = L1;
        L(2,2) = 1-2*L1;
        
        //cout << "L = " << L << endl;
        
        volume = (4./3.)*M_PI*pow(a,2)*c;
        
        alpha_0 = (3*volume*param.eps_0*param.eps_m)*(_dielectric_tensor-param.eps_m*identity)*inv(3*identity +3*L*(_dielectric_tensor - param.eps_m*identity));
        
        alpha = alpha_0*inv(identity-radiative_factor*alpha_0);
        
        //cout << "alpha = " << alpha << endl;
        
        _polarizability_tensor = alpha;
        
        this-> _polarizability_tensor = polarizability_tensor;
        this-> _volume = volume;
        
        
        //cout << "Not implemented. Program stoped." << endl;
        //exit(0);
        
    }
    
    void set_prolate_polarizability(double a, double c)
    {
        if(dielectric_controller == false)
        {
            cout << "Dielectric tensor not loaded. Program stoped." << endl;
            exit(0);
        }
        //To be implemented
        cout << "Not implemented. Program stoped." << endl;
        exit(0);
        
    }
    
    void set_polarizability_tensor(cx_mat polarizability_tensor)
    {
        
        this->_polarizability_tensor=polarizability_tensor;
        
    }
    
    void set_resonant_polarizability(Parameters param)
    {
        cx_mat alpha;
        alpha.eye(3,3);
        dcmplx factor (0,((6*M_PI*param.eps_0*param.eps_m)/pow(param.wavenumber,3)));
        double volume = 1;
        
        alpha = factor*alpha;
        
        _polarizability_tensor = alpha;
                
        this-> _polarizability_tensor = polarizability_tensor;
        this-> _volume = volume;
        
    }

    const double &x;
    const double &y;
    const double &z;
    const double &volume;
    const cx_mat &dielectric_tensor;
    const cx_mat &polarizability_tensor;
    const string &material;
};

#endif
