//
//  main.cpp
//  eris
//
//  Created by Nuno de Sousa on 27/12/13.
//  Copyright (c) 2013 Nuno de Sousa. All rights reserved.
//
//#include </usr/include/armadillo> //casa
//#include </usr/local/include/armadillo>
//
//

#define ARMA_64BIT_WORD
#include <iostream>
#include <complex>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <armadillo>
#include <cstdio>
#include <sstream>
#include <fstream>
#include <string>
#include "include/charger.h"
#include "include/particle.h"
#include "include/parameters.h"
#include "include/load_equations.h"
#include "include/power.h"
#include "include/sections.h"
#include "include/output.h"
#include "include/configurations.h"
#include "include/materials.h"
#include "include/fields.h"
#include "include/ab_coefficients.h"
#include <ctime>
#include "include/mie_polarizability/mie_polarizability.h"


using namespace std;
using namespace arma;
using namespace sp_bessel; //necessary for complex_bessel

typedef complex<double> dcmplx;


int main(int argc, const char * argv[])
{
    cout << "Armadillo version: " << arma_version::as_string() << endl;
    
    
    clock_t start;
    double duration;
    
    start = clock();
    
    // File to load the simulation parameter
    
    string filename = "simul_parameters.txt";
    
    // Load parameters of the simulation
    
    double lambda_min, lambda_max;
    int step, mod; //0 = transport, 1 = emission
    string name_load_file, output_string;
    double ux,uy,uz;
    double xsource, ysource, zsource;
    double eps_0, eps_m, E0;
    int N_part;
    string solver_method;
    int solver_Niterations;
    double solver_error;
    
    Parameters param;
    
    io_loader(&eps_0, &eps_m, &E0, &N_part, filename, &lambda_min, &lambda_max, &step, &mod, &ux, &uy, &uz, &xsource, &ysource, &zsource, &name_load_file, &output_string, &solver_method, &solver_Niterations, &solver_error);
    param.set_parameters_eps_0(eps_0);
    param.set_parameters_eps_m(eps_m);
    param.set_E0(E0);
    param.set_N_part(N_part);
    
    
    
    string filename_material = "material_list.txt";
    
    Material *mater;
    int n_materials = n_lines_material(filename_material);
    
    //cout << "Number of materials = " << n_materials << endl;
    
    mater = new Material[n_materials];
    
    for(int i = 0; i < n_materials; i ++)
    {
        mater[i].loader(i, filename_material);
    }
    
    for(double lambda=lambda_min; lambda<lambda_max; lambda=lambda+(lambda_max-lambda_min)/double(step))
    {
        cout << "lambda = " << lambda << endl;
        //double edge = 35.384615;
        Particle *part;
        part = new Particle[param.N_part];
        
        Equations eq;
        vec absorption(3);
        absorption.zeros();
        
        cx_mat epstemp;
        vec Power1(3),Power2(3);
        epstemp.eye(3,3);
        
        //dcmplx eps (real(epsilon_gold(lambda)),imag(epsilon_gold(lambda)));
        //cout << "Au_eps("<< lambda << ") = " << eps << endl;
        //dcmplx eps(8.0,0.);
        //string filemat = "dielectric_eps/Au_ext.dat";
        //dcmplx eps(real(load_refrafiles(lambda, filemat)),imag(load_refrafiles(lambda, filemat)));
        //cout << "lambda CO = " << load_refrafiles(lambda, filemat) << endl;
        //dcmplx mo(0.,0.1);
        
        
        //test to print the eps
        
        //string nameeps="eps_gold.dat";
        //print_eps_material(lambda, epsilon_gold(lambda), nameeps);
        
        //string nameeps1 = "eps_gold1.dat";
        //print_eps_material(lambda,load_refrafiles(lambda,filemat), nameeps1);
        
        //epstemp(0,0)=eps;
        //epstemp(1,1)=eps;
        //epstemp(2,2)=eps;
        //epstemp(0,1)=mo;
        //epstemp(1,0)=-mo;
        
        param.set_parameters_wavenumber(sqrt(param.eps_m)*2*M_PI/lambda);
        
        
        cx_mat pol (3,3);
        dcmplx eps_medium(param.eps_m,0);
        //cout << "polarization = " <<pol << endl;
        
        //Load file
        //load_sphere_pos(param.N_part, part, name_load_file);
        load_structure(param.N_part, part, name_load_file);
        
        for(int i = 0; i < param.N_part; i = i + 1)
        {
            //cout << "loading particle n.º = " << i << endl;
            
            for(int j = 0; j < n_materials; j = j + 1)
            {
                if(mater[j].material_nome() == part[i].material)
                {
                    part[i].set_dielectric_tensor(mater[j].return_interp_value(lambda));
                    part[i].set_sphere_polarizability(pow((3./(4.*M_PI))*part[i].volume,(1./3.)), param);
                    cx_mat test = mater[j].return_interp_value(lambda);
                    
                }
            }
            //part[i].set_dielectric_tensor(epstemp);
            //part[i].set_cube_polarizability(edge, param);
            //part[i].set_dielectric_tensor(epstemp);
            //part[i].set_cube_polarizability(edge, param);
            //part[i].set_cLDR_polarizability(edge,param);
            //part[i].set_CFD_polarizability(edge, param);
            //part[i].set_sphere_polarizability(20, param);
            //part[i].set_polarizability_tensor(Mie_Polarizability(lambda, eps_medium, eps, 20.3829));
            //cout << " polarization = " << part[i].polarizability_tensor << endl;
            
        }
        
        
        
        
        
        
        //string filemat = "Co_Palik.dat";
        //cout << "lambda CO = " << load_refrafiles(lambda, filemat) << endl;
        
        
        //generate_random_cluster(N_part, param, part, 2.5 , 54, 5, epstemp);
        
        //part[0].set_positions(0,0,40);
        //part[1].set_positions(0,0,80);
        //part[2].set_positions(0,80,0);
        //part[3].set_positions(0,0,40);
        //part[0].set_dielectric_tensor(epstemp);
        //part[1].set_dielectric_tensor(epstemp);
        //part[2].set_dielectric_tensor(epstemp);
        //part[3].set_dielectric_tensor(epstemp);
        //part[0].set_dielectric_const(eps);
        //part[1].set_dielectric_const(eps);
        //part[0].set_dielectric_tensor(epstemp);
        //part[1].set_dielectric_tensor(epstemp);
        //part[2].set_dielectric_tensor(epstemp);
        //part[3].set_dielectric_tensor(epstemp);
        //part[0].set_sphere_polarizability(20, param);
        //part[1].set_sphere_polarizability(20, param);
        //part[2].set_sphere_polarizability(20, param);
        //part[0].set_cube_polarizability(10, param);
        //part[1].set_cube_polarizability(10, param);
        //part[2].set_cube_polarizability(10, param);
        //part[3].set_cube_polarizability(10, param);
        //part[0].set_oblate_polarizability(50,5, param);
        //part[1].set_oblate_polarizability(50,5, param);
        //part[2].set_oblate_polarizability(50,5, param);
        //part[3].set_oblate_polarizability(50,5, param);
        //part[0].set_resonant_polarizability(param);
        //part[1].set_resonant_polarizability(param);
        //part[2].set_resonant_polarizability(param);
        
        
        
        //transport
        //----------
        //incident wave direction
        if(mod == 0)
        {
            
            param.set_plw_direction(ux, uy, uz);
            eq.load_equations(param, part);
            eq.load_ind_term_plw(param, part);
        }
        
        //emission
        //----------
        if(mod ==1)
        {
            param.set_source_pos(xsource,ysource,zsource);
            eq.load_equations(param, part);
            eq.load_ind_term_source(param, part);
        }
        
        //Solver
        //----------
        
        if(solver_method == "Direct")
        {
            cout << "Direct Solver init" << endl;
            eq.solver();
            cout << "Solver terminated" << endl;
        }
        //Solver using BiCStab
        
        if(solver_method == "BiCStab")
        {
            cout << "BiCStab Solver init" << endl;
            eq.solver_BiCStab(param, solver_Niterations, solver_error);
        }
        
        if(solver_method != "Direct" && solver_method != "BiCStab")
        {
            cout << "Error: The method doesn't exist" << endl;
            exit(1);
        }
        
        string namepolx = output_string + "_polarizations_x.dat";
        print_polarizations(eq, param, part, lambda, namepolx);
        string namepoly = output_string + "_polarizations_y.dat";
        print_polarizations(eq, param, part, lambda, namepoly);
        string namepolz = output_string + "_polarizations_z.dat";
        print_polarizations(eq, param, part, lambda, namepolz);
        
        //string namefield = "_fields";
        //print_fields(eq, param, lambda, namefield);
        
        //transport
        //----------
        if(mod == 0)
        {
            vec opt_theorem(3);
            vec sca_sec(3), ext_sec(3), abs_sec(3);
            
            abs_sec = abs_cross_section(&eq.X, param, part);
            ext_sec = ext_cross_section(&eq.X, param, part);
            sca_sec = sca_cross_section(&eq.X, &eq.im_parts, param, part);
            
            cout << "abs_cross_section = " << abs_sec << endl;
            cout << "ext_cross_section = " << ext_sec << endl;
            cout << "sca_cross_section = " << sca_sec << endl;
            
            
            string sectionsx = output_string + "_sections_x.dat";
            print_sections(lambda, abs_sec(0), ext_sec(0), sca_sec(0),sectionsx);
            string sectionsy = output_string + "_sections_y.dat";
            print_sections(lambda, abs_sec(1), ext_sec(1), sca_sec(1),sectionsy);
            string sectionsz = output_string + "_sections_z.dat";
            print_sections(lambda, abs_sec(2), ext_sec(2), sca_sec(2),sectionsz);
            
            for(int j = 0; j < 3; j ++)
            {
                opt_theorem = ((abs_sec + sca_sec)-(ext_sec))/((abs_sec + sca_sec)+(ext_sec));
            }
            
            cout << "Optical theorem verification -> " << opt_theorem << endl;
        }
        
        //Temporary projection over VSH basis
        
        dcmplx ae11, ao11, ae01, be11, bo11, be01;
        
        ae11 = print_ae11(&eq.X, param, part, lambda);
        ao11 = print_ao11(&eq.X, param, part, lambda);
        ae01 = print_ae01(&eq.X, param, part, lambda);
        be11 = print_be11(&eq.X, param, part, lambda);
        bo11 = print_bo11(&eq.X, param, part, lambda);
        be01 = print_be01(&eq.X, param, part, lambda);
        
        string name_aeo = output_string + "projectionaeo.dat";
        string name_beo = output_string + "projectionbeo.dat";
        
        print_aeo(lambda, ae01, ae11, ao11, name_aeo);
        print_beo(lambda, be01, be11, bo11, name_beo);
        
        
        
        //emission
        //----------
        if(mod == 1)
        {
            Power1 = power_method_1(&eq.X, param, part);
            Power2 = power_method_2(&eq.X, &eq.im_parts, param, part, &absorption);
            
            cout << "Absoption = " << absorption << endl;
            cout << "power method one = " << Power1 << endl;
            cout << "power method two = " << Power2 << endl;
            cout << "Optical theorem verification = " << (Power1-Power2)/(Power1+Power2);
            
            //Pass the name to param_sim
            string name_1 = output_string + "power1.dat";
            string name_2 = output_string + "power2.dat";
            string name_3 = output_string + "absorption.dat";
            print_Power(Power1, name_1);
            print_Power(Power2, name_2);
            print_Power(absorption, name_3);
            
        }
        //Evaluation of the far field
        
        vec point(3);
        cx_mat Efield(3,3);
        Efield.zeros(3,3);
        point(0) = 0;
        point(1) = 0;
        point(2) = 100000;
        
        Efield = E_field(point, eq, param, part);
        
        string namepfieldx = output_string + "_Epfield_xinc.dat";
        print_Efield_MO_point(lambda, point, Efield, namepfieldx, "x");
        string namepfieldy = output_string + "_Epfield_yinc.dat";
        print_Efield_MO_point(lambda, point, Efield, namepfieldy, "y");
        string namepfieldz = output_string + "_Epfield_zinc.dat";
        print_Efield_MO_point(lambda, point, Efield, namepfieldz, "z");
        
        
        
        delete []part;
    }
    
    duration = ( clock() - start ) / (double) CLOCKS_PER_SEC;
    
    std::cout<<"Duration of the execution: "<< duration <<'\n';
    
    cout << "Execution with success." << endl;
    return 0;
}
