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
#include <ctime>
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
#include "include/solver.h"
//#include "include/forces.h"
//#include "include/mie_polarizability/mie_polarizability.h"


using namespace std;
using namespace arma;
//using namespace sp_bessel; //necessary for complex_bessel

typedef complex<double> dcmplx;


int main(int argc, const char * argv[])
{
    cout << "Armadillo version: " << arma_version::as_string() << endl;
    
    //Clock operations
    clock_t start;
    double duration;
    start = clock();
    time_t prenow = time(0);

    
    
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
    string input_string_pos;
    int switch_projections;
    int switch_forces;
    int switch_sections;
    int switch_fields;
    int switch_polarizations;
    string loader_method = "static"; //can be static or dynamic
    
    Parameters param;
    
    io_loader(&eps_0, &eps_m, &E0, &N_part, filename, &lambda_min, &lambda_max, &step, &mod, &ux, &uy, &uz, &xsource, &ysource, &zsource, &name_load_file, &output_string, &switch_projections, &switch_forces, &switch_sections, &switch_fields, &switch_polarizations, &solver_method, &input_string_pos, &solver_Niterations, &solver_error);
    param.set_parameters_eps_0(eps_0);
    param.set_parameters_eps_m(eps_m);
    param.set_E0(E0);
    param.set_N_part(N_part);
    
    
    
    string filename_material = "material_list.txt";
    
    Material *mater;
    int n_materials = n_lines_material(filename_material);
    
    mater = new Material[n_materials];
    
    for(int i = 0; i < n_materials; i ++)
    {
        mater[i].loader(i, filename_material);
    }
    
    for(double lambda=lambda_min; lambda<lambda_max; lambda=lambda+(lambda_max-lambda_min)/double(step))
    {
        cout << "lambda = " << lambda << endl;
        Particle *part;
        part = new Particle[param.N_part];
        
        Equations eq;
        vec absorption(3);
        absorption.zeros();
        
        cx_mat epstemp;
        vec Power1(3),Power2(3);
        epstemp.eye(3,3);
        
        param.set_parameters_wavenumber(sqrt(param.eps_m)*2*M_PI/lambda);
        
        
        cx_mat pol (3,3);
        dcmplx eps_medium(param.eps_m,0);
        
        //Load file
        load_structure(param.N_part, part, name_load_file);
        
        for(int i = 0; i < param.N_part; i = i + 1)
        {
            //cout << "loading particle n.º = " << i << endl;
            int material_control = 0;
            
            for(int j = 0; j < n_materials; j = j + 1)
            {
                if(mater[j].material_nome() == part[i].material)
                {
                    part[i].set_dielectric_tensor(mater[j].return_interp_value(lambda));
                    part[i].set_sphere_polarizability(pow((3./(4.*M_PI))*part[i].volume,(1./3.)), param);
                    
                    cx_mat test = mater[j].return_interp_value(lambda);
                    material_control = 1;
                    
                }
                
                if("Meyer_gold" == part[i].material)
                {
                    cx_mat id33;
                    id33.eye(3,3);
                    part[i].set_dielectric_tensor(epsilon_gold(lambda)*id33);
                    part[i].set_sphere_polarizability(pow((3./(4.*M_PI))*part[i].volume,(1./3.)), param);
                    //cx_mat test = mater[j].return_interp_value(lambda);
                    
                    material_control = 1;
                }
                
                if("resonant" == part[i].material)
                {
                    part[i].set_resonant_polarizability(param);
                    
                    material_control = 1;
                }
            }
            
            if(material_control == 0)
            {
                cout << "MATERIAL DOESN'T EXIST IN THE LIST." << endl;
                exit(1);
            }
            
        }
        
        //transport
        //-------------
        //incident wave direction
        if(mod == 0)
        {
            
            param.set_plw_direction(ux, uy, uz);
            
            if(loader_method == "static")
            {
            eq.load_equations(param, part);
            }
            eq.load_ind_term_plw(param, part);
        }
        //-------------
        
        //emission
        //-------------
        if(mod ==1)
        {
            param.set_source_pos(xsource,ysource,zsource);
            eq.load_equations(param, part);
            eq.load_ind_term_source(param, part);
        }
        //-------------
        
        
        //Solver
        //-------------
        //Direct Solver
        if(solver_method == "Direct")
        {direct_solver(&eq, &eq.A, &eq.B);}
        //Biconjugate gradient stabilized method
        if(solver_method == "BiCStab")if(loader_method=="static")
            {{solver_BiCStab_static(&eq, &eq.A, &eq.B, param, solver_Niterations, solver_error);}}
        //Biconjugate gradient stabilized method with dynamic load
        if(solver_method == "BiCStab")if(loader_method=="dynamic")
        {{solver_BiCStab_dynamic(&eq, &eq.A, &eq.B, part, param, solver_Niterations, solver_error);}}
        //-------------
        //Solver Error Message
        if(solver_method != "Direct" && solver_method != "BiCStab" && loader_method != "static" && loader_method != "dynamic")
        {
            cout << "Error: The method or loader doesn't exist.\nProgram stopped." << endl;
            
            exit(1);
        }
        //-------------
        
        if(switch_polarizations == 1)
        {
        string namepolx = output_string + "_polarizations_xpol.dat";
        print_polarizations(eq, param, part, lambda, namepolx);
        string namepoly = output_string + "_polarizations_ypol.dat";
        print_polarizations(eq, param, part, lambda, namepoly);
        string namepolz = output_string + "_polarizations_zpol.dat";
        print_polarizations(eq, param, part, lambda, namepolz);
        }
        
        if(switch_fields == 1)
        {
            string namefieldx = output_string +"_fields_xpol.dat";
            print_fields(eq, param, part, lambda, namefieldx, "x");
            string namefieldy = output_string +"_fields_ypol.dat";
            print_fields(eq, param, part, lambda, namefieldy, "y");
            string namefieldz = output_string +"_fields_zpol.dat";
            print_fields(eq, param, part, lambda, namefieldz, "z");
        }
        
        
        //transport
        //----------
        if(mod == 0 && switch_sections == 1)
        {
            vec opt_theorem(3);
            vec sca_sec(3), ext_sec(3), abs_sec(3);
            
            abs_sec = abs_cross_section(&eq.X, param, part);
            ext_sec = ext_cross_section(&eq.X, param, part);
            sca_sec = sca_cross_section(&eq.X, &eq.im_parts, param, part);

            cout.precision(12);
            cout << "abs_cross_section = ";
            abs_sec.raw_print(cout);
            cout << "ext_cross_section = ";
            ext_sec.raw_print(cout);
            cout << "sca_cross_section = ";
            sca_sec.raw_print(cout);
            
            
            string sectionsx = output_string + "_sections_xpol.dat";
            print_sections(lambda, abs_sec(0), ext_sec(0), sca_sec(0),sectionsx);
            string sectionsy = output_string + "_sections_ypol.dat";
            print_sections(lambda, abs_sec(1), ext_sec(1), sca_sec(1),sectionsy);
            string sectionsz = output_string + "_sections_zpol.dat";
            print_sections(lambda, abs_sec(2), ext_sec(2), sca_sec(2),sectionsz);
            
            for(int j = 0; j < 3; j ++)
            {
                opt_theorem = ((abs_sec + sca_sec)-(ext_sec))/((abs_sec + sca_sec)+(ext_sec));
            }
            
            cout << "Optical theorem verification -> " << opt_theorem << endl;
        }
        
        //Projection over VSH basis
        if(switch_projections == 1)
        {
            dcmplx ae11, ao11, ae01, be11, bo11, be01;
            string name_aeo, name_beo;
            
            for(int pol = 0; pol < 3; pol = pol + 1)
            {
                ae11 = print_ae11(&eq.X, param, part, lambda, pol);
                ao11 = print_ao11(&eq.X, param, part, lambda, pol);
                ae01 = print_ae01(&eq.X, param, part, lambda, pol);
                be11 = print_be11(&eq.X, param, part, lambda, pol);
                bo11 = print_bo11(&eq.X, param, part, lambda, pol);
                be01 = print_be01(&eq.X, param, part, lambda, pol);
                
                if(pol == 0)
                {
                    name_aeo = output_string + "_projectionaeo_xpol.dat";
                    name_beo = output_string + "_projectionbeo_xpol.dat";
                }
                
                if(pol == 1)
                {
                    name_aeo = output_string + "_projectionaeo_ypol.dat";
                    name_beo = output_string + "_projectionbeo_ypol.dat";
                }
                
                if(pol == 2)
                {
                    name_aeo = output_string + "_projectionaeo_zpol.dat";
                    name_beo = output_string + "_projectionbeo_zpol.dat";
                }
                
                print_aeo(lambda, ae01, ae11, ao11, name_aeo);
                print_beo(lambda, be01, be11, bo11, name_beo);
            }
        }
        
        //Evaluation of the optical forces
        
        if(switch_forces == 1)
        {
            mat force(3.*param.N_part,3);
            bool numerical_verification = false; //compare with the numerical calculation of the green tensor
            //force = forces_eval(&eq.X, param, part, lambda, numerical_verification);
            //string forcex = output_string + "_force_xpol.dat";
            //print_forces(force, param, part, lambda,forcex, "x");
            //string forcey = output_string + "_force_ypol.dat";
            //print_forces(force, param, part, lambda,forcey, "y");
            //string forcez = output_string + "_force_zpol.dat";
            //print_forces(force, param, part, lambda,forcez, "z");
        }
        
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
        

        
        //Evaluation of the field in a set of points defined in a file
        
        full_field_evaluation(output_string, input_string_pos, eq, lambda, param, part);
        
        
        delete []part;
    }
    
    //Print running time
    //This is a closed function
    duration = ( clock() - start ) / (double) CLOCKS_PER_SEC;
    long duration2 = int(duration);
    int hour=duration2/3600.;
    int second=duration2 % 3600;
    int minute=second/60.;
    second %= 60;
    cout << "Computation time (s) = " << duration << "s. " << hour << "h" << minute << "min" << second << "sec" <<'\n';
    time_t now = time(0);
    duration2 = ( now - prenow);
    hour=duration2/3600.;
    second=duration2 % 3600;
    minute=second/60.;
    second %= 60;
    cout << "Lapsed time (s) = " << duration2 << "s. " << hour << "h" << minute << "min" << second << "sec" <<'\n';
    //End running time
    
    cout << "Program terminated." << endl;
    return 0;
}
