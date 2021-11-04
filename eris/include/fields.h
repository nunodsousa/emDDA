//
//  fields.h
//  eris
//
//  Created by Nuno de Sousa on 09/12/14.
//  Copyright (c) 2014 Nuno de Sousa. All rights reserved.
//

#ifndef eris_fields_h
#define eris_fields_h

#include <stdlib.h>

int n_lines_field(string filename);


//Scattered field evaluation
cx_mat E_field(vec point, Equations eq, Parameters param, Particle *part)
{
    cx_mat E(3,3);
    E.zeros(3,3);
    
    //#pragma omp parallel
    {
        cx_mat E_private(3,3);
        E_private.zeros(3,3);
        
        //#pragma omp for
        for(int i = 0; i < param.N_part; i = i + 1)
        {
            cx_mat solution(3,3), p(3,3), G(3,3);
            vec vector_dist(3),vector_dist_norm(3);
            double dist;
            solution.zeros(3,3);
            p.zeros(3,3);
            
            vector_dist_norm.zeros(3);
            
            vector_dist(0)=point(0) - part[i].x;
            vector_dist(1)=point(1) - part[i].y;
            vector_dist(2)=point(2) - part[i].z;
            
            dist = norm(vector_dist,2);
            vector_dist_norm = vector_dist/dist;
            
            G = calc_green_normv2(vector_dist_norm, param, dist);
            
            for(int m = 0; m < 3; m = m + 1)
            {
                for(int n = 0; n < 3; n = n + 1)
                {
                    solution(m,n) = eq.X(3*i+m,n);
                }
            }
            
            p = part[i].polarizability_tensor*solution;
            
            E_private = E_private + (pow(param.wavenumber,2)/(param.eps_0*param.eps_m))*G*p;
        }
        
        //#pragma omp critical
        {
            E = E_private;
        }
    }
    
    return E;
}


//----------

void full_field_evaluation(string output_string, string input_string_pos, Equations eq, double lambda, Parameters param, Particle *part)
{
    //cout << "Number of lines = " << n_lines_field(input_string_pos) << endl;
    int i = n_lines_field(input_string_pos);
    
    ifstream myFile;
    
    
    string line;
    double x, y, z;
    string temp;
    int counter = 0;
    
    myFile.open(input_string_pos.c_str());
    
    if (myFile.is_open())
    {
        
        while (!myFile.eof() && counter < i)
        {
            
            getline(myFile, line);
            
            if(!line.empty())
            {
                //If there are three elements present in each line:
                //cout << "line = " << line << endl;
                istringstream iss(line);
                int par = 0;
                string word;
                
                while(iss >> word)
                {
                    //cout << word << '\n';
                    
                    if(par == 0)
                    {
                        x = stof(word);
                    }
                    
                    if(par == 1)
                    {
                        y = stof(word);
                        
                    }
                    
                    if(par == 2)
                    {
                        z = stof(word);
                    }
                    
                    par = par + 1;
                }
                
                //j = j +1;
                
                par = 0;
            }
            
            counter = counter + 1;
            
            if(i != 0)
            {
                vec point(3);
                cx_mat Efield(3,3);
                Efield.zeros(3,3);
                point(0) = x;
                point(1) = y;
                point(2) = z;
                
                Efield = E_field(point, eq, param, part);
                
                string namepfieldx = output_string + "_Epfield_xinc.dat";
                print_Efield_MO_point(lambda, point, Efield, namepfieldx, "x");
                string namepfieldy = output_string + "_Epfield_yinc.dat";
                print_Efield_MO_point(lambda, point, Efield, namepfieldy, "y");
                string namepfieldz = output_string + "_Epfield_zinc.dat";
                print_Efield_MO_point(lambda, point, Efield, namepfieldz, "z");
            }
            
        }
    }
    
    else
    {
        cout << "Fail in the Full Field Evaluation file" << endl;
        exit(1);
    }
    
    myFile.close();
    
    
}


int n_lines_field(string filename)
{
    ifstream myFile;
    myFile.open(filename.c_str());
    string line;
    
    int counter = 0;
    
    if (myFile.is_open())
    {
        
        while (!myFile.eof())
        {
            
            getline(myFile, line);
            
            if(! line.empty())
            {
                counter = counter + 1;
            }
        }
    }
    
    myFile.close();
    
    return counter;
    
}


#endif
