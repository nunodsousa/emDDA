//
//  charger.h
//  eris
//
//  Created by Nuno de Sousa on 17/05/15.
//  Copyright (c) 2015 Nuno de Sousa. All rights reserved.
//

#ifndef eris_charger_h
#define eris_charger_h

#include "parameters.h"


using namespace std;

string input_parameter(string t);
string path(string t);

//    io_loader(&eps_0, &eps_m, &E0, &N_part, filename, &lambda_min, &lambda_max, &step, &mod, &ux, &uy, &uz, &xsource, &ysource, &zsource, &name_load_file, &output_string, &solver_method,&input_string_pos, &solver_Niterations, &solver_error);

void io_loader(double *eps_0, double *eps_m, double *E0, int *N_part, string filename, double *lambda_min, double *lambda_max, int *step, int *mod, double *ux, double *uy, double *uz, double *xsource, double *ysource, double *zsource, string *name_load_file, string *output_string, int *switch_projections, int *switch_forces, int *switch_sections, int *switch_fields, int *switch_polarizations, string *solver_method, string *input_field_pos, int *solver_Niterations, double *solver_error)
{
    ifstream myFile;
    myFile.open(filename.c_str());
    string line, loadpar;
    
    double number;
    int par = 0;
    int nn;
    
    
    
    if (myFile.is_open()) {
        
        while (!myFile.eof())
        {
            
            getline(myFile, line);
            
            if(! line.empty())
            {
                loadpar = input_parameter(line.c_str());
                
                if(par == 0)
                {
                    number = stof(loadpar);
                    *lambda_min = number;
                    cout << "lambda minimum = " << *lambda_min << endl;
                }
                
                
                if(par == 1)
                {
                    number = stof(loadpar);
                    *lambda_max = number;
                    cout << "lambda maximum = " << *lambda_max << endl;
                }
                
                if(par == 2)
                {
                    number = stoi(loadpar);
                    *step = number;
                    cout << "lambda steps = " << *step << endl;
                }
                
                if(par == 3)
                {
                    number = stof(loadpar);
                    *eps_0 = number;
                }
                
                
                if(par == 4)
                {
                    number = stof(loadpar);
                    *eps_m = number;
                }
                
                if(par == 5)
                {
                    number = stoi(loadpar);
                    *mod = number;
                }
                
                if(par == 6)
                {
                    number = stof(loadpar);
                    *ux = number;
                }
                
                if(par == 7)
                {
                    number = stof(loadpar);
                    *uy = number;
                }
                
                if(par == 8)
                {
                    number = stof(loadpar);
                    *uz = number;
                }
                
                if(par == 9)
                {
                    number = stof(loadpar);
                    *E0 = number;
                    //cout << "E0 entrance = " << number << endl;
                }
                
                if(par == 10)
                {
                    number = stof(loadpar);
                    *xsource = number;
                }
                
                if(par == 11)
                {
                    number = stof(loadpar);
                    *ysource = number;
                }
                
                
                if(par == 12)
                {
                    number = stof(loadpar);
                    *zsource = number;
                }
                
                
                if(par == 13)
                {
                    number = stoi(loadpar);
                    *N_part = number;
                }
                
                if(par == 14)
                {
                    *name_load_file = loadpar;
                }
                
                if(par == 15)
                {
                    *output_string = loadpar;
                }
                
                if(par == 16)
                {
                    *input_field_pos = loadpar;
                }
                
                if(par == 17)
                {
                    nn = stoi(loadpar);
                    *switch_sections = nn;
                }
                
                if(par == 18)
                {
                    nn = stoi(loadpar);
                    *switch_fields = nn;
                }
                
                if(par == 19)
                {
                    nn = stoi(loadpar);
                    *switch_polarizations = nn;
                }
                
                if(par == 20)
                {
                    nn = stoi(loadpar);
                    *switch_projections = nn;
                }
                
                if(par == 21)
                {
                    nn = stoi(loadpar);
                    *switch_forces = nn;
                }
                
                if(par == 22)
                {
                    *solver_method = loadpar;
                }
                
                if(par == 23)
                {
                    number = stoi(loadpar);
                    *solver_Niterations = number;
                }
                
                if(par == 24)
                {
                    number = stof(loadpar);
                    *solver_error = number;
                }
                
                par = par + 1;
            }
        }
        
    }
    
    myFile.close();
    
    
    
}


string input_parameter(string t)
{
    
    istringstream iss(t);
    string word;
    string special;
    int count = 0;
    while(iss >> word) {
        
        //cout << "word = " << word << endl;
        //cout << "iss = " << iss << endl;
        /* do stuff with word */
        
        if(count == 0)
        {
            special = word;
        }
        count ++;
    }
    
    return special;
    
}

string path(string t)
{
    
    istringstream iss(t);
    string word;
    string special;
    int count = 0;
    while(iss >> word) {
        
        //cout << "word = " << word << endl;
        //cout << "iss = " << iss << endl;
        /* do stuff with word */
        
        if(count == 1)
        {
            special = word;
        }
        count ++;
    }
    
    return special;
    
}


#endif
