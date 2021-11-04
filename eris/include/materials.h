//
//  materials.h
//  eris
//
//  Created by Nuno de Sousa on 28/11/14.
//  Copyright (c) 2014 Nuno de Sousa. All rights reserved.


#ifndef eris_materials_h
#define eris_materials_h

#include <fstream>
#include <string>
#include <sstream>
#include <iostream>

int n_lines_material(string filename);


class Material{
    
    //read only
private:
    //cx_mat _dielectric_tensor;
    string material_name;
    string material_path;
    int n_lines;
    double **matrix;
    int word_counter;
    
    
    //Material():dielectric_tensor()//_dielectric_tensor)
    //{
    //}
public:
    void loader(int i, string filename)
    {
        
        //Start by opening the "material_list.txt" and read the i line.
        
        ifstream myFile;
        myFile.open(filename.c_str());
        string line;
        int id = 0;
        
        if (myFile.is_open())
        {
            while (!myFile.eof())
            {
                getline(myFile, line);
                if(! line.empty() && i == id)
                {
                    material_name = input_parameter(line.c_str());
                    //cout << "name = " << material_name << endl;
                    material_path = path(line.c_str());
                    //cout << "path = " << path(line.c_str()) << endl;
                }
                id = id + 1;
            }
        }
        
        myFile.close();
        
        //Count the number of lines of a specific material
        
        n_lines = n_lines_material(material_path);
        //cout << "N lines = " << int(n_lines) << endl;
        
        //Load the epsilon file to the eps table
        
        matrix = new double*[n_lines];
        for(int i=0; i < n_lines; i++)
            matrix[i]=new double[19];
        
        //eps_table *eps;
        //eps = new eps_table[int(n_lines)];
        
        myFile.open(material_path.c_str());
        
        //Words per line
        
        word_counter = 0;
        
        
        if (myFile.is_open())
        {
            
            getline(myFile, line);
            istringstream iss(line);
            
            if(!line.empty())
            {
                
                string word;
                while(iss >> word)
                {
                    word_counter = word_counter +1;
                }
            }
        }
        
        myFile.close();
        
        //cout << "There are " << word_counter << " words per line" << endl;
        
        
        //-------!!------
        
        int j = 0;
        int par = 0;
        
        myFile.open(material_path.c_str());
        
        if (myFile.is_open())
        {
            
            while (!myFile.eof())
            {
                
                getline(myFile, line);
                
                if(!line.empty())
                {
                    //If there are three elements present in each line:
                    //cout << "line = " << line << endl;
                    istringstream iss(line);
                    par = 0;
                    string word;
                    if(word_counter == 3)
                    {
                        while(iss >> word)
                        {
                            //cout << word << '\n';
                            
                            if(par == 0)
                            {
                                matrix[j][0] = stof(word);
                            }
                            
                            if(par == 1)
                            {
                                matrix[j][1] = stof(word);
                                matrix[j][9] = stof(word);
                                matrix[j][17] = stof(word);
                            }
                            
                            if(par == 2)
                            {
                                matrix[j][2] = stof(word);
                                matrix[j][10] = stof(word);
                                matrix[j][18] = stof(word);
                            }
                            
                            par = par + 1;
                        }
                        
                        j = j +1;
                    }
                    
                    par = 0;
                    if(word_counter == 19)
                    {
                        
                        while(iss >> word)
                        {
                            for(int m = 0; m < 20; m ++)
                            {
                                
                                if(par == m)
                                {
                                    matrix[j][m] = stof(word);
                                    //cout << "matrix[" << j << "][" << m << "] = " << matrix[j][m] << endl;
                                }
                            }
                            
                            par = par + 1;
                            
                        }
                        
                        j = j +1;
                    }
                }
            }
        }
        
        myFile.close();
        
        
        //for(int i = 0; i < n_lines; i ++)
        //{
            
        //    cout << "matrix[" << i << "][0] = " << matrix[i][0] << endl;
        //}
        
    }
    
    
    cx_mat return_interp_value(double lambda)
    {
        double dista;
        double epsr, epsi;
        cx_mat final_result(3,3);
        cx_mat id_mat;
        id_mat.eye(3,3);
        final_result.zeros();
        
        //cout << "Very special N_lines = " << n_lines << endl;
        //cout << "Very special eps = " << matrix[0][0] << endl;
        
        if(word_counter == 3)
        {
            
            for(int j = 0; j < n_lines-1; j ++)
            {
                //modifiquei a orentação dos sinais
                if(lambda >= matrix[j][0] && lambda <= matrix[j+1][0])
                {
                    //interpola
                    dista = matrix[j+1][0] - matrix[j][0];
                    epsr = ((( matrix[j+1][0] - lambda)/dista)*matrix[j][1] + ((lambda - matrix[j][0])/dista)*matrix[j+1][1]);
                    epsi = ((( matrix[j+1][0] - lambda)/dista)*matrix[j][2] + ((lambda - matrix[j][0])/dista)*matrix[j+1][2]);
                    
                    break;
                }
                
                //Warning
                if(lambda < matrix[0][0] | lambda > matrix[n_lines-1][0])
                {
                    cout << "Wavelenth = " << lambda << " out of limits"<< endl;
                    cout << "Minimum limit = " << matrix[0][0] << "; Maximum limit = " << matrix[n_lines-1][0] << endl;
                    exit(1);
                }
                
            }
            
            //cout << "espr = " << epsr << "   epsi = " << epsi << endl;
            
            dcmplx result(epsr, epsi);
            final_result = result*id_mat;
            
            
        }
        
        //cout << "aqui" << endl;
        
        if(word_counter == 19)
        {
            
            for(int j = 0; j < n_lines-1; j ++)
            {
                //modifiquei a orentação dos sinais
                if(lambda >= matrix[j][0] && lambda <= matrix[j+1][0])
                {
                    //interpola
                    dista = matrix[j+1][0] - matrix[j][0];
                    int l = 0, m = 0;
                    
                    for(int k = 0; k < 9; k ++)
                    {
                        
                        //cout << "aqui1.5" << endl;
                        epsr = ((( matrix[j+1][0] - lambda)/dista)*matrix[j][2*k+1] + ((lambda - matrix[j][0])/dista)*matrix[j+1][2*k+1]);
                        epsi = ((( matrix[j+1][0] - lambda)/dista)*matrix[j][2*k+2] + ((lambda - matrix[j][0])/dista)*matrix[j+1][2*k+2]);
                        
                        dcmplx bla(epsr, epsi);
                        
                        //cout << "bla = " << bla << endl;
                        
                        final_result(l,m) = bla;
                        
                        //cout << "gun" << endl;
                        
                        m = m + 1;
                        
                        if(m ==3)
                        {
                            m = 0;
                            l = l + 1;
                        }
                    }
                    break;
                    
                }
                
            }
            
            //Warning
            if(lambda < matrix[0][0] | lambda > matrix[n_lines-1][0])
            {
                cout << "Wavelenth = " << lambda << " out of limits"<< endl;
                cout << "Minimum limit = " << matrix[0][0] << "; Maximum limit = " << matrix[n_lines-1][0] << endl;
                exit(1);
            }
            
        }
        
        //cout << "aqui2" << endl;
        
        
        if(word_counter != 19 && word_counter != 3)
        {
            cout << "Wrong number of entrances in the eps file. " << endl;
            cout << "file = " << material_path << endl;
            cout << "PROGRAM STOPED" << endl;
            exit(1);
            
        }
        
        
        //cout << "matrix = " << final_result << endl;
        
        return final_result;
    }
    
    string material_nome(void)
    {
        return material_name;
    }
    
};


int n_lines_material(string filename)
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

//-----------------||--------------------

dcmplx epsilon_gold(double lambda)
{
    dcmplx I (0,1.);
    dcmplx eps;
    double eps_inf = 1.54;
    double lambda_p = 143;
    double gamma_p = 14500;
    double A_1 = 1.27;
    double phi_1 = -M_PI/4;
    double lambda_1 = 470;
    double gamma_1 = 1900;
    double A_2 = 1.1;
    double phi_2 = -M_PI/4;
    double lambda_2 = 325;
    double gamma_2 = 1060;
    
    
    if(lambda < 200 | lambda > 1000)
    {
        cout << "Wavelength out of range. It should be located between 200nm and 1000nm.\nProgram stoped." << endl;
        exit(0);
        
    }
    else
    {
        eps = eps_inf - 1./(pow(lambda_p,2)*((1./pow(lambda,2))+ I/(gamma_p*lambda))) + (A_1/lambda_1)*(exp(I*phi_1)/((1./lambda_1)-(1./lambda)-(I/gamma_1)) + exp(-I*phi_1)/((1./lambda_1)+(1./lambda)+(I/gamma_1))) + (A_2/lambda_2)*(exp(I*phi_2)/((1./lambda_2)-(1./lambda)-(I/gamma_2)) + exp(-I*phi_2)/((1./lambda_2)+(1./lambda)+(I/gamma_2)));
    }
    
    
    
    return eps;
}


/*dcmplx load_refrafiles(double lambda, string name)
{
    FILE *fp;
    
    double **material;
    
    //O numero de linhas devia ser automático
    
    ifstream f(name);
    char c;
    int i = 0;
    while (f.get(c))
        if (c == '\n')
        {
            ++i;
        }
    
    int nlines = i;
    
    cout << "nlines = " << nlines << endl;
    
    char name1[50];
    for(int i=0;i<50;i++)
        name1[i]=name[i];
    
    material = new double*[nlines];
    for(int j =0; j < nlines; j++)
        material[j] = new double[3];
    
    
    double l, n, k;
    double dista;
    
    double epsr = 0, epsi = 0;
    
    fp = fopen(name1, "r");
    
    for(int i = 0; i < nlines; i ++)
    {
        fscanf(fp, "%lf\t%lf\t%lf\n", &l, &n, &k);
        material[i][0] = l;
        material[i][1] =  pow(n,2)-pow(k,2);
        material[i][2] =  2*n*k;
        
    }
    
    for(int j = 0; j < nlines-1; j ++)
    {
        //modifiquei a orentação dos sinais
        if(lambda >= material[j][0] && lambda <= material[j+1][0])
        {
            //interpola
            dista = material[j+1][0] - material[j][0];
            epsr = ((( material[j+1][0] - lambda)/dista)*material[j][1] + ((lambda - material[j][0])/dista)*material[j+1][1]);
            epsi = ((( material[j+1][0] - lambda)/dista)*material[j][2] + ((lambda - material[j][0])/dista)*material[j+1][2]);
            
            break;
        }
        
        //Aumentar a descrição do warning.
        //Colocar
        if(lambda < material[0][0] | lambda > material[nlines-1][0])
        {
            cout << "Wavelenth = " << lambda << " out of limits"<< endl;
            cout << "Minimum limit = " << material[0][0] << "; Maximum limit = " << material[nlines-1][0] << endl;
            break;
        }
        
    }
    
    dcmplx eps(epsr,epsi);
    
    return eps;
    
}*/


#endif
