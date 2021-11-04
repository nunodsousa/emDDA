//
//  output.h
//  eris
//
//  Created by Nuno de Sousa on 13/11/14.
//  Copyright (c) 2014 Nuno de Sousa. All rights reserved.
//

#ifndef eris_output_h
#define eris_output_h

#include <cstdio>
#include <cstdlib>

bool file_exists(const char * filename)
{
    if (FILE * file = fopen(filename, "r"))
    {
        fclose(file);
        return true;
    }
    return false;
}


//Print the polarization of each particle
void print_polarizations(Equations eq, Parameters param, Particle *part, double lambda, string name)
{
    FILE *file;
    cx_mat p1(3,3);
    cx_mat solution(3,3);
    solution.zeros(3,3);
    
    
    char name1[100];
    for(int i=0;i<100;i++)
        name1[i]=name[i];
    
    
    p1.zeros();
    
    for(int i = 0; i < param.N_part; i = i +1)
    {
        p1.zeros();
        
        for(int l = 0; l < 3; l = l +1)
        {
            for(int m = 0; m < 3; m = m + 1)
            {
                solution(l,m) = eq.X(3*i+l,m);
            }
        }
        
        p1 = part[i].polarizability_tensor* solution;
        
        
        if(file_exists(name1) == 0) //if doesn't exist
        {
            file = fopen(name1, "w");
            fprintf(file,"#lambda\tn\tx\ty\tz\tre{pn}_x\tim{pn}_x\tre{pn}_y\tim{pn}_y\tre{pn}_z\tim{pn}_z\n");
            fclose(file);
        }
        
        if(file_exists(name1) != 0)
        {
            file = fopen(name1, "a");
            
            fprintf(file, "%6.7e\t%d\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\n", lambda, i, part[i].x, part[i].y, part[i].z, real(p1(0,0)), imag(p1(0,0)),real(p1(1,0)), imag(p1(1,0)),real(p1(2,0)), imag(p1(2,0)));
            fclose(file);
            
        }
        
    }
}

//Used for the Problem of the dielectric magneto-optic sphere.
//Not included in the main code
void print_field_integralv(Equations eq, Parameters param, Particle *part, double lambda, string name, double edge)
{
    FILE *file;
    cx_mat p1(3,3);
    cx_mat solution(3,3);
    cx_mat Esquare(3,3);
    solution.zeros(3,3);
    Esquare.zeros(3,3);
    
    cout << name << endl;
    
    char name1[50];
    for(int i=0;i<50;i++)
        name1[i]=name[i];
    
    file = fopen(name1, "a");
    
    p1.zeros();
    
    for(int i = 0; i < param.N_part; i = i +1)
    {
        p1.zeros();
        
        for(int l = 0; l < 3; l = l +1)
        {
            for(int m = 0; m < 3; m = m + 1)
            {
                solution(l,m) = eq.X(3*i+l,m);
            }
        }
        
        
        Esquare = Esquare + (0.5)*pow(edge,3)*solution.t()*solution;
        
        //cout << "Esquare = " << Esquare << endl;
        //cout << "Solution = " << solution << endl;
        //cout << "Calc = " << solution(0,0)*conj(solution(0,0))+ solution(1,0)*conj(solution(1,0))+ solution(2,0)*conj(solution(2,0)) << endl << endl << endl;
        
        //cout << "another simultion" << endl;
        
    }
    
    fprintf(file, "%6.7e\t%6.7e\n", lambda, real(Esquare(0,0)));
    
    fclose(file);
}


//Print the structure under study
void print_structure(Parameters param, Particle *part, string name)
{
    
    FILE *file;
    
    //cout << name << endl;
    
    char name1[50];
    for(int i=0;i<50;i++)
        name1[i]=name[i];
    
    file = fopen(name1, "w");
    
    for(int i = 0; i < param.N_part; i = i +1)
    {
        fprintf(file, "%6.7e\t%6.7e\t%6.7e\n", part[i].x, part[i].y, part[i].z);
    }
    fclose(file);
    
}


//Used in the main code
//It prints the field experienced by each particle that forms the system.
void print_fields(Equations eq, Parameters param, Particle *part, double lambda, string name, string inc)
{
    FILE *file;
    cx_mat solution(3,3);
    solution.zeros(3,3);
    
    //cout << name << endl;
    
    char name1[50];
    for(int i=0;i<50;i++)
        name1[i]=name[i];
    
    if(file_exists(name1) == 0) //if doesn't exist
    {
        file = fopen(name1, "w");
        fprintf(file,"#lambda\tn\tx\ty\tz\tre{En}_x\tim{En}_x\tre{En}_y\tim{En}_y\tre{En}_z\tim{En}_z\n");
        fclose(file);
    }
    
    if(file_exists(name1) != 0)
    {
        file = fopen(name1, "a");
        
        for(int i = 0; i < param.N_part; i = i +1)
        {
            
            for(int l = 0; l < 3; l = l +1)
            {
                for(int m = 0; m < 3; m = m + 1)
                {
                    solution(l,m) = eq.X(3*i+l,m);
                }
            }
            if(inc == "x")
            fprintf(file, "%6.7e\t%d\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\n", lambda, i, part[i].x, part[i].y, part[i].z, real(solution(0,0)), imag(solution(0,0)),real(solution(1,0)), imag(solution(1,0)),real(solution(2,0)), imag(solution(2,0)));
            if(inc == "y")
            fprintf(file, "%6.7e\t%d\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\n", lambda, i, part[i].x, part[i].y, part[i].z, real(solution(0,1)), imag(solution(0,1)),real(solution(1,1)), imag(solution(1,1)),real(solution(2,1)), imag(solution(2,1)));
            if(inc == "z")
            fprintf(file, "%6.7e\t%d\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\n", lambda, i, part[i].x, part[i].y, part[i].z, real(solution(0,2)), imag(solution(0,2)),real(solution(1,2)), imag(solution(1,2)),real(solution(2,2)), imag(solution(2,2)));
        }
        
        fclose(file);
        
    }
    
    

}


//Used in the main code
//It prints the field experienced by each particle that forms the system.
void print_forces(mat forces, Parameters param, Particle *part, double lambda, string name, string inc)
{
    FILE *file;
    mat forces_bp;
    forces_bp.zeros(3, 3);
    
    char name1[50];
    for(int i=0;i<50;i++)
        name1[i]=name[i];
    
    if(file_exists(name1) == 0) //if doesn't exist
    {
        file = fopen(name1, "w");
        fprintf(file,"#lambda\tid_particle\tx\ty\tz\tF_x\tF_y\tF_z\n");
        fclose(file);
    }
    
    if(file_exists(name1) != 0)
    {
        //cout << "entra em ficheiro existente" << endl;
        file = fopen(name1, "a");
        
        for(int i = 0; i < param.N_part; i = i +1)
        {
            
            for(int l = 0; l < 3; l = l +1)
            {
                for(int m = 0; m < 3; m = m + 1)
                {
                    forces_bp(l,m) = forces(3*i+l,m);
                }
            }
            if(inc == "x")
                fprintf(file, "%6.7e\t%d\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\n", lambda, i, part[i].x, part[i].y, part[i].z, forces_bp(0,0), forces_bp(1,0), forces_bp(2,0));
            if(inc == "y")
                fprintf(file, "%6.7e\t%d\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\n", lambda, i, part[i].x, part[i].y, part[i].z, forces_bp(0,1), forces_bp(1,1), forces_bp(2,1));
            if(inc == "z")
                fprintf(file, "%6.7e\t%d\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\n", lambda, i, part[i].x, part[i].y, part[i].z, forces_bp(0,2), forces_bp(1,2), forces_bp(2,2));
        }
        
        fclose(file);
    }
    
}


//Used in the main code
//It prints the scattering, absoption and extintion cross section
void print_sections(double lambda, double abs_sec, double ext_sec, double sca_sec, string name)
{
    
    FILE *file;
    
    //cout << name << endl;
    
    char name1[50];
    for(int i=0;i<50;i++)
        name1[i]=name[i];
    
    
    if(file_exists(name1) == 0) //if doesn't exist
    {
        file = fopen(name1, "w");
        fprintf(file,"#lambda\tk\tabs_sec\text_sec\tscat_sec\n");
        fclose(file);
    }
    
    if(file_exists(name1) != 0)
    {
        //cout << "entra em ficheiro existente" << endl;
        file = fopen(name1, "a");
        fprintf(file, "%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\n", lambda, 2*M_PI/lambda, abs_sec, ext_sec, sca_sec);
        fclose(file);
    }
}

//not used in the main code
//It print the epsilon of the material. It must be given the lambda and a complex number.
//For future use must be improved to create a file if doesn't exist, etc, etc, etc.
void print_eps_material(double lambda, dcmplx eps, string name)
{
    FILE *file;
    
    cout << name << endl;
    
    char name1[50];
    for(int i=0;i<50;i++)
        name1[i]=name[i];
    
    file = fopen(name1, "a");
    
    fprintf(file,"%6.7e\t%6.7e\t%6.7e\n", lambda, real(eps), imag(eps));
    fclose(file);
    
}

//Used in the main code
//To print the radiated, absorved and total Power
//for a system with emitters
void print_Power(vec Power, string name)
{
    
    FILE *file;
    
    //cout << name << endl;
    
    char name1[50];
    for(int i=0;i<50;i++)
        name1[i]=name[i];
    
    if(file_exists(name1) == 0) //if doesn't exit
    {
        file = fopen(name1, "w");
        fprintf(file,"#x\ty\tz\n");
        fclose(file);
    }
    
    if(file_exists(name1) != 0)
    {
        //cout << "entra em ficheiro existente" << endl;
        file = fopen(name1, "a");
        fprintf(file, "%6.7e\t%6.7e\t%6.7e\n", Power(0), Power(1), Power(2));
        fclose(file);
    }
}

//Used in the main code, in the function full_field_evaluation
//Print the electromagnetic field in one point, with the rotation and elliticity
void print_Efield_MO_point(double lambda, vec point, cx_mat Efield, string name, string inc)
{
    FILE *file;
    
    char name1[50];
    for(int i=0;i<50;i++)
        name1[i]=name[i];
    
    if(inc == "x")
    {
        
        if(file_exists(name1) == 0) //if doesn't exit
        {
            file = fopen(name1, "w");
            fprintf(file,"#lambda\tx\ty\tz\tExr\tExi\tEyr\tEyi\tEzr\tEzi\trot\tellip\n");
            fclose(file);
        }
        
        
        double rot = (real(Efield(0,0))*real(Efield(1,0))+imag(Efield(0,0))*imag(Efield(1,0)))/(pow(real(Efield(0,0)),2)+pow(imag(Efield(0,0)),2));
        double phi = (real(Efield(0,0))*imag(Efield(1,0))-real(Efield(1,0))*imag(Efield(0,0)))/(pow(real(Efield(0,0)),2)+pow(imag(Efield(0,0)),2));
        
        
        if(file_exists(name1) != 0)
        {
            //cout << "entra em ficheiro existente" << endl;
            file = fopen(name1, "a");
            fprintf(file, "%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\n", lambda, point(0), point(1), point(2), real(Efield(0,0)), imag(Efield(0,0)), real(Efield(1,0)), imag(Efield(1,0)), real(Efield(2,0)) ,imag(Efield(2,0)),rot, phi);
            fclose(file);
        }
    }
    
    if(inc == "y")
    {
        
        if(file_exists(name1) == 0) //if doesn't exit
        {
            file = fopen(name1, "w");
            fprintf(file,"#lambda\tx\ty\tz\tExr\tExi\tEyr\tEyi\tEzr\tEzi\trot\tellip\n");
            fclose(file);
        }
        
        
        double rot = (real(Efield(0,1))*real(Efield(1,1))+imag(Efield(0,1))*imag(Efield(1,1)))/(pow(real(Efield(0,1)),2)+pow(imag(Efield(0,1)),2));
        double phi = (real(Efield(0,1))*imag(Efield(1,1))-real(Efield(1,1))*imag(Efield(0,1)))/(pow(real(Efield(0,1)),2)+pow(imag(Efield(0,1)),2));
        
        
        if(file_exists(name1) != 0)
        {
            //cout << "entra em ficheiro existente" << endl;
            file = fopen(name1, "a");
            fprintf(file, "%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\n", lambda, point(0), point(1), point(2), real(Efield(0,1)), imag(Efield(0,1)), real(Efield(1,1)), imag(Efield(1,1)), real(Efield(2,1)) ,imag(Efield(2,1)),rot, phi);
            fclose(file);
        }
    }
    
    if(inc == "z")
    {
        
        if(file_exists(name1) == 0) //if doesn't exit
        {
            file = fopen(name1, "w");
            fprintf(file,"#lambda\tx\ty\tz\tExr\tExi\tEyr\tEyi\tEzr\tEzi\trot\tellip\n");
            fclose(file);
        }
        
        
        double rot = (real(Efield(0,2))*real(Efield(1,2))+imag(Efield(0,2))*imag(Efield(1,2)))/(pow(real(Efield(0,2)),2)+pow(imag(Efield(0,2)),2));
        double phi = (real(Efield(0,2))*imag(Efield(1,2))-real(Efield(1,2))*imag(Efield(0,2)))/(pow(real(Efield(0,2)),2)+pow(imag(Efield(0,2)),2));
        
        
        if(file_exists(name1) != 0)
        {
            //cout << "entra em ficheiro existente" << endl;
            file = fopen(name1, "a");
            fprintf(file, "%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\n", lambda, point(0), point(1), point(2), real(Efield(0,2)), imag(Efield(0,2)), real(Efield(1,2)), imag(Efield(1,2)), real(Efield(2,2)) ,imag(Efield(2,2)),rot, phi);
            fclose(file);
        }
    }
}

void print_Efield_part(double lambda, Equations eq, Parameters param, Particle *part, string name)
{
    
    cx_mat solution(3,3);
    solution.zeros(3,3);
    FILE *file;
    
    char name1[50];
    for(int i=0;i<50;i++)
        name1[i]=name[i];
    
    file = fopen(name1, "w");
    
    fprintf(file,"#l\ti\tx\ty\tz\trEx\tiEx\trEy\tiEy\trEz\tiEz\n");
    
    for(int i = 0; i < param.N_part; i = i +1)
    {
        for(int l = 0; l < 3; l = l +1)
        {
            for(int m = 0; m < 3; m = m + 1)
            {
                solution(l,m) = eq.X(3*i+l,m);
            }
        }
        
        fprintf(file,"%6.7e\t%i\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\t%6.7e\n", lambda, i, part[i].x, part[i].y, part[i].z, real(solution(0,0)), imag(solution(0,0)), real(solution(1,0)), imag(solution(1,0)), real(solution(2,0)), imag(solution(2,0)));
        
    }
    
    fclose(file);
}



#endif

