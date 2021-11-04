//
//  configurations.h
//  eris
//
//  Created by Nuno de Sousa on 27/11/14.
//  Copyright (c) 2014 Nuno de Sousa. All rights reserved.
//

#ifndef eris_configurations_h
#define eris_configurations_h

#include <cstdio>

void load_structure(int N_part, Particle *part, string name)
{
    FILE *fp;
    
    char name1[50];
    for(int i=0;i<50;i++)
        name1[i]=name[i];
    
    double x_pos, y_pos, z_pos, volume;
    char stri [80];
    
    fp = fopen(name1, "r");
    
    if(fp != NULL){
        for(int i = 0; i < N_part; i ++)
        {
            //Não sei se devo por o \n ou não
            fscanf(fp, "%lf\t%lf\t%lf\t%lf\t%s\n", &x_pos, &y_pos, &z_pos, &volume, stri);
            part[i].set_positions(x_pos,y_pos,z_pos);
            part[i].set_volume(volume);
            string material = string(stri);
            part[i].set_material(material);
            
        }
    }
    
    else
    {
        cout << "Mierda" << endl;
        cout << "Missing structure file. Program terminated." << endl;
        exit(1);
    }
    
    fclose(fp);
    
    
    
    
}

void load_sphere_pos(int N_part, Particle *part,string name)
{
    FILE *fp;
    
    char name1[50];
    for(int i=0;i<50;i++)
        name1[i]=name[i];
    
    double x_pos, y_pos, z_pos;
    
    fp = fopen(name1, "r");
    
    for(int i = 0; i < N_part; i ++)
    {
        //Não sei se devo por o \n ou não
        fscanf(fp, "%lf\t%lf\t%lf\n", &x_pos, &y_pos, &z_pos);
        part[i].set_positions(x_pos,y_pos,z_pos);
        
        
    }
    
    fclose(fp);
    
}

void generate_random_cluster(int N_part, Parameters param, Particle *part, double particle_radius, double cluster_radius, double excl_vol_radius, cx_mat epstemp)
{
    int control;
    bool cond0 = false, cond1 = false, cond2 = false;
    double internal_radius = excl_vol_radius;
    double radius = cluster_radius;
    double x, y, z;
    
    for(int i = 0; i < N_part; i = i + 1)
    {
        control = 0;
        while(control == 0)
        {
            cond0 = true;
            cond1 = true;
            cond2 = true;
            //cout << "random = " << 10+20*(rand()/(static_cast<double>(RAND_MAX)+1.0)) << endl;
            x = (radius)*(2*(rand()/(static_cast<double>(RAND_MAX)+1.0))-1);
            y = (radius)*(2*(rand()/(static_cast<double>(RAND_MAX)+1.0))-1);
            z = (radius)*(2*(rand()/(static_cast<double>(RAND_MAX)+1.0))-1);
            
            if(sqrt(pow(x,2)+ pow(y,2) +pow(z,2)) < internal_radius)
            {
                cond0 = false;
            }
            
            
            if(sqrt(pow(x,2)+ pow(y,2) +pow(z,2)) >= (radius))
            {
                cond1= false;
            }
            
            for(int j = 0; j < i; j++)
            {
                if(sqrt(pow(x-part[j].x,2)+ pow(y-part[j].y,2) +pow(z-part[j].z,2)) < 3*particle_radius)
                {
                    cond2 = false;
                }
            }
            
            if(cond0 == true && cond1 == true && cond2 == true)
            {
                //cout << "pos = " << x << " " << y << " " << z << endl;
                part[i].set_positions(x,y,z);
                part[i].set_dielectric_tensor(epstemp);
                part[i].set_sphere_polarizability(particle_radius, param);
                control = 1;
            }
        }
    }
    
    
}

void generate_random_clusterv2(int N_part, Parameters param, Particle *part, double particle_radius, double cluster_radius, double excl_vol_radius, cx_mat epstemp)
{
    int control;
    bool cond0 = false, cond1 = false, cond2 = false;
    double internal_radius = excl_vol_radius;
    double radius = cluster_radius;
    double x, y, z, r, theta, phi;
    
    for(int i = 0; i < N_part; i = i + 1)
    {
        control = 0;
        while(control == 0)
        {
            cond0 = true;
            cond1 = true;
            cond2 = true;
            //cout << "random = " << 10+20*(rand()/(static_cast<double>(RAND_MAX)+1.0)) << endl;
            theta = 2*M_PI*(rand()/(static_cast<double>(RAND_MAX)+1.0));
            phi = acos(2*(rand()/(static_cast<double>(RAND_MAX)+1.0)) -1);
            r = (radius-internal_radius)*(rand()/(static_cast<double>(RAND_MAX)+1.0));
            
            x = r*cos(theta)*sin(phi);
            y = r*sin(phi)*sin(theta);
            z = r*cos(phi);
            
            if(r < internal_radius)
            {
                cond0 = false;
            }
            
            if(sqrt(pow(x,2)+ pow(y,2) +pow(z,2)) >= (radius))
            {
                cond1= false;
                exit(0);
            }
            
            for(int j = 0; j < i; j++)
            {
                if(sqrt(pow(x-part[j].x,2)+ pow(y-part[j].y,2) +pow(z-part[j].z,2)) < 3*particle_radius)
                {
                    cond2 = false;
                }
            }
            
            if(cond0 == true && cond1 == true && cond2 == true)
            {
                //cout << "pos = " << x << " " << y << " " << z << endl;
                part[i].set_positions(x,y,z);
                part[i].set_dielectric_tensor(epstemp);
                part[i].set_sphere_polarizability(particle_radius, param);
                control = 1;
            }
        }
    }
    
    
}



#endif
