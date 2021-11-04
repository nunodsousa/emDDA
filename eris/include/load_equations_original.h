void load_equations(Parameters param, Particle *part)
{
    cx_mat temporary; //I only use to add a diagonal identity
    _A.zeros(3*param.N_part,3*param.N_part);
    _im_parts.zeros(3*param.N_part,3*param.N_part);
    vec vector_dist (3);
    vec vector_dist_norm = zeros<vec>(3);
    double dist;
    cx_mat output_tensor(3,3);
    
    
    for(int i = 0; i < param.N_part; i = i + 1)
    {
        for(int j = 0; j < param.N_part; j = j + 1)
        {
            if(i != j)
            {
                vector_dist(0)=part[i].x - part[j].x;
                vector_dist(1)=part[i].y - part[j].y;
                vector_dist(2)=part[i].z - part[j].z;
                
                dist = norm(vector_dist,2);
                vector_dist_norm = vector_dist/dist;
                
                output_tensor = calc_green_norm(vector_dist_norm, output_tensor, param, dist);
                
                
                for(int k = 0; k < 3; k = k + 1)
                {
                    for(int l = 0; l < 3; l = l + 1)
                    {
                        _im_parts(3*i+k,3*j+l) = imag(output_tensor(k,l));
                        _im_parts(3*j+k,3*i+l) = imag(output_tensor(k,l));
                    }
                }
                output_tensor = -(pow(param.wavenumber,2)/(param.eps_0*param.eps_m))*output_tensor*part[j].polarizability_tensor;
                
                for(int k = 0; k < 3; k = k + 1)
                {
                    for(int l = 0; l < 3; l = l + 1)
                    {
                        _A(3*i+k,3*j+l) = output_tensor(k,l);
                    }
                }
                
            }
        }
    }
    
    temporary.eye(3*param.N_part,3*param.N_part);
    
    _A=_A + temporary;
    
    this-> _A = A;
    this-> _im_parts = im_parts;
}
