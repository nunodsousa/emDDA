//
//  parameters.h
//
//
//  Created by Nuno de Sousa on 04/01/14.
//
//

#ifndef _parameters_h
#define _parameters_h

typedef std::complex<double> dcmplx;
using namespace arma;

class Parameters{
    
private:
    double _wavenumber;
    double _xsource_pos;
    double _ysource_pos;
    double _zsource_pos;
    double _eps_0;
    double _eps_m;
    double _E0;
    int _N_part;
    vec _u;
    
public:
    Parameters(): wavenumber(_wavenumber), xsource_pos(_xsource_pos), ysource_pos(_ysource_pos), zsource_pos(_zsource_pos), eps_0(_eps_0), eps_m(_eps_m),N_part(_N_part), E0(_E0), u(_u)
    {
        _eps_0 = 1;
        _eps_m = 1;
        _wavenumber = 0;
        _N_part = 0;
        _xsource_pos = 0;
        _ysource_pos = 0;
        _zsource_pos = 0;
        _E0 = 1;
        _u.zeros(3);
        
    }
    
    void set_plw_direction(double ux, double uy, double uz)
    {
        vec utemp(3);
        utemp(0) = ux;
        utemp(1) = uy;
        utemp(2) = uz;
        
        this-> _u = utemp;
    }
    
    void set_E0(double E0)
    {
        this-> _E0 = E0;
    }
    
    void set_source_pos(double xsource_pos, double ysource_pos, double zsource_pos)
    {
        this->_xsource_pos=xsource_pos;
        this->_ysource_pos=ysource_pos;
        this->_zsource_pos=zsource_pos;
        
        cout << "xsource inside class = " << xsource_pos << endl;
    }

    void set_parameters_wavenumber(double wavenumber)
    {
        this->_wavenumber=wavenumber;
    }
    
    void set_parameters_eps_0(double eps_0)
    {
        this->_eps_0 = eps_0;
    }
    
    void set_parameters_eps_m(double eps_m)
    {
        
        this->_eps_m = eps_m;
    }
    
    void set_N_part(int N_part)
    {
        
        this->_N_part = N_part;
    }
    
    const double &wavenumber;
    const double &eps_0;
    const double &eps_m;
    const double &E0;
    const int &N_part;
    const double &xsource_pos;
    const double &ysource_pos;
    const double &zsource_pos;
    const vec &u;
};

#endif
