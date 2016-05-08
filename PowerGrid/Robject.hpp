/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [Robject.hpp]

    Synopsis    [Base class of implementing quadratic and non-quadratic
                    regularization.]

    Description []

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/19/2016]

 *****************************************************************************/

#ifndef PowerGrid_Robject_hpp
#define PowerGrid_Robject_hpp

using namespace arma;

template <typename T1>
class Robject
{
typedef complex<T1> CxT1;
public:
    // Robject();
    Robject() {};
    //Class members
    uword Nx;
    uword Ny;
    uword Nz;
    T1 DeltaX;
    T1 DeltaY;
    T1 DeltaZ;
    T1 Beta;


    // It was declared as type Mat<uword> and the 3D type was a cube. We need to vectorize it before it is passed to QuadPenalty.
    //Custom Class Constructor
    Robject(uword nx, uword ny, uword nz, T1 beta)
    {
        //Set Class Memebers
        this->Nx = nx;
        this->Ny = ny;
        this->Nz = nz;
        this->Beta = beta;

    }

    //Class Methods - Declared virtual so they can be implemented in the base classes. Also they are virtual so that if you try to call Robject, things crash rather than give un results.
    virtual Col<CxT1> wpot(const Col<CxT1>& d) const    {
        return ones < Col < CxT1 >> (d.n_rows);
    }

    virtual Col <CxT1> dpot(const Col <CxT1>& d) const
    {
	    return d;
    }

    virtual Col <CxT1> pot(const Col <CxT1>& d) const
    {
	    Col <T1> temp = abs(d)%abs(d)/2.0;
	    return conv_to<Col<CxT1 >> ::from(temp);
    }

    Col<CxT1>Cd(const Col<CxT1>& d, uword dim) const
    {
        /*
        Col<T1> out(Nx*Ny*Nz);
        uword offset;
        for(uword  ll = 0; ll < 2; ll++) {
            for(uword jj =0; jj < 2; jj++) {
                for (uword kk = 0; kk < 2; kk++) {
                    offset = ll + jj*Ny + kk*Nx*Ny;
                    for(uword ii = offset; ii < Ny*Nx*Nz; ii++) {
                        out(ii+((ll+2*jj+4*kk)*(Nx*Ny*Nz))) = d(ii) - d(ii - offset);
                    }
                }
            }
        }
        */
        Col<CxT1> out(Nx*Ny*Nz);
        uword ll,jj,kk;
        switch (dim)
        {
            case (uword)0:
                ll = 1;
                jj = 0;
                kk = 0;
                break;
            case (uword)1:
                ll = 0;
                jj = 1;
                kk = 0;
                break;
            case (uword)2:
                ll = 0;
                jj = 0;
                kk = 1;
                break;
            default:
                cout << "Warning regularization along dimension greater than 3! Undefined case!" << endl;
        }
        uword offset = ll + jj*Ny + kk*Nx*Ny;
        for(uword ii = offset; ii < Ny*Nx*Nz; ii++) {
            out(ii) = d(ii) - d(ii - offset);
        }

        return out;
    }

    Col<CxT1>Ctd(const Col<CxT1>& d, uword dim) const
    {
        /*
        Col<T1> out(Nx*Ny*Nz);
        uword offset;
        for(uword  ll = 0; ll < 2; ll++) {
            for( uword jj =0; jj < 2; jj++) {
                for (uword kk = 0; kk < 2; kk++) {
                    offset = ll + jj*Ny + kk*Nx*Ny;
                    for(uword ii = offset; ii < Ny*Nx*Nz; ii++) {
                        out(ii) += d(ii+((ll+2*jj+4*kk)*(Nx*Ny*Nz)));
                        out(ii-offset) -= d(ii+((ll+2*jj+4*kk)*(Nx*Ny*Nz)));
                    }
                }
            }
        }
        */
        Col<CxT1> out(Nx*Ny*Nz);

	    uword ll, jj, kk;
	    switch (dim) {
	    case (uword) 0:
		    ll = 1;
		    jj = 0;
		    kk = 0;
		    break;
	    case (uword) 1:
		    ll = 0;
		    jj = 1;
		    kk = 0;
		    break;
	    case (uword) 2:
		    ll = 0;
		    jj = 0;
		    kk = 1;
		    break;
	    default:
		    cout << "Warning regularization along dimension greater than 3! Undefined case!" << endl;
	    }

	    uword offset = ll+jj*Ny+kk*Nx*Ny;
	    for (uword ii = offset; ii<Ny*Nx*Nz; ii++) {
		    if (ii==offset-1) {
			    out(ii) = -d(ii+1);
		    }
		    else if (ii==Ny*Nx*Nz-1) {
			    out(ii-offset) = d(ii-offset);
		    }
		    else {
			    out(ii-offset) = d(ii-offset)-d(ii);

		    }


	    }

        return out;
    }
public:
    T1 Penalty(const Col <CxT1>& x) const
    {
        Col <CxT1> d = zeros< Col < T1 >>(x.n_rows);
	    T1 penal = 0;
        uword nd = 0;
        if(this->Nz == 1) {
            nd = 2;
        } else {
            nd = 3;
	        //cout << "Setting dimension to 3 in reg." << endl;
        }

        for(uword ii = 0; ii < nd; ii++) {
            d = this->Cd(x,ii);
            d = this->pot(d);

            penal = penal + abs(sum(d));
        }

        return this->Beta*penal;
    }

    Col<CxT1> Gradient(const Col<CxT1>& x) const
    {
        Col <CxT1> g = zeros< Col < CxT1 >>(x.n_rows);
        Col <CxT1> d = zeros< Col < CxT1 >>(x.n_rows);
        uword nd = 0;
        if(this->Nz == 1) {
            nd = 2;
        } else {
            nd = 3;
	        //cout << "Setting dimension to 3 in reg." << endl;

        }

        for(uword ii = 0; ii < nd; ii++) {
            d = this->Cd(x,ii);
	        d = this->dpot(d);
            d = this->Ctd(d,ii);
            g = g + d;
        }

        return this->Beta*g;
    }

    CxT1 Denom(const Col<CxT1>& ddir,const Col<CxT1>& x) const
    {
	    Col <CxT1> Cdir = zeros<Col<CxT1 >> (ddir.n_rows);
        Col <CxT1> Cx = zeros< Col < CxT1 >>(ddir.n_rows);
        CxT1 penal = 0;
        CxT1 temp;
        CxT1 cxBeta = (this->Beta, 0);
        uword nd = 0;
        if(this->Nz == 1) {
            nd = 2;
        } else {
            nd = 3;
        }

        for(uword ii = 0; ii < nd; ii++) {
	        Cdir = this->Cd(ddir, ii);
            Cx = this->wpot(this->Cd(x,ii));
	        Cx = Cx%Cdir;
	        temp = as_scalar(Cdir.t()*Cx);
            penal += temp;
        }
	    //cout << "Beta = " << Beta << endl;
        return penal * cxBeta;
    }

};

#endif //PowerGrid_Robject_hpp
