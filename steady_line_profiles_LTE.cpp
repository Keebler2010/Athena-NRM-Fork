//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file orszag_tang.cpp
//! \brief Problem generator for Orszag-Tang vortex problem.
//!
//! REFERENCE: For example, see: G. Toth,  "The div(B)=0 constraint in shock capturing
//!   MHD codes", JCP, 161, 605 (2000)
//========================================================================================


// C++ headers
#include <cmath>      // sqrt()
#include <iostream>   // endl
#include <algorithm>
#include <fstream>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <chrono>     // timing

// GSL headers
#include <gsl/gsl_histogram.h>
// #include <gsl/gsl_histogram2d.h>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

#define FROZEN_HYDRO

int Nbins,Nnu_per_bin;
Real sin_i,cos_i,c_,kb_,mbar;
Real A_abun,delta_nu0,phi_norm,a12_,h_over_c2,h_over_kb;
Real gm1,xi0,v0,nu0,y_th;
Real n_gas,T_gas,r_min,r_max;
Real p_units,length_unit;
Real y_min, y_max, bin_size, dnu;
Real theta_min,theta_max;
Real Inc;
bool use_nuD, SR, argb;

gsl_histogram* L0bins = gsl_histogram_alloc(1);
using namespace std;
using namespace std::chrono;

int write_NbinsxNnu_dump(string filename, Mesh*, int n=1);
int write_Nbins_dump(string filename, Mesh*, int n=1);
int write_TotV_dump(string filename, Mesh*, int n=1);
int write_bins(string filename);
int write_Pnu_dump(string fname, MeshBlock *pmb, int n, int ie, int je, int ke, int is, int js, int ks);
int write_nu_cmv_dump(string fname, MeshBlock *pmb, int n, int ie, int je, int ke, int is, int js, int ks, int Nbins, int Nnu_per_bin);
string default_filename(string label);

template <class T>
inline std::string to_string (const T& t)
{
   std::stringstream ss;
   ss << t;
   return ss.str();
}

void calc_velocity(Real r, Real &v_r, Real &v_t, Real &v_p) {
  // cout << "v0 = " << v0
  //      << ", r = " << r 
  //      << ", r_min = " << r_min 
  //      << endl;

  // v_r = v0*(r/r_min); // Welsh-Horne test
  v_r = v0; // Shu Vol 1 test (page 97)
  v_t = 0.;
  v_p = 0.;
}

Real calc_density(Real r, Real r_max) {
  if (r < r_max)
    return mbar*n_gas;
  return 0;
}


Real calc_temperature(Real r) {
  return T_gas;
}

Real calc_opacity(Real xi, Real T) {
  return 1.;
}



Real calc_Doppler_Shift(Real v, Real v_los, bool SR)
{
  Real LF = 1/sqrt(1 -(v*v));
  Real Dop_Shift;
    
  if (SR == 1)
    {Dop_Shift =  1.0/(LF*(1.-v_los));} // lab frame n, use equation 89.5 from Mihalas Mihalas

  else
    {Dop_Shift = (1.0 + v_los);}    

  return Dop_Shift;
}




Real calc_cmv_nu(Real nuD, Real nu)
{
  // Real dnu = nu - nuD;
  // Real nu_cmv = nu0 + dnu;
  Real nu_cmv = (nu/nuD)*nu0;

  return nu_cmv;
}

Real calc_arg(Real nu_cmv, bool argb, Real nu, Real nuD)
{
  Real arg;
  if (argb)
  {
    arg = (nu - nuD)/delta_nu0;
  }
  else 
  {
    arg = (nu_cmv - nu0)/delta_nu0;
  }


  return arg;
}

Real calc_alpha(Real n, Real T, Real nuD, Real nu, bool SR, bool argb)
{

  Real nu_cmv = calc_cmv_nu(nuD, nu);


  Real eta_ion = 1.0; // this needs to be passed in later
  Real n_1 = A_abun*eta_ion*n; // level population of lower level

  // evaluate correction factor for stimulated emission
  Real fac = 1. - exp(-h_over_kb*(nu_cmv)/T);
  // Real fac = 1. - exp(-h_over_kb*nu/T);
  // Real fac = 1. - exp(-h_over_kb*nuD/T);
  // Real fac = 1.;

  // profile function
  Real arg = calc_arg(nu_cmv, argb, nu, nuD);
  Real phi_nu = phi_norm*exp(-arg*arg);

  

    if (SR)
  {
    return (nu_cmv/nu)*n_1*a12_*fac*phi_nu;  // added relativistic correction to extinction
  }
  else 
  {
    return n_1*a12_*fac*phi_nu;
  }
  
}

Real calc_Bnu(Real T, Real nuD, Real nu, bool SR)
{
  Real nu_cmv = calc_cmv_nu(nuD, nu);

  if (SR)
  {
    return (nu/nu_cmv)*(nu/nu_cmv)*(nu/nu_cmv)*2.*h_over_c2*pow(nu_cmv,3.)/(exp(h_over_kb*nu_cmv/T) - 1.0); // added relativistic correction to source function
  }
  else 
  {
    return 2.*h_over_c2*pow(nu_cmv,3.)/(exp(h_over_kb*nu_cmv/T) - 1.0); 
  }
  
}

Real calc_jnu(Real n, Real T, Real nuD, Real nu, bool SR)
{

  Real nu_cmv = calc_cmv_nu(nuD, nu);
  Real alpha_nu = calc_alpha(n,T,nuD,nu,SR,argb);
  

  return alpha_nu*calc_Bnu(T,nuD,nu,SR); 
}

class Projections {
  public:
    Real sin_theta,cos_theta;
    Real sin_phi,cos_phi;
    Real n_r,n_t,n_p;
    Real n_x,n_y,n_z;
    Real rhat_dot_xhat,rhat_dot_yhat,rhat_dot_zhat,
         that_dot_xhat,that_dot_yhat,that_dot_zhat,
         phat_dot_xhat,phat_dot_yhat,phat_dot_zhat;
    Projections (Real,Real,Real,Real);  // constructor
    
    // member functions
    void eval_dot_products(Real, Real, Real, Real);
};

Projections::Projections(Real a, Real b, Real c, Real d)
{

  // inputs are sin(theta), cos(theta), sin(phi), cos(phi)
  sin_theta = a; 
  cos_theta = b;
  sin_phi   = c;
  cos_phi   = d;

  // calculate components of n_hat (direction toward observer)
  n_r = sin_theta*cos_phi*sin_i + cos_theta*cos_i;
  n_t = cos_theta*cos_phi*sin_i - sin_theta*cos_i;
  n_p = -sin_phi*sin_i;
  // n_x = n_r*sin_theta*cos_phi + n_t*cos_theta*cos_phi - n_p*sin_phi;
  // n_y = n_r*sin_theta*sin_phi + n_t*cos_theta*sin_phi + n_p*cos_phi;
  // n_z = n_r*cos_theta         - n_t*sin_theta;
}

void Projections::eval_dot_products(Real r, Real x, Real y, Real z)
{

  // figure out location of the point in spherical coordinates
  Real cos_theta_ = z/r;
  Real sin_theta_ = sqrt(1.-SQR(cos_theta_));  
  Real phi = atan2(y,x);   // note: atan2 is atan accounting for sign(x),sign(y)
  Real cos_phi_ = cos(phi);
  Real sin_phi_ = sin(phi);

  // set cartesian projections of spherical unit vectors
  rhat_dot_xhat = sin_theta_*cos_phi_;
  // rhat_dot_yhat = sin_theta_*sin_phi_;
  rhat_dot_zhat = cos_theta_;
  that_dot_xhat = cos_theta_*cos_phi_;
  // that_dot_yhat = cos_theta_*sin_phi_;
  that_dot_zhat = -sin_theta_;
  phat_dot_xhat = -sin_phi_;
  // phat_dot_yhat = cos_phi_;
  phat_dot_zhat = 0.;

}


void calc_LP_using_binning   (MeshBlock*, int, int, int, AthenaArray<Real>, Projections, int ctr);

void Mesh::InitUserMeshData(ParameterInput *pin) {
  // physical Constants
  Real mp_ = 1.67262178e-24; // proton mass
  Real h_ = 6.626196e-27; // Planck's constant
  c_ = 2.99792458e10; // speed of light
  kb_ = 1.380658e-16; // Boltzmann's constant
  mbar = pin->GetOrAddReal("problem","mu",0.6)*mp_;
  Real a12_over_f = 0.02654; // a_12 divided by oscillator strength (equals pi*e^2/m_e/c) 
  h_over_c2 = h_/c_/c_; 
  h_over_kb = h_/kb_;
  Real t_1day = 24.*3600.;  // seconds in a day
  length_unit = pin->GetOrAddReal("problem","length_unit",c_*t_1day); // one light-day
  //length_unit = 1; // JAK added for testing volumes and areas

  // Model parameters
  n_gas = pin->GetOrAddReal("problem","n_gas",1e9);  // gas number density
  T_gas = pin->GetOrAddReal("problem","T_gas",1e5);  // gas temperature
  r_min = pin->GetReal("mesh","x1min");
  r_max = pin->GetReal("mesh","x1max");
  Real v0_kms = pin->GetOrAddReal("problem","v0_kms",1e3);
  

  v0 = v0_kms*1e5/c_; // v0 in units of c
  
  // Real y_min = pin->GetorAddReal("problem","y_min", -v0_kms/c_)
  // Real y_max = pin->GetorAddReal("problem","y_max",  v0_kms/c_)
  

  p_units = 1./c_/c_;
  xi0 = pin->GetOrAddReal("problem","xi0",1.);
  Real i_los = pin->GetReal("problem","i_los")*PI/180.;  // input file needs to specify observer's LOS in degrees

  // Functionality parameters
  theta_min = pin->GetOrAddReal("problem","theta_min",0.);    // ignore zones with theta < theta_min
  theta_max = pin->GetOrAddReal("problem","theta_max",180.);  // ignore zones with theta > theta_max
  use_nuD = pin->GetOrAddBoolean("problem","use_nuD",false);  // set every freq. in the Gaussian to nuD
  SR = pin->GetOrAddBoolean("problem","SR",false);  // Use inertial effects of special relativity if true.  Default use non relativistic 
  argb = pin->GetOrAddBoolean("problem","argb",false);  // Use inertial effects of special relativity if true.  Default use non relativistic 


  // parameters set by the choice of ion/transition (default values are for Hbeta)
  Real atomic_number = pin->GetOrAddReal("problem","atomic_number",1.);  // number of protons + neutrons
  Real m_ion = mp_*atomic_number;
  Real v_th = sqrt(2.*kb_*T_gas/m_ion);
  y_th = v_th/c_;
  nu0 = pin->GetOrAddReal("problem","nu0",6.1668560e14);  // rest freq of line
  Real f_12 = pin->GetOrAddReal("problem","f_12",0.1193); // oscillator strength 
  A_abun = pin->GetOrAddReal("problem","A_abun",0.98); // atomic abundance
  a12_ = a12_over_f*f_12;
  delta_nu0 = nu0*v_th/c_;
  phi_norm = 1./(sqrt(PI)*delta_nu0); // assumes a Gaussian line profile

  // freq grid parameters
  Nbins = pin->GetOrAddInteger("problem","Nbins",500);  // number of bins
  Nnu_per_bin = pin->GetOrAddInteger("problem","Nnu_per_bin",10);  // number of frequencies per bin
  Inc = pin->GetOrAddReal("problem","Inc",2);  // divides the dr unit by some amount to get the dl unit


  y_min = -v0;
  y_max = v0;

  // std::cout << "y_max = " << y_max << std::endl;
  // std::cout << "v0 = " << v0 << std::endl;
  // std::cout << "4*y_th = " << 4*y_th << std::endl;
  // std::cout << "True/False = " << (v0 < 4*y_th) << std::endl;


  if (v0 < 3*y_th){   //  For cases approaching static
    y_min = -3*y_th;
    y_max = 3*y_th;}

  
  // std::cout << "y_max = " << y_max << std::endl;
  // std::cout << "y_th = " << y_th << std::endl;
  
  // initialize freq grid
  L0bins = gsl_histogram_alloc(Nbins); // reallocate histogram
  
  if (SR)
    {
    gsl_histogram_set_ranges_uniform (L0bins, (sqrt((1+y_min)/(1-y_min))-1), (sqrt((1+y_max)/(1-y_max))-1));
    } 
  else
    {
      gsl_histogram_set_ranges_uniform (L0bins, y_min, y_max);
    }

  bin_size = fabs(L0bins->range[1] - L0bins->range[0]);  // Gets absolute value of the size of the first bin
  dnu = nu0*bin_size/(Real)Nnu_per_bin;
  for (int ii=0; ii<Nbins; ii++) 
    cout << "L0bins[" << ii << "] = " << L0bins->range[ii] << endl;
  std::string file_string = default_filename("bins");
  write_bins(file_string);
  cout << "\n[InitUserMeshData]: Nbins = " << Nbins << ", Nnu_per_bin = " << Nnu_per_bin << endl;
  cout << "y_th = " << y_th << endl;
  cout << "bin_width/y_th = " << bin_size/y_th << endl;
  cout << "delta_nu0/dnu = " << delta_nu0/dnu << endl;
  cout << "T_line (h*nu0/kb) = " << h_over_kb*nu0 << endl;
  cout << "Correction for stimulated emission = " << 1. - exp(-h_over_kb*nu0/T_gas) << endl;
  cout << "a12 = " << a12_ << endl;
  cout << "delta_nu0 = " << delta_nu0 << endl;
  cout << "Bnu0 = " << calc_Bnu(T_gas,nu0,nu0,SR) << endl;
  cout << "alphanu = " << calc_alpha(n_gas,T_gas,nu0,nu0,SR,argb) << endl;
  cout << "phi_nu0 = " << phi_norm << endl; 

  // Shu test
  Real prefactor_ShuTest = A_abun*n_gas*a12_*(h_over_c2*nu0*nu0)*exp(-h_over_kb*nu0/T_gas)/v0;
  cout << "\nShu test prefactor = " << prefactor_ShuTest << endl;

  // set model-independent globals 
  Real gamma_ = pin->GetReal("hydro","gamma"); //peos->GetGamma() - 1.0;
  gm1 = gamma_ - 1.;
  sin_i = sin(i_los);
  cos_i = cos(i_los);
  
  AllocateIntUserMeshDataField(1); 
  iuser_mesh_data[0].NewAthenaArray(Nbins);

  AllocateRealUserMeshDataField(10); 
  ruser_mesh_data[0].NewAthenaArray(Nbins,Nnu_per_bin);
  ruser_mesh_data[1].NewAthenaArray(Nbins,Nnu_per_bin);  
  ruser_mesh_data[2].NewAthenaArray(Nbins,Nnu_per_bin);
  ruser_mesh_data[3].NewAthenaArray(Nbins,Nnu_per_bin);
  ruser_mesh_data[4].NewAthenaArray(Nbins);
  ruser_mesh_data[5].NewAthenaArray(Nbins);
  ruser_mesh_data[6].NewAthenaArray(Nbins,Nnu_per_bin);
  ruser_mesh_data[7].NewAthenaArray(1);
  ruser_mesh_data[8].NewAthenaArray(1);
  ruser_mesh_data[9].NewAthenaArray(Nbins,Nnu_per_bin);
  //ruser_mesh_data[10].NewAthenaArray(1);
  cout << "\n[InitUserMeshData]: Done." << endl;

  return;

}


void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{

  // user arrays
  AllocateRealUserMeshBlockDataField(7);
  AllocateUserOutputVariables(7);  

  ruser_meshblock_data[0].NewAthenaArray((ke-ks)+2*NGHOST,(je-js)+2*NGHOST,(ie-is)+2*NGHOST); 
  ruser_meshblock_data[1].NewAthenaArray((ke-ks)+2*NGHOST,(je-js)+2*NGHOST,(ie-is)+2*NGHOST);
  ruser_meshblock_data[2].NewAthenaArray((ke-ks)+2*NGHOST,(je-js)+2*NGHOST,(ie-is)+2*NGHOST); 
  ruser_meshblock_data[3].NewAthenaArray((ke-ks)+2*NGHOST,(je-js)+2*NGHOST,(ie-is)+2*NGHOST);
  ruser_meshblock_data[4].NewAthenaArray((ke-ks)+2*NGHOST,(je-js)+2*NGHOST,(ie-is)+2*NGHOST);
  ruser_meshblock_data[5].NewAthenaArray((ke-ks)+2*NGHOST,(je-js)+2*NGHOST,(ie-is)+2*NGHOST);
  ruser_meshblock_data[6].NewAthenaArray((ke-ks)+2*NGHOST,(je-js)+2*NGHOST,(ie-is)+2*NGHOST);
  // ruser_meshblock_data[7].NewAthenaArray((ke-ks)+2*NGHOST,(je-js)+2*NGHOST,(ie-is)+2*NGHOST, Nbins, Nnu_per_bin);
  // ruser_meshblock_data[8].NewAthenaArray((ke-ks)+2*NGHOST,(je-js)+2*NGHOST,(ie-is)+2*NGHOST, Nbins, Nnu_per_bin);
  // ruser_meshblock_data[9].NewAthenaArray((ke-ks)+2*NGHOST,(je-js)+2*NGHOST,(ie-is)+2*NGHOST, Nbins, Nnu_per_bin);
  // ruser_meshblock_data[10].NewAthenaArray((ke-ks)+2*NGHOST,(je-js)+2*NGHOST,(ie-is)+2*NGHOST, Nbins, Nnu_per_bin);
  // ruser_meshblock_data[11].NewAthenaArray((ke-ks)+2*NGHOST,(je-js)+2*NGHOST,(ie-is)+2*NGHOST, Nbins, Nnu_per_bin);
  // ruser_meshblock_data[12].NewAthenaArray((ke-ks)+2*NGHOST,(je-js)+2*NGHOST,(ie-is)+2*NGHOST, Nbins, Nnu_per_bin);
  // ruser_meshblock_data[13].NewAthenaArray((ke-ks)+2*NGHOST,(je-js)+2*NGHOST,(ie-is)+2*NGHOST, Nbins, Nnu_per_bin);
  // ruser_meshblock_data[14].NewAthenaArray((ke-ks)+2*NGHOST,(je-js)+2*NGHOST,(ie-is)+2*NGHOST, Nbins, Nnu_per_bin);
  // ruser_meshblock_data[15].NewAthenaArray((ke-ks)+2*NGHOST,(je-js)+2*NGHOST,(ie-is)+2*NGHOST, Nbins, Nnu_per_bin);

  return;

}


//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for computing steady line profiles.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {


  auto start = high_resolution_clock::now();
  cout << "\n[ProblemGenerator]: Timing the binning method... " << endl;

  // temporarily store the frequency and LOS velocity grids for the binning method
  // (these arrays get overwritten below)
  const size_t nx = L0bins->n;
  for (int i=0; i<nx; i++) {
    Real nu_min = nu0*(1.+L0bins->range[i]);
    Real nu1 = nu_min + 0.5*dnu;

    for (int idx=0; idx < Nnu_per_bin; idx++) {
      Real nu = nu1 + idx*dnu;
      Real y = nu/nu0 - 1.;
      pmy_mesh->ruser_mesh_data[0](i,idx) = nu;
      pmy_mesh->ruser_mesh_data[1](i,idx) = y;
    }
  }

  // write the frequency grid
  std::string dump_string = default_filename("freq_grid");
  if (write_NbinsxNnu_dump(dump_string, this->pmy_mesh, 0) == 0)
    cout << "[ProblemGenerator]: Done writing " << dump_string << "." << endl;

  // write the LOS velocity grid
  dump_string = default_filename("vlos_grid");
  if (write_NbinsxNnu_dump(dump_string, this->pmy_mesh, 1) == 0)
    cout << "[ProblemGenerator]: Done writing " << dump_string << "." << endl;

  // accessor for zone volumes
  AthenaArray<Real> vol;
  vol.NewAthenaArray(ncells1); 
  Real TVol = 0;
  Real dV = 0;
  Real ctr = 0;
  // Loop over all zones
  for (int k=ks; k<=ke; k++) {
    Real sin_phi = sin(pcoord->x3v(k));
    Real cos_phi = cos(pcoord->x3v(k));
    cout << "k = " << k
         << ", phi = " << pcoord->x3v(k)*180./PI << " deg" << endl;

    for (int j=js; j<=je; ++j)
    {
      pcoord->CellVolume(k, j, is, ie, vol);
      Real sin_theta = sin(pcoord->x2v(j));
      Real cos_theta = cos(pcoord->x2v(j));
      Real theta_deg = pcoord->x2v(j)*180./PI;
      // if (k==ks)
      //   cout << "  j = " << j
      //        << "  theta = " << theta_deg << endl;

      Projections proj(sin_theta,cos_theta,sin_phi,cos_phi);

      for (int i=is; i<=ie; i++) {
        Real r = pcoord->x1v(i);
        dV = vol(i)*pow(length_unit,3);
        TVol += dV;
        ctr += 1;

        Real rho,T,v_r,v_t,v_p;
        rho = calc_density(r, r_max);
        T = calc_temperature(r);
        calc_velocity(r,v_r,v_t,v_p);

        // pressure
        Real p = (rho/mbar)*kb_*T/p_units;

        // cout << "  i = " << i
        //      << "  r = " << r
        //      << "  rho = " << rho
        //      << ", v_r = " << v_r
        //      << ", T = " << T
        //      << ", p = " << p
        //      << endl;

        phydro->u(IDN,k,j,i) = rho;
        phydro->u(IM1,k,j,i) = rho*v_r;
        phydro->u(IM2,k,j,i) = rho*v_t;
        phydro->u(IM3,k,j,i) = rho*v_p;

        if (NON_BAROTROPIC_EOS) 
        phydro->u(IEN,k,j,i) =
              p/gm1 + 0.5*(SQR(phydro->u(IM1,k,j,i)) + 
                           SQR(phydro->u(IM2,k,j,i)) + 
                           SQR(phydro->u(IM3,k,j,i)))/phydro->u(IDN,k,j,i);


        // calculate the steady line profile using binning method 
        // over the specified theta range (default is 0 to 180 deg.)
        if (theta_deg > theta_min && theta_deg < theta_max)
          calc_LP_using_binning(this,k,j,i,vol,proj,ctr);
          //this->pmy_mesh->ruser_mesh_data[7](1) += dV;
          //this->pmy_mesh->ruser_mesh_data[8](1) += 1;
        if (i==ie && j==je && k == ke){
          cout << "  TVol = " << TVol << endl;
          cout << " Total Cells Hit = " << i*j*k << endl;}

      }
    }
  }

  // Average the luminosities within each bin
  Real delta_nu = nu0*bin_size; // freq range in each bin 

  for (int i=0; i<nx; i++) {
    Real Lnu = 0.; 
    Real Lnu_pnu1 = 0.;

    for (int idx=0; idx < Nnu_per_bin; idx++) {
      Lnu      += pmy_mesh->ruser_mesh_data[2](i,idx)*dnu;
      Lnu_pnu1 += pmy_mesh->ruser_mesh_data[3](i,idx)*dnu;
    }
    pmy_mesh->ruser_mesh_data[4](i) = Lnu/delta_nu;
    pmy_mesh->ruser_mesh_data[5](i) = Lnu_pnu1/delta_nu;
    // if (use_nuD) { // in "Sobolev mode", divide by bin count
    //   Real count = max(1.,(Real)pmy_mesh->iuser_mesh_data[0](i));
    //   pmy_mesh->ruser_mesh_data[2](i) /= count;
    //   pmy_mesh->ruser_mesh_data[3](i) /= count;
    //   pmy_mesh->ruser_mesh_data[4](i) /= count;
    //   pmy_mesh->ruser_mesh_data[5](i) /= count;
    // }
  }


  // write the steady line profile computed before averaging over frequency
  dump_string = default_filename("Lsteady");
  if (write_NbinsxNnu_dump(dump_string, this->pmy_mesh, 2) == 0){
    cout << "[ProblemGenerator]: Dump " << dump_string << " is the non-averaged line profile computed at t=0." << endl;}

  // write the (p_nu=1) steady line profile computed before averaging over frequency
  dump_string = default_filename("Lsteady-pnu1");
  if (write_NbinsxNnu_dump(dump_string, this->pmy_mesh, 3) == 0){
    cout << "[ProblemGenerator]: Dump " << dump_string << " is the non-averaged line profile (with p_nu=1) computed at t=0." << endl;}

  // // write the (p_nu=1) steady line profile computed before averaging over frequency
  // dump_string = default_filename("Bin Volume");
  // if (write_NbinsxNnu_dump(dump_string, this->pmy_mesh, 6) == 0){
  //   cout << "[ProblemGenerator]: Dump " << dump_string << " is the non-averaged volume computed at t=0." << endl;}


  // // write the (p_nu=1) steady line profile computed before averaging over frequency
  // dump_string = default_filename("Bin Count");
  // if (write_NbinsxNnu_dump(dump_string, this->pmy_mesh, 9) == 0){
  //   cout << "[ProblemGenerator]: Dump " << dump_string << " is the number of volume entries in to NbinxNu computed at t=0." << endl;}
    
  

  // write the steady line profile averaged over frequency
  dump_string = default_filename("Lavg");
  if (write_Nbins_dump(dump_string, this->pmy_mesh, 4) == 0){
    cout << "[ProblemGenerator]: Dump " << dump_string << " is the averaged line profile computed at t=0." << endl;}

  // write the (p_nu=1) steady line profile averaged over frequency
  dump_string = default_filename("Lavg-pnu1");
  if (write_Nbins_dump(dump_string, this->pmy_mesh, 5) == 0){
    cout << "[ProblemGenerator]: Dump " << dump_string << " is the averaged line profile (with p_nu=1) computed at t=0." << endl;}
 
  // // write the (p_nu=1) steady line profile averaged over frequency
  // dump_string = default_filename("Total Volume");
  // if (write_TotV_dump(dump_string, this->pmy_mesh, 7) == 0){
  //   cout << "[ProblemGenerator]: Dump " << dump_string << " is the total computed at t=0." << endl;}
  
  // // write the (p_nu=1) steady line profile averaged over frequency
  // dump_string = default_filename("Cell Count");
  // if (write_TotV_dump(dump_string, this->pmy_mesh, 8) == 0){
  //   cout << "[ProblemGenerator]: Dump " << dump_string << " is the total computed at t=0." << endl;}

  
  // // Write pnu_local_avg
  // std::string dump_string_vlos = default_filename("vlos_ijk");
  // if (write_Pnu_dump(dump_string_vlos, this, 0, ie, je, ke, is, js, ks) == 0)
  //   std::cout << "Done writing pnu_local_avg to " << dump_string_vlos << "." << std::endl;

  // // Write pnu_nonlocal_avg
  // std::string dump_string_dv = default_filename("dV");
  // if (write_Pnu_dump(dump_string_dv, this, 1, ie, je, ke, is, js, ks) == 0)
  //   std::cout << "Done writing pnu_nonlocal_avg to " << dump_string_dv << "." << std::endl;

  // // Write doppler shift non relativistic
  // std::string dump_string_nuD = default_filename("nuD");
  // if (write_Pnu_dump(dump_string_nuD, this, 2, ie, je, ke, is, js, ks) == 0)
  //   std::cout << "Done writing doppler shift non relativistic to " << dump_string_nuD << "." << std::endl;

  // // Write doppler shift relativistic
  // std::string dump_string_nuD1 = default_filename("nuD1");
  // if (write_Pnu_dump(dump_string_nuD1, this, 3, ie, je, ke, is, js, ks) == 0)
  //   std::cout << "Done writing doppler shift relativistic to " << dump_string_nuD1 << "." << std::endl;


  //   // Write comoving frequency nonrelativistic
  //   std::string dump_string_nucmv = default_filename("nucmv_nonrel");
  //   if (write_nu_cmv_dump(dump_string_nucmv, this, 7, ie, je, ke, is, js, ks, Nbins, Nnu_per_bin) == 0)
  //     std::cout << "Done writing nu_cmv to " << dump_string_nucmv << "." << std::endl;

  //   // Write comoving frequency relativistic
  //   std::string dump_string_nucmv1 = default_filename("nucmv1_rel");
  //   if (write_nu_cmv_dump(dump_string_nucmv1, this, 8, ie, je, ke, is, js, ks, Nbins, Nnu_per_bin) == 0)
  //     std::cout << "Done writing nu_cmv1 to " << dump_string_nucmv1 << "." << std::endl;
    
  //   // Write comoving frequency nonrelativistic
  //   std::string dump_string_nuD_e = default_filename("nuD_nonrel");
  //   if (write_nu_cmv_dump(dump_string_nuD_e, this, 9, ie, je, ke, is, js, ks, Nbins, Nnu_per_bin) == 0)
  //     std::cout << "Done writing nu_D to " << dump_string_nuD_e << "." << std::endl;

  //   // Write comoving frequency relativistic
  //   std::string dump_string_nuD1_e = default_filename("nuD1_rel");
  //   if (write_nu_cmv_dump(dump_string_nuD1_e, this, 10, ie, je, ke, is, js, ks, Nbins, Nnu_per_bin) == 0)
  //     std::cout << "Done writing nu_D1 to " << dump_string_nuD1_e << "." << std::endl;


  //     // pmb->ruser_meshblock_data[11](k,j,i,idx,inu) = nu;
  //     // pmb->ruser_meshblock_data[12](k,j,i,idx,inu) = v_los;
  //     // pmb->ruser_meshblock_data[13](k,j,i,idx,inu) = v_;
  //     // pmb->ruser_meshblock_data[14](k,j,i,idx,inu) = arg;
  //     // pmb->ruser_meshblock_data[15](k,j,i,idx,inu) = arg1;
  

  //   // Write comoving frequency relativistic
  //   std::string dump_string_nu_e = default_filename("nu");
  //   if (write_nu_cmv_dump(dump_string_nu_e, this, 11, ie, je, ke, is, js, ks, Nbins, Nnu_per_bin) == 0)
  //     std::cout << "Done writing nu_D1 to " << dump_string_nu_e << "." << std::endl;

  //   // Write comoving frequency relativistic
  //   std::string dump_string_vlos_e = default_filename("v_los");
  //   if (write_nu_cmv_dump(dump_string_vlos_e, this, 12, ie, je, ke, is, js, ks, Nbins, Nnu_per_bin) == 0)
  //     std::cout << "Done writing nu_D1 to " << dump_string_vlos_e << "." << std::endl;

  //         // Write comoving frequency relativistic
  //   std::string dump_string_v_e = default_filename("v_");
  //   if (write_nu_cmv_dump(dump_string_v_e, this, 13, ie, je, ke, is, js, ks, Nbins, Nnu_per_bin) == 0)
  //     std::cout << "Done writing nu_D1 to " << dump_string_v_e << "." << std::endl;

  //   // Write comoving frequency relativistic
  //   std::string dump_string_arg_e = default_filename("arg_nonrel");
  //   if (write_nu_cmv_dump(dump_string_arg_e, this, 14, ie, je, ke, is, js, ks, Nbins, Nnu_per_bin) == 0)
  //     std::cout << "Done writing nu_D1 to " << dump_string_arg_e << "." << std::endl;

  //   // Write comoving frequency relativistic
  //   std::string dump_string_arg1_e = default_filename("arg1_rel");
  //   if (write_nu_cmv_dump(dump_string_arg1_e, this, 15, ie, je, ke, is, js, ks, Nbins, Nnu_per_bin) == 0)
  //     std::cout << "Done writing nu_D1 to " << dump_string_arg1_e << "." << std::endl;


 // Added to output the escape probabliity Jan 9 2025 JAK

  // // Write pnu_local_avg
  // std::string dump_string_local = default_filename("pnu_local_avg");
  // if (write_Pnu_dump(dump_string_local, this, 0, ie, je, ke, is, js, ks) == 0)
  //   std::cout << "Done writing pnu_local_avg to " << dump_string_local << "." << std::endl;

  // // Write pnu_nonlocal_avg
  // std::string dump_string_nonlocal = default_filename("pnu_nonlocal_avg");
  // if (write_Pnu_dump(dump_string_nonlocal, this, 1, ie, je, ke, is, js, ks) == 0)
  //   std::cout << "Done writing pnu_nonlocal_avg to " << dump_string_nonlocal << "." << std::endl;

  // // Write cell match with pnu
  // std::string dump_string_count = default_filename("pnu_count_avg");
  // if (write_Pnu_dump(dump_string_count, this, 4, ie, je, ke, is, js, ks) == 0)
  //   std::cout << "Done writing pnu_nonlocal_avg to " << dump_string_nonlocal << "." << std::endl;
  
  // // Write cell match with pnu
  // std::string dump_string_tau = default_filename("pnu_tau_avg");
  // if (write_Pnu_dump(dump_string_tau, this, 5, ie, je, ke, is, js, ks) == 0)
  //   std::cout << "Done writing pnu_nonlocal_avg to " << dump_string_nonlocal << "." << std::endl;
    
  // // Write cell match with pnu
  // std::string dump_string_path = default_filename("pnu_path_path");
  // if (write_Pnu_dump(dump_string_path, this, 6, ie, je, ke, is, js, ks) == 0)
  //   std::cout << "Done writing pnu_nonlocal_avg to " << dump_string_nonlocal << "." << std::endl;
  
    

  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<seconds>(stop - start);
  cout << "[ProblemGenerator]: Binning method took " << duration.count() << " seconds." << endl;


#ifdef FROZEN_HYDRO
  this->pmy_mesh->UserWorkAfterLoop(pin);
  this->UserWorkBeforeOutput(pin);
  std::stringstream msg;
        msg << "### FROZEN_HYDRO MODE: exiting."   << endl;
          ATHENA_ERROR(msg);
#endif

  return;
}

void calc_LP_using_binning(MeshBlock* pmb, int k, int j, int i, AthenaArray<Real> vol, Projections proj, int ctr) 
{



  Hydro *phydro = pmb->phydro;
  Real r = pmb->pcoord->x1v(i);
  Real rho = phydro->u(IDN,k,j,i);
  Real n_cgs = rho/mbar;
  Real T_cgs = calc_temperature(r);

  //Real v_th = sqrt(2.*kb_*T_cgs/m_ion);
  //y_th = v_th/c_;

  // calculate LOS velocity
  Real v_r = phydro->u(IM1,k,j,i)/phydro->u(IDN,k,j,i);
  Real v_t = phydro->u(IM2,k,j,i)/phydro->u(IDN,k,j,i);
  Real v_p = phydro->u(IM3,k,j,i)/phydro->u(IDN,k,j,i);
  Real v_los = v_r*proj.n_r + v_t*proj.n_t + v_p*proj.n_p;
  Real v_ = sqrt(v_r*v_r + v_t*v_t + v_p*v_p);


  // Doppler shift of the line center at the point i,j,k in the flow. 
  Real nuD = nu0*calc_Doppler_Shift(v_,v_los,SR);

  // Unitless frequency to look up bin
  Real f = nuD/nu0 - 1;

  // Running Tally of Tau in the wing, if the wing tau is >10 then the entire band is thick and the escape probability -> 0 (skip remaining line of sight)

  Real tauW = 0;

  // verify that v_los is equivalently computed from v_x, v_z projections
  // Real x_ = r*proj.sin_theta*proj.cos_phi;
  // Real y_ = r*proj.sin_theta*proj.sin_phi;
  // Real z_ = r*proj.cos_theta;
  // proj.eval_dot_products(r,x_,y_,z_);
  // Real vx = v_r*proj.rhat_dot_xhat + v_t*proj.that_dot_xhat + v_p*proj.phat_dot_xhat;
  // Real vz = v_r*proj.rhat_dot_zhat + v_t*proj.that_dot_zhat + v_p*proj.phat_dot_zhat;
  // if (j==18 && i==10) // print along a circle in phi direction
  //   cout <<"phi = " << pmb->pcoord->x3v(k)*180./PI
  //        << ", v_los(spherical) = " << v_r*proj.n_r + v_t*proj.n_t + v_p*proj.n_p
  //        << ", v_los(cartesian) = " << vx*sin_i + vz*cos_i << endl;

  // identify velocity and time-delay bin this zone contributes to
  size_t i_y;
  gsl_histogram_find(L0bins, f, &i_y);

  // figure out ranges of frequency-loops
  Real vrange_min = y_min;
  Real vrange_max = y_max;
  int i_edge;
  int yidx_L = 0; // will only stay 0 for first bin 
  int yidx_R = 1; // will only stay 1 for last bin 
  if (i_y>0 && i_y<(Nbins-1)) {
    vrange_min = L0bins->range[i_y-1]; 
    vrange_max = L0bins->range[i_y+1] + bin_size;
    yidx_L = i_y-1;
    yidx_R = i_y+1;
    i_edge = 0;
  }
  else {
    if (i_y==0) {
      i_edge = 1;     // freqs span bins L0bins->range[0] and L0bins->range[1]
      vrange_max = L0bins->range[2];
    }
    else {
      i_edge = Nbins-1;  // freqs span bins L0bins->range[Nbins-2] and L0bins->range[Nbins-1]
      vrange_min = L0bins->range[Nbins-2];
    }
  }


  // initialize beta arrays to 1
  int idx;
  for (int yidx=yidx_L; yidx<=yidx_R; yidx++) {
    if (i_edge == 0)
      idx = yidx;
    else // for edge cases, yidx is 0 or 1
      idx = i_edge - yidx;
    for (int inu=0; inu < Nnu_per_bin; inu++) {
      pmb->pmy_mesh->ruser_mesh_data[0](idx,inu) = 1.;
      pmb->pmy_mesh->ruser_mesh_data[1](idx,inu) = 1.;

    }
  }

  // Binned method to compute the directional escape probability
  Real dl = pmb->pcoord->dx1v(i);
  Real x0 = r*proj.sin_theta*proj.cos_phi;
  Real y0 = r*proj.sin_theta*proj.sin_phi;
  Real z0 = r*proj.cos_theta;
  Real r_l = r;
  Real n_l = n_cgs;
  Real xi_l = xi0/(n_l*r_l*r_l);

  // use '_lm1' for quantities getting differentiated
  Real v_l,T_l,kappa_l,nuD_l;
  Real n_lm1 = n_l;
  Real T_lm1 = T_cgs;
  // Real kappa_lm1 = calc_opacity(xi_l,T_l);






  // Doppler shift at a point along the line of sight
  Real v_los_lm1 = v_los;
  Real v_lm1 = v_;
  Real nuD_lm1 = nu0*calc_Doppler_Shift(v_lm1,v_los_lm1,SR);

  int il = 0;

  while (r_l < r_max) {
    il += 1; 

    // define LOS variables
    Real l = dl*il;
    Real x = x0 + l*sin_i;
    Real y = y0;
    Real z = z0 + l*cos_i;
    
    // distance from origin at this point along the LOS
    r_l = sqrt(x*x + y*y + z*z);

    // evaluate v_x,v_y,v_z from v_r,v_t,v_p by projection
    calc_velocity(r_l,v_r,v_t,v_p);
    proj.eval_dot_products(r_l,x,y,z);
    Real v_x = v_r*proj.rhat_dot_xhat + v_t*proj.that_dot_xhat + v_p*proj.phat_dot_xhat;
    //Real v_y = v_r*proj.rhat_dot_yhat + v_t*proj.that_dot_yhat + v_p*proj.phat_dot_yhat;
    Real v_z = v_r*proj.rhat_dot_zhat + v_t*proj.that_dot_zhat + v_p*proj.phat_dot_zhat;

    // evaluate velocity along LOS
    v_l = v_x*sin_i + v_z*cos_i;

    // Doppler shifted line center freq
    Real v_1 = sqrt(v_r*v_r + v_t*v_t + v_p*v_p);
    Real nuD_l = nu0*calc_Doppler_Shift(v_1,v_l,SR);

    

    // calculate the Sobolev length (for comparison purposes)
    Real Q,l_sob;
    if (il == 1) {
      Q = (v_l - v_lm1)/dl;  // velocity gradient along LOS
      l_sob = y_th/fabs(Q); // definition of Sobolev length
    }
    // if (k==ks && j==js)
    //   cout <<  "r_l = " << r_l << ", l = " << l << ", l_sob = " << l_sob << endl;
    
    
    // recalculate emissivity at this position along LOS
    Real rho_l = calc_density(r_l, r_max);
    n_l = rho_l/mbar;
    xi_l = xi0/(n_l*r_l*r_l);
    T_l = calc_temperature(r_l);
    // kappa_l = calc_opacity(xi_l,T_l);

    // separately record escape probability contribution from local and non-local effects
    int idx_0or1 = (l < 1.5*l_sob) ? 0:1;

    
    
    // evaluate freq-dependent escape probability 
    if (v_l > vrange_min && v_l < vrange_max) {
      int idx;
      int ctr = 0;
      // loop over central bin as well as the neighboring bins
      for (int yidx=yidx_L; yidx<=yidx_R; yidx++) {
        if (i_edge == 0)
          idx = yidx;
        else // for edge cases, yidx is 0 or 1
          idx = i_edge - yidx;
        Real nu_min = nu0*(1.+L0bins->range[idx]);
        Real nu1 = nu_min + 0.5*dnu;
        
        // use trapezoid rule integration to evaluate p_nu 
        // for all frequencies in this bin
        for (int inu=0; inu < Nnu_per_bin; inu++) {
          Real nu = nu1 + inu*dnu;
          Real alpha_lm1 = calc_alpha(n_lm1,T_lm1,nuD_lm1,nu,SR,argb);
          Real alpha_l = calc_alpha(n_l,T_l,nuD_l,nu,SR,argb);
          Real delta_tau = (length_unit*dl)*0.5*(alpha_lm1 + alpha_l); // JAK added length unit to match alpha units correctly // Added relativisitic factor (nu/nu0) see notes
          


          // if (ctr==0 && k==ks && j==js) //delta_tau > 1e-3*dl) //k==ks && j==js
          //   cout <<  "inu = " << inu << ", arg = " << arg << ", nu/delta_nu0 = " << nu/delta_nu0
          //        <<  ", nuD/delta_nu0  = " << nuD/delta_nu0 <<  ", rho_lm1 = " << rho_lm1 
          //        <<  ", kappa_lm1 = " << kappa_lm1 << ", alpha_lm1 = " << alpha_lm1 
          //        << ", delta_tau = " << delta_tau
          //        << endl;
          // if (k==2 && j==25 && i==76 && inu==2) 
            // cout << "k = " << k << ", j = " << j << ", i = " << i << ", inu = " << inu << ", idx = " << idx << ", idx_0or1 = " << idx_0or1 << " , delta_tau = " << delta_tau << endl;
          
          
            Real p_nu = 1;

          if (r_l < r_max){
            p_nu = exp(-delta_tau);}
          
        

          pmb->pmy_mesh->ruser_mesh_data[idx_0or1](idx,inu) *= p_nu;

          // if (idx == yidx_L && inu == 0){ 
          //   tauW += delta_tau;} // Collect tau in the wing 

          // if (tauW > 5){ // if tau in the wing = 5 then escape probability is negligable for the whole frequency range. (terminate the LOS calculation from this point)
            
          //   // loop over central bin as well as the neighboring bins
          //   for (int yidx=yidx_L; yidx<=yidx_R; yidx++) {
          //     if (i_edge == 0){
          //       idx = yidx;}
          //     else { // for edge cases, yidx is 0 or 1
          //       idx = i_edge - yidx;}

          //     for (int inu=0; inu < Nnu_per_bin; inu++) {
          //         pmb->pmy_mesh->ruser_mesh_data[0](idx,inu) = 0; // Set local escape probability to zero
          //         pmb->pmy_mesh->ruser_mesh_data[1](idx,inu) = 0; // Set local escape probability to zero
          //         r_l = r_max+r_max; // and end while loop
          //         }}
                  
          //   tauW = 0;
          //   cout <<  "Terminate  i = " << i << ", j = " << j << ", k = " << k << endl;
          // }
            
          
        }
      ctr += 1;
      }
    }
    // update for next pass
    n_lm1 = n_l;
    T_lm1 = T_l; 
    nuD_lm1 = nuD_l;

    if (r_l > r_max)
      {pmb->ruser_meshblock_data[6](k,j,i) = l;}


  } // end while-loop 


  // Now that p_nu is known, compute L_nu = \int j_nu p_nu dV
  Real dV = vol(i)*pow(length_unit,3.);

  pmb->pmy_mesh->ruser_mesh_data[7](1) += dV;
  pmb->pmy_mesh->ruser_mesh_data[8](1) += 1;
  
  Real nu_sum = 0.;
  Real jnu_avg = 0.;
  Real jeff_avg = 0.;
  Real pnu_local_avg = 0.;
  Real pnu_nonlocal_avg = 0.;

  for (int yidx=yidx_L; yidx<=yidx_R; yidx++) {
    if (i_edge == 0)
      idx = yidx;
    else // for edge cases, yidx is 0 or 1
      idx = i_edge - yidx;
    Real nu_min = nu0*(1.+L0bins->range[idx]);
    Real nu1 = nu_min + 0.5*dnu;
    Real j_nu,j_nuD, j_nu1, nu_cmv, nu_cmv1, arg, arg1;
   
    if (use_nuD) // 4testing: evaluate all frequencies at line-center only
      j_nuD = calc_jnu(n_cgs,T_cgs,nuD,nuD,SR);  

    for (int inu=0; inu < Nnu_per_bin; inu++) {
      Real nu = nu1 + inu*dnu;
      
      // evaluate emissivity
      if (use_nuD) {
        j_nu = j_nuD;
        pmb->pmy_mesh->iuser_mesh_data[0](idx) += 1;
      }
      else 
        j_nu = calc_jnu(n_cgs,T_cgs,nuD,nu,SR); 
       
      // retrieve escape probability calculated in above while-loop
      Real pnu_local = pmb->pmy_mesh->ruser_mesh_data[0](idx,inu);
      Real pnu_nonlocal = pmb->pmy_mesh->ruser_mesh_data[1](idx,inu);
      Real p_nu = pnu_local*pnu_nonlocal;
      // std::cout << "p_nu = " << p_nu << std::endl;// added JAK
      // add to luminosity
      pmb->pmy_mesh->ruser_mesh_data[2](idx,inu) += j_nu*p_nu*dV;
      pmb->pmy_mesh->ruser_mesh_data[3](idx,inu) += j_nu*dV;
      pmb->pmy_mesh->ruser_mesh_data[6](idx,inu) += dV;
      pmb->pmy_mesh->ruser_mesh_data[9](idx,inu) += 1;
      // pmb->ruser_meshblock_data[7](k,j,i,idx,inu) = nu_cmv;
      // pmb->ruser_meshblock_data[8](k,j,i,idx,inu) = nu_cmv1;
      // pmb->ruser_meshblock_data[9](k,j,i,idx,inu) = nuD;
      // pmb->ruser_meshblock_data[10](k,j,i,idx,inu) = nuD1;
      // pmb->ruser_meshblock_data[11](k,j,i,idx,inu) = nu;
      // pmb->ruser_meshblock_data[12](k,j,i,idx,inu) = v_los;
      // pmb->ruser_meshblock_data[13](k,j,i,idx,inu) = v_;
      // pmb->ruser_meshblock_data[14](k,j,i,idx,inu) = arg;
      // pmb->ruser_meshblock_data[15](k,j,i,idx,inu) = arg1;
  
      // compute freq-integrated p_nus and j_nus 
      nu_sum += nu;
      jnu_avg += nu*j_nu;
      jeff_avg += nu*j_nu*p_nu;
      pnu_local_avg += nu*pnu_local;
      pnu_nonlocal_avg += nu*pnu_nonlocal;
      
    }
  }

  // finish computing the freq-averaged values and store them
  pnu_local_avg /= nu_sum;
  pnu_nonlocal_avg /= nu_sum;
  jnu_avg /= nu_sum;
  jeff_avg /= nu_sum;
  pmb->ruser_meshblock_data[0](k,j,i) = v_los;
  pmb->ruser_meshblock_data[1](k,j,i) = dV;
  pmb->ruser_meshblock_data[4](k,j,i) = ctr;
  pmb->ruser_meshblock_data[5](k,j,i) = log(pnu_local_avg*pnu_nonlocal_avg);
  pmb->ruser_meshblock_data[2](k,j,i) = jnu_avg;
  pmb->ruser_meshblock_data[3](k,j,i) = jeff_avg;


  // pmb->ruser_meshblock_data[0](k,j,i) = pnu_local_avg;
  // pmb->ruser_meshblock_data[1](k,j,i) = pnu_nonlocal_avg;
  // pmb->ruser_meshblock_data[4](k,j,i) = ctr;
  // pmb->ruser_meshblock_data[5](k,j,i) = log(pnu_local_avg*pnu_nonlocal_avg);
  // pmb->ruser_meshblock_data[2](k,j,i) = jnu_avg;
  // pmb->ruser_meshblock_data[3](k,j,i) = jeff_avg;
}


//--------------------------------------------------------------------------------------
//! \fn void MeshBlock::UserWorkInLoop(void)
//  \brief User-defined work function for every time step
void MeshBlock::UserWorkInLoop(void) {}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
{
  Real time = pmy_mesh->time;
  cout << "\n[MeshBlock::UserWorkBeforeOutput] ..." << time << "(T)" << endl;

  for(int k=ks; k<=ke; k++) {
    for(int j=js; j<=je; j++) {
      for(int i=is-NGHOST; i<=ie+NGHOST; i++) {

        user_out_var(0,k,j,i) = ruser_meshblock_data[0](k,j,i); 
        user_out_var(1,k,j,i) = ruser_meshblock_data[1](k,j,i);
        user_out_var(2,k,j,i) = ruser_meshblock_data[2](k,j,i);
        user_out_var(3,k,j,i) = ruser_meshblock_data[3](k,j,i);

      }
    }
  }

}

void Mesh::UserWorkInLoop(void){}

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  
  gsl_histogram_free(L0bins);
  cout << "\n[Mesh::UserWorkAfterLoop] Done." << endl;

  return;
}

string default_filename(string label)
{
  string name = label; 
  //append Nbins and Nnu_per_bin
  name += "_" + to_string<int>(Nbins) + "binsx" + to_string<int>(Nnu_per_bin);

  name += ".txt";
  return name;
}

int write_NbinsxNnu_dump(string fname, Mesh *pm, int n)
{
  int status;

  FILE *stream;
  std::stringstream msg;
  if ((stream = fopen(fname.c_str(),"w+")) == nullptr) {
      msg << "### FATAL ERROR in function [Mesh::UserWorkAfterLoop]"
          << std::endl << "Error output file could not be opened" <<std::endl;
      ATHENA_ERROR(msg);
    }
  // status = fprintf (stream, "%e", pm->ruser_mesh_data[1](0,0));
  
  for (int i = 0; i < Nbins; i++)
    for (int idx = 0; idx < Nnu_per_bin; idx++)
    {
      status = fprintf (stream, "%e", pm->ruser_mesh_data[n](i,idx));
      status = putc ('\n', stream);
    }

  fclose(stream);

  return 0;
}


int write_Nbins_dump(string fname, Mesh *pm, int n)
{
  int status;

  FILE *stream;
  std::stringstream msg;
  if ((stream = fopen(fname.c_str(),"w+")) == nullptr) {
      msg << "### FATAL ERROR in function [Mesh::UserWorkAfterLoop]"
          << std::endl << "Error output file could not be opened" <<std::endl;
      ATHENA_ERROR(msg);
    }
  // status = fprintf (stream, "%e", pm->ruser_mesh_data[1](0,0));
  
  for (int i = 0; i < Nbins; i++) {
    status = fprintf (stream, "%e", pm->ruser_mesh_data[n](i));
    status = putc ('\n', stream);
  }

  fclose(stream);

  return 0;
}


int write_bins(string fname)
{
  int status;

  FILE *stream;
  std::stringstream msg;
  if ((stream = fopen(fname.c_str(),"w+")) == nullptr) {
      msg << "### FATAL ERROR in function [Mesh::UserWorkAfterLoop]"
          << std::endl << "Error output file could not be opened" <<std::endl;
      ATHENA_ERROR(msg);
  }

  
  for (int idx = 0; idx < Nbins; idx++)
  {
    status = fprintf (stream, "%e", L0bins->range[idx]) ;
    status = putc ('\n', stream);
  }

  fclose(stream);

  return 0;
}





int write_Pnu_dump(string fname, MeshBlock *pmb, int n, int ie, int je, int ke, int is, int js, int ks)
{
  int status;

  FILE *stream;
  std::stringstream msg;
  if ((stream = fopen(fname.c_str(),"w+")) == nullptr) {
      msg << "### FATAL ERROR in function [write_Pnu_dump]"
          << std::endl << "Error output file could not be opened" <<std::endl;
      ATHENA_ERROR(msg);
    }
  // status = fprintf (stream, "%e", pm->ruser_mesh_data[1](0,0));
  
  for (int k = ks; k <= ke; k++)
    for (int j = js; j <= je; j++)
      for (int i =is; i <= ie; i++)
    {
      status = fprintf (stream, "%e", pmb->ruser_meshblock_data[n](k,j,i));
      status = putc ('\n', stream);
    }
  
  // pmb->ruser_meshblock_data[0](k,j,i) = pnu_local_avg;
  // pmb->ruser_meshblock_data[1](k,j,i) = pnu_nonlocal_avg;
  fclose(stream);

  return 0;
}

// ruser_meshblock_data[7](k,j,i,idx,inu)

int write_nu_cmv_dump(string fname, MeshBlock *pmb, int n, int ie, int je, int ke, int is, int js, int ks, int Nbins, int Nnu_per_bin)
{
  int status;

  FILE *stream;
  std::stringstream msg;
  if ((stream = fopen(fname.c_str(),"w+")) == nullptr) {
      msg << "### FATAL ERROR in function [write_Pnu_dump]"
          << std::endl << "Error output file could not be opened" <<std::endl;
      ATHENA_ERROR(msg);
    }
  // status = fprintf (stream, "%e", pm->ruser_mesh_data[1](0,0));
  
  for (int k = ks; k <= ke; k++)
    for (int j = js; j <= je; j++)
      for (int i =is; i <= ie; i++)
        for (int l =0; l <= Nbins; l++)
          for (int m =0; m <= Nnu_per_bin; m++)
            {
              status = fprintf (stream, "%e", pmb->ruser_meshblock_data[n](k,j,i,l,m));
              status = putc ('\n', stream);
              // if (l == Nbins)
              //   if (m == Nnu_per_bin)
              //   {
              //     cout <<  "Terminate  i = " << i << ", j = " << j << ", k = " << k << ", l =" << l << ", m =" << m << ", nu_cmv = " << pmb->ruser_meshblock_data[n](k,j,i,l,m) << endl;
              //   }
            }
  
  // pmb->ruser_meshblock_data[0](k,j,i) = pnu_local_avg;
  // pmb->ruser_meshblock_data[1](k,j,i) = pnu_nonlocal_avg;
  fclose(stream);

  return 0;
}

int write_TotV_dump(string fname, Mesh *pm, int n)
{
  int status;

  FILE *stream;
  std::stringstream msg;
  if ((stream = fopen(fname.c_str(),"w+")) == nullptr) {
      msg << "### FATAL ERROR in function [Write_TotV_dump]"
          << std::endl << "Error output file could not be opened" <<std::endl;
      ATHENA_ERROR(msg);
    }
  // status = fprintf (stream, "%e", pm->ruser_mesh_data[1](0,0));
  
  for (int i = 0; i < 2; i++) {
    status = fprintf (stream, "%e", pm->ruser_mesh_data[n](i));
    status = putc ('\n', stream);
  }

  fclose(stream);

  return 0;
}












//  Added to output the escape probabliity Jan 9 2025 JAK

// int write_ruser_meshblock_data(string fname, MeshBlock *pmb, int n)
// {
  
//   int status;

//   FILE *stream;
//   std::stringstream msg;

//   // Open the file for writing
//   if ((stream = fopen(fname.c_str(), "w+")) == nullptr) {
//     msg << "### FATAL ERROR in function [write_ruser_meshblock_data]"
//         << std::endl << "Error: output file could not be opened" << std::endl;
//     ATHENA_ERROR(msg);
//   }

//   // Loop through 3D indices (k, j, i) to write data
//   for (int k = 0; k < pmb->ke; ++k) {  // Replace `ke` with the actual size of k-dimension
//     for (int j = 0; j < pmb->je; ++j) {  // Replace `je` with the actual size of j-dimension
//       for (int i = 0; i < pmb->ie; ++i) {  // Replace `ie` with the actual size of i-dimension
//         // Write the value of pmb->ruser_meshblock_data[n](k, j, i)
//         status = fprintf(stream, "%e,%e,%e:%e\n", k,j,i,pmb->ruser_meshblock_data[n](k, j, i));
//         if (status < 0) {
//           msg << "### ERROR in function [write_ruser_meshblock_data]"
//               << std::endl << "Error: failed to write data to file" << std::endl;
//           ATHENA_ERROR(msg);
//         }
//       }
//     }
//   }

//   fclose(stream); // Close the file

//   return 0; // Return success
// }
