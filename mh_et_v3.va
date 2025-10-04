`include "constants.vams"
`include "disciplines.vams"

// MarcusHush asymmetrical surface-confined electron-transfer model
// Acts as a voltage-controlled current source: I = f(V)
// Ported for Verilog-A, using an exact integrator for unconditional stability
// Includes optional parameters like a gain scaling factor and offset for flexibility

module mh_et_v3(p, n);
  inout p, n;
  electrical p, n, cov;

  // 1) PARAMETERS & CONSTANTS (all literals)
  parameter real n_e       = 2.0;                 // electrons transferred
  parameter real F_const   = 96485.0;             // C/mol
  parameter real Rgas      = 8.314;               // J/molK
  parameter real T         = 298.15;              // K
  parameter real E0        = -0.34;               // V formal potential
  parameter real Vstart    = -0.6;                // V initial baseline
  parameter real k0_red    = 1000;                 // s preexponential
  parameter real k0_ox     = 1000;
  parameter real lambda_eV_red = 0.2;                 // eV reorganization energy
  parameter real lambda_eV_ox  = 0.2;
  parameter real Gamma_tot = 16e-9;                // mol/m total coverage
  parameter real PI        = 3.141592653589793;   // 
  parameter real gain      = 200;
  parameter real Inonfard  = 60e-9;
  //parameter real Inonfard  = 0;
  parameter real Ioffset   = 0e-9;
  parameter real Fswv      = 0;

  // 2) DERIVED CONSTANTS & STATE (all declared at module scope)
  real    lambda_red, lambda_ox, beta_red, beta_ox, A;     // computed in analog
  real    Gamma_init;          // initial coverage
  real    Vapp, eta;           // potentials
  real    k_red, k_ox;         // Marcus rates
  real    dGdt;                // d/dt
  real    sw;

  analog begin
    // 3) Compute derived constants at runtime
    lambda_red = lambda_eV_red * F_const;
    lambda_ox  = lambda_eV_ox * F_const;                
    beta_red   = 1.0/(4.0*lambda_red*Rgas*T);
    beta_ox    = 1.0/(4.0*lambda_ox*Rgas*T);            
    A          = PI*(250e-6)**2;                     

    // 4) Preequilibrate  at t=0
    @(initial_step) begin : init_blk
      real eta0, kr0, kox0;
      eta0       = Vstart - E0;
      kr0        = k0_red * exp(-((lambda_red + n_e*F_const*eta0)**2) * beta_red);
      kox0       = k0_ox * exp(-((lambda_ox - n_e*F_const*eta0)**2) * beta_ox);
      Gamma_init = Gamma_tot * kr0/(kr0 + kox0);
      V(cov, n)  <+ Gamma_init;   // set initial coverage node voltage
    end

    // 5) Compute overpotential & rates each timestep
    Vapp    = V(p, n);
    eta     = Vapp - E0;
    k_red   = k0_red * exp(-((lambda_red + n_e*F_const*eta)**2) * beta_red);
    k_ox    = k0_ox * exp(-((lambda_ox - n_e*F_const*eta)**2) * beta_ox);
    dGdt    = k_red*(Gamma_tot - V(cov,n))
            - k_ox * V(cov,n);

    // 6) Integrate  via KCL on node 'cov'
    I(cov, n) <+ ddt( V(cov, n) );  // I_C = CdV/dt, with C=1 F
    I(cov, n) <+ -dGdt;             // enforce dV/dt = dGdt

    // 7) Faradaic current source: I = -nFAd/dt
    I(p, n)   <+ -n_e * F_const * A * dGdt * gain;

    // 8) Non-faradic current
	sw = sin(2.0*`M_PI*Fswv*$abstime) > 0.0 ?  1.0 : -1.0;
	I(p,n) <+ Inonfard * sw;
    I(p,n) <+ Ioffset;
  end
endmodule
