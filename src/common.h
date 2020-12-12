/*	SET OF GLOBAL VARIABLES, CONSTANTS AND PARAMETERS COMMON
 *	BETWEEN DIFFERENT MODULES.														*/

// number of chemical species and PDEs to be solved
const int CTYPES = 23, NPDE = 21;
//const int CTYPES = 25, NPDE = 23;
const double CC_FACTOR = 1.0;	// factor for chemical reaction parameters

const double lscale = 1e-6, tscale = 2.27e-4, fscale = 4.5e-14, escale = 4.5e-20;
const double accel = 100; // time-scale acceleration for the reaction rates, 100 default
const double tscaleCC = accel*tscale;

// GPIba-vWF flex-bond parameters
const double Ftrans = 1e-11/fscale, alpha = 1.0+0.05/2.25, kf01 = 5.8e0*tscale, kf02 = 52e0*tscale,
		         sigf1 = 1.58e-9/lscale, sigf2 = 1.31e-9/lscale, dgf1 = 6.17e-21/escale, dgf2 = 6.17e-21/escale,
						 kr01 = 0.0047*tscale, kr02 = 0.0022*tscale, sigr1 = 2.52e-9/lscale, sigr2 = 1.62e-9/lscale;

// fibrin clot porosity parameters
const double fib_pow = -1.8, fib_conv = 100.0, fib_thr = 4000.0;

// solid/bound and surface reactions parameters
const double RAD = 21.0; // vessel radius
const double DX = 195.0, DY = 20.0, DZ = 20.0; // channel dimensions
const double D_ksi = 3.55e-4, D_Fick = 0.0682, TFVIIa = 2000e-6;// default TFVIIa = 50e-4, correct scaling is TFVIIa = 200e-6
const double Linj = 30.; // injury length
const int n_sur_reac = 7, nspec = 5;
const int spec[nspec] = {1, 7, 8, 10, 13};
const int sur_reac[n_sur_reac] = {0, 1, 6, 7, 12, 13, 19};
const double vals[10] = {0.034*tscaleCC, //phi_11
						 2000.,  //Phi_11M
						 32.4*tscaleCC,  //k_7,9
						 24.,  //K_7,9M
						 103.*tscaleCC, //k_7,10
						 240., //K_7,10M
						 6.52e-13*tscaleCC*1.0e12,  //k_tPA^C
						 9.27e-12*tscaleCC*1.0e12,  //k_tPA^IIa
						 5.059e-18*tscaleCC*1.0e12, //k_tPA^Ia
						 134.8*tscaleCC
						},
				cw[3] = {TFVIIa, //0.025 [TF-VIIa]
						 2.0e9*1.0e-12, //[ENDO]
						 0.0 //375 [XIIa]
					    };

// ADP release and platelet activation
const double IIa_th = 1.0, ADP_th = 1000.0;
const double reltau = 5e-2/tscaleCC, relcon = 1.2e2;
