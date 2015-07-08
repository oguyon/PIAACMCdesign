\page PIAACMC_designcodes PIAACMC design codes


# DESIGN PARAMETERS, AND HOW/WHERE TO SET THEM

## Variables driving PIAACMC optics design, excluding focal plane mask


### Non-optional parameters

Parameters with "*" appear in directory name.

Variable Name          | Description                                      | Setting location
-----------------------|--------------------------------------------------|-----------------------------------
pupil                  | pupil geometry                                   | File ./pup_1024.fits (name defined in script sim1024)
(2) coin *             | remapping input central obstruction              | argument #2 of ./runopt script -> passed to sim1024 -> PIAACMC_centobs0 in ./runPIAACMC
(3) coout *            | output central obstruction                       | argument #3 of ./runopt script -> passed to sim1024 -> PIAACMC_centobs1 in ./runPIAACMC
(4) fpmrad *           | PIAACMC design nominal focal plane mask radius   | argument #4 of ./runopt script -> passed to sim1024 -> PIAACMC_fpmradld in ./runPIAACMC
(5) lambda *           | design wavelength [nm]                           | argument #5 of ./runopt script -> passed to sim1024 -> PIAACMC_lambda in ./runPIAACMC
(6) PIAAmaterial *     | material for PIAA optics (Mirror, CaF2, etc...)  | argument #6 of ./runopt script -> passed to sim1024 -> File <confdir>/conf_PIAAmaterial_name.txt
(7) LStransm *         | Lyot stops geometric transmission                | argument #7 of ./runopt script -> passed to sim1024 -> PIAACMC_LStransm0, PIAACMC_LStransm1, PIAACMC_LStransm2 in ./runPIAACMC
(8) NBlyotstop *       | Number of Lyot Stops                             | argument #8 of ./runopt script -> passed to sim1024 -> PIAACMC_nblstop in ./runPIAACMC

### Optional parameters

The index number in the directory name is typically used to keep track of these parameters.

Variable Name          | Description                                      | Setting location
-----------------------|--------------------------------------------------|-----------------------------------
size                   | array size [pix]                                 | File ./conf_size.txt           (optional, default = 1024), read by ./runopt script
MdesignStepMax         | Maximum design step (monochromatic)              | File ./conf_MdesignStepMax.txt (optional, default = 13, max = 18), read by ./runopt script
beamrad                | Beam physical radius [m]                         | File ./conf_PIAAbeamrad.txt    (optional, default = 0.01), read by ./runPIAACMC script, PIAACMC_beamrad
Fratio                 | Focal ratio at focal plane                       | File ./conf_Fratio.txt         (optional, default = 80), read by ./runPIAACMC, PIAACMC_Fratio
r0lim                  | outer radius of first PIAA optic                 | File ./conf_PIAAr0lim.txt      (optional, default = 1.15), read by ./runPIAACMC, PIAACMC_r0lim
r1lim                  | outer radius of second PIAA optic                | File ./conf_PIAAr1lim.txt      (optional, default = 1.50), read by ./runPIAACMC, PIAACMC_r1lim
piaasep                | separation between PIAA elements [m]             | File ./conf_PIAAsep.txt        (optional, default = 1.00), read by ./runPIAACMC, PIAACMC_piaasep
piaa0pos               | conjugation of first PIAA element [m]            | File ./conf_PIAA0pos.txt       (optional, default = 1.00), read by ./runPIAACMC, PIAACMC_piaa0pos

## Variables driving focal plane mask design




Variable Name          | Description                                      | Setting location
-----------------------|--------------------------------------------------|-----------------------------------
(9)  mlambda *         | mask central wavelength [nm]                     | argument #9 of ./runopt script -> passed to sim1024 -> PIAACMC_lambda in ./runPIAACMC
(10) mlambdaB *        | mask spectral bandwidth [%]                      | argument #10 of ./runopt script -> passed to sim1024 -> PIAACMC_lambdaB in ./runPIAACMC
(11) NBrings *         | Number of focal plane mask rings                 | argument #11 of ./runopt script -> passed to sim1024 -> PIAACMC_NBrings in ./runPIAACMC
(12) maskradld *       | outer radius of focal plane mask rings           | argument #12 of ./runopt script -> passed to sim1024 -> PIAACMC_MASKRADLD in ./runPIAACMC
(13) PIAACMC_resolved *| source size (10x (-log(source rad)), 00 if pt    | argument #13 of ./runopt script -> passed to sim1024 -> PIAACMC_resolved in ./runPIAACMC
(14) PIAACMC_extmode * | how to simulate extended source (0:3pts, 1:6pts) | argument #14 of ./runopt script -> passed to sim1024 -> PIAACMC_extmode in ./runPIAACMC
(15) fpmmaterial *     | material for focal plane mask                    | argument #15 of ./runopt script -> passed to sim1024 -> File <confdir>/conf_fpmmaterial_name.txt





# FILES AND DIRECTORIES SYNTAX



PIAACMC designs are described by multipe parameters. Some are contained in the directory name in which the design is stored (mostly related to mirror design), and some are specific to focal plane masks and thus containted in the name of focal plane mask, as detailed below. 

## DIRECTORIES 

### MONOCHROMATIC SEED, CIRCULAR SYMMETRIC

All designs start from a monochromatic "seed" design, which is an idealized monochromatic PIAACMC design, with PIAA optics modeled as infinitely thin flat OPD screens.

\verbatim
piaacmcconf_<x>_coin<x.xx>_coout<x.xx>_fpmr<x.xx>

x    : focal plane scoring region (set to 1) 
coin : input central obstruction
coout: output central obstruction
fpmr : design focal plane mask radius
\endverbatim



### OPTIMIZED DESIGN, MONO- AND POLY-CHROMATIC 

\verbatim
piaacmcconf_<x>_coin<x.xxx>_coout<x.xxx>_fpmr<x.xxx>_l<xxxx>_<mmmmmm>_lt<t.tt>_ls<N>_i<xxx>

x    : focal plane scoring region (set to 1) 
coin : input central obstruction
coout: output central obstruction
fpmr : design focal plane mask radius
l    : central wavelength [nm]
mmmmm: material used for PIAA optics (corresponding code is PIAAmaterial in OPTPIAACMCDESIGN, code written as conf_PIAAmaterial.txt in directory)
lt   : Lyot stops geometric transmission [%]
ls   : Number of Lyot Stops
i    : design index (used to number designs beyond parameters above)
\endverbatim

## FOCAL PLANE MASK DESIGNS

### Design solution

\verbatim
fpm_zonez_s<s>_l<xxxx>_sr<bb>_nbr<zzz>_mr<rrr>_ssr<ss>_ssm<m>_wb<ll>.fits


s    : PIAACMC_FPMsectors (usually 1)
xxxx : central mask wavelength [nm]
bb   : spectral bandwidth [%]
zzz  : number of rings
rrr  : 100*PIAACMC_MASKRADLD
ss   : computePSF_ResolvedTarget: stellar radius used for optimization  [10x log l/D]
	radius = 10^{-0.1*ss}
	example: ss = 15  ->  radius = 0.0316 l/D
m    : computePSF_ResolvedTarget_mode
	0: 3 points
	1: 6 points
ll   : number of wavelength bins


sprintf(fname, "!%s/fpm_zonez_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_ssr%02d_ssm%d_%s_wb%02d.fits", piaacmcconfdir, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);

\endverbatim

### Linear response between focal plane mask zones sags and final focal plane complex amplitude

\verbatim
FPMresp<m>_s<s>_l<xxxx>_sr<bb>_nbr<zzz>_mr<rrr>_ssr<ss>_ssm<m>_wb<ll>.fits

m    : Scoring mask type
s    : PIAACMC_FPMsectors (usually 1)
xxxx : central mask wavelength [nm]
bb   : spectral bandwidth [%]
zzz  : number of rings
rrr  : 100*PIAACMC_MASKRADLD
ss   : computePSF_ResolvedTarget: stellar radius used for optimization  [10x log l/D]
	radius = 10^{-0.1*ss}
	example: ss = 15  ->  radius = 0.0316 l/D
m    : computePSF_ResolvedTarget_mode
	0: 3 points
	1: 6 points
ll   : number of wavelength bins


sprintf(fname, "!%s/FPMresp%d_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_ssr%02d_ssm%d_%s_wb%02d.fits", piaacmcconfdir, SCORINGMASKTYPE, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);

\endverbatim


### Evaluation files

\verbatim
Point source PSF:
psfi0_ptsr_sm<m>_s<s>_l<xxxx>_sr<bb>_nbr<zzz>_mr<rrr>_ssr<ss>_ssm<m>_wb<ll>.fits

Extended source PSF (4 points):
psfi0_extsrc<exss>_sm<m>_s<s>_l<xxxx>_sr<bb>_nbr<zzz>_mr<rrr>_ssr<ss>_ssm<m>_wb<ll>.fits

Flux: 
flux_sm<m>_s<s>_l<xxxx>_sr<bb>_nbr<zzz>_mr<rrr>_ssr<ss>_ssm<m>_wb<ll>.txt

Contrast estimate vale (point source):
contrast_sm<m>_s<s>_l<xxxx>_sr<bb>_nbr<zzz>_mr<rrr>_ssr<ss>_ssm<m>_wb<ll>.txt

Extended source contrast:
ContrastVal_extsrc<exss>_sm<m>_s<s>_l<xxxx>_sr<bb>_nbr<zzz>_mr<rrr>_ssr<ss>_ssm<m>_<material>_wb<ll>.fits
(2nd number is contrast estimate)

Extended source contrast curve:
ContrastCurve_extsrc<exss>_sm<m>_s<s>_l<xxxx>_sr<bb>_nbr<zzz>_mr<rrr>_ssr<ss>_ssm<m>_<material>_wb<ll>.fits


exss : stellar radius used for evaluation  [10x log l/D]
m    : Scoring mask type
s    : PIAACMC_FPMsectors (usually 1)
xxxx : central mask wavelength [nm]
bb   : spectral bandwidth [%]
zzz  : number of rings
rrr  : 100*PIAACMC_MASKRADLD
ss   : computePSF_ResolvedTarget: stellar radius used for optimization  [10x log l/D]
	radius = 10^{-0.1*ss}
	example: ss = 15  ->  radius = 0.0316 l/D
m    : computePSF_ResolvedTarget_mode
	0: 3 points
	1: 6 points
material: focal plane mask material
ll   : number of wavelength bins


sprintf(fname, "!%s/psfi0_sm%d_s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_ssr%02d_ssm%d_%s_wb%02d.fits", piaacmcconfdir, SCORINGMASKTYPE, PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*PIAACMC_MASKRADLD+0.1), computePSF_ResolvedTarget, computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name, piaacmc[0].nblambda);
\endverbatim






# OPERATION CODES


\verbatim



====== CODE INDICES FOR sim1024/runPIAACMC SCRIPTS =====================================================

The indices are used to instruct which computation step(s) shall be performed
Unless described otherwise, columns in this table are:
- ###     code index/step
- [##]   corresponding index in C function PIAACMCsimul_exec
- ....   Description
- [...]  approximate execution time






--------- 0-99  Monochromatic, point source steps ---------------------------------------------------

Code runs from step 0 to requested step, skipping steps that have already been completed
Note: entering step=4 will complete steps 0, 1, 2 and 3, and STOP at step 4

  0     [  0]   3.43696e-07	PREPARE - Compute on-axis PSF, monochromatic, idealized PIAACMC, first iteration  [   8mn]
		function PIAAsimul_initpiaacmcconf
		{
			create/load Cmodes_1024.fits & Fmodes_1024.fits (used to fit PIAA shapes)

			CREATE GENERALIZED PROLATE
			Uses DFT iterations
			-> APLCapo.0.800.0.300.info  (masksizeld 0.0 prolatethroughput peak)
			-> <confdir>/piaaref/APLCmaskCtransm.txt : complex amplitude mask transmission (0.51005532508) = (1.0-peak)/peak
			-> <confdir>/piaaref/apo2Drad.fits : 2D apodization (also copied on <confdir>/apo2Drad.fits)
			Assumes circular aperture with circular central obstruction
					
			FIT 2D RADIAL APODIZATION AS SUM OF 10 COSINES		
			-> <confdir>/APOmodesCos.fits : modes used for fitting (2D maps)		
			COMPUTE 2D SAG MAPS
			Uses the cosine fit above to make 1D radial profile of PIAA shapes
			-> <confdir>/PIAA_Mshape.txt (r0 z0 r1 z1)
			Make 2D sag maps
			-> <confdir>/piaam0z.fits, <confdir>/piaam1z.fits

			FIT 2D PIAA SHAPES AS SUM OF COSINES
			-> <confdir>/piaa0Cz.fits, <confdir>/piaa1Cz.fits (fitted shapes)
			-> <confdir>/piaa0Cres.fits, <confdir>/piaa1Cres.fits (residuals)
			FIT 2D RESIDUAL AS SUM OF FOURIER MODES

			NOMINAL PIAA SHAPES, FITTED TO COSINES AND FOURIER MODES
			-> <confdir>/piaaref/piaa0Cmodes.fits
			-> <confdir>/piaaref/piaa1Cmodes.fits
			-> <confdir>/piaaref/piaa0Fmodes.fits
			-> <confdir>/piaaref/piaa1Fmodes.fits
			+ copy to <confdir>/
		
			CREATE FOCAL PLANE MASK
			-> <confdir>/fpm_zonea00_00_001_001.fits	
			-> <confdir>/fpm_zonez00_00_001_001.fits	

			MAKE LYOT STOPS
			-> <confdir>/LyotStop0.fits, <confdir>/LyotStop1.fits ...
		}
		
		MAKE PIAA SHAPES FROM COEFFICIENTS (function PIAACMCsimul_makePIAAshapes())
		-> <confdir>/piaa0Cz.fits
		-> <confdir>/piaa0Fz.fits
		-> <confdir>/piaa0z.fits (sum of both terms above)
		-> <confdir>/piaa1Cz.fits
		-> <confdir>/piaa1Fz.fits
		-> <confdir>/piaa1z.fits (sum of both terms above)

		COMPUTE PSF
	       -> psfi0_step000.fits

  1     [  0]   3.43696e-07
		-> psfi0_step001.fits

  2     [   ]   5.12665e-05	Load pupil geometry                                                     [   0mn]
		-> pupa0_1024.fits

  3     [  0]   5.12665e-05	actual pupil Compute on-axis PSF                                        [   0nm]
		-> psfi0_step003.fits

  4     [  5]   3.46488e-06	Compute Lyot stops shapes and locations, 1st pass                       [   1mn]
		throughput = LStransm2 (user-specied argument to script)
		-> conjugations.txt

  5     [  2]   3.34238e-08	optimize focal plane mask transm, 1st pass -> result_fpmt.log           [   2mn]
		optimal value is written in file <confdir>/piaacmcparams.conf
		-> result_fpmt.log (ampl contrast iter range stepsize)
		

  6     [  5]   2.06295e-08	Compute Lyot stops shapes and locations, 2nd pass, 60% throughput       [   1mn]
		throughput > LStransm = LStransm0 = 0.60

  7     [ 40]    3.36149e-09    tune PIAA shapes and focal plane mask transm # modes: 10, 5             [  35mn]
		-> linoptval_step007.txt

  8     [ 40]   3.3526e-09	tune PIAA shapes and focal plane mask transm # modes: 20, 20            [  23mn]
		-> linoptval_step008.txt 
		updates piaa shapes (piaa0Cmodes, piaa0Fmodes, piaa1Cmodes, piaa1Fmodes)
		updates focal plane mask transm (<confdir>/piaacmcparams.conf)

  9     [  5]   4.3982e-09   	Compute Lyot stops shapes and locations, 2nd pass, intermediate throughput [   1mn]
		throughput = LStransm1

 10     [  1]   3.62923e-09	tune Lyot stops conjugations                                    	[   8mn]
		throughput > LStransm1 
		-> result_LMpos.log

 11     [ 40]   1.2318e-09	tune PIAA shapes and focal plane mask transm # modes: 20, 20            [  71mn]
		-> linoptval_step011.txt 

 12     [ 40]   2.9522e-10	tune PIAA shapes and focal plane mask transm, # modes: 40, 150          [ 287mn]

 13     [  5]   3.14043e-08	Compute Lyot stops shapes and locations, 3rd pass, goal throughput     	[   1mn]
 		throughput > LStransm2 
 		
 14     [  1]   1.41547e-08	Tune Lyot stops conjugations                                            [   9mn]
 		
 15     [ 40]   5.25497e-09	tune PIAA shapes and focal plane mask transm # modes: 20, 20            [  37mn]
 16     [ 40]    tune PIAA shapes and focal plane mask transm, # modes: 40, 150          	[ 169mn]
 17     [ 40]   tune PIAA shapes and focal plane mask transm, # modes: 40, 625          	[ 591mn]




--------- 100-199  Polychromatic, extended source, physical focal plane mask (consecutive) ------------------------

100	[ 11]   turn focal plane mask into zones, Compute polychromatic response to zones, store result in FPMresp
		-> FPMresp<SRE>_<xx>_<yyy>_<ll>.fits
			S:   SCORINGMASKTYPE
			R:   computePSF_ResolvedTarget
			E:   PIAACMC_FPMsectors
			xx: (long) (10.0*PIAACMC_MASKRADLD+0.1)  NOTE: this is the outer edge of the outer ring
			yyy: piaacmc[0].NBrings
			ll:  piaacmc[0].nblambda

101	[ 13]	linear piece-wise optimization - chromaticity
		<- Reads input mask from <confdir>/conf_MASKRADLD.txt <confdir>/conf_FPMsectors.txt, <confdir>/conf_NBrings.txt, <confdir>/conf_resolved.txt
		-> linoptval.txt
			local derivatives used to indentify steepest descent direction, then scan along this direction
			lines starting by "##" show values along line of steepest descent
			Best value then chosen, and new direction computed
		-> <confdir>/fpm_zonez<RE>_<rr>_<yyy>_<zzz>.fits,  <confdir>/fpm_zonea<RE>_<rr>_<yyy>_<zzz>.fits
			R:   computePSF_ResolvedTarget
			E:   PIAACMC_FPMsectors
			rr:  (long) (10.0*PIAACMC_MASKRADLD+0.1)	(=00 for idealized PIAACMC)	
			yyy: piaacmc[0].NBrings
			zzz: piaacmc[0].focmNBzone			
		-> mode13_<RE>_<yyy>_<zzz>_<ll>.txt
			R:   computePSF_ResolvedTarget
			E:   PIAACMC_FPMsectors			
			yyy: piaacmc[0].NBrings
			zzz: piaacmc[0].focmNBzone
		columns:
			1: iteration
			2: amplitude of random offset introduced at starting point (MODampl)
			3: starting point contrast (PIAACMCSIMUL_VALREF)
			4: current contrast value (PIAACMCSIMUL_VAL)
			5: best contrast value (bestval)
			6: zero starting poing flag (zeroST): 1 if starting from zero
		-> mode13_<RE>_<yyy>_<zzz>_<ll>.bestval.txt
			best contrast value
	will run until file <confdir>/stoploop13.txt detected	





----------- 200-299  Optimization (polychromatic) single steps ---------------------------------------

210	[40]	full co-optimization of PIAA shapes and focal plane mask zones
		# modes: 10, 5





--------- 500 - 599 : housekeeping -------------------------------------------------------------------

500	[N/A]	store as reference for use at other lambda / spectral band






------------ 700-799 EVALUATIONS ---------------------------------------------------------------------

A single step will be executed.

MONOCHROMATIC, PERFECT SINGLE ZONE MASK

700     [  0]   Evaluate point source mask constrat, on point source, idealized mask                   [ 12sec]
	        -> psfi0_step700.fits
		-> psfi0_%d%d_%02ld_%03ld_%03ld.fits
piaacmcconfdir, computePSF_ResolvedTarget, PIAACMC_FPMsectors, (long) (10.0*PIAACMC_MASKRADLD+0.1), piaacmc[0].NBrings, piaacmc[0].focmNBzone                                               

701     [100]   Evaluate sensitivity to pointing errors, idealized mask                                [ 25sec]
	        -> psfi0_step701.fits
		-> <confdir>/ContrastCurve%d_%03ld_%03ld_%02d_ps%03ld.txt"
PIAACMC_FPMsectors, piaacmc[0].NBrings, piaacmc[0].focmNBzone, piaacmc[0].nblambda, (long) (1000.0*ldoffset));
	

702     [   ]   Compute on-axis monochromatic PSF - extended source, using extended source focal plane mask



720	[  0]   Compute on-axis polychromatic PSF, 10 spectral values
		-> psfi0.fits, psfp0.fits, psfa0.fits		
		-> WFamp and WFpha FITS files
		-> fpm_pha.fits, fpm_amp.fits : 1 - focal plane mask complex amplitude
		-> lambdalist.txt
		-> psfi0_%d%d_%02ld_%03ld_%03ld.fits
piaacmcconfdir, computePSF_ResolvedTarget, PIAACMC_FPMsectors, (long) (10.0*PIAACMC_MASKRADLD+0.1), piaacmc[0].NBrings, piaacmc[0].focmNBzone

721	[100]   Evaluate sensitivity to pointing errors, polychromatic PSF
		-> psfi0.fits, psfp0.fits, psfa0.fits		
		-> WFamp and WFpha FITS files
		-> fpm_pha.fits, fpm_amp.fits : 1 - focal plane mask complex amplitude
		-> lambdalist.txt
		-> <confdir>/ContrastVall%d%d_%02ld_%03ld_%02d_tt%03ld.txt
computePSF_ResolvedTarget, PIAACMC_FPMsectors, (long) (10.0*PIAACMC_MASKRADLD+0.1), piaacmc[0].NBrings, piaacmc[0].nblambda, (long) (1000.0*ldoffset)
			col 1: valref = average flux over scoring region  
			col 2: aveC = averge contrast from 2 to 6 l/D
			col 3: mask rad l/D
			col 4: input Central Obs
			col 5: output Central Obs
			col 6: lambda [nm]
			col 7: bandwidth [%]
			col 8: sectors ?
			col 9: NB rings
			col 10: NB zones
			col 11: NB lambda
			col 12: source size offset
		-> <confdir>/ContrastCurve%d_%03ld_%03ld_%02d.txt"
PIAACMC_FPMsectors, piaacmc[0].NBrings, piaacmc[0].focmNBzone, piaacmc[0].nblambda);

750	[101]	Transmission curve



------------ MISC OPERATIONS -------------------------------------------------------------------------


800     [200]   Make focal plane mask OPD map -> FITS file
		-> <confdir>/fpmOPD<RE>_<yyy>_<zzz>.fits"
			R:   computePSF_ResolvedTarget
			E:   PIAACMC_FPMsectors
			yyy: piaacmc[0].NBrings
			zzz: piaacmc[0].focmNBzone



======================================================================================================

\endverbatim
