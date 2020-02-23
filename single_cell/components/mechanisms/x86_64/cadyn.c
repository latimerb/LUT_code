/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__cadyn
#define _nrn_initial _nrn_initial__cadyn
#define nrn_cur _nrn_cur__cadyn
#define _nrn_current _nrn_current__cadyn
#define nrn_jacob _nrn_jacob__cadyn
#define nrn_state _nrn_state__cadyn
#define _net_receive _net_receive__cadyn 
#define factors factors__cadyn 
#define parms parms__cadyn 
#define state state__cadyn 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define vmaxsr _p[0]
#define kpsr _p[1]
#define jmaxsr _p[2]
#define kip3 _p[3]
#define kactip3 _p[4]
#define konip3 _p[5]
#define kinhip3 _p[6]
#define kcicr _p[7]
#define ktcicr _p[8]
#define vcicr _p[9]
#define Kmer _p[10]
#define Bmer _p[11]
#define vmcu _p[12]
#define kmcu _p[13]
#define nmcu _p[14]
#define vncx _p[15]
#define kna _p[16]
#define kncx _p[17]
#define Kmmt _p[18]
#define Bmmt _p[19]
#define k1 _p[20]
#define k2 _p[21]
#define k3 _p[22]
#define k4 _p[23]
#define pump0 _p[24]
#define ica_pmp _p[25]
#define jmcu (_p + 26)
#define jmncx (_p + 38)
#define jmito (_p + 50)
#define jer (_p + 62)
#define jip3 (_p + 74)
#define jserca (_p + 86)
#define jcicr (_p + 98)
#define fmer (_p + 110)
#define fmmt (_p + 122)
#define pump _p[134]
#define pumpca _p[135]
#define ca (_p + 136)
#define caer (_p + 148)
#define camt (_p + 160)
#define hc (_p + 172)
#define ho (_p + 184)
#define Ln (_p + 196)
#define ica _p[208]
#define cai _p[209]
#define cao _p[210]
#define caeri _p[211]
#define camti _p[212]
#define last_ica_pmp _p[213]
#define parea _p[214]
#define c1 _p[215]
#define c2 _p[216]
#define c3 _p[217]
#define c4 _p[218]
#define ip3i _p[219]
#define nai _p[220]
#define Dpump _p[221]
#define Dpumpca _p[222]
#define Dca (_p + 223)
#define Dcaer (_p + 235)
#define Dcamt (_p + 247)
#define caip3ri _p[259]
#define Dcaip3ri _p[260]
#define Dhc (_p + 261)
#define Dho (_p + 273)
#define DLn (_p + 285)
#define v _p[297]
#define _g _p[298]
#define _ion_ica	*_ppvar[0]._pval
#define _ion_cai	*_ppvar[1]._pval
#define _ion_cao	*_ppvar[2]._pval
#define _ion_dicadv	*_ppvar[3]._pval
#define _style_ca	*((int*)_ppvar[4]._pvoid)
#define _ion_nai	*_ppvar[5]._pval
#define _ion_caeri	*_ppvar[6]._pval
#define _style_caer	*((int*)_ppvar[7]._pvoid)
#define _ion_camti	*_ppvar[8]._pval
#define _style_camt	*((int*)_ppvar[9]._pvoid)
#define _ion_caip3ri	*_ppvar[10]._pval
#define _style_caip3r	*((int*)_ppvar[11]._pvoid)
#define _ion_ip3i	*_ppvar[12]._pval
#define diam	*_ppvar[13]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_factors(void);
 static void _hoc_parms(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_cadyn", _hoc_setdata,
 "factors_cadyn", _hoc_factors,
 "parms_cadyn", _hoc_parms,
 0, 0
};
 #define _zfactors_done _thread[2]._pval[0]
 #define _zfrat (_thread[2]._pval + 1)
 #define _zdsq _thread[2]._pval[13]
 #define _zdsqvol _thread[2]._pval[14]
 #define _zdsqvolmt _thread[2]._pval[15]
 #define _zdsqvoler _thread[2]._pval[16]
 /* declare global and static user variables */
 static int _thread1data_inuse = 0;
static double _thread1data[12];
#define _gth 3
#define DCa DCa_cadyn
 double DCa = 0.6;
#define L L_cadyn
 double L = 0;
#define bbr bbr_cadyn
 double bbr = 370;
#define cao0 cao0_cadyn
 double cao0 = 2;
#define camti0 camti0_cadyn
 double camti0 = 0.0002;
#define caeri0 caeri0_cadyn
 double caeri0 = 0.4;
#define cai0 cai0_cadyn
 double cai0 = 0.000136;
#define vol_cadyn (_thread1data + 0)
#define vol (_thread[_gth]._pval + 0)
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "L_cadyn", "um",
 "cai0_cadyn", "mM",
 "cao0_cadyn", "mM",
 "caeri0_cadyn", "mM",
 "camti0_cadyn", "mM",
 "DCa_cadyn", "um2/ms",
 "vol_cadyn", "1",
 "vmaxsr_cadyn", "mM/ms",
 "kpsr_cadyn", "mM",
 "jmaxsr_cadyn", "mM/ms",
 "kip3_cadyn", "mM",
 "kactip3_cadyn", "mM",
 "konip3_cadyn", "/mM-ms",
 "kinhip3_cadyn", "mM",
 "kcicr_cadyn", "mM",
 "ktcicr_cadyn", "mM",
 "vcicr_cadyn", "/ms",
 "Kmer_cadyn", "mM",
 "Bmer_cadyn", "mM",
 "vmcu_cadyn", "mM/ms",
 "kmcu_cadyn", "mM",
 "nmcu_cadyn", "1",
 "vncx_cadyn", "mM/ms",
 "kna_cadyn", "mM",
 "kncx_cadyn", "mM",
 "Kmmt_cadyn", "mM",
 "Bmmt_cadyn", "mM",
 "k1_cadyn", "/mM-s",
 "k2_cadyn", "/s",
 "k3_cadyn", "/s",
 "k4_cadyn", "/mM-s",
 "pump0_cadyn", "mol/cm2",
 "pump_cadyn", "mol/cm2",
 "pumpca_cadyn", "mol/cm2",
 "ca_cadyn", "mM",
 "caer_cadyn", "mM",
 "camt_cadyn", "mM",
 "hc_cadyn", "1",
 "ho_cadyn", "1",
 "Ln_cadyn", "mM/ms",
 "ica_pmp_cadyn", "mA/cm2",
 "jmcu_cadyn", "mM/ms",
 "jmncx_cadyn", "mM /ms",
 "jmito_cadyn", "mM /ms",
 "jer_cadyn", "mM /ms",
 "jip3_cadyn", "mM /ms",
 "jserca_cadyn", "mM /ms",
 "jcicr_cadyn", "mM /ms",
 0,0
};
 static double Ln0 = 0;
 static double caip3ri0 = 0;
 static double camt0 = 0;
 static double caer0 = 0;
 static double ca0 = 0;
 static double delta_t = 0.01;
 static double ho0 = 0;
 static double hc0 = 0;
 static double pumpca0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "L_cadyn", &L_cadyn,
 "cai0_cadyn", &cai0_cadyn,
 "cao0_cadyn", &cao0_cadyn,
 "caeri0_cadyn", &caeri0_cadyn,
 "camti0_cadyn", &camti0_cadyn,
 "DCa_cadyn", &DCa_cadyn,
 "bbr_cadyn", &bbr_cadyn,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 "vol_cadyn", vol_cadyn, 12,
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[14]._i
 static void _ode_synonym(int, double**, Datum**);
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"cadyn",
 "vmaxsr_cadyn",
 "kpsr_cadyn",
 "jmaxsr_cadyn",
 "kip3_cadyn",
 "kactip3_cadyn",
 "konip3_cadyn",
 "kinhip3_cadyn",
 "kcicr_cadyn",
 "ktcicr_cadyn",
 "vcicr_cadyn",
 "Kmer_cadyn",
 "Bmer_cadyn",
 "vmcu_cadyn",
 "kmcu_cadyn",
 "nmcu_cadyn",
 "vncx_cadyn",
 "kna_cadyn",
 "kncx_cadyn",
 "Kmmt_cadyn",
 "Bmmt_cadyn",
 "k1_cadyn",
 "k2_cadyn",
 "k3_cadyn",
 "k4_cadyn",
 "pump0_cadyn",
 0,
 "ica_pmp_cadyn",
 "jmcu_cadyn[12]",
 "jmncx_cadyn[12]",
 "jmito_cadyn[12]",
 "jer_cadyn[12]",
 "jip3_cadyn[12]",
 "jserca_cadyn[12]",
 "jcicr_cadyn[12]",
 "fmer_cadyn[12]",
 "fmmt_cadyn[12]",
 0,
 "pump_cadyn",
 "pumpca_cadyn",
 "ca_cadyn[12]",
 "caer_cadyn[12]",
 "camt_cadyn[12]",
 "hc_cadyn[12]",
 "ho_cadyn[12]",
 "Ln_cadyn[12]",
 0,
 0};
 static Symbol* _morphology_sym;
 static Symbol* _ca_sym;
 static Symbol* _na_sym;
 static Symbol* _caer_sym;
 static Symbol* _camt_sym;
 static Symbol* _caip3r_sym;
 static Symbol* _ip3_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 299, _prop);
 	/*initialize range parameters*/
 	vmaxsr = 0.00027;
 	kpsr = 3.75e-06;
 	jmaxsr = 3.5e-06;
 	kip3 = 0.0008;
 	kactip3 = 0.0003;
 	konip3 = 2.7;
 	kinhip3 = 0.0002;
 	kcicr = 0.00198;
 	ktcicr = 0.0006;
 	vcicr = 5e-07;
 	Kmer = 0.5;
 	Bmer = 10;
 	vmcu = 1.4468e-06;
 	kmcu = 0.000606;
 	nmcu = 2.3;
 	vncx = 6e-05;
 	kna = 8;
 	kncx = 0.035;
 	Kmmt = 1e-05;
 	Bmmt = 0.065;
 	k1 = 3.74e+07;
 	k2 = 250000;
 	k3 = 500;
 	k4 = 5;
 	pump0 = 1.3725e-13;
 	_prop->param = _p;
 	_prop->param_size = 299;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 15, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_morphology_sym);
 	_ppvar[13]._pval = &prop_ion->param[0]; /* diam */
 prop_ion = need_memb(_ca_sym);
 nrn_check_conc_write(_prop, prop_ion, 1);
 nrn_promote(prop_ion, 3, 0);
 	_ppvar[0]._pval = &prop_ion->param[3]; /* ica */
 	_ppvar[1]._pval = &prop_ion->param[1]; /* cai */
 	_ppvar[2]._pval = &prop_ion->param[2]; /* cao */
 	_ppvar[3]._pval = &prop_ion->param[4]; /* _ion_dicadv */
 	_ppvar[4]._pvoid = (void*)(&(prop_ion->dparam[0]._i)); /* iontype for ca */
 prop_ion = need_memb(_na_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[5]._pval = &prop_ion->param[1]; /* nai */
 prop_ion = need_memb(_caer_sym);
 nrn_check_conc_write(_prop, prop_ion, 1);
 nrn_promote(prop_ion, 3, 0);
 	_ppvar[6]._pval = &prop_ion->param[1]; /* caeri */
 	_ppvar[7]._pvoid = (void*)(&(prop_ion->dparam[0]._i)); /* iontype for caer */
 prop_ion = need_memb(_camt_sym);
 nrn_check_conc_write(_prop, prop_ion, 1);
 nrn_promote(prop_ion, 3, 0);
 	_ppvar[8]._pval = &prop_ion->param[1]; /* camti */
 	_ppvar[9]._pvoid = (void*)(&(prop_ion->dparam[0]._i)); /* iontype for camt */
 prop_ion = need_memb(_caip3r_sym);
 nrn_check_conc_write(_prop, prop_ion, 1);
 nrn_promote(prop_ion, 3, 0);
 	_ppvar[10]._pval = &prop_ion->param[1]; /* caip3ri */
 	_ppvar[11]._pvoid = (void*)(&(prop_ion->dparam[0]._i)); /* iontype for caip3r */
 prop_ion = need_memb(_ip3_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[12]._pval = &prop_ion->param[1]; /* ip3i */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 "pump_cadyn", 1e-16,
 "pumpca_cadyn", 1e-16,
 "ca_cadyn", 1e-08,
 0,0
};
 static void _thread_mem_init(Datum*);
 static void _thread_cleanup(Datum*);
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _cadyn_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("ca", 2.0);
 	ion_reg("na", 1.0);
 	ion_reg("caer", 2.0);
 	ion_reg("camt", 2.0);
 	ion_reg("caip3r", 2.0);
 	ion_reg("ip3", 1.0);
 	_morphology_sym = hoc_lookup("morphology");
 	_ca_sym = hoc_lookup("ca_ion");
 	_na_sym = hoc_lookup("na_ion");
 	_caer_sym = hoc_lookup("caer_ion");
 	_camt_sym = hoc_lookup("camt_ion");
 	_caip3r_sym = hoc_lookup("caip3r_ion");
 	_ip3_sym = hoc_lookup("ip3_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 5);
  _extcall_thread = (Datum*)ecalloc(4, sizeof(Datum));
  _thread_mem_init(_extcall_thread);
  _thread1data_inuse = 0;
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 1, _thread_mem_init);
     _nrn_thread_reg(_mechtype, 0, _thread_cleanup);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 299, 15);
  hoc_register_dparam_semantics(_mechtype, 0, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "#ca_ion");
  hoc_register_dparam_semantics(_mechtype, 5, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 6, "caer_ion");
  hoc_register_dparam_semantics(_mechtype, 7, "#caer_ion");
  hoc_register_dparam_semantics(_mechtype, 8, "camt_ion");
  hoc_register_dparam_semantics(_mechtype, 9, "#camt_ion");
  hoc_register_dparam_semantics(_mechtype, 10, "caip3r_ion");
  hoc_register_dparam_semantics(_mechtype, 11, "#caip3r_ion");
  hoc_register_dparam_semantics(_mechtype, 12, "ip3_ion");
  hoc_register_dparam_semantics(_mechtype, 14, "cvodeieq");
  hoc_register_dparam_semantics(_mechtype, 13, "diam");
 	nrn_writes_conc(_mechtype, 0);
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_synonym(_mechtype, _ode_synonym);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 cadyn /home/mizzou/LUT_code/single_cell/components/mechanisms/x86_64/cadyn.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double F = 96485.3;
 static double PI = 3.14159;
 static double R = 8.3145;
 static double volo = 1e10;
 /*Top LOCAL _zfactors_done */
 /*Top LOCAL _zfrat [ 12 ] */
 /*Top LOCAL _zdsq , _zdsqvol , _zdsqvolmt , _zdsqvoler */
static int _reset;
static char *modelname = "Calcium Difffusion, Buffering, ER Mechs- SERCA, IP3R & CICR(RYR), and Mitochondrial Influx(MCU) and Eflux(MNCX) for bladder small DRG neuron soma model";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int factors(_threadargsproto_);
static int parms(_threadargsproto_);
 extern double *_nrn_thread_getelm();
 
#define _MATELM1(_row,_col) *(_nrn_thread_getelm(_so, _row + 1, _col + 1))
 
#define _RHS1(_arg) _rhs[_arg+1]
  
#define _linmat1  0
 static int _spth1 = 1;
 static int _cvspth1 = 0;
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[63], _dlist1[63]; static double *_temp1;
 static int state();
 
static int  factors ( _threadargsproto_ ) {
   double _lr , _ldr2 ;
 _lr = 1.0 / 2.0 ;
   _ldr2 = _lr / ( 12.0 - 1.0 ) / 2.0 ;
   vol [ 0 ] = 0.0 ;
   _zfrat [ 0 ] = 2.0 * _lr ;
   {int  _li ;for ( _li = 0 ; _li <= 12 - 2 ; _li ++ ) {
     vol [ _li ] = vol [ _li ] + PI * ( _lr - _ldr2 / 2.0 ) * 2.0 * _ldr2 ;
     _lr = _lr - _ldr2 ;
     _zfrat [ _li + 1 ] = 2.0 * PI * _lr / ( 2.0 * _ldr2 ) ;
     _lr = _lr - _ldr2 ;
     vol [ _li + 1 ] = PI * ( _lr + _ldr2 / 2.0 ) * 2.0 * _ldr2 ;
     } }
    return 0; }
 
static void _hoc_factors(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 factors ( _p, _ppvar, _thread, _nt );
 hoc_retpushx(_r);
}
 
static int state (void* _so, double* _rhs, double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt)
 {int _reset=0;
 {
   double b_flux, f_flux, _term; int _i;
 {int _i; double _dt1 = 1.0/dt;
for(_i=0;_i<63;_i++){
  	_RHS1(_i) = -_dt1*(_p[_slist1[_i]] - _p[_dlist1[_i]]);
	_MATELM1(_i, _i) = _dt1;
      
}  
_RHS1(0) *= ( ( 1.0 + bbr ) * diam * diam * vol [ 0 ] * 0.81) ;
_MATELM1(0, 0) *= ( ( 1.0 + bbr ) * diam * diam * vol [ 0 ] * 0.81); 
_RHS1(61) *= ( ( 1e10 ) * parea) ;
_MATELM1(61, 61) *= ( ( 1e10 ) * parea); 
_RHS1(62) *= ( ( 1e10 ) * parea) ;
_MATELM1(62, 62) *= ( ( 1e10 ) * parea);  
for (_i=0; _i < 12; _i++) {
  	_RHS1(_i + 1) *= ( ( 1.0 / fmmt [ ((int) _i ) ] ) * diam * diam * vol [ ((int) _i ) ] * 0.07) ;
_MATELM1(_i + 1, _i + 1) *= ( ( 1.0 / fmmt [ ((int) _i ) ] ) * diam * diam * vol [ ((int) _i ) ] * 0.07);  } 
for (_i=0; _i < 12; _i++) {
  	_RHS1(_i + 13) *= ( ( 1.0 / fmer [ ((int) _i ) ] ) * diam * diam * vol [ ((int) _i ) ] * 0.12) ;
_MATELM1(_i + 13, _i + 13) *= ( ( 1.0 / fmer [ ((int) _i ) ] ) * diam * diam * vol [ ((int) _i ) ] * 0.12);  } 
for (_i=0; _i < 12; _i++) {
  	_RHS1(_i + 25) *= ( ( 1.0 + bbr ) * diam * diam * vol [ ((int) _i ) ] * 0.81) ;
_MATELM1(_i + 25, _i + 25) *= ( ( 1.0 + bbr ) * diam * diam * vol [ ((int) _i ) ] * 0.81);  } }
 /* COMPARTMENT _lii , ( 1.0 + bbr ) * diam * diam * vol [ ((int) _i ) ] * 0.81 {
     ca }
   */
 /* COMPARTMENT ( 1.0 + bbr ) * diam * diam * vol [ 0 ] * 0.81 {
     caip3ri }
   */
 /* COMPARTMENT _ljj , ( 1.0 / fmer [ ((int) _i ) ] ) * diam * diam * vol [ ((int) _i ) ] * 0.12 {
     caer }
   */
 /* COMPARTMENT _lkk , ( 1.0 / fmmt [ ((int) _i ) ] ) * diam * diam * vol [ ((int) _i ) ] * 0.07 {
     camt }
   */
 /* COMPARTMENT ( 1e10 ) * parea {
     pump pumpca }
   */
 /* COMPARTMENT volo {
     }
   */
 /* ~ ca [ 0 ] < < ( - ( ica - last_ica_pmp ) * PI * diam * ( 1e4 ) * _zfrat [ 0 ] / ( 2.0 * F ) )*/
 f_flux = b_flux = 0.;
 _RHS1( 25 +  0) += (b_flux =   ( - ( ica - last_ica_pmp ) * PI * diam * ( 1e4 ) * _zfrat [ 0 ] / ( 2.0 * F ) ) );
 /*FLUX*/
  /* ~ ca [ 0 ] + pump <-> pumpca ( c1 , c2 )*/
 f_flux =  c1 * pump * ca [ 0] ;
 b_flux =  c2 * pumpca ;
 _RHS1( 62) -= (f_flux - b_flux);
 _RHS1( 25 +  0) -= (f_flux - b_flux);
 _RHS1( 61) += (f_flux - b_flux);
 
 _term =  c1 * ca [ 0] ;
 _MATELM1( 62 ,62)  += _term;
 _MATELM1( 25 +  0 ,62)  += _term;
 _MATELM1( 61 ,62)  -= _term;
 _term =  c1 * pump ;
 _MATELM1( 62 ,25 +  0)  += _term;
 _MATELM1( 25 +  0 ,25 +  0)  += _term;
 _MATELM1( 61 ,25 +  0)  -= _term;
 _term =  c2 ;
 _MATELM1( 62 ,61)  -= _term;
 _MATELM1( 25 +  0 ,61)  -= _term;
 _MATELM1( 61 ,61)  += _term;
 /*REACTION*/
  /* ~ pumpca <-> pump + cao ( c3 , c4 )*/
 f_flux =  c3 * pumpca ;
 b_flux =  c4 * cao * pump ;
 _RHS1( 61) -= (f_flux - b_flux);
 _RHS1( 62) += (f_flux - b_flux);
 
 _term =  c3 ;
 _MATELM1( 61 ,61)  += _term;
 _MATELM1( 62 ,61)  -= _term;
 _term =  c4 * cao ;
 _MATELM1( 61 ,62)  -= _term;
 _MATELM1( 62 ,62)  += _term;
 /*REACTION*/
  ica_pmp = ( 1e-4 ) * 2.0 * F * ( f_flux - b_flux ) / parea ;
   {int  _li ;for ( _li = 0 ; _li <= 12 - 2 ; _li ++ ) {
     /* ~ ca [ _li ] <-> ca [ _li + 1 ] ( DCa * _zfrat [ _li + 1 ] , DCa * _zfrat [ _li + 1 ] )*/
 f_flux =  DCa * _zfrat [ _li + 1 ] * ca [ _li] ;
 b_flux =  DCa * _zfrat [ _li + 1 ] * ca [ _li + 1] ;
 _RHS1( 25 +  _li) -= (f_flux - b_flux);
 _RHS1( 25 +  _li + 1) += (f_flux - b_flux);
 
 _term =  DCa * _zfrat [ _li + 1 ] ;
 _MATELM1( 25 +  _li ,25 +  _li)  += _term;
 _MATELM1( 25 +  _li + 1 ,25 +  _li)  -= _term;
 _term =  DCa * _zfrat [ _li + 1 ] ;
 _MATELM1( 25 +  _li ,25 +  _li + 1)  -= _term;
 _MATELM1( 25 +  _li + 1 ,25 +  _li + 1)  += _term;
 /*REACTION*/
  } }
   _zdsq = diam * diam ;
   {int  _li ;for ( _li = 0 ; _li <= 12 - 1 ; _li ++ ) {
     _zdsqvol = _zdsq * vol [ _li ] * 0.81 ;
     _zdsqvoler = _zdsq * vol [ _li ] * 0.12 ;
     jserca [ _li ] = ( ( - vmaxsr * pow( ca [ _li ] , 2.0 ) / ( pow( ca [ _li ] , 2.0 ) + pow( kpsr , 2.0 ) ) ) ) ;
     /* ~ ca [ _li ] < < ( _zdsqvol * jserca [ _li ] )*/
 f_flux = b_flux = 0.;
 _RHS1( 25 +  _li) += (b_flux =   ( _zdsqvol * jserca [ _li ] ) );
 /*FLUX*/
  /* ~ caer [ _li ] < < ( - _zdsqvoler * jserca [ _li ] )*/
 f_flux = b_flux = 0.;
 _RHS1( 13 +  _li) += (b_flux =   ( - _zdsqvoler * jserca [ _li ] ) );
 /*FLUX*/
  /* ~ hc [ _li ] <-> ho [ _li ] ( konip3 * kinhip3 , konip3 * ca [ _li ] )*/
 f_flux =  konip3 * kinhip3 * hc [ _li] ;
 b_flux =  konip3 * ca [ _li ] * ho [ _li] ;
 _RHS1( 49 +  _li) -= (f_flux - b_flux);
 _RHS1( 37 +  _li) += (f_flux - b_flux);
 
 _term =  konip3 * kinhip3 ;
 _MATELM1( 49 +  _li ,49 +  _li)  += _term;
 _MATELM1( 37 +  _li ,49 +  _li)  -= _term;
 _term =  konip3 * ca [ _li ] ;
 _MATELM1( 49 +  _li ,37 +  _li)  -= _term;
 _MATELM1( 37 +  _li ,37 +  _li)  += _term;
 /*REACTION*/
  jip3 [ _li ] = ( jmaxsr * ( 1.0 - ( ca [ _li ] / caer [ _li ] ) ) * pow( ( ( ip3i / ( ip3i + kip3 ) ) * ( ca [ _li ] / ( ca [ _li ] + kactip3 ) ) * ho [ _li ] ) , 3.0 ) ) ;
     /* ~ ca [ _li ] < < ( _zdsqvol * jip3 [ _li ] )*/
 f_flux = b_flux = 0.;
 _RHS1( 25 +  _li) += (b_flux =   ( _zdsqvol * jip3 [ _li ] ) );
 /*FLUX*/
  /* ~ caer [ _li ] < < ( - _zdsqvoler * jip3 [ _li ] )*/
 f_flux = b_flux = 0.;
 _RHS1( 13 +  _li) += (b_flux =   ( - _zdsqvoler * jip3 [ _li ] ) );
 /*FLUX*/
  if ( ((double) _li )  == 0.0 ) {
       /* ~ caip3ri < < ( _zdsqvol * jip3 [ 0 ] )*/
 f_flux = b_flux = 0.;
 _RHS1( 0) += (b_flux =   ( _zdsqvol * jip3 [ 0 ] ) );
 /*FLUX*/
  }
     if ( ca [ _li ] > ktcicr ) {
       jcicr [ _li ] = ( vcicr * ( ca [ _li ] / ( kcicr + ca [ _li ] ) ) * ( caer [ _li ] - ca [ _li ] ) ) ;
       /* ~ ca [ _li ] < < ( _zdsqvol * jcicr [ _li ] )*/
 f_flux = b_flux = 0.;
 _RHS1( 25 +  _li) += (b_flux =   ( _zdsqvol * jcicr [ _li ] ) );
 /*FLUX*/
  /* ~ caer [ _li ] < < ( - _zdsqvoler * jcicr [ _li ] )*/
 f_flux = b_flux = 0.;
 _RHS1( 13 +  _li) += (b_flux =   ( - _zdsqvoler * jcicr [ _li ] ) );
 /*FLUX*/
  }
     else {
       jcicr [ _li ] = 0.0 ;
       /* ~ ca [ _li ] < < ( _zdsqvol * jcicr [ _li ] )*/
 f_flux = b_flux = 0.;
 _RHS1( 25 +  _li) += (b_flux =   ( _zdsqvol * jcicr [ _li ] ) );
 /*FLUX*/
  /* ~ caer [ _li ] < < ( - _zdsqvoler * jcicr [ _li ] )*/
 f_flux = b_flux = 0.;
 _RHS1( 13 +  _li) += (b_flux =   ( - _zdsqvoler * jcicr [ _li ] ) );
 /*FLUX*/
  }
     /* ~ ca [ _li ] < < ( Ln [ _li ] * ( 1.0 - ca [ _li ] / caeri0 ) * _zdsqvol )*/
 f_flux = b_flux = 0.;
 _RHS1( 25 +  _li) += (b_flux =   ( Ln [ _li ] * ( 1.0 - ca [ _li ] / caeri0 ) * _zdsqvol ) );
 /*FLUX*/
  /* ~ caer [ _li ] < < ( - Ln [ _li ] * ( 1.0 - ca [ _li ] / caeri0 ) * _zdsqvoler )*/
 f_flux = b_flux = 0.;
 _RHS1( 13 +  _li) += (b_flux =   ( - Ln [ _li ] * ( 1.0 - ca [ _li ] / caeri0 ) * _zdsqvoler ) );
 /*FLUX*/
  jer [ _li ] = jserca [ _li ] + jip3 [ _li ] + jcicr [ _li ] ;
     fmer [ _li ] = 1.0 / ( 1.0 + ( Kmer * Bmer ) / pow( ( Kmer + caer [ _li ] ) , 2.0 ) ) ;
     _zdsqvolmt = _zdsq * vol [ _li ] * 0.07 ;
      jmcu [ _li ] = ( ( - vmcu * pow( ca [ _li ] , nmcu ) / ( pow( ca [ _li ] , nmcu ) + pow( kmcu , nmcu ) ) ) * 1.0 / ( camt [ _li ] * 1e3 ) ) ;
      /* ~ ca [ _li ] < < ( _zdsqvol * jmcu [ _li ] )*/
 f_flux = b_flux = 0.;
 _RHS1( 25 +  _li) += (b_flux =   ( _zdsqvol * jmcu [ _li ] ) );
 /*FLUX*/
  jmncx [ _li ] = vncx * ( pow( nai , 3.0 ) / ( pow( kna , 3.0 ) + pow( nai , 3.0 ) ) ) * ( camt [ _li ] / ( kncx + camt [ _li ] ) ) ;
     /* ~ ca [ _li ] < < ( jmncx [ _li ] * _zdsqvol )*/
 f_flux = b_flux = 0.;
 _RHS1( 25 +  _li) += (b_flux =   ( jmncx [ _li ] * _zdsqvol ) );
 /*FLUX*/
  jmito [ _li ] = jmcu [ _li ] + jmncx [ _li ] ;
     /* ~ camt [ _li ] < < ( - ( jmncx [ _li ] + jmcu [ _li ] ) * _zdsqvolmt )*/
 f_flux = b_flux = 0.;
 _RHS1( 1 +  _li) += (b_flux =   ( - ( jmncx [ _li ] + jmcu [ _li ] ) * _zdsqvolmt ) );
 /*FLUX*/
  fmmt [ _li ] = 1.0 / ( 1.0 + ( Kmmt * Bmmt ) / pow( ( Kmmt + camt [ _li ] ) , 2.0 ) ) ;
     } }
   cai = ca [ 0 ] ;
   caeri = caer [ 0 ] ;
   camti = camt [ 0 ] ;
     } return _reset;
 }
 
static int  parms ( _threadargsproto_ ) {
   parea = 2.0 * PI * ( diam / 2.0 ) ;
   c1 = ( 1e7 ) * parea * k1 ;
   c2 = ( 1e7 ) * parea * k2 ;
   c3 = ( 1e7 ) * parea * k3 ;
   c4 = ( 1e7 ) * parea * k4 ;
    return 0; }
 
static void _hoc_parms(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 parms ( _p, _ppvar, _thread, _nt );
 hoc_retpushx(_r);
}
 
/*CVODE ode begin*/
 static int _ode_spec1(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset=0;{
 double b_flux, f_flux, _term; int _i;
 {int _i; for(_i=0;_i<63;_i++) _p[_dlist1[_i]] = 0.0;}
 /* COMPARTMENT _lii , ( 1.0 + bbr ) * diam * diam * vol [ ((int) _i ) ] * 0.81 {
   ca }
 */
 /* COMPARTMENT ( 1.0 + bbr ) * diam * diam * vol [ 0 ] * 0.81 {
   caip3ri }
 */
 /* COMPARTMENT _ljj , ( 1.0 / fmer [ ((int) _i ) ] ) * diam * diam * vol [ ((int) _i ) ] * 0.12 {
   caer }
 */
 /* COMPARTMENT _lkk , ( 1.0 / fmmt [ ((int) _i ) ] ) * diam * diam * vol [ ((int) _i ) ] * 0.07 {
   camt }
 */
 /* COMPARTMENT ( 1e10 ) * parea {
   pump pumpca }
 */
 /* COMPARTMENT volo {
   }
 */
 /* ~ ca [ 0 ] < < ( - ( ica - last_ica_pmp ) * PI * diam * ( 1e4 ) * _zfrat [ 0 ] / ( 2.0 * F ) )*/
 f_flux = b_flux = 0.;
 Dca [ 0] += (b_flux =   ( - ( ica - last_ica_pmp ) * PI * diam * ( 1e4 ) * _zfrat [ 0 ] / ( 2.0 * F ) ) );
 /*FLUX*/
  /* ~ ca [ 0 ] + pump <-> pumpca ( c1 , c2 )*/
 f_flux =  c1 * pump * ca [ 0] ;
 b_flux =  c2 * pumpca ;
 Dpump -= (f_flux - b_flux);
 Dca [ 0] -= (f_flux - b_flux);
 Dpumpca += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ pumpca <-> pump + cao ( c3 , c4 )*/
 f_flux =  c3 * pumpca ;
 b_flux =  c4 * cao * pump ;
 Dpumpca -= (f_flux - b_flux);
 Dpump += (f_flux - b_flux);
 
 /*REACTION*/
  ica_pmp = ( 1e-4 ) * 2.0 * F * ( f_flux - b_flux ) / parea ;
 {int  _li ;for ( _li = 0 ; _li <= 12 - 2 ; _li ++ ) {
   /* ~ ca [ _li ] <-> ca [ _li + 1 ] ( DCa * _zfrat [ _li + 1 ] , DCa * _zfrat [ _li + 1 ] )*/
 f_flux =  DCa * _zfrat [ _li + 1 ] * ca [ _li] ;
 b_flux =  DCa * _zfrat [ _li + 1 ] * ca [ _li + 1] ;
 Dca [ _li] -= (f_flux - b_flux);
 Dca [ _li + 1] += (f_flux - b_flux);
 
 /*REACTION*/
  } }
 _zdsq = diam * diam ;
 {int  _li ;for ( _li = 0 ; _li <= 12 - 1 ; _li ++ ) {
   _zdsqvol = _zdsq * vol [ _li ] * 0.81 ;
   _zdsqvoler = _zdsq * vol [ _li ] * 0.12 ;
   jserca [ _li ] = ( ( - vmaxsr * pow( ca [ _li ] , 2.0 ) / ( pow( ca [ _li ] , 2.0 ) + pow( kpsr , 2.0 ) ) ) ) ;
   /* ~ ca [ _li ] < < ( _zdsqvol * jserca [ _li ] )*/
 f_flux = b_flux = 0.;
 Dca [ _li] += (b_flux =   ( _zdsqvol * jserca [ _li ] ) );
 /*FLUX*/
  /* ~ caer [ _li ] < < ( - _zdsqvoler * jserca [ _li ] )*/
 f_flux = b_flux = 0.;
 Dcaer [ _li] += (b_flux =   ( - _zdsqvoler * jserca [ _li ] ) );
 /*FLUX*/
  /* ~ hc [ _li ] <-> ho [ _li ] ( konip3 * kinhip3 , konip3 * ca [ _li ] )*/
 f_flux =  konip3 * kinhip3 * hc [ _li] ;
 b_flux =  konip3 * ca [ _li ] * ho [ _li] ;
 Dhc [ _li] -= (f_flux - b_flux);
 Dho [ _li] += (f_flux - b_flux);
 
 /*REACTION*/
  jip3 [ _li ] = ( jmaxsr * ( 1.0 - ( ca [ _li ] / caer [ _li ] ) ) * pow( ( ( ip3i / ( ip3i + kip3 ) ) * ( ca [ _li ] / ( ca [ _li ] + kactip3 ) ) * ho [ _li ] ) , 3.0 ) ) ;
   /* ~ ca [ _li ] < < ( _zdsqvol * jip3 [ _li ] )*/
 f_flux = b_flux = 0.;
 Dca [ _li] += (b_flux =   ( _zdsqvol * jip3 [ _li ] ) );
 /*FLUX*/
  /* ~ caer [ _li ] < < ( - _zdsqvoler * jip3 [ _li ] )*/
 f_flux = b_flux = 0.;
 Dcaer [ _li] += (b_flux =   ( - _zdsqvoler * jip3 [ _li ] ) );
 /*FLUX*/
  if ( ((double) _li )  == 0.0 ) {
     /* ~ caip3ri < < ( _zdsqvol * jip3 [ 0 ] )*/
 f_flux = b_flux = 0.;
 Dcaip3ri += (b_flux =   ( _zdsqvol * jip3 [ 0 ] ) );
 /*FLUX*/
  }
   if ( ca [ _li ] > ktcicr ) {
     jcicr [ _li ] = ( vcicr * ( ca [ _li ] / ( kcicr + ca [ _li ] ) ) * ( caer [ _li ] - ca [ _li ] ) ) ;
     /* ~ ca [ _li ] < < ( _zdsqvol * jcicr [ _li ] )*/
 f_flux = b_flux = 0.;
 Dca [ _li] += (b_flux =   ( _zdsqvol * jcicr [ _li ] ) );
 /*FLUX*/
  /* ~ caer [ _li ] < < ( - _zdsqvoler * jcicr [ _li ] )*/
 f_flux = b_flux = 0.;
 Dcaer [ _li] += (b_flux =   ( - _zdsqvoler * jcicr [ _li ] ) );
 /*FLUX*/
  }
   else {
     jcicr [ _li ] = 0.0 ;
     /* ~ ca [ _li ] < < ( _zdsqvol * jcicr [ _li ] )*/
 f_flux = b_flux = 0.;
 Dca [ _li] += (b_flux =   ( _zdsqvol * jcicr [ _li ] ) );
 /*FLUX*/
  /* ~ caer [ _li ] < < ( - _zdsqvoler * jcicr [ _li ] )*/
 f_flux = b_flux = 0.;
 Dcaer [ _li] += (b_flux =   ( - _zdsqvoler * jcicr [ _li ] ) );
 /*FLUX*/
  }
   /* ~ ca [ _li ] < < ( Ln [ _li ] * ( 1.0 - ca [ _li ] / caeri0 ) * _zdsqvol )*/
 f_flux = b_flux = 0.;
 Dca [ _li] += (b_flux =   ( Ln [ _li ] * ( 1.0 - ca [ _li ] / caeri0 ) * _zdsqvol ) );
 /*FLUX*/
  /* ~ caer [ _li ] < < ( - Ln [ _li ] * ( 1.0 - ca [ _li ] / caeri0 ) * _zdsqvoler )*/
 f_flux = b_flux = 0.;
 Dcaer [ _li] += (b_flux =   ( - Ln [ _li ] * ( 1.0 - ca [ _li ] / caeri0 ) * _zdsqvoler ) );
 /*FLUX*/
  jer [ _li ] = jserca [ _li ] + jip3 [ _li ] + jcicr [ _li ] ;
   fmer [ _li ] = 1.0 / ( 1.0 + ( Kmer * Bmer ) / pow( ( Kmer + caer [ _li ] ) , 2.0 ) ) ;
   _zdsqvolmt = _zdsq * vol [ _li ] * 0.07 ;
    jmcu [ _li ] = ( ( - vmcu * pow( ca [ _li ] , nmcu ) / ( pow( ca [ _li ] , nmcu ) + pow( kmcu , nmcu ) ) ) * 1.0 / ( camt [ _li ] * 1e3 ) ) ;
    /* ~ ca [ _li ] < < ( _zdsqvol * jmcu [ _li ] )*/
 f_flux = b_flux = 0.;
 Dca [ _li] += (b_flux =   ( _zdsqvol * jmcu [ _li ] ) );
 /*FLUX*/
  jmncx [ _li ] = vncx * ( pow( nai , 3.0 ) / ( pow( kna , 3.0 ) + pow( nai , 3.0 ) ) ) * ( camt [ _li ] / ( kncx + camt [ _li ] ) ) ;
   /* ~ ca [ _li ] < < ( jmncx [ _li ] * _zdsqvol )*/
 f_flux = b_flux = 0.;
 Dca [ _li] += (b_flux =   ( jmncx [ _li ] * _zdsqvol ) );
 /*FLUX*/
  jmito [ _li ] = jmcu [ _li ] + jmncx [ _li ] ;
   /* ~ camt [ _li ] < < ( - ( jmncx [ _li ] + jmcu [ _li ] ) * _zdsqvolmt )*/
 f_flux = b_flux = 0.;
 Dcamt [ _li] += (b_flux =   ( - ( jmncx [ _li ] + jmcu [ _li ] ) * _zdsqvolmt ) );
 /*FLUX*/
  fmmt [ _li ] = 1.0 / ( 1.0 + ( Kmmt * Bmmt ) / pow( ( Kmmt + camt [ _li ] ) , 2.0 ) ) ;
   } }
 cai = ca [ 0 ] ;
 caeri = caer [ 0 ] ;
 camti = camt [ 0 ] ;
 _p[_dlist1[0]] /= ( ( 1.0 + bbr ) * diam * diam * vol [ 0 ] * 0.81);
 for (_i=0; _i < 12; _i++) { _p[_dlist1[_i + 1]] /= ( ( 1.0 / fmmt [ ((int) _i ) ] ) * diam * diam * vol [ ((int) _i ) ] * 0.07);}
 for (_i=0; _i < 12; _i++) { _p[_dlist1[_i + 13]] /= ( ( 1.0 / fmer [ ((int) _i ) ] ) * diam * diam * vol [ ((int) _i ) ] * 0.12);}
 for (_i=0; _i < 12; _i++) { _p[_dlist1[_i + 25]] /= ( ( 1.0 + bbr ) * diam * diam * vol [ ((int) _i ) ] * 0.81);}
 _p[_dlist1[61]] /= ( ( 1e10 ) * parea);
 _p[_dlist1[62]] /= ( ( 1e10 ) * parea);
   } return _reset;
 }
 
/*CVODE matsol*/
 static int _ode_matsol1(void* _so, double* _rhs, double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset=0;{
 double b_flux, f_flux, _term; int _i;
   b_flux = f_flux = 0.;
 {int _i; double _dt1 = 1.0/dt;
for(_i=0;_i<63;_i++){
  	_RHS1(_i) = _dt1*(_p[_dlist1[_i]]);
	_MATELM1(_i, _i) = _dt1;
      
}  
_RHS1(0) *= ( ( 1.0 + bbr ) * diam * diam * vol [ 0 ] * 0.81) ;
_MATELM1(0, 0) *= ( ( 1.0 + bbr ) * diam * diam * vol [ 0 ] * 0.81); 
_RHS1(61) *= ( ( 1e10 ) * parea) ;
_MATELM1(61, 61) *= ( ( 1e10 ) * parea); 
_RHS1(62) *= ( ( 1e10 ) * parea) ;
_MATELM1(62, 62) *= ( ( 1e10 ) * parea);  
for (_i=0; _i < 12; _i++) {
  	_RHS1(_i + 1) *= ( ( 1.0 / fmmt [ ((int) _i ) ] ) * diam * diam * vol [ ((int) _i ) ] * 0.07) ;
_MATELM1(_i + 1, _i + 1) *= ( ( 1.0 / fmmt [ ((int) _i ) ] ) * diam * diam * vol [ ((int) _i ) ] * 0.07);  } 
for (_i=0; _i < 12; _i++) {
  	_RHS1(_i + 13) *= ( ( 1.0 / fmer [ ((int) _i ) ] ) * diam * diam * vol [ ((int) _i ) ] * 0.12) ;
_MATELM1(_i + 13, _i + 13) *= ( ( 1.0 / fmer [ ((int) _i ) ] ) * diam * diam * vol [ ((int) _i ) ] * 0.12);  } 
for (_i=0; _i < 12; _i++) {
  	_RHS1(_i + 25) *= ( ( 1.0 + bbr ) * diam * diam * vol [ ((int) _i ) ] * 0.81) ;
_MATELM1(_i + 25, _i + 25) *= ( ( 1.0 + bbr ) * diam * diam * vol [ ((int) _i ) ] * 0.81);  } }
 /* COMPARTMENT _lii , ( 1.0 + bbr ) * diam * diam * vol [ ((int) _i ) ] * 0.81 {
 ca }
 */
 /* COMPARTMENT ( 1.0 + bbr ) * diam * diam * vol [ 0 ] * 0.81 {
 caip3ri }
 */
 /* COMPARTMENT _ljj , ( 1.0 / fmer [ ((int) _i ) ] ) * diam * diam * vol [ ((int) _i ) ] * 0.12 {
 caer }
 */
 /* COMPARTMENT _lkk , ( 1.0 / fmmt [ ((int) _i ) ] ) * diam * diam * vol [ ((int) _i ) ] * 0.07 {
 camt }
 */
 /* COMPARTMENT ( 1e10 ) * parea {
 pump pumpca }
 */
 /* COMPARTMENT volo {
 }
 */
 /* ~ ca [ 0 ] < < ( - ( ica - last_ica_pmp ) * PI * diam * ( 1e4 ) * _zfrat [ 0 ] / ( 2.0 * F ) )*/
 /*FLUX*/
  /* ~ ca [ 0 ] + pump <-> pumpca ( c1 , c2 )*/
 _term =  c1 * ca [ 0] ;
 _MATELM1( 62 ,62)  += _term;
 _MATELM1( 25 +  0 ,62)  += _term;
 _MATELM1( 61 ,62)  -= _term;
 _term =  c1 * pump ;
 _MATELM1( 62 ,25 +  0)  += _term;
 _MATELM1( 25 +  0 ,25 +  0)  += _term;
 _MATELM1( 61 ,25 +  0)  -= _term;
 _term =  c2 ;
 _MATELM1( 62 ,61)  -= _term;
 _MATELM1( 25 +  0 ,61)  -= _term;
 _MATELM1( 61 ,61)  += _term;
 /*REACTION*/
  /* ~ pumpca <-> pump + cao ( c3 , c4 )*/
 _term =  c3 ;
 _MATELM1( 61 ,61)  += _term;
 _MATELM1( 62 ,61)  -= _term;
 _term =  c4 * cao ;
 _MATELM1( 61 ,62)  -= _term;
 _MATELM1( 62 ,62)  += _term;
 {int  _li ;for ( _li = 0 ; _li <= 12 - 2 ; _li ++ ) {
 /* ~ ca [ _li ] <-> ca [ _li + 1 ] ( DCa * _zfrat [ _li + 1 ] , DCa * _zfrat [ _li + 1 ] )*/
 _term =  DCa * _zfrat [ _li + 1 ] ;
 _MATELM1( 25 +  _li ,25 +  _li)  += _term;
 _MATELM1( 25 +  _li + 1 ,25 +  _li)  -= _term;
 _term =  DCa * _zfrat [ _li + 1 ] ;
 _MATELM1( 25 +  _li ,25 +  _li + 1)  -= _term;
 _MATELM1( 25 +  _li + 1 ,25 +  _li + 1)  += _term;
 /*REACTION*/
  } }
 _zdsq = diam * diam ;
 {int  _li ;for ( _li = 0 ; _li <= 12 - 1 ; _li ++ ) {
 _zdsqvol = _zdsq * vol [ _li ] * 0.81 ;
 _zdsqvoler = _zdsq * vol [ _li ] * 0.12 ;
 jserca [ _li ] = ( ( - vmaxsr * pow( ca [ _li ] , 2.0 ) / ( pow( ca [ _li ] , 2.0 ) + pow( kpsr , 2.0 ) ) ) ) ;
 /* ~ ca [ _li ] < < ( _zdsqvol * jserca [ _li ] )*/
 /*FLUX*/
  /* ~ caer [ _li ] < < ( - _zdsqvoler * jserca [ _li ] )*/
 /*FLUX*/
  /* ~ hc [ _li ] <-> ho [ _li ] ( konip3 * kinhip3 , konip3 * ca [ _li ] )*/
 _term =  konip3 * kinhip3 ;
 _MATELM1( 49 +  _li ,49 +  _li)  += _term;
 _MATELM1( 37 +  _li ,49 +  _li)  -= _term;
 _term =  konip3 * ca [ _li ] ;
 _MATELM1( 49 +  _li ,37 +  _li)  -= _term;
 _MATELM1( 37 +  _li ,37 +  _li)  += _term;
 /*REACTION*/
  jip3 [ _li ] = ( jmaxsr * ( 1.0 - ( ca [ _li ] / caer [ _li ] ) ) * pow( ( ( ip3i / ( ip3i + kip3 ) ) * ( ca [ _li ] / ( ca [ _li ] + kactip3 ) ) * ho [ _li ] ) , 3.0 ) ) ;
 /* ~ ca [ _li ] < < ( _zdsqvol * jip3 [ _li ] )*/
 /*FLUX*/
  /* ~ caer [ _li ] < < ( - _zdsqvoler * jip3 [ _li ] )*/
 /*FLUX*/
  if ( ((double) _li )  == 0.0 ) {
   /* ~ caip3ri < < ( _zdsqvol * jip3 [ 0 ] )*/
 /*FLUX*/
  }
 if ( ca [ _li ] > ktcicr ) {
   jcicr [ _li ] = ( vcicr * ( ca [ _li ] / ( kcicr + ca [ _li ] ) ) * ( caer [ _li ] - ca [ _li ] ) ) ;
   /* ~ ca [ _li ] < < ( _zdsqvol * jcicr [ _li ] )*/
 /*FLUX*/
  /* ~ caer [ _li ] < < ( - _zdsqvoler * jcicr [ _li ] )*/
 /*FLUX*/
  }
 else {
   jcicr [ _li ] = 0.0 ;
   /* ~ ca [ _li ] < < ( _zdsqvol * jcicr [ _li ] )*/
 /*FLUX*/
  /* ~ caer [ _li ] < < ( - _zdsqvoler * jcicr [ _li ] )*/
 /*FLUX*/
  }
 /* ~ ca [ _li ] < < ( Ln [ _li ] * ( 1.0 - ca [ _li ] / caeri0 ) * _zdsqvol )*/
 /*FLUX*/
  /* ~ caer [ _li ] < < ( - Ln [ _li ] * ( 1.0 - ca [ _li ] / caeri0 ) * _zdsqvoler )*/
 /*FLUX*/
  jer [ _li ] = jserca [ _li ] + jip3 [ _li ] + jcicr [ _li ] ;
 fmer [ _li ] = 1.0 / ( 1.0 + ( Kmer * Bmer ) / pow( ( Kmer + caer [ _li ] ) , 2.0 ) ) ;
 _zdsqvolmt = _zdsq * vol [ _li ] * 0.07 ;
  jmcu [ _li ] = ( ( - vmcu * pow( ca [ _li ] , nmcu ) / ( pow( ca [ _li ] , nmcu ) + pow( kmcu , nmcu ) ) ) * 1.0 / ( camt [ _li ] * 1e3 ) ) ;
  /* ~ ca [ _li ] < < ( _zdsqvol * jmcu [ _li ] )*/
 /*FLUX*/
  jmncx [ _li ] = vncx * ( pow( nai , 3.0 ) / ( pow( kna , 3.0 ) + pow( nai , 3.0 ) ) ) * ( camt [ _li ] / ( kncx + camt [ _li ] ) ) ;
 /* ~ ca [ _li ] < < ( jmncx [ _li ] * _zdsqvol )*/
 /*FLUX*/
  jmito [ _li ] = jmcu [ _li ] + jmncx [ _li ] ;
 /* ~ camt [ _li ] < < ( - ( jmncx [ _li ] + jmcu [ _li ] ) * _zdsqvolmt )*/
 /*FLUX*/
  fmmt [ _li ] = 1.0 / ( 1.0 + ( Kmmt * Bmmt ) / pow( ( Kmmt + camt [ _li ] ) , 2.0 ) ) ;
 } }
 cai = ca [ 0 ] ;
 caeri = caer [ 0 ] ;
 camti = camt [ 0 ] ;
   } return _reset;
 }
 
/*CVODE end*/
 
static int _ode_count(int _type){ return 63;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ica = _ion_ica;
  cai = _ion_cai;
  cao = _ion_cao;
  cai = _ion_cai;
  nai = _ion_nai;
  caeri = _ion_caeri;
  caeri = _ion_caeri;
  camti = _ion_camti;
  camti = _ion_camti;
  caip3ri = _ion_caip3ri;
  caip3ri = _ion_caip3ri;
  ip3i = _ion_ip3i;
     _ode_spec1 (_p, _ppvar, _thread, _nt);
  _ion_cai = cai;
   _ion_caeri = caeri;
  _ion_camti = camti;
  _ion_caip3ri = caip3ri;
 }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 63; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 	_pv[0] = &(_ion_caip3ri);
 }
 static void _ode_synonym(int _cnt, double** _pp, Datum** _ppd) { 
	double* _p; Datum* _ppvar;
 	int _i; 
	for (_i=0; _i < _cnt; ++_i) {_p = _pp[_i]; _ppvar = _ppd[_i];
 _ion_cai =  ca [ 0 ] ;
 _ion_caeri =  caer [ 0 ] ;
 _ion_camti =  camt [ 0 ] ;
 }}
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _cvode_sparse_thread(&_thread[_cvspth1]._pvoid, 63, _dlist1, _p, _ode_matsol1, _ppvar, _thread, _nt);
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ica = _ion_ica;
  cai = _ion_cai;
  cao = _ion_cao;
  cai = _ion_cai;
  nai = _ion_nai;
  caeri = _ion_caeri;
  caeri = _ion_caeri;
  camti = _ion_camti;
  camti = _ion_camti;
  caip3ri = _ion_caip3ri;
  caip3ri = _ion_caip3ri;
  ip3i = _ion_ip3i;
 _ode_matsol_instance1(_threadargs_);
 }}
 
static void _thread_mem_init(Datum* _thread) {
   _thread[2]._pval = (double*)ecalloc(17, sizeof(double));
  if (_thread1data_inuse) {_thread[_gth]._pval = (double*)ecalloc(12, sizeof(double));
 }else{
 _thread[_gth]._pval = _thread1data; _thread1data_inuse = 1;
 }
 }
 
static void _thread_cleanup(Datum* _thread) {
   _nrn_destroy_sparseobj_thread(_thread[_cvspth1]._pvoid);
   _nrn_destroy_sparseobj_thread(_thread[_spth1]._pvoid);
   free((void*)(_thread[2]._pval));
  if (_thread[_gth]._pval == _thread1data) {
   _thread1data_inuse = 0;
  }else{
   free((void*)_thread[_gth]._pval);
  }
 }
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_ca_sym, _ppvar, 0, 3);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 1, 1);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 2, 2);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 3, 4);
   nrn_update_ion_pointer(_na_sym, _ppvar, 5, 1);
   nrn_update_ion_pointer(_caer_sym, _ppvar, 6, 1);
   nrn_update_ion_pointer(_camt_sym, _ppvar, 8, 1);
   nrn_update_ion_pointer(_caip3r_sym, _ppvar, 10, 1);
   nrn_update_ion_pointer(_ip3_sym, _ppvar, 12, 1);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
 for (_i=0; _i<12; _i++) Ln[_i] = Ln0;
 for (_i=0; _i<12; _i++) camt[_i] = camt0;
 for (_i=0; _i<12; _i++) caer[_i] = caer0;
 for (_i=0; _i<12; _i++) ca[_i] = ca0;
 for (_i=0; _i<12; _i++) ho[_i] = ho0;
 for (_i=0; _i<12; _i++) hc[_i] = hc0;
  pumpca = pumpca0;
  pump = pump0;
 {
   double _ltotal ;
 if ( _zfactors_done  == 0.0 ) {
     _zfactors_done = 1.0 ;
     factors ( _threadargs_ ) ;
     }
   cai = cai0 ;
   parms ( _threadargs_ ) ;
   parea = PI * diam ;
   pump = pump0 ;
   pumpca = cai * pump * k1 / k2 ;
   _ltotal = pumpca + pump ;
   if ( _ltotal > 1e-9 ) {
     pump = pump * ( pump / _ltotal ) ;
     pumpca = pumpca * ( pump / _ltotal ) ;
     }
   ica_pmp = 0.0 ;
   last_ica_pmp = 0.0 ;
   {int  _li ;for ( _li = 0 ; _li <= 12 - 1 ; _li ++ ) {
     ca [ _li ] = cai0 ;
     caer [ _li ] = caeri0 ;
     camt [ _li ] = camti0 ;
     caip3ri = cai0 ;
     jserca [ _li ] = 0.0 ;
     jip3 [ _li ] = 0.0 ;
     jcicr [ _li ] = 0.0 ;
     jmcu [ _li ] = 0.0 ;
     jmncx [ _li ] = 0.0 ;
     fmmt [ _li ] = 1.0 / ( 1.0 + ( Kmmt * Bmmt ) / pow( ( Kmmt + camt [ _li ] ) , 2.0 ) ) ;
     fmer [ _li ] = 1.0 / ( 1.0 + ( Kmer * Bmer ) / pow( ( Kmer + caer [ _li ] ) , 2.0 ) ) ;
     } }
   caeri = caer [ 0 ] ;
   camti = camt [ 0 ] ;
   {int  _li ;for ( _li = 0 ; _li <= 12 - 1 ; _li ++ ) {
     ho [ _li ] = kinhip3 / ( ca [ _li ] + kinhip3 ) ;
     hc [ _li ] = 1.0 - ho [ _li ] ;
     jserca [ _li ] = ( - vmaxsr * pow( ca [ _li ] , 2.0 ) / ( pow( ca [ _li ] , 2.0 ) + pow( kpsr , 2.0 ) ) ) ;
     jip3 [ _li ] = ( jmaxsr * ( 1.0 - ( ca [ _li ] / caer [ _li ] ) ) * pow( ( ( ip3i / ( ip3i + kip3 ) ) * ( ca [ _li ] / ( ca [ _li ] + kactip3 ) ) * ho [ _li ] ) , 3.0 ) ) ;
     if ( ca [ _li ] > ktcicr ) {
       jcicr [ _li ] = ( vcicr * ( ca [ _li ] / ( kcicr + ca [ _li ] ) ) * ( caer [ _li ] - ca [ _li ] ) ) ;
       }
     else {
       jcicr [ _li ] = 0.0 ;
       }
     jer [ _li ] = jserca [ _li ] + jip3 [ _li ] + jcicr [ _li ] ;
      jmcu [ _li ] = ( - vmcu * pow( ca [ _li ] , nmcu ) / ( pow( ca [ _li ] , nmcu ) + pow( kmcu , nmcu ) ) ) ;
      jmncx [ _li ] = vncx * ( pow( nai , 3.0 ) / ( pow( kna , 3.0 ) + pow( nai , 3.0 ) ) ) * ( camt [ _li ] / ( kncx + camt [ _li ] ) ) ;
     jmito [ _li ] = jmcu [ _li ] + jmncx [ _li ] ;
     Ln [ _li ] = - ( jserca [ _li ] + jip3 [ _li ] + jcicr [ _li ] ) / ( 1.0 - ( ca [ _li ] / caeri0 ) ) ;
     } }
   }
 
}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  ica = _ion_ica;
  cai = _ion_cai;
  cao = _ion_cao;
  cai = _ion_cai;
  nai = _ion_nai;
  caeri = _ion_caeri;
  caeri = _ion_caeri;
  camti = _ion_camti;
  camti = _ion_camti;
  caip3ri = _ion_caip3ri;
  caip3ri = _ion_caip3ri;
  ip3i = _ion_ip3i;
 initmodel(_p, _ppvar, _thread, _nt);
  _ion_cai = cai;
   nrn_wrote_conc(_ca_sym, (&(_ion_cai)) - 1, _style_ca);
  _ion_caeri = caeri;
  nrn_wrote_conc(_caer_sym, (&(_ion_caeri)) - 1, _style_caer);
  _ion_camti = camti;
  nrn_wrote_conc(_camt_sym, (&(_ion_camti)) - 1, _style_camt);
  _ion_caip3ri = caip3ri;
  nrn_wrote_conc(_caip3r_sym, (&(_ion_caip3ri)) - 1, _style_caip3r);
}
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   last_ica_pmp = ica_pmp ;
   ica = ica_pmp ;
   }
 _current += ica;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
  ica = _ion_ica;
  cai = _ion_cai;
  cao = _ion_cao;
  cai = _ion_cai;
  nai = _ion_nai;
  caeri = _ion_caeri;
  caeri = _ion_caeri;
  camti = _ion_camti;
  camti = _ion_camti;
  caip3ri = _ion_caip3ri;
  caip3ri = _ion_caip3ri;
  ip3i = _ion_ip3i;
if (_nt->_vcv) { _ode_spec1(_p, _ppvar, _thread, _nt); }
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dica;
  _dica = ica;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dicadv += (_dica - ica)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_cai = cai;
  _ion_ica += ica ;
  _ion_caeri = caeri;
  _ion_camti = camti;
  _ion_caip3ri = caip3ri;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}
 
}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}
 
}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
double _dtsav = dt;
if (secondorder) { dt *= 0.5; }
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
  ica = _ion_ica;
  cai = _ion_cai;
  cao = _ion_cao;
  cai = _ion_cai;
  nai = _ion_nai;
  caeri = _ion_caeri;
  caeri = _ion_caeri;
  camti = _ion_camti;
  camti = _ion_camti;
  caip3ri = _ion_caip3ri;
  caip3ri = _ion_caip3ri;
  ip3i = _ion_ip3i;
 {  sparse_thread(&_thread[_spth1]._pvoid, 63, _slist1, _dlist1, _p, &t, dt, state, _linmat1, _ppvar, _thread, _nt);
     if (secondorder) {
    int _i;
    for (_i = 0; _i < 63; ++_i) {
      _p[_slist1[_i]] += dt*_p[_dlist1[_i]];
    }}
 }  _ion_cai = cai;
   _ion_caeri = caeri;
  _ion_camti = camti;
  _ion_caip3ri = caip3ri;
}}
 dt = _dtsav;
}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(caip3ri) - _p;  _dlist1[0] = &(Dcaip3ri) - _p;
 for(_i=0;_i<12;_i++){_slist1[1+_i] = (camt + _i) - _p;  _dlist1[1+_i] = (Dcamt + _i) - _p;}
 for(_i=0;_i<12;_i++){_slist1[13+_i] = (caer + _i) - _p;  _dlist1[13+_i] = (Dcaer + _i) - _p;}
 for(_i=0;_i<12;_i++){_slist1[25+_i] = (ca + _i) - _p;  _dlist1[25+_i] = (Dca + _i) - _p;}
 for(_i=0;_i<12;_i++){_slist1[37+_i] = (ho + _i) - _p;  _dlist1[37+_i] = (Dho + _i) - _p;}
 for(_i=0;_i<12;_i++){_slist1[49+_i] = (hc + _i) - _p;  _dlist1[49+_i] = (Dhc + _i) - _p;}
 _slist1[61] = &(pumpca) - _p;  _dlist1[61] = &(Dpumpca) - _p;
 _slist1[62] = &(pump) - _p;  _dlist1[62] = &(Dpump) - _p;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/home/mizzou/LUT_code/single_cell/components/mechanisms/modfiles/cadyn.mod";
static const char* nmodl_file_text = 
  "TITLE Calcium Difffusion, Buffering, ER Mechs- SERCA, IP3R & CICR(RYR), and Mitochondrial Influx(MCU) and Eflux(MNCX) for bladder small DRG neuron soma model\n"
  "\n"
  ": Author: Darshan Mandge (darshanmandge@iitb.ac.in)\n"
  ": Computational Neurophysiology Lab\n"
  ": Indian Institute of Technology Bombay, India \n"
  "\n"
  ": For details refer: \n"
  ": A biophysically detailed computational model of bladder small DRG neuron soma \n"
  ": Darshan Mandge and Rohit Manchanda, PLOS Computational Biology (2018)\n"
  "\n"
  "NEURON{\n"
  "	SUFFIX cadyn\n"
  "	USEION ca READ ica, cai, cao WRITE cai, ica VALENCE 2 		: writing cai and ica\n"
  "	USEION na READ nai VALENCE 1								: for mncx\n"
  "	USEION caer READ caeri WRITE caeri VALENCE 2 				: caeri = internal ER calcium. WRITING caeri\n"
  "	USEION camt READ camti WRITE camti VALENCE 2 				: camti = internal mitochondrial calcium\n"
  "	USEION caip3r READ caip3ri WRITE caip3ri VALENCE 2			: caip3ri is the ca change by ca release by outermost shell's IP3Rs. It is read by CACC\n"
  "	USEION ip3 READ ip3i VALENCE 1								: ip3i = ip3 conc. \n"
  "   \n"
  "	GLOBAL DCa, cai0, caeri0, camti0							: Diffusion constant of Ca, calcium conc. in cytoplasm, ER and mitochondria \n"
  "	RANGE k1, k2, k3, k4, ica_pmp, pump0 			 		 	: Pump Parameters: RANGE as can be different for different segments and sections\n"
  "	GLOBAL bbr													: Buffer Parameters: GLOBAL as they are props. of buffer and are constant         \n"
  "	:RANGE kmdye, Bmdye, Dbufm, bbrdye							: Dye correction: fura-2 in O'mullane-2013\n"
  "	\n"
  "	RANGE vmaxsr, kpsr		                        			: SERCA parameters\n"
  "	RANGE kactip3, konip3, kinhip3,  kip3, jmaxsr			  	: IP3 Parameters  \n"
  "	RANGE ktcicr, kcicr, vcicr                           		: CICR Parameters\n"
  "    RANGE vmcu, kmcu, nmcu, vncx, kna, kncx			    	   	: Mitochondrial Parameters: Uniporter and MNCX\n"
  "	RANGE jer, jserca, jip3, jcicr								: Flux Parameters\n"
  "	RANGE jmcu, jmncx, jmito\n"
  "\n"
  "	RANGE Kmmt, Bmmt, fmmt										: Mitochondrial Buffer Paramters\n"
  "	RANGE Kmer, Bmer, fmer										: ER Buffer Paramters\n"
  "	\n"
  "    GLOBAL vol\n"
  "	THREADSAFE\n"
  "}\n"
  "DEFINE NANN  12    :IF YOU CHANGE THIS NUMBER DONT FORGET to change the NANN value in ip3dif.mod file.\n"
  "\n"
  "UNITS{\n"
  "	(mV)	= (millivolt)\n"
  "	(um)    = (micron)\n"
  "	(mM)    = (milli/liter)\n"
  "    (nM)    = (nano/liter)\n"
  "	(mA)    = (milliamp)\n"
  "	F       = (faraday) (coulombs)\n"
  "	PI      = (pi) (1)\n"
  "	R 		= (k-mole) (joule/degC)\n"
  "    (mol)   = (1) :The term mole cannot be used here because it is already defined in NEURON's units database as 6.022169*10^23\n"
  "}\n"
  "\n"
  "PARAMETER{\n"
  "	diam	(um)\n"
  "	L		(um)\n"
  ":-------------------------------------------------------------------------------------------------------------------------------------------------------------------\n"
  "	:Ca Parameters\n"
  "	cai0 	= 136e-6	(mM)\n"
  "	cao0 	= 2 		(mM)\n"
  "	caeri0 	= 0.4		(mM)\n"
  "   	camti0 	= 2e-4 		(mM)\n"
  "	DCa 	= 0.6		(um2/ms)     : 0.6 in McHugh et al., 2004\n"
  ":-------------------------------------------------------------------------------------------------------------------------------------------------------------------\n"
  "	: Buffers\n"
  "	bbr = 370 		:Zeilhofer 1996 (J. neurophysiology) endogenous buffer binding ratio = 370 (Rat DRG)\n"
  "	\n"
  "	: Dye paramters\n"
  "	: Fura-2 O'Mullane et al., 2013 Data for bladder cai measurement\n"
  "	: kmdye = 224e-6 (mM)		: Grynkiewicz et al., 1985 224 nM in Lu and Gold 2006, 2008\n"
  "	: Bmdye = 5e-3   (mM)		: O'Mullane 2013 5 uM Fura-2\n"
  "	: Dbufm = 0.1 (um2/ms)      : Blatter and Wier, 1990 \n"
  "	\n"
  "	: Indo-1 Data for cai \n"
  "	: kmdye = 250e-6 (mM)	    : 250 nM Grynkiewicz et al., 1985		\n"
  "	: Bmdye = 100e-3 (mM)     	: 100 uM   Benham et al.,1992\n"
  "	: Dbufm =  0.1 (um2/ms)     : Blatter and Wier, 1990\n"
  "	\n"
  "	: Fura-FF shutov\n"
  "	: kmdye = 5500e-6 (mM)	    : 5.5 uM Shutov et al., 2013\n"
  "	: Bmdye = 200e-3 (mM)     	: 200 uM Shutov et al., 2013\n"
  "	: Dbufm =  0.075 (um2/ms)   : Assumed similar to Fura-2 and indo-1\n"
  "	\n"
  "	: mt-pericam\n"
  "	: kmmitodye = 1.7e-3 (mM)	: 1.7 Nagai et al., 2001\n"
  "	: Bmmitodye = 15 (mM)\n"
  "	\n"
  ":-------------------------------------------------------------------------------------------------------------------------------------------------------------------\n"
  "	:SERCA\n"
  "	vmaxsr 	= 0.00027	(mM/ms)\n"
  "	kpsr 	= 3.75e-6 	(mM) : Fink et al., 2000\n"
  "	\n"
  ":------------------------------------------------------------------------------------------------------------------------------------------------------------------\n"
  "	: IP3R\n"
  "    jmaxsr 	= 3.5e-6	(mM/ms)\n"
  "    kip3 	= 0.0008    (mM) 	: Fink et al., 2000\n"
  "	kactip3 = 0.0003	(mM) 	: Fink et al., 2000  0.3e-3\n"
  "	konip3 	= 2.7		(/mM-ms): Fink et al., 2000  2.7/mM-ms\n"
  "	kinhip3 = 0.0002	(mM)    : Fink et al., 2000  0.2 e-3 mM\n"
  "\n"
  ":-------------------------------------------------------------------------------------------------------------------------------------------------------------------\n"
  "	:CICR(RYR) Parameters\n"
  "	:Schutter and Smolen, 1998\n"
  "	kcicr 	= 0.00198	(mM) :Lokuta et al., 2002\n"
  "	ktcicr 	= 0.0006  	(mM)\n"
  "	vcicr 	= 5e-7     	(/ms)\n"
  "\n"
  "	: SER Buffer\n"
  "	Kmer = 0.5 (mM)\n"
  "	Bmer = 10 (mM)\n"
  "	\n"
  ":----------------\n"
  ":------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n"
  "	:: Mitochondrial Modelling\n"
  "	\n"
  "	vmcu = 1.4468e-6	(mM/ms)  \n"
  "	kmcu = 606e-6 	(mM) 	     : 606 nM\n"
  "    nmcu = 2.3 (1) 				 : Shutov et al., 2013\n"
  "	\n"
  "	:: Mitochondrial NCX\n"
  "	vncx = 6e-5 	(mM/ms)\n"
  "	kna = 	8		(mM)	:Boyman et al., 2013  8  	(mM)\n"
  "	kncx = 	35e-3	(mM)	:Boyman et al., 2013  13e-3	(mM)\n"
  "\n"
  ":-------------------------------------------------------------------------------------------------------------------------------------------------------------------\n"
  "	: Mito Buffer\n"
  "	Kmmt = 0.01e-3 (mM) : Faville et al., 2008\n"
  "	Bmmt = 0.065 (mM)\n"
  "	\n"
  ":------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n"
  "	:PMCA Parameters\n"
  "	k1 = 3.74e7        (/mM-s) \n"
  "	k2 = .25e6      (/s)\n"
  "	k3 = .5e3       (/s)\n"
  "	k4 = 5e0        (/mM-s)\n"
  "	pump0 = 1.3725e-13 (mol/cm2)  : set to 0 in hoc if this pump not wanted\n"
  "}\n"
  "\n"
  "ASSIGNED{\n"
  "	celsius		(degC)\n"
  "	ica			(mA/cm2)\n"
  "	\n"
  "	cai			(mM)\n"
  "	cao			(mM)\n"
  "	caeri		(mM)\n"
  "	camti 		(mM)\n"
  "	\n"
  "	vol[NANN]	(1)\n"
  ":-------------------------------------------------------------------------------------------------------------------------------------------------------------------\n"
  "	:PMCA\n"
  "	ica_pmp (mA/cm2)\n"
  "	last_ica_pmp (mA/cm2)\n"
  "	parea    (um)\n"
  "	c1      (1+8 um4/ms)\n"
  "	c2      (1-10 um/ms)\n"
  "	c3      (1-10 um/ms)\n"
  "	c4      (1+8 um4/ms)\n"
  "\n"
  ":-------------------------------------------------------------------------------------------------------------------------------------------------------------------\n"
  "	:: IP3 Ashhad et al.	\n"
  "	ip3i		(mM)\n"
  "\n"
  ":-------------------------------------------------------------------------------------------------------------------------------------------------------------------\n"
  "	: Mitochondrial NCX (MNCX)\n"
  "	nai			(mM)\n"
  ":-------------------------------------------------------------------------------------------------------------------------------------------------------------------\n"
  "	: Net mitochondrial flux into the cell	\n"
  "	jmcu[NANN] 			(mM/ms)\n"
  "	jmncx[NANN]			(mM /ms)\n"
  "	jmito[NANN]  		(mM /ms)    \n"
  ":-------------------------------------------------------------------------------------------------------------------------------------------------------------------\n"
  "	::Fluxes from ER\n"
  "	jer[NANN]		(mM /ms)\n"
  "	jip3[NANN]		(mM /ms)\n"
  "	jserca[NANN]	(mM /ms)\n"
  "	jcicr[NANN]     (mM /ms) \n"
  "	\n"
  "	: Buffer Binding Ratio of Dye\n"
  "	: bbrdye[NANN]\n"
  "	\n"
  "	: SER and Mito Buffers\n"
  "	fmer[NANN]\n"
  "	fmmt[NANN]\n"
  "}\n"
  "\n"
  "CONSTANT{\n"
  "volo = 1e10 (um2) }\n"
  "\n"
  "\n"
  "STATE{\n"
  "	:PMCA (Calcium ATPase) Pump on the membrane\n"
  "	pump            (mol/cm2) <1e-16> \n"
  "	pumpca          (mol/cm2) <1e-16>\n"
  ":-------------------------------------------------------------------------------------------------------------------------------------------------------------------\n"
  "	ca[NANN]			(mM)     <1e-8>\n"
  "	caer[NANN]			(mM)\n"
  "	camt[NANN]			(mM)\n"
  "	caip3ri				(mM)\n"
  ":-------------------------------------------------------------------------------------------------------------------------------------------------------------------\n"
  "	hc[NANN]        (1)		:IP3 channels in closed state\n"
  "	ho[NANN]		(1)     :IP3 channels in open state\n"
  ":-------------------------------------------------------------------------------------------------------------------------------------------------------------------\n"
  "	Ln[NANN]		(mM/ms)	\n"
  "}\n"
  "\n"
  "LOCAL factors_done\n"
  "\n"
  "INITIAL{ LOCAL total\n"
  "	\n"
  "	if (factors_done==0) {\n"
  "		factors_done= 1\n"
  "		factors()\n"
  "	}\n"
  ":-------------------------------------------------------------------------------------------------------------------------------------------------------------------\n"
  "	: initializing intracellular Ca concentration\n"
  "	cai = cai0	\n"
  ":-------------------------------------------------------------------------------------------------------------------------------------------------------------------	   \n"
  "	:CaATPase Pump\n"
  "	parms() \n"
  "	parea = PI*diam\n"
  "	pump = pump0\n"
  "	pumpca = cai*pump*k1/k2\n"
  "	total = pumpca + pump\n"
  "	if (total > 1e-9) {\n"
  "		pump = pump*(pump/total)\n"
  "		pumpca = pumpca*(pump/total)\n"
  "	}\n"
  "	ica_pmp = 0\n"
  "	last_ica_pmp = 0\n"
  "\n"
  ":------------------------------------------------------------------------------------------------------------------------------------------------------------------\n"
  "\n"
  "	FROM i=0 TO NANN-1{\n"
  "		ca[i]=cai0          :Intracellular      Ca Shell Concentration initialization\n"
  "		caer[i] = caeri0	:Intracellular ER   Ca Shell Concentration initialization\n"
  "		camt[i] = camti0	:Intracellular Mito Ca Shell Concentration initialization\n"
  "		caip3ri = cai0\n"
  "		\n"
  "		: ER\n"
  "		jserca[i] = 0\n"
  "		jip3[i] = 0\n"
  "		jcicr[i] = 0\n"
  "		\n"
  "		:Mito\n"
  "		jmcu[i] = 0\n"
  "		jmncx[i] = 0\n"
  "		\n"
  "		:Mitochondirial Buffer parameter\n"
  "		fmmt[i] = 1/(1+(Kmmt*Bmmt)/(Kmmt+camt[i])^2)\n"
  "		\n"
  "		:SER Buffer parameter\n"
  "		fmer[i] = 1/(1+(Kmer*Bmer)/(Kmer+caer[i])^2)\n"
  "	}\n"
  "	caeri= caer[0]\n"
  "	camti = camt[0]\n"
  ":-------------------------------------------------------------------------------------------------------------------------------------------------------------------	   \n"
  "	:Balancing SR Fluxes modified from Fink et al., 2000	\n"
  "	FROM i=0 TO NANN-1 {\n"
  "    		 ho[i] = kinhip3/(ca[i]+kinhip3) 	        : Intial open IP3 channels \n"
  "    		 hc[i] = 1-ho[i] 	   			            : Intial closed IP3 channels\n"
  "\n"
  "			   jserca[i] = (-vmaxsr*ca[i]^2 / (ca[i]^2 + kpsr^2))\n"
  "			   \n"
  "			   jip3[i] = (jmaxsr*(1-(ca[i]/caer[i])) * ( (ip3i/(ip3i+kip3)) * (ca[i]/(ca[i]+kactip3)) * ho[i] )^3 )\n"
  "			   \n"
  "			   if(ca[i] > ktcicr){			\n"
  "					jcicr[i] = (vcicr* (ca[i]/(kcicr+ca[i])) * (caer[i]-ca[i]) ) \n"
  "				} else {\n"
  "					jcicr[i] = 0\n"
  "				} 			\n"
  "		\n"
  "			   \n"
  "			   jer[i] = jserca[i]+jip3[i]+jcicr[i]\n"
  "			   \n"
  "			   UNITSOFF\n"
  "			   jmcu[i] = (-vmcu*ca[i]^nmcu / (ca[i]^nmcu + kmcu^nmcu))\n"
  "			   UNITSON\n"
  "			   jmncx[i] = vncx*(nai^3/(kna^3 + nai^3))*(camt[i]/(kncx+camt[i]))\n"
  "			   \n"
  "			   jmito[i] = jmcu[i]+jmncx[i]\n"
  "			   \n"
  "			   Ln[i] = -(jserca[i]+jip3[i]+jcicr[i])/(1 - (ca[i]/caeri0))\n"
  "    		 }\n"
  "}\n"
  "\n"
  "BREAKPOINT{\n"
  "	SOLVE state METHOD sparse\n"
  "	 last_ica_pmp = ica_pmp\n"
  "     ica = ica_pmp\n"
  "}\n"
  "\n"
  "LOCAL frat[NANN]\n"
  "\n"
  "PROCEDURE factors(){\n"
  "	LOCAL r, dr2\n"
  "	r = 1/2		        :starts at edge (half diam)\n"
  "	dr2 = r/(NANN-1)/2	:half thickness of annulus\n"
  "	vol[0] = 0\n"
  "	frat[0] = 2*r\n"
  "	FROM i=0 TO NANN-2{\n"
  "		vol[i] = vol[i] + PI*(r-dr2/2)*2*dr2 :interior half\n"
  "		r = r - dr2\n"
  "		frat[i+1] = 2*PI*r/(2*dr2)		:exterior edge of annulus\n"
  "                                        :divided by distance between centers\n"
  "		r = r - dr2\n"
  "		vol[i+1] = PI*(r+dr2/2)*2*dr2   :outer half of annulus\n"
  "		}\n"
  "}\n"
  "\n"
  "LOCAL dsq, dsqvol,dsqvolmt,dsqvoler\n"
  "\n"
  "KINETIC state {\n"
  "	COMPARTMENT ii, (1+bbr)*diam*diam*vol[ii]*0.81 {ca} : cytoplasmic volume is 0.81 of total volume in each shell. IF CHANGED, ALSO CHANGE in ip3_diff.mod\n"
  "	:COMPARTMENT ii, (1+bbr+bbrdye[ii])*diam*diam*vol[ii]*0.81 {ca} : cytoplasmic volume is 0.81 of total cytoplasmic volume in each shell. IF CHANGED, ALSO CHANGE in ip3_diff.mod file\n"
  "	COMPARTMENT     (1+bbr)*diam*diam*vol[0]*0.81 {caip3ri}: caip3ri is the ca change by ca release from ip3rs in outermost shell\n"
  "	: COMPARTMENT     (1+bbr+bbrdye[0])*diam*diam*vol[0] *0.81 {caip3ri}: caip3ri is the ca release from ip3rs in outermost shell\n"
  "	\n"
  "	COMPARTMENT jj,	(1/fmer[jj])*diam*diam*vol[jj]*0.12 {caer} :SER volume in the cell : > 10 % of total cytosolic volume Verkhratsky 2002. IF CHANGED, ALSO CHANGE in ip3_diff.mod file\n"
  "	COMPARTMENT kk, (1/fmmt[kk])*diam*diam*vol[kk]*0.07 {camt} :Mitochondrial volume in the cell : 6.96% of total cytosolic volume-Yilmaz 2017. IF CHANGED, ALSO CHANGE in ip3_diff.mod file\n"
  "	COMPARTMENT (1e10)*parea {pump pumpca}  :Calcium ATPase Pump on the membrane\n"
  "	COMPARTMENT volo {cao}\n"
  "\n"
  ":-------------------------------------------------------------------------------------------------------------------------------------------------------------------     :all currents except pump \n"
  "	~ ca[0] << (-(ica - last_ica_pmp)*PI*diam*(1e4)*frat[0]/(2*F))\n"
  "	\n"
  "	:PMCA\n"
  "	:Calcium ATPase Pump on the membrane\n"
  "	 ~ ca[0] + pump <-> pumpca  (c1,c2)  \n"
  "	 ~ pumpca <-> pump + cao    (c3,c4)\n"
  "	  \n"
  "	ica_pmp = (1e-4) * 2*F*(f_flux - b_flux)/parea \n"
  "	  \n"
  "   \n"
  "	: Diffusion\n"
  "	 FROM i=0 TO NANN-2{\n"
  "		::Diffusion of Cytoplasmic Ca\n"
  "		 ~ ca[i] <-> ca[i+1]	(DCa*frat[i+1], DCa*frat[i+1])\n"
  "}\n"
  "         \n"
  "	 dsq = diam*diam\n"
  "     	 \n"
  "	 FROM i=0 TO NANN-1{\n"
  "		 dsqvol = dsq*vol[i]*0.81\n"
  "		 dsqvoler = dsq*vol[i]*0.12\n"
  ":---------------------------------------------------------------------------------------------------------------------------------------------\n"
  "		:::SERCA pump, IP3R, CICR(RYR)\n"
  "	\n"
  "		:: SERCA pump\n"
  "		jserca[i] = ((-vmaxsr*ca[i]^2 / (ca[i]^2 + kpsr^2)))\n"
  "		\n"
  "		~ ca[i] << (dsqvol*jserca[i])\n"
  "		~ caer[i] << (-dsqvoler*jserca[i])\n"
  "		\n"
  "		:: IP3 channel\n"
  "		~ hc[i] <-> ho[i]  (konip3*kinhip3, konip3*ca[i])\n"
  "		jip3[i] = (jmaxsr*(1-(ca[i]/caer[i])) * ( (ip3i/(ip3i+kip3)) * (ca[i]/(ca[i]+kactip3)) * ho[i] )^3 )\n"
  "        \n"
  "		~ ca[i] << (dsqvol*jip3[i])\n"
  "		~ caer[i] << (-dsqvoler*jip3[i])\n"
  "		\n"
  "		:Calcium release by ip3r at the periphery. For coupling with calcium activated chloride channel (CACC). See Jin et al., 2013,2015\n"
  "		if (i==0) {\n"
  "			~ caip3ri << (dsqvol*jip3[0])\n"
  "		}\n"
  "		\n"
  "		:CICR\n"
  "		if(ca[i] > ktcicr){			\n"
  "			jcicr[i] = (vcicr* (ca[i]/(kcicr+ca[i])) * (caer[i]-ca[i]) )\n"
  "			~ ca[i] << (dsqvol * jcicr[i])\n"
  "			~ caer[i] << (-dsqvoler * jcicr[i])\n"
  "		} else {\n"
  "			jcicr[i] = 0\n"
  "			~ ca[i] << (dsqvol*jcicr[i])\n"
  "			~ caer[i] << (-dsqvoler*jcicr[i])\n"
  "		} 			\n"
  "		\n"
  "		:Leak Channels ER\n"
  "		~ ca[i] << (Ln[i]*(1-ca[i]/caeri0)*dsqvol)\n"
  "		~ caer[i] << (-Ln[i]*(1-ca[i]/caeri0)*dsqvoler)\n"
  "		\n"
  "		jer[i] = jserca[i]+jip3[i]+jcicr[i]\n"
  "		\n"
  "		:SER Buffer\n"
  "		fmer[i] = 1/(1+(Kmer*Bmer)/(Kmer+caer[i])^2) \n"
  ":---------------------------------------------------------------------------------------------------------------------------------------------------\n"
  "		: Mitochondria\n"
  "		dsqvolmt = dsq*vol[i]*0.07\n"
  "		\n"
  "		:::Influx - Jinmt  via mcu\n"
  "		UNITSOFF\n"
  "		jmcu[i] = ((-vmcu*ca[i]^nmcu / (ca[i]^nmcu + kmcu^nmcu))*1/(camt[i]*1e3))\n"
  "		UNITSON\n"
  "		~ ca[i] << (dsqvol*jmcu[i])\n"
  "		\n"
  "		::Outflux - Joutmt  via mncx\n"
  "		jmncx[i] = vncx*(nai^3/(kna^3 + nai^3))*(camt[i]/(kncx+camt[i]))\n"
  "		~ ca[i] << (jmncx[i]*dsqvol)\n"
  "		\n"
  "		:Total flux\n"
  "		jmito[i] = jmcu[i]+jmncx[i]\n"
  "		\n"
  "		:Mitochondrial Ca change\n"
  "		~ camt[i] <<  (-(jmncx[i]+jmcu[i])*dsqvolmt) \n"
  "		\n"
  "		:Mito Buffer\n"
  "		fmmt[i] = 1/(1+(Kmmt*Bmmt)/(Kmmt+camt[i])^2)\n"
  "		\n"
  "		:: Fluorescent Dye binding ratio\n"
  "		:bbrdye[i] = kmdye*Bmdye/(kmdye+ca[i])^2\n"
  "	}\n"
  "	\n"
  "	cai = ca[0]\n"
  "	caeri= caer[0]\n"
  "	camti = camt[0]\n"
  "}\n"
  "\n"
  "PROCEDURE parms() {\n"
  "	parea = 2*PI*(diam/2)\n"
  "        c1 = (1e7)*parea * k1\n"
  "        c2 = (1e7)*parea * k2\n"
  "        c3 = (1e7)*parea * k3\n"
  "        c4 = (1e7)*parea * k4\n"
  "}\n"
  "\n"
  "\n"
  "COMMENT\n"
  "The combination of voltage independent current and calcium\n"
  "accumulation is more difficult because care must be taken not to count\n"
  "the pump current twice in the computation of the change ca[0].  Hence\n"
  "the usage of last_ica_pmp to subtract the pump portion of the total\n"
  "calcium current in ica so that its effect can be calculated implicitly\n"
  "via the reaction \"pumpca <-> pump + cao\".  This artifice makes the\n"
  "pumping much more stable than the assumption of constant pump current\n"
  "during the step.  Otherwise, ca[0] is prone to become negative and that\n"
  "crashes the simulation (especially the automatic computation of eca). \n"
  "Calcium currents that are inward are generally safe to compute in\n"
  "separate models. \n"
  "\n"
  "ENDCOMMENT\n"
  ;
#endif
