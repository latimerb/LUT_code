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
 
#define nrn_init _nrn_init__ip3dif
#define _nrn_initial _nrn_initial__ip3dif
#define nrn_cur _nrn_cur__ip3dif
#define _nrn_current _nrn_current__ip3dif
#define nrn_jacob _nrn_jacob__ip3dif
#define nrn_state _nrn_state__ip3dif
#define _net_receive _net_receive__ip3dif 
#define factors factors__ip3dif 
#define state state__ip3dif 
 
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
#define ip3 (_p + 0)
#define ip3i _p[12]
#define Dip3 (_p + 13)
#define v _p[25]
#define _g _p[26]
#define _ion_ip3i	*_ppvar[0]._pval
#define _style_ip3	*((int*)_ppvar[1]._pvoid)
#define diam	*_ppvar[2]._pval
 
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
 /* declaration of user functions */
 static void _hoc_factors(void);
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
 "setdata_ip3dif", _hoc_setdata,
 "factors_ip3dif", _hoc_factors,
 0, 0
};
 #define _zfactors_done _thread[2]._pval[0]
 #define _zfrat (_thread[2]._pval + 1)
 #define _zdsq _thread[2]._pval[13]
 #define _zdsqvol _thread[2]._pval[14]
 /* declare global and static user variables */
 static int _thread1data_inuse = 0;
static double _thread1data[12];
#define _gth 3
#define DIP3 DIP3_ip3dif
 double DIP3 = 0.283;
#define ip3i0 ip3i0_ip3dif
 double ip3i0 = 0.00016;
#define kdegr kdegr_ip3dif
 double kdegr = 0.00014;
#define vol_ip3dif (_thread1data + 0)
#define vol (_thread[_gth]._pval + 0)
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "kdegr_ip3dif", "/ms",
 "DIP3_ip3dif", "um2/ms",
 "ip3i0_ip3dif", "mM",
 "ip3_ip3dif", "mM",
 0,0
};
 static double delta_t = 0.01;
 static double ip30 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "kdegr_ip3dif", &kdegr_ip3dif,
 "DIP3_ip3dif", &DIP3_ip3dif,
 "ip3i0_ip3dif", &ip3i0_ip3dif,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 "vol_ip3dif", vol_ip3dif, 12,
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
 
#define _cvode_ieq _ppvar[3]._i
 static void _ode_synonym(int, double**, Datum**);
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"ip3dif",
 0,
 0,
 "ip3_ip3dif[12]",
 0,
 0};
 static Symbol* _morphology_sym;
 static Symbol* _ip3_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 27, _prop);
 	/*initialize range parameters*/
 	_prop->param = _p;
 	_prop->param_size = 27;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_morphology_sym);
 	_ppvar[2]._pval = &prop_ion->param[0]; /* diam */
 prop_ion = need_memb(_ip3_sym);
 nrn_check_conc_write(_prop, prop_ion, 1);
 nrn_promote(prop_ion, 3, 0);
 	_ppvar[0]._pval = &prop_ion->param[1]; /* ip3i */
 	_ppvar[1]._pvoid = (void*)(&(prop_ion->dparam[0]._i)); /* iontype for ip3 */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 "ip3_ip3dif", 1e-06,
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

 void _ip3_diff_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("ip3", 1.0);
 	_morphology_sym = hoc_lookup("morphology");
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
  hoc_register_prop_size(_mechtype, 27, 4);
  hoc_register_dparam_semantics(_mechtype, 0, "ip3_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "#ip3_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "cvodeieq");
  hoc_register_dparam_semantics(_mechtype, 2, "diam");
 	nrn_writes_conc(_mechtype, 0);
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_synonym(_mechtype, _ode_synonym);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 ip3dif /home/mizzou/LUT_code/single_cell/components/mechanisms/x86_64/ip3_diff.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double FARADAY = 96485.3;
 static double PI = 3.14159;
 /*Top LOCAL _zfactors_done */
 /*Top LOCAL _zfrat [ 12 ] */
 /*Top LOCAL _zdsq , _zdsqvol */
static int _reset;
static char *modelname = "ip3 dyanmics for bladder small DRG neuron soma model";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int factors(_threadargsproto_);
 extern double *_nrn_thread_getelm();
 
#define _MATELM1(_row,_col) *(_nrn_thread_getelm(_so, _row + 1, _col + 1))
 
#define _RHS1(_arg) _rhs[_arg+1]
  
#define _linmat1  1
 static int _spth1 = 1;
 static int _cvspth1 = 0;
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[12], _dlist1[12]; static double *_temp1;
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
for(_i=0;_i<12;_i++){
  	_RHS1(_i) = -_dt1*(_p[_slist1[_i]] - _p[_dlist1[_i]]);
	_MATELM1(_i, _i) = _dt1;
      
} 
for (_i=0; _i < 12; _i++) {
  	_RHS1(_i + 0) *= ( diam * diam * vol [ ((int) _i ) ] * 0.81) ;
_MATELM1(_i + 0, _i + 0) *= ( diam * diam * vol [ ((int) _i ) ] * 0.81);  } }
 /* COMPARTMENT _li , diam * diam * vol [ ((int) _i ) ] * 0.81 {
     ip3 }
   */
 _zdsq = diam * diam ;
   {int  _li ;for ( _li = 0 ; _li <= 12 - 2 ; _li ++ ) {
     /* ~ ip3 [ _li ] <-> ip3 [ _li + 1 ] ( DIP3 * _zfrat [ _li + 1 ] , DIP3 * _zfrat [ _li + 1 ] )*/
 f_flux =  DIP3 * _zfrat [ _li + 1 ] * ip3 [ _li] ;
 b_flux =  DIP3 * _zfrat [ _li + 1 ] * ip3 [ _li + 1] ;
 _RHS1( 0 +  _li) -= (f_flux - b_flux);
 _RHS1( 0 +  _li + 1) += (f_flux - b_flux);
 
 _term =  DIP3 * _zfrat [ _li + 1 ] ;
 _MATELM1( 0 +  _li ,0 +  _li)  += _term;
 _MATELM1( 0 +  _li + 1 ,0 +  _li)  -= _term;
 _term =  DIP3 * _zfrat [ _li + 1 ] ;
 _MATELM1( 0 +  _li ,0 +  _li + 1)  -= _term;
 _MATELM1( 0 +  _li + 1 ,0 +  _li + 1)  += _term;
 /*REACTION*/
  } }
   {int  _li ;for ( _li = 0 ; _li <= 12 - 1 ; _li ++ ) {
     _zdsqvol = _zdsq * vol [ _li ] * 0.81 ;
     /* ~ ip3 [ _li ] <-> ip3i0 ( kdegr * _zdsqvol , kdegr * _zdsqvol )*/
 f_flux =  kdegr * _zdsqvol * ip3 [ _li] ;
 b_flux =  kdegr * _zdsqvol * ip3i0 ;
 _RHS1( 0 +  _li) -= (f_flux - b_flux);
 
 _term =  kdegr * _zdsqvol ;
 _MATELM1( 0 +  _li ,0 +  _li)  += _term;
 /*REACTION*/
  } }
   ip3i = ip3 [ 0 ] ;
     } return _reset;
 }
 
/*CVODE ode begin*/
 static int _ode_spec1(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset=0;{
 double b_flux, f_flux, _term; int _i;
 {int _i; for(_i=0;_i<12;_i++) _p[_dlist1[_i]] = 0.0;}
 /* COMPARTMENT _li , diam * diam * vol [ ((int) _i ) ] * 0.81 {
   ip3 }
 */
 _zdsq = diam * diam ;
 {int  _li ;for ( _li = 0 ; _li <= 12 - 2 ; _li ++ ) {
   /* ~ ip3 [ _li ] <-> ip3 [ _li + 1 ] ( DIP3 * _zfrat [ _li + 1 ] , DIP3 * _zfrat [ _li + 1 ] )*/
 f_flux =  DIP3 * _zfrat [ _li + 1 ] * ip3 [ _li] ;
 b_flux =  DIP3 * _zfrat [ _li + 1 ] * ip3 [ _li + 1] ;
 Dip3 [ _li] -= (f_flux - b_flux);
 Dip3 [ _li + 1] += (f_flux - b_flux);
 
 /*REACTION*/
  } }
 {int  _li ;for ( _li = 0 ; _li <= 12 - 1 ; _li ++ ) {
   _zdsqvol = _zdsq * vol [ _li ] * 0.81 ;
   /* ~ ip3 [ _li ] <-> ip3i0 ( kdegr * _zdsqvol , kdegr * _zdsqvol )*/
 f_flux =  kdegr * _zdsqvol * ip3 [ _li] ;
 b_flux =  kdegr * _zdsqvol * ip3i0 ;
 Dip3 [ _li] -= (f_flux - b_flux);
 
 /*REACTION*/
  } }
 ip3i = ip3 [ 0 ] ;
 for (_i=0; _i < 12; _i++) { _p[_dlist1[_i + 0]] /= ( diam * diam * vol [ ((int) _i ) ] * 0.81);}
   } return _reset;
 }
 
/*CVODE matsol*/
 static int _ode_matsol1(void* _so, double* _rhs, double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset=0;{
 double b_flux, f_flux, _term; int _i;
   b_flux = f_flux = 0.;
 {int _i; double _dt1 = 1.0/dt;
for(_i=0;_i<12;_i++){
  	_RHS1(_i) = _dt1*(_p[_dlist1[_i]]);
	_MATELM1(_i, _i) = _dt1;
      
} 
for (_i=0; _i < 12; _i++) {
  	_RHS1(_i + 0) *= ( diam * diam * vol [ ((int) _i ) ] * 0.81) ;
_MATELM1(_i + 0, _i + 0) *= ( diam * diam * vol [ ((int) _i ) ] * 0.81);  } }
 /* COMPARTMENT _li , diam * diam * vol [ ((int) _i ) ] * 0.81 {
 ip3 }
 */
 _zdsq = diam * diam ;
 {int  _li ;for ( _li = 0 ; _li <= 12 - 2 ; _li ++ ) {
 /* ~ ip3 [ _li ] <-> ip3 [ _li + 1 ] ( DIP3 * _zfrat [ _li + 1 ] , DIP3 * _zfrat [ _li + 1 ] )*/
 _term =  DIP3 * _zfrat [ _li + 1 ] ;
 _MATELM1( 0 +  _li ,0 +  _li)  += _term;
 _MATELM1( 0 +  _li + 1 ,0 +  _li)  -= _term;
 _term =  DIP3 * _zfrat [ _li + 1 ] ;
 _MATELM1( 0 +  _li ,0 +  _li + 1)  -= _term;
 _MATELM1( 0 +  _li + 1 ,0 +  _li + 1)  += _term;
 /*REACTION*/
  } }
 {int  _li ;for ( _li = 0 ; _li <= 12 - 1 ; _li ++ ) {
 _zdsqvol = _zdsq * vol [ _li ] * 0.81 ;
 /* ~ ip3 [ _li ] <-> ip3i0 ( kdegr * _zdsqvol , kdegr * _zdsqvol )*/
 _term =  kdegr * _zdsqvol ;
 _MATELM1( 0 +  _li ,0 +  _li)  += _term;
 /*REACTION*/
  } }
 ip3i = ip3 [ 0 ] ;
   } return _reset;
 }
 
/*CVODE end*/
 
static int _ode_count(int _type){ return 12;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ip3i = _ion_ip3i;
  ip3i = _ion_ip3i;
     _ode_spec1 (_p, _ppvar, _thread, _nt);
  _ion_ip3i = ip3i;
 }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 12; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 static void _ode_synonym(int _cnt, double** _pp, Datum** _ppd) { 
	double* _p; Datum* _ppvar;
 	int _i; 
	for (_i=0; _i < _cnt; ++_i) {_p = _pp[_i]; _ppvar = _ppd[_i];
 _ion_ip3i =  ip3 [ 0 ] ;
 }}
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _cvode_sparse_thread(&_thread[_cvspth1]._pvoid, 12, _dlist1, _p, _ode_matsol1, _ppvar, _thread, _nt);
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
  ip3i = _ion_ip3i;
  ip3i = _ion_ip3i;
 _ode_matsol_instance1(_threadargs_);
 }}
 
static void _thread_mem_init(Datum* _thread) {
   _thread[2]._pval = (double*)ecalloc(15, sizeof(double));
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
   nrn_update_ion_pointer(_ip3_sym, _ppvar, 0, 1);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
 for (_i=0; _i<12; _i++) ip3[_i] = ip30;
 {
   if ( _zfactors_done  == 0.0 ) {
     _zfactors_done = 1.0 ;
     factors ( _threadargs_ ) ;
     }
   ip3i = ip3i0 ;
   {int  _li ;for ( _li = 0 ; _li <= 12 - 1 ; _li ++ ) {
     ip3 [ _li ] = ip3i ;
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
  ip3i = _ion_ip3i;
  ip3i = _ion_ip3i;
 initmodel(_p, _ppvar, _thread, _nt);
  _ion_ip3i = ip3i;
  nrn_wrote_conc(_ip3_sym, (&(_ion_ip3i)) - 1, _style_ip3);
}
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{
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
  ip3i = _ion_ip3i;
  ip3i = _ion_ip3i;
 {  sparse_thread(&_thread[_spth1]._pvoid, 12, _slist1, _dlist1, _p, &t, dt, state, _linmat1, _ppvar, _thread, _nt);
     if (secondorder) {
    int _i;
    for (_i = 0; _i < 12; ++_i) {
      _p[_slist1[_i]] += dt*_p[_dlist1[_i]];
    }}
 } {
   }
  _ion_ip3i = ip3i;
}}
 dt = _dtsav;
}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 for(_i=0;_i<12;_i++){_slist1[0+_i] = (ip3 + _i) - _p;  _dlist1[0+_i] = (Dip3 + _i) - _p;}
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/home/mizzou/LUT_code/single_cell/components/mechanisms/modfiles/ip3_diff.mod";
static const char* nmodl_file_text = 
  "TITLE ip3 dyanmics for bladder small DRG neuron soma model\n"
  "\n"
  ": Adapted from Fink et al., 2000\n"
  ": Author: Darshan Mandge (darshanmandge@iitb.ac.in)\n"
  ": Computational Neurophysiology Lab\n"
  ": Indian Institute of Technology Bombay, India \n"
  "\n"
  ":For details refer: \n"
  ":A biophysically detailed computational model of bladder small DRG neuron soma \n"
  ":Darshan Mandge and Rohit Manchanda, PLOS Computational Biology (2018)\n"
  "\n"
  "NEURON {\n"
  "	 SUFFIX ip3dif\n"
  " 	 USEION ip3 READ ip3i WRITE ip3i VALENCE 1\n"
  "  	 GLOBAL vol, DIP3, ip3i0, kdegr\n"
  "     RANGE ip3i\n"
  "	 THREADSAFE\n"
  "}\n"
  "\n"
  "DEFINE NANN 12:  :This needs to be changed if NANN in cadyn.mod is changed.\n"
  "\n"
  "UNITS {\n"
  "  	(molar) = (1/liter)\n"
  "  	(mM)    = (millimolar)\n"
  "  	(uM)    = (micromolar)\n"
  "  	(um)    = (micron)\n"
  "  	(mA)    = (milliamp)\n"
  "  	FARADAY = (faraday)  (coulomb)\n"
  "  	PI      = (pi)       (1)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "  	kdegr = 0.14e-3 (/ms)  : degredation rate Fink et al.,2000\n"
  "  	DIP3 = 0.283(um2/ms)\n"
  "  	ip3i0 = 0.16e-3 (mM)   : [IP3]0  initial and resting ip3i conc\n"
  "}\n"
  "\n"
  "\n"
  "ASSIGNED {\n"
  "  	diam      (um)\n"
  " 	ip3i      (mM)\n"
  "  	vol[NANN]  		: numeric value of vol[i] equals the volume \n"
  "					: of annulus i of a 1um diameter cylinder\n"
  "					: multiply by diam^2 to get volume per um length\n"
  "}\n"
  "\n"
  "STATE {\n"
  "  	ip3[NANN]       (mM) <1e-6>\n"
  "}\n"
  "\n"
  "LOCAL factors_done\n"
  "\n"
  "BREAKPOINT { SOLVE state METHOD sparse}\n"
  "\n"
  "\n"
  "INITIAL {\n"
  "   	if (factors_done == 0) {   : flag becomes 1 in the first segment\n"
  "     	factors_done = 1       :   all subsequent segments will have\n"
  "      	factors()              :   vol = 0 unless vol is GLOBAL\n"
  "       }\n"
  "\n"
  "  	ip3i = ip3i0\n"
  "  	FROM i=0 TO NANN-1 {\n"
  "    	ip3[i] = ip3i\n"
  "	}\n"
  "}\n"
  "\n"
  "\n"
  "LOCAL frat[NANN]  : scales the rate constants for model geometry\n"
  "\n"
  "PROCEDURE factors() {\n"
  "  	LOCAL r, dr2\n"
  "  	r = 1/2                : starts at edge (half diam)\n"
  "  	dr2 = r/(NANN-1)/2     : full thickness of outermost annulus,\n"
  "						   : half thickness of all other annuli\n"
  "  	vol[0] = 0\n"
  "  	frat[0] = 2*r\n"
  " 	 FROM i=0 TO NANN-2 {\n"
  "    		vol[i] = vol[i] + PI*(r-dr2/2)*2*dr2  : interior half\n"
  "    		r = r - dr2\n"
  "    		frat[i+1] = 2*PI*r/(2*dr2)  : outer radius of annulus\n"
  "										: div by distance between centers\n"
  "   		 r = r - dr2\n"
  "    		vol[i+1] = PI*(r+dr2/2)*2*dr2  : outer half of annulus\n"
  "  	}\n"
  "}\n"
  "\n"
  "LOCAL dsq, dsqvol  : can't define local variable in KINETIC block\n"
  "\n"
  "KINETIC state {\n"
  "  	COMPARTMENT i, diam*diam*vol[i]*0.81 {ip3 ip3i0}  : cytoplasmic volume is 0.81 of total cytoplasmic volume in each shell\n"
  "	\n"
  "	dsq = diam*diam\n"
  "	\n"
  "  	FROM i=0 TO NANN-2 {\n"
  "   		 ~ ip3[i] <-> ip3[i+1]  (DIP3*frat[i+1], DIP3*frat[i+1])\n"
  "  	}\n"
  "\n"
  "  	\n"
  "  	FROM i=0 TO NANN-1 {\n"
  "    		dsqvol = dsq*vol[i]*0.81 						 : cytoplasmic volume is 0.81 of the total cell volume\n"
  "    		~ ip3[i] <-> ip3i0  (kdegr*dsqvol, kdegr*dsqvol) : ip3i0 is inactive ip3\n"
  "  	}\n"
  "\n"
  "  	ip3i = ip3[0]\n"
  "}\n"
  ;
#endif
