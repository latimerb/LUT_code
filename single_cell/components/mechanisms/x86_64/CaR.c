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
 
#define nrn_init _nrn_init__CaR
#define _nrn_initial _nrn_initial__CaR
#define nrn_cur _nrn_cur__CaR
#define _nrn_current _nrn_current__CaR
#define nrn_jacob _nrn_jacob__CaR
#define nrn_state _nrn_state__CaR
#define _net_receive _net_receive__CaR 
#define rates rates__CaR 
#define state state__CaR 
 
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
#define pmax _p[0]
#define ica _p[1]
#define mtau _p[2]
#define minf _p[3]
#define hinf _p[4]
#define htau _p[5]
#define h2tau _p[6]
#define m _p[7]
#define h _p[8]
#define h2 _p[9]
#define cai _p[10]
#define cao _p[11]
#define Dm _p[12]
#define Dh _p[13]
#define Dh2 _p[14]
#define v _p[15]
#define _g _p[16]
#define _ion_cai	*_ppvar[0]._pval
#define _ion_cao	*_ppvar[1]._pval
#define _ion_ica	*_ppvar[2]._pval
#define _ion_dicadv	*_ppvar[3]._pval
 
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
 static void _hoc_ghk(void);
 static void _hoc_rates(void);
 static void _hoc_table_tabhtau2(void);
 static void _hoc_tabhtau2(void);
 static void _hoc_table_tabhtau(void);
 static void _hoc_tabhtau(void);
 static void _hoc_table_tabmtau(void);
 static void _hoc_tabmtau(void);
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
 "setdata_CaR", _hoc_setdata,
 "ghk_CaR", _hoc_ghk,
 "rates_CaR", _hoc_rates,
 "table_tabhtau2_CaR", _hoc_table_tabhtau2,
 "tabhtau2_CaR", _hoc_tabhtau2,
 "table_tabhtau_CaR", _hoc_table_tabhtau,
 "tabhtau_CaR", _hoc_tabhtau,
 "table_tabmtau_CaR", _hoc_table_tabmtau,
 "tabmtau_CaR", _hoc_tabmtau,
 0, 0
};
#define ghk ghk_CaR
#define table_tabhtau2 table_tabhtau2_CaR
#define tabhtau2 tabhtau2_CaR
#define table_tabhtau table_tabhtau_CaR
#define tabhtau tabhtau_CaR
#define table_tabmtau table_tabmtau_CaR
#define tabmtau tabmtau_CaR
 extern double ghk( _threadargsprotocomma_ double , double , double , double );
 extern double table_tabhtau2( );
 extern double tabhtau2( _threadargsprotocomma_ double );
 extern double table_tabhtau( );
 extern double tabhtau( _threadargsprotocomma_ double );
 extern double table_tabmtau( );
 extern double tabmtau( _threadargsprotocomma_ double );
 /* declare global and static user variables */
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "pmax_CaR", "cm/s",
 "ica_CaR", "mA/cm2",
 "mtau_CaR", "ms",
 "htau_CaR", "ms",
 "h2tau_CaR", "ms",
 0,0
};
 static double delta_t = 0.01;
 static double h20 = 0;
 static double h0 = 0;
 static double m0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 0,0
};
 static DoubVec hoc_vdoub[] = {
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
 
#define _cvode_ieq _ppvar[4]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"CaR",
 "pmax_CaR",
 0,
 "ica_CaR",
 "mtau_CaR",
 "minf_CaR",
 "hinf_CaR",
 "htau_CaR",
 "h2tau_CaR",
 0,
 "m_CaR",
 "h_CaR",
 "h2_CaR",
 0,
 0};
 static Symbol* _ca_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 17, _prop);
 	/*initialize range parameters*/
 	pmax = 1e-08;
 	_prop->param = _p;
 	_prop->param_size = 17;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 5, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_ca_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[0]._pval = &prop_ion->param[1]; /* cai */
 	_ppvar[1]._pval = &prop_ion->param[2]; /* cao */
 	_ppvar[2]._pval = &prop_ion->param[3]; /* ica */
 	_ppvar[3]._pval = &prop_ion->param[4]; /* _ion_dicadv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _CaR_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("ca", -10000.);
 	_ca_sym = hoc_lookup("ca_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 17, 5);
  hoc_register_dparam_semantics(_mechtype, 0, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 CaR /home/mizzou/LUT_code/single_cell/components/mechanisms/x86_64/CaR.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double FARADAY = 96489.0;
 static double R = 8.314;
static int _reset;
static char *modelname = "R-type calcium current for bladder small DRG neuron soma model";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rates(_threadargsprotocomma_ double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 
static void* _ptable_tabhtau2 = (void*)0;
 
static void* _ptable_tabhtau = (void*)0;
 
static void* _ptable_tabmtau = (void*)0;
 static int _slist1[3], _dlist1[3];
 static int state(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset = 0; {
   rates ( _threadargscomma_ v ) ;
   Dm = ( minf - m ) / mtau ;
   Dh = ( hinf - h ) / htau ;
   Dh2 = ( hinf - h2 ) / h2tau ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
 rates ( _threadargscomma_ v ) ;
 Dm = Dm  / (1. - dt*( ( ( ( - 1.0 ) ) ) / mtau )) ;
 Dh = Dh  / (1. - dt*( ( ( ( - 1.0 ) ) ) / htau )) ;
 Dh2 = Dh2  / (1. - dt*( ( ( ( - 1.0 ) ) ) / h2tau )) ;
  return 0;
}
 /*END CVODE*/
 static int state (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) { {
   rates ( _threadargscomma_ v ) ;
    m = m + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / mtau)))*(- ( ( ( minf ) ) / mtau ) / ( ( ( ( - 1.0 ) ) ) / mtau ) - m) ;
    h = h + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / htau)))*(- ( ( ( hinf ) ) / htau ) / ( ( ( ( - 1.0 ) ) ) / htau ) - h) ;
    h2 = h2 + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / h2tau)))*(- ( ( ( hinf ) ) / h2tau ) / ( ( ( ( - 1.0 ) ) ) / h2tau ) - h2) ;
   }
  return 0;
}
 
double ghk ( _threadargsprotocomma_ double _lv , double _lci , double _lco , double _lz ) {
   double _lghk;
 double _le , _lw ;
 _lw = _lv * ( .001 ) * _lz * FARADAY / ( R * ( celsius + 273.16 ) ) ;
   if ( fabs ( _lw ) > 1e-4 ) {
     _le = _lw / ( exp ( _lw ) - 1.0 ) ;
     }
   else {
     _le = 1.0 - _lw / 2.0 ;
     }
   _lghk = - ( .001 ) * _lz * FARADAY * ( _lco - _lci * exp ( _lw ) ) * _le ;
   
return _lghk;
 }
 
static void _hoc_ghk(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  ghk ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) , *getarg(3) , *getarg(4) );
 hoc_retpushx(_r);
}
 
double tabmtau ( _threadargsprotocomma_ double _lv ) {
 double _arg[1];
 _arg[0] = _lv;
 return hoc_func_table(_ptable_tabmtau, 1, _arg);
 }
/*  }
  */
 
static void _hoc_tabmtau(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  tabmtau ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 double table_tabmtau ( ) {
	hoc_spec_table(&_ptable_tabmtau, 1);
	return 0.;
}
 
static void _hoc_table_tabmtau(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  table_tabmtau (  );
 hoc_retpushx(_r);
}
 
double tabhtau ( _threadargsprotocomma_ double _lv ) {
 double _arg[1];
 _arg[0] = _lv;
 return hoc_func_table(_ptable_tabhtau, 1, _arg);
 }
/*  }
  */
 
static void _hoc_tabhtau(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  tabhtau ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 double table_tabhtau ( ) {
	hoc_spec_table(&_ptable_tabhtau, 1);
	return 0.;
}
 
static void _hoc_table_tabhtau(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  table_tabhtau (  );
 hoc_retpushx(_r);
}
 
double tabhtau2 ( _threadargsprotocomma_ double _lv ) {
 double _arg[1];
 _arg[0] = _lv;
 return hoc_func_table(_ptable_tabhtau2, 1, _arg);
 }
/*  }
  */
 
static void _hoc_tabhtau2(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  tabhtau2 ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 double table_tabhtau2 ( ) {
	hoc_spec_table(&_ptable_tabhtau2, 1);
	return 0.;
}
 
static void _hoc_table_tabhtau2(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  table_tabhtau2 (  );
 hoc_retpushx(_r);
}
 
static int  rates ( _threadargsprotocomma_ double _lv ) {
   minf = 1.0 / ( 1.0 + exp ( ( - 5.0 - _lv ) / 5.0 ) ) ;
   mtau = tabmtau ( _threadargscomma_ _lv ) ;
   hinf = 1.0 / ( 1.0 + exp ( ( _lv + 51.0 ) / 12.0 ) ) ;
   htau = tabhtau ( _threadargscomma_ _lv ) ;
   h2tau = tabhtau2 ( _threadargscomma_ _lv ) ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 rates ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 3;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  cai = _ion_cai;
  cao = _ion_cao;
     _ode_spec1 (_p, _ppvar, _thread, _nt);
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 3; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
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
  cai = _ion_cai;
  cao = _ion_cao;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_ca_sym, _ppvar, 0, 1);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 1, 2);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 2, 3);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 3, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  h2 = h20;
  h = h0;
  m = m0;
 {
   rates ( _threadargscomma_ v ) ;
   m = minf ;
   h = hinf ;
   h2 = hinf ;
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
  cai = _ion_cai;
  cao = _ion_cao;
 initmodel(_p, _ppvar, _thread, _nt);
 }
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   ica = pmax * m * ( 0.4 * h + 0.6 * h2 ) * ghk ( _threadargscomma_ v , cai , cao , 2.0 ) ;
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
  cai = _ion_cai;
  cao = _ion_cao;
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dica;
  _dica = ica;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dicadv += (_dica - ica)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ica += ica ;
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
  cai = _ion_cai;
  cao = _ion_cao;
 {   state(_p, _ppvar, _thread, _nt);
  } }}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(m) - _p;  _dlist1[0] = &(Dm) - _p;
 _slist1[1] = &(h) - _p;  _dlist1[1] = &(Dh) - _p;
 _slist1[2] = &(h2) - _p;  _dlist1[2] = &(Dh2) - _p;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/home/mizzou/LUT_code/single_cell/components/mechanisms/modfiles/CaR.mod";
static const char* nmodl_file_text = 
  "TITLE R-type calcium current for bladder small DRG neuron soma model\n"
  "\n"
  ": Author: Darshan Mandge (darshanmandge@iitb.ac.in)\n"
  ": Computational Neurophysiology Lab\n"
  ": Indian Institute of Technology Bombay, India \n"
  "\n"
  ": For details refer: \n"
  ": A biophysically detailed computational model of bladder small DRG neuron soma \n"
  ": Darshan Mandge and Rohit Manchanda, PLOS Computational Biology (2018)\n"
  "\n"
  "NEURON {\n"
  "	SUFFIX CaR\n"
  "	USEION ca READ cai,cao WRITE ica\n"
  "	RANGE minf, mtau, hinf, htau, h2tau, ica, pmax\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(mA)	= (milliamp)\n"
  "	(mV)	= (millivolt)\n"
  "	(mM)	= (milli/liter)\n"
  "	FARADAY = 96489 (coul)\n"
  "	R       = 8.314 (volt-coul/degC)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	v		(mV)\n"
  "	celsius	(degC)\n"
  "	cai		(mM)\n"
  "	cao		(mM)\n"
  "	pmax =  1e-8	(cm/s)\n"
  "}\n"
  "\n"
  "STATE {\n"
  "	m 	:	activation\n"
  "	h 	:	slow inactivation\n"
  "	h2	: 	fast inactivation\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	ica		(mA/cm2)\n"
  "	mtau		(ms)\n"
  "	minf\n"
  "	hinf\n"
  "	htau		(ms)\n"
  "    h2tau		(ms)\n"
  "}\n"
  "\n"
  "BREAKPOINT { \n"
  "	SOLVE state METHOD cnexp\n"
  "	ica = pmax*m*(0.4*h+0.6*h2)*ghk(v,cai,cao,2)\n"
  "}\n"
  "\n"
  "DERIVATIVE state {\n"
  "	rates(v)\n"
  "	m'= (minf-m) / mtau\n"
  "	h'= (hinf-h) / htau\n"
  "    h2'= (hinf-h2) / h2tau\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "	rates(v)\n"
  "	m = minf\n"
  "	h = hinf\n"
  "    h2= hinf\n"
  "}\n"
  "\n"
  ": For GHK Current equation (adapted from Stepheen and Manchanda, 2009)\n"
  "FUNCTION ghk( v(mV), ci(mM), co(mM), z)  (millicoul/cm3) { LOCAL e, w\n"
  "        w = v * (.001) * z*FARADAY / (R*(celsius+273.16))\n"
  "        if (fabs(w)>1e-4) \n"
  "          { e = w / (exp(w)-1) }\n"
  "        else : denominator is small -> Taylor series\n"
  "          { e = 1-w/2 }\n"
  "        ghk = - (.001) * z*FARADAY * (co-ci*exp(w)) * e\n"
  "}\n"
  "\n"
  "UNITSOFF\n"
  "\n"
  "FUNCTION_TABLE tabmtau(v(mV)) (ms) \n"
  "FUNCTION_TABLE tabhtau(v(mV)) (ms)\n"
  "FUNCTION_TABLE tabhtau2(v(mV)) (ms)\n"
  "\n"
  "\n"
  "PROCEDURE rates(v(mV)) { 	\n"
  "	minf = 1 / (1+exp((-5-v)/5)) 	: Li et al., 2007\n"
  "	mtau = tabmtau(v)				: Hilaire et al., 1997\n"
  "	\n"
  "	hinf = 1 / (1+exp((v+51)/12)) 	: Hilaire et al., 1997\n"
  "	    \n"
  "    htau  = tabhtau(v)	: Hilaire et al., 1997\n"
  "	h2tau = tabhtau2(v) : Hilaire et al., 1997\n"
  "}\n"
  "\n"
  "UNITSON \n"
  ;
#endif
