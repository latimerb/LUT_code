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
 
#define nrn_init _nrn_init__cacc
#define _nrn_initial _nrn_initial__cacc
#define nrn_cur _nrn_cur__cacc
#define _nrn_current _nrn_current__cacc
#define nrn_jacob _nrn_jacob__cacc
#define nrn_state _nrn_state__cacc
#define _net_receive _net_receive__cacc 
#define rates rates__cacc 
#define states states__cacc 
 
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
#define gbar _p[0]
#define icl _p[1]
#define g _p[2]
#define ninf _p[3]
#define ntau _p[4]
#define EC50 _p[5]
#define hc _p[6]
#define n _p[7]
#define caip3ri _p[8]
#define ecl _p[9]
#define Dn _p[10]
#define v _p[11]
#define _g _p[12]
#define _ion_ecl	*_ppvar[0]._pval
#define _ion_icl	*_ppvar[1]._pval
#define _ion_dicldv	*_ppvar[2]._pval
#define _ion_caip3ri	*_ppvar[3]._pval
 
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
 static void _hoc_rates(void);
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
 "setdata_cacc", _hoc_setdata,
 "rates_cacc", _hoc_rates,
 0, 0
};
 /* declare global and static user variables */
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "gbar_cacc", "S/cm2",
 "icl_cacc", "mA/cm2",
 "g_cacc", "mho/cm2",
 "ntau_cacc", "ms",
 "EC50_cacc", "mM",
 "hc_cacc", "1",
 0,0
};
 static double delta_t = 0.01;
 static double n0 = 0;
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
"cacc",
 "gbar_cacc",
 0,
 "icl_cacc",
 "g_cacc",
 "ninf_cacc",
 "ntau_cacc",
 "EC50_cacc",
 "hc_cacc",
 0,
 "n_cacc",
 0,
 0};
 static Symbol* _cl_sym;
 static Symbol* _caip3r_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 13, _prop);
 	/*initialize range parameters*/
 	gbar = 1e-06;
 	_prop->param = _p;
 	_prop->param_size = 13;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 5, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_cl_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ecl */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* icl */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dicldv */
 prop_ion = need_memb(_caip3r_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[3]._pval = &prop_ion->param[1]; /* caip3ri */
 
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

 void _cacc_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("cl", -1.0);
 	ion_reg("caip3r", 2.0);
 	_cl_sym = hoc_lookup("cl_ion");
 	_caip3r_sym = hoc_lookup("caip3r_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 13, 5);
  hoc_register_dparam_semantics(_mechtype, 0, "cl_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "cl_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "cl_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "caip3r_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 cacc /home/mizzou/LUT_code/single_cell/components/mechanisms/x86_64/cacc.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "CACC current for bladder small DRG neuron soma model";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rates(_threadargsprotocomma_ double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[1], _dlist1[1];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset = 0; {
   rates ( _threadargscomma_ v ) ;
   Dn = ( ninf - n ) / ntau ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
 rates ( _threadargscomma_ v ) ;
 Dn = Dn  / (1. - dt*( ( ( ( - 1.0 ) ) ) / ntau )) ;
  return 0;
}
 /*END CVODE*/
 static int states (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) { {
   rates ( _threadargscomma_ v ) ;
    n = n + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / ntau)))*(- ( ( ( ninf ) ) / ntau ) / ( ( ( ( - 1.0 ) ) ) / ntau ) - n) ;
   }
  return 0;
}
 
static int  rates ( _threadargsprotocomma_ double _lv ) {
    hc = - 0.3126 * exp ( - _lv / 81.02 ) + 2.086 ;
   EC50 = 0.39175 * exp ( - _lv / 38.307 ) + 0.468 ;
   if ( _lv < - 100.0 ) {
     hc = 1.012 ;
     EC50 = 5.798 ;
     }
   ninf = ( 1.0 / ( 1.0 + pow( ( EC50 * 1e-3 / caip3ri ) , hc ) ) ) ;
   ntau = 1.0 ;
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
 
static int _ode_count(int _type){ return 1;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ecl = _ion_ecl;
  caip3ri = _ion_caip3ri;
     _ode_spec1 (_p, _ppvar, _thread, _nt);
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 1; ++_i) {
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
  ecl = _ion_ecl;
  caip3ri = _ion_caip3ri;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_cl_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_cl_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_cl_sym, _ppvar, 2, 4);
   nrn_update_ion_pointer(_caip3r_sym, _ppvar, 3, 1);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  n = n0;
 {
   rates ( _threadargscomma_ v ) ;
   n = ninf ;
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
  ecl = _ion_ecl;
  caip3ri = _ion_caip3ri;
 initmodel(_p, _ppvar, _thread, _nt);
 }
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   g = gbar * n ;
   icl = g * ( v - ecl ) ;
   }
 _current += icl;

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
  ecl = _ion_ecl;
  caip3ri = _ion_caip3ri;
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dicl;
  _dicl = icl;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dicldv += (_dicl - icl)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_icl += icl ;
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
  ecl = _ion_ecl;
  caip3ri = _ion_caip3ri;
 {   states(_p, _ppvar, _thread, _nt);
  } }}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(n) - _p;  _dlist1[0] = &(Dn) - _p;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/home/mizzou/LUT_code/single_cell/components/mechanisms/modfiles/cacc.mod";
static const char* nmodl_file_text = 
  "TITLE CACC current for bladder small DRG neuron soma model\n"
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
  "	SUFFIX cacc\n"
  "	USEION cl READ ecl WRITE icl VALENCE -1\n"
  "	USEION caip3r READ caip3ri VALENCE 2\n"
  "	RANGE g, icl  \n"
  "	RANGE gbar, ninf, ntau, vhalf, sf1, EC50, hc\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(molar) = (1/liter)\n"
  "	(mM)	= (millimolar)\n"
  "	(uM)	= (micromolar)\n"
  "	(S)  	= (siemens)\n"
  "	(mA) 	= (milliamp)\n"
  "	(mV) 	= (millivolt)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "        gbar = 1e-6 (S/cm2)\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "        v       (mV)\n"
  "        icl		(mA/cm2)\n"
  "		g		(mho/cm2)\n"
  "		ninf\n"
  "		ntau 	(ms)\n"
  "		caip3ri	(mM)\n"
  "		celsius	(degC)\n"
  "		ecl 	(mV)\n"
  "		EC50 	(mM)\n"
  "		hc 		(1)\n"
  "}\n"
  "\n"
  "STATE {\n"
  "        n\n"
  "}\n"
  " \n"
  "BREAKPOINT {\n"
  "        SOLVE states METHOD cnexp\n"
  "		g = gbar*n\n"
  "		icl = g*(v - ecl)\n"
  "}\n"
  " \n"
  "INITIAL {\n"
  "		rates(v)\n"
  "		n = ninf\n"
  "}\n"
  "\n"
  "DERIVATIVE states {  \n"
  "        rates(v)\n"
  "        n' = (ninf-n)/ntau\n"
  "}\n"
  " \n"
  "PROCEDURE rates(v(mV)) {\n"
  "		 \n"
  "		 UNITSOFF\n"
  "        :\"n\" CACC activation 		\n"
  "		hc = -0.3126*exp(-v/81.02) + 2.086\n"
  "		EC50 = 0.39175*exp(-v/38.307) + 0.468\n"
  "		\n"
  "		if(v<-100){ \n"
  "			hc = 1.012\n"
  "			EC50 = 5.798\n"
  "		}\n"
  "		\n"
  "		ninf = (1/(1+(EC50*1e-3/caip3ri)^hc))\n"
  "		ntau = 1\n"
  "}\n"
  " \n"
  "UNITSON \n"
  ;
#endif
