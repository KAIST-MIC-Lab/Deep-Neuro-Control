#ifndef __FUSION_P_H__
#define __FUSION_P_H__
#include "monty.h"
#include "mosektask_p.h"
#include "list"
#include "vector"
#include "unordered_map"
#include "fusion.h"
namespace mosek
{
  namespace fusion
  {
    struct p_Disjunction
    {
      Disjunction * _pubthis;
      static mosek::fusion::p_Disjunction* _get_impl(mosek::fusion::Disjunction * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_Disjunction * _get_impl(mosek::fusion::Disjunction::t _inst) { return _get_impl(_inst.get()); }
      p_Disjunction(Disjunction * _pubthis);
      virtual ~p_Disjunction() { /* std::cout << "~p_Disjunction" << std::endl;*/ };
      int64_t id{};

      virtual void destroy();

      static Disjunction::t _new_Disjunction(int64_t _7_id);
      void _initialize(int64_t _7_id);
    }; // struct Disjunction;

    struct p_DisjunctionTerms
    {
      DisjunctionTerms * _pubthis;
      static mosek::fusion::p_DisjunctionTerms* _get_impl(mosek::fusion::DisjunctionTerms * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_DisjunctionTerms * _get_impl(mosek::fusion::DisjunctionTerms::t _inst) { return _get_impl(_inst.get()); }
      p_DisjunctionTerms(DisjunctionTerms * _pubthis);
      virtual ~p_DisjunctionTerms() { /* std::cout << "~p_DisjunctionTerms" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Term >,1 > > terms{};

      virtual void destroy();

      static DisjunctionTerms::t _new_DisjunctionTerms(monty::rc_ptr< ::mosek::fusion::DisjunctionTerms > _8_terms1,monty::rc_ptr< ::mosek::fusion::ExprRangeDomain > _9_term);
      void _initialize(monty::rc_ptr< ::mosek::fusion::DisjunctionTerms > _8_terms1,monty::rc_ptr< ::mosek::fusion::ExprRangeDomain > _9_term);
      static DisjunctionTerms::t _new_DisjunctionTerms(monty::rc_ptr< ::mosek::fusion::DisjunctionTerms > _10_terms1,monty::rc_ptr< ::mosek::fusion::ExprLinearDomain > _11_term);
      void _initialize(monty::rc_ptr< ::mosek::fusion::DisjunctionTerms > _10_terms1,monty::rc_ptr< ::mosek::fusion::ExprLinearDomain > _11_term);
      static DisjunctionTerms::t _new_DisjunctionTerms(monty::rc_ptr< ::mosek::fusion::DisjunctionTerms > _12_terms1,monty::rc_ptr< ::mosek::fusion::DisjunctionTerms > _13_term2);
      void _initialize(monty::rc_ptr< ::mosek::fusion::DisjunctionTerms > _12_terms1,monty::rc_ptr< ::mosek::fusion::DisjunctionTerms > _13_term2);
      static DisjunctionTerms::t _new_DisjunctionTerms(monty::rc_ptr< ::mosek::fusion::DisjunctionTerms > _14_term1,std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Term >,1 > > _15_term2);
      void _initialize(monty::rc_ptr< ::mosek::fusion::DisjunctionTerms > _14_term1,std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Term >,1 > > _15_term2);
      static DisjunctionTerms::t _new_DisjunctionTerms(monty::rc_ptr< ::mosek::fusion::DisjunctionTerms > _19_term1,monty::rc_ptr< ::mosek::fusion::Term > _20_term2);
      void _initialize(monty::rc_ptr< ::mosek::fusion::DisjunctionTerms > _19_term1,monty::rc_ptr< ::mosek::fusion::Term > _20_term2);
      static DisjunctionTerms::t _new_DisjunctionTerms(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Term >,1 > > _22_terms);
      void _initialize(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Term >,1 > > _22_terms);
      static DisjunctionTerms::t _new_DisjunctionTerms(monty::rc_ptr< ::mosek::fusion::ExprRangeDomain > _24_term);
      void _initialize(monty::rc_ptr< ::mosek::fusion::ExprRangeDomain > _24_term);
      static DisjunctionTerms::t _new_DisjunctionTerms(monty::rc_ptr< ::mosek::fusion::ExprLinearDomain > _25_term);
      void _initialize(monty::rc_ptr< ::mosek::fusion::ExprLinearDomain > _25_term);
    }; // struct DisjunctionTerms;

    struct p_Term : public /*implements*/ virtual ::mosek::fusion::ExprDomain
    {
      Term * _pubthis;
      static mosek::fusion::p_Term* _get_impl(mosek::fusion::Term * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_Term * _get_impl(mosek::fusion::Term::t _inst) { return _get_impl(_inst.get()); }
      p_Term(Term * _pubthis);
      virtual ~p_Term() { /* std::cout << "~p_Term" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::DJCDomain >,1 > > domains{};
      std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Expression >,1 > > exprs{};

      virtual void destroy();

      static Term::t _new_Term(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Expression >,1 > > _26_elist,std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::DJCDomain >,1 > > _27_dlist);
      void _initialize(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Expression >,1 > > _26_elist,std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::DJCDomain >,1 > > _27_dlist);
      static Term::t _new_Term(monty::rc_ptr< ::mosek::fusion::Expression > _30_e,monty::rc_ptr< ::mosek::fusion::DJCDomain > _31_d);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _30_e,monty::rc_ptr< ::mosek::fusion::DJCDomain > _31_d);
      static Term::t _new_Term(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::SimpleTerm >,1 > > _32_t);
      void _initialize(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::SimpleTerm >,1 > > _32_t);
      virtual int32_t numaccterms() ;
      virtual int32_t numaccrows() ;
      virtual int32_t num() ;
      virtual int32_t size() ;
      virtual monty::rc_ptr< ::mosek::fusion::Term > __mosek_2fusion_2Term__toDJCTerm() ;
      virtual monty::rc_ptr< ::mosek::fusion::Term > __mosek_2fusion_2ExprDomain__toDJCTerm() { return __mosek_2fusion_2Term__toDJCTerm(); }
    }; // struct Term;

    struct p_SimpleTerm : public ::mosek::fusion::p_Term
    {
      SimpleTerm * _pubthis;
      static mosek::fusion::p_SimpleTerm* _get_impl(mosek::fusion::SimpleTerm * _inst){ return static_cast< mosek::fusion::p_SimpleTerm* >(mosek::fusion::p_Term::_get_impl(_inst)); }
      static mosek::fusion::p_SimpleTerm * _get_impl(mosek::fusion::SimpleTerm::t _inst) { return _get_impl(_inst.get()); }
      p_SimpleTerm(SimpleTerm * _pubthis);
      virtual ~p_SimpleTerm() { /* std::cout << "~p_SimpleTerm" << std::endl;*/ };

      virtual void destroy();

      static SimpleTerm::t _new_SimpleTerm(monty::rc_ptr< ::mosek::fusion::Expression > _41_e,monty::rc_ptr< ::mosek::fusion::DJCDomain > _42_d);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _41_e,monty::rc_ptr< ::mosek::fusion::DJCDomain > _42_d);
    }; // struct SimpleTerm;

    struct p_DJCDomain
    {
      DJCDomain * _pubthis;
      static mosek::fusion::p_DJCDomain* _get_impl(mosek::fusion::DJCDomain * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_DJCDomain * _get_impl(mosek::fusion::DJCDomain::t _inst) { return _get_impl(_inst.get()); }
      p_DJCDomain(DJCDomain * _pubthis);
      virtual ~p_DJCDomain() { /* std::cout << "~p_DJCDomain" << std::endl;*/ };
      mosek::fusion::DJCDomainType dom{};
      int32_t conedim{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};
      std::shared_ptr< monty::ndarray< double,1 > > par{};
      std::shared_ptr< monty::ndarray< double,1 > > b{};

      virtual void destroy();

      static DJCDomain::t _new_DJCDomain(std::shared_ptr< monty::ndarray< double,1 > > _43_b_,std::shared_ptr< monty::ndarray< double,1 > > _44_par_,std::shared_ptr< monty::ndarray< int32_t,1 > > _45_shape_,mosek::fusion::DJCDomainType _46_dom_);
      void _initialize(std::shared_ptr< monty::ndarray< double,1 > > _43_b_,std::shared_ptr< monty::ndarray< double,1 > > _44_par_,std::shared_ptr< monty::ndarray< int32_t,1 > > _45_shape_,mosek::fusion::DJCDomainType _46_dom_);
      static DJCDomain::t _new_DJCDomain(std::shared_ptr< monty::ndarray< double,1 > > _47_b_,std::shared_ptr< monty::ndarray< double,1 > > _48_par_,std::shared_ptr< monty::ndarray< int32_t,1 > > _49_shape_,int32_t _50_conedim_,mosek::fusion::DJCDomainType _51_dom_);
      void _initialize(std::shared_ptr< monty::ndarray< double,1 > > _47_b_,std::shared_ptr< monty::ndarray< double,1 > > _48_par_,std::shared_ptr< monty::ndarray< int32_t,1 > > _49_shape_,int32_t _50_conedim_,mosek::fusion::DJCDomainType _51_dom_);
      virtual int32_t numaccterms() ;
      virtual int32_t numaccrows() ;
      virtual int32_t size() ;
    }; // struct DJCDomain;

    struct p_DJC
    {
      DJC * _pubthis;
      static mosek::fusion::p_DJC* _get_impl(mosek::fusion::DJC * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_DJC * _get_impl(mosek::fusion::DJC::t _inst) { return _get_impl(_inst.get()); }
      p_DJC(DJC * _pubthis);
      virtual ~p_DJC() { /* std::cout << "~p_DJC" << std::endl;*/ };

      virtual void destroy();

      static  monty::rc_ptr< ::mosek::fusion::Term > ANDFromTerms(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Term >,1 > > _54_tlist);
      static  monty::rc_ptr< ::mosek::fusion::Term > AND(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::ExprDomain >,1 > > _61_elist);
      static  monty::rc_ptr< ::mosek::fusion::Term > AND(monty::rc_ptr< ::mosek::fusion::SimpleTerm > _63_s1,monty::rc_ptr< ::mosek::fusion::SimpleTerm > _64_s2,monty::rc_ptr< ::mosek::fusion::SimpleTerm > _65_s3);
      static  monty::rc_ptr< ::mosek::fusion::Term > AND(monty::rc_ptr< ::mosek::fusion::SimpleTerm > _66_s1,monty::rc_ptr< ::mosek::fusion::SimpleTerm > _67_s2);
      static  monty::rc_ptr< ::mosek::fusion::Term > AND(monty::rc_ptr< ::mosek::fusion::SimpleTerm > _68_s1);
      static  monty::rc_ptr< ::mosek::fusion::Term > AND(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::SimpleTerm >,1 > > _69_slist);
      static  monty::rc_ptr< ::mosek::fusion::SimpleTerm > term(monty::rc_ptr< ::mosek::fusion::Expression > _70_expr,monty::rc_ptr< ::mosek::fusion::RangeDomain > _71_dom);
      static  monty::rc_ptr< ::mosek::fusion::SimpleTerm > term(monty::rc_ptr< ::mosek::fusion::Variable > _82_x,monty::rc_ptr< ::mosek::fusion::RangeDomain > _83_dom);
      static  monty::rc_ptr< ::mosek::fusion::SimpleTerm > term(monty::rc_ptr< ::mosek::fusion::ExprRangeDomain > _84_exprdom);
      static  monty::rc_ptr< ::mosek::fusion::SimpleTerm > term(monty::rc_ptr< ::mosek::fusion::Expression > _85_expr,monty::rc_ptr< ::mosek::fusion::LinearDomain > _86_dom);
      static  monty::rc_ptr< ::mosek::fusion::SimpleTerm > term(monty::rc_ptr< ::mosek::fusion::ExprLinearDomain > _94_exprdom);
      static  monty::rc_ptr< ::mosek::fusion::SimpleTerm > term(monty::rc_ptr< ::mosek::fusion::Variable > _95_x,monty::rc_ptr< ::mosek::fusion::LinearDomain > _96_dom);
    }; // struct DJC;

    // mosek.fusion.BaseModel from file 'src/fusion/cxx/BaseModel_p.h'
    // namespace mosek::fusion
    struct p_BaseModel
    {
      p_BaseModel(BaseModel * _pubthis);
    
      void _initialize( monty::rc_ptr<BaseModel> m);
      void _initialize( const std::string & name,
                        const std::string & licpath);
    
      virtual ~p_BaseModel() { /* std::cout << "~p_BaseModel()" << std::endl;*/  }
    
      static p_BaseModel * _get_impl(Model * _inst) { return _inst->_impl; }
    
      //----------------------
    
      bool synched;
      std::string taskname;
    
      monty::rc_ptr<SolutionStruct> sol_itr;
      monty::rc_ptr<SolutionStruct> sol_itg;
      monty::rc_ptr<SolutionStruct> sol_bas;
    
      MSKsoltypee cursol;
      //---------------------
    
      std::unique_ptr<Task> task;
    
      //---------------------
      template<class T>
      using array_t = std::shared_ptr<monty::ndarray<T,1>>;
    
      virtual void clear_solutions() = 0;
      virtual void report_solution(SolutionType soltype,
                                   ProblemStatus prosta,
                                   SolutionStatus psolsta,
                                   SolutionStatus dsolsta,
                                   double pobj,
                                   double dobj,
                                   int32_t numvar,
                                   int32_t numcon,
                                   int32_t numbarelm,
                                   int32_t numacc,
                                   int32_t numaccelm,
                                   bool hasprimal,
                                   bool hasdual) = 0;
    
      void report_solution_get_xx(array_t<double> v);
      void report_solution_get_slx(array_t<double> v);
      void report_solution_get_sux(array_t<double> v);
      void report_solution_get_xc(array_t<double> v);
      void report_solution_get_slc(array_t<double> v);
      void report_solution_get_suc(array_t<double> v);
      void report_solution_get_barx(array_t<double> v);
      void report_solution_get_bars(array_t<double> v);
      void report_solution_get_accx(array_t<double> v);
      void report_solution_get_accy(array_t<double> v);
      void report_solution_get_accptr(array_t<int32_t> v);
    
      //---------------------
      void task_setLogHandler (const msghandler_t & handler);
      void task_setDataCallbackHandler (const datacbhandler_t & handler);
      void task_setCallbackHandler (const cbhandler_t & handler);
    
      int task_append_barvar(int size, int num);
    
      void task_djc_name   (int64_t index, const std::string & name);
      void task_var_name   (int index, const std::string & name);
      void task_con_name   (int index, const std::string & name);
      void task_barvar_name(int index, const std::string & name);
      void task_objectivename(         const std::string & name);
    
      void task_format_djc_names
      ( const std::shared_ptr<monty::ndarray<int64_t,1>> sub,
        const std::string                              & format,
        const std::shared_ptr<monty::ndarray<int,1>>     dims,
        const std::shared_ptr<monty::ndarray<std::shared_ptr<monty::ndarray<std::string,1>>>> names);
      void task_format_acc_names
      ( const std::shared_ptr<monty::ndarray<int64_t,1>> sub,
        const std::string                              & format,
        const std::shared_ptr<monty::ndarray<int,1>>     dims,
        const std::shared_ptr<monty::ndarray<std::shared_ptr<monty::ndarray<std::string,1>>>> names);
      void task_format_var_names
      ( const std::shared_ptr<monty::ndarray<int,1>>     subj,
        const std::string                              & format,
        const std::shared_ptr<monty::ndarray<int,1>>     dims,
        const std::shared_ptr<monty::ndarray<int64_t,1>> sp,
        const std::shared_ptr<monty::ndarray<std::shared_ptr<monty::ndarray<std::string,1>>>> names);
      void task_format_con_names
      ( const std::shared_ptr<monty::ndarray<int,1>>     subi,
        const std::string                              & format,
        const std::shared_ptr<monty::ndarray<int,1>>     dims,
        const std::shared_ptr<monty::ndarray<int64_t,1>> sp,
        const std::shared_ptr<monty::ndarray<std::shared_ptr<monty::ndarray<std::string,1>>>> names);
      void task_format_barvar_names
      ( const std::shared_ptr<monty::ndarray<int,1>>     subj,
        const std::string                              & format,
        const std::shared_ptr<monty::ndarray<int,1>>     dims,
        const std::shared_ptr<monty::ndarray<std::shared_ptr<monty::ndarray<std::string,1>>>> names);
    
      void task_break_solve();
    
      //--------------------------
    
      int task_numvar();
      int task_numcon();
      int task_numbarvar();
      int task_numacc();
      int task_numdjc();
      int task_numafe();
    
      //--------------------------
    
      void task_put_param(const std::string & name, const std::string & value);
      void task_put_param(const std::string & name, int    value);
      void task_put_param(const std::string & name, double value);
    
      double    task_get_dinf (const std::string & name);
      int       task_get_iinf (const std::string & name);
      int64_t task_get_liinf(const std::string & name);
    
      //--------------------------
    
      void task_con_putboundlist_fr(const std::shared_ptr<monty::ndarray<int,1>> idxs);
      void task_con_putboundlist_lo(const std::shared_ptr<monty::ndarray<int,1>> idxs, const std::shared_ptr<monty::ndarray<double,1>> & rhs);
      void task_con_putboundlist_up(const std::shared_ptr<monty::ndarray<int,1>> idxs, const std::shared_ptr<monty::ndarray<double,1>> & rhs);
      void task_con_putboundlist_fx(const std::shared_ptr<monty::ndarray<int,1>> idxs, const std::shared_ptr<monty::ndarray<double,1>> & rhs);
      void task_con_putboundlist_ra(const std::shared_ptr<monty::ndarray<int,1>> idxs, const std::shared_ptr<monty::ndarray<double,1>> & lb ,
                                    const std::shared_ptr<monty::ndarray<double,1>> & ub );
    
      void task_var_putboundlist_fr(const std::shared_ptr<monty::ndarray<int,1>> idxs);
      void task_var_putboundlist_lo(const std::shared_ptr<monty::ndarray<int,1>> idxs, const std::shared_ptr<monty::ndarray<double,1>> & rhs);
      void task_var_putboundlist_up(const std::shared_ptr<monty::ndarray<int,1>> idxs, const std::shared_ptr<monty::ndarray<double,1>> & rhs);
      void task_var_putboundlist_fx(const std::shared_ptr<monty::ndarray<int,1>> idxs, const std::shared_ptr<monty::ndarray<double,1>> & rhs);
      void task_var_putboundlist_ra(const std::shared_ptr<monty::ndarray<int,1>> idxs, const std::shared_ptr<monty::ndarray<double,1>> & lb ,
                                    const std::shared_ptr<monty::ndarray<double,1>> & ub );
    
      void task_var_putintlist(const std::shared_ptr<monty::ndarray<int,1>> & idxs);
      void task_var_putcontlist(const std::shared_ptr<monty::ndarray<int,1>> & idxs);
    
      //--------------------------
    
      void task_putbararowlist(const std::shared_ptr<monty::ndarray<int,1>>       subi,
                               const std::shared_ptr<monty::ndarray<int64_t,1>> ptr,
                               const std::shared_ptr<monty::ndarray<int,1>>       subj,
                               const std::shared_ptr<monty::ndarray<int64_t,1>> matidx);
    
      void task_putbaraijlist(const std::shared_ptr<monty::ndarray<int,1>> subi,
                              const std::shared_ptr<monty::ndarray<int,1>> subj,
                              std::shared_ptr<monty::ndarray<int64_t,1>> matidx);
    
      void task_putbarc(const std::shared_ptr<monty::ndarray<int,1>> subj,
                        const std::shared_ptr<monty::ndarray<int,1>> subl,
                        const std::shared_ptr<monty::ndarray<int,1>> subk,
                        const std::shared_ptr<monty::ndarray<double,1>> val);
    
      std::shared_ptr<monty::ndarray<int64_t,1>> task_appendsymmatlist (const std::shared_ptr<monty::ndarray<int,1>>       & dim,
                                                                        const std::shared_ptr<monty::ndarray<int64_t,1>> & nz,
                                                                        const std::shared_ptr<monty::ndarray<int,1>>       & subk,
                                                                        const std::shared_ptr<monty::ndarray<int,1>>       & subl,
                                                                        const std::shared_ptr<monty::ndarray<double,1>>    & val);
      int  task_barvar_dim(int j);
      void task_putbaraij (int i, int j, int k);
      void task_putbaraij (int i, int j, const std::shared_ptr<monty::ndarray<int,1>> & k);
      void task_putbarcj  (int j, int k);
      void task_putbarcj  (int j,        const std::shared_ptr<monty::ndarray<int,1>> & k);
      int  task_barvardim (int index);
    
      int task_append_var(int num);
      int task_append_con(int num);
    
      void task_cleararowlist(const std::shared_ptr<monty::ndarray<int,1>> & idxs);
      void task_clearacollist(const std::shared_ptr<monty::ndarray<int,1>> & idxs);
    
      void task_putarowlist(
        const std::shared_ptr<monty::ndarray<int,1>>       & idxs,
        const std::shared_ptr<monty::ndarray<int64_t,1>> & ptrb,
        const std::shared_ptr<monty::ndarray<int,1>>       & subj,
        const std::shared_ptr<monty::ndarray<double,1>>    & cof);
      void task_putaijlist(
        const std::shared_ptr<monty::ndarray<int,1>>       & subi,
        const std::shared_ptr<monty::ndarray<int,1>>       & subj,
        const std::shared_ptr<monty::ndarray<double,1>>    & cof,
        int64_t                           num);
    
      void task_setnumvar(int num);
      void task_cleanup(int oldnum, int oldnumcon, int oldnumcone, int oldnumbarvar);
      void task_putoptserver_host(const std::string & addr);
      void report_task_solution(MSKsoltypee st, int numvar, int numcon, int numbarelm, int64_t numacc, int64_t numaccelm);
    
      void task_solve(bool remote, const std::string & server, const std::string & port);
      void task_post_solve();
      static std::shared_ptr<monty::ndarray<SolverStatus,1>>  env_solve_batch(bool israce, 
                                                                              double timelimit, 
                                                                              int numthreads, 
                                                                              std::shared_ptr<monty::ndarray<Model::t,1>> & models);
    
      void task_putobjective(
        bool                             maximize,
        const std::shared_ptr<monty::ndarray<int,1>>    & subj    ,
        const std::shared_ptr<monty::ndarray<double,1>> & cof     ,
        double                           cfix    );
    
      void task_putclist(
        const std::shared_ptr<monty::ndarray<int,1>>    & subj,
        const std::shared_ptr<monty::ndarray<double,1>> & cof);
    
    
      void task_putobjectivename(const std::string & name);
    
      void task_write(const std::string & filename);
      void task_write_stream(const std::string & ext, std::ostream & stream);
      void task_dump (const std::string & filename);
      void task_analyze_problem(int detail);
    
      MSKtask_t task_get();
      MSKtask_t __mosek_2fusion_2BaseModel__task_get();
    
      MSKtask_t __mosek_2fusion_2BaseModel__clonetask();
      MSKtask_t clonetask();
    
    
      void dispose();
    
      void task_putxx_slice(SolutionType which, int first, int last, std::shared_ptr<monty::ndarray<double,1>> & xx);
    
      static void env_syeig (int n, std::shared_ptr<monty::ndarray<double,1>> & a, std::shared_ptr<monty::ndarray<double,1>> & w);
      static void env_potrf (int n, std::shared_ptr<monty::ndarray<double,1>> & a);
      static void env_syevd (int n, std::shared_ptr<monty::ndarray<double,1>> & a, std::shared_ptr<monty::ndarray<double,1>> & w);
    
      static void env_putlicensecode(std::shared_ptr<monty::ndarray<int,1>> code);
      static void env_putlicensepath(const std::string &licfile);
      static void env_putlicensewait(int wait);
    
      static std::string env_getversion();
    
      // void convertSolutionStatus(MSKsoltypee soltype, p_SolutionStruct * sol, MSKsolstae status, MSKprostae prosta);
    
      int64_t task_append_afes (int64_t n);
      void task_putafeflist  (array_t<int64_t> idxs, array_t<int> ptr, array_t<int>subj, array_t<double>cof, array_t<double>g);
      void task_putafebarfrowlist (array_t<int> idxs, array_t<int> ptr, array_t<int> barsubj, array_t<int64_t> symmatidx);
      void task_putafefijlist (array_t<int> &idxs, array_t<int> &subj, array_t<double> &cof);
      void task_putafefglist (array_t<int64_t> idxs, array_t<double> g);
      void task_clearafelist (array_t<int64_t>idxs);
      void task_putacclist  (array_t<int64_t>idxs, array_t<int64_t>domidxs, array_t<int64_t>afeidxs_t,array_t<double>g);
      void task_append_accs ( int64_t domidx, int numcone,array_t<int64_t> afeidxs,array_t<double> g);
    
      int64_t task_append_domain_quad     (int conesize);
      int64_t task_append_domain_rquad    (int conesize);
      int64_t task_append_domain_pexp     ();
      int64_t task_append_domain_ppow     (int conesize,array_t<double> alpha);
      int64_t task_append_domain_dexp     ();
      int64_t task_append_domain_dpow     (int conesize,array_t<double> alpha);
      /* int64_t task_append_domain_onenorm  (int conesize); */
      /* int64_t task_append_domain_infnorm  (int conesize); */
      int64_t task_append_domain_pgeomean (int conesize);
      int64_t task_append_domain_dgeomean (int conesize);
      int64_t task_append_domain_rpos     (int conesize);
      int64_t task_append_domain_rneg     (int conesize);
      int64_t task_append_domain_r        (int conesize);
      int64_t task_append_domain_rzero    (int conesize);
      int64_t task_append_domain_svec_psd (int conesize);
      int64_t task_append_domain_empty    ();
      int64_t task_append_djc             (int64_t n);
      void task_putdjcslice(int64_t first, int64_t last , array_t<int64_t> numterm_t, array_t<int64_t> termsizes, array_t<int64_t> domidxlist, array_t<int64_t> afeidxlist,  array_t<double> b);
    
    };
    // End of file 'src/fusion/cxx/BaseModel_p.h'
    struct p_Model : public ::mosek::fusion::p_BaseModel
    {
      Model * _pubthis;
      static mosek::fusion::p_Model* _get_impl(mosek::fusion::Model * _inst){ return static_cast< mosek::fusion::p_Model* >(mosek::fusion::p_BaseModel::_get_impl(_inst)); }
      static mosek::fusion::p_Model * _get_impl(mosek::fusion::Model::t _inst) { return _get_impl(_inst.get()); }
      p_Model(Model * _pubthis);
      virtual ~p_Model() { /* std::cout << "~p_Model" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::WorkStack > xs{};
      monty::rc_ptr< ::mosek::fusion::WorkStack > ws{};
      monty::rc_ptr< ::mosek::fusion::WorkStack > rs{};
      monty::rc_ptr< ::mosek::fusion::SolutionStruct > sol_itg{};
      monty::rc_ptr< ::mosek::fusion::SolutionStruct > sol_bas{};
      monty::rc_ptr< ::mosek::fusion::SolutionStruct > sol_itr{};
      monty::rc_ptr< ::mosek::fusion::Utils::StringIntMap > con_map{};
      std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::ModelConstraint >,1 > > acons{};
      std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::ModelConstraint >,1 > > cons{};
      int64_t task_numaferow{};
      std::shared_ptr< monty::ndarray< double,1 > > param_value{};
      int32_t param_num{};
      monty::rc_ptr< ::mosek::fusion::Utils::StringIntMap > par_map{};
      int32_t numparameter{};
      std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Parameter >,1 > > parameters{};
      std::shared_ptr< monty::ndarray< bool,1 > > initsol_xx_flag{};
      std::shared_ptr< monty::ndarray< double,1 > > initsol_xx{};
      monty::rc_ptr< ::mosek::fusion::Utils::StringIntMap > var_map{};
      std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::ModelVariable >,1 > > barvars{};
      std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::ModelVariable >,1 > > vars{};
      int32_t bfixidx{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > barvar_block_elm_j{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > barvar_block_elm_i{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > barvar_block_elm_barj{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > barvar_block_elm_ptr{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > barvar_block_dim{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > barvar_block_ptr{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > barvar_dim{};
      int32_t barvar_num{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > var_elm_acc_ofs{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > var_elm_acc_idx{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > var_block_acc_id{};
      monty::rc_ptr< ::mosek::fusion::LinkedBlocks > var_block_map{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > acon_elm_afe{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > acon_elm_ofs{};
      std::shared_ptr< monty::ndarray< double,1 > > acon_elm_scale{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > acon_elm_accid{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > acon_afe{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > acon_acc{};
      monty::rc_ptr< ::mosek::fusion::LinkedBlocks > acon_block_map{};
      monty::rc_ptr< ::mosek::fusion::LinkedBlocks > acc_block_map{};
      monty::rc_ptr< ::mosek::fusion::RowBlockManager > obj_blocks{};
      monty::rc_ptr< ::mosek::fusion::RowBlockManager > afe_blocks{};
      monty::rc_ptr< ::mosek::fusion::RowBlockManager > con_blocks{};
      int32_t num_task_acc{};
      int32_t num_task_afe{};
      int32_t num_task_con{};
      mosek::fusion::SolutionType solutionptr{};
      mosek::fusion::AccSolutionStatus acceptable_sol{};
      std::string model_name{};

      virtual void destroy();

      static Model::t _new_Model(monty::rc_ptr< ::mosek::fusion::Model > _697_m);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Model > _697_m);
      static Model::t _new_Model(const std::string &  _704_name,int32_t _705_basesize);
      void _initialize(const std::string &  _704_name,int32_t _705_basesize);
      static Model::t _new_Model(int32_t _713_basesize);
      void _initialize(int32_t _713_basesize);
      static Model::t _new_Model(const std::string &  _714_name);
      void _initialize(const std::string &  _714_name);
      static Model::t _new_Model();
      void _initialize();
      virtual monty::rc_ptr< ::mosek::fusion::Disjunction > __mosek_2fusion_2Model__disjunction(const std::string &  _715_name,std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Term >,1 > > _716_terms) ;
      virtual monty::rc_ptr< ::mosek::fusion::Disjunction > __mosek_2fusion_2Model__disjunction(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Term >,1 > > _762_terms) ;
      virtual monty::rc_ptr< ::mosek::fusion::Disjunction > __mosek_2fusion_2Model__disjunction(monty::rc_ptr< ::mosek::fusion::DisjunctionTerms > _763_djcterms) ;
      virtual monty::rc_ptr< ::mosek::fusion::Disjunction > __mosek_2fusion_2Model__disjunction(const std::string &  _764_name,monty::rc_ptr< ::mosek::fusion::DisjunctionTerms > _765_djcterms) ;
      virtual monty::rc_ptr< ::mosek::fusion::Disjunction > __mosek_2fusion_2Model__disjunction(const std::string &  _766_name,std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::ExprDomain >,1 > > _767_terms) ;
      virtual monty::rc_ptr< ::mosek::fusion::Disjunction > __mosek_2fusion_2Model__disjunction(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::ExprDomain >,1 > > _769_terms) ;
      virtual monty::rc_ptr< ::mosek::fusion::Disjunction > __mosek_2fusion_2Model__disjunction(const std::string &  _771_name,monty::rc_ptr< ::mosek::fusion::ExprDomain > _772_term) ;
      virtual monty::rc_ptr< ::mosek::fusion::Disjunction > __mosek_2fusion_2Model__disjunction(monty::rc_ptr< ::mosek::fusion::ExprDomain > _773_term) ;
      virtual monty::rc_ptr< ::mosek::fusion::Disjunction > __mosek_2fusion_2Model__disjunction(monty::rc_ptr< ::mosek::fusion::Term > _774_t1,monty::rc_ptr< ::mosek::fusion::Term > _775_t2,monty::rc_ptr< ::mosek::fusion::Term > _776_t3) ;
      virtual monty::rc_ptr< ::mosek::fusion::Disjunction > __mosek_2fusion_2Model__disjunction(monty::rc_ptr< ::mosek::fusion::Term > _777_t1,monty::rc_ptr< ::mosek::fusion::Term > _778_t2) ;
      virtual monty::rc_ptr< ::mosek::fusion::Disjunction > __mosek_2fusion_2Model__disjunction(monty::rc_ptr< ::mosek::fusion::Term > _779_t1) ;
      virtual monty::rc_ptr< ::mosek::fusion::Disjunction > __mosek_2fusion_2Model__disjunction(const std::string &  _780_name,monty::rc_ptr< ::mosek::fusion::Term > _781_t1,monty::rc_ptr< ::mosek::fusion::Term > _782_t2,monty::rc_ptr< ::mosek::fusion::Term > _783_t3) ;
      virtual monty::rc_ptr< ::mosek::fusion::Disjunction > __mosek_2fusion_2Model__disjunction(const std::string &  _784_name,monty::rc_ptr< ::mosek::fusion::Term > _785_t1,monty::rc_ptr< ::mosek::fusion::Term > _786_t2) ;
      virtual monty::rc_ptr< ::mosek::fusion::Disjunction > __mosek_2fusion_2Model__disjunction(const std::string &  _787_name,monty::rc_ptr< ::mosek::fusion::Term > _788_t1) ;
      virtual monty::rc_ptr< ::mosek::fusion::Utils::StringBuffer > __mosek_2fusion_2Model__formstConstr(monty::rc_ptr< ::mosek::fusion::Utils::StringBuffer > _789_sb,std::shared_ptr< monty::ndarray< int32_t,1 > > _790_shape,std::shared_ptr< monty::ndarray< int32_t,1 > > _791_idxs) ;
      virtual void acon_release(int32_t _792_id) ;
      virtual int32_t acon_allocate(int64_t _800_domidx,int32_t _801_conesize,int32_t _802_numcone,std::shared_ptr< monty::ndarray< double,1 > > _803_g,std::shared_ptr< monty::ndarray< int32_t,1 > > _804_afeidxs,std::shared_ptr< monty::ndarray< int32_t,1 > > _805_accidxs) ;
      virtual void afe_release(int32_t _831_id) ;
      virtual int32_t afe_allocate(std::shared_ptr< monty::ndarray< int32_t,1 > > _834_nativeidxs) ;
      virtual void con_release(int32_t _840_id) ;
      virtual int32_t con_allocate(std::shared_ptr< monty::ndarray< int32_t,1 > > _843_nativeidxs) ;
      virtual int32_t barvar_alloc(int32_t _851_conedim,int32_t _852_numcone,std::shared_ptr< monty::ndarray< int32_t,1 > > _853_barvaridxs,std::shared_ptr< monty::ndarray< int64_t,1 > > _854_nativeidxs) ;
      virtual int32_t conicvar_alloc(int64_t _885_domidx,int32_t _886_conesize,int32_t _887_numcone,std::shared_ptr< monty::ndarray< int32_t,1 > > _888_accidxs,std::shared_ptr< monty::ndarray< int32_t,1 > > _889_nativeidxs) ;
      virtual int32_t linearvar_alloc(int32_t _901_n,std::shared_ptr< monty::ndarray< int32_t,1 > > _902_nativeidxs) ;
      virtual void make_continuous(std::shared_ptr< monty::ndarray< int64_t,1 > > _914_idxs) ;
      virtual void make_integer(std::shared_ptr< monty::ndarray< int64_t,1 > > _920_idxs) ;
      static  void putlicensewait(bool _926_wait);
      static  void putlicensepath(const std::string &  _927_licfile);
      static  void putlicensecode(std::shared_ptr< monty::ndarray< int32_t,1 > > _928_code);
      virtual /* override */ void dispose() ;
      virtual MSKtask_t __mosek_2fusion_2Model__getTask() ;
      virtual void getConstraintDuals(bool _934_lower,std::shared_ptr< monty::ndarray< int32_t,1 > > _935_nativeidxs,std::shared_ptr< monty::ndarray< double,1 > > _936_res,int32_t _937_offset) ;
      virtual void getConstraintValues(bool _942_primal,std::shared_ptr< monty::ndarray< int32_t,1 > > _943_nativeidxs,std::shared_ptr< monty::ndarray< double,1 > > _944_res,int32_t _945_offset) ;
      virtual void getVariableDuals(bool _957_lower,std::shared_ptr< monty::ndarray< int64_t,1 > > _958_nativeidxs,std::shared_ptr< monty::ndarray< double,1 > > _959_res,int32_t _960_offset) ;
      virtual void getVariableValues(bool _966_primal,std::shared_ptr< monty::ndarray< int64_t,1 > > _967_nativeidxs,std::shared_ptr< monty::ndarray< double,1 > > _968_res,int32_t _969_offset) ;
      virtual void setVariableValues(bool _981_primal,std::shared_ptr< monty::ndarray< int64_t,1 > > _982_nativeidxs,std::shared_ptr< monty::ndarray< double,1 > > _983_values) ;
      virtual void writeTaskNoFlush(const std::string &  _990_filename) ;
      virtual void writeTaskStream(const std::string &  _991_ext,std::ostream&  _992_stream) ;
      virtual void dataReport() ;
      virtual void dataReport(int32_t _993_detail) ;
      virtual void writeTask(const std::string &  _994_filename) ;
      virtual int64_t getSolverLIntInfo(const std::string &  _995_name) ;
      virtual int32_t getSolverIntInfo(const std::string &  _996_name) ;
      virtual double getSolverDoubleInfo(const std::string &  _997_name) ;
      virtual void setCallbackHandler(mosek::cbhandler_t  _998_h) ;
      virtual void setDataCallbackHandler(mosek::datacbhandler_t  _999_h) ;
      virtual void setLogHandler(mosek::msghandler_t  _1000_h) ;
      virtual void setSolverParam(const std::string &  _1001_name,double _1002_floatval) ;
      virtual void setSolverParam(const std::string &  _1003_name,int32_t _1004_intval) ;
      virtual void setSolverParam(const std::string &  _1005_name,const std::string &  _1006_strval) ;
      virtual void breakSolver() ;
      virtual void optserverHost(const std::string &  _1007_addr) ;
      virtual /* override */ void report_solution(mosek::fusion::SolutionType _1008_soltype,mosek::fusion::ProblemStatus _1009_prosta,mosek::fusion::SolutionStatus _1010_psolsta,mosek::fusion::SolutionStatus _1011_dsolsta,double _1012_pobj,double _1013_dobj,int32_t _1014_numvar,int32_t _1015_numcon,int32_t _1016_numbarelm,int32_t _1017_numacc,int32_t _1018_numaccelm,bool _1019_hasprimal,bool _1020_hasdual) ;
      virtual /* override */ void clear_solutions() ;
      static  std::shared_ptr< monty::ndarray< mosek::fusion::SolverStatus,1 > > solveBatch(bool _1030_israce,double _1031_maxtime,int32_t _1032_numthreads,std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Model >,1 > > _1033_models);
      virtual void solve(const std::string &  _1039_addr,const std::string &  _1040_accesstoken) ;
      virtual void solve() ;
      virtual void flush_parameters() ;
      virtual void flushParameters() ;
      virtual void evaluate_parameterized(monty::rc_ptr< ::mosek::fusion::WorkStack > _1053_xs,int32_t _1054_numrow,std::shared_ptr< monty::ndarray< int32_t,1 > > _1055_rowptrb,std::shared_ptr< monty::ndarray< int32_t,1 > > _1056_rowptre,std::shared_ptr< monty::ndarray< int64_t,1 > > _1057_codenidx,std::shared_ptr< monty::ndarray< int32_t,1 > > _1058_codeptr,std::shared_ptr< monty::ndarray< int32_t,1 > > _1059_codesizes,std::shared_ptr< monty::ndarray< int32_t,1 > > _1060_code,std::shared_ptr< monty::ndarray< double,1 > > _1061_cconst,std::shared_ptr< monty::ndarray< int32_t,1 > > _1062_subj,std::shared_ptr< monty::ndarray< double,1 > > _1063_val) ;
      virtual void flushSolutions() ;
      virtual void flush_initsol(mosek::fusion::SolutionType _1074_which) ;
      virtual mosek::fusion::SolutionStatus getDualSolutionStatus() ;
      virtual mosek::fusion::ProblemStatus getProblemStatus() ;
      virtual mosek::fusion::SolutionStatus getPrimalSolutionStatus() ;
      virtual double dualObjValue() ;
      virtual double primalObjValue() ;
      virtual monty::rc_ptr< ::mosek::fusion::SolutionStruct > __mosek_2fusion_2Model__get_sol_cache(mosek::fusion::SolutionType _1081_which_,bool _1082_primal,bool _1083_nothrow) ;
      virtual monty::rc_ptr< ::mosek::fusion::SolutionStruct > __mosek_2fusion_2Model__get_sol_cache(mosek::fusion::SolutionType _1089_which_,bool _1090_primal) ;
      virtual void setSolution_xx(std::shared_ptr< monty::ndarray< int32_t,1 > > _1091_subj,std::shared_ptr< monty::ndarray< double,1 > > _1092_val) ;
      virtual void ensure_initsol_xx() ;
      virtual std::shared_ptr< monty::ndarray< int32_t,1 > > getSolution_accptr(mosek::fusion::SolutionType _1099_which) ;
      virtual std::shared_ptr< monty::ndarray< double,1 > > getSolution_accy(mosek::fusion::SolutionType _1100_which) ;
      virtual std::shared_ptr< monty::ndarray< double,1 > > getSolution_accx(mosek::fusion::SolutionType _1101_which) ;
      virtual std::shared_ptr< monty::ndarray< double,1 > > getSolution_bars(mosek::fusion::SolutionType _1102_which) ;
      virtual std::shared_ptr< monty::ndarray< double,1 > > getSolution_barx(mosek::fusion::SolutionType _1103_which) ;
      virtual std::shared_ptr< monty::ndarray< double,1 > > getSolution_y(mosek::fusion::SolutionType _1104_which) ;
      virtual std::shared_ptr< monty::ndarray< double,1 > > getSolution_xc(mosek::fusion::SolutionType _1105_which) ;
      virtual std::shared_ptr< monty::ndarray< double,1 > > getSolution_suc(mosek::fusion::SolutionType _1106_which) ;
      virtual std::shared_ptr< monty::ndarray< double,1 > > getSolution_slc(mosek::fusion::SolutionType _1107_which) ;
      virtual std::shared_ptr< monty::ndarray< double,1 > > getSolution_sux(mosek::fusion::SolutionType _1108_which) ;
      virtual std::shared_ptr< monty::ndarray< double,1 > > getSolution_slx(mosek::fusion::SolutionType _1109_which) ;
      virtual std::shared_ptr< monty::ndarray< double,1 > > getSolution_yx(mosek::fusion::SolutionType _1110_which) ;
      virtual std::shared_ptr< monty::ndarray< double,1 > > getSolution_xx(mosek::fusion::SolutionType _1111_which) ;
      virtual void selectedSolution(mosek::fusion::SolutionType _1112_soltype) ;
      virtual mosek::fusion::AccSolutionStatus getAcceptedSolutionStatus() ;
      virtual void acceptedSolutionStatus(mosek::fusion::AccSolutionStatus _1113_what) ;
      virtual mosek::fusion::ProblemStatus getProblemStatus(mosek::fusion::SolutionType _1114_which) ;
      virtual mosek::fusion::SolutionStatus getDualSolutionStatus(mosek::fusion::SolutionType _1116_which) ;
      virtual mosek::fusion::SolutionStatus getPrimalSolutionStatus(mosek::fusion::SolutionType _1117_which) ;
      virtual mosek::fusion::SolutionStatus getSolutionStatus(mosek::fusion::SolutionType _1118_which,bool _1119_primal) ;
      virtual void update(std::shared_ptr< monty::ndarray< int32_t,1 > > _1122_conidxs,monty::rc_ptr< ::mosek::fusion::Expression > _1123_expr) ;
      virtual void update(std::shared_ptr< monty::ndarray< int32_t,1 > > _1190_conidxs,monty::rc_ptr< ::mosek::fusion::Expression > _1191_expr,std::shared_ptr< monty::ndarray< int32_t,1 > > _1192_varidxs) ;
      virtual void updateObjective(monty::rc_ptr< ::mosek::fusion::Expression > _1294_expr,monty::rc_ptr< ::mosek::fusion::Variable > _1295_x) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Model__parameter_unchecked(const std::string &  _1332_name,std::shared_ptr< monty::ndarray< int32_t,1 > > _1333_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _1334_sp) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Model__parameter_(const std::string &  _1344_name,std::shared_ptr< monty::ndarray< int32_t,1 > > _1345_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _1346_sp) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Model__parameter_(const std::string &  _1351_name,std::shared_ptr< monty::ndarray< int32_t,1 > > _1352_shape,std::shared_ptr< monty::ndarray< int32_t,2 > > _1353_sparsity) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Model__parameter(const std::string &  _1361_name) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Model__parameter(const std::string &  _1363_name,int32_t _1364_d1,int32_t _1365_d2,int32_t _1366_d3) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Model__parameter(const std::string &  _1368_name,int32_t _1369_d1,int32_t _1370_d2) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Model__parameter(const std::string &  _1372_name,int32_t _1373_d1) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Model__parameter(const std::string &  _1375_name,std::shared_ptr< monty::ndarray< int32_t,1 > > _1376_shape) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Model__parameter(const std::string &  _1378_name,std::shared_ptr< monty::ndarray< int32_t,1 > > _1379_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _1380_sp) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Model__parameter(const std::string &  _1381_name,std::shared_ptr< monty::ndarray< int32_t,1 > > _1382_shape,std::shared_ptr< monty::ndarray< int32_t,2 > > _1383_sparsity) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Model__parameter() ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Model__parameter(int32_t _1385_d1,int32_t _1386_d2,int32_t _1387_d3) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Model__parameter(int32_t _1389_d1,int32_t _1390_d2) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Model__parameter(int32_t _1392_d1) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Model__parameter(std::shared_ptr< monty::ndarray< int32_t,1 > > _1394_shape) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Model__parameter(std::shared_ptr< monty::ndarray< int32_t,1 > > _1396_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _1397_sp) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Model__parameter(std::shared_ptr< monty::ndarray< int32_t,1 > > _1398_shape,std::shared_ptr< monty::ndarray< int32_t,2 > > _1399_sparsity) ;
      virtual void objective_(const std::string &  _1400_name,mosek::fusion::ObjectiveSense _1401_sense,monty::rc_ptr< ::mosek::fusion::Expression > _1402_expr) ;
      virtual void objective(double _1436_c) ;
      virtual void objective(mosek::fusion::ObjectiveSense _1437_sense,double _1438_c) ;
      virtual void objective(mosek::fusion::ObjectiveSense _1439_sense,monty::rc_ptr< ::mosek::fusion::Expression > _1440_expr) ;
      virtual void objective(const std::string &  _1441_name,double _1442_c) ;
      virtual void objective(const std::string &  _1443_name,mosek::fusion::ObjectiveSense _1444_sense,double _1445_c) ;
      virtual void objective(const std::string &  _1446_name,mosek::fusion::ObjectiveSense _1447_sense,monty::rc_ptr< ::mosek::fusion::Expression > _1448_expr) ;
      virtual monty::rc_ptr< ::mosek::fusion::RangedConstraint > __mosek_2fusion_2Model__constraint(monty::rc_ptr< ::mosek::fusion::ExprRangeDomain > _1449_exprdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::RangedConstraint > __mosek_2fusion_2Model__constraint(const std::string &  _1450_name,monty::rc_ptr< ::mosek::fusion::ExprRangeDomain > _1451_exprdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::RangedConstraint > __mosek_2fusion_2Model__constraint(monty::rc_ptr< ::mosek::fusion::Expression > _1452_expr,monty::rc_ptr< ::mosek::fusion::RangeDomain > _1453_rdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::RangedConstraint > __mosek_2fusion_2Model__constraint(const std::string &  _1454_name,monty::rc_ptr< ::mosek::fusion::Expression > _1455_expr,monty::rc_ptr< ::mosek::fusion::RangeDomain > _1456_rdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Model__constraint(monty::rc_ptr< ::mosek::fusion::ExprConicDomain > _1457_exprdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Model__constraint(const std::string &  _1458_name,monty::rc_ptr< ::mosek::fusion::ExprConicDomain > _1459_exprdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Model__constraint(monty::rc_ptr< ::mosek::fusion::Expression > _1460_expr,monty::rc_ptr< ::mosek::fusion::ConeDomain > _1461_qdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Model__constraint(const std::string &  _1462_name,monty::rc_ptr< ::mosek::fusion::Expression > _1463_expr,monty::rc_ptr< ::mosek::fusion::ConeDomain > _1464_qdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Model__constraint(monty::rc_ptr< ::mosek::fusion::ExprLinearDomain > _1465_exprdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Model__constraint(const std::string &  _1466_name,monty::rc_ptr< ::mosek::fusion::ExprLinearDomain > _1467_exprdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Model__constraint(monty::rc_ptr< ::mosek::fusion::Expression > _1468_expr,monty::rc_ptr< ::mosek::fusion::LinearDomain > _1469_ldom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Model__constraint(const std::string &  _1470_name,monty::rc_ptr< ::mosek::fusion::Expression > _1471_expr,monty::rc_ptr< ::mosek::fusion::LinearDomain > _1472_ldom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Model__constraint(monty::rc_ptr< ::mosek::fusion::ExprPSDDomain > _1473_exprdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Model__constraint(const std::string &  _1474_name,monty::rc_ptr< ::mosek::fusion::ExprPSDDomain > _1475_exprdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Model__constraint(monty::rc_ptr< ::mosek::fusion::Expression > _1476_expr,monty::rc_ptr< ::mosek::fusion::PSDDomain > _1477_psddom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Model__constraint(const std::string &  _1478_name,monty::rc_ptr< ::mosek::fusion::Expression > _1479_expr,monty::rc_ptr< ::mosek::fusion::PSDDomain > _1480_psddom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(monty::rc_ptr< ::mosek::fusion::PSDDomain > _1481_psddom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(int32_t _1482_n,int32_t _1483_m,monty::rc_ptr< ::mosek::fusion::PSDDomain > _1484_psddom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(int32_t _1485_n,monty::rc_ptr< ::mosek::fusion::PSDDomain > _1486_psddom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(const std::string &  _1487_name,monty::rc_ptr< ::mosek::fusion::PSDDomain > _1488_psddom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(const std::string &  _1489_name,int32_t _1490_n,int32_t _1491_m,monty::rc_ptr< ::mosek::fusion::PSDDomain > _1492_psddom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(const std::string &  _1493_name,int32_t _1494_n,monty::rc_ptr< ::mosek::fusion::PSDDomain > _1495_psddom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(const std::string &  _1496_name,std::shared_ptr< monty::ndarray< int32_t,1 > > _1497_shp,monty::rc_ptr< ::mosek::fusion::PSDDomain > _1498_psddom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(monty::rc_ptr< ::mosek::fusion::ConeDomain > _1499_qdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::RangedVariable > __mosek_2fusion_2Model__variable(monty::rc_ptr< ::mosek::fusion::RangeDomain > _1500_rdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(monty::rc_ptr< ::mosek::fusion::LinearDomain > _1501_ldom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(std::shared_ptr< monty::ndarray< int32_t,1 > > _1502_shp,monty::rc_ptr< ::mosek::fusion::ConeDomain > _1503_qdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::RangedVariable > __mosek_2fusion_2Model__variable(std::shared_ptr< monty::ndarray< int32_t,1 > > _1504_shp,monty::rc_ptr< ::mosek::fusion::RangeDomain > _1505_rdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(std::shared_ptr< monty::ndarray< int32_t,1 > > _1506_shp,monty::rc_ptr< ::mosek::fusion::LinearDomain > _1507_ldom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(std::shared_ptr< monty::ndarray< int32_t,1 > > _1508_shp) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(int32_t _1509_size,monty::rc_ptr< ::mosek::fusion::ConeDomain > _1510_qdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::RangedVariable > __mosek_2fusion_2Model__variable(int32_t _1511_size,monty::rc_ptr< ::mosek::fusion::RangeDomain > _1512_rdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(int32_t _1513_size,monty::rc_ptr< ::mosek::fusion::LinearDomain > _1514_ldom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(int32_t _1515_size) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable() ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(const std::string &  _1516_name,monty::rc_ptr< ::mosek::fusion::ConeDomain > _1517_qdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::RangedVariable > __mosek_2fusion_2Model__variable(const std::string &  _1518_name,monty::rc_ptr< ::mosek::fusion::RangeDomain > _1519_rdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(const std::string &  _1520_name,monty::rc_ptr< ::mosek::fusion::LinearDomain > _1521_ldom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(const std::string &  _1522_name,std::shared_ptr< monty::ndarray< int32_t,1 > > _1523_shp,monty::rc_ptr< ::mosek::fusion::ConeDomain > _1524_qdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::RangedVariable > __mosek_2fusion_2Model__variable(const std::string &  _1525_name,std::shared_ptr< monty::ndarray< int32_t,1 > > _1526_shp,monty::rc_ptr< ::mosek::fusion::RangeDomain > _1527_rdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(const std::string &  _1528_name,std::shared_ptr< monty::ndarray< int32_t,1 > > _1529_shp,monty::rc_ptr< ::mosek::fusion::LinearDomain > _1530_ldom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(const std::string &  _1531_name,std::shared_ptr< monty::ndarray< int32_t,1 > > _1532_shp) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(const std::string &  _1533_name,int32_t _1534_size,monty::rc_ptr< ::mosek::fusion::ConeDomain > _1535_qdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::RangedVariable > __mosek_2fusion_2Model__variable(const std::string &  _1536_name,int32_t _1537_size,monty::rc_ptr< ::mosek::fusion::RangeDomain > _1538_rdom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(const std::string &  _1539_name,int32_t _1540_size,monty::rc_ptr< ::mosek::fusion::LinearDomain > _1541_ldom) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(const std::string &  _1542_name,int32_t _1543_size) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable(const std::string &  _1544_name) ;
      virtual void removeConstraintBlock(int32_t _1545_conid) ;
      virtual void removeVariableBlock(int64_t _1546_varid64) ;
      virtual monty::rc_ptr< ::mosek::fusion::RangedVariable > __mosek_2fusion_2Model__ranged_variable(const std::string &  _1551_name,std::shared_ptr< monty::ndarray< int32_t,1 > > _1552_shp,monty::rc_ptr< ::mosek::fusion::RangeDomain > _1553_dom_) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable_(const std::string &  _1614_name,std::shared_ptr< monty::ndarray< int32_t,1 > > _1615_shp,monty::rc_ptr< ::mosek::fusion::ConeDomain > _1616_dom_) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable_(const std::string &  _1654_name,std::shared_ptr< monty::ndarray< int32_t,1 > > _1655_shp,monty::rc_ptr< ::mosek::fusion::LinearDomain > _1656_dom_) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__variable_(const std::string &  _1715_name,std::shared_ptr< monty::ndarray< int32_t,1 > > _1716_shp,monty::rc_ptr< ::mosek::fusion::PSDDomain > _1717_dom_) ;
      virtual void putfrows(std::shared_ptr< monty::ndarray< int32_t,1 > > _1746_nativeidxs,int32_t _1747_nativebaseptr,monty::rc_ptr< ::mosek::fusion::WorkStack > _1748_rs,int32_t _1749_nelem,int32_t _1750_nnz,int32_t _1751_ptr_base,int32_t _1752_nidxs_base,int32_t _1753_cof_base) ;
      virtual void putarows(std::shared_ptr< monty::ndarray< int32_t,1 > > _1793_nativeidxs,monty::rc_ptr< ::mosek::fusion::WorkStack > _1794_rs,int32_t _1795_nelem,int32_t _1796_nnz,int32_t _1797_ptr_base,int32_t _1798_nidxs_base,int32_t _1799_cof_base) ;
      virtual monty::rc_ptr< ::mosek::fusion::RangedConstraint > __mosek_2fusion_2Model__constraint_(const std::string &  _1836_name,monty::rc_ptr< ::mosek::fusion::Expression > _1837_expr,monty::rc_ptr< ::mosek::fusion::RangeDomain > _1838_dom_) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Model__constraint_(const std::string &  _1881_name,monty::rc_ptr< ::mosek::fusion::Expression > _1882_expr,monty::rc_ptr< ::mosek::fusion::PSDDomain > _1883_dom_) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Model__constraint_(const std::string &  _1977_name,monty::rc_ptr< ::mosek::fusion::Expression > _1978_expr,monty::rc_ptr< ::mosek::fusion::ConeDomain > _1979_dom_) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Model__constraint_(const std::string &  _2036_name,monty::rc_ptr< ::mosek::fusion::Expression > _2037_expr,monty::rc_ptr< ::mosek::fusion::LinearDomain > _2038_dom_) ;
      static  std::string getVersion();
      virtual bool hasParameter(const std::string &  _2077_name) ;
      virtual bool hasConstraint(const std::string &  _2078_name) ;
      virtual bool hasVariable(const std::string &  _2079_name) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Model__getParameter(const std::string &  _2080_name) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Model__getConstraint(int32_t _2081_index) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Model__getConstraint(const std::string &  _2083_name) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__getVariable(int32_t _2086_index) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Model__getVariable(const std::string &  _2087_name) ;
      virtual std::string getName() ;
      virtual std::shared_ptr< monty::ndarray< double,1 > > getParameterValue(std::shared_ptr< monty::ndarray< int32_t,1 > > _2089_idxs) ;
      virtual void setParameterValue(std::shared_ptr< monty::ndarray< int32_t,1 > > _2092_idxs,std::shared_ptr< monty::ndarray< double,1 > > _2093_vals) ;
      virtual monty::rc_ptr< ::mosek::fusion::Model > __mosek_2fusion_2Model__clone() ;
    }; // struct Model;

    // mosek.fusion.Debug from file 'src/fusion/cxx/Debug_p.h'
    // namespace mosek::fusion
    struct p_Debug
    {
      Debug * _pubthis;
    
      p_Debug(Debug * _pubthis) : _pubthis(_pubthis) {}
    
      static Debug::t o ()                 { return monty::rc_ptr<Debug>(new Debug()); }
      Debug::t p (const std::string & val) { std::cout << val; return Debug::t(_pubthis); }
      Debug::t p (      int val)           { std::cout << val; return Debug::t(_pubthis); }
      Debug::t p (int64_t val)           { std::cout << val; return Debug::t(_pubthis); }
      Debug::t p (   double val)           { std::cout << val; return Debug::t(_pubthis); }
      Debug::t p (     bool val)           { std::cout << val; return Debug::t(_pubthis); }
      Debug::t lf()                        { std::cout << std::endl; return Debug::t(_pubthis); }
    
      static p_Debug * _get_impl(Debug * _inst) { return _inst->ptr; }
    
      template<typename T>
      Debug::t p(const std::shared_ptr<monty::ndarray<T,1>> & val)
      {
        if (val->size() > 0)
        {
          std::cout << (*val)[0];
          for (int i = 1; i < val->size(); ++i)
            std::cout << "," << (*val)[i];
        }
        return Debug::t(_pubthis);
      }
    
      Debug::t __mosek_2fusion_2Debug__p (const std::string & val) { std::cout << val; return Debug::t(_pubthis); }
    
      template<class C>
      Debug::t __mosek_2fusion_2Debug__p (C val) { std::cout << val; return Debug::t(_pubthis); }
      Debug::t __mosek_2fusion_2Debug__lf()      { std::cout << std::endl; return Debug::t(_pubthis); }
    
    };
    // End of file 'src/fusion/cxx/Debug_p.h'
    struct p_Sort
    {
      Sort * _pubthis;
      static mosek::fusion::p_Sort* _get_impl(mosek::fusion::Sort * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_Sort * _get_impl(mosek::fusion::Sort::t _inst) { return _get_impl(_inst.get()); }
      p_Sort(Sort * _pubthis);
      virtual ~p_Sort() { /* std::cout << "~p_Sort" << std::endl;*/ };

      virtual void destroy();

      static  void argTransposeSort(std::shared_ptr< monty::ndarray< int64_t,1 > > _303_perm,std::shared_ptr< monty::ndarray< int64_t,1 > > _304_ptrb,int32_t _305_m,int32_t _306_n,int32_t _307_p,std::shared_ptr< monty::ndarray< int64_t,1 > > _308_val);
      static  void argsort(std::shared_ptr< monty::ndarray< int64_t,1 > > _316_idx,std::shared_ptr< monty::ndarray< int64_t,1 > > _317_vals1);
      static  void argsort(std::shared_ptr< monty::ndarray< int64_t,1 > > _318_idx,std::shared_ptr< monty::ndarray< int32_t,1 > > _319_vals1);
      static  void argsort(std::shared_ptr< monty::ndarray< int64_t,1 > > _320_idx,std::shared_ptr< monty::ndarray< int64_t,1 > > _321_vals1,std::shared_ptr< monty::ndarray< int64_t,1 > > _322_vals2);
      static  void argsort(std::shared_ptr< monty::ndarray< int64_t,1 > > _323_idx,std::shared_ptr< monty::ndarray< int32_t,1 > > _324_vals1,std::shared_ptr< monty::ndarray< int32_t,1 > > _325_vals2);
      static  void argsort(std::shared_ptr< monty::ndarray< int64_t,1 > > _326_idx,std::shared_ptr< monty::ndarray< int64_t,1 > > _327_vals1,int64_t _328_first,int64_t _329_last);
      static  void argsort(std::shared_ptr< monty::ndarray< int64_t,1 > > _330_idx,std::shared_ptr< monty::ndarray< int32_t,1 > > _331_vals1,int64_t _332_first,int64_t _333_last);
      static  void argsort(std::shared_ptr< monty::ndarray< int64_t,1 > > _334_idx,std::shared_ptr< monty::ndarray< int64_t,1 > > _335_vals1,std::shared_ptr< monty::ndarray< int64_t,1 > > _336_vals2,int64_t _337_first,int64_t _338_last);
      static  void argsort(std::shared_ptr< monty::ndarray< int64_t,1 > > _339_idx,std::shared_ptr< monty::ndarray< int32_t,1 > > _340_vals1,std::shared_ptr< monty::ndarray< int32_t,1 > > _341_vals2,int64_t _342_first,int64_t _343_last);
      static  void argsort(std::shared_ptr< monty::ndarray< int64_t,1 > > _344_idx,std::shared_ptr< monty::ndarray< int64_t,1 > > _345_vals1,int64_t _346_first,int64_t _347_last,bool _348_check);
      static  void argsort(std::shared_ptr< monty::ndarray< int64_t,1 > > _351_idx,std::shared_ptr< monty::ndarray< int32_t,1 > > _352_vals1,int64_t _353_first,int64_t _354_last,bool _355_check);
      static  void argsort(std::shared_ptr< monty::ndarray< int64_t,1 > > _358_idx,std::shared_ptr< monty::ndarray< int64_t,1 > > _359_vals1,std::shared_ptr< monty::ndarray< int64_t,1 > > _360_vals2,int64_t _361_first,int64_t _362_last,bool _363_check);
      static  void argsort(std::shared_ptr< monty::ndarray< int64_t,1 > > _366_idx,std::shared_ptr< monty::ndarray< int32_t,1 > > _367_vals1,std::shared_ptr< monty::ndarray< int32_t,1 > > _368_vals2,int64_t _369_first,int64_t _370_last,bool _371_check);
      static  void argbucketsort(std::shared_ptr< monty::ndarray< int64_t,1 > > _374_idx,std::shared_ptr< monty::ndarray< int64_t,1 > > _375_vals,int64_t _376_first,int64_t _377_last,int64_t _378_minv,int64_t _379_maxv);
      static  void argbucketsort(std::shared_ptr< monty::ndarray< int64_t,1 > > _380_idx,std::shared_ptr< monty::ndarray< int32_t,1 > > _381_vals,int64_t _382_first,int64_t _383_last,int32_t _384_minv,int32_t _385_maxv);
      static  void getminmax(std::shared_ptr< monty::ndarray< int64_t,1 > > _386_idx,std::shared_ptr< monty::ndarray< int64_t,1 > > _387_vals1,std::shared_ptr< monty::ndarray< int64_t,1 > > _388_vals2,int64_t _389_first,int64_t _390_last,std::shared_ptr< monty::ndarray< int64_t,1 > > _391_res);
      static  void getminmax(std::shared_ptr< monty::ndarray< int64_t,1 > > _394_idx,std::shared_ptr< monty::ndarray< int32_t,1 > > _395_vals1,std::shared_ptr< monty::ndarray< int32_t,1 > > _396_vals2,int64_t _397_first,int64_t _398_last,std::shared_ptr< monty::ndarray< int32_t,1 > > _399_res);
      static  bool issorted(std::shared_ptr< monty::ndarray< int64_t,1 > > _402_idx,std::shared_ptr< monty::ndarray< int64_t,1 > > _403_vals1,int64_t _404_first,int64_t _405_last,bool _406_check);
      static  bool issorted(std::shared_ptr< monty::ndarray< int64_t,1 > > _408_idx,std::shared_ptr< monty::ndarray< int32_t,1 > > _409_vals1,int64_t _410_first,int64_t _411_last,bool _412_check);
      static  bool issorted(std::shared_ptr< monty::ndarray< int64_t,1 > > _414_idx,std::shared_ptr< monty::ndarray< int64_t,1 > > _415_vals1,std::shared_ptr< monty::ndarray< int64_t,1 > > _416_vals2,int64_t _417_first,int64_t _418_last,bool _419_check);
      static  bool issorted(std::shared_ptr< monty::ndarray< int64_t,1 > > _421_idx,std::shared_ptr< monty::ndarray< int32_t,1 > > _422_vals1,std::shared_ptr< monty::ndarray< int32_t,1 > > _423_vals2,int64_t _424_first,int64_t _425_last,bool _426_check);
    }; // struct Sort;

    struct p_IndexCounter
    {
      IndexCounter * _pubthis;
      static mosek::fusion::p_IndexCounter* _get_impl(mosek::fusion::IndexCounter * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_IndexCounter * _get_impl(mosek::fusion::IndexCounter::t _inst) { return _get_impl(_inst.get()); }
      p_IndexCounter(IndexCounter * _pubthis);
      virtual ~p_IndexCounter() { /* std::cout << "~p_IndexCounter" << std::endl;*/ };
      int64_t start{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > dims{};
      std::shared_ptr< monty::ndarray< int64_t,1 > > strides{};
      std::shared_ptr< monty::ndarray< int64_t,1 > > st{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > ii{};
      int32_t n{};

      virtual void destroy();

      static IndexCounter::t _new_IndexCounter(std::shared_ptr< monty::ndarray< int32_t,1 > > _428_shape);
      void _initialize(std::shared_ptr< monty::ndarray< int32_t,1 > > _428_shape);
      static IndexCounter::t _new_IndexCounter(int64_t _430_start_,std::shared_ptr< monty::ndarray< int32_t,1 > > _431_dims_,std::shared_ptr< monty::ndarray< int32_t,1 > > _432_shape);
      void _initialize(int64_t _430_start_,std::shared_ptr< monty::ndarray< int32_t,1 > > _431_dims_,std::shared_ptr< monty::ndarray< int32_t,1 > > _432_shape);
      static IndexCounter::t _new_IndexCounter(int64_t _435_start_,std::shared_ptr< monty::ndarray< int32_t,1 > > _436_dims_,std::shared_ptr< monty::ndarray< int64_t,1 > > _437_strides_);
      void _initialize(int64_t _435_start_,std::shared_ptr< monty::ndarray< int32_t,1 > > _436_dims_,std::shared_ptr< monty::ndarray< int64_t,1 > > _437_strides_);
      virtual bool atEnd() ;
      virtual std::shared_ptr< monty::ndarray< int32_t,1 > > getIndex() ;
      virtual int64_t next() ;
      virtual int64_t get() ;
      virtual void inc() ;
      virtual void reset() ;
    }; // struct IndexCounter;

    struct p_CommonTools
    {
      CommonTools * _pubthis;
      static mosek::fusion::p_CommonTools* _get_impl(mosek::fusion::CommonTools * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_CommonTools * _get_impl(mosek::fusion::CommonTools::t _inst) { return _get_impl(_inst.get()); }
      p_CommonTools(CommonTools * _pubthis);
      virtual ~p_CommonTools() { /* std::cout << "~p_CommonTools" << std::endl;*/ };

      virtual void destroy();

      static  std::shared_ptr< monty::ndarray< int64_t,1 > > resize(std::shared_ptr< monty::ndarray< int64_t,1 > > _443_values,int32_t _444_newsize);
      static  std::shared_ptr< monty::ndarray< int32_t,1 > > resize(std::shared_ptr< monty::ndarray< int32_t,1 > > _446_values,int32_t _447_newsize);
      static  std::shared_ptr< monty::ndarray< double,1 > > resize(std::shared_ptr< monty::ndarray< double,1 > > _449_values,int32_t _450_newsize);
      static  int32_t binarySearch(std::shared_ptr< monty::ndarray< int32_t,1 > > _452_values,int32_t _453_target);
      static  int32_t binarySearch(std::shared_ptr< monty::ndarray< int64_t,1 > > _457_values,int64_t _458_target);
      static  int32_t binarySearchR(std::shared_ptr< monty::ndarray< int64_t,1 > > _460_values,int64_t _461_target);
      static  int32_t binarySearchL(std::shared_ptr< monty::ndarray< int64_t,1 > > _465_values,int64_t _466_target);
      static  void ndIncr(std::shared_ptr< monty::ndarray< int32_t,1 > > _470_ndidx,std::shared_ptr< monty::ndarray< int32_t,1 > > _471_first,std::shared_ptr< monty::ndarray< int32_t,1 > > _472_last);
      static  void transposeTriplets(std::shared_ptr< monty::ndarray< int32_t,1 > > _474_subi,std::shared_ptr< monty::ndarray< int32_t,1 > > _475_subj,std::shared_ptr< monty::ndarray< double,1 > > _476_val,std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< int64_t,1 > >,1 > > _477_tsubi_,std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< int64_t,1 > >,1 > > _478_tsubj_,std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< double,1 > >,1 > > _479_tval_,int64_t _480_nelm,int32_t _481_dimi,int32_t _482_dimj);
      static  void transposeTriplets(std::shared_ptr< monty::ndarray< int32_t,1 > > _495_subi,std::shared_ptr< monty::ndarray< int32_t,1 > > _496_subj,std::shared_ptr< monty::ndarray< double,1 > > _497_val,std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< int32_t,1 > >,1 > > _498_tsubi_,std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< int32_t,1 > >,1 > > _499_tsubj_,std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< double,1 > >,1 > > _500_tval_,int64_t _501_nelm,int32_t _502_dimi,int32_t _503_dimj);
      static  void tripletSort(std::shared_ptr< monty::ndarray< int32_t,1 > > _516_subi,std::shared_ptr< monty::ndarray< int32_t,1 > > _517_subj,std::shared_ptr< monty::ndarray< double,1 > > _518_val,std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< int32_t,1 > >,1 > > _519_tsubi_,std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< int32_t,1 > >,1 > > _520_tsubj_,std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< double,1 > >,1 > > _521_tval_,int64_t _522_nelm,int32_t _523_dimi,int32_t _524_dimj);
      static  void argMSort(std::shared_ptr< monty::ndarray< int32_t,1 > > _550_idx,std::shared_ptr< monty::ndarray< int32_t,1 > > _551_vals);
      static  void mergeInto(std::shared_ptr< monty::ndarray< int32_t,1 > > _556_src,std::shared_ptr< monty::ndarray< int32_t,1 > > _557_tgt,std::shared_ptr< monty::ndarray< int32_t,1 > > _558_vals,int32_t _559_si0,int32_t _560_si1_,int32_t _561_si2_);
      static  void argQsort(std::shared_ptr< monty::ndarray< int64_t,1 > > _567_idx,std::shared_ptr< monty::ndarray< int64_t,1 > > _568_vals1,std::shared_ptr< monty::ndarray< int64_t,1 > > _569_vals2,int64_t _570_first,int64_t _571_last);
      static  void argQsort(std::shared_ptr< monty::ndarray< int64_t,1 > > _572_idx,std::shared_ptr< monty::ndarray< int32_t,1 > > _573_vals1,std::shared_ptr< monty::ndarray< int32_t,1 > > _574_vals2,int64_t _575_first,int64_t _576_last);
    }; // struct CommonTools;

    struct p_SolutionStruct
    {
      SolutionStruct * _pubthis;
      static mosek::fusion::p_SolutionStruct* _get_impl(mosek::fusion::SolutionStruct * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_SolutionStruct * _get_impl(mosek::fusion::SolutionStruct::t _inst) { return _get_impl(_inst.get()); }
      p_SolutionStruct(SolutionStruct * _pubthis);
      virtual ~p_SolutionStruct() { /* std::cout << "~p_SolutionStruct" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< double,1 > > accy{};
      std::shared_ptr< monty::ndarray< double,1 > > accx{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > accptr{};
      std::shared_ptr< monty::ndarray< double,1 > > yx{};
      std::shared_ptr< monty::ndarray< double,1 > > sux{};
      std::shared_ptr< monty::ndarray< double,1 > > slx{};
      std::shared_ptr< monty::ndarray< double,1 > > bars{};
      std::shared_ptr< monty::ndarray< double,1 > > barx{};
      std::shared_ptr< monty::ndarray< double,1 > > y{};
      std::shared_ptr< monty::ndarray< double,1 > > suc{};
      std::shared_ptr< monty::ndarray< double,1 > > slc{};
      std::shared_ptr< monty::ndarray< double,1 > > xx{};
      std::shared_ptr< monty::ndarray< double,1 > > xc{};
      double dobj{};
      double pobj{};
      mosek::fusion::ProblemStatus probstatus{};
      mosek::fusion::SolutionStatus dstatus{};
      mosek::fusion::SolutionStatus pstatus{};
      int32_t sol_numbarvar{};
      int32_t sol_numaccelm{};
      int32_t sol_numacc{};
      int32_t sol_numvar{};
      int32_t sol_numcon{};

      virtual void destroy();

      static SolutionStruct::t _new_SolutionStruct(int32_t _577_numvar,int32_t _578_numcon,int32_t _579_numbarvar,int32_t _580_numacc,int32_t _581_numaccelm);
      void _initialize(int32_t _577_numvar,int32_t _578_numcon,int32_t _579_numbarvar,int32_t _580_numacc,int32_t _581_numaccelm);
      static SolutionStruct::t _new_SolutionStruct(monty::rc_ptr< ::mosek::fusion::SolutionStruct > _582_that);
      void _initialize(monty::rc_ptr< ::mosek::fusion::SolutionStruct > _582_that);
      virtual monty::rc_ptr< ::mosek::fusion::SolutionStruct > __mosek_2fusion_2SolutionStruct__clone() ;
      virtual void resize(int32_t _583_numvar,int32_t _584_numcon,int32_t _585_numbarvar,int32_t _586_numacc,int32_t _587_numaccelm) ;
      virtual bool isDualAcceptable(mosek::fusion::AccSolutionStatus _608_acceptable_sol) ;
      virtual bool isPrimalAcceptable(mosek::fusion::AccSolutionStatus _609_acceptable_sol) ;
      virtual bool isAcceptable(mosek::fusion::SolutionStatus _610_stat,mosek::fusion::AccSolutionStatus _611_accstat) ;
    }; // struct SolutionStruct;

    struct p_RowBlockManager
    {
      RowBlockManager * _pubthis;
      static mosek::fusion::p_RowBlockManager* _get_impl(mosek::fusion::RowBlockManager * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_RowBlockManager * _get_impl(mosek::fusion::RowBlockManager::t _inst) { return _get_impl(_inst.get()); }
      p_RowBlockManager(RowBlockManager * _pubthis);
      virtual ~p_RowBlockManager() { /* std::cout << "~p_RowBlockManager" << std::endl;*/ };
      int32_t varidx_used{};
      int32_t code_used{};
      std::shared_ptr< monty::ndarray< double,1 > > cconst{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > code{};
      int32_t first_free_codeitem{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > param_code_sizes{};
      std::shared_ptr< monty::ndarray< int64_t,1 > > param_varidx{};
      int32_t first_free_entry{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > row_code_ptr{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > row_param_ptre{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > row_param_ptrb{};
      monty::rc_ptr< ::mosek::fusion::LinkedBlocks > blocks{};

      virtual void destroy();

      static RowBlockManager::t _new_RowBlockManager(monty::rc_ptr< ::mosek::fusion::RowBlockManager > _612_that);
      void _initialize(monty::rc_ptr< ::mosek::fusion::RowBlockManager > _612_that);
      static RowBlockManager::t _new_RowBlockManager();
      void _initialize();
      virtual int32_t num_parameterized() ;
      virtual bool is_parameterized() ;
      virtual int32_t blocksize(int32_t _613_id) ;
      virtual int32_t block_capacity() ;
      virtual int32_t capacity() ;
      virtual void get(int32_t _614_id,std::shared_ptr< monty::ndarray< int32_t,1 > > _615_target,int32_t _616_offset) ;
      virtual void evaluate(monty::rc_ptr< ::mosek::fusion::WorkStack > _617_xs,std::shared_ptr< monty::ndarray< double,1 > > _618_param_value,std::shared_ptr< monty::ndarray< int32_t,1 > > _619_subi,std::shared_ptr< monty::ndarray< int32_t,1 > > _620_subj,std::shared_ptr< monty::ndarray< double,1 > > _621_val) ;
      virtual void replace_row_code(monty::rc_ptr< ::mosek::fusion::WorkStack > _632_rs,std::shared_ptr< monty::ndarray< int32_t,1 > > _633_nativeidxs,int32_t _634_ptr,int32_t _635_nidxs,int32_t _636_codeptr,int32_t _637_code_p,int32_t _638_cconst_p) ;
      virtual void clear_row_code(std::shared_ptr< monty::ndarray< int32_t,1 > > _661_nativeidxs) ;
      virtual void compress() ;
      virtual void ensure_code_cap(int32_t _674_nentry,int32_t _675_codesize) ;
      virtual void release(int32_t _685_id,std::shared_ptr< monty::ndarray< int32_t,1 > > _686_nativeidxs) ;
      virtual int32_t allocate(std::shared_ptr< monty::ndarray< int32_t,1 > > _690_nativeidxs) ;
      virtual bool row_is_parameterized(int32_t _696_i) ;
    }; // struct RowBlockManager;

    struct p_BaseVariable : public /*implements*/ virtual ::mosek::fusion::Variable
    {
      BaseVariable * _pubthis;
      static mosek::fusion::p_BaseVariable* _get_impl(mosek::fusion::BaseVariable * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_BaseVariable * _get_impl(mosek::fusion::BaseVariable::t _inst) { return _get_impl(_inst.get()); }
      p_BaseVariable(BaseVariable * _pubthis);
      virtual ~p_BaseVariable() { /* std::cout << "~p_BaseVariable" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< int64_t,1 > > sparsity{};
      std::shared_ptr< monty::ndarray< int64_t,1 > > basevar_nativeidxs{};
      monty::rc_ptr< ::mosek::fusion::Model > model{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};

      virtual void destroy();

      static BaseVariable::t _new_BaseVariable(monty::rc_ptr< ::mosek::fusion::BaseVariable > _2268_v,monty::rc_ptr< ::mosek::fusion::Model > _2269_m);
      void _initialize(monty::rc_ptr< ::mosek::fusion::BaseVariable > _2268_v,monty::rc_ptr< ::mosek::fusion::Model > _2269_m);
      static BaseVariable::t _new_BaseVariable(monty::rc_ptr< ::mosek::fusion::Model > _2270_m,std::shared_ptr< monty::ndarray< int32_t,1 > > _2271_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _2272_sparsity,std::shared_ptr< monty::ndarray< int64_t,1 > > _2273_basevar_nativeidxs);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Model > _2270_m,std::shared_ptr< monty::ndarray< int32_t,1 > > _2271_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _2272_sparsity,std::shared_ptr< monty::ndarray< int64_t,1 > > _2273_basevar_nativeidxs);
      virtual /* override */ std::string toString() ;
      virtual void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _2276_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _2277_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _2278_xs) ;
      virtual void remove() ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__fromTril(int32_t _2296_dim0,int32_t _2297_d) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__fromTril(int32_t _2330_d) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Variable__fromTril(int32_t _2330_d) { return __mosek_2fusion_2BaseVariable__fromTril(_2330_d); }
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__tril(int32_t _2331_dim1,int32_t _2332_dim2) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__tril() ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Variable__tril() { return __mosek_2fusion_2BaseVariable__tril(); }
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__reshape(int32_t _2386_dim0,int32_t _2387_dim1,int32_t _2388_dim2) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Variable__reshape(int32_t _2386_dim0,int32_t _2387_dim1,int32_t _2388_dim2) { return __mosek_2fusion_2BaseVariable__reshape(_2386_dim0,_2387_dim1,_2388_dim2); }
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__reshape(int32_t _2389_dim0,int32_t _2390_dim1) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Variable__reshape(int32_t _2389_dim0,int32_t _2390_dim1) { return __mosek_2fusion_2BaseVariable__reshape(_2389_dim0,_2390_dim1); }
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__reshape(int32_t _2391_dim0) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Variable__reshape(int32_t _2391_dim0) { return __mosek_2fusion_2BaseVariable__reshape(_2391_dim0); }
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__reshape(std::shared_ptr< monty::ndarray< int32_t,1 > > _2392_shape) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Variable__reshape(std::shared_ptr< monty::ndarray< int32_t,1 > > _2392_shape) { return __mosek_2fusion_2BaseVariable__reshape(_2392_shape); }
      virtual void setLevel(std::shared_ptr< monty::ndarray< double,1 > > _2396_v) ;
      virtual monty::rc_ptr< ::mosek::fusion::Model > __mosek_2fusion_2BaseVariable__getModel() ;
      virtual monty::rc_ptr< ::mosek::fusion::Model > __mosek_2fusion_2Expression__getModel() { return __mosek_2fusion_2BaseVariable__getModel(); }
      virtual int32_t getND() ;
      virtual int32_t getDim(int32_t _2399_d) ;
      virtual std::shared_ptr< monty::ndarray< int32_t,1 > > getShape() ;
      virtual int64_t getSize() ;
      virtual std::shared_ptr< monty::ndarray< double,1 > > dual() ;
      virtual std::shared_ptr< monty::ndarray< double,1 > > level() ;
      virtual void makeContinuous() ;
      virtual void makeInteger() ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__transpose() ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Variable__transpose() { return __mosek_2fusion_2BaseVariable__transpose(); }
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__index(int32_t _2420_i0,int32_t _2421_i1,int32_t _2422_i2) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Variable__index(int32_t _2420_i0,int32_t _2421_i1,int32_t _2422_i2) { return __mosek_2fusion_2BaseVariable__index(_2420_i0,_2421_i1,_2422_i2); }
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__index(int32_t _2423_i0,int32_t _2424_i1) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Variable__index(int32_t _2423_i0,int32_t _2424_i1) { return __mosek_2fusion_2BaseVariable__index(_2423_i0,_2424_i1); }
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__index(std::shared_ptr< monty::ndarray< int32_t,1 > > _2425_index) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Variable__index(std::shared_ptr< monty::ndarray< int32_t,1 > > _2425_index) { return __mosek_2fusion_2BaseVariable__index(_2425_index); }
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__index(int32_t _2428_index) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Variable__index(int32_t _2428_index) { return __mosek_2fusion_2BaseVariable__index(_2428_index); }
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__pick(std::shared_ptr< monty::ndarray< int32_t,1 > > _2429_i0,std::shared_ptr< monty::ndarray< int32_t,1 > > _2430_i1,std::shared_ptr< monty::ndarray< int32_t,1 > > _2431_i2) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Variable__pick(std::shared_ptr< monty::ndarray< int32_t,1 > > _2429_i0,std::shared_ptr< monty::ndarray< int32_t,1 > > _2430_i1,std::shared_ptr< monty::ndarray< int32_t,1 > > _2431_i2) { return __mosek_2fusion_2BaseVariable__pick(_2429_i0,_2430_i1,_2431_i2); }
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__pick(std::shared_ptr< monty::ndarray< int32_t,1 > > _2434_i0,std::shared_ptr< monty::ndarray< int32_t,1 > > _2435_i1) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Variable__pick(std::shared_ptr< monty::ndarray< int32_t,1 > > _2434_i0,std::shared_ptr< monty::ndarray< int32_t,1 > > _2435_i1) { return __mosek_2fusion_2BaseVariable__pick(_2434_i0,_2435_i1); }
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__pick(std::shared_ptr< monty::ndarray< int32_t,2 > > _2438_midxs) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Variable__pick(std::shared_ptr< monty::ndarray< int32_t,2 > > _2438_midxs) { return __mosek_2fusion_2BaseVariable__pick(_2438_midxs); }
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__pick(std::shared_ptr< monty::ndarray< int32_t,1 > > _2460_idxs) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Variable__pick(std::shared_ptr< monty::ndarray< int32_t,1 > > _2460_idxs) { return __mosek_2fusion_2BaseVariable__pick(_2460_idxs); }
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__antidiag(int32_t _2471_index) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Variable__antidiag(int32_t _2471_index) { return __mosek_2fusion_2BaseVariable__antidiag(_2471_index); }
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__antidiag() ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Variable__antidiag() { return __mosek_2fusion_2BaseVariable__antidiag(); }
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__diag(int32_t _2472_index) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Variable__diag(int32_t _2472_index) { return __mosek_2fusion_2BaseVariable__diag(_2472_index); }
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__diag() ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Variable__diag() { return __mosek_2fusion_2BaseVariable__diag(); }
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__general_diag(std::shared_ptr< monty::ndarray< int32_t,1 > > _2473_firstidx,std::shared_ptr< monty::ndarray< int32_t,1 > > _2474_step,int32_t _2475_num) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__slice(std::shared_ptr< monty::ndarray< int32_t,1 > > _2496_first,std::shared_ptr< monty::ndarray< int32_t,1 > > _2497_last) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Variable__slice(std::shared_ptr< monty::ndarray< int32_t,1 > > _2496_first,std::shared_ptr< monty::ndarray< int32_t,1 > > _2497_last) { return __mosek_2fusion_2BaseVariable__slice(_2496_first,_2497_last); }
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__slice(int32_t _2531_first,int32_t _2532_last) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2Variable__slice(int32_t _2531_first,int32_t _2532_last) { return __mosek_2fusion_2BaseVariable__slice(_2531_first,_2532_last); }
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2BaseVariable__asExpr() ;
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2Variable__asExpr() { return __mosek_2fusion_2BaseVariable__asExpr(); }
      virtual int32_t inst(int32_t _2541_spoffset,std::shared_ptr< monty::ndarray< int64_t,1 > > _2542_sparsity,int32_t _2543_nioffset,std::shared_ptr< monty::ndarray< int64_t,1 > > _2544_basevar_nativeidxs) ;
      virtual int32_t numInst() ;
      virtual void inst(int32_t _2549_offset,std::shared_ptr< monty::ndarray< int64_t,1 > > _2550_nindex) ;
      virtual void set_values(std::shared_ptr< monty::ndarray< double,1 > > _2557_values,bool _2558_primal) ;
      virtual void dual_lu(int32_t _2563_offset,std::shared_ptr< monty::ndarray< double,1 > > _2564_target,bool _2565_lower) ;
      virtual void values(int32_t _2568_offset,std::shared_ptr< monty::ndarray< double,1 > > _2569_target,bool _2570_primal) ;
      virtual void make_continuous() ;
      virtual void make_integer() ;
    }; // struct BaseVariable;

    struct p_SliceVariable : public ::mosek::fusion::p_BaseVariable
    {
      SliceVariable * _pubthis;
      static mosek::fusion::p_SliceVariable* _get_impl(mosek::fusion::SliceVariable * _inst){ return static_cast< mosek::fusion::p_SliceVariable* >(mosek::fusion::p_BaseVariable::_get_impl(_inst)); }
      static mosek::fusion::p_SliceVariable * _get_impl(mosek::fusion::SliceVariable::t _inst) { return _get_impl(_inst.get()); }
      p_SliceVariable(SliceVariable * _pubthis);
      virtual ~p_SliceVariable() { /* std::cout << "~p_SliceVariable" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};
      std::shared_ptr< monty::ndarray< int64_t,1 > > sparsity{};
      std::shared_ptr< monty::ndarray< int64_t,1 > > nativeidxs{};

      virtual void destroy();

      static SliceVariable::t _new_SliceVariable(monty::rc_ptr< ::mosek::fusion::Model > _2121_m,std::shared_ptr< monty::ndarray< int32_t,1 > > _2122_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _2123_sparsity,std::shared_ptr< monty::ndarray< int64_t,1 > > _2124_nativeidxs);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Model > _2121_m,std::shared_ptr< monty::ndarray< int32_t,1 > > _2122_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _2123_sparsity,std::shared_ptr< monty::ndarray< int64_t,1 > > _2124_nativeidxs);
      static SliceVariable::t _new_SliceVariable(monty::rc_ptr< ::mosek::fusion::SliceVariable > _2125_v);
      void _initialize(monty::rc_ptr< ::mosek::fusion::SliceVariable > _2125_v);
    }; // struct SliceVariable;

    struct p_BoundInterfaceVariable : public ::mosek::fusion::p_SliceVariable
    {
      BoundInterfaceVariable * _pubthis;
      static mosek::fusion::p_BoundInterfaceVariable* _get_impl(mosek::fusion::BoundInterfaceVariable * _inst){ return static_cast< mosek::fusion::p_BoundInterfaceVariable* >(mosek::fusion::p_SliceVariable::_get_impl(_inst)); }
      static mosek::fusion::p_BoundInterfaceVariable * _get_impl(mosek::fusion::BoundInterfaceVariable::t _inst) { return _get_impl(_inst.get()); }
      p_BoundInterfaceVariable(BoundInterfaceVariable * _pubthis);
      virtual ~p_BoundInterfaceVariable() { /* std::cout << "~p_BoundInterfaceVariable" << std::endl;*/ };
      bool islower{};

      virtual void destroy();

      static BoundInterfaceVariable::t _new_BoundInterfaceVariable(monty::rc_ptr< ::mosek::fusion::Model > _2095_m,std::shared_ptr< monty::ndarray< int32_t,1 > > _2096_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _2097_sparsity,std::shared_ptr< monty::ndarray< int64_t,1 > > _2098_nativeidxs,bool _2099_islower);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Model > _2095_m,std::shared_ptr< monty::ndarray< int32_t,1 > > _2096_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _2097_sparsity,std::shared_ptr< monty::ndarray< int64_t,1 > > _2098_nativeidxs,bool _2099_islower);
      static BoundInterfaceVariable::t _new_BoundInterfaceVariable(monty::rc_ptr< ::mosek::fusion::SliceVariable > _2100_v,bool _2101_islower);
      void _initialize(monty::rc_ptr< ::mosek::fusion::SliceVariable > _2100_v,bool _2101_islower);
      virtual /* override */ std::shared_ptr< monty::ndarray< double,1 > > dual() ;
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BoundInterfaceVariable__transpose() ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__transpose() { return __mosek_2fusion_2BoundInterfaceVariable__transpose(); }
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BoundInterfaceVariable__pick(std::shared_ptr< monty::ndarray< int32_t,1 > > _2103_i0,std::shared_ptr< monty::ndarray< int32_t,1 > > _2104_i1,std::shared_ptr< monty::ndarray< int32_t,1 > > _2105_i2) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__pick(std::shared_ptr< monty::ndarray< int32_t,1 > > _2103_i0,std::shared_ptr< monty::ndarray< int32_t,1 > > _2104_i1,std::shared_ptr< monty::ndarray< int32_t,1 > > _2105_i2) { return __mosek_2fusion_2BoundInterfaceVariable__pick(_2103_i0,_2104_i1,_2105_i2); }
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BoundInterfaceVariable__pick(std::shared_ptr< monty::ndarray< int32_t,1 > > _2106_i0,std::shared_ptr< monty::ndarray< int32_t,1 > > _2107_i1) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__pick(std::shared_ptr< monty::ndarray< int32_t,1 > > _2106_i0,std::shared_ptr< monty::ndarray< int32_t,1 > > _2107_i1) { return __mosek_2fusion_2BoundInterfaceVariable__pick(_2106_i0,_2107_i1); }
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BoundInterfaceVariable__pick(std::shared_ptr< monty::ndarray< int32_t,2 > > _2108_midxs) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__pick(std::shared_ptr< monty::ndarray< int32_t,2 > > _2108_midxs) { return __mosek_2fusion_2BoundInterfaceVariable__pick(_2108_midxs); }
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BoundInterfaceVariable__pick(std::shared_ptr< monty::ndarray< int32_t,1 > > _2109_idxs) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__pick(std::shared_ptr< monty::ndarray< int32_t,1 > > _2109_idxs) { return __mosek_2fusion_2BoundInterfaceVariable__pick(_2109_idxs); }
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BoundInterfaceVariable__antidiag(int32_t _2110_index) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__antidiag(int32_t _2110_index) { return __mosek_2fusion_2BoundInterfaceVariable__antidiag(_2110_index); }
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BoundInterfaceVariable__antidiag() ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__antidiag() { return __mosek_2fusion_2BoundInterfaceVariable__antidiag(); }
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BoundInterfaceVariable__diag(int32_t _2111_index) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__diag(int32_t _2111_index) { return __mosek_2fusion_2BoundInterfaceVariable__diag(_2111_index); }
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BoundInterfaceVariable__diag() ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__diag() { return __mosek_2fusion_2BoundInterfaceVariable__diag(); }
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BoundInterfaceVariable__slice(std::shared_ptr< monty::ndarray< int32_t,1 > > _2112_firsta,std::shared_ptr< monty::ndarray< int32_t,1 > > _2113_lasta) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__slice(std::shared_ptr< monty::ndarray< int32_t,1 > > _2112_firsta,std::shared_ptr< monty::ndarray< int32_t,1 > > _2113_lasta) { return __mosek_2fusion_2BoundInterfaceVariable__slice(_2112_firsta,_2113_lasta); }
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BoundInterfaceVariable__slice(int32_t _2114_first,int32_t _2115_last) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__slice(int32_t _2114_first,int32_t _2115_last) { return __mosek_2fusion_2BoundInterfaceVariable__slice(_2114_first,_2115_last); }
      virtual monty::rc_ptr< ::mosek::fusion::BoundInterfaceVariable > __mosek_2fusion_2BoundInterfaceVariable__from_(monty::rc_ptr< ::mosek::fusion::Variable > _2116_v) ;
    }; // struct BoundInterfaceVariable;

    struct p_ModelVariable : public ::mosek::fusion::p_BaseVariable
    {
      ModelVariable * _pubthis;
      static mosek::fusion::p_ModelVariable* _get_impl(mosek::fusion::ModelVariable * _inst){ return static_cast< mosek::fusion::p_ModelVariable* >(mosek::fusion::p_BaseVariable::_get_impl(_inst)); }
      static mosek::fusion::p_ModelVariable * _get_impl(mosek::fusion::ModelVariable::t _inst) { return _get_impl(_inst.get()); }
      p_ModelVariable(ModelVariable * _pubthis);
      virtual ~p_ModelVariable() { /* std::cout << "~p_ModelVariable" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< int64_t,1 > > sparsity{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};
      std::shared_ptr< monty::ndarray< int64_t,1 > > modelvar_nativeidxs{};
      int64_t varid{};
      std::string name{};

      virtual void destroy();

      static ModelVariable::t _new_ModelVariable(monty::rc_ptr< ::mosek::fusion::ModelVariable > _2231_v,monty::rc_ptr< ::mosek::fusion::Model > _2232_m);
      void _initialize(monty::rc_ptr< ::mosek::fusion::ModelVariable > _2231_v,monty::rc_ptr< ::mosek::fusion::Model > _2232_m);
      static ModelVariable::t _new_ModelVariable(monty::rc_ptr< ::mosek::fusion::Model > _2233_model,const std::string &  _2234_name,std::shared_ptr< monty::ndarray< int32_t,1 > > _2235_shape,int64_t _2236_varid,std::shared_ptr< monty::ndarray< int64_t,1 > > _2237_sparsity,std::shared_ptr< monty::ndarray< int64_t,1 > > _2238_modelvar_nativeidxs);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Model > _2233_model,const std::string &  _2234_name,std::shared_ptr< monty::ndarray< int32_t,1 > > _2235_shape,int64_t _2236_varid,std::shared_ptr< monty::ndarray< int64_t,1 > > _2237_sparsity,std::shared_ptr< monty::ndarray< int64_t,1 > > _2238_modelvar_nativeidxs);
      virtual void elementName(int64_t _2239_index,monty::rc_ptr< ::mosek::fusion::Utils::StringBuffer > _2240_sb) ;
      virtual /* override */ void remove() ;
      virtual monty::rc_ptr< ::mosek::fusion::ModelVariable > __mosek_2fusion_2ModelVariable__clone(monty::rc_ptr< ::mosek::fusion::Model > _2241_m) { throw monty::AbstractClassError("Call to abstract method"); }
    }; // struct ModelVariable;

    struct p_RangedVariable : public ::mosek::fusion::p_ModelVariable
    {
      RangedVariable * _pubthis;
      static mosek::fusion::p_RangedVariable* _get_impl(mosek::fusion::RangedVariable * _inst){ return static_cast< mosek::fusion::p_RangedVariable* >(mosek::fusion::p_ModelVariable::_get_impl(_inst)); }
      static mosek::fusion::p_RangedVariable * _get_impl(mosek::fusion::RangedVariable::t _inst) { return _get_impl(_inst.get()); }
      p_RangedVariable(RangedVariable * _pubthis);
      virtual ~p_RangedVariable() { /* std::cout << "~p_RangedVariable" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};
      std::string name{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > nativeidxs{};
      std::shared_ptr< monty::ndarray< int64_t,1 > > sparsity{};

      virtual void destroy();

      static RangedVariable::t _new_RangedVariable(monty::rc_ptr< ::mosek::fusion::RangedVariable > _2126_v,monty::rc_ptr< ::mosek::fusion::Model > _2127_m);
      void _initialize(monty::rc_ptr< ::mosek::fusion::RangedVariable > _2126_v,monty::rc_ptr< ::mosek::fusion::Model > _2127_m);
      static RangedVariable::t _new_RangedVariable(monty::rc_ptr< ::mosek::fusion::Model > _2128_model,const std::string &  _2129_name,int64_t _2130_varid,std::shared_ptr< monty::ndarray< int32_t,1 > > _2131_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _2132_sparsity,std::shared_ptr< monty::ndarray< int32_t,1 > > _2133_nativeidxs);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Model > _2128_model,const std::string &  _2129_name,int64_t _2130_varid,std::shared_ptr< monty::ndarray< int32_t,1 > > _2131_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _2132_sparsity,std::shared_ptr< monty::ndarray< int32_t,1 > > _2133_nativeidxs);
      virtual monty::rc_ptr< ::mosek::fusion::Utils::StringBuffer > __mosek_2fusion_2RangedVariable__elementDesc(int64_t _2134_index,monty::rc_ptr< ::mosek::fusion::Utils::StringBuffer > _2135_sb) ;
      virtual void dual_u(int32_t _2136_offset,std::shared_ptr< monty::ndarray< double,1 > > _2137_target) ;
      virtual void dual_l(int32_t _2138_offset,std::shared_ptr< monty::ndarray< double,1 > > _2139_target) ;
      virtual monty::rc_ptr< ::mosek::fusion::BoundInterfaceVariable > __mosek_2fusion_2RangedVariable__upperBoundVar() ;
      virtual monty::rc_ptr< ::mosek::fusion::BoundInterfaceVariable > __mosek_2fusion_2RangedVariable__lowerBoundVar() ;
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::ModelVariable > __mosek_2fusion_2RangedVariable__clone(monty::rc_ptr< ::mosek::fusion::Model > _2142_m) ;
      virtual monty::rc_ptr< ::mosek::fusion::ModelVariable > __mosek_2fusion_2ModelVariable__clone(monty::rc_ptr< ::mosek::fusion::Model > _2142_m) { return __mosek_2fusion_2RangedVariable__clone(_2142_m); }
      static  std::shared_ptr< monty::ndarray< int64_t,1 > > globalNativeIndexes(std::shared_ptr< monty::ndarray< int32_t,1 > > _2143_nativeidxs);
    }; // struct RangedVariable;

    struct p_LinearPSDVariable : public ::mosek::fusion::p_ModelVariable
    {
      LinearPSDVariable * _pubthis;
      static mosek::fusion::p_LinearPSDVariable* _get_impl(mosek::fusion::LinearPSDVariable * _inst){ return static_cast< mosek::fusion::p_LinearPSDVariable* >(mosek::fusion::p_ModelVariable::_get_impl(_inst)); }
      static mosek::fusion::p_LinearPSDVariable * _get_impl(mosek::fusion::LinearPSDVariable::t _inst) { return _get_impl(_inst.get()); }
      p_LinearPSDVariable(LinearPSDVariable * _pubthis);
      virtual ~p_LinearPSDVariable() { /* std::cout << "~p_LinearPSDVariable" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};
      std::string name{};
      int32_t varid{};
      std::shared_ptr< monty::ndarray< int64_t,1 > > nativeidxs{};
      int32_t conedim{};

      virtual void destroy();

      static LinearPSDVariable::t _new_LinearPSDVariable(monty::rc_ptr< ::mosek::fusion::LinearPSDVariable > _2146_v,monty::rc_ptr< ::mosek::fusion::Model > _2147_m);
      void _initialize(monty::rc_ptr< ::mosek::fusion::LinearPSDVariable > _2146_v,monty::rc_ptr< ::mosek::fusion::Model > _2147_m);
      static LinearPSDVariable::t _new_LinearPSDVariable(monty::rc_ptr< ::mosek::fusion::Model > _2148_model,const std::string &  _2149_name,int32_t _2150_varid,std::shared_ptr< monty::ndarray< int32_t,1 > > _2151_shape,int32_t _2152_conedim,std::shared_ptr< monty::ndarray< int64_t,1 > > _2153_nativeidxs);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Model > _2148_model,const std::string &  _2149_name,int32_t _2150_varid,std::shared_ptr< monty::ndarray< int32_t,1 > > _2151_shape,int32_t _2152_conedim,std::shared_ptr< monty::ndarray< int64_t,1 > > _2153_nativeidxs);
      virtual /* override */ std::string toString() ;
      virtual void make_continuous(std::shared_ptr< monty::ndarray< int64_t,1 > > _2156_idxs) ;
      virtual void make_integer(std::shared_ptr< monty::ndarray< int64_t,1 > > _2157_idxs) ;
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::ModelVariable > __mosek_2fusion_2LinearPSDVariable__clone(monty::rc_ptr< ::mosek::fusion::Model > _2158_m) ;
      virtual monty::rc_ptr< ::mosek::fusion::ModelVariable > __mosek_2fusion_2ModelVariable__clone(monty::rc_ptr< ::mosek::fusion::Model > _2158_m) { return __mosek_2fusion_2LinearPSDVariable__clone(_2158_m); }
      static  std::shared_ptr< monty::ndarray< int64_t,1 > > globalNativeIndexes(std::shared_ptr< monty::ndarray< int64_t,1 > > _2159_nativeidxs);
    }; // struct LinearPSDVariable;

    struct p_PSDVariable : public ::mosek::fusion::p_ModelVariable
    {
      PSDVariable * _pubthis;
      static mosek::fusion::p_PSDVariable* _get_impl(mosek::fusion::PSDVariable * _inst){ return static_cast< mosek::fusion::p_PSDVariable* >(mosek::fusion::p_ModelVariable::_get_impl(_inst)); }
      static mosek::fusion::p_PSDVariable * _get_impl(mosek::fusion::PSDVariable::t _inst) { return _get_impl(_inst.get()); }
      p_PSDVariable(PSDVariable * _pubthis);
      virtual ~p_PSDVariable() { /* std::cout << "~p_PSDVariable" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Model > model{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > barvaridxs{};
      int32_t conedim2{};
      int32_t conedim1{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};
      std::string name{};
      std::shared_ptr< monty::ndarray< int64_t,1 > > nativeidxs{};
      int32_t varid{};

      virtual void destroy();

      static PSDVariable::t _new_PSDVariable(monty::rc_ptr< ::mosek::fusion::PSDVariable > _2161_v,monty::rc_ptr< ::mosek::fusion::Model > _2162_m);
      void _initialize(monty::rc_ptr< ::mosek::fusion::PSDVariable > _2161_v,monty::rc_ptr< ::mosek::fusion::Model > _2162_m);
      static PSDVariable::t _new_PSDVariable(monty::rc_ptr< ::mosek::fusion::Model > _2163_model,const std::string &  _2164_name,int32_t _2165_varid,std::shared_ptr< monty::ndarray< int32_t,1 > > _2166_shape,int32_t _2167_conedim1,int32_t _2168_conedim2,std::shared_ptr< monty::ndarray< int32_t,1 > > _2169_barvaridxs,std::shared_ptr< monty::ndarray< int64_t,1 > > _2170_nativeidxs);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Model > _2163_model,const std::string &  _2164_name,int32_t _2165_varid,std::shared_ptr< monty::ndarray< int32_t,1 > > _2166_shape,int32_t _2167_conedim1,int32_t _2168_conedim2,std::shared_ptr< monty::ndarray< int32_t,1 > > _2169_barvaridxs,std::shared_ptr< monty::ndarray< int64_t,1 > > _2170_nativeidxs);
      virtual /* override */ std::string toString() ;
      virtual monty::rc_ptr< ::mosek::fusion::Utils::StringBuffer > __mosek_2fusion_2PSDVariable__elementDesc(int64_t _2173_index,monty::rc_ptr< ::mosek::fusion::Utils::StringBuffer > _2174_sb) ;
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::ModelVariable > __mosek_2fusion_2PSDVariable__clone(monty::rc_ptr< ::mosek::fusion::Model > _2175_m) ;
      virtual monty::rc_ptr< ::mosek::fusion::ModelVariable > __mosek_2fusion_2ModelVariable__clone(monty::rc_ptr< ::mosek::fusion::Model > _2175_m) { return __mosek_2fusion_2PSDVariable__clone(_2175_m); }
      static  std::shared_ptr< monty::ndarray< int64_t,1 > > fullnativeidxs(std::shared_ptr< monty::ndarray< int32_t,1 > > _2176_shape,int32_t _2177_conedim1,int32_t _2178_conedim2,std::shared_ptr< monty::ndarray< int64_t,1 > > _2179_nativeidxs);
    }; // struct PSDVariable;

    struct p_LinearVariable : public ::mosek::fusion::p_ModelVariable
    {
      LinearVariable * _pubthis;
      static mosek::fusion::p_LinearVariable* _get_impl(mosek::fusion::LinearVariable * _inst){ return static_cast< mosek::fusion::p_LinearVariable* >(mosek::fusion::p_ModelVariable::_get_impl(_inst)); }
      static mosek::fusion::p_LinearVariable * _get_impl(mosek::fusion::LinearVariable::t _inst) { return _get_impl(_inst.get()); }
      p_LinearVariable(LinearVariable * _pubthis);
      virtual ~p_LinearVariable() { /* std::cout << "~p_LinearVariable" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};
      std::shared_ptr< monty::ndarray< int64_t,1 > > sparsity{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > nativeidxs{};
      std::string name{};

      virtual void destroy();

      static LinearVariable::t _new_LinearVariable(monty::rc_ptr< ::mosek::fusion::LinearVariable > _2204_v,monty::rc_ptr< ::mosek::fusion::Model > _2205_m);
      void _initialize(monty::rc_ptr< ::mosek::fusion::LinearVariable > _2204_v,monty::rc_ptr< ::mosek::fusion::Model > _2205_m);
      static LinearVariable::t _new_LinearVariable(monty::rc_ptr< ::mosek::fusion::Model > _2206_model,const std::string &  _2207_name,int64_t _2208_varid,std::shared_ptr< monty::ndarray< int32_t,1 > > _2209_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _2210_sparsity,std::shared_ptr< monty::ndarray< int32_t,1 > > _2211_nativeidxs);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Model > _2206_model,const std::string &  _2207_name,int64_t _2208_varid,std::shared_ptr< monty::ndarray< int32_t,1 > > _2209_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _2210_sparsity,std::shared_ptr< monty::ndarray< int32_t,1 > > _2211_nativeidxs);
      virtual /* override */ std::string toString() ;
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::ModelVariable > __mosek_2fusion_2LinearVariable__clone(monty::rc_ptr< ::mosek::fusion::Model > _2214_m) ;
      virtual monty::rc_ptr< ::mosek::fusion::ModelVariable > __mosek_2fusion_2ModelVariable__clone(monty::rc_ptr< ::mosek::fusion::Model > _2214_m) { return __mosek_2fusion_2LinearVariable__clone(_2214_m); }
      static  std::shared_ptr< monty::ndarray< int64_t,1 > > globalNativeIndexes(std::shared_ptr< monty::ndarray< int32_t,1 > > _2215_nativeidxs);
    }; // struct LinearVariable;

    struct p_ConicVariable : public ::mosek::fusion::p_ModelVariable
    {
      ConicVariable * _pubthis;
      static mosek::fusion::p_ConicVariable* _get_impl(mosek::fusion::ConicVariable * _inst){ return static_cast< mosek::fusion::p_ConicVariable* >(mosek::fusion::p_ModelVariable::_get_impl(_inst)); }
      static mosek::fusion::p_ConicVariable * _get_impl(mosek::fusion::ConicVariable::t _inst) { return _get_impl(_inst.get()); }
      p_ConicVariable(ConicVariable * _pubthis);
      virtual ~p_ConicVariable() { /* std::cout << "~p_ConicVariable" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< int32_t,1 > > nativeidxs{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};
      std::string name{};
      int32_t varid{};

      virtual void destroy();

      static ConicVariable::t _new_ConicVariable(monty::rc_ptr< ::mosek::fusion::ConicVariable > _2218_v,monty::rc_ptr< ::mosek::fusion::Model > _2219_m);
      void _initialize(monty::rc_ptr< ::mosek::fusion::ConicVariable > _2218_v,monty::rc_ptr< ::mosek::fusion::Model > _2219_m);
      static ConicVariable::t _new_ConicVariable(monty::rc_ptr< ::mosek::fusion::Model > _2220_model,const std::string &  _2221_name,int32_t _2222_varid,std::shared_ptr< monty::ndarray< int32_t,1 > > _2223_shape,std::shared_ptr< monty::ndarray< int32_t,1 > > _2224_nativeidxs);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Model > _2220_model,const std::string &  _2221_name,int32_t _2222_varid,std::shared_ptr< monty::ndarray< int32_t,1 > > _2223_shape,std::shared_ptr< monty::ndarray< int32_t,1 > > _2224_nativeidxs);
      virtual /* override */ std::string toString() ;
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::ModelVariable > __mosek_2fusion_2ConicVariable__clone(monty::rc_ptr< ::mosek::fusion::Model > _2227_m) ;
      virtual monty::rc_ptr< ::mosek::fusion::ModelVariable > __mosek_2fusion_2ModelVariable__clone(monty::rc_ptr< ::mosek::fusion::Model > _2227_m) { return __mosek_2fusion_2ConicVariable__clone(_2227_m); }
      static  std::shared_ptr< monty::ndarray< int64_t,1 > > globalNativeIndexes(std::shared_ptr< monty::ndarray< int32_t,1 > > _2228_nativeidxs);
    }; // struct ConicVariable;

    struct p_NilVariable : public ::mosek::fusion::p_BaseVariable
    {
      NilVariable * _pubthis;
      static mosek::fusion::p_NilVariable* _get_impl(mosek::fusion::NilVariable * _inst){ return static_cast< mosek::fusion::p_NilVariable* >(mosek::fusion::p_BaseVariable::_get_impl(_inst)); }
      static mosek::fusion::p_NilVariable * _get_impl(mosek::fusion::NilVariable::t _inst) { return _get_impl(_inst.get()); }
      p_NilVariable(NilVariable * _pubthis);
      virtual ~p_NilVariable() { /* std::cout << "~p_NilVariable" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};

      virtual void destroy();

      static NilVariable::t _new_NilVariable(std::shared_ptr< monty::ndarray< int32_t,1 > > _2242_shape);
      void _initialize(std::shared_ptr< monty::ndarray< int32_t,1 > > _2242_shape);
      static NilVariable::t _new_NilVariable();
      void _initialize();
      virtual monty::rc_ptr< ::mosek::fusion::Utils::StringBuffer > __mosek_2fusion_2NilVariable__elementDesc(int64_t _2244_index,monty::rc_ptr< ::mosek::fusion::Utils::StringBuffer > _2245_sb) ;
      virtual void elementName(int64_t _2246_index,monty::rc_ptr< ::mosek::fusion::Utils::StringBuffer > _2247_sb) ;
      virtual /* override */ int32_t numInst() ;
      virtual int32_t inst(int32_t _2248_offset,std::shared_ptr< monty::ndarray< int64_t,1 > > _2249_sparsity,std::shared_ptr< monty::ndarray< int64_t,1 > > _2250_basevar_nativeidxs) ;
      virtual /* override */ void inst(int32_t _2251_offset,std::shared_ptr< monty::ndarray< int64_t,1 > > _2252_nindex) ;
      virtual /* override */ void set_values(std::shared_ptr< monty::ndarray< double,1 > > _2253_target,bool _2254_primal) ;
      virtual /* override */ void values(int32_t _2255_offset,std::shared_ptr< monty::ndarray< double,1 > > _2256_target,bool _2257_primal) ;
      virtual /* override */ void make_continuous() ;
      virtual /* override */ void make_integer() ;
      virtual /* override */ std::string toString() ;
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2NilVariable__index(std::shared_ptr< monty::ndarray< int32_t,1 > > _2258_first) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__index(std::shared_ptr< monty::ndarray< int32_t,1 > > _2258_first) { return __mosek_2fusion_2NilVariable__index(_2258_first); }
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2NilVariable__index(int32_t _2260_first) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__index(int32_t _2260_first) { return __mosek_2fusion_2NilVariable__index(_2260_first); }
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2NilVariable__slice(std::shared_ptr< monty::ndarray< int32_t,1 > > _2262_first,std::shared_ptr< monty::ndarray< int32_t,1 > > _2263_last) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__slice(std::shared_ptr< monty::ndarray< int32_t,1 > > _2262_first,std::shared_ptr< monty::ndarray< int32_t,1 > > _2263_last) { return __mosek_2fusion_2NilVariable__slice(_2262_first,_2263_last); }
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2NilVariable__slice(int32_t _2266_first,int32_t _2267_last) ;
      virtual monty::rc_ptr< ::mosek::fusion::Variable > __mosek_2fusion_2BaseVariable__slice(int32_t _2266_first,int32_t _2267_last) { return __mosek_2fusion_2NilVariable__slice(_2266_first,_2267_last); }
    }; // struct NilVariable;

    struct p_Var
    {
      Var * _pubthis;
      static mosek::fusion::p_Var* _get_impl(mosek::fusion::Var * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_Var * _get_impl(mosek::fusion::Var::t _inst) { return _get_impl(_inst.get()); }
      p_Var(Var * _pubthis);
      virtual ~p_Var() { /* std::cout << "~p_Var" << std::endl;*/ };

      virtual void destroy();

      static  monty::rc_ptr< ::mosek::fusion::Variable > empty(std::shared_ptr< monty::ndarray< int32_t,1 > > _2613_shape);
      static  monty::rc_ptr< ::mosek::fusion::Variable > compress(monty::rc_ptr< ::mosek::fusion::Variable > _2615_v);
      static  monty::rc_ptr< ::mosek::fusion::Variable > reshape(monty::rc_ptr< ::mosek::fusion::Variable > _2623_v,int32_t _2624_d1);
      static  monty::rc_ptr< ::mosek::fusion::Variable > reshape(monty::rc_ptr< ::mosek::fusion::Variable > _2625_v,int32_t _2626_d1,int32_t _2627_d2);
      static  monty::rc_ptr< ::mosek::fusion::Variable > flatten(monty::rc_ptr< ::mosek::fusion::Variable > _2628_v);
      static  monty::rc_ptr< ::mosek::fusion::Variable > reshape(monty::rc_ptr< ::mosek::fusion::Variable > _2629_v,std::shared_ptr< monty::ndarray< int32_t,1 > > _2630_shape);
      static  monty::rc_ptr< ::mosek::fusion::Variable > index_permute_(monty::rc_ptr< ::mosek::fusion::Variable > _2631_v,std::shared_ptr< monty::ndarray< int32_t,1 > > _2632_perm);
      static  monty::rc_ptr< ::mosek::fusion::Variable > hrepeat(monty::rc_ptr< ::mosek::fusion::Variable > _2661_v,int32_t _2662_n);
      static  monty::rc_ptr< ::mosek::fusion::Variable > vrepeat(monty::rc_ptr< ::mosek::fusion::Variable > _2663_v,int32_t _2664_n);
      static  monty::rc_ptr< ::mosek::fusion::Variable > repeat(monty::rc_ptr< ::mosek::fusion::Variable > _2665_v,int32_t _2666_n);
      static  monty::rc_ptr< ::mosek::fusion::Variable > repeat(monty::rc_ptr< ::mosek::fusion::Variable > _2667_v,int32_t _2668_n,int32_t _2669_dim);
      static  monty::rc_ptr< ::mosek::fusion::Variable > drepeat(monty::rc_ptr< ::mosek::fusion::Variable > _2670_v,int32_t _2671_dim,int32_t _2672_n);
      static  monty::rc_ptr< ::mosek::fusion::Variable > stack(std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Variable >,1 > >,1 > > _2742_vlist);
      static  monty::rc_ptr< ::mosek::fusion::Variable > vstack(monty::rc_ptr< ::mosek::fusion::Variable > _2744_v1,monty::rc_ptr< ::mosek::fusion::Variable > _2745_v2,monty::rc_ptr< ::mosek::fusion::Variable > _2746_v3);
      static  monty::rc_ptr< ::mosek::fusion::Variable > vstack(monty::rc_ptr< ::mosek::fusion::Variable > _2747_v1,monty::rc_ptr< ::mosek::fusion::Variable > _2748_v2);
      static  monty::rc_ptr< ::mosek::fusion::Variable > vstack(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Variable >,1 > > _2749_v);
      static  monty::rc_ptr< ::mosek::fusion::Variable > hstack(monty::rc_ptr< ::mosek::fusion::Variable > _2750_v1,monty::rc_ptr< ::mosek::fusion::Variable > _2751_v2,monty::rc_ptr< ::mosek::fusion::Variable > _2752_v3);
      static  monty::rc_ptr< ::mosek::fusion::Variable > hstack(monty::rc_ptr< ::mosek::fusion::Variable > _2753_v1,monty::rc_ptr< ::mosek::fusion::Variable > _2754_v2);
      static  monty::rc_ptr< ::mosek::fusion::Variable > hstack(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Variable >,1 > > _2755_v);
      static  monty::rc_ptr< ::mosek::fusion::Variable > stack(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Variable >,1 > > _2756_v,int32_t _2757_dim);
      static  monty::rc_ptr< ::mosek::fusion::Variable > stack(monty::rc_ptr< ::mosek::fusion::Variable > _2758_v1,monty::rc_ptr< ::mosek::fusion::Variable > _2759_v2,monty::rc_ptr< ::mosek::fusion::Variable > _2760_v3,int32_t _2761_dim);
      static  monty::rc_ptr< ::mosek::fusion::Variable > stack(monty::rc_ptr< ::mosek::fusion::Variable > _2762_v1,monty::rc_ptr< ::mosek::fusion::Variable > _2763_v2,int32_t _2764_dim);
      static  monty::rc_ptr< ::mosek::fusion::Variable > stack(int32_t _2765_dim,std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Variable >,1 > > _2766_v);
      static  monty::rc_ptr< ::mosek::fusion::Variable > stack(int32_t _2769_dim,monty::rc_ptr< ::mosek::fusion::Variable > _2770_v1,monty::rc_ptr< ::mosek::fusion::Variable > _2771_v2,monty::rc_ptr< ::mosek::fusion::Variable > _2772_v3);
      static  monty::rc_ptr< ::mosek::fusion::Variable > stack(int32_t _2773_dim,monty::rc_ptr< ::mosek::fusion::Variable > _2774_v1,monty::rc_ptr< ::mosek::fusion::Variable > _2775_v2);
      static  monty::rc_ptr< ::mosek::fusion::Variable > promote(monty::rc_ptr< ::mosek::fusion::Variable > _2776_v,int32_t _2777_nd);
      static  monty::rc_ptr< ::mosek::fusion::Variable > dstack(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Variable >,1 > > _2782_v,int32_t _2783_dim);
    }; // struct Var;

    struct p_Constraint
    {
      Constraint * _pubthis;
      static mosek::fusion::p_Constraint* _get_impl(mosek::fusion::Constraint * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_Constraint * _get_impl(mosek::fusion::Constraint::t _inst) { return _get_impl(_inst.get()); }
      p_Constraint(Constraint * _pubthis);
      virtual ~p_Constraint() { /* std::cout << "~p_Constraint" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< int32_t,1 > > con_nativeidxs{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};
      monty::rc_ptr< ::mosek::fusion::Model > model{};

      virtual void destroy();

      static Constraint::t _new_Constraint(monty::rc_ptr< ::mosek::fusion::Constraint > _2974_c,monty::rc_ptr< ::mosek::fusion::Model > _2975_m);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Constraint > _2974_c,monty::rc_ptr< ::mosek::fusion::Model > _2975_m);
      static Constraint::t _new_Constraint(monty::rc_ptr< ::mosek::fusion::Model > _2976_model,std::shared_ptr< monty::ndarray< int32_t,1 > > _2977_shape,std::shared_ptr< monty::ndarray< int32_t,1 > > _2978_con_nativeidxs);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Model > _2976_model,std::shared_ptr< monty::ndarray< int32_t,1 > > _2977_shape,std::shared_ptr< monty::ndarray< int32_t,1 > > _2978_con_nativeidxs);
      virtual /* override */ std::string toString() ;
      virtual void toStringArray(std::shared_ptr< monty::ndarray< int64_t,1 > > _2979_subi,int64_t _2980_dstidx,std::shared_ptr< monty::ndarray< std::string,1 > > _2981_result) ;
      virtual void dual_lu(int32_t _2982_offset,std::shared_ptr< monty::ndarray< double,1 > > _2983_target,bool _2984_islower) ;
      virtual std::shared_ptr< monty::ndarray< double,1 > > dual() ;
      virtual std::shared_ptr< monty::ndarray< double,1 > > level() ;
      virtual void values(bool _2987_primal,int32_t _2988_offset,std::shared_ptr< monty::ndarray< double,1 > > _2989_target) ;
      virtual void remove() ;
      virtual void update(std::shared_ptr< monty::ndarray< double,1 > > _2990_bfix) ;
      virtual void update(monty::rc_ptr< ::mosek::fusion::Expression > _2991_expr) ;
      virtual void update(monty::rc_ptr< ::mosek::fusion::Expression > _2995_expr,monty::rc_ptr< ::mosek::fusion::Variable > _2996_x,bool _2997_bfixupdate) ;
      virtual void update(monty::rc_ptr< ::mosek::fusion::Expression > _3017_expr,monty::rc_ptr< ::mosek::fusion::Variable > _3018_x) ;
      virtual monty::rc_ptr< ::mosek::fusion::Model > __mosek_2fusion_2Constraint__get_model() ;
      virtual int32_t get_nd() ;
      virtual int64_t size() ;
      static  monty::rc_ptr< ::mosek::fusion::Constraint > stack(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Constraint >,1 > > _3021_clist,int32_t _3022_dim);
      static  monty::rc_ptr< ::mosek::fusion::Constraint > stack(monty::rc_ptr< ::mosek::fusion::Constraint > _3023_v1,monty::rc_ptr< ::mosek::fusion::Constraint > _3024_v2,monty::rc_ptr< ::mosek::fusion::Constraint > _3025_v3,int32_t _3026_dim);
      static  monty::rc_ptr< ::mosek::fusion::Constraint > stack(monty::rc_ptr< ::mosek::fusion::Constraint > _3027_v1,monty::rc_ptr< ::mosek::fusion::Constraint > _3028_v2,int32_t _3029_dim);
      static  monty::rc_ptr< ::mosek::fusion::Constraint > hstack(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Constraint >,1 > > _3030_clist);
      static  monty::rc_ptr< ::mosek::fusion::Constraint > vstack(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Constraint >,1 > > _3031_clist);
      static  monty::rc_ptr< ::mosek::fusion::Constraint > hstack(monty::rc_ptr< ::mosek::fusion::Constraint > _3032_v1,monty::rc_ptr< ::mosek::fusion::Constraint > _3033_v2,monty::rc_ptr< ::mosek::fusion::Constraint > _3034_v3);
      static  monty::rc_ptr< ::mosek::fusion::Constraint > vstack(monty::rc_ptr< ::mosek::fusion::Constraint > _3035_v1,monty::rc_ptr< ::mosek::fusion::Constraint > _3036_v2,monty::rc_ptr< ::mosek::fusion::Constraint > _3037_v3);
      static  monty::rc_ptr< ::mosek::fusion::Constraint > hstack(monty::rc_ptr< ::mosek::fusion::Constraint > _3038_v1,monty::rc_ptr< ::mosek::fusion::Constraint > _3039_v2);
      static  monty::rc_ptr< ::mosek::fusion::Constraint > vstack(monty::rc_ptr< ::mosek::fusion::Constraint > _3040_v1,monty::rc_ptr< ::mosek::fusion::Constraint > _3041_v2);
      static  monty::rc_ptr< ::mosek::fusion::Constraint > dstack(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Constraint >,1 > > _3042_c,int32_t _3043_dim);
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Constraint__index(std::shared_ptr< monty::ndarray< int32_t,1 > > _3094_idxa) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Constraint__index(int32_t _3101_idx) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Constraint__slice(std::shared_ptr< monty::ndarray< int32_t,1 > > _3102_firsta,std::shared_ptr< monty::ndarray< int32_t,1 > > _3103_lasta) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Constraint__slice(int32_t _3122_first,int32_t _3123_last) ;
      virtual int32_t getND() ;
      virtual int32_t getSize() ;
      virtual monty::rc_ptr< ::mosek::fusion::Model > __mosek_2fusion_2Constraint__getModel() ;
      virtual std::shared_ptr< monty::ndarray< int32_t,1 > > getShape() ;
      virtual std::shared_ptr< monty::ndarray< int32_t,1 > > getNativeidxs() ;
    }; // struct Constraint;

    struct p_SliceConstraint : public ::mosek::fusion::p_Constraint
    {
      SliceConstraint * _pubthis;
      static mosek::fusion::p_SliceConstraint* _get_impl(mosek::fusion::SliceConstraint * _inst){ return static_cast< mosek::fusion::p_SliceConstraint* >(mosek::fusion::p_Constraint::_get_impl(_inst)); }
      static mosek::fusion::p_SliceConstraint * _get_impl(mosek::fusion::SliceConstraint::t _inst) { return _get_impl(_inst.get()); }
      p_SliceConstraint(SliceConstraint * _pubthis);
      virtual ~p_SliceConstraint() { /* std::cout << "~p_SliceConstraint" << std::endl;*/ };

      virtual void destroy();

      static SliceConstraint::t _new_SliceConstraint(monty::rc_ptr< ::mosek::fusion::SliceConstraint > _2926_c);
      void _initialize(monty::rc_ptr< ::mosek::fusion::SliceConstraint > _2926_c);
      static SliceConstraint::t _new_SliceConstraint(monty::rc_ptr< ::mosek::fusion::Model > _2927_model,std::shared_ptr< monty::ndarray< int32_t,1 > > _2928_shape,std::shared_ptr< monty::ndarray< int32_t,1 > > _2929_nativeidxs);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Model > _2927_model,std::shared_ptr< monty::ndarray< int32_t,1 > > _2928_shape,std::shared_ptr< monty::ndarray< int32_t,1 > > _2929_nativeidxs);
      virtual /* override */ std::string toString() ;
    }; // struct SliceConstraint;

    struct p_BoundInterfaceConstraint : public ::mosek::fusion::p_SliceConstraint
    {
      BoundInterfaceConstraint * _pubthis;
      static mosek::fusion::p_BoundInterfaceConstraint* _get_impl(mosek::fusion::BoundInterfaceConstraint * _inst){ return static_cast< mosek::fusion::p_BoundInterfaceConstraint* >(mosek::fusion::p_SliceConstraint::_get_impl(_inst)); }
      static mosek::fusion::p_BoundInterfaceConstraint * _get_impl(mosek::fusion::BoundInterfaceConstraint::t _inst) { return _get_impl(_inst.get()); }
      p_BoundInterfaceConstraint(BoundInterfaceConstraint * _pubthis);
      virtual ~p_BoundInterfaceConstraint() { /* std::cout << "~p_BoundInterfaceConstraint" << std::endl;*/ };
      bool islower{};

      virtual void destroy();

      static BoundInterfaceConstraint::t _new_BoundInterfaceConstraint(monty::rc_ptr< ::mosek::fusion::Model > _2852_m,std::shared_ptr< monty::ndarray< int32_t,1 > > _2853_shape,std::shared_ptr< monty::ndarray< int32_t,1 > > _2854_nativeidxs,bool _2855_islower);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Model > _2852_m,std::shared_ptr< monty::ndarray< int32_t,1 > > _2853_shape,std::shared_ptr< monty::ndarray< int32_t,1 > > _2854_nativeidxs,bool _2855_islower);
      static BoundInterfaceConstraint::t _new_BoundInterfaceConstraint(monty::rc_ptr< ::mosek::fusion::SliceConstraint > _2856_c,bool _2857_islower);
      void _initialize(monty::rc_ptr< ::mosek::fusion::SliceConstraint > _2856_c,bool _2857_islower);
      virtual /* override */ std::shared_ptr< monty::ndarray< double,1 > > dual() ;
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2BoundInterfaceConstraint__slice(std::shared_ptr< monty::ndarray< int32_t,1 > > _2859_firsta,std::shared_ptr< monty::ndarray< int32_t,1 > > _2860_lasta) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Constraint__slice(std::shared_ptr< monty::ndarray< int32_t,1 > > _2859_firsta,std::shared_ptr< monty::ndarray< int32_t,1 > > _2860_lasta) { return __mosek_2fusion_2BoundInterfaceConstraint__slice(_2859_firsta,_2860_lasta); }
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2BoundInterfaceConstraint__slice(int32_t _2862_first,int32_t _2863_last) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Constraint__slice(int32_t _2862_first,int32_t _2863_last) { return __mosek_2fusion_2BoundInterfaceConstraint__slice(_2862_first,_2863_last); }
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2BoundInterfaceConstraint__index(std::shared_ptr< monty::ndarray< int32_t,1 > > _2865_idxa) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Constraint__index(std::shared_ptr< monty::ndarray< int32_t,1 > > _2865_idxa) { return __mosek_2fusion_2BoundInterfaceConstraint__index(_2865_idxa); }
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2BoundInterfaceConstraint__index(int32_t _2867_idx) ;
      virtual monty::rc_ptr< ::mosek::fusion::Constraint > __mosek_2fusion_2Constraint__index(int32_t _2867_idx) { return __mosek_2fusion_2BoundInterfaceConstraint__index(_2867_idx); }
      virtual monty::rc_ptr< ::mosek::fusion::BoundInterfaceConstraint > __mosek_2fusion_2BoundInterfaceConstraint__from_(monty::rc_ptr< ::mosek::fusion::Constraint > _2869_c) ;
    }; // struct BoundInterfaceConstraint;

    struct p_ModelConstraint : public ::mosek::fusion::p_Constraint
    {
      ModelConstraint * _pubthis;
      static mosek::fusion::p_ModelConstraint* _get_impl(mosek::fusion::ModelConstraint * _inst){ return static_cast< mosek::fusion::p_ModelConstraint* >(mosek::fusion::p_Constraint::_get_impl(_inst)); }
      static mosek::fusion::p_ModelConstraint * _get_impl(mosek::fusion::ModelConstraint::t _inst) { return _get_impl(_inst.get()); }
      p_ModelConstraint(ModelConstraint * _pubthis);
      virtual ~p_ModelConstraint() { /* std::cout << "~p_ModelConstraint" << std::endl;*/ };
      int32_t conid{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > modelcon_nativeidxs{};
      std::string name{};

      virtual void destroy();

      static ModelConstraint::t _new_ModelConstraint(monty::rc_ptr< ::mosek::fusion::ModelConstraint > _2965_c,monty::rc_ptr< ::mosek::fusion::Model > _2966_m);
      void _initialize(monty::rc_ptr< ::mosek::fusion::ModelConstraint > _2965_c,monty::rc_ptr< ::mosek::fusion::Model > _2966_m);
      static ModelConstraint::t _new_ModelConstraint(monty::rc_ptr< ::mosek::fusion::Model > _2967_model,const std::string &  _2968_name,std::shared_ptr< monty::ndarray< int32_t,1 > > _2969_shape,std::shared_ptr< monty::ndarray< int32_t,1 > > _2970_nidxs,int32_t _2971_conid);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Model > _2967_model,const std::string &  _2968_name,std::shared_ptr< monty::ndarray< int32_t,1 > > _2969_shape,std::shared_ptr< monty::ndarray< int32_t,1 > > _2970_nidxs,int32_t _2971_conid);
      virtual /* override */ std::string toString() ;
      virtual monty::rc_ptr< ::mosek::fusion::ModelConstraint > __mosek_2fusion_2ModelConstraint__clone(monty::rc_ptr< ::mosek::fusion::Model > _2973_m) { throw monty::AbstractClassError("Call to abstract method"); }
      virtual /* override */ void remove() ;
    }; // struct ModelConstraint;

    struct p_LinearPSDConstraint : public ::mosek::fusion::p_ModelConstraint
    {
      LinearPSDConstraint * _pubthis;
      static mosek::fusion::p_LinearPSDConstraint* _get_impl(mosek::fusion::LinearPSDConstraint * _inst){ return static_cast< mosek::fusion::p_LinearPSDConstraint* >(mosek::fusion::p_ModelConstraint::_get_impl(_inst)); }
      static mosek::fusion::p_LinearPSDConstraint * _get_impl(mosek::fusion::LinearPSDConstraint::t _inst) { return _get_impl(_inst.get()); }
      p_LinearPSDConstraint(LinearPSDConstraint * _pubthis);
      virtual ~p_LinearPSDConstraint() { /* std::cout << "~p_LinearPSDConstraint" << std::endl;*/ };
      int32_t conedim{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};
      int32_t conid{};
      std::shared_ptr< monty::ndarray< int64_t,1 > > slackidxs{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > nativeidxs{};

      virtual void destroy();

      static LinearPSDConstraint::t _new_LinearPSDConstraint(monty::rc_ptr< ::mosek::fusion::LinearPSDConstraint > _2872_c,monty::rc_ptr< ::mosek::fusion::Model > _2873_m);
      void _initialize(monty::rc_ptr< ::mosek::fusion::LinearPSDConstraint > _2872_c,monty::rc_ptr< ::mosek::fusion::Model > _2873_m);
      static LinearPSDConstraint::t _new_LinearPSDConstraint(monty::rc_ptr< ::mosek::fusion::Model > _2874_model,const std::string &  _2875_name,int32_t _2876_conid,std::shared_ptr< monty::ndarray< int32_t,1 > > _2877_shape,int32_t _2878_conedim,std::shared_ptr< monty::ndarray< int32_t,1 > > _2879_nativeidxs,std::shared_ptr< monty::ndarray< int64_t,1 > > _2880_slackidxs);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Model > _2874_model,const std::string &  _2875_name,int32_t _2876_conid,std::shared_ptr< monty::ndarray< int32_t,1 > > _2877_shape,int32_t _2878_conedim,std::shared_ptr< monty::ndarray< int32_t,1 > > _2879_nativeidxs,std::shared_ptr< monty::ndarray< int64_t,1 > > _2880_slackidxs);
      virtual void domainToString(int64_t _2881_i,monty::rc_ptr< ::mosek::fusion::Utils::StringBuffer > _2882_sb) ;
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::ModelConstraint > __mosek_2fusion_2LinearPSDConstraint__clone(monty::rc_ptr< ::mosek::fusion::Model > _2886_m) ;
      virtual monty::rc_ptr< ::mosek::fusion::ModelConstraint > __mosek_2fusion_2ModelConstraint__clone(monty::rc_ptr< ::mosek::fusion::Model > _2886_m) { return __mosek_2fusion_2LinearPSDConstraint__clone(_2886_m); }
    }; // struct LinearPSDConstraint;

    struct p_PSDConstraint : public ::mosek::fusion::p_ModelConstraint
    {
      PSDConstraint * _pubthis;
      static mosek::fusion::p_PSDConstraint* _get_impl(mosek::fusion::PSDConstraint * _inst){ return static_cast< mosek::fusion::p_PSDConstraint* >(mosek::fusion::p_ModelConstraint::_get_impl(_inst)); }
      static mosek::fusion::p_PSDConstraint * _get_impl(mosek::fusion::PSDConstraint::t _inst) { return _get_impl(_inst.get()); }
      p_PSDConstraint(PSDConstraint * _pubthis);
      virtual ~p_PSDConstraint() { /* std::cout << "~p_PSDConstraint" << std::endl;*/ };
      int32_t conedim1{};
      int32_t conedim0{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};
      std::string name{};
      std::shared_ptr< monty::ndarray< int64_t,1 > > slackidxs{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > nativeidxs{};
      int32_t conid{};

      virtual void destroy();

      static PSDConstraint::t _new_PSDConstraint(monty::rc_ptr< ::mosek::fusion::PSDConstraint > _2887_c,monty::rc_ptr< ::mosek::fusion::Model > _2888_m);
      void _initialize(monty::rc_ptr< ::mosek::fusion::PSDConstraint > _2887_c,monty::rc_ptr< ::mosek::fusion::Model > _2888_m);
      static PSDConstraint::t _new_PSDConstraint(monty::rc_ptr< ::mosek::fusion::Model > _2889_model,const std::string &  _2890_name,int32_t _2891_conid,std::shared_ptr< monty::ndarray< int32_t,1 > > _2892_shape,int32_t _2893_conedim0,int32_t _2894_conedim1,std::shared_ptr< monty::ndarray< int64_t,1 > > _2895_slackidxs,std::shared_ptr< monty::ndarray< int32_t,1 > > _2896_nativeidxs);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Model > _2889_model,const std::string &  _2890_name,int32_t _2891_conid,std::shared_ptr< monty::ndarray< int32_t,1 > > _2892_shape,int32_t _2893_conedim0,int32_t _2894_conedim1,std::shared_ptr< monty::ndarray< int64_t,1 > > _2895_slackidxs,std::shared_ptr< monty::ndarray< int32_t,1 > > _2896_nativeidxs);
      virtual /* override */ std::string toString() ;
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::ModelConstraint > __mosek_2fusion_2PSDConstraint__clone(monty::rc_ptr< ::mosek::fusion::Model > _2897_m) ;
      virtual monty::rc_ptr< ::mosek::fusion::ModelConstraint > __mosek_2fusion_2ModelConstraint__clone(monty::rc_ptr< ::mosek::fusion::Model > _2897_m) { return __mosek_2fusion_2PSDConstraint__clone(_2897_m); }
      static  std::shared_ptr< monty::ndarray< int32_t,1 > > computenidxs(std::shared_ptr< monty::ndarray< int32_t,1 > > _2898_shape,int32_t _2899_d0,int32_t _2900_d1,std::shared_ptr< monty::ndarray< int32_t,1 > > _2901_nativeidxs);
    }; // struct PSDConstraint;

    struct p_RangedConstraint : public ::mosek::fusion::p_ModelConstraint
    {
      RangedConstraint * _pubthis;
      static mosek::fusion::p_RangedConstraint* _get_impl(mosek::fusion::RangedConstraint * _inst){ return static_cast< mosek::fusion::p_RangedConstraint* >(mosek::fusion::p_ModelConstraint::_get_impl(_inst)); }
      static mosek::fusion::p_RangedConstraint * _get_impl(mosek::fusion::RangedConstraint::t _inst) { return _get_impl(_inst.get()); }
      p_RangedConstraint(RangedConstraint * _pubthis);
      virtual ~p_RangedConstraint() { /* std::cout << "~p_RangedConstraint" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< int32_t,1 > > nativeidxs{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};

      virtual void destroy();

      static RangedConstraint::t _new_RangedConstraint(monty::rc_ptr< ::mosek::fusion::RangedConstraint > _2931_c,monty::rc_ptr< ::mosek::fusion::Model > _2932_m);
      void _initialize(monty::rc_ptr< ::mosek::fusion::RangedConstraint > _2931_c,monty::rc_ptr< ::mosek::fusion::Model > _2932_m);
      static RangedConstraint::t _new_RangedConstraint(monty::rc_ptr< ::mosek::fusion::Model > _2933_model,const std::string &  _2934_name,std::shared_ptr< monty::ndarray< int32_t,1 > > _2935_shape,std::shared_ptr< monty::ndarray< int32_t,1 > > _2936_nativeidxs,int32_t _2937_conid);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Model > _2933_model,const std::string &  _2934_name,std::shared_ptr< monty::ndarray< int32_t,1 > > _2935_shape,std::shared_ptr< monty::ndarray< int32_t,1 > > _2936_nativeidxs,int32_t _2937_conid);
      virtual monty::rc_ptr< ::mosek::fusion::BoundInterfaceConstraint > __mosek_2fusion_2RangedConstraint__upperBoundCon() ;
      virtual monty::rc_ptr< ::mosek::fusion::BoundInterfaceConstraint > __mosek_2fusion_2RangedConstraint__lowerBoundCon() ;
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::ModelConstraint > __mosek_2fusion_2RangedConstraint__clone(monty::rc_ptr< ::mosek::fusion::Model > _2938_m) ;
      virtual monty::rc_ptr< ::mosek::fusion::ModelConstraint > __mosek_2fusion_2ModelConstraint__clone(monty::rc_ptr< ::mosek::fusion::Model > _2938_m) { return __mosek_2fusion_2RangedConstraint__clone(_2938_m); }
    }; // struct RangedConstraint;

    struct p_ConicConstraint : public ::mosek::fusion::p_ModelConstraint
    {
      ConicConstraint * _pubthis;
      static mosek::fusion::p_ConicConstraint* _get_impl(mosek::fusion::ConicConstraint * _inst){ return static_cast< mosek::fusion::p_ConicConstraint* >(mosek::fusion::p_ModelConstraint::_get_impl(_inst)); }
      static mosek::fusion::p_ConicConstraint * _get_impl(mosek::fusion::ConicConstraint::t _inst) { return _get_impl(_inst.get()); }
      p_ConicConstraint(ConicConstraint * _pubthis);
      virtual ~p_ConicConstraint() { /* std::cout << "~p_ConicConstraint" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< std::string,1 > >,1 > > indexnames{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > nativeidxs{};
      bool names_flushed{};
      std::string name{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};
      monty::rc_ptr< ::mosek::fusion::ConeDomain > dom{};
      int32_t conid{};

      virtual void destroy();

      static ConicConstraint::t _new_ConicConstraint(monty::rc_ptr< ::mosek::fusion::ConicConstraint > _2939_c,monty::rc_ptr< ::mosek::fusion::Model > _2940_m);
      void _initialize(monty::rc_ptr< ::mosek::fusion::ConicConstraint > _2939_c,monty::rc_ptr< ::mosek::fusion::Model > _2940_m);
      static ConicConstraint::t _new_ConicConstraint(monty::rc_ptr< ::mosek::fusion::Model > _2941_model,const std::string &  _2942_name,monty::rc_ptr< ::mosek::fusion::ConeDomain > _2943_dom,std::shared_ptr< monty::ndarray< int32_t,1 > > _2944_shape,int32_t _2945_conid,std::shared_ptr< monty::ndarray< int32_t,1 > > _2946_nativeidxs,std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< std::string,1 > >,1 > > _2947_indexnames);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Model > _2941_model,const std::string &  _2942_name,monty::rc_ptr< ::mosek::fusion::ConeDomain > _2943_dom,std::shared_ptr< monty::ndarray< int32_t,1 > > _2944_shape,int32_t _2945_conid,std::shared_ptr< monty::ndarray< int32_t,1 > > _2946_nativeidxs,std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< std::string,1 > >,1 > > _2947_indexnames);
      virtual /* override */ std::string toString() ;
      virtual void domainToString(int64_t _2950_i,monty::rc_ptr< ::mosek::fusion::Utils::StringBuffer > _2951_sb) ;
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::ModelConstraint > __mosek_2fusion_2ConicConstraint__clone(monty::rc_ptr< ::mosek::fusion::Model > _2952_m) ;
      virtual monty::rc_ptr< ::mosek::fusion::ModelConstraint > __mosek_2fusion_2ModelConstraint__clone(monty::rc_ptr< ::mosek::fusion::Model > _2952_m) { return __mosek_2fusion_2ConicConstraint__clone(_2952_m); }
    }; // struct ConicConstraint;

    struct p_LinearConstraint : public ::mosek::fusion::p_ModelConstraint
    {
      LinearConstraint * _pubthis;
      static mosek::fusion::p_LinearConstraint* _get_impl(mosek::fusion::LinearConstraint * _inst){ return static_cast< mosek::fusion::p_LinearConstraint* >(mosek::fusion::p_ModelConstraint::_get_impl(_inst)); }
      static mosek::fusion::p_LinearConstraint * _get_impl(mosek::fusion::LinearConstraint::t _inst) { return _get_impl(_inst.get()); }
      p_LinearConstraint(LinearConstraint * _pubthis);
      virtual ~p_LinearConstraint() { /* std::cout << "~p_LinearConstraint" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< std::string,1 > >,1 > > indexnames{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > nidxs{};
      std::string name{};
      int32_t conid{};

      virtual void destroy();

      static LinearConstraint::t _new_LinearConstraint(monty::rc_ptr< ::mosek::fusion::LinearConstraint > _2953_c,monty::rc_ptr< ::mosek::fusion::Model > _2954_m);
      void _initialize(monty::rc_ptr< ::mosek::fusion::LinearConstraint > _2953_c,monty::rc_ptr< ::mosek::fusion::Model > _2954_m);
      static LinearConstraint::t _new_LinearConstraint(monty::rc_ptr< ::mosek::fusion::Model > _2955_model,const std::string &  _2956_name,int32_t _2957_conid,std::shared_ptr< monty::ndarray< int32_t,1 > > _2958_shape,std::shared_ptr< monty::ndarray< int32_t,1 > > _2959_nidxs,std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< std::string,1 > >,1 > > _2960_indexnames);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Model > _2955_model,const std::string &  _2956_name,int32_t _2957_conid,std::shared_ptr< monty::ndarray< int32_t,1 > > _2958_shape,std::shared_ptr< monty::ndarray< int32_t,1 > > _2959_nidxs,std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< std::string,1 > >,1 > > _2960_indexnames);
      virtual /* override */ std::string toString() ;
      virtual void domainToString(int64_t _2962_i,monty::rc_ptr< ::mosek::fusion::Utils::StringBuffer > _2963_sb) ;
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::ModelConstraint > __mosek_2fusion_2LinearConstraint__clone(monty::rc_ptr< ::mosek::fusion::Model > _2964_m) ;
      virtual monty::rc_ptr< ::mosek::fusion::ModelConstraint > __mosek_2fusion_2ModelConstraint__clone(monty::rc_ptr< ::mosek::fusion::Model > _2964_m) { return __mosek_2fusion_2LinearConstraint__clone(_2964_m); }
    }; // struct LinearConstraint;

    struct p_Set
    {
      Set * _pubthis;
      static mosek::fusion::p_Set* _get_impl(mosek::fusion::Set * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_Set * _get_impl(mosek::fusion::Set::t _inst) { return _get_impl(_inst.get()); }
      p_Set(Set * _pubthis);
      virtual ~p_Set() { /* std::cout << "~p_Set" << std::endl;*/ };

      virtual void destroy();

      static  int64_t size(std::shared_ptr< monty::ndarray< int32_t,1 > > _3128_shape);
      static  bool match(std::shared_ptr< monty::ndarray< int32_t,1 > > _3131_s1,std::shared_ptr< monty::ndarray< int32_t,1 > > _3132_s2);
      static  int64_t linearidx(std::shared_ptr< monty::ndarray< int32_t,1 > > _3134_shape,std::shared_ptr< monty::ndarray< int32_t,1 > > _3135_key);
      static  std::shared_ptr< monty::ndarray< int32_t,1 > > idxtokey(std::shared_ptr< monty::ndarray< int32_t,1 > > _3138_shape,int64_t _3139_idx);
      static  void idxtokey(std::shared_ptr< monty::ndarray< int32_t,1 > > _3141_shape,int64_t _3142_idx,std::shared_ptr< monty::ndarray< int32_t,1 > > _3143_dest);
      static  std::string indexToString(std::shared_ptr< monty::ndarray< int32_t,1 > > _3147_shape,int64_t _3148_key);
      static  std::string keyToString(std::shared_ptr< monty::ndarray< int32_t,1 > > _3155_key);
      static  void indexToKey(std::shared_ptr< monty::ndarray< int32_t,1 > > _3158_shape,int64_t _3159_key,std::shared_ptr< monty::ndarray< int32_t,1 > > _3160_res);
      static  std::shared_ptr< monty::ndarray< int64_t,1 > > strides(std::shared_ptr< monty::ndarray< int32_t,1 > > _3164_shape);
      static  std::shared_ptr< monty::ndarray< int32_t,1 > > make(std::shared_ptr< monty::ndarray< int32_t,1 > > _3168_set1,std::shared_ptr< monty::ndarray< int32_t,1 > > _3169_set2);
      static  std::shared_ptr< monty::ndarray< int32_t,1 > > make(std::shared_ptr< monty::ndarray< int32_t,1 > > _3173_sizes);
      static  std::shared_ptr< monty::ndarray< int32_t,1 > > make(int32_t _3175_s1,int32_t _3176_s2,int32_t _3177_s3);
      static  std::shared_ptr< monty::ndarray< int32_t,1 > > make(int32_t _3178_s1,int32_t _3179_s2);
      static  std::shared_ptr< monty::ndarray< int32_t,1 > > make(int32_t _3180_sz);
      static  std::shared_ptr< monty::ndarray< int32_t,1 > > scalar();
      static  std::shared_ptr< monty::ndarray< int32_t,1 > > make(std::shared_ptr< monty::ndarray< std::string,1 > > _3181_names);
    }; // struct Set;

    struct p_ConeDomain
    {
      ConeDomain * _pubthis;
      static mosek::fusion::p_ConeDomain* _get_impl(mosek::fusion::ConeDomain * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_ConeDomain * _get_impl(mosek::fusion::ConeDomain::t _inst) { return _get_impl(_inst.get()); }
      p_ConeDomain(ConeDomain * _pubthis);
      virtual ~p_ConeDomain() { /* std::cout << "~p_ConeDomain" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< std::string,1 > >,1 > > indexnames{};
      int64_t domsize{};
      std::shared_ptr< monty::ndarray< double,1 > > domofs{};
      std::shared_ptr< monty::ndarray< double,1 > > alpha{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};
      bool int_flag{};
      bool axisset{};
      int32_t axisidx{};
      mosek::fusion::QConeKey key{};

      virtual void destroy();

      static ConeDomain::t _new_ConeDomain(mosek::fusion::QConeKey _3182_k,std::shared_ptr< monty::ndarray< double,1 > > _3183_alpha,std::shared_ptr< monty::ndarray< int32_t,1 > > _3184_d);
      void _initialize(mosek::fusion::QConeKey _3182_k,std::shared_ptr< monty::ndarray< double,1 > > _3183_alpha,std::shared_ptr< monty::ndarray< int32_t,1 > > _3184_d);
      static ConeDomain::t _new_ConeDomain(mosek::fusion::QConeKey _3185_k,std::shared_ptr< monty::ndarray< int32_t,1 > > _3186_d);
      void _initialize(mosek::fusion::QConeKey _3185_k,std::shared_ptr< monty::ndarray< int32_t,1 > > _3186_d);
      static ConeDomain::t _new_ConeDomain(monty::rc_ptr< ::mosek::fusion::ConeDomain > _3187_other);
      void _initialize(monty::rc_ptr< ::mosek::fusion::ConeDomain > _3187_other);
      virtual bool match_shape(std::shared_ptr< monty::ndarray< int32_t,1 > > _3188_shp) ;
      virtual monty::rc_ptr< ::mosek::fusion::ConeDomain > __mosek_2fusion_2ConeDomain__integral() ;
      virtual bool axisIsSet() ;
      virtual int32_t getAxis() ;
      virtual monty::rc_ptr< ::mosek::fusion::ConeDomain > __mosek_2fusion_2ConeDomain__axis(int32_t _3189_a) ;
      virtual monty::rc_ptr< ::mosek::fusion::ConeDomain > __mosek_2fusion_2ConeDomain__withShape(int32_t _3190_dim0,int32_t _3191_dim1,int32_t _3192_dim2) ;
      virtual monty::rc_ptr< ::mosek::fusion::ConeDomain > __mosek_2fusion_2ConeDomain__withShape(int32_t _3193_dim0,int32_t _3194_dim1) ;
      virtual monty::rc_ptr< ::mosek::fusion::ConeDomain > __mosek_2fusion_2ConeDomain__withShape(int32_t _3195_dim0) ;
      virtual monty::rc_ptr< ::mosek::fusion::ConeDomain > __mosek_2fusion_2ConeDomain__withShape(std::shared_ptr< monty::ndarray< int32_t,1 > > _3196_shp) ;
      virtual monty::rc_ptr< ::mosek::fusion::ConeDomain > __mosek_2fusion_2ConeDomain__withShape_(std::shared_ptr< monty::ndarray< int32_t,1 > > _3197_shp) ;
      virtual monty::rc_ptr< ::mosek::fusion::ConeDomain > __mosek_2fusion_2ConeDomain__withNamesOnAxis(std::shared_ptr< monty::ndarray< std::string,1 > > _3198_names,int32_t _3199_axis) ;
      virtual void finalize_and_validate_inplace(std::shared_ptr< monty::ndarray< int32_t,1 > > _3204_shp) ;
      virtual monty::rc_ptr< ::mosek::fusion::ConeDomain > __mosek_2fusion_2ConeDomain__finalize_and_validate(std::shared_ptr< monty::ndarray< int32_t,1 > > _3208_shp) ;
    }; // struct ConeDomain;

    struct p_PSDDomain
    {
      PSDDomain * _pubthis;
      static mosek::fusion::p_PSDDomain* _get_impl(mosek::fusion::PSDDomain * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_PSDDomain * _get_impl(mosek::fusion::PSDDomain::t _inst) { return _get_impl(_inst.get()); }
      p_PSDDomain(PSDDomain * _pubthis);
      virtual ~p_PSDDomain() { /* std::cout << "~p_PSDDomain" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< std::string,1 > >,1 > > indexnames{};
      bool axisIsSet{};
      int32_t conedim2{};
      int32_t conedim1{};
      mosek::fusion::PSDKey key{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};

      virtual void destroy();

      static PSDDomain::t _new_PSDDomain(mosek::fusion::PSDKey _3210_k,std::shared_ptr< monty::ndarray< int32_t,1 > > _3211_shp,int32_t _3212_conedim1,int32_t _3213_conedim2);
      void _initialize(mosek::fusion::PSDKey _3210_k,std::shared_ptr< monty::ndarray< int32_t,1 > > _3211_shp,int32_t _3212_conedim1,int32_t _3213_conedim2);
      static PSDDomain::t _new_PSDDomain(mosek::fusion::PSDKey _3215_k,std::shared_ptr< monty::ndarray< int32_t,1 > > _3216_shp);
      void _initialize(mosek::fusion::PSDKey _3215_k,std::shared_ptr< monty::ndarray< int32_t,1 > > _3216_shp);
      static PSDDomain::t _new_PSDDomain(mosek::fusion::PSDKey _3217_k);
      void _initialize(mosek::fusion::PSDKey _3217_k);
      static PSDDomain::t _new_PSDDomain(monty::rc_ptr< ::mosek::fusion::PSDDomain > _3218_other);
      void _initialize(monty::rc_ptr< ::mosek::fusion::PSDDomain > _3218_other);
      virtual monty::rc_ptr< ::mosek::fusion::PSDDomain > __mosek_2fusion_2PSDDomain__axis(int32_t _3219_conedim1,int32_t _3220_conedim2) ;
      virtual monty::rc_ptr< ::mosek::fusion::PSDDomain > __mosek_2fusion_2PSDDomain__withNamesOnAxis(std::shared_ptr< monty::ndarray< std::string,1 > > _3221_names,int32_t _3222_axis) ;
      virtual void finalize_and_validate_inplace(std::shared_ptr< monty::ndarray< int32_t,1 > > _3229_shp) ;
      virtual monty::rc_ptr< ::mosek::fusion::PSDDomain > __mosek_2fusion_2PSDDomain__finalize_and_validate(std::shared_ptr< monty::ndarray< int32_t,1 > > _3232_shp) ;
    }; // struct PSDDomain;

    struct p_RangeDomain
    {
      RangeDomain * _pubthis;
      static mosek::fusion::p_RangeDomain* _get_impl(mosek::fusion::RangeDomain * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_RangeDomain * _get_impl(mosek::fusion::RangeDomain::t _inst) { return _get_impl(_inst.get()); }
      p_RangeDomain(RangeDomain * _pubthis);
      virtual ~p_RangeDomain() { /* std::cout << "~p_RangeDomain" << std::endl;*/ };
      int64_t domsize{};
      int64_t nelements{};
      std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< std::string,1 > >,1 > > indexnames{};
      bool cardinal_flag{};
      bool scalable{};
      std::shared_ptr< monty::ndarray< double,1 > > ub{};
      std::shared_ptr< monty::ndarray< double,1 > > lb{};
      std::shared_ptr< monty::ndarray< int32_t,2 > > sparsity{};
      bool empty{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};

      virtual void destroy();

      static RangeDomain::t _new_RangeDomain(bool _3234_scalable,std::shared_ptr< monty::ndarray< double,1 > > _3235_lb,std::shared_ptr< monty::ndarray< double,1 > > _3236_ub,std::shared_ptr< monty::ndarray< int32_t,1 > > _3237_dims);
      void _initialize(bool _3234_scalable,std::shared_ptr< monty::ndarray< double,1 > > _3235_lb,std::shared_ptr< monty::ndarray< double,1 > > _3236_ub,std::shared_ptr< monty::ndarray< int32_t,1 > > _3237_dims);
      static RangeDomain::t _new_RangeDomain(bool _3238_scalable,double _3239_lb,double _3240_ub,std::shared_ptr< monty::ndarray< int32_t,1 > > _3241_dims);
      void _initialize(bool _3238_scalable,double _3239_lb,double _3240_ub,std::shared_ptr< monty::ndarray< int32_t,1 > > _3241_dims);
      static RangeDomain::t _new_RangeDomain(bool _3246_scalable,std::shared_ptr< monty::ndarray< double,1 > > _3247_lb,std::shared_ptr< monty::ndarray< double,1 > > _3248_ub,std::shared_ptr< monty::ndarray< int32_t,1 > > _3249_dims,std::shared_ptr< monty::ndarray< int32_t,2 > > _3250_sp);
      void _initialize(bool _3246_scalable,std::shared_ptr< monty::ndarray< double,1 > > _3247_lb,std::shared_ptr< monty::ndarray< double,1 > > _3248_ub,std::shared_ptr< monty::ndarray< int32_t,1 > > _3249_dims,std::shared_ptr< monty::ndarray< int32_t,2 > > _3250_sp);
      static RangeDomain::t _new_RangeDomain(bool _3251_scalable,std::shared_ptr< monty::ndarray< double,1 > > _3252_lb,std::shared_ptr< monty::ndarray< double,1 > > _3253_ub,std::shared_ptr< monty::ndarray< int32_t,1 > > _3254_dims,std::shared_ptr< monty::ndarray< int32_t,2 > > _3255_sp,int32_t _3256_steal);
      void _initialize(bool _3251_scalable,std::shared_ptr< monty::ndarray< double,1 > > _3252_lb,std::shared_ptr< monty::ndarray< double,1 > > _3253_ub,std::shared_ptr< monty::ndarray< int32_t,1 > > _3254_dims,std::shared_ptr< monty::ndarray< int32_t,2 > > _3255_sp,int32_t _3256_steal);
      static RangeDomain::t _new_RangeDomain(monty::rc_ptr< ::mosek::fusion::RangeDomain > _3257_other);
      void _initialize(monty::rc_ptr< ::mosek::fusion::RangeDomain > _3257_other);
      virtual monty::rc_ptr< ::mosek::fusion::RangeDomain > __mosek_2fusion_2RangeDomain__sparse(std::shared_ptr< monty::ndarray< int32_t,2 > > _3259_sparsity) ;
      virtual monty::rc_ptr< ::mosek::fusion::RangeDomain > __mosek_2fusion_2RangeDomain__sparse(std::shared_ptr< monty::ndarray< int32_t,1 > > _3262_sparsity) ;
      virtual monty::rc_ptr< ::mosek::fusion::RangeDomain > __mosek_2fusion_2RangeDomain__sparse() ;
      virtual monty::rc_ptr< ::mosek::fusion::RangeDomain > __mosek_2fusion_2RangeDomain__integral() ;
      virtual monty::rc_ptr< ::mosek::fusion::RangeDomain > __mosek_2fusion_2RangeDomain__withShape(int32_t _3264_dim0,int32_t _3265_dim1,int32_t _3266_dim2) ;
      virtual monty::rc_ptr< ::mosek::fusion::RangeDomain > __mosek_2fusion_2RangeDomain__withShape(int32_t _3267_dim0,int32_t _3268_dim1) ;
      virtual monty::rc_ptr< ::mosek::fusion::RangeDomain > __mosek_2fusion_2RangeDomain__withShape(int32_t _3269_dim0) ;
      virtual monty::rc_ptr< ::mosek::fusion::RangeDomain > __mosek_2fusion_2RangeDomain__withShape(std::shared_ptr< monty::ndarray< int32_t,1 > > _3270_shp) ;
      virtual monty::rc_ptr< ::mosek::fusion::RangeDomain > __mosek_2fusion_2RangeDomain__withNamesOnAxis(std::shared_ptr< monty::ndarray< std::string,1 > > _3271_names,int32_t _3272_axis) ;
      virtual bool match_shape(std::shared_ptr< monty::ndarray< int32_t,1 > > _3279_shp) ;
      virtual void finalize_and_validate_inplace(std::shared_ptr< monty::ndarray< int32_t,1 > > _3281_shp) ;
      virtual monty::rc_ptr< ::mosek::fusion::RangeDomain > __mosek_2fusion_2RangeDomain__finalize_and_validate(std::shared_ptr< monty::ndarray< int32_t,1 > > _3288_shp) ;
    }; // struct RangeDomain;

    struct p_LinearDomain
    {
      LinearDomain * _pubthis;
      static mosek::fusion::p_LinearDomain* _get_impl(mosek::fusion::LinearDomain * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_LinearDomain * _get_impl(mosek::fusion::LinearDomain::t _inst) { return _get_impl(_inst.get()); }
      p_LinearDomain(LinearDomain * _pubthis);
      virtual ~p_LinearDomain() { /* std::cout << "~p_LinearDomain" << std::endl;*/ };
      int64_t nelements{};
      int64_t domsize{};
      std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< std::string,1 > >,1 > > indexnames{};
      bool empty{};
      bool scalable{};
      std::shared_ptr< monty::ndarray< int32_t,2 > > sparsity{};
      bool cardinal_flag{};
      mosek::fusion::RelationKey key{};
      std::shared_ptr< monty::ndarray< double,1 > > bnd{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};

      virtual void destroy();

      static LinearDomain::t _new_LinearDomain(mosek::fusion::RelationKey _3290_k,bool _3291_scalable,std::shared_ptr< monty::ndarray< double,1 > > _3292_rhs,std::shared_ptr< monty::ndarray< int32_t,1 > > _3293_dims);
      void _initialize(mosek::fusion::RelationKey _3290_k,bool _3291_scalable,std::shared_ptr< monty::ndarray< double,1 > > _3292_rhs,std::shared_ptr< monty::ndarray< int32_t,1 > > _3293_dims);
      static LinearDomain::t _new_LinearDomain(mosek::fusion::RelationKey _3294_k,bool _3295_scalable,std::shared_ptr< monty::ndarray< double,1 > > _3296_rhs,std::shared_ptr< monty::ndarray< int32_t,1 > > _3297_dims,std::shared_ptr< monty::ndarray< int32_t,2 > > _3298_sp,int32_t _3299_steal);
      void _initialize(mosek::fusion::RelationKey _3294_k,bool _3295_scalable,std::shared_ptr< monty::ndarray< double,1 > > _3296_rhs,std::shared_ptr< monty::ndarray< int32_t,1 > > _3297_dims,std::shared_ptr< monty::ndarray< int32_t,2 > > _3298_sp,int32_t _3299_steal);
      static LinearDomain::t _new_LinearDomain(monty::rc_ptr< ::mosek::fusion::LinearDomain > _3300_other);
      void _initialize(monty::rc_ptr< ::mosek::fusion::LinearDomain > _3300_other);
      virtual monty::rc_ptr< ::mosek::fusion::LinearDomain > __mosek_2fusion_2LinearDomain__sparse(std::shared_ptr< monty::ndarray< int32_t,2 > > _3301_sparsity) ;
      virtual monty::rc_ptr< ::mosek::fusion::LinearDomain > __mosek_2fusion_2LinearDomain__sparse(std::shared_ptr< monty::ndarray< int32_t,1 > > _3304_sparsity) ;
      virtual monty::rc_ptr< ::mosek::fusion::LinearDomain > __mosek_2fusion_2LinearDomain__sparse() ;
      virtual monty::rc_ptr< ::mosek::fusion::LinearDomain > __mosek_2fusion_2LinearDomain__integral() ;
      virtual monty::rc_ptr< ::mosek::fusion::LinearDomain > __mosek_2fusion_2LinearDomain__withShape(int32_t _3306_dim0,int32_t _3307_dim1,int32_t _3308_dim2) ;
      virtual monty::rc_ptr< ::mosek::fusion::LinearDomain > __mosek_2fusion_2LinearDomain__withShape(int32_t _3309_dim0,int32_t _3310_dim1) ;
      virtual monty::rc_ptr< ::mosek::fusion::LinearDomain > __mosek_2fusion_2LinearDomain__withShape(int32_t _3311_dim0) ;
      virtual monty::rc_ptr< ::mosek::fusion::LinearDomain > __mosek_2fusion_2LinearDomain__withShape(std::shared_ptr< monty::ndarray< int32_t,1 > > _3312_shp) ;
      virtual monty::rc_ptr< ::mosek::fusion::LinearDomain > __mosek_2fusion_2LinearDomain__withNamesOnAxis(std::shared_ptr< monty::ndarray< std::string,1 > > _3313_names,int32_t _3314_axis) ;
      virtual bool match_shape(std::shared_ptr< monty::ndarray< int32_t,1 > > _3321_shp) ;
      virtual void finalize_and_validate_inplace(std::shared_ptr< monty::ndarray< int32_t,1 > > _3323_shp) ;
      virtual monty::rc_ptr< ::mosek::fusion::LinearDomain > __mosek_2fusion_2LinearDomain__finalize_and_validate(std::shared_ptr< monty::ndarray< int32_t,1 > > _3331_shp) ;
    }; // struct LinearDomain;

    struct p_Domain
    {
      Domain * _pubthis;
      static mosek::fusion::p_Domain* _get_impl(mosek::fusion::Domain * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_Domain * _get_impl(mosek::fusion::Domain::t _inst) { return _get_impl(_inst.get()); }
      p_Domain(Domain * _pubthis);
      virtual ~p_Domain() { /* std::cout << "~p_Domain" << std::endl;*/ };

      virtual void destroy();

      static  int64_t dimsize(std::shared_ptr< monty::ndarray< int32_t,1 > > _3333_dims);
      static  monty::rc_ptr< ::mosek::fusion::RangeDomain > mkRangedDomain(monty::rc_ptr< ::mosek::fusion::Matrix > _3336_lb,monty::rc_ptr< ::mosek::fusion::Matrix > _3337_ub);
      static  monty::rc_ptr< ::mosek::fusion::RangeDomain > mkRangedDomain(std::shared_ptr< monty::ndarray< double,2 > > _3366_lb,std::shared_ptr< monty::ndarray< double,2 > > _3367_ub);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > mkLinearDomain(mosek::fusion::RelationKey _3376_k,monty::rc_ptr< ::mosek::fusion::Matrix > _3377_mx);
      static  int64_t prod(std::shared_ptr< monty::ndarray< int32_t,1 > > _3383_dim);
      static  monty::rc_ptr< ::mosek::fusion::RangeDomain > inRange(bool _3386_scalable,std::shared_ptr< monty::ndarray< double,1 > > _3387_lb,std::shared_ptr< monty::ndarray< double,1 > > _3388_ub,std::shared_ptr< monty::ndarray< int32_t,2 > > _3389_sp,std::shared_ptr< monty::ndarray< int32_t,1 > > _3390_dims);
      static  monty::rc_ptr< ::mosek::fusion::RangeDomain > sparse(monty::rc_ptr< ::mosek::fusion::RangeDomain > _3392_rd,std::shared_ptr< monty::ndarray< int32_t,2 > > _3393_sparsity);
      static  monty::rc_ptr< ::mosek::fusion::RangeDomain > sparse(monty::rc_ptr< ::mosek::fusion::RangeDomain > _3394_rd,std::shared_ptr< monty::ndarray< int32_t,1 > > _3395_sparsity);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > sparse(monty::rc_ptr< ::mosek::fusion::LinearDomain > _3396_ld,std::shared_ptr< monty::ndarray< int32_t,2 > > _3397_sparsity);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > sparse(monty::rc_ptr< ::mosek::fusion::LinearDomain > _3398_ld,std::shared_ptr< monty::ndarray< int32_t,1 > > _3399_sparsity);
      static  monty::rc_ptr< ::mosek::fusion::RangeDomain > integral(monty::rc_ptr< ::mosek::fusion::RangeDomain > _3400_rd);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > integral(monty::rc_ptr< ::mosek::fusion::LinearDomain > _3401_ld);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > integral(monty::rc_ptr< ::mosek::fusion::ConeDomain > _3402_c);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > axis(monty::rc_ptr< ::mosek::fusion::ConeDomain > _3403_c,int32_t _3404_a);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inDPowerCone(std::shared_ptr< monty::ndarray< double,1 > > _3405_alphas,std::shared_ptr< monty::ndarray< int32_t,1 > > _3406_dims);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inDPowerCone(std::shared_ptr< monty::ndarray< double,1 > > _3408_alphas,int32_t _3409_m);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inDPowerCone(std::shared_ptr< monty::ndarray< double,1 > > _3410_alphas);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inDPowerCone(double _3411_alpha,std::shared_ptr< monty::ndarray< int32_t,1 > > _3412_dims);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inDPowerCone(double _3414_alpha,int32_t _3415_m);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inDPowerCone(double _3416_alpha);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inPPowerCone(std::shared_ptr< monty::ndarray< double,1 > > _3417_alphas,std::shared_ptr< monty::ndarray< int32_t,1 > > _3418_dims);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inPPowerCone(std::shared_ptr< monty::ndarray< double,1 > > _3420_alphas,int32_t _3421_m);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inPPowerCone(std::shared_ptr< monty::ndarray< double,1 > > _3422_alphas);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inPPowerCone(double _3423_alpha,std::shared_ptr< monty::ndarray< int32_t,1 > > _3424_dims);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inPPowerCone(double _3426_alpha,int32_t _3427_m);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inPPowerCone(double _3428_alpha);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inDExpCone(std::shared_ptr< monty::ndarray< int32_t,1 > > _3429_dims);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inDExpCone(int32_t _3431_m);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inDExpCone();
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inPExpCone(std::shared_ptr< monty::ndarray< int32_t,1 > > _3432_dims);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inPExpCone(int32_t _3434_m);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inPExpCone();
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inDGeoMeanCone(std::shared_ptr< monty::ndarray< int32_t,1 > > _3435_dims);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inDGeoMeanCone(int32_t _3437_m,int32_t _3438_n);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inDGeoMeanCone(int32_t _3439_n);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inDGeoMeanCone();
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inPGeoMeanCone(std::shared_ptr< monty::ndarray< int32_t,1 > > _3440_dims);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inPGeoMeanCone(int32_t _3442_m,int32_t _3443_n);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inPGeoMeanCone(int32_t _3444_n);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inPGeoMeanCone();
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inRotatedQCone(std::shared_ptr< monty::ndarray< int32_t,1 > > _3445_dims);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inRotatedQCone(int32_t _3447_m,int32_t _3448_n);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inRotatedQCone(int32_t _3449_n);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inRotatedQCone();
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inQCone(std::shared_ptr< monty::ndarray< int32_t,1 > > _3450_dims);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inQCone(int32_t _3452_m,int32_t _3453_n);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inQCone(int32_t _3454_n);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inQCone();
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inSVecPSDCone(std::shared_ptr< monty::ndarray< int32_t,1 > > _3455_dims);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inSVecPSDCone(int32_t _3456_d1,int32_t _3457_d2);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inSVecPSDCone(int32_t _3458_n);
      static  monty::rc_ptr< ::mosek::fusion::ConeDomain > inSVecPSDCone();
      static  monty::rc_ptr< ::mosek::fusion::PSDDomain > isTrilPSD(int32_t _3459_n,int32_t _3460_m);
      static  monty::rc_ptr< ::mosek::fusion::PSDDomain > isTrilPSD(int32_t _3461_n);
      static  monty::rc_ptr< ::mosek::fusion::PSDDomain > isTrilPSD();
      static  monty::rc_ptr< ::mosek::fusion::PSDDomain > inPSDCone(int32_t _3462_n,int32_t _3463_m);
      static  monty::rc_ptr< ::mosek::fusion::PSDDomain > inPSDCone(int32_t _3464_n);
      static  monty::rc_ptr< ::mosek::fusion::PSDDomain > inPSDCone();
      static  monty::rc_ptr< ::mosek::fusion::RangeDomain > binary();
      static  monty::rc_ptr< ::mosek::fusion::RangeDomain > binary(std::shared_ptr< monty::ndarray< int32_t,1 > > _3465_dims);
      static  monty::rc_ptr< ::mosek::fusion::RangeDomain > binary(int32_t _3466_m,int32_t _3467_n);
      static  monty::rc_ptr< ::mosek::fusion::RangeDomain > binary(int32_t _3468_n);
      static  monty::rc_ptr< ::mosek::fusion::RangeDomain > inRange(monty::rc_ptr< ::mosek::fusion::Matrix > _3469_lbm,monty::rc_ptr< ::mosek::fusion::Matrix > _3470_ubm);
      static  monty::rc_ptr< ::mosek::fusion::RangeDomain > inRange(std::shared_ptr< monty::ndarray< double,2 > > _3471_lba,std::shared_ptr< monty::ndarray< double,2 > > _3472_uba);
      static  monty::rc_ptr< ::mosek::fusion::RangeDomain > inRange(std::shared_ptr< monty::ndarray< double,1 > > _3473_lba,std::shared_ptr< monty::ndarray< double,1 > > _3474_uba,std::shared_ptr< monty::ndarray< int32_t,1 > > _3475_dims);
      static  monty::rc_ptr< ::mosek::fusion::RangeDomain > inRange(std::shared_ptr< monty::ndarray< double,1 > > _3476_lba,double _3477_ub,std::shared_ptr< monty::ndarray< int32_t,1 > > _3478_dims);
      static  monty::rc_ptr< ::mosek::fusion::RangeDomain > inRange(double _3480_lb,std::shared_ptr< monty::ndarray< double,1 > > _3481_uba,std::shared_ptr< monty::ndarray< int32_t,1 > > _3482_dims);
      static  monty::rc_ptr< ::mosek::fusion::RangeDomain > inRange(double _3484_lb,double _3485_ub,std::shared_ptr< monty::ndarray< int32_t,1 > > _3486_dims);
      static  monty::rc_ptr< ::mosek::fusion::RangeDomain > inRange(std::shared_ptr< monty::ndarray< double,1 > > _3487_lba,std::shared_ptr< monty::ndarray< double,1 > > _3488_uba);
      static  monty::rc_ptr< ::mosek::fusion::RangeDomain > inRange(std::shared_ptr< monty::ndarray< double,1 > > _3489_lba,double _3490_ub);
      static  monty::rc_ptr< ::mosek::fusion::RangeDomain > inRange(double _3492_lb,std::shared_ptr< monty::ndarray< double,1 > > _3493_uba);
      static  monty::rc_ptr< ::mosek::fusion::RangeDomain > inRange(double _3495_lb,double _3496_ub);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > greaterThan(monty::rc_ptr< ::mosek::fusion::Matrix > _3497_mx);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > greaterThan(std::shared_ptr< monty::ndarray< double,1 > > _3498_a1,std::shared_ptr< monty::ndarray< int32_t,1 > > _3499_dims);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > greaterThan(std::shared_ptr< monty::ndarray< double,2 > > _3500_a2);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > greaterThan(std::shared_ptr< monty::ndarray< double,1 > > _3503_a1);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > greaterThan(double _3504_b,std::shared_ptr< monty::ndarray< int32_t,1 > > _3505_dims);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > greaterThan(double _3507_b,int32_t _3508_m,int32_t _3509_n);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > greaterThan(double _3511_b,int32_t _3512_n);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > greaterThan(double _3514_b);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > lessThan(monty::rc_ptr< ::mosek::fusion::Matrix > _3515_mx);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > lessThan(std::shared_ptr< monty::ndarray< double,1 > > _3516_a1,std::shared_ptr< monty::ndarray< int32_t,1 > > _3517_dims);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > lessThan(std::shared_ptr< monty::ndarray< double,2 > > _3518_a2);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > lessThan(std::shared_ptr< monty::ndarray< double,1 > > _3521_a1);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > lessThan(double _3522_b,std::shared_ptr< monty::ndarray< int32_t,1 > > _3523_dims);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > lessThan(double _3524_b,int32_t _3525_m,int32_t _3526_n);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > lessThan(double _3527_b,int32_t _3528_n);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > lessThan(double _3529_b);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > equalsTo(monty::rc_ptr< ::mosek::fusion::Matrix > _3530_mx);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > equalsTo(std::shared_ptr< monty::ndarray< double,1 > > _3531_a1,std::shared_ptr< monty::ndarray< int32_t,1 > > _3532_dims);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > equalsTo(std::shared_ptr< monty::ndarray< double,2 > > _3533_a2);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > equalsTo(std::shared_ptr< monty::ndarray< double,1 > > _3536_a1);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > equalsTo(double _3537_b,std::shared_ptr< monty::ndarray< int32_t,1 > > _3538_dims);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > equalsTo(double _3539_b,int32_t _3540_m,int32_t _3541_n);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > equalsTo(double _3542_b,int32_t _3543_n);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > equalsTo(double _3544_b);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > unbounded(std::shared_ptr< monty::ndarray< int32_t,1 > > _3545_dims);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > unbounded(int32_t _3547_m,int32_t _3548_n);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > unbounded(int32_t _3549_n);
      static  monty::rc_ptr< ::mosek::fusion::LinearDomain > unbounded();
    }; // struct Domain;

    struct p_Param
    {
      Param * _pubthis;
      static mosek::fusion::p_Param* _get_impl(mosek::fusion::Param * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_Param * _get_impl(mosek::fusion::Param::t _inst) { return _get_impl(_inst.get()); }
      p_Param(Param * _pubthis);
      virtual ~p_Param() { /* std::cout << "~p_Param" << std::endl;*/ };

      virtual void destroy();

      static  monty::rc_ptr< ::mosek::fusion::Parameter > repeat(monty::rc_ptr< ::mosek::fusion::Parameter > _3558_p,int32_t _3559_n,int32_t _3560_dim);
      static  monty::rc_ptr< ::mosek::fusion::Parameter > stack(int32_t _3562_dim,monty::rc_ptr< ::mosek::fusion::Parameter > _3563_p1,monty::rc_ptr< ::mosek::fusion::Parameter > _3564_p2,monty::rc_ptr< ::mosek::fusion::Parameter > _3565_p3);
      static  monty::rc_ptr< ::mosek::fusion::Parameter > stack(int32_t _3566_dim,monty::rc_ptr< ::mosek::fusion::Parameter > _3567_p1,monty::rc_ptr< ::mosek::fusion::Parameter > _3568_p2);
      static  monty::rc_ptr< ::mosek::fusion::Parameter > stack(int32_t _3569_dim,std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Parameter >,1 > > _3570_p);
      static  monty::rc_ptr< ::mosek::fusion::Parameter > stack(std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Parameter >,1 > >,1 > > _3571_p);
      static  monty::rc_ptr< ::mosek::fusion::Parameter > hstack(monty::rc_ptr< ::mosek::fusion::Parameter > _3573_p1,monty::rc_ptr< ::mosek::fusion::Parameter > _3574_p2,monty::rc_ptr< ::mosek::fusion::Parameter > _3575_p3);
      static  monty::rc_ptr< ::mosek::fusion::Parameter > hstack(monty::rc_ptr< ::mosek::fusion::Parameter > _3576_p1,monty::rc_ptr< ::mosek::fusion::Parameter > _3577_p2);
      static  monty::rc_ptr< ::mosek::fusion::Parameter > hstack(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Parameter >,1 > > _3578_p);
      static  monty::rc_ptr< ::mosek::fusion::Parameter > vstack(monty::rc_ptr< ::mosek::fusion::Parameter > _3579_p1,monty::rc_ptr< ::mosek::fusion::Parameter > _3580_p2,monty::rc_ptr< ::mosek::fusion::Parameter > _3581_p3);
      static  monty::rc_ptr< ::mosek::fusion::Parameter > vstack(monty::rc_ptr< ::mosek::fusion::Parameter > _3582_p1,monty::rc_ptr< ::mosek::fusion::Parameter > _3583_p2);
      static  monty::rc_ptr< ::mosek::fusion::Parameter > vstack(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Parameter >,1 > > _3584_p);
      static  monty::rc_ptr< ::mosek::fusion::Parameter > dstack(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Parameter >,1 > > _3585_p,int32_t _3586_dim);
    }; // struct Param;

    struct p_ParameterImpl : public /*implements*/ virtual ::mosek::fusion::Parameter
    {
      ParameterImpl * _pubthis;
      static mosek::fusion::p_ParameterImpl* _get_impl(mosek::fusion::ParameterImpl * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_ParameterImpl * _get_impl(mosek::fusion::ParameterImpl::t _inst) { return _get_impl(_inst.get()); }
      p_ParameterImpl(ParameterImpl * _pubthis);
      virtual ~p_ParameterImpl() { /* std::cout << "~p_ParameterImpl" << std::endl;*/ };
      int64_t size{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > nidxs{};
      std::shared_ptr< monty::ndarray< int64_t,1 > > sp{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};
      monty::rc_ptr< ::mosek::fusion::Model > model{};

      virtual void destroy();

      static ParameterImpl::t _new_ParameterImpl(monty::rc_ptr< ::mosek::fusion::ParameterImpl > _4375_other,monty::rc_ptr< ::mosek::fusion::Model > _4376_model);
      void _initialize(monty::rc_ptr< ::mosek::fusion::ParameterImpl > _4375_other,monty::rc_ptr< ::mosek::fusion::Model > _4376_model);
      static ParameterImpl::t _new_ParameterImpl(monty::rc_ptr< ::mosek::fusion::Model > _4377_model,std::shared_ptr< monty::ndarray< int32_t,1 > > _4378_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _4379_sp,std::shared_ptr< monty::ndarray< int32_t,1 > > _4380_nidxs);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Model > _4377_model,std::shared_ptr< monty::ndarray< int32_t,1 > > _4378_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _4379_sp,std::shared_ptr< monty::ndarray< int32_t,1 > > _4380_nidxs);
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2ParameterImpl__clone(monty::rc_ptr< ::mosek::fusion::Model > _4381_m) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Parameter__clone(monty::rc_ptr< ::mosek::fusion::Model > _4381_m) { return __mosek_2fusion_2ParameterImpl__clone(_4381_m); }
      virtual /* override */ std::string toString() ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2ParameterImpl__pick(std::shared_ptr< monty::ndarray< int32_t,2 > > _4384_midxs) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Parameter__pick(std::shared_ptr< monty::ndarray< int32_t,2 > > _4384_midxs) { return __mosek_2fusion_2ParameterImpl__pick(_4384_midxs); }
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2ParameterImpl__pick(std::shared_ptr< monty::ndarray< int32_t,1 > > _4406_idxs) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Parameter__pick(std::shared_ptr< monty::ndarray< int32_t,1 > > _4406_idxs) { return __mosek_2fusion_2ParameterImpl__pick(_4406_idxs); }
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2ParameterImpl__index(std::shared_ptr< monty::ndarray< int32_t,1 > > _4417_indexes) ;
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2Expression__index(std::shared_ptr< monty::ndarray< int32_t,1 > > _4417_indexes) { return __mosek_2fusion_2ParameterImpl__index(_4417_indexes); }
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2ParameterImpl__index(int32_t _4426_i) ;
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2Expression__index(int32_t _4426_i) { return __mosek_2fusion_2ParameterImpl__index(_4426_i); }
      virtual void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _4428_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _4429_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _4430_xs) ;
      virtual void getSp(std::shared_ptr< monty::ndarray< int64_t,1 > > _4452_dest,int32_t _4453_offset) ;
      virtual bool isSparse() ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2ParameterImpl__slice(std::shared_ptr< monty::ndarray< int32_t,1 > > _4456_astart,std::shared_ptr< monty::ndarray< int32_t,1 > > _4457_astop) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Parameter__slice(std::shared_ptr< monty::ndarray< int32_t,1 > > _4456_astart,std::shared_ptr< monty::ndarray< int32_t,1 > > _4457_astop) { return __mosek_2fusion_2ParameterImpl__slice(_4456_astart,_4457_astop); }
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2ParameterImpl__transpose() ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Parameter__transpose() { return __mosek_2fusion_2ParameterImpl__transpose(); }
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2ParameterImpl__slice(int32_t _4503_start,int32_t _4504_stop) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Parameter__slice(int32_t _4503_start,int32_t _4504_stop) { return __mosek_2fusion_2ParameterImpl__slice(_4503_start,_4504_stop); }
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2ParameterImpl__reshape(std::shared_ptr< monty::ndarray< int32_t,1 > > _4512_dims) ;
      virtual monty::rc_ptr< ::mosek::fusion::Parameter > __mosek_2fusion_2Parameter__reshape(std::shared_ptr< monty::ndarray< int32_t,1 > > _4512_dims) { return __mosek_2fusion_2ParameterImpl__reshape(_4512_dims); }
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2ParameterImpl__asExpr() ;
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2Parameter__asExpr() { return __mosek_2fusion_2ParameterImpl__asExpr(); }
      virtual int64_t getSize() ;
      virtual int32_t getNumNonzero() ;
      virtual int32_t getND() ;
      virtual std::shared_ptr< monty::ndarray< int32_t,1 > > getShape() ;
      virtual int32_t getDim(int32_t _4513_i) ;
      virtual void getAllIndexes(std::shared_ptr< monty::ndarray< int32_t,1 > > _4514_dst,int32_t _4515_ofs) ;
      virtual int32_t getIndex(int32_t _4517_i) ;
      virtual std::shared_ptr< monty::ndarray< double,1 > > getValue() ;
      virtual void setValue(std::shared_ptr< monty::ndarray< double,2 > > _4518_values2) ;
      virtual void setValue(std::shared_ptr< monty::ndarray< double,1 > > _4524_values) ;
      virtual void setValue(double _4527_value) ;
      virtual monty::rc_ptr< ::mosek::fusion::Model > __mosek_2fusion_2ParameterImpl__getModel() ;
      virtual monty::rc_ptr< ::mosek::fusion::Model > __mosek_2fusion_2Expression__getModel() { return __mosek_2fusion_2ParameterImpl__getModel(); }
    }; // struct ParameterImpl;

    struct p_ExprRangeDomain : public /*implements*/ virtual ::mosek::fusion::ExprDomain
    {
      ExprRangeDomain * _pubthis;
      static mosek::fusion::p_ExprRangeDomain* _get_impl(mosek::fusion::ExprRangeDomain * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_ExprRangeDomain * _get_impl(mosek::fusion::ExprRangeDomain::t _inst) { return _get_impl(_inst.get()); }
      p_ExprRangeDomain(ExprRangeDomain * _pubthis);
      virtual ~p_ExprRangeDomain() { /* std::cout << "~p_ExprRangeDomain" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::RangeDomain > dom{};
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};

      virtual void destroy();

      static ExprRangeDomain::t _new_ExprRangeDomain(monty::rc_ptr< ::mosek::fusion::Expression > _7199_expr,monty::rc_ptr< ::mosek::fusion::RangeDomain > _7200_dom);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _7199_expr,monty::rc_ptr< ::mosek::fusion::RangeDomain > _7200_dom);
      virtual monty::rc_ptr< ::mosek::fusion::Term > __mosek_2fusion_2ExprRangeDomain__toDJCTerm() ;
      virtual monty::rc_ptr< ::mosek::fusion::Term > __mosek_2fusion_2ExprDomain__toDJCTerm() { return __mosek_2fusion_2ExprRangeDomain__toDJCTerm(); }
    }; // struct ExprRangeDomain;

    struct p_ExprPSDDomain
    {
      ExprPSDDomain * _pubthis;
      static mosek::fusion::p_ExprPSDDomain* _get_impl(mosek::fusion::ExprPSDDomain * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_ExprPSDDomain * _get_impl(mosek::fusion::ExprPSDDomain::t _inst) { return _get_impl(_inst.get()); }
      p_ExprPSDDomain(ExprPSDDomain * _pubthis);
      virtual ~p_ExprPSDDomain() { /* std::cout << "~p_ExprPSDDomain" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::PSDDomain > dom{};
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};

      virtual void destroy();

      static ExprPSDDomain::t _new_ExprPSDDomain(monty::rc_ptr< ::mosek::fusion::Expression > _7201_expr,monty::rc_ptr< ::mosek::fusion::PSDDomain > _7202_dom);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _7201_expr,monty::rc_ptr< ::mosek::fusion::PSDDomain > _7202_dom);
    }; // struct ExprPSDDomain;

    struct p_ExprConicDomain
    {
      ExprConicDomain * _pubthis;
      static mosek::fusion::p_ExprConicDomain* _get_impl(mosek::fusion::ExprConicDomain * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_ExprConicDomain * _get_impl(mosek::fusion::ExprConicDomain::t _inst) { return _get_impl(_inst.get()); }
      p_ExprConicDomain(ExprConicDomain * _pubthis);
      virtual ~p_ExprConicDomain() { /* std::cout << "~p_ExprConicDomain" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::ConeDomain > dom{};
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};

      virtual void destroy();

      static ExprConicDomain::t _new_ExprConicDomain(monty::rc_ptr< ::mosek::fusion::Expression > _7203_expr,monty::rc_ptr< ::mosek::fusion::ConeDomain > _7204_dom);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _7203_expr,monty::rc_ptr< ::mosek::fusion::ConeDomain > _7204_dom);
    }; // struct ExprConicDomain;

    struct p_ExprLinearDomain : public /*implements*/ virtual ::mosek::fusion::ExprDomain
    {
      ExprLinearDomain * _pubthis;
      static mosek::fusion::p_ExprLinearDomain* _get_impl(mosek::fusion::ExprLinearDomain * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_ExprLinearDomain * _get_impl(mosek::fusion::ExprLinearDomain::t _inst) { return _get_impl(_inst.get()); }
      p_ExprLinearDomain(ExprLinearDomain * _pubthis);
      virtual ~p_ExprLinearDomain() { /* std::cout << "~p_ExprLinearDomain" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::LinearDomain > dom{};
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};

      virtual void destroy();

      static ExprLinearDomain::t _new_ExprLinearDomain(monty::rc_ptr< ::mosek::fusion::Expression > _7205_expr,monty::rc_ptr< ::mosek::fusion::LinearDomain > _7206_dom);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _7205_expr,monty::rc_ptr< ::mosek::fusion::LinearDomain > _7206_dom);
      virtual monty::rc_ptr< ::mosek::fusion::Term > __mosek_2fusion_2ExprLinearDomain__toDJCTerm() ;
      virtual monty::rc_ptr< ::mosek::fusion::Term > __mosek_2fusion_2ExprDomain__toDJCTerm() { return __mosek_2fusion_2ExprLinearDomain__toDJCTerm(); }
    }; // struct ExprLinearDomain;

    struct p_BaseExpression : public /*implements*/ virtual ::mosek::fusion::Expression
    {
      BaseExpression * _pubthis;
      static mosek::fusion::p_BaseExpression* _get_impl(mosek::fusion::BaseExpression * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_BaseExpression * _get_impl(mosek::fusion::BaseExpression::t _inst) { return _get_impl(_inst.get()); }
      p_BaseExpression(BaseExpression * _pubthis);
      virtual ~p_BaseExpression() { /* std::cout << "~p_BaseExpression" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Model > owner{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};

      virtual void destroy();

      static BaseExpression::t _new_BaseExpression(monty::rc_ptr< ::mosek::fusion::Model > _7207_owner,std::shared_ptr< monty::ndarray< int32_t,1 > > _7208_shape);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Model > _7207_owner,std::shared_ptr< monty::ndarray< int32_t,1 > > _7208_shape);
      virtual /* override */ std::string toString() ;
      virtual void printStack(monty::rc_ptr< ::mosek::fusion::WorkStack > _7209_rs) ;
      virtual void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _7226_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _7227_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _7228_xs) { throw monty::AbstractClassError("Call to abstract method"); }
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2BaseExpression__pick(std::shared_ptr< monty::ndarray< int32_t,2 > > _7229_indexrows) ;
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2Expression__pick(std::shared_ptr< monty::ndarray< int32_t,2 > > _7229_indexrows) { return __mosek_2fusion_2BaseExpression__pick(_7229_indexrows); }
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2BaseExpression__pick(std::shared_ptr< monty::ndarray< int32_t,1 > > _7230_indexes) ;
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2Expression__pick(std::shared_ptr< monty::ndarray< int32_t,1 > > _7230_indexes) { return __mosek_2fusion_2BaseExpression__pick(_7230_indexes); }
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2BaseExpression__index(std::shared_ptr< monty::ndarray< int32_t,1 > > _7233_indexes) ;
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2Expression__index(std::shared_ptr< monty::ndarray< int32_t,1 > > _7233_indexes) { return __mosek_2fusion_2BaseExpression__index(_7233_indexes); }
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2BaseExpression__index(int32_t _7236_i) ;
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2Expression__index(int32_t _7236_i) { return __mosek_2fusion_2BaseExpression__index(_7236_i); }
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2BaseExpression__slice(std::shared_ptr< monty::ndarray< int32_t,1 > > _7238_firsta,std::shared_ptr< monty::ndarray< int32_t,1 > > _7239_lasta) ;
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2Expression__slice(std::shared_ptr< monty::ndarray< int32_t,1 > > _7238_firsta,std::shared_ptr< monty::ndarray< int32_t,1 > > _7239_lasta) { return __mosek_2fusion_2BaseExpression__slice(_7238_firsta,_7239_lasta); }
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2BaseExpression__slice(int32_t _7240_first,int32_t _7241_last) ;
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2Expression__slice(int32_t _7240_first,int32_t _7241_last) { return __mosek_2fusion_2BaseExpression__slice(_7240_first,_7241_last); }
      virtual int64_t getSize() ;
      virtual int32_t getND() ;
      virtual int32_t getDim(int32_t _7242_d) ;
      virtual std::shared_ptr< monty::ndarray< int32_t,1 > > getShape() ;
      virtual monty::rc_ptr< ::mosek::fusion::Model > __mosek_2fusion_2BaseExpression__getModel() ;
      virtual monty::rc_ptr< ::mosek::fusion::Model > __mosek_2fusion_2Expression__getModel() { return __mosek_2fusion_2BaseExpression__getModel(); }
      static  monty::rc_ptr< ::mosek::fusion::Model > modelFrom(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Expression >,1 > > _7243_es);
      static  monty::rc_ptr< ::mosek::fusion::Model > modelFrom(monty::rc_ptr< ::mosek::fusion::Model > _7246_m1,monty::rc_ptr< ::mosek::fusion::Expression > _7247_e2);
      static  monty::rc_ptr< ::mosek::fusion::Model > modelFrom(monty::rc_ptr< ::mosek::fusion::Expression > _7248_e1,monty::rc_ptr< ::mosek::fusion::Expression > _7249_e2);
      static  monty::rc_ptr< ::mosek::fusion::Model > modelFrom(monty::rc_ptr< ::mosek::fusion::Model > _7250_m1,monty::rc_ptr< ::mosek::fusion::Model > _7251_m2);
    }; // struct BaseExpression;

    struct p_ExprParameter : public ::mosek::fusion::p_BaseExpression
    {
      ExprParameter * _pubthis;
      static mosek::fusion::p_ExprParameter* _get_impl(mosek::fusion::ExprParameter * _inst){ return static_cast< mosek::fusion::p_ExprParameter* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprParameter * _get_impl(mosek::fusion::ExprParameter::t _inst) { return _get_impl(_inst.get()); }
      p_ExprParameter(ExprParameter * _pubthis);
      virtual ~p_ExprParameter() { /* std::cout << "~p_ExprParameter" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Parameter > p{};

      virtual void destroy();

      static ExprParameter::t _new_ExprParameter(monty::rc_ptr< ::mosek::fusion::Parameter > _3550_p);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Parameter > _3550_p);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _3551_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _3552_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _3553_xs) ;
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2ExprParameter__slice(std::shared_ptr< monty::ndarray< int32_t,1 > > _3554_start,std::shared_ptr< monty::ndarray< int32_t,1 > > _3555_stop) ;
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2BaseExpression__slice(std::shared_ptr< monty::ndarray< int32_t,1 > > _3554_start,std::shared_ptr< monty::ndarray< int32_t,1 > > _3555_stop) { return __mosek_2fusion_2ExprParameter__slice(_3554_start,_3555_stop); }
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2ExprParameter__slice(int32_t _3556_start,int32_t _3557_stop) ;
      virtual monty::rc_ptr< ::mosek::fusion::Expression > __mosek_2fusion_2BaseExpression__slice(int32_t _3556_start,int32_t _3557_stop) { return __mosek_2fusion_2ExprParameter__slice(_3556_start,_3557_stop); }
      virtual /* override */ std::string toString() ;
    }; // struct ExprParameter;

    struct p_ExprMulParamScalarExpr : public ::mosek::fusion::p_BaseExpression
    {
      ExprMulParamScalarExpr * _pubthis;
      static mosek::fusion::p_ExprMulParamScalarExpr* _get_impl(mosek::fusion::ExprMulParamScalarExpr * _inst){ return static_cast< mosek::fusion::p_ExprMulParamScalarExpr* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprMulParamScalarExpr * _get_impl(mosek::fusion::ExprMulParamScalarExpr::t _inst) { return _get_impl(_inst.get()); }
      p_ExprMulParamScalarExpr(ExprMulParamScalarExpr * _pubthis);
      virtual ~p_ExprMulParamScalarExpr() { /* std::cout << "~p_ExprMulParamScalarExpr" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Expression > e{};
      monty::rc_ptr< ::mosek::fusion::Parameter > p{};

      virtual void destroy();

      static ExprMulParamScalarExpr::t _new_ExprMulParamScalarExpr(monty::rc_ptr< ::mosek::fusion::Parameter > _3645_p,monty::rc_ptr< ::mosek::fusion::Expression > _3646_e);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Parameter > _3645_p,monty::rc_ptr< ::mosek::fusion::Expression > _3646_e);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _3647_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _3648_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _3649_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprMulParamScalarExpr;

    struct p_ExprMulParamScalar : public ::mosek::fusion::p_BaseExpression
    {
      ExprMulParamScalar * _pubthis;
      static mosek::fusion::p_ExprMulParamScalar* _get_impl(mosek::fusion::ExprMulParamScalar * _inst){ return static_cast< mosek::fusion::p_ExprMulParamScalar* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprMulParamScalar * _get_impl(mosek::fusion::ExprMulParamScalar::t _inst) { return _get_impl(_inst.get()); }
      p_ExprMulParamScalar(ExprMulParamScalar * _pubthis);
      virtual ~p_ExprMulParamScalar() { /* std::cout << "~p_ExprMulParamScalar" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Expression > e{};
      monty::rc_ptr< ::mosek::fusion::Parameter > p{};

      virtual void destroy();

      static ExprMulParamScalar::t _new_ExprMulParamScalar(monty::rc_ptr< ::mosek::fusion::Parameter > _3700_p,monty::rc_ptr< ::mosek::fusion::Expression > _3701_e);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Parameter > _3700_p,monty::rc_ptr< ::mosek::fusion::Expression > _3701_e);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _3702_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _3703_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _3704_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprMulParamScalar;

    struct p_ExprMulParamDiagLeft : public ::mosek::fusion::p_BaseExpression
    {
      ExprMulParamDiagLeft * _pubthis;
      static mosek::fusion::p_ExprMulParamDiagLeft* _get_impl(mosek::fusion::ExprMulParamDiagLeft * _inst){ return static_cast< mosek::fusion::p_ExprMulParamDiagLeft* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprMulParamDiagLeft * _get_impl(mosek::fusion::ExprMulParamDiagLeft::t _inst) { return _get_impl(_inst.get()); }
      p_ExprMulParamDiagLeft(ExprMulParamDiagLeft * _pubthis);
      virtual ~p_ExprMulParamDiagLeft() { /* std::cout << "~p_ExprMulParamDiagLeft" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Expression > e{};
      monty::rc_ptr< ::mosek::fusion::Parameter > p{};

      virtual void destroy();

      static ExprMulParamDiagLeft::t _new_ExprMulParamDiagLeft(monty::rc_ptr< ::mosek::fusion::Parameter > _3747_p,monty::rc_ptr< ::mosek::fusion::Expression > _3748_e);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Parameter > _3747_p,monty::rc_ptr< ::mosek::fusion::Expression > _3748_e);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _3749_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _3750_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _3751_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprMulParamDiagLeft;

    struct p_ExprMulParamDiagRight : public ::mosek::fusion::p_BaseExpression
    {
      ExprMulParamDiagRight * _pubthis;
      static mosek::fusion::p_ExprMulParamDiagRight* _get_impl(mosek::fusion::ExprMulParamDiagRight * _inst){ return static_cast< mosek::fusion::p_ExprMulParamDiagRight* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprMulParamDiagRight * _get_impl(mosek::fusion::ExprMulParamDiagRight::t _inst) { return _get_impl(_inst.get()); }
      p_ExprMulParamDiagRight(ExprMulParamDiagRight * _pubthis);
      virtual ~p_ExprMulParamDiagRight() { /* std::cout << "~p_ExprMulParamDiagRight" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Expression > e{};
      monty::rc_ptr< ::mosek::fusion::Parameter > p{};

      virtual void destroy();

      static ExprMulParamDiagRight::t _new_ExprMulParamDiagRight(monty::rc_ptr< ::mosek::fusion::Expression > _3866_e,monty::rc_ptr< ::mosek::fusion::Parameter > _3867_p);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _3866_e,monty::rc_ptr< ::mosek::fusion::Parameter > _3867_p);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _3868_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _3869_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _3870_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprMulParamDiagRight;

    struct p_ExprDotParam : public ::mosek::fusion::p_BaseExpression
    {
      ExprDotParam * _pubthis;
      static mosek::fusion::p_ExprDotParam* _get_impl(mosek::fusion::ExprDotParam * _inst){ return static_cast< mosek::fusion::p_ExprDotParam* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprDotParam * _get_impl(mosek::fusion::ExprDotParam::t _inst) { return _get_impl(_inst.get()); }
      p_ExprDotParam(ExprDotParam * _pubthis);
      virtual ~p_ExprDotParam() { /* std::cout << "~p_ExprDotParam" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Expression > e{};
      monty::rc_ptr< ::mosek::fusion::Parameter > p{};

      virtual void destroy();

      static ExprDotParam::t _new_ExprDotParam(monty::rc_ptr< ::mosek::fusion::Parameter > _3984_p,monty::rc_ptr< ::mosek::fusion::Expression > _3985_e);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Parameter > _3984_p,monty::rc_ptr< ::mosek::fusion::Expression > _3985_e);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _3987_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _3988_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _3989_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprDotParam;

    struct p_ExprMulParamElem : public ::mosek::fusion::p_BaseExpression
    {
      ExprMulParamElem * _pubthis;
      static mosek::fusion::p_ExprMulParamElem* _get_impl(mosek::fusion::ExprMulParamElem * _inst){ return static_cast< mosek::fusion::p_ExprMulParamElem* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprMulParamElem * _get_impl(mosek::fusion::ExprMulParamElem::t _inst) { return _get_impl(_inst.get()); }
      p_ExprMulParamElem(ExprMulParamElem * _pubthis);
      virtual ~p_ExprMulParamElem() { /* std::cout << "~p_ExprMulParamElem" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Expression > e{};
      monty::rc_ptr< ::mosek::fusion::Parameter > p{};

      virtual void destroy();

      static ExprMulParamElem::t _new_ExprMulParamElem(monty::rc_ptr< ::mosek::fusion::Parameter > _4047_p,monty::rc_ptr< ::mosek::fusion::Expression > _4048_e);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Parameter > _4047_p,monty::rc_ptr< ::mosek::fusion::Expression > _4048_e);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _4050_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _4051_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _4052_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprMulParamElem;

    struct p_ExprMulParamRight : public ::mosek::fusion::p_BaseExpression
    {
      ExprMulParamRight * _pubthis;
      static mosek::fusion::p_ExprMulParamRight* _get_impl(mosek::fusion::ExprMulParamRight * _inst){ return static_cast< mosek::fusion::p_ExprMulParamRight* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprMulParamRight * _get_impl(mosek::fusion::ExprMulParamRight::t _inst) { return _get_impl(_inst.get()); }
      p_ExprMulParamRight(ExprMulParamRight * _pubthis);
      virtual ~p_ExprMulParamRight() { /* std::cout << "~p_ExprMulParamRight" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Expression > e{};
      monty::rc_ptr< ::mosek::fusion::Parameter > p{};

      virtual void destroy();

      static ExprMulParamRight::t _new_ExprMulParamRight(monty::rc_ptr< ::mosek::fusion::Expression > _4114_e,monty::rc_ptr< ::mosek::fusion::Parameter > _4115_p);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _4114_e,monty::rc_ptr< ::mosek::fusion::Parameter > _4115_p);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _4116_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _4117_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _4118_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprMulParamRight;

    struct p_ExprMulParamLeft : public ::mosek::fusion::p_BaseExpression
    {
      ExprMulParamLeft * _pubthis;
      static mosek::fusion::p_ExprMulParamLeft* _get_impl(mosek::fusion::ExprMulParamLeft * _inst){ return static_cast< mosek::fusion::p_ExprMulParamLeft* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprMulParamLeft * _get_impl(mosek::fusion::ExprMulParamLeft::t _inst) { return _get_impl(_inst.get()); }
      p_ExprMulParamLeft(ExprMulParamLeft * _pubthis);
      virtual ~p_ExprMulParamLeft() { /* std::cout << "~p_ExprMulParamLeft" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Expression > e{};
      monty::rc_ptr< ::mosek::fusion::Parameter > p{};

      virtual void destroy();

      static ExprMulParamLeft::t _new_ExprMulParamLeft(monty::rc_ptr< ::mosek::fusion::Parameter > _4218_p,monty::rc_ptr< ::mosek::fusion::Expression > _4219_e);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Parameter > _4218_p,monty::rc_ptr< ::mosek::fusion::Expression > _4219_e);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _4220_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _4221_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _4222_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprMulParamLeft;

    struct p_ExprOptimizeCode : public ::mosek::fusion::p_BaseExpression
    {
      ExprOptimizeCode * _pubthis;
      static mosek::fusion::p_ExprOptimizeCode* _get_impl(mosek::fusion::ExprOptimizeCode * _inst){ return static_cast< mosek::fusion::p_ExprOptimizeCode* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprOptimizeCode * _get_impl(mosek::fusion::ExprOptimizeCode::t _inst) { return _get_impl(_inst.get()); }
      p_ExprOptimizeCode(ExprOptimizeCode * _pubthis);
      virtual ~p_ExprOptimizeCode() { /* std::cout << "~p_ExprOptimizeCode" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};

      virtual void destroy();

      static ExprOptimizeCode::t _new_ExprOptimizeCode(monty::rc_ptr< ::mosek::fusion::Expression > _4547_expr);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _4547_expr);
      static  void compress_code(monty::rc_ptr< ::mosek::fusion::WorkStack > _4548_xs,int32_t _4549_n,std::shared_ptr< monty::ndarray< int32_t,1 > > _4550_code,int32_t _4551_code_base,std::shared_ptr< monty::ndarray< int32_t,1 > > _4552_ptr,int32_t _4553_ptr_base,std::shared_ptr< monty::ndarray< double,1 > > _4554_fixterm,int32_t _4555_fixterm_base,std::shared_ptr< monty::ndarray< double,1 > > _4556_code_consts,int32_t _4557_code_consts_base,int32_t _4558_target_code_base,int32_t _4559_target_const_base,int32_t _4560_target_ptr_base);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _4613_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _4614_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _4615_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprOptimizeCode;

    struct p_ExprCompress : public ::mosek::fusion::p_BaseExpression
    {
      ExprCompress * _pubthis;
      static mosek::fusion::p_ExprCompress* _get_impl(mosek::fusion::ExprCompress * _inst){ return static_cast< mosek::fusion::p_ExprCompress* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprCompress * _get_impl(mosek::fusion::ExprCompress::t _inst) { return _get_impl(_inst.get()); }
      p_ExprCompress(ExprCompress * _pubthis);
      virtual ~p_ExprCompress() { /* std::cout << "~p_ExprCompress" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};

      virtual void destroy();

      static ExprCompress::t _new_ExprCompress(monty::rc_ptr< ::mosek::fusion::Expression > _4682_expr);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _4682_expr);
      static  void arg_sort(monty::rc_ptr< ::mosek::fusion::WorkStack > _4683_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _4684_xs,int32_t _4685_perm,int32_t _4686_nelem,int32_t _4687_nnz,int32_t _4688_ptr,int32_t _4689_nidxs);
      static  void merge_sort(int32_t _4725_origperm1,int32_t _4726_origperm2,int32_t _4727_nelem,int32_t _4728_nnz,int32_t _4729_ptr_base,int32_t _4730_nidxs_base,std::shared_ptr< monty::ndarray< int32_t,1 > > _4731_wi32,std::shared_ptr< monty::ndarray< int64_t,1 > > _4732_wi64);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _4755_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _4756_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _4757_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprCompress;

    struct p_ExprConst : public ::mosek::fusion::p_BaseExpression
    {
      ExprConst * _pubthis;
      static mosek::fusion::p_ExprConst* _get_impl(mosek::fusion::ExprConst * _inst){ return static_cast< mosek::fusion::p_ExprConst* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprConst * _get_impl(mosek::fusion::ExprConst::t _inst) { return _get_impl(_inst.get()); }
      p_ExprConst(ExprConst * _pubthis);
      virtual ~p_ExprConst() { /* std::cout << "~p_ExprConst" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< int64_t,1 > > sparsity{};
      std::shared_ptr< monty::ndarray< double,1 > > bfix{};

      virtual void destroy();

      static ExprConst::t _new_ExprConst(std::shared_ptr< monty::ndarray< int32_t,1 > > _4843_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _4844_sparsity,std::shared_ptr< monty::ndarray< double,1 > > _4845_bfix);
      void _initialize(std::shared_ptr< monty::ndarray< int32_t,1 > > _4843_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _4844_sparsity,std::shared_ptr< monty::ndarray< double,1 > > _4845_bfix);
      static ExprConst::t _new_ExprConst(std::shared_ptr< monty::ndarray< int32_t,1 > > _4846_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _4847_sparsity,double _4848_bfix);
      void _initialize(std::shared_ptr< monty::ndarray< int32_t,1 > > _4846_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _4847_sparsity,double _4848_bfix);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _4851_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _4852_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _4853_xs) ;
      static  void validate(std::shared_ptr< monty::ndarray< int32_t,1 > > _4872_shape,std::shared_ptr< monty::ndarray< double,1 > > _4873_bfix,std::shared_ptr< monty::ndarray< int64_t,1 > > _4874_sparsity);
      virtual /* override */ std::string toString() ;
    }; // struct ExprConst;

    struct p_ExprPick : public ::mosek::fusion::p_BaseExpression
    {
      ExprPick * _pubthis;
      static mosek::fusion::p_ExprPick* _get_impl(mosek::fusion::ExprPick * _inst){ return static_cast< mosek::fusion::p_ExprPick* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprPick * _get_impl(mosek::fusion::ExprPick::t _inst) { return _get_impl(_inst.get()); }
      p_ExprPick(ExprPick * _pubthis);
      virtual ~p_ExprPick() { /* std::cout << "~p_ExprPick" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< int64_t,1 > > idxs{};
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};

      virtual void destroy();

      static ExprPick::t _new_ExprPick(monty::rc_ptr< ::mosek::fusion::Expression > _4878_expr,std::shared_ptr< monty::ndarray< int32_t,2 > > _4879_idxs);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _4878_expr,std::shared_ptr< monty::ndarray< int32_t,2 > > _4879_idxs);
      static ExprPick::t _new_ExprPick(monty::rc_ptr< ::mosek::fusion::Expression > _4891_expr,std::shared_ptr< monty::ndarray< int64_t,1 > > _4892_idxs);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _4891_expr,std::shared_ptr< monty::ndarray< int64_t,1 > > _4892_idxs);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _4897_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _4898_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _4899_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprPick;

    struct p_ExprSlice : public ::mosek::fusion::p_BaseExpression
    {
      ExprSlice * _pubthis;
      static mosek::fusion::p_ExprSlice* _get_impl(mosek::fusion::ExprSlice * _inst){ return static_cast< mosek::fusion::p_ExprSlice* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprSlice * _get_impl(mosek::fusion::ExprSlice::t _inst) { return _get_impl(_inst.get()); }
      p_ExprSlice(ExprSlice * _pubthis);
      virtual ~p_ExprSlice() { /* std::cout << "~p_ExprSlice" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< int32_t,1 > > last{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > first{};
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};

      virtual void destroy();

      static ExprSlice::t _new_ExprSlice(monty::rc_ptr< ::mosek::fusion::Expression > _4964_expr,std::shared_ptr< monty::ndarray< int32_t,1 > > _4965_first,std::shared_ptr< monty::ndarray< int32_t,1 > > _4966_last);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _4964_expr,std::shared_ptr< monty::ndarray< int32_t,1 > > _4965_first,std::shared_ptr< monty::ndarray< int32_t,1 > > _4966_last);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _4967_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _4968_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _4969_xs) ;
      static  std::shared_ptr< monty::ndarray< int32_t,1 > > makeShape(std::shared_ptr< monty::ndarray< int32_t,1 > > _5034_shape,std::shared_ptr< monty::ndarray< int32_t,1 > > _5035_first,std::shared_ptr< monty::ndarray< int32_t,1 > > _5036_last);
      virtual /* override */ std::string toString() ;
    }; // struct ExprSlice;

    struct p_ExprPermuteDims : public ::mosek::fusion::p_BaseExpression
    {
      ExprPermuteDims * _pubthis;
      static mosek::fusion::p_ExprPermuteDims* _get_impl(mosek::fusion::ExprPermuteDims * _inst){ return static_cast< mosek::fusion::p_ExprPermuteDims* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprPermuteDims * _get_impl(mosek::fusion::ExprPermuteDims::t _inst) { return _get_impl(_inst.get()); }
      p_ExprPermuteDims(ExprPermuteDims * _pubthis);
      virtual ~p_ExprPermuteDims() { /* std::cout << "~p_ExprPermuteDims" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< int32_t,1 > > dperm{};
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};

      virtual void destroy();

      static ExprPermuteDims::t _new_ExprPermuteDims(std::shared_ptr< monty::ndarray< int32_t,1 > > _5041_perm,monty::rc_ptr< ::mosek::fusion::Expression > _5042_expr);
      void _initialize(std::shared_ptr< monty::ndarray< int32_t,1 > > _5041_perm,monty::rc_ptr< ::mosek::fusion::Expression > _5042_expr);
      static ExprPermuteDims::t _new_ExprPermuteDims(std::shared_ptr< monty::ndarray< int32_t,1 > > _5048_perm,monty::rc_ptr< ::mosek::fusion::Expression > _5049_expr,int32_t _5050_validated);
      void _initialize(std::shared_ptr< monty::ndarray< int32_t,1 > > _5048_perm,monty::rc_ptr< ::mosek::fusion::Expression > _5049_expr,int32_t _5050_validated);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _5051_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _5052_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _5053_xs) ;
      static  std::shared_ptr< monty::ndarray< int32_t,1 > > computeshape(std::shared_ptr< monty::ndarray< int32_t,1 > > _5111_perm,std::shared_ptr< monty::ndarray< int32_t,1 > > _5112_shape);
    }; // struct ExprPermuteDims;

    struct p_ExprTranspose : public ::mosek::fusion::p_BaseExpression
    {
      ExprTranspose * _pubthis;
      static mosek::fusion::p_ExprTranspose* _get_impl(mosek::fusion::ExprTranspose * _inst){ return static_cast< mosek::fusion::p_ExprTranspose* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprTranspose * _get_impl(mosek::fusion::ExprTranspose::t _inst) { return _get_impl(_inst.get()); }
      p_ExprTranspose(ExprTranspose * _pubthis);
      virtual ~p_ExprTranspose() { /* std::cout << "~p_ExprTranspose" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};

      virtual void destroy();

      static ExprTranspose::t _new_ExprTranspose(monty::rc_ptr< ::mosek::fusion::Expression > _5114_expr);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _5114_expr);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _5115_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _5116_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _5117_xs) ;
      virtual /* override */ std::string toString() ;
      static  std::shared_ptr< monty::ndarray< int32_t,1 > > transposeShape(std::shared_ptr< monty::ndarray< int32_t,1 > > _5172_shape);
    }; // struct ExprTranspose;

    struct p_ExprRepeat : public ::mosek::fusion::p_BaseExpression
    {
      ExprRepeat * _pubthis;
      static mosek::fusion::p_ExprRepeat* _get_impl(mosek::fusion::ExprRepeat * _inst){ return static_cast< mosek::fusion::p_ExprRepeat* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprRepeat * _get_impl(mosek::fusion::ExprRepeat::t _inst) { return _get_impl(_inst.get()); }
      p_ExprRepeat(ExprRepeat * _pubthis);
      virtual ~p_ExprRepeat() { /* std::cout << "~p_ExprRepeat" << std::endl;*/ };
      int32_t n{};
      int32_t dim{};
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};

      virtual void destroy();

      static ExprRepeat::t _new_ExprRepeat(monty::rc_ptr< ::mosek::fusion::Expression > _5173_expr,int32_t _5174_dim,int32_t _5175_n);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _5173_expr,int32_t _5174_dim,int32_t _5175_n);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _5176_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _5177_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _5178_xs) ;
      static  std::shared_ptr< monty::ndarray< int32_t,1 > > getshape(monty::rc_ptr< ::mosek::fusion::Expression > _5243_e,int32_t _5244_dim,int32_t _5245_n);
      virtual /* override */ std::string toString() ;
    }; // struct ExprRepeat;

    struct p_ExprStack : public ::mosek::fusion::p_BaseExpression
    {
      ExprStack * _pubthis;
      static mosek::fusion::p_ExprStack* _get_impl(mosek::fusion::ExprStack * _inst){ return static_cast< mosek::fusion::p_ExprStack* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprStack * _get_impl(mosek::fusion::ExprStack::t _inst) { return _get_impl(_inst.get()); }
      p_ExprStack(ExprStack * _pubthis);
      virtual ~p_ExprStack() { /* std::cout << "~p_ExprStack" << std::endl;*/ };
      int32_t dim{};
      std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Expression >,1 > > exprs{};

      virtual void destroy();

      static ExprStack::t _new_ExprStack(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Expression >,1 > > _5250_exprs,int32_t _5251_dim);
      void _initialize(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Expression >,1 > > _5250_exprs,int32_t _5251_dim);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _5253_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _5254_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _5255_xs) ;
      static  std::shared_ptr< monty::ndarray< int32_t,1 > > getshape(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Expression >,1 > > _5399_es,int32_t _5400_dim);
      virtual /* override */ std::string toString() ;
    }; // struct ExprStack;

    struct p_ExprInner : public ::mosek::fusion::p_BaseExpression
    {
      ExprInner * _pubthis;
      static mosek::fusion::p_ExprInner* _get_impl(mosek::fusion::ExprInner * _inst){ return static_cast< mosek::fusion::p_ExprInner* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprInner * _get_impl(mosek::fusion::ExprInner::t _inst) { return _get_impl(_inst.get()); }
      p_ExprInner(ExprInner * _pubthis);
      virtual ~p_ExprInner() { /* std::cout << "~p_ExprInner" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< double,1 > > vcof{};
      std::shared_ptr< monty::ndarray< int64_t,1 > > vsub{};
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};

      virtual void destroy();

      static ExprInner::t _new_ExprInner(monty::rc_ptr< ::mosek::fusion::Expression > _5414_expr3,std::shared_ptr< monty::ndarray< int64_t,1 > > _5415_vsub3,std::shared_ptr< monty::ndarray< double,1 > > _5416_vcof3);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _5414_expr3,std::shared_ptr< monty::ndarray< int64_t,1 > > _5415_vsub3,std::shared_ptr< monty::ndarray< double,1 > > _5416_vcof3);
      static ExprInner::t _new_ExprInner(monty::rc_ptr< ::mosek::fusion::Expression > _5422_expr2,std::shared_ptr< monty::ndarray< double,1 > > _5423_vcof2);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _5422_expr2,std::shared_ptr< monty::ndarray< double,1 > > _5423_vcof2);
      static ExprInner::t _new_ExprInner(monty::rc_ptr< ::mosek::fusion::Expression > _5425_expr1,std::shared_ptr< monty::ndarray< int32_t,2 > > _5426_vsub1,std::shared_ptr< monty::ndarray< double,1 > > _5427_vcof1);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _5425_expr1,std::shared_ptr< monty::ndarray< int32_t,2 > > _5426_vsub1,std::shared_ptr< monty::ndarray< double,1 > > _5427_vcof1);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _5428_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _5429_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _5430_xs) ;
      static  std::shared_ptr< monty::ndarray< int64_t,1 > > range(int32_t _5474_n);
      static  std::shared_ptr< monty::ndarray< int64_t,1 > > convert(std::shared_ptr< monty::ndarray< int32_t,1 > > _5476_shape,std::shared_ptr< monty::ndarray< int32_t,2 > > _5477_vsub);
      virtual /* override */ std::string toString() ;
    }; // struct ExprInner;

    struct p_ExprMulDiagRight : public ::mosek::fusion::p_BaseExpression
    {
      ExprMulDiagRight * _pubthis;
      static mosek::fusion::p_ExprMulDiagRight* _get_impl(mosek::fusion::ExprMulDiagRight * _inst){ return static_cast< mosek::fusion::p_ExprMulDiagRight* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprMulDiagRight * _get_impl(mosek::fusion::ExprMulDiagRight::t _inst) { return _get_impl(_inst.get()); }
      p_ExprMulDiagRight(ExprMulDiagRight * _pubthis);
      virtual ~p_ExprMulDiagRight() { /* std::cout << "~p_ExprMulDiagRight" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};
      std::shared_ptr< monty::ndarray< double,1 > > mval{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > msubj{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > msubi{};
      int32_t mdim1{};
      int32_t mdim0{};

      virtual void destroy();

      static ExprMulDiagRight::t _new_ExprMulDiagRight(int32_t _5484_mdim0,int32_t _5485_mdim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _5486_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _5487_msubj,std::shared_ptr< monty::ndarray< double,1 > > _5488_mval,monty::rc_ptr< ::mosek::fusion::Expression > _5489_expr,int32_t _5490_validated);
      void _initialize(int32_t _5484_mdim0,int32_t _5485_mdim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _5486_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _5487_msubj,std::shared_ptr< monty::ndarray< double,1 > > _5488_mval,monty::rc_ptr< ::mosek::fusion::Expression > _5489_expr,int32_t _5490_validated);
      static ExprMulDiagRight::t _new_ExprMulDiagRight(int32_t _5491_mdim0,int32_t _5492_mdim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _5493_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _5494_msubj,std::shared_ptr< monty::ndarray< double,1 > > _5495_mval,monty::rc_ptr< ::mosek::fusion::Expression > _5496_expr);
      void _initialize(int32_t _5491_mdim0,int32_t _5492_mdim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _5493_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _5494_msubj,std::shared_ptr< monty::ndarray< double,1 > > _5495_mval,monty::rc_ptr< ::mosek::fusion::Expression > _5496_expr);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _5497_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _5498_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _5499_xs) ;
      static  int32_t validate(int32_t _5578_mdim0,int32_t _5579_mdim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _5580_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _5581_msubj,std::shared_ptr< monty::ndarray< double,1 > > _5582_mval,monty::rc_ptr< ::mosek::fusion::Expression > _5583_expr);
      virtual /* override */ std::string toString() ;
    }; // struct ExprMulDiagRight;

    struct p_ExprMulDiagLeft : public ::mosek::fusion::p_BaseExpression
    {
      ExprMulDiagLeft * _pubthis;
      static mosek::fusion::p_ExprMulDiagLeft* _get_impl(mosek::fusion::ExprMulDiagLeft * _inst){ return static_cast< mosek::fusion::p_ExprMulDiagLeft* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprMulDiagLeft * _get_impl(mosek::fusion::ExprMulDiagLeft::t _inst) { return _get_impl(_inst.get()); }
      p_ExprMulDiagLeft(ExprMulDiagLeft * _pubthis);
      virtual ~p_ExprMulDiagLeft() { /* std::cout << "~p_ExprMulDiagLeft" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};
      std::shared_ptr< monty::ndarray< double,1 > > mval{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > msubj{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > msubi{};
      int32_t mdim1{};
      int32_t mdim0{};

      virtual void destroy();

      static ExprMulDiagLeft::t _new_ExprMulDiagLeft(int32_t _5592_mdim0,int32_t _5593_mdim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _5594_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _5595_msubj,std::shared_ptr< monty::ndarray< double,1 > > _5596_mval,monty::rc_ptr< ::mosek::fusion::Expression > _5597_expr,int32_t _5598_validated);
      void _initialize(int32_t _5592_mdim0,int32_t _5593_mdim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _5594_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _5595_msubj,std::shared_ptr< monty::ndarray< double,1 > > _5596_mval,monty::rc_ptr< ::mosek::fusion::Expression > _5597_expr,int32_t _5598_validated);
      static ExprMulDiagLeft::t _new_ExprMulDiagLeft(int32_t _5599_mdim0,int32_t _5600_mdim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _5601_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _5602_msubj,std::shared_ptr< monty::ndarray< double,1 > > _5603_mval,monty::rc_ptr< ::mosek::fusion::Expression > _5604_expr);
      void _initialize(int32_t _5599_mdim0,int32_t _5600_mdim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _5601_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _5602_msubj,std::shared_ptr< monty::ndarray< double,1 > > _5603_mval,monty::rc_ptr< ::mosek::fusion::Expression > _5604_expr);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _5605_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _5606_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _5607_xs) ;
      static  int32_t validate(int32_t _5705_mdim0,int32_t _5706_mdim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _5707_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _5708_msubj,std::shared_ptr< monty::ndarray< double,1 > > _5709_mval,monty::rc_ptr< ::mosek::fusion::Expression > _5710_expr);
      virtual /* override */ std::string toString() ;
    }; // struct ExprMulDiagLeft;

    struct p_ExprMulElement : public ::mosek::fusion::p_BaseExpression
    {
      ExprMulElement * _pubthis;
      static mosek::fusion::p_ExprMulElement* _get_impl(mosek::fusion::ExprMulElement * _inst){ return static_cast< mosek::fusion::p_ExprMulElement* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprMulElement * _get_impl(mosek::fusion::ExprMulElement::t _inst) { return _get_impl(_inst.get()); }
      p_ExprMulElement(ExprMulElement * _pubthis);
      virtual ~p_ExprMulElement() { /* std::cout << "~p_ExprMulElement" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};
      std::shared_ptr< monty::ndarray< int64_t,1 > > msp{};
      std::shared_ptr< monty::ndarray< double,1 > > mcof{};

      virtual void destroy();

      static ExprMulElement::t _new_ExprMulElement(std::shared_ptr< monty::ndarray< double,1 > > _5719_mcof,std::shared_ptr< monty::ndarray< int64_t,1 > > _5720_msp,monty::rc_ptr< ::mosek::fusion::Expression > _5721_expr);
      void _initialize(std::shared_ptr< monty::ndarray< double,1 > > _5719_mcof,std::shared_ptr< monty::ndarray< int64_t,1 > > _5720_msp,monty::rc_ptr< ::mosek::fusion::Expression > _5721_expr);
      static ExprMulElement::t _new_ExprMulElement(std::shared_ptr< monty::ndarray< double,1 > > _5728_cof,std::shared_ptr< monty::ndarray< int64_t,1 > > _5729_msp,monty::rc_ptr< ::mosek::fusion::Expression > _5730_expr,int32_t _5731_validated);
      void _initialize(std::shared_ptr< monty::ndarray< double,1 > > _5728_cof,std::shared_ptr< monty::ndarray< int64_t,1 > > _5729_msp,monty::rc_ptr< ::mosek::fusion::Expression > _5730_expr,int32_t _5731_validated);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _5732_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _5733_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _5734_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprMulElement;

    struct p_ExprMulScalarConst : public ::mosek::fusion::p_BaseExpression
    {
      ExprMulScalarConst * _pubthis;
      static mosek::fusion::p_ExprMulScalarConst* _get_impl(mosek::fusion::ExprMulScalarConst * _inst){ return static_cast< mosek::fusion::p_ExprMulScalarConst* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprMulScalarConst * _get_impl(mosek::fusion::ExprMulScalarConst::t _inst) { return _get_impl(_inst.get()); }
      p_ExprMulScalarConst(ExprMulScalarConst * _pubthis);
      virtual ~p_ExprMulScalarConst() { /* std::cout << "~p_ExprMulScalarConst" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};
      double c{};

      virtual void destroy();

      static ExprMulScalarConst::t _new_ExprMulScalarConst(double _5792_c,monty::rc_ptr< ::mosek::fusion::Expression > _5793_expr);
      void _initialize(double _5792_c,monty::rc_ptr< ::mosek::fusion::Expression > _5793_expr);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _5794_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _5795_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _5796_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprMulScalarConst;

    struct p_ExprScalarMul : public ::mosek::fusion::p_BaseExpression
    {
      ExprScalarMul * _pubthis;
      static mosek::fusion::p_ExprScalarMul* _get_impl(mosek::fusion::ExprScalarMul * _inst){ return static_cast< mosek::fusion::p_ExprScalarMul* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprScalarMul * _get_impl(mosek::fusion::ExprScalarMul::t _inst) { return _get_impl(_inst.get()); }
      p_ExprScalarMul(ExprScalarMul * _pubthis);
      virtual ~p_ExprScalarMul() { /* std::cout << "~p_ExprScalarMul" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};
      std::shared_ptr< monty::ndarray< double,1 > > mval{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > msubj{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > msubi{};
      int32_t mdim1{};
      int32_t mdim0{};

      virtual void destroy();

      static ExprScalarMul::t _new_ExprScalarMul(int32_t _5834_mdim0,int32_t _5835_mdim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _5836_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _5837_msubj,std::shared_ptr< monty::ndarray< double,1 > > _5838_mval,monty::rc_ptr< ::mosek::fusion::Expression > _5839_expr,int32_t _5840_validated);
      void _initialize(int32_t _5834_mdim0,int32_t _5835_mdim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _5836_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _5837_msubj,std::shared_ptr< monty::ndarray< double,1 > > _5838_mval,monty::rc_ptr< ::mosek::fusion::Expression > _5839_expr,int32_t _5840_validated);
      static ExprScalarMul::t _new_ExprScalarMul(int32_t _5841_mdim0,int32_t _5842_mdim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _5843_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _5844_msubj,std::shared_ptr< monty::ndarray< double,1 > > _5845_mval,monty::rc_ptr< ::mosek::fusion::Expression > _5846_expr);
      void _initialize(int32_t _5841_mdim0,int32_t _5842_mdim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _5843_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _5844_msubj,std::shared_ptr< monty::ndarray< double,1 > > _5845_mval,monty::rc_ptr< ::mosek::fusion::Expression > _5846_expr);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _5847_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _5848_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _5849_xs) ;
      static  int32_t validate(int32_t _5885_mdim0,int32_t _5886_mdim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _5887_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _5888_msubj,std::shared_ptr< monty::ndarray< double,1 > > _5889_mval,monty::rc_ptr< ::mosek::fusion::Expression > _5890_expr);
      virtual /* override */ std::string toString() ;
    }; // struct ExprScalarMul;

    struct p_ExprCrossDot : public ::mosek::fusion::p_BaseExpression
    {
      ExprCrossDot * _pubthis;
      static mosek::fusion::p_ExprCrossDot* _get_impl(mosek::fusion::ExprCrossDot * _inst){ return static_cast< mosek::fusion::p_ExprCrossDot* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprCrossDot * _get_impl(mosek::fusion::ExprCrossDot::t _inst) { return _get_impl(_inst.get()); }
      p_ExprCrossDot(ExprCrossDot * _pubthis);
      virtual ~p_ExprCrossDot() { /* std::cout << "~p_ExprCrossDot" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};
      std::shared_ptr< monty::ndarray< double,1 > > mval{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > msubj{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > msubi{};
      int32_t mdim1{};
      int32_t mdim0{};

      virtual void destroy();

      static ExprCrossDot::t _new_ExprCrossDot(monty::rc_ptr< ::mosek::fusion::Matrix > _5897_m,monty::rc_ptr< ::mosek::fusion::Expression > _5898_expr);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Matrix > _5897_m,monty::rc_ptr< ::mosek::fusion::Expression > _5898_expr);
      static ExprCrossDot::t _new_ExprCrossDot(int32_t _5900_mdim0,int32_t _5901_mdim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _5902_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _5903_msubj,std::shared_ptr< monty::ndarray< double,1 > > _5904_mval,monty::rc_ptr< ::mosek::fusion::Expression > _5905_expr,int32_t _5906_validated);
      void _initialize(int32_t _5900_mdim0,int32_t _5901_mdim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _5902_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _5903_msubj,std::shared_ptr< monty::ndarray< double,1 > > _5904_mval,monty::rc_ptr< ::mosek::fusion::Expression > _5905_expr,int32_t _5906_validated);
      static ExprCrossDot::t _new_ExprCrossDot(int32_t _5907_mdim0__,int32_t _5908_mdim1__,std::shared_ptr< monty::ndarray< double,1 > > _5909_mval__,monty::rc_ptr< ::mosek::fusion::Expression > _5910_expr__);
      void _initialize(int32_t _5907_mdim0__,int32_t _5908_mdim1__,std::shared_ptr< monty::ndarray< double,1 > > _5909_mval__,monty::rc_ptr< ::mosek::fusion::Expression > _5910_expr__);
      static ExprCrossDot::t _new_ExprCrossDot(int32_t _5911_mdim0_,int32_t _5912_mdim1_,std::shared_ptr< monty::ndarray< int32_t,1 > > _5913_msubi_,std::shared_ptr< monty::ndarray< int32_t,1 > > _5914_msubj_,std::shared_ptr< monty::ndarray< double,1 > > _5915_mval_,monty::rc_ptr< ::mosek::fusion::Expression > _5916_expr_);
      void _initialize(int32_t _5911_mdim0_,int32_t _5912_mdim1_,std::shared_ptr< monty::ndarray< int32_t,1 > > _5913_msubi_,std::shared_ptr< monty::ndarray< int32_t,1 > > _5914_msubj_,std::shared_ptr< monty::ndarray< double,1 > > _5915_mval_,monty::rc_ptr< ::mosek::fusion::Expression > _5916_expr_);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _5917_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _5918_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _5919_xs) ;
      static  std::shared_ptr< monty::ndarray< int32_t,1 > > computeshape(int32_t _6062_d0,int32_t _6063_d1,std::shared_ptr< monty::ndarray< int32_t,1 > > _6064_ds);
      static  int32_t validate(int32_t _6065_mdim0,int32_t _6066_mdim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _6067_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _6068_msubj,std::shared_ptr< monty::ndarray< double,1 > > _6069_mval,monty::rc_ptr< ::mosek::fusion::Expression > _6070_expr);
      virtual /* override */ std::string toString() ;
    }; // struct ExprCrossDot;

    struct p_ExprMulVar : public ::mosek::fusion::p_BaseExpression
    {
      ExprMulVar * _pubthis;
      static mosek::fusion::p_ExprMulVar* _get_impl(mosek::fusion::ExprMulVar * _inst){ return static_cast< mosek::fusion::p_ExprMulVar* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprMulVar * _get_impl(mosek::fusion::ExprMulVar::t _inst) { return _get_impl(_inst.get()); }
      p_ExprMulVar(ExprMulVar * _pubthis);
      virtual ~p_ExprMulVar() { /* std::cout << "~p_ExprMulVar" << std::endl;*/ };
      bool left{};
      monty::rc_ptr< ::mosek::fusion::Variable > x{};
      std::shared_ptr< monty::ndarray< double,1 > > mcof{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > msubj{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > msubi{};
      int32_t mdimj{};
      int32_t mdimi{};

      virtual void destroy();

      static ExprMulVar::t _new_ExprMulVar(bool _6078_left,int32_t _6079_mdimi,int32_t _6080_mdimj,std::shared_ptr< monty::ndarray< int32_t,1 > > _6081_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _6082_msubj,std::shared_ptr< monty::ndarray< double,1 > > _6083_mcof,monty::rc_ptr< ::mosek::fusion::Variable > _6084_x);
      void _initialize(bool _6078_left,int32_t _6079_mdimi,int32_t _6080_mdimj,std::shared_ptr< monty::ndarray< int32_t,1 > > _6081_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _6082_msubj,std::shared_ptr< monty::ndarray< double,1 > > _6083_mcof,monty::rc_ptr< ::mosek::fusion::Variable > _6084_x);
      static ExprMulVar::t _new_ExprMulVar(bool _6087_left,int32_t _6088_mdimi,int32_t _6089_mdimj,std::shared_ptr< monty::ndarray< int32_t,1 > > _6090_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _6091_msubj,std::shared_ptr< monty::ndarray< double,1 > > _6092_mcof,monty::rc_ptr< ::mosek::fusion::Variable > _6093_x,int32_t _6094_unchecked_);
      void _initialize(bool _6087_left,int32_t _6088_mdimi,int32_t _6089_mdimj,std::shared_ptr< monty::ndarray< int32_t,1 > > _6090_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _6091_msubj,std::shared_ptr< monty::ndarray< double,1 > > _6092_mcof,monty::rc_ptr< ::mosek::fusion::Variable > _6093_x,int32_t _6094_unchecked_);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _6095_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _6096_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _6097_xs) ;
      virtual void eval_right(monty::rc_ptr< ::mosek::fusion::WorkStack > _6098_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _6099_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _6100_xs) ;
      virtual void eval_left(monty::rc_ptr< ::mosek::fusion::WorkStack > _6205_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _6206_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _6207_xs) ;
      virtual void validate(int32_t _6280_mdimi,int32_t _6281_mdimj,std::shared_ptr< monty::ndarray< int32_t,1 > > _6282_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _6283_msubj,std::shared_ptr< monty::ndarray< double,1 > > _6284_mcof) ;
      static  std::shared_ptr< monty::ndarray< int32_t,1 > > resshape(int32_t _6288_mdimi,int32_t _6289_mdimj,std::shared_ptr< monty::ndarray< int32_t,1 > > _6290_xshape,bool _6291_left);
      virtual /* override */ std::string toString() ;
    }; // struct ExprMulVar;

    struct p_ExprMulScalarVar : public ::mosek::fusion::p_BaseExpression
    {
      ExprMulScalarVar * _pubthis;
      static mosek::fusion::p_ExprMulScalarVar* _get_impl(mosek::fusion::ExprMulScalarVar * _inst){ return static_cast< mosek::fusion::p_ExprMulScalarVar* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprMulScalarVar * _get_impl(mosek::fusion::ExprMulScalarVar::t _inst) { return _get_impl(_inst.get()); }
      p_ExprMulScalarVar(ExprMulScalarVar * _pubthis);
      virtual ~p_ExprMulScalarVar() { /* std::cout << "~p_ExprMulScalarVar" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Variable > x{};
      std::shared_ptr< monty::ndarray< double,1 > > mcof{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > msubj{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > msubi{};
      int32_t mdimj{};
      int32_t mdimi{};

      virtual void destroy();

      static ExprMulScalarVar::t _new_ExprMulScalarVar(int32_t _6292_mdimi,int32_t _6293_mdimj,std::shared_ptr< monty::ndarray< int32_t,1 > > _6294_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _6295_msubj,std::shared_ptr< monty::ndarray< double,1 > > _6296_mcof,monty::rc_ptr< ::mosek::fusion::Variable > _6297_x);
      void _initialize(int32_t _6292_mdimi,int32_t _6293_mdimj,std::shared_ptr< monty::ndarray< int32_t,1 > > _6294_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _6295_msubj,std::shared_ptr< monty::ndarray< double,1 > > _6296_mcof,monty::rc_ptr< ::mosek::fusion::Variable > _6297_x);
      static ExprMulScalarVar::t _new_ExprMulScalarVar(int32_t _6302_mdimi,int32_t _6303_mdimj,std::shared_ptr< monty::ndarray< int32_t,1 > > _6304_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _6305_msubj,std::shared_ptr< monty::ndarray< double,1 > > _6306_mcof,monty::rc_ptr< ::mosek::fusion::Variable > _6307_x,int32_t _6308_unchecked_);
      void _initialize(int32_t _6302_mdimi,int32_t _6303_mdimj,std::shared_ptr< monty::ndarray< int32_t,1 > > _6304_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _6305_msubj,std::shared_ptr< monty::ndarray< double,1 > > _6306_mcof,monty::rc_ptr< ::mosek::fusion::Variable > _6307_x,int32_t _6308_unchecked_);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _6309_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _6310_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _6311_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprMulScalarVar;

    struct p_ExprMulVarScalarConst : public ::mosek::fusion::p_BaseExpression
    {
      ExprMulVarScalarConst * _pubthis;
      static mosek::fusion::p_ExprMulVarScalarConst* _get_impl(mosek::fusion::ExprMulVarScalarConst * _inst){ return static_cast< mosek::fusion::p_ExprMulVarScalarConst* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprMulVarScalarConst * _get_impl(mosek::fusion::ExprMulVarScalarConst::t _inst) { return _get_impl(_inst.get()); }
      p_ExprMulVarScalarConst(ExprMulVarScalarConst * _pubthis);
      virtual ~p_ExprMulVarScalarConst() { /* std::cout << "~p_ExprMulVarScalarConst" << std::endl;*/ };
      double c{};
      monty::rc_ptr< ::mosek::fusion::Variable > x{};

      virtual void destroy();

      static ExprMulVarScalarConst::t _new_ExprMulVarScalarConst(monty::rc_ptr< ::mosek::fusion::Variable > _6328_x,double _6329_c);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Variable > _6328_x,double _6329_c);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _6330_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _6331_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _6332_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprMulVarScalarConst;

    struct p_ExprSumReduceEnd : public ::mosek::fusion::p_BaseExpression
    {
      ExprSumReduceEnd * _pubthis;
      static mosek::fusion::p_ExprSumReduceEnd* _get_impl(mosek::fusion::ExprSumReduceEnd * _inst){ return static_cast< mosek::fusion::p_ExprSumReduceEnd* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprSumReduceEnd * _get_impl(mosek::fusion::ExprSumReduceEnd::t _inst) { return _get_impl(_inst.get()); }
      p_ExprSumReduceEnd(ExprSumReduceEnd * _pubthis);
      virtual ~p_ExprSumReduceEnd() { /* std::cout << "~p_ExprSumReduceEnd" << std::endl;*/ };
      int32_t ndim{};
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};

      virtual void destroy();

      static ExprSumReduceEnd::t _new_ExprSumReduceEnd(int32_t _6349_ndim,monty::rc_ptr< ::mosek::fusion::Expression > _6350_expr);
      void _initialize(int32_t _6349_ndim,monty::rc_ptr< ::mosek::fusion::Expression > _6350_expr);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _6351_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _6352_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _6353_xs) ;
      static  std::shared_ptr< monty::ndarray< int32_t,1 > > computeShape(int32_t _6419_ndim,std::shared_ptr< monty::ndarray< int32_t,1 > > _6420_shape);
      virtual /* override */ std::string toString() ;
    }; // struct ExprSumReduceEnd;

    struct p_ExprAdd : public ::mosek::fusion::p_BaseExpression
    {
      ExprAdd * _pubthis;
      static mosek::fusion::p_ExprAdd* _get_impl(mosek::fusion::ExprAdd * _inst){ return static_cast< mosek::fusion::p_ExprAdd* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprAdd * _get_impl(mosek::fusion::ExprAdd::t _inst) { return _get_impl(_inst.get()); }
      p_ExprAdd(ExprAdd * _pubthis);
      virtual ~p_ExprAdd() { /* std::cout << "~p_ExprAdd" << std::endl;*/ };
      double m2{};
      double m1{};
      monty::rc_ptr< ::mosek::fusion::Expression > e2{};
      monty::rc_ptr< ::mosek::fusion::Expression > e1{};

      virtual void destroy();

      static ExprAdd::t _new_ExprAdd(monty::rc_ptr< ::mosek::fusion::Expression > _6423_e1,monty::rc_ptr< ::mosek::fusion::Expression > _6424_e2,double _6425_m1,double _6426_m2);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _6423_e1,monty::rc_ptr< ::mosek::fusion::Expression > _6424_e2,double _6425_m1,double _6426_m2);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _6428_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _6429_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _6430_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprAdd;

    struct p_ExprWSum : public ::mosek::fusion::p_BaseExpression
    {
      ExprWSum * _pubthis;
      static mosek::fusion::p_ExprWSum* _get_impl(mosek::fusion::ExprWSum * _inst){ return static_cast< mosek::fusion::p_ExprWSum* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprWSum * _get_impl(mosek::fusion::ExprWSum::t _inst) { return _get_impl(_inst.get()); }
      p_ExprWSum(ExprWSum * _pubthis);
      virtual ~p_ExprWSum() { /* std::cout << "~p_ExprWSum" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< double,1 > > w{};
      std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Expression >,1 > > es{};

      virtual void destroy();

      static ExprWSum::t _new_ExprWSum(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Expression >,1 > > _6560_es,std::shared_ptr< monty::ndarray< double,1 > > _6561_w);
      void _initialize(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Expression >,1 > > _6560_es,std::shared_ptr< monty::ndarray< double,1 > > _6561_w);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _6568_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _6569_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _6570_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprWSum;

    struct p_ExprSumReduce : public ::mosek::fusion::p_BaseExpression
    {
      ExprSumReduce * _pubthis;
      static mosek::fusion::p_ExprSumReduce* _get_impl(mosek::fusion::ExprSumReduce * _inst){ return static_cast< mosek::fusion::p_ExprSumReduce* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprSumReduce * _get_impl(mosek::fusion::ExprSumReduce::t _inst) { return _get_impl(_inst.get()); }
      p_ExprSumReduce(ExprSumReduce * _pubthis);
      virtual ~p_ExprSumReduce() { /* std::cout << "~p_ExprSumReduce" << std::endl;*/ };
      int32_t dim{};
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};

      virtual void destroy();

      static ExprSumReduce::t _new_ExprSumReduce(int32_t _6664_dim,monty::rc_ptr< ::mosek::fusion::Expression > _6665_expr);
      void _initialize(int32_t _6664_dim,monty::rc_ptr< ::mosek::fusion::Expression > _6665_expr);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _6667_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _6668_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _6669_xs) ;
      static  std::shared_ptr< monty::ndarray< int32_t,1 > > computeShape(int32_t _6785_dim,std::shared_ptr< monty::ndarray< int32_t,1 > > _6786_shape);
      virtual /* override */ std::string toString() ;
    }; // struct ExprSumReduce;

    struct p_ExprScaleVecPSD : public ::mosek::fusion::p_BaseExpression
    {
      ExprScaleVecPSD * _pubthis;
      static mosek::fusion::p_ExprScaleVecPSD* _get_impl(mosek::fusion::ExprScaleVecPSD * _inst){ return static_cast< mosek::fusion::p_ExprScaleVecPSD* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprScaleVecPSD * _get_impl(mosek::fusion::ExprScaleVecPSD::t _inst) { return _get_impl(_inst.get()); }
      p_ExprScaleVecPSD(ExprScaleVecPSD * _pubthis);
      virtual ~p_ExprScaleVecPSD() { /* std::cout << "~p_ExprScaleVecPSD" << std::endl;*/ };
      int32_t dim1{};
      int32_t dim0{};
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};

      virtual void destroy();

      static ExprScaleVecPSD::t _new_ExprScaleVecPSD(int32_t _6790_dim0,int32_t _6791_dim1,monty::rc_ptr< ::mosek::fusion::Expression > _6792_expr);
      void _initialize(int32_t _6790_dim0,int32_t _6791_dim1,monty::rc_ptr< ::mosek::fusion::Expression > _6792_expr);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _6793_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _6794_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _6795_xs) ;
    }; // struct ExprScaleVecPSD;

    struct p_ExprDenseTril : public ::mosek::fusion::p_BaseExpression
    {
      ExprDenseTril * _pubthis;
      static mosek::fusion::p_ExprDenseTril* _get_impl(mosek::fusion::ExprDenseTril * _inst){ return static_cast< mosek::fusion::p_ExprDenseTril* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprDenseTril * _get_impl(mosek::fusion::ExprDenseTril::t _inst) { return _get_impl(_inst.get()); }
      p_ExprDenseTril(ExprDenseTril * _pubthis);
      virtual ~p_ExprDenseTril() { /* std::cout << "~p_ExprDenseTril" << std::endl;*/ };
      int32_t dim1{};
      int32_t dim0{};
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};

      virtual void destroy();

      static ExprDenseTril::t _new_ExprDenseTril(int32_t _6869_dim0,int32_t _6870_dim1,monty::rc_ptr< ::mosek::fusion::Expression > _6871_expr,int32_t _6872_unchecked_);
      void _initialize(int32_t _6869_dim0,int32_t _6870_dim1,monty::rc_ptr< ::mosek::fusion::Expression > _6871_expr,int32_t _6872_unchecked_);
      static ExprDenseTril::t _new_ExprDenseTril(int32_t _6873_dim0_,int32_t _6874_dim1_,monty::rc_ptr< ::mosek::fusion::Expression > _6875_expr);
      void _initialize(int32_t _6873_dim0_,int32_t _6874_dim1_,monty::rc_ptr< ::mosek::fusion::Expression > _6875_expr);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _6877_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _6878_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _6879_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprDenseTril;

    struct p_ExprDense : public ::mosek::fusion::p_BaseExpression
    {
      ExprDense * _pubthis;
      static mosek::fusion::p_ExprDense* _get_impl(mosek::fusion::ExprDense * _inst){ return static_cast< mosek::fusion::p_ExprDense* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprDense * _get_impl(mosek::fusion::ExprDense::t _inst) { return _get_impl(_inst.get()); }
      p_ExprDense(ExprDense * _pubthis);
      virtual ~p_ExprDense() { /* std::cout << "~p_ExprDense" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};

      virtual void destroy();

      static ExprDense::t _new_ExprDense(monty::rc_ptr< ::mosek::fusion::Expression > _6963_expr);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _6963_expr);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _6964_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _6965_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _6966_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprDense;

    struct p_ExprSymmetrize : public ::mosek::fusion::p_BaseExpression
    {
      ExprSymmetrize * _pubthis;
      static mosek::fusion::p_ExprSymmetrize* _get_impl(mosek::fusion::ExprSymmetrize * _inst){ return static_cast< mosek::fusion::p_ExprSymmetrize* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprSymmetrize * _get_impl(mosek::fusion::ExprSymmetrize::t _inst) { return _get_impl(_inst.get()); }
      p_ExprSymmetrize(ExprSymmetrize * _pubthis);
      virtual ~p_ExprSymmetrize() { /* std::cout << "~p_ExprSymmetrize" << std::endl;*/ };
      int32_t dim1{};
      int32_t dim0{};
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};

      virtual void destroy();

      static ExprSymmetrize::t _new_ExprSymmetrize(int32_t _7007_dim0,int32_t _7008_dim1,monty::rc_ptr< ::mosek::fusion::Expression > _7009_expr,int32_t _7010_unchecked_);
      void _initialize(int32_t _7007_dim0,int32_t _7008_dim1,monty::rc_ptr< ::mosek::fusion::Expression > _7009_expr,int32_t _7010_unchecked_);
      static ExprSymmetrize::t _new_ExprSymmetrize(int32_t _7011_dim0_,int32_t _7012_dim1_,monty::rc_ptr< ::mosek::fusion::Expression > _7013_expr);
      void _initialize(int32_t _7011_dim0_,int32_t _7012_dim1_,monty::rc_ptr< ::mosek::fusion::Expression > _7013_expr);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _7015_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _7016_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _7017_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprSymmetrize;

    struct p_ExprCondense : public ::mosek::fusion::p_BaseExpression
    {
      ExprCondense * _pubthis;
      static mosek::fusion::p_ExprCondense* _get_impl(mosek::fusion::ExprCondense * _inst){ return static_cast< mosek::fusion::p_ExprCondense* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprCondense * _get_impl(mosek::fusion::ExprCondense::t _inst) { return _get_impl(_inst.get()); }
      p_ExprCondense(ExprCondense * _pubthis);
      virtual ~p_ExprCondense() { /* std::cout << "~p_ExprCondense" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Expression > expr{};

      virtual void destroy();

      static ExprCondense::t _new_ExprCondense(monty::rc_ptr< ::mosek::fusion::Expression > _7159_expr);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _7159_expr);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _7160_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _7161_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _7162_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprCondense;

    struct p_ExprFromVar : public ::mosek::fusion::p_BaseExpression
    {
      ExprFromVar * _pubthis;
      static mosek::fusion::p_ExprFromVar* _get_impl(mosek::fusion::ExprFromVar * _inst){ return static_cast< mosek::fusion::p_ExprFromVar* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprFromVar * _get_impl(mosek::fusion::ExprFromVar::t _inst) { return _get_impl(_inst.get()); }
      p_ExprFromVar(ExprFromVar * _pubthis);
      virtual ~p_ExprFromVar() { /* std::cout << "~p_ExprFromVar" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Variable > x{};

      virtual void destroy();

      static ExprFromVar::t _new_ExprFromVar(monty::rc_ptr< ::mosek::fusion::Variable > _7166_x);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Variable > _7166_x);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _7167_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _7168_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _7169_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprFromVar;

    struct p_ExprReshape : public ::mosek::fusion::p_BaseExpression
    {
      ExprReshape * _pubthis;
      static mosek::fusion::p_ExprReshape* _get_impl(mosek::fusion::ExprReshape * _inst){ return static_cast< mosek::fusion::p_ExprReshape* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_ExprReshape * _get_impl(mosek::fusion::ExprReshape::t _inst) { return _get_impl(_inst.get()); }
      p_ExprReshape(ExprReshape * _pubthis);
      virtual ~p_ExprReshape() { /* std::cout << "~p_ExprReshape" << std::endl;*/ };
      monty::rc_ptr< ::mosek::fusion::Expression > e{};

      virtual void destroy();

      static ExprReshape::t _new_ExprReshape(std::shared_ptr< monty::ndarray< int32_t,1 > > _7186_shape,monty::rc_ptr< ::mosek::fusion::Expression > _7187_e);
      void _initialize(std::shared_ptr< monty::ndarray< int32_t,1 > > _7186_shape,monty::rc_ptr< ::mosek::fusion::Expression > _7187_e);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _7189_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _7190_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _7191_xs) ;
      virtual /* override */ std::string toString() ;
    }; // struct ExprReshape;

    struct p_Expr : public ::mosek::fusion::p_BaseExpression
    {
      Expr * _pubthis;
      static mosek::fusion::p_Expr* _get_impl(mosek::fusion::Expr * _inst){ return static_cast< mosek::fusion::p_Expr* >(mosek::fusion::p_BaseExpression::_get_impl(_inst)); }
      static mosek::fusion::p_Expr * _get_impl(mosek::fusion::Expr::t _inst) { return _get_impl(_inst.get()); }
      p_Expr(Expr * _pubthis);
      virtual ~p_Expr() { /* std::cout << "~p_Expr" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< int64_t,1 > > inst{};
      std::shared_ptr< monty::ndarray< double,1 > > cof_v{};
      std::shared_ptr< monty::ndarray< int64_t,1 > > subj{};
      std::shared_ptr< monty::ndarray< int64_t,1 > > ptrb{};
      std::shared_ptr< monty::ndarray< double,1 > > bfix{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > shape{};

      virtual void destroy();

      static Expr::t _new_Expr(std::shared_ptr< monty::ndarray< int64_t,1 > > _7347_ptrb,std::shared_ptr< monty::ndarray< int64_t,1 > > _7348_subj,std::shared_ptr< monty::ndarray< double,1 > > _7349_cof,std::shared_ptr< monty::ndarray< double,1 > > _7350_bfix,std::shared_ptr< monty::ndarray< int32_t,1 > > _7351_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _7352_inst);
      void _initialize(std::shared_ptr< monty::ndarray< int64_t,1 > > _7347_ptrb,std::shared_ptr< monty::ndarray< int64_t,1 > > _7348_subj,std::shared_ptr< monty::ndarray< double,1 > > _7349_cof,std::shared_ptr< monty::ndarray< double,1 > > _7350_bfix,std::shared_ptr< monty::ndarray< int32_t,1 > > _7351_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _7352_inst);
      static Expr::t _new_Expr(std::shared_ptr< monty::ndarray< int64_t,1 > > _7363_ptrb,std::shared_ptr< monty::ndarray< int64_t,1 > > _7364_subj,std::shared_ptr< monty::ndarray< double,1 > > _7365_cof,std::shared_ptr< monty::ndarray< double,1 > > _7366_bfix,std::shared_ptr< monty::ndarray< int32_t,1 > > _7367_shp,std::shared_ptr< monty::ndarray< int64_t,1 > > _7368_inst,int32_t _7369_unchecked_);
      void _initialize(std::shared_ptr< monty::ndarray< int64_t,1 > > _7363_ptrb,std::shared_ptr< monty::ndarray< int64_t,1 > > _7364_subj,std::shared_ptr< monty::ndarray< double,1 > > _7365_cof,std::shared_ptr< monty::ndarray< double,1 > > _7366_bfix,std::shared_ptr< monty::ndarray< int32_t,1 > > _7367_shp,std::shared_ptr< monty::ndarray< int64_t,1 > > _7368_inst,int32_t _7369_unchecked_);
      static Expr::t _new_Expr(monty::rc_ptr< ::mosek::fusion::Expression > _7370_e);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Expression > _7370_e);
      virtual int64_t prod(std::shared_ptr< monty::ndarray< int32_t,1 > > _7395_vals) ;
      static  std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Variable >,1 > > varstack(std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Variable >,1 > >,1 > > _7398_vs);
      static  std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Variable >,1 > > varstack(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Variable >,1 > > _7401_v1,std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Variable >,1 > > _7402_v2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > condense(monty::rc_ptr< ::mosek::fusion::Expression > _7406_e);
      static  monty::rc_ptr< ::mosek::fusion::Expression > flatten(monty::rc_ptr< ::mosek::fusion::Expression > _7407_e);
      static  monty::rc_ptr< ::mosek::fusion::Expression > reshape(monty::rc_ptr< ::mosek::fusion::Expression > _7409_e,int32_t _7410_dimi,int32_t _7411_dimj);
      static  monty::rc_ptr< ::mosek::fusion::Expression > reshape(monty::rc_ptr< ::mosek::fusion::Expression > _7412_e,int32_t _7413_size);
      static  monty::rc_ptr< ::mosek::fusion::Expression > reshape(monty::rc_ptr< ::mosek::fusion::Expression > _7414_e,std::shared_ptr< monty::ndarray< int32_t,1 > > _7415_newshape);
      static  monty::rc_ptr< ::mosek::fusion::Expression > zeros(std::shared_ptr< monty::ndarray< int32_t,1 > > _7416_shp);
      static  monty::rc_ptr< ::mosek::fusion::Expression > zeros(int32_t _7417_size);
      static  monty::rc_ptr< ::mosek::fusion::Expression > ones();
      static  monty::rc_ptr< ::mosek::fusion::Expression > ones(std::shared_ptr< monty::ndarray< int32_t,1 > > _7418_shp,std::shared_ptr< monty::ndarray< int32_t,2 > > _7419_sparsity);
      static  monty::rc_ptr< ::mosek::fusion::Expression > ones(std::shared_ptr< monty::ndarray< int32_t,1 > > _7420_shp);
      static  monty::rc_ptr< ::mosek::fusion::Expression > ones(int32_t _7421_size);
      static  monty::rc_ptr< ::mosek::fusion::Expression > constTerm(monty::rc_ptr< ::mosek::fusion::NDSparseArray > _7422_nda);
      static  monty::rc_ptr< ::mosek::fusion::Expression > constTerm(monty::rc_ptr< ::mosek::fusion::Matrix > _7423_m);
      static  monty::rc_ptr< ::mosek::fusion::Expression > constTerm(double _7432_val);
      static  monty::rc_ptr< ::mosek::fusion::Expression > constTerm(std::shared_ptr< monty::ndarray< int32_t,1 > > _7433_shp,std::shared_ptr< monty::ndarray< int32_t,2 > > _7434_sparsity,double _7435_val);
      static  monty::rc_ptr< ::mosek::fusion::Expression > constTerm(std::shared_ptr< monty::ndarray< int32_t,1 > > _7443_shp,std::shared_ptr< monty::ndarray< int32_t,2 > > _7444_sparsity,std::shared_ptr< monty::ndarray< double,1 > > _7445_vals1);
      static  monty::rc_ptr< ::mosek::fusion::Expression > constTerm(std::shared_ptr< monty::ndarray< int32_t,1 > > _7453_shp,double _7454_val);
      static  monty::rc_ptr< ::mosek::fusion::Expression > constTerm(int32_t _7455_size,double _7456_val);
      static  monty::rc_ptr< ::mosek::fusion::Expression > constTerm(std::shared_ptr< monty::ndarray< double,2 > > _7458_vals2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > constTerm(std::shared_ptr< monty::ndarray< double,1 > > _7461_vals1);
      static  monty::rc_ptr< ::mosek::fusion::Expression > sum(monty::rc_ptr< ::mosek::fusion::Expression > _7462_expr,std::shared_ptr< monty::ndarray< int32_t,1 > > _7463_dims);
      static  monty::rc_ptr< ::mosek::fusion::Expression > sum(monty::rc_ptr< ::mosek::fusion::Expression > _7472_expr,int32_t _7473_dim);
      static  monty::rc_ptr< ::mosek::fusion::Expression > sum(monty::rc_ptr< ::mosek::fusion::Expression > _7474_expr);
      static  monty::rc_ptr< ::mosek::fusion::Expression > neg(monty::rc_ptr< ::mosek::fusion::Expression > _7475_e);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulDiag(bool _7476_left,monty::rc_ptr< ::mosek::fusion::Matrix > _7477_mx,monty::rc_ptr< ::mosek::fusion::Expression > _7478_expr);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulDiag(monty::rc_ptr< ::mosek::fusion::Variable > _7485_v,monty::rc_ptr< ::mosek::fusion::Parameter > _7486_p);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulDiag(monty::rc_ptr< ::mosek::fusion::Parameter > _7487_p,monty::rc_ptr< ::mosek::fusion::Variable > _7488_v);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulDiag(monty::rc_ptr< ::mosek::fusion::Expression > _7489_expr,monty::rc_ptr< ::mosek::fusion::Parameter > _7490_p);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulDiag(monty::rc_ptr< ::mosek::fusion::Parameter > _7491_p,monty::rc_ptr< ::mosek::fusion::Expression > _7492_expr);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulDiag(monty::rc_ptr< ::mosek::fusion::Variable > _7493_v,monty::rc_ptr< ::mosek::fusion::Matrix > _7494_mx);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulDiag(monty::rc_ptr< ::mosek::fusion::Matrix > _7495_mx,monty::rc_ptr< ::mosek::fusion::Variable > _7496_v);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulDiag(monty::rc_ptr< ::mosek::fusion::Expression > _7497_expr,monty::rc_ptr< ::mosek::fusion::Matrix > _7498_mx);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulDiag(monty::rc_ptr< ::mosek::fusion::Matrix > _7499_mx,monty::rc_ptr< ::mosek::fusion::Expression > _7500_expr);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulDiag(monty::rc_ptr< ::mosek::fusion::Variable > _7501_v,std::shared_ptr< monty::ndarray< double,2 > > _7502_a);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulDiag(monty::rc_ptr< ::mosek::fusion::Expression > _7509_expr,std::shared_ptr< monty::ndarray< double,2 > > _7510_a);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulDiag(std::shared_ptr< monty::ndarray< double,2 > > _7517_a,monty::rc_ptr< ::mosek::fusion::Variable > _7518_v);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulDiag(std::shared_ptr< monty::ndarray< double,2 > > _7525_a,monty::rc_ptr< ::mosek::fusion::Expression > _7526_expr);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulElm_(monty::rc_ptr< ::mosek::fusion::Matrix > _7533_m,monty::rc_ptr< ::mosek::fusion::Expression > _7534_e);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulElm_(std::shared_ptr< monty::ndarray< double,1 > > _7543_a,monty::rc_ptr< ::mosek::fusion::Expression > _7544_expr);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulElm_(monty::rc_ptr< ::mosek::fusion::NDSparseArray > _7546_spm,monty::rc_ptr< ::mosek::fusion::Expression > _7547_expr);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mul(monty::rc_ptr< ::mosek::fusion::Expression > _7550_expr,double _7551_c);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mul(double _7552_c,monty::rc_ptr< ::mosek::fusion::Expression > _7553_expr);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mul(monty::rc_ptr< ::mosek::fusion::Expression > _7554_expr,std::shared_ptr< monty::ndarray< double,1 > > _7555_a);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mul(std::shared_ptr< monty::ndarray< double,1 > > _7556_a,monty::rc_ptr< ::mosek::fusion::Expression > _7557_expr);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mul(monty::rc_ptr< ::mosek::fusion::Expression > _7558_expr,std::shared_ptr< monty::ndarray< double,2 > > _7559_a);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mul(std::shared_ptr< monty::ndarray< double,2 > > _7560_a,monty::rc_ptr< ::mosek::fusion::Expression > _7561_expr);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mul(monty::rc_ptr< ::mosek::fusion::Expression > _7562_expr,monty::rc_ptr< ::mosek::fusion::Matrix > _7563_mx);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mul(monty::rc_ptr< ::mosek::fusion::Matrix > _7564_mx,monty::rc_ptr< ::mosek::fusion::Expression > _7565_expr);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mul(bool _7566_left,std::shared_ptr< monty::ndarray< double,1 > > _7567_mx,monty::rc_ptr< ::mosek::fusion::Expression > _7568_e);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mul(bool _7575_left,std::shared_ptr< monty::ndarray< double,2 > > _7576_mx,monty::rc_ptr< ::mosek::fusion::Expression > _7577_e);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mul(bool _7593_left,monty::rc_ptr< ::mosek::fusion::Matrix > _7594_mx,monty::rc_ptr< ::mosek::fusion::Expression > _7595_e);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mul(monty::rc_ptr< ::mosek::fusion::Variable > _7617_v,std::shared_ptr< monty::ndarray< double,2 > > _7618_mx);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mul(std::shared_ptr< monty::ndarray< double,2 > > _7628_mx,monty::rc_ptr< ::mosek::fusion::Variable > _7629_v);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mul(monty::rc_ptr< ::mosek::fusion::Variable > _7639_v,monty::rc_ptr< ::mosek::fusion::Matrix > _7640_mx);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mul(monty::rc_ptr< ::mosek::fusion::Matrix > _7646_mx,monty::rc_ptr< ::mosek::fusion::Variable > _7647_v);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mul(bool _7653_left,int32_t _7654_mdimi,int32_t _7655_mdimj,std::shared_ptr< monty::ndarray< int32_t,1 > > _7656_msubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _7657_msubj,std::shared_ptr< monty::ndarray< double,1 > > _7658_mcof,monty::rc_ptr< ::mosek::fusion::Variable > _7659_v);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mul(monty::rc_ptr< ::mosek::fusion::Expression > _7661_expr,monty::rc_ptr< ::mosek::fusion::Parameter > _7662_p);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mul(monty::rc_ptr< ::mosek::fusion::Parameter > _7663_p,monty::rc_ptr< ::mosek::fusion::Expression > _7664_expr);
      static  monty::rc_ptr< ::mosek::fusion::Expression > dot(monty::rc_ptr< ::mosek::fusion::Expression > _7665_e,monty::rc_ptr< ::mosek::fusion::Matrix > _7666_m);
      static  monty::rc_ptr< ::mosek::fusion::Expression > dot(monty::rc_ptr< ::mosek::fusion::Expression > _7675_e,std::shared_ptr< monty::ndarray< double,2 > > _7676_c2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > dot(monty::rc_ptr< ::mosek::fusion::Expression > _7681_e,monty::rc_ptr< ::mosek::fusion::NDSparseArray > _7682_nda);
      static  monty::rc_ptr< ::mosek::fusion::Expression > dot(monty::rc_ptr< ::mosek::fusion::Expression > _7685_e,std::shared_ptr< monty::ndarray< double,1 > > _7686_c1);
      static  monty::rc_ptr< ::mosek::fusion::Expression > dot(monty::rc_ptr< ::mosek::fusion::Matrix > _7691_m,monty::rc_ptr< ::mosek::fusion::Expression > _7692_e);
      static  monty::rc_ptr< ::mosek::fusion::Expression > dot(monty::rc_ptr< ::mosek::fusion::NDSparseArray > _7693_nda,monty::rc_ptr< ::mosek::fusion::Expression > _7694_e);
      static  monty::rc_ptr< ::mosek::fusion::Expression > dot(std::shared_ptr< monty::ndarray< double,2 > > _7695_c2,monty::rc_ptr< ::mosek::fusion::Expression > _7696_e);
      static  monty::rc_ptr< ::mosek::fusion::Expression > dot(std::shared_ptr< monty::ndarray< double,1 > > _7697_c1,monty::rc_ptr< ::mosek::fusion::Expression > _7698_e);
      static  monty::rc_ptr< ::mosek::fusion::Expression > dot(monty::rc_ptr< ::mosek::fusion::Expression > _7699_e,monty::rc_ptr< ::mosek::fusion::Parameter > _7700_p);
      static  monty::rc_ptr< ::mosek::fusion::Expression > dot(monty::rc_ptr< ::mosek::fusion::Parameter > _7701_p,monty::rc_ptr< ::mosek::fusion::Expression > _7702_e);
      static  monty::rc_ptr< ::mosek::fusion::Expression > outer(monty::rc_ptr< ::mosek::fusion::Parameter > _7703_p,monty::rc_ptr< ::mosek::fusion::Expression > _7704_e);
      static  monty::rc_ptr< ::mosek::fusion::Expression > outer(monty::rc_ptr< ::mosek::fusion::Expression > _7707_e,monty::rc_ptr< ::mosek::fusion::Parameter > _7708_p);
      static  monty::rc_ptr< ::mosek::fusion::Expression > outer(monty::rc_ptr< ::mosek::fusion::Matrix > _7711_m,monty::rc_ptr< ::mosek::fusion::Expression > _7712_e);
      static  monty::rc_ptr< ::mosek::fusion::Expression > outer(monty::rc_ptr< ::mosek::fusion::Expression > _7714_e,monty::rc_ptr< ::mosek::fusion::Matrix > _7715_m);
      static  monty::rc_ptr< ::mosek::fusion::Expression > outer(std::shared_ptr< monty::ndarray< double,1 > > _7717_a,monty::rc_ptr< ::mosek::fusion::Expression > _7718_e);
      static  monty::rc_ptr< ::mosek::fusion::Expression > outer(monty::rc_ptr< ::mosek::fusion::Expression > _7720_e,std::shared_ptr< monty::ndarray< double,1 > > _7721_a);
      static  monty::rc_ptr< ::mosek::fusion::Expression > outer_(int32_t _7723_edim,std::shared_ptr< monty::ndarray< int64_t,1 > > _7724_eptrb,std::shared_ptr< monty::ndarray< int64_t,1 > > _7725_esubj,std::shared_ptr< monty::ndarray< double,1 > > _7726_ecof,std::shared_ptr< monty::ndarray< double,1 > > _7727_ebfix,std::shared_ptr< monty::ndarray< int64_t,1 > > _7728_einst,std::shared_ptr< monty::ndarray< double,1 > > _7729_a,std::shared_ptr< monty::ndarray< int32_t,1 > > _7730_sub,int32_t _7731_dim,bool _7732_transpose);
      static  monty::rc_ptr< ::mosek::fusion::Expression > outer_(monty::rc_ptr< ::mosek::fusion::Variable > _7762_v,int32_t _7763_vdim,std::shared_ptr< monty::ndarray< double,1 > > _7764_a,std::shared_ptr< monty::ndarray< int32_t,1 > > _7765_sub,int32_t _7766_dim,bool _7767_transpose);
      static  monty::rc_ptr< ::mosek::fusion::Expression > stack(std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Expression >,1 > >,1 > > _7784_exprs);
      static  monty::rc_ptr< ::mosek::fusion::Expression > vstack(double _7790_a1,double _7791_a2,double _7792_a3);
      static  monty::rc_ptr< ::mosek::fusion::Expression > vstack(double _7793_a1,double _7794_a2,monty::rc_ptr< ::mosek::fusion::Expression > _7795_e3);
      static  monty::rc_ptr< ::mosek::fusion::Expression > vstack(double _7796_a1,monty::rc_ptr< ::mosek::fusion::Expression > _7797_e2,double _7798_a3);
      static  monty::rc_ptr< ::mosek::fusion::Expression > vstack(double _7799_a1,monty::rc_ptr< ::mosek::fusion::Expression > _7800_e2,monty::rc_ptr< ::mosek::fusion::Expression > _7801_e3);
      static  monty::rc_ptr< ::mosek::fusion::Expression > vstack(monty::rc_ptr< ::mosek::fusion::Expression > _7802_e1,double _7803_a2,double _7804_a3);
      static  monty::rc_ptr< ::mosek::fusion::Expression > vstack(monty::rc_ptr< ::mosek::fusion::Expression > _7805_e1,double _7806_a2,monty::rc_ptr< ::mosek::fusion::Expression > _7807_e3);
      static  monty::rc_ptr< ::mosek::fusion::Expression > vstack(monty::rc_ptr< ::mosek::fusion::Expression > _7808_e1,monty::rc_ptr< ::mosek::fusion::Expression > _7809_e2,double _7810_a3);
      static  monty::rc_ptr< ::mosek::fusion::Expression > vstack(monty::rc_ptr< ::mosek::fusion::Expression > _7811_e1,monty::rc_ptr< ::mosek::fusion::Expression > _7812_e2,monty::rc_ptr< ::mosek::fusion::Expression > _7813_e3);
      static  monty::rc_ptr< ::mosek::fusion::Expression > vstack(double _7814_a1,monty::rc_ptr< ::mosek::fusion::Expression > _7815_e2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > vstack(monty::rc_ptr< ::mosek::fusion::Expression > _7816_e1,double _7817_a2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > vstack(monty::rc_ptr< ::mosek::fusion::Expression > _7818_e1,monty::rc_ptr< ::mosek::fusion::Expression > _7819_e2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > vstack(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Expression >,1 > > _7820_exprs);
      static  monty::rc_ptr< ::mosek::fusion::Expression > hstack(monty::rc_ptr< ::mosek::fusion::Expression > _7822_e1,monty::rc_ptr< ::mosek::fusion::Expression > _7823_e2,monty::rc_ptr< ::mosek::fusion::Expression > _7824_e3);
      static  monty::rc_ptr< ::mosek::fusion::Expression > hstack(monty::rc_ptr< ::mosek::fusion::Expression > _7825_e1,monty::rc_ptr< ::mosek::fusion::Expression > _7826_e2,double _7827_a3);
      static  monty::rc_ptr< ::mosek::fusion::Expression > hstack(monty::rc_ptr< ::mosek::fusion::Expression > _7828_e1,double _7829_a2,monty::rc_ptr< ::mosek::fusion::Expression > _7830_e3);
      static  monty::rc_ptr< ::mosek::fusion::Expression > hstack(monty::rc_ptr< ::mosek::fusion::Expression > _7831_e1,double _7832_a2,double _7833_a3);
      static  monty::rc_ptr< ::mosek::fusion::Expression > hstack(double _7834_a1,monty::rc_ptr< ::mosek::fusion::Expression > _7835_e2,monty::rc_ptr< ::mosek::fusion::Expression > _7836_e3);
      static  monty::rc_ptr< ::mosek::fusion::Expression > hstack(double _7837_a1,monty::rc_ptr< ::mosek::fusion::Expression > _7838_e2,double _7839_a3);
      static  monty::rc_ptr< ::mosek::fusion::Expression > hstack(double _7840_a1,double _7841_a2,monty::rc_ptr< ::mosek::fusion::Expression > _7842_e3);
      static  monty::rc_ptr< ::mosek::fusion::Expression > hstack(double _7843_a1,monty::rc_ptr< ::mosek::fusion::Expression > _7844_e2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > hstack(monty::rc_ptr< ::mosek::fusion::Expression > _7845_e1,double _7846_a2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > hstack(monty::rc_ptr< ::mosek::fusion::Expression > _7847_e1,monty::rc_ptr< ::mosek::fusion::Expression > _7848_e2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > hstack(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Expression >,1 > > _7849_exprs);
      static  monty::rc_ptr< ::mosek::fusion::Expression > stack(int32_t _7851_dim,monty::rc_ptr< ::mosek::fusion::Expression > _7852_e1,monty::rc_ptr< ::mosek::fusion::Expression > _7853_e2,monty::rc_ptr< ::mosek::fusion::Expression > _7854_e3);
      static  monty::rc_ptr< ::mosek::fusion::Expression > stack(int32_t _7855_dim,monty::rc_ptr< ::mosek::fusion::Expression > _7856_e1,monty::rc_ptr< ::mosek::fusion::Expression > _7857_e2,double _7858_a3);
      static  monty::rc_ptr< ::mosek::fusion::Expression > stack(int32_t _7859_dim,monty::rc_ptr< ::mosek::fusion::Expression > _7860_e1,double _7861_a2,monty::rc_ptr< ::mosek::fusion::Expression > _7862_e3);
      static  monty::rc_ptr< ::mosek::fusion::Expression > stack(int32_t _7863_dim,monty::rc_ptr< ::mosek::fusion::Expression > _7864_e1,double _7865_a2,double _7866_a3);
      static  monty::rc_ptr< ::mosek::fusion::Expression > stack(int32_t _7867_dim,double _7868_a1,monty::rc_ptr< ::mosek::fusion::Expression > _7869_e2,monty::rc_ptr< ::mosek::fusion::Expression > _7870_e3);
      static  monty::rc_ptr< ::mosek::fusion::Expression > stack(int32_t _7871_dim,double _7872_a1,monty::rc_ptr< ::mosek::fusion::Expression > _7873_e2,double _7874_a3);
      static  monty::rc_ptr< ::mosek::fusion::Expression > stack(int32_t _7875_dim,double _7876_a1,double _7877_a2,monty::rc_ptr< ::mosek::fusion::Expression > _7878_e1);
      static  monty::rc_ptr< ::mosek::fusion::Expression > stack(int32_t _7879_dim,double _7880_a1,monty::rc_ptr< ::mosek::fusion::Expression > _7881_e2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > stack(int32_t _7882_dim,monty::rc_ptr< ::mosek::fusion::Expression > _7883_e1,double _7884_a2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > stack(int32_t _7885_dim,monty::rc_ptr< ::mosek::fusion::Expression > _7886_e1,monty::rc_ptr< ::mosek::fusion::Expression > _7887_e2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > stack(int32_t _7888_dim,std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Expression >,1 > > _7889_exprs);
      static  monty::rc_ptr< ::mosek::fusion::Expression > stack_(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Expression >,1 > > _7890_exprs,int32_t _7891_dim);
      static  std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Expression >,1 > > promote(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Expression >,1 > > _7892_exprs,int32_t _7893_dim);
      static  monty::rc_ptr< ::mosek::fusion::Expression > repeat(monty::rc_ptr< ::mosek::fusion::Variable > _7906_x,int32_t _7907_n,int32_t _7908_d);
      static  monty::rc_ptr< ::mosek::fusion::Expression > repeat(monty::rc_ptr< ::mosek::fusion::Expression > _7909_e,int32_t _7910_n,int32_t _7911_d);
      static  monty::rc_ptr< ::mosek::fusion::Expression > add(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Expression >,1 > > _7914_exps);
      static  monty::rc_ptr< ::mosek::fusion::Expression > add(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Variable >,1 > > _7916_vs);
      static  monty::rc_ptr< ::mosek::fusion::Expression > add_(monty::rc_ptr< ::mosek::fusion::Expression > _7949_e1,double _7950_m1,monty::rc_ptr< ::mosek::fusion::Expression > _7951_e2,double _7952_m2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > transpose(monty::rc_ptr< ::mosek::fusion::Expression > _7963_e);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulElm(monty::rc_ptr< ::mosek::fusion::Matrix > _7964_m,monty::rc_ptr< ::mosek::fusion::Expression > _7965_expr);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulElm(monty::rc_ptr< ::mosek::fusion::NDSparseArray > _7966_spm,monty::rc_ptr< ::mosek::fusion::Expression > _7967_expr);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulElm(std::shared_ptr< monty::ndarray< double,2 > > _7968_a2,monty::rc_ptr< ::mosek::fusion::Expression > _7969_expr);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulElm(std::shared_ptr< monty::ndarray< double,1 > > _7970_a1,monty::rc_ptr< ::mosek::fusion::Expression > _7971_expr);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulElm(monty::rc_ptr< ::mosek::fusion::Expression > _7972_expr,monty::rc_ptr< ::mosek::fusion::Matrix > _7973_m);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulElm(monty::rc_ptr< ::mosek::fusion::Expression > _7974_expr,std::shared_ptr< monty::ndarray< double,2 > > _7975_a2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulElm(monty::rc_ptr< ::mosek::fusion::Expression > _7976_expr,std::shared_ptr< monty::ndarray< double,1 > > _7977_a1);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulElm(monty::rc_ptr< ::mosek::fusion::Expression > _7978_expr,monty::rc_ptr< ::mosek::fusion::NDSparseArray > _7979_spm);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulElm(monty::rc_ptr< ::mosek::fusion::Parameter > _7980_p,monty::rc_ptr< ::mosek::fusion::Expression > _7981_expr);
      static  monty::rc_ptr< ::mosek::fusion::Expression > mulElm(monty::rc_ptr< ::mosek::fusion::Expression > _7982_expr,monty::rc_ptr< ::mosek::fusion::Parameter > _7983_p);
      static  monty::rc_ptr< ::mosek::fusion::Expression > sub(monty::rc_ptr< ::mosek::fusion::NDSparseArray > _7984_n,monty::rc_ptr< ::mosek::fusion::Expression > _7985_e2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > sub(monty::rc_ptr< ::mosek::fusion::Expression > _7986_e1,monty::rc_ptr< ::mosek::fusion::NDSparseArray > _7987_n);
      static  monty::rc_ptr< ::mosek::fusion::Expression > sub(monty::rc_ptr< ::mosek::fusion::Matrix > _7988_m,monty::rc_ptr< ::mosek::fusion::Expression > _7989_e2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > sub(monty::rc_ptr< ::mosek::fusion::Expression > _7990_e1,monty::rc_ptr< ::mosek::fusion::Matrix > _7991_m);
      static  monty::rc_ptr< ::mosek::fusion::Expression > sub(double _7992_c,monty::rc_ptr< ::mosek::fusion::Expression > _7993_e2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > sub(monty::rc_ptr< ::mosek::fusion::Expression > _7994_e1,double _7995_c);
      static  monty::rc_ptr< ::mosek::fusion::Expression > sub(std::shared_ptr< monty::ndarray< double,2 > > _7996_a2,monty::rc_ptr< ::mosek::fusion::Expression > _7997_e2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > sub(std::shared_ptr< monty::ndarray< double,1 > > _7998_a1,monty::rc_ptr< ::mosek::fusion::Expression > _7999_e2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > sub(monty::rc_ptr< ::mosek::fusion::Expression > _8000_e1,std::shared_ptr< monty::ndarray< double,2 > > _8001_a2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > sub(monty::rc_ptr< ::mosek::fusion::Expression > _8002_e1,std::shared_ptr< monty::ndarray< double,1 > > _8003_a1);
      static  monty::rc_ptr< ::mosek::fusion::Expression > sub(monty::rc_ptr< ::mosek::fusion::Expression > _8004_e1,monty::rc_ptr< ::mosek::fusion::Expression > _8005_e2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > add(monty::rc_ptr< ::mosek::fusion::NDSparseArray > _8006_n,monty::rc_ptr< ::mosek::fusion::Expression > _8007_e2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > add(monty::rc_ptr< ::mosek::fusion::Expression > _8008_e1,monty::rc_ptr< ::mosek::fusion::NDSparseArray > _8009_n);
      static  monty::rc_ptr< ::mosek::fusion::Expression > add(monty::rc_ptr< ::mosek::fusion::Matrix > _8010_m,monty::rc_ptr< ::mosek::fusion::Expression > _8011_e2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > add(monty::rc_ptr< ::mosek::fusion::Expression > _8012_e1,monty::rc_ptr< ::mosek::fusion::Matrix > _8013_m);
      static  monty::rc_ptr< ::mosek::fusion::Expression > add(double _8014_c,monty::rc_ptr< ::mosek::fusion::Expression > _8015_e2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > add(monty::rc_ptr< ::mosek::fusion::Expression > _8016_e1,double _8017_c);
      static  monty::rc_ptr< ::mosek::fusion::Expression > add(std::shared_ptr< monty::ndarray< double,2 > > _8018_a2,monty::rc_ptr< ::mosek::fusion::Expression > _8019_e2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > add(std::shared_ptr< monty::ndarray< double,1 > > _8020_a1,monty::rc_ptr< ::mosek::fusion::Expression > _8021_e2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > add(monty::rc_ptr< ::mosek::fusion::Expression > _8022_e1,std::shared_ptr< monty::ndarray< double,2 > > _8023_a2);
      static  monty::rc_ptr< ::mosek::fusion::Expression > add(monty::rc_ptr< ::mosek::fusion::Expression > _8024_e1,std::shared_ptr< monty::ndarray< double,1 > > _8025_a1);
      static  monty::rc_ptr< ::mosek::fusion::Expression > add(monty::rc_ptr< ::mosek::fusion::Expression > _8026_e1,monty::rc_ptr< ::mosek::fusion::Expression > _8027_e2);
      virtual /* override */ void eval(monty::rc_ptr< ::mosek::fusion::WorkStack > _8028_rs,monty::rc_ptr< ::mosek::fusion::WorkStack > _8029_ws,monty::rc_ptr< ::mosek::fusion::WorkStack > _8030_xs) ;
      static  void validateData(std::shared_ptr< monty::ndarray< int64_t,1 > > _8048_ptrb,std::shared_ptr< monty::ndarray< int64_t,1 > > _8049_subj,std::shared_ptr< monty::ndarray< double,1 > > _8050_cof,std::shared_ptr< monty::ndarray< double,1 > > _8051_bfix,std::shared_ptr< monty::ndarray< int32_t,1 > > _8052_shape,std::shared_ptr< monty::ndarray< int64_t,1 > > _8053_inst);
      static  monty::rc_ptr< ::mosek::fusion::Model > extractModel(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Variable >,1 > > _8066_v);
    }; // struct Expr;

    struct p_WorkStack
    {
      WorkStack * _pubthis;
      static mosek::fusion::p_WorkStack* _get_impl(mosek::fusion::WorkStack * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_WorkStack * _get_impl(mosek::fusion::WorkStack::t _inst) { return _get_impl(_inst.get()); }
      p_WorkStack(WorkStack * _pubthis);
      virtual ~p_WorkStack() { /* std::cout << "~p_WorkStack" << std::endl;*/ };
      int32_t code_base{};
      int32_t cconst_base{};
      int32_t codeptr_base{};
      int32_t cof_base{};
      int32_t nidxs_base{};
      int32_t sp_base{};
      int32_t shape_base{};
      int32_t ptr_base{};
      bool hassp{};
      int32_t ncodeatom{};
      int32_t nelem{};
      int32_t nnz{};
      int32_t nd{};
      int32_t pf64{};
      int32_t pi64{};
      int32_t pi32{};
      std::shared_ptr< monty::ndarray< double,1 > > f64{};
      std::shared_ptr< monty::ndarray< int64_t,1 > > i64{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > i32{};

      virtual void destroy();

      static WorkStack::t _new_WorkStack();
      void _initialize();
      virtual void validate() ;
      virtual void validate_all() ;
      virtual void validate_from(int32_t _7252_fi32,int32_t _7253_fi64,int32_t _7254_ff64,bool _7255_rec) ;
      virtual std::string formatCurrent() ;
      virtual bool peek_hassp() ;
      virtual int32_t peek_nnz() ;
      virtual int32_t peek_nelem() ;
      virtual int32_t peek_dim(int32_t _7289_i) ;
      virtual int32_t peek_nd() ;
      virtual void alloc_expr(int32_t _7290_nd,int32_t _7291_nelem,int32_t _7292_nnz,bool _7293_hassp) ;
      virtual void alloc_expr(int32_t _7294_nd,int32_t _7295_nelem,int32_t _7296_nnz,bool _7297_hassp,int32_t _7298_ncodeatom) ;
      virtual void pop_expr() ;
      virtual void move_expr(monty::rc_ptr< ::mosek::fusion::WorkStack > _7299_to) ;
      virtual void peek_expr() ;
      virtual void ensure_sparsity() ;
      virtual void clear() ;
      virtual int32_t allocf64(int32_t _7314_n) ;
      virtual int32_t alloci64(int32_t _7316_n) ;
      virtual int32_t alloci32(int32_t _7318_n) ;
      virtual void pushf64(double _7320_v) ;
      virtual void pushi64(int64_t _7321_v) ;
      virtual void pushi32(int32_t _7322_v) ;
      virtual void ensuref64(int32_t _7323_n) ;
      virtual void ensurei64(int32_t _7326_n) ;
      virtual void ensurei32(int32_t _7329_n) ;
      virtual int32_t popf64(int32_t _7332_n) ;
      virtual int32_t popi64(int32_t _7333_n) ;
      virtual int32_t popi32(int32_t _7334_n) ;
      virtual void popf64(int32_t _7335_n,std::shared_ptr< monty::ndarray< double,1 > > _7336_r,int32_t _7337_ofs) ;
      virtual void popi64(int32_t _7338_n,std::shared_ptr< monty::ndarray< int64_t,1 > > _7339_r,int32_t _7340_ofs) ;
      virtual void popi32(int32_t _7341_n,std::shared_ptr< monty::ndarray< int32_t,1 > > _7342_r,int32_t _7343_ofs) ;
      virtual double popf64() ;
      virtual int64_t popi64() ;
      virtual int32_t popi32() ;
      virtual double peekf64() ;
      virtual int64_t peeki64() ;
      virtual int32_t peeki32() ;
      virtual double peekf64(int32_t _7344_i) ;
      virtual int64_t peeki64(int32_t _7345_i) ;
      virtual int32_t peeki32(int32_t _7346_i) ;
    }; // struct WorkStack;

    struct p_SymmetricMatrix
    {
      SymmetricMatrix * _pubthis;
      static mosek::fusion::p_SymmetricMatrix* _get_impl(mosek::fusion::SymmetricMatrix * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_SymmetricMatrix * _get_impl(mosek::fusion::SymmetricMatrix::t _inst) { return _get_impl(_inst.get()); }
      p_SymmetricMatrix(SymmetricMatrix * _pubthis);
      virtual ~p_SymmetricMatrix() { /* std::cout << "~p_SymmetricMatrix" << std::endl;*/ };
      int32_t nnz{};
      double scale{};
      std::shared_ptr< monty::ndarray< double,1 > > vval{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > vsubj{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > vsubi{};
      std::shared_ptr< monty::ndarray< double,1 > > uval{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > usubj{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > usubi{};
      int32_t d1{};
      int32_t d0{};

      virtual void destroy();

      static SymmetricMatrix::t _new_SymmetricMatrix(int32_t _8079_dim0,int32_t _8080_dim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _8081_usubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _8082_usubj,std::shared_ptr< monty::ndarray< double,1 > > _8083_uval,std::shared_ptr< monty::ndarray< int32_t,1 > > _8084_vsubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _8085_vsubj,std::shared_ptr< monty::ndarray< double,1 > > _8086_vval,double _8087_scale);
      void _initialize(int32_t _8079_dim0,int32_t _8080_dim1,std::shared_ptr< monty::ndarray< int32_t,1 > > _8081_usubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _8082_usubj,std::shared_ptr< monty::ndarray< double,1 > > _8083_uval,std::shared_ptr< monty::ndarray< int32_t,1 > > _8084_vsubi,std::shared_ptr< monty::ndarray< int32_t,1 > > _8085_vsubj,std::shared_ptr< monty::ndarray< double,1 > > _8086_vval,double _8087_scale);
      static  monty::rc_ptr< ::mosek::fusion::SymmetricMatrix > rankOne(int32_t _8088_n,std::shared_ptr< monty::ndarray< int32_t,1 > > _8089_sub,std::shared_ptr< monty::ndarray< double,1 > > _8090_v);
      static  monty::rc_ptr< ::mosek::fusion::SymmetricMatrix > rankOne(std::shared_ptr< monty::ndarray< double,1 > > _8098_v);
      static  monty::rc_ptr< ::mosek::fusion::SymmetricMatrix > antiDiag(std::shared_ptr< monty::ndarray< double,1 > > _8106_vals);
      static  monty::rc_ptr< ::mosek::fusion::SymmetricMatrix > diag(std::shared_ptr< monty::ndarray< double,1 > > _8113_vals);
      virtual monty::rc_ptr< ::mosek::fusion::SymmetricMatrix > __mosek_2fusion_2SymmetricMatrix__add(monty::rc_ptr< ::mosek::fusion::SymmetricMatrix > _8119_m) ;
      virtual monty::rc_ptr< ::mosek::fusion::SymmetricMatrix > __mosek_2fusion_2SymmetricMatrix__sub(monty::rc_ptr< ::mosek::fusion::SymmetricMatrix > _8139_m) ;
      virtual monty::rc_ptr< ::mosek::fusion::SymmetricMatrix > __mosek_2fusion_2SymmetricMatrix__mul(double _8140_v) ;
      virtual int32_t getdim() ;
    }; // struct SymmetricMatrix;

    struct p_NDSparseArray
    {
      NDSparseArray * _pubthis;
      static mosek::fusion::p_NDSparseArray* _get_impl(mosek::fusion::NDSparseArray * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_NDSparseArray * _get_impl(mosek::fusion::NDSparseArray::t _inst) { return _get_impl(_inst.get()); }
      p_NDSparseArray(NDSparseArray * _pubthis);
      virtual ~p_NDSparseArray() { /* std::cout << "~p_NDSparseArray" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< double,1 > > cof{};
      std::shared_ptr< monty::ndarray< int64_t,1 > > inst{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > dims{};
      int64_t size{};

      virtual void destroy();

      static NDSparseArray::t _new_NDSparseArray(std::shared_ptr< monty::ndarray< int32_t,1 > > _8141_dims_,std::shared_ptr< monty::ndarray< int32_t,2 > > _8142_sub,std::shared_ptr< monty::ndarray< double,1 > > _8143_cof_);
      void _initialize(std::shared_ptr< monty::ndarray< int32_t,1 > > _8141_dims_,std::shared_ptr< monty::ndarray< int32_t,2 > > _8142_sub,std::shared_ptr< monty::ndarray< double,1 > > _8143_cof_);
      static NDSparseArray::t _new_NDSparseArray(std::shared_ptr< monty::ndarray< int32_t,1 > > _8164_dims_,std::shared_ptr< monty::ndarray< int64_t,1 > > _8165_inst_,std::shared_ptr< monty::ndarray< double,1 > > _8166_cof_);
      void _initialize(std::shared_ptr< monty::ndarray< int32_t,1 > > _8164_dims_,std::shared_ptr< monty::ndarray< int64_t,1 > > _8165_inst_,std::shared_ptr< monty::ndarray< double,1 > > _8166_cof_);
      static NDSparseArray::t _new_NDSparseArray(monty::rc_ptr< ::mosek::fusion::Matrix > _8182_m);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Matrix > _8182_m);
      static  monty::rc_ptr< ::mosek::fusion::NDSparseArray > make(monty::rc_ptr< ::mosek::fusion::Matrix > _8190_m);
      static  monty::rc_ptr< ::mosek::fusion::NDSparseArray > make(std::shared_ptr< monty::ndarray< int32_t,1 > > _8191_dims,std::shared_ptr< monty::ndarray< int64_t,1 > > _8192_inst,std::shared_ptr< monty::ndarray< double,1 > > _8193_cof);
      static  monty::rc_ptr< ::mosek::fusion::NDSparseArray > make(std::shared_ptr< monty::ndarray< int32_t,1 > > _8194_dims,std::shared_ptr< monty::ndarray< int32_t,2 > > _8195_sub,std::shared_ptr< monty::ndarray< double,1 > > _8196_cof);
    }; // struct NDSparseArray;

    struct p_Matrix
    {
      Matrix * _pubthis;
      static mosek::fusion::p_Matrix* _get_impl(mosek::fusion::Matrix * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_Matrix * _get_impl(mosek::fusion::Matrix::t _inst) { return _get_impl(_inst.get()); }
      p_Matrix(Matrix * _pubthis);
      virtual ~p_Matrix() { /* std::cout << "~p_Matrix" << std::endl;*/ };
      int32_t dimj{};
      int32_t dimi{};

      virtual void destroy();

      static Matrix::t _new_Matrix(int32_t _8266_di,int32_t _8267_dj);
      void _initialize(int32_t _8266_di,int32_t _8267_dj);
      virtual /* override */ std::string toString() ;
      virtual void switchDims() ;
      static  monty::rc_ptr< ::mosek::fusion::Matrix > diag(int32_t _8269_num,monty::rc_ptr< ::mosek::fusion::Matrix > _8270_mv);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > diag(std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Matrix >,1 > > _8272_md);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > antidiag(int32_t _8290_n,double _8291_val,int32_t _8292_k);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > antidiag(int32_t _8293_n,double _8294_val);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > diag(int32_t _8295_n,double _8296_val,int32_t _8297_k);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > diag(int32_t _8298_n,double _8299_val);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > antidiag(std::shared_ptr< monty::ndarray< double,1 > > _8300_d,int32_t _8301_k);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > antidiag(std::shared_ptr< monty::ndarray< double,1 > > _8311_d);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > diag(std::shared_ptr< monty::ndarray< double,1 > > _8312_d,int32_t _8313_k);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > diag(std::shared_ptr< monty::ndarray< double,1 > > _8321_d);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > ones(int32_t _8322_n,int32_t _8323_m);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > eye(int32_t _8324_n);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > dense(monty::rc_ptr< ::mosek::fusion::Matrix > _8326_other);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > dense(int32_t _8327_dimi,int32_t _8328_dimj,double _8329_value);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > dense(int32_t _8330_dimi,int32_t _8331_dimj,std::shared_ptr< monty::ndarray< double,1 > > _8332_data);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > dense(std::shared_ptr< monty::ndarray< double,2 > > _8333_data);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > sparse(monty::rc_ptr< ::mosek::fusion::Matrix > _8334_mx);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > sparse(std::shared_ptr< monty::ndarray< std::shared_ptr< monty::ndarray< monty::rc_ptr< ::mosek::fusion::Matrix >,1 > >,1 > > _8338_blocks);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > sparse(std::shared_ptr< monty::ndarray< double,2 > > _8369_data);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > sparse(int32_t _8379_nrow,int32_t _8380_ncol);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > sparse(int32_t _8381_nrow,int32_t _8382_ncol,std::shared_ptr< monty::ndarray< int32_t,1 > > _8383_subi,std::shared_ptr< monty::ndarray< int32_t,1 > > _8384_subj,double _8385_val);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > sparse(std::shared_ptr< monty::ndarray< int32_t,1 > > _8387_subi,std::shared_ptr< monty::ndarray< int32_t,1 > > _8388_subj,double _8389_val);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > sparse(std::shared_ptr< monty::ndarray< int32_t,1 > > _8394_subi,std::shared_ptr< monty::ndarray< int32_t,1 > > _8395_subj,std::shared_ptr< monty::ndarray< double,1 > > _8396_val);
      static  monty::rc_ptr< ::mosek::fusion::Matrix > sparse(int32_t _8401_nrow,int32_t _8402_ncol,std::shared_ptr< monty::ndarray< int32_t,1 > > _8403_subi,std::shared_ptr< monty::ndarray< int32_t,1 > > _8404_subj,std::shared_ptr< monty::ndarray< double,1 > > _8405_val);
      virtual double get(int32_t _8410_i,int32_t _8411_j) { throw monty::AbstractClassError("Call to abstract method"); }
      virtual bool isSparse() { throw monty::AbstractClassError("Call to abstract method"); }
      virtual std::shared_ptr< monty::ndarray< double,1 > > getDataAsArray() { throw monty::AbstractClassError("Call to abstract method"); }
      virtual void getDataAsTriplets(std::shared_ptr< monty::ndarray< int32_t,1 > > _8412_subi,std::shared_ptr< monty::ndarray< int32_t,1 > > _8413_subj,std::shared_ptr< monty::ndarray< double,1 > > _8414_val) { throw monty::AbstractClassError("Call to abstract method"); }
      virtual monty::rc_ptr< ::mosek::fusion::Matrix > __mosek_2fusion_2Matrix__transpose() { throw monty::AbstractClassError("Call to abstract method"); }
      virtual int64_t numNonzeros() { throw monty::AbstractClassError("Call to abstract method"); }
      virtual int32_t numColumns() ;
      virtual int32_t numRows() ;
    }; // struct Matrix;

    struct p_DenseMatrix : public ::mosek::fusion::p_Matrix
    {
      DenseMatrix * _pubthis;
      static mosek::fusion::p_DenseMatrix* _get_impl(mosek::fusion::DenseMatrix * _inst){ return static_cast< mosek::fusion::p_DenseMatrix* >(mosek::fusion::p_Matrix::_get_impl(_inst)); }
      static mosek::fusion::p_DenseMatrix * _get_impl(mosek::fusion::DenseMatrix::t _inst) { return _get_impl(_inst.get()); }
      p_DenseMatrix(DenseMatrix * _pubthis);
      virtual ~p_DenseMatrix() { /* std::cout << "~p_DenseMatrix" << std::endl;*/ };
      int64_t nnz{};
      std::shared_ptr< monty::ndarray< double,1 > > data{};

      virtual void destroy();

      static DenseMatrix::t _new_DenseMatrix(int32_t _8197_dimi_,int32_t _8198_dimj_,std::shared_ptr< monty::ndarray< double,1 > > _8199_cof);
      void _initialize(int32_t _8197_dimi_,int32_t _8198_dimj_,std::shared_ptr< monty::ndarray< double,1 > > _8199_cof);
      static DenseMatrix::t _new_DenseMatrix(monty::rc_ptr< ::mosek::fusion::Matrix > _8200_m_);
      void _initialize(monty::rc_ptr< ::mosek::fusion::Matrix > _8200_m_);
      static DenseMatrix::t _new_DenseMatrix(std::shared_ptr< monty::ndarray< double,2 > > _8205_d);
      void _initialize(std::shared_ptr< monty::ndarray< double,2 > > _8205_d);
      static DenseMatrix::t _new_DenseMatrix(int32_t _8208_dimi_,int32_t _8209_dimj_,double _8210_value_);
      void _initialize(int32_t _8208_dimi_,int32_t _8209_dimj_,double _8210_value_);
      virtual /* override */ std::string toString() ;
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Matrix > __mosek_2fusion_2DenseMatrix__transpose() ;
      virtual monty::rc_ptr< ::mosek::fusion::Matrix > __mosek_2fusion_2Matrix__transpose() { return __mosek_2fusion_2DenseMatrix__transpose(); }
      virtual /* override */ bool isSparse() ;
      virtual /* override */ std::shared_ptr< monty::ndarray< double,1 > > getDataAsArray() ;
      virtual /* override */ void getDataAsTriplets(std::shared_ptr< monty::ndarray< int32_t,1 > > _8223_subi,std::shared_ptr< monty::ndarray< int32_t,1 > > _8224_subj,std::shared_ptr< monty::ndarray< double,1 > > _8225_cof) ;
      virtual /* override */ double get(int32_t _8229_i,int32_t _8230_j) ;
      virtual /* override */ int64_t numNonzeros() ;
    }; // struct DenseMatrix;

    struct p_SparseMatrix : public ::mosek::fusion::p_Matrix
    {
      SparseMatrix * _pubthis;
      static mosek::fusion::p_SparseMatrix* _get_impl(mosek::fusion::SparseMatrix * _inst){ return static_cast< mosek::fusion::p_SparseMatrix* >(mosek::fusion::p_Matrix::_get_impl(_inst)); }
      static mosek::fusion::p_SparseMatrix * _get_impl(mosek::fusion::SparseMatrix::t _inst) { return _get_impl(_inst.get()); }
      p_SparseMatrix(SparseMatrix * _pubthis);
      virtual ~p_SparseMatrix() { /* std::cout << "~p_SparseMatrix" << std::endl;*/ };
      int64_t nnz{};
      std::shared_ptr< monty::ndarray< double,1 > > val{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > subj{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > subi{};

      virtual void destroy();

      static SparseMatrix::t _new_SparseMatrix(int32_t _8231_dimi_,int32_t _8232_dimj_,std::shared_ptr< monty::ndarray< int32_t,1 > > _8233_subi_,std::shared_ptr< monty::ndarray< int32_t,1 > > _8234_subj_,std::shared_ptr< monty::ndarray< double,1 > > _8235_val_,int64_t _8236_nelm);
      void _initialize(int32_t _8231_dimi_,int32_t _8232_dimj_,std::shared_ptr< monty::ndarray< int32_t,1 > > _8233_subi_,std::shared_ptr< monty::ndarray< int32_t,1 > > _8234_subj_,std::shared_ptr< monty::ndarray< double,1 > > _8235_val_,int64_t _8236_nelm);
      static SparseMatrix::t _new_SparseMatrix(int32_t _8242_dimi_,int32_t _8243_dimj_,std::shared_ptr< monty::ndarray< int32_t,1 > > _8244_subi_,std::shared_ptr< monty::ndarray< int32_t,1 > > _8245_subj_,std::shared_ptr< monty::ndarray< double,1 > > _8246_val_);
      void _initialize(int32_t _8242_dimi_,int32_t _8243_dimj_,std::shared_ptr< monty::ndarray< int32_t,1 > > _8244_subi_,std::shared_ptr< monty::ndarray< int32_t,1 > > _8245_subj_,std::shared_ptr< monty::ndarray< double,1 > > _8246_val_);
      virtual std::shared_ptr< monty::ndarray< int64_t,1 > > formPtrb() ;
      virtual /* override */ std::string toString() ;
      virtual /* override */ int64_t numNonzeros() ;
      virtual /* override */ monty::rc_ptr< ::mosek::fusion::Matrix > __mosek_2fusion_2SparseMatrix__transpose() ;
      virtual monty::rc_ptr< ::mosek::fusion::Matrix > __mosek_2fusion_2Matrix__transpose() { return __mosek_2fusion_2SparseMatrix__transpose(); }
      virtual /* override */ bool isSparse() ;
      virtual /* override */ std::shared_ptr< monty::ndarray< double,1 > > getDataAsArray() ;
      virtual /* override */ void getDataAsTriplets(std::shared_ptr< monty::ndarray< int32_t,1 > > _8258_subi_,std::shared_ptr< monty::ndarray< int32_t,1 > > _8259_subj_,std::shared_ptr< monty::ndarray< double,1 > > _8260_cof_) ;
      virtual /* override */ double get(int32_t _8261_i,int32_t _8262_j) ;
    }; // struct SparseMatrix;

    struct p_LinkedBlocks
    {
      LinkedBlocks * _pubthis;
      static mosek::fusion::p_LinkedBlocks* _get_impl(mosek::fusion::LinkedBlocks * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_LinkedBlocks * _get_impl(mosek::fusion::LinkedBlocks::t _inst) { return _get_impl(_inst.get()); }
      p_LinkedBlocks(LinkedBlocks * _pubthis);
      virtual ~p_LinkedBlocks() { /* std::cout << "~p_LinkedBlocks" << std::endl;*/ };
      std::shared_ptr< monty::ndarray< int32_t,1 > > bfirst{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > bsize{};
      monty::rc_ptr< ::mosek::fusion::LinkedInts > blocks{};
      monty::rc_ptr< ::mosek::fusion::LinkedInts > ints{};

      virtual void destroy();

      static LinkedBlocks::t _new_LinkedBlocks();
      void _initialize();
      static LinkedBlocks::t _new_LinkedBlocks(int32_t _8440_n);
      void _initialize(int32_t _8440_n);
      static LinkedBlocks::t _new_LinkedBlocks(monty::rc_ptr< ::mosek::fusion::LinkedBlocks > _8442_other);
      void _initialize(monty::rc_ptr< ::mosek::fusion::LinkedBlocks > _8442_other);
      virtual void free(int32_t _8443_bkey) ;
      virtual int32_t alloc(int32_t _8445_size) ;
      virtual int32_t maxidx(int32_t _8450_bkey) ;
      virtual int32_t numallocated() ;
      virtual void get(int32_t _8451_bkey,std::shared_ptr< monty::ndarray< int32_t,1 > > _8452_target,int32_t _8453_offset) ;
      virtual int32_t numblocks() ;
      virtual int32_t blocksize(int32_t _8454_bkey) ;
      virtual int32_t block_capacity() ;
      virtual int32_t capacity() ;
      virtual bool validate() ;
    }; // struct LinkedBlocks;

    struct p_LinkedInts
    {
      LinkedInts * _pubthis;
      static mosek::fusion::p_LinkedInts* _get_impl(mosek::fusion::LinkedInts * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_LinkedInts * _get_impl(mosek::fusion::LinkedInts::t _inst) { return _get_impl(_inst.get()); }
      p_LinkedInts(LinkedInts * _pubthis);
      virtual ~p_LinkedInts() { /* std::cout << "~p_LinkedInts" << std::endl;*/ };
      int32_t nfree{};
      int32_t last_free{};
      int32_t first_free{};
      int32_t first_used{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > prev{};
      std::shared_ptr< monty::ndarray< int32_t,1 > > next{};

      virtual void destroy();

      static LinkedInts::t _new_LinkedInts(int32_t _8455_cap_);
      void _initialize(int32_t _8455_cap_);
      static LinkedInts::t _new_LinkedInts();
      void _initialize();
      static LinkedInts::t _new_LinkedInts(monty::rc_ptr< ::mosek::fusion::LinkedInts > _8458_other);
      void _initialize(monty::rc_ptr< ::mosek::fusion::LinkedInts > _8458_other);
      virtual void free(int32_t _8459_i,int32_t _8460_num) ;
      virtual int32_t alloc() ;
      virtual int32_t alloc(int32_t _8466_n) ;
      virtual void alloc(int32_t _8467_num,std::shared_ptr< monty::ndarray< int32_t,1 > > _8468_target,int32_t _8469_offset) ;
      virtual void get(int32_t _8472_i,int32_t _8473_num,std::shared_ptr< monty::ndarray< int32_t,1 > > _8474_target,int32_t _8475_offset) ;
      virtual int32_t numallocated() ;
      virtual int32_t maxidx(int32_t _8478_i,int32_t _8479_num) ;
      virtual int32_t allocblock(int32_t _8483_num) ;
      virtual void recap(int32_t _8489_ncap) ;
      virtual int32_t capacity() ;
      virtual bool validate() ;
    }; // struct LinkedInts;

    struct p_Parameters
    {
      Parameters * _pubthis;
      static mosek::fusion::p_Parameters* _get_impl(mosek::fusion::Parameters * _inst){ assert(_inst); assert(_inst->_impl); return _inst->_impl; }
      static mosek::fusion::p_Parameters * _get_impl(mosek::fusion::Parameters::t _inst) { return _get_impl(_inst.get()); }
      p_Parameters(Parameters * _pubthis);
      virtual ~p_Parameters() { /* std::cout << "~p_Parameters" << std::endl;*/ };

      virtual void destroy();

      static  void setParameter(monty::rc_ptr< ::mosek::fusion::Model > _8498_M,const std::string &  _8499_name,double _8500_value);
      static  void setParameter(monty::rc_ptr< ::mosek::fusion::Model > _8608_M,const std::string &  _8609_name,int32_t _8610_value);
      static  void setParameter(monty::rc_ptr< ::mosek::fusion::Model > _8718_M,const std::string &  _8719_name,const std::string &  _8720_value);
      static  int32_t string_to_variabletype_value(const std::string &  _9002_v);
      static  int32_t string_to_value_value(const std::string &  _9003_v);
      static  int32_t string_to_streamtype_value(const std::string &  _9004_v);
      static  int32_t string_to_startpointtype_value(const std::string &  _9005_v);
      static  int32_t string_to_stakey_value(const std::string &  _9006_v);
      static  int32_t string_to_sparam_value(const std::string &  _9007_v);
      static  int32_t string_to_solveform_value(const std::string &  _9008_v);
      static  int32_t string_to_soltype_value(const std::string &  _9009_v);
      static  int32_t string_to_solsta_value(const std::string &  _9010_v);
      static  int32_t string_to_solitem_value(const std::string &  _9011_v);
      static  int32_t string_to_simseltype_value(const std::string &  _9012_v);
      static  int32_t string_to_sensitivitytype_value(const std::string &  _9013_v);
      static  int32_t string_to_scalingmethod_value(const std::string &  _9014_v);
      static  int32_t string_to_scalingtype_value(const std::string &  _9015_v);
      static  int32_t string_to_rescodetype_value(const std::string &  _9016_v);
      static  int32_t string_to_rescode_value(const std::string &  _9017_v);
      static  int32_t string_to_prosta_value(const std::string &  _9018_v);
      static  int32_t string_to_problemtype_value(const std::string &  _9019_v);
      static  int32_t string_to_problemitem_value(const std::string &  _9020_v);
      static  int32_t string_to_parametertype_value(const std::string &  _9021_v);
      static  int32_t string_to_foldingmode_value(const std::string &  _9022_v);
      static  int32_t string_to_presolvemode_value(const std::string &  _9023_v);
      static  int32_t string_to_orderingtype_value(const std::string &  _9024_v);
      static  int32_t string_to_optimizertype_value(const std::string &  _9025_v);
      static  int32_t string_to_onoffkey_value(const std::string &  _9026_v);
      static  int32_t string_to_objsense_value(const std::string &  _9027_v);
      static  int32_t string_to_mpsformat_value(const std::string &  _9028_v);
      static  int32_t string_to_miovarseltype_value(const std::string &  _9029_v);
      static  int32_t string_to_mionodeseltype_value(const std::string &  _9030_v);
      static  int32_t string_to_miomode_value(const std::string &  _9031_v);
      static  int32_t string_to_miocontsoltype_value(const std::string &  _9032_v);
      static  int32_t string_to_miodatapermmethod_value(const std::string &  _9033_v);
      static  int32_t string_to_miqcqoreformmethod_value(const std::string &  _9034_v);
      static  int32_t string_to_branchdir_value(const std::string &  _9035_v);
      static  int32_t string_to_iparam_value(const std::string &  _9036_v);
      static  int32_t string_to_iomode_value(const std::string &  _9037_v);
      static  int32_t string_to_internal_iinf_value(const std::string &  _9038_v);
      static  int32_t string_to_internal_dinf_value(const std::string &  _9039_v);
      static  int32_t string_to_inftype_value(const std::string &  _9040_v);
      static  int32_t string_to_iinfitem_value(const std::string &  _9041_v);
      static  int32_t string_to_internal_liinf_value(const std::string &  _9042_v);
      static  int32_t string_to_liinfitem_value(const std::string &  _9043_v);
      static  int32_t string_to_dparam_value(const std::string &  _9044_v);
      static  int32_t string_to_feature_value(const std::string &  _9045_v);
      static  int32_t string_to_dinfitem_value(const std::string &  _9046_v);
      static  int32_t string_to_solformat_value(const std::string &  _9047_v);
      static  int32_t string_to_dataformat_value(const std::string &  _9048_v);
      static  int32_t string_to_symmattype_value(const std::string &  _9049_v);
      static  int32_t string_to_nametype_value(const std::string &  _9050_v);
      static  int32_t string_to_domaintype_value(const std::string &  _9051_v);
      static  int32_t string_to_conetype_value(const std::string &  _9052_v);
      static  int32_t string_to_compresstype_value(const std::string &  _9053_v);
      static  int32_t string_to_callbackcode_value(const std::string &  _9054_v);
      static  int32_t string_to_intpnthotstart_value(const std::string &  _9055_v);
      static  int32_t string_to_simhotstart_value(const std::string &  _9056_v);
      static  int32_t string_to_simdupvec_value(const std::string &  _9057_v);
      static  int32_t string_to_simreform_value(const std::string &  _9058_v);
      static  int32_t string_to_uplo_value(const std::string &  _9059_v);
      static  int32_t string_to_transpose_value(const std::string &  _9060_v);
      static  int32_t string_to_simdegen_value(const std::string &  _9061_v);
      static  int32_t string_to_simprecision_value(const std::string &  _9062_v);
      static  int32_t string_to_mark_value(const std::string &  _9063_v);
      static  int32_t string_to_boundkey_value(const std::string &  _9064_v);
      static  int32_t string_to_basindtype_value(const std::string &  _9065_v);
      static  int32_t string_to_language_value(const std::string &  _9066_v);
    }; // struct Parameters;

  }
}
namespace mosek
{
  namespace fusion
  {
    namespace Utils
    {
      // mosek.fusion.Utils.IntMap from file 'src/fusion/cxx/IntMap_p.h'
      struct p_IntMap
      {
        IntMap * _pubself;
      
        static p_IntMap * _get_impl(IntMap * _inst) { return _inst->_impl.get(); }
      
        p_IntMap(IntMap * _pubself) : _pubself(_pubself) {}
      
        static IntMap::t _new_IntMap() { return new IntMap(); }
      
        ::std::unordered_map<int64_t,int> m;
      
        bool hasItem (int64_t key) { return m.find(key) != m.end(); }
        int  getItem (int64_t key) { return m.find(key)->second; } // will probably throw something or crash of no such key
        void setItem (int64_t key, int val) { m[key] = val; }
      
        std::shared_ptr<monty::ndarray<int64_t,1>> keys()
        {
          size_t size = m.size();
          auto res = std::shared_ptr<monty::ndarray<int64_t,1>>(new monty::ndarray<int64_t,1>(monty::shape((int)size)));
      
          ptrdiff_t i = 0;
          for (auto it = m.begin(); it != m.end(); ++it)
            (*res)[i++] = it->first;
      
          return res;
        }
      
        std::shared_ptr<monty::ndarray<int,1>> values()
        {
          size_t size = m.size();
          auto res = std::shared_ptr<monty::ndarray<int,1>>(new monty::ndarray<int,1>(monty::shape((int)size)));
      
          ptrdiff_t i = 0;
          for (auto it = m.begin(); it != m.end(); ++it)
            (*res)[i++] = it->second;
      
          return res;
        }
      
        IntMap::t clone();
        IntMap::t __mosek_2fusion_2Utils_2IntMap__clone();
      };
      
      
      
      struct p_StringIntMap
      {
        StringIntMap * _pubself;
      
        static p_StringIntMap * _get_impl(StringIntMap * _inst) { return _inst->_impl.get(); }
      
        p_StringIntMap(StringIntMap * _pubself) : _pubself(_pubself) {}
      
        static StringIntMap::t _new_StringIntMap() { return new StringIntMap(); }
      
        ::std::unordered_map<std::string,int> m;
      
        bool hasItem (const std::string & key) { return m.find(key) != m.end(); }
        int  getItem (const std::string & key) { return m.find(key)->second; } // will probably throw something or crash of no such key
        void setItem (const std::string & key, int val) { m[key] = val; }
      
        std::shared_ptr<monty::ndarray<std::string,1>> keys()
        {
          size_t size = m.size();
          auto res = std::shared_ptr<monty::ndarray<std::string,1>>(new monty::ndarray<std::string,1>(monty::shape((int)size)));
      
          ptrdiff_t i = 0;
          for (auto it = m.begin(); it != m.end(); ++it)
            (*res)[i++] = it->first;
      
          return res;
        }
      
        std::shared_ptr<monty::ndarray<int,1>> values()
        {
          size_t size = m.size();
          auto res = std::shared_ptr<monty::ndarray<int,1>>(new monty::ndarray<int,1>(monty::shape((int)size)));
      
          ptrdiff_t i = 0;
          for (auto it = m.begin(); it != m.end(); ++it)
            (*res)[i++] = it->second;
      
          return res;
        }
      
        StringIntMap::t clone();
        StringIntMap::t __mosek_2fusion_2Utils_2StringIntMap__clone();
      };
      // End of file 'src/fusion/cxx/IntMap_p.h'
      // mosek.fusion.Utils.StringBuffer from file 'src/fusion/cxx/StringBuffer_p.h'
      // namespace mosek::fusion::Utils
      struct p_StringBuffer
      {
        StringBuffer * _pubthis;
        std::stringstream ss;
      
        p_StringBuffer(StringBuffer * _pubthis) : _pubthis(_pubthis) {}
      
        static p_StringBuffer * _get_impl(StringBuffer::t ptr) { return ptr->_impl.get(); }
        static p_StringBuffer * _get_impl(StringBuffer * ptr) { return ptr->_impl.get(); }
      
        static StringBuffer::t _new_StringBuffer() { return new StringBuffer(); }
      
        StringBuffer::t clear ();
      
        StringBuffer::t a (const monty::ndarray<std::string,1> & val);
        StringBuffer::t a (const monty::ndarray<int,1> & val);
        StringBuffer::t a (const monty::ndarray<int64_t,1> & val);
        StringBuffer::t a (const monty::ndarray<double,1> & val);
      
      
        StringBuffer::t a (const int & val);
        StringBuffer::t a (const int64_t & val);
        StringBuffer::t a (const double & val);
        StringBuffer::t a (const bool & val);
        StringBuffer::t a (const std::string & val);
      
        StringBuffer::t lf ();
        StringBuffer::t __mosek_2fusion_2Utils_2StringBuffer__clear ();
      
        StringBuffer::t __mosek_2fusion_2Utils_2StringBuffer__a (const monty::ndarray<std::string,1> & val);
        StringBuffer::t __mosek_2fusion_2Utils_2StringBuffer__a (const monty::ndarray<int,1> & val);
        StringBuffer::t __mosek_2fusion_2Utils_2StringBuffer__a (const monty::ndarray<int64_t,1> & val);
        StringBuffer::t __mosek_2fusion_2Utils_2StringBuffer__a (const monty::ndarray<double,1> & val);
      
        StringBuffer::t __mosek_2fusion_2Utils_2StringBuffer__a (const int & val);
        StringBuffer::t __mosek_2fusion_2Utils_2StringBuffer__a (const int64_t & val);
        StringBuffer::t __mosek_2fusion_2Utils_2StringBuffer__a (const double & val);
        StringBuffer::t __mosek_2fusion_2Utils_2StringBuffer__a (const bool & val);
        StringBuffer::t __mosek_2fusion_2Utils_2StringBuffer__a (const std::string & val);
      
        StringBuffer::t __mosek_2fusion_2Utils_2StringBuffer__lf ();
      
        std::string toString () const;
      };
      // End of file 'src/fusion/cxx/StringBuffer_p.h'
    }
  }
}
#endif
