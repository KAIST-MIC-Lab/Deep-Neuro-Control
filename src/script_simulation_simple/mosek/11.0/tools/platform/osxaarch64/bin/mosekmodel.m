classdef mosekmodel < handle
    properties
        modelName    string

        % objective
        objSense     string 
        c            (:,1) double
        cfix         double
        
        % Variables
        % Lengths must be the same and indicates the number of scalar variables
        numVar        uint32
        varType       (:,1) uint8
        varNames      (:,1) string

        % Constraints
        % Ft'x+g defines constraint matrix rows
        % Keep F in transposed for to allow cheaper row and columns adding
        Ft            (:,:) double
        g             (:,1) double

        % Ft as components. These are set immediately before call to solve
        % because accessing sparse matrix components from c++ is super
        % inefficient!
        F_subi (:,1) uint64
        F_subj (:,1) uint64
        F_cofij (:,1) double

        % For efficiency reasons we over-allocate Ft and g. F_ncols and and F_nrows keeps track of how much we are actually using
        F_nrows       uint64

        % Continuous conic constraints.
        %
        %   Ax+b in K1 x ... x Kn + RHS
        %  
        conIndexes    (:,1) uint64 % Constraint rows from `F` and `g`
        % We overallocate conIndexes, and the actual number of used elements in conIndexes is indicated by conIndexes_n
        conIndexes_n  uint64
        conDomain     (:,1) mosekdomain % Domains for each constraint block
        conDomain_n   uint64
        numBlock      uint64
        conBlockSize  (:,1) uint64 % For each constraint block, this is the number of domain entries in the block
        conBlockNames (:,1) string % 

        % DJCs
        djcIndexes    (:,1) uint64 % Constraint rows from `F` and `g`
        djcClauseNTerm(:,1) uint64 % Number of terms in each clause
        djcNumClauses (:,1) uint64 % Number of clauses in each DJC
        djcDomain     (:,1) mosekdomain % Domains for each constraint block
        djcName       (:,1) string % Array of strings
    
        % Which solution to use as init sol:
        % 0 or undefined: None
        % 1 : Use init_sol_...
        % 2 : Use current solution, if available, in order of priority:
        %   1 : integer
        %   2 : basis
        %   3 : interior
        use_init_sol int32

        % initial solution
        init_sol_x   (:,1) double
        init_sol_sc  (:,1) double

        % current solutions:
        sol_itr_available  logical
        sol_itr_prosta     string
        sol_itr_solsta     string
        sol_itr_x          (:,1) double
        sol_itr_y          (:,1) double
 
        sol_bas_available  logical
        sol_bas_prosta     string
        sol_bas_solsta     string
        sol_bas_x          (:,1) double
        sol_bas_y          (:,1) double

        sol_itg_available  logical
        sol_itg_prosta     string
        sol_itg_solsta     string
        sol_itg_x          (:,1) double

        % latest info items
        info
    end % properties

    methods
        function model = mosekmodel(options)
            % MOSEKMODEL  Create a new model object. All model data can be specified at 
            %    creation time with optional arguments or be added subsequently.
            %
            %    model = MOSEKMODEL(options...)
            %
            %    options:
            %    - name string (default: "model")
            %      Name of the model.
            %    - numVar uint32 (default: 0)
            %      Number of scalar variables in the model.
            %    - objsense string (default: "min")
            %      Objective sense, either "min", "minimize", "max" or "maximize"
            %    - objective (:,1) double
            %      Objective vector, this must match the given numVar.
            %    - objfixterm double     
            %      Objective fixed scalar term.
            %    - intvars (:,1) uint64
            %      List of indexes of variables that are integer.
            %    - F (:,:) double
            %      The constraint matrix.
            %    - g (:,1) double (default: zeros(size(F,1),1))
            %      The constraint left-hand constant vector
            %    - conIndexes (:,1) uint64 (default: [1:size(f,1))]
            %      This is the list of row-indexes of F that
            %      define the actual constraints. This allows is to reuse the
            %      same F row in multiple constraint without repeating the row.
            %    - conDomain (:,1) mosekdomain
            %      List of constraint domains.
            %    - solution_x (:,1) double
            %      Initial x solution.
            %    - solution_y (:,1) double
            %      Initial y solution.
            %
            %    See also: APPENDVARS, APPENDCONS, APPENDROWS, OBJECTIVE
            arguments
                options.name          string       = "model"
                options.objsense      string       = "min"
                options.objective     (:,1) double
                options.objfixterm    double

                options.numvar        uint32       = 0
                options.intvars       (:,1) uint64

                options.F             (:,:) double
                options.g             (:,1) double

                options.conindexes    (:,1) uint64 % optional, default is 1..numrows
                options.domain        (:,1) mosekdomain

                % Initial solution 
                options.solution_x    (:,1) double
                options.solution_y    (:,1) double
                %options.solution_sk   (:,1) int
            end % arguments

            model.modelName = options.name;
            switch (options.objsense)
                case {"min","minimize"}
                    model.objSense  = "min";
                case {"max","maximize"}
                    model.objSense  = "max";
                otherwise
                    err = MException("mosekmodel:DataError","Invalid objective sense");
                    throw(err);
            end

            if isfield(options,"objfixterm")
                model.cfix = options.objfixterm;
            end

            model.sol_itr_available = false;
            model.sol_bas_available = false;
            model.sol_itg_available = false;

            numvar = options.numvar;
            model.numVar = numvar;
            model.varType = zeros(numvar,1,'uint8');
            if isfield(options,"intvars")
                if max(options.intvars) > numvar
                    err = MException("mosekmodel:DataError","Integer variable index out of bounds");
                    throw(err);
                end
                model.varType(options.intvars) = char(1);
            end

            if isfield(options,"objective")
                if size(options.objective,1) > numvar
                    err = MException("mosekmodel:DataError","Objective vector exceeds number of variables");
                    throw(err);
                end
                model.c = options.objective;
            else
                model.c = zeros(numvar,1);
            end

            if isfield(options,"objfixterm")
                model.cfix = options.objfixterm;
            else
                model.cfix = 0;
            end

            nrow = 0;
            if isfield(options,"F")
                nrow = size(options.F,1);
            elseif isfield(options,"g")
                nrow = size(options.g,1);
            end

            if isfield(options,"F")
                if size(options.F,2) > numvar 
                    err = MException("mosekmodel:DataError","Coefficient matrix width exceeds number of variables");
                    throw(err);
                elseif isfield(options,"g") && size(options.F,1) ~= size(options.g,1)
                    err = MException("mosekmodel:DataError","Mismatching F and g sizes");
                    throw(err);
                end
                model.Ft = sparse(options.F');
            else 
                model.Ft = sparse(numvar,nrow);
            end 
            model.F_nrows = nrow;

            if isfield(options,"g")
                if size(options.g,1) ~= nrow
                    err = MException("mosekmodel:DataError","Mismatching F and g sizes");
                    throw(err);
                end
                model.g = options.g;
            else 
                model.g = zeros(nrow,1);
            end

            % Constarint blocks

            if isfield(options,"domain") || isfield(options,"conindexes")
                if ~isfield(options,"domain")
                    err = MException("mosekmodel:DataError","Missing conDomain argument");
                    throw(err);
                end
                
                nblocks = size(options.domain,1);
                nrows = 0;
                for i = 1:nblocks
                    nrows = nrows + options.domain(i,1).numrow();
                end

                if isfield(options,"conindexes") 
                    if max(options.conindexes) > size(model.Ft,2)
                        err = MException("mosekmodel:DataError","Argument 'conindexes' contains entry that is out of bounds");
                        throw(err);
                    elseif size(options.conindexes,1) ~= nrows
                        err = MException("mosekmodel:DataError","Length of 'conindexes' does not match the domains");
                        throw(err);
                    end

                    model.conIndexes = options.conindexes;
                else
                    model.conIndexes = uint64([1:nrows]);
                end

                model.conDomain    = options.domain;
            else
                model.conIndexes   = uint64.empty(0,1);
                model.conDomain    = mosekdomain.empty(0,1);
            end

            model.numBlock = 0;
            modek.conBlockSize = zeros(256,1,'uint64');
            modek.conBlockNames = [];

            if isfield(options,"solution_x")
                model.use_init_sol = 1;
                model.init_sol_x  = options.solution_x;
            end
            if isfield(options,"solution_y")
                model.use_init_sol = 1;
                model.init_sol_y  = options.solution_y;
            end
            %if isfield(options,"solution_sk")
            %    model.init_sol_sk = options.solution_sk;
            %end

            model.conIndexes_n = size(model.conIndexes,1);
            model.conDomain_n  = size(model.conDomain,1);

            model.djcIndexes     = 0;
            model.djcClauseNTerm = [];
            model.djcNumClauses  = [];
            %model.djcDomain      = 
            model.djcName        = [];
        end % function MosekMatlab


        function varname(model,indexes,names)
            % VARNAME  Set variable names
            %
            %    VARNAME(model,indexes,names)
            %
            %    Set the names of a subset of variables.
            %
            %    options:
            %    - indexes (:,1) uint64
            %      List of variables indexes to name
            %    - names (:,1) string
            %      List of names.
            arguments
                model mosekmodel
                indexes (:,1) uint64
                names (:,1) string
            end
            
            if max(indexes) > model.numVar || min(indexes) == 0
                err = MException("mosekmodel:DataError","Invalid variable indexes");
                throw(err);
            elseif size(names,1) ~= size(indexes,1)
                err = MException("mosekmodel:DataError","Mismatching lengths of indexes and names");
                throw(err);
            end

            model.varNames(indexes,1) = names;
        end

        function result = appendvars(model,num,options)
            % APPENDVARS  Add new scalar variables to the model.
            %
            %    vars = APPENDVARS(model,num,options...)
            %
            %    Append num scalar variables to the model, returning the indexes of the 
            %    added variables.
            % 
            %    options:
            %    - intVars (:,1) uint64
            %      List if indexes of the new variables that are integer.
            %
            %    See also: APPENDCONS, APPENDROWS, OBJECTIVE
            arguments
                model
                num      uint32 = 0
                options.intvars  (:,1) uint64
            end

            oldnumvar = model.numVar;
            if num > 0
                model.numVar = model.numVar + num;
                model.c = [model.c; zeros(num, 1)];
            end
            addnumvar = num;

            if isfield(options,"intVars")
                if max(options.intvars) > num
                    err = MException("mosekmodel:DataError","Number of constraint row right-hand-sides does not match the number of rows");
                    throw(err);
                end

                newVarType = zeros(num,1,'uint8');
                newVarType(options.intvars) = uint8(1);
                if addnumvar > 0 
                    model.varType = [ model.varType; newVarType ];
                end
            elseif addnumvar > 0
                model.varType = [ model.varType; zeros(num,1,'uint8') ];
            end
            
            if addnumvar > 0
                uint32([ oldnumvar+1 : oldnumvar+num ]);
            else
                uint32([]);
            end
        end % appendvars

        function result = appendrows(model, options)
            % APPENDROWS  Append unused matrix rows to the model. These can subsequently be 
            %    used by APPENDCONS.
            %
            %    rows = APPENDROWS(model,options...)
            
            %    Add new matrix rows, returning the indexes of the added rows.
            %
            %    options:
            %    - F (:,:) double
            %      Matrix block to append
            %    - g (:,1) double
            %      Left-hand side constant term
            arguments
                model
                options.F             (:,:) double
                options.g             (:,1) double
            end

            oldnumrow = model.F_nrows;
            addnrow = 0;
            if isfield(options,"F") && isfield(options,"g")
                if size(options.F,1) ~= size(options.g,1)
                    err = MException("mosekmodel:DataError","Number of F rows do not match length of g");
                    throw(err);
                end
                addnrow = size(options.F,1);
                if size(options.F,2) == model.numVar
                    addFt = options.F';
                else
                    addFt = sparse([],[],[],model.numVar,addnrow,nnz(options.F));
                    addFt(1:size(options.F,2),:) = options.F';
                end
                addg  = options.g;
            elseif isfield(options,"F")
                addnrow = size(options.F,1);
                if size(options.F,2) == model.numVar
                    addFt = options.F';
                else
                    addFt = sparse([],[],[],model.numVar,addnrow,nnz(options.F));
                    addFt(1:size(options.F,2),:) = options.F';
                end
                addg  = zeros(addnrow,1);
            elseif isfield(options,"g")
                addnrow = size(options.g,1);
                addFt = sparse(model.numVar,addnrow);
                addg  = options.g;
            else
                error("F or g must be supplied")
            end

            nrow = model.F_nrows;
            curnnz = nnz(model.Ft);
            capnnz = nzmax(model.Ft);
            rowcap = size(model.g,1); % == size(model.Ft,2)

            addnnz = nnz(addFt);

            % Inrease Ft and g capacity
            if nrow+addnrow > rowcap && curnnz+addnnz > capnnz
                newnnzcap = max(capnnz*2,curnnz+addnnz);
                newrowcap = max(rowcap*2,nrow+addnrow);

                model.g(rowcap+1:newrowcap,1) = zeros(newrowcap-rowcap,1);
                Ft = sparse([],[],[],model.numVar,newrowcap,newnnzcap);
                Ft(:,1:nrow) = model.Ft(:,1:nrow);
                model.Ft = Ft;
            elseif nrow+addnrow > rowcap 
                newrowcap = max(rowcap*2,nrow+addnrow);
                model.Ft(:,rowcap+1:newrowcap) = sparse(model.numVar,newrowcap-rowcap);
                model.g(rowcap+1:newrowcap,1) = zeros(newrowcap-rowcap,1);
            elseif curnnz+addnnz > capnnz
                newnnzcap = max(capnnz*2,curnnz+addnnz);
                Ft = sparse([],[],[],model.numVar,rowcap,newnnzcap);
                Ft(:,1:nrow) = model.Ft(:,1:nrow);
                model.Ft = Ft;
            end

            if isfield(options,"F")
                model.Ft(:,nrow+1:nrow+addnrow) = addFt;
            end

            if isfield(options,"g")
                model.g(nrow+1:nrow+addnrow,1) = addg;
            end
            
            model.F_nrows = model.F_nrows + addnrow;

            if addnrow > 0 
                result = uint64([oldnumrow+1 : oldnumrow+addnrow]');
            else 
                result = uint64.empty(0,1);
            end
        end % appendrows

        function result = appendcons(model, options)
            % APPENDCONS  Append a new constraint block to the model.
            %
            %    cons = APPENDCONS(model,options...)
            %
            %    Append a block of constraints, returning the indexes of the new constraints.
            %
            %    options:
            %    - domain  (:,1) mosekdomain 
            %      Domain or list of domains
            %    - F       (:,:) double
            %      Constraint matrix. If argument x is given, the number of columns in F 
            %      must be equal the number of columns in x.
            %    - g       (:,1) double
            %      Constraint left-hand side constant term.
            %    - x       (:,1) uint32
            %      Optional. Only valid if F is specified. Specifies the variable vector 
            %      for the constraints.
            %    - indexes (:,1) uint64
            %      Refers to constraint indexes already in model. Only allowed if neither 
            %      F or g are specified.
            %    - name string
            %      Name of the constraint block
            arguments
                model
                options.domain  (:,1) mosekdomain % domain or list of domains
                options.F       (:,:) double % 
                options.g       (:,1) double % 
                options.x       (:,1) uint32 % x vector
                options.indexes (:,1) uint64 % refer to constraint indexes already in model, mutually exclusive with F,g
                options.name    string
            end

            oldnumcon = size(model.conIndexes,1);
        
            if ~isfield(options,"domain")
                err = MException("mosekmodel:DataError","Missing 'domain' argument");
                throw(err);
            else
                nblocks = size(options.domain,1);
                conDomain    = [ model.conDomain; options.domain ];
                nrows = 0;
                for i = 1:nblocks
                    nrows = nrows + options.domain(i).numrow();
                end
            end

            
            if isfield(options,"F") || isfield(options,"g") 
                if isfield(options,"indexes")
                    err = MException("mosekmodel:DataError","Arguments 'indexes' and ('F','g') are mutially exclusive");
                    throw(err);
                end

                if isfield(options,"F") && isfield(options,"g")
                    indexes = appendrows(model,F=options.F,g=options.g);
                elseif isfield(options,"F")
                    indexes = appendrows(model,F=options.F);
                else
                    indexes = appendrows(model,g=options.g);
                end
            elseif isfield(options,"indexes")
                if max(options.indexes) > size(model.Ft,2)
                    err = MException("mosekmodel:DataError","Entry in 'indexes' is out of bounds");
                    throw(err);
                end
                indexes = options.indexes;
            else
                err = MException("mosekmodel:DataError","Entry in 'indexes' or 'F' and 'g' must be specified");
                throw(err);
            end

            if size(indexes,1) ~= nrows
                err = MException("mosekmodel:DataError","Mismatching sizes of expression and domain");
                throw(err);
            end


            n = size(indexes,1);
            idxcap = size(model.conIndexes,1);
            if model.conIndexes_n+n > idxcap
                newnum = max(idxcap*2,model.conIndexes_n+n);
                model.conIndexes(idxcap+1:newnum,1) = zeros(newnum-idxcap,1);
            end           
            model.conIndexes(model.conIndexes_n+1:model.conIndexes_n+n,1) = indexes;
            model.conIndexes_n = model.conIndexes_n + n;

            n = size(options.domain,1);
            domcap = size(model.conDomain,1);
            if model.conDomain_n+n > domcap
                newnum = max(domcap*2,model.conDomain_n+n);
                model.conDomain(domcap+1:newnum,1) = repelem(mosekdomain(),newnum-domcap,1); 
            end
            model.conDomain(model.conDomain_n+1:model.conDomain_n+n) = options.domain;
            model.conDomain_n = model.conDomain_n + n;

            model.numBlock = model.numBlock + 1;
            model.conBlockSize(model.numBlock,1) = n;
            if isfield(options,"name")
                model.conBlockNames(model.numBlock,1) = options.name;
            end

            n = size(model.conIndexes,1);

            if nrows > 0 
                result = uint64([oldnumcon+1:oldnumcon+nrows]);
            else 
                result = uint64.empty(0,1);
            end
        end % appendcons

        function result = clause(model,options)
            % CLAUSE  Creates a conic term for a disjunctive constraint.
            %
            %    dom = TERM(F,g,dom)
            %    
            %    Create a conic block representing (Fx+g in dom).
            %    Optional arguments:
            %    - F (:,:) double
            %      The linear coefficient matrix
            %    - g (:,1) double
            %      The left-hand side constant term
            %    - dom (:,1) mosekdomain
            %      An array of domains defining the domain of the term.
            %    - idxs (:,1) uint64
            %      Array of row indexes referring to F rows. This argument is mutually 
            %      exclusive with F and g.
            %
            %    The number of rows in F and g must be the same, and the total size
            %    of dom must be the same as the number of rows of F and g.
            arguments
                model mosekmodel
                options.domain (:,1) mosekdomain
                options.F      (:,:) double
                options.g      (:,1) double
                options.indexes(:,1) uint64
            end


            if ~ isfield(options,"domain") 
                error("Domain argument required")
            else 
                dom = options.domain;
            end

            if size(dom,1) == 1 
                totdomsize = dom.n * dom.domdim;
            else
                tottomsize = sum([dom.n] .* [dom.domdim]);
            end

            if isfield(options,"indexes")
                if isfield(options,"F") || isfield(options,"g")
                    error("Argument indexes cannot be specified with F or g");
                end
                
                idxs = options.indexes;
                if min(idxs) < 1 || max(idxs) > size(model.F_nrows)
                    error("Idxs element out of bounds");
                elseif totdomsize ~= size(idxs,1) 
                    error("Length of idxs does not match domain");
                end
            elseif isfield(options,"F") || isfield(options,"g")
                if isfield(options,"F") && isfield(options,"g")
                    [am,an] = size(options.F);
                    [bm] = size(options.g,1)
                    if bm ~= am
                        error('Mismatching F and g in term');
                    end

                    if an == model.numVar 
                        idxs = model.appendrows(F=options.F,g=options.g);
                    elseif an < model.numVar
                        idxs = model.appendrows(F=[options.F sparse(am,model.numVar-an)],g=g);
                    else
                        error('a matrix is too wide');
                    end
                elseif isfield(options,"F")
                    [am,an] = size(options.F);
                    if an == model.numVar 
                        idxs = model.appendrows(F=options.F);
                    elseif an < model.numVar
                        idxs = model.appendrows(F=[options.F sparse(am,model.numVar-an)]);
                    else
                        error('a matrix is too wide');
                    end
                else
                    [bm] = size(options.g,1)
                    idxs = model.appendrows(g=options.g);
                end
            else
                error("Either indexes or F and g must be specified")
            end

            result = mosekterm(idxs,dom);
        end

        function appenddisjunction(model,terms,options)
            % DISJUNCTION  Append one or more disjunctive constraints
            % 
            %    OBJECTIVE(model,varargin)
            %
            %    Arguments:
            %    - terms (:,1) mosekterm
            %      Disjunctive constraint terms.
            arguments
                model mosekmodel
                terms (:,1) mosekterm
                options.name string
            end
      
            idxs = vertcat(terms.rowidx);
            dom  = vertcat(terms.dom);
            domnum = vertcat(terms.numdom);

            curnidx = size(model.djcIndexes,1);
            curndom = size(model.djcDomain,1);
            curclnterm = size(model.djcClauseNTerm,1);

            model.djcIndexes(curnidx+1:curnidx+size(idxs,1))              = idxs;
            model.djcDomain(curndom+1:curndom+size(dom,1))                = dom;
            model.djcNumClauses(size(model.djcNumClauses,1)+1,1)          = size(domnum,1);
            model.djcClauseNTerm(curclnterm+1:curclnterm+size(domnum,1),1) = domnum;
            if isfield(options,"name")
                model.djcName(size(model.djcName,1)+1,1) = options.name;
            end
        end


        function objective(model, objsense, c, options)
            % OBJECTIVE  Set the objective function
            % 
            %    OBJECTIVE(model,objsense,objective,options)
            %
            %    Arguments:
            %    - objsense string
            %      The objective sense, "min", "minimize", "max", "maximize"
            %    - c (:,1) double
            %      The objective vector. Note that this must match the number of variables 
            %      in the model.
            %    Options:
            %    - objfixterm double
            %      Objective fixed term
            arguments
                model mosekmodel
                objsense string
                c (:,1) double
                options.objfixterm double
            end
            if size(c,1) > model.numVar
                err = MException("mosekmodel:DataError","Objective vector exceeds number of variables in the model");
                throw(err);
            end
            switch objsense
                case {"min","minimize"}
                    model.objSense = "min";
                case {"max","maximize"}
                    model.objSense = "max";
                otherwise
                    err = MException("mosekmodel:DataError","Objective sense is invalid");
                    throw(err);
            end

            if isfield(options,"objfixterm")
                model.cfix = options.objfixterm;
            end

            model.c = c;
        end

        function write(model,filename,options)
            % WRITE  Write problem to a file
            %
            %    WRITE(model,filename)
            %
            %    - filename string
            %      File name to write the problem to. The format is determined 
            %      by the file name extension.
            arguments
                model mosekmodel
                filename string
                options.param  (1,:) string
            end

            if isfield(options,"param")
                mosekenv("write",model,filename,options.param);
            else
                mosekenv("write",model,filename);
            end
        end

        function solve(model,options)
            % SOLVE  Solve the model. 
            %
            %    SOLVE(model,options...)
            % 
            %    options:
            %    - licfile string
            %      Override default license file.
            %    - liccode (1,:) char
            %      Provide a license key as a string
            %    - write_to_file string
            %      Write the problem to this file before solving. The extension specified 
            %      the file format.
            %    - logfile string
            %      Write log to this file.
            %    - nosolve logical
            %      Format and load problem, but do not solve.
            %    - quiet logical
            %      Do not print log.
            %    - optserver string
            %      Optimization server host name. 
            %    - analyze logical
            %      Enable the problem analyzer. This will run after optimization unless 
            %      nosolve=true.
            %    - param  (1,:) string            
            %      Array of parameters in the form [name, value, name, value,...]. Names
            %      are normal MOSEK parameter names, values are parameter values as strings.
            %
            %    See also: HASSOLUTION, GETSOLUTION
            arguments
                model           mosekmodel
                options.quiet   logical
                options.nosolve logical
                options.logfile string
                options.licfile string
                options.liccode (:,1) int32                
                options.write_to_file string
                options.analyze logical
                options.optserver string
                options.param  (1,:) string
            end

            if isfield(options,"licfile")
                licfile = options.licfile;
            else 
                licfile = "";
            end

            if isfield(options,"write_to_file")
                write_to_file = options.write_to_file;
            else 
                write_to_file = "";
            end
            
            if isfield(options,"logfile")
                logfile = options.logfile;
            else 
                logfile = "";
            end
            
            if isfield(options,"liccode")
                liccode = options.liccode;
            else 
                liccode = [];
            end

            if isfield(options,"quiet")
                quiet = options.quiet;
            else
                quiet = false;
            end
            
            if isfield(options,"analyze")
                analyze = options.analyze;
            else
                analyze = false;
            end

            if isfield(options,"optserver")
                optserver = options.optserver;
            else 
                optserver = "";
            end

            nosolve = isfield(options,"nosolve") && options.nosolve;
            
            [F_subi,F_subj,F_cofij] = find(model.Ft');
            model.F_subj = F_subj-1;
            model.F_subi = F_subi-1;
            model.F_cofij = F_cofij;

            if isfield(options,"param")
                [itr_sol,bas_sol,itg_sol,info] = mosekenv("solve",model,nosolve,write_to_file,licfile,liccode,logfile,quiet,optserver,analyze,options.param); 
            else 
                [itr_sol,bas_sol,itg_sol,info] = mosekenv("solve",model,nosolve,write_to_file,licfile,liccode,logfile,quiet,optserver,analyze);
            end
            
           
            if ~ nosolve 
                model.sol_itr_available = false;
                model.sol_bas_available = false;
                model.sol_itg_available = false;
                model.info = info;

                if size(itr_sol,1) == 5
                    model.use_init_sol = 2;
                    model.sol_itr_available = true;
                    itr_prosta = itr_sol{1};
                    itr_solsta = itr_sol{2};
                    switch itr_prosta
                        case 0
                            model.sol_itr_prosta = "UNKNOWN";
                        case 1
                            model.sol_itr_prosta = "PRIM_AND_DUAL_FEAS";
                        case 2
                            model.sol_itr_prosta = "PRIM_FEAS";
                        case 3
                            model.sol_itr_prosta = "DUAL_FEAS";
                        case 4 
                            model.sol_itr_prosta = "PRIM_INFEAS";
                        case 5 
                            model.sol_itr_prosta = "DUAL_INFEAS";
                        case 6 
                            model.sol_itr_prosta = "PRIM_AND_DUAL_INFEAS";
                        case 7 
                            model.sol_itr_prosta = "ILL_POSED";
                        case 8 
                            model.sol_itr_prosta = "PRIM_INFEAS_OR_UNBOUNDED";
                        otherwise 
                            model.sol_itr_available = false;
                    end

                    switch itr_solsta
                        case 0 
                            model.sol_itr_solsta = "UNKNOWN";
                        case 1 
                            model.sol_itr_solsta = "OPTIMAL";
                        case 2 
                            model.sol_itr_solsta = "PRIM_FEAS";
                        case 3 
                            model.sol_itr_solsta = "DUAL_FEAS";
                        case 4 
                            model.sol_itr_solsta = "PRIM_AND_DUAL_FEAS";
                        case 5 
                            model.sol_itr_solsta = "PRIM_INFEAS_CER";
                        case 6 
                            model.sol_itr_solsta = "DUAL_INFEAS_CER";
                        case 7 
                            model.sol_itr_solsta = "PRIM_ILLPOSED_CER";
                        case 8 
                            model.sol_itr_solsta = "DUAL_ILLPOSED_CER";
                        case 9 
                            model.sol_itr_solsta = "INTEGER_OPTIMAL";
                        otherwise 
                            model.sol_itr_available = false;
                    end

                    model.sol_itr_x  = itr_sol{3};
                    model.sol_itr_y  = itr_sol{4};
                    %model.sol_itr_sk = (itr_sol{5} ~= 0);
                end
                
                if size(bas_sol,1) == 5 
                    model.sol_bas_available = true;
                    bas_prosta = bas_sol{1};
                    bas_solsta = bas_sol{2};
                    model.use_init_sol = 2;
                    switch bas_prosta 
                        case 0
                            model.sol_bas_prosta = "UNKNOWN";
                        case 1
                            model.sol_bas_prosta = "PRIM_AND_DUAL_FEAS";
                        case 2
                            model.sol_bas_prosta = "PRIM_FEAS";
                        case 3
                            model.sol_bas_prosta = "DUAL_FEAS";
                        case 4 
                            model.sol_bas_prosta = "PRIM_INFEAS";
                        case 5 
                            model.sol_bas_prosta = "DUAL_INFEAS";
                        case 6 
                            model.sol_bas_prosta = "PRIM_AND_DUAL_INFEAS";
                        case 7 
                            model.sol_bas_prosta = "ILL_POSED";
                        case 8 
                            model.sol_bas_prosta = "PRIM_INFEAS_OR_UNBOUNDED";
                        otherwise 
                            model.sol_bas_available = false;
                    end

                    switch bas_solsta
                        case 0 
                            model.sol_bas_solsta = "UNKNOWN";
                        case 1 
                            model.sol_bas_solsta = "OPTIMAL";
                        case 2 
                            model.sol_bas_solsta = "PRIM_FEAS";
                        case 3 
                            model.sol_bas_solsta = "DUAL_FEAS";
                        case 4 
                            model.sol_bas_solsta = "PRIM_AND_DUAL_FEAS";
                        case 5 
                            model.sol_bas_solsta = "PRIM_INFEAS_CER";
                        case 6 
                            model.sol_bas_solsta = "DUAL_INFEAS_CER";
                        case 7 
                            model.sol_bas_solsta = "PRIM_ILLPOSED_CER";
                        case 8 
                            model.sol_bas_solsta = "DUAL_ILLPOSED_CER";
                        case 9 
                            model.sol_bas_solsta = "INTEGER_OPTIMAL";
                        otherwise 
                            model.sol_bas_available = false;
                    end

                    model.sol_bas_x  = bas_sol{3};
                    model.sol_bas_y  = bas_sol{4};
                    %model.sol_bas_sk = bas_sol{5} ~= 0;
                end
                
                if size(itg_sol,1) == 3
                    itg_prosta = itg_sol{1};
                    itg_solsta = itg_sol{2};
                    model.sol_itg_available = true;
                    model.use_init_sol = 2;
                    switch itg_prosta 
                        case 0
                            model.sol_itg_prosta = "UNKNOWN";
                        case 1
                            model.sol_itg_prosta = "PRIM_AND_DUAL_FEAS";
                        case 2
                            model.sol_itg_prosta = "PRIM_FEAS";
                        case 3
                            model.sol_itg_prosta = "DUAL_FEAS";
                        case 4 
                            model.sol_itg_prosta = "PRIM_INFEAS";
                        case 5 
                            model.sol_itg_prosta = "DUAL_INFEAS";
                        case 6 
                            model.sol_itg_prosta = "PRIM_AND_DUAL_INFEAS";
                        case 7 
                            model.sol_itg_prosta = "ILL_POSED";
                        case 8 
                            model.sol_itg_prosta = "PRIM_INFEAS_OR_UNBOUNDED";
                        otherwise 
                            model.sol_itg_available = false;
                    end

                    switch itg_solsta
                        case 0 
                            model.sol_itg_solsta = "UNKNOWN";
                        case 1 
                            model.sol_itg_solsta = "OPTIMAL";
                        case 2 
                            model.sol_itg_solsta = "PRIM_FEAS";
                        case 3 
                            model.sol_itg_solsta = "DUAL_FEAS";
                        case 4 
                            model.sol_itg_solsta = "PRIM_AND_DUAL_FEAS";
                        case 5 
                            model.sol_itg_solsta = "PRIM_INFEAS_CER";
                        case 6 
                            model.sol_itg_solsta = "DUAL_INFEAS_CER";
                        case 7 
                            model.sol_itg_solsta = "PRIM_ILLPOSED_CER";
                        case 8 
                            model.sol_itg_solsta = "DUAL_ILLPOSED_CER";
                        case 9 
                            model.sol_itg_solsta = "INTEGER_OPTIMAL";
                        otherwise 
                            model.sol_itg_available = false;
                    end

                    model.sol_itg_x = itg_sol{3};
                end
            end
        end
 
        function [result,prosta,solsta] = hassolution(model,which)
            % HASSOLUTION  Check is a specific solution is available and returns it status.
            % 
            %    exists = HASSOLUTION(model,which)
            %    
            %    Only check if the solution is available.
            %
            %    [exists,prosta,solsta] = hassolution(model,which)
            %             
            %    Checks of the solution is available and returns problem and solutions
            %    status if it is.
            %
            %    The which argument is either "interior", "basic" or "integer".
            %
            %    See also: GETSOLUTION, SOLVE
            switch which 
                case {"itr","interior"}
                    result = model.sol_itr_available;
                    prosta = model.sol_itr_prosta;
                    solsta = model.sol_itr_solsta;
                case {"bas","basic","basis"}
                    result = model.sol_bas_available;
                    prosta = model.sol_bas_prosta;
                    solsta = model.sol_bas_solsta;
                case {"itg","integer"}
                    result = model.sol_itg_available;
                    prosta = model.sol_itg_prosta;
                    solsta = model.sol_itg_solsta;
                otherwise
                    err = MException("mosekmodel:KeyError","Invalid solution selector key");
                    throw(err); 
            end
        end

        function setsolution(model,part,value)
            % SETSOLUTION  Set an initial solution. 
            %
            %   SETSOLUTION(model,part,value)
            % 
            %   Sets an initial solution that will be passed to the solver. Note 
            %   that the solution vector lengths are only checked when solve is
            %   called, so they do not need to match the current number of
            %   variables or constraints, only the final numbers.
            %
            %   Arguments:
            %   - part string
            %     Indicates which part of the solution is being provided. One of 
            %     the values "x" (primal solution values), "y" (dual solution 
            %     values) or "sk" (constraint status).
            %   - value (:,1)
            %     Vector of solution values. For part="x" and part="y" this is 
            %     a vector of solution values, for part="sk" this is a vector of 
            %     logicals indicating the the individual constraints are binding.
            arguments
                model mosekmodel
                part string
                value (:,1) 
            end

            model.use_init_sol = 1;
            switch part 
                case "x"
                    model.init_sol_x = value;
                case "y"
                    model.init_sol_y = value;
                %case "sk"
                %    model.init_sol_sk = value
                otherwise
                    err = MException("mosekmodel:KeyError","Invalid solution selector key");
                    throw(err); 
            end
        
        end

        function [result,prosta,solsta] = getsolution(model,which,part) 
            % GETSOLUTION  Fetch a solution.
            %
            %   [result,prosta,solsta] = GETSOLUTION(model,which,part)
            % 
            %   Fetch the solution values, problem status and solution status.
            %   Not that if the solution is a certificate, this will be
            %   indicated by the solsta value.
            %   
            %   Arguments:
            %   - which string (default: "any")
            %     One of the values "any", "interior", "basic" or "integer"
            %   - part string (default : "x")
            %     One of the values 
            %     - "x" indicating primal variable values, and 
            %     - "y" indicating dual constraint values
            % 
            %   See also: HASSOLUTION, SOLVE
            arguments
                model mosekmodel
                which string = "any"
                part  string = "x"
            end
            if which == "any"
                if model.sol_itg_available
                    which = "itg";
                elseif model.sol_bas_available
                    which = "bas";
                else
                    which = "itr";
                end
            end

            switch which 
                case {"itr","interior"}
                    if model.sol_itr_available
                        switch part
                            case "x"
                                result = model.sol_itr_x;
                            case "y"
                                result = model.sol_itr_y;
                            otherwise
                               err = MException("mosekmodel:KeyError","Invalid solution part selector key");
                               throw(err); 
                        end
                        prosta = model.sol_itr_prosta;
                        solsta = model.sol_itr_solsta;
                    else
                        err = MException("mosekmodel:SolutionError","The requested solution is not available");
                        throw(err); 
                    end
                case {"bas","basic","basis"}
                    if model.sol_bas_available
                        switch part
                            case "x"
                                result = model.sol_bas_x;
                            case "y"
                                result = model.sol_bas_y;
                            otherwise
                               err = MException("mosekmodel:KeyError","Invalid solution part selector key");
                               throw(err); 
                        end
                        prosta = model.sol_bas_prosta;
                        solsta = model.sol_bas_solsta;
                    else
                        err = MException("mosekmodel:SolutionError","The requested solution is not available");
                        throw(err); 
                    end
                case {"itg","integer"}
                    if model.sol_itg_available
                        if part == "x"
                            result = model.sol_itg_x;
                            prosta = model.sol_itg_prosta;
                            solsta = model.sol_itg_solsta;
                        else
                            err = MException("mosekmodel:KeyError","Invalid solution part selector key");
                            throw(err); 
                        end
                    else
                        err = MException("mosekmodel:SolutionError","The requested solution is not available");
                        throw(err); 
                    end
                otherwise
                    err = MException("mosekmodel:KeyError","Invalid solution selector key");
                    throw(err); 
            end
        end


        function result = struct(model) 
            arguments 
                model mosekmodel
            end
                
            result.modelName = convertStringsToChars(model.modelName);
            result.objSense = convertStringsToChars(model.objSense); 
            result.c = model.c;                                      % (:,1) double
            result.cfix = model.cfix;                                % double
            result.numVar = model.numVar;                            % uint32
            result.varType = model.varType;                          % (:,1) uint8
            result.varNames = convertStringsToChars(model.varNames); % (:,1) string
            result.Ft = model.Ft;                                    % (:,:) double
            result.g = model.g;                                      % (:,1) double
            result.F_subi = model.F_subi;                            % (:,1) uint64
            result.F_subj = model.F_subj;                            % (:,1) uint64
            result.F_cofij = model.F_cofij;                          % (:,1) double
            result.F_nrows = model.F_nrows;                         % uint64
            result.conIndexes = model.conIndexes;                   % (:,1) uint64 % Constraint rows from `F` and `g`
            result.conIndexes_n = model.conIndexes_n;               % uint64

            result.conDomain = model.conDomain.struct();                     % (:,1) mosekdomain % Domains for each constraint block
            result.conDomain_n = model.conDomain_n;                 % uint64
            result.numBlock = model.numBlock;                       % uint64
            result.conBlockSize = model.conBlockSize;               % (:,1) uint64 % For each constraint block, this is the number of domain entries in the block
            result.conBlockNames = convertStringsToChars(model.conBlockNames);              % (:,1) string % 

            result.djcIndexes = model.djcIndexes;                    % (:,1) uint64 % Constraint rows from `F` and `g`
            result.djcClauseNTerm = model.djcClauseNTerm;  % uint64 % Number of terms in each clause
            result.djcNumClauses = model.djcNumClauses;              % (:,1) uint64 % Number of clauses in each DJC
            result.djcDomain = model.djcDomain.struct();                      % (:,1) mosekdomain % Domains for each constraint block
            result.djcName = convertStringsToChars(model.djcName);
        
            result.use_init_sol = model.use_init_sol;

            % initial solution
            if isfield(model,"init_sol_x")
                result.init_sol_x = model.init_sol_x;
            end
            if isfield(model,"init_sol_sc")
                result.init_sol_sc = model.init_sol_sc;
            end

            % current solutions:
            result.sol_itr_available = model.sol_itr_available;
            if model.sol_itr_available
                result.sol_itr_prosta = model.sol_itr_prosta;
                result.sol_itr_solsta = model.sol_itr_solsta;
                result.sol_itr_x = model.sol_itr_x;
                result.sol_itr_y = model.sol_itr_y;
            end
     
            result.sol_bas_available = model.sol_bas_available;
            if model.sol_bas_available
                result.sol_bas_prosta = model.sol_bas_prosta;
                result.sol_bas_solsta = model.sol_bas_solsta;
                result.sol_bas_x = model.sol_bas_x;
                result.sol_bas_y = model.sol_bas_y;
            end

            result.sol_itg_available = model.sol_itg_available;
            if model.sol_itg_available
                result.sol_itg_prosta = model.sol_itg_prosta;
                result.sol_itg_solsta = model.sol_itg_solsta;
                result.sol_itg_x = model.sol_itg_x;
            end
        end
      end % methods
 end  
