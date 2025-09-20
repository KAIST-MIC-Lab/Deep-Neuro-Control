classdef moseklinmodel < handle
    % Encapsulates following problem:
    %
    %  min/max   c'x + cfix
    %  such that A'x = b
    %            l < x < u
    %
    properties
        modelName string

        % objective
        objSense     string
        c            (:,1) double
        cfix         double

        % Variables
        % Lengths must be the same and indicates the number of scalar variables
        numvar        uint32
        varnames      (:,1) string
        bx            (:,2) double

        % Constraints
        % At'x = b defines constraint matrix rows
        % Keep A in transposed for to allow cheaper row and columns adding
        connames      (:,1) string
        At            (:,:) double
        b             (:,1) double

        % At as components. These are set immediately before call to solve
        % because accessing sparse matrix components from c++ is super
        % inefficient!
        A_nrows  uint32
        A_ncols  uint32
        A_subi  (:,1) uint64
        A_subj  (:,1) uint64
        A_cofij (:,1) double

        % Initial solution
        init_sol_available logical
        init_sol_x   (:,1) double
        init_sol_slx (:,1) double
        init_sol_sux (:,1) double
        init_sol_y   (:,1) double
        init_sol_skx (:,1) int32
        init_sol_skc (:,1) int32

        % Result solutions
        sol_bas_available logical
        sol_bas_prosta    string
        sol_bas_solsta    string
        sol_bas_x         (:,1) double
        sol_bas_slx       (:,1) double
        sol_bas_sux       (:,1) double
        sol_bas_y         (:,1) double
        sol_bas_skx       (:,1) int32
        sol_bas_skc       (:,1) int32

        % latest information items
        info
    end

    methods
        function model = moseklinmodel(options)
            % MOSEKLINMODELCreate a new linear model object.
            %
            %    model = MOSEKLINMODEL(options...)
            %
            %    options:
            %    - name string (default: "model")
            %      Name of the model.
            %    - numvar uint32 (default: 0)
            %      Number of scalar variables.
            %    - varnames (:,1) string
            %      Variable names
            %    - blx (:,1) double
            %      Variable lower bound
            %    - bux (:,1) double
            %      Variable upper bound
            %    - objsense string (default: "min")
            %      Objective sense, "min", "minimize", "max" or "maximize"
            %    - objective (:,1) double
            %      Objective vector, this must match the given numvar.
            %    - objfixterm double
            %      Objective fixed scalar term.
            %    - A (:,:) double
            %      The constraint matrix.
            %    - b (:,1) double 
            %      The constraint right-hand constant vector
            %    - connames (:,1) string
            %      Constraint names
            %    - solution_x (:,1) double
            %      Initial x solution.
            %    - solution_y (:,1) double
            %      Initial y solution.
            %    - solution_slx (:,1) double
            %      Initial slx solution.
            %    - solution_sux (:,1) double
            %      Initial sux solution.
            %    - solution_skx (:,1) int32
            %      Initial skx solution.
            %    - solution_skc (:,1) int32
            %      Initial skc solution.
            %
            %    See also: APPENDVARS, APPENDCONS, OBJECTIVE
            arguments
                options.name          string       = "model"
                options.objsense      string       = "min"
                options.objective     (:,1) double
                options.objfixterm    double

                options.numvar        uint32
                options.varnames      (:,1) string
                options.blx           (:,1) double
                options.bux           (:,1) double

                options.A             (:,:) double
                options.b             (:,1) double
                options.connames      (:,1) string

                % Initial solution
                options.solution_x    (:,1) double
                options.solution_y    (:,1) double
                options.solution_slx  (:,1) double
                options.solution_sux  (:,1) double
                options.solution_skx  (:,1) string
                options.solution_skc  (:,1) string
            end

            model.sol_bas_available = false;
            model.init_sol_available = false;

            model.modelName = options.name;
            switch (options.objsense)
                case {"min","minimize"}
                    model.objSense  = "min";
                case {"max","maximize"}
                    model.objSense  = "max";
                otherwise
                    err = MException("moseklinmodel:DataError","Invalid objective sense");
                    throw(err);
            end

            if isfield(options,"numvar")
                numvar = options.numvar;
            elseif isfield(options,"blx")
                numvar = size(options.blx,1);
            elseif isfield(options,"bux")
                numvar = size(options.bux,1);
            elseif isfield(options,"varnames")
                numvar = size(options.varnames,1);
            else
                numvar = 0;
            end

            model.numvar = numvar;

            if isfield(options,"blx")
                if numvar ~= size(options.blx,1)
                    err = MException("moseklinmodel:DataError","Invalid blx size");
                    throw(err);
                end
            elseif isfield(options,"bux")
                if numvar ~= size(options.bux,1)
                    err = MException("moseklinmodel:DataError","Invalid bux size");
                    throw(err);
                end
            elseif isfield(options,"varnames")
                if numvar ~= size(options.varnames,1)
                    err = MException("moseklinmodel:DataError","Invalid varnames size");
                    throw(err);
                end
            end

            model.bx = zeros(numvar,2);
            model.bx(:,1) = -inf;
            model.bx(:,2) = inf;

            if isfield(options,"blx")
                model.bx(:,1) = options.blx;
            end
            if isfield(options,"bux")
                model.bx(:,2) = options.bux;
            end

            if isfield(options,"objective")
                if size(options.objective,1) > numvar
                    err = MException("moseklinmodel:DataError","Objective vector exceeds number of variables");
                    throw(err);
                end
                model.c = options.objective;
            else
                model.c = zeros(numvar,1);
            end

            if isfield(options,"varnames")
                model.varnames = options.varnames;
            end

            if isfield(options,"objfixterm")
                model.cfix = options.objfixterm;
            else
                model.cfix = 0.0;
            end

            nrow = 0;
            if isfield(options,"A")
                nrow = size(options.A,1);
            elseif isfield(options,"b")
                nrow = size(options.b,1);
            end

            if isfield(options,"A")
                if size(options.A,2) > numvar
                    err = MException("moseklinmodel:DataError","Coefficient matrix width exceeds number of variables");
                    throw(err);
                elseif isfield(options,"b") && size(options.A,1) ~= size(options.b,1)
                    err = MException("moseklinmodel:DataError","Mismatching A and b sizes");
                    throw(err);
                end
                model.At = sparse(options.A');
            else
                model.At = sparse(numvar,nrow);
            end
            model.A_nrows = nrow;


            if isfield(options,"b")
                if size(options.b,1) ~= nrow
                    err = MException("moseklinmodel:DataError","Mismatching A and b sizes");
                    throw(err);
                end
                model.b = options.b;
            else
                model.b = zeros(nrow,1);
            end

            if isfield(options,"connames")
                if size(options.connames,1) ~= model.A_nrows
                    err = MException("moseklinmodel:DataError","Invalid connames size");
                    throw(err);
                end
                model.connames = options.connames;
            end

            if isfield(options,"solution_x")
                model.init_sol_x   = options.solution_x;
                model.init_sol_available = true;
            end
            if isfield(options,"solution_y")
                model.init_sol_y   = options.solution_y;
                model.init_sol_available = true;
            end
            if isfield(options,"solution_slx")
                model.init_sol_slx = options.solution_slx;
                model.init_sol_available = true;
            end
            if isfield(options,"solution_sux")
                model.init_sol_sux = options.solution_sux;
                model.init_sol_available = true;
            end
            if isfield(options,"solution_skx")
                model.init_sol_skx = arrayfun(@str_to_sk,options.solution_skx);
                model.init_sol_available = true;
            end
            if isfield(options,"solution_skc")
                model.init_sol_skc = arrayfun(@str_to_sk,options.solution_skc);
                model.init_sol_available = true;
            end
        end % constructor



        function varname(model,indexes,names)
            % VARNAME  Set variables names
            %
            %    VARNAME(indexes,names)
            %
            %    Arguments:
            %    - indexes (:,1) int32
            %      Indexes of variables to set names for.
            %    - names (:,1) string
            %      List of names. The length must match indexes.
            arguments
                model moseklinmodel
                indexes (:,1) uint32
                names   (:,1) string
            end

            if max(indexes) > model.numvar || min(indexes) == 0
                err = MException("moseklinmodel:DataError","Invalid variable indexes");
                throw(err);
            elseif size(names,1) ~= size(indexes,1)
                err = MException("moseklinmodel:DataError","Mismatching lengths of indexes and names");
                throw(err);
            end

            model.varnames(indexes,1) = names;
        end % varname

        function conname(model,indexes,names)
            % CONNAME  Set constraint names
            %
            %    CONNAME(indexes,names)
            %
            %    Arguments:
            %    - indexes (:,1) int32
            %      Indexes of constraints to set names for.
            %    - names (:,1) string
            %      List of names. The length must match indexes.
            arguments
                model moseklinmodel
                indexes (:,1) uint32
                names   (:,1) string
            end

            numcon = size(model.b,1);

            if max(indexes) > numcon || min(indexes) == 0
                err = MException("moseklinmodel:DataError","Invalid constraint indexes");
                throw(err);
            elseif size(names,1) ~= size(indexes,1)
                err = MException("moseklinmodel:DataError","Mismatching lengths of indexes and names");
                throw(err);
            end

            model.connames(indexes,1) = names;
        end % conname

        function result = appendvars(model,num,options)
            % APPENDVARS  Add new scalar variables to the model.
            %
            %    vars = APPENDVARS(model,num,options...)
            %
            %    Append num scalar variables to the model, returning the indexes of the
            %    added variables.
            %
            %    Options:
            %    - bl (:,1) double
            %      Vector of variable lower bounds
            %    - bu (:,1) double
            %      Vector of variable upper bounds
            %    - c  (:,1) doble
            %      Vector of objective coefficients for the new variables
            %    - A  (:,:) double
            %      Constraint matrix nonzeros for the new columns
            %    - names (:,1) string
            %      Names of the new variables
            %
            %    See also: APPENDCONS, Objective
            arguments
                model
                num              uint32 = 0
                options.bl (:,1) double
                options.bu (:,1) double
                options.c  (:,1) double
                options.A  (:,:) double
                options.names (:,1) string
            end

            if isfield(options,"bl")
                if size(options.bl,1) ~= num
                    err = MException("moseklinmodel:DataError","Mismatching num and bl");
                    throw(err);
                end
            end
            if isfield(options,"bu")
                if size(options.bu,1) ~= num
                    err = MException("moseklinmodel:DataError","Mismatching num and bu");
                    throw(err);
                end
            end

            if isfield(options,"A") 
                if size(options.A,2) ~= num
                    err = MException("moseklinmodel:DataError","Mismatching num and A");
                    throw(err);
                end
            end
            
            if isfield(options,"c") 
                if size(options.c,1) ~= num
                    err = MException("moseklinmodel:DataError","Mismatching num and c");
                    throw(err);
                end
            end
            
            if isfield(options,"names")
                if size(options.names,1) ~= num
                    err = MException("moseklinmodel:DataError","Mismatching num and names");
                    throw(err);
                end
            end

            oldnumvar = model.numvar;
            if num > 0
                model.numvar = model.numvar + num;
                model.bx(oldnumvar+1:model.numvar,:) = 0.0;
            end
    
            if isfield(options,"bl")
                bl = options.bl;
            else 
                bl = -inf(num,1);
            end
            if isfield(options,"bu")
                bu = options.bu;
            else 
                bu = inf(num,1);
            end

            if isfield(options,"A")
                A = options.A';
            else
                A = sparse(num, size(model.At,2));
            end

            if isfield(options,"c")
                c = options.c;
            else
                c = zeros(num,1);
            end

            model.bx(oldnumvar+1:model.numvar,1)  = bl;
            model.bx(oldnumvar+1:model.numvar,2)  = bu;
            model.At(oldnumvar+1:oldnumvar+num,:) = A;
            model.c(oldnumvar+1:oldnumvar+num,1)  = c;


            if num > 0
                result = uint32([ oldnumvar+1 : oldnumvar+num ]);
            else
                result = uint32([]);
            end            
            if isfield(options,"names")
                n = size(model.varnames,1);
                model.varnames(n+1:n+num) = options.names;
            end
        end % appendvars

        function appendcons(model,A,b,options)
            % APPENDCONS  Append a new constraint block to the model.
            %
            %    APPENDCONS(model,options...)
            %
            %    Append a block of constraints.
            %
            %    Arguments:
            %    - A       (:,:) double
            %      Constraint matrix. If argument x is given, the number of columns in A
            %      must be equal the number of columns in x.
            %    - b       (:,1) double
            %      Constraint left-hand side constant term.
            %    Options:
            %    - x       (:,1) uint32
            %      Optional. Specifies the variable vector for the constraints.
            %    - names (:,1) string
            %      Names of the constraints
            arguments
                model
                A               (:,:) double %
                b               (:,1) double %
                options.x       (:,1) uint32 % x vector
                options.names   (:,1) string
            end

            oldnumcon = size(model.At,2);

            [nrows,ncols] = size(A);
            if isfield(options,"names")
                if size(options.names,1) ~= nrows
                    err = MException("moseklinmodel:DataError","Mismatching A and names");
                    throw(err);
                end
            end

            if nrows ~= size(b,1)
                err = MException("moseklinmodel:DataError","Mismatching A and b");
                throw(err);
            end

            if isfield(options,"x")
                if ncols > size(options.x,1)
                    err = MException("moseklinmodel:DataError","Mismatching A and x");
                    throw(err);
                elseif max(options.x) > model.numvar
                    err = MException("moseklinmodel:DataError","Element in x is out of bounds");
                    throw(err);
                end

                [Ai,Aj,Ac] = find(A);
                Aj = x(Aj);
                At = sparse(ncols,nrows,Aj,Ai,Ac);
            else
                if ncols > model.numvar
                    err = MException("moseklinmodel:DataError","Invalid number of columns in A");
                    throw(err);
                end

                At = A';
            end

            model.At(:,oldnumcon+1:oldnumcon+nrows) = At;
            model.b(oldnumcon+1:oldnumcon+nrows,1)  = b;
            model.A_nrows = model.A_nrows + nrows;

            if nrows > 0
                result = uint64([oldnumcon+1:oldnumcon+nrows]);
            else
                result = uint64.empty(0,1);
            end

            if isfield(options,"names")
                n = size(model.connames,1);
                model.connames(n+1:n+size(options.names,1)) = options.names;
            end
        end % appendcons

        function numcon = getnumcon(model)
            % GET_NUMCON  Get number of linear constraints.
            %
            %    numcon = GET_NUMCON(model)
            %    
            numcon = size(model.At,2);
        end

        function numvar = getnumvar(model)
            % GET_NUMVAR  Get number of linear variables.
            %
            %    numvar = GET_NUMVAR(model)
            %    
            numvar = model.numvar;
        end

        function setb(model,b,options)
            % SET_B  Update the constraint b vector or a part of it.
            %
            %    SET_B(model,b,options)
            %
            %    Update the constraint b vector.
            %
            %    Arguments:
            %    - b (:,1) double 
            %      New values to put into the constraint b vector.
            %    Options:
            %    - first uint32 = 1
            %      Index of the first element if the b vector to update.
            arguments
                model moseklinmodel
                b (:,1) double
                options.first uint32
            end

            if isfield(options,"first")
                first = options.first;
            else
                first = 1;
            end

            n = size(b,1);
            if first+n-1 > size(model.b,1)
                err = MException("moseklinmodel:DataError","Range out of bounds");
                throw(err);
            end
            model.b(first:first+n-1,1) = b;
        end

        function setvarbounds(model,options)
            % SET_VAR_BOUNDS  Update the variable bounds.
            % 
            %    SET_BAR_BOUNDS(model,options)
            %
            %    Options:
            %    - bl (:,1) double
            %      Lower bound values.
            %    - bu (:,1) double
            %      Upper bound values.
            %    - first uint32
            %      Index of the first variable to update.
            arguments
                model moseklinmodel
                options.bl (:,1) double
                options.bu (:,1) double
                options.first uint32
            end
            
            if isfield(options,"first")
                first = options.first;
            else
                first = 1;
            end

            if isfield(options,"bl")
                b = options.bl;
                n = size(b, 1);
                if first+n-1 > size(model.bx,1)
                    err = MException("moseklinmodel:DataError","Range of bl out of bounds");
                    throw(err);
                end
                model.bx(first:first+n-1,1) = b;
            end

            if isfield(options,"bu")
                b = options.bl;
                n = size(b, 1);
                if first+n-1 > size(model.bx,1)
                    err = MException("moseklinmodel:DataError","Range of bu out of bounds");
                    throw(err);
                end
                model.bx(first:first+n-1,2) = b;
            end

        end


        function setrows(model,A,options) 
            % SET_ROWS  Set row non-zeros of the model constraint matrix.
            %
            %    SET_ROWS(model,A,options)
            %    - model moseklinmodel
            %    - A (:,:) double 
            %      A matrix row values. The number of columns in A must match
            %      the model's number of columns.
            %    Options:
            %    - first uint32
            %      First row to update.
            arguments
                model moseklinmodel
                A     (:,:) double
                options.first uint32
            end

            if isfield(options,"first")
                first = options.first;
            else 
                first = 1;
            end

            numvar = size(model.At,1);  
            numcon = size(model.At,2);
        

            if numvar ~= size(A,2)
                err = MException("moseklinmodel:DataError","Incorrect number of columns in A");
                throw(err);
            end

            n = size(A,1);
            if first+n-1 > numvar 
                err = MException("moseklinmodel:DataError","Range out of bounds");
                throw(err);
            end

            model.At(:,first:first+n-1) = A';
        end

        function setcolumns(model,A,options)
            % SET_COLUMNS  Set column non-zeros of the model constraint matrix.
            %
            %    SET_COLUMNS(model,A,options)
            %    - model moseklinmodel
            %    - A (:,:) double
            %      A matrix column values. The number of rows in A must match
            %      the model's number of rows.
            %    Options:
            %    - first uint32
            %      First column to update.
            arguments
                model moseklinmodel
                A     (:,:) double
                options.first uint32
            end

            if isfield(options,"first")
                first = options.first;
            else 
                first = 1;
            end

            numvar = size(model.At,1);
            numcon = size(model.At,2);

            if numcon ~= size(A,1)
                err = MException("moseklinmodel:DataError","Incorrect number of columns in A");
                throw(err);
            end
           
            n = size(A,2);
            if first+n-1 > numcon 
                err = MException("moseklinmodel:DataError","Range out of bounds");
                throw(err);
            end

            model.At(first:first+n-1,:) = A';
        end

        function setc(model,c,options)
            % SETC  Set the objective function coefficients
            % 
            %    SETC(model,c,options)
            %
            %    - c (:,1) double
            %      The objective vector or parital vector. Note that this must
            %      have length at most the number of variables in the model.
            %      Only the coefficients in the range (first:first-1+n), where
            %      n = size(c,1).
            %    Options:
            %    - first uint32
            %      Index of the first variable whose coefficient to change
            arguments
                model moseklinmodel
                c (:,1) double
                options.first   uint32
            end

            if isfield(options,"first")
                first = options.first;
            else
                first = 1;
            end

            if size(c,1) + first - 1 > model.numvar
                err = MException("moseklinmodel:DataError","Objective vector exceeds number of variables in the model");
                throw(err);
            end
            model.c(first:first-1+size(c,1),1) = c;
        end

        function objective(model, objsense, c, options)
            % OBJECTIVE  Set the objective function
            % 
            %    OBJECTIVE(model,objsense,c,options)
            %    - objsense string
            %      The objective sense, "min", "minimize", "max" or "maximize"
            %    - c (:,1) double
            %      The objective vector. Note that this must match the number of variables in the model.
            %    Options:
            %    - objfixterm double
            %      Objective fixed term
            arguments
                model moseklinmodel
                objsense string
                c (:,1) double
                options.objfixterm double
            end

            if size(c,1) ~= model.numvar
                err = MException("moseklinmodel:DataError","Objective vector has wrong length");
                throw(err);
            end
            model.c = c;

            switch objsense
                case {"min","minimize"}
                    model.objSense = "min";
                case {"max","maximize"}
                    model.objSense = "max";
                otherwise
                    err = MException("moseklinmodel:DataError","Objective sense is invalid");
                    throw(err);
            end

            if isfield(options,"objfixterm")
                model.cfix = options.objfixterm;
            end
        end % objective

        function write(model,filename,options)
            % WRITE  Write problem to a file
            %
            %    WRITE(model,filename,options...)
            %
            %    - filename string
            %      File name to write the problem to. The format is determined 
            %      by the file name extension.
            arguments
                model moseklinmodel
                filename string
                options.param  (1,:) string
            end

            if isfield(options,"param")
                mosekenv("writelinear",model,filename,options.param);
            else
                mosekenv("writelinear",model,filename);
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
                model           moseklinmodel
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
            
            [A_subj,A_subi,A_cofij] = find(model.At);
            [model.A_ncols,model.A_nrows] = size(model.At);
            model.A_subj  = A_subj-1;
            model.A_subi  = A_subi-1;
            model.A_cofij = A_cofij;

            if isfield(options,"param")
                [sol,info] = mosekenv("solvelinear",model,nosolve,write_to_file,licfile,liccode,logfile,quiet,optserver,analyze,options.param); 
            else 
                [sol,info] = mosekenv("solvelinear",model,nosolve,write_to_file,licfile,liccode,logfile,quiet,optserver,analyze);
            end

            if ~ nosolve 
                model.init_sol_available = false;

                % erase old user-provided solution
                if isfield(model,"init_sol_x")
                    rmfield(model,"init_sol_x");
                end
                if isfield(model,"init_sol_y")
                    rmfield(model,"init_sol_y");
                end
                if isfield(model,"init_sol_slx");
                    rmfield(model,"init_sol_slx");
                end
                if isfield(model,"init_sol_sux");
                    rmfield(model,"init_sol_sux");
                end
                if isfield(model,"init_sol_skx");
                    rmfield(model,"init_sol_skx");
                end
                if isfield(model,"init_sol_skc");
                    rmfield(model,"init_sol_skc");
                end
                
                model.sol_bas_available = false;
           
                model.info = info;
                if size(sol,1) == 8
                    model.sol_bas_available = true;
                    bas_prosta = sol{1};
                    bas_solsta = sol{2};
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

                    model.sol_bas_x   = sol{3};
                    model.sol_bas_slx = sol{4};
                    model.sol_bas_sux = sol{5};
                    model.sol_bas_y   = sol{6};
                    model.sol_bas_skx = sol{7}; 
                    model.sol_bas_skc = sol{8};
                end
            end
        end % solve

        function [result,prosta,solsta] = hassolution(model)
            % HASSOLUTION  Check is a specific solution is available and returns it status.
            % 
            %    exists = HASSOLUTION(model)
            %    
            %    Only check if the solution is available.
            %
            %    [exists,prosta,solsta] = HASSOLUTION(model)
            %             
            %    Checks of the solution is available and returns problem and solutions
            %    status if it is.
            %
            %    See also: GETSOLUTION, SOLVE
            result = model.sol_bas_available;
            prosta = model.sol_bas_prosta;
            solsta = model.sol_bas_solsta;
        end % hassolution

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
            %     values), "slx", "sux" (variable bound dual values)
            %     "skc" (constraint status) or "skx" (variable status).
            %   - value (:,1)
            %     Vector of solution values. For numerical values this is 
            %     a vector of solution values, for part="sk" this is a vector of 
            %     logicals indicating the the individual constraints are binding.
            arguments
                model moseklinmodel
                part string
                value (:,1) 
            end

            switch part 
                case "x"
                    model.init_sol_x = value;
                    model.init_sol_available = true;
                case "y"
                    model.init_sol_y = value;
                    model.init_sol_available = true;
                case "slx"
                    model.init_sol_slx = value;
                    model.init_sol_available = true;
                case "sux"
                    model.init_sol_sux = value;
                    model.init_sol_available = true;
                case "skx"
                    model.init_sol_skx = arrayfun(@str_to_sk,value);
                    model.init_sol_available = true;
                case "skc"
                    model.init_sol_skc = arrayfun(@str_to_sk,value);
                    model.init_sol_available = true;
                otherwise
                    err = MException("mosekmodel:KeyError","Invalid solution part selector key");
                    throw(err); 
            end
        end % setsolution

        function [result,prosta,solsta] = getsolution(model,part) 
            % GETSOLUTION  Fetch a solution.
            %
            %   [result,prosta,solsta] = GETSOLUTION(model,part)
            % 
            %   Fetch the solution values, problem status and solution status.
            %   Not that if the solution is a certificate, this will be
            %   indicated by the solsta value.
            %   
            %   Arguments:
            %   - part string (default : "x")
            %     One of the values 
            %     - "x" indicating primal variable values, and 
            %     - "y" indicating dual constraint values
            %     - "slx" indicating dual slack variable solution related to lower bound on variable
            %     - "sux" indicating dual slack variable solution related to upper bound on variable
            %     - "skx" indicating variable status
            %     - "skc" indicating constraint status
            %   Returns:
            %   - result (:,1) double or (:.1) string
            %     The actual solution part value. For part="skx" or part="skc"
            %     this is an array of strings, otherwise an array of doubles.
            %   - solsta string
            %     String indicating the solution status
            %   - prosta string
            %     String indicating the problem status
            %
            %   See also: HASSOLUTION, SOLVE
            arguments
                model moseklinmodel
                part  string = "x"
            end

            if model.sol_bas_available                
                keys = [ "UNK"; "BAS"; "SUPBAS"; "LOW"; "UPR"; "FIX"; "INF" ];
                switch part
                    case "x"
                        result = model.sol_bas_x;
                    case "slx"
                        result = model.sol_bas_slx;
                    case "sux"
                        result = model.sol_bas_sux;
                    case "y"
                        result = model.sol_bas_y;
                    case "skx"
                        result = keys(model.sol_bas_skx+1);
                    case "skc"
                        result = keys(model.sol_bas_skc+1);
                    otherwise
                        err = MException("mosekmodel:KeyError","Invalid solution part selector key");
                        throw(err); 
                end
                prosta = model.sol_bas_prosta;
                solsta = model.sol_bas_solsta;
            else
                err = MException("mosekmodel:SolutionError","Basic solution is not available");
                throw(err); 
            end
        end
    end
end


function key = str_to_sk(s)
    switch s
        case {"unk", "UNK"}
            key = 0;
        case {"bas", "BAS"}
            key = 1;
        case {"supbas", "SUPBAS"}
            key = 2;
        case {"low", "LOW"}
            key = 3;
        case {"upr", "UPR"}
            key = 4;
        case {"fix", "FIX"}
            key = 5;
        case {"inf", "INF"}
            key = 6;
        otherwise
            key = 0;
    end
end
