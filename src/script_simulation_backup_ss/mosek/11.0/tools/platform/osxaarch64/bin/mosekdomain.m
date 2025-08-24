classdef mosekdomain < handle
    properties
        domtype string
        domdim  uint64
        n       uint64
        domofs  (:,1) double
        dompar  (:,1) double

        nelm    uint64
    end

    methods
        function result = mosekdomain(which,options)
            % MOSEKDOMAIN  Create a new constraint domain
            %
            %    dom = MOSEKDOMAIN(which,options...)
            % 
            %    Create a new domain of the type indicated by which.
            %    Arguments:
            %    - which string
            %      Defines the domain type, one of the strings:
            %      - "rminus" ( or "less than", "lt", "nonpositive", "r-", "<=")
            %        Nonpositive values. By default dim=1.
            %      - "rplus" (or "greater than", "gt", "nonnegative", "r+", ">=")
            %        Nonnegative values. By default dim=1.
            %      - "r" (or "unbounded")
            %        Unbounded values. By default dim=1.
            %      - "zero" (or "equal", "eq", "fixed", "zero", "==")
            %        Zero. By default dim=1.
            %      - "qcone" (or "quadratic cone")
            %        Second order cone. The dim argument must be given and dim >= 2.
            %      - "rqcone" (or "rotated quadratic cone")
            %        Rotated second order cone. The dim argument must be given and dim >= 3.
            %      - "exp" (or "exponential cone")
            %        Exponential cone. Dimension is always exactly 3.
            %      - "dexp" (or "dual exponential cone")
            %        Dual exponential cone. Dimension is always exactly 3.
            %      - "pow" (or "power cone")
            %        The power cone. The alpha and dim arguments are required.
            %      - "dpow" (or "dual power cone")
            %        The dual power cone. The alpha and dim arguments are required.
            %      - "geomean" (or "geometric mean cone")
            %        Geometric mean cone. The dim argument is required.
            %      - "dgeomean" (or "dual geometric mean cone")
            %        Geometric mean cone. The dim argument is required.
            %    Options:
            %    - rhs (:,1) double
            %      Domain offset. This is a constant vector added to the domain, so for 
            %      example, for which="equal", this would simply be the right-hand-side of 
            %      the constraint.
            %    - dim uint64
            %      Dimension of the cone. For "exp" and "dexp" it is implied. For "rminus",
            %      "rplus","r", and "zero" it is by default 1, but if a value larger than 1 
            %      is given, it will converted from n dim-sized domains into n*dim 1-sized 
            %      domains. For the remaining cones it is required.
            %    - n uint64 (default: 1)            
            %      Number of cones of the given type and dimension. 
            %    - alpha (:,1) double
            %      Required by "pow" and "dpow" cones, unused for all others. These are the 
            %      powers for the cone. The length must be at most dim-1. If the length 
            %      is 1, it is interpreted as [alpha, 1.0-alpha]. The vector is normalized 
            %      since elements should sum to 1.0, so the vector we actually use 
            %      is alpha/sum(alpha).
            arguments
                which string = ""
                options.rhs (:,1) double
                options.dim uint64
                options.n   uint64 = 1
                options.alpha (:,1) double
            end

            if which == "" 
                return
            end

            if isfield(options,"n")
                result.n = options.n;
            else
                result.n = 1;
            end

            switch which
                case {"less","less than","lt","nonpositive","greater","greater than","gt","nonnegative","unbounded","equal","equals","eq","fixed","r","rplus","rminus","zero","r+","r-","==","<=",">="}
                    switch which
                        case {"less", "less than", "lt", "nonpositive","rminus","r-", "<="}
                            result.domtype = "rminus";
                        case {"greater", "greater than", "gt", "nonnegative", "rplus", "r+", ">=" }
                            result.domtype = "rplus";
                        case {"unbounded","r"}
                            result.domtype = "r";
                        case {"fixed","zero","eq","equal","equals","=="}
                            result.domtype = "zero";
                    end

                    if isfield(options,"dim") && isfield(options,"n")
                        result.n = options.dim * result.n;
                    elseif isfield(options,"dim")
                        result.n = options.dim;
                    elseif isfield(options,"n")
                        result.n = options.n;
                    elseif isfield(options,"rhs")
                        result.n = size(options.rhs,1);
                    else
                        result.n = 1;
                    end
                    result.domdim = 1;

                case {"qcone","quadratic cone"}
                    result.domtype = "qcone";

                    if ~ isfield(options,"dim")
                        err = MException("mosekdomain:MissingArgument","Missing 'dim' argument for domain");
                        throw(err);
                    elseif options.dim < 2
                        err = MException("mosekdomain:InvalidArgument","Invalid 'dim' argument for quadratic cone domain");
                        throw(err);
                    end

                    result.domdim = options.dim;
                case {"rqcone", "rotated qcone","rotated quadratic cone"}
                    result.domtype = "rqcone";
                    if ~ isfield(options,"dim")
                        err = MException("mosekdomain:MissingArgument","Missing 'dim' argument for domain");
                        throw(err);
                    elseif options.dim < 3
                        err = MException("mosekdomain:InvalidArgument","Invalid 'dim' argument for rotated quadratic cone domain");
                        throw(err);
                    end
                    result.domdim = options.dim;
                case {"exp","exponential cone"}
                    result.domtype = "exp";
                    if isfield(options,"dim") && options.dim ~= 3
                        err = MException("mosekdomain:InvalidArgument","Invalid 'dim' argument for exponential cone domain");
                        throw(err);
                    end
                    result.domdim = 3;
                case {"dexp","dual exponential cone"}
                    result.domtype = "dexp";
                    if isfield(options,"dim") && options.dim ~= 3
                        err = MException("mosekdomain:InvalidArgument","Invalid 'dim' argument for exponential cone domain");
                        throw(err);
                    end
                    result.domdim = 3;
                case {"pow","dpow","power cone","dual power cone","dual pow"}
                    switch which
                        case {"pow","power cone"}
                            result.domtype = "pow";
                        case {"dpow","dual power cone","dual pow"}
                            result.domtype = "dpow";
                    end

                    if ~ isfield(options,"alpha")
                        err = MException("mosekdomain:InvalidArgument","Missing argument 'alpha' for power cone domain");
                        throw(err);
                    elseif size(options.alpha,1) == 1 
                        alpha = [ options.alpha 1.0-options.alpha ];
                    else
                        alpha = options.alpha;
                    end
                    if ~ isfield(options,"dim")
                        err = MException("mosekdomain:MissingArgument","Missing argument 'dim' for power cone domain");
                        throw(err);
                    end
                   
                    if size(alpha,1) >= options.dim 
                        err = MException("mosekdomain:InvalidArgument","Invalid argument 'alpha' for power cone domain");
                        throw(err);
                    end
                    
                    result.domdim = options.dim;
                    result.dompar = alpha;
                case {"geomean","dual geomean","dgeomean","geometric mean cone","dual geometric mean cone"}
                    switch which 
                        case {"geomean","geometric mean cone"}
                            result.domtype = "geomean";
                        case {"dual geomean","dgeomean","dual geometric mean cone"}
                            result.domtype = "dgeomean";
                    end

                    if ~ isfield(options,"dim")
                        err = MException("mosekdomain:MissingArgument","Missing argument 'dim' for geometric mean cone domain");
                        throw(err);
                    end
                    result.domdim = options.dim;
                case {"svecpsd","scaled vector psd"}
                    if ~ isfield(options,"dim")
                        err = MException("mosekdomain:MissingArgument","Missing argument 'dim' for scaled vector PSD cone domain");
                        throw(err);
                    end
                    d = floor((sqrt(1+ 8*options.dim)-1)/2);
                    if options.dim ~= d*(d+1)/2
                        err = MException("mosekdomain:InvalidArgument","Invalid value of argument 'dim' for scaled vector PSD cone domain");
                        throw(err);
                    end

                    result.domdim = options.dim;
                otherwise
                    if size(which,1) == 0
                        err = MException("mosekmodel:InvalidDomain","<missing>");
                    else
                        err = MException("mosekmodel:InvalidDomain",which);
                    end
                    throw(err);
            end % switch

            if isfield(options,"rhs")
                if size(options.rhs,1) ~= result.domdim * result.n
                    err = MException("mosekdomain:InvalidArgument","Right-hand side length does not match domain size");
                    throw(err);
                end
                result.domofs = options.rhs;
                result.nelm = size(options.rhs,1);
            else 
                result.domofs = zeros(result.domdim * result.n,1);
                result.nelm = result.domdim * result.n;
            end
        end
        function result = numrow(obj)
            result = obj.n * obj.domdim;
        end

        function result = struct(obj)
            [h,w] = size(obj);
            if h ~= 1 | w ~= 1 
                result = arrayfun(@struct,obj);
            else
                result.domtype = convertStringsToChars(obj.domtype);
                result.domdim  = obj.domdim;
                result.n       = obj.n;    
                result.domofs  = obj.domofs;
                result.dompar  = obj.dompar;
                result.nelm    = obj.nelm;
            end
        end
    end % methods
end
