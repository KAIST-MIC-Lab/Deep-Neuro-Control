classdef mosekterm < handle
    properties
        rowidx (:,1) uint64
        dom (:,1) mosekdomain
        numdom uint64
    end

    methods 
        function term = mosekterm(rowidx,dom)
            arguments
                rowidx (:,1) uint64
                dom    (:,1) mosekdomain
            end
            
            if size(dom,1) == 1
                if dom.n * dom.domdim ~= size(rowidx,1)
                    error('Size of domain does not match expression Ax+b')
                end
            elseif sum([dom.n] .* [dom.domdim]) ~= size(rowidx,1)
                error('Size of domain does not match expression Ax+b')
            end

            term.rowidx = rowidx;
            term.dom    = dom;
            term.numdom = size(dom,1);
        end
    end
end
