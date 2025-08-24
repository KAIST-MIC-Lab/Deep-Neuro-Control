%%
%  Copyright: Copyright (c) MOSEK ApS, Denmark. All rights reserved.
%
%  File:      facility_location.m
%
%  Purpose: Demonstrates a small one-facility location problem.
%
%  Given 10 customers placed in a grid we wish to place a facility
%  somewhere so that the total sum of distances to customers is
%  minimized.
%
%  The problem is formulated as a conic optimization problem as follows.
%  Let f=(fx,fy) be the (unknown) location of the facility, and let
%  c_i=(cx_i,cy_i) be the (known) customer locations; then we wish to
%  minimize
%      sum_i || f - c_i ||
%  where
%      ||.||
%  denotes the euclidian norm.
%  This is formulated as
%
%  minimize   sum(d_i)
%  such that  d_i ^ 2 > tx_i ^ 2 + ty_i ^ 2, for all i (a.k.a. (d_i,tx_i,ty_i) in C^3_r)
%             tx_i = cx_i - fx, for all i
%             ty_i = cy_i - fy, for all i
%             d_i > 0, for all i
%%
function [xx,prosta,solsta] = facility_location()
    customerloc = [12.0  2.0 ; ... 
                   15.0 13.0 ; ... 
                   10.0  8.0 ; ... 
                    0.0 10.0 ; ... 
                    6.0 13.0 ; ... 
                    5.0  8.0 ; ... 
                   10.0 12.0 ; ... 
                    4.0  6.0 ; ... 
                    5.0  2.0 ; ... 
                    1.0 10.0];
    N = size(customerloc,1);
    model = mosekmodel(...
        name="FacilityLocation",...
        numVar=2+N);
    % variable(1:N) are the distances from the facility to each customerloc
    % variable(N+1:N+2) is the location of the facility

    % distance per customer
    appendcons(model,...
               F = sparse([1:N*3], reshape([ [1:N] ; repmat([N+1;N+2],1,N)],N*3,1)', ones(1,N*3)),...
               g = reshape([ zeros(1,N) ; -customerloc' ], 3*N,1),...
               domain=mosekdomain("qcone",n=N,dim=3));
    % Summed distances
    objective(model,...
              "minimize", ...
              [ 0; 0; ones(N,1)] );

    solve(model);
    if hassolution(model,"interior") 
        [xx,prosta,solsta] = getsolution(model,"interior","x");
    end
end
