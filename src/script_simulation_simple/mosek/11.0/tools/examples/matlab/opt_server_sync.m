%
%  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
%
%  File :      opt_server_sync.m
%
%  Purpose : Demonstrates solving a problem remotely using the OptServer.
%
%            url should be 'http://server:port' or 'https://server:port'
%            cert is the path to a TLS certificate, if using https
%

function opt_server_sync(url, cert)

% Here we can set up a model
model = mosekmodel(...
              objsense = "min", ...
              objective = [ 1, 2 ]', ...
              numvar = 2, ...
              F = [ 1 0 ; ...
                    0 1 ; ...
                    3 -1 ], ...
              domain = [ mosekdomain("greater than", rhs = 2), ... 
                         mosekdomain("greater than", rhs = 3), ... 
                         mosekdomain("less than", rhs = 1) ]);     


% Set up the certificate path, if using TLS, otherwise ignore
if exist('cert','var')
	param = ["MSK_SPAR_REMOTE_TLS_CERT_PATH", cert ];
else
	param = [];
end

% Optimize using the remote server
model.solve(param = param, optserver = url);

% Use the optimal solution
if model.hassolution("interior")
	xx = model.getsolution("interior", "x");
	disp(xx);
end

