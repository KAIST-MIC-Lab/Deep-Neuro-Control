function mosekcheck(options)
%
% MOSEK diagnostics script for MATLAB.
%
arguments
    options.debug logical = false;
end

%
% Basic system information
%
fprintf("Matlab version   : %s\n", version());
fprintf("Architecture     : %s\n", computer());


%
% Check if version is not too small
%
s = sscanf(version(),'%d.%d');
if s(1) < 9 || (s(1) == 9 && s(2) < 10)
    error(["MOSEK supports MATLAB 9.10 or newer. Installed version: %d.%d"], ...
          s(1), s(2))
end


%
% See if mosekenv is available.
%
c = which("mosekenv");
if length(c) == 0
    error("mosekenv is not included in the MATLAB path.");
else
    fprintf("MOSEK MeX path   : %s\n", c);
end


%
% Get version
%
try
    [major, minor, revision] = mosekenv("version");
catch ME 
    disp("*******************************************************");    
    disp("Error: could not load MOSEK MeX file.");
    disp("Very likely the folder with shared libraries is not in the environment");
    disp("variable PATH or you missed some operating-system-specific installation step.");    
    disp("Below is a more detailed debug info to help you identify the problem.");
    disp("See the interface manual:");
    disp("   https://docs.mosek.com/latest/matlabapi/install-interface.html#troubleshooting");
    disp("for typical troubleshooting tips.");
    disp("Share the full output of mosekcheck(debug = True) if contacting support.");
    disp("************ EXTENDED ERROR LOG BELOW ******************");
    rethrow(ME);
end

fprintf("MOSEK version    : %d.%d.%d\n", major, minor, revision);

%
% Attempt to check out licenses
%
pts = false;
try
    mosekenv("checkoutlicense", "pts");
    pts = true;
    fprintf("PTS license      : Yes.\n");
catch ME
    fprintf("PTS license      : No.\n");
end

pton = false;
try
    mosekenv("checkoutlicense", "pton");
    pton = true;
    fprintf("PTON license     : Yes.\n");
catch ME
    fprintf("PTON license     : No.\n");
end

if ~pts
    fprintf("Unable to check out any license.\n");
else
    fprintf("MOSEK works OK   : You can use MOSEK in MATLAB.\n");
end

if ~options.debug
    fprintf("\nTo obtain more information and debugging tips run:\n")
    fprintf("   mosekcheck(debug = true)\n");
end

mosekenv("checkinlicense");


if options.debug

%
% Testing a simple optimization problem
%
model = mosekmodel();
try
    disp(newline);
    disp("*************** SOLVER LOG OUTPUT: ******************");
    disp(newline);

    model.solve(quiet = false);
    fprintf("Test solve       : Success.\n");
    fprintf("MOSEK works OK   : You can use MOSEK in MATLAB.\n");
catch ME
    fprintf("Test solve       : Error.\n");
    fprintf("Error message    : %s\n", ME.message);

    disp(newline);
    disp("************* ERROR INFORMATION: ********************");    
    disp(newline);
    disp("Above is a detailed log output to help you identify the problem.");
    disp("See the interface manual:");
    disp("   https://docs.mosek.com/latest/matlabapi/install-interface.html#troubleshooting");
    disp("for typical troubleshooting tips.");
    disp("Share the full output of mosekcheck(debug = true) if contacting support.");

    if length(strfind(ME.message, "LICENSE")) > 0 || length(strfind(ME.message, "license")) > 0 || length(strfind(ME.message, "License")) > 0
        disp("This error is related to licensing.");
        disp("The line with `License path` below shows where MOSEK is expecting the license.");
        disp("Make sure the license is placed exactly there and is readable and up to date.");
        disp("Try restarting MATLAB if you installed a new license file.");
    end

end

disp(newline);
disp("************ EXTENDED MOSEK CONFIG ***************");
disp(newline);
%
% Prints the intro
%
mosekenv("intro");

mosekenv("checkinlicense");

end % if debug

end % mosekcheck
