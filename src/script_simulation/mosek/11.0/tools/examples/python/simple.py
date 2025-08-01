##
#  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
#
#  File :      simple.py
#
#  Purpose :   Demonstrates a very simple example using MOSEK by
#              reading a problem file, solving the problem and
#              writing the solution to a file.
##
import mosek
import sys

def streamprinter(msg):
    sys.stdout.write(msg)
    sys.stdout.flush()

if len(sys.argv) <= 1:
    print("Missing argument, syntax is:")
    print("  simple inputfile [ solutionfile ]")
else:
    with mosek.Task() as task:
        task.set_Stream(mosek.streamtype.log, streamprinter)

        # We assume that a problem file was given as the first command
        # line argument (received in `argv')

        task.readdata(sys.argv[1])

        # Solve the problem
        task.optimize()

        # Print a summary of the solution
        task.solutionsummary(mosek.streamtype.log)

        # If an output file was specified, save problem to a file
        if len(sys.argv) >= 3:
            # If using OPF format, these parameters will specify what to include in output
            task.putintparam(mosek.iparam.opf_write_solutions, mosek.onoffkey.on)
            task.putintparam(mosek.iparam.opf_write_problem, mosek.onoffkey.on)
            task.putintparam(mosek.iparam.opf_write_hints, mosek.onoffkey.off)
            task.putintparam(mosek.iparam.opf_write_parameters, mosek.onoffkey.off)

            task.writedata(sys.argv[2])
