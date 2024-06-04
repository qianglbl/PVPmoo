# PVPmoo
PVPmoo: A Parallel Variable Population Multi-Objective Optimizer

****************************

*** Copyright Notice ***

Parallel Variable Population Multi-Objective Optimizer (pvpmoo) 
Copyright (c) 2024, The Regents of the University of California,
through Lawrence Berkeley National Laboratory (subject to receipt of
any required approvals from the U.S. Dept. of Energy). All rights reserved.

If you have questions about your rights to use or distribute this software,
please contact Berkeley Lab's Intellectual Property Office at
IPO@lbl.gov.

NOTICE.  This Software was developed under funding from the U.S. Department
of Energy and the U.S. Government consequently retains certain rights.  As
such, the U.S. Government has been granted for itself and others acting on
its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
Software to reproduce, distribute copies to the public, prepare derivative 
works, and perform publicly and display publicly, and to permit others to do so.


****************************

! This program was developed by Ji Qiang (jqiang@lbl.gov) of the Lawrence Berkeley National Lab.

!.........................................................................

Some features include:

1) The population size varies from generation to generation.

2) The population is uniformly distributed to a number of parallel processors

3) The objective functions are attained from the ouptputs of an external serial program.

Ref: J. Qiang, "A parallel variable population multi-objective optimizer for
accelerator beam dynamics optimization," Nuclear Inst. and Methods in Physics Research,
A 1054 (2023) 168402.

The optimization input parameters are in pvpmoo.in.
The objective function file "objfunc.f90" should be modified according to
user's application program.

The outputs are obj.out (objective function values) and pop.out (control parameter values).

!.........................................................................

