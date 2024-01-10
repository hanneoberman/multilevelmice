## 0README.txt

This file simply describes a number of examples and test scripts
related to the optimx package. 

Note that even methods that are generally quite good still show
failures in some of these examples.

Most examples can be run in sequence by starting R IN THIS DIRECTORY
and then issuing the command 

source("runex.R")

The examples fhop.R and HessQuality.R need user input and should be
source'd separately.

John C. Nash 2023-06-13

The files are:

3Rosen.R -- this illustrates three different extensions of the famous
            Rosenbrock banana-shaped valley problem which has just 2
            parameters. The three extensions to more parameters give
            different function, gradient and hessian values.

argclash.R -- Illustration of some issues when the objective function
            uses information from dot arguments and the name is the
            same as one that is an argument of a called function, 
            e.g., gradient approximation

axsearch-Ex.R -- a simple illustration of the use and output of the
            axsearch() function to get "tilt" and "radius of curvature"
            (roc) in each parameter direction at a point on a function 
            surface

brown_test -- A quite nasty function that has the name "Brown" but does
            not seem related to the Brown Badly Scaled, Brown and Dennis,
            or Brown Almost Linear functions. The current code only
            has the function, without gradient or Hessian, and a number 
            of methods fail. A reference to its origin, as well as code
            for gradient and Hessian would be welcome. 

broydt_test  --  Test opm on the Broyden tridiagonal function (taken
            from package funconstrain)

chen_test    --  Test opm on the Chen function. It would be helpful to
            have a reference to the origin of this function.

chenlog_test --  Test opm on the Chen function (log-transformed parameters)

cyq_test     --  Test opm on the Fletcher Chebyquad problem of different sizes

dropnmbk.R -- A little script to show that method nmbk fails when 
            started with a parameter on a bound.

fhop.R -- To run many methods on some or all of the problems in the
            funconstrain package with upper, lower, both, or no 
            bounds.

froth_test -- Test opm on the Freudenstein and Roth function

genrose_test -- Test opm on a variant of the generalized Rosenbrock function
            ?? seems to be a problem with the dotarg passing

hessapptest.R -- to show how to change the source of the Hessian and
           Jacobian approximations used e.g., pracma vs numDeriv

hessian-used.R -- Illustration of use of hessian in optimization with
            several methods from optimx.

HessQuality.R -- tests on funconstrain set to check symmetry and other
           characteristics of the returned Hessian

hobbs.R -- Various applications of methods to the Hobbs Weed problem
           from Nash, Compact numerical methods for computers, 
           Adam Hilger:Bristol, 1979 and Second Edition, IOP: Bristol,
           1990. This is a problem that has some nasty features yet
           is very simple to state.

jonesrun.R -- tests using a function from Owen Jones, Robert Maillardet, 
           and Andrew Robinson, Scientific Programming and Simulation Using 
           R, Second Edition, Chapman & Hall/CRC, Boca Raton Page 230

maxtestJN -- Test opm with maximizing rather than minimizing

ncgtests.R -- various tests of the revised R conjugate gradients minimizer
	   ncg by J C Nash. Note that ncg() is in ongoing development.

onepar_test -- Test opm with a function of one parameter

optimrgrapprox.R -- Examples of calls to gradient approximations, showing
           failure if non-existent method specified

ox        -- Interactive examples using the legacy optimx() function.

poissmix_test -- Test opm on a Poisson mixture problem

rosbkext_test --    Test opm on the Rosenbrock function (extended)

rosenbrock.R -- diverse examples using the genrose extended Rosenbrock
           test function. genrose and rosbkext are DIFFERENT!

sc2_test    -- Tests with the relatively simple SC2 function. It would be
           nice to have a reference to the origin of this relatively
           simple but variable sized problem.

scalechk-Ex.R -- examples of calls to scalechk() function for checking the
           parameter scaling in optimization problems

simplefuntst.R -- Extremely simple example used to show syntax of calls.

snewtonbtest.R -- bounded tests using snewton methods

specctrlhobbs.R -- tests passing of special controls to bobyqa using Hobbs.

ssqbtest.R -- a simple sum of squares test showing differences between
              opm() and optimx()

trig1507.R -- More' Garbow and Hillstrom, 1981, problem no. 26 tests.
           Note that this function is also part of the funconstrain
           set and fhop.R will use it also.

trystarttests.R -- testing the starttests option in older optimx()
	   function (deprecated).

valley_test -- Test opm on the valley function. Again it would be nice to
           have a good reference for this function. The example shows the
           effectiveness of the ncg method for large numbers of parameters.

vmmix_test  -- Test opm on the mixture function

woodtest.R -- Some tests using the well-known 4-parameter wood function.

