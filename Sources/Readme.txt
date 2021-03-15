*********************************************************************************
* USER NOTES
* - Sudip Seal, ORNL, April 26, 2013
*********************************************************************************
ENVIRONMENT VARIABLES

Environment variables that can be set by users are :
1. PARSOLVER (Default:TRUE or NSCALED): It can be set to FALSE (or ScaLAPACK)
using, for example, the following command at bash prompt:
  export PARSOLVER=FALSE or 
  export PARSOLVER=ScaLAPACK

2. PARFUNCTISL (DEFAULT:TRUE): This variable is automatically set FALSE when
PARSOLVER=FALSE/ScaLAPACK. When PARSOLVER=TRUE/NSCALED, SIESTA will execute
faster with PARFUNCTISL=TRUE than with PARFUNCTISL=FALSE for most cases. 

Several other environment variables for code development purposes may be found
in nscalingtools.f90 but they should not be changed by the user. Modifying these
variables maye result in indeterminate code behavior in terms of correctness as
well as runtime performance.
*********************************************************************************
