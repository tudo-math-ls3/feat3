/// \file FEAST user configuration header for Visual Studio 2010
/// \author Peter Zajac

/* ********************************************************************************************* */
// GENERAL INFORMATION
/* ********************************************************************************************* */
// This file is the user-configuration header template for the hand-made Visual Studio 2010
// FEAST project files.

// WARNING:
// If the file you are currently reading resides in the '[FEAST]/kernel/vc10' directory, then do
// *NOT* modify it. Instead, copy this file into your FEAST root-directory and modify that copy;
// Visual Studio will then automatically select your modified file instead of the original one.


/* ********************************************************************************************* */
// GENERAL USER-DEFINED SETTINGS
/* ********************************************************************************************* */

// disable the context stack - in VS we can use a *real* debugger for that...
#define FEAST_NO_CONTEXT 1


/* ********************************************************************************************* */
// INTEL MKL SUPPORT
/* ********************************************************************************************* */

// To enable Intel MKL support, uncomment the following macro-definition

// #define FEAST_BACKENDS_MKL 1


// Furthermore, define the FEAST_VC10_MKL_ROOT_DIR macro to your MKL installation directory
// without quotations marks and using forwards slashes, e.g.

// #define FEAST_VC10_MKL_ROOT_DIR C:/Program Files (x86)/Intel/mkl

/* END-OF-FILE */
