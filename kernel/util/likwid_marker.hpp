// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#ifndef FEAT_LIKWID_MARKER_HPP
#define FEAT_LIKWID_MARKER_HPP 1


/** \file
 *  \author Maximilian Esser
 */


// This block enables to compile the code with and without the likwid header.
// Since these are pre-compile macros, these are "ignored" if the code is build without
// the likwid-library
#ifdef LIKWID_PERFMON
#include <likwid-marker.h>
#else
// cpu markers
#define LIKWID_MARKER_INIT
#define LIKWID_MARKER_THREADINIT
#define LIKWID_MARKER_SWITCH
#define LIKWID_MARKER_REGISTER(regionTag)
#define LIKWID_MARKER_START(regionTag)
#define LIKWID_MARKER_STOP(regionTag)
#define LIKWID_MARKER_CLOSE
#define LIKWID_MARKER_GET(regionTag, nevents, events, time, count)
#endif

// special case if we did not activate cuda support, this logic requires
// that LIKWID_NVMON is only defined if LIKWID_PERFMON is also defined
#ifndef LIKWID_NVMARKER_INIT
//gpu markers
#define LIKWID_NVMARKER_INIT
#define LIKWID_NVMARKER_SWITCH
#define LIKWID_NVMARKER_REGISTER(regionTag)
#define LIKWID_NVMARKER_START(regionTag)
#define LIKWID_NVMARKER_STOP(regionTag)
#define LIKWID_NVMARKER_RESET(regionTag)
#define LIKWID_NVMARKER_CLOSE
#define LIKWID_NVMARKER_GET(name, ngpu, nevents, eventlist, time, count)
#endif

//these are simple aliases which should only be called once per programm, i.e. in FEAT runtime guard
/// Init the marker api
#define FEAT_MARKER_INIT LIKWID_MARKER_INIT
/// Only needed in non pthread based thread creation
#define FEAT_MARKER_THREADINIT LIKWID_MARKER_THREADINIT
/// Switch to the next performance group
#define FEAT_MARKER_SWITCH LIKWID_MARKER_SWITCH
/// Finalize marker api, writes data to file, only call once
#define FEAT_MARKER_CLOSE LIKWID_MARKER_CLOSE


// Now we define three different aliases of these markers, which essentially
// do the same thing, but can be activated seperatly...

#ifdef FEAT_KERNEL_MARKER_ACTIVATED

  #define FEAT_KERNEL_MARKER_REGISTER(regionTag) LIKWID_MARKER_REGISTER(regionTag)
  #define FEAT_KERNEL_MARKER_START(regionTag) LIKWID_MARKER_START(regionTag)
  #define FEAT_KERNEL_MARKER_STOP(regionTag) LIKWID_MARKER_STOP(regionTag)
  #define FEAT_KERNEL_MARKER_GET(regionTag, nevents, events, time, count) LIKWID_MARKER_GET(regionTag, nevents, events, time, count)
#else
  #define FEAT_KERNEL_MARKER_REGISTER(regionTag)
  #define FEAT_KERNEL_MARKER_START(regionTag)
  #define FEAT_KERNEL_MARKER_STOP(regionTag)
  #define FEAT_KERNEL_MARKER_GET(regionTag, nevents, events, time, count)

#endif //FEAT_KERNEL_MARKER_ACTIVATED

#ifdef FEAT_APPLICATION_MARKER_ACTIVATED
  #define FEAT_APPLICATION_MARKER_REGISTER(regionTag) LIKWID_MARKER_REGISTER(regionTag)
  #define FEAT_APPLICATION_MARKER_START(regionTag) LIKWID_MARKER_START(regionTag)
  #define FEAT_APPLICATION_MARKER_STOP(regionTag) LIKWID_MARKER_STOP(regionTag)
  #define FEAT_APPLICATION_MARKER_GET(regionTag, nevents, events, time, count) LIKWID_MARKER_GET(regionTag, nevents, events, time, count)
#else
  #define FEAT_APPLICATION_MARKER_REGISTER(regionTag)
  #define FEAT_APPLICATION_MARKER_START(regionTag)
  #define FEAT_APPLICATION_MARKER_STOP(regionTag)
  #define FEAT_APPLICATION_MARKER_GET(regionTag, nevents, events, time, count)
#endif //FEAT_APPLICATION_MARKER_ACTIVATED

#ifdef FEAT_SPECIAL_MARKER_ACTIVATED
  #define FEAT_SPECIAL_MARKER_REGISTER(regionTag) LIKWID_MARKER_REGISTER(regionTag)
  #define FEAT_SPECIAL_MARKER_START(regionTag) LIKWID_MARKER_START(regionTag)
  #define FEAT_SPECIAL_MARKER_STOP(regionTag) LIKWID_MARKER_STOP(regionTag)
  #define FEAT_SPECIAL_MARKER_GET(regionTag, nevents, events, time, count) LIKWID_MARKER_GET(regionTag, nevents, events, time, count)
#else
  #define FEAT_SPECIAL_MARKER_REGISTER(regionTag)
  #define FEAT_SPECIAL_MARKER_START(regionTag)
  #define FEAT_SPECIAL_MARKER_STOP(regionTag)
  #define FEAT_SPECIAL_MARKER_GET(regionTag, nevents, events, time, count)
#endif //FEAT_KERNEL_MARKER_ACTIVATED

// FEAT CUDA marker

//these are simple aliases which should only be called once per programm, i.e. in FEAT runtime guard
/// Init the marker api
#define FEAT_NVMARKER_INIT LIKWID_NVMARKER_INIT
/// Only needed in non pthread based thread creation
#define FEAT_NVMARKER_THREADINIT LIKWID_NVMARKER_THREADINIT
/// Switch to the next performance group
#define FEAT_NVMARKER_SWITCH LIKWID_NVMARKER_SWITCH
/// Finalize marker api, writes data to file, only call once
#define FEAT_NVMARKER_CLOSE LIKWID_NVMARKER_CLOSE

#define FEAT_NVMARKER_REGISTER(regionTag) LIKWID_NVMARKER_REGISTER(regionTag)
#define FEAT_NVMARKER_START(regionTag) LIKWID_NVMARKER_START(regionTag)
#define FEAT_NVMARKER_STOP(regionTag) LIKWID_NVMARKER_STOP(regionTag)
#define FEAT_NVMARKER_RESET(regionTag) LIKWID_NVMARKER_RESET(regionTag)
#define FEAT_NVMARKER_GET(name, ngpu, nevents, eventlist, time, count) LIKWID_NVMARKER_GET(name, ngpu, nevents, eventlist, time, count)

#endif //FEAT_LIKWID_MARKER_HPP