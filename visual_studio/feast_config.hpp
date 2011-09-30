#pragma once
#ifndef VISUAL_STUDIO_FEAST_CONFIG_HPP
#define VISUAL_STUDIO_FEAST_CONFIG_HPP 1

/**
 * \file
 * \brief FEAST configuration file for Visual Studio projects.
 *
 * \author Peter Zajac
 */

// disable the context stack - in VS we can use a *real* debugger for that...
#define FEAST_NO_CONTEXT 1

#endif // VISUAL_STUDIO_FEAST_CONFIG_HPP
