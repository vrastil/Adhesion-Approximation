/**
 * @file:	stdafx.h
 * @brief:	system include files and for project-specific include files that are used frequently but are changed infrequently
 */

#pragma once

/**********************
 * STANDARD LIBRARIES *
 **********************/

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <random>
#include <omp.h>
#include <sstream>
#include <vector>

/*******************
 * OTHER LIBRARIES *
 *******************/

#include <ccl.h>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include <fftw3.h>

#include <gsl/gsl_deriv.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_spline.h>

#include <json.hpp>