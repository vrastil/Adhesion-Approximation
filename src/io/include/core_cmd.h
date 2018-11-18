/**
 * @brief command line arguments manipulation
 * 
 * @file core_cmd.h
 * @author Michal Vrastil
 * @date 2018-07-11
 */

#pragma once

class Sim_Param;

/**
 * @brief parse command line arguments
 * 
 * Take command line arguments and store parameters in \p sim.
 * If config file is specified in arguments (or the default one is present),
 * parse also this file.
 * 
 * @param ac number of command line arguments
 * @param av command line arguments
 * @param sim where to store simulation parameters
 */
void handle_cmd_line(int ac, const char* const av[], Sim_Param& sim);