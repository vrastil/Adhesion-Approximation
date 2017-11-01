/**
 * @file:	emu.h
 * @brief:	Cosmic Emu : https://github.com/lanl/CosmicEmu
 */

#pragma once
#include "core.h"

const int nmode = 351;
void emu(double *xstar, double *ystar);
Data_x_y<double> init_emu(const Sim_Param &sim, double z);