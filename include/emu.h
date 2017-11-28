/**
 * @file:	emu.h
 * @brief:	Cosmic Emu : https://github.com/lanl/CosmicEmu
 */

#pragma once
#include "core.h"
namespace emu
{
    extern const int nmode;
    void emu(double *xstar, double *ystar);
    Data_Vec<double, 2> init_emu(const Sim_Param &sim, double z);
}
