/**
 * @file:	emu.h
 * @brief:	Cosmic Emu : https://github.com/lanl/CosmicEmu
 */

#pragma once
#include "core.h"
namespace emu
{
    extern const int nmode;
    void emu(FTYPE *xstar, FTYPE *ystar);
    Data_Vec<FTYPE, 2> init_emu(const Sim_Param &sim, FTYPE z);
}
