/**
 * @file:	emu.h
 * @brief:	Cosmic Emu : https://github.com/lanl/CosmicEmu
 */

#pragma once
#include "params.hpp"
#include "precision.hpp"
#include "templates/class_data_vec.hpp"

namespace emu
{
    extern const int nmode;
    void emu(FTYPE_t *xstar, FTYPE_t *ystar);
    Data_Vec<FTYPE_t, 2> init_emu(const Sim_Param &sim, FTYPE_t z);
}
