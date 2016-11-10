#!/bin/bash

DIR="/cygdrive/e/Dropbox/Argonne/Adhesion_Approximation/output/datastar"
DIR_anl="/home/ac.vrastil/Adhesion_Approximation/output/run"

scp anl:$DIR_anl/* $DIR/
