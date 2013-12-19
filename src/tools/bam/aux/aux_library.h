/**
* Copyright (C) 2013 Raúl Moreno Galdón
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef AUX_LIBRARY_H_
#define AUX_LIBRARY_H_

//System libs
#include <math.h>
#include <float.h>
#include <fenv.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>

//Bioinfo libs
#include "bioformats/bam/samtools/bam.h"
#include "bioformats/bam/bam_file.h"

//Common
#include "aux_common.h"

//Additional headers
#include "recalibrate/recal_config.h"

//Auxiliar headers
#include "aux_bam.h"
#include "aux_cigar.h"
#include "aux_math.h"
#include "aux_misc.h"
#include "aux_nucleotide.h"
#include "aux_quality.h"
#include "aux_vector.h"


#endif /* AUX_LIBRARY_H_ */
