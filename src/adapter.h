#ifdef _ADAPTER_H
#define _ADAPTER_H

#include <string.h>
#include <stdlib.h>

#include "bioformats/fastq/fastq_read.h"

//--------------------------------------------------------------------

void cut_adapter(char *adapter, int adapter_length, fastq_read_t *read);

//--------------------------------------------------------------------
//--------------------------------------------------------------------
#endif // _ADAPTER_H
