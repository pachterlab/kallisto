#pragma once

// MetalReadLoader is a thin alias for HostReadLoader defined in MetalPipeline.h.
// CPU gzip decompression via kseq/ProcessReads; result lives in shared MetalBuffers
// via DeviceKmerLoader::load().

#include "MetalPipeline.h"

// (No additional declarations needed beyond MetalPipeline.h)
