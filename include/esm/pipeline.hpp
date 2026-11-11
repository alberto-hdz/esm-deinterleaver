#pragma once

#include <cstddef>
#include <vector>

#include "esm/clusterer.hpp"
#include "esm/pdw.hpp"
#include "esm/pri_analyzer.hpp"

namespace esm {

// One emitter the pipeline thinks it found.
struct Emitter {
    int id = 0;
    int cluster = -1;         // which parameter cluster it came out of
    double pri_us = 0.0;      // 0 means no PRI was found (pulses grouped by parameters only)
    PriType pri_type = PriType::Constant;
    double jitter = 0.0;              // +/- fraction of PRI, jitter only
    std::vector<double> stagger_us;   // intervals in one frame, stagger only (pri_us = frame)
    double freq_mhz = 0.0;
    double pw_us = 0.0;
    double aoa_deg = 0.0;
    std::vector<std::size_t> pulses;  // indices into the input PDW vector
};

struct PipelineResult {
    std::vector<Emitter> emitters;
    std::vector<int> pulse_label;  // emitter id per input pulse, -1 = unassigned
    int clusters = 0;
};

// Cluster on RF/AOA/PW, then run PRI deinterleaving inside each cluster.
// Pulses left over in a cluster after PRI analysis are kept as one more emitter
// with pri_us = 0 if there are enough of them, since they still matched on RF/AOA/PW.
PipelineResult run_pipeline(const std::vector<Pdw>& pdws, const ClusterParams& cp,
                            const PriParams& pp);

}  // namespace esm
