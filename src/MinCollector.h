#ifndef KALLISTO_MINCOLLECTOR_H
#define KALLISTO_MINCOLLECTOR_H

#include "common.h"
#include <iostream>
#include <vector>

template <typename Index>
struct MinCollector {
	
MinCollector(Index &ind, const ProgramOptions& opt) : index(ind), counts(index.ecmap.size(), 0) {}


	
	void collect(const std::vector<int>& v) {
		if (v.empty()) {
			return;
		}

		std::vector<int> u = index.ecmap[v[0]];
		for (int i = 0; i < v.size(); i++) {
			u = index.intersect(v[i],u);
			if (u.empty()) {
				break;
			}
		}
		// if u is empty do nothing
		if (u.empty()) {
			return;
		}
		
		auto search = index.ecmapinv.find(u);
		if (search != index.ecmapinv.end()) {
			// ec class already exists, update count
			++counts[search->second];
		} else {
			// new ec class, update the index and count
			auto necs = counts.size();
			index.ecmap.insert({necs,u});
			index.ecmapinv.insert({u,necs});
			counts.push_back(1);
		}
	}
	
	void write(std::ostream& o) {
		for (int id = 0; id < counts.size(); id++) {
			o << id << "\t" << counts[id] << "\n";
		}
	}

	Index &index;
	std::vector<int> counts;
	
};

#endif // KALLISTO_MINCOLLECTOR_H
