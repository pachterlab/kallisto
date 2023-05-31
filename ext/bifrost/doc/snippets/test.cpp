#include <iostream>
#include <stack>
#include <queue>
#include <unordered_set>

#include <bifrost/ColoredCDBG.hpp>

using namespace std;

// Boolean class indicating if its associated unitig was "visited", "seen" or none of those two.
// The class inherits from CCDBG_Data_t<MyBool> to be used with ColoredCDBG and from
// CDBG_Data_t<MyBool> to be used with CompactedDBG.

class MyBool : public CCDBG_Data_t<MyBool>, CDBG_Data_t<MyBool> {

    public:

        MyBool() : b(NOT_VISITED_SEEN) {} // Initiate the boolean to "not visited"

        // Clear method for CompactedDBG
        void clear(const UnitigMap<MyBool>& um_dest){

        	set_not_seen_visited(); // Set the new unitig to "not seen nor visited"
        }

         // Clear method for CompactedDBG
        void clear(const UnitigColorMap<MyBool>& um_dest){

        	set_not_seen_visited(); // Set the new unitig to "not seen nor visited"
        }

        // Concatenation method for ColoredCDBG
        void concat(const UnitigColorMap<MyBool>& um_dest, const UnitigColorMap<MyBool>& um_src){

            // When concatenating the reference unitig of um_src to the reference unitig of um_dest,
            // we set the boolean of the new unitig to "not visited" because it will be a new unitig.

            set_not_seen_visited(); // Set the new unitig to "not seen nor visited"
        }

        // Concatenation method for ColoredCDBG
        void concat(const UnitigMap<MyBool>& um_dest, const UnitigMap<MyBool>& um_src){

            // When concatenating the reference unitig of um_src to the reference unitig of um_dest,
            // we set the boolean of the new unitig to "not visited" because it will be a new unitig.

            set_not_seen_visited(); // Set the new unitig to "not seen nor visited"
        }

        // Extraction method for ColoredCDBG
        void extract(const UnitigColors* uc_dest, const UnitigColorMap<MyBool>& um_src, const bool last_extraction) {

            // This function creates a new unitig which is a sub-unitig from the reference unitig of um_src.
            // The new unitig created is set to "not seen nor visited" as a measure of precaution (it is already
            // initiated by default to "not visited" in the constructor)

            set_not_seen_visited(); // Set the new unitig to "not seen nor visited"
        }

        // Extraction method for ColoredCDBG
        void extract(const UnitigMap<MyBool>& um_src, bool last_extraction) {

            // This function creates a new unitig which is a sub-unitig from the reference unitig of um_src.
            // The new unitig created is set to "not seen nor visited" as a measure of precaution (it is already
            // initiated by default to "not visited" in the constructor)

            set_not_seen_visited(); // Set the new unitig to "not seen nor visited"
        }

        // Methods MyBool::merge() and MyBool

        string toString() const {

            if (is_visited()) return string("visited");
            if (is_seen()) return string("seen");

            return string("Not seen nor visited");
        }

        inline void set_visited() { b = VISITED; } // Set the boolean to "visited"
        inline void set_seen() { b = SEEN; } // Set the boolean to "seen"
        inline void set_not_seen_visited() { b = NOT_VISITED_SEEN; } // Set the boolean to "not seen and not visited"

        inline bool is_visited() const { return (b == VISITED); } // return if the boolean is "visited"
        inline bool is_not_visited() const { return !is_visited(); } // return if the boolean is "not visited"

        inline bool is_seen() const { return (b == SEEN); } // return if the boolean is "seen"
        inline bool is_not_seen() const { return !is_seen(); } // return if the boolean is "not seen"

    private:

        const static uint8_t NOT_VISITED_SEEN = 0x0;
        const static uint8_t VISITED = 0x1;
        const static uint8_t SEEN = 0x2;

        uint8_t b;
};

void clearMarking(const unordered_set<UnitigColorMap<MyBool>, UnitigMapHash<DataAccessor<MyBool>, DataStorage<MyBool>>>& set_km_seen){

    for (const auto& ucm : set_km_seen){

        DataAccessor<MyBool>* da_ucm = ucm.getData();
        MyBool* data_ucm = da_ucm->getData(ucm);

        data_ucm->clear(ucm);
    }
}

void clearMarking(ColoredCDBG<MyBool>& ccdbg){

    for (const auto& unitig : ccdbg){

        DataAccessor<MyBool>* da_ucm = unitig.getData();
        MyBool* data_ucm = da_ucm->getData(unitig);

        data_ucm->clear(unitig);
    }
}

void BFS_Recursive(const UnitigColorMap<MyBool>& ucm){

    DataAccessor<MyBool>* da = ucm.getData(); // Get DataAccessor from unitig
    MyBool* data = da->getData(ucm); // Get boolean from DataAccessor

    if (data->is_not_visited() && data->is_not_seen()){

	    data->set_visited(); // Set boolean to indicate unitig was visited

	    for (auto& successor : ucm.getSuccessors()){ // Iterate over successors of a unitig

	        DataAccessor<MyBool>* da_succ = successor.getData(); // Get DataAccessor from unitig successor
	        MyBool* data_succ = da_succ->getData(successor); // Get boolean from DataAccessor

	        if (data_succ->is_not_visited()) data_succ->set_seen(); // Set boolean to indicate unitig was seen (visited but not traversed yet)
	    }

	    for (auto& predecessor : ucm.getPredecessors()){ // Iterate over predecessors of a unitig

	        DataAccessor<MyBool>* da_pred = predecessor.getData(); // Get DataAccessor from unitig predecessor
	        MyBool* data_pred = da_pred->getData(predecessor); // Get boolean from DataAccessor

	        if (data_pred->is_not_visited()) data_pred->set_seen(); // Set boolean to indicate unitig was visited (visited but not traversed yet)
	    }

	    // Traverse successors
	    for (auto& successor : ucm.getSuccessors()){

	        DataAccessor<MyBool>* da_succ = successor.getData(); // Get DataAccessor from unitig successor
	        MyBool* data_succ = da_succ->getData(successor); // Get boolean from DataAccessor

	    	if (data_succ->is_seen()){

	    		data_succ->set_not_seen_visited();

	    		BFS_Recursive(successor);
	    	}
	    }

	    // Traverse predecessors
	    for (auto& predecessor : ucm.getSuccessors()){

	        DataAccessor<MyBool>* da_pred = predecessor.getData(); // Get DataAccessor from unitig predecessor
	        MyBool* data_pred = da_pred->getData(predecessor); // Get boolean from DataAccessor

	    	if (data_pred->is_seen()){

	    		data_pred->set_not_seen_visited();

	    		BFS_Recursive(predecessor);
	    	}
	    }
	}
}

void BFS_Iterative(const UnitigColorMap<MyBool>& ucm){

    queue<UnitigColorMap<MyBool>> q; // Create queue of unitig to traverse
    UnitigColorMap<MyBool> ucm_tmp(ucm); // Create a non-const local copy of unitig given in parameter

    DataAccessor<MyBool>* da = ucm_tmp.getData(); // Get DataAccessor from unitig
    MyBool* data = da->getData(ucm_tmp); // Get boolean from DataAccessor

    data->set_visited(); // Set boolean to indicate unitig was visited

    q.push(ucm_tmp); // Push unitig to traverse on the stack

    while (!q.empty()){ // While they are unitigs to traverse in the stack

        ucm_tmp = q.front(); // Get unitig at the front of the queue

        q.pop(); // Delete unitig at the front of the queue

        for (auto& successor : ucm_tmp.getSuccessors()){ // Traverse successors

            DataAccessor<MyBool>* da_succ = successor.getData(); // Get DataAccessor from successor
            MyBool* data_succ = da_succ->getData(successor); // Get boolean from DataAccessor

            if (data_succ->is_not_visited()){ // If boolean indicates the successor was not visited

                data_succ->set_visited(); // Set boolean to indicate successor was visited

                q.push(successor); // Traverse neighbors of successor
            }
        }

        // Traverse predecessors
        for (auto& predecessor : ucm_tmp.getPredecessors()){

            DataAccessor<MyBool>* da_pred = predecessor.getData(); // Get DataAccessor from predecessor
            MyBool* data_pred = da_pred->getData(predecessor); // Get boolean from DataAccessor

            if (data_pred->is_not_visited()){ // If boolean indicates the predecessor was not visited

                data_pred->set_visited(); // Set boolean to indicate predecessor was visited

                q.push(predecessor); // Traverse neighbors of predecessor
            }
        }
    }
}

void DFS_Recursive(const UnitigColorMap<MyBool>& ucm){

    DataAccessor<MyBool>* da = ucm.getData(); // Get DataAccessor from unitig
    MyBool* data = da->getData(ucm); // Get boolean from DataAccessor

    if (data->is_not_visited()){

	    data->set_visited(); // Set boolean to indicate unitig was visited

	    for (auto& successor : ucm.getSuccessors()) DFS_Recursive(successor); // Traverse neighbors of successor
	    for (auto& predecessor : ucm.getPredecessors()) DFS_Recursive(predecessor); // Traverse neighbors of predecessor
	}
}

void DFS_Iterative(const UnitigColorMap<MyBool>& ucm){

    stack<UnitigColorMap<MyBool>> stck; // Create stack of unitig to traverse
    UnitigColorMap<MyBool> ucm_tmp(ucm); // Create a non-const local copy of unitig given in parameter

    stck.push(ucm_tmp); // Push first unitig to traverse on the stack

    while (!stck.empty()){ // While they are unitigs to traverse in the stack

        ucm_tmp = stck.top(); // Get the unitig on top of the stack

        stck.pop(); // Delete unitig on the top of the stack

        DataAccessor<MyBool>* da = ucm_tmp.getData(); // Get DataAccessor from unitig
        MyBool* data = da->getData(ucm_tmp); // Get boolean from DataAccessor

        if (data->is_not_visited()){ // If boolean indicates the unitig was not visited

            data->set_visited(); // Set boolean to indicate unitig was visited

            // Add successors to stack of unitigs to traverse
            for (auto& successor : ucm_tmp.getSuccessors()) stck.push(successor);
            // Add predecessors to stack of unitigs to traverse
            for (auto& predecessor : ucm_tmp.getPredecessors()) stck.push(predecessor);
        }
    }
}

// - Parameter "iterative" indicates if you want to traverse the graph in a iterative (true)
// or recursive (false) way. Recursive version is limited to small components because of
// the stack size (recursive calls are piling up variables on your stack) which is not
// the case of the iterative version. Default is iterative.
void traverse(ColoredCDBG<MyBool>& ccdbg, const bool DFS = true, const bool iterative = true){

    for (auto& unitig : ccdbg){ // Iterate over unitigs of a colored de Bruijn graph

        if (iterative) DFS ? DFS_Iterative(unitig) : BFS_Iterative(unitig); // Traverse neighbors of unitig in an iterative manner
        else DFS ? DFS_Recursive(unitig) : BFS_Recursive(unitig); // Traverse neighbors of unitig in a recursive manner
    }

    clearMarking(ccdbg);
}

size_t getNbConnectedComponent(ColoredCDBG<MyBool>& ccdbg){

    size_t nb_cc = 0; // Number of connected components

    for (auto& unitig : ccdbg){ // Iterate over unitigs of a colored de Bruijn graph

        DataAccessor<MyBool>* da = unitig.getData(); // Get DataAccessor from unitig
        MyBool* data = da->getData(unitig); // Get boolean from DataAccessor

        if (data->is_not_visited()){ // If boolean indicates the unitig was not visited

            ++nb_cc; // It's a new connected components

            DFS_Iterative(unitig); // Traverse neighbors of unitig in DFS
        }
    }

    clearMarking(ccdbg);

    return nb_cc;
}

pair<UnitigColorMap<MyBool>, UnitigColorMap<MyBool>> extractSuperBubble(const UnitigColorMap<MyBool>& s){

    vector<UnitigColorMap<MyBool>> vertices_visit;

    pair<UnitigColorMap<MyBool>, UnitigColorMap<MyBool>> p;

    unordered_set<UnitigColorMap<MyBool>, UnitigMapHash<DataAccessor<MyBool>, DataStorage<MyBool>>> set_km_seen;

    UnitigColorMap<MyBool> v(s);

    vertices_visit.push_back(v);

    while (!vertices_visit.empty()){ // While there are vertices to visit

        v = vertices_visit.back(); // Pick arbitrary vertex in the set

        vertices_visit.pop_back(); // Delete that vertex from the set

        DataAccessor<MyBool>* da = v.getData(); // Get DataAccessor from unitig
        MyBool* data = da->getData(v); // Get boolean from DataAccessor

        //cout << v.getUnitigHead().toString() << " - " << data->toString() << endl;

        data->set_visited(); // Set boolean to indicate unitig was visited
        set_km_seen.insert(v);

        if (!v.getSuccessors().hasSuccessors()){ // Unitig is no successor, it is a tip

            clearMarking(set_km_seen);
            return p;
        }

        for (auto& u : v.getSuccessors()){

            if (u == s){ // Cycle including s

                clearMarking(set_km_seen);
                return p;
            }

            DataAccessor<MyBool>* da_u = u.getData(); // Get DataAccessor from successor
            MyBool* bool_u = da_u->getData(u); // Get boolean from DataAccessor

            if (bool_u->is_not_visited()){

	            bool_u->set_seen(); // Set boolean to indicate unitig was seen
	            set_km_seen.insert(u);

	            bool all_predecessor_visited = true;

	            for (const auto& predecessor : u.getPredecessors()){

	                const DataAccessor<MyBool>* da_pred = predecessor.getData();
	                const MyBool* data_pred = da_pred->getData(predecessor);

	                if (data_pred->is_not_visited()){

	                    all_predecessor_visited = false;
	                    break;
	                }
	            }

	            if (all_predecessor_visited) vertices_visit.push_back(u);
        	}
        	else {

                clearMarking(set_km_seen);
                return p;
        	}
        }

        if (vertices_visit.size() == 1){

            bool not_seen = true;

            for (const auto& cucm : set_km_seen){

                if (cucm != vertices_visit[0]){

                    const DataAccessor<MyBool>* da_cucm = cucm.getData();
                    const MyBool* data_cucm = da_cucm->getData(cucm);

                    if (data_cucm->is_seen()){

                        not_seen = false;
                        break;
                    }
                }
            }

            if (not_seen){

                for (const auto& successor : vertices_visit[0].getSuccessors()){

                    if (successor == s){ // cycle

                        clearMarking(set_km_seen);
                        return p;
                    }
                }

                p.first = s;
                p.second = vertices_visit[0];

                clearMarking(set_km_seen);
                return p;
            }
        }
    }

    clearMarking(set_km_seen);
    return p;
}

int main(int argc, char *argv[])
{

    if (argc != 1){

    	ColoredCDBG<MyBool> ccdbg;

        // Read input filenames
        const string filename_prefix(argv[1]);

        cout << "=== Reading graph ===" << endl;

        if (ccdbg.read(filename_prefix, 4)){

	        cout << "=== Traversing graph: Depth First Search ===" << endl;

	        traverse(ccdbg, true);

	        cout << "=== Traversing graph: Breadth First Search ===" << endl;

	        traverse(ccdbg, false);

	        cout << "=== Computing number of connected components ===" << endl;

	        const size_t nb_cc = getNbConnectedComponent(ccdbg);

	        cout << nb_cc << " connected components found" << endl;

	        cout << "=== Computing super bubbles ===" << endl;

	        size_t nb_super_bubble = 0;
	        size_t nb_unitig_processed = 0;

	        for (const auto& unitig : ccdbg){

	            const pair<UnitigColorMap<MyBool>, UnitigColorMap<MyBool>> p = extractSuperBubble(unitig);

	            if (!p.first.isEmpty && !p.second.isEmpty) ++nb_super_bubble;

	            ++nb_unitig_processed;

	            if (nb_unitig_processed % 100000 == 0){

	                cout << "Processed " << nb_unitig_processed << " unitigs: " << nb_super_bubble << " super bubbles so far" << endl;
	            }
	        }

	        cout << nb_super_bubble << " super bubbles found" << endl;
    	}
    	else cerr << "Invalid graph filename prefix provided" << endl;
    }
    else cerr << "No input graph filename prefix provided" << endl;
}
