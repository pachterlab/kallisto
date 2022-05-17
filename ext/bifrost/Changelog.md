### Changelog

API only.

* **04-28-2022**
	* Color files generated prior to version 1.0.6.2 are **not** compatible with version 1.0.6.2 and onward.
	* `CompactedDBG::simplify()` and `ColoredCDBG::simplify()` now return true even if no simplification was performed ("null-simplification" in case all input parameters are set to false). The goal is to only return false if the graph is invalid or in case of unexpected behavior. 
* **08-29-2018**
	* `UnitigColors::const_iterator` only considers now the k-mer positions of the unitig mapping provided in the `UnitigMap`/`UnitigColorMap` parameter of `UnitigColors::begin()`.
* **08-28-2018**
	* Functions `CDBG_Data_t::serialize()` and `CCDBG_Data_t::serialize()` have now one parameter which is a `const_UnitigMap`/`const_UnitigColorMap` reference representing the (reference) unitig to which the data are associated to.
* **08-23-2018**
	* Function `UnitigMap::toString()` was ambiguous as if it would generate the string of the mapped sequence or the string of the reference unitig used in the mapping. It has been replaced by two functions: `UnitigMap::mappedSequenceToString()` and `UnitigMap::referenceUnitigToString()`.
* **08-07-2018**
	* Add de Bruijn graphs merging functions (`CompactedDBG::merge()` and `ColoredCDBG::merge()`) and addition assignment operators (`CompactedDBG::operator+=()` and `ColoredCDBG::operator+=()`, same as `merge()` but uses only one thread).
	* Add de Bruijn graphs comparison functions `CompactedDBG::operator==()`, `CompactedDBG::operator!=()`, `ColoredCDBG::operator==()` and `ColoredCDBG::operator!=()`.
	* Delete `CompactedDBG::empty()` and `ColoredCDBG::empty()` to be consistent with STD containers (those functions were emptying the graph of its content while `empty()` of STD containers returns whether the container is empty). Now, to empty the graph, use `CompactedDBG::clear()` and `ColoredCDBG::clear()`.
    * Major changes in the abstract class `CDBG_Data_t` and `CCDBG_Data_t`:
    	* All the functions are now **non**-static.
    	* Function `join()` is renamed `concat()` and works a bit differently (have a look at the doc). Quickly, `join()` was concatenating two unitigs A and B such that the result was A=AB and B was deleted from the graph. Now, `concat()` deletes A and B from the graph and adds a new unitig C=AB.
    	* Function `sub()` is renamed `extract()`.
    	* Add the functions `merge()` and `clear()` which **must** be overloaded too in the derived class of `CDBG_Data_t` and `CCDBG_Data_t`.
