#ifndef MINHASHITERATOR_H
#define MINHASHITERATOR_H

#include <stdint.h>
#include <deque>

//----
#include <iostream>
using namespace std;
//---

//using std::vector;

struct minHashResult {
	minHashResult() : hash((uint64_t) -1),pos(-1) {}
	minHashResult(const uint64_t h, const int p) : hash(h), pos(p) {}
	uint64_t hash;
	int pos;
};

template<class HF>
struct minHashResultIterator;


// TODO: return current position as well as minimum										
// templated?
template<class HF>
class minHashIterator {
public:
	minHashIterator(int _k, int _g, HF _h) : s(NULL), n(0), k(_k), g(_g), hf(_h), v(k-g+1), p(-1),invalid(true) {
		hf.setK(g);
  }
	
	minHashIterator(const char* _s, int _length, int _k, int _g, HF _h) : k(_k), g(_g), hf(_h), v(k-g+1), p(-1), invalid(true) {
		hf.setK(g);
		initString(_s,_length);
	}

  minHashIterator() : s(NULL),n(0), k(0), g(0),hf(HF(0)),invalid(true) {}

  void seed(HF &ohf) {
    hf = ohf;
    invalid=true;
  }
  
	void initString(const char* _s, int _length) {
		n = _length;
		s = _s;
		p = -1;// reset position
		invalid = false;
		v.clear();
		
		// advance to first position or set to invalid
		if (n < k || k < g) {
			invalid = true;
		} else {
			operator++();
		}

	}
		
	
	bool operator==(const minHashIterator& o) {
		if(invalid || o.invalid) {
			return invalid && o.invalid;
		} else {
			// same string, position and parameters, rest is deterministic
			return s==o.s && n==o.n && g==o.g && k==o.k; 
		}
	}
	
	bool operator!=(const minHashIterator& o) {return !this->operator==(o);}

	// invariant:
	//   v[0..sz] contains only ascending elements in s from [p,p+k-g+1)
	//
	//   there exists no a<b s.t. v[a] > v[b].
	minHashIterator& operator++() {
		//cerr << "operator++";
		if (invalid) {
			return *this;
		}

		++p; // advance to next k-mer
		if (p >= n-k+1 || s[p+k-1] == 0) {
			// out of bounds
			invalid = true;
			return *this;
		}
		//cerr << ", p = " << p << endl;

		if (p==0) {
			//cerr << "first " << endl;
			// for first position, find all min values
			hf.init(s);
			v.push_back(minHashResult(hf.hash(),0));
			//for (int i = 0; i < v.size(); i++) { cerr << "(" << v[i].hash << "," << v[i].pos << ")" << ", ";} cerr << "]" << endl;

			//cerr << "  " << hf.hash() << endl;
			
			for (int j = 1; j < k-g+1; j++) {
				// invariant: v[0..sz) contains ascending values in s[p..p+j)
				hf.update(s[j-1],s[j+g-1]);
				uint64_t h = hf.hash();
				//cerr << "  " << h << " - v = [";

        int t = ((int)v.size())-1; 
				while (t >= 0 && v[t].hash > h) { 
					// h is strictly smaller than last element, remove it
					v.pop_back();
					t--;
				}
				v.push_back(minHashResult(h,j));
				//for (int i = 0; i < v.size(); i++) { cerr << "(" << v[i].hash << "," << v[i].pos << ")" << ", ";} cerr << "]" << endl;
			}
		} else {
			//cerr << "update" << endl;
			if (v[0].pos < p) {
				v.pop_front(); // remove first element, fell outside of window
			}
			hf.update(s[p+k-g-1],s[p+k-1]);
			uint64_t h = hf.hash();
			//cerr << "  " << h << " - v = [";
      int t = ((int) v.size())-1;
			while (t >= 0 && v[t].hash > h) {
				v.pop_back();
				t--;
			}
			v.push_back(minHashResult(h,p+k-g));
			//for (int i = 0; i < v.size(); i++) { cerr << "(" << v[i].hash << "," << v[i].pos << ")" << ", ";} cerr << "]" << endl;
		}
    return *this;
	}
	
	minHashIterator operator++(int) {
		minHashIterator tmp(*this);
		operator++();
		return tmp;
	}

	
	minHashResultIterator<HF> operator*() const {return minHashResultIterator<HF>(this);}
	//minHashResultIterator<HF>* operator->() { return &(operator*());}

	const char *s;
	int n;
	int k;
	int g;
	HF hf;
	deque<minHashResult> v;
	int p;
	bool invalid;


	// private copy constructor
  minHashIterator(const minHashIterator& o) : s(o.s), n(o.n), k(o.k), g(o.g), hf(o.hf), v(o.v), p(o.p), invalid(o.invalid) {}
};


template <class HF>
struct minHashResultIterator {
  minHashResultIterator(const minHashIterator<HF> *p) : p(p), invalid(false), pos(0), p_pos(p->p), p_s(p->s) {}
  minHashResultIterator() : invalid(true), p_pos(-1), p_s(NULL) {}
  
  
  bool operator==(const minHashResultIterator& o) {
    if (o.invalid || invalid) {
      return o.invalid && invalid;
    } else {
      return p_pos == o.p_pos && p_s == o.p_s && pos==o.pos;
    }
  }
  
  bool operator!=(const minHashResultIterator& o) {
    return !this->operator==(o);
  }
  
  minHashResultIterator& operator++() {
    if (invalid) {
      return *this;
    }
    
    // check if parent iterator has moved
    if ((p_pos != p->p || p_s != p->s)
        || (pos>=p->v.size()-1) // or if we advance past the end
        || (p->v[pos+1].hash != p->v[pos].hash))
      // or if the next position doesn't share the hash value
    {
      // invalidate the iterator
      p_s = NULL;
      invalid=true;
      return *this;
    } else {
      pos++; // advance to next equal hash value
    }
    return *this;
  }
  
  
  minHashResultIterator operator++(int) {
    minHashResultIterator tmp(*this);
    operator++();
    return tmp;
  }
  
  const minHashResult& operator*() const {return p->v[pos];}
  const minHashResult* operator->() const {return &(p->v[pos]);}
  
  
  // pos points to a minHashIterator, all the values from p.v[0] to p.v[pos] have the
  // same (minimum) hash value or the this.invalid is true
  // at the time this was created p.s==p_s and p.p=p_pos
  const minHashIterator<HF> *p;
  bool invalid;
  int pos;
  const int p_pos;
  const char *p_s;
};





#endif // MINHASHITERATOR_H
