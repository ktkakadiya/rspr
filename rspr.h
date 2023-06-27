/*******************************************************************************
rspr.h

Calculate approximate and exact Subtree Prune and Regraft (rSPR)
distances and the associated maximum agreement forests (MAFs) between pairs
of rooted binary trees.
Supports arbitrary labels. See the
README for more information.

Copyright 2009-2014 Chris Whidden
whidden@cs.dal.ca
http://kiwi.cs.dal.ca/Software/RSPR
April 29, 2014
Version 1.2.2

This file is part of rspr.

rspr is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

rspr is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with rspr.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************/

#define RSPR
//#define DEBUG 1
//#define DEBUG_CONTRACTED 1
//#define DEBUG_APPROX 1
//#define DEBUG_CLUSTERS 1
//#define DEBUG_SYNC 1
//#define DEBUG_UNDO 1
//#define DEBUG_DEPTHS 1
//#define DEBUG_CASE_COUNTER 1
//#define MUL-T_PICK_LARGEST_GROUP 1
//#define DEBUG_LGT_EVENTS 1

#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <climits>
#include <vector>
#include <map>
#include <set>
#include <list>
#include <unordered_set>
#include <algorithm>
#include "Forest.h"
#include "ClusterForest.h"
#include "LCA.h"
#include "ClusterInstance.h"
#include "SiblingPair.h"
#include "UndoMachine.h"
#include "ClusterMergeForest.h"

using namespace std;

enum RELAXATION {STRICT, NEGATIVE_RELAXED, ALL_RELAXED};

const string whitespaces = " \t\f\v\n\r";

// note: not using undo
int rSPR_worse_3_mult_approx_hlpr(Forest *T1, Forest *T2, list<Node *> *singletons, list<Node *> *sibling_pairs, Forest **F1, Forest **F2, bool save_forests);
int rSPR_worse_3_mult_approx(Forest *T1, Forest *T2);
int rSPR_worse_3_mult_approx(Forest *T1, Forest *T2, bool sync);

int rSPR_branch_and_bound_mult_range(Forest *T1, Forest *T2, int start_k);
int rSPR_branch_and_bound_mult_range(Forest *T1, Forest *T2, int start_k, int end_k);
int rSPR_branch_and_bound_mult(Forest *T1, Forest *T2, int k);
int rSPR_branch_and_bound_mult_hlpr(Forest *T1, Forest *T2, int k, list<Node*> *sibling_groups, list<Node*> *singletons, Node *protected_stack, list<pair<Forest,Forest>> *AFs, int* num_ties);



int rSPR_3_approx_hlpr(Forest *T1, Forest *T2, list<Node *> *singletons,
		list<Node *> *sibling_pairs);
int rSPR_3_approx(Forest *T1, Forest *T2);
int rSPR_worse_3_approx_hlpr(Forest *T1, Forest *T2, list<Node *> *singletons, list<Node *> *sibling_pairs, Forest **F1, Forest **F2, bool save_forests);
int rSPR_worse_3_approx(Forest *T1, Forest *T2);
int rSPR_worse_3_approx(Forest *T1, Forest *T2, bool sync);
int rSPR_worse_3_approx(Node *subtree, Forest *T1, Forest *T2);
int rSPR_worse_3_approx(Node *subtree, Forest *T1, Forest *T2, bool sync);
int rSPR_worse_3_approx_binary_hlpr(Forest *T1, Forest *T2, list<Node *> *singletons, list<Node *> *sibling_pairs, Forest **F1, Forest **F2, bool save_forests);
int rSPR_worse_3_approx_binary(Forest *T1, Forest *T2, bool sync);
int rSPR_worse_3_approx_binary(Forest *T1, Forest *T2);
int rSPR_branch_and_bound(Forest *T1, Forest *T2);
int rSPR_branch_and_bound(Forest *T1, Forest *T2, int k);
int rSPR_branch_and_bound(Forest *T1, Forest *T2, int k, list<pair<Forest,Forest>> *extAFs);
int rSPR_branch_and_bound(Forest *T1, Forest *T2, int k,
		map<string, int> *label_map,
		map<int, string> *reverse_label_map,
		list<pair<Forest,Forest>> *extAFs);

int rSPR_branch_and_bound_range(Forest *T1, Forest *T2, int end_k);
int rSPR_branch_and_bound_range(Forest *T1, Forest *T2, int start_k,
		int end_k);
int rSPR_branch_and_bound_hlpr(Forest *T1, Forest *T2, int k,
		set<SiblingPair> *sibling_pairs, list<Node *> *singletons, bool cut_b_only,
		list<pair<Forest,Forest> > *AFs, list<Node *> *protected_stack,
		int *num_ties);
int rSPR_branch_and_bound_hlpr(Forest *T1, Forest *T2, int k,
		set<SiblingPair> *sibling_pairs, list<Node *> *singletons, bool cut_b_only,
		list<pair<Forest,Forest> > *AFs, list<Node *> *protected_stack,
		int *num_ties, Node *prev_T1_a, Node *prev_T1_c);
void rspr_branch_and_bound_cut_a_hlpr(Forest *T1, Forest *T2, int k,
	set<SiblingPair> *sibling_pairs, list<Node *> *singletons, list<pair<Forest,Forest> > *AFs,
	list<Node *> *protected_stack, int *num_ties, Node *T1_c, Node *T2_a, Node *T2_b, Node* T2_c,
	bool cut_a_only, bool cut_b_only, bool cut_c_only, int path_length, UndoMachine *um,
	Node *T2_ab, bool balanced, bool multi_b1, bool multi_b2, Node *T2_d, int &best_k);
void rspr_branch_and_bound_cut_c_hlpr(Forest *T1, Forest *T2, int k,
	set<SiblingPair> *sibling_pairs, list<Node *> *singletons, list<pair<Forest,Forest> > *AFs,
	list<Node *> *protected_stack, int *num_ties, Node *T1_a, Node *T2_a, Node *T2_b, Node* T2_c,
	bool cut_a_only, bool cut_b_only, bool cut_c_only, int path_length, UndoMachine *um,
	bool balanced, bool multi_b1, bool multi_b2, bool cut_a_or_merge_ac, bool cut_ab_only,
	Node *T2_d, int lca_depth, int &best_k);
int rSPR_total_approx_distance(Node *T1, vector<Node *> &gene_trees);
int rSPR_total_approx_distance(Node *T1, vector<Node *> &gene_trees,
		int threshold);
int rSPR_total_distance(Node *T1, vector<Node *> &gene_trees);
int rSPR_total_distance(Node *T1, vector<Node *> &gene_trees,
		vector<int> *original_scores);
void rSPR_pairwise_distance(Node *T1, vector<Node *> &gene_trees);
void rSPR_pairwise_distance(Node *T1, vector<Node *> &gene_trees, bool approx);
void rSPR_pairwise_distance(Node *T1, vector<Node *> &gene_trees, int start, int end);
void rSPR_pairwise_distance(Node *T1, vector<Node *> &gene_trees, int start, int end, bool approx);
void rSPR_pairwise_distance(Node *T1, vector<Node *> &gene_trees, int max_spr);
void rSPR_pairwise_distance(Node *T1, vector<Node *> &gene_trees, int max_spr, int start, int end);
void rSPR_pairwise_distance_unrooted(Node *T1, vector<Node *> &gene_trees);
void rSPR_pairwise_distance_unrooted(Node *T1, vector<Node *> &gene_trees, int start, int end);
void rSPR_pairwise_distance_unrooted(Node *T1, vector<Node *> &gene_trees, bool approx);
void rSPR_pairwise_distance_unrooted(Node *T1, vector<Node *> &gene_trees, int start, int end, bool approx);
int rf_total_distance(Node *T1, vector<Node *> &gene_trees);
int rf_total_distance_unrooted(Node *T1, vector<Node *> &gene_trees);
void rf_pairwise_distance(Node *T1, vector<Node *> &gene_trees);
void rf_pairwise_distance(Node *T1, vector<Node *> &gene_trees, int start, int end);
void rf_pairwise_distance_unrooted(Node *T1, vector<Node *> &gene_trees);
void rf_pairwise_distance_unrooted(Node *T1, vector<Node *> &gene_trees, int start, int end);
int rSPR_total_distance_unrooted(Node *T1, vector<Node *> &gene_trees, int threshold);
int rSPR_total_distance_unrooted(Node *T1, vector<Node *> &gene_trees, int threshold, vector<int> *original_scores);
int rSPR_branch_and_bound_simple_clustering(Node *T1, Node *T2, Forest **out_F1, Forest **out_F2);
int rSPR_branch_and_bound_simple_clustering(Node *T1, Node *T2, bool verbose, map<string, int> *label_map, map<int, string> *reverse_label_map);
int rSPR_branch_and_bound_simple_clustering(Node *T1, Node *T2, bool verbose, map<string, int> *label_map, map<int, string> *reverse_label_map, int min_k, int max_k);
int rSPR_branch_and_bound_simple_clustering(Node *T1, Node *T2);
int rSPR_branch_and_bound_simple_clustering(Node *T1, Node *T2, bool verbose);
int rSPR_branch_and_bound_simple_clustering(Node *T1, Node *T2, bool verbose, int min_k, int max_k);
int rSPR_branch_and_bound_simple_clustering(Node *T1, Node *T2, bool verbose, map<string, int> *label_map, map<int, string> *reverse_label_map, int min_k, int max_k, Forest **out_F1, Forest **out_F2);
void reduction_leaf_mult(Forest *T1, Forest *T2);
void reduction_leaf(Forest *T1, Forest *T2);
void reduction_leaf(Forest *T1, Forest *T2, UndoMachine *um);
bool chain_match(Node *T1_node, Node *T2_node, Node *T2_node_end);
Node *find_subtree_of_approx_distance(Node *n, Forest *F1, Forest *F2, int target_size);
Node *find_best_root(Node *T1, Node *T2);
double find_best_root_acc(Node *T1, Node *T2);
void find_best_root_hlpr(Node *T2, int pre_separator, int group_1_total,
		int group_2_total, Node **best_root, double *best_root_b_acc);
void find_best_root_hlpr(Node *n, int pre_separator, int group_1_total,
		int group_2_total, Node **best_root, double *best_root_b_acc,
		int *p_group_1_descendants, int *p_group_2_descendants, int *num_ties);
int rf_distance(Node *T1, Node *T2);
int count_differing_bipartitions(Node *n);
bool contains_bipartition(Node *n, int pre_start, int pre_end,
		int group_1_total, int group_2_total, int *p_group_1_descendants,
		int *p_group_2_descendants);
void modify_bipartition_support(Node *T1, Node *T2, enum RELAXATION relaxed);
void modify_bipartition_support(Node *n, Forest *F1, Forest *F2,
		Node *T1, Node *T2, vector<int> *F1_descendant_counts, enum RELAXATION);
void modify_bipartition_support(Forest *F1, Forest *F2, Node *n1);
bool is_nonbranching(Forest *T1, Forest *T2, Node *T1_a, Node *T1_c, Node *T2_a, Node *T2_c);
bool outgroup_root(Node *T, set<string, StringCompare> outgroup);
bool outgroup_root(Node *n, vector<int> &num_in, vector<int> &num_out);
bool outgroup_reroot(Node *n, vector<int> &num_in, vector<int> &num_out);
void count_in_out(Node *n, vector<int> &num_in, vector<int> &num_out,
		set<string, StringCompare> &outgroup);
void randomize_tree_with_spr(Node* T1, Node* T2, int count);
/*Joel's part*/
int rSPR_branch_and_bound_simple_clustering(Forest *T1, Forest *T2, bool verbose, map<string, int> *label_map, map<int, string> *reverse_label_map);
int rSPR_branch_and_bound_simple_clustering(Forest *T1, Forest *T2);
int rSPR_branch_and_bound_simple_clustering(Forest *T1, Forest *T2, bool verbose);
int rSPR_total_distance(Forest *T1, vector<Node *> &gene_trees);

bool BB = false;
bool APPROX_CHECK_COMPONENT = false;
bool APPROX_REVERSE_CUT_ONE_B = false;
bool APPROX_REVERSE_CUT_ONE_B_2 = false;
bool APPROX_CUT_ONE_B = false;
bool APPROX_CUT_TWO_B = false;
bool APPROX_CUT_TWO_B_ROOT = false;
bool APPROX_EDGE_PROTECTION = false;
bool CUT_ONE_B = false;
bool REVERSE_CUT_ONE_B = false;
bool REVERSE_CUT_ONE_B_2 = false;
bool REVERSE_CUT_ONE_B_3 = false;
bool CUT_TWO_B = false;
bool CUT_TWO_B_ROOT = false;
bool CUT_ALL_B = false;
bool CUT_AC_SEPARATE_COMPONENTS = false;
bool CUT_ONE_AB = false;
bool CLUSTER_REDUCTION = false;
bool PREFER_RHO = false;
bool MAIN_CALL = true;
bool MEMOIZE = false;
bool MULTIFURCATING = false;
bool MULT_4_BRANCH = false;
bool USE_CASE_7 = true;
bool ALL_MAFS = true;
int NUM_CLUSTERS = 0;
int MAX_CLUSTERS = -1;
bool UNROOTED_MIN_APPROX = false;
bool VERBOSE = false;
bool CLAMP = false;
int MAX_SPR = 1000;
int CLUSTER_MAX_SPR = MAX_SPR;
int MIN_SPR = 0;
bool FIND_RATE = false;
bool EDGE_PROTECTION = false;
bool EDGE_PROTECTION_TWO_B = false;
bool ABORT_AT_FIRST_SOLUTION = false;
bool PREORDER_SIBLING_PAIRS = false;
bool DEEPEST_ORDER = false;
bool DEEPEST_PROTECTED_ORDER = false;
bool NEAR_PREORDER_SIBLING_PAIRS = false;
bool LEAF_REDUCTION = false;
bool LEAF_REDUCTION2 = false;
bool SPLIT_APPROX = false;
bool IN_SPLIT_APPROX = false;
int SPLIT_APPROX_THRESHOLD = 25;
float INITIAL_TREE_FRACTION = 0.4;
bool COUNT_LOSSES = false;
bool CUT_LOST = false;
bool CHECK_MERGE_DEPTH = false;
bool check_all_pairs = true;
bool PREFER_NONBRANCHING = false;
int CLUSTER_TUNE = -1;
int SIMPLE_UNROOTED_LEAF = 0;
bool SHOW_PERCENT_LGT_EVENTS = false;
string ALL_MAFS_CASE = "2358";
bool ALL_MERGED_MAFS = true;

class ProblemSolution {
public:
string T1;
string T2;
int k;

ProblemSolution(Forest *t1, Forest *t2, int new_k) {
	T1 = t1->str();
	T2 = t2->str();
	k = new_k;
}
	};

	map<string, ProblemSolution> memoized_clusters = map<string, ProblemSolution>();

/*******************************************************************************
	RSPR WORSE_3_MULT_APPROX
*******************************************************************************/

/* rSPR_worse_3_mult_approx
 * Calculate an approximate maximum agreement forest and SPR distance for two multifurcating trees
 * RETURN At most 3 times the rSPR distance
 * NOTE: destructive. The computed forests replace T1 and T2.
 * T1 and T2 can be  multifurcating forests.
 */
int rSPR_worse_3_mult_approx(Forest *T1, Forest *T2) {
	return rSPR_worse_3_mult_approx(T1, T2, true);
}

int rSPR_worse_3_mult_approx(Forest *T1, Forest *T2, bool sync) {
	// match up nodes of T1 and T2
	if (sync) {
if (!sync_twins(T1, T2))
	return 0;
	}

	if (LEAF_REDUCTION) {
	  reduction_leaf_mult(T1, T2);
	}
//	cout << "T1: "; T1->print_components();
//	cout << "T2: "; T2->print_components();
	// find sibling pairs of T1
	list<Node *> *sibling_groups = T1->find_sibling_groups();


	/*
	list<Node *>::iterator c;
	list<Node *>::iterator ch;
	
	for (c = sibling_groups->begin(); c != sibling_groups->end(); c++) {
	  cout << "Sibling group parent: " << (*c)->str() << endl;
		list<Node *> children = (*c)->get_children();
		for (ch = children.begin(); ch != children.end(); ch++) {
		  cout << "\tChild: " << (*ch)->str() <<endl;
		}
		list<list<Node *>> *identical_sibling_groups = (*c)->find_identical_sibling_groups();
		if (identical_sibling_groups->size() > 0) {
		  cout << "Printing identical sibling group:" << endl;
		  for (auto g = identical_sibling_groups->begin(); g != identical_sibling_groups->end(); g++) {
			cout << "Sibling identical group:" << endl;
			for (auto ge = g->begin(); ge != g->end(); ge++) {
			  cout << "\tchild:" << (*ge)->get_name() << endl;
			}			  
		  }
		}
		else {		  
		  cout << "No identical sibling groups" <<endl;
		}
	}	
	//return 0;
	*/
	
	
	// find singletons of T2
	list<Node *> singletons = T2->find_singletons();
	//list<pair<Forest,Forest> > AFs = list<pair<Forest,Forest> >();

	Forest *F1;
	Forest *F2;

	T2->max_preorder = T2->components[0]->get_max_preorder_number(0);//preorder_number(0);
	int ans = rSPR_worse_3_mult_approx_hlpr(T1, T2, &singletons, sibling_groups, &F1, &F2, true);

	F1->swap(T1);
	F2->swap(T2);
	sync_twins(T1,T2);


	delete sibling_groups;
	delete F1;
	delete F2;
	return ans;
}


// rSPR_worse_3_mult_approx recursive helper function
int rSPR_worse_3_mult_approx_hlpr(Forest *T1, Forest *T2, list<Node *> *singletons, list<Node *> *sibling_groups, Forest **F1, Forest **F2, bool save_forests) {

  int num_cut = 0;
  Node* previous_group = NULL;
  while(!singletons->empty() || !sibling_groups->empty()) {
	  
    // Case 1 - Remove singletons
    while(!singletons->empty()) {
      
      Node *T2_a = singletons->back();
      singletons->pop_back();
      #ifdef DEBUG_APPROX
      cout << "Handling singleton: " << T2_a->str() << endl;
      #endif
      // find twin in T1
      Node *T1_a = T2_a->get_twin();
      // if this is in the first component of T_2 then
      // it is not really a singleton.
      // TODO: problem when we cluster and have a singleton as the
      //		first comp of T2
      //    NEED TO MODIFY CUTTING?
      // 		HERE AND IN BB?
      Node *T1_a_p = T1_a->parent();
      if (T1_a_p == NULL)
	continue;

      
      if (T2_a == T2->get_component(0)){
	if (!T1->contains_rho()) {
	  T1->add_rho();
	  T2->add_rho();
	  num_cut++;
	}
      }
	
	//continue;

      bool is_sibling_group = T1_a_p->is_sibling_group();
      // cut the edge above T1_a
      T1_a->cut_parent();      
      if (!T1_a->is_leaf()) {
	T1_a_p->decrement_non_leaf_children();//although would this ever be a non leaf?
      }
      T1->add_component(T1_a);

      //only contract if one node
      if (T1_a_p->get_children().size() == 1) {
	if (is_sibling_group) {
	  sibling_groups->remove(T1_a_p);
          #ifdef DEBUG_APPROX	  
	  cout << "Removed ";
	  for (list<Node*>::iterator i = T1_a_p->get_children().begin(); i != T1_a_p->get_children().end(); i++) {
	    cout << (*i)->str();
	  }
	  cout  << " from sibling groups" << endl;
	  #endif
	}
	Node *possible_previous_sibling = T1_a_p->get_children().front();
	bool was_sibling_group = possible_previous_sibling->is_sibling_group();
	Node *node = T1_a_p->contract(true);
	if (node != NULL) {
	  node->recalculate_non_leaf_children(); //can we tell what this would be instead of recalculating?

	  if (node->is_sibling_group()) {
	    if (was_sibling_group) {
	      list<Node*>::iterator i = find(sibling_groups->begin(), sibling_groups->end(), possible_previous_sibling);
	      *i = node;
	    }
	    else{
	      sibling_groups->push_front(node);
	    }
	  }
	}
      }
    }//!singletons->empty()

    
    if(!sibling_groups->empty()) {
      //Get the first group that has identical sibling groups, otherwise default to the group on the back
      list<Node*>::iterator i = sibling_groups->end();
      i--;
      Node *T1_sibling_group = sibling_groups->back();
      list<list<Node*>> identical_sibling_groups = list<list<Node*>>();
      T1_sibling_group->find_identical_sibling_groups(&identical_sibling_groups);
      for (; i != sibling_groups->begin(); i-- ){
	(*i)->find_identical_sibling_groups(&identical_sibling_groups);
	if (identical_sibling_groups.size() > 0) {
	  T1_sibling_group = (*i);
	  break;
	}
      }
      
      #ifdef DEBUG_APPROX
      cout << "F2: ";
      T2->print_components();
      cout << endl;
      cout << "F1: ";
      T1->print_components();
      cout << endl;
      #endif
      
      /* 
	 Case where a subset of the group have the same parent both in T1 and T2
	 Step 5 in paper
      */
      // Case 2 - Contract identical sibling pair
      if (identical_sibling_groups.size() > 0) {	  		
	list<list<Node *>>::iterator i;
	for (i = identical_sibling_groups.begin(); i != identical_sibling_groups.end(); i++) {
	  list<Node *> T2_group = (*i);
	  Node *T2_p = T2_group.front()->parent();
	  #ifdef DEBUG_APPROX
	  cout << "Contracting T1... " << endl;
	  #endif
	  Node *T1_group_new = T1_sibling_group->contract_twin_group(&T2_group);
	  #ifdef DEBUG_APPROX
	  cout << "Contracting T2... " << endl;
	  #endif
	  Node *T2_group_new = T2_p->contract_sibling_group(&T2_group);
	  
	  T1_group_new->set_twin(T2_group_new);
	  T2_group_new->set_twin(T1_group_new);			

	  // check if T2_p is a singleton after the contraction
	  if (T2_p->is_singleton() && T2_p != T2->get_component(0)) {
	    //(T2_p->is_singleton() && T1_sibling_group != T1->get_component(0) && T2_p != T2->get_component(0)) {
	    singletons->push_front(T2_p);
	  }
	  if (T1_sibling_group->parent() != NULL) {
	    //Check if the contraction made a new sibling group
	    T1_sibling_group->parent()->recalculate_non_leaf_children();
	    if (T1_sibling_group->parent()->is_sibling_group()) {
	      sibling_groups->push_front(T1_sibling_group->parent());
	    }
	  }
	  if (!T1_sibling_group->is_sibling_group()) {
	    sibling_groups->remove(T1_sibling_group);
	  }	  
	}
      }

      /*
	4 branching case
	Step 6-8 in paper
	Cut above a1, b1, a2, b2
	Part 1: Get the LCA, this is the numbering to the top part
	Part 2: Get subset of sibling group that is descendant of this LCA
	Part 3: Sort them based on depth
	Part 4: Since this is the approximation, we cut deepest 2
      */

      // Case 3
      else {
	#ifdef DEBUG_APPROX
	cout << "Sibling group to be cutting: " << T1_sibling_group->str_subtree() << endl;
	#endif
	//vector of ints describing how many of the siblings are in its descendants, indexed by preorder number
	vector<int> descendant_count = T1_sibling_group->find_pseudo_lca_descendant_count(T2->max_preorder + 1);
	Node* arbitrary_lca = T1_sibling_group->find_arbitrary_lca(T2->components, descendant_count);
	vector<Node *> deepest_siblings;
	
	//All siblings are in different components, ie no path between them
	//Get depth of siblings from root of each component
	if (arbitrary_lca == NULL) {	
	  vector<vector<Node *>> siblings_by_depth = vector<vector<Node *>>(10);
	  for (int i = 0; i != T2->components.size(); i++) {
	    T2->components[i]->get_deepest_siblings(descendant_count, siblings_by_depth);
	  }
	  deepest_siblings = contract_deepest_siblings(siblings_by_depth);
	}
	//Otherwise they share an LCA
	else {
	  vector<vector<Node *>> siblings_by_depth = arbitrary_lca->get_deepest_siblings(descendant_count);
	  deepest_siblings = contract_deepest_siblings(siblings_by_depth);
	}
      
	// Should assert here
	if (deepest_siblings.size() < 2) { cout << "improper length" << endl; }

	//Get the deepest two of the siblings
	Node *T2_a1 = deepest_siblings[0];
	Node *T2_a2 = deepest_siblings[1];
	#ifdef DEBUG_APPROX
	cout << "a1: " << T2_a1->str() << " a2: " << T2_a2->str() << endl;
	#endif
	
	bool cut_a1   = false;
	bool cut_a1_p = false;
	bool cut_a2   = false;
	bool cut_a2_p = false;
      
	if (T1_sibling_group->get_children().size() == 2) {
	  /*
	    7.1 case
	    Cut a1, pa1, a2 in F2, add 3 to num_cut
	    Consider adding 3 regardless if we actually cut 3, 
	  */	
	  cut_a1   = true;
	  cut_a1_p = true;
	  cut_a2   = true;
	  //num_cut += 3;
	  /*
	    7.2 case
	    Cut a1, a2, pa1, pa2, add 4 to num_cuts
	  */
	  if (previous_group == T1_sibling_group) {
	    cut_a2_p = true;
	    #ifdef DEBUG_APPROX
	    cout <<"Case 7.2" << endl;
	    #endif
	    //num_cut += 1;
	  }
	  else {
	    #if DEBUG_APPROX
	    cout << "Case 7.1" << endl;
	    #endif
	  }
	} // size == 2
      
	else if (T1_sibling_group->get_children().size() > 2) {
	  /*
	    7.3 case
	    If a2's parent's only sibling is part of the sibling group, 
	    and a1's parent is a root or has a sibling that is not part of the sibling group
	    then cut a2 and a2_p otherwise a1 and a1_p
	  */
	  if (previous_group != T1_sibling_group) {
	    Node* T2_a2_p = T2_a2->parent();
	    bool x_2 = false;
	    bool a2_p_one_sibling = (T2_a2_p != NULL) &&
	      (T2_a2_p->parent() != NULL) &&
	      (T2_a2_p->parent()->get_children().size() == 2);	  
	    if (a2_p_one_sibling) {
	      list<Node *> group = T1_sibling_group->get_children();
	      //get the other one
	      Node *a2_p_sibling = T2_a2_p->parent()->get_children().front() == T2_a2_p ?
		T2_a2_p->parent()->get_children().back() :
		T2_a2_p->parent()->get_children().front();	
	      //check if it is part of sibling group
	      bool a2_p_sibling_in_group = descendant_count[a2_p_sibling->get_preorder_number()] == -1;
	      if (a2_p_sibling_in_group) {
		Node* T2_a1_p = T2_a1->parent();
		bool a1_p_is_root = T2_a1_p->parent() == NULL;
		bool a1_p_sibling_not_in_group = false;
		if (!a1_p_is_root) {
		  list<Node*> a1_p_siblings = T2_a1_p->parent()->get_children();
		  for (list<Node*>::iterator i = a1_p_siblings.begin(); i != a1_p_siblings.end(); i++) {
		    if (descendant_count[(*i)->get_preorder_number()] != -1) {
		      a1_p_sibling_not_in_group = true;
		      break;
		    }
		  }
		}
		if (a2_p_one_sibling && a2_p_sibling_in_group && (a1_p_is_root || a1_p_sibling_not_in_group)) {
		  x_2 = true;

		}
	      }
	    }
	    #ifdef DEBUG_APPROX
	    cout << "Case 7.3" << endl;
	    #endif
	    //num_cut += 2;
	    if (x_2){
	    
	      cut_a2   = true;
	      cut_a2_p = true;
	    }
	    else {
	      cut_a1   = true;
	      cut_a1_p = true;
	    }
	  }
	  /*7.4 case
	    cut a1 and a1_p
	  */
	  else if (previous_group == T1_sibling_group) {
	    #ifdef DEBUG_APPROX
	    cout << "Case 7.4" << endl;
	    #endif
	    cut_a1   = true;
	    cut_a1_p = true;
	    //num_cut += 2;
	  }
	}
	/*
	  Cutting section
	*/
	Node* T2_a2_p = NULL;
	if (cut_a1) {
	  Node *T2_a1_p = T2_a1->parent();
	  //singletons? components?
	  if (T2_a1_p != NULL) {
	    //Cut connections
	    T2_a1->cut_parent();
	    num_cut++;
	    //add as components
	    T2->add_component(T2_a1);
	    //Just cut 1 of two children of a1 parent
	    if (T2_a1_p->get_children().size() == 1) {
	      if (T2_a1_p->parent() == NULL) {
		T2_a1_p->contract(true);
		if (T2_a1_p->is_singleton() && T2_a1_p != T2->get_component(0)) {
		  singletons->push_front(T2_a1_p);
		}
	      }
	      else {
		Node* T2_b1 = T2_a1_p->get_children().front();
		T2_a1_p->contract(true);
		T2_a1_p = T2_b1; // <------------------------------------
	      }
	    }
	    //Check for singletons
	    if (T2_a1->is_singleton()) // wont ever be C0?
	      singletons->push_front(T2_a1);
	    

	    if (cut_a1_p) {
	      if (T2_a1_p->parent() != NULL) {
		bool aborted_a2 = false;
		if (T2_a1_p->parent() == T2_a2->parent() &&
		    T2_a2->parent()->get_children().size() == 2)
		  {
		    cut_a2 = false; //cutting a1_p will cause a2 to get contracted up, so there is no more a2 to cut
		    cut_a2_p = true; //Instead we cut a2_p
		    aborted_a2 = true;
		  }
		Node *T2_a1_gp = T2_a1_p->parent();
		//Cut connections
		T2_a1_p->cut_parent();
		num_cut++;
		//add as components
		T2->add_component(T2_a1_p);
		//Just cut 1 of two children of a1 parent
		if (T2_a1_gp->get_children().size() == 1) {
		  if (T2_a1_gp->parent() == NULL) {
		    T2_a1_gp->contract(true);		    
		  }
		  else {
		    Node* T2_b1_p = T2_a1_gp->get_children().front();
		    T2_a1_gp->contract(true);
		    T2_a1_gp = T2_b1_p;
		  }
		  if (aborted_a2) { T2_a2_p = T2_a1_gp; }
		  if (T2_a1_gp != NULL && T2_a1_gp->is_singleton() && (T2_a1_gp != T2->get_component(0) || aborted_a2)) {
		    singletons->push_front(T2_a1_gp);
		  }
		}
		//Check for singletons
		if (T2_a1_p->is_singleton() && T2_a1_p != T2->get_component(0)) {
		  singletons->push_front(T2_a1_p);	      
		}	
	      }
	    }
	  }//T2_a1_p() != NULL
	}
	if (cut_a2){
	  T2_a2_p = T2_a2->parent();
	  //could have cut a2's parent in previous steps, so could be singleton now
	  if (T2_a2->is_leaf() && T2_a2_p == NULL) {
	    singletons->push_front(T2_a2);
	  }
	  else if (T2_a2_p != NULL) {
	    //Cut connections
	    T2_a2->cut_parent();
	    num_cut++;
	    //add as components
	    T2->add_component(T2_a2);
	    //Just cut 1 of two children of a2 parent
	    if (T2_a2_p->get_children().size() == 1) {
	      if (T2_a2_p->parent() == NULL) {
		T2_a2_p->contract(true);
		if (T2_a2_p->is_singleton() && T2_a2_p != T2->get_component(0)) {
		  singletons->push_front(T2_a2_p);
		}
	      }
	      else {
		Node* T2_b2 = T2_a2_p->get_children().front();
		T2_a2_p->contract(true);
		T2_a2_p = T2_b2;
	      }
	    }
	    //Check for singletons
	    if (T2_a2->is_singleton())
	      singletons->push_front(T2_a2);
	    
	  }
	}
	if (cut_a2_p) {
	  if (T2_a2_p != NULL && T2_a2_p->parent() != NULL) {
	    Node *T2_a2_gp = T2_a2_p->parent();
	    //Cut connections
	    T2_a2_p->cut_parent();
	    	    num_cut++;
	    //add as components
	    T2->add_component(T2_a2_p);
	    //Just cut 1 of two children of a2 parent
	    if (T2_a2_gp->get_children().size() == 1) {
	      if (T2_a2_gp->parent() == NULL) {
		T2_a2_gp->contract(true);
	      }
	      else {
		Node* T2_b2_p = T2_a2_gp->get_children().front();
		T2_a2_gp->contract(true);
		T2_a2_gp = T2_b2_p;
	      }

	      if (T2_a2_gp->is_singleton() && T2_a2_gp != T2->get_component(0)) {
		singletons->push_front(T2_a2_gp);
	      }
	    }
	    //Check for singletons
	    if (T2_a2_p->is_singleton()) {
	      singletons->push_front(T2_a2_p);	      
	    }

	  }	 	 
	}
      

      }//else
      //delete identical_sibling_groups;
      previous_group = T1_sibling_group;
    }//!sibling_groups->empty()

  } //while(!sibling_groups->empty() && !singletons->empty()
  // if the first component of the forests differ then we have cut p
  /*
  if (T1->get_component(0)->get_twin() != T2->get_component(0)) {
    if (!T1->contains_rho()) {
      T1->add_rho();
      T2->add_rho();
    }
    else
      // hack to ignore rho when it shouldn't be in a cluster
      num_cut -=3;
  }
  */
  if (save_forests) {
    *F1 = new Forest(T1);
    *F2 = new Forest(T2);
  }
  //if (num_cut) num_cut--;
  return num_cut;
}

#ifdef DEBUG_CASE_COUNTER
struct mult_case_counter {
  int case_71, case_72, case_73, case_74;
  int case_81, case_82, case_83, case_84;
  int case_85, case_86, case_87;
};
static mult_case_counter case_counter = {};
void print_mult_case_count() {
  cout << "Stats:" << endl;
  cout << "\tCase 7.1:" << case_counter.case_71 << endl;
  cout << "\tCase 7.2:" << case_counter.case_72 << endl;
  cout << "\tCase 7.3:" << case_counter.case_73 << endl;
  cout << "\tCase 7.4:" << case_counter.case_74 << endl;
  cout << "\tCase 8.1:" << case_counter.case_81 << endl;
  cout << "\tCase 8.2:" << case_counter.case_82 << endl;
  cout << "\tCase 8.3:" << case_counter.case_83 << endl;
  cout << "\tCase 8.4:" << case_counter.case_84 << endl;
  cout << "\tCase 8.5:" << case_counter.case_85 << endl;
  cout << "\tCase 8.6:" << case_counter.case_86 << endl;
  cout << "\tCase 8.7:" << case_counter.case_87 << endl;
}
  
#endif

int rSPR_branch_and_bound_mult_range(Forest *T1, Forest *T2, int start_k){
  return rSPR_branch_and_bound_mult_range(T1, T2, start_k, MAX_SPR);
}
int rSPR_branch_and_bound_mult_range(Forest *T1, Forest *T2, int start_k, int end_k){
  int exact_spr = -1;
  int k;
  for (k = start_k; k <= end_k; k++) {
    Forest *F1 = new Forest(T1);
    Forest *F2 = new Forest(T2);
    if (!sync_twins(F1,F2)) {
      exact_spr = 0;
      continue;
    }
    if (LEAF_REDUCTION) {
      reduction_leaf_mult(F1, F2);
    }
    #ifdef DEBUG
    cout << "Trying K = " << k << endl << "------------------" << endl;
    #else
    cout << k << " " << endl;
    #endif
    exact_spr = rSPR_branch_and_bound_mult(F1, F2, k);
    #ifdef DEBUG
    cout << "------------------" << endl;
    cout << "Finished K = " << k << " return value : " << exact_spr << endl;   
    #endif
    if (exact_spr >= 0) {
      F1->swap(T1);
      F2->swap(T2);
    }
    delete F1;
    delete F2;

    if (exact_spr >= 0) {    
      break;
    }
  }
  #ifdef DEBUG_CASE_COUNTER
  print_mult_case_count();
  #endif
  if (k > end_k) {
    k = -1;
  }
  return k;
}
int rSPR_branch_and_bound_mult(Forest *T1, Forest *T2, int k){

  if (!sync_twins(T1,T2)) {
    return 0;      
  }
  T2->max_preorder = T2->components[0]->get_max_preorder_number(0);//preorder_number(0);          
  list<Node *> *sibling_groups = T1->find_sibling_groups();
  //sibling_groups->push_front(sibling_groups->back());
  //sibling_groups->pop_back();
  list<Node *> singletons     = T1->find_singletons();
  
  list<pair<Forest,Forest>> AFs = list<pair<Forest,Forest>>();
  //list<Node *> protected_stack = list<Node*>();
  int num_ties = 2;
  int final_k = rSPR_branch_and_bound_mult_hlpr(T1, T2, k, sibling_groups, &singletons, NULL, &AFs, &num_ties);  

  //print AFs
  if (!AFs.empty() && final_k > -1) {
    if (ALL_MAFS) {
      cout << endl << endl << "FOUND ANSWER" << endl;
      // TODO: this is a cheap hack
      for (list<pair<Forest,Forest> >::iterator x = AFs.begin(); x != AFs.end(); x++) {
	cout << "\tT1: ";
	x->first.print_components();
	cout << "\tT2: ";
	x->second.print_components();
      }
    }
      AFs.front().first.swap(T1);
  AFs.front().second.swap(T2);
  sync_twins(T1,T2);

  }

  
  delete sibling_groups;
  if (final_k >= 0)
  return k - final_k;
  else
    return final_k;
    
}


class bb_mult_recurse_data {
public:
  Forest *T1;
  Forest *T2;
  list<Node*> *sibling_groups;
  list<Node*> *singletons;
  map<Node *, Node*> node_map;
};

/* Generates a copy of the data used to recurse on rSPR_branch_and_bound_mult_hlpr so 
   that original forests aren't clobbered,
   updates pointers accordingly */
//TODO: 1 new not 5
bb_mult_recurse_data *generate_mult_recurse_data(Forest *T1, Forest *T2, list<Node*> *sibling_groups, list<Node*> *singletons) {

  bb_mult_recurse_data *data = new bb_mult_recurse_data();

  data->node_map = map<Node*, Node*>();
  data->T1 = new Forest(T1, &data->node_map);
  data->T2 = new Forest(T2, &data->node_map);
  //TODO: smarter way of syncing
  //for (auto n = node_map.begin(); n != node_map.end(); n++) {
    //(*n).first->set_twin(node_map[(*n).first->get_twin()]);
    /*
      if ((*n).second->is_leaf()) {
      (*n).second->set_twin(node_map[(*n).first->get_twin()]);
      node_map[(*n).first->get_twin()]->set_twin((*n).second);
      cout << "Synced twin : " << (*n).second->str() << " -> " << (*n).second->get_twin()->str() << endl;}*/
  //}
  sync_twins(data->T1, data->T2);
  data->sibling_groups = new list<Node*>();
  for (auto n = sibling_groups->begin(); n != sibling_groups->end(); n++) {
    data->sibling_groups->push_back(data->node_map[*n]);
  }
  data->singletons = new list<Node*>();
  for (auto n = singletons->begin(); n != singletons->end(); n++) {
    data->singletons->push_back(data->node_map[*n]);
  }
  return data;
}
//Cuts to_cut, adds to components, conditionally adds to singletons
//Assumes parent is not null and no Null parameters
void mult_cut_and_cleanup(Node* to_cut, Forest *T2, list<Node*> *singletons) {
  //Node* T2_a1_next = next_data->node_map[T2_a1];
  Node* to_cut_p = to_cut->parent();
  //Cut connections  
  to_cut->cut_parent();
  //add as components
  T2->add_component(to_cut);
  //Just cut 1 of two children of a1 parent
  if (to_cut_p->get_children().size() == 1) {
    if (to_cut_p->parent() == NULL) {
      to_cut_p->contract(true);
      if (to_cut_p->is_singleton() && to_cut_p != T2->get_component(0)) {
      singletons->push_front(to_cut_p);
      }
    }
    else {
      Node* to_cut_b = to_cut_p->get_children().front();
      to_cut_p->contract(true);
      if (to_cut_b->is_singleton() && to_cut_b != T2->get_component(0)) {
	singletons->push_front(to_cut_b);
      }
    }
  }
  //Check for singletons
  if (to_cut->is_singleton())
    singletons->push_front(to_cut);  
}
//cuts everything except node, possibly expanding, adds to components, conditionally adds to singletons
//Assumes parent is not null and no Null parameters
//NOTE: potentially unsafe for preorder numbers
void mult_cut_all_except_and_cleanup(Node* T2_a1, Forest *T2, list<Node*> *singletons) {
  Node* T2_b1;
  Node* parent = T2_a1->parent();
  if (parent->get_children().size() == 2) {
    T2_b1 = parent->get_children().front() == T2_a1 ?
      parent->get_children().back() :
      parent->get_children().front();		
  }
  else {
    list<Node*> all_but_a1 = list<Node*>(parent->get_children());
    all_but_a1.remove(T2_a1);	    
    T2_b1 = parent->expand_children_out(all_but_a1);
    T2_b1->set_preorder_number(parent->get_preorder_number());
  }
  //hack for now
  //We immediately contract a1 up so we know the parent's preorder number is available
  //Cut connections

  T2_b1->cut_parent();

  //add as components
  T2->add_component(T2_b1);
	    
  //Just cut 1 of two children of a2 parent (will always be in this case?)
  if (parent->get_children().size() == 1) {
    if (parent->parent() == NULL) {
      parent->contract(true);
      if (parent->is_singleton() && parent != T2->get_component(0)) {
	singletons->push_front(parent);
      }
    }
    else {
      parent = parent->contract(true);
      if (T2_a1->is_singleton() && T2_a1 != T2->get_component(0))
	singletons->push_front(T2_a1);
    }
  }
  if (T2_b1->is_singleton())
    singletons->push_front(T2_b1);
}

//Adds rho and a certain singleton to cut. Basically its for realising when we have cut rho but it is already a singleton so we explicitly continue on 
#define MULT_RHO_CUT_AND_RESOLVE(rho_to_singleton, node_to_protect) {			\
  bb_mult_recurse_data *next_data = generate_mult_recurse_data(T1, T2, sibling_groups, singletons); \
  Node* protect = NULL;							\
  next_data->T1->add_rho();						\
  next_data->T2->add_rho();						\
  next_data->singletons->push_front(next_data->node_map[rho_to_singleton]);	\
  int result_k = rSPR_branch_and_bound_mult_hlpr(next_data->T1, next_data->T2, k - 1, next_data->sibling_groups, next_data->singletons, protect, AFs, num_ties); \
  delete next_data->T1;							\
  delete next_data->T2;						\
  delete next_data->sibling_groups;				\
  delete next_data->singletons;					\
  delete next_data;						\
  if (result_k > best_k) {					\
    best_k = result_k;						\
  }								\
  }								\

//TODO (Ben): try using a routine instead of a macro, I suspect it will slow down because of
//            stack and parameter passing though
#define MULT_BB_CUT_AND_RESOLVE(nodes_to_cut, nodes_to_exclude_cutting, node_to_protect) { \
  bb_mult_recurse_data *next_data = generate_mult_recurse_data(T1, T2, sibling_groups, singletons);\
  int num_cuts = 0;							\
  for (int i = 0; i < nodes_to_cut.size(); i++) {			\
    Node* T2_ax_next = next_data->node_map[nodes_to_cut[i]];		\
    mult_cut_and_cleanup(T2_ax_next, next_data->T2, next_data->singletons);\
    num_cuts++;								\
  }									\
  for (int i = 0; i < nodes_to_exclude_cutting.size(); i++) {		\
    Node* T2_ax_next = next_data->node_map[nodes_to_exclude_cutting[i]]; \
    mult_cut_all_except_and_cleanup(T2_ax_next, next_data->T2, next_data->singletons); \
    num_cuts++;								\
  }									\
  Node* protect = next_data->node_map[node_to_protect];			\
  int result_k = rSPR_branch_and_bound_mult_hlpr(next_data->T1, next_data->T2, k - num_cuts, next_data->sibling_groups, next_data->singletons, protect, AFs, num_ties); \
  delete next_data->T1;							\
  delete next_data->T2;							\
  delete next_data->sibling_groups;					\
  delete next_data->singletons;						\
  delete next_data;							\
  if (result_k > best_k) {						\
    best_k = result_k;							\
    if (!ALL_MAFS && best_k > -1) {					\
    }									\
  }									\
}

//TODO: UndoMachine, then cleanup all constructors, bb_mult_recurse_data relying on copies of trees
int rSPR_branch_and_bound_mult_hlpr(Forest *T1, Forest *T2,
				    int k,
				    list<Node*> *sibling_groups, list<Node*> *singletons,
				    Node *protected_node, list<pair<Forest,Forest>> *AFs,
				    int* num_ties) {
  //run out of cuts if k < 0
  if (k < 0) {
    return k;
  }
  Node* previous_group = sibling_groups->back();
  int best_k = -1;
  while(!singletons->empty() || !sibling_groups->empty()) {
	  
    // Case 1 - Remove singletons
    while(!singletons->empty()) {
      
      Node *T2_a = singletons->back();
      singletons->pop_back();
      #ifdef DEBUG
      cout << "Handling singleton: " << T2_a->str() << endl;
      #endif
      
      Node *T1_a = T2_a->get_twin();
      Node *T1_a_p = T1_a->parent();
      
      if (T1_a_p == NULL)
	continue;      

      // find twin in T1
      //If we have added component 0, then we have made a B cut and added rho
      if (T2_a == T2->get_component(0)){// && T1_a != T1->get_component(0)) {
	if (!T1->contains_rho()) {
	  T1->add_rho();
	  T2->add_rho();
	  k--;
	  //continue;
	}
      }

      bool is_sibling_group = T1_a_p->is_sibling_group();
      // cut the edge above T1_a
      T1_a->cut_parent();
      if (!T1_a->is_leaf()) {
	T1_a_p->decrement_non_leaf_children();
      }
      T1->add_component(T1_a);
      
      //only contract if one node
      if (T1_a_p->get_children().size() == 1) {
	//If we contract this node, and it used to be a sibling group
	//then it is not a sibling group anymore
	if (is_sibling_group) {
	  sibling_groups->remove(T1_a_p);
          #ifdef DEBUG
	  cout << "Removed " << T1_a_p->str() << " from sibling groups" << endl;
	  #endif
	}
	
	Node *possible_previous_sibling = T1_a_p->get_children().front();
	bool was_sibling_group = possible_previous_sibling->is_sibling_group();
	
	Node *T1_a_gp = T1_a_p->parent();
	Node *T1_new_a_p = T1_a_p;
	
	//If the child is a sibling group, there is a possibility contract()
	//will delete it, so we need to update it in the sibling groups
	if (T1_a_p->parent() == NULL) {
	  if (was_sibling_group){
	    list<Node*>::iterator i = find(sibling_groups->begin(), sibling_groups->end(), possible_previous_sibling);
	    *i = T1_new_a_p;	    
	  }
	}
	else {
	  T1_new_a_p = possible_previous_sibling;
	}
       
        T1_a_p->contract(true);
	T1_new_a_p->recalculate_non_leaf_children();
	//After contracting this, the grandparent may be a sibling group now
	if (T1_a_gp != NULL) {
	  T1_a_gp->recalculate_non_leaf_children(); //can we tell what this would be instead of recalculating?
	  if (T1_a_gp->is_sibling_group()) {	 
	    sibling_groups->push_front(T1_a_gp);	 
	  }
	}
      }

    }//!singletons->empty()

    //NOTE: we know there are no singletons left here
    if(!sibling_groups->empty()) {
      //Find identical sibling groups
      //Get the first group that has identical sibling groups, otherwise default to the group on the back
      list<Node*>::reverse_iterator i = sibling_groups->rbegin();
      Node *T1_sibling_group = sibling_groups->back();
      list<list<Node*>> identical_sibling_groups;
      T1_sibling_group->find_identical_sibling_groups(&identical_sibling_groups);
      bool found_identical = false;
      for (; i != sibling_groups->rend(); i++ ){
	(*i)->find_identical_sibling_groups(&identical_sibling_groups);
	if (identical_sibling_groups.size() > 0) {
	  T1_sibling_group = (*i);
	  found_identical = true;
	  break;
	}
      }
      #ifdef MULT_PICK_LARGEST_GROUP
      if (!found_identical) {
	int max_size = 0;
	Node* largest_group = NULL;
	for (auto i = sibling_groups->begin(); i != sibling_groups->end(); i++) {
	  if ((*i)->get_children().size() > max_size) {
	    max_size = (*i)->get_children().size();
	    largest_group = (*i);
	  }
	}
	T1_sibling_group = largest_group;
      }
      #endif
      #ifdef DEBUG
      cout << "K = " << k << endl;
      cout << "F2: ";
      T2->print_components();
      cout << endl;
      cout << "F1: ";
      T1->print_components();
      cout << endl;
      #endif

      /* 
	 Case where a subset of the group have the same parent both in T1 and T2
      */
      // Case 2 - Contract identical sibling pair
      if (identical_sibling_groups.size() > 0) {	  		
	list<list<Node *>>::iterator i;
	for (i = identical_sibling_groups.begin(); i != identical_sibling_groups.end(); i++) {
	  //Contract the groups
	  list<Node *> T2_group = (*i);
	  Node *T2_p = T2_group.front()->parent();
	  #ifdef DEBUG
	  //cout << "Contracting T1... " << endl;
	  #endif
	  Node *T1_group_new = T1_sibling_group->contract_twin_group(&T2_group);
	  #ifdef DEBUG
	  //cout << "Contracting T2... " << endl;
	  #endif
	  Node *T2_group_new = T2_p->contract_sibling_group(&T2_group);

	  //Maintain twins
	  T1_group_new->set_twin(T2_group_new);
	  T2_group_new->set_twin(T1_group_new);			

	  // check if T2_p is a singleton after the contraction
	  if (T2_p->is_singleton() && T2_p != T2->components[0]){
	    singletons->push_front(T2_p);
	  }
	  if (T1_sibling_group->parent() != NULL) {
	    //Check if the contraction made a new sibling group
	    T1_sibling_group->parent()->recalculate_non_leaf_children();
	    if (T1_sibling_group->parent()->is_sibling_group()) {
	      sibling_groups->push_front(T1_sibling_group->parent());
	      #ifdef DEBUG
	      //cout << "Added new sibling group after contraction: " << T1_sibling_group->parent()->str_subtree() << endl;
	      #endif
	    }
	  }
	  if (!T1_sibling_group->is_sibling_group()) {
	    sibling_groups->remove(T1_sibling_group);
	      #ifdef DEBUG
	    //cout << "Removed new sibling group after contraction: " << T1_sibling_group->str() << endl;
	      #endif

	  }	  
	}
      }

      /*
	4 branching case
	Step 6-8 in paper
	Cut above a1, b1, a2, b2
	Part 1: Get the LCA, this is the numbering to the top part
	Part 2: Get subset of sibling group that is descendant of this LCA
	Part 3: Sort them based on depth
	Part 4: Since this is the approximation, we cut deepest 2
      */

      // Case 3
      else {

	//Check branch and bound
	if (BB) {
	  //copies for the approx so we dont clobber this tree
	  map<Node*, Node*> approx_map = map<Node*, Node*>();
	  Forest T1_approx = Forest(T1, &approx_map);
	  Forest T2_approx = Forest(T2, &approx_map);
	  sync_twins(&T1_approx, &T2_approx);
	  list<Node*> sibling_group_approx = list<Node*>();
	  for (auto n = sibling_groups->begin(); n != sibling_groups->end(); n++) {
	    sibling_group_approx.push_back(approx_map[*n]);
	  }
	  list<Node*> singletons_approx = list<Node*>();
	  for (auto n = singletons->begin(); n != singletons->end(); n++) {
	    singletons_approx.push_back(approx_map[*n]);
	  }
	  int approx_spr = rSPR_worse_3_mult_approx_hlpr(&T1_approx, &T2_approx, &singletons_approx, &sibling_group_approx, NULL, NULL, false);
	  //TODO: Sometimes approx returns 1 over the right amount. For example 4 for a 1 cut tree or 16 for a 3 cut tree
	  //If approx_spr > 3k then we will not have enough cuts
	  if (approx_spr > 3*k + 1) {
#ifdef DEBUG
	    cout << "approx failed approx k = " << approx_spr  <<  endl;
#endif
	    return -1;
	  }
	}

	//Finding deepest siblings
	#ifdef DEBUG
	cout << "Sibling group to be cutting: " << T1_sibling_group->str_subtree() << endl;
	#endif
	//vector of ints describing how many of the siblings are in its descendants, indexed by preorder number
	vector<int> descendant_count = T1_sibling_group->find_pseudo_lca_descendant_count(T2->max_preorder + 1);
	Node* arbitrary_lca = T1_sibling_group->find_arbitrary_lca(T2->components, descendant_count);
	vector<Node *> deepest_siblings;
	vector<vector<Node*>> siblings_by_depth;
	//Get depth of siblings from root of each component
	//If the lca is null, all siblings are in different components, ie no path between them
	if (arbitrary_lca == NULL) {	
	  siblings_by_depth = vector<vector<Node *>>(10);
	  for (int i = 0; i != T2->components.size(); i++) {
	    //if preorder is -1 then this is rho
	    if (T2->components[i]->get_preorder_number() == -1) { continue; }
	    T2->components[i]->get_deepest_siblings(descendant_count, siblings_by_depth);
	  }
	}
	//Otherwise they share an LCA
	else {
	  siblings_by_depth = arbitrary_lca->get_deepest_siblings(descendant_count);
	}
	map<Node*, int> s_map = map<Node*, int>();
	//Get deepest siblings returns sparse vector, this compacts it
	deepest_siblings = contract_deepest_siblings(siblings_by_depth, &s_map);      
	// Should assert here
	if (deepest_siblings.size() < 2) { cout << "improper length" << endl; }

	//Get the deepest two of the siblings
	Node *T2_a1 = deepest_siblings[0];
	Node *T2_a2 = deepest_siblings[1];
	#ifdef DEBUG
	cout << "a1: " << T2_a1->str() << " a2: " << T2_a2->str() << endl;
	#endif
	best_k = -1;

	//MULT_4_BRANCH is naive approach of making 4 cuts every time
	if (!MULT_4_BRANCH) {
	//step 7
	/*
	  if a1 == prot, x = 2 else x = 1
	*/
	if (protected_node != NULL) {
	  Node* T2_ax;
	  if (T2_a1 != protected_node) {
	    T2_ax = T2_a1;
	  }
	  else {
	    T2_ax = T2_a2;
	  }
	  
	  bool a0_descendant_of_lca = false;
	  Node* stepper = protected_node;
	  while (stepper->parent() != NULL && stepper->parent() != arbitrary_lca) {
	    stepper = stepper->parent();
	  }
	  if (stepper->parent() == arbitrary_lca) {
	    a0_descendant_of_lca = true;
	  }

	  //TODO (Ben) do one iteration of pl parent for has_aj_child and has_non_aj_child,
	  //           as well as any other flags, before step 7 check
	  bool pl_has_aj_child = false;
	  if (arbitrary_lca != NULL) {
	    Node* pl = arbitrary_lca->parent();	  
	    if (pl != NULL) {
	      for (auto i = pl->get_children().begin(); i != pl->get_children().end(); i++) {
		if (*i != arbitrary_lca && descendant_count[(*i)->get_preorder_number()] == -1) {
		  pl_has_aj_child = true;
		  break;
		}
	      }
	    }
	  }
	  //step 7.1
	  if (arbitrary_lca == NULL) {
	    #ifdef DEBUG
	    cout << "Case 7.1 T2_ax: " << T2_ax->str() << endl;
	    #endif
	    #ifdef DEBUG_CASE_COUNTER
	    case_counter.case_71++;
	    #endif
	    /* 
	       recurse on cut ax prot protected node
	    */
	    //7.1 : cut ax
	    if (T2_ax->parent() != NULL) {
	      vector<Node*> to_cut = {T2_ax};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, protected_node);
	    }
	  }

	  //step 7.2
	  else if (!a0_descendant_of_lca) {
	    /*
	      recurse on cut a1's B's up to LCA prot protected node

	    */
	    #ifdef DEBUG
	    cout << "Case 7.2a cut all B1's Protected Node: " << protected_node->str() <<  endl;
	    #endif
	    #ifdef DEBUG_CASE_COUNTER
	    case_counter.case_72++;
	    #endif

	    {
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {};
	      Node* stepper = T2_a1;
		
	      while(stepper->parent() != arbitrary_lca) { 
		to_cut_except.push_back(stepper);
		stepper = stepper->parent();
	      }
	      //TODO: use list to push_front or figure out how to add from top to bottom 
	      reverse(to_cut_except.begin(), to_cut_except.end());
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, protected_node);
	    }
	    
	    if (s_map[T2_a1] > 1 ||
		s_map[deepest_siblings[deepest_siblings.size()-1]] > 0){
	      /*
		recurse on cut a1 prot protected node
	      */
#ifdef DEBUG
	    cout << "Case 7.2b Cut a1 Protected Node: " << protected_node->str() <<  endl;
#endif
	      //7.2b cut a1
	      if (T2_a1->parent() != NULL) {
		vector<Node*> to_cut = {T2_a1};
		vector<Node*> to_cut_except = {};
		MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, protected_node);
	      }     
	    }
	     /*

	      Not part of the outlined special cases

	     */
	    //if (T2_a1->parent()->get_children().size() == 2)
	    {
	      #ifdef DEBUG
	      cout << "Case 7.2c Cut a2 T2_ax: " << T2_ax->str() << " Protected Node: " << protected_node->str() << endl;
	      #endif
	      vector<Node*> to_cut = {T2_a2};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, protected_node);	      
	    }
	  }	
	    

	  //TODO: Is component 0 considered a root?
	  //If we do get component 0 can we consider this if rho has been added?
	  //step 7.3
	  else if (deepest_siblings.size() == 2 &&
		   T2_ax->parent() == arbitrary_lca &&
		   a0_descendant_of_lca &&
		   (
		    (arbitrary_lca->parent() == NULL && arbitrary_lca != T2->get_component(0)) ||
		    pl_has_aj_child
		    )) {
	    /*
	      recurse on cut Bx prot protected node
	    */
	    //7.3 cut Bx

	    #ifdef DEBUG
	    cout << "Case 7.3 Cut T2_bx T2_ax: " << T2_ax->str() << " Protected Node: " << protected_node->str() << endl;
	    #endif
	    #ifdef DEBUG_CASE_COUNTER
	    case_counter.case_73++;
	    #endif

	    {
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {T2_ax};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, protected_node);
	    }
	    /*

	      This is NOT in the cases. Figure out whats happening with 8.2 with r = 2

	     */
	    /*
	    {
	      cout << "Case 7.3 Cut T2_ax: " << T2_ax->str() << " Protected Node: " << protected_node->str() << endl;
	      vector<Node*> to_cut = {T2_ax};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, protected_node);
	      }
	    */
	  }


	  //step 7.4
	  else if (a0_descendant_of_lca) {
	    /*
	      recurse on cut ax prot protected_node
	                 if m > 2 and r > 2 recurse on 
			    cut all Bx's then B`x
	     */
	    //7.4 cut ax
	    #ifdef DEBUG
	    cout << "Case 7.4 T2_ax: " << T2_ax->str() << " Protected Node: " << protected_node->str() << endl;
	    #endif
	    #ifdef DEBUG_CASE_COUNTER
	    case_counter.case_74++;
	    #endif

	    {
	      vector<Node*> to_cut = {T2_ax};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, protected_node);
	    }
	    if (T1_sibling_group->get_children().size() > 2 && deepest_siblings.size() > 2) {
	      #ifdef DEBUG
	      cout << "Case 7.4b T2_ax: " << T2_ax->str() << " Protected Node: " << protected_node->str() << endl;
	      #endif
	      //7.4 cut all Bx B`x
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {};
	      Node* stepper = T2_ax;
		
	      while(stepper->parent() != arbitrary_lca) { //Do we know ax is descendant of lca?
		to_cut_except.push_back(stepper);
		stepper = stepper->parent();
	      }
	      //TODO: use list to push_front or figure out how to add from top to bottom 
	      reverse(to_cut_except.begin(), to_cut_except.end());
	      to_cut_except.push_back(T2_ax);
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, protected_node);	      
	    }
	    /*

	      Not part of the outlined special cases

	     */
	    {
	      vector<Node*> to_cut = {T2_a2};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, protected_node);	     
	    }
	    {
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {T2_a1};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, protected_node);	     
	    }

	    /*
	    if (T1_sibling_group->get_children().size() == 2 &&
		deepest_siblings.size() == 2 &&
		T2_a1 == T2_ax &&
		T2_a2 != protected_node) {
	      vector<Node*> to_cut = {T2_a2};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, protected_node);
	      }*/
	  }
	}

	//step 8
	else {
	  bool all_but_ar_s1 = true;
	  for (int i = 0; i < deepest_siblings.size() - 1; i++) {
	    if (s_map[deepest_siblings[i]] != 1) {
	      all_but_ar_s1 = false;
	      break;
	    }
	  }
	  bool all_s1 = s_map[deepest_siblings[deepest_siblings.size()-1]] == 1 &&
	    all_but_ar_s1;

	  bool lca_p_contains_sibling = false;	    
	  if (arbitrary_lca != NULL && arbitrary_lca->parent() != NULL) {
	    for (list<Node*>::iterator i = arbitrary_lca->parent()->get_children().begin();
		 i != arbitrary_lca->parent()->get_children().end();
		 i++) {
	      if (descendant_count[(*i)->get_preorder_number()] == -1) {
		lca_p_contains_sibling = true;
		break;
	      }
	    }
	  }

	  //step 8.1
	  if (arbitrary_lca == NULL) {
	    /*
	      recurse on cut a1 no prot,
	      cut a2 no prot
	    */
#ifdef DEBUG
	    cout << "Case 8.1a cut a1" << endl;
#endif
#ifdef DEBUG_CASE_COUNTER
	    case_counter.case_81++;
#endif

	    //8.1 : Cut a1
	    if (T2_a1->parent() != NULL) {
	      vector<Node*> to_cut = {T2_a1};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }
	    else if (T2_a1 == T2->get_component(0)) {
	      if (!T1->contains_rho()) {
		MULT_RHO_CUT_AND_RESOLVE(T2_a1, NULL);
	      }
	    }
#ifdef DEBUG
	    cout << "Case 8.1b cut a2" << endl;
#endif
	    //8.1 : Cut a2
	    if (T2_a2->parent() != NULL) {
	      vector<Node*> to_cut = {T2_a2};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }
	    else if (T2_a2 == T2->get_component(0)) {
	      if (!T1->contains_rho()) {
		MULT_RHO_CUT_AND_RESOLVE(T2_a2, NULL);
	      }
	    }
	  
	  }
	  //step 8.2
	  else if (all_but_ar_s1 &&
		   deepest_siblings[deepest_siblings.size()-1]->parent() == arbitrary_lca) {
	    /*
	      recurse on cut all B's except shallowest
	      for each sibling except for shallowest,
	      cut all other B's protect ai
		           
	    */
	    //8.2 B's
	    {
#ifdef DEBUG
	    cout << "Case 8.2a cut B1 - B(r-1)" << endl;
#endif
#ifdef DEBUG_CASE_COUNTER
	    case_counter.case_82++;
#endif

	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {};
	      for (int i = 0; i < deepest_siblings.size() - 1; i++) {
		to_cut_except.push_back(deepest_siblings[i]);
	      }
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }
	    //all other B's except ai
	    //Doesnt explicitly say this in the paper, but
	    //this would only happen if r > 2


	    if (deepest_siblings.size() > 2){
	      for (int i = 0; i < deepest_siblings.size() - 1; i++) {
#ifdef DEBUG
	      cout << "Case 8.2b for all ai cut all other B" << endl;
#endif

		vector<Node*> to_cut = {};
		vector<Node*> to_cut_except = {};
		for (int j = 0; j < deepest_siblings.size() - 1; j++) {
		  if (i != j) {
		    to_cut_except.push_back(deepest_siblings[j]);
		  }
		}
		MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, deepest_siblings[i]);
	      }
	      }
	    /*
	    else if (T1_sibling_group->parent() != NULL && T1_sibling_group->parent()->get_children().size() == 2) {
	      Node* gp = T1_sibling_group->parent();
	      Node* aunt = gp->get_children().front() == T1_sibling_group ?
		gp->get_children().back() :
		gp->get_children().front();
	      if (aunt->is_leaf()) {
		Node* a1_p = T2_a1->parent();
		vector<Node*> a1_p_leaves = a1_p->find_leaves();
		Node* T2_aunt = aunt->get_twin();
		for (int i = 0; i < a1_p_leaves.size(); i++) {
		  if (a1_p_leaves[i] == T2_aunt) {
#ifdef DEBUG
		    cout << "Case 8.2c cut a2" << endl;
#endif
		    vector<Node*> to_cut = {T2_a2};
		    vector<Node*> to_cut_except = {};
		    MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
		    break;
		  }
		}

	      }
	      }*/


	      else if (deepest_siblings.size() == 2) {
		  vector<Node*> to_cut = {T2_a2};
		  vector<Node*> to_cut_except = {};
		  MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	      }

	  }
	  

	  //step 8.3
	  else if (T1_sibling_group->get_children().size() == 2 &&
		   s_map[T2_a1] + s_map[T2_a2] >= 2) {
	    /*recurse on cut a1 no prot,
	      cut a2 no prot,
	      cut all along a1 to a2 no prot
	    */
	    //8.3 : Cut a1
	    if (T2_a1->parent() != NULL) {
#ifdef DEBUG
	      cout << "Case 8.3a cut a1" << endl;
#endif
#ifdef DEBUG_CASE_COUNTER
	    case_counter.case_83++;
#endif

	      vector<Node*> to_cut = {T2_a1};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }
	    //8.3 : Cut a2
	    if (T2_a2->parent() != NULL) {
#ifdef DEBUG
	    cout << "Case 8.3b cut a2" << endl;
#endif
	      vector<Node*> to_cut = {T2_a2};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }
	    //8.3 : All along path from a1, a2
	    {
#ifdef DEBUG
	    cout << "Case 8.3c cut all B1's and B2's" << endl;
#endif
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {};
	      Node* stepper = T2_a1;
		
	      while(stepper->parent() != arbitrary_lca) { //no need to check for null! we know theres a path
		to_cut_except.push_back(stepper);
		stepper = stepper->parent();
	      }
	      stepper = T2_a2;
	      while(stepper->parent() != arbitrary_lca) { //no need to check for null! we know theres a path
		to_cut_except.push_back(stepper);
		stepper = stepper->parent();
	      }
	      /*
		Cut from top to bottom,
		Cutting from bottom to top causes the contraction to invalidate
		the node, since contraction is implemented by cutting the parent, then giving
		the child to the parent
	      */
	      //TODO: use list to push_front or figure out how to add from top to bottom 
	      reverse(to_cut_except.begin(), to_cut_except.end());
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }	      
	  }

	  //step 8.4
	  else if (T1_sibling_group->get_children().size() > 2 &&
		   deepest_siblings.size() == 2 &&
		   s_map[T2_a1] == 1 &&
		   s_map[T2_a2] == 1) {
	    /* 
	       recurse on cut a1 and a2 no prot
	       cut b1 and b2 no prot
	       cut all except b1 off of lca, b1
	       cut all except b2 off of lca, b2
	    */
	    //8.4 : Cut a1 and a2
	    {
#ifdef DEBUG
	    cout << "Case 8.4a cut a1 & a2" << endl;
#endif
#ifdef DEBUG_CASE_COUNTER
	    case_counter.case_84++;
#endif

	      vector<Node*> to_cut = {T2_a1, T2_a2};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }
	    //8.4 : Cut a1 and a2's B's
	    {
#ifdef DEBUG
	    cout << "Case 8.4b cut B1 & B2" << endl;
#endif
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {T2_a1, T2_a2};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }

	    //8.4 : All except a1->p, then b1
	    {
#ifdef DEBUG
	    cout << "Case 8.4c cut B1 & B`1" << endl;
#endif
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {T2_a1->parent(), T2_a1};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, T2_a2);
	    }
	    //8.4 : All except a2->p, then b2
	    {
#ifdef DEBUG
	    cout << "Case 8.4d cut B2 & B`2" << endl;
#endif
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {T2_a2->parent(), T2_a2};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, T2_a1);
	    }
	    Node* pl = arbitrary_lca->parent();
	    bool pl_has_non_aj_child = false;
	    if (pl != NULL) {
	      for (auto i = pl->get_children().begin(); i != pl->get_children().end(); i++) {
		if (*i != arbitrary_lca && descendant_count[(*i)->get_preorder_number()] != -1) {
		  pl_has_non_aj_child = true;
		  break;
		}
	      }
	    }
	    bool gpl_has_non_aj_child = false;
	    //place this conditional so we do not unnecessarily iterate
	    //If pl_has_non_aj_child then the if will evaluate to true later on anyways
	    if (!pl_has_non_aj_child) {
	      if (pl != NULL) {
		Node* gpl = pl->parent();		    		    
		if (gpl != NULL) {
		  for (auto i = gpl->get_children().begin(); i != gpl->get_children().end(); i++) {
		    if (*i != pl && descendant_count[(*i)->get_preorder_number()] != -1) {
		      gpl_has_non_aj_child = true;
		      break;
		    }
		  }
		}
	      }
	    }
	    //8.4 special case
	    if (arbitrary_lca->parent() == NULL ||
		pl_has_non_aj_child ||
		gpl_has_non_aj_child) {
	      //8.4 Cut a1 prot a2
	      {
#ifdef DEBUG
	      cout << "Case 8.4e cut a1" << endl;
#endif
		vector<Node*> to_cut = {T2_a1};
		vector<Node*> to_cut_except = {};
		MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, T2_a2);	
	      }
	      //8.4 Cut a2 prot a1
	      {
#ifdef DEBUG
	      cout << "Case 8.4f a2" << endl;
#endif
		vector<Node*> to_cut = {T2_a2};
		vector<Node*> to_cut_except = {};
		MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, T2_a1);	
	      }

	    }
	  }

	  //step 8.5
	  else if (T1_sibling_group->get_children().size() > 2 &&
		   deepest_siblings.size() > 2 &&
		   all_s1) {
	    /*
	      recurse on cut all a1 through ar no prot,
	      cut all b1 through br no prot,
	      for each ai, cut all a's except ai prot ai
	      for each ai, cut all B's except for ai's B prot ai
	    */
	    //8.5 all ai
	    {
#ifdef DEBUG
	    cout << "Case 8.5a cut all a1-ar" << endl;
#endif
#ifdef DEBUG_CASE_COUNTER
	    case_counter.case_85++;
#endif	    
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {};
	      for (int i = 0; i < deepest_siblings.size(); i++) {
		to_cut.push_back(deepest_siblings[i]);
	      }
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }
	    //8.5 cut all bi
	    {
#ifdef DEBUG
	      cout << "Case 8.5b cut all B1-Br" << endl;
#endif
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {};
	      for (int i = 0; i < deepest_siblings.size(); i++) {
		to_cut_except.push_back(deepest_siblings[i]);
	      }
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }
	    //8.5 for each ai cut all a's except ai
	    {	      
	      for (int i = 0; i < deepest_siblings.size(); i++) {
#ifdef DEBUG
		cout << "Case 8.5c cut all a1-ar except ai" << endl;
#endif
		vector<Node*> to_cut = {};
		vector<Node*> to_cut_except = {};
		for (int j = 0; j < deepest_siblings.size(); j++) {
		  if (i != j) {
		    to_cut.push_back(deepest_siblings[j]);
		  }
		}
		MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, deepest_siblings[i]);//TODO protect ai
	      }
	    }
	    //8.5 for each ai cut all b's except ai's
	    {	      
	      for (int i = 0; i < deepest_siblings.size(); i++) {
#ifdef DEBUG
		cout << "Case 8.5d cut all B1-Br except Bi" << endl;
#endif
		vector<Node*> to_cut = {};
		vector<Node*> to_cut_except = {};
		for (int j = 0; j < deepest_siblings.size(); j++) {
		  if (i != j) {
		    to_cut_except.push_back(deepest_siblings[j]);
		  }
		}
		MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, deepest_siblings[i]);//TODO protect ai
	      }
	    }
	  }

	  //step 8.6
	  else if (T1_sibling_group->get_children().size() > 2 &&
		   deepest_siblings.size() == 2 &&
		   s_map[T2_a1] >= 2 &&
		   T2_a2->parent() == arbitrary_lca && //equivalent to s_map[T2_a2] == 0
		   (
		    arbitrary_lca->parent() == NULL ||
		    lca_p_contains_sibling
		    )) {
	    /*
	      recurse on cut a1 no prot,
	      cut all B's leading up to LCA from a1, no prot,
	      cut B2 (essentially a1's whole branch) no prot
	    */
	    //if (T2_a1->parent() != NULL) {
	    //8.6 cut a1
	    {
#ifdef DEBUG	     
	      cout << "Case 8.6a cut a1" << endl;
#endif
#ifdef DEBUG_CASE_COUNTER
	    case_counter.case_86++;
#endif
	      vector<Node*> to_cut = {T2_a1};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }
	    //8.6 cut all B1s
	    {
#ifdef DEBUG	     
	      cout << "Case 8.6b cut all B1's" << endl;
#endif
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {};
	      Node* stepper = T2_a1;
		
	      while(stepper->parent() != arbitrary_lca) { 
		to_cut_except.push_back(stepper);
		stepper = stepper->parent();
	      }
	      //TODO: use list to push_front or figure out how to add from top to bottom 
	      reverse(to_cut_except.begin(), to_cut_except.end());

	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }
	    //8.6 cut B2
	    {
#ifdef DEBUG	     
	      cout << "Case 8.6c cut B2" << endl;
#endif
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {T2_a2};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }	      
	  }

	  //step 8.7
	  else if (T1_sibling_group->get_children().size() > 2 &&
		   s_map[T2_a1] >= 2) {
	    /*
	      recurse on cut a1 no prot,
	      cut a2 prot a1,
	      cut all B1's leading up to LCA no prot,
			   
	    */
	    //8.7 cut a1
	    {
#ifdef DEBUG	     
	      cout << "Case 8.7a cut a1" << endl;
#endif
#ifdef DEBUG_CASE_COUNTER
	    case_counter.case_87++;
#endif
	      vector<Node*> to_cut = {T2_a1};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }	      
	    //8.7 cut a2
	    {
#ifdef DEBUG	     
	      cout << "Case 8.7b cut a2" << endl;
#endif
	      vector<Node*> to_cut = {T2_a2};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, T2_a1);
	    }
	    //8.7 cut all B1s
	    {
#ifdef DEBUG	     
	      cout << "Case 8.7c cut all B1's" << endl;
#endif
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {};
	      Node* stepper = T2_a1;
		
	      while(stepper->parent() != arbitrary_lca) { 
		to_cut_except.push_back(stepper);
		stepper = stepper->parent();
	      }
	      //TODO: use list to push_front or figure out how to add from top to bottom 
	      reverse(to_cut_except.begin(), to_cut_except.end());

	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);		
	    }
	    //8.7 cut all B2s
	    {
#ifdef DEBUG
	    cout << "Case 8.7d";
	    if (deepest_siblings.size() == 2) {
	      cout << "2" << endl;
	    }
	    else {
	      cout << "n" << endl;
	    }
#endif

	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {};
	      Node* stepper = T2_a2;
		
	      while(stepper->parent() != arbitrary_lca) { 
		to_cut_except.push_back(stepper);
		stepper = stepper->parent();
	      }
	      //TODO: use list to push_front or figure out how to add from top to bottom 
	      reverse(to_cut_except.begin(), to_cut_except.end());
	      //8.7 cut all B2's, cut B`2, otherwise r > 2
	      if(deepest_siblings.size() == 2) {
		to_cut_except.push_back(T2_a2); // cut B`2 at the end
	      }
	      //This exception should be caught by 8.6 case.
	      /*
		if (to_cut_except.size() == 0) {
		cout << "\n\n\nto_cut_except empty in 8.7n. Parent is lca ERROR \n\n\n"; 
		}*/

	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, T2_a1);	
	    }
	  }

	}

	}
	else { //4_MULT_BRANCH

	  //To debug: set all to true. Make all cuts all the time
	  //Otherwise uses same logic as approximation algorithm
	  bool cut_a1 = true;
	  bool cut_b1 = true;
	  bool cut_a2 = true;
	  bool cut_b2 = true;

	  if (T1_sibling_group->get_children().size() == 2) {
	    /*
	      7.1 case
	      Cut a1, pa1, a2 in F2, add 3 to num_cut
	    */	
	    cut_a1   = true;
	    cut_b1 = true;
	    cut_a2   = true;
	    /*
	      7.2 case
	      Cut a1, a2, pa1, pa2, add 4 to num_cuts
	    */
	    if (previous_group == T1_sibling_group) {
	      cut_b2 = true;
	    }
	  } // size == 2
      
	  else if (T1_sibling_group->get_children().size() > 2) {
	    /*
	      7.3 case
	      If a2's parent's only sibling is part of the sibling group, 
	      and a1's parent is a root or has a sibling that is not part of the sibling group
	      then cut a2 and a2_p otherwise a1 and a1_p
	    */
	    if (previous_group != T1_sibling_group) {
	      Node* T2_a2_p = T2_a2->parent();
	      bool x_2 = false;
	      bool a2_p_one_sibling = (T2_a2_p != NULL) &&
		(T2_a2_p->parent() != NULL) &&
		(T2_a2_p->parent()->get_children().size() == 2);	  
	      if (a2_p_one_sibling) {
		list<Node *> group = T1_sibling_group->get_children();
		//get the other one
		Node *a2_p_sibling = T2_a2_p->parent()->get_children().front() == T2_a2_p ?
		  T2_a2_p->parent()->get_children().back() :
		  T2_a2_p->parent()->get_children().front();	
		//check if it is part of sibling group
		bool a2_p_sibling_in_group = descendant_count[a2_p_sibling->get_preorder_number()] == -1;
		if (a2_p_sibling_in_group) {
		  Node* T2_a1_p = T2_a1->parent();
		  bool a1_p_is_root = T2_a1_p->parent() == NULL;
		  bool a1_p_sibling_not_in_group = false;
		  if (!a1_p_is_root) {
		    list<Node*> a1_p_siblings = T2_a1_p->parent()->get_children();
		    for (list<Node*>::iterator i = a1_p_siblings.begin(); i != a1_p_siblings.end(); i++) {
		      if (descendant_count[(*i)->get_preorder_number()] != -1) {
			a1_p_sibling_not_in_group = true;
			break;
		      }
		    }
		  }
		  if (a2_p_one_sibling && a2_p_sibling_in_group && (a1_p_is_root || a1_p_sibling_not_in_group)) {
		    x_2 = true;
		  }
		}
	      }
	      if (x_2){
		cut_a2 = true;
		cut_b2 = true;
		//cut_a1 = true;
		//cut_b1 = true;
	      }
	      else {
		cut_a1 = true;
		cut_b1 = true;
		cut_a2 = true;
		//cut_b2 = true;
	      }
	    }
	    /*7.4 case
	      cut a1 and a1_p
	    */
	    else if (previous_group == T1_sibling_group) {
	      cut_a1   = true;
	      cut_b1 = true;
	    }
	  }
	  /*
	    Cutting section
	  */
	  Node *T2_a1_p = T2_a1->parent();

	  if (cut_a1) {	  
	    if (T2_a1_p != NULL) {
#ifdef DEBUG
	      cout << "Case Cut a1" << endl;
#endif
	      vector<Node*> to_cut = {T2_a1};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }
	    //mimics case 8.1, should add rho
	    else if (T2_a1 == T2->get_component(0) && arbitrary_lca == NULL) {
	      if (!T1->contains_rho()) {
		MULT_RHO_CUT_AND_RESOLVE(T2_a1, NULL);
	      }
	    }
	  }
	  if (cut_b1) {
	    if (T2_a1_p != NULL) {
#ifdef DEBUG
	      cout << "Case Cut b1" << endl;
#endif
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {T2_a1};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }
	  }
	  Node *T2_a2_p = T2_a2->parent();
	  if (cut_a2){
#ifdef DEBUG
	      cout << "Case Cut a2" << endl;
#endif
	    if (T2_a2_p != NULL) {
	      vector<Node*> to_cut = {T2_a2};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }
	    else if (T2_a2 == T2->get_component(0) && arbitrary_lca == NULL) {
	      
	      if (!T1->contains_rho()) {
		MULT_RHO_CUT_AND_RESOLVE(T2_a2, NULL);
	      }
	    }

	  }
	  if (cut_b2) {
	    if (T2_a2_p != NULL) {
#ifdef DEBUG
	      cout << "Case Cut b2" << endl;
#endif
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {T2_a2};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }
	  }
	}
	sibling_groups->pop_back();
	return best_k;
      }//else (cutting)
      previous_group = T1_sibling_group;
    }//!sibling_groups->empty()

  } //while(!sibling_groups->empty() && !singletons->empty()

  //Made it to end, so add to results
  if (k >= 0) {
    //if (PREFER_RHO && !AFs->empty() && !AFs->front().first.contains_rho() && T1->contains_rho()) {
    if (true && !AFs->empty() && !AFs->front().first.contains_rho() && T1->contains_rho()) {
      if (!ALL_MAFS)
	AFs->clear();
      AFs->push_front(make_pair(Forest(T1),Forest(T2)));
      *num_ties = 2;
    }
    else if (ALL_MAFS || AFs->empty()) {
      AFs->push_back(make_pair(Forest(T1),Forest(T2)));
    }
    //else if (!PREFER_RHO || AFs->front().first.contains_rho() == T1->contains_rho()) {
    else if (false || AFs->front().first.contains_rho() == T1->contains_rho()) {
      if (rand() < RAND_MAX/ *num_ties) {
	AFs->clear();
	AFs->push_back(make_pair(Forest(T1),Forest(T2)));
      }
      (*num_ties)++;
    }
  }
  return k;
}


/* rSPR_3_approx
 * Calculate an approximate maximum agreement forest and SPR distance
 * RETURN At most 3 times the rSPR distance
 * NOTE: destructive. The computed forests replace T1 and T2.
 */
int rSPR_3_approx(Forest *T1, Forest *T2) {
	// find sibling pairs of T1
	// match up nodes of T1 and T2
	if (!sync_twins(T1, T2))
return 0;
	// find singletons of T2
	list<Node *> *sibling_pairs = T1->find_sibling_pairs();
	list<Node *> singletons = T2->find_singletons();
	int ans = rSPR_3_approx_hlpr(T1, T2, &singletons, sibling_pairs);
	delete sibling_pairs;
	return ans;
}

// rSPR_3_approx recursive helper function
int rSPR_3_approx_hlpr(Forest *T1, Forest *T2, list<Node *> *singletons,
list<Node *> *sibling_pairs) {
	int num_cut = 0;
	while(!singletons->empty() || !sibling_pairs->empty()) {
// Case 1 - Remove singletons
while(!singletons->empty()) {


	Node *T2_a = singletons->back();
	singletons->pop_back();
  // find twin in T1
	Node *T1_a = T2_a->get_twin();
	// if this is in the first component of T_2 then
	// it is not really a singleton.
	if (T2_a == T2->get_component(0))
		continue;

	Node *T1_a_parent = T1_a->parent();
	if (T1_a_parent == NULL)
		continue;
	bool potential_new_sibling_pair = T1_a_parent->is_sibling_pair();
	// cut the edge above T1_a
	T1_a->cut_parent();
	T1->add_component(T1_a);
	if (T1_a->get_sibling_pair_status() > 0)
		T1_a->clear_sibling_pair(sibling_pairs);
	//delete(T1_a);

	Node *node = T1_a_parent->contract();
	if (node != NULL && potential_new_sibling_pair && node->is_sibling_pair()){
		node->rchild()->add_to_front_sibling_pairs(sibling_pairs, 2);
		node->lchild()->add_to_front_sibling_pairs(sibling_pairs, 1);
	}

}
if(!sibling_pairs->empty()) {
	Node *T1_a = sibling_pairs->back();
	sibling_pairs->pop_back();
	Node *T1_c = sibling_pairs->back();
	sibling_pairs->pop_back();
	T1_a->clear_sibling_pair_status();
	T1_c->clear_sibling_pair_status();
	if (T1_a->parent() == NULL || T1_a->parent() != T1_c->parent()) {
		continue;
	}
	Node *T1_ac = T1_a->parent();
	// lookup in T2 and determine the case
	Node *T2_a = T1_a->get_twin();
	Node *T2_c = T1_c->get_twin();

	// Case 2 - Contract identical sibling pair
	if (T2_a->parent() != NULL && T2_a->parent() == T2_c->parent()) {
		Node *T2_ac = T2_a->parent();
		T1_ac->contract_sibling_pair();
		T2_ac->contract_sibling_pair();
		T1_ac->set_twin(T2_ac);
		T2_ac->set_twin(T1_ac);
		T1->add_deleted_node(T1_a);
		T1->add_deleted_node(T1_c);
		T2->add_deleted_node(T2_a);
		T2->add_deleted_node(T2_c);

		// check if T2_ac is a singleton
		if (T2_ac->is_singleton() && !T1_ac->is_singleton() && T2_ac != T2->get_component(0))
			singletons->push_back(T2_ac);
		// check if T1_ac is part of a sibling pair
		if (T1_ac->parent() != NULL && T1_ac->parent()->is_sibling_pair()) {
			T1_ac->parent()->lchild()->add_to_sibling_pairs(sibling_pairs, 1);
			T1_ac->parent()->rchild()->add_to_sibling_pairs(sibling_pairs, 2);
		}
	}
	// Case 3
	else {
		
		//  ensure T2_a is below T2_c
		if (T2_a->get_depth() < T2_c->get_depth()) {
			swap(&T1_a, &T1_c);
			swap(&T2_a, &T2_c);
		}
		else if (T2_a->get_depth() == T2_c->get_depth()) {
			if (T2_a->parent() && T2_c->parent() &&
					(T2_a->parent()->get_depth() <
					T2_c->parent()->get_depth())) {
			swap(&T1_a, &T1_c);
			swap(&T2_a, &T2_c);
			}
		}

		// get T2_b
		Node *T2_ab = T2_a->parent();
		Node *T2_b = T2_ab->rchild();
		if (T2_b == T2_a)
			T2_b = T2_ab->lchild();
		// cut T1_a, T1_c, T2_a, T2_b, T2_c

		bool cut_b_only = false;
		if (T2_a->parent() != NULL && T2_a->parent()->parent() != NULL && T2_a->parent()->parent() == T2_c->parent()) {
			cut_b_only = true;
			T1_a->add_to_sibling_pairs(sibling_pairs,1);
			T1_c->add_to_sibling_pairs(sibling_pairs,2);
		}

		if (!cut_b_only) {
			T1_a->cut_parent();
			T1_c->cut_parent();
			// contract parents
			Node *node = T1_ac->contract();
			// check for T1_ac sibling pair
			if (node != NULL && node && node->is_sibling_pair()){
				node->lchild()->add_to_sibling_pairs(sibling_pairs,1);
				node->rchild()->add_to_sibling_pairs(sibling_pairs,2);
			}
		}

		bool same_component = true;
		if (APPROX_CHECK_COMPONENT)
			same_component = (T2_a->find_root() == T2_c->find_root());

		if (!cut_b_only) {
			T2_a->cut_parent();
			num_cut++;
		}
		bool cut_b = false;
		if (same_component && T2_ab->parent() != NULL) {
			T2_b->cut_parent();
			num_cut++;
			cut_b = true;
		}
		// T2_b will move up after contraction
		else {
			T2_b = T2_b->parent();
		}
		// check for T2 parents as singletons
		Node *node = T2_ab->contract();
		if (node != NULL && node->is_singleton()
				&& node != T2->get_component(0))
			singletons->push_back(node);

		// if T2_c is gone then its replacement is in singleton list
		// contract might delete old T2_c, see where it is
		bool add_T2_c = true;
		T2_c = T1_c->get_twin();
		// ignore T2_c if it is a singleton
		if (T2_c != node && T2_c->parent() != NULL && !cut_b_only) {

			Node *T2_c_parent = T2_c->parent();
			T2_c->cut_parent();
			num_cut++;
			node = T2_c_parent->contract();
			if (node != NULL && node->is_singleton()
					&& node != T2->get_component(0))
				singletons->push_back(node);
		}
		else {
			add_T2_c = false;
		}

		
		if (!cut_b_only)
			T1->add_component(T1_a);
		if (!cut_b_only)
			T1->add_component(T1_c);
		// put T2 cut parts into T2
		if (!cut_b_only) {
			T2->add_component(T2_a);
		}
		// may have already been added
		if (cut_b) {
			T2->add_component(T2_b);
		}
		// problem if c is deleted
		if (add_T2_c) {
			T2->add_component(T2_c);
		}

		// may have already been added
		if (T2_b->is_leaf())
			singletons->push_back(T2_b);

	}
}
	}
// if the first component of the forests differ then we have to cut p
if (T1->get_component(0)->get_twin() != T2->get_component(0)) {
	num_cut++;
	T1->add_rho();
	T2->add_rho();
}
return num_cut;
}
/*******************************************************************************
	RSPR WORSE_3_APPROX
*******************************************************************************/

/* rSPR_worse_3_approx
 * Calculate an approximate maximum agreement forest and SPR distance
 * RETURN At most 3 times the rSPR distance
 * NOTE: destructive. The computed forests replace T1 and T2.
 * T1 must be a binary tree. T2 can be a multifurcating forest.
 */
int rSPR_worse_3_approx(Forest *T1, Forest *T2) {
	return rSPR_worse_3_approx(T1, T2, true);
}

int rSPR_worse_3_approx(Forest *T1, Forest *T2, bool sync) {
	// match up nodes of T1 and T2
	if (sync) {
if (!sync_twins(T1, T2))
	return 0;
	}
//	cout << "T1: "; T1->print_components();
//	cout << "T2: "; T2->print_components();
	// find sibling pairs of T1
	list<Node *> *sibling_pairs = T1->find_sibling_pairs();
	// find singletons of T2
	list<Node *> singletons = T2->find_singletons();
	list<pair<Forest,Forest> > AFs = list<pair<Forest,Forest> >();

	Forest *F1;
	Forest *F2;

	int ans = rSPR_worse_3_approx_hlpr(T1, T2, &singletons, sibling_pairs, &F1, &F2, true);

	F1->swap(T1);
	F2->swap(T2);
	sync_twins(T1,T2);


	delete sibling_pairs;
	delete F1;
	delete F2;
	return ans;
}

int rSPR_worse_3_approx_distance_only(Forest *T1, Forest *T2) {
if (!sync_twins(T1, T2))
	return 0;
	list<Node *> *sibling_pairs = T1->find_sibling_pairs();
	list<Node *> singletons = T2->find_singletons();
	list<pair<Forest,Forest> > AFs = list<pair<Forest,Forest> >();

	int ans = rSPR_worse_3_approx_hlpr(T1, T2, &singletons, sibling_pairs, NULL, NULL, false);

	delete sibling_pairs;
	return ans;
}

int rSPR_worse_3_approx(Node *subtree, Forest *T1, Forest *T2) {
	return rSPR_worse_3_approx(subtree, T1, T2, true);
}

int rSPR_worse_3_approx(Node *subtree, Forest *T1, Forest *T2, bool sync) {
	// match up nodes of T1 and T2
	if (sync) {
if (!sync_twins(T1, T2))
	return 0;
	}
//	cout << "T1: "; T1->print_components();
//	cout << "T2: "; T2->print_components();
	// find sibling pairs of T1
	list<Node *> *sibling_pairs = subtree->find_sibling_pairs();
	// find singletons of T2
	list<Node *> singletons = T2->find_singletons();
	list<pair<Forest,Forest> > AFs = list<pair<Forest,Forest> >();

	Forest *F1;
	Forest *F2;

	int ans = rSPR_worse_3_approx_hlpr(T1, T2, &singletons, sibling_pairs, &F1, &F2, true);

	F1->swap(T1);
	F2->swap(T2);
	sync_twins(T1,T2);


	delete sibling_pairs;
	delete F1;
	delete F2;
	return ans;
}

// rSPR_worse_3_approx recursive helper function
int rSPR_worse_3_approx_hlpr(Forest *T1, Forest *T2, list<Node *> *singletons, list<Node *> *sibling_pairs, Forest **F1, Forest **F2, bool save_forests) {
	#ifdef DEBUG_APPROX
cout << "rSPR_worse_3_approx_hlpr" << endl;
			cout << "\tT1: ";
			T1->print_components_with_twins();
			cout << "\tT2: ";
			T2->print_components_with_twins();
			cout << "sibling pairs:";
			for (list<Node *>::iterator i = sibling_pairs->begin(); i != sibling_pairs->end(); i++) {
				cout << "  ";
				(*i)->print_subtree_hlpr();
			}
			cout << endl;
	#endif
	int num_cut = 0;
	UndoMachine um = UndoMachine();
	while(!singletons->empty() || !sibling_pairs->empty()) {
// Case 1 - Remove singletons
while(!singletons->empty()) {
	#ifdef DEBUG_APPROX
		cout << "Case 1" << endl;
	#endif

	Node *T2_a = singletons->back();
	singletons->pop_back();
	// find twin in T1
	Node *T1_a = T2_a->get_twin();
	// if this is in the first component of T_2 then
	// it is not really a singleton.
	// TODO: problem when we cluster and have a singleton as the
	//		first comp of T2
	//    NEED TO MODIFY CUTTING?
	// 		HERE AND IN BB?
	if (T2_a == T2->get_component(0))
		continue;

	Node *T1_a_parent = T1_a->parent();
	if (T1_a_parent == NULL)
		continue;
	bool potential_new_sibling_pair = T1_a_parent->is_sibling_pair();
	// cut the edge above T1_a
	um.add_event(new CutParent(T1_a));
	T1_a->cut_parent();
	um.add_event(new AddComponent(T1));
	T1->add_component(T1_a);
	//if (T1_a->get_sibling_pair_status() > 0)
	//	T1_a->clear_sibling_pair(sibling_pairs);
	//delete(T1_a);

	ContractEvent(&um, T1_a_parent);
	Node *node = T1_a_parent->contract();
	if (node != NULL && potential_new_sibling_pair &&
			node->is_sibling_pair()){
		um.add_event(new AddToFrontSiblingPairs(sibling_pairs));
		sibling_pairs->push_front(node->rchild());
		sibling_pairs->push_front(node->lchild());
	}

	#ifdef DEBUG_APPROX
			cout << "\tT1: ";
			T1->print_components();
			cout << "\tT2: ";
			T2->print_components();
	#endif
}
if(!sibling_pairs->empty()) {
	/*
	if (PREORDER_SIBLING_PAIRS) {
		T1->get_component(0)->preorder_number();
		list<Node *>::iterator c;
		list<Node *>::iterator best_sib = sibling_pairs->end();
		int best_prenum = INT_MAX;
		list<Node *>::iterator T1_a_i;
		list<Node *>::iterator T1_c_i;
		for(c = sibling_pairs->begin(); c != sibling_pairs->end(); ) {
				T1_c_i = c;
				T1_c = *c;
				c++;
				T1_a_i = c;
				T1_a = *c;
				c++;
				cout << T1_a->str_subtree() << endl;
				cout << T1_c->str_subtree() << endl;
//					if (T1_a->parent() == NULL || T1_a->parent() != T1_c->parent()) {
//						cout << "invalid" << endl;
//						sibling_pairs->erase(T1_c_i);
//						sibling_pairs->erase(T1_a_i);
//						um.add_event(new PopSiblingPair(T1_a, T1_c, sibling_pairs));
//						continue;
//					}
//					else {
				//int prenum = T1_a->parent()->get_preorder_number();
				int prenum = T1_a->get_preorder_number();
				cout << "prenum=" << prenum << endl;
				cout << "old_prenum=" << best_prenum << endl;
				if (prenum < best_prenum) {
					best_sib = T1_c_i; 
					best_prenum = prenum;
				}
				cout << "new_prenum=" << best_prenum << endl;
//					}
		}
		cout << endl;
		if (best_prenum == INT_MAX)
			continue;
		else {
			T1_c_i = best_sib;
			T1_c = *T1_c_i;
			best_sib++;
			T1_a_i = best_sib;
			T1_a = *T1_a_i;
			sibling_pairs->erase(T1_a_i);
			sibling_pairs->erase(T1_c_i);
		}
	}
	else {
	*/
	Node *T1_a = sibling_pairs->back();
	sibling_pairs->pop_back();
	Node *T1_c = sibling_pairs->back();
	sibling_pairs->pop_back();
	um.add_event(new PopSiblingPair(T1_a, T1_c, sibling_pairs));

	//if (T1_a->get_sibling_pair_status() == 0 ||
	//		T1_c->get_sibling_pair_status() == 0) {
	//	continue;
//			}

	//T1_a->clear_sibling_pair_status();
	//T1_c->clear_sibling_pair_status();
	if (T1_a->parent() == NULL || T1_c->parent() == NULL || T1_a->parent() != T1_c->parent()) {
		continue;
	}
	if (!T1_a->can_be_sibling() || !T1_c->can_be_sibling()
	|| num_cut >= INT_MAX - 3) {
		continue;
	}
	Node *T1_ac = T1_a->parent();
	// lookup in T2 and determine the case
	Node *T2_a = T1_a->get_twin();
	Node *T2_c = T1_c->get_twin();

	#ifdef DEBUG_APPROX
		cout << "Fetching sibling pair" << endl;
		T1_ac->print_subtree();
		cout << "T2_a" << ": ";
		cout << " d=" << T2_a->get_depth() << " ";
		T2_a->print_subtree();
		cout << "T1_c" << ": ";
		T1_c->print_subtree();
		cout << "T2_c" << ": ";
		cout << " d=" << T2_c->get_depth() << " ";
		T2_c->print_subtree();
	#endif

	// Case 2 - Contract identical sibling pair
	if (T2_a->parent() != NULL && T2_a->parent() == T2_c->parent()) {
		#ifdef DEBUG_APPROX
			cout << "Case 2" << endl;
			T1->print_components();
			T2->print_components();
		#endif
		Node *T2_ac = T2_a->parent();
		um.add_event(new ContractSiblingPair(T1_ac));
		T1_ac->contract_sibling_pair_undoable();
		um.add_event(new ContractSiblingPair(T2_ac, T2_a, T2_c, &um));
		Node *T2_ac_new = T2_ac->contract_sibling_pair_undoable(T2_a, T2_c);
		if (T2_ac_new != NULL && T2_ac_new != T2_ac) {
			T2_ac = T2_ac_new;
			um.add_event(new CreateNode(T2_ac));
			um.add_event(new ContractSiblingPair(T2_ac));
			T2_ac->contract_sibling_pair_undoable();
		}
		um.add_event(new SetTwin(T1_ac));
		um.add_event(new SetTwin(T2_ac));
		T1_ac->set_twin(T2_ac);
		T2_ac->set_twin(T1_ac);
		//T2_ac->fix_contracted_order();
		//T1->add_deleted_node(T1_a);
		//T1->add_deleted_node(T1_c);
		//T2->add_deleted_node(T2_a);
		//T2->add_deleted_node(T2_c);

		// check if T2_ac is a singleton
		//if (T2_ac->is_singleton() && !T1_ac->is_singleton() && T2_ac != T2->get_component(0))
		if (T2_ac->is_singleton() && T1_ac != T1->get_component(0) && T2_ac != T2->get_component(0))
			singletons->push_back(T2_ac);
		// check if T1_ac is part of a sibling pair
		if (T1_ac->parent() != NULL && T1_ac->parent()->is_sibling_pair()) {
			um.add_event(new AddToSiblingPairs(sibling_pairs));
			sibling_pairs->push_back(T1_ac->parent()->lchild());
			sibling_pairs->push_back(T1_ac->parent()->rchild());
		}
	}
	// Case 3
	else {
		#ifdef DEBUG_APPROX
			cout << "Case 3" << endl;
		#endif
		
		//  ensure T2_a is below T2_c
		if ((T2_a->get_depth() < T2_c->get_depth()
				&& T2_c->parent() != NULL)
				|| T2_a->parent() == NULL) {
			#ifdef DEBUG_APPROX
				cout << "swapping" << endl;
			#endif
		
			swap(&T1_a, &T1_c);
			swap(&T2_a, &T2_c);

		}
		else if (T2_a->get_depth() == T2_c->get_depth()) {
			if (T2_a->parent() && T2_c->parent() &&
					(T2_a->parent()->get_depth() <
					T2_c->parent()->get_depth())) {
			swap(&T1_a, &T1_c);
			swap(&T2_a, &T2_c);
			}
		}

		// get T2_b
		bool multi_node = false;
		Node *T2_ab = T2_a->parent();
		Node *T2_b = T2_ab;
		if (T2_ab->get_children().size() > 2) {
			multi_node = true;
		}
		else {
			T2_b = T2_ab->rchild();
			if (T2_b == T2_a)
				T2_b = T2_ab->lchild();
		}

		#ifdef DEBUG_APPROX
		cout << "T2_b" << ": ";
		cout.flush();
		T2_b->print_subtree();
	#endif
		// cut T1_a, T1_c, T2_a, T2_b, T2_c

		bool cut_a_only = false;
		bool cut_b_only = false;
		bool cut_c_only = false;
		bool cut_b_only_if_not_a_or_c = false;
		if (APPROX_CUT_ONE_B && T2_a->parent() != NULL && T2_a->parent()->parent() != NULL && T2_a->parent()->parent() == T2_c->parent() && !multi_node
						&& (!APPROX_EDGE_PROTECTION || !T2_b->is_protected())) {
			cut_b_only = true;
			um.add_event(new AddToSiblingPairs(sibling_pairs));
			sibling_pairs->push_back(T1_c);
			sibling_pairs->push_back(T1_a);
		}
	if (APPROX_CUT_TWO_B && !cut_b_only && T1_ac->parent() != NULL
						&& (!APPROX_EDGE_PROTECTION || !T2_b->is_protected())) {
		Node *T1_s = T1_ac->get_sibling();
		if (T1_s->is_leaf()) {
			Node *T2_l = T2_a->parent()->parent();
			if (T2_l != NULL) {
				if (T2_c->parent() != NULL && T2_c->parent()->parent() == T2_l
						&& T2_a->parent()->get_children().size() > 2
						&& T2_c->parent()->get_children().size() > 2) {
					if (T2_l->get_sibling() == T1_s->get_twin()) {
						cut_b_only=true;
					}
					else if (T2_l->parent() == NULL &&
							(T2->contains_rho() ||
							 T2->get_component(0) != T2_l)) {
						cut_b_only_if_not_a_or_c=true;
					}
				}
				else if ((T2_l = T2_l->parent()) != NULL
						&& T2_c->parent() == T2_l
						&& T2_a->parent()->get_children().size() > 2
						&& T2_a->parent()->parent()->get_children().size() > 2) {
					if (T2_l->get_sibling() == T1_s->get_twin()) {
						cut_b_only=true;
					}
					else if (T2_l->parent() == NULL &&
							(T2->contains_rho() ||
							 T2->get_component(0) != T2_l)) {
						cut_b_only_if_not_a_or_c=true;
					}
				}
			}
		}
	}
	if (APPROX_REVERSE_CUT_ONE_B && !cut_b_only && T1_ac->parent() != NULL) {
		Node *T1_s = T1_ac->get_sibling();
		if (T1_s->is_leaf()) {
			if (T1_s->get_twin()->parent() == T2_a->parent()//) {
						&& (!APPROX_EDGE_PROTECTION || !T2_c->is_protected())) {
				cut_c_only=true;
			}
			else if (T1_s->get_twin()->parent() == T2_c->parent()//) {
						&& (!APPROX_EDGE_PROTECTION || !T2_a->is_protected())
							&& T2_c->parent()->get_children().size() <= 2) {
				cut_a_only=true;
			}
		}
		else if (APPROX_REVERSE_CUT_ONE_B_2) {
			if (T2_c->parent() != NULL
				&& chain_match(T1_s, T2_c->get_sibling(), T2_a) //)
						&& (!APPROX_EDGE_PROTECTION || !T2_a->is_protected()))
			cut_a_only = true;
		}
	}
	if (APPROX_CUT_TWO_B_ROOT && cut_a_only == false && cut_c_only == false
			&& cut_b_only_if_not_a_or_c == true) {
		cut_b_only = true;
	}
	/*
	if (CUT_LOST) {
		if (T1_a->num_lost_children() > 0
				|| T2_a->num_lost_children() > 0) {
			cut_a_only = true;
			cut_b_only = false;
			cut_c_only = false;
			num_cut-=3;
		}
		else if (T1_c->num_lost_children() > 0
				|| T2_c->num_lost_children() > 0) {
			cut_a_only = false;
			cut_b_only = false;
			cut_c_only = true;
			num_cut-=3;
		}
		else if (T2_b->is_leaf()) {
			if (T2_b->num_lost_children() > 0
				|| T2_b->get_twin()->num_lost_children() > 0) {
			cut_a_only = false;
			cut_b_only = true;
			cut_c_only = false;
			num_cut-=3;
			}
		}
	}
	*/


		Node *node;

		bool cut_a = false;
		bool cut_c = false;
		if (!cut_b_only || T2_a->parent()->get_children().size() > 2) {
			if (!cut_c_only &&
					(!APPROX_EDGE_PROTECTION
					 	|| (!T2_a->is_protected()
							&& (T2_a->parent()->parent() != NULL
								|| !T2_b->is_protected()
								|| T2_a->parent()->get_children().size() > 2)))) {
//					|| cut_a_only)) {
				um.add_event(new CutParent(T1_a));
				T1_a->cut_parent();
				cut_a = true;

				ContractEvent(&um, T1_ac);
				node = T1_ac->contract();
			}
			else
				node = T1_ac;
			if (!cut_a_only &&
					(!APPROX_EDGE_PROTECTION
					 	|| (!T2_c->is_protected()
						&& (T2_c->parent() == NULL
								|| T2_c->parent()->parent() != NULL
								|| !T2_c->get_sibling()->is_protected()
								|| T2_c->parent()->get_children().size() > 2)))) {// &&
//					|| cut_c_only)) {
				um.add_event(new CutParent(T1_c));
				T1_c->cut_parent();
				cut_c = true;

				if (node) {
					ContractEvent(&um, node);
					node = node->contract();
				}
			}

			// contract parents
			// check for T1_ac sibling pair
			if (node && node->is_sibling_pair()){
				um.add_event(new AddToSiblingPairs(sibling_pairs));
				sibling_pairs->push_back(node->lchild());
				sibling_pairs->push_back(node->rchild());
			}
		}

		bool same_component = true;
		if (APPROX_CHECK_COMPONENT && !cut_a_only && !cut_c_only)
			same_component = (T2_a->find_root() == T2_c->find_root());

		Node *T2_ab_parent = T2_ab->parent();
		node = T2_ab;
		if (cut_a) {
			um.add_event(new CutParent(T2_a));
			T2_a->cut_parent();

			//ContractEvent(&um, T2_ab);
			//node = T2_ab->contract();
		}
		bool cut_b = false;
		if (same_component && T2_ab_parent != NULL
				&& !cut_a_only && !cut_c_only
				&& (!APPROX_EDGE_PROTECTION
					|| (!T2_b->is_protected() ))) {
//							&& (T2_b->parentT2_a->parent()->parent() != NULL
//								|| !T2_a->is_protected())))) {
//					|| cut_b_only)) {
			if (multi_node) {
				T2_b = T2_ab;
				um.add_event(new CutParent(T2_ab));
				T2_ab->cut_parent();
				if (T2_a->parent() != NULL) {
					um.add_event(new CutParent(T2_a));
					T2_a->cut_parent();
					um.add_event(new AddChild(T2_a));
					T2_ab_parent->add_child(T2_a);
				}
				else
					node = T2_ab_parent;
			}
			else {
				um.add_event(new CutParent(T2_b));
				T2_b->cut_parent();
				//ContractEvent(&um, node);
				//node = node->contract();
			}
			cut_b = true;
		}
		// T2_b will move up after contraction
		else if (!multi_node) {
			T2_b = T2_b->parent();
		}
		if (node != NULL) {
			ContractEvent(&um, node);
			node = node->contract();
		// check for T2 parents as singletons
		if (node != NULL && node->is_singleton()
				&& node != T2->get_component(0))
			singletons->push_back(node);
		}

		// if T2_c is gone then its replacement is in singleton list
		// contract might delete old T2_c, see where it is
		bool add_T2_c = true;
		T2_c = T1_c->get_twin();
		// ignore T2_c if it is a singleton
		if (cut_c && T2_c != node && T2_c->parent() != NULL) {
			Node *T2_c_parent = T2_c->parent();
			um.add_event(new CutParent(T2_c));
			T2_c->cut_parent();
			ContractEvent(&um, T2_c_parent);
			node = T2_c_parent->contract();
			if (node != NULL && node->is_singleton()
					&& node != T2->get_component(0))
				singletons->push_back(node);
		}
		else {
			add_T2_c = false;
		}

		
		if (cut_a) {
			um.add_event(new AddComponent(T1));
			T1->add_component(T1_a);
			um.add_event(new AddComponent(T2));
			T2->add_component(T2_a);
		}
		if (cut_c) {
			um.add_event(new AddComponent(T1));
			T1->add_component(T1_c);
		}
		if (cut_b) {
			um.add_event(new AddComponent(T2));
			T2->add_component(T2_b);
		}
		// problem if c is deleted
		if (add_T2_c) {
			um.add_event(new AddComponent(T2));
			T2->add_component(T2_c);
		}

		// may have already been added
		if (T2_b->is_leaf() && cut_b)
			singletons->push_back(T2_b);

		num_cut+=3;

		if (cut_a == false && cut_b == false && cut_c == false) {
			num_cut = INT_MAX-3;
		}

	}
}
	}
// if the first component of the forests differ then we have cut p
if (T1->get_component(0)->get_twin() != T2->get_component(0)) {
	if (!T1->contains_rho()) {
		um.add_event(new AddRho(T1));
		um.add_event(new AddRho(T2));
		T1->add_rho();
		T2->add_rho();
	}
	else
		// hack to ignore rho when it shouldn't be in a cluster
		num_cut -=3;
}
if (save_forests) {
	*F1 = new Forest(T1);
	*F2 = new Forest(T2);
}
 
#ifdef DEBUG_APPROX
#ifdef DEBUG_UNDO
 while(um.num_events() > 0) {
		cout << "Undo step " << um.num_events() << endl;
		cout << "T1: ";
		T1->print_components();
		cout << "T2: ";
		T2->print_components();
			cout << "sibling pairs:";
			for (list<Node *>::iterator i = sibling_pairs->begin(); i != sibling_pairs->end(); i++) {
				cout << "  ";
				(*i)->print_subtree_hlpr();
			}
			cout << endl;
	 um.undo();
	cout << endl;
 }
#else
 um.undo_all();
#endif
#else
 um.undo_all();
#endif

 
//		 for(int i = 0; i < T1->num_components(); i++)
//		 	T1->get_component(i)->fix_parents();
//		 for(int i = 0; i < T2->num_components(); i++)
//		 	T2->get_component(i)->fix_parents();
return num_cut;
}

/*******************************************************************************
	RSPR WORSE_3_APPROX_BINARY
*******************************************************************************/

int rSPR_worse_3_approx_binary(Forest *T1, Forest *T2, bool sync) {
	// match up nodes of T1 and T2
	if (sync) {
if (!sync_twins(T1, T2))
	return 0;
	}
	// find sibling pairs of T1
	list<Node *> *sibling_pairs = T1->find_sibling_pairs();
	// find singletons of T2
	list<Node *> singletons = T2->find_singletons();
	list<pair<Forest,Forest> > AFs = list<pair<Forest,Forest> >();

	Forest *F1;
	Forest *F2;

	int ans = rSPR_worse_3_approx_binary_hlpr(T1, T2, &singletons, sibling_pairs, &F1, &F2, true);

	F1->swap(T1);
	F2->swap(T2);
	sync_twins(T1,T2);


	delete sibling_pairs;
	delete F1;
	delete F2;
	return ans;
}

// rSPR_worse_3_approx_binary recursive helper function
int rSPR_worse_3_approx_binary_hlpr(Forest *T1, Forest *T2, list<Node *> *singletons, list<Node *> *sibling_pairs, Forest **F1, Forest **F2, bool save_forests) {
	#ifdef DEBUG_APPROX
cout << "rSPR_worse_3_approx_binary_hlpr" << endl;
			cout << "\tT1: ";
			T1->print_components_with_twins();
			cout << "\tT2: ";
			T2->print_components_with_twins();
			cout << "sibling pairs:";
			for (list<Node *>::iterator i = sibling_pairs->begin(); i != sibling_pairs->end(); i++) {
				cout << "  ";
				(*i)->print_subtree_hlpr();
			}
			cout << endl;
	#endif
	int num_cut = 0;
	UndoMachine um = UndoMachine();
	while(!singletons->empty() || !sibling_pairs->empty()) {
// Case 1 - Remove singletons
while(!singletons->empty()) {
	#ifdef DEBUG_APPROX
		cout << "Case 1" << endl;
	#endif

	Node *T2_a = singletons->back();
	singletons->pop_back();
	// find twin in T1
	Node *T1_a = T2_a->get_twin();
	// if this is in the first component of T_2 then
	// it is not really a singleton.
	// TODO: problem when we cluster and have a singleton as the
	//		first comp of T2
	//    NEED TO MODIFY CUTTING?
	// 		HERE AND IN BB?
	if (T2_a == T2->get_component(0))
		continue;

	Node *T1_a_parent = T1_a->parent();
	if (T1_a_parent == NULL)
		continue;
	bool potential_new_sibling_pair = T1_a_parent->is_sibling_pair();
	// cut the edge above T1_a
	um.add_event(new CutParent(T1_a));
	T1_a->cut_parent();
	um.add_event(new AddComponent(T1));
	T1->add_component(T1_a);
	//if (T1_a->get_sibling_pair_status() > 0)
	//	T1_a->clear_sibling_pair(sibling_pairs);
	//delete(T1_a);

	ContractEvent(&um, T1_a_parent);
	Node *node = T1_a_parent->contract();
	if (node != NULL && potential_new_sibling_pair && node->is_sibling_pair()){
		um.add_event(new AddToFrontSiblingPairs(sibling_pairs));
		sibling_pairs->push_front(node->rchild());
		sibling_pairs->push_front(node->lchild());
	}

	#ifdef DEBUG_APPROX
			cout << "\tT1: ";
			T1->print_components();
			cout << "\tT2: ";
			T2->print_components();
	#endif
}
if(!sibling_pairs->empty()) {
	Node *T1_a = sibling_pairs->back();
	sibling_pairs->pop_back();
	Node *T1_c = sibling_pairs->back();
	sibling_pairs->pop_back();
	um.add_event(new PopSiblingPair(T1_a, T1_c, sibling_pairs));

	//if (T1_a->get_sibling_pair_status() == 0 ||
	//		T1_c->get_sibling_pair_status() == 0) {
	//	continue;
	//}

	//T1_a->clear_sibling_pair_status();
	//T1_c->clear_sibling_pair_status();
	if (T1_a->parent() == NULL || T1_c->parent() == NULL || T1_a->parent() != T1_c->parent()) {
		continue;
	}
	Node *T1_ac = T1_a->parent();
	// lookup in T2 and determine the case
	Node *T2_a = T1_a->get_twin();
	Node *T2_c = T1_c->get_twin();

	#ifdef DEBUG_APPROX
		cout << "Fetching sibling pair" << endl;
		T1_ac->print_subtree();
		cout << "T2_a" << ": ";
		cout << " d=" << T2_a->get_depth() << " ";
		T2_a->print_subtree();
		cout << "T1_c" << ": ";
		T1_c->print_subtree();
		cout << "T2_c" << ": ";
		cout << " d=" << T2_c->get_depth() << " ";
		T2_c->print_subtree();
	#endif

	// Case 2 - Contract identical sibling pair
	if (T2_a->parent() != NULL && T2_a->parent() == T2_c->parent()) {
		#ifdef DEBUG_APPROX
			cout << "Case 2" << endl;
			T1->print_components();
			T2->print_components();
		#endif
		Node *T2_ac = T2_a->parent();
		um.add_event(new ContractSiblingPair(T1_ac));
		um.add_event(new ContractSiblingPair(T2_ac));
		T1_ac->contract_sibling_pair_undoable();
		T2_ac->contract_sibling_pair_undoable();
		um.add_event(new SetTwin(T1_ac));
		um.add_event(new SetTwin(T2_ac));
		T1_ac->set_twin(T2_ac);
		T2_ac->set_twin(T1_ac);
		//T1->add_deleted_node(T1_a);
		//T1->add_deleted_node(T1_c);
		//T2->add_deleted_node(T2_a);
		//T2->add_deleted_node(T2_c);

		// check if T2_ac is a singleton
		if (T2_ac->is_singleton() && !T1_ac->is_singleton() && T2_ac != T2->get_component(0))
			singletons->push_back(T2_ac);
		// check if T1_ac is part of a sibling pair
		if (T1_ac->parent() != NULL && T1_ac->parent()->is_sibling_pair()) {
			um.add_event(new AddToSiblingPairs(sibling_pairs));
			sibling_pairs->push_back(T1_ac->parent()->lchild());
			sibling_pairs->push_back(T1_ac->parent()->rchild());
		}
	}
	// Case 3
	else {
		#ifdef DEBUG_APPROX
			cout << "Case 3" << endl;
		#endif
		
		//  ensure T2_a is below T2_c
		if ((T2_a->get_depth() < T2_c->get_depth()
				&& T2_c->parent() != NULL)
				|| T2_a->parent() == NULL) {
			#ifdef DEBUG_APPROX
				cout << "swapping" << endl;
			#endif
		
			swap(&T1_a, &T1_c);
			swap(&T2_a, &T2_c);

		}
		else if (T2_a->get_depth() == T2_c->get_depth()) {
			if (T2_a->parent() && T2_c->parent() &&
					(T2_a->parent()->get_depth() <
					T2_c->parent()->get_depth())) {
			swap(&T1_a, &T1_c);
			swap(&T2_a, &T2_c);
			}
		}

		// get T2_b
		Node *T2_ab = T2_a->parent();
		Node *T2_b = T2_ab->rchild();
		if (T2_b == T2_a)
			T2_b = T2_ab->lchild();

		#ifdef DEBUG_APPROX
		cout << "T2_b" << ": ";
		cout.flush();
		T2_b->print_subtree();
	#endif
		// cut T1_a, T1_c, T2_a, T2_b, T2_c

		bool cut_b_only = false;
		if (T2_a->parent() != NULL && T2_a->parent()->parent() != NULL && T2_a->parent()->parent() == T2_c->parent()) {
			cut_b_only = true;
			um.add_event(new AddToSiblingPairs(sibling_pairs));
			sibling_pairs->push_back(T1_c);
			sibling_pairs->push_back(T1_a);
		}

		Node *node;

		if (!cut_b_only) {
			um.add_event(new CutParent(T1_a));
			T1_a->cut_parent();

			ContractEvent(&um, T1_ac);
			node = T1_ac->contract();

			um.add_event(new CutParent(T1_c));
			T1_c->cut_parent();


			ContractEvent(&um, node);
			node = node->contract();

			// contract parents
			// check for T1_ac sibling pair
			if (node && node->is_sibling_pair()){
				um.add_event(new AddToSiblingPairs(sibling_pairs));
				sibling_pairs->push_back(node->lchild());
				sibling_pairs->push_back(node->rchild());
			}
		}

		bool same_component = true;
		if (APPROX_CHECK_COMPONENT)
			same_component = (T2_a->find_root() == T2_c->find_root());

		Node *T2_ab_parent = T2_ab->parent();
		node = T2_ab;
		if (!cut_b_only) {
			um.add_event(new CutParent(T2_a));
			T2_a->cut_parent();

			//ContractEvent(&um, T2_ab);
			//node = T2_ab->contract();
		}
		bool cut_b = false;
		if (same_component && T2_ab_parent != NULL) {
			um.add_event(new CutParent(T2_b));
			T2_b->cut_parent();
			//ContractEvent(&um, node);
			//node = node->contract();
			cut_b = true;
		}
		// T2_b will move up after contraction
		else {
			T2_b = T2_b->parent();
		}
			ContractEvent(&um, node);
			node = node->contract();
		// check for T2 parents as singletons
		if (node != NULL && node->is_singleton()
				&& node != T2->get_component(0))
			singletons->push_back(node);

		// if T2_c is gone then its replacement is in singleton list
		// contract might delete old T2_c, see where it is
		bool add_T2_c = true;
		T2_c = T1_c->get_twin();
		// ignore T2_c if it is a singleton
		if (T2_c != node && T2_c->parent() != NULL && !cut_b_only) {

			Node *T2_c_parent = T2_c->parent();
			um.add_event(new CutParent(T2_c));
			T2_c->cut_parent();
			ContractEvent(&um, T2_c_parent);
			node = T2_c_parent->contract();
			if (node != NULL && node->is_singleton()
					&& node != T2->get_component(0))
				singletons->push_back(node);
		}
		else {
			add_T2_c = false;
		}

		
		if (!cut_b_only) {
			um.add_event(new AddComponent(T1));
			T1->add_component(T1_a);
			um.add_event(new AddComponent(T1));
			T1->add_component(T1_c);
			// put T2 cut parts into T2
			um.add_event(new AddComponent(T2));
			T2->add_component(T2_a);
			// may have already been added
		}
		if (cut_b) {
			um.add_event(new AddComponent(T2));
			T2->add_component(T2_b);
		}
		// problem if c is deleted
		if (add_T2_c) {
			um.add_event(new AddComponent(T2));
			T2->add_component(T2_c);
		}

		// may have already been added
		if (T2_b->is_leaf())
			singletons->push_back(T2_b);

		num_cut+=3;

	}
}
	}
// if the first component of the forests differ then we have cut p
if (T1->get_component(0)->get_twin() != T2->get_component(0)) {
	if (!T1->contains_rho()) {
		um.add_event(new AddRho(T1));
		um.add_event(new AddRho(T2));
		T1->add_rho();
		T2->add_rho();
	}
	else
		// hack to ignore rho when it shouldn't be in a cluster
		num_cut -=3;
}
if (save_forests) {
	*F1 = new Forest(T1);
	*F2 = new Forest(T2);
}
 

#ifdef DEBUG_UNDO
 while(um.num_events() > 0) {
		cout << "Undo step " << um.num_events() << endl;
		cout << "T1: ";
		T1->print_components();
		cout << "T2: ";
		T2->print_components();
			cout << "sibling pairs:";
			for (list<Node *>::iterator i = sibling_pairs->begin(); i != sibling_pairs->end(); i++) {
				cout << "  ";
				(*i)->print_subtree_hlpr();
			}
			cout << endl;
	 um.undo();
	cout << endl;
 }
#else
 um.undo_all();
#endif

 
//		 for(int i = 0; i < T1->num_components(); i++)
//		 	T1->get_component(i)->fix_parents();
//		 for(int i = 0; i < T2->num_components(); i++)
//		 	T2->get_component(i)->fix_parents();
return num_cut;
}


int rSPR_branch_and_bound(Forest *T1, Forest *T2) {
	return rSPR_branch_and_bound_range(T1, T2, MAX_SPR);
}


int rSPR_branch_and_bound_range(Forest *T1, Forest *T2, int end_k) {
	string problem_key;
	map<string,ProblemSolution>::iterator i;

	if (MEMOIZE) {
problem_key = T1->str() + ":" + T2->str();
i = memoized_clusters.find(problem_key);
if (i != memoized_clusters.end()) {
	//cout << "already solved: " << endl;
	//cout << problem_key << endl;
	//cout << i->second.T2 << endl;
	//cout << "start" << endl;
	Forest *new_T1 = build_finished_forest(i->second.T1);
	//cout << "middle" << endl;
	Forest *new_T2 = build_finished_forest(i->second.T2);
	//cout << "end" << endl;
	T1->swap(new_T1);
	T2->swap(new_T2);
	sync_twins(T1, T2);
	delete new_T1;
	delete new_T2;
	return i->second.k;
}
	}
	Forest F1 = Forest(T1);
	Forest F2 = Forest(T2);
	int approx_spr = rSPR_worse_3_approx(&F1, &F2);
	int min_spr = approx_spr / 3;
	int exact_spr = rSPR_branch_and_bound_range(T1, T2, min_spr, end_k);
	if (MEMOIZE && exact_spr >= 0 && i == memoized_clusters.end()) {
//string solution_key = T1->str() + ":" + T2->str();
memoized_clusters.insert(make_pair(problem_key,
		ProblemSolution(T1,T2,exact_spr)));
	}

	return exact_spr;
}
	
int rSPR_branch_and_bound_range(Forest *T1, Forest *T2, int start_k,
int end_k) {
	int exact_spr = -1;
	bool in_main = MAIN_CALL;
	MAIN_CALL = false;
	int k;
	for(k = start_k; k <= end_k; k++) {
if (in_main) {
	cout << " " << k;
	cout.flush();
}
//Forest F1 = Forest(T1);
//Forest F2 = Forest(T2);
//exact_spr = rSPR_branch_and_bound(&F1, &F2, k);
exact_spr = rSPR_branch_and_bound(T1,T2, k);
//if (exact_spr >= 0 || k == end_k) {
if (exact_spr >= 0) {
//			F1.swap(T1);
//			F2.swap(T2);
	break;
}
	}
	if (in_main)
cout << endl;
	if (k > end_k)
k = -1;
	return k;
}

int rSPR_branch_and_bound(Forest *T1, Forest *T2, int k) {
	return rSPR_branch_and_bound(T1, T2, k, NULL, NULL, NULL);
}

int rSPR_branch_and_bound(Forest *T1, Forest *T2, int k, list<pair<Forest,Forest>> *extAFs) {
	return rSPR_branch_and_bound(T1, T2, k, NULL, NULL, extAFs);
}

/* rSPR_branch_and_bound
 * Calculate a maximum agreement forest and SPR distance
 * Uses a branch and bound optimization to not explore paths
 * guaranteed to be incorrect based on rspr_3_approx
 * RETURN The rSPR distance
 * NOTE: destructive. The computed forests replace T1 and T2.
 */
int rSPR_branch_and_bound(Forest *T1, Forest *T2, int k,
		map<string, int> *label_map,
		map<int, string> *reverse_label_map,
		list<pair<Forest,Forest>> *extAFs) {
	// find sibling pairs of T1
//	cout << "foo1" << endl;
	if (!sync_twins(T1, T2))
return 0;
	if (PREORDER_SIBLING_PAIRS &&
			T1->get_component(0)->get_preorder_number() == -1) {
		T1->get_component(0)->preorder_number();
		T2->get_component(0)->preorder_number();
}
	if (DEEPEST_PROTECTED_ORDER
			&& T1->get_component(0)->get_edge_pre_start() == -1) {
		T1->get_component(0)->edge_preorder_interval();
		T2->get_component(0)->edge_preorder_interval();
	}

	set<SiblingPair> *sibling_pairs;
	list<Node *> singletons;
	
	list<pair<Forest,Forest> > AFs = list<pair<Forest,Forest> >();
	sibling_pairs = find_sibling_pairs_set(T1);
	singletons = T2->find_singletons();
	list<Node *> protected_stack = list<Node *>();
	int num_ties = 2;


	int final_k = 
rSPR_branch_and_bound_hlpr(T1, T2, k, sibling_pairs, &singletons, false, &AFs, &protected_stack, &num_ties);

//		cout << "foo" << endl;
	// TODO: this is a cheap hack
	if (!AFs.empty()) {
if (ALL_MAFS
#ifdef DEBUG
		|| true
#endif
		) {

	vector<string> vecComponents;
	unordered_set<string> setMAFs;
	cout << endl << endl << "FOUND ANSWERS " << AFs.size() << endl;
	// TODO: this is a cheap hack
	for (list<pair<Forest,Forest> >::iterator x = AFs.begin(); x != AFs.end(); x++) {
		if (label_map != NULL && reverse_label_map != NULL) {
			x->first.numbers_to_labels(reverse_label_map);
			x->second.numbers_to_labels(reverse_label_map);
		}
		
		string strMAF = (x->first).add_vec_components(&vecComponents);
		std::pair<std::unordered_set<string>::iterator, bool> insertResult = setMAFs.insert(strMAF);
		if(ALL_MERGED_MAFS && insertResult.second){
			extAFs->push_back(make_pair(x->first, x->second));
		}

		cout << "\tT1: ";
		x->first.print_components();
		cout << "\tT2: ";
		x->second.print_components();
		if (label_map != NULL && reverse_label_map != NULL) {
			x->first.labels_to_numbers(label_map, reverse_label_map);
			x->second.labels_to_numbers(label_map, reverse_label_map);
		}
	}
	cout << "FOUND UNIQUE ANSWERS " << setMAFs.size() << endl;
}

	if(SHOW_PERCENT_LGT_EVENTS){
		map<int, int> cut_edge_count;
		int totalAFs = AFs.size();
		for (list<pair<Forest,Forest> >::iterator x = AFs.begin(); x != AFs.end(); x++) {
			vector<Node *> components = x->second.components;
			vector<Node *>::iterator it;
			for(it = components.begin(); it != components.end(); it++) {
				Node *root = *it;
				if (root != NULL &&  root->get_edge_pre_start() != 0){
					cut_edge_count[root->get_preorder_number()]++;
				}
			}
		}

		std::map<int, int>::iterator itr = cut_edge_count.begin();
		while (itr != cut_edge_count.end())
		{
			int count = itr->second;
			float percent = (count * 100) / totalAFs;
			std::cout << "Key: " << itr->first << ", Value: " << percent << std::endl;
			++itr;
		}
	}

AFs.front().first.swap(T1);
AFs.front().second.swap(T2);
sync_twins(T1,T2);
	}
	if (final_k >= 0)
final_k = k - final_k;
	delete sibling_pairs;
	return final_k;
}

void add_sibling_pair(set<SiblingPair> *sibling_pairs, Node *a, Node *c, UndoMachine *um) {
	SiblingPair sp = SiblingPair(a,c);
	pair< set<SiblingPair>::iterator, bool> ins = 
	sibling_pairs->insert(sp);
	if (ins.second == false) {
um->add_event(new RemoveSetSiblingPairs(sibling_pairs, *(ins.first)));
sibling_pairs->erase(ins.first);
ins = sibling_pairs->insert(sp);
	}
	um->add_event(new AddToSetSiblingPairs(sibling_pairs, *(ins.first)));
}

SiblingPair pop_sibling_pair(set<SiblingPair> *sibling_pairs, UndoMachine *um) {
	set<SiblingPair>::iterator s = sibling_pairs->begin();
	SiblingPair spair = SiblingPair(*s); 
	um->add_event(new RemoveSetSiblingPairs(sibling_pairs, spair));
	sibling_pairs->erase(s);
	return spair;
}

SiblingPair pop_sibling_pair(set<SiblingPair>::iterator s, set<SiblingPair> *sibling_pairs, UndoMachine *um) {
	SiblingPair spair = SiblingPair(*s); 
	um->add_event(new RemoveSetSiblingPairs(sibling_pairs, spair));
	sibling_pairs->erase(s);
	return spair;
}

inline int rSPR_branch_and_bound_hlpr(Forest *T1, Forest *T2, int k,
set<SiblingPair> *sibling_pairs, list<Node *> *singletons,
bool cut_b_only, list<pair<Forest,Forest> > *AFs,
list<Node *> *protected_stack, int *num_ties) {
	return rSPR_branch_and_bound_hlpr(T1, T2, k, sibling_pairs,
			singletons, cut_b_only, AFs, protected_stack, num_ties, NULL, NULL);
}

// rSPR_branch_and_bound recursive helper function
int rSPR_branch_and_bound_hlpr(Forest *T1, Forest *T2, int k,
set<SiblingPair> *sibling_pairs, list<Node *> *singletons,
bool cut_b_only, list<pair<Forest,Forest> > *AFs,
list<Node *> *protected_stack, int *num_ties, Node *prev_T1_a, Node *prev_T1_c) {
	cout << "rSPR_branch_and_bound_hlpr()" << endl;
	#ifdef DEBUG
	cout << "\tT1: ";
	T1->print_components();
	cout << "\tT2: ";
	T2->print_components();
	cout << "K=" << k << endl;
	cout << "sibling pairs:";
	for (set<SiblingPair>::iterator i = sibling_pairs->begin(); i != sibling_pairs->end(); i++) {
cout << "  ";
(*i).a->print_subtree_hlpr();
cout << ",";
(*i).c->print_subtree_hlpr();
	}
	cout << endl;
	cout << "protected_stack:";
	for (list<Node *>::iterator i = protected_stack->begin(); i != protected_stack->end(); i++) {
cout << "  ";
(*i)->print_subtree_hlpr();
	}
	cout << endl;
	#endif
	UndoMachine um = UndoMachine();
	#ifdef DEBUG_LGT_EVENTS
		cout << "T1111" << endl;
		T1->print_components();
		cout << "T2222" << endl;
		T2->print_components();
	#endif

	while(!singletons->empty() || !sibling_pairs->empty()) {
		// Case 1 - Remove singletons
		while(!singletons->empty()) {
			Node *T2_a = singletons->back();
			#ifdef DEBUG
				cout << "Case 1" << endl;
				cout << "a " << T2_a->str_subtree() << endl;
			#endif
		
			singletons->pop_back();
			// find twin in T1
			Node *T1_a = T2_a->get_twin();
		
			//if (T1_a->get_sibling_pair_status() > 0)
			//	cout << T1_a->get_sibling(sibling_pairs) << endl;
			// if this is in the first component of T_2 then
			// it is not really a singleton.
			Node *T1_a_parent = T1_a->parent();
			if (T1_a_parent == NULL)
				continue;
			bool potential_new_sibling_pair = T1_a_parent->is_sibling_pair();
			if (T2_a == T2->get_component(0)) {
				// TODO: should we do this when it happens?
				if (!T1->contains_rho()) {
					um.add_event(new AddRho(T1));
					um.add_event(new AddRho(T2));
					T1->add_rho();
					T2->add_rho();
					k--;
					#ifdef DEBUG
					cout << "adding p element, k=" << k << endl;
					#endif
				}
			}
		
			// cut the edge above T1_a
			um.add_event(new CutParent(T1_a));
			T1_a->cut_parent();
		
			um.add_event(new AddComponent(T1));
			T1->add_component(T1_a);
			ContractEvent(&um, T1_a_parent);
		
		
			Node *node = T1_a_parent->contract();
		
			if (node != NULL && potential_new_sibling_pair && node->is_sibling_pair()){
				add_sibling_pair(sibling_pairs, node->lchild(), node->rchild(),
						&um);
			}
			#ifdef DEBUG
				cout << "\tT1: ";
				T1->print_components();
				cout << "\tT2: ";
				T2->print_components();
			#endif
		//			sp_i = sibling_pairs->begin();
		}
		if(!sibling_pairs->empty()) {
			Node *T1_a;
			Node *T1_c;
			set<SiblingPair>::iterator deepest_valid = sibling_pairs->end();
			int deepest_depth = INT_MAX;
			int deepest_depth_2 = INT_MAX;
			Node *best_a = NULL;
			Node *best_c = NULL;
			/* pop protected_stack when out of order sibling pairs
			 have already contracted it */
			if(!protected_stack->empty()
					&& protected_stack->back()->get_twin()->parent() == NULL
					&& protected_stack->back()->get_twin() != T1->get_component(0)) {
				um.undo_all();
				return -1;
			}
			while(!protected_stack->empty()
					&& (protected_stack->back()->is_contracted()
					// this shouldn't happen
						|| protected_stack->back()->get_twin()->parent() == NULL)) {
				um.add_event(new ListPopBack(protected_stack));
				protected_stack->pop_back();
			}
			if (LEAF_REDUCTION && !cut_b_only) {
				bool found = false;
				set<SiblingPair>::iterator sp_i = sibling_pairs->begin();
				// correct in case sibling pair involves previous
				/*				if (sp_i != sibling_pairs->begin()) {
					if (check_all_pairs)
						sp_i = sibling_pairs->begin();
					else
						sp_i--;
				}
				*/
				while (sp_i != sibling_pairs->end()) {
					T1_a = (*sp_i).a;
					T1_c = (*sp_i).c;
					if (T1_a->parent() == NULL || T1_a->parent() != T1_c->parent()) {
						um.add_event(new RemoveSetSiblingPairs(sibling_pairs,
									SiblingPair(T1_a, T1_c)));
						set<SiblingPair>::iterator rem = sp_i;
						sp_i++;
						sibling_pairs->erase(rem);
						continue;
					}
					Node *T2_a = T1_a->get_twin();
					Node *T2_c = T1_c->get_twin();
					// select if this is a Case 2 or, optionally, nonbranching
					if ((T2_a->parent() != NULL && T2_a->parent() == T2_c->parent())
							|| (!cut_b_only && PREFER_NONBRANCHING
									&& is_nonbranching(T1, T2, T1_a, T1_c, T2_a, T2_c))) {
						um.add_event(new RemoveSetSiblingPairs(sibling_pairs,
									SiblingPair(T1_a, T1_c)));
						set<SiblingPair>::iterator rem = sp_i;
						sp_i++;
						sibling_pairs->erase(rem);
						found = true;
						break;
					}
					if (DEEPEST_ORDER) {
						int depth;
						int depth2;
						if (T1_a->get_depth() < T1_c->get_depth())
							depth = T1_c->get_depth();
						else
							depth = T1_a->get_depth();
						if (T2_a->get_depth() < T2_c->get_depth())
							depth2 = T2_c->get_depth();
						else
							depth2 = T2_a->get_depth();
						if (deepest_valid == sibling_pairs->end()
								|| deepest_depth < depth
								|| (deepest_depth == depth && deepest_depth_2 < depth2) ) {
							// TODO: this crashes on bigtest2
							// Why can we end up cutting the protected node?
							if (!DEEPEST_PROTECTED_ORDER
									|| protected_stack->empty()
									|| (protected_stack->back()->get_twin()->parent()->get_edge_pre_start()
											<= T1_a->get_preorder_number()
										&& protected_stack->back()->get_twin()->parent()->get_edge_pre_end()
											>= T1_a->get_preorder_number())) {
								deepest_valid = sp_i;
								deepest_depth = depth;
								deepest_depth_2 = depth2;
							}
						}
					}
					/* TODO: create a stack of protected nodes (intervals?) and only
						 accept the deepest sibling pair within the interval
					*/

					/* TODO: remember to pop the stack when we include the protected
					   node */
					sp_i++;
				}
				if (!found) {
					if (sibling_pairs->empty())
						continue;
					else {
						SiblingPair spair;
//						cout << "depth: " << deepest_depth << endl;
						if (DEEPEST_ORDER && deepest_valid != sibling_pairs->end())
							spair = pop_sibling_pair(deepest_valid, sibling_pairs, &um);
						else
							spair = pop_sibling_pair(sibling_pairs, &um);
						T1_a = spair.a;
						T1_c = spair.c;
					}
				}
			}
			else {
				if (prev_T1_a != NULL && prev_T1_c != NULL) {
					T1_a = prev_T1_a;
					T1_c = prev_T1_c;
					prev_T1_a = NULL;
					prev_T1_c = NULL;
				}
				else {
					SiblingPair spair = pop_sibling_pair(sibling_pairs, &um);
					T1_a = spair.a;
					T1_c = spair.c;
				}
			}
//			if (T1_a->parent() != NULL)
//				cout << "a_p: " << T1_a->parent()->str_subtree() << endl;
//			if (T1_c->parent() != NULL)
//				cout << "c_p: " << T1_c->parent()->str_subtree() << endl;
			if (T1_a->parent() == NULL || T1_a->parent() != T1_c->parent()) {
				continue;
			}
			if (!T1_a->can_be_sibling() || !T1_c->can_be_sibling()) {
				continue;
			}
			Node *T1_ac = T1_a->parent();
			// lookup in T2 and determine the case
			Node *T2_a = T1_a->get_twin();
			Node *T2_c = T1_c->get_twin();

			if (T2_a->parent() != NULL && T2_a->parent() == T2_c->parent()) {
				#ifdef DEBUG
					cout << "Case 2" << endl;
					T1_ac->print_subtree();
				#endif
				Node *T2_ac = T2_a->parent();

				if (CHECK_MERGE_DEPTH &&
						(T2_a->get_max_merge_depth() > T2_ac->get_depth()
							|| T2_c->get_max_merge_depth() > T2_ac->get_depth())) {
					um.undo_all();
					return -1;
				}

				if (!protected_stack->empty() &&
						(T2_a == protected_stack->back()
						 	|| T2_c == protected_stack->back())) {
					um.add_event(new ListPopBack(protected_stack));
					protected_stack->pop_back();
				}
				// CAN THIS HAPPEN TWICE?
				if (!protected_stack->empty() &&
						(T2_a == protected_stack->back()
						 	|| T2_c == protected_stack->back())) {
					um.add_event(new ListPopBack(protected_stack));
					protected_stack->pop_back();
				}


				um.add_event(new ContractSiblingPair(T1_ac));
				T1_ac->contract_sibling_pair_undoable();
				um.add_event(new ContractSiblingPair(T2_ac, T2_a, T2_c, &um));
				Node *T2_ac_new = T2_ac->contract_sibling_pair_undoable(T2_a, T2_c);
				if (T2_ac_new != NULL && T2_ac_new != T2_ac) {
					T2_ac = T2_ac_new;
					um.add_event(new CreateNode(T2_ac));
					um.add_event(new ContractSiblingPair(T2_ac));
					T2_ac->contract_sibling_pair_undoable();
				}

				um.add_event(new SetTwin(T1_ac));
				um.add_event(new SetTwin(T2_ac));
				T1_ac->set_twin(T2_ac);
				T2_ac->set_twin(T1_ac);
				//T1->add_deleted_node(T1_a);
				//T1->add_deleted_node(T1_c);
				//T2->add_deleted_node(T2_a);
				//T2->add_deleted_node(T2_c);

				// check if T2_ac is a singleton
				if (T2_ac->is_singleton() && !T1_ac->is_singleton() && T2_ac != T2->get_component(0))
					singletons->push_back(T2_ac);
				// check if T1_ac is part of a sibling pair
				if (T1_ac->parent() != NULL && T1_ac->parent()->is_sibling_pair()) {
				add_sibling_pair(sibling_pairs, T1_ac->parent()->lchild(), T1_ac->parent()->rchild(),
						&um);
				}
				#ifdef DEBUG
					cout << "\tT1: ";
					T1->print_components();
					cout << "\tT2: ";
					T2->print_components();
				#endif
			}
			/* need to copy trees and lists for branching
			 * use forest copy constructor for T1 and T2 giving T1' and T2'
			 * T1' twins are in T2, and same for T2' and T1.
			 * singleton list will be empty except for maybe above the cut,
			 * so this can be created.
			 * fix one set of twins (T2->T1' or T1->T2' not sure)
			 * exploit chained twin relationship to copy sibling pair list
			 * fix other set of twins
			 * swap T2 and T2' root nodes
			 * now do the cut
			 *
			 * note: don't copy for 3rd cut, is a waste
			 */

			// Case 3
			// note: guaranteed that singleton list is empty
			else {
				if (k <= 0) {
					if ((!CUT_LOST || k < 0 ||
								(T1_a->num_lost_children() == 0 &&
								 T1_c->num_lost_children() == 0))
							&& ((T2_c->parent() != NULL && T2_a->parent() != NULL)|| !T2->contains_rho())) {
						singletons->clear();
						um.undo_all();
						return k-1;
					}
				}
				Forest *best_T1;
				Forest *best_T2;
				int best_k = -1;
				int answer_a = -1;
				int answer_b = -1;
				int answer_c = -1;
				bool cut_ab_only = false;
				bool cut_a_only = false;
				bool cut_c_only = false;
				bool cut_a_or_merge_ac = false;
				bool same_component = true;
				int lca_depth = -1;
				int path_length = -1;
				bool cut_b_only_if_not_a_or_c = false;
				bool cob = false;
				int undo_state = um.num_events();
				//  ensure T2_a is below T2_c
				if ((T2_a->get_depth() < T2_c->get_depth()
						&& T2_c->parent() != NULL)
						|| T2_a->parent() == NULL) {
					swap(&T1_a, &T1_c);
					swap(&T2_a, &T2_c);
				}
				else if (T2_a->get_depth() == T2_c->get_depth()) {
					if (T2_a->parent() && T2_c->parent() &&
							(T2_a->parent()->get_depth() <
							T2_c->parent()->get_depth()
							//|| (T2_c->parent()->parent()
							//&& T2_c->parent()->parent() == T2_a->parent())
							)) {
					swap(&T1_a, &T1_c);
					swap(&T2_a, &T2_c);
					}
				}
				Node *T2_b = T2_a->parent()->rchild();
				if (T2_b == T2_a)
					T2_b = T2_a->parent()->lchild();
				bool multi_node = false;
				if (T2_a->parent()->get_children().size() > 2)
					multi_node = true;

			if (CUT_ONE_B) {
				if (T2_a->parent()->parent() == T2_c->parent()
					&& T2_c->parent() != NULL && !cut_b_only)
					cut_b_only=true;
					cob = true;
			}
			else if (CUT_ONE_AB) {
				if (T2_a->parent()->parent() == T2_c->parent()
					&& T2_c->parent() != NULL)
					cut_ab_only=true;
			}
			if (CUT_TWO_B && !cut_b_only && T1_ac->parent() != NULL) {
				Node *T1_s = T1_ac->get_sibling();
				if (T1_s->is_leaf()) {
					Node *T2_l = T2_a->parent()->parent();
					// Note: is this too harsh? If T2_l is nonbinary then can we do cut_b_only_if_not_a_or_c ?
					if (T2_l != NULL && T2_l->get_children().size() <= 2) {
						if (T2_c->parent() != NULL && T2_c->parent()->parent() == T2_l
								&& ((T2_a->parent()->get_children().size() <= 2
										&& T2_c->parent()->get_children().size() <= 2)
									|| T1_s->get_twin()->is_protected())) {
							if (T2_l->get_sibling() == T1_s->get_twin()) {
								cut_b_only=true;
							}
							else if (T2_l->parent() == NULL &&
									(T2->contains_rho() ||
									 T2->get_component(0) != T2_l)) {
								cut_b_only_if_not_a_or_c=true;
							}
						}
						else if ((T2_l = T2_l->parent()) != NULL
								&& T2_c->parent() == T2_l
								&& ((T2_a->parent()->get_children().size() <= 2
										&& T2_a->parent()->parent()->get_children().size() <= 2
										&& T2_l->get_children().size() <= 2)
									|| T1_s->get_twin()->is_protected())) {
							if (T2_l->get_sibling() == T1_s->get_twin()) {
								cut_b_only=true;
							}
							else if (T2_l->parent() == NULL &&
									(T2->contains_rho() ||
									 T2->get_component(0) != T2_l)) {
								cut_b_only_if_not_a_or_c=true;
							}
						}
					}
				}
			}
			if (REVERSE_CUT_ONE_B && (!cut_b_only || (cob && multi_node)) &&
					T1_ac->parent() != NULL) {
				Node *T1_s = T1_ac->get_sibling();
				if (T1_s->is_leaf()) {
					Node *T2_s = T1_s->get_twin();
					if (T2_s->parent() == T2_a->parent()) {
						cut_c_only=true;
						cut_b_only=false;
						cob=false;
					}
					else if (T2_s->parent() == T2_c->parent()) {
						if (T2_c->parent()->get_children().size() <= 2) {
							cut_a_only=true;
							cut_b_only=false;
							cob=false;
						}
						else {
							cut_a_or_merge_ac=true;
							cut_b_only=false;
							cob=false;
						}
					}
					else if (REVERSE_CUT_ONE_B_3
							// TODO: there is a chance for an additional optimization
							// here. If T2_s is not protected then we can cut c or
							// have to cut a (and s?)
//							&& T2_c->parent()->parent() != NULL
							// TODO: buggy? Do we also have to consider cutting b_1 through the other b's?
							&& T2_s->is_protected()
							&& T2_s->parent() != NULL
							&& T2_s->parent()->parent() == T2_a->parent()
							&& T2_s->parent()->get_children().size() <= 2) {
						//cut_c_only = true;
						cut_b_only=false;
						cob=false;
						if (!T2_a->is_protected()) {
							um.add_event(new ProtectEdge(T2_a));
							T2_a->protect_edge();
						}
					}
				}
				else if (REVERSE_CUT_ONE_B_2 && T2_c->parent() != NULL
						&& chain_match(T1_s, T2_c->get_sibling(), T2_a)) {
					cut_a_only = true;
					cut_b_only=false;
					cob=false;
				}
			}
			if (REVERSE_CUT_ONE_B_3 && (!cut_b_only || (cob && multi_node)) &&
					T1_ac->parent() != NULL && T1_ac->parent()->parent() != NULL) {
				Node *T1_s = T1_ac->parent()->get_sibling();
				Node *T2_s = T1_s->get_sibling();
				if (T1_s->is_leaf()
						&& T2_s->is_protected()
						&& T2_s->parent() != NULL
						&& T2_s->parent()->parent() == T2_a->parent()
						&& T2_s->parent()->get_children().size() <= 2) {
					//cut_c_only = true;
					cut_b_only=false;
					cob=false;
					if (!T2_a->is_protected()) {
						um.add_event(new ProtectEdge(T2_a));
						T2_a->protect_edge();
					}
				}
			}
			if (CUT_TWO_B_ROOT && cut_a_only == false && cut_c_only == false
					&& cut_b_only_if_not_a_or_c == true) {
				cut_b_only = true;
			}
/*			if (CUT_LOST) {
				if (T1_a->num_lost_children() > 0
						|| T2_a->num_lost_children() > 0) {
					cut_a_only = true;
					cut_b_only = false;
					cut_c_only = false;
					k++;
				}
				else if (T1_c->num_lost_children() > 0
						|| T2_c->num_lost_children() > 0) {
					cut_a_only = false;
					cut_b_only = false;
					cut_c_only = true;
					k++;
				}
				else if (T2_b->is_leaf()) {
					if (T2_b->num_lost_children() > 0
						|| T2_b->get_twin()->num_lost_children() > 0) {
					cut_a_only = false;
					cut_b_only = true;
					cut_c_only = false;
					k++;
					}
				}
			}
*/

			#ifdef DEBUG
					cout << "Case 3" << endl;
					cout << "\tT1: ";
					T1->print_components();
					cout << "\tT2: ";
					T2->print_components();
					cout << "\tK=" << k << endl;
					cout << "\tsibling pairs:";
					for (set<SiblingPair>::iterator i = sibling_pairs->begin(); i != sibling_pairs->end(); i++) {
						cout << "  ";
						(*i).a->print_subtree_hlpr();
						cout << ",";
						(*i).c->print_subtree_hlpr();
					}
					cout << endl;
					cout << "\tprotected_stack:";
					for (list<Node *>::iterator i = protected_stack->begin(); i != protected_stack->end(); i++) {
						cout << "  ";
						(*i)->print_subtree_hlpr();
					}
					cout << endl;
					cout << "\tcut_a_only=" << cut_a_only << endl;
					cout << "\tcut_b_only=" << cut_b_only << endl;
					cout << "\tcut_c_only=" << cut_c_only << endl;
					cout << "\tT2_a " << T2_a->str() << " "
						<< T2_a->get_depth() << endl;
					cout << "\tT2_c " << T2_c->str() << " "
						<< T2_c->get_depth() << endl;
					cout << "\tT2_b " << T2_b->str_subtree() << " "
						<< T2_b->get_depth() << endl;
				#endif
			


				// copy elements
					/*
				Forest *T1_copy;
				Forest *T2_copy;
				list<Node *> *sibling_pairs_copy;
				Node *T1_a_copy;
				Node *T1_c_copy;
				Node *T2_a_copy;
				Node *T2_c_copy;
				*/
				//list<Node *> *singletons_copy = new list<Node *>();

				// make copies for the approx
				// be careful we do not kill real T1 and T2
				// ie use the copies
				if (BB && !cut_a_only && !cut_b_only && !cut_c_only) {
					list<Node *> *spairs;
					spairs = new list<Node *>();
					spairs->push_back(T1_c);
					spairs->push_back(T1_a);
					for (set<SiblingPair>::iterator i = sibling_pairs->begin(); i != sibling_pairs->end(); i++) {
						spairs->push_back((*i).a);
						spairs->push_back((*i).c);
					}
					int approx_spr = rSPR_worse_3_approx_hlpr(T1, T2,
							singletons, spairs, NULL, NULL, false);
					delete spairs;
					#ifdef DEBUG
						cout << "\tT1: ";
						T1->print_components();
						cout << "\tT2: ";
						T2->print_components();
						cout << "approx =" << approx_spr << endl;
					#endif
					if (approx_spr  >  3*k){
						#ifdef DEBUG
							cout << "approx failed" << endl;
						#endif
						um.undo_all();
						return -1;
					}
					um.undo_to(undo_state);
				}

				if (!cut_a_only && !cut_c_only && !cut_b_only)
					same_component = T2_a->same_component(T2_c, lca_depth, path_length);
				if (cut_b_only)
					same_component = true;
				
			if (CLUSTER_REDUCTION && (MAX_CLUSTERS < 0 || NUM_CLUSTERS < MAX_CLUSTERS)) {
				// clean up singletons
				// TODO: this is duplication
				/*
			while(!singletons->empty()) {
				Node *T2_a = singletons->back();
				singletons->pop_back();
				// find twin in T1
				Node *T1_a = T2_a->get_twin();
				// if this is in the first component of T_2 then
				// it is not really a singleton.
				Node *T1_a_parent = T1_a->parent();
				if (T1_a_parent == NULL)
					continue;
				bool potential_new_sibling_pair = T1_a_parent->is_sibling_pair();
				if (T2_a == T2->get_component(0)) {
					T1->add_rho();
					T2->add_rho();
					k--;
				}
	
				// cut the edge above T1_a
				T1_a->cut_parent();
				T1->add_component(T1_a);
				Node *node = T1_a_parent->contract();
				if (node != NULL && potential_new_sibling_pair && node->is_sibling_pair()){
					sibling_pairs->push_front(node->lchild());
					sibling_pairs->push_front(node->rchild());
				}
			}
			*/
	//			cout << "foo" << endl;
	//			cout << "foo2" << endl;
//				cout << "\tT1: ";
//				T1->print_components();
//				cout << "\tT2: ";
//				T2->print_components();
				sync_interior_twins_real(T1, T2);
				list<Node *> *cluster_points = find_cluster_points(T1, T2);
				//cluster_points->erase(++cluster_points->begin(),cluster_points->end());
	
				// TODO: could this be faster by using the approx to allocate
				// a certain amount of the k to different clusters?
				// TODO: write pseudocode for what we need
				// TODO: then implement it
				// NOTE: need to make a list of ClusterInstances and then
				// solve each.
				// TODO: where should we do this? Just before we would
				// normally branch?
	
//				cout << "k=" << k << endl;
//				cout << "cp=" << cluster_points->size() << endl;
				if (!cluster_points->empty()) {
					NUM_CLUSTERS++;
					sibling_pairs->clear();
#ifdef DEBUG_CLUSTERS
					cout << "CLUSTERS" << endl;
					for(int j = 0; j < 70; j++) {
						cout << "*";
					}
					cout << endl;
					for(list<Node *>::iterator i = cluster_points->begin();
							i != cluster_points->end(); i++) {
						cout << (*i)->str_subtree() << endl;
						cout << (*i)->get_twin()->str_subtree() << endl;
						for(int j = 0; j < 70; j++) {
							cout << "*";
						}
						cout << endl;
					}
					cout << endl;
#endif
	
					list<ClusterInstance> clusters =
						cluster_reduction(T1, T2, cluster_points);
	
					// TODO: make it so we don't need this?
					T1->unsync_interior();
					T2->unsync_interior();
					while(!clusters.empty()) {
						ClusterInstance cluster = clusters.front();
						clusters.pop_front();
						cluster.F1->unsync_interior();
						cluster.F2->unsync_interior();
#ifdef DEBUG_CLUSTERS
						cout << "CLUSTER_START" << endl;
						cout << &(*cluster.F1->get_component(0)) << endl;
						if (!clusters.empty())
							cout << &(*clusters.front().F1->get_component(0)) << endl;
						cout << "\tF1: ";
						cluster.F1->print_components();
						cout << "\tF2: ";
						cluster.F2->print_components();
						cout << "K=" << k << endl;
#endif
						int cluster_spr = -1;
						if (k >= 0) {
							// hack for clusters with no rho
//							cout << __LINE__ << endl;
//							cout << cluster.F2_cluster_node << endl;
//							cout << cluster.F2_has_component_zero << endl;
							if ((cluster.F2_cluster_node == NULL
										|| (cluster.F2_cluster_node->is_leaf()
												&& cluster.F2_cluster_node->parent() == NULL
												&& cluster.F2_cluster_node->
												get_num_clustered_children() <= 1
												&& (cluster.F2_cluster_node !=
														cluster.F2_cluster_node->get_forest()->
															get_component(0))))
										&& cluster.F2_has_component_zero == false) {
//							cout << __LINE__ << endl;
								cluster.F1->add_rho();
								cluster.F2->add_rho();
							}
							cluster_spr = rSPR_branch_and_bound_range(cluster.F1,
									cluster.F2, k);
							if (cluster_spr >= 0) {
//							cout << "cluster k=" << cluster_spr << endl;
//							cout << "\tF1: ";
//							cluster.F1->print_components();
//							cout << "\tF2: ";
//							cluster.F2->print_components();
								k -= cluster_spr;
							}
							else {
								k = -1;
							}
						}
						if (k > -1) {
//							cout << "\tF1: ";
//							T1->print_components();
//							cout << "\tF2: ";
//							T2->print_components();
								if (!cluster.is_original()) {
									int adjustment = cluster.join_cluster(T1, T2);
									k -= adjustment;
									delete cluster.F1;
									delete cluster.F2;
	
	//						cout << cluster.F1_cluster_node->str_subtree() << endl;
	//						cout << cluster.F2_cluster_node->str_subtree() << endl;
	//						Node *p = cluster.F1_cluster_node;
	//						while (p->parent() != NULL)
	//							p = p->parent();
	//						cout << &(*p) << endl;
	//						cout << p->str_subtree() << endl;
	
	
							//cout << "\tF1: ";
							//T1->print_components();
							//cout << "\tF2: ";
							//T2->print_components();
								}
//						else {
//							cout << "original" << endl;
//						}
							}
							else {
								if (!cluster.is_original()) {
									//if (cluster.F1_cluster_node != NULL)
									//	cluster.F1_cluster_node->contract();
									//if (cluster.F2_cluster_node != NULL)
									//	cluster.F2_cluster_node->contract();
									delete cluster.F1;
									delete cluster.F2;
								}
							}
					}
					delete cluster_points;
//					cout << "returning k=" << k << endl;
					NUM_CLUSTERS--;
					return k;
				}
				else {
					T1->unsync_interior();
					T2->unsync_interior();
					// HACK to allow only initial clusters
					// TODO: use UndoMachine in this clustering section
					// and update this clustering to not require copying
					NUM_CLUSTERS++;
				}
				delete cluster_points;
	
	//			cout << "done" << endl;
			}
				 // make copies for the branching
			/*
				copy_trees(&T1, &T2, &sibling_pairs, &T1_a, &T1_c, &T2_a, &T2_c,
						&T1_copy, &T2_copy, &sibling_pairs_copy,
					&T1_a_copy, &T1_c_copy, &T2_a_copy, &T2_c_copy);
					*/

				Node *node;
				Node *T2_ab = T2_a->parent();

				bool balanced = false;
				bool multi_b1 = T2_ab->get_children().size() > 2;
				bool multi_b2 = false;
				Node *T2_d = NULL;
				if (path_length == 4) {
					if (T2_a->parent()->parent() == T2_c->parent()->parent())
						balanced = true;
					if (balanced)
						multi_b2 = T2_c->parent()->get_children().size() > 2;
					else
						multi_b2 = T2_ab->parent()->get_children().size() > 2;
					if (balanced)
						T2_d = T2_c->get_sibling();
				}


				rspr_branch_and_bound_cut_a_hlpr(T1, T2, k, sibling_pairs,
					singletons, AFs, protected_stack, num_ties, T1_c, T2_a, T2_b, T2_c,
					cut_a_only, cut_b_only, cut_c_only, path_length, &um, 
					T2_ab, balanced, multi_b1, multi_b2, T2_d, best_k);

				best_T1 = T1;
				best_T2 = T2;

				um.undo_to(undo_state);

/*				#ifdef DEBUG
					cout << "Case 3 CHECK" << endl;
					cout << "\tT1: ";
					T1->print_components();
					cout << "\tT2: ";
					T2->print_components();
					cout << "\tK=" << k << endl;
					cout << "\tsibling pairs:";
					for (list<Node *>::iterator i = sibling_pairs->begin(); i != sibling_pairs->end(); i++) {
						cout << "  ";
						(*i)->print_subtree_hlpr();
					}
					cout << endl;
					cout << "\tcut_b_only=" << cut_b_only << endl;
					cout << "\tT2_a " << T2_a->str() << " "
						<< T2_a->get_depth() << endl;
					cout << "\tT2_c " << T2_c->str() << " "
						<< T2_c->get_depth() << endl;
					cout << "\tT2_b " << T2_b->str_subtree() << " "
						<< T2_b->get_depth() << endl;
				#endif
*/

				//load the copy
				/*
				T1 = T1_copy;
				T2 = T2_copy;
				T1_a = T1_a_copy;
				T1_c = T1_c_copy;
				T2_a = T2_a_copy;
				T2_c = T2_c_copy;
				sibling_pairs = sibling_pairs_copy;
				singletons = new list<Node *>();
				*/


				// make copies for the branching
				/*
				copy_trees(&T1, &T2, &sibling_pairs, &T1_a, &T1_c, &T2_a, &T2_c,
						&T1_copy, &T2_copy, &sibling_pairs_copy,
						&T1_a_copy, &T1_c_copy, &T2_a_copy, &T2_c_copy);
						*/


				// get T2_b
				T2_ab = T2_a->parent();
				T2_b = T2_ab->rchild();

				if (T2_b == T2_a)
					T2_b = T2_ab->lchild();


//				if ((!CUT_AC_SEPARATE_COMPONENTS ||
//							T2_a->find_root() == T2_c->find_root())
//						&& (!multi_node && T2_b->is_protected()))
//					cout << "protected k=" << k << endl;

				// BUG: we need to cut b or c when b is a multifurcating COB
				// TODO: if the LCA of B_1's leaves in T1 is not an ancestor of 
				// a then this is safe
				if (multi_node && cob) {
					vector<Node *> B_1 = T2_a->parent()->find_leaves();
					LCA T1_LCA = LCA(T1->get_component(0));
					Node *B_1_lca = NULL;
					for(int i = 0; i < B_1.size(); i++) {
						if (B_1[i] == T2_a) {
							continue;
						}
						if (B_1_lca == NULL) {
							B_1_lca = B_1[i]->get_twin();
						}
						else {
							B_1_lca = T1_LCA.get_lca(B_1_lca, B_1[i]->get_twin());
						}
					}
					Node *T1_a_ancestor = T1_a;
					while(T1_a_ancestor != NULL && T1_a_ancestor != B_1_lca) {
						T1_a_ancestor = T1_a_ancestor->parent();
					}
					if (T1_a_ancestor != NULL) {
						#ifdef DEBUG
							cout << "\tmultifurcating cut_b_only so will also cut c" << endl;
						#endif
						cut_b_only=false;
					}
				}

				// cut T2_b
				if ((!CUT_AC_SEPARATE_COMPONENTS || same_component)
//						&& ((!T2_b->parent()->is_protected()
								&& (((multi_node || !T2_b->is_protected())))
						&& (!ABORT_AT_FIRST_SOLUTION || best_k < 0
							|| !PREFER_RHO || !AFs->front().first.contains_rho() )
						&& !cut_a_only && !cut_c_only
						&& (T2_a->parent()->parent() != NULL
								|| !T2_a->is_protected()
								|| (T2_a->parent() == T2->get_component(0)
										&& !T2->contains_rho()))) {
					if (multi_node) {
						um.add_event(new ChangeEdgePreInterval(T2_a));
						T2_a->copy_edge_pre_interval(T2_ab);
						um.add_event(new CutParent(T2_a));
						T2_a->cut_parent();
						um.add_event(new ChangeEdgePreInterval(T2_ab));
						T2_ab->set_edge_pre_start(-1);
						T2_ab->set_edge_pre_end(-1);
						Node *T2_ab_parent = T2_ab->parent();
						if (T2_ab_parent != NULL) {
							um.add_event(new CutParent(T2_ab));
							T2_ab->cut_parent();
							um.add_event(new AddChild(T2_a));
							T2_ab_parent->add_child(T2_a);
							um.add_event(new AddComponent(T2));
							T2->add_component(T2_ab);
						}
						else {
							if (T2->get_component(0) == T2_ab) {
								um.add_event(new AddComponentToFront(T2));
								T2->add_component(0, T2_a);
							}
							else {
								um.add_event(new AddComponent(T2));
								T2->add_component(T2_a);
								singletons->push_back(T2_a);
							}
						}
					}
					else {
						um.add_event(new CutParent(T2_b));
						T2_b->cut_parent();
						ContractEvent(&um, T2_ab);
						node = T2_ab->contract();
						if (node != NULL && node->is_singleton()
								&& node != T2->get_component(0))
								singletons->push_back(node);
						um.add_event(new AddComponent(T2));
						T2->add_component(T2_b);
						if (T2_b->is_leaf())
							singletons->push_back(T2_b);
					}
				add_sibling_pair(sibling_pairs, T1_a, T1_c,
						&um);

					// TODO: check carefully

					if (cut_a_or_merge_ac) {
						if (!T2_a->is_protected()) {
							um.add_event(new ProtectEdge(T2_a));
							T2_a->protect_edge();
							um.add_event(new ListPushBack(protected_stack));
							protected_stack->push_back(T2_a);
						}
						if (!T2_c->is_protected()) {
							um.add_event(new ProtectEdge(T2_c));
							T2_c->protect_edge();
						}
					}

					if (CUT_ALL_B) {
						answer_b =
							rSPR_branch_and_bound_hlpr(T1, T2, k-1,
									sibling_pairs, singletons, true, AFs, protected_stack,
									num_ties, T1_a, T1_c);
					}
					else {
						answer_b =
							rSPR_branch_and_bound_hlpr(T1, T2, k-1,
									sibling_pairs, singletons, false, AFs, protected_stack,
									num_ties, T1_a, T1_c);
					}
				}

				if (answer_b > best_k
						|| (answer_b == best_k
							&& PREFER_RHO
							&& T2->contains_rho() )) {
					best_k = answer_b;
					//swap(&best_T1, &T1);
					//swap(&best_T2, &T2);
				}

				um.undo_to(undo_state);

				/*
				delete T1;
				delete T2;
				delete sibling_pairs;
				delete singletons;
				*/
				rspr_branch_and_bound_cut_c_hlpr(T1, T2, k, sibling_pairs, 
					singletons, AFs, protected_stack, num_ties, T1_a, T2_a, T2_b, T2_c,
					cut_a_only, cut_b_only, cut_c_only, path_length, &um, balanced, 
					multi_b1, multi_b2, cut_a_or_merge_ac, cut_ab_only, T2_d, lca_depth, best_k);

				/*
				delete T1;
				delete T2;
				delete sibling_pairs;
				delete singletons;
				*/

				um.undo_to(undo_state);

				//T1 = best_T1;
				//T2 = best_T2;

#ifdef DEBUG_UNDO
		 while(um.num_events() > 0) {
				cout << "Undo step " << um.num_events() << endl;
				cout << "T1: ";
				T1->print_components();
				cout << "T2: ";
				T2->print_components();
					cout << "sibling pairs:";
					for (set<SiblingPair>::iterator i = sibling_pairs->begin(); i != sibling_pairs->end(); i++) {
						cout << "  ";
						(*i).a->print_subtree_hlpr();
						cout << ",";
						(*i).c->print_subtree_hlpr();
					}
					cout << endl;
			 um.undo();
			cout << endl;
		 }
#else
		 um.undo_all();
#endif
				singletons->clear();
				return best_k;
			} // 3rd case
			cut_b_only = false;
		} // sibling pairs
	} // both singletons and sibling pairs

	if (k >= 0) {
		if (PREFER_RHO && !AFs->empty() && !AFs->front().first.contains_rho() && T1->contains_rho()) {
			if (!ALL_MAFS)
				AFs->clear();
			AFs->push_front(make_pair(Forest(T1),Forest(T2)));
			*num_ties = 2;
		}
		else if (ALL_MAFS || AFs->empty()) {
			AFs->push_back(make_pair(Forest(T1),Forest(T2)));
		}
		else if (!PREFER_RHO || AFs->front().first.contains_rho() == T1->contains_rho()) {
			if (rand() < RAND_MAX/ *num_ties) {
				AFs->clear();
				AFs->push_back(make_pair(Forest(T1),Forest(T2)));
			}
			(*num_ties)++;
		}
	}

#ifdef DEBUG_UNDO
		 while(um.num_events() > 0) {
				cout << "Undo step " << um.num_events() << endl;
				cout << "T1: ";
				T1->print_components();
				cout << "T2: ";
				T2->print_components();
					cout << "sibling pairs:";
					for (set<SiblingPair>::iterator i = sibling_pairs->begin(); i != sibling_pairs->end(); i++) {
						cout << "  ";
						(*i).a->print_subtree_hlpr();
						cout << ",";
						(*i).c->print_subtree_hlpr();
					}
					cout << endl;
			 um.undo();
			cout << endl;
		 }
#else
		 um.undo_all();
#endif

	return k;
}

void rspr_branch_and_bound_cut_a_hlpr(Forest *T1, Forest *T2, int k,
	set<SiblingPair> *sibling_pairs, list<Node *> *singletons, list<pair<Forest,Forest> > *AFs,
	list<Node *> *protected_stack, int *num_ties, Node *T1_c, Node *T2_a, Node *T2_b, Node* T2_c,
	bool cut_a_only, bool cut_b_only, bool cut_c_only, int path_length, UndoMachine *um,
	Node *T2_ab, bool balanced, bool multi_b1, bool multi_b2, Node *T2_d, int &best_k)
{
	int answer_a = -1;
	Node *node;

	// cut T2_a
	//				if (cut_b_only == false && T2_a->is_protected())
	//					cout << "protected k=" << k << endl;
	if (cut_b_only == false 
		&& cut_c_only == false &&
			!T2_a->is_protected()
			&& (T2_a->parent()->parent() != NULL
					|| (T2_a->parent() == T2->get_component(0)
							&& !T2->contains_rho())
					|| !T2_b->is_protected()
					|| T2_a->parent()->get_children().size() > 2)) {// &&
					//(!T2_a->parent()->is_protected() ||
					//T2_a->parent()->get_children().size() > 2)) { }
		um->add_event(new CutParent(T2_a));
		T2_a->cut_parent();
		ContractEvent(um, T2_ab);
		node = T2_ab->contract();
		if (node != NULL && node->is_singleton() &&
				node != T2->get_component(0))
			singletons->push_back(node);
		um->add_event(new AddComponent(T2));
		T2->add_component(T2_a);
		singletons->push_back(T2_a);

		// also require !cut_a_only ?
		//					if (EDGE_PROTECTION_TWO_B && T2_c->is_protected() && !cut_a_only) {
		//				}

		if (EDGE_PROTECTION_TWO_B && T2_c->is_protected() && !cut_a_only){
			if (path_length == 4) {
				if (!multi_b1 && !multi_b2 && !T2_b->is_protected()) {
					um->add_event(new ProtectEdge(T2_b));
					T2_b->protect_edge();
				}
				if (!multi_b2 && !multi_b1) {
					Node *T2_b2 = T2_b->parent()->get_sibling();
					if (balanced)
						T2_b2 = T2_d;
					if (!T2_b2->is_protected()) {
						um->add_event(new ProtectEdge(T2_b2));
						T2_b2->protect_edge();
					}
				}
			}
		}

		if (cut_a_only) {
			answer_a =
				rSPR_branch_and_bound_hlpr(T1, T2, k-1, sibling_pairs,
						singletons, false, AFs, protected_stack, num_ties, T1_c, T1_c->get_sibling());
		}
		else {
			answer_a =
				rSPR_branch_and_bound_hlpr(T1, T2, k-1, sibling_pairs,
						singletons, false, AFs, protected_stack, num_ties);
		}
	}

	if (answer_a > best_k
			|| (answer_a == best_k
				&& PREFER_RHO
				&& T2->contains_rho() )) {
		best_k = answer_a;
	}
}

void rspr_branch_and_bound_cut_c_hlpr(Forest *T1, Forest *T2, int k,
	set<SiblingPair> *sibling_pairs, list<Node *> *singletons, list<pair<Forest,Forest> > *AFs,
	list<Node *> *protected_stack, int *num_ties, Node *T1_a, Node *T2_a, Node *T2_b, Node* T2_c,
	bool cut_a_only, bool cut_b_only, bool cut_c_only, int path_length, UndoMachine *um,
	bool balanced, bool multi_b1, bool multi_b2, bool cut_a_or_merge_ac, bool cut_ab_only, 
	Node *T2_d, int lca_depth, int &best_k) 
{	
	int answer_c = -1;
	Node *node;
	
	// load the copy
	/*
	T1 = T1_copy;
	T2 = T2_copy;
	T1_a = T1_a_copy;
	T1_c = T1_c_copy;
	T2_a = T2_a_copy;
	T2_c = T2_c_copy;
	sibling_pairs = sibling_pairs_copy;
	singletons = new list<Node *>();
	*/
	//				if (T2_c->is_protected())
	//					cout << "protected k=" << k << endl;
	if (!T2_c->is_protected() &&
			!cut_a_or_merge_ac &&
	//					(T2_c->parent() == NULL || !T2_c->parent()->is_protected() ||
	//						T2_c->parent()->get_children().size() > 2) &&
			(!ABORT_AT_FIRST_SOLUTION || best_k < 0
				|| !PREFER_RHO || !AFs->front().first.contains_rho() )
			&& cut_b_only == false && cut_ab_only == false
			&& cut_a_only == false
			// TODO: do we allow this if T2_c has no parent?
			// it has to be under rho, right?
			&& (T2_c->parent() == NULL
					|| T2_c->parent()->parent() != NULL
					|| (T2_c->parent() == T2->get_component(0)
							&& !T2->contains_rho())
					|| !T2_c->get_sibling()->is_protected()
					|| T2_c->parent()->get_children().size() > 2)) {// &&


		if (T2_c->parent() != NULL) {
			Node *T2_c_parent = T2_c->parent();
			um->add_event(new CutParent(T2_c));
			T2_c->cut_parent();
			ContractEvent(um, T2_c_parent);
			node = T2_c_parent->contract();
			if (node != NULL && node->is_singleton()
					&& node != T2->get_component(0))
				singletons->push_back(node);
			um->add_event(new AddComponent(T2));
			T2->add_component(T2_c);
		}
		else {
			// don't decrease k
			k++;
		}

		if (EDGE_PROTECTION && !cut_c_only) {
			if (!T2_a->is_protected()) {
				um->add_event(new ProtectEdge(T2_a));
				T2_a->protect_edge();
//							if (DEEPEST_PROTECTED_ORDER && !cut_c_only) {
				if (DEEPEST_PROTECTED_ORDER) {
					um->add_event(new ListPushBack(protected_stack));
					protected_stack->push_back(T2_a);
				}
				// TODO: add to protected list
			}
//						if (EDGE_PROTECTION_TWO_B && !cut_c_only) {
//					}
			// TODO: problem here :(
			if (EDGE_PROTECTION_TWO_B) {
				if (path_length == 4) {
					if (!multi_b1 && !multi_b2 && !T2_b->is_protected()) {
						um->add_event(new ProtectEdge(T2_b));
						T2_b->protect_edge();
					}
					if (!multi_b2 && !multi_b1) {
						Node *T2_b2 = T2_b->parent()->get_sibling();
						if (balanced)
							T2_b2 = T2_d;
						if (!T2_b2->is_protected()) {
							um->add_event(new ProtectEdge(T2_b2));
							T2_b2->protect_edge();
						}
					}
				}
			}
			if (path_length == 5) 
				T2_a->set_max_merge_depth(lca_depth);
		}

		singletons->push_back(T2_c);
		if (cut_c_only) {
			answer_c =
				rSPR_branch_and_bound_hlpr(T1, T2, k-1, sibling_pairs,
						singletons, false, AFs, protected_stack, num_ties, T1_a, T1_a->get_sibling());
		}
		else {
			answer_c =
				rSPR_branch_and_bound_hlpr(T1, T2, k-1, sibling_pairs,
						singletons, false, AFs, protected_stack, num_ties);
		}
		if (answer_c > best_k
					|| (answer_c == best_k
					&& PREFER_RHO
					&& T2->contains_rho() )) {
			best_k = answer_c;
			//swap(&best_T1, &T1);
			//swap(&best_T2, &T2);
		}
	}
}

int rSPR_branch_and_bound_simple_clustering(Node *T1, Node *T2, bool verbose, map<string, int> *label_map, map<int, string> *reverse_label_map) {
	return rSPR_branch_and_bound_simple_clustering(T1,T2, verbose, label_map, reverse_label_map, -1, -1, NULL, NULL);
}

int rSPR_branch_and_bound_simple_clustering(Node *T1, Node *T2, bool verbose, map<string, int> *label_map, map<int, string> *reverse_label_map, int min_k, int max_k) {
	return rSPR_branch_and_bound_simple_clustering(T1,T2, verbose, label_map, reverse_label_map, min_k, max_k, NULL, NULL);
}

int rSPR_branch_and_bound_simple_clustering(Node *T1, Node *T2, bool verbose, map<string, int> *label_map, map<int, string> *reverse_label_map, int min_k, int max_k, Forest **out_F1, Forest **out_F2) {
	bool do_cluster = true;
	if (max_k > MAX_SPR)
		max_k = MAX_SPR;
	else if (max_k == -1)
		max_k = INT_MAX;

	if (T2->get_preorder_number() == -1) {
	  T2->preorder_number();
	}

	ClusterForest F1 = ClusterForest(T1);
	ClusterForest F2 = ClusterForest(T2);
	Forest F3 = Forest(F1);
	Forest F4 = Forest(F2);


//	bool old_rho = PREFER_RHO;
	PREFER_RHO = true;
	if (verbose) {
		cout << "T1: ";
		F1.print_components();
		cout << "T2: ";
		F2.print_components();
	}
	int full_approx_spr;
	if (MULTIFURCATING) {
	  full_approx_spr = rSPR_worse_3_mult_approx(&F3, &F4);
	}
	else {
	  full_approx_spr = rSPR_worse_3_approx(&F3, &F4);
	}
	/*if (full_approx_spr < CLUSTER_TUNE) {
		do_cluster = false;
	}*/
	if (verbose) {

		cout << "approx F1: ";
		F3.print_components();
		cout << "approx F2: ";
		F4.print_components();
		// what the AF shows
		cout << "approx drSPR=" << F4.num_components()-1 << endl;
		/* what we use to get the lower bound: 3 * the number of cutting rounds in
			 the approx algorithm
		*/
		//cout << "approx drSPR=" << full_approx_spr << endl;
		cout << "\n";
	}
//	if (F1.get_component(0)->get_preorder_number() == -1)
//		F1.get_component(0)->preorder_number();
//	if (F2.get_component(0)->get_preorder_number() == -1)
//		F2.get_component(0)->preorder_number();

	if (!sync_twins(&F1, &F2))
		return 0;
	if (F1.get_component(0)->is_leaf())
		return 0;
	if (F1.get_component(0)->get_preorder_number() == -1) {
		F1.get_component(0)->preorder_number();
		F2.get_component(0)->preorder_number();
	}
	int loss = 0;
	list<Node *> *cluster_points;
	if (F1.get_component(0)->get_edge_pre_start() == -1) {
		F1.get_component(0)->edge_preorder_interval();
		F2.get_component(0)->edge_preorder_interval();
	}
	if (LEAF_REDUCTION2) {
	  
	  if (MULTIFURCATING) {
	    reduction_leaf_mult(&F1, &F2);
	  }
	  else {	    
	    reduction_leaf(&F1, &F2);
	    }
//		F1.get_component(0)->preorder_number();
//		F2.get_component(0)->preorder_number();
//		F1.get_component(0)->edge_preorder_interval();
//		F2.get_component(0)->edge_preorder_interval();
	}
	if (COUNT_LOSSES) {
		loss += F1.get_component(0)->count_lost_subtree();
		loss += F2.get_component(0)->count_lost_subtree();
	}
		//F1.print_components();
		//F2.print_components();
	
	ClusterMergeForest* cluster_merge_forest = NULL;
	int cur_cluster_index = 0;
	if (do_cluster) {
		sync_interior_twins(&F1, &F2);
		cluster_points = find_cluster_points(&F1, &F2);

		if(ALL_MERGED_MAFS){
			cluster_merge_forest = new ClusterMergeForest(cluster_points->size()+1);
		}

		//	list<Node *> *cluster_points = new list<Node *>();
		for(list<Node *>::iterator i = cluster_points->begin();
				i != cluster_points->end(); i++) {
			string cluster_name = "X";
			/*
			if (verbose) {
					stringstream ss;
					ss << F1.size();
					cluster_name += ss.str();
					//int num_labels = label_map.size();
					//label_map.insert(make_pair(cluster_name,num_labels));
					//reverse_label_map.insert(
					//		make_pair(num_labels,cluster_name));
					//ss.str("");
					//ss << num_labels;
					//cluster_name = ss.str();
			}
			*/
	
			Node *n = *i;
			/*if (n->parent()->parent() == NULL
					&& n->get_sibling() != NULL &&
					n->get_sibling()->get_name() == "X")
				continue;*/
			Node *n_parent = n->parent();
			Node *twin = n->get_twin();
			Node *twin_parent = twin->parent();
	
			if (twin_parent == NULL)
				continue;

			if(ALL_MERGED_MAFS && cluster_merge_forest != NULL){
				cur_cluster_index++;
				cluster_merge_forest->assign_cluster_parent(n, twin, cur_cluster_index);
			}

			/*cout << "Cluster prenums " << endl;
			for(int j=0; j<cluster_prenums.size(); j++){
				cout << cluster_prenums[j] << endl;
			}
			cout << "Cluster Parents " << endl;
			for(int j=0; j<cluster_parent_index.size(); j++){
				cout << cluster_parent_index[j] << endl;
			}*/
	
			F1.add_cluster(n,cluster_name);
	
			F2.add_cluster(twin,cluster_name);
	
			Node *n_cluster =
					F1.get_cluster_node(F1.num_clusters()-1);
			Node *twin_cluster =
					F2.get_cluster_node(F2.num_clusters()-1);
			n_cluster->set_twin(twin_cluster);
			twin_cluster->set_twin(n_cluster);
	
		}
		if (verbose)
			cout << endl << "CLUSTERS" << endl;

	}

	// component 0 needs to be done last
	F1.add_component(F1.get_component(0));
	F2.add_component(F2.get_component(0));

	int k;
	int num_clusters = F1.num_components();
	int total_k = 0;

	for(int i = 1; i < num_clusters; i++) {
		if (i == num_clusters - 1) {
			PREFER_RHO = false;
		}
		int exact_spr = -1;
		//vector<Node *> comps = vector<Node *>();
		//comps.push_back(F1.get_component(i));
		Forest f1 = Forest(F1.get_component(i));
		//Forest f1 = Forest(comps);
		//comps.clear();

		//comps.push_back(F1.get_component(i));
		Forest f2 = Forest(F2.get_component(i));
		//Forest f2 = Forest(comps);
		//comps.clear();
		Forest f1a = Forest(f1);
		Forest f2a = Forest(f2);
		Forest *f1_cluster;
		Forest *f2_cluster;

		if (verbose) {
			cout << "C" << i << "_1: ";
			f1.print_components();
			cout << "C" << i << "_2: ";
			f2.print_components();
		}
		int approx_spr;
		if (MULTIFURCATING) {
		  approx_spr = rSPR_worse_3_mult_approx(&f1a, &f2a);
		}
		else {
		  approx_spr = rSPR_worse_3_approx(&f1a, &f2a);
		}
		if (verbose) {
			cout << "cluster approx drSPR=" << f2a.num_components()-1 << endl;
			//cout << "cluster approx drSPR=" << approx_spr << endl;

			cout << endl;
		}

		int min_spr = approx_spr / 3;
		if (min_spr < MIN_SPR - total_k)
			min_spr = MIN_SPR - total_k;
		int total_split_k = 0;

		bool done_cluster = false;
		bool done_split = false;

		double tree_fraction = INITIAL_TREE_FRACTION;

		if (min_spr < min_k)
			min_spr = min_k;

		while(!done_cluster) {
			done_cluster = true;

			for(k = min_spr - total_split_k; true; k++) {
				if (k < 0)
					k = 0;
				if (SPLIT_APPROX && !done_split && k >= SPLIT_APPROX_THRESHOLD) {
					done_cluster = false;
					break;
				}
				Forest f1t = Forest(f1);
//				Forest f1t = f1;
				Forest f2t = Forest(f2);
//				Forest f2t = f2;
				
				#ifdef DEBUG_LGT_EVENTS
					f2t.print_protected_edge_list();
				#endif

				f1t.unsync();
				f2t.unsync();
				exact_spr = -1;
				if (verbose) {
					cout << k << " ";
  				cout.flush();
				}

				list<pair<Forest,Forest>> extAFs = list<pair<Forest,Forest>>();
				if (k + total_k <= max_k && k <= CLUSTER_MAX_SPR) {
					if (f1t.get_component(0)->get_name() == DEAD_COMPONENT) {
						f1t.add_rho();
						f2t.add_rho();
					}

					if (MULTIFURCATING) {
					  exact_spr = rSPR_branch_and_bound_mult(&f1t, &f2t, k);
					}
					else {
						if(ALL_MERGED_MAFS)
							exact_spr = rSPR_branch_and_bound(&f1t, &f2t, k, &extAFs);
						else
							exact_spr = rSPR_branch_and_bound(&f1t, &f2t, k);
					}
				}
				if (exact_spr >= 0 || k + total_k > max_k ||
						k > CLUSTER_MAX_SPR) {
					if (k > CLUSTER_MAX_SPR) {
						f1t.swap(&f1a);
						f2t.swap(&f2a);
						cout << "foo" << endl;
					}
					if (exact_spr >= 0) {
						exact_spr += total_split_k;
						if (verbose) {
	  					cout << endl;
	  					cout << "F" << i << "_1: ";
	  					f1t.print_components();
	  					cout << "F" << i << "_2: ";
	  					f2t.print_components();
	  					cout << "cluster exact drSPR=" << exact_spr << endl;
	  					cout << endl;
						}
						total_k += exact_spr;
					}
					else {
						// TODO: don't just the MAX_SPR here
						// incorporate extra information
						// toggle?
						if (verbose) {
							cout << "cluster exact drSPR=?  " << "k=" << k << " too large"
								<< endl;
							cout << "\n";
						}
						if (false && k > CLUSTER_MAX_SPR) {
							// TODO: this should be an approx of the remaining forest
//							total_k += approx_spr;
						}
						else if (CLAMP) {
							total_k = max_k;
						}
						else {
							Forest f1a = Forest(f1);
							Forest f2a = Forest(f2);
							
							int approx_spr;
							if (MULTIFURCATING) {
							  approx_spr = rSPR_worse_3_mult_approx(&f1a, &f2a);
							}
							else{
							  approx_spr = rSPR_worse_3_approx(&f1a, &f2a);
							}
								//total_k += min_spr;
								total_k += approx_spr / 3;
						}
					}
					if(ALL_MERGED_MAFS && cluster_merge_forest) {
						bool has_cluster_node = true;
						if (f1t.contains_rho())
							has_cluster_node = false;

						vector<int> rho_indexes = cluster_merge_forest->get_forests_with_rho(extAFs);
						if(rho_indexes.size() == 0 && i < num_clusters - 1){
							//Run with one more distance
							Forest f1t1 = Forest(f1);
							Forest f2t1 = Forest(f2);
							f1t1.unsync();
							f2t1.unsync();
							f1t1.add_rho();
							f2t1.add_rho();

							list<pair<Forest,Forest>> extRhoAFs = list<pair<Forest,Forest>>();
							int new_k = rSPR_branch_and_bound(&f1t1, &f2t1, k, &extRhoAFs);
							cluster_merge_forest->add_non_duplicate_maf(&extAFs, &extRhoAFs);
						}
						cluster_merge_forest->update_merged_forests(i, num_clusters, &extAFs, has_cluster_node, exact_spr);
					}
					if ( i < num_clusters - 1) {
						F1.join_cluster(i, &f1t);
						F2.join_cluster(i, &f2t);
					}
					else {
						F1.join_cluster(&f1t);
						F2.join_cluster(&f2t);
						if(ALL_MERGED_MAFS && cluster_merge_forest) {
							cluster_merge_forest->print_merged_agreement_forest(i, label_map, reverse_label_map);
						}
					}
					break;
				}
			}
			done_split = done_cluster;
			int num_splits = 0;
			while (SPLIT_APPROX && !done_split) {
				//IN_SPLIT_APPROX = true;
				Node *original_split_node = find_subtree_of_approx_distance(
						f1.get_component(0), &f1, &f2, SPLIT_APPROX_THRESHOLD*2);
				if (original_split_node == f1.get_component(0) &&
						num_splits > 0)
					done_split = true;
				else {
					Forest f1a = Forest(f1);
					Forest f2a = Forest(f2);
					Node *a_split_node =
					f1a.find_by_prenum(original_split_node->get_preorder_number());
					f1a.get_component(0)->disallow_siblings_subtree();
						a_split_node->allow_siblings_subtree();
//					if (a_split_node->lchild() != NULL)
//						a_split_node->lchild()->allow_siblings_subtree();
//					if (a_split_node->rchild() != NULL)
//						a_split_node->rchild()->allow_siblings_subtree();
					// something odd going on here
					int start = rSPR_worse_3_approx(a_split_node, &f1a, &f2a);
					if (start == INT_MAX)
						start = 0;
					start /= 3;
					int end = f1.get_component(0)->size();
					for(k = start; true; k++) {
						// TODO: figure out the bug here
						if (k > end) {
							k = 0;
							done_split = true;
							break;
						}
				/*	if (k > SPLIT_APPROX_THRESHOLD) {
						k = 0;
						tree_fraction *= 0.75;
						if (verbose)
							cout << "tree_fraction: " << tree_fraction << endl;
						continue;
					}*/
						Forest f1s = Forest(f1);
						Forest f2s = Forest(f2);
						if (!sync_twins(&f1s, &f2s)) {
							k = 0;
							done_split = true;
							break;
						}
						if (verbose) {
							cout << k << " ";
		  				cout.flush();
						}
						Node *split_node = f1s.find_by_prenum(original_split_node->get_preorder_number());
						f1s.get_component(0)->disallow_siblings_subtree();
							split_node->allow_siblings_subtree();
//						if (split_node->lchild() != NULL)
//							split_node->lchild()->allow_siblings_subtree();
//						if (split_node->rchild() != NULL)
//							split_node->rchild()->allow_siblings_subtree();
							//f1s.get_component(0)->find_subtree_of_size(tree_fraction);
							set<SiblingPair > *sibling_pairs =
								find_sibling_pairs_set(split_node);
							list<Node *> singletons = f2s.find_singletons();
							list<pair<Forest,Forest> > AFs = list<pair<Forest,Forest> >();
							list<Node *> protected_stack = list<Node *>();

							int num_ties = 2;

							int split_k = rSPR_branch_and_bound_hlpr(&f1s, &f2s, k,
									sibling_pairs, &singletons, false, &AFs,
									&protected_stack, &num_ties);
							delete sibling_pairs;
							if (!AFs.empty()) {
								AFs.front().first.swap(&f1);
								AFs.front().second.swap(&f2);
								f2.unprotect_edges();
								f1.get_component(0)->allow_siblings_subtree();
								AFs.clear();
								total_split_k += k - split_k;
		//						if (k < SPLIT_APPROX_THRESHOLD * 0.75) {
		//							tree_fraction *= 2;
		//							if (tree_fraction > INITIAL_TREE_FRACTION)
		//								tree_fraction = INITIAL_TREE_FRACTION;
		//						}
								if (verbose)
									cout << "split_k: " << k << endl;
								break;
							}
					}
				}
				//IN_SPLIT_APPROX = false;
				num_splits++;
			}
	
			// TODO: approx again? seperate approxes ?
		}
	}

		if (F1.contains_rho()) {
			F1.get_component(0)->delete_tree();
			F2.get_component(0)->delete_tree();
			F1.erase_components(0, num_clusters);
			F2.erase_components(0, num_clusters);
		}
		else {
			F1.get_component(num_clusters)->delete_tree();
			F2.get_component(num_clusters)->delete_tree();
			F1.erase_components(1, num_clusters+1);
			F2.erase_components(1, num_clusters+1);
		}
		// fix hanging roots
		for(int i = 0; i < F1.num_components(); i++) {
			F1.get_component(i)->contract(true);
			F2.get_component(i)->contract(true);
		}

		if (verbose) {
			F1.numbers_to_labels(reverse_label_map);
			F2.numbers_to_labels(reverse_label_map);
			cout << "F1: ";
			F1.print_components();
			cout << "F2: ";
			F2.print_components();
			cout << "total exact drSPR=" << total_k << endl;
		}
		if (out_F1 != NULL)
			*out_F1 = new Forest(&F1);
			//F1.swap(out_F1);
		if (out_F2 != NULL)
			*out_F2 = new Forest(&F2);
//			F2.swap(out_F2);
		if (out_F1 != NULL && out_F2 != NULL) {
//			out_F1->resync();
			sync_twins(*out_F1, *out_F2);
	//		sync_interior_twins_real(out_F1, out_F2);
		}

	if (do_cluster) {
		delete cluster_points;
	}
//	PREFER_RHO = old_rho;
	total_k += loss;
/*	cout << "F1: ";
	for (int i = 0; i < F1.num_components(); i++) {
		if (i > 0)
			cout << " ";
		F1.get_component(i)->expand_contracted_nodes();
		cout << F1.get_component(i)->str_edge_pre_interval_subtree();
	}
	cout << endl;
	cout << "F2: ";
	for (int i = 0; i < F2.num_components(); i++) {
		if (i > 0)
			cout << " ";
		F2.get_component(i)->expand_contracted_nodes();
		cout << F2.get_component(i)->str_edge_pre_interval_subtree();
	}
	cout << endl << endl;
*/
	#ifdef DEBUG_CASE_COUNTER
	print_mult_case_count();
        #endif

	return total_k;
}

int rSPR_branch_and_bound_simple_clustering(Forest *T1, Forest *T2, bool verbose, map<string, int> *label_map, map<int, string> *reverse_label_map) {
	Forest F1 = *T1;//Forest(T1);
	Forest F2 = *T2;//Forest(T2);
	Forest F3 = Forest(F1);
	Forest F4 = Forest(F2);

	bool do_cluster = true;

//	bool old_rho = PREFER_RHO;
	PREFER_RHO = true;
	if (verbose) {
		cout << "T1: ";
		F1.print_components();
		cout << "T2: ";
		F2.print_components();
	}

	int full_approx_spr;
	if (MULTIFURCATING) {
	  full_approx_spr = rSPR_worse_3_mult_approx(&F3, &F4);
	}
	else {
	  full_approx_spr = rSPR_worse_3_approx(&F3, &F4);
	}
	if (full_approx_spr <= CLUSTER_TUNE) {
		do_cluster = false;
	}
	if (verbose) {

		cout << "approx F1: ";
		F3.print_components();
		cout << "approx F2: ";
		F4.print_components();
		// what the AF shows
		cout << "approx drSPR=" << F4.num_components()-1 << endl;
		/* what we use to get the lower bound: 3 * the number of cutting rounds in
			 the approx algorithm
		*/
		//cout << "approx drSPR=" << full_approx_spr << endl;
		cout << "\n";
	}

	if (!sync_twins(&F1, &F2))
		return 0;
	if (F1.get_component(0)->is_leaf())
		return 0;
	list<Node *> *cluster_points;
	sync_interior_twins_real(&F1, &F2);
	if (do_cluster) {
		cluster_points = find_cluster_points(&F1, &F2);
	}

		int total_k = 0;

		if (do_cluster && !cluster_points->empty()) {
	
			list<ClusterInstance> clusters =
				cluster_reduction(&F1, &F2, cluster_points);
	
			F1.unsync_interior();
			F2.unsync_interior();
		int k = 0;
		int i = 0;
		while(!clusters.empty()) {
			i++;
			ClusterInstance cluster = clusters.front();
			clusters.pop_front();
			cluster.F1->unsync_interior();
			cluster.F2->unsync_interior();
			if (verbose) {
				cout << "C" << i << "_1: ";
				cluster.F1->print_components();
				cout << "C" << i << "_2: ";
				cluster.F2->print_components();
			}
			Forest f1 = Forest(cluster.F1);
			Forest f2 = Forest(cluster.F2);
			
			int min_spr;
			if (MULTIFURCATING) {
			  min_spr = rSPR_worse_3_mult_approx(&f1, &f2);
			}
			else {
			  min_spr = rSPR_worse_3_approx(&f1, &f2);
			}
			min_spr /= 3;

			if (verbose) {
				cout << "cluster approx drSPR=" << f1.num_components()-1 << endl;
				//cout << "cluster approx drSPR=" << min_spr << endl;
				cout << endl;
			}

			int cluster_spr = -1;
			k = MAX_SPR - total_k;
			if (k >= 0) {
				// hack for clusters with no rho
				if ((cluster.F2_cluster_node == NULL
							|| (cluster.F2_cluster_node->is_leaf()
									&& cluster.F2_cluster_node->parent() == NULL
									&& cluster.F2_cluster_node->
									get_num_clustered_children() <= 1
									&& (cluster.F2_cluster_node !=
											cluster.F2_cluster_node->get_forest()->
												get_component(0))))
							&& cluster.F2_has_component_zero == false) {
					cluster.F1->add_rho();
					cluster.F2->add_rho();
				}
				if (MULTIFURCATING) {
				  cluster_spr = rSPR_branch_and_bound_mult_range(cluster.F1,
						cluster.F2, min_spr, MAX_SPR - total_k);
				}
				else {
				  cluster_spr = rSPR_branch_and_bound_range(cluster.F1,
						cluster.F2, min_spr, MAX_SPR - total_k);
				}
				if (cluster_spr >= 0) {
					if (verbose) {
	  				cout << endl;
	  				cout << "F" << i << "_1: ";
	  				cluster.F1->print_components();
	  				cout << "F" << i << "_2: ";
	  				cluster.F2->print_components();
	  				cout << "cluster exact drSPR=" << cluster_spr << endl;
	  				cout << endl;
					}
					total_k += cluster_spr;
					k -= cluster_spr;
				}
				else {
					k = -1;
					if (verbose) {
						cout << "cluster exact drSPR=?  " << "k=" << k << " too large"
							<< endl;
						cout << "\n";
					}
					if (CLAMP) {
						total_k = MAX_SPR + 1;
					}
					else {
							total_k += min_spr;
					}
				}
			}
			if (k > -1) {
				if (!cluster.is_original()) {
					int adjustment = cluster.join_cluster(&F1, &F2);
					total_k += adjustment;
					delete cluster.F1;
					delete cluster.F2;
				}
			}
			else {
				if (!cluster.is_original()) {
					//if (cluster.F1_cluster_node != NULL)
					//	cluster.F1_cluster_node->contract();
					//if (cluster.F2_cluster_node != NULL)
					//	cluster.F2_cluster_node->contract();
					delete cluster.F1;
					delete cluster.F2;
				}
			}
		}
		delete cluster_points;
	}
	else {
		if (do_cluster) {
			delete cluster_points;
		}
		full_approx_spr /= 3;
		if (MULTIFURCATING) {
		  total_k = rSPR_branch_and_bound_mult_range(&F1, &F2, full_approx_spr, MAX_SPR);
		}
		else {
		  total_k = rSPR_branch_and_bound_range(&F1, &F2, full_approx_spr, MAX_SPR);
		}
		int i = 1;
		if (total_k < 0) {
			if (CLAMP) {
				total_k = MAX_SPR;
			}
			else {
		 	 total_k = full_approx_spr;
			}
		}

	}

	if (verbose) {
		F1.numbers_to_labels(reverse_label_map);
		F2.numbers_to_labels(reverse_label_map);
	 	cout << endl;
	 	cout << "F1: ";
	 	F1.print_components();
	 	cout << "F2: ";
	 	F2.print_components();
	 	cout << "total exact drSPR=" << total_k << endl;
	 	cout << endl;
	}
	return total_k;
}

/*Joel's part*/
int rSPR_branch_and_bound_simple_clustering(Forest *T1, Forest *T2, bool verbose){
	return rSPR_branch_and_bound_simple_clustering(T1, T2, false, NULL, NULL);
}

int rSPR_branch_and_bound_simple_clustering(Forest *T1, Forest *T2){
	return rSPR_branch_and_bound_simple_clustering(T1, T2, false, NULL, NULL);
}

int rSPR_branch_and_bound_simple_clustering(Node *T1, Node *T2, bool verbose) {
	return rSPR_branch_and_bound_simple_clustering(T1, T2, false, NULL, NULL, -1, -1);
}

int rSPR_branch_and_bound_simple_clustering(Node *T1, Node *T2, bool verbose, int min_k, int max_k) {
	return rSPR_branch_and_bound_simple_clustering(T1, T2, false, NULL, NULL, min_k, max_k);
}

int rSPR_branch_and_bound_simple_clustering(Node *T1, Node *T2) {
	return rSPR_branch_and_bound_simple_clustering(T1, T2, false);
}
int rSPR_branch_and_bound_simple_clustering(Node *T1, Node *T2, Forest **out_F1, Forest **out_F2) {
	return rSPR_branch_and_bound_simple_clustering(T1,T2, false, NULL, NULL, -1, -1, out_F1, out_F2);
}

// T1 and T2 are assumed to already be synced
void reduction_leaf_mult(Forest *T1, Forest* T2) {
  
  list<Node *> *sibling_groups = T1->find_sibling_groups();
  while (!sibling_groups->empty()) {
    //Get a sibling group with identical siblings
    list<Node*>::reverse_iterator i = sibling_groups->rbegin();
    Node *T1_sibling_group = sibling_groups->back();
    list<list<Node*>> identical_sibling_groups;
    T1_sibling_group->find_identical_sibling_groups(&identical_sibling_groups);
    for (; i != sibling_groups->rend(); i++ ){
      (*i)->find_identical_sibling_groups(&identical_sibling_groups);
      if (identical_sibling_groups.size() > 0) {
	T1_sibling_group = (*i);
	break;
      }
    }

    // Checked all sibling groups, found none with identical groups in T2
    // Therefore there are no more contractions to be made

    if (identical_sibling_groups.size() == 0) {
      delete sibling_groups;
      return;
      }
    // Contract them
    else {
      list<list<Node *>>::iterator i;
      for (i = identical_sibling_groups.begin(); i != identical_sibling_groups.end(); i++) {
	list<Node *> T2_group = (*i);
	Node *T2_p = T2_group.front()->parent();
	Node *T1_group_new = T1_sibling_group->contract_twin_group(&T2_group);
	Node *T2_group_new = T2_p->contract_sibling_group(&T2_group);
	  
	T1_group_new->set_twin(T2_group_new);
	T2_group_new->set_twin(T1_group_new);			

	if (T1_sibling_group->parent() != NULL) {
	  //Check if the contraction made a new sibling group
	  T1_sibling_group->parent()->recalculate_non_leaf_children();
	  if (T1_sibling_group->parent()->is_sibling_group()) {
	    sibling_groups->push_front(T1_sibling_group->parent());
	  }
	}
	//Check if this contraction removed a sibling group
	if (!T1_sibling_group->is_sibling_group()) {
	  sibling_groups->remove(T1_sibling_group);
	}	  
      }
    }
  }
  delete sibling_groups;
}
  
// T1 and T2 are assumed to already be synced
void reduction_leaf(Forest *T1, Forest *T2) {
	reduction_leaf(T1, T2, NULL);
}

void reduction_leaf(Forest *T1, Forest *T2, UndoMachine *um) {
	list<Node *> *sibling_pairs = T1->find_sibling_pairs();
	Node *T1_a;
	Node *T1_c;
	while (!sibling_pairs->empty()) {
		T1_a = sibling_pairs->front();
		sibling_pairs->pop_front();
		T1_c = sibling_pairs->front();
		sibling_pairs->pop_front();
		// shouldn't happen here
		if (T1_a->parent() == NULL || T1_a->parent() != T1_c->parent()) {
				continue;
		}
		Node *T2_a = T1_a->get_twin();
		Node *T2_c = T1_c->get_twin();
		if (T2_a->parent() != NULL && T2_a->parent() == T2_c->parent()) {
			Node *T1_ac = T1_a->parent();
			Node *T2_ac = T2_a->parent();
			T1_ac->contract_sibling_pair_undoable();
			Node *T2_ac_new = T2_ac->contract_sibling_pair_undoable(T2_a, T2_c);
			if (T2_ac_new != NULL && T2_ac_new != T2_ac) {
				T2_ac = T2_ac_new;
				T2_ac->contract_sibling_pair_undoable();
			}

			T1_ac->set_twin(T2_ac);
			T2_ac->set_twin(T1_ac);

			// check if T2_ac is a singleton
			// also shouldn't happen
//				if (T2_ac->is_singleton() && !T1_ac->is_singleton() && T2_ac != T2->get_component(0))

			// check if T1_ac is part of a sibling pair
			if (T1_ac->parent() != NULL &&
					T1_ac->parent()->is_sibling_pair()) {
				sibling_pairs->push_back(T1_ac->parent()->lchild());
				sibling_pairs->push_back(T1_ac->parent()->rchild());
			}
			#ifdef DEBUG
				cout << "\tT1: ";
				T1->print_components();
				cout << "\tT2: ";
				T2->print_components();
			#endif
		}
	}
	delete sibling_pairs;
}

/* return true if T1_node matches the chain between T2_node and
	 T2_node_end
*/
bool chain_match(Node *T1_node, Node *T2_node, Node *T2_node_end) {
	Node *T1_pendant;
	Node *T2_pendant;
	bool pendant_found = false;
	if (T2_node->is_leaf())
		return false;
	// T1_node is a leaf
	if (T1_node->is_leaf()) {
		T1_pendant = T1_node;
		if (T1_pendant->get_twin() == T2_node->lchild()) {
			if (T2_node->rchild() == T2_node_end)
				return true;
		}
		else if (T1_pendant->get_twin() == T2_node->rchild()) {
			if (T2_node->lchild() == T2_node_end)
				return true;
		}
		return false;
	}
	// T1_pendant is T1_node->lchild()
	T1_pendant = T1_node->lchild();
	if (T1_pendant->is_leaf()) {
		T2_pendant = T2_node->lchild();
		if (T2_pendant->is_leaf() && T1_pendant->get_twin() == T2_pendant) {
			return chain_match(T1_pendant->get_sibling(),
					T2_pendant->get_sibling(), T2_node_end);
		}
		T2_pendant = T2_node->rchild();
		if (T2_pendant->is_leaf() && T1_pendant->get_twin() == T2_pendant) {
			return chain_match(T1_pendant->get_sibling(),
					T2_pendant->get_sibling(), T2_node_end);
		}
	}
	// T1_pendant is T1_node->rchild()
	if (T1_pendant->is_leaf()) {
		T2_pendant = T2_node->lchild();
		if (T2_pendant->is_leaf() && T1_pendant->get_twin() == T2_pendant) {
			return chain_match(T1_pendant->get_sibling(),
					T2_pendant->get_sibling(), T2_node_end);
		}
		T2_pendant = T2_node->rchild();
		if (T2_pendant->is_leaf() && T1_pendant->get_twin() == T2_pendant) {
			return chain_match(T1_pendant->get_sibling(),
					T2_pendant->get_sibling(), T2_node_end);
		}
	}
	return false;
}

int rSPR_total_distance(Node *T1, vector<Node *> &gene_trees) {
	return rSPR_total_distance(T1, gene_trees, NULL);
}

int rSPR_total_distance(Node *T1, vector<Node *> &gene_trees,
		vector<int> *original_scores) {
	int total = 0;
	MAIN_CALL = false;
	int end = gene_trees.size();
//	T1->preorder_number();
	#pragma omp parallel for reduction(+ : total) firstprivate(PREFER_RHO)  // firstprivate(IN_SPLIT_APPROX)
//	for(int j = 0; j < 10; j++)
//	cout << "T1: " << T1->str_subtree() << endl;
	for(int i = 0; i < end; i++) {
			//		cout << i << endl;
	  cout << "Trying tree #" << i << " : " << gene_trees[i]->str_subtree() << endl;
	  int k = rSPR_branch_and_bound_simple_clustering(T1, gene_trees[i], VERBOSE);

	  MULTIFURCATING = true;
		//MULT_4_BRANCH = true;
		int mult_k = rSPR_branch_and_bound_simple_clustering(T1, gene_trees[i], VERBOSE);
		MULTIFURCATING = false;
		//MULT_4_BRANCH = false;
		
		if (k != mult_k) {
		  cout << "BINARY DOES NOT MATCH MULT" << endl;
		  cout << "T1: " << T1->str_subtree() << endl;;
		  cout << "BINARY k = " << k << " mult_k = " << mult_k << endl;
		  break;
		  }
		else {
		  cout << "\tMATCHES: k = " << k << endl;
		}
		
		//cout << "\t: k = " << k << endl;
//		k *= mylog2(gene_trees[i]->size());

		if (original_scores != NULL)
			(*original_scores)[i] = k;
		total += k;
//		cout << "T2: " << gene_trees[i]->str_subtree() << endl;
//		cout << " k: " << k << endl;
		if (FIND_RATE) {
			if (k > 0) {
				int size = gene_trees[i]->find_leaves().size();
//				cout << k << endl;
//				cout << size << endl;
				cout << "rate=" << (float)k / size << endl;
//				cout << T1->str_subtree() << endl;
//				cout << gene_trees[i]->str_subtree() << endl;
			}
		}
//		Forest F1 = Forest(T1);
//		Forest F2 = Forest(gene_trees[i]);
//		total += rSPR_branch_and_bound(&F1, &F2);
	}
	return total;
}

void rSPR_pairwise_distance(Node *T1, vector<Node *> &gene_trees) {
	rSPR_pairwise_distance(T1, gene_trees, 0, gene_trees.size());
}

void rSPR_pairwise_distance(Node *T1, vector<Node *> &gene_trees, bool APPROX) {
	rSPR_pairwise_distance(T1, gene_trees, 0, gene_trees.size(), APPROX);
}

void rSPR_pairwise_distance(Node *T1, vector<Node *> &gene_trees, int start, int end) {
	rSPR_pairwise_distance(T1, gene_trees, start, end, false);
}

void rSPR_pairwise_distance(Node *T1, vector<Node *> &gene_trees, int start, int end, bool approx) {
	MAIN_CALL = false;
//	T1->preorder_number();
	vector<int> distances = vector<int>(end-start);
	#pragma omp parallel for shared(distances) firstprivate(PREFER_RHO)
	for(int i = start; i < end; i++) {
		int k;
		if (approx) {
			Forest F1 = Forest(T1);
			Forest F2 = Forest(gene_trees[i]);
			k = rSPR_worse_3_approx_distance_only(&F1, &F2)/3;
		}
		else {
			k = rSPR_branch_and_bound_simple_clustering(T1, gene_trees[i]);
		}
		distances[i-start] = k;
	}

	cout << distances[0];
	for(int i = 1; i < end-start; i++) {
		cout << "," << distances[i];
	}
	cout << "\n";
}


void rSPR_pairwise_distance(Node *T1, vector<Node *> &gene_trees, int max_spr) {
	rSPR_pairwise_distance(T1, gene_trees, max_spr, 0, (int)gene_trees.size());
}

void rSPR_pairwise_distance(Node *T1, vector<Node *> &gene_trees, int max_spr, int start, int end) {
	MAIN_CALL = false;
//	T1->preorder_number();
	vector<int> distances = vector<int>(end-start);
	#pragma omp parallel for shared(distances) firstprivate(PREFER_RHO)
	for(int i = start; i < end; i++) {
		Forest F1 = Forest(T1);
		Forest F2 = Forest(gene_trees[i]);
		int k = rSPR_branch_and_bound_range(&F1, &F2, 0, max_spr);
		distances[i-start] = k;
	}

	cout << distances[0];
	for(int i = 1; i < end-start; i++) {
		cout << "," << distances[i];
	}
	cout << "\n";
}

void rSPR_pairwise_distance_unrooted(Node *T1, vector<Node *> &gene_trees) {
	rSPR_pairwise_distance_unrooted(T1, gene_trees, 0, gene_trees.size());
}

void rSPR_pairwise_distance_unrooted(Node *T1, vector<Node *> &gene_trees, bool approx) {
	rSPR_pairwise_distance_unrooted(T1, gene_trees, 0, gene_trees.size(), approx);
}

void rSPR_pairwise_distance_unrooted(Node *T1, vector<Node *> &gene_trees, int start, int end) {
	rSPR_pairwise_distance_unrooted(T1, gene_trees, 0, gene_trees.size(), false);
}

void rSPR_pairwise_distance_unrooted(Node *T1, vector<Node *> &gene_trees, int start, int end, bool approx) {
	MAIN_CALL = false;
	T1->preorder_number();
	vector<int> distances = vector<int>(end-start);
	#pragma omp parallel for shared(distances) firstprivate(PREFER_RHO)
	for(int i = start; i < end; i++) {
		int best_k = INT_MAX;
		Node *T2_copy = new Node(*(gene_trees[i]));
		vector<Node *> descendants = 
				T2_copy->find_descendants();
		for(int j = 0; j < descendants.size(); j++) {
			T2_copy->reroot(descendants[j]);
			T2_copy->set_depth(0);
			T2_copy->fix_depths();
			T2_copy->preorder_number();
	//				cout << i << "," << j << endl;
	//				cout << T1->str_subtree() << endl;
	//				cout << gene_trees[i]->str_subtree() << endl;
			int k;
			if (approx) {
				Forest F1 = Forest(T1);
				Forest F2 = Forest(T2_copy);
				k = rSPR_worse_3_approx(&F1, &F2) / 3;
			}
			else {
				k = rSPR_branch_and_bound_simple_clustering(T1, T2_copy, false);
			}
			if (k < best_k) {
				best_k = k;
			}
		}
		distances[i-start] = best_k;
		T2_copy->delete_tree();
	}

	cout << distances[0];
	for(int i = 1; i < end-start; i++) {
		cout << "," << distances[i];
	}
	cout << "\n";
}

void rSPR_pairwise_distance_unrooted(Node *T1, vector<Node *> &gene_trees, int max_spr) {
	rSPR_pairwise_distance_unrooted(T1, gene_trees, max_spr, 0, (int)gene_trees.size());
}

void rSPR_pairwise_distance_unrooted(Node *T1, vector<Node *> &gene_trees, int max_spr, int start, int end) {
	MAIN_CALL = false;
	T1->preorder_number();
	vector<int> distances = vector<int>(end-start);
	#pragma omp parallel for shared(distances) firstprivate(PREFER_RHO)
	for(int i = start; i < end; i++) {
		int best_k = -1;
		Node *T2_copy = new Node(*(gene_trees[i]));
		vector<Node *> descendants = 
				T2_copy->find_descendants();
		for(int j = 0; j < descendants.size(); j++) {
			T2_copy->reroot(descendants[j]);
			T2_copy->set_depth(0);
			T2_copy->fix_depths();
			T2_copy->preorder_number();
	//				cout << i << "," << j << endl;
	//				cout << T1->str_subtree() << endl;
	//				cout << gene_trees[i]->str_subtree() << endl;
			Forest F1 = Forest(T1);
			Forest F2 = Forest(T2_copy);
			int k = rSPR_branch_and_bound_range(&F1, &F2, 0, max_spr);
			if ((best_k == -1) || (k < best_k && k >= 0)) {
				best_k = k;
			}
		}
		distances[i-start] = best_k;
		T2_copy->delete_tree();
	}

	cout << distances[0];
	for(int i = 1; i < end-start; i++) {
		cout << "," << distances[i];
	}
	cout << "\n";
}

int rSPR_total_distance_precomputed(Node *T1, vector<Node *> &gene_trees,
		vector<int> *original_scores, vector<int> *new_original_scores, Node *old_T1) {
	int total = 0;
	MAIN_CALL = false;
	int end = gene_trees.size();
//	T1->preorder_number();
	#pragma omp parallel for reduction(+ : total) firstprivate(PREFER_RHO)  // firstprivate(IN_SPLIT_APPROX)
	for(int i = 0; i < end; i++) {
		// check that the SPR move affects the projection of T1
		Forest F1 = Forest(T1);
		Forest F2 = Forest(gene_trees[i]);
		Forest F1_old = Forest(old_T1);
		sync_twins(&F1, &F2);
		sync_twins(&F1, &F1_old);
		int k = 0;
		if (original_scores == NULL
				|| rSPR_worse_3_approx(&F1, &F1_old) > 0) {
			k = rSPR_branch_and_bound_simple_clustering(T1, gene_trees[i], VERBOSE);
		}
		else {
			k = (*original_scores)[i];
		}
		if (new_original_scores != NULL) {
			(*new_original_scores)[i] = k;
		}

		total += k;
	}
	return total;
}


int rf_total_distance(Node *T1, vector<Node *> &gene_trees) {
	int total = 0;
	int end = gene_trees.size();
	#pragma omp parallel for reduction(+ : total) firstprivate(PREFER_RHO)  // firstprivate(IN_SPLIT_APPROX)
	for(int i = 0; i < end; i++) {
			//		cout << i << endl;
		int k = rf_distance(T1, gene_trees[i]);
		total += k;
	}
	return total;
}

int rf_total_distance_unrooted(Node *T1, vector<Node *> &gene_trees) {
	int total = 0;
	int end = gene_trees.size();
	#pragma omp parallel for reduction(+ : total) firstprivate(PREFER_RHO)  // firstprivate(IN_SPLIT_APPROX)
	for(int i = 0; i < end; i++) {
		int best_k = INT_MAX;
		Node T2_copy = Node(*(gene_trees[i]));
		vector<Node *> descendants = 
				T2_copy.find_descendants();
		for(int j = 0; j < descendants.size(); j++) {
			T2_copy.reroot(descendants[j]);
			T2_copy.set_depth(0);
			T2_copy.fix_depths();
			T2_copy.preorder_number();
			int k = rf_distance(T1, &T2_copy);
			if (k < best_k) {
				best_k = k;
			}
		}
		total += best_k;
	}
	return total;
}



void rf_pairwise_distance(Node *T1, vector<Node *> &gene_trees) {
	rf_pairwise_distance(T1, gene_trees, 0, gene_trees.size());
}

void rf_pairwise_distance(Node *T1, vector<Node *> &gene_trees, int start, int end) {
	MAIN_CALL = false;
//	T1->preorder_number();
	vector<int> distances = vector<int>(end-start);
	#pragma omp parallel for shared(distances) firstprivate(PREFER_RHO)
	for(int i = start; i < end; i++) {
		int k = rf_distance(T1, gene_trees[i]);
		distances[i-start] = k;
	}

	cout << distances[0];
	for(int i = 1; i < end-start; i++) {
		cout << "," << distances[i];
	}
	cout << "\n";
}

void rf_pairwise_distance_unrooted(Node *T1, vector<Node *> &gene_trees) {
	rf_pairwise_distance_unrooted(T1, gene_trees, 0, gene_trees.size());
}

void rf_pairwise_distance_unrooted(Node *T1, vector<Node *> &gene_trees, int start, int end) {
	MAIN_CALL = false;
	T1->preorder_number();
	vector<int> distances = vector<int>(end-start);
	#pragma omp parallel for shared(distances) firstprivate(PREFER_RHO)
	for(int i = start; i < end; i++) {
		int best_k = INT_MAX;
		Node T2_copy = Node(*(gene_trees[i]));
		vector<Node *> descendants = 
				T2_copy.find_descendants();
		for(int j = 0; j < descendants.size(); j++) {
			T2_copy.reroot(descendants[j]);
			T2_copy.set_depth(0);
			T2_copy.fix_depths();
			T2_copy.preorder_number();
	//				cout << i << "," << j << endl;
	//				cout << T1->str_subtree() << endl;
	//				cout << gene_trees[i]->str_subtree() << endl;
			int k = rf_distance(T1, &T2_copy);
			if (k < best_k) {
				best_k = k;
			}
		}
		distances[i-start] = best_k;
	}

	cout << distances[0];
	for(int i = 1; i < end-start; i++) {
		cout << "," << distances[i];
	}
	cout << "\n";
}

int rSPR_total_distance(Node *T1, vector<Node *> &gene_trees, int threshold) {
	int total = 0;
	MAIN_CALL = false;
	int end = gene_trees.size();
	T1->preorder_number();
	#pragma omp parallel for reduction(+ : total) firstprivate(PREFER_RHO)  // firstprivate(IN_SPLIT_APPROX)
	for(int i = 0; i < end; i++) {
		int k = rSPR_branch_and_bound_simple_clustering(T1, gene_trees[i], VERBOSE);
//		k *= mylog2(gene_trees[i]->size());
		total += k;
//		if (total > threshold) {
//			break;
//		}
	}
	return total;
}

/*Joel's part*/
int rSPR_total_distance(Forest *T1, vector<Node *> &gene_trees){
	int total = 0;
	#pragma omp parallel for reduction(+ : total) firstprivate(PREFER_RHO)
	for(int i = 0; i < gene_trees.size(); i++) {
		Forest T2 = Forest(gene_trees[i]);
		total += rSPR_branch_and_bound_simple_clustering(&T2, T1, VERBOSE);
	}
	return total;
}

int rSPR_total_approx_distance(Forest *T1, vector<Node *> &gene_trees) {
	int total = 0;
	#pragma omp parallel for reduction(+ : total)
	for(int i = 0; i < gene_trees.size(); i++) {
		Forest F1 = Forest(T1);
		Forest F2 = Forest(gene_trees[i]);
//		cout << i << endl;
//		cout << T1->str_subtree() << endl;
//		cout << gene_trees[i]->str_subtree() << endl;
		//total += rSPR_worse_3_approx(&F2, &F1)/3;
		total += rSPR_worse_3_approx(&F2, &F1)/3;
	}
	return total;
}

int rSPR_total_distance_unrooted(Node *T1, vector<Node *> &gene_trees) {
	return rSPR_total_distance_unrooted(T1, gene_trees, INT_MAX, NULL);
}

int rSPR_total_distance_unrooted(Node *T1, vector<Node *> &gene_trees,
		int threshold) {
	return rSPR_total_distance_unrooted(T1, gene_trees, threshold, NULL);
}

int rSPR_total_distance_unrooted(Node *T1, vector<Node *> &gene_trees, int threshold, vector<int> *original_scores) {
	//cout << "rSPR_total_distance_unrooted" << endl;
	int total = 0;
	MAIN_CALL = false;
	T1->preorder_number();
	#pragma omp parallel for reduction(+ : total) firstprivate(PREFER_RHO) firstprivate(MAX_SPR) firstprivate(MIN_SPR)
	for(int i = 0; i < gene_trees.size(); i++) {
//		cout << "T1: " << T1->str_subtree() << endl;
//		cout << "T2: " << gene_trees[i]->str_subtree() << endl;
		Forest f1 = Forest(T1);
		//f1.print_components();
		Forest f2 = Forest(gene_trees[i]);
		//f2.print_components();
		if (!sync_twins(&f1, &f2))
			continue;
		if (f2.get_component(0)->get_children().size() > 2) {
			f2.get_component(0)->fixroot();
			f2.get_component(0)->set_depth(0);
			f2.get_component(0)->fix_depths();
			f2.get_component(0)->preorder_number();
		}
		//f1.print_components();
		//f2.print_components();
		int size = f2.get_component(0)->size();
		int best_distance = INT_MAX;
		int old_max = MAX_SPR;
		bool done = false;
		int NO_CLUSTER_ROUNDS=15;
//		cout << "boo" << endl;
		if (!UNROOTED_MIN_APPROX) {
//			for(int k = 0; !done; k++) {
			int best_min_spr = INT_MAX;
			vector<Node *> descendants = 
				f2.get_component(0)->find_descendants();
/*
			for(int j = 0; j < descendants.size(); j++) {
				f2.get_component(0)->reroot(descendants[j]);
				f2.get_component(0)->set_depth(0);
				f2.get_component(0)->fix_depths();
				f2.get_component(0)->preorder_number();
				Forest F1 = Forest(f1);
				Forest F2 = Forest(f2);
				int distance = rSPR_worse_3_approx(&F1, &F2)/3;
				if (distance < best_min_spr)
					best_min_spr = distance;
			}
			*/
		
			int min_spr = 0;
			if (best_min_spr < INT_MAX)
				min_spr = best_min_spr;
			for(int k = min_spr; k <= NO_CLUSTER_ROUNDS; k++) {
//			for(int k = min_spr; !done; k++) {
////				cout << k << endl;
				MIN_SPR=k;
				MAX_SPR=k;
//				Node *original_lc = f2.get_component(0)->lchild();
////					f2.print_components();
////					cout << endl;
//				vector<Node *> descendants = 
//				f2.get_component(0)->find_descendants();
				for(int j = 0; j < descendants.size(); j++) {
//					cout << "J=" << j << endl;
//					cout << i << "," << k << "," << j << endl;
//					cout << "rooting at: " << descendants[j]->str_subtree() << endl;
					//f2.get_component(0)->reroot(original_lc);
					f2.get_component(0)->reroot(descendants[j]);
					f2.get_component(0)->set_depth(0);
					f2.get_component(0)->fix_depths();
					f2.get_component(0)->preorder_number();
////					f2.print_components();
////					cout << endl;
	//			cout << T1->str_subtree() << endl;
	//			cout << gene_trees[i]->str_subtree() << endl;
					int distance;
					Forest *F1 = new Forest(f1);
					Forest *F2 = new Forest(f2);
					if (k <= NO_CLUSTER_ROUNDS)
						distance = rSPR_branch_and_bound_range(F1, F2, MIN_SPR, MAX_SPR);
//					else
//						break;
//						distance = rSPR_branch_and_bound_simple_clustering(F1->get_component(0), F2->get_component(0), VERBOSE, k, k);
					if (distance < 0)
						distance = k+1;
					delete F1;
					delete F2;
					if (distance <= k) {
						best_distance = distance;
						k=NO_CLUSTER_ROUNDS+1;
						done = true;
						break;
					}
				}
////				cout << endl;
				//f2.get_component(0)->reroot(original_lc);
//				cout << endl;
			}
			MAX_SPR=old_max;
			MIN_SPR=0;
			if (!done) {
				vector<Node *> descendants = 
					f2.get_component(0)->find_descendants();
				for(int j = 0; j < descendants.size(); j++) {
					f2.get_component(0)->reroot(descendants[j]);
					f2.get_component(0)->set_depth(0);
					f2.get_component(0)->fix_depths();
					f2.get_component(0)->preorder_number();
	//				cout << i << "," << j << endl;
	//				cout << T1->str_subtree() << endl;
	//				cout << gene_trees[i]->str_subtree() << endl;
					int distance = rSPR_branch_and_bound_simple_clustering(f1.get_component(0), f2.get_component(0), VERBOSE);
					if (distance <= best_distance) {
							best_distance = distance;
					}
				}
			}
	//		cout << "best_distance: " << best_distance << endl;
			if (best_distance == INT_MAX)
				best_distance = 0;
			total += best_distance;
			if (original_scores != NULL)
				(*original_scores)[i] = best_distance;
	//		cout << "total: " << total << endl;
		}
		else {
			int best_approx = INT_MAX;
			Node *best_rooting = f2.get_component(0)->lchild();
			int num_ties = 2;
			vector<Node *> descendants = 
				f2.get_component(0)->find_descendants();
			int NUM_ROOTINGS = 0;
//			int NUM_ROOTINGS = 6;
//			if (descendants.size() > 70)
//				NUM_ROOTINGS = descendants.size() / 10;
//			NUM_ROOTINGS = sqrt(descendants.size());
			if (NUM_ROOTINGS > 0 && descendants.size() < NUM_ROOTINGS + 1) {
				vector<Node *> rand_descendants =
					random_select(descendants, NUM_ROOTINGS);
				descendants = rand_descendants;
				descendants.push_back(f2.get_component(0)->lchild());
			}
			for(int j = 0; j < descendants.size(); j++) {
				f2.get_component(0)->reroot(descendants[j]);
					f2.get_component(0)->set_depth(0);
					f2.get_component(0)->fix_depths();
					f2.get_component(0)->preorder_number();
				//Forest F1 = Forest(f1);
				//Forest F2 = Forest(f2);
				int distance = rSPR_worse_3_approx_distance_only(&f1, &f2)/3;
				if (distance < best_approx) {
					best_approx = distance;
					best_rooting = descendants[j];
					num_ties = 2;
				}
				else if (distance == best_approx) {
					int r = rand();
					if (r < RAND_MAX/num_ties) {
						best_approx = distance;
						best_rooting = descendants[j];
					}
					num_ties++;
				}
			}
			f2.get_component(0)->reroot(best_rooting);
					f2.get_component(0)->set_depth(0);
					f2.get_component(0)->fix_depths();
					f2.get_component(0)->preorder_number();
			int k;
			if (best_approx > 20)
				k = rSPR_branch_and_bound_simple_clustering(f1.get_component(0), f2.get_component(0), VERBOSE);
			else
					k = rSPR_branch_and_bound_range(&f1, &f2, best_approx/3, best_approx);
			total += k;
		if (original_scores != NULL)
			(*original_scores)[i] = k;
		}
//		if (total > threshold)
//			break;
	}
	return total;
}

int rSPR_total_approx_distance_unrooted(Node *T1, vector<Node *> &gene_trees) {
	int total = 0;
	MAIN_CALL = false;
	#pragma omp parallel for reduction(+: total)
	for(int i = 0; i < gene_trees.size(); i++) {
		Forest f1 = Forest(T1);
		Forest f2 = Forest(gene_trees[i]);
		if (!sync_twins(&f1, &f2))
			continue;
		if (f2.get_component(0)->get_children().size() > 2) {
			f2.get_component(0)->fixroot();
			f2.get_component(0)->set_depth(0);
			f2.get_component(0)->fix_depths();
			f2.get_component(0)->preorder_number();
		}
		int size = f2.get_component(0)->size();
		int best_distance = INT_MAX;
		vector<Node *> descendants = 
			f2.get_component(0)->find_descendants();
		for(int j = 0; j < descendants.size(); j++) {
			f2.get_component(0)->reroot(descendants[j]);
					f2.get_component(0)->set_depth(0);
					f2.get_component(0)->fix_depths();
					f2.get_component(0)->preorder_number();
			Forest F1 = Forest(f1);
			Forest F2 = Forest(f2);

			int distance = rSPR_worse_3_approx(&F1, &F2)/3;
			if (distance < best_distance)
				best_distance = distance;
		}
		if (best_distance == INT_MAX)
			best_distance = 0;
		total += best_distance;
	}
	return total;
}

int rSPR_total_approx_distance(Node *T1, vector<Node *> &gene_trees) {
	return rSPR_total_approx_distance(T1, gene_trees, INT_MAX);
}

int rSPR_total_approx_distance(Node *T1, vector<Node *> &gene_trees,
		int threshold) {
	int total = 0;
	MAIN_CALL = false;
	#pragma omp parallel for reduction(+ : total)
	for(int i = 0; i < gene_trees.size(); i++) {
		Forest F1 = Forest(T1);
		Forest F2 = Forest(gene_trees[i]);
//		cout << i << endl;
//		cout << T1->str_subtree() << endl;
//		cout << gene_trees[i]->str_subtree() << endl;
		total += rSPR_worse_3_approx(&F1, &F2)/3;
//		if (total > threshold)
//			break;
	}
	return total;
}


string itos(int i) {
	stringstream ss;
	string a;
	ss << i;
	a = ss.str();
	return a;
}

Node *find_subtree_of_approx_distance_hlpr(Node *n, Forest *F1, Forest *F2, int target_size) {
	Node *largest_child_subtree = NULL;
	int lcs_size = 0;
	list<Node *>::iterator c;
	for(c = n->get_children().begin(); c != n->get_children().end(); c++) {
		Forest f1 = Forest(F1);
		Forest f2 = Forest(F2);
		Node *subtree = f1.find_by_prenum((*c)->get_preorder_number());
		f1.get_component(0)->disallow_siblings_subtree();
		subtree->allow_siblings_subtree();
//		if (subtree->lchild() != NULL)
//			subtree->lchild()->allow_siblings_subtree();
//		if (subtree->rchild() != NULL)
//			subtree->rchild()->allow_siblings_subtree();

		int cs_size = rSPR_worse_3_approx(subtree, &f1, &f2);
		if (cs_size > lcs_size) {
			largest_child_subtree = *c;
			lcs_size = cs_size;
		}
	}
	if (lcs_size <= 0)
		return n;
	else if (lcs_size < target_size)
		return largest_child_subtree;
	else
		return find_subtree_of_approx_distance_hlpr(largest_child_subtree,
				F1, F2, target_size);
}

Node *find_subtree_of_approx_distance(Node *n, Forest *F1, Forest *F2, int target_size) {
		Forest f1 = Forest(F1);
		Forest f2 = Forest(F2);
		Node *subtree = f1.find_by_prenum(n->get_preorder_number());
		f1.get_component(0)->disallow_siblings_subtree();
		subtree->allow_siblings_subtree();
//		if (subtree->lchild() != NULL)
//			subtree->lchild()->allow_siblings_subtree();
//		if (subtree->rchild() != NULL)
//			subtree->rchild()->allow_siblings_subtree();
		int size = rSPR_worse_3_approx(subtree, &f1, &f2);
		if (size > target_size)
			return find_subtree_of_approx_distance_hlpr(n, F1, F2, target_size);
		else 
			return n;
}

Node *find_best_root(Node *T1, Node *T2, double *best_root_b_acc) {
	Forest F1 = Forest(T1);
	Forest F2 = Forest(T2);
	Node *t1 = F1.get_component(0);
	Node *t2 = F2.get_component(0);
	//F1->preorder_number();
	int lchild_pre = t1->lchild()->get_preorder_number();
	int rchild_pre = t1->rchild()->get_preorder_number();
	if (!sync_twins(&F1, &F2)) {
		return NULL;
	}
//	if (t1->lchild()->get_preorder_number() != lchild_pre ||
//			t1->rchild()->get_preorder_number() != rchild_pre)
//		return NULL;
	if (F2.get_component(0)->get_children().size() > 2)
		F2.get_component(0)->fixroot();
	// TODO: maybe stop if F2 is too small?
	int pre_separator = t1->lchild()->get_preorder_number();
	int group_1_total;
	int group_2_total;
	if (t1->rchild()->get_preorder_number() > pre_separator) {
		pre_separator = t1->rchild()->get_preorder_number();
		group_1_total = t1->lchild()->find_leaves().size();
		group_2_total = t1->rchild()->find_leaves().size();
	}
	else {
		group_1_total = t1->rchild()->find_leaves().size();
		group_2_total = t1->lchild()->find_leaves().size();
	}
//	cout << "g1_total: " << group_1_total << endl;
//	cout << "g2_total: " << group_2_total << endl;
	Node *best_root = t2->lchild();
	find_best_root_hlpr(t2, pre_separator, group_1_total, group_2_total,
			&best_root, best_root_b_acc);
//	cout << "t1: " << t1->str_subtree() << endl;
//	cout << "t2: " << t2->str_subtree() << endl;
//	cout << "best_root: " << best_root->str_subtree() << endl;
	best_root = T2->find_by_prenum(best_root->get_preorder_number());
//	cout << "T1: " << T1->str_subtree() << endl;
//	cout << "best_root: " << best_root->str_subtree() << endl;
//	T2->reroot(best_root);
//	cout << "T2: " << T2->str_subtree() << endl;
	return best_root;
}

Node *find_best_root(Node *T1, Node *T2) {
	double best_root_b_acc = 0;
	Forest f1 = Forest(T1);
	Forest f2 = Forest(T2);
	sync_twins(&f1, &f2);
	Node *new_root =
		find_best_root(f1.get_component(0), f2.get_component(0), &best_root_b_acc);
	if (new_root != NULL)
		new_root = T2->find_by_prenum(new_root->get_preorder_number());
	return new_root;
}

double find_best_root_acc(Node *T1, Node *T2) {
	double best_root_b_acc = -1;
	Forest f1 = Forest(T1);
	Forest f2 = Forest(T2);
	sync_twins(&f1, &f2);
	find_best_root(f1.get_component(0), f2.get_component(0), &best_root_b_acc);
	return best_root_b_acc;
}

void find_best_root_hlpr(Node *T2, int pre_separator, int group_1_total,
		int group_2_total, Node **best_root, double *best_root_b_acc) {
	list<Node*>::iterator c;
	int group_1_descendants = 0;
	int group_2_descendants = 0;
	int num_ties = 2;
	for(c = T2->get_children().begin(); c != T2->get_children().end(); c++) {
		find_best_root_hlpr(*c, pre_separator, group_1_total,
				group_2_total, best_root, best_root_b_acc,
				&group_1_descendants, &group_2_descendants, &num_ties);
	}
}

/*	class child_ba_comp {
		private:
			int g_1_total;
			int g_2_total;
			bool d;
		public:
			child_ba_comp(int group_1_total, int group_2_total, bool direction) {
				g_1_total = group_1_total;
				g_2_total = group_2_total;
				d = direction;
			}
			bool operator()(const pair<int,int> x,const pair<int,int> y) {
				if (d) {
					return ((x.first / x.second * (x.first + x.second)) - 
							(y.first / y.second * (y.first + y.second)));
				}
				else {
					return ((x.second / x.first * (x.first + x.second)) - 
							(y.second / y.first * (y.first + y.second)));
				}
			}
	};
	*/

void find_best_root_hlpr(Node *n, int pre_separator, int group_1_total,
		int group_2_total, Node **best_root, double *best_root_b_acc,
		int *p_group_1_descendants, int *p_group_2_descendants, int *num_ties) {
	list<Node*>::iterator c;
	int group_1_descendants = 0;
	int group_2_descendants = 0;
//	vector<pair<int,int> > children_splits = vector<pair<int,int>>();
	for(c = n->get_children().begin(); c != n->get_children().end(); c++) {
//		int g_1_desc = 0;
//		int g_2_desc = 0;
		find_best_root_hlpr(*c, pre_separator, group_1_total,
				group_2_total, best_root, best_root_b_acc,
				&group_1_descendants, &group_2_descendants, num_ties);
//				&g_1_desc, &g_2_desc, num_ties);
//		children_splits.push_back(make_pair(g_1_desc,g_2_desc));
//		group_1_descendants += g_1_desc;
//		group_2_descendants += g_2_desc;
	}
	if (n->is_leaf()) {
		int pre = n->get_twin()->get_preorder_number();
		if (pre < pre_separator)
			group_1_descendants++;
		else
			group_2_descendants++;
	}
//	else if (n->get_children.size() > 2) {
//		if (group_1_descendants / (double) group_1_total >= group_2_descendants / (double) group_2_total) {
//			sort(children_splits.begin(),children_splits.end(), g_1_comp(group_1_total, group_2_total, true));
//		}
//	}
	// balanced accuracy
	// don't bother averaging since we only directly compare them
	double tpos = group_1_descendants;
	double fpos = group_2_descendants;
	double fneg = (group_1_total - group_1_descendants);
	double tneg = (group_2_total - group_2_descendants);
	// balanced accuracy for this bipartition
	double b_acc =  tpos / (tpos + fneg)
			+ tneg / (tneg + fpos);
	// balanced accuracy for the opposite bipartition
	double b_acc_opp = fpos / (fpos + tneg)
			+ fneg / (fneg + tpos);
	// use the better bipartition
//	cout << "n: " << n->str_subtree() << endl;
//	cout << "b_acc: " << b_acc << endl;
//	cout << "b_acc_opp: " << b_acc_opp << endl;
//	cout << n->str_subtree() << endl;
//	cout << b_acc << "\t" << b_acc_opp << endl;
	if (b_acc_opp > b_acc)
		b_acc = b_acc_opp;
	if (b_acc > *best_root_b_acc) {
		*best_root = n;
		*best_root_b_acc = b_acc;
		*num_ties = 2;
	}
	else if(b_acc == *best_root_b_acc) {
		int r = rand();
		if (r < RAND_MAX/ *num_ties) {
			*best_root = n;
			*best_root_b_acc = b_acc;
		}
	}
	*p_group_1_descendants += group_1_descendants;
	*p_group_2_descendants += group_2_descendants;
}

Node *find_random_root(Node *T1, Node *T2) {
	vector<Node *> rroots = T2->find_descendants();
	int r = rand() % rroots.size();
	return rroots[r];
}
Node *find_best_root_rspr(Node *T1, Node *T2) {
	Node *t1 = new Node(*T1);
//	t1->preorder_number();
	Node *t2  = new Node(*T2);
	int new_prenum = T2->lchild()->get_preorder_number();
	vector<Node *> roots = t2->find_descendants();
	vector<int> root_prenums = vector<int>(roots.size());
	for(int i = 0; i < roots.size(); i++) {
		root_prenums[i] = roots[i]->get_preorder_number();
	}
	int best_distance = INT_MAX;
	int num_ties = 2;
//	cout << endl;
//	cout << "find_best_root_rspr" << endl;
//	cout << "T1: " << T1->str_subtree() << endl;
//	cout << "T2: " << T2->str_subtree() << endl;
//	cout << roots.size() << " Rootings" << endl;
	for(int i = 0; i < roots.size(); i++) {
		Node *root = roots[i];
		t2->reroot(root);
//		cout << "\t" << t2->str_subtree() << endl;
		t2->set_depth(0);
		t2->fix_depths();
		t2->preorder_number();
		int distance = rSPR_branch_and_bound_simple_clustering(t1, t2);
//		int distance = rf_distance(t1, t2);
		if (distance < best_distance) {
			best_distance = distance;
			new_prenum = root_prenums[i];
			num_ties = 2;
		}
		else if (distance == best_distance) {
			int r = rand();
			if (r < RAND_MAX / num_ties) {
				best_distance = distance;
				new_prenum = root_prenums[i];
			}
			num_ties++;
		}
	}
	Node *new_root = T2->find_by_prenum(new_prenum);
	t1->delete_tree();
	t2->delete_tree();
	return new_root;
}

// assume already sync_twins and preorder numbered
bool contains_bipartition(Node *n, int pre_start, int pre_end,
		int group_1_total, int group_2_total, int *p_group_1_descendants,
		int *p_group_2_descendants) {
	list<Node*>::iterator c;
	int group_1_descendants = 0;
	int group_2_descendants = 0;
	bool found = false;
	bool proper_split = true;
	for(c = n->get_children().begin(); c != n->get_children().end(); c++) {
		int c_group_1_descendants = 0;
		int c_group_2_descendants = 0;
		found = contains_bipartition(*c, pre_start, pre_end, group_1_total,
				group_2_total, &c_group_1_descendants, &c_group_2_descendants);
		if (found)
			return true;
		group_1_descendants += c_group_1_descendants;
		group_2_descendants += c_group_2_descendants;
		if (c_group_1_descendants > 0 && c_group_2_descendants > 0)
			proper_split = false;
	}

	if (n->is_leaf()) {
		int pre = n->get_twin()->get_preorder_number();
		if (pre >= pre_start && pre <= pre_end)
			group_1_descendants++;
		else
			group_2_descendants++;
	}
	else if (proper_split) {
		if (group_1_descendants == group_1_total
				|| group_2_descendants == group_2_total)
			return true;
	}
	if (p_group_1_descendants != NULL)
		*p_group_1_descendants += group_1_descendants;
	if (p_group_2_descendants != NULL)
		*p_group_2_descendants += group_2_descendants;
	return false;
}

// root the tree based on an outgroup
// returns false if the outgroup is not found or not a clade
bool outgroup_root(Node *T, set<string, StringCompare> outgroup) {
	vector<int> num_in = vector<int>();
	vector<int> num_out = vector<int>();
	count_in_out(T, num_in, num_out, outgroup);
	list<Node *>::iterator c;
	int pre = T->get_preorder_number();
	bool clean_split = true;
	if (num_out[pre] == 0)
		return false;
	for(c = T->get_children().begin(); c != T->get_children().end(); c++) {
		int c_pre = (*c)->get_preorder_number();
		if (num_out[pre] == num_out[c_pre]) {
			// split the out_group
			return outgroup_root(*c, num_in, num_out);
		}
		else if (num_in[pre] == num_in[c_pre]) {
			// split the in_group
			return outgroup_root(*c, num_out, num_in);
		}
		else if (num_out[c_pre] > 0 && num_in[c_pre] > 0)
			clean_split = false;
	}
	if (clean_split)
		return outgroup_reroot(T, num_in, num_out);
	else
		return false;
}

bool outgroup_root(Node *n, vector<int> &num_in, vector<int> &num_out) {
	list<Node *>::iterator c;
	int pre = n->get_preorder_number();
	bool clean_split = true;
	for(c = n->get_children().begin(); c != n->get_children().end(); c++) {
		int c_pre = (*c)->get_preorder_number();
		if (num_out[pre] == num_out[c_pre]) {
			// split the out_group
			return outgroup_root(*c, num_in, num_out);
		}
		else if (num_out[c_pre] > 0 && num_in[c_pre] > 0)
			clean_split = false;
	}
	if (clean_split)
		return outgroup_reroot(n, num_in, num_out);
	else
		return false;
}

bool outgroup_reroot(Node *n, vector<int> &num_in, vector<int> &num_out) {
	Node *T = n->find_root();
	if (num_in[n->get_preorder_number()] == 0) {
		if (!n->is_leaf()) {
			T->reroot(n);
			T->set_depth(0);
			T->fix_depths();
			T->preorder_number();
		}
		return true;
	}
	Node *new_split = new Node("");
	n->add_child(new_split);
	list<Node *>::iterator c = n->get_children().begin();
	while(c != n->get_children().end()) {
		Node *child = *c;
		c++;
		int c_pre = child->get_preorder_number();
		if (c_pre >= 0 && num_in[c_pre] == 0) {
			//child->cut_parent();
			new_split->add_child(child);
		}
	}
	T->reroot(new_split);
	T->set_depth(0);
	T->fix_depths();
	T->preorder_number();
	return true;
}

void count_in_out(Node *n, vector<int> &num_in, vector<int> &num_out,
		set<string, StringCompare> &outgroup) {
	list<Node *>::iterator c;
	int pre = n->get_preorder_number();
	if (num_in.size() <= pre)
		num_in.resize(pre + 1, 0);
	if (num_out.size() <= pre)
		num_out.resize(pre + 1, 0);
	if (n->is_leaf()) {
		if (outgroup.find(n->get_name()) != outgroup.end()) {
				num_out[pre] = 1;
				num_in[pre] = 0;
		}
		else {
				num_out[pre] = 0;
				num_in[pre] = 1;
		}
	}
	else {
		for(c = n->get_children().begin(); c != n->get_children().end(); c++) {
			count_in_out(*c, num_in, num_out, outgroup);
			num_in[pre] += num_in[(*c)->get_preorder_number()];
			num_out[pre] += num_out[(*c)->get_preorder_number()];
		}
	}
}

void modify_bipartition_support(Node *T1, Node *T2, enum RELAXATION relaxed) {
	Forest F1 = Forest(T1);
	Forest F2 = Forest(T2);
	if (!sync_twins(&F1, &F2)) {
		return;
	}
	Node *n = F1.get_component(0);
	vector<int> *F1_descendant_counts = n->find_leaf_counts();
	modify_bipartition_support(n, &F1, &F2, T1, T2, F1_descendant_counts,
			relaxed);
	delete F1_descendant_counts;
}

void modify_bipartition_support(Node *n, Forest *F1, Forest *F2,
		Node *T1, Node *T2, vector<int> *F1_descendant_counts, enum RELAXATION relaxed) {
	if (n->is_leaf())
		return;
	list<Node *>::iterator c;
	for(c = n->get_children().begin(); c != n->get_children().end(); c++) {
		modify_bipartition_support(*c, F1, F2, T1, T2, F1_descendant_counts,
				relaxed);
	}
	Node *t = T1->find_by_prenum(n->get_preorder_number());
	if (t->parent() == NULL)
		return;
	int pre_start = n->get_edge_pre_start();
	int pre_end = n->get_edge_pre_end();
	int group_1_total = (*F1_descendant_counts)[n->get_preorder_number()];
	int group_2_total = (*F1_descendant_counts)[F1->get_component(0)->get_preorder_number()] - group_1_total;
	if (group_2_total >= 2) {
		if (contains_bipartition(F2->get_component(0), pre_start, pre_end,
				group_1_total, group_2_total, NULL, NULL)) {
			t->a_inc_support();
			t->a_inc_support_normalization();
			// relaxed
			if (false && relaxed == ALL_RELAXED) {
				int stop_pre = 0;
				if (n->parent() != NULL) {
					stop_pre = n->parent()->get_preorder_number();
				while ((t = t->parent()) != NULL
						&& t->get_preorder_number() != stop_pre)
					t->a_inc_support();
					t->a_inc_support_normalization();
				}
			}
		}
		else {
//			t->a_dec_support();
			t->a_inc_support_normalization();
			// relaxed
			if (relaxed == ALL_RELAXED || relaxed == NEGATIVE_RELAXED) {
				int stop_pre = 0;
				if (n->parent() != NULL) {
					stop_pre = n->parent()->get_preorder_number();
				while ((t = t->parent()) != NULL
						&& t->get_preorder_number() != stop_pre)
//					t->a_dec_support();
					t->a_inc_support_normalization();
				}
			}
		}
	}
}

int rf_distance(Node *T1, Node *T2) {
	Forest F1 = Forest(T1);
	Forest F2 = Forest(T2);
	if (!sync_twins(&F1, &F2))
		return 0;
	if (F1.get_component(0)->is_leaf())
		return 0;
	sync_interior_twins(&F1, &F2);
	int rf_d = 0;
	rf_d += count_differing_bipartitions(F1.get_component(0));
	rf_d += count_differing_bipartitions(F2.get_component(0));
	return rf_d;
}

int count_differing_bipartitions(Node *n) {
	//cout << "Start: " << n->str_subtree() << endl;
	int count = 0;
	list<Node *>::iterator c;
	for(c = n->get_children().begin(); c != n->get_children().end(); c++) {
		count += count_differing_bipartitions(*c);
	}
	if (n->get_twin() == NULL ||
//			n->get_depth() > n->get_twin()->get_twin()->get_depth())
			n != n->get_twin()->get_twin()) {
		count++;
			}
	return count;
}

bool is_nonbranching(Forest *T1, Forest *T2, Node *T1_a, Node *T1_c, Node *T2_a, Node *T2_c) {
	if ((T2_a->get_depth() < T2_c->get_depth()
			&& T2_c->parent() != NULL)
			|| T2_a->parent() == NULL) {
		swap(&T1_a, &T1_c);
		swap(&T2_a, &T2_c);
	}
	else if (T2_a->get_depth() == T2_c->get_depth()) {
		if (T2_a->parent() && T2_c->parent() &&
				(T2_a->parent()->get_depth() <
				T2_c->parent()->get_depth()
				//|| (T2_c->parent()->parent()
				//&& T2_c->parent()->parent() == T2_a->parent())
				)) {
		swap(&T1_a, &T1_c);
		swap(&T2_a, &T2_c);
		}
	}
	int num_protected = T2_a->is_protected() + T2_c->is_protected();
	if (T2_a->parent()->get_children().size() == 2)
		num_protected += T2_a->get_sibling()->is_protected();
	if (num_protected >= 2)
		return true;
	if (CUT_ONE_B) {
		if (T2_a->parent()->parent() == T2_c->parent()
			&& T2_c->parent() != NULL
			&& T2_a->parent()->get_children().size() <= 2)
			return true;
	}
	if ((CUT_TWO_B) && T1_a->parent()->parent() != NULL) {
		Node *T1_s = T1_a->parent()->get_sibling();
		if (T1_s->is_leaf()) {
			Node *T2_l = T2_a->parent()->parent();
			if (T2_l != NULL && T2_l->get_children().size() <= 2) {
				if (T2_c->parent() != NULL && T2_c->parent()->parent() == T2_l
						&& ((T2_a->parent()->get_children().size() <= 2
						&& T2_c->parent()->get_children().size() <= 2)
						|| T1_s->get_twin()->is_protected())){
					if (T2_l->get_sibling() == T1_s->get_twin()) {
						return true;
					}
					else if (CUT_TWO_B_ROOT && T2_l->parent() == NULL &&
							(T2->contains_rho() ||
							 T2->get_component(0) != T2_l)) {
						return true;
					}
				}
				else if ((T2_l = T2_l->parent()) != NULL
						&& T2_c->parent() == T2_l
						&& ((T2_a->parent()->get_children().size() <= 2
						&& T2_a->parent()->parent()->get_children().size() <= 2
						&& T2_l->get_children().size() <= 2)
						|| T1_s->get_twin()->is_protected())){
					if (T2_l->get_sibling() == T1_s->get_twin()) {
						return true;
					}
					else if (CUT_TWO_B_ROOT && T2_l->parent() == NULL &&
							(T2->contains_rho() ||
							 T2->get_component(0) != T2_l)) {
						return true;
					}
				}
			}
		}
	}
	if (REVERSE_CUT_ONE_B && T1_a->parent()->parent() != NULL) {
		Node *T1_s = T1_a->parent()->get_sibling();
		Node *T2_s = T1_s->get_twin();
		if (T1_s->is_leaf()) {
			if (T2_s->parent() == T2_a->parent()) {
				return true;
			}
			else if (T2_s->parent() == T2_c->parent()
					&& T2_c->parent()->get_children().size() <= 2) {
				return true;
			}
//			else if (REVERSE_CUT_ONE_B_3
//							&& T2_s->is_protected()
//							&& T2_s->parent() != NULL
//							&& T2_s->parent()->parent() == T2_a->parent()
//							&& T2_s->parent()->get_children().size() <= 2) {
//				return true;
//			}
		}
		else if (REVERSE_CUT_ONE_B_2 && T2_c->parent() != NULL
				&& chain_match(T1_s, T2_c->get_sibling(), T2_a))
			return true;
	}
	return false;
}

void strip_whitespace(string &str) {
	std::string::iterator end_pos = std::remove_if(str.begin(), str.end(), ::isspace);
	str.erase(end_pos, str.end());
}

void strip_trailing_whitespace(string &str) {
	size_t start_pos = str.find_first_not_of(whitespaces);
	str.erase(0, start_pos);
	size_t end_pos = str.find_last_not_of(whitespaces)+1;
	str.erase(end_pos, str.size()-end_pos);
}

//randomizes T2 with count number of sprs. If T1 equals T2 at the beginning,
//then the spr distance between T1 and T2 is equal to (or possibly less than) count
//Assumes they are already synced, assumes there are valid sprs to be made
void randomize_tree_with_spr(Forest* T1, Forest* T2, int count) {
  for (int spr = 0; spr < count; spr++) {
    vector<Node*> all_nodes = T2->get_component(0)->find_nodes_in_subtree();
    Node* source = all_nodes[rand() % all_nodes.size()];
    Node* target = all_nodes[rand() % all_nodes.size()];
    bool target_in_subtree = false;
    Node* first_leaf = target->find_leaves()[0];
    vector<Node*> source_leaves = source->find_leaves();
    for (int i = 0; i < source_leaves.size(); i++) {
      if (source_leaves[i] == first_leaf) {
	target_in_subtree = true;
	break;
      }
    }

    //Get random node
    //if target is in source's subtree repick
    //spr
    while (source == target ||
	   source->is_sibling_of(target) ||
	   source->parent() == target || 
	   target_in_subtree) {	   
      source = all_nodes[rand() % all_nodes.size()];
      target = all_nodes[rand() % all_nodes.size()];
      //cout << "Trying: " << source->str_subtree() << " and "<<  target->str_subtree() << endl;
      target_in_subtree = false;
      Node* first_leaf = target->find_leaves()[0];
      
      int child_count = source->get_children().size();
      if (child_count > 2) {
	int rand_count = rand() % (child_count + 1);
	if (rand_count == child_count || rand_count == 0){
	  //cout << "moving whole tree" << endl;
	}
	else if (rand_count == 1) {
	  source = source->get_children().front();
	  //cout << "moving first child" << endl;
	}
	else {
	  list<Node*> to_expand = list<Node*>();
	  list<Node*>::iterator c = source->get_children().begin();
	  for (int i = 0; i < rand_count; i++) {
	    to_expand.push_back(*c);
	    c++;
	  }	
	  source = source->expand_children_out(to_expand);
	  //cout << "Moving part " << rand_count<< endl;
	}
      }
      
      vector<Node*> source_leaves = source->find_leaves();
      for (int i = 0; i < source_leaves.size(); i++) {
	if (source_leaves[i] == first_leaf) {
	  //cout << "target in subtree" << endl;
	  target_in_subtree = true;
	  break;
	}
      }
    } 

    //cout << "Moving : " << source->str_subtree() << " to " << target->str_subtree() << endl;
    
    Node* parent = source->parent();      
    if (parent != NULL) {
      source->cut_parent();
      if (parent->get_children().size() == 1) {
	parent->contract(true);
      }
    }
    else {
      continue;
    }
    //regraft
    if (target->parent() != NULL) {
      target->parent()->add_child(source);
    }
    else {
      Node* new_parent = target;//new Node();
      Node* replace = new Node(*target);
      new_parent->get_children().clear();
      new_parent->add_child(replace);
      new_parent->add_child(source);
    }

    //T1->print_components();
    //T2->print_components();
  }
}
