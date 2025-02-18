/*******************************************************************************
Node.h

Data structure for a node of a binary tree
Contains methods to recursively work on a node's subtree

Copyright 2009-2014 Chris Whidden
cwhidden@dal.ca
http://kiwi.cs.dal.ca/Software/RSPR
March 3, 2014
Version 1.2.1

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

#ifndef INCLUDE_NODE

#define INCLUDE_NODE

#define COPY_CONTRACTED
#define DEBUG_SUPPORT 0

#include <cstdio>
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <map>
#include <set>
#include "Forest.h"

using namespace std;

bool IGNORE_MULTI = false;
double REQUIRED_SUPPORT = 0.0;
double SUPPORT_BY_VALUE = 1.0;
double MIN_LENGTH = -1.0;

struct StringCompare {
	bool operator() (const string &a, const string &b) const {
		return strcmp(a.c_str(), b.c_str()) < 0;
	}
};


// representation of a component with no leaves
#define DEAD_COMPONENT "*"
//#define DEBUG_PROTECTED "@"
/*void find_sibling_pairs_hlpr(Node *node, list<Node *> &sibling_pairs);
void find_leaves_hlpr(Node *node, vector<Node *> &leaves);
void str_subtree_hlpr(string *s);
list<Node *> find_sibling_pairs(Node *node);
vector<Node *> find_leaves(Node *node);
*/
int stomini(string s);


class Forest;

class Node {
	private:
	//Node *lc;			// left child
	//Node *rc;			// right child
	list<Node *> children; // children
	Node *p;			// parent
	list<Node *>:: iterator p_link;		// location in parents list
	Node *twin;			// counterpart in another tree
	string name;		// label
	int depth;			//distance from root
	int pre_num;	// preorder number
	int edge_pre_start;
	int edge_pre_end;

	int component_number;
	list <Node *> active_descendants;
	list <Node *> root_lcas;
	list<list <Node *>::iterator> removable_descendants;
	list <Node *>::iterator sibling_pair_loc;
	int sibling_pair_status;
	int num_clustered_children;
	Forest *forest;
        list <Node *> contracted_children;
	Node *contracted_lc;
	Node *contracted_rc;
	bool contracted;
	bool edge_protected;
	bool edge_protected_forcefully;
	int max_merge_depth;
	bool allow_sibling;
	int lost_children;
	double support;
	double support_normalization;
        // TODO (Ben): put initialization in correct constructor
        int non_leaf_children = 0;

	public:
	Node() {
		init(NULL, NULL, NULL, "", 0);
	}
	Node(string n) {
		init(NULL, NULL, NULL, n, 0);
	}
	Node(string n, int d) {
		init(NULL, NULL, NULL, n, d);
	}
	Node(Node *lc, Node *rc, Node *p, string n, int d) {
		init(lc, rc, p, n, d);
	}
	void init(Node *lc, Node *rc, Node *p, string n, int d) {
//		this->lc = lc;
//		this->rc = rc;
		this->p = p;
		this->name = string(n);
		this->twin = NULL;
		this->depth = d;
		this->pre_num = -1;
		this->edge_pre_start = -1;
		this->edge_pre_end = -1;
		this->component_number = -2;
		this->active_descendants = list <Node *>();
		this->root_lcas = list <Node *>();
		this->removable_descendants = list< list<Node *>::iterator>();
		this->sibling_pair_loc = list<Node *>::iterator(); 
		this->sibling_pair_status = 0;
		this->num_clustered_children = 0;
		this->forest = NULL;
		this->contracted_lc = NULL;
		this->contracted_rc = NULL;
		this->contracted = false;
		this->edge_protected = false;
		this->edge_protected_forcefully = false;
		this->allow_sibling = true;
		this->lost_children = 0;
		this->max_merge_depth = -1;
		this->support = -1;
		this->support_normalization = -1;
		this->children = list<Node *>();
		if (lc != NULL)
			add_child(lc);
		if (rc != NULL)
			add_child(rc);
	}
	// copy constructor
	Node(const Node &n) {
		p = NULL;
		name = n.name.c_str();
		twin = n.twin;
		depth = n.depth;
//		depth = 0;
		pre_num = n.pre_num;
		edge_pre_start = n.edge_pre_start;
		edge_pre_end = n.edge_pre_end;
		component_number = n.component_number;
		this->active_descendants = list <Node *>();
		this->removable_descendants = list< list<Node *>::iterator>();
		this->root_lcas = list <Node *>();
		//sibling_pair_loc = n.sibling_pair_loc;
		//sibling_pair_status = n.sibling_pair_status;
		this->sibling_pair_loc = list<Node *>::iterator(); 
		this->sibling_pair_status = 0;
		this->num_clustered_children = 0;
		this->forest = NULL;
		list<Node *>::const_iterator c;
		this->children = list<Node *>();
		for(c = n.children.begin(); c != n.children.end(); c++) {
			add_child(new Node(**c));
		}
#ifdef COPY_CONTRACTED
		if (n.contracted_lc == NULL)
			contracted_lc = NULL;
		else
			contracted_lc = new Node(*(n.contracted_lc), this);
		if (n.contracted_rc == NULL)
			contracted_rc = NULL;
		else
			contracted_rc = new Node(*(n.contracted_rc), this);
		this->contracted = n.contracted;
		for (auto i = n.contracted_children.begin(); i != n.contracted_children.end(); i++) {
		  Node* contracted_node = new Node(**i, this);
		  add_contracted_child(contracted_node);
		  contracted_node->set_parent(this);
		} 
#else
		this->contracted_lc = n.contracted_lc;
		this->contracted_rc = n.contracted_rc;
		this->contracted = n.contracted;
		for (auto i = n.contracted_children.begin(); i != n.contracted_children.end(); i++) {
		  add_contracted_child(**i);
		  contracted_node->set_parent(this);
		} 

#endif
		this->edge_protected = n.edge_protected;
		this->edge_protected_forcefully = n.edge_protected_forcefully;
		this->allow_sibling = n.allow_sibling;
		this->lost_children = n.lost_children;
		this->max_merge_depth = n.max_merge_depth;
		this->support = n.support;
		this->support_normalization = n.support_normalization;
		this->non_leaf_children = n.non_leaf_children;
	}

        /* Temporary until rSPR_branch_and_bound_mult_hlpr uses the UndoMachine */
        // copy constructor
        Node(const Node &n, map<Node*, Node*> *node_map) {
		p = NULL;
		name = n.name.c_str();
//		twin = n.twin;
		depth = n.depth;
//		depth = 0;
		pre_num = n.pre_num;
		edge_pre_start = n.edge_pre_start;
		edge_pre_end = n.edge_pre_end;
		component_number = n.component_number;
		this->active_descendants = list <Node *>();
		this->removable_descendants = list< list<Node *>::iterator>();
		this->root_lcas = list <Node *>();
		//sibling_pair_loc = n.sibling_pair_loc;
		//sibling_pair_status = n.sibling_pair_status;
		this->sibling_pair_loc = list<Node *>::iterator(); 
		this->sibling_pair_status = 0;
		this->num_clustered_children = 0;
		this->forest = NULL;
		list<Node *>::const_iterator c;
		this->children = list<Node *>();
		for(c = n.children.begin(); c != n.children.end(); c++) {
		        Node* new_node = new Node(**c, node_map);
			add_child(new_node);
			node_map->emplace(pair<Node*,Node*>(*c, new_node));
		}
#ifdef COPY_CONTRACTED
		if (n.contracted_lc == NULL)
			contracted_lc = NULL;
		else {
		  contracted_lc = new Node(*(n.contracted_lc), this);// node_map);
		  // contracted_lc->set_parent(this);
		  }
		if (n.contracted_rc == NULL)
			contracted_rc = NULL;
		else {
		  contracted_rc = new Node(*(n.contracted_rc), this);//node_map);
			//contracted_rc->set_parent(this);
		}
		this->contracted = n.contracted;
		for (auto i = n.contracted_children.begin(); i != n.contracted_children.end(); i++) {
		  Node* contracted_node = new Node(**i, node_map);
		  contracted_children.push_back(contracted_node);
		  contracted_node->set_parent(this);
		  //TODO: why does this change distance?
		  //node_map->emplace(pair<Node*,Node*>(*i, contracted_node));
		} 
#else
		this->contracted_lc = n.contracted_lc;
		this->contracted_rc = n.contracted_rc;
		this->contracted = n.contracted;
#endif
		this->edge_protected = n.edge_protected;
		this->edge_protected_forcefully = n.edge_protected_forcefully;
		this->allow_sibling = n.allow_sibling;
		this->lost_children = n.lost_children;
		this->max_merge_depth = n.max_merge_depth;
		this->support = n.support;
		this->support_normalization = n.support_normalization;
		this->non_leaf_children = n.non_leaf_children;
	}

  
	Node(const Node &n, Node *parent) {
		p = parent;
		name = n.name.c_str();
		twin = n.twin;
		if (p != NULL)
			depth = p->depth+1;
		else
		depth = n.depth;
		pre_num = n.pre_num;
		edge_pre_start = n.edge_pre_start;
		edge_pre_end = n.edge_pre_end;
		component_number = n.component_number;
		this->active_descendants = list <Node *>();
		this->removable_descendants = list< list<Node *>::iterator>();
		this->root_lcas = list <Node *>();
		//sibling_pair_loc = n.sibling_pair_loc;
		//sibling_pair_status = n.sibling_pair_status;
		this->sibling_pair_loc = list<Node *>::iterator(); 
		this->sibling_pair_status = 0;
		this->num_clustered_children = 0;
		this->forest = NULL;
		this->children = list<Node *>();
		list<Node *>::const_iterator c;
		for(c = n.children.begin(); c != n.children.end(); c++) {
			add_child(new Node(**c));
		}
#ifdef COPY_CONTRACTED
		if (n.contracted_lc == NULL)
			contracted_lc = NULL;
		else
			contracted_lc = new Node(*(n.contracted_lc), this);
		if (n.contracted_rc == NULL)
			contracted_rc = NULL;
		else
			contracted_rc = new Node(*(n.contracted_rc), this);
		this->contracted = n.contracted;
		for (auto i = n.contracted_children.begin(); i != n.contracted_children.end(); i++) {
		  Node* contracted_node = new Node(**i, this);
		  add_contracted_child(contracted_node);
		  contracted_node->set_parent(this);
		} 

#else
		this->contracted_lc = n.contracted_lc;
		this->contracted_rc = n.contracted_rc;
		this->contracted = n.contracted;
		for (auto i = n.contracted_children.begin(); i != n.contracted_children.end(); i++) {
		  add_contracted_child(**i);
		  contracted_node->set_parent(this);
		} 

#endif
		this->edge_protected = n.edge_protected;
		this->edge_protected_forcefully = n.edge_protected_forcefully;
		this->allow_sibling = n.allow_sibling;
		this->lost_children = n.lost_children;
		this->max_merge_depth = n.max_merge_depth;
		this->support = n.support;
		this->support_normalization = n.support_normalization;
	}
	// TODO: clear_parent function
	~Node() {
		list<Node *>::iterator c = children.begin();
		while(c!= children.end()) {
			Node *n = *c;
			c++;
			n->cut_parent();
		}
		cut_parent();
		active_descendants.clear();
		root_lcas.clear();
		removable_descendants.clear();
#ifdef COPY_CONTRACTED
		if (contracted_lc != NULL) {
			contracted_lc->delete_tree();
		}
		contracted_lc = NULL;
		if (contracted_rc != NULL) {
			contracted_rc->delete_tree();
		}
		contracted_rc = NULL;
#endif
	}
	// TODO: is this still useful?
	/*
	void fake_delete() {
		if (lc != NULL) {
			lc->p = NULL;
			//lc = NULL;
		}
		if (rc != NULL) {
			rc->p = NULL;
			//rc = NULL;
		}
		if (p != NULL) {
			p->delete_child(this);
			//p = NULL;
		}
		active_descendants.clear();
		root_lcas.clear();
		removable_descendants.clear();
	}
	*/

	// delete a subtree
	void delete_tree() {
		list<Node *>::iterator c = children.begin();
		while(c!= children.end()) {
			Node *n = *c;
			c++;
			n->delete_tree();
		}

		cut_parent();
		/*
		children.clear();
		*/

#ifdef COPY_CONTRACTED
		if (contracted_lc != NULL) {
			contracted_lc->delete_tree();
			contracted_lc = NULL;
		}
		if (contracted_rc != NULL) {
			contracted_rc->delete_tree();
			contracted_rc = NULL;
		}
		for (auto i = contracted_children.begin(); i != contracted_children.end(); i++) {
		  (*i)->delete_tree();
		}
		contracted_children.clear();
#endif
		delete this;
	}

	// cut edge between parent and child
	// should really be cut_child, no deleting occurs
	// TODO: is this useful? The only reason for this would be to
	// be sure it works when the child is not correctly set
	void delete_child(Node *n) {
		n->cut_parent();
	}

	// TODO: make sure this doesn't break things with >2 children
	// add a child
	void add_child(Node *n) {
		if (n->p != NULL)
			n->cut_parent();
		n->p_link = children.insert(children.end(),n);
		n->p = this;
		n->depth = depth+1;
		n->contracted = false;
	}

	// TODO: make sure this doesn't break things with >2 children
	// add a child
	void add_child_keep_depth(Node *n) {
		if (n->p != NULL)
			n->cut_parent();
		n->p_link = children.insert(children.end(),n);
		n->p = this;
		n->contracted = false;
	}


	/* TODO: need new method of putting a child in a specific spot
	   in the children list

		 maybe adjacent to an iterator?
	*/

	// insert a child before the given sibling
	 void insert_child(Node *sibling, Node *n) {
		if (n->p != NULL)
			n->cut_parent();
		n->p_link = children.insert(sibling->p_link, n);
		n->depth = depth+1;
		n->p = this;
		n->contracted = false;
	 }

	// insert a child before the given sibling
	 void insert_child_keep_depth(Node *sibling, Node *n) {
		if (n->p != NULL)
			n->cut_parent();
	 	n->p_link = children.insert(sibling->p_link, n);
		n->p = this;
		n->contracted = false;
	 }


	/* dangerous, not relevant to multifurcating
	Node *set_lchild(Node *n) {
		this->lc = n;
		if (n != NULL) {
			n->p = this;
			n->depth = depth+1;
		}
		return lc;
	}

	Node *set_lchild_keep_depth(Node *n) {
		this->lc = n;
		if (n != NULL) {
			n->p = this;
		}
		return lc;
	}
	Node *set_rchild(Node *n) {
		this->rc = n;
		if (n != NULL) {
			n->p = this;
			n->depth = depth+1;
		}
		return rc;
	}
	Node *set_rchild_keep_depth(Node *n) {
		this->rc = n;
		if (n != NULL) {
			n->p = this;
		}
		return rc;
	}
*/
	// potentially dangerous
	Node *set_parent(Node *n) {
		p = n;
		return p;
	}
  
        void add_children(list<Node *> *nodes) {
	    list<Node *>::iterator i;
	    for (i = nodes->begin(); i != nodes->end(); i++) {
	      //children.push_back((*i));
	      add_child(*i);
	    }
	}
	Node *set_twin(Node *n) {
		twin = n;
		return twin;
	}
	void set_name(string n) {
		name = string(n);
	}
	int set_depth(int d) {
		depth = d;
		return depth;
	}
	void fix_depths() {
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->depth = depth+1;
			(*c)->fix_depths();
		}
	}
	int set_preorder_number(int p) {
		pre_num = p;
		return pre_num;
	}

	int set_edge_pre_start(int p) {
		edge_pre_start= p;
		return edge_pre_start;
	}

	int set_edge_pre_end(int p) {
		edge_pre_end= p;
		return edge_pre_end;
	}
        void set_non_leaf_children(int p) {
	        non_leaf_children = p;
	}

        void decrement_non_leaf_children() {
	        non_leaf_children--;
	}
        void increment_non_leaf_children() {
	        non_leaf_children++;
	}


	void copy_edge_pre_interval(Node *n) {
		if (n->edge_pre_start > -1) {
			edge_pre_start = n->edge_pre_start;
		}
		if (n->edge_pre_end > -1) {
			edge_pre_end = n->edge_pre_end;
		}
	}
	int set_component_number(int c) {
		component_number = c;
		return c;
	}
	list<Node *>& get_children() {
		return children;
	}
	Node *get_contracted_lc() {
		return contracted_lc;
	}
	Node *get_contracted_rc() {
		return contracted_rc;
	}

	void set_contracted_lc(Node *n) {
		contracted_lc = n;
	}
	void set_contracted_rc(Node *n) {
		contracted_rc = n;
	}


	bool is_protected() {
		return edge_protected;
	}
	bool is_contracted() {
		return contracted;
	}

	void protect_edge() {
		edge_protected = true;
		#ifdef DEBUG_LGT_EVENTS
			cout << "Protecting edge " << this->get_preorder_number() << endl;
			if(is_protected_forcefully()){
				cout << "Already protected forcefully "<< this->get_preorder_number() << endl;
			}
		#endif
	}

	void unprotect_edge() {
		#ifdef DEBUG_LGT_EVENTS
			cout << "Unrotecting edge " << this->get_preorder_number() << endl;
		#endif
		if(!is_protected_forcefully()){
			edge_protected = false;
		}
		else {
			#ifdef DEBUG_LGT_EVENTS
				cout << "Unprotection failed " << this->get_preorder_number() << endl;
			#endif
		}
	}

	void unprotect_subtree() {
		unprotect_edge();
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->unprotect_edge();
		}
	}

	void unprotect_tree() {
		unprotect_edge();
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->unprotect_tree();
		}
	}

	bool is_protected_forcefully() {
		return edge_protected_forcefully;
	}

	void protect_edge_forcefully() {
		#ifdef DEBUG_LGT_EVENTS
			cout << "Frocefully protecting edge " << this->get_preorder_number() << endl;
		#endif
		edge_protected = true;
		edge_protected_forcefully = true;
	}

	void unprotect_edge_forcefully() {
		edge_protected_forcefully = false;
		edge_protected = false;
	}

	void print_protected_edge_list_hlpr() {
		if(is_protected_forcefully()){
			cout << "Edge " << this->get_preorder_number() << " protected forcefully " << endl;
		}
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->print_protected_edge_list_hlpr();
		}
	}

	void protect_supported_edges() {
		if (support > 0)
			edge_protected = true;
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->protect_supported_edges();
		}
	}

	bool can_be_sibling() {
		return allow_sibling;
	}

	void disallow_siblings() {
		allow_sibling = false;
	}

	void disallow_siblings_subtree() {
		allow_sibling = false;
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->disallow_siblings_subtree();
		}
	}

	void allow_siblings() {
		allow_sibling = true;
	}

	void allow_siblings_subtree() {
		allow_sibling = true;
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->allow_siblings_subtree();
		}
	}

        bool is_sibling_of(Node* other) {
	  if (parent() != NULL && other->parent() != NULL) {	    
	    return parent() == other->parent();
	  }
	  return false;
	}


	int num_lost_children() {
		return lost_children;
	}

	int get_max_merge_depth() {
		return max_merge_depth;
	}
	void set_max_merge_depth(int d) {
		max_merge_depth = d;
	}


	int count_lost_children_subtree() {
		int lost_children_count = lost_children;
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			lost_children_count += (*c)->count_lost_children_subtree();
		}
		return lost_children_count;
	}
	int count_lost_subtree() {
		int lost_children_count = (lost_children > 0) ? 1 : 0;
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			lost_children_count += (*c)->count_lost_subtree();
		}
		return lost_children_count;
	}


	void lost_child() {
		lost_children++;
	}

	void no_lost_children() {
		lost_children = 0;
	}

	void no_lost_children_subtree() {
		lost_children = 0;
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->no_lost_children_subtree();
		}
	}

	double get_support() {
		return support;
	}
	void set_support(double s) {
		support = s;
	}
	void a_inc_support() {
#pragma omp atomic
		support += 1;
	}
	void a_dec_support() {
#pragma omp atomic
		support -= 1;
	}
	double get_support_normalization() {
		return support_normalization;
	}
	void set_support_normalization(double s) {
		support_normalization = s;
	}
	void a_inc_support_normalization() {
#pragma omp atomic
		support_normalization += 1;
	}
	void a_dec_support_normalization() {
#pragma omp atomic
		support_normalization -= 1;
	}

	void normalize_support() {
		if (support_normalization != 0)
			support /= support_normalization;
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->normalize_support();
		}
	}

	int get_component_number() {
		return component_number;
	}
	void increase_clustered_children() {
		num_clustered_children++;
	}
	void decrease_clustered_children() {
		num_clustered_children--;
	}
	void set_num_clustered_children(int c) {
		num_clustered_children = c;
	}
	int get_num_clustered_children() {
		return num_clustered_children;
	}
	list <Node *> *get_active_descendants() {
		return &active_descendants;
	}
	list <Node *> *get_root_lcas(){
		return &root_lcas;
	}
	int get_sibling_pair_status(){
		return sibling_pair_status;
	}
	void set_sibling_pair_status(int s){
		sibling_pair_status = s;
	}
	void set_forest(Forest *f) {
		forest = f;
	}
	void set_forest_rec(Forest *f) {
		forest = f;
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->set_forest_rec(f);
		}
	}
	Forest *get_forest() {
		return forest;
	}
	list<list <Node *>::iterator> *get_removable_descendants() {
		return &removable_descendants;
	}
	void initialize_component_number(int value) {
		component_number = value;
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->initialize_component_number(value);
		}
	}
	void initialize_active_descendants(list <Node *> value) {
		active_descendants = value;
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->initialize_active_descendants(value);
		}
	}
	void initialize_root_lcas(list <Node *> value) {
		root_lcas = value;
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->initialize_root_lcas(value);
		}
	}
	void initialize_removable_descendants(list<list <Node *>::iterator> value) {
		removable_descendants = value;
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->initialize_removable_descendants(value);
		}
	}


	/* contract:
	 * if this node has a parent and one child then contract it out
	 * if this node has a parent and no child then contract it out
	 *   parent will have one child so contract it as well.
	 * if this node has no parent and one child then take the
	 *	child's children and contract it out
	 * return the first degree two parent found or NULL if there
	 * was no contraction
	 */
	// TODO: check handling twins for interior nodes
	Node *contract(bool remove) {
		Node *parent = p;
		Node *child;
		Node *ret = NULL;
		// contract out this node and give child to parent
		if (parent != NULL) {
		//cout << p->str_subtree() << endl;
			if (children.size() == 1) {
					child = children.front();
					if (this == parent->children.back()) {
						parent->add_child_keep_depth(child);
						child->set_depth(depth);
					}
					else {
						list<Node *>::iterator sib = p_link;
						sib++;
						Node *sibling = *sib;
						parent->insert_child_keep_depth(sibling, child);
						child->set_depth(depth);
					}
					child->copy_edge_pre_interval(this);
					if (edge_protected && !child->is_protected())
						child->protect_edge();
					cut_parent();
					if (remove)
						delete this;
					ret = parent;
			}
			else if (children.empty()) {
				cut_parent();
				ret = parent->contract(remove);
				if (remove)
					delete this;
				//this->fake_delete();
			}
			else
				ret = this;

		}
		// if no parent then take children of single child and remove it
		else {

			// dead component or singleton, will be cleaned up by the forest
			if (children.empty()) {
				if (str() == "")
					name = DEAD_COMPONENT;
			}
			if (children.size() == 1) {
					child = children.front();
					child->cut_parent();

				/* cluster hack - if we delete a cluster node then
				 * we may try to use it later. This only happens once
				 * per cluster so we can spend linear time to update
				 * the forest
				 */
				if (child->num_clustered_children > 0) {
					// TODO: redundant?
					delete_child(child);
					if (remove)
						delete this;
					//this->fake_delete();
					ret = child;
				}
				else {
					// if child is a leaf then get rid of this so we don't lose refs
					// problem: if the child is not c, then we want to copy
					// otherwise we don't
					// copy other parameters and join the twin
					//to this if the child is a label
					Node *new_lc = child->lchild();
					Node *new_rc = child->rchild();
					if (child->is_leaf()) {
						if (child->get_twin() != NULL) {
							set_twin(child->get_twin());
							child->get_twin()->set_twin(this);
						}
						name = child->get_name().c_str();
//						name = child->str();
					}
					// TODO: redundant?
					child->cut_parent();
					list<Node *>::iterator c = child->children.begin();
					while(c!= child->children.end()) {
						Node *new_child = *c;
						c++;
						new_child->cut_parent();
						add_child(new_child);
					}
					if (child->contracted_lc != NULL)
						contracted_lc = child->contracted_lc;
					if (child->contracted_rc != NULL)
						contracted_rc = child->contracted_rc;
					pre_num = child->get_preorder_number();

					c = child->contracted_children.begin();
					while(c!= child->contracted_children.end()) {
						Node *new_child = *c;
						c++;
						new_child->set_parent(this);
						add_contracted_child(new_child);
						contracted = true;
					}
					child->contracted_children.clear();
					if (remove) {
						child->contracted_lc = NULL;
						child->contracted_rc = NULL;
						delete child;
					}
					ret = this;
				}
			}
		}

		return ret;
	}

	Node *contract() {
		return contract(false);
	}

	// TODO: binary only

	/* contract_sibling_pair:
	 * if this node has two child leaves then contract them out
	 * return true if contracted, otherwise false
	 */
	bool contract_sibling_pair() {
		if (lchild() != NULL && lchild()->is_leaf()
				&& rchild() != NULL && rchild()->is_leaf()) {
			#ifdef DEBUG
				string new_name = "<" + lchild()->str() + "," + rchild()->str() + ">";
			#else
				string new_name = "(" + lchild()->str() + "," + rchild()->str() + ")";
			#endif
			set_name(new_name);
			lchild()->cut_parent();
			rchild()->cut_parent();
			return true;
		}
		return false;
	}

	// TODO: binary only
	bool contract_sibling_pair_undoable() {
		if (lchild() != NULL && lchild()->is_leaf()
				&& rchild() != NULL && rchild()->is_leaf()) {
			/*
			#ifdef DEBUG
				string new_name = "<" + lc->str() + "," + rc->str() + ">";
			#else
				string new_name = "(" + lc->str() + "," + rc->str() + ")";
			#endif
			set_name(new_name);
			*/
			Node *lc = lchild();
			Node *rc = rchild();
			contracted_lc = lc;
			contracted_rc = rc;
			rc->cut_parent();
			lc->cut_parent();
			contracted_lc->contracted = true;
			contracted_rc->contracted = true;
			edge_protected = false;
			return true;
		}
		return false;
	}

	/* contract_sibling_pair_undoable
	 * works with multifurcating trees
	 * returns NULL for no contract, otherwise returns the contracted
	 * parent of the nodes
	 */
	Node *contract_sibling_pair_undoable(Node *child1, Node *child2) {

		if (child1->parent() != this ||
				child2->parent() != this)
			return NULL;
		if (children.size() == 2) {
			contract_sibling_pair_undoable();
			return this;
		}
		else {
			Node *new_child = new Node();
			// buggy
//			new_child->set_preorder_number(pre_num);
			if (child1->get_preorder_number() < child2->get_preorder_number()) {
				new_child->set_preorder_number(child1->get_preorder_number());
			}
			else {
				new_child->set_preorder_number(child2->get_preorder_number());
			}
			add_child(new_child);
			new_child->add_child(child1);
			new_child->add_child(child2);
			edge_protected = false;
			//new_child->contract_sibling_pair_undoable();
			return new_child;
		}
	}

	// TODO: binary only
	void undo_contract_sibling_pair() {
		// hacky, might hide problems
		if (contracted_lc != NULL)
			add_child(contracted_lc);
		if (contracted_rc != NULL)
			add_child(contracted_rc);
		contracted_lc = NULL;
		contracted_rc = NULL;
	}

        int recalculate_non_leaf_children() {
	        int total = 0;
		list<Node *>::iterator i;
		for (i = children.begin(); i != children.end(); i++) {
		  if (!(*i)->is_leaf()) {
			total++;
		  }
		}
		non_leaf_children = total;
		return total;
	}

        void add_contracted_child(Node *node) {
                contracted_children.push_back(node);
		node->contracted = true;
		node->set_parent(this);
	}

        list<Node *> *get_combined_children(list<Node *> *nodes) {
	        list<Node *> *combined = new list<Node *>();
		list<Node *>::iterator i;
		for (i = nodes->begin(); i != nodes->end(); i++) {
		  list<Node *> children = (*i)->get_children();
		  combined->insert(combined->end(), children.begin(), children.end());
		}
		return combined;
	}

  
        Node *contract_sibling_group(list<Node *> *nodes) {
	  //check to make sure all nodes are siblings
	  list<Node *>::iterator i;
	  for (i = nodes->begin(); i != nodes->end(); i++) {
	    if ((*i)->parent() != this) {
	      return NULL;
	    }
	  }
	  //contract all children into this
	  if (nodes->size() == children.size()) {
	    //TODO: Both have references to children, double free?
	    list<Node*> *new_children = get_combined_children(nodes);
	    children.clear();
	    for (i = new_children->begin(); i != new_children->end(); i++) {
	      children.push_back(*i);
	    }
	    delete new_children;
	    int min_pre_num = 0x7FFFFFFF; //max value constant?
	    for (i = nodes->begin(); i != nodes->end(); i++) {
	      //children.remove(*i); already removed in setting children
	      add_contracted_child(*i);
	      if ((*i)->get_preorder_number() < min_pre_num) {
		  min_pre_num = (*i)->get_preorder_number();
	      }
	    }
	    recalculate_non_leaf_children();
	    return this;
	  }
	  //Otherwise we create a new node and contract them into that
	  else {
	    Node *new_child = new Node();
	    add_child(new_child);
	    int min_pre_num = 0x7FFFFFFF; //max value constant?
	    for (i = nodes->begin(); i != nodes->end(); i++) {
	      new_child->add_contracted_child((*i));
	      children.remove((*i));
	      if ((*i)->get_preorder_number() < min_pre_num) {
		  min_pre_num = (*i)->get_preorder_number();
		}
	    }
	    new_child->set_preorder_number(min_pre_num);

	    
	    list<Node*> *new_children = get_combined_children(nodes);
	     new_child->add_children(new_children);
	    for (i = new_children->begin(); i != new_children->end(); i++) {
	      (*i)->set_parent(new_child);

	    }
	    delete new_children;
	    new_child->set_parent(this);
	    recalculate_non_leaf_children();
	    new_child->recalculate_non_leaf_children();
	    return new_child;
	  }
	}

    //Takes a list of nodes in the other tree to contract in this tree
    Node *contract_twin_group(list<Node *> *TO_nodes) {
	  list<Node *> TT_nodes;
	    list<Node *>::iterator i;
		for( i = TO_nodes->begin(); i != TO_nodes->end(); i++) {
		    TT_nodes.push_back((*i)->get_twin());			
		}
		return contract_sibling_group(&TT_nodes);		
    }
  
	void fix_contracted_order() {
		if (twin != NULL && twin->contracted_lc->twin != contracted_lc) {
			Node *swap = contracted_lc;
			contracted_lc = contracted_rc;
			contracted_rc = swap;
		}
	}


	// cut the edge between this node and its parent
	void cut_parent() {
		if (p != NULL) {
			// TODO hacky: fix this to use a multi list for contractions
			if (!contracted) {
				p->children.erase(p_link);
			}
			else {
				if (p->contracted_lc == this)
					p->contracted_lc = NULL;
				if (p->contracted_rc == this)
					p->contracted_rc = NULL;
			}
			p = NULL;
			p_link = children.end();
		}
	}

	// caution: destructive
	void contract_node() {
		if (p == NULL || is_leaf())
			return;
		list<Node *>::iterator c = children.begin();
		while(c!= children.end()) {
			Node *n = *c;
			c++;
			p->add_child(n);
		}

#ifdef COPY_CONTRACTED
		if (contracted_lc != NULL) {
			p->add_child(contracted_lc);
		}
		if (contracted_rc != NULL) {
			p->add_child(contracted_rc);
		}
#endif
		delete this;
	}

	Node *parent() {
		return p;
	}
	inline Node *lchild() {
		if (children.empty())
			return NULL;
		else
			return children.front();
	}

	Node *rchild() {
		if (children.empty())
				return NULL;
		list<Node *>::iterator c = ++(children.begin());
		if (c == children.end())
			return NULL;
		else
			return *c;
	}
	Node *get_twin() {
		return twin;
	}
	int get_depth() {
		return depth;
	}
	int get_preorder_number() {
		return pre_num;
	}
	int get_edge_pre_start() {
		return edge_pre_start;
	}
	int get_edge_pre_end() {
		return edge_pre_end;
	}
	string str() {
		string s = "";
		str_hlpr(&s);
		return s;
	}
	string get_name() {
		return name;
	}
  
	void str_hlpr(string *s) {
		if (!name.empty())
			*s += name.c_str();
		if (contracted_lc != NULL || contracted_rc != NULL || !contracted_children.empty()) {
			#ifdef DEBUG_CONTRACTED
				*s += "<";
			#else
				*s += "(";
			#endif
			if (contracted_lc != NULL) {
				contracted_lc->str_c_subtree_hlpr(s);
				*s += ",";
			}
			for (auto c = contracted_children.begin(); c != contracted_children.end(); c++) {
			  (*c)->str_c_subtree_hlpr(s);
			  *s += ",";
			}

			if (contracted_rc != NULL) {
				contracted_rc->str_c_subtree_hlpr(s);
			}
			else {
			  s->pop_back();
			}
			#ifdef DEBUG_CONTRACTED
				*s += ">";
			#else
				*s += ")";
			#endif
		}
	}

	string str_subtree() {
		string s = "";
		str_subtree_hlpr(&s);
		return s;
	}

	void str_subtree_hlpr(string *s) {
		str_hlpr(s);
		if (!is_leaf()) {
			*s += "(";
			list<Node *>::iterator c;
			for(c = children.begin(); c != children.end(); c++) {
				if (c != children.begin())
					*s += ",";
				(*c)->str_subtree_hlpr(s);
				if ((*c)->parent() != this)
					cout << "#";
			}
			*s += ")";
		}
#ifdef DEBUG_DEPTHS
		*s+= ":";
	stringstream ss;
	string a;
	ss << depth;
	a = ss.str();
		*s+= a;
#endif
#ifdef DEBUG_PROTECTED
		if (edge_protected)
			*s += DEBUG_PROTECTED;
#endif
	}

	string str_support_subtree(bool allow_negative) {
		string s = "";
		str_support_subtree_hlpr(&s, allow_negative);
		return s;
	}

	string str_support_subtree() {
		return str_support_subtree(false);
	}

	void str_support_subtree_hlpr(string *s, bool allow_negative) {
		str_hlpr(s);
		if (!is_leaf()) {
			*s += "(";
			list<Node *>::iterator c;
			for(c = children.begin(); c != children.end(); c++) {
				if (c != children.begin())
					*s += ",";
				(*c)->str_support_subtree_hlpr(s, allow_negative);
				if ((*c)->parent() != this)
					cout << "#";
			}
			*s += ")";
			if (get_support() > -1 || allow_negative) {
				stringstream ss;
				ss << setprecision (2) << get_support();
				*s+= ss.str();
				if (get_support_normalization() > -1 || allow_negative) {
					stringstream ss;
					ss << "#";
					ss << get_support_normalization();
					*s+= ss.str();
				}
			}
		}
	}

	string str_edge_pre_interval_subtree() {
		string s = "";
		str_edge_pre_interval_subtree_hlpr(&s);
		return s;
	}

	void str_edge_pre_interval_subtree_hlpr(string *s) {
		str_hlpr(s);
		if (!is_leaf()) {
			*s += "(";
			list<Node *>::iterator c;
			for(c = children.begin(); c != children.end(); c++) {
				if (c != children.begin())
					*s += ",";
				(*c)->str_edge_pre_interval_subtree_hlpr(s);
				if ((*c)->parent() != this)
					cout << "#";
			}
			*s += ")";
		}
		if (get_preorder_number() > -1) {
			stringstream ss;
			ss << ":";
			if (get_edge_pre_start() > -1) {
				ss << get_edge_pre_start();
				ss << ";";
			}
			ss << get_preorder_number();
			if (get_edge_pre_end() > -1) {
				ss << ";";
				ss << get_edge_pre_end();
			}
			*s+= ss.str();
		}
	}

	void str_c_subtree_hlpr(string *s) {
		str_hlpr(s);
		if (!is_leaf()) {
			#ifdef DEBUG_CONTRACTED
				*s += "<";
			#else
				*s += "(";
			#endif
			list<Node *>::iterator c;
			for(c = children.begin(); c != children.end(); c++) {
				if (c != children.begin())
					*s += ",";
				(*c)->str_c_subtree_hlpr(s);
				if ((*c)->parent() != this)
					cout << "#";
			}
			#ifdef DEBUG_CONTRACTED
				*s += ">";
			#else
				*s += ")";
			#endif
		}
	}

	string str_subtree_twin() {
		string s = "";
		str_subtree_twin_hlpr(&s);
		return s;
	}

	void str_subtree_twin_hlpr(string *s) {
		*s += name.c_str();//str_hlpr(s);
		if (twin != NULL) {
			*s += "{";
				twin->str_subtree_hlpr(s);
			*s += "}";
		}
		if (!is_leaf()) {
			*s += "(";
			list<Node *>::iterator c;
			for(c = children.begin(); c != children.end(); c++) {
				if (c != children.begin())
					*s += ",";
				(*c)->str_subtree_hlpr(s);
				if ((*c)->parent() != this)
					cout << "#";
			}
			*s += ")";
		}
	}


	void print() {
		cout << name;
	}
	void print_subtree() {
		cout << str_subtree();
		cout << endl;
	}

  	void print_subtree_hlpr() {
		cout << str_subtree();
	}
	void print_subtree_twin_hlpr() {
		cout << str_subtree_twin();
	}

	bool is_leaf() {
		return children.empty();
	}

	// TODO: binary only
	bool is_sibling_pair() {
		return (lchild() != NULL && lchild()->is_leaf()
				&& rchild() != NULL && rchild()->is_leaf());
	}

	bool is_temp_parent(){
		return (non_leaf_children == 0 && !is_leaf() 
				&& lchild() != NULL && pre_num == lchild()->get_preorder_number());
	}

	bool is_temp_target_parent(){
		return (non_leaf_children == 0 && !is_leaf() 
				&& pre_num == get_child_min_preorder());
	}
	
	bool is_dead_component() {
		if(this->name == DEAD_COMPONENT)
			return true;
		return false;
	}

	int get_child_min_preorder(){
		int min_pre_num = INT_MAX;
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			if(min_pre_num > (*c)->get_preorder_number()){
				min_pre_num = (*c)->get_preorder_number();
			}
		}
		return min_pre_num;
	}

	void get_all_children_preorder(list<int> *children_preorders) {
	    list<Node *>::iterator c;
		if(is_leaf()){
			children_preorders->push_back(pre_num);
			return;
		}
		for (c = children.begin(); c != children.end(); c++) {
		  	if (!(*c)->is_leaf()) {
				(*c)->get_all_children_preorder(children_preorders);
		  	}
		}		
    }

    bool is_sibling_group() {
	  return non_leaf_children == 0 && !is_leaf();
    }

	bool is_singleton() {
		return (parent() == NULL && is_leaf());
	}

	// TODO: binary only
	void find_sibling_pairs_hlpr(list<Node *> *sibling_pairs) {
		Node *lchild = this->lchild();
		Node *rchild = this->rchild();
		bool lchild_leaf = false;
		bool rchild_leaf = false;
		if (lchild != NULL) {
			if (lchild->is_leaf())
				lchild_leaf = true;
			else
				lchild->find_sibling_pairs_hlpr(sibling_pairs);
		}
		if (rchild != NULL) {
			if (rchild->is_leaf())
				rchild_leaf = true;
			else
				rchild->find_sibling_pairs_hlpr(sibling_pairs);
		}
		if (lchild_leaf && rchild_leaf) {
			sibling_pairs->push_back(lchild);
			sibling_pairs->push_back(rchild);
			//lchild->add_to_sibling_pairs(sibling_pairs, 1);
			//rchild->add_to_sibling_pairs(sibling_pairs, 2);
		}
	}
	
	// find the sibling pairs in this node's subtree
	void append_sibling_pairs(list<Node *> *sibling_pairs) {
		find_sibling_pairs_hlpr(sibling_pairs);
	}

	// find the sibling pairs in this node's subtree
	list<Node *> *find_sibling_pairs() {
		list<Node *> *sibling_pairs = new list<Node *>();
		find_sibling_pairs_hlpr(sibling_pairs);
		return sibling_pairs;
	}
  
	void find_leaves_hlpr(vector<Node *> &leaves) {
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			if ((*c)->is_leaf())
				leaves.push_back(*c);
			else
				(*c)->find_leaves_hlpr(leaves);
		}

	}
      // finds the parents of the sibling groups in this node's subtree
    void append_sibling_groups(list<Node *> *sibling_groups) {
	  find_sibling_groups_hlpr(sibling_groups);	
    }

    // finds the parents of the sibling groups in this node's subtree
    list<Node *> *find_sibling_groups() {
	    list<Node *> *sibling_groups = new list<Node *>();
		find_sibling_groups_hlpr(sibling_groups);
		return sibling_groups;
    }
    void find_sibling_groups_hlpr(list<Node*> *sibling_groups) {
	    list<Node *>::iterator c;
		int total_non_leaves = 0;
		for (c = children.begin(); c != children.end(); c++) {
		  if (!(*c)->is_leaf()) {
			total_non_leaves++;
		    (*c)->find_sibling_groups_hlpr(sibling_groups);
		  }
		}
		non_leaf_children = total_non_leaves;
		if (is_sibling_group())
		  sibling_groups->push_back(this);		
    }

    // Assumes this is the parent of the sibling group in T1
    //outputs a list of list of nodes of siblings in T2 that have the same parents in T2
    void find_identical_sibling_groups(list<list<Node *>> *identical_sibling_groups) {
                identical_sibling_groups->clear();
		list<Node *>::iterator c;
		map<Node *, list<Node*>> parent_to_children = map<Node *, list<Node *>>();
		for (c = children.begin(); c != children.end(); c++) {
		  Node *T1_a = *c;
		  Node *T2_a = T1_a->get_twin();
		  if (T2_a->parent() == NULL)	{
			continue;
		  }
		  auto search = parent_to_children.find(T2_a->parent());
		  //if it does not equal the end, then we already have this parent in the map
		  if (search != parent_to_children.end()) {
			search->second.push_back(T2_a);
		  }
		  //otherwise we add it
		  else {
			parent_to_children.insert({T2_a->parent(), {T2_a}});
		  }
		}
		//Check which parents have more than one identical sibling
		map<Node *, list<Node*>>::iterator parents;
		for (parents = parent_to_children.begin(); parents != parent_to_children.end(); parents++) {
		  if (parents->second.size() > 1) {
			identical_sibling_groups->push_back(parents->second);
		  }
		}
    }
  
	// find the leaves in this node's subtree
	vector<Node *> find_leaves() {
		vector<Node *> leaves = vector<Node *>();
		if (is_leaf())
			leaves.push_back(this);
		else
			find_leaves_hlpr(leaves);
		return leaves;
	}
        void find_nodes_in_subtree_hlpr(vector<Node*> *nodes) {
	  for  (auto i = children.begin(); i != children.end(); i++) {
	    (*i)->find_nodes_in_subtree_hlpr(nodes);
	  }
	  nodes->push_back(this);
	}
        vector<Node*> find_nodes_in_subtree() {
                vector<Node*> nodes = vector<Node*>();
		for (auto i = children.begin(); i != children.end(); i++) {
		  (*i)->find_nodes_in_subtree_hlpr(&nodes);
		}
		return nodes;
	}

      
  //find the preorder numbers of the leaves in this node's subtree
        vector<int> find_leaf_preorders() {
	  vector<Node*> leaves = find_leaves();
	  vector<int> preorders = vector<int>();
	  for (int i = 0; i < leaves.size(); i++) {
	    preorders.push_back(leaves[i]->get_preorder_number());
	  }
	  
	  return preorders;
	}

	void find_interior_hlpr(vector<Node *> &interior) {
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			if (!(*c)->is_leaf()) {
				interior.push_back(*c);
				(*c)->find_interior_hlpr(interior);
			}
		}

	}
  
	// find the interior nodes in this node's subtree
	// does not include this node
	vector<Node *> find_interior() {
		vector<Node *> interior = vector<Node *>();
		if (!is_leaf())
			find_interior_hlpr(interior);
		return interior;
	}
  /*
        int find_arbitrary_lca_hlpr(vector<Node *> &sibling_descendants_counter) {
		list<Node *>::iterator i;
		int total = 0;
		for (i = children->begin(); i != children->end(); i++) {
		  total += (*i)->find_arbitrary_lca_hlpr(sibling_descendants_counter);
		}
		sibling_descendants_counter[get_preorder_number()] += total;
		return sibling_descendants_counter[get_preorder_number()];
	}
  
        // Finds one of the lca's in T2 based on a sibling group parent of T1
        vector<Node *> find_arbitrary_lca(Node* T1_parent) {
	        vector<Node *> sibling_descendants_counter = vector<Node *>();
		list<Node *>::iterator i;
		list<Node *> siblings = T1_parent->get_children();
		for (i = siblings->begin(); i != siblings->end(); i++)
		  {
		    Node* twin = (*i)->get_twin();
		    if (twin != NULL) {
		      //1 it
		    }
		    else {
		      //0 it
		    }
		  }
		for () {
		     //i in components of t2 forest)
		  find_arbitrary_lca_hlpr(sibling_descendants_counter);
		}

	}
  */

        int get_deepest_siblings_hlpr(vector<int> &descendants, Node** to_be_placed) {
	  list<Node*>::iterator i;
	  for (i = children.begin(); i != children.end(); i++) {
	    int descendant_value = descendants[(*i)->get_preorder_number()];
	    if (descendant_value == -1) {
	      *to_be_placed = (*i);
	      return 2;
	    }
	    else if (descendant_value == 1) {
	      return 1 + (*i)->get_deepest_siblings_hlpr(descendants, to_be_placed);
	    }
	  }
	  return -300; //shouldn't get here, it should break when it tries to index -300.
	  //should really assert
	}

        void get_deepest_siblings(vector<int> &descendants, vector<vector<Node*>> &siblings_by_depth) {
	  get_deepest_siblings_strt(descendants, siblings_by_depth);
	}
        vector<vector<Node*>> get_deepest_siblings(vector<int> &descendants) {
	  vector<vector<Node *>> siblings_by_depth = vector<vector<Node *>>(10);
	  get_deepest_siblings_strt(descendants, siblings_by_depth);
	  return siblings_by_depth;
	}
  
        /*
	  Takes a vector of ints with the number of siblings that descend 
	  from that node indexed by their pre-order number
	  Assumes siblings are marked with -1
	  TODO (Ben): proper counting sort
	 */
         void get_deepest_siblings_strt(vector<int> &descendants, vector<vector<Node*>> &siblings_by_depth) {
	   if (get_preorder_number() == -1) { return; } //if -1 then this is rho
	   //This is the sibling, and is a component
	   if (descendants[get_preorder_number()] == -1) {
	     siblings_by_depth[0].push_back(this);
	     return;
	   }
	   list<Node *>::iterator i;
	  //Go down each path of 1s and take note of their path length to this node
	  //There should be the same amount of 1 paths as nodes.size()
	  for (i = children.begin(); i != children.end(); i++) {

	    if (descendants[(*i)->get_preorder_number()] == 1) {
	      Node* placed_node;
	      int depth = (*i)->get_deepest_siblings_hlpr(descendants, &placed_node);
	      while (siblings_by_depth.size() - 1 < depth)
		{
		  siblings_by_depth.resize(siblings_by_depth.size() * 2);
		}
	      if (placed_node != NULL){
		siblings_by_depth[depth].push_back(placed_node);
	      }
	    }	
	    //if it is -1 then this child is the node itself
	    else if (descendants[(*i)->get_preorder_number()] == -1) {
	      siblings_by_depth[1].push_back(*i);
	    }
	  }
	}

        
  
        Node *find_arbitrary_lca_hlpr(vector<int> &descendants) {
	        if (descendants[get_preorder_number()] < 2) { return NULL; }
		list<Node*>::iterator i;
	        for (i = children.begin(); i != children.end(); i++) {
		  if (descendants[(*i)->get_preorder_number()] > 1) {
		    return (*i)->find_arbitrary_lca_hlpr(descendants);
		  }
		}
		return this;
	}
        /*
	  Takes a vector of components in T2, 
	  and a count of how many sibling descendants the node 
	  has in the index of its preorder number
	  Returns lca of the nodes
	  Returns NULL if there is no path between any of the two siblings
	  
	  NOTE: this can be a function
        */
        Node *find_arbitrary_lca(vector<Node *> components, vector<int> &descendants) {	        
		for (int i = 0; i < components.size(); i++) {
		  Node* root = components[i];
		  if (root->get_preorder_number() == -1) {continue;}
		  if (descendants[root->get_preorder_number()] > 1) {
		    return root->find_arbitrary_lca_hlpr(descendants);
		  }
		}
		return NULL;
	}
  /*This is potentially slower than just doing a traversal of the whole tree, 
    since we could be visiting the same path many times. 
    It does mean if this is a short branch and the others are
    long that we dont visit those branches though
    TODO: traverse tree and sum on return. Label sibling group with 1s rest with 0
  */
        void find_pseudo_lca_descendant_count_hlpr(vector <int> &desc_counter) {
		if (p != NULL) {
		  desc_counter[p->get_preorder_number()]++;
		  p->find_pseudo_lca_descendant_count_hlpr(desc_counter);
		}
	}
        vector<int> find_pseudo_lca_descendant_count(int max_preorder) {
	        vector<int> sibling_descendants_counter = vector<int>(max_preorder, 0);
		list<Node *>::iterator i;
		for (i = children.begin(); i != children.end(); i++)
		  {
		    Node* twin = (*i)->get_twin();
		    if (twin != NULL) {
		      //set the twin itself to -1
		      //This will be used later when traversing down to notify we have hit the node
		      sibling_descendants_counter[twin->get_preorder_number()] = -1;
		      twin->find_pseudo_lca_descendant_count_hlpr(sibling_descendants_counter);
		    }
		  }
		return sibling_descendants_counter;
	}
  
	void find_descendants_hlpr(vector<Node *> &descendants) {
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
				descendants.push_back(*c);
			if (!(*c)->is_leaf()) {
				(*c)->find_descendants_hlpr(descendants);
			}
		}

	}
	
	// find the descendants nodes in this node's subtree
	// does not include this node
	vector<Node *> find_descendants() {
		vector<Node *> descendants = vector<Node *>();
		if (!is_leaf())
			find_descendants_hlpr(descendants);
		return descendants;
	}

	bool contains_leaf(int number) {
		if (stomini(get_name()) == number)
			return true;
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			bool ans = (*c)->contains_leaf(number);
			if (ans == true)
				return true;
		}
		vector<Node *> descendants = vector<Node *>();
		return false;
	}

	// make twins point to this tree in this node's subtree
	void resync() {
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->resync();
		}
		if (twin != NULL)
			twin->set_twin(this);
	}

	// remove all twins
	void unsync() {
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->unsync();
		}
		twin = NULL;
	}

	// remove all interior twins
	void unsync_interior() {
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->unsync();
		}
		if(!is_leaf())
			twin = NULL;
	}

	void sync_af_twins() {
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->sync_af_twins();
		}
		if (!is_leaf()) {
			twin = lchild()->get_twin()->parent();
			twin->set_twin(this);
		}
	}

	// find the root of this node's tree
	Node *find_root() {
		Node *root = this;
		while (root->parent() != NULL)
			root = root->parent();
		return root;
	}

	bool same_component(Node *n, int &lca_depth, int &path_length) {
		Node *a = this;
		Node *b = n;
		path_length = 0;
		while(a != b) {
			if ((b == NULL) || (a != NULL && a->get_depth() > b->get_depth()))
				a = a->parent();
			else
				b = b->parent();
			path_length++;
		}
		if (a == NULL) {
			path_length = -1;
			return false;
		}
		lca_depth = a->get_depth();
		return true;
	}

	bool same_component(Node *n) {
		int a,b;
		return same_component(n, a, b);
	}

	bool same_component(Node *n, int &lca_depth) {
		int a;
		return same_component(n, lca_depth, a);
	}

 
	void labels_to_numbers(map<string, int> *label_map, map<int, string> *reverse_label_map) {
		if (name != "") {
			map<string, int>::iterator i = label_map->find(name);
			stringstream ss;
			if (i != label_map->end()) {
				ss << i->second;
				name = ss.str();
			}
			else {
				int num = label_map->size();
				ss << num;
				label_map->insert(make_pair(name, num));
				reverse_label_map->insert(make_pair(num, name));
				name = ss.str();
			}
		}
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->labels_to_numbers(label_map, reverse_label_map);
		}
		if (contracted_lc != NULL)
			contracted_lc->labels_to_numbers(label_map, reverse_label_map);
		if (contracted_rc != NULL)
			contracted_rc->labels_to_numbers(label_map, reverse_label_map);
		for(c = contracted_children.begin(); c != contracted_children.end(); c++) {
			(*c)->labels_to_numbers(label_map, reverse_label_map);
		}
	}


  	void numbers_to_labels(map<int, string> *reverse_label_map) {
		if (name != "") {
			string converted_name = "";
			string current_num = "";
			string::iterator i = name.begin();
			size_t old_loc = 0;
			size_t loc = 0;
			while ((loc = name.find_first_of("0123456789", old_loc)) != string::npos) {
				converted_name.append(name.substr(old_loc, loc - old_loc)); 
				old_loc = loc;
				loc = name.find_first_not_of("0123456789", old_loc);
				string label = "";
				if (loc == string::npos)
					loc = name.size();
				label = name.substr(old_loc, loc - old_loc);
				map<int, string>::iterator j = reverse_label_map->find(atoi(label.c_str()));
				if (j != reverse_label_map->end()) {
					stringstream ss;
					ss << j->second;
					label = ss.str();
				}
				converted_name.append(label);
				old_loc = loc;
			}
			converted_name.append(name.substr(old_loc, name.size() - old_loc)); 
			name = converted_name;



		}
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->numbers_to_labels(reverse_label_map);
		}
		if (contracted_lc != NULL)
			contracted_lc->numbers_to_labels(reverse_label_map);
		if (contracted_rc != NULL)
			contracted_rc->numbers_to_labels(reverse_label_map);
		for(c = contracted_children.begin(); c != contracted_children.end(); c++) {
			(*c)->numbers_to_labels(reverse_label_map);
		}

	}

	void build_name_to_pre_map(map<string, int> *name_to_pre) {
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->build_name_to_pre_map(name_to_pre);
		}
		if (is_leaf())
			name_to_pre->insert(make_pair(get_name(), get_preorder_number()));
	}

	void count_numbered_labels(vector<int> *label_counts) {
		if (name != "") {
			int label = stomini(name);
			if (label_counts->size() <= label)
				label_counts->resize(label+1,0);
			(*label_counts)[label]++;
		}

		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->count_numbered_labels(label_counts);
		}
		if (contracted_lc != NULL)
			contracted_lc->count_numbered_labels(label_counts);
		if (contracted_rc != NULL)
			contracted_rc->count_numbered_labels(label_counts);
	}

  int get_max_preorder_number(int next) {
    if (next < get_preorder_number()) {
      next = get_preorder_number();
    }
    for (auto c = children.begin(); c != children.end(); c++) {
      next = (*c)->get_max_preorder_number(next);
    }
    return next;
  }
  
	void print_preorder_number() {
		cout << name << " - " << pre_num << " - " << edge_pre_start << " - " << edge_pre_end << endl;
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			(*c)->print_preorder_number();
		}
	}

	void preorder_number() {
		preorder_number(0);
	}
  

	int preorder_number(int next) {
		set_preorder_number(next);
		next++;
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			next = (*c)->preorder_number(next);
		}
		return next;
	}

	void edge_preorder_interval() {
		edge_pre_start = pre_num;
		if (is_leaf()) {
			edge_pre_end = pre_num;
		}
		else {
			list<Node *>::iterator c;
			for(c = children.begin(); c != children.end(); c++) {
				(*c)->edge_preorder_interval();
				if (edge_pre_end == -1 || (*c)->edge_pre_end > edge_pre_end)
					edge_pre_end = (*c)->edge_pre_end;
			}
		}
	}

	int size() {
		int s  = 1;
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			s += (*c)->size();
		}
		return s;
	}

	int size_using_prenum() {
		if (is_leaf())
			return get_preorder_number();
		else
			return children.back()->size_using_prenum();
	}

	int max_depth() {
		int d  = 0;
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			int c_d = (*c)->max_depth();
			if (c_d > d)
				d = c_d;
		}
		return d+1;
	}

	int max_degree() {
		int d = children.size();
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			int c_d = (*c)->max_degree();
			if (c_d > d)
				d = c_d;
		}
		return d;
	}

	int get_sackin_index()
	{
		return this->get_sackin_index_helper(0);
	}

	int get_sackin_index_helper(int height)
	{
		if(this->is_leaf()){
			return height;
		}

		int index = 0;
		for(list<Node *>::iterator c = children.begin(); c != children.end(); c++) 
		{
			index += (*c)->get_sackin_index_helper(height+1);
		}
		return index;
	}

	// TODO: binary only
	// these will potentially be removed

void add_to_front_sibling_pairs(list<Node *> *sibling_pairs, int status) {
	sibling_pairs->push_front(this);
	clear_sibling_pair(sibling_pairs);
	sibling_pair_status = status;
	sibling_pair_loc = sibling_pairs->begin();
}

void add_to_sibling_pairs(list<Node *> *sibling_pairs, int status) {
	sibling_pairs->push_back(this);
	clear_sibling_pair(sibling_pairs);
	sibling_pair_status = status;
	sibling_pair_loc = sibling_pairs->end();
	sibling_pair_loc--;
}

void remove_sibling_pair(list<Node *> *sibling_pairs) {
	if (sibling_pair_status > 0) {
		list<Node *>::iterator loc = sibling_pair_loc;
		list<Node *>::iterator sibling_loc = loc;
		if (sibling_pair_status == 1)
			sibling_loc++;
		else if (sibling_pair_status == 2)
			sibling_loc--;

		if (sibling_loc != sibling_pairs->end()) {
			Node *old_sibling = *sibling_loc;
			old_sibling->sibling_pair_status = 0;
			sibling_pairs->erase(sibling_loc);
		}
		sibling_pairs->erase(loc);
		sibling_pair_status = 0;
	}
}

void clear_sibling_pair(list<Node *> *sibling_pairs) {
	if (sibling_pair_status > 0) {
		list<Node *>::iterator loc = sibling_pair_loc;
		list<Node *>::iterator sibling_loc = loc;
		if (sibling_pair_status == 1)
			sibling_loc++;
		else if (sibling_pair_status == 2)
			sibling_loc--;

		if (sibling_loc != sibling_pairs->end()) {
			Node *old_sibling = *sibling_loc;
			old_sibling->sibling_pair_status = 0;
		}
		sibling_pair_status = 0;
	}
}

Node *get_sibling() {
	Node *ret = NULL;
	if (p != NULL && p->children.size() > 1) {
		list<Node *>::iterator s = p_link;
		if (p_link == p->children.begin())
			s++;
		else
			s--;
		ret = *s;
	}
	return ret;
}

Node *get_right_sibling() {
	Node *ret = NULL;
	if (p != NULL && p->children.size() > 1) {
		list<Node *>::iterator s = p_link;
		s++;
		if (s != children.end())
			ret = *s;
	}
	return ret;
}

// TODO: binary only
Node *get_sibling(list<Node *> *sibling_pairs) {
	if (sibling_pair_status > 0) {
		list<Node *>::iterator loc = sibling_pair_loc;
		list<Node *>::iterator sibling_loc = loc;
		if (sibling_pair_status == 1)
			sibling_loc++;
		else if (sibling_pair_status == 2)
			sibling_loc--;
		return *sibling_loc;
	}
	else
		return NULL;
}

void set_sibling(Node *sibling) {
	if (sibling->sibling_pair_status > 0) {
		sibling_pair_loc = sibling->sibling_pair_loc;
		if (sibling->sibling_pair_status == 1)
			sibling_pair_loc++;
		else if (sibling->sibling_pair_status == 2)
			sibling_pair_loc--;
	}
}

void clear_sibling_pair_status() {
	sibling_pair_status = 0;
}

// fix parents
void fix_parents() {
		list<Node *>::iterator c;
		for(c = children.begin(); c != children.end(); c++) {
			if ((*c)->parent() != this) {
				(*c)->p = this;
				(*c)->p_link = c;
			}
			(*c)->fix_parents();
		}
}

// TODO: binary only
/* TODO: think of a better way of trying all rootings. Maybe create
	 a function that copies a tree and then reroots it with a given
	 subtree as a child of the root.
*/

void left_rotate() {
	if (lchild()->lchild() != NULL) {
		Node *new_lc = lchild()->lchild();
		Node *new_rc = lchild();
		Node *new_rc_lc = lchild()->rchild();
		Node *new_rc_rc = rchild();

		new_lc->cut_parent();
		new_rc->cut_parent();
		new_rc_lc->cut_parent();
		new_rc_rc->cut_parent();
		add_child(new_lc);
		add_child(new_rc);
		rchild()->add_child(new_rc_lc);
		rchild()->add_child(new_rc_rc);
	}
}

void right_rotate() {
	if (rchild()->lchild() != NULL) {
		Node *new_lc = rchild()->lchild();
		Node *new_rc = rchild();
		Node *new_rc_lc = rchild()->rchild();
		Node *new_rc_rc = lchild();

		new_lc->cut_parent();
		new_rc->cut_parent();
		new_rc_lc->cut_parent();
		new_rc_rc->cut_parent();
		add_child(new_lc);
		add_child(new_rc);
		rchild()->add_child(new_rc_lc);
		rchild()->add_child(new_rc_rc);
	}
}

void next_rooting() {
	if (lchild()->lchild() != NULL)
		left_rotate();
	else if (rchild()->lchild() != NULL)
			right_rotate();
	else
		return;
	if (lchild()->pre_num < rchild()->pre_num && (lchild()->pre_num != 1 || rchild()->pre_num != 2 || !lchild()->is_leaf()))
		next_rooting();
}

/* reroot
 * rerooots the tree between new_lc and its parent. Maintains this
 * Node as the root
 * assumes that new_lc is a descendant of this Node
 * assumes this node is the root and has two children
 * other nodes may be multifurcating
*/
void reroot(Node *new_lc) {
	Node *new_rc = new_lc->parent();
	if (new_rc == this || new_rc == NULL) {
		return;
	}
	Node *prev = new_rc;
	Node *next = new_rc->parent();
//	Node *old_lc = lchild();
//	Node *old_rc_rc = rchild();
	new_lc->cut_parent();
	new_rc->cut_parent();
	while(next != NULL) {
		Node *current = next;
		next = current->parent();
		prev->add_child(current);
		prev = current;
	}
	Node *root = prev;
	while(root->get_children().size() > 0)
		root->parent()->add_child(root->lchild());
//	root->parent()->add_child(root->lchild());
	root->cut_parent();
	root->add_child(new_lc);
	root->add_child(new_rc);
}

void reroot_clean(Node *new_lc) {
	this->reroot(new_lc);
	this->set_depth(0);
	this->fix_depths();
	this->preorder_number();
	this->edge_preorder_interval();
}

// make the root binay again
void fixroot() {
	if (get_children().size() > 2) {
		Node *new_lc = new Node();
		while(get_children().size() > 1) {
			new_lc->add_child(get_children().back());
		}
		add_child(new_lc);
	}
}

  //Currently only used on T2, if you want to use it on T1, update non_leaf_children
Node *expand_children_out(list<Node*> nodes) {
  Node *new_child = new Node();
  //new_child->set_parent(this);
  int min_preorder = 0x7FFFFFFF;
  list<Node *>::iterator i;
  for (i = nodes.begin(); i != nodes.end(); i++) {
    new_child->add_child(*i);
    (*i)->set_parent(new_child);//not needed?
    children.remove(*i);
    if ((*i)->get_preorder_number() < min_preorder) {
      min_preorder = (*i)->get_preorder_number();
    }
  }
  //super hacky, may overlap with others. Currently this is being used only when the other node is being contracted anyways, but be careful
  new_child->set_preorder_number(min_preorder);
  add_child(new_child);
  return new_child;
}

Node *expand_parent_edge(Node *n) {
	if (p != NULL) {
		Node *old_p = p;
		cut_parent();
		Node *new_p = new Node();
		new_p->add_child(n);
		old_p->add_child(new_p);
		return p;
	}
	else {
		Node *new_child = new Node(name);
		Node *old_lc = lchild();
		Node *old_rc = rchild();
		old_lc->cut_parent();
		old_rc->cut_parent();

		new_child->add_child(old_lc);
		new_child->add_child(old_rc);
		new_child->contracted_lc = contracted_lc;
		if (contracted_lc != NULL)
			contracted_lc->p = new_child;
		new_child->contracted_rc = contracted_rc;
		if (contracted_rc != NULL)
		contracted_rc->p = new_child;
		name = "";
		contracted_lc = NULL;
		contracted_rc = NULL;
		add_child(new_child);
		return this;
	}
}

Node *undo_expand_parent_edge() {
	if (p != NULL) {
		Node *old_p = p;
		cut_parent();
		Node *child = children.front();
		child->cut_parent();
		old_p->add_child(child);
		return this;
	}
	else {
		Node *child = children.front();
		name = child->name;
		Node *new_lc = child->lchild();
		Node *new_rc = child->rchild();
		new_lc->cut_parent();
		new_rc->cut_parent();
		contracted_lc = child->contracted_lc;
		if (contracted_lc != NULL)
			contracted_lc->p = this;
		contracted_rc = child->contracted_rc;
		if (contracted_rc != NULL)
			contracted_rc->p = this;
		child->contracted_lc = NULL;
		child->contracted_rc = NULL;
		child->name = "";
		add_child(new_lc);
		add_child(new_rc);
		return child;
	}
}

/*  apply an SPR operation to move this subtree to be a
 *	sibling of new_sibling
 *
 *  Note: problems will occur if new_sibling is a descendant of this
 *  Note: does nothing if this is the root
 *
 *	Returns a node that will reverse the spr (old_sibling unless it
 *		was moved to maintain a root, in which case the root is returned,
 *		NULL if no spr occured)
 *
 * Note: binary only
 * TODO (ben): multifurcating
 */
Node *spr(Node *new_sibling, int &which_child) {
	Node *reverse;
	int prev_child_loc = 0;
	if (p == NULL || new_sibling == NULL)
		return NULL;
	Node *old_sibling = get_sibling();
	if (old_sibling == new_sibling)
		return NULL;
	Node *grandparent = p->p;
	if (p->lchild() == this)
		prev_child_loc = 1;
	else
		prev_child_loc = 2;
	// Prune
	if (grandparent != NULL) {
		bool leftc = false;
		if (old_sibling->parent() == grandparent->lchild())
			leftc = true;

		old_sibling->cut_parent();
		//p->delete_child(old_sibling);
		p->cut_parent(); //LEAK?
//		grandparent->delete_child(p);
		if (leftc && !grandparent->is_leaf()) {
				Node *ns = grandparent->children.front();
				grandparent->insert_child(ns,old_sibling);
		}
		else
			grandparent->add_child(old_sibling);

		reverse = old_sibling;
	}
	else {
		if (old_sibling->is_leaf())
			return NULL;
		Node *root = p;
		bool leftc = false;
		if (root->lchild() == this)
			leftc = true;
		Node *lc = old_sibling->lchild();
		Node *rc = old_sibling->rchild();
		root->delete_child(this);
		root->delete_child(old_sibling) ;
		if (lc != NULL)
			root->add_child(lc);
		if (rc != NULL)
			root->add_child(rc);

		//old_sibling->delete_child(old_sibling->lchild());
		//old_sibling->delete_child(old_sibling->rchild());
		if (leftc) {
			if (old_sibling->is_leaf())
				old_sibling->add_child(this);
			else {
				Node *new_sibling = old_sibling->children.front();
				old_sibling->insert_child(new_sibling,this);
			}
		}
		else
			old_sibling->add_child(this);
		reverse = root;
	}


	// Regraft
	if (new_sibling->p != NULL) {
		grandparent = new_sibling->p;
//		grandparent->delete_child(new_sibling);
		bool leftc = false;
		if (new_sibling == grandparent->lchild())
			leftc = true;
		if (leftc && !grandparent->is_leaf()) {
				Node *new_sibling = grandparent->children.front();
				grandparent->insert_child(new_sibling,p);
		}
		else {
			grandparent->add_child(p);
		}

		if (which_child == 1)
			p->add_child(new_sibling);
		else {
				Node *ns = p->children.front();
				p->insert_child(ns,new_sibling);
		}

	}
	else {
		Node *root = new_sibling;
		new_sibling = p;
		// TODO: still broken
		Node *lc = root->lchild();
		Node *rc = root->rchild();
		p->add_child(lc);
		p->add_child(rc);
//		p->delete_child(this);
		// problem here
		if (which_child == 0)
			which_child = prev_child_loc;
		if (which_child == 1) {
			root->add_child(this);
			root->add_child(new_sibling);
		}
		else {
			root->add_child(new_sibling);
			root->add_child(this);
		}


	}

	which_child = prev_child_loc;
	return reverse;
}

Node *spr(Node *new_sibling) {
	int na = 0;
	return spr(new_sibling, na);
}

  /*
    prune this subtree and attach as a sibling of new_sibiling
   */

void spr_mult(Node *new_sibling) {
  //prune
  Node* parent = p;
  if (parent != NULL) {
    cut_parent();
    if (parent->get_children().size() == 1) {
      parent->contract(true);
    }
  }
  else {
    return;
  }
  //regraft
  if (new_sibling->parent() != NULL) {
    Node* new_parent = new Node();
    new_sibling->parent()->add_child(new_parent);
    //new_parent->set_parent(new_sibling->parent());
    new_parent->add_child(new_sibling);
    new_parent->add_child(this);
    }
  else {
    Node* new_parent = new_sibling;//new Node();
    Node* replace = new Node(*new_sibling);
    new_parent->get_children().clear();
    new_parent->add_child(replace);
    new_parent->add_child(this);
  }
}

  
void find_descendant_counts_hlpr(vector<int> *dc) {
	int num_descendants = 0;
	list<Node *>::iterator c;
	for(c = children.begin(); c != children.end(); c++) {
		(*c)->find_descendant_counts_hlpr(dc);
		num_descendants += (*dc)[(*c)->get_preorder_number()];
		num_descendants += 1;
	}
	if (dc->size() <= get_preorder_number())
		dc->resize(get_preorder_number() + 1, -1);
	(*dc)[get_preorder_number()] = num_descendants;
}

vector<int> *find_descendant_counts() {
	vector<int> *dc = new vector<int>();
	find_descendant_counts_hlpr(dc);
	return dc;
}

void find_leaf_counts_hlpr(vector<int> *lc) {
	int num_leaves = 0;
	list<Node *>::iterator c;
	for(c = children.begin(); c != children.end(); c++) {
		(*c)->find_leaf_counts_hlpr(lc);
		num_leaves += (*lc)[(*c)->get_preorder_number()];
	}
	if (is_leaf())
		num_leaves++;
	if (lc->size() <= get_preorder_number())
		lc->resize(get_preorder_number() + 1, -1);
	(*lc)[get_preorder_number()] = num_leaves;
}

vector<int> *find_leaf_counts() {
	vector<int> *lc = new vector<int>();
	find_leaf_counts_hlpr(lc);
	return lc;
}

Node *find_median() {
	vector<int> *dc = find_descendant_counts();
	return find_median_hlpr(dc, (*dc)[get_preorder_number()] / 2);
}

Node *find_subtree_of_size(double percentage) {
	vector<int> *dc = find_descendant_counts();
	return find_median_hlpr(dc, (int)((*dc)[get_preorder_number()] * percentage));
}

Node *find_median_hlpr(vector<int> *dc, int target_size) {
	Node *largest_child_subtree = NULL;
	int lcs_size = 0;
	list<Node *>::iterator c;
	for(c = children.begin(); c != children.end(); c++) {
		int cs_size = (*dc)[(*c)->get_preorder_number()] + 1; 
		if (cs_size > lcs_size) {
			largest_child_subtree = *c;
			lcs_size = cs_size;
		}
	}
	if (lcs_size > target_size)
		return largest_child_subtree->find_median_hlpr(dc, target_size);
	else
		return this;
}

int any_leaf_preorder_number() {
	if (is_leaf()) {
		if (contracted_lc != NULL)
			return contracted_lc->any_leaf_preorder_number();
		else return get_preorder_number();
	}
	else
		return (*children.begin())->any_leaf_preorder_number();
}

Node *find_by_prenum(int prenum) {
//	cout << "find_by_prenum: " << str_subtree() << endl;
	if (prenum == get_preorder_number())
		return this;
	Node *search_child = NULL;
	int best_prenum = -1;
	list<Node *>::iterator c;
	for(c = children.begin(); c != children.end(); c++) {
//		cout << "\tchild: " << (*c)->str_subtree() << endl;
//		cout << "\tpre: " << (*c)->get_preorder_number() << endl;
		int p = (*c)->get_preorder_number();
		if (p == prenum)
			return *c;
		else if (p < prenum && p > best_prenum) {
			best_prenum = p;
			search_child = *c;
		}
	}
	if (search_child != NULL)
		return search_child->find_by_prenum(prenum);
	else
		return NULL;
}

Node *find_by_prenum_full(int prenum) {
	if (prenum == get_preorder_number())
		return this;
	list<Node *>::iterator c;
	for(c = children.begin(); c != children.end(); c++) {
		Node *test = (*c)->find_by_prenum_full(prenum);
		if (test != NULL)
			return test;
	}
	return NULL;
}

// expand all contracted nodes of a subtree starting at n
void expand_contracted_nodes() {
  	list<Node *>::iterator c;
	if (is_leaf()) {
		if (contracted_lc != NULL) {
			add_child(contracted_lc);
			contracted_lc = NULL;
		}
		if (contracted_rc != NULL) {
			add_child(contracted_rc);
			contracted_rc = NULL;
		}
		for(c = contracted_children.begin(); c != contracted_children.end(); c++) {
		        add_child(*c);
		}
		contracted_children.clear();

	}
	for(c = get_children().begin(); c != get_children().end(); c++) {
		(*c)->expand_contracted_nodes();
	}
	/*
	//if -1 then this is a node we created for a contraction, and it should go up one
	if (get_edge_pre_end() == -1) {
	  if (p != NULL) {
	    while (!children.empty()) {
	      p->add_child(children.front());
	    }
	    delete this;
	    return;
	  }
	  }*/
}

int get_name_num() {
	return atoi(get_name().c_str());
}


};

// function prototypes
Node *build_tree(string s);
Node *build_tree(string s, set<string, StringCompare> *include_only);
Node *build_tree(string s, int start_depth);
Node *build_tree(string s, int start_depth, set<string, StringCompare> *include_only);

int build_tree_helper(int start, const string& s, Node *parent,
		bool &valid, set<string, StringCompare> *include_only);
//void preorder_number(Node *node);
//int preorder_number(Node *node, int next);
string strip_newick_name(string &T);
vector<Node *> contract_deepest_siblings(vector<vector<Node *>> &siblings_by_depth);


// build a tree from a newick string
Node *build_tree(string s) {
	return build_tree(s, 0, NULL);
}
Node *build_tree(string s, int start_depth) {
	return build_tree(s, start_depth, NULL);
}
Node *build_tree(string s, set<string, StringCompare> *include_only) {
	return build_tree(s, 0, include_only);
}
Node *build_tree(string s, int start_depth, set<string, StringCompare> *include_only) {
	if (s == "")
		return new Node();
	Node *dummy_head = new Node("p", start_depth-1);
	bool valid = true;
	build_tree_helper(0, s, dummy_head, valid, include_only);
	Node *head = dummy_head->lchild();
	if (valid && head != NULL) {
		delete dummy_head;
		return head;
	}
	else {
		if (head != NULL)
			head->delete_tree();
		return dummy_head;
	}

}

// build_tree recursive helper function
int build_tree_helper(int start, const string& s, Node *parent,
		bool &valid, set<string, StringCompare> *include_only) {
	int loc = s.find_first_of("(,)", start);
	if (loc == string::npos) {
		string name = s.substr(start, s.size() - start);
		int name_end = name.find(':');
		if (name_end != string::npos)
			name = name.substr(0, name_end);
		if (include_only == NULL ||
				include_only->find(name) != include_only->end()) {
			Node *node = new Node(name);
			parent->add_child(node);
		}
		loc = s.size()-1;
		return loc;
	}
	while(s[start] == ' ' || s[start] == '\t')
		start++;
	int end = loc;
	while(s[end] == ' ' || s[end] == '\t')
		end--;
	string name = s.substr(start, end - start);
	int name_end = name.find(':');
	if (name_end != string::npos)
		name = name.substr(0, name_end);
	Node *node = NULL;
	if (include_only == NULL ||
			include_only->find(name) != include_only->end()) {
		node = new Node(name);
		parent->add_child(node);
	}

	int count = 1;
	if (s[loc] == '(') {
			loc = build_tree_helper(loc + 1, s, node, valid, include_only);
			while(s[loc] == ',') {
				loc = build_tree_helper(loc + 1, s, node, valid, include_only);
				count++;
			}
//			int loc_check = s.find_first_of("(,)", loc);
//			if (loc_check != string::npos &&
//					s[loc_check] == ','
			if (s[loc] != ')'
					|| (IGNORE_MULTI && count > 2)) {
				valid = false;
					return s.size()-1;
			}
			// TODO: get the support values here (and branch lengths?)
			// contract_node() if support is less than a threshold
			loc++;
			if (s[loc-1] == ')') {
				int numc = node->get_children().size();
				bool contracted = false;
				int next = s.find_first_of(",)", loc);
				if (next != string::npos) {
					if (next > loc && (REQUIRED_SUPPORT > 0 || MIN_LENGTH >= 0)) {
						string info = s.substr(loc, next - loc);
						// support values (a number after a close bracket)
						if (info[0] != ':') {
							double support = atof(info.c_str()) / SUPPORT_BY_VALUE;
							#if DEBUG_SUPPORT
							cout << "support=" << support << endl;
							#endif
							if (REQUIRED_SUPPORT > 0 && support < REQUIRED_SUPPORT && numc > 0) {
							  #if DEBUG_SUPPORT
 								cout << "contracting (support)" << endl;
								#endif
								node->contract_node();
								contracted = true;
							}
						}
						// branch lengths (a number after a colon)
						int next_colon = info.find_first_of(":");
						if (next_colon != string::npos) {
							string length_str = info.substr(next_colon+1);
							double length = atof(length_str.c_str());
#if DEBUG_SUPPORT
							cout << "length=" << length << endl;
							#endif
							if (MIN_LENGTH >= 0 && !contracted && length <= MIN_LENGTH && numc > 0) {
#if DEBUG_SUPPORT
								cout << "contracting (length)" << endl;
								#endif
								node->contract_node();
								contracted = true;
							}
						}
						// support values (a number in square brackets)
						int square_bracket = info.find_first_of("[");
						if (square_bracket != string::npos) {
							string support_str = info.substr(square_bracket+1);
							double support = atof(support_str.c_str()) / 100.0;
							#if DEBUG_SUPPORT
							cout << "support=" << support << endl;
							#endif
							if (REQUIRED_SUPPORT > 0 && !contracted && support < REQUIRED_SUPPORT && numc > 0) {
							  #if DEBUG_SUPPORT
 								cout << "contracting (support)" << endl;
								#endif
								node->contract_node();
								contracted = true;
							}
						}

					}
					loc=next;
				}
				if (!contracted) {
					if (numc == 1)
						node->contract_node();
					else if (numc == 0 && name == "") {
						node->cut_parent();
						delete node;
					}
				}
			}
	}
	return loc;
}

// swap two nodes
void swap(Node **a, Node **b) {
	Node *temp = *a;
	*a = *b;
	*b = temp;
}


/*
void preorder_number(Node *node) {
	preorder_number(node, 0);
}
int preorder_number(Node *node, int next) {
	node->set_preorder_number(next);
	next++;
	if(node->lchild() != NULL) {
		next = preorder_number(node->lchild(), next);
	}
	if(node->rchild() != NULL) {
		next = preorder_number(node->rchild(), next);
	}
	return next;
}
*/
// return the smallest number in s

int stomini(string s) {
//	cout << "stomini" << endl;
	string number_characters = "+-0123456789";
	int i = 0;
	int min = INT_MAX;
	string current = "";
	for(int i = 0; i < s.size(); i++) {
		if (number_characters.find(s[i]) != string::npos) {
			current += s[i];
		}
		else if (current.size() > 0) {
			int num = atoi(current.c_str());
			if (num < min)
				min = num;
			current = "";
		}
	}
	if (current.size() > 0) {
		int num = atoi(current.c_str());
		if (num < min)
			min = num;
		current = "";
	}
//	cout << "returning " << min << endl;
	return min;
}

// assumes that an unrooted tree is represented with a 3-way multifurcation
string root(string s) {
//	cout << "root(string s)" << endl;
//	cout << s << endl;
	string r = "";
	int i = 0;
	int depth = 0;
	int first_c = -1;
	int second_c = -1;
	int last_bracket = -1;
	for(int i = 0; i < s.size(); i++) {
		if (s[i] == '(')
			depth++;
		else if (s[i] == ')') {
			depth--;
			last_bracket = i;
		}
		else if (depth == 1 && s[i] == ',') {
			if (first_c == -1)
				first_c = i;
			else if  (second_c == -1)
				second_c = i;
		}
	}
	if (second_c == -1 || last_bracket == -1)
		return s;
	else {
		r.append(s.substr(0,first_c+1));
		r.append("(");
		r.append(s.substr(first_c+1,last_bracket-first_c));
		r.append(")");
		r.append(s.substr(last_bracket+1,string::npos));
	}
//	cout << r << endl;
	return r;
}

template <typename T> vector<T> &random_select(vector <T> &V, int n) {
	vector<T> *ret = new vector<T>;
	int end = V.size();
	for(int i = 0; i < n && end > 0; i++) {
		int x = rand() % end;
		ret->push_back(V[x]);
		V[x] = V[end-1];
		// Is it better to remove them or leave them?
		V.pop_back();
		end--;
	}
	return *ret;
}

template <typename T> void print_vector(vector <T> V) {
	int end = V.size();
	for(int i = 0; i < end; i++) {
		if (i > 0)
			cout << ",";
		cout << V[i];
	}
	cout << endl;
}

string strip_newick_name(string &line) {
	string name;
	size_t loc = line.find_first_of("(");
	if (loc != string::npos) {
		name = "";
		if (loc != 0) {
			name = line.substr(0,loc);
			line.erase(0,loc);
		}
	}
	return name;
}

void reroot_safe(Node **n, Node *new_lc) {
	(*n)->reroot(new_lc);
	(*n)->set_depth(0);
	(*n)->fix_depths();
	(*n)->preorder_number();
	(*n)->edge_preorder_interval();
	// hack TODO fix
	string str = (*n)->str_subtree();
	(*n)->delete_tree();
	(*n) = build_tree(str);
	(*n)->preorder_number();
	(*n)->edge_preorder_interval();
}

//TODO: remove references to this, replace with map version
vector<Node *> contract_deepest_siblings(vector<vector<Node *>> &siblings_by_depth) {	  
	  vector<Node*> deepest_siblings = vector<Node*>();
	  for (int j = siblings_by_depth.size() - 1; j >= 0; j--) {
	    if (siblings_by_depth[j].size() != 0) {
	      for (int k = 0; k < siblings_by_depth[j].size(); k++) {
		deepest_siblings.push_back(siblings_by_depth[j][k]);
	      }
	    }
	  }
	  return deepest_siblings; 
	}
vector<Node *> contract_deepest_siblings(vector<vector<Node *>> &siblings_by_depth, map<Node*, int> *s_map) {	  
	  vector<Node*> deepest_siblings = vector<Node*>();
	  for (int j = siblings_by_depth.size() - 1; j >= 0; j--) {
	    if (siblings_by_depth[j].size() != 0) {
	      for (int k = 0; k < siblings_by_depth[j].size(); k++) {
		deepest_siblings.push_back(siblings_by_depth[j][k]);
		s_map->insert({siblings_by_depth[j][k], j - 1});
	      }
	    }
	  }
	  
	  return deepest_siblings; 
	}





#endif
