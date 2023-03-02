/*******************************************************************************
lgt.h

Utilities for inferring LGT events.

Copyright 2013-2014 Chris Whidden
whidden@cs.dal.ca
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

#ifndef LGT
#define LGT

//#define DEBUG_LGT

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
#include <algorithm>
#include "Forest.h"
#include "ClusterForest.h"
#include "LCA.h"
#include "ClusterInstance.h"
#include "SiblingPair.h"
#include "UndoMachine.h"

bool LGT_MOVE_PARENT = false;
bool LGT_MOVE_INDIVIDUAL_NODE = false;
bool LGT_MAINTAIN_LIST = false;

using namespace std;

class transfer {
	public:
	int source_pre;
	int target_pre;
	int final_source_depth;
        vector<int> source_child_pre_nums;
        vector<int> target_child_pre_nums;

	transfer(int s, int t, int f_d) {
		source_pre = s;
		target_pre = t;
		final_source_depth = f_d;
	}
  transfer(int s, int t, int f_d, vector<int> scpn, vector<int> tcpn) {
		source_pre = s;
		target_pre = t;
		final_source_depth = f_d;
		source_child_pre_nums = scpn;
		target_child_pre_nums = tcpn;
	}

};

void add_transfers(vector<vector<int> > *transfer_counts, Node *super_tree,
		vector<Node *> *gene_trees);
void add_transfers(vector<vector<int> > *transfer_counts, Forest *F1,
		Forest *F2, Forest *MAF1, Forest *MAF2);
void add_transfers(vector<vector<int> > *transfer_counts, Forest *F1,
		Forest *F2, Forest *MAF1, Forest *MAF2, map<int, string> *reverse_label_map);
void add_transfers(vector<vector<int> > *transfer_counts, Node *super_tree,
		vector<Node *> *gene_trees, map<int, string> *reverse_label_map);
void print_transfers(Node *super_tree, Forest *F1, Forest *F2, Forest *MAF1,
		Forest *MAF2, map<int, string> *reverse_label_map);
void print_leaf_list(Node *F1_source, map<int, string> *reverse_label_map);
bool map_transfer(Node *F2_source, Forest *F1, Forest *MAF2,
		Node **F1_source_out, Node **F1_target_out);
Node *find_best_target(Node *source, Forest *AF);
Node *find_best_target(Node *source, Node *target, Node **best_target);
void show_moves(Node *T1, Node *T2, map<string, int> *label_map,
		map<int, string> *reverse_label_map);


void add_transfers(vector<vector<int> > *transfer_counts, Node *super_tree,
		vector<Node *> *gene_trees, map<int, string> *reverse_label_map) {
	#pragma omp parallel for
	for(int i = 0; i < gene_trees->size(); i++) {
		Forest *MAF1 = NULL;
		Forest *MAF2 = NULL;
		Forest F1 = Forest(super_tree);
		Forest F2 = Forest((*gene_trees)[i]);
		if (sync_twins(&F1,&F2)) {
			int distance = rSPR_branch_and_bound_simple_clustering(F1.get_component(0), F2.get_component(0), &MAF1, &MAF2);
			expand_contracted_nodes(MAF1);
			expand_contracted_nodes(MAF2);
#ifdef DEBUG_LGT			
			cout << i << ": " << distance << endl;
			//cout << "\tT1: "; F1.print_components();
			//cout << "\tT2: "; F2.print_components();
			cout << "\tT1: "; F1.print_components_with_edge_pre_interval();
			cout << "\tT2: "; F2.print_components_with_edge_pre_interval();
			cout << "\tF1: "; MAF1->print_components_with_edge_pre_interval();
			cout << "\tF2: "; MAF2->print_components_with_edge_pre_interval();
#endif
			sync_af_twins(MAF1, MAF2);
			add_transfers(transfer_counts, &F1, &F2, MAF1, MAF2, reverse_label_map);
		}
		if (MAF1 != NULL)
			delete MAF1;
		if (MAF2 != NULL)
			delete MAF2;
	}
}

void print_transfer(Node *F1_source, Node *F1_target, map<int, string> *reverse_label_map){
	cout << F1_source->get_name() << "(" << F1_source->get_preorder_number() << ")" << " -> " 
		 << F1_target->get_name() << "(" << F1_target->get_preorder_number() << ")" << endl;
	cout << "F1 source " << endl;
	F1_source->print_subtree();
	print_leaf_list(F1_source, reverse_label_map);
	cout << "F1 target " << endl;
	F1_target->print_subtree();
	print_leaf_list(F1_target, reverse_label_map);
}

void add_transfers(vector<vector<int> > *transfer_counts, Forest *F1,
		Forest *F2, Forest *MAF1, Forest *MAF2, map<int, string> *reverse_label_map) {
	int start = 1;
	if (MAF2->contains_rho())
		start = 0;

#ifdef DEBUG
	cout << "MAF1 : ";
	MAF1->print_components();
	MAF1->print_components_preorder();
	cout << "MAF2 : ";
	MAF2->print_components();
	MAF2->print_components_preorder();
#endif

	cout << "LGT events" << endl;
	int transfer_count = 0;
	map<string, list<list<int>>> map_transfer_sibling;			

	for(int i = start; i < MAF2->num_components(); i++) {
		Node *F2_source = MAF2->get_component(i);
		Node *F1_source;
		Node *F1_target;
		if (!map_transfer(F2_source, F1, MAF2, &F1_source,
					&F1_target)) {
			continue;
		}

		#pragma omp atomic
		if(!IGNORE_MULTI && MULTIFURCATING){
			Node *F1_source_new = F1_source;
			Node *F1_target_new = F1_target;
			if(LGT_MOVE_PARENT){
				if(F1_target->is_temp_target_parent()){
					Node* f1_child = F1->find_by_prenum(
									F1_target->get_preorder_number());
					if(f1_child != NULL && f1_child->parent() != NULL){
						F1_target_new = f1_child->parent();
					}
				}
				if(F1_source->is_temp_parent()){
					Node* f1_lchild = F1->find_by_prenum(
										F1_source->get_preorder_number());
					if(f1_lchild != NULL && f1_lchild->parent() != NULL){
						F1_source_new = f1_lchild->parent();
					}
				}
				(*transfer_counts)[F1_source_new->get_preorder_number()][F1_target_new->get_preorder_number()]++;
				print_transfer(F1_source_new, F1_target_new, reverse_label_map);
				transfer_count++;	
			}
			else if(LGT_MOVE_INDIVIDUAL_NODE){
				list<Node*> lstSource = {F1_source_new};
				list<Node*> lstTarget = {F1_target_new};
				if(F1_target->is_temp_target_parent()){
					/*Node* f1_child = F1->find_by_prenum(
									F1_target->get_preorder_number());

					if(f1_child != NULL && f1_child->parent() != NULL){
						F1_target_new = f1_child->parent();
						lstTarget = F1_target_new->get_children();
					}*/
					lstTarget = F1_target->get_children();
				}
				if(F1_source->is_temp_parent()){
					lstSource = F1_source->get_children();
				}
				list<Node *>::const_iterator c, ct;
				for(c = lstSource.begin(); c != lstSource.end(); c++) {
					for(ct = lstTarget.begin(); ct != lstTarget.end(); ct++) {
						(*transfer_counts)[(*c)->get_preorder_number()][(*ct)->get_preorder_number()]++;
						print_transfer(*c, *ct, reverse_label_map);
						transfer_count++;
					}
				}
			}
			else if(LGT_MAINTAIN_LIST){
				list<int> F1_source_sibling;
				list<int> F1_target_sibling;
				if(F1_target->is_temp_parent()){
					Node* f1_child = F1->find_by_prenum(
										F1_target->get_preorder_number());
					if(f1_child != NULL && f1_child->parent() != NULL){
						F1_target_new = f1_child->parent();
					}
					list<Node *>::const_iterator c;
					for(c = F1_target->get_children().begin(); 
								c != F1_target->get_children().end(); c++) {
						F1_target_sibling.push_back((*c)->get_preorder_number());
						F1_target->get_all_children_preorder(&F1_target_sibling);
					}
				}
				if(F1_source->is_temp_parent()){
					Node* f1_lchild = F1->find_by_prenum(
										F1_source->get_preorder_number());
					if(f1_lchild != NULL && f1_lchild->parent() != NULL){
						F1_source_new = f1_lchild->parent();
					}
					list<Node *>::const_iterator c;
					for(c = F1_source->get_children().begin(); 
								c != F1_source->get_children().end(); c++) {
						F1_source_sibling.push_back((*c)->get_preorder_number());
						F1_source->get_all_children_preorder(&F1_source_sibling);
					}
				}
				string key = to_string(F1_source_new->get_preorder_number()) 
								+ "_" + to_string(F1_target_new->get_preorder_number());
				map_transfer_sibling[key] = {F1_source_sibling, F1_target_sibling};
				(*transfer_counts)[F1_source_new->get_preorder_number()][F1_target_new->get_preorder_number()]++;
				print_transfer(F1_source_new, F1_target_new, reverse_label_map);
				transfer_count++;			
			}
			else{
				(*transfer_counts)[F1_source->get_preorder_number()][F1_target->get_preorder_number()]++;
				print_transfer(F1_source, F1_target, reverse_label_map);
				transfer_count++;
			}
		}
		else{
			(*transfer_counts)[F1_source->get_preorder_number()][F1_target->get_preorder_number()]++;
			print_transfer(F1_source, F1_target, reverse_label_map);
			transfer_count++;
		}

		// do we want to check that the move is valid here?
	}
	cout << "Transfer count " << transfer_count << endl;
	//Print transfer map
	if(LGT_MAINTAIN_LIST){
		cout << "Transfer map" << endl;
		for(auto it = map_transfer_sibling.cbegin(); it != map_transfer_sibling.cend(); ++it){
			std::cout << "Key : " << it->first << endl;
			for(auto &lst : it->second){
				for (auto const &i: lst) {
					cout << i << " ";
				}
				cout << endl;
			}
		}
	}
	cout << endl;
					// TODO: identify a valid move (if any) for each component
					// loop over the components of MAF2
							// map source to destination in MAF2
							// map source to MAF1
							// map dest to MAF1
							// check if move is valid
							// ? map to full trees? automatic with the pre numbers?
							// store in a list/vector? counts?
					// TODO: make this functiony (reusable)

				// output ideas
				// 1 list of transfers with freqs
				// 2 some kind of processing to make a network
				// 3 using input groups of interest map highways between those
				// 4 groups

}

void print_transfers(Node *super_tree, vector<Node *> *gene_trees,
		vector<string> *gene_tree_names, map<int, string> *reverse_label_map) {
	#pragma omp parallel for
	for(int i = 0; i < gene_trees->size(); i++) {
		Forest *MAF1 = NULL;
		Forest *MAF2 = NULL;
		Forest F1 = Forest(super_tree);
		Forest F2 = Forest((*gene_trees)[i]);
		cout << (*gene_tree_names)[i] << endl;
		if (sync_twins(&F1,&F2)) {
			int distance = rSPR_branch_and_bound_simple_clustering(F1.get_component(0), F2.get_component(0), &MAF1, &MAF2);
			expand_contracted_nodes(MAF1);
			expand_contracted_nodes(MAF2);
			/*
			cout << i << ": " << distance << endl;
			cout << "\tT1: "; F1.print_components();
			cout << "\tT2: "; F2.print_components();
			cout << "\tF1: "; MAF1->print_components_with_edge_pre_interval();
			cout << "\tF2: "; MAF2->print_components_with_edge_pre_interval();
			*/
			sync_af_twins(MAF1, MAF2);
			print_transfers(super_tree, &F1, &F2, MAF1, MAF2, reverse_label_map);
		}
		if (MAF1 != NULL)
			delete MAF1;
		if (MAF2 != NULL)
			delete MAF2;
	}
}

void print_transfers(Node *super_tree, Forest *F1, Forest *F2, Forest *MAF1, Forest *MAF2,
		map<int, string> *reverse_label_map) {
	int start = 1;
	if (MAF2->contains_rho())
		start = 0;
	for(int i = start; i < MAF2->num_components(); i++) {
		Node *F2_source = MAF2->get_component(i);
		Node *F1_source;
		Node *F1_target;
		if (!map_transfer(F2_source, F1, MAF2, &F1_source,
					&F1_target)) {
			continue;
		}

//		cout << F1_source->str_subtree() << endl;
//		cout << super_tree->find_by_prenum(F1_source->get_preorder_number())->str_subtree() << endl;
//		print_leaf_list(F1_source, reverse_label_map);
		print_leaf_list(super_tree->find_by_prenum(F1_source->get_preorder_number()), reverse_label_map);

	}
}
void add_transfers(list<transfer> *transfer_list, Forest *F1,
		Forest *F2, Forest *MAF1, Forest *MAF2);

bool transfer_compare (const transfer &a, const transfer &b) {
	return a.source_pre > b.source_pre;
}

void list_transfers(list<transfer> *transfer_list, Node *super_tree,
		Node *gene_tree) {
		Forest *MAF1 = NULL;
		Forest *MAF2 = NULL;
		Forest F1 = Forest(super_tree);
		Forest F2 = Forest(gene_tree);
		if (sync_twins(&F1,&F2)) {
			int distance = rSPR_branch_and_bound_simple_clustering(F1.get_component(0), F2.get_component(0), &MAF1, &MAF2);
			expand_contracted_nodes(MAF1);
			expand_contracted_nodes(MAF2);
#ifdef DEBUG_LGT
			cout << "\tT1: "; F1.print_components();
			cout << "\tT2: "; F2.print_components();
			//cout << "\tT1: "; F1.print_components_with_edge_pre_interval();
			//cout << "\tT2: "; F2.print_components_with_edge_pre_interval();
			cout << "\tF1: "; MAF1->print_components_with_edge_pre_interval();
			cout << "\tF2: "; MAF2->print_components_with_edge_pre_interval();
#endif
			sync_af_twins(MAF1, MAF2);
			add_transfers(transfer_list, &F1, &F2, MAF1, MAF2);
			transfer_list->sort(transfer_compare);
		}
		if (MAF1 != NULL)
			delete MAF1;
		if (MAF2 != NULL)
			delete MAF2;
}

void add_transfers(list<transfer> *transfer_list, Forest *F1,
		Forest *F2, Forest *MAF1, Forest *MAF2) {
	int start = 1;
	if (MAF2->contains_rho())
		start = 0;
	for(int i = start; i < MAF2->num_components(); i++) {
		Node *F2_source = MAF2->get_component(i);
		Node *F1_source;
		Node *F1_target;
		if (!map_transfer(F2_source, F1, MAF2, &F1_source,
					&F1_target)) {
			continue;
		}

		vector<int> source_child_pre_nums = vector<int>();
		vector<int> target_child_pre_nums = vector<int>();
		//if (F1_source->get_edge_pre_start() == -1) {
		
		
		for (auto c = F1_source->get_children().begin();
		     c != F1_source->get_children().end();
		     c++) {
		  if ((*c)->get_preorder_number() == F1_source->get_preorder_number()) {
		    /*
		    for (auto i = F1_source->get_children().begin();
			 i != F1_source->get_children().end();
			 i++) {
		      source_child_pre_nums.push_back((*i)->get_preorder_number());
		      }*/
		    vector<Node*> leaves = F1_source->find_leaves();
		    for (int i = 0; i < leaves.size(); i++) {
		      source_child_pre_nums.push_back(leaves[i]->get_preorder_number());
		    }
		    break;
		  }
		}
		for (auto c = F1_target->get_children().begin();
		     c != F1_target->get_children().end();
		     c++) {
		  if ((*c)->get_preorder_number() == F1_target->get_preorder_number()) {
		    /*
		    for (auto i = F1_target->get_children().begin();
			 i != F1_target->get_children().end();
			 i++) {
		      target_child_pre_nums.push_back((*i)->get_preorder_number());
		    }
		    */
		    vector<Node*> leaves = F1_target->find_leaves();
		    for (int i = 0; i < leaves.size(); i++) {
		      target_child_pre_nums.push_back(leaves[i]->get_preorder_number());
		    }

		    break;
		  }
		}
		
		/*

		else if (F2_source->get_edge_pre_start() == -1) {
		  for (list<Node*>::iterator i = F2_source->get_children().begin(); i != F2_source->get_children().end(); i++) {
		    source_child_pre_nums.push_back((*i)->get_twin()->get_preorder_number());
		  }
		}
		*/
		/*
		if (F1_target->get_edge_pre_start() == -1) {
		  for (list<Node*>::iterator i = F1_target->get_children().begin(); i != F1_target->get_children().end(); i++) {
		    target_child_pre_nums.push_back((*i)->get_preorder_number());
		  }
		}
		*/
		transfer_list->push_back(transfer(F1_source->get_preorder_number(),
					F1_target->get_preorder_number(),
						  F2_source->get_depth(),
						  source_child_pre_nums,
						  target_child_pre_nums));
	}
}

void print_leaf_list(Node *T, map<int, string> *reverse_label_map) {
	vector<Node *> leaves = T->find_leaves();
	vector<string> leaf_labels = vector<string>();
	if (!leaves.empty()) {
		for(int i = 0; i < leaves.size(); i++) {
			map<int, string>::iterator j = reverse_label_map->find(atoi(leaves[i]->str().c_str()));
			if (j != reverse_label_map->end()) {
				stringstream ss;
				ss << j->second;
				leaf_labels.push_back(ss.str());
			}
		}
		sort(leaf_labels.begin(), leaf_labels.end());
		cout << "\t";
		cout << leaf_labels[0];
		for(int i = 1; i < leaf_labels.size(); i++) {
			cout << "," << leaf_labels[i];
		}
		cout << endl;
	}
}

void print_leaf_list(Node *T) {
	vector<Node *> leaves = T->find_leaves();
	vector<string> leaf_labels = vector<string>();
	if (!leaves.empty()) {
		for(int i = 0; i < leaves.size(); i++) {
			leaf_labels.push_back(leaves[i]->str());
		}
		sort(leaf_labels.begin(), leaf_labels.end());
		cout << leaf_labels[0];
		for(int i = 1; i < leaf_labels.size(); i++) {
			cout << "," << leaf_labels[i];
		}
	}
}

bool map_transfer(Node *F2_source, Forest *F1, Forest *MAF2,
		Node **F1_source_out, Node **F1_target_out) {
	bool ret_val = false;
	if (F2_source->str() == "p")
		return ret_val;
	Node *F2_target = find_best_target(F2_source, MAF2);
	#ifdef DEBUG_LGT			
		cout << "\tF2s: " << F2_source->str_edge_pre_interval_subtree()
			<< endl;
		if (F2_target != NULL)
			cout << "\tF2t: " << F2_target->str_edge_pre_interval_subtree()
				<< endl;
		else
			cout << "\tF2t: p" << endl;
	#endif

	Node *F1_source = F2_source->get_twin();
			//save preorder numbers
	Node *F1_target;
	if (F2_target != NULL)
		F1_target = F2_target->get_twin();
	else
		F1_target = F1->get_component(0);
	#ifdef DEBUG_LGT			
		cout << "\tF1s: " << F1_source->str_edge_pre_interval_subtree()
			<< endl;
		if (F2_target != NULL)
			cout << "\tF1t: " << F1_target->str_edge_pre_interval_subtree()
				<< endl;
		else
			cout << "\tF1t: p" << endl;
		cout << endl;
	#endif
	string a = F1_source->get_name();
	string b = F1_target->get_name();

	*F1_source_out = F1_source;
	*F1_target_out = F1_target;

	ret_val = true;
	return ret_val;
}

Node *find_best_target(Node *source, Forest *AF) {
	Node *best_target = NULL;
	for(int i = 0; i < AF->num_components(); i++) {
		Node *target = AF->get_component(i);
		if (source != target)
			find_best_target(source, target, &best_target);
	}
	return best_target;
}

Node* find_best_target(Node *source, Node *target, Node **best_target) {
	if (target->get_edge_pre_start() <= source->get_preorder_number()
			&& target->get_edge_pre_end() >= source->get_preorder_number()
			&& (*best_target == NULL || target->get_edge_pre_start() >
					(*best_target)->get_edge_pre_start())) {
		*best_target = target;
	}
	list<Node *>::const_iterator c;
	for(c = target->get_children().begin();
			c != target->get_children().end(); c++) {
		find_best_target(source, *c, best_target);
	}
	return *best_target;
}

void add_lcas_to_groups(vector<int> *pre_to_group, Node *subtree) {
	list<Node *>::const_iterator c;
	int group = -1;
	for(c = subtree->get_children().begin();
			c != subtree->get_children().end(); c++) {
		add_lcas_to_groups(pre_to_group, *c);
		int child_group = (*pre_to_group)[(*c)->get_preorder_number()];
		if (group == -1)
			group = child_group;
		else if (group != child_group)
			group = -2;
	}
	if (group >= 0)
		(*pre_to_group)[subtree->get_preorder_number()] = group;
}

bool contains_any(vector<int>& A, vector<int>& B)
{
  for (int i = 0; i < B.size(); i++) {
    for (int j = 0; j < A.size(); j++) {
      if (B[i] == A[j]) {
	return true;
      }
    }
  }
  return false;
}

bool is_subset(vector<int>& A, vector<int>& B)
{
  if (B.size() == 0 && A.size() > 0) { return false; }
  for (int i = 0; i < A.size(); i++) {
    bool found = false;
    for (int j = 0; j < B.size(); j++) {
      if (B[j] == A[i]) {
	found = true;
	break;
      }
    }
    if (!found) {
      return false;
    }
  }
  return true;
}

void show_moves(Node *T1, Node *T2, map<string, int> *label_map,
		map<int, string> *reverse_label_map) {
	T1->preorder_number();
	T2->preorder_number();
	T1->edge_preorder_interval();
	T2->edge_preorder_interval();
	int num_nodes = T1->size();
	int distance = rSPR_branch_and_bound_simple_clustering(T1, T2);
	int current_distance = distance;
	T1->numbers_to_labels(reverse_label_map);
	cout << T1->str_subtree() << endl;
	T1->labels_to_numbers(label_map, reverse_label_map);
	while (current_distance > 0) {
		list<transfer> transfer_list = list<transfer>();
		list_transfers(&transfer_list, T1, T2);
		T1->numbers_to_labels(reverse_label_map);
	choose_different_tranfer:
		transfer trans = transfer_list.front();
		transfer_list.pop_front();
		Node *s = T1->find_by_prenum_full(trans.source_pre);
		Node *t = T1->find_by_prenum_full(trans.target_pre);

#ifdef DEBUG_LGT
		cout << trans.source_pre << endl;
		cout << trans.target_pre << endl;
		print_leaf_list(s);
		cout << " : ";
		print_leaf_list(t);
		cout << endl;
#endif

		//TODO (ben): try updating the preorder interval when contracting

		if (MULTIFURCATING) {
		  vector<Node*> source_nodes_to_expand = vector<Node*>();
		  
		  if (trans.source_child_pre_nums.size() > 0) {
		    //move up until all the preorders are in s's subtree
		    s = s->parent();
		    vector<int> leaf_preorders = s->find_leaf_preorders();
		    bool leaf_is_subset = is_subset(trans.source_child_pre_nums, leaf_preorders);
		    while (!leaf_is_subset){
			s = s->parent();
			leaf_preorders = s->find_leaf_preorders();
			leaf_is_subset = is_subset(trans.source_child_pre_nums, leaf_preorders);
		    }
		    //cout << "New S " << s->get_preorder_number() << endl;
		    //get the children that have the preorders in their subtree
		    for (auto i = s->get_children().begin(); i != s->get_children().end(); i++) {
		      vector<int> this_childs_preorder_list = (*i)->find_leaf_preorders();
		      if (contains_any(trans.source_child_pre_nums, this_childs_preorder_list)){
			source_nodes_to_expand.push_back(*i);
		      }
		    }

		    //before we expand, we need to check if the target is in any of the nodes_to_expand's subtrees
		  }

		  
		  //want to get the subset of children
		  //the prenum will be the child though, so move up one
		  vector<Node*> target_nodes_to_expand = vector<Node*>();
		  if (trans.target_child_pre_nums.size() > 0) {
		    //move up until all the preorders are in t's subtree
		    t = t->parent();
		    vector<int> leaf_preorders = t->find_leaf_preorders();
		    bool leaf_is_subset = is_subset(trans.target_child_pre_nums, leaf_preorders);
		    while (!leaf_is_subset){
			t = t->parent();
			leaf_preorders = t->find_leaf_preorders();
			leaf_is_subset = is_subset(trans.target_child_pre_nums, leaf_preorders);
		      }
		    
		    //get the children that have the preorders in their subtree
		    for (auto i = t->get_children().begin(); i != t->get_children().end(); i++) {
		      vector<int> this_childs_preorder_list = (*i)->find_leaf_preorders();
		      if (contains_any(trans.target_child_pre_nums, this_childs_preorder_list)){
			target_nodes_to_expand.push_back(*i);
		      }
		    }		   
		  }

		  //check if the target will be in any of the source nodes subtrees
		    //if it is then we are attaching it to itself and it is an invalid spr
		  if (source_nodes_to_expand.size() > 0) {		    
		    for (int i = 0; i < source_nodes_to_expand.size(); i++) {
		      if (source_nodes_to_expand[i]->find_by_prenum_full(t->get_preorder_number()) != NULL) {
			cout << "INVALID MOVE, RETRYING..." << endl;
			goto choose_different_tranfer;
		      }
		    }
		  }
		  else {
		    if (s->find_by_prenum_full(t->get_preorder_number()) != NULL) {
		      cout << "INVALID MOVE, RETRYING..." << endl;
		      goto choose_different_tranfer;
		    }
		  }
		  //cout << "Big Decision T0.0 ";
			//T1->print_subtree();
			if (trans.source_child_pre_nums.size() > 0) {
		    //if there are children who don't have the correct preorders, expand the others
		    if (source_nodes_to_expand.size() < s->get_children().size()) {
		      Node* new_s = new Node();
		      s->add_child(new_s);
		      for (int i = 0; i < source_nodes_to_expand.size(); i++) {
			new_s->add_child(source_nodes_to_expand[i]);
		      }
		      s = new_s;
		    }
		  }
		  if (trans.target_child_pre_nums.size() > 0) {
			//if there are children who don't have the correct preorders, expand the others
		    if (target_nodes_to_expand.size() < t->get_children().size()){
		      Node* new_t = new Node();
		      t->add_child(new_t);
		      for (int i = 0; i < target_nodes_to_expand.size(); i++) {
			new_t->add_child(target_nodes_to_expand[i]);
		      }
		      t = new_t;
		    }
		  }
		  if (s->parent() != NULL &&
		      s->parent()->get_children().size() > 2) {
		    Node* new_s = new Node();
		    s->parent()->add_child(new_s);		    
		    new_s->add_child(s);
		    new_s->add_child(s->parent()->parent()->lchild());
		    //s = new_s;
		  }
			#ifdef DEBUG_LGT
		  print_leaf_list(s);
		  cout << " : ";
		  print_leaf_list(t);
		  cout << endl;
			#endif

		  s->spr_mult(t);
		}
		else {
		  s->spr(t);
		}

		print_leaf_list(s);
		cout << " : ";
		print_leaf_list(t);
		cout << endl;
		cout << T1->str_subtree() << endl;
		T1->labels_to_numbers(label_map, reverse_label_map);

		// cleanup T1
		// TODO: the SPR does something odd to the tree
		// find a way to remove this hack of creating a new tree each time
		string str = T1->str_subtree();
		T1 = build_tree(str);
		T1->preorder_number();
		T1->edge_preorder_interval();

		current_distance--;
		int distance = rSPR_branch_and_bound_simple_clustering(T1, T2);

	}

	Forest F1 = Forest(T1);
	Forest F2 = Forest(T2);
	int approx_spr;
	if (MULTIFURCATING) {
	  approx_spr = rSPR_worse_3_mult_approx(&F1, &F2);
	}
	else {
	  approx_spr = rSPR_worse_3_approx(&F1, &F2);
	}
	if (approx_spr > 0) {
		cout << "WARNING: final tree does not match T2" << endl;
	}
	return;
}

#endif
