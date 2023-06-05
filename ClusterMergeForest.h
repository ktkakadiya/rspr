/*******************************************************************************
ClusterMergeForest.h


Copyright 2011-2014 Chris Whidden
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
#ifndef INCLUDE_CLUSTERMERGEFOREST

#define INCLUDE_CLUSTERMERGEFOREST

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

using namespace std;

class ClusterMergeForest {

    private:
        int cur_cluster_index;
        vector<pair<int, int>> unsolved_cluster_prenums;
        vector<int> cluster_parent_index;
        vector<int> cluster_prenums_f1;
        vector<int> cluster_prenums_f2;
	    map<int, list<pair<Forest,Forest>>> map_cluster_mafs;

	public:
	ClusterMergeForest(){
		init(0);
	}

    ClusterMergeForest(int cluster_size){
		init(cluster_size);
	}

	void init(int cluster_size) {
        cur_cluster_index = 0;
		cluster_parent_index.assign(cluster_size, -1);
        cluster_prenums_f1.assign(cluster_size, -1);
        cluster_prenums_f2.assign(cluster_size, -1);
	}

	~ClusterMergeForest() {
		unsolved_cluster_prenums.clear();
		cluster_parent_index.clear();
		cluster_prenums_f1.clear();
		cluster_prenums_f2.clear();
        map_cluster_mafs.clear();
	}

    /**
     * @brief Assign which cluster is parent of which other clusters and save cluster points prenums
     * 
     * @param cluster_node_prenum 
     * @param cluster_node_twin_prenum 
     */
    void assign_cluster_parent(int cluster_node_prenum, int cluster_node_twin_prenum) {
        cur_cluster_index++;
        if(unsolved_cluster_prenums.size() > 0){
            while(unsolved_cluster_prenums.back().second > cluster_node_prenum){
                cluster_parent_index[unsolved_cluster_prenums.back().first] = cur_cluster_index;
                unsolved_cluster_prenums.pop_back();
            }
        }
        unsolved_cluster_prenums.push_back(make_pair(cur_cluster_index, cluster_node_prenum));
        cluster_prenums_f1[cur_cluster_index] = cluster_node_prenum;
        cluster_prenums_f2[cur_cluster_index] = cluster_node_twin_prenum;	
    }

    Forest* merge_agreement_forests(Forest *upper_forest, Forest *lower_forest, 
                                        bool has_sep_comp, int lower_cluster_prenum) {
        Forest* merged_forest = new Forest(upper_forest);
        Node* cluster_node = merged_forest->get_contracted_node_with_prenum(lower_cluster_prenum);
        if(!cluster_node){
            cout << "Cluster node not found!" << endl;
            return merged_forest;
        }

        vector<Node *> components = lower_forest->components;
        bool has_rho = lower_forest->contains_rho();
        
        if(has_sep_comp){
            merged_forest->remove_component_with_prenum(lower_cluster_prenum);
            vector<Node *>::iterator it;	
            for(it = components.begin(); it != components.end(); it++) {
                if ((*it)->str() != "p"){
                    merged_forest->add_component(new Node(**it));
                }
            }
        }
        else {
            vector<Node *>::iterator it;
            int i=0;
            for(it = components.begin(); it != components.end(); it++) {
                if(!has_rho){
                    if(i == 0){		
                        Node* lower_clstr_parent = cluster_node->parent();
                        Node* lower_comp = new Node(**it, lower_clstr_parent);
                        if(lower_clstr_parent != NULL){
                            lower_clstr_parent->replace_contracted_child(cluster_node, lower_comp);
                        }
                    }
                    else{
                        merged_forest->add_component(new Node(**it));
                    }
                }
                i++;
            }
        }
        return merged_forest;
    }

    list<pair<Forest,Forest>> merge_cluster_forests(list<pair<Forest,Forest>> upper_cluster_mafs,
                                                    list<pair<Forest,Forest>> lower_cluster_mafs,
                                                    int lower_cluster_prenum_f1, int lower_cluster_prenum_f2){
        list<pair<Forest,Forest>> merged_mafs = list<pair<Forest,Forest>>();
        for (const auto& upper_maf_pair : upper_cluster_mafs) {
            Forest uF1 = upper_maf_pair.first;
            Forest uF2 = upper_maf_pair.second;			

            bool has_sep_comp = uF1.has_component_with_prenum(lower_cluster_prenum_f1);
            for (const auto& lower_maf_pair : lower_cluster_mafs) {
                Forest lF1 = lower_maf_pair.first;
                Forest lF2 = lower_maf_pair.second;		
                Forest* merged_maf1 = merge_agreement_forests(&uF1, &lF1, has_sep_comp, lower_cluster_prenum_f1);
                Forest* merged_maf2 = merge_agreement_forests(&uF2, &lF2, has_sep_comp, lower_cluster_prenum_f2);
                merged_maf1->print_components();
                merged_maf2->print_components();
                merged_mafs.push_back(make_pair(Forest(merged_maf1), Forest(merged_maf2)));
            }
        }
        return merged_mafs;
    }

    /**
     * @brief Update merged agreement forests
     * 
     * @param cluster_idx 
     * @param num_clusters 
     * @param extAFs 
     */
    void update_merged_afs(int cluster_idx, int num_clusters, list<pair<Forest,Forest>> *extAFs){
        int child_cluster_index = cluster_idx - 1;
        int cur_parent_idx = cluster_idx;
        if (cluster_idx >= num_clusters - 1) {
            cur_parent_idx = -1;
        }
        map_cluster_mafs[cluster_idx] = *extAFs;
        while(child_cluster_index > 0){
            if(cluster_parent_index[child_cluster_index] == cur_parent_idx){
                map_cluster_mafs[cluster_idx] = merge_cluster_forests(map_cluster_mafs[cluster_idx],
                                    map_cluster_mafs[child_cluster_index],
                                    cluster_prenums_f1[child_cluster_index],
                                    cluster_prenums_f2[child_cluster_index]);
            }
            child_cluster_index--;
        }
    }

    void print_merged_agreement_forest(int cluster_index, map<string, int> *label_map,
                                                map<int, string> *reverse_label_map){
        vector<string> vecComponents;
        unordered_set<string> setMAFs;

        cout << endl << endl << "FOUND ANSWERS " << map_cluster_mafs[cluster_index].size() << endl;

        list<pair<Forest,Forest> >::iterator x = map_cluster_mafs[cluster_index].begin();
        for (; x != map_cluster_mafs[cluster_index].end(); x++) {
            if (label_map != NULL && reverse_label_map != NULL) {
                x->first.numbers_to_labels(reverse_label_map);
                x->second.numbers_to_labels(reverse_label_map);
            }
            
            string strMAF = (x->first).add_vec_components(&vecComponents);
            setMAFs.insert(strMAF);

            cout << "\tMerged T1: ";
            x->first.print_components();
            cout << "\tMerged T2: ";
            x->second.print_components();
            if (label_map != NULL && reverse_label_map != NULL) {
                x->first.labels_to_numbers(label_map, reverse_label_map);
                x->second.labels_to_numbers(label_map, reverse_label_map);
            }
        }
        cout << "FOUND MERGED UNIQUE ANSWERS " << setMAFs.size() << endl;
    }
};

#endif