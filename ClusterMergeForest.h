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

class MergedSolution {
    public:
    Forest *forest;
    Node* merging_point;

    MergedSolution(Forest* forest, Node* merging_point) {
        this->forest = forest;
        this->merging_point = merging_point;   
    }

    void set_merging_point(Node* merging_point){
        this->merging_point = merging_point;
    }
};

class ClusterMergeForest {

    private:
        //temporary stack used to store currently unsolved cluster to get parent cluster assigned
        vector<pair<int, int>> unsolved_cluster_prenums;

        //Indexes of parent cluster
        vector<int> cluster_parent_index;

        //Preorder numbers of clusters in F1
        vector<int> cluster_prenums_f1;

        //Preorder numbers of clusters in F1
        vector<pair<int, int>> cluster_edge_interval_f1;

        //Preorder numbers of clusters in F2
        vector<int> cluster_prenums_f2;

        //Preorder numbers of clusters in F2
        vector<pair<int, int>> cluster_edge_interval_f2;

        //Start index of clusters extra MAFs that are added with k+1 distance
        vector<int> cluster_exact_spr;

        //Map of merged MAFs of clusters at given index
	    map<int, list<pair<Forest,Forest>>> map_cluster_mafs;

        //Map used to store whether cluster at given index has any node in upper cluster
        vector<bool> map_has_cluster_node;

	public:
	ClusterMergeForest(){
		init(0);
	}

    ClusterMergeForest(int cluster_size){
		init(cluster_size);
	}

	void init(int cluster_size) {
		cluster_parent_index.assign(cluster_size, -1);
        cluster_prenums_f1.assign(cluster_size, -1);
        cluster_edge_interval_f1.assign(cluster_size, make_pair(-1, -1));
        cluster_prenums_f2.assign(cluster_size, -1);
        cluster_edge_interval_f2.assign(cluster_size, make_pair(-1, -1));
        cluster_exact_spr.assign(cluster_size, 0);
        map_has_cluster_node.assign(cluster_size, true);
	}

	~ClusterMergeForest() {
		unsolved_cluster_prenums.clear();
		cluster_parent_index.clear();
		cluster_prenums_f1.clear();
		cluster_edge_interval_f1.clear();
		cluster_prenums_f2.clear();
		cluster_edge_interval_f2.clear();
		cluster_exact_spr.clear();
        map_cluster_mafs.clear();
        map_has_cluster_node.clear();
	}

    /**
     * @brief Assign which cluster is parent of which other clusters and save cluster points prenums
     * 
     * @param cluster_node_prenum 
     * @param cluster_node_twin_prenum 
     * @param cur_cluster_index
     */
    void assign_cluster_parent(Node* node, Node* twin, int cur_cluster_index) {
        int cluster_node_prenum = node->get_preorder_number();
        int cluster_node_twin_prenum = twin->get_preorder_number();
        while(unsolved_cluster_prenums.size() > 0 &&
                            unsolved_cluster_prenums.back().second > cluster_node_prenum){
            cluster_parent_index[unsolved_cluster_prenums.back().first] = cur_cluster_index;
            unsolved_cluster_prenums.pop_back();
        }
        unsolved_cluster_prenums.push_back(make_pair(cur_cluster_index, cluster_node_prenum));
        cluster_prenums_f1[cur_cluster_index] = cluster_node_prenum;
        cluster_prenums_f2[cur_cluster_index] = cluster_node_twin_prenum;
        cluster_edge_interval_f1[cur_cluster_index] = make_pair(node->get_edge_pre_start(),
                                                                node->get_edge_pre_end());
        cluster_edge_interval_f2[cur_cluster_index] = make_pair(twin->get_edge_pre_start(),
                                                                twin->get_edge_pre_end());

        /*cout << "Cluster F1 prenums " << endl;
        for(auto it=cluster_prenums_f1.begin(); it != cluster_prenums_f1.end(); it++){
            cout << " " << *it;
        }
        cout << endl;
        cout << "Cluster F2 prenums " << endl;
        for(auto it=cluster_prenums_f2.begin(); it != cluster_prenums_f2.end(); it++){
            cout << " " << *it;
        }
        cout << endl;
        cout << "Cluster Parent index " << endl;
        for(auto it=cluster_parent_index.begin(); it != cluster_parent_index.end(); it++){
            cout << " " << *it;
        }
        cout << endl;*/
    }

    MergedSolution* merge_agreement_forests(Forest *upper_forest, Forest *lower_forest, 
                                        bool has_sep_comp, int lower_cluster_prenum,
                                        pair<int, int> lower_cluster_range,
                                        bool has_cluster_node) {

        Forest* merged_forest = new Forest(upper_forest);
        vector<Node *> components = lower_forest->components;
        bool has_rho = lower_forest->contains_rho();

        MergedSolution *soln = new MergedSolution(merged_forest, NULL);

        if(has_sep_comp){
            
            //Update lower cluster root component preorder range
            Node* upper_component = merged_forest->get_component_with_prenum(lower_cluster_prenum);
            if(upper_component && components.size() > 0 && components[0]->str() != "p"){
                components[0]->set_preorder_number(upper_component->get_preorder_number());
                components[0]->copy_edge_pre_interval(upper_component);
            }

            merged_forest->remove_component_with_prenum(lower_cluster_prenum);
            vector<Node *>::iterator it;	
            for(it = components.begin(); it != components.end(); it++) {
                if ((*it)->str() != "p"){
                    merged_forest->add_component(new Node(**it));
                }
            }
        }
        else if(!has_cluster_node){
            if(has_rho){
                vector<Node *>::iterator it;	
                for(it = components.begin(); it != components.end(); it++) {
                    if ((*it)->str() != "p"){
                        merged_forest->add_component(new Node(**it));
                    }
                }
            }
            else{
                Node* cluster_node = merged_forest->get_best_node_with_prenum(lower_cluster_prenum);
                if(!cluster_node){
                    cout << "Cluster node not found 1!" << endl;
                    return soln;
                }

                soln->set_merging_point(cluster_node);
                vector<Node *>::iterator it;
                int i=0;
                for(it = components.begin(); it != components.end(); it++) {
                    if(i == 0){
                        cluster_node->place_contracted_child(new Node(**it));
                    }
                    else{
                        merged_forest->add_component(new Node(**it));
                    }
                    i++;
                }
            }
        }
        else if(has_cluster_node){
            Node* cluster_node = merged_forest->get_contracted_node_with_prenum_range(lower_cluster_range);
            if(!cluster_node){
                cout << "Cluster node not found 2!" << endl;
                return soln;
            }
            
            vector<Node *>::iterator it;
            int i=0;
            for(it = components.begin(); it != components.end(); it++) {
                if(!has_rho){
                    if(i == 0){		
                        Node* lower_clstr_parent = cluster_node->parent();
                        if(lower_clstr_parent != NULL){
                            if(!cluster_node->is_contracted()){
                                lower_clstr_parent->contract_sibling_pair_undoable();
                            }
                            Node* lower_comp = new Node(**it, lower_clstr_parent);
                            lower_clstr_parent->replace_contracted_child(cluster_node, lower_comp);
                        }
                    }
                    else{
                        merged_forest->add_component(new Node(**it));
                    }
                }
                else{
                    if((*it)->str() == "p"){		
                        Node* lower_clstr_parent = cluster_node->parent();
                        if(lower_clstr_parent != NULL){
                            lower_clstr_parent->replace_contracted_child(cluster_node, NULL);
                            lower_clstr_parent->contract(true);
                        }
                    }
                    else{
                        merged_forest->add_component(new Node(**it));
                    }
                }
                i++;
            }
        }
        return soln;
    }

    list<pair<Forest,Forest>> merge_cluster_forests(list<pair<Forest,Forest>> upper_cluster_mafs,
            list<pair<Forest,Forest>> lower_cluster_mafs, int lower_cluster_prenum_f1, 
            pair<int, int> lower_cluster_range_f1, int lower_cluster_prenum_f2,
            pair<int, int> lower_cluster_range_f2, bool has_cluster_node, int total_spr){

        list<pair<Forest,Forest>> merged_mafs = list<pair<Forest,Forest>>();
        for (const auto& upper_maf_pair : upper_cluster_mafs) {
            Forest uF1 = upper_maf_pair.first;
            Forest uF2 = upper_maf_pair.second;
            
            bool has_sep_comp = false;
            if(has_cluster_node)
                has_sep_comp = uF1.has_component_with_prenum(lower_cluster_prenum_f1);

            for (const auto& lower_maf_pair : lower_cluster_mafs) {
                Forest lF1 = lower_maf_pair.first;
                Forest lF2 = lower_maf_pair.second;

                MergedSolution *merged_soln1 = merge_agreement_forests(&uF1, &lF1, has_sep_comp, 
                                                    lower_cluster_prenum_f1, lower_cluster_range_f1,
                                                    has_cluster_node);
                MergedSolution *merged_soln2 = merge_agreement_forests(&uF2, &lF2, has_sep_comp, 
                                                    lower_cluster_prenum_f2, lower_cluster_range_f2,
                                                    has_cluster_node);

                bool bValid = false;
                if(!merged_soln1->merging_point && !merged_soln2->merging_point){
                    bValid = true;
                }
                else if(merged_soln1->merging_point && merged_soln2->merging_point){
                    Forest mf1 = Forest(merged_soln1->merging_point);
                    Forest mf2 = Forest(merged_soln2->merging_point);
                    sync_twins(&mf1, &mf2);
                    Node* merging_point1 = mf1.get_component(0);
                    Node* merging_point2 = mf2.get_component(0);
                    if(merging_point1->get_twin() == merging_point2){
                        bValid = true;
                    }
                }  

                if(merged_soln1->forest->num_components() != total_spr+1 || 
                        merged_soln2->forest->num_components() != total_spr+1){
                    bValid = false;
                }

                if(bValid){
                    Forest *merged_maf1 = merged_soln1->forest;
                    Forest *merged_maf2 = merged_soln2->forest;
                    //merged_maf1->print_components();
                    //merged_maf2->print_components();
                    merged_mafs.push_back(make_pair(Forest(merged_maf1), Forest(merged_maf2)));
                }
            }
        }
        return merged_mafs;
    }

    vector<int> get_forests_with_rho(list<pair<Forest,Forest>> lower_cluster_mafs){
        vector<int> indexes;
        int index=0;
        for (const auto& lower_maf_pair : lower_cluster_mafs) {
            Forest uF1 = lower_maf_pair.first;
            if(uF1.contains_rho()){
                indexes.push_back(index);
                index++;
            }
        }
        return indexes;
    }

    /**
     * @brief Update merged agreement forests
     * 
     * @param cluster_idx 
     * @param num_clusters 
     * @param extAFs 
     */
    void update_merged_forests(int cluster_idx, int num_clusters, list<pair<Forest,Forest>> *extAFs,
                                                                bool has_cluster_node, int exact_spr){
        int child_cluster_index = cluster_idx - 1;
        int cur_parent_idx = cluster_idx;
        if (cluster_idx >= num_clusters - 1) {
            cur_parent_idx = -1;
        }
        map_cluster_mafs[cluster_idx] = *extAFs;
        map_has_cluster_node[cluster_idx] = has_cluster_node;
        cluster_exact_spr[cluster_idx] = exact_spr;
        while(child_cluster_index > 0){
            if(cluster_parent_index[child_cluster_index] == cur_parent_idx){
                int total_spr = cluster_exact_spr[cluster_idx] + cluster_exact_spr[child_cluster_index];
                cluster_exact_spr[cluster_idx] = total_spr;
                map_cluster_mafs[cluster_idx] = merge_cluster_forests(map_cluster_mafs[cluster_idx],
                                    map_cluster_mafs[child_cluster_index],
                                    cluster_prenums_f1[child_cluster_index],
                                    cluster_edge_interval_f1[child_cluster_index],
                                    cluster_prenums_f2[child_cluster_index],
                                    cluster_edge_interval_f2[child_cluster_index],
                                    map_has_cluster_node[child_cluster_index], total_spr);
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

    void add_non_duplicate_maf(list<pair<Forest,Forest>> *extAFs, list<pair<Forest,Forest>> *extRhoAFs){
        vector<string> vecComponents;
        unordered_set<string> setMAFs;
        for (list<pair<Forest,Forest> >::iterator x = extAFs->begin(); x != extAFs->end(); x++) {
            string strMAF = (x->first).add_vec_components(&vecComponents);
            setMAFs.insert(strMAF);
        }

        for (list<pair<Forest,Forest> >::iterator x = extRhoAFs->begin(); x != extRhoAFs->end(); x++) {
            string strMAF = (x->first).add_vec_components(&vecComponents, false);
            std::pair<std::unordered_set<string>::iterator, bool> insertResult = setMAFs.insert(strMAF);
            if(insertResult.second){
                extAFs->push_back(make_pair(x->first, x->second));
            }
        }
    }
};

#endif