#include <limits>
#include <TProcessID.h>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <algorithm>

#include "CCMAnalysis/CCMIO/CCMRootHandle.h"

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

CCMRootHandle::CCMRootHandle(TFile * input_file_ptr) :
    input_file_ptr(input_file_ptr), tree_size(0), current_index(0), tree_name_specified(false) {
    if(input_file_ptr == nullptr)
        throw std::runtime_error("Invalid TFile pointer passed to CCMRootHandle()");

    size_t tsize = 150 * 1000 * 1000 * 1000; //set the limit of the event tree to be large (GB)
    TTree::SetMaxTreeSize(tsize);
    LoadAllBranches();
}

CCMRootHandle::CCMRootHandle(TFile * input_file_ptr, std::vector<std::string> branch_names) :
    input_file_ptr(input_file_ptr), tree_size(0), current_index(0), tree_name_specified(false) {
    if(input_file_ptr == nullptr)
        throw std::runtime_error("Invalid TFile pointer passed to CCMRootHandle()");

    size_t tsize = 150 * 1000 * 1000 * 1000; //set the limit of the event tree to be large (GB)
    TTree::SetMaxTreeSize(tsize);
    LoadBranches(branch_names);
}

CCMRootHandle::CCMRootHandle(TFile * input_file_ptr, std::string const & tree_name, std::vector<std::string> branch_names) :
    input_file_ptr(input_file_ptr), tree_size(0), current_index(0), tree_name_specified(true), specified_tree_name(tree_name) {
    if(input_file_ptr == nullptr)
        throw std::runtime_error("Invalid TFile pointer passed to CCMRootHandle()");

    size_t tsize = 150 * 1000 * 1000 * 1000; //set the limit of the event tree to be large (GB)
    TTree::SetMaxTreeSize(tsize);
    LoadBranches(tree_name, branch_names);
}

std::set<std::string> CCMRootHandle::GetTreeNames() {
    TObjArray * key_array = input_file_ptr->GetListOfKeys();
    std::set<std::string> tree_names;
    for(TObjArray * obj_ptr : *key_array) {
        TKey * key_ptr = (TKey *) obj_ptr;
        std::string branch_name = key_ptr->GetClassName();
        if(branch_name == "TTree") {
            tree_names.insert(branch_name);
        }
    }
    return tree_names;
}

CCMRootHandle::CCMRootHandle(TFile * input_file_ptr, std::string const & tree_name) : CCMRootHandle(tree_name, std::vector<std::string>()) {};

void CCMRootHandle::LoadAllBranches() {
    std::set<std::string> tree_names = GetTreeNames();

    std::set<std::string> branch_names;
    for(std::string const & tree_name : tree_names) {
        std::set<std::string> tree_branch_names;
        TTree * tree_ptr;
        input_file_ptr->GetObject(tree_name.c_str(), tree_ptr);
        TObjArray * branch_array = input_file_ptr->GetListOfBranches();
        for(TObjArray * obj_ptr : *branch_array) {
            TBranch * branch_ptr = (TBranch *) obj_ptr;
            std::string branch_name = branch_ptr->GetName();
            tree_branch_names.insert(branch_name);
        }
        for(std::string const & branch_name : tree_branch_names) {
            if(branch_names.count(branch_name) != 0)
                throw std::runtime_error("Multiple trees in the root file have branches with the same name");
            TBranch * branch_ptr = tree_ptr->GetBranch(branch_name.c_str());
            branches.insert({branch_name, branch_ptr});
            num_entries.insert({branch_name, branch_ptr->GetEntries()});
            branch_names.insert(branch_name);
            if(tree_size < branch_ptr->GetEntries())
                tree_size = branch_ptr->GetEntries();
        }
    }
}

void CCMRootHandle::LoadBranches(std::vector<std::string> const & branch_name_vector) {
    std::set<std::string> tree_names = GetTreeNames();

    std::set<std::string> branch_names(branch_name_vector.begin(), branch_name_vector.end());
    std::set<std::string> found_branch_names;
    for(std::string const & tree_name : tree_names) {
        TTree * tree_ptr;
        input_file_ptr->GetObject(tree_name.c_str(), tree_ptr);
        for(std::string const & branch_name : branch_names) {
            TBranch * branch_ptr = tree_ptr->GetBranch(branch_name.c_str());
            if(branch_ptr == nullptr)
                continue;
            if(found_branch_names.count(branch_name) != 0)
                throw std::runtime_error("Multiple trees in the root file have branches with the same name");
            branches.insert({branch_name, branch_ptr});
            num_entries.insert({branch_name, branch_ptr->GetEntries()});
            if(tree_size < branch_ptr->GetEntries())
                tree_size = branch_ptr->GetEntries();
        }
    }
}

void CCMRootHandle::LoadBranches(std::string const & tree_name, std::vector<std::string> const & branch_name_vector) {
    std::set<std::string> tree_names = GetTreeNames();
    if(tree_names.count(tree_names) == 0)
        throw std::runtime_error("Tree name \"" + tree_name + "\" not found in root file");

    std::set<std::string> branch_names(branch_name_vector.begin(), branch_name_vector.end());
    TTree * tree_ptr;
    input_file_ptr->GetObject(tree_name.c_str(), tree_ptr);
    for(std::string const & branch_name : branch_names) {
        TBranch * branch_ptr = tree_ptr->GetBranch(branch_name.c_str());
        if(branch_ptr == nullptr)
            throw std::runtime_error("Could not find branch named \"" + branch_name + "\" in tree named \"" + tree_name + "\"");
        branches.insert({branch_name, branch_ptr});
        num_entries.insert({branch_name, branch_ptr->GetEntries()});
        if(tree_size < branch_ptr->GetEntries())
            tree_size = branch_ptr->GetEntries();
    }
}

void CCMRootHandle::LoadBranch(std::string const & branch_name) {
    if(tree_name_specified) {
        LoadBranch(tree_name, branch_name);
        return;
    }
    std::set<std::string> tree_names = GetTreeNames();

    std::set<std::string> found_branch_names;
    for(std::string const & tree_name : tree_names) {
        TTree * tree_ptr;
        input_file_ptr->GetObject(tree_name.c_str(), tree_ptr);
        TBranch * branch_ptr = tree_ptr->GetBranch(branch_name.c_str());
        if(branch_ptr == nullptr)
            continue;
        if(found_branch_names.count(branch_name) != 0)
            throw std::runtime_error("Multiple trees in the root file have branches with the same name");
        branches.insert({branch_name, branch_ptr});
        num_entries.insert({branch_name, branch_ptr->GetEntries()});
        if(tree_size < branch_ptr->GetEntries())
            tree_size = branch_ptr->GetEntries();
    }
}

void CCMRootHandle::Advance(size_t n) {
    current_index += n;
}

size_t CCMRootHandle::Rewind(size_t n) {
    if(current_index >= n)
        current_index -= n;
    else
        current_index = 0;
}

bool CCMRootHandle::More() {
    return current_index < (tree_size - 1);
}

bool CCMRootHandle::Valid() {
    return current_index < tree_size;
}
