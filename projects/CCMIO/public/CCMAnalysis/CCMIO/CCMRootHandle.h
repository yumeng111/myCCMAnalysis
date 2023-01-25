#ifndef CCM_CCMRootHandle_H
#define CCM_CCMRootHandle_H

#include <map>
#include <set>
#include <list>
#include <memory>
#include <vector>
#include <stdint.h>
#include <iostream>
#include <boost/shared_ptr.hpp>

#include <icetray/I3FrameObject.h>
#include <icetray/serialization.h>

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

class CCMRootHandle {
protected:
    TFile * input_file_ptr;
    size_t tree_size;
    size_t current_index;
    bool tree_name_specified;
    std::string specified_tree_name;

    std::map<std::string, TBranch *> branches;
    std::map<std::string, size_t> num_entries;

    std::set<std::string> GetTreeNames();
    void LoadBranch(std::string const & branch_name);
    void LoadBranches(std::vector<std::string> const & branch_name);
    void LoadBranches(std::string const & tree_name, std::vector<std::string> const & branch_name);
    void LoadAllBranches();
public:
    CCMRootHandle(TFile * input_file_ptr);
    CCMRootHandle(TFile * input_file_ptr, std::vector<std::string> const & branch_names);
    CCMRootHandle(TFile * input_file_ptr, std::string const & tree_name, std::vector<std::string> const & branch_names);
    CCMRootHandle(TFile * input_file_ptr, std::string const & tree_name);

    template<typename T>
    void Get(std::string const & branch_name, boost::shared_ptr<T> data_ptr) {
        std::cerr << "Getting entry number " << current_index << " on branch " << branch_name << std::endl;
        if(branches.count(branch_name) == 0)
            LoadBranch(branch_name);
        TBranch * branch_ptr = branches[branch_name];
        if(branch_ptr == nullptr or current_index >= num_entries[branch_name]) {
            data_ptr = nullptr;
            return;
        }
        std::cerr << "Branch has " << branch_ptr->GetEntries() << " entries" << std::endl;
        T * data = new T();
        branch_ptr->SetAddress(&data);
        branch_ptr->GetEntry(current_index);
        data_ptr = boost::shared_ptr<T>(data);
    }

    template<typename T>
    boost::shared_ptr<T> Get(std::string const & branch_name) {
        std::cerr << "Getting entry number " << current_index << " on branch " << branch_name << std::endl;
        if(branches.count(branch_name) == 0)
            LoadBranch(branch_name);
        TBranch * branch_ptr = branches[branch_name];
        if(branch_ptr == nullptr or current_index >= num_entries[branch_name]) {
            return nullptr;
        }
        T * data = new T();
        std::cerr << "Branch has " << branch_ptr->GetEntries() << " entries" << std::endl;

        TTree * tree_ptr = nullptr;
        input_file_ptr->GetObject("EventTree", tree_ptr);
        TBranch * alt_branch_ptr = nullptr;
        alt_branch_ptr = tree_ptr->GetBranch("rawData");
        std::cerr << "stored branch_ptr: " << branch_ptr << std::endl;
        std::cerr << "new branch_ptr: " << alt_branch_ptr << std::endl;
        alt_branch_ptr->SetAddress(&data);
        alt_branch_ptr->GetEntry(current_index);
        branch_ptr->SetAddress(&data);
        branch_ptr->GetEntry(current_index);
        boost::shared_ptr<T> data_ptr(data);
        return data_ptr;
    }

    void Advance(size_t n = 1);
    void Rewind(size_t n = 1);
    bool More();
    bool Valid();
};

I3_POINTER_TYPEDEFS(CCMRootHandle);

#endif // CCM_CCMRootHandle_H

