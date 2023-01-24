#ifndef CCM_CCMRootHandle_H
#define CCM_CCMRootHandle_H

#include <list>
#include <memory>
#include <vector>
#include <stdint.h>
#include <boost/shared_ptr.hpp>

#include "TBranch.h"

class TFile;
class TTree;
class TBranch;

class CCMRootHandle {
protected:
    std::map<std::string, TBranch *> branches;
    std::map<std::string, size_t> num_entries;
    size_t tree_size;
    size_t current_index;
    TFile * input_file_ptr;

    bool tree_name_specified;
    std::string specified_tree_name;

    std::set<std::string> GetTreeNames();
    void LoadBranch(std::string const & branch_name);
    void LoadBranches(std::vector<std::string> const & branch_name);
    void LoadBranches(std::string const & tree_name, std::vector<std::string> const & branch_name);
    void LoadAllBranches();
    std::set<std::string> GetTreeNames();
public:
    CCMRootHandle(TFile * input_file_ptr);
    CCMRootHandle(TFile * input_file_ptr, std::vector<std::string> const & branch_names);
    CCMRootHandle(TFile * input_file_ptr, std::string const & tree_name, std::vector<std::string> const & branch_names);
    CCMRootHandle(TFile * input_file_ptr, std::string const & tree_name);

    template<typename T>
    void Get(std::string const & branch_name, boost::shared_ptr<T> data_ptr) {
        if(branches.count(branch_name) == 0)
            LoadBranch(tree_name)
        TBranch * branch_ptr = branches[branch_name];
        if(current_index >= num_entries[branch_name]) {
            data_ptr = nullptr;
            return;
        }
        branch_ptr->SetAddress(&(data_ptr.get()));
        branch_ptr->GetEntry(current_index);
    }

    template<typename T>
    boost::shared_ptr<T> Get(std::string const & branch_name) {
        boost::shared_ptr<T> data_ptr(nullptr);
        this->Get<T>(branch_name, data_ptr)
        return data_ptr;
    }

    size_t Advance(size_t n = 1);
    size_t Rewind(size_t n = 1);
    bool More();
    bool Valid();
};

//////////////////////////////////////////////////////////////
// Template functions must be inlines
//////////////////////////////////////////////////////////////

#endif // CCM_CCMRootHandle_H

