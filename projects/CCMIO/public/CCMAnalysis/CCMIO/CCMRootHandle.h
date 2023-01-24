#ifndef CCM_CCMRootHandle_H
#define CCM_CCMRootHandle_H

#include <map>
#include <set>
#include <list>
#include <memory>
#include <vector>
#include <stdint.h>
#include <boost/shared_ptr.hpp>

#include <icetray/I3FrameObject.h>
#include <icetray/serialization.h>

#include <TBranch.h>

class TFile;

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
        if(branches.count(branch_name) == 0)
            LoadBranch(branch_name);
        TBranch * branch_ptr = branches[branch_name];
        if(branch_ptr == nullptr or current_index >= num_entries[branch_name]) {
            data_ptr = nullptr;
            return;
        }
        T * ptr_destination = data_ptr.get();
        bool new_allocation = ptr_destination == nullptr;
        if(new_allocation)
            ptr_destination = new T();
        branch_ptr->SetAddress(&ptr_destination);
        branch_ptr->GetEntry(current_index);

        if(new_allocation)
            data_ptr = boost::shared_ptr<T>(ptr_destination);
    }

    template<typename T>
    boost::shared_ptr<T> Get(std::string const & branch_name) {
        boost::shared_ptr<T> data_ptr(nullptr);
        this->Get<T>(branch_name, data_ptr);
        return data_ptr;
    }

    void Advance(size_t n = 1);
    void Rewind(size_t n = 1);
    bool More();
    bool Valid();
};

I3_POINTER_TYPEDEFS(CCMRootHandle);

#endif // CCM_CCMRootHandle_H

