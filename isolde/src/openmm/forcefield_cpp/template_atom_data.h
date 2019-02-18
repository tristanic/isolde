#pragma once

#include <map>
#include <string>

namespace OpenMM_FF {

class TemplateAtomData
{
public:
    TemplateAtomData(const std::string& name, const std::string& type,
        const std::string& element, const std::map<std::string, double> params)
        : name_(name), type_(type), element_(element), params_(params) {}

    inline const std::string& name() const { return name_; }
    inline const std::string& element() const { return element_; }
    inline const std::string& type() const { return type_; }
    inline double get_param(const std::string& name) const { return params_.at(name); }
    inline const std::vector<int>& bondedTo() const { return bonded_to_; }
    inline addBond(int index) { external_bonds_.push_back(index); }
    inline const int& externalBonds() const { return external_bonds_; }
    inline int& externalBonds() { return external_bonds_; }

private:
    std::string name_;
    std::string type_;
    std::string element_;
    std::map<std::string, double> params_;
    std::vector<int> bonded_to_;
    int external_bonds_ = 0;

}; // class TemplateAtomData

} // namespace OpenMM_FF
