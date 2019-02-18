#pragma once

#include <map>
#include <string>

#include "template_atom_data.h"

namespace OpenMM_FF
{

class TemplateData
{
public:
    TemplateData(const std::string& name )
        : name_(name);

private:
    std::string name_;
    std::vector<TemplateAtomData> atoms_;
    


}; // class TemplateData

} // namespace OpenMM_FF
