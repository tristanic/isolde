/**
 * @Author: Tristan Croll <tic20>
 * @Date:   11-Jun-2019
 * @Email:  tcroll@altoslabs.com
 * @Last modified by:   tic20
 * @Last modified time: 11-Jun-2019
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright: 2016-2019 Tristan Croll
 */



#pragma once

#include <map>
#include <string>

#include "template_atom_data.h"

namespace OpenMM_FF
{

class ResidueTemplate
{
public:
    ResidueTemplate(const std::string& name )
        : name_(name);

private:
    std::string name_;
    std::vector<TemplateAtomData> atoms_;



}; // class TemplateData

} // namespace OpenMM_FF
