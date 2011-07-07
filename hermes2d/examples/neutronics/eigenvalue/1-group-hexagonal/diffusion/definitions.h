#include "hermes2d.h"
#include <sstream>

using namespace WeakFormsNeutronics::Multigroup;
using namespace CompleteWeakForms::Diffusion;
using namespace ElementaryForms::Diffusion;

class CustomWeakForm : public DefaultWeakFormSourceIteration
{
  public:
    CustomWeakForm(const MaterialPropertyMaps& matprop,
                   const Hermes::vector<Solution*>& iterates,
                   const Hermes::vector<std::string>& fission_regions,
                   double init_keff, std::string bdy_vacuum);
};

void report_num_dof(const std::string& msg, const Hermes::vector<Space *> spaces);
                    