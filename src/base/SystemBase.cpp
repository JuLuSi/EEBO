#include "SystemBase.h"

using namespace EEBO;
using namespace libMesh;

SystemBase::~SystemBase() {

}

SystemBase::SystemBase(libMesh::EquationSystems &eqs, const std::string &name, const unsigned int number) :
    TransientNonlinearImplicitSystem(eqs, name, number)
{

}
