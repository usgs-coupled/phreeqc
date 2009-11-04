#include "Phreeqc.h"
#include <algorithm>			// std::replace

Phreeqc::Phreeqc(void)
{
}
Phreeqc::Phreeqc(CParser & parser)
{
  this->clean_up();
}

Phreeqc::~Phreeqc(void)
{
}
