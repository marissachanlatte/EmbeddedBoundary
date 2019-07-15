#ifndef NORMALS_H_
#define NORMALS_H_

#include "normals/NormalBase.h"

class Normal : public NormalBase{
  public:
    ~Normal() = default;
    int makeNormal();
};

#endif
