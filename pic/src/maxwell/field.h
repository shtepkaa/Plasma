#ifndef FIELD_H
#define FIELD_H

#include<types.h>

namespace maxwell {

template<Dimension dim, typename Type>
struct ElectroMagneticField
{
    Tuple<dim, Type> E;
    Tuple<dim, Type> B;
    
    Field(const Tuple<dim, Type> &, const Tuple<dim, Type> &);
    ~Field() {}
};

} // namespace maxwell

#endif // FIELD_H
