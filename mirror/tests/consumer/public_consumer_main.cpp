#include <openvfife/openvfife.hpp>

int openvfife_consumer_particle_id();

int main()
{
    openvfife::LinearElastic material(1);
    material.setE(200.0);
    openvfife::Rectangle section(2, 1.0, 2.0);
    return (openvfife_consumer_particle_id() == 42 && material.E() == 200.0 && section.A() == 2.0) ? 0 : 1;
}
