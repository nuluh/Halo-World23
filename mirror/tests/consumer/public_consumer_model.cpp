#include <openvfife/openvfife.hpp>

int openvfife_consumer_particle_id()
{
    openvfife::Particle particle(42, 1.0, 2.0, 3.0);
    return particle.id();
}
