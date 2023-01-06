// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include <cmath>
#include "sim.h"
#include "my_treadmilling_fiber_prop.h"
#include "my_treadmilling_fiber.h"
#include "property_list.h"
#include "simul_prop.h"
#include "exceptions.h"
#include "glossary.h"
#include "messages.h"
#include "simul.h"


Fiber* MyTreadmillingFiberProp::newFiber() const
{
    return new MyTreadmillingFiber(this);
}

void MyTreadmillingFiberProp::clear()
{
    FiberProp::clear();
    
    for ( int i = 0; i < 2; ++i )
    {
        growing_speed[i]	    = 0;
        growing_force[i]  	    = INFINITY;
        shrinking_speed[i] 	    = 0;
    }
}

void MyTreadmillingFiberProp::read(Glossary& glos)
{
    FiberProp::read(glos);
    
    glos.set(growing_speed,	2,	    "growing_speed");
    glos.set(growing_force,	2, 	    "growing_force");
    glos.set(shrinking_speed, 2, 	"shrinking_speed");

}


void MyTreadmillingFiberProp::complete(Simul const& sim)
{
    FiberProp::complete(sim);

    for ( int i = 0; i < 2; ++i )
    {
        if ( growing_force[i] <= 0 )
            throw InvalidParameter("fiber:growing_force should be > 0");
        if ( growing_speed[i] < 0 )
            throw InvalidParameter("fiber:growing_speed should be >= 0");
        
        growing_speed_dt[i]   = growing_speed[i] * sim.time_step();
//        growing_rate_dt[i] = sim.time_step() * fabs(growing_speed[i]) / unit_length;
        
        if ( shrinking_speed[i] > 0 )
            throw InvalidParameter("fiber:shrinking_speed should be <= 0");
        shrinking_speed_dt[i] = shrinking_speed[i] * sim.time_step();
//        shrinking_rate_dt[i] = sim.time_step() * fabs(shrinking_speed[i]) / unit_length;
    }
}


void MyTreadmillingFiberProp::write_values(std::ostream& os) const
{
    FiberProp::write_values(os);
    write_value(os, "unit_length",          unit_length);
    write_value(os, "growing_speed",	    growing_speed, 2);
    write_value(os, "growing_force",	    growing_force, 2);
    write_value(os, "shrinking_speed", 	    shrinking_speed, 2);
    write_value(os, "growing_speed_dt",     growing_speed_dt, 2);
    write_value(os, "shrinking_speed_dt",   shrinking_speed_dt, 2);
}

