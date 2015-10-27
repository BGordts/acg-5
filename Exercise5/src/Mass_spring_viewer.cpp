//=============================================================================
//
//   Exercise code for the lecture
//   "Advanced Computer Graphics"
//
//   Adapted from Prof. Dr. Mario Botsch, Bielefeld University
//
//   Copyright (C) 2013 LGG, epfl
//
//   DO NOT REDISTRIBUTE
//=============================================================================


//== INCLUDES =================================================================

#include "Mass_spring_viewer.h"
#include "utils/gl_text.h"
#include <sstream>
#include <cmath>


//== IMPLEMENTATION ==========================================================


Mass_spring_viewer::
Mass_spring_viewer(const char* _title, int _width, int _height)
: Viewer_2D(_title, _width, _height)
{
    integration_         = Euler;
    collisions_          = Force_based;
    external_force_      = None;
    animate_             = false;
    area_forces_         = true;
    show_forces_         = false;

    time_step_           = 0.001;
    particle_radius_     = 0.03;

    particle_mass_       = 0.1;
    damping_             = 0.1;
    collision_stiffness_ = 1000.0;
    spring_stiffness_    = 1000.0;
    spring_damping_      = 1.0;
    area_stiffness_      = 100000.0;

    mouse_spring_.active = false;
}


//-----------------------------------------------------------------------------


void Mass_spring_viewer::keyboard(int key, int x, int y)
{
    switch (key)
    {
        // setup problem 1
        case '1':
        {
            body_.clear();
            body_.add_particle( vec2(-0.5, -0.5), vec2(14.0, -2.0), particle_mass_, false );
            glutPostRedisplay();
            break;
        }

        // setup problem 5
        case '2':
        {
            body_.clear();
            for (int i=0; i<100; ++i)
            {
                body_.add_particle( vec2(0.9* cos(i/50.0*M_PI), 0.9*sin(i/50.0*M_PI)), vec2(-sin(i/50.0*M_PI), cos(i/50.0*M_PI)), particle_mass_, false );
            }

            glutPostRedisplay();
            break;
        }

        // setup problem 4
        case '3':
        {
            body_.clear();

            for (int i=0; i<10; ++i)
            {
                body_.add_particle( vec2(i*0.1, 0.8), vec2(0.0, 0.0), particle_mass_, i==0 );
            }

            for (unsigned int i=0; i<9; ++i)
            {
                body_.add_spring(i, i+1);
            }

            glutPostRedisplay();
            break;
        }

        // setup problem 2
        case '4':
        {
            body_.clear();

            body_.add_particle( vec2(-0.1, 0.7), vec2(0.0, 0.0), particle_mass_, false );
            body_.add_particle( vec2( 0.0, 0.6), vec2(0.0, 0.0), particle_mass_, false );
            body_.add_particle( vec2( 0.1, 0.7), vec2(0.0, 0.0), particle_mass_, false );

            body_.add_spring(0, 1);
            body_.add_spring(0, 2);
            body_.add_spring(1, 2);

            body_.add_triangle(0, 1, 2);

            glutPostRedisplay();
            break;
        }

        // setup problem 3
        case '5':
        {
            body_.clear();

            for (int i=0; i<8; ++i)
            {
                body_.add_particle( vec2(-0.5+0.2*cos(0.25*i*M_PI), -0.5+0.2*sin(0.25*i*M_PI)), vec2(5.0, 5.0), particle_mass_, false );
            }

            body_.add_particle( vec2(-0.5, -0.5), vec2(5.0, 5.0), particle_mass_, false );

            for (unsigned int i=0; i<8; ++i)
            {
                body_.add_spring( i, (i+1) % 8 );
                body_.add_spring(i, 8);
                body_.add_triangle(i, (i+1)%8, 8);
            }

            glutPostRedisplay();
            break;
        }

        // switch between time integration methods
        case 'i':
        {
            switch (integration_)
            {
                case Euler:    integration_ = Midpoint; break;
                case Midpoint: integration_ = Verlet;   break;
                case Verlet:   integration_ = Euler;    break;
            }
            glutPostRedisplay();
            break;
        }

        // switch between center and gravitation force
        case 'f':
        {
            switch (external_force_)
            {
                case None:        external_force_ = Center; break;
                case Center:      external_force_ = Gravitation; break;
                case Gravitation: external_force_ = None;      break;
            }
            glutPostRedisplay();
            break;
        }

        // switch between force-based and impulse-based collisions
        case 'c':
        {
            switch (collisions_)
            {
                case Force_based:   collisions_ = Impulse_based; break;
                case Impulse_based: collisions_ = Force_based;   break;
            }
            glutPostRedisplay();
            break;
        }


        // toggle area forces on/off
        case 'a':
        {
            area_forces_ = !area_forces_;
            glutPostRedisplay();
            break;
        }


        // visualization of particle forces on/off
        case 'v':
        {
            show_forces_ = !show_forces_;
            glutPostRedisplay();
            break;
        }


        // let parent class do the work
        default:
        {
            Viewer_2D::keyboard(key, x, y);
            break;
        }
    }
}


//-----------------------------------------------------------------------------


void Mass_spring_viewer::draw()
{
    // parent's status text
    Viewer_2D::draw();


    // draw some status text
    glDisable(GL_LIGHTING);
    glColor3f(1,0,0);
    std::ostringstream oss;

    oss.str("");
    oss << "Integration: ";
    switch (integration_)
    {
        case Euler:    oss << "Euler";    break;
        case Midpoint: oss << "Midpoint"; break;
        case Verlet:   oss << "Verlet";   break;
    }
    glText(20, height_-40, oss.str());

    oss.str("");
    oss << "#Particles: " << body_.particles.size();
    glText(20, height_-60, oss.str());

    oss.str("");
    oss << "#Springs: " << body_.springs.size();
    glText(20, height_-80, oss.str());

    oss.str("");
    oss << "#Triangles: " << body_.triangles.size();
    glText(20, height_-100, oss.str());

    oss.str("");
    oss << "Area Forces: " << (area_forces_ ? "on" : "off");
    glText(20, height_-120, oss.str());

    oss.str("");
    oss << "Collisions: " << (collisions_ == Force_based ? "force" : "impulse");
    glText(20, height_-140, oss.str());

    oss.str("");
    oss << "External force: ";
    switch (external_force_)
    {
        case None:        oss << "None";        break;
        case Center:      oss << "Center";      break;
        case Gravitation: oss << "Gravitation"; break;
    }
    glText(20, height_-160, oss.str());

    oss.str("");
    oss << "Visualize forces: " << (show_forces_ ? "on" : "off");
    glText(20, height_-180, oss.str());


    // draw walls
    glDisable(GL_LIGHTING);
    glLineWidth(1.0);
    glColor3f(0.5,0.5,0.5);
    glBegin(GL_LINE_STRIP);
    glVertex2f( -1.0,  1.0 );
    glVertex2f( -1.0, -1.0 );
    glVertex2f(  1.0, -1.0 );
    glVertex2f(  1.0,  1.0 );
    glVertex2f( -1.0,  1.0 );
    glEnd();


    // draw mouse spring
    if (mouse_spring_.active)
    {
        glDisable(GL_LIGHTING);
        glLineWidth(5.0);
        glColor3f(1,0,0);
        glBegin(GL_LINES);
        glVertex2fv( body_.particles[mouse_spring_.particle_index].position.data() );
        glVertex2fv( mouse_spring_.mouse_position.data() );
        glEnd();
    }


    // draw particles, springs, triangles
    body_.draw(particle_radius_, show_forces_);
}


//-----------------------------------------------------------------------------


void Mass_spring_viewer::mouse(int _button, int _state, int _x, int _y)
{
    // need particles to do interaction
    if (!body_.particles.empty())
    {
        // mouse button release destroys current mouse spring
        if (_state == GLUT_UP)
        {
            mouse_spring_.active = false;
        }

        // mouse button press generates new mouse spring
        else if (_state == GLUT_DOWN)
        {
            // get point under mouse cursor
            vec2 p = pick(_x, _y);

            // find closest particle
            int   pidx = -1;
            float dmin = FLT_MAX;
            for (unsigned int i=0; i<body_.particles.size(); ++i)
            {
                float d = norm(p - body_.particles[i].position);
                if (d < dmin)
                {
                    dmin = d;
                    pidx = i;
                }
            }

            // construct mouse spring
            mouse_spring_.mouse_position = p;
            mouse_spring_.particle_index = pidx;
            mouse_spring_.active = true;
        }
    }

    glutPostRedisplay();
}


//-----------------------------------------------------------------------------


void Mass_spring_viewer::motion(int _x, int _y)
{
    if (mouse_spring_.active)
    {
        // update mouse positions
        mouse_spring_.mouse_position = pick(_x, _y);
    }

    glutPostRedisplay();
}


//-----------------------------------------------------------------------------


void Mass_spring_viewer::time_integration(float dt)
{
    switch (integration_)
    {
        case Euler:
        {
            /** \todo (Part 1) Implement Euler integration scheme
             \li Hint: compute_forces() computes all forces for the current positions and velocities.
             */

            compute_forces();

            for (unsigned int i=0; i<body_.particles.size(); ++i)
            {
                body_.particles[i].acceleration = body_.particles[i].force / body_.particles[i].mass;
                body_.particles[i].velocity += dt*body_.particles[i].acceleration;
                body_.particles[i].position += dt*body_.particles[i].velocity; // Use new velocity
            }

            break;
        }

        case Midpoint:
        {
            /** \todo (Part 2) Implement the Midpoint time integration scheme
             \li The Particle class has variables position_t and velocity_t to store current values
             \li Hint: compute_forces() computes all forces for the current positions and velocities.
             */

            break;
        }


        case Verlet:
        {
            /** \todo (Part 2) Implement the Verlet time integration scheme
             \li The Particle class has a variable acceleration to remember the previous values
             \li Hint: compute_forces() computes all forces for the current positions and velocities.
             */

            break;
        }
    }


    // impulse-based collision handling
    if (collisions_ == Impulse_based)
    {
        impulse_based_collisions();
    }


    glutPostRedisplay();
}


//-----------------------------------------------------------------------------


void
Mass_spring_viewer::compute_forces()
{
    // clear forces
    for (unsigned int i=0; i<body_.particles.size(); ++i)
        body_.particles[i].force = vec2(0,0);


    /** \todo (Part 1) Implement center force
     */
    if (external_force_ == Center)
    {
        for (unsigned int i=0; i<body_.particles.size(); ++i)
        {
            body_.particles[i].force += 20*(vec2(0,0) - body_.particles[i].position);
        }
    }


    /** \todo (Part 1) Implement damping force
     \li The damping coefficient is given as member damping_
     */
    for (unsigned int i=0; i<body_.particles.size(); ++i)
    {
        body_.particles[i].force += -damping_*body_.particles[i].velocity;
    }



    /** \todo (Part 1) Implement gravitation force
     \li Particle mass available as particle_mass_
     */
    if (external_force_ == Gravitation)
    {
        for (unsigned int i=0; i<body_.particles.size(); ++i)
        {
            body_.particles[i].force += vec2(0.f,-9.81)*body_.particles[i].mass;
            // body_.particles[i].force += vec2(0.f,-9.81)*particle_mass_;
        }
    }


    /** \todo (Part 1) Implement force based boundary collisions
     \li Collision coefficient given as collision_stiffness_
     */
    // collision forces
    if (collisions_ == Force_based)
    {
        const int nbPanes = 4;
        float planes[nbPanes][3] = {
            {  0.0,  1.0, 1.0 },
            {  0.0, -1.0, 1.0 },
            {  1.0,  0.0, 1.0 },
            { -1.0,  0.0, 1.0 }/*,
            {  -0.2, 1.0, 0.9 }*/ // Test with an oblique pane (Warning: change "nbPanes" to 5)
        };

        vec2 planesNorms[nbPanes]; // Norms of the plane
        float planesP[nbPanes]; // Distance from the origin

        for(unsigned int i=0; i<nbPanes; ++i)
        {
            // According to http://mathworld.wolfram.com/HessianNormalForm.html

            float sqrtCoef = std::sqrt(planes[i][0]*planes[i][0] +
                                       planes[i][1]*planes[i][1] +
                                       planes[i][2]*planes[i][2]);

            planesNorms[i] = vec2(planes[i][0]/sqrtCoef,
                                  planes[i][1]/sqrtCoef);

            planesP[i] = planes[i][2]/sqrtCoef;
        }

        for (unsigned int i=0; i<body_.particles.size(); ++i)
        {
            vec2 boudaryVec(0,0);
            bool hasCollisions = false;

            for (unsigned int j=0 ; j<nbPanes ; ++j)
            {
                float distPointPlanes = dot(planesNorms[j], body_.particles[i].position) + planesP[j];
                if(distPointPlanes < particle_radius_) // Collision
                {
                    // std::cout << distPointPlanes << std::endl;
                    boudaryVec += collision_stiffness_ * (particle_radius_-distPointPlanes) * planesNorms[j];
                    hasCollisions = true;
                }
            }

            if(hasCollisions)
            {
                body_.particles[i].force += boudaryVec;
            }
        }

    }


    /** \todo (Part 1) Compute force of the interactive mouse spring
     \li Required coefficients are given as spring_stiffness_ and spring_damping_
     */
    if (mouse_spring_.active)
    {
        Particle& p0 = body_.particles[ mouse_spring_.particle_index ];

        vec2 pos0 = p0.position;
        vec2 pos1 = mouse_spring_.mouse_position;

        vec2 vel0 = p0.velocity;
        vec2 vel1 = vec2(0,0);

        float d = norm(pos0 - pos1);

        p0.force += -1 * (spring_stiffness_*d + spring_damping_*dot(vel0-vel1, pos0-pos1)/d) * ((pos0-pos1)/d);

    }


    /** \todo (Part 1) Compute spring forces
     \li Required information about springs in the scene are found in body_.springs
     \li Required coefficients are given as spring_stiffness_ and spring_damping_
     */
    for (unsigned int i=0; i<body_.springs.size(); ++i){
        Spring &nextSpring = body_.springs[i];

        vec2 pos0 = nextSpring.particle0->position;
        vec2 pos1 = nextSpring.particle1->position;

        vec2 vel0 = nextSpring.particle0->velocity;
        vec2 vel1 = nextSpring.particle1->velocity;

        float d = norm(pos0 - pos1);

        vec2 f0 = -1 * (spring_stiffness_*(d-nextSpring.rest_length) + spring_damping_*dot(vel0-vel1, pos0-pos1)/d) * ((pos0-pos1)/d);
        nextSpring.particle0->force += f0;
        nextSpring.particle1->force -= f0;

    }


    /** \todo (Part 2) Compute more forces in part 2 of the exercise: triangle-area forces, binding forces, etc.
     */
}


//-----------------------------------------------------------------------------


void Mass_spring_viewer::impulse_based_collisions()
{
    /** \todo (Part 2) Handle collisions based on impulses
     */
    // planes for which we compute collisions
    float planes[4][3] = {
        {  0.0,  1.0, 1.0 },
        {  0.0, -1.0, 1.0 },
        {  1.0,  0.0, 1.0 },
        { -1.0,  0.0, 1.0 }
    };
}


//=============================================================================
