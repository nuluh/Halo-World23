/* ************************************************************************* **
**    OpenVFIFE - Open System for Vector Form Instrinsic                     **
**                Finite Element Method (VFIFE)                              **
**                GinkGo(Tan Biao)                                           **
**                                                                           **
**                                                                           **
** (C) Copyright 2021, The GinkGo(Tan Biao). All Rights Reserved.            **
**                                                                           **
** Commercial use of this program without express permission of              **
** GinkGo(Tan Biao), is strictly prohibited.  See                            **
** file 'COPYRIGHT'  in main directory for information on usage and          **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.                  **
**                                                                           **
** Developed by:                                                             **
**      Tan Biao (ginkgoltd@outlook.com)                                     **
**                                                                           **
** ************************************************************************* */

// $Date: 2020-05-11 $
// Written: Tan Biao
// Revised:
//
// Purpose: This file contains the class definition for StructSystem

// The interface:
//

#include <openvfife/solver/struct_system.hpp>
#include "filesystem_helpers.hpp"

#include <filesystem>
#include <stdexcept>
#include <system_error>
#include <fstream>
#include <iostream>
using namespace std;

namespace openvfife {


/*****************************************************************************/
/*                            other functions                                */
/*****************************************************************************/
bool isFileExists(const string &fname)
{
    std::error_code ec;
    return std::filesystem::is_regular_file(fname, ec);
}

bool isDirExists(const string &path)
{
    std::error_code ec;
    return std::filesystem::is_directory(path, ec);
}

bool makedir(const string &path)
{
    std::error_code ec;
    if (std::filesystem::create_directories(path, ec))
    {
        return true;
    }
    return !ec && std::filesystem::is_directory(path);
}


/*****************************************************************************/
/*                            helper functions                               */
/*****************************************************************************/
inline void clearCsv(const string &fname)
{
    if (isFileExists(fname))
    {
        fstream fout;
        fout.open(fname, ios::out | ios::trunc);
        fout.close();
    }
}

inline void writeData(fstream &fout, const int ID, const StdArray6d& arr)
{
    fout << ID << ", ";
    for (int i = 0; i < 5; i++)
    {
        fout << arr[i] << ", ";
    }
    fout << arr[5] << endl;
}

/*****************************************************************************/
/*                               StructSystem                                */
/*****************************************************************************/
StructSystem::StructSystem(int id)
{
    cout << "########## StructSystem Established! ##########" << endl;
    id_ = id;
    h_ = 1e-3;
    alpha_ = 0.2;
    beta_ = 0.0;
    critical_time_step_ = 1e-2;
    workdir_ = std::filesystem::current_path().string();
    jobname_ = "filename";
    elem_max_nodes_ = 0;
    elem_types_.push_back("Link");
    elem_types_.push_back("Beam");
    elem_types_.push_back("Plane");
    elem_types_.push_back("Solid");
}

StructSystem::~StructSystem()
{
    clearOwnedObjects();
}

string StructSystem::jobname() const
{
    return jobname_;
}

string StructSystem::workdir() const
{
    return workdir_;
}

double StructSystem::timestep() const
{
    return h_;
}

double StructSystem::critical_time_step() const
{
    return critical_time_step_;
}

const StructSystem::ParticleMap& StructSystem::particles() const noexcept
{
    return particles_;
}

const StructSystem::MaterialMap& StructSystem::materials() const noexcept
{
    return materials_;
}

const StructSystem::SectionMap& StructSystem::sections() const noexcept
{
    return sections_;
}

const StructSystem::ElementMap& StructSystem::elements() const noexcept
{
    return elements_;
}

const map<int, DOF>& StructSystem::constraints() const noexcept
{
    return constraints_;
}

void StructSystem::info() const
{
    for (int i = 0; i < 80; i++)
    {
        cout << "#";
    }
    cout << endl << "#";
    for (int i = 0; i < 30; i++)
    {
        cout << " ";
    }
    cout << "System Information";
    for (int i = 0; i < 30; i++)
    {
        cout << " ";
    }
    cout << endl;
    for (int i = 0; i < 80; i++)
    {
        cout << "#";
    }
    cout << endl;
    cout << "1. Working Directory: " << workdir_ << endl;
    cout << "2. Job Name: " << jobname_ << endl;
    cout << "3. Number of Particles: " << particles_.size() << endl;
    cout << "4. Number of Elements: " << elements_.size() << endl;
    cout << endl;
    // print paritcle details
    for (int i = 0; i < 31; i++)
    {
        cout << "#";
    }
    cout << " Particle Details ";
    for (int i = 0; i < 31; i++)
    {
        cout << "#";
    }
    cout << endl;
    for (auto iter = particles_.begin(); iter != particles_.end(); ++iter)
    {
        iter->second->info();
    }
    // print element details
    cout << endl;
    for (int i = 0; i < 31; i++)
    {
        cout << "#";
    }
    cout << " Element Details ";
    for (int i = 0; i < 32; i++)
    {
        cout << "#";
    }
    cout << endl;
    for (auto iter = elements_.begin(); iter != elements_.end(); ++iter)
    {
        iter->second->info();
    }
    // print conststaint details
    cout << endl;
    for (int i = 0; i < 30; i++)
    {
        cout << "#";
    }
    cout << " Constraint Details ";
    for (int i = 0; i < 30; i++)
    {
        cout << "#";
    }
    cout << endl;
    cout << "PID |    Ux    |    Uy    |    Uz    |" <<
                 "   Rotx   |   Roty   |   Rotz   |" << endl;
    cout << "    | key, val | key, val | key, val |" <<
                 " key, val | key, val | key, val |" << endl;
    for (auto iter = constraints_.begin(); iter != constraints_.end(); ++iter)
    {
        cout << iter->first << " | ";
        for (int i = 0; i < 6; i++)
        {
            cout << iter->second.key[i] << ", " << iter->second.val[i] << " | ";
        }
        cout << endl;
    }
    cout << endl;
}

void StructSystem::setJobname(const string &job)
{
    jobname_ = job;
}

void StructSystem::setWorkdir(const string &dir)
{
    workdir_ = dir;
}

void StructSystem::setTimeStep(double h)
{
    if (!(h > 0.0))
    {
        throw std::invalid_argument("time step must be positive");
    }
    h_ = h;
}

void StructSystem::setDampCoeff(double alpha, double beta)
{
    if (alpha < 0.0 || beta < 0.0)
    {
        throw std::invalid_argument("damping coefficients must be non-negative");
    }
    alpha_ = alpha;
    beta_ = beta;
}

void StructSystem::addParticle(std::unique_ptr<Particle> p)
{
    if (!p)
    {
        throw std::invalid_argument("particle ownership transfer requires a non-null object");
    }
    const int id = p->id();
    if (particles_.contains(id))
    {
        throw std::invalid_argument("duplicate particle id");
    }
    particles_.emplace(id, std::move(p));
}

void StructSystem::addParticles(std::vector<std::unique_ptr<Particle>> particles)
{
    for (auto& particle : particles)
    {
        addParticle(std::move(particle));
    }
}

void StructSystem::addElement(std::unique_ptr<BaseElement> e)
{
    if (!e)
    {
        throw std::invalid_argument("element ownership transfer requires a non-null object");
    }
    const int id = e->id();
    if (elements_.contains(id))
    {
        throw std::invalid_argument("duplicate element id");
    }
    BaseElement* raw = e.get();
    elements_.emplace(id, std::move(e));
    // update critical time step
    if (critical_time_step_ > raw->calcCriticalTimeStep())
    {
        critical_time_step_ = raw->calcCriticalTimeStep();
    }
    // update max element nodes
    if (elem_max_nodes_ <= raw->num_of_particles())
    {
        elem_max_nodes_ = raw->num_of_particles();
    }
    // update particle_connected_elems_
    const auto& ps = raw->particles();
    for (auto p = ps.begin(); p != ps.end(); p++)
    {
       if (particle_connected_elems_.count(p->first))
       {
           particle_connected_elems_[p->first] += 1;
       }
       else
       {
           particle_connected_elems_[p->first] = 1;
       }
    }
}

void StructSystem::addElements(std::vector<std::unique_ptr<BaseElement>> elements)
{
    for (auto& element : elements)
    {
        addElement(std::move(element));
    }
}

void StructSystem::addMaterial(std::unique_ptr<BaseMaterial> mat)
{
    if (!mat)
    {
        throw std::invalid_argument("material ownership transfer requires a non-null object");
    }
    const int id = mat->id();
    if (materials_.contains(id))
    {
        throw std::invalid_argument("duplicate material id");
    }
    materials_.emplace(id, std::move(mat));
}

void StructSystem::addMaterials(std::vector<std::unique_ptr<BaseMaterial>> materials)
{
    for (auto& material : materials)
    {
        addMaterial(std::move(material));
    }
}

void StructSystem::addSection(std::unique_ptr<BaseSection> sect)
{
    if (!sect)
    {
        throw std::invalid_argument("section ownership transfer requires a non-null object");
    }
    const int id = sect->id();
    if (sections_.contains(id))
    {
        throw std::invalid_argument("duplicate section id");
    }
    sections_.emplace(id, std::move(sect));
}

void StructSystem::addSections(std::vector<std::unique_ptr<BaseSection>> sections)
{
    for (auto& section : sections)
    {
        addSection(std::move(section));
    }
}

void StructSystem::addConstraint(int pid, const DOF &dof)
{
    if (!particles_.contains(pid))
    {
        throw std::out_of_range("cannot constrain a particle that is not in the system");
    }
    if (constraints_.contains(pid))
    {
        throw std::invalid_argument("particle is already constrained");
    }
    constraints_.emplace(pid, dof);
}

bool StructSystem::deleteConstraint(int pid)
{
    return constraints_.erase(pid) != 0;
}

void StructSystem::changeConstraint(int pid, const DOF &dof)
{
    if (!particles_.contains(pid))
    {
        throw std::out_of_range("cannot change a constraint for a missing particle");
    }
    auto it = constraints_.find(pid);
    if (it == constraints_.end())
    {
        throw std::out_of_range("particle does not have an existing constraint");
    }
    it->second = dof;
}

void StructSystem::clearConstraint()
{
    constraints_.clear();
}

void StructSystem::constraint()
{
    for (auto iter = constraints_.begin(); iter != constraints_.end(); iter++)
    {
        int pid = iter->first;
        particles_[pid]->constraintDof(iter->second);
    }
}

void StructSystem::setParticleDisplace(int pid, const StdArray6d &val)
{
    particles_[pid]->setDisplace(val);
}

void StructSystem::addParticleDisplace(int pid, const StdArray6d &val)
{
    particles_[pid]->addDisplace(val);
}

void StructSystem::setAccelerate(double ax, double ay, double az)
{
    StdArray6d acc = {ax, ay, az, 0, 0, 0};
    StdArray6d force = {0, 0, 0, 0, 0, 0};
    for (auto iter = particles_.begin(); iter != particles_.end(); iter++)
    {
        const Mass& m = iter->second->mass();
        for (int i = 0; i < 6; i++)
        {
            // Caution: the orient of force and acclerate are opposite
            force[i] = -m.lumped_mass[i] * acc[i];
        }
        iter->second->addForce(force);
    }
}

void StructSystem::setExternalForce(int pid, const StdArray6d &val)
{
    if (particles_.count(pid))
    {
        particles_[pid]->setForce(val);
    }
    else
    {
        cerr << "WARNING: Particle " << pid << " not found!" << endl;
    }
}

void StructSystem::addExternalForce(int pid, const StdArray6d &val)
{
    if (particles_.count(pid))
    {
        particles_[pid]->addForce(val);
    }
    else
    {
        cerr << "WARNING: Particle " << pid << " not found!" << endl;
    }
}

void StructSystem::setInternalForce()
{
    // loop for elements
    for (auto iter = elements_.begin(); iter != elements_.end(); iter++)
    {
        iter->second->setParticleForce();
    }
}

void StructSystem::clearParticleForce()
{
    // clear all force
    for (auto iter = particles_.begin(); iter != particles_.end(); iter++)
    {
        iter->second->clearForce();
    }
}

double StructSystem::autoTimeStep()
{
    h_ = pow(10, floor(log10(critical_time_step_)));
    return h_;
}

void StructSystem::solve(double h, double zeta, bool firstStep)
{
    if (!(h > 0.0))
    {
        throw std::invalid_argument("time step must be positive");
    }
    // loop for particles
    for (auto iter = particles_.begin(); iter != particles_.end(); iter++)
    {
        // cdm(iter->second, h, zeta, firstStep);
        if (firstStep)
        {
            iter->second->solveFirstStep(h, zeta);
        }
        else
        {
            iter->second->solve(h, zeta);
        }
    }
    constraint();
}

void StructSystem::solve(bool firstStep)
{
    // loop for particles
    for (auto iter = particles_.begin(); iter != particles_.end(); iter++)
    {
        // cdm(iter->second, h, zeta, firstStep);
        if (firstStep)
        {
            iter->second->solveFirstStep(h_, alpha_);
        }
        else
        {
            iter->second->solve(h_, alpha_);
        }
    }
    constraint();
}

void StructSystem::saveModel(const string &path)
{
    // if directory does not exist, create it
    if (!isDirExists(path))
    {
        makedir(path);
    }

    // save particles information
    {
        string fname = path + "/coordinates.csv";
        clearCsv(fname);
        fstream fout;
        fout.open(fname, ios::out | ios::app);
        fout << "PID, x, y, z, thx, thy, thz" << endl;
        for (auto iter = particles_.begin(); iter != particles_.end(); iter++)
        {
            writeData(fout, iter->first, iter->second->coordinate());
        }
        fout.close();
    }

    // save element information
    {
        string fname = path + "/elements.csv";
        clearCsv(fname);
        fstream fout;
        fout.open(fname, ios::out | ios::app);
        fout << "EID, Type, ";
        for (int i = 1; i < elem_max_nodes_; i++)
        {
            fout << "p" << i << ", ";
        }
        fout << "p" << elem_max_nodes_ << endl;
        for (auto iter = elements_.begin(); iter != elements_.end(); iter++)
        {
            const auto& ps = iter->second->particles();
            fout << iter->first << ", " << iter->second->type() << ", ";
            for (auto p = ps.begin(); p != ps.end(); p++)
            {
                fout << p->first << ", ";
            }
            fout << "\n";
        }
        fout.close();
    }

    // save constraint information
    {
        string fname = path + "/constraints.csv";
        clearCsv(fname);
        fstream fout;
        fout.open(fname, ios::out | ios::app);
        fout << "PID, Ux, Uy, Uz, Rotx, Roty, Rotz" << endl;

        for (auto dof = constraints_.begin(); dof != constraints_.end(); dof++)
        {
            StdArray6d array {};
            for (int i = 0; i < 6; i++)
            {
                if (dof->second.key[i])
                {
                    array[i] = dof->second.val[i];
                }
                else
                {
                    array[i] = 0.0 / 0.0;
                }
            }
            writeData(fout, dof->first, array);
        }
        fout.close();
    }
}

void StructSystem::saveParticleResult(const string &path)
{
    // if directory does not exist, create it
    if (!isDirExists(path))
    {
        makedir(path);
    }

    // create displace.csv
    vector<string> items = {"displace", "velocity", "accelerate"};
    for (int i = 0; i < 3; i++)
    {
        // create file
        string fname = path + "/" + items[i] + ".csv";
        clearCsv(fname);
        fstream fout;
        fout.open(fname, ios::out | ios::app);
        fout << "PID, Ux, Uy, Uz, Rotx, Roty, Rotz" << endl;

        // write data
        for (auto iter = particles_.begin(); iter != particles_.end(); iter++)
        {
            if (i == 0)
            {
                iter->second->outputDisplace(fout);
            }
            else if (i == 1)
            {
                iter->second->outputVelocity(fout);
            }
            else
            {
                iter->second->outputAccelarate(fout);
            }
        }

        // close file
        fout.close();
    }
}

void StructSystem::saveElementResult(const std::string &path)
{
    // if directory does not exist, create it
    if (!isDirExists(path))
    {
        makedir(path);
    }

    // create link_element_force.csv
    string link_fname = path + "/link_element_force.csv";
    clearCsv(link_fname);
    fstream link_out;
    link_out.open(link_fname, ios::out | ios::app);
    string link_header = "EID, N";
    link_out << link_header << endl;

    // create cable_element_force.csv
    string cable_fname = path + "/cable_element_force.csv";
    clearCsv(cable_fname);
    fstream cable_out;
    cable_out.open(cable_fname, ios::out | ios::app);
    string cable_header = "EID, N";
    cable_out << cable_header << endl;

    // create beam_element_force.csv
    string beam_fname = path + "/beam_element_force.csv";
    clearCsv(beam_fname);
    fstream beam_out;
    beam_out.open(beam_fname, ios::out | ios::app);
    string beam_header = "EID, PID, Fx, SFx, SFy, Mx, My, Mz";
    beam_out << beam_header << endl;

    // write data
    for (auto iter = elements_.begin(); iter != elements_.end(); iter++)
    {
        if (iter->second->type().find("Link") != string::npos)
        {
            iter->second->outputForce(link_out);
        }
        else if (iter->second->type().find("Beam") != string::npos)
        {
            iter->second->outputForce(beam_out);
        }
        else if (iter->second->type().find("Cable") != string::npos)
        {
            iter->second->outputForce(cable_out);
        }
    }

    // close file
    link_out.close();
    beam_out.close();
    cable_out.close();
}

void StructSystem::saveSupportReact(const std::string &path)
{
    // if directory does not exist, create it
    if (!isDirExists(path))
    {
        makedir(path);
    }

    // create file
    string fname = path + "/" + "support_reaction.csv";
    clearCsv(fname);
    fstream fout;
    fout.open(fname, ios::out | ios::app);

    // write header
    string header = "PID, Fx, Fy, Fz, Mx, My, Mz";
    fout << header << endl;

    // write data
    for (auto iter = constraints_.begin(); iter != constraints_.end(); iter++)
    {
        particles_[iter->first]->outputReactionForce(fout);
    }

    // close file
    fout.close();
}

void StructSystem::saveResult(const string &path)
{
    saveModel(path);
    saveParticleResult(path);
    saveElementResult(path);
    saveSupportReact(path);
}

void StructSystem::clearOwnedObjects() noexcept
{
    // Elements observe particles/materials/sections, so destroy them first.
    failure_elems_.clear();
    elements_.clear();
    free_particles_.clear();
    particles_.clear();
    sections_.clear();
    materials_.clear();
    particle_connected_elems_.clear();
}


void StructSystem::killElement(const int id)
{
    if (!elements_.count(id))
    {
        cout << "Warning: Element " << id << " not exists!" << endl;
    }
    else
    {
        // change particle_connected_elems_
        const auto& ps = elements_[id]->particles();
        for (auto p = ps.begin(); p != ps.end(); p++)
        {
            if (particle_connected_elems_.count(p->first))
            {
                particle_connected_elems_[p->first] -= 1;
            }
            else
            {
                cout << "Warning: Particle " << p->first << " not found!\n";
            }
        }
        failure_elems_[id] = std::move(elements_.at(id));
        elements_.erase(id);
    }
}

void StructSystem::autoKillElement()
{
    // WARNING: erase data in a loop might cause collapse of program
    for (auto iter = elements_.begin(); iter != elements_.end();)
    {
        if (iter->second->is_failed())
        {
            // change particle_connected_elems_
            const auto& ps = iter->second->particles();
            for (auto p = ps.begin(); p != ps.end(); p++)
            {
                if (particle_connected_elems_.count(p->first))
                {
                    particle_connected_elems_[p->first] -= 1;
                }
                else
                {
                    cout << "Warning: Particle " << p->first << " not found!\n";
                }
            }
            failure_elems_[iter->first] = std::move(iter->second);
            iter = elements_.erase(iter);
        }
        else
        {
            iter++;
        }
    }
}

void StructSystem::removeFreeParticle()
{
    for (auto iter = particle_connected_elems_.begin();
        iter != particle_connected_elems_.end(); iter++)
    {
        if ((particle_connected_elems_[iter->first] <= 0) &&
            (!constraints_.count(iter->first)))
        {
            free_particles_[iter->first] = std::move(particles_.at(iter->first));
            particles_.erase(iter->first);
        }
    }
}

} // namespace openvfife
