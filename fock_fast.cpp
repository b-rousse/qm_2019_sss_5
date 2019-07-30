#include <Eigen/Dense>
#include <map>
#include <string>
#include <vector>
#include <iostream>

int atom(int ao_index, int orbitals_per_atom);
std::string orb(int ao_index, std::vector<std::string> orbital_types, int orbitals_per_atom);
int ao_index(int atom_p, std::string orb_p, std::vector<std::string> orbital_types, int orbitals_per_atom);
double chi_on_atom(std::string orb1, std::string orb2, std::string orb3, std::map<std::string, double> model_parameters);
Eigen::MatrixXd calc_fock_matrix_fast(Eigen::MatrixXd hamiltonian_matrix, Eigen::MatrixXd interaction_matrix, Eigen::MatrixXd density_matrix, std::map<std::string, double> model_parameters, int orbitals_per_atom);

int main()
{
    std::map<std::string, double> model_parameters;
    model_parameters["r_hop"] = 3.1810226927827516;
    model_parameters["t_ss"] = 0.03365982238611262;
    model_parameters["t_sp"] = -0.029154833035109226;
    model_parameters["t_pp1"] = -0.0804163845390335;
    model_parameters["t_pp2"] = -0.01393611496959445;
    model_parameters["r_pseudo"] = 2.60342991362958;
    model_parameters["v_pseudo"] = 0.022972992186364977;
    model_parameters["dipole"] = 2.781629275106456;
    model_parameters["energy_s"] = 3.1659446174413004;
    model_parameters["energy_p"] = -2.3926873325346554;
    model_parameters["coulomb_s"] = 0.3603533286088998;
    model_parameters["coulomb_p"] = -0.003267991835806299;

    Eigen::MatrixXd one_particle_hamiltonian(8, 8);
    Eigen::MatrixXd interaction_matrix(8, 8);
    Eigen::MatrixXd density_matrix(8, 8);
    one_particle_hamiltonian << 2.317455540100693, -0.141367040480457, -0.18848938730727602, -0.23561173413409503, 0.0006538083460047041, 0.0005340767672242702, 0.0007121023562990268, 0.0008901279453737837, -0.141367040480457, -3.241176409875263, 0.0, 0.0, -0.0005340767672242702, 0.00029247941230302376, 0.0021734004930512737, 0.002716750616314092, -0.18848938730727602, 0.0, -3.241176409875263, 0.0, -0.0007121023562990268, 0.0021734004930512737, 0.0015602963665829331, 0.0036223341550854563, -0.23561173413409503, 0.0, 0.0, -3.241176409875263, -0.0008901279453737837, 0.002716750616314092, 0.0036223341550854563, 0.0031903467363713894, 0.0006538083460047041, -0.0005340767672242702, -0.0007121023562990268, -0.0008901279453737837, 2.317455540100693, 0.141367040480457, 0.18848938730727602, 0.23561173413409503, 0.0005340767672242702, 0.00029247941230302376, 0.0021734004930512737, 0.002716750616314092, 0.141367040480457, -3.241176409875263, 0.0, 0.0, 0.0007121023562990268, 0.0021734004930512737, 0.0015602963665829331, 0.0036223341550854563, 0.18848938730727602, 0.0, -3.241176409875263, 0.0, 0.0008901279453737837, 0.002716750616314092, 0.0036223341550854563, 0.0031903467363713894, 0.23561173413409503, 0.0, 0.0, -3.241176409875263;
    interaction_matrix << 0.3603533286088998, 0.0, 0.0, 0.0, 0.1414213562373095, -0.00848528137423857, -0.01131370849898476, -0.01414213562373095, 0.0, -0.003267991835806299, 0.0, 0.0, 0.00848528137423857, 0.0013010764773832475, -0.0020364675298172566, -0.0025455844122715707, 0.0, 0.0, -0.003267991835806299, 0.0, 0.01131370849898476, -0.0020364675298172566, 0.00011313708498984776, -0.003394112549695428, 0.0, 0.0, 0.0, -0.003267991835806299, 0.01414213562373095, -0.0025455844122715707, -0.003394112549695428, -0.0014142135623730948, 0.1414213562373095, 0.00848528137423857, 0.01131370849898476, 0.01414213562373095, 0.3603533286088998, 0.0, 0.0, 0.0, -0.00848528137423857, 0.0013010764773832475, -0.0020364675298172566, -0.0025455844122715707, 0.0, -0.003267991835806299, 0.0, 0.0, -0.01131370849898476, -0.0020364675298172566, 0.00011313708498984776, -0.003394112549695428, 0.0, 0.0, -0.003267991835806299, 0.0, -0.01414213562373095, -0.0025455844122715707, -0.003394112549695428, -0.0014142135623730948, 0.0, 0.0, 0.0, -0.003267991835806299;
    density_matrix << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0;

    //one_particle_hamiltonian = Eigen::MatrixXd::Zero(8, 8);
    std::cout << "one_particle_hamiltonian" << std::endl << one_particle_hamiltonian << std::endl;
    std::cout << "interaction_matrix" << std::endl << interaction_matrix << std::endl;
    std::cout << "density_matrix" << std::endl << density_matrix << std::endl;
    //interaction_matrix = Eigen::MatrixXd::Zero(8, 8);
    
    Eigen::MatrixXd fock = calc_fock_matrix_fast(one_particle_hamiltonian, interaction_matrix, density_matrix, model_parameters, 4);
    //std::cout << fock << std::endl;

    //std::cout << "Testing atom" << std::endl;
    /* for (int i = 0; i < 8; i++)
    {
        std::cout << atom(i, 4) << std::endl;
    }*/
    std::vector<std::string> orbital_types;
    orbital_types.push_back("s");
    orbital_types.push_back("px");
    orbital_types.push_back("py");
    orbital_types.push_back("pz");
    //std::cout << "Testing orb" << std::endl;
    /*for (int i = 0; i < 8; i++)
    {
        std::cout << orb(i, orbital_types, 4) << std::endl;
    }*/

    /* std::cout << "Testing ao_index" << std::endl;
    for (int i = 0; i < 2; i++)
    {
        for (auto it : orbital_types)
        {
            std::cout << ao_index(i, it, orbital_types, 4) << std::endl;
        }
    }*/

    /*std::cout << "Testing chi_on_atom" << std::endl;
    for (auto i : orbital_types)
    {
        for (auto j : orbital_types)
        {
            for (auto k : orbital_types)
            {
                std::cout << chi_on_atom(i, j, k, model_parameters) << std::endl;
            }
        }
    }*/

    return 0;
}


int atom(int ao_index, int orbitals_per_atom)
{
    return floor(ao_index / orbitals_per_atom);
}

std::string orb(int ao_index, std::vector<std::string> orbital_types, int orbitals_per_atom)
{
    int orb_index = ao_index % orbitals_per_atom;
    return orbital_types[orb_index];
}

int ao_index(int atom_p, std::string orb_p, std::vector<std::string> orbital_types, int orbitals_per_atom)
{
    int p = atom_p * orbitals_per_atom;
    std::vector<std::string>::iterator it = std::find(orbital_types.begin(), orbital_types.end(), orb_p);
    int orb_index = std::distance(orbital_types.begin(), it);
    return p + orb_index;
}

double chi_on_atom(std::string orb1, std::string orb2, std::string orb3, std::map<std::string, double> model_parameters)
{
    if (orb1 == orb2 && orb3 == "s")
    {
        return 1.0;
    }
    if (orb1 == orb3 && (orb3 == "px" || orb3 == "py" || orb3 == "pz") && orb2 == "s")
    {
        //std::cout << model_parameters["dipole"] << std::endl;
        return model_parameters["dipole"];
    }
    if (orb2 == orb3 && (orb3 == "px" || orb3 == "py" || orb3 == "pz") && orb1 == "s")
    {
        return model_parameters["dipole"];
    }
    /*else {
        std::cout << "ERROR" << std::endl;
    }*/
    return 0.0;
}

Eigen::MatrixXd calc_fock_matrix_fast(Eigen::MatrixXd hamiltonian_matrix, Eigen::MatrixXd interaction_matrix, Eigen::MatrixXd density_matrix, std::map<std::string, double> model_parameters, int orbitals_per_atom)
{
    std::vector<std::string> orbital_types;
    orbital_types.push_back("s");
    orbital_types.push_back("px");
    orbital_types.push_back("py");
    orbital_types.push_back("pz");
    int ndof = hamiltonian_matrix.rows();
    Eigen::MatrixXd fock_matrix = hamiltonian_matrix;

    // Hartree potential
    for (int p = 0; p < ndof; p++)
    {
        int atom_p = atom(p, orbitals_per_atom);
        std::string orb_p = orb(p, orbital_types, orbitals_per_atom);
        std::cout << "atom_p: " << atom_p << "  orbital: " << orb_p << std::endl;
        for (auto orb_q : orbital_types)
        {
            int atom_p = atom(p, orbitals_per_atom);
            int q = ao_index(atom_p, orb_q, orbital_types, orbitals_per_atom);
            for (auto orb_t : orbital_types)
            {
                int t = ao_index(atom_p, orb_t, orbital_types, orbitals_per_atom);
                double chi_pqt = chi_on_atom(orb_p, orb_q, orb_t, model_parameters);
                for (int r = 0; r < ndof; r++)
                {
                    int atom_r = atom(r, orbitals_per_atom);
                    std::string orb_r = orb(r, orbital_types, orbitals_per_atom);
                    for (auto orb_s : orbital_types)
                    {
                        int s = ao_index(atom_r, orb_s, orbital_types, orbitals_per_atom);
                        for (auto orb_u : orbital_types)
                        {
                            int u = ao_index(atom_r, orb_u, orbital_types, orbitals_per_atom);
                            //std::cout << orb_r << " " << orb_s << " " << orb_u << " " << r << " " << s << " " << t << " " << u << std::endl;
                            double chi_rsu = chi_on_atom(orb_r, orb_s, orb_u, model_parameters);
                            fock_matrix(p, q) += 2.0 * chi_pqt * chi_rsu * interaction_matrix.coeffRef(t, u) * density_matrix.coeffRef(r, s);
                        }
                    }
                }
            }
        }
    }

    // Fock exchange term
    for (int p = 0; p < ndof; p++)
    {
        int atom_p = atom(p, orbitals_per_atom);
        std::string orb_p = orb(p, orbital_types, orbitals_per_atom);
        for (auto orb_s : orbital_types)
        {
            int s = ao_index(atom_p, orb_s, orbital_types, orbitals_per_atom);
            for (auto orb_u : orbital_types)
            {
                int u = ao_index(atom_p, orb_u, orbital_types, orbitals_per_atom);
                double chi_psu = chi_on_atom(orb_p, orb_s, orb_u, model_parameters);
                for (int q = 0; q < ndof; q++)
                {
                    int atom_q = atom(q, orbitals_per_atom);
                    std::string orb_q = orb(q, orbital_types, orbitals_per_atom);
                    for (auto orb_r : orbital_types)
                    {
                        int r = ao_index(atom_q, orb_r, orbital_types, orbitals_per_atom);
                        int atom_r = atom(r, orbitals_per_atom);
                        for (auto orb_t : orbital_types)
                        {
                            int t = ao_index(atom_q, orb_t, orbital_types, orbitals_per_atom);
                            double chi_rqt = chi_on_atom(orb_r, orb_q, orb_t, model_parameters);
                            fock_matrix(p, q) -= chi_rqt * chi_psu * interaction_matrix.coeffRef(t, u) * density_matrix.coeffRef(r, s);
                        }
                    }
                }
            }
        }
    }
    std::cout << "fock matrix:" << std::endl << fock_matrix << std::endl;
    return fock_matrix;
}
