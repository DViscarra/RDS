#include <classes.h>
#include <helpers.h>
#include <iostream>
#include <string>
#include <math.h>
#include <stdexcept>
#include <algorithm>

Chemical::Chemical() {
    default_id++;
    id = char(default_id);
    diff_const = 1.0;
}

Chemical::Chemical(double diff_const_param, std::string id_param) {
    id = id_param;
    diff_const = diff_const_param;
}

size_t Chemical::ChemicalHashFunction::operator()(const Chemical& chemical) const {
    return std::hash<std::string>{}(chemical.id);
}

bool Chemical::ChemicalEqualFunction::operator()(const Chemical& lhs, const Chemical& rhs) const {
    return lhs.id == rhs.id && lhs.diff_const == rhs.diff_const;
}

bool Chemical::operator==(const Chemical& other) const {
    return id == other.id;
}

bool Chemical::operator!=(const Chemical& other) const {
    return !(*this == other);
}

int Chemical::default_id = 64;

Reaction::Reaction(ChemicalMap chems_in_param, ChemicalMap chems_out_param, double reaction_const_param) {
    chems_in = chems_in_param;
    chems_out = chems_out_param;

    for (auto chem : chems_in) {
        chems_list.push_back(chem.first);
    }
    for (auto chem : chems_out) {
        for (auto chem_in : chems_list) {
            if (!(chem.first.id == chem_in.id && chem.first.diff_const == chem_in.diff_const)) {
                chems_list.push_back(chem.first);
            }
        }
    }

    reaction_const = reaction_const_param;

    id = "";
    bool first = true;
    if(!chems_in.empty()) {
        for (auto chem : chems_in) {
            Chemical c = chem.first;
            if (first) {
                id = id + double_to_string(chems_in[c], 2) + c.id + " ";
                first = false;
            } else id = id + "+ " + double_to_string(chems_in[c], 2) + c.id + " ";
        }
    } else id = id + "Ø ";

    id = id + "-(" + double_to_string(reaction_const, 2) + ")> ";

    first = true;
    if (!chems_out.empty()) {
        for (auto chem : chems_out) {
            Chemical c = chem.first;
            if (first) {
                id = id + double_to_string(chems_out[c], 2) + c.id + "";
                first = false;
            } else id = id + "+ " + double_to_string(chems_out[c], 2) + c.id + " ";
        }
    } else id = id + "Ø ";
}

std::pair<ChemicalMap, ChemicalMap> Reaction::pred_prey(ChemicalMap chems_val) {
    double c = reaction_const;
    for (std::pair<const Chemical, double> chem_in : chems_in) {
        Chemical chem = chem_in.first;
        c = c * pow(chems_val[chem], chems_in[chem]);
    }

    //consumption
    ChemicalMap consumption;
    double val;
    for (Chemical chem : chems_list) {
        if (chems_in.find(chem) != chems_in.end()) {
            if (chems_out.find(chem) != chems_out.end()) {
                if (chems_in[chem] > chems_out[chem]) val = chems_in[chem] - chems_out[chem];
                else val = 0.0;
            } else val = chems_in[chem];
        } else val = 0.0;
        consumption[chem] = c * val;
    }

    //creation
    ChemicalMap creation;
    val = 0;
    if (chems_in.empty()) {
        for (std::pair<const Chemical, double> chem_out : chems_out) {
            Chemical chem = chem_out.first;
            creation[chem] = reaction_const * chems_out[chem];
        }
    }
    for (Chemical chem : chems_list) {
        if (chems_out.find(chem) != chems_out.end()) {
            if (chems_in.find(chem) != chems_in.end()) {
                if (chems_in[chem] < chems_out[chem]) val = chems_out[chem] - chems_in[chem];
                else val = 0.0;
            } else val = chems_out[chem];
        } else val = 0.0;
        creation[chem] = c * val;
    }

    std::pair<ChemicalMap, ChemicalMap> ret;
    ret.first = consumption; 
    ret.second = creation;
    return ret;
}

Rod::Rod(double dt_param, double dx_param, double length_param, double temp_param) {
    dx = dx_param;
    dt = dt_param;
    length = length_param;
    nx = int(length/dx);
    temp = temp_param;
    coef = dt/(dx*dx);
    for (int i = 0; i < nx; i++) {
        x.push_back(i*dx);
    }
}

void Rod::add_reaction(Reaction react) {
    reactions.push_back(react);
}

void Rod::add_reaction(std::vector<Reaction> reacts) {
    reactions.insert(reactions.end(), reacts.begin(), reacts.end());
}

void Rod::add_chem(Chemical chem, std::vector<int> region, std::vector<double> vals) {
    if (region.size() != vals.size()) throw std::invalid_argument("Provide region and vals of the same size");

    if (std::find(chems_list.begin(), chems_list.end(), chem) == chems_list.end()) {
        chems_list.push_back(chem);
        std::vector<double> v(nx, 0.0), vn(nx, 0.0);
        chems[chem] = v;
        chems[chem] = vn;
    }

    for (int i = 0; i < region.size(); i++) {
        chems[chem][region[i]] += vals[i];
    }
}

void Rod::transfer_to(Rod destination, Chemical chem, std::vector<int> range_in, std::vector<int> range_out) {
    if (range_in.empty()) {
        for (int i = 0; i < nx; i++) range_in.push_back(i); 
    }

    if (range_out.empty()) {
        for (int i = 0; i < destination.nx; i++) range_out.push_back(i); 
    }

    double transfer_vol = 0;
    for (int i = 0; i < range_in.size(); i++) {
        transfer_vol += chems[chem][i];
        chems[chem][range_in[i]] = 0;
    }

    // Different from transfer_vol
    std::vector<double> transfer_vals;
    for (int i = 0; i < range_out.size(); i++) transfer_vals[i] = transfer_vol/range_out.size();
    destination.add_chem(chem, range_out, transfer_vals);
}

ChemicalMap Rod::integrate_chems(Chemical chem) {
    ChemicalMap integrated;
    integrated[chem] = 0;
    for (double value : chems[chem]) integrated[chem] += value*dx;

    return integrated;
}

ChemicalMap Rod::integrate_chems(std::vector<Chemical> chems_int) {
    ChemicalMap integrated;
    if (chems_int.empty()) chems_int = chems_list;

    for(Chemical chem : chems_int) {
        integrated[chem] = 0;
        for (double value : chems[chem]) integrated[chem] += value*dx;
    }

    return integrated;
}

void Rod::update() {
    for (int i = 1; i < nx-1; i++) {
        ChemicalMap cons, cre;
        for (Reaction react : reactions) {
            if (is_subset(react.chems_list, chems_list)) {
                ChemicalMap chems_val;
                std::vector<Chemical> c_i;
                for (std::pair<Chemical, double> p : react.chems_in) c_i.push_back(p.first);
                for (Chemical chem : c_i) chems_val[chem] = chems[chem][i];
                std::pair<ChemicalMap, ChemicalMap> temp = react.pred_prey(chems_val);
                ChemicalMap cons_temp = temp.first, cre_temp = temp.second;
                for (std::pair<Chemical, double> j : cons_temp) {
                    cons[j.first] += cons_temp[j.first];
                    cre[j.first] += cre_temp[j.first];
                }
            } else throw std::runtime_error("Missing chemicals in rod for reaction");
        }
        for (Chemical chem : chems_list) {
            if (cons.find(chem) == cons.end()) cons[chem] = 0.0;
            if (cre.find(chem) == cre.end()) cre[chem] = 0.0;

            chems_next[chem][i] = coef*chem.diff_const*(chems[chem][i-1] - 2*chems[chem][i] + chems[chem][i+1])
                                  + chems[chem][i]
                                  + dt*(cre[chem] - cons[chem]);
        }
    }

    // Boundary conditions
    for (Chemical chem : chems_list) {
        chems_next[chem][0] = chems_next[chem][1];
        chems_next[chem][nx-1] = chems_next[chem][nx-2];
    }

    // Move chems_next into chems
    for (Chemical chem : chems_list) for (int i = 0; i < chems.size(); i++) chems[chem][i] = chems_next[chem][i];
}