#include <string>
#include <unordered_map>
#include <vector>

#ifndef CLASSES_H
#define CLASSES_H

class Chemical {
    private:
        static int default_id;
        
    public:
        double diff_const;
        std::string id;
        Chemical();
        Chemical(double diff_const, std::string id);

        struct ChemicalHashFunction {
            size_t operator()(const Chemical& chemical) const;
        };

        struct ChemicalEqualFunction {
            bool operator()(const Chemical& lhs, const Chemical& rhs) const;
        };

        bool operator==(const Chemical& other) const;
        bool operator!=(const Chemical& other) const;
};

using ChemicalMap = std::unordered_map<Chemical, double, Chemical::ChemicalHashFunction, Chemical::ChemicalEqualFunction>;

class Reaction {
    private:
        
        
    public:
        double reaction_const;
        std::string id;
        std::vector<Chemical> chems_list;
        ChemicalMap chems_in, chems_out;
        
        Reaction(ChemicalMap chems_in_param, ChemicalMap chems_out_param, double reaction_const = 1.0);

        std::pair<ChemicalMap, ChemicalMap> pred_prey(ChemicalMap chems_val);

};

class Rod {
    private:
        std::vector<Chemical> chems_list;
        std::vector<Reaction> reactions;
        std::vector<double> x;
        std::unordered_map<Chemical, std::vector<double>, Chemical::ChemicalHashFunction, Chemical::ChemicalEqualFunction> chems, chems_next;
        std::string id;
        double dx, dt, length, temp, coef;

    public:

        int nx;

        // Rod constructor
        Rod(double dt, double dx, double length = 1, double temp = 1);

        // add a single reaction to Rod
        void add_reaction(Reaction react);
        // add multiple reactions to Rod through vector
        void add_reaction(std::vector<Reaction> reacts);
        // add chemical, pass chemical, region (position in x vector where to place chemical), vals (respective values for each place in x)
        void add_chem(Chemical chem, std::vector<int> region, std::vector<double> vals);
        // transfer chemicals in region to destination, if range_in is empty, take all chemicals in current Rod, if range_out empty, transfer evenly to all space in target rod
        void transfer_to(Rod destination, Chemical chem, std::vector<int> range_in, std::vector<int> range_out);
        // integrate the total volume of a single chemical
        ChemicalMap integrate_chems(Chemical chem);
        // integrate the total volume of a vector of chemicals, pass empty vector for all chemicals in rod
        ChemicalMap integrate_chems(std::vector<Chemical> chems);
        // iterate through the solution
        void update();

};

#endif 