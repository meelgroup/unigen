#include <unigen/unigen.h>
#include <approxmc/approxmc.h>
#include <arjun/arjun.h>
#include <iostream>
#include <memory>
#include <vector>
using std::vector;
using CMSat::Lit;

void print_sample(const vector<int>& sample, void*) {
    for (const int lit: sample) std::cout << lit << " ";
    std::cout << "0" << std::endl;
}

int main() {
    std::unique_ptr<CMSat::FieldGen> fg = std::make_unique<ArjunNS::FGenMpz>();
    ArjunNS::SimplifiedCNF cnf(fg);
    cnf.set_need_aig();
    cnf.new_vars(4);
    cnf.add_clause({Lit(0, false), Lit(1, false)}); // x1 v x2
    cnf.add_clause({Lit(2, true), Lit(3, false)});  // -x3 v x4
    cnf.add_clause({Lit(3, true), Lit(0, false)});  // -x4 v x1
    cnf.set_sampl_vars(vector<uint32_t>{0, 1, 2});

    ApproxMC::AppMC appmc(fg);
    UniGen::UniG unigen(&appmc);
    unigen.set_callback(print_sample, nullptr);
    if (!unigen.simplify_with_arjun(cnf)) {
        std::cout << "UNSAT" << std::endl;
        return 1;
    }
    // x4 is not sampled: its value comes from Arjun's definition of it
    unigen.set_full_sampling_vars({0, 1, 2, 3});
    const auto sol_count = appmc.count();
    unigen.sample(&sol_count, 5);
    return 0;
}
