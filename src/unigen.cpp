/*
 UniGen

 Copyright (c) 2019-2020, Mate Soos and Kuldeep S. Meel. All rights reserved
 Copyright (c) 2009-2018, Mate Soos. All rights reserved.
 Copyright (c) 2015, Supratik Chakraborty, Daniel J. Fremont,
 Kuldeep S. Meel, Sanjit A. Seshia, Moshe Y. Vardi
 Copyright (c) 2014, Supratik Chakraborty, Kuldeep S. Meel, Moshe Y. Vardi

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
 */

#include <ctime>
#include <cstring>
#include <errno.h>
#include <algorithm>
#include <string.h>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <set>
#include <fstream>
#include <sys/stat.h>
#include <string.h>
#include <list>
#include <array>
#include <cmath>
#include <complex>
//#include <coz.h>

#include "approxmc/approxmc.h"
#include "approxmc/approxmcconfig.h"
#include "time_mem.h"
#include "cryptominisat5/cryptominisat.h"
#include "cryptominisat5/solvertypesmini.h"
#include "GitSHA1.h"

using std::cout;
using std::cerr;
using std::endl;
using std::list;
using std::map;

Hash UniGen::add_hash(uint32_t hash_index)
{
    const string randomBits =
        gen_rnd_bits(conf.sampling_set.size(), hash_index);

    vector<uint32_t> vars;
    for (uint32_t j = 0; j < conf.sampling_set.size(); j++) {
        if (randomBits[j] == '1') {
            vars.push_back(conf.sampling_set[j]);
        }
    }

    solver->new_var();
    const uint32_t act_var = solver->nVars()-1;
    const bool rhs = gen_rhs();
    Hash h(act_var, vars, rhs);

    vars.push_back(act_var);
    solver->add_xor_clause(vars, rhs);
    if (conf.verb_UniGen_cls) {
        print_xor(vars, rhs);
    }

    return h;
}

void UniGen::ban_one(const uint32_t act_var, const vector<lbool>& model)
{
    vector<Lit> lits;
    lits.push_back(Lit(act_var, false));
    for (const uint32_t var: conf.sampling_set) {
        lits.push_back(Lit(var, model[var] == l_True));
    }
    solver->add_clause(lits);
}

///adding banning clauses for repeating solutions
uint64_t UniGen::add_glob_banning_cls(
    const HashesModels* hm
    , const uint32_t act_var
    , const uint32_t num_hashes)
{
    if (hm == NULL)
        return 0;

    assert(act_var != std::numeric_limits<uint32_t>::max());
    assert(num_hashes != std::numeric_limits<uint32_t>::max());

    uint64_t repeat = 0;
    vector<Lit> lits;
    for (uint32_t i = 0; i < hm->glob_model.size(); i++) {
        const SavedModel& sm = hm->glob_model[i];
        //Model was generated with 'sm.hash_num' active
        //We will have 'num_hashes' hashes active

        if (sm.hash_num >= num_hashes) {
            ban_one(act_var, sm.model);
            repeat++;
        } else if ((int)num_hashes - (int)sm.hash_num < 9) {
            //Model has to fit all hashes
            bool ok = true;
            uint32_t checked = 0;
            for(const auto& h: hm->hashes) {
                //This hash is number: h.first
                //Only has to match hashes below current need
                //note that "h.first" is numbered from 0, so this is a "<" not "<="
                if (h.first < num_hashes) {
                    checked++;
                    ok &= check_model_against_hash(h.second, sm.model);
                    if (!ok) break;
                }
            }
            if (ok) {
                //cout << "Found repeat model, had to check " << checked << " hashes" << endl;
                ban_one(act_var, sm.model);
                repeat++;
            }
        }
    }
    return repeat;
}

SolNum UniGen::bounded_sol_count(
        uint32_t maxSolutions,
        const vector<Lit>* assumps,
        const uint32_t hashCount,
        uint32_t minSolutions,
        HashesModels* hm,
        vector<string>* out_solutions
) {
    cout << "c [UniGen] "
    "[ " << std::setw(7) << std::setprecision(2) << std::fixed
    << (cpuTimeTotal()-startTime)
    << " ]"
    << " bounded_sol_count looking for " << std::setw(4) << maxSolutions << " solutions"
    << " -- hashes active: " << hashCount << endl;

    //Will we need to extend the solution?
    bool only_indep_sol = true;
    if (out_solutions != NULL) {
        only_indep_sol = conf.only_indep_samples;
    }

    //Turn off improvement from ApproxMC4 research paper
    if (conf.force_sol_extension) {
        only_indep_sol = false;
    }

    //Set up things for adding clauses that can later be removed
    vector<Lit> new_assumps;
    if (assumps) {
        assert(assumps->size() == hashCount);
        new_assumps = *assumps;
    } else {
        assert(hashCount == 0);
    }
    solver->new_var();
    const uint32_t sol_ban_var = solver->nVars()-1;
    new_assumps.push_back(Lit(sol_ban_var, true));

    if (conf.simplify >= 2) {
        if (conf.verb >= 2) {
            cout << "c [UniGen] inter-simplifying" << endl;
        }
        double myTime = cpuTime();
        solver->simplify(&new_assumps);
        solver->set_verbosity(0);
        total_inter_simp_time += cpuTime() - myTime;
        if (conf.verb >= 1) {
            cout << "c [UniGen] inter-simp finished, total simp time: "
            << total_inter_simp_time << endl;
        }

    }

    const uint64_t repeat = add_glob_banning_cls(hm, sol_ban_var, hashCount);
    uint64_t solutions = repeat;
    double last_found_time = cpuTimeTotal();
    vector<vector<lbool>> models;
    while (solutions < maxSolutions) {
        lbool ret = solver->solve(&new_assumps, only_indep_sol);
        //COZ_PROGRESS_NAMED("one solution")
        assert(ret == l_False || ret == l_True);

        if (conf.verb >= 2) {
            cout << "c [UniGen] bounded_sol_count ret: " << std::setw(7) << ret;
            if (ret == l_True) {
                cout << " sol no.  " << std::setw(3) << solutions;
            } else {
                cout << " No more. " << std::setw(3) << "";
            }
            cout << " T: "
            << std::setw(7) << std::setprecision(2) << std::fixed
            << (cpuTimeTotal()-startTime)
            << " -- hashes act: " << hashCount
            << " -- T since last: "
            << std::setw(7) << std::setprecision(2) << std::fixed
            << (cpuTimeTotal()-last_found_time)
            << endl;
            if (conf.verb >= 3) {
                solver->print_stats();
            }
            last_found_time = cpuTimeTotal();
        }

        if (ret != l_True) {
            break;
        }

        //Add solution to set
        solutions++;
        const vector<lbool> model = solver->get_model();
        //#ifdef SLOW_DEBUG
        check_model(model, hm, hashCount);
        //#endif
        models.push_back(model);
        if (out_solutions) {
            out_solutions->push_back(get_solution_str(model));
        }

        //ban solution
        vector<Lit> lits;
        lits.push_back(Lit(sol_ban_var, false));
        for (const uint32_t var: conf.sampling_set) {
            assert(solver->get_model()[var] != l_Undef);
            lits.push_back(Lit(var, solver->get_model()[var] == l_True));
        }
        if (conf.verb_UniGen_cls) {
            cout << "c [UniGen] Adding banning clause: " << lits << endl;
        }
        solver->add_clause(lits);
    }

    if (solutions < maxSolutions) {
        //Sampling -- output a random sample of N solutions
        if (solutions >= minSolutions && samples_out != NULL) {
            assert(minSolutions > 0);
            vector<size_t> modelIndices;
            for (uint32_t i = 0; i < models.size(); i++) {
                modelIndices.push_back(i);
            }
            std::shuffle(modelIndices.begin(), modelIndices.end(), randomEngine);

            for (uint32_t i = 0; i < sols_to_return(solutions); i++) {
                const auto& model = models.at(modelIndices.at(i));
                (*samples_out) << get_solution_str(model) << endl << std::flush;
            }
        }
    }

    //Save global models
    if (hm && conf.reuse_models) {
        for (const auto& model: models) {
            hm->glob_model.push_back(SavedModel(hashCount, model));
        }
    }

    //Remove solution banning
    vector<Lit> cl_that_removes;
    cl_that_removes.push_back(Lit(sol_ban_var, false));
    solver->add_clause(cl_that_removes);

    return SolNum(solutions, repeat);
}

void UniGen::print_final_count_stats(SATCount solCount)
{
    cout << "c [UniGen] FINISHED UniGen T: "
    << (cpuTimeTotal() - startTime) << " s"
    << endl;

    if (solCount.hashCount == 0 && solCount.cellSolCount == 0) {
        cout << "c [UniGen] Formula was UNSAT " << endl;
    }

    if (conf.verb > 2) {
        solver->print_stats();
    }
    solCount.print_num_solutions();
}

int UniGen::solve(UniGenConfig _conf)
{
    conf = _conf;
    orig_num_vars = solver->nVars();
    startTime = cpuTimeTotal();

    openLogFile();
    randomEngine.seed(conf.seed);
    if (conf.startiter > conf.sampling_set.size()) {
        cerr << "ERROR: Manually-specified startiter for gen_n_samples"
             "is larger than the size of the independent set.\n" << endl;
        exit(-1);
    }

    /* Compute threshold via formula from TACAS-15 paper */
    threshold_UniGengen = ceil(4.03 * (1 + (1/conf.kappa)) * (1 + (1/conf.kappa)));

    //No startiter, we have to figure it out
    if (conf.startiter == 0) {
        std::ostream* backup = samples_out;
        samples_out = NULL;
        AppMC counter = new AppMC;
        AppMCConfig appconf; //should set the default parameters
        SATCount solCount = counter->count(appconf);
        cout << "c [UniGen] finished counting solutions in "
        << (cpuTimeTotal() - startTime) << " s" << endl;

        if (solCount.hashCount == 0 && solCount.cellSolCount == 0) {
            cout << "c [UniGen] The input formula is unsatisfiable." << endl;
            return 0;
        }

        if (conf.verb) {
            solver->print_stats();
        }
        solCount.print_num_solutions();

        double si = round(solCount.hashCount + log2(solCount.cellSolCount)
            + log2(1.8) - log2(threshold_UniGengen)) - 2;
        if (si > 0) {
            conf.startiter = si;
        } else {
            conf.startiter = 0;   /* Indicate ideal sampling case */
        }
        samples_out = backup;
    } else {
        cout << "c Using manually-specified startiter for sample generation" << endl;
    }
    generate_samples();


    return 0;
}

vector<Lit> UniGen::set_num_hashes(
    uint32_t num_wanted,
    map<uint64_t, Hash>& hashes
) {
    vector<Lit> assumps;
    for(uint32_t i = 0; i < num_wanted; i++) {
        if (hashes.find(i) != hashes.end()) {
            assumps.push_back(Lit(hashes[i].act_var, true));
        } else {
            Hash h = add_hash(i);
            assumps.push_back(Lit(h.act_var, true));
            hashes[i] = h;
        }
    }
    assert(num_wanted == assumps.size());

    return assumps;
}

void UniGen::simplify()
{
    if (conf.verb >= 1) {
        cout << "c [UniGen] simplifying" << endl;
    }

    solver->set_sls(1);
    solver->set_intree_probe(1);
    solver->set_full_bve_iter_ratio(conf.var_elim_ratio);
    solver->set_full_bve(1);
    solver->set_bva(1);
    solver->set_distill(1);
    solver->set_scc(1);

    solver->simplify();

    solver->set_sls(0);
    solver->set_intree_probe(0);
    solver->set_full_bve(0);
    solver->set_bva(0);
    solver->set_distill(0);
    //solver->set_scc(0);
}


void UniGen::generate_samples()
{
    assert(samples_out != NULL);
    double genStartTime = cpuTimeTotal();

    hiThresh = ceil(1 + (1.4142136 * (1 + conf.kappa) * threshold_UniGengen));
    loThresh = floor(threshold_UniGengen / (1.4142136 * (1 + conf.kappa)));
    const uint32_t samplesPerCall = sols_to_return(conf.samples);
    const uint32_t callsNeeded =
        conf.samples / samplesPerCall + (bool)(conf.samples % samplesPerCall);
    cout << "Samples requested: " << conf.samples << endl;
    cout << "samples per XOR set:" << samplesPerCall << endl;
    cout << "-> calls needed: " << callsNeeded << endl;

    //TODO WARNING what is this 14???????????????????
    uint32_t callsPerLoop = std::min(solver->nVars() / 14, callsNeeded);
    callsPerLoop = std::max(callsPerLoop, 1U);
    cout << "callsPerLoop:" << callsPerLoop << endl;

    cout << "c [UniGen] starting sample generation."
    << " loThresh: " << loThresh
    << ", hiThresh: " << hiThresh
    << ", startiter: " << conf.startiter << endl;

    uint32_t samples = 0;
    if (conf.startiter > 0) {
        uint32_t lastSuccessfulHashOffset = 0;
        while(samples < conf.samples) {
            samples += gen_n_samples(
                callsPerLoop,
                &lastSuccessfulHashOffset);
        }
    } else {
        /* Ideal sampling case; enumerate all solutions */
        vector<string> out_solutions;
        const uint32_t count = bounded_sol_count(
            std::numeric_limits<uint32_t>::max() //max no. solutions
            , NULL //assumps is empty
            , 0 //number of hashes (information only)
            , 1 //min num. solutions
            , NULL //gobal model (would be banned)
            , &out_solutions
        ).solutions;
        assert(count > 0);

        std::uniform_int_distribution<unsigned> uid {0, count-1};
        for (uint32_t i = 0; i < conf.samples; ++i) {
            vector<string>::iterator it = out_solutions.begin();
            for (uint32_t j = uid(randomEngine); j > 0; --j)    // TODO improve hack
            {
                ++it;
            }
            samples++;
            (*samples_out) << *it << endl << std::flush;
        }
    }

    cout
    << "[UniGen]"
    << " Time to sample: " << cpuTimeTotal() - genStartTime << " s"
    << " -- Time count+samples: " << cpuTimeTotal() << " s"
    << endl;

    cout << "c [UniGen] Samples generated: " << samples << endl;
}

uint32_t UniGen::gen_n_samples(
    const uint32_t num_calls
    , uint32_t* lastSuccessfulHashOffset)
{
    SparseData sparse_data(-1);
    uint32_t num_samples = 0;
    uint32_t i = 0;
    while(i < num_calls) {
        uint32_t hashOffsets[3];
        hashOffsets[0] = *lastSuccessfulHashOffset;

        //Specific values
        if (hashOffsets[0] == 0) { // Starting at q-2; go to q-1 then q
            hashOffsets[1] = 1;
            hashOffsets[2] = 2;
        }
        if (hashOffsets[0] == 2) { // Starting at q; go to q-1 then q-2
            hashOffsets[1] = 1;
            hashOffsets[2] = 0;
        }

        map<uint64_t, Hash> hashes;
        bool ok;
        for (uint32_t j = 0; j < 3; j++) {
            uint32_t currentHashOffset = hashOffsets[j];
            uint32_t currentHashCount = currentHashOffset + conf.startiter;
            const vector<Lit> assumps = set_num_hashes(currentHashCount, hashes, sparse_data);

            double myTime = cpuTime();
            const uint64_t solutionCount = bounded_sol_count(
                hiThresh // max num solutions
                , &assumps //assumptions to use
                , currentHashCount
                , loThresh //min number of solutions (samples not output otherwise)
            ).solutions;
            ok = (solutionCount < hiThresh && solutionCount >= loThresh);
            write_log(
                true, //sampling
                i, currentHashCount, solutionCount == hiThresh,
                      solutionCount, 0, cpuTime()-myTime);

            if (ok) {
                num_samples += sols_to_return(conf.samples);
                *lastSuccessfulHashOffset = currentHashOffset;
                break;
            }
            // Number of solutions too small or too large

            // At q-1, and need to pick next hash count
            if (j == 0 && currentHashOffset == 1) {
                if (solutionCount < loThresh) {
                    // Go to q-2; next will be q
                    hashOffsets[1] = 0;
                    hashOffsets[2] = 2;
                } else {
                    // Go to q; next will be q-2
                    hashOffsets[1] = 2;
                    hashOffsets[2] = 0;
                }
            }
        }

        if (ok) {
            i++;
        }
        if (conf.simplify >= 1) {
            simplify();
        }
    }
    return num_samples;
}

////////////////////
//Helper functions
////////////////////

std::string UniGen::get_solution_str(const vector<lbool>& model)
{
    assert(samples_out != NULL);

    std::stringstream  solution;
    if (conf.only_indep_samples) {
        for (uint32_t j = 0; j < conf.sampling_set.size(); j++) {
            uint32_t var = conf.sampling_set[j];
            assert(model[var] != l_Undef);
            solution << ((model[var] != l_True) ? "-":"") << var + 1 << " ";
        }
    } else {
        for(uint32_t var = 0; var < orig_num_vars; var++) {
            assert(model[var] != l_Undef);
            solution << ((model[var] != l_True) ? "-":"") << var + 1 << " ";
        }
    }
    solution << "0";
    return solution.str();
}

bool UniGen::gen_rhs()
{
    std::uniform_int_distribution<uint32_t> dist{0, 1};
    bool rhs = dist(randomEngine);
    //cout << "rnd rhs:" << (int)rhs << endl;
    return rhs;
}

string UniGen::gen_rnd_bits(
    const uint32_t size,
    // The name of parameter was changed to indicate that this is the index of hash function
    const uint32_t hash_index)
{
    string randomBits;
    std::uniform_int_distribution<uint32_t> dist{0, 1000};
    uint32_t cutoff = 500;

    while (randomBits.size() < size) {
        bool val = dist(randomEngine) < cutoff;
        randomBits += '0' + val;
    }
    assert(randomBits.size() >= size);

    //cout << "rnd bits: " << randomBits << endl;
    return randomBits;
}

void UniGen::print_xor(const vector<uint32_t>& vars, const uint32_t rhs)
{
    cout << "c [UniGen] Added XOR ";
    for (size_t i = 0; i < vars.size(); i++) {
        cout << vars[i]+1;
        if (i < vars.size()-1) {
            cout << " + ";
        }
    }
    cout << " = " << (rhs ? "True" : "False") << endl;
}


void printVersionInfoUniGen()
{
    cout << "c UniGen SHA revision " << ::get_version_sha1() << endl;
    cout << "c UniGen version " << ::get_version_tag() << endl;
    cout << "c UniGen compilation env " << ::get_compilation_env() << endl;
    #ifdef __GNUC__
    cout << "c UniGen compiled with gcc version " << __VERSION__ << endl;
    #else
    cout << "c UniGen compiled with non-gcc compiler" << endl;
    #endif
}

void UniGen::printVersionInfo() const
{
    ::printVersionInfoUniGen();
    cout << solver->get_text_version_info();
}

/* Number of solutions to return from one invocation of gen_n_samples. */
uint32_t UniGen::sols_to_return(uint32_t numSolutions)
{
    if (conf.startiter == 0)   // TODO improve hack for ideal sampling case?
        return numSolutions;
    else if (conf.multisample)
        return loThresh;
    else
        return 1;
}

void UniGen::openLogFile()
{
    if (!conf.logfilename.empty()) {
        logfile.open(conf.logfilename.c_str());
        if (!logfile.is_open()) {
            cout << "[UniGen] Cannot open UniGen log file '" << conf.logfilename
                 << "' for writing." << endl;
            exit(1);
        }

        logfile << std::left
        << std::setw(5) << "sampl"
        << " " << std::setw(4) << "iter"
        << " " << std::setw(4) << "hash"
        << " " << std::setw(4) << "full"
        << " " << std::setw(4) << "sols"
        << " " << std::setw(4) << "rep"
        << " " << std::setw(7) << "T"
        << " " << std::setw(7) << "total T"
        << endl;

    }
}

void UniGen::write_log(
    bool sampling,
    int iter,
    uint32_t hashCount,
    int found_full,
    uint32_t num_sols,
    uint32_t repeat_sols,
    double used_time
)
{
    if (!conf.logfilename.empty()) {
        logfile
        << std::left
        << std::setw(5) << (int)sampling
        << " " << std::setw(4) << iter
        << " " << std::setw(4) << hashCount
        << " " << std::setw(4) << found_full
        << " " << std::setw(4) << num_sols
        << " " << std::setw(4) << repeat_sols
        << " " << std::setw(7) << std::fixed << std::setprecision(2) << used_time
        << " " << std::setw(7) << std::fixed << std::setprecision(2) << (cpuTimeTotal() - startTime)
        << endl;
    }
}


void UniGen::check_model(
    const vector<lbool>& model,
    const HashesModels* const hm,
    const uint32_t hashCount
)
{
    for(uint32_t var: conf.sampling_set) {
        assert(model[var] != l_Undef);
    }

    if (!hm)
        return;

    uint32_t checked = 0;
    bool ok = true;
    for(const auto& h: hm->hashes) {
        //This hash is number: h.first
        //Only has to match hashes at & below
        //Notice that "h.first" is numbered from 0, so it's a "<" not "<="
        if (h.first < hashCount) {
            //cout << "Checking model against hash" << h.first << endl;
            checked++;
            ok &= check_model_against_hash(h.second, model);
            if (!ok) break;
        }
    }
    assert(ok);
}

bool UniGen::check_model_against_hash(const Hash& h, const vector<lbool>& model)
{
    bool rhs = h.rhs;
    for (const uint32_t var: h.hash_vars) {
        assert(model[var] != l_Undef);
        rhs ^= model[var] == l_True;
    }

    //If we started with rhs=FALSE and we XOR-ed in only FALSE
    //rhs is FALSE but we should return TRUE

    //If we started with rhs=TRUE and we XOR-ed in only one TRUE
    //rhs is FALSE but we should return TRUE

    //hence return !rhs
    return !rhs;
}
