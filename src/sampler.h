/*
Sampler

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


#ifndef Sampler_H_
#define Sampler_H_

#include <fstream>
#include <random>
#include <map>
#include <cstdint>
#ifdef CMS_LOCAL_BUILD
#include "cryptominisat.h"
#include "approxmc.h"
#else
#include <cryptominisat5/cryptominisat.h>
#include <approxmc/approxmc.h>
#endif
#include "unigen.h"
#include "config.h"

using std::string;
using std::vector;
using std::map;
using std::cout;
using std::endl;
using ApproxMC::SolCount;
using namespace CMSat;
using namespace ApproxMC;


struct SavedModel
{
    SavedModel(uint32_t _hash_num, const vector<lbool>& _model) :
        model(_model),
        hash_num(_hash_num)
    {
    }

    vector<lbool> model;
    uint32_t hash_num;
};

struct Hash {
    Hash(uint32_t _act_var, vector<uint32_t>& _hash_vars, bool _rhs) :
        act_var(_act_var),
        hash_vars(_hash_vars),
        rhs(_rhs)
    {}

    Hash()
    {}

    uint32_t act_var;
    vector<uint32_t> hash_vars;
    bool rhs;
};

struct HashesModels {
    map<uint64_t, Hash> hashes;
    vector<SavedModel> glob_model; //global table storing models
};

struct SolNum {
    SolNum(uint64_t _solutions, uint64_t _repeated) :
        solutions(_solutions),
        repeated(_repeated)
    {}
    uint64_t solutions = 0;
    uint64_t repeated = 0;
};

struct SparseData {
    explicit SparseData(int _table_no) :
        table_no(_table_no)
    {}

    uint32_t next_index = 0;
    double sparseprob = 0.5;
    int table_no = -1;
};

class Sampler {
public:
    void sample(
        const Config conf,
        const SolCount sol_count,
        const uint32_t num_samples);
    AppMC* appmc;
    SATSolver* solver = NULL;
    string get_version_info() const;

    ///What to call on samples
    UniGen::callback callback_func = NULL;
    void* callback_func_data = NULL;

private:
    uint32_t loThresh;
    uint32_t hiThresh;
    uint32_t threshold_Samplergen;

    Config conf;
    string gen_rnd_bits(const uint32_t size,
                        const uint32_t numhashes);
    uint32_t sols_to_return(uint32_t numSolutions);
    void add_Sampler_options();
    bool gen_rhs();
    uint32_t gen_n_samples(
        const uint32_t samples
        , uint32_t* lastSuccessfulHashOffset
        , const uint32_t num_samples_needed
    );
    Hash add_hash(uint32_t total_num_hashes);
    string binary(const uint32_t x, const uint32_t length);
    void generate_samples(const uint32_t num_samples);
    SolNum bounded_sol_count(
        uint32_t maxSolutions,
        const vector<Lit>* assumps,
        const uint32_t hashCount,
        uint32_t minSolutions = 1,
        HashesModels* hm = NULL,
        vector<vector<int>>* out_solutions = NULL
    );
    vector<Lit> set_num_hashes(
        uint32_t num_wanted,
        map<uint64_t, Hash>& hashes
    );
    void simplify();

    ////////////////
    //Helper functions
    ////////////////
    void print_xor(const vector<uint32_t>& vars, const uint32_t rhs);
    vector<int> get_solution_ints(const vector<lbool>& model);
    void write_log(
        bool sampling,
        int iter,
        uint32_t hashCount,
        int found_full,
        uint32_t num_sols,
        uint32_t repeat_sols,
        double used_time
    );
    void openLogFile();
    void call_after_parse();
    void ban_one(const uint32_t act_var, const vector<lbool>& model);
    void check_model(
        const vector<lbool>& model,
        const HashesModels* const hm,
        const uint32_t hashCount
    );
    bool check_model_against_hash(const Hash& h, const vector<lbool>& model);
    uint64_t add_glob_banning_cls(
        const HashesModels* glob_model = NULL
        , const uint32_t act_var = std::numeric_limits<uint32_t>::max()
        , const uint32_t num_hashes = std::numeric_limits<uint32_t>::max()
    );

    void readInAFile(SATSolver* solver2, const string& filename);
    void readInStandardInput(SATSolver* solver2);


    //Data so we can output temporary count when catching the signal
    vector<uint64_t> numHashList;
    vector<int64_t> numCountList;
    template<class T> T findMedian(vector<T>& numList);
    template<class T> T findMin(vector<T>& numList);

    ////////////////
    // internal data
    ////////////////
    double startTime;
    std::mt19937 randomEngine;
    uint32_t orig_num_vars;
    double total_inter_simp_time = 0;
    uint32_t threshold; //precision, it's computed
};


#endif //Sampler_H_
