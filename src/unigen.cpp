/*
 ApproxMC

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

#include "sampler.h"
#include "GitSHA1.h"
#include "config.h"
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <memory>
#include <cassert>
using namespace ApproxMC;

using std::cout;
using std::endl;

#if defined _WIN32
    #define DLL_PUBLIC __declspec(dllexport)
#else
    #define DLL_PUBLIC __attribute__ ((visibility ("default")))
    #define DLL_LOCAL  __attribute__ ((visibility ("hidden")))
#endif

#define set_get_macro(TYPE, NAME) \
DLL_PUBLIC void UniG::set_##NAME (TYPE NAME) \
{ \
    data->conf.NAME = NAME; \
} \
DLL_PUBLIC TYPE UniG::get_##NAME () const \
{ \
    return data->conf.NAME; \
} \

namespace UniGen {
    // How to map a sample of Arjun's simplified CNF back to the input's vars
    struct ArjunMap {
        uint32_t orig_nvars;
        std::vector<uint32_t> orig_sampl_vars;
        std::vector<CMSat::Lit> orig_to_new;
        std::vector<ArjunNS::aig_lit> defs;
        std::vector<std::vector<CMSat::Lit>> input_cls;
        std::mt19937 rnd;
    };

    struct UniGenPrivateData {
        Sampler sampler;
        Config conf;
        std::unique_ptr<ArjunMap> amap;
    };
}

using namespace UniGen;

DLL_PUBLIC UniG::UniG(AppMC* appmc)
{
    data = new UniGenPrivateData;
    data->sampler.appmc = appmc;
}

DLL_PUBLIC UniG::~UniG()
{
    delete data;
}

DLL_PUBLIC void UniG::set_callback(
    UniGen::callback _callback_func,
    void* _callback_func_data)
{
    data->sampler.callback_func = _callback_func;
    data->sampler.callback_func_data = _callback_func_data;
}

DLL_PUBLIC bool UniG::simplify_with_arjun(ArjunNS::SimplifiedCNF& cnf, const ArjunConf& aconf)
{
    using namespace CMSat;
    if (!cnf.get_need_aig()) {
        cout << "ERROR: the CNF given to simplify_with_arjun() must be set_need_aig()" << endl;
        exit(-1);
    }
    if (!cnf.get_sampl_vars_set()) {
        vector<uint32_t> all;
        for(uint32_t i = 0; i < cnf.nVars(); i++) all.push_back(i);
        cnf.set_sampl_vars(all);
    }
    if (!cnf.get_opt_sampl_vars_set()) cnf.set_opt_sampl_vars(cnf.get_sampl_vars());
    const bool all_indep = cnf.get_sampl_vars().size() == cnf.nVars();

    auto m = std::make_unique<ArjunMap>();
    m->orig_nvars = cnf.nVars();
    m->orig_sampl_vars = cnf.get_sampl_vars();
    m->input_cls = cnf.get_clauses();
    m->rnd.seed(data->sampler.appmc->get_seed());
    {
        // Arjun's definitions need a satisfiable input
        SATSolver s;
        s.new_vars(cnf.nVars());
        for(const auto& c: cnf.get_clauses()) s.add_clause(c);
        if (s.solve() == l_False) return false;
    }

    ArjunNS::Arjun arjun;
    arjun.set_verb(data->conf.verb);
    arjun.set_seed(data->sampler.appmc->get_seed());
    arjun.set_or_gate_based(aconf.gates);
    arjun.set_xor_gates_based(aconf.gates);
    arjun.set_ite_gate_based(aconf.gates);
    arjun.set_irreg_gate_based(aconf.gates);
    if (aconf.backbone) arjun.standalone_backbone(cnf);
    const auto info = arjun.standalone_minimize_indep_info(cnf, aconf.iconf, all_indep);
    // Free in the projection, so Arjun multiplies them out, but they are
    // still in the clauses: their values must come from the sampler too
    std::set<uint32_t> sampl(cnf.get_sampl_vars().begin(), cnf.get_sampl_vars().end());
    for(const auto& v: info.free_vars) {
        if (v >= cnf.nVars() || sampl.count(v)) {
            cout << "ERROR: Arjun gave a wrong free var: " << v+1 << endl;
            exit(-1);
        }
        sampl.insert(v);
    }
    if (!info.free_vars.empty()) cnf.set_multiplier_weight(cnf.fg->one());
    // The independent support is the input of Arjun's definitions, not the
    // input's sampling set: Arjun only keeps defs of vars the inputs define,
    // and needs every opt sampling var to be an input or defined
    cnf.set_sampl_vars(sampl, true);
    cnf.set_opt_sampl_vars(sampl);
    cnf.set_orig_clauses(cnf.get_clauses());
    cnf.set_orig_sampl_vars(cnf.get_sampl_vars());
    if (aconf.elim) {
        auto etof = aconf.etof;
        etof.all_indep = false;
        arjun.standalone_elim_to_file(cnf, etof, aconf.simp, aconf.iconf);
    }

    AppMC* appmc = data->sampler.appmc;
    appmc->new_vars(cnf.nVars());
    appmc->set_sampl_vars(cnf.get_sampl_vars());
    for(const auto& c: cnf.get_clauses()) appmc->add_clause(c);
    for(const auto& c: cnf.get_red_clauses()) appmc->add_red_clause(c);
    appmc->set_multiplier_weight(cnf.get_multiplier_weight());

    m->defs.resize(cnf.num_defs());
    for(uint32_t v = 0; v < cnf.num_defs(); v++) m->defs[v] = cnf.get_def(v);
    m->orig_to_new.assign(cnf.num_defs(), lit_Undef);
    for(const auto& [o, n]: cnf.get_orig_to_new_var()) m->orig_to_new[o] = n;
    data->amap = std::move(m);
    return true;
}

// solution is over all vars of the simplified CNF, in order
static vector<CMSat::lbool> extend_to_orig(const vector<int>& solution, ArjunMap& m) {
    using namespace CMSat;
    vector<lbool> vals(m.defs.size(), l_Undef);
    for(uint32_t o = 0; o < vals.size(); o++) {
        const Lit n = m.orig_to_new[o];
        if (n != lit_Undef) {
            assert(n.var() < solution.size() && std::abs(solution[n.var()]) == (int)n.var()+1);
            vals[o] = boolToLBool(solution[n.var()] > 0) ^ n.sign();
        } else if (m.defs[o] == nullptr) {
            // in no clause of the input, or a sampling var Arjun took out into the
            // multiplier: either way, every value of it is part of a solution
            vals[o] = boolToLBool(m.rnd() & 1);
        }
    }
    std::map<ArjunNS::aig_lit, lbool> cache;
    for(uint32_t o = 0; o < vals.size(); o++) {
        if (m.orig_to_new[o] != lit_Undef || m.defs[o] == nullptr) continue;
        vals[o] = ArjunNS::AIG::evaluate(vals, m.defs[o], m.defs, cache);
        assert(vals[o] != l_Undef);
    }
    // Arjun's definitions are not always right: never give a sample that is not a solution
    for(const auto& cl: m.input_cls) {
        bool sat = false;
        for(const auto& l: cl) if ((vals[l.var()] ^ l.sign()) == l_True) { sat = true; break; }
        if (!sat) {
            cout << "ERROR: extended sample falsifies input clause " << cl
                << ". Arjun gave a wrong definition, this is a bug. Try --withe 0" << endl;
            exit(-1);
        }
    }
    return vals;
}

DLL_PUBLIC void UniG::sample(
    const SolCount* sol_count,
    uint32_t num_samples)
{
    if (data->sampler.callback_func == NULL) {
        cout << "ERROR! You must set the callback function or your samples will be lost" << endl;
        exit(-1);
    }
    if (!data->amap) {
        data->sampler.sample(data->conf, *sol_count, num_samples);
        return;
    }

    ArjunMap& m = *data->amap;
    const vector<uint32_t> out_vars = data->conf.full_sampling_vars.empty()
        ? m.orig_sampl_vars : data->conf.full_sampling_vars;
    for(const auto& v: out_vars) if (v >= m.orig_nvars) {
        cout << "ERROR: sampling var " << v+1 << " is not a var of the CNF" << endl;
        exit(-1);
    }
    Config conf = data->conf;
    conf.full_sampling_vars.clear();
    for(uint32_t v = 0; v < data->sampler.appmc->nVars(); v++) conf.full_sampling_vars.push_back(v);
    const auto user_func = data->sampler.callback_func;
    void* const user_data = data->sampler.callback_func_data;
    data->sampler.callback_func = [&](const vector<int>& solution, void*) {
        const auto vals = extend_to_orig(solution, m);
        vector<int> out;
        for(const uint32_t v: out_vars) out.push_back(vals[v] == CMSat::l_True ? (int)v+1 : -(int)v-1);
        user_func(out, user_data);
    };
    data->sampler.sample(conf, *sol_count, num_samples);
    data->sampler.callback_func = user_func;
}

DLL_PUBLIC string UniG::get_version_sha1()
{
    return UnigenIntNS::get_version_sha1();
}

set_get_macro(double, kappa)
set_get_macro(bool, multisample)
set_get_macro(bool, force_sol_extension)
set_get_macro(const std::vector<uint32_t>&, full_sampling_vars)
set_get_macro(bool, verb_sampler_cls)

DLL_PUBLIC void UniG::set_verbosity(uint32_t verb)
{
    data->conf.verb = verb;
}
