/*
 UniGen

 Copyright (c) 2019-2020, Mate Soos and Kuldeep S. Meel. All rights reserved
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

#include "argparse.hpp"
using std::string;
using std::vector;

#if defined(__GNUC__) && defined(__linux__)
#include <fenv.h>
#endif

#include "unigen.h"
#include "time_mem.h"
#include <approxmc/approxmc.h>
#include <fstream>
#include <set>
#include <random>

#include <cryptominisat5/dimacsparser.h>
#include <cryptominisat5/streambuffer.h>
#include <arjun/arjun.h>

using namespace CMSat;
using namespace UniGen;
using std::cout;
using std::endl;
using std::string;
using std::set;
ApproxMC::AppMC* appmc = NULL;
UniG* unigen = NULL;
argparse::ArgumentParser program = argparse::ArgumentParser("unigen",
        UniGen::UniG::get_version_sha1(),
        argparse::default_arguments::help);
std::unique_ptr<CMSat::FieldGen> fg;

uint32_t verb = 1;
uint32_t seed;
double epsilon;
double delta;
uint32_t verb_banning_cls = 0;
uint32_t simplify;
double var_elim_ratio;
uint32_t reuse_models = 1;
uint32_t sparse;

// Arjun
int do_arjun = 1;
UniGen::ArjunConf aconf;

// UniGen
uint32_t num_samples = 500;
int multisample;
string sample_fname;
double kappa = 0.638;  /* Corresponds to UniGen's epsilon=16 in the TACAS-15 paper */
bool verb_sampler_cls ;
int full_sol = 0;

#define myopt(name, var, fun, hhelp) \
    program.add_argument(name) \
        .action([&](const auto& a) {var = std::fun(a.c_str());}) \
        .default_value(var) \
        .help(hhelp)
#define myopt2(name1, name2, var, fun, hhelp) \
    program.add_argument(name1, name2) \
        .action([&](const auto& a) {var = std::fun(a.c_str());}) \
        .default_value(var) \
        .help(hhelp)

void print_version() {
    cout << "c o UniGen SHA1: " << UniGen::UniG::get_version_sha1() << endl;
    cout << "c o ApproxMC SHA1: " << ApproxMC::AppMC::get_version_sha1() << endl;
    cout << "c o Arjun SHA1: " << ArjunNS::Arjun::get_version_sha1() << endl;
    cout << "c o SBVA SHA1: " << ArjunNS::Arjun::get_sbva_version_sha1() << endl;
    cout << "c o CMS SHA1: " << CMSat::SATSolver::get_version_sha1() << endl;
    cout << "c o CaDiCaL SHA1: " << CMSat::SATSolver::get_cadical_version_sha1() << endl;
    cout << "c o CadiBack SHA1: " << CMSat::SATSolver::get_cadiback_version_sha1() << endl;
    cout << CMSat::SATSolver::get_thanks_info("c o ");
    cout << ArjunNS::Arjun::get_thanks_info("c o ");
}

void SIGINT_handler(int) {
    return; //Perhaps we should output all the generated samples so far
}

void add_unigen_options() {
    ApproxMC::AppMC tmp(fg);
    epsilon = tmp.get_epsilon();
    delta = tmp.get_delta();
    simplify = tmp.get_simplify();
    var_elim_ratio = tmp.get_var_elim_ratio();
    sparse = tmp.get_sparse();
    seed = tmp.get_seed();

    myopt2("-v", "--verb", verb, atoi, "Verbosity");
    program.add_argument("--version") \
        .action([&](const auto&) {print_version(); exit(0);}) \
        .flag()
        .help("Print version and exit");
    myopt2("-s", "--seed", seed, atoi, "Seed");
    myopt2("-e", "--epsilon", epsilon, stod,
            "Tolerance parameter, i.e. how close is the count from the correct count? "
            "Count output is within bounds of (exact_count/(1+e)) < count < (exact_count*(1+e)). "
            "So e=0.8 means we'll output at most 180%% of exact count and at least 55%% of exact count. "
            "Lower value means more precise.");
    myopt2("-d", "--delta", delta, stod, "Confidence parameter, i.e. how sure are we of the result? "
            "(1-d) = probability the count is within range as per epsilon parameter. "
            "So d=0.2 means we are 80%% sure the count is within range as specified by epsilon. "
            "The lower, the higher confidence we have in the count.");
    myopt("--kappa", kappa, atof, "Uniformity parameter (see TACAS-15 paper)");

    myopt("--arjun", do_arjun, atoi, "Use arjun to minimize sampling set");
    myopt("--withe", aconf.elim, atoi, "With --arjun: also simplify the CNF, eliminating vars. "
            "Samples are mapped back to the input via Arjun's definitions of the removed vars");
    myopt("--backbone", aconf.backbone, atoi, "With --arjun: run backbone analysis first");
    myopt("--sbva", aconf.etof.num_sbva_steps, atoi, "With --withe: SBVA steps (in 1000s)");
    myopt("--eappmc", aconf.simp.appmc, atoi, "With --withe: simplify for ApproxMC's XOR hashing");
    myopt("--evivif", aconf.simp.oracle_vivify, atoi, "With --withe: oracle vivification");
    myopt("--eiter1", aconf.simp.iter1, atoi, "With --withe: simplification iterations, 1st round");
    myopt("--eiter2", aconf.simp.iter2, atoi, "With --withe: simplification iterations, 2nd round");
    myopt("--autarky", aconf.etof.do_autarky, atoi, "With --withe: autarky");
    myopt("--fullsol", full_sol, atoi, "Print a value for every var of the input, not only the sampling set");
    myopt("--sparse", sparse, atoi, "Generate sparse XORs when possible");
    myopt("--reusemodels", reuse_models, atoi, "Reuse models while counting solutions");
    myopt("--verbanbcls", verb_banning_cls, atoi ,"Print banning clause + xor clauses. Highly verbose.");
    myopt("--simplify", simplify, atoi , "Simplify agressiveness");
    myopt("--velimratio", var_elim_ratio, atof, "Variable elimination ratio for each simplify run");
    myopt("--samples", num_samples, atoi, "Number of random samples to generate");
    myopt("--multisample", multisample, atoi, "Return multiple samples from each call");
    myopt("--sampleout", sample_fname, string, "Write samples to this file");
    myopt("--verbsamplercls", verb_sampler_cls, atoi, "Print XOR constraints added for sampling");

    program.add_argument("inputfile").remaining().help("input CNF");
}

void parse_supported_options(int argc, char** argv) {
    add_unigen_options();
    try {
        program.parse_args(argc, argv);
        if (program.is_used("--help")) {
            cout << "Probilistic Uniform Sample Generator" << endl << endl
            << "approxmc [options] inputfile" << endl;
            cout << program << endl;
            exit(0);
        }
    }
    catch (const std::exception& err) {
        std::cerr << err.what() << std::endl;
        exit(-1);
    }
}

template<class T> void parse_file(const std::string& filename, T* reader) {
  #ifndef USE_ZLIB
  FILE * in = fopen(filename.c_str(), "rb");
  DimacsParser<StreamBuffer<FILE*, CMSat::FN>, T> parser(reader, nullptr, 0, fg);
  #else
  gzFile in = gzopen(filename.c_str(), "rb");
  DimacsParser<StreamBuffer<gzFile, CMSat::GZ>, T> parser(reader, nullptr, 0, fg);
  #endif
  if (in == nullptr) {
      std::cout << "ERROR! Could not open file '" << filename
      << "' for reading: " << strerror(errno) << endl;
      std::exit(-1);
  }
  if (!parser.parse_DIMACS(in, true)) exit(-1);
  #ifndef USE_ZLIB
  fclose(in);
  #else
  gzclose(in);
  #endif

  if (!reader->get_sampl_vars_set()) {
    vector<uint32_t> tmp;
    for(uint32_t i = 0; i < reader->nVars(); i++) tmp.push_back(i);
    reader->set_sampl_vars(tmp);
  } else {
    for(auto const& s: reader->get_sampl_vars()) {
      if (s >= reader->nVars()) {
        cout << "ERROR: Sampling var " << s+1 << " is larger than number of vars in formula: "
          << reader->nVars() << endl;
        exit(-1);
      }
    }
    if (!reader->get_opt_sampl_vars_set()) {
      reader->set_opt_sampl_vars(reader->get_sampl_vars());
    }
  }
}

void mycallback(const std::vector<int>& solution, void *data)
{
    std::ostream& os = *(std::ostream*)data;
    for(const int l: solution) os << l << " ";
    os << "0" << endl;
}

inline double stats_line_percent(double num, double total)
{
    if (total == 0) {
        return 0;
    } else {
        return num/total*100.0;
    }
}

void print_final_indep_set(const vector<uint32_t>& indep_set, uint32_t orig_sampling_set_size)
{
    cout << "c o ind ";
    for(const uint32_t s: indep_set) cout << s+1 << " ";
    cout << "0" << endl;

    cout
    << "c o [arjun] final set size:      " << std::setw(7) << indep_set.size()
    << " percent of original: "
    <<  std::setw(6) << std::setprecision(4)
    << stats_line_percent(indep_set.size(), orig_sampling_set_size)
    << " %" << endl;
}

template <class T>
void check_sanity_sampling_vars(T vars, const uint32_t nvars)
{
    for(const auto& v: vars) if (v >= nvars) {
        cout << "ERROR: sampling set provided is incorrect, it has a variable in it: " << v+1 << " that is larger than the total number of variables: " << nvars << endl;
        exit(-1);
    }
}

int main(int argc, char** argv) {
    #if defined(__GNUC__) && defined(__linux__)
    feenableexcept(FE_INVALID   |
                   FE_DIVBYZERO |
                   FE_OVERFLOW
                  );
    #endif

    //Reconstruct the command line so we can emit it later if needed
    string command_line;
    for(int i = 0; i < argc; i++) {
        command_line += string(argv[i]);
        if (i+1 < argc) {
            command_line += " ";
        }
    }

    fg = std::make_unique<ArjunNS::FGenMpz>();
    appmc = new ApproxMC::AppMC(fg);
    unigen = new UniG(appmc);
    parse_supported_options(argc, argv);
    if (verb) {
        print_version();
        cout << "c o executed with command line: " << command_line << endl;
    }

    appmc->set_verbosity(verb);
    appmc->set_seed(seed);
    appmc->set_reuse_models(reuse_models);
    appmc->set_sparse(sparse);
    appmc->set_simplify(simplify);
    appmc->set_var_elim_ratio(var_elim_ratio);
    vector<uint32_t> sampling_vars_orig;

    if (!program.is_used("inputfile")) {
        cout << "ERROR: you need to provide an input file. Exiting." << endl;
        exit(-1);
    }
    const auto& files = program.get<std::vector<std::string>>("inputfile");
    if (files.empty()) {
      cout << "ERROR: you provided --inputfile but no file. Strange. Exiting. " << endl;
      exit(-1);
    }
    const string fname(files[0]);
    ArjunNS::SimplifiedCNF cnf(fg);
    unigen->set_verbosity(verb);
    if (do_arjun) {
        cnf.set_need_aig();
        parse_file(fname, &cnf);
        const uint32_t orig_nvars = cnf.nVars();
        sampling_vars_orig = cnf.get_sampl_vars();
        const double my_time = cpuTime();
        if (!unigen->simplify_with_arjun(cnf, aconf)) {
            cout << "c o [unig] The input formula is unsatisfiable." << endl;
            exit(-1);
        }
        print_final_indep_set(appmc->get_sampl_vars(), sampling_vars_orig.size());
        cout << "c o [arjun] Arjun finished. T: " << (cpuTime() - my_time) << endl;
        if (full_sol) {
            vector<uint32_t> all;
            for(uint32_t v = 0; v < orig_nvars; v++) all.push_back(v);
            unigen->set_full_sampling_vars(all);
        }
    } else {
        parse_file(fname, appmc);
        sampling_vars_orig = appmc->get_sampl_vars();
        print_final_indep_set(appmc->get_sampl_vars(), appmc->get_sampl_vars().size());
        vector<uint32_t> print_vars = sampling_vars_orig;
        if (full_sol) {
            print_vars.clear();
            for(uint32_t v = 0; v < appmc->nVars(); v++) print_vars.push_back(v);
        }
        unigen->set_full_sampling_vars(print_vars);
    }

    auto sol_count_unig = appmc->count();
    unigen->set_verb_sampler_cls(verb_banning_cls);
    unigen->set_kappa(kappa);
    unigen->set_multisample(multisample);

    std::ostream* os = &std::cout;
    std::ofstream sample_out;
    if (sample_fname != "") {
        sample_out.open(sample_fname.c_str());
        if (!sample_out.is_open()) {
            cout << "[Sampler] Cannot open samples file '" << sample_fname
                 << "' for writing." << endl;
            exit(-1);
        }
        os = &sample_out;
    }
    unigen->set_callback(mycallback, os);
    unigen->sample(&sol_count_unig, num_samples);

    delete unigen;
    delete appmc;

    return 0;
}
