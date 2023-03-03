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

#include <boost/program_options.hpp>
namespace po = boost::program_options;
using std::string;
using std::vector;

#if defined(__GNUC__) && defined(__linux__)
#include <fenv.h>
#endif

#include "config.h"
#include "unigen.h"
#include "time_mem.h"
#include <approxmc/approxmc.h>
#include <fstream>

#include <cryptominisat5/dimacsparser.h>
#include <cryptominisat5/streambuffer.h>
#include <arjun/arjun.h>

using namespace CMSat;
using namespace UniGen;
using std::cout;
using std::cerr;
using std::endl;
ApproxMC::AppMC* appmc = NULL;
ArjunNS::Arjun* arjun = NULL;
UniG* unigen = NULL;

po::options_description main_options = po::options_description("Main options");
po::options_description arjun_options = po::options_description("Arjun options");
po::options_description improvement_options = po::options_description("Improvement options");
po::options_description misc_options = po::options_description("Misc options");
po::options_description sampling_options = po::options_description("Sampling options");
po::options_description help_options;
po::variables_map vm;
po::positional_options_description p;

uint32_t verbosity = 1;
uint32_t seed;
double epsilon;
double delta;
string logfilename;
uint32_t verb_banning_cls = 0;
uint32_t simplify;
double var_elim_ratio;
uint32_t detach_xors = 1;
uint32_t reuse_models = 1;
uint32_t force_sol_extension = 0;
uint32_t sparse;

//Arjun
vector<uint32_t> sampling_vars;
bool sampling_vars_found = false;
int ignore_sampl_set = 0;
int do_arjun = 1;
int debug_arjun = 0;

//sampling
uint32_t num_samples = 500;
int only_indep_samples ;
int multisample;
std::string sample_fname;
double kappa;
bool verb_sampler_cls ;

//signal code
void SIGINT_handler(int)
{
    return; //Perhaps we should output all the generated samples so far
}


void add_UniGen_options()
{
    ApproxMC::AppMC tmp;
    epsilon = tmp.get_epsilon();
    delta = tmp.get_delta();
    simplify = tmp.get_simplify();
    var_elim_ratio = tmp.get_var_elim_ratio();
    sparse = tmp.get_sparse();
    seed = tmp.get_seed();

    UniG tmp2(NULL);
    kappa = tmp2.get_kappa();
    multisample = tmp2.get_multisample();
    only_indep_samples = tmp2.get_only_indep_samples();
    force_sol_extension = tmp2.get_force_sol_extension();
    verb_sampler_cls = tmp2. get_verb_sampler_cls();

    std::ostringstream my_epsilon;
    std::ostringstream my_delta;
    std::ostringstream my_var_elim_ratio;
    std::ostringstream my_kappa;

    my_epsilon << std::setprecision(8) << epsilon;
    my_delta << std::setprecision(8) << delta;
    my_var_elim_ratio << std::setprecision(8) << var_elim_ratio;
    my_kappa << std::setprecision(8) << kappa;

    main_options.add_options()
    ("help,h", "Prints help")
    ("input", po::value< vector<string> >(), "file(s) to read")
    ("verb,v", po::value(&verbosity)->default_value(1), "verbosity")
    ("seed,s", po::value(&seed)->default_value(seed), "Seed")
    ("version", "Print version info")

    ("epsilon", po::value(&epsilon)->default_value(epsilon, my_epsilon.str())
        , "epsilon parameter as per PAC guarantees")
    ("delta", po::value(&delta)->default_value(delta, my_delta.str())
        , "delta parameter as per PAC guarantees; 1-delta is the confidence")
    ("log", po::value(&logfilename),
         "Logs of ApproxMC execution")
    ;

    ArjunNS::Arjun tmpa;
    arjun_options.add_options()
    ("arjun", po::value(&do_arjun)->default_value(do_arjun)
        , "Use arjun to minimize sampling set")
    ("debugarjun", po::value(&debug_arjun)->default_value(debug_arjun)
        , "Use CNF from Arjun, but use sampling set from CNF")
    ;

    improvement_options.add_options()
    ("sparse", po::value(&sparse)->default_value(sparse)
        , "Generate sparse XORs when possible")
    ("detachxor", po::value(&detach_xors)->default_value(detach_xors)
        , "Detach XORs in CMS")
    ("reusemodels", po::value(&reuse_models)->default_value(reuse_models)
        , "Reuse models while counting solutions")
    ("forcesolextension", po::value(&force_sol_extension)->default_value(force_sol_extension)
        , "Use trick of not extending solutions in the SAT solver to full solution")
    ;

    misc_options.add_options()
    ("verbanbcls", po::value(&verb_banning_cls)->default_value(verb_banning_cls)
        ,"Print banning clause + xor clauses. Highly verbose.")
    ("simplify", po::value(&simplify)->default_value(simplify)
        , "Simplify agressiveness")
    ("velimratio", po::value(&var_elim_ratio)->default_value(var_elim_ratio, my_var_elim_ratio.str())
        , "Variable elimination ratio for each simplify run")
    ;

    sampling_options.add_options()
    ("samples", po::value(&num_samples)->default_value(num_samples)
        , "Number of random samples to generate. CAV2020 paper has 500, the default")
    ("nosolext", po::value(&only_indep_samples)->default_value(only_indep_samples)
        , "Should only output the independent vars from the samples")
    ("multisample", po::value(&multisample)->default_value(multisample)
        , "Return multiple samples from each call")
    ("sampleout", po::value(&sample_fname)
        , "Write samples to this file")
    ("kappa", po::value(&kappa)->default_value(kappa, my_kappa.str())
        , "Uniformity parameter (see TACAS-15 paper)")
    ("verbsamplercls", po::value(&verb_sampler_cls)->default_value(verb_sampler_cls)
        , "Print XOR constraints added for sampling")
    ;

    help_options.add(main_options);
    help_options.add(arjun_options);
    help_options.add(sampling_options);
    help_options.add(improvement_options);
    help_options.add(misc_options);
}

void add_supported_options(int argc, char** argv)
{
    add_UniGen_options();
    p.add("input", 1);

    try {
        po::store(po::command_line_parser(argc, argv).options(help_options).positional(p).run(), vm);
        if (vm.count("help"))
        {
            cout
            << "Probably Approximate counter" << endl;

            cout
            << "approxmc [options] inputfile" << endl << endl;

            cout << help_options << endl;
            std::exit(0);
        }

        if (vm.count("version")) {
            ApproxMC::AppMC tmp;
            UniG tmp2(&tmp);
            cout << tmp2.get_version_info();
            std::exit(0);
        }

        po::notify(vm);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::unknown_option> >& c
    ) {
        cerr
        << "ERROR: Some option you gave was wrong. Please give '--help' to get help" << endl
        << "       Unkown option: " << c.what() << endl;
        std::exit(-1);
    } catch (boost::bad_any_cast &e) {
        std::cerr
        << "ERROR! You probably gave a wrong argument type" << endl
        << "       Bad cast: " << e.what()
        << endl;

        std::exit(-1);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::invalid_option_value> >& what
    ) {
        cerr
        << "ERROR: Invalid value '" << what.what() << "'" << endl
        << "       given to option '" << what.get_option_name() << "'"
        << endl;

        std::exit(-1);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::multiple_occurrences> >& what
    ) {
        cerr
        << "ERROR: " << what.what() << " of option '"
        << what.get_option_name() << "'"
        << endl;

        std::exit(-1);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::required_option> >& what
    ) {
        cerr
        << "ERROR: You forgot to give a required option '"
        << what.get_option_name() << "'"
        << endl;

        std::exit(-1);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::too_many_positional_options_error> >& what
    ) {
        cerr
        << "ERROR: You gave too many positional arguments. Only the input CNF can be given as a positional option." << endl;
        std::exit(-1);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::ambiguous_option> >& what
    ) {
        cerr
        << "ERROR: The option you gave was not fully written and matches" << endl
        << "       more than one option. Please give the full option name." << endl
        << "       The option you gave: '" << what.get_option_name() << "'" <<endl
        << "       The alternatives are: ";
        for(size_t i = 0; i < what.alternatives().size(); i++) {
            cout << what.alternatives()[i];
            if (i+1 < what.alternatives().size()) {
                cout << ", ";
            }
        }
        cout << endl;

        std::exit(-1);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::invalid_command_line_syntax> >& what
    ) {
        cerr
        << "ERROR: The option you gave is missing the argument or the" << endl
        << "       argument is given with space between the equal sign." << endl
        << "       detailed error message: " << what.what() << endl
        ;
        std::exit(-1);
    }
}

template<class T>
void read_in_file(const string& filename, T* myreader)
{
    #ifndef USE_ZLIB
    FILE * in = fopen(filename.c_str(), "rb");
    DimacsParser<StreamBuffer<FILE*, FN>, T> parser(myreader, NULL, verbosity);
    #else
    gzFile in = gzopen(filename.c_str(), "rb");
    DimacsParser<StreamBuffer<gzFile, GZ>, T> parser(myreader, NULL, verbosity);
    #endif

    if (in == NULL) {
        std::cerr
        << "ERROR! Could not open file '"
        << filename
        << "' for reading: " << strerror(errno) << endl;

        std::exit(-1);
    }

    if (!parser.parse_DIMACS(in, true)) {
        exit(-1);
    }

    sampling_vars = parser.sampling_vars;
    sampling_vars_found = parser.sampling_vars_found;

    #ifndef USE_ZLIB
    fclose(in);
    #else
    gzclose(in);
    #endif
}

template<class T>
void read_stdin(T* myreader)
{
    cout
    << "c Reading from standard input... Use '-h' or '--help' for help."
    << endl;

    #ifndef USE_ZLIB
    FILE * in = stdin;
    #else
    gzFile in = gzdopen(0, "rb"); //opens stdin, which is 0
    #endif

    if (in == NULL) {
        std::cerr << "ERROR! Could not open standard input for reading" << endl;
        std::exit(1);
    }

    #ifndef USE_ZLIB
    DimacsParser<StreamBuffer<FILE*, FN>, T> parser(myreader, NULL, verbosity);
    #else
    DimacsParser<StreamBuffer<gzFile, GZ>, T> parser(myreader, NULL, verbosity);
    #endif

    if (!parser.parse_DIMACS(in, false)) {
        exit(-1);
    }

    sampling_vars = parser.sampling_vars;
    sampling_vars_found = parser.sampling_vars_found;

    #ifdef USE_ZLIB
    gzclose(in);
    #endif
}

void mycallback(const std::vector<int>& solution, void *file)
{
    std::ostream* os = (std::ostream*)file;
    for(uint32_t i = 0; i < solution.size(); i++) {
        (*os) << solution[i] <<  " ";
    }
    (*os) << "0" << endl;
}

template<class T>
void read_input_cnf(T* reader)
{
    //Init Arjun, read in file, get minimal indep set
    if (vm.count("input") != 0) {
        vector<string> inp = vm["input"].as<vector<string> >();
        if (inp.size() > 1) {
            cout << "[appmc] ERROR: you must only give one CNF as input" << endl;
            exit(-1);
        }

        read_in_file(inp[0].c_str(), reader);
    } else {
        read_stdin(reader);
    }
}

void print_orig_sampling_vars(const vector<uint32_t>& orig_sampling_vars)
{
    cout << "c [unig] Original sampling vars: ";
    for(auto v: orig_sampling_vars) cout << v << " ";
    cout << endl;
    cout << "c [unig] Orig sampling vars size: " << orig_sampling_vars.size() << endl;
}

void set_up_sampling_set()
{
    if (!sampling_vars_found || ignore_sampl_set) {
        if (!ignore_sampl_set) assert(sampling_vars.empty());
        sampling_vars.clear();
        for(uint32_t i = 0; i < arjun->nVars(); i++) sampling_vars.push_back(i);
    }
    arjun->set_starting_sampling_set(sampling_vars);
}

void get_cnf_from_arjun()
{
    bool ret = true;
    const uint32_t orig_num_vars = arjun->get_orig_num_vars();
    appmc->new_vars(orig_num_vars);
    arjun->start_getting_small_clauses(
        std::numeric_limits<uint32_t>::max(),
        std::numeric_limits<uint32_t>::max(),
        false);
    vector<Lit> clause;
    while (ret) {
        ret = arjun->get_next_small_clause(clause);
        if (!ret) {
            break;
        }

        bool ok = true;
        for(auto l: clause) {
            if (l.var() >= orig_num_vars) {
                ok = false;
                break;
            }
        }

        if (ok) {
            appmc->add_clause(clause);
        }

    }
    arjun->end_getting_small_clauses();
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
    cout << "c [arjun] final indep set: ";
    for(const uint32_t s: indep_set) cout << s+1 << " ";
    cout << "0" << endl;

    cout << "c [arjun] final set size: " << std::setw(8)
    << indep_set.size()
    << " percent of original: "
    <<  std::setw(6) << std::setprecision(4)
    << stats_line_percent(indep_set.size(), orig_sampling_set_size)
    << " %" << endl << std::flush;
}

int main(int argc, char** argv)
{
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

    appmc = new ApproxMC::AppMC;
    unigen = new UniG(appmc);
    add_supported_options(argc, argv);
    if (verbosity) {
        cout << unigen->get_version_info();
        cout << "c executed with command line: " << command_line << endl;
    }
    if (!only_indep_samples) {
        force_sol_extension = true;
        if (verbosity) {
            cout << "c Setting '--forcesolextension' to 1 since '--nosolext' is set to 0" << endl;
        }
    }

    //Main options
    appmc->set_verbosity(verbosity);
    if (verbosity) {
        appmc->set_detach_warning();
    }
    appmc->set_seed(seed);

    //Improvement options
    appmc->set_detach_xors(detach_xors);
    appmc->set_reuse_models(reuse_models);
    appmc->set_force_sol_extension(force_sol_extension);
    appmc->set_sparse(sparse);

    //Misc options
    appmc->set_simplify(simplify);
    appmc->set_var_elim_ratio(var_elim_ratio);

    if (logfilename != "") {
        appmc->set_up_log(logfilename);
        cout << "c [appmc] Logfile set " << logfilename << endl;
    }

    vector<uint32_t> empty_occ_sampl_vars;
    vector<uint32_t> sampling_vars_orig;
    if (do_arjun) {
        arjun = new ArjunNS::Arjun;
        arjun->set_seed(seed);
        arjun->set_verbosity(verbosity);
        if (verbosity) {
            cout << "c Arjun SHA revision " <<  arjun->get_version_info() << endl;
        }

        read_input_cnf(arjun);
        set_up_sampling_set();
        sampling_vars_orig = sampling_vars;
        check_sanity_sampling_vars(sampling_vars, arjun->get_orig_num_vars());
        print_orig_sampling_vars(sampling_vars_orig);
        get_cnf_from_arjun();
        sampling_vars = arjun->get_indep_set();
        print_final_indep_set(sampling_vars , sampling_vars_orig.size());
        if (debug_arjun) sampling_vars = sampling_vars_orig;
        delete arjun;
    } else {
        read_input_cnf(appmc);
        if (!sampling_vars_found || ignore_sampl_set) {
            sampling_vars.clear();
            for(uint32_t i = 0; i < appmc->nVars(); i++) sampling_vars.push_back(i);
        }
        sampling_vars_orig = sampling_vars;
        //print_orig_sampling_vars(sampling_vars, appmc);
    }

    appmc->set_projection_set(sampling_vars);
    auto sol_count = appmc->count();

    appmc->set_projection_set(sampling_vars_orig);
    unigen->set_verbosity(verbosity);
    unigen->set_verb_banning_cls(verb_banning_cls);
    unigen->set_kappa(kappa);
    unigen->set_multisample(multisample);
    unigen->set_only_indep_samples(only_indep_samples);
    unigen->set_force_sol_extension(force_sol_extension);

    std::ofstream logfile;
    if (logfilename != "") {
        logfile.open(logfilename.c_str());
        if (!logfile.is_open()) {
            cout << "[Sampler] Cannot open Sampler log file '" << logfilename
                 << "' for writing." << endl;
            exit(1);
        }
        unigen->set_logfile(&logfile);
    }

    void *myfile = &std::cout;
    std::ofstream sample_out;
    if (vm.count("sampleout") != 0) {
        sample_out.open(sample_fname.c_str());
        if (!sample_out.is_open()) {
            cout << "[Sampler] Cannot open samples file '" << sample_fname
                 << "' for writing." << endl;
            exit(-1);
        }
        myfile = &sample_out;
    }

    unigen->set_callback(mycallback, myfile);
    unigen->sample(&sol_count, num_samples);

    delete unigen;
    delete appmc;

    return 0;
}
