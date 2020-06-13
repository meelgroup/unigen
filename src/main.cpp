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
using boost::lexical_cast;
namespace po = boost::program_options;
using std::string;
using std::vector;

#if defined(__GNUC__) && defined(__linux__)
#include <fenv.h>
#endif
#include <signal.h>

#include "unigenconfig.h"
#include "unigen.h"
#include "time_mem.h"
#include "approxmc/approxmc.h"

#include <cryptominisat5/cryptominisat.h>
#include "cryptominisat5/dimacsparser.h"
#include "cryptominisat5/streambuffer.h"

using namespace CMSat;
using std::cout;
using std::cerr;
using std::endl;
string command_line;
UniGen* unigen = NULL;

UniGenConfig conf;
po::options_description UniGen_options = po::options_description("Main options");
po::options_description UniGengen_options = po::options_description("Sampling options");
po::options_description UniGen4_options = po::options_description("ApproxMC4 paper options");
po::options_description help_options;
po::variables_map vm;
po::positional_options_description p;

//signal code
void SIGINT_handler(int)
{
    return; //Perhaps we should output all the generated samples so far
}


void add_UniGen_options()
{
    std::ostringstream my_kappa;

    my_kappa << std::setprecision(8) << conf.kappa;

    UniGen_options.add_options()
    ("help,h", "Prints help")
    ("version", "Print version info")
    ("input", po::value< vector<string> >(), "file(s) to read")
    ("verb,v", po::value(&conf.verb)->default_value(conf.verb), "verbosity")
    ("seed,s", po::value(&conf.seed)->default_value(conf.seed), "Seed")
    ("start", po::value(&conf.startiter)->default_value(conf.startiter),
         "Start at this many XORs")
    ("log", po::value(&conf.logfilename)->default_value(conf.logfilename),
         "Logs of ApproxMC execution")
    ("th", po::value(&conf.num_threads)->default_value(conf.num_threads),
         "How many solving threads to use per solver call")
    ("vcl", po::value(&conf.verb_UniGen_cls)->default_value(conf.verb_UniGen_cls)
        ,"Print banning clause + xor clauses. Highly verbose.")
    ("simplify", po::value(&conf.simplify)->default_value(conf.simplify)
        , "Simplify agressiveness")
    //blasted_TR_ptb_1_linear.cnf.gz.no_w.cnf.gz is sensitive to below.
    //1.0 will mess it up. 0.3 will work.
    ("velimratio", po::value(&conf.var_elim_ratio)->default_value(conf.var_elim_ratio)
        , "Variable elimination ratio for each simplify run")
    ;

    //Improvements from ApproxMC4 paper
    UniGen4_options.add_options()
    ("detachxor", po::value(&conf.cms_detach_xor)->default_value(conf.cms_detach_xor)
        , "Detach XORs in CMS")
    ("reusemodels", po::value(&conf.reuse_models)->default_value(conf.reuse_models)
        , "Reuse models while counting solutions")
    ("forcesolextension", po::value(&conf.force_sol_extension)->default_value(conf.force_sol_extension)
        , "Use trick of not extending solutions in the SAT solver to full solution");

    UniGengen_options.add_options()
    ("samples", po::value(&conf.samples)->default_value(conf.samples)
        , "Number of random samples to generate")
    ("indepsamples", po::value(&conf.only_indep_samples)->default_value(conf.only_indep_samples)
        , "Should only output the independent vars from the samples")
    ("multisample", po::value(&conf.multisample)->default_value(conf.multisample)
        , "Return multiple samples from each call")
    ("sampleout", po::value(&conf.sampleFilename)
        , "Write samples to this file")
    ("kappa", po::value(&conf.kappa)->default_value(conf.kappa, my_kappa.str())
        , "Uniformity parameter (see TACAS-15 paper)")

    ;

    help_options.add(UniGen_options);
    help_options.add(UniGengen_options);
    help_options.add(UniGen4_options);
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
            unigen->printVersionInfo();
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

     if (conf.samples == 0 && vm.count("sampleout")){
        cerr << "ERROR: You did not give the '--samples N' option, but you gave the '--sampleout FNAME' option." << endl;
        cout << "ERROR: This is confusing. Please give '--samples N' if you give '--sampleout FNAME'" << endl;
        exit(-1);
     }
}

void readInAFile(SATSolver* solver2, const string& filename)
{
    solver2->add_sql_tag("filename", filename);
    #ifndef USE_ZLIB
    FILE * in = fopen(filename.c_str(), "rb");
    DimacsParser<StreamBuffer<FILE*, FN> > parser(solver, NULL, 2);
    #else
    gzFile in = gzopen(filename.c_str(), "rb");
    DimacsParser<StreamBuffer<gzFile, GZ> > parser(unigen->solver, NULL, 2);
    #endif

    if (in == NULL) {
        std::cerr
        << "ERROR! Could not open file '"
        << filename
        << "' for reading: " << strerror(errno) << endl;

        std::exit(-1);
    }

    if (!parser.parse_DIMACS(in, false)) {
        exit(-1);
    }

    conf.sampling_set.swap(parser.sampling_vars);

    #ifndef USE_ZLIB
        fclose(in);
    #else
        gzclose(in);
    #endif
}

void readInStandardInput(SATSolver* solver2)
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
    DimacsParser<StreamBuffer<FILE*, FN> > parser(solver2, NULL, 2);
    #else
    DimacsParser<StreamBuffer<gzFile, GZ> > parser(solver2, NULL, 2);
    #endif

    if (!parser.parse_DIMACS(in, false)) {
        exit(-1);
    }

    conf.sampling_set.swap(parser.sampling_vars);

    #ifdef USE_ZLIB
        gzclose(in);
    #endif
}

void set_sampling_vars()
{
    if (conf.sampling_set.empty()) {
        cout
        << "c [UniGen] WARNING! Sampling set was not declared with 'c ind var1 [var2 var3 ..] 0'"
        " notation in the CNF." << endl
        << "c [UniGen] we may work substantially worse!" << endl;
        for (size_t i = 0; i < unigen->solver->nVars(); i++) {
            conf.sampling_set.push_back(i);
        }
    }
    cout << "c [UniGen] Sampling set size: " << conf.sampling_set.size() << endl;
    if (conf.sampling_set.size() > 100) {
        cout
        << "c [UniGen] Sampling var set contains over 100 variables, not displaying"
        << endl;
    } else {
        cout << "c [UniGen] Sampling set: ";
        for (auto v: conf.sampling_set) {
            cout << v+1 << ", ";
        }
        cout << endl;
    }
    unigen->solver->set_sampling_vars(&conf.sampling_set);
}

std::ostream* open_samples_file()
{
    std::ostream* os;
    std::ofstream* sampleFile = NULL;
    if (conf.sampleFilename.length() != 0)
    {
        sampleFile = new std::ofstream;
        sampleFile->open(conf.sampleFilename.c_str());
        if (!(*sampleFile)) {
            cout
            << "ERROR: Couldn't open sample file '"
            << conf.sampleFilename
            << "' for writing!"
            << endl;
            std::exit(-1);
        }
        os = sampleFile;
    } else {
        os = &cout;
    }

    return os;
}

int main(int argc, char** argv)
{
    #if defined(__GNUC__) && defined(__linux__)
    feenableexcept(FE_INVALID   |
                   FE_DIVBYZERO |
                   FE_OVERFLOW
                  );
    #endif
    signal(SIGINT, SIGINT_handler);
    signal(SIGALRM, SIGINT_handler);
    signal(SIGTERM, SIGINT_handler);
    signal(SIGKILL, SIGINT_handler);
    //Reconstruct the command line so we can emit it later if needed
    for(int i = 0; i < argc; i++) {
        command_line += string(argv[i]);
        if (i+1 < argc) {
            command_line += " ";
        }
    }

    UniGen *unigen = new UniGen;
    add_supported_options(argc, argv);
    unigen->printVersionInfo();
    cout
    << "c executed with command line: "
    << command_line
    << endl;

    cout << "c [UniGen] using seed: " << conf.seed << endl;

    if (vm.count("log") == 0) {
        if (vm.count("input") != 0) {
            conf.logfilename = vm["input"].as<vector<string> >()[0] + ".log";
            cout << "c [UniGen] Logfile name not given, assumed to be " << conf.logfilename << endl;
        } else {
            std::cerr << "[UniGen] ERROR: You must provide the logfile name" << endl;
            exit(-1);
        }
    }

    if (!conf.only_indep_samples) {
        cout << "ERROR: You requested samples with full solutions but '--cmpindeponly 1' is set. Set it to false: '--indep 0'" << endl;
        exit(-1);
    }

    //startTime = cpuTimeTotal();
    unigen->solver = new SATSolver();
    unigen->solver->set_up_for_scalmc();

    if (conf.verb > 2) {
        unigen->solver->set_verbosity(conf.verb-2);
    }
    unigen->solver->set_allow_otf_gauss();
    unigen->solver->set_xor_detach(conf.cms_detach_xor);

    if (conf.num_threads > 1) {
        unigen->solver->set_num_threads(conf.num_threads);
    }

    
    if (conf.samples < 0) {
        cout << "[UniGen] ERROR: The number of samples should be greater than zero" << endl;
    }

    //parsing the input
    if (vm.count("input") != 0) {
        vector<string> inp = vm["input"].as<vector<string> >();
        if (inp.size() > 1) {
            cout << "[UniGen] ERROR: can only parse in one file" << endl;
        }
        readInAFile(unigen->solver, inp[0].c_str());
    } else {
        readInStandardInput(unigen->solver);
    }
    set_sampling_vars();

    std::ostream* out = open_samples_file();
    unigen->set_samples_file(out);

    //Counting
    if (conf.samples > 0 && conf.startiter == 0) {
        cout << "c [UniGen] Using UniGen to compute startiter for gen_n_samples" << endl;
    }

    if (conf.startiter > conf.sampling_set.size()) {
        cout << "c [UniGen] ERROR: Manually-specified start_iter"
             "is larger than the size of the sampling set.\n" << endl;
        exit(-1);
    }

    auto ret = unigen->solve(conf);
    if (out != &cout) {
        delete out;
    }

    return ret;
}
