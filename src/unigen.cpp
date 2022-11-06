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
#include "config.h"
#include <iostream>
#include <fstream>
using namespace ApproxMC;

using std::cout;
using std::endl;

#if defined _WIN32
    #define DLL_PUBLIC __declspec(dllexport)
#else
    #define DLL_PUBLIC __attribute__ ((visibility ("default")))
    #define DLL_LOCAL  __attribute__ ((visibility ("hidden")))
#endif

namespace UniGen {
    struct UniGenPrivateData {
        Sampler sampler;
        AppMC* appmc;
        Config conf;
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

DLL_PUBLIC void UniG::sample(
    const SolCount* sol_count,
    uint32_t num_samples)
{
    if (data->sampler.callback_func == NULL) {
        std::cout << "ERROR! You must set the callback function or your samples will be lost" << endl;
        exit(-1);
    }
    data->sampler.sample(data->conf, *sol_count, num_samples);
}

DLL_PUBLIC string UniG::get_version_info()
{
    return data->sampler.get_version_info();
}

DLL_PUBLIC double UniG::get_epsilon()
{
    return data->conf.epsilon;
}

DLL_PUBLIC bool UniG::get_multisample()
{
    return data->conf.multisample;
}

DLL_PUBLIC bool UniG::get_only_indep_samples()
{
    return data->conf.only_indep_samples;
}

DLL_PUBLIC bool UniG::get_verb_sampler_cls()
{
    return data->conf.verb_banning_cls;
}

DLL_PUBLIC void UniG::set_epsilon(double epsilon)
{
    data->conf.epsilon = epsilon;
}

DLL_PUBLIC void UniG::set_multisample(bool multisample)
{
    data->conf.multisample = multisample;
}

DLL_PUBLIC void UniG::set_only_indep_samples(bool only_indep_samples)
{
    data->conf.only_indep_samples = only_indep_samples;
}

DLL_PUBLIC void UniG::set_verb_banning_cls(bool verb_banning_cls)
{
    data->conf.verb_banning_cls = verb_banning_cls;
}

DLL_PUBLIC bool UniG::get_force_sol_extension()
{
    return data->conf.force_sol_extension;
}

DLL_PUBLIC void UniG::set_force_sol_extension(bool force_sol_extension)
{
    data->conf.force_sol_extension = force_sol_extension;
}

DLL_PUBLIC void UniG::set_logfile(std::ostream* logfile)
{
    data->conf.logfile = logfile;
}

DLL_PUBLIC void UniG::set_verbosity(uint32_t verb)
{
    data->conf.verb = verb;
}
