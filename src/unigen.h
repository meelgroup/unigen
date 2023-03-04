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

#ifndef UNIGEN_H__
#define UNIGEN_H__

#include <cstdint>
#include <string>
#include <vector>
#include <fstream>
#include <functional>
#ifdef CMS_LOCAL_BUILD
#include "cryptominisat.h"
#else
#include <cryptominisat5/cryptominisat.h>
#endif

namespace ApproxMC {
    class AppMC;
    class SolCount;
}

namespace UniGen {

typedef std::function<void(const std::vector<int>& solution, void* data)> callback;

struct UniGenPrivateData;
#ifdef _WIN32
class __declspec(dllexport) UniG
#else
class UniG
#endif
{
public:
    UniG(ApproxMC::AppMC* appmc);
    ~UniG();
    std::string get_version_info();
    void sample(
        const ApproxMC::SolCount* sol_count,
        uint32_t num_samples);

    //Misc options -- do NOT to change unless you know what you are doing!
    void set_kappa(double kappa);
    void set_multisample(bool multisample);
    void set_only_indep_samples(bool only_indep_samples);
    void set_verb_sampler_cls(bool verb_sampler_cls);
    void set_force_sol_extension(bool force_sol_extension);
    void set_logfile(std::ostream* logfile);
    void set_verbosity(uint32_t verb);
    void set_callback(UniGen::callback f, void* data);
    void set_full_sampling_vars(const std::vector<uint32_t>& full_sampl_vars);

    //Querying default values
    double get_kappa() const;
    bool get_multisample() const;
    bool get_only_indep_samples() const;
    bool get_verb_sampler_cls() const;
    bool get_force_sol_extension() const;
    bool get_verb_banning_cls() const;
    std::ostream* get_logfile() const;
    const std::vector<uint32_t>& get_full_sampling_vars() const;

private:
    ////////////////////////////
    // Do not bother with this, it's private
    ////////////////////////////
    UniGenPrivateData* data;
};

}

#endif
