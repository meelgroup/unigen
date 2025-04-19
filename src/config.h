/*
 ApproxMC and gen_n_samples

 Copyright (c) 2009-2018, Mate Soos. All rights reserved.
 Copyright (c) 2014, Supratik Chakraborty, Kuldeep S. Meel, Moshe Y. Vardi
 Copyright (c) 2015, Supratik Chakraborty, Daniel J. Fremont,
 Kuldeep S. Meel, Sanjit A. Seshia, Moshe Y. Vardi

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


#pragma once

#include <vector>
#include <cstdint>

struct Config {
    unsigned verb = 1;
    unsigned verb_sampler_cls = 0;
    uint32_t multisample = 1;
    int force_sol_extension = 0;
    double kappa = 0.638;  /* Corresponds to UniGen's epsilon=16 in the TACAS-15 paper */
    std::vector<uint32_t> full_sampling_vars;
    std::vector<uint32_t> empty_sampling_vars;
};
