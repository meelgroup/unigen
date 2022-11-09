/*************
Python bindings for Unigen, heavily based on the Python bindings written for CryptoMiniSat

Copyright (c) 2021 Eric Vin
              2022 Mate Soos

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
**********************************/

#include <Python.h>
#include "../cryptominisat/src/cryptominisat.h"
#include "../approxmc/src/approxmc.h"
#include "../../src/unigen.h"

#include <limits>

#define MODULE_NAME "pyunigen"
#define MODULE_DOC "Unigen almost uniform sampler."

struct PySampler {
    PyObject_HEAD
    UniGen::UniG* unig;
    PyObject* sample_list = NULL;
    ApproxMC::AppMC* appmc;

    // config
    int verbosity = 0;
    uint32_t seed = 1;
    double kappa;
    double epsilon;
    double delta;
    bool multisample = false;

    // internal
    std::vector<CMSat::Lit> tmp_cl_lits;
    std::vector<uint32_t> sampling_set;
    bool called_already = false;
    uint32_t samples_needed = 5;
    uint32_t samples_generated = 0;
};

static const char sampler_create_docstring[] = \
"Sampler(verbosity=0, seed=1)\n\
Create Sampler object.\n\
\n\
:param verbosity: Verbosity level: 0: nothing printed; 15: very verbose.\n\
:param seed: Random seed\n\
:param delta: \n\
:param epsilon: \n\
:param kappa: Uniformity parameter (see TACAS-15 paper)\n\
";

/********** Internal Functions **********/

/* Helper functions */

void pybinding_callback(const std::vector<int>& solution, void *self_in)
{
    PySampler* self = (PySampler*) self_in;
    if (self->samples_generated >= self->samples_needed) return;

    PyObject* sample = PyList_New(solution.size());
    if (sample == NULL) {
        PyErr_SetString(PyExc_SystemError, "failed to create a list");
        return;
    }

    for (unsigned int i = 0; i < solution.size(); i++) {
        PyObject *lit = PyLong_FromLong((long)solution[i]);
        if (lit == NULL) {
            PyErr_SetString(PyExc_SystemError, "failed to create a list");
            return;
        }
        PyList_SET_ITEM(sample, i, lit);
    }
    PyList_Append(self->sample_list, sample);
    self->samples_generated++;
}

static int parse_sampling_set(PySampler *self, PyObject *sampling_set_obj)
{
    PyObject *iterator = PyObject_GetIter(sampling_set_obj);
    if (iterator == NULL) {
        PyErr_SetString(PyExc_TypeError, "iterable object expected");
        return NULL;
    }

    PyObject *lit;
    while ((lit = PyIter_Next(iterator)) != NULL) {
        long val = PyLong_AsLong(lit);
        if (val <= 0) {
            PyErr_SetString(PyExc_ValueError, "Sampling set must be positive numbers");
            return NULL;
        }
        if (val > std::numeric_limits<int>::max()/2
            || val < std::numeric_limits<int>::min()/2
        ) {
            PyErr_Format(PyExc_ValueError, "integer %ld is too small or too large", val);
            return NULL;
        }

        self->sampling_set.push_back(val - 1);
        Py_DECREF(lit);
    }
    Py_DECREF(iterator);

    return 1;
}

static void setup_sampler(PySampler *self, PyObject *args, PyObject *kwds)
{
    self->verbosity = 0;
    self->seed = 1;
    self->multisample = false;
    self->called_already = false;
    self->samples_needed = 5;
    self->samples_generated = 0;
    self->appmc = new ApproxMC::AppMC;
    self->unig = new UniGen::UniG(self->appmc);
    self->epsilon = self->appmc->get_epsilon();
    self->delta = self->appmc->get_delta();
    self->kappa = self->unig->get_kappa();

    static char* kwlist[] = {"verbosity", "seed", "epsilon", "delta", "kappa", "multisample", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|iIdddp", kwlist,
        &self->verbosity, &self->seed, &self->epsilon, &self->delta, &self->kappa, &self->multisample))
    {
        return;
    }

    if (self->verbosity < 0) {
        PyErr_SetString(PyExc_ValueError, "verbosity must be at least 0");
        return;
    }
    if (self->epsilon <= 0 || self->epsilon >= 1) {
        PyErr_SetString(PyExc_ValueError, "epsilon must be greater than 0");
        return;
    }
    if (self->delta <= 0 || self->delta >= 1) {
        PyErr_SetString(PyExc_ValueError, "delta must be greater than 0");
        return;
    }

    if (self->kappa <= 0 || self->kappa >= 1) {
        PyErr_SetString(PyExc_ValueError, "kappa must be greater than 0");
        return;
    }

    self->appmc->set_verbosity(self->verbosity);
    self->appmc->set_seed(self->seed);
    self->appmc->set_epsilon(self->epsilon);
    self->appmc->set_delta(self->delta);

    self->unig->set_verbosity(self->verbosity);
    self->unig->set_kappa(self->kappa);
    self->unig->set_multisample(self->multisample);

    self->unig->set_callback(pybinding_callback, self);

    return;
}

static int convert_lit_to_sign_and_var(PyObject* lit, long& var, bool& sign)
{
    if (!PyLong_Check(lit))  {
        PyErr_SetString(PyExc_TypeError, "integer expected !");
        return 0;
    }

    long val = PyLong_AsLong(lit);
    if (val == 0) {
        PyErr_SetString(PyExc_ValueError, "non-zero integer expected");
        return 0;
    }
    if (val > std::numeric_limits<int>::max()/2
        || val < std::numeric_limits<int>::min()/2
    ) {
        PyErr_Format(PyExc_ValueError, "integer %ld is too small or too large", val);
        return 0;
    }

    sign = (val < 0);
    var = std::abs(val) - 1;

    return 1;
}

static int parse_clause(PySampler *self, PyObject *clause, std::vector<CMSat::Lit>& lits)
{
    PyObject *iterator = PyObject_GetIter(clause);
    if (iterator == NULL) {
        PyErr_SetString(PyExc_TypeError, "iterable object expected");
        return 0;
    }

    PyObject *lit;
    long int max_var = 0;
    while ((lit = PyIter_Next(iterator)) != NULL) {
        long var;
        bool sign;
        int ret = convert_lit_to_sign_and_var(lit, var, sign);
        Py_DECREF(lit);
        if (!ret) {
            Py_DECREF(iterator);
            return 0;
        }
        max_var = std::max(var, max_var);

        lits.push_back(CMSat::Lit(var, sign));
    }

    if (!lits.empty() && max_var >= (long int)self->appmc->nVars()) {
        self->appmc->new_vars(max_var-(long int)self->appmc->nVars()+1);
    }

    Py_DECREF(iterator);
    if (PyErr_Occurred()) {
        return 0;
    }

    return 1;
}

static int _add_clause(PySampler *self, PyObject *clause)
{
    self->tmp_cl_lits.clear();
    if (!parse_clause(self, clause, self->tmp_cl_lits)) {
        return 0;
    }
    self->appmc->add_clause(self->tmp_cl_lits);

    return 1;
}

/* add_clause function */

PyDoc_STRVAR(add_clause_doc,
"add_clause(clause)\n\
Add a clause to the solver.\n\
\n\
:param clause: An iterator contains literals (ints)"
);

static PyObject* add_clause(PySampler *self, PyObject *args, PyObject *kwds)
{
    static char* kwlist[] = {"clause", NULL};
    PyObject *clause;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O", kwlist, &clause)) {
        return NULL;
    }

    if (_add_clause(self, clause) == 0 ) {
        return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;

}

static int parse_cell_hash_count(PyObject* obj, ApproxMC::SolCount& sol_count)
{
    if (!PyTuple_Check(obj)) {
        PyErr_SetString(PyExc_SystemError, "The hash&cell count must be a tuple");
        return NULL;
    }
    if (PyTuple_Size(obj) != 2) {
        PyErr_SetString(PyExc_SystemError, "The hash&cell count must be a tuple of size 2");
        return NULL;
    }

    PyObject* tmp = PyTuple_GetItem(obj, 0);
    if (!PyLong_Check(tmp))  {
        PyErr_SetString(PyExc_TypeError, "integer expected for cell count");
        return 0;
    }
    long cell_count = PyLong_AsLong(tmp);
    if (cell_count < 0) {
        PyErr_SetString(PyExc_TypeError, ">=0 integer expected for cell count");
        return 0;
    }

    tmp = PyTuple_GetItem(obj, 1);
    if (!PyLong_Check(tmp))  {
        PyErr_SetString(PyExc_TypeError, "integer expected for hash count");
        return 0;
    }
    long hash_count = PyLong_AsLong(tmp);
    if (hash_count < 0) {
        PyErr_SetString(PyExc_TypeError, ">=0 integer expected for hash count");
        return 0;
    }

    sol_count.cellSolCount = cell_count;
    sol_count.hashCount = hash_count;
    return 1;
}

PyDoc_STRVAR(sample_doc,
"sample()\n\
Sample almost uniformly from the solutions for the clauses that have been \n\
added with add_clause().\n\
\n\
:param sampling_set: (Optional) If provided, the solutions are sampled almost uniformly\n\
    over the variables in sampling_set.\n\
:param cell_hash_count: (Optional) Pair of hash_count and cell_count. If provided, they are used instead of\n\
    calling ApproxMC internally.\n\
:return: The Python tuple (cell count, hash count, list of samples)\n\
    where the first two elements are from ApproxMC's count, and the last\n\
    element is a list of samples"
);

static PyObject* sample(PySampler *self, PyObject *args, PyObject *kwds)
{
    if (self->called_already) {
        PyErr_SetString(PyExc_SystemError, "You can only call sample() once");
        return NULL;
    }
    self->called_already = true;

    PyObject* sampling_set_obj = NULL;
    PyObject* cell_hash_count_obj = NULL;
    static char* kwlist[] = {"num", "sampling_set", "cell_hash_count", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|IOO", kwlist,
    	&self->samples_needed, &sampling_set_obj, &cell_hash_count_obj))
    {
        return NULL;
    }

    if (sampling_set_obj != NULL && !parse_sampling_set(self, sampling_set_obj)) {
        return NULL;
    }

    ApproxMC::SolCount sol_count;
    if (cell_hash_count_obj != NULL && !parse_cell_hash_count(cell_hash_count_obj, sol_count)) {
        return NULL;
    }

    if (sampling_set_obj == NULL) {
        assert(self->sampling_set.empty());
        for(uint32_t i = 0; i < self->appmc->nVars(); i++) self->sampling_set.push_back(i);
    }

    self->appmc->set_projection_set(self->sampling_set);

    if (cell_hash_count_obj == NULL) {
    	sol_count = self->appmc->count();
    }

    self->sample_list = PyList_New(0);
    if (self->sample_list == NULL) {
        PyErr_SetString(PyExc_SystemError, "failed to create a list");
        return NULL;
    }
    self->unig->sample(&sol_count, self->samples_needed);

    PyObject *result = PyTuple_New((Py_ssize_t) 3);
    if (result == NULL) {
        PyErr_SetString(PyExc_SystemError, "failed to create a tuple");
        return NULL;
    }
    PyTuple_SET_ITEM(result, 0, PyLong_FromLong((long)sol_count.cellSolCount));
    PyTuple_SET_ITEM(result, 1, PyLong_FromLong((long)sol_count.hashCount));
    PyTuple_SET_ITEM(result, 2, self->sample_list);

    return result;
}

/********** Python Bindings **********/
static PyMethodDef PySampler_methods[] = {
    {"sample",    (PyCFunction) sample,      METH_VARARGS | METH_KEYWORDS, sample_doc},
    {"add_clause",(PyCFunction) add_clause,  METH_VARARGS | METH_KEYWORDS, add_clause_doc},
    {NULL, NULL}  // Sentinel
};

static void PySampler_dealloc(PySampler* self)
{
    delete self->unig;
    delete self->appmc;
    Py_TYPE(self)->tp_free ((PyObject*) self);
}

static int PySampler_init(PySampler *self, PyObject *args, PyObject *kwds)
{
    if (self->unig != NULL) {
        delete self->unig;
    }

    if (self->appmc != NULL) {
        delete self->appmc;
    }

    setup_sampler(self, args, kwds);

    if (!self->unig) {
        return -1;
    }

    if (!self->appmc) {
        return -1;
    }

    return 0;
}

static PyTypeObject pyunigen_PySamplerType =
{
    PyVarObject_HEAD_INIT(NULL, 0)  /*ob_size*/
    "pyunigen.Sampler",           /*tp_name*/
    sizeof(PySampler),                /*tp_basicsize*/
    0,                              /*tp_itemsize*/
    (destructor)PySampler_dealloc,    /*tp_dealloc*/
    0,                              /*tp_print*/
    0,                              /*tp_getattr*/
    0,                              /*tp_setattr*/
    0,                              /*tp_compare*/
    0,                              /*tp_repr*/
    0,                              /*tp_as_number*/
    0,                              /*tp_as_sequence*/
    0,                              /*tp_as_mapping*/
    0,                              /*tp_hash */
    0,                              /*tp_call*/
    0,                              /*tp_str*/
    0,                              /*tp_getattro*/
    0,                              /*tp_setattro*/
    0,                              /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    sampler_create_docstring,       /* tp_doc */
    0,                              /* tp_traverse */
    0,                              /* tp_clear */
    0,                              /* tp_richcompare */
    0,                              /* tp_weaklistoffset */
    0,                              /* tp_iter */
    0,                              /* tp_iternext */
    PySampler_methods,                /* tp_methods */
    0,                              /* tp_members */
    0,                              /* tp_getset */
    0,                              /* tp_base */
    0,                              /* tp_dict */
    0,                              /* tp_descr_get */
    0,                              /* tp_descr_set */
    0,                              /* tp_dictoffset */
    (initproc)PySampler_init,         /* tp_init */
};

PyMODINIT_FUNC PyInit_pyunigen(void)
{
    PyObject* m;

    pyunigen_PySamplerType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&pyunigen_PySamplerType) < 0) {
        // Return NULL on Python3 and on Python2 with MODULE_INIT_FUNC macro
        // In pure Python2: return nothing.
        return NULL;
    }

    static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,  /* m_base */
        MODULE_NAME,            /* m_name */
        MODULE_DOC,             /* m_doc */
        -1,                     /* m_size */
        NULL,                   /* m_methods */
        NULL,                   /* m_reload */
        NULL,                   /* m_traverse */
        NULL,                   /* m_clear */
        NULL,                   /* m_free */
    };

    m = PyModule_Create(&moduledef);

    // Add the version string so users know what version of UniGen
    // they're using.
    if (PyModule_AddStringConstant(m, "__version__", UNIGEN_FULL_VERSION) == -1) {
        Py_DECREF(m);
        return NULL;
    }
    if (PyModule_AddStringConstant(m, "VERSION", UNIGEN_FULL_VERSION) == -1) {
        Py_DECREF(m);
        return NULL;
    }

    if (!m) {
        return NULL;
    }

    // Add the PySampler type
    Py_INCREF(&pyunigen_PySamplerType);
    if (PyModule_AddObject(m, "Sampler", (PyObject *)&pyunigen_PySamplerType)) {
        Py_DECREF(m);
        return NULL;
    }

    return m;
}

